// Plotの数学的解析IRと軽量symbolic prepass
#include "plot_analysis.hpp"

#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "error/error_message.hpp"
#include "expression/exact_value.hpp"
#include "mathematics/exact_algebra.hpp"
#include "mathematics/exact_trigonometry.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "solver/polynomial_solver.hpp"
#include "solver/real_function_analysis.hpp"
#include "solver/solve_constraints.hpp"
#include "solver/solver_support.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::plot {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::RelationKind;
using mathematics::TruthValue;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

inline constexpr std::size_t maxSymbolicLandmarks = 256;

[[nodiscard]] Expr integerExpr(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] std::optional<int> compareExact(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (lhs == rhs)
        return 0;
    if (const auto left = expression::exact::realRational(lhs)) {
        if (const auto right = expression::exact::realRational(rhs))
            return *left < *right ? -1 : (*right < *left ? 1 : 0);
    }

    // Plotの定数順序はexactな代数数同士の比較を毎回構成する必要がない。
    // certified enclosureが分離すればそれ自体が証明なので，安価な数値証明を先に使う。
    // equality等でenclosureだけでは決まらない場合に限りsymbolic proofへ戻す。
    const approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    for (const std::size_t bits : {96U, 192U, 384U}) {
        try {
            const auto left = certified.enclose(lhs, bits);
            const auto right = certified.enclose(rhs, bits);
            if (!left || !right || !left->isReal() || !right->isReal()) {
                break;
            }
            const auto& l = left->asReal();
            const auto& r = right->asReal();
            if (l.upper() < r.lower())
                return -1;
            if (l.lower() > r.upper())
                return 1;
            if (l.isPoint() && r.isPoint() && l.lower() == r.lower())
                return 0;
        }
        catch (const approximation::PrecisionInsufficient& exception) {
            if (!exception.refinable()) {
                break;
            }
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            break;
        }
    }

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::relation(RelationKind::Less, lhs, rhs))
        == TruthValue::True)
        return -1;
    if (knowledge.prove(mathematics::relation(RelationKind::Greater, lhs, rhs))
        == TruthValue::True)
        return 1;
    return std::nullopt;
}

[[nodiscard]] bool insideRequestedRange(
    const Expr& point,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto lower = compareExact(
        point, request.lower, builtins, mathematics, angles, assumptions);
    const auto upper = compareExact(
        point, request.upper, builtins, mathematics, angles, assumptions);
    return lower && upper && *lower >= 0 && *upper <= 0;
}

void addLandmark(
    PlotAnalysis& analysis,
    Expr position,
    PlotLandmarkKind kind,
    PlotLandmarkConfidence confidence = PlotLandmarkConfidence::Proven,
    std::optional<Expr> finiteLimit = std::nullopt) {
    const PlotLandmark landmark{
        std::move(position), kind, confidence, std::move(finiteLimit)};
    if (std::find(analysis.landmarks.begin(), analysis.landmarks.end(), landmark)
        == analysis.landmarks.end())
        analysis.landmarks.push_back(landmark);
}

[[nodiscard]] PlotEndpointInclusion inclusion(bool inclusive) noexcept {
    return inclusive ? PlotEndpointInclusion::Closed : PlotEndpointInclusion::Open;
}

[[nodiscard]] bool clipInterval(
    const solver::RealDomainInterval& source,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::optional<PlotInterval>& clipped) {
    PlotInterval result{
        request.lower,
        PlotEndpointInclusion::Closed,
        request.upper,
        PlotEndpointInclusion::Closed};

    if (source.lower) {
        const auto comparison = compareExact(
            *source.lower, result.lower, builtins, mathematics, angles, assumptions);
        if (!comparison)
            return false;
        if (*comparison > 0) {
            result.lower = *source.lower;
            result.lowerInclusion = inclusion(source.lowerInclusive);
        }
        else if (*comparison == 0 && !source.lowerInclusive)
            result.lowerInclusion = PlotEndpointInclusion::Open;
    }

    if (source.upper) {
        const auto comparison = compareExact(
            *source.upper, result.upper, builtins, mathematics, angles, assumptions);
        if (!comparison)
            return false;
        if (*comparison < 0) {
            result.upper = *source.upper;
            result.upperInclusion = inclusion(source.upperInclusive);
        }
        else if (*comparison == 0 && !source.upperInclusive)
            result.upperInclusion = PlotEndpointInclusion::Open;
    }

    const auto ordering = compareExact(
        result.lower, result.upper, builtins, mathematics, angles, assumptions);
    if (!ordering)
        return false;
    if (*ordering > 0
        || (*ordering == 0
            && (result.lowerInclusion == PlotEndpointInclusion::Open
                || result.upperInclusion == PlotEndpointInclusion::Open))) {
        clipped.reset();
        return true;
    }
    clipped = std::move(result);
    return true;
}

[[nodiscard]] std::optional<Expr> simplifiedAtPoint(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& point,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    try {
        Expr value = symbolic::substituteSymbol(expression, variable, point);
        value = solver::simplifyForSolve(
            std::move(value), builtins, mathematics, angles, assumptions);
        return value;
    }
    catch (const error::CalcError& exception) {
        // 境界点そのものが未定義でもdomain分割は既に証明済み。
        // semantic分類だけ諦め，ResourceLimit等は従来どおり上位へ伝播する。
        if (exception.type() == error::CalcErrorType::Domain)
            return std::nullopt;
        throw;
    }
}

[[nodiscard]] bool exactZero(const Expr& expression) {
    const auto value = expression::exact::realRational(expression);
    return value && value->isZero();
}

[[nodiscard]] bool exactMinusOne(const Expr& expression) {
    const auto value = expression::exact::realRational(expression);
    return value && *value == Rational{BigInt{-1}};
}

// Plotのdomain分割は既にexactに証明済みなので，境界の意味分類のためだけに
// 一般Limitを走らせない。深い合成函数ではLimitが描画本体より桁違いに重くなるため，
// 構造からcheapに証明できるvertical asymptote / poleだけを追加する。
void classifyBoundaryCheap(
    PlotAnalysis& analysis,
    const PlotRequest& request,
    const Expr& point,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!request.expression.isCall())
        return;
    const auto* definition = builtins.find(request.expression.asCall().head);
    if (!definition)
        return;
    const auto& arguments = request.expression.asCall().arguments;

    if ((definition->id == BuiltinId::Log
            || definition->id == BuiltinId::Log2
            || definition->id == BuiltinId::Log10)
        && arguments.size() == 1) {
        const auto argument = simplifiedAtPoint(
            arguments.front(), request.variable, point,
            builtins, mathematics, angles, assumptions);
        if (argument && exactZero(*argument))
            addLandmark(analysis, point, PlotLandmarkKind::VerticalAsymptote);
        return;
    }

    if (definition->id == BuiltinId::Log1p && arguments.size() == 1) {
        const auto argument = simplifiedAtPoint(
            arguments.front(), request.variable, point,
            builtins, mathematics, angles, assumptions);
        if (argument && exactMinusOne(*argument))
            addLandmark(analysis, point, PlotLandmarkKind::VerticalAsymptote);
        return;
    }

    if ((definition->id == BuiltinId::ExponentialIntegralEi
            || definition->id == BuiltinId::CosineIntegralCi)
        && arguments.size() == 1) {
        const auto argument = simplifiedAtPoint(
            arguments.front(), request.variable, point,
            builtins, mathematics, angles, assumptions);
        if (argument && exactZero(*argument))
            addLandmark(analysis, point, PlotLandmarkKind::VerticalAsymptote);
        return;
    }

    if (definition->id == BuiltinId::Power && arguments.size() == 2) {
        // x^x は実Plotでは x>0 で定義され，x->0+ で 1 へ収束する。
        // requestが非負側へclipされている場合，x=0 は removable hole として
        // 明示してよい。負側を含む request は上流で unsafe として弾かれる。
        if (point.isNumber() && point.asNumber().isReal() && point.asNumber().asReal().isInteger()
            && point.asNumber().asReal().asInteger().isZero()
            && arguments[0].isSymbol() && arguments[1].isSymbol()
            && arguments[0].asSymbol().sameIdentity(request.variable)
            && arguments[1].asSymbol().sameIdentity(request.variable)) {
            addLandmark(analysis, point, PlotLandmarkKind::RemovableSingularity,
                PlotLandmarkConfidence::Proven, integerExpr(1));
        }
        return;
    }

    if (definition->id != BuiltinId::Divide || arguments.size() != 2)
        return;

    const auto numerator = simplifiedAtPoint(
        arguments[0], request.variable, point,
        builtins, mathematics, angles, assumptions);
    const auto denominator = simplifiedAtPoint(
        arguments[1], request.variable, point,
        builtins, mathematics, angles, assumptions);
    if (!numerator || !denominator || !exactZero(*denominator))
        return;

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::relation(
            RelationKind::NotEqual, *numerator, integerExpr(0))) != TruthValue::True)
        return;

    addLandmark(analysis, point, PlotLandmarkKind::VerticalAsymptote);

    // Poleはmeromorphicな零点と証明できる場合だけ付ける。
    // sqrt/abs等の零点はvertical asymptoteでもpoleではない。
    const auto polynomial = symbolic::toExpressionPolynomial(
        arguments[1], request.variable, builtins, mathematics, angles);
    if (polynomial && polynomial->degree() >= 1)
        addLandmark(analysis, point, PlotLandmarkKind::Pole);
}

[[nodiscard]] std::optional<Expr> directBranchPoint(
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins) {
    if (!request.expression.isCall()
        || request.expression.asCall().arguments.size() != 1)
        return std::nullopt;
    const Expr& argument = request.expression.asCall().arguments.front();
    if (!argument.isSymbol() || !argument.asSymbol().sameIdentity(request.variable))
        return std::nullopt;

    if (builtins.isCallTo(request.expression, BuiltinId::Sqrt)
        || builtins.isCallTo(request.expression, BuiltinId::Log)
        || builtins.isCallTo(request.expression, BuiltinId::Log2)
        || builtins.isCallTo(request.expression, BuiltinId::Log10)
        || builtins.isCallTo(request.expression, BuiltinId::ExponentialIntegralEi)
        || builtins.isCallTo(request.expression, BuiltinId::CosineIntegralCi))
        return integerExpr(0);
    if (builtins.isCallTo(request.expression, BuiltinId::Log1p))
        return integerExpr(-1);
    return std::nullopt;
}

[[nodiscard]] BigInt floorRational(const Rational& value) {
    auto result = numeric::divmod(value.numerator(), value.denominator());
    if (value.numerator().isNegative() && !result.remainder.isZero())
        result.quotient -= BigInt{1};
    return result.quotient;
}

[[nodiscard]] BigInt ceilRational(const Rational& value) {
    auto result = numeric::divmod(value.numerator(), value.denominator());
    if (value.numerator().isPositive() && !result.remainder.isZero())
        result.quotient += BigInt{1};
    return result.quotient;
}

[[nodiscard]] Expr angleFromTurns(
    const Rational& turns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Degree:
        return Expr{Number{turns * Rational{BigInt{360}}}};
    case mathematics::AngleUnit::Gradian:
        return Expr{Number{turns * Rational{BigInt{400}}}};
    case mathematics::AngleUnit::Radian: {
        const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
        if (!pi)
            return integerExpr(0);
        return mathematics::scaleExactExpression(
            turns * Rational{BigInt{2}}, Expr{pi->symbol}, builtins);
    }
    }
    return integerExpr(0);
}


struct StepDiscontinuity final {
    Expr position;
    bool attachLeft = false;
    bool attachRight = false;
    PlotLandmarkConfidence confidence = PlotLandmarkConfidence::Proven;
};

[[nodiscard]] std::optional<std::pair<Rational, Rational>> certifiedConstantBounds(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    for (const std::size_t bits : {96U, 192U}) {
        try {
            const auto value = certified.enclose(expression, bits);
            if (!value || !value->isReal())
                return std::nullopt;
            const auto& interval = value->asReal();
            return std::pair<Rational, Rational>{
                interval.lower().toRational(), interval.upper().toRational()};
        }
        catch (const approximation::PrecisionInsufficient& exception) {
            if (!exception.refinable())
                return std::nullopt;
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

[[nodiscard]] bool strictlyInsideRequestedRange(
    const Expr& point,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto lower = compareExact(
        point, request.lower, builtins, mathematics, angles, assumptions);
    const auto upper = compareExact(
        point, request.upper, builtins, mathematics, angles, assumptions);
    return lower && upper && *lower > 0 && *upper < 0;
}

[[nodiscard]] bool mergeStepDiscontinuity(
    std::vector<StepDiscontinuity>& cuts,
    StepDiscontinuity cut,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    for (auto& existing : cuts) {
        const auto comparison = compareExact(
            existing.position, cut.position,
            builtins, mathematics, angles, assumptions);
        if (!comparison)
            continue;
        if (*comparison == 0) {
            existing.attachLeft = existing.attachLeft && cut.attachLeft;
            existing.attachRight = existing.attachRight && cut.attachRight;
            if (cut.confidence == PlotLandmarkConfidence::Candidate)
                existing.confidence = PlotLandmarkConfidence::Candidate;
            return true;
        }
    }
    if (cuts.size() >= maxSymbolicLandmarks)
        return false;
    cuts.push_back(std::move(cut));
    return true;
}


[[nodiscard]] std::optional<std::vector<Expr>> finitePolynomialLevelRoots(
    const Expr& expression,
    const Expr& level,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    bool includeRequestEndpoints = false) {
    Expr relation = Expr::call(
        builtins.symbol(BuiltinId::Equal), {expression, level});
    relation = solver::simplifyForSolve(
        std::move(relation), builtins, mathematics, angles, assumptions);

    solver::SolveConstraints realConstraint;
    realConstraint.domain = mathematics::NumericDomain::Real;
    solver::SolutionSet solutions = solver::applySolveConstraints(
        solver::solveUnivariatePolynomialRelation(
            relation, request.variable, builtins, mathematics, angles),
        realConstraint, builtins, mathematics, angles);
    if (solutions.kind() == solver::SolutionSetKind::Unresolved) {
        if (auto algebraic = solver::solveRealAlgebraicPolynomialEquation(
                relation, request.variable, builtins, mathematics, angles))
            solutions = std::move(*algebraic);
    }
    if (solutions.kind() == solver::SolutionSetKind::Empty)
        return std::vector<Expr>{};
    if (solutions.kind() != solver::SolutionSetKind::Finite
        || !solutions.conditions().empty()
        || solutions.branches().size() > maxSymbolicLandmarks)
        return std::nullopt;

    std::vector<Expr> roots;
    for (const auto& branch : solutions.branches()) {
        if (!branch.unconditional() || !branch.freeVariables.empty()
            || branch.bindings.size() != 1
            || !branch.bindings.front().variable.sameIdentity(request.variable))
            return std::nullopt;
        Expr root = branch.bindings.front().value;
        const bool inRange = includeRequestEndpoints
            ? insideRequestedRange(
                root, request, builtins, mathematics, angles, assumptions)
            : strictlyInsideRequestedRange(
                root, request, builtins, mathematics, angles, assumptions);
        if (!inRange)
            continue;
        bool duplicate = false;
        for (const Expr& existing : roots) {
            const auto comparison = compareExact(
                existing, root, builtins, mathematics, angles, assumptions);
            if (comparison && *comparison == 0) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate)
            roots.push_back(std::move(root));
    }
    return roots;
}

[[nodiscard]] std::optional<int> derivativeDirectionAt(
    const Expr& derivative,
    const Expr& point,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr value = symbolic::substituteSymbol(derivative, request.variable, point);
    value = solver::simplifyForSolve(
        std::move(value), builtins, mathematics, angles, assumptions);
    if (const auto exact = expression::exact::realRational(value)) {
        if (exact->isZero()) return 0;
        return exact->numerator().isNegative() ? -1 : 1;
    }
    const auto bounds = certifiedConstantBounds(value, builtins, mathematics, angles);
    if (!bounds)
        return std::nullopt;
    const Rational zero{BigInt{0}};
    if (bounds->first > zero)
        return 1;
    if (bounds->second < zero)
        return -1;
    return 0;
}

[[nodiscard]] std::optional<bool> collectAffinePeriodicStepDiscontinuities(
    BuiltinId stepFunction,
    const Expr& argument,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::vector<StepDiscontinuity>& cuts,
    PlotLandmarkConfidence confidence) {
    if (!argument.isCall() || argument.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* trig = builtins.find(argument.asCall().head);
    if (!trig || (trig->id != BuiltinId::Sin && trig->id != BuiltinId::Cos))
        return std::nullopt;

    const Expr& phase = argument.asCall().arguments.front();
    const auto affine = symbolic::toExpressionPolynomial(
        phase, request.variable, builtins, mathematics, angles);
    if (!affine || affine->degree() != 1)
        return false;
    const auto slope = expression::exact::realRational(affine->coefficient(1));
    if (!slope || slope->isZero())
        return false;

    auto endpointAngle = [&](const Expr& endpoint) {
        Expr value = symbolic::substituteSymbol(phase, request.variable, endpoint);
        value = solver::simplifyForSolve(
            std::move(value), builtins, mathematics, angles, assumptions);
        return mathematics::extractExactAngle(value, builtins, mathematics, angles);
    };
    const auto firstAngle = endpointAngle(request.lower);
    const auto secondAngle = endpointAngle(request.upper);
    if (!firstAngle || !secondAngle)
        return false;

    const Rational lowerTurns = std::min(firstAngle->turns, secondAngle->turns);
    const Rational upperTurns = std::max(firstAngle->turns, secondAngle->turns);
    const Expr derivative = symbolic::differentiateExpression(
        argument, request.variable, builtins, mathematics, angles);

    auto addSequence = [&](Rational baseTurns, Rational spacingTurns, bool extremum) {
        const BigInt firstK = ceilRational((lowerTurns - baseTurns) / spacingTurns);
        const BigInt lastK = floorRational((upperTurns - baseTurns) / spacingTurns);
        for (BigInt k = firstK; k <= lastK; k += BigInt{1}) {
            if (cuts.size() >= maxSymbolicLandmarks)
                return false;
            const Rational turns = baseTurns + spacingTurns * Rational{k};
            const Rational deltaTurns = turns - firstAngle->turns;
            Expr delta = angleFromTurns(deltaTurns, builtins, mathematics, angles);
            delta = mathematics::scaleExactExpression(
                Rational{BigInt{1}} / *slope, std::move(delta), builtins);
            Expr position = Expr::call(
                builtins.symbol(BuiltinId::Add),
                std::vector<Expr>{request.lower, std::move(delta)});
            position = solver::simplifyForSolve(
                std::move(position), builtins, mathematics, angles, assumptions);
            if (!insideRequestedRange(
                    position, request, builtins, mathematics, angles, assumptions))
                continue;

            bool attachLeft = false;
            bool attachRight = false;
            if (stepFunction != BuiltinId::Sign && !extremum) {
                const auto direction = derivativeDirectionAt(
                    derivative, position, request,
                    builtins, mathematics, angles, assumptions);
                if (!direction || *direction == 0)
                    return false;
                const bool increasing = *direction > 0;
                if (stepFunction == BuiltinId::Floor) {
                    attachLeft = !increasing;
                    attachRight = increasing;
                }
                else {
                    attachLeft = increasing;
                    attachRight = !increasing;
                }
            }

            if (!mergeStepDiscontinuity(
                    cuts, StepDiscontinuity{
                        std::move(position), attachLeft, attachRight, confidence},
                    builtins, mathematics, angles, assumptions))
                return false;
        }
        return true;
    };

    // sin/cosの値域は[-1,1]。floor/ceilで本当に不連続になる整数levelだけを
    // 列挙し，極値levelは周囲と異なる一点値なのでsingletonとして分離する。
    if (trig->id == BuiltinId::Sin) {
        if (stepFunction == BuiltinId::Floor) {
            if (!addSequence(Rational{}, Rational{BigInt{1}, BigInt{2}}, false)) return false;
            if (!addSequence(Rational{BigInt{1}, BigInt{4}}, Rational{BigInt{1}}, true)) return false;
        }
        else if (stepFunction == BuiltinId::Ceil) {
            if (!addSequence(Rational{BigInt{3}, BigInt{4}}, Rational{BigInt{1}}, true)) return false;
            if (!addSequence(Rational{}, Rational{BigInt{1}, BigInt{2}}, false)) return false;
        }
        else if (!addSequence(Rational{}, Rational{BigInt{1}, BigInt{2}}, false))
            return false;
    }
    else {
        if (stepFunction == BuiltinId::Floor) {
            if (!addSequence(Rational{BigInt{1}, BigInt{4}}, Rational{BigInt{1}, BigInt{2}}, false)) return false;
            if (!addSequence(Rational{}, Rational{BigInt{1}}, true)) return false;
        }
        else if (stepFunction == BuiltinId::Ceil) {
            if (!addSequence(Rational{BigInt{1}, BigInt{2}}, Rational{BigInt{1}}, true)) return false;
            if (!addSequence(Rational{BigInt{1}, BigInt{4}}, Rational{BigInt{1}, BigInt{2}}, false)) return false;
        }
        else if (!addSequence(
                Rational{BigInt{1}, BigInt{4}}, Rational{BigInt{1}, BigInt{2}}, false))
            return false;
    }
    return true;
}

[[nodiscard]] bool collectStepDiscontinuities(
    const Expr& expression,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::vector<StepDiscontinuity>& cuts,
    bool topLevel = true) {
    if (!expression.isCall())
        return true;

    for (const Expr& argument : expression.asCall().arguments)
        if (!collectStepDiscontinuities(
                argument, request, builtins, mathematics, angles,
                assumptions, cuts, false))
            return false;

    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || expression.asCall().arguments.size() != 1)
        return true;
    if (definition->id != BuiltinId::Floor
        && definition->id != BuiltinId::Ceil
        && definition->id != BuiltinId::Trunc
        && definition->id != BuiltinId::Round
        && definition->id != BuiltinId::Frac
        && definition->id != BuiltinId::Sign)
        return true;

    const Expr& argument = expression.asCall().arguments.front();
    const auto polynomial = symbolic::toExpressionPolynomial(
        argument, request.variable, builtins, mathematics, angles);
    if (!polynomial) {
        const auto periodic = collectAffinePeriodicStepDiscontinuities(
            definition->id, argument, request,
            builtins, mathematics, angles, assumptions, cuts,
            topLevel ? PlotLandmarkConfidence::Proven
                     : PlotLandmarkConfidence::Candidate);
        return periodic ? *periodic : false;
    }
    if (polynomial->degree() == 0)
        return true;

    const PlotLandmarkConfidence confidence = topLevel
        ? PlotLandmarkConfidence::Proven
        : PlotLandmarkConfidence::Candidate;

    // affineはlevel preimageが一意なのでsolverを通さずexactに戻す。
    if (polynomial->degree() == 1) {
        const auto slope = expression::exact::realRational(polynomial->coefficient(1));
        const auto intercept = expression::exact::realRational(polynomial->coefficient(0));
        if (!slope || slope->isZero() || !intercept)
            return false;
        const bool increasing = slope->numerator().isPositive();

        auto addAtLevel = [&](const Rational& level, bool attachLeft, bool attachRight) {
            const Rational position = (level - *intercept) / *slope;
            Expr point{Number{position}};
            if (!insideRequestedRange(
                    point, request, builtins, mathematics, angles, assumptions))
                return true;
            return mergeStepDiscontinuity(
                cuts, StepDiscontinuity{
                    std::move(point), attachLeft, attachRight, confidence},
                builtins, mathematics, angles, assumptions);
        };

        if (definition->id == BuiltinId::Sign)
            return addAtLevel(Rational{}, false, false);

        Expr lowerValue = symbolic::substituteSymbol(argument, request.variable, request.lower);
        Expr upperValue = symbolic::substituteSymbol(argument, request.variable, request.upper);
        lowerValue = solver::simplifyForSolve(
            std::move(lowerValue), builtins, mathematics, angles, assumptions);
        upperValue = solver::simplifyForSolve(
            std::move(upperValue), builtins, mathematics, angles, assumptions);
        const auto lowerBounds = certifiedConstantBounds(
            lowerValue, builtins, mathematics, angles);
        const auto upperBounds = certifiedConstantBounds(
            upperValue, builtins, mathematics, angles);
        if (!lowerBounds || !upperBounds)
            return false;

        const Rational rangeLower = std::min(lowerBounds->first, upperBounds->first);
        const Rational rangeUpper = std::max(lowerBounds->second, upperBounds->second);
        auto addIntegerBoundary = [&](const BigInt& level) {
            bool attachArgumentLow = false;
            bool attachArgumentHigh = false;
            switch (definition->id) {
            case BuiltinId::Floor:
            case BuiltinId::Frac:
                attachArgumentHigh = true;
                break;
            case BuiltinId::Ceil:
                attachArgumentLow = true;
                break;
            case BuiltinId::Trunc:
                if (level.isZero())
                    return true; // trunc is continuous at zero.
                if (level.isNegative())
                    attachArgumentLow = true; // negative side behaves like ceil.
                else
                    attachArgumentHigh = true; // positive side behaves like floor.
                break;
            default:
                return false;
            }
            const bool attachLeft = increasing ? attachArgumentLow : attachArgumentHigh;
            const bool attachRight = increasing ? attachArgumentHigh : attachArgumentLow;
            return addAtLevel(Rational{level}, attachLeft, attachRight);
        };

        if (definition->id == BuiltinId::Round) {
            const Rational half{BigInt{1}, BigInt{2}};
            BigInt first = ceilRational(rangeLower - half);
            BigInt last = floorRational(rangeUpper - half);
            if (last < first)
                return true;
            const BigInt span = last - first;
            if (span.bitLength() > 8)
                return false;
            std::size_t visited = 0;
            for (BigInt lowerInteger = first;
                 lowerInteger <= last; lowerInteger += BigInt{1}) {
                if (++visited > maxSymbolicLandmarks + 2)
                    return false;
                // ties-to-even: at n+1/2 the point belongs to the side whose
                // integer value is even.
                const bool lowerEven = (lowerInteger % BigInt{2}).isZero();
                const bool attachArgumentLow = lowerEven;
                const bool attachArgumentHigh = !lowerEven;
                const bool attachLeft = increasing ? attachArgumentLow : attachArgumentHigh;
                const bool attachRight = increasing ? attachArgumentHigh : attachArgumentLow;
                if (!addAtLevel(
                        Rational{lowerInteger} + half, attachLeft, attachRight))
                    return false;
            }
            return true;
        }

        BigInt first = floorRational(rangeLower);
        BigInt last = ceilRational(rangeUpper);
        if (last < first)
            std::swap(first, last);

        const BigInt span = last - first;
        if (span.bitLength() > 8)
            return false;

        std::size_t visited = 0;
        for (BigInt level = first; level <= last; level += BigInt{1}) {
            if (++visited > maxSymbolicLandmarks + 2)
                return false;
            if (!addIntegerBoundary(level))
                return false;
        }
        return true;
    }

    // 非線形多項式は臨界点でrequest区間を分ければ，各piece上で単調になる。
    // endpointと臨界点の函数値から整数levelの有限範囲をcertifiedに求め，
    // 各levelのexact real rootだけをjump候補として列挙する。
    const Expr derivative = symbolic::differentiateExpression(
        argument, request.variable, builtins, mathematics, angles);
    const auto critical = finitePolynomialLevelRoots(
        derivative, integerExpr(0), request,
        builtins, mathematics, angles, assumptions);
    if (!critical)
        return false;

    if (definition->id == BuiltinId::Sign) {
        const auto zeros = finitePolynomialLevelRoots(
            argument, integerExpr(0), request,
            builtins, mathematics, angles, assumptions, true);
        if (!zeros)
            return false;
        for (const Expr& root : *zeros) {
            if (!mergeStepDiscontinuity(
                    cuts, StepDiscontinuity{
                        root, false, false, confidence},
                    builtins, mathematics, angles, assumptions))
                return false;
        }
        return true;
    }

    std::vector<Expr> rangePoints{request.lower, request.upper};
    rangePoints.insert(rangePoints.end(), critical->begin(), critical->end());
    std::optional<Rational> rangeLower;
    std::optional<Rational> rangeUpper;
    for (const Expr& point : rangePoints) {
        Expr value = symbolic::substituteSymbol(argument, request.variable, point);
        value = solver::simplifyForSolve(
            std::move(value), builtins, mathematics, angles, assumptions);
        const auto bounds = certifiedConstantBounds(value, builtins, mathematics, angles);
        if (!bounds)
            return false;
        if (!rangeLower || bounds->first < *rangeLower)
            rangeLower = bounds->first;
        if (!rangeUpper || bounds->second > *rangeUpper)
            rangeUpper = bounds->second;
    }
    if (!rangeLower || !rangeUpper)
        return false;

    const bool roundFunction = definition->id == BuiltinId::Round;
    const Rational half{BigInt{1}, BigInt{2}};
    BigInt first = roundFunction
        ? ceilRational(*rangeLower - half)
        : floorRational(*rangeLower);
    BigInt last = roundFunction
        ? floorRational(*rangeUpper - half)
        : ceilRational(*rangeUpper);
    const BigInt span = last - first;
    if (span.isNegative())
        return true;
    if (span.bitLength() > 8)
        return false;

    std::size_t visitedLevels = 0;
    for (BigInt levelIndex = first; levelIndex <= last; levelIndex += BigInt{1}) {
        if (++visitedLevels > maxSymbolicLandmarks + 2)
            return false;
        if (definition->id == BuiltinId::Trunc && levelIndex.isZero())
            continue;
        const Rational level = roundFunction
            ? Rational{levelIndex} + half
            : Rational{levelIndex};
        const auto roots = finitePolynomialLevelRoots(
            argument, Expr{Number{level}}, request,
            builtins, mathematics, angles, assumptions, true);
        if (!roots)
            return false;
        for (const Expr& root : *roots) {
            const auto direction = derivativeDirectionAt(
                derivative, root, request,
                builtins, mathematics, angles, assumptions);
            bool attachLeft = false;
            bool attachRight = false;
            if (direction && *direction != 0) {
                const bool increasing = *direction > 0;
                bool attachArgumentLow = false;
                bool attachArgumentHigh = false;
                switch (definition->id) {
                case BuiltinId::Floor:
                case BuiltinId::Frac:
                    attachArgumentHigh = true;
                    break;
                case BuiltinId::Ceil:
                    attachArgumentLow = true;
                    break;
                case BuiltinId::Trunc:
                    if (levelIndex.isNegative())
                        attachArgumentLow = true;
                    else
                        attachArgumentHigh = true;
                    break;
                case BuiltinId::Round: {
                    const bool lowerEven = (levelIndex % BigInt{2}).isZero();
                    attachArgumentLow = lowerEven;
                    attachArgumentHigh = !lowerEven;
                    break;
                }
                default:
                    break;
                }
                attachLeft = increasing ? attachArgumentLow : attachArgumentHigh;
                attachRight = increasing ? attachArgumentHigh : attachArgumentLow;
            }
            // derivative==0ではlevelに接するだけでjumpが無い場合もある。
            // その判定を推測せずsingletonへ分離すれば，偽の接続だけは作らない。
            if (!mergeStepDiscontinuity(
                    cuts, StepDiscontinuity{
                        root, attachLeft, attachRight, confidence},
                    builtins, mathematics, angles, assumptions))
                return false;
        }
    }
    return true;
}

[[nodiscard]] bool splitDomainAtStepDiscontinuities(
    PlotAnalysis& analysis,
    std::vector<StepDiscontinuity> cuts,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (cuts.empty())
        return true;

    std::stable_sort(cuts.begin(), cuts.end(), [&](const auto& lhs, const auto& rhs) {
        const auto comparison = compareExact(
            lhs.position, rhs.position, builtins, mathematics, angles, assumptions);
        return comparison && *comparison < 0;
    });
    for (std::size_t i = 1; i < cuts.size(); ++i) {
        const auto comparison = compareExact(
            cuts[i - 1].position, cuts[i].position,
            builtins, mathematics, angles, assumptions);
        if (!comparison || *comparison >= 0)
            return false;
    }

    std::vector<PlotInterval> split;
    for (const PlotInterval& source : analysis.domain.intervals) {
        const auto sourceSpan = compareExact(
            source.lower, source.upper, builtins, mathematics, angles, assumptions);
        if (!sourceSpan)
            return false;
        if (*sourceSpan == 0) {
            split.push_back(source);
            continue;
        }

        Expr currentLower = source.lower;
        PlotEndpointInclusion currentLowerInclusion = source.lowerInclusion;
        PlotEndpointInclusion finalUpperInclusion = source.upperInclusion;
        bool appendUpperSingleton = false;

        for (const auto& cut : cuts) {
            const auto toLower = compareExact(
                cut.position, source.lower, builtins, mathematics, angles, assumptions);
            const auto toUpper = compareExact(
                cut.position, source.upper, builtins, mathematics, angles, assumptions);
            if (!toLower || !toUpper)
                return false;
            if (*toLower < 0 || *toUpper > 0)
                continue;

            // request/domain端そのものがstep jumpなら，単なるPlot終端とは区別する。
            // visible側と函数値が一致しない場合だけsingletonへ分離し，斜め/縦の偽線を防ぐ。
            if (*toLower == 0) {
                if (source.lowerInclusion == PlotEndpointInclusion::Closed && !cut.attachRight) {
                    split.push_back(PlotInterval{
                        cut.position, PlotEndpointInclusion::Closed,
                        cut.position, PlotEndpointInclusion::Closed});
                    currentLower = cut.position;
                    currentLowerInclusion = PlotEndpointInclusion::Open;
                }
                if (source.lowerInclusion == PlotEndpointInclusion::Closed)
                    addLandmark(
                        analysis, cut.position,
                        PlotLandmarkKind::JumpDiscontinuity, cut.confidence);
                continue;
            }
            if (*toUpper == 0) {
                if (source.upperInclusion == PlotEndpointInclusion::Closed && !cut.attachLeft) {
                    finalUpperInclusion = PlotEndpointInclusion::Open;
                    appendUpperSingleton = true;
                }
                if (source.upperInclusion == PlotEndpointInclusion::Closed)
                    addLandmark(
                        analysis, cut.position,
                        PlotLandmarkKind::JumpDiscontinuity, cut.confidence);
                continue;
            }

            split.push_back(PlotInterval{
                currentLower,
                currentLowerInclusion,
                cut.position,
                cut.attachLeft ? PlotEndpointInclusion::Closed : PlotEndpointInclusion::Open});
            if (!cut.attachLeft && !cut.attachRight) {
                split.push_back(PlotInterval{
                    cut.position, PlotEndpointInclusion::Closed,
                    cut.position, PlotEndpointInclusion::Closed});
            }
            currentLower = cut.position;
            currentLowerInclusion = cut.attachRight
                ? PlotEndpointInclusion::Closed
                : PlotEndpointInclusion::Open;

            addLandmark(
                analysis, cut.position,
                PlotLandmarkKind::JumpDiscontinuity, cut.confidence);
        }
        split.push_back(PlotInterval{
            std::move(currentLower), currentLowerInclusion,
            source.upper, finalUpperInclusion});
        if (appendUpperSingleton) {
            split.push_back(PlotInterval{
                source.upper, PlotEndpointInclusion::Closed,
                source.upper, PlotEndpointInclusion::Closed});
        }
    }
    if (split.size() > maxSymbolicLandmarks * 3 + analysis.domain.intervals.size())
        return false;
    analysis.domain.intervals = std::move(split);
    return true;
}


enum class PeriodicRealConstraint {
    NonNegative,
    Positive,
    NonZero,
    Defined,
    UnitClosed,
    UnitOpen
};

// 外側函数が周期函数のzero/poleへ近づいたときの実軸上の極限挙動。
// domain制約とは分離して保持し，例えばtan[x]^(-1/2)ではtanのpoleが
// 外側Powerによって0へ写ることを，元のtanのpoleと取り違えない。
enum class PeriodicBoundaryLimitBehavior {
    None,
    FiniteZero,
    Diverges
};

enum class PeriodicBoundaryKind {
    Zero,
    Pole,
    UnitMagnitude
};

struct PeriodicConstraintBoundary final {
    Rational turns;
    Expr position;
    PeriodicBoundaryKind kind = PeriodicBoundaryKind::Zero;
};

[[nodiscard]] Rational normalizedUnitTurn(const Rational& turns) {
    return turns - Rational{floorRational(turns)};
}

[[nodiscard]] int sinSignAtTurn(const Rational& turns) {
    const Rational normalized = normalizedUnitTurn(turns);
    const Rational half{BigInt{1}, BigInt{2}};
    if (normalized.isZero() || normalized == half)
        return 0;
    return normalized < half ? 1 : -1;
}

[[nodiscard]] int cosSignAtTurn(const Rational& turns) {
    const Rational normalized = normalizedUnitTurn(turns);
    const Rational quarter{BigInt{1}, BigInt{4}};
    const Rational threeQuarters{BigInt{3}, BigInt{4}};
    if (normalized == quarter || normalized == threeQuarters)
        return 0;
    return normalized < quarter || normalized > threeQuarters ? 1 : -1;
}

[[nodiscard]] int periodicTrigSignAtTurn(BuiltinId function, const Rational& turns) {
    const int sinSign = sinSignAtTurn(turns);
    const int cosSign = cosSignAtTurn(turns);
    switch (function) {
    case BuiltinId::Sin: return sinSign;
    case BuiltinId::Cos: return cosSign;
    case BuiltinId::Tan:
    case BuiltinId::Cot: return sinSign * cosSign;
    case BuiltinId::Sec: return cosSign;
    case BuiltinId::Csc: return sinSign;
    default: return 0;
    }
}

[[nodiscard]] bool periodicAbsLessThanOneAtTurn(
    BuiltinId function,
    const Rational& turns) {
    const Rational half{BigInt{1}, BigInt{2}};
    const Rational eighth{BigInt{1}, BigInt{8}};
    Rational withinHalf = turns - Rational{floorRational(turns / half)} * half;
    if (function == BuiltinId::Tan)
        return withinHalf < eighth || withinHalf > half - eighth;
    if (function == BuiltinId::Cot) {
        const Rational distance = withinHalf >= Rational{BigInt{1}, BigInt{4}}
            ? withinHalf - Rational{BigInt{1}, BigInt{4}}
            : Rational{BigInt{1}, BigInt{4}} - withinHalf;
        return distance < eighth;
    }
    if (function == BuiltinId::Sin || function == BuiltinId::Cos)
        return true;
    return false;
}

[[nodiscard]] bool constraintAcceptsMidpoint(
    PeriodicRealConstraint constraint,
    BuiltinId function,
    const Rational& turns) {
    const int sign = periodicTrigSignAtTurn(function, turns);
    switch (constraint) {
    case PeriodicRealConstraint::NonNegative: return sign >= 0;
    case PeriodicRealConstraint::Positive: return sign > 0;
    case PeriodicRealConstraint::NonZero: return sign != 0;
    case PeriodicRealConstraint::Defined: return true;
    case PeriodicRealConstraint::UnitClosed:
    case PeriodicRealConstraint::UnitOpen:
        return periodicAbsLessThanOneAtTurn(function, turns);
    }
    return false;
}

[[nodiscard]] PlotEndpointInclusion boundaryInclusion(
    PeriodicRealConstraint constraint,
    PeriodicBoundaryKind kind) noexcept {
    if (kind == PeriodicBoundaryKind::Zero
        && constraint == PeriodicRealConstraint::NonNegative)
        return PlotEndpointInclusion::Closed;
    if (kind == PeriodicBoundaryKind::UnitMagnitude
        && constraint == PeriodicRealConstraint::UnitClosed)
        return PlotEndpointInclusion::Closed;
    return PlotEndpointInclusion::Open;
}

// 実Plotでは，周期函数の複素continuationではなく実軸上の符号区間を直接使う。
// sqrt/log/reciprocalのdomain制約をturn空間でexactに解き，pole/zeroを跨ぐsamplingを防ぐ。
[[nodiscard]] std::optional<PlotAnalysis> analyzeAffinePeriodicRealConstraint(
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const Expr* periodic = nullptr;
    PeriodicRealConstraint constraint = PeriodicRealConstraint::NonNegative;
    PeriodicBoundaryLimitBehavior zeroBehavior = PeriodicBoundaryLimitBehavior::None;
    PeriodicBoundaryLimitBehavior poleBehavior = PeriodicBoundaryLimitBehavior::None;
    bool reciprocal = false;

    if (request.expression.isCall() && request.expression.asCall().arguments.size() == 1) {
        const auto* outer = builtins.find(request.expression.asCall().head);
        if (outer && outer->id == BuiltinId::Sqrt) {
            periodic = &request.expression.asCall().arguments.front();
            constraint = PeriodicRealConstraint::NonNegative;
            zeroBehavior = PeriodicBoundaryLimitBehavior::FiniteZero;
            poleBehavior = PeriodicBoundaryLimitBehavior::Diverges;
        }
        else if (outer && (outer->id == BuiltinId::Log
                || outer->id == BuiltinId::Log2
                || outer->id == BuiltinId::Log10)) {
            periodic = &request.expression.asCall().arguments.front();
            constraint = PeriodicRealConstraint::Positive;
            zeroBehavior = PeriodicBoundaryLimitBehavior::Diverges;
            poleBehavior = PeriodicBoundaryLimitBehavior::Diverges;
        }
        else if (outer && (outer->id == BuiltinId::Asin || outer->id == BuiltinId::Acos)) {
            periodic = &request.expression.asCall().arguments.front();
            constraint = PeriodicRealConstraint::UnitClosed;
        }
        else if (outer && outer->id == BuiltinId::Atanh) {
            periodic = &request.expression.asCall().arguments.front();
            constraint = PeriodicRealConstraint::UnitOpen;
        }
    }
    else if (builtins.isCallTo(request.expression, BuiltinId::Power)
        && request.expression.asCall().arguments.size() == 2) {
        const auto& arguments = request.expression.asCall().arguments;
        const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
        const auto exponentFacts = knowledge.facts(arguments[1]);
        const bool exponentNonInteger = knowledge.prove(mathematics::elementOf(
            arguments[1], mathematics::NumericDomain::Integer)) == TruthValue::False;
        // principal Powerの実非整数指数は，周期baseの符号区間をsqrt/logと同じ
        // turn空間でexactに切れる。指数が整数になり得る場合は負側の離散実点を
        // 落とし得るため，この特殊化を使わない。
        if (exponentFacts.isProvablyReal() && exponentNonInteger) {
            if (exponentFacts.sign == mathematics::RealSign::Positive
                || exponentFacts.sign == mathematics::RealSign::NonNegative) {
                periodic = &arguments[0];
                constraint = PeriodicRealConstraint::NonNegative;
                zeroBehavior = PeriodicBoundaryLimitBehavior::FiniteZero;
                poleBehavior = PeriodicBoundaryLimitBehavior::Diverges;
            }
            else if (exponentFacts.sign == mathematics::RealSign::Negative
                || exponentFacts.sign == mathematics::RealSign::NonPositive) {
                periodic = &arguments[0];
                constraint = PeriodicRealConstraint::Positive;
                zeroBehavior = PeriodicBoundaryLimitBehavior::Diverges;
                poleBehavior = PeriodicBoundaryLimitBehavior::FiniteZero;
            }
        }
    }
    else if (builtins.isCallTo(request.expression, BuiltinId::Divide)
        && request.expression.asCall().arguments.size() == 2) {
        const auto numerator = expression::exact::realRational(
            request.expression.asCall().arguments[0]);
        if (numerator && !numerator->isZero()) {
            periodic = &request.expression.asCall().arguments[1];
            constraint = PeriodicRealConstraint::NonZero;
            zeroBehavior = PeriodicBoundaryLimitBehavior::Diverges;
            poleBehavior = PeriodicBoundaryLimitBehavior::FiniteZero;
            reciprocal = true;
        }
    }
    if (!periodic)
        return std::nullopt;

    if ((constraint == PeriodicRealConstraint::NonNegative
            || constraint == PeriodicRealConstraint::Positive)
        && builtins.isCallTo(*periodic, BuiltinId::Abs)
        && periodic->asCall().arguments.size() == 1) {
        periodic = &periodic->asCall().arguments.front();
        constraint = constraint == PeriodicRealConstraint::NonNegative
            ? PeriodicRealConstraint::Defined
            : PeriodicRealConstraint::NonZero;
    }
    else if ((constraint == PeriodicRealConstraint::NonNegative
            || constraint == PeriodicRealConstraint::Positive)
        && builtins.isCallTo(*periodic, BuiltinId::Power)
        && periodic->asCall().arguments.size() == 2) {
        const Expr& exponent = periodic->asCall().arguments[1];
        if (exponent.isNumber() && exponent.asNumber().isReal()
            && exponent.asNumber().asReal().isInteger()) {
            const BigInt& power = exponent.asNumber().asReal().asInteger();
            if (power.isPositive() && (power % BigInt{2}).isZero()) {
                periodic = &periodic->asCall().arguments[0];
                constraint = constraint == PeriodicRealConstraint::NonNegative
                    ? PeriodicRealConstraint::Defined
                    : PeriodicRealConstraint::NonZero;
            }
        }
    }

    if (!periodic->isCall() || periodic->asCall().arguments.size() != 1)
        return std::nullopt;

    const auto* trig = builtins.find(periodic->asCall().head);
    if (!trig || (trig->id != BuiltinId::Sin && trig->id != BuiltinId::Cos
        && trig->id != BuiltinId::Tan && trig->id != BuiltinId::Cot
        && trig->id != BuiltinId::Sec && trig->id != BuiltinId::Csc))
        return std::nullopt;

    const Expr& phase = periodic->asCall().arguments.front();
    const auto affine = symbolic::toExpressionPolynomial(
        phase, request.variable, builtins, mathematics, angles);
    if (!affine || affine->degree() != 1)
        return std::nullopt;
    const auto slope = expression::exact::realRational(affine->coefficient(1));
    if (!slope || slope->isZero())
        return std::nullopt;

    auto endpointAngle = [&](const Expr& endpoint) {
        Expr value = symbolic::substituteSymbol(phase, request.variable, endpoint);
        value = solver::simplifyForSolve(
            std::move(value), builtins, mathematics, angles, assumptions);
        return mathematics::extractExactAngle(value, builtins, mathematics, angles);
    };
    const auto firstAngle = endpointAngle(request.lower);
    const auto secondAngle = endpointAngle(request.upper);
    if (!firstAngle || !secondAngle)
        return std::nullopt;

    if ((constraint == PeriodicRealConstraint::Defined
            || constraint == PeriodicRealConstraint::UnitClosed)
        && (trig->id == BuiltinId::Sin || trig->id == BuiltinId::Cos)) {
        PlotAnalysis result;
        result.domain.coverage = PlotDomainCoverage::Complete;
        result.domain.intervals.push_back(PlotInterval{
            request.lower, PlotEndpointInclusion::Closed,
            request.upper, PlotEndpointInclusion::Closed});
        return result;
    }

    auto positionFromTurn = [&](const Rational& turns) {
        const Rational deltaTurns = turns - firstAngle->turns;
        Expr delta = angleFromTurns(deltaTurns, builtins, mathematics, angles);
        delta = mathematics::scaleExactExpression(
            Rational{BigInt{1}} / *slope, std::move(delta), builtins);
        Expr position = Expr::call(
            builtins.symbol(BuiltinId::Add),
            std::vector<Expr>{request.lower, std::move(delta)});
        return solver::simplifyForSolve(
            std::move(position), builtins, mathematics, angles, assumptions);
    };

    const Rational lowerTurns = std::min(firstAngle->turns, secondAngle->turns);
    const Rational upperTurns = std::max(firstAngle->turns, secondAngle->turns);
    std::vector<PeriodicConstraintBoundary> boundaries;

    auto addSequence = [&](Rational base, Rational spacing, PeriodicBoundaryKind kind) {
        const BigInt firstK = ceilRational((lowerTurns - base) / spacing);
        const BigInt lastK = floorRational((upperTurns - base) / spacing);
        for (BigInt k = firstK; k <= lastK; k += BigInt{1}) {
            if (boundaries.size() >= maxSymbolicLandmarks)
                return false;
            const Rational turns = base + spacing * Rational{k};
            boundaries.push_back(PeriodicConstraintBoundary{
                turns, positionFromTurn(turns), kind});
        }
        return true;
    };

    const Rational half{BigInt{1}, BigInt{2}};
    const Rational quarter{BigInt{1}, BigInt{4}};
    if (constraint == PeriodicRealConstraint::Defined) {
        if ((trig->id == BuiltinId::Tan || trig->id == BuiltinId::Sec)
            && !addSequence(quarter, half, PeriodicBoundaryKind::Pole))
            return std::nullopt;
        if ((trig->id == BuiltinId::Cot || trig->id == BuiltinId::Csc)
            && !addSequence(Rational{}, half, PeriodicBoundaryKind::Pole))
            return std::nullopt;
    }
    else if (constraint == PeriodicRealConstraint::UnitClosed
        || constraint == PeriodicRealConstraint::UnitOpen) {
        switch (trig->id) {
        case BuiltinId::Sin:
            if (!addSequence(quarter, half, PeriodicBoundaryKind::UnitMagnitude)) return std::nullopt;
            break;
        case BuiltinId::Cos:
            if (!addSequence(Rational{}, half, PeriodicBoundaryKind::UnitMagnitude)) return std::nullopt;
            break;
        case BuiltinId::Tan:
        case BuiltinId::Cot:
            if (!addSequence(
                    Rational{BigInt{1}, BigInt{8}},
                    Rational{BigInt{1}, BigInt{4}},
                    PeriodicBoundaryKind::UnitMagnitude)) return std::nullopt;
            break;
        case BuiltinId::Sec:
            if (!addSequence(Rational{}, half, PeriodicBoundaryKind::UnitMagnitude)) return std::nullopt;
            break;
        case BuiltinId::Csc:
            if (!addSequence(quarter, half, PeriodicBoundaryKind::UnitMagnitude)) return std::nullopt;
            break;
        default:
            return std::nullopt;
        }
    }
    else {
        switch (trig->id) {
        case BuiltinId::Sin:
            if (!addSequence(Rational{}, half, PeriodicBoundaryKind::Zero)) return std::nullopt;
            break;
        case BuiltinId::Cos:
            if (!addSequence(quarter, half, PeriodicBoundaryKind::Zero)) return std::nullopt;
            break;
        case BuiltinId::Tan:
            if (!addSequence(Rational{}, half, PeriodicBoundaryKind::Zero)
                || !addSequence(quarter, half, PeriodicBoundaryKind::Pole)) return std::nullopt;
            break;
        case BuiltinId::Cot:
            if (!addSequence(quarter, half, PeriodicBoundaryKind::Zero)
                || !addSequence(Rational{}, half, PeriodicBoundaryKind::Pole)) return std::nullopt;
            break;
        case BuiltinId::Sec:
            if (!addSequence(quarter, half, PeriodicBoundaryKind::Pole)) return std::nullopt;
            break;
        case BuiltinId::Csc:
            if (!addSequence(Rational{}, half, PeriodicBoundaryKind::Pole)) return std::nullopt;
            break;
        default:
            return std::nullopt;
        }
    }

    std::sort(boundaries.begin(), boundaries.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.turns < rhs.turns;
    });
    if (slope->numerator().isNegative())
        std::reverse(boundaries.begin(), boundaries.end());

    struct Node final {
        Rational turns;
        Expr position;
        std::optional<PeriodicBoundaryKind> kind;
    };
    std::vector<Node> nodes;
    nodes.push_back(Node{firstAngle->turns, request.lower, std::nullopt});
    for (const auto& boundary : boundaries) {
        if (boundary.turns == firstAngle->turns) {
            nodes.front().kind = boundary.kind;
            continue;
        }
        if (boundary.turns == secondAngle->turns)
            continue;
        nodes.push_back(Node{boundary.turns, boundary.position, boundary.kind});
    }
    nodes.push_back(Node{secondAngle->turns, request.upper, std::nullopt});
    for (const auto& boundary : boundaries) {
        if (boundary.turns == secondAngle->turns) {
            nodes.back().kind = boundary.kind;
            break;
        }
    }

    PlotAnalysis result;
    result.domain.coverage = PlotDomainCoverage::Complete;
    for (std::size_t i = 0; i + 1 < nodes.size(); ++i) {
        const Rational midpointTurns = (nodes[i].turns + nodes[i + 1].turns)
            / Rational{BigInt{2}};
        if (!constraintAcceptsMidpoint(constraint, trig->id, midpointTurns))
            continue;

        const PlotEndpointInclusion lowerInclusion = nodes[i].kind
            ? boundaryInclusion(constraint, *nodes[i].kind)
            : PlotEndpointInclusion::Closed;
        const PlotEndpointInclusion upperInclusion = nodes[i + 1].kind
            ? boundaryInclusion(constraint, *nodes[i + 1].kind)
            : PlotEndpointInclusion::Closed;
        result.domain.intervals.push_back(PlotInterval{
            nodes[i].position, lowerInclusion,
            nodes[i + 1].position, upperInclusion});
    }

    if (constraint == PeriodicRealConstraint::NonNegative
        || constraint == PeriodicRealConstraint::UnitClosed) {
        // request端点zeroや|f|=1だけがdomainに残る場合はsingletonとして保持する。
        const PeriodicBoundaryKind singletonKind = constraint == PeriodicRealConstraint::NonNegative
            ? PeriodicBoundaryKind::Zero
            : PeriodicBoundaryKind::UnitMagnitude;
        for (const Node& node : nodes) {
            if (!node.kind || *node.kind != singletonKind)
                continue;
            bool included = false;
            for (const PlotInterval& interval : result.domain.intervals) {
                if ((interval.lower == node.position
                        && interval.lowerInclusion == PlotEndpointInclusion::Closed)
                    || (interval.upper == node.position
                        && interval.upperInclusion == PlotEndpointInclusion::Closed)) {
                    included = true;
                    break;
                }
            }
            if (!included)
                result.domain.intervals.push_back(PlotInterval{
                    node.position, PlotEndpointInclusion::Closed,
                    node.position, PlotEndpointInclusion::Closed});
        }
    }

    std::sort(result.domain.intervals.begin(), result.domain.intervals.end(),
        [&](const PlotInterval& lhs, const PlotInterval& rhs) {
            const auto comparison = compareExact(
                lhs.lower, rhs.lower, builtins, mathematics, angles, assumptions);
            return comparison && *comparison < 0;
        });

    for (const auto& boundary : boundaries) {
        addLandmark(result, boundary.position, PlotLandmarkKind::DomainBoundary);
        const bool undefinedBoundary = boundary.kind == PeriodicBoundaryKind::Pole
            || (boundary.kind == PeriodicBoundaryKind::Zero
                && constraint != PeriodicRealConstraint::NonNegative)
            || (boundary.kind == PeriodicBoundaryKind::UnitMagnitude
                && constraint == PeriodicRealConstraint::UnitOpen);
        if (undefinedBoundary)
            addLandmark(result, boundary.position, PlotLandmarkKind::UndefinedPoint);

        if (constraint == PeriodicRealConstraint::NonNegative
            && boundary.kind == PeriodicBoundaryKind::Zero)
            addLandmark(result, boundary.position, PlotLandmarkKind::BranchPoint);
        if (boundary.kind == PeriodicBoundaryKind::Zero
            && zeroBehavior == PeriodicBoundaryLimitBehavior::Diverges) {
            if (reciprocal)
                addLandmark(result, boundary.position, PlotLandmarkKind::Pole);
            addLandmark(result, boundary.position, PlotLandmarkKind::VerticalAsymptote);
        }
        else if (boundary.kind == PeriodicBoundaryKind::Pole
            && poleBehavior == PeriodicBoundaryLimitBehavior::Diverges)
            addLandmark(result, boundary.position, PlotLandmarkKind::VerticalAsymptote);
        else if (boundary.kind == PeriodicBoundaryKind::Pole
            && poleBehavior == PeriodicBoundaryLimitBehavior::FiniteZero)
            addLandmark(
                result, boundary.position, PlotLandmarkKind::RemovableSingularity,
                PlotLandmarkConfidence::Proven, integerExpr(0));
        else if (constraint == PeriodicRealConstraint::UnitOpen
            && boundary.kind == PeriodicBoundaryKind::UnitMagnitude)
            addLandmark(result, boundary.position, PlotLandmarkKind::VerticalAsymptote);
    }
    return result;
}


struct PeriodicUndefinedPoint final {
    Expr position;
};

[[nodiscard]] bool mergePeriodicUndefinedPoint(
    std::vector<PeriodicUndefinedPoint>& points,
    Expr point,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    for (const auto& existing : points) {
        const auto comparison = compareExact(
            existing.position, point, builtins, mathematics, angles, assumptions);
        if (comparison && *comparison == 0)
            return true;
    }
    if (points.size() >= maxSymbolicLandmarks)
        return false;
    points.push_back(PeriodicUndefinedPoint{std::move(point)});
    return true;
}

// tan/cot/sec/cscが算術やreal-safe函数の内側にいても，そのpoleでは元式自体が未定義。
// outer函数の発散性までは推測せず，samplingを跨がせないholeだけを再帰的に伝播する。
[[nodiscard]] bool collectNestedAffinePeriodicUndefinedPoints(
    const Expr& expression,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::vector<PeriodicUndefinedPoint>& points) {
    if (!expression.isCall())
        return true;

    for (const Expr& argument : expression.asCall().arguments)
        if (!collectNestedAffinePeriodicUndefinedPoints(
                argument, request, builtins, mathematics, angles, assumptions, points))
            return false;

    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || expression.asCall().arguments.size() != 1)
        return true;

    Rational poleBaseTurns;
    if (definition->id == BuiltinId::Tan || definition->id == BuiltinId::Sec
        || definition->id == BuiltinId::Tanc)
        poleBaseTurns = Rational{BigInt{1}, BigInt{4}};
    else if (definition->id == BuiltinId::Cot || definition->id == BuiltinId::Csc)
        poleBaseTurns = Rational{};
    else
        return true;

    const Expr& phase = expression.asCall().arguments.front();
    const auto affine = symbolic::toExpressionPolynomial(
        phase, request.variable, builtins, mathematics, angles);
    if (!affine || affine->degree() != 1)
        return false;
    const auto slope = expression::exact::realRational(affine->coefficient(1));
    if (!slope || slope->isZero())
        return false;

    auto endpointAngle = [&](const Expr& endpoint) {
        Expr value = symbolic::substituteSymbol(phase, request.variable, endpoint);
        value = solver::simplifyForSolve(
            std::move(value), builtins, mathematics, angles, assumptions);
        return mathematics::extractExactAngle(value, builtins, mathematics, angles);
    };
    const auto firstAngle = endpointAngle(request.lower);
    const auto secondAngle = endpointAngle(request.upper);
    if (!firstAngle || !secondAngle)
        return false;

    const Rational lowerTurns = std::min(firstAngle->turns, secondAngle->turns);
    const Rational upperTurns = std::max(firstAngle->turns, secondAngle->turns);
    const Rational firstKValue = (lowerTurns - poleBaseTurns) * Rational{BigInt{2}};
    const Rational lastKValue = (upperTurns - poleBaseTurns) * Rational{BigInt{2}};
    const BigInt firstK = ceilRational(firstKValue);
    const BigInt lastK = floorRational(lastKValue);

    for (BigInt k = firstK; k <= lastK; k += BigInt{1}) {
        const Rational turns = poleBaseTurns + Rational{k, BigInt{2}};
        const Rational deltaTurns = turns - firstAngle->turns;
        Expr delta = angleFromTurns(deltaTurns, builtins, mathematics, angles);
        delta = mathematics::scaleExactExpression(
            Rational{BigInt{1}} / *slope, std::move(delta), builtins);
        Expr position = Expr::call(
            builtins.symbol(BuiltinId::Add),
            std::vector<Expr>{request.lower, std::move(delta)});
        position = solver::simplifyForSolve(
            std::move(position), builtins, mathematics, angles, assumptions);
        if (!insideRequestedRange(
                position, request, builtins, mathematics, angles, assumptions))
            continue;
        if (!mergePeriodicUndefinedPoint(
                points, std::move(position), builtins, mathematics, angles, assumptions))
            return false;
    }
    return true;
}

[[nodiscard]] bool splitDomainAtOpenCuts(
    PlotAnalysis& analysis,
    std::vector<PeriodicUndefinedPoint> cuts,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (cuts.empty())
        return true;
    std::stable_sort(cuts.begin(), cuts.end(), [&](const auto& lhs, const auto& rhs) {
        const auto comparison = compareExact(
            lhs.position, rhs.position, builtins, mathematics, angles, assumptions);
        return comparison && *comparison < 0;
    });

    std::vector<PlotInterval> split;
    for (const PlotInterval& source : analysis.domain.intervals) {
        Expr currentLower = source.lower;
        PlotEndpointInclusion currentLowerInclusion = source.lowerInclusion;
        PlotEndpointInclusion finalUpperInclusion = source.upperInclusion;

        for (const auto& cut : cuts) {
            const auto toLower = compareExact(
                cut.position, source.lower, builtins, mathematics, angles, assumptions);
            const auto toUpper = compareExact(
                cut.position, source.upper, builtins, mathematics, angles, assumptions);
            if (!toLower || !toUpper)
                return false;
            if (*toLower < 0 || *toUpper > 0)
                continue;
            if (*toLower == 0) {
                currentLowerInclusion = PlotEndpointInclusion::Open;
                addLandmark(analysis, cut.position, PlotLandmarkKind::DomainBoundary);
                addLandmark(analysis, cut.position, PlotLandmarkKind::UndefinedPoint);
                continue;
            }
            if (*toUpper == 0) {
                finalUpperInclusion = PlotEndpointInclusion::Open;
                addLandmark(analysis, cut.position, PlotLandmarkKind::DomainBoundary);
                addLandmark(analysis, cut.position, PlotLandmarkKind::UndefinedPoint);
                continue;
            }

            split.push_back(PlotInterval{
                currentLower, currentLowerInclusion,
                cut.position, PlotEndpointInclusion::Open});
            currentLower = cut.position;
            currentLowerInclusion = PlotEndpointInclusion::Open;
            addLandmark(analysis, cut.position, PlotLandmarkKind::DomainBoundary);
            addLandmark(analysis, cut.position, PlotLandmarkKind::UndefinedPoint);
        }
        split.push_back(PlotInterval{
            std::move(currentLower), currentLowerInclusion,
            source.upper, finalUpperInclusion});
    }
    if (split.size() > maxSymbolicLandmarks + analysis.domain.intervals.size())
        return false;
    analysis.domain.intervals = std::move(split);
    return true;
}

[[nodiscard]] bool containsPrincipalPowerNeedingBranchProof(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    if (!expression.isCall())
        return false;

    if (builtins.isCallTo(expression, BuiltinId::Power)
        && expression.asCall().arguments.size() == 2) {
        const auto& arguments = expression.asCall().arguments;
        const Expr& exponent = arguments[1];
        const bool exactInteger = exponent.isNumber()
            && exponent.asNumber().isReal()
            && exponent.asNumber().asReal().isInteger();
        if (!exactInteger) {
            const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
            const auto baseFacts = knowledge.facts(arguments[0]);
            // 正の実数底ではprincipal Powerは任意の実指数で実数なので，
            // domain不完全性があっても原因は子式側にあり，既存のcut解析へ任せられる。
            if (baseFacts.sign != mathematics::RealSign::Positive)
                return true;
        }
    }

    for (const Expr& argument : expression.asCall().arguments)
        if (containsPrincipalPowerNeedingBranchProof(
                argument, builtins, mathematics, assumptions))
            return true;
    return false;
}

[[nodiscard]] bool containsRestrictedRealWrapperOverPeriodic(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall())
        return false;

    const auto containsPeriodic = [&](const auto& self, const Expr& node) -> bool {
        if (!node.isCall())
            return false;
        if (const auto* definition = builtins.find(node.asCall().head)) {
            if (definition->id == BuiltinId::Tan || definition->id == BuiltinId::Cot
                || definition->id == BuiltinId::Sec || definition->id == BuiltinId::Csc)
                return true;
        }
        for (const Expr& argument : node.asCall().arguments)
            if (self(self, argument))
                return true;
        return false;
    };

    if (const auto* definition = builtins.find(expression.asCall().head)) {
        switch (definition->id) {
        case BuiltinId::Sqrt:
        case BuiltinId::Log:
        case BuiltinId::Log1p:
        case BuiltinId::Log2:
        case BuiltinId::Log10:
        case BuiltinId::ExponentialIntegralEi:
        case BuiltinId::CosineIntegralCi:
        case BuiltinId::Asin:
        case BuiltinId::Acos:
        case BuiltinId::Acosh:
        case BuiltinId::Atanh:
            if (!expression.asCall().arguments.empty()
                && containsPeriodic(containsPeriodic, expression.asCall().arguments.front()))
                return true;
            break;
        case BuiltinId::Power:
            if (expression.asCall().arguments.size() == 2
                && containsPeriodic(containsPeriodic, expression.asCall().arguments.front())) {
                const Expr& exponent = expression.asCall().arguments[1];
                const bool exactInteger = exponent.isNumber()
                    && exponent.asNumber().isReal()
                    && exponent.asNumber().asReal().isInteger();
                if (!exactInteger)
                    return true;
            }
            break;
        default:
            break;
        }
    }
    for (const Expr& argument : expression.asCall().arguments)
        if (containsRestrictedRealWrapperOverPeriodic(argument, builtins))
            return true;
    return false;
}

[[nodiscard]] std::optional<PlotAnalysis> analyzeAffinePeriodicPoles(
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions) {
    (void)infinitySymbol;
    if (!request.expression.isCall())
        return std::nullopt;

    const Expr* argument = nullptr;
    Rational poleBaseTurns;
    bool directPeriodicPole = false;

    if (request.expression.asCall().arguments.size() == 1) {
        if (const auto* definition = builtins.find(request.expression.asCall().head)) {
            switch (definition->id) {
            case BuiltinId::Tan:
            case BuiltinId::Sec:
            case BuiltinId::Tanc:
                poleBaseTurns = Rational{BigInt{1}, BigInt{4}};
                argument = &request.expression.asCall().arguments.front();
                directPeriodicPole = true;
                break;
            case BuiltinId::Cot:
            case BuiltinId::Csc:
                poleBaseTurns = Rational{};
                argument = &request.expression.asCall().arguments.front();
                directPeriodicPole = true;
                break;
            default:
                break;
            }
        }
    }

    // 1/sin[ax+b], p(x)/cos[ax+b]等はglobal real-domain解析では無限個のzeroを
    // 有限intervalへ表現できない。Plot要求区間内だけexactに列挙し，元式が未定義な点を跨がない。
    if (!argument && builtins.isCallTo(request.expression, BuiltinId::Divide)
        && request.expression.asCall().arguments.size() == 2) {
        const Expr& denominator = request.expression.asCall().arguments[1];
        if (denominator.isCall() && denominator.asCall().arguments.size() == 1) {
            if (const auto* definition = builtins.find(denominator.asCall().head)) {
                if (definition->id == BuiltinId::Sin) {
                    poleBaseTurns = Rational{};
                    argument = &denominator.asCall().arguments.front();
                }
                else if (definition->id == BuiltinId::Cos) {
                    poleBaseTurns = Rational{BigInt{1}, BigInt{4}};
                    argument = &denominator.asCall().arguments.front();
                }
            }
        }
    }
    if (!argument)
        return std::nullopt;

    const Expr& periodicArgument = *argument;
    const auto affine = symbolic::toExpressionPolynomial(
        periodicArgument, request.variable, builtins, mathematics, angles);
    if (!affine || affine->degree() != 1)
        return std::nullopt;
    const auto slope = expression::exact::realRational(affine->coefficient(1));
    if (!slope || slope->isZero())
        return std::nullopt;

    // request端点をargument-spaceへ写し，exact angleとして証明できる場合だけ
    // 周期poleを有限列挙する。証明できない場合はcomplete domainを捏造しない。
    auto endpointAngle = [&](const Expr& endpoint) {
        Expr value = symbolic::substituteSymbol(periodicArgument, request.variable, endpoint);
        value = solver::simplifyForSolve(
            std::move(value), builtins, mathematics, angles, assumptions);
        return mathematics::extractExactAngle(value, builtins, mathematics, angles);
    };
    const auto firstAngle = endpointAngle(request.lower);
    const auto secondAngle = endpointAngle(request.upper);
    if (!firstAngle || !secondAngle)
        return std::nullopt;

    const Rational lowerTurns = std::min(firstAngle->turns, secondAngle->turns);
    const Rational upperTurns = std::max(firstAngle->turns, secondAngle->turns);
    const Rational lowerKValue = (lowerTurns - poleBaseTurns) * Rational{BigInt{2}};
    const Rational upperKValue = (upperTurns - poleBaseTurns) * Rational{BigInt{2}};
    const BigInt firstK = ceilRational(lowerKValue);
    const BigInt lastK = floorRational(upperKValue);

    struct PeriodicPole final {
        Rational turns;
        Expr position;
    };
    std::vector<PeriodicPole> poles;
    for (BigInt k = firstK; k <= lastK; k += BigInt{1}) {
        if (poles.size() >= maxSymbolicLandmarks)
            return std::nullopt;
        const Rational turns = poleBaseTurns + Rational{k, BigInt{2}};

        // affineの定数項やPi式のcanonical形に依存せず，request.lowerを基準に
        // argument-spaceのturn差をx-spaceへexactに逆写像する。
        const Rational deltaTurns = turns - firstAngle->turns;
        Expr delta = angleFromTurns(deltaTurns, builtins, mathematics, angles);
        delta = mathematics::scaleExactExpression(
            Rational{BigInt{1}} / *slope, std::move(delta), builtins);
        Expr position = Expr::call(
            builtins.symbol(BuiltinId::Add),
            std::vector<Expr>{request.lower, std::move(delta)});
        position = solver::simplifyForSolve(
            std::move(position), builtins, mathematics, angles, assumptions);
        poles.push_back(PeriodicPole{turns, std::move(position)});
    }
    if (slope->numerator().isNegative())
        std::reverse(poles.begin(), poles.end());

    // poleの順序はargument-spaceのexact turnから既に確定している。
    // Piを含むx座標を一般KnowledgeContextで再比較するとUnknownへ落ち得るため，
    // endpoint判定もturnのRational equalityだけで行う。
    const bool lowerIsPole = !poles.empty()
        && poles.front().turns == firstAngle->turns;
    const bool upperIsPole = !poles.empty()
        && poles.back().turns == secondAngle->turns;

    PlotAnalysis result;
    result.domain.coverage = PlotDomainCoverage::Complete;
    Expr current = request.lower;
    PlotEndpointInclusion currentInclusion = lowerIsPole
        ? PlotEndpointInclusion::Open
        : PlotEndpointInclusion::Closed;

    for (std::size_t i = 0; i < poles.size(); ++i) {
        const auto& pole = poles[i];
        const bool atLower = i == 0 && lowerIsPole;
        const bool atUpper = i + 1 == poles.size() && upperIsPole;

        if (!atLower) {
            result.domain.intervals.push_back(PlotInterval{
                current, currentInclusion, pole.position, PlotEndpointInclusion::Open});
        }
        current = pole.position;
        currentInclusion = PlotEndpointInclusion::Open;

        addLandmark(result, pole.position, PlotLandmarkKind::DomainBoundary);
        addLandmark(result, pole.position, PlotLandmarkKind::UndefinedPoint);
        if (directPeriodicPole) {
            addLandmark(result, pole.position, PlotLandmarkKind::Pole);
            addLandmark(result, pole.position, PlotLandmarkKind::VerticalAsymptote);
        }
        else {
            // reciprocal trigは分子zeroとの相殺可能性がある。pole点を直接評価せず，
            // 分子がその点で非零と証明できる場合だけpole/asymptoteを付ける。
            const Expr& numerator = request.expression.asCall().arguments[0];
            Expr numeratorAtPole = symbolic::substituteSymbol(
                numerator, request.variable, pole.position);
            numeratorAtPole = solver::simplifyForSolve(
                std::move(numeratorAtPole), builtins, mathematics, angles, assumptions);
            const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
            if (knowledge.prove(mathematics::relation(
                    RelationKind::NotEqual, numeratorAtPole, integerExpr(0)))
                == TruthValue::True) {
                addLandmark(result, pole.position, PlotLandmarkKind::Pole);
                addLandmark(result, pole.position, PlotLandmarkKind::VerticalAsymptote);
            }
        }

        if (atUpper)
            break;
    }

    if (!upperIsPole) {
        result.domain.intervals.push_back(PlotInterval{
            current, currentInclusion, request.upper, PlotEndpointInclusion::Closed});
    }
    return result;
}

} // namespace

PlotAnalysis makeInitialPlotAnalysis(const PlotRequest& request) {
    PlotDomain domain;
    domain.coverage = PlotDomainCoverage::Unknown;
    domain.intervals.push_back(PlotInterval{
        request.lower,
        PlotEndpointInclusion::Closed,
        request.upper,
        PlotEndpointInclusion::Closed});
    return PlotAnalysis{std::move(domain), {}};
}

PlotAnalysis analyzePlotRequest(
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions) {
    mathematics::AssumptionSet realAssumptions = solver::withRealVariable(
        assumptions, request.variable);
    // Plot要求区間そのものもlocal assumptionとして使う。これによりx^x等でも，
    // requestが正の区間に完全に入る場合はprincipal Powerの実性をexactに証明できる。
    // 解析結果は最終的に同じrequest区間へclipするため，区間外を誤って主張しない。
    realAssumptions.add(mathematics::relation(
        RelationKind::GreaterEqual, Expr{request.variable}, request.lower));
    realAssumptions.add(mathematics::relation(
        RelationKind::LessEqual, Expr{request.variable}, request.upper));

    if (auto constrained = analyzeAffinePeriodicRealConstraint(
            request, builtins, mathematics, angles, realAssumptions))
        return *constrained;

    if (auto periodic = analyzeAffinePeriodicPoles(
            request, builtins, mathematics, angles, infinitySymbol, realAssumptions))
        return *periodic;

    PlotAnalysis result = makeInitialPlotAnalysis(request);
    const solver::RealDomainAnalysis realDomain = solver::analyzeRealDomain(
        request.expression, request.variable,
        builtins, mathematics, angles, realAssumptions);
    if (!realDomain.complete) {
        std::vector<PeriodicUndefinedPoint> periodicCuts;
        if (!collectNestedAffinePeriodicUndefinedPoints(
                request.expression, request, builtins, mathematics, angles,
                realAssumptions, periodicCuts)
            || !splitDomainAtOpenCuts(
                result, std::move(periodicCuts), builtins, mathematics, angles, realAssumptions)) {
            result.samplingSafety = PlotSamplingSafety::UnsupportedDiscontinuity;
            return result;
        }
        if (containsRestrictedRealWrapperOverPeriodic(request.expression, builtins)
            || containsPrincipalPowerNeedingBranchProof(
                request.expression, builtins, mathematics, realAssumptions)) {
            result.samplingSafety = PlotSamplingSafety::UnsupportedDiscontinuity;
            return result;
        }

        std::vector<StepDiscontinuity> stepCuts;
        if (!collectStepDiscontinuities(
                request.expression, request, builtins, mathematics, angles,
                realAssumptions, stepCuts)
            || !splitDomainAtStepDiscontinuities(
                result, std::move(stepCuts), builtins, mathematics, angles, realAssumptions))
            result.samplingSafety = PlotSamplingSafety::UnsupportedDiscontinuity;
        return result;
    }

    PlotDomain clipped;
    clipped.coverage = PlotDomainCoverage::Complete;
    for (const solver::RealDomainInterval& interval : realDomain.intervals) {
        std::optional<PlotInterval> clippedInterval;
        if (!clipInterval(
                interval, request, builtins, mathematics, angles, realAssumptions, clippedInterval)) {
            // exactにrequest範囲と比較できない場合は，不完全なintervalを捏造しない。
            return result;
        }
        if (clippedInterval)
            clipped.intervals.push_back(std::move(*clippedInterval));
    }
    result.domain = std::move(clipped);

    std::vector<Expr> boundaries;
    for (const solver::RealDomainInterval& interval : realDomain.intervals) {
        const auto collect = [&](const std::optional<Expr>& endpoint, bool inclusiveEndpoint) {
            if (!endpoint || !insideRequestedRange(
                    *endpoint, request, builtins, mathematics, angles, realAssumptions))
                return;
            addLandmark(result, *endpoint, PlotLandmarkKind::DomainBoundary);
            if (!inclusiveEndpoint)
                addLandmark(result, *endpoint, PlotLandmarkKind::UndefinedPoint);
            if (std::find(boundaries.begin(), boundaries.end(), *endpoint) == boundaries.end())
                boundaries.push_back(*endpoint);
        };
        collect(interval.lower, interval.lowerInclusive);
        collect(interval.upper, interval.upperInclusive);
    }

    if (const auto branchPoint = directBranchPoint(request, builtins);
        branchPoint && insideRequestedRange(
            *branchPoint, request, builtins, mathematics, angles, realAssumptions))
        addLandmark(result, *branchPoint, PlotLandmarkKind::BranchPoint);

    for (const Expr& boundary : boundaries)
        classifyBoundaryCheap(
            result, request, boundary,
            builtins, mathematics, angles, realAssumptions);

    std::vector<PeriodicUndefinedPoint> periodicCuts;
    if (!collectNestedAffinePeriodicUndefinedPoints(
            request.expression, request, builtins, mathematics, angles,
            realAssumptions, periodicCuts)
        || !splitDomainAtOpenCuts(
            result, std::move(periodicCuts), builtins, mathematics, angles, realAssumptions)) {
        result.samplingSafety = PlotSamplingSafety::UnsupportedDiscontinuity;
        return result;
    }
    if (!realDomain.complete
        && containsRestrictedRealWrapperOverPeriodic(request.expression, builtins)) {
        result.samplingSafety = PlotSamplingSafety::UnsupportedDiscontinuity;
        return result;
    }

    std::vector<StepDiscontinuity> stepCuts;
    if (!collectStepDiscontinuities(
            request.expression, request, builtins, mathematics, angles,
            realAssumptions, stepCuts)
        || !splitDomainAtStepDiscontinuities(
            result, std::move(stepCuts), builtins, mathematics, angles, realAssumptions))
        result.samplingSafety = PlotSamplingSafety::UnsupportedDiscontinuity;
    return result;
}

} // namespace mmcal::plot
