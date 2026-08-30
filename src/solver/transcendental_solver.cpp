// 実軸で安全な超越方程式反転
#include "transcendental_solver.hpp"

#include "mathematics/knowledge_context.hpp"
#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "polynomial_solver.hpp"
#include "solve_constraints.hpp"
#include "solver_support.hpp"
#include "symbolic/polynomial.hpp"

#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mmcal::solver {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::FunctionDefinition;
using mathematics::RealRangeRule;
using mathematics::RelationKind;
using mathematics::TruthValue;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return builtins.isCallTo(expression, id);
}

[[nodiscard]] bool isEqualRelation(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 2)
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    return definition && definition->id == BuiltinId::Equal;
}

[[nodiscard]] bool containsVariable(
    const Expr& expression,
    const expression::Symbol& variable) {
    return symbolic::containsSymbol(expression, variable);
}

struct FunctionSide final {
    const FunctionDefinition* definition = nullptr;
    Expr argument;
    Expr rhs;
};

[[nodiscard]] std::optional<FunctionSide> matchFunctionSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const mathematics::MathRegistry& mathematics) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 1
        || !containsVariable(lhs.asCall().arguments[0], variable)
        || containsVariable(rhs, variable))
        return std::nullopt;
    const auto* definition = mathematics.findFunction(lhs.asCall().head);
    if (!definition || !definition->inverseFunction || !definition->realGloballyInjective)
        return std::nullopt;
    return FunctionSide{definition, lhs.asCall().arguments[0], rhs};
}

[[nodiscard]] mathematics::AssumptionSet rangeConditions(
    const FunctionDefinition& definition,
    const Expr& rhs) {
    mathematics::AssumptionSet result;
    result.add(mathematics::elementOf(rhs, mathematics::NumericDomain::Real));
    switch (definition.realRangeRule) {
    case RealRangeRule::Positive:
        result.add(mathematics::relation(RelationKind::Greater, rhs, integerExpr(0)));
        break;
    case RealRangeRule::NonNegative:
        result.add(mathematics::relation(RelationKind::GreaterEqual, rhs, integerExpr(0)));
        break;
    case RealRangeRule::OpenMinusOneToOne:
        result.add(mathematics::relation(RelationKind::Greater, rhs, integerExpr(-1)));
        result.add(mathematics::relation(RelationKind::Less, rhs, integerExpr(1)));
        break;
    case RealRangeRule::OpenZeroToTwo:
        result.add(mathematics::relation(RelationKind::Greater, rhs, integerExpr(0)));
        result.add(mathematics::relation(RelationKind::Less, rhs, integerExpr(2)));
        break;
    case RealRangeRule::ClosedMinusOneToOne:
        result.add(mathematics::relation(RelationKind::GreaterEqual, rhs, integerExpr(-1)));
        result.add(mathematics::relation(RelationKind::LessEqual, rhs, integerExpr(1)));
        break;
    case RealRangeRule::OneToInfinity:
        result.add(mathematics::relation(RelationKind::GreaterEqual, rhs, integerExpr(1)));
        break;
    case RealRangeRule::AllReal:
    case RealRangeRule::Unknown:
        break;
    }
    return result;
}

[[nodiscard]] std::optional<mathematics::AssumptionSet> proveRangeConditions(
    const FunctionDefinition& definition,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    mathematics::AssumptionSet conditions = rangeConditions(definition, rhs);
    mathematics::AssumptionSet remaining;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    for (const auto& predicate : conditions.predicates()) {
        const TruthValue truth = knowledge.prove(predicate);
        if (truth == TruthValue::False)
            return std::nullopt;
        if (truth == TruthValue::Unknown)
            remaining.add(predicate);
    }
    return remaining;
}


[[nodiscard]] std::optional<FunctionSide> matchPeriodicFunctionSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const mathematics::MathRegistry& mathematics) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 1
        || !containsVariable(lhs.asCall().arguments[0], variable)
        || containsVariable(rhs, variable))
        return std::nullopt;

    const auto* definition = mathematics.findFunction(lhs.asCall().head);
    if (!definition || !definition->periodTurns || !definition->inverseFunction
        || definition->realGloballyInjective)
        return std::nullopt;
    if (definition->id != mathematics::FunctionId::Sin
        && definition->id != mathematics::FunctionId::Cos
        && definition->id != mathematics::FunctionId::Tan)
        return std::nullopt;
    return FunctionSide{definition, lhs.asCall().arguments[0], rhs};
}

[[nodiscard]] Expr angleValueFromTurns(
    const numeric::Rational& turns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Degree:
        return Expr{Number{turns * numeric::Rational{BigInt{360}}}};
    case mathematics::AngleUnit::Gradian:
        return Expr{Number{turns * numeric::Rational{BigInt{400}}}};
    case mathematics::AngleUnit::Radian: {
        const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
        if (!pi)
            error::throwCalcError(error::CalcErrorType::Internal, "Pi is not registered");
        const numeric::Rational coefficient = turns * numeric::Rational{BigInt{2}};
        if (coefficient == numeric::Rational{BigInt{1}})
            return Expr{pi->symbol};
        return Expr::call(
            builtins.symbol(BuiltinId::Multiply),
            {Expr{Number{coefficient}}, Expr{pi->symbol}});
    }
    }
    return Expr{Number{BigInt{0}}};
}

void collectSymbolNames(const Expr& expression, std::unordered_set<std::string>& names) {
    if (expression.isSymbol()) {
        names.insert(expression.asSymbol().name());
        return;
    }
    if (expression.isCall()) {
        for (const Expr& argument : expression.asCall().arguments)
            collectSymbolNames(argument, names);
        return;
    }
    if (expression.isArray()) {
        for (std::size_t i = 0; i < expression.asArray().size(); ++i)
            collectSymbolNames(expression.asArray().element(i), names);
        return;
    }
    if (expression.isList())
        for (const Expr& element : expression.asList().elements)
            collectSymbolNames(element, names);
}

[[nodiscard]] expression::Symbol freshIntegerParameter(
    const Expr& relation,
    const expression::Symbol& variable) {
    std::unordered_set<std::string> names;
    collectSymbolNames(relation, names);
    names.insert(variable.name());
    for (std::size_t index = 0;; ++index) {
        const std::string name = index == 0 ? "k" : "k" + std::to_string(index);
        if (!names.contains(name))
            return expression::Symbol{name};
    }
}

[[nodiscard]] Expr periodicTarget(
    Expr base,
    const Expr& period,
    const expression::Symbol& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr multiple = Expr::call(
        builtins.symbol(BuiltinId::Multiply), {period, Expr{parameter}});
    return simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Add), {std::move(base), std::move(multiple)}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] bool exactEndpointOne(const Expr& value) {
    if (!value.isNumber() || !value.asNumber().isReal())
        return false;
    const numeric::Rational rational = value.asNumber().asReal().toRational();
    return rational == numeric::Rational{BigInt{1}}
        || rational == numeric::Rational{BigInt{-1}};
}

[[nodiscard]] bool exactZero(const Expr& value) {
    return value.isNumber() && value.asNumber().isReal() && value.asNumber().isZero();
}


enum class CertifiedOrder {
    Less,
    Equal,
    Greater,
    Unknown
};

[[nodiscard]] CertifiedOrder certifiedConstantOrder(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs == rhs)
        return CertifiedOrder::Equal;

    const mathematics::AssumptionSet noAssumptions;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, noAssumptions};
    if (knowledge.prove(mathematics::relation(RelationKind::Less, lhs, rhs)) == TruthValue::True)
        return CertifiedOrder::Less;
    if (knowledge.prove(mathematics::relation(RelationKind::Greater, lhs, rhs)) == TruthValue::True)
        return CertifiedOrder::Greater;
    if (knowledge.prove(mathematics::relation(RelationKind::Equal, lhs, rhs)) == TruthValue::True)
        return CertifiedOrder::Equal;

    // 定数だけからなる超越式の大小は、guessではなくcertified enclosureが分離した場合だけ採用する。
    const approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    for (const std::size_t bits : {96U, 192U, 384U}) {
        try {
            const auto left = certified.enclose(lhs, bits);
            const auto right = certified.enclose(rhs, bits);
            if (!left || !right || !left->isReal() || !right->isReal())
                return CertifiedOrder::Unknown;
            const auto& l = left->asReal();
            const auto& r = right->asReal();
            if (l.upper() < r.lower())
                return CertifiedOrder::Less;
            if (l.lower() > r.upper())
                return CertifiedOrder::Greater;
            if (l.isPoint() && r.isPoint() && l.lower() == r.lower())
                return CertifiedOrder::Equal;
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            return CertifiedOrder::Unknown;
        }
    }
    return CertifiedOrder::Unknown;
}

[[nodiscard]] bool isVariableSquare(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 2)
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::Power)
        return false;
    const auto& arguments = expression.asCall().arguments;
    if (!arguments[0].isSymbol() || arguments[0].asSymbol() != variable
        || !arguments[1].isNumber() || !arguments[1].asNumber().isReal()
        || !arguments[1].asNumber().asReal().isInteger())
        return false;
    return arguments[1].asNumber().asReal().asInteger() == BigInt{2};
}



struct LambertNormalSide final {
    Expr argument;
    Expr rhs;
};

[[nodiscard]] std::optional<LambertNormalSide> matchLambertNormalSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (containsVariable(rhs, variable)
        || !isHead(lhs, builtins, BuiltinId::Multiply)
        || lhs.asCall().arguments.size() < 2)
        return std::nullopt;

    const auto& factors = lhs.asCall().arguments;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        const Expr& exponential = factors[i];
        if (!isHead(exponential, builtins, BuiltinId::Exp)
            || exponential.asCall().arguments.size() != 1)
            continue;

        std::vector<Expr> remaining;
        remaining.reserve(factors.size() - 1);
        for (std::size_t j = 0; j < factors.size(); ++j)
            if (j != i)
                remaining.push_back(factors[j]);
        Expr multiplier = remaining.size() == 1
            ? remaining.front()
            : Expr::call(builtins.symbol(BuiltinId::Multiply), std::move(remaining));
        multiplier = simplifyForSolve(
            std::move(multiplier), builtins, mathematics, angles, assumptions);
        const Expr& argument = exponential.asCall().arguments[0];
        if (!containsVariable(argument, variable))
            continue;
        Expr difference = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Add), {
                argument, Expr::call(builtins.symbol(BuiltinId::Negate), {multiplier})}),
            builtins, mathematics, angles, assumptions);
        if (exactZero(difference))
            return LambertNormalSide{argument, rhs};
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> matchFunctionPlusArgumentZeroSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    mathematics::FunctionId function) {
    if (!exactZero(rhs)
        || !isHead(lhs, builtins, BuiltinId::Add)
        || lhs.asCall().arguments.size() < 2)
        return std::nullopt;

    const auto& terms = lhs.asCall().arguments;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        const Expr& call = terms[i];
        if (!call.isCall() || call.asCall().arguments.size() != 1)
            continue;
        const auto* definition = mathematics.findFunction(call.asCall().head);
        if (!definition || definition->id != function)
            continue;
        const Expr& argument = call.asCall().arguments[0];
        if (!containsVariable(argument, variable))
            continue;

        std::vector<Expr> remaining;
        remaining.reserve(terms.size() - 1);
        for (std::size_t j = 0; j < terms.size(); ++j)
            if (j != i)
                remaining.push_back(terms[j]);
        Expr rest = remaining.size() == 1
            ? remaining.front()
            : Expr::call(builtins.symbol(BuiltinId::Add), std::move(remaining));
        rest = simplifyForSolve(
            std::move(rest), builtins, mathematics, angles, assumptions);
        Expr difference = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Subtract), {rest, argument}),
            builtins, mathematics, angles, assumptions);
        if (exactZero(difference))
            return argument;
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> matchFunctionNegativeSelfSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    mathematics::FunctionId function) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = mathematics.findFunction(lhs.asCall().head);
    if (!definition || definition->id != function)
        return std::nullopt;
    const Expr& argument = lhs.asCall().arguments[0];
    if (!containsVariable(argument, variable) || !containsVariable(rhs, variable))
        return std::nullopt;
    Expr sum = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Add), {rhs, argument}),
        builtins, mathematics, angles, assumptions);
    if (!exactZero(sum))
        return std::nullopt;
    return argument;
}

[[nodiscard]] std::optional<LambertNormalSide> matchNegativeExponentialFixedPoint(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isHead(lhs, builtins, BuiltinId::Exp)
        || lhs.asCall().arguments.size() != 1
        || !containsVariable(rhs, variable))
        return std::nullopt;

    const Expr& exponent = lhs.asCall().arguments[0];
    if (!containsVariable(exponent, variable))
        return std::nullopt;
    if (isHead(exponent, builtins, BuiltinId::Negate)
        && exponent.asCall().arguments.size() == 1
        && exponent.asCall().arguments[0] == rhs)
        return LambertNormalSide{rhs, integerExpr(1)};
    Expr sum = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Add), {exponent, rhs}),
        builtins, mathematics, angles, assumptions);
    if (!exactZero(sum))
        return std::nullopt;
    return LambertNormalSide{rhs, integerExpr(1)};
}

[[nodiscard]] std::optional<FunctionSide> matchNamedFunctionSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const mathematics::MathRegistry& mathematics,
    mathematics::FunctionId function) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 1
        || !containsVariable(lhs.asCall().arguments[0], variable)
        || containsVariable(rhs, variable))
        return std::nullopt;
    const auto* definition = mathematics.findFunction(lhs.asCall().head);
    if (!definition || definition->id != function)
        return std::nullopt;
    return FunctionSide{definition, lhs.asCall().arguments[0], rhs};
}

[[nodiscard]] std::optional<Expr> matchUnarySelfRelationSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const mathematics::MathRegistry& mathematics,
    mathematics::FunctionId function) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 1)
        return std::nullopt;
    const Expr& argument = lhs.asCall().arguments[0];
    if (!containsVariable(argument, variable) || argument != rhs)
        return std::nullopt;
    const auto* definition = mathematics.findFunction(lhs.asCall().head);
    if (!definition || definition->id != function)
        return std::nullopt;
    return argument;
}

[[nodiscard]] bool provablyRealForRealVariable(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    mathematics::AssumptionSet proofAssumptions = assumptions;
    proofAssumptions.add(mathematics::elementOf(
        Expr{variable}, mathematics::NumericDomain::Real));
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, proofAssumptions};
    return knowledge.prove(mathematics::elementOf(
        expression, mathematics::NumericDomain::Real)) == TruthValue::True;
}

struct ExponentialPowerSide final {
    Expr base;
    Expr exponent;
    Expr rhs;
};

struct RealAffineForm final {
    Expr slope;
    Expr intercept;
};

struct BaseLogSide final {
    Expr base;
    Expr argument;
    Expr rhs;
};

[[nodiscard]] std::optional<BaseLogSide> matchBaseLogSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 2
        || containsVariable(rhs, variable))
        return std::nullopt;
    const auto* definition = builtins.find(lhs.asCall().head);
    if (!definition || definition->id != BuiltinId::Log)
        return std::nullopt;
    const auto& arguments = lhs.asCall().arguments;
    if (containsVariable(arguments[0], variable)
        || !containsVariable(arguments[1], variable))
        return std::nullopt;
    return BaseLogSide{arguments[0], arguments[1], rhs};
}

[[nodiscard]] std::optional<ExponentialPowerSide> matchExponentialPowerSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    // 右辺はx^2等で変数依存してよい。ここでは左辺が「定数base ^ 変数依存指数」かだけを見る。
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto* definition = builtins.find(lhs.asCall().head);
    if (!definition || definition->id != BuiltinId::Power)
        return std::nullopt;
    const auto& arguments = lhs.asCall().arguments;
    if (containsVariable(arguments[0], variable)
        || !containsVariable(arguments[1], variable))
        return std::nullopt;
    return ExponentialPowerSide{arguments[0], arguments[1], rhs};
}

[[nodiscard]] std::optional<RealAffineForm> realAffineForm(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto polynomial = symbolic::toExpressionPolynomial(
        expression, variable, builtins, mathematics, angles);
    if (!polynomial || polynomial->degree() != 1)
        return std::nullopt;

    const Expr& slope = polynomial->coefficient(1);
    const Expr& intercept = polynomial->coefficient(0);
    if (containsVariable(slope, variable) || containsVariable(intercept, variable))
        return std::nullopt;

    mathematics::AssumptionSet proofAssumptions = assumptions;
    proofAssumptions.add(mathematics::elementOf(
        Expr{variable}, mathematics::NumericDomain::Real));
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, proofAssumptions};
    if (knowledge.prove(mathematics::elementOf(
            slope, mathematics::NumericDomain::Real)) != TruthValue::True
        || knowledge.prove(mathematics::elementOf(
            intercept, mathematics::NumericDomain::Real)) != TruthValue::True
        || knowledge.prove(mathematics::relation(
            RelationKind::NotEqual, slope, integerExpr(0))) != TruthValue::True)
        return std::nullopt;

    return RealAffineForm{slope, intercept};
}

[[nodiscard]] std::optional<Expr> expArgument(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::Exp)
        return std::nullopt;
    return expression.asCall().arguments[0];
}

[[nodiscard]] std::optional<Expr> selfPowerConstantTarget(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (containsVariable(rhs, variable)
        || !lhs.isCall() || lhs.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto* definition = builtins.find(lhs.asCall().head);
    if (!definition || definition->id != BuiltinId::Power)
        return std::nullopt;
    const auto& arguments = lhs.asCall().arguments;
    if (!arguments[0].isSymbol() || !arguments[1].isSymbol()
        || !arguments[0].asSymbol().sameIdentity(variable)
        || !arguments[1].asSymbol().sameIdentity(variable))
        return std::nullopt;
    return rhs;
}

[[nodiscard]] Expr lambertW(
    int branch,
    Expr argument,
    const evaluation::BuiltinRegistry& builtins) {
    if (branch == 0)
        return Expr::call(builtins.symbol(BuiltinId::LambertW), {std::move(argument)});
    return Expr::call(
        builtins.symbol(BuiltinId::LambertW),
        {integerExpr(branch), std::move(argument)});
}

[[nodiscard]] SolutionBranch realBinding(
    const expression::Symbol& variable,
    Expr value,
    mathematics::AssumptionSet conditions = {}) {
    SolutionBranch branch;
    branch.bindings.push_back(SolutionBinding{variable, std::move(value)});
    branch.conditions = std::move(conditions);
    branch.bindingsCertifiedDomain = mathematics::NumericDomain::Real;
    return branch;
}

[[nodiscard]] std::optional<SolutionSet> solveRealAffineExponentialEquation(
    const RealAffineForm& exponent,
    const RealAffineForm& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    mathematics::AssumptionSet proofAssumptions = assumptions;
    proofAssumptions.add(mathematics::elementOf(
        Expr{variable}, mathematics::NumericDomain::Real));

    // exp[p x+q] == c x+d を y=c x+d へ移し，
    // -(p/c)y Exp[-(p/c)y] == -(p/c)Exp[q-pd/c] としてLambert Wへ落とす。
    Expr pdOverC = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Divide), {
            Expr::call(builtins.symbol(BuiltinId::Multiply),
                {exponent.slope, rhs.intercept}),
            rhs.slope}),
        builtins, mathematics, angles, proofAssumptions);
    Expr shiftedExponent = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Subtract), {
            exponent.intercept, std::move(pdOverC)}),
        builtins, mathematics, angles, proofAssumptions);
    Expr pOverC = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Divide), {
            exponent.slope, rhs.slope}),
        builtins, mathematics, angles, proofAssumptions);
    Expr z = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Negate), {
            Expr::call(builtins.symbol(BuiltinId::Multiply), {
                pOverC,
                Expr::call(builtins.symbol(BuiltinId::Exp), {std::move(shiftedExponent)})})}),
        builtins, mathematics, angles, proofAssumptions);

    const CertifiedOrder zToZero = certifiedConstantOrder(
        z, integerExpr(0), builtins, mathematics, angles);
    if (zToZero == CertifiedOrder::Equal || zToZero == CertifiedOrder::Unknown)
        return std::nullopt;

    const auto* e = mathematics.findConstant(mathematics::ConstantId::E);
    if (!e)
        error::throwCalcError(error::CalcErrorType::Internal, "E is not registered");
    Expr branchPoint = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Negate), {
            Expr::call(builtins.symbol(BuiltinId::Divide), {integerExpr(1), Expr{e->symbol}})}),
        builtins, mathematics, angles, proofAssumptions);

    auto rootFromBranch = [&](int branch) {
        Expr first = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Divide), {
                Expr::call(builtins.symbol(BuiltinId::Negate), {
                    lambertW(branch, z, builtins)}),
                exponent.slope}),
            builtins, mathematics, angles, proofAssumptions);
        Expr shift = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Divide), {
                rhs.intercept, rhs.slope}),
            builtins, mathematics, angles, proofAssumptions);
        return simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Subtract), {
                std::move(first), std::move(shift)}),
            builtins, mathematics, angles, proofAssumptions);
    };

    if (zToZero == CertifiedOrder::Greater)
        return SolutionSet::finite(
            {SolverVariable{variable, mathematics::NumericDomain::Real}},
            {realBinding(variable, rootFromBranch(0))});

    const CertifiedOrder zToBranchPoint = certifiedConstantOrder(
        z, branchPoint, builtins, mathematics, angles);
    if (zToBranchPoint == CertifiedOrder::Less)
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});
    if (zToBranchPoint == CertifiedOrder::Unknown)
        return std::nullopt;
    if (zToBranchPoint == CertifiedOrder::Equal) {
        // W_0(-1/e)=W_-1(-1/e)=-1なので，branch pointでは重複を作らず直接値へ戻す。
        Expr first = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Divide), {
                integerExpr(1), exponent.slope}),
            builtins, mathematics, angles, proofAssumptions);
        Expr shift = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Divide), {
                rhs.intercept, rhs.slope}),
            builtins, mathematics, angles, proofAssumptions);
        Expr value = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Subtract), {
                std::move(first), std::move(shift)}),
            builtins, mathematics, angles, proofAssumptions);
        return SolutionSet::finite(
            {SolverVariable{variable, mathematics::NumericDomain::Real}},
            {realBinding(variable, std::move(value))});
    }

    return SolutionSet::finite(
        {SolverVariable{variable, mathematics::NumericDomain::Real}},
        {realBinding(variable, rootFromBranch(0)),
         realBinding(variable, rootFromBranch(-1))});
}

[[nodiscard]] std::optional<SolutionSet> solveRealSelfPowerConstantEquation(
    const Expr& target,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    mathematics::AssumptionSet proofAssumptions = assumptions;
    proofAssumptions.add(mathematics::elementOf(
        Expr{variable}, mathematics::NumericDomain::Real));
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, proofAssumptions};
    if (knowledge.prove(mathematics::elementOf(
            target, mathematics::NumericDomain::Real)) != TruthValue::True)
        return std::nullopt;

    const CertifiedOrder toZero = certifiedConstantOrder(
        target, integerExpr(0), builtins, mathematics, angles);
    const CertifiedOrder toOne = certifiedConstantOrder(
        target, integerExpr(1), builtins, mathematics, angles);
    const CertifiedOrder toMinusOne = certifiedConstantOrder(
        target, integerExpr(-1), builtins, mathematics, angles);

    if (toZero == CertifiedOrder::Equal)
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});
    if (toOne == CertifiedOrder::Equal)
        return SolutionSet::finite(
            {SolverVariable{variable, mathematics::NumericDomain::Real}},
            {realBinding(variable, integerExpr(1))});
    if (toMinusOne == CertifiedOrder::Equal)
        return SolutionSet::finite(
            {SolverVariable{variable, mathematics::NumericDomain::Real}},
            {realBinding(variable, integerExpr(-1))});
    if (toMinusOne == CertifiedOrder::Less)
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    // principal x^xでx<0が正実数になるにはxが負の偶整数でなければならず，
    // その値は|x|^x<=1/4。したがってtarget>1では負実根は存在せず，
    // x>0上の log変換 x log[x]==log[target] が完全である。
    if (toOne != CertifiedOrder::Greater)
        return std::nullopt;

    Expr logTarget = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Log), {target}),
        builtins, mathematics, angles, proofAssumptions);
    Expr value = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Exp), {
            lambertW(0, std::move(logTarget), builtins)}),
        builtins, mathematics, angles, proofAssumptions);
    return SolutionSet::finite(
        {SolverVariable{variable, mathematics::NumericDomain::Real}},
        {realBinding(variable, std::move(value))});
}

[[nodiscard]] std::optional<SolutionSet> solveRealTargetEquation(
    const Expr& argument,
    Expr target,
    const expression::Symbol& variable,
    const mathematics::AssumptionSet& branchConditions,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (argument.isSymbol() && argument.asSymbol().sameIdentity(variable)) {
        return SolutionSet::finite(
            {SolverVariable{variable, mathematics::NumericDomain::Real}},
            {realBinding(variable, std::move(target), branchConditions)});
    }

    // Polynomial solverはLambertW等の一般x-free symbolic targetを係数として
    // 受理しない経路がある。argument自体がexact affineならここで直接反転し，
    // inverse-function/Lambertの結果を不要にUnresolvedへ落とさない。
    if (const auto affine = symbolic::toExpressionPolynomial(
            argument, variable, builtins, mathematics, angles);
        affine && affine->degree() == 1) {
        const Expr& slope = affine->coefficient(1);
        if (slope.isNumber() && slope.asNumber().isReal() && !slope.asNumber().isZero()) {
            Expr numerator = simplifyForSolve(
                Expr::call(builtins.symbol(BuiltinId::Subtract),
                    {std::move(target), affine->coefficient(0)}),
                builtins, mathematics, angles, {});
            Expr value = simplifyForSolve(
                Expr::call(builtins.symbol(BuiltinId::Divide),
                    {std::move(numerator), slope}),
                builtins, mathematics, angles, {});
            return SolutionSet::finite(
                {SolverVariable{variable, mathematics::NumericDomain::Real}},
                {realBinding(variable, std::move(value), branchConditions)});
        }
    }

    Expr relation = Expr::call(
        builtins.symbol(BuiltinId::Equal), {argument, std::move(target)});
    SolutionSet solved = solveUnivariatePolynomialRelation(
        relation, variable, builtins, mathematics, angles);
    if (solved.kind() != SolutionSetKind::Finite)
        return std::nullopt;

    std::vector<SolutionBranch> branches(solved.branches().begin(), solved.branches().end());
    for (SolutionBranch& branch : branches) {
        branch.bindingsCertifiedDomain = mathematics::NumericDomain::Real;
        for (const auto& predicate : branchConditions.predicates())
            branch.conditions.add(predicate);
    }
    return SolutionSet::finite(
        std::vector<SolverVariable>{solved.variables().begin(), solved.variables().end()},
        std::move(branches)).withAdditionalConditions(solved.conditions());
}

[[nodiscard]] std::optional<SolutionSet> mergeFiniteTargets(
    std::vector<SolutionSet> sets,
    const expression::Symbol& variable) {
    std::vector<SolutionBranch> branches;
    mathematics::AssumptionSet globalConditions;
    for (const SolutionSet& set : sets) {
        if (set.kind() != SolutionSetKind::Finite)
            return std::nullopt;
        branches.insert(branches.end(), set.branches().begin(), set.branches().end());
        for (const auto& predicate : set.conditions().predicates())
            globalConditions.add(predicate);
    }
    SolutionSet result = SolutionSet::finite(
        {SolverVariable{variable, mathematics::NumericDomain::Real}}, std::move(branches));
    return result.withAdditionalConditions(globalConditions);
}

[[nodiscard]] std::optional<SolutionSet> solveLambertNormalForm(
    const LambertNormalSide& matched,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    mathematics::AssumptionSet proofAssumptions = assumptions;
    proofAssumptions.add(mathematics::elementOf(Expr{variable}, mathematics::NumericDomain::Real));
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, proofAssumptions};

    const TruthValue rhsReal = knowledge.prove(mathematics::elementOf(
        matched.rhs, mathematics::NumericDomain::Real));
    if (rhsReal == TruthValue::False)
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    const auto* e = mathematics.findConstant(mathematics::ConstantId::E);
    if (!e)
        error::throwCalcError(error::CalcErrorType::Internal, "E is not registered");
    Expr branchPoint = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Negate), {
            Expr::call(builtins.symbol(BuiltinId::Divide), {integerExpr(1), Expr{e->symbol}})}),
        builtins, mathematics, angles, proofAssumptions);

    const CertifiedOrder toBranchPoint = certifiedConstantOrder(
        matched.rhs, branchPoint, builtins, mathematics, angles);
    const CertifiedOrder toZero = certifiedConstantOrder(
        matched.rhs, integerExpr(0), builtins, mathematics, angles);
    if (toBranchPoint == CertifiedOrder::Less)
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    std::vector<SolutionSet> targets;
    mathematics::AssumptionSet principalConditions;
    if (rhsReal == TruthValue::Unknown)
        principalConditions.add(mathematics::elementOf(
            matched.rhs, mathematics::NumericDomain::Real));
    if (toBranchPoint == CertifiedOrder::Unknown)
        principalConditions.add(mathematics::relation(
            RelationKind::GreaterEqual, matched.rhs, branchPoint));

    auto principal = solveRealTargetEquation(
        matched.argument, lambertW(0, matched.rhs, builtins), variable,
        principalConditions, builtins, mathematics, angles);
    if (!principal)
        return std::nullopt;
    targets.push_back(std::move(*principal));

    // W_-1が第2の実逆分岐を与えるのは(-1/e,0)だけである。
    // At -1/e it coincides with W_0, so strict lower-bound avoids duplicate roots.
    if (toZero == CertifiedOrder::Less
        && toBranchPoint != CertifiedOrder::Equal) {
        mathematics::AssumptionSet lowerConditions;
        if (rhsReal == TruthValue::Unknown)
            lowerConditions.add(mathematics::elementOf(
                matched.rhs, mathematics::NumericDomain::Real));
        if (toBranchPoint == CertifiedOrder::Unknown)
            lowerConditions.add(mathematics::relation(
                RelationKind::Greater, matched.rhs, branchPoint));
        auto lower = solveRealTargetEquation(
            matched.argument, lambertW(-1, matched.rhs, builtins), variable,
            lowerConditions, builtins, mathematics, angles);
        if (!lower)
            return std::nullopt;
        targets.push_back(std::move(*lower));
    }
    else if (toZero == CertifiedOrder::Unknown) {
        mathematics::AssumptionSet lowerConditions;
        lowerConditions.add(mathematics::elementOf(
            matched.rhs, mathematics::NumericDomain::Real));
        lowerConditions.add(mathematics::relation(
            RelationKind::Greater, matched.rhs, branchPoint));
        lowerConditions.add(mathematics::relation(
            RelationKind::Less, matched.rhs, integerExpr(0)));
        auto lower = solveRealTargetEquation(
            matched.argument, lambertW(-1, matched.rhs, builtins), variable,
            lowerConditions, builtins, mathematics, angles);
        if (!lower)
            return std::nullopt;
        targets.push_back(std::move(*lower));
    }

    return mergeFiniteTargets(std::move(targets), variable);
}

[[nodiscard]] std::optional<SolutionSet> solvePeriodicTarget(
    const Expr& argument,
    const Expr& target,
    const expression::Symbol& variable,
    const expression::Symbol& parameter,
    const mathematics::AssumptionSet& remainingConditions,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr relation = Expr::call(
        builtins.symbol(BuiltinId::Equal), {argument, target});
    SolutionSet result = solveUnivariatePolynomialRelation(
        relation, variable, builtins, mathematics, angles);
    if (result.kind() != SolutionSetKind::Finite)
        return std::nullopt;

    std::vector<SolutionBranch> branches(result.branches().begin(), result.branches().end());
    for (SolutionBranch& branch : branches) {
        branch.freeVariables.push_back(
            SolverVariable{parameter, mathematics::NumericDomain::Integer});
        for (const auto& predicate : remainingConditions.predicates())
            branch.conditions.add(predicate);
    }
    return SolutionSet::finite(
        std::vector<SolverVariable>{result.variables().begin(), result.variables().end()},
        std::move(branches)).withAdditionalConditions(result.conditions());
}

[[nodiscard]] mathematics::Predicate simplifyPredicate(
    const mathematics::Predicate& predicate,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate))
        return mathematics::relation(
            relation->relation,
            simplifyForSolve(relation->lhs, builtins, mathematics, angles, assumptions),
            simplifyForSolve(relation->rhs, builtins, mathematics, angles, assumptions));
    const auto& domain = std::get<mathematics::DomainPredicate>(predicate);
    return mathematics::elementOf(
        simplifyForSolve(domain.expression, builtins, mathematics, angles, assumptions), domain.domain);
}

} // namespace

std::optional<SolutionSet> solvePrincipalLambertRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isEqualRelation(relation, builtins))
        return std::nullopt;

    const auto& sides = relation.asCall().arguments;
    auto matched = matchNamedFunctionSide(
        sides[0], sides[1], variable, mathematics, mathematics::FunctionId::LambertW);
    if (!matched)
        matched = matchNamedFunctionSide(
            sides[1], sides[0], variable, mathematics, mathematics::FunctionId::LambertW);
    if (!matched
        || !matched->argument.isSymbol()
        || !matched->argument.asSymbol().sameIdentity(variable))
        return std::nullopt;

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::elementOf(
            matched->rhs, mathematics::NumericDomain::Real)) != TruthValue::True
        || knowledge.prove(mathematics::relation(
            RelationKind::GreaterEqual, matched->rhs, integerExpr(-1))) != TruthValue::True)
        return std::nullopt;

    Expr target = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Multiply), {
            matched->rhs,
            Expr::call(builtins.symbol(BuiltinId::Exp), {matched->rhs})}),
        builtins, mathematics, angles, assumptions);
    SolutionBranch branch;
    branch.bindings.push_back(SolutionBinding{variable, std::move(target)});
    branch.bindingsCertifiedDomain = mathematics::NumericDomain::Real;
    return SolutionSet::finite(
        {SolverVariable{variable, mathematics::NumericDomain::Complex}},
        {std::move(branch)});
}

std::optional<SolutionSet> solveRealExponentialRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isEqualRelation(relation, builtins))
        return std::nullopt;

    const auto& sides = relation.asCall().arguments;

    // Canonical Lambert normal form u Exp[u] == a.
    // This is a family rule, not an equation-specific rewrite: the argument u may itself be
    // any polynomial expression that the exact polynomial solver can invert afterwards.
    auto lambertNormal = matchLambertNormalSide(
        sides[0], sides[1], variable, builtins, mathematics, angles, assumptions);
    if (!lambertNormal)
        lambertNormal = matchLambertNormalSide(
            sides[1], sides[0], variable, builtins, mathematics, angles, assumptions);
    if (lambertNormal)
        return solveLambertNormalForm(
            *lambertNormal, variable, builtins, mathematics, angles, assumptions);

    // log[u]+u==0 は log[u]==-u，したがって u Exp[u]==1 と同値。
    // 変換後の実根W(1)>0がLogの実定義域も自動的に満たすため，余分な候補を作らない。
    auto logNegativeSelf = matchFunctionNegativeSelfSide(
        sides[0], sides[1], variable, builtins, mathematics, angles, assumptions,
        mathematics::FunctionId::Log);
    if (!logNegativeSelf)
        logNegativeSelf = matchFunctionNegativeSelfSide(
            sides[1], sides[0], variable, builtins, mathematics, angles, assumptions,
            mathematics::FunctionId::Log);
    if (!logNegativeSelf)
        logNegativeSelf = matchFunctionPlusArgumentZeroSide(
            sides[0], sides[1], variable, builtins, mathematics, angles, assumptions,
            mathematics::FunctionId::Log);
    if (!logNegativeSelf)
        logNegativeSelf = matchFunctionPlusArgumentZeroSide(
            sides[1], sides[0], variable, builtins, mathematics, angles, assumptions,
            mathematics::FunctionId::Log);
    if (logNegativeSelf && provablyRealForRealVariable(
            *logNegativeSelf, variable, builtins, mathematics, assumptions))
        return solveLambertNormalForm(
            LambertNormalSide{*logNegativeSelf, integerExpr(1)},
            variable, builtins, mathematics, angles, assumptions);

    // exp[u]+u==0 はv=-uと置けば exp[-v]==v，すなわちv Exp[v]==1。
    auto expNegativeSelf = matchFunctionNegativeSelfSide(
        sides[0], sides[1], variable, builtins, mathematics, angles, assumptions,
        mathematics::FunctionId::Exp);
    if (!expNegativeSelf)
        expNegativeSelf = matchFunctionNegativeSelfSide(
            sides[1], sides[0], variable, builtins, mathematics, angles, assumptions,
            mathematics::FunctionId::Exp);
    if (!expNegativeSelf)
        expNegativeSelf = matchFunctionPlusArgumentZeroSide(
            sides[0], sides[1], variable, builtins, mathematics, angles, assumptions,
            mathematics::FunctionId::Exp);
    if (!expNegativeSelf)
        expNegativeSelf = matchFunctionPlusArgumentZeroSide(
            sides[1], sides[0], variable, builtins, mathematics, angles, assumptions,
            mathematics::FunctionId::Exp);
    if (expNegativeSelf && provablyRealForRealVariable(
            *expNegativeSelf, variable, builtins, mathematics, assumptions)) {
        Expr negated = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Negate), {*expNegativeSelf}),
            builtins, mathematics, angles, assumptions);
        return solveLambertNormalForm(
            LambertNormalSide{std::move(negated), integerExpr(1)},
            variable, builtins, mathematics, angles, assumptions);
    }

    // Exp[-u] == u is algebraically equivalent on the real axis to u Exp[u] == 1.
    // Matching is structural after safe Solve normalization, so exp[-x]==x and the same
    // pattern with a polynomial u share one rule.
    auto fixedPoint = matchNegativeExponentialFixedPoint(
        sides[0], sides[1], variable, builtins, mathematics, angles, assumptions);
    if (!fixedPoint)
        fixedPoint = matchNegativeExponentialFixedPoint(
            sides[1], sides[0], variable, builtins, mathematics, angles, assumptions);
    if (fixedPoint)
        return solveLambertNormalForm(
            *fixedPoint, variable, builtins, mathematics, angles, assumptions);

    // Exp[p x+q] == c x+d はLambert Wで完全に反転できる。
    // branch point判定までexact/certifiedに閉じる場合だけ有限実解を返す。
    auto expSide = expArgument(sides[0], builtins);
    const Expr* expRhs = &sides[1];
    if (!expSide) {
        expSide = expArgument(sides[1], builtins);
        expRhs = &sides[0];
    }
    if (expSide) {
        const auto exponentAffine = realAffineForm(
            *expSide, variable, builtins, mathematics, angles, assumptions);
        const auto rhsAffine = realAffineForm(
            *expRhs, variable, builtins, mathematics, angles, assumptions);
        if (exponentAffine && rhsAffine) {
            if (auto solved = solveRealAffineExponentialEquation(
                    *exponentAffine, *rhsAffine, variable,
                    builtins, mathematics, angles, assumptions))
                return solved;
        }
    }

    // principal x^xは負実軸で一般に複素値になるため，正のbranchだけを
    // 機械的に採用しない。target>=1等，負実根を完全に排除できる範囲だけ扱う。
    auto selfPower = selfPowerConstantTarget(
        sides[0], sides[1], variable, builtins);
    if (!selfPower)
        selfPower = selfPowerConstantTarget(
            sides[1], sides[0], variable, builtins);
    if (selfPower) {
        if (auto solved = solveRealSelfPowerConstantEquation(
                *selfPower, variable, builtins, mathematics, angles, assumptions))
            return solved;
    }

    auto matched = matchExponentialPowerSide(sides[0], sides[1], variable, builtins);
    if (!matched)
        matched = matchExponentialPowerSide(sides[1], sides[0], variable, builtins);
    if (!matched)
        return std::nullopt;

    mathematics::AssumptionSet proofAssumptions = assumptions;
    proofAssumptions.add(mathematics::elementOf(Expr{variable}, mathematics::NumericDomain::Real));
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, proofAssumptions};

    // a>0 かつ実指数なら principal Power[a,u] は exp[u log[a]] と一致し常に正。
    // 0^uや負baseを同じ規則で消すとdomain/branchを壊すため、証明できる場合だけ使う。
    if (knowledge.prove(mathematics::relation(
            RelationKind::Greater, matched->base, integerExpr(0))) != TruthValue::True
        || knowledge.prove(mathematics::elementOf(
            matched->exponent, mathematics::NumericDomain::Real)) != TruthValue::True)
        return std::nullopt;

    if (exactZero(matched->rhs))
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    // positive constant baseのaffine指数はExpへexactに書き換えられるので，
    // a^(m x+n)==c x+d も同じLambert W classifierへ送る。
    if (const auto exponentAffine = realAffineForm(
            matched->exponent, variable, builtins, mathematics, angles, assumptions);
        exponentAffine) {
        if (const auto rhsAffine = realAffineForm(
                matched->rhs, variable, builtins, mathematics, angles, assumptions);
            rhsAffine) {
            const CertifiedOrder baseToOne = certifiedConstantOrder(
                matched->base, integerExpr(1), builtins, mathematics, angles);
            if (baseToOne == CertifiedOrder::Equal) {
                Expr transformed = Expr::call(
                    builtins.symbol(BuiltinId::Equal), {integerExpr(1), matched->rhs});
                return solveUnivariatePolynomialRelation(
                    transformed, variable, builtins, mathematics, angles);
            }
            if (baseToOne != CertifiedOrder::Unknown) {
                Expr logBase = simplifyForSolve(
                    Expr::call(builtins.symbol(BuiltinId::Log), {matched->base}),
                    builtins, mathematics, angles, proofAssumptions);
                RealAffineForm scaledExponent{
                    simplifyForSolve(
                        Expr::call(builtins.symbol(BuiltinId::Multiply), {
                            exponentAffine->slope, logBase}),
                        builtins, mathematics, angles, proofAssumptions),
                    simplifyForSolve(
                        Expr::call(builtins.symbol(BuiltinId::Multiply), {
                            exponentAffine->intercept, logBase}),
                        builtins, mathematics, angles, proofAssumptions)};
                if (auto solved = solveRealAffineExponentialEquation(
                        scaledExponent, *rhsAffine, variable,
                        builtins, mathematics, angles, assumptions))
                    return solved;
            }
        }
    }

    // a^u == d でdが変数に依存しない場合は，a>0かつa!=1の実軸上で
    // u == log[a,d] と完全に同値。Powerを一般Evaluatorへ通さず，HoldAllのSolve内で
    // 証明済みの単射性だけを使って指数側をpolynomial solverへ渡す。
    if (!containsVariable(matched->rhs, variable)) {
        const TruthValue rhsReal = knowledge.prove(mathematics::elementOf(
            matched->rhs, mathematics::NumericDomain::Real));
        const TruthValue rhsPositive = knowledge.prove(mathematics::relation(
            RelationKind::Greater, matched->rhs, integerExpr(0)));
        if (rhsReal == TruthValue::True
            && knowledge.prove(mathematics::relation(
                RelationKind::LessEqual, matched->rhs, integerExpr(0))) == TruthValue::True)
            return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

        const CertifiedOrder baseToOne = certifiedConstantOrder(
            matched->base, integerExpr(1), builtins, mathematics, angles);
        if (baseToOne == CertifiedOrder::Equal) {
            const TruthValue rhsOne = knowledge.prove(mathematics::relation(
                RelationKind::Equal, matched->rhs, integerExpr(1)));
            if (rhsOne == TruthValue::True)
                return SolutionSet::universal(
                    {SolverVariable{variable, mathematics::NumericDomain::Real}});
            if (rhsOne == TruthValue::False)
                return SolutionSet::empty(
                    {SolverVariable{variable, mathematics::NumericDomain::Real}});
            return std::nullopt;
        }

        if (rhsReal == TruthValue::True && rhsPositive == TruthValue::True
            && baseToOne != CertifiedOrder::Unknown) {
            Expr target = simplifyForSolve(
                Expr::call(builtins.symbol(BuiltinId::Log),
                    {matched->base, matched->rhs}),
                builtins, mathematics, angles, proofAssumptions);
            Expr transformed = Expr::call(
                builtins.symbol(BuiltinId::Equal), {matched->exponent, std::move(target)});
            return solveUnivariatePolynomialRelation(
                transformed, variable, builtins, mathematics, angles);
        }
    }

    // variable-dependent RHSをLambert Wへ落とす初版は a^x == x^2 に限定する。
    // affine指数＋constant RHSは上で処理済み。一般のnonconstant P(x)はbranch/domain条件が増えるため別段階へ送る。
    if (!matched->exponent.isSymbol() || matched->exponent.asSymbol() != variable
        || !isVariableSquare(matched->rhs, variable, builtins))
        return std::nullopt;

    const CertifiedOrder baseToOne = certifiedConstantOrder(
        matched->base, integerExpr(1), builtins, mathematics, angles);
    if (baseToOne == CertifiedOrder::Unknown)
        return std::nullopt;
    if (baseToOne == CertifiedOrder::Equal) {
        Expr transformed = Expr::call(
            builtins.symbol(BuiltinId::Equal), {integerExpr(1), matched->rhs});
        return solveUnivariatePolynomialRelation(
            transformed, variable, builtins, mathematics, angles);
    }

    Expr logBase = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Log), {matched->base}),
        builtins, mathematics, angles, proofAssumptions);
    Expr magnitudeLog = baseToOne == CertifiedOrder::Greater
        ? logBase
        : simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Negate), {logBase}),
            builtins, mathematics, angles, proofAssumptions);
    Expr halfMagnitude = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Divide), {magnitudeLog, integerExpr(2)}),
        builtins, mathematics, angles, proofAssumptions);
    // base<1ではmagnitudeLog=-logBaseなので、-(magnitude/2)を機械的に作ると
    // --log[a]/2という非canonical形が残り得る。符号証明済みのlogBaseから直接構成する。
    Expr negativeHalfMagnitude = simplifyForSolve(
        Expr::call(
            builtins.symbol(BuiltinId::Divide),
            {baseToOne == CertifiedOrder::Greater
                ? Expr::call(builtins.symbol(BuiltinId::Negate), {logBase})
                : logBase,
             integerExpr(2)}),
        builtins, mathematics, angles, proofAssumptions);
    Expr scale = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Divide), {integerExpr(-2), logBase}),
        builtins, mathematics, angles, proofAssumptions);

    auto scaledW = [&](int branch, Expr argument) {
        return simplifyForSolve(
            Expr::call(
                builtins.symbol(BuiltinId::Multiply),
                {scale, lambertW(branch, std::move(argument), builtins)}),
            builtins, mathematics, angles, proofAssumptions);
    };

    std::vector<SolutionBranch> branches;
    // 符号がbase-1と反対側の根は、W_0(|log a|/2) が正実数上で常に存在するため無条件。
    branches.push_back(realBinding(variable, scaledW(0, halfMagnitude)));

    // 同符号側の2根は -|log a|/2 >= -1/e、すなわち |log a| <= 2/e のときだけ実在する。
    const auto* e = mathematics.findConstant(mathematics::ConstantId::E);
    if (!e)
        error::throwCalcError(error::CalcErrorType::Internal, "E is not registered");
    Expr threshold = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Divide), {integerExpr(2), Expr{e->symbol}}),
        builtins, mathematics, angles, proofAssumptions);
    const CertifiedOrder branchCondition = certifiedConstantOrder(
        magnitudeLog, threshold, builtins, mathematics, angles);

    if (branchCondition == CertifiedOrder::Less) {
        branches.push_back(realBinding(variable, scaledW(0, negativeHalfMagnitude)));
        branches.push_back(realBinding(variable, scaledW(-1, negativeHalfMagnitude)));
    }
    else if (branchCondition == CertifiedOrder::Equal) {
        // branch point -1/E では W_0 = W_-1 = -1。同じ実根を二重に返さない。
        branches.push_back(realBinding(variable, scaledW(0, negativeHalfMagnitude)));
    }
    else if (branchCondition == CertifiedOrder::Unknown) {
        mathematics::AssumptionSet principalCondition;
        principalCondition.add(mathematics::relation(
            RelationKind::LessEqual, magnitudeLog, threshold));
        branches.push_back(realBinding(
            variable, scaledW(0, negativeHalfMagnitude), std::move(principalCondition)));

        // lower branch はbranch pointでprincipal branchと一致するので、重複回避のためstrict条件にする。
        mathematics::AssumptionSet lowerCondition;
        lowerCondition.add(mathematics::relation(
            RelationKind::Less, magnitudeLog, threshold));
        branches.push_back(realBinding(
            variable, scaledW(-1, negativeHalfMagnitude), std::move(lowerCondition)));
    }

    return SolutionSet::finite(
        {SolverVariable{variable, mathematics::NumericDomain::Real}}, std::move(branches));
}

std::optional<SolutionSet> solveRealPeriodicFunctionRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isEqualRelation(relation, builtins))
        return std::nullopt;

    const auto& sides = relation.asCall().arguments;
    auto matched = matchPeriodicFunctionSide(sides[0], sides[1], variable, mathematics);
    if (!matched)
        matched = matchPeriodicFunctionSide(sides[1], sides[0], variable, mathematics);
    if (!matched)
        return std::nullopt;

    // 初版はargumentが変数についてaffineで、一次係数がexact非零の場合に限定する。
    // symbolic coefficientや非線形argumentのparameter条件を不完全に返さない。
    const auto polynomial = symbolic::toExpressionPolynomial(
        matched->argument, variable, builtins, mathematics, angles);
    if (!polynomial || polynomial->degree() != 1)
        return std::nullopt;
    const Expr& linearCoefficient = polynomial->coefficient(1);
    if (!linearCoefficient.isNumber() || !linearCoefficient.asNumber().isReal()
        || linearCoefficient.asNumber().isZero())
        return std::nullopt;

    matched->rhs = simplifyForSolve(
        matched->rhs, builtins, mathematics, angles, assumptions);
    const auto remainingConditions = proveRangeConditions(
        *matched->definition, matched->rhs, builtins, mathematics, assumptions);
    if (!remainingConditions)
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    const auto* inverse = mathematics.findFunction(*matched->definition->inverseFunction);
    if (!inverse)
        return std::nullopt;

    const expression::Symbol parameter = freshIntegerParameter(relation, variable);
    Expr period = simplifyForSolve(
        angleValueFromTurns(*matched->definition->periodTurns, builtins, mathematics, angles),
        builtins, mathematics, angles, assumptions);
    const Expr inverseValue = simplifyForSolve(
        Expr::call(inverse->symbol, {matched->rhs}),
        builtins, mathematics, angles, assumptions);

    std::vector<Expr> bases;
    switch (matched->definition->id) {
    case mathematics::FunctionId::Tan:
        bases.push_back(inverseValue);
        break;
    case mathematics::FunctionId::Cos:
        if (exactZero(matched->rhs)) {
            // cos[u]==0 は ±quarter-turn + full-turn*k の2branchを
            // quarter-turn + half-turn*k へ一意にまとめられる。
            bases.push_back(angleValueFromTurns(
                numeric::Rational{BigInt{1}, BigInt{4}}, builtins, mathematics, angles));
            period = angleValueFromTurns(
                numeric::Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles);
        }
        else {
            bases.push_back(inverseValue);
            if (!exactEndpointOne(matched->rhs))
                bases.push_back(simplifyForSolve(
                    Expr::call(builtins.symbol(BuiltinId::Negate), {inverseValue}),
                    builtins, mathematics, angles, assumptions));
        }
        break;
    case mathematics::FunctionId::Sin:
        if (exactZero(matched->rhs)) {
            // sin[u]==0 は full-turnごとの2branchではなくhalf-turn*kがcanonical。
            bases.push_back(integerExpr(0));
            period = angleValueFromTurns(
                numeric::Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles);
        }
        else {
            bases.push_back(inverseValue);
            if (!exactEndpointOne(matched->rhs)) {
                Expr halfTurn = angleValueFromTurns(
                    numeric::Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles);
                bases.push_back(simplifyForSolve(
                    Expr::call(builtins.symbol(BuiltinId::Subtract),
                        {std::move(halfTurn), inverseValue}),
                    builtins, mathematics, angles, assumptions));
            }
        }
        break;
    default:
        return std::nullopt;
    }

    std::vector<SolutionBranch> combined;
    std::vector<SolverVariable> resultVariables;
    mathematics::AssumptionSet globalConditions;
    for (Expr& base : bases) {
        const Expr target = periodicTarget(
            std::move(base), period, parameter,
            builtins, mathematics, angles, assumptions);
        auto family = solvePeriodicTarget(
            matched->argument, target, variable, parameter, *remainingConditions,
            builtins, mathematics, angles);
        if (!family)
            return std::nullopt;
        if (resultVariables.empty())
            resultVariables.assign(family->variables().begin(), family->variables().end());
        for (const auto& predicate : family->conditions().predicates())
            globalConditions.add(predicate);
        combined.insert(combined.end(), family->branches().begin(), family->branches().end());
    }

    if (combined.empty())
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});
    return SolutionSet::finite(std::move(resultVariables), std::move(combined))
        .withAdditionalConditions(globalConditions);
}

std::optional<SolutionSet> solveRealInjectiveFunctionRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isEqualRelation(relation, builtins))
        return std::nullopt;

    const auto& sides = relation.asCall().arguments;

    mathematics::AssumptionSet proofAssumptions = assumptions;
    proofAssumptions.add(mathematics::elementOf(
        Expr{variable}, mathematics::NumericDomain::Real));
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, proofAssumptions};

    // principal real Logの定義域u>0では log[u] <= u-1 < u。
    // よって log[u]==u は実固定点を持たない。複素固定点を誤って除外しないよう，
    // uが実数と証明できる場合だけこの完全性を使う。
    auto logarithmFixedPoint = matchUnarySelfRelationSide(
        sides[0], sides[1], variable, mathematics, mathematics::FunctionId::Log);
    if (!logarithmFixedPoint)
        logarithmFixedPoint = matchUnarySelfRelationSide(
            sides[1], sides[0], variable, mathematics, mathematics::FunctionId::Log);
    if (logarithmFixedPoint && provablyRealForRealVariable(
            *logarithmFixedPoint, variable, builtins, mathematics, assumptions))
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    // Principal Lambert W is a one-to-one real map [-1/e,inf) -> [-1,inf).
    // W_0(u)==r therefore reduces exactly to u==r Exp[r] once r>=-1 is proved;
    // 証明できない記号的値域は推測せず，側条件として保持する。
    auto lambertSide = matchNamedFunctionSide(
        sides[0], sides[1], variable, mathematics, mathematics::FunctionId::LambertW);
    if (!lambertSide)
        lambertSide = matchNamedFunctionSide(
            sides[1], sides[0], variable, mathematics, mathematics::FunctionId::LambertW);
    if (lambertSide) {
        const TruthValue rhsReal = knowledge.prove(mathematics::elementOf(
            lambertSide->rhs, mathematics::NumericDomain::Real));
        const TruthValue rhsBelowRange = knowledge.prove(mathematics::relation(
            RelationKind::Less, lambertSide->rhs, integerExpr(-1)));
        if (rhsReal == TruthValue::False || rhsBelowRange == TruthValue::True)
            return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

        mathematics::AssumptionSet conditions;
        if (rhsReal == TruthValue::Unknown)
            conditions.add(mathematics::elementOf(
                lambertSide->rhs, mathematics::NumericDomain::Real));
        if (knowledge.prove(mathematics::relation(
                RelationKind::GreaterEqual, lambertSide->rhs, integerExpr(-1))) != TruthValue::True)
            conditions.add(mathematics::relation(
                RelationKind::GreaterEqual, lambertSide->rhs, integerExpr(-1)));
        Expr target = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Multiply), {
                lambertSide->rhs,
                Expr::call(builtins.symbol(BuiltinId::Exp), {lambertSide->rhs})}),
            builtins, mathematics, angles, proofAssumptions);
        return solveRealTargetEquation(
            lambertSide->argument, std::move(target), variable, conditions,
            builtins, mathematics, angles);
    }

    // cosh is even and strictly increasing on [0,inf).  For real target r>=1,
    // cosh(u)==r is exactly u==+/-acosh(r); at r==1 the two branches coalesce.
    auto coshSide = matchNamedFunctionSide(
        sides[0], sides[1], variable, mathematics, mathematics::FunctionId::Cosh);
    if (!coshSide)
        coshSide = matchNamedFunctionSide(
            sides[1], sides[0], variable, mathematics, mathematics::FunctionId::Cosh);
    if (coshSide) {
        const TruthValue rhsReal = knowledge.prove(mathematics::elementOf(
            coshSide->rhs, mathematics::NumericDomain::Real));
        const TruthValue belowRange = knowledge.prove(mathematics::relation(
            RelationKind::Less, coshSide->rhs, integerExpr(1)));
        if (rhsReal == TruthValue::False || belowRange == TruthValue::True)
            return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

        mathematics::AssumptionSet conditions;
        if (rhsReal == TruthValue::Unknown)
            conditions.add(mathematics::elementOf(
                coshSide->rhs, mathematics::NumericDomain::Real));
        if (knowledge.prove(mathematics::relation(
                RelationKind::GreaterEqual, coshSide->rhs, integerExpr(1))) != TruthValue::True)
            conditions.add(mathematics::relation(
                RelationKind::GreaterEqual, coshSide->rhs, integerExpr(1)));

        Expr positive = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Acosh), {coshSide->rhs}),
            builtins, mathematics, angles, proofAssumptions);
        auto plus = solveRealTargetEquation(
            coshSide->argument, positive, variable, conditions,
            builtins, mathematics, angles);
        if (!plus)
            return std::nullopt;
        if (knowledge.prove(mathematics::relation(
                RelationKind::Equal, coshSide->rhs, integerExpr(1))) == TruthValue::True)
            return plus;

        Expr negative = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Negate), {std::move(positive)}),
            builtins, mathematics, angles, proofAssumptions);
        auto minus = solveRealTargetEquation(
            coshSide->argument, std::move(negative), variable, conditions,
            builtins, mathematics, angles);
        if (!minus)
            return std::nullopt;
        return mergeFiniteTargets({std::move(*plus), std::move(*minus)}, variable);
    }

    // log[b,u] == r は b>0, b!=1, r∈Real が証明できる場合，
    // u == b^r へ安全に反転できる。log2/log10はSolve normalizationで
    // このcanonical 2引数Logへ寄るため，同じknowledge pathを使う。
    auto baseLog = matchBaseLogSide(sides[0], sides[1], variable, builtins);
    if (!baseLog)
        baseLog = matchBaseLogSide(sides[1], sides[0], variable, builtins);
    if (baseLog) {
        const CertifiedOrder baseToOne = certifiedConstantOrder(
            baseLog->base, integerExpr(1), builtins, mathematics, angles);
        if (baseToOne != CertifiedOrder::Unknown && baseToOne != CertifiedOrder::Equal
            && knowledge.prove(mathematics::relation(
                RelationKind::Greater, baseLog->base, integerExpr(0))) == TruthValue::True
            && knowledge.prove(mathematics::elementOf(
                baseLog->rhs, mathematics::NumericDomain::Real)) == TruthValue::True) {
            Expr target = simplifyForSolve(
                Expr::call(builtins.symbol(BuiltinId::Power),
                    {baseLog->base, baseLog->rhs}),
                builtins, mathematics, angles, proofAssumptions);
            Expr transformed = Expr::call(
                builtins.symbol(BuiltinId::Equal), {baseLog->argument, std::move(target)});
            return solveUnivariatePolynomialRelation(
                transformed, variable, builtins, mathematics, angles);
        }
    }

    auto matched = matchFunctionSide(sides[0], sides[1], variable, mathematics);
    if (!matched)
        matched = matchFunctionSide(sides[1], sides[0], variable, mathematics);
    if (!matched)
        return std::nullopt;

    matched->rhs = simplifyForSolve(
        matched->rhs, builtins, mathematics, angles, assumptions);

    const mathematics::AssumptionSet allRangeConditions =
        rangeConditions(*matched->definition, matched->rhs);
    const auto remainingConditions = proveRangeConditions(
        *matched->definition, matched->rhs, builtins, mathematics, assumptions);
    if (!remainingConditions)
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    const auto* inverse = mathematics.findFunction(*matched->definition->inverseFunction);
    if (!inverse)
        return std::nullopt;

    Expr inverseValue = Expr::call(inverse->symbol, {matched->rhs});
    Expr transformed = Expr::call(
        builtins.symbol(BuiltinId::Equal), {matched->argument, std::move(inverseValue)});
    SolutionSet result = solveUnivariatePolynomialRelation(
        transformed, variable, builtins, mathematics, angles);

    // 逆函数の値域条件は、変換後に現れたLog/Atanh等のdefinednessを含意する。
    // polynomial solverが安全のため収集した条件を、この追加知識で再検証して a>0 と a!=0 のような冗長条件を残さない。
    mathematics::AssumptionSet rangeProofAssumptions = assumptions;
    for (const auto& predicate : allRangeConditions.predicates())
        rangeProofAssumptions.add(predicate);
    const mathematics::KnowledgeContext proofKnowledge{
        builtins, mathematics, rangeProofAssumptions};
    mathematics::AssumptionSet filteredGlobal;
    for (const auto& predicate : result.conditions().predicates()) {
        const mathematics::Predicate simplifiedPredicate = simplifyPredicate(
            predicate, builtins, mathematics, angles, rangeProofAssumptions);
        const TruthValue truth = proofKnowledge.prove(simplifiedPredicate);
        if (truth == TruthValue::False)
            return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});
        if (truth == TruthValue::Unknown)
            filteredGlobal.add(simplifiedPredicate);
    }

    if (result.kind() == SolutionSetKind::Finite) {
        std::vector<SolutionBranch> branches(result.branches().begin(), result.branches().end());
        for (SolutionBranch& branch : branches)
            for (const auto& predicate : remainingConditions->predicates())
                branch.conditions.add(predicate);
        SolutionSet cleaned = SolutionSet::finite(
            std::vector<SolverVariable>{result.variables().begin(), result.variables().end()},
            std::move(branches));
        return cleaned.withAdditionalConditions(filteredGlobal);
    }

    SolutionSet cleaned = result.withAdditionalConditions(filteredGlobal);
    return cleaned.withAdditionalConditions(*remainingConditions);
}

} // namespace mmcal::solver
