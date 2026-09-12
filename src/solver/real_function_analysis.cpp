// 一変数実函数の定義域・単調性・値域解析
#include "real_function_analysis.hpp"

#include "expression/exact_value.hpp"

#include "mathematics/knowledge_context.hpp"
#include "polynomial_solver.hpp"
#include "radical_solver.hpp"
#include "solve_constraints.hpp"
#include "solve_normalization.hpp"
#include "solver_limits.hpp"
#include "solver_support.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/limit.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"
#include "transcendental_solver.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace mmcal::solver {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::NumericDomain;
using mathematics::RealSign;
using mathematics::RelationKind;
using mathematics::TruthValue;
using numeric::BigInt;

[[nodiscard]] Expr negativeInfinity(
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    return builtinCall(builtins, BuiltinId::Negate, {Expr{infinity}});
}



[[nodiscard]] bool unresolvedLimit(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    return builtins.isCallTo(expression, BuiltinId::Limit);
}

struct DomainRequirements final {
    enum class State { Supported, Empty, Unknown } state = State::Supported;
    std::vector<mathematics::Predicate> predicates;
    // trueなら，predicates外でも式自体は定義されるが非実値を取り得る。
    // これを持つ部分式をexp等へ無条件に合成すると，再び実値へ戻る点を
    // 見落とし得るため，real domainの完全性証明に使わない。
    bool hasNonRealContinuation = false;
};

void appendPredicate(DomainRequirements& requirements, mathematics::Predicate predicate) {
    requirements.predicates.push_back(std::move(predicate));
}

[[nodiscard]] DomainRequirements mergeRequirements(
    DomainRequirements lhs,
    DomainRequirements rhs) {
    if (lhs.state == DomainRequirements::State::Empty
        || rhs.state == DomainRequirements::State::Empty)
        return DomainRequirements{DomainRequirements::State::Empty, {}, false};
    if (lhs.state == DomainRequirements::State::Unknown
        || rhs.state == DomainRequirements::State::Unknown)
        return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
    lhs.predicates.insert(
        lhs.predicates.end(),
        std::make_move_iterator(rhs.predicates.begin()),
        std::make_move_iterator(rhs.predicates.end()));
    lhs.hasNonRealContinuation =
        lhs.hasNonRealContinuation || rhs.hasNonRealContinuation;
    return lhs;
}

[[nodiscard]] DomainRequirements collectRealDomainRequirements(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions);

[[nodiscard]] DomainRequirements collectAllArguments(
    const std::vector<Expr>& arguments,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    DomainRequirements result;
    for (const Expr& argument : arguments) {
        result = mergeRequirements(
            std::move(result),
            collectRealDomainRequirements(
                argument, variable, builtins, mathematics, assumptions));
        if (result.state != DomainRequirements::State::Supported)
            return result;
    }
    return result;
}

[[nodiscard]] DomainRequirements constantRealRequirement(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::elementOf(expression, NumericDomain::Real))
        == TruthValue::True)
        return {};
    const auto facts = knowledge.facts(expression);
    if (facts.provablyNonReal)
        return DomainRequirements{DomainRequirements::State::Empty, {}, false};
    return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
}

[[nodiscard]] DomainRequirements collectRealDomainRequirements(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    if (!symbolic::containsSymbol(expression, variable))
        return constantRealRequirement(expression, builtins, mathematics, assumptions);

    if (expression.isSymbol())
        return expression.asSymbol().sameIdentity(variable)
            ? DomainRequirements{}
            : DomainRequirements{DomainRequirements::State::Unknown, {}};
    if (expression.isNumber())
        return expression.asNumber().isReal()
            ? DomainRequirements{}
            : DomainRequirements{DomainRequirements::State::Empty, {}};
    if (!expression.isCall())
        return DomainRequirements{DomainRequirements::State::Unknown, {}, false};

    const auto* builtin = builtins.find(expression.asCall().head);
    if (!builtin)
        return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
    const auto& arguments = expression.asCall().arguments;

    switch (builtin->id) {
    case BuiltinId::Negate:
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        return collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);

    case BuiltinId::Add:
    case BuiltinId::Subtract: {
        DomainRequirements result;
        std::size_t nonRealContinuations = 0;
        for (const Expr& argument : arguments) {
            DomainRequirements current = collectRealDomainRequirements(
                argument, variable, builtins, mathematics, assumptions);
            if (current.state != DomainRequirements::State::Supported)
                return current;
            if (current.hasNonRealContinuation)
                ++nonRealContinuations;
            result = mergeRequirements(std::move(result), std::move(current));
        }
        // 非実部を持ち得る項が2つ以上あると，real domain外でのimaginary cancellationを
        // 除外できない。1項だけなら実数項の加減算では非実部を打ち消せない。
        if (nonRealContinuations > 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        return result;
    }

    case BuiltinId::Multiply: {
        DomainRequirements result = collectAllArguments(
            arguments, variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported)
            return result;
        // 非実factorに実factorのzeroが掛かると，real domain外の孤立点で0へ戻り得る。
        // そのzero setまで同時に解析しない初版では完全domainを主張しない。
        if (result.hasNonRealContinuation)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        return result;
    }

    case BuiltinId::Divide: {
        if (arguments.size() != 2)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements numerator = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        DomainRequirements denominator = collectRealDomainRequirements(
            arguments[1], variable, builtins, mathematics, assumptions);
        if (numerator.state != DomainRequirements::State::Supported)
            return numerator;
        if (denominator.state != DomainRequirements::State::Supported)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        if (denominator.hasNonRealContinuation) {
            // 実かつ非零の分子を複素分母で割った値が実になるのは，分母自身が実のときだけ。
            // したがってこの場合は分母の実domainをそのまま採用し，zeroだけ除外できる。
            // 分子が0になり得る場合は，複素側でも0へ戻る孤立点を除外できないため保守的にUnknownとする。
            if (!numerator.hasNonRealContinuation) {
                const mathematics::KnowledgeContext knowledge{
                    builtins, mathematics, assumptions};
                if (knowledge.prove(mathematics::relation(
                        RelationKind::NotEqual, arguments[0], integerExpr(0)))
                    == TruthValue::True) {
                    DomainRequirements result = mergeRequirements(
                        std::move(numerator), std::move(denominator));
                    result.hasNonRealContinuation = false;
                    appendPredicate(result, mathematics::relation(
                        RelationKind::NotEqual, arguments[1], integerExpr(0)));
                    return result;
                }
            }
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        }
        DomainRequirements result = mergeRequirements(
            std::move(numerator), std::move(denominator));
        appendPredicate(result, mathematics::relation(
            RelationKind::NotEqual, arguments[1], integerExpr(0)));
        return result;
    }

    case BuiltinId::Power: {
        if (arguments.size() != 2)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};

        DomainRequirements baseRequirements = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        DomainRequirements exponentRequirements = collectRealDomainRequirements(
            arguments[1], variable, builtins, mathematics, assumptions);
        if (baseRequirements.state != DomainRequirements::State::Supported)
            return baseRequirements;
        if (exponentRequirements.state != DomainRequirements::State::Supported)
            return exponentRequirements;

        // 部分式が実domain外でも複素値として継続し，その後Powerで実軸へ戻る点まで
        // 同時に証明する一般算法はまだ持たない。不完全なdomainを捏造しない。
        if (baseRequirements.hasNonRealContinuation
            || exponentRequirements.hasNonRealContinuation)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};

        DomainRequirements result = mergeRequirements(
            std::move(baseRequirements), std::move(exponentRequirements));
        const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
        const mathematics::ValueFacts baseFacts = knowledge.facts(arguments[0]);
        const mathematics::ValueFacts exponentFacts = knowledge.facts(arguments[1]);

        const auto signIsNonZero = [](RealSign sign) {
            return sign == RealSign::Positive || sign == RealSign::Negative
                || sign == RealSign::NonZero;
        };
        const auto signIsNonPositive = [](RealSign sign) {
            return sign == RealSign::Negative || sign == RealSign::Zero
                || sign == RealSign::NonPositive;
        };

        const TruthValue exponentIntegerKnowledge = knowledge.prove(
            mathematics::elementOf(arguments[1], NumericDomain::Integer));
        const bool exponentInteger = exponentIntegerKnowledge == TruthValue::True;
        const bool exponentNonInteger = exponentIntegerKnowledge == TruthValue::False;

        // 実整数指数なら負の底も常に実数。0だけは指数の符号で定義性が変わる。
        if (exponentInteger) {
            if (exponentFacts.sign == RealSign::Positive)
                return result;
            if (signIsNonZero(baseFacts.sign))
                return result;
            if (signIsNonPositive(exponentFacts.sign)) {
                appendPredicate(result, mathematics::relation(
                    RelationKind::NotEqual, arguments[0], integerExpr(0)));
                return result;
            }
            // integer-valuedだが符号未知，かつbaseが0を取り得る場合は
            // 0^positiveだけを残す相関条件が必要になるため，現表現では完全化しない。
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        }

        // 正の実数底ならprincipal Logが実数なので，任意の実指数を安全に扱える。
        if (baseFacts.sign == RealSign::Positive)
            return result;

        // identically zeroな底は0^qの規則だけで完全に記述できる。
        if (baseFacts.sign == RealSign::Zero) {
            appendPredicate(result, mathematics::relation(
                RelationKind::Greater, arguments[1], integerExpr(0)));
            return result;
        }

        // 非負底まで証明できる場合は，指数の符号だけでzeroを含めるか決まる。
        if (baseFacts.sign == RealSign::NonNegative) {
            if (exponentFacts.sign == RealSign::Positive)
                return result;
            if (signIsNonPositive(exponentFacts.sign)) {
                appendPredicate(result, mathematics::relation(
                    RelationKind::Greater, arguments[0], integerExpr(0)));
                return result;
            }
            if (exponentNonInteger
                && exponentFacts.sign == RealSign::NonNegative)
                return result; // 非整数かつ>=0なら0ではないので実際には正。
            if (arguments[0] == arguments[1]
                && exponentFacts.sign == RealSign::NonNegative) {
                // x^x型ではbase=0なら指数も必ず0なので，0^0だけを除けばよい。
                appendPredicate(result, mathematics::relation(
                    RelationKind::Greater, arguments[0], integerExpr(0)));
                return result;
            }
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        }

        // principal Power(b,q)=Exp(q principal Log(b))では，実非整数qに対し
        // b<0は複素値になる。q>0ならb=0を含み，q<0ならstrict positiveだけ。
        // ここでprovablyNonIntegerを要求するのは，負側の孤立した整数指数点を
        // 誤って落とさないためである。
        if (exponentNonInteger) {
            if (exponentFacts.sign == RealSign::Positive
                || exponentFacts.sign == RealSign::NonNegative) {
                appendPredicate(result, mathematics::relation(
                    RelationKind::GreaterEqual, arguments[0], integerExpr(0)));
                result.hasNonRealContinuation = true;
                return result;
            }
            if (exponentFacts.sign == RealSign::Negative
                || exponentFacts.sign == RealSign::NonPositive) {
                appendPredicate(result, mathematics::relation(
                    RelationKind::Greater, arguments[0], integerExpr(0)));
                result.hasNonRealContinuation = true;
                return result;
            }
        }

        // x^x等は負側に離散的な実値点を持つ。intervalだけで完全domainを
        // 表せない現段階では，正側だけを描いて「完全」とは主張しない。
        return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
    }

    case BuiltinId::Sqrt: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported)
            return result;
        appendPredicate(result, mathematics::relation(
            RelationKind::GreaterEqual, arguments[0], integerExpr(0)));
        result.hasNonRealContinuation = true;
        return result;
    }

    case BuiltinId::Log:
    case BuiltinId::Log2:
    case BuiltinId::Log10: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported)
            return result;
        appendPredicate(result, mathematics::relation(
            RelationKind::Greater, arguments[0], integerExpr(0)));
        result.hasNonRealContinuation = true;
        return result;
    }

    case BuiltinId::Log1p: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported)
            return result;
        appendPredicate(result, mathematics::relation(
            RelationKind::Greater, arguments[0], integerExpr(-1)));
        result.hasNonRealContinuation = true;
        return result;
    }

    case BuiltinId::ExponentialIntegralEi: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported
            || result.hasNonRealContinuation)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        // Principal Ei is real on both real half-axes and singular only at zero.
        appendPredicate(result, mathematics::relation(
            RelationKind::NotEqual, arguments[0], integerExpr(0)));
        return result;
    }

    case BuiltinId::CosineIntegralCi: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported)
            return result;
        // Principal Ci(x) is real only for x>0. On x<0 it acquires the
        // principal-log branch offset +i Pi, so a real Plot must not draw it.
        appendPredicate(result, mathematics::relation(
            RelationKind::Greater, arguments[0], integerExpr(0)));
        result.hasNonRealContinuation = true;
        return result;
    }

    case BuiltinId::Asin:
    case BuiltinId::Acos: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported)
            return result;
        appendPredicate(result, mathematics::relation(
            RelationKind::GreaterEqual, arguments[0], integerExpr(-1)));
        appendPredicate(result, mathematics::relation(
            RelationKind::LessEqual, arguments[0], integerExpr(1)));
        result.hasNonRealContinuation = true;
        return result;
    }

    case BuiltinId::Acosh: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported)
            return result;
        appendPredicate(result, mathematics::relation(
            RelationKind::GreaterEqual, arguments[0], integerExpr(1)));
        result.hasNonRealContinuation = true;
        return result;
    }

    case BuiltinId::Atanh: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported)
            return result;
        appendPredicate(result, mathematics::relation(
            RelationKind::Greater, arguments[0], integerExpr(-1)));
        appendPredicate(result, mathematics::relation(
            RelationKind::Less, arguments[0], integerExpr(1)));
        result.hasNonRealContinuation = true;
        return result;
    }

    case BuiltinId::Cbrt: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported)
            return result;
        // cbrtはmmCalではReal専用。非実argumentでは未定義なのでnon-real continuationを持たない。
        result.hasNonRealContinuation = false;
        return result;
    }

    case BuiltinId::Floor:
    case BuiltinId::Ceil:
    case BuiltinId::Sign: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported
            || result.hasNonRealContinuation)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        return result;
    }

    case BuiltinId::Exp:
    case BuiltinId::Expm1:
    case BuiltinId::Sin:
    case BuiltinId::Cos:
    case BuiltinId::Sinh:
    case BuiltinId::Cosh:
    case BuiltinId::Tanh:
    case BuiltinId::Atan:
    case BuiltinId::Asinh:
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS: {
        if (arguments.size() != 1)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported
            || result.hasNonRealContinuation)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        return result;
    }

    // tan/sec/cot/csc等は周期的に無限個のpoleを持つ。有限interval列を返す初版では
    // 不完全なdomainを捏造せず，周期domain analyzerの後続拡張へ残す。
    case BuiltinId::Tan:
    case BuiltinId::Cot:
    case BuiltinId::Sec:
    case BuiltinId::Csc:
    case BuiltinId::Csch:
    case BuiltinId::Coth:
    case BuiltinId::Gamma:
    case BuiltinId::LogGamma:
    case BuiltinId::LogarithmicIntegralLi:
    case BuiltinId::Polylog:
    case BuiltinId::Hypergeometric1F1:
    case BuiltinId::Hypergeometric2F1:
    case BuiltinId::EllipticF:
    case BuiltinId::EllipticE:
    case BuiltinId::EllipticPi:
        return DomainRequirements{DomainRequirements::State::Unknown, {}, false};

    default:
        break;
    }

    // Registry上RealToRealかComplexToComplexRealPreservingでeverywhereなら，
    // 実引数の実domainをそのまま持ち上げられる。特殊definednessは上で明示したものだけ扱う。
    const auto* definition = mathematics.findFunction(expression.asCall().head);
    if (definition && arguments.size() == 1
        && definition->definednessRule == mathematics::FunctionDefinednessRule::Everywhere
        && (definition->domainRule == mathematics::FunctionDomainRule::RealToReal
            || definition->domainRule
                == mathematics::FunctionDomainRule::ComplexToComplexRealPreserving))
    {
        DomainRequirements result = collectRealDomainRequirements(
            arguments[0], variable, builtins, mathematics, assumptions);
        if (result.state != DomainRequirements::State::Supported
            || result.hasNonRealContinuation)
            return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
        return result;
    }

    return DomainRequirements{DomainRequirements::State::Unknown, {}, false};
}

[[nodiscard]] std::optional<int> compareFinite(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    if (lhs == rhs)
        return 0;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::relation(RelationKind::Less, lhs, rhs))
        == TruthValue::True)
        return -1;
    if (knowledge.prove(mathematics::relation(RelationKind::Greater, lhs, rhs))
        == TruthValue::True)
        return 1;

    // 同じglobal monotone函数の値同士なら，引数の順序へexactに戻す。
    // asinh[-1] < asinh[1]のようなdomain境界を数値近似へ落とさず比較できる。
    if (lhs.isCall() && rhs.isCall()
        && lhs.asCall().head.sameIdentity(rhs.asCall().head)
        && lhs.asCall().arguments.size() == 1
        && rhs.asCall().arguments.size() == 1) {
        if (const auto* function = mathematics.findFunction(lhs.asCall().head);
            function && function->realGloballyInjective
            && function->realMonotonicity != mathematics::RealMonotonicity::Unknown) {
            const auto inner = compareFinite(
                lhs.asCall().arguments.front(), rhs.asCall().arguments.front(),
                builtins, mathematics, assumptions);
            if (!inner)
                return std::nullopt;
            return function->realMonotonicity == mathematics::RealMonotonicity::Decreasing
                ? -*inner : *inner;
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<RealDomainInterval> intersectIntervals(
    const RealDomainInterval& lhs,
    const RealDomainInterval& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    RealDomainInterval result;

    if (!lhs.lower) {
        result.lower = rhs.lower;
        result.lowerInclusive = rhs.lowerInclusive;
    }
    else if (!rhs.lower) {
        result.lower = lhs.lower;
        result.lowerInclusive = lhs.lowerInclusive;
    }
    else {
        const auto comparison = compareFinite(
            *lhs.lower, *rhs.lower, builtins, mathematics, assumptions);
        if (!comparison)
            return std::nullopt;
        if (*comparison > 0) {
            result.lower = lhs.lower;
            result.lowerInclusive = lhs.lowerInclusive;
        }
        else if (*comparison < 0) {
            result.lower = rhs.lower;
            result.lowerInclusive = rhs.lowerInclusive;
        }
        else {
            result.lower = lhs.lower;
            result.lowerInclusive = lhs.lowerInclusive && rhs.lowerInclusive;
        }
    }

    if (!lhs.upper) {
        result.upper = rhs.upper;
        result.upperInclusive = rhs.upperInclusive;
    }
    else if (!rhs.upper) {
        result.upper = lhs.upper;
        result.upperInclusive = lhs.upperInclusive;
    }
    else {
        const auto comparison = compareFinite(
            *lhs.upper, *rhs.upper, builtins, mathematics, assumptions);
        if (!comparison)
            return std::nullopt;
        if (*comparison < 0) {
            result.upper = lhs.upper;
            result.upperInclusive = lhs.upperInclusive;
        }
        else if (*comparison > 0) {
            result.upper = rhs.upper;
            result.upperInclusive = rhs.upperInclusive;
        }
        else {
            result.upper = lhs.upper;
            result.upperInclusive = lhs.upperInclusive && rhs.upperInclusive;
        }
    }

    if (result.lower && result.upper) {
        const auto comparison = compareFinite(
            *result.lower, *result.upper, builtins, mathematics, assumptions);
        if (!comparison)
            return std::nullopt;
        if (*comparison > 0)
            return RealDomainInterval{integerExpr(0), false, integerExpr(0), false};
        if (*comparison == 0 && !(result.lowerInclusive && result.upperInclusive))
            return RealDomainInterval{integerExpr(0), false, integerExpr(0), false};
    }
    return result;
}

[[nodiscard]] bool intervalIsEmpty(const RealDomainInterval& interval) {
    return interval.lower && interval.upper
        && *interval.lower == integerExpr(0) && *interval.upper == integerExpr(0)
        && !interval.lowerInclusive && !interval.upperInclusive;
}

[[nodiscard]] std::optional<std::vector<RealDomainInterval>> intersectIntervalSets(
    const std::vector<RealDomainInterval>& lhs,
    const std::vector<RealDomainInterval>& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    std::vector<RealDomainInterval> result;
    for (const RealDomainInterval& left : lhs) {
        for (const RealDomainInterval& right : rhs) {
            const auto intersection = intersectIntervals(
                left, right, builtins, mathematics, assumptions);
            if (!intersection)
                return std::nullopt;
            if (!intervalIsEmpty(*intersection))
                result.push_back(*intersection);
            if (result.size() > limits::realDomainPieces)
                return std::nullopt;
        }
    }
    return result;
}

[[nodiscard]] bool parseBoundPredicate(
    const mathematics::RelationPredicate& predicate,
    const expression::Symbol& variable,
    RealDomainInterval& interval) {
    RelationKind relation = predicate.relation;
    Expr bound = predicate.rhs;
    if (!(predicate.lhs.isSymbol()
            && predicate.lhs.asSymbol().sameIdentity(variable))) {
        if (!(predicate.rhs.isSymbol()
                && predicate.rhs.asSymbol().sameIdentity(variable)))
            return false;
        relation = reversedRelation(relation);
        bound = predicate.lhs;
    }

    switch (relation) {
    case RelationKind::Greater:
        interval.lower = std::move(bound);
        interval.lowerInclusive = false;
        return true;
    case RelationKind::GreaterEqual:
        interval.lower = std::move(bound);
        interval.lowerInclusive = true;
        return true;
    case RelationKind::Less:
        interval.upper = std::move(bound);
        interval.upperInclusive = false;
        return true;
    case RelationKind::LessEqual:
        interval.upper = std::move(bound);
        interval.upperInclusive = true;
        return true;
    case RelationKind::Equal:
        interval.lower = bound;
        interval.lowerInclusive = true;
        interval.upper = std::move(bound);
        interval.upperInclusive = true;
        return true;
    case RelationKind::NotEqual:
        return false;
    }
    return false;
}

[[nodiscard]] std::optional<std::vector<RealDomainInterval>> intervalsFromSolutionSet(
    const SolutionSet& solutions,
    const expression::Symbol& variable) {
    if (solutions.kind() == SolutionSetKind::Empty)
        return std::vector<RealDomainInterval>{};
    if (solutions.kind() == SolutionSetKind::Universal)
        return std::vector<RealDomainInterval>{RealDomainInterval{}};
    if (solutions.kind() != SolutionSetKind::Finite)
        return std::nullopt;

    std::vector<RealDomainInterval> result;
    for (const SolutionBranch& branch : solutions.branches()) {
        if (branch.bindings.size() == 1
            && branch.bindings.front().variable.sameIdentity(variable)) {
            const Expr point = branch.bindings.front().value;
            result.push_back(RealDomainInterval{point, true, point, true});
            continue;
        }
        if (!branch.bindings.empty())
            return std::nullopt;

        RealDomainInterval interval;
        for (const mathematics::Predicate& predicate : branch.conditions.predicates()) {
            const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate);
            if (!relation || !parseBoundPredicate(*relation, variable, interval))
                return std::nullopt;
        }
        result.push_back(std::move(interval));
    }
    return result;
}

[[nodiscard]] SolutionSet solveDomainRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    SolveConstraints realConstraint;
    realConstraint.domain = NumericDomain::Real;
    const auto constrain = [&](SolutionSet solutions) {
        return applySolveConstraints(
            std::move(solutions), realConstraint, builtins, mathematics, angles);
    };

    // 既存のdomain解析は多項式を最優先にする。ここで解けないrelationだけ，
    // solver側で既にexactな実反転を持つ限定familyへfallbackする。
    SolutionSet polynomial = constrain(solveUnivariatePolynomialRelation(
        relation, variable, builtins, mathematics, angles));
    if (polynomial.kind() != SolutionSetKind::Unresolved)
        return polynomial;

    if (auto radical = solveRadicalRelation(
            relation, variable, builtins, mathematics, angles, assumptions)) {
        SolutionSet constrained = constrain(std::move(*radical));
        if (constrained.kind() != SolutionSetKind::Unresolved)
            return constrained;
    }
    if (auto exponential = solveRealExponentialRelation(
            relation, variable, builtins, mathematics, angles, assumptions)) {
        SolutionSet constrained = constrain(std::move(*exponential));
        if (constrained.kind() != SolutionSetKind::Unresolved)
            return constrained;
    }
    if (auto injective = solveRealInjectiveFunctionRelation(
            relation, variable, builtins, mathematics, angles, assumptions)) {
        SolutionSet constrained = constrain(std::move(*injective));
        if (constrained.kind() != SolutionSetKind::Unresolved)
            return constrained;
    }
    return polynomial;
}

[[nodiscard]] std::optional<std::vector<Expr>> finiteRealRoots(
    Expr equality,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    SolutionSet roots = solveDomainRelation(
        equality, variable, builtins, mathematics, angles, assumptions);
    if (roots.kind() == SolutionSetKind::Empty)
        return std::vector<Expr>{};
    if (roots.kind() != SolutionSetKind::Finite)
        return std::nullopt;

    std::vector<Expr> result;
    for (const SolutionBranch& branch : roots.branches()) {
        if (!branch.unconditional() || branch.bindings.size() != 1
            || !branch.bindings.front().variable.sameIdentity(variable))
            return std::nullopt;
        result.push_back(branch.bindings.front().value);
    }
    return result;
}

[[nodiscard]] std::optional<std::vector<RealDomainInterval>> intervalsForNotEqual(
    const mathematics::RelationPredicate& predicate,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr equality = normalizeForSolve(
        relationExpr(RelationKind::Equal, predicate.lhs, predicate.rhs, builtins),
        builtins, mathematics, angles, assumptions);
    auto roots = finiteRealRoots(
        std::move(equality), variable, builtins, mathematics, angles, assumptions);
    if (!roots)
        return std::nullopt;
    if (roots->empty())
        return std::vector<RealDomainInterval>{RealDomainInterval{}};
    if (roots->size() > limits::realDomainPieces - 1)
        return std::nullopt;

    std::stable_sort(roots->begin(), roots->end(), [&](const Expr& lhs, const Expr& rhs) {
        const auto comparison = compareFinite(
            lhs, rhs, builtins, mathematics, assumptions);
        return comparison && *comparison < 0;
    });
    for (std::size_t i = 1; i < roots->size(); ++i) {
        const auto comparison = compareFinite(
            (*roots)[i - 1], (*roots)[i], builtins, mathematics, assumptions);
        if (!comparison)
            return std::nullopt;
    }

    std::vector<RealDomainInterval> result;
    result.reserve(roots->size() + 1);
    for (std::size_t i = 0; i <= roots->size(); ++i) {
        RealDomainInterval interval;
        if (i != 0) {
            interval.lower = (*roots)[i - 1];
            interval.lowerInclusive = false;
        }
        if (i != roots->size()) {
            interval.upper = (*roots)[i];
            interval.upperInclusive = false;
        }
        result.push_back(std::move(interval));
    }
    return result;
}


[[nodiscard]] bool exactRationalEquals(const Expr& expression, std::int64_t value) {
    const auto rational = expression::exact::realRational(expression);
    return rational && *rational == numeric::Rational{BigInt{value}};
}

[[nodiscard]] std::optional<TruthValue> relationTruthFromRealRange(
    const mathematics::FunctionDefinition& definition,
    RelationKind relation,
    const Expr& rhs) {
    const auto target = expression::exact::realRational(rhs);
    if (!target)
        return std::nullopt;

    struct Bound final {
        std::optional<numeric::Rational> lower;
        bool lowerInclusive = false;
        std::optional<numeric::Rational> upper;
        bool upperInclusive = false;
    } range;

    switch (definition.realRangeRule) {
    case mathematics::RealRangeRule::Unknown:
        return std::nullopt;
    case mathematics::RealRangeRule::AllReal:
        return std::nullopt;
    case mathematics::RealRangeRule::Positive:
        range.lower = numeric::Rational{BigInt{0}};
        range.lowerInclusive = false;
        break;
    case mathematics::RealRangeRule::NonNegative:
        range.lower = numeric::Rational{BigInt{0}};
        range.lowerInclusive = true;
        break;
    case mathematics::RealRangeRule::OpenMinusOneToOne:
        range.lower = numeric::Rational{BigInt{-1}};
        range.upper = numeric::Rational{BigInt{1}};
        break;
    case mathematics::RealRangeRule::OpenZeroToTwo:
        range.lower = numeric::Rational{BigInt{0}};
        range.upper = numeric::Rational{BigInt{2}};
        break;
    case mathematics::RealRangeRule::ClosedMinusOneToOne:
        range.lower = numeric::Rational{BigInt{-1}};
        range.lowerInclusive = true;
        range.upper = numeric::Rational{BigInt{1}};
        range.upperInclusive = true;
        break;
    case mathematics::RealRangeRule::OneToInfinity:
        range.lower = numeric::Rational{BigInt{1}};
        range.lowerInclusive = true;
        break;
    }

    const auto belowLower = [&] {
        return range.lower && (*target < *range.lower
            || (*target == *range.lower && !range.lowerInclusive));
    };
    const auto aboveUpper = [&] {
        return range.upper && (*target > *range.upper
            || (*target == *range.upper && !range.upperInclusive));
    };
    const auto allGreater = [&] {
        return range.lower && (*range.lower > *target
            || (*range.lower == *target && !range.lowerInclusive));
    };
    const auto allGreaterEqual = [&] {
        return range.lower && *range.lower >= *target;
    };
    const auto allLess = [&] {
        return range.upper && (*range.upper < *target
            || (*range.upper == *target && !range.upperInclusive));
    };
    const auto allLessEqual = [&] {
        return range.upper && *range.upper <= *target;
    };
    const auto noneGreater = [&] {
        return range.upper && *range.upper <= *target;
    };
    const auto noneGreaterEqual = [&] {
        return range.upper && (*range.upper < *target
            || (*range.upper == *target && !range.upperInclusive));
    };
    const auto noneLess = [&] {
        return range.lower && *range.lower >= *target;
    };
    const auto noneLessEqual = [&] {
        return range.lower && (*range.lower > *target
            || (*range.lower == *target && !range.lowerInclusive));
    };

    switch (relation) {
    case RelationKind::Greater:
        if (allGreater()) return TruthValue::True;
        if (noneGreater()) return TruthValue::False;
        break;
    case RelationKind::GreaterEqual:
        if (allGreaterEqual()) return TruthValue::True;
        if (noneGreaterEqual()) return TruthValue::False;
        break;
    case RelationKind::Less:
        if (allLess()) return TruthValue::True;
        if (noneLess()) return TruthValue::False;
        break;
    case RelationKind::LessEqual:
        if (allLessEqual()) return TruthValue::True;
        if (noneLessEqual()) return TruthValue::False;
        break;
    case RelationKind::Equal:
        if (belowLower() || aboveUpper()) return TruthValue::False;
        break;
    case RelationKind::NotEqual:
        if (belowLower() || aboveUpper()) return TruthValue::True;
        break;
    }
    return std::nullopt;
}

[[nodiscard]] bool exactRationalInsideRealRange(
    const mathematics::FunctionDefinition& definition,
    const numeric::Rational& value) {
    const numeric::Rational zero{BigInt{0}};
    const numeric::Rational one{BigInt{1}};
    const numeric::Rational minusOne{BigInt{-1}};
    const numeric::Rational two{BigInt{2}};

    switch (definition.realRangeRule) {
    case mathematics::RealRangeRule::Unknown:
        return false;
    case mathematics::RealRangeRule::AllReal:
        return true;
    case mathematics::RealRangeRule::Positive:
        return value > zero;
    case mathematics::RealRangeRule::NonNegative:
        return value >= zero;
    case mathematics::RealRangeRule::OpenMinusOneToOne:
        return value > minusOne && value < one;
    case mathematics::RealRangeRule::OpenZeroToTwo:
        return value > zero && value < two;
    case mathematics::RealRangeRule::ClosedMinusOneToOne:
        return value >= minusOne && value <= one;
    case mathematics::RealRangeRule::OneToInfinity:
        return value >= one;
    }
    return false;
}

[[nodiscard]] std::optional<mathematics::RelationPredicate> reduceKnownDomainRelation(
    const mathematics::RelationPredicate& predicate,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    mathematics::RelationPredicate reduced = predicate;

    const auto constantTruth = [&](bool truth) {
        return mathematics::RelationPredicate{
            RelationKind::Equal, integerExpr(0), integerExpr(truth ? 0 : 1)};
    };

    // c/f(x) R t を，divisionのdefinedness（f!=0）と組み合わせて分母側へ戻す。
    // 特にsqrt[1/log[x]]やasin[1/sqrt[x]]のような入子で，
    // 一般rational-transcendental inequality solverを起動せずexactにdomainを絞れる。
    if (builtins.isCallTo(reduced.lhs, BuiltinId::Divide)
        && reduced.lhs.asCall().arguments.size() == 2) {
        const Expr& numeratorExpr = reduced.lhs.asCall().arguments[0];
        const Expr& denominator = reduced.lhs.asCall().arguments[1];
        const auto numerator = expression::exact::realRational(numeratorExpr);
        const auto target = expression::exact::realRational(reduced.rhs);
        if (numerator && target && !numerator->isZero()) {
            const bool numeratorPositive = numerator->numerator().isPositive();

            if (target->isZero()) {
                switch (reduced.relation) {
                case RelationKind::Greater:
                case RelationKind::GreaterEqual:
                    return mathematics::RelationPredicate{
                        numeratorPositive ? RelationKind::Greater : RelationKind::Less,
                        denominator, integerExpr(0)};
                case RelationKind::Less:
                case RelationKind::LessEqual:
                    return mathematics::RelationPredicate{
                        numeratorPositive ? RelationKind::Less : RelationKind::Greater,
                        denominator, integerExpr(0)};
                case RelationKind::Equal:
                    return constantTruth(false);
                case RelationKind::NotEqual:
                    return constantTruth(true);
                }
            }

            bool denominatorPositiveWhenDefined = false;
            if (denominator.isCall() && denominator.asCall().arguments.size() == 1) {
                if (const auto* function = mathematics.findFunction(denominator.asCall().head)) {
                    switch (function->realRangeRule) {
                    case mathematics::RealRangeRule::Positive:
                    case mathematics::RealRangeRule::NonNegative:
                    case mathematics::RealRangeRule::OpenZeroToTwo:
                    case mathematics::RealRangeRule::OneToInfinity:
                        denominatorPositiveWhenDefined = true;
                        break;
                    default:
                        break;
                    }
                }
            }

            if (denominatorPositiveWhenDefined) {
                const bool targetPositive = target->numerator().isPositive();
                const bool targetNegative = target->numerator().isNegative();
                if (numeratorPositive && targetNegative) {
                    switch (reduced.relation) {
                    case RelationKind::Greater:
                    case RelationKind::GreaterEqual:
                    case RelationKind::NotEqual:
                        return constantTruth(true);
                    case RelationKind::Less:
                    case RelationKind::LessEqual:
                    case RelationKind::Equal:
                        return constantTruth(false);
                    }
                }
                if (!numeratorPositive && targetPositive) {
                    switch (reduced.relation) {
                    case RelationKind::Less:
                    case RelationKind::LessEqual:
                    case RelationKind::NotEqual:
                        return constantTruth(true);
                    case RelationKind::Greater:
                    case RelationKind::GreaterEqual:
                    case RelationKind::Equal:
                        return constantTruth(false);
                    }
                }

                if ((numeratorPositive && targetPositive)
                    || (!numeratorPositive && targetNegative)) {
                    const numeric::Rational boundary = *numerator / *target;
                    RelationKind denominatorRelation = reduced.relation;
                    switch (reduced.relation) {
                    case RelationKind::Greater: denominatorRelation = RelationKind::Less; break;
                    case RelationKind::GreaterEqual: denominatorRelation = RelationKind::LessEqual; break;
                    case RelationKind::Less: denominatorRelation = RelationKind::Greater; break;
                    case RelationKind::LessEqual: denominatorRelation = RelationKind::GreaterEqual; break;
                    case RelationKind::Equal: break;
                    case RelationKind::NotEqual: break;
                    }
                    return mathematics::RelationPredicate{
                        denominatorRelation, denominator, Expr{numeric::Number{boundary}}};
                }
            }
        }
    }

    // f-g R 0 は f R g へ戻す。solverの個別familyはこのrelation形を契約にしている。
    if (exactRationalEquals(reduced.rhs, 0)
        && builtins.isCallTo(reduced.lhs, BuiltinId::Subtract)
        && reduced.lhs.asCall().arguments.size() == 2) {
        reduced.rhs = reduced.lhs.asCall().arguments[1];
        reduced.lhs = reduced.lhs.asCall().arguments[0];
        return reduced;
    }
    if (exactRationalEquals(reduced.lhs, 0)
        && builtins.isCallTo(reduced.rhs, BuiltinId::Subtract)
        && reduced.rhs.asCall().arguments.size() == 2) {
        reduced.lhs = reduced.rhs.asCall().arguments[1];
        reduced.rhs = reduced.rhs.asCall().arguments[0];
        return reduced;
    }

    if (!reduced.lhs.isCall() || reduced.lhs.asCall().arguments.size() != 1
        || symbolic::containsSymbol(reduced.rhs, variable))
        return std::nullopt;
    const auto* definition = builtins.find(reduced.lhs.asCall().head);
    if (!definition)
        return std::nullopt;
    const Expr& argument = reduced.lhs.asCall().arguments.front();

    std::optional<Expr> target;
    RelationKind relation = reduced.relation;
    switch (definition->id) {
    case BuiltinId::Log:
    case BuiltinId::Log2:
    case BuiltinId::Log10:
        if (exactRationalEquals(reduced.rhs, 0)) {
            target = integerExpr(1);
        }
        break;
    case BuiltinId::Cbrt: {
        const auto rhs = expression::exact::realRational(reduced.rhs);
        if (rhs)
            target = Expr{numeric::Number{*rhs * *rhs * *rhs}};
        break;
    }
    case BuiltinId::Log1p:
    case BuiltinId::Expm1:
    case BuiltinId::Sinh:
    case BuiltinId::Tanh:
    case BuiltinId::Asinh:
    case BuiltinId::Atanh:
    case BuiltinId::Asin:
    case BuiltinId::Atan:
    case BuiltinId::Erf:
        if (exactRationalEquals(reduced.rhs, 0))
            target = integerExpr(0);
        break;
    case BuiltinId::Sqrt: {
        const auto rhs = expression::exact::realRational(reduced.rhs);
        if (rhs && !rhs->numerator().isNegative()) {
            // sqrt[g]は実domainで非負かつ単調なので，非負exact境界との比較は
            // squaringしてg側へexactに戻せる。domain条件g>=0は別predicateで保持される。
            target = Expr{numeric::Number{*rhs * *rhs}};
        }
        break;
    }
    case BuiltinId::Exp:
        if (exactRationalEquals(reduced.rhs, 1)) {
            target = integerExpr(0);
        }
        break;
    case BuiltinId::Acos:
        if (exactRationalEquals(reduced.rhs, 0)) {
            target = integerExpr(1);
            relation = reversedRelation(relation);
        }
        break;
    case BuiltinId::Acosh:
        if (exactRationalEquals(reduced.rhs, 0)) {
            target = integerExpr(1);
        }
        break;
    case BuiltinId::Erfc:
        if (exactRationalEquals(reduced.rhs, 1)) {
            target = integerExpr(0);
            relation = reversedRelation(relation);
        }
        break;
    default:
        break;
    }

    // MathRegistryで実軸上のglobal injectivity・単調性・inverseが証明済みなら，
    // exact rational境界をinverse側へ戻す。周期函数やbranchごとのinverseには使わない。
    if (!target) {
        const auto* function = mathematics.findFunction(reduced.lhs.asCall().head);
        const auto rhs = expression::exact::realRational(reduced.rhs);
        if (function && rhs && function->realGloballyInjective
            && function->inverseFunction
            && function->realMonotonicity != mathematics::RealMonotonicity::Unknown
            && exactRationalInsideRealRange(*function, *rhs)) {
            const auto* inverse = mathematics.findFunction(*function->inverseFunction);
            if (inverse) {
                target = Expr::call(inverse->symbol, {reduced.rhs});
                if (function->realMonotonicity == mathematics::RealMonotonicity::Decreasing)
                    relation = reversedRelation(relation);
            }
        }
    }

    if (!target)
        return std::nullopt;
    return mathematics::RelationPredicate{relation, argument, std::move(*target)};
}

[[nodiscard]] std::optional<std::vector<RealDomainInterval>> intervalsForPredicate(
    const mathematics::Predicate& predicate,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate);
    if (!relation) {
        const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
        if (knowledge.prove(predicate) == TruthValue::True)
            return std::vector<RealDomainInterval>{RealDomainInterval{}};
        if (knowledge.prove(predicate) == TruthValue::False)
            return std::vector<RealDomainInterval>{};
        return std::nullopt;
    }

    if (relation->lhs.isCall() && relation->lhs.asCall().arguments.size() == 1) {
        if (const auto* definition = mathematics.findFunction(relation->lhs.asCall().head)) {
            if (const auto truth = relationTruthFromRealRange(
                    *definition, relation->relation, relation->rhs)) {
                if (*truth == TruthValue::True)
                    return std::vector<RealDomainInterval>{RealDomainInterval{}};
                if (*truth == TruthValue::False)
                    return std::vector<RealDomainInterval>{};
            }
        }
    }

    if (const auto reduced = reduceKnownDomainRelation(
            *relation, variable, builtins, mathematics);
        reduced && !(*reduced == *relation)) {
        return intervalsForPredicate(
            mathematics::Predicate{*reduced}, variable,
            builtins, mathematics, angles, assumptions);
    }

    auto directVariableBound = [&](
        const Expr& variableSide, const Expr& bound, RelationKind relationKind)
        -> std::optional<std::vector<RealDomainInterval>> {
        if (!variableSide.isSymbol() || !variableSide.asSymbol().sameIdentity(variable)
            || symbolic::containsSymbol(bound, variable))
            return std::nullopt;
        const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
        if (knowledge.prove(mathematics::elementOf(bound, NumericDomain::Real))
            != TruthValue::True)
            return std::nullopt;

        switch (relationKind) {
        case RelationKind::Greater:
            return std::vector<RealDomainInterval>{
                RealDomainInterval{bound, false, std::nullopt, false}};
        case RelationKind::GreaterEqual:
            return std::vector<RealDomainInterval>{
                RealDomainInterval{bound, true, std::nullopt, false}};
        case RelationKind::Less:
            return std::vector<RealDomainInterval>{
                RealDomainInterval{std::nullopt, false, bound, false}};
        case RelationKind::LessEqual:
            return std::vector<RealDomainInterval>{
                RealDomainInterval{std::nullopt, false, bound, true}};
        case RelationKind::Equal:
            return std::vector<RealDomainInterval>{
                RealDomainInterval{bound, true, bound, true}};
        case RelationKind::NotEqual:
            return std::vector<RealDomainInterval>{
                RealDomainInterval{std::nullopt, false, bound, false},
                RealDomainInterval{bound, false, std::nullopt, false}};
        }
        return std::nullopt;
    };

    if (auto direct = directVariableBound(
            relation->lhs, relation->rhs, relation->relation))
        return direct;
    if (auto direct = directVariableBound(
            relation->rhs, relation->lhs, reversedRelation(relation->relation)))
        return direct;

    if (!symbolic::containsSymbol(relation->lhs, variable)
        && !symbolic::containsSymbol(relation->rhs, variable)) {
        const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
        if (knowledge.prove(predicate) == TruthValue::True)
            return std::vector<RealDomainInterval>{RealDomainInterval{}};
        if (knowledge.prove(predicate) == TruthValue::False)
            return std::vector<RealDomainInterval>{};
        return std::nullopt;
    }

    if (relation->relation == RelationKind::NotEqual)
        return intervalsForNotEqual(
            *relation, variable, builtins, mathematics, angles, assumptions);

    Expr relationExpression = normalizeForSolve(
        relationExpr(relation->relation, relation->lhs, relation->rhs, builtins),
        builtins, mathematics, angles, assumptions);
    SolutionSet solutions = solveDomainRelation(
        relationExpression, variable, builtins, mathematics, angles, assumptions);
    return intervalsFromSolutionSet(solutions, variable);
}

[[nodiscard]] std::optional<std::vector<RealDomainInterval>> realDomainIntervals(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    DomainRequirements requirements = collectRealDomainRequirements(
        expression, variable, builtins, mathematics, assumptions);
    if (requirements.state == DomainRequirements::State::Empty)
        return std::vector<RealDomainInterval>{};
    if (requirements.state == DomainRequirements::State::Unknown)
        return std::nullopt;

    std::vector<RealDomainInterval> domain{RealDomainInterval{}};
    for (const mathematics::Predicate& predicate : requirements.predicates) {
        auto constraint = intervalsForPredicate(
            predicate, variable, builtins, mathematics, angles, assumptions);
        if (!constraint)
            return std::nullopt;
        auto intersection = intersectIntervalSets(
            domain, *constraint, builtins, mathematics, assumptions);
        if (!intersection)
            return std::nullopt;
        domain = std::move(*intersection);
        if (domain.empty())
            break;
    }
    return domain;
}

[[nodiscard]] mathematics::AssumptionSet interiorAssumptions(
    const RealDomainInterval& interval,
    const mathematics::AssumptionSet& assumptions,
    const expression::Symbol& variable) {
    mathematics::AssumptionSet result = withRealVariable(assumptions, variable);
    if (interval.lower)
        result.add(mathematics::relation(
            RelationKind::Greater, Expr{variable}, *interval.lower));
    if (interval.upper)
        result.add(mathematics::relation(
            RelationKind::Less, Expr{variable}, *interval.upper));
    return result;
}

[[nodiscard]] bool pointInsideInterval(
    const Expr& point,
    const RealDomainInterval& interval,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions,
    bool strict) {
    if (interval.lower) {
        const auto comparison = compareFinite(
            point, *interval.lower, builtins, mathematics, assumptions);
        if (!comparison || *comparison < 0
            || (*comparison == 0 && (strict || !interval.lowerInclusive)))
            return false;
    }
    if (interval.upper) {
        const auto comparison = compareFinite(
            point, *interval.upper, builtins, mathematics, assumptions);
        if (!comparison || *comparison > 0
            || (*comparison == 0 && (strict || !interval.upperInclusive)))
            return false;
    }
    return true;
}

[[nodiscard]] std::optional<SolutionSet> solveCriticalEquation(
    const Expr& derivative,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    std::optional<Expr> normalizedRelation;
    if (derivative.isCall()) {
        const auto* definition = builtins.find(derivative.asCall().head);
        if (definition && definition->id == BuiltinId::Subtract
            && derivative.asCall().arguments.size() == 2) {
            normalizedRelation = relationExpr(
                RelationKind::Equal,
                derivative.asCall().arguments[0], derivative.asCall().arguments[1], builtins);
        }
    }
    Expr relation = normalizeForSolve(
        normalizedRelation.value_or(
            relationExpr(RelationKind::Equal, derivative, integerExpr(0), builtins)),
        builtins, mathematics, angles, assumptions);

    if (auto radical = solveRadicalRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return radical;
    if (auto exponential = solveRealExponentialRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return exponential;
    if (auto periodic = solveRealPeriodicFunctionRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return periodic;
    if (auto injective = solveRealInjectiveFunctionRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return injective;

    SolveConstraints realConstraint;
    realConstraint.domain = NumericDomain::Real;
    SolutionSet polynomial = applySolveConstraints(
        solveUnivariatePolynomialRelation(
            relation, variable, builtins, mathematics, angles),
        realConstraint, builtins, mathematics, angles);
    if (polynomial.kind() != SolutionSetKind::Unresolved)
        return polynomial;
    if (auto algebraic = solveRealAlgebraicPolynomialEquation(
            relation, variable, builtins, mathematics, angles))
        return algebraic;
    return std::nullopt;
}

[[nodiscard]] std::optional<std::vector<Expr>> criticalPoints(
    const Expr& derivative,
    const std::vector<RealDomainInterval>& domain,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto solutions = solveCriticalEquation(
        derivative, variable, builtins, mathematics, angles, assumptions);
    if (!solutions)
        return std::vector<Expr>{};
    if (solutions->kind() == SolutionSetKind::Empty)
        return std::vector<Expr>{};
    if (solutions->kind() != SolutionSetKind::Finite
        || !solutions->conditions().empty()
        || solutions->branches().size() > limits::realCriticalPoints)
        return std::nullopt;

    std::vector<Expr> points;
    for (const SolutionBranch& branch : solutions->branches()) {
        if (!branch.unconditional() || !branch.freeVariables.empty()
            || branch.bindings.size() != 1
            || !branch.bindings.front().variable.sameIdentity(variable))
            return std::nullopt;
        const Expr& point = branch.bindings.front().value;
        bool inside = false;
        for (const RealDomainInterval& interval : domain) {
            if (pointInsideInterval(
                    point, interval, builtins, mathematics, assumptions, true)) {
                inside = true;
                break;
            }
        }
        if (inside)
            points.push_back(point);
    }
    return points;
}

[[nodiscard]] std::optional<std::vector<RealDomainInterval>> splitAtCriticalPoints(
    const std::vector<RealDomainInterval>& domain,
    std::vector<Expr> points,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    std::vector<RealDomainInterval> result;
    for (const RealDomainInterval& interval : domain) {
        std::vector<Expr> local;
        for (const Expr& point : points)
            if (pointInsideInterval(
                    point, interval, builtins, mathematics, assumptions, true))
                local.push_back(point);

        std::stable_sort(local.begin(), local.end(), [&](const Expr& lhs, const Expr& rhs) {
            const auto comparison = compareFinite(
                lhs, rhs, builtins, mathematics, assumptions);
            return comparison && *comparison < 0;
        });
        for (std::size_t i = 1; i < local.size(); ++i) {
            const auto comparison = compareFinite(
                local[i - 1], local[i], builtins, mathematics, assumptions);
            if (!comparison)
                return std::nullopt;
        }

        if (local.empty()) {
            result.push_back(interval);
            continue;
        }

        std::optional<Expr> lower = interval.lower;
        bool lowerInclusive = interval.lowerInclusive;
        for (const Expr& point : local) {
            result.push_back(RealDomainInterval{
                lower, lowerInclusive, point, true});
            lower = point;
            lowerInclusive = true;
        }
        result.push_back(RealDomainInterval{
            lower, lowerInclusive, interval.upper, interval.upperInclusive});
    }
    if (result.size() > limits::realDomainPieces)
        return std::nullopt;
    return result;
}

[[nodiscard]] bool openIntervalContainedIn(
    const RealDomainInterval& interval,
    const RealDomainInterval& container,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    if (!interval.lower) {
        if (container.lower)
            return false;
    }
    else if (container.lower) {
        const auto comparison = compareFinite(
            *container.lower, *interval.lower, builtins, mathematics, assumptions);
        if (!comparison || *comparison > 0)
            return false;
    }

    if (!interval.upper) {
        if (container.upper)
            return false;
    }
    else if (container.upper) {
        const auto comparison = compareFinite(
            *container.upper, *interval.upper, builtins, mathematics, assumptions);
        if (!comparison || *comparison < 0)
            return false;
    }
    return true;
}

[[nodiscard]] bool derivativeHasSignOnInterior(
    const Expr& derivative,
    RelationKind signRelation,
    const RealDomainInterval& interval,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::AssumptionSet local = interiorAssumptions(
        interval, assumptions, variable);
    const mathematics::Predicate predicate = mathematics::relation(
        signRelation, derivative, integerExpr(0));
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, local};
    if (knowledge.prove(predicate) == TruthValue::True)
        return true;

    const auto signSet = intervalsForPredicate(
        predicate, variable, builtins, mathematics, angles, assumptions);
    if (!signSet)
        return false;
    for (const RealDomainInterval& candidate : *signSet)
        if (openIntervalContainedIn(
                interval, candidate, builtins, mathematics, assumptions))
            return true;
    return false;
}

[[nodiscard]] bool zeroSetHasEmptyInterior(
    const Expr& derivative,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto zeros = solveCriticalEquation(
        derivative, variable, builtins, mathematics, angles, assumptions);
    if (!zeros)
        return false;
    if (zeros->kind() == SolutionSetKind::Empty)
        return true;
    if (zeros->kind() != SolutionSetKind::Finite
        || zeros->branches().size() > limits::realCriticalPoints)
        return false;

    // finite個のbranchで，各自由parameterがIntegerならzero setは高々可算であり，
    // 実区間を含み得ない。これによりf'>=0かつ零点が孤立/周期離散族なら，
    // 非減少ではなくstrict increasingまで証明できる。
    for (const SolutionBranch& branch : zeros->branches()) {
        if (branch.bindings.size() != 1
            || !branch.bindings.front().variable.sameIdentity(variable))
            return false;
        for (const SolverVariable& parameter : branch.freeVariables)
            if (parameter.domain != NumericDomain::Integer)
                return false;
    }
    return true;
}

[[nodiscard]] RealIntervalMonotonicity monotonicityOn(
    const Expr& derivative,
    const RealDomainInterval& interval,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::AssumptionSet local = interiorAssumptions(
        interval, assumptions, variable);
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, local};
    const auto facts = knowledge.facts(derivative);
    if (facts.isProvablyReal()) {
        if (facts.sign == RealSign::Positive)
            return RealIntervalMonotonicity::Increasing;
        if (facts.sign == RealSign::Negative)
            return RealIntervalMonotonicity::Decreasing;
        if (facts.sign == RealSign::Zero)
            return RealIntervalMonotonicity::Constant;
    }

    // ValueFactsだけでsignが閉じないrational derivativeは，既存のexact sign-chartを
    // 再利用する。臨界点で分割済みなので，各区間内部全体が
    // derivative>0/<0の解区間へ含まれることを証明できればstrict monotoneである。
    if (derivativeHasSignOnInterior(
            derivative, RelationKind::Greater, interval, variable,
            builtins, mathematics, angles, assumptions))
        return RealIntervalMonotonicity::Increasing;
    if (derivativeHasSignOnInterior(
            derivative, RelationKind::Less, interval, variable,
            builtins, mathematics, angles, assumptions))
        return RealIntervalMonotonicity::Decreasing;

    // f'>=0 / <=0だけではconstant区間を除外できない。zero setの完全解が
    // finiteまたはInteger parameter族で高々可算と証明できた場合だけstrict化する。
    const bool thinZeroSet = zeroSetHasEmptyInterior(
        derivative, variable, builtins, mathematics, angles, assumptions);
    if (thinZeroSet && derivativeHasSignOnInterior(
            derivative, RelationKind::GreaterEqual, interval, variable,
            builtins, mathematics, angles, assumptions))
        return RealIntervalMonotonicity::Increasing;
    if (thinZeroSet && derivativeHasSignOnInterior(
            derivative, RelationKind::LessEqual, interval, variable,
            builtins, mathematics, angles, assumptions))
        return RealIntervalMonotonicity::Decreasing;
    return RealIntervalMonotonicity::Unknown;
}

[[nodiscard]] std::optional<Expr> endpointLimit(
    const Expr& expression,
    const RealDomainInterval& interval,
    bool lower,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    Expr point = lower
        ? (interval.lower ? *interval.lower : negativeInfinity(builtins, infinity))
        : (interval.upper ? *interval.upper : Expr{infinity});
    symbolic::LimitDirection direction = symbolic::LimitDirection::TwoSided;
    if (lower && interval.lower)
        direction = symbolic::LimitDirection::Right;
    else if (!lower && interval.upper)
        direction = symbolic::LimitDirection::Left;

    Expr result = symbolic::limitExpression(
        expression, variable, point, direction,
        builtins, mathematics, angles, infinity, assumptions);
    if (unresolvedLimit(result, builtins))
        return std::nullopt;
    return result;
}

[[nodiscard]] bool endpointIncludedInRange(
    const RealDomainInterval& interval,
    bool lower,
    const std::optional<Expr>& limit,
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const bool included = lower ? interval.lowerInclusive : interval.upperInclusive;
    const std::optional<Expr>& point = lower ? interval.lower : interval.upper;
    if (!included || !point || !limit)
        return false;
    Expr value = simplifyForSolve(
        symbolic::substituteSymbol(expression, variable, *point),
        builtins, mathematics, angles, assumptions);
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    return value == *limit
        || knowledge.prove(mathematics::relation(
            RelationKind::Equal, value, *limit)) == TruthValue::True;
}

[[nodiscard]] std::optional<RealValueRange> rangeFromMonotonicity(
    const Expr& expression,
    const RealDomainInterval& interval,
    RealIntervalMonotonicity monotonicity,
    const std::optional<Expr>& lowerLimit,
    const std::optional<Expr>& upperLimit,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!lowerLimit || !upperLimit
        || monotonicity == RealIntervalMonotonicity::Unknown)
        return std::nullopt;

    const bool lowerAttained = endpointIncludedInRange(
        interval, true, lowerLimit, expression, variable,
        builtins, mathematics, angles, assumptions);
    const bool upperAttained = endpointIncludedInRange(
        interval, false, upperLimit, expression, variable,
        builtins, mathematics, angles, assumptions);

    if (monotonicity == RealIntervalMonotonicity::Increasing)
        return RealValueRange{*lowerLimit, lowerAttained, *upperLimit, upperAttained};
    if (monotonicity == RealIntervalMonotonicity::Decreasing)
        return RealValueRange{*upperLimit, upperAttained, *lowerLimit, lowerAttained};
    if (monotonicity == RealIntervalMonotonicity::Constant) {
        const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
        if (*lowerLimit == *upperLimit
            || knowledge.prove(mathematics::relation(
                RelationKind::Equal, *lowerLimit, *upperLimit)) == TruthValue::True)
            return RealValueRange{*lowerLimit, true, *upperLimit, true};
    }
    return std::nullopt;
}

} // namespace

RealDomainAnalysis analyzeRealDomain(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    RealDomainAnalysis result;
    if (expressionNodeCount(expression, limits::realAnalysisNodes) > limits::realAnalysisNodes)
        return result;

    const mathematics::AssumptionSet realAssumptions = withRealVariable(
        assumptions, variable);
    auto domain = realDomainIntervals(
        expression, variable, builtins, mathematics, angles, realAssumptions);
    if (!domain)
        return result;
    result.complete = true;
    result.intervals = std::move(*domain);
    return result;
}

RealFunctionAnalysis analyzeRealFunction(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions) {
    RealFunctionAnalysis result;
    const RealDomainAnalysis domainAnalysis = analyzeRealDomain(
        expression, variable, builtins, mathematics, angles, assumptions);
    if (!domainAnalysis.complete)
        return result;

    const mathematics::AssumptionSet realAssumptions = withRealVariable(
        assumptions, variable);
    result.domainComplete = true;
    result.domain = domainAnalysis.intervals;
    if (result.domain.empty())
        return result;
    const auto& domain = result.domain;

    Expr derivative = symbolic::differentiateExpression(
        expression, variable, builtins, mathematics, angles);
    if (containsBuiltinCall(derivative, BuiltinId::Derivative, builtins)
        || expressionNodeCount(derivative, limits::realAnalysisNodes) > limits::realAnalysisNodes) {
        // domainだけは完全に分かっているため，pieceをmonotonicity unknownで返す。
        for (const RealDomainInterval& interval : domain)
            result.pieces.push_back(RealIntervalFunctionAnalysis{
                interval, RealIntervalMonotonicity::Unknown,
                std::nullopt, std::nullopt, std::nullopt});
        return result;
    }
    derivative = simplifyForSolve(
        std::move(derivative), builtins, mathematics, angles, realAssumptions);
    result.derivative = derivative;

    auto points = criticalPoints(
        derivative, domain, variable,
        builtins, mathematics, angles, realAssumptions);
    std::vector<RealDomainInterval> pieces = domain;
    if (points) {
        if (auto split = splitAtCriticalPoints(
                domain, std::move(*points), builtins, mathematics, realAssumptions))
            pieces = std::move(*split);
    }

    result.pieces.reserve(pieces.size());
    for (const RealDomainInterval& interval : pieces) {
        RealIntervalFunctionAnalysis analysis;
        analysis.domain = interval;
        analysis.monotonicity = monotonicityOn(
            derivative, interval, variable, builtins, mathematics,
            angles, realAssumptions);
        analysis.lowerLimit = endpointLimit(
            expression, interval, true, variable,
            builtins, mathematics, angles, infinitySymbol, realAssumptions);
        analysis.upperLimit = endpointLimit(
            expression, interval, false, variable,
            builtins, mathematics, angles, infinitySymbol, realAssumptions);
        analysis.range = rangeFromMonotonicity(
            expression, interval, analysis.monotonicity,
            analysis.lowerLimit, analysis.upperLimit,
            variable, builtins, mathematics, angles,
            realAssumptions);
        result.pieces.push_back(std::move(analysis));
    }
    return result;
}

} // namespace mmcal::solver
