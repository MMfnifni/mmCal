// 極限limit
#include "limit.hpp"

#include "mathematics/definedness.hpp"
#include "mathematics/knowledge_context.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "simplification/full_simplifier.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::symbolic {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

constexpr std::size_t maximumLimitDepth = 24;
constexpr std::size_t maximumLHopitalSteps = 12;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr rational(const Rational& value) {
    return Expr{Number{value}};
}

[[nodiscard]] bool isZero(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return builtins.isCallTo(expression, id);
}

[[nodiscard]] Expr call(
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id,
    std::vector<Expr> arguments) {
    return Expr::call(builtins.symbol(id), std::move(arguments));
}

[[nodiscard]] Expr simplify(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
}

[[nodiscard]] Expr fullSimplify(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplification::fullSimplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions},
        simplification::FullSimplificationOptions{48});
}

[[nodiscard]] Expr negate(
    Expr value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplify(call(builtins, BuiltinId::Negate, {std::move(value)}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr divide(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplify(call(builtins, BuiltinId::Divide, {std::move(lhs), std::move(rhs)}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr unresolved(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& point,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Expr> arguments{expression, Expr{variable}, point};
    if (direction == LimitDirection::Left)
        arguments.push_back(integer(-1));
    else if (direction == LimitDirection::Right)
        arguments.push_back(integer(1));
    return call(builtins, BuiltinId::Limit, std::move(arguments));
}

[[nodiscard]] bool isInfinity(const Expr& expression, const expression::Symbol& infinity) {
    return expression.isSymbol() && expression.asSymbol().sameIdentity(infinity);
}

[[nodiscard]] bool isNegativeInfinity(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    return isHead(expression, builtins, BuiltinId::Negate)
        && expression.asCall().arguments.size() == 1
        && isInfinity(expression.asCall().arguments[0], infinity);
}

[[nodiscard]] Expr positiveInfinity(const expression::Symbol& infinity) {
    return Expr{infinity};
}

[[nodiscard]] Expr signedInfinity(
    int sign,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    if (sign >= 0)
        return positiveInfinity(infinity);
    return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<Rational> exactRealRational(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    return expression.asNumber().asReal().toRational();
}

[[nodiscard]] bool domainConditionsHold(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const auto conditions = mathematics::expressionDomainConditions(expression, builtins, mathematics);
    if (!conditions)
        return true;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    for (const auto& predicate : conditions->predicates()) {
        if (knowledge.prove(predicate) == mathematics::TruthValue::False)
            return false;
    }
    return true;
}

struct LocalPolynomialBehavior final {
    std::size_t zeroOrder = 0;
    Rational leading{};
};

[[nodiscard]] std::optional<LocalPolynomialBehavior> localBehavior(
    RationalPolynomial polynomial,
    const Rational& point) {
    if (polynomial.isZero())
        return std::nullopt;

    std::size_t order = 0;
    while (evaluatePolynomial(polynomial, point).isZero()) {
        const auto quotient = divideByLinearFactor(polynomial, point);
        if (!quotient)
            return std::nullopt;
        polynomial = *quotient;
        ++order;
    }
    return LocalPolynomialBehavior{order, evaluatePolynomial(polynomial, point)};
}

[[nodiscard]] int rationalSign(const Rational& value) {
    if (value.isZero()) return 0;
    return value.numerator().isNegative() ? -1 : 1;
}

[[nodiscard]] std::optional<Expr> rationalFunctionFiniteLimit(
    const Expr& expression,
    const expression::Symbol& variable,
    const Rational& point,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    Expr numerator = expression;
    Expr denominator = integer(1);
    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        numerator = expression.asCall().arguments[0];
        denominator = expression.asCall().arguments[1];
    }

    const auto np = toRationalPolynomial(numerator, variable, builtins, {256, 2048});
    const auto dp = toRationalPolynomial(denominator, variable, builtins, {256, 2048});
    if (!np || !dp || dp->isZero())
        return std::nullopt;

    const auto nb = localBehavior(*np, point);
    const auto db = localBehavior(*dp, point);
    if (!nb || !db)
        return std::nullopt;

    if (nb->zeroOrder > db->zeroOrder)
        return integer(0);

    if (nb->zeroOrder == db->zeroOrder)
        return rational(nb->leading / db->leading);

    const std::size_t poleOrder = db->zeroOrder - nb->zeroOrder;
    int sign = rationalSign(nb->leading / db->leading);
    if (direction == LimitDirection::TwoSided && (poleOrder % 2) != 0)
        return std::nullopt;
    if (direction == LimitDirection::Left && (poleOrder % 2) != 0)
        sign = -sign;
    return signedInfinity(sign, builtins, mathematics, angles, infinity, assumptions);
}

[[nodiscard]] std::optional<Expr> rationalFunctionInfiniteLimit(
    const Expr& expression,
    const expression::Symbol& variable,
    bool negativeInfinityPoint,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    Expr numerator = expression;
    Expr denominator = integer(1);
    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        numerator = expression.asCall().arguments[0];
        denominator = expression.asCall().arguments[1];
    }

    const auto np = toRationalPolynomial(numerator, variable, builtins, {256, 2048});
    const auto dp = toRationalPolynomial(denominator, variable, builtins, {256, 2048});
    if (!np || !dp || np->isZero() || dp->isZero())
        return std::nullopt;

    if (np->degree() < dp->degree())
        return integer(0);

    const Rational ratio = np->coefficient(np->degree()) / dp->coefficient(dp->degree());
    if (np->degree() == dp->degree())
        return rational(ratio);

    const std::size_t degreeDifference = np->degree() - dp->degree();
    int sign = rationalSign(ratio);
    if (negativeInfinityPoint && (degreeDifference % 2) != 0)
        sign = -sign;
    return signedInfinity(sign, builtins, mathematics, angles, infinity, assumptions);
}

[[nodiscard]] Expr inverseHalfTurn(
    bool negative,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr value = [&]() -> Expr {
        switch (angles.defaultUnit()) {
        case mathematics::AngleUnit::Degree:
            return integer(90);
        case mathematics::AngleUnit::Gradian:
            return integer(100);
        case mathematics::AngleUnit::Radian: {
            const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
            return divide(Expr{pi->symbol}, integer(2),
                builtins, mathematics, angles, assumptions);
        }
        }
        return integer(0);
    }();
    return negative
        ? negate(std::move(value), builtins, mathematics, angles, assumptions)
        : value;
}

[[nodiscard]] std::optional<Rational> affineSlopeAtInfinity(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    const auto polynomial = toRationalPolynomial(expression, variable, builtins, {1, 4});
    if (!polynomial || polynomial->degree() != 1)
        return std::nullopt;
    return polynomial->coefficient(1);
}

[[nodiscard]] Expr limitCore(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& point,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions,
    std::size_t depth) {
    if (depth > maximumLimitDepth)
        return unresolved(expression, variable, point, direction, builtins);
    if (!containsSymbol(expression, variable))
        return expression;

    const bool atPositiveInfinity = isInfinity(point, infinity);
    const bool atNegativeInfinity = isNegativeInfinity(point, builtins, infinity);

    // 線形な外側構造は先に分解する。improper integralの原始函数で
    // -exp[-x] や atan[x]-atan[0] を無用に未評価へ落とさない。
    if (expression.isCall()) {
        const auto* outer = builtins.find(expression.asCall().head);
        const auto& arguments = expression.asCall().arguments;
        if (outer && outer->id == BuiltinId::Negate && arguments.size() == 1) {
            Expr inner = limitCore(arguments[0], variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, depth + 1);
            if (!isHead(inner, builtins, BuiltinId::Limit))
                return negate(std::move(inner), builtins, mathematics, angles, assumptions);
        }
        if (outer && (outer->id == BuiltinId::Add || outer->id == BuiltinId::Subtract)) {
            std::vector<Expr> values;
            values.reserve(arguments.size());
            bool positiveInfinitySeen = false;
            bool negativeInfinitySeen = false;
            for (const Expr& argument : arguments) {
                Expr value = limitCore(argument, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, depth + 1);
                if (isHead(value, builtins, BuiltinId::Limit))
                    return unresolved(expression, variable, point, direction, builtins);
                if (outer->id == BuiltinId::Subtract && values.size() == 1)
                    value = negate(std::move(value), builtins, mathematics, angles, assumptions);
                positiveInfinitySeen |= isInfinity(value, infinity);
                negativeInfinitySeen |= isNegativeInfinity(value, builtins, infinity);
                values.push_back(std::move(value));
            }
            if (positiveInfinitySeen && negativeInfinitySeen)
                return unresolved(expression, variable, point, direction, builtins);
            if (positiveInfinitySeen)
                return Expr{infinity};
            if (negativeInfinitySeen)
                return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);
            const BuiltinId combinedId = outer->id == BuiltinId::Subtract
                ? BuiltinId::Add : outer->id;
            return simplify(call(builtins, combinedId, std::move(values)),
                builtins, mathematics, angles, assumptions);
        }
    }

    if (atPositiveInfinity || atNegativeInfinity) {
        if (auto rationalLimit = rationalFunctionInfiniteLimit(
                expression, variable, atNegativeInfinity,
                builtins, mathematics, angles, infinity, assumptions))
            return *rationalLimit;

        if (expression.isCall() && expression.asCall().arguments.size() == 1) {
            const auto* definition = builtins.find(expression.asCall().head);
            const Expr& argument = expression.asCall().arguments[0];
            if (definition) {
                const auto slope = affineSlopeAtInfinity(argument, variable, builtins);
                const int infinitySign = atNegativeInfinity ? -1 : 1;
                if (slope && !slope->isZero()) {
                    const int argumentSign = rationalSign(*slope) * infinitySign;
                    switch (definition->id) {
                    case BuiltinId::Exp:
                        return argumentSign > 0 ? Expr{infinity} : integer(0);
                    case BuiltinId::Log:
                        if (argumentSign > 0)
                            return Expr{infinity};
                        break;
                    case BuiltinId::Atan:
                        return inverseHalfTurn(argumentSign < 0,
                            builtins, mathematics, angles, assumptions);
                    case BuiltinId::Tanh:
                    case BuiltinId::Erf:
                        return integer(argumentSign > 0 ? 1 : -1);
                    case BuiltinId::Erfc:
                        return integer(argumentSign > 0 ? 0 : 2);
                    case BuiltinId::Asinh:
                        return signedInfinity(argumentSign,
                            builtins, mathematics, angles, infinity, assumptions);
                    default:
                        break;
                    }
                }
            }
        }
        return unresolved(expression, variable, point, direction, builtins);
    }

    const auto rationalPoint = exactRealRational(point);
    if (rationalPoint) {
        if (auto rationalLimit = rationalFunctionFiniteLimit(
                expression, variable, *rationalPoint, direction,
                builtins, mathematics, angles, infinity, assumptions))
            return *rationalLimit;

        // principal Logの実軸右極限。左側は -Infinity + I Pi となるため、
        // 現在のInfinity sentinelで複素無限大を捏造せず未評価に残す。
        if (rationalPoint->isZero() && direction == LimitDirection::Right
            && isHead(expression, builtins, BuiltinId::Log)
            && expression.asCall().arguments.size() == 1
            && expression.asCall().arguments[0].isSymbol()
            && expression.asCall().arguments[0].asSymbol().sameIdentity(variable))
            return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);

        // x^a log[x] -> 0 (x->0+, a>0)。improper integralのendpointで頻出し、
        // 実軸上の標準極限としてbranch安全に扱える。
        if (rationalPoint->isZero() && direction == LimitDirection::Right
            && isHead(expression, builtins, BuiltinId::Multiply)) {
            const auto& factors = expression.asCall().arguments;
            if (factors.size() == 2) {
                for (std::size_t logIndex = 0; logIndex < 2; ++logIndex) {
                    const Expr& logarithm = factors[logIndex];
                    const Expr& vanishing = factors[1 - logIndex];
                    if (!isHead(logarithm, builtins, BuiltinId::Log)
                        || logarithm.asCall().arguments.size() != 1
                        || !logarithm.asCall().arguments[0].isSymbol()
                        || !logarithm.asCall().arguments[0].asSymbol().sameIdentity(variable))
                        continue;
                    Rational exponent{BigInt{0}};
                    bool matched = false;
                    if (vanishing.isSymbol() && vanishing.asSymbol().sameIdentity(variable)) {
                        exponent = Rational{BigInt{1}};
                        matched = true;
                    }
                    else if (isHead(vanishing, builtins, BuiltinId::Power)
                        && vanishing.asCall().arguments.size() == 2
                        && vanishing.asCall().arguments[0].isSymbol()
                        && vanishing.asCall().arguments[0].asSymbol().sameIdentity(variable)) {
                        if (const auto value = exactRealRational(vanishing.asCall().arguments[1])) {
                            exponent = *value;
                            matched = !exponent.isZero() && !exponent.numerator().isNegative();
                        }
                    }
                    if (matched)
                        return integer(0);
                }
            }
        }
    }

    Expr rawSubstituted = substituteSymbol(expression, variable, point);
    if (domainConditionsHold(rawSubstituted, builtins, mathematics, assumptions)) {
        Expr substituted = simplify(
            std::move(rawSubstituted), builtins, mathematics, angles, assumptions);
        if (!containsSymbol(substituted, variable))
            return fullSimplify(std::move(substituted), builtins, mathematics, angles, assumptions);
    }

    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        Expr numerator = expression.asCall().arguments[0];
        Expr denominator = expression.asCall().arguments[1];
        for (std::size_t step = 0; step < maximumLHopitalSteps; ++step) {
            Expr numeratorAt = simplify(substituteSymbol(numerator, variable, point),
                builtins, mathematics, angles, assumptions);
            Expr denominatorAt = simplify(substituteSymbol(denominator, variable, point),
                builtins, mathematics, angles, assumptions);
            if (!(isZero(numeratorAt) && isZero(denominatorAt)))
                break;
            numerator = differentiateExpression(
                numerator, variable, builtins, mathematics, angles);
            denominator = differentiateExpression(
                denominator, variable, builtins, mathematics, angles);
            if (isHead(numerator, builtins, BuiltinId::Derivative)
                || isHead(denominator, builtins, BuiltinId::Derivative))
                break;
            Expr quotient = divide(numerator, denominator,
                builtins, mathematics, angles, assumptions);
            Expr result = limitCore(quotient, variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, depth + 1);
            if (!isHead(result, builtins, BuiltinId::Limit))
                return result;
        }
    }

    if (expression.isCall()) {
        const auto* definition = builtins.find(expression.asCall().head);
        const auto& arguments = expression.asCall().arguments;
        if (definition && definition->id == BuiltinId::Multiply && arguments.size() == 2) {
            // x*log(x) -> 0 (x->0+) など、0*Infinity型の代表的な形は商へ直してL'Hopitalへ送る。
            for (std::size_t logIndex = 0; logIndex < 2; ++logIndex) {
                const Expr& logarithm = arguments[logIndex];
                const Expr& factor = arguments[1 - logIndex];
                if (!isHead(logarithm, builtins, BuiltinId::Log)
                    || logarithm.asCall().arguments.size() != 1)
                    continue;
                Expr rewritten = divide(logarithm,
                    divide(integer(1), factor, builtins, mathematics, angles, assumptions),
                    builtins, mathematics, angles, assumptions);
                Expr result = limitCore(rewritten, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, depth + 1);
                if (!isHead(result, builtins, BuiltinId::Limit))
                    return result;
            }
        }
    }

    return unresolved(expression, variable, point, direction, builtins);
}

} // namespace

Expr limitExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& point,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions) {
    mathematics::AssumptionSet local = assumptions;
    if (!isInfinity(point, infinitySymbol)
        && !isNegativeInfinity(point, builtins, infinitySymbol)
        && direction != LimitDirection::TwoSided) {
        local.add(mathematics::elementOf(Expr{variable}, mathematics::NumericDomain::Real));
        local.add(mathematics::relation(
            direction == LimitDirection::Right
                ? mathematics::RelationKind::Greater
                : mathematics::RelationKind::Less,
            Expr{variable}, point));
    }
    const Expr prepared = simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles, local});
    return limitCore(prepared, variable, point, direction,
        builtins, mathematics, angles, infinitySymbol, local, 0);
}

} // namespace mmcal::symbolic
