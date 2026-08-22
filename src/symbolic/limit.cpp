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

[[nodiscard]] Expr caseBranch(
    const evaluation::BuiltinRegistry& builtins,
    Expr value,
    std::optional<Expr> condition = std::nullopt) {
    std::vector<Expr> arguments{std::move(value)};
    if (condition)
        arguments.push_back(std::move(*condition));
    return call(builtins, BuiltinId::CaseBranch, std::move(arguments));
}

[[nodiscard]] Expr cases(
    const evaluation::BuiltinRegistry& builtins,
    std::vector<Expr> branches) {
    return call(builtins, BuiltinId::Cases, std::move(branches));
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

[[nodiscard]] bool isUnaryFunctionOfVariable(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return isHead(expression, builtins, id)
        && expression.asCall().arguments.size() == 1
        && expression.asCall().arguments[0].isSymbol()
        && expression.asCall().arguments[0].asSymbol().sameIdentity(variable);
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

[[nodiscard]] Expr imaginaryPi(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
    Expr imaginaryUnit{Number::complex(
        numeric::RealNumber{BigInt{0}}, numeric::RealNumber{BigInt{1}})};
    return simplify(call(builtins, BuiltinId::Multiply, {
        std::move(imaginaryUnit), Expr{pi->symbol}}),
        builtins, mathematics, angles, assumptions);
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


[[nodiscard]] bool isSignedInfinity(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    return isInfinity(expression, infinity)
        || isNegativeInfinity(expression, builtins, infinity);
}

[[nodiscard]] bool isOscillatoryPeriodicBuiltin(BuiltinId id) {
    return id == BuiltinId::Sin || id == BuiltinId::Cos || id == BuiltinId::Tan;
}

[[nodiscard]] bool isExactRealRationalFunction(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        const auto numerator = toRationalPolynomial(
            expression.asCall().arguments[0], variable, builtins, {256, 2048});
        const auto denominator = toRationalPolynomial(
            expression.asCall().arguments[1], variable, builtins, {256, 2048});
        return numerator && denominator && !denominator->isZero();
    }
    return toRationalPolynomial(expression, variable, builtins, {256, 2048}).has_value();
}

[[nodiscard]] bool isBoundedRealTrigFactor(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 1)
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || (definition->id != BuiltinId::Sin && definition->id != BuiltinId::Cos))
        return false;
    return isExactRealRationalFunction(expression.asCall().arguments[0], variable, builtins);
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
    const expression::Symbol* complexInfinity,
    const expression::Symbol* indeterminate,
    std::size_t depth) {
    if (depth > maximumLimitDepth)
        return unresolved(expression, variable, point, direction, builtins);
    if (!containsSymbol(expression, variable))
        return expression;

    // casesの条件が極限変数に依存しない場合だけ，各branchの極限へ分配できる。
    // 変数依存条件ではapproach directionとbranch境界の解析が必要であり，
    // 点代入でbranch値を先に評価すると0/0等を誤ってDomainErrorへ落とすため未解決で保持する。
    if (isHead(expression, builtins, BuiltinId::Cases)) {
        const auto& sourceBranches = expression.asCall().arguments;
        std::vector<Expr> branches;
        branches.reserve(sourceBranches.size());
        for (const Expr& branchExpression : sourceBranches) {
            if (!isHead(branchExpression, builtins, BuiltinId::CaseBranch)
                || branchExpression.asCall().arguments.empty()
                || branchExpression.asCall().arguments.size() > 2)
                return unresolved(expression, variable, point, direction, builtins);
            const auto& branch = branchExpression.asCall().arguments;
            if (branch.size() == 2 && containsSymbol(branch[1], variable))
                return unresolved(expression, variable, point, direction, builtins);
            Expr value = limitCore(
                branch[0], variable, point, direction, builtins, mathematics, angles,
                infinity, assumptions, complexInfinity, indeterminate, depth + 1);
            if (branch.size() == 2)
                branches.push_back(caseBranch(builtins, std::move(value), branch[1]));
            else
                branches.push_back(caseBranch(builtins, std::move(value)));
        }
        return cases(builtins, std::move(branches));
    }

    // 実引数が無限大へ走る周期三角函数は単一の極限値を持たない。
    // これは「未実装」ではなく不存在を証明できる場合なので，Evaluatorから
    // Indeterminate atomが渡されていればそれを返す。有限点の二側極限では，
    // どちらか一方で無限振動を証明できれば二側極限の不存在も確定する。
    if (expression.isCall() && expression.asCall().arguments.size() == 1) {
        const auto* definition = builtins.find(expression.asCall().head);
        if (definition && isOscillatoryPeriodicBuiltin(definition->id) && indeterminate) {
            const Expr& argument = expression.asCall().arguments[0];
            if (direction == LimitDirection::TwoSided
                && !isInfinity(point, infinity)
                && !isNegativeInfinity(point, builtins, infinity)) {
                Expr left = limitCore(argument, variable, point, LimitDirection::Left,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity,
                    indeterminate, depth + 1);
                if (isSignedInfinity(left, builtins, infinity))
                    return Expr{*indeterminate};
                Expr right = limitCore(argument, variable, point, LimitDirection::Right,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity,
                    indeterminate, depth + 1);
                if (isSignedInfinity(right, builtins, infinity))
                    return Expr{*indeterminate};
            }
            else {
                Expr inner = limitCore(argument, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity,
                    indeterminate, depth + 1);
                if (isSignedInfinity(inner, builtins, infinity))
                    return Expr{*indeterminate};
            }
        }
    }

    const bool atPositiveInfinity = isInfinity(point, infinity);
    const bool atNegativeInfinity = isNegativeInfinity(point, builtins, infinity);

    // 線形な外側構造は先に分解する。improper integralの原始函数で
    // -exp[-x] や atan[x]-atan[0] を無用に未評価へ落とさない。
    if (expression.isCall()) {
        const auto* outer = builtins.find(expression.asCall().head);
        const auto& arguments = expression.asCall().arguments;
        if (outer && outer->id == BuiltinId::Negate && arguments.size() == 1) {
            Expr inner = limitCore(arguments[0], variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
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
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
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
        if (outer && outer->id == BuiltinId::Divide && arguments.size() == 2) {
            Expr numerator = limitCore(arguments[0], variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
            Expr denominator = limitCore(arguments[1], variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
            const bool numeratorFinite = !isHead(numerator, builtins, BuiltinId::Limit)
                && !isInfinity(numerator, infinity)
                && !isNegativeInfinity(numerator, builtins, infinity);
            const bool denominatorFinite = !isHead(denominator, builtins, BuiltinId::Limit)
                && !isInfinity(denominator, infinity)
                && !isNegativeInfinity(denominator, builtins, infinity);
            if (numeratorFinite && denominatorFinite && !isZero(denominator))
                return divide(std::move(numerator), std::move(denominator),
                    builtins, mathematics, angles, assumptions);
        }
        if (outer && outer->id == BuiltinId::Multiply) {
            // 実有理函数を引数に取るsin/cosは実軸上で絶対値1以下なので，
            // もう一方のfactorが0へ収束する2因子積はsqueeze theoremで0となる。
            if (arguments.size() == 2) {
                for (std::size_t boundedIndex = 0; boundedIndex < 2; ++boundedIndex) {
                    if (!isBoundedRealTrigFactor(arguments[boundedIndex], variable, builtins))
                        continue;
                    Expr vanishing = limitCore(arguments[1 - boundedIndex], variable, point, direction,
                        builtins, mathematics, angles, infinity, assumptions, complexInfinity,
                        indeterminate, depth + 1);
                    if (isZero(vanishing))
                        return integer(0);
                }
            }

            // すべてのfactorが有限極限へ収束する場合だけ積を合成する。
            // 0*Infinity等の不定形はここで決めず、既存の専用ruleへ残す。
            std::vector<Expr> values;
            values.reserve(arguments.size());
            bool finite = true;
            for (const Expr& argument : arguments) {
                Expr value = limitCore(argument, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
                if (isHead(value, builtins, BuiltinId::Limit)
                    || isInfinity(value, infinity)
                    || isNegativeInfinity(value, builtins, infinity)) {
                    finite = false;
                    break;
                }
                values.push_back(std::move(value));
            }
            if (finite)
                return simplify(call(builtins, BuiltinId::Multiply, std::move(values)),
                    builtins, mathematics, angles, assumptions);
        }
    }

    if (atPositiveInfinity || atNegativeInfinity) {
        if (auto rationalLimit = rationalFunctionInfiniteLimit(
                expression, variable, atNegativeInfinity,
                builtins, mathematics, angles, infinity, assumptions))
            return *rationalLimit;

        // Classical integral functions have branch-sensitive but exact real-axis asymptotics.
        // Ei(x) ~ exp(x)/x for x->+Infinity and tends to zero for x->-Infinity.
        if (isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::ExponentialIntegralEi))
            return atNegativeInfinity ? integer(0) : Expr{infinity};

        // Principal Ci(x) tends to zero on the positive real axis. On the negative real
        // branch cut Ci(-r)=Ci(r)+I Pi (r>0), hence x->-Infinity tends exactly to I Pi.
        if (isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::CosineIntegralCi))
            return atNegativeInfinity
                ? imaginaryPi(builtins, mathematics, angles, assumptions)
                : integer(0);

        // li(x)=Ei(Log(x)). On the positive real axis li(x)->+Infinity.
        // Along x->-Infinity on the principal branch, Log(x)=log|x|+I Pi and
        // |Ei(Log(x))| ~ |x|/|Log(x)| -> Infinity while the value is complex.
        // mmCal has no DirectedInfinity yet, so ComplexInfinity is the exact conservative
        // extended-complex result: it records unbounded magnitude without inventing a real value.
        if (isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::LogarithmicIntegralLi)) {
            if (!atNegativeInfinity)
                return Expr{infinity};
            if (complexInfinity)
                return Expr{*complexInfinity};
        }

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
        if (rationalPoint->isZero()) {
            // Ei(x)=gamma+log|x|+O(x) on the real axis, so both real one-sided
            // limits at its logarithmic singularity are -Infinity.
            if (isUnaryFunctionOfVariable(
                    expression, variable, builtins, BuiltinId::ExponentialIntegralEi))
                return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);

            // Ci(z)=gamma+Log(z)+O(z^2). From the right this is real -Infinity;
            // from the left the principal branch adds the bounded term I Pi. In either
            // case the directed limit is -Infinity, so the two-sided real limit agrees.
            if (isUnaryFunctionOfVariable(
                    expression, variable, builtins, BuiltinId::CosineIntegralCi))
                return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);

            // li(x)=Ei(Log(x)) tends to zero at the origin from either real side,
            // despite x=0 being a branch point of the complex principal function.
            if (isUnaryFunctionOfVariable(
                    expression, variable, builtins, BuiltinId::LogarithmicIntegralLi))
                return integer(0);
        }

        // li(x)=Ei(Log(x)) and Ei(t)->-Infinity as t->0 from either real side.
        // Hence both one-sided limits, and therefore the two-sided real limit, at x=1
        // are -Infinity even though li(1) itself is undefined.
        if (*rationalPoint == Rational{BigInt{1}}
            && isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::LogarithmicIntegralLi))
            return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);

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
                builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
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
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
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
    const mathematics::AssumptionSet& assumptions,
    const expression::Symbol* complexInfinitySymbol,
    const expression::Symbol* indeterminateSymbol) {
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
        builtins, mathematics, angles, infinitySymbol, local, complexInfinitySymbol,
        indeterminateSymbol, 0);
}

} // namespace mmcal::symbolic
