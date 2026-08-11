// 記号積分integrate
#include "integration.hpp"

#include "approximation/certified_evaluator.hpp"
#include "approximation/real_interval.hpp"
#include "mathematics/definedness.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "simplification/full_simplifier.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/limit.hpp"
#include "symbolic/algebra_transforms.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::symbolic {
namespace {

using approximation::CertifiedBinding;
using approximation::CertifiedEvaluator;
using approximation::CertifiedValue;
using approximation::RealInterval;
using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

constexpr std::size_t maximumIntegrationDepth = 24;
constexpr std::size_t maximumSubstitutionCandidates = 32;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr rational(const Rational& value) {
    return Expr{Number{value}};
}

[[nodiscard]] bool isZero(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] bool isOne(const Expr& expression) {
    return expression.isNumber()
        && expression.asNumber().isReal()
        && expression.asNumber().asReal().toRational() == Rational{BigInt{1}};
}

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return expression.isCall()
        && expression.asCall().head.sameIdentity(builtins.symbol(id));
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
    const mathematics::AngleSemantics& angles) {
    return simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles});
}

[[nodiscard]] Expr fullSimplify(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplification::fullSimplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles},
        simplification::FullSimplificationOptions{48});
}

[[nodiscard]] Expr add(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::vector<Expr> arguments) {
    return simplify(call(builtins, BuiltinId::Add, std::move(arguments)), builtins, mathematics, angles);
}

[[nodiscard]] Expr multiply(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::vector<Expr> arguments) {
    return simplify(call(builtins, BuiltinId::Multiply, std::move(arguments)), builtins, mathematics, angles);
}

[[nodiscard]] Expr divide(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    Expr numerator,
    Expr denominator) {
    if (denominator.isNumber() && denominator.asNumber().isReal()) {
        const Rational value = denominator.asNumber().asReal().toRational();
        if (value == Rational{BigInt{1}})
            return numerator;
        if (value == Rational{BigInt{-1}})
            return simplify(call(builtins, BuiltinId::Negate, {std::move(numerator)}),
                builtins, mathematics, angles);
    }
    return simplify(call(builtins, BuiltinId::Divide,
        {std::move(numerator), std::move(denominator)}), builtins, mathematics, angles);
}

[[nodiscard]] Expr subtract(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    Expr lhs,
    Expr rhs) {
    return simplify(call(builtins, BuiltinId::Subtract,
        {std::move(lhs), std::move(rhs)}), builtins, mathematics, angles);
}

[[nodiscard]] Expr negate(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    Expr value) {
    return simplify(call(builtins, BuiltinId::Negate, {std::move(value)}), builtins, mathematics, angles);
}

[[nodiscard]] Expr power(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    Expr base,
    Expr exponent) {
    return simplify(call(builtins, BuiltinId::Power,
        {std::move(base), std::move(exponent)}), builtins, mathematics, angles);
}

[[nodiscard]] Expr unresolved(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    return call(builtins, BuiltinId::SymbolicIntegral, {expression, Expr{variable}});
}

[[nodiscard]] bool containsVariable(
    const Expr& expression,
    const expression::Symbol& variable) {
    return containsSymbol(expression, variable);
}

[[nodiscard]] std::optional<Rational> exactRealRational(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    return expression.asNumber().asReal().toRational();
}

[[nodiscard]] bool provablyNonZero(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    if (expression.isNumber())
        return !expression.asNumber().isZero();
    const mathematics::ValueFacts facts = mathematics::inferValueFacts(
        expression, builtins, mathematics);
    return facts.sign == mathematics::RealSign::Positive
        || facts.sign == mathematics::RealSign::Negative
        || facts.provablyNonReal;
}

struct FactorSplit final {
    Expr constant;
    Expr dependent;
};

[[nodiscard]] FactorSplit splitConstantFactor(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!containsVariable(expression, variable))
        return FactorSplit{expression, integer(1)};

    std::vector<Expr> constants;
    std::vector<Expr> dependent;
    if (isHead(expression, builtins, BuiltinId::Multiply)) {
        for (const Expr& factor : expression.asCall().arguments) {
            if (containsVariable(factor, variable))
                dependent.push_back(factor);
            else
                constants.push_back(factor);
        }
    }
    else if (isHead(expression, builtins, BuiltinId::Negate)
        && expression.asCall().arguments.size() == 1) {
        FactorSplit inner = splitConstantFactor(
            expression.asCall().arguments[0], variable, builtins, mathematics, angles);
        inner.constant = negate(builtins, mathematics, angles, std::move(inner.constant));
        return inner;
    }
    else if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2
        && !containsVariable(expression.asCall().arguments[1], variable)) {
        FactorSplit inner = splitConstantFactor(
            expression.asCall().arguments[0], variable, builtins, mathematics, angles);
        inner.constant = divide(
            builtins, mathematics, angles,
            std::move(inner.constant), expression.asCall().arguments[1]);
        return inner;
    }
    else {
        dependent.push_back(expression);
    }

    return FactorSplit{
        constants.empty() ? integer(1)
            : multiply(builtins, mathematics, angles, std::move(constants)),
        dependent.empty() ? integer(1)
            : multiply(builtins, mathematics, angles, std::move(dependent))};
}

[[nodiscard]] std::optional<Expr> polynomialProportionalFactor(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto left = toExpressionPolynomial(
        lhs, variable, builtins, mathematics, angles, PolynomialConversionOptions{64, 256});
    auto right = toExpressionPolynomial(
        rhs, variable, builtins, mathematics, angles, PolynomialConversionOptions{64, 256});
    if (!left || !right || left->degree() != right->degree())
        return std::nullopt;

    std::optional<Expr> ratio;
    const std::size_t count = std::max(left->coefficients().size(), right->coefficients().size());
    for (std::size_t i = 0; i < count; ++i) {
        const Expr lc = left->coefficient(i);
        const Expr rc = right->coefficient(i);
        if (isZero(lc) && isZero(rc))
            continue;
        if (isZero(rc))
            return std::nullopt;
        if (!ratio) {
            ratio = divide(builtins, mathematics, angles, lc, rc);
            if (containsVariable(*ratio, variable))
                return std::nullopt;
            continue;
        }
        Expr difference = subtract(builtins, mathematics, angles,
            lc,
            multiply(builtins, mathematics, angles, {*ratio, rc}));
        if (!isZero(difference))
            return std::nullopt;
    }
    return ratio;
}

[[nodiscard]] std::optional<Expr> denominatorProportionalFactor(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs == rhs)
        return integer(1);
    if (const auto polynomialRatio = polynomialProportionalFactor(
            lhs, rhs, variable, builtins, mathematics, angles))
        return polynomialRatio;

    // principal sqrtについて sqrt(k z)=sqrt(k)sqrt(z) は k>0 の実数なら大域的に安全。
    // 二次式の平方完成で生じる正のRational scaleだけをここで認める。
    if (isHead(lhs, builtins, BuiltinId::Sqrt)
        && isHead(rhs, builtins, BuiltinId::Sqrt)
        && lhs.asCall().arguments.size() == 1
        && rhs.asCall().arguments.size() == 1) {
        const auto radicandRatio = polynomialProportionalFactor(
            lhs.asCall().arguments[0], rhs.asCall().arguments[0], variable,
            builtins, mathematics, angles);
        if (radicandRatio) {
            const auto exactRatio = exactRealRational(*radicandRatio);
            if (exactRatio && *exactRatio > Rational{BigInt{0}}) {
                return simplify(call(builtins, BuiltinId::Sqrt, {*radicandRatio}),
                    builtins, mathematics, angles);
            }
        }
    }

    FactorSplit left = splitConstantFactor(lhs, variable, builtins, mathematics, angles);
    FactorSplit right = splitConstantFactor(rhs, variable, builtins, mathematics, angles);
    if (left.dependent == right.dependent)
        return divide(builtins, mathematics, angles, left.constant, right.constant);
    if (!isOne(left.constant) || !isOne(right.constant)) {
        if (const auto dependentRatio = denominatorProportionalFactor(
                left.dependent, right.dependent, variable,
                builtins, mathematics, angles)) {
            return multiply(builtins, mathematics, angles, {
                divide(builtins, mathematics, angles, left.constant, right.constant),
                *dependentRatio});
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> proportionalFactor(
    const Expr& target,
    const Expr& reference,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Expr left = simplify(target, builtins, mathematics, angles);
    const Expr right = simplify(reference, builtins, mathematics, angles);
    if (left == right)
        return integer(1);
    if (isZero(right))
        return std::nullopt;

    if (isHead(left, builtins, BuiltinId::Divide)
        && isHead(right, builtins, BuiltinId::Divide)
        && left.asCall().arguments.size() == 2
        && right.asCall().arguments.size() == 2) {
        const auto numeratorRatio = proportionalFactor(
            left.asCall().arguments[0], right.asCall().arguments[0], variable,
            builtins, mathematics, angles);
        if (!numeratorRatio)
            return std::nullopt;

        if (left.asCall().arguments[1] == right.asCall().arguments[1])
            return numeratorRatio;
        if (const auto denominatorRatio = denominatorProportionalFactor(
                left.asCall().arguments[1], right.asCall().arguments[1], variable,
                builtins, mathematics, angles)) {
            return divide(builtins, mathematics, angles,
                *numeratorRatio, *denominatorRatio);
        }
    }

    FactorSplit lhs = splitConstantFactor(left, variable, builtins, mathematics, angles);
    FactorSplit rhs = splitConstantFactor(right, variable, builtins, mathematics, angles);
    Expr constantRatio = divide(
        builtins, mathematics, angles, lhs.constant, rhs.constant);
    if (containsVariable(constantRatio, variable))
        return std::nullopt;

    if (lhs.dependent == rhs.dependent)
        return constantRatio;

    // 定数因子を一段剥がしたことで式形が変わった場合だけ再帰する。
    // 1/(1+(x/2)^2)/2 のようなnested constant divisionを正規化し、
    // exact polynomial denominator ratioへ到達させる。
    if (!isOne(lhs.constant) || !isOne(rhs.constant)) {
        if (const auto dependentRatio = proportionalFactor(
                lhs.dependent, rhs.dependent, variable,
                builtins, mathematics, angles)) {
            return multiply(builtins, mathematics, angles,
                {std::move(constantRatio), *dependentRatio});
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> verifiedScaledPrimitive(
    const Expr& integrand,
    Expr candidate,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    candidate = simplify(std::move(candidate), builtins, mathematics, angles);
    Expr derivative = differentiateExpression(candidate, variable, builtins, mathematics, angles);
    const auto ratio = proportionalFactor(
        integrand, derivative, variable, builtins, mathematics, angles);
    if (!ratio)
        return std::nullopt;

    // derivativeとintegrandの比例関係は structural factor / exact polynomial ratio で
    // 既に証明済み。domain-safe Simplifierが x/x を消さない場合でも、この証明は有効。
    return multiply(builtins, mathematics, angles, {*ratio, std::move(candidate)});
}

[[nodiscard]] Expr pi(
    const mathematics::MathRegistry& mathematics) {
    const auto* definition = mathematics.findConstant(mathematics::ConstantId::Pi);
    return definition ? Expr{definition->symbol} : integer(0);
}

[[nodiscard]] Expr inverseAngleScale(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Radian:
        return integer(1);
    case mathematics::AngleUnit::Degree:
        return divide(builtins, mathematics, angles, integer(180), pi(mathematics));
    case mathematics::AngleUnit::Gradian:
        return divide(builtins, mathematics, angles, integer(200), pi(mathematics));
    }
    return integer(1);
}

[[nodiscard]] std::vector<Expr> primitiveTemplates(
    const Expr& u,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> candidates;
    candidates.reserve(40);

    const auto unary = [&](BuiltinId id) {
        return call(builtins, id, {u});
    };
    const Expr one = integer(1);
    const Expr two = integer(2);

    candidates.push_back(unary(BuiltinId::Log));
    candidates.push_back(unary(BuiltinId::Exp));
    candidates.push_back(unary(BuiltinId::Sin));
    candidates.push_back(unary(BuiltinId::Cos));
    candidates.push_back(unary(BuiltinId::Tan));
    candidates.push_back(unary(BuiltinId::Cot));
    candidates.push_back(unary(BuiltinId::Sinh));
    candidates.push_back(unary(BuiltinId::Cosh));
    candidates.push_back(unary(BuiltinId::Tanh));
    candidates.push_back(unary(BuiltinId::Coth));
    candidates.push_back(unary(BuiltinId::Sec));
    candidates.push_back(unary(BuiltinId::Csc));
    candidates.push_back(unary(BuiltinId::Sech));
    candidates.push_back(unary(BuiltinId::Csch));
    candidates.push_back(unary(BuiltinId::Sqrt));
    candidates.push_back(unary(BuiltinId::Cbrt));
    candidates.push_back(unary(BuiltinId::Log1p));
    candidates.push_back(unary(BuiltinId::Expm1));
    candidates.push_back(unary(BuiltinId::Erf));
    candidates.push_back(unary(BuiltinId::Erfc));
    candidates.push_back(unary(BuiltinId::Asin));
    candidates.push_back(unary(BuiltinId::Acos));
    candidates.push_back(unary(BuiltinId::Atan));
    candidates.push_back(unary(BuiltinId::Asinh));
    candidates.push_back(unary(BuiltinId::Acosh));
    candidates.push_back(unary(BuiltinId::Atanh));

    // sec/cscの標準原始函数。branchを含むため必ずDで検証してから採用する。
    candidates.push_back(call(builtins, BuiltinId::Log, {
        add(builtins, mathematics, angles,
            {unary(BuiltinId::Sec), unary(BuiltinId::Tan)})}));
    candidates.push_back(negate(builtins, mathematics, angles,
        call(builtins, BuiltinId::Log, {
            add(builtins, mathematics, angles,
                {unary(BuiltinId::Csc), unary(BuiltinId::Cot)})})));

    // erf/erfcの原始函数（u自身を変数とした標準形）。合成時は比例係数をDから決める。
    const Expr sqrtPi = call(builtins, BuiltinId::Sqrt, {pi(mathematics)});
    const Expr gaussian = call(builtins, BuiltinId::Exp, {
        negate(builtins, mathematics, angles,
            power(builtins, mathematics, angles, u, two))});
    candidates.push_back(add(builtins, mathematics, angles, {
        multiply(builtins, mathematics, angles, {u, unary(BuiltinId::Erf)}),
        divide(builtins, mathematics, angles, gaussian, sqrtPi)}));
    candidates.push_back(subtract(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles, {u, unary(BuiltinId::Erfc)}),
        divide(builtins, mathematics, angles, gaussian, sqrtPi)));

    // log(u)の原始函数。u'が定数の場合に比例検出で採用される。
    candidates.push_back(subtract(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles, {u, unary(BuiltinId::Log)}), u));

    // log(1+u) の標準原始函数。
    Expr onePlusU = add(builtins, mathematics, angles, {one, u});
    candidates.push_back(subtract(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles,
            {onePlusU, call(builtins, BuiltinId::Log, {onePlusU})}), u));

    // sqrt/cbrt/expm1 の標準原始函数。u'が定数または全体がchain-rule形ならD検証で採用する。
    candidates.push_back(multiply(builtins, mathematics, angles, {
        rational(Rational{BigInt{2}, BigInt{3}}), u, unary(BuiltinId::Sqrt)}));
    candidates.push_back(multiply(builtins, mathematics, angles, {
        rational(Rational{BigInt{3}, BigInt{4}}), u, unary(BuiltinId::Cbrt)}));
    candidates.push_back(subtract(builtins, mathematics, angles, unary(BuiltinId::Exp), u));

    // 逆三角函数の原始函数。出力角度単位のscaleを明示してRadian以外でも成立させる。
    const Expr angleScale = inverseAngleScale(builtins, mathematics, angles);
    const Expr uSquared = power(builtins, mathematics, angles, u, two);
    const Expr sqrtOneMinusSquare = call(builtins, BuiltinId::Sqrt, {
        subtract(builtins, mathematics, angles, one, uSquared)});
    candidates.push_back(add(builtins, mathematics, angles, {
        multiply(builtins, mathematics, angles, {u, unary(BuiltinId::Asin)}),
        multiply(builtins, mathematics, angles, {angleScale, sqrtOneMinusSquare})}));
    candidates.push_back(subtract(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles, {u, unary(BuiltinId::Acos)}),
        multiply(builtins, mathematics, angles, {angleScale, sqrtOneMinusSquare})));
    candidates.push_back(subtract(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles, {u, unary(BuiltinId::Atan)}),
        multiply(builtins, mathematics, angles, {
            divide(builtins, mathematics, angles, angleScale, two),
            call(builtins, BuiltinId::Log, {
                add(builtins, mathematics, angles, {one, uSquared})})})));

    // 逆双曲線函数。
    candidates.push_back(subtract(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles, {u, unary(BuiltinId::Asinh)}),
        call(builtins, BuiltinId::Sqrt, {
            add(builtins, mathematics, angles, {uSquared, one})})));
    candidates.push_back(add(builtins, mathematics, angles, {
        multiply(builtins, mathematics, angles, {u, unary(BuiltinId::Atanh)}),
        multiply(builtins, mathematics, angles, {
            rational(Rational{BigInt{1}, BigInt{2}}),
            call(builtins, BuiltinId::Log, {
                subtract(builtins, mathematics, angles, one, uSquared)})})}));

    // 積のreverse chain用。係数はD側から求めるので1/2等をここでhard-codeしない。
    candidates.push_back(power(builtins, mathematics, angles, unary(BuiltinId::Sin), two));
    candidates.push_back(power(builtins, mathematics, angles, unary(BuiltinId::Cos), two));
    candidates.push_back(power(builtins, mathematics, angles, unary(BuiltinId::Sinh), two));
    candidates.push_back(power(builtins, mathematics, angles, unary(BuiltinId::Cosh), two));

    return candidates;
}

[[nodiscard]] std::vector<Expr> dependentSubexpressions(
    const Expr& expression,
    const expression::Symbol& variable) {
    std::vector<Expr> result;
    std::vector<Expr> pending{expression};
    while (!pending.empty() && result.size() < maximumSubstitutionCandidates) {
        Expr current = std::move(pending.back());
        pending.pop_back();
        if (!containsVariable(current, variable))
            continue;

        const auto duplicate = std::find(result.begin(), result.end(), current);
        if (duplicate == result.end())
            result.push_back(current);

        if (current.isCall()) {
            for (const Expr& argument : current.asCall().arguments)
                pending.push_back(argument);
        }
    }

    const Expr variableExpr{variable};
    if (std::find(result.begin(), result.end(), variableExpr) == result.end())
        result.push_back(variableExpr);
    return result;
}

[[nodiscard]] std::optional<Expr> tryReverseChainRule(
    const Expr& integrand,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    for (const Expr& u : dependentSubexpressions(integrand, variable)) {
        for (Expr candidate : primitiveTemplates(u, builtins, mathematics, angles)) {
            if (auto result = verifiedScaledPrimitive(
                    integrand, std::move(candidate), variable,
                    builtins, mathematics, angles))
                return result;
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> integratePolynomial(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto polynomial = toExpressionPolynomial(
        expression, variable, builtins, mathematics, angles,
        PolynomialConversionOptions{256, 1024});
    if (!polynomial)
        return std::nullopt;

    std::vector<Expr> terms;
    const auto coefficients = polynomial->coefficients();
    terms.reserve(coefficients.size());
    for (std::size_t exponent = 0; exponent < coefficients.size(); ++exponent) {
        if (isZero(coefficients[exponent]))
            continue;
        const std::size_t resultExponent = exponent + 1;
        Expr variableTerm{variable};
        if (resultExponent != 1)
            variableTerm = power(builtins, mathematics, angles,
                std::move(variableTerm), rational(Rational{BigInt{static_cast<std::int64_t>(resultExponent)}}));
        Expr coefficient = divide(builtins, mathematics, angles,
            coefficients[exponent],
            rational(Rational{BigInt{static_cast<std::int64_t>(resultExponent)}}));
        terms.push_back(multiply(
            builtins, mathematics, angles,
            {std::move(coefficient), std::move(variableTerm)}));
    }
    if (terms.empty())
        return integer(0);
    return add(builtins, mathematics, angles, std::move(terms));
}

[[nodiscard]] std::optional<Expr> integratePowerRule(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Power)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const auto& a = expression.asCall().arguments;
    Expr exponentExpression = simplify(a[1], builtins, mathematics, angles);
    std::optional<Rational> exponent = exactRealRational(exponentExpression);
    if (!exponent && isHead(exponentExpression, builtins, BuiltinId::Negate)
        && exponentExpression.asCall().arguments.size() == 1) {
        if (const auto positive = exactRealRational(exponentExpression.asCall().arguments[0]))
            exponent = -*positive;
    }
    if (!exponent || containsVariable(a[1], variable))
        return std::nullopt;

    Expr baseDerivative = differentiateExpression(a[0], variable, builtins, mathematics, angles);
    if (containsVariable(baseDerivative, variable)
        || !provablyNonZero(baseDerivative, builtins, mathematics))
        return std::nullopt;

    if (*exponent == Rational{BigInt{-1}})
        return divide(builtins, mathematics, angles,
            call(builtins, BuiltinId::Log, {a[0]}), baseDerivative);

    const Rational next = *exponent + Rational{BigInt{1}};
    if (next.isZero())
        return std::nullopt;
    Expr denominator = multiply(
        builtins, mathematics, angles,
        {baseDerivative, rational(next)});
    return divide(builtins, mathematics, angles,
        power(builtins, mathematics, angles, a[0], rational(next)),
        std::move(denominator));
}

[[nodiscard]] Expr inverseAngleScale(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] Expr radiansPerInverseAngleUnit(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Radian:
        return integer(1);
    case mathematics::AngleUnit::Degree:
        return divide(builtins, mathematics, angles, pi(mathematics), integer(180));
    case mathematics::AngleUnit::Gradian:
        return divide(builtins, mathematics, angles, pi(mathematics), integer(200));
    }
    return integer(1);
}

[[nodiscard]] std::optional<std::array<Rational, 3>> exactQuadraticCoefficients(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    const auto polynomial = toRationalPolynomial(
        expression, variable, builtins, PolynomialConversionOptions{32, 128});
    if (!polynomial || polynomial->degree() != 2)
        return std::nullopt;
    return std::array<Rational, 3>{
        polynomial->coefficient(0), polynomial->coefficient(1), polynomial->coefficient(2)};
}

[[nodiscard]] Expr quadraticLinearNumerator(
    const std::array<Rational, 3>& coefficients,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return add(builtins, mathematics, angles, {
        multiply(builtins, mathematics, angles, {
            rational(Rational{BigInt{2}} * coefficients[2]), Expr{variable}}),
        rational(coefficients[1])});
}

[[nodiscard]] std::optional<Expr> tryQuadraticReciprocal(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr numerator = integer(1);
    std::optional<Expr> denominator;
    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2
        && !containsVariable(expression.asCall().arguments[0], variable)) {
        numerator = expression.asCall().arguments[0];
        denominator = expression.asCall().arguments[1];
    }
    else if (isHead(expression, builtins, BuiltinId::Power)
        && expression.asCall().arguments.size() == 2
        && exactRealRational(expression.asCall().arguments[1]) == Rational{BigInt{-1}}) {
        denominator = expression.asCall().arguments[0];
    }
    else {
        return std::nullopt;
    }

    const auto coefficients = exactQuadraticCoefficients(*denominator, variable, builtins);
    if (!coefficients || (*coefficients)[2].isZero())
        return std::nullopt;

    const Rational discriminant = (*coefficients)[1] * (*coefficients)[1]
        - Rational{BigInt{4}} * (*coefficients)[2] * (*coefficients)[0];
    Expr linear = quadraticLinearNumerator(
        *coefficients, variable, builtins, mathematics, angles);

    if (discriminant < Rational{BigInt{0}}) {
        Expr scale = call(builtins, BuiltinId::Sqrt, {rational(-discriminant)});
        Expr argument = divide(
            builtins, mathematics, angles, std::move(linear), scale);
        Expr coefficient = divide(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {
                integer(2), radiansPerInverseAngleUnit(builtins, mathematics, angles)}),
            scale);
        return multiply(builtins, mathematics, angles, {
            std::move(numerator), std::move(coefficient),
            call(builtins, BuiltinId::Atan, {std::move(argument)})});
    }
    if (discriminant > Rational{BigInt{0}}) {
        Expr scale = call(builtins, BuiltinId::Sqrt, {rational(discriminant)});
        const bool reverse = (*coefficients)[2] < Rational{BigInt{0}};
        if (reverse)
            linear = negate(builtins, mathematics, angles, std::move(linear));
        Expr argument = divide(
            builtins, mathematics, angles, std::move(linear), scale);
        Expr coefficient = divide(builtins, mathematics, angles,
            reverse ? integer(2) : integer(-2), scale);
        return multiply(builtins, mathematics, angles, {
            std::move(numerator), std::move(coefficient),
            call(builtins, BuiltinId::Atanh, {std::move(argument)})});
    }

    Expr halfLinear = add(builtins, mathematics, angles, {
        multiply(builtins, mathematics, angles, {
            rational((*coefficients)[2]), Expr{variable}}),
        rational((*coefficients)[1] / Rational{BigInt{2}})});
    return divide(builtins, mathematics, angles,
        negate(builtins, mathematics, angles, std::move(numerator)),
        std::move(halfLinear));
}

struct PolynomialDivision final {
    RationalPolynomial quotient;
    RationalPolynomial remainder;
};

[[nodiscard]] PolynomialDivision dividePolynomials(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator) {
    std::vector<Rational> remainder = numerator.coefficients();
    std::vector<Rational> quotient(
        numerator.degree() >= denominator.degree()
            ? numerator.degree() - denominator.degree() + 1
            : 1,
        Rational{BigInt{0}});

    const auto trim = [](std::vector<Rational>& coefficients) {
        while (coefficients.size() > 1 && coefficients.back().isZero())
            coefficients.pop_back();
        if (coefficients.empty())
            coefficients.push_back(Rational{BigInt{0}});
    };
    trim(remainder);

    const Rational leading = denominator.coefficient(denominator.degree());
    while (!(remainder.size() == 1 && remainder[0].isZero())
        && remainder.size() - 1 >= denominator.degree()) {
        const std::size_t degree = remainder.size() - 1;
        const std::size_t shift = degree - denominator.degree();
        const Rational factor = remainder.back() / leading;
        quotient[shift] += factor;
        for (std::size_t i = 0; i <= denominator.degree(); ++i)
            remainder[i + shift] -= factor * denominator.coefficient(i);
        trim(remainder);
    }

    return PolynomialDivision{
        RationalPolynomial{std::move(quotient)},
        RationalPolynomial{std::move(remainder)}};
}


[[nodiscard]] RationalPolynomial multiplyPolynomials(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    if (lhs.isZero() || rhs.isZero())
        return RationalPolynomial{};
    std::vector<Rational> coefficients(
        lhs.degree() + rhs.degree() + 1, Rational{BigInt{0}});
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        for (std::size_t j = 0; j <= rhs.degree(); ++j)
            coefficients[i + j] += lhs.coefficient(i) * rhs.coefficient(j);
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] RationalPolynomial powerPolynomial(
    RationalPolynomial base,
    std::size_t exponent) {
    RationalPolynomial result{{Rational{BigInt{1}}}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = multiplyPolynomials(result, base);
        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyPolynomials(base, base);
    }
    return result;
}

[[nodiscard]] RationalPolynomial multiplyPolynomialByMonomial(
    const RationalPolynomial& polynomial,
    std::size_t exponent) {
    std::vector<Rational> coefficients(exponent, Rational{BigInt{0}});
    const auto source = polynomial.coefficients();
    coefficients.insert(coefficients.end(), source.begin(), source.end());
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] std::optional<Rational> perfectRationalSquareRoot(const Rational& value) {
    if (value < Rational{BigInt{0}})
        return std::nullopt;
    const auto numerator = numeric::integerSqrt(value.numerator());
    const auto denominator = numeric::integerSqrt(value.denominator());
    if (!numerator.remainder.isZero() || !denominator.remainder.isZero())
        return std::nullopt;
    return Rational{numerator.root, denominator.root};
}

struct RationalFactor final {
    RationalPolynomial polynomial;
    std::size_t multiplicity = 1;
};

[[nodiscard]] std::optional<std::vector<RationalFactor>> factorIntoLinearQuadratic(
    RationalPolynomial polynomial) {
    std::vector<RationalFactor> factors;
    while (polynomial.degree() > 0) {
        std::optional<Rational> root;
        if (polynomial.degree() == 1) {
            root = -polynomial.coefficient(0) / polynomial.coefficient(1);
        }
        else if (polynomial.degree() == 2) {
            const Rational a = polynomial.coefficient(2);
            const Rational b = polynomial.coefficient(1);
            const Rational c = polynomial.coefficient(0);
            const Rational discriminant = b * b - Rational{BigInt{4}} * a * c;
            if (const auto squareRoot = perfectRationalSquareRoot(discriminant))
                root = (-b + *squareRoot) / (Rational{BigInt{2}} * a);
            else {
                factors.push_back(RationalFactor{polynomial, 1});
                break;
            }
        }
        else {
            const RationalRootSearchResult search = findRationalRoot(
                polynomial, RationalRootSearchOptions{1'000'000, 200'000});
            if (search.root)
                root = *search.root;
            else
                return std::nullopt;
        }

        RationalPolynomial linear{{-*root, Rational{BigInt{1}}}};
        std::size_t multiplicity = 0;
        while (true) {
            const auto quotient = divideByLinearFactor(polynomial, *root);
            if (!quotient)
                break;
            polynomial = *quotient;
            ++multiplicity;
            if (polynomial.degree() == 0)
                break;
        }
        if (multiplicity == 0)
            return std::nullopt;
        factors.push_back(RationalFactor{std::move(linear), multiplicity});
    }
    return factors;
}

struct PartialFractionBasis final {
    std::size_t factorIndex = 0;
    std::size_t denominatorPower = 1;
    std::size_t numeratorExponent = 0;
};

[[nodiscard]] std::optional<std::vector<Rational>> solveRationalLinearSystem(
    std::vector<std::vector<Rational>> matrix,
    std::vector<Rational> rhs) {
    const std::size_t n = matrix.size();
    if (rhs.size() != n)
        return std::nullopt;
    for (const auto& row : matrix)
        if (row.size() != n)
            return std::nullopt;

    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        while (pivot < n && matrix[pivot][column].isZero())
            ++pivot;
        if (pivot == n)
            return std::nullopt;
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            std::swap(rhs[pivot], rhs[column]);
        }

        const Rational pivotValue = matrix[column][column];
        for (std::size_t j = column; j < n; ++j)
            matrix[column][j] /= pivotValue;
        rhs[column] /= pivotValue;

        for (std::size_t row = 0; row < n; ++row) {
            if (row == column || matrix[row][column].isZero())
                continue;
            const Rational factor = matrix[row][column];
            for (std::size_t j = column; j < n; ++j)
                matrix[row][j] -= factor * matrix[column][j];
            rhs[row] -= factor * rhs[column];
        }
    }
    return rhs;
}

[[nodiscard]] std::optional<Expr> tryRationalLowDegree(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] std::optional<Expr> tryRationalFactored(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const Expr& numeratorExpr = expression.asCall().arguments[0];
    const Expr& denominatorExpr = expression.asCall().arguments[1];
    const auto numerator = toRationalPolynomial(
        numeratorExpr, variable, builtins, PolynomialConversionOptions{256, 2048});
    const auto denominator = toRationalPolynomial(
        denominatorExpr, variable, builtins, PolynomialConversionOptions{64, 1024});
    if (!numerator || !denominator || denominator->degree() <= 2 || denominator->degree() > 12)
        return std::nullopt;

    PolynomialDivision division = dividePolynomials(*numerator, *denominator);
    const auto factors = factorIntoLinearQuadratic(*denominator);
    if (!factors)
        return std::nullopt;

    std::vector<PartialFractionBasis> basis;
    for (std::size_t i = 0; i < factors->size(); ++i) {
        const std::size_t degree = (*factors)[i].polynomial.degree();
        if (degree == 0 || degree > 2)
            return std::nullopt;
        if (degree == 2 && (*factors)[i].multiplicity != 1)
            return std::nullopt;
        for (std::size_t k = 1; k <= (*factors)[i].multiplicity; ++k)
            for (std::size_t j = 0; j < degree; ++j)
                basis.push_back(PartialFractionBasis{i, k, j});
    }

    const std::size_t unknowns = basis.size();
    if (unknowns != denominator->degree() || unknowns == 0)
        return std::nullopt;

    std::vector<std::vector<Rational>> matrix(
        unknowns, std::vector<Rational>(unknowns, Rational{BigInt{0}}));
    for (std::size_t column = 0; column < unknowns; ++column) {
        const auto& descriptor = basis[column];
        const RationalPolynomial divisor = powerPolynomial(
            (*factors)[descriptor.factorIndex].polynomial,
            descriptor.denominatorPower);
        PolynomialDivision quotient = dividePolynomials(*denominator, divisor);
        if (!quotient.remainder.isZero())
            return std::nullopt;
        RationalPolynomial contribution = multiplyPolynomialByMonomial(
            quotient.quotient, descriptor.numeratorExponent);
        for (std::size_t row = 0; row < unknowns; ++row)
            matrix[row][column] = contribution.coefficient(row);
    }

    std::vector<Rational> rhs(unknowns, Rational{BigInt{0}});
    for (std::size_t row = 0; row < unknowns; ++row)
        rhs[row] = division.remainder.coefficient(row);
    const auto coefficients = solveRationalLinearSystem(std::move(matrix), std::move(rhs));
    if (!coefficients)
        return std::nullopt;

    std::vector<Expr> primitiveTerms;
    if (!division.quotient.isZero()) {
        Expr quotientExpr = polynomialToExpandedExpr(division.quotient, variable, builtins);
        const auto primitive = integratePolynomial(
            quotientExpr, variable, builtins, mathematics, angles);
        if (!primitive)
            return std::nullopt;
        primitiveTerms.push_back(*primitive);
    }

    for (std::size_t factorIndex = 0; factorIndex < factors->size(); ++factorIndex) {
        const RationalFactor& factor = (*factors)[factorIndex];
        const Expr factorExpr = polynomialToExpandedExpr(factor.polynomial, variable, builtins);
        if (factor.polynomial.degree() == 1) {
            for (std::size_t k = 1; k <= factor.multiplicity; ++k) {
                Rational coefficient{BigInt{0}};
                for (std::size_t column = 0; column < basis.size(); ++column)
                    if (basis[column].factorIndex == factorIndex
                        && basis[column].denominatorPower == k)
                        coefficient += (*coefficients)[column];
                if (coefficient.isZero())
                    continue;
                if (k == 1) {
                    primitiveTerms.push_back(multiply(builtins, mathematics, angles, {
                        rational(coefficient), call(builtins, BuiltinId::Log, {factorExpr})}));
                }
                else {
                    const Rational exponent{BigInt{1 - static_cast<std::int64_t>(k)}};
                    primitiveTerms.push_back(multiply(builtins, mathematics, angles, {
                        rational(coefficient / exponent),
                        power(builtins, mathematics, angles, factorExpr, rational(exponent))}));
                }
            }
            continue;
        }

        std::vector<Rational> numeratorCoefficients(2, Rational{BigInt{0}});
        for (std::size_t column = 0; column < basis.size(); ++column) {
            if (basis[column].factorIndex == factorIndex
                && basis[column].denominatorPower == 1)
                numeratorCoefficients[basis[column].numeratorExponent] += (*coefficients)[column];
        }
        RationalPolynomial partialNumerator{std::move(numeratorCoefficients)};
        if (partialNumerator.isZero())
            continue;
        Expr partial = divide(
            builtins, mathematics, angles,
            polynomialToExpandedExpr(partialNumerator, variable, builtins), factorExpr);
        const auto primitive = tryRationalLowDegree(
            partial, variable, builtins, mathematics, angles);
        if (!primitive)
            return std::nullopt;
        primitiveTerms.push_back(*primitive);
    }

    if (primitiveTerms.empty())
        return integer(0);
    // 係数は denominator * partialFraction == remainder となるexact Rational
    // linear systemを解いて得ている。これは数値推測ではなく代数恒等式なので、そのまま採用する。
    return add(builtins, mathematics, angles, std::move(primitiveTerms));
}

[[nodiscard]] std::optional<Expr> tryRationalLowDegree(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const Expr& numeratorExpr = expression.asCall().arguments[0];
    const Expr& denominatorExpr = expression.asCall().arguments[1];
    const auto numerator = toRationalPolynomial(
        numeratorExpr, variable, builtins, PolynomialConversionOptions{128, 512});
    const auto denominator = toRationalPolynomial(
        denominatorExpr, variable, builtins, PolynomialConversionOptions{32, 128});
    if (!numerator || !denominator
        || denominator->degree() == 0 || denominator->degree() > 2)
        return std::nullopt;

    PolynomialDivision division = dividePolynomials(*numerator, *denominator);
    std::vector<Expr> primitiveTerms;

    if (!division.quotient.isZero()) {
        Expr quotientExpr = polynomialToExpandedExpr(
            division.quotient, variable, builtins);
        const auto quotientPrimitive = integratePolynomial(
            quotientExpr, variable, builtins, mathematics, angles);
        if (!quotientPrimitive)
            return std::nullopt;
        primitiveTerms.push_back(*quotientPrimitive);
    }

    if (denominator->degree() == 1) {
        const Rational slope = denominator->coefficient(1);
        if (slope.isZero())
            return std::nullopt;
        const Rational remainder = division.remainder.coefficient(0);
        if (!remainder.isZero()) {
            primitiveTerms.push_back(multiply(builtins, mathematics, angles, {
                rational(remainder / slope),
                call(builtins, BuiltinId::Log, {denominatorExpr})}));
        }
    }
    else {
        const Rational a = denominator->coefficient(2);
        if (a.isZero())
            return std::nullopt;
        const Rational linear = division.remainder.coefficient(1);
        const Rational constant = division.remainder.coefficient(0);
        const Rational logCoefficient = linear / (Rational{BigInt{2}} * a);
        const Rational reciprocalCoefficient = constant
            - logCoefficient * denominator->coefficient(1);

        if (!logCoefficient.isZero()) {
            primitiveTerms.push_back(multiply(builtins, mathematics, angles, {
                rational(logCoefficient),
                call(builtins, BuiltinId::Log, {denominatorExpr})}));
        }

        if (!reciprocalCoefficient.isZero()) {
            Expr reciprocal = divide(
                builtins, mathematics, angles, integer(1), denominatorExpr);
            const auto primitive = tryQuadraticReciprocal(
                reciprocal, variable, builtins, mathematics, angles);
            if (!primitive)
                return std::nullopt;
            primitiveTerms.push_back(multiply(builtins, mathematics, angles, {
                rational(reciprocalCoefficient), *primitive}));
        }
    }

    if (primitiveTerms.empty())
        return integer(0);
    return add(builtins, mathematics, angles, std::move(primitiveTerms));
}


[[nodiscard]] std::optional<Expr> rewriteForSquareRootSubstitution(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    bool& sawPrincipalRoot) {
    if (expression.isSymbol() && expression.asSymbol().sameIdentity(variable))
        return power(builtins, mathematics, angles, Expr{variable}, integer(2));

    if (!expression.isCall())
        return expression;

    if (isHead(expression, builtins, BuiltinId::Sqrt)
        && expression.asCall().arguments.size() == 1
        && expression.asCall().arguments[0].isSymbol()
        && expression.asCall().arguments[0].asSymbol().sameIdentity(variable)) {
        sawPrincipalRoot = true;
        return Expr{variable};
    }

    std::vector<Expr> arguments;
    arguments.reserve(expression.asCall().arguments.size());
    for (const Expr& argument : expression.asCall().arguments) {
        auto rewritten = rewriteForSquareRootSubstitution(
            argument, variable, builtins, mathematics, angles, sawPrincipalRoot);
        if (!rewritten)
            return std::nullopt;
        arguments.push_back(std::move(*rewritten));
    }
    return simplify(Expr::call(expression.asCall().head, std::move(arguments)),
        builtins, mathematics, angles);
}

// sqrt[q(sqrt[x])] で q がexact Rational係数の2次式なら、t=sqrt[x] により
// 2 t sqrt[q(t)] dt へ落とす。一般的な置換探索ではなく、この代数的class全体を
// 一つの局所置換則として扱う。principal sqrt/logのglobal identityは追加しない。
[[nodiscard]] std::optional<Expr> tryNestedQuadraticSquareRoot(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Sqrt)
        || expression.asCall().arguments.size() != 1)
        return std::nullopt;

    bool sawPrincipalRoot = false;
    auto transformed = rewriteForSquareRootSubstitution(
        expression.asCall().arguments[0], variable,
        builtins, mathematics, angles, sawPrincipalRoot);
    if (!transformed || !sawPrincipalRoot)
        return std::nullopt;

    const auto coefficients = exactQuadraticCoefficients(*transformed, variable, builtins);
    if (!coefficients)
        return std::nullopt;

    const Rational c = (*coefficients)[0];
    const Rational b = (*coefficients)[1];
    const Rational a = (*coefficients)[2];
    if (!(a > Rational{BigInt{0}}))
        return std::nullopt;

    const Expr t{variable};
    const Expr q = polynomialToExpandedExpr(
        RationalPolynomial{std::vector<Rational>{c, b, a}}, variable, builtins);
    const Expr sqrtQ = call(builtins, BuiltinId::Sqrt, {q});
    const Expr sqrtA = call(builtins, BuiltinId::Sqrt, {rational(a)});
    const Expr linear = add(builtins, mathematics, angles, {
        multiply(builtins, mathematics, angles,
            {rational(Rational{BigInt{2}} * a), t}),
        rational(b)});

    // J = integral sqrt[a t^2+b t+c] dt
    //   = (2at+b)sqrt[q]/(4a)
    //     - (b^2-4ac)/(8 a sqrt[a]) Log[2sqrt[a]sqrt[q]+2at+b]
    // はprincipal branches上の局所原始函数。定数差は表示しない。
    const Rational discriminant = b * b - Rational{BigInt{4}} * a * c;
    Expr j = divide(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles, {linear, sqrtQ}),
        rational(Rational{BigInt{4}} * a));
    if (!discriminant.isZero()) {
        Expr logArgument = add(builtins, mathematics, angles, {
            multiply(builtins, mathematics, angles,
                {integer(2), sqrtA, sqrtQ}),
            linear});
        Expr logCoefficient = divide(builtins, mathematics, angles,
            rational(-discriminant),
            multiply(builtins, mathematics, angles,
                {integer(8), rational(a), sqrtA}));
        j = add(builtins, mathematics, angles, {
            std::move(j),
            multiply(builtins, mathematics, angles,
                {std::move(logCoefficient), call(builtins, BuiltinId::Log, {std::move(logArgument)})})});
    }

    // 2t = q'(t)/a - b/a。
    Expr primitiveInT = add(builtins, mathematics, angles, {
        multiply(builtins, mathematics, angles, {
            rational(Rational{BigInt{2}} / (Rational{BigInt{3}} * a)),
            q, sqrtQ}),
        multiply(builtins, mathematics, angles, {
            rational(-b / a), std::move(j)})});

    const Expr rootX = call(builtins, BuiltinId::Sqrt, {Expr{variable}});
    Expr result = substituteSymbol(primitiveInT, variable, rootX);
    return fullSimplify(std::move(result), builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> tryQuadraticInverseSqrt(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr numerator = integer(1);
    std::optional<Expr> radicand;
    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2
        && !containsVariable(expression.asCall().arguments[0], variable)
        && isHead(expression.asCall().arguments[1], builtins, BuiltinId::Sqrt)
        && expression.asCall().arguments[1].asCall().arguments.size() == 1) {
        numerator = expression.asCall().arguments[0];
        radicand = expression.asCall().arguments[1].asCall().arguments[0];
    }
    else if (isHead(expression, builtins, BuiltinId::Power)
        && expression.asCall().arguments.size() == 2
        && exactRealRational(expression.asCall().arguments[1])
            == Rational{BigInt{-1}, BigInt{2}}) {
        radicand = expression.asCall().arguments[0];
    }
    else {
        return std::nullopt;
    }

    const auto coefficients = exactQuadraticCoefficients(*radicand, variable, builtins);
    if (!coefficients || (*coefficients)[2].isZero())
        return std::nullopt;

    const Rational a = (*coefficients)[2];
    const Rational discriminant = (*coefficients)[1] * (*coefficients)[1]
        - Rational{BigInt{4}} * a * (*coefficients)[0];
    Expr linear = quadraticLinearNumerator(
        *coefficients, variable, builtins, mathematics, angles);

    if (a > Rational{BigInt{0}} && discriminant < Rational{BigInt{0}}) {
        Expr scale = call(builtins, BuiltinId::Sqrt, {rational(-discriminant)});
        Expr argument = divide(
            builtins, mathematics, angles, std::move(linear), std::move(scale));
        Expr leadingRoot = call(builtins, BuiltinId::Sqrt, {rational(a)});
        return multiply(builtins, mathematics, angles, {
            std::move(numerator),
            divide(builtins, mathematics, angles,
                call(builtins, BuiltinId::Asinh, {std::move(argument)}),
                std::move(leadingRoot))});
    }
    if (a > Rational{BigInt{0}} && discriminant > Rational{BigInt{0}}) {
        // 実根を持つ正leading quadratic。大域的なsqrt因数分解は行わず、
        // principal式の局所原始函数候補を直接作り、Dによるexact検証に通す。
        Expr sqrtA = call(builtins, BuiltinId::Sqrt, {rational(a)});
        Expr inner = add(builtins, mathematics, angles, {
            multiply(builtins, mathematics, angles, {
                sqrtA, call(builtins, BuiltinId::Sqrt, {*radicand})}),
            multiply(builtins, mathematics, angles, {rational(a), Expr{variable}}),
            rational((*coefficients)[1] / Rational{BigInt{2}})});
        Expr candidate = multiply(builtins, mathematics, angles, {
            numerator,
            divide(builtins, mathematics, angles,
                call(builtins, BuiltinId::Log, {std::move(inner)}), sqrtA)});
        // これはsqrt因数分解ではなく、quadratic completionから直接得られる
        // principal式の局所原始函数。共通解析領域上の代表元として採用する。
        return candidate;
    }
    if (a < Rational{BigInt{0}} && discriminant > Rational{BigInt{0}}) {
        Expr scale = call(builtins, BuiltinId::Sqrt, {rational(discriminant)});
        linear = negate(builtins, mathematics, angles, std::move(linear));
        Expr argument = divide(
            builtins, mathematics, angles, std::move(linear), std::move(scale));
        Expr leadingRoot = call(builtins, BuiltinId::Sqrt, {rational(-a)});
        Expr coefficient = divide(builtins, mathematics, angles,
            radiansPerInverseAngleUnit(builtins, mathematics, angles),
            std::move(leadingRoot));
        return multiply(builtins, mathematics, angles, {
            std::move(numerator), std::move(coefficient),
            call(builtins, BuiltinId::Asin, {std::move(argument)})});
    }
    return std::nullopt;
}

struct TrigArgument final {
    Expr argument;
    Expr scale;
    Expr inverseScale;
};

[[nodiscard]] std::optional<mathematics::AngleUnit> explicitAngleUnit(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::UnitApplied))
        return std::nullopt;
    const auto& a = expression.asCall().arguments;
    if (a.size() != 2 || !a[1].isString())
        return std::nullopt;
    return mathematics::AngleSemantics::parseUnit(a[1].asString());
}

[[nodiscard]] Expr directTrigScale(
    mathematics::AngleUnit unit,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (unit) {
    case mathematics::AngleUnit::Radian:
        return integer(1);
    case mathematics::AngleUnit::Degree:
        return divide(builtins, mathematics, angles, pi(mathematics), integer(180));
    case mathematics::AngleUnit::Gradian:
        return divide(builtins, mathematics, angles, pi(mathematics), integer(200));
    }
    return integer(1);
}

[[nodiscard]] Expr directTrigInverseScale(
    mathematics::AngleUnit unit,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (unit) {
    case mathematics::AngleUnit::Radian:
        return integer(1);
    case mathematics::AngleUnit::Degree:
        return divide(builtins, mathematics, angles, integer(180), pi(mathematics));
    case mathematics::AngleUnit::Gradian:
        return divide(builtins, mathematics, angles, integer(200), pi(mathematics));
    }
    return integer(1);
}

[[nodiscard]] TrigArgument trigArgument(
    const Expr& source,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    mathematics::AngleUnit unit = angles.defaultUnit();
    Expr argument = source;
    if (const auto explicitUnit = explicitAngleUnit(source, builtins)) {
        unit = *explicitUnit;
        argument = source.asCall().arguments[0];
    }
    return TrigArgument{
        std::move(argument),
        directTrigScale(unit, builtins, mathematics, angles),
        directTrigInverseScale(unit, builtins, mathematics, angles)};
}

[[nodiscard]] Expr doubledTrigArgument(
    const Expr& source,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (isHead(source, builtins, BuiltinId::UnitApplied)
        && source.asCall().arguments.size() == 2) {
        return call(builtins, BuiltinId::UnitApplied, {
            multiply(builtins, mathematics, angles,
                {integer(2), source.asCall().arguments[0]}),
            source.asCall().arguments[1]});
    }
    return multiply(builtins, mathematics, angles, {integer(2), source});
}

[[nodiscard]] std::optional<Expr> rewriteSquareIdentity(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Power)
        || expression.asCall().arguments.size() != 2
        || exactRealRational(expression.asCall().arguments[1]) != Rational{BigInt{2}})
        return std::nullopt;

    const Expr& base = expression.asCall().arguments[0];
    if (!base.isCall() || base.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = builtins.find(base.asCall().head);
    if (!definition)
        return std::nullopt;

    const Expr& u = base.asCall().arguments[0];
    const Expr one = integer(1);
    const Expr two = integer(2);
    const auto squareOf = [&](BuiltinId id) {
        return power(builtins, mathematics, angles,
            call(builtins, id, {u}), two);
    };

    switch (definition->id) {
    case BuiltinId::Sin:
        return divide(builtins, mathematics, angles,
            subtract(builtins, mathematics, angles, one,
                call(builtins, BuiltinId::Cos, {
                    doubledTrigArgument(u, builtins, mathematics, angles)})), two);
    case BuiltinId::Cos:
        return divide(builtins, mathematics, angles,
            add(builtins, mathematics, angles, {one,
                call(builtins, BuiltinId::Cos, {
                    doubledTrigArgument(u, builtins, mathematics, angles)})}), two);
    case BuiltinId::Tan:
        return subtract(builtins, mathematics, angles, squareOf(BuiltinId::Sec), one);
    case BuiltinId::Cot:
        return subtract(builtins, mathematics, angles, squareOf(BuiltinId::Csc), one);
    case BuiltinId::Sinh:
        return divide(builtins, mathematics, angles,
            subtract(builtins, mathematics, angles,
                call(builtins, BuiltinId::Cosh, {
                    multiply(builtins, mathematics, angles, {two, u})}), one), two);
    case BuiltinId::Cosh:
        return divide(builtins, mathematics, angles,
            add(builtins, mathematics, angles, {
                call(builtins, BuiltinId::Cosh, {
                    multiply(builtins, mathematics, angles, {two, u})}), one}), two);
    case BuiltinId::Tanh:
        return subtract(builtins, mathematics, angles, one, squareOf(BuiltinId::Sech));
    case BuiltinId::Coth:
        return add(builtins, mathematics, angles, {one, squareOf(BuiltinId::Csch)});
    default:
        return std::nullopt;
    }
}

[[nodiscard]] std::optional<Expr> integrateElementaryUnary(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return std::nullopt;

    const Expr& sourceArgument = expression.asCall().arguments[0];
    Expr differentiationArgument = sourceArgument;
    Expr derivativeScale = integer(1);
    Expr inverseDerivativeScale = integer(1);

    switch (definition->id) {
    case BuiltinId::Sin:
    case BuiltinId::Cos:
    case BuiltinId::Tan:
    case BuiltinId::Cot:
    case BuiltinId::Sec:
    case BuiltinId::Csc: {
        TrigArgument info = trigArgument(
            sourceArgument, builtins, mathematics, angles);
        differentiationArgument = std::move(info.argument);
        derivativeScale = std::move(info.scale);
        inverseDerivativeScale = std::move(info.inverseScale);
        break;
    }
    default:
        break;
    }

    Expr du = differentiateExpression(
        differentiationArgument, variable, builtins, mathematics, angles);
    if (containsVariable(du, variable) || !provablyNonZero(du, builtins, mathematics))
        return std::nullopt;
    Expr totalScale = multiply(
        builtins, mathematics, angles, {derivativeScale, du});
    if (!provablyNonZero(totalScale, builtins, mathematics))
        return std::nullopt;
    inverseDerivativeScale = divide(
        builtins, mathematics, angles, std::move(inverseDerivativeScale), du);

    const auto scaledTrig = [&](Expr primitive) {
        return multiply(builtins, mathematics, angles,
            {inverseDerivativeScale, std::move(primitive)});
    };
    const auto same = [&](BuiltinId id) {
        return call(builtins, id, {sourceArgument});
    };
    switch (definition->id) {
    case BuiltinId::Exp:
        return divide(builtins, mathematics, angles, expression, du);
    case BuiltinId::Sin:
        return scaledTrig(negate(builtins, mathematics, angles, same(BuiltinId::Cos)));
    case BuiltinId::Cos:
        return scaledTrig(same(BuiltinId::Sin));
    case BuiltinId::Tan:
        return scaledTrig(negate(builtins, mathematics, angles,
            call(builtins, BuiltinId::Log, {same(BuiltinId::Cos)})));
    case BuiltinId::Cot:
        return scaledTrig(call(builtins, BuiltinId::Log, {same(BuiltinId::Sin)}));
    case BuiltinId::Sec:
        return scaledTrig(call(builtins, BuiltinId::Log, {
            add(builtins, mathematics, angles, {same(BuiltinId::Sec), same(BuiltinId::Tan)})}));
    case BuiltinId::Csc:
        return scaledTrig(negate(builtins, mathematics, angles,
            call(builtins, BuiltinId::Log, {
                add(builtins, mathematics, angles, {same(BuiltinId::Csc), same(BuiltinId::Cot)})})));
    case BuiltinId::Sinh:
        return divide(builtins, mathematics, angles, same(BuiltinId::Cosh), du);
    case BuiltinId::Cosh:
        return divide(builtins, mathematics, angles, same(BuiltinId::Sinh), du);
    case BuiltinId::Tanh:
        return divide(builtins, mathematics, angles,
            call(builtins, BuiltinId::Log, {same(BuiltinId::Cosh)}), du);
    case BuiltinId::Coth:
        return divide(builtins, mathematics, angles,
            call(builtins, BuiltinId::Log, {same(BuiltinId::Sinh)}), du);
    case BuiltinId::Sech:
        return divide(builtins, mathematics, angles,
            call(builtins, BuiltinId::Atan, {same(BuiltinId::Sinh)}),
            multiply(builtins, mathematics, angles, {
                inverseAngleScale(builtins, mathematics, angles), du}));
    case BuiltinId::Csch:
        return divide(builtins, mathematics, angles,
            call(builtins, BuiltinId::Log, {
                call(builtins, BuiltinId::Tanh, {
                    divide(builtins, mathematics, angles, sourceArgument, integer(2))})}),
            du);
    default:
        return std::nullopt;
    }
}

[[nodiscard]] std::optional<Expr> integrateStandardUnary(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return std::nullopt;

    const Expr& u = expression.asCall().arguments[0];
    Expr du = differentiateExpression(u, variable, builtins, mathematics, angles);
    if (containsVariable(du, variable) || !provablyNonZero(du, builtins, mathematics))
        return std::nullopt;

    const Expr one = integer(1);
    const Expr two = integer(2);
    const Expr u2 = power(builtins, mathematics, angles, u, two);
    const auto divideByDu = [&](Expr numerator) {
        return divide(builtins, mathematics, angles, std::move(numerator), du);
    };

    switch (definition->id) {
    case BuiltinId::Log:
        return divideByDu(subtract(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {u, expression}), u));

    case BuiltinId::Log1p: {
        Expr onePlusU = add(builtins, mathematics, angles, {one, u});
        return divideByDu(subtract(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {
                onePlusU, call(builtins, BuiltinId::Log, {onePlusU})}), u));
    }

    case BuiltinId::Sqrt:
        return divideByDu(multiply(builtins, mathematics, angles, {
            rational(Rational{BigInt{2}, BigInt{3}}), u, expression}));

    case BuiltinId::Cbrt:
        return divideByDu(multiply(builtins, mathematics, angles, {
            rational(Rational{BigInt{3}, BigInt{4}}), u, expression}));

    case BuiltinId::Expm1:
        return divideByDu(subtract(builtins, mathematics, angles,
            call(builtins, BuiltinId::Exp, {u}), u));

    case BuiltinId::Erf:
    case BuiltinId::Erfc: {
        Expr gaussianOverSqrtPi = divide(builtins, mathematics, angles,
            call(builtins, BuiltinId::Exp, {
                negate(builtins, mathematics, angles, u2)}),
            call(builtins, BuiltinId::Sqrt, {pi(mathematics)}));
        Expr main = multiply(builtins, mathematics, angles, {u, expression});
        return divideByDu(definition->id == BuiltinId::Erf
            ? add(builtins, mathematics, angles, {std::move(main), std::move(gaussianOverSqrtPi)})
            : subtract(builtins, mathematics, angles, std::move(main), std::move(gaussianOverSqrtPi)));
    }

    case BuiltinId::Asin:
    case BuiltinId::Acos:
    case BuiltinId::Atan: {
        const Expr scale = inverseAngleScale(builtins, mathematics, angles);
        Expr main = multiply(builtins, mathematics, angles, {u, expression});
        Expr correction = integer(0);
        if (definition->id == BuiltinId::Atan) {
            correction = multiply(builtins, mathematics, angles, {
                divide(builtins, mathematics, angles, scale, two),
                call(builtins, BuiltinId::Log, {
                    add(builtins, mathematics, angles, {one, u2})})});
            return divideByDu(subtract(
                builtins, mathematics, angles, std::move(main), std::move(correction)));
        }
        correction = multiply(builtins, mathematics, angles, {
            scale,
            call(builtins, BuiltinId::Sqrt, {
                subtract(builtins, mathematics, angles, one, u2)})});
        return divideByDu(definition->id == BuiltinId::Asin
            ? add(builtins, mathematics, angles, {std::move(main), std::move(correction)})
            : subtract(builtins, mathematics, angles, std::move(main), std::move(correction)));
    }

    case BuiltinId::Asinh:
        return divideByDu(subtract(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {u, expression}),
            call(builtins, BuiltinId::Sqrt, {
                add(builtins, mathematics, angles, {u2, one})})));

    case BuiltinId::Acosh:
        return divideByDu(subtract(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {u, expression}),
            multiply(builtins, mathematics, angles, {
                call(builtins, BuiltinId::Sqrt, {
                    subtract(builtins, mathematics, angles, u, one)}),
                call(builtins, BuiltinId::Sqrt, {
                    add(builtins, mathematics, angles, {u, one})})})));

    case BuiltinId::Atanh:
        return divideByDu(add(builtins, mathematics, angles, {
            multiply(builtins, mathematics, angles, {u, expression}),
            multiply(builtins, mathematics, angles, {
                rational(Rational{BigInt{1}, BigInt{2}}),
                call(builtins, BuiltinId::Log, {
                    subtract(builtins, mathematics, angles, one, u2)})})}));

    default:
        return std::nullopt;
    }
}

[[nodiscard]] Expr integrateCore(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t depth);

[[nodiscard]] std::optional<Expr> integrateMonomialTimesLog(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Multiply))
        return std::nullopt;
    const auto& factors = expression.asCall().arguments;
    if (factors.size() != 2)
        return std::nullopt;

    for (std::size_t logIndex = 0; logIndex < 2; ++logIndex) {
        const Expr& logFactor = factors[logIndex];
        const Expr& monomial = factors[1 - logIndex];
        if (!isHead(logFactor, builtins, BuiltinId::Log)
            || logFactor.asCall().arguments.size() != 1
            || !(logFactor.asCall().arguments[0].isSymbol()
                && logFactor.asCall().arguments[0].asSymbol().sameIdentity(variable)))
            continue;

        std::size_t n = 0;
        if (monomial.isSymbol() && monomial.asSymbol().sameIdentity(variable))
            n = 1;
        else if (isHead(monomial, builtins, BuiltinId::Power)
            && monomial.asCall().arguments.size() == 2
            && monomial.asCall().arguments[0].isSymbol()
            && monomial.asCall().arguments[0].asSymbol().sameIdentity(variable)) {
            const auto exponent = exactRealRational(monomial.asCall().arguments[1]);
            if (!exponent || !exponent->isInteger() || exponent->numerator().isNegative())
                continue;
            const auto value = numeric::tryToUint64(exponent->numerator());
            if (!value || *value > 256)
                continue;
            n = static_cast<std::size_t>(*value);
        }
        else {
            continue;
        }

        const Rational m{BigInt{static_cast<std::int64_t>(n + 1)}};
        Expr xPower = power(builtins, mathematics, angles,
            Expr{variable}, rational(m));
        Expr first = multiply(builtins, mathematics, angles, {
            divide(builtins, mathematics, angles, xPower, rational(m)), logFactor});
        Expr second = divide(builtins, mathematics, angles,
            xPower, rational(m * m));
        return subtract(builtins, mathematics, angles, std::move(first), std::move(second));
    }
    return std::nullopt;
}


[[nodiscard]] std::optional<Expr> tryExponentialTrigProduct(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Multiply))
        return std::nullopt;
    const auto& factors = expression.asCall().arguments;
    if (factors.size() != 2)
        return std::nullopt;

    const Expr* exponential = nullptr;
    const Expr* trigonometric = nullptr;
    BuiltinId trigId = BuiltinId::Sin;
    for (const Expr& factor : factors) {
        if (isHead(factor, builtins, BuiltinId::Exp)
            && factor.asCall().arguments.size() == 1) {
            exponential = &factor;
            continue;
        }
        if ((isHead(factor, builtins, BuiltinId::Sin)
                || isHead(factor, builtins, BuiltinId::Cos))
            && factor.asCall().arguments.size() == 1) {
            trigonometric = &factor;
            trigId = isHead(factor, builtins, BuiltinId::Sin)
                ? BuiltinId::Sin : BuiltinId::Cos;
        }
    }
    if (!exponential || !trigonometric)
        return std::nullopt;

    const Expr& exponent = exponential->asCall().arguments[0];
    const Expr& trigSource = trigonometric->asCall().arguments[0];
    Expr exponentRate = simplify(
        differentiateExpression(exponent, variable, builtins, mathematics, angles),
        builtins, mathematics, angles);
    const TrigArgument info = trigArgument(trigSource, builtins, mathematics, angles);
    Expr trigRate = multiply(builtins, mathematics, angles, {
        info.scale,
        differentiateExpression(info.argument, variable, builtins, mathematics, angles)});
    if (containsVariable(exponentRate, variable) || containsVariable(trigRate, variable))
        return std::nullopt;

    Expr denominator = add(builtins, mathematics, angles, {
        multiply(builtins, mathematics, angles, {exponentRate, exponentRate}),
        multiply(builtins, mathematics, angles, {trigRate, trigRate})});
    if (!provablyNonZero(denominator, builtins, mathematics))
        return std::nullopt;

    Expr sine = call(builtins, BuiltinId::Sin, {trigSource});
    Expr cosine = call(builtins, BuiltinId::Cos, {trigSource});
    Expr numerator = trigId == BuiltinId::Cos
        ? add(builtins, mathematics, angles, {
            multiply(builtins, mathematics, angles, {exponentRate, cosine}),
            multiply(builtins, mathematics, angles, {trigRate, sine})})
        : subtract(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {exponentRate, sine}),
            multiply(builtins, mathematics, angles, {trigRate, cosine}));
    // a=u', b=d(theta)/dx が定数なら、部分積分2回から得られる連立一次方程式のexact解。
    // a^2+b^2 != 0 は上で証明済みなので、D側の簡約能力に依存せず採用できる。
    return multiply(builtins, mathematics, angles, {
        *exponential,
        divide(builtins, mathematics, angles, std::move(numerator), std::move(denominator))});
}

[[nodiscard]] std::optional<Expr> tryIntegrationByParts(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t depth) {
    if (depth >= maximumIntegrationDepth
        || !isHead(expression, builtins, BuiltinId::Multiply))
        return std::nullopt;

    const auto& factors = expression.asCall().arguments;
    if (factors.size() < 2 || factors.size() > 8)
        return std::nullopt;

    for (std::size_t selected = 0; selected < factors.size(); ++selected) {
        if (!containsVariable(factors[selected], variable))
            continue;
        if (!factors[selected].isCall())
            continue;
        const auto* selectedDefinition = builtins.find(factors[selected].asCall().head);
        if (!selectedDefinition)
            continue;
        switch (selectedDefinition->id) {
        case BuiltinId::Exp:
        case BuiltinId::Sin:
        case BuiltinId::Cos:
        case BuiltinId::Sinh:
        case BuiltinId::Cosh:
            break;
        default:
            continue;
        }

        std::vector<Expr> polynomialFactors;
        polynomialFactors.reserve(factors.size() - 1);
        for (std::size_t i = 0; i < factors.size(); ++i)
            if (i != selected)
                polynomialFactors.push_back(factors[i]);
        Expr polynomialPart = multiply(
            builtins, mathematics, angles, std::move(polynomialFactors));

        auto polynomial = toExpressionPolynomial(
            polynomialPart, variable, builtins, mathematics, angles,
            PolynomialConversionOptions{64, 256});
        if (!polynomial || polynomial->degree() == 0)
            continue;

        Expr primitiveOfSelected = integrateCore(
            factors[selected], variable, builtins, mathematics, angles, depth + 1);
        if (isHead(primitiveOfSelected, builtins, BuiltinId::SymbolicIntegral))
            continue;

        Expr derivativePolynomial = differentiateExpression(
            polynomialPart, variable, builtins, mathematics, angles);
        Expr remainderIntegrand = multiply(
            builtins, mathematics, angles,
            {std::move(derivativePolynomial), primitiveOfSelected});
        Expr remainder = integrateCore(
            remainderIntegrand, variable, builtins, mathematics, angles, depth + 1);
        if (isHead(remainder, builtins, BuiltinId::SymbolicIntegral))
            continue;

        // ∫p g = p G - ∫p' G。pは多項式で次数が下がり、Gは既に得られた原始函数。
        // この再帰自体がexact identityなので、Simplifierのdomain-safe cancellation能力には依存しない。
        Expr candidate = subtract(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles,
                {polynomialPart, primitiveOfSelected}),
            std::move(remainder));
        return simplify(expandExpression(
            candidate, builtins, mathematics, angles, AlgebraTransformOptions{512}),
            builtins, mathematics, angles);
    }
    return std::nullopt;
}

[[nodiscard]] Expr integrateCore(
    const Expr& original,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t depth) {
    if (depth > maximumIntegrationDepth)
        return unresolved(original, variable, builtins);

    const Expr expression = simplify(original, builtins, mathematics, angles);
    if (!containsVariable(expression, variable))
        return multiply(builtins, mathematics, angles, {expression, Expr{variable}});

    if (expression.isSymbol() && expression.asSymbol().sameIdentity(variable))
        return divide(builtins, mathematics, angles,
            power(builtins, mathematics, angles, Expr{variable}, integer(2)), integer(2));

    if (isHead(expression, builtins, BuiltinId::Power)) {
        if (auto result = integratePowerRule(
                expression, variable, builtins, mathematics, angles))
            return *result;
        if (const auto rewritten = rewriteSquareIdentity(
                expression, builtins, mathematics, angles)) {
            Expr result = integrateCore(
                *rewritten, variable, builtins, mathematics, angles, depth + 1);
            if (!isHead(result, builtins, BuiltinId::SymbolicIntegral))
                return result;
        }
    }

    if (auto polynomial = integratePolynomial(
            expression, variable, builtins, mathematics, angles))
        return *polynomial;

    if (auto rationalQuadratic = tryRationalLowDegree(
            expression, variable, builtins, mathematics, angles))
        return *rationalQuadratic;
    if (auto rationalFactored = tryRationalFactored(
            expression, variable, builtins, mathematics, angles))
        return *rationalFactored;
    if (auto quadratic = tryQuadraticReciprocal(
            expression, variable, builtins, mathematics, angles))
        return *quadratic;
    if (auto quadraticRoot = tryQuadraticInverseSqrt(
            expression, variable, builtins, mathematics, angles))
        return *quadraticRoot;
    if (auto nestedRoot = tryNestedQuadraticSquareRoot(
            expression, variable, builtins, mathematics, angles))
        return *nestedRoot;

    if (expression.isCall()) {
        const auto* definition = builtins.find(expression.asCall().head);
        const auto& a = expression.asCall().arguments;
        if (definition) {
            switch (definition->id) {
            case BuiltinId::Log:
                if (a.size() == 2 && !containsVariable(a[0], variable)) {
                    Expr unaryLog = call(builtins, BuiltinId::Log, {a[1]});
                    Expr primitive = integrateCore(
                        unaryLog, variable, builtins, mathematics, angles, depth + 1);
                    if (!isHead(primitive, builtins, BuiltinId::SymbolicIntegral)) {
                        return divide(builtins, mathematics, angles,
                            std::move(primitive), call(builtins, BuiltinId::Log, {a[0]}));
                    }
                }
                break;
            case BuiltinId::Add: {
                // 線形性は未評価項が残っても厳密。既知項まで巻き戻さず、
                // integrate[unknown,x] をその項だけに保持して部分結果を返す。
                std::vector<Expr> terms;
                terms.reserve(a.size());
                for (const Expr& term : a)
                    terms.push_back(integrateCore(
                        term, variable, builtins, mathematics, angles, depth + 1));
                return add(builtins, mathematics, angles, std::move(terms));
            }
            case BuiltinId::Subtract:
                if (a.size() == 2) {
                    Expr lhs = integrateCore(a[0], variable, builtins, mathematics, angles, depth + 1);
                    Expr rhs = integrateCore(a[1], variable, builtins, mathematics, angles, depth + 1);
                    return subtract(builtins, mathematics, angles, std::move(lhs), std::move(rhs));
                }
                break;
            case BuiltinId::Negate:
                if (a.size() == 1)
                    return negate(builtins, mathematics, angles, integrateCore(
                        a[0], variable, builtins, mathematics, angles, depth + 1));
                break;
            case BuiltinId::Multiply: {
                FactorSplit split = splitConstantFactor(
                    expression, variable, builtins, mathematics, angles);
                if (!isOne(split.constant) && !isOne(split.dependent)) {
                    Expr integrated = integrateCore(
                        split.dependent, variable, builtins, mathematics, angles, depth + 1);
                    return multiply(
                        builtins, mathematics, angles,
                        {split.constant, std::move(integrated)});
                }
                if (auto logProduct = integrateMonomialTimesLog(
                        expression, variable, builtins, mathematics, angles))
                    return *logProduct;
                if (auto expTrig = tryExponentialTrigProduct(
                        expression, variable, builtins, mathematics, angles))
                    return *expTrig;
                if (auto byParts = tryIntegrationByParts(
                        expression, variable, builtins, mathematics, angles, depth))
                    return *byParts;
                break;
            }
            case BuiltinId::Divide:
                if (a.size() == 2 && !containsVariable(a[0], variable)) {
                    Expr denominatorDerivative = differentiateExpression(
                        a[1], variable, builtins, mathematics, angles);
                    if (!containsVariable(denominatorDerivative, variable)
                        && provablyNonZero(denominatorDerivative, builtins, mathematics)) {
                        return multiply(builtins, mathematics, angles, {
                            a[0],
                            divide(builtins, mathematics, angles,
                                call(builtins, BuiltinId::Log, {a[1]}),
                                std::move(denominatorDerivative))});
                    }
                }
                if (a.size() == 2 && !containsVariable(a[1], variable)) {
                    Expr numerator = integrateCore(
                        a[0], variable, builtins, mathematics, angles, depth + 1);
                    return divide(
                        builtins, mathematics, angles, std::move(numerator), a[1]);
                }
                break;
            case BuiltinId::Power:
                if (auto result = integratePowerRule(
                        expression, variable, builtins, mathematics, angles))
                    return *result;
                break;
            default:
                break;
            }
        }
    }

    if (auto elementary = integrateElementaryUnary(
            expression, variable, builtins, mathematics, angles))
        return *elementary;

    if (auto standard = integrateStandardUnary(
            expression, variable, builtins, mathematics, angles))
        return *standard;

    if (auto reverse = tryReverseChainRule(
            expression, variable, builtins, mathematics, angles))
        return *reverse;

    if (auto byParts = tryIntegrationByParts(
            expression, variable, builtins, mathematics, angles, depth))
        return *byParts;

    return unresolved(original, variable, builtins);
}

[[nodiscard]] bool isInfinityExpr(
    const Expr& expression,
    const expression::Symbol& infinity) {
    return expression.isSymbol() && expression.asSymbol().sameIdentity(infinity);
}

[[nodiscard]] bool isNegativeInfinityExpr(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    return isHead(expression, builtins, BuiltinId::Negate)
        && expression.asCall().arguments.size() == 1
        && isInfinityExpr(expression.asCall().arguments[0], infinity);
}

[[nodiscard]] bool boundsAreProvablyReal(
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    return knowledge.facts(lower).isProvablyReal()
        && knowledge.facts(upper).isProvablyReal();
}

[[nodiscard]] bool domainConditionsProven(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const auto conditions = mathematics::expressionDomainConditions(expression, builtins, mathematics);
    if (!conditions)
        return false;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    for (const mathematics::Predicate& predicate : conditions->predicates()) {
        if (knowledge.prove(predicate) != mathematics::TruthValue::True)
            return false;
    }
    return true;
}

[[nodiscard]] bool integrandHasNoDomainConditions(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    const auto conditions = mathematics::expressionDomainConditions(expression, builtins, mathematics);
    return conditions && conditions->empty();
}

[[nodiscard]] bool certifyOnNumericRealInterval(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    try {
        constexpr std::size_t bits = 192;
        CertifiedEvaluator evaluator{builtins, mathematics, angles};
        const auto lowerValue = evaluator.enclose(lower, bits);
        const auto upperValue = evaluator.enclose(upper, bits);
        if (!lowerValue || !upperValue || !lowerValue->isReal() || !upperValue->isReal())
            return false;

        const RealInterval interval = approximation::hull(
            lowerValue->asReal(), upperValue->asReal());
        const CertifiedBinding binding{variable, CertifiedValue{interval}};
        return evaluator.enclose(
            expression, bits, std::span<const CertifiedBinding>{&binding, 1}).has_value();
    }
    catch (...) {
        return false;
    }
}

[[nodiscard]] bool safeForDefiniteIntegral(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!boundsAreProvablyReal(lower, upper, builtins, mathematics, assumptions))
        return false;
    if (integrandHasNoDomainConditions(expression, builtins, mathematics))
        return true;
    if (domainConditionsProven(expression, builtins, mathematics, assumptions))
        return true;
    return certifyOnNumericRealInterval(
        expression, variable, lower, upper, builtins, mathematics, angles);
}

[[nodiscard]] bool rationalDenominatorHasNoRealPoleOnImproperPath(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return false;
    const auto denominator = toRationalPolynomial(
        expression.asCall().arguments[1], variable, builtins, {16, 128});
    if (!denominator || denominator->isZero())
        return false;
    if (denominator->degree() == 0)
        return true;

    const auto lowerValue = exactRealRational(lower);
    const auto upperValue = exactRealRational(upper);
    const bool lowerInfinite = isNegativeInfinityExpr(lower, builtins, infinity);
    const bool upperInfinite = isInfinityExpr(upper, infinity);
    if ((!lowerValue && !lowerInfinite) || (!upperValue && !upperInfinite))
        return false;

    auto inOpenPath = [&](const Rational& root) {
        if (!lowerInfinite && !(*lowerValue < root))
            return false;
        if (!upperInfinite && !(root < *upperValue))
            return false;
        return true;
    };

    RationalPolynomial remaining = *denominator;
    for (std::size_t count = 0; count < 16 && remaining.degree() > 0; ++count) {
        const RationalRootSearchResult search = findRationalRoot(
            remaining, RationalRootSearchOptions{200000, 10000});
        if (search.root) {
            if (inOpenPath(*search.root))
                return false;
            const auto quotient = divideByLinearFactor(remaining, *search.root);
            if (!quotient)
                return false;
            remaining = *quotient;
            continue;
        }
        if (remaining.degree() == 2) {
            const Rational a = remaining.coefficient(2);
            const Rational b = remaining.coefficient(1);
            const Rational c = remaining.coefficient(0);
            const Rational discriminant = b * b - Rational{BigInt{4}} * a * c;
            // 負の判別式なら実根なし。0は重根でRationalなので上の探索で見つかる。
            if (discriminant.numerator().isNegative())
                return true;
            return false;
        }
        return search.complete && remaining.degree() == 0;
    }
    return remaining.degree() == 0;
}

[[nodiscard]] bool inverseSqrtEndpointClass(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    if (isInfinityExpr(upper, infinity) || isNegativeInfinityExpr(lower, builtins, infinity))
        return false;
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return false;
    const Expr& numerator = expression.asCall().arguments[0];
    const Expr& denominator = expression.asCall().arguments[1];
    if (!numerator.isNumber() || !isHead(denominator, builtins, BuiltinId::Sqrt)
        || denominator.asCall().arguments.size() != 1)
        return false;
    const auto affine = toRationalPolynomial(
        denominator.asCall().arguments[0], variable, builtins, {1, 4});
    const auto lo = exactRealRational(lower);
    const auto hi = exactRealRational(upper);
    if (!affine || affine->degree() != 1 || !lo || !hi)
        return false;
    const Rational slope = affine->coefficient(1);
    const Rational atLo = evaluatePolynomial(*affine, *lo);
    const Rational atHi = evaluatePolynomial(*affine, *hi);
    if (slope.numerator().isNegative())
        return atHi.isZero() && !atLo.numerator().isNegative();
    return atLo.isZero() && !atHi.numerator().isNegative();
}

[[nodiscard]] bool logarithmEndpointClass(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Log)
        || expression.asCall().arguments.size() != 1)
        return false;
    const auto affine = toRationalPolynomial(
        expression.asCall().arguments[0], variable, builtins, {1, 4});
    const auto lo = exactRealRational(lower);
    const auto hi = exactRealRational(upper);
    if (!affine || affine->degree() != 1 || !lo || !hi)
        return false;
    const Rational atLo = evaluatePolynomial(*affine, *lo);
    const Rational atHi = evaluatePolynomial(*affine, *hi);
    const Rational slope = affine->coefficient(1);
    if (slope.numerator().isNegative())
        return atHi.isZero() && !atLo.numerator().isNegative();
    return atLo.isZero() && !atHi.numerator().isNegative();
}

[[nodiscard]] bool safeForImproperIntegral(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const expression::Symbol& infinity) {
    if (integrandHasNoDomainConditions(expression, builtins, mathematics))
        return true;
    if (rationalDenominatorHasNoRealPoleOnImproperPath(
            expression, variable, lower, upper, builtins, infinity))
        return true;
    if (inverseSqrtEndpointClass(
            expression, variable, lower, upper, builtins, infinity))
        return true;
    return logarithmEndpointClass(expression, variable, lower, upper, builtins);
}

[[nodiscard]] Expr endpointPrimitiveValue(
    const Expr& primitive,
    const expression::Symbol& variable,
    const Expr& endpoint,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    if (isInfinityExpr(endpoint, infinity) || isNegativeInfinityExpr(endpoint, builtins, infinity))
        return limitExpression(
            primitive, variable, endpoint, LimitDirection::TwoSided,
            builtins, mathematics, angles, infinity, assumptions);

    Expr raw = substituteSymbol(primitive, variable, endpoint);
    const auto rawConditions = mathematics::expressionDomainConditions(raw, builtins, mathematics);
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    bool falseCondition = false;
    if (rawConditions)
        for (const auto& predicate : rawConditions->predicates())
            falseCondition |= knowledge.prove(predicate) == mathematics::TruthValue::False;
    if (!falseCondition) {
        Expr direct = simplification::Simplifier{}.simplify(
            std::move(raw),
            simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
        return direct;
    }

    return limitExpression(
        primitive, variable, endpoint, direction,
        builtins, mathematics, angles, infinity, assumptions);
}

[[nodiscard]] Expr subtractEndpointValues(
    Expr high,
    Expr low,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions,
    const Expr& unresolvedExpression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper) {
    const bool highPosInf = isInfinityExpr(high, infinity);
    const bool highNegInf = isNegativeInfinityExpr(high, builtins, infinity);
    const bool lowPosInf = isInfinityExpr(low, infinity);
    const bool lowNegInf = isNegativeInfinityExpr(low, builtins, infinity);
    if ((highPosInf && lowPosInf) || (highNegInf && lowNegInf)
        || (highPosInf && lowNegInf) || (highNegInf && lowPosInf)) {
        // 両端が無限大になる場合、符号だけではCauchy cancellationを定義しない。
        return call(builtins, BuiltinId::SymbolicIntegral, {
            unresolvedExpression, Expr::array({3}, {Expr{variable}, lower, upper})});
    }
    if (highPosInf || lowNegInf)
        return Expr{infinity};
    if (highNegInf)
        return call(builtins, BuiltinId::Negate, {Expr{infinity}});
    if (lowPosInf)
        return call(builtins, BuiltinId::Negate, {Expr{infinity}});

    Expr difference = subtract(
        builtins, mathematics, angles, std::move(high), std::move(low));
    difference = expandExpression(
        difference, builtins, mathematics, angles, AlgebraTransformOptions{1024});
    return simplification::fullSimplify(
        std::move(difference),
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions},
        simplification::FullSimplificationOptions{48});
}

} // namespace

Expr integrateExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const Expr prepared = simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
    return integrateCore(prepared, variable, builtins, mathematics, angles, 0);
}

Expr integrateExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions) {
    const Expr prepared = simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
    Expr primitive = integrateCore(prepared, variable, builtins, mathematics, angles, 0);
    if (isHead(primitive, builtins, BuiltinId::SymbolicIntegral)) {
        return call(builtins, BuiltinId::SymbolicIntegral, {
            expression,
            Expr::array({3}, {Expr{variable}, lower, upper})});
    }

    const bool improper = isInfinityExpr(upper, infinitySymbol)
        || isNegativeInfinityExpr(lower, builtins, infinitySymbol)
        || !safeForDefiniteIntegral(
            prepared, variable, lower, upper,
            builtins, mathematics, angles, assumptions);

    if (improper) {
        if (!safeForImproperIntegral(
                prepared, variable, lower, upper,
                builtins, mathematics, infinitySymbol)) {
            return call(builtins, BuiltinId::SymbolicIntegral, {
                expression,
                Expr::array({3}, {Expr{variable}, lower, upper})});
        }
        Expr high = endpointPrimitiveValue(
            primitive, variable, upper, LimitDirection::Left,
            builtins, mathematics, angles, infinitySymbol, assumptions);
        Expr low = endpointPrimitiveValue(
            primitive, variable, lower, LimitDirection::Right,
            builtins, mathematics, angles, infinitySymbol, assumptions);
        if (isHead(high, builtins, BuiltinId::Limit)
            || isHead(low, builtins, BuiltinId::Limit)) {
            return call(builtins, BuiltinId::SymbolicIntegral, {
                expression,
                Expr::array({3}, {Expr{variable}, lower, upper})});
        }
        return subtractEndpointValues(
            std::move(high), std::move(low), builtins, mathematics, angles,
            infinitySymbol, assumptions, expression, variable, lower, upper);
    }

    Expr high = simplification::Simplifier{}.simplify(
        substituteSymbol(primitive, variable, upper),
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
    Expr low = simplification::Simplifier{}.simplify(
        substituteSymbol(primitive, variable, lower),
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
    return subtractEndpointValues(
        std::move(high), std::move(low), builtins, mathematics, angles,
        infinitySymbol, assumptions, expression, variable, lower, upper);
}

} // namespace mmcal::symbolic
