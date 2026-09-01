// 記号積分integrate
#include "integration.hpp"

#include "approximation/certified_evaluator.hpp"
#include "approximation/real_interval.hpp"
#include "mathematics/definedness.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/trigonometric_polynomial.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "simplification/full_simplifier.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/algebraic_expression.hpp"
#include "symbolic/limit.hpp"
#include "symbolic/algebra_transforms.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"
#include "symbolic/series.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
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
    return builtins.isCallTo(expression, id);
}

[[nodiscard]] bool containsHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    if (isHead(expression, builtins, id))
        return true;
    if (!expression.isCall())
        return false;
    for (const Expr& argument : expression.asCall().arguments)
        if (containsHead(argument, builtins, id))
            return true;
    return false;
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

[[nodiscard]] std::optional<Expr> integrateGenericSymbolicPowerRule(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Power)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const auto& arguments = expression.asCall().arguments;
    if (containsVariable(arguments[1], variable)
        || exactRealRational(arguments[1]))
        return std::nullopt;

    Expr baseDerivative = simplify(
        differentiateExpression(arguments[0], variable, builtins, mathematics, angles),
        builtins, mathematics, angles);
    if (containsVariable(baseDerivative, variable)
        || !provablyNonZero(baseDerivative, builtins, mathematics))
        return std::nullopt;

    // Generic parameter rule。明示的なexponent=-1は上のexact ruleでlogへ送る。
    // 記号parameterが後に-1へ特殊化される点はこのgeneric代表式の例外であり，
    // integrate[x^-1,x]自体は別規則を持つ。
    Expr nextExponent = simplify(
        add(builtins, mathematics, angles, {arguments[1], integer(1)}),
        builtins, mathematics, angles);
    Expr denominator = multiply(
        builtins, mathematics, angles, {std::move(baseDerivative), nextExponent});
    return divide(
        builtins, mathematics, angles,
        power(builtins, mathematics, angles, arguments[0], std::move(nextExponent)),
        std::move(denominator));
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
    candidates.reserve(48);

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
    candidates.push_back(unary(BuiltinId::ExponentialIntegralEi));
    candidates.push_back(unary(BuiltinId::SineIntegralSi));
    candidates.push_back(unary(BuiltinId::CosineIntegralCi));
    candidates.push_back(unary(BuiltinId::LogarithmicIntegralLi));
    // D[Gamma]=Gamma Digamma, D[LogGamma]=Digamma, D[Digamma]=Trigamma。
    // reverse-chain候補へ共有して積分側でも既存の微分知識を再利用する。
    candidates.push_back(unary(BuiltinId::Gamma));
    candidates.push_back(unary(BuiltinId::LogGamma));
    candidates.push_back(unary(BuiltinId::Digamma));
    candidates.push_back(call(builtins, BuiltinId::Polylog, {integer(2), u}));
    // Li_2(-u) は log(1+u)/u 系の逆chainを一括で拾う。
    // 個別に log[1+x]/x, log[1+x^2]/x を表登録せず、Dで比例係数を決める。
    candidates.push_back(negate(builtins, mathematics, angles,
        call(builtins, BuiltinId::Polylog, {
            integer(2), negate(builtins, mathematics, angles, u)})));
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

    // 対数微分の基本Knowledge。候補をDで厳密検証してから採用するため、
    // branchを無視したglobal rewriteにはせず、局所primitiveとしてだけ使う。
    const Expr logU = unary(BuiltinId::Log);
    candidates.push_back(power(builtins, mathematics, angles, logU, two));
    candidates.push_back(call(builtins, BuiltinId::Log, {logU}));

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

    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::IntegrationCandidate,
        candidates.size());
    return candidates;
}

[[nodiscard]] std::vector<Expr> dependentSubexpressions(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Expr> result;
    std::vector<Expr> pending{expression};
    while (!pending.empty() && result.size() < maximumSubstitutionCandidates) {
        Expr current = std::move(pending.back());
        pending.pop_back();
        if (!containsVariable(current, variable))
            continue;

        const auto duplicate = std::find(result.begin(), result.end(), current);
        if (duplicate == result.end()) {
            evaluation::consumeEvaluationBudget(
                evaluation::EvaluationResource::IntegrationCandidate);
            result.push_back(current);
        }

        if (current.isCall()) {
            for (const Expr& argument : current.asCall().arguments)
                pending.push_back(argument);
        }
    }

    // x^4 の内部ASTには x^2 が存在しないため、u=x^2 で一発の
    // x/(1+x^4), x/sqrt[1-x^4] を従来のsubexpression探索だけでは見落としていた。
    // 変数の正整数冪が現れた場合、その真の約数冪も少数だけsubstitution候補へ加える。
    const std::vector<Expr> snapshot = result;
    for (const Expr& candidate : snapshot) {
        if (result.size() >= maximumSubstitutionCandidates
            || !isHead(candidate, builtins, BuiltinId::Power)
            || candidate.asCall().arguments.size() != 2)
            continue;

        const auto& arguments = candidate.asCall().arguments;
        if (!arguments[0].isSymbol()
            || !arguments[0].asSymbol().sameIdentity(variable))
            continue;
        const auto exponent = exactRealRational(arguments[1]);
        if (!exponent || !exponent->isInteger() || exponent->numerator().isNegative())
            continue;
        const auto magnitude = numeric::tryToUint64(exponent->numerator());
        if (!magnitude || *magnitude < 4 || *magnitude > 4096)
            continue;

        for (std::uint64_t divisor = 2;
             divisor * divisor <= *magnitude && result.size() < maximumSubstitutionCandidates;
             ++divisor) {
            if ((*magnitude % divisor) != 0)
                continue;
            const auto appendPower = [&](std::uint64_t powerValue) {
                Expr generated = Expr::call(builtins.symbol(BuiltinId::Power), {
                    Expr{variable}, integer(static_cast<std::int64_t>(powerValue))});
                if (std::find(result.begin(), result.end(), generated) == result.end())
                    result.push_back(std::move(generated));
            };
            appendPower(divisor);
            if (divisor * divisor != *magnitude)
                appendPower(*magnitude / divisor);
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
    for (const Expr& u : dependentSubexpressions(integrand, variable, builtins)) {
        // まずsubexpression自身を原始函数候補としてDで検証する。
        // Gamma*DigammaやBetaのlogarithmic derivativeなど，既にDが知る函数族を
        // integrate側へ個別に二重登録せず再利用できる。
        if (auto result = verifiedScaledPrimitive(
                integrand, u, variable, builtins, mathematics, angles))
            return result;

        for (Expr candidate : primitiveTemplates(u, builtins, mathematics, angles)) {
            if (auto result = verifiedScaledPrimitive(
                    integrand, std::move(candidate), variable,
                    builtins, mathematics, angles))
                return result;
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> trySpecialFunctionShiftPrimitive(
    const Expr& integrand,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    // 微分側に既にあるparameter-shift公式を逆向きに使う。
    // candidateは必ずDで元のintegrandとの比例関係をexact検証してから採用する。
    for (const Expr& node : dependentSubexpressions(integrand, variable, builtins)) {
        if (!node.isCall())
            continue;
        const auto* definition = builtins.find(node.asCall().head);
        if (!definition)
            continue;
        const auto& a = node.asCall().arguments;
        std::optional<Expr> candidate;

        switch (definition->id) {
        case BuiltinId::Polylog:
            if (a.size() == 2 && !containsVariable(a[0], variable)) {
                Expr nextOrder = simplify(
                    add(builtins, mathematics, angles, {a[0], integer(1)}),
                    builtins, mathematics, angles);
                Expr primitive = call(builtins, BuiltinId::Polylog, {
                    std::move(nextOrder), a[1]});

                // D側はz=0のremovable singularityをifで埋めるが，元integrandが
                // Li_s(z) z'/z と明示している場合はそのpunctured-domain公式を直接照合する。
                Expr dz = simplify(
                    differentiateExpression(a[1], variable, builtins, mathematics, angles),
                    builtins, mathematics, angles);
                Expr kernel = multiply(builtins, mathematics, angles, {
                    node, divide(builtins, mathematics, angles, std::move(dz), a[1])});
                if (const auto ratio = proportionalFactor(
                        integrand, kernel, variable, builtins, mathematics, angles))
                    return multiply(builtins, mathematics, angles, {*ratio, std::move(primitive)});
                candidate = std::move(primitive);
            }
            break;

        case BuiltinId::Hypergeometric1F1:
            if (a.size() == 3
                && !containsVariable(a[0], variable)
                && !containsVariable(a[1], variable)) {
                Expr previousA = simplify(
                    subtract(builtins, mathematics, angles, a[0], integer(1)),
                    builtins, mathematics, angles);
                Expr previousB = simplify(
                    subtract(builtins, mathematics, angles, a[1], integer(1)),
                    builtins, mathematics, angles);
                // M(a-1,b-1,z)がparameter poleへ落ちる場合は使わない。
                if (provablyNonZero(previousA, builtins, mathematics)
                    && provablyNonZero(previousB, builtins, mathematics)) {
                    candidate = call(builtins, BuiltinId::Hypergeometric1F1, {
                        std::move(previousA), std::move(previousB), a[2]});
                }
            }
            break;

        case BuiltinId::Hypergeometric2F1:
            if (a.size() == 4
                && !containsVariable(a[0], variable)
                && !containsVariable(a[1], variable)
                && !containsVariable(a[2], variable)) {
                Expr previousA = simplify(
                    subtract(builtins, mathematics, angles, a[0], integer(1)),
                    builtins, mathematics, angles);
                Expr previousB = simplify(
                    subtract(builtins, mathematics, angles, a[1], integer(1)),
                    builtins, mathematics, angles);
                Expr previousC = simplify(
                    subtract(builtins, mathematics, angles, a[2], integer(1)),
                    builtins, mathematics, angles);
                if (provablyNonZero(previousA, builtins, mathematics)
                    && provablyNonZero(previousB, builtins, mathematics)
                    && provablyNonZero(previousC, builtins, mathematics)) {
                    candidate = call(builtins, BuiltinId::Hypergeometric2F1, {
                        std::move(previousA), std::move(previousB),
                        std::move(previousC), a[3]});
                }
            }
            break;

        case BuiltinId::IncompleteBeta:
            if (a.size() == 3
                && !containsVariable(a[0], variable)
                && !containsVariable(a[1], variable)) {
                Expr parameterSum = simplify(
                    add(builtins, mathematics, angles, {a[0], a[1]}),
                    builtins, mathematics, angles);
                if (provablyNonZero(parameterSum, builtins, mathematics)) {
                    Expr shiftedA = simplify(
                        add(builtins, mathematics, angles, {a[0], integer(1)}),
                        builtins, mathematics, angles);
                    Expr correction = multiply(builtins, mathematics, angles, {
                        divide(builtins, mathematics, angles, a[0], parameterSum),
                        call(builtins, BuiltinId::IncompleteBeta, {
                            std::move(shiftedA), a[1], a[2]})});
                    candidate = subtract(builtins, mathematics, angles,
                        multiply(builtins, mathematics, angles, {a[2], node}),
                        std::move(correction));
                }
            }
            break;

        case BuiltinId::LambertW:
            if (a.size() == 1 || (a.size() == 2 && !containsVariable(a[0], variable))) {
                // d(W+W^2/2) = W(z) z'/z。branch indexを含むW_kでも同じ。
                candidate = add(builtins, mathematics, angles, {
                    node,
                    divide(builtins, mathematics, angles,
                        power(builtins, mathematics, angles, node, integer(2)), integer(2))});
            }
            break;

        default:
            break;
        }

        if (candidate) {
            if (auto result = verifiedScaledPrimitive(
                    integrand, std::move(*candidate), variable,
                    builtins, mathematics, angles))
                return result;
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> integrateIncompleteBetaAffine(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::IncompleteBeta)
        || expression.asCall().arguments.size() != 3)
        return std::nullopt;
    const auto& a = expression.asCall().arguments;
    if (containsVariable(a[0], variable) || containsVariable(a[1], variable))
        return std::nullopt;

    Expr du = simplify(
        differentiateExpression(a[2], variable, builtins, mathematics, angles),
        builtins, mathematics, angles);
    if (containsVariable(du, variable) || !provablyNonZero(du, builtins, mathematics))
        return std::nullopt;

    Expr parameterSum = simplify(
        add(builtins, mathematics, angles, {a[0], a[1]}),
        builtins, mathematics, angles);
    if (!provablyNonZero(parameterSum, builtins, mathematics))
        return std::nullopt;

    // ∫ I_u(a,b) du = u I_u(a,b) - a/(a+b) I_u(a+1,b)。
    // 現在のibetaはregularized incomplete betaなのでこの係数を使う。
    Expr shiftedA = simplify(
        add(builtins, mathematics, angles, {a[0], integer(1)}),
        builtins, mathematics, angles);
    Expr correction = multiply(builtins, mathematics, angles, {
        divide(builtins, mathematics, angles, a[0], parameterSum),
        call(builtins, BuiltinId::IncompleteBeta, {
            std::move(shiftedA), a[1], a[2]})});
    Expr primitive = subtract(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles, {a[2], expression}),
        std::move(correction));
    return divide(builtins, mathematics, angles, std::move(primitive), std::move(du));
}

[[nodiscard]] std::optional<Expr> integrateLambertWLogDerivative(
    const Expr& integrand,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    for (const Expr& node : dependentSubexpressions(integrand, variable, builtins)) {
        if (!isHead(node, builtins, BuiltinId::LambertW))
            continue;
        const auto& a = node.asCall().arguments;
        if ((a.size() != 1 && a.size() != 2)
            || (a.size() == 2 && containsVariable(a[0], variable)))
            continue;
        const Expr& z = a.back();
        Expr dz = simplify(
            differentiateExpression(z, variable, builtins, mathematics, angles),
            builtins, mathematics, angles);
        if (containsVariable(dz, variable) && isZero(dz))
            continue;
        Expr kernel = multiply(builtins, mathematics, angles, {
            node, divide(builtins, mathematics, angles, dz, z)});
        const auto ratio = proportionalFactor(
            integrand, kernel, variable, builtins, mathematics, angles);
        if (!ratio)
            continue;

        Expr primitive = add(builtins, mathematics, angles, {
            node,
            divide(builtins, mathematics, angles,
                power(builtins, mathematics, angles, node, integer(2)), integer(2))});
        return multiply(builtins, mathematics, angles, {*ratio, std::move(primitive)});
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

    if (*exponent == Rational{BigInt{-1}}) {
        Expr logArgument = a[0];
        if (const auto affine = toRationalPolynomial(
                a[0], variable, builtins, PolynomialConversionOptions{2, 16});
            affine && affine->degree() == 1)
            logArgument = polynomialToExpandedExpr(*affine, variable, builtins);
        return divide(builtins, mathematics, angles,
            call(builtins, BuiltinId::Log, {std::move(logArgument)}), baseDerivative);
    }

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

[[nodiscard]] RationalPolynomial addPolynomials(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    const std::size_t count = std::max(
        lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<Rational> coefficients(count, Rational{BigInt{0}});
    for (std::size_t i = 0; i < count; ++i)
        coefficients[i] = lhs.coefficient(i) + rhs.coefficient(i);
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] RationalPolynomial negatePolynomial(
    const RationalPolynomial& polynomial) {
    std::vector<Rational> coefficients = polynomial.coefficients();
    for (Rational& coefficient : coefficients)
        coefficient = -coefficient;
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] RationalPolynomial scalePolynomial(
    const RationalPolynomial& polynomial,
    const Rational& scale) {
    std::vector<Rational> coefficients = polynomial.coefficients();
    for (Rational& coefficient : coefficients)
        coefficient *= scale;
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] RationalPolynomial polynomialGcdMonic(
    RationalPolynomial lhs,
    RationalPolynomial rhs) {
    while (!rhs.isZero()) {
        PolynomialDivision division = dividePolynomials(lhs, rhs);
        lhs = std::move(rhs);
        rhs = std::move(division.remainder);
    }
    if (lhs.isZero())
        return RationalPolynomial{{Rational{BigInt{1}}}};
    const Rational leading = lhs.coefficient(lhs.degree());
    return scalePolynomial(lhs, Rational{BigInt{1}} / leading);
}


[[nodiscard]] RationalPolynomial differentiatePolynomial(
    const RationalPolynomial& polynomial) {
    if (polynomial.degree() == 0)
        return RationalPolynomial{};
    std::vector<Rational> coefficients(polynomial.degree(), Rational{BigInt{0}});
    for (std::size_t exponent = 1; exponent <= polynomial.degree(); ++exponent)
        coefficients[exponent - 1] = polynomial.coefficient(exponent)
            * Rational{BigInt::fromUnsigned(exponent)};
    return RationalPolynomial{std::move(coefficients)};
}

struct PolynomialExtendedGcd final {
    RationalPolynomial gcd;
    RationalPolynomial lhsCoefficient;
    RationalPolynomial rhsCoefficient;
};

[[nodiscard]] PolynomialExtendedGcd polynomialExtendedGcd(
    RationalPolynomial lhs,
    RationalPolynomial rhs) {
    RationalPolynomial oldS{{Rational{BigInt{1}}}};
    RationalPolynomial s;
    RationalPolynomial oldT;
    RationalPolynomial t{{Rational{BigInt{1}}}};

    while (!rhs.isZero()) {
        PolynomialDivision division = dividePolynomials(lhs, rhs);
        RationalPolynomial nextS = addPolynomials(
            oldS, negatePolynomial(multiplyPolynomials(division.quotient, s)));
        RationalPolynomial nextT = addPolynomials(
            oldT, negatePolynomial(multiplyPolynomials(division.quotient, t)));
        lhs = std::move(rhs);
        rhs = std::move(division.remainder);
        oldS = std::move(s);
        s = std::move(nextS);
        oldT = std::move(t);
        t = std::move(nextT);
    }

    if (lhs.isZero())
        return PolynomialExtendedGcd{
            RationalPolynomial{}, RationalPolynomial{}, RationalPolynomial{}};
    const Rational inverseLeading = Rational{BigInt{1}} / lhs.coefficient(lhs.degree());
    return PolynomialExtendedGcd{
        scalePolynomial(lhs, inverseLeading),
        scalePolynomial(oldS, inverseLeading),
        scalePolynomial(oldT, inverseLeading)};
}

[[nodiscard]] std::optional<RationalPolynomial> polynomialInverseModulo(
    const RationalPolynomial& value,
    const RationalPolynomial& modulus) {
    PolynomialExtendedGcd extended = polynomialExtendedGcd(value, modulus);
    if (extended.gcd.degree() != 0 || extended.gcd.coefficient(0).isZero())
        return std::nullopt;
    PolynomialDivision reduced = dividePolynomials(extended.lhsCoefficient, modulus);
    return reduced.remainder;
}

[[nodiscard]] std::optional<AlgebraicNumber> evaluatePolynomialAtAlgebraic(
    const RationalPolynomial& polynomial,
    const AlgebraicNumber& point) {
    auto value = AlgebraicNumber::fromRational(
        polynomial.coefficient(polynomial.degree()));
    if (!value)
        return std::nullopt;
    for (std::size_t exponent = polynomial.degree(); exponent-- > 0;) {
        auto multiplied = AlgebraicNumber::combine(
            *value, point, AlgebraicBinaryOperation::Multiply);
        if (!multiplied)
            return std::nullopt;
        auto coefficient = AlgebraicNumber::fromRational(polynomial.coefficient(exponent));
        if (!coefficient)
            return std::nullopt;
        auto added = AlgebraicNumber::combine(
            *multiplied, *coefficient, AlgebraicBinaryOperation::Add);
        if (!added)
            return std::nullopt;
        value = std::move(added);
    }
    return value;
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


[[nodiscard]] bool isNonzeroConstantPolynomial(const RationalPolynomial& polynomial) {
    return polynomial.degree() == 0 && !polynomial.coefficient(0).isZero();
}

[[nodiscard]] std::optional<std::vector<RationalFactor>> factorIntoSquareFreePowers(
    RationalPolynomial polynomial) {
    if (polynomial.degree() == 0)
        return std::nullopt;

    // Yunのsquare-free decompositionをQ[x]上で行う。irreducible factorizationは不要で，
    // f = product f_i^i の各f_iだけをgcd(f,f')からexactに回収する。
    const Rational leading = polynomial.coefficient(polynomial.degree());
    polynomial = scalePolynomial(polynomial, Rational{BigInt{1}} / leading);
    RationalPolynomial repeated = polynomialGcdMonic(
        polynomial, differentiatePolynomial(polynomial));
    PolynomialDivision initial = dividePolynomials(polynomial, repeated);
    if (!initial.remainder.isZero())
        return std::nullopt;
    RationalPolynomial remaining = std::move(initial.quotient);

    std::vector<RationalFactor> factors;
    for (std::size_t multiplicity = 1;
         !isNonzeroConstantPolynomial(remaining);
         ++multiplicity) {
        if (multiplicity > polynomial.degree())
            return std::nullopt;
        RationalPolynomial shared = polynomialGcdMonic(remaining, repeated);
        PolynomialDivision distinct = dividePolynomials(remaining, shared);
        if (!distinct.remainder.isZero())
            return std::nullopt;
        if (!isNonzeroConstantPolynomial(distinct.quotient))
            factors.push_back(RationalFactor{
                std::move(distinct.quotient), multiplicity});

        remaining = std::move(shared);
        PolynomialDivision nextRepeated = dividePolynomials(repeated, remaining);
        if (!nextRepeated.remainder.isZero())
            return std::nullopt;
        repeated = std::move(nextRepeated.quotient);
    }

    if (factors.empty())
        return std::nullopt;

    RationalPolynomial reconstructed{{Rational{BigInt{1}}}};
    for (const RationalFactor& factor : factors)
        reconstructed = multiplyPolynomials(
            reconstructed, powerPolynomial(factor.polynomial, factor.multiplicity));
    if (reconstructed.degree() != polynomial.degree())
        return std::nullopt;
    for (std::size_t exponent = 0; exponent <= polynomial.degree(); ++exponent)
        if (!(reconstructed.coefficient(exponent) == polynomial.coefficient(exponent)))
            return std::nullopt;
    return factors;
}

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


[[nodiscard]] bool samePolynomial(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    if (lhs.degree() != rhs.degree())
        return false;
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        if (!(lhs.coefficient(i) == rhs.coefficient(i)))
            return false;
    return true;
}

[[nodiscard]] std::optional<std::vector<RationalFactor>> factorIntoRationalFactorsViaFactorEngine(
    const RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Expr source = polynomialToExpandedExpr(polynomial, variable, builtins);
    const Expr factored = factorExpression(source, builtins, mathematics, angles);

    std::vector<RationalFactor> factors;
    Rational scalar{BigInt{1}};
    const auto appendFactor = [&](const Expr& factor, std::size_t multiplicity,
                                  auto&& appendFactorRef) -> bool {
        if (multiplicity == 0)
            return true;
        if (isHead(factor, builtins, BuiltinId::Multiply)) {
            if (multiplicity != 1)
                return false;
            for (const Expr& nested : factor.asCall().arguments)
                if (!appendFactorRef(nested, 1, appendFactorRef))
                    return false;
            return true;
        }
        if (isHead(factor, builtins, BuiltinId::Power)
            && factor.asCall().arguments.size() == 2) {
            const auto exponent = exactRealRational(factor.asCall().arguments[1]);
            if (!exponent || !exponent->isInteger() || exponent->numerator().isNegative())
                return false;
            const auto count = numeric::tryToUint64(exponent->numerator());
            if (!count || *count == 0 || *count > 64)
                return false;
            return appendFactorRef(
                factor.asCall().arguments[0],
                multiplicity * static_cast<std::size_t>(*count),
                appendFactorRef);
        }
        if (!containsVariable(factor, variable)) {
            const auto value = exactRealRational(factor);
            if (!value)
                return false;
            for (std::size_t i = 0; i < multiplicity; ++i)
                scalar *= *value;
            return true;
        }

        const auto part = toRationalPolynomial(
            factor, variable, builtins, PolynomialConversionOptions{64, 1024});
        if (!part || part->degree() == 0 || part->degree() > 12)
            return false;
        for (RationalFactor& existing : factors) {
            if (samePolynomial(existing.polynomial, *part)) {
                existing.multiplicity += multiplicity;
                return true;
            }
        }
        factors.push_back(RationalFactor{*part, multiplicity});
        return true;
    };

    if (!appendFactor(factored, 1, appendFactor) || factors.empty())
        return std::nullopt;

    RationalPolynomial reconstructed{{scalar}};
    for (const RationalFactor& factor : factors)
        reconstructed = multiplyPolynomials(
            reconstructed, powerPolynomial(factor.polynomial, factor.multiplicity));
    if (!samePolynomial(reconstructed, polynomial))
        return std::nullopt;
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

[[nodiscard]] std::optional<Expr> integrateInverseQuadraticPower(
    const RationalPolynomial& quadratic,
    std::size_t powerIndex,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (quadratic.degree() != 2 || powerIndex == 0)
        return std::nullopt;
    const Rational a = quadratic.coefficient(2);
    const Rational b = quadratic.coefficient(1);
    const Rational c = quadratic.coefficient(0);
    const Rational delta = Rational{BigInt{4}} * a * c - b * b;
    if (a.isZero() || delta.isZero())
        return std::nullopt;

    const Expr q = polynomialToExpandedExpr(quadratic, variable, builtins);
    Expr currentFraction = call(builtins, BuiltinId::Divide, {integer(1), q});
    auto current = tryRationalLowDegree(
        currentFraction, variable, builtins, mathematics, angles);
    if (!current)
        return std::nullopt;
    if (powerIndex == 1)
        return current;

    // I_k = (2ax+b)/(delta(k-1) q^(k-1))
    //     + 2a(2k-3)/(delta(k-1)) I_(k-1).
    // Completing the square proves this recurrence over Q without numerical roots.
    for (std::size_t k = 2; k <= powerIndex; ++k) {
        const Rational km1{BigInt::fromUnsigned(k - 1)};
        Expr linear = add(builtins, mathematics, angles, {
            multiply(builtins, mathematics, angles, {
                rational(Rational{BigInt{2}} * a), Expr{variable}}),
            rational(b)});
        Expr first = divide(builtins, mathematics, angles,
            std::move(linear),
            multiply(builtins, mathematics, angles, {
                rational(delta * km1),
                power(builtins, mathematics, angles, q,
                    integer(static_cast<std::int64_t>(k - 1)))}));
        const Rational recurrence = Rational{BigInt{2}} * a
            * Rational{BigInt::fromUnsigned(2 * k - 3)} / (delta * km1);
        *current = add(builtins, mathematics, angles, {
            std::move(first),
            multiply(builtins, mathematics, angles, {
                rational(recurrence), std::move(*current)})});
    }
    return current;
}

[[nodiscard]] std::optional<Expr> integrateQuadraticPartialFraction(
    const RationalPolynomial& numerator,
    const RationalPolynomial& quadratic,
    std::size_t denominatorPower,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (quadratic.degree() != 2 || numerator.degree() > 1 || denominatorPower == 0)
        return std::nullopt;
    const Rational a = quadratic.coefficient(2);
    const Rational b = quadratic.coefficient(1);
    if (a.isZero())
        return std::nullopt;
    const Rational linear = numerator.coefficient(1);
    const Rational constant = numerator.coefficient(0);
    const Rational alpha = linear / (Rational{BigInt{2}} * a);
    const Rational beta = constant - alpha * b;
    const Expr q = polynomialToExpandedExpr(quadratic, variable, builtins);

    std::vector<Expr> terms;
    if (!alpha.isZero()) {
        if (denominatorPower == 1)
            terms.push_back(multiply(builtins, mathematics, angles, {
                rational(alpha), call(builtins, BuiltinId::Log, {q})}));
        else {
            const Rational exponent{BigInt{1 - static_cast<std::int64_t>(denominatorPower)}};
            terms.push_back(multiply(builtins, mathematics, angles, {
                rational(alpha / exponent),
                power(builtins, mathematics, angles, q, rational(exponent))}));
        }
    }
    if (!beta.isZero()) {
        auto inversePower = integrateInverseQuadraticPower(
            quadratic, denominatorPower, variable, builtins, mathematics, angles);
        if (!inversePower)
            return std::nullopt;
        terms.push_back(multiply(builtins, mathematics, angles, {
            rational(beta), std::move(*inversePower)}));
    }
    if (terms.empty())
        return integer(0);
    return add(builtins, mathematics, angles, std::move(terms));
}

struct HermitePowerReduction final {
    std::vector<std::pair<RationalPolynomial, std::size_t>> rationalTerms;
    RationalPolynomial squareFreeNumerator;
};

[[nodiscard]] std::optional<HermitePowerReduction> hermiteReduceSquareFreePower(
    RationalPolynomial numerator,
    const RationalPolynomial& squareFreeFactor,
    std::size_t denominatorPower,
    const RationalPolynomial& inverseDerivativeModuloFactor) {
    if (denominatorPower == 0)
        return std::nullopt;

    HermitePowerReduction result;
    const RationalPolynomial derivative = differentiatePolynomial(squareFreeFactor);
    while (denominatorPower > 1 && !numerator.isZero()) {
        PolynomialDivision inverseProduct = dividePolynomials(
            multiplyPolynomials(numerator, inverseDerivativeModuloFactor),
            squareFreeFactor);
        const Rational scale = -Rational{BigInt{1}}
            / Rational{BigInt::fromUnsigned(denominatorPower - 1)};
        RationalPolynomial correction = scalePolynomial(inverseProduct.remainder, scale);

        // d(B/f^(k-1)) = (B' f - (k-1) B f') / f^k.
        // B = -A (f')^-1/(k-1) mod f と取ると残差はexactにfで割れ，
        // 分母冪を1段下げられる。数値rootや因数分解は使わないHermite stepである。
        RationalPolynomial residual = addPolynomials(numerator,
            addPolynomials(
                negatePolynomial(multiplyPolynomials(
                    differentiatePolynomial(correction), squareFreeFactor)),
                scalePolynomial(
                    multiplyPolynomials(correction, derivative),
                    Rational{BigInt::fromUnsigned(denominatorPower - 1)})));
        PolynomialDivision lowered = dividePolynomials(residual, squareFreeFactor);
        if (!lowered.remainder.isZero())
            return std::nullopt;

        if (!correction.isZero())
            result.rationalTerms.emplace_back(
                std::move(correction), denominatorPower - 1);
        numerator = std::move(lowered.quotient);
        --denominatorPower;
    }
    result.squareFreeNumerator = std::move(numerator);
    return result;
}

[[nodiscard]] std::optional<Expr> integrateSquareFreeAlgebraicLog(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (numerator.isZero())
        return integer(0);
    if (denominator.degree() < 3 || denominator.degree() > 12)
        return std::nullopt;

    const RationalPolynomial derivative = differentiatePolynomial(denominator);
    if (polynomialGcdMonic(denominator, derivative).degree() != 0)
        return std::nullopt;

    auto roots = ComplexAlgebraicNumber::isolateAll(denominator.coefficients());
    if (!roots || roots->size() != denominator.degree())
        return std::nullopt;
    if (auto canonical = ComplexAlgebraicNumber::canonicalizeAll(*roots))
        roots = std::move(canonical);

    std::vector<Expr> terms;
    terms.reserve(roots->size());
    for (const ComplexAlgebraicNumber& root : *roots) {
        const AlgebraicNumber point =
            AlgebraicNumber::fromComplexRoot(root).withGeneratorField();
        const auto numeratorAtRoot = evaluatePolynomialAtAlgebraic(numerator, point);
        const auto derivativeAtRoot = evaluatePolynomialAtAlgebraic(derivative, point);
        if (!numeratorAtRoot || !derivativeAtRoot)
            return std::nullopt;
        const auto residue = AlgebraicNumber::combine(
            *numeratorAtRoot, *derivativeAtRoot, AlgebraicBinaryOperation::Divide);
        if (!residue)
            return std::nullopt;
        const auto exactResidue = residue->exactRationalParts();
        if (exactResidue && exactResidue->first.isZero() && exactResidue->second.isZero())
            continue;

        Expr residueExpr = makeCanonicalAlgebraicExpression(*residue, builtins);
        Expr rootExpr = makeCanonicalRootExpression(root, builtins);
        Expr logarithm = call(builtins, BuiltinId::Log, {
            subtract(builtins, mathematics, angles, Expr{variable}, std::move(rootExpr))});
        terms.push_back(multiply(builtins, mathematics, angles, {
            std::move(residueExpr), std::move(logarithm)}));
    }

    if (terms.empty())
        return integer(0);
    // square-free Q に対する P/Q = sum P(r)/(Q'(r)(x-r)) はexactな留数分解。
    // Rootは各rをcertifiedに識別するため，近似候補を恒等式の根拠には使わない。
    return add(builtins, mathematics, angles, std::move(terms));
}

[[nodiscard]] std::optional<Expr> tryRationalFactored(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    bool allowAlgebraicFactors = false) {
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
    auto factors = factorIntoLinearQuadratic(*denominator);
    if (!factors)
        factors = factorIntoRationalFactorsViaFactorEngine(
            *denominator, variable, builtins, mathematics, angles);
    if (factors && allowAlgebraicFactors) {
        for (const RationalFactor& factor : *factors) {
            if (polynomialGcdMonic(
                    factor.polynomial, differentiatePolynomial(factor.polynomial)).degree() != 0) {
                factors.reset();
                break;
            }
        }
    }
    if (!factors && allowAlgebraicFactors)
        factors = factorIntoSquareFreePowers(*denominator);
    if (!factors)
        return std::nullopt;

    std::vector<PartialFractionBasis> basis;
    for (std::size_t i = 0; i < factors->size(); ++i) {
        const std::size_t degree = (*factors)[i].polynomial.degree();
        if (degree == 0 || degree > 12 || (!allowAlgebraicFactors && degree > 2))
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

        if (factor.polynomial.degree() == 2) {
            for (std::size_t k = 1; k <= factor.multiplicity; ++k) {
                std::vector<Rational> numeratorCoefficients(2, Rational{BigInt{0}});
                for (std::size_t column = 0; column < basis.size(); ++column) {
                    if (basis[column].factorIndex == factorIndex
                        && basis[column].denominatorPower == k)
                        numeratorCoefficients[basis[column].numeratorExponent]
                            += (*coefficients)[column];
                }
                RationalPolynomial partialNumerator{std::move(numeratorCoefficients)};
                if (partialNumerator.isZero())
                    continue;
                const auto primitive = integrateQuadraticPartialFraction(
                    partialNumerator, factor.polynomial, k,
                    variable, builtins, mathematics, angles);
                if (!primitive)
                    return std::nullopt;
                primitiveTerms.push_back(*primitive);
            }
            continue;
        }

        const RationalPolynomial derivative = differentiatePolynomial(factor.polynomial);
        if (polynomialGcdMonic(factor.polynomial, derivative).degree() != 0)
            return std::nullopt;
        const auto inverseDerivative = polynomialInverseModulo(derivative, factor.polynomial);
        if (!inverseDerivative)
            return std::nullopt;

        RationalPolynomial squareFreeNumerator;
        std::vector<RationalPolynomial> rationalNumerators(
            factor.multiplicity, RationalPolynomial{});
        for (std::size_t k = 1; k <= factor.multiplicity; ++k) {
            std::vector<Rational> numeratorCoefficients(
                factor.polynomial.degree(), Rational{BigInt{0}});
            for (std::size_t column = 0; column < basis.size(); ++column) {
                if (basis[column].factorIndex == factorIndex
                    && basis[column].denominatorPower == k)
                    numeratorCoefficients[basis[column].numeratorExponent]
                        += (*coefficients)[column];
            }
            RationalPolynomial partialNumerator{std::move(numeratorCoefficients)};
            if (partialNumerator.isZero())
                continue;
            if (k == 1) {
                squareFreeNumerator = addPolynomials(squareFreeNumerator, partialNumerator);
                continue;
            }

            const auto reduced = hermiteReduceSquareFreePower(
                partialNumerator, factor.polynomial, k, *inverseDerivative);
            if (!reduced)
                return std::nullopt;
            squareFreeNumerator = addPolynomials(
                squareFreeNumerator, reduced->squareFreeNumerator);
            for (const auto& [rationalNumerator, denominatorPower] : reduced->rationalTerms) {
                if (denominatorPower == 0 || denominatorPower >= rationalNumerators.size())
                    return std::nullopt;
                rationalNumerators[denominatorPower] = addPolynomials(
                    rationalNumerators[denominatorPower], rationalNumerator);
            }
        }

        for (std::size_t denominatorPower = 1;
             denominatorPower < rationalNumerators.size(); ++denominatorPower) {
            if (rationalNumerators[denominatorPower].isZero())
                continue;
            Expr rationalNumerator = polynomialToExpandedExpr(
                rationalNumerators[denominatorPower], variable, builtins);
            Expr rationalDenominator = power(
                builtins, mathematics, angles, factorExpr,
                integer(static_cast<std::int64_t>(denominatorPower)));
            primitiveTerms.push_back(divide(
                builtins, mathematics, angles,
                std::move(rationalNumerator), std::move(rationalDenominator)));
        }

        if (!squareFreeNumerator.isZero()) {
            const auto logarithmic = integrateSquareFreeAlgebraicLog(
                squareFreeNumerator, factor.polynomial,
                variable, builtins, mathematics, angles);
            if (!logarithmic)
                return std::nullopt;
            primitiveTerms.push_back(*logarithmic);
        }
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
            // normalizeRationalFunctionがcontentを外へ出した場合でも，
            // Logの引数は元の一次多項式の展開形を使う。
            // 定数倍の違いは積分定数へ吸収できるが，2(x+3/2)より2x+3をcanonicalとする。
            Expr logArgument = polynomialToExpandedExpr(*denominator, variable, builtins);
            primitiveTerms.push_back(multiply(builtins, mathematics, angles, {
                rational(remainder / slope),
                call(builtins, BuiltinId::Log, {std::move(logArgument)})}));
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
    return simplify(Expr::rebuildCall(expression.asCall(), std::move(arguments)),
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

// exact Rational係数の二次式 q(x)=a x^2+b x+c に対する sqrt[q(x)] の局所原始函数。
// 逆平方根は以前から扱えていたが、sqrt[q] 本体が未実装だったため
// sqrt[4-x^2], sqrt[x^2+4], sqrt[x^2-4] が不自然に未評価で残っていた。
// principal sqrt/log/asin の大域的な恒等変形は行わず、標準の平方完成公式から
// primitiveを直接構成する。
[[nodiscard]] std::optional<Expr> tryQuadraticSquareRoot(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Sqrt)
        || expression.asCall().arguments.size() != 1)
        return std::nullopt;

    const Expr& radicand = expression.asCall().arguments[0];
    const auto coefficients = exactQuadraticCoefficients(radicand, variable, builtins);
    if (!coefficients || (*coefficients)[2].isZero())
        return std::nullopt;

    const Rational c = (*coefficients)[0];
    const Rational b = (*coefficients)[1];
    const Rational a = (*coefficients)[2];
    const Rational discriminant = b * b - Rational{BigInt{4}} * a * c;
    Expr sqrtQ = call(builtins, BuiltinId::Sqrt, {radicand});
    Expr linear = quadraticLinearNumerator(
        *coefficients, variable, builtins, mathematics, angles);
    Expr algebraic = divide(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles, {linear, sqrtQ}),
        rational(Rational{BigInt{4}} * a));

    if (a > Rational{BigInt{0}}) {
        if (discriminant.isZero())
            return algebraic;

        Expr sqrtA = call(builtins, BuiltinId::Sqrt, {rational(a)});
        Expr logArgument = add(builtins, mathematics, angles, {
            multiply(builtins, mathematics, angles,
                {integer(2), sqrtA, sqrtQ}),
            linear});
        Expr logCoefficient = divide(builtins, mathematics, angles,
            rational(-discriminant),
            multiply(builtins, mathematics, angles,
                {integer(8), rational(a), sqrtA}));
        return add(builtins, mathematics, angles, {
            std::move(algebraic),
            multiply(builtins, mathematics, angles, {
                std::move(logCoefficient),
                call(builtins, BuiltinId::Log, {std::move(logArgument)})})});
    }

    if (a < Rational{BigInt{0}} && discriminant > Rational{BigInt{0}}) {
        Expr scale = call(builtins, BuiltinId::Sqrt, {rational(discriminant)});
        Expr argument = divide(builtins, mathematics, angles,
            negate(builtins, mathematics, angles, linear), scale);
        Expr coefficient = divide(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {
                rational(-discriminant),
                radiansPerInverseAngleUnit(builtins, mathematics, angles)}),
            multiply(builtins, mathematics, angles, {
                integer(8), rational(a),
                call(builtins, BuiltinId::Sqrt, {rational(-a)})}));
        return add(builtins, mathematics, angles, {
            std::move(algebraic),
            multiply(builtins, mathematics, angles, {
                std::move(coefficient), call(builtins, BuiltinId::Asin, {std::move(argument)})})});
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



[[nodiscard]] std::optional<Expr> integrateSineCosineOverArgument(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const auto& quotient = expression.asCall().arguments;
    if (!quotient[0].isCall() || quotient[0].asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = builtins.find(quotient[0].asCall().head);
    if (!definition
        || (definition->id != BuiltinId::Sin && definition->id != BuiltinId::Cos))
        return std::nullopt;

    // sin(u)/u と cos(u)/u の核は，角度単位をRadianへ正規化したSi/Ciで閉じる。
    // u'が定数の場合だけ採用し，非線形uでは余分な1/u'を記号的に持ち込まない。
    const TrigArgument info = trigArgument(
        quotient[0].asCall().arguments[0], builtins, mathematics, angles);
    if (quotient[1] != info.argument)
        return std::nullopt;

    Expr du = simplify(
        differentiateExpression(info.argument, variable, builtins, mathematics, angles),
        builtins, mathematics, angles);
    if (containsVariable(du, variable) || !provablyNonZero(du, builtins, mathematics))
        return std::nullopt;

    Expr radians = multiply(
        builtins, mathematics, angles, {info.scale, info.argument});
    Expr primitive = call(builtins,
        definition->id == BuiltinId::Sin
            ? BuiltinId::SineIntegralSi
            : BuiltinId::CosineIntegralCi,
        {std::move(radians)});
    return divide(
        builtins, mathematics, angles, std::move(primitive), std::move(du));
}

struct EllipticTrigKernel final {
    Expr parameter;
    Expr sourceArgument;
};

[[nodiscard]] std::optional<EllipticTrigKernel> matchOneMinusParameterSinSquared(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Subtract)
        || expression.asCall().arguments.size() != 2
        || !isOne(expression.asCall().arguments[0]))
        return std::nullopt;

    const Expr& term = expression.asCall().arguments[1];
    std::vector<std::pair<Expr, Expr>> candidates;
    if (isHead(term, builtins, BuiltinId::Multiply)
        && term.asCall().arguments.size() == 2) {
        const auto& factors = term.asCall().arguments;
        candidates.emplace_back(factors[0], factors[1]);
        candidates.emplace_back(factors[1], factors[0]);
    }
    else if (isHead(term, builtins, BuiltinId::Divide)
        && term.asCall().arguments.size() == 2
        && !containsVariable(term.asCall().arguments[1], variable)) {
        candidates.emplace_back(
            Expr::call(builtins.symbol(BuiltinId::Divide), {
                integer(1), term.asCall().arguments[1]}),
            term.asCall().arguments[0]);
    }
    else {
        candidates.emplace_back(integer(1), term);
    }

    for (const auto& [parameter, sinePower] : candidates) {
        if (containsVariable(parameter, variable)
            || !isHead(sinePower, builtins, BuiltinId::Power)
            || sinePower.asCall().arguments.size() != 2)
            continue;
        const auto exponent = exactRealRational(sinePower.asCall().arguments[1]);
        if (!exponent || *exponent != Rational{BigInt{2}})
            continue;
        const Expr& sine = sinePower.asCall().arguments[0];
        if (!isHead(sine, builtins, BuiltinId::Sin)
            || sine.asCall().arguments.size() != 1)
            continue;
        return EllipticTrigKernel{parameter, sine.asCall().arguments[0]};
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> integrateEllipticTrigKernel(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    BuiltinId target = BuiltinId::EllipticF;
    Expr characteristic = integer(0);
    const Expr* radicand = nullptr;
    std::optional<EllipticTrigKernel> matched;

    if (isHead(expression, builtins, BuiltinId::Sqrt)
        && expression.asCall().arguments.size() == 1) {
        target = BuiltinId::EllipticE;
        radicand = &expression.asCall().arguments[0];
        matched = matchOneMinusParameterSinSquared(*radicand, variable, builtins);
    }
    else if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2
        && isOne(expression.asCall().arguments[0])) {
        const Expr& denominator = expression.asCall().arguments[1];
        if (isHead(denominator, builtins, BuiltinId::Sqrt)
            && denominator.asCall().arguments.size() == 1) {
            target = BuiltinId::EllipticF;
            radicand = &denominator.asCall().arguments[0];
            matched = matchOneMinusParameterSinSquared(*radicand, variable, builtins);
        }
        else if (isHead(denominator, builtins, BuiltinId::Multiply)) {
            const auto& factors = denominator.asCall().arguments;
            if (factors.size() == 2) {
                for (std::size_t rootIndex = 0; rootIndex < 2; ++rootIndex) {
                    const Expr& root = factors[rootIndex];
                    if (!isHead(root, builtins, BuiltinId::Sqrt)
                        || root.asCall().arguments.size() != 1)
                        continue;
                    auto mKernel = matchOneMinusParameterSinSquared(
                        root.asCall().arguments[0], variable, builtins);
                    auto nKernel = matchOneMinusParameterSinSquared(
                        factors[1 - rootIndex], variable, builtins);
                    if (mKernel && nKernel
                        && mKernel->sourceArgument == nKernel->sourceArgument) {
                        target = BuiltinId::EllipticPi;
                        matched = std::move(mKernel);
                        characteristic = std::move(nKernel->parameter);
                        break;
                    }
                }
            }
        }
    }
    if (!matched)
        return std::nullopt;

    TrigArgument info = trigArgument(
        matched->sourceArgument, builtins, mathematics, angles);
    Expr du = differentiateExpression(info.argument, variable, builtins, mathematics, angles);
    if (containsVariable(du, variable) || !provablyNonZero(du, builtins, mathematics))
        return std::nullopt;

    Expr amplitude = isOne(info.scale)
        ? info.argument
        : multiply(builtins, mathematics, angles, {info.scale, info.argument});
    Expr primitive = target == BuiltinId::EllipticPi
        ? call(builtins, target, {
            std::move(characteristic), std::move(amplitude), matched->parameter})
        : call(builtins, target, {std::move(amplitude), matched->parameter});
    return multiply(builtins, mathematics, angles, {
        std::move(info.inverseScale),
        divide(builtins, mathematics, angles, std::move(primitive), std::move(du))});
}

[[nodiscard]] std::optional<std::uint64_t> negativeIntegerMagnitude(const Expr& expression) {
    const auto value = exactRealRational(expression);
    if (!value || !value->isInteger() || !value->numerator().isNegative())
        return std::nullopt;
    return numeric::tryToUint64(-value->numerator());
}

[[nodiscard]] std::optional<std::uint64_t> positiveIntegerMagnitude(const Expr& expression) {
    const auto value = exactRealRational(expression);
    if (!value || !value->isInteger() || !value->numerator().isPositive())
        return std::nullopt;
    return numeric::tryToUint64(value->numerator());
}

[[nodiscard]] std::optional<Expr> integrateReciprocalTrigPower(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Power)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const auto& powerArguments = expression.asCall().arguments;
    const auto order = negativeIntegerMagnitude(powerArguments[1]);
    if (!order || *order == 0 || *order > 256)
        return std::nullopt;

    const Expr& base = powerArguments[0];
    if (!base.isCall() || base.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = builtins.find(base.asCall().head);
    if (!definition || (definition->id != BuiltinId::Sin && definition->id != BuiltinId::Cos))
        return std::nullopt;

    const Expr& sourceArgument = base.asCall().arguments[0];
    TrigArgument info = trigArgument(sourceArgument, builtins, mathematics, angles);
    Expr du = differentiateExpression(info.argument, variable, builtins, mathematics, angles);
    if (containsVariable(du, variable) || !provablyNonZero(du, builtins, mathematics))
        return std::nullopt;
    Expr inverseScale = divide(
        builtins, mathematics, angles, std::move(info.inverseScale), std::move(du));

    const bool sine = definition->id == BuiltinId::Sin;
    const BuiltinId reciprocalId = sine ? BuiltinId::Csc : BuiltinId::Sec;
    const BuiltinId companionId = sine ? BuiltinId::Cot : BuiltinId::Tan;
    const auto same = [&](BuiltinId id) { return call(builtins, id, {sourceArgument}); };

    // 旧実装ではPowerの一般ruleが base'=const の場合しか扱えず、sin[2x]^-2 は
    // base'=2 cos[2x] が変数依存なので未評価だった。負整数冪はcsc/secの標準漸化式として
    // MathKnowledge化し、-2だけの個別hackではなく -1..-256 を同じ規則で扱う。
    std::vector<std::optional<Expr>> primitives(static_cast<std::size_t>(*order) + 1);
    if (sine) {
        primitives[1] = negate(builtins, mathematics, angles,
            call(builtins, BuiltinId::Log, {
                add(builtins, mathematics, angles, {same(BuiltinId::Csc), same(BuiltinId::Cot)})}));
        if (*order >= 2)
            primitives[2] = negate(builtins, mathematics, angles, same(BuiltinId::Cot));
    }
    else {
        primitives[1] = call(builtins, BuiltinId::Log, {
            add(builtins, mathematics, angles, {same(BuiltinId::Sec), same(BuiltinId::Tan)})});
        if (*order >= 2)
            primitives[2] = same(BuiltinId::Tan);
    }

    for (std::uint64_t n = 3; n <= *order; ++n) {
        Expr reciprocalPower = n == 2
            ? same(reciprocalId)
            : power(builtins, mathematics, angles, same(reciprocalId),
                integer(static_cast<std::int64_t>(n - 2)));
        Expr boundary = multiply(builtins, mathematics, angles,
            {std::move(reciprocalPower), same(companionId)});
        if (sine)
            boundary = negate(builtins, mathematics, angles, std::move(boundary));
        boundary = divide(builtins, mathematics, angles,
            std::move(boundary), integer(static_cast<std::int64_t>(n - 1)));
        Expr recurrence = multiply(builtins, mathematics, angles, {
            rational(Rational{BigInt::fromUnsigned(n - 2), BigInt::fromUnsigned(n - 1)}),
            *primitives[static_cast<std::size_t>(n - 2)]});
        primitives[static_cast<std::size_t>(n)] = add(
            builtins, mathematics, angles, {std::move(boundary), std::move(recurrence)});
    }

    return multiply(builtins, mathematics, angles, {
        std::move(inverseScale), *primitives[static_cast<std::size_t>(*order)]});
}


[[nodiscard]] std::optional<Expr> integrateDirectTrigPower(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Power)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const auto& a = expression.asCall().arguments;
    const auto order = positiveIntegerMagnitude(a[1]);
    if (!order || *order < 2 || *order > 256)
        return std::nullopt;
    const Expr& base = a[0];
    if (!base.isCall() || base.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = builtins.find(base.asCall().head);
    if (!definition)
        return std::nullopt;

    const BuiltinId id = definition->id;
    if (id != BuiltinId::Tan && id != BuiltinId::Cot
        && id != BuiltinId::Sec && id != BuiltinId::Csc)
        return std::nullopt;

    const Expr& sourceArgument = base.asCall().arguments[0];
    TrigArgument info = trigArgument(sourceArgument, builtins, mathematics, angles);
    Expr du = differentiateExpression(info.argument, variable, builtins, mathematics, angles);
    if (containsVariable(du, variable) || !provablyNonZero(du, builtins, mathematics))
        return std::nullopt;
    Expr inverseScale = divide(
        builtins, mathematics, angles, std::move(info.inverseScale), std::move(du));

    const auto same = [&](BuiltinId function) {
        return call(builtins, function, {sourceArgument});
    };

    std::vector<std::optional<Expr>> primitives(static_cast<std::size_t>(*order) + 1);
    if (id == BuiltinId::Tan || id == BuiltinId::Cot) {
        primitives[0] = sourceArgument;
        primitives[1] = id == BuiltinId::Tan
            ? negate(builtins, mathematics, angles,
                call(builtins, BuiltinId::Log, {same(BuiltinId::Cos)}))
            : call(builtins, BuiltinId::Log, {same(BuiltinId::Sin)});

        // tan^n = tan^(n-2)(sec^2-1),
        // cot^n = cot^(n-2)(csc^2-1) を使う標準漸化式。
        for (std::uint64_t n = 2; n <= *order; ++n) {
            Expr boundary = power(builtins, mathematics, angles, same(id),
                integer(static_cast<std::int64_t>(n - 1)));
            boundary = divide(builtins, mathematics, angles,
                std::move(boundary), integer(static_cast<std::int64_t>(n - 1)));
            if (id == BuiltinId::Cot)
                boundary = negate(builtins, mathematics, angles, std::move(boundary));
            primitives[static_cast<std::size_t>(n)] = subtract(
                builtins, mathematics, angles,
                std::move(boundary), *primitives[static_cast<std::size_t>(n - 2)]);
        }
    }
    else {
        const bool secant = id == BuiltinId::Sec;
        const BuiltinId companion = secant ? BuiltinId::Tan : BuiltinId::Cot;
        primitives[1] = secant
            ? call(builtins, BuiltinId::Log, {
                add(builtins, mathematics, angles, {same(BuiltinId::Sec), same(BuiltinId::Tan)})})
            : negate(builtins, mathematics, angles,
                call(builtins, BuiltinId::Log, {
                    add(builtins, mathematics, angles, {same(BuiltinId::Csc), same(BuiltinId::Cot)})}));
        primitives[2] = secant ? same(BuiltinId::Tan)
            : negate(builtins, mathematics, angles, same(BuiltinId::Cot));

        // sec/cscの標準 reduction formula。既存 sin^-n/cos^-n と同じ数学Knowledgeを
        // 直接函数表記にも適用し、sec^3/csc^3 を個別表で持たない。
        for (std::uint64_t n = 3; n <= *order; ++n) {
            Expr reciprocalPower = power(builtins, mathematics, angles, same(id),
                integer(static_cast<std::int64_t>(n - 2)));
            Expr boundary = multiply(builtins, mathematics, angles,
                {std::move(reciprocalPower), same(companion)});
            if (!secant)
                boundary = negate(builtins, mathematics, angles, std::move(boundary));
            boundary = divide(builtins, mathematics, angles,
                std::move(boundary), integer(static_cast<std::int64_t>(n - 1)));
            Expr recurrence = multiply(builtins, mathematics, angles, {
                rational(Rational{BigInt::fromUnsigned(n - 2), BigInt::fromUnsigned(n - 1)}),
                *primitives[static_cast<std::size_t>(n - 2)]});
            primitives[static_cast<std::size_t>(n)] = add(
                builtins, mathematics, angles, {std::move(boundary), std::move(recurrence)});
        }
    }

    return multiply(builtins, mathematics, angles, {
        std::move(inverseScale), *primitives[static_cast<std::size_t>(*order)]});
}


[[nodiscard]] std::optional<Expr> integrateQuarticEllipticF(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Expr* radicand = nullptr;
    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2
        && isOne(expression.asCall().arguments[0])
        && isHead(expression.asCall().arguments[1], builtins, BuiltinId::Sqrt)
        && expression.asCall().arguments[1].asCall().arguments.size() == 1)
        radicand = &expression.asCall().arguments[1].asCall().arguments[0];
    if (!radicand)
        return std::nullopt;

    const auto polynomial = toRationalPolynomial(
        *radicand, variable, builtins, PolynomialConversionOptions{4, 16});
    if (!polynomial || polynomial->degree() != 4
        || polynomial->coefficient(0) != Rational{BigInt{1}}
        || !polynomial->coefficient(1).isZero()
        || !polynomial->coefficient(2).isZero()
        || !polynomial->coefficient(3).isZero()
        || polynomial->coefficient(4) != Rational{BigInt{-1}})
        return std::nullopt;

    // ellipticFのamplitudeはRadian固定。asinのsession角度単位をRadianへ戻す。
    Expr amplitude = call(builtins, BuiltinId::Asin, {Expr{variable}});
    Expr radianScale = radiansPerInverseAngleUnit(builtins, mathematics, angles);
    if (!isOne(radianScale))
        amplitude = multiply(builtins, mathematics, angles, {
            std::move(radianScale), std::move(amplitude)});
    return call(builtins, BuiltinId::EllipticF, {
        std::move(amplitude), integer(-1)});
}

[[nodiscard]] std::optional<Expr> integrateBinomialPower2F1(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    // ∫x^m(1+beta*x^n)^p dx
    // = x^(m+1)/(m+1) 2F1(-p,(m+1)/n;1+(m+1)/n;-beta*x^n)。
    // ここではmを非負整数monomialに限定し、principal branchの既存binomial ruleを
    // numerator付きへexactに拡張する。
    const Expr* base = nullptr;
    std::optional<Rational> exponent;
    Rational outside{BigInt{1}};
    std::size_t numeratorDegree = 0;

    if (isHead(expression, builtins, BuiltinId::Sqrt)
        && expression.asCall().arguments.size() == 1) {
        base = &expression.asCall().arguments[0];
        exponent = Rational{BigInt{1}, BigInt{2}};
    }
    else if (isHead(expression, builtins, BuiltinId::Power)
        && expression.asCall().arguments.size() == 2) {
        base = &expression.asCall().arguments[0];
        exponent = exactRealRational(expression.asCall().arguments[1]);
    }
    else if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        const Expr& numeratorExpr = expression.asCall().arguments[0];
        const Expr& denominator = expression.asCall().arguments[1];
        const auto numerator = toRationalPolynomial(
            numeratorExpr, variable, builtins, PolynomialConversionOptions{4096, 8192});
        if (!numerator)
            return std::nullopt;
        std::optional<std::size_t> monomialDegree;
        for (std::size_t i = 0; i <= numerator->degree(); ++i) {
            if (numerator->coefficient(i).isZero())
                continue;
            if (monomialDegree)
                return std::nullopt;
            monomialDegree = i;
            outside = numerator->coefficient(i);
        }
        if (!monomialDegree)
            return integer(0);
        numeratorDegree = *monomialDegree;

        if (isHead(denominator, builtins, BuiltinId::Sqrt)
            && denominator.asCall().arguments.size() == 1) {
            base = &denominator.asCall().arguments[0];
            exponent = Rational{BigInt{-1}, BigInt{2}};
        }
        else {
            base = &denominator;
            exponent = Rational{BigInt{-1}};
        }
    }
    if (!base || !exponent)
        return std::nullopt;

    const auto polynomial = toRationalPolynomial(
        *base, variable, builtins, PolynomialConversionOptions{4096, 8192});
    if (!polynomial || polynomial->degree() < 2
        || polynomial->coefficient(0) != Rational{BigInt{1}})
        return std::nullopt;
    const std::size_t degree = polynomial->degree();
    for (std::size_t i = 1; i < degree; ++i)
        if (!polynomial->coefficient(i).isZero())
            return std::nullopt;
    const Rational beta = polynomial->coefficient(degree);
    if (beta.isZero() || degree > 4096)
        return std::nullopt;

    const BigInt mPlusOne = BigInt::fromUnsigned(numeratorDegree + 1);
    const Rational b{mPlusOne, BigInt::fromUnsigned(degree)};
    const Rational c = Rational{BigInt{1}} + b;
    Expr xPower = power(builtins, mathematics, angles,
        Expr{variable}, integer(static_cast<std::int64_t>(degree)));
    Expr z = multiply(builtins, mathematics, angles, {
        rational(-beta), std::move(xPower)});
    Expr leading = divide(builtins, mathematics, angles,
        multiply(builtins, mathematics, angles, {
            rational(outside),
            power(builtins, mathematics, angles, Expr{variable}, Expr{Number{mPlusOne}})}),
        Expr{Number{mPlusOne}});
    return multiply(builtins, mathematics, angles, {
        std::move(leading),
        call(builtins, BuiltinId::Hypergeometric2F1, {
            rational(-*exponent), rational(b), rational(c), std::move(z)})});
}

[[nodiscard]] std::optional<Expr> integrateExponentialMonomial1F1(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Exp)
        || expression.asCall().arguments.size() != 1)
        return std::nullopt;

    const Expr& exponent = expression.asCall().arguments[0];
    // Multiplyだけを特別扱いしていた旧実装では -x^4 が Negate[Power[...]] のため
    // exp[-x^4] を1F1へ還元できなかった。定数因子分離はNegate/Divideも既に扱えるので、
    // exponent全体へ一律に適用する。
    FactorSplit split = splitConstantFactor(
        exponent, variable, builtins, mathematics, angles);
    if (isOne(split.dependent))
        return std::nullopt;
    Expr coefficient = std::move(split.constant);
    Expr powerExpression = std::move(split.dependent);

    // 旧実装はローカルFactorSplit::dependentへのpointerをifブロック外へ保持していたため、
    // exp[2x^3]のような係数付き指数でdangling pointerとなりsegfaultした。
    // dependent式を値として所有し、以後の解析中に寿命が切れないようにする。
    if (!isHead(powerExpression, builtins, BuiltinId::Power)
        || powerExpression.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& powerArguments = powerExpression.asCall().arguments;
    if (!powerArguments[0].isSymbol()
        || !powerArguments[0].asSymbol().sameIdentity(variable))
        return std::nullopt;
    const auto order = positiveIntegerMagnitude(powerArguments[1]);
    if (!order || *order < 2 || *order > 4096)
        return std::nullopt;

    // Gaussian exp[-a x^2] (a>0 exact Rational) は一般1F1より erf を preferred form とする。
    // sqrt[a] はprincipal exact rootで，Dにより元のintegrandへ戻る局所原始函数を構成する。
    if (*order == 2) {
        const auto exactCoefficient = exactRealRational(coefficient);
        if (exactCoefficient && *exactCoefficient < Rational{BigInt{0}}) {
            const Rational positiveScale = -*exactCoefficient;
            Expr sqrtScale = call(builtins, BuiltinId::Sqrt, {rational(positiveScale)});
            Expr argument = multiply(builtins, mathematics, angles,
                {sqrtScale, Expr{variable}});
            Expr numerator = multiply(builtins, mathematics, angles, {
                call(builtins, BuiltinId::Sqrt, {pi(mathematics)}),
                call(builtins, BuiltinId::Erf, {std::move(argument)})});
            Expr denominator = multiply(builtins, mathematics, angles,
                {integer(2), std::move(sqrtScale)});
            return divide(
                builtins, mathematics, angles, std::move(numerator), std::move(denominator));
        }
    }

    // ∫ exp(c x^n) dx = x 1F1(1/n;1+1/n;c x^n)。
    // 右辺はx=0でもentireで、incomplete-Gamma表現の見かけのbranch/holeを持ち込まない。
    const Rational a{BigInt{1}, BigInt::fromUnsigned(*order)};
    const Rational b = Rational{BigInt{1}} + a;
    Expr argument = isOne(coefficient)
        ? powerExpression
        : multiply(builtins, mathematics, angles, {coefficient, powerExpression});
    return multiply(builtins, mathematics, angles, {
        Expr{variable},
        call(builtins, BuiltinId::Hypergeometric1F1,
            {rational(a), rational(b), std::move(argument)})});
}

[[nodiscard]] std::optional<Expr> integrateLogOnePlusMonomialOverX(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (!arguments[1].isSymbol()
        || !arguments[1].asSymbol().sameIdentity(variable)
        || !isHead(arguments[0], builtins, BuiltinId::Log)
        || arguments[0].asCall().arguments.size() != 1)
        return std::nullopt;

    const auto polynomial = toRationalPolynomial(
        arguments[0].asCall().arguments[0], variable, builtins,
        PolynomialConversionOptions{4096, 4096});
    if (!polynomial || polynomial->degree() < 1
        || polynomial->coefficient(0) != Rational{BigInt{1}})
        return std::nullopt;
    const std::size_t degree = polynomial->degree();
    for (std::size_t i = 1; i < degree; ++i)
        if (!polynomial->coefficient(i).isZero())
            return std::nullopt;
    const Rational beta = polynomial->coefficient(degree);
    if (beta.isZero() || degree > 4096)
        return std::nullopt;

    Expr xPower = power(builtins, mathematics, angles,
        Expr{variable}, integer(static_cast<std::int64_t>(degree)));
    Expr z = multiply(builtins, mathematics, angles, {
        rational(-beta), std::move(xPower)});
    Expr dilogarithm = call(builtins, BuiltinId::Polylog, {
        integer(2), std::move(z)});
    return divide(builtins, mathematics, angles,
        negate(builtins, mathematics, angles, std::move(dilogarithm)),
        integer(static_cast<std::int64_t>(degree)));
}

[[nodiscard]] std::optional<Expr> integrateQuadraticFresnel(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || (definition->id != BuiltinId::Sin && definition->id != BuiltinId::Cos))
        return std::nullopt;

    const Expr& sourceArgument = expression.asCall().arguments[0];
    TrigArgument info = trigArgument(sourceArgument, builtins, mathematics, angles);

    // Pi*x^2/2 のように二次位相全体へsymbolicな定数scaleが掛かる場合も、
    // polynomial部分だけをRational係数として解析する。旧実装はargument全体を
    // toRationalPolynomialへ渡していたためFresnelの定義形そのものを認識できなかった。
    FactorSplit phaseSplit = splitConstantFactor(
        info.argument, variable, builtins, mathematics, angles);
    Expr phaseScale = std::move(phaseSplit.constant);
    Expr polynomialArgument = std::move(phaseSplit.dependent);
    const mathematics::ValueFacts scaleFacts = mathematics::inferValueFacts(
        phaseScale, builtins, mathematics);
    if (scaleFacts.sign == mathematics::RealSign::Negative) {
        phaseScale = negate(builtins, mathematics, angles, std::move(phaseScale));
        polynomialArgument = negate(
            builtins, mathematics, angles, std::move(polynomialArgument));
    }
    else if (scaleFacts.sign != mathematics::RealSign::Positive && !isOne(phaseScale)) {
        return std::nullopt;
    }

    const auto polynomial = toRationalPolynomial(
        polynomialArgument, variable, builtins, PolynomialConversionOptions{2, 8});
    if (!polynomial || polynomial->degree() != 2 || polynomial->coefficient(2).isZero())
        return std::nullopt;

    const Rational quadratic = polynomial->coefficient(2);
    const Rational linear = polynomial->coefficient(1);
    const Rational constant = polynomial->coefficient(0);
    const bool negative = quadratic < Rational{BigInt{0}};
    const Rational magnitude = negative ? -quadratic : quadratic;

    // q x^2+l x+c = q(x+l/(2q))^2 + c-l^2/(4q)。
    // 旧実装はl=c=0だけを認識していたため、一般的なchirp位相をFresnelへ落とせなかった。
    const Rational shift = linear / (Rational{BigInt{2}} * quadratic);
    const Rational phaseOffset = constant
        - linear * linear / (Rational{BigInt{4}} * quadratic);
    Expr shiftedVariable = shift.isZero()
        ? Expr{variable}
        : add(builtins, mathematics, angles, {Expr{variable}, rational(shift)});

    Expr effectiveMagnitude = multiply(builtins, mathematics, angles,
        {info.scale, phaseScale, rational(magnitude)});
    Expr twoEffective = multiply(builtins, mathematics, angles,
        {integer(2), effectiveMagnitude});
    Expr fresnelScale = call(builtins, BuiltinId::Sqrt, {
        divide(builtins, mathematics, angles, twoEffective, pi(mathematics))});
    Expr fresnelArgument = multiply(builtins, mathematics, angles,
        {fresnelScale, shiftedVariable});
    Expr prefactor = divide(builtins, mathematics, angles, integer(1), fresnelScale);

    Expr cPart = multiply(builtins, mathematics, angles, {
        prefactor, call(builtins, BuiltinId::FresnelC, {fresnelArgument})});
    Expr sPart = multiply(builtins, mathematics, angles, {
        prefactor, call(builtins, BuiltinId::FresnelS, {fresnelArgument})});

    if (phaseOffset.isZero()) {
        if (definition->id == BuiltinId::Cos)
            return cPart;
        return negative
            ? negate(builtins, mathematics, angles, std::move(sPart))
            : std::optional<Expr>{std::move(sPart)};
    }

    // phaseOffsetはinfo.scale適用後のRadian量。session angle modeに再解釈させない。
    Expr delta = multiply(builtins, mathematics, angles,
        {info.scale, phaseScale, rational(phaseOffset)});
    Expr radianDelta = call(builtins, BuiltinId::UnitApplied, {
        std::move(delta), Expr{std::string{"Rad"}}});
    Expr cosDelta = call(builtins, BuiltinId::Cos, {radianDelta});
    Expr sinDelta = call(builtins, BuiltinId::Sin, {radianDelta});

    if (definition->id == BuiltinId::Cos) {
        Expr sineContribution = multiply(
            builtins, mathematics, angles, {std::move(sinDelta), std::move(sPart)});
        if (!negative)
            sineContribution = negate(
                builtins, mathematics, angles, std::move(sineContribution));
        return add(builtins, mathematics, angles, {
            multiply(builtins, mathematics, angles, {std::move(cosDelta), std::move(cPart)}),
            std::move(sineContribution)});
    }

    Expr sineSquareContribution = multiply(
        builtins, mathematics, angles, {std::move(cosDelta), std::move(sPart)});
    if (negative)
        sineSquareContribution = negate(
            builtins, mathematics, angles, std::move(sineSquareContribution));
    return add(builtins, mathematics, angles, {
        multiply(builtins, mathematics, angles, {std::move(sinDelta), std::move(cPart)}),
        std::move(sineSquareContribution)});
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

[[nodiscard]] std::optional<Expr> integrateTrigFourierExpansion(
    mathematics::TrigFourierExpansion expansion,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    TrigArgument info = trigArgument(expansion.argument, builtins, mathematics, angles);
    Expr du = differentiateExpression(
        info.argument, variable, builtins, mathematics, angles);
    if (containsVariable(du, variable) || !provablyNonZero(du, builtins, mathematics))
        return std::nullopt;

    Expr inverseDerivativeScale = divide(
        builtins, mathematics, angles, std::move(info.inverseScale), std::move(du));
    std::vector<Expr> primitives;
    primitives.reserve(expansion.terms.size());

    for (const mathematics::TrigFourierTerm& term : expansion.terms) {
        if (term.frequency == 0) {
            primitives.push_back(call(builtins, BuiltinId::Multiply, {
                rational(term.coefficient), Expr{variable}}));
            continue;
        }

        Expr atom = call(builtins, term.sine ? BuiltinId::Cos : BuiltinId::Sin, {
            mathematics::scaledTrigArgumentForFrequency(
                expansion.argument, term.frequency, builtins)});
        Rational coefficient = term.coefficient
            / Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(term.frequency))};
        if (term.sine)
            coefficient = -coefficient;

        std::vector<Expr> factors;
        if (coefficient != Rational{BigInt{1}})
            factors.push_back(rational(coefficient));
        if (!isOne(inverseDerivativeScale))
            factors.push_back(inverseDerivativeScale);
        factors.push_back(std::move(atom));
        primitives.push_back(factors.size() == 1
            ? std::move(factors.front())
            : call(builtins, BuiltinId::Multiply, std::move(factors)));
    }

    if (primitives.empty())
        return integer(0);
    Expr result = primitives.size() == 1
        ? std::move(primitives.front())
        : call(builtins, BuiltinId::Add, std::move(primitives));
    return result;
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

    // ∫log(log(u)) dx = (u log(log(u)) - li(u))/u' for constant nonzero u'.
    // li' = 1/log(u) なので微分でexactに検証できる局所原始函数である。
    if (definition->id == BuiltinId::Log
        && isHead(u, builtins, BuiltinId::Log)
        && u.asCall().arguments.size() == 1) {
        const Expr& inner = u.asCall().arguments[0];
        Expr innerDerivative = simplify(
            differentiateExpression(inner, variable, builtins, mathematics, angles),
            builtins, mathematics, angles);
        if (!containsVariable(innerDerivative, variable)
            && provablyNonZero(innerDerivative, builtins, mathematics)) {
            Expr numerator = subtract(builtins, mathematics, angles,
                multiply(builtins, mathematics, angles, {inner, expression}),
                call(builtins, BuiltinId::LogarithmicIntegralLi, {inner}));
            return divide(builtins, mathematics, angles,
                std::move(numerator), std::move(innerDerivative));
        }
    }

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

    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS: {
        Expr phase = divide(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {pi(mathematics), u2}), integer(2));
        Expr radianPhase = call(builtins, BuiltinId::UnitApplied, {
            std::move(phase), Expr{std::string{"Rad"}}});
        Expr oscillation = call(builtins,
            definition->id == BuiltinId::FresnelC ? BuiltinId::Sin : BuiltinId::Cos,
            {std::move(radianPhase)});
        Expr correction = divide(builtins, mathematics, angles,
            std::move(oscillation), pi(mathematics));
        Expr main = multiply(builtins, mathematics, angles, {u, expression});
        return divideByDu(definition->id == BuiltinId::FresnelC
            ? subtract(builtins, mathematics, angles, std::move(main), std::move(correction))
            : add(builtins, mathematics, angles, {std::move(main), std::move(correction)}));
    }

    case BuiltinId::Sinc:
    case BuiltinId::Cosc: {
        // cardinal trigは内部でangle operandをRadianへ変換してから sin(y)/y 等を定義する。
        // uに明示UnitAppliedがある場合は上のdu判定で未解決になるため，ここではsession既定単位だけ扱う。
        const Expr angleScale = directTrigScale(
            angles.defaultUnit(), builtins, mathematics, angles);
        Expr radians = multiply(builtins, mathematics, angles, {angleScale, u});
        Expr denominator = multiply(builtins, mathematics, angles, {angleScale, du});
        if (definition->id == BuiltinId::Sinc)
            return divide(builtins, mathematics, angles,
                call(builtins, BuiltinId::SineIntegralSi, {std::move(radians)}),
                std::move(denominator));
        Expr logRadians = call(builtins, BuiltinId::Log, {radians});
        Expr ciRadians = call(builtins, BuiltinId::CosineIntegralCi, {std::move(radians)});
        Expr primitive = subtract(builtins, mathematics, angles,
            std::move(logRadians), std::move(ciRadians));
        return divide(builtins, mathematics, angles,
            std::move(primitive), std::move(denominator));
    }

    case BuiltinId::Expc:
        // expc(u)=(Exp(u)-1)/u。Ei(u)-Log(u)は0で個別に特異だが，局所原始函数としてDで一致する。
        return divideByDu(subtract(builtins, mathematics, angles,
            call(builtins, BuiltinId::ExponentialIntegralEi, {u}),
            call(builtins, BuiltinId::Log, {u})));

    case BuiltinId::ExponentialIntegralEi:
        // ∫Ei(u)du = u Ei(u)-Exp(u)。
        return divideByDu(subtract(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {u, expression}),
            call(builtins, BuiltinId::Exp, {u})));

    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi: {
        // Si/Ciの定義核はsession angle modeではなく常にRadian。
        Expr radianU = call(builtins, BuiltinId::UnitApplied, {
            u, Expr{std::string{"Rad"}}});
        Expr correction = call(builtins,
            definition->id == BuiltinId::SineIntegralSi ? BuiltinId::Cos : BuiltinId::Sin,
            {std::move(radianU)});
        Expr main = multiply(builtins, mathematics, angles, {u, expression});
        return divideByDu(definition->id == BuiltinId::SineIntegralSi
            ? add(builtins, mathematics, angles, {std::move(main), std::move(correction)})
            : subtract(builtins, mathematics, angles, std::move(main), std::move(correction)));
    }

    case BuiltinId::LogarithmicIntegralLi: {
        // principal li(u)=Ei(Log(u))に対するbranch-safeな局所原始函数。
        // li(u^2)のようなglobal branch rewriteは使わない。
        Expr doubledLog = multiply(builtins, mathematics, angles, {
            integer(2), call(builtins, BuiltinId::Log, {u})});
        return divideByDu(subtract(builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {u, expression}),
            call(builtins, BuiltinId::ExponentialIntegralEi, {std::move(doubledLog)})));
    }

    case BuiltinId::Digamma:
        return divideByDu(call(builtins, BuiltinId::LogGamma, {u}));

    case BuiltinId::Trigamma:
        return divideByDu(call(builtins, BuiltinId::Digamma, {u}));

    case BuiltinId::LambertW: {
        // DLMF 4.13.12と同値だが z/W(z) を使わず，W=0で偽のholeを作らない形。
        // ∫W(u)du = Exp(W(u)) (W(u)^2-W(u)+1)。
        Expr w = expression;
        Expr polynomial = add(builtins, mathematics, angles, {
            power(builtins, mathematics, angles, w, integer(2)),
            negate(builtins, mathematics, angles, w),
            one});
        return divideByDu(multiply(builtins, mathematics, angles, {
            call(builtins, BuiltinId::Exp, {w}), std::move(polynomial)}));
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

struct PowerLikeFactor final {
    Expr base;
    Rational exponent;
};

[[nodiscard]] std::optional<PowerLikeFactor> powerLikeFactor(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (isHead(expression, builtins, BuiltinId::Power)
        && expression.asCall().arguments.size() == 2) {
        const auto exponent = exactRealRational(expression.asCall().arguments[1]);
        if (exponent)
            return PowerLikeFactor{expression.asCall().arguments[0], *exponent};
    }
    if (isHead(expression, builtins, BuiltinId::Sqrt)
        && expression.asCall().arguments.size() == 1)
        return PowerLikeFactor{
            expression.asCall().arguments[0], Rational{BigInt{1}, BigInt{2}}};
    if (isHead(expression, builtins, BuiltinId::Cbrt)
        && expression.asCall().arguments.size() == 1)
        return PowerLikeFactor{
            expression.asCall().arguments[0], Rational{BigInt{1}, BigInt{3}}};
    return std::nullopt;
}

// f'(x) f(x)^p を式形の簡約能力に依存せず直接認識する。
// 旧reverse-chainは候補primitiveを一度Dして全体比例を調べるため、
// D[(2/3)u sqrt[u]] が複数項へ展開された場合に x sqrt[1+x^2] を見落とした。
[[nodiscard]] std::optional<Expr> tryPowerChainProduct(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Multiply)
        && !isHead(expression, builtins, BuiltinId::Divide))
        return std::nullopt;
    for (const Expr& subexpression : dependentSubexpressions(expression, variable, builtins)) {
        const auto factor = powerLikeFactor(subexpression, builtins);
        if (!factor || !containsVariable(factor->base, variable))
            continue;

        Expr derivative = differentiateExpression(
            factor->base, variable, builtins, mathematics, angles);
        Expr target = multiply(builtins, mathematics, angles, {
            std::move(derivative), subexpression});
        const auto ratio = proportionalFactor(
            expression, target, variable, builtins, mathematics, angles);
        if (!ratio)
            continue;

        if (factor->exponent == Rational{BigInt{-1}}) {
            Expr primitive = call(builtins, BuiltinId::Log, {factor->base});
            return multiply(builtins, mathematics, angles, {*ratio, std::move(primitive)});
        }
        const Rational next = factor->exponent + Rational{BigInt{1}};
        if (next.isZero())
            continue;
        Expr primitive = divide(builtins, mathematics, angles,
            power(builtins, mathematics, angles, factor->base, rational(next)),
            rational(next));
        return multiply(builtins, mathematics, angles, {*ratio, std::move(primitive)});
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> tryPrincipalSquareRootSubstitution(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t depth) {
    bool sawPrincipalRoot = false;
    auto transformed = rewriteForSquareRootSubstitution(
        expression, variable, builtins, mathematics, angles, sawPrincipalRoot);
    if (!transformed || !sawPrincipalRoot)
        return std::nullopt;

    // t=sqrt[x], x=t^2, dx=2t dt。rewrite helperでは作業変数に同じSymbolを再利用する。
    Expr transformedIntegrand = multiply(builtins, mathematics, angles, {
        integer(2), Expr{variable}, std::move(*transformed)});
    Expr primitive = integrateCore(
        transformedIntegrand, variable, builtins, mathematics, angles, depth + 1);
    if (containsHead(primitive, builtins, BuiltinId::SymbolicIntegral))
        return std::nullopt;

    Expr root = call(builtins, BuiltinId::Sqrt, {Expr{variable}});
    return simplify(
        substituteSymbol(primitive, variable, root), builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> tryDistributeProductOverSum(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t depth) {
    if (!isHead(expression, builtins, BuiltinId::Multiply))
        return std::nullopt;
    const auto& factors = expression.asCall().arguments;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        if (!isHead(factors[i], builtins, BuiltinId::Add)
            || factors[i].asCall().arguments.size() > 8)
            continue;

        std::vector<Expr> expandedTerms;
        expandedTerms.reserve(factors[i].asCall().arguments.size());
        for (const Expr& term : factors[i].asCall().arguments) {
            std::vector<Expr> product = factors;
            product[i] = term;
            expandedTerms.push_back(multiply(
                builtins, mathematics, angles, std::move(product)));
        }
        Expr expanded = add(builtins, mathematics, angles, std::move(expandedTerms));
        Expr primitive = integrateCore(
            expanded, variable, builtins, mathematics, angles, depth + 1);
        if (!isHead(primitive, builtins, BuiltinId::SymbolicIntegral))
            return primitive;
    }
    return std::nullopt;
}

struct RationalFunctionForm final {
    RationalPolynomial numerator;
    RationalPolynomial denominator;
};

[[nodiscard]] bool rationalFunctionBudgetOkay(const RationalFunctionForm& value) {
    return value.numerator.degree() <= 64 && value.denominator.degree() <= 64
        && value.numerator.coefficients().size() <= 128
        && value.denominator.coefficients().size() <= 128;
}

[[nodiscard]] std::optional<RationalFunctionForm> toRationalFunctionForm(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    const RationalPolynomial one{{Rational{BigInt{1}}}};
    if (expression.isNumber() && expression.asNumber().isReal())
        return RationalFunctionForm{
            RationalPolynomial{{expression.asNumber().asReal().toRational()}}, one};
    if (expression.isSymbol() && expression.asSymbol().sameIdentity(variable))
        return RationalFunctionForm{
            RationalPolynomial{{Rational{BigInt{0}}, Rational{BigInt{1}}}}, one};
    if (!expression.isCall())
        return std::nullopt;

    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;

    const auto combine = [&](RationalFunctionForm lhs,
                             const RationalFunctionForm& rhs,
                             BuiltinId operation) -> std::optional<RationalFunctionForm> {
        RationalFunctionForm result{one, one};
        if (operation == BuiltinId::Add || operation == BuiltinId::Subtract) {
            RationalPolynomial left = multiplyPolynomials(lhs.numerator, rhs.denominator);
            RationalPolynomial right = multiplyPolynomials(rhs.numerator, lhs.denominator);
            if (operation == BuiltinId::Subtract)
                right = negatePolynomial(right);
            result.numerator = addPolynomials(left, right);
            result.denominator = multiplyPolynomials(lhs.denominator, rhs.denominator);
        }
        else if (operation == BuiltinId::Multiply) {
            result.numerator = multiplyPolynomials(lhs.numerator, rhs.numerator);
            result.denominator = multiplyPolynomials(lhs.denominator, rhs.denominator);
        }
        else if (operation == BuiltinId::Divide) {
            if (rhs.numerator.isZero())
                return std::nullopt;
            result.numerator = multiplyPolynomials(lhs.numerator, rhs.denominator);
            result.denominator = multiplyPolynomials(lhs.denominator, rhs.numerator);
        }
        else {
            return std::nullopt;
        }
        if (!rationalFunctionBudgetOkay(result))
            return std::nullopt;
        return result;
    };

    switch (definition->id) {
    case BuiltinId::Negate:
        if (arguments.size() == 1) {
            auto inner = toRationalFunctionForm(arguments[0], variable, builtins);
            if (!inner)
                return std::nullopt;
            inner->numerator = negatePolynomial(inner->numerator);
            return inner;
        }
        return std::nullopt;

    case BuiltinId::Add:
    case BuiltinId::Multiply: {
        RationalFunctionForm accumulated = definition->id == BuiltinId::Add
            ? RationalFunctionForm{RationalPolynomial{}, one}
            : RationalFunctionForm{one, one};
        for (const Expr& argument : arguments) {
            auto part = toRationalFunctionForm(argument, variable, builtins);
            if (!part)
                return std::nullopt;
            auto next = combine(
                std::move(accumulated), *part, definition->id);
            if (!next)
                return std::nullopt;
            accumulated = std::move(*next);
        }
        return accumulated;
    }

    case BuiltinId::Subtract:
    case BuiltinId::Divide:
        if (arguments.size() != 2)
            return std::nullopt;
        if (auto lhs = toRationalFunctionForm(arguments[0], variable, builtins)) {
            if (auto rhs = toRationalFunctionForm(arguments[1], variable, builtins))
                return combine(std::move(*lhs), *rhs, definition->id);
        }
        return std::nullopt;

    case BuiltinId::Power:
        if (arguments.size() == 2) {
            auto base = toRationalFunctionForm(arguments[0], variable, builtins);
            const auto exponent = exactRealRational(arguments[1]);
            if (!base || !exponent || !exponent->isInteger())
                return std::nullopt;
            const BigInt& integerExponent = exponent->numerator();
            const auto magnitude = numeric::tryToUint64(
                integerExponent.isNegative() ? -integerExponent : integerExponent);
            if (!magnitude || *magnitude > 32)
                return std::nullopt;
            RationalFunctionForm result{
                powerPolynomial(base->numerator, *magnitude),
                powerPolynomial(base->denominator, *magnitude)};
            if (integerExponent.isNegative())
                std::swap(result.numerator, result.denominator);
            if (!rationalFunctionBudgetOkay(result) || result.denominator.isZero())
                return std::nullopt;
            return result;
        }
        return std::nullopt;

    default:
        return std::nullopt;
    }
}

[[nodiscard]] std::optional<Expr> normalizeRationalFunction(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto form = toRationalFunctionForm(expression, variable, builtins);
    if (!form)
        return std::nullopt;
    RationalPolynomial numeratorPolynomial = form->numerator;
    RationalPolynomial denominatorPolynomial = form->denominator;
    if (!numeratorPolynomial.isZero()) {
        const RationalPolynomial gcd = polynomialGcdMonic(
            numeratorPolynomial, denominatorPolynomial);
        if (gcd.degree() > 0) {
            PolynomialDivision numeratorDivision = dividePolynomials(numeratorPolynomial, gcd);
            PolynomialDivision denominatorDivision = dividePolynomials(denominatorPolynomial, gcd);
            if (numeratorDivision.remainder.isZero() && denominatorDivision.remainder.isZero()) {
                numeratorPolynomial = std::move(numeratorDivision.quotient);
                denominatorPolynomial = std::move(denominatorDivision.quotient);
            }
        }
    }
    // Q[x]上の有理函数は分母leading coefficientを1へ正規化する。
    // これにより 216x/(216x^3-216) のような共通scalarもexactに消える。
    if (!denominatorPolynomial.isZero()) {
        const Rational scale = denominatorPolynomial.coefficient(denominatorPolynomial.degree());
        if (!(scale == Rational{BigInt{1}})) {
            std::vector<Rational> numeratorCoefficients(
                numeratorPolynomial.coefficients().begin(), numeratorPolynomial.coefficients().end());
            std::vector<Rational> denominatorCoefficients(
                denominatorPolynomial.coefficients().begin(), denominatorPolynomial.coefficients().end());
            for (Rational& coefficient : numeratorCoefficients)
                coefficient /= scale;
            for (Rational& coefficient : denominatorCoefficients)
                coefficient /= scale;
            numeratorPolynomial = RationalPolynomial{std::move(numeratorCoefficients)};
            denominatorPolynomial = RationalPolynomial{std::move(denominatorCoefficients)};
        }
    }

    Expr numerator = polynomialToExpandedExpr(numeratorPolynomial, variable, builtins);
    Expr denominator = polynomialToExpandedExpr(denominatorPolynomial, variable, builtins);
    return simplify(
        divide(builtins, mathematics, angles, std::move(numerator), std::move(denominator)),
        builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> rewriteWeierstrassRational(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::optional<Expr>& commonRawArgument,
    std::optional<Expr>& commonScale,
    std::optional<mathematics::AngleUnit>& commonUnit,
    bool& sawTrig) {
    if (!containsVariable(expression, variable))
        return expression;
    if (expression.isSymbol())
        return std::nullopt;
    if (!expression.isCall())
        return std::nullopt;

    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return std::nullopt;
    if ((definition->id == BuiltinId::Sin || definition->id == BuiltinId::Cos)
        && expression.asCall().arguments.size() == 1) {
        const Expr& source = expression.asCall().arguments[0];
        TrigArgument info = trigArgument(source, builtins, mathematics, angles);
        const auto unit = explicitAngleUnit(source, builtins);
        if (commonRawArgument) {
            if (*commonRawArgument != info.argument || *commonScale != info.scale
                || commonUnit != unit)
                return std::nullopt;
        }
        else {
            commonRawArgument = info.argument;
            commonScale = info.scale;
            commonUnit = unit;
        }
        sawTrig = true;

        Expr t{variable};
        Expr t2 = power(builtins, mathematics, angles, t, integer(2));
        Expr denominator = add(builtins, mathematics, angles, {integer(1), t2});
        if (definition->id == BuiltinId::Sin)
            return divide(builtins, mathematics, angles,
                multiply(builtins, mathematics, angles, {integer(2), t}), denominator);
        return divide(builtins, mathematics, angles,
            subtract(builtins, mathematics, angles, integer(1), t2), denominator);
    }

    switch (definition->id) {
    case BuiltinId::Add:
    case BuiltinId::Subtract:
    case BuiltinId::Negate:
    case BuiltinId::Multiply:
    case BuiltinId::Divide:
    case BuiltinId::Power:
        break;
    default:
        return std::nullopt;
    }

    if (definition->id == BuiltinId::Power
        && expression.asCall().arguments.size() == 2
        && containsVariable(expression.asCall().arguments[1], variable))
        return std::nullopt;

    std::vector<Expr> arguments;
    arguments.reserve(expression.asCall().arguments.size());
    for (const Expr& argument : expression.asCall().arguments) {
        auto rewritten = rewriteWeierstrassRational(
            argument, variable, builtins, mathematics, angles,
            commonRawArgument, commonScale, commonUnit, sawTrig);
        if (!rewritten)
            return std::nullopt;
        arguments.push_back(std::move(*rewritten));
    }
    return simplify(Expr::rebuildCall(expression.asCall(), std::move(arguments)),
        builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> tryWeierstrassSubstitution(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t depth) {
    std::optional<Expr> rawArgument;
    std::optional<Expr> angleScale;
    std::optional<mathematics::AngleUnit> explicitUnit;
    bool sawTrig = false;
    auto rewritten = rewriteWeierstrassRational(
        expression, variable, builtins, mathematics, angles,
        rawArgument, angleScale, explicitUnit, sawTrig);
    if (!rewritten || !sawTrig || !rawArgument || !angleScale)
        return std::nullopt;

    Expr rawDerivative = differentiateExpression(
        *rawArgument, variable, builtins, mathematics, angles);
    if (containsVariable(rawDerivative, variable))
        return std::nullopt;
    Expr thetaDerivative = multiply(
        builtins, mathematics, angles, {*angleScale, std::move(rawDerivative)});
    if (!provablyNonZero(thetaDerivative, builtins, mathematics))
        return std::nullopt;

    Expr t{variable};
    Expr measure = divide(builtins, mathematics, angles, integer(2),
        multiply(builtins, mathematics, angles, {
            thetaDerivative,
            add(builtins, mathematics, angles, {
                integer(1), power(builtins, mathematics, angles, t, integer(2))})}));
    Expr rationalIntegrand = simplify(
        multiply(builtins, mathematics, angles, {std::move(*rewritten), std::move(measure)}),
        builtins, mathematics, angles);
    if (const auto normalized = normalizeRationalFunction(
            rationalIntegrand, variable, builtins, mathematics, angles))
        rationalIntegrand = *normalized;
    Expr primitive = integrateCore(
        rationalIntegrand, variable, builtins, mathematics, angles, depth + 1);
    if (containsHead(primitive, builtins, BuiltinId::SymbolicIntegral))
        return std::nullopt;

    Expr halfRaw = divide(
        builtins, mathematics, angles, *rawArgument, integer(2));
    Expr halfAngle = halfRaw;
    if (explicitUnit) {
        const char* unitName = *explicitUnit == mathematics::AngleUnit::Degree ? "Deg"
            : (*explicitUnit == mathematics::AngleUnit::Gradian ? "Grad" : "Rad");
        halfAngle = call(builtins, BuiltinId::UnitApplied, {
            std::move(halfRaw), Expr{std::string{unitName}}});
    }
    Expr tangent = call(builtins, BuiltinId::Tan, {std::move(halfAngle)});
    return simplify(
        substituteSymbol(primitive, variable, tangent), builtins, mathematics, angles);
}

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
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::IntegrationCandidate);
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
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::IntegrationCandidate);
    if (depth > maximumIntegrationDepth)
        return unresolved(original, variable, builtins);

    const Expr expression = simplify(original, builtins, mathematics, angles);
    if (!containsVariable(expression, variable))
        return multiply(builtins, mathematics, angles, {expression, Expr{variable}});

    if (expression.isSymbol() && expression.asSymbol().sameIdentity(variable))
        return divide(builtins, mathematics, angles,
            power(builtins, mathematics, angles, Expr{variable}, integer(2)), integer(2));

    // principal li(x) = Ei(Log(x)) の標準原始函数。
    // d/dx Ei(2 Log(x)) = Exp(2 Log(x))/(x Log(x)) = x/Log(x)
    // はprincipal Logの定義域で成立するため、
    // d/dx [x li(x) - Ei(2 Log(x))] = li(x)。
    // li(x^2) と書くと principal branch 上で Log(x^2)=2 Log(x) が
    // 大域的には成立しないため、branch-safeな Ei(2 Log(x)) 形を保持する。
    if (isHead(expression, builtins, BuiltinId::LogarithmicIntegralLi)
        && expression.asCall().arguments.size() == 1
        && expression.asCall().arguments[0].isSymbol()
        && expression.asCall().arguments[0].asSymbol().sameIdentity(variable)) {
        Expr x{variable};
        Expr logarithm = call(builtins, BuiltinId::Log, {x});
        Expr doubledLog = multiply(
            builtins, mathematics, angles, {integer(2), std::move(logarithm)});
        return subtract(
            builtins, mathematics, angles,
            multiply(builtins, mathematics, angles, {x, expression}),
            call(builtins, BuiltinId::ExponentialIntegralEi, {std::move(doubledLog)}));
    }

    if (isHead(expression, builtins, BuiltinId::Power)) {
        if (auto result = integrateReciprocalTrigPower(
                expression, variable, builtins, mathematics, angles))
            return *result;
        if (auto result = integrateDirectTrigPower(
                expression, variable, builtins, mathematics, angles))
            return *result;
        if (auto result = integratePowerRule(
                expression, variable, builtins, mathematics, angles))
            return *result;
        if (auto result = integrateGenericSymbolicPowerRule(
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

    if (auto powerChain = tryPowerChainProduct(
            expression, variable, builtins, mathematics, angles))
        return *powerChain;

    if (auto polynomial = integratePolynomial(
            expression, variable, builtins, mathematics, angles))
        return *polynomial;

    if (auto normalizedRational = normalizeRationalFunction(
            expression, variable, builtins, mathematics, angles);
        normalizedRational && *normalizedRational != expression) {
        if (auto rationalQuadratic = tryRationalLowDegree(
                *normalizedRational, variable, builtins, mathematics, angles))
            return *rationalQuadratic;
        if (auto rationalFactored = tryRationalFactored(
                *normalizedRational, variable, builtins, mathematics, angles))
            return *rationalFactored;
    }

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
    if (auto quadraticSqrt = tryQuadraticSquareRoot(
            expression, variable, builtins, mathematics, angles))
        return *quadraticSqrt;
    if (auto nestedRoot = tryNestedQuadraticSquareRoot(
            expression, variable, builtins, mathematics, angles))
        return *nestedRoot;

    if (auto sineCosineIntegral = integrateSineCosineOverArgument(
            expression, variable, builtins, mathematics, angles))
        return *sineCosineIntegral;

    if (expression.isCall()) {
        const auto* definition = builtins.find(expression.asCall().head);
        const auto& a = expression.asCall().arguments;
        if (definition) {
            switch (definition->id) {
            case BuiltinId::Cases: {
                std::vector<Expr> branches;
                branches.reserve(a.size());
                for (const Expr& branchExpression : a) {
                    if (!isHead(branchExpression, builtins, BuiltinId::CaseBranch)
                        || branchExpression.asCall().arguments.empty()
                        || branchExpression.asCall().arguments.size() > 2)
                        return unresolved(expression, variable, builtins);
                    const auto& branch = branchExpression.asCall().arguments;
                    // 条件が積分変数へ依存しない場合だけ線形にbranchへ分配する。
                    // x依存の境界を持つpiecewise積分では積分定数の整合や境界連続性が別問題になるため，
                    // 現段階では元のcasesを未解決のまま保持する。
                    if (branch.size() == 2 && containsVariable(branch[1], variable))
                        return unresolved(expression, variable, builtins);
                    Expr primitive = integrateCore(
                        branch[0], variable, builtins, mathematics, angles, depth + 1);
                    if (branch.size() == 2)
                        branches.push_back(caseBranch(builtins, std::move(primitive), branch[1]));
                    else
                        branches.push_back(caseBranch(builtins, std::move(primitive)));
                }
                return cases(builtins, std::move(branches));
            }
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
                // chain-ruleで一発に閉じる式を積への分配より先に試す。
                // 例えば (3x^2+1)exp[x^3+x] は全体でexp[x^3+x]の微分だが、
                // 先に分配すると各項が未評価になり、既に持っていた積分能力を失う。
                if (auto reverse = tryReverseChainRule(
                        expression, variable, builtins, mathematics, angles))
                    return *reverse;
                if (auto distributed = tryDistributeProductOverSum(
                        expression, variable, builtins, mathematics, angles, depth))
                    return *distributed;
                // Fresnel等の特殊函数primitiveが増えるとpartsが先に成功して式を肥大化し得る。
                // exact chain-ruleとbounded distributionの双方をpartsより先に採る。
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
                        Expr logArgument = a[1];
                        if (const auto affine = toRationalPolynomial(
                                a[1], variable, builtins, PolynomialConversionOptions{2, 16});
                            affine && affine->degree() == 1)
                            logArgument = polynomialToExpandedExpr(*affine, variable, builtins);
                        return multiply(builtins, mathematics, angles, {
                            a[0],
                            divide(builtins, mathematics, angles,
                                call(builtins, BuiltinId::Log, {std::move(logArgument)}),
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
                if (auto result = integrateGenericSymbolicPowerRule(
                        expression, variable, builtins, mathematics, angles))
                    return *result;
                break;
            default:
                break;
            }
        }
    }

    if (auto ellipticTrig = integrateEllipticTrigKernel(
            expression, variable, builtins, mathematics, angles))
        return *ellipticTrig;

    if (auto fresnel = integrateQuadraticFresnel(
            expression, variable, builtins, mathematics, angles))
        return *fresnel;

    if (auto elliptic = integrateQuarticEllipticF(
            expression, variable, builtins, mathematics, angles))
        return *elliptic;

    if (auto hypergeometric = integrateExponentialMonomial1F1(
            expression, variable, builtins, mathematics, angles))
        return *hypergeometric;

    if (auto dilogarithm = integrateLogOnePlusMonomialOverX(
            expression, variable, builtins, mathematics, angles))
        return *dilogarithm;

    if (auto elementary = integrateElementaryUnary(
            expression, variable, builtins, mathematics, angles))
        return *elementary;

    if (auto standard = integrateStandardUnary(
            expression, variable, builtins, mathematics, angles))
        return *standard;

    if (auto incompleteBeta = integrateIncompleteBetaAffine(
            expression, variable, builtins, mathematics, angles))
        return *incompleteBeta;

    if (auto lambertLogDerivative = integrateLambertWLogDerivative(
            expression, variable, builtins, mathematics, angles))
        return *lambertLogDerivative;

    if (auto shifted = trySpecialFunctionShiftPrimitive(
            expression, variable, builtins, mathematics, angles))
        return *shifted;

    // 同一argumentの高次sin/cos整数冪は有限Fourier和へ厳密還元できる。
    // 汎用reverse-chain / Weierstrass / parts探索より先に処理し、
    // 明らかに有限なfamilyで候補探索コストを払わない。平方とsin*cosは既存のcompact解を優先する。
    if (auto expansion = mathematics::expandTrigMonomial(expression, builtins, 256);
        expansion && expansion->totalDegree >= 3) {
        if (auto result = integrateTrigFourierExpansion(
                std::move(*expansion), variable, builtins, mathematics, angles))
            return *result;
    }

    if (auto reverse = tryReverseChainRule(
            expression, variable, builtins, mathematics, angles))
        return *reverse;

    if (auto rootSubstitution = tryPrincipalSquareRootSubstitution(
            expression, variable, builtins, mathematics, angles, depth))
        return *rootSubstitution;

    // R(sin(theta),cos(theta)) を t=tan(theta/2) で有理函数へ落とす。
    // 局所的な変数変換なのでtanのpoleを跨ぐ大域同値は主張しない。
    if (auto weierstrass = tryWeierstrassSubstitution(
            expression, variable, builtins, mathematics, angles, depth))
        return *weierstrass;

    if (auto byParts = tryIntegrationByParts(
            expression, variable, builtins, mathematics, angles, depth))
        return *byParts;

    // 旧実装ではsin/cosの二乗だけをrewriteSquareIdentityで個別処理していたため、
    // sin[2x]^6のような有限Fourier多項式へ厳密還元できる整数冪が未評価で残った。
    // 既存のcompactなreverse-chain/parts解を優先した後のfallbackとして共有TrigKnowledgeを使い、
    // sin[u]^m cos[u]^nを有限Fourier和へ落として既存の線形積分規則へ再投入する。
    // 積分器は大きな出力を明示的に要求された場合に限り256次まで展開する。
    // FullSimplify側の既定64次上限は維持し、通常の簡約で128項級の式爆発を起こさない。
    if (auto expansion = mathematics::expandTrigMonomial(expression, builtins, 256)) {
        if (auto result = integrateTrigFourierExpansion(
                std::move(*expansion), variable, builtins, mathematics, angles))
            return *result;

        // 非線形argument等で共通duを定数として抜けない場合は，従来どおり
        // Fourier式を汎用integratorへ戻す。Fresnel等への後続変換能力を失わない。
        if (const auto reduced = mathematics::reduceTrigMonomial(expression, builtins, 256)) {
            Expr result = integrateCore(
                simplify(*reduced, builtins, mathematics, angles),
                variable, builtins, mathematics, angles, depth + 1);
            if (!isHead(result, builtins, BuiltinId::SymbolicIntegral))
                return result;
        }
    }

    // 異なる引数の一次sin/cos積は上のmonomial規則では扱えない。
    // 積和恒等式は複素引数でも大域的に成立するため、安全な最後の代数fallbackとして使う。
    if (const auto reduced = mathematics::reduceTrigProduct(expression, builtins)) {
        Expr result = integrateCore(
            fullSimplify(*reduced, builtins, mathematics, angles),
            variable, builtins, mathematics, angles, depth + 1);
        if (!isHead(result, builtins, BuiltinId::SymbolicIntegral))
            return result;
    }

    // 一般binomial 2F1はexactだが，atan/asinや単純reverse-chainより表現が重い。
    // したがって特殊化されたelementary substitutionをすべて試した後のfallbackに置く。
    // algebraic Root/Log全分解よりはcompactなので，該当するbinomial familyでは2F1を優先する。
    if (auto hypergeometric2F1 = integrateBinomialPower2F1(
            expression, variable, builtins, mathematics, angles))
        return *hypergeometric2F1;

    // 高次有理因子のHermite reduction / algebraic-log分解はexactだが，
    // elementary/2F1より表現が重い最終fallbackである。通常の有理部分分数では
    // 一次・二次因子だけを先に処理し，次数3以上はここまで閉じなかった場合だけ展開する。
    if (auto normalizedRational = normalizeRationalFunction(
            expression, variable, builtins, mathematics, angles);
        normalizedRational) {
        if (auto rationalAlgebraic = tryRationalFactored(
                *normalizedRational, variable, builtins, mathematics, angles, true))
            return *rationalAlgebraic;
    }
    if (auto rationalAlgebraic = tryRationalFactored(
            expression, variable, builtins, mathematics, angles, true))
        return *rationalAlgebraic;

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
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
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
    catch (const error::CalcError&) {
        // Unified EvaluationBudgetの超過を「intervalで証明できなかった」へ
        // 読み替えない。資源diagnosticはtop-levelまでそのまま伝播させる。
        throw;
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

[[nodiscard]] bool nestedLogarithmEndpointClass(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    if (!isHead(expression, builtins, BuiltinId::Log)
        || expression.asCall().arguments.size() != 1)
        return false;
    const Expr& innerLog = expression.asCall().arguments[0];
    if (!isHead(innerLog, builtins, BuiltinId::Log)
        || innerLog.asCall().arguments.size() != 1)
        return false;
    const Expr& inner = innerLog.asCall().arguments[0];
    if (!inner.isSymbol() || !inner.asSymbol().sameIdentity(variable))
        return false;

    const auto lo = exactRealRational(lower);
    if (!lo || *lo != Rational{BigInt{1}})
        return false;
    if (const auto hi = exactRealRational(upper))
        return *hi > Rational{BigInt{1}};
    if (!upper.isSymbol())
        return false;
    const auto* constant = mathematics.findConstant(upper.asSymbol());
    return constant && (constant->id == mathematics::ConstantId::E
        || constant->id == mathematics::ConstantId::Pi
        || constant->id == mathematics::ConstantId::Phi);
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
    if (logarithmEndpointClass(expression, variable, lower, upper, builtins))
        return true;
    return nestedLogarithmEndpointClass(
        expression, variable, lower, upper, builtins, mathematics);
}

[[nodiscard]] bool isNestedLogOfVariable(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Log)
        || expression.asCall().arguments.size() != 1)
        return false;
    const Expr& inner = expression.asCall().arguments[0];
    return isHead(inner, builtins, BuiltinId::Log)
        && inner.asCall().arguments.size() == 1
        && inner.asCall().arguments[0].isSymbol()
        && inner.asCall().arguments[0].asSymbol().sameIdentity(variable);
}

[[nodiscard]] bool isVariableTimesNestedLog(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Multiply)
        || expression.asCall().arguments.size() != 2)
        return false;
    const auto& factors = expression.asCall().arguments;
    return (factors[0].isSymbol() && factors[0].asSymbol().sameIdentity(variable)
            && isNestedLogOfVariable(factors[1], variable, builtins))
        || (factors[1].isSymbol() && factors[1].asSymbol().sameIdentity(variable)
            && isNestedLogOfVariable(factors[0], variable, builtins));
}

[[nodiscard]] bool isLiOfVariable(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    return isHead(expression, builtins, BuiltinId::LogarithmicIntegralLi)
        && expression.asCall().arguments.size() == 1
        && expression.asCall().arguments[0].isSymbol()
        && expression.asCall().arguments[0].asSymbol().sameIdentity(variable);
}

[[nodiscard]] bool isLogLogPrimitiveAtUnitEndpoint(
    const Expr& primitive,
    const expression::Symbol& variable,
    const Expr& endpoint,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins) {
    const auto point = exactRealRational(endpoint);
    if (!point || *point != Rational{BigInt{1}} || direction != LimitDirection::Right)
        return false;
    if (!isHead(primitive, builtins, BuiltinId::Subtract)
        || primitive.asCall().arguments.size() != 2)
        return false;
    const auto& terms = primitive.asCall().arguments;
    return isVariableTimesNestedLog(terms[0], variable, builtins)
        && isLiOfVariable(terms[1], variable, builtins);
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
    // lim_{x->1+}(x log(log(x)) - li(x)) = -EulerGamma = digamma(1).
    // 各項は個別に発散するため、通常の項別endpoint評価より先に組として処理する。
    if (isLogLogPrimitiveAtUnitEndpoint(primitive, variable, endpoint, direction, builtins))
        return call(builtins, BuiltinId::Digamma, {integer(1)});

    if (const auto endpointValue = exactRealRational(endpoint);
        endpointValue && endpointValue->isZero()
        && (isHead(primitive, builtins, BuiltinId::SineIntegralSi)
            || isHead(primitive, builtins, BuiltinId::FresnelC)
            || isHead(primitive, builtins, BuiltinId::FresnelS))
        && primitive.asCall().arguments.size() == 1
        && primitive.asCall().arguments[0].isSymbol()
        && primitive.asCall().arguments[0].asSymbol().sameIdentity(variable))
        return integer(0);

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

namespace {

[[nodiscard]] bool isRecognizedNoFiniteClosedFormFamily(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    // ここは「一般に積分不能」を証明する判定器ではない。
    // カタログで意図的に未評価とする代表familyだけを明示認識し、
    // 実装漏れによるUnsupportedとユーザー表示上で区別する。
    if (isHead(expression, builtins, BuiltinId::Power)
        && expression.asCall().arguments.size() == 2
        && expression.asCall().arguments[0].isSymbol()
        && expression.asCall().arguments[0].asSymbol().sameIdentity(variable)
        && expression.asCall().arguments[1].isSymbol()
        && expression.asCall().arguments[1].asSymbol().sameIdentity(variable))
        return true;

    if (isHead(expression, builtins, BuiltinId::Sin)
        && expression.asCall().arguments.size() == 1) {
        const Expr& inner = expression.asCall().arguments[0];
        if (isHead(inner, builtins, BuiltinId::Sin)
            && inner.asCall().arguments.size() == 1
            && inner.asCall().arguments[0].isSymbol()
            && inner.asCall().arguments[0].asSymbol().sameIdentity(variable))
            return true;
    }

    if (isHead(expression, builtins, BuiltinId::Exp)
        && expression.asCall().arguments.size() == 1) {
        const Expr& inner = expression.asCall().arguments[0];
        if (isHead(inner, builtins, BuiltinId::Sin)
            && inner.asCall().arguments.size() == 1
            && inner.asCall().arguments[0].isSymbol()
            && inner.asCall().arguments[0].asSymbol().sameIdentity(variable))
            return true;
    }
    return false;
}


[[nodiscard]] bool isExactZeroExpr(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] bool isExactOneExpr(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isReal()
        && expression.asNumber().asReal() == numeric::RealNumber{BigInt{1}};
}

[[nodiscard]] bool proveGreaterThan(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    return mathematics::KnowledgeContext{builtins, mathematics, assumptions}.prove(
        mathematics::relation(mathematics::RelationKind::Greater, lhs, rhs))
        == mathematics::TruthValue::True;
}

[[nodiscard]] std::optional<Expr> variablePowerExponent(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (expression.isSymbol() && expression.asSymbol().sameIdentity(variable))
        return integer(1);
    if (isHead(expression, builtins, BuiltinId::Power)
        && expression.asCall().arguments.size() == 2
        && expression.asCall().arguments[0].isSymbol()
        && expression.asCall().arguments[0].asSymbol().sameIdentity(variable)
        && !containsVariable(expression.asCall().arguments[1], variable))
        return expression.asCall().arguments[1];
    return std::nullopt;
}

struct GammaKernelMatch final {
    Expr outside;
    Expr nu;
    Rational mu;
    Expr scale;
};

[[nodiscard]] std::optional<GammaKernelMatch> matchPositiveHalfLineGammaKernel(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    FactorSplit outer = splitConstantFactor(
        expression, variable, builtins, mathematics, angles);

    std::vector<Expr> factors;
    if (isHead(outer.dependent, builtins, BuiltinId::Multiply))
        factors.assign(
            outer.dependent.asCall().arguments.begin(),
            outer.dependent.asCall().arguments.end());
    else
        factors.push_back(outer.dependent);

    const Expr* exponential = nullptr;
    std::optional<Expr> monomialExponent;
    for (const Expr& factor : factors) {
        if (isHead(factor, builtins, BuiltinId::Exp)
            && factor.asCall().arguments.size() == 1) {
            if (exponential)
                return std::nullopt;
            exponential = &factor;
            continue;
        }
        if (isOne(factor))
            continue;
        if (monomialExponent)
            return std::nullopt;
        monomialExponent = variablePowerExponent(factor, variable, builtins);
        if (!monomialExponent)
            return std::nullopt;
    }
    if (!exponential)
        return std::nullopt;

    const Expr& exponent = exponential->asCall().arguments[0];
    FactorSplit exponentSplit = splitConstantFactor(
        exponent, variable, builtins, mathematics, angles);
    auto muExpr = variablePowerExponent(exponentSplit.dependent, variable, builtins);
    if (!muExpr)
        return std::nullopt;
    const auto mu = exactRealRational(*muExpr);
    if (!mu || !(*mu > Rational{BigInt{0}}))
        return std::nullopt;

    Expr scale = simplify(
        negate(builtins, mathematics, angles, exponentSplit.constant),
        builtins, mathematics, angles);
    Expr nu = simplify(
        add(builtins, mathematics, angles,
            {monomialExponent ? *monomialExponent : integer(0), integer(1)}),
        builtins, mathematics, angles);
    return GammaKernelMatch{std::move(outer.constant), std::move(nu), *mu, std::move(scale)};
}

[[nodiscard]] std::optional<std::pair<Expr, Expr>> matchEulerBetaKernel(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Multiply)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const Expr oneMinusX = simplify(
        subtract(builtins, mathematics, angles, integer(1), Expr{variable}),
        builtins, mathematics, angles);
    std::optional<Expr> a;
    std::optional<Expr> b;
    for (const Expr& factor : expression.asCall().arguments) {
        if (!isHead(factor, builtins, BuiltinId::Power)
            || factor.asCall().arguments.size() != 2
            || containsVariable(factor.asCall().arguments[1], variable))
            return std::nullopt;
        const auto& args = factor.asCall().arguments;
        Expr parameter = simplify(
            add(builtins, mathematics, angles, {args[1], integer(1)}),
            builtins, mathematics, angles);
        if (args[0].isSymbol() && args[0].asSymbol().sameIdentity(variable))
            a = std::move(parameter);
        else if (simplify(args[0], builtins, mathematics, angles) == oneMinusX)
            b = std::move(parameter);
        else
            return std::nullopt;
    }
    if (!a || !b)
        return std::nullopt;
    return std::pair<Expr, Expr>{std::move(*a), std::move(*b)};
}

[[nodiscard]] bool isVariableSquareExpr(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Power)
        || expression.asCall().arguments.size() != 2)
        return false;
    const auto& args = expression.asCall().arguments;
    return args[0].isSymbol() && args[0].asSymbol().sameIdentity(variable)
        && args[1].isNumber() && args[1].asNumber().isReal()
        && args[1].asNumber().asReal() == numeric::RealNumber{BigInt{2}};
}

[[nodiscard]] std::optional<Expr> squaredParameter(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Power)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& args = expression.asCall().arguments;
    if (containsVariable(args[0], variable)
        || !args[1].isNumber() || !args[1].asNumber().isReal()
        || !(args[1].asNumber().asReal() == numeric::RealNumber{BigInt{2}}))
        return std::nullopt;
    return args[0];
}

[[nodiscard]] std::optional<Expr> symmetricCircleParameter(
    const Expr& radicand,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(radicand, builtins, BuiltinId::Subtract)
        || radicand.asCall().arguments.size() != 2
        || !isVariableSquareExpr(radicand.asCall().arguments[1], variable, builtins))
        return std::nullopt;
    return squaredParameter(radicand.asCall().arguments[0], variable, builtins);
}

[[nodiscard]] std::optional<Expr> cauchyScaleParameter(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2
        || !isExactOneExpr(expression.asCall().arguments[0]))
        return std::nullopt;
    const Expr& denominator = expression.asCall().arguments[1];
    if (!isHead(denominator, builtins, BuiltinId::Add)
        || denominator.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& terms = denominator.asCall().arguments;
    if (isVariableSquareExpr(terms[0], variable, builtins))
        return squaredParameter(terms[1], variable, builtins);
    if (isVariableSquareExpr(terms[1], variable, builtins))
        return squaredParameter(terms[0], variable, builtins);
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> zetaMellinParameter(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& args = expression.asCall().arguments;
    const auto exponent = variablePowerExponent(args[0], variable, builtins);
    if (!exponent)
        return std::nullopt;
    if (!isHead(args[1], builtins, BuiltinId::Subtract)
        || args[1].asCall().arguments.size() != 2
        || !isExactOneExpr(args[1].asCall().arguments[1]))
        return std::nullopt;
    const Expr& exponential = args[1].asCall().arguments[0];
    if (!isHead(exponential, builtins, BuiltinId::Exp)
        || exponential.asCall().arguments.size() != 1
        || !(exponential.asCall().arguments[0].isSymbol()
            && exponential.asCall().arguments[0].asSymbol().sameIdentity(variable)))
        return std::nullopt;
    return simplify(add(builtins, mathematics, angles, {*exponent, integer(1)}),
        builtins, mathematics, angles);
}

struct BetaReflectionKernel final {
    Expr a;
    Rational q;
};

[[nodiscard]] std::optional<BetaReflectionKernel> betaReflectionKernel(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& args = expression.asCall().arguments;

    Expr numeratorExponent = integer(0);
    if (!isExactOneExpr(args[0])) {
        auto exponent = variablePowerExponent(args[0], variable, builtins);
        if (!exponent)
            return std::nullopt;
        numeratorExponent = std::move(*exponent);
    }

    if (!isHead(args[1], builtins, BuiltinId::Add)
        || args[1].asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& denominatorTerms = args[1].asCall().arguments;
    const Expr* variablePower = nullptr;
    if (isExactOneExpr(denominatorTerms[0]))
        variablePower = &denominatorTerms[1];
    else if (isExactOneExpr(denominatorTerms[1]))
        variablePower = &denominatorTerms[0];
    if (!variablePower)
        return std::nullopt;

    auto qExpr = variablePowerExponent(*variablePower, variable, builtins);
    if (!qExpr)
        return std::nullopt;
    const auto q = exactRealRational(*qExpr);
    if (!q || !(*q > Rational{BigInt{0}}))
        return std::nullopt;

    Expr a = simplify(add(builtins, mathematics, angles, {
        std::move(numeratorExponent), integer(1)}), builtins, mathematics, angles);
    return BetaReflectionKernel{std::move(a), *q};
}

struct UnitPoleKernel final {
    Expr pole;
    bool negated = false;
};

[[nodiscard]] std::optional<UnitPoleKernel> unitPoleKernel(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2
        || !isExactOneExpr(expression.asCall().arguments[0]))
        return std::nullopt;
    const Expr& denominator = expression.asCall().arguments[1];
    if (!isHead(denominator, builtins, BuiltinId::Subtract)
        || denominator.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& terms = denominator.asCall().arguments;
    if (terms[0].isSymbol() && terms[0].asSymbol().sameIdentity(variable)
        && !containsVariable(terms[1], variable))
        return UnitPoleKernel{terms[1], false};
    if (terms[1].isSymbol() && terms[1].asSymbol().sameIdentity(variable)
        && !containsVariable(terms[0], variable))
        return UnitPoleKernel{terms[0], true};
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> fermiMellinParameter(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& args = expression.asCall().arguments;
    const auto exponent = variablePowerExponent(args[0], variable, builtins);
    if (!exponent)
        return std::nullopt;
    if (!isHead(args[1], builtins, BuiltinId::Add)
        || args[1].asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& denominatorTerms = args[1].asCall().arguments;
    const Expr* exponential = nullptr;
    if (isExactOneExpr(denominatorTerms[0]))
        exponential = &denominatorTerms[1];
    else if (isExactOneExpr(denominatorTerms[1]))
        exponential = &denominatorTerms[0];
    if (!exponential || !isHead(*exponential, builtins, BuiltinId::Exp)
        || exponential->asCall().arguments.size() != 1
        || !(exponential->asCall().arguments[0].isSymbol()
            && exponential->asCall().arguments[0].asSymbol().sameIdentity(variable)))
        return std::nullopt;
    return simplify(add(builtins, mathematics, angles, {*exponent, integer(1)}),
        builtins, mathematics, angles);
}

[[nodiscard]] std::optional<std::pair<Expr, Expr>> frullaniScales(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2
        || !(expression.asCall().arguments[1].isSymbol()
            && expression.asCall().arguments[1].asSymbol().sameIdentity(variable)))
        return std::nullopt;
    const Expr& numerator = expression.asCall().arguments[0];
    if (!isHead(numerator, builtins, BuiltinId::Subtract)
        || numerator.asCall().arguments.size() != 2)
        return std::nullopt;

    const auto scaleOf = [&](const Expr& term) -> std::optional<Expr> {
        if (!isHead(term, builtins, BuiltinId::Exp)
            || term.asCall().arguments.size() != 1)
            return std::nullopt;
        FactorSplit split = splitConstantFactor(
            term.asCall().arguments[0], variable, builtins, mathematics, angles);
        if (!(split.dependent.isSymbol()
            && split.dependent.asSymbol().sameIdentity(variable)))
            return std::nullopt;
        return simplify(negate(builtins, mathematics, angles, split.constant),
            builtins, mathematics, angles);
    };
    auto first = scaleOf(numerator.asCall().arguments[0]);
    auto second = scaleOf(numerator.asCall().arguments[1]);
    if (!first || !second)
        return std::nullopt;
    return std::pair<Expr, Expr>{std::move(*first), std::move(*second)};
}

struct LogMomentKernel final {
    Expr exponent;
    std::uint64_t logPower = 0;
};

[[nodiscard]] std::optional<LogMomentKernel> logMomentKernel(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    auto logPowerOf = [&](const Expr& factor) -> std::optional<std::uint64_t> {
        if (isHead(factor, builtins, BuiltinId::Log)
            && factor.asCall().arguments.size() == 1
            && factor.asCall().arguments[0].isSymbol()
            && factor.asCall().arguments[0].asSymbol().sameIdentity(variable))
            return 1;
        if (!isHead(factor, builtins, BuiltinId::Power)
            || factor.asCall().arguments.size() != 2
            || !isHead(factor.asCall().arguments[0], builtins, BuiltinId::Log)
            || factor.asCall().arguments[0].asCall().arguments.size() != 1
            || !factor.asCall().arguments[0].asCall().arguments[0].isSymbol()
            || !factor.asCall().arguments[0].asCall().arguments[0].asSymbol().sameIdentity(variable))
            return std::nullopt;
        const auto power = exactRealRational(factor.asCall().arguments[1]);
        if (!power || !power->isInteger() || !(*power >= Rational{BigInt{1}}))
            return std::nullopt;
        return numeric::tryToUint64(power->numerator());
    };

    if (auto m = logPowerOf(expression))
        return LogMomentKernel{integer(0), *m};
    if (!isHead(expression, builtins, BuiltinId::Multiply)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const auto& factors = expression.asCall().arguments;
    for (std::size_t logIndex = 0; logIndex < 2; ++logIndex) {
        auto m = logPowerOf(factors[logIndex]);
        if (!m)
            continue;
        auto exponent = variablePowerExponent(factors[1 - logIndex], variable, builtins);
        if (exponent)
            return LogMomentKernel{std::move(*exponent), *m};
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Rational> reciprocalSqrtOneMinusPower(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2
        || !isExactOneExpr(expression.asCall().arguments[0]))
        return std::nullopt;
    const Expr& denominator = expression.asCall().arguments[1];
    if (!isHead(denominator, builtins, BuiltinId::Sqrt)
        || denominator.asCall().arguments.size() != 1)
        return std::nullopt;
    const Expr& radicand = denominator.asCall().arguments[0];
    if (!isHead(radicand, builtins, BuiltinId::Subtract)
        || radicand.asCall().arguments.size() != 2
        || !isExactOneExpr(radicand.asCall().arguments[0]))
        return std::nullopt;
    auto qExpr = variablePowerExponent(radicand.asCall().arguments[1], variable, builtins);
    if (!qExpr)
        return std::nullopt;
    auto q = exactRealRational(*qExpr);
    if (!q || !(*q > Rational{BigInt{0}}))
        return std::nullopt;
    return q;
}

[[nodiscard]] Expr gammaExactKernelValue(
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    if (const auto value = exactRealRational(argument)) {
        if (*value == Rational{BigInt{1}})
            return integer(1);
        if (*value == Rational{BigInt{1}, BigInt{2}})
            return call(builtins, BuiltinId::Sqrt, {pi(mathematics)});
    }
    return call(builtins, BuiltinId::Gamma, {argument});
}

[[nodiscard]] std::optional<Expr> tryAssumptionAwareDefiniteIntegral(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& lower,
    const Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions) {
    const simplification::SimplificationContext context{
        builtins, mathematics, angles, assumptions};
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};

    const bool zeroToPositiveInfinity = isExactZeroExpr(lower)
        && isInfinityExpr(upper, infinitySymbol);
    const bool wholeReal = isNegativeInfinityExpr(lower, builtins, infinitySymbol)
        && isInfinityExpr(upper, infinitySymbol);

    // DLMF 25.5.1: Mellin representation of zeta, Re(s)>1.
    if (zeroToPositiveInfinity) {
        if (auto s = zetaMellinParameter(
                expression, variable, builtins, mathematics, angles)) {
            if (proveGreaterThan(*s, integer(1), builtins, mathematics, assumptions))
                return simplification::Simplifier{}.simplify(
                    multiply(builtins, mathematics, angles, {
                        gammaExactKernelValue(*s, builtins, mathematics),
                        call(builtins, BuiltinId::Zeta, {*s})}), context);
        }

        // DLMF 25.5.3.  The product (1-2^(1-s)) zeta[s] has a removable
        // singularity at s=1, so keep that point as an explicit exact branch.
        if (auto s = fermiMellinParameter(
                expression, variable, builtins, mathematics, angles)) {
            if (proveGreaterThan(*s, integer(0), builtins, mathematics, assumptions)) {
                Expr generic = multiply(builtins, mathematics, angles, {
                    gammaExactKernelValue(*s, builtins, mathematics),
                    subtract(builtins, mathematics, angles, integer(1),
                        power(builtins, mathematics, angles, integer(2),
                            subtract(builtins, mathematics, angles, integer(1), *s))),
                    call(builtins, BuiltinId::Zeta, {*s})});
                return call(builtins, BuiltinId::Cases, {
                    call(builtins, BuiltinId::CaseBranch, {
                        call(builtins, BuiltinId::Log, {integer(2)}),
                        call(builtins, BuiltinId::Equal, {*s, integer(1)})}),
                    call(builtins, BuiltinId::CaseBranch, {
                        std::move(generic),
                        call(builtins, BuiltinId::NotEqual, {*s, integer(1)})})});
            }
        }

        // Euler beta integral after t=x^q:
        // int_0^inf x^(a-1)/(1+x^q) dx = Pi/(q sin(Pi a/q)), 0<a<q.
        if (auto kernel = betaReflectionKernel(
                expression, variable, builtins, mathematics, angles)) {
            const Expr q = rational(kernel->q);
            if (proveGreaterThan(kernel->a, integer(0), builtins, mathematics, assumptions)
                && proveGreaterThan(q, kernel->a, builtins, mathematics, assumptions)) {
                Expr ratio = simplify(
                    divide(builtins, mathematics, angles, kernel->a, q),
                    builtins, mathematics, angles);
                Expr phase = call(builtins, BuiltinId::UnitApplied, {
                    multiply(builtins, mathematics, angles, {pi(mathematics), ratio}),
                    Expr{std::string{"Rad"}}});
                return divide(builtins, mathematics, angles, pi(mathematics),
                    multiply(builtins, mathematics, angles, {
                        q, call(builtins, BuiltinId::Sin, {std::move(phase)})}));
            }
        }

        // Frullani: int_0^inf (exp(-a x)-exp(-b x))/x dx = log(b/a), a,b>0.
        if (auto scales = frullaniScales(
                expression, variable, builtins, mathematics, angles)) {
            if (proveGreaterThan(scales->first, integer(0), builtins, mathematics, assumptions)
                && proveGreaterThan(scales->second, integer(0), builtins, mathematics, assumptions))
                return call(builtins, BuiltinId::Log, {
                    divide(builtins, mathematics, angles, scales->second, scales->first)});
        }
    }

    // A single symbolic pole outside a finite real interval is safe.  Using the
    // positive endpoint ratio avoids principal-Log branch constants on the negative side.
    if (!isInfinityExpr(lower, infinitySymbol)
        && !isInfinityExpr(upper, infinitySymbol)
        && !isNegativeInfinityExpr(lower, builtins, infinitySymbol)
        && proveGreaterThan(upper, lower, builtins, mathematics, assumptions)) {
        if (auto pole = unitPoleKernel(expression, variable, builtins)) {
            const bool below = proveGreaterThan(
                lower, pole->pole, builtins, mathematics, assumptions);
            const bool above = proveGreaterThan(
                pole->pole, upper, builtins, mathematics, assumptions);
            if (below || above) {
                Expr ratio = divide(builtins, mathematics, angles,
                    subtract(builtins, mathematics, angles, upper, pole->pole),
                    subtract(builtins, mathematics, angles, lower, pole->pole));
                Expr result = call(builtins, BuiltinId::Log, {
                    simplify(std::move(ratio), builtins, mathematics, angles)});
                if (pole->negated)
                    result = negate(builtins, mathematics, angles, std::move(result));
                return simplification::Simplifier{}.simplify(std::move(result), context);
            }
        }
    }

    // Cauchy kernel on the whole real axis.
    if (wholeReal) {
        if (auto scale = cauchyScaleParameter(expression, variable, builtins)) {
            if (proveGreaterThan(*scale, integer(0), builtins, mathematics, assumptions))
                return divide(builtins, mathematics, angles, pi(mathematics), *scale);
        }
    }

    // Semicircle / arcsine kernels on [-a,a], a>0.  These are endpoint-improper
    // but convergent; matching the exact radicand avoids relying on generic endpoint safety.
    if (isHead(lower, builtins, BuiltinId::Negate)
        && lower.asCall().arguments.size() == 1
        && upper == lower.asCall().arguments[0]
        && proveGreaterThan(upper, integer(0), builtins, mathematics, assumptions)) {
        const Expr* radicand = nullptr;
        bool reciprocalSqrt = false;
        if (isHead(expression, builtins, BuiltinId::Sqrt)
            && expression.asCall().arguments.size() == 1)
            radicand = &expression.asCall().arguments[0];
        else if (isHead(expression, builtins, BuiltinId::Divide)
            && expression.asCall().arguments.size() == 2
            && isExactOneExpr(expression.asCall().arguments[0])
            && isHead(expression.asCall().arguments[1], builtins, BuiltinId::Sqrt)
            && expression.asCall().arguments[1].asCall().arguments.size() == 1) {
            reciprocalSqrt = true;
            radicand = &expression.asCall().arguments[1].asCall().arguments[0];
        }
        if (radicand) {
            if (auto scale = symmetricCircleParameter(*radicand, variable, builtins);
                scale && *scale == upper) {
                if (reciprocalSqrt)
                    return pi(mathematics);
                return divide(builtins, mathematics, angles,
                    multiply(builtins, mathematics, angles, {
                        pi(mathematics), power(builtins, mathematics, angles, *scale, integer(2))}),
                    integer(2));
            }
        }
    }

    // DLMF 5.9.1: int_0^inf exp(-z x^mu) x^(nu-1) dx
    // = Gamma(nu/mu)/(mu z^(nu/mu)), for mu>0, Re(nu)>0, Re(z)>0.
    if (zeroToPositiveInfinity || wholeReal) {
        if (auto gammaKernel = matchPositiveHalfLineGammaKernel(
                expression, variable, builtins, mathematics, angles)) {
            if (proveGreaterThan(gammaKernel->nu, integer(0), builtins, mathematics, assumptions)
                && proveGreaterThan(gammaKernel->scale, integer(0), builtins, mathematics, assumptions)) {
                Expr ratio = simplification::Simplifier{}.simplify(
                    divide(builtins, mathematics, angles,
                        gammaKernel->nu, rational(gammaKernel->mu)), context);
                Expr scalePower = isExactOneExpr(gammaKernel->scale)
                    ? integer(1)
                    : power(builtins, mathematics, angles, gammaKernel->scale, ratio);
                Expr denominator = multiply(builtins, mathematics, angles, {
                    rational(gammaKernel->mu), std::move(scalePower)});
                Expr result = multiply(builtins, mathematics, angles, {
                    gammaKernel->outside,
                    divide(builtins, mathematics, angles,
                        gammaExactKernelValue(ratio, builtins, mathematics), denominator)});

                if (wholeReal) {
                    // exp(-a x^m) is even for positive even integer m.  General
                    // x^(nu-1) still needs branch/parity analysis on the negative half-line,
                    // so the whole-line shortcut remains restricted to nu=1.
                    if (!(gammaKernel->nu == integer(1))
                        || !gammaKernel->mu.isInteger()
                        || !(gammaKernel->mu > Rational{BigInt{0}})
                        || !(gammaKernel->mu.numerator() % BigInt{2}).isZero())
                        return std::nullopt;
                    result = multiply(builtins, mathematics, angles, {integer(2), result});
                }
                return simplification::Simplifier{}.simplify(std::move(result), context);
            }
        }
    }

    // Log moments on [0,1]: d^m/da^m int_0^1 x^a dx
    // = (-1)^m m!/(a+1)^(m+1), for a>-1.
    if (isExactZeroExpr(lower) && isExactOneExpr(upper)) {
        if (auto moment = logMomentKernel(expression, variable, builtins)) {
            if (moment->logPower <= 64
                && proveGreaterThan(moment->exponent, integer(-1),
                    builtins, mathematics, assumptions)) {
                BigInt factorial{1};
                for (std::uint64_t k = 2; k <= moment->logPower; ++k)
                    factorial *= BigInt::fromUnsigned(k);
                if (moment->logPower % 2 != 0)
                    factorial = -factorial;
                Expr next = simplify(add(builtins, mathematics, angles, {
                    moment->exponent, integer(1)}), builtins, mathematics, angles);
                return simplification::Simplifier{}.simplify(
                    divide(builtins, mathematics, angles,
                        Expr{Number{factorial}},
                        power(builtins, mathematics, angles, next,
                            integer(static_cast<std::int64_t>(moment->logPower + 1)))),
                    context);
            }
        }

        // t=x^q maps the endpoint algebraic singularity to Euler's beta integral.
        if (auto q = reciprocalSqrtOneMinusPower(expression, variable, builtins)) {
            Expr qExpr = rational(*q);
            Expr inverseQ = rational(Rational{BigInt{1}} / *q);
            return simplification::Simplifier{}.simplify(
                divide(builtins, mathematics, angles,
                    call(builtins, BuiltinId::Beta, {
                        inverseQ, rational(Rational{BigInt{1}, BigInt{2}})}),
                    qExpr), context);
        }
    }

    // Euler beta integral on [0,1].
    if (isExactZeroExpr(lower) && isExactOneExpr(upper)) {
        if (auto betaKernel = matchEulerBetaKernel(
                expression, variable, builtins, mathematics, angles)) {
            if (proveGreaterThan(betaKernel->first, integer(0), builtins, mathematics, assumptions)
                && proveGreaterThan(betaKernel->second, integer(0), builtins, mathematics, assumptions))
                return call(builtins, BuiltinId::Beta,
                    {std::move(betaKernel->first), std::move(betaKernel->second)});
        }
    }

    // int_0^a x^n dx = a^(n+1)/(n+1), for a>0 and n>-1.
    if (isExactZeroExpr(lower) && !isInfinityExpr(upper, infinitySymbol)) {
        if (auto exponent = variablePowerExponent(expression, variable, builtins)) {
            Expr next = simplification::Simplifier{}.simplify(
                add(builtins, mathematics, angles, {*exponent, integer(1)}), context);
            if (proveGreaterThan(upper, integer(0), builtins, mathematics, assumptions)
                && proveGreaterThan(*exponent, integer(-1), builtins, mathematics, assumptions)) {
                Expr numerator = isExactOneExpr(upper)
                    ? integer(1)
                    : power(builtins, mathematics, angles, upper, next);
                return simplification::Simplifier{}.simplify(
                    divide(builtins, mathematics, angles, std::move(numerator), next), context);
            }
        }
    }

    return std::nullopt;
}

[[nodiscard]] bool unresolvedNeedsDomainInformation(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    if (!containsHead(expression, builtins, BuiltinId::Abs)
        && !containsHead(expression, builtins, BuiltinId::Sign))
        return false;
    const mathematics::ValueFacts facts = mathematics::inferValueFacts(
        Expr{variable}, builtins, mathematics, assumptions);
    return !facts.isProvablyReal();
}

} // namespace

std::optional<Expr> normalizeRationalExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return normalizeRationalFunction(expression, variable, builtins, mathematics, angles);
}

IntegrationResult integrateExpressionDetailed(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const simplification::SimplificationContext context{
        builtins, mathematics, angles, assumptions};
    const Expr prepared = simplification::Simplifier{}.simplify(expression, context);
    if (const auto series = parseSeriesData(prepared, builtins)) {
        if (auto primitive = integrateSeriesExpression(
                *series, variable, builtins, mathematics, angles))
            return IntegrationResult{std::move(*primitive), IntegrationDisposition::Solved};
    }
    Expr primitive = integrateCore(prepared, variable, builtins, mathematics, angles, 0);

    // rule採用後のprimitiveだけをcanonical化する。derivative-backをruntime gateへ
    // 昇格させず、積分能力を維持したままformat/reparse時のAST順序を安定させる。
    primitive = simplification::Simplifier{}.simplify(primitive, context);
    if (!containsHead(primitive, builtins, BuiltinId::SymbolicIntegral))
        return IntegrationResult{std::move(primitive), IntegrationDisposition::Solved};

    const Expr directUnresolved = unresolved(prepared, variable, builtins);
    if (primitive != directUnresolved)
        return IntegrationResult{std::move(primitive), IntegrationDisposition::Partial};
    if (unresolvedNeedsDomainInformation(
            prepared, variable, builtins, mathematics, assumptions))
        return IntegrationResult{
            std::move(primitive), IntegrationDisposition::ConditionsRequired};
    if (isRecognizedNoFiniteClosedFormFamily(prepared, variable, builtins))
        return IntegrationResult{
            std::move(primitive), IntegrationDisposition::KnownNoFiniteClosedForm};
    return IntegrationResult{
        std::move(primitive), IntegrationDisposition::UnsupportedByEngine};
}

Expr integrateExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return integrateExpressionDetailed(
        expression, variable, builtins, mathematics, angles, assumptions).expression;
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
    const simplification::SimplificationContext context{
        builtins, mathematics, angles, assumptions};
    const Expr prepared = simplification::Simplifier{}.simplify(expression, context);
    if (auto known = tryAssumptionAwareDefiniteIntegral(
            prepared, variable, lower, upper, builtins, mathematics, angles,
            infinitySymbol, assumptions))
        return *known;

    Expr primitive = integrateCore(prepared, variable, builtins, mathematics, angles, 0);
    primitive = simplification::Simplifier{}.simplify(primitive, context);
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
