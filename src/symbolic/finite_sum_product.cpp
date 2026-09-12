// 記号有限和・積の一般kernel
#include "finite_sum_product.hpp"

#include "builtins/exact_operations.hpp"
#include "expression/exact_value.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "symbolic/algebra_transforms.hpp"
#include "symbolic/cases.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
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

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr unsignedInteger(std::size_t value) {
    return Expr{Number{BigInt::fromUnsigned(value)}};
}

[[nodiscard]] bool isExactIntegerValue(const Expr& expression, std::int64_t value) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return false;
    return expression.asNumber().asReal().asInteger() == BigInt{value};
}

[[nodiscard]] bool unitStepIterator(const evaluation::TableIteratorSpec& iterator) {
    if (iterator.rangeArguments.size() == 1 || iterator.rangeArguments.size() == 2)
        return true;
    return iterator.rangeArguments.size() == 3
        && isExactIntegerValue(iterator.rangeArguments[2], 1);
}

struct Bounds final {
    Expr lower;
    Expr upper;
};

[[nodiscard]] std::optional<Bounds> symbolicBounds(
    const evaluation::TableIteratorSpec& iterator) {
    if (!unitStepIterator(iterator))
        return std::nullopt;
    if (iterator.rangeArguments.size() == 1)
        return Bounds{integer(1), iterator.rangeArguments[0]};
    if (iterator.rangeArguments.size() == 2 || iterator.rangeArguments.size() == 3)
        return Bounds{iterator.rangeArguments[0], iterator.rangeArguments[1]};
    return std::nullopt;
}

[[nodiscard]] Expr countExpression(
    const Bounds& bounds,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return builtins::exact::add({
        builtins::exact::subtract(bounds.upper, bounds.lower, builtins, mathematics, angles),
        integer(1)}, builtins, mathematics, angles);
}

[[nodiscard]] Expr power(
    Expr base,
    Expr exponent,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return builtins::exact::call(
        BuiltinId::Power, {std::move(base), std::move(exponent)},
        builtins, mathematics, angles);
}

[[nodiscard]] Expr comparison(
    BuiltinId id,
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins) {
    return Expr::call(builtins.symbol(id), {std::move(lhs), std::move(rhs)});
}

[[nodiscard]] Expr guardedForNonemptyRange(
    Expr value,
    const Bounds& bounds,
    Expr emptyIdentity,
    const evaluation::BuiltinRegistry& builtins) {
    // finite iteratorのempty semanticsをsymbolic boundsでも失わない。
    // condition未確定ならCasesがそのまま保持されるため，proof-only方針を維持する。
    Expr condition = comparison(
        BuiltinId::LessEqual, bounds.lower, bounds.upper, builtins);
    return detail::makeCases(builtins, {
        detail::makeCaseBranch(builtins, std::move(value), std::move(condition)),
        detail::makeCaseBranch(builtins, std::move(emptyIdentity))});
}

[[nodiscard]] BigInt binomialBig(std::size_t n, std::size_t k) {
    if (k > n)
        return BigInt{0};
    k = std::min(k, n - k);
    BigInt result{1};
    for (std::size_t i = 1; i <= k; ++i) {
        result *= BigInt::fromUnsigned(n - k + i);
        result /= BigInt::fromUnsigned(i);
    }
    return result;
}

[[nodiscard]] std::optional<Expr> polynomialAntidifference(
    const Expr& body,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    constexpr std::size_t maximumDegree = 128;
    const auto polynomial = toExpressionPolynomial(
        body, variable, builtins, mathematics, angles,
        PolynomialConversionOptions{maximumDegree, 2048});
    if (!polynomial)
        return std::nullopt;

    const std::size_t degree = polynomial->degree();
    std::vector<Expr> q(degree + 2, integer(0));

    // Q(k+1)-Q(k)=P(k) を高次係数から三角的に解く。
    // Bernoulli/Faulhaberを個別表にせず，同じexact差分kernelで任意次数を扱う。
    for (std::size_t reverse = 0; reverse <= degree; ++reverse) {
        const std::size_t j = degree - reverse;
        std::vector<Expr> contributionTerms;
        for (std::size_t m = j + 2; m <= degree + 1; ++m) {
            const BigInt coefficient = binomialBig(m, j);
            if (coefficient.isZero())
                continue;
            contributionTerms.push_back(builtins::exact::multiply(
                {Expr{Number{coefficient}}, q[m]}, builtins, mathematics, angles));
        }
        Expr contribution = contributionTerms.empty()
            ? integer(0)
            : builtins::exact::add(
                std::move(contributionTerms), builtins, mathematics, angles);
        Expr numerator = builtins::exact::subtract(
            polynomial->coefficient(j), std::move(contribution),
            builtins, mathematics, angles);
        q[j + 1] = builtins::exact::divide(
            std::move(numerator), unsignedInteger(j + 1), builtins, mathematics, angles);
    }

    return expressionPolynomialToCollectedExpr(
        ExpressionPolynomial{std::move(q)}, variable,
        builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> polynomialSum(
    const Expr& body,
    const expression::Symbol& variable,
    const Bounds& bounds,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto antidifference = polynomialAntidifference(
        body, variable, builtins, mathematics, angles);
    if (!antidifference)
        return std::nullopt;

    Expr upperPlusOne = builtins::exact::add(
        {bounds.upper, integer(1)}, builtins, mathematics, angles);
    Expr atUpper = substituteSymbol(*antidifference, variable, upperPlusOne);
    Expr atLower = substituteSymbol(*antidifference, variable, bounds.lower);
    return builtins::exact::subtract(
        std::move(atUpper), std::move(atLower), builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> normalizeRationalTelescopingBody(
    const Expr& body,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!builtins.isCallTo(body, BuiltinId::Divide)
        || body.asCall().arguments.size() != 2
        || containsSymbol(body.asCall().arguments[0], variable))
        return std::nullopt;

    const Expr& numerator = body.asCall().arguments[0];
    Expr denominator = body.asCall().arguments[1];
    // 入力が既に積なら記号係数を含むfactor境界を壊さない。展開済みの
    // i^2+i等だけ一般Q[x] factorへ送り，線形因子を回収する。
    if (!builtins.isCallTo(denominator, BuiltinId::Multiply))
        denominator = factorExpression(denominator, builtins, mathematics, angles);
    std::vector<Expr> factors;
    if (builtins.isCallTo(denominator, BuiltinId::Multiply))
        factors.assign(
            denominator.asCall().arguments.begin(), denominator.asCall().arguments.end());
    else
        factors.push_back(std::move(denominator));

    std::vector<Expr> independent;
    std::vector<Expr> varying;
    for (const Expr& factor : factors) {
        if (containsSymbol(factor, variable))
            varying.push_back(factor);
        else
            independent.push_back(factor);
    }
    if (varying.size() != 2)
        return std::nullopt;

    const auto left = toExpressionPolynomial(
        varying[0], variable, builtins, mathematics, angles,
        PolynomialConversionOptions{1, 16});
    const auto right = toExpressionPolynomial(
        varying[1], variable, builtins, mathematics, angles,
        PolynomialConversionOptions{1, 16});
    if (!left || !right || left->degree() != 1 || right->degree() != 1)
        return std::nullopt;
    const auto leftSlope = expression::exact::realRational(left->coefficient(1));
    const auto rightSlope = expression::exact::realRational(right->coefficient(1));
    if (!leftSlope || !rightSlope || leftSlope->isZero()
        || !(*leftSlope == *rightSlope))
        return std::nullopt;

    // B-Aがexact定数なら 1/(A B)=(1/A-1/B)/(B-A) である。
    // shift幅の整数性は後段の構造照合で確認し，ここでは候補変形だけを作る。
    Expr delta = builtins::exact::subtract(
        right->coefficient(0), left->coefficient(0),
        builtins, mathematics, angles);
    const auto deltaValue = expression::exact::realRational(delta);
    if (!deltaValue || deltaValue->isZero())
        return std::nullopt;

    Expr constantDenominator = independent.empty()
        ? integer(1)
        : builtins::exact::multiply(
            std::move(independent), builtins, mathematics, angles);
    constantDenominator = builtins::exact::multiply(
        {std::move(constantDenominator), std::move(delta)},
        builtins, mathematics, angles);
    Expr coefficient = builtins::exact::divide(
        numerator, std::move(constantDenominator),
        builtins, mathematics, angles);
    Expr lhs = builtins::exact::divide(
        coefficient, varying[0], builtins, mathematics, angles);
    Expr rhs = builtins::exact::divide(
        coefficient, varying[1], builtins, mathematics, angles);
    return builtins::exact::subtract(
        std::move(lhs), std::move(rhs), builtins, mathematics, angles);
}

[[nodiscard]] Expr shiftedIteratorValue(
    const Expr& endpoint,
    std::size_t offset,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (offset == 0)
        return endpoint;
    return builtins::exact::add(
        {endpoint, unsignedInteger(offset)}, builtins, mathematics, angles);
}

[[nodiscard]] Expr telescopingBoundaryBlock(
    const Expr& function,
    const expression::Symbol& variable,
    const Expr& endpoint,
    std::size_t firstOffset,
    std::size_t count,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> terms;
    terms.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        Expr point = shiftedIteratorValue(
            endpoint, firstOffset + i, builtins, mathematics, angles);
        terms.push_back(substituteSymbol(function, variable, point));
    }
    return terms.size() == 1
        ? terms.front()
        : builtins::exact::add(std::move(terms), builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> telescopingSum(
    const Expr& body,
    const expression::Symbol& variable,
    const Bounds& bounds,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!builtins.isCallTo(body, BuiltinId::Subtract)
        || body.asCall().arguments.size() != 2)
        return std::nullopt;

    const Expr& lhs = body.asCall().arguments[0];
    const Expr& rhs = body.asCall().arguments[1];
    const Expr normalizedLhs = builtins::exact::simplify(lhs, builtins, mathematics, angles);
    const Expr normalizedRhs = builtins::exact::simplify(rhs, builtins, mathematics, angles);

    constexpr std::size_t maximumShift = 64;
    for (std::size_t shift = 1; shift <= maximumShift; ++shift) {
        Expr variablePlusShift = shiftedIteratorValue(
            Expr{variable}, shift, builtins, mathematics, angles);
        const Expr shiftedRhs = builtins::exact::simplify(
            substituteSymbol(rhs, variable, variablePlusShift),
            builtins, mathematics, angles);
        if (shiftedRhs == normalizedLhs) {
            Expr upperBlock = telescopingBoundaryBlock(
                rhs, variable, bounds.upper, 1, shift,
                builtins, mathematics, angles);
            Expr lowerBlock = telescopingBoundaryBlock(
                rhs, variable, bounds.lower, 0, shift,
                builtins, mathematics, angles);
            return builtins::exact::subtract(
                std::move(upperBlock), std::move(lowerBlock),
                builtins, mathematics, angles);
        }

        const Expr shiftedLhs = builtins::exact::simplify(
            substituteSymbol(lhs, variable, variablePlusShift),
            builtins, mathematics, angles);
        if (shiftedLhs == normalizedRhs) {
            Expr lowerBlock = telescopingBoundaryBlock(
                lhs, variable, bounds.lower, 0, shift,
                builtins, mathematics, angles);
            Expr upperBlock = telescopingBoundaryBlock(
                lhs, variable, bounds.upper, 1, shift,
                builtins, mathematics, angles);
            return builtins::exact::subtract(
                std::move(lowerBlock), std::move(upperBlock),
                builtins, mathematics, angles);
        }
    }
    return std::nullopt;
}

struct GeometricTerm final {
    Expr coefficient;
    Expr base;
    Expr exponentOffset;
    Expr exponentSlope;
};

[[nodiscard]] std::optional<GeometricTerm> geometricTerm(
    const Expr& body,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> factors;
    if (builtins.isCallTo(body, BuiltinId::Multiply))
        factors.assign(body.asCall().arguments.begin(), body.asCall().arguments.end());
    else
        factors.push_back(body);

    std::optional<Expr> varyingPower;
    std::vector<Expr> constants;
    for (const Expr& factor : factors) {
        if (!containsSymbol(factor, variable)) {
            constants.push_back(factor);
            continue;
        }
        if (varyingPower || !builtins.isCallTo(factor, BuiltinId::Power)
            || factor.asCall().arguments.size() != 2)
            return std::nullopt;
        varyingPower = factor;
    }
    if (!varyingPower)
        return std::nullopt;

    const auto& powerArgs = varyingPower->asCall().arguments;
    if (containsSymbol(powerArgs[0], variable))
        return std::nullopt;
    const auto exponent = toExpressionPolynomial(
        powerArgs[1], variable, builtins, mathematics, angles,
        PolynomialConversionOptions{1, 16});
    if (!exponent || exponent->degree() != 1)
        return std::nullopt;
    if (containsSymbol(exponent->coefficient(0), variable)
        || containsSymbol(exponent->coefficient(1), variable))
        return std::nullopt;

    Expr coefficient = constants.empty()
        ? integer(1)
        : builtins::exact::multiply(
            std::move(constants), builtins, mathematics, angles);
    return GeometricTerm{
        std::move(coefficient), powerArgs[0],
        exponent->coefficient(0), exponent->coefficient(1)};
}

[[nodiscard]] std::optional<Expr> geometricSum(
    const Expr& body,
    const expression::Symbol& variable,
    const Bounds& bounds,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto term = geometricTerm(body, variable, builtins, mathematics, angles);
    if (!term)
        return std::nullopt;

    const Expr count = countExpression(bounds, builtins, mathematics, angles);
    const Expr exponentAtLower = builtins::exact::add({
        builtins::exact::multiply(
            {term->exponentSlope, bounds.lower}, builtins, mathematics, angles),
        term->exponentOffset}, builtins, mathematics, angles);
    Expr first = builtins::exact::multiply({
        term->coefficient,
        power(term->base, exponentAtLower, builtins, mathematics, angles)},
        builtins, mathematics, angles);
    Expr ratio = power(
        term->base, term->exponentSlope, builtins, mathematics, angles);

    Expr equalValue = builtins::exact::multiply(
        {first, count}, builtins, mathematics, angles);
    Expr numerator = builtins::exact::subtract(
        integer(1), power(ratio, count, builtins, mathematics, angles),
        builtins, mathematics, angles);
    Expr denominator = builtins::exact::subtract(
        integer(1), ratio, builtins, mathematics, angles);
    Expr genericValue = builtins::exact::multiply({
        std::move(first),
        builtins::exact::divide(
            std::move(numerator), std::move(denominator),
            builtins, mathematics, angles)},
        builtins, mathematics, angles);

    if (ratio == integer(1))
        return equalValue;
    if (const auto exactRatio = expression::exact::realRational(ratio);
        exactRatio && !(*exactRatio == Rational{BigInt{1}}))
        return genericValue;

    Expr ratioEqualOne = comparison(BuiltinId::Equal, ratio, integer(1), builtins);
    return detail::makeCases(builtins, {
        detail::makeCaseBranch(builtins, std::move(equalValue), std::move(ratioEqualOne)),
        detail::makeCaseBranch(builtins, std::move(genericValue))});
}

[[nodiscard]] std::optional<Expr> binomialPowerSum(
    const Expr& body,
    const expression::Symbol& variable,
    const Bounds& bounds,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    // comb[n,k] z^k -> (1+z)^n。hypergeometric recognizerを実際の有限和へ
    // 接続する最小の基礎familyとして扱う。
    std::optional<Expr> combination;
    Expr weight = integer(1);
    if (builtins.isCallTo(body, BuiltinId::Combination)) {
        combination = body;
    }
    else if (builtins.isCallTo(body, BuiltinId::Multiply)) {
        std::vector<Expr> other;
        for (const Expr& factor : body.asCall().arguments) {
            if (!combination && builtins.isCallTo(factor, BuiltinId::Combination)) {
                combination = factor;
                continue;
            }
            other.push_back(factor);
        }
        if (!combination)
            return std::nullopt;
        weight = other.empty() ? integer(1)
            : builtins::exact::multiply(std::move(other), builtins, mathematics, angles);
    }
    else {
        return std::nullopt;
    }

    if (!combination || combination->asCall().arguments.size() != 2)
        return std::nullopt;
    const Expr& n = combination->asCall().arguments[0];
    if (!(combination->asCall().arguments[1].isSymbol()
        && combination->asCall().arguments[1].asSymbol().sameIdentity(variable)))
        return std::nullopt;
    if (!isExactIntegerValue(bounds.lower, 0) || !(bounds.upper == n))
        return std::nullopt;

    if (isExactIntegerValue(weight, 1))
        return power(integer(2), n, builtins, mathematics, angles);

    // weightはz^kのみを受理する。
    if (!builtins.isCallTo(weight, BuiltinId::Power)
        || weight.asCall().arguments.size() != 2
        || !(weight.asCall().arguments[1].isSymbol()
            && weight.asCall().arguments[1].asSymbol().sameIdentity(variable))
        || containsSymbol(weight.asCall().arguments[0], variable))
        return std::nullopt;
    Expr onePlusZ = builtins::exact::add(
        {integer(1), weight.asCall().arguments[0]}, builtins, mathematics, angles);
    return power(std::move(onePlusZ), n, builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> telescopingProduct(
    const Expr& body,
    const expression::Symbol& variable,
    const Bounds& bounds,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!builtins.isCallTo(body, BuiltinId::Divide)
        || body.asCall().arguments.size() != 2)
        return std::nullopt;
    const Expr& numerator = body.asCall().arguments[0];
    const Expr& denominator = body.asCall().arguments[1];
    Expr variablePlusOne = builtins::exact::add(
        {Expr{variable}, integer(1)}, builtins, mathematics, angles);
    const Expr shiftedDenominator = builtins::exact::simplify(
        substituteSymbol(denominator, variable, variablePlusOne),
        builtins, mathematics, angles);
    if (shiftedDenominator == builtins::exact::simplify(
            numerator, builtins, mathematics, angles)) {
        Expr upperPlusOne = builtins::exact::add(
            {bounds.upper, integer(1)}, builtins, mathematics, angles);
        return builtins::exact::divide(
            substituteSymbol(denominator, variable, upperPlusOne),
            substituteSymbol(denominator, variable, bounds.lower),
            builtins, mathematics, angles);
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> affineProduct(
    const Expr& body,
    const expression::Symbol& variable,
    const Bounds& bounds,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto polynomial = toExpressionPolynomial(
        body, variable, builtins, mathematics, angles,
        PolynomialConversionOptions{1, 16});
    if (!polynomial || polynomial->degree() > 1)
        return std::nullopt;

    const Expr count = countExpression(bounds, builtins, mathematics, angles);
    if (polynomial->degree() == 0)
        return power(polynomial->coefficient(0), count, builtins, mathematics, angles);

    const Expr& intercept = polynomial->coefficient(0);
    const Expr& slope = polynomial->coefficient(1);
    // 記号slopeの除算でdefinednessを捏造しない。exact非零実数係数だけ一般化する。
    const auto slopeValue = expression::exact::realRational(slope);
    if (!slopeValue || slopeValue->isZero())
        return std::nullopt;

    if (*slopeValue == Rational{BigInt{1}}
        && isExactIntegerValue(intercept, 0)
        && isExactIntegerValue(bounds.lower, 1))
        return Expr::call(
            builtins.symbol(BuiltinId::Factorial), {bounds.upper});

    Expr start = builtins::exact::add({
        bounds.lower,
        builtins::exact::divide(intercept, slope, builtins, mathematics, angles)},
        builtins, mathematics, angles);
    Expr rising = Expr::call(
        builtins.symbol(BuiltinId::RisingFactorial), {std::move(start), count});
    if (*slopeValue == Rational{BigInt{1}})
        return rising;
    return builtins::exact::multiply({
        power(slope, count, builtins, mathematics, angles),
        std::move(rising)}, builtins, mathematics, angles);
}

[[nodiscard]] RationalPolynomial rpMonic(const RationalPolynomial& polynomial) {
    if (polynomial.isZero())
        return polynomial;
    const Rational leading = polynomial.coefficient(polynomial.degree());
    std::vector<Rational> coefficients(
        polynomial.coefficients().begin(), polynomial.coefficients().end());
    for (Rational& coefficient : coefficients)
        coefficient /= leading;
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] RationalPolynomial rpAdd(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs,
    bool subtract = false) {
    std::vector<Rational> coefficients(
        std::max(lhs.degree(), rhs.degree()) + 1, Rational{BigInt{0}});
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        coefficients[i] += lhs.coefficient(i);
    for (std::size_t i = 0; i <= rhs.degree(); ++i) {
        if (subtract)
            coefficients[i] -= rhs.coefficient(i);
        else
            coefficients[i] += rhs.coefficient(i);
    }
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] std::optional<RationalPolynomial> rpMultiply(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs,
    std::size_t maximumDegree = 96) {
    if (lhs.isZero() || rhs.isZero())
        return RationalPolynomial{};
    if (lhs.degree() > maximumDegree - std::min(rhs.degree(), maximumDegree)
        || lhs.degree() + rhs.degree() > maximumDegree)
        return std::nullopt;
    std::vector<Rational> coefficients(
        lhs.degree() + rhs.degree() + 1, Rational{BigInt{0}});
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        for (std::size_t j = 0; j <= rhs.degree(); ++j)
            coefficients[i + j] += lhs.coefficient(i) * rhs.coefficient(j);
    return RationalPolynomial{std::move(coefficients)};
}

struct RpDivision final {
    RationalPolynomial quotient;
    RationalPolynomial remainder;
};

[[nodiscard]] std::optional<RpDivision> rpDivide(
    const RationalPolynomial& dividend,
    const RationalPolynomial& divisor) {
    if (divisor.isZero())
        return std::nullopt;
    if (dividend.degree() < divisor.degree())
        return RpDivision{RationalPolynomial{}, dividend};
    std::vector<Rational> remainder(
        dividend.coefficients().begin(), dividend.coefficients().end());
    std::vector<Rational> quotient(
        dividend.degree() - divisor.degree() + 1, Rational{BigInt{0}});
    const Rational leading = divisor.coefficient(divisor.degree());
    for (std::size_t degree = dividend.degree() + 1; degree-- > divisor.degree();) {
        const Rational amount = remainder[degree] / leading;
        const std::size_t shift = degree - divisor.degree();
        quotient[shift] += amount;
        for (std::size_t i = 0; i <= divisor.degree(); ++i)
            remainder[i + shift] -= amount * divisor.coefficient(i);
    }
    remainder.resize(divisor.degree());
    return RpDivision{
        RationalPolynomial{std::move(quotient)},
        RationalPolynomial{std::move(remainder)}};
}

[[nodiscard]] std::optional<RationalPolynomial> rpExactQuotient(
    const RationalPolynomial& dividend,
    const RationalPolynomial& divisor) {
    auto division = rpDivide(dividend, divisor);
    if (!division || !division->remainder.isZero())
        return std::nullopt;
    return std::move(division->quotient);
}

[[nodiscard]] RationalPolynomial rpGcd(
    RationalPolynomial lhs,
    RationalPolynomial rhs) {
    while (!rhs.isZero()) {
        auto division = rpDivide(lhs, rhs);
        if (!division)
            return RationalPolynomial{{Rational{BigInt{1}}}};
        lhs = std::move(rhs);
        rhs = std::move(division->remainder);
    }
    return rpMonic(lhs);
}

[[nodiscard]] RationalPolynomial rpShift(
    const RationalPolynomial& polynomial,
    std::int64_t shift) {
    std::vector<Rational> coefficients(
        polynomial.degree() + 1, Rational{BigInt{0}});
    const Rational offset{BigInt{shift}};
    for (std::size_t exponent = 0; exponent <= polynomial.degree(); ++exponent) {
        Rational offsetPower{BigInt{1}};
        for (std::size_t reverse = 0; reverse <= exponent; ++reverse) {
            const std::size_t outputExponent = exponent - reverse;
            coefficients[outputExponent] += polynomial.coefficient(exponent)
                * Rational{binomialBig(exponent, outputExponent)} * offsetPower;
            offsetPower *= offset;
        }
    }
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] std::optional<RationalPolynomial> rpLcm(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs,
    std::size_t maximumDegree) {
    const RationalPolynomial common = rpGcd(lhs, rhs);
    auto quotient = rpExactQuotient(lhs, common);
    if (!quotient)
        return std::nullopt;
    auto product = rpMultiply(*quotient, rhs, maximumDegree);
    if (!product)
        return std::nullopt;
    return rpMonic(*product);
}

[[nodiscard]] std::optional<RationalPolynomial> rpPower(
    RationalPolynomial base,
    std::uint64_t exponent,
    std::size_t maximumDegree) {
    RationalPolynomial result{{Rational{BigInt{1}}}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0) {
            auto product = rpMultiply(result, base, maximumDegree);
            if (!product)
                return std::nullopt;
            result = std::move(*product);
        }
        exponent >>= 1U;
        if (exponent != 0) {
            auto square = rpMultiply(base, base, maximumDegree);
            if (!square)
                return std::nullopt;
            base = std::move(*square);
        }
    }
    return result;
}

struct ExactRationalFunction final {
    RationalPolynomial numerator;
    RationalPolynomial denominator{{Rational{BigInt{1}}}};
};

[[nodiscard]] std::optional<ExactRationalFunction> normalizeRationalFunction(
    ExactRationalFunction function) {
    if (function.denominator.isZero())
        return std::nullopt;
    if (function.numerator.isZero())
        return ExactRationalFunction{};
    const RationalPolynomial common = rpGcd(
        function.numerator, function.denominator);
    auto numerator = rpExactQuotient(function.numerator, common);
    auto denominator = rpExactQuotient(function.denominator, common);
    if (!numerator || !denominator)
        return std::nullopt;
    const Rational leading = denominator->coefficient(denominator->degree());
    std::vector<Rational> numeratorCoefficients(
        numerator->coefficients().begin(), numerator->coefficients().end());
    std::vector<Rational> denominatorCoefficients(
        denominator->coefficients().begin(), denominator->coefficients().end());
    for (Rational& coefficient : numeratorCoefficients)
        coefficient /= leading;
    for (Rational& coefficient : denominatorCoefficients)
        coefficient /= leading;
    return ExactRationalFunction{
        RationalPolynomial{std::move(numeratorCoefficients)},
        RationalPolynomial{std::move(denominatorCoefficients)}};
}

[[nodiscard]] std::optional<ExactRationalFunction> combineRationalFunctions(
    const ExactRationalFunction& lhs,
    const ExactRationalFunction& rhs,
    BuiltinId operation,
    std::size_t maximumDegree) {
    ExactRationalFunction result;
    if (operation == BuiltinId::Add || operation == BuiltinId::Subtract) {
        auto left = rpMultiply(lhs.numerator, rhs.denominator, maximumDegree);
        auto right = rpMultiply(rhs.numerator, lhs.denominator, maximumDegree);
        auto denominator = rpMultiply(lhs.denominator, rhs.denominator, maximumDegree);
        if (!left || !right || !denominator)
            return std::nullopt;
        result.numerator = rpAdd(*left, *right, operation == BuiltinId::Subtract);
        result.denominator = std::move(*denominator);
    }
    else if (operation == BuiltinId::Multiply) {
        auto numerator = rpMultiply(lhs.numerator, rhs.numerator, maximumDegree);
        auto denominator = rpMultiply(lhs.denominator, rhs.denominator, maximumDegree);
        if (!numerator || !denominator)
            return std::nullopt;
        result = ExactRationalFunction{std::move(*numerator), std::move(*denominator)};
    }
    else if (operation == BuiltinId::Divide) {
        if (rhs.numerator.isZero())
            return std::nullopt;
        auto numerator = rpMultiply(lhs.numerator, rhs.denominator, maximumDegree);
        auto denominator = rpMultiply(lhs.denominator, rhs.numerator, maximumDegree);
        if (!numerator || !denominator)
            return std::nullopt;
        result = ExactRationalFunction{std::move(*numerator), std::move(*denominator)};
    }
    else {
        return std::nullopt;
    }
    return normalizeRationalFunction(std::move(result));
}

[[nodiscard]] std::optional<std::int64_t> boundedIntegerExponent(
    const Expr& expression,
    std::int64_t magnitudeLimit) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    const BigInt& integerValue = expression.asNumber().asReal().asInteger();
    const auto magnitude = numeric::tryToUint64(integerValue.abs());
    if (!magnitude || *magnitude > static_cast<std::uint64_t>(magnitudeLimit))
        return std::nullopt;
    const auto value = static_cast<std::int64_t>(*magnitude);
    return integerValue.isNegative() ? -value : value;
}

[[nodiscard]] std::optional<ExactRationalFunction> extractRationalFunction(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    std::size_t maximumDegree,
    std::size_t depth = 0) {
    if (depth > 64)
        return std::nullopt;
    if (auto polynomial = toRationalPolynomial(
            expression, variable, builtins,
            PolynomialConversionOptions{maximumDegree, 2048}))
        return ExactRationalFunction{std::move(*polynomial),
            RationalPolynomial{{Rational{BigInt{1}}}}};
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;

    if (definition->id == BuiltinId::Negate && arguments.size() == 1) {
        auto value = extractRationalFunction(
            arguments[0], variable, builtins, maximumDegree, depth + 1);
        if (!value)
            return std::nullopt;
        std::vector<Rational> coefficients(
            value->numerator.coefficients().begin(),
            value->numerator.coefficients().end());
        for (Rational& coefficient : coefficients)
            coefficient = -coefficient;
        value->numerator = RationalPolynomial{std::move(coefficients)};
        return value;
    }
    if ((definition->id == BuiltinId::Add
            || definition->id == BuiltinId::Multiply)
        && !arguments.empty()) {
        ExactRationalFunction result = definition->id == BuiltinId::Add
            ? ExactRationalFunction{}
            : ExactRationalFunction{
                RationalPolynomial{{Rational{BigInt{1}}}},
                RationalPolynomial{{Rational{BigInt{1}}}}};
        for (const Expr& argument : arguments) {
            auto value = extractRationalFunction(
                argument, variable, builtins, maximumDegree, depth + 1);
            if (!value)
                return std::nullopt;
            auto combined = combineRationalFunctions(
                result, *value, definition->id, maximumDegree);
            if (!combined)
                return std::nullopt;
            result = std::move(*combined);
        }
        return result;
    }
    if ((definition->id == BuiltinId::Subtract
            || definition->id == BuiltinId::Divide)
        && arguments.size() == 2) {
        auto lhs = extractRationalFunction(
            arguments[0], variable, builtins, maximumDegree, depth + 1);
        auto rhs = extractRationalFunction(
            arguments[1], variable, builtins, maximumDegree, depth + 1);
        if (!lhs || !rhs)
            return std::nullopt;
        return combineRationalFunctions(
            *lhs, *rhs, definition->id, maximumDegree);
    }
    if (definition->id == BuiltinId::Power && arguments.size() == 2) {
        const auto exponent = boundedIntegerExponent(arguments[1], 16);
        if (!exponent)
            return std::nullopt;
        auto base = extractRationalFunction(
            arguments[0], variable, builtins, maximumDegree, depth + 1);
        if (!base || (*exponent < 0 && base->numerator.isZero()))
            return std::nullopt;
        const std::uint64_t magnitude = static_cast<std::uint64_t>(
            *exponent < 0 ? -*exponent : *exponent);
        auto numerator = rpPower(
            *exponent < 0 ? base->denominator : base->numerator,
            magnitude, maximumDegree);
        auto denominator = rpPower(
            *exponent < 0 ? base->numerator : base->denominator,
            magnitude, maximumDegree);
        if (!numerator || !denominator)
            return std::nullopt;
        return normalizeRationalFunction(
            ExactRationalFunction{std::move(*numerator), std::move(*denominator)});
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<std::vector<Rational>> solveExactLinearSystem(
    std::vector<std::vector<Rational>> matrix,
    std::size_t unknownCount) {
    std::size_t pivotRow = 0;
    std::vector<std::size_t> pivotColumns;
    for (std::size_t column = 0;
         column < unknownCount && pivotRow < matrix.size(); ++column) {
        std::size_t selected = pivotRow;
        while (selected < matrix.size() && matrix[selected][column].isZero())
            ++selected;
        if (selected == matrix.size())
            continue;
        std::swap(matrix[pivotRow], matrix[selected]);
        const Rational pivot = matrix[pivotRow][column];
        for (std::size_t j = column; j <= unknownCount; ++j)
            matrix[pivotRow][j] /= pivot;
        for (std::size_t row = 0; row < matrix.size(); ++row) {
            if (row == pivotRow || matrix[row][column].isZero())
                continue;
            const Rational amount = matrix[row][column];
            for (std::size_t j = column; j <= unknownCount; ++j)
                matrix[row][j] -= amount * matrix[pivotRow][j];
        }
        pivotColumns.push_back(column);
        ++pivotRow;
    }
    for (const auto& row : matrix) {
        bool allZero = true;
        for (std::size_t column = 0; column < unknownCount; ++column)
            allZero = allZero && row[column].isZero();
        if (allZero && !row[unknownCount].isZero())
            return std::nullopt;
    }
    std::vector<Rational> solution(unknownCount, Rational{BigInt{0}});
    for (std::size_t row = 0; row < pivotColumns.size(); ++row)
        solution[pivotColumns[row]] = matrix[row][unknownCount];
    return solution;
}

[[nodiscard]] bool samePolynomial(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    return lhs.coefficients() == rhs.coefficients();
}

[[nodiscard]] std::optional<ExactRationalFunction> abramovProperAntidifference(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator) {
    constexpr std::size_t maximumDispersion = 32;
    constexpr std::size_t maximumDegree = 96;
    constexpr std::size_t maximumEquations = 256;
    if (numerator.isZero() || denominator.degree() == 0
        || denominator.degree() > maximumDegree)
        return std::nullopt;

    std::size_t dispersion = 0;
    for (std::size_t shift = 1; shift <= maximumDispersion; ++shift) {
        const RationalPolynomial shifted = rpShift(
            denominator, static_cast<std::int64_t>(shift));
        if (rpGcd(denominator, shifted).degree() > 0)
            dispersion = shift;
    }

    RationalPolynomial universal = rpMonic(denominator);
    for (std::int64_t shift = -static_cast<std::int64_t>(dispersion);
         shift <= static_cast<std::int64_t>(dispersion); ++shift) {
        auto next = rpLcm(
            universal, rpShift(denominator, shift), maximumDegree);
        if (!next)
            return std::nullopt;
        universal = std::move(*next);
    }
    if (universal.degree() == 0 || universal.degree() > maximumDegree)
        return std::nullopt;

    const RationalPolynomial shiftedUniversal = rpShift(universal, 1);
    auto universalTimesDenominator = rpMultiply(
        universal, denominator, maximumEquations);
    auto shiftedTimesDenominator = rpMultiply(
        shiftedUniversal, denominator, maximumEquations);
    auto rhsFirst = rpMultiply(numerator, universal, maximumEquations);
    if (!universalTimesDenominator || !shiftedTimesDenominator || !rhsFirst)
        return std::nullopt;
    auto rhs = rpMultiply(*rhsFirst, shiftedUniversal, maximumEquations);
    if (!rhs)
        return std::nullopt;

    const std::size_t unknownCount = universal.degree();
    std::vector<RationalPolynomial> columns;
    columns.reserve(unknownCount);
    std::size_t equationCount = rhs->degree() + 1;
    for (std::size_t exponent = 0; exponent < unknownCount; ++exponent) {
        std::vector<Rational> monomialCoefficients(
            exponent + 1, Rational{BigInt{0}});
        monomialCoefficients.back() = Rational{BigInt{1}};
        RationalPolynomial monomial{std::move(monomialCoefficients)};
        auto left = rpMultiply(
            rpShift(monomial, 1), *universalTimesDenominator,
            maximumEquations);
        auto right = rpMultiply(
            monomial, *shiftedTimesDenominator, maximumEquations);
        if (!left || !right)
            return std::nullopt;
        columns.push_back(rpAdd(*left, *right, true));
        equationCount = std::max(equationCount, columns.back().degree() + 1);
    }
    if (equationCount > maximumEquations)
        return std::nullopt;
    std::vector<std::vector<Rational>> matrix(
        equationCount,
        std::vector<Rational>(unknownCount + 1, Rational{BigInt{0}}));
    for (std::size_t row = 0; row < equationCount; ++row) {
        for (std::size_t column = 0; column < unknownCount; ++column)
            matrix[row][column] = columns[column].coefficient(row);
        matrix[row][unknownCount] = rhs->coefficient(row);
    }
    auto solution = solveExactLinearSystem(std::move(matrix), unknownCount);
    if (!solution)
        return std::nullopt;
    RationalPolynomial antidifferenceNumerator{std::move(*solution)};

    // N(x+1)D(x)Q(x)-N(x)D(x+1)Q(x)=R(x)D(x)D(x+1)
    // をもう一度exactに構成し、線形系実装から独立したcertificateにする。
    auto certifiedLeftA = rpMultiply(
        rpShift(antidifferenceNumerator, 1), *universalTimesDenominator,
        maximumEquations);
    auto certifiedLeftB = rpMultiply(
        antidifferenceNumerator, *shiftedTimesDenominator,
        maximumEquations);
    if (!certifiedLeftA || !certifiedLeftB
        || !samePolynomial(rpAdd(*certifiedLeftA, *certifiedLeftB, true), *rhs))
        return std::nullopt;

    return normalizeRationalFunction(ExactRationalFunction{
        std::move(antidifferenceNumerator), std::move(universal)});
}

[[nodiscard]] std::optional<Expr> abramovRationalAntidifference(
    const Expr& body,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    constexpr std::size_t maximumDegree = 96;
    auto function = extractRationalFunction(
        body, variable, builtins, maximumDegree);
    if (!function || function->denominator.degree() == 0)
        return std::nullopt;
    auto division = rpDivide(function->numerator, function->denominator);
    if (!division || division->remainder.isZero())
        return std::nullopt;
    auto proper = abramovProperAntidifference(
        division->remainder, function->denominator);
    if (!proper)
        return std::nullopt;

    std::vector<Expr> parts;
    if (!division->quotient.isZero()) {
        Expr polynomialBody = polynomialToExpandedExpr(
            division->quotient, variable, builtins);
        auto polynomialPart = polynomialAntidifference(
            polynomialBody, variable, builtins, mathematics, angles);
        if (!polynomialPart)
            return std::nullopt;
        parts.push_back(std::move(*polynomialPart));
    }
    Expr rationalNumerator = polynomialToExpandedExpr(
        proper->numerator, variable, builtins);
    Expr rationalDenominator = polynomialToExpandedExpr(
        proper->denominator, variable, builtins);
    parts.push_back(builtins::exact::divide(
        std::move(rationalNumerator), std::move(rationalDenominator),
        builtins, mathematics, angles));
    return parts.size() == 1
        ? std::optional<Expr>{std::move(parts.front())}
        : std::optional<Expr>{builtins::exact::add(
            std::move(parts), builtins, mathematics, angles)};
}

[[nodiscard]] std::optional<Expr> abramovRationalSum(
    const Expr& body,
    const expression::Symbol& variable,
    const Bounds& bounds,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto antidifference = abramovRationalAntidifference(
        body, variable, builtins, mathematics, angles);
    if (!antidifference)
        return std::nullopt;
    Expr upperPlusOne = builtins::exact::add(
        {bounds.upper, integer(1)}, builtins, mathematics, angles);
    Expr atUpper = substituteSymbol(*antidifference, variable, upperPlusOne);
    Expr atLower = substituteSymbol(*antidifference, variable, bounds.lower);
    return builtins::exact::subtract(
        std::move(atUpper), std::move(atLower),
        builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> sumIndependentCases(
    const Expr& body,
    const evaluation::TableIteratorSpec& iterator,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!builtins.isCallTo(body, BuiltinId::Cases))
        return std::nullopt;

    std::vector<Expr> branches;
    branches.reserve(body.asCall().arguments.size());
    for (const Expr& branchExpression : body.asCall().arguments) {
        if (!builtins.isCallTo(branchExpression, BuiltinId::CaseBranch)
            || branchExpression.asCall().arguments.empty()
            || branchExpression.asCall().arguments.size() > 2)
            return std::nullopt;
        const auto& branch = branchExpression.asCall().arguments;
        if (branch.size() == 2 && containsSymbol(branch[1], iterator.variable))
            return std::nullopt;
        auto value = finiteSymbolicSum(
            branch[0], iterator, builtins, mathematics, angles);
        if (!value)
            return std::nullopt;
        branches.push_back(detail::makeCaseBranch(
            builtins, std::move(*value),
            branch.size() == 2 ? std::optional<Expr>{branch[1]} : std::nullopt));
    }
    return detail::makeCases(builtins, std::move(branches));
}

[[nodiscard]] std::optional<Expr> productIndependentCases(
    const Expr& body,
    const evaluation::TableIteratorSpec& iterator,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!builtins.isCallTo(body, BuiltinId::Cases))
        return std::nullopt;

    std::vector<Expr> branches;
    branches.reserve(body.asCall().arguments.size());
    for (const Expr& branchExpression : body.asCall().arguments) {
        if (!builtins.isCallTo(branchExpression, BuiltinId::CaseBranch)
            || branchExpression.asCall().arguments.empty()
            || branchExpression.asCall().arguments.size() > 2)
            return std::nullopt;
        const auto& branch = branchExpression.asCall().arguments;
        if (branch.size() == 2 && containsSymbol(branch[1], iterator.variable))
            return std::nullopt;
        auto value = finiteSymbolicProduct(
            branch[0], iterator, builtins, mathematics, angles);
        if (!value)
            return std::nullopt;
        branches.push_back(detail::makeCaseBranch(
            builtins, std::move(*value),
            branch.size() == 2 ? std::optional<Expr>{branch[1]} : std::nullopt));
    }
    return detail::makeCases(builtins, std::move(branches));
}

[[nodiscard]] bool rationalFunctionInVariable(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (toRationalPolynomial(expression, variable, builtins,
            PolynomialConversionOptions{128, 1024}))
        return true;
    if (!builtins.isCallTo(expression, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return false;
    return toRationalPolynomial(expression.asCall().arguments[0], variable, builtins,
               PolynomialConversionOptions{128, 1024}).has_value()
        && toRationalPolynomial(expression.asCall().arguments[1], variable, builtins,
               PolynomialConversionOptions{128, 1024}).has_value();
}

} // namespace

std::optional<Expr> finiteSymbolicSum(
    const Expr& body,
    const evaluation::TableIteratorSpec& iterator,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto bounds = symbolicBounds(iterator);
    if (!bounds)
        return std::nullopt;

    if (auto distributed = sumIndependentCases(
            body, iterator, builtins, mathematics, angles))
        return distributed;

    if (!containsSymbol(body, iterator.variable)) {
        Expr value = builtins::exact::multiply({
            body, countExpression(*bounds, builtins, mathematics, angles)},
            builtins, mathematics, angles);
        return guardedForNonemptyRange(
            std::move(value), *bounds, integer(0), builtins);
    }

    if (auto result = binomialPowerSum(
            body, iterator.variable, *bounds, builtins, mathematics, angles))
        return guardedForNonemptyRange(
            std::move(*result), *bounds, integer(0), builtins);
    Expr telescopingBody = body;
    if (auto normalized = normalizeRationalTelescopingBody(
            body, iterator.variable, builtins, mathematics, angles))
        telescopingBody = std::move(*normalized);
    if (auto result = telescopingSum(
            telescopingBody, iterator.variable, *bounds, builtins, mathematics, angles))
        return guardedForNonemptyRange(
            std::move(*result), *bounds, integer(0), builtins);
    if (auto result = abramovRationalSum(
            body, iterator.variable, *bounds, builtins, mathematics, angles))
        return guardedForNonemptyRange(
            std::move(*result), *bounds, integer(0), builtins);
    if (auto result = polynomialSum(
            body, iterator.variable, *bounds, builtins, mathematics, angles))
        return guardedForNonemptyRange(
            std::move(*result), *bounds, integer(0), builtins);
    if (auto result = geometricSum(
            body, iterator.variable, *bounds, builtins, mathematics, angles))
        return guardedForNonemptyRange(
            std::move(*result), *bounds, integer(0), builtins);
    return std::nullopt;
}

std::optional<Expr> finiteSymbolicProduct(
    const Expr& body,
    const evaluation::TableIteratorSpec& iterator,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto bounds = symbolicBounds(iterator);
    if (!bounds)
        return std::nullopt;

    if (auto distributed = productIndependentCases(
            body, iterator, builtins, mathematics, angles))
        return distributed;

    if (!containsSymbol(body, iterator.variable)) {
        Expr value = power(
            body, countExpression(*bounds, builtins, mathematics, angles),
            builtins, mathematics, angles);
        return guardedForNonemptyRange(
            std::move(value), *bounds, integer(1), builtins);
    }
    if (auto result = telescopingProduct(
            body, iterator.variable, *bounds, builtins, mathematics, angles))
        return guardedForNonemptyRange(
            std::move(*result), *bounds, integer(1), builtins);
    if (auto result = affineProduct(
            body, iterator.variable, *bounds, builtins, mathematics, angles))
        return guardedForNonemptyRange(
            std::move(*result), *bounds, integer(1), builtins);
    return std::nullopt;
}

std::optional<HypergeometricTermRecognition> recognizeHypergeometricTerm(
    const Expr& body,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr variablePlusOne = builtins::exact::add(
        {Expr{variable}, integer(1)}, builtins, mathematics, angles);
    Expr shifted = substituteSymbol(body, variable, variablePlusOne);
    Expr ratio = builtins::exact::divide(
        std::move(shifted), body, builtins, mathematics, angles);
    ratio = builtins::exact::simplify(std::move(ratio), builtins, mathematics, angles);
    if (!rationalFunctionInVariable(ratio, variable, builtins))
        return std::nullopt;
    return HypergeometricTermRecognition{std::move(ratio)};
}

} // namespace mmcal::symbolic
