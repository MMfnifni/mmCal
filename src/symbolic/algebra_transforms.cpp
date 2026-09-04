// 展開・因数分解・collect
#include "algebra_transforms.hpp"

#include "mathematics/exact_algebra.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "symbolic/polynomial.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <optional>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace mmcal::symbolic {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

[[nodiscard]] Rational rational(std::int64_t value) { return Rational{BigInt{value}}; }

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return builtins.isCallTo(expression, id);
}

[[nodiscard]] std::optional<std::size_t> smallNonNegativeInteger(
    const Expr& expression,
    std::size_t maximum) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    const BigInt& value = expression.asNumber().asReal().asInteger();
    if (value.isNegative())
        return std::nullopt;
    const std::string text = value.toString();
    std::size_t result = 0;
    const auto conversion = std::from_chars(text.data(), text.data() + text.size(), result);
    if (conversion.ec != std::errc{} || conversion.ptr != text.data() + text.size()
        || result > maximum)
        return std::nullopt;
    return result;
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

[[nodiscard]] std::vector<Expr> addTerms(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (isHead(expression, builtins, BuiltinId::Add))
        return expression.asCall().arguments;
    return {expression};
}

[[nodiscard]] std::optional<Expr> distributeProduct(
    const std::vector<Expr>& factors,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t maximumTerms) {
    std::vector<std::vector<Expr>> choices;
    choices.reserve(factors.size());
    std::size_t termCount = 1;
    for (const Expr& factor : factors) {
        auto terms = addTerms(factor, builtins);
        if (!terms.empty() && termCount > maximumTerms / terms.size())
            return std::nullopt;
        termCount *= terms.size();
        if (termCount > maximumTerms)
            return std::nullopt;
        choices.push_back(std::move(terms));
    }

    std::vector<Expr> products;
    products.reserve(termCount);
    std::vector<std::size_t> indices(choices.size(), 0);
    for (std::size_t generated = 0; generated < termCount; ++generated) {
        std::vector<Expr> productFactors;
        productFactors.reserve(choices.size());
        for (std::size_t i = 0; i < choices.size(); ++i)
            productFactors.push_back(choices[i][indices[i]]);
        products.push_back(productFactors.size() == 1
            ? productFactors.front()
            : Expr::call(builtins.symbol(BuiltinId::Multiply), std::move(productFactors)));

        for (std::size_t position = choices.size(); position-- > 0;) {
            if (++indices[position] < choices[position].size())
                break;
            indices[position] = 0;
        }
    }

    Expr result = products.size() == 1
        ? products.front()
        : Expr::call(builtins.symbol(BuiltinId::Add), std::move(products));
    return simplify(std::move(result), builtins, mathematics, angles);
}

[[nodiscard]] Expr expandRecursive(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AlgebraTransformOptions& options) {
    if (!expression.isCall())
        return expression;

    const auto& call = expression.asCall();
    std::vector<Expr> arguments;
    arguments.reserve(call.arguments.size());
    for (const Expr& argument : call.arguments)
        arguments.push_back(expandRecursive(argument, builtins, mathematics, angles, options));

    Expr rebuilt = Expr::rebuildCall(call, arguments);
    if (isHead(rebuilt, builtins, BuiltinId::Subtract) && arguments.size() == 2) {
        rebuilt = Expr::call(
            builtins.symbol(BuiltinId::Add),
            {arguments[0], Expr::call(builtins.symbol(BuiltinId::Negate), {arguments[1]})});
    }
    else if (isHead(rebuilt, builtins, BuiltinId::Power) && arguments.size() == 2) {
        const auto exponent = smallNonNegativeInteger(arguments[1], options.maximumExpandedTerms);
        if (exponent && *exponent > 0 && isHead(arguments[0], builtins, BuiltinId::Add)) {
            std::vector<Expr> factors(*exponent, arguments[0]);
            if (const auto distributed = distributeProduct(
                factors, builtins, mathematics, angles, options.maximumExpandedTerms))
                rebuilt = *distributed;
        }
    }

    if (isHead(rebuilt, builtins, BuiltinId::Multiply)) {
        if (const auto distributed = distributeProduct(
            rebuilt.asCall().arguments,
            builtins,
            mathematics,
            angles,
            options.maximumExpandedTerms))
            rebuilt = *distributed;
    }

    return simplify(std::move(rebuilt), builtins, mathematics, angles);
}

[[nodiscard]] std::optional<Rational> perfectRationalSquareRoot(const Rational& value) {
    if (value.numerator().isNegative())
        return std::nullopt;
    const auto numerator = numeric::integerSqrt(value.numerator());
    const auto denominator = numeric::integerSqrt(value.denominator());
    if (!numerator.remainder.isZero() || !denominator.remainder.isZero())
        return std::nullopt;
    return Rational{numerator.root, denominator.root};
}

[[nodiscard]] std::optional<BigInt> exactIntegerCubeRoot(const BigInt& value) {
    if (value.isZero())
        return BigInt{0};
    const bool negative = value.isNegative();
    const BigInt magnitude = value.abs();
    const std::size_t highBit = (magnitude.bitLength() + 2) / 3 + 1;
    BigInt low{0};
    BigInt high = BigInt{1} << highBit;
    const BigInt one{1};
    while (high - low > one) {
        const BigInt middle = (low + high) / BigInt{2};
        const BigInt cube = middle * middle * middle;
        if (cube <= magnitude)
            low = middle;
        else
            high = middle;
    }
    if (low * low * low != magnitude)
        return std::nullopt;
    return negative ? -low : low;
}

[[nodiscard]] std::optional<Rational> perfectRationalCubeRoot(const Rational& value) {
    const auto numerator = exactIntegerCubeRoot(value.numerator());
    const auto denominator = exactIntegerCubeRoot(value.denominator());
    if (!numerator || !denominator)
        return std::nullopt;
    return Rational{*numerator, *denominator};
}

[[nodiscard]] std::optional<BigInt> exactIntegerNthRoot(
    const BigInt& value,
    std::uint64_t degree) {
    if (degree == 0)
        return std::nullopt;
    if (degree == 1 || value.isZero())
        return value;

    const bool negative = value.isNegative();
    if (negative && (degree % 2U) == 0U)
        return std::nullopt;
    const BigInt magnitude = value.abs();
    const std::size_t rootBits = (magnitude.bitLength() + degree - 1) / degree;
    BigInt low{0};
    BigInt high = BigInt{1} << (rootBits + 1);
    const BigInt one{1};
    while (high - low > one) {
        const BigInt middle = (low + high) / BigInt{2};
        if (numeric::pow(middle, degree) <= magnitude)
            low = middle;
        else
            high = middle;
    }
    if (numeric::pow(low, degree) != magnitude)
        return std::nullopt;
    return negative ? -low : low;
}

[[nodiscard]] std::optional<Rational> perfectRationalNthRoot(
    const Rational& value,
    std::uint64_t degree) {
    const auto numerator = exactIntegerNthRoot(value.numerator(), degree);
    const auto denominator = exactIntegerNthRoot(value.denominator(), degree);
    if (!numerator || !denominator)
        return std::nullopt;
    return Rational{*numerator, *denominator};
}

[[nodiscard]] RationalPolynomial multiplyRationalPolynomials(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    if (lhs.isZero() || rhs.isZero())
        return RationalPolynomial{};
    std::vector<Rational> coefficients(lhs.degree() + rhs.degree() + 1, rational(0));
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        for (std::size_t j = 0; j <= rhs.degree(); ++j)
            coefficients[i + j] += lhs.coefficient(i) * rhs.coefficient(j);
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] RationalPolynomial powerRationalPolynomial(
    RationalPolynomial base,
    std::uint64_t exponent) {
    RationalPolynomial result{{rational(1)}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0U)
            result = multiplyRationalPolynomials(result, base);
        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyRationalPolynomials(base, base);
    }
    return result;
}

[[nodiscard]] std::optional<RationalPolynomial> exactPolynomialPowerRoot(
    const RationalPolynomial& polynomial,
    std::uint64_t degree) {
    if (degree < 2 || polynomial.isZero() || polynomial.degree() % degree != 0)
        return std::nullopt;

    const std::size_t rootDegree = polynomial.degree() / degree;
    const auto leadingRoot = perfectRationalNthRoot(
        polynomial.coefficient(polynomial.degree()), degree);
    if (!leadingRoot || leadingRoot->isZero())
        return std::nullopt;

    std::vector<Rational> coefficients(rootDegree + 1, rational(0));
    coefficients[rootDegree] = *leadingRoot;
    const Rational linearScale = Rational{BigInt::fromUnsigned(degree)}
        * Rational{numeric::pow(leadingRoot->numerator(), degree - 1),
            numeric::pow(leadingRoot->denominator(), degree - 1)};

    // q(x)^k の上位係数を順に一致させる。x^((k-1)m+j) の係数では
    // 未知q_jは k*q_m^(k-1)*q_j と一次にしか現れないため，Q上でexactに決定できる。
    for (std::size_t offset = 0; offset < rootDegree; ++offset) {
        const std::size_t j = rootDegree - 1 - offset;
        RationalPolynomial partial{coefficients};
        const RationalPolynomial powered = powerRationalPolynomial(partial, degree);
        const std::size_t targetExponent = (degree - 1) * rootDegree + j;
        coefficients[j] = (polynomial.coefficient(targetExponent)
            - powered.coefficient(targetExponent)) / linearScale;
    }

    RationalPolynomial root{std::move(coefficients)};
    if (powerRationalPolynomial(root, degree).coefficients() != polynomial.coefficients())
        return std::nullopt;
    return root;
}

[[nodiscard]] std::optional<Expr> factorPerfectUnivariatePower(
    const RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    constexpr std::size_t maximumExponent = 64;
    const std::size_t upper = std::min(polynomial.degree(), maximumExponent);
    for (std::size_t exponent = upper; exponent >= 2; --exponent) {
        if (polynomial.degree() % exponent != 0)
            continue;
        const auto root = exactPolynomialPowerRoot(polynomial, exponent);
        if (!root)
            continue;
        return Expr::call(
            builtins.symbol(BuiltinId::Power),
            {polynomialToExpandedExpr(*root, variable, builtins),
                Expr{Number{BigInt::fromUnsigned(exponent)}}});
    }
    return std::nullopt;
}

[[nodiscard]] Monomial divideMonomial(
    const Monomial& value,
    const Monomial& divisor) {
    std::vector<MonomialFactor> factors;
    for (const MonomialFactor& factor : value.factors()) {
        const std::size_t divisorExponent = divisor.exponentOf(factor.variable);
        if (factor.exponent > divisorExponent)
            factors.push_back(MonomialFactor{factor.variable, factor.exponent - divisorExponent});
    }
    return Monomial{std::move(factors)};
}

[[nodiscard]] Monomial commonMonomial(const MultivariateRationalPolynomial& polynomial) {
    if (polynomial.isZero())
        return Monomial{};
    std::vector<MonomialFactor> common(
        polynomial.terms().front().monomial.factors().begin(),
        polynomial.terms().front().monomial.factors().end());
    for (const PolynomialTerm& term : polynomial.terms()) {
        for (MonomialFactor& factor : common)
            factor.exponent = std::min(
                factor.exponent, term.monomial.exponentOf(factor.variable));
    }
    return Monomial{std::move(common)};
}

[[nodiscard]] Rational rationalContent(const MultivariateRationalPolynomial& polynomial) {
    if (polynomial.isZero())
        return rational(1);
    BigInt numeratorGcd{0};
    BigInt denominatorLcm{1};
    for (const PolynomialTerm& term : polynomial.terms()) {
        numeratorGcd = numeratorGcd.isZero()
            ? term.coefficient.numerator().abs()
            : numeric::gcd(numeratorGcd, term.coefficient.numerator().abs());
        denominatorLcm = numeric::lcm(denominatorLcm, term.coefficient.denominator());
    }
    Rational content{numeratorGcd, denominatorLcm};
    if (polynomial.terms().front().coefficient.numerator().isNegative())
        content = -content;
    return content;
}

[[nodiscard]] MultivariateRationalPolynomial divideCommon(
    const MultivariateRationalPolynomial& polynomial,
    const Rational& content,
    const Monomial& monomial) {
    std::vector<PolynomialTerm> terms;
    terms.reserve(polynomial.termCount());
    for (const PolynomialTerm& term : polynomial.terms()) {
        terms.push_back(PolynomialTerm{
            divideMonomial(term.monomial, monomial),
            term.coefficient / content});
    }
    return MultivariateRationalPolynomial{std::move(terms)};
}

[[nodiscard]] Expr monomialExpr(
    const Monomial& monomial,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Expr> factors;
    for (const MonomialFactor& factor : monomial.factors()) {
        Expr value{factor.variable};
        if (factor.exponent > 1) {
            value = Expr::call(
                builtins.symbol(BuiltinId::Power),
                {std::move(value), Expr{Number{BigInt::parse(std::to_string(factor.exponent))}}});
        }
        factors.push_back(std::move(value));
    }
    if (factors.empty())
        return Expr{Number{BigInt{1}}};
    if (factors.size() == 1)
        return factors.front();
    return Expr::call(builtins.symbol(BuiltinId::Multiply), std::move(factors));
}

[[nodiscard]] std::optional<Expr> monomialRootExpr(
    const PolynomialTerm& term,
    unsigned degree,
    const evaluation::BuiltinRegistry& builtins) {
    std::optional<Rational> coefficientRoot;
    if (degree == 2)
        coefficientRoot = perfectRationalSquareRoot(term.coefficient);
    else if (degree == 3)
        coefficientRoot = perfectRationalCubeRoot(term.coefficient);
    else
        return std::nullopt;
    if (!coefficientRoot)
        return std::nullopt;

    std::vector<MonomialFactor> factors;
    for (const MonomialFactor& factor : term.monomial.factors()) {
        if (factor.exponent % degree != 0)
            return std::nullopt;
        if (factor.exponent / degree != 0)
            factors.push_back(MonomialFactor{factor.variable, factor.exponent / degree});
    }
    Expr atom = monomialExpr(Monomial{std::move(factors)}, builtins);
    return mathematics::scaleExactExpression(*coefficientRoot, atom, builtins);
}

[[nodiscard]] Expr addExpr(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins) {
    return Expr::call(builtins.symbol(BuiltinId::Add), {std::move(lhs), std::move(rhs)});
}

[[nodiscard]] Expr subtractExpr(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins) {
    return Expr::call(builtins.symbol(BuiltinId::Subtract), {std::move(lhs), std::move(rhs)});
}

[[nodiscard]] Expr multiplyExpr(
    std::vector<Expr> factors,
    const evaluation::BuiltinRegistry& builtins) {
    if (factors.empty())
        return Expr{Number{BigInt{1}}};
    if (factors.size() == 1)
        return factors.front();
    return Expr::call(builtins.symbol(BuiltinId::Multiply), std::move(factors));
}


struct ProductParts final {
    bool negative = false;
    std::vector<Expr> factors;
};

[[nodiscard]] ProductParts decomposeProduct(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    ProductParts result;
    std::vector<Expr> pending{expression};
    while (!pending.empty()) {
        Expr current = std::move(pending.back());
        pending.pop_back();
        if (isHead(current, builtins, BuiltinId::Negate)
            && current.asCall().arguments.size() == 1) {
            result.negative = !result.negative;
            pending.push_back(current.asCall().arguments.front());
            continue;
        }
        if (isHead(current, builtins, BuiltinId::Multiply)) {
            const auto& arguments = current.asCall().arguments;
            for (auto iterator = arguments.rbegin(); iterator != arguments.rend(); ++iterator)
                pending.push_back(*iterator);
            continue;
        }
        result.factors.push_back(std::move(current));
    }
    return result;
}

[[nodiscard]] std::vector<Expr> structuralCommonFactors(
    std::span<const Expr> terms,
    const evaluation::BuiltinRegistry& builtins) {
    if (terms.size() < 2)
        return {};

    ProductParts first = decomposeProduct(terms.front(), builtins);
    std::vector<Expr> common;
    for (const Expr& factor : first.factors)
        if (!factor.isNumber())
            common.push_back(factor);

    for (std::size_t termIndex = 1; termIndex < terms.size() && !common.empty(); ++termIndex) {
        ProductParts parts = decomposeProduct(terms[termIndex], builtins);
        std::vector<bool> used(parts.factors.size(), false);
        std::vector<Expr> retained;
        retained.reserve(common.size());
        for (const Expr& candidate : common) {
            for (std::size_t i = 0; i < parts.factors.size(); ++i) {
                if (used[i] || parts.factors[i].isNumber() || !(parts.factors[i] == candidate))
                    continue;
                used[i] = true;
                retained.push_back(candidate);
                break;
            }
        }
        common = std::move(retained);
    }
    return common;
}

[[nodiscard]] Expr removeStructuralFactors(
    const Expr& term,
    std::span<const Expr> common,
    const evaluation::BuiltinRegistry& builtins) {
    ProductParts parts = decomposeProduct(term, builtins);
    std::vector<bool> removed(parts.factors.size(), false);
    for (const Expr& candidate : common) {
        for (std::size_t i = 0; i < parts.factors.size(); ++i) {
            if (!removed[i] && parts.factors[i] == candidate) {
                removed[i] = true;
                break;
            }
        }
    }

    std::vector<Expr> remaining;
    for (std::size_t i = 0; i < parts.factors.size(); ++i)
        if (!removed[i])
            remaining.push_back(parts.factors[i]);
    Expr result = multiplyExpr(std::move(remaining), builtins);
    if (parts.negative)
        result = Expr::call(builtins.symbol(BuiltinId::Negate), {std::move(result)});
    return result;
}

[[nodiscard]] std::optional<Expr> factorStructuralCommonExpression(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Add))
        return std::nullopt;
    const auto& terms = expression.asCall().arguments;
    const std::vector<Expr> common = structuralCommonFactors(terms, builtins);
    if (common.empty())
        return std::nullopt;

    std::vector<Expr> reducedTerms;
    reducedTerms.reserve(terms.size());
    for (const Expr& term : terms)
        reducedTerms.push_back(removeStructuralFactors(term, common, builtins));
    Expr inner = simplify(
        Expr::call(builtins.symbol(BuiltinId::Add), std::move(reducedTerms)),
        builtins, mathematics, angles);
    inner = factorExpression(inner, builtins, mathematics, angles);

    std::vector<Expr> factors = common;
    factors.push_back(std::move(inner));
    return simplify(multiplyExpr(std::move(factors), builtins), builtins, mathematics, angles);
}


[[nodiscard]] std::optional<Expr> factorCommonExpressionCoefficient(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::vector<expression::Symbol> symbols = collectSymbols(expression);
    for (const expression::Symbol& variable : symbols) {
        const auto polynomial = toExpressionPolynomial(
            expression, variable, builtins, mathematics, angles);
        if (!polynomial || polynomial->degree() == 0)
            continue;

        std::vector<Expr> factoredCoefficients;
        factoredCoefficients.reserve(polynomial->coefficients().size());
        std::vector<Expr> nonZeroCoefficients;
        for (const Expr& coefficient : polynomial->coefficients()) {
            if (coefficient.isNumber() && coefficient.asNumber().isZero()) {
                factoredCoefficients.push_back(coefficient);
                continue;
            }
            Expr factored = factorExpression(coefficient, builtins, mathematics, angles);
            nonZeroCoefficients.push_back(factored);
            factoredCoefficients.push_back(std::move(factored));
        }
        if (nonZeroCoefficients.size() < 2)
            continue;

        const std::vector<Expr> common = structuralCommonFactors(
            nonZeroCoefficients, builtins);
        if (common.empty())
            continue;

        std::vector<Expr> reducedCoefficients;
        reducedCoefficients.reserve(factoredCoefficients.size());
        for (const Expr& coefficient : factoredCoefficients) {
            if (coefficient.isNumber() && coefficient.asNumber().isZero())
                reducedCoefficients.push_back(coefficient);
            else
                reducedCoefficients.push_back(simplify(
                    removeStructuralFactors(coefficient, common, builtins),
                    builtins, mathematics, angles));
        }

        Expr core = expressionPolynomialToCollectedExpr(
            ExpressionPolynomial{std::move(reducedCoefficients)},
            variable, builtins, mathematics, angles);
        core = factorExpression(core, builtins, mathematics, angles);
        std::vector<Expr> factors = common;
        factors.push_back(std::move(core));
        return simplify(
            multiplyExpr(std::move(factors), builtins),
            builtins, mathematics, angles);
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> factorDifferenceOfSquares(
    const MultivariateRationalPolynomial& polynomial,
    const evaluation::BuiltinRegistry& builtins) {
    if (polynomial.termCount() != 2)
        return std::nullopt;
    const PolynomialTerm* positive = nullptr;
    const PolynomialTerm* negative = nullptr;
    for (const PolynomialTerm& term : polynomial.terms()) {
        if (term.coefficient > rational(0))
            positive = &term;
        else if (term.coefficient < rational(0))
            negative = &term;
    }
    if (!positive || !negative)
        return std::nullopt;

    const auto left = monomialRootExpr(*positive, 2, builtins);
    PolynomialTerm positiveNegative{negative->monomial, -negative->coefficient};
    const auto right = monomialRootExpr(positiveNegative, 2, builtins);
    if (!left || !right)
        return std::nullopt;
    return multiplyExpr(
        {subtractExpr(*left, *right, builtins), addExpr(*left, *right, builtins)},
        builtins);
}

[[nodiscard]] bool sameMonomialProduct(
    const Monomial& lhs,
    const Monomial& rhs,
    const Monomial& expected) {
    std::vector<MonomialFactor> factors(lhs.factors().begin(), lhs.factors().end());
    factors.insert(factors.end(), rhs.factors().begin(), rhs.factors().end());
    return Monomial{std::move(factors)} == expected;
}

[[nodiscard]] std::optional<Expr> factorPerfectSquareTrinomial(
    const MultivariateRationalPolynomial& polynomial,
    const evaluation::BuiltinRegistry& builtins) {
    if (polynomial.termCount() != 3)
        return std::nullopt;

    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = i + 1; j < 3; ++j) {
            const auto a = monomialRootExpr(polynomial.terms()[i], 2, builtins);
            const auto b = monomialRootExpr(polynomial.terms()[j], 2, builtins);
            if (!a || !b)
                continue;
            const std::size_t k = 3 - i - j;
            const PolynomialTerm& cross = polynomial.terms()[k];
            const auto ai = perfectRationalSquareRoot(polynomial.terms()[i].coefficient);
            const auto bj = perfectRationalSquareRoot(polynomial.terms()[j].coefficient);
            if (!ai || !bj)
                continue;

            std::vector<MonomialFactor> af;
            for (const MonomialFactor& factor : polynomial.terms()[i].monomial.factors())
                af.push_back(MonomialFactor{factor.variable, factor.exponent / 2});
            std::vector<MonomialFactor> bf;
            for (const MonomialFactor& factor : polynomial.terms()[j].monomial.factors())
                bf.push_back(MonomialFactor{factor.variable, factor.exponent / 2});
            if (!sameMonomialProduct(Monomial{af}, Monomial{bf}, cross.monomial))
                continue;

            const Rational expected = rational(2) * *ai * *bj;
            std::optional<Expr> binomial;
            if (cross.coefficient == expected)
                binomial = addExpr(*a, *b, builtins);
            else if (cross.coefficient == -expected)
                binomial = subtractExpr(*a, *b, builtins);
            else
                continue;
            return Expr::call(
                builtins.symbol(BuiltinId::Power),
                {std::move(*binomial), Expr{Number{BigInt{2}}}});
        }
    }
    return std::nullopt;
}

[[nodiscard]] Expr squareExpr(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins) {
    if (value.isNumber() && value.asNumber().isReal()
        && value.asNumber().asReal().toRational() == rational(1))
        return value;
    return Expr::call(
        builtins.symbol(BuiltinId::Power),
        {value, Expr{Number{BigInt{2}}}});
}

[[nodiscard]] std::optional<Expr> factorSumOfCubes(
    const MultivariateRationalPolynomial& polynomial,
    const evaluation::BuiltinRegistry& builtins) {
    if (polynomial.termCount() != 2)
        return std::nullopt;

    PolynomialTerm leftTerm = polynomial.terms()[0];
    PolynomialTerm rightTerm = polynomial.terms()[1];
    const bool difference = (leftTerm.coefficient < rational(0))
        != (rightTerm.coefficient < rational(0));

    // 読みやすい a^3-b^3 形にするため、差の場合は正項をa、負項の絶対値をbにする。
    if (difference && leftTerm.coefficient < rational(0))
        std::swap(leftTerm, rightTerm);
    if (difference)
        rightTerm.coefficient = -rightTerm.coefficient;

    const auto a = monomialRootExpr(leftTerm, 3, builtins);
    const auto b = monomialRootExpr(rightTerm, 3, builtins);
    if (!a || !b)
        return std::nullopt;

    Expr first = difference
        ? subtractExpr(*a, *b, builtins)
        : addExpr(*a, *b, builtins);
    Expr a2 = squareExpr(*a, builtins);
    Expr b2 = squareExpr(*b, builtins);
    Expr ab = Expr::call(builtins.symbol(BuiltinId::Multiply), {*a, *b});
    Expr second = difference
        ? addExpr(addExpr(std::move(a2), std::move(ab), builtins), std::move(b2), builtins)
        : addExpr(subtractExpr(std::move(a2), std::move(ab), builtins), std::move(b2), builtins);
    return multiplyExpr({std::move(first), std::move(second)}, builtins);
}

[[nodiscard]] Expr primitiveLinearFactor(
    const expression::Symbol& variable,
    const Rational& root,
    const evaluation::BuiltinRegistry& builtins) {
    // root = p/q (gcd(p,q)=1, q>0) に対して q*x-p を返す。
    // x-p/q のまま因子化すると最終的なleading coefficientと分母が
    // 打ち消し合う表示になりやすいため、Q[x]上でもprimitiveな整数係数因子を優先する。
    const BigInt& numerator = root.numerator();
    const BigInt& denominator = root.denominator();
    Expr variableTerm{variable};
    if (!(denominator == BigInt{1})) {
        variableTerm = Expr::call(
            builtins.symbol(BuiltinId::Multiply),
            {Expr{Number{denominator}}, std::move(variableTerm)});
    }
    if (numerator.isZero())
        return variableTerm;
    if (numerator.isNegative())
        return Expr::call(
            builtins.symbol(BuiltinId::Add),
            {std::move(variableTerm), Expr{Number{-numerator}}});
    return Expr::call(
        builtins.symbol(BuiltinId::Subtract),
        {std::move(variableTerm), Expr{Number{numerator}}});
}

[[nodiscard]] std::optional<Expr> factorQuadraticInPower(
    const RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    const std::size_t degree = polynomial.degree();
    if (degree < 4 || degree % 2 != 0
        || polynomial.coefficient(degree) != rational(1))
        return std::nullopt;

    const std::size_t middle = degree / 2;
    for (std::size_t i = 1; i < degree; ++i)
        if (i != middle && !polynomial.coefficient(i).isZero())
            return std::nullopt;

    const Rational b = polynomial.coefficient(middle);
    const Rational c = polynomial.coefficient(0);
    const Rational discriminant = b * b - rational(4) * c;
    if (discriminant < rational(0))
        return std::nullopt;
    const auto numeratorRoot = numeric::integerSqrt(discriminant.numerator());
    const auto denominatorRoot = numeric::integerSqrt(discriminant.denominator());
    if (!numeratorRoot.remainder.isZero() || !denominatorRoot.remainder.isZero())
        return std::nullopt;

    const Rational sqrtDiscriminant{numeratorRoot.root, denominatorRoot.root};
    const Rational y1 = (-b - sqrtDiscriminant) / rational(2);
    const Rational y2 = (-b + sqrtDiscriminant) / rational(2);
    if (y1 == y2)
        return std::nullopt;

    auto power = [&] {
        if (middle == 1)
            return Expr{variable};
        return Expr::call(builtins.symbol(BuiltinId::Power), {
            Expr{variable}, Expr{Number{BigInt::fromUnsigned(middle)}}});
    };
    auto factor = [&](const Rational& root) {
        if (root.isZero())
            return power();
        return subtractExpr(power(), Expr{Number{root}}, builtins);
    };
    return multiplyExpr({factor(y1), factor(y2)}, builtins);
}

[[nodiscard]] std::optional<Expr> factorSparseCyclotomicTrinomial(
    const RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    const std::size_t degree = polynomial.degree();
    if (degree < 4 || degree % 4 != 0)
        return std::nullopt;
    const std::size_t middle = degree / 2;
    const std::size_t inner = degree / 4;
    if (!(polynomial.coefficient(degree) == rational(1))
        || !(polynomial.coefficient(middle) == rational(1))
        || !(polynomial.coefficient(0) == rational(1)))
        return std::nullopt;
    for (std::size_t i = 1; i < degree; ++i)
        if (i != middle && !polynomial.coefficient(i).isZero())
            return std::nullopt;

    auto powerExpr = [&](std::size_t exponent) {
        if (exponent == 1)
            return Expr{variable};
        return Expr::call(builtins.symbol(BuiltinId::Power), {
            Expr{variable}, Expr{Number{BigInt::fromUnsigned(exponent)}}});
    };
    Expr high = powerExpr(2 * inner);
    Expr low = powerExpr(inner);
    Expr plus = addExpr(addExpr(high, low, builtins), Expr{Number{BigInt{1}}}, builtins);
    Expr minus = addExpr(subtractExpr(powerExpr(2 * inner), powerExpr(inner), builtins),
        Expr{Number{BigInt{1}}}, builtins);
    return multiplyExpr({std::move(plus), std::move(minus)}, builtins);
}

[[nodiscard]] std::optional<Expr> factorUnivariate(
    RationalPolynomial polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> factors;
    bool foundFactor = false;
    Rational primitiveScale = rational(1);

    while (polynomial.degree() > 0) {
        const RationalRootSearchResult search = findRationalRoot(polynomial);
        if (!search.root)
            break;

        const Rational root = *search.root;
        std::size_t multiplicity = 0;
        while (true) {
            auto quotient = divideByLinearFactor(polynomial, root);
            if (!quotient)
                break;
            polynomial = std::move(*quotient);
            ++multiplicity;
            if (polynomial.degree() == 0)
                break;
        }
        if (multiplicity == 0)
            break;
        foundFactor = true;
        Expr factor = primitiveLinearFactor(variable, root, builtins);
        for (std::size_t i = 0; i < multiplicity; ++i)
            primitiveScale = primitiveScale / Rational{root.denominator()};
        if (multiplicity > 1) {
            factor = Expr::call(
                builtins.symbol(BuiltinId::Power),
                {std::move(factor), Expr{Number{BigInt::parse(std::to_string(multiplicity))}}});
        }
        factors.push_back(std::move(factor));
    }

    if (!foundFactor)
        return std::nullopt;

    if (polynomial.degree() == 0) {
        const Rational scalar = polynomial.coefficient(0) * primitiveScale;
        if (!(scalar == rational(1)))
            factors.insert(factors.begin(), Expr{Number{scalar}});
    }
    else {
        if (!(primitiveScale == rational(1)))
            factors.insert(factors.begin(), Expr{Number{primitiveScale}});
        factors.push_back(polynomialToExpandedExpr(polynomial, variable, builtins));
    }

    static_cast<void>(mathematics);
    static_cast<void>(angles);
    return multiplyExpr(std::move(factors), builtins);
}

[[nodiscard]] Expr collectExpressionRecursive(
    const Expr& expression,
    std::span<const expression::Symbol> variables,
    std::size_t variableIndex,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (variableIndex >= variables.size())
        return expression;

    const expression::Symbol& variable = variables[variableIndex];
    const auto polynomial = toExpressionPolynomial(
        expression, variable, builtins, mathematics, angles);
    if (!polynomial)
        return expression;

    std::vector<Expr> coefficients;
    coefficients.reserve(polynomial->coefficients().size());
    for (const Expr& coefficient : polynomial->coefficients())
        coefficients.push_back(collectExpressionRecursive(
            coefficient, variables, variableIndex + 1,
            builtins, mathematics, angles));
    return expressionPolynomialToCollectedExpr(
        ExpressionPolynomial{std::move(coefficients)},
        variable, builtins, mathematics, angles);
}

} // namespace

Expr expandExpression(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    AlgebraTransformOptions options) {
    const Expr expanded = expandRecursive(expression, builtins, mathematics, angles, options);
    // 純粋な有理係数多項式なら多変数Polynomialを一度通して項を標準順へ揃える。
    if (const auto polynomial = toMultivariateRationalPolynomial(expanded, builtins))
        return polynomialToExpandedExpr(*polynomial, builtins);
    return expanded;
}

Expr collectExpression(
    const Expr& expression,
    std::span<const expression::Symbol> variables,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Expr expanded = expandExpression(expression, builtins, mathematics, angles);
    if (variables.empty())
        return expanded;

    // 指定順に一変数ExpressionPolynomialを再帰的に作る。各段の係数は
    // 次の変数や sin/log などを含む任意のexact Exprでよいため、
    // collect[sin[z] x y + log[z] x + E, {x,y}] のような式も
    // 「xの係数がyの式」という自然な塔構造として保持できる。
    return collectExpressionRecursive(
        expanded, variables, 0, builtins, mathematics, angles);
}

Expr collectExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::vector<expression::Symbol> variables{variable};
    return collectExpression(expression, variables, builtins, mathematics, angles);
}

Expr factorExpression(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr simplified = simplify(expression, builtins, mathematics, angles);

    // Rational polynomialへ落ちない函数係数でも、各項に同じexact部分式が
    // 構造的に掛かっているなら安全に括り出せる。
    // 例: sin[y] x^2 + 2 sin[y] x + sin[y] -> sin[y] (x+1)^2。
    if (const auto structural = factorStructuralCommonExpression(
        simplified, builtins, mathematics, angles))
        return *structural;

    // 各項で係数の形が一度展開されていても、指定変数に関する
    // ExpressionPolynomialへ見直すと共通係数を回収できる。
    // 例: a x^2+b x^2+2 a x+2 b x+a+b -> (a+b)(x+1)^2。
    if (const auto coefficient = factorCommonExpressionCoefficient(
        simplified, builtins, mathematics, angles))
        return *coefficient;

    const auto polynomial = toMultivariateRationalPolynomial(simplified, builtins);
    if (!polynomial || polynomial->isZero() || polynomial->termCount() < 2)
        return simplified;

    const Rational content = rationalContent(*polynomial);
    const Monomial common = commonMonomial(*polynomial);
    MultivariateRationalPolynomial primitive = divideCommon(*polynomial, content, common);

    Expr core = polynomialToExpandedExpr(primitive, builtins);
    bool factoredCore = false;

    if (const auto square = factorPerfectSquareTrinomial(primitive, builtins)) {
        core = *square;
        factoredCore = true;
    }
    else if (const auto difference = factorDifferenceOfSquares(primitive, builtins)) {
        core = *difference;
        factoredCore = true;
    }
    else if (const auto cubes = factorSumOfCubes(primitive, builtins)) {
        core = *cubes;
        factoredCore = true;
    }
    else {
        const auto variables = primitive.variables();
        if (variables.size() == 1) {
            if (const auto univariate = toRationalPolynomial(core, variables.front(), builtins);
                univariate && univariate->degree() > 1) {
                if (const auto quadraticPower = factorQuadraticInPower(
                    *univariate, variables.front(), builtins)) {
                    core = *quadraticPower;
                    factoredCore = true;
                }
                else if (const auto sparse = factorSparseCyclotomicTrinomial(
                    *univariate, variables.front(), builtins)) {
                    core = *sparse;
                    factoredCore = true;
                }
                else if (const auto perfectPower = factorPerfectUnivariatePower(
                    *univariate, variables.front(), builtins)) {
                    core = *perfectPower;
                    factoredCore = true;
                }
                else if (const auto result = factorUnivariate(
                    *univariate, variables.front(), builtins, mathematics, angles)) {
                    core = *result;
                    factoredCore = true;
                }
            }
        }
    }

    // difference of squares/cubes等で積へ分かれた後は、各因子も同じfactor engineへ
    // 再帰的に渡す。次数が下がる因子だけを辿るため、x^6-1 のような式を
    // (x-1)(x+1)(x^2+x+1)(x^2-x+1) まで段階的に分解できる。
    if (factoredCore && isHead(core, builtins, BuiltinId::Multiply)) {
        std::vector<Expr> recursiveFactors;
        recursiveFactors.reserve(core.asCall().arguments.size());
        for (const Expr& factor : core.asCall().arguments) {
            // 一次因子は Q 上ですでに既約であり、さらに factorExpression へ渡すと
            // rational content の正規化と linear-root 分解が互いに形を変え、
            // 例: x-3/2 <-> (1/2)(2x-3) のような循環を作り得る。
            // 再帰は次数を厳密に下げられる非線形因子だけに限定する。
            const auto factorPolynomial = toMultivariateRationalPolynomial(factor, builtins);
            if (factorPolynomial && factorPolynomial->totalDegree() <= 1) {
                recursiveFactors.push_back(factor);
                continue;
            }
            recursiveFactors.push_back(factorExpression(
                factor, builtins, mathematics, angles));
        }
        core = multiplyExpr(std::move(recursiveFactors), builtins);
    }

    const bool hasContent = !(content == rational(1));
    const bool hasCommonMonomial = !common.isOne();
    if (!factoredCore && !hasContent && !hasCommonMonomial)
        return polynomialToExpandedExpr(*polynomial, builtins);

    std::vector<Expr> factors;
    if (hasContent)
        factors.emplace_back(Number{content});
    if (hasCommonMonomial)
        factors.push_back(monomialExpr(common, builtins));
    factors.push_back(std::move(core));
    return multiplyExpr(std::move(factors), builtins);
}

} // namespace mmcal::symbolic
