// 多項式表現と操作
#include "polynomial.hpp"

#include "linear_algebra/exact_lll.hpp"
#include "mathematics/exact_algebra.hpp"
#include "mathematics/definedness.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
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

[[nodiscard]] Rational zero() { return Rational{BigInt{0}}; }
[[nodiscard]] Rational one() { return Rational{BigInt{1}}; }

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return builtins.isCallTo(expression, id);
}

[[nodiscard]] std::optional<std::size_t> nonNegativeIntegerExponent(
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
    if (conversion.ec != std::errc{}
        || conversion.ptr != text.data() + text.size()
        || result > maximum)
        return std::nullopt;
    return result;
}

[[nodiscard]] bool symbolLess(
    const expression::Symbol& lhs,
    const expression::Symbol& rhs) noexcept {
    return lhs < rhs;
}

[[nodiscard]] bool monomialLess(const Monomial& lhs, const Monomial& rhs) {
    // 表示に都合のよい graded lexicographic order の逆順を作るため、ここでは「より高次数・名前の若い変数の指数が大きい方」を先とする。
    if (lhs.totalDegree() != rhs.totalDegree())
        return lhs.totalDegree() > rhs.totalDegree();

    const auto left = lhs.factors();
    const auto right = rhs.factors();
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < left.size() || j < right.size()) {
        if (i == left.size())
            return false;
        if (j == right.size())
            return true;

        if (left[i].variable == right[j].variable) {
            if (left[i].exponent != right[j].exponent)
                return left[i].exponent > right[j].exponent;
            ++i;
            ++j;
            continue;
        }
        if (symbolLess(left[i].variable, right[j].variable))
            return true;
        return false;
    }
    return false;
}

[[nodiscard]] Monomial multiplyMonomial(const Monomial& lhs, const Monomial& rhs) {
    std::vector<MonomialFactor> factors;
    factors.reserve(lhs.factors().size() + rhs.factors().size());
    factors.insert(factors.end(), lhs.factors().begin(), lhs.factors().end());
    factors.insert(factors.end(), rhs.factors().begin(), rhs.factors().end());
    return Monomial{std::move(factors)};
}

[[nodiscard]] MultivariateRationalPolynomial add(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs) {
    std::vector<PolynomialTerm> terms;
    terms.reserve(lhs.termCount() + rhs.termCount());
    terms.insert(terms.end(), lhs.terms().begin(), lhs.terms().end());
    terms.insert(terms.end(), rhs.terms().begin(), rhs.terms().end());
    return MultivariateRationalPolynomial{std::move(terms)};
}

[[nodiscard]] MultivariateRationalPolynomial negate(
    const MultivariateRationalPolynomial& value) {
    std::vector<PolynomialTerm> terms(value.terms().begin(), value.terms().end());
    for (PolynomialTerm& term : terms)
        term.coefficient = -term.coefficient;
    return MultivariateRationalPolynomial{std::move(terms)};
}

[[nodiscard]] std::optional<MultivariateRationalPolynomial> multiply(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs,
    const PolynomialConversionOptions& options) {
    if (lhs.isZero() || rhs.isZero())
        return MultivariateRationalPolynomial{};
    if (lhs.termCount() > options.maximumTerms / rhs.termCount())
        return std::nullopt;

    std::vector<PolynomialTerm> terms;
    terms.reserve(lhs.termCount() * rhs.termCount());
    for (const PolynomialTerm& left : lhs.terms()) {
        for (const PolynomialTerm& right : rhs.terms()) {
            Monomial monomial = multiplyMonomial(left.monomial, right.monomial);
            if (monomial.totalDegree() > options.maximumDegree)
                return std::nullopt;
            terms.push_back(PolynomialTerm{
                std::move(monomial), left.coefficient * right.coefficient});
        }
    }
    MultivariateRationalPolynomial result{std::move(terms)};
    if (result.termCount() > options.maximumTerms)
        return std::nullopt;
    return result;
}

[[nodiscard]] std::optional<MultivariateRationalPolynomial> power(
    MultivariateRationalPolynomial base,
    std::size_t exponent,
    const PolynomialConversionOptions& options) {
    MultivariateRationalPolynomial result{{PolynomialTerm{Monomial{}, one()}}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0U) {
            auto product = multiply(result, base, options);
            if (!product)
                return std::nullopt;
            result = std::move(*product);
        }
        exponent >>= 1U;
        if (exponent == 0)
            break;
        auto square = multiply(base, base, options);
        if (!square)
            return std::nullopt;
        base = std::move(*square);
    }
    return result;
}

[[nodiscard]] std::optional<MultivariateRationalPolynomial> convert(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const PolynomialConversionOptions& options) {
    if (expression.isNumber()) {
        if (!expression.asNumber().isReal())
            return std::nullopt;
        return MultivariateRationalPolynomial{{PolynomialTerm{
            Monomial{}, expression.asNumber().asReal().toRational()}}};
    }

    if (expression.isSymbol()) {
        return MultivariateRationalPolynomial{{PolynomialTerm{
            Monomial{{MonomialFactor{expression.asSymbol(), 1}}}, one()}}};
    }

    if (!expression.isCall())
        return std::nullopt;

    const auto& arguments = expression.asCall().arguments;
    if (isHead(expression, builtins, BuiltinId::Add)) {
        MultivariateRationalPolynomial result;
        for (const Expr& argument : arguments) {
            auto term = convert(argument, builtins, options);
            if (!term)
                return std::nullopt;
            result = add(result, *term);
            if (result.termCount() > options.maximumTerms)
                return std::nullopt;
        }
        return result;
    }

    if (isHead(expression, builtins, BuiltinId::Subtract) && arguments.size() == 2) {
        auto lhs = convert(arguments[0], builtins, options);
        auto rhs = convert(arguments[1], builtins, options);
        if (!lhs || !rhs)
            return std::nullopt;
        return add(*lhs, negate(*rhs));
    }

    if (isHead(expression, builtins, BuiltinId::Negate) && arguments.size() == 1) {
        auto operand = convert(arguments[0], builtins, options);
        return operand
            ? std::optional<MultivariateRationalPolynomial>{negate(*operand)}
            : std::nullopt;
    }

    if (isHead(expression, builtins, BuiltinId::Multiply)) {
        MultivariateRationalPolynomial result{{PolynomialTerm{Monomial{}, one()}}};
        for (const Expr& argument : arguments) {
            auto factor = convert(argument, builtins, options);
            if (!factor)
                return std::nullopt;
            auto product = multiply(result, *factor, options);
            if (!product)
                return std::nullopt;
            result = std::move(*product);
        }
        return result;
    }

    if (isHead(expression, builtins, BuiltinId::Divide) && arguments.size() == 2) {
        auto numerator = convert(arguments[0], builtins, options);
        if (!numerator || !arguments[1].isNumber()
            || !arguments[1].asNumber().isReal())
            return std::nullopt;
        const Rational denominator = arguments[1].asNumber().asReal().toRational();
        if (denominator.isZero())
            return std::nullopt;
        std::vector<PolynomialTerm> terms(numerator->terms().begin(), numerator->terms().end());
        for (PolynomialTerm& term : terms)
            term.coefficient /= denominator;
        return MultivariateRationalPolynomial{std::move(terms)};
    }

    if (isHead(expression, builtins, BuiltinId::Power) && arguments.size() == 2) {
        auto base = convert(arguments[0], builtins, options);
        const auto exponent = nonNegativeIntegerExponent(arguments[1], options.maximumDegree);
        if (!base || !exponent)
            return std::nullopt;
        if (!base->isZero() && *exponent != 0
            && base->totalDegree() > options.maximumDegree / *exponent)
            return std::nullopt;
        return power(std::move(*base), *exponent, options);
    }

    return std::nullopt;
}

[[nodiscard]] Expr integerExpr(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr rationalExpr(const Rational& value) {
    return Expr{Number{value}};
}

[[nodiscard]] Expr monomialExpr(
    const Monomial& monomial,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Expr> factors;
    factors.reserve(monomial.factors().size());
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
        return integerExpr(1);
    if (factors.size() == 1)
        return factors.front();
    return Expr::call(builtins.symbol(BuiltinId::Multiply), std::move(factors));
}

[[nodiscard]] Expr termExpr(
    const PolynomialTerm& term,
    const evaluation::BuiltinRegistry& builtins) {
    if (term.monomial.isOne())
        return rationalExpr(term.coefficient);
    return mathematics::scaleExactExpression(
        term.coefficient, monomialExpr(term.monomial, builtins), builtins);
}

[[nodiscard]] Expr termsToExpr(
    std::vector<PolynomialTerm> terms,
    const evaluation::BuiltinRegistry& builtins) {
    if (terms.empty())
        return integerExpr(0);
    std::sort(terms.begin(), terms.end(), [](const PolynomialTerm& lhs, const PolynomialTerm& rhs) {
        return monomialLess(lhs.monomial, rhs.monomial);
    });
    std::vector<Expr> expressions;
    expressions.reserve(terms.size());
    for (const PolynomialTerm& term : terms)
        expressions.push_back(termExpr(term, builtins));
    if (expressions.size() == 1)
        return expressions.front();
    return Expr::call(builtins.symbol(BuiltinId::Add), std::move(expressions));
}

[[nodiscard]] Expr variablePower(
    const expression::Symbol& variable,
    std::size_t exponent,
    const evaluation::BuiltinRegistry& builtins) {
    if (exponent == 0)
        return integerExpr(1);
    Expr result{variable};
    if (exponent > 1) {
        result = Expr::call(
            builtins.symbol(BuiltinId::Power),
            {std::move(result), Expr{Number{BigInt::parse(std::to_string(exponent))}}});
    }
    return result;
}

[[nodiscard]] Expr collectRecursive(
    std::span<const PolynomialTerm> terms,
    std::span<const expression::Symbol> variables,
    std::size_t variableIndex,
    const evaluation::BuiltinRegistry& builtins) {
    if (terms.empty())
        return integerExpr(0);
    if (variableIndex >= variables.size())
        return termsToExpr(std::vector<PolynomialTerm>{terms.begin(), terms.end()}, builtins);

    const expression::Symbol& variable = variables[variableIndex];
    std::size_t maximumExponent = 0;
    for (const PolynomialTerm& term : terms)
        maximumExponent = std::max(maximumExponent, term.monomial.exponentOf(variable));

    std::vector<Expr> resultTerms;
    for (std::size_t reverse = 0; reverse <= maximumExponent; ++reverse) {
        const std::size_t exponent = maximumExponent - reverse;
        std::vector<PolynomialTerm> coefficientTerms;
        for (const PolynomialTerm& term : terms) {
            if (term.monomial.exponentOf(variable) != exponent)
                continue;
            coefficientTerms.push_back(PolynomialTerm{
                term.monomial.without(variable), term.coefficient});
        }
        if (coefficientTerms.empty())
            continue;

        const Expr coefficient = collectRecursive(
            coefficientTerms, variables, variableIndex + 1, builtins);
        if (exponent == 0) {
            resultTerms.push_back(coefficient);
            continue;
        }

        Expr powerExpr = variablePower(variable, exponent, builtins);
        if (coefficient.isNumber() && coefficient.asNumber().isReal()
            && coefficient.asNumber().asReal().toRational() == one()) {
            resultTerms.push_back(std::move(powerExpr));
            continue;
        }
        resultTerms.push_back(Expr::call(
            builtins.symbol(BuiltinId::Multiply),
            {coefficient, std::move(powerExpr)}));
    }

    if (resultTerms.empty())
        return integerExpr(0);
    if (resultTerms.size() == 1)
        return resultTerms.front();
    return Expr::call(builtins.symbol(BuiltinId::Add), std::move(resultTerms));
}


[[nodiscard]] Expr zeroExpr() { return Expr{Number{BigInt{0}}}; }
[[nodiscard]] Expr oneExpr() { return Expr{Number{BigInt{1}}}; }

[[nodiscard]] bool isExactZeroExpr(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] bool isExactOneExpr(const Expr& expression) {
    return expression.isNumber()
        && expression.asNumber().isReal()
        && expression.asNumber().asReal().isInteger()
        && expression.asNumber().asReal().asInteger() == BigInt{1};
}

[[nodiscard]] bool containsSymbolImpl(
    const Expr& expression,
    const expression::Symbol& symbol) {
    std::vector<Expr> pending{expression};
    while (!pending.empty()) {
        Expr current = std::move(pending.back());
        pending.pop_back();
        if (current.isSymbol() && current.asSymbol() == symbol)
            return true;
        if (current.isCall()) {
            for (const Expr& argument : current.asCall().arguments)
                pending.push_back(argument);
        }
        else if (current.isArray()
            && current.asArray().storageKind() == expression::ArrayStorageKind::Generic) {
            for (const Expr& element : current.asArray().storedExpressions())
                pending.push_back(element);
        }
    }
    return false;
}

[[nodiscard]] bool scalarCoefficientCandidate(const Expr& expression) noexcept {
    return expression.isNumber() || expression.isSymbol() || expression.isCall();
}

[[nodiscard]] Expr simplifyCoefficient(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles});
}

[[nodiscard]] Expr coefficientNegate(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (isExactZeroExpr(value))
        return zeroExpr();
    return simplifyCoefficient(
        Expr::call(builtins.symbol(BuiltinId::Negate), {value}),
        builtins, mathematics, angles);
}

[[nodiscard]] Expr coefficientAdd(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (isExactZeroExpr(lhs))
        return rhs;
    if (isExactZeroExpr(rhs))
        return lhs;
    return simplifyCoefficient(
        Expr::call(builtins.symbol(BuiltinId::Add), {lhs, rhs}),
        builtins, mathematics, angles);
}

[[nodiscard]] Expr coefficientMultiply(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (isExactZeroExpr(lhs) || isExactZeroExpr(rhs))
        return zeroExpr();
    if (isExactOneExpr(lhs))
        return rhs;
    if (isExactOneExpr(rhs))
        return lhs;
    return simplifyCoefficient(
        Expr::call(builtins.symbol(BuiltinId::Multiply), {lhs, rhs}),
        builtins, mathematics, angles);
}

[[nodiscard]] Expr coefficientDivide(
    const Expr& numerator,
    const Expr& denominator,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (isExactZeroExpr(numerator))
        return zeroExpr();
    if (isExactOneExpr(denominator))
        return numerator;
    return simplifyCoefficient(
        Expr::call(builtins.symbol(BuiltinId::Divide), {numerator, denominator}),
        builtins, mathematics, angles);
}

[[nodiscard]] ExpressionPolynomial expressionPolynomialAdd(
    const ExpressionPolynomial& lhs,
    const ExpressionPolynomial& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::size_t size = std::max(lhs.degree(), rhs.degree()) + 1;
    std::vector<Expr> coefficients;
    coefficients.reserve(size);
    for (std::size_t exponent = 0; exponent < size; ++exponent)
        coefficients.push_back(coefficientAdd(
            lhs.coefficient(exponent), rhs.coefficient(exponent),
            builtins, mathematics, angles));
    return ExpressionPolynomial{std::move(coefficients)};
}

[[nodiscard]] ExpressionPolynomial expressionPolynomialNegate(
    const ExpressionPolynomial& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> coefficients;
    coefficients.reserve(value.coefficients().size());
    for (const Expr& coefficient : value.coefficients())
        coefficients.push_back(coefficientNegate(
            coefficient, builtins, mathematics, angles));
    return ExpressionPolynomial{std::move(coefficients)};
}

[[nodiscard]] std::optional<ExpressionPolynomial> expressionPolynomialMultiply(
    const ExpressionPolynomial& lhs,
    const ExpressionPolynomial& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const PolynomialConversionOptions& options) {
    if (lhs.isZero() || rhs.isZero())
        return ExpressionPolynomial{{zeroExpr()}};
    if (lhs.degree() > options.maximumDegree - std::min(rhs.degree(), options.maximumDegree))
        return std::nullopt;
    const std::size_t degree = lhs.degree() + rhs.degree();
    if (degree > options.maximumDegree)
        return std::nullopt;

    std::vector<Expr> coefficients(degree + 1, zeroExpr());
    for (std::size_t i = 0; i <= lhs.degree(); ++i) {
        for (std::size_t j = 0; j <= rhs.degree(); ++j) {
            Expr product = coefficientMultiply(
                lhs.coefficient(i), rhs.coefficient(j), builtins, mathematics, angles);
            coefficients[i + j] = coefficientAdd(
                coefficients[i + j], product, builtins, mathematics, angles);
        }
    }
    return ExpressionPolynomial{std::move(coefficients)};
}

[[nodiscard]] std::optional<ExpressionPolynomial> expressionPolynomialPower(
    ExpressionPolynomial base,
    std::size_t exponent,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const PolynomialConversionOptions& options) {
    ExpressionPolynomial result{{oneExpr()}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0) {
            auto product = expressionPolynomialMultiply(
                result, base, builtins, mathematics, angles, options);
            if (!product)
                return std::nullopt;
            result = std::move(*product);
        }
        exponent >>= 1U;
        if (exponent == 0)
            break;
        auto square = expressionPolynomialMultiply(
            base, base, builtins, mathematics, angles, options);
        if (!square)
            return std::nullopt;
        base = std::move(*square);
    }
    return result;
}


[[nodiscard]] std::optional<ExpressionPolynomial> convertExpressionPolynomial(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const PolynomialConversionOptions& options,
    mathematics::AssumptionSet* domainConditions) {
    if (!containsSymbolImpl(expression, variable)) {
        if (!scalarCoefficientCandidate(expression))
            return std::nullopt;
        if (domainConditions) {
            const auto required = mathematics::expressionDomainConditions(
                expression, builtins, mathematics);
            if (!required)
                return std::nullopt;
            for (const mathematics::Predicate& predicate : required->predicates())
                domainConditions->add(predicate);
        }
        return ExpressionPolynomial{{simplifyCoefficient(
            expression, builtins, mathematics, angles)}};
    }

    if (expression.isSymbol()) {
        if (!(expression.asSymbol() == variable))
            return std::nullopt;
        return ExpressionPolynomial{{zeroExpr(), oneExpr()}};
    }
    if (!expression.isCall())
        return std::nullopt;

    const auto& arguments = expression.asCall().arguments;
    if (isHead(expression, builtins, BuiltinId::Add)) {
        ExpressionPolynomial result{{zeroExpr()}};
        for (const Expr& argument : arguments) {
            auto term = convertExpressionPolynomial(
                argument, variable, builtins, mathematics, angles, options, domainConditions);
            if (!term)
                return std::nullopt;
            result = expressionPolynomialAdd(result, *term, builtins, mathematics, angles);
        }
        return result;
    }

    if (isHead(expression, builtins, BuiltinId::Subtract) && arguments.size() == 2) {
        auto lhs = convertExpressionPolynomial(
            arguments[0], variable, builtins, mathematics, angles, options, domainConditions);
        auto rhs = convertExpressionPolynomial(
            arguments[1], variable, builtins, mathematics, angles, options, domainConditions);
        if (!lhs || !rhs)
            return std::nullopt;
        return expressionPolynomialAdd(
            *lhs,
            expressionPolynomialNegate(*rhs, builtins, mathematics, angles),
            builtins, mathematics, angles);
    }

    if (isHead(expression, builtins, BuiltinId::Negate) && arguments.size() == 1) {
        auto operand = convertExpressionPolynomial(
            arguments[0], variable, builtins, mathematics, angles, options, domainConditions);
        if (!operand)
            return std::nullopt;
        return expressionPolynomialNegate(*operand, builtins, mathematics, angles);
    }

    if (isHead(expression, builtins, BuiltinId::Multiply)) {
        ExpressionPolynomial result{{oneExpr()}};
        for (const Expr& argument : arguments) {
            auto factor = convertExpressionPolynomial(
                argument, variable, builtins, mathematics, angles, options, domainConditions);
            if (!factor)
                return std::nullopt;
            auto product = expressionPolynomialMultiply(
                result, *factor, builtins, mathematics, angles, options);
            if (!product)
                return std::nullopt;
            result = std::move(*product);
        }
        return result;
    }

    if (isHead(expression, builtins, BuiltinId::Divide) && arguments.size() == 2) {
        if (containsSymbolImpl(arguments[1], variable)
            || !scalarCoefficientCandidate(arguments[1]))
            return std::nullopt;
        auto numerator = convertExpressionPolynomial(
            arguments[0], variable, builtins, mathematics, angles, options, domainConditions);
        if (!numerator)
            return std::nullopt;
        const Expr denominator = simplifyCoefficient(
            arguments[1], builtins, mathematics, angles);
        if (isExactZeroExpr(denominator))
            return std::nullopt;

        const mathematics::AssumptionSet assumptions;
        const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
        if (knowledge.prove(mathematics::relation(
                mathematics::RelationKind::NotEqual, denominator, zeroExpr()))
            != mathematics::TruthValue::True) {
            if (!domainConditions)
                return std::nullopt;
            domainConditions->add(mathematics::relation(
                mathematics::RelationKind::NotEqual, denominator, zeroExpr()));
        }
        std::vector<Expr> coefficients;
        coefficients.reserve(numerator->coefficients().size());
        for (const Expr& coefficient : numerator->coefficients())
            coefficients.push_back(coefficientDivide(
                coefficient, denominator, builtins, mathematics, angles));
        return ExpressionPolynomial{std::move(coefficients)};
    }

    if (isHead(expression, builtins, BuiltinId::Power) && arguments.size() == 2) {
        const auto exponent = nonNegativeIntegerExponent(arguments[1], options.maximumDegree);
        if (!exponent)
            return std::nullopt;
        auto base = convertExpressionPolynomial(
            arguments[0], variable, builtins, mathematics, angles, options, domainConditions);
        if (!base)
            return std::nullopt;
        if (!base->isZero() && *exponent != 0
            && base->degree() > options.maximumDegree / *exponent)
            return std::nullopt;
        return expressionPolynomialPower(
            std::move(*base), *exponent, builtins, mathematics, angles, options);
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<std::uint64_t> toUint64(const BigInt& value) {
    if (value.isNegative())
        return std::nullopt;
    const std::string text = value.toString();
    std::uint64_t result = 0;
    const auto conversion = std::from_chars(text.data(), text.data() + text.size(), result);
    if (conversion.ec != std::errc{} || conversion.ptr != text.data() + text.size())
        return std::nullopt;
    return result;
}

struct DivisorList final {
    std::vector<std::uint64_t> values;
    bool complete = false;
};

[[nodiscard]] DivisorList positiveDivisors(
    const BigInt& value,
    std::size_t maximumTrialDivisor) {
    const BigInt magnitude = value.abs();
    if (magnitude.isZero())
        return DivisorList{{0}, true};
    const auto converted = toUint64(magnitude);
    if (!converted)
        return {};

    const std::uint64_t n = *converted;
    std::vector<std::uint64_t> low;
    std::vector<std::uint64_t> high;
    bool complete = true;
    std::uint64_t divisor = 1;
    for (; divisor <= n / divisor; ++divisor) {
        if (divisor > maximumTrialDivisor) {
            complete = false;
            break;
        }
        if (n % divisor != 0)
            continue;
        low.push_back(divisor);
        const std::uint64_t pair = n / divisor;
        if (pair != divisor)
            high.push_back(pair);
    }
    std::reverse(high.begin(), high.end());
    low.insert(low.end(), high.begin(), high.end());
    return DivisorList{std::move(low), complete};
}

[[nodiscard]] std::vector<BigInt> primitiveIntegerCoefficients(
    const RationalPolynomial& polynomial) {
    BigInt commonDenominator{1};
    for (const Rational& coefficient : polynomial.coefficients())
        commonDenominator = numeric::lcm(commonDenominator, coefficient.denominator());

    std::vector<BigInt> integers;
    integers.reserve(polynomial.coefficients().size());
    BigInt commonNumerator{0};
    for (const Rational& coefficient : polynomial.coefficients()) {
        BigInt value = coefficient.numerator()
            * (commonDenominator / coefficient.denominator());
        integers.push_back(value);
        commonNumerator = commonNumerator.isZero()
            ? value.abs()
            : numeric::gcd(commonNumerator, value.abs());
    }
    if (!commonNumerator.isZero() && !(commonNumerator == BigInt{1}))
        for (BigInt& value : integers)
            value /= commonNumerator;
    return integers;
}

} // namespace

RationalPolynomial::RationalPolynomial()
    : coefficients_{zero()} {}

RationalPolynomial::RationalPolynomial(std::vector<Rational> coefficients)
    : coefficients_(std::move(coefficients)) {
    if (coefficients_.empty())
        coefficients_.push_back(zero());
    normalize();
}

void RationalPolynomial::normalize() {
    while (coefficients_.size() > 1 && coefficients_.back().isZero())
        coefficients_.pop_back();
}

bool RationalPolynomial::isZero() const noexcept {
    return coefficients_.size() == 1 && coefficients_.front().isZero();
}

std::size_t RationalPolynomial::degree() const noexcept {
    return coefficients_.size() - 1;
}

const Rational& RationalPolynomial::coefficient(std::size_t exponent) const noexcept {
    static const Rational zeroValue{BigInt{0}};
    return exponent < coefficients_.size() ? coefficients_[exponent] : zeroValue;
}

const std::vector<Rational>& RationalPolynomial::coefficients() const noexcept {
    return coefficients_;
}


ExpressionPolynomial::ExpressionPolynomial(std::vector<Expr> coefficients)
    : coefficients_(std::move(coefficients)) {
    if (coefficients_.empty())
        coefficients_.push_back(zeroExpr());
    normalize();
}

void ExpressionPolynomial::normalize() {
    while (coefficients_.size() > 1 && isExactZeroExpr(coefficients_.back()))
        coefficients_.pop_back();
}

bool ExpressionPolynomial::isZero() const noexcept {
    return coefficients_.size() == 1 && isExactZeroExpr(coefficients_.front());
}

std::size_t ExpressionPolynomial::degree() const noexcept {
    return coefficients_.size() - 1;
}

const Expr& ExpressionPolynomial::coefficient(std::size_t exponent) const noexcept {
    static const Expr zeroValue{Number{BigInt{0}}};
    return exponent < coefficients_.size() ? coefficients_[exponent] : zeroValue;
}

std::span<const Expr> ExpressionPolynomial::coefficients() const noexcept {
    return coefficients_;
}

Monomial::Monomial(std::vector<MonomialFactor> factors) {
    factors.erase(
        std::remove_if(factors.begin(), factors.end(), [](const MonomialFactor& factor) {
            return !factor.variable.valid() || factor.exponent == 0;
        }),
        factors.end());
    std::sort(factors.begin(), factors.end(), [](const MonomialFactor& lhs, const MonomialFactor& rhs) {
        return symbolLess(lhs.variable, rhs.variable);
    });

    for (const MonomialFactor& factor : factors) {
        if (!factors_.empty() && factors_.back().variable == factor.variable) {
            if (factors_.back().exponent > std::numeric_limits<std::size_t>::max() - factor.exponent)
                throw std::overflow_error("Monomial exponent overflow");
            factors_.back().exponent += factor.exponent;
        }
        else
            factors_.push_back(factor);
    }
}

bool Monomial::isOne() const noexcept { return factors_.empty(); }

std::size_t Monomial::totalDegree() const noexcept {
    std::size_t result = 0;
    for (const MonomialFactor& factor : factors_) {
        if (result > std::numeric_limits<std::size_t>::max() - factor.exponent)
            return std::numeric_limits<std::size_t>::max();
        result += factor.exponent;
    }
    return result;
}

std::size_t Monomial::exponentOf(const expression::Symbol& variable) const noexcept {
    for (const MonomialFactor& factor : factors_)
        if (factor.variable == variable)
            return factor.exponent;
    return 0;
}

std::span<const MonomialFactor> Monomial::factors() const noexcept { return factors_; }

Monomial Monomial::without(const expression::Symbol& variable) const {
    std::vector<MonomialFactor> factors;
    factors.reserve(factors_.size());
    for (const MonomialFactor& factor : factors_)
        if (!(factor.variable == variable))
            factors.push_back(factor);
    return Monomial{std::move(factors)};
}

std::optional<MonomialOrder> parseMonomialOrder(std::string_view name) noexcept {
    if (name == "Lex" || name == "lex")
        return MonomialOrder::Lex;
    if (name == "GrLex" || name == "grlex")
        return MonomialOrder::GrLex;
    if (name == "GrevLex" || name == "grevlex")
        return MonomialOrder::GrevLex;
    return std::nullopt;
}

std::string_view monomialOrderName(MonomialOrder order) noexcept {
    switch (order) {
    case MonomialOrder::Lex: return "Lex";
    case MonomialOrder::GrLex: return "GrLex";
    case MonomialOrder::GrevLex: return "GrevLex";
    }
    return "GrevLex";
}

PolynomialRing::PolynomialRing(
    std::vector<expression::Symbol> variables,
    MonomialOrder order)
    : variables_(std::move(variables)), order_(order) {
    for (std::size_t i = 0; i < variables_.size(); ++i) {
        if (!variables_[i].valid())
            throw std::invalid_argument("Polynomial ring contains an invalid variable");
        for (std::size_t j = 0; j < i; ++j)
            if (variables_[i] == variables_[j])
                throw std::invalid_argument("Polynomial ring variables must be unique");
    }
}

std::span<const expression::Symbol> PolynomialRing::variables() const noexcept {
    return variables_;
}

MonomialOrder PolynomialRing::order() const noexcept { return order_; }

bool PolynomialRing::contains(const expression::Symbol& variable) const noexcept {
    return std::find(variables_.begin(), variables_.end(), variable) != variables_.end();
}

bool PolynomialRing::contains(const Monomial& monomial) const noexcept {
    for (const MonomialFactor& factor : monomial.factors())
        if (!contains(factor.variable))
            return false;
    return true;
}

int PolynomialRing::compare(const Monomial& lhs, const Monomial& rhs) const {
    if (!contains(lhs) || !contains(rhs))
        throw std::invalid_argument("Monomial does not belong to the polynomial ring");

    const auto lexCompare = [&]() -> int {
        for (const expression::Symbol& variable : variables_) {
            const std::size_t left = lhs.exponentOf(variable);
            const std::size_t right = rhs.exponentOf(variable);
            if (left != right)
                return left < right ? -1 : 1;
        }
        return 0;
    };

    if (order_ == MonomialOrder::Lex)
        return lexCompare();

    const std::size_t leftDegree = lhs.totalDegree();
    const std::size_t rightDegree = rhs.totalDegree();
    if (leftDegree != rightDegree)
        return leftDegree < rightDegree ? -1 : 1;
    if (order_ == MonomialOrder::GrLex)
        return lexCompare();

    // graded reverse lexicographic: total degree equalなら、後ろから最初に
    // 異なる指数が小さい側を大きいmonomialとする。
    for (auto iterator = variables_.rbegin(); iterator != variables_.rend(); ++iterator) {
        const std::size_t left = lhs.exponentOf(*iterator);
        const std::size_t right = rhs.exponentOf(*iterator);
        if (left != right)
            return left < right ? 1 : -1;
    }
    return 0;
}

bool monomialDivides(const Monomial& divisor, const Monomial& dividend) noexcept {
    for (const MonomialFactor& factor : divisor.factors())
        if (factor.exponent > dividend.exponentOf(factor.variable))
            return false;
    return true;
}

Monomial multiplyMonomials(const Monomial& lhs, const Monomial& rhs) {
    return multiplyMonomial(lhs, rhs);
}

Monomial leastCommonMultiple(const Monomial& lhs, const Monomial& rhs) {
    std::vector<MonomialFactor> factors;
    factors.reserve(lhs.factors().size() + rhs.factors().size());
    for (const MonomialFactor& factor : lhs.factors())
        factors.push_back(factor);
    for (const MonomialFactor& factor : rhs.factors()) {
        const auto existing = std::find_if(factors.begin(), factors.end(), [&](const MonomialFactor& current) {
            return current.variable == factor.variable;
        });
        if (existing == factors.end())
            factors.push_back(factor);
        else
            existing->exponent = std::max(existing->exponent, factor.exponent);
    }
    return Monomial{std::move(factors)};
}

std::optional<Monomial> divideMonomials(
    const Monomial& dividend,
    const Monomial& divisor) {
    if (!monomialDivides(divisor, dividend))
        return std::nullopt;
    std::vector<MonomialFactor> factors;
    factors.reserve(dividend.factors().size());
    for (const MonomialFactor& factor : dividend.factors()) {
        const std::size_t divisorExponent = divisor.exponentOf(factor.variable);
        if (factor.exponent > divisorExponent)
            factors.push_back(MonomialFactor{
                factor.variable, factor.exponent - divisorExponent});
    }
    return Monomial{std::move(factors)};
}

MultivariateRationalPolynomial::MultivariateRationalPolynomial() = default;

MultivariateRationalPolynomial::MultivariateRationalPolynomial(std::vector<PolynomialTerm> terms)
    : terms_(std::move(terms)) {
    normalize();
}

void MultivariateRationalPolynomial::normalize() {
    std::vector<PolynomialTerm> normalized;
    for (PolynomialTerm& term : terms_) {
        if (term.coefficient.isZero())
            continue;
        auto existing = std::find_if(normalized.begin(), normalized.end(), [&](const PolynomialTerm& candidate) {
            return candidate.monomial == term.monomial;
        });
        if (existing == normalized.end())
            normalized.push_back(std::move(term));
        else
            existing->coefficient += term.coefficient;
    }
    normalized.erase(
        std::remove_if(normalized.begin(), normalized.end(), [](const PolynomialTerm& term) {
            return term.coefficient.isZero();
        }),
        normalized.end());
    std::sort(normalized.begin(), normalized.end(), [](const PolynomialTerm& lhs, const PolynomialTerm& rhs) {
        return monomialLess(lhs.monomial, rhs.monomial);
    });
    terms_ = std::move(normalized);
}

bool MultivariateRationalPolynomial::isZero() const noexcept { return terms_.empty(); }
std::size_t MultivariateRationalPolynomial::termCount() const noexcept { return terms_.size(); }

std::size_t MultivariateRationalPolynomial::totalDegree() const noexcept {
    std::size_t result = 0;
    for (const PolynomialTerm& term : terms_)
        result = std::max(result, term.monomial.totalDegree());
    return result;
}

std::size_t MultivariateRationalPolynomial::degree(const expression::Symbol& variable) const noexcept {
    std::size_t result = 0;
    for (const PolynomialTerm& term : terms_)
        result = std::max(result, term.monomial.exponentOf(variable));
    return result;
}

std::span<const PolynomialTerm> MultivariateRationalPolynomial::terms() const noexcept { return terms_; }

std::vector<expression::Symbol> MultivariateRationalPolynomial::variables() const {
    std::vector<expression::Symbol> result;
    for (const PolynomialTerm& term : terms_) {
        for (const MonomialFactor& factor : term.monomial.factors()) {
            if (std::find(result.begin(), result.end(), factor.variable) == result.end())
                result.push_back(factor.variable);
        }
    }
    std::sort(result.begin(), result.end(), symbolLess);
    return result;
}

bool MultivariateRationalPolynomial::belongsTo(const PolynomialRing& ring) const noexcept {
    for (const PolynomialTerm& term : terms_)
        if (!ring.contains(term.monomial))
            return false;
    return true;
}

std::optional<PolynomialTerm> MultivariateRationalPolynomial::leadingTerm(
    const PolynomialRing& ring) const {
    if (isZero())
        return std::nullopt;
    if (!belongsTo(ring))
        throw std::invalid_argument("Polynomial does not belong to the polynomial ring");
    const PolynomialTerm* leading = &terms_.front();
    for (const PolynomialTerm& term : terms_)
        if (ring.compare(term.monomial, leading->monomial) > 0)
            leading = &term;
    return *leading;
}

MultivariateRationalPolynomial negatePolynomial(
    const MultivariateRationalPolynomial& value) {
    return negate(value);
}

MultivariateRationalPolynomial addPolynomials(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs) {
    return add(lhs, rhs);
}

MultivariateRationalPolynomial subtractPolynomials(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs) {
    return add(lhs, negate(rhs));
}

std::optional<MultivariateRationalPolynomial> multiplyPolynomials(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs,
    PolynomialConversionOptions options) {
    return multiply(lhs, rhs, options);
}

MultivariateRationalPolynomial multiplyPolynomialByTerm(
    const MultivariateRationalPolynomial& polynomial,
    const PolynomialTerm& multiplier) {
    if (polynomial.isZero() || multiplier.coefficient.isZero())
        return MultivariateRationalPolynomial{};
    std::vector<PolynomialTerm> terms;
    terms.reserve(polynomial.termCount());
    for (const PolynomialTerm& term : polynomial.terms())
        terms.push_back(PolynomialTerm{
            multiplyMonomial(term.monomial, multiplier.monomial),
            term.coefficient * multiplier.coefficient});
    return MultivariateRationalPolynomial{std::move(terms)};
}

MultivariateRationalPolynomial monicPolynomial(
    const MultivariateRationalPolynomial& polynomial,
    const PolynomialRing& ring) {
    const auto leading = polynomial.leadingTerm(ring);
    if (!leading)
        return polynomial;
    std::vector<PolynomialTerm> terms(polynomial.terms().begin(), polynomial.terms().end());
    for (PolynomialTerm& term : terms)
        term.coefficient /= leading->coefficient;
    return MultivariateRationalPolynomial{std::move(terms)};
}

std::optional<MultivariateRationalPolynomial> toMultivariateRationalPolynomial(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    PolynomialConversionOptions options) {
    return convert(expression, builtins, options);
}

std::optional<MultivariateRationalPolynomial> toMultivariateRationalPolynomial(
    const Expr& expression,
    const PolynomialRing& ring,
    const evaluation::BuiltinRegistry& builtins,
    PolynomialConversionOptions options) {
    auto polynomial = convert(expression, builtins, options);
    if (!polynomial || !polynomial->belongsTo(ring))
        return std::nullopt;
    return polynomial;
}

std::optional<RationalPolynomial> toRationalPolynomial(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    PolynomialConversionOptions options) {
    if (!variable.valid())
        return std::nullopt;
    const auto multivariate = convert(expression, builtins, options);
    if (!multivariate)
        return std::nullopt;

    std::vector<Rational> coefficients(multivariate->degree(variable) + 1, zero());
    for (const PolynomialTerm& term : multivariate->terms()) {
        for (const MonomialFactor& factor : term.monomial.factors())
            if (!(factor.variable == variable))
                return std::nullopt;
        coefficients[term.monomial.exponentOf(variable)] += term.coefficient;
    }
    return RationalPolynomial{std::move(coefficients)};
}


bool containsSymbol(
    const Expr& expression,
    const expression::Symbol& symbol) {
    return symbol.valid() && containsSymbolImpl(expression, symbol);
}

std::optional<ExpressionPolynomial> toExpressionPolynomial(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    PolynomialConversionOptions options) {
    if (!variable.valid())
        return std::nullopt;
    return convertExpressionPolynomial(
        expression, variable, builtins, mathematics, angles, options, nullptr);
}

std::optional<ExpressionPolynomialConversionResult>
toExpressionPolynomialWithConditions(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    PolynomialConversionOptions options) {
    if (!variable.valid())
        return std::nullopt;
    mathematics::AssumptionSet domainConditions;
    auto polynomial = convertExpressionPolynomial(
        expression, variable, builtins, mathematics, angles, options, &domainConditions);
    if (!polynomial)
        return std::nullopt;
    return ExpressionPolynomialConversionResult{
        std::move(*polynomial), std::move(domainConditions)};
}

std::optional<mathematics::AssumptionSet> scalarExpressionDomainConditions(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    if (!scalarCoefficientCandidate(expression))
        return std::nullopt;
    return mathematics::expressionDomainConditions(expression, builtins, mathematics);
}

Expr expressionPolynomialToCollectedExpr(
    const ExpressionPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (polynomial.isZero())
        return zeroExpr();

    std::vector<Expr> terms;
    for (std::size_t reverse = 0; reverse <= polynomial.degree(); ++reverse) {
        const std::size_t exponent = polynomial.degree() - reverse;
        const Expr& coefficient = polynomial.coefficient(exponent);
        if (isExactZeroExpr(coefficient))
            continue;
        if (exponent == 0) {
            terms.push_back(coefficient);
            continue;
        }
        Expr powerExpr = variablePower(variable, exponent, builtins);
        if (isExactOneExpr(coefficient)) {
            terms.push_back(std::move(powerExpr));
            continue;
        }
        terms.push_back(Expr::call(
            builtins.symbol(BuiltinId::Multiply), {coefficient, std::move(powerExpr)}));
    }
    if (terms.empty())
        return zeroExpr();
    if (terms.size() == 1)
        return terms.front();
    // collectは「指定変数の降冪順」という表示上の構造自体に意味がある。
    // ここで通常Simplifierへ戻すと可換順序で定数項が先頭へ移るため、
    // 係数は既に正規化済みとしてAddの項順を保持する。
    static_cast<void>(mathematics);
    static_cast<void>(angles);
    return Expr::call(builtins.symbol(BuiltinId::Add), std::move(terms));
}

Expr polynomialToExpandedExpr(
    const RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (polynomial.isZero())
        return integerExpr(0);

    std::vector<PolynomialTerm> terms;
    for (std::size_t reverse = 0; reverse <= polynomial.degree(); ++reverse) {
        const std::size_t exponent = polynomial.degree() - reverse;
        const Rational& coefficient = polynomial.coefficient(exponent);
        if (coefficient.isZero())
            continue;
        Monomial monomial = exponent == 0
            ? Monomial{}
            : Monomial{{MonomialFactor{variable, exponent}}};
        terms.push_back(PolynomialTerm{std::move(monomial), coefficient});
    }
    return termsToExpr(std::move(terms), builtins);
}

Expr polynomialToExpandedExpr(
    const MultivariateRationalPolynomial& polynomial,
    const evaluation::BuiltinRegistry& builtins) {
    return termsToExpr(
        std::vector<PolynomialTerm>{polynomial.terms().begin(), polynomial.terms().end()},
        builtins);
}

Expr polynomialToCollectedExpr(
    const MultivariateRationalPolynomial& polynomial,
    std::span<const expression::Symbol> variables,
    const evaluation::BuiltinRegistry& builtins) {
    if (polynomial.isZero())
        return integerExpr(0);
    if (variables.empty())
        return polynomialToExpandedExpr(polynomial, builtins);
    return collectRecursive(polynomial.terms(), variables, 0, builtins);
}

Rational evaluatePolynomial(const RationalPolynomial& polynomial, const Rational& value) {
    Rational result{BigInt{0}};
    for (std::size_t i = polynomial.degree() + 1; i-- > 0;)
        result = result * value + polynomial.coefficient(i);
    return result;
}

std::optional<RationalPolynomial> divideByLinearFactor(
    const RationalPolynomial& polynomial,
    const Rational& root) {
    if (polynomial.degree() == 0)
        return std::nullopt;

    std::vector<Rational> quotient(polynomial.degree(), zero());
    quotient.back() = polynomial.coefficient(polynomial.degree());
    for (std::size_t exponent = polynomial.degree() - 1; exponent > 0; --exponent)
        quotient[exponent - 1] = polynomial.coefficient(exponent) + root * quotient[exponent];
    const Rational remainder = polynomial.coefficient(0) + root * quotient[0];
    if (!remainder.isZero())
        return std::nullopt;
    return RationalPolynomial{std::move(quotient)};
}

RationalRootSearchResult findRationalRoot(
    const RationalPolynomial& polynomial,
    RationalRootSearchOptions options) {
    if (polynomial.degree() == 0)
        return RationalRootSearchResult{std::nullopt, true};
    if (polynomial.coefficient(0).isZero())
        return RationalRootSearchResult{Rational{BigInt{0}}, true};

    const std::vector<BigInt> integerCoefficients = primitiveIntegerCoefficients(polynomial);
    const DivisorList numerators = positiveDivisors(
        integerCoefficients.front(), options.maximumTrialDivisor);
    const DivisorList denominators = positiveDivisors(
        integerCoefficients.back(), options.maximumTrialDivisor);
    if (numerators.values.empty() || denominators.values.empty())
        return RationalRootSearchResult{std::nullopt, false};

    std::size_t tested = 0;
    for (const std::uint64_t p : numerators.values) {
        for (const std::uint64_t q : denominators.values) {
            if (q == 0)
                continue;
            if (++tested > options.maximumCandidates)
                return RationalRootSearchResult{std::nullopt, false};
            Rational positive{BigInt::parse(std::to_string(p)), BigInt::parse(std::to_string(q))};
            if (evaluatePolynomial(polynomial, positive).isZero())
                return RationalRootSearchResult{std::move(positive), true};
            Rational negative = -positive;
            if (evaluatePolynomial(polynomial, negative).isZero())
                return RationalRootSearchResult{std::move(negative), true};
        }
    }

    return RationalRootSearchResult{
        std::nullopt,
        numerators.complete && denominators.complete};
}

namespace {

[[nodiscard]] RationalPolynomial monicRationalPolynomial(
    const RationalPolynomial& polynomial) {
    if (polynomial.isZero())
        return polynomial;
    const Rational leading = polynomial.coefficient(polynomial.degree());
    std::vector<Rational> coefficients(polynomial.coefficients().begin(),
        polynomial.coefficients().end());
    for (Rational& coefficient : coefficients)
        coefficient /= leading;
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] std::optional<RationalPolynomial> divideRationalPolynomialsExact(
    const RationalPolynomial& dividend,
    const RationalPolynomial& divisor) {
    if (divisor.isZero() || dividend.degree() < divisor.degree())
        return std::nullopt;

    std::vector<Rational> remainder(dividend.coefficients().begin(),
        dividend.coefficients().end());
    std::vector<Rational> quotient(
        dividend.degree() - divisor.degree() + 1, zero());
    const Rational leading = divisor.coefficient(divisor.degree());
    for (std::size_t degree = dividend.degree() + 1; degree-- > divisor.degree();) {
        const Rational amount = remainder[degree] / leading;
        const std::size_t shift = degree - divisor.degree();
        quotient[shift] += amount;
        for (std::size_t i = 0; i <= divisor.degree(); ++i)
            remainder[i + shift] -= amount * divisor.coefficient(i);
    }
    for (std::size_t i = 0; i < divisor.degree(); ++i)
        if (!remainder[i].isZero())
            return std::nullopt;
    return RationalPolynomial{std::move(quotient)};
}

[[nodiscard]] RationalPolynomial multiplyRationalPolynomialValues(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    if (lhs.isZero() || rhs.isZero())
        return RationalPolynomial{};
    std::vector<Rational> result(lhs.degree() + rhs.degree() + 1, zero());
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        for (std::size_t j = 0; j <= rhs.degree(); ++j)
            result[i + j] += lhs.coefficient(i) * rhs.coefficient(j);
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial rationalPolynomialRemainder(
    const RationalPolynomial& dividend,
    const RationalPolynomial& divisor) {
    if (divisor.isZero())
        return dividend;
    if (dividend.degree() < divisor.degree())
        return dividend;

    std::vector<Rational> remainder(
        dividend.coefficients().begin(), dividend.coefficients().end());
    const Rational leading = divisor.coefficient(divisor.degree());
    for (std::size_t degree = dividend.degree() + 1; degree-- > divisor.degree();) {
        const Rational amount = remainder[degree] / leading;
        const std::size_t shift = degree - divisor.degree();
        for (std::size_t i = 0; i <= divisor.degree(); ++i)
            remainder[i + shift] -= amount * divisor.coefficient(i);
    }
    remainder.resize(divisor.degree());
    return RationalPolynomial{std::move(remainder)};
}

[[nodiscard]] RationalPolynomial rationalPolynomialDerivative(
    const RationalPolynomial& polynomial) {
    if (polynomial.degree() == 0)
        return RationalPolynomial{};
    std::vector<Rational> result(polynomial.degree(), zero());
    for (std::size_t exponent = 1; exponent <= polynomial.degree(); ++exponent)
        result[exponent - 1] = polynomial.coefficient(exponent)
            * Rational{BigInt::fromUnsigned(exponent)};
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial rationalPolynomialGcdMonic(
    RationalPolynomial lhs,
    RationalPolynomial rhs) {
    while (!rhs.isZero()) {
        RationalPolynomial remainder = rationalPolynomialRemainder(lhs, rhs);
        lhs = std::move(rhs);
        rhs = std::move(remainder);
    }
    return monicRationalPolynomial(lhs);
}

[[nodiscard]] bool isOnePolynomial(const RationalPolynomial& polynomial) {
    return polynomial.degree() == 0 && polynomial.coefficient(0) == one();
}

struct SquareFreeRationalFactor final {
    RationalPolynomial polynomial;
    std::size_t multiplicity = 1;
};

[[nodiscard]] std::optional<std::vector<SquareFreeRationalFactor>>
squareFreeRationalFactorization(const RationalPolynomial& input) {
    const RationalPolynomial polynomial = monicRationalPolynomial(input);
    if (polynomial.degree() <= 1)
        return std::vector<SquareFreeRationalFactor>{{polynomial, 1}};

    const RationalPolynomial derivative = rationalPolynomialDerivative(polynomial);
    RationalPolynomial repeated = rationalPolynomialGcdMonic(polynomial, derivative);
    auto squareFree = divideRationalPolynomialsExact(polynomial, repeated);
    if (!squareFree)
        return std::nullopt;

    std::vector<SquareFreeRationalFactor> result;
    RationalPolynomial current = monicRationalPolynomial(*squareFree);
    std::size_t multiplicity = 1;
    while (!isOnePolynomial(current)) {
        if (multiplicity > polynomial.degree())
            return std::nullopt;
        RationalPolynomial shared = rationalPolynomialGcdMonic(current, repeated);
        auto component = divideRationalPolynomialsExact(current, shared);
        auto nextRepeated = divideRationalPolynomialsExact(repeated, shared);
        if (!component || !nextRepeated)
            return std::nullopt;
        RationalPolynomial normalizedComponent = monicRationalPolynomial(*component);
        if (!isOnePolynomial(normalizedComponent))
            result.push_back(SquareFreeRationalFactor{
                std::move(normalizedComponent), multiplicity});
        current = monicRationalPolynomial(shared);
        repeated = monicRationalPolynomial(*nextRepeated);
        ++multiplicity;
    }
    return result;
}

using ModPolynomial = std::vector<std::uint32_t>;
using IntegerPolynomial = std::vector<BigInt>;

void normalizeModPolynomial(ModPolynomial& polynomial) {
    while (!polynomial.empty() && polynomial.back() == 0)
        polynomial.pop_back();
}

void normalizeIntegerPolynomial(IntegerPolynomial& polynomial) {
    while (polynomial.size() > 1 && polynomial.back().isZero())
        polynomial.pop_back();
    if (polynomial.empty())
        polynomial.emplace_back(0);
}

[[nodiscard]] std::uint32_t modularMultiply(
    std::uint32_t lhs,
    std::uint32_t rhs,
    std::uint32_t prime) {
    return static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(lhs) * rhs) % prime);
}

[[nodiscard]] std::uint32_t modularPower(
    std::uint32_t base,
    std::uint64_t exponent,
    std::uint32_t prime) {
    std::uint32_t result = 1;
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = modularMultiply(result, base, prime);
        exponent >>= 1U;
        if (exponent != 0)
            base = modularMultiply(base, base, prime);
    }
    return result;
}

[[nodiscard]] ModPolynomial subtractModPolynomials(
    const ModPolynomial& lhs,
    const ModPolynomial& rhs,
    std::uint32_t prime) {
    ModPolynomial result(std::max(lhs.size(), rhs.size()), 0);
    for (std::size_t i = 0; i < result.size(); ++i) {
        const std::uint32_t left = i < lhs.size() ? lhs[i] : 0;
        const std::uint32_t right = i < rhs.size() ? rhs[i] : 0;
        result[i] = left >= right ? left - right : left + prime - right;
    }
    normalizeModPolynomial(result);
    return result;
}

[[nodiscard]] ModPolynomial scaleModPolynomial(
    ModPolynomial polynomial,
    std::uint32_t scalar,
    std::uint32_t prime) {
    for (std::uint32_t& coefficient : polynomial)
        coefficient = modularMultiply(coefficient, scalar, prime);
    normalizeModPolynomial(polynomial);
    return polynomial;
}

[[nodiscard]] ModPolynomial multiplyModPolynomials(
    const ModPolynomial& lhs,
    const ModPolynomial& rhs,
    std::uint32_t prime) {
    if (lhs.empty() || rhs.empty())
        return {};
    ModPolynomial result(lhs.size() + rhs.size() - 1, 0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result[i + j] = static_cast<std::uint32_t>(
                (static_cast<std::uint64_t>(result[i + j])
                    + modularMultiply(lhs[i], rhs[j], prime)) % prime);
        }
    }
    normalizeModPolynomial(result);
    return result;
}

[[nodiscard]] std::pair<ModPolynomial, ModPolynomial> divideModPolynomials(
    ModPolynomial dividend,
    const ModPolynomial& divisor,
    std::uint32_t prime) {
    normalizeModPolynomial(dividend);
    if (divisor.empty() || dividend.size() < divisor.size())
        return {{}, std::move(dividend)};
    ModPolynomial quotient(dividend.size() - divisor.size() + 1, 0);
    const std::uint32_t inverseLeading = modularPower(
        divisor.back(), prime - 2, prime);
    while (!dividend.empty() && dividend.size() >= divisor.size()) {
        const std::size_t shift = dividend.size() - divisor.size();
        const std::uint32_t amount = modularMultiply(
            dividend.back(), inverseLeading, prime);
        quotient[shift] = amount;
        for (std::size_t i = 0; i < divisor.size(); ++i) {
            const std::uint32_t product = modularMultiply(amount, divisor[i], prime);
            const std::size_t index = i + shift;
            dividend[index] = dividend[index] >= product
                ? dividend[index] - product
                : dividend[index] + prime - product;
        }
        normalizeModPolynomial(dividend);
    }
    normalizeModPolynomial(quotient);
    return {std::move(quotient), std::move(dividend)};
}

[[nodiscard]] ModPolynomial remainderModPolynomial(
    ModPolynomial dividend,
    const ModPolynomial& divisor,
    std::uint32_t prime) {
    return divideModPolynomials(std::move(dividend), divisor, prime).second;
}

[[nodiscard]] ModPolynomial gcdModPolynomial(
    ModPolynomial lhs,
    ModPolynomial rhs,
    std::uint32_t prime) {
    while (!rhs.empty()) {
        ModPolynomial remainder = remainderModPolynomial(
            std::move(lhs), rhs, prime);
        lhs = std::move(rhs);
        rhs = std::move(remainder);
    }
    if (lhs.empty())
        return lhs;
    const std::uint32_t inverseLeading = modularPower(
        lhs.back(), prime - 2, prime);
    return scaleModPolynomial(
        std::move(lhs), inverseLeading, prime);
}

[[nodiscard]] ModPolynomial derivativeModPolynomial(
    const ModPolynomial& polynomial,
    std::uint32_t prime) {
    if (polynomial.size() <= 1)
        return {};
    ModPolynomial result(polynomial.size() - 1, 0);
    for (std::size_t i = 1; i < polynomial.size(); ++i)
        result[i - 1] = modularMultiply(
            polynomial[i], static_cast<std::uint32_t>(i % prime), prime);
    normalizeModPolynomial(result);
    return result;
}

[[nodiscard]] ModPolynomial multiplyReduceModPolynomial(
    const ModPolynomial& lhs,
    const ModPolynomial& rhs,
    const ModPolynomial& modulus,
    std::uint32_t prime) {
    return remainderModPolynomial(
        multiplyModPolynomials(lhs, rhs, prime), modulus, prime);
}

[[nodiscard]] ModPolynomial powerModPolynomial(
    ModPolynomial base,
    std::uint64_t exponent,
    const ModPolynomial& modulus,
    std::uint32_t prime) {
    ModPolynomial result{1};
    base = remainderModPolynomial(std::move(base), modulus, prime);
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = multiplyReduceModPolynomial(result, base, modulus, prime);
        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyReduceModPolynomial(base, base, modulus, prime);
    }
    return result;
}

[[nodiscard]] ModPolynomial powerModPolynomial(
    ModPolynomial base,
    BigInt exponent,
    const ModPolynomial& modulus,
    std::uint32_t prime) {
    ModPolynomial result{1};
    base = remainderModPolynomial(std::move(base), modulus, prime);
    while (!exponent.isZero()) {
        if (exponent.testBit(0))
            result = multiplyReduceModPolynomial(result, base, modulus, prime);
        exponent >>= 1;
        if (!exponent.isZero())
            base = multiplyReduceModPolynomial(base, base, modulus, prime);
    }
    return result;
}

[[nodiscard]] std::optional<ModPolynomial> inverseModPolynomial(
    ModPolynomial value,
    const ModPolynomial& modulus,
    std::uint32_t prime) {
    value = remainderModPolynomial(std::move(value), modulus, prime);
    ModPolynomial oldRemainder = modulus;
    ModPolynomial remainder = std::move(value);
    ModPolynomial oldCoefficient;
    ModPolynomial coefficient{1};
    while (!remainder.empty()) {
        auto [quotient, nextRemainder] = divideModPolynomials(
            oldRemainder, remainder, prime);
        ModPolynomial nextCoefficient = subtractModPolynomials(
            oldCoefficient,
            multiplyModPolynomials(quotient, coefficient, prime),
            prime);
        oldRemainder = std::move(remainder);
        remainder = std::move(nextRemainder);
        oldCoefficient = std::move(coefficient);
        coefficient = std::move(nextCoefficient);
    }
    if (oldRemainder.size() != 1)
        return std::nullopt;
    const std::uint32_t inverse = modularPower(
        oldRemainder.front(), prime - 2, prime);
    return remainderModPolynomial(
        scaleModPolynomial(std::move(oldCoefficient), inverse, prime),
        modulus, prime);
}

[[nodiscard]] ModPolynomial integerPolynomialModulo(
    const IntegerPolynomial& polynomial,
    std::uint32_t prime) {
    ModPolynomial result;
    result.reserve(polynomial.size());
    for (const BigInt& coefficient : polynomial)
        result.push_back(coefficient.modulo(prime));
    normalizeModPolynomial(result);
    return result;
}

[[nodiscard]] bool hasModularSquareFreeCertificate(
    const RationalPolynomial& polynomial) {
    if (polynomial.degree() <= 1)
        return true;
    constexpr std::array<std::uint32_t, 8> primes{
        3, 5, 7, 11, 13, 17, 19, 23};
    for (const std::uint32_t prime : primes) {
        ModPolynomial reduced(polynomial.degree() + 1, 0);
        bool usable = true;
        for (std::size_t i = 0; i <= polynomial.degree(); ++i) {
            const Rational& coefficient = polynomial.coefficient(i);
            const std::uint32_t denominator = coefficient.denominator().modulo(prime);
            if (denominator == 0) {
                usable = false;
                break;
            }
            reduced[i] = modularMultiply(
                coefficient.numerator().modulo(prime),
                modularPower(denominator, prime - 2, prime), prime);
        }
        if (!usable)
            continue;
        normalizeModPolynomial(reduced);
        if (reduced.size() != polynomial.degree() + 1)
            continue;
        const ModPolynomial derivative = derivativeModPolynomial(reduced, prime);
        if (!derivative.empty()
            && gcdModPolynomial(reduced, derivative, prime).size() == 1)
            return true;
    }
    return false;
}

[[nodiscard]] std::vector<ModPolynomial> nullspaceModPrime(
    std::vector<std::vector<std::uint32_t>> matrix,
    std::uint32_t prime) {
    if (matrix.empty())
        return {};
    const std::size_t rowCount = matrix.size();
    const std::size_t columnCount = matrix.front().size();
    std::vector<std::size_t> pivotColumns;
    std::size_t pivotRow = 0;
    for (std::size_t column = 0;
         column < columnCount && pivotRow < rowCount; ++column) {
        std::size_t selected = pivotRow;
        while (selected < rowCount && matrix[selected][column] == 0)
            ++selected;
        if (selected == rowCount)
            continue;
        std::swap(matrix[pivotRow], matrix[selected]);
        const std::uint32_t inverse = modularPower(
            matrix[pivotRow][column], prime - 2, prime);
        for (std::size_t j = 0; j < columnCount; ++j)
            matrix[pivotRow][j] = modularMultiply(
                matrix[pivotRow][j], inverse, prime);
        for (std::size_t row = 0; row < rowCount; ++row) {
            if (row == pivotRow || matrix[row][column] == 0)
                continue;
            const std::uint32_t amount = matrix[row][column];
            for (std::size_t j = 0; j < columnCount; ++j) {
                const std::uint32_t product = modularMultiply(
                    amount, matrix[pivotRow][j], prime);
                matrix[row][j] = matrix[row][j] >= product
                    ? matrix[row][j] - product
                    : matrix[row][j] + prime - product;
            }
        }
        pivotColumns.push_back(column);
        ++pivotRow;
    }

    std::vector<bool> isPivot(columnCount, false);
    for (const std::size_t column : pivotColumns)
        isPivot[column] = true;
    std::vector<ModPolynomial> basis;
    for (std::size_t freeColumn = 0; freeColumn < columnCount; ++freeColumn) {
        if (isPivot[freeColumn])
            continue;
        ModPolynomial vector(columnCount, 0);
        vector[freeColumn] = 1;
        for (std::size_t row = 0; row < pivotColumns.size(); ++row) {
            const std::uint32_t coefficient = matrix[row][freeColumn];
            vector[pivotColumns[row]] = coefficient == 0 ? 0 : prime - coefficient;
        }
        normalizeModPolynomial(vector);
        basis.push_back(std::move(vector));
    }
    return basis;
}

[[nodiscard]] std::optional<std::vector<ModPolynomial>> berlekampFactorization(
    const ModPolynomial& polynomial,
    std::uint32_t prime) {
    if (polynomial.size() <= 2)
        return std::vector<ModPolynomial>{polynomial};
    const std::size_t degree = polynomial.size() - 1;
    std::vector<std::vector<std::uint32_t>> matrix(
        degree, std::vector<std::uint32_t>(degree, 0));
    const ModPolynomial x{0, 1};
    for (std::size_t column = 0; column < degree; ++column) {
        const ModPolynomial power = powerModPolynomial(
            x, static_cast<std::uint64_t>(prime) * column,
            polynomial, prime);
        for (std::size_t row = 0; row < power.size(); ++row)
            matrix[row][column] = power[row];
        matrix[column][column] = matrix[column][column] == 0
            ? prime - 1
            : matrix[column][column] - 1;
    }
    const std::vector<ModPolynomial> basis = nullspaceModPrime(
        std::move(matrix), prime);
    if (basis.empty())
        return std::nullopt;
    if (basis.size() == 1)
        return std::vector<ModPolynomial>{polynomial};

    std::vector<ModPolynomial> factors{polynomial};
    for (const ModPolynomial& splitter : basis) {
        if (splitter.size() <= 1)
            continue;
        for (std::uint32_t value = 0;
             value < prime && factors.size() < basis.size(); ++value) {
            std::vector<ModPolynomial> next;
            for (const ModPolynomial& factor : factors) {
                if (factor.size() <= 2) {
                    next.push_back(factor);
                    continue;
                }
                ModPolynomial shifted = remainderModPolynomial(
                    splitter, factor, prime);
                if (shifted.empty())
                    shifted.resize(1, 0);
                shifted[0] = shifted[0] >= value
                    ? shifted[0] - value
                    : shifted[0] + prime - value;
                normalizeModPolynomial(shifted);
                ModPolynomial common = gcdModPolynomial(factor, shifted, prime);
                if (common.size() <= 1 || common.size() == factor.size()) {
                    next.push_back(factor);
                    continue;
                }
                auto [quotient, remainder] = divideModPolynomials(
                    factor, common, prime);
                if (!remainder.empty() || quotient.empty())
                    return std::nullopt;
                next.push_back(std::move(common));
                next.push_back(std::move(quotient));
            }
            factors = std::move(next);
        }
        if (factors.size() == basis.size())
            break;
    }
    if (factors.size() != basis.size())
        return std::nullopt;
    std::sort(factors.begin(), factors.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.size() != rhs.size())
            return lhs.size() < rhs.size();
        return lhs < rhs;
    });
    return factors;
}

struct DistinctDegreeFactor final {
    ModPolynomial polynomial;
    std::size_t irreducibleDegree = 0;
};

[[nodiscard]] std::optional<std::vector<DistinctDegreeFactor>>
distinctDegreeFactorization(
    const ModPolynomial& polynomial,
    std::uint32_t prime) {
    if (polynomial.size() <= 1)
        return std::vector<DistinctDegreeFactor>{};

    std::vector<DistinctDegreeFactor> result;
    ModPolynomial remaining = polynomial;
    ModPolynomial frobenius{0, 1};
    const ModPolynomial x{0, 1};
    for (std::size_t degree = 1;
         remaining.size() > 1 && 2 * degree <= remaining.size() - 1;
         ++degree) {
        frobenius = powerModPolynomial(
            std::move(frobenius), static_cast<std::uint64_t>(prime),
            remaining, prime);
        ModPolynomial common = gcdModPolynomial(
            remaining,
            subtractModPolynomials(frobenius, x, prime),
            prime);
        if (common.size() <= 1)
            continue;

        auto [quotient, remainder] = divideModPolynomials(
            remaining, common, prime);
        if (!remainder.empty() || quotient.empty())
            return std::nullopt;
        result.push_back(DistinctDegreeFactor{
            std::move(common), degree});
        remaining = std::move(quotient);
        if (remaining.size() > 1)
            frobenius = remainderModPolynomial(
                std::move(frobenius), remaining, prime);
    }
    if (remaining.size() > 1) {
        const std::size_t remainingDegree = remaining.size() - 1;
        result.push_back(DistinctDegreeFactor{
            std::move(remaining), remainingDegree});
    }
    return result;
}

[[nodiscard]] std::uint64_t nextFactorRandom(std::uint64_t& state) noexcept {
    // 再現可能なLas Vegas試行にする。乱数は候補生成にしか使わず，
    // 採用は有限体上のexact gcd/divisionで検証する。
    state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t value = state;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

[[nodiscard]] std::uint64_t modPolynomialSeed(
    const ModPolynomial& polynomial,
    std::uint32_t prime,
    std::size_t equalDegree) noexcept {
    std::uint64_t seed = 0xcbf29ce484222325ULL
        ^ static_cast<std::uint64_t>(prime)
        ^ (static_cast<std::uint64_t>(equalDegree) << 32U);
    for (const std::uint32_t coefficient : polynomial) {
        seed ^= coefficient;
        seed *= 0x100000001b3ULL;
    }
    return seed;
}

[[nodiscard]] std::optional<std::vector<ModPolynomial>>
equalDegreeFactorization(
    const ModPolynomial& polynomial,
    std::size_t irreducibleDegree,
    std::uint32_t prime,
    std::size_t maximumAttempts) {
    if (polynomial.size() <= 1 || irreducibleDegree == 0
        || (polynomial.size() - 1) % irreducibleDegree != 0)
        return std::nullopt;
    if (polynomial.size() - 1 == irreducibleDegree)
        return std::vector<ModPolynomial>{polynomial};

    const BigInt exponent = (
        numeric::pow(BigInt::fromUnsigned(prime),
            static_cast<std::uint64_t>(irreducibleDegree)) - BigInt{1})
        / BigInt{2};
    std::uint64_t randomState = modPolynomialSeed(
        polynomial, prime, irreducibleDegree);
    std::vector<ModPolynomial> factors;

    std::function<bool(const ModPolynomial&)> split =
        [&](const ModPolynomial& current) {
            const std::size_t degree = current.size() - 1;
            if (degree == irreducibleDegree) {
                factors.push_back(current);
                return true;
            }
            if (degree < irreducibleDegree
                || degree % irreducibleDegree != 0)
                return false;

            for (std::size_t attempt = 0; attempt < maximumAttempts; ++attempt) {
                ModPolynomial candidate(degree, 0);
                bool nonConstant = false;
                for (std::size_t i = 0; i < candidate.size(); ++i) {
                    candidate[i] = static_cast<std::uint32_t>(
                        nextFactorRandom(randomState) % prime);
                    nonConstant = nonConstant || (i != 0 && candidate[i] != 0);
                }
                normalizeModPolynomial(candidate);
                if (!nonConstant || candidate.empty())
                    continue;

                ModPolynomial common = gcdModPolynomial(
                    current, candidate, prime);
                if (common.size() <= 1 || common.size() == current.size()) {
                    ModPolynomial powered = powerModPolynomial(
                        candidate, exponent, current, prime);
                    common = gcdModPolynomial(
                        current,
                        subtractModPolynomials(powered, ModPolynomial{1}, prime),
                        prime);
                }
                if (common.size() <= 1 || common.size() == current.size())
                    continue;

                auto [quotient, remainder] = divideModPolynomials(
                    current, common, prime);
                if (!remainder.empty() || quotient.empty())
                    return false;
                return split(common) && split(quotient);
            }
            return false;
        };

    if (!split(polynomial))
        return std::nullopt;
    std::sort(factors.begin(), factors.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.size() != rhs.size())
            return lhs.size() < rhs.size();
        return lhs < rhs;
    });
    return factors;
}

[[nodiscard]] std::optional<std::vector<ModPolynomial>>
cantorZassenhausFactorization(
    const ModPolynomial& polynomial,
    std::uint32_t prime,
    std::size_t maximumAttempts) {
    if (polynomial.size() <= 2)
        return std::vector<ModPolynomial>{polynomial};
    if (prime == 2)
        return std::nullopt;

    const auto distinct = distinctDegreeFactorization(polynomial, prime);
    if (!distinct)
        return std::nullopt;
    std::vector<ModPolynomial> result;
    for (const DistinctDegreeFactor& component : *distinct) {
        auto equal = equalDegreeFactorization(
            component.polynomial, component.irreducibleDegree,
            prime, maximumAttempts);
        if (!equal)
            return std::nullopt;
        result.insert(
            result.end(),
            std::make_move_iterator(equal->begin()),
            std::make_move_iterator(equal->end()));
    }

    ModPolynomial product{1};
    for (const ModPolynomial& factor : result)
        product = multiplyModPolynomials(product, factor, prime);
    if (product != polynomial)
        return std::nullopt;
    std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.size() != rhs.size())
            return lhs.size() < rhs.size();
        return lhs < rhs;
    });
    return result;
}

[[nodiscard]] std::optional<std::vector<ModPolynomial>>
factorFiniteFieldPolynomial(
    const ModPolynomial& polynomial,
    std::uint32_t prime,
    const RationalPolynomialFactorOptions& options) {
    const std::size_t degree = polynomial.size() - 1;
    const bool berlekampFits = degree == 0
        || degree <= options.maximumBerlekampMatrixEntries / degree;
    if (berlekampFits)
        return berlekampFactorization(polynomial, prime);
    return cantorZassenhausFactorization(
        polynomial, prime, options.maximumCantorZassenhausAttempts);
}

[[nodiscard]] IntegerPolynomial multiplyIntegerPolynomials(
    const IntegerPolynomial& lhs,
    const IntegerPolynomial& rhs) {
    if ((lhs.size() == 1 && lhs.front().isZero())
        || (rhs.size() == 1 && rhs.front().isZero()))
        return IntegerPolynomial{BigInt{0}};
    IntegerPolynomial result(lhs.size() + rhs.size() - 1, BigInt{0});
    for (std::size_t i = 0; i < lhs.size(); ++i)
        for (std::size_t j = 0; j < rhs.size(); ++j)
            result[i + j] += lhs[i] * rhs[j];
    normalizeIntegerPolynomial(result);
    return result;
}

[[nodiscard]] BigInt positiveResidue(
    BigInt value,
    const BigInt& modulus) {
    value %= modulus;
    if (value.isNegative())
        value += modulus;
    return value;
}

[[nodiscard]] IntegerPolynomial reduceIntegerPolynomialModulo(
    IntegerPolynomial polynomial,
    const BigInt& modulus) {
    for (BigInt& coefficient : polynomial)
        coefficient = positiveResidue(std::move(coefficient), modulus);
    normalizeIntegerPolynomial(polynomial);
    return polynomial;
}

[[nodiscard]] bool isZeroIntegerPolynomial(
    const IntegerPolynomial& polynomial) noexcept {
    return polynomial.empty()
        || (polynomial.size() == 1 && polynomial.front().isZero());
}

[[nodiscard]] IntegerPolynomial multiplyIntegerPolynomialsModulo(
    const IntegerPolynomial& lhs,
    const IntegerPolynomial& rhs,
    const BigInt& modulus) {
    IntegerPolynomial result(lhs.size() + rhs.size() - 1, BigInt{0});
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result[i + j] += lhs[i] * rhs[j];
            result[i + j] = positiveResidue(
                std::move(result[i + j]), modulus);
        }
    }
    normalizeIntegerPolynomial(result);
    return result;
}

[[nodiscard]] IntegerPolynomial subtractIntegerPolynomialsModulo(
    const IntegerPolynomial& lhs,
    const IntegerPolynomial& rhs,
    const BigInt& modulus) {
    IntegerPolynomial result(std::max(lhs.size(), rhs.size()), BigInt{0});
    for (std::size_t i = 0; i < result.size(); ++i) {
        const BigInt left = i < lhs.size() ? lhs[i] : BigInt{0};
        const BigInt right = i < rhs.size() ? rhs[i] : BigInt{0};
        result[i] = positiveResidue(left - right, modulus);
    }
    normalizeIntegerPolynomial(result);
    return result;
}

[[nodiscard]] std::pair<IntegerPolynomial, IntegerPolynomial>
divideMonicIntegerPolynomialsModulo(
    IntegerPolynomial dividend,
    const IntegerPolynomial& divisor,
    const BigInt& modulus) {
    dividend = reduceIntegerPolynomialModulo(std::move(dividend), modulus);
    if (divisor.empty() || divisor.back() != BigInt{1}
        || dividend.size() < divisor.size())
        return {{BigInt{0}}, std::move(dividend)};

    IntegerPolynomial quotient(
        dividend.size() - divisor.size() + 1, BigInt{0});
    while (dividend.size() >= divisor.size()
        && !isZeroIntegerPolynomial(dividend)) {
        const std::size_t shift = dividend.size() - divisor.size();
        const BigInt amount = dividend.back();
        quotient[shift] = amount;
        for (std::size_t i = 0; i < divisor.size(); ++i) {
            const std::size_t index = i + shift;
            dividend[index] = positiveResidue(
                dividend[index] - amount * divisor[i], modulus);
        }
        normalizeIntegerPolynomial(dividend);
    }
    normalizeIntegerPolynomial(quotient);
    return {std::move(quotient), std::move(dividend)};
}

[[nodiscard]] IntegerPolynomial remainderMonicIntegerPolynomialModulo(
    IntegerPolynomial dividend,
    const IntegerPolynomial& divisor,
    const BigInt& modulus) {
    return divideMonicIntegerPolynomialsModulo(
        std::move(dividend), divisor, modulus).second;
}

[[nodiscard]] bool integerPolynomialsCongruent(
    const IntegerPolynomial& lhs,
    const IntegerPolynomial& rhs,
    const BigInt& modulus) {
    const std::size_t size = std::max(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < size; ++i) {
        const BigInt left = i < lhs.size() ? lhs[i] : BigInt{0};
        const BigInt right = i < rhs.size() ? rhs[i] : BigInt{0};
        if (!((left - right) % modulus).isZero())
            return false;
    }
    return true;
}

[[nodiscard]] std::optional<IntegerPolynomial> divideIntegerPolynomialsExact(
    const IntegerPolynomial& dividend,
    const IntegerPolynomial& divisor) {
    if (divisor.empty() || divisor.back() != BigInt{1}
        || dividend.size() < divisor.size())
        return std::nullopt;
    IntegerPolynomial remainder = dividend;
    IntegerPolynomial quotient(dividend.size() - divisor.size() + 1, BigInt{0});
    for (std::size_t degree = dividend.size(); degree-- >= divisor.size() - 1;) {
        const std::size_t shift = degree - (divisor.size() - 1);
        const BigInt amount = remainder[degree];
        quotient[shift] = amount;
        for (std::size_t i = 0; i < divisor.size(); ++i)
            remainder[i + shift] -= amount * divisor[i];
        if (degree == divisor.size() - 1)
            break;
    }
    for (std::size_t i = 0; i + 1 < divisor.size(); ++i)
        if (!remainder[i].isZero())
            return std::nullopt;
    normalizeIntegerPolynomial(quotient);
    return quotient;
}

[[nodiscard]] std::optional<IntegerPolynomial> monicIntegerTransform(
    const RationalPolynomial& polynomial,
    BigInt& scale,
    std::size_t maximumBits) {
    const IntegerPolynomial primitive = primitiveIntegerCoefficients(polynomial);
    if (primitive.empty() || primitive.back().isZero())
        return std::nullopt;
    scale = primitive.back();
    if (scale.isNegative())
        scale = -scale;
    const std::size_t degree = polynomial.degree();
    IntegerPolynomial transformed(degree + 1, BigInt{0});
    for (std::size_t exponent = 0; exponent < degree; ++exponent) {
        transformed[exponent] = primitive[exponent]
            * numeric::pow(scale, static_cast<std::uint64_t>(degree - 1 - exponent));
        if (transformed[exponent].abs().bitLength() > maximumBits)
            return std::nullopt;
    }
    transformed[degree] = BigInt{1};
    return transformed;
}

[[nodiscard]] RationalPolynomial inverseMonicIntegerTransform(
    const IntegerPolynomial& polynomial,
    const BigInt& scale) {
    const std::size_t degree = polynomial.size() - 1;
    std::vector<Rational> coefficients(degree + 1, zero());
    for (std::size_t exponent = 0; exponent <= degree; ++exponent) {
        coefficients[exponent] = exponent == degree
            ? one()
            : Rational{polynomial[exponent], numeric::pow(
                scale, static_cast<std::uint64_t>(degree - exponent))};
    }
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] std::size_t mignotteCoefficientBoundBits(
    const IntegerPolynomial& polynomial) {
    std::size_t maximum = 0;
    for (const BigInt& coefficient : polynomial)
        maximum = std::max(maximum, coefficient.abs().bitLength());
    std::size_t rootSizeBits = 0;
    for (std::size_t n = polynomial.size(); n > 1; n = (n + 1) / 2)
        ++rootSizeBits;
    return maximum + polynomial.size() - 1 + (rootSizeBits + 1) / 2 + 2;
}

struct HenselLiftResult final {
    std::vector<IntegerPolynomial> factors;
    BigInt modulus;
};

[[nodiscard]] std::optional<HenselLiftResult> henselLiftFactorsLinear(
    const IntegerPolynomial& polynomial,
    const std::vector<ModPolynomial>& modularFactors,
    std::uint32_t prime,
    std::size_t targetBits,
    std::size_t maximumWorkingBits) {
    std::vector<IntegerPolynomial> lifted;
    lifted.reserve(modularFactors.size());
    for (const ModPolynomial& factor : modularFactors) {
        IntegerPolynomial value;
        value.reserve(factor.size());
        for (const std::uint32_t coefficient : factor)
            value.push_back(BigInt::fromUnsigned(coefficient));
        lifted.push_back(std::move(value));
    }

    std::vector<ModPolynomial> inverseCofactors;
    inverseCofactors.reserve(modularFactors.size());
    for (std::size_t i = 0; i < modularFactors.size(); ++i) {
        ModPolynomial cofactor{1};
        for (std::size_t j = 0; j < modularFactors.size(); ++j)
            if (j != i)
                cofactor = multiplyModPolynomials(cofactor, modularFactors[j], prime);
        auto inverse = inverseModPolynomial(
            remainderModPolynomial(cofactor, modularFactors[i], prime),
            modularFactors[i], prime);
        if (!inverse)
            return std::nullopt;
        inverseCofactors.push_back(std::move(*inverse));
    }

    BigInt modulus = BigInt::fromUnsigned(prime);
    while (modulus.bitLength() <= targetBits) {
        IntegerPolynomial product{BigInt{1}};
        for (const IntegerPolynomial& factor : lifted)
            product = multiplyIntegerPolynomials(product, factor);
        const std::size_t size = std::max(polynomial.size(), product.size());
        ModPolynomial error(size, 0);
        for (std::size_t i = 0; i < size; ++i) {
            const BigInt target = i < polynomial.size() ? polynomial[i] : BigInt{0};
            const BigInt actual = i < product.size() ? product[i] : BigInt{0};
            const BigInt difference = target - actual;
            if (!(difference % modulus).isZero())
                return std::nullopt;
            error[i] = (difference / modulus).modulo(prime);
        }
        normalizeModPolynomial(error);

        for (std::size_t i = 0; i < lifted.size(); ++i) {
            ModPolynomial correction = remainderModPolynomial(
                multiplyModPolynomials(error, inverseCofactors[i], prime),
                modularFactors[i], prime);
            if (lifted[i].size() < modularFactors[i].size())
                lifted[i].resize(modularFactors[i].size(), BigInt{0});
            for (std::size_t j = 0; j < correction.size(); ++j)
                lifted[i][j] += modulus * BigInt::fromUnsigned(correction[j]);
        }
        const BigInt nextModulus = modulus * BigInt::fromUnsigned(prime);
        if (nextModulus.bitLength() > maximumWorkingBits)
            return std::nullopt;
        modulus = nextModulus;
    }
    return HenselLiftResult{std::move(lifted), std::move(modulus)};
}

[[nodiscard]] IntegerPolynomial modPolynomialToInteger(
    const ModPolynomial& polynomial) {
    IntegerPolynomial result;
    result.reserve(polynomial.size());
    for (const std::uint32_t coefficient : polynomial)
        result.push_back(BigInt::fromUnsigned(coefficient));
    normalizeIntegerPolynomial(result);
    return result;
}

[[nodiscard]] ModPolynomial productModularFactors(
    const std::vector<ModPolynomial>& factors,
    std::size_t begin,
    std::size_t end,
    std::uint32_t prime) {
    ModPolynomial product{1};
    for (std::size_t i = begin; i < end; ++i)
        product = multiplyModPolynomials(product, factors[i], prime);
    return product;
}

struct QuadraticPairLift final {
    IntegerPolynomial left;
    IntegerPolynomial right;
};

[[nodiscard]] std::optional<QuadraticPairLift> liftFactorPairQuadratically(
    const IntegerPolynomial& target,
    const ModPolynomial& modularLeft,
    const ModPolynomial& modularRight,
    std::uint32_t prime,
    const BigInt& targetModulus) {
    IntegerPolynomial left = modPolynomialToInteger(modularLeft);
    IntegerPolynomial right = modPolynomialToInteger(modularRight);
    auto modularInverse = inverseModPolynomial(
        remainderModPolynomial(modularRight, modularLeft, prime),
        modularLeft, prime);
    if (!modularInverse)
        return std::nullopt;
    IntegerPolynomial inverse = modPolynomialToInteger(*modularInverse);

    BigInt modulus = BigInt::fromUnsigned(prime);
    while (modulus < targetModulus) {
        const BigInt nextModulus = modulus * modulus;
        if (nextModulus > targetModulus)
            return std::nullopt;

        const IntegerPolynomial product = multiplyIntegerPolynomials(left, right);
        const std::size_t errorSize = std::max(target.size(), product.size());
        IntegerPolynomial error(errorSize, BigInt{0});
        for (std::size_t i = 0; i < errorSize; ++i) {
            const BigInt expected = i < target.size() ? target[i] : BigInt{0};
            const BigInt actual = i < product.size() ? product[i] : BigInt{0};
            const BigInt difference = expected - actual;
            if (!(difference % modulus).isZero())
                return std::nullopt;
            error[i] = positiveResidue(difference / modulus, modulus);
        }
        normalizeIntegerPolynomial(error);

        IntegerPolynomial deltaLeft = remainderMonicIntegerPolynomialModulo(
            multiplyIntegerPolynomialsModulo(error, inverse, modulus),
            left, modulus);
        IntegerPolynomial residual = subtractIntegerPolynomialsModulo(
            error,
            multiplyIntegerPolynomialsModulo(deltaLeft, right, modulus),
            modulus);
        auto [deltaRight, residualRemainder] =
            divideMonicIntegerPolynomialsModulo(
                std::move(residual), left, modulus);
        if (!isZeroIntegerPolynomial(residualRemainder)
            || deltaLeft.size() >= left.size()
            || deltaRight.size() >= right.size())
            return std::nullopt;

        IntegerPolynomial nextLeft = left;
        nextLeft.resize(std::max(nextLeft.size(), deltaLeft.size()), BigInt{0});
        for (std::size_t i = 0; i < deltaLeft.size(); ++i)
            nextLeft[i] += modulus * deltaLeft[i];
        nextLeft = reduceIntegerPolynomialModulo(
            std::move(nextLeft), nextModulus);

        IntegerPolynomial nextRight = right;
        nextRight.resize(std::max(nextRight.size(), deltaRight.size()), BigInt{0});
        for (std::size_t i = 0; i < deltaRight.size(); ++i)
            nextRight[i] += modulus * deltaRight[i];
        nextRight = reduceIntegerPolynomialModulo(
            std::move(nextRight), nextModulus);
        if (nextLeft.size() != left.size() || nextRight.size() != right.size()
            || nextLeft.back() != BigInt{1} || nextRight.back() != BigInt{1}
            || !integerPolynomialsCongruent(
                target,
                multiplyIntegerPolynomialsModulo(
                    nextLeft, nextRight, nextModulus),
                nextModulus))
            return std::nullopt;

        // inverseもNewton補正し，次のdoublingで再利用する。
        IntegerPolynomial inverseError = subtractIntegerPolynomialsModulo(
            IntegerPolynomial{BigInt{1}},
            multiplyIntegerPolynomialsModulo(inverse, nextRight, nextModulus),
            nextModulus);
        inverseError = remainderMonicIntegerPolynomialModulo(
            std::move(inverseError), nextLeft, nextModulus);
        IntegerPolynomial scaledInverseError(
            inverseError.size(), BigInt{0});
        for (std::size_t i = 0; i < inverseError.size(); ++i) {
            if (!(inverseError[i] % modulus).isZero())
                return std::nullopt;
            scaledInverseError[i] = positiveResidue(
                inverseError[i] / modulus, modulus);
        }
        normalizeIntegerPolynomial(scaledInverseError);
        IntegerPolynomial deltaInverse = remainderMonicIntegerPolynomialModulo(
            multiplyIntegerPolynomialsModulo(
                inverse, scaledInverseError, modulus),
            left, modulus);
        IntegerPolynomial nextInverse = inverse;
        nextInverse.resize(
            std::max(nextInverse.size(), deltaInverse.size()), BigInt{0});
        for (std::size_t i = 0; i < deltaInverse.size(); ++i)
            nextInverse[i] += modulus * deltaInverse[i];
        nextInverse = reduceIntegerPolynomialModulo(
            std::move(nextInverse), nextModulus);

        IntegerPolynomial verification = subtractIntegerPolynomialsModulo(
            multiplyIntegerPolynomialsModulo(
                nextInverse, nextRight, nextModulus),
            IntegerPolynomial{BigInt{1}}, nextModulus);
        verification = remainderMonicIntegerPolynomialModulo(
            std::move(verification), nextLeft, nextModulus);
        if (!isZeroIntegerPolynomial(verification))
            return std::nullopt;

        left = std::move(nextLeft);
        right = std::move(nextRight);
        inverse = std::move(nextInverse);
        modulus = nextModulus;
    }
    return QuadraticPairLift{std::move(left), std::move(right)};
}

[[nodiscard]] std::optional<std::vector<IntegerPolynomial>>
liftFactorTreeQuadratically(
    const IntegerPolynomial& target,
    const std::vector<ModPolynomial>& modularFactors,
    std::size_t begin,
    std::size_t end,
    std::uint32_t prime,
    const BigInt& targetModulus) {
    if (end - begin == 1) {
        IntegerPolynomial leaf = reduceIntegerPolynomialModulo(
            target, targetModulus);
        if (integerPolynomialModulo(leaf, prime) != modularFactors[begin]
            || leaf.back() != BigInt{1})
            return std::nullopt;
        return std::vector<IntegerPolynomial>{std::move(leaf)};
    }

    const std::size_t middle = begin + (end - begin) / 2;
    const ModPolynomial modularLeft = productModularFactors(
        modularFactors, begin, middle, prime);
    const ModPolynomial modularRight = productModularFactors(
        modularFactors, middle, end, prime);
    auto pair = liftFactorPairQuadratically(
        target, modularLeft, modularRight, prime, targetModulus);
    if (!pair)
        return std::nullopt;

    auto left = liftFactorTreeQuadratically(
        pair->left, modularFactors, begin, middle, prime, targetModulus);
    auto right = liftFactorTreeQuadratically(
        pair->right, modularFactors, middle, end, prime, targetModulus);
    if (!left || !right)
        return std::nullopt;
    left->insert(
        left->end(),
        std::make_move_iterator(right->begin()),
        std::make_move_iterator(right->end()));
    return left;
}

[[nodiscard]] std::optional<HenselLiftResult> henselLiftFactorsQuadratic(
    const IntegerPolynomial& polynomial,
    const std::vector<ModPolynomial>& modularFactors,
    std::uint32_t prime,
    std::size_t targetBits,
    std::size_t maximumWorkingBits) {
    if (modularFactors.empty())
        return std::nullopt;
    BigInt modulus = BigInt::fromUnsigned(prime);
    while (modulus.bitLength() <= targetBits) {
        const BigInt next = modulus * modulus;
        if (next.bitLength() > maximumWorkingBits)
            return std::nullopt;
        modulus = next;
    }

    auto lifted = liftFactorTreeQuadratically(
        polynomial, modularFactors, 0, modularFactors.size(),
        prime, modulus);
    if (!lifted)
        return std::nullopt;
    IntegerPolynomial product{BigInt{1}};
    for (const IntegerPolynomial& factor : *lifted)
        product = multiplyIntegerPolynomialsModulo(product, factor, modulus);
    if (!integerPolynomialsCongruent(product, polynomial, modulus))
        return std::nullopt;
    return HenselLiftResult{std::move(*lifted), std::move(modulus)};
}

[[nodiscard]] std::optional<HenselLiftResult> henselLiftFactors(
    const IntegerPolynomial& polynomial,
    const std::vector<ModPolynomial>& modularFactors,
    std::uint32_t prime,
    std::size_t targetBits,
    std::size_t maximumWorkingBits) {
    if (auto quadratic = henselLiftFactorsQuadratic(
            polynomial, modularFactors, prime,
            targetBits, maximumWorkingBits))
        return quadratic;
    return henselLiftFactorsLinear(
        polynomial, modularFactors, prime, targetBits, maximumWorkingBits);
}

[[nodiscard]] IntegerPolynomial centeredLiftedProduct(
    const std::vector<std::size_t>& indices,
    const std::vector<IntegerPolynomial>& lifted,
    const BigInt& modulus) {
    IntegerPolynomial product{BigInt{1}};
    for (const std::size_t index : indices)
        product = multiplyIntegerPolynomialsModulo(product, lifted[index], modulus);
    const BigInt half = modulus / BigInt{2};
    for (BigInt& coefficient : product)
        if (coefficient > half)
            coefficient -= modulus;
    normalizeIntegerPolynomial(product);
    return product;
}

using FactorDegreeMask = std::vector<bool>;

[[nodiscard]] FactorDegreeMask modularFactorDegreeMask(
    const std::vector<ModPolynomial>& factors,
    std::size_t totalDegree) {
    FactorDegreeMask possible(totalDegree + 1, false);
    possible[0] = true;
    for (const ModPolynomial& factor : factors) {
        const std::size_t degree = factor.size() - 1;
        for (std::size_t current = totalDegree + 1; current-- > degree;)
            possible[current] = possible[current] || possible[current - degree];
    }
    return possible;
}

[[nodiscard]] bool hasProperFactorDegree(
    const FactorDegreeMask& possible) noexcept {
    if (possible.size() <= 2)
        return false;
    const std::size_t totalDegree = possible.size() - 1;
    for (std::size_t degree = 1; degree <= totalDegree / 2; ++degree)
        if (possible[degree])
            return true;
    return false;
}

[[nodiscard]] std::size_t estimatedAllowedSubsets(
    const std::vector<ModPolynomial>& factors,
    const FactorDegreeMask& allowed,
    std::size_t saturation) {
    std::vector<std::size_t> counts(allowed.size(), 0);
    counts[0] = 1;
    for (const ModPolynomial& factor : factors) {
        const std::size_t degree = factor.size() - 1;
        for (std::size_t current = counts.size(); current-- > degree;) {
            const std::size_t addition = counts[current - degree];
            counts[current] = std::min(
                saturation,
                counts[current] > saturation - std::min(addition, saturation)
                    ? saturation
                    : counts[current] + addition);
        }
    }
    std::size_t result = 0;
    const std::size_t totalDegree = allowed.size() - 1;
    for (std::size_t degree = 1; degree <= totalDegree / 2; ++degree) {
        if (!allowed[degree])
            continue;
        result = std::min(
            saturation,
            result > saturation - std::min(counts[degree], saturation)
                ? saturation
                : result + counts[degree]);
    }
    return result;
}

[[nodiscard]] std::size_t ceilLog2Size(std::size_t value) noexcept {
    if (value <= 1)
        return 0;
    --value;
    std::size_t result = 0;
    while (value != 0) {
        ++result;
        value >>= 1;
    }
    return result;
}

[[nodiscard]] std::optional<std::size_t> cldCoefficientBoundBits(
    std::size_t factorCoefficientBoundBits,
    std::size_t degree) noexcept {
    // g,hの各係数をBで押さえると，coeff(h*g')は高々degree^2*B^2。
    // Mignotte boundは余裕を含むため，さらに2 bitをguardへ足す。
    const std::size_t degreeBits = ceilLog2Size(degree + 1);
    constexpr std::size_t guardBits = 2;
    if (factorCoefficientBoundBits
        > (std::numeric_limits<std::size_t>::max()
            - 2 * degreeBits - guardBits) / 2)
        return std::nullopt;
    return 2 * factorCoefficientBoundBits
        + 2 * degreeBits + guardBits;
}

[[nodiscard]] IntegerPolynomial integerPolynomialDerivativeModulo(
    const IntegerPolynomial& polynomial,
    const BigInt& modulus) {
    if (polynomial.size() <= 1)
        return IntegerPolynomial{BigInt{0}};
    IntegerPolynomial result(polynomial.size() - 1, BigInt{0});
    for (std::size_t exponent = 1; exponent < polynomial.size(); ++exponent) {
        result[exponent - 1] = positiveResidue(
            polynomial[exponent]
                * BigInt::fromUnsigned(static_cast<std::uint64_t>(exponent)),
            modulus);
    }
    normalizeIntegerPolynomial(result);
    return result;
}

[[nodiscard]] BigInt centeredResidue(
    BigInt value,
    const BigInt& modulus) {
    value = positiveResidue(std::move(value), modulus);
    if (value > modulus / BigInt{2})
        value -= modulus;
    return value;
}

[[nodiscard]] std::optional<std::vector<IntegerPolynomial>>
coefficientLogarithmicDerivativeRows(
    const IntegerPolynomial& polynomial,
    const std::vector<std::size_t>& indices,
    const std::vector<IntegerPolynomial>& lifted,
    const BigInt& modulus) {
    std::vector<IntegerPolynomial> result;
    result.reserve(indices.size());
    for (const std::size_t index : indices) {
        if (index >= lifted.size() || lifted[index].back() != BigInt{1})
            return std::nullopt;
        auto [cofactor, remainder] = divideMonicIntegerPolynomialsModulo(
            polynomial, lifted[index], modulus);
        if (!isZeroIntegerPolynomial(remainder))
            return std::nullopt;
        IntegerPolynomial derivative = integerPolynomialDerivativeModulo(
            lifted[index], modulus);
        IntegerPolynomial cld = multiplyIntegerPolynomialsModulo(
            cofactor, derivative, modulus);
        cld.resize(polynomial.size() - 1, BigInt{0});
        result.push_back(std::move(cld));
    }
    return result;
}

[[nodiscard]] std::vector<std::size_t> cldCoefficientWindow(
    std::size_t degree,
    std::size_t count,
    std::size_t pass) {
    std::vector<std::size_t> result;
    result.reserve(count);
    if (pass == 0) {
        for (std::size_t index = 0; index < count; ++index)
            result.push_back(degree - 1 - index);
        return result;
    }
    if (pass == 1) {
        for (std::size_t index = 0; index < count; ++index)
            result.push_back(index);
        return result;
    }

    const std::size_t shift = (pass - 2) % degree;
    for (std::size_t index = 0; index < count; ++index)
        result.push_back(((index * degree) / count + shift) % degree);
    return result;
}

struct ExactRecombinationSplit final {
    IntegerPolynomial factor;
    IntegerPolynomial quotient;
    std::vector<std::size_t> selected;
    std::vector<std::size_t> complement;
};

struct CldSearchBudget final {
    std::size_t lattices = 0;
    std::size_t candidates = 0;
};

[[nodiscard]] bool sameCldSignature(
    std::size_t lhs,
    std::size_t rhs,
    const linear_algebra::IntegerLatticeBasis& reduced,
    std::span<const std::size_t> shortRows) {
    for (const std::size_t row : shortRows)
        if (reduced[row][lhs] != reduced[row][rhs])
            return false;
    return true;
}

void appendUniqueIndexSet(
    std::vector<std::vector<std::size_t>>& candidates,
    std::vector<std::size_t> candidate,
    std::size_t factorCount) {
    std::sort(candidate.begin(), candidate.end());
    candidate.erase(
        std::unique(candidate.begin(), candidate.end()), candidate.end());
    if (candidate.empty() || candidate.size() == factorCount)
        return;
    if (std::find(candidates.begin(), candidates.end(), candidate)
        == candidates.end())
        candidates.push_back(std::move(candidate));
}

[[nodiscard]] std::vector<std::vector<std::size_t>>
cldCandidateIndexSets(
    const linear_algebra::IntegerLatticeBasis& reduced,
    std::size_t factorCount,
    const std::vector<IntegerPolynomial>& cldRows,
    const BigInt& identityWeight,
    const BigInt& modulus) {
    std::vector<std::size_t> shortRows;
    const BigInt modulusSquared = modulus * modulus;
    const BigInt normScale = BigInt::fromUnsigned(
        static_cast<std::uint64_t>(
            16 * (factorCount + cldRows.front().size())));
    for (std::size_t row = 0; row < reduced.size(); ++row) {
        bool hasFactorCoordinate = false;
        bool validFactorCoordinates = true;
        std::vector<BigInt> combination(factorCount);
        BigInt normSquared;
        for (std::size_t factor = 0; factor < factorCount; ++factor) {
            if (!reduced[row][factor].isZero())
                hasFactorCoordinate = true;
            normSquared += reduced[row][factor] * reduced[row][factor];
            auto divided = numeric::divmod(
                reduced[row][factor], identityWeight);
            if (!divided.remainder.isZero()) {
                validFactorCoordinates = false;
                break;
            }
            combination[factor] = std::move(divided.quotient);
        }
        if (!hasFactorCoordinate || !validFactorCoordinates)
            continue;

        // 選んだ係数窓だけの偶然の短関係をblock情報へ混ぜない。
        // 同じcombinationを全CLD係数で再評価し，Mに比べて全体が短いものだけを使う。
        bool globallyShort = normScale * normSquared < modulusSquared;
        for (std::size_t coefficient = 0;
             globallyShort && coefficient < cldRows.front().size();
             ++coefficient) {
            BigInt value;
            for (std::size_t factor = 0; factor < factorCount; ++factor) {
                if (combination[factor].isZero())
                    continue;
                value += positiveResidue(
                    combination[factor], modulus)
                    * cldRows[factor][coefficient];
                value = positiveResidue(std::move(value), modulus);
            }
            value = centeredResidue(std::move(value), modulus);
            normSquared += value * value;
            globallyShort = normScale * normSquared < modulusSquared;
        }
        if (globallyShort)
            shortRows.push_back(row);
    }
    if (shortRows.empty())
        return {};

    std::vector<std::vector<std::size_t>> candidates;
    std::vector<bool> grouped(factorCount, false);
    for (std::size_t factor = 0; factor < factorCount; ++factor) {
        if (grouped[factor])
            continue;
        std::vector<std::size_t> group{factor};
        grouped[factor] = true;
        for (std::size_t other = factor + 1; other < factorCount; ++other) {
            if (!grouped[other]
                && sameCldSignature(
                    factor, other, reduced, shortRows)) {
                grouped[other] = true;
                group.push_back(other);
            }
        }
        appendUniqueIndexSet(candidates, std::move(group), factorCount);
    }

    // 短基底がblock indicatorそのものを含む場合と，二blockの差を含む場合の
    // 両方を拾う。false positiveは後段のexact divisionで棄却する。
    for (const std::size_t row : shortRows) {
        std::vector<std::size_t> positive;
        std::vector<std::size_t> negative;
        std::vector<std::size_t> support;
        for (std::size_t factor = 0; factor < factorCount; ++factor) {
            if (reduced[row][factor].isPositive())
                positive.push_back(factor);
            if (reduced[row][factor].isNegative())
                negative.push_back(factor);
            if (!reduced[row][factor].isZero())
                support.push_back(factor);
        }
        appendUniqueIndexSet(candidates, std::move(positive), factorCount);
        appendUniqueIndexSet(candidates, std::move(negative), factorCount);
        appendUniqueIndexSet(candidates, std::move(support), factorCount);
    }
    return candidates;
}

[[nodiscard]] std::optional<ExactRecombinationSplit>
verifyLiftedIndexSet(
    const IntegerPolynomial& polynomial,
    const std::vector<std::size_t>& indices,
    std::span<const std::size_t> localSelection,
    const std::vector<IntegerPolynomial>& lifted,
    const BigInt& modulus,
    const FactorDegreeMask& allowedDegrees,
    const RationalPolynomialFactorOptions& options,
    CldSearchBudget& budget) {
    if (budget.candidates >= options.maximumCldCandidates)
        return std::nullopt;
    ++budget.candidates;

    std::vector<bool> selectedLocally(indices.size(), false);
    std::vector<std::size_t> selected;
    selected.reserve(localSelection.size());
    for (const std::size_t local : localSelection) {
        if (local >= indices.size())
            return std::nullopt;
        selectedLocally[local] = true;
        selected.push_back(indices[local]);
    }
    std::vector<std::size_t> complement;
    for (std::size_t local = 0; local < indices.size(); ++local)
        if (!selectedLocally[local])
            complement.push_back(indices[local]);
    if (selected.empty() || complement.empty())
        return std::nullopt;

    std::size_t selectedDegree = 0;
    for (const std::size_t index : selected)
        selectedDegree += lifted[index].size() - 1;
    const std::size_t polynomialDegree = polynomial.size() - 1;
    if (selectedDegree == 0 || selectedDegree >= polynomialDegree
        || selectedDegree >= allowedDegrees.size()
        || polynomialDegree - selectedDegree >= allowedDegrees.size()
        || !allowedDegrees[selectedDegree]
        || !allowedDegrees[polynomialDegree - selectedDegree])
        return std::nullopt;

    IntegerPolynomial candidate = centeredLiftedProduct(
        selected, lifted, modulus);
    if (candidate.size() <= 1 || candidate.size() >= polynomial.size())
        return std::nullopt;
    auto quotient = divideIntegerPolynomialsExact(polynomial, candidate);
    if (!quotient || quotient->size() <= 1)
        return std::nullopt;
    return ExactRecombinationSplit{
        std::move(candidate), std::move(*quotient),
        std::move(selected), std::move(complement)};
}

[[nodiscard]] std::optional<ExactRecombinationSplit>
findCldRecombinationSplit(
    const IntegerPolynomial& polynomial,
    const std::vector<std::size_t>& indices,
    const std::vector<IntegerPolynomial>& lifted,
    const BigInt& modulus,
    const FactorDegreeMask& allowedDegrees,
    const RationalPolynomialFactorOptions& options,
    CldSearchBudget& budget) {
    const std::size_t factorCount = indices.size();
    const std::size_t degree = polynomial.size() - 1;
    if (factorCount < options.minimumCldModularFactors
        || budget.lattices >= options.maximumCldLattices
        || budget.candidates >= options.maximumCldCandidates)
        return std::nullopt;

    const std::size_t dimensionLimit = std::min(
        options.maximumLllRank, options.maximumLllColumns);
    if (factorCount >= dimensionLimit)
        return std::nullopt;
    const std::size_t coefficientCount = std::min({
        degree,
        dimensionLimit - factorCount,
        options.maximumCldCoefficientColumns});
    if (coefficientCount == 0)
        return std::nullopt;

    auto cldRows = coefficientLogarithmicDerivativeRows(
        polynomial, indices, lifted, modulus);
    if (!cldRows)
        return std::nullopt;

    // indicator座標を1のまま置くと，M-bitのCLD列とのscale差だけで
    // LLL swapが増える。M/(64*dimension)以下の2冪へ持ち上げ，真のblock
    // vectorがshortである余地を保ちながらbasisをpreconditionする。
    const std::size_t identityGuardBits = ceilLog2Size(
        factorCount + degree) + 6;
    const std::size_t identityBits = modulus.bitLength() > identityGuardBits
        ? modulus.bitLength() - identityGuardBits
        : 0;
    const BigInt identityWeight = BigInt{1} << identityBits;

    const std::size_t passLimit = options.maximumCldLattices - budget.lattices;
    for (std::size_t pass = 0; pass < passLimit; ++pass) {
        const std::vector<std::size_t> coefficientIndices =
            cldCoefficientWindow(degree, coefficientCount, pass);
        const std::size_t dimension = factorCount + coefficientCount;
        linear_algebra::IntegerLatticeBasis lattice(
            dimension,
            linear_algebra::IntegerLatticeRow(dimension, BigInt{0}));
        for (std::size_t factor = 0; factor < factorCount; ++factor) {
            lattice[factor][factor] = identityWeight;
            for (std::size_t column = 0;
                 column < coefficientCount; ++column) {
                lattice[factor][factorCount + column] = centeredResidue(
                    (*cldRows)[factor][coefficientIndices[column]], modulus);
            }
        }
        for (std::size_t column = 0; column < coefficientCount; ++column)
            lattice[factorCount + column][factorCount + column] = modulus;

        linear_algebra::ExactLllOptions lllOptions;
        lllOptions.maximumRank = options.maximumLllRank;
        lllOptions.maximumColumns = options.maximumLllColumns;
        lllOptions.maximumIntermediateBits =
            options.maximumLllIntermediateBits;
        lllOptions.maximumSwaps = options.maximumLllSwaps;
        lllOptions.maximumSizeReductions =
            options.maximumLllSizeReductions;
        ++budget.lattices;
        const auto reduced = linear_algebra::exactLllReduceRows(
            std::move(lattice), lllOptions);
        const bool usableBudgetStop =
            reduced.status == linear_algebra::ExactLllStatus::SwapLimitExceeded
            || reduced.status
                == linear_algebra::ExactLllStatus::SizeReductionLimitExceeded;
        if (!reduced.reduced() && !usableBudgetStop)
            return std::nullopt;

        auto candidates = cldCandidateIndexSets(
            reduced.basis, factorCount, *cldRows,
            identityWeight, modulus);
        std::stable_sort(
            candidates.begin(), candidates.end(),
            [&](const auto& lhs, const auto& rhs) {
                const auto degreeOf = [&](const auto& selection) {
                    std::size_t result = 0;
                    for (const std::size_t local : selection)
                        result += lifted[indices[local]].size() - 1;
                    return result;
                };
                return degreeOf(lhs) < degreeOf(rhs);
            });
        for (const auto& candidate : candidates) {
            if (auto split = verifyLiftedIndexSet(
                    polynomial, indices, candidate, lifted, modulus,
                    allowedDegrees, options, budget))
                return split;
            if (budget.candidates >= options.maximumCldCandidates)
                return std::nullopt;
        }
        if (!reduced.reduced())
            return std::nullopt;
    }
    return std::nullopt;
}

struct RecombineResult final {
    std::vector<IntegerPolynomial> factors;
    bool complete = false;
};

[[nodiscard]] RecombineResult recombineLiftedFactorsRecursive(
    const IntegerPolynomial& polynomial,
    const std::vector<std::size_t>& indices,
    const std::vector<IntegerPolynomial>& lifted,
    const BigInt& modulus,
    const FactorDegreeMask& allowedDegrees,
    const RationalPolynomialFactorOptions& options,
    bool useCld,
    CldSearchBudget& cldBudget,
    std::size_t& combinations) {
    if (indices.size() <= 1 || polynomial.size() <= 2) {
        return RecombineResult{{polynomial}, true};
    }

    if (useCld) {
        if (auto split = findCldRecombinationSplit(
                polynomial, indices, lifted, modulus,
                allowedDegrees, options, cldBudget)) {
            RecombineResult left = recombineLiftedFactorsRecursive(
                split->factor, split->selected, lifted, modulus,
                allowedDegrees, options, useCld, cldBudget, combinations);
            RecombineResult right = recombineLiftedFactorsRecursive(
                split->quotient, split->complement, lifted, modulus,
                allowedDegrees, options, useCld, cldBudget, combinations);
            left.factors.insert(
                left.factors.end(),
                std::make_move_iterator(right.factors.begin()),
                std::make_move_iterator(right.factors.end()));
            left.complete = left.complete && right.complete;
            return left;
        }
    }

    for (std::size_t count = 1; count <= indices.size() / 2; ++count) {
        std::vector<std::size_t> selected;
        std::optional<RecombineResult> splitResult;
        bool budgetExceeded = false;
        std::function<void(std::size_t, std::size_t)> search =
            [&](std::size_t position, std::size_t remaining) {
                if (splitResult || budgetExceeded)
                    return;
                if (remaining == 0) {
                    if (++combinations > options.maximumCombinations) {
                        budgetExceeded = true;
                        return;
                    }
                    std::size_t selectedDegree = 0;
                    for (const std::size_t index : selected)
                        selectedDegree += lifted[index].size() - 1;
                    const std::size_t polynomialDegree = polynomial.size() - 1;
                    if (selectedDegree == 0 || selectedDegree >= polynomialDegree
                        || selectedDegree >= allowedDegrees.size()
                        || polynomialDegree - selectedDegree >= allowedDegrees.size()
                        || !allowedDegrees[selectedDegree]
                        || !allowedDegrees[polynomialDegree - selectedDegree])
                        return;
                    IntegerPolynomial candidate = centeredLiftedProduct(
                        selected, lifted, modulus);
                    if (candidate.size() <= 1 || candidate.size() >= polynomial.size())
                        return;
                    auto quotient = divideIntegerPolynomialsExact(polynomial, candidate);
                    if (!quotient || quotient->size() <= 1)
                        return;

                    std::vector<bool> chosen(lifted.size(), false);
                    for (const std::size_t index : selected)
                        chosen[index] = true;
                    std::vector<std::size_t> complement;
                    for (const std::size_t index : indices)
                        if (!chosen[index])
                            complement.push_back(index);
                    RecombineResult left = recombineLiftedFactorsRecursive(
                        candidate, selected, lifted, modulus, allowedDegrees,
                        options, useCld, cldBudget, combinations);
                    RecombineResult right = recombineLiftedFactorsRecursive(
                        *quotient, complement, lifted, modulus, allowedDegrees,
                        options, useCld, cldBudget, combinations);
                    left.factors.insert(
                        left.factors.end(),
                        std::make_move_iterator(right.factors.begin()),
                        std::make_move_iterator(right.factors.end()));
                    left.complete = left.complete && right.complete;
                    splitResult = std::move(left);
                    return;
                }
                if (indices.size() - position < remaining)
                    return;
                selected.push_back(indices[position]);
                search(position + 1, remaining - 1);
                selected.pop_back();
                search(position + 1, remaining);
            };
        search(0, count);
        if (splitResult)
            return std::move(*splitResult);
        if (budgetExceeded)
            return RecombineResult{{polynomial}, false};
    }

    // lifted irreducible factorsのどのproper subsetもexact因子を与えなかった。
    return RecombineResult{{polynomial}, true};
}

struct ModularFactorizationResult final {
    std::vector<RationalPolynomial> factors;
    bool complete = false;
};

[[nodiscard]] std::optional<ModularFactorizationResult>
factorSquareFreeRationalPolynomialModular(
    const RationalPolynomial& polynomial,
    const RationalPolynomialFactorOptions& options) {
    if (polynomial.degree() <= 1)
        return ModularFactorizationResult{
            {monicRationalPolynomial(polynomial)}, true};
    if (polynomial.degree() > options.maximumFiniteFieldDegree)
        return std::nullopt;

    BigInt scale;
    auto transformed = monicIntegerTransform(
        polynomial, scale, options.maximumLiftBits);
    if (!transformed)
        return std::nullopt;
    const std::size_t coefficientBoundBits = mignotteCoefficientBoundBits(*transformed);
    if (coefficientBoundBits + 1 > options.maximumLiftBits)
        return std::nullopt;

    constexpr std::array<std::uint32_t, 24> primes{
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
        43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    struct GoodPrimeFactorization final {
        std::uint32_t prime = 0;
        std::vector<ModPolynomial> factors;
    };
    std::vector<GoodPrimeFactorization> goodPrimes;
    FactorDegreeMask allowedDegrees(polynomial.degree() + 1, true);
    const std::size_t trialCount = std::min(options.maximumPrimeTrials, primes.size());
    for (std::size_t trial = 0;
         trial < trialCount && goodPrimes.size() < options.maximumGoodPrimes;
         ++trial) {
        const std::uint32_t prime = primes[trial];
        ModPolynomial reduced = integerPolynomialModulo(*transformed, prime);
        if (reduced.size() != transformed->size())
            continue;
        ModPolynomial derivative = derivativeModPolynomial(reduced, prime);
        if (derivative.empty()
            || gcdModPolynomial(reduced, derivative, prime).size() != 1)
            continue;
        auto modularFactors = factorFiniteFieldPolynomial(
            reduced, prime, options);
        if (!modularFactors)
            continue;
        if (modularFactors->size() == 1)
            return ModularFactorizationResult{
                {monicRationalPolynomial(polynomial)}, true};

        const FactorDegreeMask primeDegrees = modularFactorDegreeMask(
            *modularFactors, polynomial.degree());
        for (std::size_t degree = 0; degree < allowedDegrees.size(); ++degree)
            allowedDegrees[degree] = allowedDegrees[degree] && primeDegrees[degree];
        goodPrimes.push_back(GoodPrimeFactorization{
            prime, std::move(*modularFactors)});
        if (!hasProperFactorDegree(allowedDegrees)) {
            return ModularFactorizationResult{
                {monicRationalPolynomial(polynomial)}, true};
        }
        if (goodPrimes.size() >= 2) {
            const std::size_t scanSaturation =
                options.preferredRecombinationCandidates
                    == std::numeric_limits<std::size_t>::max()
                ? options.preferredRecombinationCandidates
                : options.preferredRecombinationCandidates + 1;
            std::size_t bestCandidateCount = scanSaturation;
            for (const GoodPrimeFactorization& candidate : goodPrimes) {
                bestCandidateCount = std::min(
                    bestCandidateCount,
                    estimatedAllowedSubsets(
                        candidate.factors, allowedDegrees, scanSaturation));
            }
            // 追加primeの価値が小さい領域ではlift/recombinationへ進む。
            // exact divisionが最終certificateなので，これは性能上の停止条件であり
            // 正しさやcomplete判定を弱めない。
            if (bestCandidateCount
                <= options.preferredRecombinationCandidates)
                break;
        }
    }
    if (goodPrimes.empty())
        return std::nullopt;

    const std::size_t saturation = options.maximumCombinations == 0
        ? 1
        : options.maximumCombinations;
    const auto best = std::min_element(
        goodPrimes.begin(), goodPrimes.end(),
        [&](const GoodPrimeFactorization& lhs,
            const GoodPrimeFactorization& rhs) {
            const std::size_t leftCandidates = estimatedAllowedSubsets(
                lhs.factors, allowedDegrees, saturation);
            const std::size_t rightCandidates = estimatedAllowedSubsets(
                rhs.factors, allowedDegrees, saturation);
            if (leftCandidates != rightCandidates)
                return leftCandidates < rightCandidates;
            if (lhs.factors.size() != rhs.factors.size())
                return lhs.factors.size() < rhs.factors.size();
            return lhs.prime < rhs.prime;
        });

    auto lifted = henselLiftFactors(
        *transformed, best->factors, best->prime,
        coefficientBoundBits + 1, options.maximumHenselWorkingBits);
    if (!lifted)
        return std::nullopt;

    const std::size_t cldSaturation =
        options.preferredRecombinationCandidates
            == std::numeric_limits<std::size_t>::max()
        ? options.preferredRecombinationCandidates
        : options.preferredRecombinationCandidates + 1;
    const std::size_t estimatedCandidates = estimatedAllowedSubsets(
        best->factors, allowedDegrees, cldSaturation);
    const bool cldWorthwhile =
        best->factors.size() >= options.minimumCldModularFactors
        && estimatedCandidates > options.preferredRecombinationCandidates;
    bool useCld = cldWorthwhile
        && options.maximumCldLattices != 0
        && options.maximumCldCandidates != 0;
    const std::size_t cldCoefficientColumns = std::min({
        polynomial.degree(),
        options.maximumCldCoefficientColumns,
        options.maximumLllRank > best->factors.size()
            ? options.maximumLllRank - best->factors.size()
            : std::size_t{0},
        options.maximumLllColumns > best->factors.size()
            ? options.maximumLllColumns - best->factors.size()
            : std::size_t{0}});
    const std::size_t cldDimension =
        best->factors.size() + cldCoefficientColumns;
    // 小dimensionではrigorous CLD boundまで上げても安い。大きいlatticeは，
    // 既存Hensel精度で先にprobeし，exact divisionによってsoundnessを閉じる。
    // 高精度化でLLL swap数まで増幅する逆効果を避ける。
    if (useCld && cldCoefficientColumns != 0 && cldDimension <= 12) {
        const auto coefficientBits = cldCoefficientBoundBits(
            coefficientBoundBits, polynomial.degree());
        // coefficient windowだけでなく全CLD vectorを短関係として再検証するため，
        // full degreeのEuclidean normにも余裕を持たせる。
        const std::size_t extraBits = ceilLog2Size(std::max(
            polynomial.degree(), options.maximumLllRank)) + 6;
        if (coefficientBits
            && *coefficientBits <= options.maximumCldLiftBits
            && extraBits <= options.maximumCldLiftBits - *coefficientBits) {
            const std::size_t cldTargetBits = *coefficientBits + extraBits;
            if (lifted->modulus.bitLength() > cldTargetBits) {
                useCld = true;
            }
            else if (auto cldLifted = henselLiftFactors(
                    *transformed, best->factors, best->prime,
                    cldTargetBits, options.maximumHenselWorkingBits)) {
                lifted = std::move(cldLifted);
                useCld = true;
            }
        }
    }
    else if (cldCoefficientColumns == 0) {
        useCld = false;
    }

    std::vector<std::size_t> indices(lifted->factors.size());
    for (std::size_t i = 0; i < indices.size(); ++i)
        indices[i] = i;
    std::size_t combinations = 0;
    CldSearchBudget cldBudget;
    RecombineResult recombined = recombineLiftedFactorsRecursive(
        *transformed, indices, lifted->factors, lifted->modulus,
        allowedDegrees, options, useCld, cldBudget, combinations);

    std::vector<RationalPolynomial> result;
    result.reserve(recombined.factors.size());
    RationalPolynomial product{{one()}};
    for (const IntegerPolynomial& factor : recombined.factors) {
        RationalPolynomial rational = inverseMonicIntegerTransform(factor, scale);
        product = multiplyRationalPolynomialValues(product, rational);
        result.push_back(std::move(rational));
    }
    if (product.coefficients()
        != monicRationalPolynomial(polynomial).coefficients())
        return std::nullopt;
    return ModularFactorizationResult{
        std::move(result), recombined.complete};
}

[[nodiscard]] BigInt evaluateIntegerCoefficients(
    const std::vector<BigInt>& coefficients,
    std::int64_t value) {
    BigInt result{0};
    const BigInt point{value};
    for (auto iterator = coefficients.rbegin(); iterator != coefficients.rend(); ++iterator)
        result = result * point + *iterator;
    return result;
}

[[nodiscard]] std::optional<std::vector<std::uint64_t>> factoredPositiveDivisors(
    const BigInt& value,
    std::size_t maximumDivisors) {
    const auto magnitude = numeric::tryToUint64(value.abs());
    if (!magnitude || *magnitude == 0)
        return std::nullopt;

    std::vector<std::uint64_t> primes;
    if (!numeric::factorUint64(*magnitude, primes))
        return std::nullopt;
    std::sort(primes.begin(), primes.end());
    std::vector<std::pair<std::uint64_t, std::size_t>> groups;
    for (const std::uint64_t prime : primes) {
        if (!groups.empty() && groups.back().first == prime)
            ++groups.back().second;
        else
            groups.emplace_back(prime, 1);
    }

    std::vector<std::uint64_t> divisors{1};
    for (const auto& [prime, exponent] : groups) {
        const std::size_t oldSize = divisors.size();
        std::uint64_t primePower = 1;
        for (std::size_t e = 1; e <= exponent; ++e) {
            if (primePower > std::numeric_limits<std::uint64_t>::max() / prime)
                return std::nullopt;
            primePower *= prime;
            if (oldSize > maximumDivisors - divisors.size())
                return std::nullopt;
            for (std::size_t i = 0; i < oldSize; ++i) {
                if (divisors[i] > std::numeric_limits<std::uint64_t>::max() / primePower)
                    return std::nullopt;
                divisors.push_back(divisors[i] * primePower);
            }
        }
    }
    std::sort(divisors.begin(), divisors.end());
    return divisors;
}

[[nodiscard]] std::optional<BigInt> exactIntegerNthRoot(
    const BigInt& value,
    std::size_t degree) {
    if (degree == 0)
        return std::nullopt;
    if (degree == 1 || value.isZero())
        return value;
    const bool negative = value.isNegative();
    if (negative && degree % 2 == 0)
        return std::nullopt;
    const BigInt magnitude = value.abs();
    const std::size_t rootBits = (magnitude.bitLength() + degree - 1) / degree;
    BigInt low{0};
    BigInt high = BigInt{1} << (rootBits + 1);
    while (high - low > BigInt{1}) {
        const BigInt middle = (low + high) / BigInt{2};
        if (numeric::pow(middle, static_cast<std::uint64_t>(degree)) <= magnitude)
            low = middle;
        else
            high = middle;
    }
    if (numeric::pow(low, static_cast<std::uint64_t>(degree)) != magnitude)
        return std::nullopt;
    return negative ? -low : low;
}

[[nodiscard]] std::optional<Rational> exactRationalNthRoot(
    const Rational& value,
    std::size_t degree) {
    const auto numerator = exactIntegerNthRoot(value.numerator(), degree);
    const auto denominator = exactIntegerNthRoot(value.denominator(), degree);
    if (!numerator || !denominator)
        return std::nullopt;
    return Rational{*numerator, *denominator};
}

struct CertifiedPolynomialSplit final {
    std::optional<std::pair<RationalPolynomial, RationalPolynomial>> factors;
    bool irreducible = false;
};

[[nodiscard]] std::vector<std::size_t> distinctPrimeDivisors(
    std::size_t value) {
    std::vector<std::size_t> result;
    for (std::size_t divisor = 2; divisor <= value / divisor; ++divisor) {
        if (value % divisor != 0)
            continue;
        result.push_back(divisor);
        while (value % divisor == 0)
            value /= divisor;
    }
    if (value > 1)
        result.push_back(value);
    return result;
}

[[nodiscard]] CertifiedPolynomialSplit splitRationalBinomial(
    const RationalPolynomial& input) {
    const RationalPolynomial polynomial = monicRationalPolynomial(input);
    if (polynomial.degree() < 2 || polynomial.coefficient(0).isZero())
        return {};
    for (std::size_t exponent = 1; exponent < polynomial.degree(); ++exponent)
        if (!polynomial.coefficient(exponent).isZero())
            return {};

    // Capelliの判定をQ上でexactに適用する。x^n-aは，nの素因数qに
    // 対してaがq乗，または4|nかつa=-4b^4の場合に限り可約である。
    const Rational radicand = -polynomial.coefficient(0);
    for (const std::size_t prime : distinctPrimeDivisors(polynomial.degree())) {
        const auto root = exactRationalNthRoot(radicand, prime);
        if (!root)
            continue;
        const std::size_t reducedDegree = polynomial.degree() / prime;
        std::vector<Rational> coefficients(reducedDegree + 1, zero());
        coefficients[0] = -*root;
        coefficients[reducedDegree] = one();
        RationalPolynomial factor{std::move(coefficients)};
        auto quotient = divideRationalPolynomialsExact(polynomial, factor);
        if (quotient && quotient->degree() > 0)
            return CertifiedPolynomialSplit{
                std::pair{
                    monicRationalPolynomial(factor),
                    monicRationalPolynomial(*quotient)},
                false};
    }

    if (polynomial.degree() % 4 == 0) {
        const auto fourthRoot = exactRationalNthRoot(
            (-radicand) / Rational{BigInt{4}}, 4);
        if (fourthRoot) {
            const std::size_t quarter = polynomial.degree() / 4;
            std::vector<Rational> coefficients(2 * quarter + 1, zero());
            const Rational two{BigInt{2}};
            coefficients[0] = two * *fourthRoot * *fourthRoot;
            coefficients[quarter] = -two * *fourthRoot;
            coefficients[2 * quarter] = one();
            RationalPolynomial factor{std::move(coefficients)};
            auto quotient = divideRationalPolynomialsExact(polynomial, factor);
            if (quotient && quotient->degree() > 0)
                return CertifiedPolynomialSplit{
                    std::pair{
                        monicRationalPolynomial(factor),
                        monicRationalPolynomial(*quotient)},
                    false};
        }
    }
    return CertifiedPolynomialSplit{std::nullopt, true};
}

[[nodiscard]] std::optional<std::pair<RationalPolynomial, RationalPolynomial>>
boundedSparseBinomialSplit(
    const RationalPolynomial& input,
    const RationalPolynomialFactorOptions& options) {
    const RationalPolynomial polynomial = monicRationalPolynomial(input);
    std::vector<std::size_t> support;
    for (std::size_t exponent = 0; exponent <= polynomial.degree(); ++exponent)
        if (!polynomial.coefficient(exponent).isZero())
            support.push_back(exponent);
    if (support.size() < 3 || support.size() > options.maximumSparseTerms)
        return std::nullopt;

    std::vector<bool> seenDegree(polynomial.degree() / 2 + 1, false);
    std::vector<std::size_t> candidateDegrees;
    for (std::size_t i = 0; i < support.size(); ++i) {
        for (std::size_t j = i + 1; j < support.size(); ++j) {
            const std::size_t difference = support[j] - support[i];
            if (difference == 0 || difference >= seenDegree.size()
                || seenDegree[difference])
                continue;
            seenDegree[difference] = true;
            candidateDegrees.push_back(difference);
        }
    }
    std::sort(candidateDegrees.begin(), candidateDegrees.end());

    // supportがHとH+dの重ならない2層で，対応係数比が一定なら
    // f=(1+r*x^d)hである。一般候補の多項式除算を繰り返す前にO(|S|^2)
    // のsupport検査と1回のexact divisionだけで回収する。
    std::vector<bool> present(polynomial.degree() + 1, false);
    for (const std::size_t exponent : support)
        present[exponent] = true;
    for (const std::size_t degree : candidateDegrees) {
        std::vector<bool> consumed(polynomial.degree() + 1, false);
        std::optional<Rational> ratio;
        bool paired = true;
        for (const std::size_t exponent : support) {
            if (consumed[exponent])
                continue;
            if (exponent > polynomial.degree() - degree
                || !present[exponent + degree]) {
                paired = false;
                break;
            }
            const Rational currentRatio = polynomial.coefficient(exponent + degree)
                / polynomial.coefficient(exponent);
            if (ratio && currentRatio != *ratio) {
                paired = false;
                break;
            }
            ratio = currentRatio;
            consumed[exponent] = true;
            consumed[exponent + degree] = true;
        }
        if (!paired || !ratio || ratio->isZero())
            continue;
        std::vector<Rational> coefficients(degree + 1, zero());
        coefficients[0] = one();
        coefficients[degree] = *ratio;
        RationalPolynomial factor = monicRationalPolynomial(
            RationalPolynomial{std::move(coefficients)});
        auto quotient = divideRationalPolynomialsExact(polynomial, factor);
        if (quotient && quotient->degree() > 0)
            return std::pair{
                std::move(factor), monicRationalPolynomial(*quotient)};
    }

    const std::vector<BigInt> primitive = primitiveIntegerCoefficients(polynomial);
    const auto leadingDivisors = factoredPositiveDivisors(
        primitive.back(), options.maximumDivisorsPerSample);
    const auto constantDivisors = factoredPositiveDivisors(
        primitive.front(), options.maximumDivisorsPerSample);
    if (!leadingDivisors || !constantDivisors)
        return std::nullopt;

    std::size_t candidates = 0;
    for (const std::size_t degree : candidateDegrees) {
        for (const std::uint64_t leading : *leadingDivisors) {
            for (const std::uint64_t constant : *constantDivisors) {
                if (std::gcd(leading, constant) != 1)
                    continue;
                for (const int sign : std::array<int, 2>{-1, 1}) {
                    if (++candidates > options.maximumSparseCandidates)
                        return std::nullopt;
                    std::vector<Rational> coefficients(degree + 1, zero());
                    coefficients[0] = Rational{
                        sign < 0
                            ? -BigInt::fromUnsigned(constant)
                            : BigInt::fromUnsigned(constant)};
                    coefficients[degree] = Rational{
                        BigInt::fromUnsigned(leading)};
                    RationalPolynomial factor = monicRationalPolynomial(
                        RationalPolynomial{std::move(coefficients)});
                    auto quotient = divideRationalPolynomialsExact(
                        polynomial, factor);
                    if (!quotient || quotient->degree() == 0)
                        continue;
                    return std::pair{
                        std::move(factor),
                        monicRationalPolynomial(*quotient)};
                }
            }
        }
    }
    return std::nullopt;
}

struct DeflatedPolynomial final {
    RationalPolynomial polynomial;
    std::size_t exponentGcd = 1;
};

[[nodiscard]] std::optional<DeflatedPolynomial> deflateRationalPolynomial(
    const RationalPolynomial& polynomial) {
    std::size_t exponentGcd = 0;
    for (std::size_t exponent = 1; exponent <= polynomial.degree(); ++exponent) {
        if (polynomial.coefficient(exponent).isZero())
            continue;
        exponentGcd = std::gcd(exponentGcd, exponent);
    }
    if (exponentGcd <= 1)
        return std::nullopt;
    std::vector<Rational> coefficients(
        polynomial.degree() / exponentGcd + 1, zero());
    for (std::size_t exponent = 0; exponent <= polynomial.degree(); ++exponent) {
        if (!polynomial.coefficient(exponent).isZero())
            coefficients[exponent / exponentGcd] = polynomial.coefficient(exponent);
    }
    return DeflatedPolynomial{
        RationalPolynomial{std::move(coefficients)}, exponentGcd};
}

[[nodiscard]] RationalPolynomial inflateRationalPolynomial(
    const RationalPolynomial& polynomial,
    std::size_t exponentMultiplier) {
    std::vector<Rational> coefficients(
        polynomial.degree() * exponentMultiplier + 1, zero());
    for (std::size_t exponent = 0; exponent <= polynomial.degree(); ++exponent)
        coefficients[exponent * exponentMultiplier] = polynomial.coefficient(exponent);
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] RationalPolynomial multiplyByMonicLinear(
    const RationalPolynomial& polynomial,
    std::int64_t root) {
    std::vector<Rational> result(polynomial.degree() + 2, zero());
    const Rational rationalRoot{BigInt{root}};
    for (std::size_t i = 0; i <= polynomial.degree(); ++i) {
        result[i] -= polynomial.coefficient(i) * rationalRoot;
        result[i + 1] += polynomial.coefficient(i);
    }
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial interpolateIntegerSamples(
    const std::vector<std::int64_t>& points,
    const std::vector<BigInt>& values) {
    std::vector<Rational> result(points.size(), zero());
    for (std::size_t i = 0; i < points.size(); ++i) {
        RationalPolynomial basis{{one()}};
        BigInt denominator{1};
        for (std::size_t j = 0; j < points.size(); ++j) {
            if (i == j)
                continue;
            basis = multiplyByMonicLinear(basis, points[j]);
            denominator *= BigInt{points[i] - points[j]};
        }
        const Rational scale{values[i], denominator};
        for (std::size_t exponent = 0; exponent <= basis.degree(); ++exponent)
            result[exponent] += scale * basis.coefficient(exponent);
    }
    return RationalPolynomial{std::move(result)};
}

struct KroneckerSample final {
    std::int64_t point = 0;
    std::vector<std::uint64_t> divisors;
};

struct RationalPolynomialSplit final {
    std::optional<std::pair<RationalPolynomial, RationalPolynomial>> factors;
    bool exhaustive = false;
};

[[nodiscard]] RationalPolynomialSplit kroneckerSplitRationalPolynomial(
    const RationalPolynomial& input,
    const RationalPolynomialFactorOptions& options) {
    const RationalPolynomial polynomial = monicRationalPolynomial(input);
    if (polynomial.degree() < 2
        || polynomial.degree() > options.maximumKroneckerDegree)
        return {};

    // 一次因子はRational Root Theoremで先に除き，補間の組合せを減らす。
    if (const auto root = findRationalRoot(polynomial).root) {
        RationalPolynomial factor{{-*root, one()}};
        if (auto quotient = divideRationalPolynomialsExact(polynomial, factor))
            return RationalPolynomialSplit{
                std::pair{std::move(factor), monicRationalPolynomial(*quotient)}, true};
    }

    const std::vector<BigInt> integerPolynomial = primitiveIntegerCoefficients(polynomial);
    std::vector<KroneckerSample> samples;
    bool sampleEnumerationComplete = true;
    const std::size_t radiusLimit = std::min<std::size_t>(
        options.maximumSampleRadius,
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()));
    for (std::size_t radius = 0; radius <= radiusLimit; ++radius) {
        const auto signedRadius = static_cast<std::int64_t>(radius);
        const std::array<std::int64_t, 2> candidates{signedRadius, -signedRadius};
        for (std::size_t candidateIndex = 0; candidateIndex < candidates.size(); ++candidateIndex) {
            if (radius == 0 && candidateIndex != 0)
                continue;
            const std::int64_t point = candidates[candidateIndex];
            const BigInt value = evaluateIntegerCoefficients(integerPolynomial, point);
            if (value.isZero())
                continue;
            auto divisors = factoredPositiveDivisors(
                value, options.maximumDivisorsPerSample);
            if (!divisors) {
                sampleEnumerationComplete = false;
                continue;
            }
            samples.push_back(KroneckerSample{point, std::move(*divisors)});
        }
    }
    std::sort(samples.begin(), samples.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.divisors.size() != rhs.divisors.size())
            return lhs.divisors.size() < rhs.divisors.size();
        return lhs.point < rhs.point;
    });

    std::size_t totalCombinations = 0;
    bool exhaustive = sampleEnumerationComplete;
    for (std::size_t factorDegree = 1;
         factorDegree <= polynomial.degree() / 2; ++factorDegree) {
        if (samples.size() < factorDegree + 1) {
            exhaustive = false;
            break;
        }
        const std::vector<KroneckerSample> selected(
            samples.begin(),
            samples.begin() + static_cast<std::ptrdiff_t>(factorDegree + 1));
        std::vector<std::int64_t> points;
        points.reserve(selected.size());
        for (const auto& sample : selected)
            points.push_back(sample.point);
        std::vector<BigInt> values(selected.size());
        std::optional<std::pair<RationalPolynomial, RationalPolynomial>> found;

        std::function<void(std::size_t)> search = [&](std::size_t index) {
            if (found || totalCombinations >= options.maximumCombinations)
                return;
            if (index == selected.size()) {
                ++totalCombinations;
                RationalPolynomial candidate = interpolateIntegerSamples(points, values);
                if (candidate.degree() != factorDegree)
                    return;
                for (const Rational& coefficient : candidate.coefficients())
                    if (!coefficient.isInteger())
                        return;
                candidate = monicRationalPolynomial(candidate);
                if (candidate.degree() == 0 || candidate.degree() >= polynomial.degree())
                    return;
                auto quotient = divideRationalPolynomialsExact(polynomial, candidate);
                if (!quotient || quotient->degree() == 0)
                    return;
                found = std::pair{
                    std::move(candidate), monicRationalPolynomial(*quotient)};
                return;
            }

            for (const std::uint64_t divisor : selected[index].divisors) {
                const BigInt positive = BigInt::fromUnsigned(divisor);
                values[index] = positive;
                search(index + 1);
                if (found || totalCombinations >= options.maximumCombinations)
                    return;
                // 全値の符号を同時に反転した候補は同じmonic因子になるため，
                // 最初のsampleだけ正に固定して重複を半減する。
                if (index != 0) {
                    values[index] = -positive;
                    search(index + 1);
                    if (found || totalCombinations >= options.maximumCombinations)
                        return;
                }
            }
        };
        search(0);
        if (found)
            return RationalPolynomialSplit{std::move(found), true};
        if (totalCombinations >= options.maximumCombinations) {
            exhaustive = false;
            break;
        }
    }
    return RationalPolynomialSplit{std::nullopt, exhaustive};
}

void factorRationalPolynomialRecursive(
    const RationalPolynomial& polynomial,
    const RationalPolynomialFactorOptions& options,
    std::vector<RationalPolynomial>& factors,
    bool& complete) {
    if (polynomial.degree() <= 1) {
        factors.push_back(monicRationalPolynomial(polynomial));
        return;
    }
    RationalPolynomialSplit split = kroneckerSplitRationalPolynomial(polynomial, options);
    if (!split.factors) {
        factors.push_back(monicRationalPolynomial(polynomial));
        complete = complete && split.exhaustive;
        return;
    }
    factorRationalPolynomialRecursive(
        split.factors->first, options, factors, complete);
    factorRationalPolynomialRecursive(
        split.factors->second, options, factors, complete);
}

struct AdaptiveFactorizationResult final {
    std::vector<RationalPolynomial> factors;
    bool complete = false;
};

[[nodiscard]] AdaptiveFactorizationResult factorSquareFreePolynomialAdaptive(
    const RationalPolynomial& polynomial,
    const RationalPolynomialFactorOptions& options);

[[nodiscard]] AdaptiveFactorizationResult factorAdaptiveSplit(
    const std::pair<RationalPolynomial, RationalPolynomial>& split,
    const RationalPolynomialFactorOptions& options) {
    AdaptiveFactorizationResult left = factorSquareFreePolynomialAdaptive(
        split.first, options);
    AdaptiveFactorizationResult right = factorSquareFreePolynomialAdaptive(
        split.second, options);
    left.factors.insert(
        left.factors.end(),
        std::make_move_iterator(right.factors.begin()),
        std::make_move_iterator(right.factors.end()));
    left.complete = left.complete && right.complete;
    return left;
}

[[nodiscard]] AdaptiveFactorizationResult factorSquareFreePolynomialAdaptive(
    const RationalPolynomial& input,
    const RationalPolynomialFactorOptions& options) {
    const RationalPolynomial polynomial = monicRationalPolynomial(input);
    if (polynomial.degree() <= 1)
        return AdaptiveFactorizationResult{{polynomial}, true};

    const CertifiedPolynomialSplit binomial = splitRationalBinomial(polynomial);
    if (binomial.factors)
        return factorAdaptiveSplit(*binomial.factors, options);
    if (binomial.irreducible)
        return AdaptiveFactorizationResult{{polynomial}, true};

    if (const auto sparse = boundedSparseBinomialSplit(polynomial, options))
        return factorAdaptiveSplit(*sparse, options);

    // g(x^d)を先にg(y)へ縮める。gの既約性をg(x^d)へは継承せず，
    // gが実際に複数因子へ分かれた場合だけinflateして再帰する。
    if (const auto deflated = deflateRationalPolynomial(polynomial)) {
        const RationalPolynomialFactorization base =
            factorRationalPolynomialOverQ(deflated->polynomial, options);
        if (base.factors.size() > 1 && base.scalar == one()) {
            AdaptiveFactorizationResult inflatedResult{{}, true};
            RationalPolynomial product{{one()}};
            for (const RationalPolynomial& factor : base.factors) {
                RationalPolynomial inflated = inflateRationalPolynomial(
                    factor, deflated->exponentGcd);
                const RationalPolynomialFactorization decomposition =
                    factorRationalPolynomialOverQ(inflated, options);
                if (!(decomposition.scalar == one())
                    || decomposition.factors.empty()) {
                    inflatedResult.factors.clear();
                    break;
                }
                inflatedResult.complete =
                    inflatedResult.complete && decomposition.complete;
                for (const RationalPolynomial& child : decomposition.factors) {
                    product = multiplyRationalPolynomialValues(product, child);
                    inflatedResult.factors.push_back(child);
                }
            }
            if (!inflatedResult.factors.empty()
                && product.coefficients() == polynomial.coefficients())
                return inflatedResult;
        }
    }

    if (auto modular = factorSquareFreeRationalPolynomialModular(
            polynomial, options))
        return AdaptiveFactorizationResult{
            std::move(modular->factors), modular->complete};

    AdaptiveFactorizationResult fallback{{}, true};
    factorRationalPolynomialRecursive(
        polynomial, options, fallback.factors, fallback.complete);
    return fallback;
}

void sortRationalPolynomialFactors(
    std::vector<RationalPolynomial>& factors) {
    std::sort(factors.begin(), factors.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.degree() != rhs.degree())
            return lhs.degree() < rhs.degree();
        const auto& left = lhs.coefficients();
        const auto& right = rhs.coefficients();
        return std::lexicographical_compare(
            left.rbegin(), left.rend(), right.rbegin(), right.rend());
    });
}

} // namespace

RationalPolynomialFactorization factorRationalPolynomialOverQ(
    const RationalPolynomial& polynomial,
    RationalPolynomialFactorOptions options) {
    if (polynomial.isZero())
        return RationalPolynomialFactorization{zero(), {}, true};
    if (polynomial.degree() == 0)
        return RationalPolynomialFactorization{polynomial.coefficient(0), {}, true};

    const Rational scalar = polynomial.coefficient(polynomial.degree());
    const RationalPolynomial normalized = monicRationalPolynomial(polynomial);

    const auto factorCertifiedSplit = [&](
        const std::pair<RationalPolynomial, RationalPolynomial>& split) {
        RationalPolynomialFactorization left = factorRationalPolynomialOverQ(
            split.first, options);
        RationalPolynomialFactorization right = factorRationalPolynomialOverQ(
            split.second, options);
        left.scalar = scalar * left.scalar * right.scalar;
        left.factors.insert(
            left.factors.end(),
            std::make_move_iterator(right.factors.begin()),
            std::make_move_iterator(right.factors.end()));
        left.complete = left.complete && right.complete;
        sortRationalPolynomialFactors(left.factors);
        return left;
    };

    // sparse/exact structureはQ[x]のEuclidean square-free分解より先に試す。
    // 高次数の疎積ではRational係数のgcdが膨張し得る一方，ここで採るsplitは
    // exact division済みなので，先行しても完全性を損なわない。
    const CertifiedPolynomialSplit binomial = splitRationalBinomial(normalized);
    if (binomial.factors)
        return factorCertifiedSplit(*binomial.factors);
    if (binomial.irreducible)
        return RationalPolynomialFactorization{
            scalar, {normalized}, true};
    if (const auto sparse = boundedSparseBinomialSplit(normalized, options))
        return factorCertifiedSplit(*sparse);

    if (const auto deflated = deflateRationalPolynomial(normalized)) {
        const RationalPolynomialFactorization base =
            factorRationalPolynomialOverQ(deflated->polynomial, options);
        if (base.factors.size() > 1 && base.scalar == one()) {
            RationalPolynomialFactorization inflated{scalar, {}, true};
            for (const RationalPolynomial& factor : base.factors) {
                RationalPolynomialFactorization child = factorRationalPolynomialOverQ(
                    inflateRationalPolynomial(factor, deflated->exponentGcd), options);
                inflated.scalar *= child.scalar;
                inflated.complete = inflated.complete && child.complete;
                inflated.factors.insert(
                    inflated.factors.end(),
                    std::make_move_iterator(child.factors.begin()),
                    std::make_move_iterator(child.factors.end()));
            }
            if (!inflated.factors.empty()) {
                sortRationalPolynomialFactors(inflated.factors);
                return inflated;
            }
        }
    }

    std::vector<RationalPolynomial> factors;
    bool complete = true;
    if (hasModularSquareFreeCertificate(normalized)) {
        AdaptiveFactorizationResult adaptive =
            factorSquareFreePolynomialAdaptive(normalized, options);
        factors = std::move(adaptive.factors);
        complete = adaptive.complete;
    }
    else {
        const auto squareFree = squareFreeRationalFactorization(normalized);
        if (!squareFree) {
            factorRationalPolynomialRecursive(
                normalized, options, factors, complete);
        }
        else {
            for (const SquareFreeRationalFactor& component : *squareFree) {
                AdaptiveFactorizationResult adaptive =
                    factorSquareFreePolynomialAdaptive(
                        component.polynomial, options);
                complete = complete && adaptive.complete;
                for (std::size_t multiplicity = 0;
                     multiplicity < component.multiplicity; ++multiplicity)
                    factors.insert(
                        factors.end(), adaptive.factors.begin(), adaptive.factors.end());
            }
        }
    }
    sortRationalPolynomialFactors(factors);
    return RationalPolynomialFactorization{scalar, std::move(factors), complete};
}

std::vector<expression::Symbol> collectSymbols(const Expr& expression) {
    std::vector<expression::Symbol> result;
    std::vector<Expr> stack{expression};
    while (!stack.empty()) {
        Expr current = std::move(stack.back());
        stack.pop_back();
        if (current.isSymbol()) {
            const auto symbol = current.asSymbol();
            if (std::find(result.begin(), result.end(), symbol) == result.end())
                result.push_back(symbol);
            continue;
        }
        if (current.isArray()) {
            if (current.asArray().storageKind() == expression::ArrayStorageKind::Generic)
                for (const Expr& element : current.asArray().storedExpressions())
                    stack.push_back(element);
            continue;
        }
        if (current.isCall())
            for (const Expr& argument : current.asCall().arguments)
                stack.push_back(argument);
    }
    std::sort(result.begin(), result.end(), symbolLess);
    return result;
}

} // namespace mmcal::symbolic
