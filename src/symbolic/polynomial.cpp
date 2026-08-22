// 多項式表現と操作
#include "polynomial.hpp"

#include "mathematics/exact_algebra.hpp"
#include "mathematics/definedness.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <limits>
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
