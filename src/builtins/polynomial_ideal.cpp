// Gröbner basis / multivariate polynomial reductionの公開builtin境界
#include "polynomial_ideal.hpp"

#include "error/error_message.hpp"
#include "expression/array_utils.hpp"
#include "symbolic/groebner.hpp"
#include "symbolic/polynomial.hpp"

#include <optional>
#include <string_view>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using expression::Expr;

[[nodiscard]] std::optional<std::vector<Expr>> braceElements(const Expr& expression) {
    if (expression.isArray()) {
        const auto& array = expression.asArray();
        if (array.rank() != 1)
            return std::nullopt;
        return array.materialize();
    }
    if (expression.isList())
        return expression.asList().elements;
    return std::nullopt;
}

[[nodiscard]] std::vector<expression::Symbol> parseVariables(const Expr& expression) {
    const auto values = braceElements(expression);
    if (!values || values->empty())
        error::throwCalcError(
            error::CalcErrorType::Type,
            "polynomial variables must be a non-empty brace value of symbols");
    std::vector<expression::Symbol> variables;
    variables.reserve(values->size());
    for (const Expr& value : *values) {
        if (!value.isSymbol())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "polynomial variables must be symbols");
        variables.push_back(value.asSymbol());
    }
    return variables;
}

[[nodiscard]] symbolic::MonomialOrder parseOrder(
    std::span<const Expr> arguments,
    std::size_t index,
    symbolic::MonomialOrder defaultOrder) {
    if (arguments.size() <= index)
        return defaultOrder;
    std::string_view name;
    if (arguments[index].isSymbol())
        name = arguments[index].asSymbol().view();
    else if (arguments[index].isString())
        name = arguments[index].asString();
    else
        error::throwCalcError(
            error::CalcErrorType::Type,
            "monomial order must be Lex, GrLex, or GrevLex");
    const auto order = symbolic::parseMonomialOrder(name);
    if (!order)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "monomial order must be Lex, GrLex, or GrevLex");
    return *order;
}

[[nodiscard]] std::vector<symbolic::MultivariateRationalPolynomial> parsePolynomials(
    const Expr& expression,
    const symbolic::PolynomialRing& ring,
    const evaluation::BuiltinRegistry& builtins) {
    const auto values = braceElements(expression);
    if (!values)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "polynomials must be supplied as a brace value");
    std::vector<symbolic::MultivariateRationalPolynomial> polynomials;
    polynomials.reserve(values->size());
    for (const Expr& value : *values) {
        const auto polynomial = symbolic::toMultivariateRationalPolynomial(
            value, ring, builtins,
            symbolic::PolynomialConversionOptions{4096, 100'000});
        if (!polynomial)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Groebner arithmetic requires exact Rational-coefficient polynomials in the declared variables");
        polynomials.push_back(*polynomial);
    }
    return polynomials;
}

[[nodiscard]] Expr polynomialList(
    std::span<const symbolic::MultivariateRationalPolynomial> polynomials,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Expr> values;
    values.reserve(polynomials.size());
    for (const auto& polynomial : polynomials)
        values.push_back(symbolic::polynomialToExpandedExpr(polynomial, builtins));
    return expression::braceValue(std::move(values));
}

} // namespace

Expr evaluateGroebnerBasis(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins) {
    const auto variables = parseVariables(arguments[1]);
    const symbolic::MonomialOrder order = parseOrder(
        arguments, 2, symbolic::MonomialOrder::GrevLex);
    const symbolic::PolynomialRing ring{variables, order};

    // HoldAllは多項式変数を既存のsession定義から保護するために必要である。
    // ただしgroebnerBasis[groebnerBasis[...], ...]の直接合成まで拒否する理由はないため，
    // 同じsymbolic境界内の明示的なnested Gröbner callだけは安全に先に解決する。
    Expr generatorExpression = arguments[0];
    if (builtins.isCallTo(generatorExpression, evaluation::BuiltinId::GroebnerBasis))
        generatorExpression = evaluateGroebnerBasis(
            generatorExpression.asCall().arguments, builtins);
    const auto generators = parsePolynomials(generatorExpression, ring, builtins);
    const symbolic::GroebnerComputation result = symbolic::groebnerBasis(
        generators, ring);
    return polynomialList(result.basis, builtins);
}

Expr evaluatePolynomialReduce(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins) {
    const auto variables = parseVariables(arguments[2]);
    const symbolic::MonomialOrder order = parseOrder(
        arguments, 3, symbolic::MonomialOrder::GrevLex);
    const symbolic::PolynomialRing ring{variables, order};

    Expr divisorExpression = arguments[1];
    if (builtins.isCallTo(divisorExpression, evaluation::BuiltinId::GroebnerBasis))
        divisorExpression = evaluateGroebnerBasis(
            divisorExpression.asCall().arguments, builtins);
    const auto divisors = parsePolynomials(divisorExpression, ring, builtins);
    const auto dividend = symbolic::toMultivariateRationalPolynomial(
        arguments[0], ring, builtins,
        symbolic::PolynomialConversionOptions{4096, 100'000});
    if (!dividend)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "polynomialReduce requires an exact Rational-coefficient polynomial in the declared variables");

    const symbolic::PolynomialDivisionResult result = symbolic::multivariateDivide(
        *dividend, divisors, ring);
    return expression::braceValue({
        polynomialList(result.quotients, builtins),
        symbolic::polynomialToExpandedExpr(result.remainder, builtins)});
}

} // namespace mmcal::builtins
