#include "risch_expression.hpp"

#include "expression/exact_value.hpp"
#include "numeric/big_int.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "symbolic/algebraic_expression.hpp"
#include "symbolic/algebraic_number.hpp"
#include "symbolic/polynomial.hpp"

#include <utility>
#include <vector>

namespace mmcal::symbolic::risch {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using expression::exact::integer;
using expression::exact::rational;

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

[[nodiscard]] Expr add(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplify(
        call(builtins, BuiltinId::Add, {std::move(lhs), std::move(rhs)}),
        builtins, mathematics, angles);
}

[[nodiscard]] Expr multiply(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplify(
        call(builtins, BuiltinId::Multiply, {std::move(lhs), std::move(rhs)}),
        builtins, mathematics, angles);
}

[[nodiscard]] Expr evaluateResiduePolynomial(
    const RationalPolynomial& polynomial,
    const Expr& root,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr value = rational(polynomial.coefficient(polynomial.degree()));
    for (std::size_t exponent = polynomial.degree(); exponent != 0; --exponent)
        value = add(
            multiply(value, root, builtins, mathematics, angles),
            rational(polynomial.coefficient(exponent - 1)),
            builtins, mathematics, angles);
    return value;
}

[[nodiscard]] Expr evaluateBivariateAtResidue(
    const BivariateRationalPolynomial& polynomial,
    const Expr& residue,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr value = evaluateResiduePolynomial(
        polynomial.coefficientInX(polynomial.degreeInX()), residue,
        builtins, mathematics, angles);
    const Expr variableExpression{variable};
    for (std::size_t exponent = polynomial.degreeInX(); exponent != 0; --exponent)
        value = add(
            multiply(value, variableExpression, builtins, mathematics, angles),
            evaluateResiduePolynomial(
                polynomial.coefficientInX(exponent - 1), residue,
                builtins, mathematics, angles),
            builtins, mathematics, angles);
    return value;
}

[[nodiscard]] Expr rationalFunctionExpression(
    const RationalFunction& value,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (value.numerator.isZero())
        return integer(0);
    Expr numerator = polynomialToExpandedExpr(value.numerator, variable, builtins);
    if (value.denominator.degree() == 0) {
        const numeric::Rational scale = value.denominator.coefficient(0);
        if (scale == numeric::Rational{numeric::BigInt{1}})
            return numerator;
    }
    Expr denominator = polynomialToExpandedExpr(value.denominator, variable, builtins);
    return simplify(
        call(builtins, BuiltinId::Divide, {
            std::move(numerator), std::move(denominator)}),
        builtins, mathematics, angles);
}

} // namespace

RischStageResult<expression::Expr> materializeLrtLogarithms(
    const RationalFunction& properSquareFreePart,
    const LrtResult& result,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const RischOptions& options) {
    if (!result.exactVerified
        || !verifyLrtResult(properSquareFreePart, result, options))
        return {std::nullopt, RischFailure::CertificateFailed};

    std::vector<Expr> logarithms;
    for (const AlgebraicResidueLogTerm& term : result.logarithmicTerms) {
        if (!term.exactVerified || term.residuePolynomial.degree() == 0
            || term.residuePolynomial.degree() > options.maximumResidueDegree)
            return {std::nullopt, RischFailure::CertificateFailed};
        auto roots = ComplexAlgebraicNumber::isolateAll(
            term.residuePolynomial.coefficients());
        if (!roots || roots->size() != term.residuePolynomial.degree())
            return {std::nullopt, RischFailure::UnsupportedExtension};
        if (auto canonical = ComplexAlgebraicNumber::canonicalizeAll(*roots))
            roots = std::move(canonical);
        logarithms.reserve(logarithms.size() + roots->size());
        for (const ComplexAlgebraicNumber& root : *roots) {
            Expr residue = makeCanonicalRootExpression(root, builtins);
            Expr argument = evaluateBivariateAtResidue(
                term.logArgument, residue, variable,
                builtins, mathematics, angles);
            Expr logarithm = call(
                builtins, BuiltinId::Log, {std::move(argument)});
            logarithms.push_back(multiply(
                residue, std::move(logarithm),
                builtins, mathematics, angles));
        }
    }
    if (logarithms.empty())
        return {integer(0), RischFailure::None};
    Expr sum = std::move(logarithms.front());
    for (std::size_t i = 1; i < logarithms.size(); ++i)
        sum = add(
            std::move(sum), std::move(logarithms[i]),
            builtins, mathematics, angles);
    return {std::move(sum), RischFailure::None};
}

RischStageResult<expression::Expr> materializeRationalRischResult(
    const RischResult& result,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const RischOptions& options) {
    if (result.status != RischResultStatus::Elementary
        || !result.exactVerified || !result.rational)
        return {std::nullopt, RischFailure::CertificateFailed};
    auto logarithmic = materializeLrtLogarithms(
        result.rational->hermite.squareFreePart,
        result.rational->logarithmicPart, variable,
        builtins, mathematics, angles, options);
    if (!logarithmic)
        return logarithmic;
    Expr rationalPart = rationalFunctionExpression(
        result.rational->hermite.rationalPart, variable,
        builtins, mathematics, angles);
    return {
        add(std::move(rationalPart), std::move(*logarithmic.value),
            builtins, mathematics, angles),
        RischFailure::None};
}

} // namespace mmcal::symbolic::risch
