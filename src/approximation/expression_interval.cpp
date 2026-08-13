// Exprとcertified interval backendの共通変換
#include "expression_interval.hpp"

#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/real_number.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcal::approximation {
namespace {

[[nodiscard]] RealInterval intervalFromDecimal(
    const numeric::DecimalApproximation& value,
    std::size_t precisionBits) {
    return RealInterval::fromRationalBounds(
        value.certifiedLower(), value.certifiedUpper(), precisionBits);
}

[[nodiscard]] bool exactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

void collectApproximationDigits(
    const expression::Expr& root,
    std::optional<std::size_t>& digits) {
    std::vector<const expression::Expr*> pending{&root};
    while (!pending.empty()) {
        const expression::Expr& current = *pending.back();
        pending.pop_back();

        std::optional<std::size_t> currentDigits;
        if (current.isDecimalApproximation())
            currentDigits = current.asDecimalApproximation().requestedFractionalDigits();
        else if (current.isComplexDecimalApproximation()) {
            const auto& complex = current.asComplexDecimalApproximation();
            currentDigits = std::min(
                complex.real().requestedFractionalDigits(),
                complex.imaginary().requestedFractionalDigits());
        }

        if (currentDigits && *currentDigits != 0)
            digits = digits ? std::min(*digits, *currentDigits) : currentDigits;

        if (current.isArray())
            for (const expression::Expr& element : current.asArray().elements)
                pending.push_back(&element);
        else if (current.isList())
            for (const expression::Expr& element : current.asList().elements)
                pending.push_back(&element);
        else if (current.isCall())
            for (const expression::Expr& argument : current.asCall().arguments)
                pending.push_back(&argument);
    }
}

} // namespace

std::optional<ComplexInterval> encloseComplexExpression(
    const expression::Expr& expression,
    std::size_t precisionBits,
    const CertifiedEvaluator& certified) {
    if (expression.isDecimalApproximation())
        return ComplexInterval::fromReal(
            intervalFromDecimal(expression.asDecimalApproximation(), precisionBits));

    if (expression.isComplexDecimalApproximation()) {
        const auto& value = expression.asComplexDecimalApproximation();
        return ComplexInterval{
            intervalFromDecimal(value.real(), precisionBits),
            intervalFromDecimal(value.imaginary(), precisionBits)};
    }

    const auto enclosed = certified.enclose(expression, precisionBits);
    if (!enclosed)
        return std::nullopt;
    return enclosed->toComplex();
}

std::optional<expression::Expr> decimalExpression(
    const RealInterval& value,
    std::size_t fractionalDigits) {
    if (value.isPoint())
        return expression::Expr{numeric::DecimalApproximation::fromReal(
            numeric::RealNumber{value.lower().toRational()}, fractionalDigits)};

    const auto decimal = numeric::DecimalApproximation::fromCertifiedInterval(
        value.lower().toRational(), value.upper().toRational(), fractionalDigits);
    return decimal ? std::optional<expression::Expr>{expression::Expr{*decimal}}
                   : std::nullopt;
}

std::optional<expression::Expr> decimalExpression(
    const ComplexInterval& value,
    std::size_t fractionalDigits) {
    if (value.real().isPoint() && value.imaginary().isPoint()) {
        const auto real = numeric::DecimalApproximation::fromReal(
            numeric::RealNumber{value.real().lower().toRational()}, fractionalDigits);
        const auto imaginary = numeric::DecimalApproximation::fromReal(
            numeric::RealNumber{value.imaginary().lower().toRational()}, fractionalDigits);
        if (exactZero(value.imaginary()))
            return expression::Expr{real};
        return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
            real, imaginary, exactZero(value.real()), false)};
    }

    const auto real = numeric::DecimalApproximation::fromCertifiedInterval(
        value.real().lower().toRational(), value.real().upper().toRational(), fractionalDigits);
    const auto imaginary = numeric::DecimalApproximation::fromCertifiedInterval(
        value.imaginary().lower().toRational(), value.imaginary().upper().toRational(), fractionalDigits);
    if (!real || !imaginary)
        return std::nullopt;
    if (exactZero(value.imaginary()))
        return expression::Expr{*real};
    return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
        *real, *imaginary, exactZero(value.real()), false)};
}

std::optional<ApproximationContext> inferredApproximationContext(
    std::span<const expression::Expr> expressions) {
    std::optional<std::size_t> digits;
    for (const expression::Expr& expression : expressions)
        collectApproximationDigits(expression, digits);
    if (!digits)
        return std::nullopt;
    return ApproximationContext{*digits};
}

std::size_t nextGuardDigits(std::size_t current) {
    const std::size_t growth = std::max<std::size_t>(8, current / 2);
    if (growth > std::numeric_limits<std::size_t>::max() - current)
        throw std::overflow_error("Approximation precision is too large");
    return current + growth;
}

} // namespace mmcal::approximation
