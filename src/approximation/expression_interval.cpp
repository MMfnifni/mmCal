// Exprとcertified interval backendの共通変換
#include "expression_interval.hpp"

#include "certification_error.hpp"

#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/real_number.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
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

[[nodiscard]] RealInterval informationIntervalFromDecimal(
    const numeric::DecimalApproximation& value,
    std::size_t precisionBits) {
    return RealInterval::fromRationalBounds(
        value.informationLower(), value.informationUpper(), precisionBits);
}

[[nodiscard]] bool exactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

[[nodiscard]] std::optional<CertifiedValue> storedNumericInterval(
    const expression::Expr& value,
    std::size_t precisionBits,
    bool information) {
    if (value.isNumber()) {
        const auto& number = value.asNumber();
        if (number.isReal())
            return CertifiedValue{RealInterval::fromRational(
                number.asReal().toRational(), precisionBits)};

        const auto& complex = number.asComplex();
        return CertifiedValue{ComplexInterval{
            RealInterval::fromRational(complex.real.toRational(), precisionBits),
            RealInterval::fromRational(complex.imaginary.toRational(), precisionBits)}};
    }

    if (value.isDecimalApproximation())
        return CertifiedValue{information
            ? informationIntervalFromDecimal(value.asDecimalApproximation(), precisionBits)
            : intervalFromDecimal(value.asDecimalApproximation(), precisionBits)};

    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        return CertifiedValue{ComplexInterval{
            information ? informationIntervalFromDecimal(complex.real(), precisionBits)
                        : intervalFromDecimal(complex.real(), precisionBits),
            information ? informationIntervalFromDecimal(complex.imaginary(), precisionBits)
                        : intervalFromDecimal(complex.imaginary(), precisionBits)}};
    }

    return std::nullopt;
}

[[nodiscard]] CertifiedValue normalizeComplex(ComplexInterval value) {
    if (value.isProvablyReal())
        return CertifiedValue{value.real()};
    return CertifiedValue{std::move(value)};
}

[[nodiscard]] CertifiedValue addValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{add(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeComplex(add(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue subtractValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{subtract(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeComplex(subtract(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue multiplyValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{multiply(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeComplex(multiply(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue divideValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{divide(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeComplex(divide(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue negateValue(const CertifiedValue& value) {
    if (value.isReal())
        return CertifiedValue{negate(value.asReal())};
    return normalizeComplex(negate(value.asComplex()));
}

[[nodiscard]] std::size_t guaranteedDigits(
    const numeric::Rational& error,
    std::size_t cap);

[[nodiscard]] std::size_t intervalDigitHint(
    const RealInterval& value,
    std::size_t maximumDigits) {
    if (value.isPoint())
        return maximumDigits;

    const numeric::Rational lower = value.lower().toRational();
    const numeric::Rational upper = value.upper().toRational();
    if (lower <= numeric::Rational{} && upper >= numeric::Rational{})
        return 1;
    const numeric::Rational minimumMagnitude = lower > numeric::Rational{} ? lower : -upper;
    if (minimumMagnitude.isZero())
        return 1;
    const numeric::Rational width = upper - lower;
    return std::max<std::size_t>(1, guaranteedDigits(width / minimumMagnitude, maximumDigits));
}

[[nodiscard]] std::size_t intervalDigitHint(
    const CertifiedValue& value,
    std::size_t maximumFractionalDigits) {
    if (value.isReal())
        return intervalDigitHint(value.asReal(), maximumFractionalDigits);

    return std::min(
        intervalDigitHint(value.asComplex().real(), maximumFractionalDigits),
        intervalDigitHint(value.asComplex().imaginary(), maximumFractionalDigits));
}

[[nodiscard]] std::optional<expression::Expr> decimalAtDigits(
    const CertifiedValue& value,
    std::size_t significantDigits) {
    return value.isReal()
        ? decimalExpression(value.asReal(), significantDigits)
        : decimalExpression(value.asComplex(), significantDigits);
}

[[nodiscard]] std::optional<expression::Expr> bestDecimalExpression(
    const CertifiedValue& value,
    std::size_t maximumDigits) {
    if (maximumDigits == 0)
        return std::nullopt;
    if (const auto decimal = decimalAtDigits(value, maximumDigits))
        return decimal;

    // certified幅に対する相対桁数から候補へ跳ぶ。巨大precisionでも1桁ずつ降りない。
    const std::size_t hint = intervalDigitHint(value, maximumDigits);
    std::size_t digits = std::min(maximumDigits - 1, hint + 2);

    for (; digits > 1; --digits)
        if (const auto decimal = decimalAtDigits(value, digits))
            return decimal;
    return decimalAtDigits(value, 1);
}

[[nodiscard]] numeric::Rational absRational(const numeric::Rational& value) {
    return value.numerator().isNegative() ? -value : value;
}

[[nodiscard]] numeric::Rational maximum(
    const numeric::Rational& lhs,
    const numeric::Rational& rhs) {
    return lhs < rhs ? rhs : lhs;
}

[[nodiscard]] std::size_t guaranteedDigits(
    const numeric::Rational& error,
    std::size_t cap) {
    if (error.isZero())
        return cap;
    if (cap == 0 || error >= numeric::Rational{numeric::BigInt{1}})
        return 0;

    const std::size_t numeratorDigits = error.numerator().abs().toString().size();
    const std::size_t denominatorDigits = error.denominator().toString().size();
    if (denominatorDigits <= numeratorDigits)
        return 0;

    std::size_t candidate = std::min(cap, denominatorDigits - numeratorDigits);
    const auto exponent = static_cast<std::uint64_t>(candidate);
    const numeric::BigInt scaledNumerator = error.numerator().abs()
        * numeric::pow(numeric::BigInt{10}, exponent);
    if (scaledNumerator < error.denominator())
        return candidate;
    return candidate == 0 ? 0 : candidate - 1;
}

[[nodiscard]] std::size_t componentPrecisionCap(
    const numeric::DecimalApproximation& displayed,
    const RealInterval& information,
    std::size_t cap) {
    const numeric::Rational lowerError = absRational(
        displayed.displayedValue() - information.lower().toRational());
    const numeric::Rational upperError = absRational(
        displayed.displayedValue() - information.upper().toRational());
    const numeric::Rational error = maximum(lowerError, upperError);
    if (error.isZero())
        return cap;

    const numeric::Rational lower = information.lower().toRational();
    const numeric::Rational upper = information.upper().toRational();
    if (lower <= numeric::Rational{} && upper >= numeric::Rational{})
        return 0;

    const numeric::Rational minimumMagnitude = lower > numeric::Rational{}
        ? lower
        : -upper;
    if (minimumMagnitude.isZero())
        return 0;
    return guaranteedDigits(error / minimumMagnitude, cap);
}

[[nodiscard]] std::size_t informationPrecisionCap(
    const expression::Expr& displayed,
    const CertifiedValue& information,
    std::size_t cap) {
    if (displayed.isDecimalApproximation()) {
        const RealInterval& interval = information.isReal()
            ? information.asReal() : information.asComplex().real();
        return componentPrecisionCap(displayed.asDecimalApproximation(), interval, cap);
    }

    if (!displayed.isComplexDecimalApproximation())
        return cap;
    const auto& complex = displayed.asComplexDecimalApproximation();
    const ComplexInterval interval = information.toComplex();
    return std::min(
        componentPrecisionCap(complex.real(), interval.real(), cap),
        componentPrecisionCap(complex.imaginary(), interval.imaginary(), cap));
}

[[nodiscard]] std::optional<expression::Expr> decimalExpressionWithInformation(
    const RealInterval& certified,
    const RealInterval& information,
    std::size_t significantDigits) {
    const auto decimal = numeric::DecimalApproximation::fromCertifiedIntervalWithInformationSignificant(
        certified.lower().toRational(),
        certified.upper().toRational(),
        information.lower().toRational(),
        information.upper().toRational(),
        significantDigits);
    return decimal ? std::optional<expression::Expr>{expression::Expr{*decimal}}
                   : std::nullopt;
}

[[nodiscard]] std::optional<expression::Expr> decimalExpressionWithInformation(
    const ComplexInterval& certified,
    const ComplexInterval& information,
    std::size_t significantDigits) {
    const auto real = numeric::DecimalApproximation::fromCertifiedIntervalWithInformationSignificant(
        certified.real().lower().toRational(),
        certified.real().upper().toRational(),
        information.real().lower().toRational(),
        information.real().upper().toRational(),
        significantDigits);
    const auto imaginary = numeric::DecimalApproximation::fromCertifiedIntervalWithInformationSignificant(
        certified.imaginary().lower().toRational(),
        certified.imaginary().upper().toRational(),
        information.imaginary().lower().toRational(),
        information.imaginary().upper().toRational(),
        significantDigits);
    if (!real || !imaginary)
        return std::nullopt;
    if (exactZero(certified.imaginary()))
        return expression::Expr{*real};
    return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
        *real, *imaginary, exactZero(certified.real()), false)};
}

[[nodiscard]] std::optional<expression::Expr> decimalExpressionWithInformation(
    const CertifiedValue& certified,
    const CertifiedValue& information,
    std::size_t significantDigits) {
    if (certified.isReal()) {
        const RealInterval& info = information.isReal()
            ? information.asReal() : information.asComplex().real();
        return decimalExpressionWithInformation(certified.asReal(), info, significantDigits);
    }
    return decimalExpressionWithInformation(
        certified.asComplex(), information.toComplex(), significantDigits);
}

[[nodiscard]] std::optional<expression::Expr> finalizeApproximateOperation(
    const CertifiedValue& certified,
    const CertifiedValue& information,
    std::size_t sourceDigits) {
    const auto preliminary = bestDecimalExpression(certified, sourceDigits);
    if (!preliminary)
        return std::nullopt;
    const std::size_t cap = informationPrecisionCap(*preliminary, information, sourceDigits);
    if (cap == 0) {
        const bool crossesZero = information.isReal()
            ? information.asReal().containsZero()
            : information.asComplex().containsZero();
        if (crossesZero)
            return decimalExpressionWithInformation(certified, information, sourceDigits);
    }
    return decimalExpressionWithInformation(certified, information, std::max<std::size_t>(cap, 1));
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
            currentDigits = current.asDecimalApproximation().requestedSignificantDigits();
        else if (current.isComplexDecimalApproximation()) {
            const auto& complex = current.asComplexDecimalApproximation();
            currentDigits = std::min(
                complex.real().requestedSignificantDigits(),
                complex.imaginary().requestedSignificantDigits());
        }

        if (currentDigits && *currentDigits != 0)
            digits = digits ? std::min(*digits, *currentDigits) : currentDigits;

        if (current.isArray()) {
            const auto& array = current.asArray();
            for (std::size_t i = 0; i < array.size(); ++i) {
                switch (array.storedKindAt(i)) {
                case expression::ArrayStorageKind::DecimalApproximation: {
                    const auto valueDigits = array.decimalAt(i).requestedSignificantDigits();
                    if (valueDigits != 0)
                        digits = digits ? std::min(*digits, valueDigits) : valueDigits;
                    break;
                }
                case expression::ArrayStorageKind::ComplexDecimalApproximation: {
                    const auto& value = array.complexDecimalAt(i);
                    const auto valueDigits = std::min(
                        value.real().requestedSignificantDigits(),
                        value.imaginary().requestedSignificantDigits());
                    if (valueDigits != 0)
                        digits = digits ? std::min(*digits, valueDigits) : valueDigits;
                    break;
                }
                case expression::ArrayStorageKind::Generic:
                    pending.push_back(&array.expressionAt(i));
                    break;
                default:
                    break;
                }
            }
        }
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
    std::size_t significantDigits) {
    if (value.isPoint()) {
        const numeric::RealNumber exact{value.lower().toRational()};
        return expression::Expr{significantDigits == 0
            ? numeric::DecimalApproximation::fromRealFixed(exact, 0)
            : numeric::DecimalApproximation::fromRealSignificant(exact, significantDigits)};
    }

    const auto decimal = numeric::DecimalApproximation::fromCertifiedIntervalSignificant(
        value.lower().toRational(), value.upper().toRational(), significantDigits);
    return decimal ? std::optional<expression::Expr>{expression::Expr{*decimal}}
                   : std::nullopt;
}

std::optional<expression::Expr> decimalExpression(
    const ComplexInterval& value,
    std::size_t significantDigits) {
    if (value.real().isPoint() && value.imaginary().isPoint()) {
        const numeric::RealNumber exactReal{value.real().lower().toRational()};
        const numeric::RealNumber exactImaginary{value.imaginary().lower().toRational()};
        const auto real = significantDigits == 0
            ? numeric::DecimalApproximation::fromRealFixed(exactReal, 0)
            : numeric::DecimalApproximation::fromRealSignificant(exactReal, significantDigits);
        const auto imaginary = significantDigits == 0
            ? numeric::DecimalApproximation::fromRealFixed(exactImaginary, 0)
            : numeric::DecimalApproximation::fromRealSignificant(exactImaginary, significantDigits);
        if (exactZero(value.imaginary()))
            return expression::Expr{real};
        return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
            real, imaginary, exactZero(value.real()), false)};
    }

    const auto real = numeric::DecimalApproximation::fromCertifiedIntervalSignificant(
        value.real().lower().toRational(), value.real().upper().toRational(), significantDigits);
    const auto imaginary = numeric::DecimalApproximation::fromCertifiedIntervalSignificant(
        value.imaginary().lower().toRational(), value.imaginary().upper().toRational(), significantDigits);
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

std::optional<expression::Expr> evaluateApproximateExpression(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::array<expression::Expr, 1> expressions{expression};
    auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    CertifiedEvaluator evaluator{builtins, mathematics, angles};
    for (std::size_t attempt = 0; attempt < 12; ++attempt) {
        try {
            const std::size_t bits = context->workingBinaryBits();
            const auto certified = evaluator.enclose(
                expression, bits, CertifiedEvaluator::EnclosureKind::Certified);
            const auto information = evaluator.enclose(
                expression, bits, CertifiedEvaluator::EnclosureKind::Information);
            if (!certified || !information)
                return std::nullopt;
            if (const auto result = finalizeApproximateOperation(
                    *certified, *information, context->decimalDigits()))
                return result;
        }
        catch (const PrecisionInsufficient&) {
            // backendのguard不足なら再試行する。InformationEnclosureそのものが
            // branch/domain境界を跨ぐ場合は上限回数で保守的に未評価へ戻る。
        }
        catch (const std::domain_error&) {
            return std::nullopt;
        }
        context->setGuardDigits(nextGuardDigits(context->guardDigits()));
    }
    return std::nullopt;
}

std::optional<expression::Expr> addApproximateScalars(
    std::span<const expression::Expr> expressions) {
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    CertifiedValue result{RealInterval::fromRational(numeric::Rational{}, precisionBits)};
    CertifiedValue information{RealInterval::fromRational(numeric::Rational{}, precisionBits)};
    for (const auto& expression : expressions) {
        const auto enclosed = storedNumericInterval(expression, precisionBits, false);
        const auto informationValue = storedNumericInterval(expression, precisionBits, true);
        if (!enclosed || !informationValue)
            return std::nullopt;
        result = addValues(result, *enclosed, precisionBits);
        information = addValues(information, *informationValue, precisionBits);
    }
    return finalizeApproximateOperation(result, information, context->decimalDigits());
}

std::optional<expression::Expr> subtractApproximateScalars(
    const expression::Expr& lhs,
    const expression::Expr& rhs) {
    const std::array<expression::Expr, 2> expressions{lhs, rhs};
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    const auto left = storedNumericInterval(lhs, precisionBits, false);
    const auto right = storedNumericInterval(rhs, precisionBits, false);
    const auto informationLeft = storedNumericInterval(lhs, precisionBits, true);
    const auto informationRight = storedNumericInterval(rhs, precisionBits, true);
    if (!left || !right || !informationLeft || !informationRight)
        return std::nullopt;
    return finalizeApproximateOperation(
        subtractValues(*left, *right, precisionBits),
        subtractValues(*informationLeft, *informationRight, precisionBits),
        context->decimalDigits());
}

std::optional<expression::Expr> multiplyApproximateScalars(
    std::span<const expression::Expr> expressions) {
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    CertifiedValue result{RealInterval::fromRational(
        numeric::Rational{numeric::BigInt{1}}, precisionBits)};
    CertifiedValue information{RealInterval::fromRational(
        numeric::Rational{numeric::BigInt{1}}, precisionBits)};
    for (const auto& expression : expressions) {
        const auto enclosed = storedNumericInterval(expression, precisionBits, false);
        const auto informationValue = storedNumericInterval(expression, precisionBits, true);
        if (!enclosed || !informationValue)
            return std::nullopt;
        result = multiplyValues(result, *enclosed, precisionBits);
        information = multiplyValues(information, *informationValue, precisionBits);
    }
    return finalizeApproximateOperation(result, information, context->decimalDigits());
}

std::optional<expression::Expr> divideApproximateScalars(
    const expression::Expr& lhs,
    const expression::Expr& rhs) {
    const std::array<expression::Expr, 2> expressions{lhs, rhs};
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    const auto left = storedNumericInterval(lhs, precisionBits, false);
    const auto right = storedNumericInterval(rhs, precisionBits, false);
    const auto informationLeft = storedNumericInterval(lhs, precisionBits, true);
    const auto informationRight = storedNumericInterval(rhs, precisionBits, true);
    if (!left || !right || !informationLeft || !informationRight)
        return std::nullopt;

    try {
        return finalizeApproximateOperation(
            divideValues(*left, *right, precisionBits),
            divideValues(*informationLeft, *informationRight, precisionBits),
            context->decimalDigits());
    }
    catch (const std::domain_error&) {
        return std::nullopt;
    }
}

std::optional<expression::Expr> negateApproximateScalar(
    const expression::Expr& value) {
    const std::array<expression::Expr, 1> expressions{value};
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    const auto enclosed = storedNumericInterval(value, precisionBits, false);
    const auto informationValue = storedNumericInterval(value, precisionBits, true);
    if (!enclosed || !informationValue)
        return std::nullopt;
    return finalizeApproximateOperation(
        negateValue(*enclosed), negateValue(*informationValue), context->decimalDigits());
}

std::size_t nextGuardDigits(std::size_t current) {
    const std::size_t growth = std::max<std::size_t>(8, current / 2);
    if (growth > std::numeric_limits<std::size_t>::max() - current)
        throw std::overflow_error("Approximation precision is too large");
    return current + growth;
}

} // namespace mmcal::approximation
