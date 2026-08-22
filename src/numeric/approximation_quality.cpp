// DecimalApproximation / ComplexDecimalApproximationの品質計算
#include "approximation_quality.hpp"

#include "integer_algorithms.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>

namespace mmcal::numeric {
namespace {

[[nodiscard]] Rational absRational(const Rational& value) {
    return value.numerator().isNegative() ? -value : value;
}

[[nodiscard]] Rational maximum(const Rational& lhs, const Rational& rhs) {
    return lhs < rhs ? rhs : lhs;
}

[[nodiscard]] Rational absoluteErrorForBounds(
    const Rational& displayed,
    const Rational& lower,
    const Rational& upper) {
    const Rational lowerError = absRational(displayed - lower);
    const Rational upperError = absRational(displayed - upper);
    return maximum(lowerError, upperError);
}

[[nodiscard]] Rational minimumMagnitudeForBounds(
    const Rational& lower,
    const Rational& upper) {
    if (lower <= Rational{} && upper >= Rational{})
        return Rational{};
    return lower > Rational{} ? lower : -upper;
}

[[nodiscard]] std::size_t decimalIntegerDigits(const Rational& value) {
    const BigInt integerPart = value.numerator().abs() / value.denominator();
    if (integerPart.isZero())
        return 1;
    return integerPart.toString().size();
}

[[nodiscard]] std::size_t guaranteedDigits(const Rational& error, std::size_t cap) {
    if (error.isZero())
        return cap;
    if (cap == 0 || error >= Rational{BigInt{1}})
        return 0;

    const std::size_t numeratorDigits = error.numerator().abs().toString().size();
    const std::size_t denominatorDigits = error.denominator().toString().size();
    if (denominatorDigits <= numeratorDigits)
        return 0;

    std::size_t candidate = std::min(cap, denominatorDigits - numeratorDigits);
    const auto exponent = static_cast<std::uint64_t>(candidate);
    const BigInt scaledNumerator = error.numerator().abs()
        * pow(BigInt{10}, exponent);
    if (scaledNumerator < error.denominator())
        return candidate;
    return candidate == 0 ? 0 : candidate - 1;
}

[[nodiscard]] std::size_t guaranteedDigits(const Rational& error) {
    if (error.isZero() || error >= Rational{BigInt{1}})
        return 0;

    const std::size_t numeratorDigits = error.numerator().abs().toString().size();
    const std::size_t denominatorDigits = error.denominator().toString().size();
    if (denominatorDigits <= numeratorDigits)
        return 0;

    std::size_t candidate = denominatorDigits - numeratorDigits;
    const BigInt scaledNumerator = error.numerator().abs()
        * pow(BigInt{10}, static_cast<std::uint64_t>(candidate));
    if (scaledNumerator < error.denominator())
        return candidate;
    return candidate == 0 ? 0 : candidate - 1;
}

[[nodiscard]] Rational square(const Rational& value) {
    return value * value;
}

[[nodiscard]] std::size_t activeComplexCap(const ComplexDecimalApproximation& value) {
    std::size_t cap = std::numeric_limits<std::size_t>::max();
    bool any = false;

    if (!value.realInformationExactlyZero()) {
        cap = requestedDigitCap(value.real());
        any = true;
    }
    if (!value.imaginaryInformationExactlyZero()) {
        const std::size_t imaginaryCap = requestedDigitCap(value.imaginary());
        cap = any ? std::min(cap, imaginaryCap) : imaginaryCap;
        any = true;
    }

    if (any)
        return cap;
    return std::min(requestedDigitCap(value.real()), requestedDigitCap(value.imaginary()));
}

} // namespace

Rational informationAbsoluteError(const DecimalApproximation& value) {
    return absoluteErrorForBounds(
        value.displayedValue(), value.informationLower(), value.informationUpper());
}

Rational minimumInformationMagnitude(const DecimalApproximation& value) {
    return minimumMagnitudeForBounds(value.informationLower(), value.informationUpper());
}

std::size_t requestedDigitCap(const DecimalApproximation& value) {
    if (value.requestedSignificantDigits() != 0)
        return value.requestedSignificantDigits();
    return value.requestedFractionalDigits()
        + decimalIntegerDigits(value.displayedValue()) + 2;
}

std::size_t accuracyDigits(const DecimalApproximation& value) {
    const Rational error = informationAbsoluteError(value);
    if (error.isZero())
        return value.requestedSignificantDigits() != 0
            ? value.requestedSignificantDigits()
            : value.requestedFractionalDigits();
    return guaranteedDigits(error);
}

std::size_t precisionDigitsForBounds(
    const Rational& displayed,
    const Rational& lower,
    const Rational& upper,
    std::size_t cap) {
    const Rational minimumMagnitude = minimumMagnitudeForBounds(lower, upper);
    if (minimumMagnitude.isZero())
        return 0;

    const Rational error = absoluteErrorForBounds(displayed, lower, upper);
    if (error.isZero())
        return cap;
    return guaranteedDigits(error / minimumMagnitude, cap);
}

std::size_t precisionDigits(const DecimalApproximation& value) {
    return precisionDigitsForBounds(
        value.displayedValue(), value.informationLower(), value.informationUpper(),
        requestedDigitCap(value));
}

std::size_t accuracyDigits(const ComplexDecimalApproximation& value) {
    const Rational realError = value.realInformationExactlyZero()
        ? Rational{} : informationAbsoluteError(value.real());
    const Rational imaginaryError = value.imaginaryInformationExactlyZero()
        ? Rational{} : informationAbsoluteError(value.imaginary());
    const Rational errorSquared = square(realError) + square(imaginaryError);
    if (errorSquared.isZero()) {
        const std::size_t cap = activeComplexCap(value);
        return cap;
    }

    return guaranteedDigits(errorSquared) / 2;
}

std::size_t complexPrecisionDigitsForBounds(
    const Rational& realDisplayed,
    const Rational& realLower,
    const Rational& realUpper,
    bool realInformationExactlyZero,
    const Rational& imaginaryDisplayed,
    const Rational& imaginaryLower,
    const Rational& imaginaryUpper,
    bool imaginaryInformationExactlyZero,
    std::size_t cap) {
    const Rational realError = realInformationExactlyZero
        ? Rational{} : absoluteErrorForBounds(realDisplayed, realLower, realUpper);
    const Rational imaginaryError = imaginaryInformationExactlyZero
        ? Rational{} : absoluteErrorForBounds(imaginaryDisplayed, imaginaryLower, imaginaryUpper);
    const Rational errorSquared = square(realError) + square(imaginaryError);

    const Rational realMagnitude = realInformationExactlyZero
        ? Rational{} : minimumMagnitudeForBounds(realLower, realUpper);
    const Rational imaginaryMagnitude = imaginaryInformationExactlyZero
        ? Rational{} : minimumMagnitudeForBounds(imaginaryLower, imaginaryUpper);
    const Rational magnitudeSquared = square(realMagnitude) + square(imaginaryMagnitude);
    if (magnitudeSquared.isZero())
        return 0;
    if (errorSquared.isZero())
        return cap;

    const Rational relativeErrorSquared = errorSquared / magnitudeSquared;
    const std::size_t doubledCap = cap > std::numeric_limits<std::size_t>::max() / 2
        ? std::numeric_limits<std::size_t>::max()
        : cap * 2;
    return std::min(cap, guaranteedDigits(relativeErrorSquared, doubledCap) / 2);
}

std::size_t precisionDigits(const ComplexDecimalApproximation& value) {
    return complexPrecisionDigitsForBounds(
        value.real().displayedValue(),
        value.realInformationLower(), value.realInformationUpper(),
        value.realInformationExactlyZero(),
        value.imaginary().displayedValue(),
        value.imaginaryInformationLower(), value.imaginaryInformationUpper(),
        value.imaginaryInformationExactlyZero(),
        activeComplexCap(value));
}

} // namespace mmcal::numeric
