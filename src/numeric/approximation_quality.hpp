#pragma once

#include "complex_decimal_approximation.hpp"
#include "decimal_approximation.hpp"

#include <cstddef>

namespace mmcal::numeric {

// InformationEnclosureから導出する近似値の品質指標を一元化する。
// requested precisionは現行仕様どおり上限として扱う。表示規則とは独立である。
[[nodiscard]] Rational informationAbsoluteError(const DecimalApproximation& value);
[[nodiscard]] Rational minimumInformationMagnitude(const DecimalApproximation& value);
[[nodiscard]] std::size_t requestedDigitCap(const DecimalApproximation& value);
[[nodiscard]] std::size_t accuracyDigits(const DecimalApproximation& value);
[[nodiscard]] std::size_t precisionDigits(const DecimalApproximation& value);
[[nodiscard]] std::size_t accuracyDigits(const ComplexDecimalApproximation& value);
[[nodiscard]] std::size_t precisionDigits(const ComplexDecimalApproximation& value);

// stored DecimalApproximation以外のcertified/information intervalにも同じ品質基準を適用する。
[[nodiscard]] std::size_t precisionDigitsForBounds(
    const Rational& displayed,
    const Rational& lower,
    const Rational& upper,
    std::size_t cap);

[[nodiscard]] std::size_t complexPrecisionDigitsForBounds(
    const Rational& realDisplayed,
    const Rational& realLower,
    const Rational& realUpper,
    bool realInformationExactlyZero,
    const Rational& imaginaryDisplayed,
    const Rational& imaginaryLower,
    const Rational& imaginaryUpper,
    bool imaginaryInformationExactlyZero,
    std::size_t cap);

} // namespace mmcal::numeric
