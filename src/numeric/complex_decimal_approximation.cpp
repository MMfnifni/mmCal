// 複素近似値metadata
#include "complex_decimal_approximation.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::numeric {
namespace {

[[nodiscard]] bool negativeText(std::string_view text) noexcept {
    return !text.empty() && text.front() == '-';
}

[[nodiscard]] std::string_view magnitudeText(std::string_view text) noexcept {
    return negativeText(text) ? text.substr(1) : text;
}

[[nodiscard]] std::string imaginaryText(const DecimalApproximation& value) {
    const std::string_view text = value.text();
    const std::string_view magnitude = magnitudeText(text);

    std::string result;
    if (negativeText(text))
        result.push_back('-');

    // exactな1なら従来のNumberと同様に係数を省略する。
    // certified interval由来の1.0は近似値であることを示すため係数を省略しない。
    if (magnitude != "1")
        result += magnitude;
    result.push_back('I');
    return result;
}

} // namespace

ComplexDecimalApproximation::ComplexDecimalApproximation(
    std::string text,
    DecimalApproximation real,
    DecimalApproximation imaginary,
    bool realExactlyZero,
    bool imaginaryExactlyZero)
    : text_(std::move(text)),
      real_(std::move(real)),
      imaginary_(std::move(imaginary)),
      realExactlyZero_(realExactlyZero),
      imaginaryExactlyZero_(imaginaryExactlyZero) {
    if (text_.empty())
        throw std::invalid_argument("Complex decimal approximation text cannot be empty");
}

ComplexDecimalApproximation ComplexDecimalApproximation::fromComponents(
    DecimalApproximation real,
    DecimalApproximation imaginary,
    bool realExactlyZero,
    bool imaginaryExactlyZero) {
    std::string text;
    if (imaginaryExactlyZero)
        text = std::string{real.text()};
    else if (realExactlyZero)
        text = imaginaryText(imaginary);
    else {
        const std::string_view imaginaryValue = imaginary.text();
        const bool negativeImaginary = negativeText(imaginaryValue);
        const std::string_view magnitude = magnitudeText(imaginaryValue);
        text = std::string{real.text()};
        text += negativeImaginary ? "-" : "+";
        if (magnitude != "1")
            text += magnitude;
        text.push_back('I');
    }

    return ComplexDecimalApproximation{
        std::move(text), std::move(real), std::move(imaginary),
        realExactlyZero, imaginaryExactlyZero};
}

std::string_view ComplexDecimalApproximation::text() const noexcept {
    return text_;
}

const DecimalApproximation& ComplexDecimalApproximation::real() const noexcept {
    return real_;
}

const DecimalApproximation& ComplexDecimalApproximation::imaginary() const noexcept {
    return imaginary_;
}

bool ComplexDecimalApproximation::realExactlyZero() const noexcept {
    return realExactlyZero_;
}

bool ComplexDecimalApproximation::imaginaryExactlyZero() const noexcept {
    return imaginaryExactlyZero_;
}

} // namespace mmcal::numeric
