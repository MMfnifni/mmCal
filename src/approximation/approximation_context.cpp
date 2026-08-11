// 近似精度とガード桁の管理
#include "approximation_context.hpp"
#include "precision.hpp"

#include <limits>
#include <stdexcept>

namespace mmcal::approximation {

ApproximationContext::ApproximationContext(
    std::size_t decimalDigits,
    std::size_t guardDigits,
    RoundingMode roundingMode)
    : guardDigits_(guardDigits), roundingMode_(roundingMode) {
    setDecimalDigits(decimalDigits);
}

std::size_t ApproximationContext::decimalDigits() const noexcept {
    return decimalDigits_;
}

std::size_t ApproximationContext::guardDigits() const noexcept {
    return guardDigits_;
}

std::size_t ApproximationContext::workingDecimalDigits() const {
    if (guardDigits_ > std::numeric_limits<std::size_t>::max() - decimalDigits_)
        throw std::overflow_error("Approximation working precision is too large");

    return decimalDigits_ + guardDigits_;
}

std::size_t ApproximationContext::workingBinaryBits() const {
    return decimalDigitsToBinaryBits(workingDecimalDigits());
}

RoundingMode ApproximationContext::roundingMode() const noexcept {
    return roundingMode_;
}

void ApproximationContext::setDecimalDigits(std::size_t decimalDigits) {
    if (decimalDigits == 0)
        throw std::invalid_argument("Approximation precision must be greater than zero");

    decimalDigits_ = decimalDigits;
}

void ApproximationContext::setGuardDigits(std::size_t guardDigits) noexcept {
    guardDigits_ = guardDigits;
}

void ApproximationContext::setRoundingMode(RoundingMode roundingMode) noexcept {
    roundingMode_ = roundingMode;
}

} // namespace mmcal::approximation
