// 精度・誤差計算
#include "precision.hpp"

#include <limits>
#include <stdexcept>

namespace mmcal::approximation {

std::size_t decimalDigitsToBinaryBits(std::size_t decimalDigits) {
    if (decimalDigits == 0)
        return 1;

    constexpr std::size_t numerator = 3322;
    constexpr std::size_t denominator = 1000;

    const std::size_t whole = decimalDigits / denominator;
    const std::size_t remainder = decimalDigits % denominator;

    if (whole > std::numeric_limits<std::size_t>::max() / numerator)
        throw std::overflow_error("Requested precision is too large");

    std::size_t bits = whole * numerator;
    const std::size_t remainderProduct = remainder * numerator;
    const std::size_t remainderBits =
        remainderProduct / denominator
        + (remainderProduct % denominator != 0 ? 1 : 0);

    if (bits > std::numeric_limits<std::size_t>::max() - remainderBits)
        throw std::overflow_error("Requested precision is too large");
    bits += remainderBits;
    return bits == 0 ? 1 : bits;
}

} // namespace mmcal::approximation
