#pragma once

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace mmcal::approximation {

[[nodiscard]] inline std::size_t checkedPrecisionAdd(
    std::size_t lhs,
    std::size_t rhs,
    const char* message) {
    if (rhs > std::numeric_limits<std::size_t>::max() - lhs)
        throw std::overflow_error(message);
    return lhs + rhs;
}

[[nodiscard]] inline numeric::Rational binaryPrecisionThreshold(std::size_t bits) {
    numeric::BigInt denominator{1};
    denominator <<= bits;
    return numeric::Rational{numeric::BigInt{1}, std::move(denominator)};
}

} // namespace mmcal::approximation
