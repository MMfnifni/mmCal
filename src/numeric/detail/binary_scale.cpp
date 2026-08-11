// 二進scale補助
#include "binary_scale.hpp"

#include <limits>
#include <stdexcept>

namespace mmcal::numeric::detail {

std::int64_t floorLog2PositiveRatio(
    const BigInt& numerator,
    const BigInt& denominator) {
    if (!numerator.isPositive() || !denominator.isPositive())
        throw std::invalid_argument("Binary scale requires a positive ratio");

    const std::size_t numeratorBits = numerator.bitLength();
    const std::size_t denominatorBits = denominator.bitLength();

    if (numeratorBits >= denominatorBits) {
        const std::size_t difference = numeratorBits - denominatorBits;
        if (difference > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
            throw std::overflow_error("Binary scale exponent overflow");

        const BigInt boundary = denominator << difference;
        auto result = static_cast<std::int64_t>(difference);
        if (numerator < boundary)
            --result;
        return result;
    }

    const std::size_t difference = denominatorBits - numeratorBits;
    if (difference > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::overflow_error("Binary scale exponent underflow");

    const BigInt scaledNumerator = numerator << difference;
    auto result = -static_cast<std::int64_t>(difference);
    if (scaledNumerator < denominator)
        --result;
    return result;
}

} // namespace mmcal::numeric::detail
