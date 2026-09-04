#pragma once

#include "approximation/real_interval.hpp"
#include "numeric/big_float.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <cstddef>

namespace mmcal::linear_algebra::detail {

// SVD/Eigenの中心値計算は保証区間ではなくNearestEvenのBigFloat点を使う。
// 丸めmodeを各算法へ散らすと同じ演算でも中心値がずれるため，薄い点演算をここへ固定する。
[[nodiscard]] inline numeric::BigFloat zero(std::size_t bits) {
    return numeric::BigFloat::fromBigInt(
        numeric::BigInt{}, bits, numeric::RoundingMode::NearestEven);
}

[[nodiscard]] inline numeric::BigFloat one(std::size_t bits) {
    return numeric::BigFloat::fromBigInt(
        numeric::BigInt{1}, bits, numeric::RoundingMode::NearestEven);
}

[[nodiscard]] inline numeric::BigFloat two(std::size_t bits) {
    return numeric::BigFloat::fromBigInt(
        numeric::BigInt{2}, bits, numeric::RoundingMode::NearestEven);
}

[[nodiscard]] inline numeric::BigFloat add(
    const numeric::BigFloat& lhs,
    const numeric::BigFloat& rhs,
    std::size_t bits) {
    return numeric::add(lhs, rhs, bits, numeric::RoundingMode::NearestEven);
}

[[nodiscard]] inline numeric::BigFloat subtract(
    const numeric::BigFloat& lhs,
    const numeric::BigFloat& rhs,
    std::size_t bits) {
    return numeric::subtract(lhs, rhs, bits, numeric::RoundingMode::NearestEven);
}

[[nodiscard]] inline numeric::BigFloat multiply(
    const numeric::BigFloat& lhs,
    const numeric::BigFloat& rhs,
    std::size_t bits) {
    return numeric::multiply(lhs, rhs, bits, numeric::RoundingMode::NearestEven);
}

[[nodiscard]] inline numeric::BigFloat divide(
    const numeric::BigFloat& lhs,
    const numeric::BigFloat& rhs,
    std::size_t bits) {
    return numeric::divide(lhs, rhs, bits, numeric::RoundingMode::NearestEven);
}

[[nodiscard]] inline numeric::BigFloat absolute(const numeric::BigFloat& value) {
    return value.isNegative() ? -value : value;
}

[[nodiscard]] inline numeric::Rational midpointRational(
    const approximation::RealInterval& interval) {
    return (interval.lower().toRational() + interval.upper().toRational())
        / numeric::Rational{numeric::BigInt{2}};
}

[[nodiscard]] inline numeric::BigFloat midpoint(
    const approximation::RealInterval& interval,
    std::size_t bits) {
    return numeric::BigFloat::fromRational(
        midpointRational(interval), bits, numeric::RoundingMode::NearestEven);
}

} // namespace mmcal::linear_algebra::detail
