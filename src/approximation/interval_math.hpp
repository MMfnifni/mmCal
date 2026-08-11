#pragma once

#include "real_interval.hpp"

#include <cstddef>

namespace mmcal::approximation {

// 同一変数の二乗 x^2 を通常の区間積 x*x より鋭く囲う。特に0を跨ぐ区間 [a,b] では下限を厳密に0へ固定できる。
[[nodiscard]] RealInterval squareInterval(
    const RealInterval& value,
    std::size_t precisionBits);

// 区間 |x| の値域を囲う。0を跨ぐ場合は [0,max(|a|,|b|)]。
[[nodiscard]] RealInterval absoluteInterval(
    const RealInterval& value,
    std::size_t precisionBits);

} // namespace mmcal::approximation
