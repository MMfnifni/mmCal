#pragma once

#include "../big_int.hpp"

#include <cstdint>

namespace mmcal::numeric::detail {

// 正の整数比 numerator / denominator に対して floor(log2()) を厳密に求める。
//
// 浮動小数点のlog2は一切使わず、BigIntのbit長と整数比較だけで決定する。
// BigFloatのRational変換だけでなく、平方根など「2進桁位置を数学的に決める」
// 算法で共通利用する。
[[nodiscard]] std::int64_t floorLog2PositiveRatio(
    const BigInt& numerator,
    const BigInt& denominator);

} // namespace mmcal::numeric::detail
