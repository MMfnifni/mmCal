#pragma once

#include <cstddef>

namespace mmcal::approximation {

// 10進n桁を保持するのに十分な2進bit数の安全な上界を返す。
// log2(10) ~= 3.321928... に対して 3.322 = 3322/1000 を使うため、戻り値は必要bit数を下回らない。
[[nodiscard]] std::size_t decimalDigitsToBinaryBits(std::size_t decimalDigits);

} // namespace mmcal::approximation
