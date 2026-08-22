#pragma once

#include "big_int.hpp"
#include "rational.hpp"

namespace mmcal::numeric {

// Rationalを0方向・床・天井・最近接偶数丸めでBigIntへ写す。
// C++の整数除算は0方向へ丸めるため，負値のfloor/ceilは剰余を見て補正する。
[[nodiscard]] BigInt truncateToInteger(const Rational& value);
[[nodiscard]] BigInt floorToInteger(const Rational& value);
[[nodiscard]] BigInt ceilToInteger(const Rational& value);
[[nodiscard]] BigInt roundToNearestEvenInteger(const Rational& value);

} // namespace mmcal::numeric
