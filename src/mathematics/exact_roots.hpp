#pragma once

#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "numeric/real_number.hpp"

#include <optional>

namespace mmcal::mathematics {

// 正の有理数 q に対して
//     sqrt(q) = coefficient * sqrt(radicand)
// となる形へ、証明できた平方因子だけを外へ出す。radicandは正の整数。
// 分母の平方根は有理化してradicand側へ移すため、表示上もsqrt(2/3) -> sqrt(6)/3 のような形を作れる。
//
// 巨大整数の完全素因数分解は行わないため、非常に大きな未知因子については
// 「完全なsquare-free化」ではなく安全な部分正規化になる。等式自体は常にexact。
struct RationalSquareRootDecomposition final {
    numeric::Rational coefficient;
    numeric::BigInt radicand;
};

[[nodiscard]] RationalSquareRootDecomposition decomposePositiveRationalSquareRoot(
    const numeric::Rational& value);


// exactな有理実数について、平方根が現在のRealNumber domainで閉じる場合だけ返す。
// 例: 4 -> 2, 9/16 -> 3/4, 2 -> nullopt。
[[nodiscard]] std::optional<numeric::RealNumber> exactSquareRoot(
    const numeric::RealNumber& value);

// exactな有理実数について、実立方根が現在のRealNumber domainで閉じる場合だけ返す。
// 負数では実根を選ぶため principal Power[x,1/3] とは別意味を持つ。
[[nodiscard]] std::optional<numeric::RealNumber> exactRealCubeRoot(
    const numeric::RealNumber& value);

// principal square root が現在のNumber domain（有理実部・有理虚部）で閉じる場合だけ返す。
// principal sqrt は常に実部 >= 0 を選び、実部が0なら虚部 >= 0 を選ぶ。
//
// 例:
//   sqrt(3 + 4I)   = 2 + I
//   sqrt(-3 + 4I)  = 1 + 2I
//   sqrt(3 - 4I)   = 2 - I
//   sqrt(-4)       = 2I
//   sqrt(1 + I)    = nullopt  // 係数にsqrt(2)が必要でNumberだけでは閉じない
[[nodiscard]] std::optional<numeric::Number> exactPrincipalSquareRoot(
    const numeric::Number& value);

} // namespace mmcal::mathematics
