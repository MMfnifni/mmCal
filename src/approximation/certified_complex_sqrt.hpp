#pragma once

#include "complex_interval.hpp"

#include <cstddef>

namespace mmcal::approximation {

// ComplexIntervalが表す任意の z に対するprincipal sqrt(z)を必ず含む長方形区間を返す。
// principal branchは
//   Re(sqrt(z)) >= 0
// を採り、負の実軸上では虚部が正の値を選ぶ。
//
// 入力区間がbranch cut（負実軸）を跨ぐ場合、principal sqrtは上下半平面からの
// 極限で虚部の符号が反転するため、虚部区間を両符号へ広げて包含を保証する。
[[nodiscard]] ComplexInterval enclosePrincipalComplexSqrt(
    const ComplexInterval& value,
    std::size_t precisionBits);

} // namespace mmcal::approximation
