#pragma once

#include "real_interval.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/rational.hpp"
#include "numeric/real_number.hpp"

#include <cstddef>

namespace mmcal::approximation {

struct CertifiedSqrtEnclosure final {
    RealInterval interval;
    std::size_t precisionBits = 0;
};

// 非負のexact Rational xに対して sqrt(x) を必ず含むdyadic区間を作る。
// Taylor/Newtonの「収束したように見える」判定は使わず、整数平方根から直接上下界を作る。
[[nodiscard]] CertifiedSqrtEnclosure encloseSqrt(
    const numeric::Rational& value,
    std::size_t precisionBits);

// 非負実数区間 [a,b] に対し、sqrtの単調性から [sqrt(a),sqrt(b)] を外向きに囲う。
[[nodiscard]] CertifiedSqrtEnclosure encloseSqrt(
    const RealInterval& value,
    std::size_t precisionBits);

// 要求小数桁への丸めが上下端で一致するまで作業精度を増やす。
[[nodiscard]] numeric::DecimalApproximation approximateSqrt(
    const numeric::RealNumber& value,
    std::size_t fractionalDigits);

} // namespace mmcal::approximation
