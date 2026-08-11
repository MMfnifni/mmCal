#pragma once

#include "real_interval.hpp"

#include <cstddef>

namespace mmcal::approximation {

struct CertifiedAtanResult final {
    RealInterval interval;
    std::size_t termsUsed = 0;
};

// 実数区間上のprincipal atanをcertifiedに囲う。
// atanは実軸全体で単調増加なので、入力区間の両端を独立に囲えばよい。
[[nodiscard]] CertifiedAtanResult encloseAtan(
    const RealInterval& input,
    std::size_t precisionBits);

} // namespace mmcal::approximation
