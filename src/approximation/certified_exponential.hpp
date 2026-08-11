#pragma once

#include "real_interval.hpp"

#include <cstddef>

namespace mmcal::approximation {

struct CertifiedExponentialResult final {
    RealInterval interval;
    std::size_t termsUsed = 0;
    std::size_t squarings = 0;
};

// 実区間上のexpを包含保証付きで評価する。expは単調増加なので、区間端点をそれぞれcertified評価して外側端点を採る。
[[nodiscard]] CertifiedExponentialResult encloseExp(
    const RealInterval& input,
    std::size_t precisionBits);

} // namespace mmcal::approximation
