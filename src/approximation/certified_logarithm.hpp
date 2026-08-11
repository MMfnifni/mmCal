#pragma once

#include "real_interval.hpp"

#include <cstddef>

namespace mmcal::approximation {

struct CertifiedLogarithmResult final {
    RealInterval interval;
    std::size_t termsUsed = 0;
};

// 正の実区間上の自然対数を包含保証付きで評価する。0以下を含む区間はreal Logの定義域外なので拒否する。
[[nodiscard]] CertifiedLogarithmResult encloseLogPositive(
    const RealInterval& input,
    std::size_t precisionBits);

} // namespace mmcal::approximation
