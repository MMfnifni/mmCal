#pragma once

#include "certified_value.hpp"

#include <cstddef>

namespace mmcal::approximation {

// CertifiedValueのReal/Complex昇格規則を一箇所に集約する。
// Complex -> Realの縮約は虚部が厳密なpoint zeroと証明できる場合だけ行う。
[[nodiscard]] CertifiedValue normalizeCertifiedComplex(ComplexInterval value);
[[nodiscard]] CertifiedValue addCertifiedValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits);
[[nodiscard]] CertifiedValue subtractCertifiedValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits);
[[nodiscard]] CertifiedValue multiplyCertifiedValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits);
[[nodiscard]] CertifiedValue divideCertifiedValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits);
[[nodiscard]] CertifiedValue negateCertifiedValue(const CertifiedValue& value);

} // namespace mmcal::approximation
