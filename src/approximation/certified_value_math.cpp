#include "certified_value_math.hpp"

#include <utility>

namespace mmcal::approximation {

CertifiedValue normalizeCertifiedComplex(ComplexInterval value) {
    if (value.isProvablyReal())
        return CertifiedValue{value.real()};
    return CertifiedValue{std::move(value)};
}

CertifiedValue addCertifiedValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{add(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeCertifiedComplex(add(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

CertifiedValue subtractCertifiedValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{subtract(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeCertifiedComplex(subtract(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

CertifiedValue multiplyCertifiedValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{multiply(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeCertifiedComplex(multiply(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

CertifiedValue divideCertifiedValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{divide(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeCertifiedComplex(divide(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

CertifiedValue negateCertifiedValue(const CertifiedValue& value) {
    if (value.isReal())
        return CertifiedValue{negate(value.asReal())};
    return normalizeCertifiedComplex(negate(value.asComplex()));
}

} // namespace mmcal::approximation
