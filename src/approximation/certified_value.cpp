#include "certified_value.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::approximation {

CertifiedValue::CertifiedValue(RealInterval real)
    : value_(std::move(real)) {}

CertifiedValue::CertifiedValue(ComplexInterval complex)
    : value_(std::move(complex)) {}

bool CertifiedValue::isReal() const noexcept {
    return std::holds_alternative<RealInterval>(value_);
}

bool CertifiedValue::isComplex() const noexcept {
    return std::holds_alternative<ComplexInterval>(value_);
}

const RealInterval& CertifiedValue::asReal() const {
    if (!isReal())
        throw std::logic_error("CertifiedValue does not contain a real interval");
    return std::get<RealInterval>(value_);
}

const ComplexInterval& CertifiedValue::asComplex() const {
    if (!isComplex())
        throw std::logic_error("CertifiedValue does not contain a complex interval");
    return std::get<ComplexInterval>(value_);
}

ComplexInterval CertifiedValue::toComplex() const {
    if (isComplex())
        return asComplex();
    return ComplexInterval::fromReal(asReal());
}

} // namespace mmcal::approximation
