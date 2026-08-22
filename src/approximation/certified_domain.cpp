#include "certified_domain.hpp"

#include <cstdint>

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "numeric/rational_rounding.hpp"

namespace mmcal::approximation {
namespace {

using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational rational(std::int64_t value) {
    return Rational{BigInt{value}};
}

[[nodiscard]] bool realIntervalMayContainNonPositiveInteger(
    const RealInterval& value) {
    const Rational lower = value.lower().toRational();
    const Rational upper = value.upper().toRational();
    if (lower > rational(0))
        return false;

    const Rational clippedUpper = upper < rational(0) ? upper : rational(0);
    return numeric::ceilToInteger(lower) <= numeric::floorToInteger(clippedUpper);
}

} // namespace

bool intervalIsExactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

bool intervalIsExactZero(const ComplexInterval& value) noexcept {
    return intervalIsExactZero(value.real()) && intervalIsExactZero(value.imaginary());
}

SingularityRelation classifyNonPositiveIntegerPole(const RealInterval& value) {
    if (value.isPoint()) {
        const Rational point = value.lower().toRational();
        if (point.isInteger() && point <= rational(0))
            return SingularityRelation::ExactSingularity;
    }
    return realIntervalMayContainNonPositiveInteger(value)
        ? SingularityRelation::MayContainSingularity
        : SingularityRelation::Clear;
}

SingularityRelation classifyNonPositiveIntegerPole(const ComplexInterval& value) {
    if (!value.imaginary().containsZero())
        return SingularityRelation::Clear;

    if (intervalIsExactZero(value.imaginary())) {
        const SingularityRelation realRelation = classifyNonPositiveIntegerPole(value.real());
        if (realRelation == SingularityRelation::ExactSingularity)
            return SingularityRelation::ExactSingularity;
        if (realRelation == SingularityRelation::Clear)
            return SingularityRelation::Clear;
    }

    return realIntervalMayContainNonPositiveInteger(value.real())
        ? SingularityRelation::MayContainSingularity
        : SingularityRelation::Clear;
}

bool informationCrossesPrincipalNegativeRealCut(const ComplexInterval& value) {
    const numeric::BigFloat zero;
    return value.real().lower() < zero
        && value.imaginary().lower() < zero
        && value.imaginary().upper() >= zero;
}

bool informationIsAmbiguousAtOuterRealCuts(const ComplexInterval& value) {
    const numeric::BigFloat zero;
    const Rational lower = value.real().lower().toRational();
    const Rational upper = value.real().upper().toRational();
    const RealInterval& imaginary = value.imaginary();
    const bool positiveCutWrongSide = upper >= rational(1)
        && imaginary.lower() <= zero && imaginary.upper() > zero;
    const bool negativeCutWrongSide = lower <= rational(-1)
        && imaginary.lower() < zero && imaginary.upper() >= zero;
    return positiveCutWrongSide || negativeCutWrongSide;
}

bool informationIsAmbiguousAtAcoshCut(const ComplexInterval& value) {
    const numeric::BigFloat zero;
    const RealInterval& imaginary = value.imaginary();
    return value.real().lower().toRational() <= rational(1)
        && imaginary.lower() < zero && imaginary.upper() >= zero;
}

bool informationIsAmbiguousAtOuterImaginaryCuts(const ComplexInterval& value) {
    const numeric::BigFloat zero;
    const Rational lower = value.imaginary().lower().toRational();
    const Rational upper = value.imaginary().upper().toRational();
    const RealInterval& real = value.real();
    const bool positiveCutWrongSide = upper >= rational(1)
        && real.lower() < zero && real.upper() >= zero;
    const bool negativeCutWrongSide = lower <= rational(-1)
        && real.lower() <= zero && real.upper() > zero;
    return positiveCutWrongSide || negativeCutWrongSide;
}

bool informationIsAmbiguousAtPositiveRealCut(const ComplexInterval& value) {
    const numeric::BigFloat zero;
    const RealInterval& imaginary = value.imaginary();
    return value.real().upper().toRational() >= rational(1)
        && imaginary.lower() <= zero && imaginary.upper() > zero;
}

bool informationMayContainComplexZero(const ComplexInterval& value) noexcept {
    return value.containsZero() && !intervalIsExactZero(value);
}

} // namespace mmcal::approximation
