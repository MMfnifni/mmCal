// 線形代数内部のBigFloat複素点演算
#include "complex_point.hpp"
#include "point_arithmetic.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "approximation/certified_sqrt.hpp"
#include "approximation/real_interval.hpp"

namespace mmcal::linear_algebra::detail {
namespace {

using numeric::BigFloat;
using numeric::BigInt;

[[nodiscard]] BigFloat squareRoot(const BigFloat& value, std::size_t bits) {
    const auto enclosure = approximation::encloseSqrt(value.toRational(), bits + 8).interval;
    return midpoint(enclosure, bits);
}

} // namespace

ComplexPoint complexZero(std::size_t bits) {
    return ComplexPoint{zero(bits), zero(bits)};
}

ComplexPoint complexOne(std::size_t bits) {
    return ComplexPoint{one(bits), zero(bits)};
}

ComplexPoint complexConjugate(const ComplexPoint& value) {
    return ComplexPoint{value.real, -value.imaginary};
}

ComplexPoint complexNegate(const ComplexPoint& value) {
    return ComplexPoint{-value.real, -value.imaginary};
}

ComplexPoint complexAdd(
    const ComplexPoint& lhs, const ComplexPoint& rhs, std::size_t bits) {
    return ComplexPoint{
        add(lhs.real, rhs.real, bits),
        add(lhs.imaginary, rhs.imaginary, bits)
    };
}

ComplexPoint complexSubtract(
    const ComplexPoint& lhs, const ComplexPoint& rhs, std::size_t bits) {
    return ComplexPoint{
        subtract(lhs.real, rhs.real, bits),
        subtract(lhs.imaginary, rhs.imaginary, bits)
    };
}

ComplexPoint complexMultiply(
    const ComplexPoint& lhs, const ComplexPoint& rhs, std::size_t bits) {
    const BigFloat ac = multiply(lhs.real, rhs.real, bits);
    const BigFloat bd = multiply(lhs.imaginary, rhs.imaginary, bits);
    const BigFloat ad = multiply(lhs.real, rhs.imaginary, bits);
    const BigFloat bc = multiply(lhs.imaginary, rhs.real, bits);
    return ComplexPoint{
        subtract(ac, bd, bits),
        add(ad, bc, bits)
    };
}

ComplexPoint complexScale(
    const ComplexPoint& value, const BigFloat& factor, std::size_t bits) {
    return ComplexPoint{
        multiply(value.real, factor, bits),
        multiply(value.imaginary, factor, bits)
    };
}

ComplexPoint complexDivideReal(
    const ComplexPoint& value, const BigFloat& divisor, std::size_t bits) {
    return ComplexPoint{
        divide(value.real, divisor, bits),
        divide(value.imaginary, divisor, bits)
    };
}

ComplexPoint complexDivide(
    const ComplexPoint& numerator, const ComplexPoint& denominator, std::size_t bits) {
    const BigFloat denominatorSquared = complexMagnitudeSquared(denominator, bits);
    const ComplexPoint quotientNumerator = complexMultiply(
        numerator, complexConjugate(denominator), bits);
    return complexDivideReal(quotientNumerator, denominatorSquared, bits);
}

BigFloat complexMagnitude(
    const ComplexPoint& value, std::size_t bits) {
    return squareRoot(complexMagnitudeSquared(value, bits), bits);
}

BigFloat complexMagnitudeSquared(const ComplexPoint& value, std::size_t bits) {
    return add(
        multiply(value.real, value.real, bits),
        multiply(value.imaginary, value.imaginary, bits), bits);
}

bool complexIsZero(const ComplexPoint& value) noexcept {
    return value.real.isZero() && value.imaginary.isZero();
}

} // namespace mmcal::linear_algebra::detail
