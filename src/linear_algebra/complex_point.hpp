#pragma once

#include "numeric/big_float.hpp"

#include <cstddef>

namespace mmcal::linear_algebra::detail {

// BigFloat二成分だけを持つ線形代数内部用の複素点。
// certified enclosureそのものではなく，反復算法の中心値計算にだけ使う。
struct ComplexPoint final {
    numeric::BigFloat real;
    numeric::BigFloat imaginary;
};

[[nodiscard]] ComplexPoint complexZero(std::size_t bits);
[[nodiscard]] ComplexPoint complexOne(std::size_t bits);
[[nodiscard]] ComplexPoint complexConjugate(const ComplexPoint& value);
[[nodiscard]] ComplexPoint complexNegate(const ComplexPoint& value);
[[nodiscard]] ComplexPoint complexAdd(
    const ComplexPoint& lhs, const ComplexPoint& rhs, std::size_t bits);
[[nodiscard]] ComplexPoint complexSubtract(
    const ComplexPoint& lhs, const ComplexPoint& rhs, std::size_t bits);
[[nodiscard]] ComplexPoint complexMultiply(
    const ComplexPoint& lhs, const ComplexPoint& rhs, std::size_t bits);
[[nodiscard]] ComplexPoint complexScale(
    const ComplexPoint& value, const numeric::BigFloat& factor, std::size_t bits);
[[nodiscard]] ComplexPoint complexDivideReal(
    const ComplexPoint& value, const numeric::BigFloat& divisor, std::size_t bits);
[[nodiscard]] ComplexPoint complexDivide(
    const ComplexPoint& numerator, const ComplexPoint& denominator, std::size_t bits);
[[nodiscard]] numeric::BigFloat complexMagnitude(
    const ComplexPoint& value, std::size_t bits);
[[nodiscard]] numeric::BigFloat complexMagnitudeSquared(
    const ComplexPoint& value, std::size_t bits);
[[nodiscard]] bool complexIsZero(const ComplexPoint& value) noexcept;

} // namespace mmcal::linear_algebra::detail
