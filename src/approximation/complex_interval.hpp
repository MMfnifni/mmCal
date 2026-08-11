#pragma once

#include "real_interval.hpp"

#include <cstddef>

namespace mmcal::approximation {

// 真の複素数 z = x + yI を、実部・虚部それぞれのcertified区間で囲う。
//
// RealIntervalと同様に「中心値と誤差」ではなく包含区間そのものを持つ。
// したがって、例えば虚部区間が厳密に [0,0] のときだけ「実数へ縮約できる」と証明できる。
// 単に虚部が小さいという理由で実数扱いすることはしない。
class ComplexInterval final {
public:
    ComplexInterval(RealInterval real, RealInterval imaginary);

    [[nodiscard]] static ComplexInterval fromReal(RealInterval real);

    [[nodiscard]] const RealInterval& real() const noexcept;
    [[nodiscard]] const RealInterval& imaginary() const noexcept;

    [[nodiscard]] bool isProvablyReal() const noexcept;
    [[nodiscard]] bool containsZero() const noexcept;

    [[nodiscard]] ComplexInterval roundedOutward(std::size_t precisionBits) const;

private:
    RealInterval real_;
    RealInterval imaginary_;
};

[[nodiscard]] ComplexInterval negate(const ComplexInterval& value);

[[nodiscard]] ComplexInterval add(
    const ComplexInterval& lhs,
    const ComplexInterval& rhs,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval subtract(
    const ComplexInterval& lhs,
    const ComplexInterval& rhs,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval multiply(
    const ComplexInterval& lhs,
    const ComplexInterval& rhs,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval divide(
    const ComplexInterval& lhs,
    const ComplexInterval& rhs,
    std::size_t precisionBits);

} // namespace mmcal::approximation
