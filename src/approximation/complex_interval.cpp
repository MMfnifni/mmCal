// 複素区間
#include "complex_interval.hpp"

#include "interval_math.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::approximation {
namespace {

[[nodiscard]] RealInterval zeroInterval(std::size_t precisionBits) {
    return RealInterval::fromRational(
        numeric::Rational{numeric::BigInt{0}}, precisionBits);
}

[[nodiscard]] bool isExactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

} // namespace

ComplexInterval::ComplexInterval(RealInterval real, RealInterval imaginary)
    : real_(std::move(real)), imaginary_(std::move(imaginary)) {}

ComplexInterval ComplexInterval::fromReal(RealInterval real) {
    const std::size_t precisionBits = real.lower().precisionBits();
    return ComplexInterval{std::move(real), zeroInterval(precisionBits)};
}

const RealInterval& ComplexInterval::real() const noexcept {
    return real_;
}

const RealInterval& ComplexInterval::imaginary() const noexcept {
    return imaginary_;
}

bool ComplexInterval::isProvablyReal() const noexcept {
    return isExactZero(imaginary_);
}

bool ComplexInterval::containsZero() const noexcept {
    return real_.containsZero() && imaginary_.containsZero();
}

ComplexInterval ComplexInterval::roundedOutward(std::size_t precisionBits) const {
    return ComplexInterval{
        real_.roundedOutward(precisionBits),
        imaginary_.roundedOutward(precisionBits)
    };
}

ComplexInterval negate(const ComplexInterval& value) {
    return ComplexInterval{
        approximation::negate(value.real()),
        approximation::negate(value.imaginary())
    };
}

ComplexInterval add(
    const ComplexInterval& lhs,
    const ComplexInterval& rhs,
    std::size_t precisionBits) {
    return ComplexInterval{
        approximation::add(lhs.real(), rhs.real(), precisionBits),
        approximation::add(lhs.imaginary(), rhs.imaginary(), precisionBits)
    };
}

ComplexInterval subtract(
    const ComplexInterval& lhs,
    const ComplexInterval& rhs,
    std::size_t precisionBits) {
    return ComplexInterval{
        approximation::subtract(lhs.real(), rhs.real(), precisionBits),
        approximation::subtract(lhs.imaginary(), rhs.imaginary(), precisionBits)
    };
}

ComplexInterval multiply(
    const ComplexInterval& lhs,
    const ComplexInterval& rhs,
    std::size_t precisionBits) {
    // (a + bI)(c + dI) = (ac - bd) + (ad + bc)I。
    // 各a,b,c,dは単一値ではなく区間だが、RealIntervalの四則演算が常に外向きに
    // 丸めるため、この代数式をそのまま区間演算へ持ち上げれば真の積を失わない。
    const RealInterval ac = approximation::multiply(
        lhs.real(), rhs.real(), precisionBits);
    const RealInterval bd = approximation::multiply(
        lhs.imaginary(), rhs.imaginary(), precisionBits);
    const RealInterval ad = approximation::multiply(
        lhs.real(), rhs.imaginary(), precisionBits);
    const RealInterval bc = approximation::multiply(
        lhs.imaginary(), rhs.real(), precisionBits);

    return ComplexInterval{
        approximation::subtract(ac, bd, precisionBits),
        approximation::add(ad, bc, precisionBits)
    };
}

ComplexInterval divide(
    const ComplexInterval& lhs,
    const ComplexInterval& rhs,
    std::size_t precisionBits) {
    // (a+bI)/(c+dI)
    //   = ((ac+bd) + (bc-ad)I) / (c^2+d^2)
    //
    // denominatorの区間が0を含む場合は、現precisionでは「分母が0でない」と
    // 証明できない。この関数は無理に巨大値を返さずdomain_errorとする。
    // 上位のcertified evaluatorは必要なら作業precisionを増やして再試行できる。
    const RealInterval cSquared = squareInterval(rhs.real(), precisionBits);
    const RealInterval dSquared = squareInterval(rhs.imaginary(), precisionBits);
    const RealInterval denominator = approximation::add(
        cSquared, dSquared, precisionBits);

    if (denominator.containsZero())
        throw std::domain_error(
            "ComplexInterval division denominator may contain zero");

    const RealInterval ac = approximation::multiply(
        lhs.real(), rhs.real(), precisionBits);
    const RealInterval bd = approximation::multiply(
        lhs.imaginary(), rhs.imaginary(), precisionBits);
    const RealInterval bc = approximation::multiply(
        lhs.imaginary(), rhs.real(), precisionBits);
    const RealInterval ad = approximation::multiply(
        lhs.real(), rhs.imaginary(), precisionBits);

    return ComplexInterval{
        approximation::divide(
            approximation::add(ac, bd, precisionBits),
            denominator,
            precisionBits),
        approximation::divide(
            approximation::subtract(bc, ad, precisionBits),
            denominator,
            precisionBits)
    };
}

} // namespace mmcal::approximation
