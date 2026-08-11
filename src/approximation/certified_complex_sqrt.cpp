#include "certified_complex_sqrt.hpp"

#include "certified_sqrt.hpp"
#include "interval_math.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <stdexcept>

namespace mmcal::approximation {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

[[nodiscard]] RealInterval zeroInterval(std::size_t precisionBits) {
    return RealInterval::fromRational(rational(0), precisionBits);
}

[[nodiscard]] RealInterval twoInterval(std::size_t precisionBits) {
    return RealInterval::fromRational(rational(2), precisionBits);
}

[[nodiscard]] bool isExactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

[[nodiscard]] RealInterval clampNonNegative(
    const RealInterval& value,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (value.upper() < zero) {
        // この関数へ渡す量は数学的には必ず非負。
        // 上端まで負なら、包含保証をどこかで破った内部バグを意味する。
        throw std::logic_error(
            "Certified complex sqrt produced an impossible negative radicand interval");
    }

    if (value.lower() >= zero)
        return value;

    return RealInterval{zeroInterval(precisionBits).lower(), value.upper()};
}

} // namespace

ComplexInterval enclosePrincipalComplexSqrt(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument(
            "Certified complex sqrt precision must be at least one bit");

    // z = x + yI, r = |z| = sqrt(x^2+y^2) とするとprincipal sqrtは
    //
    //   u = sqrt((r+x)/2) >= 0
    //   |v| = sqrt((r-x)/2)
    //
    // で表せる。y != 0ならvの符号はyと同じ。
    // 各演算はRealIntervalの外向き丸めなので、途中の丸め誤差もすべて包含される。
    const RealInterval xSquared = squareInterval(value.real(), precisionBits);
    const RealInterval ySquared = squareInterval(value.imaginary(), precisionBits);
    const RealInterval magnitudeSquared = add(xSquared, ySquared, precisionBits);
    const RealInterval magnitude = encloseSqrt(magnitudeSquared, precisionBits).interval;

    const RealInterval two = twoInterval(precisionBits);
    const RealInterval uSquared = clampNonNegative(
        divide(add(magnitude, value.real(), precisionBits), two, precisionBits),
        precisionBits);
    const RealInterval vSquared = clampNonNegative(
        divide(subtract(magnitude, value.real(), precisionBits), two, precisionBits),
        precisionBits);

    const RealInterval u = encloseSqrt(uSquared, precisionBits).interval;
    const RealInterval vMagnitude = encloseSqrt(vSquared, precisionBits).interval;
    const BigFloat zero;

    // Im(z)の符号が証明できるならprincipal sqrtの虚部符号も一意。
    if (value.imaginary().lower() > zero)
        return ComplexInterval{u, vMagnitude};
    if (value.imaginary().upper() < zero)
        return ComplexInterval{u, negate(vMagnitude)};

    // Im(z)が厳密に0なら負実軸上だけ+側の平方根を採る。
    if (isExactZero(value.imaginary())) {
        if (value.real().lower() >= zero)
            return ComplexInterval{u, zeroInterval(precisionBits)};
        if (value.real().upper() < zero)
            return ComplexInterval{u, vMagnitude};

        // x自体の符号がまだ確定しない区間。真値が正なら虚部0、負なら+sqrt(-x)。
        // [0, upper] とすれば両方を包含できる。
        return ComplexInterval{
            u,
            RealInterval{zeroInterval(precisionBits).lower(), vMagnitude.upper()}
        };
    }

    // 入力長方形が実軸を跨ぐ場合。特に負実軸（principal sqrtのbranch cut）を跨ぐと、上半平面からは +v、下半平面からは -v へ近づく。
    // どちらが真値か現precisionで証明できない以上、両方を含む対称区間へ広げる。
    return ComplexInterval{
        u,
        RealInterval{-vMagnitude.upper(), vMagnitude.upper()}
    };
}

} // namespace mmcal::approximation
