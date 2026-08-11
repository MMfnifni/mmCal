// 区間演算
#include "interval_math.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

namespace mmcal::approximation {
namespace {

[[nodiscard]] RealInterval zeroInterval(std::size_t precisionBits) {
    return RealInterval::fromRational(
        numeric::Rational{numeric::BigInt{0}}, precisionBits);
}

[[nodiscard]] numeric::BigFloat absoluteValue(const numeric::BigFloat& value) {
    return value.isNegative() ? -value : value;
}

} // namespace

RealInterval squareInterval(
    const RealInterval& value,
    std::size_t precisionBits) {
    const numeric::BigFloat zero;

    // x >= 0 なら x^2 は単調増加なので通常の区間積で十分。
    if (value.lower() >= zero)
        return multiply(value, value, precisionBits);

    // x <= 0 なら -x >= 0 へ反転してから二乗する。
    if (value.upper() <= zero) {
        const RealInterval positive = negate(value);
        return multiply(positive, positive, precisionBits);
    }

    // 0を跨ぐ場合、通常の interval multiply(value,value) では負の下限まで含んでしまうdependency problemがある。
    // しかし同一変数の平方は必ず非負。上限だけ両端の平方の大きい方を採れば、正しい値域 [0,max(a^2,b^2)] になる。
    const RealInterval lowerPoint = RealInterval::point(value.lower());
    const RealInterval upperPoint = RealInterval::point(value.upper());
    const RealInterval lowerSquared = multiply(lowerPoint, lowerPoint, precisionBits);
    const RealInterval upperSquared = multiply(upperPoint, upperPoint, precisionBits);
    const numeric::BigFloat upper = lowerSquared.upper() > upperSquared.upper()
        ? lowerSquared.upper()
        : upperSquared.upper();

    return RealInterval{zeroInterval(precisionBits).lower(), upper};
}

RealInterval absoluteInterval(
    const RealInterval& value,
    std::size_t precisionBits) {
    const numeric::BigFloat zero;

    if (value.lower() >= zero)
        return value.roundedOutward(precisionBits);
    if (value.upper() <= zero)
        return negate(value).roundedOutward(precisionBits);

    const numeric::BigFloat lowerMagnitude = absoluteValue(value.lower());
    const numeric::BigFloat upperMagnitude = absoluteValue(value.upper());
    const numeric::BigFloat upper = lowerMagnitude > upperMagnitude
        ? lowerMagnitude
        : upperMagnitude;
    return RealInterval{zeroInterval(precisionBits).lower(), upper};
}

} // namespace mmcal::approximation
