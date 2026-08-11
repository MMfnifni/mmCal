// 任意精度二進浮動小数
#include "big_float.hpp"
#include "detail/binary_scale.hpp"

#include <limits>
#include <stdexcept>
#include <utility>

namespace mmcal::numeric {
namespace {

void validatePrecision(std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("BigFloat precision must be at least one bit");
}

BigInt powerOfTwo(std::size_t exponent) {
    BigInt result{1};
    result <<= exponent;
    return result;
}

std::size_t toShiftCount(std::uint64_t value) {
    if (value > std::numeric_limits<std::size_t>::max())
        throw std::length_error("BigFloat shift distance is too large");
    return static_cast<std::size_t>(value);
}

// a >= b が既知のとき、符号付き64bitをオーバーフローさせずに a-b を求める。
std::uint64_t exponentDistance(
    BigFloat::exponent_type a,
    BigFloat::exponent_type b) noexcept {
    return static_cast<std::uint64_t>(a) - static_cast<std::uint64_t>(b);
}

BigFloat::exponent_type checkedAdd(
    BigFloat::exponent_type lhs,
    BigFloat::exponent_type rhs) {
    constexpr auto min = std::numeric_limits<BigFloat::exponent_type>::min();
    constexpr auto max = std::numeric_limits<BigFloat::exponent_type>::max();

    if (rhs > 0 && lhs > max - rhs)
        throw std::overflow_error("BigFloat exponent overflow");
    if (rhs < 0 && lhs < min - rhs)
        throw std::overflow_error("BigFloat exponent underflow");
    return static_cast<BigFloat::exponent_type>(lhs + rhs);
}

BigFloat::exponent_type checkedSubtract(
    BigFloat::exponent_type lhs,
    BigFloat::exponent_type rhs) {
    if (rhs == std::numeric_limits<BigFloat::exponent_type>::min()) {
        if (lhs >= 0)
            throw std::overflow_error("BigFloat exponent overflow");

        // lhs - INT64_MIN = lhs + 2^63。lhs<0なら結果は表現可能。
        const auto magnitude = static_cast<std::uint64_t>(-(lhs + 1)) + 1;
        const auto value = (std::uint64_t{1} << 63) - magnitude;
        return static_cast<BigFloat::exponent_type>(value);
    }
    return checkedAdd(lhs, static_cast<BigFloat::exponent_type>(-rhs));
}

BigFloat::exponent_type checkedAddShift(
    BigFloat::exponent_type exponent,
    std::size_t shift) {
    if (shift > static_cast<std::size_t>(
            std::numeric_limits<BigFloat::exponent_type>::max()))
        throw std::overflow_error("BigFloat exponent overflow");
    return checkedAdd(exponent, static_cast<BigFloat::exponent_type>(shift));
}

BigFloat::exponent_type checkedSubtractShift(
    BigFloat::exponent_type exponent,
    std::size_t shift) {
    constexpr auto max = std::numeric_limits<BigFloat::exponent_type>::max();
    if (shift > static_cast<std::size_t>(max))
        throw std::overflow_error("BigFloat exponent underflow");
    return checkedAdd(exponent, -static_cast<BigFloat::exponent_type>(shift));
}

bool shouldIncrement(
    const BigInt& quotient,
    const BigInt& remainder,
    const BigInt& divisor,
    bool negative,
    RoundingMode mode) {
    if (remainder.isZero())
        return false;

    switch (mode) {
    case RoundingMode::TowardZero:
        return false;
    case RoundingMode::TowardPositive:
        return !negative;
    case RoundingMode::TowardNegative:
        return negative;
    case RoundingMode::NearestEven:
        break;
    }

    const BigInt twiceRemainder = remainder * BigInt{2};
    const auto order = twiceRemainder <=> divisor;
    if (order == std::strong_ordering::greater)
        return true;
    if (order == std::strong_ordering::less)
        return false;

    // ちょうど中点なら、保持側の最下位bitが0になる方を選ぶ。
    return !quotient.isZero() && quotient.trailingZeroBits() == 0;
}

BigFloat fromPositiveRatio(
    const BigInt& numerator,
    const BigInt& denominator,
    bool negative,
    BigFloat::exponent_type exponentOffset,
    std::size_t precisionBits,
    RoundingMode roundingMode) {
    validatePrecision(precisionBits);

    if (numerator.isZero())
        return BigFloat::fromDyadic(BigInt{}, 0, precisionBits, roundingMode);
    if (denominator.isZero())
        throw std::domain_error("BigFloat division by zero");

    const auto ratioExponent = detail::floorLog2PositiveRatio(numerator, denominator);

    // m = (numerator / denominator) * 2^(p-1-ratioExponent)
    // とすると、floor(m) は原則p bitの仮数になる。
    // ここでは割り算の余りを捨てず、その余りから指定丸めを厳密に決定する。
    BigInt scaledNumerator = numerator;
    BigInt scaledDenominator = denominator;

    const std::size_t pMinusOne = precisionBits - 1;
    if (ratioExponent >= 0) {
        const auto ratioExp = static_cast<std::uint64_t>(ratioExponent);
        if (ratioExp <= pMinusOne)
            scaledNumerator <<= pMinusOne - static_cast<std::size_t>(ratioExp);
        else
            scaledDenominator <<= toShiftCount(ratioExp - pMinusOne);
    }
    else {
        const auto magnitude = std::uint64_t{0}
            - static_cast<std::uint64_t>(ratioExponent);
        if (magnitude > std::numeric_limits<std::size_t>::max() - pMinusOne)
            throw std::length_error("BigFloat conversion shift is too large");
        scaledNumerator <<= pMinusOne + static_cast<std::size_t>(magnitude);
    }

    auto quotientAndRemainder = divmod(scaledNumerator, scaledDenominator);
    BigInt quotient = std::move(quotientAndRemainder.quotient);
    const BigInt& remainder = quotientAndRemainder.remainder;

    if (shouldIncrement(
            quotient,
            remainder,
            scaledDenominator,
            negative,
            roundingMode))
        quotient += BigInt{1};

    if (negative)
        quotient = -quotient;

    auto outputExponent = checkedAdd(exponentOffset, ratioExponent);
    outputExponent = checkedSubtractShift(outputExponent, pMinusOne);
    return BigFloat::fromDyadic(
        std::move(quotient),
        outputExponent,
        precisionBits,
        RoundingMode::TowardZero);
}

std::strong_ordering compareMagnitude(
    const BigFloat& lhs,
    const BigFloat& rhs) {
    const BigInt lhsMagnitude = lhs.significand().abs();
    const BigInt rhsMagnitude = rhs.significand().abs();

    if (lhs.exponent() == rhs.exponent())
        return lhsMagnitude <=> rhsMagnitude;

    if (lhs.exponent() > rhs.exponent()) {
        const auto gap = exponentDistance(lhs.exponent(), rhs.exponent());
        if (gap >= rhsMagnitude.bitLength())
            return std::strong_ordering::greater;

        const BigInt scaledLhs = lhsMagnitude << toShiftCount(gap);
        return scaledLhs <=> rhsMagnitude;
    }

    const auto gap = exponentDistance(rhs.exponent(), lhs.exponent());
    if (gap >= lhsMagnitude.bitLength())
        return std::strong_ordering::less;

    const BigInt scaledRhs = rhsMagnitude << toShiftCount(gap);
    return lhsMagnitude <=> scaledRhs;
}

} // namespace

BigFloat::BigFloat() = default;

BigFloat::BigFloat(
    BigInt significand,
    exponent_type exponent,
    std::size_t precisionBits)
    : significand_(std::move(significand)),
      exponent_(exponent),
      precisionBits_(precisionBits) {
    validatePrecision(precisionBits_);
    normalize();
}

BigFloat BigFloat::fromBigInt(
    const BigInt& value,
    std::size_t precisionBits,
    RoundingMode roundingMode) {
    return fromDyadic(value, 0, precisionBits, roundingMode);
}

BigFloat BigFloat::fromRational(
    const Rational& value,
    std::size_t precisionBits,
    RoundingMode roundingMode) {
    const bool negative = value.numerator().isNegative();
    return fromPositiveRatio(
        value.numerator().abs(),
        value.denominator(),
        negative,
        0,
        precisionBits,
        roundingMode);
}

BigFloat BigFloat::fromDyadic(
    BigInt exactSignificand,
    exponent_type exponent,
    std::size_t precisionBits,
    RoundingMode roundingMode) {
    validatePrecision(precisionBits);

    if (exactSignificand.isZero())
        return BigFloat{BigInt{}, 0, precisionBits};

    const bool negative = exactSignificand.isNegative();
    BigInt magnitude = exactSignificand.abs();
    const std::size_t bitLength = magnitude.bitLength();

    if (bitLength <= precisionBits)
        return BigFloat{std::move(exactSignificand), exponent, precisionBits};

    const std::size_t discardedBits = bitLength - precisionBits;
    BigInt quotient = magnitude >> discardedBits;
    const BigInt retainedPart = quotient << discardedBits;
    const BigInt remainder = magnitude - retainedPart;
    const BigInt divisor = powerOfTwo(discardedBits);

    if (shouldIncrement(
            quotient,
            remainder,
            divisor,
            negative,
            roundingMode))
        quotient += BigInt{1};

    if (negative)
        quotient = -quotient;

    exponent = checkedAddShift(exponent, discardedBits);
    return BigFloat{std::move(quotient), exponent, precisionBits};
}

bool BigFloat::isZero() const noexcept {
    return significand_.isZero();
}

bool BigFloat::isNegative() const noexcept {
    return significand_.isNegative();
}

bool BigFloat::isPositive() const noexcept {
    return significand_.isPositive();
}

std::size_t BigFloat::precisionBits() const noexcept {
    return precisionBits_;
}

const BigInt& BigFloat::significand() const noexcept {
    return significand_;
}

BigFloat::exponent_type BigFloat::exponent() const noexcept {
    return exponent_;
}

BigFloat BigFloat::operator-() const {
    return BigFloat{-significand_, exponent_, precisionBits_};
}

BigFloat BigFloat::rounded(
    std::size_t precisionBits,
    RoundingMode roundingMode) const {
    return fromDyadic(significand_, exponent_, precisionBits, roundingMode);
}

Rational BigFloat::toRational() const {
    if (isZero())
        return Rational{};

    if (exponent_ >= 0) {
        const auto shift = toShiftCount(static_cast<std::uint64_t>(exponent_));
        return Rational{significand_ << shift};
    }

    const auto shift = toShiftCount(
        std::uint64_t{0} - static_cast<std::uint64_t>(exponent_));
    return Rational{significand_, powerOfTwo(shift)};
}

std::string BigFloat::toString() const {
    if (isZero())
        return "0";
    if (exponent_ == 0)
        return significand_.toString();
    return significand_.toString() + " * 2^" + std::to_string(exponent_);
}

std::strong_ordering BigFloat::operator<=>(const BigFloat& rhs) const {
    if (isZero() && rhs.isZero())
        return std::strong_ordering::equal;
    if (isZero())
        return rhs.isNegative()
            ? std::strong_ordering::greater
            : std::strong_ordering::less;
    if (rhs.isZero())
        return isNegative()
            ? std::strong_ordering::less
            : std::strong_ordering::greater;
    if (isNegative() != rhs.isNegative())
        return isNegative()
            ? std::strong_ordering::less
            : std::strong_ordering::greater;

    const auto magnitudeOrder = compareMagnitude(*this, rhs);
    if (!isNegative())
        return magnitudeOrder;

    if (magnitudeOrder == std::strong_ordering::less)
        return std::strong_ordering::greater;
    if (magnitudeOrder == std::strong_ordering::greater)
        return std::strong_ordering::less;
    return std::strong_ordering::equal;
}

bool BigFloat::operator==(const BigFloat& rhs) const {
    return (*this <=> rhs) == std::strong_ordering::equal;
}

void BigFloat::normalize() {
    if (significand_.isZero()) {
        exponent_ = 0;
        return;
    }

    // 同じdyadic値に複数の表現を作らないため、仮数から2の因子を全て除く。
    // 例: 12 * 2^-3 -> 3 * 2^-1。
    const std::size_t zeros = significand_.trailingZeroBits();
    if (zeros == 0)
        return;

    significand_ >>= zeros;
    exponent_ = checkedAddShift(exponent_, zeros);
}

BigFloat add(
    const BigFloat& lhs,
    const BigFloat& rhs,
    std::size_t precisionBits,
    RoundingMode roundingMode) {
    validatePrecision(precisionBits);

    if (lhs.isZero())
        return rhs.rounded(precisionBits, roundingMode);
    if (rhs.isZero())
        return lhs.rounded(precisionBits, roundingMode);

    const auto commonExponent = lhs.exponent() < rhs.exponent()
        ? lhs.exponent()
        : rhs.exponent();

    BigInt lhsSignificand = lhs.significand();
    BigInt rhsSignificand = rhs.significand();

    if (lhs.exponent() > commonExponent)
        lhsSignificand <<= toShiftCount(
            exponentDistance(lhs.exponent(), commonExponent));
    if (rhs.exponent() > commonExponent)
        rhsSignificand <<= toShiftCount(
            exponentDistance(rhs.exponent(), commonExponent));

    return BigFloat::fromDyadic(
        lhsSignificand + rhsSignificand,
        commonExponent,
        precisionBits,
        roundingMode);
}

BigFloat subtract(
    const BigFloat& lhs,
    const BigFloat& rhs,
    std::size_t precisionBits,
    RoundingMode roundingMode) {
    return add(lhs, -rhs, precisionBits, roundingMode);
}

BigFloat multiply(
    const BigFloat& lhs,
    const BigFloat& rhs,
    std::size_t precisionBits,
    RoundingMode roundingMode) {
    validatePrecision(precisionBits);

    if (lhs.isZero() || rhs.isZero())
        return BigFloat::fromDyadic(BigInt{}, 0, precisionBits, roundingMode);

    return BigFloat::fromDyadic(
        lhs.significand() * rhs.significand(),
        checkedAdd(lhs.exponent(), rhs.exponent()),
        precisionBits,
        roundingMode);
}

BigFloat divide(
    const BigFloat& lhs,
    const BigFloat& rhs,
    std::size_t precisionBits,
    RoundingMode roundingMode) {
    validatePrecision(precisionBits);
    if (rhs.isZero())
        throw std::domain_error("BigFloat division by zero");
    if (lhs.isZero())
        return BigFloat::fromDyadic(BigInt{}, 0, precisionBits, roundingMode);

    const bool negative = lhs.isNegative() != rhs.isNegative();
    return fromPositiveRatio(
        lhs.significand().abs(),
        rhs.significand().abs(),
        negative,
        checkedSubtract(lhs.exponent(), rhs.exponent()),
        precisionBits,
        roundingMode);
}

} // namespace mmcal::numeric
