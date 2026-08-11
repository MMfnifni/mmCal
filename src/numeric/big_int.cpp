// 任意精度整数BigInt
#include "big_int.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::numeric {

BigInt::BigInt(std::int64_t value) {
    if (value >= 0) {
        magnitude_ = detail::BigUInt{static_cast<std::uint64_t>(value)};
        return;
    }

    // INT64_MIN を直接反転すると int64_t がオーバーフローするため、1ずらして絶対値を求める。
    const auto magnitude = static_cast<std::uint64_t>(-(value + 1)) + 1;
    magnitude_ = detail::BigUInt{magnitude};
    negative_ = true;
}

BigInt::BigInt(detail::BigUInt magnitude, bool negative) noexcept
    : negative_(negative), magnitude_(std::move(magnitude)) {
    normalizeSign();
}

BigInt BigInt::parse(std::string_view text, unsigned radix) {
    if (text.empty())
        throw std::invalid_argument("BigInt cannot parse an empty string");

    bool negative = false;
    std::size_t position = 0;

    if (text.front() == '+' || text.front() == '-') {
        negative = text.front() == '-';
        position = 1;
    }

    if (position == text.size())
        throw std::invalid_argument("BigInt requires at least one digit");

    auto magnitude = detail::BigUInt::parse(text.substr(position), radix);
    return BigInt{std::move(magnitude), negative};
}

bool BigInt::isZero() const noexcept {
    return magnitude_.isZero();
}

bool BigInt::isNegative() const noexcept {
    return negative_;
}

bool BigInt::isPositive() const noexcept {
    return !negative_ && !isZero();
}

std::size_t BigInt::bitLength() const noexcept {
    return magnitude_.bitLength();
}

std::size_t BigInt::trailingZeroBits() const noexcept {
    return magnitude_.trailingZeroBits();
}

BigInt BigInt::abs() const {
    return BigInt{magnitude_, false};
}

std::string BigInt::toString(unsigned radix) const {
    auto result = magnitude_.toString(radix);
    if (negative_)
        result.insert(result.begin(), '-');
    return result;
}

BigInt BigInt::operator-() const {
    BigInt result = *this;
    if (!result.isZero())
        result.negative_ = !result.negative_;
    return result;
}

BigInt& BigInt::operator+=(const BigInt& rhs) {
    if (negative_ == rhs.negative_) {
        magnitude_ += rhs.magnitude_;
        return *this;
    }

    if (magnitude_ == rhs.magnitude_) {
        magnitude_ = detail::BigUInt{};
        negative_ = false;
        return *this;
    }

    if (magnitude_ > rhs.magnitude_)
        magnitude_ -= rhs.magnitude_;
    else {
        auto resultMagnitude = rhs.magnitude_ - magnitude_;
        magnitude_ = std::move(resultMagnitude);
        negative_ = rhs.negative_;
    }

    normalizeSign();
    return *this;
}

BigInt& BigInt::operator-=(const BigInt& rhs) {
    return *this += -rhs;
}

BigInt& BigInt::operator*=(const BigInt& rhs) {
    const bool resultNegative = negative_ != rhs.negative_;
    magnitude_ *= rhs.magnitude_;
    negative_ = resultNegative;
    normalizeSign();
    return *this;
}

BigInt& BigInt::operator/=(const BigInt& rhs) {
    auto result = divmod(*this, rhs);
    *this = std::move(result.quotient);
    return *this;
}

BigInt& BigInt::operator%=(const BigInt& rhs) {
    auto result = divmod(*this, rhs);
    *this = std::move(result.remainder);
    return *this;
}

BigInt& BigInt::operator<<=(std::size_t bits) {
    magnitude_ <<= bits;
    return *this;
}

BigInt& BigInt::operator>>=(std::size_t bits) {
    // 符号と絶対値を分離しているため、右シフトは絶対値を切り捨てる。
    // 負数でも算術シフトではなく 0 方向への切り捨てになる。
    magnitude_ >>= bits;
    normalizeSign();
    return *this;
}

std::strong_ordering BigInt::operator<=>(const BigInt& rhs) const noexcept {
    if (negative_ != rhs.negative_)
        return negative_ ? std::strong_ordering::less : std::strong_ordering::greater;

    const auto magnitudeOrder = magnitude_ <=> rhs.magnitude_;
    if (!negative_)
        return magnitudeOrder;

    if (magnitudeOrder == std::strong_ordering::less)
        return std::strong_ordering::greater;
    if (magnitudeOrder == std::strong_ordering::greater)
        return std::strong_ordering::less;
    return std::strong_ordering::equal;
}

bool BigInt::operator==(const BigInt& rhs) const noexcept {
    return negative_ == rhs.negative_ && magnitude_ == rhs.magnitude_;
}

void BigInt::normalizeSign() noexcept {
    if (magnitude_.isZero())
        negative_ = false;
}

BigIntDivModResult divmod(const BigInt& dividend, const BigInt& divisor) {
    if (divisor.isZero())
        throw std::domain_error("BigInt division by zero");

    auto unsignedResult = detail::divmod(dividend.magnitude_, divisor.magnitude_);

    BigInt quotient{
        std::move(unsignedResult.quotient),
        dividend.negative_ != divisor.negative_};
    BigInt remainder{
        std::move(unsignedResult.remainder),
        dividend.negative_};

    return {std::move(quotient), std::move(remainder)};
}

BigInt operator+(BigInt lhs, const BigInt& rhs) {
    lhs += rhs;
    return lhs;
}

BigInt operator-(BigInt lhs, const BigInt& rhs) {
    lhs -= rhs;
    return lhs;
}

BigInt operator*(BigInt lhs, const BigInt& rhs) {
    lhs *= rhs;
    return lhs;
}

BigInt operator/(BigInt lhs, const BigInt& rhs) {
    lhs /= rhs;
    return lhs;
}

BigInt operator%(BigInt lhs, const BigInt& rhs) {
    lhs %= rhs;
    return lhs;
}

BigInt operator<<(BigInt value, std::size_t bits) {
    value <<= bits;
    return value;
}

BigInt operator>>(BigInt value, std::size_t bits) {
    value >>= bits;
    return value;
}

} // namespace mmcal::numeric
