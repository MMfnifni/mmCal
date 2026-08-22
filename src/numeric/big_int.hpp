#pragma once

#include "detail/big_uint.hpp"

#include <compare>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace mmcal::numeric {

class BigInt;
struct BigIntDivModResult;
struct IntegerSqrtResult;

[[nodiscard]] BigIntDivModResult divmod(const BigInt& dividend, const BigInt& divisor);
[[nodiscard]] IntegerSqrtResult integerSqrt(const BigInt& value);

class BigInt final {
public:
    BigInt() = default;
    explicit BigInt(std::int64_t value);

    [[nodiscard]] static BigInt fromUnsigned(std::uint64_t value);
    [[nodiscard]] static BigInt parse(std::string_view text, unsigned radix = 10);

    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] bool isNegative() const noexcept;
    [[nodiscard]] bool isPositive() const noexcept;
    [[nodiscard]] std::size_t bitLength() const noexcept;
    [[nodiscard]] std::size_t trailingZeroBits() const noexcept;
    [[nodiscard]] std::size_t populationCount() const noexcept;
    [[nodiscard]] bool testBit(std::size_t index) const;
    [[nodiscard]] std::uint32_t modulo(std::uint32_t divisor) const;
    [[nodiscard]] BigInt abs() const;
    [[nodiscard]] std::string toString(unsigned radix = 10) const;

    [[nodiscard]] BigInt operator-() const;

    BigInt& operator+=(const BigInt& rhs);
    BigInt& operator-=(const BigInt& rhs);
    BigInt& operator*=(const BigInt& rhs);
    BigInt& operator/=(const BigInt& rhs);
    BigInt& operator%=(const BigInt& rhs);
    BigInt& operator&=(const BigInt& rhs);
    BigInt& operator|=(const BigInt& rhs);
    BigInt& operator^=(const BigInt& rhs);
    BigInt& operator<<=(std::size_t bits);
    BigInt& operator>>=(std::size_t bits);

    [[nodiscard]] std::strong_ordering operator<=>(const BigInt& rhs) const noexcept;
    [[nodiscard]] bool operator==(const BigInt& rhs) const noexcept;

private:
    bool negative_ = false;
    detail::BigUInt magnitude_;

    explicit BigInt(detail::BigUInt magnitude, bool negative) noexcept;
    void normalizeSign() noexcept;

    friend BigIntDivModResult divmod(const BigInt& dividend, const BigInt& divisor);
    friend IntegerSqrtResult integerSqrt(const BigInt& value);
};

struct BigIntDivModResult final {
    BigInt quotient;
    BigInt remainder;
};

[[nodiscard]] BigInt operator+(BigInt lhs, const BigInt& rhs);
[[nodiscard]] BigInt operator-(BigInt lhs, const BigInt& rhs);
[[nodiscard]] BigInt operator*(BigInt lhs, const BigInt& rhs);
[[nodiscard]] BigInt operator/(BigInt lhs, const BigInt& rhs);
[[nodiscard]] BigInt operator%(BigInt lhs, const BigInt& rhs);
[[nodiscard]] BigInt operator&(BigInt lhs, const BigInt& rhs);
[[nodiscard]] BigInt operator|(BigInt lhs, const BigInt& rhs);
[[nodiscard]] BigInt operator^(BigInt lhs, const BigInt& rhs);
[[nodiscard]] BigInt operator~(const BigInt& value);
[[nodiscard]] BigInt operator<<(BigInt value, std::size_t bits);
[[nodiscard]] BigInt operator>>(BigInt value, std::size_t bits);

} // namespace mmcal::numeric
