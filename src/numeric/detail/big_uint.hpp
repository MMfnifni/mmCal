#pragma once

#include <compare>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace mmcal::numeric::detail {

class BigUInt;
struct BigUIntDivModResult;

class BigUInt final {
public:
    using limb_type = std::uint32_t;
    using double_limb_type = std::uint64_t;

    BigUInt() = default;
    explicit BigUInt(std::uint64_t value);

    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] std::size_t limbCount() const noexcept;
    [[nodiscard]] std::size_t bitLength() const noexcept;
    [[nodiscard]] std::size_t trailingZeroBits() const noexcept;
    [[nodiscard]] std::size_t populationCount() const noexcept;
    [[nodiscard]] bool testBit(std::size_t index) const noexcept;

    [[nodiscard]] static BigUInt parse(std::string_view text, unsigned radix = 10);
    [[nodiscard]] std::string toString(unsigned radix = 10) const;

    BigUInt& operator+=(const BigUInt& rhs);
    BigUInt& operator-=(const BigUInt& rhs);
    BigUInt& operator*=(const BigUInt& rhs);
    BigUInt& operator/=(const BigUInt& rhs);
    BigUInt& operator%=(const BigUInt& rhs);
    BigUInt& operator&=(const BigUInt& rhs);
    BigUInt& operator|=(const BigUInt& rhs);
    BigUInt& operator^=(const BigUInt& rhs);

    // 任意ビット数だけシフトする。右シフトでは下位ビットを切り捨てる。
    BigUInt& operator<<=(std::size_t bits);
    BigUInt& operator>>=(std::size_t bits);

    [[nodiscard]] std::strong_ordering operator<=>(const BigUInt& rhs) const noexcept;
    [[nodiscard]] bool operator==(const BigUInt& rhs) const noexcept;

private:
    // 内部基数は 2^32。下位limbから並べるリトルエンディアン形式で保持する。
    // limbs_[0] が最下位32bitを保持する。
    std::vector<limb_type> limbs_;

    void normalize() noexcept;
    void addSmall(limb_type value);
    void multiplySmall(limb_type value);
    [[nodiscard]] limb_type divideSmall(limb_type divisor);

    [[nodiscard]] static unsigned digitValue(char ch) noexcept;
    [[nodiscard]] static char digitCharacter(unsigned value) noexcept;

    // 巨大除算の内部block操作。公開APIへlimb表現を漏らさず、
    // Burnikel-Ziegler再帰だけがlittle-endian limb列を直接扱う。
    [[nodiscard]] static BigUInt limbSlice(
        const BigUInt& value, std::size_t start, std::size_t count);
    [[nodiscard]] static BigUInt concatenateLimbs(
        const BigUInt& high, const BigUInt& low, std::size_t lowWidth);
    [[nodiscard]] static BigUInt allOnes(std::size_t limbCount);
    [[nodiscard]] static BigUIntDivModResult divideKnuth(
        const BigUInt& dividend, const BigUInt& divisor);
    [[nodiscard]] static BigUIntDivModResult divide2n1n(
        const BigUInt& dividend, const BigUInt& divisor);
    [[nodiscard]] static BigUIntDivModResult divide3n2n(
        const BigUInt& dividend, const BigUInt& divisor);
    [[nodiscard]] static BigUIntDivModResult divideBurnikelZiegler(
        const BigUInt& dividend, const BigUInt& divisor);

    friend BigUIntDivModResult divmod(const BigUInt& dividend, const BigUInt& divisor);
};

struct BigUIntDivModResult final {
    BigUInt quotient;
    BigUInt remainder;
};

[[nodiscard]] BigUIntDivModResult divmod(
    const BigUInt& dividend,
    const BigUInt& divisor);

[[nodiscard]] BigUInt operator+(BigUInt lhs, const BigUInt& rhs);
[[nodiscard]] BigUInt operator-(BigUInt lhs, const BigUInt& rhs);
[[nodiscard]] BigUInt operator*(BigUInt lhs, const BigUInt& rhs);
[[nodiscard]] BigUInt operator/(BigUInt lhs, const BigUInt& rhs);
[[nodiscard]] BigUInt operator%(BigUInt lhs, const BigUInt& rhs);
[[nodiscard]] BigUInt operator&(BigUInt lhs, const BigUInt& rhs);
[[nodiscard]] BigUInt operator|(BigUInt lhs, const BigUInt& rhs);
[[nodiscard]] BigUInt operator^(BigUInt lhs, const BigUInt& rhs);
[[nodiscard]] BigUInt operator<<(BigUInt value, std::size_t bits);
[[nodiscard]] BigUInt operator>>(BigUInt value, std::size_t bits);

} // namespace mmcal::numeric::detail
