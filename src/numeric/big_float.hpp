#pragma once

#include "big_int.hpp"
#include "rational.hpp"
#include "rounding_mode.hpp"

#include <compare>
#include <cstddef>
#include <cstdint>
#include <string>

namespace mmcal::numeric {

// 有限な2進任意精度浮動小数点数。
//
// 値は常に
//     significand * 2^exponent
// として保持する。significand は符号付きBigIntで、0以外では末尾の2の因子を取り除いた奇数に正規化する。
// precisionBits は、この値を生成したときの目標有効ビット数を表す。
//
// NaN / ±Infinity はまだ持たない。数学的domain errorは上位層で明示的に扱う。
class BigFloat final {
public:
    using exponent_type = std::int64_t;

    BigFloat();

    [[nodiscard]] static BigFloat fromBigInt(
        const BigInt& value,
        std::size_t precisionBits,
        RoundingMode roundingMode = RoundingMode::NearestEven);

    [[nodiscard]] static BigFloat fromRational(
        const Rational& value,
        std::size_t precisionBits,
        RoundingMode roundingMode = RoundingMode::NearestEven);

    // exactSignificand * 2^exponent を precisionBits へ丸めて構築する。
    [[nodiscard]] static BigFloat fromDyadic(
        BigInt exactSignificand,
        exponent_type exponent,
        std::size_t precisionBits,
        RoundingMode roundingMode = RoundingMode::NearestEven);

    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] bool isNegative() const noexcept;
    [[nodiscard]] bool isPositive() const noexcept;
    [[nodiscard]] std::size_t precisionBits() const noexcept;
    [[nodiscard]] const BigInt& significand() const noexcept;
    [[nodiscard]] exponent_type exponent() const noexcept;

    [[nodiscard]] BigFloat operator-() const;
    [[nodiscard]] BigFloat rounded(
        std::size_t precisionBits,
        RoundingMode roundingMode = RoundingMode::NearestEven) const;

    // 現在のdyadic値を厳密なRationalへ戻す。表示や検証用にも使える。
    [[nodiscard]] Rational toRational() const;
    [[nodiscard]] std::string toString() const;

    [[nodiscard]] std::strong_ordering operator<=>(const BigFloat& rhs) const;
    [[nodiscard]] bool operator==(const BigFloat& rhs) const;

private:
    BigInt significand_;
    exponent_type exponent_ = 0;
    std::size_t precisionBits_ = 1;

    BigFloat(BigInt significand, exponent_type exponent, std::size_t precisionBits);
    void normalize();
};

// 四則演算は丸め方向を暗黙にしない。
// RealIntervalで下端/上端を別方向へ丸めるため、必ず呼出側が指定する。
[[nodiscard]] BigFloat add(
    const BigFloat& lhs,
    const BigFloat& rhs,
    std::size_t precisionBits,
    RoundingMode roundingMode);

[[nodiscard]] BigFloat subtract(
    const BigFloat& lhs,
    const BigFloat& rhs,
    std::size_t precisionBits,
    RoundingMode roundingMode);

[[nodiscard]] BigFloat multiply(
    const BigFloat& lhs,
    const BigFloat& rhs,
    std::size_t precisionBits,
    RoundingMode roundingMode);

[[nodiscard]] BigFloat divide(
    const BigFloat& lhs,
    const BigFloat& rhs,
    std::size_t precisionBits,
    RoundingMode roundingMode);

} // namespace mmcal::numeric
