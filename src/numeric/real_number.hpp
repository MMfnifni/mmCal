#pragma once

#include "big_int.hpp"
#include "rational.hpp"

#include <compare>
#include <string>
#include <variant>

namespace mmcal::numeric {

// 現段階で厳密に表せる実数。整数化できる分数は常にBigIntへ縮約する。
class RealNumber final {
public:
    RealNumber() = default;
    RealNumber(BigInt integer);
    RealNumber(Rational rational);

    [[nodiscard]] bool isInteger() const noexcept;
    [[nodiscard]] bool isRational() const noexcept;
    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] bool isNegative() const noexcept;

    [[nodiscard]] const BigInt& asInteger() const;
    [[nodiscard]] const Rational& asRational() const;
    [[nodiscard]] Rational toRational() const;
    [[nodiscard]] RealNumber abs() const;
    [[nodiscard]] std::string toString(unsigned radix = 10) const;

    [[nodiscard]] RealNumber operator-() const;

    RealNumber& operator+=(const RealNumber& rhs);
    RealNumber& operator-=(const RealNumber& rhs);
    RealNumber& operator*=(const RealNumber& rhs);
    RealNumber& operator/=(const RealNumber& rhs);

    [[nodiscard]] std::strong_ordering operator<=>(const RealNumber& rhs) const;
    [[nodiscard]] bool operator==(const RealNumber& rhs) const;

private:
    std::variant<BigInt, Rational> value_;

    void normalize();
};

[[nodiscard]] RealNumber operator+(RealNumber lhs, const RealNumber& rhs);
[[nodiscard]] RealNumber operator-(RealNumber lhs, const RealNumber& rhs);
[[nodiscard]] RealNumber operator*(RealNumber lhs, const RealNumber& rhs);
[[nodiscard]] RealNumber operator/(RealNumber lhs, const RealNumber& rhs);

} // namespace mmcal::numeric
