#pragma once

#include "big_int.hpp"

#include <compare>
#include <string>
#include <string_view>

namespace mmcal::numeric {

class Rational final {
public:
    Rational();
    explicit Rational(BigInt integer);
    Rational(BigInt numerator, BigInt denominator);

    [[nodiscard]] static Rational parse(std::string_view text, unsigned radix = 10);

    [[nodiscard]] const BigInt& numerator() const noexcept;
    [[nodiscard]] const BigInt& denominator() const noexcept;
    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] bool isInteger() const noexcept;
    [[nodiscard]] std::string toString(unsigned radix = 10) const;

    [[nodiscard]] Rational operator-() const;

    Rational& operator+=(const Rational& rhs);
    Rational& operator-=(const Rational& rhs);
    Rational& operator*=(const Rational& rhs);
    Rational& operator/=(const Rational& rhs);

    [[nodiscard]] std::strong_ordering operator<=>(const Rational& rhs) const;
    [[nodiscard]] bool operator==(const Rational& rhs) const noexcept;

private:
    BigInt numerator_;
    BigInt denominator_{1};

    void normalize();
    void normalizeSignAndZero();
};

[[nodiscard]] Rational operator+(Rational lhs, const Rational& rhs);
[[nodiscard]] Rational operator-(Rational lhs, const Rational& rhs);
[[nodiscard]] Rational operator*(Rational lhs, const Rational& rhs);
[[nodiscard]] Rational operator/(Rational lhs, const Rational& rhs);

} // namespace mmcal::numeric
