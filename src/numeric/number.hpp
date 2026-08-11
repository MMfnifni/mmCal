#pragma once

#include "real_number.hpp"

#include <cstdint>
#include <string>
#include <variant>

namespace mmcal::numeric {

struct ComplexNumber final {
    RealNumber real;
    RealNumber imaginary;

    [[nodiscard]] bool operator==(const ComplexNumber&) const = default;
};

// 厳密実数と厳密複素数を統一して扱う数値型。
class Number final {
public:
    Number() = default;
    Number(BigInt integer);
    Number(Rational rational);
    Number(RealNumber real);

    [[nodiscard]] static Number complex(RealNumber real, RealNumber imaginary);

    [[nodiscard]] bool isReal() const noexcept;
    [[nodiscard]] bool isComplex() const noexcept;
    [[nodiscard]] bool isZero() const noexcept;

    [[nodiscard]] const RealNumber& asReal() const;
    [[nodiscard]] const ComplexNumber& asComplex() const;
    [[nodiscard]] RealNumber realPart() const;
    [[nodiscard]] RealNumber imaginaryPart() const;
    [[nodiscard]] Number conjugate() const;
    [[nodiscard]] std::string toString(unsigned radix = 10) const;

    [[nodiscard]] Number operator-() const;

    Number& operator+=(const Number& rhs);
    Number& operator-=(const Number& rhs);
    Number& operator*=(const Number& rhs);
    Number& operator/=(const Number& rhs);

    [[nodiscard]] bool operator==(const Number& rhs) const;

private:
    std::variant<RealNumber, ComplexNumber> value_;

    explicit Number(ComplexNumber complex);
    void normalize();
};

[[nodiscard]] Number operator+(Number lhs, const Number& rhs);
[[nodiscard]] Number operator-(Number lhs, const Number& rhs);
[[nodiscard]] Number operator*(Number lhs, const Number& rhs);
[[nodiscard]] Number operator/(Number lhs, const Number& rhs);
[[nodiscard]] Number integerPower(Number base, std::uint64_t exponent);

} // namespace mmcal::numeric
