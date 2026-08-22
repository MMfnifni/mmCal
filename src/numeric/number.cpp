// 実数・複素数を統合するNumber
#include "number.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::numeric {
namespace {

[[nodiscard]] std::string formatImaginaryTerm(
    const RealNumber& magnitude,
    unsigned radix) {
    if (magnitude.isInteger()) {
        if (magnitude == RealNumber{BigInt{1}})
            return "I";
        return magnitude.toString(radix) + "I";
    }

    // a/b I は a/(b I) とも読めるため，複素有理数は aI/b と表示する。
    // parser上も (a*I)/b となり，元のexact値へ再parseできる。
    const Rational& value = magnitude.asRational();
    std::string numerator = value.numerator() == BigInt{1}
        ? std::string{"I"}
        : value.numerator().toString(radix) + "I";
    return numerator + "/" + value.denominator().toString(radix);
}

} // namespace

Number::Number(BigInt integer)
    : value_(RealNumber{std::move(integer)}) {}

Number::Number(Rational rational)
    : value_(RealNumber{std::move(rational)}) {}

Number::Number(RealNumber real)
    : value_(std::move(real)) {}

Number::Number(ComplexNumber complex)
    : value_(std::move(complex)) {
    normalize();
}

Number Number::complex(RealNumber real, RealNumber imaginary) {
    return Number{ComplexNumber{std::move(real), std::move(imaginary)}};
}

bool Number::isReal() const noexcept {
    return std::holds_alternative<RealNumber>(value_);
}

bool Number::isComplex() const noexcept {
    return std::holds_alternative<ComplexNumber>(value_);
}

bool Number::isZero() const noexcept {
    if (const auto* real = std::get_if<RealNumber>(&value_))
        return real->isZero();

    const auto& complex = std::get<ComplexNumber>(value_);
    return complex.real.isZero() && complex.imaginary.isZero();
}

const RealNumber& Number::asReal() const {
    if (!isReal())
        throw std::logic_error("Number does not contain a real value");

    return std::get<RealNumber>(value_);
}

const ComplexNumber& Number::asComplex() const {
    if (!isComplex())
        throw std::logic_error("Number does not contain a complex value");

    return std::get<ComplexNumber>(value_);
}

RealNumber Number::realPart() const {
    if (const auto* real = std::get_if<RealNumber>(&value_))
        return *real;

    return std::get<ComplexNumber>(value_).real;
}

RealNumber Number::imaginaryPart() const {
    if (isReal())
        return RealNumber{};

    return std::get<ComplexNumber>(value_).imaginary;
}

Number Number::conjugate() const {
    if (isReal())
        return *this;
    const auto& value = asComplex();
    return complex(value.real, -value.imaginary);
}

std::string Number::toString(unsigned radix) const {
    if (const auto* real = std::get_if<RealNumber>(&value_))
        return real->toString(radix);

    const auto& value = std::get<ComplexNumber>(value_);
    const bool hasReal = !value.real.isZero();
    const bool negativeImaginary = value.imaginary.isNegative();
    const RealNumber magnitude = value.imaginary.abs();
    const std::string imaginary = formatImaginaryTerm(magnitude, radix);

    if (!hasReal)
        return negativeImaginary ? "-" + imaginary : imaginary;

    return value.real.toString(radix)
        + (negativeImaginary ? "-" : "+")
        + imaginary;
}

Number Number::operator-() const {
    if (isReal())
        return Number{-asReal()};
    const auto& value = asComplex();
    return complex(-value.real, -value.imaginary);
}

Number& Number::operator+=(const Number& rhs) {
    if (isReal() && rhs.isReal()) {
        std::get<RealNumber>(value_) += rhs.asReal();
        return *this;
    }

    if (auto* lhs = std::get_if<ComplexNumber>(&value_)) {
        if (rhs.isReal()) {
            lhs->real += rhs.asReal();
            normalize();
            return *this;
        }
        const auto& other = rhs.asComplex();
        lhs->real += other.real;
        lhs->imaginary += other.imaginary;
        normalize();
        return *this;
    }

    const RealNumber lhsReal = asReal();
    const auto& other = rhs.asComplex();
    value_ = ComplexNumber{lhsReal + other.real, other.imaginary};
    normalize();
    return *this;
}

Number& Number::operator-=(const Number& rhs) {
    if (isReal() && rhs.isReal()) {
        std::get<RealNumber>(value_) -= rhs.asReal();
        return *this;
    }

    if (auto* lhs = std::get_if<ComplexNumber>(&value_)) {
        if (rhs.isReal()) {
            lhs->real -= rhs.asReal();
            normalize();
            return *this;
        }
        const auto& other = rhs.asComplex();
        lhs->real -= other.real;
        lhs->imaginary -= other.imaginary;
        normalize();
        return *this;
    }

    const RealNumber lhsReal = asReal();
    const auto& other = rhs.asComplex();
    value_ = ComplexNumber{lhsReal - other.real, -other.imaginary};
    normalize();
    return *this;
}

Number& Number::operator*=(const Number& rhs) {
    if (isReal() && rhs.isReal()) {
        std::get<RealNumber>(value_) *= rhs.asReal();
        return *this;
    }

    if (auto* lhs = std::get_if<ComplexNumber>(&value_)) {
        if (rhs.isReal()) {
            lhs->real *= rhs.asReal();
            lhs->imaginary *= rhs.asReal();
            normalize();
            return *this;
        }

        const RealNumber a = lhs->real;
        const RealNumber b = lhs->imaginary;
        const auto& other = rhs.asComplex();
        // rhsが*this自身でも後半式が前半代入の影響を受けないよう、両成分を先に退避する。
        const RealNumber c = other.real;
        const RealNumber d = other.imaginary;
        lhs->real = a * c - b * d;
        lhs->imaginary = a * d + b * c;
        normalize();
        return *this;
    }

    const RealNumber lhsReal = asReal();
    const auto& other = rhs.asComplex();
    value_ = ComplexNumber{lhsReal * other.real, lhsReal * other.imaginary};
    normalize();
    return *this;
}

Number& Number::operator/=(const Number& rhs) {
    if (rhs.isZero())
        throw std::domain_error("Division by zero");

    if (isReal() && rhs.isReal()) {
        std::get<RealNumber>(value_) /= rhs.asReal();
        return *this;
    }

    if (auto* lhs = std::get_if<ComplexNumber>(&value_); lhs && rhs.isReal()) {
        lhs->real /= rhs.asReal();
        lhs->imaginary /= rhs.asReal();
        normalize();
        return *this;
    }

    const RealNumber a = realPart();
    const RealNumber b = imaginaryPart();
    const auto& other = rhs.asComplex();
    const RealNumber denominator = other.real * other.real + other.imaginary * other.imaginary;

    value_ = ComplexNumber{
        (a * other.real + b * other.imaginary) / denominator,
        (b * other.real - a * other.imaginary) / denominator};
    normalize();
    return *this;
}

bool Number::operator==(const Number& rhs) const {
    if (isReal() && rhs.isReal())
        return asReal() == rhs.asReal();
    if (isComplex() && rhs.isComplex())
        return asComplex() == rhs.asComplex();
    return false;
}

void Number::normalize() {
    auto* complex = std::get_if<ComplexNumber>(&value_);
    if (complex && complex->imaginary.isZero())
        value_ = complex->real;
}

Number operator+(Number lhs, const Number& rhs) {
    lhs += rhs;
    return lhs;
}

Number operator-(Number lhs, const Number& rhs) {
    lhs -= rhs;
    return lhs;
}

Number operator*(Number lhs, const Number& rhs) {
    lhs *= rhs;
    return lhs;
}

Number operator/(Number lhs, const Number& rhs) {
    lhs /= rhs;
    return lhs;
}

Number integerPower(Number base, std::uint64_t exponent) {
    Number result{BigInt{1}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result *= base;
        exponent >>= 1U;
        if (exponent != 0)
            base *= base;
    }
    return result;
}

} // namespace mmcal::numeric
