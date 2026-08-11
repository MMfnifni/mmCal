// 整数・有理数を統合するRealNumber
#include "real_number.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::numeric {

RealNumber::RealNumber(BigInt integer)
    : value_(std::move(integer)) {}

RealNumber::RealNumber(Rational rational)
    : value_(std::move(rational)) {
    normalize();
}

bool RealNumber::isInteger() const noexcept {
    return std::holds_alternative<BigInt>(value_);
}

bool RealNumber::isRational() const noexcept {
    return std::holds_alternative<Rational>(value_);
}

bool RealNumber::isZero() const noexcept {
    if (const auto* integer = std::get_if<BigInt>(&value_))
        return integer->isZero();

    return std::get<Rational>(value_).isZero();
}

bool RealNumber::isNegative() const noexcept {
    if (const auto* integer = std::get_if<BigInt>(&value_))
        return integer->isNegative();

    return std::get<Rational>(value_).numerator().isNegative();
}

const BigInt& RealNumber::asInteger() const {
    if (!isInteger())
        throw std::logic_error("RealNumber does not contain an integer");

    return std::get<BigInt>(value_);
}

const Rational& RealNumber::asRational() const {
    if (!isRational())
        throw std::logic_error("RealNumber does not contain a rational value");

    return std::get<Rational>(value_);
}

Rational RealNumber::toRational() const {
    if (const auto* integer = std::get_if<BigInt>(&value_))
        return Rational{*integer};

    return std::get<Rational>(value_);
}

RealNumber RealNumber::abs() const {
    return isNegative() ? -*this : *this;
}

std::string RealNumber::toString(unsigned radix) const {
    if (const auto* integer = std::get_if<BigInt>(&value_))
        return integer->toString(radix);

    return std::get<Rational>(value_).toString(radix);
}

RealNumber RealNumber::operator-() const {
    if (const auto* integer = std::get_if<BigInt>(&value_))
        return RealNumber{-*integer};

    return RealNumber{-std::get<Rational>(value_)};
}

RealNumber& RealNumber::operator+=(const RealNumber& rhs) {
    if (isInteger() && rhs.isInteger()) {
        std::get<BigInt>(value_) += rhs.asInteger();
        return *this;
    }

    if (auto* lhs = std::get_if<Rational>(&value_)) {
        if (rhs.isInteger())
            *lhs += Rational{rhs.asInteger()};
        else
            *lhs += rhs.asRational();
    }
    else {
        Rational result{asInteger()};
        result += rhs.asRational();
        value_ = std::move(result);
    }
    normalize();
    return *this;
}

RealNumber& RealNumber::operator-=(const RealNumber& rhs) {
    if (isInteger() && rhs.isInteger()) {
        std::get<BigInt>(value_) -= rhs.asInteger();
        return *this;
    }

    if (auto* lhs = std::get_if<Rational>(&value_)) {
        if (rhs.isInteger())
            *lhs -= Rational{rhs.asInteger()};
        else
            *lhs -= rhs.asRational();
    }
    else {
        Rational result{asInteger()};
        result -= rhs.asRational();
        value_ = std::move(result);
    }
    normalize();
    return *this;
}

RealNumber& RealNumber::operator*=(const RealNumber& rhs) {
    if (isInteger() && rhs.isInteger()) {
        std::get<BigInt>(value_) *= rhs.asInteger();
        return *this;
    }

    if (auto* lhs = std::get_if<Rational>(&value_)) {
        if (rhs.isInteger())
            *lhs *= Rational{rhs.asInteger()};
        else
            *lhs *= rhs.asRational();
    }
    else {
        Rational result{asInteger()};
        result *= rhs.asRational();
        value_ = std::move(result);
    }
    normalize();
    return *this;
}

RealNumber& RealNumber::operator/=(const RealNumber& rhs) {
    if (rhs.isZero())
        throw std::domain_error("Division by zero");

    if (isInteger() && rhs.isInteger()) {
        value_ = Rational{asInteger(), rhs.asInteger()};
        normalize();
        return *this;
    }

    if (auto* lhs = std::get_if<Rational>(&value_)) {
        if (rhs.isInteger())
            *lhs /= Rational{rhs.asInteger()};
        else
            *lhs /= rhs.asRational();
    }
    else {
        Rational result{asInteger()};
        result /= rhs.asRational();
        value_ = std::move(result);
    }
    normalize();
    return *this;
}

std::strong_ordering RealNumber::operator<=>(const RealNumber& rhs) const {
    if (isInteger() && rhs.isInteger())
        return asInteger() <=> rhs.asInteger();

    return toRational() <=> rhs.toRational();
}

bool RealNumber::operator==(const RealNumber& rhs) const {
    if (isInteger() && rhs.isInteger())
        return asInteger() == rhs.asInteger();

    return toRational() == rhs.toRational();
}

void RealNumber::normalize() {
    auto* rational = std::get_if<Rational>(&value_);
    if (!rational || !rational->isInteger())
        return;

    value_ = rational->numerator();
}

RealNumber operator+(RealNumber lhs, const RealNumber& rhs) {
    lhs += rhs;
    return lhs;
}

RealNumber operator-(RealNumber lhs, const RealNumber& rhs) {
    lhs -= rhs;
    return lhs;
}

RealNumber operator*(RealNumber lhs, const RealNumber& rhs) {
    lhs *= rhs;
    return lhs;
}

RealNumber operator/(RealNumber lhs, const RealNumber& rhs) {
    lhs /= rhs;
    return lhs;
}

} // namespace mmcal::numeric
