// 任意精度有理数Rational
#include "rational.hpp"

#include "integer_algorithms.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::numeric {
namespace {

void validateRadix(unsigned radix) {
    if (radix < 2 || radix > 36)
        throw std::invalid_argument("Rational radix must be in the range 2..36");
}

} // namespace

Rational::Rational() = default;

Rational::Rational(BigInt integer)
    : numerator_(std::move(integer)) {}

Rational::Rational(BigInt numerator, BigInt denominator)
    : numerator_(std::move(numerator)), denominator_(std::move(denominator)) {
    normalize();
}

Rational Rational::parse(std::string_view text, unsigned radix) {
    validateRadix(radix);

    if (text.empty())
        throw std::invalid_argument("Rational cannot parse an empty string");

    bool negative = false;
    std::size_t position = 0;

    if (text.front() == '+' || text.front() == '-') {
        negative = text.front() == '-';
        position = 1;
    }

    if (position == text.size())
        throw std::invalid_argument("Rational requires at least one digit");

    const auto unsignedText = text.substr(position);
    const auto point = unsignedText.find('.');

    if (point == std::string_view::npos) {
        BigInt integer = BigInt::parse(unsignedText, radix);
        if (negative)
            integer = -integer;
        return Rational{std::move(integer)};
    }

    if (unsignedText.find('.', point + 1) != std::string_view::npos)
        throw std::invalid_argument("Rational contains multiple radix points");

    const auto whole = unsignedText.substr(0, point);
    const auto fraction = unsignedText.substr(point + 1);
    if (whole.empty() && fraction.empty())
        throw std::invalid_argument("Rational requires at least one digit");

    std::string digits;
    digits.reserve(whole.size() + fraction.size());
    digits.append(whole);
    digits.append(fraction);

    BigInt numerator = BigInt::parse(digits.empty() ? "0" : digits, radix);
    if (negative)
        numerator = -numerator;

    const BigInt denominator = pow(
        BigInt{static_cast<std::int64_t>(radix)},
        static_cast<std::uint64_t>(fraction.size()));

    return Rational{std::move(numerator), denominator};
}

const BigInt& Rational::numerator() const noexcept {
    return numerator_;
}

const BigInt& Rational::denominator() const noexcept {
    return denominator_;
}

bool Rational::isZero() const noexcept {
    return numerator_.isZero();
}

bool Rational::isInteger() const noexcept {
    return denominator_ == BigInt{1};
}

std::string Rational::toString(unsigned radix) const {
    if (isInteger())
        return numerator_.toString(radix);

    return numerator_.toString(radix) + '/' + denominator_.toString(radix);
}

Rational Rational::operator-() const {
    return Rational{-numerator_, denominator_};
}

Rational& Rational::operator+=(const Rational& rhs) {
    if (rhs.isZero())
        return *this;
    if (isZero()) {
        *this = rhs;
        return *this;
    }

    // a/b + c/d。入力は既約なので、和の分子と最小公倍分母に残り得る共通因子は g=gcd(b,d) の約数だけである。
    // 巨大な最終分母全体とのGCDを取り直さず、gとのGCDだけで完全に既約化できる。
    const BigInt common = gcd(denominator_, rhs.denominator_);
    const BigInt lhsScale = rhs.denominator_ / common;
    const BigInt rhsScale = denominator_ / common;
    BigInt resultNumerator = numerator_ * lhsScale + rhs.numerator_ * rhsScale;
    const BigInt reduction = gcd(resultNumerator.abs(), common);
    resultNumerator /= reduction;
    BigInt resultDenominator = (denominator_ / reduction) * lhsScale;

    numerator_ = std::move(resultNumerator);
    denominator_ = std::move(resultDenominator);
    normalizeSignAndZero();
    return *this;
}

Rational& Rational::operator-=(const Rational& rhs) {
    return *this += -rhs;
}

Rational& Rational::operator*=(const Rational& rhs) {
    if (isZero() || rhs.isZero()) {
        numerator_ = BigInt{};
        denominator_ = BigInt{1};
        return *this;
    }

    // 乗算前に交差約分する。入力Rationalが既約であるため、この2回の交差約分後は結果も既約であり、巨大な積に対する3回目のGCDは不要。
    const BigInt leftCancel = gcd(numerator_.abs(), rhs.denominator_);
    const BigInt rightCancel = gcd(rhs.numerator_.abs(), denominator_);
    BigInt resultNumerator =
        (numerator_ / leftCancel) * (rhs.numerator_ / rightCancel);
    BigInt resultDenominator =
        (denominator_ / rightCancel) * (rhs.denominator_ / leftCancel);

    numerator_ = std::move(resultNumerator);
    denominator_ = std::move(resultDenominator);
    normalizeSignAndZero();
    return *this;
}

Rational& Rational::operator/=(const Rational& rhs) {
    if (rhs.numerator_.isZero())
        throw std::domain_error("Rational division by zero");
    if (isZero())
        return *this;

    // a/b ÷ c/d = ad/bc。分子同士・分母同士を先に約分すれば、入力が既約である限り結果も既約なので、積を作った後のGCD再計算は不要。
    const BigInt numeratorCancel = gcd(numerator_.abs(), rhs.numerator_.abs());
    const BigInt denominatorCancel = gcd(denominator_, rhs.denominator_);
    BigInt resultNumerator = (numerator_ / numeratorCancel)
        * (rhs.denominator_ / denominatorCancel);
    BigInt resultDenominator = (denominator_ / denominatorCancel)
        * (rhs.numerator_ / numeratorCancel);

    numerator_ = std::move(resultNumerator);
    denominator_ = std::move(resultDenominator);
    normalizeSignAndZero();
    return *this;
}

std::strong_ordering Rational::operator<=>(const Rational& rhs) const {
    return numerator_ * rhs.denominator_ <=> rhs.numerator_ * denominator_;
}

bool Rational::operator==(const Rational& rhs) const noexcept {
    return numerator_ == rhs.numerator_ && denominator_ == rhs.denominator_;
}

void Rational::normalize() {
    if (denominator_.isZero())
        throw std::domain_error("Rational denominator cannot be zero");

    normalizeSignAndZero();
    if (numerator_.isZero())
        return;

    const BigInt factor = gcd(numerator_.abs(), denominator_);
    numerator_ /= factor;
    denominator_ /= factor;
}

void Rational::normalizeSignAndZero() {
    if (numerator_.isZero()) {
        denominator_ = BigInt{1};
        return;
    }
    if (denominator_.isNegative()) {
        numerator_ = -numerator_;
        denominator_ = -denominator_;
    }
}

Rational operator+(Rational lhs, const Rational& rhs) {
    lhs += rhs;
    return lhs;
}

Rational operator-(Rational lhs, const Rational& rhs) {
    lhs -= rhs;
    return lhs;
}

Rational operator*(Rational lhs, const Rational& rhs) {
    lhs *= rhs;
    return lhs;
}

Rational operator/(Rational lhs, const Rational& rhs) {
    lhs /= rhs;
    return lhs;
}

} // namespace mmcal::numeric
