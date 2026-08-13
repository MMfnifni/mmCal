// 十進近似値metadata
#include "decimal_approximation.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace mmcal::numeric {
namespace {

[[nodiscard]] bool isDivisibleBy(const BigInt& value, std::int64_t divisor) {
    return (value % BigInt{divisor}).isZero();
}

[[nodiscard]] bool hasTerminatingDecimal(BigInt denominator) {
    while (isDivisibleBy(denominator, 2))
        denominator /= BigInt{2};
    while (isDivisibleBy(denominator, 5))
        denominator /= BigInt{5};

    return denominator == BigInt{1};
}

[[nodiscard]] char decimalDigit(const BigInt& value) {
    const std::string text = value.toString();
    if (text.size() != 1)
        throw std::logic_error("Decimal digit is outside the expected range");

    return text.front();
}

void incrementDecimal(std::string& digits, BigInt& integerPart) {
    for (auto iterator = digits.rbegin(); iterator != digits.rend(); ++iterator) {
        if (*iterator != '9') {
            ++*iterator;
            return;
        }

        *iterator = '0';
    }

    integerPart += BigInt{1};
}


[[nodiscard]] Rational exactDecimalValue(std::string_view text) {
    bool negative = false;
    if (!text.empty() && text.front() == '-') {
        negative = true;
        text.remove_prefix(1);
    }

    const std::size_t point = text.find('.');
    if (point == std::string_view::npos) {
        BigInt value = BigInt::parse(std::string{text});
        if (negative)
            value = -value;
        return Rational{std::move(value)};
    }

    std::string digits{text};
    digits.erase(point, 1);
    BigInt numerator = BigInt::parse(digits);
    if (negative)
        numerator = -numerator;

    BigInt denominator{1};
    for (std::size_t i = point; i + 1 < text.size(); ++i)
        denominator *= BigInt{10};
    return Rational{std::move(numerator), std::move(denominator)};
}

[[nodiscard]] bool isZeroText(const BigInt& integerPart, std::string_view digits) {
    if (!integerPart.isZero())
        return false;

    for (const char digit : digits)
        if (digit != '0')
            return false;

    return true;
}

struct FixedDecimal final {
    std::string text;
    std::size_t fractionalDigits = 0;
    bool exact = false;
};

// certified区間由来の固定桁表示では，要求桁まで並んだ末尾0をすべて見せる必要はない。
// ただし近似値であることと，最後に観測された非零桁より一段下まで保証があることを
// 視覚的に残すため，末尾0が複数ある場合は1個だけ保持する。
// exact値由来の有限小数はfromReal側で必要最小桁表示のままとする。
[[nodiscard]] std::size_t compactCertifiedDecimal(std::string& text) {
    const std::size_t point = text.find('.');
    if (point == std::string::npos)
        return 0;

    const std::size_t fractionalDigits = text.size() - point - 1;
    std::size_t trailingZeros = 0;
    while (trailingZeros < fractionalDigits
        && text[text.size() - 1 - trailingZeros] == '0')
        ++trailingZeros;

    if (trailingZeros <= 1)
        return fractionalDigits;

    const std::size_t removed = trailingZeros - 1;
    text.erase(text.size() - removed);
    return fractionalDigits - removed;
}

[[nodiscard]] FixedDecimal roundFixed(
    const Rational& rational,
    std::size_t fractionalDigits) {
    if (fractionalDigits == std::numeric_limits<std::size_t>::max())
        throw std::length_error("Decimal precision is too large");

    const bool negative = rational.numerator().isNegative();
    BigInt numerator = rational.numerator().abs();
    const BigInt denominator = rational.denominator();

    BigInt integerPart = numerator / denominator;
    BigInt remainder = numerator % denominator;

    std::string digits;
    digits.reserve(fractionalDigits + 1);

    // 要求桁 + 1桁（guard digit）まで必ず生成する。途中で割り切れた場合も0を
    // 補うことで、区間端点同士を同一桁数で比較できるようにする。
    for (std::size_t i = 0; i < fractionalDigits + 1; ++i) {
        remainder *= BigInt{10};
        const BigInt digit = remainder / denominator;
        remainder %= denominator;
        digits.push_back(decimalDigit(digit));
    }

    const int guardDigit = digits.back() - '0';
    digits.pop_back();

    // 最近接・偶数丸め。guardが5より大きい、または5でその後に非0が続くなら上へ。
    // 厳密に...5000...で終わるtieだけは最後に保持する桁の偶奇で決める。
    const bool hasDiscardedNonZero = !remainder.isZero();
    const bool lastKeptIsOdd = fractionalDigits == 0
        ? (integerPart % BigInt{2}) != BigInt{}
        : (digits.back() - '0') % 2 != 0;
    const bool roundUp = guardDigit > 5
        || (guardDigit == 5 && (hasDiscardedNonZero || lastKeptIsOdd));
    if (roundUp)
        incrementDecimal(digits, integerPart);

    std::string text;
    if (negative && !isZeroText(integerPart, digits))
        text.push_back('-');

    text += integerPart.toString();
    if (fractionalDigits != 0) {
        text.push_back('.');
        text += digits;
    }

    return FixedDecimal{
        std::move(text),
        fractionalDigits,
        remainder.isZero() && guardDigit == 0
    };
}

} // namespace

DecimalApproximation::DecimalApproximation(
    std::string text,
    std::size_t fractionalDigits,
    std::size_t requestedFractionalDigits,
    bool rounded,
    ApproximationOrigin origin,
    Rational displayedValue,
    Rational certifiedLower,
    Rational certifiedUpper)
    : text_(std::move(text)),
      fractionalDigits_(fractionalDigits),
      requestedFractionalDigits_(requestedFractionalDigits),
      rounded_(rounded),
      origin_(origin),
      displayedValue_(std::move(displayedValue)),
      certifiedLower_(std::move(certifiedLower)),
      certifiedUpper_(std::move(certifiedUpper)) {
    if (text_.empty())
        throw std::invalid_argument("Decimal approximation text cannot be empty");
    if (certifiedLower_ > certifiedUpper_)
        throw std::invalid_argument("Decimal approximation enclosure is reversed");
}

DecimalApproximation DecimalApproximation::fromReal(
    const RealNumber& value,
    std::size_t repeatingFractionalDigits) {
    if (repeatingFractionalDigits == 0)
        throw std::invalid_argument("Decimal precision must be greater than zero");

    const Rational rational = value.toRational();
    const bool negative = rational.numerator().isNegative();
    BigInt numerator = rational.numerator().abs();
    const BigInt denominator = rational.denominator();

    BigInt integerPart = numerator / denominator;
    BigInt remainder = numerator % denominator;
    if (remainder.isZero())
        return DecimalApproximation{
            rational.numerator().toString(), 0, repeatingFractionalDigits, false,
            ApproximationOrigin::ExactValue, rational, rational, rational};

    const bool terminating = hasTerminatingDecimal(denominator);
    std::string digits;

    if (!terminating) {
        if (repeatingFractionalDigits == std::numeric_limits<std::size_t>::max())
            throw std::length_error("Decimal precision is too large");
        digits.reserve(repeatingFractionalDigits + 1);
    }

    do {
        remainder *= BigInt{10};
        const BigInt digit = remainder / denominator;
        remainder %= denominator;
        digits.push_back(decimalDigit(digit));
    } while (terminating
        ? !remainder.isZero()
        : digits.size() < repeatingFractionalDigits + 1);

    bool rounded = false;
    if (!terminating) {
        const int guardDigit = digits.back() - '0';
        digits.pop_back();

        // 近似値は最近接へ丸め、厳密な中間値では偶数丸めを用いる。
        const bool hasDiscardedNonZero = !remainder.isZero();
        const bool lastKeptIsOdd = (digits.back() - '0') % 2 != 0;
        const bool roundUp = guardDigit > 5
            || (guardDigit == 5 && (hasDiscardedNonZero || lastKeptIsOdd));
        if (roundUp)
            incrementDecimal(digits, integerPart);

        rounded = true;
    }

    std::string text;
    if (negative && !isZeroText(integerPart, digits))
        text.push_back('-');

    text += integerPart.toString();
    text.push_back('.');
    text += digits;

    return DecimalApproximation{
        text, digits.size(), repeatingFractionalDigits, rounded,
        ApproximationOrigin::ExactValue, exactDecimalValue(text), rational, rational};
}

DecimalApproximation DecimalApproximation::fromRealFixed(
    const RealNumber& value,
    std::size_t fractionalDigits) {
    const FixedDecimal rounded = roundFixed(value.toRational(), fractionalDigits);
    const Rational rational = value.toRational();
    return DecimalApproximation{
        rounded.text,
        rounded.fractionalDigits,
        fractionalDigits,
        !rounded.exact,
        ApproximationOrigin::ExactValue,
        exactDecimalValue(rounded.text),
        rational,
        rational
    };
}

std::optional<DecimalApproximation> DecimalApproximation::fromCertifiedInterval(
    const Rational& lower,
    const Rational& upper,
    std::size_t fractionalDigits) {
    if ((lower <=> upper) == std::strong_ordering::greater)
        throw std::invalid_argument("Certified decimal interval is reversed");

    const FixedDecimal lowerRounded = roundFixed(lower, fractionalDigits);
    const FixedDecimal upperRounded = roundFixed(upper, fractionalDigits);
    if (lowerRounded.text != upperRounded.text)
        return std::nullopt;

    // 区間両端が同じ丸め結果を持つため、その間にある真値も必ず同じ結果になる。
    // 値自体が厳密に10進有限であるとは限らないので、certified interval由来は常に「丸められた近似値」として扱う。
    // 内部保証は要求桁数のまま保持し，表示だけ末尾0を1桁まで圧縮する。
    std::string text = lowerRounded.text;
    const Rational displayed = exactDecimalValue(text);
    const std::size_t displayedFractionalDigits = compactCertifiedDecimal(text);
    return DecimalApproximation{
        std::move(text),
        displayedFractionalDigits,
        fractionalDigits,
        true,
        ApproximationOrigin::CertifiedInterval,
        displayed,
        lower,
        upper
    };
}

std::string_view DecimalApproximation::text() const noexcept {
    return text_;
}

std::size_t DecimalApproximation::fractionalDigits() const noexcept {
    return fractionalDigits_;
}

std::size_t DecimalApproximation::requestedFractionalDigits() const noexcept {
    return requestedFractionalDigits_;
}

bool DecimalApproximation::isRounded() const noexcept {
    return rounded_;
}

ApproximationOrigin DecimalApproximation::origin() const noexcept {
    return origin_;
}

const Rational& DecimalApproximation::displayedValue() const noexcept {
    return displayedValue_;
}

const Rational& DecimalApproximation::certifiedLower() const noexcept {
    return certifiedLower_;
}

const Rational& DecimalApproximation::certifiedUpper() const noexcept {
    return certifiedUpper_;
}

bool DecimalApproximation::enclosureIsPoint() const noexcept {
    return certifiedLower_ == certifiedUpper_;
}

} // namespace mmcal::numeric
