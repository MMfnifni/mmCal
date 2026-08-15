// 十進近似値metadata
#include "decimal_approximation.hpp"

#include "integer_algorithms.hpp"

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


struct InformationBounds final {
    Rational lower;
    Rational upper;
};

struct FixedDecimal final {
    std::string text;
    std::size_t fractionalDigits = 0;
    bool exact = false;
};

[[nodiscard]] FixedDecimal roundFixed(
    const Rational& rational,
    std::size_t fractionalDigits);

struct SignificantDecimal final {
    std::string text;
    std::size_t fractionalDigits = 0;
    std::size_t roundingFractionalDigits = 0;
    Rational displayed;
    Rational halfQuantum;
    bool exact = false;
};

[[nodiscard]] Rational informationHalfQuantum(std::size_t fractionalDigits) {
    const auto exponent = static_cast<std::uint64_t>(fractionalDigits);
    BigInt denominator = pow(BigInt{10}, exponent);
    denominator *= BigInt{2};
    return Rational{BigInt{1}, std::move(denominator)};
}

[[nodiscard]] std::int64_t checkedSignedSize(std::size_t value) {
    if (value > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error("Decimal exponent is too large");
    return static_cast<std::int64_t>(value);
}

[[nodiscard]] std::int64_t decimalExponent(const Rational& value) {
    if (value.numerator().isZero())
        return 0;

    const BigInt numerator = value.numerator().abs();
    const BigInt denominator = value.denominator();
    const std::size_t numeratorDigits = numerator.toString().size();
    const std::size_t denominatorDigits = denominator.toString().size();

    if (numeratorDigits >= denominatorDigits) {
        const std::size_t gap = numeratorDigits - denominatorDigits;
        const BigInt threshold = denominator * pow(BigInt{10}, static_cast<std::uint64_t>(gap));
        const std::int64_t candidate = checkedSignedSize(gap);
        return numerator < threshold ? candidate - 1 : candidate;
    }

    const std::size_t gap = denominatorDigits - numeratorDigits;
    const BigInt scaled = numerator * pow(BigInt{10}, static_cast<std::uint64_t>(gap));
    const std::int64_t candidate = -checkedSignedSize(gap);
    return scaled < denominator ? candidate - 1 : candidate;
}

[[nodiscard]] Rational powerOfTen(std::int64_t exponent) {
    if (exponent >= 0)
        return Rational{pow(BigInt{10}, static_cast<std::uint64_t>(exponent))};
    if (exponent == std::numeric_limits<std::int64_t>::min())
        throw std::length_error("Decimal exponent is too large");
    return Rational{BigInt{1}, pow(BigInt{10}, static_cast<std::uint64_t>(-exponent))};
}

[[nodiscard]] std::size_t trimTrailingFractionalZeros(std::string& text) {
    const std::size_t point = text.find('.');
    if (point == std::string::npos)
        return 0;

    while (text.size() > point + 1 && text.back() == '0')
        text.pop_back();
    if (text.back() == '.')
        text.pop_back();
    return text.find('.') == std::string::npos ? 0 : text.size() - text.find('.') - 1;
}

[[nodiscard]] SignificantDecimal roundSignificant(
    const Rational& value,
    std::size_t significantDigits) {
    if (significantDigits == 0)
        throw std::invalid_argument("Decimal precision must be greater than zero");

    const std::int64_t exponent = value.numerator().isZero()
        ? -checkedSignedSize(significantDigits)
        : decimalExponent(value) - checkedSignedSize(significantDigits) + 1;
    const Rational quantum = powerOfTen(exponent);
    const FixedDecimal scaledRounded = roundFixed(value / quantum, 0);
    const BigInt roundedInteger = BigInt::parse(scaledRounded.text);
    const Rational displayed = Rational{roundedInteger} * quantum;

    const std::int64_t displayedQuantumExponent = displayed.numerator().isZero()
        ? -checkedSignedSize(significantDigits)
        : decimalExponent(displayed) - checkedSignedSize(significantDigits) + 1;
    const Rational displayedQuantum = powerOfTen(displayedQuantumExponent);

    std::size_t fractionalDigits = 0;
    if (displayedQuantumExponent < 0) {
        if (displayedQuantumExponent == std::numeric_limits<std::int64_t>::min())
            throw std::length_error("Decimal precision is too large");
        fractionalDigits = static_cast<std::size_t>(-displayedQuantumExponent);
    }
    const FixedDecimal displayedFixed = roundFixed(displayed, fractionalDigits);
    std::string text = displayedFixed.text;
    const std::size_t compactDigits = trimTrailingFractionalZeros(text);

    return SignificantDecimal{
        std::move(text),
        compactDigits,
        fractionalDigits,
        displayed,
        displayedQuantum / Rational{BigInt{2}},
        displayed == value
    };
}


struct ZeroCenteredDecimal final {
    std::string text;
    std::size_t fractionalDigits = 0;
    Rational halfQuantum;
};

[[nodiscard]] std::optional<ZeroCenteredDecimal> zeroCenteredDecimal(
    const Rational& lower,
    const Rational& upper,
    std::size_t significantDigits) {
    if (lower > Rational{} || upper < Rational{})
        return std::nullopt;

    const Rational lowerMagnitude = lower.numerator().isNegative() ? -lower : lower;
    const Rational upperMagnitude = upper.numerator().isNegative() ? -upper : upper;
    const Rational maximumMagnitude = lowerMagnitude < upperMagnitude
        ? upperMagnitude : lowerMagnitude;

    std::size_t fractionalDigits = significantDigits;
    if (!maximumMagnitude.isZero()) {
        const Rational twiceMagnitude = maximumMagnitude * Rational{BigInt{2}};
        const std::int64_t floorExponent = decimalExponent(twiceMagnitude);
        const std::int64_t ceilExponent = twiceMagnitude == powerOfTen(floorExponent)
            ? floorExponent : floorExponent + 1;
        if (ceilExponent > 0)
            return std::nullopt;
        fractionalDigits = ceilExponent < 0
            ? static_cast<std::size_t>(-ceilExponent)
            : 0;
    }

    const FixedDecimal zero = roundFixed(Rational{}, fractionalDigits);
    return ZeroCenteredDecimal{
        zero.text,
        fractionalDigits,
        informationHalfQuantum(fractionalDigits)
    };
}

[[nodiscard]] InformationBounds defaultInformationBoundsSignificant(
    const SignificantDecimal& rounded,
    const Rational& certifiedLower,
    const Rational& certifiedUpper) {
    const Rational roundedLower = rounded.displayed - rounded.halfQuantum;
    const Rational roundedUpper = rounded.displayed + rounded.halfQuantum;
    return InformationBounds{
        certifiedLower < roundedLower ? certifiedLower : roundedLower,
        certifiedUpper > roundedUpper ? certifiedUpper : roundedUpper
    };
}

[[nodiscard]] InformationBounds explicitInformationBoundsSignificant(
    const SignificantDecimal& rounded,
    const Rational& certifiedLower,
    const Rational& certifiedUpper,
    const Rational& informationLower,
    const Rational& informationUpper) {
    if (informationLower > informationUpper)
        throw std::invalid_argument("Decimal approximation information enclosure is reversed");
    if (informationLower > certifiedLower || informationUpper < certifiedUpper)
        throw std::invalid_argument("Decimal approximation information enclosure must contain the certified enclosure");

    const Rational roundedLower = rounded.displayed - rounded.halfQuantum;
    const Rational roundedUpper = rounded.displayed + rounded.halfQuantum;
    return InformationBounds{
        informationLower < roundedLower ? informationLower : roundedLower,
        informationUpper > roundedUpper ? informationUpper : roundedUpper
    };
}

[[nodiscard]] InformationBounds defaultInformationBounds(
    const Rational& displayed,
    const Rational& certifiedLower,
    const Rational& certifiedUpper,
    std::size_t fractionalDigits) {
    const Rational halfQuantum = informationHalfQuantum(fractionalDigits);
    const Rational roundedLower = displayed - halfQuantum;
    const Rational roundedUpper = displayed + halfQuantum;
    return InformationBounds{
        certifiedLower < roundedLower ? certifiedLower : roundedLower,
        certifiedUpper > roundedUpper ? certifiedUpper : roundedUpper
    };
}

[[nodiscard]] InformationBounds explicitInformationBounds(
    const Rational& displayed,
    const Rational& certifiedLower,
    const Rational& certifiedUpper,
    const Rational& informationLower,
    const Rational& informationUpper,
    std::size_t fractionalDigits) {
    if (informationLower > informationUpper)
        throw std::invalid_argument("Decimal approximation information enclosure is reversed");
    if (informationLower > certifiedLower || informationUpper < certifiedUpper)
        throw std::invalid_argument("Decimal approximation information enclosure must contain the certified enclosure");

    const Rational halfQuantum = informationHalfQuantum(fractionalDigits);
    const Rational roundedLower = displayed - halfQuantum;
    const Rational roundedUpper = displayed + halfQuantum;
    return InformationBounds{
        informationLower < roundedLower ? informationLower : roundedLower,
        informationUpper > roundedUpper ? informationUpper : roundedUpper
    };
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
    std::size_t requestedSignificantDigits,
    bool rounded,
    ApproximationOrigin origin,
    Rational displayedValue,
    Rational certifiedLower,
    Rational certifiedUpper,
    Rational informationLower,
    Rational informationUpper)
    : text_(std::move(text)),
      fractionalDigits_(fractionalDigits),
      requestedFractionalDigits_(requestedFractionalDigits),
      requestedSignificantDigits_(requestedSignificantDigits),
      rounded_(rounded),
      origin_(origin),
      displayedValue_(std::move(displayedValue)),
      certifiedLower_(std::move(certifiedLower)),
      certifiedUpper_(std::move(certifiedUpper)),
      informationLower_(std::move(informationLower)),
      informationUpper_(std::move(informationUpper)) {
    if (text_.empty())
        throw std::invalid_argument("Decimal approximation text cannot be empty");
    if (certifiedLower_ > certifiedUpper_)
        throw std::invalid_argument("Decimal approximation enclosure is reversed");
    if (informationLower_ > informationUpper_)
        throw std::invalid_argument("Decimal approximation information enclosure is reversed");
    if (informationLower_ > certifiedLower_ || informationUpper_ < certifiedUpper_)
        throw std::invalid_argument("Decimal approximation information enclosure must contain the certified enclosure");
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
    if (remainder.isZero()) {
        const std::string text = rational.numerator().toString();
        const InformationBounds information = defaultInformationBounds(
            rational, rational, rational, repeatingFractionalDigits);
        return DecimalApproximation{
            text, 0, repeatingFractionalDigits, 0, false, ApproximationOrigin::ExactValue,
            rational, rational, rational, information.lower, information.upper};
    }

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

    const Rational displayed = exactDecimalValue(text);
    const InformationBounds information = defaultInformationBounds(
        displayed, rational, rational, repeatingFractionalDigits);
    return DecimalApproximation{
        text, digits.size(), repeatingFractionalDigits, 0, rounded,
        ApproximationOrigin::ExactValue, displayed, rational, rational,
        information.lower, information.upper};
}

DecimalApproximation DecimalApproximation::fromRealSignificant(
    const RealNumber& value,
    std::size_t significantDigits) {
    const Rational rational = value.toRational();
    const SignificantDecimal rounded = roundSignificant(rational, significantDigits);
    const InformationBounds information = defaultInformationBoundsSignificant(
        rounded, rational, rational);
    return DecimalApproximation{
        rounded.text,
        rounded.fractionalDigits,
        rounded.roundingFractionalDigits,
        significantDigits,
        !rounded.exact,
        ApproximationOrigin::ExactValue,
        rounded.displayed,
        rational,
        rational,
        information.lower,
        information.upper
    };
}

DecimalApproximation DecimalApproximation::fromRealFixed(
    const RealNumber& value,
    std::size_t fractionalDigits) {
    const FixedDecimal rounded = roundFixed(value.toRational(), fractionalDigits);
    const Rational rational = value.toRational();
    const Rational displayed = exactDecimalValue(rounded.text);
    const InformationBounds information = defaultInformationBounds(
        displayed, rational, rational, fractionalDigits);
    return DecimalApproximation{
        rounded.text,
        rounded.fractionalDigits,
        fractionalDigits,
        0,
        !rounded.exact,
        ApproximationOrigin::ExactValue,
        displayed,
        rational,
        rational,
        information.lower,
        information.upper
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

    std::string text = lowerRounded.text;
    const Rational displayed = exactDecimalValue(text);
    const InformationBounds information = defaultInformationBounds(
        displayed, lower, upper, fractionalDigits);
    const std::size_t displayedFractionalDigits = compactCertifiedDecimal(text);
    return DecimalApproximation{
        std::move(text),
        displayedFractionalDigits,
        fractionalDigits,
        0,
        true,
        ApproximationOrigin::CertifiedInterval,
        displayed,
        lower,
        upper,
        information.lower,
        information.upper
    };
}

std::optional<DecimalApproximation> DecimalApproximation::fromCertifiedIntervalSignificant(
    const Rational& lower,
    const Rational& upper,
    std::size_t significantDigits) {
    if (lower > upper)
        throw std::invalid_argument("Certified decimal interval is reversed");

    const SignificantDecimal lowerRounded = roundSignificant(lower, significantDigits);
    const SignificantDecimal upperRounded = roundSignificant(upper, significantDigits);
    if (lowerRounded.displayed != upperRounded.displayed) {
        const auto zero = zeroCenteredDecimal(lower, upper, significantDigits);
        if (!zero)
            return std::nullopt;
        const Rational roundedLower = -zero->halfQuantum;
        const Rational roundedUpper = zero->halfQuantum;
        const InformationBounds information{
            lower < roundedLower ? lower : roundedLower,
            upper > roundedUpper ? upper : roundedUpper
        };
        return DecimalApproximation{
            zero->text,
            zero->fractionalDigits,
            zero->fractionalDigits,
            significantDigits,
            true,
            ApproximationOrigin::CertifiedInterval,
            Rational{},
            lower,
            upper,
            information.lower,
            information.upper
        };
    }

    const InformationBounds information = defaultInformationBoundsSignificant(
        lowerRounded, lower, upper);
    return DecimalApproximation{
        lowerRounded.text,
        lowerRounded.fractionalDigits,
        lowerRounded.roundingFractionalDigits,
        significantDigits,
        true,
        ApproximationOrigin::CertifiedInterval,
        lowerRounded.displayed,
        lower,
        upper,
        information.lower,
        information.upper
    };
}

std::optional<DecimalApproximation> DecimalApproximation::fromCertifiedIntervalWithInformation(
    const Rational& certifiedLower,
    const Rational& certifiedUpper,
    const Rational& informationLower,
    const Rational& informationUpper,
    std::size_t fractionalDigits) {
    if (certifiedLower > certifiedUpper)
        throw std::invalid_argument("Certified decimal interval is reversed");

    // certified truthがpointなら，既存Nと同じく有限10進値の不要な末尾0を増やさない。
    // InformationEnclosureは別metadataとして保持するため，最小表示にしても情報量を回収しない。
    if (certifiedLower == certifiedUpper) {
        const RealNumber exact{certifiedLower};
        const DecimalApproximation base = fractionalDigits == 0
            ? fromRealFixed(exact, 0)
            : fromReal(exact, fractionalDigits);
        const InformationBounds information = explicitInformationBounds(
            base.displayedValue(), certifiedLower, certifiedUpper,
            informationLower, informationUpper, fractionalDigits);
        return DecimalApproximation{
            std::string{base.text()},
            base.fractionalDigits(),
            fractionalDigits,
            0,
            base.isRounded(),
            ApproximationOrigin::CertifiedInterval,
            base.displayedValue(),
            certifiedLower,
            certifiedUpper,
            information.lower,
            information.upper
        };
    }

    const FixedDecimal lowerRounded = roundFixed(certifiedLower, fractionalDigits);
    const FixedDecimal upperRounded = roundFixed(certifiedUpper, fractionalDigits);
    if (lowerRounded.text != upperRounded.text)
        return std::nullopt;

    std::string text = lowerRounded.text;
    const Rational displayed = exactDecimalValue(text);
    const InformationBounds information = explicitInformationBounds(
        displayed, certifiedLower, certifiedUpper, informationLower, informationUpper,
        fractionalDigits);
    const std::size_t displayedFractionalDigits = compactCertifiedDecimal(text);
    return DecimalApproximation{
        std::move(text),
        displayedFractionalDigits,
        fractionalDigits,
        0,
        true,
        ApproximationOrigin::CertifiedInterval,
        displayed,
        certifiedLower,
        certifiedUpper,
        information.lower,
        information.upper
    };
}

std::optional<DecimalApproximation> DecimalApproximation::fromCertifiedIntervalWithInformationSignificant(
    const Rational& certifiedLower,
    const Rational& certifiedUpper,
    const Rational& informationLower,
    const Rational& informationUpper,
    std::size_t significantDigits) {
    if (certifiedLower > certifiedUpper)
        throw std::invalid_argument("Certified decimal interval is reversed");

    const SignificantDecimal lowerRounded = roundSignificant(certifiedLower, significantDigits);
    const SignificantDecimal upperRounded = roundSignificant(certifiedUpper, significantDigits);
    if (lowerRounded.displayed != upperRounded.displayed) {
        if (informationLower > informationUpper)
            throw std::invalid_argument("Decimal approximation information enclosure is reversed");
        if (informationLower > certifiedLower || informationUpper < certifiedUpper)
            throw std::invalid_argument("Decimal approximation information enclosure must contain the certified enclosure");
        const auto zero = zeroCenteredDecimal(informationLower, informationUpper, significantDigits);
        if (!zero)
            return std::nullopt;
        const Rational roundedLower = -zero->halfQuantum;
        const Rational roundedUpper = zero->halfQuantum;
        const InformationBounds information{
            informationLower < roundedLower ? informationLower : roundedLower,
            informationUpper > roundedUpper ? informationUpper : roundedUpper
        };
        return DecimalApproximation{
            zero->text,
            zero->fractionalDigits,
            zero->fractionalDigits,
            significantDigits,
            true,
            ApproximationOrigin::CertifiedInterval,
            Rational{},
            certifiedLower,
            certifiedUpper,
            information.lower,
            information.upper
        };
    }

    const InformationBounds information = explicitInformationBoundsSignificant(
        lowerRounded,
        certifiedLower,
        certifiedUpper,
        informationLower,
        informationUpper);
    return DecimalApproximation{
        lowerRounded.text,
        lowerRounded.fractionalDigits,
        lowerRounded.roundingFractionalDigits,
        significantDigits,
        true,
        ApproximationOrigin::CertifiedInterval,
        lowerRounded.displayed,
        certifiedLower,
        certifiedUpper,
        information.lower,
        information.upper
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

std::size_t DecimalApproximation::requestedSignificantDigits() const noexcept {
    return requestedSignificantDigits_;
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

const Rational& DecimalApproximation::informationLower() const noexcept {
    return informationLower_;
}

const Rational& DecimalApproximation::informationUpper() const noexcept {
    return informationUpper_;
}

bool DecimalApproximation::certifiedEnclosureIsPoint() const noexcept {
    return certifiedLower_ == certifiedUpper_;
}

bool DecimalApproximation::informationEnclosureIsPoint() const noexcept {
    return informationLower_ == informationUpper_;
}

} // namespace mmcal::numeric
