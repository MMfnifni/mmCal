// Exprとcertified interval backendの共通変換
#include "expression_interval.hpp"

#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/real_number.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcal::approximation {
namespace {

[[nodiscard]] numeric::Rational semanticHalfQuantum(
    const numeric::DecimalApproximation& value) {
    const auto digits = static_cast<std::uint64_t>(value.requestedFractionalDigits());
    numeric::BigInt denominator = numeric::pow(numeric::BigInt{10}, digits);
    denominator *= numeric::BigInt{2};
    return numeric::Rational{numeric::BigInt{1}, std::move(denominator)};
}

[[nodiscard]] RealInterval intervalFromDecimal(
    const numeric::DecimalApproximation& value,
    std::size_t precisionBits) {
    return RealInterval::fromRationalBounds(
        value.certifiedLower(), value.certifiedUpper(), precisionBits);
}

[[nodiscard]] RealInterval semanticIntervalFromDecimal(
    const numeric::DecimalApproximation& value,
    std::size_t precisionBits) {
    // certified source enclosureはguard bits分だけ要求桁より狭いことがある。
    // semantic側では±0.5*10^-nも包含し，後続演算が宣言済みaccuracy以上の
    // 情報を回収しないための上限としてだけ利用する。実際の包含保証は別のintervalで保持する。
    const numeric::Rational halfQuantum = semanticHalfQuantum(value);
    const numeric::Rational semanticLower = value.displayedValue() - halfQuantum;
    const numeric::Rational semanticUpper = value.displayedValue() + halfQuantum;
    const numeric::Rational lower = value.certifiedLower() < semanticLower
        ? value.certifiedLower() : semanticLower;
    const numeric::Rational upper = value.certifiedUpper() > semanticUpper
        ? value.certifiedUpper() : semanticUpper;
    return RealInterval::fromRationalBounds(lower, upper, precisionBits);
}

[[nodiscard]] bool exactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

[[nodiscard]] std::optional<CertifiedValue> storedNumericInterval(
    const expression::Expr& value,
    std::size_t precisionBits,
    bool semantic) {
    if (value.isNumber()) {
        const auto& number = value.asNumber();
        if (number.isReal())
            return CertifiedValue{RealInterval::fromRational(
                number.asReal().toRational(), precisionBits)};

        const auto& complex = number.asComplex();
        return CertifiedValue{ComplexInterval{
            RealInterval::fromRational(complex.real.toRational(), precisionBits),
            RealInterval::fromRational(complex.imaginary.toRational(), precisionBits)}};
    }

    if (value.isDecimalApproximation())
        return CertifiedValue{semantic
            ? semanticIntervalFromDecimal(value.asDecimalApproximation(), precisionBits)
            : intervalFromDecimal(value.asDecimalApproximation(), precisionBits)};

    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        return CertifiedValue{ComplexInterval{
            semantic ? semanticIntervalFromDecimal(complex.real(), precisionBits)
                     : intervalFromDecimal(complex.real(), precisionBits),
            semantic ? semanticIntervalFromDecimal(complex.imaginary(), precisionBits)
                     : intervalFromDecimal(complex.imaginary(), precisionBits)}};
    }

    return std::nullopt;
}

[[nodiscard]] CertifiedValue normalizeComplex(ComplexInterval value) {
    if (value.isProvablyReal())
        return CertifiedValue{value.real()};
    return CertifiedValue{std::move(value)};
}

[[nodiscard]] CertifiedValue addValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{add(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeComplex(add(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue subtractValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{subtract(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeComplex(subtract(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue multiplyValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{multiply(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeComplex(multiply(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue divideValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{divide(lhs.asReal(), rhs.asReal(), precisionBits)};
    return normalizeComplex(divide(lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue negateValue(const CertifiedValue& value) {
    if (value.isReal())
        return CertifiedValue{negate(value.asReal())};
    return normalizeComplex(negate(value.asComplex()));
}

[[nodiscard]] std::size_t intervalDigitHint(
    const RealInterval& value,
    std::size_t maximumFractionalDigits) {
    if (value.isPoint())
        return maximumFractionalDigits;

    const numeric::Rational width = value.upper().toRational() - value.lower().toRational();
    const std::size_t numeratorDigits = width.numerator().abs().toString().size();
    const std::size_t denominatorDigits = width.denominator().toString().size();
    if (denominatorDigits <= numeratorDigits + 1)
        return 1;

    return std::min(maximumFractionalDigits, denominatorDigits - numeratorDigits - 1);
}

[[nodiscard]] std::size_t intervalDigitHint(
    const CertifiedValue& value,
    std::size_t maximumFractionalDigits) {
    if (value.isReal())
        return intervalDigitHint(value.asReal(), maximumFractionalDigits);

    return std::min(
        intervalDigitHint(value.asComplex().real(), maximumFractionalDigits),
        intervalDigitHint(value.asComplex().imaginary(), maximumFractionalDigits));
}

[[nodiscard]] std::optional<expression::Expr> decimalAtDigits(
    const CertifiedValue& value,
    std::size_t fractionalDigits) {
    return value.isReal()
        ? decimalExpression(value.asReal(), fractionalDigits)
        : decimalExpression(value.asComplex(), fractionalDigits);
}

[[nodiscard]] std::optional<expression::Expr> bestDecimalExpression(
    const CertifiedValue& value,
    std::size_t maximumFractionalDigits) {
    if (const auto decimal = decimalAtDigits(value, maximumFractionalDigits))
        return decimal;
    if (maximumFractionalDigits == 0)
        return std::nullopt;

    // 大きなscale変更でaccuracyが数万桁落ち得るため，要求桁から1桁ずつ降りない。
    // 区間幅の10進桁数から候補へ跳び，境界付近の丸めを考慮して2桁だけ上から確認する。
    const std::size_t hint = intervalDigitHint(value, maximumFractionalDigits);
    std::size_t digits = hint;
    if (digits < maximumFractionalDigits)
        digits = std::min(maximumFractionalDigits - 1, digits + 2);

    for (; digits != 0; --digits)
        if (const auto decimal = decimalAtDigits(value, digits))
            return decimal;
    return decimalAtDigits(value, 0);
}

[[nodiscard]] numeric::Rational absRational(const numeric::Rational& value) {
    return value.numerator().isNegative() ? -value : value;
}

[[nodiscard]] numeric::Rational maximum(
    const numeric::Rational& lhs,
    const numeric::Rational& rhs) {
    return lhs < rhs ? rhs : lhs;
}

[[nodiscard]] std::size_t guaranteedDigits(
    const numeric::Rational& error,
    std::size_t cap) {
    if (error.isZero())
        return cap;
    if (cap == 0 || error >= numeric::Rational{numeric::BigInt{1}})
        return 0;

    const std::size_t numeratorDigits = error.numerator().abs().toString().size();
    const std::size_t denominatorDigits = error.denominator().toString().size();
    if (denominatorDigits <= numeratorDigits)
        return 0;

    std::size_t candidate = std::min(cap, denominatorDigits - numeratorDigits);
    const auto exponent = static_cast<std::uint64_t>(candidate);
    const numeric::BigInt scaledNumerator = error.numerator().abs()
        * numeric::pow(numeric::BigInt{10}, exponent);
    if (scaledNumerator < error.denominator())
        return candidate;
    return candidate == 0 ? 0 : candidate - 1;
}

[[nodiscard]] std::size_t componentAccuracyCap(
    const numeric::DecimalApproximation& displayed,
    const RealInterval& semantic,
    std::size_t cap) {
    const numeric::Rational lowerError = absRational(
        displayed.displayedValue() - semantic.lower().toRational());
    const numeric::Rational upperError = absRational(
        displayed.displayedValue() - semantic.upper().toRational());
    return guaranteedDigits(maximum(lowerError, upperError), cap);
}

[[nodiscard]] std::size_t semanticAccuracyCap(
    const expression::Expr& displayed,
    const CertifiedValue& semantic,
    std::size_t cap) {
    if (displayed.isDecimalApproximation()) {
        const RealInterval& interval = semantic.isReal()
            ? semantic.asReal() : semantic.asComplex().real();
        return componentAccuracyCap(displayed.asDecimalApproximation(), interval, cap);
    }

    if (!displayed.isComplexDecimalApproximation())
        return cap;
    const auto& complex = displayed.asComplexDecimalApproximation();
    const ComplexInterval interval = semantic.toComplex();
    return std::min(
        componentAccuracyCap(complex.real(), interval.real(), cap),
        componentAccuracyCap(complex.imaginary(), interval.imaginary(), cap));
}

[[nodiscard]] std::optional<expression::Expr> finalizeApproximateOperation(
    const CertifiedValue& certified,
    const CertifiedValue& semantic,
    std::size_t sourceDigits) {
    const auto preliminary = bestDecimalExpression(certified, sourceDigits);
    if (!preliminary)
        return std::nullopt;
    const std::size_t cap = semanticAccuracyCap(*preliminary, semantic, sourceDigits);
    return bestDecimalExpression(certified, cap);
}

void collectApproximationDigits(
    const expression::Expr& root,
    std::optional<std::size_t>& digits) {
    std::vector<const expression::Expr*> pending{&root};
    while (!pending.empty()) {
        const expression::Expr& current = *pending.back();
        pending.pop_back();

        std::optional<std::size_t> currentDigits;
        if (current.isDecimalApproximation())
            currentDigits = current.asDecimalApproximation().requestedFractionalDigits();
        else if (current.isComplexDecimalApproximation()) {
            const auto& complex = current.asComplexDecimalApproximation();
            currentDigits = std::min(
                complex.real().requestedFractionalDigits(),
                complex.imaginary().requestedFractionalDigits());
        }

        if (currentDigits && *currentDigits != 0)
            digits = digits ? std::min(*digits, *currentDigits) : currentDigits;

        if (current.isArray()) {
            const auto& array = current.asArray();
            for (std::size_t i = 0; i < array.size(); ++i) {
                switch (array.storedKindAt(i)) {
                case expression::ArrayStorageKind::DecimalApproximation: {
                    const auto valueDigits = array.decimalAt(i).requestedFractionalDigits();
                    if (valueDigits != 0)
                        digits = digits ? std::min(*digits, valueDigits) : valueDigits;
                    break;
                }
                case expression::ArrayStorageKind::ComplexDecimalApproximation: {
                    const auto& value = array.complexDecimalAt(i);
                    const auto valueDigits = std::min(
                        value.real().requestedFractionalDigits(),
                        value.imaginary().requestedFractionalDigits());
                    if (valueDigits != 0)
                        digits = digits ? std::min(*digits, valueDigits) : valueDigits;
                    break;
                }
                case expression::ArrayStorageKind::Generic:
                    pending.push_back(&array.expressionAt(i));
                    break;
                default:
                    break;
                }
            }
        }
        else if (current.isList())
            for (const expression::Expr& element : current.asList().elements)
                pending.push_back(&element);
        else if (current.isCall())
            for (const expression::Expr& argument : current.asCall().arguments)
                pending.push_back(&argument);
    }
}

} // namespace

std::optional<ComplexInterval> encloseComplexExpression(
    const expression::Expr& expression,
    std::size_t precisionBits,
    const CertifiedEvaluator& certified) {
    if (expression.isDecimalApproximation())
        return ComplexInterval::fromReal(
            intervalFromDecimal(expression.asDecimalApproximation(), precisionBits));

    if (expression.isComplexDecimalApproximation()) {
        const auto& value = expression.asComplexDecimalApproximation();
        return ComplexInterval{
            intervalFromDecimal(value.real(), precisionBits),
            intervalFromDecimal(value.imaginary(), precisionBits)};
    }

    const auto enclosed = certified.enclose(expression, precisionBits);
    if (!enclosed)
        return std::nullopt;
    return enclosed->toComplex();
}

std::optional<expression::Expr> decimalExpression(
    const RealInterval& value,
    std::size_t fractionalDigits) {
    if (value.isPoint()) {
        const numeric::RealNumber exact{value.lower().toRational()};
        return expression::Expr{fractionalDigits == 0
            ? numeric::DecimalApproximation::fromRealFixed(exact, 0)
            : numeric::DecimalApproximation::fromReal(exact, fractionalDigits)};
    }

    const auto decimal = numeric::DecimalApproximation::fromCertifiedInterval(
        value.lower().toRational(), value.upper().toRational(), fractionalDigits);
    return decimal ? std::optional<expression::Expr>{expression::Expr{*decimal}}
                   : std::nullopt;
}

std::optional<expression::Expr> decimalExpression(
    const ComplexInterval& value,
    std::size_t fractionalDigits) {
    if (value.real().isPoint() && value.imaginary().isPoint()) {
        const numeric::RealNumber exactReal{value.real().lower().toRational()};
        const numeric::RealNumber exactImaginary{value.imaginary().lower().toRational()};
        const auto real = fractionalDigits == 0
            ? numeric::DecimalApproximation::fromRealFixed(exactReal, 0)
            : numeric::DecimalApproximation::fromReal(exactReal, fractionalDigits);
        const auto imaginary = fractionalDigits == 0
            ? numeric::DecimalApproximation::fromRealFixed(exactImaginary, 0)
            : numeric::DecimalApproximation::fromReal(exactImaginary, fractionalDigits);
        if (exactZero(value.imaginary()))
            return expression::Expr{real};
        return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
            real, imaginary, exactZero(value.real()), false)};
    }

    const auto real = numeric::DecimalApproximation::fromCertifiedInterval(
        value.real().lower().toRational(), value.real().upper().toRational(), fractionalDigits);
    const auto imaginary = numeric::DecimalApproximation::fromCertifiedInterval(
        value.imaginary().lower().toRational(), value.imaginary().upper().toRational(), fractionalDigits);
    if (!real || !imaginary)
        return std::nullopt;
    if (exactZero(value.imaginary()))
        return expression::Expr{*real};
    return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
        *real, *imaginary, exactZero(value.real()), false)};
}

std::optional<ApproximationContext> inferredApproximationContext(
    std::span<const expression::Expr> expressions) {
    std::optional<std::size_t> digits;
    for (const expression::Expr& expression : expressions)
        collectApproximationDigits(expression, digits);
    if (!digits)
        return std::nullopt;
    return ApproximationContext{*digits};
}

std::optional<expression::Expr> addApproximateScalars(
    std::span<const expression::Expr> expressions) {
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    CertifiedValue result{RealInterval::fromRational(numeric::Rational{}, precisionBits)};
    CertifiedValue semantic{RealInterval::fromRational(numeric::Rational{}, precisionBits)};
    for (const auto& expression : expressions) {
        const auto enclosed = storedNumericInterval(expression, precisionBits, false);
        const auto semanticValue = storedNumericInterval(expression, precisionBits, true);
        if (!enclosed || !semanticValue)
            return std::nullopt;
        result = addValues(result, *enclosed, precisionBits);
        semantic = addValues(semantic, *semanticValue, precisionBits);
    }
    return finalizeApproximateOperation(result, semantic, context->decimalDigits());
}

std::optional<expression::Expr> subtractApproximateScalars(
    const expression::Expr& lhs,
    const expression::Expr& rhs) {
    const std::array<expression::Expr, 2> expressions{lhs, rhs};
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    const auto left = storedNumericInterval(lhs, precisionBits, false);
    const auto right = storedNumericInterval(rhs, precisionBits, false);
    const auto semanticLeft = storedNumericInterval(lhs, precisionBits, true);
    const auto semanticRight = storedNumericInterval(rhs, precisionBits, true);
    if (!left || !right || !semanticLeft || !semanticRight)
        return std::nullopt;
    return finalizeApproximateOperation(
        subtractValues(*left, *right, precisionBits),
        subtractValues(*semanticLeft, *semanticRight, precisionBits),
        context->decimalDigits());
}

std::optional<expression::Expr> multiplyApproximateScalars(
    std::span<const expression::Expr> expressions) {
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    CertifiedValue result{RealInterval::fromRational(
        numeric::Rational{numeric::BigInt{1}}, precisionBits)};
    CertifiedValue semantic{RealInterval::fromRational(
        numeric::Rational{numeric::BigInt{1}}, precisionBits)};
    for (const auto& expression : expressions) {
        const auto enclosed = storedNumericInterval(expression, precisionBits, false);
        const auto semanticValue = storedNumericInterval(expression, precisionBits, true);
        if (!enclosed || !semanticValue)
            return std::nullopt;
        result = multiplyValues(result, *enclosed, precisionBits);
        semantic = multiplyValues(semantic, *semanticValue, precisionBits);
    }
    return finalizeApproximateOperation(result, semantic, context->decimalDigits());
}

std::optional<expression::Expr> divideApproximateScalars(
    const expression::Expr& lhs,
    const expression::Expr& rhs) {
    const std::array<expression::Expr, 2> expressions{lhs, rhs};
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    const auto left = storedNumericInterval(lhs, precisionBits, false);
    const auto right = storedNumericInterval(rhs, precisionBits, false);
    const auto semanticLeft = storedNumericInterval(lhs, precisionBits, true);
    const auto semanticRight = storedNumericInterval(rhs, precisionBits, true);
    if (!left || !right || !semanticLeft || !semanticRight)
        return std::nullopt;

    try {
        return finalizeApproximateOperation(
            divideValues(*left, *right, precisionBits),
            divideValues(*semanticLeft, *semanticRight, precisionBits),
            context->decimalDigits());
    }
    catch (const std::domain_error&) {
        return std::nullopt;
    }
}

std::optional<expression::Expr> negateApproximateScalar(
    const expression::Expr& value) {
    const std::array<expression::Expr, 1> expressions{value};
    const auto context = inferredApproximationContext(expressions);
    if (!context)
        return std::nullopt;

    const std::size_t precisionBits = context->workingBinaryBits();
    const auto enclosed = storedNumericInterval(value, precisionBits, false);
    const auto semanticValue = storedNumericInterval(value, precisionBits, true);
    if (!enclosed || !semanticValue)
        return std::nullopt;
    return finalizeApproximateOperation(
        negateValue(*enclosed), negateValue(*semanticValue), context->decimalDigits());
}

std::size_t nextGuardDigits(std::size_t current) {
    const std::size_t growth = std::max<std::size_t>(8, current / 2);
    if (growth > std::numeric_limits<std::size_t>::max() - current)
        throw std::overflow_error("Approximation precision is too large");
    return current + growth;
}

} // namespace mmcal::approximation
