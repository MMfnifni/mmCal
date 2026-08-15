// precision・accuracy・rationalize
#include "approximation_utilities.hpp"
#include "expression/array_utils.hpp"

#include "error/error_message.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/number.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using expression::Expr;
using numeric::BigInt;
using numeric::DecimalApproximation;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] Rational absRational(const Rational& value) {
    return value.numerator().isNegative() ? -value : value;
}

[[nodiscard]] BigInt floorRational(const Rational& value) {
    BigInt quotient = value.numerator() / value.denominator();
    const BigInt remainder = value.numerator() % value.denominator();
    if (value.numerator().isNegative() && !remainder.isZero())
        quotient -= BigInt{1};
    return quotient;
}

[[nodiscard]] Rational maximum(const Rational& lhs, const Rational& rhs) {
    return lhs < rhs ? rhs : lhs;
}

[[nodiscard]] Rational informationAbsoluteError(const DecimalApproximation& value) {
    const Rational lowerError = absRational(value.displayedValue() - value.informationLower());
    const Rational upperError = absRational(value.displayedValue() - value.informationUpper());
    return maximum(lowerError, upperError);
}

[[nodiscard]] std::size_t decimalIntegerDigits(const Rational& value) {
    const BigInt integerPart = value.numerator().abs() / value.denominator();
    if (integerPart.isZero())
        return 1;
    return integerPart.toString().size();
}

[[nodiscard]] std::size_t guaranteedDigits(const Rational& error, std::size_t cap) {
    if (error.isZero())
        return cap;

    Rational threshold{BigInt{1}};
    std::size_t digits = 0;
    for (std::size_t next = 1; next <= cap; ++next) {
        threshold /= Rational{BigInt{10}};
        if (!(error < threshold))
            break;
        digits = next;
    }
    return digits;
}

[[nodiscard]] std::size_t accuracyDigits(const DecimalApproximation& value) {
    return guaranteedDigits(
        informationAbsoluteError(value),
        value.requestedFractionalDigits());
}

[[nodiscard]] std::size_t precisionDigits(const DecimalApproximation& value) {
    const Rational& lower = value.informationLower();
    const Rational& upper = value.informationUpper();
    if (lower <= Rational{} && upper >= Rational{})
        return 0;

    const Rational minimumMagnitude = lower > Rational{}
        ? lower
        : -upper;
    if (minimumMagnitude.isZero())
        return 0;

    const Rational relativeError = informationAbsoluteError(value) / minimumMagnitude;
    const std::size_t cap = value.requestedFractionalDigits()
        + decimalIntegerDigits(value.displayedValue()) + 2;
    return guaranteedDigits(relativeError, cap);
}

[[nodiscard]] Expr integerResult(std::size_t value) {
    return Expr{Number{BigInt::parse(std::to_string(value))}};
}

[[nodiscard]] bool containsApproximation(const Expr& root) {
    std::vector<Expr> stack{root};
    while (!stack.empty()) {
        const Expr current = stack.back();
        stack.pop_back();
        if (current.isDecimalApproximation() || current.isComplexDecimalApproximation())
            return true;
        if (current.isCall())
            for (const Expr& argument : current.asCall().arguments)
                stack.push_back(argument);
        else if (current.isArray()) {
            const auto& array = current.asArray();
            for (std::size_t i = 0; i < array.size(); ++i) {
                const auto kind = array.storedKindAt(i);
                if (kind == expression::ArrayStorageKind::DecimalApproximation
                    || kind == expression::ArrayStorageKind::ComplexDecimalApproximation)
                    return true;
                if (kind == expression::ArrayStorageKind::Generic)
                    stack.push_back(array.expressionAt(i));
            }
        }
        else if (current.isList())
            for (const Expr& element : current.asList().elements)
                stack.push_back(element);
    }
    return false;
}

[[nodiscard]] bool isExactNumericLike(const Expr& value) {
    return value.isNumber() || value.isSymbol() || value.isCall()
        || value.isArray() || value.isList();
}

[[nodiscard]] std::optional<Expr> accuracyOf(const Expr& value, const expression::Symbol& infinity) {
    if (value.isNumber())
        return Expr{infinity};
    if (value.isDecimalApproximation())
        return integerResult(accuracyDigits(value.asDecimalApproximation()));
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        return integerResult(std::min(
            accuracyDigits(complex.real()),
            accuracyDigits(complex.imaginary())));
    }
    if (isExactNumericLike(value) && !containsApproximation(value))
        return Expr{infinity};
    return {};
}

[[nodiscard]] std::optional<Expr> precisionOf(const Expr& value, const expression::Symbol& infinity) {
    if (value.isNumber())
        return Expr{infinity};
    if (value.isDecimalApproximation())
        return integerResult(precisionDigits(value.asDecimalApproximation()));
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        return integerResult(std::min(
            precisionDigits(complex.real()),
            precisionDigits(complex.imaginary())));
    }
    if (isExactNumericLike(value) && !containsApproximation(value))
        return Expr{infinity};
    return {};
}

[[nodiscard]] Rational simplestPositive(const Rational& lower, const Rational& upper) {
    if (!(Rational{} < lower) || upper < lower)
        error::throwCalcError(error::CalcErrorType::Internal, "Invalid positive rational interval");

    if (lower.isInteger())
        return lower;

    const BigInt lowerFloor = floorRational(lower);
    const BigInt upperFloor = floorRational(upper);
    if (lowerFloor != upperFloor)
        return Rational{lowerFloor + BigInt{1}};

    const Rational integerPart{lowerFloor};
    const Rational lowFraction = lower - integerPart;
    const Rational highFraction = upper - integerPart;
    if (lowFraction.isZero())
        return lower;

    const Rational reciprocalLow = Rational{BigInt{1}} / highFraction;
    const Rational reciprocalHigh = Rational{BigInt{1}} / lowFraction;
    const Rational inner = simplestPositive(reciprocalLow, reciprocalHigh);
    return integerPart + Rational{BigInt{1}} / inner;
}

[[nodiscard]] Rational simplestInInterval(Rational lower, Rational upper) {
    if (upper < lower)
        std::swap(lower, upper);
    if (lower <= Rational{} && upper >= Rational{})
        return Rational{};
    if (upper < Rational{})
        return -simplestPositive(-upper, -lower);
    return simplestPositive(lower, upper);
}

[[nodiscard]] std::optional<Rational> exactNonnegativeTolerance(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    const Rational value = expression.asNumber().asReal().toRational();
    if (value < Rational{})
        return std::nullopt;
    return value;
}

[[nodiscard]] Rational rationalizeDecimal(
    const DecimalApproximation& value,
    const std::optional<Rational>& tolerance) {
    if (tolerance) {
        if (tolerance->isZero())
            return value.displayedValue();
        return simplestInInterval(
            value.displayedValue() - *tolerance,
            value.displayedValue() + *tolerance);
    }
    return simplestInInterval(value.informationLower(), value.informationUpper());
}

[[nodiscard]] std::optional<Expr> rationalizeValue(
    const Expr& value,
    const std::optional<Rational>& tolerance) {
    if (value.isNumber())
        return value;
    if (value.isDecimalApproximation())
        return Expr{Number{rationalizeDecimal(value.asDecimalApproximation(), tolerance)}};
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        const Rational real = rationalizeDecimal(complex.real(), tolerance);
        const Rational imaginary = rationalizeDecimal(complex.imaginary(), tolerance);
        return Expr{Number::complex(RealNumber{real}, RealNumber{imaginary})};
    }
    if (value.isArray()) {
        const auto& array = value.asArray();
        if (array.hasExactNumberStorage())
            return value;
        if (array.storageKind() == expression::ArrayStorageKind::DecimalApproximation) {
            std::vector<Rational> elements;
            elements.reserve(array.size());
            for (std::size_t i = 0; i < array.size(); ++i)
                elements.push_back(rationalizeDecimal(array.decimalAt(i), tolerance));
            return Expr::rationalArray(array.shape, std::move(elements));
        }
        if (array.storageKind() == expression::ArrayStorageKind::ComplexDecimalApproximation) {
            std::vector<Number> elements;
            elements.reserve(array.size());
            for (std::size_t i = 0; i < array.size(); ++i) {
                const auto& element = array.complexDecimalAt(i);
                const Rational real = rationalizeDecimal(element.real(), tolerance);
                const Rational imaginary = rationalizeDecimal(element.imaginary(), tolerance);
                elements.push_back(Number::complex(RealNumber{real}, RealNumber{imaginary}));
            }
            return Expr::numberArray(array.shape, std::move(elements));
        }
        expression::ArrayBuilder builder;
        builder.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i) {
            const Expr element = array.element(i);
            if (const auto rationalized = rationalizeValue(element, tolerance))
                builder.append(*rationalized);
            else
                builder.append(element);
        }
        return Expr::array(builder.finish(array.shape));
    }
    if (value.isList()) {
        const auto& list = value.asList();
        std::vector<Expr> elements;
        elements.reserve(list.elements.size());
        for (const Expr& element : list.elements)
            if (const auto rationalized = rationalizeValue(element, tolerance))
                elements.push_back(*rationalized);
            else
                elements.push_back(element);
        return expression::braceValue(std::move(elements));
    }
    if (value.isCall()) {
        std::vector<Expr> arguments;
        arguments.reserve(value.asCall().arguments.size());
        for (const Expr& argument : value.asCall().arguments)
            if (const auto rationalized = rationalizeValue(argument, tolerance))
                arguments.push_back(*rationalized);
            else
                arguments.push_back(argument);
        return Expr::call(value.asCall().head, std::move(arguments));
    }
    if (value.isSymbol())
        return value;
    return {};
}

} // namespace

std::optional<expression::Expr> evaluatePrecision(
    std::span<const expression::Expr> arguments,
    const expression::Symbol& infinity) {
    return precisionOf(arguments.front(), infinity);
}

std::optional<expression::Expr> evaluateAccuracy(
    std::span<const expression::Expr> arguments,
    const expression::Symbol& infinity) {
    return accuracyOf(arguments.front(), infinity);
}

std::optional<expression::Expr> evaluateRationalize(
    std::span<const expression::Expr> arguments) {
    std::optional<Rational> tolerance;
    if (arguments.size() == 2) {
        tolerance = exactNonnegativeTolerance(arguments[1]);
        if (!tolerance)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "rationalize tolerance must be a nonnegative exact real number");
    }

    return rationalizeValue(arguments.front(), tolerance);
}

} // namespace mmcal::builtins
