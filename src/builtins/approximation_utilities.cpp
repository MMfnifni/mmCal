// precision・accuracy・rationalize
#include "approximation_utilities.hpp"
#include "expression/array_utils.hpp"

#include "error/error_message.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/approximation_quality.hpp"
#include "numeric/number.hpp"
#include "numeric/rational_rounding.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
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
        return integerResult(numeric::accuracyDigits(value.asDecimalApproximation()));
    if (value.isComplexDecimalApproximation())
        return integerResult(numeric::accuracyDigits(value.asComplexDecimalApproximation()));
    if (isExactNumericLike(value) && !containsApproximation(value))
        return Expr{infinity};
    return {};
}

[[nodiscard]] std::optional<Expr> precisionOf(const Expr& value, const expression::Symbol& infinity) {
    if (value.isNumber())
        return Expr{infinity};
    if (value.isDecimalApproximation())
        return integerResult(numeric::precisionDigits(value.asDecimalApproximation()));
    if (value.isComplexDecimalApproximation())
        return integerResult(numeric::precisionDigits(value.asComplexDecimalApproximation()));
    if (isExactNumericLike(value) && !containsApproximation(value))
        return Expr{infinity};
    return {};
}

[[nodiscard]] Rational simplestPositive(const Rational& lower, const Rational& upper) {
    if (!(Rational{} < lower) || upper < lower)
        error::throwCalcError(error::CalcErrorType::Internal, "Invalid positive rational interval");

    if (lower.isInteger())
        return lower;

    const BigInt lowerFloor = numeric::floorToInteger(lower);
    const BigInt upperFloor = numeric::floorToInteger(upper);
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
        return Expr::rebuildCall(value.asCall(), std::move(arguments));
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
