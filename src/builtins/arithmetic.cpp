// 四則演算、冪
#include "arithmetic.hpp"

#include "approximation/expression_interval.hpp"
#include "error/error_message.hpp"
#include "names.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/integer_algorithms.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] Expr numberExpr(Number number) {
    return Expr{std::move(number)};
}

[[nodiscard]] Expr integerExpr(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

void requireArity(std::span<const Expr> arguments, std::size_t expected, std::string_view name) {
    if (arguments.size() != expected)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects " + std::to_string(expected) + " argument(s)");
}

[[nodiscard]] std::optional<std::uint64_t> toUint64(const BigInt& value) {
    if (value.isNegative())
        return std::nullopt;

    const std::string text = value.toString();
    std::uint64_t result = 0;
    const auto conversion = std::from_chars(text.data(), text.data() + text.size(), result);
    if (conversion.ec != std::errc{} || conversion.ptr != text.data() + text.size())
        return std::nullopt;

    return result;
}


[[nodiscard]] Expr scalarAdd(
    const std::vector<Expr>& arguments,
    const evaluation::BuiltinRegistry& registry) {
    Number sum{BigInt{0}};
    bool numeric = true;
    for (const Expr& argument : arguments) {
        if (!argument.isNumber()) {
            numeric = false;
            break;
        }
        sum += argument.asNumber();
    }
    if (numeric)
        return numberExpr(std::move(sum));
    if (const auto approximate = approximation::addApproximateScalars(arguments))
        return *approximate;
    return Expr::call(registry.symbol(evaluation::BuiltinId::Add), arguments);
}

[[nodiscard]] Expr scalarMultiply(
    const std::vector<Expr>& arguments,
    const evaluation::BuiltinRegistry& registry) {
    Number product{BigInt{1}};
    bool numeric = true;
    for (const Expr& argument : arguments) {
        if (!argument.isNumber()) {
            numeric = false;
            break;
        }
        product *= argument.asNumber();
    }
    if (numeric)
        return numberExpr(std::move(product));
    if (const auto approximate = approximation::multiplyApproximateScalars(arguments))
        return *approximate;
    return Expr::call(registry.symbol(evaluation::BuiltinId::Multiply), arguments);
}

[[nodiscard]] Expr scalarSubtract(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& registry) {
    if (lhs.isNumber() && rhs.isNumber())
        return numberExpr(lhs.asNumber() - rhs.asNumber());
    if (const auto approximate = approximation::subtractApproximateScalars(lhs, rhs))
        return *approximate;
    return Expr::call(registry.symbol(evaluation::BuiltinId::Subtract), {lhs, rhs});
}

[[nodiscard]] Expr scalarNegate(
    const Expr& value,
    const evaluation::BuiltinRegistry& registry) {
    if (value.isNumber())
        return numberExpr(-value.asNumber());
    if (const auto approximate = approximation::negateApproximateScalar(value))
        return *approximate;
    return Expr::call(registry.symbol(evaluation::BuiltinId::Negate), {value});
}

[[nodiscard]] bool isCertifiedExactZero(const Expr& value) noexcept {
    if (value.isNumber())
        return value.asNumber().isZero();
    if (value.isDecimalApproximation()) {
        const auto& decimal = value.asDecimalApproximation();
        return decimal.enclosureIsPoint() && decimal.certifiedLower().isZero();
    }
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        return complex.real().enclosureIsPoint() && complex.real().certifiedLower().isZero()
            && complex.imaginary().enclosureIsPoint()
            && complex.imaginary().certifiedLower().isZero();
    }
    return false;
}

[[noreturn]] void arrayArithmeticError(std::string message) {
    error::throwCalcError(error::CalcErrorType::Type, std::move(message));
}

} // namespace

Expr evaluateAdd(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    const bool hasArray = std::any_of(arguments.begin(), arguments.end(),
        [](const Expr& argument) { return argument.isArray(); });
    if (!hasArray)
        return scalarAdd(std::vector<Expr>{arguments.begin(), arguments.end()}, registry);

    for (const Expr& argument : arguments)
        if (!argument.isArray())
            arrayArithmeticError("Array addition requires arrays with identical shapes");

    const auto& first = arguments.front().asArray();
    for (const Expr& argument : arguments)
        if (argument.asArray().shape != first.shape)
            error::throwCalcError(error::CalcErrorType::Domain,
                "Array addition requires identical shapes");

    const bool allExact = std::all_of(arguments.begin(), arguments.end(),
        [](const Expr& argument) { return argument.asArray().hasExactNumberStorage(); });
    if (allExact) {
        std::vector<Number> values(first.size(), Number{BigInt{0}});
        for (const Expr& argument : arguments) {
            const auto& array = argument.asArray();
            for (std::size_t i = 0; i < values.size(); ++i)
                values[i] += array.exactNumber(i);
        }
        return Expr::numberArray(first.shape, std::move(values));
    }

    std::vector<Expr> elements;
    elements.reserve(first.size());
    for (std::size_t i = 0; i < first.size(); ++i) {
        std::vector<Expr> terms;
        terms.reserve(arguments.size());
        for (const Expr& argument : arguments)
            terms.push_back(argument.asArray().element(i));
        elements.push_back(scalarAdd(terms, registry));
    }
    return Expr::array(first.shape, std::move(elements));
}

Expr evaluateSubtract(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::subtract);
    const Expr& lhs = arguments[0];
    const Expr& rhs = arguments[1];

    if (lhs.isArray() || rhs.isArray()) {
        if (!lhs.isArray() || !rhs.isArray())
            arrayArithmeticError("Array subtraction requires two arrays with identical shapes");
        if (lhs.asArray().shape != rhs.asArray().shape)
            error::throwCalcError(error::CalcErrorType::Domain,
                "Array subtraction requires identical shapes");

        if (lhs.asArray().hasExactNumberStorage() && rhs.asArray().hasExactNumberStorage()) {
            std::vector<Number> values;
            values.reserve(lhs.asArray().size());
            for (std::size_t i = 0; i < lhs.asArray().size(); ++i)
                values.push_back(lhs.asArray().exactNumber(i) - rhs.asArray().exactNumber(i));
            return Expr::numberArray(lhs.asArray().shape, std::move(values));
        }

        std::vector<Expr> elements;
        elements.reserve(lhs.asArray().size());
        for (std::size_t i = 0; i < lhs.asArray().size(); ++i) {
            const Expr left = lhs.asArray().element(i);
            const Expr right = rhs.asArray().element(i);
            elements.push_back(scalarSubtract(left, right, registry));
        }
        return Expr::array(lhs.asArray().shape, std::move(elements));
    }

    return scalarSubtract(lhs, rhs, registry);
}

Expr evaluateMultiply(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    std::size_t arrayCount = 0;
    std::size_t arrayIndex = 0;
    for (std::size_t i = 0; i < arguments.size(); ++i)
        if (arguments[i].isArray()) {
            ++arrayCount;
            arrayIndex = i;
        }

    if (arrayCount == 0)
        return scalarMultiply(std::vector<Expr>{arguments.begin(), arguments.end()}, registry);
    if (arrayCount > 1)
        arrayArithmeticError(
            "Array multiplication is scalar-only; use dot[...] for vector or matrix contraction");

    const auto& array = arguments[arrayIndex].asArray();
    std::vector<Expr> scalarFactors;
    scalarFactors.reserve(arguments.size() - 1);
    for (std::size_t i = 0; i < arguments.size(); ++i)
        if (i != arrayIndex)
            scalarFactors.push_back(arguments[i]);

    if (array.hasExactNumberStorage()
        && std::all_of(scalarFactors.begin(), scalarFactors.end(),
            [](const Expr& value) { return value.isNumber(); })) {
        Number scalar{BigInt{1}};
        for (const Expr& factor : scalarFactors)
            scalar *= factor.asNumber();
        std::vector<Number> values;
        values.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i)
            values.push_back(array.exactNumber(i) * scalar);
        return Expr::numberArray(array.shape, std::move(values));
    }

    std::vector<Expr> elements;
    elements.reserve(array.size());
    for (std::size_t i = 0; i < array.size(); ++i) {
        const Expr element = array.element(i);
        std::vector<Expr> factors;
        factors.reserve(scalarFactors.size() + 1);
        factors.push_back(element);
        factors.insert(factors.end(), scalarFactors.begin(), scalarFactors.end());
        elements.push_back(scalarMultiply(factors, registry));
    }
    return Expr::array(array.shape, std::move(elements));
}

Expr evaluateDivide(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::divide);

    const Expr& numerator = arguments[0];
    const Expr& denominator = arguments[1];
    if (isCertifiedExactZero(denominator))
        error::throwCalcError(error::CalcErrorType::Domain, "Division by zero");
    if (numerator.isNumber() && denominator.isNumber())
        return numberExpr(numerator.asNumber() / denominator.asNumber());
    if (const auto approximate = approximation::divideApproximateScalars(
        numerator, denominator))
        return *approximate;

    return Expr::call(
        registry.symbol(evaluation::BuiltinId::Divide),
        {numerator, denominator});
}

Expr evaluatePower(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireArity(arguments, 2, names::power);

    const Expr& base = arguments[0];
    const Expr& exponent = arguments[1];

    if (base.isNumber() && base.asNumber().isZero() && exponent.isNumber()
        && exponent.asNumber().isReal()) {
        const RealNumber& realExponent = exponent.asNumber().asReal();
        if (realExponent.isZero())
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "Zero to the zero power is indeterminate");
        if (realExponent.isNegative())
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "Zero cannot be raised to a negative power");

        return integerExpr(0);
    }

    if (base.isNumber() && base.asNumber() == Number{BigInt{1}})
        return integerExpr(1);

    // mmCalでは有限小数もexact Rationalとして読むため、0.5 は厳密に 1/2。
    // Power[x, 1/2] はprincipal square rootの定義と一致させる。これにより
    //     (-2)^0.5  -> I * sqrt[2]
    //     4^0.5     -> 2
    //     x^0.5     -> sqrt[x]
    // となり、「小数で書いたからmachine-real powerへ落ちる」という別意味を作らない。
    if (exponent.isNumber() && exponent.asNumber().isReal()
        && exponent.asNumber().asReal().toRational()
            == Rational{BigInt{1}, BigInt{2}}) {
        return evaluateSqrt(std::span<const Expr>{&base, 1}, registry, mathematics);
    }

    if (!exponent.isNumber() || !exponent.asNumber().isReal()
        || !exponent.asNumber().asReal().isInteger())
        return Expr::call(registry.symbol(evaluation::BuiltinId::Power), {base, exponent});

    const BigInt& integerExponent = exponent.asNumber().asReal().asInteger();
    if (integerExponent.isZero()) {
        // 0^0をDomainErrorとしている以上、未知のsymbolic baseに対してx^0 -> 1 と無条件簡約するのは安全ではない。
        // baseが非零と証明できる場合だけ1へ畳み込み、未知ならPower式を保持する。
        if (base.isNumber())
            return integerExpr(1);

        const mathematics::ValueFacts baseFacts = mathematics::inferValueFacts(
            base, registry, mathematics);
        if (baseFacts.sign == mathematics::RealSign::Positive
            || baseFacts.sign == mathematics::RealSign::Negative
            || baseFacts.sign == mathematics::RealSign::NonZero)
            return integerExpr(1);

        return Expr::call(registry.symbol(evaluation::BuiltinId::Power), {base, exponent});
    }
    if (integerExponent == BigInt{1})
        return base;

    if (!base.isNumber())
        return Expr::call(registry.symbol(evaluation::BuiltinId::Power), {base, exponent});

    const bool negativeExponent = integerExponent.isNegative();
    const auto magnitude = toUint64(integerExponent.abs());
    if (!magnitude)
        error::throwCalcError(
            error::CalcErrorType::Overflow,
            "Exponent is too large for exact evaluation");

    if (negativeExponent && base.asNumber().isZero())
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "Zero cannot be raised to a negative power");

    Number result = numeric::integerPower(base.asNumber(), *magnitude);
    if (negativeExponent)
        result = Number{BigInt{1}} / result;

    return numberExpr(std::move(result));
}

Expr evaluateNegate(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::negate);
    if (arguments.front().isArray()) {
        const auto& array = arguments.front().asArray();
        if (array.hasExactNumberStorage()) {
            std::vector<Number> values;
            values.reserve(array.size());
            for (std::size_t i = 0; i < array.size(); ++i)
                values.push_back(-array.exactNumber(i));
            return Expr::numberArray(array.shape, std::move(values));
        }
        std::vector<Expr> elements;
        elements.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i) {
            elements.push_back(scalarNegate(array.element(i), registry));
        }
        return Expr::array(array.shape, std::move(elements));
    }
    return scalarNegate(arguments.front(), registry);
}

Expr evaluateFactorial(std::span<const Expr> arguments) {
    requireArity(arguments, 1, names::factorial);

    const Expr& argument = arguments.front();
    if (!argument.isNumber() || !argument.asNumber().isReal()
        || !argument.asNumber().asReal().isInteger())
        error::throwCalcError(
            error::CalcErrorType::Type,
            "Factorial requires a non-negative integer");

    const BigInt& value = argument.asNumber().asReal().asInteger();
    if (value.isNegative())
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "Factorial is undefined for negative integers");

    const auto count = toUint64(value);
    if (!count)
        error::throwCalcError(
            error::CalcErrorType::Overflow,
            "Factorial argument is too large for exact evaluation");

    return Expr{Number{numeric::factorial(*count)}};
}

Expr evaluateSqrt(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireArity(arguments, 1, names::sqrt);
    static_cast<void>(mathematics);

    // principal sqrtのexact/conditional簡約は共通Simplifierに集約する。
    return Expr::call(
        registry.symbol(evaluation::BuiltinId::Sqrt),
        {arguments.front()});
}

} // namespace mmcal::builtins
