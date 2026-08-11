// 四則演算、冪
#include "arithmetic.hpp"

#include "error/error_message.hpp"
#include "names.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/integer_algorithms.hpp"

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

} // namespace

Expr evaluateAdd(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    bool allNumeric = true;
    Number sum{BigInt{0}};
    for (const Expr& argument : arguments) {
        if (!argument.isNumber()) {
            allNumeric = false;
            break;
        }
        sum += argument.asNumber();
    }

    if (allNumeric)
        return numberExpr(std::move(sum));

    // 記号項のflatten、係数収集、0除去は共通Simplifierの責務。
    return Expr::call(
        registry.symbol(evaluation::BuiltinId::Add),
        std::vector<Expr>{arguments.begin(), arguments.end()});
}

Expr evaluateSubtract(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::subtract);
    if (arguments[0].isNumber() && arguments[1].isNumber())
        return numberExpr(arguments[0].asNumber() - arguments[1].asNumber());
    return Expr::call(
        registry.symbol(evaluation::BuiltinId::Subtract),
        {arguments[0], arguments[1]});
}

Expr evaluateMultiply(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    bool allNumeric = true;
    Number product{BigInt{1}};
    for (const Expr& argument : arguments) {
        if (!argument.isNumber()) {
            allNumeric = false;
            break;
        }
        product *= argument.asNumber();
    }

    if (allNumeric)
        return numberExpr(std::move(product));

    // 記号因子のflatten、0/1除去、係数整理は共通Simplifierへ集約する。
    return Expr::call(
        registry.symbol(evaluation::BuiltinId::Multiply),
        std::vector<Expr>{arguments.begin(), arguments.end()});
}

Expr evaluateDivide(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::divide);

    const Expr& numerator = arguments[0];
    const Expr& denominator = arguments[1];
    if (denominator.isNumber() && denominator.asNumber().isZero())
        error::throwCalcError(error::CalcErrorType::Domain, "Division by zero");
    if (numerator.isNumber() && denominator.isNumber())
        return numberExpr(numerator.asNumber() / denominator.asNumber());

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
    if (arguments.front().isNumber())
        return numberExpr(-arguments.front().asNumber());
    return Expr::call(
        registry.symbol(evaluation::BuiltinId::Negate),
        {arguments.front()});
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
