// 丸め・剰余・整数演算
#include "discrete_math.hpp"
#include "expression/exact_value.hpp"

#include "builtins/names.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/certification_error.hpp"
#include "error/error_message.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"

#include <compare>
#include <algorithm>
#include <cstdint>
#include <charconv>
#include <limits>
#include <numeric>
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

void requireArity(std::span<const Expr> arguments, std::size_t arity, std::string_view name) {
    if (arguments.size() != arity)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects " + std::to_string(arity) + " argument(s)");
}

[[nodiscard]] Expr holdUnary(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    evaluation::BuiltinId id,
    std::string_view name) {
    requireArity(arguments, 1, name);
    return Expr::call(registry.symbol(id), {arguments.front()});
}

[[nodiscard]] const BigInt* exactInteger(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return nullptr;
    return &expression.asNumber().asReal().asInteger();
}

[[nodiscard]] BigInt truncRational(const Rational& value) {
    return value.numerator() / value.denominator();
}

[[nodiscard]] BigInt floorRational(const Rational& value) {
    auto result = numeric::divmod(value.numerator(), value.denominator());
    if (!result.remainder.isZero() && value.numerator().isNegative())
        result.quotient -= BigInt{1};
    return result.quotient;
}

[[nodiscard]] BigInt ceilRational(const Rational& value) {
    auto result = numeric::divmod(value.numerator(), value.denominator());
    if (!result.remainder.isZero() && value.numerator().isPositive())
        result.quotient += BigInt{1};
    return result.quotient;
}

[[nodiscard]] BigInt roundNearestEven(const Rational& value) {
    auto result = numeric::divmod(value.numerator(), value.denominator());
    if (result.remainder.isZero())
        return result.quotient;

    const BigInt twiceRemainder = result.remainder.abs() * BigInt{2};
    const auto comparison = twiceRemainder <=> value.denominator();
    bool awayFromZero = comparison == std::strong_ordering::greater;
    if (comparison == std::strong_ordering::equal) {
        const BigInt parity = result.quotient.abs() % BigInt{2};
        awayFromZero = !parity.isZero();
    }

    if (awayFromZero)
        result.quotient += value.numerator().isNegative() ? BigInt{-1} : BigInt{1};
    return result.quotient;
}

[[nodiscard]] Expr integerResult(BigInt value) {
    return Expr{Number{std::move(value)}};
}


[[nodiscard]] std::int64_t requireSignedSmallInteger(
    const Expr& expression, std::string_view name, std::int64_t limit = 100'000) {
    const BigInt* value = exactInteger(expression);
    if (!value)
        error::throwCalcError(error::CalcErrorType::Type,
            std::string{name} + " requires an exact integer argument");
    const std::string text = value->toString();
    std::int64_t result = 0;
    const auto converted = std::from_chars(text.data(), text.data() + text.size(), result);
    if (converted.ec != std::errc{} || converted.ptr != text.data() + text.size()
        || result < -limit || result > limit)
        error::throwCalcError(error::CalcErrorType::Domain,
            std::string{name} + " integer argument exceeds the current limit");
    return result;
}

[[nodiscard]] Rational roundDecimal(const Rational& value, std::int64_t digits) {
    const std::size_t magnitude = static_cast<std::size_t>(digits < 0 ? -digits : digits);
    const BigInt scale = numeric::pow(BigInt{10}, static_cast<std::uint64_t>(magnitude));
    if (digits >= 0) {
        const Rational scaled = value * Rational{scale};
        return Rational{roundNearestEven(scaled), scale};
    }
    const Rational scaled = value / Rational{scale};
    return Rational{roundNearestEven(scaled) * scale};
}

[[nodiscard]] std::optional<Rational> certifyDecimalRound(
    const Expr& expression,
    std::int64_t digits,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const approximation::CertifiedEvaluator evaluator{registry, mathematics, angles};
    for (std::size_t bits = 64; bits <= 4096; bits *= 2) {
        try {
            const auto enclosed = evaluator.enclose(
                expression, bits,
                approximation::CertifiedEvaluator::EnclosureKind::Information);
            if (!enclosed)
                return std::nullopt;
            const approximation::RealInterval* real = nullptr;
            if (enclosed->isReal())
                real = &enclosed->asReal();
            else if (enclosed->asComplex().isProvablyReal())
                real = &enclosed->asComplex().real();
            else
                error::throwCalcError(error::CalcErrorType::Type, "round requires a real argument");
            const Rational lower = roundDecimal(real->lower().toRational(), digits);
            const Rational upper = roundDecimal(real->upper().toRational(), digits);
            if (lower == upper)
                return lower;
        }
        catch (const approximation::PrecisionInsufficient&) {}
        catch (const approximation::CertifiedBackendUnsupported&) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::size_t requireBitIndex(const Expr& expression, std::string_view name) {
    const BigInt* value = exactInteger(expression);
    if (!value || value->isNegative())
        error::throwCalcError(error::CalcErrorType::Type,
            std::string{name} + " requires a nonnegative integer bit index");
    const auto converted = numeric::tryToUint64(*value);
    if (!converted || *converted > std::numeric_limits<std::size_t>::max())
        error::throwCalcError(error::CalcErrorType::Overflow,
            std::string{name} + " bit index is too large");
    return static_cast<std::size_t>(*converted);
}

[[nodiscard]] BigInt arithmeticShiftRight(BigInt value, std::size_t bits) {
    if (!value.isNegative())
        return value >> bits;
    if (value.isZero())
        return value;
    BigInt magnitude = value.abs();
    if (bits >= magnitude.bitLength())
        return BigInt{-1};
    const BigInt bias = (BigInt{1} << bits) - BigInt{1};
    magnitude += bias;
    magnitude >>= bits;
    return -magnitude;
}

enum class IntegralRoundingOperation {
    Floor,
    Ceil,
    Trunc,
    Round
};

[[nodiscard]] BigInt applyIntegralRounding(
    const Rational& value,
    IntegralRoundingOperation operation) {
    switch (operation) {
    case IntegralRoundingOperation::Floor: return floorRational(value);
    case IntegralRoundingOperation::Ceil: return ceilRational(value);
    case IntegralRoundingOperation::Trunc: return truncRational(value);
    case IntegralRoundingOperation::Round: return roundNearestEven(value);
    }
    return BigInt{};
}

[[nodiscard]] std::optional<BigInt> certifyIntegralRounding(
    const Expr& expression,
    IntegralRoundingOperation operation,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const approximation::CertifiedEvaluator evaluator{registry, mathematics, angles};
    for (std::size_t bits = 64; bits <= 4096; bits *= 2) {
        try {
            const auto enclosed = evaluator.enclose(
                expression, bits,
                approximation::CertifiedEvaluator::EnclosureKind::Information);
            if (!enclosed)
                return std::nullopt;

            const approximation::RealInterval* real = nullptr;
            if (enclosed->isReal())
                real = &enclosed->asReal();
            else if (enclosed->asComplex().isProvablyReal())
                real = &enclosed->asComplex().real();
            else
                error::throwCalcError(
                    error::CalcErrorType::Type,
                    "integer-part function requires a real argument");

            const BigInt lower = applyIntegralRounding(
                real->lower().toRational(), operation);
            const BigInt upper = applyIntegralRounding(
                real->upper().toRational(), operation);
            if (lower == upper)
                return lower;
        }
        catch (const approximation::PrecisionInsufficient&) {
            // 分岐や符号をまだ証明できない場合だけ精度を上げて再試行する。
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::int64_t checkedSignedShift(std::size_t magnitude) {
    if (magnitude > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        error::throwCalcError(
            error::CalcErrorType::Overflow,
            "nextpow2 exponent is too large");
    return static_cast<std::int64_t>(magnitude);
}

[[nodiscard]] std::int64_t nextPow2Exponent(const Rational& value) {
    if (value <= Rational{BigInt{0}})
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "nextpow2 requires a positive real argument");

    const BigInt& numerator = value.numerator();
    const BigInt& denominator = value.denominator();
    const std::size_t numeratorBits = numerator.bitLength();
    const std::size_t denominatorBits = denominator.bitLength();

    std::int64_t floorExponent = 0;
    bool exactPower = false;
    if (numeratorBits >= denominatorBits) {
        const std::size_t shift = numeratorBits - denominatorBits;
        const std::int64_t signedShift = checkedSignedShift(shift);
        const BigInt scaledDenominator = denominator << shift;
        if (numerator < scaledDenominator)
            floorExponent = signedShift - 1;
        else {
            floorExponent = signedShift;
            exactPower = numerator == scaledDenominator;
        }
    }
    else {
        const std::size_t shift = denominatorBits - numeratorBits;
        const std::int64_t signedShift = checkedSignedShift(shift);
        const BigInt scaledNumerator = numerator << shift;
        if (scaledNumerator < denominator)
            floorExponent = -signedShift - 1;
        else {
            floorExponent = -signedShift;
            exactPower = scaledNumerator == denominator;
        }
    }

    if (exactPower)
        return floorExponent;
    if (floorExponent == std::numeric_limits<std::int64_t>::max())
        error::throwCalcError(
            error::CalcErrorType::Overflow,
            "nextpow2 exponent is too large");
    return floorExponent + 1;
}

[[nodiscard]] std::optional<std::int64_t> certifyNextPow2(
    const Expr& expression,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const approximation::CertifiedEvaluator evaluator{registry, mathematics, angles};
    for (std::size_t bits = 64; bits <= 4096; bits *= 2) {
        try {
            const auto enclosed = evaluator.enclose(
                expression, bits,
                approximation::CertifiedEvaluator::EnclosureKind::Information);
            if (!enclosed)
                return std::nullopt;
            if (!enclosed->isReal())
                error::throwCalcError(
                    error::CalcErrorType::Type,
                    "nextpow2 requires a real argument");

            const Rational lower = enclosed->asReal().lower().toRational();
            const Rational upper = enclosed->asReal().upper().toRational();
            if (upper <= Rational{BigInt{0}})
                error::throwCalcError(
                    error::CalcErrorType::Domain,
                    "nextpow2 requires a positive real argument");
            if (lower <= Rational{BigInt{0}})
                continue;

            const std::int64_t lowExponent = nextPow2Exponent(lower);
            const std::int64_t highExponent = nextPow2Exponent(upper);
            if (lowExponent == highExponent)
                return lowExponent;
        }
        catch (const approximation::PrecisionInsufficient&) {
            // power-of-two境界や符号が未確定なら精度を増やして再試行する。
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

[[nodiscard]] BigInt requireInteger(const Expr& expression, std::string_view name) {
    const BigInt* value = exactInteger(expression);
    if (!value)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " requires integer arguments");
    return *value;
}


[[nodiscard]] Expr factorList(std::vector<std::uint64_t> factors, bool negative) {
    std::sort(factors.begin(), factors.end());
    std::vector<Expr> elements;
    elements.reserve(factors.size() + (negative ? 1U : 0U));
    if (negative)
        elements.emplace_back(Number{BigInt{-1}});
    for (const std::uint64_t factor : factors)
        elements.emplace_back(Number{BigInt::fromUnsigned(factor)});
    return Expr::list(std::move(elements));
}

} // namespace

Expr evaluateFloor(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::floor);
    const auto value = expression::exact::realRational(arguments.front());
    if (!value) {
        if (arguments.front().isNumber())
            error::throwCalcError(error::CalcErrorType::Type, "floor requires a real argument");
        if (const auto certified = certifyIntegralRounding(
                arguments.front(), IntegralRoundingOperation::Floor,
                registry, mathematics, angles))
            return integerResult(*certified);
        return holdUnary(arguments, registry, evaluation::BuiltinId::Floor, names::floor);
    }
    return integerResult(floorRational(*value));
}

Expr evaluateCeil(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::ceil);
    const auto value = expression::exact::realRational(arguments.front());
    if (!value) {
        if (arguments.front().isNumber())
            error::throwCalcError(error::CalcErrorType::Type, "ceil requires a real argument");
        if (const auto certified = certifyIntegralRounding(
                arguments.front(), IntegralRoundingOperation::Ceil,
                registry, mathematics, angles))
            return integerResult(*certified);
        return holdUnary(arguments, registry, evaluation::BuiltinId::Ceil, names::ceil);
    }
    return integerResult(ceilRational(*value));
}

Expr evaluateTrunc(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::trunc);
    const auto value = expression::exact::realRational(arguments.front());
    if (!value) {
        if (arguments.front().isNumber())
            error::throwCalcError(error::CalcErrorType::Type, "trunc requires a real argument");
        if (const auto certified = certifyIntegralRounding(
                arguments.front(), IntegralRoundingOperation::Trunc,
                registry, mathematics, angles))
            return integerResult(*certified);
        return holdUnary(arguments, registry, evaluation::BuiltinId::Trunc, names::trunc);
    }
    return integerResult(truncRational(*value));
}

Expr evaluateRound(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() < 1 || arguments.size() > 2)
        error::throwCalcError(error::CalcErrorType::Type, "round expects 1 or 2 arguments");
    if (arguments.size() == 1) {
        const auto value = expression::exact::realRational(arguments.front());
        if (!value) {
            if (arguments.front().isNumber())
                error::throwCalcError(error::CalcErrorType::Type, "round requires a real argument");
            if (const auto certified = certifyIntegralRounding(
                    arguments.front(), IntegralRoundingOperation::Round,
                    registry, mathematics, angles))
                return integerResult(*certified);
            return holdUnary(arguments, registry, evaluation::BuiltinId::Round, names::round);
        }
        return integerResult(roundNearestEven(*value));
    }

    const std::int64_t digits = requireSignedSmallInteger(arguments[1], names::round);
    if (const auto value = expression::exact::realRational(arguments[0]))
        return Expr{Number{roundDecimal(*value, digits)}};
    if (arguments[0].isNumber())
        error::throwCalcError(error::CalcErrorType::Type, "round requires a real argument");
    if (const auto certified = certifyDecimalRound(
            arguments[0], digits, registry, mathematics, angles))
        return Expr{Number{*certified}};
    return Expr::call(registry.symbol(evaluation::BuiltinId::Round), {arguments[0], arguments[1]});
}

Expr evaluateFrac(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::frac);
    const auto value = expression::exact::realRational(arguments.front());
    if (!value) {
        if (arguments.front().isNumber())
            error::throwCalcError(error::CalcErrorType::Type, "frac requires a real argument");
        if (const auto certified = certifyIntegralRounding(
                arguments.front(), IntegralRoundingOperation::Floor,
                registry, mathematics, angles))
            return Expr::call(
                registry.symbol(evaluation::BuiltinId::Subtract),
                {arguments.front(), integerResult(*certified)});
        return holdUnary(arguments, registry, evaluation::BuiltinId::Frac, names::frac);
    }
    return Expr{Number{*value - Rational{floorRational(*value)}}};
}

Expr evaluateBitAnd(std::span<const Expr> arguments, const evaluation::BuiltinRegistry&) {
    BigInt result = requireInteger(arguments.front(), names::bitAnd);
    for (std::size_t i = 1; i < arguments.size(); ++i)
        result &= requireInteger(arguments[i], names::bitAnd);
    return integerResult(std::move(result));
}

Expr evaluateBitOr(std::span<const Expr> arguments, const evaluation::BuiltinRegistry&) {
    BigInt result = requireInteger(arguments.front(), names::bitOr);
    for (std::size_t i = 1; i < arguments.size(); ++i)
        result |= requireInteger(arguments[i], names::bitOr);
    return integerResult(std::move(result));
}

Expr evaluateBitXor(std::span<const Expr> arguments, const evaluation::BuiltinRegistry&) {
    BigInt result = requireInteger(arguments.front(), names::bitXor);
    for (std::size_t i = 1; i < arguments.size(); ++i)
        result ^= requireInteger(arguments[i], names::bitXor);
    return integerResult(std::move(result));
}

Expr evaluateBitNot(std::span<const Expr> arguments, const evaluation::BuiltinRegistry&) {
    requireArity(arguments, 1, names::bitNot);
    return integerResult(~requireInteger(arguments[0], names::bitNot));
}

Expr evaluateBitShiftLeft(std::span<const Expr> arguments, const evaluation::BuiltinRegistry&) {
    requireArity(arguments, 2, names::bitShiftLeft);
    BigInt value = requireInteger(arguments[0], names::bitShiftLeft);
    const std::int64_t shift = requireSignedSmallInteger(arguments[1], names::bitShiftLeft, 10'000'000);
    if (shift >= 0)
        value <<= static_cast<std::size_t>(shift);
    else
        value = arithmeticShiftRight(std::move(value), static_cast<std::size_t>(-shift));
    return integerResult(std::move(value));
}

Expr evaluateBitShiftRight(std::span<const Expr> arguments, const evaluation::BuiltinRegistry&) {
    requireArity(arguments, 2, names::bitShiftRight);
    BigInt value = requireInteger(arguments[0], names::bitShiftRight);
    const std::int64_t shift = requireSignedSmallInteger(arguments[1], names::bitShiftRight, 10'000'000);
    if (shift >= 0)
        value = arithmeticShiftRight(std::move(value), static_cast<std::size_t>(shift));
    else
        value <<= static_cast<std::size_t>(-shift);
    return integerResult(std::move(value));
}

Expr evaluateBitLength(std::span<const Expr> arguments, const evaluation::BuiltinRegistry&) {
    requireArity(arguments, 1, names::bitLength);
    const BigInt value = requireInteger(arguments[0], names::bitLength);
    return integerResult(BigInt::fromUnsigned(static_cast<std::uint64_t>(value.bitLength())));
}

Expr evaluateBitCount(std::span<const Expr> arguments, const evaluation::BuiltinRegistry&) {
    requireArity(arguments, 1, names::bitCount);
    const BigInt value = requireInteger(arguments[0], names::bitCount);
    if (value.isNegative())
        error::throwCalcError(error::CalcErrorType::Domain,
            "bitcount is defined only for nonnegative integers");
    return integerResult(BigInt::fromUnsigned(static_cast<std::uint64_t>(value.populationCount())));
}

Expr evaluateBitGet(std::span<const Expr> arguments, const evaluation::BuiltinRegistry&) {
    requireArity(arguments, 2, names::bitGet);
    const BigInt value = requireInteger(arguments[0], names::bitGet);
    const std::size_t index = requireBitIndex(arguments[1], names::bitGet);
    return integerResult(BigInt{value.testBit(index) ? 1 : 0});
}


Expr evaluateGcd(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    static_cast<void>(registry);
    BigInt result;
    bool first = true;
    for (const Expr& argument : arguments) {
        const BigInt value = requireInteger(argument, names::gcd);
        result = first ? value.abs() : numeric::gcd(std::move(result), value);
        first = false;
    }
    return integerResult(std::move(result));
}

Expr evaluateLcm(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    static_cast<void>(registry);
    BigInt result{1};
    for (const Expr& argument : arguments)
        result = numeric::lcm(result, requireInteger(argument, names::lcm));
    return integerResult(std::move(result));
}

Expr evaluateQuotient(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    static_cast<void>(registry);
    requireArity(arguments, 2, names::quotient);
    const BigInt lhs = requireInteger(arguments[0], names::quotient);
    const BigInt rhs = requireInteger(arguments[1], names::quotient);
    if (rhs.isZero())
        error::throwCalcError(error::CalcErrorType::Domain, "quotient divisor must be nonzero");
    return integerResult(lhs / rhs);
}

Expr evaluateIsPrime(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::isPrime);
    const BigInt* integerValue = exactInteger(arguments.front());
    if (!integerValue)
        error::throwCalcError(error::CalcErrorType::Type, "isprime requires an integer argument");
    if (integerValue->isNegative())
        return Expr{false};
    const auto value = numeric::tryToUint64(*integerValue);
    if (!value)
        return holdUnary(arguments, registry, evaluation::BuiltinId::IsPrime, names::isPrime);
    return Expr{numeric::isPrimeUint64(*value)};
}

Expr evaluateNextPrime(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::nextPrime);
    const BigInt* integerValue = exactInteger(arguments.front());
    if (!integerValue)
        error::throwCalcError(error::CalcErrorType::Type, "nextprime requires an integer argument");
    if (integerValue->isNegative() || *integerValue < BigInt{2})
        return integerResult(BigInt{2});
    const auto value = numeric::tryToUint64(*integerValue);
    if (!value || *value >= std::numeric_limits<std::uint64_t>::max() - 2)
        return holdUnary(arguments, registry, evaluation::BuiltinId::NextPrime, names::nextPrime);

    std::uint64_t candidate = *value + 1;
    if (candidate > 2 && (candidate & 1U) == 0)
        ++candidate;
    while (candidate > *value) {
        if (numeric::isPrimeUint64(candidate))
            return integerResult(BigInt::fromUnsigned(candidate));
        if (candidate > std::numeric_limits<std::uint64_t>::max() - 2)
            break;
        candidate += candidate == 2 ? 1 : 2;
    }
    return holdUnary(arguments, registry, evaluation::BuiltinId::NextPrime, names::nextPrime);
}

Expr evaluatePreviousPrime(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::previousPrime);
    const BigInt* integerValue = exactInteger(arguments.front());
    if (!integerValue)
        error::throwCalcError(error::CalcErrorType::Type, "prevprime requires an integer argument");
    if (*integerValue <= BigInt{2})
        error::throwCalcError(error::CalcErrorType::Domain, "prevprime has no positive prime below 2");
    const auto value = numeric::tryToUint64(*integerValue);
    if (!value)
        return holdUnary(arguments, registry, evaluation::BuiltinId::PreviousPrime, names::previousPrime);

    std::uint64_t candidate = *value - 1;
    if (candidate > 2 && (candidate & 1U) == 0)
        --candidate;
    for (;;) {
        if (numeric::isPrimeUint64(candidate))
            return integerResult(BigInt::fromUnsigned(candidate));
        if (candidate <= 3)
            return integerResult(BigInt{2});
        candidate -= 2;
    }
}

Expr evaluateFactorInteger(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::factorInteger);
    const BigInt* integerValue = exactInteger(arguments.front());
    if (!integerValue)
        error::throwCalcError(error::CalcErrorType::Type, "factorint requires an integer argument");
    if (integerValue->isZero())
        error::throwCalcError(error::CalcErrorType::Domain, "factorint is undefined for zero");

    const bool negative = integerValue->isNegative();
    const BigInt magnitude = integerValue->abs();
    const auto value = numeric::tryToUint64(magnitude);
    if (!value)
        return holdUnary(arguments, registry, evaluation::BuiltinId::FactorInteger, names::factorInteger);
    if (*value == 1)
        return Expr::list({Expr{Number{negative ? BigInt{-1} : BigInt{1}}}});

    std::vector<std::uint64_t> factors;
    if (!numeric::factorUint64(*value, factors))
        return holdUnary(arguments, registry, evaluation::BuiltinId::FactorInteger, names::factorInteger);
    return factorList(std::move(factors), negative);
}

Expr evaluateTotient(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::totient);
    const BigInt* integerValue = exactInteger(arguments.front());
    if (!integerValue)
        error::throwCalcError(error::CalcErrorType::Type, "totient requires an integer argument");
    if (!integerValue->isPositive())
        error::throwCalcError(error::CalcErrorType::Domain, "totient requires a positive integer");
    const auto value = numeric::tryToUint64(*integerValue);
    if (!value)
        return holdUnary(arguments, registry, evaluation::BuiltinId::Totient, names::totient);
    if (*value == 1)
        return integerResult(BigInt{1});

    std::vector<std::uint64_t> factors;
    if (!numeric::factorUint64(*value, factors))
        return holdUnary(arguments, registry, evaluation::BuiltinId::Totient, names::totient);
    std::sort(factors.begin(), factors.end());
    std::uint64_t result = *value;
    std::uint64_t previous = 0;
    for (const std::uint64_t prime : factors) {
        if (prime == previous)
            continue;
        result = (result / prime) * (prime - 1);
        previous = prime;
    }
    return integerResult(BigInt::fromUnsigned(result));
}

Expr evaluateRem(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    static_cast<void>(registry);
    requireArity(arguments, 2, names::rem);
    const BigInt lhs = requireInteger(arguments[0], names::rem);
    const BigInt rhs = requireInteger(arguments[1], names::rem);
    if (rhs.isZero())
        error::throwCalcError(error::CalcErrorType::Domain, "rem divisor must be nonzero");
    return integerResult(lhs % rhs);
}

Expr evaluateMod(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    static_cast<void>(registry);
    requireArity(arguments, 2, names::mod);
    const BigInt lhs = requireInteger(arguments[0], names::mod);
    const BigInt rhs = requireInteger(arguments[1], names::mod);
    if (rhs.isZero())
        error::throwCalcError(error::CalcErrorType::Domain, "mod modulus must be nonzero");

    // mod[a,m] = a - m floor(a/m)。結果の符号はmに従う。
    auto division = numeric::divmod(lhs, rhs);
    if (!division.remainder.isZero()
        && (lhs.isNegative() != rhs.isNegative()))
        division.quotient -= BigInt{1};
    return integerResult(lhs - rhs * division.quotient);
}

Expr evaluateNextPow2(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::nextPow2);
    if (const auto exact = expression::exact::realRational(arguments.front()))
        return integerResult(BigInt{nextPow2Exponent(*exact)});
    if (arguments.front().isNumber())
        error::throwCalcError(
            error::CalcErrorType::Type,
            "nextpow2 requires a real argument");
    if (const auto certified = certifyNextPow2(
            arguments.front(), registry, mathematics, angles))
        return integerResult(BigInt{*certified});
    return holdUnary(arguments, registry, evaluation::BuiltinId::NextPow2, names::nextPow2);
}

} // namespace mmcal::builtins
