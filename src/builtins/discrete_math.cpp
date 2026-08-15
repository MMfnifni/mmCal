// 丸め・剰余・整数演算
#include "discrete_math.hpp"

#include "builtins/names.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/certification_error.hpp"
#include "error/error_message.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"

#include <compare>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

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

[[nodiscard]] std::optional<Rational> exactRealRational(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    return expression.asNumber().asReal().toRational();
}

[[nodiscard]] Expr integerResult(BigInt value) {
    return Expr{Number{std::move(value)}};
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

} // namespace

Expr evaluateFloor(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::floor);
    const auto value = exactRealRational(arguments.front());
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
    const auto value = exactRealRational(arguments.front());
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
    const auto value = exactRealRational(arguments.front());
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
    requireArity(arguments, 1, names::round);
    const auto value = exactRealRational(arguments.front());
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

Expr evaluateFrac(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::frac);
    const auto value = exactRealRational(arguments.front());
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
    if (const auto exact = exactRealRational(arguments.front()))
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
