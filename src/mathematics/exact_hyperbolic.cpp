// 双曲線函数
#include "exact_hyperbolic.hpp"

#include "angle.hpp"
#include "error/error_message.hpp"
#include "exact_trigonometry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <utility>

namespace mmcal::mathematics {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

[[nodiscard]] Expr integerExpr(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] bool isExactReal(const Expr& expression, std::int64_t value) {
    return expression.isNumber()
        && expression.asNumber().isReal()
        && expression.asNumber().asReal() == RealNumber{BigInt{value}};
}

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return expression.isCall()
        && expression.asCall().head.sameIdentity(builtins.symbol(id));
}

[[nodiscard]] std::optional<Rational> pureImaginaryCoefficient(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isComplex())
        return std::nullopt;
    const auto& complex = expression.asNumber().asComplex();
    if (!complex.real.isZero())
        return std::nullopt;
    return complex.imaginary.toRational();
}

// z = I*q*Pi の形だけをexactに認識する。一般複素式を極形式へ推測しない。
[[nodiscard]] std::optional<Rational> extractImaginaryPiMultiple(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (isHead(expression, builtins, BuiltinId::Negate)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 1)
            return std::nullopt;
        auto coefficient = extractImaginaryPiMultiple(arguments[0], builtins, mathematics);
        return coefficient ? std::optional<Rational>{-*coefficient} : std::nullopt;
    }

    if (isHead(expression, builtins, BuiltinId::Divide)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 2 || !arguments[1].isNumber()
            || !arguments[1].asNumber().isReal())
            return std::nullopt;
        const Rational divisor = arguments[1].asNumber().asReal().toRational();
        if (divisor.isZero())
            return std::nullopt;
        auto numerator = extractImaginaryPiMultiple(arguments[0], builtins, mathematics);
        return numerator ? std::optional<Rational>{*numerator / divisor} : std::nullopt;
    }

    if (!isHead(expression, builtins, BuiltinId::Multiply))
        return std::nullopt;

    Rational scale = rational(1);
    bool foundImaginary = false;
    bool foundPi = false;
    for (const Expr& factor : expression.asCall().arguments) {
        if (const auto imaginary = pureImaginaryCoefficient(factor)) {
            if (foundImaginary)
                return std::nullopt;
            scale *= *imaginary;
            foundImaginary = true;
            continue;
        }
        if (factor.isNumber() && factor.asNumber().isReal()) {
            scale *= factor.asNumber().asReal().toRational();
            continue;
        }
        const auto pi = extractRationalPiMultiple(factor, builtins, mathematics);
        if (!pi || foundPi)
            return std::nullopt;
        scale *= *pi;
        foundPi = true;
    }
    return foundImaginary && foundPi ? std::optional<Rational>{scale} : std::nullopt;
}

[[nodiscard]] Expr piMultipleExpr(
    const Rational& coefficient,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    const auto* pi = mathematics.findConstant(ConstantId::Pi);
    if (!pi)
        return integerExpr(0);
    if (coefficient.isZero())
        return integerExpr(0);

    Expr numerator{pi->symbol};
    const bool negative = coefficient.numerator().isNegative();
    const BigInt magnitude = coefficient.numerator().abs();
    if (!(magnitude == BigInt{1}))
        numerator = Expr::call(
            builtins.symbol(BuiltinId::Multiply),
            {Expr{Number{magnitude}}, std::move(numerator)});
    if (negative)
        numerator = Expr::call(builtins.symbol(BuiltinId::Negate), {std::move(numerator)});
    if (coefficient.denominator() == BigInt{1})
        return numerator;
    return Expr::call(
        builtins.symbol(BuiltinId::Divide),
        {std::move(numerator), Expr{Number{coefficient.denominator()}}});
}

[[nodiscard]] Expr radianAngle(
    const Rational& piCoefficient,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    return Expr::call(
        builtins.symbol(BuiltinId::UnitApplied),
        {piMultipleExpr(piCoefficient, builtins, mathematics), Expr{std::string{"Rad"}}});
}

[[nodiscard]] Expr imaginaryUnitExpr() {
    return Expr{Number::complex(RealNumber{}, RealNumber{BigInt{1}})};
}

[[nodiscard]] Expr multiplyExpr(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins) {
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() * rhs.asNumber()};
    return Expr::call(builtins.symbol(BuiltinId::Multiply), {std::move(lhs), std::move(rhs)});
}

[[nodiscard]] Expr divideExpr(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins,
    const char* domainMessage) {
    if (rhs.isNumber() && rhs.asNumber().isZero())
        error::throwCalcError(error::CalcErrorType::Domain, domainMessage);
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() / rhs.asNumber()};
    return Expr::call(builtins.symbol(BuiltinId::Divide), {std::move(lhs), std::move(rhs)});
}

[[nodiscard]] std::optional<Expr> simplifyImaginaryPiHyperbolic(
    FunctionId function,
    const Rational& coefficient,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    const Expr angle = radianAngle(coefficient, builtins, mathematics);
    auto sine = simplifyExactTrig(
        FunctionId::Sin, angle, builtins, mathematics, AngleSemantics{AngleUnit::Degree});
    auto cosine = simplifyExactTrig(
        FunctionId::Cos, angle, builtins, mathematics, AngleSemantics{AngleUnit::Degree});
    if (!sine || !cosine)
        return std::nullopt;

    Expr hyperbolicSine = multiplyExpr(imaginaryUnitExpr(), *sine, builtins);
    Expr hyperbolicCosine = *cosine;

    switch (function) {
    case FunctionId::Sinh:
        return hyperbolicSine;
    case FunctionId::Cosh:
        return hyperbolicCosine;
    case FunctionId::Tanh:
        return divideExpr(
            std::move(hyperbolicSine), std::move(hyperbolicCosine), builtins,
            "tanh is undefined where cosh is zero");
    case FunctionId::Csch:
        return divideExpr(
            integerExpr(1), std::move(hyperbolicSine), builtins,
            "csch is undefined where sinh is zero");
    case FunctionId::Sech:
        return divideExpr(
            integerExpr(1), std::move(hyperbolicCosine), builtins,
            "sech is undefined where cosh is zero");
    case FunctionId::Coth:
        return divideExpr(
            std::move(hyperbolicCosine), std::move(hyperbolicSine), builtins,
            "coth is undefined where sinh is zero");
    default:
        return std::nullopt;
    }
}

} // namespace

std::optional<Expr> simplifyExactHyperbolic(
    FunctionId function,
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (isExactReal(argument, 0)) {
        switch (function) {
        case FunctionId::Sinh:
        case FunctionId::Tanh:
        case FunctionId::Asinh:
        case FunctionId::Atanh:
            return integerExpr(0);
        case FunctionId::Cosh:
        case FunctionId::Sech:
            return integerExpr(1);
        case FunctionId::Csch:
            error::throwCalcError(error::CalcErrorType::Domain, "csch is undefined at zero");
        case FunctionId::Coth:
            error::throwCalcError(error::CalcErrorType::Domain, "coth is undefined at zero");
        default:
            break;
        }
    }

    if (function == FunctionId::Acosh && isExactReal(argument, 1))
        return integerExpr(0);
    if (function == FunctionId::Atanh
        && (isExactReal(argument, 1) || isExactReal(argument, -1)))
        error::throwCalcError(error::CalcErrorType::Domain, "atanh is undefined at +/-1");

    if (const auto coefficient = extractImaginaryPiMultiple(argument, builtins, mathematics))
        return simplifyImaginaryPiHyperbolic(function, *coefficient, builtins, mathematics);

    return std::nullopt;
}

} // namespace mmcal::mathematics
