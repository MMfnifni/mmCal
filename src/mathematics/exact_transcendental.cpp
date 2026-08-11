// 指数・対数・特殊値
#include "exact_transcendental.hpp"

#include "error/error_message.hpp"
#include "exact_trigonometry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "value_facts.hpp"

#include <cstdint>
#include <optional>
#include <stdexcept>
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

[[nodiscard]] Expr numberExpr(Number value) {
    return Expr{std::move(value)};
}

[[nodiscard]] Expr integerExpr(std::int64_t value) {
    return numberExpr(Number{BigInt{value}});
}

[[nodiscard]] const expression::Symbol& constantSymbol(
    const MathRegistry& mathematics,
    ConstantId id) {
    const ConstantDefinition* definition = mathematics.findConstant(id);
    if (!definition)
        throw std::logic_error("Required mathematical constant is not registered");
    return definition->symbol;
}

[[nodiscard]] Expr piExpr(const MathRegistry& mathematics) {
    return Expr{constantSymbol(mathematics, ConstantId::Pi)};
}

[[nodiscard]] Expr eExpr(const MathRegistry& mathematics) {
    return Expr{constantSymbol(mathematics, ConstantId::E)};
}

[[nodiscard]] Expr imaginaryUnitExpr() {
    return numberExpr(Number::complex(RealNumber{BigInt{0}}, RealNumber{BigInt{1}}));
}

[[nodiscard]] bool isConstant(
    const Expr& expression,
    const MathRegistry& mathematics,
    ConstantId id) {
    if (!expression.isSymbol())
        return false;
    const ConstantDefinition* definition = mathematics.findConstant(expression.asSymbol());
    return definition && definition->id == id;
}

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return expression.isCall()
        && expression.asCall().head.sameIdentity(builtins.symbol(id));
}

[[nodiscard]] bool isExactZero(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] bool isExactOne(const Expr& expression) {
    return expression.isNumber()
        && expression.asNumber().isReal()
        && expression.asNumber().asReal() == RealNumber{BigInt{1}};
}


[[nodiscard]] bool isExactValidLogBase(const Expr& expression) {
    if (!expression.isNumber())
        return false;
    const Number& number = expression.asNumber();
    return !number.isZero() && !(number == Number{BigInt{1}});
}

struct IntegerLikeRational final {
    BigInt magnitude;
    bool reciprocal = false;
};

// n または 1/n (nは2以上の正整数) だけを抽出する。
// 任意Rationalの素因数分解をここへ持ち込まず、安価に証明できるexact冪だけを扱う。
[[nodiscard]] std::optional<IntegerLikeRational> integerLikeRational(const Rational& value) {
    if (value.numerator().isNegative() || value.numerator().isZero())
        return std::nullopt;

    if (value.denominator() == BigInt{1} && value.numerator() > BigInt{1})
        return IntegerLikeRational{value.numerator(), false};
    if (value.numerator() == BigInt{1} && value.denominator() > BigInt{1})
        return IntegerLikeRational{value.denominator(), true};
    return std::nullopt;
}

[[nodiscard]] std::optional<BigInt> exactIntegerPowerExponent(
    const BigInt& base,
    BigInt value) {
    if (base <= BigInt{1} || value.isNegative() || value.isZero())
        return std::nullopt;
    if (value == BigInt{1})
        return BigInt{0};

    BigInt exponent{0};
    while (value > BigInt{1}) {
        auto division = numeric::divmod(value, base);
        if (!division.remainder.isZero())
            return std::nullopt;
        value = std::move(division.quotient);
        exponent += BigInt{1};
    }
    return exponent;
}

[[nodiscard]] std::optional<Rational> exactRationalLogExponent(
    const Rational& base,
    const Rational& value) {
    if (value == Rational{BigInt{1}})
        return Rational{BigInt{0}};

    const auto normalizedBase = integerLikeRational(base);
    const auto normalizedValue = integerLikeRational(value);
    if (!normalizedBase || !normalizedValue)
        return std::nullopt;

    const bool negative = normalizedBase->reciprocal != normalizedValue->reciprocal;
    if (auto exponent = exactIntegerPowerExponent(
        normalizedBase->magnitude, normalizedValue->magnitude)) {
        Rational result{std::move(*exponent)};
        return negative ? -result : result;
    }

    // 逆向きに base=value^k なら log_base(value)=1/k。
    // これで log[4,2]=1/2, log[8,2]=1/3 も因数分解なしでexact化できる。
    if (auto inverseExponent = exactIntegerPowerExponent(
        normalizedValue->magnitude, normalizedBase->magnitude)) {
        if (inverseExponent->isZero())
            return std::nullopt;
        Rational result{BigInt{1}, std::move(*inverseExponent)};
        return negative ? -result : result;
    }
    return std::nullopt;
}

[[nodiscard]] Expr negateExpr(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (expression.isNumber())
        return numberExpr(-expression.asNumber());
    return Expr::call(builtins.symbol(BuiltinId::Negate), {std::move(expression)});
}

[[nodiscard]] Expr multiplyExpr(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins) {
    if (isExactZero(lhs) || isExactZero(rhs))
        return integerExpr(0);
    if (isExactOne(lhs))
        return rhs;
    if (isExactOne(rhs))
        return lhs;
    if (lhs.isNumber() && rhs.isNumber())
        return numberExpr(lhs.asNumber() * rhs.asNumber());
    return Expr::call(
        builtins.symbol(BuiltinId::Multiply),
        {std::move(lhs), std::move(rhs)});
}

[[nodiscard]] Expr addExpr(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins) {
    if (isExactZero(lhs))
        return rhs;
    if (isExactZero(rhs))
        return lhs;
    if (lhs.isNumber() && rhs.isNumber())
        return numberExpr(lhs.asNumber() + rhs.asNumber());
    return Expr::call(
        builtins.symbol(BuiltinId::Add),
        {std::move(lhs), std::move(rhs)});
}

[[nodiscard]] Expr divideExpr(
    Expr numerator,
    Expr denominator,
    const evaluation::BuiltinRegistry& builtins) {
    if (isExactOne(denominator))
        return numerator;
    if (numerator.isNumber() && denominator.isNumber())
        return numberExpr(numerator.asNumber() / denominator.asNumber());
    return Expr::call(
        builtins.symbol(BuiltinId::Divide),
        {std::move(numerator), std::move(denominator)});
}

// q*Pi を読みやすい式へ戻す。内部的にはRational qを保ったままなので、
// Piの数値近似は一切発生しない。
[[nodiscard]] Expr piMultipleExpr(
    const Rational& coefficient,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (coefficient.isZero())
        return integerExpr(0);

    const BigInt& numerator = coefficient.numerator();
    const BigInt& denominator = coefficient.denominator();
    const bool negative = numerator.isNegative();
    const BigInt magnitude = numerator.abs();

    Expr numeratorExpr = piExpr(mathematics);
    if (!(magnitude == BigInt{1}))
        numeratorExpr = multiplyExpr(numberExpr(Number{magnitude}), std::move(numeratorExpr), builtins);

    if (negative)
        numeratorExpr = negateExpr(std::move(numeratorExpr), builtins);
    return denominator == BigInt{1}
        ? std::move(numeratorExpr)
        : divideExpr(std::move(numeratorExpr), numberExpr(Number{denominator}), builtins);
}

[[nodiscard]] std::optional<Rational> pureImaginaryCoefficient(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isComplex())
        return std::nullopt;

    const auto& complex = expression.asNumber().asComplex();
    if (!complex.real.isZero())
        return std::nullopt;
    return complex.imaginary.toRational();
}

// expression = I * q * Pi の形だけを厳密に認識する。
// ExpのEuler特殊値に使うが、一般の複素式を無理に極形式へ変換しない。
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

    Rational imaginaryScale = rational(1);
    bool foundImaginary = false;
    bool foundPi = false;

    for (const Expr& factor : expression.asCall().arguments) {
        if (const auto imaginary = pureImaginaryCoefficient(factor)) {
            if (foundImaginary)
                return std::nullopt;
            imaginaryScale *= *imaginary;
            foundImaginary = true;
            continue;
        }

        if (factor.isNumber() && factor.asNumber().isReal()) {
            imaginaryScale *= factor.asNumber().asReal().toRational();
            continue;
        }

        const auto pi = extractRationalPiMultiple(factor, builtins, mathematics);
        if (!pi || foundPi)
            return std::nullopt;
        imaginaryScale *= *pi;
        foundPi = true;
    }

    return foundImaginary && foundPi
        ? std::optional<Rational>{std::move(imaginaryScale)}
        : std::nullopt;
}

[[nodiscard]] Expr explicitRadianAngle(
    const Rational& piCoefficient,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    return Expr::call(
        builtins.symbol(BuiltinId::UnitApplied),
        {piMultipleExpr(piCoefficient, builtins, mathematics), Expr{std::string{"Rad"}}});
}

[[nodiscard]] std::optional<Expr> expImaginaryPi(
    const Rational& coefficient,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    const Expr angle = explicitRadianAngle(coefficient, builtins, mathematics);
    auto cosine = simplifyExactTrig(
        FunctionId::Cos, angle, builtins, mathematics, AngleSemantics{AngleUnit::Degree});
    auto sine = simplifyExactTrig(
        FunctionId::Sin, angle, builtins, mathematics, AngleSemantics{AngleUnit::Degree});
    if (!cosine || !sine)
        return std::nullopt;

    Expr imaginaryPart = multiplyExpr(imaginaryUnitExpr(), std::move(*sine), builtins);
    return addExpr(std::move(*cosine), std::move(imaginaryPart), builtins);
}

[[nodiscard]] std::optional<Rational> exactArgForNumber(const Number& number) {
    if (number.isZero())
        return std::nullopt;

    if (number.isReal())
        return number.asReal().isNegative() ? rational(1) : rational(0);

    const auto& complex = number.asComplex();
    const RealNumber& x = complex.real;
    const RealNumber& y = complex.imaginary;

    if (x.isZero())
        return y.isNegative() ? rational(-1, 2) : rational(1, 2);

    // x,yが同じ絶対値なら45度刻みの偏角をexactに決定できる。
    // 一般のatan2は数値/記号函数の責務なのでここでは推測しない。
    if (x.abs() == y.abs()) {
        if (!x.isNegative() && !y.isNegative()) return rational(1, 4);
        if (x.isNegative() && !y.isNegative()) return rational(3, 4);
        if (x.isNegative() && y.isNegative()) return rational(-3, 4);
        return rational(-1, 4);
    }

    return std::nullopt;
}

[[nodiscard]] Expr principalLogNegativeReal(
    Expr magnitude,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    Expr realPart = isExactOne(magnitude)
        ? integerExpr(0)
        : (isConstant(magnitude, mathematics, ConstantId::E)
            ? integerExpr(1)
            : Expr::call(builtins.symbol(BuiltinId::Log), {std::move(magnitude)}));
    Expr imaginaryPart = multiplyExpr(imaginaryUnitExpr(), piExpr(mathematics), builtins);
    return addExpr(std::move(realPart), std::move(imaginaryPart), builtins);
}

} // namespace

std::optional<Expr> simplifyExactArg(
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (argument.isNumber()) {
        if (argument.asNumber().isZero())
            error::throwCalcError(error::CalcErrorType::Domain, "Arg is undefined at zero");
        if (const auto coefficient = exactArgForNumber(argument.asNumber()))
            // Argはsession既定単位から独立したprincipal angleを返すため、出力には必ずRadを明示する。
            // これにより既定単位を変更しても sin[arg[I]] は sin[(Pi/2) Rad] として正しく1へ簡約できる。
            return explicitRadianAngle(*coefficient, builtins, mathematics);
        return std::nullopt;
    }

    const ValueFacts facts = inferValueFacts(argument, builtins, mathematics);
    if (facts.isProvablyReal()) {
        if (facts.sign == RealSign::Positive)
            return explicitRadianAngle(rational(0), builtins, mathematics);
        if (facts.sign == RealSign::Negative)
            return explicitRadianAngle(rational(1), builtins, mathematics);
        if (facts.sign == RealSign::Zero)
            error::throwCalcError(error::CalcErrorType::Domain, "Arg is undefined at zero");
    }

    return std::nullopt;
}

std::optional<Expr> simplifyExactLog(
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (argument.isNumber()) {
        const Number& number = argument.asNumber();
        if (number.isZero())
            error::throwCalcError(error::CalcErrorType::Domain, "Log is undefined at zero");

        if (number.isReal()) {
            const RealNumber& real = number.asReal();
            if (real == RealNumber{BigInt{1}})
                return integerExpr(0);
            if (real.isNegative())
                return principalLogNegativeReal(
                    numberExpr(Number{real.abs()}), builtins, mathematics);
            return std::nullopt;
        }

        // |I|=|-I|=1なので実部ln|z|は0。Argだけでprincipal Logがexactに決まる。
        const auto& complex = number.asComplex();
        if (complex.real.isZero()
            && complex.imaginary.abs() == RealNumber{BigInt{1}}) {
            const Rational arg = complex.imaginary.isNegative()
                ? rational(-1, 2)
                : rational(1, 2);
            return multiplyExpr(
                imaginaryUnitExpr(),
                piMultipleExpr(arg, builtins, mathematics),
                builtins);
        }
        return std::nullopt;
    }

    if (isConstant(argument, mathematics, ConstantId::E))
        return integerExpr(1);

    const ValueFacts facts = inferValueFacts(argument, builtins, mathematics);
    if (facts.isProvablyNegativeReal()) {
        Expr magnitude = isHead(argument, builtins, BuiltinId::Negate)
                && argument.asCall().arguments.size() == 1
            ? argument.asCall().arguments.front()
            : Expr::call(builtins.symbol(BuiltinId::Negate), {argument});
        return principalLogNegativeReal(std::move(magnitude), builtins, mathematics);
    }
    if (facts.isProvablyReal() && facts.sign == RealSign::Zero)
        error::throwCalcError(error::CalcErrorType::Domain, "Log is undefined at zero");

    return std::nullopt;
}


std::optional<Expr> simplifyExactLog(
    const Expr& base,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (base.isNumber()) {
        const Number& baseNumber = base.asNumber();
        if (baseNumber.isZero())
            error::throwCalcError(error::CalcErrorType::Domain, "Logarithm base cannot be zero");
        if (baseNumber == Number{BigInt{1}})
            error::throwCalcError(error::CalcErrorType::Domain, "Logarithm base cannot be one");
    }
    if (value.isNumber() && value.asNumber().isZero())
        error::throwCalcError(error::CalcErrorType::Domain, "Log is undefined at zero");

    // E底は自然対数そのもの。E!=0,1 は定数メタデータから既知。
    if (isConstant(base, mathematics, ConstantId::E))
        return Expr::call(builtins.symbol(BuiltinId::Log), {value});

    // exactなbaseについてのみ、Log[base,1]=0 / Log[base,base]=1 を安全に畳む。
    if (isExactValidLogBase(base)) {
        if (isExactOne(value))
            return integerExpr(0);
        if (base == value)
            return integerExpr(1);
    }

    // 正の整数/その逆数について、整数指数の冪関係を因数分解なしで検出する。
    if (base.isNumber() && value.isNumber()
        && base.asNumber().isReal() && value.asNumber().isReal()) {
        const RealNumber& baseReal = base.asNumber().asReal();
        const RealNumber& valueReal = value.asNumber().asReal();
        if (!baseReal.isNegative() && !baseReal.isZero()
            && !valueReal.isNegative() && !valueReal.isZero()) {
            if (auto exponent = exactRationalLogExponent(
                baseReal.toRational(), valueReal.toRational()))
                return numberExpr(Number{std::move(*exponent)});
        }
    }

    return std::nullopt;
}

std::optional<Expr> simplifyExactExp(
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (argument.isNumber()) {
        const Number& number = argument.asNumber();
        if (number.isReal()) {
            if (number.asReal().isZero())
                return integerExpr(1);
            if (number.asReal() == RealNumber{BigInt{1}})
                return eExpr(mathematics);
        }
    }

    // Exp[Log[z]] = z は principal Log の定義域 z != 0 では常に成立する。
    // z=0の可能性が残る記号式には適用しない。
    if (isHead(argument, builtins, BuiltinId::Log)) {
        const auto& logArguments = argument.asCall().arguments;
        if (logArguments.size() == 1) {
            const ValueFacts facts = inferValueFacts(logArguments[0], builtins, mathematics);
            const bool provablyNonZero = facts.sign == RealSign::Positive
                || facts.sign == RealSign::Negative
                || facts.sign == RealSign::NonZero
                || facts.provablyNonReal;
            if (provablyNonZero)
                return logArguments[0];
        }
    }

    if (const auto coefficient = extractImaginaryPiMultiple(argument, builtins, mathematics))
        return expImaginaryPi(*coefficient, builtins, mathematics);

    return std::nullopt;
}

} // namespace mmcal::mathematics
