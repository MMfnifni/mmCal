// 三角函数の特殊値
#include "exact_trigonometry.hpp"
#include "expression/exact_value.hpp"

#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "trigonometric_reduction.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>

namespace mmcal::mathematics {
namespace {

using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

[[nodiscard]] Rational absRational(const Rational& value) {
    return value.numerator().isNegative() ? -value : value;
}

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    evaluation::BuiltinId id) {
    return expression.isCall()
        && expression.asCall().head.sameIdentity(builtins.symbol(id));
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

[[nodiscard]] std::optional<std::pair<Expr, AngleUnit>> explicitAngleUnit(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, evaluation::BuiltinId::UnitApplied))
        return std::nullopt;

    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() != 2 || !arguments[1].isString())
        return std::nullopt;

    const auto unit = AngleSemantics::parseUnit(arguments[1].asString());
    if (!unit)
        return std::nullopt;

    return std::pair<Expr, AngleUnit>{arguments[0], *unit};
}

[[nodiscard]] std::optional<Rational> extractExplicitAngleTurns(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    // UnitAppliedが角度式のどこかに明示されている場合、その単位をturnへ変換する。
    // turnを中間表現にすることで、Degree/Radian/Gradianを混ぜた後の加減乗除もPiの10進値へ落とさず厳密なRationalだけで扱える。
    if (const auto explicitUnit = explicitAngleUnit(expression, builtins)) {
        const Expr& value = explicitUnit->first;
        switch (explicitUnit->second) {
        case AngleUnit::Degree:
            if (const auto numeric = expression::exact::realRational(value))
                return *numeric / rational(360);
            return std::nullopt;

        case AngleUnit::Gradian:
            if (const auto numeric = expression::exact::realRational(value))
                return *numeric / rational(400);
            return std::nullopt;

        case AngleUnit::Radian:
            if (const auto numeric = expression::exact::realRational(value); numeric && numeric->isZero())
                return rational(0);
            if (const auto piMultiple = extractRationalPiMultiple(value, builtins, mathematics))
                return *piMultiple / rational(2);
            return std::nullopt;
        }
    }

    if (isHead(expression, builtins, evaluation::BuiltinId::Negate)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 1)
            return std::nullopt;
        auto turns = extractExplicitAngleTurns(arguments[0], builtins, mathematics);
        return turns ? std::optional<Rational>{-*turns} : std::nullopt;
    }

    if (isHead(expression, builtins, evaluation::BuiltinId::Add)
        || isHead(expression, builtins, evaluation::BuiltinId::Subtract)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 2)
            return std::nullopt;
        auto left = extractExplicitAngleTurns(arguments[0], builtins, mathematics);
        auto right = extractExplicitAngleTurns(arguments[1], builtins, mathematics);
        if (!left || !right)
            return std::nullopt;
        if (isHead(expression, builtins, evaluation::BuiltinId::Add))
            return *left + *right;
        return *left - *right;
    }

    if (isHead(expression, builtins, evaluation::BuiltinId::Multiply)) {
        Rational scale = rational(1);
        std::optional<Rational> angle;
        for (const Expr& factor : expression.asCall().arguments) {
            if (const auto numeric = expression::exact::realRational(factor)) {
                scale *= *numeric;
                continue;
            }

            if (angle)
                return std::nullopt;
            angle = extractExplicitAngleTurns(factor, builtins, mathematics);
            if (!angle)
                return std::nullopt;
        }
        return angle ? std::optional<Rational>{*angle * scale} : std::nullopt;
    }

    if (isHead(expression, builtins, evaluation::BuiltinId::Divide)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 2)
            return std::nullopt;
        auto numerator = extractExplicitAngleTurns(arguments[0], builtins, mathematics);
        const auto denominator = expression::exact::realRational(arguments[1]);
        if (!numerator || !denominator || denominator->isZero())
            return std::nullopt;
        return *numerator / *denominator;
    }

    return std::nullopt;
}

[[nodiscard]] Expr integerExpr(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr rationalExpr(std::int64_t numerator, std::int64_t denominator) {
    return Expr{Number{Rational{BigInt{numerator}, BigInt{denominator}}}};
}

[[nodiscard]] Expr sqrtHalfExpr(
    std::int64_t radicand,
    const evaluation::BuiltinRegistry& builtins) {
    Expr root = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sqrt),
        {integerExpr(radicand)});
    return Expr::call(
        builtins.symbol(evaluation::BuiltinId::Divide),
        {std::move(root), integerExpr(2)});
}

[[nodiscard]] Expr sqrtExpr(
    std::int64_t radicand,
    const evaluation::BuiltinRegistry& builtins) {
    return Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sqrt),
        {integerExpr(radicand)});
}

[[nodiscard]] Expr radicalPairQuarterExpr(
    bool subtract,
    const evaluation::BuiltinRegistry& builtins) {
    // 15度と75度は、45度±30度の加法定理からexactな根号式を得られる。
    // sin(15°) = sin(45°-30°)
    //          = (sqrt(6) - sqrt(2)) / 4
    // cos(15°) = cos(45°-30°)
    //          = (sqrt(6) + sqrt(2)) / 4
    // ここでも根号を小数近似へ落とさない。
    // sqrt[2], sqrt[6] はexactな記号式として残り、後段の代数簡約器が同類項をまとめられる。
    Expr left = sqrtExpr(6, builtins);
    Expr right = sqrtExpr(2, builtins);
    Expr numerator = Expr::call(
        builtins.symbol(subtract
            ? evaluation::BuiltinId::Subtract
            : evaluation::BuiltinId::Add),
        {std::move(left), std::move(right)});
    return Expr::call(
        builtins.symbol(evaluation::BuiltinId::Divide),
        {std::move(numerator), integerExpr(4)});
}

[[nodiscard]] Expr negateExpr(Expr value, const evaluation::BuiltinRegistry& builtins) {
    if (value.isNumber())
        return Expr{-value.asNumber()};
    return Expr::call(
        builtins.symbol(evaluation::BuiltinId::Negate),
        {std::move(value)});
}

[[nodiscard]] Expr rationalExpr(const Rational& value) {
    return Expr{Number{value}};
}

[[nodiscard]] Expr multiplyExpr(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins) {
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() * rhs.asNumber()};
    return Expr::call(
        builtins.symbol(evaluation::BuiltinId::Multiply),
        {std::move(lhs), std::move(rhs)});
}

[[nodiscard]] Expr divideExpr(
    Expr numerator,
    Expr denominator,
    const evaluation::BuiltinRegistry& builtins) {
    if (denominator.isNumber() && denominator.asNumber().isZero())
        error::throwCalcError(error::CalcErrorType::Domain, "Reciprocal trigonometric function is undefined at this angle");
    if (numerator.isNumber() && denominator.isNumber())
        return Expr{numerator.asNumber() / denominator.asNumber()};
    return Expr::call(
        builtins.symbol(evaluation::BuiltinId::Divide),
        {std::move(numerator), std::move(denominator)});
}

[[nodiscard]] Expr piExpr(const MathRegistry& mathematics) {
    const ConstantDefinition* pi = mathematics.findConstant(ConstantId::Pi);
    if (!pi)
        throw std::logic_error("Pi is not registered in MathRegistry");
    return Expr{pi->symbol};
}

[[nodiscard]] Expr piMultipleExpr(
    const Rational& coefficient,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (coefficient.isZero())
        return integerExpr(0);

    const BigInt numerator = coefficient.numerator();
    const BigInt denominator = coefficient.denominator();
    const bool negative = numerator.isNegative();
    const BigInt magnitude = numerator.abs();

    Expr value = piExpr(mathematics);
    if (!(magnitude == BigInt{1}))
        value = multiplyExpr(Expr{Number{magnitude}}, std::move(value), builtins);
    if (!(denominator == BigInt{1}))
        value = divideExpr(std::move(value), Expr{Number{denominator}}, builtins);
    return negative ? negateExpr(std::move(value), builtins) : std::move(value);
}

[[nodiscard]] Expr angleValueFromTurns(
    const Rational& turns,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AngleSemantics& angleSemantics) {
    switch (angleSemantics.defaultUnit()) {
    case AngleUnit::Degree:
        return rationalExpr(turns * rational(360));
    case AngleUnit::Gradian:
        return rationalExpr(turns * rational(400));
    case AngleUnit::Radian:
        return piMultipleExpr(turns * rational(2), builtins, mathematics);
    }
    return rationalExpr(turns * rational(360));
}

struct ScaledSqrt final {
    Rational scale{BigInt{1}};
    std::int64_t radicand = 0;
};

[[nodiscard]] std::optional<ScaledSqrt> extractScaledSqrt(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (isHead(expression, builtins, evaluation::BuiltinId::Negate)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 1)
            return std::nullopt;
        auto result = extractScaledSqrt(arguments[0], builtins);
        if (result)
            result->scale = -result->scale;
        return result;
    }

    if (isHead(expression, builtins, evaluation::BuiltinId::Divide)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 2)
            return std::nullopt;
        auto result = extractScaledSqrt(arguments[0], builtins);
        const auto denominator = expression::exact::realRational(arguments[1]);
        if (!result || !denominator || denominator->isZero())
            return std::nullopt;
        result->scale /= *denominator;
        return result;
    }

    if (isHead(expression, builtins, evaluation::BuiltinId::Multiply)) {
        Rational scale = rational(1);
        std::optional<ScaledSqrt> radical;
        for (const Expr& factor : expression.asCall().arguments) {
            if (const auto numeric = expression::exact::realRational(factor)) {
                scale *= *numeric;
                continue;
            }
            if (radical)
                return std::nullopt;
            radical = extractScaledSqrt(factor, builtins);
            if (!radical)
                return std::nullopt;
        }
        if (!radical)
            return std::nullopt;
        radical->scale *= scale;
        return radical;
    }

    if (!isHead(expression, builtins, evaluation::BuiltinId::Sqrt))
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() != 1 || !arguments[0].isNumber()
        || !arguments[0].asNumber().isReal()
        || !arguments[0].asNumber().asReal().isInteger())
        return std::nullopt;

    const BigInt& radicand = arguments[0].asNumber().asReal().asInteger();
    if (radicand == BigInt{2})
        return ScaledSqrt{rational(1), 2};
    if (radicand == BigInt{3})
        return ScaledSqrt{rational(1), 3};
    return std::nullopt;
}

[[nodiscard]] std::optional<Rational> inverseTrigTurns(
    FunctionId function,
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins) {
    if (const auto value = expression::exact::realRational(argument)) {
        if (function == FunctionId::Asin) {
            if (*value == rational(-1)) return rational(-1, 4);
            if (*value == rational(-1, 2)) return rational(-1, 12);
            if (value->isZero()) return rational(0);
            if (*value == rational(1, 2)) return rational(1, 12);
            if (*value == rational(1)) return rational(1, 4);
        }
        if (function == FunctionId::Acos) {
            if (*value == rational(-1)) return rational(1, 2);
            if (*value == rational(-1, 2)) return rational(1, 3);
            if (value->isZero()) return rational(1, 4);
            if (*value == rational(1, 2)) return rational(1, 6);
            if (*value == rational(1)) return rational(0);
        }
        if (function == FunctionId::Atan) {
            if (*value == rational(-1)) return rational(-1, 8);
            if (value->isZero()) return rational(0);
            if (*value == rational(1)) return rational(1, 8);
        }
    }

    const auto radical = extractScaledSqrt(argument, builtins);
    if (!radical)
        return std::nullopt;

    if (function == FunctionId::Asin) {
        if (radical->radicand == 2 && radical->scale == rational(1, 2)) return rational(1, 8);
        if (radical->radicand == 2 && radical->scale == rational(-1, 2)) return rational(-1, 8);
        if (radical->radicand == 3 && radical->scale == rational(1, 2)) return rational(1, 6);
        if (radical->radicand == 3 && radical->scale == rational(-1, 2)) return rational(-1, 6);
    }
    if (function == FunctionId::Acos) {
        if (radical->radicand == 2 && radical->scale == rational(1, 2)) return rational(1, 8);
        if (radical->radicand == 2 && radical->scale == rational(-1, 2)) return rational(3, 8);
        if (radical->radicand == 3 && radical->scale == rational(1, 2)) return rational(1, 12);
        if (radical->radicand == 3 && radical->scale == rational(-1, 2)) return rational(5, 12);
    }
    if (function == FunctionId::Atan && radical->radicand == 3) {
        if (radical->scale == rational(1, 3)) return rational(1, 12);
        if (radical->scale == rational(-1, 3)) return rational(-1, 12);
        if (radical->scale == rational(1)) return rational(1, 6);
        if (radical->scale == rational(-1)) return rational(-1, 6);
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> sineFirstQuadrant(
    const Rational& reference,
    const evaluation::BuiltinRegistry& builtins) {
    if (reference == rational(0))
        return integerExpr(0);
    if (reference == rational(1, 24))
        return radicalPairQuarterExpr(true, builtins);
    if (reference == rational(1, 12))
        return rationalExpr(1, 2);
    if (reference == rational(1, 8))
        return sqrtHalfExpr(2, builtins);
    if (reference == rational(1, 6))
        return sqrtHalfExpr(3, builtins);
    if (reference == rational(5, 24))
        return radicalPairQuarterExpr(false, builtins);
    if (reference == rational(1, 4))
        return integerExpr(1);
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> cosineFirstQuadrant(
    const Rational& reference,
    const evaluation::BuiltinRegistry& builtins) {
    if (reference == rational(0))
        return integerExpr(1);
    if (reference == rational(1, 24))
        return radicalPairQuarterExpr(false, builtins);
    if (reference == rational(1, 12))
        return sqrtHalfExpr(3, builtins);
    if (reference == rational(1, 8))
        return sqrtHalfExpr(2, builtins);
    if (reference == rational(1, 6))
        return rationalExpr(1, 2);
    if (reference == rational(5, 24))
        return radicalPairQuarterExpr(true, builtins);
    if (reference == rational(1, 4))
        return integerExpr(0);
    return std::nullopt;
}

[[nodiscard]] Expr addOrSubtractIntegerAndSqrt3(
    bool subtract,
    const evaluation::BuiltinRegistry& builtins) {
    return Expr::call(
        builtins.symbol(subtract
            ? evaluation::BuiltinId::Subtract
            : evaluation::BuiltinId::Add),
        {integerExpr(2), sqrtExpr(3, builtins)});
}

[[nodiscard]] std::optional<Expr> tangentFirstQuadrant(
    const Rational& reference,
    const evaluation::BuiltinRegistry& builtins) {
    if (reference == rational(0))
        return integerExpr(0);
    if (reference == rational(1, 24))
        return addOrSubtractIntegerAndSqrt3(true, builtins);
    if (reference == rational(1, 12))
        return Expr::call(
            builtins.symbol(evaluation::BuiltinId::Divide),
            {sqrtExpr(3, builtins), integerExpr(3)});
    if (reference == rational(1, 8))
        return integerExpr(1);
    if (reference == rational(1, 6))
        return sqrtExpr(3, builtins);
    if (reference == rational(5, 24))
        return addOrSubtractIntegerAndSqrt3(false, builtins);
    if (reference == rational(1, 4))
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "tan is undefined where cos is zero");
    return std::nullopt;
}


[[nodiscard]] Expr sqrt6PlusOrMinusSqrt2(
    bool subtract,
    const evaluation::BuiltinRegistry& builtins) {
    return Expr::call(
        builtins.symbol(subtract
            ? evaluation::BuiltinId::Subtract
            : evaluation::BuiltinId::Add),
        {sqrtExpr(6, builtins), sqrtExpr(2, builtins)});
}

[[nodiscard]] Expr twoSqrt3OverThree(
    const evaluation::BuiltinRegistry& builtins) {
    return Expr::call(
        builtins.symbol(evaluation::BuiltinId::Divide),
        {Expr::call(
            builtins.symbol(evaluation::BuiltinId::Multiply),
            {integerExpr(2), sqrtExpr(3, builtins)}),
         integerExpr(3)});
}

[[nodiscard]] std::optional<Expr> reciprocalTrigFirstQuadrant(
    FunctionId function,
    const Rational& reference,
    const evaluation::BuiltinRegistry& builtins) {
    switch (function) {
    case FunctionId::Cot:
        if (reference == rational(0))
            error::throwCalcError(error::CalcErrorType::Domain, "cot is undefined where sin is zero");
        if (reference == rational(1, 24)) return addOrSubtractIntegerAndSqrt3(false, builtins);
        if (reference == rational(1, 12)) return sqrtExpr(3, builtins);
        if (reference == rational(1, 8)) return integerExpr(1);
        if (reference == rational(1, 6))
            return Expr::call(
                builtins.symbol(evaluation::BuiltinId::Divide),
                {sqrtExpr(3, builtins), integerExpr(3)});
        if (reference == rational(5, 24)) return addOrSubtractIntegerAndSqrt3(true, builtins);
        if (reference == rational(1, 4)) return integerExpr(0);
        return std::nullopt;

    case FunctionId::Sec:
        if (reference == rational(0)) return integerExpr(1);
        if (reference == rational(1, 24)) return sqrt6PlusOrMinusSqrt2(true, builtins);
        if (reference == rational(1, 12)) return twoSqrt3OverThree(builtins);
        if (reference == rational(1, 8)) return sqrtExpr(2, builtins);
        if (reference == rational(1, 6)) return integerExpr(2);
        if (reference == rational(5, 24)) return sqrt6PlusOrMinusSqrt2(false, builtins);
        if (reference == rational(1, 4))
            error::throwCalcError(error::CalcErrorType::Domain, "sec is undefined where cos is zero");
        return std::nullopt;

    case FunctionId::Csc:
        if (reference == rational(0))
            error::throwCalcError(error::CalcErrorType::Domain, "csc is undefined where sin is zero");
        if (reference == rational(1, 24)) return sqrt6PlusOrMinusSqrt2(false, builtins);
        if (reference == rational(1, 12)) return integerExpr(2);
        if (reference == rational(1, 8)) return sqrtExpr(2, builtins);
        if (reference == rational(1, 6)) return twoSqrt3OverThree(builtins);
        if (reference == rational(5, 24)) return sqrt6PlusOrMinusSqrt2(true, builtins);
        if (reference == rational(1, 4)) return integerExpr(1);
        return std::nullopt;

    default:
        return std::nullopt;
    }
}

[[nodiscard]] std::optional<Expr> simplifyReducedReciprocalTrig(
    FunctionId function,
    Rational turns,
    const evaluation::BuiltinRegistry& builtins) {
    const FunctionId reductionFunction = function == FunctionId::Cot
        ? FunctionId::Tan
        : function == FunctionId::Sec ? FunctionId::Cos : FunctionId::Sin;
    const ReducedTrigAngle reduced = reduceTrigTurns(reductionFunction, std::move(turns));
    auto value = reciprocalTrigFirstQuadrant(function, reduced.referenceTurns, builtins);
    if (!value)
        return std::nullopt;
    return reduced.negative ? negateExpr(std::move(*value), builtins) : std::move(*value);
}

[[nodiscard]] std::optional<Expr> simplifyReducedTrig(
    FunctionId function,
    Rational turns,
    const evaluation::BuiltinRegistry& builtins) {
    const ReducedTrigAngle reduced = reduceTrigTurns(function, std::move(turns));

    std::optional<Expr> value;
    switch (function) {
    case FunctionId::Expm1:
    case FunctionId::Log1p:
    case FunctionId::Sinc:
    case FunctionId::Cosc:
    case FunctionId::Tanc:
    case FunctionId::Sinhc:
    case FunctionId::Tanhc:
    case FunctionId::Expc:
    case FunctionId::Gamma:
    case FunctionId::LogGamma:
    case FunctionId::LambertW:
    case FunctionId::Erf:
    case FunctionId::Erfc:
    case FunctionId::FresnelC:
    case FunctionId::FresnelS:
    case FunctionId::Hypergeometric1F1:
    case FunctionId::Hypergeometric2F1:
    case FunctionId::EllipticF:
    case FunctionId::EllipticE:
    case FunctionId::EllipticPi:
    case FunctionId::ExponentialIntegralEi:
    case FunctionId::SineIntegralSi:
    case FunctionId::CosineIntegralCi:
    case FunctionId::LogarithmicIntegralLi:
    case FunctionId::Polylog:
    case FunctionId::Beta:
    case FunctionId::BetaLog:
    case FunctionId::Zeta:
    case FunctionId::Digamma:
    case FunctionId::Trigamma:
    case FunctionId::IncompleteBeta:
    case FunctionId::Cbrt:
    case FunctionId::Hypot:
    case FunctionId::Cis:
    case FunctionId::Polar:
    case FunctionId::DegreeToRadian:
    case FunctionId::DegreeToGradian:
    case FunctionId::RadianToDegree:
    case FunctionId::RadianToGradian:
    case FunctionId::GradianToDegree:
    case FunctionId::GradianToRadian:
    case FunctionId::Sqrt:
    case FunctionId::Abs:
    case FunctionId::Sign:
    case FunctionId::Re:
    case FunctionId::Im:
    case FunctionId::Conj:
    case FunctionId::Cot:
    case FunctionId::Sec:
    case FunctionId::Csc:
    case FunctionId::Asin:
    case FunctionId::Acos:
    case FunctionId::Atan:
    case FunctionId::Atan2:
    case FunctionId::Sinh:
    case FunctionId::Cosh:
    case FunctionId::Tanh:
    case FunctionId::Asinh:
    case FunctionId::Acosh:
    case FunctionId::Atanh:
    case FunctionId::Csch:
    case FunctionId::Sech:
    case FunctionId::Coth:
    case FunctionId::Arg:
    case FunctionId::Log:
    case FunctionId::Exp:
    case FunctionId::Power:
        return std::nullopt;
    case FunctionId::Sin:
        value = sineFirstQuadrant(reduced.referenceTurns, builtins);
        break;
    case FunctionId::Cos:
        value = cosineFirstQuadrant(reduced.referenceTurns, builtins);
        break;
    case FunctionId::Tan:
        value = tangentFirstQuadrant(reduced.referenceTurns, builtins);
        break;
    }

    if (!value)
        return std::nullopt;
    return reduced.negative ? negateExpr(std::move(*value), builtins) : std::move(*value);
}

} // namespace

std::optional<Rational> extractRationalPiMultiple(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (isConstant(expression, mathematics, ConstantId::Pi))
        return rational(1);
    if (isHead(expression, builtins, evaluation::BuiltinId::Negate)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 1)
            return std::nullopt;
        auto coefficient = extractRationalPiMultiple(arguments[0], builtins, mathematics);
        return coefficient ? std::optional<Rational>{-*coefficient} : std::nullopt;
    }

    if (isHead(expression, builtins, evaluation::BuiltinId::Add)
        || isHead(expression, builtins, evaluation::BuiltinId::Subtract)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 2)
            return std::nullopt;
        const auto left = extractRationalPiMultiple(arguments[0], builtins, mathematics);
        const auto right = extractRationalPiMultiple(arguments[1], builtins, mathematics);
        if (!left || !right)
            return std::nullopt;
        return isHead(expression, builtins, evaluation::BuiltinId::Add)
            ? std::optional<Rational>{*left + *right}
            : std::optional<Rational>{*left - *right};
    }

    if (isHead(expression, builtins, evaluation::BuiltinId::Multiply)) {
        Rational coefficient = rational(1);
        bool foundPi = false;
        for (const Expr& factor : expression.asCall().arguments) {
            if (const auto numeric = expression::exact::realRational(factor)) {
                coefficient *= *numeric;
                continue;
            }

            if (foundPi)
                return std::nullopt;

            const auto piFactor = extractRationalPiMultiple(factor, builtins, mathematics);
            if (!piFactor)
                return std::nullopt;
            coefficient *= *piFactor;
            foundPi = true;
        }
        return foundPi ? std::optional<Rational>{std::move(coefficient)} : std::nullopt;
    }

    if (isHead(expression, builtins, evaluation::BuiltinId::Divide)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 2)
            return std::nullopt;
        auto numerator = extractRationalPiMultiple(arguments[0], builtins, mathematics);
        const auto denominator = expression::exact::realRational(arguments[1]);
        if (!numerator || !denominator || denominator->isZero())
            return std::nullopt;
        return *numerator / *denominator;
    }

    return std::nullopt;
}

std::optional<ExactAngle> extractExactAngle(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AngleSemantics& angleSemantics) {
    // 明示単位が式内にある場合は、まず単位付き角度算術としてturnへ変換する。
    // 例: (Pi Rad) / 6, 2 * (90 Deg), (100 Grad) + (90 Deg)。
    if (const auto explicitTurns = extractExplicitAngleTurns(expression, builtins, mathematics))
        return ExactAngle{*explicitTurns};

    // 明示単位が無い式はセッション既定単位として解釈する。
    switch (angleSemantics.defaultUnit()) {
    case AngleUnit::Degree:
        if (const auto numeric = expression::exact::realRational(expression))
            return ExactAngle{*numeric / rational(360)};
        return std::nullopt;

    case AngleUnit::Gradian:
        if (const auto numeric = expression::exact::realRational(expression))
            return ExactAngle{*numeric / rational(400)};
        return std::nullopt;

    case AngleUnit::Radian:
        if (const auto numeric = expression::exact::realRational(expression); numeric && numeric->isZero())
            return ExactAngle{rational(0)};
        if (const auto piMultiple = extractRationalPiMultiple(expression, builtins, mathematics))
            return ExactAngle{*piMultiple / rational(2)};
        return std::nullopt;
    }

    return std::nullopt;
}

std::optional<Expr> simplifyExactTrig(
    FunctionId function,
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AngleSemantics& angleSemantics) {
    const auto angle = extractExactAngle(argument, builtins, mathematics, angleSemantics);
    if (!angle)
        return std::nullopt;

    if (function == FunctionId::Cot || function == FunctionId::Sec || function == FunctionId::Csc)
        return simplifyReducedReciprocalTrig(function, angle->turns, builtins);

    const FunctionDefinition* definition = mathematics.findFunction(function);
    if (!definition || !definition->periodTurns)
        return std::nullopt;

    switch (function) {
    case FunctionId::Sin:
    case FunctionId::Cos:
    case FunctionId::Tan:
        return simplifyReducedTrig(function, angle->turns, builtins);
    default:
        return std::nullopt;
    }
}

std::optional<Expr> simplifyExactInverseTrig(
    FunctionId function,
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AngleSemantics& angleSemantics) {
    if (function != FunctionId::Asin && function != FunctionId::Acos
        && function != FunctionId::Atan)
        return std::nullopt;

    const auto turns = inverseTrigTurns(function, argument, builtins);
    if (!turns)
        return std::nullopt;
    return angleValueFromTurns(*turns, builtins, mathematics, angleSemantics);
}

std::optional<Expr> simplifyExactAtan2(
    const Expr& y,
    const Expr& x,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AngleSemantics& angleSemantics) {
    // atan2は実数座標専用。exactな複素数が渡された時点でdomain違反を確定できる。
    if ((y.isNumber() && !y.asNumber().isReal())
        || (x.isNumber() && !x.asNumber().isReal()))
        error::throwCalcError(error::CalcErrorType::Domain, "atan2 expects real arguments");

    const auto yValue = expression::exact::realRational(y);
    const auto xValue = expression::exact::realRational(x);
    if (!yValue || !xValue)
        return std::nullopt;

    if (xValue->isZero() && yValue->isZero())
        error::throwCalcError(error::CalcErrorType::Domain, "atan2 is undefined at (0, 0)");

    Rational turns;
    if (yValue->isZero())
        turns = xValue->numerator().isNegative() ? rational(1, 2) : rational(0);
    else if (xValue->isZero())
        turns = yValue->numerator().isNegative() ? rational(-1, 4) : rational(1, 4);
    else if (absRational(*yValue) == absRational(*xValue)) {
        if (!xValue->numerator().isNegative() && !yValue->numerator().isNegative()) turns = rational(1, 8);
        else if (xValue->numerator().isNegative() && !yValue->numerator().isNegative()) turns = rational(3, 8);
        else if (xValue->numerator().isNegative() && yValue->numerator().isNegative()) turns = rational(-3, 8);
        else turns = rational(-1, 8);
    }
    else
        return std::nullopt;

    return angleValueFromTurns(turns, builtins, mathematics, angleSemantics);
}

std::optional<RealNumber> extractExactRadianValue(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const AngleSemantics& angleSemantics) {
    Expr value = expression;
    AngleUnit unit = angleSemantics.defaultUnit();

    if (const auto explicitUnit = explicitAngleUnit(expression, builtins)) {
        value = explicitUnit->first;
        unit = explicitUnit->second;
    }

    // 一般のRationalラジアン値（例: 1 Rad）は、turnへexact変換するにはPiで割る必要がある。
    // そのためDegree/Gradや q*Pi Rad とは分け、Rationalラジアンのreference Taylor backendへ直接渡すためにここで抽出する。
    if (unit != AngleUnit::Radian)
        return std::nullopt;
    if (!value.isNumber() || !value.asNumber().isReal())
        return std::nullopt;
    return value.asNumber().asReal();
}

} // namespace mmcal::mathematics
