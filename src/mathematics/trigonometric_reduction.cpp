// 角度還元
#include "trigonometric_reduction.hpp"

#include "numeric/big_int.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::mathematics {
namespace {

using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

} // namespace

Rational normalizeTurns(Rational turns) {
    // Rational は denominator > 0 に正規化されている。
    // したがって分子を分母で剰余化し、負なら1周期分だけ戻せば、浮動小数を使わず厳密に [0, 1) turn の代表値を得られる。
    const BigInt denominator = turns.denominator();
    BigInt remainder = turns.numerator() % denominator;
    if (remainder.isNegative())
        remainder += denominator;
    return Rational{std::move(remainder), denominator};
}

ReducedTrigAngle reduceTrigTurns(FunctionId function, Rational turns) {
    turns = normalizeTurns(std::move(turns));

    const Rational quarter = rational(1, 4);
    const Rational half = rational(1, 2);
    const Rational threeQuarters = rational(3, 4);
    const Rational one = rational(1);

    ReducedTrigAngle result{rational(0), false};

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
    case FunctionId::Erf:
    case FunctionId::Erfc:
    case FunctionId::FresnelC:
    case FunctionId::FresnelS:
    case FunctionId::Hypergeometric1F1:
    case FunctionId::Hypergeometric2F1:
    case FunctionId::EllipticF:
    case FunctionId::EllipticE:
    case FunctionId::EllipticPi:
    case FunctionId::Beta:
    case FunctionId::BetaLog:
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
        throw std::invalid_argument("Trigonometric reduction requires sin, cos, or tan");
    case FunctionId::Sin:
        // sin は各象限を次のように第1象限へ写せる。
        //   Q1:  sin(t)
        //   Q2:  sin(1/2 - t)
        //   Q3: -sin(t - 1/2)
        //   Q4: -sin(1 - t)
        if (turns <= quarter)
            result.referenceTurns = turns;
        else if (turns <= half)
            result.referenceTurns = half - turns;
        else if (turns <= threeQuarters) {
            result.referenceTurns = turns - half;
            result.negative = true;
        }
        else {
            result.referenceTurns = one - turns;
            result.negative = true;
        }
        return result;

    case FunctionId::Cos:
        // cos は
        //   Q1:  cos(t)
        //   Q2: -cos(1/2 - t)
        //   Q3: -cos(t - 1/2)
        //   Q4:  cos(1 - t)
        // として同じ第1象限の基準角へ写せる。
        if (turns <= quarter)
            result.referenceTurns = turns;
        else if (turns <= half) {
            result.referenceTurns = half - turns;
            result.negative = true;
        }
        else if (turns <= threeQuarters) {
            result.referenceTurns = turns - half;
            result.negative = true;
        }
        else
            result.referenceTurns = one - turns;
        return result;

    case FunctionId::Tan: {
        // tan の周期は1/2 turn。まず [0,1/2) に正規化し、
        //   [0,1/4]   :  tan(t)
        //   (1/4,1/2): -tan(1/2-t)
        // として第1象限へ写す。1/4 turn はpoleなので、呼出側が厳密に判定する。
        Rational halfTurn = normalizeTurns(turns * rational(2)) / rational(2);
        if (halfTurn <= quarter)
            result.referenceTurns = std::move(halfTurn);
        else {
            result.referenceTurns = half - halfTurn;
            result.negative = true;
        }
        return result;
    }
    }

    throw std::invalid_argument("Unsupported trigonometric function ID");
}

} // namespace mmcal::mathematics
