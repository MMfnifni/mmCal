// 複素超越函数の保証付き評価
#include "certified_complex_transcendental.hpp"
#include "certified_precision.hpp"

#include "certified_atan.hpp"
#include "certification_error.hpp"
#include "certified_constants.hpp"
#include "certified_exponential.hpp"
#include "certified_logarithm.hpp"
#include "certified_sqrt.hpp"
#include "certified_trigonometry.hpp"
#include "interval_math.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <limits>
#include <stdexcept>

namespace mmcal::approximation {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}


[[nodiscard]] RealInterval exactReal(
    std::int64_t value,
    std::size_t precisionBits) {
    return RealInterval::fromRational(rational(value), precisionBits);
}

[[nodiscard]] RealInterval exactHalf(std::size_t precisionBits) {
    return RealInterval::fromRational(rational(1, 2), precisionBits);
}

[[nodiscard]] bool isExactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

[[nodiscard]] RealInterval piInterval(std::size_t precisionBits) {
    return enclosePi(precisionBits).interval;
}

[[nodiscard]] RealInterval halfPiInterval(std::size_t precisionBits) {
    return multiply(piInterval(precisionBits), exactHalf(precisionBits), precisionBits);
}

[[nodiscard]] RealInterval fullPrincipalArgumentRange(std::size_t precisionBits) {
    const RealInterval pi = piInterval(precisionBits);
    return RealInterval{-pi.upper(), pi.upper()};
}

[[nodiscard]] RealInterval zeroInterval(std::size_t precisionBits) {
    return exactReal(0, precisionBits);
}

} // namespace

CertifiedArgumentResult enclosePrincipalArgument(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified Arg precision must be at least one bit");

    const BigFloat zero;
    const RealInterval& x = value.real();
    const RealInterval& y = value.imaginary();
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 16, "Certified Arg working precision is too large");

    // z=0はArgの定義域外。区間が単に0を含むだけなら、真値が非零である可能性もあるため即DomainErrorにはしない。上位Nはprecisionを増やせる。
    if (isExactZero(x) && isExactZero(y))
        throw std::domain_error("Arg is undefined at zero");

    // 実軸上の点はatanを経由せずexactなbranch規約を使う。principal Argの負実軸は +Pi 側であり -Pi ではない。
    if (isExactZero(y)) {
        if (x.lower() > zero)
            return CertifiedArgumentResult{zeroInterval(precisionBits), 0};
        if (x.upper() < zero)
            return CertifiedArgumentResult{piInterval(precisionBits), 0};

        // xの符号が現precisionでは未確定。実軸上なので可能値は0または+Pi。
        const RealInterval pi = piInterval(precisionBits);
        return CertifiedArgumentResult{
            RealInterval{zeroInterval(precisionBits).lower(), pi.upper()},
            0
        };
    }

    // 上半平面 y>0 では、xの符号に関係なく
    //   Arg(x+iy) = Pi/2 - atan(x/y)
    // が成立する。分母yが正なのでinterval divisionも0を跨がない。
    if (y.lower() > zero) {
        const RealInterval ratio = divide(x, y, workBits);
        const CertifiedAtanResult atan = encloseAtan(ratio, workBits);
        const RealInterval result = subtract(
            halfPiInterval(workBits), atan.interval, workBits);
        return CertifiedArgumentResult{
            result.roundedOutward(precisionBits), atan.termsUsed};
    }

    // 下半平面 y<0 では
    //   Arg(x+iy) = -Pi/2 + atan(x/(-y))
    // とすれば同じく分母を正にできる。
    if (y.upper() < zero) {
        const RealInterval positiveY = negate(y);
        const RealInterval ratio = divide(x, positiveY, workBits);
        const CertifiedAtanResult atan = encloseAtan(ratio, workBits);
        const RealInterval result = add(
            negate(halfPiInterval(workBits)), atan.interval, workBits);
        return CertifiedArgumentResult{
            result.roundedOutward(precisionBits), atan.termsUsed};
    }

    // yが0を跨いでも右半平面 x>0 ならbranch cutには触れない。
    // この場合は通常の atan(y/x) が [-Pi/2,Pi/2] 内で連続に使える。
    if (x.lower() > zero) {
        const RealInterval ratio = divide(y, x, workBits);
        const CertifiedAtanResult atan = encloseAtan(ratio, workBits);
        return CertifiedArgumentResult{
            atan.interval.roundedOutward(precisionBits), atan.termsUsed};
    }

    // 左半平面でyが0を跨ぐ長方形はprincipal Argのbranch cut（負実軸）を跨ぐ。
    // 上側では +Pi に、下側では -Pi に近づくため、一つの通常区間で両方を包含するには [-Pi,+Pi] まで広げる必要がある。これは粗いが数学的に正しい。
    // 真値が実際にはcutの片側なら、入力区間が狭まった段階で上の分岐へ移る。
    return CertifiedArgumentResult{fullPrincipalArgumentRange(precisionBits), 0};
}

CertifiedComplexTranscendentalResult enclosePrincipalComplexLog(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified complex Log precision must be at least one bit");

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified complex Log working precision is too large");

    // |z| = sqrt(x^2+y^2)。同一変数の平方にはsquareIntervalを使い、[-a,a]*[-a,a] のようなdependencyによる偽の負値を作らない。
    const RealInterval xSquared = squareInterval(value.real(), workBits);
    const RealInterval ySquared = squareInterval(value.imaginary(), workBits);
    const RealInterval magnitudeSquared = add(xSquared, ySquared, workBits);

    if (magnitudeSquared.upper().isZero())
        throw std::domain_error("Log is undefined at zero");

    const RealInterval magnitude = encloseSqrt(magnitudeSquared, workBits).interval;
    if (magnitude.lower() <= BigFloat{}) {
        // 真値が非零でも、現在の区間幅が広すぎると|z|の正の下限を証明できない。
        // ここで偽のlog下限を作るより、上位の作業precision増加に委ねる。
        throw PrecisionInsufficient{
            "Certified complex Log cannot yet prove that the argument is nonzero"};
    }

    const CertifiedLogarithmResult realPart = encloseLogPositive(magnitude, workBits);
    const CertifiedArgumentResult imaginaryPart = enclosePrincipalArgument(value, workBits);

    return CertifiedComplexTranscendentalResult{
        ComplexInterval{
            realPart.interval.roundedOutward(precisionBits),
            imaginaryPart.interval.roundedOutward(precisionBits)},
        realPart.termsUsed + imaginaryPart.termsUsed
    };
}

CertifiedComplexTranscendentalResult encloseComplexExp(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified complex Exp precision must be at least one bit");

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified complex Exp working precision is too large");

    // 複素指数函数の定義
    //   Exp(x+yI) = exp(x) (cos(y) + I sin(y))
    // で評価する。ここでyは複素数の虚部という数学的なラジアン量であり、ユーザー向けAngleSemanticsとは独立に、複素指数内部の角度は常に数学上のRadianとして扱う。
    const CertifiedExponentialResult magnitude = encloseExp(value.real(), workBits);
    const CertifiedTrigEnclosure cosine = encloseCosRadianInterval(
        value.imaginary(), workBits);
    const CertifiedTrigEnclosure sine = encloseSinRadianInterval(
        value.imaginary(), workBits);

    const RealInterval realPart = multiply(
        magnitude.interval, cosine.interval, workBits);
    const RealInterval imaginaryPart = multiply(
        magnitude.interval, sine.interval, workBits);

    return CertifiedComplexTranscendentalResult{
        ComplexInterval{
            realPart.roundedOutward(precisionBits),
            imaginaryPart.roundedOutward(precisionBits)},
        magnitude.termsUsed + cosine.termsUsed + sine.termsUsed
    };
}

CertifiedComplexTranscendentalResult enclosePrincipalPower(
    const ComplexInterval& base,
    const ComplexInterval& exponent,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified Power precision must be at least one bit");

    // principal Power の定義を一箇所に固定する。
    //   Power(z,w) := Exp(w * Log(z))
    // ここでLogはprincipal Log、Argのbranchは (-Pi,Pi]。
    // したがって例えば (-8)^(1/3) は「実立方根 -2」ではなくprincipal値 1 + I*sqrt(3) を選ぶ。方程式 z^3=-8 の全3解はSolverのSolutionSetの責務。
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified Power working precision is too large");
    const CertifiedComplexTranscendentalResult logarithm =
        enclosePrincipalComplexLog(base, workBits);
    const ComplexInterval exponentTimesLog = multiply(
        exponent, logarithm.interval, workBits);
    CertifiedComplexTranscendentalResult result = encloseComplexExp(
        exponentTimesLog, workBits);
    result.interval = result.interval.roundedOutward(precisionBits);
    result.termsUsed += logarithm.termsUsed;
    return result;
}

} // namespace mmcal::approximation
