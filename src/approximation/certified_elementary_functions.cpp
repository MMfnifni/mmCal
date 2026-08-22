// 安定初等函数の保証付き評価
#include "certified_elementary_functions.hpp"
#include "certified_precision.hpp"

#include "certification_error.hpp"
#include "certified_atan.hpp"
#include "certified_complex_sqrt.hpp"
#include "certified_complex_transcendental.hpp"
#include "certified_constants.hpp"
#include "certified_exponential.hpp"
#include "certified_logarithm.hpp"
#include "certified_sqrt.hpp"
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
    const Rational& value,
    std::size_t precisionBits) {
    return RealInterval::fromRational(value, precisionBits);
}

[[nodiscard]] RealInterval exactReal(
    std::int64_t value,
    std::size_t precisionBits) {
    return exactReal(rational(value), precisionBits);
}

[[nodiscard]] ComplexInterval exactComplex(
    std::int64_t real,
    std::int64_t imaginary,
    std::size_t precisionBits) {
    return ComplexInterval{exactReal(real, precisionBits), exactReal(imaginary, precisionBits)};
}

[[nodiscard]] ComplexInterval multiplyByI(const ComplexInterval& value) {
    return ComplexInterval{negate(value.imaginary()), value.real()};
}

[[nodiscard]] ComplexInterval multiplyByNegativeI(const ComplexInterval& value) {
    return ComplexInterval{value.imaginary(), negate(value.real())};
}

[[nodiscard]] ComplexInterval scaleComplex(
    const ComplexInterval& value,
    const Rational& scale,
    std::size_t precisionBits) {
    const RealInterval factor = exactReal(scale, precisionBits);
    return ComplexInterval{
        multiply(value.real(), factor, precisionBits),
        multiply(value.imaginary(), factor, precisionBits)};
}

[[nodiscard]] ComplexInterval oneComplex(std::size_t precisionBits) {
    return exactComplex(1, 0, precisionBits);
}

[[nodiscard]] RealInterval piInterval(std::size_t precisionBits) {
    return enclosePi(precisionBits).interval;
}

[[nodiscard]] RealInterval halfPiInterval(std::size_t precisionBits) {
    return multiply(piInterval(precisionBits), exactReal(rational(1, 2), precisionBits), precisionBits);
}

[[nodiscard]] RealInterval pointCubeRoot(
    const Rational& x,
    std::size_t precisionBits) {
    if (x.isZero())
        return exactReal(0, precisionBits);

    const bool negative = x.numerator().isNegative();
    const Rational magnitude = negative ? -x : x;
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified cube-root working precision is too large");
    const RealInterval logarithm = encloseLogPositive(
        exactReal(magnitude, workBits), workBits).interval;
    const RealInterval divided = multiply(
        logarithm, exactReal(rational(1, 3), workBits), workBits);
    RealInterval root = encloseExp(divided, workBits).interval.roundedOutward(precisionBits);
    return negative ? negate(root) : root;
}

[[nodiscard]] RealInterval pointAsin(
    const Rational& x,
    std::size_t precisionBits) {
    if (x < rational(-1) || x > rational(1))
        throw std::domain_error("Real asin input must be in [-1, 1]");
    if (x == rational(1))
        return halfPiInterval(precisionBits);
    if (x == rational(-1))
        return negate(halfPiInterval(precisionBits));
    if (x.isZero())
        return exactReal(0, precisionBits);
    if (x == rational(1, 2))
        return multiply(
            piInterval(precisionBits), exactReal(rational(1, 6), precisionBits), precisionBits);
    if (x == rational(-1, 2))
        return negate(multiply(
            piInterval(precisionBits), exactReal(rational(1, 6), precisionBits), precisionBits));

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified asin working precision is too large");
    const Rational radicand = rational(1) - x * x;
    const RealInterval root = encloseSqrt(exactReal(radicand, workBits), workBits).interval;
    const RealInterval ratio = divide(exactReal(x, workBits), root, workBits);
    return encloseAtan(ratio, workBits).interval.roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval pointAsinhPositive(
    const Rational& x,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified asinh working precision is too large");
    const RealInterval root = encloseSqrt(
        exactReal(x * x + rational(1), workBits), workBits).interval;
    const RealInterval sum = add(exactReal(x, workBits), root, workBits);
    return encloseLogPositive(sum, workBits).interval.roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval pointAsinh(
    const Rational& x,
    std::size_t precisionBits) {
    if (x.isZero())
        return exactReal(0, precisionBits);
    if (x.numerator().isNegative())
        return negate(pointAsinhPositive(-x, precisionBits));
    return pointAsinhPositive(x, precisionBits);
}

[[nodiscard]] RealInterval pointAcosh(
    const Rational& x,
    std::size_t precisionBits) {
    if (x < rational(1))
        throw std::domain_error("Real acosh input must be at least one");
    if (x == rational(1))
        return exactReal(0, precisionBits);

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified acosh working precision is too large");
    const RealInterval root = encloseSqrt(
        exactReal(x * x - rational(1), workBits), workBits).interval;
    const RealInterval sum = add(exactReal(x, workBits), root, workBits);
    return encloseLogPositive(sum, workBits).interval.roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval pointAtanh(
    const Rational& x,
    std::size_t precisionBits) {
    if (x <= rational(-1) || x >= rational(1))
        throw std::domain_error("Real atanh input must be in (-1, 1)");
    if (x.isZero())
        return exactReal(0, precisionBits);

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified atanh working precision is too large");
    const RealInterval numerator = exactReal(rational(1) + x, workBits);
    const RealInterval denominator = exactReal(rational(1) - x, workBits);
    const RealInterval ratio = divide(numerator, denominator, workBits);
    const RealInterval logarithm = encloseLogPositive(ratio, workBits).interval;
    return multiply(logarithm, exactReal(rational(1, 2), workBits), workBits)
        .roundedOutward(precisionBits);
}

[[nodiscard]] ComplexInterval quotientOrRetry(
    const ComplexInterval& numerator,
    const ComplexInterval& denominator,
    std::size_t precisionBits,
    const char* message) {
    if (denominator.containsZero())
        throw PrecisionInsufficient{message};
    return divide(numerator, denominator, precisionBits);
}

enum class PrincipalInverseKind {
    Asin,
    Acos,
    Atan,
    Asinh,
    Acosh,
    Atanh
};

[[nodiscard]] bool excludesZero(const RealInterval& value) {
    const BigFloat zero;
    return value.upper() < zero || value.lower() > zero;
}

[[nodiscard]] bool isAnalyticOnRectangle(
    PrincipalInverseKind kind,
    const ComplexInterval& value) {
    const Rational realLower = value.real().lower().toRational();
    const Rational realUpper = value.real().upper().toRational();
    const Rational imagLower = value.imaginary().lower().toRational();
    const Rational imagUpper = value.imaginary().upper().toRational();

    switch (kind) {
    case PrincipalInverseKind::Asin:
    case PrincipalInverseKind::Acos:
    case PrincipalInverseKind::Atanh:
        // principal cut: (-infinity,-1] U [1,infinity)
        return excludesZero(value.imaginary())
            || (realLower > rational(-1) && realUpper < rational(1));
    case PrincipalInverseKind::Acosh:
        // principal cut: (-infinity,1]
        return excludesZero(value.imaginary()) || realLower > rational(1);
    case PrincipalInverseKind::Atan:
    case PrincipalInverseKind::Asinh:
        // principal cut: -I[infinity,1] U I[1,infinity)
        return excludesZero(value.real())
            || (imagLower > rational(-1) && imagUpper < rational(1));
    }
    return false;
}

[[nodiscard]] Rational axisDistanceFromZero(const RealInterval& value) {
    const BigFloat zero;
    if (value.lower() > zero)
        return value.lower().toRational();
    if (value.upper() < zero)
        return -value.upper().toRational();
    return Rational{};
}

// 複素矩形から原点までの距離のcertified下界を返す。
// 0は「正の下界を証明できない」を表す。
[[nodiscard]] Rational complexMagnitudeLower(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const Rational dx = axisDistanceFromZero(value.real());
    const Rational dy = axisDistanceFromZero(value.imaginary());
    const Rational squared = dx * dx + dy * dy;
    if (squared.isZero())
        return Rational{};

    const RealInterval distance = encloseSqrt(
        exactReal(squared, precisionBits), precisionBits).interval;
    const Rational lower = distance.lower().toRational();
    return lower > Rational{} ? lower : Rational{};
}

[[nodiscard]] Rational complexRectangleRadiusUpper(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const Rational half = rational(1, 2);
    const Rational dx = (value.real().upper().toRational()
        - value.real().lower().toRational()) * half;
    const Rational dy = (value.imaginary().upper().toRational()
        - value.imaginary().lower().toRational()) * half;
    const Rational squared = dx * dx + dy * dy;
    if (squared.isZero())
        return Rational{};
    return encloseSqrt(exactReal(squared, precisionBits), precisionBits)
        .interval.upper().toRational();
}

[[nodiscard]] ComplexInterval rawPrincipalComplexAsinh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified complex asinh working precision is too large");
    const ComplexInterval square = multiply(value, value, workBits);
    const ComplexInterval radicand = add(square, oneComplex(workBits), workBits);
    const ComplexInterval root = enclosePrincipalComplexSqrt(radicand, workBits);
    return enclosePrincipalComplexLog(add(value, root, workBits), workBits)
        .interval.roundedOutward(precisionBits);
}

[[nodiscard]] std::optional<ComplexInterval> stableAcoshAcrossInteriorCut(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const Rational xLower = value.real().lower().toRational();
    const Rational xUpper = value.real().upper().toRational();
    if (xLower <= rational(-1) || xUpper >= rational(1))
        return std::nullopt;

    const BigFloat zero;
    const bool upperHalfPlane = value.imaginary().lower() > zero;
    const bool lowerHalfPlane = value.imaginary().upper() < zero;
    if (!upperHalfPlane && !lowerHalfPlane)
        return std::nullopt;

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified complex acosh working precision is too large");
    const RealInterval one = exactReal(1, workBits);
    const RealInterval four = exactReal(4, workBits);
    const RealInterval half = exactReal(rational(1, 2), workBits);
    const RealInterval x = value.real().roundedOutward(workBits);
    const RealInterval y = (upperHalfPlane ? value.imaginary() : negate(value.imaginary()))
        .roundedOutward(workBits);
    const RealInterval ySquared = multiply(y, y, workBits);

    const RealInterval xPlusOne = add(x, one, workBits);
    const RealInterval oneMinusX = subtract(one, x, workBits);
    const RealInterval rPlus = encloseSqrt(add(
        multiply(xPlusOne, xPlusOne, workBits), ySquared, workBits), workBits).interval;
    const RealInterval rMinus = encloseSqrt(add(
        multiply(oneMinusX, oneMinusX, workBits), ySquared, workBits), workBits).interval;

    // -1<x<1 では rPlus-(x+1), rMinus-(1-x) を有理化して加える。
    // rPlus+rMinus-2 を直接減算しないため，虚部が極小でも有効桁を失わない。
    const RealInterval deltaPlus = divide(
        ySquared, add(rPlus, xPlusOne, workBits), workBits);
    const RealInterval deltaMinus = divide(
        ySquared, add(rMinus, oneMinusX, workBits), workBits);
    const RealInterval delta = add(deltaPlus, deltaMinus, workBits);

    const RealInterval sinhA = multiply(
        encloseSqrt(multiply(delta, add(four, delta, workBits), workBits), workBits).interval,
        half, workBits);
    const RealInterval coshA = add(one, multiply(delta, half, workBits), workBits);
    const RealInterval cosB = divide(x, coshA, workBits);
    const Rational cosLower = cosB.lower().toRational();
    const Rational cosUpper = cosB.upper().toRational();
    if (cosLower < rational(-1) || cosUpper > rational(1))
        throw PrecisionInsufficient{
            "Acosh interior-cut angle cannot yet be certified at the current precision"};

    // Complex precisionは成分ごとの相対精度ではなく，複素数全体の相対精度である。
    // cutへ極端に近い点ではRe[acosh]とIm[acosh]の補正が要求精度より小さい。
    // この場合はacos(x/cosh(a))を直接高精度化せず，acos(x)からの変化量を
    // |acos'(t)|<=1/(1-q), q=max|t|<1 で厳密に包む。巨大な有理数へ成長した
    // x/cosh(a)をacosへ渡さないため，枝の安全性を保ったまま性能の崖を除く。
    constexpr std::size_t complexComponentGuardBits = 27;
    const RealInterval baseImaginary = encloseAcosRealRadian(x, workBits);
    const Rational baseImaginaryLower = baseImaginary.lower().toRational();
    const Rational sinhUpper = sinhA.upper().toRational();
    const bool negligibleRealPart = precisionBits > complexComponentGuardBits
        && baseImaginaryLower > Rational{}
        && sinhUpper <= baseImaginaryLower
            * binaryPrecisionThreshold(precisionBits - complexComponentGuardBits);

    std::optional<RealInterval> realPart;
    RealInterval imaginaryPart = baseImaginary;
    if (negligibleRealPart) {
        const auto absolute = [](const Rational& value) {
            return value.numerator().isNegative() ? -value : value;
        };
        const auto maximum = [](const Rational& lhs, const Rational& rhs) {
            return lhs < rhs ? rhs : lhs;
        };
        const Rational xLow = x.lower().toRational();
        const Rational xHigh = x.upper().toRational();
        const Rational q = maximum(
            maximum(absolute(xLow), absolute(xHigh)),
            maximum(absolute(cosLower), absolute(cosUpper)));

        if (q < rational(1)) {
            const Rational difference = maximum(
                absolute(cosLower - xHigh),
                absolute(cosUpper - xLow));
            const Rational angleError = difference / (rational(1) - q);
            imaginaryPart = RealInterval::fromRationalBounds(
                baseImaginary.lower().toRational() - angleError,
                baseImaginary.upper().toRational() + angleError,
                workBits);
            // 0<=asinh(t)<=t (t>=0) なので，微小実部はこの粗い区間で十分である。
            realPart = RealInterval::fromRationalBounds(Rational{}, sinhUpper, workBits);
        }
        else {
            imaginaryPart = encloseAcosRealRadian(cosB, workBits);
        }
    }
    else {
        imaginaryPart = encloseAcosRealRadian(cosB, workBits);
    }

    if (!realPart)
        realPart = encloseAsinhReal(sinhA, workBits);
    if (lowerHalfPlane)
        imaginaryPart = negate(imaginaryPart);

    return ComplexInterval{*realPart, imaginaryPart}.roundedOutward(precisionBits);
}

[[nodiscard]] ComplexInterval rawPrincipalComplexAcosh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified complex acosh working precision is too large");
    const ComplexInterval one = oneComplex(workBits);
    const ComplexInterval rootPlus = enclosePrincipalComplexSqrt(
        add(value, one, workBits), workBits);
    const ComplexInterval rootMinus = enclosePrincipalComplexSqrt(
        subtract(value, one, workBits), workBits);
    const ComplexInterval product = multiply(rootPlus, rootMinus, workBits);
    return enclosePrincipalComplexLog(add(value, product, workBits), workBits)
        .interval.roundedOutward(precisionBits);
}

[[nodiscard]] ComplexInterval rawPrincipalComplexAtanh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified complex atanh working precision is too large");
    const ComplexInterval one = oneComplex(workBits);
    const ComplexInterval numerator = add(one, value, workBits);
    const ComplexInterval denominator = subtract(one, value, workBits);
    const ComplexInterval left = enclosePrincipalComplexLog(numerator, workBits).interval;
    const ComplexInterval right = enclosePrincipalComplexLog(denominator, workBits).interval;
    return scaleComplex(subtract(left, right, workBits), rational(1, 2), workBits)
        .roundedOutward(precisionBits);
}

[[nodiscard]] ComplexInterval rawPrincipalInverse(
    PrincipalInverseKind kind,
    const ComplexInterval& value,
    std::size_t precisionBits) {
    switch (kind) {
    case PrincipalInverseKind::Asinh:
        return rawPrincipalComplexAsinh(value, precisionBits);
    case PrincipalInverseKind::Acosh:
        return rawPrincipalComplexAcosh(value, precisionBits);
    case PrincipalInverseKind::Atanh:
        return rawPrincipalComplexAtanh(value, precisionBits);
    case PrincipalInverseKind::Asin:
        return multiplyByNegativeI(
            rawPrincipalComplexAsinh(multiplyByI(value), precisionBits));
    case PrincipalInverseKind::Acos: {
        const std::size_t workBits = checkedPrecisionAdd(
            precisionBits, 16, "Certified complex acos working precision is too large");
        const ComplexInterval asinValue = rawPrincipalInverse(
            PrincipalInverseKind::Asin, value, workBits);
        const ComplexInterval halfPi = ComplexInterval::fromReal(halfPiInterval(workBits));
        return subtract(halfPi, asinValue, workBits).roundedOutward(precisionBits);
    }
    case PrincipalInverseKind::Atan:
        return multiplyByNegativeI(
            rawPrincipalComplexAtanh(multiplyByI(value), precisionBits));
    }
    throw std::logic_error("Unknown principal inverse function");
}

[[nodiscard]] Rational derivativeMagnitudeUpper(
    PrincipalInverseKind kind,
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const ComplexInterval one = oneComplex(precisionBits);

    Rational denominatorLower;
    switch (kind) {
    case PrincipalInverseKind::Asin:
    case PrincipalInverseKind::Acos:
    case PrincipalInverseKind::Atanh: {
        const ComplexInterval denominator = subtract(
            one, multiply(value, value, precisionBits), precisionBits);
        denominatorLower = complexMagnitudeLower(denominator, precisionBits);
        if ((kind == PrincipalInverseKind::Asin || kind == PrincipalInverseKind::Acos)
            && !denominatorLower.isZero()) {
            denominatorLower = encloseSqrt(
                exactReal(denominatorLower, precisionBits), precisionBits)
                .interval.lower().toRational();
        }
        break;
    }
    case PrincipalInverseKind::Atan:
    case PrincipalInverseKind::Asinh: {
        const ComplexInterval denominator = add(
            one, multiply(value, value, precisionBits), precisionBits);
        denominatorLower = complexMagnitudeLower(denominator, precisionBits);
        if (kind == PrincipalInverseKind::Asinh && !denominatorLower.isZero()) {
            denominatorLower = encloseSqrt(
                exactReal(denominatorLower, precisionBits), precisionBits)
                .interval.lower().toRational();
        }
        break;
    }
    case PrincipalInverseKind::Acosh: {
        const Rational plusLower = complexMagnitudeLower(
            add(value, one, precisionBits), precisionBits);
        const Rational minusLower = complexMagnitudeLower(
            subtract(value, one, precisionBits), precisionBits);
        if (plusLower.isZero() || minusLower.isZero())
            return Rational{};
        denominatorLower = encloseSqrt(
            exactReal(plusLower * minusLower, precisionBits), precisionBits)
            .interval.lower().toRational();
        break;
    }
    }

    if (denominatorLower <= Rational{})
        return Rational{};
    return rational(1) / denominatorLower;
}

[[nodiscard]] RealInterval intersectIntervals(
    const RealInterval& lhs,
    const RealInterval& rhs) {
    const BigFloat lower = lhs.lower() > rhs.lower() ? lhs.lower() : rhs.lower();
    const BigFloat upper = lhs.upper() < rhs.upper() ? lhs.upper() : rhs.upper();
    // 双方ともcertified enclosureなら本来必ず交差する。万一丸め境界で交差を
    // 証明できない場合は、狭めず元のenclosureを優先する。
    return lower <= upper ? RealInterval{lower, upper} : lhs;
}

[[nodiscard]] ComplexInterval tightenPrincipalInverse(
    PrincipalInverseKind kind,
    const ComplexInterval& value,
    const ComplexInterval& raw,
    std::size_t precisionBits) {
    if (!isAnalyticOnRectangle(kind, value))
        return raw;

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified inverse-function tightening precision is too large");
    const Rational radius = complexRectangleRadiusUpper(value, workBits);
    if (radius.isZero())
        return raw;

    const Rational derivative = derivativeMagnitudeUpper(kind, value, workBits);
    if (derivative.isZero())
        return raw;

    const Rational errorRadius = radius * derivative;
    const Rational candidateDiameter = rational(2) * errorRadius;
    const Rational rawRealWidth = raw.real().upper().toRational()
        - raw.real().lower().toRational();
    const Rational rawImaginaryWidth = raw.imaginary().upper().toRational()
        - raw.imaginary().lower().toRational();
    if (candidateDiameter >= rawRealWidth && candidateDiameter >= rawImaginaryWidth)
        return raw;

    const Rational centerReal = (value.real().lower().toRational()
        + value.real().upper().toRational()) * rational(1, 2);
    const Rational centerImaginary = (value.imaginary().lower().toRational()
        + value.imaginary().upper().toRational()) * rational(1, 2);
    const ComplexInterval centerInput{
        exactReal(centerReal, workBits),
        exactReal(centerImaginary, workBits)};
    const ComplexInterval centerValue = rawPrincipalInverse(kind, centerInput, workBits);

    // 凸な入力矩形上で |f(z)-f(c)| <= sup|f'| |z-c| を使う。
    // branch cutから離れ、導関数分母の正の下界を証明できた場合だけ適用する。
    const RealInterval error = exactReal(errorRadius, workBits);
    const RealInterval realLower = subtract(centerValue.real(), error, workBits);
    const RealInterval realUpper = add(centerValue.real(), error, workBits);
    const RealInterval imagLower = subtract(centerValue.imaginary(), error, workBits);
    const RealInterval imagUpper = add(centerValue.imaginary(), error, workBits);
    const ComplexInterval meanValue{
        RealInterval{realLower.lower(), realUpper.upper()},
        RealInterval{imagLower.lower(), imagUpper.upper()}};

    return ComplexInterval{
        intersectIntervals(raw.real(), meanValue.real()),
        intersectIntervals(raw.imaginary(), meanValue.imaginary())}
        .roundedOutward(precisionBits);
}

} // namespace

RealInterval encloseRealCubeRoot(
    const RealInterval& value,
    std::size_t precisionBits) {
    // cbrtは実軸全体で単調増加なので、両端を独立に囲えば区間像を得られる。
    const RealInterval lower = pointCubeRoot(value.lower().toRational(), precisionBits);
    const RealInterval upper = pointCubeRoot(value.upper().toRational(), precisionBits);
    return RealInterval{lower.lower(), upper.upper()};
}

ComplexInterval encloseComplexSinh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified complex sinh working precision is too large");
    const ComplexInterval positive = encloseComplexExp(value, workBits).interval;
    const ComplexInterval negative = encloseComplexExp(negate(value), workBits).interval;
    return scaleComplex(subtract(positive, negative, workBits), rational(1, 2), workBits)
        .roundedOutward(precisionBits);
}

ComplexInterval encloseComplexCosh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified complex cosh working precision is too large");
    const ComplexInterval positive = encloseComplexExp(value, workBits).interval;
    const ComplexInterval negative = encloseComplexExp(negate(value), workBits).interval;
    return scaleComplex(add(positive, negative, workBits), rational(1, 2), workBits)
        .roundedOutward(precisionBits);
}

ComplexInterval encloseComplexTanh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified complex tanh working precision is too large");
    const ComplexInterval numerator = encloseComplexSinh(value, workBits);
    const ComplexInterval denominator = encloseComplexCosh(value, workBits);
    return quotientOrRetry(
        numerator, denominator, workBits,
        "Complex tanh denominator cannot yet be proven nonzero")
        .roundedOutward(precisionBits);
}

ComplexInterval encloseComplexSinRadian(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    return multiplyByNegativeI(
        encloseComplexSinh(multiplyByI(value), precisionBits));
}

ComplexInterval encloseComplexCosRadian(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    return encloseComplexCosh(multiplyByI(value), precisionBits);
}

ComplexInterval encloseComplexTanRadian(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    return multiplyByNegativeI(
        encloseComplexTanh(multiplyByI(value), precisionBits));
}

RealInterval encloseSinhReal(
    const RealInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified sinh working precision is too large");
    const RealInterval positive = encloseExp(value, workBits).interval;
    const RealInterval negative = encloseExp(negate(value), workBits).interval;
    return multiply(
        subtract(positive, negative, workBits),
        exactReal(rational(1, 2), workBits), workBits).roundedOutward(precisionBits);
}

RealInterval encloseCoshReal(
    const RealInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified cosh working precision is too large");
    const RealInterval positive = encloseExp(value, workBits).interval;
    const RealInterval negative = encloseExp(negate(value), workBits).interval;
    return multiply(
        add(positive, negative, workBits),
        exactReal(rational(1, 2), workBits), workBits).roundedOutward(precisionBits);
}

RealInterval encloseTanhReal(
    const RealInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Certified tanh working precision is too large");
    return divide(
        encloseSinhReal(value, workBits),
        encloseCoshReal(value, workBits), workBits).roundedOutward(precisionBits);
}

RealInterval encloseAsinRealRadian(
    const RealInterval& value,
    std::size_t precisionBits) {
    const Rational lower = value.lower().toRational();
    const Rational upper = value.upper().toRational();
    if (lower < rational(-1) || upper > rational(1))
        throw std::domain_error("Real asin interval is outside [-1, 1]");

    const RealInterval lowerValue = pointAsin(lower, precisionBits);
    const RealInterval upperValue = pointAsin(upper, precisionBits);
    return RealInterval{lowerValue.lower(), upperValue.upper()};
}

RealInterval encloseAcosRealRadian(
    const RealInterval& value,
    std::size_t precisionBits) {
    // acos(x)=Pi/2-asin(x) かつ単調減少。
    const RealInterval asinValue = encloseAsinRealRadian(value, precisionBits);
    return subtract(halfPiInterval(precisionBits), asinValue, precisionBits);
}

ComplexInterval enclosePrincipalComplexAsinh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const ComplexInterval raw = rawPrincipalInverse(
        PrincipalInverseKind::Asinh, value, precisionBits);
    return tightenPrincipalInverse(
        PrincipalInverseKind::Asinh, value, raw, precisionBits);
}

ComplexInterval enclosePrincipalComplexAcosh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    if (const auto stable = stableAcoshAcrossInteriorCut(value, precisionBits))
        return *stable;
    const ComplexInterval raw = rawPrincipalInverse(
        PrincipalInverseKind::Acosh, value, precisionBits);
    return tightenPrincipalInverse(
        PrincipalInverseKind::Acosh, value, raw, precisionBits);
}

ComplexInterval enclosePrincipalComplexAtanh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const ComplexInterval raw = rawPrincipalInverse(
        PrincipalInverseKind::Atanh, value, precisionBits);
    return tightenPrincipalInverse(
        PrincipalInverseKind::Atanh, value, raw, precisionBits);
}

ComplexInterval enclosePrincipalComplexAsinRadian(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const ComplexInterval raw = rawPrincipalInverse(
        PrincipalInverseKind::Asin, value, precisionBits);
    return tightenPrincipalInverse(
        PrincipalInverseKind::Asin, value, raw, precisionBits);
}

ComplexInterval enclosePrincipalComplexAcosRadian(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const ComplexInterval raw = rawPrincipalInverse(
        PrincipalInverseKind::Acos, value, precisionBits);
    return tightenPrincipalInverse(
        PrincipalInverseKind::Acos, value, raw, precisionBits);
}

ComplexInterval enclosePrincipalComplexAtanRadian(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const ComplexInterval raw = rawPrincipalInverse(
        PrincipalInverseKind::Atan, value, precisionBits);
    return tightenPrincipalInverse(
        PrincipalInverseKind::Atan, value, raw, precisionBits);
}

RealInterval encloseAsinhReal(
    const RealInterval& value,
    std::size_t precisionBits) {
    const RealInterval lower = pointAsinh(value.lower().toRational(), precisionBits);
    const RealInterval upper = pointAsinh(value.upper().toRational(), precisionBits);
    return RealInterval{lower.lower(), upper.upper()};
}

RealInterval encloseAcoshReal(
    const RealInterval& value,
    std::size_t precisionBits) {
    if (value.lower().toRational() < rational(1))
        throw std::domain_error("Real acosh interval is below one");
    const RealInterval lower = pointAcosh(value.lower().toRational(), precisionBits);
    const RealInterval upper = pointAcosh(value.upper().toRational(), precisionBits);
    return RealInterval{lower.lower(), upper.upper()};
}

RealInterval encloseAtanhReal(
    const RealInterval& value,
    std::size_t precisionBits) {
    if (value.lower().toRational() <= rational(-1)
        || value.upper().toRational() >= rational(1))
        throw std::domain_error("Real atanh interval is outside (-1, 1)");
    const RealInterval lower = pointAtanh(value.lower().toRational(), precisionBits);
    const RealInterval upper = pointAtanh(value.upper().toRational(), precisionBits);
    return RealInterval{lower.lower(), upper.upper()};
}

} // namespace mmcal::approximation
