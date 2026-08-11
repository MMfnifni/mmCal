// 安定初等函数の保証付き評価
#include "certified_elementary_functions.hpp"

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

[[nodiscard]] std::size_t checkedAdd(
    std::size_t lhs,
    std::size_t rhs,
    const char* message) {
    if (rhs > std::numeric_limits<std::size_t>::max() - lhs)
        throw std::overflow_error(message);
    return lhs + rhs;
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
    const std::size_t workBits = checkedAdd(
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

    const std::size_t workBits = checkedAdd(
        precisionBits, 24, "Certified asin working precision is too large");
    const Rational radicand = rational(1) - x * x;
    const RealInterval root = encloseSqrt(exactReal(radicand, workBits), workBits).interval;
    const RealInterval ratio = divide(exactReal(x, workBits), root, workBits);
    return encloseAtan(ratio, workBits).interval.roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval pointAsinhPositive(
    const Rational& x,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
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

    const std::size_t workBits = checkedAdd(
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

    const std::size_t workBits = checkedAdd(
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
    const std::size_t workBits = checkedAdd(
        precisionBits, 24, "Certified complex sinh working precision is too large");
    const ComplexInterval positive = encloseComplexExp(value, workBits).interval;
    const ComplexInterval negative = encloseComplexExp(negate(value), workBits).interval;
    return scaleComplex(subtract(positive, negative, workBits), rational(1, 2), workBits)
        .roundedOutward(precisionBits);
}

ComplexInterval encloseComplexCosh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 24, "Certified complex cosh working precision is too large");
    const ComplexInterval positive = encloseComplexExp(value, workBits).interval;
    const ComplexInterval negative = encloseComplexExp(negate(value), workBits).interval;
    return scaleComplex(add(positive, negative, workBits), rational(1, 2), workBits)
        .roundedOutward(precisionBits);
}

ComplexInterval encloseComplexTanh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
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
    const std::size_t workBits = checkedAdd(
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
    const std::size_t workBits = checkedAdd(
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
    const std::size_t workBits = checkedAdd(
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
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Certified complex asinh working precision is too large");
    const ComplexInterval square = multiply(value, value, workBits);
    const ComplexInterval radicand = add(square, oneComplex(workBits), workBits);
    const ComplexInterval root = enclosePrincipalComplexSqrt(radicand, workBits);
    return enclosePrincipalComplexLog(add(value, root, workBits), workBits)
        .interval.roundedOutward(precisionBits);
}

ComplexInterval enclosePrincipalComplexAcosh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
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

ComplexInterval enclosePrincipalComplexAtanh(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Certified complex atanh working precision is too large");
    const ComplexInterval one = oneComplex(workBits);
    const ComplexInterval numerator = add(one, value, workBits);
    const ComplexInterval denominator = subtract(one, value, workBits);
    const ComplexInterval left = enclosePrincipalComplexLog(numerator, workBits).interval;
    const ComplexInterval right = enclosePrincipalComplexLog(denominator, workBits).interval;
    return scaleComplex(subtract(left, right, workBits), rational(1, 2), workBits)
        .roundedOutward(precisionBits);
}

ComplexInterval enclosePrincipalComplexAsinRadian(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    return multiplyByNegativeI(
        enclosePrincipalComplexAsinh(multiplyByI(value), precisionBits));
}

ComplexInterval enclosePrincipalComplexAcosRadian(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 16, "Certified complex acos working precision is too large");
    const ComplexInterval asinValue = enclosePrincipalComplexAsinRadian(value, workBits);
    const ComplexInterval halfPi = ComplexInterval::fromReal(halfPiInterval(workBits));
    return subtract(halfPi, asinValue, workBits).roundedOutward(precisionBits);
}

ComplexInterval enclosePrincipalComplexAtanRadian(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    return multiplyByNegativeI(
        enclosePrincipalComplexAtanh(multiplyByI(value), precisionBits));
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
