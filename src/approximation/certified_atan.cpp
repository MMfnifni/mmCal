#include "certified_atan.hpp"
#include "certified_precision.hpp"

#include "certified_constants.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace mmcal::approximation {
namespace {

using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}



struct PointAtanResult final {
    Rational lower;
    Rational upper;
    std::size_t termsUsed = 0;
};

[[nodiscard]] PointAtanResult negateResult(PointAtanResult value) {
    return PointAtanResult{-value.upper, -value.lower, value.termsUsed};
}

[[nodiscard]] PointAtanResult addResults(
    const PointAtanResult& lhs,
    const PointAtanResult& rhs) {
    return PointAtanResult{
        lhs.lower + rhs.lower,
        lhs.upper + rhs.upper,
        lhs.termsUsed + rhs.termsUsed
    };
}

// |x| <= 1/2 に縮約済みの一点 atan(x) を、交代級数だけで厳密に囲う。
//
//   atan(x) = x - x^3/3 + x^5/5 - ...
//
// 0 <= x <= 1/2 では項の絶対値が単調減少するので、交代級数の基本定理から
// 真値は「現在の部分和」と「次項まで加えた部分和」の間に必ず存在する。
//
// 入力xはRationalとしてexactだが、級数の全中間値までnormalized Rationalで
// 保持する必要はない。高precision dyadic xでは分母が項ごとに巨大化し、Arg/Log
// 全体のperformance cliffになる。固定working precisionのoutward intervalで
// 隣接部分和を包含し、そのhullを返せば交代級数の証明はそのまま保てる。
[[nodiscard]] PointAtanResult encloseAtanSeriesPositive(
    const Rational& x,
    std::size_t precisionBits) {
    if (x < rational(0) || x > rational(1, 2))
        throw std::invalid_argument("Certified atan series argument must be in [0, 1/2]");
    if (x.isZero())
        return PointAtanResult{rational(0), rational(0), 1};

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified atan working precision is too large");
    const Rational threshold = binaryPrecisionThreshold(checkedPrecisionAdd(
        precisionBits, 16, "Certified atan precision is too large"));
    const RealInterval xInterval = RealInterval::fromRational(x, workBits);
    const RealInterval xSquared = multiply(xInterval, xInterval, workBits);

    RealInterval power = xInterval; // x^(2k+1)
    RealInterval sum = xInterval;   // k=0 の部分和
    bool nextNegative = true;
    std::uint64_t odd = 1;
    std::size_t termsUsed = 1;

    for (;;) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        if (odd > std::numeric_limits<std::uint64_t>::max() - 2)
            throw std::overflow_error("Certified atan series index overflow");
        odd += 2;
        power = multiply(power, xSquared, workBits);

        if (odd > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
            throw std::overflow_error("Certified atan series index exceeds BigInt small-integer range");

        const Rational inverseOdd{BigInt{1}, BigInt{static_cast<std::int64_t>(odd)}};
        const RealInterval magnitude = multiply(
            power, RealInterval::fromRational(inverseOdd, workBits), workBits);
        const RealInterval nextSum = nextNegative
            ? subtract(sum, magnitude, workBits)
            : add(sum, magnitude, workBits);

        // magnitudeの上端が閾値以下なら、真の次項も必ず閾値以下である。
        // 真値はexactな隣接部分和の間にあり、それぞれをsum/nextSumが包含するため、
        // interval hullは丸めを含めても厳密なatan(x) enclosureになる。
        if (magnitude.upper().toRational() <= threshold) {
            const RealInterval enclosure = hull(sum, nextSum).roundedOutward(precisionBits);
            return PointAtanResult{
                enclosure.lower().toRational(),
                enclosure.upper().toRational(),
                termsUsed};
        }

        sum = nextSum;
        nextNegative = !nextNegative;
        if (termsUsed != std::numeric_limits<std::size_t>::max())
            ++termsUsed;
    }
}

[[nodiscard]] PointAtanResult encloseAtanPointPositive(
    const Rational& x,
    std::size_t precisionBits) {
    if (x < rational(0))
        throw std::invalid_argument("Positive atan reducer received a negative argument");

    // x > 1 では atan(x) = Pi/2 - atan(1/x)。
    // これで残る引数は (0,1]。Piもcertified intervalなので、変換自体に
    // machine constantを一切持ち込まない。
    if (x > rational(1)) {
        const auto reciprocal = encloseAtanPointPositive(rational(1) / x, precisionBits);
        const CertifiedConstantResult pi = enclosePi(checkedPrecisionAdd(
            precisionBits, 16, "Certified atan Pi precision is too large"));
        const Rational piLowerHalf = pi.interval.lower().toRational() / rational(2);
        const Rational piUpperHalf = pi.interval.upper().toRational() / rational(2);
        return PointAtanResult{
            piLowerHalf - reciprocal.upper,
            piUpperHalf - reciprocal.lower,
            reciprocal.termsUsed + pi.termsUsed
        };
    }

    // 0 <= x <= 1/2 ならそのまま高速な交代級数へ入れる。
    if (x <= rational(1, 2))
        return encloseAtanSeriesPositive(x, precisionBits);

    // 1/2 < x <= 1 では加法定理
    //
    //   atan(x) = atan(1/2) + atan((x-1/2)/(1+x/2))
    //
    // を使う。変換後 t=(2x-1)/(x+2) は 0 < t <= 1/3 なので、両方のatanを
    // |argument|<=1/2 の同じreference級数で評価できる。
    const Rational half = rational(1, 2);
    const Rational transformed = (x - half) / (rational(1) + x * half);
    return addResults(
        encloseAtanSeriesPositive(half, precisionBits),
        encloseAtanSeriesPositive(transformed, precisionBits));
}

[[nodiscard]] PointAtanResult encloseAtanPoint(
    const Rational& x,
    std::size_t precisionBits) {
    if (x < rational(0))
        return negateResult(encloseAtanPointPositive(-x, precisionBits));
    return encloseAtanPointPositive(x, precisionBits);
}

} // namespace

CertifiedAtanResult encloseAtan(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified atan precision must be at least one bit");

    // atanは実数全体で厳密に単調増加。
    // よって atan([a,b]) = [atan(a),atan(b)] を端点ごとのcertified enclosureで作れる。
    const PointAtanResult lower = encloseAtanPoint(
        input.lower().toRational(), precisionBits);
    const PointAtanResult upper = encloseAtanPoint(
        input.upper().toRational(), precisionBits);

    return CertifiedAtanResult{
        RealInterval::fromRationalBounds(lower.lower, upper.upper, precisionBits),
        lower.termsUsed + upper.termsUsed
    };
}

} // namespace mmcal::approximation
