#include "certified_atan.hpp"

#include "certified_constants.hpp"
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

[[nodiscard]] Rational binaryThreshold(std::size_t bits) {
    BigInt denominator{1};
    denominator <<= bits;
    return Rational{BigInt{1}, std::move(denominator)};
}

[[nodiscard]] std::size_t checkedAdd(
    std::size_t lhs,
    std::size_t rhs,
    const char* message) {
    if (rhs > std::numeric_limits<std::size_t>::max() - lhs)
        throw std::overflow_error(message);
    return lhs + rhs;
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
// つまり誤差の大きさを経験的に推測する必要がない。
//
// ここはreference実装としてRationalを完全exactに保つ。atanはArg/Log/Powerの
// branchを決める根幹なので、まず証明が最も単純な算法を採用する。
[[nodiscard]] PointAtanResult encloseAtanSeriesPositive(
    const Rational& x,
    std::size_t precisionBits) {
    if (x < rational(0) || x > rational(1, 2))
        throw std::invalid_argument("Certified atan series argument must be in [0, 1/2]");
    if (x.isZero())
        return PointAtanResult{rational(0), rational(0), 1};

    const Rational threshold = binaryThreshold(checkedAdd(
        precisionBits, 16, "Certified atan precision is too large"));
    const Rational xSquared = x * x;

    Rational power = x;   // x^(2k+1)
    Rational sum = x;     // k=0 の部分和
    bool nextNegative = true;
    std::uint64_t odd = 1;
    std::size_t termsUsed = 1;

    for (;;) {
        if (odd > std::numeric_limits<std::uint64_t>::max() - 2)
            throw std::overflow_error("Certified atan series index overflow");
        odd += 2;
        power *= xSquared;

        if (odd > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
            throw std::overflow_error("Certified atan series index exceeds BigInt small-integer range");

        const Rational magnitude = power / Rational{BigInt{static_cast<std::int64_t>(odd)}};
        const Rational nextSum = nextNegative ? sum - magnitude : sum + magnitude;

        // 交代級数では真値が隣接する二つの部分和の間にある。
        // magnitude <= threshold まで狭まれば、この区間は要求binary precisionより
        // 十分細かい。最終10進丸めの一意性はさらに上位層が検査する。
        if (magnitude <= threshold) {
            const Rational lower = sum < nextSum ? sum : nextSum;
            const Rational upper = sum < nextSum ? nextSum : sum;
            return PointAtanResult{lower, upper, termsUsed};
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
        const CertifiedConstantResult pi = enclosePi(checkedAdd(
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
