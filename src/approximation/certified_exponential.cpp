// 指数函数の保証付き評価
#include "certified_exponential.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstdint>
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

[[nodiscard]] RealInterval exactIntegerInterval(
    std::int64_t value,
    std::size_t precisionBits) {
    return RealInterval::fromRational(rational(value), precisionBits);
}

struct PointExpResult final {
    RealInterval interval;
    std::size_t termsUsed = 0;
    std::size_t squarings = 0;
};

// 0 <= x <= 1/2 の一点をTaylor級数で囲う。
//   exp(x) = sum_{n=0}^N x^n/n! + R_N
// x>=0なので全項は非負。次項 t_(N+1) 以降の比は
//   t_(m+1)/t_m = x/(m+1) <= x/(N+2) = r < 1
// だから残り全部は幾何級数で
//   0 <= R_N <= t_(N+1)/(1-r)
// と厳密に上から押さえられる。部分和の各演算はRealIntervalで外向き丸めするため、打切り誤差とBigFloat丸め誤差の両方を同じ包含区間へ入れられる。
[[nodiscard]] PointExpResult encloseExpSmallNonNegative(
    const Rational& x,
    std::size_t precisionBits) {
    if (x < rational(0) || x > rational(1, 2))
        throw std::invalid_argument("Reduced exponential argument must be in [0, 1/2]");

    const RealInterval one = exactIntegerInterval(1, precisionBits);
    if (x.isZero())
        return PointExpResult{one, 1, 0};

    const RealInterval xInterval = RealInterval::fromRational(x, precisionBits);
    RealInterval term = one;
    RealInterval sum = one;
    std::size_t termsUsed = 1;

    const Rational threshold = binaryThreshold(checkedAdd(
        precisionBits, 24, "Certified exp precision is too large"));

    for (std::uint64_t n = 1;; ++n) {
        // 級数indexは通常ごく小さいが、BigIntの小整数constructorはint64_t。
        // 極端なprecision指定でも符号付きcastをwrapさせず、明示的に停止する。
        if (n > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max() - 2))
            throw std::overflow_error("Certified exp series index exceeds BigInt small-integer range");
        const RealInterval divisor = RealInterval::fromRational(
            Rational{BigInt{static_cast<std::int64_t>(n)}}, precisionBits);
        term = divide(multiply(term, xInterval, precisionBits), divisor, precisionBits);
        sum = add(sum, term, precisionBits);
        ++termsUsed;

        const std::uint64_t nextIndex = n + 1;
        const RealInterval nextDivisor = RealInterval::fromRational(
            Rational{BigInt{static_cast<std::int64_t>(nextIndex)}}, precisionBits);
        const RealInterval nextTerm = divide(
            multiply(term, xInterval, precisionBits), nextDivisor, precisionBits);

        // tail <= nextTerm / (1 - x/(n+2))。
        // nextTerm.upper() は外向き丸め済みなので、そこから作るtailBoundも真の剰余を必ず上から押さえる。
        const Rational ratio = x / Rational{BigInt{static_cast<std::int64_t>(n + 2)}};
        const Rational geometricFactor = rational(1) / (rational(1) - ratio);
        const Rational tailBound = nextTerm.upper().toRational() * geometricFactor;

        if (tailBound <= threshold) {
            const RealInterval tail = RealInterval::fromRationalBounds(
                rational(0), tailBound, precisionBits);
            return PointExpResult{add(sum, tail, precisionBits), termsUsed, 0};
        }

        if (n == std::numeric_limits<std::uint64_t>::max())
            throw std::overflow_error("Certified exp series iteration overflow");
    }
}

[[nodiscard]] PointExpResult encloseExpPoint(
    Rational x,
    std::size_t precisionBits) {
    const bool negative = x < rational(0);
    if (negative)
        x = -x;

    // exp(x) = exp(x / 2^k)^(2^k)。
    // |x/2^k| <= 1/2 までexact Rationalで縮約し、Taylorの収束を安定化する。最後の復元は区間二乗なので、途中丸めも包含保証を失わない。
    std::size_t squarings = 0;
    while (x > rational(1, 2)) {
        x /= rational(2);
        ++squarings;
    }

    PointExpResult result = encloseExpSmallNonNegative(x, precisionBits);
    result.squarings = squarings;
    for (std::size_t i = 0; i < squarings; ++i)
        result.interval = multiply(result.interval, result.interval, precisionBits);

    if (!negative)
        return result;

    // exp(-x)=1/exp(x)。正数区間なので0を跨がず、安全に逆数区間を作れる。
    result.interval = divide(exactIntegerInterval(1, precisionBits), result.interval, precisionBits);
    return result;
}

} // namespace

CertifiedExponentialResult encloseExp(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified exp precision must be at least one bit");

    // expは実軸上で厳密に単調増加するため、
    //   exp([a,b]) = [exp(a), exp(b)]
    // の包含を端点別のcertified計算だけで構成できる。
    const PointExpResult lower = encloseExpPoint(input.lower().toRational(), precisionBits);

    /*
    旧実装ではpoint intervalでもlower/upperを同じ値から二度Taylor評価していた。
    N[E,n]やN[exp[1],n]はexact pointなので完全な重複計算になる。
    pointなら一度のcertified enclosureをそのまま返し、区間入力だけ上端を別評価する。
    */
    if (input.isPoint())
        return CertifiedExponentialResult{
            lower.interval, lower.termsUsed, lower.squarings};

    const PointExpResult upper = encloseExpPoint(input.upper().toRational(), precisionBits);
    return CertifiedExponentialResult{
        RealInterval{lower.interval.lower(), upper.interval.upper()},
        lower.termsUsed + upper.termsUsed,
        lower.squarings + upper.squarings
    };
}

} // namespace mmcal::approximation
