// 対数函数の保証付き評価
#include "certified_logarithm.hpp"

#include "numeric/detail/binary_scale.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

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

[[nodiscard]] std::size_t checkedShiftCount(std::uint64_t bits) {
    if (bits > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()))
        throw std::overflow_error("Binary scale exceeds addressable shift range");
    return static_cast<std::size_t>(bits);
}

[[nodiscard]] Rational scalePowerOfTwo(Rational value, std::int64_t exponent) {
    // value * 2^exponent をexact Rationalで作る。
    //
    // BigIntのshift countはsize_tなので、int64_tから無条件にcastしない。
    // x64では通常同じ幅だが、型の意味としては別物であり、32-bit buildでもsilent truncationを起こさないよう明示的に範囲検査する。
    BigInt numerator = value.numerator();
    BigInt denominator = value.denominator();
    if (exponent >= 0) {
        numerator <<= checkedShiftCount(static_cast<std::uint64_t>(exponent));
    }
    else {
        if (exponent == std::numeric_limits<std::int64_t>::min())
            throw std::overflow_error("Binary scale exponent is too small");
        denominator <<= checkedShiftCount(static_cast<std::uint64_t>(-exponent));
    }
    return Rational{std::move(numerator), std::move(denominator)};
}

[[nodiscard]] RealInterval exactIntegerInterval(
    const BigInt& value,
    std::size_t precisionBits) {
    return RealInterval::fromRational(Rational{value}, precisionBits);
}

struct SeriesResult final {
    RealInterval interval;
    std::size_t termsUsed = 0;
};

// 1 <= m <= 2 に対し、
//   log(m) = 2 atanh(t),  t=(m-1)/(m+1)
//          = 2 (t + t^3/3 + t^5/5 + ...)
// を使う。この範囲なら 0 <= t <= 1/3 なので非常に安定して収束する。全項が非負で、次項以降の比は t^2 以下だから、次項Aから先のtailは
//   tail <= A / (1 - t^2)
// と厳密に上から押さえられる。
[[nodiscard]] SeriesResult encloseLogMantissa(
    const Rational& m,
    std::size_t precisionBits) {
    if (m < rational(1) || m > rational(2))
        throw std::invalid_argument("Log mantissa must be in [1, 2]");
    if (m == rational(1))
        return SeriesResult{RealInterval::fromRational(rational(0), precisionBits), 0};

    const Rational tExact = (m - rational(1)) / (m + rational(1));
    const Rational tSquaredExact = tExact * tExact;
    const RealInterval t = RealInterval::fromRational(tExact, precisionBits);
    const RealInterval tSquared = RealInterval::fromRational(tSquaredExact, precisionBits);

    RealInterval power = t; // t^(2n+1)
    RealInterval sum = RealInterval::fromRational(rational(0), precisionBits);
    std::uint64_t odd = 1;
    std::size_t termsUsed = 0;

    const Rational threshold = binaryThreshold(checkedAdd(
        precisionBits, 24, "Certified log precision is too large"));
    const Rational geometricFactor = rational(1) / (rational(1) - tSquaredExact);

    for (;;) {
        // BigIntの小整数constructorはint64_tなので、理論上そこを越えるほどの項数を要求された場合はwrapさせず、資源上限として明示的に止める。
        if (odd > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
            throw std::overflow_error("Certified log series index exceeds BigInt small-integer range");
        const RealInterval divisor = RealInterval::fromRational(
            Rational{BigInt{static_cast<std::int64_t>(odd)}}, precisionBits);
        const RealInterval term = divide(power, divisor, precisionBits);
        sum = add(sum, term, precisionBits);
        ++termsUsed;

        if (odd > std::numeric_limits<std::uint64_t>::max() - 2)
            throw std::overflow_error("Certified log series iteration overflow");
        const std::uint64_t nextOdd = odd + 2;
        if (nextOdd > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
            throw std::overflow_error("Certified log series index exceeds BigInt small-integer range");
        const RealInterval nextPower = multiply(power, tSquared, precisionBits);
        const RealInterval nextDivisor = RealInterval::fromRational(
            Rational{BigInt{static_cast<std::int64_t>(nextOdd)}}, precisionBits);
        const RealInterval nextTerm = divide(nextPower, nextDivisor, precisionBits);
        const Rational tailBound = nextTerm.upper().toRational() * geometricFactor;

        if (tailBound <= threshold) {
            const RealInterval tail = RealInterval::fromRationalBounds(
                rational(0), tailBound, precisionBits);
            const RealInterval two = RealInterval::fromRational(rational(2), precisionBits);
            return SeriesResult{
                multiply(add(sum, tail, precisionBits), two, precisionBits),
                termsUsed
            };
        }

        power = nextPower;
        odd = nextOdd;
    }
}

[[nodiscard]] SeriesResult encloseLogPoint(
    const Rational& x,
    std::size_t precisionBits,
    const SeriesResult* cachedLog2) {
    if (x <= rational(0))
        throw std::domain_error("Real Log requires a positive argument");

    std::int64_t k = numeric::detail::floorLog2PositiveRatio(
        x.numerator(), x.denominator());

    // 通常のbinary normalizationでは 1 <= m < 2 だが、atanh級数の |t|=(m-1)/(m+1) を小さくするため、
    // sqrt(2)に近いexact rational threshold 99/70 を使って上半分をさらに2で割る。恒等式 x=m*2^k はexactのままで、近似thresholdは値そのものには入らない。
    Rational m = scalePowerOfTwo(x, -k);
    if (m > Rational{BigInt{99}, BigInt{70}}) {
        m /= rational(2);
        ++k;
    }

    // m<1なら log(m)=-log(1/m)。上のrange reductionにより1/mも約sqrt(2)以下なので、級数の最大|t|は約0.172に抑えられる。
    bool negateMantissa = false;
    if (m < rational(1)) {
        m = rational(1) / m;
        negateMantissa = true;
    }
    SeriesResult logM = encloseLogMantissa(m, precisionBits);
    if (negateMantissa)
        logM.interval = approximation::negate(logM.interval);
    if (k == 0)
        return logM;

    if (!cachedLog2)
        throw std::logic_error("Certified log2 cache is missing");
    const RealInterval kInterval = exactIntegerInterval(BigInt{k}, precisionBits);
    const RealInterval scaledLog2 = multiply(cachedLog2->interval, kInterval, precisionBits);
    logM.interval = add(logM.interval, scaledLog2, precisionBits);
    logM.termsUsed += cachedLog2->termsUsed;
    return logM;
}

} // namespace

CertifiedLogarithmResult encloseLogPositive(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified Log precision must be at least one bit");

    const BigFloat zero;
    if (input.lower() <= zero)
        throw std::domain_error("Real Log interval must be strictly positive");

    // log(2)は同一enclosure内で一度だけ計算する。旧実装はlower/upperの各point評価で重複していた。
    const SeriesResult log2 = encloseLogMantissa(rational(2), precisionBits);

    // point intervalなら同じ級数を二度評価しない。
    const SeriesResult lower = encloseLogPoint(
        input.lower().toRational(), precisionBits, &log2);
    if (input.isPoint())
        return CertifiedLogarithmResult{lower.interval, lower.termsUsed};

    // logは正実軸上で単調増加。
    const SeriesResult upper = encloseLogPoint(
        input.upper().toRational(), precisionBits, &log2);
    return CertifiedLogarithmResult{
        RealInterval{lower.interval.lower(), upper.interval.upper()},
        lower.termsUsed + upper.termsUsed
    };
}

} // namespace mmcal::approximation
