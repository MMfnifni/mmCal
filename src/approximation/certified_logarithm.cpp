// 対数函数の保証付き評価
#include "certified_logarithm.hpp"
#include "certified_precision.hpp"
#include "certified_sqrt.hpp"
#include "evaluation/evaluation_budget.hpp"

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


[[nodiscard]] std::size_t checkedShiftCount(std::uint64_t bits) {
    if (bits > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()))
        throw std::overflow_error("Binary scale exceeds addressable shift range");
    return static_cast<std::size_t>(bits);
}

[[nodiscard]] Rational scalePowerOfTwo(Rational value, std::int64_t exponent) {
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

struct LogSeriesSplit final {
    BigInt p{1};
    BigInt q{1};
    BigInt t{};
};

struct IntervalSeriesSplit final {
    RealInterval p;
    RealInterval t;
};

[[nodiscard]] bool preferExactSeries(const Rational& value) {
    constexpr std::size_t exactOperandBits = 16;
    return value.numerator().abs().bitLength() <= exactOperandBits
        && value.denominator().bitLength() <= exactOperandBits;
}

[[nodiscard]] IntervalSeriesSplit splitLogSeriesInterval(
    const RealInterval& tSquared,
    std::uint64_t begin,
    std::uint64_t end,
    std::size_t precisionBits) {
    if (end == begin + 1) {
        if (begin > (std::numeric_limits<std::uint64_t>::max() - 1) / 2)
            throw std::overflow_error("Certified log series index is too large");
        const Rational oddRatio{
            BigInt::fromUnsigned(2 * begin - 1),
            BigInt::fromUnsigned(2 * begin + 1)};
        const RealInterval ratio = multiply(
            tSquared, RealInterval::fromRational(oddRatio, precisionBits), precisionBits);
        return IntervalSeriesSplit{ratio, ratio};
    }

    const std::uint64_t middle = begin + (end - begin) / 2;
    IntervalSeriesSplit left = splitLogSeriesInterval(
        tSquared, begin, middle, precisionBits);
    IntervalSeriesSplit right = splitLogSeriesInterval(
        tSquared, middle, end, precisionBits);
    return IntervalSeriesSplit{
        multiply(left.p, right.p, precisionBits),
        add(left.t, multiply(left.p, right.t, precisionBits), precisionBits)
    };
}

// u_n=t^(2n+1)/(2n+1) は
//   u_n/u_(n-1) = t^2 * (2n-1)/(2n+1)
// を満たす。t=a/bとして、このratio列をbinary splittingする。
[[nodiscard]] LogSeriesSplit splitLogSeries(
    const BigInt& numeratorSquared,
    const BigInt& denominatorSquared,
    std::uint64_t begin,
    std::uint64_t end) {
    if (end == begin + 1) {
        if (begin > (std::numeric_limits<std::uint64_t>::max() - 1) / 2)
            throw std::overflow_error("Certified log series index is too large");
        const std::uint64_t previousOdd = 2 * begin - 1;
        const std::uint64_t nextOdd = 2 * begin + 1;
        const BigInt p = numeratorSquared * BigInt::fromUnsigned(previousOdd);
        const BigInt q = denominatorSquared * BigInt::fromUnsigned(nextOdd);
        return LogSeriesSplit{p, q, p};
    }

    const std::uint64_t middle = begin + (end - begin) / 2;
    LogSeriesSplit left = splitLogSeries(
        numeratorSquared, denominatorSquared, begin, middle);
    LogSeriesSplit right = splitLogSeries(
        numeratorSquared, denominatorSquared, middle, end);

    LogSeriesSplit result;
    result.p = left.p * right.p;
    result.q = left.q * right.q;
    result.t = left.t * right.q + left.p * right.t;
    return result;
}

[[nodiscard]] SeriesResult encloseLogMantissaReducedInterval(
    const Rational& m,
    std::size_t precisionBits) {
    constexpr std::size_t sqrtReductions = 16;
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified log precision is too large");

    RealInterval reduced = RealInterval::fromRational(m, workBits);
    for (std::size_t i = 0; i < sqrtReductions; ++i) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        reduced = encloseSqrt(reduced, workBits).interval;
    }

    const RealInterval one = RealInterval::fromRational(rational(1), workBits);
    const RealInterval t = divide(
        subtract(reduced, one, workBits),
        add(reduced, one, workBits),
        workBits);
    const RealInterval tSquared = multiply(t, t, workBits);

    // 1<=m<=2なら、sqrtをs回適用した後は
    //   t=(root-1)/(root+1) < 2^-(s+1)
    // と保守的に押さえられる。元のlogへ戻す2^(s+1)倍も含め、
    // tailが要求bitを下回る項数を固定回数ではなく精度から決める。
    const std::size_t targetBits = checkedPrecisionAdd(
        precisionBits, 28, "Certified log precision is too large");
    constexpr std::size_t reductionDenominator = 2 * (sqrtReductions + 1);
    if (targetBits > std::numeric_limits<std::size_t>::max() - reductionDenominator)
        throw std::overflow_error("Certified log precision is too large");
    const std::size_t terms = (targetBits + reductionDenominator - 1)
        / reductionDenominator + 2;
    if (terms == 0 || terms >= std::numeric_limits<std::uint64_t>::max())
        throw std::overflow_error("Certified log series index is too large");

    RealInterval normalizedSum = one;
    RealInterval lastRatioProduct = one;
    if (terms > 1) {
        const IntervalSeriesSplit split = splitLogSeriesInterval(
            tSquared, 1, static_cast<std::uint64_t>(terms), workBits);
        normalizedSum = add(normalizedSum, split.t, workBits);
        lastRatioProduct = split.p;
    }

    const RealInterval partialHalf = multiply(t, normalizedSum, workBits);
    const std::uint64_t nextIndex = static_cast<std::uint64_t>(terms);
    if (nextIndex > (std::numeric_limits<std::uint64_t>::max() - 1) / 2)
        throw std::overflow_error("Certified log series index is too large");
    const Rational nextOddRatio{
        BigInt::fromUnsigned(2 * nextIndex - 1),
        BigInt::fromUnsigned(2 * nextIndex + 1)};
    const RealInterval nextRatio = multiply(
        tSquared,
        RealInterval::fromRational(nextOddRatio, workBits),
        workBits);
    const RealInterval nextTerm = multiply(
        multiply(t, lastRatioProduct, workBits), nextRatio, workBits);
    const RealInterval tailHalf = divide(
        nextTerm,
        subtract(one, tSquared, workBits),
        workBits);
    const RealInterval nonNegativeTail = RealInterval::fromRationalBounds(
        rational(0), tailHalf.upper().toRational(), workBits);

    BigInt scaleInteger{1};
    scaleInteger <<= sqrtReductions + 1;
    const RealInterval scale = RealInterval::fromRational(
        Rational{std::move(scaleInteger)}, workBits);
    const RealInterval lower = multiply(partialHalf, scale, workBits);
    const RealInterval upper = multiply(
        add(partialHalf, nonNegativeTail, workBits), scale, workBits);
    return SeriesResult{
        RealInterval{lower.lower(), upper.upper()}.roundedOutward(precisionBits),
        terms
    };
}

[[nodiscard]] std::size_t requiredLogTerms(std::size_t precisionBits) {
    // 1<=m<=2では t=(m-1)/(m+1)<=1/3。
    // M項採用後、log(m)=2*sum u_n のtailは
    //   2*R_M <= 2*u_M/(1-t^2) < 2^(-3M)
    // と保守的に抑えられる。したがってprecision+guardを3で割るだけで
    // correctnessに依存しない十分な項数を決められる。
    const std::size_t targetBits = checkedPrecisionAdd(
        precisionBits, 28, "Certified log precision is too large");
    if (targetBits > std::numeric_limits<std::size_t>::max() - 2)
        throw std::overflow_error("Certified log precision is too large");
    return (targetBits + 2) / 3;
}


// 1 <= m <= 2 に対して log(m)=2*atanh((m-1)/(m+1)) を保証付き評価する。
[[nodiscard]] SeriesResult encloseLogMantissa(
    const Rational& m,
    std::size_t precisionBits) {
    if (m < rational(1) || m > rational(2))
        throw std::invalid_argument("Log mantissa must be in [1, 2]");
    if (m == rational(1))
        return SeriesResult{RealInterval::fromRational(rational(0), precisionBits), 0};

    const Rational t = (m - rational(1)) / (m + rational(1));
    const Rational tSquared = t * t;
    const std::size_t terms = requiredLogTerms(precisionBits);
    if (terms == 0)
        throw std::logic_error("Certified log term count is zero");
    if (terms >= std::numeric_limits<std::uint64_t>::max())
        throw std::overflow_error("Certified log series index is too large");

    const std::uint64_t nextIndex = static_cast<std::uint64_t>(terms);
    if (nextIndex > (std::numeric_limits<std::uint64_t>::max() - 1) / 2)
        throw std::overflow_error("Certified log series index is too large");

    if (preferExactSeries(t)) {
        Rational normalizedSum = rational(1);
        Rational lastRatioProduct = rational(1);
        if (terms > 1) {
            const LogSeriesSplit split = splitLogSeries(
                tSquared.numerator(), tSquared.denominator(),
                1, static_cast<std::uint64_t>(terms));
            normalizedSum += Rational{split.t, split.q};
            lastRatioProduct = Rational{split.p, split.q};
        }

        const Rational partialHalf = t * normalizedSum;
        const BigInt nextNumerator = tSquared.numerator()
            * BigInt::fromUnsigned(2 * nextIndex - 1);
        const BigInt nextDenominator = tSquared.denominator()
            * BigInt::fromUnsigned(2 * nextIndex + 1);
        const Rational nextTerm = t * lastRatioProduct
            * Rational{nextNumerator, nextDenominator};
        const Rational tailHalf = nextTerm / (rational(1) - tSquared);

        const Rational lower = partialHalf * rational(2);
        const Rational upper = (partialHalf + tailHalf) * rational(2);
        return SeriesResult{
            RealInterval::fromRationalBounds(lower, upper, precisionBits),
            terms
        };
    }

    // 巨大Rationalではexact splitのa^(2N), b^(2N)が必要precisionを超えて
    // 巨大化する。旧fallbackは同じtのまま固定precision interval splitしていたが、
    // 項数自体は減らないため5000桁級で秒単位まで伸びた。
    // 16回のcertified sqrtでmantissaを1へ近づけてから級数を評価し、
    // 最後に2^16倍してlog(m)へ戻す。
    return encloseLogMantissaReducedInterval(m, precisionBits);
}

[[nodiscard]] SeriesResult encloseLogPoint(
    const Rational& x,
    std::size_t precisionBits,
    const SeriesResult* cachedLog2) {
    if (x <= rational(0))
        throw std::domain_error("Real Log requires a positive argument");

    std::int64_t k = numeric::detail::floorLog2PositiveRatio(
        x.numerator(), x.denominator());

    // 1<=m<2へbinary normalizationし、sqrt(2)近傍の99/70を境にさらに2で割る。
    Rational m = scalePowerOfTwo(x, -k);
    if (m > Rational{BigInt{99}, BigInt{70}}) {
        m /= rational(2);
        ++k;
    }

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

    // log(2)は同一enclosure内で一度だけ計算する。
    const SeriesResult log2 = encloseLogMantissa(rational(2), precisionBits);

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
