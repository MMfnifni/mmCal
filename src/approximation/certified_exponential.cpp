// 指数函数の保証付き評価
#include "certified_exponential.hpp"
#include "certified_precision.hpp"

#include "evaluation/evaluation_budget.hpp"
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

struct ExpSeriesSplit final {
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

[[nodiscard]] IntervalSeriesSplit splitExpSeriesInterval(
    const RealInterval& x,
    std::uint64_t begin,
    std::uint64_t end,
    std::size_t precisionBits) {
    if (end == begin + 1) {
        const RealInterval divisor = RealInterval::fromRational(
            Rational{BigInt::fromUnsigned(begin)}, precisionBits);
        const RealInterval ratio = divide(x, divisor, precisionBits);
        return IntervalSeriesSplit{ratio, ratio};
    }

    const std::uint64_t middle = begin + (end - begin) / 2;
    IntervalSeriesSplit left = splitExpSeriesInterval(
        x, begin, middle, precisionBits);
    IntervalSeriesSplit right = splitExpSeriesInterval(
        x, middle, end, precisionBits);
    return IntervalSeriesSplit{
        multiply(left.p, right.p, precisionBits),
        add(left.t, multiply(left.p, right.t, precisionBits), precisionBits)
    };
}

// exp(x) のTaylor項 t_n=x^n/n! は
//   t_n/t_(n-1) = numerator / (denominator*n)
// というhypergeometric recurrenceを持つ。
// [begin,end) をbinary splittingし、T/Qとしてその区間のprefix積和を返す。
[[nodiscard]] ExpSeriesSplit splitExpSeries(
    const BigInt& numerator,
    const BigInt& denominator,
    std::uint64_t begin,
    std::uint64_t end) {
    if (end == begin + 1) {
        const BigInt p = numerator;
        const BigInt q = denominator * BigInt::fromUnsigned(begin);
        return ExpSeriesSplit{p, q, p};
    }

    const std::uint64_t middle = begin + (end - begin) / 2;
    ExpSeriesSplit left = splitExpSeries(numerator, denominator, begin, middle);
    ExpSeriesSplit right = splitExpSeries(numerator, denominator, middle, end);

    ExpSeriesSplit result;
    result.p = left.p * right.p;
    result.q = left.q * right.q;
    result.t = left.t * right.q + left.p * right.t;
    return result;
}

[[nodiscard]] std::size_t requiredExpPower(
    std::size_t precisionBits,
    std::size_t reductionBits) {
    // |x|<=2^-r では、N次まで採った後のtailは
    //   R_N <= 2*t_(N+1) <= 2 / (2^(r(N+1)) (N+1)!)
    // と抑えられる。factorialのbit長だけでこの上界を保証し、
    // 項数決定のために浮動小数点logや経験的停止条件を使わない。
    const std::size_t targetBits = checkedPrecisionAdd(
        precisionBits, 28, "Certified exp precision is too large");

    BigInt factorial{1};
    for (std::uint64_t m = 1;; ++m) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        factorial *= BigInt::fromUnsigned(m);
        const std::size_t factorialBits = factorial.bitLength();
        const std::uint64_t n = m - 1;
        if (n <= static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            const std::size_t nBits = static_cast<std::size_t>(n);
            if (reductionBits != 0
                && m <= std::numeric_limits<std::size_t>::max() / reductionBits) {
                const std::size_t reducedBits = static_cast<std::size_t>(m) * reductionBits;
                if (factorialBits != 0
                    && factorialBits <= std::numeric_limits<std::size_t>::max() - reducedBits
                    && factorialBits + reducedBits >= targetBits + 2)
                    return nBits;
            }
        }

        if (m == std::numeric_limits<std::uint64_t>::max())
            throw std::overflow_error("Certified exp series index is too large");
    }
}


// 0 <= x <= 1/2 の一点をbinary-splitting Taylor級数で囲う。
[[nodiscard]] PointExpResult encloseExpSmallNonNegative(
    const Rational& x,
    std::size_t precisionBits,
    std::size_t reductionBits) {
    if (x < rational(0) || x > rational(1, 2))
        throw std::invalid_argument("Reduced exponential argument must be in [0, 1/2]");

    const RealInterval one = exactIntegerInterval(1, precisionBits);
    if (x.isZero())
        return PointExpResult{one, 1, 0};

    const std::size_t maxPower = requiredExpPower(precisionBits, reductionBits);
    if (maxPower == 0)
        return PointExpResult{one, 1, 0};
    if (maxPower >= std::numeric_limits<std::uint64_t>::max())
        throw std::overflow_error("Certified exp series index is too large");

    const std::uint64_t end = static_cast<std::uint64_t>(maxPower) + 1;

    if (preferExactSeries(x)) {
        const ExpSeriesSplit split = splitExpSeries(
            x.numerator(), x.denominator(), 1, end);

        // T/Q = x/1! + ... + x^N/N!、P/Q = x^N/N!。
        const Rational partial = rational(1) + Rational{split.t, split.q};
        const Rational lastTerm{split.p, split.q};
        const Rational nextTerm = lastTerm * x
            / Rational{BigInt::fromUnsigned(end)};
        const Rational ratio = x / Rational{BigInt::fromUnsigned(end + 1)};
        const Rational tailBound = nextTerm / (rational(1) - ratio);
        return PointExpResult{
            RealInterval::fromRationalBounds(
                partial, partial + tailBound, precisionBits),
            maxPower + 1,
            0
        };
    }

    // 巨大な分子・分母を持つRationalをexact binary splittingすると、
    // q^Nの中間整数が要求precisionを大幅に超えて膨らむ。
    // この場合だけ固定precisionのinterval binary splittingへ切り替え、
    // 証明付き外向き丸めを保ちながら中間整数の無制限成長を避ける。
    const RealInterval xInterval = RealInterval::fromRational(x, precisionBits);
    const IntervalSeriesSplit split = splitExpSeriesInterval(
        xInterval, 1, end, precisionBits);
    const RealInterval partial = add(
        exactIntegerInterval(1, precisionBits), split.t, precisionBits);
    const RealInterval nextIndex = RealInterval::fromRational(
        Rational{BigInt::fromUnsigned(end)}, precisionBits);
    const RealInterval nextTerm = divide(
        multiply(split.p, xInterval, precisionBits), nextIndex, precisionBits);
    const RealInterval ratio = divide(
        xInterval,
        RealInterval::fromRational(
            Rational{BigInt::fromUnsigned(end + 1)}, precisionBits),
        precisionBits);
    const RealInterval tail = divide(
        nextTerm,
        subtract(exactIntegerInterval(1, precisionBits), ratio, precisionBits),
        precisionBits);
    const RealInterval nonNegativeTail = RealInterval::fromRationalBounds(
        rational(0), tail.upper().toRational(), precisionBits);
    return PointExpResult{
        add(partial, nonNegativeTail, precisionBits),
        maxPower + 1,
        0
    };
}

[[nodiscard]] PointExpResult encloseExpPoint(
    Rational x,
    std::size_t precisionBits) {
    const bool negative = x < rational(0);
    if (negative)
        x = -x;

    // exp(x) = exp(x / 2^k)^(2^k)。
    // 小さいRationalは1/2までの縮約でexact binary splittingが最速だが、
    // 巨大な分子・分母を持つ値は固定precision interval演算へ落ちるため、
    // さらに2^-24まで縮約して級数項数を減らした方が速い。
    const bool exactSeries = preferExactSeries(x);
    const std::size_t reductionBits = exactSeries ? 1 : 24;
    BigInt reductionDenominator{1};
    reductionDenominator <<= reductionBits;
    const Rational reductionTarget{BigInt{1}, std::move(reductionDenominator)};

    std::size_t squarings = 0;
    while (x > reductionTarget) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        x /= rational(2);
        ++squarings;
    }

    PointExpResult result = encloseExpSmallNonNegative(
        x, precisionBits, reductionBits);
    result.squarings = squarings;
    for (std::size_t i = 0; i < squarings; ++i) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        result.interval = multiply(result.interval, result.interval, precisionBits);
    }

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

    // expは実軸上で厳密に単調増加するため、exp([a,b])=[exp(a),exp(b)]。
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
