// 特殊函数の保証付き評価
#include "certified_special_functions.hpp"

#include "certification_error.hpp"
#include "certified_constants.hpp"
#include "certified_exponential.hpp"
#include "certified_logarithm.hpp"
#include "certified_sqrt.hpp"
#include "certified_trigonometry.hpp"
#include "interval_math.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <stdexcept>
#include <vector>

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

[[nodiscard]] Rational binaryThreshold(std::size_t bits) {
    BigInt denominator{1};
    denominator <<= bits;
    return Rational{BigInt{1}, std::move(denominator)};
}

[[nodiscard]] Rational absRational(Rational value) {
    return value.numerator().isNegative() ? -value : value;
}

[[nodiscard]] RealInterval exactInterval(
    const Rational& value,
    std::size_t bits) {
    return RealInterval::fromRational(value, bits);
}

[[nodiscard]] RealInterval exactInterval(std::int64_t value, std::size_t bits) {
    return exactInterval(rational(value), bits);
}

// Akiyama-Tanigawa法でBernoulli数をexact Rationalとして一度だけ生成する。GammaのStirling剰余評価では偶数添字だけを使う。
[[nodiscard]] const std::vector<Rational>& bernoulliNumbers() {
    static const std::vector<Rational> values = [] {
        constexpr std::size_t maximum = 128;
        std::vector<Rational> a(maximum + 1);
        std::vector<Rational> b(maximum + 1);
        for (std::size_t m = 0; m <= maximum; ++m) {
            a[m] = Rational{BigInt{1}, BigInt::parse(std::to_string(m + 1))};
            for (std::size_t j = m; j >= 1; --j) {
                a[j - 1] = Rational{BigInt::parse(std::to_string(j))} * (a[j - 1] - a[j]);
                if (j == 1)
                    break;
            }
            b[m] = a[0];
        }
        return b;
    }();
    return values;
}

[[nodiscard]] Rational stirlingCoefficient(std::size_t k) {
    const auto& b = bernoulliNumbers();
    const std::size_t n = 2 * k;
    const BigInt denominator = BigInt::parse(std::to_string(n))
        * BigInt::parse(std::to_string(n - 1));
    return b[n] / Rational{denominator};
}

struct StirlingPlan final {
    std::size_t shift = 0;
    std::size_t omittedK = 0;
    Rational remainderBound;
};

[[nodiscard]] StirlingPlan chooseStirlingPlan(
    const Rational& inputLower,
    std::size_t precisionBits) {
    if (inputLower <= rational(0))
        throw std::domain_error("LogGamma requires a positive interval in the Stirling backend");

    const Rational threshold = binaryThreshold(checkedAdd(
        precisionBits, 20, "Gamma precision is too large"));
    constexpr std::size_t maximumK = 64;

    for (std::size_t shift = 0; shift <= 1'000'000; shift += 8) {
        const Rational x = inputLower + Rational{BigInt::parse(std::to_string(shift))};
        if (x < rational(4))
            continue;

        const Rational inverseSquare = rational(1) / (x * x);
        Rational inverseOdd = rational(1) / x;
        for (std::size_t k = 1; k <= maximumK; ++k) {
            const Rational bound = absRational(stirlingCoefficient(k)) * inverseOdd;
            if (bound <= threshold)
                return StirlingPlan{shift, k, bound};
            inverseOdd *= inverseSquare;
        }
    }
    throw std::overflow_error("Gamma precision requires an excessive recurrence shift");
}

[[nodiscard]] RealInterval encloseLogGammaPositive(
    const RealInterval& input,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (input.lower() <= zero)
        throw std::domain_error("LogGamma positive backend requires x > 0");

    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Gamma working precision is too large");
    const StirlingPlan plan = chooseStirlingPlan(
        input.lower().toRational(), workBits);

    const RealInterval shift = exactInterval(
        Rational{BigInt::parse(std::to_string(plan.shift))}, workBits);
    const RealInterval x = add(input.roundedOutward(workBits), shift, workBits);
    const RealInterval logX = encloseLogPositive(x, workBits).interval;
    const RealInterval xMinusHalf = subtract(x, exactInterval(rational(1, 2), workBits), workBits);

    RealInterval result = subtract(
        multiply(xMinusHalf, logX, workBits), x, workBits);

    // 1/2 log(2 Pi)
    const RealInterval twoPi = multiply(
        enclosePi(workBits).interval, exactInterval(2, workBits), workBits);
    const RealInterval halfLogTwoPi = multiply(
        encloseLogPositive(twoPi, workBits).interval,
        exactInterval(rational(1, 2), workBits), workBits);
    result = add(result, halfLogTwoPi, workBits);

    const RealInterval inverseX = divide(exactInterval(1, workBits), x, workBits);
    const RealInterval inverseSquare = multiply(inverseX, inverseX, workBits);
    RealInterval inverseOdd = inverseX;
    for (std::size_t k = 1; k < plan.omittedK; ++k) {
        const RealInterval coefficient = exactInterval(stirlingCoefficient(k), workBits);
        result = add(result, multiply(coefficient, inverseOdd, workBits), workBits);
        inverseOdd = multiply(inverseOdd, inverseSquare, workBits);
    }

    // 正実軸上のStirling級数の剰余は最初の省略項と同符号で、絶対値はその項を超えない。入力区間ではlower endpointが最大絶対値を与える。
    const Rational omittedCoefficient = stirlingCoefficient(plan.omittedK);
    const Rational bound = plan.remainderBound;
    const RealInterval remainder = omittedCoefficient.numerator().isNegative()
        ? RealInterval::fromRationalBounds(-bound, rational(0), workBits)
        : RealInterval::fromRationalBounds(rational(0), bound, workBits);
    result = add(result, remainder, workBits);

    if (plan.shift != 0) {
        // Gamma(x+s)=Gamma(x) product_{j=0}^{s-1}(x+j)。各logを個別に取らずproductを区間で作ってから1回だけlogを取る。
        RealInterval product = exactInterval(1, workBits);
        for (std::size_t j = 0; j < plan.shift; ++j) {
            const RealInterval factor = add(
                input.roundedOutward(workBits),
                exactInterval(Rational{BigInt::parse(std::to_string(j))}, workBits),
                workBits);
            product = multiply(product, factor, workBits);
        }
        result = subtract(result, encloseLogPositive(product, workBits).interval, workBits);
    }

    return result.roundedOutward(precisionBits);
}

[[nodiscard]] bool exactNonPositiveIntegerPoint(const RealInterval& input) {
    if (!input.isPoint())
        return false;
    const Rational value = input.lower().toRational();
    return value.isInteger() && value.numerator() <= BigInt{0};
}

[[nodiscard]] RealInterval gammaNegativeByReflection(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (exactNonPositiveIntegerPoint(input))
        throw std::domain_error("gamma is undefined at non-positive integers");

    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "Gamma reflection precision is too large");
    const RealInterval oneMinusX = subtract(
        exactInterval(1, workBits), input.roundedOutward(workBits), workBits);
    const RealInterval gammaComplement = encloseExp(
        encloseLogGammaPositive(oneMinusX, workBits), workBits).interval;
    const RealInterval pi = enclosePi(workBits).interval;
    const RealInterval piX = multiply(pi, input.roundedOutward(workBits), workBits);
    const RealInterval sine = encloseSinRadianInterval(piX, workBits).interval;
    if (sine.containsZero())
        throw PrecisionInsufficient{"Gamma reflection cannot yet exclude a non-positive-integer pole"};
    const RealInterval denominator = multiply(sine, gammaComplement, workBits);
    if (denominator.containsZero())
        throw PrecisionInsufficient{"Gamma reflection denominator cannot yet be proven nonzero"};
    return divide(pi, denominator, workBits).roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval twoOverSqrtPi(std::size_t bits) {
    const RealInterval rootPi = encloseSqrt(enclosePi(bits).interval, bits).interval;
    return divide(exactInterval(2, bits), rootPi, bits);
}


[[nodiscard]] Rational unsignedRational(std::size_t value) {
    if (value > static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max()))
        throw std::overflow_error("Fresnel series index is too large");
    return Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(value))};
}

[[nodiscard]] RealInterval symmetricError(
    const Rational& bound,
    std::size_t precisionBits) {
    return RealInterval::fromRationalBounds(-bound, bound, precisionBits);
}

[[nodiscard]] Rational intervalAbsUpper(
    const RealInterval& interval,
    std::size_t precisionBits) {
    return absoluteInterval(interval, precisionBits).upper().toRational();
}

[[nodiscard]] RealInterval pointFresnelSeriesPositive(
    const Rational& x,
    bool cosineIntegral,
    std::size_t precisionBits) {
    // FresnelのMaclaurin級数は大きいxでは巨大な中間項が相殺する。
    // 旧+48bit固定guardではx=4程度でも高精度時に包含幅が縮まらないため、
    // x^2に比例したguardを追加して相殺分を明示的に吸収する。
    const Rational x2ForGuard = x * x;
    const BigInt guardQuotient = x2ForGuard.numerator() / x2ForGuard.denominator();
    const auto guardMagnitude = numeric::tryToUint64(guardQuotient);
    const std::size_t cancellationGuard = guardMagnitude
        ? static_cast<std::size_t>(std::min<std::uint64_t>(*guardMagnitude, 100'000ULL)) * 4U
        : 400'000U;
    const std::size_t workBits = checkedAdd(
        precisionBits,
        checkedAdd(64, cancellationGuard, "Fresnel cancellation guard is too large"),
        "Fresnel working precision is too large");
    const RealInterval pi = enclosePi(workBits).interval;
    const Rational piUpper = pi.upper().toRational();
    const Rational x2 = x * x;
    const Rational x4 = x2 * x2;
    const Rational commonUpper = piUpper * piUpper * x4 / rational(4);
    const RealInterval common = divide(
        multiply(multiply(pi, pi, workBits), exactInterval(x4, workBits), workBits),
        exactInterval(4, workBits), workBits);

    RealInterval term = cosineIntegral
        ? exactInterval(x, workBits)
        : divide(
            multiply(pi, exactInterval(x * x2, workBits), workBits),
            exactInterval(6, workBits), workBits);
    RealInterval sum = term;
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 12, "Fresnel target precision is too large"));

    constexpr std::size_t maximumTerms = 1'000'000;
    for (std::size_t n = 0; n < maximumTerms; ++n) {
        const std::size_t numeratorIndex = cosineIntegral ? 4 * n + 1 : 4 * n + 3;
        const std::size_t d0 = cosineIntegral ? 2 * n + 1 : 2 * n + 2;
        const std::size_t d1 = cosineIntegral ? 2 * n + 2 : 2 * n + 3;
        const std::size_t d2 = cosineIntegral ? 4 * n + 5 : 4 * n + 7;
        const Rational ratioUpper = commonUpper * unsignedRational(numeratorIndex)
            / (unsignedRational(d0) * unsignedRational(d1) * unsignedRational(d2));

        // 現項より後の比が1未満に入れば以後は単調減少する。
        // 最初の未加算項を等比級数で上から押さえ、Taylor剰余を明示的に区間へ足す。
        if (ratioUpper < rational(1)) {
            const Rational nextBound = intervalAbsUpper(term, workBits) * ratioUpper;
            const Rational tailBound = nextBound / (rational(1) - ratioUpper);
            if (tailBound <= target) {
                sum = add(sum, symmetricError(tailBound, workBits), workBits);
                return sum.roundedOutward(precisionBits);
            }
        }

        const RealInterval ratio = divide(
            multiply(common, exactInterval(unsignedRational(numeratorIndex), workBits), workBits),
            exactInterval(unsignedRational(d0) * unsignedRational(d1) * unsignedRational(d2), workBits),
            workBits);
        term = negate(multiply(term, ratio, workBits));
        sum = add(sum, term, workBits);
    }
    throw std::overflow_error("Fresnel series requires too many terms");
}

struct FresnelPair final {
    RealInterval c;
    RealInterval s;
};

[[nodiscard]] FresnelPair pointFresnelAsymptoticPositive(
    const Rational& x,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 64, "Fresnel asymptotic precision is too large");
    const RealInterval pi = enclosePi(workBits).interval;
    const RealInterval xInterval = exactInterval(x, workBits);
    const RealInterval x2 = exactInterval(x * x, workBits);
    const RealInterval phase = divide(
        multiply(pi, x2, workBits), exactInterval(2, workBits), workBits);
    const RealInterval sine = encloseSinRadianInterval(phase, workBits).interval;
    const RealInterval cosine = encloseCosRadianInterval(phase, workBits).interval;

    RealInterval amplitude = divide(
        exactInterval(1, workBits),
        multiply(pi, xInterval, workBits), workBits);
    RealInterval tailReal = exactInterval(0, workBits);
    RealInterval tailImag = exactInterval(0, workBits);
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 14, "Fresnel asymptotic target precision is too large"));

    constexpr std::size_t maximumTerms = 4096;
    Rational remainderBound;
    for (std::size_t m = 0; m < maximumTerms; ++m) {
        RealInterval realPart = exactInterval(0, workBits);
        RealInterval imagPart = exactInterval(0, workBits);
        switch (m & 3U) {
        case 0: // i A_m
            realPart = negate(multiply(sine, amplitude, workBits));
            imagPart = multiply(cosine, amplitude, workBits);
            break;
        case 1: // +A_m
            realPart = multiply(cosine, amplitude, workBits);
            imagPart = multiply(sine, amplitude, workBits);
            break;
        case 2: // -i A_m
            realPart = multiply(sine, amplitude, workBits);
            imagPart = negate(multiply(cosine, amplitude, workBits));
            break;
        case 3: // -A_m
            realPart = negate(multiply(cosine, amplitude, workBits));
            imagPart = negate(multiply(sine, amplitude, workBits));
            break;
        }
        tailReal = add(tailReal, realPart, workBits);
        tailImag = add(tailImag, imagPart, workBits);

        // m+1項まで展開した部分積分公式の剰余は、最後に加えたA_m以下。
        // x>=4ではA_mが必要精度まで減少する範囲で打ち切るため、発散域へ進まない。
        remainderBound = intervalAbsUpper(amplitude, workBits);
        if (remainderBound <= target) {
            const RealInterval error = symmetricError(remainderBound, workBits);
            const RealInterval half = exactInterval(rational(1, 2), workBits);
            return FresnelPair{
                add(subtract(half, tailReal, workBits), error, workBits).roundedOutward(precisionBits),
                add(subtract(half, tailImag, workBits), error, workBits).roundedOutward(precisionBits)};
        }

        const Rational odd = unsignedRational(2 * m + 1);
        const RealInterval scale = divide(
            exactInterval(odd, workBits),
            multiply(pi, x2, workBits), workBits);
        const RealInterval nextAmplitude = multiply(amplitude, scale, workBits);
        if (intervalAbsUpper(nextAmplitude, workBits) >= remainderBound)
            break;
        amplitude = nextAmplitude;
    }
    throw PrecisionInsufficient{"Fresnel asymptotic expansion did not reach the requested precision"};
}

[[nodiscard]] RealInterval pointFresnel(
    Rational x,
    bool cosineIntegral,
    std::size_t precisionBits) {
    if (x.isZero())
        return exactInterval(0, precisionBits);
    if (x.numerator().isNegative())
        return negate(pointFresnel(-x, cosineIntegral, precisionBits));

    if (x < rational(8))
        return pointFresnelSeriesPositive(x, cosineIntegral, precisionBits);

    // 漸近級数は固定xで任意精度まで収束する級数ではない。
    // 要求精度に届かない場合は正則なMaclaurin級数へ戻し、速度のために保証を捨てない。
    try {
        const FresnelPair pair = pointFresnelAsymptoticPositive(x, precisionBits);
        return cosineIntegral ? pair.c : pair.s;
    }
    catch (const PrecisionInsufficient&) {
        return pointFresnelSeriesPositive(x, cosineIntegral, precisionBits);
    }
}

[[nodiscard]] RealInterval encloseFresnelRealImpl(
    const RealInterval& input,
    bool cosineIntegral,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Fresnel interval precision is too large");
    const Rational lower = input.lower().toRational();
    const Rational upper = input.upper().toRational();
    RealInterval value = pointFresnel(lower, cosineIntegral, workBits);

    // |C'(x)|=|cos(pi x^2/2)|<=1, |S'(x)|<=1。
    // 入力が丸め区間でもlower endpointからの距離だけ膨らませれば真値を必ず包含できる。
    const Rational width = upper - lower;
    if (!width.isZero())
        value = add(value, symmetricError(width, workBits), workBits);
    return value.roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval pointErfSeriesPositive(
    const Rational& x,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "erf working precision is too large");
    const Rational threshold = binaryThreshold(checkedAdd(
        workBits, 12, "erf precision is too large"));
    const Rational x2 = x * x;

    Rational term = x;
    Rational sum = term;
    Rational tailBound;
    std::uint64_t n = 0;
    for (;;) {
        const std::uint64_t next = n + 1;
        const Rational ratio = -x2
            * Rational{BigInt::parse(std::to_string(2 * n + 1))}
            / Rational{
                BigInt::parse(std::to_string(next))
                * BigInt::parse(std::to_string(2 * n + 3))};
        const Rational nextTerm = term * ratio;

        // 次項以降の絶対比は x^2/(n+2) より小さい。
        const Rational q = x2 / Rational{BigInt::parse(std::to_string(n + 2))};
        if (q < rational(1)) {
            tailBound = absRational(nextTerm) / (rational(1) - q);
            if (tailBound <= threshold)
                break;
        }

        sum += nextTerm;
        term = nextTerm;
        n = next;
        if (n > 1'000'000)
            throw std::overflow_error("erf series did not converge within the iteration limit");
    }

    const RealInterval sumWithTail = RealInterval::fromRationalBounds(
        sum - tailBound, sum + tailBound, workBits);
    return multiply(sumWithTail, twoOverSqrtPi(workBits), workBits)
        .roundedOutward(precisionBits);
}

[[nodiscard]] std::optional<RealInterval> pointErfcAsymptoticPositive(
    const Rational& x,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "erfc working precision is too large");
    const Rational threshold = binaryThreshold(checkedAdd(
        workBits, 12, "erfc precision is too large"));
    const Rational x2 = x * x;

    const RealInterval xInterval = exactInterval(x, workBits);
    const RealInterval exponential = encloseExp(
        exactInterval(-x2, workBits), workBits).interval;
    const RealInterval rootPi = encloseSqrt(enclosePi(workBits).interval, workBits).interval;
    const RealInterval prefactor = divide(
        exponential, multiply(xInterval, rootPi, workBits), workBits);

    Rational term{BigInt{1}};
    Rational sum = term;
    Rational previousMagnitude = absRational(term);
    for (std::uint64_t n = 0; n < 1'000'000; ++n) {
        const Rational nextTerm = -term
            * Rational{BigInt::parse(std::to_string(2 * n + 1))}
            / (rational(2) * x2);
        const Rational nextMagnitude = absRational(nextTerm);
        const RealInterval error = multiply(
            prefactor, exactInterval(nextMagnitude, workBits), workBits);
        if (error.upper().toRational() <= threshold) {
            const RealInterval series = RealInterval::fromRationalBounds(
                sum - nextMagnitude, sum + nextMagnitude, workBits);
            return multiply(prefactor, series, workBits).roundedOutward(precisionBits);
        }

        // 漸近級数は最小項を越えると発散する。要求精度へ届く前に項が増加へ転じた場合はMaclaurin側へfallbackする。
        if (nextMagnitude >= previousMagnitude)
            return std::nullopt;

        sum += nextTerm;
        term = nextTerm;
        previousMagnitude = nextMagnitude;
    }
    return std::nullopt;
}

[[nodiscard]] RealInterval pointErf(
    const Rational& x,
    std::size_t precisionBits) {
    if (x.isZero())
        return exactInterval(0, precisionBits);
    if (x.numerator().isNegative())
        return negate(pointErf(-x, precisionBits));

    if (x < rational(4))
        return pointErfSeriesPositive(x, precisionBits);

    if (const auto erfc = pointErfcAsymptoticPositive(x, precisionBits))
        return subtract(exactInterval(1, precisionBits), *erfc, precisionBits);
    return pointErfSeriesPositive(x, precisionBits);
}

[[nodiscard]] RealInterval pointErfc(
    const Rational& x,
    std::size_t precisionBits) {
    if (x.isZero())
        return exactInterval(1, precisionBits);
    if (x.numerator().isNegative())
        return subtract(exactInterval(2, precisionBits), pointErfc(-x, precisionBits), precisionBits);
    if (x < rational(4))
        return subtract(exactInterval(1, precisionBits), pointErfSeriesPositive(x, precisionBits), precisionBits);
    if (const auto asymptotic = pointErfcAsymptoticPositive(x, precisionBits))
        return *asymptotic;
    return subtract(exactInterval(1, precisionBits), pointErfSeriesPositive(x, precisionBits), precisionBits);
}

} // namespace


[[nodiscard]] std::optional<std::uint64_t> ceilAbsToUint64(const Rational& value) {
    const BigInt numerator = value.numerator().abs();
    const BigInt& denominator = value.denominator();
    if (numerator.isZero())
        return 0;
    const BigInt quotient = (numerator + denominator - BigInt{1}) / denominator;
    return numeric::tryToUint64(quotient);
}

[[nodiscard]] bool nonPositiveInteger(const Rational& value) {
    return value.isInteger() && !value.numerator().isPositive();
}

[[nodiscard]] RealInterval pointHypergeometric1F1(
    const Rational& a,
    const Rational& b,
    const Rational& z,
    std::size_t precisionBits) {
    if (nonPositiveInteger(b))
        throw std::domain_error("hypergeometric1F1 has a pole at a non-positive integer b");
    if (z.isZero())
        return exactInterval(1, precisionBits);

    const auto absA = ceilAbsToUint64(a);
    const auto absB = ceilAbsToUint64(b);
    const auto absZ = ceilAbsToUint64(z);
    if (!absA || !absB || !absZ)
        throw PrecisionInsufficient{"hypergeometric1F1 argument is too large for the series backend"};

    // j>=2|a|,2|b|,6|z| なら
    // |t_{j+1}/t_j| = |z||a+j|/(|b+j|(j+1)) <= 1/2。
    // 以後のtailは次項の2倍で厳密に上から押さえられる。
    const std::uint64_t ratioStart = std::max({
        2U * *absA + 2U,
        2U * *absB + 2U,
        6U * *absZ + 2U});
    constexpr std::uint64_t maximumTerms = 200000;
    if (ratioStart > maximumTerms)
        throw PrecisionInsufficient{"hypergeometric1F1 requires too many series terms"};

    Rational term{BigInt{1}};
    Rational sum{BigInt{1}};
    const Rational tolerance = binaryThreshold(checkedAdd(
        precisionBits, 16, "hypergeometric1F1 precision is too large"));

    for (std::uint64_t n = 0; n < maximumTerms; ++n) {
        const Rational numeratorFactor = a + Rational{BigInt::fromUnsigned(n)};
        const Rational denominatorFactor = b + Rational{BigInt::fromUnsigned(n)};
        if (denominatorFactor.isZero())
            throw std::domain_error("hypergeometric1F1 denominator parameter reaches a pole");

        const Rational next = term * numeratorFactor * z
            / (denominatorFactor * Rational{BigInt::fromUnsigned(n + 1)});

        if (n >= ratioStart) {
            const Rational tailBound = rational(2) * absRational(next);
            if (tailBound <= tolerance)
                return RealInterval::fromRationalBounds(
                    sum - tailBound, sum + tailBound, precisionBits);
        }

        term = next;
        sum += term;
        if (term.isZero())
            return exactInterval(sum, precisionBits);
    }

    throw PrecisionInsufficient{"hypergeometric1F1 series did not converge within the term limit"};
}

RealInterval encloseGammaReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Gamma precision must be at least one bit");
    if (exactNonPositiveIntegerPoint(input))
        throw std::domain_error("gamma is undefined at non-positive integers");

    const BigFloat zero;
    if (input.lower() > zero)
        return encloseExp(encloseLogGammaPositive(input, precisionBits), precisionBits).interval;
    if (input.upper() < zero)
        return gammaNegativeByReflection(input, precisionBits);

    throw PrecisionInsufficient{"Gamma interval straddles zero or a possible pole"};
}

RealInterval encloseLogGammaReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (exactNonPositiveIntegerPoint(input))
        throw std::domain_error("lgamma is undefined at non-positive integers");

    const BigFloat zero;
    if (input.lower() > zero)
        return encloseLogGammaPositive(input, precisionBits);

    const RealInterval gamma = encloseGammaReal(input, checkedAdd(
        precisionBits, 24, "lgamma working precision is too large"));
    const RealInterval magnitude = absoluteInterval(gamma, checkedAdd(
        precisionBits, 16, "lgamma absolute-value precision is too large"));
    if (magnitude.containsZero())
        throw PrecisionInsufficient{"lgamma could not yet prove Gamma away from zero"};
    return encloseLogPositive(magnitude, precisionBits).interval;
}

RealInterval encloseErfReal(
    const RealInterval& input,
    std::size_t precisionBits) { // erfは実軸上で単調増加。
    const RealInterval lower = pointErf(input.lower().toRational(), precisionBits);
    const RealInterval upper = pointErf(input.upper().toRational(), precisionBits);
    return RealInterval{lower.lower(), upper.upper()};
}

RealInterval encloseErfcReal(
    const RealInterval& input,
    std::size_t precisionBits) { // erfcは実軸上で単調減少。
    const RealInterval lower = pointErfc(input.upper().toRational(), precisionBits);
    const RealInterval upper = pointErfc(input.lower().toRational(), precisionBits);
    return RealInterval{lower.lower(), upper.upper()};
}

RealInterval encloseBetaLogPositive(
    const RealInterval& a,
    const RealInterval& b,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (a.lower() <= zero || b.lower() <= zero)
        throw std::domain_error("betaln certified evaluation currently requires a > 0 and b > 0");
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Beta working precision is too large");
    const RealInterval sum = add(a, b, workBits);
    return subtract(
        add(encloseLogGammaPositive(a, workBits), encloseLogGammaPositive(b, workBits), workBits),
        encloseLogGammaPositive(sum, workBits), workBits)
        .roundedOutward(precisionBits);
}


RealInterval encloseFresnelCReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Fresnel precision must be at least one bit");
    return encloseFresnelRealImpl(input, true, precisionBits);
}

RealInterval encloseFresnelSReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Fresnel precision must be at least one bit");
    return encloseFresnelRealImpl(input, false, precisionBits);
}


RealInterval encloseHypergeometric1F1Real(
    const Rational& a,
    const Rational& b,
    const Rational& z,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("hypergeometric1F1 precision must be at least one bit");
    return pointHypergeometric1F1(a, b, z, precisionBits);
}

RealInterval encloseBetaPositive(
    const RealInterval& a,
    const RealInterval& b,
    std::size_t precisionBits) {
    return encloseExp(encloseBetaLogPositive(a, b, precisionBits), precisionBits).interval;
}

} // namespace mmcal::approximation
