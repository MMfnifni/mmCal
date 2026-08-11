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

RealInterval encloseBetaPositive(
    const RealInterval& a,
    const RealInterval& b,
    std::size_t precisionBits) {
    return encloseExp(encloseBetaLogPositive(a, b, precisionBits), precisionBits).interval;
}

} // namespace mmcal::approximation
