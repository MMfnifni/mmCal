// 実代数数のSturm分離
#include "algebraic_number.hpp"

#include "numeric/big_int.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace mmcal::symbolic {
namespace {

using numeric::BigInt;
using numeric::Rational;

constexpr std::size_t maximumAlgebraicDegree = 64;
constexpr std::size_t maximumIsolationSplits = 1'000'000;

using Polynomial = std::vector<Rational>;

void normalize(Polynomial& polynomial) {
    while (!polynomial.empty() && polynomial.back().isZero())
        polynomial.pop_back();
}

[[nodiscard]] Polynomial normalized(std::span<const Rational> coefficients) {
    Polynomial result(coefficients.begin(), coefficients.end());
    normalize(result);
    if (result.size() > 1) {
        const Rational leading = result.back();
        for (Rational& coefficient : result)
            coefficient /= leading;
    }
    return result;
}

[[nodiscard]] Rational absolute(Rational value) {
    return value < Rational{} ? -value : value;
}

[[nodiscard]] Rational evaluate(const Polynomial& polynomial, const Rational& x) {
    Rational result;
    for (auto iterator = polynomial.rbegin(); iterator != polynomial.rend(); ++iterator)
        result = result * x + *iterator;
    return result;
}

[[nodiscard]] Polynomial derivative(const Polynomial& polynomial) {
    if (polynomial.size() <= 1)
        return {};
    Polynomial result(polynomial.size() - 1);
    for (std::size_t exponent = 1; exponent < polynomial.size(); ++exponent)
        result[exponent - 1] = polynomial[exponent]
            * Rational{BigInt::fromUnsigned(exponent)};
    normalize(result);
    return result;
}

[[nodiscard]] Polynomial remainder(Polynomial dividend, const Polynomial& divisor) {
    if (divisor.empty())
        throw std::logic_error("Polynomial remainder divisor is zero");

    normalize(dividend);
    const std::size_t divisorDegree = divisor.size() - 1;
    const Rational divisorLeading = divisor.back();
    while (!dividend.empty() && dividend.size() - 1 >= divisorDegree) {
        const std::size_t shift = dividend.size() - 1 - divisorDegree;
        const Rational factor = dividend.back() / divisorLeading;
        for (std::size_t i = 0; i <= divisorDegree; ++i)
            dividend[i + shift] -= factor * divisor[i];
        normalize(dividend);
    }
    return dividend;
}

[[nodiscard]] Polynomial polynomialGcd(Polynomial lhs, Polynomial rhs) {
    normalize(lhs);
    normalize(rhs);
    while (!rhs.empty()) {
        Polynomial next = remainder(std::move(lhs), rhs);
        lhs = std::move(rhs);
        rhs = std::move(next);
    }
    if (lhs.empty())
        return lhs;
    const Rational leading = lhs.back();
    for (Rational& coefficient : lhs)
        coefficient /= leading;
    return lhs;
}

[[nodiscard]] Polynomial divideExact(Polynomial dividend, const Polynomial& divisor) {
    if (divisor.empty())
        throw std::logic_error("Polynomial exact divisor is zero");
    normalize(dividend);
    if (dividend.size() < divisor.size())
        throw std::logic_error("Polynomial exact division is not divisible");

    const std::size_t divisorDegree = divisor.size() - 1;
    const Rational divisorLeading = divisor.back();
    Polynomial quotient(dividend.size() - divisorDegree);
    while (!dividend.empty() && dividend.size() - 1 >= divisorDegree) {
        const std::size_t shift = dividend.size() - 1 - divisorDegree;
        const Rational factor = dividend.back() / divisorLeading;
        quotient[shift] += factor;
        for (std::size_t i = 0; i <= divisorDegree; ++i)
            dividend[i + shift] -= factor * divisor[i];
        normalize(dividend);
    }
    if (!dividend.empty())
        throw std::logic_error("Polynomial exact division left a remainder");
    normalize(quotient);
    return quotient;
}

[[nodiscard]] Polynomial canonicalPolynomial(std::span<const Rational> coefficients) {
    Polynomial result = normalized(coefficients);
    if (result.size() <= 2)
        return result;

    Polynomial d = derivative(result);
    Polynomial common = polynomialGcd(result, std::move(d));
    if (common.size() > 1)
        result = divideExact(std::move(result), common);

    if (result.size() > 1) {
        const Rational leading = result.back();
        for (Rational& coefficient : result)
            coefficient /= leading;
    }
    return result;
}

[[nodiscard]] std::vector<Polynomial> sturmSequence(const Polynomial& polynomial) {
    std::vector<Polynomial> sequence;
    sequence.push_back(polynomial);
    Polynomial d = derivative(polynomial);
    if (d.empty())
        return sequence;
    sequence.push_back(std::move(d));

    while (true) {
        Polynomial next = remainder(sequence[sequence.size() - 2], sequence.back());
        if (next.empty())
            break;
        for (Rational& coefficient : next)
            coefficient = -coefficient;
        normalize(next);
        sequence.push_back(std::move(next));
    }
    return sequence;
}

[[nodiscard]] int sign(const Rational& value) noexcept {
    if (value.isZero())
        return 0;
    return value < Rational{} ? -1 : 1;
}

[[nodiscard]] std::size_t variationsAt(
    const std::vector<Polynomial>& sequence,
    const Rational& x) {
    std::size_t changes = 0;
    int previous = 0;
    for (const Polynomial& polynomial : sequence) {
        const int current = sign(evaluate(polynomial, x));
        if (current == 0)
            continue;
        if (previous != 0 && previous != current)
            ++changes;
        previous = current;
    }
    return changes;
}

[[nodiscard]] std::size_t rootsBetween(
    const std::vector<Polynomial>& sequence,
    const Rational& lower,
    const Rational& upper) {
    const std::size_t left = variationsAt(sequence, lower);
    const std::size_t right = variationsAt(sequence, upper);
    return left >= right ? left - right : 0;
}

[[nodiscard]] Rational cauchyBound(const Polynomial& polynomial) {
    const Rational leadingMagnitude = absolute(polynomial.back());
    Rational maximum;
    for (std::size_t i = 0; i + 1 < polynomial.size(); ++i) {
        const Rational ratio = absolute(polynomial[i]) / leadingMagnitude;
        if (maximum < ratio)
            maximum = ratio;
    }
    return maximum + Rational{BigInt{1}};
}

[[nodiscard]] Rational nonRootSplit(
    const Polynomial& polynomial,
    const Rational& lower,
    const Rational& upper) {
    Rational split = (lower + upper) / Rational{BigInt{2}};
    if (!evaluate(polynomial, split).isZero())
        return split;

    // midpointが偶然有理根でも，分割点を少しずらせばSturm端点を根にせずに済む。
    // 区間内の根は有限個なので1/3,2/3,1/4,...を順に試せば必ず見つかる。
    for (std::uint64_t denominator = 3; denominator < 1024; ++denominator) {
        for (std::uint64_t numerator = 1; numerator < denominator; ++numerator) {
            Rational candidate = lower
                + (upper - lower) * Rational{
                    BigInt::fromUnsigned(numerator), BigInt::fromUnsigned(denominator)};
            if (!evaluate(polynomial, candidate).isZero())
                return candidate;
        }
    }
    throw std::runtime_error("Could not choose a non-root Sturm split point");
}

struct PendingInterval final {
    Rational lower;
    Rational upper;
    std::size_t rootCount = 0;
};

[[nodiscard]] std::optional<std::vector<RationalRootInterval>> isolateIntervals(
    const Polynomial& polynomial) {
    if (polynomial.size() <= 1)
        return std::vector<RationalRootInterval>{};
    if (polynomial.size() - 1 > maximumAlgebraicDegree)
        return std::nullopt;

    const auto sturm = sturmSequence(polynomial);
    Rational bound = cauchyBound(polynomial);
    Rational lower = -bound;
    Rational upper = bound;
    while (evaluate(polynomial, lower).isZero() || evaluate(polynomial, upper).isZero()) {
        bound += Rational{BigInt{1}};
        lower = -bound;
        upper = bound;
    }

    const std::size_t total = rootsBetween(sturm, lower, upper);
    if (total == 0)
        return std::vector<RationalRootInterval>{};

    std::vector<RationalRootInterval> result;
    result.reserve(total);
    std::vector<PendingInterval> stack;
    stack.push_back(PendingInterval{lower, upper, total});
    std::size_t splits = 0;

    while (!stack.empty()) {
        PendingInterval current = std::move(stack.back());
        stack.pop_back();
        if (current.rootCount == 0)
            continue;
        if (current.rootCount == 1) {
            result.push_back(RationalRootInterval{
                std::move(current.lower), std::move(current.upper)});
            continue;
        }
        if (++splits > maximumIsolationSplits)
            return std::nullopt;

        const Rational split = nonRootSplit(polynomial, current.lower, current.upper);
        const std::size_t leftCount = rootsBetween(sturm, current.lower, split);
        const std::size_t rightCount = current.rootCount - leftCount;
        // stackはLIFOなので右を先に積み，最終resultを昇順にする。
        if (rightCount != 0)
            stack.push_back(PendingInterval{split, current.upper, rightCount});
        if (leftCount != 0)
            stack.push_back(PendingInterval{current.lower, split, leftCount});
    }

    return result;
}

[[nodiscard]] Rational powerOfTwo(std::size_t exponent) {
    BigInt value{1};
    value <<= exponent;
    return Rational{std::move(value)};
}

[[nodiscard]] bool narrowEnough(
    const RationalRootInterval& interval,
    std::size_t precisionBits) {
    if (interval.isPoint())
        return true;
    const Rational width = interval.upper - interval.lower;
    const Rational absLower = absolute(interval.lower);
    const Rational absUpper = absolute(interval.upper);
    const Rational minimumMagnitude = absLower < absUpper ? absLower : absUpper;
    if (interval.lower <= Rational{} && interval.upper >= Rational{})
        return false;
    return width * powerOfTwo(precisionBits + 8) <= minimumMagnitude;
}

} // namespace

RealAlgebraicNumber::RealAlgebraicNumber(
    std::vector<Rational> polynomial,
    std::size_t rootIndex,
    RationalRootInterval interval)
    : polynomial_(std::move(polynomial)),
      rootIndex_(rootIndex),
      interval_(std::move(interval)) {}

std::optional<RealAlgebraicNumber> RealAlgebraicNumber::create(
    std::span<const Rational> polynomial,
    std::size_t rootIndex) {
    if (rootIndex == 0)
        return std::nullopt;
    Polynomial normalizedPolynomial = normalized(polynomial);
    if (normalizedPolynomial.size() <= 1)
        return std::nullopt;
    if (normalizedPolynomial.size() - 1 > maximumAlgebraicDegree)
        return std::nullopt;
    normalizedPolynomial = canonicalPolynomial(normalizedPolynomial);
    const auto intervals = isolateIntervals(normalizedPolynomial);
    if (!intervals || rootIndex > intervals->size())
        return std::nullopt;
    return RealAlgebraicNumber{
        std::move(normalizedPolynomial), rootIndex, (*intervals)[rootIndex - 1]};
}

std::optional<std::vector<RealAlgebraicNumber>> RealAlgebraicNumber::isolateAll(
    std::span<const Rational> polynomial) {
    Polynomial normalizedPolynomial = normalized(polynomial);
    if (normalizedPolynomial.size() <= 1)
        return std::vector<RealAlgebraicNumber>{};
    if (normalizedPolynomial.size() - 1 > maximumAlgebraicDegree)
        return std::nullopt;
    normalizedPolynomial = canonicalPolynomial(normalizedPolynomial);
    const auto intervals = isolateIntervals(normalizedPolynomial);
    if (!intervals)
        return std::nullopt;

    std::vector<RealAlgebraicNumber> roots;
    roots.reserve(intervals->size());
    for (std::size_t i = 0; i < intervals->size(); ++i)
        roots.push_back(RealAlgebraicNumber{normalizedPolynomial, i + 1, (*intervals)[i]});
    return roots;
}

std::span<const Rational> RealAlgebraicNumber::polynomial() const noexcept {
    return polynomial_;
}

std::size_t RealAlgebraicNumber::degree() const noexcept {
    return polynomial_.size() - 1;
}

std::size_t RealAlgebraicNumber::rootIndex() const noexcept { return rootIndex_; }

const RationalRootInterval& RealAlgebraicNumber::isolatingInterval() const noexcept {
    return interval_;
}

RationalRootInterval RealAlgebraicNumber::refined(std::size_t precisionBits) const {
    RationalRootInterval result = interval_;
    if (result.isPoint())
        return result;
    if (result.lower <= Rational{} && result.upper >= Rational{}
        && evaluate(polynomial_, Rational{}).isZero())
        return RationalRootInterval{Rational{}, Rational{}};

    const auto sturm = sturmSequence(polynomial_);
    std::size_t iterations = 0;
    const std::size_t maximumIterations = precisionBits > (std::numeric_limits<std::size_t>::max() - 4096) / 4
        ? std::numeric_limits<std::size_t>::max()
        : precisionBits * 4 + 4096;

    while (!narrowEnough(result, precisionBits)) {
        if (++iterations > maximumIterations)
            throw std::runtime_error("Algebraic root refinement did not converge");

        const Rational midpoint = (result.lower + result.upper) / Rational{BigInt{2}};
        if (evaluate(polynomial_, midpoint).isZero())
            return RationalRootInterval{midpoint, midpoint};

        const std::size_t leftRoots = rootsBetween(sturm, result.lower, midpoint);
        if (leftRoots == 1)
            result.upper = midpoint;
        else
            result.lower = midpoint;
    }
    return result;
}

} // namespace mmcal::symbolic
