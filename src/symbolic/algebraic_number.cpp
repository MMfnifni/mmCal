// 実代数数のSturm分離
#include "algebraic_number.hpp"
#include "number_field.hpp"
#include "rational_linear_basis.hpp"

#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/big_float.hpp"
#include "approximation/certified_constants.hpp"
#include "approximation/certified_trigonometry.hpp"
#include "approximation/real_interval.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace mmcal::symbolic {
namespace {

using numeric::BigInt;
using numeric::Rational;

constexpr std::size_t maximumAlgebraicDegree = 64;
constexpr std::size_t maximumIsolationSplits = 1'000'000;
constexpr std::size_t maximumComparisonRefinementBits = 4096;

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


// minimal polynomial縮約はresultant次数budgetと同程度の小～中次数へ限定する。
// 高次数を無理にfactorしてRoot生成自体を重くしない。
constexpr std::size_t maximumMinimalPolynomialDegree = 16;
constexpr std::size_t maximumKroneckerDivisors = 4096;
constexpr std::size_t maximumKroneckerCombinations = 200'000;

using ModPolynomial = std::vector<std::uint32_t>;

[[nodiscard]] std::uint32_t modMultiply(
    std::uint32_t lhs,
    std::uint32_t rhs,
    std::uint32_t modulus) noexcept {
    return static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(lhs) * rhs) % modulus);
}

[[nodiscard]] std::uint32_t modPower(
    std::uint32_t base,
    std::uint64_t exponent,
    std::uint32_t modulus) noexcept {
    std::uint32_t result = 1 % modulus;
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = modMultiply(result, base, modulus);
        exponent >>= 1U;
        if (exponent != 0)
            base = modMultiply(base, base, modulus);
    }
    return result;
}

void normalizeMod(ModPolynomial& polynomial) {
    while (!polynomial.empty() && polynomial.back() == 0)
        polynomial.pop_back();
}

[[nodiscard]] std::optional<std::uint32_t> rationalModulo(
    const Rational& value,
    std::uint32_t prime) {
    const BigInt p = BigInt::fromUnsigned(prime);
    BigInt numerator = value.numerator() % p;
    BigInt denominator = value.denominator() % p;
    if (numerator.isNegative())
        numerator += p;
    if (denominator.isNegative())
        denominator += p;
    const auto n = numeric::tryToUint64(numerator);
    const auto d = numeric::tryToUint64(denominator);
    if (!n || !d || *d == 0)
        return std::nullopt;
    const std::uint32_t inverse = modPower(
        static_cast<std::uint32_t>(*d), prime - 2, prime);
    return modMultiply(static_cast<std::uint32_t>(*n), inverse, prime);
}

[[nodiscard]] std::optional<ModPolynomial> polynomialModulo(
    const Polynomial& polynomial,
    std::uint32_t prime) {
    ModPolynomial result;
    result.reserve(polynomial.size());
    for (const Rational& coefficient : polynomial) {
        const auto value = rationalModulo(coefficient, prime);
        if (!value)
            return std::nullopt;
        result.push_back(*value);
    }
    normalizeMod(result);
    if (result.size() != polynomial.size())
        return std::nullopt;
    return result;
}

[[nodiscard]] ModPolynomial modRemainder(
    ModPolynomial dividend,
    const ModPolynomial& divisor,
    std::uint32_t prime) {
    normalizeMod(dividend);
    if (divisor.empty())
        throw std::logic_error("Finite-field polynomial divisor is zero");
    const std::size_t divisorDegree = divisor.size() - 1;
    const std::uint32_t inverseLeading = modPower(divisor.back(), prime - 2, prime);
    while (!dividend.empty() && dividend.size() - 1 >= divisorDegree) {
        const std::size_t shift = dividend.size() - 1 - divisorDegree;
        const std::uint32_t factor = modMultiply(dividend.back(), inverseLeading, prime);
        for (std::size_t i = 0; i <= divisorDegree; ++i) {
            const std::uint32_t amount = modMultiply(factor, divisor[i], prime);
            dividend[i + shift] = dividend[i + shift] >= amount
                ? dividend[i + shift] - amount
                : static_cast<std::uint32_t>(dividend[i + shift] + prime - amount);
        }
        normalizeMod(dividend);
    }
    return dividend;
}

[[nodiscard]] ModPolynomial modGcd(
    ModPolynomial lhs,
    ModPolynomial rhs,
    std::uint32_t prime) {
    normalizeMod(lhs);
    normalizeMod(rhs);
    while (!rhs.empty()) {
        ModPolynomial next = modRemainder(std::move(lhs), rhs, prime);
        lhs = std::move(rhs);
        rhs = std::move(next);
    }
    if (lhs.empty())
        return lhs;
    const std::uint32_t inverse = modPower(lhs.back(), prime - 2, prime);
    for (std::uint32_t& coefficient : lhs)
        coefficient = modMultiply(coefficient, inverse, prime);
    return lhs;
}

[[nodiscard]] ModPolynomial modMultiplyReduce(
    const ModPolynomial& lhs,
    const ModPolynomial& rhs,
    const ModPolynomial& modulusPolynomial,
    std::uint32_t prime) {
    if (lhs.empty() || rhs.empty())
        return {};
    ModPolynomial product(lhs.size() + rhs.size() - 1);
    for (std::size_t i = 0; i < lhs.size(); ++i)
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            const std::uint32_t amount = modMultiply(lhs[i], rhs[j], prime);
            product[i + j] = static_cast<std::uint32_t>(
                (static_cast<std::uint64_t>(product[i + j]) + amount) % prime);
        }
    return modRemainder(std::move(product), modulusPolynomial, prime);
}

[[nodiscard]] ModPolynomial modPowerPolynomial(
    ModPolynomial base,
    std::uint64_t exponent,
    const ModPolynomial& modulusPolynomial,
    std::uint32_t prime) {
    ModPolynomial result{1};
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = modMultiplyReduce(result, base, modulusPolynomial, prime);
        exponent >>= 1U;
        if (exponent != 0)
            base = modMultiplyReduce(base, base, modulusPolynomial, prime);
    }
    return result;
}

[[nodiscard]] std::vector<std::size_t> distinctPrimeFactors(std::size_t value) {
    std::vector<std::size_t> result;
    for (std::size_t divisor = 2; divisor <= value / divisor; ++divisor) {
        if (value % divisor != 0)
            continue;
        result.push_back(divisor);
        while (value % divisor == 0)
            value /= divisor;
    }
    if (value > 1)
        result.push_back(value);
    return result;
}

[[nodiscard]] bool irreducibleModuloPrime(
    const Polynomial& polynomial,
    std::uint32_t prime) {
    const auto reduced = polynomialModulo(polynomial, prime);
    if (!reduced || reduced->size() <= 1)
        return false;
    const std::size_t degree = reduced->size() - 1;
    ModPolynomial x{0, 1};
    ModPolynomial frobenius = x;
    std::vector<ModPolynomial> powers(degree + 1);
    powers[0] = x;
    for (std::size_t i = 1; i <= degree; ++i) {
        frobenius = modPowerPolynomial(
            std::move(frobenius), prime, *reduced, prime);
        powers[i] = frobenius;
    }
    if (powers[degree] != x)
        return false;

    for (const std::size_t divisor : distinctPrimeFactors(degree)) {
        ModPolynomial difference = powers[degree / divisor];
        if (difference.size() < 2)
            difference.resize(2);
        difference[1] = difference[1] == 0 ? prime - 1 : difference[1] - 1;
        normalizeMod(difference);
        const ModPolynomial common = modGcd(*reduced, std::move(difference), prime);
        if (common.size() > 1)
            return false;
    }
    return true;
}

[[nodiscard]] bool provenIrreducibleOverQ(const Polynomial& polynomial) {
    if (polynomial.size() <= 2)
        return polynomial.size() == 2;
    constexpr std::uint32_t primes[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
    };
    for (const std::uint32_t prime : primes)
        if (irreducibleModuloPrime(polynomial, prime))
            return true;
    return false;
}

struct IntegerPolynomial final {
    std::vector<BigInt> coefficients;
};

[[nodiscard]] IntegerPolynomial primitiveIntegerPolynomial(const Polynomial& polynomial) {
    BigInt denominatorLcm{1};
    for (const Rational& coefficient : polynomial)
        denominatorLcm = numeric::lcm(denominatorLcm, coefficient.denominator());

    std::vector<BigInt> coefficients;
    coefficients.reserve(polynomial.size());
    BigInt common;
    for (const Rational& coefficient : polynomial) {
        BigInt value = coefficient.numerator()
            * (denominatorLcm / coefficient.denominator());
        common = common.isZero() ? value.abs() : numeric::gcd(common, value.abs());
        coefficients.push_back(std::move(value));
    }
    if (!common.isZero() && common != BigInt{1})
        for (BigInt& coefficient : coefficients)
            coefficient /= common;
    if (!coefficients.empty() && coefficients.back().isNegative())
        for (BigInt& coefficient : coefficients)
            coefficient = -coefficient;
    return IntegerPolynomial{std::move(coefficients)};
}

[[nodiscard]] BigInt evaluateIntegerPolynomial(
    const IntegerPolynomial& polynomial,
    std::int64_t value) {
    BigInt result;
    const BigInt x{value};
    for (auto iterator = polynomial.coefficients.rbegin();
         iterator != polynomial.coefficients.rend(); ++iterator)
        result = result * x + *iterator;
    return result;
}

[[nodiscard]] std::optional<std::vector<std::uint64_t>> positiveDivisors(
    const BigInt& value) {
    const auto magnitude = numeric::tryToUint64(value.abs());
    if (!magnitude || *magnitude == 0)
        return std::nullopt;
    std::vector<std::uint64_t> factors;
    if (!numeric::factorUint64(*magnitude, factors))
        return std::nullopt;
    std::sort(factors.begin(), factors.end());
    std::vector<std::pair<std::uint64_t, std::size_t>> groups;
    for (const std::uint64_t factor : factors) {
        if (!groups.empty() && groups.back().first == factor)
            ++groups.back().second;
        else
            groups.emplace_back(factor, 1);
    }
    std::vector<std::uint64_t> divisors{1};
    for (const auto& [prime, exponent] : groups) {
        const std::size_t existing = divisors.size();
        std::uint64_t power = 1;
        for (std::size_t e = 1; e <= exponent; ++e) {
            if (power > std::numeric_limits<std::uint64_t>::max() / prime)
                return std::nullopt;
            power *= prime;
            if (divisors.size() + existing > maximumKroneckerDivisors)
                return std::nullopt;
            for (std::size_t i = 0; i < existing; ++i) {
                if (divisors[i] > std::numeric_limits<std::uint64_t>::max() / power)
                    return std::nullopt;
                divisors.push_back(divisors[i] * power);
            }
        }
    }
    std::sort(divisors.begin(), divisors.end());
    return divisors;
}

[[nodiscard]] Polynomial multiplyByLinear(
    const Polynomial& polynomial,
    const Rational& root) {
    Polynomial result(polynomial.size() + 1);
    for (std::size_t i = 0; i < polynomial.size(); ++i) {
        result[i] -= polynomial[i] * root;
        result[i + 1] += polynomial[i];
    }
    normalize(result);
    return result;
}

[[nodiscard]] Polynomial interpolatePoints(
    const std::vector<std::int64_t>& points,
    const std::vector<BigInt>& values) {
    Polynomial result(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        Polynomial basis{Rational{BigInt{1}}};
        BigInt denominator{1};
        for (std::size_t j = 0; j < points.size(); ++j) {
            if (i == j)
                continue;
            basis = multiplyByLinear(basis, Rational{BigInt{points[j]}});
            denominator *= BigInt{points[i] - points[j]};
        }
        const Rational scale{values[i], denominator};
        for (std::size_t k = 0; k < basis.size(); ++k)
            result[k] += basis[k] * scale;
    }
    normalize(result);
    return result;
}

struct KroneckerSample final {
    std::int64_t x = 0;
    std::vector<std::uint64_t> divisors;
};

[[nodiscard]] std::optional<std::pair<Polynomial, Polynomial>> kroneckerSplit(
    const Polynomial& polynomial) {
    const std::size_t degree = polynomial.size() - 1;
    if (degree < 2 || degree > maximumMinimalPolynomialDegree)
        return std::nullopt;
    const IntegerPolynomial integerPolynomial = primitiveIntegerPolynomial(polynomial);

    for (std::int64_t x = -8; x <= 8; ++x) {
        if (evaluateIntegerPolynomial(integerPolynomial, x).isZero()) {
            Polynomial factor{-Rational{BigInt{x}}, Rational{BigInt{1}}};
            Polynomial quotient = divideExact(polynomial, factor);
            return std::pair<Polynomial, Polynomial>{
                canonicalPolynomial(factor), canonicalPolynomial(quotient)};
        }
    }

    std::vector<KroneckerSample> pool;
    for (std::int64_t radius = 0; radius <= 24; ++radius) {
        const std::array<std::int64_t, 2> candidates{radius, -radius};
        for (std::size_t candidateIndex = 0; candidateIndex < candidates.size(); ++candidateIndex) {
            if (radius == 0 && candidateIndex != 0)
                continue;
            const std::int64_t x = candidates[candidateIndex];
            const BigInt value = evaluateIntegerPolynomial(integerPolynomial, x);
            if (value.isZero())
                continue;
            const auto divisors = positiveDivisors(value);
            if (divisors)
                pool.push_back(KroneckerSample{x, *divisors});
        }
    }
    std::sort(pool.begin(), pool.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.divisors.size() < rhs.divisors.size();
    });

    for (std::size_t factorDegree = 1; factorDegree <= degree / 2; ++factorDegree) {
        if (pool.size() < factorDegree + 1)
            break;
        std::vector<KroneckerSample> samples(
            pool.begin(), pool.begin() + static_cast<std::ptrdiff_t>(factorDegree + 1));
        std::vector<std::int64_t> points;
        points.reserve(samples.size());
        for (const auto& sample : samples)
            points.push_back(sample.x);
        std::vector<BigInt> values(samples.size());
        std::size_t combinations = 0;
        std::optional<std::pair<Polynomial, Polynomial>> found;

        std::function<void(std::size_t)> search = [&](std::size_t index) {
            if (found || combinations >= maximumKroneckerCombinations)
                return;
            if (index == samples.size()) {
                ++combinations;
                Polynomial candidate = interpolatePoints(points, values);
                if (candidate.size() <= 1 || candidate.size() >= polynomial.size())
                    return;
                for (const Rational& coefficient : candidate)
                    if (!coefficient.isInteger())
                        return;
                candidate = canonicalPolynomial(candidate);
                if (candidate.size() <= 1 || candidate.size() >= polynomial.size())
                    return;
                if (!remainder(polynomial, candidate).empty())
                    return;
                Polynomial quotient = divideExact(polynomial, candidate);
                if (quotient.size() <= 1)
                    return;
                found = std::pair<Polynomial, Polynomial>{
                    std::move(candidate), canonicalPolynomial(quotient)};
                return;
            }

            for (const std::uint64_t divisor : samples[index].divisors) {
                const BigInt positive = BigInt::fromUnsigned(divisor);
                values[index] = positive;
                search(index + 1);
                if (found || combinations >= maximumKroneckerCombinations)
                    return;
                if (index != 0) {
                    values[index] = -positive;
                    search(index + 1);
                    if (found || combinations >= maximumKroneckerCombinations)
                        return;
                }
            }
        };
        search(0);
        if (found)
            return found;
    }
    return std::nullopt;
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


struct ReducedRealRoot final {
    Polynomial polynomial;
    std::size_t rootIndex = 0;
    RationalRootInterval interval;
    bool minimalPolynomialProven = false;
};

[[nodiscard]] std::optional<std::pair<std::size_t, RationalRootInterval>>
findRealRootInFactor(
    const Polynomial& factor,
    const RationalRootInterval& target) {
    if (target.isPoint()) {
        if (!evaluate(factor, target.lower).isZero())
            return std::nullopt;
    }
    else {
        const auto sturm = sturmSequence(factor);
        if (rootsBetween(sturm, target.lower, target.upper) != 1)
            return std::nullopt;
    }

    const auto intervals = isolateIntervals(factor);
    if (!intervals)
        return std::nullopt;
    const auto sturm = sturmSequence(factor);
    for (std::size_t i = 0; i < intervals->size(); ++i) {
        const RationalRootInterval& interval = (*intervals)[i];
        if (target.isPoint()) {
            if (interval.lower <= target.lower && target.lower <= interval.upper)
                return std::pair<std::size_t, RationalRootInterval>{i + 1, interval};
            continue;
        }
        const Rational lower = target.lower < interval.lower ? interval.lower : target.lower;
        const Rational upper = target.upper < interval.upper ? target.upper : interval.upper;
        if (lower < upper && rootsBetween(sturm, lower, upper) == 1)
            return std::pair<std::size_t, RationalRootInterval>{i + 1, interval};
    }
    return std::nullopt;
}

[[nodiscard]] ReducedRealRoot reduceRealRootPolynomial(
    Polynomial polynomial,
    std::size_t rootIndex,
    RationalRootInterval interval) {
    while (polynomial.size() > 2 && polynomial.size() - 1 <= maximumMinimalPolynomialDegree) {
        if (provenIrreducibleOverQ(polynomial))
            return ReducedRealRoot{
                std::move(polynomial), rootIndex, std::move(interval), true};
        const auto split = kroneckerSplit(polynomial);
        if (!split)
            break;
        const auto left = findRealRootInFactor(split->first, interval);
        const auto right = findRealRootInFactor(split->second, interval);
        if (left && !right) {
            polynomial = split->first;
            rootIndex = left->first;
            interval = left->second;
            continue;
        }
        if (right && !left) {
            polynomial = split->second;
            rootIndex = right->first;
            interval = right->second;
            continue;
        }
        break;
    }
    const bool proven = polynomial.size() == 2 || provenIrreducibleOverQ(polynomial);
    return ReducedRealRoot{
        std::move(polynomial), rootIndex, std::move(interval), proven};
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
    ReducedRealRoot reduced = reduceRealRootPolynomial(
        std::move(normalizedPolynomial), rootIndex, (*intervals)[rootIndex - 1]);
    return RealAlgebraicNumber{
        std::move(reduced.polynomial), reduced.rootIndex, std::move(reduced.interval)};
}

std::optional<RealAlgebraicNumber> RealAlgebraicNumber::createFromMinimalPolynomialInterval(
    std::span<const Rational> polynomial,
    RationalRootInterval interval) {
    if (interval.upper < interval.lower)
        return std::nullopt;

    Polynomial normalizedPolynomial = normalized(polynomial);
    if (normalizedPolynomial.size() <= 1
        || normalizedPolynomial.size() - 1 > maximumAlgebraicDegree)
        return std::nullopt;
    normalizedPolynomial = canonicalPolynomial(normalizedPolynomial);
    if (normalizedPolynomial.size() > 2 && !provenIrreducibleOverQ(normalizedPolynomial))
        return std::nullopt;

    if (interval.isPoint()) {
        if (!evaluate(normalizedPolynomial, interval.lower).isZero())
            return std::nullopt;
        if (normalizedPolynomial.size() != 2)
            return std::nullopt;
        return RealAlgebraicNumber{
            std::move(normalizedPolynomial), 1, std::move(interval)};
    }
    if (evaluate(normalizedPolynomial, interval.lower).isZero()
        || evaluate(normalizedPolynomial, interval.upper).isZero())
        return std::nullopt;

    const auto sturm = sturmSequence(normalizedPolynomial);
    if (rootsBetween(sturm, interval.lower, interval.upper) != 1)
        return std::nullopt;

    Rational bound = cauchyBound(normalizedPolynomial);
    Rational lowerBound = -bound;
    while (evaluate(normalizedPolynomial, lowerBound).isZero()) {
        bound += Rational{BigInt{1}};
        lowerBound = -bound;
    }
    const std::size_t rootIndex = rootsBetween(sturm, lowerBound, interval.lower) + 1;
    return RealAlgebraicNumber{
        std::move(normalizedPolynomial), rootIndex, std::move(interval)};
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


namespace {

using numeric::BigFloat;
using numeric::RoundingMode;

struct RationalComplex final {
    Rational real;
    Rational imaginary;
};

[[nodiscard]] RationalComplex rcAdd(const RationalComplex& a, const RationalComplex& b) {
    return {a.real + b.real, a.imaginary + b.imaginary};
}
[[nodiscard]] RationalComplex rcMultiply(const RationalComplex& a, const RationalComplex& b) {
    return {
        a.real * b.real - a.imaginary * b.imaginary,
        a.real * b.imaginary + a.imaginary * b.real};
}
[[nodiscard]] RationalComplex rcScale(const RationalComplex& a, const Rational& s) {
    return {a.real * s, a.imaginary * s};
}
[[nodiscard]] Rational rcL1(const RationalComplex& a) {
    return absolute(a.real) + absolute(a.imaginary);
}
[[nodiscard]] Rational rcMagnitudeSquared(const RationalComplex& a) {
    return a.real * a.real + a.imaginary * a.imaginary;
}
[[nodiscard]] RationalComplex rcPower(RationalComplex base, std::size_t exponent) {
    RationalComplex result{Rational{BigInt{1}}, Rational{}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = rcMultiply(result, base);
        exponent >>= 1U;
        if (exponent != 0)
            base = rcMultiply(base, base);
    }
    return result;
}

[[nodiscard]] std::uint64_t binomialSmall(std::size_t n, std::size_t k) {
    k = std::min(k, n - k);
    std::uint64_t result = 1;
    for (std::size_t i = 1; i <= k; ++i)
        result = (result * static_cast<std::uint64_t>(n - k + i)) / static_cast<std::uint64_t>(i);
    return result;
}

[[nodiscard]] std::vector<RationalComplex> taylorAt(
    const Polynomial& polynomial,
    const RationalComplex& center) {
    const std::size_t degree = polynomial.size() - 1;
    std::vector<RationalComplex> result(degree + 1);
    for (std::size_t k = 0; k <= degree; ++k) {
        RationalComplex sum{};
        for (std::size_t j = k; j <= degree; ++j) {
            const Rational choose{BigInt::fromUnsigned(binomialSmall(j, k))};
            const RationalComplex term = rcScale(
                rcPower(center, j - k), polynomial[j] * choose);
            sum = rcAdd(sum, term);
        }
        result[k] = std::move(sum);
    }
    return result;
}

[[nodiscard]] bool certifiesUniqueRoot(
    const Polynomial& polynomial,
    const RationalComplexDisk& disk) {
    if (disk.radius <= Rational{})
        return false;
    const auto taylor = taylorAt(polynomial, {disk.real, disk.imaginary});
    if (taylor.size() < 2)
        return false;

    Rational remainderBound = rcL1(taylor[0]);
    Rational radiusPower = disk.radius * disk.radius;
    for (std::size_t k = 2; k < taylor.size(); ++k) {
        remainderBound += rcL1(taylor[k]) * radiusPower;
        radiusPower *= disk.radius;
    }
    const Rational derivativeSquared = rcMagnitudeSquared(taylor[1]);
    if (derivativeSquared.isZero())
        return false;
    // |q0| + sum_{k>=2}|qk|r^k < |q1|r なら，境界上で線形項が残差を支配する。
    // Roucheによりdisk内のpの零点数は線形多項式と同じ1個になる。
    return remainderBound * remainderBound
        < derivativeSquared * disk.radius * disk.radius;
}

[[nodiscard]] bool disksDisjoint(
    const RationalComplexDisk& a,
    const RationalComplexDisk& b) {
    const Rational dx = a.real - b.real;
    const Rational dy = a.imaginary - b.imaginary;
    const Rational radius = a.radius + b.radius;
    return dx * dx + dy * dy > radius * radius;
}

struct ApproxComplex final {
    BigFloat real;
    BigFloat imaginary;
};

[[nodiscard]] BigFloat bf(const Rational& value, std::size_t bits) {
    return BigFloat::fromRational(value, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat bfAdd(const BigFloat& a, const BigFloat& b, std::size_t bits) {
    return numeric::add(a, b, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat bfSub(const BigFloat& a, const BigFloat& b, std::size_t bits) {
    return numeric::subtract(a, b, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat bfMul(const BigFloat& a, const BigFloat& b, std::size_t bits) {
    return numeric::multiply(a, b, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat bfDiv(const BigFloat& a, const BigFloat& b, std::size_t bits) {
    return numeric::divide(a, b, bits, RoundingMode::NearestEven);
}
[[nodiscard]] ApproxComplex acSub(const ApproxComplex& a, const ApproxComplex& b, std::size_t bits) {
    return {bfSub(a.real, b.real, bits), bfSub(a.imaginary, b.imaginary, bits)};
}
[[nodiscard]] ApproxComplex acMul(const ApproxComplex& a, const ApproxComplex& b, std::size_t bits) {
    const BigFloat ac = bfMul(a.real, b.real, bits);
    const BigFloat bd = bfMul(a.imaginary, b.imaginary, bits);
    const BigFloat ad = bfMul(a.real, b.imaginary, bits);
    const BigFloat bc = bfMul(a.imaginary, b.real, bits);
    return {bfSub(ac, bd, bits), bfAdd(ad, bc, bits)};
}
[[nodiscard]] ApproxComplex acDiv(const ApproxComplex& a, const ApproxComplex& b, std::size_t bits) {
    const BigFloat denominator = bfAdd(
        bfMul(b.real, b.real, bits), bfMul(b.imaginary, b.imaginary, bits), bits);
    if (denominator.isZero())
        throw std::domain_error("Complex algebraic root iteration encountered a zero denominator");
    return {
        bfDiv(bfAdd(bfMul(a.real, b.real, bits), bfMul(a.imaginary, b.imaginary, bits), bits), denominator, bits),
        bfDiv(bfSub(bfMul(a.imaginary, b.real, bits), bfMul(a.real, b.imaginary, bits), bits), denominator, bits)};
}
[[nodiscard]] ApproxComplex acScale(const ApproxComplex& a, const BigFloat& s, std::size_t bits) {
    return {bfMul(a.real, s, bits), bfMul(a.imaginary, s, bits)};
}

[[nodiscard]] ApproxComplex evaluateApprox(
    const Polynomial& polynomial,
    const ApproxComplex& z,
    std::size_t bits) {
    ApproxComplex result{bf(Rational{}, bits), bf(Rational{}, bits)};
    for (auto it = polynomial.rbegin(); it != polynomial.rend(); ++it) {
        result = acMul(result, z, bits);
        result.real = bfAdd(result.real, bf(*it, bits), bits);
    }
    return result;
}

[[nodiscard]] ApproxComplex evaluateDerivativeApprox(
    const Polynomial& polynomial,
    const ApproxComplex& z,
    std::size_t bits) {
    ApproxComplex result{bf(Rational{}, bits), bf(Rational{}, bits)};
    for (std::size_t exponent = polynomial.size() - 1; exponent > 0; --exponent) {
        result = acMul(result, z, bits);
        const Rational coefficient = polynomial[exponent]
            * Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(exponent))};
        result.real = bfAdd(result.real, bf(coefficient, bits), bits);
    }
    return result;
}

[[nodiscard]] BigFloat intervalMidpoint(
    const approximation::RealInterval& interval,
    std::size_t bits) {
    const Rational midpoint = (interval.lower().toRational() + interval.upper().toRational())
        / Rational{BigInt{2}};
    return bf(midpoint, bits);
}

[[nodiscard]] std::vector<ApproxComplex> durandKernerCandidates(
    const Polynomial& polynomial,
    std::size_t bits) {
    const std::size_t degree = polynomial.size() - 1;
    const Rational radiusRational = cauchyBound(polynomial) + Rational{BigInt{1}};
    const BigFloat radius = bf(radiusRational, bits);
    std::vector<ApproxComplex> roots;
    roots.reserve(degree);
    for (std::size_t k = 0; k < degree; ++k) {
        const Rational turns{
            BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k + 1)),
            BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * degree))};
        const auto cosine = approximation::encloseCosTurns(turns, bits + 16).interval;
        const auto sine = approximation::encloseSinTurns(turns, bits + 16).interval;
        roots.push_back(acScale(
            {intervalMidpoint(cosine, bits), intervalMidpoint(sine, bits)}, radius, bits));
    }

    const std::size_t iterations = std::min<std::size_t>(4096, 128 + bits * 2);
    for (std::size_t iteration = 0; iteration < iterations; ++iteration) {
        std::vector<ApproxComplex> next = roots;
        bool tinyCorrection = true;
        const Rational threshold = Rational{BigInt{1}, BigInt{1} << std::min<std::size_t>(bits / 2, 4096)};
        for (std::size_t i = 0; i < degree; ++i) {
            ApproxComplex denominator{bf(Rational{BigInt{1}}, bits), bf(Rational{}, bits)};
            for (std::size_t j = 0; j < degree; ++j)
                if (i != j)
                    denominator = acMul(denominator, acSub(roots[i], roots[j], bits), bits);
            if (denominator.real.isZero() && denominator.imaginary.isZero()) {
                tinyCorrection = false;
                continue;
            }
            const ApproxComplex correction = acDiv(evaluateApprox(polynomial, roots[i], bits), denominator, bits);
            next[i] = acSub(roots[i], correction, bits);
            const Rational correctionL1 = absolute(correction.real.toRational())
                + absolute(correction.imaginary.toRational());
            if (correctionL1 > threshold)
                tinyCorrection = false;
        }
        roots = std::move(next);
        if (tinyCorrection && iteration > degree)
            break;
    }
    return roots;
}

[[nodiscard]] std::optional<RationalComplexDisk> certifiedDiskAround(
    const Polynomial& polynomial,
    ApproxComplex center,
    std::size_t bits) {
    // Newtonで中心を数回磨き，その後はexact Rational中心へ固定してRouche判定する。
    for (std::size_t i = 0; i < 12; ++i) {
        const ApproxComplex derivative = evaluateDerivativeApprox(polynomial, center, bits);
        if (derivative.real.isZero() && derivative.imaginary.isZero())
            break;
        center = acSub(center,
            acDiv(evaluateApprox(polynomial, center, bits), derivative, bits), bits);
    }
    const RationalComplex rationalCenter{center.real.toRational(), center.imaginary.toRational()};
    const ApproxComplex derivative = evaluateDerivativeApprox(polynomial, center, bits);
    Rational correctionL1;
    if (!(derivative.real.isZero() && derivative.imaginary.isZero())) {
        const ApproxComplex correction = acDiv(evaluateApprox(polynomial, center, bits), derivative, bits);
        correctionL1 = absolute(correction.real.toRational())
            + absolute(correction.imaginary.toRational());
    }
    const std::size_t floorBits = std::min<std::size_t>(bits / 2, 4096);
    const Rational floorRadius{BigInt{1}, BigInt{1} << floorBits};
    Rational base = correctionL1 * Rational{BigInt{2}};
    if (base < floorRadius)
        base = floorRadius;

    for (const std::int64_t multiplier : std::array<std::int64_t, 8>{1,2,4,8,16,32,64,128}) {
        RationalComplexDisk disk{
            rationalCenter.real, rationalCenter.imaginary,
            base * Rational{BigInt{multiplier}}};
        if (certifiesUniqueRoot(polynomial, disk))
            return disk;
    }
    return std::nullopt;
}

struct RationalIntervalPair final { Rational lower; Rational upper; };

[[nodiscard]] RationalIntervalPair multiplyIntervals(
    const RationalIntervalPair& a,
    const RationalIntervalPair& b) {
    const std::array<Rational,4> values{
        a.lower*b.lower, a.lower*b.upper, a.upper*b.lower, a.upper*b.upper};
    auto [minIt, maxIt] = std::minmax_element(values.begin(), values.end());
    return {*minIt, *maxIt};
}

[[nodiscard]] RationalIntervalPair orderingKeyInterval(
    const RationalComplexDisk& disk,
    const RationalIntervalPair& pi) {
    const RationalIntervalPair imaginary{
        disk.imaginary - disk.radius, disk.imaginary + disk.radius};
    const auto piImag = multiplyIntervals(pi, imaginary);
    return {
        disk.real - disk.radius + piImag.lower,
        disk.real + disk.radius + piImag.upper};
}

[[nodiscard]] std::optional<std::vector<RationalComplexDisk>> isolateComplexDisks(
    const Polynomial& polynomial,
    std::size_t minimumBits = 128) {
    if (polynomial.size() <= 1)
        return std::vector<RationalComplexDisk>{};
    if (polynomial.size() - 1 > maximumAlgebraicDegree)
        return std::nullopt;
    const std::size_t degree = polynomial.size() - 1;

    for (std::size_t bits = std::max<std::size_t>(128, minimumBits); bits <= 4096; bits *= 2) {
        std::vector<ApproxComplex> candidates;
        try {
            candidates = durandKernerCandidates(polynomial, bits);
        }
        catch (const std::exception&) {
            continue;
        }
        if (candidates.size() != degree)
            continue;

        std::vector<RationalComplexDisk> disks;
        disks.reserve(degree);
        bool failed = false;
        for (ApproxComplex& candidate : candidates) {
            const auto disk = certifiedDiskAround(polynomial, std::move(candidate), bits);
            if (!disk) {
                failed = true;
                break;
            }
            disks.push_back(*disk);
        }
        if (failed)
            continue;

        for (std::size_t i = 0; i < disks.size() && !failed; ++i)
            for (std::size_t j = i + 1; j < disks.size(); ++j)
                if (!disksDisjoint(disks[i], disks[j])) {
                    failed = true;
                    break;
                }
        if (failed)
            continue;

        const auto piBox = approximation::enclosePi(bits + 32).interval;
        const RationalIntervalPair pi{
            piBox.lower().toRational(), piBox.upper().toRational()};
        std::sort(disks.begin(), disks.end(), [&](const RationalComplexDisk& lhs, const RationalComplexDisk& rhs) {
            const Rational lhsMid = lhs.real + ((pi.lower + pi.upper) / Rational{BigInt{2}}) * lhs.imaginary;
            const Rational rhsMid = rhs.real + ((pi.lower + pi.upper) / Rational{BigInt{2}}) * rhs.imaginary;
            return lhsMid < rhsMid;
        });
        bool ordered = true;
        for (std::size_t i = 1; i < disks.size(); ++i) {
            const auto left = orderingKeyInterval(disks[i - 1], pi);
            const auto right = orderingKeyInterval(disks[i], pi);
            if (!(left.upper < right.lower)) {
                ordered = false;
                break;
            }
        }
        if (ordered)
            return disks;
    }
    return std::nullopt;
}

[[nodiscard]] bool complexDiskNarrowEnough(
    const RationalComplexDisk& disk,
    std::size_t precisionBits) {
    const Rational scale = Rational{BigInt{1}}
        + absolute(disk.real) + absolute(disk.imaginary);
    return disk.radius * powerOfTwo(precisionBits + 8) <= scale;
}


struct ReducedComplexRoot final {
    Polynomial polynomial;
    std::size_t rootIndex = 0;
    RationalComplexDisk disk;
    bool minimalPolynomialProven = false;
};

[[nodiscard]] std::optional<std::pair<std::size_t, RationalComplexDisk>>
findComplexRootInFactor(
    const Polynomial& factor,
    const RationalComplexDisk& target) {
    for (std::size_t bits : std::array<std::size_t, 5>{128,256,512,1024,2048}) {
        const auto disks = isolateComplexDisks(factor, bits);
        if (!disks)
            return std::nullopt;
        std::size_t matchCount = 0;
        std::size_t matchIndex = 0;
        RationalComplexDisk matchDisk{};
        for (std::size_t i = 0; i < disks->size(); ++i) {
            if (!disksDisjoint((*disks)[i], target)) {
                ++matchCount;
                matchIndex = i + 1;
                matchDisk = (*disks)[i];
            }
        }
        if (matchCount == 1)
            return std::pair<std::size_t, RationalComplexDisk>{matchIndex, matchDisk};
        if (matchCount == 0)
            return std::nullopt;
    }
    return std::nullopt;
}

[[nodiscard]] ReducedComplexRoot reduceComplexRootPolynomial(
    Polynomial polynomial,
    std::size_t rootIndex,
    RationalComplexDisk disk) {
    while (polynomial.size() > 2 && polynomial.size() - 1 <= maximumMinimalPolynomialDegree) {
        if (provenIrreducibleOverQ(polynomial))
            return ReducedComplexRoot{
                std::move(polynomial), rootIndex, std::move(disk), true};
        const auto split = kroneckerSplit(polynomial);
        if (!split)
            break;
        const auto left = findComplexRootInFactor(split->first, disk);
        const auto right = findComplexRootInFactor(split->second, disk);
        if (left && !right) {
            polynomial = split->first;
            rootIndex = left->first;
            disk = left->second;
            continue;
        }
        if (right && !left) {
            polynomial = split->second;
            rootIndex = right->first;
            disk = right->second;
            continue;
        }
        break;
    }
    const bool proven = polynomial.size() == 2 || provenIrreducibleOverQ(polynomial);
    return ReducedComplexRoot{
        std::move(polynomial), rootIndex, std::move(disk), proven};
}

} // namespace

ComplexAlgebraicNumber::ComplexAlgebraicNumber(
    std::vector<Rational> polynomial,
    std::size_t rootIndex,
    RationalComplexDisk disk)
    : polynomial_(std::move(polynomial)), rootIndex_(rootIndex), disk_(std::move(disk)) {}

std::optional<ComplexAlgebraicNumber> ComplexAlgebraicNumber::create(
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
    const auto disks = isolateComplexDisks(normalizedPolynomial);
    if (!disks || rootIndex > disks->size())
        return std::nullopt;
    ReducedComplexRoot reduced = reduceComplexRootPolynomial(
        std::move(normalizedPolynomial), rootIndex, (*disks)[rootIndex - 1]);
    return ComplexAlgebraicNumber{
        std::move(reduced.polynomial), reduced.rootIndex, std::move(reduced.disk)};
}

std::optional<std::vector<ComplexAlgebraicNumber>> ComplexAlgebraicNumber::isolateAll(
    std::span<const Rational> polynomial) {
    Polynomial normalizedPolynomial = normalized(polynomial);
    if (normalizedPolynomial.size() <= 1)
        return std::vector<ComplexAlgebraicNumber>{};
    if (normalizedPolynomial.size() - 1 > maximumAlgebraicDegree)
        return std::nullopt;
    normalizedPolynomial = canonicalPolynomial(normalizedPolynomial);
    const auto disks = isolateComplexDisks(normalizedPolynomial);
    if (!disks)
        return std::nullopt;
    std::vector<ComplexAlgebraicNumber> roots;
    roots.reserve(disks->size());
    for (std::size_t i = 0; i < disks->size(); ++i)
        roots.push_back(ComplexAlgebraicNumber{normalizedPolynomial, i + 1, (*disks)[i]});
    return roots;
}

std::span<const Rational> ComplexAlgebraicNumber::polynomial() const noexcept { return polynomial_; }
std::size_t ComplexAlgebraicNumber::degree() const noexcept { return polynomial_.size() - 1; }
std::size_t ComplexAlgebraicNumber::rootIndex() const noexcept { return rootIndex_; }
const RationalComplexDisk& ComplexAlgebraicNumber::isolatingDisk() const noexcept { return disk_; }

RationalComplexDisk ComplexAlgebraicNumber::refined(std::size_t precisionBits) const {
    if (complexDiskNarrowEnough(disk_, precisionBits))
        return disk_;
    const auto disks = isolateComplexDisks(polynomial_, std::min<std::size_t>(4096, precisionBits + 64));
    if (!disks || rootIndex_ > disks->size())
        throw std::runtime_error("Complex algebraic root refinement did not converge");
    return (*disks)[rootIndex_ - 1];
}


namespace {

constexpr std::size_t maximumAlgebraicArithmeticDegree = 16;

[[nodiscard]] Rational determinant(std::vector<std::vector<Rational>> matrix) {
    if (matrix.empty())
        return Rational{BigInt{1}};
    const std::size_t size = matrix.size();
    Rational result{BigInt{1}};
    bool negative = false;
    for (std::size_t column = 0; column < size; ++column) {
        std::size_t pivot = column;
        while (pivot < size && matrix[pivot][column].isZero())
            ++pivot;
        if (pivot == size)
            return Rational{};
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            negative = !negative;
        }
        const Rational pivotValue = matrix[column][column];
        result *= pivotValue;
        for (std::size_t row = column + 1; row < size; ++row) {
            if (matrix[row][column].isZero())
                continue;
            const Rational factor = matrix[row][column] / pivotValue;
            for (std::size_t k = column + 1; k < size; ++k)
                matrix[row][k] -= factor * matrix[column][k];
            matrix[row][column] = Rational{};
        }
    }
    return negative ? -result : result;
}

[[nodiscard]] Rational resultant(const Polynomial& lhs, const Polynomial& rhs) {
    if (lhs.size() <= 1 || rhs.size() <= 1)
        return Rational{};
    const std::size_t m = lhs.size() - 1;
    const std::size_t n = rhs.size() - 1;
    const std::size_t size = m + n;
    std::vector<std::vector<Rational>> matrix(size, std::vector<Rational>(size));
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t j = 0; j <= m; ++j)
            matrix[row][row + j] = lhs[j];
    for (std::size_t row = 0; row < m; ++row)
        for (std::size_t j = 0; j <= n; ++j)
            matrix[n + row][row + j] = rhs[j];
    return determinant(std::move(matrix));
}

[[nodiscard]] BigInt binomialInteger(std::size_t n, std::size_t k) {
    if (k > n)
        return BigInt{};
    k = std::min(k, n - k);
    BigInt result{1};
    for (std::size_t i = 1; i <= k; ++i) {
        result *= BigInt::fromUnsigned(n - k + i);
        result /= BigInt::fromUnsigned(i);
    }
    return result;
}

[[nodiscard]] Rational integerPowerRational(Rational base, std::size_t exponent) {
    Rational result{BigInt{1}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result *= base;
        exponent >>= 1U;
        if (exponent != 0)
            base *= base;
    }
    return result;
}

[[nodiscard]] Polynomial shiftedNegatedArgument(const Polynomial& polynomial, const Rational& shift) {
    // q(shift-y) をyの昇冪係数へ展開する。
    Polynomial result(polynomial.size());
    for (std::size_t j = 0; j < polynomial.size(); ++j) {
        for (std::size_t k = 0; k <= j; ++k) {
            Rational coefficient = polynomial[j]
                * Rational{binomialInteger(j, k)}
                * integerPowerRational(shift, j - k);
            if ((k & 1U) != 0)
                coefficient = -coefficient;
            result[k] += coefficient;
        }
    }
    normalize(result);
    return result;
}

[[nodiscard]] Polynomial multiplicationTransform(const Polynomial& polynomial, const Rational& value) {
    // y^n q(value/y) = sum_j q_j value^j y^(n-j)
    const std::size_t degree = polynomial.size() - 1;
    Polynomial result(degree + 1);
    Rational valuePower{BigInt{1}};
    for (std::size_t j = 0; j <= degree; ++j) {
        result[degree - j] += polynomial[j] * valuePower;
        valuePower *= value;
    }
    normalize(result);
    return result;
}

[[nodiscard]] Polynomial multiplyByLinearFactor(
    const Polynomial& polynomial,
    const Rational& constant) {
    // polynomial * (x-constant)
    Polynomial result(polynomial.size() + 1);
    for (std::size_t i = 0; i < polynomial.size(); ++i) {
        result[i] -= polynomial[i] * constant;
        result[i + 1] += polynomial[i];
    }
    normalize(result);
    return result;
}

[[nodiscard]] Polynomial interpolateIntegerGrid(const std::vector<Rational>& values) {
    if (values.empty())
        return {};
    std::vector<Rational> differences = values;
    std::vector<Rational> forward;
    forward.reserve(values.size());
    for (std::size_t order = 0; order < values.size(); ++order) {
        forward.push_back(differences.front());
        for (std::size_t i = 0; i + 1 < differences.size(); ++i)
            differences[i] = differences[i + 1] - differences[i];
        differences.pop_back();
    }

    Polynomial result(values.size());
    Polynomial falling{Rational{BigInt{1}}};
    BigInt factorial{1};
    for (std::size_t k = 0; k < forward.size(); ++k) {
        const Rational scale = forward[k] / Rational{factorial};
        for (std::size_t i = 0; i < falling.size(); ++i)
            result[i] += falling[i] * scale;
        if (k + 1 < forward.size()) {
            falling = multiplyByLinearFactor(falling, Rational{BigInt::fromUnsigned(k)});
            factorial *= BigInt::fromUnsigned(k + 1);
        }
    }
    normalize(result);
    return result;
}

[[nodiscard]] Polynomial resultantPolynomial(
    const Polynomial& lhs,
    const Polynomial& rhs,
    AlgebraicBinaryOperation operation) {
    Polynomial right = rhs;
    if (operation == AlgebraicBinaryOperation::Subtract) {
        for (std::size_t i = 1; i < right.size(); i += 2)
            right[i] = -right[i];
        operation = AlgebraicBinaryOperation::Add;
    }
    else if (operation == AlgebraicBinaryOperation::Divide) {
        std::reverse(right.begin(), right.end());
        normalize(right); // zero rootは逆数を持たないため自然に次数から落ちる。
        if (right.size() <= 1)
            return {};
        operation = AlgebraicBinaryOperation::Multiply;
    }

    const std::size_t lhsDegree = lhs.size() - 1;
    const std::size_t rhsDegree = right.size() - 1;
    const std::size_t degreeBound = lhsDegree * rhsDegree;
    if (degreeBound == 0 || degreeBound > maximumAlgebraicArithmeticDegree)
        return {};

    std::vector<Rational> values;
    values.reserve(degreeBound + 1);
    for (std::size_t sample = 0; sample <= degreeBound; ++sample) {
        const Rational x{BigInt::fromUnsigned(sample)};
        Polynomial transformed = operation == AlgebraicBinaryOperation::Add
            ? shiftedNegatedArgument(right, x)
            : multiplicationTransform(right, x);
        values.push_back(resultant(lhs, transformed));
    }
    return canonicalPolynomial(interpolateIntegerGrid(values));
}

[[nodiscard]] RationalComplexDisk asDisk(const AlgebraicNumber& value, std::size_t bits) {
    if (const auto* real = value.asReal()) {
        const RationalRootInterval interval = real->refined(bits);
        const Rational two{BigInt{2}};
        return RationalComplexDisk{
            (interval.lower + interval.upper) / two,
            Rational{},
            (interval.upper - interval.lower) / two};
    }
    return value.asComplex()->refined(bits);
}

[[nodiscard]] Rational complexL1(const RationalComplexDisk& disk) {
    return absolute(disk.real) + absolute(disk.imaginary);
}

[[nodiscard]] RationalComplexDisk addDisks(
    const RationalComplexDisk& lhs,
    const RationalComplexDisk& rhs,
    bool subtract) {
    return RationalComplexDisk{
        subtract ? lhs.real - rhs.real : lhs.real + rhs.real,
        subtract ? lhs.imaginary - rhs.imaginary : lhs.imaginary + rhs.imaginary,
        lhs.radius + rhs.radius};
}

[[nodiscard]] RationalComplexDisk scaleDisk(
    const RationalComplexDisk& value,
    const Rational& scale) {
    return RationalComplexDisk{
        value.real * scale,
        value.imaginary * scale,
        value.radius * absolute(scale)};
}

[[nodiscard]] RationalComplexDisk multiplyDisks(
    const RationalComplexDisk& lhs,
    const RationalComplexDisk& rhs) {
    const Rational real = lhs.real * rhs.real - lhs.imaginary * rhs.imaginary;
    const Rational imaginary = lhs.real * rhs.imaginary + lhs.imaginary * rhs.real;
    const Rational radius = complexL1(lhs) * rhs.radius
        + complexL1(rhs) * lhs.radius + lhs.radius * rhs.radius;
    return RationalComplexDisk{real, imaginary, radius};
}

[[nodiscard]] std::optional<RationalComplexDisk> reciprocalDisk(
    const RationalComplexDisk& value) {
    const Rational absReal = absolute(value.real);
    const Rational absImaginary = absolute(value.imaginary);
    const Rational magnitudeLower = absReal < absImaginary ? absImaginary : absReal;
    if (magnitudeLower <= value.radius)
        return std::nullopt;
    const Rational denominator = value.real * value.real + value.imaginary * value.imaginary;
    if (denominator.isZero())
        return std::nullopt;
    const Rational real = value.real / denominator;
    const Rational imaginary = -value.imaginary / denominator;
    const Rational distanceLower = magnitudeLower - value.radius;
    const Rational radius = value.radius / (magnitudeLower * distanceLower);
    return RationalComplexDisk{real, imaginary, radius};
}

[[nodiscard]] std::optional<RationalComplexDisk> imageDisk(
    const AlgebraicNumber& lhs,
    const AlgebraicNumber& rhs,
    AlgebraicBinaryOperation operation,
    std::size_t bits) {
    const RationalComplexDisk left = asDisk(lhs, bits);
    const RationalComplexDisk right = asDisk(rhs, bits);
    switch (operation) {
    case AlgebraicBinaryOperation::Add:
        return addDisks(left, right, false);
    case AlgebraicBinaryOperation::Subtract:
        return addDisks(left, right, true);
    case AlgebraicBinaryOperation::Multiply:
        return multiplyDisks(left, right);
    case AlgebraicBinaryOperation::Divide: {
        const auto reciprocal = reciprocalDisk(right);
        if (!reciprocal)
            return std::nullopt;
        return multiplyDisks(left, *reciprocal);
    }
    }
    return std::nullopt;
}

[[nodiscard]] bool disksIntersect(
    const RationalComplexDisk& lhs,
    const RationalComplexDisk& rhs) {
    const Rational dx = lhs.real - rhs.real;
    const Rational dy = lhs.imaginary - rhs.imaginary;
    const Rational radius = lhs.radius + rhs.radius;
    return dx * dx + dy * dy <= radius * radius;
}


using RationalVector = std::vector<Rational>;
using RationalMatrix = std::vector<RationalVector>;

[[nodiscard]] RationalVector tensorBasisElement(
    std::size_t lhsDegree,
    std::size_t rhsDegree,
    std::size_t lhsExponent,
    std::size_t rhsExponent) {
    RationalVector value(lhsDegree * rhsDegree);
    value[lhsExponent * rhsDegree + rhsExponent] = Rational{BigInt{1}};
    return value;
}

[[nodiscard]] RationalVector tensorMultiply(
    const RationalVector& lhs,
    const RationalVector& rhs,
    const Polynomial& lhsPolynomial,
    const Polynomial& rhsPolynomial) {
    const std::size_t m = lhsPolynomial.size() - 1;
    const std::size_t n = rhsPolynomial.size() - 1;
    if (lhs.size() != m * n || rhs.size() != m * n)
        throw std::logic_error("Primitive-element tensor size mismatch");

    const std::size_t alphaCount = 2 * m - 1;
    const std::size_t betaCount = 2 * n - 1;
    std::vector<Rational> temporary(alphaCount * betaCount);
    auto at = [&](std::size_t alpha, std::size_t beta) -> Rational& {
        return temporary[alpha * betaCount + beta];
    };
    for (std::size_t ai = 0; ai < m; ++ai)
        for (std::size_t bi = 0; bi < n; ++bi) {
            const Rational& left = lhs[ai * n + bi];
            if (left.isZero())
                continue;
            for (std::size_t aj = 0; aj < m; ++aj)
                for (std::size_t bj = 0; bj < n; ++bj) {
                    const Rational& right = rhs[aj * n + bj];
                    if (!right.isZero())
                        at(ai + aj, bi + bj) += left * right;
                }
        }

    for (std::int64_t alpha = static_cast<std::int64_t>(alphaCount) - 1;
         alpha >= static_cast<std::int64_t>(m); --alpha) {
        for (std::size_t beta = 0; beta < betaCount; ++beta) {
            Rational coefficient = at(static_cast<std::size_t>(alpha), beta);
            if (coefficient.isZero())
                continue;
            at(static_cast<std::size_t>(alpha), beta) = Rational{};
            const std::size_t shift = static_cast<std::size_t>(alpha) - m;
            for (std::size_t k = 0; k < m; ++k)
                at(shift + k, beta) -= coefficient * lhsPolynomial[k];
        }
    }
    for (std::int64_t beta = static_cast<std::int64_t>(betaCount) - 1;
         beta >= static_cast<std::int64_t>(n); --beta) {
        for (std::size_t alpha = 0; alpha < m; ++alpha) {
            Rational coefficient = at(alpha, static_cast<std::size_t>(beta));
            if (coefficient.isZero())
                continue;
            at(alpha, static_cast<std::size_t>(beta)) = Rational{};
            const std::size_t shift = static_cast<std::size_t>(beta) - n;
            for (std::size_t k = 0; k < n; ++k)
                at(alpha, shift + k) -= coefficient * rhsPolynomial[k];
        }
    }

    RationalVector result(m * n);
    for (std::size_t alpha = 0; alpha < m; ++alpha)
        for (std::size_t beta = 0; beta < n; ++beta)
            result[alpha * n + beta] = at(alpha, beta);
    return result;
}

[[nodiscard]] std::optional<Polynomial> tensorMinimalPolynomial(
    const RationalVector& element,
    const Polynomial& lhsPolynomial,
    const Polynomial& rhsPolynomial) {
    const std::size_t dimension = element.size();
    RationalVector current(dimension);
    current.front() = Rational{BigInt{1}};

    detail::RationalLinearBasis krylovBasis(dimension);
    if (krylovBasis.append(current))
        return std::nullopt;

    for (std::size_t degree = 1; degree <= dimension; ++degree) {
        current = tensorMultiply(current, element, lhsPolynomial, rhsPolynomial);
        if (auto relation = krylovBasis.append(current))
            return canonicalPolynomial(*relation);
    }
    return std::nullopt;
}

struct PrimitiveElementReduction final {
    std::shared_ptr<const NumberFieldContext> field;
    AlgebraicElement alpha;
    AlgebraicElement beta;
};

constexpr std::size_t maximumCachedPrimitiveElementReductions = 64;

class PrimitiveElementReductionCache final {
public:
    [[nodiscard]] std::optional<PrimitiveElementReduction> find(
        const AlgebraicNumber& lhs,
        const AlgebraicNumber& rhs) {
        std::lock_guard lock(mutex_);
        pruneExpired();
        for (std::size_t i = 0; i < entries_.size(); ++i) {
            Entry& entry = entries_[i];
            const bool direct = entry.lhs.hasSameRootIdentity(lhs)
                && entry.rhs.hasSameRootIdentity(rhs);
            const bool reversed = entry.lhs.hasSameRootIdentity(rhs)
                && entry.rhs.hasSameRootIdentity(lhs);
            if (!direct && !reversed)
                continue;

            auto field = entry.field.lock();
            if (!field)
                continue;
            auto alpha = AlgebraicElement::create(
                field, direct ? entry.alphaCoordinates : entry.betaCoordinates);
            auto beta = AlgebraicElement::create(
                field, direct ? entry.betaCoordinates : entry.alphaCoordinates);
            if (!alpha || !beta)
                continue;

            if (i + 1 != entries_.size()) {
                Entry hit = std::move(entry);
                entries_.erase(entries_.begin() + static_cast<std::ptrdiff_t>(i));
                entries_.push_back(std::move(hit));
            }
            return PrimitiveElementReduction{
                std::move(field), std::move(*alpha), std::move(*beta)};
        }
        return std::nullopt;
    }

    [[nodiscard]] PrimitiveElementReduction publish(
        const AlgebraicNumber& lhs,
        const AlgebraicNumber& rhs,
        PrimitiveElementReduction reduction) {
        std::lock_guard lock(mutex_);
        pruneExpired();
        for (std::size_t i = 0; i < entries_.size(); ++i) {
            Entry& entry = entries_[i];
            const bool direct = entry.lhs.hasSameRootIdentity(lhs)
                && entry.rhs.hasSameRootIdentity(rhs);
            const bool reversed = entry.lhs.hasSameRootIdentity(rhs)
                && entry.rhs.hasSameRootIdentity(lhs);
            if (!direct && !reversed)
                continue;
            auto field = entry.field.lock();
            if (!field)
                continue;
            auto alpha = AlgebraicElement::create(
                field, direct ? entry.alphaCoordinates : entry.betaCoordinates);
            auto beta = AlgebraicElement::create(
                field, direct ? entry.betaCoordinates : entry.alphaCoordinates);
            if (alpha && beta)
                return PrimitiveElementReduction{
                    std::move(field), std::move(*alpha), std::move(*beta)};
        }

        if (entries_.size() >= maximumCachedPrimitiveElementReductions)
            entries_.erase(entries_.begin());
        entries_.push_back(Entry{
            lhs.withArithmeticElement({}),
            rhs.withArithmeticElement({}),
            reduction.field,
            std::vector<Rational>{
                reduction.alpha.coefficients().begin(), reduction.alpha.coefficients().end()},
            std::vector<Rational>{
                reduction.beta.coefficients().begin(), reduction.beta.coefficients().end()}});
        return reduction;
    }

private:
    struct Entry final {
        AlgebraicNumber lhs;
        AlgebraicNumber rhs;
        std::weak_ptr<const NumberFieldContext> field;
        RationalVector alphaCoordinates;
        RationalVector betaCoordinates;
    };

    std::mutex mutex_;
    std::vector<Entry> entries_;

    void pruneExpired() {
        entries_.erase(
            std::remove_if(entries_.begin(), entries_.end(),
                [](const Entry& entry) { return entry.field.expired(); }),
            entries_.end());
    }
};

[[nodiscard]] PrimitiveElementReductionCache& primitiveElementReductionCache() {
    static PrimitiveElementReductionCache cache;
    return cache;
}

[[nodiscard]] std::optional<AlgebraicNumber> selectPrimitiveGenerator(
    const Polynomial& polynomial,
    const AlgebraicNumber& lhs,
    const AlgebraicNumber& rhs,
    std::int64_t multiplier) {
    const AlgebraicRootDomain domain =
        lhs.domain() == AlgebraicRootDomain::Real && rhs.domain() == AlgebraicRootDomain::Real
        ? AlgebraicRootDomain::Real : AlgebraicRootDomain::Complex;
    const Rational scale{BigInt{multiplier}};

    if (domain == AlgebraicRootDomain::Real) {
        const auto roots = RealAlgebraicNumber::isolateAll(polynomial);
        if (!roots)
            return std::nullopt;
        for (std::size_t bits : std::array<std::size_t, 5>{128,256,512,1024,2048}) {
            const RationalComplexDisk target = addDisks(
                asDisk(lhs, bits), scaleDisk(asDisk(rhs, bits), scale), false);
            std::size_t match = 0;
            std::size_t matchIndex = 0;
            for (const RealAlgebraicNumber& root : *roots) {
                const RationalRootInterval interval = root.refined(bits);
                const Rational two{BigInt{2}};
                const RationalComplexDisk disk{
                    (interval.lower + interval.upper) / two,
                    Rational{},
                    (interval.upper - interval.lower) / two};
                if (disksIntersect(disk, target)) {
                    ++match;
                    matchIndex = root.rootIndex();
                }
            }
            if (match == 1)
                return AlgebraicNumber::create(polynomial, matchIndex, domain);
        }
        return std::nullopt;
    }

    const auto roots = ComplexAlgebraicNumber::isolateAll(polynomial);
    if (!roots)
        return std::nullopt;
    for (std::size_t bits : std::array<std::size_t, 5>{128,256,512,1024,2048}) {
        const RationalComplexDisk target = addDisks(
            asDisk(lhs, bits), scaleDisk(asDisk(rhs, bits), scale), false);
        std::size_t match = 0;
        std::size_t matchIndex = 0;
        for (const ComplexAlgebraicNumber& root : *roots) {
            if (disksIntersect(root.refined(bits), target)) {
                ++match;
                matchIndex = root.rootIndex();
            }
        }
        if (match == 1)
            return AlgebraicNumber::create(polynomial, matchIndex, domain);
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<PrimitiveElementReduction> primitiveElementReduction(
    const RationalVector& alpha,
    const RationalVector& beta,
    const Polynomial& lhsPolynomial,
    const Polynomial& rhsPolynomial,
    const AlgebraicNumber& lhs,
    const AlgebraicNumber& rhs) {
    if (auto cached = primitiveElementReductionCache().find(lhs, rhs))
        return cached;

    const std::size_t dimension = alpha.size();
    constexpr std::int64_t candidates[] = {1, 2, -1, 3, -2, 4, -3};
    for (const std::int64_t multiplier : candidates) {
        RationalVector theta(dimension);
        const Rational scale{BigInt{multiplier}};
        for (std::size_t i = 0; i < dimension; ++i)
            theta[i] = alpha[i] + scale * beta[i];
        const auto polynomial = tensorMinimalPolynomial(theta, lhsPolynomial, rhsPolynomial);
        if (!polynomial || polynomial->size() - 1 != dimension
            || !provenIrreducibleOverQ(*polynomial))
            continue;

        detail::RationalLinearBasis powerBasis(dimension);
        RationalVector current(dimension);
        current.front() = Rational{BigInt{1}};
        bool basisComplete = true;
        for (std::size_t exponent = 0; exponent < dimension; ++exponent) {
            if (powerBasis.append(current)) {
                basisComplete = false;
                break;
            }
            current = tensorMultiply(current, theta, lhsPolynomial, rhsPolynomial);
        }
        if (!basisComplete || powerBasis.rank() != dimension)
            continue;

        // degree==dimensionの既約relationがあるので1,theta,...,theta^(d-1)はbasis。
        auto alphaCoordinates = powerBasis.coordinates(alpha);
        auto betaCoordinates = powerBasis.coordinates(beta);
        if (!alphaCoordinates || !betaCoordinates)
            continue;

        auto generator = selectPrimitiveGenerator(*polynomial, lhs, rhs, multiplier);
        if (!generator)
            continue;
        auto field = NumberFieldContext::create(std::move(*generator));
        if (!field)
            continue;

        auto alphaElement = AlgebraicElement::create(field, std::move(*alphaCoordinates));
        auto betaElement = AlgebraicElement::create(field, std::move(*betaCoordinates));
        if (!alphaElement || !betaElement)
            continue;
        return primitiveElementReductionCache().publish(
            lhs, rhs, PrimitiveElementReduction{
                std::move(field), std::move(*alphaElement), std::move(*betaElement)});
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<AlgebraicElement> rationalFieldElement(
    const AlgebraicElement& reference,
    const Rational& value) {
    RationalVector coefficients(reference.field()->degree());
    coefficients.front() = value;
    return AlgebraicElement::create(reference.field(), std::move(coefficients));
}

[[nodiscard]] std::optional<AlgebraicElement> fieldOperationElement(
    const AlgebraicElement& lhs,
    const AlgebraicElement& rhs,
    AlgebraicBinaryOperation operation) {
    switch (operation) {
    case AlgebraicBinaryOperation::Add:
        return lhs.add(rhs);
    case AlgebraicBinaryOperation::Subtract:
        return lhs.subtract(rhs);
    case AlgebraicBinaryOperation::Multiply:
        return lhs.multiply(rhs);
    case AlgebraicBinaryOperation::Divide:
        return lhs.divide(rhs);
    }
    return std::nullopt;
}

[[nodiscard]] RationalComplexDisk intervalDisk(const RationalRootInterval& interval) {
    const Rational two{BigInt{2}};
    return RationalComplexDisk{
        (interval.lower + interval.upper) / two,
        Rational{},
        (interval.upper - interval.lower) / two};
}

[[nodiscard]] std::optional<AlgebraicNumber> selectResultRoot(
    const Polynomial& polynomial,
    AlgebraicRootDomain domain,
    const AlgebraicNumber& lhs,
    const AlgebraicNumber& rhs,
    AlgebraicBinaryOperation operation) {
    if (polynomial.size() <= 1 || polynomial.size() - 1 > maximumAlgebraicDegree)
        return std::nullopt;

    if (domain == AlgebraicRootDomain::Real) {
        const auto roots = RealAlgebraicNumber::isolateAll(polynomial);
        if (!roots)
            return std::nullopt;
        for (std::size_t bits : std::array<std::size_t, 5>{128,256,512,1024,2048}) {
            const auto target = imageDisk(lhs, rhs, operation, bits);
            if (!target)
                continue;
            std::size_t match = 0;
            std::size_t matchIndex = 0;
            for (const RealAlgebraicNumber& root : *roots) {
                if (disksIntersect(intervalDisk(root.refined(bits)), *target)) {
                    ++match;
                    matchIndex = root.rootIndex();
                }
            }
            if (match == 1)
                return AlgebraicNumber::create(polynomial, matchIndex, domain);
        }
        return std::nullopt;
    }

    const auto roots = ComplexAlgebraicNumber::isolateAll(polynomial);
    if (!roots)
        return std::nullopt;
    for (std::size_t bits : std::array<std::size_t, 5>{128,256,512,1024,2048}) {
        const auto target = imageDisk(lhs, rhs, operation, bits);
        if (!target)
            continue;
        std::size_t match = 0;
        std::size_t matchIndex = 0;
        for (const ComplexAlgebraicNumber& root : *roots) {
            if (disksIntersect(root.refined(bits), *target)) {
                ++match;
                matchIndex = root.rootIndex();
            }
        }
        if (match == 1)
            return AlgebraicNumber::create(polynomial, matchIndex, domain);
    }
    return std::nullopt;
}


[[nodiscard]] std::optional<AlgebraicNumber> materializeFieldElement(
    AlgebraicElement element,
    const AlgebraicNumber& lhs,
    const AlgebraicNumber& rhs,
    AlgebraicBinaryOperation operation) {
    const auto minimal = element.minimalPolynomial();
    if (!minimal || minimal->size() <= 1)
        return std::nullopt;
    if (!provenIrreducibleOverQ(*minimal) && minimal->size() > 2)
        return std::nullopt;

    const AlgebraicRootDomain resultDomain =
        lhs.domain() == AlgebraicRootDomain::Real && rhs.domain() == AlgebraicRootDomain::Real
        ? AlgebraicRootDomain::Real : AlgebraicRootDomain::Complex;

    if (resultDomain == AlgebraicRootDomain::Real) {
        for (const std::size_t bits : std::array<std::size_t, 7>{32,64,128,256,512,1024,2048}) {
            auto interval = element.refinedRealInterval(bits);
            if (!interval)
                break;
            auto root = RealAlgebraicNumber::createFromMinimalPolynomialInterval(
                *minimal, std::move(*interval));
            if (!root)
                continue;
            AlgebraicNumber result = AlgebraicNumber::fromRealRoot(std::move(*root));
            return result.withArithmeticElement(
                std::make_shared<const AlgebraicElement>(std::move(element)));
        }
    }

    auto result = selectResultRoot(*minimal, resultDomain, lhs, rhs, operation);
    if (!result)
        return std::nullopt;
    return result->withArithmeticElement(
        std::make_shared<const AlgebraicElement>(std::move(element)));
}

[[nodiscard]] std::optional<AlgebraicNumber> primitiveElementCombine(
    const AlgebraicNumber& lhs,
    const AlgebraicNumber& rhs,
    AlgebraicBinaryOperation operation) {
    Polynomial left(lhs.polynomial().begin(), lhs.polynomial().end());
    Polynomial right(rhs.polynomial().begin(), rhs.polynomial().end());
    if (left.size() <= 2 || right.size() <= 2)
        return std::nullopt;
    if (!provenIrreducibleOverQ(left) || !provenIrreducibleOverQ(right))
        return std::nullopt;

    const std::size_t m = left.size() - 1;
    const std::size_t n = right.size() - 1;
    const std::size_t dimension = m * n;
    if (dimension == 0 || dimension > maximumAlgebraicArithmeticDegree)
        return std::nullopt;

    const RationalVector alpha = tensorBasisElement(m, n, 1, 0);
    const RationalVector beta = tensorBasisElement(m, n, 0, 1);
    const auto primitive = primitiveElementReduction(
        alpha, beta, left, right, lhs, rhs);
    if (!primitive)
        return std::nullopt;

    auto resultElement = fieldOperationElement(
        primitive->alpha, primitive->beta, operation);
    if (!resultElement)
        return std::nullopt;
    return materializeFieldElement(
        std::move(*resultElement), lhs, rhs, operation);
}


[[nodiscard]] BigInt floorRational(const Rational& value) {
    BigInt quotient = value.numerator() / value.denominator();
    if (value.numerator().isNegative()
        && !(value.numerator() % value.denominator()).isZero())
        quotient -= BigInt{1};
    return quotient;
}

[[nodiscard]] Rational simplestPositive(const Rational& lower, const Rational& upper) {
    if (!(Rational{} < lower) || upper < lower)
        throw std::logic_error("Invalid positive rational interval");
    if (lower.isInteger())
        return lower;
    const BigInt lowerFloor = floorRational(lower);
    const BigInt upperFloor = floorRational(upper);
    if (lowerFloor != upperFloor)
        return Rational{lowerFloor + BigInt{1}};
    const Rational integerPart{lowerFloor};
    const Rational lowFraction = lower - integerPart;
    const Rational highFraction = upper - integerPart;
    if (lowFraction.isZero())
        return lower;
    return integerPart + Rational{BigInt{1}} /
        simplestPositive(Rational{BigInt{1}} / highFraction,
            Rational{BigInt{1}} / lowFraction);
}

[[nodiscard]] Rational simplestInInterval(Rational lower, Rational upper) {
    if (upper < lower)
        std::swap(lower, upper);
    if (lower <= Rational{} && Rational{} <= upper)
        return Rational{};
    if (upper < Rational{})
        return -simplestPositive(-upper, -lower);
    return simplestPositive(lower, upper);
}

[[nodiscard]] RationalComplex evaluateComplex(
    const Polynomial& polynomial,
    const RationalComplex& value) {
    RationalComplex result{};
    for (auto iterator = polynomial.rbegin(); iterator != polynomial.rend(); ++iterator)
        result = rcAdd(rcMultiply(result, value), RationalComplex{*iterator, Rational{}});
    return result;
}

[[nodiscard]] bool hasSamePolynomial(
    const AlgebraicNumber& lhs,
    const AlgebraicNumber& rhs) noexcept {
    const auto left = lhs.polynomial();
    const auto right = rhs.polynomial();
    return left.size() == right.size()
        && std::equal(left.begin(), left.end(), right.begin());
}

[[nodiscard]] bool isExactlyZero(const AlgebraicNumber& value) {
    if (const AlgebraicElement* element = value.arithmeticElement())
        return element->isZero();

    const auto polynomial = value.polynomial();
    if (polynomial.empty() || !polynomial.front().isZero())
        return false;

    if (const RealAlgebraicNumber* real = value.asReal()) {
        const RationalRootInterval& interval = real->isolatingInterval();
        return interval.lower <= Rational{} && Rational{} <= interval.upper;
    }

    const RationalComplexDisk& disk = value.asComplex()->isolatingDisk();
    return disk.real * disk.real + disk.imaginary * disk.imaginary
        <= disk.radius * disk.radius;
}

[[nodiscard]] std::optional<AlgebraicSign> exactRealSign(
    const AlgebraicNumber& value) {
    if (value.domain() != AlgebraicRootDomain::Real)
        return std::nullopt;

    if (const AlgebraicElement* element = value.arithmeticElement()) {
        if (const auto sign = element->exactSign())
            return sign;
    }
    if (isExactlyZero(value))
        return AlgebraicSign::Zero;

    const RealAlgebraicNumber* real = value.asReal();
    if (!real)
        return std::nullopt;
    for (std::size_t bits = 32; bits <= maximumComparisonRefinementBits; bits *= 2) {
        const RationalRootInterval interval = real->refined(bits);
        if (interval.upper < Rational{})
            return AlgebraicSign::Negative;
        if (Rational{} < interval.lower)
            return AlgebraicSign::Positive;
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<AlgebraicElement> embedRational(
    const AlgebraicElement& exemplar,
    const Rational& value) {
    std::vector<Rational> coordinates(exemplar.field()->degree());
    if (coordinates.empty())
        return std::nullopt;
    coordinates.front() = value;
    return AlgebraicElement::create(exemplar.field(), std::move(coordinates));
}

} // namespace

AlgebraicNumber::AlgebraicNumber(RealAlgebraicNumber value) : value_(std::move(value)) {}
AlgebraicNumber::AlgebraicNumber(ComplexAlgebraicNumber value) : value_(std::move(value)) {}

std::optional<AlgebraicNumber> AlgebraicNumber::create(
    std::span<const Rational> polynomial,
    std::size_t rootIndex,
    AlgebraicRootDomain domain) {
    if (domain == AlgebraicRootDomain::Real) {
        auto value = RealAlgebraicNumber::create(polynomial, rootIndex);
        if (!value)
            return std::nullopt;
        return AlgebraicNumber{std::move(*value)};
    }
    auto value = ComplexAlgebraicNumber::create(polynomial, rootIndex);
    if (!value)
        return std::nullopt;
    return AlgebraicNumber{std::move(*value)};
}

AlgebraicNumber AlgebraicNumber::fromRealRoot(RealAlgebraicNumber value) {
    return AlgebraicNumber{std::move(value)};
}

AlgebraicNumber AlgebraicNumber::fromComplexRoot(ComplexAlgebraicNumber value) {
    return AlgebraicNumber{std::move(value)};
}

std::optional<AlgebraicNumber> AlgebraicNumber::fromRational(const Rational& value) {
    const std::array<Rational, 2> polynomial{-value, Rational{BigInt{1}}};
    return create(polynomial, 1, AlgebraicRootDomain::Real);
}

std::optional<AlgebraicNumber> AlgebraicNumber::fromComplexRational(
    const Rational& real,
    const Rational& imaginary) {
    if (imaginary.isZero())
        return fromRational(real);
    const std::array<Rational, 3> polynomial{
        real * real + imaginary * imaginary,
        -Rational{BigInt{2}} * real,
        Rational{BigInt{1}}};
    const auto roots = ComplexAlgebraicNumber::isolateAll(polynomial);
    if (!roots)
        return std::nullopt;
    const RationalComplexDisk point{real, imaginary, Rational{}};
    for (const ComplexAlgebraicNumber& root : *roots)
        if (disksIntersect(root.isolatingDisk(), point))
            return AlgebraicNumber{root};
    return std::nullopt;
}

std::optional<AlgebraicNumber> AlgebraicNumber::combine(
    const AlgebraicNumber& lhs,
    const AlgebraicNumber& rhs,
    AlgebraicBinaryOperation operation) {
    std::shared_ptr<const AlgebraicElement> leftElement = lhs.arithmeticElement_;
    std::shared_ptr<const AlgebraicElement> rightElement = rhs.arithmeticElement_;

    // 同一canonical Root identityなら，片側が保持するfield座標をそのまま共有できる。
    // pointer identityが異なるContextを数学的に同一視する一般mergeはここでは行わない。
    if (lhs.hasSameRootIdentity(rhs)) {
        if (leftElement)
            rightElement = leftElement;
        else if (rightElement)
            leftElement = rightElement;
    }

    // Qは任意のQ(theta)へ定数座標として埋め込める。Rational側が別の内部表現を
    // 持っていても，既存fieldを優先して同一Context演算へ落とす。
    if (leftElement) {
        const auto exact = rhs.exactRationalParts();
        if (exact && exact->second.isZero()) {
            auto embedded = rationalFieldElement(*leftElement, exact->first);
            auto fieldResult = embedded
                ? fieldOperationElement(*leftElement, *embedded, operation)
                : std::nullopt;
            if (fieldResult) {
                if (auto materialized = materializeFieldElement(
                        std::move(*fieldResult), lhs, rhs, operation))
                    return materialized;
            }
        }
    }
    if (rightElement) {
        const auto exact = lhs.exactRationalParts();
        if (exact && exact->second.isZero()) {
            auto embedded = rationalFieldElement(*rightElement, exact->first);
            auto fieldResult = embedded
                ? fieldOperationElement(*embedded, *rightElement, operation)
                : std::nullopt;
            if (fieldResult) {
                if (auto materialized = materializeFieldElement(
                        std::move(*fieldResult), lhs, rhs, operation))
                    return materialized;
            }
        }
    }

    if (leftElement && rightElement) {
        auto fieldResult = fieldOperationElement(*leftElement, *rightElement, operation);
        if (fieldResult) {
            if (auto materialized = materializeFieldElement(
                    std::move(*fieldResult), lhs, rhs, operation))
                return materialized;
        }
    }

    if (const auto reduced = primitiveElementCombine(lhs, rhs, operation))
        return reduced;

    Polynomial left(lhs.polynomial().begin(), lhs.polynomial().end());
    Polynomial right(rhs.polynomial().begin(), rhs.polynomial().end());
    if (left.size() <= 1 || right.size() <= 1)
        return std::nullopt;
    if ((left.size() - 1) * (right.size() - 1) > maximumAlgebraicArithmeticDegree)
        return std::nullopt;

    Polynomial result = resultantPolynomial(left, right, operation);
    if (result.size() <= 1)
        return std::nullopt;
    const AlgebraicRootDomain resultDomain =
        lhs.domain() == AlgebraicRootDomain::Real && rhs.domain() == AlgebraicRootDomain::Real
        ? AlgebraicRootDomain::Real : AlgebraicRootDomain::Complex;
    return selectResultRoot(result, resultDomain, lhs, rhs, operation);
}

AlgebraicRootDomain AlgebraicNumber::domain() const noexcept {
    return std::holds_alternative<RealAlgebraicNumber>(value_)
        ? AlgebraicRootDomain::Real : AlgebraicRootDomain::Complex;
}
std::span<const Rational> AlgebraicNumber::polynomial() const noexcept {
    if (const auto* value = std::get_if<RealAlgebraicNumber>(&value_))
        return value->polynomial();
    return std::get<ComplexAlgebraicNumber>(value_).polynomial();
}
std::size_t AlgebraicNumber::rootIndex() const noexcept {
    if (const auto* value = std::get_if<RealAlgebraicNumber>(&value_))
        return value->rootIndex();
    return std::get<ComplexAlgebraicNumber>(value_).rootIndex();
}
const RealAlgebraicNumber* AlgebraicNumber::asReal() const noexcept {
    return std::get_if<RealAlgebraicNumber>(&value_);
}
const ComplexAlgebraicNumber* AlgebraicNumber::asComplex() const noexcept {
    return std::get_if<ComplexAlgebraicNumber>(&value_);
}

bool AlgebraicNumber::hasSameRootIdentity(const AlgebraicNumber& rhs) const noexcept {
    if (domain() != rhs.domain() || rootIndex() != rhs.rootIndex())
        return false;
    return hasSamePolynomial(*this, rhs);
}

std::optional<bool> AlgebraicNumber::exactEquals(const AlgebraicNumber& rhs) const {
    if (hasSameRootIdentity(rhs))
        return true;

    const bool samePolynomial = hasSamePolynomial(*this, rhs);
    // 同じ根集合を同じdomainで列挙しているなら異なるindexは異なる根である。
    if (domain() == rhs.domain() && samePolynomial)
        return false;

    // 異なるQ上既約minimal polynomialは共通根を持たない。
    if (!samePolynomial) {
        Polynomial leftPolynomial(polynomial().begin(), polynomial().end());
        Polynomial rightPolynomial(rhs.polynomial().begin(), rhs.polynomial().end());
        if (provenIrreducibleOverQ(leftPolynomial)
            && provenIrreducibleOverQ(rightPolynomial))
            return false;
    }

    const AlgebraicElement* leftElement = arithmeticElement();
    const AlgebraicElement* rightElement = rhs.arithmeticElement();
    if (leftElement && rightElement) {
        if (const auto equal = leftElement->exactEquals(*rightElement))
            return equal;
    }

    if (leftElement) {
        const auto right = rhs.exactRationalParts();
        if (right && right->second.isZero()) {
            const auto embedded = embedRational(*leftElement, right->first);
            if (embedded)
                if (const auto equal = leftElement->exactEquals(*embedded))
                    return equal;
        }
    }
    if (rightElement) {
        const auto left = exactRationalParts();
        if (left && left->second.isZero()) {
            const auto embedded = embedRational(*rightElement, left->first);
            if (embedded)
                if (const auto equal = embedded->exactEquals(*rightElement))
                    return equal;
        }
    }

    if (const auto left = exactRationalParts()) {
        if (const auto right = rhs.exactRationalParts())
            return *left == *right;
    }

    // 分離領域が既に交わらなければ，common fieldを構成せずFalseを証明できる。
    if (const RealAlgebraicNumber* realLeft = asReal()) {
        if (const RealAlgebraicNumber* realRight = rhs.asReal()) {
            const RationalRootInterval& a = realLeft->isolatingInterval();
            const RationalRootInterval& b = realRight->isolatingInterval();
            if (a.upper < b.lower || b.upper < a.lower)
                return false;
        }
    }
    else if (const ComplexAlgebraicNumber* complexLeft = asComplex()) {
        if (const ComplexAlgebraicNumber* complexRight = rhs.asComplex())
            if (disksDisjoint(complexLeft->isolatingDisk(), complexRight->isolatingDisk()))
                return false;
    }

    // 既存のbounded primitive-element/resultant machineryで差をexactに構成する。
    // 0判定はfield座標またはrootのcertified isolationだけを使う。
    const auto difference = combine(*this, rhs, AlgebraicBinaryOperation::Subtract);
    if (!difference)
        return std::nullopt;
    return isExactlyZero(*difference);
}

std::optional<AlgebraicOrder> AlgebraicNumber::exactRealCompare(
    const AlgebraicNumber& rhs) const {
    if (domain() != AlgebraicRootDomain::Real
        || rhs.domain() != AlgebraicRootDomain::Real)
        return std::nullopt;
    if (hasSameRootIdentity(rhs))
        return AlgebraicOrder::Equal;

    const AlgebraicElement* leftElement = arithmeticElement();
    const AlgebraicElement* rightElement = rhs.arithmeticElement();
    if (leftElement && rightElement
        && leftElement->field().get() == rightElement->field().get()) {
        const auto difference = leftElement->subtract(*rightElement);
        const auto sign = difference ? difference->exactSign() : std::nullopt;
        if (sign == AlgebraicSign::Negative) return AlgebraicOrder::Less;
        if (sign == AlgebraicSign::Zero) return AlgebraicOrder::Equal;
        if (sign == AlgebraicSign::Positive) return AlgebraicOrder::Greater;
    }

    // 異なるfieldでも実根のcertified intervalsが分離すれば即比較できる。
    const RealAlgebraicNumber* left = asReal();
    const RealAlgebraicNumber* right = rhs.asReal();
    if (!left || !right)
        return std::nullopt;
    for (std::size_t bits = 32; bits <= maximumComparisonRefinementBits; bits *= 2) {
        const RationalRootInterval a = left->refined(bits);
        const RationalRootInterval b = right->refined(bits);
        if (a.upper < b.lower) return AlgebraicOrder::Less;
        if (b.upper < a.lower) return AlgebraicOrder::Greater;
    }

    // intervalが重なり続ける場合は差をexact algebraic valueとして構成する。
    const auto difference = combine(*this, rhs, AlgebraicBinaryOperation::Subtract);
    if (!difference)
        return std::nullopt;
    const auto sign = exactRealSign(*difference);
    if (sign == AlgebraicSign::Negative) return AlgebraicOrder::Less;
    if (sign == AlgebraicSign::Zero) return AlgebraicOrder::Equal;
    if (sign == AlgebraicSign::Positive) return AlgebraicOrder::Greater;
    return std::nullopt;
}

const AlgebraicElement* AlgebraicNumber::arithmeticElement() const noexcept {
    return arithmeticElement_.get();
}

AlgebraicNumber AlgebraicNumber::withArithmeticElement(
    std::shared_ptr<const AlgebraicElement> element) const {
    AlgebraicNumber result = *this;
    result.arithmeticElement_ = std::move(element);
    return result;
}

AlgebraicNumber AlgebraicNumber::withGeneratorField() const {
    if (arithmeticElement_)
        return *this;

    Polynomial definingPolynomial(polynomial().begin(), polynomial().end());
    if (definingPolynomial.size() <= 1
        || !provenIrreducibleOverQ(definingPolynomial))
        return *this;

    auto field = NumberFieldContext::create(*this);
    auto element = field ? AlgebraicElement::generator(std::move(field)) : std::nullopt;
    if (!element)
        return *this;
    return withArithmeticElement(
        std::make_shared<const AlgebraicElement>(std::move(*element)));
}

std::optional<std::pair<Rational, Rational>> AlgebraicNumber::exactRationalParts() const {
    // arithmeticElement_はminimal polynomialの既約性を証明できたRootにだけ付与する。
    // したがって実次数>1はRationalではなく，複素次数>2はQ+iQへ退化しない。
    // 高精度refinementで毎回それを再確認する固定費を避ける。
    if (arithmeticElement_) {
        const std::size_t degree = polynomial().size() - 1;
        if ((domain() == AlgebraicRootDomain::Real && degree > 1)
            || (domain() == AlgebraicRootDomain::Complex && degree > 2))
            return std::nullopt;
    }

    if (const auto* real = asReal()) {
        const RationalRootInterval interval = real->refined(192);
        const Rational candidate = simplestInInterval(interval.lower, interval.upper);
        if (evaluate(Polynomial{real->polynomial().begin(), real->polynomial().end()}, candidate).isZero())
            return std::pair<Rational, Rational>{candidate, Rational{}};
        return std::nullopt;
    }

    const auto* complex = asComplex();
    const RationalComplexDisk disk = complex->refined(192);
    const Rational realCandidate = simplestInInterval(
        disk.real - disk.radius, disk.real + disk.radius);
    const Rational imaginaryCandidate = simplestInInterval(
        disk.imaginary - disk.radius, disk.imaginary + disk.radius);
    const RationalComplex candidate{realCandidate, imaginaryCandidate};
    const RationalComplex value = evaluateComplex(
        Polynomial{complex->polynomial().begin(), complex->polynomial().end()}, candidate);
    if (value.real.isZero() && value.imaginary.isZero())
        return std::pair<Rational, Rational>{realCandidate, imaginaryCandidate};
    return std::nullopt;
}

} // namespace mmcal::symbolic
