// 多変数多項式idealとGröbner basis
#include "groebner.hpp"

#include "evaluation/evaluation_budget.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace mmcal::symbolic {
namespace {

using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational one() { return Rational{BigInt{1}}; }

[[nodiscard]] MultivariateRationalPolynomial termPolynomial(const PolynomialTerm& term) {
    return MultivariateRationalPolynomial{{term}};
}

[[nodiscard]] PolynomialTerm quotientTerm(
    const PolynomialTerm& numerator,
    const PolynomialTerm& denominator) {
    const auto monomial = divideMonomials(numerator.monomial, denominator.monomial);
    if (!monomial)
        throw std::invalid_argument("Leading monomial is not divisible");
    return PolynomialTerm{*monomial, numerator.coefficient / denominator.coefficient};
}

void checkPolynomialSize(
    const MultivariateRationalPolynomial& polynomial,
    std::size_t maximumTerms) {
    if (polynomial.termCount() > maximumTerms)
        throw std::length_error("Groebner polynomial term limit exceeded");
}

[[nodiscard]] bool relativelyPrime(const Monomial& lhs, const Monomial& rhs) noexcept {
    for (const MonomialFactor& factor : lhs.factors())
        if (rhs.exponentOf(factor.variable) != 0)
            return false;
    return true;
}

struct PairKey final {
    std::size_t first = 0;
    std::size_t second = 0;

    [[nodiscard]] bool operator<(const PairKey& rhs) const noexcept {
        return std::tie(first, second) < std::tie(rhs.first, rhs.second);
    }
};

[[nodiscard]] PairKey pairKey(std::size_t lhs, std::size_t rhs) noexcept {
    return lhs < rhs ? PairKey{lhs, rhs} : PairKey{rhs, lhs};
}

struct CriticalPair final {
    std::size_t first = 0;
    std::size_t second = 0;
    Monomial lcm;
    std::size_t sugar = 0;
};

[[nodiscard]] bool pairLess(
    const CriticalPair& lhs,
    const CriticalPair& rhs,
    const PolynomialRing& ring) {
    if (lhs.sugar != rhs.sugar)
        return lhs.sugar < rhs.sugar;
    const int comparison = ring.compare(lhs.lcm, rhs.lcm);
    if (comparison != 0)
        return comparison < 0;
    return pairKey(lhs.first, lhs.second) < pairKey(rhs.first, rhs.second);
}

[[nodiscard]] bool chainCriterion(
    const CriticalPair& pair,
    std::span<const MultivariateRationalPolynomial> basis,
    const PolynomialRing& ring,
    const std::set<PairKey>& zeroPairs) {
    static_cast<void>(ring);
    for (std::size_t k = 0; k < basis.size(); ++k) {
        if (k == pair.first || k == pair.second)
            continue;
        const auto leading = basis[k].leadingTerm(ring);
        if (!leading || !monomialDivides(leading->monomial, pair.lcm))
            continue;
        if (zeroPairs.contains(pairKey(pair.first, k))
            && zeroPairs.contains(pairKey(k, pair.second)))
            return true;
    }
    return false;
}

void appendCriticalPair(
    std::vector<CriticalPair>& pairs,
    const std::vector<MultivariateRationalPolynomial>& basis,
    std::size_t first,
    std::size_t second,
    const PolynomialRing& ring,
    std::set<PairKey>& zeroPairs,
    std::size_t maximumPairs) {
    const auto left = basis[first].leadingTerm(ring);
    const auto right = basis[second].leadingTerm(ring);
    if (!left || !right)
        return;

    // Buchberger product criterion。leading monomialが互いに素ならS-polynomialは0へ還元される。
    if (relativelyPrime(left->monomial, right->monomial)) {
        zeroPairs.insert(pairKey(first, second));
        return;
    }

    if (pairs.size() >= maximumPairs)
        throw std::length_error("Groebner critical-pair limit exceeded");
    Monomial lcm = leastCommonMultiple(left->monomial, right->monomial);
    pairs.push_back(CriticalPair{
        first, second, std::move(lcm),
        std::max({basis[first].totalDegree(), basis[second].totalDegree(),
            leastCommonMultiple(left->monomial, right->monomial).totalDegree()})});
}

[[nodiscard]] bool samePolynomial(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs) {
    const auto left = lhs.terms();
    const auto right = rhs.terms();
    return std::equal(left.begin(), left.end(), right.begin(), right.end());
}

} // namespace

PolynomialDivisionResult multivariateDivide(
    const MultivariateRationalPolynomial& dividend,
    std::span<const MultivariateRationalPolynomial> divisors,
    const PolynomialRing& ring,
    std::size_t maximumReductionSteps) {
    if (!dividend.belongsTo(ring))
        throw std::invalid_argument("Dividend does not belong to the polynomial ring");
    for (const auto& divisor : divisors)
        if (!divisor.belongsTo(ring))
            throw std::invalid_argument("Divisor does not belong to the polynomial ring");

    PolynomialDivisionResult result;
    result.quotients.resize(divisors.size());
    MultivariateRationalPolynomial current = dividend;
    std::size_t steps = 0;

    while (!current.isZero()) {
        evaluation::checkEvaluationCancellation();
        if (++steps > maximumReductionSteps)
            throw std::length_error("Multivariate polynomial reduction step limit exceeded");
        const PolynomialTerm leading = *current.leadingTerm(ring);
        bool reduced = false;

        for (std::size_t i = 0; i < divisors.size(); ++i) {
            if (divisors[i].isZero())
                continue;
            const PolynomialTerm divisorLeading = *divisors[i].leadingTerm(ring);
            if (!monomialDivides(divisorLeading.monomial, leading.monomial))
                continue;

            const PolynomialTerm factor = quotientTerm(leading, divisorLeading);
            result.quotients[i] = addPolynomials(
                result.quotients[i], termPolynomial(factor));
            current = subtractPolynomials(
                current, multiplyPolynomialByTerm(divisors[i], factor));
            reduced = true;
            break;
        }

        if (reduced)
            continue;

        const MultivariateRationalPolynomial leadingPolynomial = termPolynomial(leading);
        result.remainder = addPolynomials(result.remainder, leadingPolynomial);
        current = subtractPolynomials(current, leadingPolynomial);
    }

    return result;
}

MultivariateRationalPolynomial normalForm(
    const MultivariateRationalPolynomial& polynomial,
    std::span<const MultivariateRationalPolynomial> basis,
    const PolynomialRing& ring,
    std::size_t maximumReductionSteps) {
    return multivariateDivide(polynomial, basis, ring, maximumReductionSteps).remainder;
}

MultivariateRationalPolynomial sPolynomial(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs,
    const PolynomialRing& ring) {
    const auto left = lhs.leadingTerm(ring);
    const auto right = rhs.leadingTerm(ring);
    if (!left || !right)
        return MultivariateRationalPolynomial{};

    const Monomial lcm = leastCommonMultiple(left->monomial, right->monomial);
    const Monomial leftMonomial = *divideMonomials(lcm, left->monomial);
    const Monomial rightMonomial = *divideMonomials(lcm, right->monomial);
    const PolynomialTerm leftFactor{leftMonomial, one() / left->coefficient};
    const PolynomialTerm rightFactor{rightMonomial, one() / right->coefficient};
    return subtractPolynomials(
        multiplyPolynomialByTerm(lhs, leftFactor),
        multiplyPolynomialByTerm(rhs, rightFactor));
}

std::vector<MultivariateRationalPolynomial> reducedGroebnerBasis(
    std::span<const MultivariateRationalPolynomial> source,
    const PolynomialRing& ring,
    std::size_t maximumReductionSteps) {
    std::vector<MultivariateRationalPolynomial> basis;
    basis.reserve(source.size());
    for (const auto& polynomial : source)
        if (!polynomial.isZero())
            basis.push_back(monicPolynomial(polynomial, ring));

    // leading idealに冗長なgeneratorを除く。
    std::vector<bool> remove(basis.size(), false);
    for (std::size_t i = 0; i < basis.size(); ++i) {
        const auto leadingI = basis[i].leadingTerm(ring);
        for (std::size_t j = 0; j < basis.size(); ++j) {
            if (i == j)
                continue;
            const auto leadingJ = basis[j].leadingTerm(ring);
            if (!leadingI || !leadingJ)
                continue;
            if (monomialDivides(leadingJ->monomial, leadingI->monomial)) {
                const bool same = leadingJ->monomial == leadingI->monomial;
                if (!same || j < i) {
                    remove[i] = true;
                    break;
                }
            }
        }
    }

    std::vector<MultivariateRationalPolynomial> minimal;
    for (std::size_t i = 0; i < basis.size(); ++i)
        if (!remove[i])
            minimal.push_back(std::move(basis[i]));

    // 各generatorを他のgeneratorで完全還元しmonic化する。
    std::vector<MultivariateRationalPolynomial> reduced;
    reduced.reserve(minimal.size());
    for (std::size_t i = 0; i < minimal.size(); ++i) {
        std::vector<MultivariateRationalPolynomial> others;
        others.reserve(minimal.size() - 1);
        for (std::size_t j = 0; j < minimal.size(); ++j)
            if (j != i)
                others.push_back(minimal[j]);
        MultivariateRationalPolynomial value = normalForm(
            minimal[i], others, ring, maximumReductionSteps);
        if (!value.isZero())
            reduced.push_back(monicPolynomial(value, ring));
    }

    // reductionで同一generatorが生じた場合を除去。
    std::vector<MultivariateRationalPolynomial> unique;
    for (auto& polynomial : reduced) {
        const bool duplicate = std::any_of(unique.begin(), unique.end(), [&](const auto& existing) {
            return samePolynomial(existing, polynomial);
        });
        if (!duplicate)
            unique.push_back(std::move(polynomial));
    }

    std::sort(unique.begin(), unique.end(), [&](const auto& lhs, const auto& rhs) {
        const auto left = lhs.leadingTerm(ring);
        const auto right = rhs.leadingTerm(ring);
        if (!left) return false;
        if (!right) return true;
        return ring.compare(left->monomial, right->monomial) > 0;
    });
    return unique;
}

GroebnerComputation groebnerBasis(
    std::span<const MultivariateRationalPolynomial> generators,
    const PolynomialRing& ring,
    GroebnerOptions options) {
    GroebnerComputation result;
    std::set<PairKey> zeroPairs;

    // generator自身を順次normal formへ落とし、初期basisの重複・冗長性を抑える。
    for (const auto& generator : generators) {
        if (!generator.belongsTo(ring))
            throw std::invalid_argument("Groebner generator does not belong to the polynomial ring");
        if (generator.isZero())
            continue;
        MultivariateRationalPolynomial value = normalForm(
            generator, result.basis, ring, options.maximumReductionSteps);
        if (value.isZero())
            continue;
        checkPolynomialSize(value, options.maximumTermsPerPolynomial);
        value = monicPolynomial(value, ring);
        if (result.basis.size() >= options.maximumBasisSize)
            throw std::length_error("Groebner basis size limit exceeded");
        result.basis.push_back(std::move(value));
    }

    // nonzero constantがidealへ入ればbasisは{1}で終了。
    for (const auto& polynomial : result.basis) {
        const auto leading = polynomial.leadingTerm(ring);
        if (leading && leading->monomial.isOne()) {
            result.basis = {MultivariateRationalPolynomial{{PolynomialTerm{Monomial{}, one()}}}};
            return result;
        }
    }

    std::vector<CriticalPair> pairs;
    for (std::size_t j = 0; j < result.basis.size(); ++j)
        for (std::size_t i = 0; i < j; ++i)
            appendCriticalPair(
                pairs, result.basis, i, j, ring, zeroPairs,
                options.maximumCriticalPairs);

    while (!pairs.empty()) {
        evaluation::consumeEvaluationBudget(evaluation::EvaluationResource::SolverBranch);
        const auto selected = std::min_element(
            pairs.begin(), pairs.end(), [&](const CriticalPair& lhs, const CriticalPair& rhs) {
                return pairLess(lhs, rhs, ring);
            });
        CriticalPair pair = std::move(*selected);
        pairs.erase(selected);
        if (++result.criticalPairsProcessed > options.maximumCriticalPairs)
            throw std::length_error("Groebner critical-pair processing limit exceeded");

        if (chainCriterion(pair, result.basis, ring, zeroPairs)) {
            zeroPairs.insert(pairKey(pair.first, pair.second));
            continue;
        }

        MultivariateRationalPolynomial s = sPolynomial(
            result.basis[pair.first], result.basis[pair.second], ring);
        MultivariateRationalPolynomial remainder = normalForm(
            s, result.basis, ring, options.maximumReductionSteps);
        ++result.reductions;
        if (remainder.isZero()) {
            zeroPairs.insert(pairKey(pair.first, pair.second));
            continue;
        }

        checkPolynomialSize(remainder, options.maximumTermsPerPolynomial);
        remainder = monicPolynomial(remainder, ring);
        if (result.basis.size() >= options.maximumBasisSize)
            throw std::length_error("Groebner basis size limit exceeded");
        const std::size_t newIndex = result.basis.size();
        result.basis.push_back(std::move(remainder));

        const auto leading = result.basis.back().leadingTerm(ring);
        if (leading && leading->monomial.isOne()) {
            result.basis = {MultivariateRationalPolynomial{{PolynomialTerm{Monomial{}, one()}}}};
            return result;
        }

        for (std::size_t i = 0; i < newIndex; ++i)
            appendCriticalPair(
                pairs, result.basis, i, newIndex, ring, zeroPairs,
                options.maximumCriticalPairs);
    }

    result.basis = reducedGroebnerBasis(
        result.basis, ring, options.maximumReductionSteps);
    return result;
}

} // namespace mmcal::symbolic
