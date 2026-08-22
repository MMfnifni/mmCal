#pragma once

#include "polynomial.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace mmcal::symbolic {

struct PolynomialDivisionResult final {
    std::vector<MultivariateRationalPolynomial> quotients;
    MultivariateRationalPolynomial remainder;
};

struct GroebnerOptions final {
    // General算法ではあるが、対話CASとして資源爆発を無制限には許さない。
    std::size_t maximumBasisSize = 256;
    std::size_t maximumCriticalPairs = 50'000;
    std::size_t maximumReductionSteps = 500'000;
    std::size_t maximumTermsPerPolynomial = 100'000;
};

struct GroebnerComputation final {
    std::vector<MultivariateRationalPolynomial> basis;
    std::size_t criticalPairsProcessed = 0;
    std::size_t reductions = 0;
};

[[nodiscard]] PolynomialDivisionResult multivariateDivide(
    const MultivariateRationalPolynomial& dividend,
    std::span<const MultivariateRationalPolynomial> divisors,
    const PolynomialRing& ring,
    std::size_t maximumReductionSteps = 500'000);

[[nodiscard]] MultivariateRationalPolynomial normalForm(
    const MultivariateRationalPolynomial& polynomial,
    std::span<const MultivariateRationalPolynomial> basis,
    const PolynomialRing& ring,
    std::size_t maximumReductionSteps = 500'000);

[[nodiscard]] MultivariateRationalPolynomial sPolynomial(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs,
    const PolynomialRing& ring);

// exact Q[x1,...,xn] 上のreduced Gröbner basis。
// pair selectionはsugar degree、Buchberger product criterionとzero-pair chain criterionを使う。
[[nodiscard]] GroebnerComputation groebnerBasis(
    std::span<const MultivariateRationalPolynomial> generators,
    const PolynomialRing& ring,
    GroebnerOptions options = {});

[[nodiscard]] std::vector<MultivariateRationalPolynomial> reducedGroebnerBasis(
    std::span<const MultivariateRationalPolynomial> basis,
    const PolynomialRing& ring,
    std::size_t maximumReductionSteps = 500'000);

} // namespace mmcal::symbolic
