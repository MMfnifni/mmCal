#pragma once

#include "symbolic/risch_differential_reduction.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace mmcal::symbolic::risch {

struct RationalRdeSolution final {
    // Dy + coefficient*y = rightHandSide のQ(x)解。
    RationalFunction solution;
    RationalPolynomial weakNormalizer;
    RationalPolynomial normalDenominator;
    std::size_t polynomialDegreeBound = 0;
    std::size_t linearRank = 0;
    std::size_t steps = 0;
    bool exactVerified = false;
};

// 基礎微分体(Q(x),D=d/dx)の一階Risch微分方程式をexactに解く。
// 解なしを返すのは，weak normalization後のnormal-denominator条件が破れるか，
// 証明済み次数境界の有限線形系がinconsistentである場合だけである。
[[nodiscard]] RischStageResult<RationalRdeSolution>
solveRationalDifferentialEquation(
    const RationalFunction& coefficient,
    const RationalFunction& rightHandSide,
    const RischOptions& options = {});

[[nodiscard]] bool verifyRationalDifferentialEquation(
    const RationalFunction& coefficient,
    const RationalFunction& rightHandSide,
    const RationalRdeSolution& solution);

struct LimitedRationalIntegral final {
    // input = D(rationalPart) + constantMultiple*differentialCoefficient。
    RationalFunction rationalPart;
    numeric::Rational constantMultiple;
    bool exactVerified = false;
};

// Primitive extensionのleading coefficientに必要なlimited integration。
[[nodiscard]] RischStageResult<LimitedRationalIntegral>
limitedIntegrateRationalFunction(
    const RationalFunction& input,
    const RationalFunction& differentialCoefficient,
    const RischOptions& options = {});

struct PrimitivePolynomialReduction final {
    // input = D(polynomialPart) + lowerFieldRemainder。
    DifferentialPolynomial polynomialPart;
    RationalFunction lowerFieldRemainder;
    std::size_t steps = 0;
    bool exactVerified = false;
};

[[nodiscard]] RischStageResult<PrimitivePolynomialReduction>
reducePrimitivePolynomial(
    const DifferentialPolynomial& input,
    const RationalFunction& differentialCoefficient,
    const RischOptions& options = {});

struct DifferentialRationalFunction final {
    DifferentialPolynomial numerator;
    DifferentialPolynomial denominator;
};

struct PrimitiveLogarithmicTerm final {
    numeric::Rational coefficient;
    DifferentialPolynomial argument;
};

struct PrimitiveRationalReduction final {
    // input = D(polynomialPart + sum rationalPart)
    //       + sum coefficient*D(argument)/argument
    //       + lowerFieldRemainder + sum residualPart。
    DifferentialPolynomial polynomialPart;
    std::vector<DifferentialRationalFunction> rationalPart;
    std::vector<PrimitiveLogarithmicTerm> logarithmicPart;
    RationalFunction lowerFieldRemainder;
    std::vector<DifferentialRationalFunction> residualPart;
    std::size_t steps = 0;
    bool exactVerified = false;
};

// Primitive extension K(t)，Dt=differentialCoefficientにおけるK(t)の
// factor-free Hermite reduction。通常の多項式square-free decompositionと
// differential normalityを分離し，証明できないspecial residueは残差に保つ。
[[nodiscard]] RischStageResult<PrimitiveRationalReduction>
reducePrimitiveRationalFunction(
    const DifferentialRationalFunction& input,
    const RationalFunction& differentialCoefficient,
    const RischOptions& options = {});

[[nodiscard]] bool verifyPrimitiveRationalReduction(
    const DifferentialRationalFunction& input,
    const RationalFunction& differentialCoefficient,
    const PrimitiveRationalReduction& reduction,
    const RischOptions& options = {});

struct ExponentialLaurentTerm final {
    std::int64_t exponent = 0;
    RationalFunction coefficient;
};

struct ExponentialLaurentReduction final {
    // input = D(sum coefficient*t^exponent) + lowerFieldRemainder。
    std::vector<ExponentialLaurentTerm> laurentPart;
    RationalFunction lowerFieldRemainder;
    std::size_t steps = 0;
    bool exactVerified = false;
};

// Dt/t=differentialCoefficientであるexponential extensionのLaurent partを，
// 各次数の基礎体RDEへ分解してexactに還元する。
[[nodiscard]] RischStageResult<ExponentialLaurentReduction>
reduceExponentialLaurentPolynomial(
    const std::vector<ExponentialLaurentTerm>& input,
    const RationalFunction& differentialCoefficient,
    const RischOptions& options = {});

[[nodiscard]] bool verifyExponentialLaurentReduction(
    const std::vector<ExponentialLaurentTerm>& input,
    const RationalFunction& differentialCoefficient,
    const ExponentialLaurentReduction& reduction);

} // namespace mmcal::symbolic::risch
