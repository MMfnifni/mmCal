#pragma once

#include "symbolic/differential_tower.hpp"
#include "symbolic/polynomial.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::symbolic::risch {

struct RischOptions final {
    std::size_t maximumTowerDepth = 16;
    std::size_t maximumHermiteDegree = 512;
    std::size_t maximumLrtDegree = 32;
    std::size_t maximumSubresultantSteps = 64;
    std::size_t maximumResultantMatrixEntries = 4096;
    std::size_t maximumBivariateCoefficients = 65536;
    std::size_t maximumIntermediateBits = 8192;
    std::size_t maximumResidueDegree = 96;
    std::size_t maximumDifferentialOperations = 32768;
    std::size_t maximumRdeDegree = 128;
    std::size_t maximumRdeMatrixEntries = 65536;
    std::size_t maximumRdeSteps = 4096;
    std::size_t maximumRdeResidueCandidates = 200000;
    std::size_t maximumRdeIntegerResidue = 64;
};

enum class RischFailure {
    None,
    InvalidRationalFunction,
    ZeroDenominator,
    NonSquareFreeDenominator,
    DegreeLimit,
    TowerDepthLimit,
    SubresultantStepLimit,
    MatrixSizeLimit,
    CoefficientCountLimit,
    IntermediateBitLimit,
    DifferentialOperationLimit,
    RdeMatrixSizeLimit,
    RdeStepLimit,
    ResidueSearchLimit,
    ExactDivisionFailed,
    NonInvertibleLeadingCoefficient,
    CertificateFailed,
    NonNormalDifferentialFactor,
    UnsupportedExtension,
    NoRationalSolution
};

template <class Value>
struct RischStageResult final {
    std::optional<Value> value;
    RischFailure failure = RischFailure::None;

    [[nodiscard]] explicit operator bool() const noexcept {
        return value.has_value() && failure == RischFailure::None;
    }
};

struct RationalFunction final {
    RationalPolynomial numerator;
    RationalPolynomial denominator;

    RationalFunction();
    RationalFunction(RationalPolynomial numerator, RationalPolynomial denominator);
};

struct RationalPolynomialDivision final {
    RationalPolynomial quotient;
    RationalPolynomial remainder;
};

// Risch各段が共有するQ[x]のexact演算。名前を限定し，通常の
// MultivariateRationalPolynomial演算とのoverload衝突を避ける。
[[nodiscard]] RationalPolynomial addRationalPolynomialsExact(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs);
[[nodiscard]] RationalPolynomial subtractRationalPolynomialsExact(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs);
[[nodiscard]] RationalPolynomial negateRationalPolynomialExact(
    const RationalPolynomial& value);
[[nodiscard]] RationalPolynomial scaleRationalPolynomialExact(
    const RationalPolynomial& value,
    const numeric::Rational& scale);
[[nodiscard]] RationalPolynomial multiplyRationalPolynomialsExact(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs);
[[nodiscard]] RationalPolynomial powerRationalPolynomialExact(
    RationalPolynomial base,
    std::size_t exponent);
[[nodiscard]] RationalPolynomial differentiateRationalPolynomialExact(
    const RationalPolynomial& value);
[[nodiscard]] RationalPolynomial integrateRationalPolynomialExact(
    const RationalPolynomial& value);
[[nodiscard]] RationalPolynomial monicRationalPolynomialExact(
    const RationalPolynomial& value);
[[nodiscard]] RationalPolynomial gcdRationalPolynomialsMonic(
    RationalPolynomial lhs,
    RationalPolynomial rhs);
[[nodiscard]] RationalPolynomialDivision divideRationalPolynomials(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator);
[[nodiscard]] std::optional<RationalPolynomial>
divideRationalPolynomialsExactly(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator);

// Differential extension の係数体 K=Q(x) が共有する exact 四則演算。
// すべて既約化し，等値判定はcross multiplicationで行う。
[[nodiscard]] RationalFunction canonicalizeRationalFunction(
    RationalFunction value);
[[nodiscard]] RationalFunction addRationalFunctionsExact(
    const RationalFunction& lhs,
    const RationalFunction& rhs);
[[nodiscard]] RationalFunction subtractRationalFunctionsExact(
    const RationalFunction& lhs,
    const RationalFunction& rhs);
[[nodiscard]] RationalFunction multiplyRationalFunctionsExact(
    const RationalFunction& lhs,
    const RationalFunction& rhs);
[[nodiscard]] std::optional<RationalFunction> divideRationalFunctionsExact(
    const RationalFunction& numerator,
    const RationalFunction& denominator);
[[nodiscard]] RationalFunction differentiateRationalFunctionExact(
    const RationalFunction& value);
[[nodiscard]] bool equivalentRationalFunctions(
    const RationalFunction& lhs,
    const RationalFunction& rhs);

struct HermitePowerReduction final {
    std::vector<std::pair<RationalPolynomial, std::size_t>> rationalTerms;
    RationalPolynomial squareFreeNumerator;
    std::size_t steps = 0;
    bool exactVerified = false;
};

struct RationalHermiteReduction final {
    // input = D(rationalPart) + squareFreePart。
    RationalFunction rationalPart;
    RationalFunction squareFreePart;
    bool exactVerified = false;
};

// sum_i coefficientInX(i)(z) x^i を表す Q[z][x] 多項式。
class BivariateRationalPolynomial final {
public:
    BivariateRationalPolynomial();
    explicit BivariateRationalPolynomial(
        std::vector<RationalPolynomial> coefficientsInX);

    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] std::size_t degreeInX() const noexcept;
    [[nodiscard]] const RationalPolynomial& coefficientInX(
        std::size_t exponent) const noexcept;
    [[nodiscard]] const std::vector<RationalPolynomial>& coefficientsInX() const noexcept;
    [[nodiscard]] bool operator==(const BivariateRationalPolynomial& rhs) const noexcept;

private:
    std::vector<RationalPolynomial> coefficientsInX_;
    void normalize();
};

struct AlgebraicResidueLogTerm final {
    // residuePolynomial(z)=0 の全根 a について
    // sum_a a log(logArgument(a,x)) を表す compact な RootSum certificate。
    RationalPolynomial residuePolynomial;
    BivariateRationalPolynomial logArgument;
    std::size_t poleMultiplicity = 0;
    bool exactVerified = false;
};

struct LrtResult final {
    RationalPolynomial residueResultant;
    std::vector<BivariateRationalPolynomial> subresultantSequence;
    std::vector<AlgebraicResidueLogTerm> logarithmicTerms;
    bool exactVerified = false;
};

struct RationalRischDecomposition final {
    RationalHermiteReduction hermite;
    LrtResult logarithmicPart;
};

enum class RischResultStatus {
    Elementary,
    ProvenNonElementary,
    ConditionsRequired,
    ResourceLimit,
    UnsupportedExtension
};

struct RischResult final {
    RischResultStatus status = RischResultStatus::UnsupportedExtension;
    RischFailure failure = RischFailure::None;
    std::optional<RationalRischDecomposition> rational;
    bool exactVerified = false;
};

// 現行積分器から移植した square-free power Hermite step。
[[nodiscard]] RischStageResult<HermitePowerReduction> hermiteReduceSquareFreePower(
    RationalPolynomial numerator,
    const RationalPolynomial& squareFreeFactor,
    std::size_t denominatorPower,
    const RationalPolynomial& inverseDerivativeModuloFactor,
    const RischOptions& options = {});

[[nodiscard]] RischStageResult<RationalHermiteReduction> hermiteReduceRationalFunction(
    const RationalFunction& input,
    const RischOptions& options = {});

[[nodiscard]] bool verifyHermiteReduction(
    const RationalFunction& input,
    const RationalHermiteReduction& reduction);

// Lazard-Rioboo-Trager: resultant と subresultant PRS から algebraic residue を
// compact に構成し、Q[z]/(residuePolynomial) 上の整除で証明する。
[[nodiscard]] RischStageResult<LrtResult> lazardRiobooTrager(
    const RationalFunction& properSquareFreePart,
    const RischOptions& options = {});

[[nodiscard]] bool verifyLrtResult(
    const RationalFunction& properSquareFreePart,
    const LrtResult& result,
    const RischOptions& options = {});

[[nodiscard]] RischResult integrateRationalRisch(
    const RationalFunction& input,
    const RischOptions& options = {});

} // namespace mmcal::symbolic::risch
