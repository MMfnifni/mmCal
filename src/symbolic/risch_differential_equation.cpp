#include "risch_differential_equation.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "symbolic/polynomial.hpp"

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::symbolic::risch {
namespace {

using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational zero() {
    return Rational{BigInt{0}};
}

[[nodiscard]] Rational one() {
    return Rational{BigInt{1}};
}

[[nodiscard]] RationalPolynomial constantPolynomial(const Rational& value) {
    return RationalPolynomial{{value}};
}

[[nodiscard]] RationalPolynomial onePolynomial() {
    return constantPolynomial(one());
}

[[nodiscard]] RationalFunction constantFunction(const Rational& value) {
    return RationalFunction{constantPolynomial(value), onePolynomial()};
}

[[nodiscard]] bool rationalWithinBudget(
    const Rational& value,
    const RischOptions& options) {
    return value.numerator().bitLength() <= options.maximumIntermediateBits
        && value.denominator().bitLength() <= options.maximumIntermediateBits;
}

[[nodiscard]] RischFailure polynomialBudgetFailure(
    const RationalPolynomial& value,
    const RischOptions& options) {
    if (value.degree() > options.maximumRdeDegree)
        return RischFailure::DegreeLimit;
    for (const Rational& coefficient : value.coefficients())
        if (!rationalWithinBudget(coefficient, options))
            return RischFailure::IntermediateBitLimit;
    return RischFailure::None;
}

[[nodiscard]] RischFailure functionBudgetFailure(
    const RationalFunction& value,
    const RischOptions& options) {
    if (value.denominator.isZero())
        return RischFailure::ZeroDenominator;
    const RischFailure numeratorFailure = polynomialBudgetFailure(
        value.numerator, options);
    if (numeratorFailure != RischFailure::None)
        return numeratorFailure;
    return polynomialBudgetFailure(value.denominator, options);
}

class RdeBudget final {
public:
    explicit RdeBudget(const RischOptions& options)
        : options_(options) {}

    [[nodiscard]] bool charge(std::size_t amount = 1) {
        if (failure_ != RischFailure::None)
            return false;
        if (amount > options_.maximumRdeSteps - std::min(
                steps_, options_.maximumRdeSteps)) {
            failure_ = RischFailure::RdeStepLimit;
            return false;
        }
        steps_ += amount;
        return true;
    }

    [[nodiscard]] bool check(const Rational& value) {
        if (!charge())
            return false;
        if (!rationalWithinBudget(value, options_)) {
            failure_ = RischFailure::IntermediateBitLimit;
            return false;
        }
        return true;
    }

    [[nodiscard]] bool check(const RationalPolynomial& value) {
        if (!charge())
            return false;
        const RischFailure failure = polynomialBudgetFailure(value, options_);
        if (failure != RischFailure::None) {
            failure_ = failure;
            return false;
        }
        return true;
    }

    [[nodiscard]] bool check(const RationalFunction& value) {
        if (!charge())
            return false;
        const RischFailure failure = functionBudgetFailure(value, options_);
        if (failure != RischFailure::None) {
            failure_ = failure;
            return false;
        }
        return true;
    }

    [[nodiscard]] RischFailure failure() const noexcept {
        return failure_;
    }

    [[nodiscard]] std::size_t steps() const noexcept {
        return steps_;
    }

private:
    const RischOptions& options_;
    std::size_t steps_ = 0;
    RischFailure failure_ = RischFailure::None;
};

[[nodiscard]] std::optional<std::uint64_t> positiveIntegerValue(
    const Rational& value) {
    if (!value.isInteger() || !value.numerator().isPositive())
        return std::nullopt;
    const std::string text = value.numerator().toString();
    std::uint64_t result = 0;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), result);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
        return std::nullopt;
    return result;
}

[[nodiscard]] RationalPolynomial evaluateAtResidue(
    const BivariateRationalPolynomial& value,
    const Rational& residue) {
    std::vector<Rational> coefficients(value.degreeInX() + 1, zero());
    for (std::size_t exponent = 0; exponent <= value.degreeInX(); ++exponent)
        coefficients[exponent] = evaluatePolynomial(
            value.coefficientInX(exponent), residue);
    return RationalPolynomial{std::move(coefficients)};
}

struct WeakNormalization final {
    RationalFunction coefficient;
    RationalFunction rightHandSide;
    RationalPolynomial normalizer = onePolynomial();
};

[[nodiscard]] RischStageResult<WeakNormalization> weakNormalize(
    const RationalFunction& coefficient,
    const RationalFunction& rightHandSide,
    const RischOptions& options,
    RdeBudget& budget) {
    auto hermite = hermiteReduceRationalFunction(coefficient, options);
    if (!hermite)
        return {std::nullopt, hermite.failure};
    auto residues = lazardRiobooTrager(
        hermite.value->squareFreePart, options);
    if (!residues)
        return {std::nullopt, residues.failure};

    RationalPolynomial normalizer = onePolynomial();
    for (const AlgebraicResidueLogTerm& term
         : residues.value->logarithmicTerms) {
        RationalPolynomial remaining = term.residuePolynomial;
        while (remaining.degree() != 0) {
            if (!budget.charge())
                return {std::nullopt, budget.failure()};
            const RationalRootSearchResult root = findRationalRoot(
                remaining,
                RationalRootSearchOptions{
                    options.maximumRdeResidueCandidates,
                    options.maximumRdeResidueCandidates});
            if (!root.root) {
                if (!root.complete)
                    return {std::nullopt, RischFailure::ResidueSearchLimit};
                break;
            }
            const auto quotient = divideByLinearFactor(remaining, *root.root);
            if (!quotient)
                return {std::nullopt, RischFailure::CertificateFailed};
            if (root.root->isInteger()
                && root.root->numerator().isPositive()) {
                const auto integer = positiveIntegerValue(*root.root);
                if (!integer)
                    return {std::nullopt, RischFailure::DegreeLimit};
                if (*integer > options.maximumRdeIntegerResidue)
                    return {std::nullopt, RischFailure::DegreeLimit};
                RationalPolynomial factor = evaluateAtResidue(
                    term.logArgument, *root.root);
                if (factor.isZero() || factor.degree() == 0)
                    return {std::nullopt, RischFailure::CertificateFailed};
                factor = powerRationalPolynomialExact(
                    std::move(factor), static_cast<std::size_t>(*integer));
                normalizer = multiplyRationalPolynomialsExact(
                    normalizer, factor);
                if (!budget.check(normalizer))
                    return {std::nullopt, budget.failure()};
            }
            remaining = *quotient;
        }
    }

    const RationalFunction logarithmicDerivative{
        differentiateRationalPolynomialExact(normalizer), normalizer};
    WeakNormalization result{
        subtractRationalFunctionsExact(coefficient, logarithmicDerivative),
        multiplyRationalFunctionsExact(
            rightHandSide,
            RationalFunction{normalizer, onePolynomial()}),
        std::move(normalizer)};
    if (!budget.check(result.coefficient)
        || !budget.check(result.rightHandSide))
        return {std::nullopt, budget.failure()};
    return {std::move(result), RischFailure::None};
}

struct PolynomialRde final {
    RationalPolynomial derivativeCoefficient;
    RationalPolynomial valueCoefficient;
    RationalPolynomial rightHandSide;
    RationalPolynomial denominator;
};

[[nodiscard]] RischStageResult<PolynomialRde> normalDenominator(
    const RationalFunction& coefficient,
    const RationalFunction& rightHandSide,
    RdeBudget& budget) {
    const RationalFunction normalizedCoefficient =
        canonicalizeRationalFunction(coefficient);
    const RationalFunction normalizedRightHandSide =
        canonicalizeRationalFunction(rightHandSide);
    if (normalizedCoefficient.denominator.isZero()
        || normalizedRightHandSide.denominator.isZero())
        return {std::nullopt, RischFailure::ZeroDenominator};

    const RationalPolynomial sharedDenominator = gcdRationalPolynomialsMonic(
        normalizedCoefficient.denominator,
        normalizedRightHandSide.denominator);
    const RationalPolynomial repeatedRight = gcdRationalPolynomialsMonic(
        normalizedRightHandSide.denominator,
        differentiateRationalPolynomialExact(
            normalizedRightHandSide.denominator));
    const RationalPolynomial repeatedShared = gcdRationalPolynomialsMonic(
        sharedDenominator,
        differentiateRationalPolynomialExact(sharedDenominator));
    auto denominator = divideRationalPolynomialsExactly(
        repeatedRight, repeatedShared);
    if (!denominator)
        return {std::nullopt, RischFailure::ExactDivisionFailed};

    PolynomialRde result;
    result.denominator = monicRationalPolynomialExact(*denominator);
    result.derivativeCoefficient = multiplyRationalPolynomialsExact(
        normalizedCoefficient.denominator, result.denominator);
    result.valueCoefficient = subtractRationalPolynomialsExact(
        multiplyRationalPolynomialsExact(
            normalizedCoefficient.numerator, result.denominator),
        multiplyRationalPolynomialsExact(
            normalizedCoefficient.denominator,
            differentiateRationalPolynomialExact(result.denominator)));
    const RationalPolynomial rightNumerator =
        multiplyRationalPolynomialsExact(
            multiplyRationalPolynomialsExact(
                normalizedCoefficient.denominator,
                powerRationalPolynomialExact(result.denominator, 2)),
            normalizedRightHandSide.numerator);
    auto polynomialRight = divideRationalPolynomialsExactly(
        rightNumerator, normalizedRightHandSide.denominator);
    if (!polynomialRight)
        return {std::nullopt, RischFailure::NoRationalSolution};
    result.rightHandSide = std::move(*polynomialRight);
    for (const RationalPolynomial* polynomial : {
             &result.derivativeCoefficient,
             &result.valueCoefficient,
             &result.rightHandSide,
             &result.denominator})
        if (!budget.check(*polynomial))
            return {std::nullopt, budget.failure()};
    return {std::move(result), RischFailure::None};
}

[[nodiscard]] std::ptrdiff_t degreeOrMinusOne(
    const RationalPolynomial& value) {
    return value.isZero()
        ? std::ptrdiff_t{-1}
        : static_cast<std::ptrdiff_t>(value.degree());
}

[[nodiscard]] RischStageResult<std::size_t> polynomialDegreeBound(
    const PolynomialRde& equation,
    const RischOptions& options) {
    if (equation.derivativeCoefficient.isZero())
        return {std::nullopt, RischFailure::InvalidRationalFunction};
    if (equation.rightHandSide.isZero())
        return {std::size_t{0}, RischFailure::None};
    const std::ptrdiff_t derivativeDegree =
        degreeOrMinusOne(equation.derivativeCoefficient) - 1;
    const std::ptrdiff_t valueDegree =
        degreeOrMinusOne(equation.valueCoefficient);
    const std::ptrdiff_t rightDegree =
        degreeOrMinusOne(equation.rightHandSide);
    const std::ptrdiff_t operatorDegree = std::max(
        derivativeDegree, valueDegree);
    std::ptrdiff_t bound = std::max<std::ptrdiff_t>(
        0, rightDegree - operatorDegree);

    if (!equation.valueCoefficient.isZero()
        && valueDegree == derivativeDegree) {
        const Rational resonance = -equation.valueCoefficient.coefficient(
            equation.valueCoefficient.degree())
            / equation.derivativeCoefficient.coefficient(
                equation.derivativeCoefficient.degree());
        if (resonance.isInteger() && !resonance.numerator().isNegative()) {
            const std::string text = resonance.numerator().toString();
            std::uint64_t parsedValue = 0;
            const auto parsed = std::from_chars(
                text.data(), text.data() + text.size(), parsedValue);
            if (parsed.ec != std::errc{}
                || parsed.ptr != text.data() + text.size()
                || parsedValue > options.maximumRdeDegree)
                return {std::nullopt, RischFailure::DegreeLimit};
            const std::size_t exceptional =
                static_cast<std::size_t>(parsedValue);
            bound = std::max(
                bound, static_cast<std::ptrdiff_t>(exceptional));
        }
    }
    if (bound > static_cast<std::ptrdiff_t>(options.maximumRdeDegree))
        return {std::nullopt, RischFailure::DegreeLimit};
    return {static_cast<std::size_t>(bound), RischFailure::None};
}

struct LinearPolynomialSolution final {
    RationalPolynomial polynomial;
    std::size_t rank = 0;
};

[[nodiscard]] RischStageResult<LinearPolynomialSolution>
solvePolynomialEquation(
    const PolynomialRde& equation,
    std::size_t degreeBound,
    const RischOptions& options,
    RdeBudget& budget) {
    const std::size_t columns = degreeBound + 1;
    const std::size_t derivativeOutput = equation.derivativeCoefficient.degree()
        + (degreeBound == 0 ? 0 : degreeBound - 1);
    const std::size_t valueOutput = equation.valueCoefficient.isZero()
        ? 0
        : equation.valueCoefficient.degree() + degreeBound;
    const std::size_t rows = std::max({
        derivativeOutput,
        valueOutput,
        equation.rightHandSide.degree()}) + 1;
    if (columns > std::numeric_limits<std::size_t>::max() / (rows + 1)
        || rows * (columns + 1) > options.maximumRdeMatrixEntries)
        return {std::nullopt, RischFailure::RdeMatrixSizeLimit};

    std::vector<std::vector<Rational>> matrix(
        rows, std::vector<Rational>(columns + 1, zero()));
    for (std::size_t unknown = 0; unknown < columns; ++unknown) {
        if (unknown != 0) {
            const Rational multiplier{BigInt::fromUnsigned(unknown)};
            for (std::size_t i = 0;
                 i <= equation.derivativeCoefficient.degree(); ++i)
                matrix[i + unknown - 1][unknown] +=
                    equation.derivativeCoefficient.coefficient(i) * multiplier;
        }
        if (!equation.valueCoefficient.isZero())
            for (std::size_t i = 0;
                 i <= equation.valueCoefficient.degree(); ++i)
                matrix[i + unknown][unknown] +=
                    equation.valueCoefficient.coefficient(i);
    }
    for (std::size_t row = 0; row <= equation.rightHandSide.degree(); ++row)
        matrix[row][columns] = equation.rightHandSide.coefficient(row);
    if (!budget.charge(rows * (columns + 1)))
        return {std::nullopt, budget.failure()};

    std::vector<std::size_t> pivotColumns;
    std::size_t pivotRow = 0;
    for (std::size_t column = 0;
         column < columns && pivotRow < rows; ++column) {
        std::size_t selected = pivotRow;
        while (selected < rows && matrix[selected][column].isZero())
            ++selected;
        if (selected == rows)
            continue;
        if (selected != pivotRow)
            std::swap(matrix[selected], matrix[pivotRow]);
        const Rational pivot = matrix[pivotRow][column];
        for (std::size_t entry = column; entry <= columns; ++entry) {
            matrix[pivotRow][entry] /= pivot;
            if (!budget.check(matrix[pivotRow][entry]))
                return {std::nullopt, budget.failure()};
        }
        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivotRow || matrix[row][column].isZero())
                continue;
            const Rational scale = matrix[row][column];
            for (std::size_t entry = column; entry <= columns; ++entry) {
                matrix[row][entry] -= scale * matrix[pivotRow][entry];
                if (!budget.check(matrix[row][entry]))
                    return {std::nullopt, budget.failure()};
            }
        }
        pivotColumns.push_back(column);
        ++pivotRow;
    }

    for (std::size_t row = 0; row < rows; ++row) {
        bool zeroLeft = true;
        for (std::size_t column = 0; column < columns; ++column)
            zeroLeft = zeroLeft && matrix[row][column].isZero();
        if (zeroLeft && !matrix[row][columns].isZero())
            return {std::nullopt, RischFailure::NoRationalSolution};
    }

    std::vector<Rational> coefficients(columns, zero());
    for (std::size_t row = 0; row < pivotColumns.size(); ++row)
        coefficients[pivotColumns[row]] = matrix[row][columns];
    RationalPolynomial solution{std::move(coefficients)};
    const RationalPolynomial reconstructed = addRationalPolynomialsExact(
        multiplyRationalPolynomialsExact(
            equation.derivativeCoefficient,
            differentiateRationalPolynomialExact(solution)),
        multiplyRationalPolynomialsExact(
            equation.valueCoefficient, solution));
    if (reconstructed.coefficients()
        != equation.rightHandSide.coefficients())
        return {std::nullopt, RischFailure::CertificateFailed};
    return {
        LinearPolynomialSolution{std::move(solution), pivotColumns.size()},
        RischFailure::None};
}

[[nodiscard]] bool validFunction(const RationalFunction& value) {
    return !value.denominator.isZero();
}

[[nodiscard]] std::optional<Rational> constantValue(
    const RationalFunction& value) {
    const RationalFunction normalized = canonicalizeRationalFunction(value);
    if (normalized.denominator.isZero()
        || normalized.numerator.degree() != 0
        || normalized.denominator.degree() != 0)
        return std::nullopt;
    return normalized.numerator.coefficient(0)
        / normalized.denominator.coefficient(0);
}

[[nodiscard]] DifferentialPolynomial addDifferentialPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs) {
    const std::size_t count = std::max(
        lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<RationalFunction> coefficients(count);
    for (std::size_t i = 0; i < count; ++i)
        coefficients[i] = addRationalFunctionsExact(
            lhs.coefficient(i), rhs.coefficient(i));
    return DifferentialPolynomial{std::move(coefficients)};
}

[[nodiscard]] DifferentialPolynomial subtractDifferentialPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs) {
    const std::size_t count = std::max(
        lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<RationalFunction> coefficients(count);
    for (std::size_t i = 0; i < count; ++i)
        coefficients[i] = subtractRationalFunctionsExact(
            lhs.coefficient(i), rhs.coefficient(i));
    return DifferentialPolynomial{std::move(coefficients)};
}

[[nodiscard]] bool differentialPolynomialWithinBudget(
    const DifferentialPolynomial& value,
    const RischOptions& options) {
    if (value.degree() > options.maximumRdeDegree)
        return false;
    for (const RationalFunction& coefficient : value.coefficients())
        if (functionBudgetFailure(coefficient, options) != RischFailure::None)
            return false;
    return true;
}

using LaurentMap = std::map<std::int64_t, RationalFunction>;

[[nodiscard]] LaurentMap normalizeLaurent(
    const std::vector<ExponentialLaurentTerm>& terms) {
    LaurentMap result;
    for (const ExponentialLaurentTerm& term : terms) {
        auto [iterator, inserted] = result.try_emplace(
            term.exponent, canonicalizeRationalFunction(term.coefficient));
        if (!inserted)
            iterator->second = addRationalFunctionsExact(
                iterator->second, term.coefficient);
    }
    for (auto iterator = result.begin(); iterator != result.end();) {
        if (iterator->second.numerator.isZero())
            iterator = result.erase(iterator);
        else
            ++iterator;
    }
    return result;
}

[[nodiscard]] std::vector<ExponentialLaurentTerm> laurentTerms(
    LaurentMap value) {
    std::vector<ExponentialLaurentTerm> result;
    result.reserve(value.size());
    for (auto& [exponent, coefficient] : value)
        result.push_back({exponent, std::move(coefficient)});
    return result;
}

} // namespace

RischStageResult<RationalRdeSolution> solveRationalDifferentialEquation(
    const RationalFunction& coefficient,
    const RationalFunction& rightHandSide,
    const RischOptions& options) {
    if (!validFunction(coefficient) || !validFunction(rightHandSide))
        return {std::nullopt, RischFailure::ZeroDenominator};
    const RationalFunction normalizedCoefficient =
        canonicalizeRationalFunction(coefficient);
    const RationalFunction normalizedRightHandSide =
        canonicalizeRationalFunction(rightHandSide);
    if (normalizedRightHandSide.numerator.isZero()) {
        RationalRdeSolution zeroSolution;
        zeroSolution.weakNormalizer = onePolynomial();
        zeroSolution.normalDenominator = onePolynomial();
        zeroSolution.exactVerified = true;
        return {std::move(zeroSolution), RischFailure::None};
    }
    const RischFailure coefficientFailure = functionBudgetFailure(
        normalizedCoefficient, options);
    if (coefficientFailure != RischFailure::None)
        return {std::nullopt, coefficientFailure};
    const RischFailure rightFailure = functionBudgetFailure(
        normalizedRightHandSide, options);
    if (rightFailure != RischFailure::None)
        return {std::nullopt, rightFailure};

    RdeBudget budget{options};
    auto weak = weakNormalize(
        normalizedCoefficient, normalizedRightHandSide, options, budget);
    if (!weak)
        return {std::nullopt, weak.failure};
    auto polynomialEquation = normalDenominator(
        weak.value->coefficient, weak.value->rightHandSide, budget);
    if (!polynomialEquation)
        return {std::nullopt, polynomialEquation.failure};
    auto bound = polynomialDegreeBound(*polynomialEquation.value, options);
    if (!bound)
        return {std::nullopt, bound.failure};
    auto polynomialSolution = solvePolynomialEquation(
        *polynomialEquation.value, *bound.value, options, budget);
    if (!polynomialSolution)
        return {std::nullopt, polynomialSolution.failure};

    RationalFunction solution{
        polynomialSolution.value->polynomial,
        multiplyRationalPolynomialsExact(
            polynomialEquation.value->denominator,
            weak.value->normalizer)};
    solution = canonicalizeRationalFunction(std::move(solution));
    if (!budget.check(solution))
        return {std::nullopt, budget.failure()};
    RationalRdeSolution result{
        std::move(solution),
        weak.value->normalizer,
        polynomialEquation.value->denominator,
        *bound.value,
        polynomialSolution.value->rank,
        budget.steps(),
        false};
    result.exactVerified = verifyRationalDifferentialEquation(
        normalizedCoefficient, normalizedRightHandSide, result);
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

bool verifyRationalDifferentialEquation(
    const RationalFunction& coefficient,
    const RationalFunction& rightHandSide,
    const RationalRdeSolution& solution) {
    if (!validFunction(coefficient) || !validFunction(rightHandSide)
        || !validFunction(solution.solution))
        return false;
    const RationalFunction reconstructed = addRationalFunctionsExact(
        differentiateRationalFunctionExact(solution.solution),
        multiplyRationalFunctionsExact(coefficient, solution.solution));
    return equivalentRationalFunctions(reconstructed, rightHandSide);
}

RischStageResult<LimitedRationalIntegral>
limitedIntegrateRationalFunction(
    const RationalFunction& input,
    const RationalFunction& differentialCoefficient,
    const RischOptions& options) {
    if (!validFunction(input) || !validFunction(differentialCoefficient))
        return {std::nullopt, RischFailure::ZeroDenominator};
    auto inputHermite = hermiteReduceRationalFunction(input, options);
    if (!inputHermite)
        return {std::nullopt, inputHermite.failure};
    auto differentialHermite = hermiteReduceRationalFunction(
        differentialCoefficient, options);
    if (!differentialHermite)
        return {std::nullopt, differentialHermite.failure};

    Rational multiplier = zero();
    const RationalFunction& inputRemainder =
        inputHermite.value->squareFreePart;
    const RationalFunction& differentialRemainder =
        differentialHermite.value->squareFreePart;
    if (!inputRemainder.numerator.isZero()) {
        if (differentialRemainder.numerator.isZero())
            return {std::nullopt, RischFailure::NoRationalSolution};
        const auto quotient = divideRationalFunctionsExact(
            inputRemainder, differentialRemainder);
        if (!quotient)
            return {std::nullopt, RischFailure::ExactDivisionFailed};
        const auto constant = constantValue(*quotient);
        if (!constant)
            return {std::nullopt, RischFailure::NoRationalSolution};
        multiplier = *constant;
    }

    const RationalFunction rationalPart = subtractRationalFunctionsExact(
        inputHermite.value->rationalPart,
        multiplyRationalFunctionsExact(
            constantFunction(multiplier),
            differentialHermite.value->rationalPart));
    LimitedRationalIntegral result{
        rationalPart, multiplier, false};
    const RationalFunction reconstructed = addRationalFunctionsExact(
        differentiateRationalFunctionExact(result.rationalPart),
        multiplyRationalFunctionsExact(
            constantFunction(result.constantMultiple),
            differentialCoefficient));
    result.exactVerified = equivalentRationalFunctions(reconstructed, input);
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

RischStageResult<PrimitivePolynomialReduction> reducePrimitivePolynomial(
    const DifferentialPolynomial& input,
    const RationalFunction& differentialCoefficient,
    const RischOptions& options) {
    if (!validFunction(differentialCoefficient))
        return {std::nullopt, RischFailure::ZeroDenominator};
    if (input.degree() > options.maximumRdeDegree)
        return {std::nullopt, RischFailure::DegreeLimit};
    const DifferentialDerivation derivation{
        DifferentialExtensionKind::Primitive, differentialCoefficient};
    DifferentialPolynomial remainder = input;
    DifferentialPolynomial antiderivative;
    std::size_t steps = 0;
    while (!remainder.isZero() && remainder.degree() != 0) {
        if (++steps > options.maximumRdeSteps)
            return {std::nullopt, RischFailure::RdeStepLimit};
        const std::size_t degree = remainder.degree();
        auto leadingIntegral = limitedIntegrateRationalFunction(
            remainder.coefficient(degree), differentialCoefficient, options);
        if (!leadingIntegral)
            return {std::nullopt, leadingIntegral.failure};
        std::vector<RationalFunction> coefficients(degree + 2);
        coefficients[degree] = leadingIntegral.value->rationalPart;
        coefficients[degree + 1] = constantFunction(
            leadingIntegral.value->constantMultiple
            / Rational{BigInt::fromUnsigned(degree + 1)});
        DifferentialPolynomial correction{std::move(coefficients)};
        auto derivative = differentiateDifferentialPolynomial(
            correction, derivation, options);
        if (!derivative)
            return {std::nullopt, derivative.failure};
        remainder = subtractDifferentialPolynomials(
            remainder, *derivative.value);
        antiderivative = addDifferentialPolynomials(
            antiderivative, correction);
        if (!differentialPolynomialWithinBudget(remainder, options)
            || !differentialPolynomialWithinBudget(antiderivative, options))
            return {std::nullopt, RischFailure::IntermediateBitLimit};
        if (!remainder.isZero() && remainder.degree() >= degree)
            return {std::nullopt, RischFailure::CertificateFailed};
    }

    PrimitivePolynomialReduction result{
        std::move(antiderivative),
        remainder.coefficient(0),
        steps,
        false};
    auto derivative = differentiateDifferentialPolynomial(
        result.polynomialPart, derivation, options);
    if (!derivative)
        return {std::nullopt, derivative.failure};
    const DifferentialPolynomial reconstructed = addDifferentialPolynomials(
        *derivative.value,
        DifferentialPolynomial{{result.lowerFieldRemainder}});
    result.exactVerified = equivalentDifferentialPolynomials(
        reconstructed, input);
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

RischStageResult<ExponentialLaurentReduction>
reduceExponentialLaurentPolynomial(
    const std::vector<ExponentialLaurentTerm>& input,
    const RationalFunction& differentialCoefficient,
    const RischOptions& options) {
    if (!validFunction(differentialCoefficient))
        return {std::nullopt, RischFailure::ZeroDenominator};
    const LaurentMap normalized = normalizeLaurent(input);
    LaurentMap antiderivative;
    RationalFunction lowerFieldRemainder;
    std::size_t steps = 0;
    for (const auto& [exponent, coefficient] : normalized) {
        if (exponent == 0) {
            lowerFieldRemainder = addRationalFunctionsExact(
                lowerFieldRemainder, coefficient);
            continue;
        }
        if (++steps > options.maximumRdeSteps)
            return {std::nullopt, RischFailure::RdeStepLimit};
        const RationalFunction rdeCoefficient =
            multiplyRationalFunctionsExact(
                constantFunction(Rational{BigInt{exponent}}),
                differentialCoefficient);
        auto solved = solveRationalDifferentialEquation(
            rdeCoefficient, coefficient, options);
        if (!solved)
            return {std::nullopt, solved.failure};
        antiderivative.emplace(exponent, solved.value->solution);
    }
    ExponentialLaurentReduction result{
        laurentTerms(std::move(antiderivative)),
        std::move(lowerFieldRemainder),
        steps,
        false};
    result.exactVerified = verifyExponentialLaurentReduction(
        input, differentialCoefficient, result);
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

bool verifyExponentialLaurentReduction(
    const std::vector<ExponentialLaurentTerm>& input,
    const RationalFunction& differentialCoefficient,
    const ExponentialLaurentReduction& reduction) {
    if (!validFunction(differentialCoefficient)
        || !validFunction(reduction.lowerFieldRemainder))
        return false;
    const LaurentMap expected = normalizeLaurent(input);
    LaurentMap reconstructed;
    for (const ExponentialLaurentTerm& term : reduction.laurentPart) {
        const RationalFunction differentiated = addRationalFunctionsExact(
            differentiateRationalFunctionExact(term.coefficient),
            multiplyRationalFunctionsExact(
                constantFunction(Rational{BigInt{term.exponent}}),
                multiplyRationalFunctionsExact(
                    differentialCoefficient, term.coefficient)));
        auto [iterator, inserted] = reconstructed.try_emplace(
            term.exponent, differentiated);
        if (!inserted)
            iterator->second = addRationalFunctionsExact(
                iterator->second, differentiated);
    }
    reconstructed[0] = addRationalFunctionsExact(
        reconstructed[0], reduction.lowerFieldRemainder);
    for (auto iterator = reconstructed.begin(); iterator != reconstructed.end();) {
        if (iterator->second.numerator.isZero())
            iterator = reconstructed.erase(iterator);
        else
            ++iterator;
    }
    if (reconstructed.size() != expected.size())
        return false;
    auto lhs = reconstructed.begin();
    auto rhs = expected.begin();
    for (; lhs != reconstructed.end(); ++lhs, ++rhs)
        if (lhs->first != rhs->first
            || !equivalentRationalFunctions(lhs->second, rhs->second))
            return false;
    return true;
}

} // namespace mmcal::symbolic::risch
