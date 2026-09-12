#include "risch_differential_equation.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::symbolic::risch {
namespace {

using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational one() {
    return Rational{BigInt{1}};
}

[[nodiscard]] RationalFunction constantFunction(const Rational& value) {
    return RationalFunction{
        RationalPolynomial{{value}},
        RationalPolynomial{{one()}}};
}

[[nodiscard]] RationalFunction oneFunction() {
    return constantFunction(one());
}

[[nodiscard]] bool functionWithinBudget(
    const RationalFunction& value,
    const RischOptions& options) {
    if (value.denominator.isZero()
        || value.numerator.degree() > options.maximumRdeDegree
        || value.denominator.degree() > options.maximumRdeDegree)
        return false;
    for (const RationalPolynomial* polynomial
         : {&value.numerator, &value.denominator})
        for (const Rational& coefficient : polynomial->coefficients())
            if (coefficient.numerator().bitLength() > options.maximumIntermediateBits
                || coefficient.denominator().bitLength()
                    > options.maximumIntermediateBits)
                return false;
    return true;
}

class PrimitiveBudget final {
public:
    explicit PrimitiveBudget(const RischOptions& options)
        : options_(options) {}

    bool check(
        const DifferentialPolynomial& value,
        std::size_t charge = 1) {
        if (!consume(charge))
            return false;
        if (value.degree() > options_.maximumHermiteDegree) {
            failure_ = RischFailure::DegreeLimit;
            return false;
        }
        for (const RationalFunction& coefficient : value.coefficients())
            if (!functionWithinBudget(coefficient, options_)) {
                failure_ = RischFailure::IntermediateBitLimit;
                return false;
            }
        return true;
    }

    [[nodiscard]] bool consume(std::size_t amount = 1) {
        if (failure_ != RischFailure::None)
            return false;
        if (amount > options_.maximumDifferentialOperations
                - std::min(operations_, options_.maximumDifferentialOperations)) {
            failure_ = RischFailure::DifferentialOperationLimit;
            return false;
        }
        operations_ += amount;
        return true;
    }

    void fail(RischFailure failure) {
        if (failure_ == RischFailure::None)
            failure_ = failure;
    }

    [[nodiscard]] RischFailure failure() const noexcept {
        return failure_;
    }

    [[nodiscard]] std::size_t operations() const noexcept {
        return operations_;
    }

private:
    const RischOptions& options_;
    std::size_t operations_ = 0;
    RischFailure failure_ = RischFailure::None;
};

[[nodiscard]] DifferentialPolynomial constantPolynomial(
    const RationalFunction& value) {
    return DifferentialPolynomial{{value}};
}

[[nodiscard]] DifferentialPolynomial onePolynomial() {
    return constantPolynomial(oneFunction());
}

[[nodiscard]] bool isOnePolynomial(const DifferentialPolynomial& value) {
    return value.degree() == 0
        && equivalentRationalFunctions(value.coefficient(0), oneFunction());
}

[[nodiscard]] DifferentialPolynomial addPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs,
    PrimitiveBudget& budget) {
    const std::size_t count = std::max(
        lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<RationalFunction> coefficients(count);
    for (std::size_t i = 0; i < count; ++i)
        coefficients[i] = addRationalFunctionsExact(
            lhs.coefficient(i), rhs.coefficient(i));
    DifferentialPolynomial result{std::move(coefficients)};
    budget.check(result, count);
    return result;
}

[[nodiscard]] DifferentialPolynomial subtractPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs,
    PrimitiveBudget& budget) {
    const std::size_t count = std::max(
        lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<RationalFunction> coefficients(count);
    for (std::size_t i = 0; i < count; ++i)
        coefficients[i] = subtractRationalFunctionsExact(
            lhs.coefficient(i), rhs.coefficient(i));
    DifferentialPolynomial result{std::move(coefficients)};
    budget.check(result, count);
    return result;
}

[[nodiscard]] DifferentialPolynomial multiplyPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs,
    PrimitiveBudget& budget) {
    if (lhs.isZero() || rhs.isZero())
        return {};
    const std::size_t coefficientCount = lhs.degree() + rhs.degree() + 1;
    std::vector<RationalFunction> coefficients(coefficientCount);
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        for (std::size_t j = 0; j <= rhs.degree(); ++j)
            coefficients[i + j] = addRationalFunctionsExact(
                coefficients[i + j],
                multiplyRationalFunctionsExact(
                    lhs.coefficient(i), rhs.coefficient(j)));
    DifferentialPolynomial result{std::move(coefficients)};
    budget.check(result, (lhs.degree() + 1) * (rhs.degree() + 1));
    return result;
}

[[nodiscard]] DifferentialPolynomial scalePolynomial(
    const DifferentialPolynomial& value,
    const RationalFunction& scale,
    PrimitiveBudget& budget) {
    std::vector<RationalFunction> coefficients;
    coefficients.reserve(value.coefficients().size());
    for (const RationalFunction& coefficient : value.coefficients())
        coefficients.push_back(multiplyRationalFunctionsExact(
            coefficient, scale));
    DifferentialPolynomial result{std::move(coefficients)};
    budget.check(result, value.coefficients().size());
    return result;
}

struct PolynomialDivision final {
    DifferentialPolynomial quotient;
    DifferentialPolynomial remainder;
    bool exact = false;
};

[[nodiscard]] PolynomialDivision dividePolynomials(
    const DifferentialPolynomial& numerator,
    const DifferentialPolynomial& denominator,
    PrimitiveBudget& budget) {
    if (denominator.isZero()) {
        budget.fail(RischFailure::ZeroDenominator);
        return {};
    }
    DifferentialPolynomial remainder = numerator;
    std::vector<RationalFunction> quotient(
        numerator.degree() >= denominator.degree()
            ? numerator.degree() - denominator.degree() + 1
            : 1);
    while (!remainder.isZero()
        && remainder.degree() >= denominator.degree()
        && budget.failure() == RischFailure::None) {
        const std::size_t shift = remainder.degree() - denominator.degree();
        const auto scale = divideRationalFunctionsExact(
            remainder.coefficient(remainder.degree()),
            denominator.coefficient(denominator.degree()));
        if (!scale) {
            budget.fail(RischFailure::ExactDivisionFailed);
            return {};
        }
        quotient[shift] = addRationalFunctionsExact(quotient[shift], *scale);
        std::vector<RationalFunction> product(shift);
        product.reserve(shift + denominator.coefficients().size());
        for (const RationalFunction& coefficient : denominator.coefficients())
            product.push_back(multiplyRationalFunctionsExact(
                coefficient, *scale));
        remainder = subtractPolynomials(
            remainder, DifferentialPolynomial{std::move(product)}, budget);
        if (!budget.consume())
            return {};
    }
    DifferentialPolynomial resultQuotient{std::move(quotient)};
    budget.check(resultQuotient);
    return {
        std::move(resultQuotient), std::move(remainder),
        budget.failure() == RischFailure::None};
}

[[nodiscard]] std::optional<DifferentialPolynomial> exactQuotient(
    const DifferentialPolynomial& numerator,
    const DifferentialPolynomial& denominator,
    PrimitiveBudget& budget) {
    PolynomialDivision division = dividePolynomials(
        numerator, denominator, budget);
    if (!division.exact || !division.remainder.isZero()) {
        budget.fail(RischFailure::ExactDivisionFailed);
        return std::nullopt;
    }
    return std::move(division.quotient);
}

[[nodiscard]] DifferentialPolynomial remainderPolynomial(
    const DifferentialPolynomial& value,
    const DifferentialPolynomial& modulus,
    PrimitiveBudget& budget) {
    return dividePolynomials(value, modulus, budget).remainder;
}

[[nodiscard]] DifferentialPolynomial monicPolynomial(
    const DifferentialPolynomial& value,
    PrimitiveBudget& budget) {
    if (value.isZero())
        return value;
    const auto inverse = divideRationalFunctionsExact(
        oneFunction(), value.coefficient(value.degree()));
    if (!inverse) {
        budget.fail(RischFailure::ExactDivisionFailed);
        return {};
    }
    return scalePolynomial(value, *inverse, budget);
}

[[nodiscard]] DifferentialPolynomial gcdPolynomial(
    DifferentialPolynomial lhs,
    DifferentialPolynomial rhs,
    PrimitiveBudget& budget) {
    while (!rhs.isZero() && budget.failure() == RischFailure::None) {
        PolynomialDivision division = dividePolynomials(lhs, rhs, budget);
        if (!division.exact)
            return {};
        lhs = std::move(rhs);
        rhs = std::move(division.remainder);
    }
    return monicPolynomial(lhs, budget);
}

struct ExtendedGcd final {
    DifferentialPolynomial gcd;
    DifferentialPolynomial lhsCoefficient;
};

[[nodiscard]] ExtendedGcd extendedGcd(
    DifferentialPolynomial lhs,
    DifferentialPolynomial rhs,
    PrimitiveBudget& budget) {
    DifferentialPolynomial oldS = onePolynomial();
    DifferentialPolynomial s;
    while (!rhs.isZero() && budget.failure() == RischFailure::None) {
        PolynomialDivision division = dividePolynomials(lhs, rhs, budget);
        if (!division.exact)
            return {};
        DifferentialPolynomial nextS = subtractPolynomials(
            oldS,
            multiplyPolynomials(division.quotient, s, budget),
            budget);
        lhs = std::move(rhs);
        rhs = std::move(division.remainder);
        oldS = std::move(s);
        s = std::move(nextS);
    }
    if (lhs.isZero() || budget.failure() != RischFailure::None)
        return {};
    const auto inverse = divideRationalFunctionsExact(
        oneFunction(), lhs.coefficient(lhs.degree()));
    if (!inverse) {
        budget.fail(RischFailure::ExactDivisionFailed);
        return {};
    }
    return {
        scalePolynomial(lhs, *inverse, budget),
        scalePolynomial(oldS, *inverse, budget)};
}

[[nodiscard]] DifferentialPolynomial powerPolynomial(
    DifferentialPolynomial base,
    std::size_t exponent,
    PrimitiveBudget& budget) {
    DifferentialPolynomial result = onePolynomial();
    while (exponent != 0 && budget.failure() == RischFailure::None) {
        if ((exponent & 1U) != 0)
            result = multiplyPolynomials(result, base, budget);
        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyPolynomials(base, base, budget);
    }
    return result;
}

[[nodiscard]] DifferentialPolynomial formalDerivative(
    const DifferentialPolynomial& value,
    PrimitiveBudget& budget) {
    if (value.degree() == 0)
        return {};
    std::vector<RationalFunction> coefficients(value.degree());
    for (std::size_t i = 1; i <= value.degree(); ++i)
        coefficients[i - 1] = multiplyRationalFunctionsExact(
            value.coefficient(i),
            constantFunction(Rational{BigInt::fromUnsigned(i)}));
    DifferentialPolynomial result{std::move(coefficients)};
    budget.check(result, value.degree());
    return result;
}

[[nodiscard]] DifferentialPolynomial totalDerivative(
    const DifferentialPolynomial& value,
    const RationalFunction& differentialCoefficient,
    PrimitiveBudget& budget) {
    std::vector<RationalFunction> coefficients(value.coefficients().size());
    for (std::size_t i = 0; i <= value.degree(); ++i) {
        coefficients[i] = addRationalFunctionsExact(
            coefficients[i],
            differentiateRationalFunctionExact(value.coefficient(i)));
        if (i != 0)
            coefficients[i - 1] = addRationalFunctionsExact(
                coefficients[i - 1],
                multiplyRationalFunctionsExact(
                    constantFunction(Rational{BigInt::fromUnsigned(i)}),
                    multiplyRationalFunctionsExact(
                        value.coefficient(i), differentialCoefficient)));
    }
    DifferentialPolynomial result{std::move(coefficients)};
    budget.check(result, value.coefficients().size() * 2);
    return result;
}

[[nodiscard]] bool samePolynomial(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs) {
    if (lhs.degree() != rhs.degree())
        return false;
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        if (!equivalentRationalFunctions(
                lhs.coefficient(i), rhs.coefficient(i)))
            return false;
    return true;
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

[[nodiscard]] std::optional<Rational> constantMultiple(
    const DifferentialPolynomial& numerator,
    const DifferentialPolynomial& denominator) {
    std::optional<Rational> multiplier;
    const std::size_t degree = std::max(
        numerator.degree(), denominator.degree());
    for (std::size_t i = 0; i <= degree; ++i) {
        const RationalFunction& lhs = numerator.coefficient(i);
        const RationalFunction& rhs = denominator.coefficient(i);
        if (rhs.numerator.isZero()) {
            if (!lhs.numerator.isZero())
                return std::nullopt;
            continue;
        }
        const auto quotient = divideRationalFunctionsExact(lhs, rhs);
        if (!quotient)
            return std::nullopt;
        const auto value = constantValue(*quotient);
        if (!value || (multiplier && *multiplier != *value))
            return std::nullopt;
        multiplier = *value;
    }
    return multiplier;
}

struct SquareFreeFactor final {
    DifferentialPolynomial factor;
    std::size_t multiplicity = 0;
};

[[nodiscard]] std::optional<RationalPolynomial> constantCoefficientPolynomial(
    const DifferentialPolynomial& value) {
    std::vector<Rational> coefficients;
    coefficients.reserve(value.coefficients().size());
    for (const RationalFunction& coefficient : value.coefficients()) {
        const auto constant = constantValue(coefficient);
        if (!constant)
            return std::nullopt;
        coefficients.push_back(*constant);
    }
    return RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] DifferentialPolynomial differentialPolynomial(
    const RationalPolynomial& value) {
    std::vector<RationalFunction> coefficients;
    coefficients.reserve(value.coefficients().size());
    for (const Rational& coefficient : value.coefficients())
        coefficients.push_back(constantFunction(coefficient));
    return DifferentialPolynomial{std::move(coefficients)};
}

void appendConstantCoefficientFactors(
    const DifferentialPolynomial& factor,
    std::size_t multiplicity,
    std::vector<SquareFreeFactor>& result,
    PrimitiveBudget& budget) {
    const auto polynomial = constantCoefficientPolynomial(factor);
    if (!polynomial || polynomial->degree() > 64) {
        result.push_back({factor, multiplicity});
        return;
    }
    RationalPolynomialFactorOptions options;
    options.maximumFiniteFieldDegree = 64;
    options.maximumKroneckerDegree = 32;
    options.maximumCombinations = 32768;
    const RationalPolynomialFactorization factorization =
        factorRationalPolynomialOverQ(*polynomial, options);
    if (factorization.factors.size() <= 1
        || factorization.scalar != one()) {
        result.push_back({factor, multiplicity});
        return;
    }

    DifferentialPolynomial reconstructed = onePolynomial();
    std::vector<DifferentialPolynomial> exactFactors;
    exactFactors.reserve(factorization.factors.size());
    for (const RationalPolynomial& exactFactor : factorization.factors) {
        exactFactors.push_back(differentialPolynomial(exactFactor));
        reconstructed = multiplyPolynomials(
            reconstructed, exactFactors.back(), budget);
    }
    if (budget.failure() != RischFailure::None
        || !samePolynomial(reconstructed, factor)) {
        if (budget.failure() == RischFailure::None)
            budget.fail(RischFailure::CertificateFailed);
        return;
    }
    for (DifferentialPolynomial& exactFactor : exactFactors)
        result.push_back({std::move(exactFactor), multiplicity});
}

[[nodiscard]] std::optional<std::vector<SquareFreeFactor>>
squareFreeDecomposition(
    const DifferentialPolynomial& input,
    PrimitiveBudget& budget) {
    DifferentialPolynomial polynomial = monicPolynomial(input, budget);
    if (polynomial.isZero() || budget.failure() != RischFailure::None)
        return std::nullopt;
    if (polynomial.degree() == 0)
        return std::vector<SquareFreeFactor>{};
    DifferentialPolynomial repeated = gcdPolynomial(
        polynomial, formalDerivative(polynomial, budget), budget);
    auto squareFree = exactQuotient(polynomial, repeated, budget);
    if (!squareFree)
        return std::nullopt;

    std::vector<SquareFreeFactor> result;
    for (std::size_t multiplicity = 1;
         !isOnePolynomial(*squareFree)
            && budget.failure() == RischFailure::None;
         ++multiplicity) {
        DifferentialPolynomial shared = gcdPolynomial(
            *squareFree, repeated, budget);
        auto factor = exactQuotient(*squareFree, shared, budget);
        if (!factor)
            return std::nullopt;
        if (!isOnePolynomial(*factor))
            appendConstantCoefficientFactors(
                *factor, multiplicity, result, budget);
        auto nextRepeated = exactQuotient(repeated, shared, budget);
        if (!nextRepeated)
            return std::nullopt;
        *squareFree = std::move(shared);
        repeated = std::move(*nextRepeated);
        if (multiplicity >= input.degree() + 1) {
            budget.fail(RischFailure::CertificateFailed);
            return std::nullopt;
        }
    }
    if (!isOnePolynomial(repeated)) {
        budget.fail(RischFailure::CertificateFailed);
        return std::nullopt;
    }
    return result;
}

[[nodiscard]] DifferentialRationalFunction normalizedFunction(
    const DifferentialRationalFunction& input,
    PrimitiveBudget& budget) {
    if (input.denominator.isZero()) {
        budget.fail(RischFailure::ZeroDenominator);
        return {};
    }
    DifferentialPolynomial common = gcdPolynomial(
        input.numerator, input.denominator, budget);
    auto numerator = exactQuotient(input.numerator, common, budget);
    auto denominator = exactQuotient(input.denominator, common, budget);
    if (!numerator || !denominator)
        return {};
    const RationalFunction leading = denominator->coefficient(
        denominator->degree());
    const auto inverse = divideRationalFunctionsExact(oneFunction(), leading);
    if (!inverse) {
        budget.fail(RischFailure::ExactDivisionFailed);
        return {};
    }
    return {
        scalePolynomial(*numerator, *inverse, budget),
        scalePolynomial(*denominator, *inverse, budget)};
}

[[nodiscard]] DifferentialRationalFunction addFunctions(
    const DifferentialRationalFunction& lhs,
    const DifferentialRationalFunction& rhs,
    PrimitiveBudget& budget) {
    return {
        addPolynomials(
            multiplyPolynomials(lhs.numerator, rhs.denominator, budget),
            multiplyPolynomials(rhs.numerator, lhs.denominator, budget),
            budget),
        multiplyPolynomials(lhs.denominator, rhs.denominator, budget)};
}

[[nodiscard]] DifferentialRationalFunction derivativeFunction(
    const DifferentialRationalFunction& value,
    const RationalFunction& differentialCoefficient,
    PrimitiveBudget& budget) {
    const DifferentialPolynomial numeratorDerivative = totalDerivative(
        value.numerator, differentialCoefficient, budget);
    const DifferentialPolynomial denominatorDerivative = totalDerivative(
        value.denominator, differentialCoefficient, budget);
    return {
        subtractPolynomials(
            multiplyPolynomials(
                numeratorDerivative, value.denominator, budget),
            multiplyPolynomials(
                value.numerator, denominatorDerivative, budget),
            budget),
        multiplyPolynomials(value.denominator, value.denominator, budget)};
}

[[nodiscard]] bool sameFunction(
    const DifferentialRationalFunction& lhs,
    const DifferentialRationalFunction& rhs,
    PrimitiveBudget& budget) {
    if (lhs.denominator.isZero() || rhs.denominator.isZero())
        return false;
    return samePolynomial(
        multiplyPolynomials(lhs.numerator, rhs.denominator, budget),
        multiplyPolynomials(rhs.numerator, lhs.denominator, budget));
}

[[nodiscard]] bool verifyReduction(
    const DifferentialRationalFunction& input,
    const RationalFunction& differentialCoefficient,
    const PrimitiveRationalReduction& reduction,
    const RischOptions& options) {
    PrimitiveBudget budget{options};
    DifferentialRationalFunction reconstructed{
        totalDerivative(
            reduction.polynomialPart, differentialCoefficient, budget),
        onePolynomial()};
    reconstructed = addFunctions(
        reconstructed,
        DifferentialRationalFunction{
            constantPolynomial(reduction.lowerFieldRemainder),
            onePolynomial()},
        budget);
    for (const DifferentialRationalFunction& term : reduction.rationalPart)
        reconstructed = addFunctions(
            reconstructed,
            derivativeFunction(term, differentialCoefficient, budget),
            budget);
    for (const PrimitiveLogarithmicTerm& term : reduction.logarithmicPart) {
        DifferentialPolynomial numerator = totalDerivative(
            term.argument, differentialCoefficient, budget);
        numerator = scalePolynomial(
            numerator, constantFunction(term.coefficient), budget);
        reconstructed = addFunctions(
            reconstructed,
            DifferentialRationalFunction{
                std::move(numerator), term.argument},
            budget);
    }
    for (const DifferentialRationalFunction& term : reduction.residualPart)
        reconstructed = addFunctions(reconstructed, term, budget);
    return budget.failure() == RischFailure::None
        && sameFunction(input, reconstructed, budget)
        && budget.failure() == RischFailure::None;
}

} // namespace

RischStageResult<PrimitiveRationalReduction>
reducePrimitiveRationalFunction(
    const DifferentialRationalFunction& input,
    const RationalFunction& differentialCoefficient,
    const RischOptions& options) {
    if (input.denominator.isZero()
        || differentialCoefficient.denominator.isZero())
        return {std::nullopt, RischFailure::ZeroDenominator};
    PrimitiveBudget budget{options};
    DifferentialRationalFunction normalized = normalizedFunction(input, budget);
    if (budget.failure() != RischFailure::None)
        return {std::nullopt, budget.failure()};

    PolynomialDivision polynomialDivision = dividePolynomials(
        normalized.numerator, normalized.denominator, budget);
    if (!polynomialDivision.exact)
        return {std::nullopt, budget.failure()};
    DifferentialPolynomial polynomialIntegrand =
        std::move(polynomialDivision.quotient);
    const DifferentialPolynomial properNumerator =
        std::move(polynomialDivision.remainder);

    PrimitiveRationalReduction result;
    if (!properNumerator.isZero()) {
        auto factors = squareFreeDecomposition(
            normalized.denominator, budget);
        if (!factors)
            return {std::nullopt, budget.failure()};
        for (const SquareFreeFactor& entry : *factors) {
            DifferentialPolynomial blockDenominator = powerPolynomial(
                entry.factor, entry.multiplicity, budget);
            auto cofactor = exactQuotient(
                normalized.denominator, blockDenominator, budget);
            if (!cofactor)
                return {std::nullopt, budget.failure()};
            const ExtendedGcd inverse = extendedGcd(
                *cofactor, blockDenominator, budget);
            if (budget.failure() != RischFailure::None
                || inverse.gcd.degree() != 0)
                return {std::nullopt,
                    budget.failure() == RischFailure::None
                        ? RischFailure::CertificateFailed
                        : budget.failure()};
            DifferentialPolynomial blockNumerator = remainderPolynomial(
                multiplyPolynomials(
                    properNumerator, inverse.lhsCoefficient, budget),
                blockDenominator, budget);

            const DifferentialDerivation derivation{
                DifferentialExtensionKind::Primitive,
                differentialCoefficient};
            const auto classification = classifyDifferentialPolynomial(
                entry.factor, derivation, options);
            if (!classification)
                return {std::nullopt, classification.failure};
            if (classification.value->classification
                != DifferentialPolynomialClass::Normal) {
                result.residualPart.push_back({
                    std::move(blockNumerator),
                    std::move(blockDenominator)});
                continue;
            }

            auto reduced = hermiteReduceNormalDifferentialPower(
                std::move(blockNumerator), entry.factor,
                entry.multiplicity, derivation, options);
            if (!reduced)
                return {std::nullopt, reduced.failure};
            result.steps += reduced.value->steps;
            for (auto& [numerator, power] : reduced.value->rationalTerms)
                result.rationalPart.push_back({
                    std::move(numerator),
                    powerPolynomial(entry.factor, power, budget)});

            PolynomialDivision squareFreeDivision = dividePolynomials(
                reduced.value->squareFreeNumerator, entry.factor, budget);
            if (!squareFreeDivision.exact)
                return {std::nullopt, budget.failure()};
            polynomialIntegrand = addPolynomials(
                polynomialIntegrand, squareFreeDivision.quotient, budget);
            if (squareFreeDivision.remainder.isZero())
                continue;
            DifferentialPolynomial logarithmicNumerator = remainderPolynomial(
                totalDerivative(
                    entry.factor, differentialCoefficient, budget),
                entry.factor, budget);
            const auto multiplier = constantMultiple(
                squareFreeDivision.remainder, logarithmicNumerator);
            if (multiplier)
                result.logarithmicPart.push_back({
                    *multiplier, entry.factor});
            else
                result.residualPart.push_back({
                    std::move(squareFreeDivision.remainder),
                    entry.factor});
        }
    }

    if (budget.failure() != RischFailure::None)
        return {std::nullopt, budget.failure()};
    auto polynomial = reducePrimitivePolynomial(
        polynomialIntegrand, differentialCoefficient, options);
    if (!polynomial)
        return {std::nullopt, polynomial.failure};
    result.polynomialPart = std::move(polynomial.value->polynomialPart);
    result.lowerFieldRemainder =
        std::move(polynomial.value->lowerFieldRemainder);
    result.steps += polynomial.value->steps + budget.operations();
    result.exactVerified = verifyReduction(
        input, differentialCoefficient, result, options);
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

bool verifyPrimitiveRationalReduction(
    const DifferentialRationalFunction& input,
    const RationalFunction& differentialCoefficient,
    const PrimitiveRationalReduction& reduction,
    const RischOptions& options) {
    return verifyReduction(
        input, differentialCoefficient, reduction, options);
}

} // namespace mmcal::symbolic::risch
