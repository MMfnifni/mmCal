#include "risch_differential_reduction.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::symbolic::risch {
namespace {

using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] RationalPolynomial constantPolynomial(const Rational& value) {
    return RationalPolynomial{{value}};
}

[[nodiscard]] RationalFunction oneRationalFunction() {
    return RationalFunction{
        constantPolynomial(Rational{BigInt{1}}),
        constantPolynomial(Rational{BigInt{1}})};
}

[[nodiscard]] RationalFunction integerRationalFunction(std::size_t value) {
    return RationalFunction{
        constantPolynomial(Rational{BigInt::fromUnsigned(value)}),
        constantPolynomial(Rational{BigInt{1}})};
}

[[nodiscard]] bool rationalFunctionWithinBitBudget(
    const RationalFunction& value,
    const RischOptions& options) {
    if (value.denominator.isZero())
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

class DifferentialArithmetic final {
public:
    explicit DifferentialArithmetic(const RischOptions& options)
        : options_(options) {}

    [[nodiscard]] RationalFunction add(
        const RationalFunction& lhs,
        const RationalFunction& rhs) {
        if (!charge())
            return {};
        return checked(addRationalFunctionsExact(lhs, rhs));
    }

    [[nodiscard]] RationalFunction subtract(
        const RationalFunction& lhs,
        const RationalFunction& rhs) {
        if (!charge())
            return {};
        return checked(subtractRationalFunctionsExact(lhs, rhs));
    }

    [[nodiscard]] RationalFunction multiply(
        const RationalFunction& lhs,
        const RationalFunction& rhs) {
        if (!charge())
            return {};
        return checked(multiplyRationalFunctionsExact(lhs, rhs));
    }

    [[nodiscard]] std::optional<RationalFunction> divide(
        const RationalFunction& lhs,
        const RationalFunction& rhs) {
        if (!charge())
            return std::nullopt;
        auto result = divideRationalFunctionsExact(lhs, rhs);
        if (!result) {
            if (failure_ == RischFailure::None)
                failure_ = RischFailure::ExactDivisionFailed;
            return std::nullopt;
        }
        *result = checked(std::move(*result));
        if (failure_ != RischFailure::None)
            return std::nullopt;
        return result;
    }

    [[nodiscard]] RationalFunction derivative(const RationalFunction& value) {
        if (!charge())
            return {};
        return checked(differentiateRationalFunctionExact(value));
    }

    [[nodiscard]] RischFailure failure() const noexcept {
        return failure_;
    }

private:
    [[nodiscard]] bool charge() {
        if (failure_ != RischFailure::None)
            return false;
        if (operations_ >= options_.maximumDifferentialOperations) {
            failure_ = RischFailure::DifferentialOperationLimit;
            return false;
        }
        ++operations_;
        return true;
    }

    [[nodiscard]] RationalFunction checked(RationalFunction value) {
        if (value.denominator.isZero())
            failure_ = RischFailure::ZeroDenominator;
        else if (!rationalFunctionWithinBitBudget(value, options_))
            failure_ = RischFailure::IntermediateBitLimit;
        return value;
    }

    const RischOptions& options_;
    std::size_t operations_ = 0;
    RischFailure failure_ = RischFailure::None;
};

[[nodiscard]] DifferentialPolynomial addPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs,
    DifferentialArithmetic& arithmetic) {
    const std::size_t count = std::max(
        lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<RationalFunction> result(count);
    for (std::size_t i = 0; i < count; ++i)
        result[i] = arithmetic.add(lhs.coefficient(i), rhs.coefficient(i));
    return DifferentialPolynomial{std::move(result)};
}

[[nodiscard]] DifferentialPolynomial subtractPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs,
    DifferentialArithmetic& arithmetic) {
    const std::size_t count = std::max(
        lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<RationalFunction> result(count);
    for (std::size_t i = 0; i < count; ++i)
        result[i] = arithmetic.subtract(lhs.coefficient(i), rhs.coefficient(i));
    return DifferentialPolynomial{std::move(result)};
}

[[nodiscard]] DifferentialPolynomial multiplyPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs,
    DifferentialArithmetic& arithmetic) {
    if (lhs.isZero() || rhs.isZero())
        return {};
    std::vector<RationalFunction> result(
        lhs.degree() + rhs.degree() + 1);
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        for (std::size_t j = 0; j <= rhs.degree(); ++j)
            result[i + j] = arithmetic.add(
                result[i + j], arithmetic.multiply(
                    lhs.coefficient(i), rhs.coefficient(j)));
    return DifferentialPolynomial{std::move(result)};
}

[[nodiscard]] DifferentialPolynomial scalePolynomial(
    const DifferentialPolynomial& value,
    const RationalFunction& scale,
    DifferentialArithmetic& arithmetic) {
    std::vector<RationalFunction> result;
    result.reserve(value.coefficients().size());
    for (const RationalFunction& coefficient : value.coefficients())
        result.push_back(arithmetic.multiply(coefficient, scale));
    return DifferentialPolynomial{std::move(result)};
}

[[nodiscard]] DifferentialPolynomial multiplyByMonomial(
    const DifferentialPolynomial& value,
    const RationalFunction& coefficient,
    std::size_t exponent,
    DifferentialArithmetic& arithmetic) {
    std::vector<RationalFunction> result(exponent);
    result.reserve(exponent + value.coefficients().size());
    for (const RationalFunction& current : value.coefficients())
        result.push_back(arithmetic.multiply(current, coefficient));
    return DifferentialPolynomial{std::move(result)};
}

struct PolynomialDivision final {
    DifferentialPolynomial quotient;
    DifferentialPolynomial remainder;
    bool exact = false;
};

[[nodiscard]] PolynomialDivision dividePolynomials(
    const DifferentialPolynomial& numerator,
    const DifferentialPolynomial& denominator,
    DifferentialArithmetic& arithmetic) {
    if (denominator.isZero())
        return {{}, numerator, false};
    DifferentialPolynomial remainder = numerator;
    std::vector<RationalFunction> quotient(
        numerator.degree() >= denominator.degree()
            ? numerator.degree() - denominator.degree() + 1
            : 1);
    while (!remainder.isZero()
        && remainder.degree() >= denominator.degree()
        && arithmetic.failure() == RischFailure::None) {
        const std::size_t shift = remainder.degree() - denominator.degree();
        auto scale = arithmetic.divide(
            remainder.coefficient(remainder.degree()),
            denominator.coefficient(denominator.degree()));
        if (!scale)
            return {DifferentialPolynomial{std::move(quotient)}, remainder, false};
        quotient[shift] = arithmetic.add(quotient[shift], *scale);
        remainder = subtractPolynomials(
            remainder,
            multiplyByMonomial(denominator, *scale, shift, arithmetic),
            arithmetic);
    }
    return {
        DifferentialPolynomial{std::move(quotient)},
        std::move(remainder),
        arithmetic.failure() == RischFailure::None};
}

[[nodiscard]] DifferentialPolynomial polynomialRemainder(
    const DifferentialPolynomial& value,
    const DifferentialPolynomial& modulus,
    DifferentialArithmetic& arithmetic) {
    return dividePolynomials(value, modulus, arithmetic).remainder;
}

[[nodiscard]] DifferentialPolynomial powerPolynomial(
    DifferentialPolynomial base,
    std::size_t exponent,
    DifferentialArithmetic& arithmetic) {
    DifferentialPolynomial result{{oneRationalFunction()}};
    while (exponent != 0 && arithmetic.failure() == RischFailure::None) {
        if ((exponent & 1U) != 0)
            result = multiplyPolynomials(result, base, arithmetic);
        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyPolynomials(base, base, arithmetic);
    }
    return result;
}

struct ExtendedGcd final {
    DifferentialPolynomial gcd;
    DifferentialPolynomial lhsCoefficient;
};

[[nodiscard]] ExtendedGcd extendedGcd(
    DifferentialPolynomial lhs,
    DifferentialPolynomial rhs,
    DifferentialArithmetic& arithmetic) {
    DifferentialPolynomial oldS{{oneRationalFunction()}};
    DifferentialPolynomial s;
    while (!rhs.isZero() && arithmetic.failure() == RischFailure::None) {
        PolynomialDivision division = dividePolynomials(lhs, rhs, arithmetic);
        if (!division.exact)
            break;
        DifferentialPolynomial nextS = subtractPolynomials(
            oldS, multiplyPolynomials(division.quotient, s, arithmetic), arithmetic);
        lhs = std::move(rhs);
        rhs = std::move(division.remainder);
        oldS = std::move(s);
        s = std::move(nextS);
    }
    if (lhs.isZero() || arithmetic.failure() != RischFailure::None)
        return {};
    auto inverseLeading = arithmetic.divide(
        oneRationalFunction(), lhs.coefficient(lhs.degree()));
    if (!inverseLeading)
        return {};
    DifferentialPolynomial normalizedGcd = scalePolynomial(
        lhs, *inverseLeading, arithmetic);
    DifferentialPolynomial normalizedLhsCoefficient = scalePolynomial(
        oldS, *inverseLeading, arithmetic);
    return {std::move(normalizedGcd), std::move(normalizedLhsCoefficient)};
}

[[nodiscard]] DifferentialPolynomial differentiatePolynomial(
    const DifferentialPolynomial& polynomial,
    const DifferentialDerivation& derivation,
    DifferentialArithmetic& arithmetic) {
    std::vector<RationalFunction> result(polynomial.coefficients().size());
    for (std::size_t i = 0; i <= polynomial.degree(); ++i) {
        result[i] = arithmetic.add(
            result[i], arithmetic.derivative(polynomial.coefficient(i)));
        if (i == 0)
            continue;
        RationalFunction generatorTerm = arithmetic.multiply(
            polynomial.coefficient(i), derivation.differentialCoefficient);
        generatorTerm = arithmetic.multiply(
            generatorTerm, integerRationalFunction(i));
        const std::size_t target =
            derivation.kind == DifferentialExtensionKind::Primitive ? i - 1 : i;
        result[target] = arithmetic.add(result[target], generatorTerm);
    }
    return DifferentialPolynomial{std::move(result)};
}

[[nodiscard]] bool samePolynomial(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs) {
    if (lhs.degree() != rhs.degree())
        return false;
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        if (!equivalentRationalFunctions(lhs.coefficient(i), rhs.coefficient(i)))
            return false;
    return true;
}

[[nodiscard]] bool validPolynomial(const DifferentialPolynomial& polynomial) {
    for (const RationalFunction& coefficient : polynomial.coefficients())
        if (coefficient.denominator.isZero())
            return false;
    return true;
}

[[nodiscard]] bool validDerivation(const DifferentialDerivation& derivation) {
    return !derivation.differentialCoefficient.denominator.isZero();
}

[[nodiscard]] RischFailure arithmeticFailureOr(
    const DifferentialArithmetic& arithmetic,
    RischFailure fallback) {
    return arithmetic.failure() == RischFailure::None
        ? fallback
        : arithmetic.failure();
}

[[nodiscard]] DifferentialPolynomialClassification classify(
    const DifferentialPolynomial& polynomial,
    const DifferentialDerivation& derivation,
    DifferentialArithmetic& arithmetic) {
    if (polynomial.degree() == 0)
        return {
            DifferentialPolynomialClass::Constant,
            DifferentialPolynomial{{oneRationalFunction()}},
            true};
    const DifferentialPolynomial derivative = differentiatePolynomial(
        polynomial, derivation, arithmetic);
    const ExtendedGcd gcd = extendedGcd(polynomial, derivative, arithmetic);
    if (arithmetic.failure() != RischFailure::None || gcd.gcd.isZero())
        return {};
    DifferentialPolynomialClass kind = DifferentialPolynomialClass::Mixed;
    if (gcd.gcd.degree() == 0)
        kind = DifferentialPolynomialClass::Normal;
    else if (gcd.gcd.degree() == polynomial.degree())
        kind = DifferentialPolynomialClass::Special;
    const bool verified = dividePolynomials(
        polynomial, gcd.gcd, arithmetic).remainder.isZero()
        && dividePolynomials(derivative, gcd.gcd, arithmetic).remainder.isZero();
    return {kind, gcd.gcd, verified};
}

struct ReductionVerification final {
    bool verified = false;
    RischFailure failure = RischFailure::None;
};

[[nodiscard]] ReductionVerification verifyReduction(
    const DifferentialPolynomial& originalNumerator,
    const DifferentialPolynomial& factor,
    std::size_t originalPower,
    const DifferentialDerivation& derivation,
    const DifferentialNormalReduction& reduction,
    const RischOptions& options) {
    if (factor.isZero() || factor.degree() == 0 || originalPower == 0
        || !validPolynomial(originalNumerator) || !validPolynomial(factor)
        || !validDerivation(derivation))
        return {};
    DifferentialArithmetic arithmetic{options};
    const DifferentialPolynomial factorDerivative = differentiatePolynomial(
        factor, derivation, arithmetic);
    DifferentialPolynomial reconstructed = multiplyPolynomials(
        reduction.squareFreeNumerator,
        powerPolynomial(factor, originalPower - 1, arithmetic),
        arithmetic);
    for (const auto& [numerator, power] : reduction.rationalTerms) {
        if (power == 0 || power >= originalPower)
            return {};
        const DifferentialPolynomial derivative = differentiatePolynomial(
            numerator, derivation, arithmetic);
        const DifferentialPolynomial derivativeContribution =
            multiplyPolynomials(derivative, factor, arithmetic);
        const DifferentialPolynomial factorDerivativeContribution =
            multiplyPolynomials(numerator, factorDerivative, arithmetic);
        const DifferentialPolynomial scaledFactorDerivativeContribution =
            scalePolynomial(
                factorDerivativeContribution,
                integerRationalFunction(power), arithmetic);
        DifferentialPolynomial differentiatedNumerator = subtractPolynomials(
            derivativeContribution,
            scaledFactorDerivativeContribution, arithmetic);
        differentiatedNumerator = multiplyPolynomials(
            differentiatedNumerator,
            powerPolynomial(factor, originalPower - power - 1, arithmetic),
            arithmetic);
        reconstructed = addPolynomials(
            reconstructed, differentiatedNumerator, arithmetic);
    }
    if (arithmetic.failure() != RischFailure::None)
        return {false, arithmetic.failure()};
    return {
        samePolynomial(originalNumerator, reconstructed),
        RischFailure::None};
}

} // namespace

DifferentialPolynomial::DifferentialPolynomial()
    : coefficients_(1) {}

DifferentialPolynomial::DifferentialPolynomial(
    std::vector<RationalFunction> coefficients)
    : coefficients_(std::move(coefficients)) {
    normalize();
}

bool DifferentialPolynomial::isZero() const noexcept {
    return coefficients_.size() == 1 && coefficients_.front().numerator.isZero();
}

std::size_t DifferentialPolynomial::degree() const noexcept {
    return coefficients_.size() - 1;
}

const RationalFunction& DifferentialPolynomial::coefficient(
    std::size_t exponent) const noexcept {
    static const RationalFunction zero;
    return exponent < coefficients_.size() ? coefficients_[exponent] : zero;
}

const std::vector<RationalFunction>&
DifferentialPolynomial::coefficients() const noexcept {
    return coefficients_;
}

void DifferentialPolynomial::normalize() {
    for (RationalFunction& coefficient : coefficients_)
        coefficient = canonicalizeRationalFunction(std::move(coefficient));
    while (coefficients_.size() > 1 && coefficients_.back().numerator.isZero())
        coefficients_.pop_back();
    if (coefficients_.empty())
        coefficients_.emplace_back();
}

RischStageResult<DifferentialPolynomial> differentiateDifferentialPolynomial(
    const DifferentialPolynomial& polynomial,
    const DifferentialDerivation& derivation,
    const RischOptions& options) {
    if (!validPolynomial(polynomial) || !validDerivation(derivation))
        return {std::nullopt, RischFailure::ZeroDenominator};
    if (polynomial.degree() > options.maximumHermiteDegree)
        return {std::nullopt, RischFailure::DegreeLimit};
    DifferentialArithmetic arithmetic{options};
    DifferentialPolynomial result = differentiatePolynomial(
        polynomial, derivation, arithmetic);
    if (arithmetic.failure() != RischFailure::None)
        return {std::nullopt, arithmetic.failure()};
    return {std::move(result), RischFailure::None};
}

RischStageResult<DifferentialPolynomialClassification>
classifyDifferentialPolynomial(
    const DifferentialPolynomial& polynomial,
    const DifferentialDerivation& derivation,
    const RischOptions& options) {
    if (polynomial.isZero())
        return {std::nullopt, RischFailure::InvalidRationalFunction};
    if (!validPolynomial(polynomial) || !validDerivation(derivation))
        return {std::nullopt, RischFailure::ZeroDenominator};
    if (polynomial.degree() > options.maximumHermiteDegree)
        return {std::nullopt, RischFailure::DegreeLimit};
    DifferentialArithmetic arithmetic{options};
    DifferentialPolynomialClassification result = classify(
        polynomial, derivation, arithmetic);
    if (arithmetic.failure() != RischFailure::None)
        return {std::nullopt, arithmetic.failure()};
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

RischStageResult<DifferentialNormalReduction>
hermiteReduceNormalDifferentialPower(
    DifferentialPolynomial numerator,
    const DifferentialPolynomial& factor,
    std::size_t denominatorPower,
    const DifferentialDerivation& derivation,
    const RischOptions& options) {
    if (denominatorPower == 0 || factor.isZero() || factor.degree() == 0)
        return {std::nullopt, RischFailure::InvalidRationalFunction};
    if (!validPolynomial(numerator) || !validPolynomial(factor)
        || !validDerivation(derivation))
        return {std::nullopt, RischFailure::ZeroDenominator};
    if (factor.degree() > options.maximumHermiteDegree)
        return {std::nullopt, RischFailure::DegreeLimit};

    const DifferentialPolynomial originalNumerator = numerator;
    const std::size_t originalPower = denominatorPower;
    DifferentialArithmetic arithmetic{options};
    const DifferentialPolynomialClassification classification = classify(
        factor, derivation, arithmetic);
    if (arithmetic.failure() != RischFailure::None)
        return {std::nullopt, arithmetic.failure()};
    if (!classification.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    if (classification.classification != DifferentialPolynomialClass::Normal)
        return {std::nullopt, RischFailure::NonNormalDifferentialFactor};

    const DifferentialPolynomial factorDerivative = differentiatePolynomial(
        factor, derivation, arithmetic);
    const ExtendedGcd inverseData = extendedGcd(
        factorDerivative, factor, arithmetic);
    if (arithmetic.failure() != RischFailure::None)
        return {std::nullopt, arithmetic.failure()};
    if (inverseData.gcd.degree() != 0)
        return {std::nullopt, RischFailure::NonNormalDifferentialFactor};
    auto inverseConstant = arithmetic.divide(
        oneRationalFunction(), inverseData.gcd.coefficient(0));
    if (!inverseConstant)
        return {std::nullopt, arithmeticFailureOr(
            arithmetic, RischFailure::ExactDivisionFailed)};
    const DifferentialPolynomial inverseDerivative = polynomialRemainder(
        scalePolynomial(inverseData.lhsCoefficient, *inverseConstant, arithmetic),
        factor, arithmetic);

    DifferentialNormalReduction result;
    while (denominatorPower > 1 && !numerator.isZero()
        && arithmetic.failure() == RischFailure::None) {
        DifferentialPolynomial correction = polynomialRemainder(
            multiplyPolynomials(numerator, inverseDerivative, arithmetic),
            factor, arithmetic);
        auto negativeScale = arithmetic.divide(
            RationalFunction{
                constantPolynomial(Rational{BigInt{-1}}),
                constantPolynomial(Rational{BigInt{1}})},
            integerRationalFunction(denominatorPower - 1));
        if (!negativeScale)
            return {std::nullopt, arithmeticFailureOr(
                arithmetic, RischFailure::ExactDivisionFailed)};
        correction = scalePolynomial(correction, *negativeScale, arithmetic);

        const DifferentialPolynomial correctionDerivative =
            differentiatePolynomial(correction, derivation, arithmetic);
        const DifferentialPolynomial derivativeContribution =
            multiplyPolynomials(correctionDerivative, factor, arithmetic);
        const DifferentialPolynomial factorDerivativeContribution =
            multiplyPolynomials(correction, factorDerivative, arithmetic);
        const DifferentialPolynomial scaledFactorDerivativeContribution =
            scalePolynomial(
                factorDerivativeContribution,
                integerRationalFunction(denominatorPower - 1), arithmetic);
        const DifferentialPolynomial residual = addPolynomials(
            subtractPolynomials(
                numerator, derivativeContribution, arithmetic),
            scaledFactorDerivativeContribution, arithmetic);
        PolynomialDivision lowered = dividePolynomials(
            residual, factor, arithmetic);
        if (arithmetic.failure() != RischFailure::None)
            return {std::nullopt, arithmetic.failure()};
        if (!lowered.exact || !lowered.remainder.isZero())
            return {std::nullopt, RischFailure::ExactDivisionFailed};
        result.rationalTerms.emplace_back(
            std::move(correction), denominatorPower - 1);
        numerator = std::move(lowered.quotient);
        --denominatorPower;
        ++result.steps;
    }
    result.squareFreeNumerator = std::move(numerator);
    if (arithmetic.failure() != RischFailure::None)
        return {std::nullopt, arithmetic.failure()};
    const ReductionVerification verification = verifyReduction(
        originalNumerator, factor, originalPower,
        derivation, result, options);
    if (verification.failure != RischFailure::None)
        return {std::nullopt, verification.failure};
    result.exactVerified = verification.verified;
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

bool verifyNormalDifferentialReduction(
    const DifferentialPolynomial& originalNumerator,
    const DifferentialPolynomial& factor,
    std::size_t originalPower,
    const DifferentialDerivation& derivation,
    const DifferentialNormalReduction& reduction,
    const RischOptions& options) {
    return verifyReduction(
        originalNumerator, factor, originalPower,
        derivation, reduction, options).verified;
}

bool equivalentDifferentialPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs) {
    return samePolynomial(lhs, rhs);
}

} // namespace mmcal::symbolic::risch
