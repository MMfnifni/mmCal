#include "risch_core.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <limits>
#include <utility>

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

[[nodiscard]] bool samePolynomial(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    return lhs.coefficients() == rhs.coefficients();
}

struct PolynomialDivision final {
    RationalPolynomial quotient;
    RationalPolynomial remainder;
};

[[nodiscard]] PolynomialDivision dividePolynomials(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator) {
    if (denominator.isZero())
        return {RationalPolynomial{}, numerator};

    std::vector<Rational> remainder = numerator.coefficients();
    std::vector<Rational> quotient(
        numerator.degree() >= denominator.degree()
            ? numerator.degree() - denominator.degree() + 1
            : 1,
        zero());
    const auto trim = [](std::vector<Rational>& coefficients) {
        while (coefficients.size() > 1 && coefficients.back().isZero())
            coefficients.pop_back();
        if (coefficients.empty())
            coefficients.push_back(Rational{BigInt{0}});
    };
    trim(remainder);

    const Rational leading = denominator.coefficient(denominator.degree());
    while (!(remainder.size() == 1 && remainder.front().isZero())
        && remainder.size() - 1 >= denominator.degree()) {
        const std::size_t shift = remainder.size() - 1 - denominator.degree();
        const Rational scale = remainder.back() / leading;
        quotient[shift] += scale;
        for (std::size_t i = 0; i <= denominator.degree(); ++i)
            remainder[i + shift] -= scale * denominator.coefficient(i);
        trim(remainder);
    }
    return {
        RationalPolynomial{std::move(quotient)},
        RationalPolynomial{std::move(remainder)}};
}

[[nodiscard]] RationalPolynomial addPolynomials(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    const std::size_t count = std::max(
        lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<Rational> result(count, zero());
    for (std::size_t i = 0; i < count; ++i)
        result[i] = lhs.coefficient(i) + rhs.coefficient(i);
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial negatePolynomial(
    const RationalPolynomial& value) {
    std::vector<Rational> result = value.coefficients();
    for (Rational& coefficient : result)
        coefficient = -coefficient;
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial subtractPolynomials(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    return addPolynomials(lhs, negatePolynomial(rhs));
}

[[nodiscard]] RationalPolynomial scalePolynomial(
    const RationalPolynomial& value,
    const Rational& scale) {
    std::vector<Rational> result = value.coefficients();
    for (Rational& coefficient : result)
        coefficient *= scale;
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial multiplyPolynomials(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    if (lhs.isZero() || rhs.isZero())
        return RationalPolynomial{};
    std::vector<Rational> result(
        lhs.degree() + rhs.degree() + 1, zero());
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        for (std::size_t j = 0; j <= rhs.degree(); ++j)
            result[i + j] += lhs.coefficient(i) * rhs.coefficient(j);
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial powerPolynomial(
    RationalPolynomial base,
    std::size_t exponent) {
    RationalPolynomial result = onePolynomial();
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = multiplyPolynomials(result, base);
        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyPolynomials(base, base);
    }
    return result;
}

[[nodiscard]] RationalPolynomial monicPolynomial(
    const RationalPolynomial& value) {
    if (value.isZero())
        return value;
    return scalePolynomial(
        value, one() / value.coefficient(value.degree()));
}

[[nodiscard]] RationalPolynomial differentiatePolynomial(
    const RationalPolynomial& value) {
    if (value.degree() == 0)
        return RationalPolynomial{};
    std::vector<Rational> result(value.degree(), zero());
    for (std::size_t i = 1; i <= value.degree(); ++i)
        result[i - 1] = value.coefficient(i)
            * Rational{BigInt::fromUnsigned(i)};
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial integratePolynomial(
    const RationalPolynomial& value) {
    if (value.isZero())
        return RationalPolynomial{};
    std::vector<Rational> result(value.degree() + 2, zero());
    for (std::size_t i = 0; i <= value.degree(); ++i)
        result[i + 1] = value.coefficient(i)
            / Rational{BigInt::fromUnsigned(i + 1)};
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial polynomialGcdMonic(
    RationalPolynomial lhs,
    RationalPolynomial rhs) {
    while (!rhs.isZero()) {
        PolynomialDivision division = dividePolynomials(lhs, rhs);
        lhs = std::move(rhs);
        rhs = std::move(division.remainder);
    }
    if (lhs.isZero())
        return onePolynomial();
    return monicPolynomial(lhs);
}

struct ExtendedGcd final {
    RationalPolynomial gcd;
    RationalPolynomial lhsCoefficient;
    RationalPolynomial rhsCoefficient;
};

[[nodiscard]] ExtendedGcd extendedGcd(
    RationalPolynomial lhs,
    RationalPolynomial rhs) {
    RationalPolynomial oldS = onePolynomial();
    RationalPolynomial s;
    RationalPolynomial oldT;
    RationalPolynomial t = onePolynomial();
    while (!rhs.isZero()) {
        PolynomialDivision division = dividePolynomials(lhs, rhs);
        RationalPolynomial nextS = subtractPolynomials(
            oldS, multiplyPolynomials(division.quotient, s));
        RationalPolynomial nextT = subtractPolynomials(
            oldT, multiplyPolynomials(division.quotient, t));
        lhs = std::move(rhs);
        rhs = std::move(division.remainder);
        oldS = std::move(s);
        s = std::move(nextS);
        oldT = std::move(t);
        t = std::move(nextT);
    }
    if (lhs.isZero())
        return {};
    const Rational scale = one() / lhs.coefficient(lhs.degree());
    return {
        scalePolynomial(lhs, scale),
        scalePolynomial(oldS, scale),
        scalePolynomial(oldT, scale)};
}

[[nodiscard]] std::optional<RationalPolynomial> inverseModulo(
    const RationalPolynomial& value,
    const RationalPolynomial& modulus) {
    ExtendedGcd result = extendedGcd(value, modulus);
    if (result.gcd.degree() != 0 || result.gcd.coefficient(0).isZero())
        return std::nullopt;
    return dividePolynomials(result.lhsCoefficient, modulus).remainder;
}

[[nodiscard]] bool exactPolynomialQuotient(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator,
    RationalPolynomial& quotient) {
    if (denominator.isZero())
        return false;
    PolynomialDivision division = dividePolynomials(numerator, denominator);
    if (!division.remainder.isZero())
        return false;
    quotient = std::move(division.quotient);
    return true;
}

[[nodiscard]] bool polynomialWithinBudget(
    const RationalPolynomial& polynomial,
    const RischOptions& options) {
    if (polynomial.coefficients().size() > options.maximumBivariateCoefficients)
        return false;
    for (const Rational& coefficient : polynomial.coefficients())
        if (coefficient.numerator().bitLength() > options.maximumIntermediateBits
            || coefficient.denominator().bitLength() > options.maximumIntermediateBits)
            return false;
    return true;
}

[[nodiscard]] RischFailure polynomialBudgetFailure(
    const RationalPolynomial& polynomial,
    const RischOptions& options) {
    if (polynomial.coefficients().size() > options.maximumBivariateCoefficients)
        return RischFailure::CoefficientCountLimit;
    if (!polynomialWithinBudget(polynomial, options))
        return RischFailure::IntermediateBitLimit;
    return RischFailure::None;
}

[[nodiscard]] RationalFunction normalizeRationalFunction(
    RationalFunction value) {
    if (value.denominator.isZero())
        return value;
    if (value.numerator.isZero())
        return RationalFunction{};
    RationalPolynomial gcd = polynomialGcdMonic(
        value.numerator, value.denominator);
    RationalPolynomial numerator;
    RationalPolynomial denominator;
    if (!exactPolynomialQuotient(value.numerator, gcd, numerator)
        || !exactPolynomialQuotient(value.denominator, gcd, denominator))
        return value;
    const Rational scale = one() / denominator.coefficient(denominator.degree());
    return RationalFunction{
        scalePolynomial(numerator, scale),
        scalePolynomial(denominator, scale)};
}

[[nodiscard]] RationalFunction addRationalFunctions(
    const RationalFunction& lhs,
    const RationalFunction& rhs) {
    return normalizeRationalFunction(RationalFunction{
        addPolynomials(
            multiplyPolynomials(lhs.numerator, rhs.denominator),
            multiplyPolynomials(rhs.numerator, lhs.denominator)),
        multiplyPolynomials(lhs.denominator, rhs.denominator)});
}

[[nodiscard]] RationalFunction subtractRationalFunctions(
    const RationalFunction& lhs,
    const RationalFunction& rhs) {
    return normalizeRationalFunction(RationalFunction{
        subtractPolynomials(
            multiplyPolynomials(lhs.numerator, rhs.denominator),
            multiplyPolynomials(rhs.numerator, lhs.denominator)),
        multiplyPolynomials(lhs.denominator, rhs.denominator)});
}

[[nodiscard]] RationalFunction multiplyRationalFunctions(
    const RationalFunction& lhs,
    const RationalFunction& rhs) {
    return normalizeRationalFunction(RationalFunction{
        multiplyPolynomials(lhs.numerator, rhs.numerator),
        multiplyPolynomials(lhs.denominator, rhs.denominator)});
}

[[nodiscard]] std::optional<RationalFunction> divideRationalFunctions(
    const RationalFunction& numerator,
    const RationalFunction& denominator) {
    if (numerator.denominator.isZero()
        || denominator.denominator.isZero() || denominator.numerator.isZero())
        return std::nullopt;
    return normalizeRationalFunction(RationalFunction{
        multiplyPolynomials(numerator.numerator, denominator.denominator),
        multiplyPolynomials(numerator.denominator, denominator.numerator)});
}

[[nodiscard]] RationalFunction differentiateRationalFunction(
    const RationalFunction& value) {
    return normalizeRationalFunction(RationalFunction{
        subtractPolynomials(
            multiplyPolynomials(
                differentiatePolynomial(value.numerator), value.denominator),
            multiplyPolynomials(
                value.numerator, differentiatePolynomial(value.denominator))),
        multiplyPolynomials(value.denominator, value.denominator)});
}

[[nodiscard]] bool sameRationalFunction(
    const RationalFunction& lhs,
    const RationalFunction& rhs) {
    if (lhs.denominator.isZero() || rhs.denominator.isZero())
        return false;
    return samePolynomial(
        multiplyPolynomials(lhs.numerator, rhs.denominator),
        multiplyPolynomials(rhs.numerator, lhs.denominator));
}

struct SquareFreeFactor final {
    RationalPolynomial polynomial;
    std::size_t multiplicity = 0;
};

[[nodiscard]] RischStageResult<std::vector<SquareFreeFactor>> squareFreeDecomposition(
    RationalPolynomial polynomial,
    const RischOptions& options) {
    if (polynomial.isZero() || polynomial.degree() == 0)
        return {std::nullopt, RischFailure::InvalidRationalFunction};
    polynomial = monicPolynomial(polynomial);
    RationalPolynomial repeated = polynomialGcdMonic(
        polynomial, differentiatePolynomial(polynomial));
    RationalPolynomial remaining;
    if (!exactPolynomialQuotient(polynomial, repeated, remaining))
        return {std::nullopt, RischFailure::ExactDivisionFailed};

    std::vector<SquareFreeFactor> factors;
    for (std::size_t multiplicity = 1;
         !(remaining.degree() == 0 && !remaining.coefficient(0).isZero());
         ++multiplicity) {
        if (multiplicity > polynomial.degree())
            return {std::nullopt, RischFailure::ExactDivisionFailed};
        RationalPolynomial shared = polynomialGcdMonic(remaining, repeated);
        RationalPolynomial distinct;
        if (!exactPolynomialQuotient(remaining, shared, distinct))
            return {std::nullopt, RischFailure::ExactDivisionFailed};
        if (distinct.degree() != 0)
            factors.push_back({std::move(distinct), multiplicity});
        remaining = std::move(shared);
        RationalPolynomial nextRepeated;
        if (!exactPolynomialQuotient(repeated, remaining, nextRepeated))
            return {std::nullopt, RischFailure::ExactDivisionFailed};
        repeated = std::move(nextRepeated);
        if (factors.size() > options.maximumSubresultantSteps)
            return {std::nullopt, RischFailure::SubresultantStepLimit};
    }
    return {std::move(factors), RischFailure::None};
}

[[nodiscard]] BivariateRationalPolynomial bivariateScale(
    const BivariateRationalPolynomial& value,
    const RationalPolynomial& scale) {
    std::vector<RationalPolynomial> result;
    result.reserve(value.coefficientsInX().size());
    for (const RationalPolynomial& coefficient : value.coefficientsInX())
        result.push_back(multiplyPolynomials(coefficient, scale));
    return BivariateRationalPolynomial{std::move(result)};
}

[[nodiscard]] BivariateRationalPolynomial bivariateSubtract(
    const BivariateRationalPolynomial& lhs,
    const BivariateRationalPolynomial& rhs) {
    const std::size_t count = std::max(
        lhs.coefficientsInX().size(), rhs.coefficientsInX().size());
    std::vector<RationalPolynomial> result(count);
    for (std::size_t i = 0; i < count; ++i)
        result[i] = subtractPolynomials(
            lhs.coefficientInX(i), rhs.coefficientInX(i));
    return BivariateRationalPolynomial{std::move(result)};
}

[[nodiscard]] BivariateRationalPolynomial bivariateMultiplyMonomial(
    const BivariateRationalPolynomial& value,
    const RationalPolynomial& coefficient,
    std::size_t xExponent) {
    std::vector<RationalPolynomial> result(xExponent);
    result.reserve(xExponent + value.coefficientsInX().size());
    for (const RationalPolynomial& current : value.coefficientsInX())
        result.push_back(multiplyPolynomials(current, coefficient));
    return BivariateRationalPolynomial{std::move(result)};
}

[[nodiscard]] std::size_t bivariateCoefficientCount(
    const BivariateRationalPolynomial& value) {
    std::size_t count = 0;
    for (const RationalPolynomial& coefficient : value.coefficientsInX()) {
        if (count > std::numeric_limits<std::size_t>::max()
                - coefficient.coefficients().size())
            return std::numeric_limits<std::size_t>::max();
        count += coefficient.coefficients().size();
    }
    return count;
}

[[nodiscard]] RischFailure bivariateBudgetFailure(
    const BivariateRationalPolynomial& value,
    const RischOptions& options) {
    if (bivariateCoefficientCount(value) > options.maximumBivariateCoefficients)
        return RischFailure::CoefficientCountLimit;
    for (const RationalPolynomial& coefficient : value.coefficientsInX()) {
        const RischFailure failure = polynomialBudgetFailure(coefficient, options);
        if (failure != RischFailure::None)
            return failure;
    }
    return RischFailure::None;
}

[[nodiscard]] BivariateRationalPolynomial pseudoRemainder(
    const BivariateRationalPolynomial& dividend,
    const BivariateRationalPolynomial& divisor) {
    if (divisor.isZero() || dividend.isZero()
        || dividend.degreeInX() < divisor.degreeInX())
        return dividend;
    BivariateRationalPolynomial remainder = dividend;
    const RationalPolynomial divisorLeading = divisor.coefficientInX(divisor.degreeInX());
    std::size_t remainingScales = dividend.degreeInX() - divisor.degreeInX() + 1;
    while (!remainder.isZero() && remainder.degreeInX() >= divisor.degreeInX()) {
        const std::size_t shift = remainder.degreeInX() - divisor.degreeInX();
        const RationalPolynomial leading = remainder.coefficientInX(remainder.degreeInX());
        remainder = bivariateSubtract(
            bivariateScale(remainder, divisorLeading),
            bivariateMultiplyMonomial(divisor, leading, shift));
        --remainingScales;
    }
    if (remainingScales != 0)
        remainder = bivariateScale(
            remainder, powerPolynomial(divisorLeading, remainingScales));
    return remainder;
}

[[nodiscard]] std::optional<BivariateRationalPolynomial> exactBivariateCoefficientQuotient(
    const BivariateRationalPolynomial& value,
    const RationalPolynomial& divisor) {
    std::vector<RationalPolynomial> result;
    result.reserve(value.coefficientsInX().size());
    for (const RationalPolynomial& coefficient : value.coefficientsInX()) {
        RationalPolynomial quotient;
        if (!exactPolynomialQuotient(coefficient, divisor, quotient))
            return std::nullopt;
        result.push_back(std::move(quotient));
    }
    return BivariateRationalPolynomial{std::move(result)};
}

[[nodiscard]] RischStageResult<std::vector<BivariateRationalPolynomial>>
subresultantPolynomialRemainderSequence(
    BivariateRationalPolynomial first,
    BivariateRationalPolynomial second,
    const RischOptions& options) {
    if (first.isZero() || second.isZero())
        return {std::nullopt, RischFailure::InvalidRationalFunction};
    if (first.degreeInX() < second.degreeInX())
        std::swap(first, second);

    std::vector<BivariateRationalPolynomial> sequence{first, second};
    std::size_t m = second.degreeInX();
    std::size_t delta = first.degreeInX() - m;
    RationalPolynomial beta = constantPolynomial(
        ((delta + 1) & 1U) == 0 ? one() : -one());
    BivariateRationalPolynomial remainder = bivariateScale(
        pseudoRemainder(first, second), beta);
    RationalPolynomial leading = second.coefficientInX(second.degreeInX());
    RationalPolynomial psi = negatePolynomial(powerPolynomial(leading, delta));

    while (!remainder.isZero()) {
        if (sequence.size() >= options.maximumSubresultantSteps)
            return {std::nullopt, RischFailure::SubresultantStepLimit};
        const RischFailure budgetFailure = bivariateBudgetFailure(remainder, options);
        if (budgetFailure != RischFailure::None)
            return {std::nullopt, budgetFailure};
        const std::size_t k = remainder.degreeInX();
        sequence.push_back(remainder);

        leading = second.coefficientInX(second.degreeInX());
        RationalPolynomial denominator = powerPolynomial(psi, delta - 1);
        RationalPolynomial nextPsi;
        const RationalPolynomial psiNumerator =
            powerPolynomial(negatePolynomial(leading), delta);
        if (!exactPolynomialQuotient(
                psiNumerator, denominator, nextPsi))
            return {std::nullopt, RischFailure::ExactDivisionFailed};

        first = std::move(second);
        second = std::move(remainder);
        delta = m - k;
        m = k;
        beta = negatePolynomial(multiplyPolynomials(
            leading, powerPolynomial(nextPsi, delta)));
        remainder = pseudoRemainder(first, second);
        auto divided = exactBivariateCoefficientQuotient(remainder, beta);
        if (!divided)
            return {std::nullopt, RischFailure::ExactDivisionFailed};
        remainder = std::move(*divided);
        psi = std::move(nextPsi);
    }
    return {std::move(sequence), RischFailure::None};
}

[[nodiscard]] RischStageResult<RationalPolynomial> resultantByBareiss(
    const BivariateRationalPolynomial& first,
    const BivariateRationalPolynomial& second,
    const RischOptions& options) {
    const std::size_t m = first.degreeInX();
    const std::size_t n = second.degreeInX();
    if (first.isZero() || second.isZero())
        return {RationalPolynomial{}, RischFailure::None};
    if (n == 0)
        return {
            powerPolynomial(second.coefficientInX(0), m),
            RischFailure::None};
    if (m == 0)
        return {
            powerPolynomial(first.coefficientInX(0), n),
            RischFailure::None};
    const std::size_t dimension = m + n;
    if (dimension != 0
        && dimension > options.maximumResultantMatrixEntries / dimension)
        return {std::nullopt, RischFailure::MatrixSizeLimit};

    std::vector<std::vector<RationalPolynomial>> matrix(
        dimension, std::vector<RationalPolynomial>(dimension));
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t j = 0; j <= m; ++j)
            matrix[row][row + j] = first.coefficientInX(m - j);
    for (std::size_t row = 0; row < m; ++row)
        for (std::size_t j = 0; j <= n; ++j)
            matrix[n + row][row + j] = second.coefficientInX(n - j);

    RationalPolynomial previousPivot = onePolynomial();
    bool negateResult = false;
    for (std::size_t k = 0; k + 1 < dimension; ++k) {
        std::size_t pivotRow = k;
        while (pivotRow < dimension && matrix[pivotRow][k].isZero())
            ++pivotRow;
        if (pivotRow == dimension)
            return {RationalPolynomial{}, RischFailure::None};
        if (pivotRow != k) {
            std::swap(matrix[pivotRow], matrix[k]);
            negateResult = !negateResult;
        }
        const RationalPolynomial pivot = matrix[k][k];
        for (std::size_t i = k + 1; i < dimension; ++i) {
            for (std::size_t j = k + 1; j < dimension; ++j) {
                RationalPolynomial numerator = subtractPolynomials(
                    multiplyPolynomials(matrix[i][j], pivot),
                    multiplyPolynomials(matrix[i][k], matrix[k][j]));
                RationalPolynomial quotient;
                if (!exactPolynomialQuotient(
                        numerator, previousPivot, quotient))
                    return {std::nullopt, RischFailure::ExactDivisionFailed};
                const RischFailure failure = polynomialBudgetFailure(quotient, options);
                if (failure != RischFailure::None)
                    return {std::nullopt, failure};
                matrix[i][j] = std::move(quotient);
            }
        }
        previousPivot = pivot;
    }
    RationalPolynomial determinant = matrix.back().back();
    if (negateResult)
        determinant = negatePolynomial(determinant);
    return {std::move(determinant), RischFailure::None};
}

[[nodiscard]] RationalPolynomial reduceModulo(
    const RationalPolynomial& value,
    const RationalPolynomial& modulus) {
    return dividePolynomials(value, modulus).remainder;
}

[[nodiscard]] BivariateRationalPolynomial reduceModulo(
    const BivariateRationalPolynomial& value,
    const RationalPolynomial& modulus) {
    std::vector<RationalPolynomial> result;
    result.reserve(value.coefficientsInX().size());
    for (const RationalPolynomial& coefficient : value.coefficientsInX())
        result.push_back(reduceModulo(coefficient, modulus));
    return BivariateRationalPolynomial{std::move(result)};
}

[[nodiscard]] BivariateRationalPolynomial bivariateOverResidueField(
    const RationalPolynomial& value) {
    std::vector<RationalPolynomial> coefficients;
    coefficients.reserve(value.coefficients().size());
    for (const Rational& coefficient : value.coefficients())
        coefficients.push_back(constantPolynomial(coefficient));
    return BivariateRationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] std::optional<BivariateRationalPolynomial> makeMonicModulo(
    const BivariateRationalPolynomial& value,
    const RationalPolynomial& modulus) {
    BivariateRationalPolynomial reduced = reduceModulo(value, modulus);
    if (reduced.isZero())
        return std::nullopt;
    const auto inverse = inverseModulo(
        reduced.coefficientInX(reduced.degreeInX()), modulus);
    if (!inverse)
        return std::nullopt;
    return reduceModulo(bivariateScale(reduced, *inverse), modulus);
}

[[nodiscard]] bool dividesModulo(
    BivariateRationalPolynomial dividend,
    const BivariateRationalPolynomial& monicDivisor,
    const RationalPolynomial& modulus) {
    dividend = reduceModulo(dividend, modulus);
    if (monicDivisor.isZero())
        return false;
    while (!dividend.isZero()
        && dividend.degreeInX() >= monicDivisor.degreeInX()) {
        const std::size_t shift = dividend.degreeInX() - monicDivisor.degreeInX();
        const RationalPolynomial leading = dividend.coefficientInX(dividend.degreeInX());
        dividend = reduceModulo(
            bivariateSubtract(
                dividend,
                bivariateMultiplyMonomial(monicDivisor, leading, shift)),
            modulus);
    }
    return dividend.isZero();
}

[[nodiscard]] bool isZeroResiduePolynomial(const RationalPolynomial& value) {
    if (value.degree() != 1)
        return false;
    return value.coefficient(0).isZero()
        && value.coefficient(1) == one();
}

[[nodiscard]] RischResult resourceFailure(RischFailure failure) {
    return {RischResultStatus::ResourceLimit, failure, std::nullopt, false};
}

[[nodiscard]] bool isResourceFailure(RischFailure failure) {
    switch (failure) {
    case RischFailure::DegreeLimit:
    case RischFailure::TowerDepthLimit:
    case RischFailure::SubresultantStepLimit:
    case RischFailure::MatrixSizeLimit:
    case RischFailure::CoefficientCountLimit:
    case RischFailure::IntermediateBitLimit:
    case RischFailure::DifferentialOperationLimit:
    case RischFailure::RdeMatrixSizeLimit:
    case RischFailure::RdeStepLimit:
    case RischFailure::ResidueSearchLimit:
        return true;
    default:
        return false;
    }
}

} // namespace

RationalFunction::RationalFunction()
    : numerator(), denominator(onePolynomial()) {}

RationalFunction::RationalFunction(
    RationalPolynomial numeratorValue,
    RationalPolynomial denominatorValue)
    : numerator(std::move(numeratorValue)),
      denominator(std::move(denominatorValue)) {}

RationalPolynomial addRationalPolynomialsExact(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    return addPolynomials(lhs, rhs);
}

RationalPolynomial subtractRationalPolynomialsExact(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    return subtractPolynomials(lhs, rhs);
}

RationalPolynomial negateRationalPolynomialExact(
    const RationalPolynomial& value) {
    return negatePolynomial(value);
}

RationalPolynomial scaleRationalPolynomialExact(
    const RationalPolynomial& value,
    const Rational& scale) {
    return scalePolynomial(value, scale);
}

RationalPolynomial multiplyRationalPolynomialsExact(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    return multiplyPolynomials(lhs, rhs);
}

RationalPolynomial powerRationalPolynomialExact(
    RationalPolynomial base,
    std::size_t exponent) {
    return powerPolynomial(std::move(base), exponent);
}

RationalPolynomial differentiateRationalPolynomialExact(
    const RationalPolynomial& value) {
    return differentiatePolynomial(value);
}

RationalPolynomial integrateRationalPolynomialExact(
    const RationalPolynomial& value) {
    return integratePolynomial(value);
}

RationalPolynomial shiftRationalPolynomialExact(
    const RationalPolynomial& value,
    std::int64_t shift) {
    const RationalPolynomial linear{{Rational{BigInt{shift}}, one()}};
    RationalPolynomial result;
    for (std::size_t index = value.degree() + 1; index-- > 0;) {
        result = addPolynomials(
            multiplyPolynomials(result, linear),
            constantPolynomial(value.coefficient(index)));
    }
    return result;
}

RationalPolynomial monicRationalPolynomialExact(
    const RationalPolynomial& value) {
    return monicPolynomial(value);
}

RationalPolynomial gcdRationalPolynomialsMonic(
    RationalPolynomial lhs,
    RationalPolynomial rhs) {
    return polynomialGcdMonic(std::move(lhs), std::move(rhs));
}

RationalPolynomialDivision divideRationalPolynomials(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator) {
    PolynomialDivision result = dividePolynomials(numerator, denominator);
    return {std::move(result.quotient), std::move(result.remainder)};
}

std::optional<RationalPolynomial> divideRationalPolynomialsExactly(
    const RationalPolynomial& numerator,
    const RationalPolynomial& denominator) {
    RationalPolynomial quotient;
    if (!exactPolynomialQuotient(numerator, denominator, quotient))
        return std::nullopt;
    return quotient;
}

BivariateRationalPolynomial::BivariateRationalPolynomial()
    : coefficientsInX_(1) {}

BivariateRationalPolynomial::BivariateRationalPolynomial(
    std::vector<RationalPolynomial> coefficientsInX)
    : coefficientsInX_(std::move(coefficientsInX)) {
    normalize();
}

bool BivariateRationalPolynomial::isZero() const noexcept {
    return coefficientsInX_.size() == 1 && coefficientsInX_.front().isZero();
}

std::size_t BivariateRationalPolynomial::degreeInX() const noexcept {
    return coefficientsInX_.size() - 1;
}

const RationalPolynomial& BivariateRationalPolynomial::coefficientInX(
    std::size_t exponent) const noexcept {
    static const RationalPolynomial zeroPolynomial;
    return exponent < coefficientsInX_.size()
        ? coefficientsInX_[exponent]
        : zeroPolynomial;
}

const std::vector<RationalPolynomial>&
BivariateRationalPolynomial::coefficientsInX() const noexcept {
    return coefficientsInX_;
}

bool BivariateRationalPolynomial::operator==(
    const BivariateRationalPolynomial& rhs) const noexcept {
    if (coefficientsInX_.size() != rhs.coefficientsInX_.size())
        return false;
    for (std::size_t i = 0; i < coefficientsInX_.size(); ++i)
        if (!samePolynomial(coefficientsInX_[i], rhs.coefficientsInX_[i]))
            return false;
    return true;
}

void BivariateRationalPolynomial::normalize() {
    while (coefficientsInX_.size() > 1 && coefficientsInX_.back().isZero())
        coefficientsInX_.pop_back();
    if (coefficientsInX_.empty())
        coefficientsInX_.emplace_back();
}

RischStageResult<HermitePowerReduction> hermiteReduceSquareFreePower(
    RationalPolynomial numerator,
    const RationalPolynomial& squareFreeFactor,
    std::size_t denominatorPower,
    const RationalPolynomial& inverseDerivativeModuloFactor,
    const RischOptions& options) {
    if (denominatorPower == 0 || squareFreeFactor.degree() == 0)
        return {std::nullopt, RischFailure::InvalidRationalFunction};
    if (squareFreeFactor.degree() > options.maximumHermiteDegree)
        return {std::nullopt, RischFailure::DegreeLimit};
    if (polynomialGcdMonic(
            squareFreeFactor, differentiatePolynomial(squareFreeFactor)).degree() != 0)
        return {std::nullopt, RischFailure::NonSquareFreeDenominator};

    const RationalPolynomial originalNumerator = numerator;
    const std::size_t originalPower = denominatorPower;
    HermitePowerReduction result;
    const RationalPolynomial derivative = differentiatePolynomial(squareFreeFactor);
    while (denominatorPower > 1 && !numerator.isZero()) {
        PolynomialDivision inverseProduct = dividePolynomials(
            multiplyPolynomials(numerator, inverseDerivativeModuloFactor),
            squareFreeFactor);
        const Rational scale = -one()
            / Rational{BigInt::fromUnsigned(denominatorPower - 1)};
        RationalPolynomial correction = scalePolynomial(
            inverseProduct.remainder, scale);

        RationalPolynomial residual = addPolynomials(numerator,
            addPolynomials(
                negatePolynomial(multiplyPolynomials(
                    differentiatePolynomial(correction), squareFreeFactor)),
                scalePolynomial(
                    multiplyPolynomials(correction, derivative),
                    Rational{BigInt::fromUnsigned(denominatorPower - 1)})));
        PolynomialDivision lowered = dividePolynomials(residual, squareFreeFactor);
        if (!lowered.remainder.isZero())
            return {std::nullopt, RischFailure::ExactDivisionFailed};
        const RischFailure failure = polynomialBudgetFailure(lowered.quotient, options);
        if (failure != RischFailure::None)
            return {std::nullopt, failure};

        if (!correction.isZero())
            result.rationalTerms.emplace_back(
                std::move(correction), denominatorPower - 1);
        numerator = std::move(lowered.quotient);
        --denominatorPower;
        ++result.steps;
    }
    result.squareFreeNumerator = std::move(numerator);

    RationalFunction reconstructed{
        result.squareFreeNumerator, squareFreeFactor};
    for (const auto& [correction, power] : result.rationalTerms) {
        const RationalFunction term{
            correction, powerPolynomial(squareFreeFactor, power)};
        reconstructed = addRationalFunctions(
            reconstructed, differentiateRationalFunction(term));
    }
    const RationalFunction original{
        originalNumerator, powerPolynomial(squareFreeFactor, originalPower)};
    result.exactVerified = sameRationalFunction(original, reconstructed);
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

RischStageResult<RationalHermiteReduction> hermiteReduceRationalFunction(
    const RationalFunction& input,
    const RischOptions& options) {
    if (input.denominator.isZero())
        return {std::nullopt, RischFailure::ZeroDenominator};
    RationalFunction normalized = normalizeRationalFunction(input);
    if (normalized.denominator.degree() > options.maximumHermiteDegree)
        return {std::nullopt, RischFailure::DegreeLimit};

    PolynomialDivision proper = dividePolynomials(
        normalized.numerator, normalized.denominator);
    RationalFunction rationalPart{integratePolynomial(proper.quotient), onePolynomial()};
    RationalFunction squareFreePart;
    if (!proper.remainder.isZero()) {
        auto factors = squareFreeDecomposition(normalized.denominator, options);
        if (!factors)
            return {std::nullopt, factors.failure};
        for (const SquareFreeFactor& factor : *factors.value) {
            const RationalPolynomial powered = powerPolynomial(
                factor.polynomial, factor.multiplicity);
            RationalPolynomial cofactor;
            if (!exactPolynomialQuotient(
                    normalized.denominator, powered, cofactor))
                return {std::nullopt, RischFailure::ExactDivisionFailed};
            const auto inverseCofactor = inverseModulo(cofactor, powered);
            if (!inverseCofactor)
                return {std::nullopt, RischFailure::ExactDivisionFailed};
            const RationalPolynomial localNumerator = dividePolynomials(
                multiplyPolynomials(proper.remainder, *inverseCofactor),
                powered).remainder;
            const auto inverseDerivative = inverseModulo(
                differentiatePolynomial(factor.polynomial), factor.polynomial);
            if (!inverseDerivative)
                return {std::nullopt, RischFailure::NonSquareFreeDenominator};
            auto reduced = hermiteReduceSquareFreePower(
                localNumerator, factor.polynomial, factor.multiplicity,
                *inverseDerivative, options);
            if (!reduced)
                return {std::nullopt, reduced.failure};
            for (const auto& [correction, power] : reduced.value->rationalTerms)
                rationalPart = addRationalFunctions(
                    rationalPart,
                    RationalFunction{
                        correction, powerPolynomial(factor.polynomial, power)});
            squareFreePart = addRationalFunctions(
                squareFreePart,
                RationalFunction{
                    reduced.value->squareFreeNumerator, factor.polynomial});
        }
    }

    RationalHermiteReduction result{
        normalizeRationalFunction(std::move(rationalPart)),
        normalizeRationalFunction(std::move(squareFreePart)),
        false};
    result.exactVerified = verifyHermiteReduction(normalized, result);
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

bool verifyHermiteReduction(
    const RationalFunction& input,
    const RationalHermiteReduction& reduction) {
    if (input.denominator.isZero()
        || reduction.rationalPart.denominator.isZero()
        || reduction.squareFreePart.denominator.isZero())
        return false;
    const RationalFunction reconstructed = addRationalFunctions(
        differentiateRationalFunction(reduction.rationalPart),
        reduction.squareFreePart);
    if (!sameRationalFunction(input, reconstructed))
        return false;
    if (!reduction.squareFreePart.numerator.isZero()) {
        if (reduction.squareFreePart.numerator.degree()
            >= reduction.squareFreePart.denominator.degree())
            return false;
        if (polynomialGcdMonic(
                reduction.squareFreePart.numerator,
                reduction.squareFreePart.denominator).degree() != 0)
            return false;
        if (polynomialGcdMonic(
                reduction.squareFreePart.denominator,
                differentiatePolynomial(
                    reduction.squareFreePart.denominator)).degree() != 0)
            return false;
    }
    return true;
}

RischStageResult<LrtResult> lazardRiobooTrager(
    const RationalFunction& properSquareFreePart,
    const RischOptions& options) {
    if (properSquareFreePart.denominator.isZero())
        return {std::nullopt, RischFailure::ZeroDenominator};
    const RationalFunction input = normalizeRationalFunction(properSquareFreePart);
    if (input.numerator.isZero()) {
        LrtResult zeroResult;
        zeroResult.residueResultant = onePolynomial();
        zeroResult.exactVerified = true;
        return {std::move(zeroResult), RischFailure::None};
    }
    if (input.denominator.degree() == 0
        || input.numerator.degree() >= input.denominator.degree())
        return {std::nullopt, RischFailure::InvalidRationalFunction};
    if (input.denominator.degree() > options.maximumLrtDegree)
        return {std::nullopt, RischFailure::DegreeLimit};
    const RationalPolynomial denominatorDerivative =
        differentiatePolynomial(input.denominator);
    if (polynomialGcdMonic(input.denominator, denominatorDerivative).degree() != 0)
        return {std::nullopt, RischFailure::NonSquareFreeDenominator};

    const BivariateRationalPolynomial denominator =
        bivariateOverResidueField(input.denominator);
    std::vector<RationalPolynomial> residueEquationCoefficients(
        input.denominator.degree());
    for (std::size_t i = 0; i < input.denominator.degree(); ++i)
        residueEquationCoefficients[i] = RationalPolynomial{{
            input.numerator.coefficient(i),
            -denominatorDerivative.coefficient(i)}};
    const BivariateRationalPolynomial residueEquation{
        std::move(residueEquationCoefficients)};

    auto sequence = subresultantPolynomialRemainderSequence(
        denominator, residueEquation, options);
    if (!sequence)
        return {std::nullopt, sequence.failure};
    auto rawResultant = resultantByBareiss(
        denominator, residueEquation, options);
    if (!rawResultant)
        return {std::nullopt, rawResultant.failure};
    if (rawResultant.value->isZero())
        return {std::nullopt, RischFailure::CertificateFailed};
    RationalPolynomial resultant = monicPolynomial(*rawResultant.value);
    if (resultant.degree() > options.maximumResidueDegree)
        return {std::nullopt, RischFailure::DegreeLimit};
    const BivariateRationalPolynomial& lastSubresultant =
        sequence.value->back();
    if (lastSubresultant.degreeInX() != 0
        || !samePolynomial(
            monicPolynomial(lastSubresultant.coefficientInX(0)), resultant))
        return {std::nullopt, RischFailure::CertificateFailed};

    auto factors = squareFreeDecomposition(resultant, options);
    if (!factors)
        return {std::nullopt, factors.failure};

    LrtResult result;
    result.residueResultant = resultant;
    result.subresultantSequence = std::move(*sequence.value);
    for (const SquareFreeFactor& factor : *factors.value) {
        if (isZeroResiduePolynomial(factor.polynomial))
            continue;
        const BivariateRationalPolynomial* candidate = nullptr;
        if (factor.multiplicity == input.denominator.degree())
            candidate = &denominator;
        else {
            for (const BivariateRationalPolynomial& member
                 : result.subresultantSequence)
                if (!member.isZero()
                    && member.degreeInX() == factor.multiplicity)
                    candidate = &member;
        }
        if (!candidate)
            return {std::nullopt, RischFailure::CertificateFailed};
        auto monic = makeMonicModulo(*candidate, factor.polynomial);
        if (!monic)
            return {std::nullopt, RischFailure::NonInvertibleLeadingCoefficient};
        AlgebraicResidueLogTerm term{
            factor.polynomial,
            std::move(*monic),
            factor.multiplicity,
            false};
        term.exactVerified = term.logArgument.degreeInX() == term.poleMultiplicity
            && dividesModulo(denominator, term.logArgument, term.residuePolynomial)
            && dividesModulo(residueEquation, term.logArgument, term.residuePolynomial);
        if (!term.exactVerified)
            return {std::nullopt, RischFailure::CertificateFailed};
        result.logarithmicTerms.push_back(std::move(term));
    }
    result.exactVerified = verifyLrtResult(input, result, options);
    if (!result.exactVerified)
        return {std::nullopt, RischFailure::CertificateFailed};
    return {std::move(result), RischFailure::None};
}

bool verifyLrtResult(
    const RationalFunction& properSquareFreePart,
    const LrtResult& result,
    const RischOptions& options) {
    if (properSquareFreePart.denominator.isZero())
        return false;
    const RationalFunction input = normalizeRationalFunction(properSquareFreePart);
    if (input.numerator.isZero())
        return result.logarithmicTerms.empty()
            && samePolynomial(result.residueResultant, onePolynomial());
    if (input.denominator.degree() == 0
        || input.numerator.degree() >= input.denominator.degree()
        || input.denominator.degree() > options.maximumLrtDegree)
        return false;

    const RationalPolynomial derivative = differentiatePolynomial(input.denominator);
    if (polynomialGcdMonic(input.denominator, derivative).degree() != 0)
        return false;
    const BivariateRationalPolynomial denominator =
        bivariateOverResidueField(input.denominator);
    std::vector<RationalPolynomial> residueCoefficients(input.denominator.degree());
    for (std::size_t i = 0; i < input.denominator.degree(); ++i)
        residueCoefficients[i] = RationalPolynomial{{
            input.numerator.coefficient(i), -derivative.coefficient(i)}};
    const BivariateRationalPolynomial residueEquation{
        std::move(residueCoefficients)};
    auto rawResultant = resultantByBareiss(denominator, residueEquation, options);
    if (!rawResultant || rawResultant.value->isZero()
        || !samePolynomial(
            monicPolynomial(*rawResultant.value), result.residueResultant))
        return false;
    if (result.subresultantSequence.empty()
        || result.subresultantSequence.back().degreeInX() != 0
        || !samePolynomial(
            monicPolynomial(
                result.subresultantSequence.back().coefficientInX(0)),
            result.residueResultant))
        return false;

    RationalPolynomial represented = onePolynomial();
    for (const AlgebraicResidueLogTerm& term : result.logarithmicTerms) {
        if (!term.exactVerified || term.poleMultiplicity == 0
            || term.logArgument.degreeInX() != term.poleMultiplicity)
            return false;
        if (polynomialGcdMonic(
                term.residuePolynomial,
                differentiatePolynomial(term.residuePolynomial)).degree() != 0)
            return false;
        const auto monic = makeMonicModulo(
            term.logArgument, term.residuePolynomial);
        if (!monic || !(*monic == term.logArgument)
            || !dividesModulo(denominator, term.logArgument, term.residuePolynomial)
            || !dividesModulo(residueEquation, term.logArgument, term.residuePolynomial))
            return false;
        represented = multiplyPolynomials(
            represented,
            powerPolynomial(term.residuePolynomial, term.poleMultiplicity));
    }
    // 0 residue は integral に寄与しないので certificate 集合から省く。
    PolynomialDivision missing = dividePolynomials(result.residueResultant, represented);
    if (!missing.remainder.isZero())
        return false;
    if (missing.quotient.degree() == 0)
        return true;
    return isZeroResiduePolynomial(monicPolynomial(missing.quotient));
}

RischResult integrateRationalRisch(
    const RationalFunction& input,
    const RischOptions& options) {
    auto hermite = hermiteReduceRationalFunction(input, options);
    if (!hermite) {
        if (isResourceFailure(hermite.failure))
            return resourceFailure(hermite.failure);
        return {
            RischResultStatus::UnsupportedExtension,
            hermite.failure, std::nullopt, false};
    }
    auto logarithmic = lazardRiobooTrager(
        hermite.value->squareFreePart, options);
    if (!logarithmic) {
        if (isResourceFailure(logarithmic.failure))
            return resourceFailure(logarithmic.failure);
        return {
            RischResultStatus::UnsupportedExtension,
            logarithmic.failure, std::nullopt, false};
    }
    RationalRischDecomposition decomposition{
        std::move(*hermite.value), std::move(*logarithmic.value)};
    const bool verified = decomposition.hermite.exactVerified
        && decomposition.logarithmicPart.exactVerified;
    return {
        RischResultStatus::Elementary,
        RischFailure::None,
        std::move(decomposition),
        verified};
}

RationalFunction canonicalizeRationalFunction(RationalFunction value) {
    return normalizeRationalFunction(std::move(value));
}

RationalFunction addRationalFunctionsExact(
    const RationalFunction& lhs,
    const RationalFunction& rhs) {
    return addRationalFunctions(lhs, rhs);
}

RationalFunction subtractRationalFunctionsExact(
    const RationalFunction& lhs,
    const RationalFunction& rhs) {
    return subtractRationalFunctions(lhs, rhs);
}

RationalFunction multiplyRationalFunctionsExact(
    const RationalFunction& lhs,
    const RationalFunction& rhs) {
    return multiplyRationalFunctions(lhs, rhs);
}

std::optional<RationalFunction> divideRationalFunctionsExact(
    const RationalFunction& numerator,
    const RationalFunction& denominator) {
    return divideRationalFunctions(numerator, denominator);
}

RationalFunction differentiateRationalFunctionExact(
    const RationalFunction& value) {
    return differentiateRationalFunction(value);
}

bool equivalentRationalFunctions(
    const RationalFunction& lhs,
    const RationalFunction& rhs) {
    return sameRationalFunction(lhs, rhs);
}

} // namespace mmcal::symbolic::risch
