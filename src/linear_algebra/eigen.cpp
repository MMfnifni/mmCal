// eigenvalue/eigenvector。exact小規模caseと，BigFloat complex Schur numerical backend。
#include "eigen.hpp"

#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/certified_sqrt.hpp"
#include "approximation/expression_interval.hpp"
#include "approximation/interval_math.hpp"
#include "approximation/precision.hpp"
#include "builtins/exact_operations.hpp"
#include "expression/array_utils.hpp"
#include "linear_algebra/complex_point.hpp"
#include "numeric/big_float.hpp"
#include "numeric/big_int.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/number.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra {
namespace {

using approximation::ComplexInterval;
using approximation::RealInterval;
using detail::ComplexPoint;
using expression::ArrayExpr;
using expression::Expr;
using numeric::BigFloat;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RoundingMode;

constexpr std::size_t maximumPrecisionRetries = 10;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr add(Expr lhs, Expr rhs, const ExactMatrixContext& context) {
    return builtins::exact::add({std::move(lhs), std::move(rhs)},
        context.builtins, context.mathematics, context.angles);
}
[[nodiscard]] Expr subtract(Expr lhs, Expr rhs, const ExactMatrixContext& context) {
    return builtins::exact::subtract(std::move(lhs), std::move(rhs),
        context.builtins, context.mathematics, context.angles);
}
[[nodiscard]] Expr multiply(Expr lhs, Expr rhs, const ExactMatrixContext& context) {
    return builtins::exact::multiply({std::move(lhs), std::move(rhs)},
        context.builtins, context.mathematics, context.angles);
}
[[nodiscard]] Expr divide(Expr lhs, Expr rhs, const ExactMatrixContext& context) {
    return builtins::exact::divide(std::move(lhs), std::move(rhs),
        context.builtins, context.mathematics, context.angles);
}
[[nodiscard]] Expr square(Expr value, const ExactMatrixContext& context) {
    return multiply(value, value, context);
}
[[nodiscard]] Expr squareRoot(Expr value, const ExactMatrixContext& context) {
    return builtins::exact::sqrt(std::move(value),
        context.builtins, context.mathematics, context.angles);
}

[[nodiscard]] bool isUpperTriangular(const MatrixView& matrix) {
    for (std::size_t row = 1; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < row; ++column) {
            const Expr& value = matrix(row, column);
            if (!value.isNumber() || !value.asNumber().isZero())
                return false;
        }
    return true;
}

[[nodiscard]] bool isDiagonal(const MatrixView& matrix) {
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < matrix.columns(); ++column) {
            if (row == column)
                continue;
            const Expr& value = matrix(row, column);
            if (!value.isNumber() || !value.asNumber().isZero())
                return false;
        }
    return true;
}

[[nodiscard]] Expr diagonalValues(const MatrixView& matrix) {
    std::vector<Expr> values;
    values.reserve(matrix.rows());
    for (std::size_t i = 0; i < matrix.rows(); ++i)
        values.push_back(matrix(i, i));
    return Expr::array({matrix.rows()}, std::move(values));
}

[[nodiscard]] Expr identityExpr(std::size_t size) {
    std::vector<Expr> values(size * size, integer(0));
    for (std::size_t i = 0; i < size; ++i)
        values[i * size + i] = integer(1);
    return Expr::array({size, size}, std::move(values));
}

[[nodiscard]] Expr exactTwoByTwoEigenvalues(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    const Expr trace = add(matrix(0, 0), matrix(1, 1), context);
    const Expr diagonalDifference = subtract(matrix(0, 0), matrix(1, 1), context);
    const Expr discriminant = add(
        square(diagonalDifference, context),
        multiply(integer(4), multiply(matrix(0, 1), matrix(1, 0), context), context),
        context);
    const Expr root = squareRoot(discriminant, context);
    const Expr two = integer(2);
    return Expr::array({2}, {
        divide(add(trace, root, context), two, context),
        divide(subtract(trace, root, context), two, context)
    });
}

[[nodiscard]] bool allNumbers(const MatrixView& matrix) {
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < matrix.columns(); ++column)
            if (!matrix(row, column).isNumber())
                return false;
    return true;
}

[[nodiscard]] std::optional<Expr> exactTwoByTwoEigenvectors(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    if (!allNumbers(matrix))
        return std::nullopt;
    const Expr values = exactTwoByTwoEigenvalues(matrix, context);
    const auto& lambdas = values.asArray().elements;
    if (lambdas[0] == lambdas[1])
        return std::nullopt;

    const bool useUpper = !matrix(0, 1).asNumber().isZero();
    const bool useLower = !matrix(1, 0).asNumber().isZero();
    if (!useUpper && !useLower)
        return identityExpr(2);

    std::vector<Expr> result;
    result.reserve(4);
    if (useUpper) {
        const Expr b = matrix(0, 1);
        result.push_back(b);
        result.push_back(b);
        result.push_back(subtract(lambdas[0], matrix(0, 0), context));
        result.push_back(subtract(lambdas[1], matrix(0, 0), context));
    }
    else {
        result.push_back(subtract(lambdas[0], matrix(1, 1), context));
        result.push_back(subtract(lambdas[1], matrix(1, 1), context));
        const Expr c = matrix(1, 0);
        result.push_back(c);
        result.push_back(c);
    }
    return Expr::array({2, 2}, std::move(result));
}

[[nodiscard]] BigFloat zero(std::size_t bits) {
    return BigFloat::fromBigInt(BigInt{}, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat one(std::size_t bits) {
    return BigFloat::fromBigInt(BigInt{1}, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat two(std::size_t bits) {
    return BigFloat::fromBigInt(BigInt{2}, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat add(const BigFloat& lhs, const BigFloat& rhs, std::size_t bits) {
    return numeric::add(lhs, rhs, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat subtract(const BigFloat& lhs, const BigFloat& rhs, std::size_t bits) {
    return numeric::subtract(lhs, rhs, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat multiply(const BigFloat& lhs, const BigFloat& rhs, std::size_t bits) {
    return numeric::multiply(lhs, rhs, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat divide(const BigFloat& lhs, const BigFloat& rhs, std::size_t bits) {
    return numeric::divide(lhs, rhs, bits, RoundingMode::NearestEven);
}
[[nodiscard]] BigFloat absolute(const BigFloat& value) {
    return value.isNegative() ? -value : value;
}

[[nodiscard]] Rational midpointRational(const RealInterval& interval) {
    return (interval.lower().toRational() + interval.upper().toRational()) / Rational{BigInt{2}};
}
[[nodiscard]] BigFloat midpoint(const RealInterval& interval, std::size_t bits) {
    return BigFloat::fromRational(midpointRational(interval), bits, RoundingMode::NearestEven);
}
[[nodiscard]] ComplexPoint midpoint(const ComplexInterval& interval, std::size_t bits) {
    return ComplexPoint{midpoint(interval.real(), bits), midpoint(interval.imaginary(), bits)};
}

[[nodiscard]] BigFloat pointSqrt(const BigFloat& value, std::size_t bits) {
    if (value.isNegative())
        throw approximation::PrecisionInsufficient("Eigen iteration encountered a negative rounded norm");
    const auto enclosure = approximation::encloseSqrt(value.toRational(), bits + 8).interval;
    return midpoint(enclosure, bits);
}

[[nodiscard]] ComplexPoint complexSquareRoot(const ComplexPoint& value, std::size_t bits) {
    if (detail::complexIsZero(value))
        return detail::complexZero(bits);
    const BigFloat magnitude = detail::complexMagnitude(value, bits);
    BigFloat realPart = divide(add(magnitude, value.real, bits), two(bits), bits);
    BigFloat imaginaryPart = divide(subtract(magnitude, value.real, bits), two(bits), bits);
    if (realPart.isNegative())
        realPart = zero(bits);
    if (imaginaryPart.isNegative())
        imaginaryPart = zero(bits);
    BigFloat x = pointSqrt(realPart, bits);
    BigFloat y = pointSqrt(imaginaryPart, bits);
    if (value.imaginary.isNegative())
        y = -y;
    return ComplexPoint{std::move(x), std::move(y)};
}

class ComplexMatrix final {
public:
    ComplexMatrix(std::size_t rows, std::size_t columns, std::size_t bits)
        : rows_(rows), columns_(columns), values_(rows * columns, detail::complexZero(bits)) {}
    ComplexMatrix(std::size_t rows, std::size_t columns, std::vector<ComplexPoint> values)
        : rows_(rows), columns_(columns), values_(std::move(values)) {
        if (values_.size() != rows_ * columns_)
            throw std::invalid_argument("ComplexMatrix shape does not match element count");
    }
    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] ComplexPoint& operator()(std::size_t row, std::size_t column) noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const ComplexPoint& operator()(std::size_t row, std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const std::vector<ComplexPoint>& values() const noexcept { return values_; }
private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<ComplexPoint> values_;
};

[[nodiscard]] ComplexMatrix identity(std::size_t size, std::size_t bits) {
    ComplexMatrix result{size, size, bits};
    for (std::size_t i = 0; i < size; ++i)
        result(i, i) = detail::complexOne(bits);
    return result;
}

struct Reflector final {
    std::size_t first = 0;
    std::vector<ComplexPoint> vector;
    BigFloat beta;
};

[[nodiscard]] std::optional<Reflector> reflectorForValues(
    std::vector<ComplexPoint> values,
    std::size_t first,
    std::size_t bits) {
    BigFloat normSquared = zero(bits);
    for (const ComplexPoint& value : values)
        normSquared = add(normSquared, detail::complexMagnitudeSquared(value, bits), bits);
    if (normSquared.isZero())
        return std::nullopt;

    const BigFloat norm = pointSqrt(normSquared, bits);
    ComplexPoint phase = detail::complexOne(bits);
    const BigFloat firstMagnitudeSquared = detail::complexMagnitudeSquared(values.front(), bits);
    if (!firstMagnitudeSquared.isZero())
        phase = detail::complexDivideReal(
            values.front(), pointSqrt(firstMagnitudeSquared, bits), bits);
    const ComplexPoint alpha = detail::complexNegate(detail::complexScale(phase, norm, bits));
    values.front() = detail::complexSubtract(values.front(), alpha, bits);

    BigFloat vNormSquared = zero(bits);
    for (const ComplexPoint& value : values)
        vNormSquared = add(vNormSquared, detail::complexMagnitudeSquared(value, bits), bits);
    if (vNormSquared.isZero())
        return std::nullopt;
    return Reflector{first, std::move(values), divide(two(bits), vNormSquared, bits)};
}

void applyLeft(
    ComplexMatrix& matrix,
    const Reflector& reflector,
    std::size_t firstColumn,
    std::size_t lastColumn,
    std::size_t bits) {
    for (std::size_t column = firstColumn; column <= lastColumn; ++column) {
        ComplexPoint dot = detail::complexZero(bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            dot = detail::complexAdd(dot,
                detail::complexMultiply(detail::complexConjugate(reflector.vector[i]),
                    matrix(reflector.first + i, column), bits), bits);
        const ComplexPoint scaled = detail::complexScale(dot, reflector.beta, bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            matrix(reflector.first + i, column) = detail::complexSubtract(
                matrix(reflector.first + i, column),
                detail::complexMultiply(reflector.vector[i], scaled, bits), bits);
    }
}

void applyRight(
    ComplexMatrix& matrix,
    const Reflector& reflector,
    std::size_t firstRow,
    std::size_t lastRow,
    std::size_t bits) {
    for (std::size_t row = firstRow; row <= lastRow; ++row) {
        ComplexPoint dot = detail::complexZero(bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            dot = detail::complexAdd(dot,
                detail::complexMultiply(matrix(row, reflector.first + i), reflector.vector[i], bits), bits);
        const ComplexPoint scaled = detail::complexScale(dot, reflector.beta, bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            matrix(row, reflector.first + i) = detail::complexSubtract(
                matrix(row, reflector.first + i),
                detail::complexMultiply(scaled, detail::complexConjugate(reflector.vector[i]), bits), bits);
    }
}

void hessenbergReduce(
    ComplexMatrix& matrix,
    ComplexMatrix* vectors,
    std::size_t bits) {
    const std::size_t size = matrix.rows();
    if (size < 3)
        return;
    for (std::size_t column = 0; column + 2 < size; ++column) {
        std::vector<ComplexPoint> values;
        values.reserve(size - column - 1);
        for (std::size_t row = column + 1; row < size; ++row)
            values.push_back(matrix(row, column));
        const auto reflector = reflectorForValues(std::move(values), column + 1, bits);
        if (!reflector)
            continue;
        applyLeft(matrix, *reflector, column, size - 1, bits);
        applyRight(matrix, *reflector, 0, size - 1, bits);
        if (vectors)
            applyRight(*vectors, *reflector, 0, size - 1, bits);
        for (std::size_t row = column + 2; row < size; ++row)
            matrix(row, column) = detail::complexZero(bits);
    }
}

[[nodiscard]] BigFloat tolerance(std::size_t targetBits, std::size_t bits) {
    const auto exponent = -static_cast<BigFloat::exponent_type>(
        std::min<std::size_t>(targetBits,
            static_cast<std::size_t>(std::numeric_limits<BigFloat::exponent_type>::max() / 2)));
    return BigFloat::fromDyadic(BigInt{1}, exponent, bits, RoundingMode::TowardPositive);
}

[[nodiscard]] bool smallSubdiagonal(
    const ComplexMatrix& matrix,
    std::size_t row,
    const BigFloat& epsilon,
    std::size_t bits) {
    const BigFloat sub = detail::complexMagnitude(matrix(row, row - 1), bits);
    BigFloat scale = add(
        detail::complexMagnitude(matrix(row - 1, row - 1), bits),
        detail::complexMagnitude(matrix(row, row), bits), bits);
    scale = add(scale, one(bits), bits);
    return sub <= multiply(epsilon, scale, bits);
}

[[nodiscard]] ComplexPoint wilkinsonShift(
    const ComplexMatrix& matrix,
    std::size_t high,
    std::size_t bits) {
    if (high == 0)
        return matrix(0, 0);
    const ComplexPoint a = matrix(high - 1, high - 1);
    const ComplexPoint b = matrix(high - 1, high);
    const ComplexPoint c = matrix(high, high - 1);
    const ComplexPoint d = matrix(high, high);
    const ComplexPoint halfTrace = detail::complexScale(
        detail::complexAdd(a, d, bits), divide(one(bits), two(bits), bits), bits);
    const ComplexPoint halfDifference = detail::complexScale(
        detail::complexSubtract(a, d, bits), divide(one(bits), two(bits), bits), bits);
    const ComplexPoint radicand = detail::complexAdd(
        detail::complexMultiply(halfDifference, halfDifference, bits),
        detail::complexMultiply(b, c, bits), bits);
    const ComplexPoint root = complexSquareRoot(radicand, bits);
    const ComplexPoint first = detail::complexAdd(halfTrace, root, bits);
    const ComplexPoint second = detail::complexSubtract(halfTrace, root, bits);
    const BigFloat firstDistance = detail::complexMagnitude(
        detail::complexSubtract(first, d, bits), bits);
    const BigFloat secondDistance = detail::complexMagnitude(
        detail::complexSubtract(second, d, bits), bits);
    return firstDistance <= secondDistance ? first : second;
}

void implicitShiftSweep(
    ComplexMatrix& matrix,
    ComplexMatrix* vectors,
    std::size_t low,
    std::size_t high,
    const ComplexPoint& shift,
    std::size_t bits) {
    for (std::size_t k = low; k < high; ++k) {
        ComplexPoint x;
        ComplexPoint y;
        if (k == low) {
            x = detail::complexSubtract(matrix(k, k), shift, bits);
            y = matrix(k + 1, k);
        }
        else {
            x = matrix(k, k - 1);
            y = matrix(k + 1, k - 1);
        }
        std::vector<ComplexPoint> pair{x, y};
        const auto reflector = reflectorForValues(std::move(pair), k, bits);
        if (!reflector)
            continue;
        const std::size_t firstColumn = k == low ? low : k - 1;
        applyLeft(matrix, *reflector, firstColumn, matrix.columns() - 1, bits);
        applyRight(matrix, *reflector, 0, std::min(matrix.rows() - 1, k + 2), bits);
        if (vectors)
            applyRight(*vectors, *reflector, 0, vectors->rows() - 1, bits);
        if (k > low)
            matrix(k + 1, k - 1) = detail::complexZero(bits);
    }
}

struct SchurResult final {
    ComplexMatrix triangular;
    std::optional<ComplexMatrix> vectors;
};

[[nodiscard]] SchurResult complexSchur(
    ComplexMatrix source,
    std::size_t targetBits,
    std::size_t bits,
    bool needVectors) {
    const std::size_t size = source.rows();
    std::optional<ComplexMatrix> vectors;
    if (needVectors)
        vectors.emplace(identity(size, bits));
    hessenbergReduce(source, vectors ? &*vectors : nullptr, bits);
    if (size < 2)
        return SchurResult{std::move(source), std::move(vectors)};

    const BigFloat epsilon = tolerance(targetBits, bits);
    std::size_t high = size - 1;
    std::size_t sweepsWithoutDeflation = 0;
    const std::size_t maximumSweeps = std::max<std::size_t>(256, 96 * size);
    while (high > 0) {
        if (smallSubdiagonal(source, high, epsilon, bits)) {
            source(high, high - 1) = detail::complexZero(bits);
            --high;
            sweepsWithoutDeflation = 0;
            continue;
        }
        std::size_t low = high;
        while (low > 0 && !smallSubdiagonal(source, low, epsilon, bits))
            --low;
        if (low > 0)
            source(low, low - 1) = detail::complexZero(bits);

        ComplexPoint shift = wilkinsonShift(source, high, bits);
        if (sweepsWithoutDeflation != 0 && sweepsWithoutDeflation % 24 == 0) {
            // exceptional shiftで稀な停滞を外す。doubleへ落とさずBigFloatのまま。
            shift = detail::complexAdd(shift,
                detail::complexScale(source(high, high - 1), divide(one(bits), two(bits), bits), bits), bits);
        }
        implicitShiftSweep(source, vectors ? &*vectors : nullptr, low, high, shift, bits);
        ++sweepsWithoutDeflation;
        if (sweepsWithoutDeflation > maximumSweeps)
            throw approximation::PrecisionInsufficient(
                "Eigen QR iteration did not converge at the current precision");
    }
    return SchurResult{std::move(source), std::move(vectors)};
}

[[nodiscard]] std::optional<std::vector<ComplexInterval>> encloseElements(
    const ArrayExpr& array,
    std::size_t bits,
    const approximation::CertifiedEvaluator& certified) {
    std::vector<ComplexInterval> values;
    values.reserve(array.elements.size());
    for (const Expr& element : array.elements) {
        const auto enclosed = approximation::encloseComplexExpression(element, bits, certified);
        if (!enclosed)
            return std::nullopt;
        values.push_back(*enclosed);
    }
    return values;
}

[[nodiscard]] ComplexMatrix midpointMatrix(
    std::size_t rows,
    std::size_t columns,
    const std::vector<ComplexInterval>& values,
    std::size_t bits) {
    std::vector<ComplexPoint> points;
    points.reserve(values.size());
    for (const ComplexInterval& value : values)
        points.push_back(midpoint(value, bits));
    return ComplexMatrix{rows, columns, std::move(points)};
}

[[nodiscard]] ComplexMatrix multiplyMatrices(
    const ComplexMatrix& lhs,
    const ComplexMatrix& rhs,
    std::size_t bits) {
    ComplexMatrix result{lhs.rows(), rhs.columns(), bits};
    for (std::size_t row = 0; row < lhs.rows(); ++row)
        for (std::size_t column = 0; column < rhs.columns(); ++column) {
            ComplexPoint sum = detail::complexZero(bits);
            for (std::size_t k = 0; k < lhs.columns(); ++k)
                sum = detail::complexAdd(sum,
                    detail::complexMultiply(lhs(row, k), rhs(k, column), bits), bits);
            result(row, column) = std::move(sum);
        }
    return result;
}

[[nodiscard]] ComplexMatrix conjugateTranspose(const ComplexMatrix& matrix, std::size_t bits) {
    ComplexMatrix result{matrix.columns(), matrix.rows(), bits};
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < matrix.columns(); ++column)
            result(column, row) = detail::complexConjugate(matrix(row, column));
    return result;
}

[[nodiscard]] BigFloat matrixDifferenceMax(
    const ComplexMatrix& lhs,
    const ComplexMatrix& rhs,
    std::size_t bits) {
    BigFloat result = zero(bits);
    for (std::size_t row = 0; row < lhs.rows(); ++row)
        for (std::size_t column = 0; column < lhs.columns(); ++column) {
            const BigFloat magnitude = detail::complexMagnitude(
                detail::complexSubtract(lhs(row, column), rhs(row, column), bits), bits);
            if (magnitude > result)
                result = magnitude;
        }
    return result;
}

[[nodiscard]] ComplexInterval pointInterval(
    const ComplexPoint& value) {
    return ComplexInterval{
        RealInterval::point(value.real),
        RealInterval::point(value.imaginary)
    };
}

[[nodiscard]] BigFloat intervalComponentMagnitudeUpper(
    const ComplexInterval& value,
    std::size_t bits) {
    const RealInterval realAbsolute = approximation::absoluteInterval(value.real(), bits);
    const RealInterval imaginaryAbsolute = approximation::absoluteInterval(value.imaginary(), bits);
    return realAbsolute.upper() >= imaginaryAbsolute.upper()
        ? realAbsolute.upper()
        : imaginaryAbsolute.upper();
}

[[nodiscard]] BigFloat schurRelationResidualUpper(
    const std::vector<ComplexInterval>& source,
    std::size_t size,
    const ComplexMatrix& vectors,
    const ComplexMatrix& triangular,
    std::size_t bits) {
    BigFloat residual = zero(bits);
    for (std::size_t row = 0; row < size; ++row)
        for (std::size_t column = 0; column < size; ++column) {
            ComplexInterval aq = ComplexInterval::fromReal(
                RealInterval::fromRational(Rational{BigInt{0}}, bits));
            ComplexInterval qt = aq;
            for (std::size_t k = 0; k < size; ++k) {
                aq = approximation::add(aq,
                    approximation::multiply(source[row * size + k],
                        pointInterval(vectors(k, column)), bits), bits);
                qt = approximation::add(qt,
                    approximation::multiply(pointInterval(vectors(row, k)),
                        pointInterval(triangular(k, column)), bits), bits);
            }
            const ComplexInterval difference = approximation::subtract(aq, qt, bits);
            const BigFloat current = intervalComponentMagnitudeUpper(difference, bits);
            if (current > residual)
                residual = current;
        }
    return residual;
}

[[nodiscard]] BigFloat decimalTolerance(std::size_t digits, std::size_t bits) {
    if (digits > std::numeric_limits<std::uint64_t>::max() - 3)
        throw std::overflow_error("Eigen requested decimal precision is too large");
    const BigInt denominator = numeric::pow(BigInt{10}, static_cast<std::uint64_t>(digits + 3));
    return BigFloat::fromRational(Rational{BigInt{1}, denominator}, bits, RoundingMode::TowardPositive);
}

[[nodiscard]] bool verifySchur(
    const std::vector<ComplexInterval>& source,
    const SchurResult& result,
    std::size_t digits,
    std::size_t bits) {
    if (!result.vectors)
        return true;
    const BigFloat residual = schurRelationResidualUpper(
        source, result.triangular.rows(), *result.vectors, result.triangular, bits);
    const ComplexMatrix gram = multiplyMatrices(
        conjugateTranspose(*result.vectors, bits), *result.vectors, bits);
    const BigFloat orthogonality = matrixDifferenceMax(gram, identity(gram.rows(), bits), bits);
    const BigFloat allowed = decimalTolerance(digits, bits);
    return residual <= allowed && orthogonality <= allowed;
}

[[nodiscard]] bool closeToZero(
    const ComplexPoint& value,
    const ComplexPoint& scale,
    std::size_t targetBits,
    std::size_t bits) {
    BigFloat threshold = multiply(tolerance(targetBits, bits),
        add(detail::complexMagnitude(scale, bits), one(bits), bits), bits);
    return detail::complexMagnitude(value, bits) <= threshold;
}

[[nodiscard]] std::optional<ComplexMatrix> schurEigenvectors(
    const SchurResult& schur,
    std::size_t targetBits,
    std::size_t bits) {
    if (!schur.vectors)
        return std::nullopt;
    const std::size_t size = schur.triangular.rows();
    ComplexMatrix vectors{size, size, bits};
    for (std::size_t eigenIndex = 0; eigenIndex < size; ++eigenIndex) {
        const ComplexPoint lambda = schur.triangular(eigenIndex, eigenIndex);
        std::vector<ComplexPoint> x(size, detail::complexZero(bits));
        x[eigenIndex] = detail::complexOne(bits);
        for (std::size_t offset = 0; offset < eigenIndex; ++offset) {
            const std::size_t row = eigenIndex - 1 - offset;
            ComplexPoint sum = detail::complexZero(bits);
            for (std::size_t column = row + 1; column <= eigenIndex; ++column)
                sum = detail::complexAdd(sum,
                    detail::complexMultiply(schur.triangular(row, column), x[column], bits), bits);
            const ComplexPoint denominator = detail::complexSubtract(
                schur.triangular(row, row), lambda, bits);
            if (closeToZero(denominator, lambda, targetBits / 2 + 4, bits))
                return std::nullopt; // repeated/near-defective caseは無理にvectorを捏造しない。
            x[row] = detail::complexNegate(detail::complexDivide(sum, denominator, bits));
        }

        BigFloat normSquared = zero(bits);
        for (const ComplexPoint& value : x)
            normSquared = add(normSquared, detail::complexMagnitudeSquared(value, bits), bits);
        if (normSquared.isZero())
            return std::nullopt;
        const BigFloat norm = pointSqrt(normSquared, bits);
        for (ComplexPoint& value : x)
            value = detail::complexDivideReal(value, norm, bits);

        std::vector<ComplexPoint> v(size, detail::complexZero(bits));
        for (std::size_t row = 0; row < size; ++row)
            for (std::size_t column = 0; column < size; ++column)
                v[row] = detail::complexAdd(v[row],
                    detail::complexMultiply((*schur.vectors)(row, column), x[column], bits), bits);

        std::size_t pivot = 0;
        BigFloat pivotMagnitude = zero(bits);
        for (std::size_t row = 0; row < size; ++row) {
            const BigFloat magnitude = detail::complexMagnitudeSquared(v[row], bits);
            if (magnitude > pivotMagnitude) {
                pivotMagnitude = magnitude;
                pivot = row;
            }
        }
        if (!pivotMagnitude.isZero()) {
            const ComplexPoint factor = detail::complexDivideReal(
                detail::complexConjugate(v[pivot]), pointSqrt(pivotMagnitude, bits), bits);
            for (ComplexPoint& value : v)
                value = detail::complexMultiply(value, factor, bits);
        }
        for (std::size_t row = 0; row < size; ++row)
            vectors(row, eigenIndex) = std::move(v[row]);
    }
    return vectors;
}

[[nodiscard]] bool verifyEigenpairs(
    const std::vector<ComplexInterval>& source,
    std::size_t size,
    const ComplexMatrix& vectors,
    const ComplexMatrix& triangular,
    std::size_t digits,
    std::size_t bits) {
    BigFloat residual = zero(bits);
    for (std::size_t column = 0; column < size; ++column) {
        const ComplexPoint lambda = triangular(column, column);
        BigFloat normSquared = zero(bits);
        for (std::size_t row = 0; row < size; ++row) {
            ComplexInterval av = ComplexInterval::fromReal(
                RealInterval::fromRational(Rational{BigInt{0}}, bits));
            for (std::size_t k = 0; k < size; ++k)
                av = approximation::add(av,
                    approximation::multiply(source[row * size + k],
                        pointInterval(vectors(k, column)), bits), bits);
            const ComplexInterval lv = approximation::multiply(
                pointInterval(lambda), pointInterval(vectors(row, column)), bits);
            const BigFloat current = intervalComponentMagnitudeUpper(
                approximation::subtract(av, lv, bits), bits);
            if (current > residual)
                residual = current;
            normSquared = add(normSquared,
                detail::complexMagnitudeSquared(vectors(row, column), bits), bits);
        }
        const BigFloat normResidual = absolute(subtract(normSquared, one(bits), bits));
        if (normResidual > residual)
            residual = normResidual;
    }
    return residual <= decimalTolerance(digits, bits);
}

[[nodiscard]] std::optional<Expr> decimalPoint(const ComplexPoint& value, std::size_t digits) {
    const Rational real = value.real.toRational();
    const Rational imaginary = value.imaginary.toRational();
    const auto realDecimal = numeric::DecimalApproximation::fromCertifiedInterval(real, real, digits);
    const auto imaginaryDecimal = numeric::DecimalApproximation::fromCertifiedInterval(imaginary, imaginary, digits);
    if (!realDecimal || !imaginaryDecimal)
        return std::nullopt;
    if (imaginary.isZero())
        return Expr{*realDecimal};
    return Expr{numeric::ComplexDecimalApproximation::fromComponents(
        *realDecimal, *imaginaryDecimal, real.isZero(), false)};
}

[[nodiscard]] std::optional<Expr> decimalEigenvalues(
    const ComplexMatrix& triangular,
    std::size_t digits) {
    std::vector<Expr> values;
    values.reserve(triangular.rows());
    for (std::size_t i = 0; i < triangular.rows(); ++i) {
        const auto value = decimalPoint(triangular(i, i), digits);
        if (!value)
            return std::nullopt;
        values.push_back(*value);
    }
    return Expr::array({triangular.rows()}, std::move(values));
}

[[nodiscard]] std::optional<Expr> decimalMatrix(
    const ComplexMatrix& matrix,
    std::size_t digits) {
    std::vector<Expr> values;
    values.reserve(matrix.values().size());
    for (const ComplexPoint& point : matrix.values()) {
        const auto value = decimalPoint(point, digits);
        if (!value)
            return std::nullopt;
        values.push_back(*value);
    }
    return Expr::array({matrix.rows(), matrix.columns()}, std::move(values));
}

enum class ApproximateEigenOutput {
    Values,
    Vectors,
    System
};

[[nodiscard]] std::optional<Expr> approximateAtPrecision(
    const ArrayExpr& source,
    std::size_t bits,
    std::size_t digits,
    const approximation::CertifiedEvaluator& certified,
    ApproximateEigenOutput output) {
    const auto enclosed = encloseElements(source, bits, certified);
    if (!enclosed)
        return std::nullopt;
    const std::size_t size = source.shape[0];
    const ComplexMatrix midpointSource = midpointMatrix(size, size, *enclosed, bits);
    const std::size_t targetBits = approximation::decimalDigitsToBinaryBits(digits + 10);
    SchurResult schur = complexSchur(midpointSource, targetBits, bits, true);
    if (!verifySchur(*enclosed, schur, digits, bits))
        throw approximation::PrecisionInsufficient(
            "Eigen Schur relation is not certified at the current precision");
    const auto values = decimalEigenvalues(schur.triangular, digits);
    if (!values)
        return std::nullopt;
    if (output == ApproximateEigenOutput::Values)
        return values;

    const auto vectors = schurEigenvectors(schur, targetBits, bits);
    if (!vectors)
        return std::nullopt;
    if (!verifyEigenpairs(*enclosed, size, *vectors, schur.triangular, digits, bits))
        throw approximation::PrecisionInsufficient(
            "Eigenpair residual is not certified at the current precision");
    const auto vectorExpr = decimalMatrix(*vectors, digits);
    if (!vectorExpr)
        return std::nullopt;
    if (output == ApproximateEigenOutput::Vectors)
        return vectorExpr;
    return expression::braceValue({*values, *vectorExpr});
}

[[nodiscard]] std::optional<Expr> approximateEigen(
    const ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context,
    ApproximateEigenOutput output) {
    if (!matrix.isMatrix() || matrix.shape[0] != matrix.shape[1])
        return std::nullopt;
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    for (std::size_t attempt = 0; attempt < maximumPrecisionRetries; ++attempt) {
        try {
            if (const auto result = approximateAtPrecision(
                matrix, context.workingBinaryBits(), context.decimalDigits(), certified, output))
                return result;
        }
        catch (const approximation::PrecisionInsufficient&) {
        }
        context.setGuardDigits(approximation::nextGuardDigits(context.guardDigits()));
    }
    return std::nullopt;
}

} // namespace

std::optional<Expr> eigenvalues(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    if (matrix.rows() != matrix.columns())
        return std::nullopt;
    if (isUpperTriangular(matrix))
        return diagonalValues(matrix);
    if (matrix.rows() == 2)
        return exactTwoByTwoEigenvalues(matrix, context);
    return std::nullopt;
}

std::optional<Expr> eigenvectors(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    if (matrix.rows() != matrix.columns())
        return std::nullopt;
    if (isDiagonal(matrix))
        return identityExpr(matrix.rows());
    if (matrix.rows() == 2)
        return exactTwoByTwoEigenvectors(matrix, context);
    return std::nullopt;
}

std::optional<Expr> eigensystem(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    if (matrix.rows() != matrix.columns())
        return std::nullopt;
    if (isDiagonal(matrix)) {
        const Expr values = diagonalValues(matrix);
        const Expr vectors = identityExpr(matrix.rows());
        return expression::braceValue({values, vectors});
    }
    if (matrix.rows() == 2) {
        const auto vectors = exactTwoByTwoEigenvectors(matrix, context);
        if (!vectors)
            return std::nullopt;
        return expression::braceValue({exactTwoByTwoEigenvalues(matrix, context), *vectors});
    }
    return std::nullopt;
}

std::optional<Expr> approximateEigenvalues(
    const ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return approximateEigen(matrix, builtins, mathematics, angles, std::move(context),
        ApproximateEigenOutput::Values);
}

std::optional<Expr> approximateEigenvectors(
    const ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return approximateEigen(matrix, builtins, mathematics, angles, std::move(context),
        ApproximateEigenOutput::Vectors);
}

std::optional<Expr> approximateEigensystem(
    const ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return approximateEigen(matrix, builtins, mathematics, angles, std::move(context),
        ApproximateEigenOutput::System);
}

} // namespace mmcal::linear_algebra
