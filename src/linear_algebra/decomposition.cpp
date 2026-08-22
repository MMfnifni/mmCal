// LU / Householder QR。exactとcertified approximateの分解kernelを一箇所にまとめる。
#include "decomposition.hpp"

#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/certified_sqrt.hpp"
#include "approximation/expression_interval.hpp"
#include "approximation/interval_math.hpp"
#include "builtins/exact_operations.hpp"
#include "expression/array_utils.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "mathematics/exact_roots.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <span>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra {
namespace {

using approximation::ComplexInterval;
using approximation::RealInterval;
using expression::ArrayExpr;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

constexpr std::size_t maximumPrecisionRetries = 12;
constexpr std::size_t defaultQrBlockColumns = 1;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

class NumberMatrix final {
public:
    NumberMatrix(std::size_t rows, std::size_t columns)
        : rows_(rows), columns_(columns) {
        const std::size_t shape[] = {rows, columns};
        values_.assign(expression::arrayElementCount(shape), Number{BigInt{0}});
    }

    explicit NumberMatrix(const MatrixView& source)
        : NumberMatrix(source.rows(), source.columns()) {
        for (std::size_t i = 0; i < values_.size(); ++i)
            values_[i] = source.array().exactNumber(i);
    }

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] Number& operator()(std::size_t row, std::size_t column) noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const Number& operator()(std::size_t row, std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }
    void swapRows(std::size_t lhs, std::size_t rhs) noexcept {
        if (lhs == rhs)
            return;
        for (std::size_t column = 0; column < columns_; ++column)
            std::swap((*this)(lhs, column), (*this)(rhs, column));
    }

    [[nodiscard]] const std::vector<Number>& values() const noexcept { return values_; }

private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<Number> values_;
};

[[nodiscard]] bool allExactRealNumbers(const MatrixView& matrix) noexcept {
    return matrix.array().hasExactRealStorage();
}

[[nodiscard]] Expr packedNumberMatrices(
    std::span<const NumberMatrix> matrices,
    std::size_t rows,
    std::size_t columns) {
    std::vector<Expr> values;
    const std::size_t shape[] = {matrices.size(), rows, columns};
    values.reserve(expression::arrayElementCount(shape));
    for (const NumberMatrix& matrix : matrices)
        for (const Number& value : matrix.values())
            values.emplace_back(value);
    return Expr::array({matrices.size(), rows, columns}, std::move(values));
}

[[nodiscard]] std::optional<std::size_t> exactPivotRow(
    const NumberMatrix& matrix,
    std::size_t firstRow,
    std::size_t column) {
    for (std::size_t row = firstRow; row < matrix.rows(); ++row)
        if (!matrix(row, column).isZero())
            return row;
    return std::nullopt;
}

[[nodiscard]] Expr exactLuOfNumbers(const MatrixView& source) {
    const std::size_t n = source.rows();
    NumberMatrix upper{source};
    NumberMatrix lower{n, n};
    NumberMatrix permutation{n, n};
    for (std::size_t i = 0; i < n; ++i) {
        lower(i, i) = Number{BigInt{1}};
        permutation(i, i) = Number{BigInt{1}};
    }

    for (std::size_t column = 0; column < n; ++column) {
        const auto pivot = exactPivotRow(upper, column, column);
        if (!pivot)
            continue;
        if (*pivot != column) {
            upper.swapRows(*pivot, column);
            permutation.swapRows(*pivot, column);
            for (std::size_t c = 0; c < column; ++c)
                std::swap(lower(*pivot, c), lower(column, c));
        }

        const Number pivotValue = upper(column, column);
        for (std::size_t row = column + 1; row < n; ++row) {
            if (upper(row, column).isZero())
                continue;
            const Number factor = upper(row, column) / pivotValue;
            lower(row, column) = factor;
            upper(row, column) = Number{BigInt{0}};
            for (std::size_t c = column + 1; c < n; ++c)
                upper(row, c) -= factor * upper(column, c);
        }
    }

    const NumberMatrix factors[] = {
        std::move(permutation), std::move(lower), std::move(upper)};
    return packedNumberMatrices(factors, n, n);
}

[[nodiscard]] bool exactZero(const Expr& value) {
    return value.isNumber() && value.asNumber().isZero();
}

[[nodiscard]] bool provablyNonZero(const Expr& value, const ExactMatrixContext& context) {
    if (value.isNumber())
        return !value.asNumber().isZero();
    const auto facts = mathematics::inferValueFacts(value, context.builtins, context.mathematics);
    return facts.sign == mathematics::RealSign::Positive
        || facts.sign == mathematics::RealSign::Negative
        || facts.sign == mathematics::RealSign::NonZero;
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

[[nodiscard]] Expr packedExprMatrices(
    std::span<const MatrixBuffer* const> matrices,
    std::size_t rows,
    std::size_t columns) {
    std::vector<Expr> output;
    const std::size_t shape[] = {matrices.size(), rows, columns};
    output.reserve(expression::arrayElementCount(shape));
    for (const MatrixBuffer* matrix : matrices)
        for (const Expr& value : matrix->elements())
            output.push_back(value);
    return Expr::array({matrices.size(), rows, columns}, std::move(output));
}

[[nodiscard]] Expr matrixExpr(const MatrixBuffer& matrix) {
    return Expr::array({matrix.rows(), matrix.columns()}, matrix.elements());
}

[[nodiscard]] Expr braceExprMatrices(std::span<const MatrixBuffer* const> matrices) {
    std::vector<Expr> factors;
    factors.reserve(matrices.size());
    for (const MatrixBuffer* matrix : matrices)
        factors.push_back(matrixExpr(*matrix));
    return expression::braceValue(std::move(factors));
}

[[nodiscard]] MatrixBuffer identityExprMatrix(std::size_t size) {
    MatrixBuffer result{size, size, integer(0)};
    for (std::size_t i = 0; i < size; ++i)
        result(i, i) = integer(1);
    return result;
}

[[nodiscard]] std::optional<Expr> symbolicLu(MatrixBuffer upper, const ExactMatrixContext& context) {
    const std::size_t n = upper.rows();
    MatrixBuffer lower{n, n, integer(0)};
    MatrixBuffer permutation{n, n, integer(0)};
    for (std::size_t i = 0; i < n; ++i) {
        lower(i, i) = integer(1);
        permutation(i, i) = integer(1);
    }

    for (std::size_t column = 0; column < n; ++column) {
        std::optional<std::size_t> pivot;
        bool undecidable = false;
        for (std::size_t row = column; row < n; ++row) {
            if (exactZero(upper(row, column)))
                continue;
            if (provablyNonZero(upper(row, column), context)) {
                pivot = row;
                break;
            }
            undecidable = true;
        }
        if (!pivot) {
            if (undecidable)
                return std::nullopt;
            continue;
        }
        if (*pivot != column) {
            upper.swapRows(*pivot, column);
            permutation.swapRows(*pivot, column);
            for (std::size_t c = 0; c < column; ++c)
                std::swap(lower(*pivot, c), lower(column, c));
        }

        const Expr pivotValue = upper(column, column);
        for (std::size_t row = column + 1; row < n; ++row) {
            if (exactZero(upper(row, column)))
                continue;
            const Expr factor = divide(upper(row, column), pivotValue, context);
            lower(row, column) = factor;
            upper(row, column) = integer(0);
            for (std::size_t c = column + 1; c < n; ++c)
                upper(row, c) = subtract(upper(row, c),
                    multiply(factor, upper(column, c), context), context);
        }
    }

    const MatrixBuffer* factors[] = {&permutation, &lower, &upper};
    return packedExprMatrices(factors, n, n);
}

[[nodiscard]] bool upperTriangular(const MatrixView& matrix) {
    for (std::size_t row = 1; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < std::min(row, matrix.columns()); ++column)
            if (!exactZero(matrix(row, column)))
                return false;
    return true;
}

using IntegerVector = std::vector<BigInt>;
using RationalVector = std::vector<Rational>;

[[nodiscard]] Rational exactRational(const MatrixView& matrix,
    std::size_t row,
    std::size_t column) {
    return matrix(row, column).asNumber().asReal().toRational();
}

[[nodiscard]] BigInt lcmPositive(const BigInt& lhs, const BigInt& rhs) {
    if (lhs.isZero() || rhs.isZero())
        return BigInt{};
    return (lhs / numeric::gcd(lhs, rhs)) * rhs;
}

void primitiveReduce(IntegerVector& vector) {
    BigInt content{};
    for (const BigInt& value : vector) {
        if (value.isZero())
            continue;
        content = content.isZero()
            ? value.abs()
            : numeric::gcd(std::move(content), value.abs());
        if (content == BigInt{1})
            return;
    }
    if (content.isZero() || content == BigInt{1})
        return;
    for (BigInt& value : vector)
        value /= content;
}

[[nodiscard]] bool zeroVector(const IntegerVector& vector) noexcept {
    return std::all_of(vector.begin(), vector.end(),
        [](const BigInt& value) { return value.isZero(); });
}

[[nodiscard]] BigInt integerDot(const IntegerVector& lhs, const IntegerVector& rhs) {
    BigInt result{};
    for (std::size_t i = 0; i < lhs.size(); ++i)
        result += lhs[i] * rhs[i];
    return result;
}

[[nodiscard]] Rational mixedDot(const IntegerVector& lhs, const RationalVector& rhs) {
    Rational result{};
    for (std::size_t i = 0; i < lhs.size(); ++i)
        result += Rational{lhs[i]} * rhs[i];
    return result;
}

[[nodiscard]] IntegerVector primitiveIntegerDirection(const RationalVector& source) {
    BigInt commonDenominator{1};
    for (const Rational& value : source)
        commonDenominator = lcmPositive(commonDenominator, value.denominator());

    IntegerVector result;
    result.reserve(source.size());
    for (const Rational& value : source)
        result.push_back(value.numerator() * (commonDenominator / value.denominator()));
    primitiveReduce(result);
    return result;
}

// p への射影を除く際に除算を作らず，整数vectorの方向だけを更新する。
void fractionFreeProjectAway(IntegerVector& vector,
    const IntegerVector& p,
    const BigInt& normSquared) {
    BigInt projection = integerDot(p, vector);
    if (projection.isZero())
        return;

    const BigInt common = numeric::gcd(normSquared, projection.abs());
    const BigInt normFactor = normSquared / common;
    projection /= common;
    for (std::size_t row = 0; row < vector.size(); ++row)
        vector[row] = normFactor * vector[row] - projection * p[row];
    primitiveReduce(vector);
}

[[nodiscard]] IntegerVector orthogonalizeFractionFree(
    IntegerVector vector,
    std::span<const IntegerVector> basis,
    std::span<const BigInt> normSquared) {
    for (std::size_t i = 0; i < basis.size(); ++i) {
        fractionFreeProjectAway(vector, basis[i], normSquared[i]);
        if (zeroVector(vector))
            break;
    }
    return vector;
}

[[nodiscard]] IntegerVector completionDirection(
    std::size_t rows,
    std::span<const IntegerVector> basis,
    std::span<const BigInt> normSquared) {
    for (std::size_t candidate = 0; candidate < rows; ++candidate) {
        IntegerVector vector(rows, BigInt{});
        vector[candidate] = BigInt{1};
        vector = orthogonalizeFractionFree(
            std::move(vector), basis, normSquared);
        if (!zeroVector(vector))
            return vector;
    }
    return {};
}

[[nodiscard]] std::optional<std::pair<std::vector<IntegerVector>, std::vector<BigInt>>>
bareissGramOrthogonalBasis(std::span<const IntegerVector> columns) {
    const std::size_t size = columns.size();
    if (size == 0)
        return std::pair{std::vector<IntegerVector>{}, std::vector<BigInt>{}};

    std::vector<IntegerVector> work(size, IntegerVector(size));
    std::vector<IntegerVector> lower(size, IntegerVector(size));
    for (std::size_t row = 0; row < size; ++row)
        for (std::size_t column = 0; column <= row; ++column)
            work[row][column] = work[column][row] = integerDot(columns[row], columns[column]);

    std::vector<BigInt> determinants(size + 1);
    determinants[0] = BigInt{1};
    BigInt previousPivot{1};
    for (std::size_t k = 0; k < size; ++k) {
        const BigInt pivot = work[k][k];
        if (pivot.isZero())
            return std::nullopt;
        determinants[k + 1] = pivot;
        lower[k][k] = pivot;
        for (std::size_t row = k + 1; row < size; ++row)
            lower[row][k] = work[row][k];

        for (std::size_t row = k + 1; row < size; ++row) {
            for (std::size_t column = row; column < size; ++column) {
                BigInt numerator = pivot * work[row][column]
                    - work[row][k] * work[k][column];
                const auto division = numeric::divmod(numerator, previousPivot);
                if (!division.remainder.isZero())
                    return std::nullopt;
                work[row][column] = work[column][row] = division.quotient;
            }
        }
        previousPivot = pivot;
    }

    std::vector<IntegerVector> basis;
    std::vector<BigInt> normSquared;
    basis.reserve(size);
    normSquared.reserve(size);
    for (std::size_t column = 0; column < size; ++column) {
        const BigInt diagonalScale = determinants[column] * determinants[column + 1];
        IntegerVector coefficients(size);
        for (std::size_t i = column + 1; i-- > 0;) {
            BigInt rhs = i == column ? diagonalScale : BigInt{};
            for (std::size_t j = i + 1; j <= column; ++j)
                rhs -= lower[j][i] * coefficients[j];
            const auto division = numeric::divmod(rhs, lower[i][i]);
            if (!division.remainder.isZero())
                return std::nullopt;
            coefficients[i] = division.quotient;
        }

        IntegerVector vector(columns.front().size());
        for (std::size_t i = 0; i <= column; ++i) {
            if (coefficients[i].isZero())
                continue;
            for (std::size_t row = 0; row < vector.size(); ++row)
                vector[row] += columns[i][row] * coefficients[i];
        }
        primitiveReduce(vector);
        if (zeroVector(vector))
            return std::nullopt;
        basis.push_back(std::move(vector));
        normSquared.push_back(integerDot(basis.back(), basis.back()));
    }
    return std::pair{std::move(basis), std::move(normSquared)};
}

[[nodiscard]] Expr rationalExpr(Rational value) {
    return Expr{Number{std::move(value)}};
}

[[nodiscard]] Expr scaledExactFactor(
    Rational coefficient,
    const Expr& factor,
    const ExactMatrixContext& context) {
    if (coefficient.isZero())
        return integer(0);
    if (coefficient == Rational{BigInt{1}})
        return factor;
    return Expr::call(context.builtins.symbol(evaluation::BuiltinId::Multiply),
        {rationalExpr(std::move(coefficient)), factor});
}

struct ExactInverseNorm final {
    Rational scale;
    std::optional<Expr> radical;
};

[[nodiscard]] ExactInverseNorm makeInverseIntegerNorm(
    const BigInt& normSquared,
    const ExactMatrixContext& context) {
    const auto decomposition = mathematics::decomposePositiveRationalSquareRoot(
        Rational{normSquared});
    if (decomposition.radicand == BigInt{1})
        return ExactInverseNorm{
            Rational{BigInt{1}} / decomposition.coefficient, std::nullopt};

    Rational scale = Rational{BigInt{1}}
        / (decomposition.coefficient * Rational{decomposition.radicand});
    Expr radical = Expr::call(context.builtins.symbol(evaluation::BuiltinId::Sqrt),
        {rationalExpr(Rational{decomposition.radicand})});
    return ExactInverseNorm{std::move(scale), std::move(radical)};
}

[[nodiscard]] Expr applyInverseNorm(
    Rational coefficient,
    const ExactInverseNorm& inverseNorm,
    const ExactMatrixContext& context) {
    coefficient *= inverseNorm.scale;
    if (!inverseNorm.radical)
        return rationalExpr(std::move(coefficient));
    return scaledExactFactor(std::move(coefficient), *inverseNorm.radical, context);
}

[[nodiscard]] std::optional<Expr> exactFractionFreeQrReal(
    const MatrixView& source,
    const ExactMatrixContext& context) {
    const std::size_t rows = source.rows();
    const std::size_t columns = source.columns();
    const std::size_t qColumns = std::min(rows, columns);

    std::vector<RationalVector> sourceColumns(columns, RationalVector(rows));
    for (std::size_t column = 0; column < columns; ++column)
        for (std::size_t row = 0; row < rows; ++row)
            sourceColumns[column][row] = exactRational(source, row, column);

    std::vector<IntegerVector> sourceDirections;
    sourceDirections.reserve(qColumns);
    for (std::size_t column = 0; column < qColumns; ++column)
        sourceDirections.push_back(primitiveIntegerDirection(sourceColumns[column]));

    std::vector<IntegerVector> basis;
    std::vector<BigInt> normSquared;
    if (auto gram = bareissGramOrthogonalBasis(sourceDirections)) {
        basis = std::move(gram->first);
        normSquared = std::move(gram->second);
    } else {
        // rank落ちや先頭principal minorが消えるcaseは，直接fraction-free直交化で基底を補完する。
        basis.reserve(qColumns);
        normSquared.reserve(qColumns);
        for (std::size_t column = 0; column < qColumns; ++column) {
            IntegerVector vector = orthogonalizeFractionFree(
                std::move(sourceDirections[column]), basis, normSquared);
            if (zeroVector(vector))
                vector = completionDirection(rows, basis, normSquared);
            if (zeroVector(vector))
                return std::nullopt;

            BigInt squared = integerDot(vector, vector);
            if (squared.isZero())
                return std::nullopt;
            basis.push_back(std::move(vector));
            normSquared.push_back(std::move(squared));
        }
    }

    MatrixBuffer q{rows, qColumns, integer(0)};
    MatrixBuffer r{qColumns, columns, integer(0)};
    for (std::size_t i = 0; i < qColumns; ++i) {
        const ExactInverseNorm inverseNorm = makeInverseIntegerNorm(normSquared[i], context);
        for (std::size_t row = 0; row < rows; ++row) {
            if (basis[i][row].isZero())
                continue;
            q(row, i) = applyInverseNorm(
                Rational{basis[i][row]}, inverseNorm, context);
        }
        for (std::size_t column = 0; column < columns; ++column) {
            const Rational projection = mixedDot(basis[i], sourceColumns[column]);
            if (projection.isZero())
                continue;
            r(i, column) = applyInverseNorm(projection, inverseNorm, context);
        }
    }

    const MatrixBuffer* factors[] = {&q, &r};
    return braceExprMatrices(factors);
}

[[nodiscard]] ComplexInterval exactComplex(std::int64_t real, std::size_t bits) {
    return ComplexInterval::fromReal(RealInterval::fromRational(
        Rational{BigInt{real}}, bits));
}
[[nodiscard]] bool exactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}
[[nodiscard]] bool exactZero(const ComplexInterval& value) noexcept {
    return exactZero(value.real()) && exactZero(value.imaginary());
}
[[nodiscard]] ComplexInterval conjugate(const ComplexInterval& value) {
    return ComplexInterval{value.real(), approximation::negate(value.imaginary())};
}

class IntervalMatrix final {
public:
    IntervalMatrix(std::size_t rows, std::size_t columns, std::vector<ComplexInterval> values)
        : rows_(rows), columns_(columns), values_(std::move(values)) {
        const std::size_t shape[] = {rows, columns};
        if (values_.size() != expression::arrayElementCount(shape))
            throw std::invalid_argument("IntervalMatrix shape does not match element count");
    }
    IntervalMatrix(std::size_t rows, std::size_t columns, ComplexInterval initial)
        : rows_(rows), columns_(columns) {
        const std::size_t shape[] = {rows, columns};
        values_.assign(expression::arrayElementCount(shape), std::move(initial));
    }
    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] ComplexInterval& operator()(std::size_t row, std::size_t column) noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const ComplexInterval& operator()(std::size_t row, std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }
    void swapRows(std::size_t lhs, std::size_t rhs) noexcept {
        if (lhs == rhs)
            return;
        for (std::size_t column = 0; column < columns_; ++column)
            std::swap((*this)(lhs, column), (*this)(rhs, column));
    }
    [[nodiscard]] const std::vector<ComplexInterval>& values() const noexcept { return values_; }
private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<ComplexInterval> values_;
};

[[nodiscard]] std::optional<std::vector<ComplexInterval>> encloseElements(
    const ArrayExpr& array,
    std::size_t bits,
    const approximation::CertifiedEvaluator& certified) {
    std::vector<ComplexInterval> values;
    values.reserve(array.size());
    for (std::size_t i = 0; i < array.size(); ++i) {
        const Expr element = array.element(i);
        const auto value = approximation::encloseComplexExpression(element, bits, certified);
        if (!value)
            return std::nullopt;
        values.push_back(*value);
    }
    return values;
}

[[nodiscard]] RealInterval magnitudeSquared(const ComplexInterval& value, std::size_t bits) {
    return approximation::add(
        approximation::squareInterval(value.real(), bits),
        approximation::squareInterval(value.imaginary(), bits), bits);
}

[[nodiscard]] std::optional<std::size_t> intervalPivotRow(
    const IntervalMatrix& matrix,
    std::size_t firstRow,
    std::size_t column,
    std::size_t bits) {
    bool uncertain = false;
    std::optional<std::size_t> bestRow;
    std::optional<numeric::BigFloat> bestLowerMagnitudeSquared;
    for (std::size_t row = firstRow; row < matrix.rows(); ++row) {
        const auto& value = matrix(row, column);
        if (value.containsZero()) {
            if (!exactZero(value))
                uncertain = true;
            continue;
        }

        const auto magnitude = magnitudeSquared(value, bits);
        const auto lower = magnitude.lower();
        if (!bestRow || lower > *bestLowerMagnitudeSquared) {
            bestRow = row;
            bestLowerMagnitudeSquared = lower;
        }
    }
    if (bestRow)
        return bestRow;
    if (uncertain)
        throw approximation::PrecisionInsufficient(
            "LU pivot nonzero status is not certified at the current precision");
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> decimalPacked(
    std::vector<std::size_t> shape,
    const std::vector<ComplexInterval>& values,
    std::size_t digits) {
    std::vector<Expr> output;
    output.reserve(values.size());
    for (const ComplexInterval& value : values) {
        const auto converted = approximation::decimalExpression(value, digits);
        if (!converted)
            return std::nullopt;
        output.push_back(*converted);
    }
    return Expr::array(std::move(shape), std::move(output));
}

[[nodiscard]] std::optional<Expr> decimalMatrix(
    std::size_t rows,
    std::size_t columns,
    const std::vector<ComplexInterval>& values,
    std::size_t digits) {
    if (values.size() != rows * columns)
        return std::nullopt;
    std::vector<Expr> output;
    output.reserve(values.size());
    for (const ComplexInterval& value : values) {
        const auto converted = approximation::decimalExpression(value, digits);
        if (!converted)
            return std::nullopt;
        output.push_back(*converted);
    }
    return Expr::array({rows, columns}, std::move(output));
}

template <class Operation>
[[nodiscard]] std::optional<Expr> retryApproximation(
    approximation::ApproximationContext context,
    Operation&& operation) {
    for (std::size_t attempt = 0; attempt < maximumPrecisionRetries; ++attempt) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        try {
            if (const auto result = operation(context))
                return result;
        }
        catch (const approximation::PrecisionInsufficient&) {
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            return std::nullopt;
        }
        context.setGuardDigits(approximation::nextGuardDigits(context.guardDigits()));
    }
    return std::nullopt;
}

struct IntervalReflector final {
    std::size_t firstRow = 0;
    std::vector<ComplexInterval> vector;
    ComplexInterval beta;
};

void applyIntervalHouseholderLeft(
    IntervalMatrix& matrix,
    const IntervalReflector& reflector,
    std::size_t firstColumn,
    std::size_t bits,
    std::size_t blockColumns) {
    const std::size_t width = std::max<std::size_t>(1, blockColumns);
    for (std::size_t block = firstColumn; block < matrix.columns(); block += width) {
        const std::size_t end = std::min(matrix.columns(), block + width);
        std::vector<ComplexInterval> dots(end - block, exactComplex(0, bits));
        for (std::size_t i = 0; i < reflector.vector.size(); ++i) {
            const ComplexInterval cv = conjugate(reflector.vector[i]);
            const std::size_t row = reflector.firstRow + i;
            for (std::size_t column = block; column < end; ++column)
                dots[column - block] = approximation::add(dots[column - block],
                    approximation::multiply(cv, matrix(row, column), bits), bits);
        }
        for (std::size_t i = 0; i < reflector.vector.size(); ++i) {
            const std::size_t row = reflector.firstRow + i;
            for (std::size_t column = block; column < end; ++column) {
                const ComplexInterval correction = approximation::multiply(
                    reflector.vector[i],
                    approximation::multiply(reflector.beta, dots[column - block], bits), bits);
                matrix(row, column) = approximation::subtract(
                    matrix(row, column), correction, bits);
            }
        }
    }
}

[[nodiscard]] std::optional<Expr> approximateLuAtPrecision(
    const ArrayExpr& source,
    std::size_t bits,
    std::size_t digits,
    const approximation::CertifiedEvaluator& certified) {
    const auto enclosed = encloseElements(source, bits, certified);
    if (!enclosed)
        return std::nullopt;
    const std::size_t n = source.shape[0];
    IntervalMatrix upper{n, n, *enclosed};
    IntervalMatrix lower{n, n, exactComplex(0, bits)};
    IntervalMatrix permutation{n, n, exactComplex(0, bits)};
    for (std::size_t i = 0; i < n; ++i) {
        lower(i, i) = exactComplex(1, bits);
        permutation(i, i) = exactComplex(1, bits);
    }

    for (std::size_t column = 0; column < n; ++column) {
        const auto pivot = intervalPivotRow(upper, column, column, bits);
        if (!pivot)
            continue;
        if (*pivot != column) {
            upper.swapRows(*pivot, column);
            permutation.swapRows(*pivot, column);
            for (std::size_t c = 0; c < column; ++c)
                std::swap(lower(*pivot, c), lower(column, c));
        }
        const ComplexInterval pivotValue = upper(column, column);
        for (std::size_t row = column + 1; row < n; ++row) {
            if (exactZero(upper(row, column)))
                continue;
            const ComplexInterval factor = approximation::divide(
                upper(row, column), pivotValue, bits);
            lower(row, column) = factor;
            upper(row, column) = exactComplex(0, bits);
            for (std::size_t c = column + 1; c < n; ++c)
                upper(row, c) = approximation::subtract(upper(row, c),
                    approximation::multiply(factor, upper(column, c), bits), bits);
        }
    }

    std::vector<ComplexInterval> output;
    const std::size_t shape[] = {3, n, n};
    output.reserve(expression::arrayElementCount(shape));
    for (const IntervalMatrix* matrix : {&permutation, &lower, &upper})
        output.insert(output.end(), matrix->values().begin(), matrix->values().end());
    return decimalPacked({3, n, n}, output, digits);
}

[[nodiscard]] std::optional<Expr> approximateQrAtPrecision(
    const ArrayExpr& source,
    std::size_t bits,
    std::size_t digits,
    const approximation::CertifiedEvaluator& certified,
    std::size_t blockColumns) {
    const auto enclosed = encloseElements(source, bits, certified);
    if (!enclosed)
        return std::nullopt;
    const std::size_t rows = source.shape[0];
    const std::size_t columns = source.shape[1];
    const std::size_t qColumns = std::min(rows, columns);
    const std::size_t reflectorSteps = rows == 0 ? 0 : std::min(columns, rows - 1);
    IntervalMatrix r{rows, columns, *enclosed};
    std::vector<IntervalReflector> reflectors;
    reflectors.reserve(reflectorSteps);

    for (std::size_t k = 0; k < reflectorSteps; ++k) {
        RealInterval normSq = RealInterval::fromRational(Rational{}, bits);
        for (std::size_t row = k; row < rows; ++row)
            normSq = approximation::add(normSq, magnitudeSquared(r(row, k), bits), bits);
        if (exactZero(normSq))
            continue;
        if (normSq.containsZero())
            throw approximation::PrecisionInsufficient(
                "QR column norm is not certified nonzero at the current precision");
        const RealInterval norm = approximation::encloseSqrt(normSq, bits).interval;

        ComplexInterval phase = exactComplex(1, bits);
        const ComplexInterval x0 = r(k, k);
        if (!exactZero(x0)) {
            const RealInterval x0Sq = magnitudeSquared(x0, bits);
            if (x0Sq.containsZero())
                throw approximation::PrecisionInsufficient(
                    "QR Householder phase is not certified at the current precision");
            const RealInterval magnitude = approximation::encloseSqrt(x0Sq, bits).interval;
            phase = approximation::divide(x0, ComplexInterval::fromReal(magnitude), bits);
        }
        const ComplexInterval alpha = approximation::negate(
            approximation::multiply(phase, ComplexInterval::fromReal(norm), bits));

        std::vector<ComplexInterval> v;
        v.reserve(rows - k);
        v.push_back(approximation::subtract(x0, alpha, bits));
        for (std::size_t row = k + 1; row < rows; ++row)
            v.push_back(r(row, k));

        RealInterval vNormSq = RealInterval::fromRational(Rational{}, bits);
        for (const ComplexInterval& item : v)
            vNormSq = approximation::add(vNormSq, magnitudeSquared(item, bits), bits);
        if (vNormSq.containsZero())
            throw approximation::PrecisionInsufficient(
                "QR Householder vector norm is not certified at the current precision");
        const RealInterval betaReal = approximation::divide(
            RealInterval::fromRational(Rational{BigInt{2}}, bits), vNormSq, bits);
        IntervalReflector reflector{k, std::move(v), ComplexInterval::fromReal(betaReal)};
        applyIntervalHouseholderLeft(r, reflector, k, bits, blockColumns);
        for (std::size_t row = k + 1; row < rows; ++row)
            r(row, k) = exactComplex(0, bits);
        reflectors.push_back(std::move(reflector));
    }

    IntervalMatrix q{rows, qColumns, exactComplex(0, bits)};
    for (std::size_t i = 0; i < qColumns; ++i)
        q(i, i) = exactComplex(1, bits);
    for (std::size_t i = reflectors.size(); i-- > 0;)
        applyIntervalHouseholderLeft(q, reflectors[i], 0, bits, blockColumns);

    std::vector<ComplexInterval> reducedR;
    reducedR.reserve(qColumns * columns);
    for (std::size_t row = 0; row < qColumns; ++row)
        for (std::size_t column = 0; column < columns; ++column)
            reducedR.push_back(r(row, column));

    const auto qExpr = decimalMatrix(rows, qColumns, q.values(), digits);
    const auto rExpr = decimalMatrix(qColumns, columns, reducedR, digits);
    if (!qExpr || !rExpr)
        return std::nullopt;
    return expression::braceValue({*qExpr, *rExpr});
}

} // namespace

std::optional<Expr> luDecomposition(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    if (matrix.rows() != matrix.columns())
        throw std::invalid_argument("luDecomposition currently requires a square matrix");
    if (upperTriangular(matrix)) {
        MatrixBuffer permutation = identityExprMatrix(matrix.rows());
        MatrixBuffer lower = identityExprMatrix(matrix.rows());
        MatrixBuffer upper{matrix};
        const MatrixBuffer* factors[] = {&permutation, &lower, &upper};
        return packedExprMatrices(factors, matrix.rows(), matrix.columns());
    }
    if (allExactNumbers(matrix))
        return exactLuOfNumbers(matrix);
    return symbolicLu(MatrixBuffer{matrix}, context);
}

std::optional<Expr> qrDecomposition(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    if (upperTriangular(matrix)) {
        const std::size_t k = std::min(matrix.rows(), matrix.columns());
        MatrixBuffer q{matrix.rows(), k, integer(0)};
        for (std::size_t i = 0; i < k; ++i)
            q(i, i) = integer(1);
        MatrixBuffer r{k, matrix.columns(), integer(0)};
        for (std::size_t row = 0; row < k; ++row)
            for (std::size_t column = 0; column < matrix.columns(); ++column)
                r(row, column) = matrix(row, column);
        const MatrixBuffer* factors[] = {&q, &r};
        return braceExprMatrices(factors);
    }
    if (!allExactRealNumbers(matrix))
        return std::nullopt;
    return exactFractionFreeQrReal(matrix, context);
}

std::optional<Expr> approximateLuDecomposition(
    const ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (!matrix.isMatrix() || matrix.shape[0] != matrix.shape[1])
        return std::nullopt;
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) {
        return approximateLuAtPrecision(matrix, current.workingBinaryBits(),
            current.decimalDigits(), certified);
    });
}

std::optional<Expr> approximateQrDecompositionWithBlockSize(
    const ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context,
    std::size_t blockColumns) {
    if (!matrix.isMatrix())
        return std::nullopt;
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) {
        return approximateQrAtPrecision(matrix, current.workingBinaryBits(),
            current.decimalDigits(), certified, blockColumns);
    });
}

std::optional<Expr> approximateQrDecomposition(
    const ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return approximateQrDecompositionWithBlockSize(
        matrix, builtins, mathematics, angles, context, defaultQrBlockColumns);
}

} // namespace mmcal::linear_algebra
