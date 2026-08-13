// precision-aware certified Matrix backend
#include "approximate_matrix.hpp"

#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/certified_sqrt.hpp"
#include "approximation/expression_interval.hpp"
#include "approximation/interval_math.hpp"
#include "builtins/array_helpers.hpp"
#include "error/error_message.hpp"
#include "expression/array_utils.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra {
namespace {

using approximation::ComplexInterval;
using approximation::RealInterval;
using expression::ArrayExpr;
using expression::Expr;
using numeric::BigInt;
using numeric::Rational;

constexpr std::size_t maximumPrecisionRetries = 12;

[[nodiscard]] ComplexInterval exactComplex(
    std::int64_t real,
    std::size_t precisionBits) {
    return ComplexInterval::fromReal(RealInterval::fromRational(
        Rational{BigInt{real}}, precisionBits));
}

[[nodiscard]] bool exactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

[[nodiscard]] bool exactZero(const ComplexInterval& value) noexcept {
    return exactZero(value.real()) && exactZero(value.imaginary());
}

[[nodiscard]] std::size_t augmentedColumns(std::size_t columns) {
    if (columns > std::numeric_limits<std::size_t>::max() / 2)
        throw std::length_error("Augmented matrix column count exceeds the size_t range");
    return columns * 2;
}

class IntervalMatrix final {
public:
    IntervalMatrix(std::size_t rows, std::size_t columns, std::vector<ComplexInterval> values)
        : rows_(rows), columns_(columns), values_(std::move(values)) {
        const std::size_t shape[] = {rows_, columns_};
        if (values_.size() != expression::arrayElementCount(shape))
            throw std::invalid_argument("IntervalMatrix shape does not match the element count");
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

private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<ComplexInterval> values_;
};

[[nodiscard]] std::optional<std::vector<ComplexInterval>> encloseElements(
    const ArrayExpr& array,
    std::size_t precisionBits,
    const approximation::CertifiedEvaluator& certified) {
    std::vector<ComplexInterval> values;
    values.reserve(array.elements.size());
    for (const Expr& element : array.elements) {
        const auto enclosed = approximation::encloseComplexExpression(
            element, precisionBits, certified);
        if (!enclosed)
            return std::nullopt;
        values.push_back(*enclosed);
    }
    return values;
}

[[nodiscard]] std::optional<IntervalMatrix> encloseMatrix(
    const ArrayExpr& array,
    std::size_t precisionBits,
    const approximation::CertifiedEvaluator& certified) {
    const auto values = encloseElements(array, precisionBits, certified);
    if (!values)
        return std::nullopt;
    return IntervalMatrix{array.shape[0], array.shape[1], *values};
}

[[nodiscard]] std::optional<std::size_t> selectPivot(
    const IntervalMatrix& matrix,
    std::size_t firstRow,
    std::size_t column) {
    bool uncertain = false;
    for (std::size_t row = firstRow; row < matrix.rows(); ++row) {
        const ComplexInterval& value = matrix(row, column);
        if (!value.containsZero())
            return row;
        if (!exactZero(value))
            uncertain = true;
    }
    if (uncertain)
        throw approximation::PrecisionInsufficient(
            "Matrix pivot nonzero status is not certified at the current precision");
    return std::nullopt;
}

struct RrefResult final {
    IntervalMatrix matrix;
    std::size_t rank = 0;
    std::vector<std::size_t> pivotColumns;
};

[[nodiscard]] RrefResult intervalRref(
    IntervalMatrix matrix,
    std::size_t precisionBits,
    std::size_t pivotColumnLimit) {
    std::size_t pivotRow = 0;
    const std::size_t limit = std::min(pivotColumnLimit, matrix.columns());
    std::vector<std::size_t> pivotColumns;
    pivotColumns.reserve(std::min(matrix.rows(), limit));

    for (std::size_t column = 0; column < limit && pivotRow < matrix.rows(); ++column) {
        const auto selected = selectPivot(matrix, pivotRow, column);
        if (!selected)
            continue;

        matrix.swapRows(*selected, pivotRow);
        const ComplexInterval pivot = matrix(pivotRow, column);
        for (std::size_t c = 0; c < matrix.columns(); ++c)
            matrix(pivotRow, c) = approximation::divide(
                matrix(pivotRow, c), pivot, precisionBits);

        for (std::size_t row = 0; row < matrix.rows(); ++row) {
            if (row == pivotRow || exactZero(matrix(row, column)))
                continue;
            const ComplexInterval factor = matrix(row, column);
            for (std::size_t c = 0; c < matrix.columns(); ++c)
                matrix(row, c) = approximation::subtract(
                    matrix(row, c),
                    approximation::multiply(factor, matrix(pivotRow, c), precisionBits),
                    precisionBits);
        }
        pivotColumns.push_back(column);
        ++pivotRow;
    }
    return RrefResult{std::move(matrix), pivotRow, std::move(pivotColumns)};
}

[[nodiscard]] std::optional<Expr> decimalArray(
    std::vector<std::size_t> shape,
    const std::vector<ComplexInterval>& values,
    std::size_t digits) {
    std::vector<Expr> output;
    output.reserve(values.size());
    for (const ComplexInterval& value : values) {
        const auto decimal = approximation::decimalExpression(value, digits);
        if (!decimal)
            return std::nullopt;
        output.push_back(*decimal);
    }
    return Expr::array(std::move(shape), std::move(output));
}

template <class Operation>
[[nodiscard]] std::optional<Expr> retryApproximation(
    approximation::ApproximationContext context,
    Operation&& operation) {
    for (std::size_t attempt = 0; attempt < maximumPrecisionRetries; ++attempt) {
        try {
            if (const auto result = operation(context))
                return result;
        }
        catch (const approximation::PrecisionInsufficient&) {
        }
        context.setGuardDigits(approximation::nextGuardDigits(context.guardDigits()));
    }
    return std::nullopt;
}

[[nodiscard]] RealInterval normSquared(
    const std::vector<ComplexInterval>& values,
    std::size_t bits) {
    RealInterval result = RealInterval::fromRational(Rational{}, bits);
    for (const ComplexInterval& value : values) {
        const RealInterval magnitudeSquared = approximation::add(
            approximation::squareInterval(value.real(), bits),
            approximation::squareInterval(value.imaginary(), bits),
            bits);
        result = approximation::add(result, magnitudeSquared, bits);
    }
    return result;
}

} // namespace

std::optional<Expr> approximateDot(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 2 || !arguments[0].isArray() || !arguments[1].isArray())
        return std::nullopt;
    const ArrayExpr& lhs = arguments[0].asArray();
    const ArrayExpr& rhs = arguments[1].asArray();
    if (lhs.rank() < 1 || lhs.rank() > 2 || rhs.rank() < 1 || rhs.rank() > 2)
        return std::nullopt;

    const std::size_t lhsRows = lhs.rank() == 1 ? 1 : lhs.shape[0];
    const std::size_t lhsColumns = lhs.rank() == 1 ? lhs.shape[0] : lhs.shape[1];
    const std::size_t rhsRows = rhs.shape[0];
    const std::size_t rhsColumns = rhs.rank() == 1 ? 1 : rhs.shape[1];
    if (lhsColumns != rhsRows)
        return std::nullopt;

    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto left = encloseElements(lhs, bits, certified);
        const auto right = encloseElements(rhs, bits, certified);
        if (!left || !right)
            return std::nullopt;

        auto lhsAt = [&](std::size_t row, std::size_t column) -> const ComplexInterval& {
            return lhs.rank() == 1 ? (*left)[column] : (*left)[row * lhsColumns + column];
        };
        auto rhsAt = [&](std::size_t row, std::size_t column) -> const ComplexInterval& {
            return rhs.rank() == 1 ? (*right)[row] : (*right)[row * rhsColumns + column];
        };

        std::vector<ComplexInterval> output;
        const std::size_t outputShape[] = {lhsRows, rhsColumns};
        output.reserve(expression::arrayElementCount(outputShape));
        for (std::size_t row = 0; row < lhsRows; ++row) {
            for (std::size_t column = 0; column < rhsColumns; ++column) {
                ComplexInterval sum = exactComplex(0, bits);
                for (std::size_t k = 0; k < lhsColumns; ++k)
                    sum = approximation::add(sum,
                        approximation::multiply(lhsAt(row, k), rhsAt(k, column), bits), bits);
                output.push_back(std::move(sum));
            }
        }

        if (lhs.rank() == 1 && rhs.rank() == 1)
            return approximation::decimalExpression(output.front(), current.decimalDigits());
        const std::vector<std::size_t> shape = lhs.rank() == 1 || rhs.rank() == 1
            ? std::vector<std::size_t>{lhs.rank() == 1 ? rhsColumns : lhsRows}
            : std::vector<std::size_t>{lhsRows, rhsColumns};
        return decimalArray(shape, output, current.decimalDigits());
    });
}

std::optional<Expr> approximateDeterminant(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& array = arguments[0].asArray();
    if (array.shape[0] != array.shape[1])
        return std::nullopt;

    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        auto matrix = encloseMatrix(array, bits, certified);
        if (!matrix)
            return std::nullopt;

        ComplexInterval determinant = exactComplex(1, bits);
        bool negative = false;
        for (std::size_t column = 0; column < matrix->columns(); ++column) {
            const auto pivot = selectPivot(*matrix, column, column);
            if (!pivot)
                return approximation::decimalExpression(exactComplex(0, bits), current.decimalDigits());
            if (*pivot != column) {
                matrix->swapRows(*pivot, column);
                negative = !negative;
            }

            const ComplexInterval pivotValue = (*matrix)(column, column);
            determinant = approximation::multiply(determinant, pivotValue, bits);
            for (std::size_t row = column + 1; row < matrix->rows(); ++row) {
                if (exactZero((*matrix)(row, column)))
                    continue;
                const ComplexInterval factor = approximation::divide(
                    (*matrix)(row, column), pivotValue, bits);
                (*matrix)(row, column) = exactComplex(0, bits);
                for (std::size_t c = column + 1; c < matrix->columns(); ++c)
                    (*matrix)(row, c) = approximation::subtract(
                        (*matrix)(row, c),
                        approximation::multiply(factor, (*matrix)(column, c), bits), bits);
            }
        }
        if (negative)
            determinant = approximation::negate(determinant);
        return approximation::decimalExpression(determinant, current.decimalDigits());
    });
}

std::optional<Expr> approximateInverse(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    if (source.shape[0] != source.shape[1])
        return std::nullopt;

    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto input = encloseElements(source, bits, certified);
        if (!input)
            return std::nullopt;
        const std::size_t n = source.shape[0];
        const std::size_t columns = augmentedColumns(n);
        const std::size_t augmentedShape[] = {n, columns};

        std::vector<ComplexInterval> values;
        values.reserve(expression::arrayElementCount(augmentedShape));
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t column = 0; column < n; ++column)
                values.push_back((*input)[row * n + column]);
            for (std::size_t column = 0; column < n; ++column)
                values.push_back(exactComplex(row == column ? 1 : 0, bits));
        }
        IntervalMatrix augmented{n, columns, std::move(values)};
        const RrefResult reduced = intervalRref(std::move(augmented), bits, n);
        if (reduced.rank != n)
            throw std::domain_error("Matrix is singular");

        std::vector<ComplexInterval> output;
        output.reserve(source.elements.size());
        for (std::size_t row = 0; row < n; ++row)
            for (std::size_t column = 0; column < n; ++column)
                output.push_back(reduced.matrix(row, n + column));
        return decimalArray({n, n}, output, current.decimalDigits());
    });
}

std::optional<Expr> approximateRref(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        auto matrix = encloseMatrix(source, bits, certified);
        if (!matrix)
            return std::nullopt;
        const RrefResult reduced = intervalRref(std::move(*matrix), bits, source.shape[1]);

        std::vector<ComplexInterval> output;
        output.reserve(source.elements.size());
        for (std::size_t row = 0; row < reduced.matrix.rows(); ++row)
            for (std::size_t column = 0; column < reduced.matrix.columns(); ++column)
                output.push_back(reduced.matrix(row, column));
        return decimalArray(source.shape, output, current.decimalDigits());
    });
}

std::optional<Expr> approximateMatrixRank(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        auto matrix = encloseMatrix(source, bits, certified);
        if (!matrix)
            return std::nullopt;
        const RrefResult reduced = intervalRref(std::move(*matrix), bits, source.shape[1]);
        return Expr{numeric::Number{BigInt::parse(std::to_string(reduced.rank))}};
    });
}

std::optional<Expr> approximateSolveLinear(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 2 || !arguments[0].isArray() || !arguments[1].isArray()
        || !arguments[0].asArray().isMatrix() || !arguments[1].asArray().isVector())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    const ArrayExpr& rhs = arguments[1].asArray();
    if (rhs.shape[0] != source.shape[0])
        return std::nullopt;
    if (source.shape[1] == std::numeric_limits<std::size_t>::max())
        throw std::length_error("Linear system augmented column count exceeds the size_t range");

    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto coefficients = encloseElements(source, bits, certified);
        const auto right = encloseElements(rhs, bits, certified);
        if (!coefficients || !right)
            return std::nullopt;

        const std::size_t rows = source.shape[0];
        const std::size_t variables = source.shape[1];
        const std::size_t augmentedColumns = variables + 1;
        const std::size_t augmentedShape[] = {rows, augmentedColumns};
        std::vector<ComplexInterval> values;
        values.reserve(expression::arrayElementCount(augmentedShape));
        for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t column = 0; column < variables; ++column)
                values.push_back((*coefficients)[row * variables + column]);
            values.push_back((*right)[row]);
        }

        RrefResult reduced = intervalRref(
            IntervalMatrix{rows, augmentedColumns, std::move(values)}, bits, variables);
        for (std::size_t row = reduced.rank; row < rows; ++row) {
            const ComplexInterval& residual = reduced.matrix(row, variables);
            if (!residual.containsZero())
                throw std::domain_error("Linear system is inconsistent");
            if (!exactZero(residual))
                throw approximation::PrecisionInsufficient(
                    "Linear system consistency is not certified at the current precision");
        }
        if (reduced.rank != variables)
            throw std::domain_error("Linear system does not have a unique solution");

        std::vector<ComplexInterval> solution;
        solution.reserve(variables);
        for (std::size_t variable = 0; variable < variables; ++variable)
            solution.push_back(reduced.matrix(variable, variables));
        return decimalArray({variables}, solution, current.decimalDigits());
    });
}

std::optional<Expr> approximateNullSpace(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        auto matrix = encloseMatrix(source, bits, certified);
        if (!matrix)
            return std::nullopt;
        RrefResult reduced = intervalRref(std::move(*matrix), bits, source.shape[1]);

        const std::size_t variables = source.shape[1];
        const std::size_t noPivot = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> pivotRowByColumn(variables, noPivot);
        for (std::size_t row = 0; row < reduced.pivotColumns.size(); ++row)
            pivotRowByColumn[reduced.pivotColumns[row]] = row;

        const std::size_t nullity = variables - reduced.pivotColumns.size();
        std::vector<ComplexInterval> basis;
        const std::size_t shape[] = {nullity, variables};
        basis.reserve(expression::arrayElementCount(shape));
        for (std::size_t freeColumn = 0; freeColumn < variables; ++freeColumn) {
            if (pivotRowByColumn[freeColumn] != noPivot)
                continue;
            for (std::size_t variable = 0; variable < variables; ++variable) {
                if (variable == freeColumn) {
                    basis.push_back(exactComplex(1, bits));
                    continue;
                }

                const std::size_t pivotRow = pivotRowByColumn[variable];
                if (pivotRow == noPivot) {
                    basis.push_back(exactComplex(0, bits));
                    continue;
                }
                basis.push_back(approximation::negate(
                    reduced.matrix(pivotRow, freeColumn)));
            }
        }
        return decimalArray({nullity, variables}, basis, current.decimalDigits());
    });
}

std::optional<Expr> approximateNorm(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isVector())
        return std::nullopt;
    const ArrayExpr& vector = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto values = encloseElements(vector, bits, certified);
        if (!values)
            return std::nullopt;
        const RealInterval squared = normSquared(*values, bits);
        const RealInterval norm = approximation::encloseSqrt(squared, bits).interval;
        return approximation::decimalExpression(norm, current.decimalDigits());
    });
}

std::optional<Expr> approximateNormalize(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isVector())
        return std::nullopt;
    const ArrayExpr& vector = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto values = encloseElements(vector, bits, certified);
        if (!values)
            return std::nullopt;
        const RealInterval squared = normSquared(*values, bits);
        const RealInterval norm = approximation::encloseSqrt(squared, bits).interval;
        if (exactZero(norm))
            throw std::domain_error("normalize requires a nonzero vector");
        if (norm.containsZero())
            throw approximation::PrecisionInsufficient(
                "Vector norm nonzero status is not certified at the current precision");

        const ComplexInterval denominator = ComplexInterval::fromReal(norm);
        std::vector<ComplexInterval> output;
        output.reserve(values->size());
        for (const ComplexInterval& value : *values)
            output.push_back(approximation::divide(value, denominator, bits));
        return decimalArray({vector.shape[0]}, output, current.decimalDigits());
    });
}

std::optional<Expr> approximateTrace(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& matrix = arguments[0].asArray();
    if (matrix.shape[0] != matrix.shape[1])
        return std::nullopt;
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto values = encloseElements(matrix, bits, certified);
        if (!values)
            return std::nullopt;
        ComplexInterval sum = exactComplex(0, bits);
        for (std::size_t i = 0; i < matrix.shape[0]; ++i)
            sum = approximation::add(sum, (*values)[i * matrix.shape[1] + i], bits);
        return approximation::decimalExpression(sum, current.decimalDigits());
    });
}

} // namespace mmcal::linear_algebra
