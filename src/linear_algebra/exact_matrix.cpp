// exact Matrix elimination。Number行列はExpr/Simplifierを経由せず直接処理する。
#include "exact_matrix.hpp"
#include "linear_algebra/exact_matrix_detail.hpp"
#include "linear_algebra/number_matrix.hpp"
#include "evaluation/evaluation_budget.hpp"

#include "builtins/exact_operations.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "expression/array_utils.hpp"
#include "linear_algebra/fraction_free_elimination.hpp"
#include "linear_algebra/modular_linear_algebra.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra {
namespace {

using expression::Expr;
using numeric::BigInt;
using numeric::Number;

using detail::NumberMatrix;
using detail::add;
using detail::divide;
using detail::exactZero;
using detail::integer;
using detail::multiply;
using detail::negate;
using detail::provablyNonZero;
using detail::simplify;
using detail::subtract;


[[nodiscard]] std::size_t augmentedColumns(std::size_t columns);

struct IntegerLift final {
    IntegerMatrixBuffer matrix;
    std::vector<BigInt> rowScales;
};

[[nodiscard]] bool allExactRealNumbers(const MatrixView& matrix) noexcept {
    return matrix.array().hasExactRealStorage();
}

template <class RationalAt>
[[nodiscard]] IntegerLift liftRealRows(
    std::size_t rows,
    std::size_t columns,
    RationalAt&& rationalAt) {
    IntegerMatrixBuffer lifted{rows, columns};
    std::vector<BigInt> rowScales(rows, BigInt{1});

    for (std::size_t row = 0; row < rows; ++row) {
        BigInt scale{1};
        for (std::size_t column = 0; column < columns; ++column) {
            const auto rational = rationalAt(row, column);
            scale = numeric::lcm(scale, rational.denominator());
        }
        rowScales[row] = scale;

        for (std::size_t column = 0; column < columns; ++column) {
            const auto rational = rationalAt(row, column);
            lifted(row, column) = rational.numerator()
                * (scale / rational.denominator());
        }
    }
    return {std::move(lifted), std::move(rowScales)};
}

[[nodiscard]] IntegerLift liftRealMatrix(const MatrixView& source) {
    return liftRealRows(source.rows(), source.columns(),
        [&](std::size_t row, std::size_t column) {
            return source.array().exactNumber(row * source.columns() + column)
                .asReal().toRational();
        });
}

[[nodiscard]] BigInt rowScaleProduct(const std::vector<BigInt>& scales) {
    BigInt result{1};
    for (const BigInt& scale : scales)
        result *= scale;
    return result;
}

[[nodiscard]] std::vector<numeric::Rational> bareissBackSubstituteUnique(
    const BareissEchelonResult& echelon,
    std::size_t variables,
    std::size_t rhsFirstColumn) {
    if (variables == 0)
        return {};
    if (echelon.pivotColumns.size() != variables)
        throw std::invalid_argument("Bareiss back substitution requires full column rank");
    for (std::size_t i = 0; i < variables; ++i)
        if (echelon.pivotColumns[i] != i)
            throw std::logic_error("Bareiss full-rank pivot columns are not canonical");
    if (rhsFirstColumn < variables || rhsFirstColumn > echelon.matrix.columns())
        throw std::invalid_argument("Bareiss back substitution RHS column is invalid");

    const std::size_t rhsColumns = echelon.matrix.columns() - rhsFirstColumn;
    const BigInt commonDenominator = echelon.matrix(variables - 1, variables - 1);
    if (commonDenominator.isZero())
        throw std::logic_error("Bareiss full-rank system has a zero final pivot");

    std::vector<numeric::Rational> result(variables * rhsColumns);
    std::vector<BigInt> numerators(variables);
    for (std::size_t rhs = 0; rhs < rhsColumns; ++rhs) {
        std::fill(numerators.begin(), numerators.end(), BigInt{});
        for (std::size_t row = variables; row-- > 0;) {
            BigInt numerator = echelon.matrix(row, rhsFirstColumn + rhs)
                * commonDenominator;
            for (std::size_t column = row + 1; column < variables; ++column)
                numerator -= echelon.matrix(row, column) * numerators[column];

            auto division = numeric::divmod(numerator, echelon.matrix(row, row));
            if (!division.remainder.isZero())
                throw std::logic_error("Bareiss back substitution division was not exact");
            numerators[row] = std::move(division.quotient);
        }
        for (std::size_t row = 0; row < variables; ++row)
            result[row * rhsColumns + rhs] = numeric::Rational{
                std::move(numerators[row]), commonDenominator};
    }
    return result;
}

[[nodiscard]] NumberMatrix rrefFromBareiss(BareissEchelonResult result) {
    NumberMatrix matrix{result.matrix.rows(), result.matrix.columns()};
    for (std::size_t row = 0; row < result.matrix.rows(); ++row)
        for (std::size_t column = 0; column < result.matrix.columns(); ++column)
            matrix(row, column) = Number{result.matrix(row, column)};

    for (std::size_t pivotIndex = result.pivotColumns.size(); pivotIndex-- > 0;) {
        const std::size_t row = pivotIndex;
        const std::size_t column = result.pivotColumns[pivotIndex];
        const Number pivot = matrix(row, column);
        matrix(row, column) = Number{BigInt{1}};
        for (std::size_t c = column + 1; c < matrix.columns(); ++c)
            matrix(row, c) /= pivot;

        for (std::size_t upper = 0; upper < row; ++upper) {
            const Number factor = matrix(upper, column);
            if (factor.isZero())
                continue;
            matrix(upper, column) = Number{BigInt{0}};
            for (std::size_t c = column + 1; c < matrix.columns(); ++c)
                matrix(upper, c) -= factor * matrix(row, c);
        }
    }
    return matrix;
}

[[nodiscard]] Number exactDeterminantOfRealMatrix(const MatrixView& source) {
    IntegerLift lift = liftRealMatrix(source);
    BigInt determinant = preferModularDeterminant(lift.matrix)
        ? modularDeterminant(lift.matrix)
        : bareissDeterminant(std::move(lift.matrix));
    return Number{numeric::Rational{
        std::move(determinant), rowScaleProduct(lift.rowScales)}};
}

[[nodiscard]] MatrixBuffer bareissRrefOfRealMatrix(const MatrixView& source) {
    IntegerLift lift = liftRealMatrix(source);
    auto echelon = bareissEchelon(std::move(lift.matrix), source.columns());

    // full column rankならRREFは上側の単位行列と余剰zero rowで確定する。
    // 高価なRational後退消去を行う必要はない。
    if (echelon.pivotColumns.size() == source.columns()) {
        MatrixBuffer result{source.rows(), source.columns(), integer(0)};
        for (std::size_t i = 0; i < source.columns(); ++i)
            result(i, i) = integer(1);
        return result;
    }
    return rrefFromBareiss(std::move(echelon)).toExprBuffer();
}

[[nodiscard]] std::size_t bareissRankOfRealMatrix(const MatrixView& source) {
    IntegerLift lift = liftRealMatrix(source);
    return bareissEchelon(std::move(lift.matrix), source.columns()).pivotColumns.size();
}

[[nodiscard]] Expr numericNullSpaceBasis(
    const NumberMatrix& reduced,
    const std::vector<std::size_t>& pivotColumns,
    std::size_t variables) {
    std::vector<bool> isPivot(variables, false);
    for (const std::size_t column : pivotColumns)
        isPivot[column] = true;

    const std::size_t nullity = variables - pivotColumns.size();
    std::vector<Expr> elements;
    const std::size_t shape[] = {nullity, variables};
    elements.reserve(expression::arrayElementCount(shape));

    for (std::size_t freeColumn = 0; freeColumn < variables; ++freeColumn) {
        if (isPivot[freeColumn])
            continue;

        std::vector<Number> basis(variables, Number{BigInt{0}});
        basis[freeColumn] = Number{BigInt{1}};
        for (std::size_t pivotRow = 0; pivotRow < pivotColumns.size(); ++pivotRow)
            basis[pivotColumns[pivotRow]] = -reduced(pivotRow, freeColumn);
        for (const Number& value : basis)
            elements.emplace_back(value);
    }
    return Expr::array({nullity, variables}, std::move(elements));
}

[[nodiscard]] std::vector<std::size_t> numericPivotColumns(
    const NumberMatrix& reduced,
    std::size_t variables) {
    std::vector<std::size_t> pivots;
    pivots.reserve(std::min(reduced.rows(), variables));
    for (std::size_t row = 0; row < reduced.rows(); ++row)
        for (std::size_t column = 0; column < variables; ++column)
            if (!reduced(row, column).isZero()) {
                pivots.push_back(column);
                break;
            }
    return pivots;
}

[[nodiscard]] Expr bareissNullSpaceOfRealMatrix(const MatrixView& source) {
    IntegerLift lift = liftRealMatrix(source);
    auto echelon = bareissEchelon(std::move(lift.matrix), source.columns());
    const std::vector<std::size_t> pivots = echelon.pivotColumns;
    if (pivots.size() == source.columns())
        return Expr::array({0, source.columns()}, {});
    return numericNullSpaceBasis(
        rrefFromBareiss(std::move(echelon)), pivots, source.columns());
}

[[nodiscard]] MatrixBuffer exactInverseOfRealMatrix(const MatrixView& source) {
    const std::size_t n = source.rows();
    IntegerLift lift = liftRealMatrix(source);

    const std::size_t columns = augmentedColumns(n);
    IntegerMatrixBuffer augmented{n, columns};
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column)
            augmented(row, column) = lift.matrix(row, column);
        augmented(row, n + row) = lift.rowScales[row];
    }

    auto echelon = bareissEchelon(std::move(augmented), n);
    if (echelon.pivotColumns.size() != n)
        throw std::domain_error("Matrix is singular");

    // Bareissの最終pivot（determinant up to row swaps）を全列の共通分母に使い，
    // 後退代入もBigIntだけで行う。generic Rational RREFをn本分作るより大幅に軽い。
    const auto inverse = bareissBackSubstituteUnique(echelon, n, n);
    std::vector<Expr> elements;
    const std::size_t shape[] = {n, n};
    elements.reserve(expression::arrayElementCount(shape));
    for (const numeric::Rational& value : inverse)
        elements.emplace_back(Number{value});
    return MatrixBuffer{n, n, std::move(elements)};
}

// exact complexは整数lift対象外なので、専用環を導入するまでは従来Gaussianを保持する。
[[nodiscard]] Number numericDeterminantGaussian(NumberMatrix matrix) {
    const std::size_t n = matrix.rows();
    Number result{BigInt{1}};
    bool negative = false;

    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        while (pivot < n && matrix(pivot, column).isZero())
            ++pivot;
        if (pivot == n)
            return Number{BigInt{0}};
        if (pivot != column) {
            matrix.swapRows(pivot, column);
            negative = !negative;
        }

        const Number pivotValue = matrix(column, column);
        result *= pivotValue;
        for (std::size_t row = column + 1; row < n; ++row) {
            if (matrix(row, column).isZero())
                continue;
            const Number factor = matrix(row, column) / pivotValue;
            matrix(row, column) = Number{BigInt{0}};
            for (std::size_t c = column + 1; c < n; ++c)
                matrix(row, c) -= factor * matrix(column, c);
        }
    }

    return negative ? -result : result;
}

// 同じfallbackをinverse/rref/rankで共有する。実Rational主経路はBareissへ送る。
[[nodiscard]] NumberMatrix numericRrefGaussian(
    NumberMatrix matrix,
    std::size_t pivotColumnLimit) {
    const std::size_t rows = matrix.rows();
    const std::size_t columns = matrix.columns();
    const std::size_t limit = std::min(columns, pivotColumnLimit);
    std::size_t pivotRow = 0;

    for (std::size_t column = 0; column < limit && pivotRow < rows; ++column) {
        std::size_t selected = pivotRow;
        while (selected < rows && matrix(selected, column).isZero())
            ++selected;
        if (selected == rows)
            continue;

        matrix.swapRows(selected, pivotRow);
        const Number pivot = matrix(pivotRow, column);
        for (std::size_t c = 0; c < columns; ++c)
            matrix(pivotRow, c) /= pivot;

        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivotRow || matrix(row, column).isZero())
                continue;
            const Number factor = matrix(row, column);
            matrix(row, column) = Number{BigInt{0}};
            for (std::size_t c = 0; c < columns; ++c) {
                if (c == column)
                    continue;
                matrix(row, c) -= factor * matrix(pivotRow, c);
            }
        }
        ++pivotRow;
    }
    return matrix;
}

[[nodiscard]] std::size_t numericRank(const NumberMatrix& reduced) {
    std::size_t rank = 0;
    for (std::size_t row = 0; row < reduced.rows(); ++row) {
        bool nonZero = false;
        for (std::size_t column = 0; column < reduced.columns(); ++column)
            if (!reduced(row, column).isZero()) {
                nonZero = true;
                break;
            }
        rank += nonZero ? 1 : 0;
    }
    return rank;
}

[[nodiscard]] std::size_t augmentedColumns(std::size_t columns) {
    if (columns > std::numeric_limits<std::size_t>::max() / 2)
        throw std::length_error("Augmented matrix column count exceeds the size_t range");
    return columns * 2;
}

[[nodiscard]] MatrixBuffer numericInverseGaussian(const MatrixView& source) {
    const std::size_t n = source.rows();
    const std::size_t columns = augmentedColumns(n);
    NumberMatrix augmented{n, columns};
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column)
            augmented(row, column) = source(row, column).asNumber();
        augmented(row, n + row) = Number{BigInt{1}};
    }

    NumberMatrix reduced = numericRrefGaussian(std::move(augmented), n);
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t column = 0; column < n; ++column) {
            const Number expected{BigInt{row == column ? 1 : 0}};
            if (!(reduced(row, column) == expected))
                throw std::domain_error("Matrix is singular");
        }

    std::vector<Expr> elements;
    elements.reserve(n * n);
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t column = 0; column < n; ++column)
            elements.emplace_back(reduced(row, n + column));
    return MatrixBuffer{n, n, std::move(elements)};
}

[[nodiscard]] std::optional<MatrixBuffer> symbolicRref(
    MatrixBuffer matrix,
    const ExactMatrixContext& context,
    std::size_t pivotColumnLimit);

[[nodiscard]] bool allExactRealNumbers(const expression::ArrayExpr& array) noexcept {
    return array.hasExactRealStorage();
}

[[nodiscard]] bool allExactNumbers(const expression::ArrayExpr& array) noexcept {
    return array.hasExactNumberStorage();
}

[[nodiscard]] Expr exactSolveLinearOfRealMatrix(
    const MatrixView& source,
    const expression::ArrayExpr& rhs) {
    const std::size_t rows = source.rows();
    const std::size_t variables = source.columns();
    if (variables == std::numeric_limits<std::size_t>::max())
        throw std::length_error("Linear system augmented column count exceeds the size_t range");

    IntegerLift lift = liftRealRows(rows, variables + 1,
        [&](std::size_t row, std::size_t column) {
            return column == variables
                ? rhs.exactNumber(row).asReal().toRational()
                : source.array().exactNumber(row * variables + column).asReal().toRational();
        });

    if (preferModularSolve(lift.matrix, variables)) {
        if (const auto modular = modularSolve(lift.matrix, variables)) {
            std::vector<Expr> solution;
            solution.reserve(variables);
            for (const numeric::Rational& value : *modular)
                solution.emplace_back(Number{value});
            return Expr::array({variables}, std::move(solution));
        }
    }

    auto echelon = bareissEchelon(std::move(lift.matrix), variables);
    for (std::size_t row = 0; row < rows; ++row) {
        bool coefficientNonZero = false;
        for (std::size_t column = 0; column < variables; ++column)
            if (!echelon.matrix(row, column).isZero()) {
                coefficientNonZero = true;
                break;
            }
        if (!coefficientNonZero && !echelon.matrix(row, variables).isZero())
            throw std::domain_error("Linear system is inconsistent");
    }
    if (echelon.pivotColumns.size() != variables)
        throw std::domain_error("Linear system does not have a unique solution");

    const auto solved = bareissBackSubstituteUnique(echelon, variables, variables);
    std::vector<Expr> solution;
    solution.reserve(variables);
    for (const numeric::Rational& value : solved)
        solution.emplace_back(Number{value});
    return Expr::array({variables}, std::move(solution));
}

[[nodiscard]] Expr numericSolveLinearGaussian(
    const MatrixView& source,
    const expression::ArrayExpr& rhs) {
    const std::size_t rows = source.rows();
    const std::size_t variables = source.columns();
    if (variables == std::numeric_limits<std::size_t>::max())
        throw std::length_error("Linear system augmented column count exceeds the size_t range");

    NumberMatrix augmented{rows, variables + 1};
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t column = 0; column < variables; ++column)
            augmented(row, column) = source(row, column).asNumber();
        augmented(row, variables) = rhs.exactNumber(row);
    }
    NumberMatrix reduced = numericRrefGaussian(std::move(augmented), variables);

    std::size_t rank = 0;
    for (std::size_t row = 0; row < rows; ++row) {
        bool coefficientNonZero = false;
        for (std::size_t column = 0; column < variables; ++column)
            if (!reduced(row, column).isZero()) {
                coefficientNonZero = true;
                break;
            }
        if (coefficientNonZero)
            ++rank;
        else if (!reduced(row, variables).isZero())
            throw std::domain_error("Linear system is inconsistent");
    }
    if (rank != variables)
        throw std::domain_error("Linear system does not have a unique solution");

    std::vector<Expr> solution;
    solution.reserve(variables);
    for (std::size_t variable = 0; variable < variables; ++variable)
        solution.emplace_back(reduced(variable, variables));
    return Expr::array({variables}, std::move(solution));
}

[[nodiscard]] std::optional<Expr> symbolicSolveLinear(
    const MatrixView& source,
    const expression::ArrayExpr& rhs,
    const ExactMatrixContext& context) {
    const std::size_t rows = source.rows();
    const std::size_t variables = source.columns();
    if (variables == std::numeric_limits<std::size_t>::max())
        throw std::length_error("Linear system augmented column count exceeds the size_t range");

    std::vector<Expr> elements;
    const std::size_t shape[] = {rows, variables + 1};
    elements.reserve(expression::arrayElementCount(shape));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t column = 0; column < variables; ++column)
            elements.push_back(source(row, column));
        elements.push_back(rhs.element(row));
    }
    auto reduced = symbolicRref(
        MatrixBuffer{rows, variables + 1, std::move(elements)}, context, variables);
    if (!reduced)
        return std::nullopt;

    std::size_t rank = 0;
    for (std::size_t row = 0; row < rows; ++row) {
        bool coefficientNonZero = false;
        bool coefficientUndecidable = false;
        for (std::size_t column = 0; column < variables; ++column) {
            const Expr& value = (*reduced)(row, column);
            if (exactZero(value))
                continue;
            if (provablyNonZero(value, context)) {
                coefficientNonZero = true;
                break;
            }
            coefficientUndecidable = true;
        }
        if (coefficientUndecidable && !coefficientNonZero)
            return std::nullopt;
        if (coefficientNonZero) {
            ++rank;
            continue;
        }

        const Expr& residual = (*reduced)(row, variables);
        if (exactZero(residual))
            continue;
        if (provablyNonZero(residual, context))
            throw std::domain_error("Linear system is inconsistent");
        return std::nullopt;
    }
    if (rank != variables)
        throw std::domain_error("Linear system does not have a unique solution");

    std::vector<Expr> solution;
    solution.reserve(variables);
    for (std::size_t variable = 0; variable < variables; ++variable)
        solution.push_back((*reduced)(variable, variables));
    return Expr::array({variables}, std::move(solution));
}

[[nodiscard]] MatrixBuffer minorMatrix(
    const MatrixBuffer& matrix,
    std::size_t removedRow,
    std::size_t removedColumn) {
    const std::size_t n = matrix.rows();
    std::vector<Expr> elements;
    elements.reserve((n - 1) * (n - 1));
    for (std::size_t row = 0; row < n; ++row) {
        if (row == removedRow)
            continue;
        for (std::size_t column = 0; column < n; ++column)
            if (column != removedColumn)
                elements.push_back(matrix(row, column));
    }
    return MatrixBuffer{n - 1, n - 1, std::move(elements)};
}

[[nodiscard]] std::optional<Expr> triangularDeterminant(
    const MatrixBuffer& matrix,
    const ExactMatrixContext& context) {
    const std::size_t n = matrix.rows();
    bool upper = true;
    bool lower = true;
    for (std::size_t row = 0; row < n && (upper || lower); ++row)
        for (std::size_t column = 0; column < n; ++column) {
            if (row > column && !exactZero(matrix(row, column)))
                upper = false;
            if (row < column && !exactZero(matrix(row, column)))
                lower = false;
        }
    if (!upper && !lower)
        return std::nullopt;

    Expr result = integer(1);
    for (std::size_t i = 0; i < n; ++i)
        result = multiply(std::move(result), matrix(i, i), context);
    return simplify(std::move(result), context);
}

[[nodiscard]] std::optional<Expr> symbolicDeterminant(
    const MatrixBuffer& matrix,
    const ExactMatrixContext& context,
    std::size_t& expansionBudget) {
    const std::size_t n = matrix.rows();
    if (n == 0)
        return integer(1);
    if (n == 1)
        return matrix(0, 0);
    if (n == 2)
        return subtract(
            multiply(matrix(0, 0), matrix(1, 1), context),
            multiply(matrix(0, 1), matrix(1, 0), context),
            context);
    if (const auto triangular = triangularDeterminant(matrix, context))
        return triangular;

    // Laplace展開は疎行列には有効だが、一般symbolic行列では階乗級に膨張する。
    // 全再帰で共有するbudgetを消費し、上限を超える場合は未評価のまま返す。
    std::size_t expansionRow = 0;
    std::size_t bestZeros = 0;
    for (std::size_t row = 0; row < n; ++row) {
        std::size_t zeros = 0;
        for (std::size_t column = 0; column < n; ++column)
            zeros += exactZero(matrix(row, column)) ? 1 : 0;
        if (zeros > bestZeros) {
            bestZeros = zeros;
            expansionRow = row;
        }
    }

    Expr result = integer(0);
    for (std::size_t column = 0; column < n; ++column) {
        if (exactZero(matrix(expansionRow, column)))
            continue;
        if (expansionBudget == 0)
            return std::nullopt;
        --expansionBudget;
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::TemporaryMatrixElement);

        const auto minor = symbolicDeterminant(
            minorMatrix(matrix, expansionRow, column), context, expansionBudget);
        if (!minor)
            return std::nullopt;
        Expr term = multiply(matrix(expansionRow, column), *minor, context);
        if (((expansionRow + column) & 1U) != 0)
            term = negate(std::move(term), context);
        result = add(std::move(result), std::move(term), context);
    }
    return simplify(std::move(result), context);
}

[[nodiscard]] std::optional<MatrixBuffer> symbolicRref(
    MatrixBuffer matrix,
    const ExactMatrixContext& context,
    std::size_t pivotColumnLimit) {
    const std::size_t rows = matrix.rows();
    const std::size_t columns = matrix.columns();
    const std::size_t limit = std::min(columns, pivotColumnLimit);
    std::size_t pivotRow = 0;

    for (std::size_t column = 0; column < limit && pivotRow < rows; ++column) {
        std::optional<std::size_t> selected;
        bool hasUndecidable = false;
        for (std::size_t row = pivotRow; row < rows; ++row) {
            if (exactZero(matrix(row, column)))
                continue;
            if (provablyNonZero(matrix(row, column), context)) {
                selected = row;
                break;
            }
            hasUndecidable = true;
        }
        if (!selected) {
            if (hasUndecidable)
                return std::nullopt;
            continue;
        }

        matrix.swapRows(*selected, pivotRow);
        const Expr pivot = matrix(pivotRow, column);
        for (std::size_t c = 0; c < columns; ++c)
            matrix(pivotRow, c) = divide(matrix(pivotRow, c), pivot, context);

        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivotRow || exactZero(matrix(row, column)))
                continue;
            const Expr factor = matrix(row, column);
            for (std::size_t c = 0; c < columns; ++c)
                matrix(row, c) = subtract(
                    matrix(row, c), multiply(factor, matrix(pivotRow, c), context), context);
        }
        ++pivotRow;
    }
    return matrix;
}

[[nodiscard]] std::optional<std::size_t> symbolicRank(
    const MatrixBuffer& matrix,
    const ExactMatrixContext& context) {
    std::size_t rank = 0;
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        bool nonZero = false;
        bool undecidable = false;
        for (std::size_t column = 0; column < matrix.columns(); ++column) {
            const Expr& item = matrix(row, column);
            if (exactZero(item))
                continue;
            if (provablyNonZero(item, context)) {
                nonZero = true;
                break;
            }
            undecidable = true;
        }
        if (nonZero)
            ++rank;
        else if (undecidable)
            return std::nullopt;
    }
    return rank;
}

[[nodiscard]] std::optional<std::vector<std::size_t>> symbolicPivotColumns(
    const MatrixBuffer& reduced,
    const ExactMatrixContext& context) {
    std::vector<std::size_t> pivots;
    pivots.reserve(std::min(reduced.rows(), reduced.columns()));
    for (std::size_t row = 0; row < reduced.rows(); ++row) {
        bool found = false;
        for (std::size_t column = 0; column < reduced.columns(); ++column) {
            const Expr& value = reduced(row, column);
            if (exactZero(value))
                continue;
            if (!provablyNonZero(value, context))
                return std::nullopt;
            pivots.push_back(column);
            found = true;
            break;
        }
        if (!found)
            continue;
    }
    return pivots;
}

[[nodiscard]] std::optional<Expr> symbolicNullSpace(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    auto reduced = symbolicRref(MatrixBuffer{matrix}, context, matrix.columns());
    if (!reduced)
        return std::nullopt;
    const auto pivots = symbolicPivotColumns(*reduced, context);
    if (!pivots)
        return std::nullopt;

    const std::size_t variables = matrix.columns();
    std::vector<bool> isPivot(variables, false);
    for (const std::size_t column : *pivots)
        isPivot[column] = true;

    const std::size_t nullity = variables - pivots->size();
    std::vector<Expr> elements;
    const std::size_t shape[] = {nullity, variables};
    elements.reserve(expression::arrayElementCount(shape));
    for (std::size_t freeColumn = 0; freeColumn < variables; ++freeColumn) {
        if (isPivot[freeColumn])
            continue;
        std::vector<Expr> basis(variables, integer(0));
        basis[freeColumn] = integer(1);
        for (std::size_t pivotRow = 0; pivotRow < pivots->size(); ++pivotRow)
            basis[(*pivots)[pivotRow]] = simplify(
                negate((*reduced)(pivotRow, freeColumn), context), context);
        for (Expr& value : basis)
            elements.push_back(std::move(value));
    }
    return Expr::array({nullity, variables}, std::move(elements));
}

} // namespace

bool allExactNumbers(const MatrixView& matrix) noexcept {
    return matrix.array().hasExactNumberStorage();
}

std::optional<Expr> determinant(const MatrixView& matrix, const ExactMatrixContext& context) {
    if (matrix.rows() != matrix.columns())
        throw std::invalid_argument("determinant requires a square matrix");
    if (allExactRealNumbers(matrix))
        return Expr{exactDeterminantOfRealMatrix(matrix)};
    if (allExactNumbers(matrix))
        return Expr{numericDeterminantGaussian(NumberMatrix{matrix})};

    constexpr std::size_t symbolicExpansionBudget = 512;
    std::size_t budget = symbolicExpansionBudget;
    return symbolicDeterminant(MatrixBuffer{matrix}, context, budget);
}

std::optional<Expr> inverse(const MatrixView& matrix, const ExactMatrixContext& context) {
    if (matrix.rows() != matrix.columns())
        throw std::invalid_argument("inverse requires a square matrix");
    if (matrix.rows() == 0)
        return MatrixBuffer{matrix}.toExpr();
    if (allExactRealNumbers(matrix))
        return exactInverseOfRealMatrix(matrix).toExpr();
    if (allExactNumbers(matrix))
        return numericInverseGaussian(matrix).toExpr();

    constexpr std::size_t symbolicExpansionBudget = 256;
    std::size_t budget = symbolicExpansionBudget;
    MatrixBuffer source{matrix};
    const auto det = symbolicDeterminant(source, context, budget);
    if (!det)
        return std::nullopt;
    if (exactZero(*det))
        throw std::domain_error("Matrix is singular");

    const std::size_t n = matrix.rows();
    MatrixBuffer output{n, n, integer(0)};
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) {
            const auto minor = symbolicDeterminant(
                minorMatrix(source, column, row), context, budget);
            if (!minor)
                return std::nullopt;
            Expr cofactor = *minor;
            if (((row + column) & 1U) != 0)
                cofactor = negate(std::move(cofactor), context);
            output(row, column) = divide(std::move(cofactor), *det, context);
        }
    }
    return std::move(output).toExpr();
}

std::optional<MatrixBuffer> rref(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    if (allExactRealNumbers(matrix))
        return bareissRrefOfRealMatrix(matrix);
    if (allExactNumbers(matrix))
        return numericRrefGaussian(NumberMatrix{matrix}, matrix.columns()).toExprBuffer();
    return symbolicRref(MatrixBuffer{matrix}, context, matrix.columns());
}

std::optional<std::size_t> matrixRank(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    if (allExactRealNumbers(matrix))
        return bareissRankOfRealMatrix(matrix);
    if (allExactNumbers(matrix))
        return numericRank(numericRrefGaussian(NumberMatrix{matrix}, matrix.columns()));
    const auto reduced = symbolicRref(MatrixBuffer{matrix}, context, matrix.columns());
    if (!reduced)
        return std::nullopt;
    return symbolicRank(*reduced, context);
}

std::optional<Expr> solveLinear(
    const MatrixView& matrix,
    const expression::ArrayExpr& rhs,
    const ExactMatrixContext& context) {
    if (!rhs.isVector() || rhs.shape[0] != matrix.rows())
        throw std::invalid_argument("solveLinear right-hand side size does not match matrix rows");
    if (allExactRealNumbers(matrix) && allExactRealNumbers(rhs))
        return exactSolveLinearOfRealMatrix(matrix, rhs);
    if (allExactNumbers(matrix) && allExactNumbers(rhs))
        return numericSolveLinearGaussian(matrix, rhs);
    return symbolicSolveLinear(matrix, rhs, context);
}

std::optional<Expr> nullSpace(
    const MatrixView& matrix,
    const ExactMatrixContext& context) {
    if (allExactRealNumbers(matrix))
        return bareissNullSpaceOfRealMatrix(matrix);
    if (allExactNumbers(matrix)) {
        NumberMatrix reduced = numericRrefGaussian(NumberMatrix{matrix}, matrix.columns());
        const auto pivots = numericPivotColumns(reduced, matrix.columns());
        return numericNullSpaceBasis(reduced, pivots, matrix.columns());
    }
    return symbolicNullSpace(matrix, context);
}

} // namespace mmcal::linear_algebra
