// 行列・線形代数
#include "linear_algebra.hpp"

#include "builtins/array_helpers.hpp"
#include "builtins/exact_operations.hpp"
#include "error/error_message.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using evaluation::BuiltinId;
using expression::ArrayExpr;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr simplify(
    Expr expression,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return exact::simplify(std::move(expression), registry, mathematics, angles);
}

[[nodiscard]] Expr add(
    Expr lhs, Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return exact::add({std::move(lhs), std::move(rhs)}, registry, mathematics, angles);
}

[[nodiscard]] Expr subtract(
    Expr lhs, Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return exact::subtract(std::move(lhs), std::move(rhs), registry, mathematics, angles);
}

[[nodiscard]] Expr multiply(
    Expr lhs, Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return exact::multiply({std::move(lhs), std::move(rhs)}, registry, mathematics, angles);
}

[[nodiscard]] Expr divide(
    Expr lhs, Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return exact::divide(std::move(lhs), std::move(rhs), registry, mathematics, angles);
}

[[nodiscard]] Expr negate(
    Expr value,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return exact::negate(std::move(value), registry, mathematics, angles);
}

[[nodiscard]] bool exactZero(const Expr& value) {
    return value.isNumber() && value.asNumber().isZero();
}

[[nodiscard]] bool provablyNonZero(
    const Expr& value,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    if (value.isNumber())
        return !value.asNumber().isZero();
    const auto facts = mathematics::inferValueFacts(value, registry, mathematics);
    return facts.sign == mathematics::RealSign::Positive
        || facts.sign == mathematics::RealSign::Negative
        || facts.sign == mathematics::RealSign::NonZero;
}

[[nodiscard]] std::vector<std::vector<Expr>> rowsOf(const ArrayExpr& matrix) {
    const std::size_t rows = matrix.shape[0];
    const std::size_t columns = matrix.shape[1];
    std::vector<std::vector<Expr>> result(rows, std::vector<Expr>{});
    for (std::size_t r = 0; r < rows; ++r) {
        result[r].reserve(columns);
        for (std::size_t c = 0; c < columns; ++c)
            result[r].push_back(matrix.elements[detail::matrixIndex(r, c, columns)]);
    }
    return result;
}

[[nodiscard]] Expr matrixExpr(const std::vector<std::vector<Expr>>& rows) {
    const std::size_t rowCount = rows.size();
    const std::size_t columnCount = rowCount == 0 ? 0 : rows.front().size();
    std::vector<Expr> elements;
    elements.reserve(rowCount * columnCount);
    for (const auto& row : rows) {
        if (row.size() != columnCount)
            throw std::logic_error("Internal matrix rows have inconsistent lengths");
        elements.insert(elements.end(), row.begin(), row.end());
    }
    return Expr::array({rowCount, columnCount}, std::move(elements));
}

[[nodiscard]] bool allNumbers(const std::vector<std::vector<Expr>>& matrix) {
    for (const auto& row : matrix)
        for (const Expr& item : row)
            if (!item.isNumber())
                return false;
    return true;
}

[[nodiscard]] Expr numericDeterminant(std::vector<std::vector<Expr>> matrix) {
    const std::size_t n = matrix.size();
    if (n == 0)
        return integer(1);

    Number determinant{BigInt{1}};
    bool negative = false;
    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        while (pivot < n && matrix[pivot][column].asNumber().isZero())
            ++pivot;
        if (pivot == n)
            return integer(0);
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            negative = !negative;
        }

        const Number pivotValue = matrix[column][column].asNumber();
        determinant *= pivotValue;
        for (std::size_t row = column + 1; row < n; ++row) {
            if (matrix[row][column].asNumber().isZero())
                continue;
            const Number factor = matrix[row][column].asNumber() / pivotValue;
            for (std::size_t c = column + 1; c < n; ++c)
                matrix[row][c] = Expr{matrix[row][c].asNumber() - factor * matrix[column][c].asNumber()};
            matrix[row][column] = integer(0);
        }
    }
    if (negative)
        determinant = -determinant;
    return Expr{std::move(determinant)};
}

[[nodiscard]] std::vector<std::vector<Expr>> minorMatrix(
    const std::vector<std::vector<Expr>>& matrix,
    std::size_t removedRow,
    std::size_t removedColumn) {
    std::vector<std::vector<Expr>> result;
    result.reserve(matrix.size() - 1);
    for (std::size_t r = 0; r < matrix.size(); ++r) {
        if (r == removedRow)
            continue;
        std::vector<Expr> row;
        row.reserve(matrix.size() - 1);
        for (std::size_t c = 0; c < matrix.size(); ++c)
            if (c != removedColumn)
                row.push_back(matrix[r][c]);
        result.push_back(std::move(row));
    }
    return result;
}

[[nodiscard]] Expr determinantOf(
    const std::vector<std::vector<Expr>>& matrix,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::size_t n = matrix.size();
    if (n == 0)
        return integer(1);
    if (n == 1)
        return matrix[0][0];
    if (allNumbers(matrix))
        return numericDeterminant(matrix);
    if (n == 2)
        return subtract(
            multiply(matrix[0][0], matrix[1][1], registry, mathematics, angles),
            multiply(matrix[0][1], matrix[1][0], registry, mathematics, angles),
            registry, mathematics, angles);

    // Laplace展開では、exact zeroを最も多く含む行を選んで式膨張を抑える。
    std::size_t expansionRow = 0;
    std::size_t bestZeros = 0;
    for (std::size_t r = 0; r < n; ++r) {
        const std::size_t zeros = static_cast<std::size_t>(std::count_if(
            matrix[r].begin(), matrix[r].end(), exactZero));
        if (zeros > bestZeros) {
            bestZeros = zeros;
            expansionRow = r;
        }
    }

    Expr result = integer(0);
    for (std::size_t c = 0; c < n; ++c) {
        if (exactZero(matrix[expansionRow][c]))
            continue;
        Expr term = multiply(
            matrix[expansionRow][c],
            determinantOf(minorMatrix(matrix, expansionRow, c), registry, mathematics, angles),
            registry, mathematics, angles);
        if (((expansionRow + c) & 1U) != 0)
            term = negate(std::move(term), registry, mathematics, angles);
        result = add(std::move(result), std::move(term), registry, mathematics, angles);
    }
    return simplify(std::move(result), registry, mathematics, angles);
}

[[nodiscard]] std::optional<std::vector<std::vector<Expr>>> gaussianRref(
    std::vector<std::vector<Expr>> matrix,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (matrix.empty())
        return matrix;
    const std::size_t rows = matrix.size();
    const std::size_t columns = matrix.front().size();
    std::size_t pivotRow = 0;

    for (std::size_t column = 0; column < columns && pivotRow < rows; ++column) {
        std::optional<std::size_t> selected;
        bool hasUndecidable = false;
        for (std::size_t row = pivotRow; row < rows; ++row) {
            if (exactZero(matrix[row][column]))
                continue;
            if (provablyNonZero(matrix[row][column], registry, mathematics)) {
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

        if (*selected != pivotRow)
            std::swap(matrix[*selected], matrix[pivotRow]);

        const Expr pivot = matrix[pivotRow][column];
        for (std::size_t c = 0; c < columns; ++c)
            matrix[pivotRow][c] = divide(matrix[pivotRow][c], pivot, registry, mathematics, angles);

        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivotRow || exactZero(matrix[row][column]))
                continue;
            const Expr factor = matrix[row][column];
            for (std::size_t c = 0; c < columns; ++c) {
                matrix[row][c] = subtract(
                    matrix[row][c],
                    multiply(factor, matrix[pivotRow][c], registry, mathematics, angles),
                    registry, mathematics, angles);
            }
        }
        ++pivotRow;
    }
    return matrix;
}

[[nodiscard]] std::optional<std::size_t> rankOfRref(
    const std::vector<std::vector<Expr>>& matrix,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    std::size_t rank = 0;
    for (const auto& row : matrix) {
        bool nonZero = false;
        bool undecidable = false;
        for (const Expr& item : row) {
            if (exactZero(item))
                continue;
            if (provablyNonZero(item, registry, mathematics)) {
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

} // namespace

Expr evaluateTranspose(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry&) {
    if (arguments.size() != 1)
        detail::arrayTypeError("transpose expects one array");
    const ArrayExpr& array = detail::requireArray(arguments.front(), "transpose");
    if (array.rank() == 1)
        return arguments.front();
    if (array.rank() != 2)
        detail::arrayTypeError("transpose currently supports rank-1 or rank-2 arrays");

    const std::size_t rows = array.shape[0];
    const std::size_t columns = array.shape[1];
    std::vector<Expr> elements;
    elements.reserve(array.elements.size());
    for (std::size_t c = 0; c < columns; ++c)
        for (std::size_t r = 0; r < rows; ++r)
            elements.push_back(array.elements[detail::matrixIndex(r, c, columns)]);
    return Expr::array({columns, rows}, std::move(elements));
}

Expr evaluateMatrixAdd(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() < 2)
        detail::arrayTypeError("madd expects at least two arrays");
    const ArrayExpr& first = detail::requireArray(arguments.front(), "madd");
    std::vector<Expr> result = first.elements;
    for (std::size_t a = 1; a < arguments.size(); ++a) {
        const ArrayExpr& next = detail::requireArray(arguments[a], "madd");
        if (next.shape != first.shape)
            error::throwCalcError(error::CalcErrorType::Domain, "madd requires identical array shapes");
        for (std::size_t i = 0; i < result.size(); ++i)
            result[i] = add(result[i], next.elements[i], registry, mathematics, angles);
    }
    return Expr::array(first.shape, std::move(result));
}

Expr evaluateMatrixMultiply(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() != 2)
        detail::arrayTypeError("matmul expects two arrays");
    const ArrayExpr& lhs = detail::requireArray(arguments[0], "matmul");
    const ArrayExpr& rhs = detail::requireArray(arguments[1], "matmul");
    if (lhs.rank() < 1 || lhs.rank() > 2 || rhs.rank() < 1 || rhs.rank() > 2)
        detail::arrayTypeError("matmul currently supports vectors and matrices");

    const std::size_t lhsRows = lhs.rank() == 1 ? 1 : lhs.shape[0];
    const std::size_t lhsColumns = lhs.rank() == 1 ? lhs.shape[0] : lhs.shape[1];
    const std::size_t rhsRows = rhs.rank() == 1 ? rhs.shape[0] : rhs.shape[0];
    const std::size_t rhsColumns = rhs.rank() == 1 ? 1 : rhs.shape[1];
    if (lhsColumns != rhsRows)
        error::throwCalcError(error::CalcErrorType::Domain, "matmul inner dimensions do not agree");

    const auto lhsAt = [&](std::size_t r, std::size_t c) -> const Expr& {
        return lhs.rank() == 1 ? lhs.elements[c] : lhs.elements[detail::matrixIndex(r, c, lhsColumns)];
    };
    const auto rhsAt = [&](std::size_t r, std::size_t c) -> const Expr& {
        return rhs.rank() == 1 ? rhs.elements[r] : rhs.elements[detail::matrixIndex(r, c, rhsColumns)];
    };

    std::vector<Expr> elements;
    elements.reserve(lhsRows * rhsColumns);
    for (std::size_t r = 0; r < lhsRows; ++r) {
        for (std::size_t c = 0; c < rhsColumns; ++c) {
            Expr sum = integer(0);
            for (std::size_t k = 0; k < lhsColumns; ++k)
                sum = add(std::move(sum), multiply(lhsAt(r, k), rhsAt(k, c), registry, mathematics, angles),
                    registry, mathematics, angles);
            elements.push_back(std::move(sum));
        }
    }

    if (lhs.rank() == 1 && rhs.rank() == 1)
        return elements.front();
    if (lhs.rank() == 1 || rhs.rank() == 1)
        return Expr::array({lhs.rank() == 1 ? rhsColumns : lhsRows}, std::move(elements));
    return Expr::array({lhsRows, rhsColumns}, std::move(elements));
}

Expr evaluateDeterminant(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() != 1)
        detail::arrayTypeError("det expects one matrix");
    const ArrayExpr& matrix = detail::requireMatrix(arguments.front(), "det");
    if (matrix.shape[0] != matrix.shape[1])
        error::throwCalcError(error::CalcErrorType::Domain, "det requires a square matrix");
    return determinantOf(rowsOf(matrix), registry, mathematics, angles);
}

Expr evaluateMatrixInverse(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() != 1)
        detail::arrayTypeError("inverse expects one matrix");
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "inverse");
    if (array.shape[0] != array.shape[1])
        error::throwCalcError(error::CalcErrorType::Domain, "inverse requires a square matrix");
    const std::size_t n = array.shape[0];
    if (n == 0)
        return arguments.front();

    const auto matrix = rowsOf(array);
    const Expr determinant = determinantOf(matrix, registry, mathematics, angles);
    if (exactZero(determinant))
        error::throwCalcError(error::CalcErrorType::Domain, "Matrix is singular");

    // 数値行列はGauss-JordanでO(n^3)。symbolic行列はadjugate/detにして、determinant != 0 というinverse自身の定義域を各Divideへ保持する。
    if (allNumbers(matrix)) {
        std::vector<std::vector<Expr>> augmented(n, std::vector<Expr>(2 * n, integer(0)));
        for (std::size_t r = 0; r < n; ++r) {
            for (std::size_t c = 0; c < n; ++c)
                augmented[r][c] = matrix[r][c];
            augmented[r][n + r] = integer(1);
        }
        const auto reduced = gaussianRref(std::move(augmented), registry, mathematics, angles);
        if (!reduced)
            error::throwCalcError(error::CalcErrorType::Internal, "Exact numeric inverse pivoting failed");
        std::vector<std::vector<Expr>> output(n, std::vector<Expr>(n, integer(0)));
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t c = 0; c < n; ++c)
                output[r][c] = (*reduced)[r][n + c];
        return matrixExpr(output);
    }

    std::vector<std::vector<Expr>> output(n, std::vector<Expr>(n, integer(0)));
    for (std::size_t r = 0; r < n; ++r) {
        for (std::size_t c = 0; c < n; ++c) {
            Expr cofactor = determinantOf(minorMatrix(matrix, c, r), registry, mathematics, angles);
            if (((r + c) & 1U) != 0)
                cofactor = negate(std::move(cofactor), registry, mathematics, angles);
            output[r][c] = divide(std::move(cofactor), determinant, registry, mathematics, angles);
        }
    }
    return matrixExpr(output);
}

Expr evaluateRref(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() != 1)
        detail::arrayTypeError("rref expects one matrix");
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "rref");
    const auto reduced = gaussianRref(rowsOf(array), registry, mathematics, angles);
    if (!reduced)
        return Expr::call(registry.symbol(BuiltinId::Rref), {arguments.front()});
    return matrixExpr(*reduced);
}

Expr evaluateMatrixRank(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() != 1)
        detail::arrayTypeError("rank expects one matrix");
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "rank");
    const auto reduced = gaussianRref(rowsOf(array), registry, mathematics, angles);
    if (!reduced)
        return Expr::call(registry.symbol(BuiltinId::Rank), {arguments.front()});
    const auto rank = rankOfRref(*reduced, registry, mathematics);
    if (!rank)
        return Expr::call(registry.symbol(BuiltinId::Rank), {arguments.front()});
    return Expr{Number{BigInt::parse(std::to_string(*rank))}};
}

} // namespace mmcal::builtins
