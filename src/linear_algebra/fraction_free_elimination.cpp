// Bareiss fraction-free elimination。中間Rational生成を避け、BigInt上で厳密消去する。
#include "fraction_free_elimination.hpp"

#include "expression/array_utils.hpp"
#include "evaluation/evaluation_budget.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace mmcal::linear_algebra {
namespace {

using numeric::BigInt;

[[nodiscard]] BigInt exactDivide(BigInt numerator, const BigInt& denominator) {
    if (denominator == BigInt{1})
        return numerator;
    if (denominator == BigInt{-1})
        return -numerator;

    auto result = numeric::divmod(numerator, denominator);
    if (!result.remainder.isZero())
        throw std::logic_error("Bareiss division was not exact");
    return std::move(result.quotient);
}

[[nodiscard]] std::size_t selectPivot(
    const IntegerMatrixBuffer& matrix,
    std::size_t firstRow,
    std::size_t column) {
    std::size_t selected = matrix.rows();
    std::size_t bestBits = 0;
    for (std::size_t row = firstRow; row < matrix.rows(); ++row) {
        const BigInt& value = matrix(row, column);
        if (value.isZero())
            continue;
        const std::size_t bits = value.bitLength();
        if (selected == matrix.rows() || bits < bestBits) {
            selected = row;
            bestBits = bits;
        }
    }
    return selected;
}

} // namespace

IntegerMatrixBuffer::IntegerMatrixBuffer(std::size_t rows, std::size_t columns)
    : rows_(rows), columns_(columns) {
    const std::size_t shape[] = {rows, columns};
    elements_.assign(expression::arrayElementCount(shape), BigInt{});
}

IntegerMatrixBuffer::IntegerMatrixBuffer(
    std::size_t rows,
    std::size_t columns,
    std::vector<BigInt> elements)
    : rows_(rows), columns_(columns), elements_(std::move(elements)) {
    const std::size_t shape[] = {rows, columns};
    if (elements_.size() != expression::arrayElementCount(shape))
        throw std::invalid_argument("Integer matrix element count does not match shape");
}

std::size_t IntegerMatrixBuffer::rows() const noexcept { return rows_; }
std::size_t IntegerMatrixBuffer::columns() const noexcept { return columns_; }
std::size_t IntegerMatrixBuffer::size() const noexcept { return elements_.size(); }

BigInt& IntegerMatrixBuffer::operator()(std::size_t row, std::size_t column) noexcept {
    return elements_[row * columns_ + column];
}

const BigInt& IntegerMatrixBuffer::operator()(
    std::size_t row, std::size_t column) const noexcept {
    return elements_[row * columns_ + column];
}

void IntegerMatrixBuffer::swapRows(std::size_t lhs, std::size_t rhs) noexcept {
    if (lhs == rhs)
        return;
    for (std::size_t column = 0; column < columns_; ++column)
        std::swap((*this)(lhs, column), (*this)(rhs, column));
}

std::vector<BigInt>& IntegerMatrixBuffer::elements() noexcept { return elements_; }
const std::vector<BigInt>& IntegerMatrixBuffer::elements() const noexcept { return elements_; }

BareissEchelonResult bareissEchelon(
    IntegerMatrixBuffer matrix,
    std::size_t pivotColumnLimit) {
    pivotColumnLimit = std::min(pivotColumnLimit, matrix.columns());
    std::vector<std::size_t> pivotColumns;
    pivotColumns.reserve(std::min(matrix.rows(), pivotColumnLimit));

    std::size_t pivotRow = 0;
    std::size_t rowSwaps = 0;
    BigInt previousPivot{1};

    for (std::size_t column = 0;
         column < pivotColumnLimit && pivotRow < matrix.rows();
         ++column) {
        evaluation::checkEvaluationCancellation();
        const std::size_t selected = selectPivot(matrix, pivotRow, column);
        if (selected == matrix.rows())
            continue;
        if (selected != pivotRow) {
            matrix.swapRows(selected, pivotRow);
            ++rowSwaps;
        }

        const BigInt pivot = matrix(pivotRow, column);
        for (std::size_t row = pivotRow + 1; row < matrix.rows(); ++row) {
            evaluation::checkEvaluationCancellation();
            const BigInt factor = matrix(row, column);
            for (std::size_t c = column + 1; c < matrix.columns(); ++c) {
                BigInt numerator = pivot * matrix(row, c);
                if (!factor.isZero())
                    numerator -= factor * matrix(pivotRow, c);
                matrix(row, c) = exactDivide(std::move(numerator), previousPivot);
            }
            matrix(row, column) = BigInt{};
        }

        previousPivot = pivot;
        pivotColumns.push_back(column);
        ++pivotRow;
    }

    return {std::move(matrix), std::move(pivotColumns), rowSwaps};
}

BigInt bareissDeterminant(IntegerMatrixBuffer matrix) {
    if (matrix.rows() != matrix.columns())
        throw std::invalid_argument("Bareiss determinant requires a square matrix");
    const std::size_t n = matrix.rows();
    if (n == 0)
        return BigInt{1};
    if (n == 1)
        return matrix(0, 0);

    auto result = bareissEchelon(std::move(matrix), n);
    if (result.pivotColumns.size() != n)
        return BigInt{};

    BigInt determinant = result.matrix(n - 1, n - 1);
    if ((result.rowSwaps & 1U) != 0)
        determinant = -determinant;
    return determinant;
}

} // namespace mmcal::linear_algebra
