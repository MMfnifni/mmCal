#pragma once

#include "expression/array_utils.hpp"
#include "linear_algebra/matrix.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra::detail {

// exact Numberだけを密に保持する内部行列。
// elimination/decompositionごとに同じ行列bufferを持つとshape計算やswap実装が分岐するため，共通化する。
class NumberMatrix final {
public:
    NumberMatrix(std::size_t rows, std::size_t columns)
        : rows_(rows), columns_(columns) {
        const std::size_t shape[] = {rows, columns};
        values_.assign(expression::arrayElementCount(shape), numeric::Number{numeric::BigInt{0}});
    }

    explicit NumberMatrix(const MatrixView& source)
        : NumberMatrix(source.rows(), source.columns()) {
        for (std::size_t i = 0; i < values_.size(); ++i)
            values_[i] = source.array().exactNumber(i);
    }

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }

    [[nodiscard]] numeric::Number& operator()(std::size_t row, std::size_t column) noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const numeric::Number& operator()(
        std::size_t row,
        std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }

    void swapRows(std::size_t lhs, std::size_t rhs) noexcept {
        if (lhs == rhs)
            return;
        for (std::size_t column = 0; column < columns_; ++column)
            std::swap((*this)(lhs, column), (*this)(rhs, column));
    }

    [[nodiscard]] const std::vector<numeric::Number>& values() const noexcept {
        return values_;
    }

    [[nodiscard]] MatrixBuffer toExprBuffer() const {
        std::vector<expression::Expr> elements;
        elements.reserve(values_.size());
        for (const numeric::Number& value : values_)
            elements.emplace_back(value);
        return MatrixBuffer{rows_, columns_, std::move(elements)};
    }

private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<numeric::Number> values_;
};

} // namespace mmcal::linear_algebra::detail
