#pragma once

#include "expression/expr.hpp"

#include <cstddef>
#include <vector>

namespace mmcal::linear_algebra {

// ArrayExprのrank-2 row-major storageをコピーせず参照する軽量view。
class MatrixView final {
public:
    explicit MatrixView(const expression::ArrayExpr& array);

    [[nodiscard]] std::size_t rows() const noexcept;
    [[nodiscard]] std::size_t columns() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] expression::Expr operator()(
        std::size_t row, std::size_t column) const;
    [[nodiscard]] const expression::ArrayExpr& array() const noexcept;

private:
    const expression::ArrayExpr* array_ = nullptr;
};

// elimination等で書換えが必要な場合だけ一度flat copyする作業領域。
class MatrixBuffer final {
public:
    MatrixBuffer(std::size_t rows, std::size_t columns, expression::Expr initialValue);
    explicit MatrixBuffer(const MatrixView& source);
    MatrixBuffer(std::size_t rows, std::size_t columns, std::vector<expression::Expr> elements);

    [[nodiscard]] std::size_t rows() const noexcept;
    [[nodiscard]] std::size_t columns() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;

    [[nodiscard]] expression::Expr& operator()(std::size_t row, std::size_t column) noexcept;
    [[nodiscard]] const expression::Expr& operator()(std::size_t row, std::size_t column) const noexcept;

    void swapRows(std::size_t lhs, std::size_t rhs) noexcept;

    [[nodiscard]] std::vector<expression::Expr>& elements() noexcept;
    [[nodiscard]] const std::vector<expression::Expr>& elements() const noexcept;
    [[nodiscard]] expression::Expr toExpr() const &;
    [[nodiscard]] expression::Expr toExpr() &&;

private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<expression::Expr> elements_;
};

} // namespace mmcal::linear_algebra
