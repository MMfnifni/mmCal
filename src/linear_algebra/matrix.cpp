// rank-2 Arrayのrow-major viewと作業buffer
#include "matrix.hpp"

#include "expression/array_utils.hpp"
#include "evaluation/evaluation_budget.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace mmcal::linear_algebra {

MatrixView::MatrixView(const expression::ArrayExpr& array)
    : array_(&array) {
    if (!array.isMatrix())
        throw std::invalid_argument("MatrixView requires a rank-2 array");
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::DenseArrayElement, array.size());
}

std::size_t MatrixView::rows() const noexcept {
    return array_->shape[0];
}

std::size_t MatrixView::columns() const noexcept {
    return array_->shape[1];
}

std::size_t MatrixView::size() const noexcept {
    return array_->size();
}

expression::Expr MatrixView::operator()(
    std::size_t row, std::size_t column) const {
    return array_->element(row * columns() + column);
}

const expression::ArrayExpr& MatrixView::array() const noexcept {
    return *array_;
}

MatrixBuffer::MatrixBuffer(
    std::size_t rows,
    std::size_t columns,
    expression::Expr initialValue)
    : rows_(rows), columns_(columns) {
    const std::size_t shape[] = {rows, columns};
    const std::size_t count = expression::arrayElementCount(shape);
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::TemporaryMatrixElement, count);
    elements_.assign(count, std::move(initialValue));
}

MatrixBuffer::MatrixBuffer(const MatrixView& source)
    : rows_(source.rows()), columns_(source.columns()) {
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::TemporaryMatrixElement, source.size());
    elements_ = source.array().materialize();
}

MatrixBuffer::MatrixBuffer(
    std::size_t rows,
    std::size_t columns,
    std::vector<expression::Expr> elements)
    : rows_(rows), columns_(columns), elements_(std::move(elements)) {
    const std::size_t shape[] = {rows, columns};
    if (expression::arrayElementCount(shape) != elements_.size())
        throw std::invalid_argument("MatrixBuffer shape does not match the element count");
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::TemporaryMatrixElement, elements_.size());
    for (const expression::Expr& element : elements_)
        if (element.isArray())
            throw std::invalid_argument("MatrixBuffer elements must be scalar expressions");
}

std::size_t MatrixBuffer::rows() const noexcept {
    return rows_;
}

std::size_t MatrixBuffer::columns() const noexcept {
    return columns_;
}

std::size_t MatrixBuffer::size() const noexcept {
    return elements_.size();
}

expression::Expr& MatrixBuffer::operator()(
    std::size_t row, std::size_t column) noexcept {
    return elements_[row * columns_ + column];
}

const expression::Expr& MatrixBuffer::operator()(
    std::size_t row, std::size_t column) const noexcept {
    return elements_[row * columns_ + column];
}

void MatrixBuffer::swapRows(std::size_t lhs, std::size_t rhs) noexcept {
    if (lhs == rhs)
        return;
    for (std::size_t column = 0; column < columns_; ++column)
        std::swap((*this)(lhs, column), (*this)(rhs, column));
}

std::vector<expression::Expr>& MatrixBuffer::elements() noexcept {
    return elements_;
}

const std::vector<expression::Expr>& MatrixBuffer::elements() const noexcept {
    return elements_;
}

expression::Expr MatrixBuffer::toExpr() const & {
    return expression::Expr::array({rows_, columns_}, elements_);
}

expression::Expr MatrixBuffer::toExpr() && {
    return expression::Expr::array({rows_, columns_}, std::move(elements_));
}

} // namespace mmcal::linear_algebra
