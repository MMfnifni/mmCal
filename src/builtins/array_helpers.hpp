// 配列shape検査の共通補助
#pragma once

#include "error/error_message.hpp"
#include "expression/expr.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>

namespace mmcal::builtins::detail {

[[noreturn]] inline void arrayTypeError(std::string message) {
    error::throwCalcError(error::CalcErrorType::Type, std::move(message));
}

[[nodiscard]] inline const expression::ArrayExpr& requireArray(
    const expression::Expr& expression, std::string_view name) {
    if (!expression.isArray())
        arrayTypeError(std::string{name} + " expects an array");
    return expression.asArray();
}

[[nodiscard]] inline const expression::ArrayExpr& requireRank(
    const expression::Expr& expression, std::size_t rank, std::string_view name) {
    const expression::ArrayExpr& array = requireArray(expression, name);
    if (array.rank() != rank)
        arrayTypeError(std::string{name} + " expects a rank-" + std::to_string(rank) + " array");
    return array;
}

[[nodiscard]] inline const expression::ArrayExpr& requireVector(
    const expression::Expr& expression, std::string_view name) {
    return requireRank(expression, 1, name);
}

[[nodiscard]] inline const expression::ArrayExpr& requireMatrix(
    const expression::Expr& expression, std::string_view name) {
    return requireRank(expression, 2, name);
}

[[nodiscard]] constexpr std::size_t matrixIndex(
    std::size_t row, std::size_t column, std::size_t columns) noexcept {
    return row * columns + column;
}

} // namespace mmcal::builtins::detail
