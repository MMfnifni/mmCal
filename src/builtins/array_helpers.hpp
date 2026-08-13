// 配列shape・index検査の共通補助
#pragma once

#include "error/error_message.hpp"
#include "expression/array_utils.hpp"
#include "expression/expr.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <charconv>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
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

[[nodiscard]] inline std::optional<std::size_t> exactSizeValue(
    const expression::Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;

    const numeric::BigInt& value = expression.asNumber().asReal().asInteger();
    if (value.isNegative())
        return std::nullopt;

    const std::string text = value.toString();
    std::size_t result = 0;
    const auto converted = std::from_chars(text.data(), text.data() + text.size(), result);
    if (converted.ec != std::errc{} || converted.ptr != text.data() + text.size())
        return std::nullopt;
    return result;
}

[[nodiscard]] inline std::size_t requireSize(
    const expression::Expr& expression, std::string_view name) {
    const auto value = exactSizeValue(expression);
    if (!value)
        arrayTypeError(std::string{name} + " requires a nonnegative machine-size integer");
    return *value;
}

[[nodiscard]] inline expression::Expr sizeExpr(std::size_t value) {
    return expression::Expr{numeric::Number{
        numeric::BigInt::parse(std::to_string(value))}};
}

[[nodiscard]] constexpr std::size_t matrixIndex(
    std::size_t row, std::size_t column, std::size_t columns) noexcept {
    return row * columns + column;
}

} // namespace mmcal::builtins::detail
