#pragma once

#include "expression/expr.hpp"
#include "expression/symbol.hpp"

#include <optional>
#include <vector>

namespace mmcal::evaluation {

// {variable, lower, upper} 形式の有限区間iterator。
// integrate/nintegrate/sum/product等のbinderで同じ構文規約を共有する。
struct RangeIteratorSpec final {
    expression::Symbol variable;
    expression::Expr lower;
    expression::Expr upper;
};

[[nodiscard]] inline const expression::ArrayExpr* rangeIteratorArray(
    const expression::Expr& expression) noexcept {
    if (!expression.isArray())
        return nullptr;

    const expression::ArrayExpr& array = expression.asArray();
    if (array.shape.size() != 1 || array.shape[0] != 3 || array.size() != 3)
        return nullptr;
    return &array;
}

[[nodiscard]] inline std::optional<RangeIteratorSpec> parseRangeIteratorSpec(
    const expression::Expr& expression) {
    const expression::ArrayExpr* array = rangeIteratorArray(expression);
    if (!array)
        return std::nullopt;

    const expression::Expr variable = array->element(0);
    if (!variable.isSymbol())
        return std::nullopt;

    return RangeIteratorSpec{
        variable.asSymbol(),
        array->element(1),
        array->element(2)
    };
}

// table専用iterator。{i,n} / {i,lower,upper} / {i,lower,upper,step}を受ける。
struct TableIteratorSpec final {
    expression::Symbol variable;
    std::vector<expression::Expr> rangeArguments;
};

[[nodiscard]] inline const expression::ArrayExpr* tableIteratorArray(
    const expression::Expr& expression) noexcept {
    if (!expression.isArray())
        return nullptr;

    const expression::ArrayExpr& array = expression.asArray();
    if (array.rank() != 1 || array.size() < 2 || array.size() > 4)
        return nullptr;
    return &array;
}

[[nodiscard]] inline std::optional<TableIteratorSpec> parseTableIteratorSpec(
    const expression::Expr& expression) {
    const expression::ArrayExpr* array = tableIteratorArray(expression);
    if (!array)
        return std::nullopt;

    const expression::Expr variable = array->element(0);
    if (!variable.isSymbol())
        return std::nullopt;

    TableIteratorSpec result;
    result.variable = variable.asSymbol();
    result.rangeArguments.reserve(array->size() - 1);
    for (std::size_t i = 1; i < array->size(); ++i)
        result.rangeArguments.push_back(array->element(i));
    return result;
}

} // namespace mmcal::evaluation
