#pragma once

#include "expression/expr.hpp"
#include "expression/symbol.hpp"

#include <optional>

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
    if (array.shape.size() != 1 || array.shape[0] != 3 || array.elements.size() != 3)
        return nullptr;
    return &array;
}

[[nodiscard]] inline std::optional<RangeIteratorSpec> parseRangeIteratorSpec(
    const expression::Expr& expression) {
    const expression::ArrayExpr* array = rangeIteratorArray(expression);
    if (!array || !array->elements[0].isSymbol())
        return std::nullopt;

    return RangeIteratorSpec{
        array->elements[0].asSymbol(),
        array->elements[1],
        array->elements[2]
    };
}

} // namespace mmcal::evaluation
