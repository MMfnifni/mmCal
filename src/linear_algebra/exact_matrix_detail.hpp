#pragma once

#include "builtins/exact_operations.hpp"
#include "expression/exact_value.hpp"
#include "linear_algebra/exact_matrix.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <cstdint>
#include <utility>

namespace mmcal::linear_algebra::detail {

[[nodiscard]] inline expression::Expr integer(std::int64_t value) {
    return expression::exact::integer(value);
}

[[nodiscard]] inline bool exactZero(const expression::Expr& value) {
    return expression::exact::isZero(value);
}

[[nodiscard]] inline bool provablyNonZero(
    const expression::Expr& value,
    const ExactMatrixContext& context) {
    if (value.isNumber())
        return !value.asNumber().isZero();
    const auto facts = mathematics::inferValueFacts(
        value, context.builtins, context.mathematics);
    return facts.sign == mathematics::RealSign::Positive
        || facts.sign == mathematics::RealSign::Negative
        || facts.sign == mathematics::RealSign::NonZero;
}

[[nodiscard]] inline expression::Expr add(
    expression::Expr lhs,
    expression::Expr rhs,
    const ExactMatrixContext& context) {
    return builtins::exact::add({std::move(lhs), std::move(rhs)},
        context.builtins, context.mathematics, context.angles);
}

[[nodiscard]] inline expression::Expr subtract(
    expression::Expr lhs,
    expression::Expr rhs,
    const ExactMatrixContext& context) {
    return builtins::exact::subtract(std::move(lhs), std::move(rhs),
        context.builtins, context.mathematics, context.angles);
}

[[nodiscard]] inline expression::Expr multiply(
    expression::Expr lhs,
    expression::Expr rhs,
    const ExactMatrixContext& context) {
    return builtins::exact::multiply({std::move(lhs), std::move(rhs)},
        context.builtins, context.mathematics, context.angles);
}

[[nodiscard]] inline expression::Expr divide(
    expression::Expr lhs,
    expression::Expr rhs,
    const ExactMatrixContext& context) {
    return builtins::exact::divide(std::move(lhs), std::move(rhs),
        context.builtins, context.mathematics, context.angles);
}

[[nodiscard]] inline expression::Expr negate(
    expression::Expr value,
    const ExactMatrixContext& context) {
    return builtins::exact::negate(std::move(value),
        context.builtins, context.mathematics, context.angles);
}

[[nodiscard]] inline expression::Expr simplify(
    expression::Expr value,
    const ExactMatrixContext& context) {
    return builtins::exact::simplify(std::move(value),
        context.builtins, context.mathematics, context.angles);
}

[[nodiscard]] inline expression::Expr squareRoot(
    expression::Expr value,
    const ExactMatrixContext& context) {
    return builtins::exact::sqrt(std::move(value),
        context.builtins, context.mathematics, context.angles);
}

} // namespace mmcal::linear_algebra::detail
