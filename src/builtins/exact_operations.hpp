// 組込み函数向けexact演算補助
#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"

#include <utility>
#include <vector>

namespace mmcal::builtins::exact {

[[nodiscard]] inline expression::Expr simplify(
    expression::Expr expression,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{registry, mathematics, angles});
}

[[nodiscard]] inline expression::Expr call(
    evaluation::BuiltinId id,
    std::vector<expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplify(expression::Expr::call(registry.symbol(id), std::move(arguments)),
        registry, mathematics, angles);
}

[[nodiscard]] inline expression::Expr add(
    std::vector<expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return call(evaluation::BuiltinId::Add, std::move(arguments), registry, mathematics, angles);
}

[[nodiscard]] inline expression::Expr multiply(
    std::vector<expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return call(evaluation::BuiltinId::Multiply, std::move(arguments), registry, mathematics, angles);
}

[[nodiscard]] inline expression::Expr subtract(
    expression::Expr lhs,
    expression::Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return call(evaluation::BuiltinId::Subtract, {std::move(lhs), std::move(rhs)},
        registry, mathematics, angles);
}

[[nodiscard]] inline expression::Expr divide(
    expression::Expr lhs,
    expression::Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return call(evaluation::BuiltinId::Divide, {std::move(lhs), std::move(rhs)},
        registry, mathematics, angles);
}

[[nodiscard]] inline expression::Expr negate(
    expression::Expr value,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return call(evaluation::BuiltinId::Negate, {std::move(value)},
        registry, mathematics, angles);
}

[[nodiscard]] inline expression::Expr sqrt(
    expression::Expr value,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return call(evaluation::BuiltinId::Sqrt, {std::move(value)},
        registry, mathematics, angles);
}

} // namespace mmcal::builtins::exact
