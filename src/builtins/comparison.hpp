#pragma once

#include "expression/expr.hpp"
#include "evaluation/builtin_registry.hpp"
#include "mathematics/math_registry.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluateComparison(
    const expression::Symbol& head,
    std::span<const expression::Expr> arguments);
[[nodiscard]] expression::Expr evaluateLogicalAnd(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateElement(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);

} // namespace mmcal::builtins
