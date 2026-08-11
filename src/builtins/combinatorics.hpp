#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluatePermutation(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateCombination(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateFibonacci(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);

} // namespace mmcal::builtins
