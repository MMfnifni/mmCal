#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluateGroebnerBasis(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins);

[[nodiscard]] expression::Expr evaluatePolynomialReduce(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins);

} // namespace mmcal::builtins
