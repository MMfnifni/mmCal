#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/math_registry.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluateAbs(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);
[[nodiscard]] expression::Expr evaluateSign(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);
[[nodiscard]] expression::Expr evaluateRe(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);
[[nodiscard]] expression::Expr evaluateIm(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);
[[nodiscard]] expression::Expr evaluateConj(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);

} // namespace mmcal::builtins
