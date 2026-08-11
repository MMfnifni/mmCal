#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/math_registry.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluateSinh(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&);
[[nodiscard]] expression::Expr evaluateCosh(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&);
[[nodiscard]] expression::Expr evaluateTanh(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&);
[[nodiscard]] expression::Expr evaluateAsinh(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&);
[[nodiscard]] expression::Expr evaluateAcosh(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&);
[[nodiscard]] expression::Expr evaluateAtanh(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&);
[[nodiscard]] expression::Expr evaluateCsch(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&);
[[nodiscard]] expression::Expr evaluateSech(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&);
[[nodiscard]] expression::Expr evaluateCoth(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&);

} // namespace mmcal::builtins
