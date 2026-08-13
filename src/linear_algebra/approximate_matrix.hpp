#pragma once

#include "approximation/approximation_context.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <optional>
#include <span>

namespace mmcal::linear_algebra {

[[nodiscard]] std::optional<expression::Expr> approximateDot(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateDeterminant(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateInverse(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateRref(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateMatrixRank(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateSolveLinear(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateNullSpace(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateNorm(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateNormalize(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateTrace(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

} // namespace mmcal::linear_algebra
