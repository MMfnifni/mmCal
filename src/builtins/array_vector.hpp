#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluateIdentity(std::span<const expression::Expr> arguments);
[[nodiscard]] expression::Expr evaluateZeros(std::span<const expression::Expr> arguments);
[[nodiscard]] expression::Expr evaluateRows(std::span<const expression::Expr> arguments);
[[nodiscard]] expression::Expr evaluateCols(std::span<const expression::Expr> arguments);
[[nodiscard]] expression::Expr evaluateDiag(std::span<const expression::Expr> arguments);

[[nodiscard]] expression::Expr evaluateVectorAdd(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorSubtract(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorScale(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorCross(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorNorm(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorManhattan(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorEuclidean(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorNormalize(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorProject(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorAngle(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorReflect(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorReflectAxis(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateVectorSum(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::builtins
