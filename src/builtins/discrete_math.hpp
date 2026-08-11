#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluateFloor(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateCeil(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateTrunc(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateRound(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateFrac(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateGcd(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateLcm(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateMod(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateRem(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateQuotient(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateNextPow2(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::builtins
