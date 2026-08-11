#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/math_registry.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluateCbrt(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);
[[nodiscard]] expression::Expr evaluateHypot(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateCis(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluatePolar(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);

[[nodiscard]] expression::Expr evaluateDegreeToRadian(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);
[[nodiscard]] expression::Expr evaluateDegreeToGradian(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateRadianToDegree(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);
[[nodiscard]] expression::Expr evaluateRadianToGradian(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);
[[nodiscard]] expression::Expr evaluateGradianToDegree(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry);
[[nodiscard]] expression::Expr evaluateGradianToRadian(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics);

} // namespace mmcal::builtins
