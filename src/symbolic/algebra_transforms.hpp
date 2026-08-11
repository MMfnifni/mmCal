#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <cstddef>
#include <span>

namespace mmcal::symbolic {

struct AlgebraTransformOptions final {
    // expandの項数爆発を防ぐ。上限を超える場合は元の部分式を保持する。
    std::size_t maximumExpandedTerms = 4096;
};

[[nodiscard]] expression::Expr expandExpression(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    AlgebraTransformOptions options = {});

[[nodiscard]] expression::Expr collectExpression(
    const expression::Expr& expression,
    std::span<const expression::Symbol> variables,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] expression::Expr collectExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] expression::Expr factorExpression(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::symbolic
