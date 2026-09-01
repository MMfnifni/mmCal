#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <cstdint>
#include <optional>

namespace mmcal::symbolic {

// principal branch と現在の角度意味論を保った記号微分。
// 安全な微分則を持たない非正則函数は D[...] のまま保持する。
[[nodiscard]] expression::Expr differentiateExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// D[..., {x,n}]の各段で生じる入れ子の加減算を，exactな有理係数の
// 線形結合としてだけ平坦化する。一般expandは行わず，積やbranch構造は保つ。
[[nodiscard]] expression::Expr canonicalizeDerivativeOutput(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// 高階微分を閉形式で構成できる既知family用のfast path。
// 現在はdirect-variable Polylogを扱い，非該当ならnulloptを返す。
[[nodiscard]] std::optional<expression::Expr> differentiateKnownRepeatedExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    std::uint64_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::symbolic
