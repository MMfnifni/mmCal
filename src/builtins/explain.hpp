#pragma once

#include "expression/expr.hpp"
#include "evaluation/builtin_registry.hpp"
#include "mathematics/math_registry.hpp"
#include "symbols/symbol_registry.hpp"

#include <span>

namespace mmcal::builtins {

// 評価済みExprが既に保持している安価なmetadataだけを返す。
// 数学的な追加計算や全要素走査は行わない。
[[nodiscard]] expression::Expr evaluateExplain(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const symbols::SymbolRegistry& symbols,
    const mathematics::MathRegistry& mathematics);

} // namespace mmcal::builtins
