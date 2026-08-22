#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"

namespace mmcal::symbolic {

enum class LimitDirection {
    TwoSided,
    Left,
    Right
};

// 実軸上の記号極限。証明できない場合は limit[...] を保持する。
// Infinity は言語上のsentinel Symbolなので、MathRegistryの数学定数とは分離して受け取る。
[[nodiscard]] expression::Expr limitExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const expression::Expr& point,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const expression::Symbol* complexInfinitySymbol = nullptr,
    const expression::Symbol* indeterminateSymbol = nullptr);

} // namespace mmcal::symbolic
