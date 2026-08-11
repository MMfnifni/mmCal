#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

namespace mmcal::symbolic {

// principal branch と現在の角度意味論を保った記号微分。
// 安全な微分則を持たない非正則函数は D[...] のまま保持する。
[[nodiscard]] expression::Expr differentiateExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::symbolic
