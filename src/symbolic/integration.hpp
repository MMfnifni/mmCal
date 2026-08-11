#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"

namespace mmcal::symbolic {

// 不定積分の原始函数代表元を返す。積分定数は表示しない。
// assumptionsは積分前の安全な簡約・definedness証明に利用する。
// 未対応の場合は integrate[expression, variable] をそのまま返す。
[[nodiscard]] expression::Expr integrateExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// 実区間の記号定積分。Infinity endpointと証明可能なimproper integralも扱う。
// 安全性を証明できない場合は integrate[...] を保持する。
[[nodiscard]] expression::Expr integrateExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const expression::Expr& lower,
    const expression::Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::symbolic
