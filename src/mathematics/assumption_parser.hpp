#pragma once

#include "assumption_set.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "math_registry.hpp"

namespace mmcal::mathematics {

// ユーザーが与えた比較式・domain所属を内部Predicateへ変換する。
// Simplifyと将来のRefine/Assumingで同じ解釈を共有するため、Evaluatorへ埋め込まない。
[[nodiscard]] AssumptionSet parseAssumptions(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics);

} // namespace mmcal::mathematics
