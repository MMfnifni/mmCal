#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <span>

namespace mmcal::builtins {

// diff[expr, x, at] / diff[expr, x, at, digits]
// 記号微分Dの結果をcertified evaluatorで点評価する数値微分。
[[nodiscard]] expression::Expr evaluateNumericDerivative(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// nintegrate[expr, {x, a, b}] / nintegrate[expr, {x, a, b}, digits]
// certifiableな有限実区間をt∈[0,1]へ写像し、9階導函数のinterval boundとexact Rationalな9点closed Newton-Cotes誤差包絡を用いたcertified enclosureを返す。
[[nodiscard]] expression::Expr evaluateNumericIntegral(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::builtins
