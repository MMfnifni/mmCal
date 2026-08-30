#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"
#include "solution_set.hpp"

#include <optional>

namespace mmcal::solver {

// 主値Lambert Wの値が実主値域[-1,inf)にあると証明できる場合だけ，
// W_0(x)==r をdomain非依存に x==r Exp[r] へ反転する。
// 複素r一般のprincipal-range判定はここでは推測しない。
[[nodiscard]] std::optional<SolutionSet> solvePrincipalLambertRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// 実軸上の周期函数sin/cos/tanについて、affine引数に限り整数parameterを持つ全解族へ反転する。
[[nodiscard]] std::optional<SolutionSet> solveRealPeriodicFunctionRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// 実指数方程式を分類する。正の定数baseの対数反転，
// affine exponential == affine，boundedなx^x constant target，a^x==x^2をLambert Wへ落とす。
[[nodiscard]] std::optional<SolutionSet> solveRealExponentialRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// MathRegistryで「実軸上global injective」と証明済みの函数だけを
// 主値逆函数で反転する。
[[nodiscard]] std::optional<SolutionSet> solveRealInjectiveFunctionRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::solver
