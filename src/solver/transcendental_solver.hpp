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

// Lambert Wの明示branchについて，実値像をexactに証明できる範囲だけ
// W_k(x)==r を x==r Exp[r] へ反転する。W_0は[-1,inf)，W_-1は(-inf,-1]。
// 一般の複素targetに対するbranch imageは推測せず未解決を保つ。
[[nodiscard]] std::optional<SolutionSet> solveLambertBranchRelation(
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

// Complex上のbranch-aware transcendental inversion。Exp / Sin / Cos / Tanに加え，
// Sinh / Cosh / Tanhの虚周期解族，u Exp[u]==aのLambert W全整数branchを扱う。
// principal Log[u]==cはprincipal-log像(-Pi,Pi]をexact/certifiedに検査して反転する。
[[nodiscard]] std::optional<SolutionSet> solveComplexTranscendentalRelation(
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
