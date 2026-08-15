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

// 実軸上の周期函数sin/cos/tanについて、affine引数に限り整数parameterを持つ全解族へ反転する。
[[nodiscard]] std::optional<SolutionSet> solveRealPeriodicFunctionRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// MathRegistryで「実軸上global injective」と証明済みの函数だけを
// principal inverseで反転する。
[[nodiscard]] std::optional<SolutionSet> solveRealInjectiveFunctionRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::solver
