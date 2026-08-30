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

// 主値sqrt / 実cbrtの等式を，冪乗後の多項式候補と値域条件の組としてexactに解く。
// sqrt[A]==B は A==B^2 だけでは不十分なので、Bがprincipal sqrtの像に入ることまで証明する。
// cbrt[A]==B は A==B^3 と B in Real が同値条件になる。
[[nodiscard]] std::optional<SolutionSet> solveRadicalRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::solver
