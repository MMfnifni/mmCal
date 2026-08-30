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

// 一変数Real等式について，明示的な逆函数を持たなくても
// 符号・函数値域・単調性・凸性から完全性を証明できる場合だけ解集合を返す。
// 数値scanや「根が見つからなかった」ことはNoSolutionの根拠にしない。
[[nodiscard]] std::optional<SolutionSet> solveRealEquationByProof(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::solver
