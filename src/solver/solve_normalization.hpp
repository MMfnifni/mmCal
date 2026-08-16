#pragma once

#include "expression/expr.hpp"
#include "mathematics/assumption_set.hpp"

namespace mmcal::evaluation { class BuiltinRegistry; }
namespace mmcal::mathematics { class MathRegistry; class AngleSemantics; }

namespace mmcal::solver {

// SolveはHoldAllなので，一般Evaluatorを通すとuser definitionや副作用まで実行し得る。
// ここではbuiltin aliasのcanonical head化と，Simplifierが証明付きで行える数学的書換えだけを適用する。
[[nodiscard]] expression::Expr normalizeForSolve(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::solver
