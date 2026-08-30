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

// Real variableに対するabs等式・不等式をexactに変換する。
// equalityはComplex上では円/曲線になり得るため，explicitRealDomain=falseでは扱わない。
[[nodiscard]] std::optional<SolutionSet> solveRealAbsoluteValueRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    bool explicitRealDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::solver
