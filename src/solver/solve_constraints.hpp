#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "mathematics/numeric_domain.hpp"
#include "solution_set.hpp"

#include <optional>
#include <span>

namespace mmcal::solver {

// solveの第3引数を、探索domainと追加条件へ正規化したもの。
// domainは全solver変数へ適用し、assumptionsは x!=0 等の局所条件を保持する。
struct SolveConstraints final {
    // 未指定と「Complexを明示指定」を区別する。等式solverの既定はComplexだが、
    // 順序不等式はそれ自体がReal domainを要求するため、solverが決めたambient domainを
    // 第3引数なしで上書きしてはいけない。
    std::optional<mathematics::NumericDomain> domain;
    mathematics::AssumptionSet assumptions;
};

// spec は Real / Complex / Integer / Rational、比較式、またはそれらの一次元配列を受ける。
[[nodiscard]] SolveConstraints parseSolveConstraints(
    const expression::Expr& spec,
    std::span<const expression::Symbol> variables,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// 解候補をdomain/条件で絞り込む。証明できない条件は捨てずにSolutionBranchへ残す。
[[nodiscard]] SolutionSet applySolveConstraints(
    SolutionSet solutions,
    const SolveConstraints& constraints,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::solver
