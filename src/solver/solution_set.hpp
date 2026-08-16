#pragma once

#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/numeric_domain.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace mmcal::solver {

// Solverが何を変数として解いているかを、ambient domainと一緒に保持する。
// 同じ方程式でもReal上とComplex上では解集合が違うため、domainを結果から落とさない。
struct SolverVariable final {
    expression::Symbol symbol;
    mathematics::NumericDomain domain = mathematics::NumericDomain::Complex;

    [[nodiscard]] bool operator==(const SolverVariable&) const = default;
};

struct SolutionBinding final {
    expression::Symbol variable;
    expression::Expr value;

    [[nodiscard]] bool operator==(const SolutionBinding&) const = default;
};

// 1つの解の枝。
// conditionsは「この式になるための追加条件」を保持する。
// multiplicityは多項式solverなど、重根情報を証明できるsolverだけが設定する。
struct SolutionBranch final {
    std::vector<SolutionBinding> bindings;
    mathematics::AssumptionSet conditions;
    std::optional<std::size_t> multiplicity;
    // 解族をparameter化する自由変数。連立一次方程式の未拘束solver変数に加え、
    // 周期解の整数kのようなformal parameterも保持できる。
    std::vector<SolverVariable> freeVariables;
    // 専用solverが各bindingのdomain所属を数学的に証明済みなら、そのdomain以上への
    // 再検査を省ける。Real証明をRational/Integer証明へ誤って強めないためboolにはしない。
    std::optional<mathematics::NumericDomain> bindingsCertifiedDomain;

    [[nodiscard]] bool unconditional() const noexcept { return conditions.empty(); }
    [[nodiscard]] bool operator==(const SolutionBranch&) const = default;
};

enum class SolutionSetKind {
    Empty,
    Finite,
    Universal,
    Conditional,
    Unresolved
};

// 係数に未確定なSymbolを含む方程式では、解集合そのものが条件で分岐する。
// 例: a x + b == 0 は a!=0 / a==0,b==0 / a==0,b!=0 で
// Finite / Universal / Empty に分かれる。branchのside conditionだけでは
// Universal/Empty側を表せないため、集合レベルのcaseを独立して保持する。
struct SolutionCase final {
    mathematics::AssumptionSet conditions;
    SolutionSetKind outcome = SolutionSetKind::Unresolved;
    std::vector<SolutionBranch> branches;

    [[nodiscard]] bool operator==(const SolutionCase&) const = default;
};

// 「函数の値」と「方程式の全解」を混ぜないためのSolver専用結果型。
//
// Power[-8, 1/3] は将来principal Powerとして一価に定義する一方、
// solve[z^3 == -8, z] は3本のSolutionBranchを持つFinite集合になる。
class SolutionSet final {
public:
    [[nodiscard]] static SolutionSet empty(std::vector<SolverVariable> variables);
    [[nodiscard]] static SolutionSet finite(
        std::vector<SolverVariable> variables,
        std::vector<SolutionBranch> branches);
    [[nodiscard]] static SolutionSet universal(
        std::vector<SolverVariable> variables,
        mathematics::AssumptionSet conditions = {});
    [[nodiscard]] static SolutionSet conditional(
        std::vector<SolverVariable> variables,
        std::vector<SolutionCase> cases);
    [[nodiscard]] static SolutionSet unresolved(
        std::vector<SolverVariable> variables,
        mathematics::AssumptionSet conditions = {});

    [[nodiscard]] SolutionSetKind kind() const noexcept;
    [[nodiscard]] std::span<const SolverVariable> variables() const noexcept;
    [[nodiscard]] std::span<const SolutionBranch> branches() const noexcept;
    [[nodiscard]] std::span<const SolutionCase> cases() const noexcept;
    [[nodiscard]] const mathematics::AssumptionSet& conditions() const noexcept;
    [[nodiscard]] SolutionSet withAdditionalConditions(
        const mathematics::AssumptionSet& conditions) const;
    [[nodiscard]] bool operator==(const SolutionSet&) const = default;

private:
    SolutionSetKind kind_ = SolutionSetKind::Unresolved;
    std::vector<SolverVariable> variables_;
    std::vector<SolutionBranch> branches_;
    std::vector<SolutionCase> cases_;
    mathematics::AssumptionSet conditions_;

    SolutionSet(
        SolutionSetKind kind,
        std::vector<SolverVariable> variables,
        std::vector<SolutionBranch> branches,
        std::vector<SolutionCase> cases,
        mathematics::AssumptionSet conditions);
};

} // namespace mmcal::solver
