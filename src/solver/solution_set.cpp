// 解集合表現
#include "solution_set.hpp"

#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace mmcal::solver {
namespace {

void validateVariables(std::span<const SolverVariable> variables) {
    std::unordered_set<symbols::SymbolId, symbols::SymbolIdHash> seen;
    for (const SolverVariable& variable : variables) {
        if (!variable.symbol.valid())
            throw std::invalid_argument("Solver variable is invalid");
        if (variable.domain == mathematics::NumericDomain::Unknown)
            throw std::invalid_argument("Solver variable domain cannot be Unknown");
        if (!seen.insert(variable.symbol.id()).second)
            throw std::invalid_argument("Solver variable is duplicated");
    }
}

void validateBranch(
    const SolutionBranch& branch,
    std::span<const SolverVariable> variables) {
    if (branch.multiplicity && *branch.multiplicity == 0)
        throw std::invalid_argument("Solution multiplicity must be positive");

    std::unordered_set<symbols::SymbolId, symbols::SymbolIdHash> allowed;
    for (const SolverVariable& variable : variables)
        allowed.insert(variable.symbol.id());

    std::unordered_set<symbols::SymbolId, symbols::SymbolIdHash> bound;
    for (const SolutionBinding& binding : branch.bindings) {
        if (!allowed.contains(binding.variable.id()))
            throw std::invalid_argument("Solution binds a symbol that is not a solver variable");
        if (!bound.insert(binding.variable.id()).second)
            throw std::invalid_argument("Solution binds the same variable more than once");
    }

    std::unordered_set<symbols::SymbolId, symbols::SymbolIdHash> free;
    for (const SolverVariable& variable : branch.freeVariables) {
        if (!allowed.contains(variable.symbol.id()))
            throw std::invalid_argument("Free solution parameter is not a solver variable");
        if (bound.contains(variable.symbol.id()))
            throw std::invalid_argument("A solver variable cannot be both bound and free");
        if (!free.insert(variable.symbol.id()).second)
            throw std::invalid_argument("Free solution parameter is duplicated");
    }

}

void validateCase(
    const SolutionCase& solutionCase,
    std::span<const SolverVariable> variables) {
    if (solutionCase.outcome == SolutionSetKind::Conditional)
        throw std::invalid_argument("A solution case cannot contain another conditional solution set");
    if (solutionCase.outcome == SolutionSetKind::Finite && solutionCase.branches.empty())
        throw std::invalid_argument("A finite solution case must contain at least one branch");
    if (solutionCase.outcome != SolutionSetKind::Finite && !solutionCase.branches.empty())
        throw std::invalid_argument("Only a finite solution case may contain explicit branches");
    for (const SolutionBranch& branch : solutionCase.branches)
        validateBranch(branch, variables);
}

} // namespace

SolutionSet::SolutionSet(
    SolutionSetKind kind,
    std::vector<SolverVariable> variables,
    std::vector<SolutionBranch> branches,
    std::vector<SolutionCase> cases,
    mathematics::AssumptionSet conditions)
    : kind_(kind),
      variables_(std::move(variables)),
      branches_(std::move(branches)),
      cases_(std::move(cases)),
      conditions_(std::move(conditions)) {
    validateVariables(variables_);

    if (kind_ != SolutionSetKind::Finite && !branches_.empty())
        throw std::invalid_argument("Only a finite solution set may contain explicit branches");
    if (kind_ == SolutionSetKind::Finite && branches_.empty())
        throw std::invalid_argument("Use SolutionSet::empty for an empty finite solution set");
    if (kind_ != SolutionSetKind::Conditional && !cases_.empty())
        throw std::invalid_argument("Only a conditional solution set may contain cases");
    if (kind_ == SolutionSetKind::Conditional && cases_.empty())
        throw std::invalid_argument("A conditional solution set must contain at least one case");

    for (const SolutionBranch& branch : branches_)
        validateBranch(branch, variables_);
    for (const SolutionCase& solutionCase : cases_)
        validateCase(solutionCase, variables_);
}

SolutionSet SolutionSet::empty(std::vector<SolverVariable> variables) {
    return SolutionSet{SolutionSetKind::Empty, std::move(variables), {}, {}, {}};
}

SolutionSet SolutionSet::finite(
    std::vector<SolverVariable> variables,
    std::vector<SolutionBranch> branches) {
    return SolutionSet{
        SolutionSetKind::Finite,
        std::move(variables),
        std::move(branches),
        {},
        {}};
}

SolutionSet SolutionSet::universal(
    std::vector<SolverVariable> variables,
    mathematics::AssumptionSet conditions) {
    return SolutionSet{
        SolutionSetKind::Universal,
        std::move(variables),
        {},
        {},
        std::move(conditions)};
}

SolutionSet SolutionSet::conditional(
    std::vector<SolverVariable> variables,
    std::vector<SolutionCase> cases) {
    return SolutionSet{
        SolutionSetKind::Conditional,
        std::move(variables),
        {},
        std::move(cases),
        {}};
}

SolutionSet SolutionSet::unresolved(
    std::vector<SolverVariable> variables,
    mathematics::AssumptionSet conditions) {
    return SolutionSet{
        SolutionSetKind::Unresolved, std::move(variables), {}, {}, std::move(conditions)};
}

SolutionSetKind SolutionSet::kind() const noexcept {
    return kind_;
}

std::span<const SolverVariable> SolutionSet::variables() const noexcept {
    return variables_;
}

std::span<const SolutionBranch> SolutionSet::branches() const noexcept {
    return branches_;
}

std::span<const SolutionCase> SolutionSet::cases() const noexcept {
    return cases_;
}

const mathematics::AssumptionSet& SolutionSet::conditions() const noexcept {
    return conditions_;
}

SolutionSet SolutionSet::withAdditionalConditions(
    const mathematics::AssumptionSet& conditions) const {
    SolutionSet result = *this;
    for (const mathematics::Predicate& predicate : conditions.predicates())
        result.conditions_.add(predicate);
    return result;
}

} // namespace mmcal::solver
