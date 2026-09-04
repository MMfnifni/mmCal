// Real絶対値関係式solver
#include "absolute_value_solver.hpp"

#include "mathematics/definedness.hpp"
#include "mathematics/knowledge_context.hpp"
#include "polynomial_solver.hpp"
#include "solve_constraints.hpp"
#include "solve_normalization.hpp"
#include "solver_support.hpp"
#include "symbolic/polynomial.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::solver {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::NumericDomain;
using mathematics::RelationKind;
using mathematics::TruthValue;

enum class RhsSign {
    Negative,
    Zero,
    Positive,
    Unknown
};

struct AbsoluteSide final {
    Expr argument;
    Expr rhs;
    RelationKind relation = RelationKind::Equal;
};

[[nodiscard]] std::optional<AbsoluteSide> matchAbsoluteSide(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    const auto kind = relationKindOf(relation, builtins);
    if (!kind)
        return std::nullopt;

    const auto& sides = relation.asCall().arguments;
    const auto match = [&](const Expr& candidate, const Expr& rhs, RelationKind relationKind)
        -> std::optional<AbsoluteSide> {
        if (!builtins.isCallTo(candidate, BuiltinId::Abs)
            || candidate.asCall().arguments.size() != 1)
            return std::nullopt;
        if (!symbolic::containsSymbol(candidate.asCall().arguments[0], variable)
            && !symbolic::containsSymbol(rhs, variable))
            return std::nullopt;
        return AbsoluteSide{candidate.asCall().arguments[0], rhs, relationKind};
    };

    if (auto lhs = match(sides[0], sides[1], *kind))
        return lhs;
    return match(sides[1], sides[0], reversedRelation(*kind));
}

[[nodiscard]] RhsSign classifyRhsSign(
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::relation(
            RelationKind::Less, rhs, integerExpr(0))) == TruthValue::True)
        return RhsSign::Negative;
    if (knowledge.prove(mathematics::relation(
            RelationKind::Equal, rhs, integerExpr(0))) == TruthValue::True)
        return RhsSign::Zero;
    if (knowledge.prove(mathematics::relation(
            RelationKind::Greater, rhs, integerExpr(0))) == TruthValue::True)
        return RhsSign::Positive;
    return RhsSign::Unknown;
}

[[nodiscard]] SolutionSet solveRealTransformed(
    Expr relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    relation = normalizeForSolve(relation, builtins, mathematics, angles);
    SolveConstraints real;
    real.domain = NumericDomain::Real;
    SolutionSet result = applySolveConstraints(
        solveUnivariatePolynomialRelation(
            relation, variable, builtins, mathematics, angles),
        real, builtins, mathematics, angles);
    if (result.kind() == SolutionSetKind::Unresolved) {
        if (auto algebraic = solveRealAlgebraicPolynomialEquation(
                relation, variable, builtins, mathematics, angles))
            result = *algebraic;
    }
    return result;
}

void clearMultiplicities(SolutionSet& solutions) {
    if (solutions.kind() != SolutionSetKind::Finite)
        return;
    std::vector<SolutionBranch> branches(solutions.branches().begin(), solutions.branches().end());
    for (SolutionBranch& branch : branches)
        branch.multiplicity.reset();
    solutions = SolutionSet::finite(
        std::vector<SolverVariable>{solutions.variables().begin(), solutions.variables().end()},
        std::move(branches)).withAdditionalConditions(solutions.conditions());
}

[[nodiscard]] std::optional<SolutionSet> unionFinite(
    SolutionSet lhs,
    SolutionSet rhs,
    const expression::Symbol& variable) {
    const std::vector<SolverVariable> variables{{variable, NumericDomain::Real}};
    if (lhs.kind() == SolutionSetKind::Universal || rhs.kind() == SolutionSetKind::Universal)
        return SolutionSet::universal(variables);
    if (lhs.kind() == SolutionSetKind::Empty)
        return rhs;
    if (rhs.kind() == SolutionSetKind::Empty)
        return lhs;
    if (lhs.kind() != SolutionSetKind::Finite || rhs.kind() != SolutionSetKind::Finite
        || !lhs.conditions().empty() || !rhs.conditions().empty())
        return std::nullopt;

    std::vector<SolutionBranch> branches(lhs.branches().begin(), lhs.branches().end());
    for (const SolutionBranch& branch : rhs.branches())
        if (std::find(branches.begin(), branches.end(), branch) == branches.end())
            branches.push_back(branch);
    if (branches.empty())
        return SolutionSet::empty(variables);
    return SolutionSet::finite(variables, std::move(branches));
}


[[nodiscard]] std::optional<SolutionSet> universalOnDefinedDomain(
    const Expr& argument,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto argumentDomain = mathematics::expressionDomainConditions(
        argument, builtins, mathematics);
    const auto rhsDomain = mathematics::expressionDomainConditions(
        rhs, builtins, mathematics);
    if (!argumentDomain || !rhsDomain)
        return std::nullopt;

    SolveConstraints domain;
    domain.domain = NumericDomain::Real;
    for (const auto& predicate : argumentDomain->predicates())
        domain.assumptions.add(predicate);
    for (const auto& predicate : rhsDomain->predicates())
        domain.assumptions.add(predicate);

    const std::vector<SolverVariable> variables{{variable, NumericDomain::Real}};
    return applySolveConstraints(
        SolutionSet::universal(variables), domain, builtins, mathematics, angles);
}

} // namespace

std::optional<SolutionSet> solveRealAbsoluteValueRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    bool explicitRealDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    auto matched = matchAbsoluteSide(relation, variable, builtins);
    if (!matched)
        return std::nullopt;

    const bool ordered = matched->relation == RelationKind::Less
        || matched->relation == RelationKind::LessEqual
        || matched->relation == RelationKind::Greater
        || matched->relation == RelationKind::GreaterEqual;
    if (!ordered && !explicitRealDomain)
        return std::nullopt;

    const mathematics::AssumptionSet local = withRealVariable(assumptions, variable);
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, local};
    const bool argumentReal = knowledge.prove(mathematics::elementOf(
        matched->argument, NumericDomain::Real)) == TruthValue::True;
    if (knowledge.prove(mathematics::elementOf(
            matched->rhs, NumericDomain::Real)) != TruthValue::True)
        return std::nullopt;

    matched->argument = simplifyForSolve(
        std::move(matched->argument), builtins, mathematics, angles, local);
    matched->rhs = simplifyForSolve(
        std::move(matched->rhs), builtins, mathematics, angles, local);
    const RhsSign sign = classifyRhsSign(
        matched->rhs, builtins, mathematics, local);
    const std::vector<SolverVariable> variables{{variable, NumericDomain::Real}};

    if (matched->relation == RelationKind::Equal && argumentReal) {
        // |u|=u / |u|=-u はrhsの符号を先に推測する必要がなく，
        // それぞれu>=0 / u<=0とexactに同値である。symbolic rhsを
        // Unknown扱いしてUnresolvedへ落とさず，既存のReal relation solverへ渡す。
        if (matched->rhs == matched->argument)
            return solveRealTransformed(
                relationExpr(RelationKind::GreaterEqual, matched->argument, integerExpr(0), builtins),
                variable, builtins, mathematics, angles);

        Expr negativeArgument = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Negate), {matched->argument}),
            builtins, mathematics, angles, local);
        if (matched->rhs == negativeArgument)
            return solveRealTransformed(
                relationExpr(RelationKind::LessEqual, matched->argument, integerExpr(0), builtins),
                variable, builtins, mathematics, angles);
    }

    if (matched->relation == RelationKind::Equal) {
        if (!argumentReal)
            return std::nullopt;
        if (sign == RhsSign::Negative)
            return SolutionSet::empty(variables);
        if (sign == RhsSign::Unknown)
            return std::nullopt;

        if (sign == RhsSign::Zero) {
            Expr transformed = Expr::call(
                builtins.symbol(BuiltinId::Equal), {matched->argument, integerExpr(0)});
            SolutionSet result = solveRealTransformed(
                std::move(transformed), variable, builtins, mathematics, angles);
            clearMultiplicities(result);
            return result;
        }

        Expr negativeRhs = simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Negate), {matched->rhs}),
            builtins, mathematics, angles, local);
        SolutionSet positive = solveRealTransformed(
            Expr::call(builtins.symbol(BuiltinId::Equal), {matched->argument, matched->rhs}),
            variable, builtins, mathematics, angles);
        SolutionSet negative = solveRealTransformed(
            Expr::call(builtins.symbol(BuiltinId::Equal), {matched->argument, std::move(negativeRhs)}),
            variable, builtins, mathematics, angles);
        clearMultiplicities(positive);
        clearMultiplicities(negative);
        return unionFinite(std::move(positive), std::move(negative), variable);
    }

    if (matched->relation == RelationKind::NotEqual) {
        if (sign == RhsSign::Negative)
            return universalOnDefinedDomain(
                matched->argument, matched->rhs, variable,
                builtins, mathematics, angles);
        if (sign == RhsSign::Unknown || !argumentReal)
            return std::nullopt;
    }

    switch (matched->relation) {
    case RelationKind::Less:
        if (sign == RhsSign::Negative || sign == RhsSign::Zero)
            return SolutionSet::empty(variables);
        if (sign == RhsSign::Unknown)
            return std::nullopt;
        break;
    case RelationKind::LessEqual:
        if (sign == RhsSign::Negative)
            return SolutionSet::empty(variables);
        if (sign == RhsSign::Unknown)
            return std::nullopt;
        break;
    case RelationKind::Greater:
        if (sign == RhsSign::Negative)
            return universalOnDefinedDomain(
                matched->argument, matched->rhs, variable,
                builtins, mathematics, angles);
        if (sign == RhsSign::Unknown)
            return std::nullopt;
        break;
    case RelationKind::GreaterEqual:
        if (sign == RhsSign::Negative || sign == RhsSign::Zero)
            return universalOnDefinedDomain(
                matched->argument, matched->rhs, variable,
                builtins, mathematics, angles);
        if (sign == RhsSign::Unknown)
            return std::nullopt;
        break;
    case RelationKind::Equal:
        return std::nullopt;
    case RelationKind::NotEqual:
        break;
    }

    if (!argumentReal)
        return std::nullopt;

    Expr lhsSquared = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Power), {matched->argument, integerExpr(2)}),
        builtins, mathematics, angles, local);
    Expr rhsSquared = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Power), {matched->rhs, integerExpr(2)}),
        builtins, mathematics, angles, local);
    Expr transformed = Expr::call(
        builtins.symbol(builtinForRelation(matched->relation)),
        {std::move(lhsSquared), std::move(rhsSquared)});
    return solveRealTransformed(
        std::move(transformed), variable, builtins, mathematics, angles);
}

} // namespace mmcal::solver
