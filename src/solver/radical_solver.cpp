// 主値radical等式のexact solve
#include "radical_solver.hpp"

#include "mathematics/knowledge_context.hpp"
#include "numeric/number.hpp"
#include "polynomial_solver.hpp"
#include "solve_constraints.hpp"
#include "solver_support.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::solver {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::RelationKind;
using mathematics::TruthValue;
using numeric::Number;

enum class RadicalKind {
    Sqrt,
    Cbrt
};

struct RadicalSide final {
    RadicalKind kind = RadicalKind::Sqrt;
    Expr argument;
    Expr rhs;
};

enum class RangeDecision {
    Accept,
    Reject,
    RequireNonNegativeReal,
    Unknown
};

[[nodiscard]] bool isEqualRelation(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 2)
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    return definition && definition->id == BuiltinId::Equal;
}

[[nodiscard]] std::optional<RadicalSide> matchRadicalSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = builtins.find(lhs.asCall().head);
    if (!definition)
        return std::nullopt;

    RadicalKind kind;
    if (definition->id == BuiltinId::Sqrt)
        kind = RadicalKind::Sqrt;
    else if (definition->id == BuiltinId::Cbrt)
        kind = RadicalKind::Cbrt;
    else
        return std::nullopt;

    if (!symbolic::containsSymbol(lhs.asCall().arguments[0], variable)
        && !symbolic::containsSymbol(rhs, variable))
        return std::nullopt;
    return RadicalSide{kind, lhs.asCall().arguments[0], rhs};
}

void appendAssumptions(
    mathematics::AssumptionSet& target,
    const mathematics::AssumptionSet& source) {
    for (const auto& predicate : source.predicates())
        target.add(predicate);
}

[[nodiscard]] RangeDecision principalSqrtRange(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    // 主値sqrt自身の値は定義上，主値sqrtの像に入る。
    if (builtins.isCallTo(value, BuiltinId::Sqrt))
        return RangeDecision::Accept;

    if (value.isNumber()) {
        const Number& number = value.asNumber();
        if (number.isReal())
            return number.asReal().isNegative()
                ? RangeDecision::Reject
                : RangeDecision::Accept;

        // 主値sqrtの像は Re(w)>0 または Re(w)==0 && Im(w)>=0。
        const auto real = number.realPart();
        const auto imaginary = number.imaginaryPart();
        if (!real.isNegative() && !real.isZero())
            return RangeDecision::Accept;
        if (real.isNegative())
            return RangeDecision::Reject;
        return imaginary.isNegative()
            ? RangeDecision::Reject
            : RangeDecision::Accept;
    }

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto realPredicate = mathematics::elementOf(value, mathematics::NumericDomain::Real);
    const TruthValue real = knowledge.prove(realPredicate);
    if (real != TruthValue::True)
        return RangeDecision::Unknown;

    const auto nonnegative = mathematics::relation(
        RelationKind::GreaterEqual, value, integerExpr(0));
    const TruthValue sign = knowledge.prove(nonnegative);
    if (sign == TruthValue::True)
        return RangeDecision::Accept;
    if (sign == TruthValue::False)
        return RangeDecision::Reject;
    return RangeDecision::RequireNonNegativeReal;
}


[[nodiscard]] bool containsRootCall(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (builtins.isCallTo(expression, BuiltinId::Root))
        return true;
    if (expression.isCall())
        for (const Expr& argument : expression.asCall().arguments)
            if (containsRootCall(argument, builtins))
                return true;
    if (expression.isArray())
        for (std::size_t i = 0; i < expression.asArray().size(); ++i)
            if (containsRootCall(expression.asArray().element(i), builtins))
                return true;
    if (expression.isList())
        for (const Expr& element : expression.asList().elements)
            if (containsRootCall(element, builtins))
                return true;
    return false;
}

[[nodiscard]] bool containsRootBinding(
    const SolutionSet& solutions,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (solutions.kind() != SolutionSetKind::Finite)
        return false;
    for (const SolutionBranch& branch : solutions.branches())
        for (const SolutionBinding& binding : branch.bindings)
            if (binding.variable == variable && containsRootCall(binding.value, builtins))
                return true;
    return false;
}
[[nodiscard]] bool affineRealRhsForcesVariableReal(
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    const auto polynomial = symbolic::toRationalPolynomial(rhs, variable, builtins);
    return polynomial && polynomial->degree() == 1
        && !polynomial->coefficient(1).isZero();
}

[[nodiscard]] SolutionSet withAmbientComplexDomain(
    SolutionSet solutions,
    const expression::Symbol& variable) {
    const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Complex}};
    switch (solutions.kind()) {
    case SolutionSetKind::Empty:
        return SolutionSet::empty(variables).withAdditionalConditions(solutions.conditions());
    case SolutionSetKind::Finite: {
        std::vector<SolutionBranch> branches(solutions.branches().begin(), solutions.branches().end());
        for (SolutionBranch& branch : branches)
            branch.bindingsCertifiedDomain = mathematics::NumericDomain::Real;
        return SolutionSet::finite(variables, std::move(branches))
            .withAdditionalConditions(solutions.conditions());
    }
    default:
        return SolutionSet::unresolved(variables).withAdditionalConditions(solutions.conditions());
    }
}

[[nodiscard]] SolutionSet solveTransformedRelation(
    const RadicalSide& matched,
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    // cbrt[A]==a*x+b (a,b in Q, a!=0) なら、cbrtの値がRealであることから
    // 任意の解でxもRealと従う。高次多項式を最初からReal isolationへ渡せる。
    if (matched.kind == RadicalKind::Cbrt
        && affineRealRhsForcesVariableReal(matched.rhs, variable, builtins)) {
        SolutionSet ordinary = solveUnivariatePolynomialRelation(
            relation, variable, builtins, mathematics, angles);

        // Rational deflation / quadratic radicalsで閉じる場合は、そのcanonical outputを
        // 先に活かす。Complex Root fallbackまで到達した高次式だけReal isolationへ切り替える。
        if (!containsRootBinding(ordinary, variable, builtins)) {
            SolveConstraints realConstraint;
            realConstraint.domain = mathematics::NumericDomain::Real;
            ordinary = applySolveConstraints(
                std::move(ordinary), realConstraint, builtins, mathematics, angles);
            return withAmbientComplexDomain(std::move(ordinary), variable);
        }

        if (auto real = solveRealAlgebraicPolynomialEquation(
                relation, variable, builtins, mathematics, angles))
            return withAmbientComplexDomain(std::move(*real), variable);
        return ordinary;
    }

    return solveUnivariatePolynomialRelation(
        relation, variable, builtins, mathematics, angles);
}

[[nodiscard]] std::optional<SolutionSet> filterCandidates(
    SolutionSet candidates,
    const RadicalSide& matched,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Complex}};
    if (candidates.kind() == SolutionSetKind::Empty)
        return SolutionSet::empty(variables).withAdditionalConditions(candidates.conditions());
    if (candidates.kind() != SolutionSetKind::Finite)
        return std::nullopt;

    std::vector<SolutionBranch> accepted;
    accepted.reserve(candidates.branches().size());
    for (const SolutionBranch& candidate : candidates.branches()) {
        const auto binding = std::find_if(
            candidate.bindings.begin(), candidate.bindings.end(),
            [&](const SolutionBinding& item) { return item.variable == variable; });
        if (binding == candidate.bindings.end() || !candidate.freeVariables.empty())
            return std::nullopt;

        mathematics::AssumptionSet localAssumptions = assumptions;
        appendAssumptions(localAssumptions, candidate.conditions);
        Expr rhs = symbolic::substituteSymbol(matched.rhs, variable, binding->value);
        rhs = simplifyForSolve(std::move(rhs), builtins, mathematics, angles, localAssumptions);

        SolutionBranch branch = candidate;
        branch.multiplicity.reset();

        if (matched.kind == RadicalKind::Sqrt) {
            switch (principalSqrtRange(rhs, builtins, mathematics, localAssumptions)) {
            case RangeDecision::Reject:
                continue;
            case RangeDecision::Unknown:
                // 複素principal sqrtの像は半平面境界を含むdisjunctionになる。
                // 現AssumptionSetでそのORを弱めて表現すると解を落とすため、推測しない。
                return std::nullopt;
            case RangeDecision::RequireNonNegativeReal:
                branch.conditions.add(mathematics::relation(
                    RelationKind::GreaterEqual, rhs, integerExpr(0)));
                break;
            case RangeDecision::Accept:
                break;
            }
        }
        else {
            const mathematics::KnowledgeContext knowledge{
                builtins, mathematics, localAssumptions};
            const auto realPredicate = mathematics::elementOf(
                rhs, mathematics::NumericDomain::Real);
            const TruthValue real = knowledge.prove(realPredicate);
            if (real == TruthValue::False)
                continue;
            if (real == TruthValue::Unknown)
                branch.conditions.add(realPredicate);
        }
        accepted.push_back(std::move(branch));
    }

    if (accepted.empty())
        return SolutionSet::empty(variables).withAdditionalConditions(candidates.conditions());
    return SolutionSet::finite(variables, std::move(accepted))
        .withAdditionalConditions(candidates.conditions());
}

} // namespace

std::optional<SolutionSet> solveRadicalRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isEqualRelation(relation, builtins))
        return std::nullopt;

    const auto& sides = relation.asCall().arguments;
    auto matched = matchRadicalSide(sides[0], sides[1], variable, builtins);
    if (!matched)
        matched = matchRadicalSide(sides[1], sides[0], variable, builtins);
    if (!matched)
        return std::nullopt;

    const std::int64_t exponent = matched->kind == RadicalKind::Sqrt ? 2 : 3;
    Expr rhsPower = Expr::call(
        builtins.symbol(BuiltinId::Power), {matched->rhs, integerExpr(exponent)});
    Expr transformed = Expr::call(
        builtins.symbol(BuiltinId::Equal), {matched->argument, std::move(rhsPower)});
    transformed = simplifyForSolve(
        std::move(transformed), builtins, mathematics, angles, assumptions);

    SolutionSet candidates = solveTransformedRelation(
        *matched, transformed, variable, builtins, mathematics, angles);
    return filterCandidates(
        std::move(candidates), *matched, variable,
        builtins, mathematics, angles, assumptions);
}

} // namespace mmcal::solver
