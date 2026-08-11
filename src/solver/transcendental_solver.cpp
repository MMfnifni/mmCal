// 実軸で安全な超越方程式反転
#include "transcendental_solver.hpp"

#include "mathematics/knowledge_context.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "polynomial_solver.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "symbolic/polynomial.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace mmcal::solver {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::FunctionDefinition;
using mathematics::RealRangeRule;
using mathematics::RelationKind;
using mathematics::TruthValue;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] bool isEqualRelation(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 2)
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    return definition && definition->id == BuiltinId::Equal;
}

[[nodiscard]] bool containsVariable(
    const Expr& expression,
    const expression::Symbol& variable) {
    return symbolic::containsSymbol(expression, variable);
}

struct FunctionSide final {
    const FunctionDefinition* definition = nullptr;
    Expr argument;
    Expr rhs;
};

[[nodiscard]] std::optional<FunctionSide> matchFunctionSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const mathematics::MathRegistry& mathematics) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 1
        || !containsVariable(lhs.asCall().arguments[0], variable)
        || containsVariable(rhs, variable))
        return std::nullopt;
    const auto* definition = mathematics.findFunction(lhs.asCall().head);
    if (!definition || !definition->inverseFunction || !definition->realGloballyInjective)
        return std::nullopt;
    return FunctionSide{definition, lhs.asCall().arguments[0], rhs};
}

[[nodiscard]] mathematics::AssumptionSet rangeConditions(
    const FunctionDefinition& definition,
    const Expr& rhs) {
    mathematics::AssumptionSet result;
    result.add(mathematics::elementOf(rhs, mathematics::NumericDomain::Real));
    switch (definition.realRangeRule) {
    case RealRangeRule::Positive:
        result.add(mathematics::relation(RelationKind::Greater, rhs, integer(0)));
        break;
    case RealRangeRule::NonNegative:
        result.add(mathematics::relation(RelationKind::GreaterEqual, rhs, integer(0)));
        break;
    case RealRangeRule::OpenMinusOneToOne:
        result.add(mathematics::relation(RelationKind::Greater, rhs, integer(-1)));
        result.add(mathematics::relation(RelationKind::Less, rhs, integer(1)));
        break;
    case RealRangeRule::ClosedMinusOneToOne:
        result.add(mathematics::relation(RelationKind::GreaterEqual, rhs, integer(-1)));
        result.add(mathematics::relation(RelationKind::LessEqual, rhs, integer(1)));
        break;
    case RealRangeRule::OneToInfinity:
        result.add(mathematics::relation(RelationKind::GreaterEqual, rhs, integer(1)));
        break;
    case RealRangeRule::AllReal:
    case RealRangeRule::Unknown:
        break;
    }
    return result;
}

[[nodiscard]] std::optional<mathematics::AssumptionSet> proveRangeConditions(
    const FunctionDefinition& definition,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    mathematics::AssumptionSet conditions = rangeConditions(definition, rhs);
    mathematics::AssumptionSet remaining;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    for (const auto& predicate : conditions.predicates()) {
        const TruthValue truth = knowledge.prove(predicate);
        if (truth == TruthValue::False)
            return std::nullopt;
        if (truth == TruthValue::Unknown)
            remaining.add(predicate);
    }
    return remaining;
}


[[nodiscard]] mathematics::Predicate simplifyPredicate(
    const mathematics::Predicate& predicate,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const simplification::SimplificationContext context{
        builtins, mathematics, angles, assumptions};
    if (const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate))
        return mathematics::relation(
            relation->relation,
            simplification::Simplifier{}.simplify(relation->lhs, context),
            simplification::Simplifier{}.simplify(relation->rhs, context));
    const auto& domain = std::get<mathematics::DomainPredicate>(predicate);
    return mathematics::elementOf(
        simplification::Simplifier{}.simplify(domain.expression, context), domain.domain);
}

} // namespace

std::optional<SolutionSet> solveRealInjectiveFunctionRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isEqualRelation(relation, builtins))
        return std::nullopt;

    const auto& sides = relation.asCall().arguments;
    auto matched = matchFunctionSide(sides[0], sides[1], variable, mathematics);
    if (!matched)
        matched = matchFunctionSide(sides[1], sides[0], variable, mathematics);
    if (!matched)
        return std::nullopt;

    matched->rhs = simplification::Simplifier{}.simplify(
        matched->rhs,
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});

    const mathematics::AssumptionSet allRangeConditions =
        rangeConditions(*matched->definition, matched->rhs);
    const auto remainingConditions = proveRangeConditions(
        *matched->definition, matched->rhs, builtins, mathematics, assumptions);
    if (!remainingConditions)
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    const auto* inverse = mathematics.findFunction(*matched->definition->inverseFunction);
    if (!inverse)
        return std::nullopt;

    Expr inverseValue = Expr::call(inverse->symbol, {matched->rhs});
    Expr transformed = Expr::call(
        builtins.symbol(BuiltinId::Equal), {matched->argument, std::move(inverseValue)});
    SolutionSet result = solveUnivariatePolynomialRelation(
        transformed, variable, builtins, mathematics, angles);

    // 逆函数の値域条件は、変換後に現れたLog/Atanh等のdefinednessを含意する。
    // polynomial solverが安全のため収集した条件を、この追加知識で再検証して a>0 と a!=0 のような冗長条件を残さない。
    mathematics::AssumptionSet proofAssumptions = assumptions;
    for (const auto& predicate : allRangeConditions.predicates())
        proofAssumptions.add(predicate);
    const mathematics::KnowledgeContext proofKnowledge{builtins, mathematics, proofAssumptions};
    mathematics::AssumptionSet filteredGlobal;
    for (const auto& predicate : result.conditions().predicates()) {
        const mathematics::Predicate simplifiedPredicate = simplifyPredicate(
            predicate, builtins, mathematics, angles, proofAssumptions);
        const TruthValue truth = proofKnowledge.prove(simplifiedPredicate);
        if (truth == TruthValue::False)
            return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});
        if (truth == TruthValue::Unknown)
            filteredGlobal.add(simplifiedPredicate);
    }

    if (result.kind() == SolutionSetKind::Finite) {
        std::vector<SolutionBranch> branches(result.branches().begin(), result.branches().end());
        for (SolutionBranch& branch : branches)
            for (const auto& predicate : remainingConditions->predicates())
                branch.conditions.add(predicate);
        SolutionSet cleaned = SolutionSet::finite(
            std::vector<SolverVariable>{result.variables().begin(), result.variables().end()},
            std::move(branches));
        return cleaned.withAdditionalConditions(filteredGlobal);
    }

    SolutionSet cleaned = result.withAdditionalConditions(filteredGlobal);
    return cleaned.withAdditionalConditions(*remainingConditions);
}

} // namespace mmcal::solver
