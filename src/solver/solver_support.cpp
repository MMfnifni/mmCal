// solver間で共有する小さな副作用なし補助処理
#include "solver_support.hpp"

#include "evaluation/builtin_registry.hpp"
#include "mathematics/numeric_domain.hpp"
#include "mathematics/definedness.hpp"
#include "mathematics/relation_builtin.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"

#include <utility>

namespace mmcal::solver {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::RelationKind;
using numeric::BigInt;
using numeric::Number;

Expr integerExpr(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

Expr builtinCall(
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id,
    std::vector<Expr> arguments) {
    return Expr::call(builtins.symbol(id), std::move(arguments));
}

Expr simplifyForSolve(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{
            builtins, mathematics, angles, assumptions});
}

std::optional<RelationKind> relationKindOf(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return std::nullopt;

    return mathematics::relationKindForBuiltin(definition->id);
}

RelationKind reversedRelation(RelationKind relation) noexcept {
    return mathematics::reverseRelation(relation);
}

BuiltinId builtinForRelation(RelationKind relation) noexcept {
    return mathematics::builtinForRelation(relation);
}

Expr relationExpr(
    RelationKind relation,
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins) {
    return builtinCall(
        builtins, solver::builtinForRelation(relation), {std::move(lhs), std::move(rhs)});
}

std::optional<SolutionSet> solveIdenticalEquality(
    const Expr& relation,
    std::vector<SolverVariable> variables,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto kind = relationKindOf(relation, builtins);
    if (!kind || (*kind != RelationKind::Equal && *kind != RelationKind::NotEqual))
        return std::nullopt;

    const Expr& lhs = relation.asCall().arguments[0];
    const Expr& rhs = relation.asCall().arguments[1];
    const auto lhsDomain = mathematics::expressionDomainConditions(lhs, builtins, mathematics);
    const auto rhsDomain = mathematics::expressionDomainConditions(rhs, builtins, mathematics);
    if (!lhsDomain || !rhsDomain)
        return std::nullopt;

    mathematics::AssumptionSet definedness;
    for (const auto& predicate : lhsDomain->predicates())
        definedness.add(predicate);
    for (const auto& predicate : rhsDomain->predicates())
        definedness.add(predicate);

    // syntactically同一でなくても，両辺が定義される領域ではinverse composition等が
    // 安全に同一式へ落ちることがある。definednessとsolver変数domainだけを仮定し，
    // generic Simplifierが証明できた場合に限って恒等式として扱う。
    mathematics::AssumptionSet proofAssumptions = definedness;
    for (const SolverVariable& variable : variables) {
        if (mathematics::isSubdomainOf(variable.domain, mathematics::NumericDomain::Real))
            proofAssumptions.add(mathematics::elementOf(
                Expr{variable.symbol}, variable.domain));
    }

    const Expr normalizedLhs = simplifyForSolve(
        lhs, builtins, mathematics, angles, proofAssumptions);
    const Expr normalizedRhs = simplifyForSolve(
        rhs, builtins, mathematics, angles, proofAssumptions);
    if (normalizedLhs != normalizedRhs)
        return std::nullopt;

    if (*kind == RelationKind::NotEqual)
        return SolutionSet::empty(std::move(variables));
    return SolutionSet::universal(std::move(variables), definedness);
}

mathematics::AssumptionSet withRealVariable(
    const mathematics::AssumptionSet& assumptions,
    const expression::Symbol& variable) {
    mathematics::AssumptionSet result = assumptions;
    result.add(mathematics::elementOf(
        Expr{variable}, mathematics::NumericDomain::Real));
    return result;
}

std::size_t expressionNodeCount(const Expr& expression, std::size_t limit) {
    std::size_t count = 1;
    if (count > limit)
        return count;

    const auto addChild = [&](const Expr& child, std::size_t& total) {
        total += expressionNodeCount(child, limit > total ? limit - total : 0);
        return total <= limit;
    };

    if (expression.isCall()) {
        for (const Expr& argument : expression.asCall().arguments)
            if (!addChild(argument, count))
                return count;
    }
    else if (expression.isArray()) {
        for (std::size_t i = 0; i < expression.asArray().size(); ++i)
            if (!addChild(expression.asArray().element(i), count))
                return count;
    }
    else if (expression.isList()) {
        for (const Expr& element : expression.asList().elements)
            if (!addChild(element, count))
                return count;
    }
    return count;
}

bool containsBuiltinCall(
    const Expr& expression,
    BuiltinId id,
    const evaluation::BuiltinRegistry& builtins) {
    if (builtins.isCallTo(expression, id))
        return true;
    if (expression.isCall()) {
        for (const Expr& argument : expression.asCall().arguments)
            if (containsBuiltinCall(argument, id, builtins))
                return true;
    }
    else if (expression.isArray()) {
        for (std::size_t i = 0; i < expression.asArray().size(); ++i)
            if (containsBuiltinCall(expression.asArray().element(i), id, builtins))
                return true;
    }
    else if (expression.isList()) {
        for (const Expr& element : expression.asList().elements)
            if (containsBuiltinCall(element, id, builtins))
                return true;
    }
    return false;
}

} // namespace mmcal::solver
