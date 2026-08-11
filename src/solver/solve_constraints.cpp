// 定義域・制約付きsolve
#include "solve_constraints.hpp"

#include "error/error_message.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/predicate.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"

#include <algorithm>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

namespace mmcal::solver {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::NumericDomain;
using mathematics::Predicate;
using mathematics::RelationKind;
using mathematics::TruthValue;

[[nodiscard]] std::optional<NumericDomain> domainFromSymbol(const expression::Symbol& symbol) {
    const std::string_view name = symbol.view();
    if (name == "Integer") return NumericDomain::Integer;
    if (name == "Rational") return NumericDomain::Rational;
    if (name == "Real") return NumericDomain::Real;
    if (name == "Complex") return NumericDomain::Complex;
    return std::nullopt;
}

[[nodiscard]] std::optional<RelationKind> relationKind(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    switch (definition->id) {
    case BuiltinId::Equal: return RelationKind::Equal;
    case BuiltinId::NotEqual: return RelationKind::NotEqual;
    case BuiltinId::Less: return RelationKind::Less;
    case BuiltinId::LessEqual: return RelationKind::LessEqual;
    case BuiltinId::Greater: return RelationKind::Greater;
    case BuiltinId::GreaterEqual: return RelationKind::GreaterEqual;
    default: return std::nullopt;
    }
}

[[nodiscard]] Expr zeroExpr() {
    return Expr{numeric::Number{numeric::BigInt{0}}};
}

[[nodiscard]] Predicate canonicalRelation(
    RelationKind kind,
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (rhs.isNumber() && rhs.asNumber().isZero())
        return mathematics::relation(kind, lhs, rhs);

    Expr difference = simplification::Simplifier{}.simplify(
        Expr::call(builtins.symbol(BuiltinId::Subtract), {lhs, rhs}),
        simplification::SimplificationContext{builtins, mathematics, angles});
    return mathematics::relation(kind, std::move(difference), zeroExpr());
}

void mergeDomain(std::optional<NumericDomain>& current, NumericDomain candidate) {
    // Integer ⊂ Rational ⊂ Real ⊂ Complex なので、複数指定は交差=狭い方。
    if (!current)
        current = candidate;
    else if (mathematics::isSubdomainOf(candidate, *current))
        current = candidate;
}

void parseSpec(
    const Expr& spec,
    SolveConstraints& result,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (spec.isSymbol()) {
        if (const auto domain = domainFromSymbol(spec.asSymbol())) {
            mergeDomain(result.domain, *domain);
            return;
        }
        error::throwCalcError(
            error::CalcErrorType::Type,
            "solve constraint symbol must be Integer, Rational, Real, or Complex");
    }

    if (spec.isArray()) {
        if (spec.asArray().rank() != 1)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "solve constraints must be a one-dimensional array");
        for (const Expr& item : spec.asArray().elements)
            parseSpec(item, result, builtins, mathematics, angles);
        return;
    }

    if (spec.isCall()) {
        if (const auto relation = relationKind(spec, builtins)) {
            const auto& arguments = spec.asCall().arguments;
            if (*relation == RelationKind::Less
                || *relation == RelationKind::LessEqual
                || *relation == RelationKind::Greater
                || *relation == RelationKind::GreaterEqual) {
                // 順序比較そのものがReal（またはその部分domain）を要求する。
                // x>0を単なるBoolean条件としてComplex ambient domainに残さない。
                mergeDomain(result.domain, NumericDomain::Real);
            }
            result.assumptions.add(canonicalRelation(
                *relation, arguments[0], arguments[1], builtins, mathematics, angles));
            return;
        }

        const auto* definition = builtins.find(spec.asCall().head);
        if (definition && definition->id == BuiltinId::LogicalAnd) {
            for (const Expr& argument : spec.asCall().arguments)
                parseSpec(argument, result, builtins, mathematics, angles);
            return;
        }
    }

    error::throwCalcError(
        error::CalcErrorType::Type,
        "solve constraints must be a domain, comparison, or one-dimensional array");
}

[[nodiscard]] Expr simplifyExpr(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
}

[[nodiscard]] Expr substituteExpr(
    const Expr& expression,
    std::span<const SolutionBinding> bindings) {
    if (expression.isSymbol()) {
        for (const SolutionBinding& binding : bindings)
            if (expression.asSymbol() == binding.variable)
                return binding.value;
        return expression;
    }
    if (expression.isCall()) {
        std::vector<Expr> arguments;
        arguments.reserve(expression.asCall().arguments.size());
        for (const Expr& argument : expression.asCall().arguments)
            arguments.push_back(substituteExpr(argument, bindings));
        return Expr::call(expression.asCall().head, std::move(arguments));
    }
    if (expression.isArray()) {
        std::vector<Expr> elements;
        elements.reserve(expression.asArray().elements.size());
        for (const Expr& element : expression.asArray().elements)
            elements.push_back(substituteExpr(element, bindings));
        return Expr::array(expression.asArray().shape, std::move(elements));
    }
    return expression;
}

[[nodiscard]] Predicate substitutePredicate(
    const Predicate& predicate,
    std::span<const SolutionBinding> bindings,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate)) {
        return mathematics::relation(
            relation->relation,
            simplifyExpr(substituteExpr(relation->lhs, bindings), builtins, mathematics, angles, assumptions),
            simplifyExpr(substituteExpr(relation->rhs, bindings), builtins, mathematics, angles, assumptions));
    }
    const auto& domain = std::get<mathematics::DomainPredicate>(predicate);
    return mathematics::elementOf(
        simplifyExpr(substituteExpr(domain.expression, bindings), builtins, mathematics, angles, assumptions),
        domain.domain);
}

void appendUnique(mathematics::AssumptionSet& target, const mathematics::AssumptionSet& source) {
    for (const Predicate& predicate : source.predicates())
        target.add(predicate);
}

[[nodiscard]] TruthValue proveAll(
    const mathematics::AssumptionSet& predicates,
    const mathematics::AssumptionSet& assumptions,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    bool unknown = false;
    for (const Predicate& predicate : predicates.predicates()) {
        const TruthValue value = knowledge.prove(predicate);
        if (value == TruthValue::False)
            return TruthValue::False;
        unknown = unknown || value == TruthValue::Unknown;
    }
    return unknown ? TruthValue::Unknown : TruthValue::True;
}

[[nodiscard]] std::optional<SolutionBranch> refineBranch(
    SolutionBranch branch,
    std::span<const SolverVariable> variables,
    const SolveConstraints& constraints,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    for (SolverVariable& freeVariable : branch.freeVariables) {
        const auto it = std::find_if(variables.begin(), variables.end(), [&](const SolverVariable& candidate) {
            return candidate.symbol == freeVariable.symbol;
        });
        if (it != variables.end())
            freeVariable.domain = it->domain;
    }

    // まずambient domainだけを前提に、solver自身が付けたbranch conditionを再検証する。
    // branch conditionを最初からAssumptionへ入れると「自分自身を根拠にTrue」となり、
    // Real制約によって x!=I が自明になった場合などを除去できない。
    mathematics::AssumptionSet processedAssumptions;
    for (const SolverVariable& variable : variables)
        processedAssumptions.add(mathematics::elementOf(Expr{variable.symbol}, variable.domain));

    mathematics::AssumptionSet pending;
    for (const Predicate& condition : branch.conditions.predicates()) {
        const mathematics::KnowledgeContext knowledge{
            builtins, mathematics, processedAssumptions};
        const TruthValue truth = knowledge.prove(condition);
        if (truth == TruthValue::False)
            return std::nullopt;
        if (truth == TruthValue::Unknown) {
            pending.add(condition);
            processedAssumptions.add(condition);
        }
    }

    mathematics::AssumptionSet domainKnowledgeAssumptions = processedAssumptions;
    appendUnique(domainKnowledgeAssumptions, constraints.assumptions);
    const mathematics::KnowledgeContext domainKnowledge{
        builtins, mathematics, domainKnowledgeAssumptions};

    for (const SolverVariable& variable : variables) {
        const auto binding = std::find_if(
            branch.bindings.begin(), branch.bindings.end(), [&](const SolutionBinding& candidate) {
                return candidate.variable == variable.symbol;
            });
        if (binding == branch.bindings.end())
            continue;

        // Complexはsolveの既定ambient domainであり、solverが構成した数値式を
        // さらに「Complexと証明できるか」で絞る必要はない。Real以下だけ検査する。
        if (variable.domain != NumericDomain::Complex) {
            const Predicate domainPredicate = mathematics::elementOf(binding->value, variable.domain);
            const TruthValue domainTruth = domainKnowledge.prove(domainPredicate);
            if (domainTruth == TruthValue::False)
                return std::nullopt;
            if (domainTruth == TruthValue::Unknown) {
                pending.add(domainPredicate);
                processedAssumptions.add(domainPredicate);
            }
        }
    }

    // user constraintはbindingを代入してから判定する。Unknownだけを結果へ残す。
    for (const Predicate& constraint : constraints.assumptions.predicates()) {
        const Predicate substituted = substitutePredicate(
            constraint, branch.bindings, builtins, mathematics, angles, processedAssumptions);
        const mathematics::KnowledgeContext branchKnowledge{
            builtins, mathematics, processedAssumptions};
        const TruthValue truth = branchKnowledge.prove(substituted);
        if (truth == TruthValue::False)
            return std::nullopt;
        if (truth == TruthValue::Unknown) {
            pending.add(substituted);
            processedAssumptions.add(substituted);
        }
    }

    branch.conditions = std::move(pending);
    return branch;
}

[[nodiscard]] std::vector<SolutionBranch> refineBranches(
    std::span<const SolutionBranch> branches,
    std::span<const SolverVariable> variables,
    const SolveConstraints& constraints,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<SolutionBranch> result;
    result.reserve(branches.size());
    for (const SolutionBranch& branch : branches)
        if (auto refined = refineBranch(branch, variables, constraints, builtins, mathematics, angles))
            result.push_back(std::move(*refined));
    return result;
}

[[nodiscard]] bool isUniversalFreeBranch(
    const SolutionBranch& branch,
    std::span<const SolverVariable> variables) {
    if (!branch.bindings.empty() || !branch.conditions.empty()
        || branch.freeVariables.size() != variables.size())
        return false;
    for (const SolverVariable& variable : variables) {
        const auto iterator = std::find_if(
            branch.freeVariables.begin(), branch.freeVariables.end(),
            [&](const SolverVariable& freeVariable) {
                return freeVariable.symbol == variable.symbol
                    && freeVariable.domain == variable.domain;
            });
        if (iterator == branch.freeVariables.end())
            return false;
    }
    return true;
}

[[nodiscard]] std::vector<SolverVariable> constrainedVariables(
    std::span<const SolverVariable> variables,
    const std::optional<NumericDomain>& requestedDomain) {
    std::vector<SolverVariable> result(variables.begin(), variables.end());
    if (!requestedDomain)
        return result;

    for (SolverVariable& variable : result) {
        // domain指定は既存ambient domainとの交差として解釈する。
        // 例: 不等式solverがRealを要求しているとき、Complex指定でRealをComplexへ
        // 広げることはしない。Integer指定ならRealの部分domainなのでIntegerへ狭める。
        if (mathematics::isSubdomainOf(*requestedDomain, variable.domain))
            variable.domain = *requestedDomain;
        else if (!mathematics::isSubdomainOf(variable.domain, *requestedDomain))
            variable.domain = *requestedDomain;
    }
    return result;
}

[[nodiscard]] SolutionSet outcomeFromCase(
    std::vector<SolverVariable> variables,
    const SolutionCase& item,
    const SolveConstraints& constraints,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (item.outcome) {
    case SolutionSetKind::Empty:
        return SolutionSet::empty(std::move(variables));
    case SolutionSetKind::Universal: {
        mathematics::AssumptionSet conditions = item.conditions;
        appendUnique(conditions, constraints.assumptions);
        return SolutionSet::universal(std::move(variables), std::move(conditions));
    }
    case SolutionSetKind::Finite: {
        auto branches = refineBranches(
            item.branches, variables, constraints, builtins, mathematics, angles);
        return branches.empty()
            ? SolutionSet::empty(std::move(variables))
            : SolutionSet::finite(std::move(variables), std::move(branches));
    }
    case SolutionSetKind::Unresolved: {
        mathematics::AssumptionSet conditions = item.conditions;
        appendUnique(conditions, constraints.assumptions);
        return SolutionSet::unresolved(std::move(variables), std::move(conditions));
    }
    case SolutionSetKind::Conditional:
        break;
    }
    return SolutionSet::unresolved(std::move(variables), constraints.assumptions);
}

} // namespace

SolveConstraints parseSolveConstraints(
    const Expr& spec,
    std::span<const expression::Symbol> variables,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (variables.empty())
        error::throwCalcError(error::CalcErrorType::Type, "solve requires at least one variable");
    SolveConstraints result;
    parseSpec(spec, result, builtins, mathematics, angles);
    return result;
}

SolutionSet applySolveConstraints(
    SolutionSet solutions,
    const SolveConstraints& constraints,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const mathematics::AssumptionSet originalGlobal = solutions.conditions();
    std::vector<SolverVariable> variables = constrainedVariables(
        solutions.variables(), constraints.domain);

    auto preserveGlobal = [&](SolutionSet result) {
        return result.withAdditionalConditions(originalGlobal);
    };

    switch (solutions.kind()) {
    case SolutionSetKind::Empty:
        return preserveGlobal(SolutionSet::empty(std::move(variables)));
    case SolutionSetKind::Finite: {
        auto branches = refineBranches(
            solutions.branches(), variables, constraints, builtins, mathematics, angles);
        if (branches.empty())
            return preserveGlobal(SolutionSet::empty(variables));
        if (branches.size() == 1 && isUniversalFreeBranch(branches.front(), variables))
            return preserveGlobal(SolutionSet::universal(variables));
        return preserveGlobal(SolutionSet::finite(variables, std::move(branches)));
    }
    case SolutionSetKind::Universal: {
        mathematics::AssumptionSet conditions = originalGlobal;
        appendUnique(conditions, constraints.assumptions);
        return SolutionSet::universal(std::move(variables), std::move(conditions));
    }
    case SolutionSetKind::Unresolved: {
        mathematics::AssumptionSet conditions = originalGlobal;
        appendUnique(conditions, constraints.assumptions);
        return SolutionSet::unresolved(std::move(variables), std::move(conditions));
    }
    case SolutionSetKind::Conditional:
        break;
    }

    std::vector<SolutionCase> remaining;
    for (const SolutionCase& item : solutions.cases()) {
        mathematics::AssumptionSet caseKnowledge = constraints.assumptions;
        appendUnique(caseKnowledge, originalGlobal);
        const TruthValue truth = proveAll(item.conditions, caseKnowledge, builtins, mathematics);
        if (truth == TruthValue::False)
            continue;
        if (truth == TruthValue::True) {
            SolutionSet selected = outcomeFromCase(
                variables, item, constraints, builtins, mathematics, angles);
            return preserveGlobal(std::move(selected));
        }

        SolutionCase refined = item;
        if (refined.outcome == SolutionSetKind::Finite) {
            refined.branches = refineBranches(
                item.branches, variables, constraints, builtins, mathematics, angles);
            if (refined.branches.empty())
                refined.outcome = SolutionSetKind::Empty;
        }
        remaining.push_back(std::move(refined));
    }

    SolutionSet result = remaining.empty()
        ? SolutionSet::empty(variables)
        : SolutionSet::conditional(variables, std::move(remaining));
    return preserveGlobal(std::move(result));
}

} // namespace mmcal::solver
