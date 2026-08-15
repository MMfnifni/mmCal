// 実軸で安全な超越方程式反転
#include "transcendental_solver.hpp"

#include "mathematics/knowledge_context.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "polynomial_solver.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "symbolic/polynomial.hpp"

#include <optional>
#include <string>
#include <unordered_set>
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


[[nodiscard]] std::optional<FunctionSide> matchPeriodicFunctionSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const mathematics::MathRegistry& mathematics) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 1
        || !containsVariable(lhs.asCall().arguments[0], variable)
        || containsVariable(rhs, variable))
        return std::nullopt;

    const auto* definition = mathematics.findFunction(lhs.asCall().head);
    if (!definition || !definition->periodTurns || !definition->inverseFunction
        || definition->realGloballyInjective)
        return std::nullopt;
    if (definition->id != mathematics::FunctionId::Sin
        && definition->id != mathematics::FunctionId::Cos
        && definition->id != mathematics::FunctionId::Tan)
        return std::nullopt;
    return FunctionSide{definition, lhs.asCall().arguments[0], rhs};
}

[[nodiscard]] Expr angleValueFromTurns(
    const numeric::Rational& turns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Degree:
        return Expr{Number{turns * numeric::Rational{BigInt{360}}}};
    case mathematics::AngleUnit::Gradian:
        return Expr{Number{turns * numeric::Rational{BigInt{400}}}};
    case mathematics::AngleUnit::Radian: {
        const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
        if (!pi)
            error::throwCalcError(error::CalcErrorType::Internal, "Pi is not registered");
        const numeric::Rational coefficient = turns * numeric::Rational{BigInt{2}};
        if (coefficient == numeric::Rational{BigInt{1}})
            return Expr{pi->symbol};
        return Expr::call(
            builtins.symbol(BuiltinId::Multiply),
            {Expr{Number{coefficient}}, Expr{pi->symbol}});
    }
    }
    return Expr{Number{BigInt{0}}};
}

void collectSymbolNames(const Expr& expression, std::unordered_set<std::string>& names) {
    if (expression.isSymbol()) {
        names.insert(expression.asSymbol().name());
        return;
    }
    if (expression.isCall()) {
        for (const Expr& argument : expression.asCall().arguments)
            collectSymbolNames(argument, names);
        return;
    }
    if (expression.isArray()) {
        for (std::size_t i = 0; i < expression.asArray().size(); ++i)
            collectSymbolNames(expression.asArray().element(i), names);
        return;
    }
    if (expression.isList())
        for (const Expr& element : expression.asList().elements)
            collectSymbolNames(element, names);
}

[[nodiscard]] expression::Symbol freshIntegerParameter(
    const Expr& relation,
    const expression::Symbol& variable) {
    std::unordered_set<std::string> names;
    collectSymbolNames(relation, names);
    names.insert(variable.name());
    for (std::size_t index = 0;; ++index) {
        const std::string name = index == 0 ? "k" : "k" + std::to_string(index);
        if (!names.contains(name))
            return expression::Symbol{name};
    }
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

[[nodiscard]] Expr periodicTarget(
    Expr base,
    const Expr& period,
    const expression::Symbol& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr multiple = Expr::call(
        builtins.symbol(BuiltinId::Multiply), {period, Expr{parameter}});
    return simplifyExpr(
        Expr::call(builtins.symbol(BuiltinId::Add), {std::move(base), std::move(multiple)}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] bool exactEndpointOne(const Expr& value) {
    if (!value.isNumber() || !value.asNumber().isReal())
        return false;
    const numeric::Rational rational = value.asNumber().asReal().toRational();
    return rational == numeric::Rational{BigInt{1}}
        || rational == numeric::Rational{BigInt{-1}};
}

[[nodiscard]] bool exactZero(const Expr& value) {
    return value.isNumber() && value.asNumber().isReal() && value.asNumber().isZero();
}

[[nodiscard]] std::optional<SolutionSet> solvePeriodicTarget(
    const Expr& argument,
    const Expr& target,
    const expression::Symbol& variable,
    const expression::Symbol& parameter,
    const mathematics::AssumptionSet& remainingConditions,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr relation = Expr::call(
        builtins.symbol(BuiltinId::Equal), {argument, target});
    SolutionSet result = solveUnivariatePolynomialRelation(
        relation, variable, builtins, mathematics, angles);
    if (result.kind() != SolutionSetKind::Finite)
        return std::nullopt;

    std::vector<SolutionBranch> branches(result.branches().begin(), result.branches().end());
    for (SolutionBranch& branch : branches) {
        branch.freeVariables.push_back(
            SolverVariable{parameter, mathematics::NumericDomain::Integer});
        for (const auto& predicate : remainingConditions.predicates())
            branch.conditions.add(predicate);
    }
    return SolutionSet::finite(
        std::vector<SolverVariable>{result.variables().begin(), result.variables().end()},
        std::move(branches)).withAdditionalConditions(result.conditions());
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

std::optional<SolutionSet> solveRealPeriodicFunctionRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isEqualRelation(relation, builtins))
        return std::nullopt;

    const auto& sides = relation.asCall().arguments;
    auto matched = matchPeriodicFunctionSide(sides[0], sides[1], variable, mathematics);
    if (!matched)
        matched = matchPeriodicFunctionSide(sides[1], sides[0], variable, mathematics);
    if (!matched)
        return std::nullopt;

    // 初版はargumentが変数についてaffineで、一次係数がexact非零の場合に限定する。
    // symbolic coefficientや非線形argumentのparameter条件を不完全に返さない。
    const auto polynomial = symbolic::toExpressionPolynomial(
        matched->argument, variable, builtins, mathematics, angles);
    if (!polynomial || polynomial->degree() != 1)
        return std::nullopt;
    const Expr& linearCoefficient = polynomial->coefficient(1);
    if (!linearCoefficient.isNumber() || !linearCoefficient.asNumber().isReal()
        || linearCoefficient.asNumber().isZero())
        return std::nullopt;

    matched->rhs = simplifyExpr(
        matched->rhs, builtins, mathematics, angles, assumptions);
    const auto remainingConditions = proveRangeConditions(
        *matched->definition, matched->rhs, builtins, mathematics, assumptions);
    if (!remainingConditions)
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    const auto* inverse = mathematics.findFunction(*matched->definition->inverseFunction);
    if (!inverse)
        return std::nullopt;

    const expression::Symbol parameter = freshIntegerParameter(relation, variable);
    Expr period = simplifyExpr(
        angleValueFromTurns(*matched->definition->periodTurns, builtins, mathematics, angles),
        builtins, mathematics, angles, assumptions);
    const Expr inverseValue = simplifyExpr(
        Expr::call(inverse->symbol, {matched->rhs}),
        builtins, mathematics, angles, assumptions);

    std::vector<Expr> bases;
    switch (matched->definition->id) {
    case mathematics::FunctionId::Tan:
        bases.push_back(inverseValue);
        break;
    case mathematics::FunctionId::Cos:
        if (exactZero(matched->rhs)) {
            // cos[u]==0 は ±quarter-turn + full-turn*k の2branchを
            // quarter-turn + half-turn*k へ一意にまとめられる。
            bases.push_back(angleValueFromTurns(
                numeric::Rational{BigInt{1}, BigInt{4}}, builtins, mathematics, angles));
            period = angleValueFromTurns(
                numeric::Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles);
        }
        else {
            bases.push_back(inverseValue);
            if (!exactEndpointOne(matched->rhs))
                bases.push_back(simplifyExpr(
                    Expr::call(builtins.symbol(BuiltinId::Negate), {inverseValue}),
                    builtins, mathematics, angles, assumptions));
        }
        break;
    case mathematics::FunctionId::Sin:
        if (exactZero(matched->rhs)) {
            // sin[u]==0 は full-turnごとの2branchではなくhalf-turn*kがcanonical。
            bases.push_back(integer(0));
            period = angleValueFromTurns(
                numeric::Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles);
        }
        else {
            bases.push_back(inverseValue);
            if (!exactEndpointOne(matched->rhs)) {
                Expr halfTurn = angleValueFromTurns(
                    numeric::Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles);
                bases.push_back(simplifyExpr(
                    Expr::call(builtins.symbol(BuiltinId::Subtract),
                        {std::move(halfTurn), inverseValue}),
                    builtins, mathematics, angles, assumptions));
            }
        }
        break;
    default:
        return std::nullopt;
    }

    std::vector<SolutionBranch> combined;
    std::vector<SolverVariable> resultVariables;
    mathematics::AssumptionSet globalConditions;
    for (Expr& base : bases) {
        const Expr target = periodicTarget(
            std::move(base), period, parameter,
            builtins, mathematics, angles, assumptions);
        auto family = solvePeriodicTarget(
            matched->argument, target, variable, parameter, *remainingConditions,
            builtins, mathematics, angles);
        if (!family)
            return std::nullopt;
        if (resultVariables.empty())
            resultVariables.assign(family->variables().begin(), family->variables().end());
        for (const auto& predicate : family->conditions().predicates())
            globalConditions.add(predicate);
        combined.insert(combined.end(), family->branches().begin(), family->branches().end());
    }

    if (combined.empty())
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});
    return SolutionSet::finite(std::move(resultVariables), std::move(combined))
        .withAdditionalConditions(globalConditions);
}

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
