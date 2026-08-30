// 実方程式の非存在・一意性証明
#include "real_equation_prover.hpp"

#include "mathematics/definedness.hpp"
#include "mathematics/knowledge_context.hpp"
#include "polynomial_solver.hpp"
#include "real_function_analysis.hpp"
#include "radical_solver.hpp"
#include "solve_constraints.hpp"
#include "solve_normalization.hpp"
#include "solver_limits.hpp"
#include "solver_support.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"
#include "transcendental_solver.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::solver {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::RealRangeRule;
using mathematics::RealSign;
using mathematics::RelationKind;
using mathematics::TruthValue;

[[nodiscard]] bool exactZero(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] bool isEqualRelation(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    return expression.isCall()
        && expression.asCall().arguments.size() == 2
        && builtins.isCallTo(expression, BuiltinId::Equal);
}

[[nodiscard]] Expr residualExpression(
    const Expr& relation,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto& sides = relation.asCall().arguments;
    return simplifyForSolve(
        Expr::call(
            builtins.symbol(BuiltinId::Subtract),
            {sides[0], sides[1]}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] bool provablyRealAndDefinedEverywhere(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::elementOf(
            expression, mathematics::NumericDomain::Real)) != TruthValue::True)
        return false;

    const auto conditions = mathematics::expressionDomainConditions(
        expression, builtins, mathematics);
    if (!conditions)
        return false;
    for (const auto& predicate : conditions->predicates())
        if (knowledge.prove(predicate) != TruthValue::True)
            return false;
    return true;
}

[[nodiscard]] std::optional<RealSign> strictSign(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(expression);
    if (!facts.isProvablyReal())
        return std::nullopt;
    if (facts.sign == RealSign::Positive || facts.sign == RealSign::Negative)
        return facts.sign;
    return std::nullopt;
}

[[nodiscard]] SolutionBranch realBinding(
    const expression::Symbol& variable,
    Expr value) {
    SolutionBranch branch;
    branch.bindings.push_back(SolutionBinding{variable, std::move(value)});
    branch.bindingsCertifiedDomain = mathematics::NumericDomain::Real;
    return branch;
}

[[nodiscard]] std::optional<SolutionSet> solveRealTarget(
    const Expr& argument,
    Expr target,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (argument.isSymbol() && argument.asSymbol().sameIdentity(variable))
        return SolutionSet::finite(
            {SolverVariable{variable, mathematics::NumericDomain::Real}},
            {realBinding(variable, std::move(target))});

    Expr relation = normalizeForSolve(
        Expr::call(
            builtins.symbol(BuiltinId::Equal), {argument, std::move(target)}),
        builtins, mathematics, angles, assumptions);

    if (auto radical = solveRadicalRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return radical;
    if (auto exponential = solveRealExponentialRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return exponential;
    if (auto injective = solveRealInjectiveFunctionRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return injective;

    SolveConstraints realConstraint;
    realConstraint.domain = mathematics::NumericDomain::Real;
    SolutionSet polynomial = applySolveConstraints(
        solveUnivariatePolynomialRelation(
            relation, variable, builtins, mathematics, angles),
        realConstraint, builtins, mathematics, angles);
    if (polynomial.kind() != SolutionSetKind::Unresolved)
        return polynomial;

    if (auto algebraic = solveRealAlgebraicPolynomialEquation(
            relation, variable, builtins, mathematics, angles))
        return algebraic;
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> singleUnconditionalRoot(
    const SolutionSet& solutions,
    const expression::Symbol& variable) {
    if (solutions.kind() != SolutionSetKind::Finite
        || !solutions.conditions().empty()
        || solutions.branches().size() != 1)
        return std::nullopt;
    const SolutionBranch& branch = solutions.branches().front();
    if (!branch.unconditional() || !branch.freeVariables.empty()
        || branch.bindings.size() != 1
        || !branch.bindings.front().variable.sameIdentity(variable))
        return std::nullopt;
    return branch.bindings.front().value;
}

[[nodiscard]] std::optional<SolutionSet> solveCriticalEquation(
    const Expr& derivative,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr relation = normalizeForSolve(
        Expr::call(
            builtins.symbol(BuiltinId::Equal), {derivative, integerExpr(0)}),
        builtins, mathematics, angles, assumptions);

    if (auto exponential = solveRealExponentialRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return exponential;
    if (auto periodic = solveRealPeriodicFunctionRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return periodic;
    if (auto injective = solveRealInjectiveFunctionRelation(
            relation, variable, builtins, mathematics, angles, assumptions))
        return injective;

    SolveConstraints realConstraint;
    realConstraint.domain = mathematics::NumericDomain::Real;
    SolutionSet polynomial = applySolveConstraints(
        solveUnivariatePolynomialRelation(
            relation, variable, builtins, mathematics, angles),
        realConstraint, builtins, mathematics, angles);
    if (polynomial.kind() != SolutionSetKind::Unresolved)
        return polynomial;
    if (auto algebraic = solveRealAlgebraicPolynomialEquation(
            relation, variable, builtins, mathematics, angles))
        return algebraic;
    return std::nullopt;
}

[[nodiscard]] mathematics::AssumptionSet realRangeConditions(
    const mathematics::FunctionDefinition& definition,
    const Expr& rhs) {
    mathematics::AssumptionSet result;
    result.add(mathematics::elementOf(rhs, mathematics::NumericDomain::Real));
    switch (definition.realRangeRule) {
    case RealRangeRule::Positive:
        result.add(mathematics::relation(RelationKind::Greater, rhs, integerExpr(0)));
        break;
    case RealRangeRule::NonNegative:
        result.add(mathematics::relation(RelationKind::GreaterEqual, rhs, integerExpr(0)));
        break;
    case RealRangeRule::OpenMinusOneToOne:
        result.add(mathematics::relation(RelationKind::Greater, rhs, integerExpr(-1)));
        result.add(mathematics::relation(RelationKind::Less, rhs, integerExpr(1)));
        break;
    case RealRangeRule::OpenZeroToTwo:
        result.add(mathematics::relation(RelationKind::Greater, rhs, integerExpr(0)));
        result.add(mathematics::relation(RelationKind::Less, rhs, integerExpr(2)));
        break;
    case RealRangeRule::ClosedMinusOneToOne:
        result.add(mathematics::relation(RelationKind::GreaterEqual, rhs, integerExpr(-1)));
        result.add(mathematics::relation(RelationKind::LessEqual, rhs, integerExpr(1)));
        break;
    case RealRangeRule::OneToInfinity:
        result.add(mathematics::relation(RelationKind::GreaterEqual, rhs, integerExpr(1)));
        break;
    case RealRangeRule::AllReal:
    case RealRangeRule::Unknown:
        break;
    }
    return result;
}

struct UnaryFunctionSide final {
    const mathematics::FunctionDefinition* definition = nullptr;
    Expr argument;
    Expr rhs;
};

[[nodiscard]] std::optional<UnaryFunctionSide> matchUnaryFunctionSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const mathematics::MathRegistry& mathematics) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 1
        || !symbolic::containsSymbol(lhs.asCall().arguments[0], variable)
        || symbolic::containsSymbol(rhs, variable))
        return std::nullopt;
    const auto* definition = mathematics.findFunction(lhs.asCall().head);
    if (!definition)
        return std::nullopt;
    return UnaryFunctionSide{definition, lhs.asCall().arguments[0], rhs};
}

[[nodiscard]] std::optional<SolutionSet> proveUnaryRangeOrAnchor(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto& sides = relation.asCall().arguments;
    auto matched = matchUnaryFunctionSide(
        sides[0], sides[1], variable, mathematics);
    if (!matched)
        matched = matchUnaryFunctionSide(
            sides[1], sides[0], variable, mathematics);
    if (!matched)
        return std::nullopt;

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::elementOf(
            matched->argument, mathematics::NumericDomain::Real)) != TruthValue::True)
        return std::nullopt;

    if (matched->definition->realRangeRule != RealRangeRule::Unknown) {
        const auto conditions = realRangeConditions(*matched->definition, matched->rhs);
        for (const auto& predicate : conditions.predicates())
            if (knowledge.prove(predicate) == TruthValue::False)
                return SolutionSet::empty(
                    {SolverVariable{variable, mathematics::NumericDomain::Real}});
    }

    // inverse builtinを持たない函数でも，実軸上global injectiveで既知のexact値が
    // targetと一致するなら，f(u)==f(c)をu==cへ完全に帰着できる。
    if (!matched->definition->realGloballyInjective)
        return std::nullopt;

    const std::array<Expr, 3> anchors{integerExpr(0), integerExpr(1), integerExpr(-1)};
    for (const Expr& anchor : anchors) {
        Expr value = simplifyForSolve(
            Expr::call(matched->definition->symbol, {anchor}),
            builtins, mathematics, angles, assumptions);
        if (knowledge.prove(mathematics::relation(
                RelationKind::Equal, value, matched->rhs)) != TruthValue::True)
            continue;
        return solveRealTarget(
            matched->argument, anchor, variable,
            builtins, mathematics, angles, assumptions);
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<SolutionSet> proveGlobalSignExclusion(
    const Expr& residual,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(residual);
    if (!facts.isProvablyReal())
        return std::nullopt;
    if (facts.sign == RealSign::Positive
        || facts.sign == RealSign::Negative
        || facts.sign == RealSign::NonZero)
        return SolutionSet::empty(
            {SolverVariable{variable, mathematics::NumericDomain::Real}});
    return std::nullopt;
}

[[nodiscard]] bool exactRootAt(
    const Expr& residual,
    const Expr& candidate,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr value = simplifyForSolve(
        symbolic::substituteSymbol(residual, variable, candidate),
        builtins, mathematics, angles, assumptions);
    if (exactZero(value))
        return true;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    return knowledge.prove(mathematics::relation(
        RelationKind::Equal, value, integerExpr(0))) == TruthValue::True;
}

[[nodiscard]] std::optional<SolutionSet> proveStrictMonotoneAnchor(
    const Expr& residual,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!provablyRealAndDefinedEverywhere(
            residual, builtins, mathematics, assumptions))
        return std::nullopt;

    Expr derivative = symbolic::differentiateExpression(
        residual, variable, builtins, mathematics, angles);
    if (containsBuiltinCall(derivative, BuiltinId::Derivative, builtins)
        || expressionNodeCount(derivative, limits::realProofNodes) > limits::realProofNodes)
        return std::nullopt;
    derivative = simplifyForSolve(
        std::move(derivative), builtins, mathematics, angles, assumptions);
    if (!provablyRealAndDefinedEverywhere(
            derivative, builtins, mathematics, assumptions)
        || !strictSign(derivative, builtins, mathematics, assumptions))
        return std::nullopt;

    const std::array<Expr, 3> candidates{integerExpr(0), integerExpr(1), integerExpr(-1)};
    for (const Expr& candidate : candidates)
        if (exactRootAt(
                residual, candidate, variable,
                builtins, mathematics, angles, assumptions))
            return SolutionSet::finite(
                {SolverVariable{variable, mathematics::NumericDomain::Real}},
                {realBinding(variable, candidate)});
    return std::nullopt;
}

[[nodiscard]] std::optional<SolutionSet> proveStrictConvexExtremum(
    const Expr& residual,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!provablyRealAndDefinedEverywhere(
            residual, builtins, mathematics, assumptions))
        return std::nullopt;

    Expr first = symbolic::differentiateExpression(
        residual, variable, builtins, mathematics, angles);
    if (containsBuiltinCall(first, BuiltinId::Derivative, builtins)
        || expressionNodeCount(first, limits::realProofNodes) > limits::realProofNodes)
        return std::nullopt;
    first = simplifyForSolve(
        std::move(first), builtins, mathematics, angles, assumptions);
    if (!provablyRealAndDefinedEverywhere(
            first, builtins, mathematics, assumptions))
        return std::nullopt;

    Expr second = symbolic::differentiateExpression(
        first, variable, builtins, mathematics, angles);
    if (containsBuiltinCall(second, BuiltinId::Derivative, builtins)
        || expressionNodeCount(second, limits::realProofNodes) > limits::realProofNodes)
        return std::nullopt;
    second = simplifyForSolve(
        std::move(second), builtins, mathematics, angles, assumptions);
    if (!provablyRealAndDefinedEverywhere(
            second, builtins, mathematics, assumptions))
        return std::nullopt;

    const auto curvature = strictSign(
        second, builtins, mathematics, assumptions);
    if (!curvature)
        return std::nullopt;

    std::optional<Expr> point;
    // strict curvatureが既にf'の単調性を証明しているので，exact anchorでf'=0を
    // 1点見つければcritical pointはそれだけで唯一である。一般solverを呼ぶ前に
    // 低コストな0,±1を試し，見つからない場合だけ既存exact solverへ委譲する。
    const std::array<Expr, 3> criticalAnchors{integerExpr(0), integerExpr(1), integerExpr(-1)};
    for (const Expr& candidate : criticalAnchors) {
        if (exactRootAt(
                first, candidate, variable,
                builtins, mathematics, angles, assumptions)) {
            point = candidate;
            break;
        }
    }
    if (!point) {
        const auto critical = solveCriticalEquation(
            first, variable, builtins, mathematics, angles, assumptions);
        if (!critical)
            return std::nullopt;
        point = singleUnconditionalRoot(*critical, variable);
        if (!point)
            return std::nullopt;
    }

    Expr value = simplifyForSolve(
        symbolic::substituteSymbol(residual, variable, *point),
        builtins, mathematics, angles, assumptions);
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto valueFacts = knowledge.facts(value);
    const bool zero = exactZero(value)
        || knowledge.prove(mathematics::relation(
            RelationKind::Equal, value, integerExpr(0))) == TruthValue::True;

    if (*curvature == RealSign::Positive) {
        if (valueFacts.isProvablyReal() && valueFacts.sign == RealSign::Positive)
            return SolutionSet::empty(
                {SolverVariable{variable, mathematics::NumericDomain::Real}});
        if (zero)
            return SolutionSet::finite(
                {SolverVariable{variable, mathematics::NumericDomain::Real}},
                {realBinding(variable, *point)});
    }
    else {
        if (valueFacts.isProvablyReal() && valueFacts.sign == RealSign::Negative)
            return SolutionSet::empty(
                {SolverVariable{variable, mathematics::NumericDomain::Real}});
        if (zero)
            return SolutionSet::finite(
                {SolverVariable{variable, mathematics::NumericDomain::Real}},
                {realBinding(variable, *point)});
    }
    return std::nullopt;
}


[[nodiscard]] bool rangeExcludesZero(
    const RealValueRange& range,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    if (range.lower.isSymbol() && range.lower.asSymbol().sameIdentity(infinity))
        return true;
    if (builtins.isCallTo(range.upper, BuiltinId::Negate)
        && range.upper.asCall().arguments.size() == 1
        && range.upper.asCall().arguments[0].isSymbol()
        && range.upper.asCall().arguments[0].asSymbol().sameIdentity(infinity))
        return true;

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const Expr zero = integerExpr(0);
    const TruthValue lowerPositive = knowledge.prove(mathematics::relation(
        RelationKind::Greater, range.lower, zero));
    if (lowerPositive == TruthValue::True)
        return true;
    const TruthValue lowerZero = knowledge.prove(mathematics::relation(
        RelationKind::Equal, range.lower, zero));
    if (lowerZero == TruthValue::True && !range.lowerInclusive)
        return true;

    const TruthValue upperNegative = knowledge.prove(mathematics::relation(
        RelationKind::Less, range.upper, zero));
    if (upperNegative == TruthValue::True)
        return true;
    const TruthValue upperZero = knowledge.prove(mathematics::relation(
        RelationKind::Equal, range.upper, zero));
    return upperZero == TruthValue::True && !range.upperInclusive;
}

[[nodiscard]] std::optional<int> realSignOfValue(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    if (value.isSymbol() && value.asSymbol().sameIdentity(infinity))
        return 1;
    if (builtins.isCallTo(value, BuiltinId::Negate)
        && value.asCall().arguments.size() == 1
        && value.asCall().arguments[0].isSymbol()
        && value.asCall().arguments[0].asSymbol().sameIdentity(infinity))
        return -1;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(value);
    if (!facts.isProvablyReal())
        return std::nullopt;
    if (facts.sign == RealSign::Positive)
        return 1;
    if (facts.sign == RealSign::Negative)
        return -1;
    if (facts.sign == RealSign::Zero)
        return 0;
    return std::nullopt;
}

[[nodiscard]] bool monotoneEndpointExcludesZero(
    const RealIntervalFunctionAnalysis& piece,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    if (piece.monotonicity != RealIntervalMonotonicity::Increasing
        && piece.monotonicity != RealIntervalMonotonicity::Decreasing)
        return false;

    const auto lowerSign = piece.lowerLimit
        ? realSignOfValue(*piece.lowerLimit, builtins, mathematics, infinity, assumptions)
        : std::nullopt;
    const auto upperSign = piece.upperLimit
        ? realSignOfValue(*piece.upperLimit, builtins, mathematics, infinity, assumptions)
        : std::nullopt;
    const bool lowerAttained = piece.domain.lower && piece.domain.lowerInclusive;
    const bool upperAttained = piece.domain.upper && piece.domain.upperInclusive;
    const bool lowerKnown = lowerSign.has_value();
    const bool upperKnown = upperSign.has_value();
    const int lower = lowerSign.value_or(0);
    const int upper = upperSign.value_or(0);

    if (piece.monotonicity == RealIntervalMonotonicity::Increasing) {
        if (lowerKnown && (lower > 0 || (lower == 0 && !lowerAttained)))
            return true;
        if (upperKnown && (upper < 0 || (upper == 0 && !upperAttained)))
            return true;
    }
    else {
        if (lowerKnown && (lower < 0 || (lower == 0 && !lowerAttained)))
            return true;
        if (upperKnown && (upper > 0 || (upper == 0 && !upperAttained)))
            return true;
    }
    return false;
}

[[nodiscard]] bool pointInAnalyzedInterval(
    const Expr& point,
    const RealDomainInterval& interval,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (interval.lower) {
        const TruthValue relation = knowledge.prove(mathematics::relation(
            interval.lowerInclusive ? RelationKind::GreaterEqual : RelationKind::Greater,
            point, *interval.lower));
        if (relation != TruthValue::True)
            return false;
    }
    if (interval.upper) {
        const TruthValue relation = knowledge.prove(mathematics::relation(
            interval.upperInclusive ? RelationKind::LessEqual : RelationKind::Less,
            point, *interval.upper));
        if (relation != TruthValue::True)
            return false;
    }
    return true;
}

void appendUniqueRoot(std::vector<Expr>& roots, const Expr& root) {
    if (std::find(roots.begin(), roots.end(), root) == roots.end())
        roots.push_back(root);
}

[[nodiscard]] std::optional<SolutionSet> proveByIntervalAnalysis(
    const Expr& residual,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    const RealFunctionAnalysis analysis = analyzeRealFunction(
        residual, variable, builtins, mathematics, angles, infinity, assumptions);
    if (!analysis.domainComplete)
        return std::nullopt;
    if (analysis.pieces.empty())
        return SolutionSet::empty(
            {SolverVariable{variable, mathematics::NumericDomain::Real}});

    std::vector<Expr> roots;
    const std::array<Expr, 3> simpleAnchors{integerExpr(0), integerExpr(1), integerExpr(-1)};
    for (const RealIntervalFunctionAnalysis& piece : analysis.pieces) {
        if ((piece.range && rangeExcludesZero(
                *piece.range, builtins, mathematics, infinity, assumptions))
            || monotoneEndpointExcludesZero(
                piece, builtins, mathematics, infinity, assumptions))
            continue;

        if (piece.monotonicity != RealIntervalMonotonicity::Increasing
            && piece.monotonicity != RealIntervalMonotonicity::Decreasing)
            return std::nullopt;

        std::optional<Expr> root;
        const auto tryCandidate = [&](const Expr& candidate) {
            if (root || !pointInAnalyzedInterval(
                    candidate, piece.domain, builtins, mathematics, assumptions))
                return;
            if (exactRootAt(
                    residual, candidate, variable,
                    builtins, mathematics, angles, assumptions))
                root = candidate;
        };

        if (piece.domain.lower && piece.domain.lowerInclusive)
            tryCandidate(*piece.domain.lower);
        if (piece.domain.upper && piece.domain.upperInclusive)
            tryCandidate(*piece.domain.upper);
        for (const Expr& candidate : simpleAnchors)
            tryCandidate(candidate);

        if (!root)
            return std::nullopt;
        appendUniqueRoot(roots, *root);
    }

    if (roots.empty())
        return SolutionSet::empty(
            {SolverVariable{variable, mathematics::NumericDomain::Real}});
    std::vector<SolutionBranch> branches;
    branches.reserve(roots.size());
    for (Expr& root : roots)
        branches.push_back(realBinding(variable, std::move(root)));
    return SolutionSet::finite(
        {SolverVariable{variable, mathematics::NumericDomain::Real}},
        std::move(branches));
}

} // namespace

std::optional<SolutionSet> solveRealEquationByProof(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions) {
    if (!isEqualRelation(relation, builtins)
        || expressionNodeCount(relation, limits::realProofNodes) > limits::realProofNodes)
        return std::nullopt;

    const mathematics::AssumptionSet proofAssumptions =
        withRealVariable(assumptions, variable);
    Expr residual = residualExpression(
        relation, builtins, mathematics, angles, proofAssumptions);

    // 最廉価の証明から順に試す。ValueFactsが全域で0を排除できれば，
    // 微分やcritical point探索を行う必要はない。
    if (auto sign = proveGlobalSignExclusion(
            residual, variable, builtins, mathematics, proofAssumptions))
        return sign;

    // Wolfram FunctionRange/Reduceと同様に，既知の実値域からtarget不在を
    // 証明できる場合は解集合を空にする。値域内でもinjectivityとexact anchorが
    // 揃えば逆函数builtinなしで唯一解まで閉じられる。
    if (auto range = proveUnaryRangeOrAnchor(
            relation, variable, builtins, mathematics, angles, proofAssumptions))
        return range;

    // 実函数の定義域を連結区間へ分解し，臨界点でさらに分割する。
    // 各pieceのexact endpoint limitとstrict monotonicityから0がrange外なら排除し，
    // range内でもexact anchorがあればそのpieceの唯一解として証明する。
    if (auto intervalProof = proveByIntervalAnalysis(
            residual, variable, builtins, mathematics, angles,
            infinitySymbol, proofAssumptions))
        return intervalProof;

    if (auto monotone = proveStrictMonotoneAnchor(
            residual, variable, builtins, mathematics, angles, proofAssumptions))
        return monotone;

    // f''が全実軸でstrict signを持ち，f'=0の唯一の点を既存exact solverで
    // 構成できる場合，そこはglobal extremumである。extremum値の符号だけで
    // NoSolution，0なら接する唯一解を証明する。
    return proveStrictConvexExtremum(
        residual, variable, builtins, mathematics, angles, proofAssumptions);
}

} // namespace mmcal::solver
