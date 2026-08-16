// 実軸で安全な超越方程式反転
#include "transcendental_solver.hpp"

#include "mathematics/knowledge_context.hpp"
#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
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


enum class CertifiedOrder {
    Less,
    Equal,
    Greater,
    Unknown
};

[[nodiscard]] CertifiedOrder certifiedConstantOrder(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs == rhs)
        return CertifiedOrder::Equal;

    const mathematics::AssumptionSet noAssumptions;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, noAssumptions};
    if (knowledge.prove(mathematics::relation(RelationKind::Less, lhs, rhs)) == TruthValue::True)
        return CertifiedOrder::Less;
    if (knowledge.prove(mathematics::relation(RelationKind::Greater, lhs, rhs)) == TruthValue::True)
        return CertifiedOrder::Greater;
    if (knowledge.prove(mathematics::relation(RelationKind::Equal, lhs, rhs)) == TruthValue::True)
        return CertifiedOrder::Equal;

    // 定数だけからなる超越式の大小は、guessではなくcertified enclosureが分離した場合だけ採用する。
    const approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    for (const std::size_t bits : {96U, 192U, 384U}) {
        try {
            const auto left = certified.enclose(lhs, bits);
            const auto right = certified.enclose(rhs, bits);
            if (!left || !right || !left->isReal() || !right->isReal())
                return CertifiedOrder::Unknown;
            const auto& l = left->asReal();
            const auto& r = right->asReal();
            if (l.upper() < r.lower())
                return CertifiedOrder::Less;
            if (l.lower() > r.upper())
                return CertifiedOrder::Greater;
            if (l.isPoint() && r.isPoint() && l.lower() == r.lower())
                return CertifiedOrder::Equal;
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            return CertifiedOrder::Unknown;
        }
    }
    return CertifiedOrder::Unknown;
}

[[nodiscard]] bool isVariableSquare(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 2)
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::Power)
        return false;
    const auto& arguments = expression.asCall().arguments;
    if (!arguments[0].isSymbol() || arguments[0].asSymbol() != variable
        || !arguments[1].isNumber() || !arguments[1].asNumber().isReal()
        || !arguments[1].asNumber().asReal().isInteger())
        return false;
    return arguments[1].asNumber().asReal().asInteger() == BigInt{2};
}

struct ExponentialPowerSide final {
    Expr base;
    Expr exponent;
    Expr rhs;
};

struct BaseLogSide final {
    Expr base;
    Expr argument;
    Expr rhs;
};

[[nodiscard]] std::optional<BaseLogSide> matchBaseLogSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 2
        || containsVariable(rhs, variable))
        return std::nullopt;
    const auto* definition = builtins.find(lhs.asCall().head);
    if (!definition || definition->id != BuiltinId::Log)
        return std::nullopt;
    const auto& arguments = lhs.asCall().arguments;
    if (containsVariable(arguments[0], variable)
        || !containsVariable(arguments[1], variable))
        return std::nullopt;
    return BaseLogSide{arguments[0], arguments[1], rhs};
}

[[nodiscard]] std::optional<ExponentialPowerSide> matchExponentialPowerSide(
    const Expr& lhs,
    const Expr& rhs,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    // 右辺はx^2等で変数依存してよい。ここでは左辺が「定数base ^ 変数依存指数」かだけを見る。
    if (!lhs.isCall() || lhs.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto* definition = builtins.find(lhs.asCall().head);
    if (!definition || definition->id != BuiltinId::Power)
        return std::nullopt;
    const auto& arguments = lhs.asCall().arguments;
    if (containsVariable(arguments[0], variable)
        || !containsVariable(arguments[1], variable))
        return std::nullopt;
    return ExponentialPowerSide{arguments[0], arguments[1], rhs};
}

[[nodiscard]] Expr lambertW(
    int branch,
    Expr argument,
    const evaluation::BuiltinRegistry& builtins) {
    if (branch == 0)
        return Expr::call(builtins.symbol(BuiltinId::LambertW), {std::move(argument)});
    return Expr::call(
        builtins.symbol(BuiltinId::LambertW),
        {integer(branch), std::move(argument)});
}

[[nodiscard]] SolutionBranch realBinding(
    const expression::Symbol& variable,
    Expr value,
    mathematics::AssumptionSet conditions = {}) {
    SolutionBranch branch;
    branch.bindings.push_back(SolutionBinding{variable, std::move(value)});
    branch.conditions = std::move(conditions);
    branch.bindingsCertifiedDomain = mathematics::NumericDomain::Real;
    return branch;
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

std::optional<SolutionSet> solveRealExponentialRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isEqualRelation(relation, builtins))
        return std::nullopt;

    const auto& sides = relation.asCall().arguments;
    auto matched = matchExponentialPowerSide(sides[0], sides[1], variable, builtins);
    if (!matched)
        matched = matchExponentialPowerSide(sides[1], sides[0], variable, builtins);
    if (!matched)
        return std::nullopt;

    mathematics::AssumptionSet proofAssumptions = assumptions;
    proofAssumptions.add(mathematics::elementOf(Expr{variable}, mathematics::NumericDomain::Real));
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, proofAssumptions};

    // a>0 かつ実指数なら principal Power[a,u] は exp[u log[a]] と一致し常に正。
    // 0^uや負baseを同じ規則で消すとdomain/branchを壊すため、証明できる場合だけ使う。
    if (knowledge.prove(mathematics::relation(
            RelationKind::Greater, matched->base, integer(0))) != TruthValue::True
        || knowledge.prove(mathematics::elementOf(
            matched->exponent, mathematics::NumericDomain::Real)) != TruthValue::True)
        return std::nullopt;

    if (exactZero(matched->rhs))
        return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

    // a^u == d でdが変数に依存しない場合は，a>0かつa!=1の実軸上で
    // u == log[a,d] と完全に同値。Powerを一般Evaluatorへ通さず，HoldAllのSolve内で
    // 証明済みの単射性だけを使って指数側をpolynomial solverへ渡す。
    if (!containsVariable(matched->rhs, variable)) {
        const TruthValue rhsReal = knowledge.prove(mathematics::elementOf(
            matched->rhs, mathematics::NumericDomain::Real));
        const TruthValue rhsPositive = knowledge.prove(mathematics::relation(
            RelationKind::Greater, matched->rhs, integer(0)));
        if (rhsReal == TruthValue::True
            && knowledge.prove(mathematics::relation(
                RelationKind::LessEqual, matched->rhs, integer(0))) == TruthValue::True)
            return SolutionSet::empty({SolverVariable{variable, mathematics::NumericDomain::Real}});

        const CertifiedOrder baseToOne = certifiedConstantOrder(
            matched->base, integer(1), builtins, mathematics, angles);
        if (baseToOne == CertifiedOrder::Equal) {
            const TruthValue rhsOne = knowledge.prove(mathematics::relation(
                RelationKind::Equal, matched->rhs, integer(1)));
            if (rhsOne == TruthValue::True)
                return SolutionSet::universal(
                    {SolverVariable{variable, mathematics::NumericDomain::Real}});
            if (rhsOne == TruthValue::False)
                return SolutionSet::empty(
                    {SolverVariable{variable, mathematics::NumericDomain::Real}});
            return std::nullopt;
        }

        if (rhsReal == TruthValue::True && rhsPositive == TruthValue::True
            && baseToOne != CertifiedOrder::Unknown) {
            Expr target = simplifyExpr(
                Expr::call(builtins.symbol(BuiltinId::Log),
                    {matched->base, matched->rhs}),
                builtins, mathematics, angles, proofAssumptions);
            Expr transformed = Expr::call(
                builtins.symbol(BuiltinId::Equal), {matched->exponent, std::move(target)});
            return solveUnivariatePolynomialRelation(
                transformed, variable, builtins, mathematics, angles);
        }
    }

    // variable-dependent RHSをLambert Wへ落とす初版は a^x == x^2 に限定する。
    // affine指数＋constant RHSは上で処理済み。一般のnonconstant P(x)はbranch/domain条件が増えるため別段階へ送る。
    if (!matched->exponent.isSymbol() || matched->exponent.asSymbol() != variable
        || !isVariableSquare(matched->rhs, variable, builtins))
        return std::nullopt;

    const CertifiedOrder baseToOne = certifiedConstantOrder(
        matched->base, integer(1), builtins, mathematics, angles);
    if (baseToOne == CertifiedOrder::Unknown)
        return std::nullopt;
    if (baseToOne == CertifiedOrder::Equal) {
        Expr transformed = Expr::call(
            builtins.symbol(BuiltinId::Equal), {integer(1), matched->rhs});
        return solveUnivariatePolynomialRelation(
            transformed, variable, builtins, mathematics, angles);
    }

    Expr logBase = simplifyExpr(
        Expr::call(builtins.symbol(BuiltinId::Log), {matched->base}),
        builtins, mathematics, angles, proofAssumptions);
    Expr magnitudeLog = baseToOne == CertifiedOrder::Greater
        ? logBase
        : simplifyExpr(
            Expr::call(builtins.symbol(BuiltinId::Negate), {logBase}),
            builtins, mathematics, angles, proofAssumptions);
    Expr halfMagnitude = simplifyExpr(
        Expr::call(builtins.symbol(BuiltinId::Divide), {magnitudeLog, integer(2)}),
        builtins, mathematics, angles, proofAssumptions);
    // base<1ではmagnitudeLog=-logBaseなので、-(magnitude/2)を機械的に作ると
    // --log[a]/2という非canonical形が残り得る。符号証明済みのlogBaseから直接構成する。
    Expr negativeHalfMagnitude = simplifyExpr(
        Expr::call(
            builtins.symbol(BuiltinId::Divide),
            {baseToOne == CertifiedOrder::Greater
                ? Expr::call(builtins.symbol(BuiltinId::Negate), {logBase})
                : logBase,
             integer(2)}),
        builtins, mathematics, angles, proofAssumptions);
    Expr scale = simplifyExpr(
        Expr::call(builtins.symbol(BuiltinId::Divide), {integer(-2), logBase}),
        builtins, mathematics, angles, proofAssumptions);

    auto scaledW = [&](int branch, Expr argument) {
        return simplifyExpr(
            Expr::call(
                builtins.symbol(BuiltinId::Multiply),
                {scale, lambertW(branch, std::move(argument), builtins)}),
            builtins, mathematics, angles, proofAssumptions);
    };

    std::vector<SolutionBranch> branches;
    // 符号がbase-1と反対側の根は、W_0(|log a|/2) が正実数上で常に存在するため無条件。
    branches.push_back(realBinding(variable, scaledW(0, halfMagnitude)));

    // 同符号側の2根は -|log a|/2 >= -1/e、すなわち |log a| <= 2/e のときだけ実在する。
    const auto* e = mathematics.findConstant(mathematics::ConstantId::E);
    if (!e)
        error::throwCalcError(error::CalcErrorType::Internal, "E is not registered");
    Expr threshold = simplifyExpr(
        Expr::call(builtins.symbol(BuiltinId::Divide), {integer(2), Expr{e->symbol}}),
        builtins, mathematics, angles, proofAssumptions);
    const CertifiedOrder branchCondition = certifiedConstantOrder(
        magnitudeLog, threshold, builtins, mathematics, angles);

    if (branchCondition == CertifiedOrder::Less) {
        branches.push_back(realBinding(variable, scaledW(0, negativeHalfMagnitude)));
        branches.push_back(realBinding(variable, scaledW(-1, negativeHalfMagnitude)));
    }
    else if (branchCondition == CertifiedOrder::Equal) {
        // branch point -1/E では W_0 = W_-1 = -1。同じ実根を二重に返さない。
        branches.push_back(realBinding(variable, scaledW(0, negativeHalfMagnitude)));
    }
    else if (branchCondition == CertifiedOrder::Unknown) {
        mathematics::AssumptionSet principalCondition;
        principalCondition.add(mathematics::relation(
            RelationKind::LessEqual, magnitudeLog, threshold));
        branches.push_back(realBinding(
            variable, scaledW(0, negativeHalfMagnitude), std::move(principalCondition)));

        // lower branch はbranch pointでprincipal branchと一致するので、重複回避のためstrict条件にする。
        mathematics::AssumptionSet lowerCondition;
        lowerCondition.add(mathematics::relation(
            RelationKind::Less, magnitudeLog, threshold));
        branches.push_back(realBinding(
            variable, scaledW(-1, negativeHalfMagnitude), std::move(lowerCondition)));
    }

    return SolutionSet::finite(
        {SolverVariable{variable, mathematics::NumericDomain::Real}}, std::move(branches));
}

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

    // log[b,u] == r は b>0, b!=1, r∈Real が証明できる場合，
    // u == b^r へ安全に反転できる。log2/log10はSolve normalizationで
    // このcanonical 2引数Logへ寄るため，同じknowledge pathを使う。
    auto baseLog = matchBaseLogSide(sides[0], sides[1], variable, builtins);
    if (!baseLog)
        baseLog = matchBaseLogSide(sides[1], sides[0], variable, builtins);
    if (baseLog) {
        mathematics::AssumptionSet proofAssumptions = assumptions;
        proofAssumptions.add(mathematics::elementOf(
            Expr{variable}, mathematics::NumericDomain::Real));
        const mathematics::KnowledgeContext knowledge{builtins, mathematics, proofAssumptions};
        const CertifiedOrder baseToOne = certifiedConstantOrder(
            baseLog->base, integer(1), builtins, mathematics, angles);
        if (baseToOne != CertifiedOrder::Unknown && baseToOne != CertifiedOrder::Equal
            && knowledge.prove(mathematics::relation(
                RelationKind::Greater, baseLog->base, integer(0))) == TruthValue::True
            && knowledge.prove(mathematics::elementOf(
                baseLog->rhs, mathematics::NumericDomain::Real)) == TruthValue::True) {
            Expr target = simplifyExpr(
                Expr::call(builtins.symbol(BuiltinId::Power),
                    {baseLog->base, baseLog->rhs}),
                builtins, mathematics, angles, proofAssumptions);
            Expr transformed = Expr::call(
                builtins.symbol(BuiltinId::Equal), {baseLog->argument, std::move(target)});
            return solveUnivariatePolynomialRelation(
                transformed, variable, builtins, mathematics, angles);
        }
    }

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
