// 局所仮定を含む数学知識コンテキスト
#include "knowledge_context.hpp"

#include "definedness.hpp"
#include "numeric/number.hpp"
#include "symbolic/algebraic_expression.hpp"

#include <optional>

namespace mmcal::mathematics {
namespace {

using expression::Expr;
using numeric::Number;

[[nodiscard]] TruthValue boolTruth(bool value) noexcept {
    return value ? TruthValue::True : TruthValue::False;
}

[[nodiscard]] TruthValue proveNumberRelation(
    RelationKind relation,
    const Number& lhs,
    const Number& rhs) {
    if (relation == RelationKind::Equal)
        return boolTruth(lhs == rhs);
    if (relation == RelationKind::NotEqual)
        return boolTruth(!(lhs == rhs));

    // 複素数には順序関係を定義しない。
    if (!lhs.isReal() || !rhs.isReal())
        return TruthValue::Unknown;

    const auto& l = lhs.asReal();
    const auto& r = rhs.asReal();
    switch (relation) {
    case RelationKind::Less: return boolTruth(l < r);
    case RelationKind::LessEqual: return boolTruth(l <= r);
    case RelationKind::Greater: return boolTruth(l > r);
    case RelationKind::GreaterEqual: return boolTruth(l >= r);
    case RelationKind::Equal:
    case RelationKind::NotEqual:
        break;
    }
    return TruthValue::Unknown;
}


[[nodiscard]] TruthValue proveAlgebraicRelation(
    RelationKind relation,
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    const auto left = symbolic::exactAlgebraicValue(lhs, builtins, mathematics);
    const auto right = symbolic::exactAlgebraicValue(rhs, builtins, mathematics);
    if (!left || !right)
        return TruthValue::Unknown;

    if (relation == RelationKind::Equal || relation == RelationKind::NotEqual) {
        const auto equal = left->exactEquals(*right);
        if (!equal)
            return TruthValue::Unknown;
        return boolTruth(relation == RelationKind::Equal ? *equal : !*equal);
    }

    const auto order = left->exactRealCompare(*right);
    if (!order)
        return TruthValue::Unknown;
    switch (relation) {
    case RelationKind::Less:
        return boolTruth(*order == symbolic::AlgebraicOrder::Less);
    case RelationKind::LessEqual:
        return boolTruth(*order != symbolic::AlgebraicOrder::Greater);
    case RelationKind::Greater:
        return boolTruth(*order == symbolic::AlgebraicOrder::Greater);
    case RelationKind::GreaterEqual:
        return boolTruth(*order != symbolic::AlgebraicOrder::Less);
    case RelationKind::Equal:
    case RelationKind::NotEqual:
        break;
    }
    return TruthValue::Unknown;
}

[[nodiscard]] TruthValue proveAlgebraicDomain(
    const Expr& expression,
    NumericDomain domain,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    const auto algebraic = symbolic::exactAlgebraicValue(expression, builtins, mathematics);
    if (!algebraic)
        return TruthValue::Unknown;

    if (domain == NumericDomain::Complex)
        return TruthValue::True;
    if (domain == NumericDomain::Real
        && algebraic->domain() == symbolic::AlgebraicRootDomain::Real)
        return TruthValue::True;

    if (domain != NumericDomain::Integer && domain != NumericDomain::Rational)
        return TruthValue::Unknown;

    if (const auto parts = algebraic->exactRationalParts()) {
        if (!parts->second.isZero())
            return TruthValue::False;
        if (domain == NumericDomain::Rational)
            return TruthValue::True;
        return boolTruth(parts->first.denominator() == numeric::BigInt{1});
    }

    // degree>1のcanonical AlgebraicNumberはQ上既約minimal polynomialを持つ。
    if (algebraic->polynomial().size() > 2)
        return TruthValue::False;
    return TruthValue::Unknown;
}
[[nodiscard]] bool isZero(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] TruthValue proveSignRelation(
    RelationKind relation,
    const ValueFacts& facts,
    bool expressionOnLeft) noexcept {
    if (!expressionOnLeft) {
        switch (relation) {
        case RelationKind::Less: relation = RelationKind::Greater; break;
        case RelationKind::LessEqual: relation = RelationKind::GreaterEqual; break;
        case RelationKind::Greater: relation = RelationKind::Less; break;
        case RelationKind::GreaterEqual: relation = RelationKind::LessEqual; break;
        case RelationKind::Equal:
        case RelationKind::NotEqual:
            break;
        }
    }

    if (!facts.isProvablyReal())
        return facts.provablyNonReal ? TruthValue::False : TruthValue::Unknown;

    switch (relation) {
    case RelationKind::Equal:
        if (facts.sign == RealSign::Zero)
            return TruthValue::True;
        if (facts.sign == RealSign::Positive
            || facts.sign == RealSign::Negative
            || facts.sign == RealSign::NonZero)
            return TruthValue::False;
        return TruthValue::Unknown;
    case RelationKind::NotEqual:
        if (facts.sign == RealSign::Zero)
            return TruthValue::False;
        if (facts.sign == RealSign::Positive
            || facts.sign == RealSign::Negative
            || facts.sign == RealSign::NonZero)
            return TruthValue::True;
        return TruthValue::Unknown;
    case RelationKind::Less:
        if (facts.sign == RealSign::Negative)
            return TruthValue::True;
        if (facts.sign == RealSign::Zero
            || facts.sign == RealSign::Positive
            || facts.sign == RealSign::NonNegative)
            return TruthValue::False;
        return TruthValue::Unknown;
    case RelationKind::LessEqual:
        if (facts.sign == RealSign::Negative
            || facts.sign == RealSign::Zero
            || facts.sign == RealSign::NonPositive)
            return TruthValue::True;
        if (facts.sign == RealSign::Positive)
            return TruthValue::False;
        return TruthValue::Unknown;
    case RelationKind::Greater:
        if (facts.sign == RealSign::Positive)
            return TruthValue::True;
        if (facts.sign == RealSign::Zero
            || facts.sign == RealSign::Negative
            || facts.sign == RealSign::NonPositive)
            return TruthValue::False;
        return TruthValue::Unknown;
    case RelationKind::GreaterEqual:
        if (facts.sign == RealSign::Positive
            || facts.sign == RealSign::Zero
            || facts.sign == RealSign::NonNegative)
            return TruthValue::True;
        if (facts.sign == RealSign::Negative)
            return TruthValue::False;
        return TruthValue::Unknown;
    }
    return TruthValue::Unknown;
}


[[nodiscard]] RelationKind reversedRelation(RelationKind relation) noexcept {
    switch (relation) {
    case RelationKind::Less: return RelationKind::Greater;
    case RelationKind::LessEqual: return RelationKind::GreaterEqual;
    case RelationKind::Greater: return RelationKind::Less;
    case RelationKind::GreaterEqual: return RelationKind::LessEqual;
    case RelationKind::Equal: return RelationKind::Equal;
    case RelationKind::NotEqual: return RelationKind::NotEqual;
    }
    return relation;
}

struct RationalRangeBound final {
    numeric::Rational value;
    bool inclusive = false;
};

struct RationalRange final {
    std::optional<RationalRangeBound> lower;
    std::optional<RationalRangeBound> upper;
};

[[nodiscard]] std::optional<RationalRange> knownRealRange(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AssumptionSet& assumptions) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* function = mathematics.findFunction(expression.asCall().head);
    if (!function || function->realRangeRule == RealRangeRule::Unknown
        || function->definednessRule != FunctionDefinednessRule::Everywhere)
        return std::nullopt;

    const ValueFacts argument = inferValueFacts(
        expression.asCall().arguments[0], builtins, mathematics, assumptions);
    if (!argument.isProvablyReal())
        return std::nullopt;

    using numeric::BigInt;
    using numeric::Rational;
    const auto q = [](std::int64_t value) { return Rational{BigInt{value}}; };
    switch (function->realRangeRule) {
    case RealRangeRule::AllReal:
        return RationalRange{};
    case RealRangeRule::Positive:
        return RationalRange{RationalRangeBound{q(0), false}, std::nullopt};
    case RealRangeRule::NonNegative:
        return RationalRange{RationalRangeBound{q(0), true}, std::nullopt};
    case RealRangeRule::OpenMinusOneToOne:
        return RationalRange{
            RationalRangeBound{q(-1), false}, RationalRangeBound{q(1), false}};
    case RealRangeRule::OpenZeroToTwo:
        return RationalRange{
            RationalRangeBound{q(0), false}, RationalRangeBound{q(2), false}};
    case RealRangeRule::ClosedMinusOneToOne:
        return RationalRange{
            RationalRangeBound{q(-1), true}, RationalRangeBound{q(1), true}};
    case RealRangeRule::OneToInfinity:
        return RationalRange{RationalRangeBound{q(1), true}, std::nullopt};
    case RealRangeRule::Unknown:
        break;
    }
    return std::nullopt;
}

[[nodiscard]] TruthValue proveRealRangeRelation(
    RelationKind relation,
    const Expr& functionExpression,
    const Expr& boundExpression,
    bool functionOnLeft,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AssumptionSet& assumptions) {
    if (!boundExpression.isNumber() || !boundExpression.asNumber().isReal())
        return TruthValue::Unknown;
    const auto range = knownRealRange(
        functionExpression, builtins, mathematics, assumptions);
    if (!range)
        return TruthValue::Unknown;

    if (!functionOnLeft)
        relation = reversedRelation(relation);

    const numeric::Rational bound = boundExpression.asNumber().asReal().toRational();
    const auto lowerComparison = range->lower
        ? std::optional<int>{range->lower->value < bound ? -1
            : range->lower->value > bound ? 1 : 0}
        : std::nullopt;
    const auto upperComparison = range->upper
        ? std::optional<int>{range->upper->value < bound ? -1
            : range->upper->value > bound ? 1 : 0}
        : std::nullopt;

    switch (relation) {
    case RelationKind::Less:
        if (upperComparison && (*upperComparison < 0
                || (*upperComparison == 0 && !range->upper->inclusive)))
            return TruthValue::True;
        if (lowerComparison && *lowerComparison >= 0)
            return TruthValue::False;
        return TruthValue::Unknown;
    case RelationKind::LessEqual:
        if (upperComparison && *upperComparison <= 0)
            return TruthValue::True;
        if (lowerComparison && (*lowerComparison > 0
                || (*lowerComparison == 0 && !range->lower->inclusive)))
            return TruthValue::False;
        return TruthValue::Unknown;
    case RelationKind::Greater:
        if (lowerComparison && (*lowerComparison > 0
                || (*lowerComparison == 0 && !range->lower->inclusive)))
            return TruthValue::True;
        if (upperComparison && *upperComparison <= 0)
            return TruthValue::False;
        return TruthValue::Unknown;
    case RelationKind::GreaterEqual:
        if (lowerComparison && *lowerComparison >= 0)
            return TruthValue::True;
        if (upperComparison && (*upperComparison < 0
                || (*upperComparison == 0 && !range->upper->inclusive)))
            return TruthValue::False;
        return TruthValue::Unknown;
    case RelationKind::Equal: {
        const bool below = lowerComparison && (*lowerComparison > 0
            || (*lowerComparison == 0 && !range->lower->inclusive));
        const bool above = upperComparison && (*upperComparison < 0
            || (*upperComparison == 0 && !range->upper->inclusive));
        return below || above ? TruthValue::False : TruthValue::Unknown;
    }
    case RelationKind::NotEqual: {
        const bool below = lowerComparison && (*lowerComparison > 0
            || (*lowerComparison == 0 && !range->lower->inclusive));
        const bool above = upperComparison && (*upperComparison < 0
            || (*upperComparison == 0 && !range->upper->inclusive));
        return below || above ? TruthValue::True : TruthValue::Unknown;
    }
    }
    return TruthValue::Unknown;
}

} // namespace

KnowledgeContext::KnowledgeContext(
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AssumptionSet& assumptions)
    : builtins_(builtins), mathematics_(mathematics), assumptions_(assumptions) {}

ValueFacts KnowledgeContext::facts(const Expr& expression) const {
    return inferValueFacts(expression, builtins_, mathematics_, assumptions_);
}

TruthValue KnowledgeContext::prove(const Predicate& predicate) const {
    if (const auto* domain = std::get_if<DomainPredicate>(&predicate)) {
        // 明示Assumptionより恒久的な数学知識を優先する。 例えば I ∈ Real という矛盾したAssumptionで、Iの非実性を上書きしない。
        const ValueFacts permanentFacts = inferValueFacts(
            domain->expression, builtins_, mathematics_);
        if (isSubdomainOf(permanentFacts.domain, domain->domain))
            return TruthValue::True;
        if ((domain->domain == NumericDomain::Integer && permanentFacts.provablyNonInteger)
            || (domain->domain == NumericDomain::Rational && permanentFacts.provablyNonRational)
            || (permanentFacts.provablyNonReal && domain->domain != NumericDomain::Complex))
            return TruthValue::False;

        if (const TruthValue algebraic = proveAlgebraicDomain(
                domain->expression, domain->domain, builtins_, mathematics_);
            algebraic != TruthValue::Unknown)
            return algebraic;

        if (assumptions_.contains(predicate))
            return TruthValue::True;

        const ValueFacts valueFacts = facts(domain->expression);
        if (isSubdomainOf(valueFacts.domain, domain->domain))
            return TruthValue::True;
        if ((domain->domain == NumericDomain::Integer && valueFacts.provablyNonInteger)
            || (domain->domain == NumericDomain::Rational && valueFacts.provablyNonRational)
            || (valueFacts.provablyNonReal && domain->domain != NumericDomain::Complex))
            return TruthValue::False;

        return TruthValue::Unknown;
    }

    const auto& relationPredicate = std::get<RelationPredicate>(predicate);
    if (relationPredicate.lhs == relationPredicate.rhs) {
        if (relationPredicate.relation == RelationKind::Equal
            || relationPredicate.relation == RelationKind::LessEqual
            || relationPredicate.relation == RelationKind::GreaterEqual)
            return TruthValue::True;
        if (relationPredicate.relation == RelationKind::NotEqual
            || relationPredicate.relation == RelationKind::Less
            || relationPredicate.relation == RelationKind::Greater)
            return TruthValue::False;
    }

    if (relationPredicate.lhs.isNumber() && relationPredicate.rhs.isNumber()) {
        const TruthValue numeric = proveNumberRelation(
            relationPredicate.relation,
            relationPredicate.lhs.asNumber(),
            relationPredicate.rhs.asNumber());
        if (numeric != TruthValue::Unknown)
            return numeric;
    }

    if (const TruthValue algebraic = proveAlgebraicRelation(
            relationPredicate.relation, relationPredicate.lhs, relationPredicate.rhs,
            builtins_, mathematics_);
        algebraic != TruthValue::Unknown)
        return algebraic;

    // MathRegistryのreal rangeが有限境界を持つ函数は，実引数・everywhere-definedが
    // 証明できる場合だけ，sin[x] <= 1 のような全域boundを利用する。
    if (const TruthValue range = proveRealRangeRelation(
            relationPredicate.relation, relationPredicate.lhs, relationPredicate.rhs, true,
            builtins_, mathematics_, assumptions_);
        range != TruthValue::Unknown)
        return range;
    if (const TruthValue range = proveRealRangeRelation(
            relationPredicate.relation, relationPredicate.rhs, relationPredicate.lhs, false,
            builtins_, mathematics_, assumptions_);
        range != TruthValue::Unknown)
        return range;

    // MathRegistryに「定義域内では決して0にならない」と登録された函数は、
    // exp[z]!=0 のような非零性を局所実装へ重複記述せず証明できる。
    if ((relationPredicate.relation == RelationKind::Equal
            || relationPredicate.relation == RelationKind::NotEqual)
        && (isZero(relationPredicate.lhs) || isZero(relationPredicate.rhs))) {
        const Expr& candidate = isZero(relationPredicate.lhs)
            ? relationPredicate.rhs : relationPredicate.lhs;
        if (candidate.isCall()) {
            const auto* function = mathematics_.findFunction(candidate.asCall().head);
            if (function && function->zeroRule == FunctionZeroRule::NeverZero) {
                // zeroRuleは「定義される点では非零」であってdefinednessそのものではない。
                // gamma[x]/gamma[x]のpole等を消さないため、現在の仮定から候補式の
                // 定義条件をすべて証明できる場合に限って非零性を利用する。
                const auto conditions = expressionDomainConditions(candidate, builtins_, mathematics_);
                if (conditions) {
                    bool defined = true;
                    for (const Predicate& condition : conditions->predicates()) {
                        if (prove(condition) != TruthValue::True) {
                            defined = false;
                            break;
                        }
                    }
                    if (defined)
                        return relationPredicate.relation == RelationKind::NotEqual
                            ? TruthValue::True : TruthValue::False;
                }
            }
        }
    }

    // Realと「非実であることが証明済み」のComplexは同値になり得ない。
    // これにより x in Real の下で x != I をTrueと証明できる。
    // 近似的な虚部判定ではなくValueFactsのexact domain知識だけを使う。
    if (relationPredicate.relation == RelationKind::Equal
        || relationPredicate.relation == RelationKind::NotEqual) {
        const ValueFacts lhsFacts = facts(relationPredicate.lhs);
        const ValueFacts rhsFacts = facts(relationPredicate.rhs);
        const bool disjointRealComplex =
            (lhsFacts.isProvablyReal() && rhsFacts.provablyNonReal)
            || (rhsFacts.isProvablyReal() && lhsFacts.provablyNonReal);
        if (disjointRealComplex)
            return relationPredicate.relation == RelationKind::NotEqual
                ? TruthValue::True
                : TruthValue::False;
    }

    // 数値・構造だけで確定できない命題に限って、明示Assumptionを採用する。
    if (assumptions_.contains(predicate))
        return TruthValue::True;
    if (assumptions_.contains(relation(
            reversedRelation(relationPredicate.relation),
            relationPredicate.rhs, relationPredicate.lhs)))
        return TruthValue::True;

    // a-b relation 0 は a relation b へexactに正規化できる。ordered relationでは
    // 再帰先がReal性まで証明できた場合だけTrue/Falseになるため，Complex順序を捏造しない。
    auto differenceRelation = [&](const Expr& difference) -> std::optional<TruthValue> {
        if (!difference.isCall())
            return std::nullopt;
        const auto* definition = builtins_.find(difference.asCall().head);
        if (!definition || definition->id != evaluation::BuiltinId::Subtract
            || difference.asCall().arguments.size() != 2)
            return std::nullopt;
        const auto& arguments = difference.asCall().arguments;
        return prove(relation(
            relationPredicate.relation, arguments[0], arguments[1]));
    };

    if (isZero(relationPredicate.rhs)) {
        if (const auto normalized = differenceRelation(relationPredicate.lhs))
            if (*normalized != TruthValue::Unknown)
                return *normalized;
    }
    if (isZero(relationPredicate.lhs)) {
        if (const auto normalized = differenceRelation(relationPredicate.rhs))
            if (*normalized != TruthValue::Unknown)
                return *normalized;
    }

    // 明示前提の論理的な補集合も最低限認識する。solverのcase分岐で
    // a!=0 が与えられているのに a==0 側をUnknownとして残さないため。
    RelationKind complement = relationPredicate.relation;
    switch (relationPredicate.relation) {
    case RelationKind::Equal: complement = RelationKind::NotEqual; break;
    case RelationKind::NotEqual: complement = RelationKind::Equal; break;
    case RelationKind::Less: complement = RelationKind::GreaterEqual; break;
    case RelationKind::LessEqual: complement = RelationKind::Greater; break;
    case RelationKind::Greater: complement = RelationKind::LessEqual; break;
    case RelationKind::GreaterEqual: complement = RelationKind::Less; break;
    }
    if (assumptions_.contains(relation(
            complement, relationPredicate.lhs, relationPredicate.rhs))
        || assumptions_.contains(relation(
            reversedRelation(complement), relationPredicate.rhs, relationPredicate.lhs)))
        return TruthValue::False;

    if (isZero(relationPredicate.rhs))
        return proveSignRelation(
            relationPredicate.relation,
            facts(relationPredicate.lhs),
            true);
    if (isZero(relationPredicate.lhs))
        return proveSignRelation(
            relationPredicate.relation,
            facts(relationPredicate.rhs),
            false);

    return TruthValue::Unknown;
}

} // namespace mmcal::mathematics
