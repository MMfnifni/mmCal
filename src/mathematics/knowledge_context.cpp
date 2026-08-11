// 局所仮定を含む数学知識コンテキスト
#include "knowledge_context.hpp"

#include "numeric/number.hpp"

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
        if (permanentFacts.provablyNonReal
            && domain->domain != NumericDomain::Complex)
            return TruthValue::False;

        if (assumptions_.contains(predicate))
            return TruthValue::True;

        const ValueFacts valueFacts = facts(domain->expression);
        if (isSubdomainOf(valueFacts.domain, domain->domain))
            return TruthValue::True;
        if (valueFacts.provablyNonReal
            && domain->domain != NumericDomain::Complex)
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

    // a-b == 0 / != 0 は a == b / != b と完全に同値。
    // FullSimplify等が因子 (x-1) の非零性を x!=1 から証明できるよう、数値展開に依存しないこの最小限のrelation正規化を共有知識層で行う。
    if (relationPredicate.relation == RelationKind::Equal
        || relationPredicate.relation == RelationKind::NotEqual) {
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
                return *normalized;
        }
        if (isZero(relationPredicate.lhs)) {
            if (const auto normalized = differenceRelation(relationPredicate.rhs))
                return *normalized;
        }
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
            complement, relationPredicate.lhs, relationPredicate.rhs)))
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
