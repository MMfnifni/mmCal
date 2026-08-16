// 式の符号・実数性などの事実推論
#include "value_facts.hpp"

#include "assumption_set.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "symbolic/algebraic_number.hpp"

#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mmcal::mathematics {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] bool isAtMostReal(NumericDomain domain) noexcept {
    return domain == NumericDomain::Integer
        || domain == NumericDomain::Rational
        || domain == NumericDomain::Real;
}

[[nodiscard]] NumericDomain widestArithmeticDomain(
    NumericDomain lhs,
    NumericDomain rhs) noexcept {
    if (lhs == NumericDomain::Unknown || rhs == NumericDomain::Unknown)
        return NumericDomain::Unknown;
    if (lhs == NumericDomain::Complex || rhs == NumericDomain::Complex)
        return NumericDomain::Complex;
    if (lhs == NumericDomain::Real || rhs == NumericDomain::Real)
        return NumericDomain::Real;
    if (lhs == NumericDomain::Rational || rhs == NumericDomain::Rational)
        return NumericDomain::Rational;
    return NumericDomain::Integer;
}

[[nodiscard]] RealSign negateSign(RealSign sign) noexcept {
    switch (sign) {
    case RealSign::Negative: return RealSign::Positive;
    case RealSign::Zero: return RealSign::Zero;
    case RealSign::Positive: return RealSign::Negative;
    case RealSign::NonPositive: return RealSign::NonNegative;
    case RealSign::NonNegative: return RealSign::NonPositive;
    case RealSign::NonZero: return RealSign::NonZero;
    case RealSign::Unknown: return RealSign::Unknown;
    }
    return RealSign::Unknown;
}

[[nodiscard]] RealSign multiplySigns(RealSign lhs, RealSign rhs) noexcept {
    if (lhs == RealSign::Zero || rhs == RealSign::Zero)
        return RealSign::Zero;

    const auto exactSign = [](RealSign sign) -> int {
        if (sign == RealSign::Positive)
            return 1;
        if (sign == RealSign::Negative)
            return -1;
        return 0;
    };

    const int l = exactSign(lhs);
    const int r = exactSign(rhs);
    if (l != 0 && r != 0)
        return l == r ? RealSign::Positive : RealSign::Negative;

    if ((lhs == RealSign::NonNegative && rhs == RealSign::NonNegative)
        || (lhs == RealSign::NonPositive && rhs == RealSign::NonPositive))
        return RealSign::NonNegative;
    if ((lhs == RealSign::NonNegative && rhs == RealSign::NonPositive)
        || (lhs == RealSign::NonPositive && rhs == RealSign::NonNegative))
        return RealSign::NonPositive;

    return RealSign::Unknown;
}

[[nodiscard]] RealSign addSigns(RealSign lhs, RealSign rhs) noexcept {
    if (lhs == RealSign::Zero)
        return rhs;
    if (rhs == RealSign::Zero)
        return lhs;
    if (lhs == RealSign::Positive && rhs == RealSign::Positive)
        return RealSign::Positive;
    if (lhs == RealSign::Negative && rhs == RealSign::Negative)
        return RealSign::Negative;
    if ((lhs == RealSign::Positive || lhs == RealSign::NonNegative)
        && (rhs == RealSign::Positive || rhs == RealSign::NonNegative))
        return lhs == RealSign::Positive || rhs == RealSign::Positive
            ? RealSign::Positive
            : RealSign::NonNegative;
    if ((lhs == RealSign::Negative || lhs == RealSign::NonPositive)
        && (rhs == RealSign::Negative || rhs == RealSign::NonPositive))
        return lhs == RealSign::Negative || rhs == RealSign::Negative
            ? RealSign::Negative
            : RealSign::NonPositive;
    return RealSign::Unknown;
}

[[nodiscard]] ValueFacts factsForNumber(const Number& number) {
    if (number.isComplex()) {
        ValueFacts facts{NumericDomain::Complex, RealSign::Unknown, true, true};
        facts.provablyNonInteger = true;
        facts.provablyNonRational = true;
        return facts;
    }

    const auto& real = number.asReal();
    const bool integer = real.isInteger();
    NumericDomain domain = integer ? NumericDomain::Integer : NumericDomain::Rational;
    RealSign sign = RealSign::Zero;
    if (real.isNegative())
        sign = RealSign::Negative;
    else if (!real.isZero())
        sign = RealSign::Positive;
    ValueFacts facts{domain, sign, true, false};
    facts.provablyNonInteger = !integer;
    return facts;
}

[[nodiscard]] ValueFacts factsForSymbol(
    const expression::Symbol& symbol,
    const MathRegistry& mathematics) {
    if (const ConstantDefinition* constant = mathematics.findConstant(symbol)) {
        RealSign sign = RealSign::Unknown;
        if (constant->properties.positive)
            sign = RealSign::Positive;
        ValueFacts facts{
            constant->properties.real ? NumericDomain::Real : NumericDomain::Complex,
            sign,
            constant->properties.exact,
            !constant->properties.real
        };
        const bool irrational = constant->properties.irrational
            || constant->properties.arithmeticClass == ArithmeticClass::Transcendental;
        facts.provablyNonInteger = irrational || !constant->properties.real;
        facts.provablyNonRational = irrational || !constant->properties.real;
        return facts;
    }

    // 未束縛Symbolを無条件にRealと仮定すると、将来のassumption systemや複素変数を不当に狭める。
    // mmCalの通常入力では未定義SymbolはEvaluatorがNameErrorにするため、数学知識層ではdomainを証明できないものを正直にUnknownとして残す。
    // ただしSymbolという式そのものはmachine approximationではないのでexact性は保つ。
    return ValueFacts{NumericDomain::Unknown, RealSign::Unknown, true, false};
}


[[nodiscard]] bool isZeroExpression(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] RealSign signFromRelation(
    RelationKind relation,
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

    switch (relation) {
    case RelationKind::Equal: return RealSign::Zero;
    case RelationKind::NotEqual: return RealSign::NonZero;
    case RelationKind::Less: return RealSign::Negative;
    case RelationKind::LessEqual: return RealSign::NonPositive;
    case RelationKind::Greater: return RealSign::Positive;
    case RelationKind::GreaterEqual: return RealSign::NonNegative;
    }
    return RealSign::Unknown;
}

[[nodiscard]] ValueFacts applyAssumptions(
    const Expr& expression,
    ValueFacts facts,
    const AssumptionSet* assumptions) {
    if (!assumptions)
        return facts;

    for (const Predicate& predicate : assumptions->predicates()) {
        if (const auto* domain = std::get_if<DomainPredicate>(&predicate)) {
            if (!(domain->expression == expression))
                continue;

            // 前提は「少なくともこの集合に属する」という知識ではなく、solver変数のambient domainを指定するものとして扱う。
            // 既により狭いdomainを証明できているなら、その情報を失わない。
            if (facts.provablyNonReal
                && domain->domain != NumericDomain::Complex)
                continue;

            if (facts.domain == NumericDomain::Unknown
                || isSubdomainOf(domain->domain, facts.domain))
                facts.domain = domain->domain;
            continue;
        }

        const auto* relationPredicate = std::get_if<RelationPredicate>(&predicate);
        if (!relationPredicate)
            continue;

        bool expressionOnLeft = false;
        if (relationPredicate->lhs == expression
            && isZeroExpression(relationPredicate->rhs))
            expressionOnLeft = true;
        else if (!(relationPredicate->rhs == expression)
            || !isZeroExpression(relationPredicate->lhs))
            continue;

        const RealSign assumedSign = signFromRelation(
            relationPredicate->relation,
            expressionOnLeft);
        if (assumedSign == RealSign::Unknown)
            continue;

        // x != 0 はComplexでも意味を持ち、xがRealであることを意味しない。
        // 一方、順序比較は実数上の命題であり、x == 0 もxが実数(実際は0)と分かる。
        const bool notEqualOnly = relationPredicate->relation == RelationKind::NotEqual;
        if (notEqualOnly) {
            if (facts.isProvablyReal())
                facts.sign = RealSign::NonZero;
            continue;
        }

        // 順序関係や0との等値比較は実数に対する命題。
        // 恒久知識で非実数と証明済みなら矛盾した前提で上書きせず、知識側を保つ。
        if (facts.provablyNonReal)
            continue;
        if (facts.domain == NumericDomain::Unknown
            || facts.domain == NumericDomain::Complex)
            facts.domain = NumericDomain::Real;
        facts.sign = assumedSign;
        facts.provablyNonReal = false;
    }

    return facts;
}

[[nodiscard]] ValueFacts sqrtFacts(const ValueFacts& input) {
    if (!input.isNumeric())
        return {};
    if (input.isProvablyNegativeReal())
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, input.exact, true};
    if (input.isProvablyReal()) {
        RealSign sign = RealSign::Unknown;
        if (input.sign == RealSign::Zero)
            sign = RealSign::Zero;
        else if (input.sign == RealSign::Positive)
            sign = RealSign::Positive;
        else if (input.sign == RealSign::NonNegative)
            sign = RealSign::NonNegative;
        return ValueFacts{NumericDomain::Real, sign, input.exact, false};
    }
    return ValueFacts{NumericDomain::Complex, RealSign::Unknown, input.exact, input.provablyNonReal};
}

[[nodiscard]] const ValueFacts& childFacts(
    const std::unordered_map<const void*, ValueFacts>& facts,
    const Expr& expression) {
    return facts.at(expression.identity());
}

[[nodiscard]] ValueFacts inferCallFacts(
    const expression::CallExpr& call,
    const std::unordered_map<const void*, ValueFacts>& facts,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    const auto* builtin = builtins.find(call.head);
    if (!builtin)
        return {};

    const auto argument = [&](std::size_t index) -> const ValueFacts& {
        return childFacts(facts, call.arguments[index]);
    };

    switch (builtin->id) {
    case BuiltinId::Negate:
        if (call.arguments.size() != 1)
            return {};
        {
            ValueFacts result = argument(0);
            if (result.isProvablyReal())
                result.sign = negateSign(result.sign);
            return result;
        }

    case BuiltinId::Add:
    case BuiltinId::Subtract:
        if (call.arguments.empty())
            return {};
        {
            ValueFacts result = argument(0);
            if (builtin->id == BuiltinId::Subtract && call.arguments.size() != 2)
                return {};
            for (std::size_t i = 1; i < call.arguments.size(); ++i) {
                ValueFacts rhs = argument(i);
                if (builtin->id == BuiltinId::Subtract)
                    rhs.sign = negateSign(rhs.sign);
                const bool lhsWasNonReal = result.provablyNonReal;
                const bool lhsWasReal = result.isProvablyReal();
                const bool rhsIsReal = rhs.isProvablyReal();
                result.domain = widestArithmeticDomain(result.domain, rhs.domain);
                result.exact = result.exact && rhs.exact;

                // 非実数 + 実数 は虚部を消せないので非実数のまま。
                // ただし非実数同士の和は I + (-I) = 0 のように実数へ戻り得るため、「どちらかが非実数なら非実数」という危険な推論はしない。
                result.provablyNonReal =
                    (lhsWasNonReal && rhsIsReal)
                    || (rhs.provablyNonReal && lhsWasReal);
                result.sign = isAtMostReal(result.domain)
                    ? addSigns(result.sign, rhs.sign)
                    : RealSign::Unknown;
            }
            return result;
        }

    case BuiltinId::Multiply:
        if (call.arguments.empty())
            return {};
        {
            ValueFacts result = argument(0);
            for (std::size_t i = 1; i < call.arguments.size(); ++i) {
                const ValueFacts& rhs = argument(i);
                const bool lhsWasNonReal = result.provablyNonReal;
                const bool lhsWasReal = result.isProvablyReal();
                const RealSign lhsSign = result.sign;
                const bool rhsIsReal = rhs.isProvablyReal();
                result.domain = widestArithmeticDomain(result.domain, rhs.domain);
                result.exact = result.exact && rhs.exact;

                // 非実数 * 非0実数 は非実数。ただし I*I=-1 のように複素同士では実数へ戻ることがある。
                // また0を掛ければ当然0になるので、実数因子が非0だと証明できる場合だけ非実性を伝播する。
                const bool lhsRealNonZero = lhsWasReal
                    && (lhsSign == RealSign::Positive
                        || lhsSign == RealSign::Negative
                        || lhsSign == RealSign::NonZero);
                const bool rhsRealNonZero = rhsIsReal
                    && (rhs.sign == RealSign::Positive
                        || rhs.sign == RealSign::Negative
                        || rhs.sign == RealSign::NonZero);
                result.provablyNonReal =
                    (lhsWasNonReal && rhsRealNonZero)
                    || (rhs.provablyNonReal && lhsRealNonZero);
                result.sign = isAtMostReal(result.domain)
                    ? multiplySigns(lhsSign, rhs.sign)
                    : RealSign::Unknown;
            }
            return result;
        }

    case BuiltinId::Divide:
        if (call.arguments.size() != 2)
            return {};
        {
            const ValueFacts& lhs = argument(0);
            const ValueFacts& rhs = argument(1);
            NumericDomain domain = widestArithmeticDomain(lhs.domain, rhs.domain);
            if (domain == NumericDomain::Integer)
                domain = NumericDomain::Rational;
            const bool lhsRealNonZero = lhs.isProvablyReal()
                && (lhs.sign == RealSign::Positive
                    || lhs.sign == RealSign::Negative
                    || lhs.sign == RealSign::NonZero);
            const bool rhsRealNonZero = rhs.isProvablyReal()
                && (rhs.sign == RealSign::Positive
                    || rhs.sign == RealSign::Negative
                    || rhs.sign == RealSign::NonZero);
            const bool provablyNonReal =
                (lhs.provablyNonReal && rhsRealNonZero)
                || (rhs.provablyNonReal && lhsRealNonZero);
            return ValueFacts{
                domain,
                isAtMostReal(domain) ? multiplySigns(lhs.sign, rhs.sign) : RealSign::Unknown,
                lhs.exact && rhs.exact,
                provablyNonReal
            };
        }

    case BuiltinId::Power:
        if (call.arguments.size() != 2)
            return {};
        {
            const ValueFacts& base = argument(0);
            const Expr& exponentExpression = call.arguments[1];
            if (!exponentExpression.isNumber()
                || !exponentExpression.asNumber().isReal())
                return ValueFacts{NumericDomain::Complex, RealSign::Unknown, base.exact && argument(1).exact, false};

            const auto exponentRational = exponentExpression.asNumber().asReal().toRational();
            if (exponentRational == numeric::Rational{BigInt{1}, BigInt{2}}) {
                ValueFacts result = sqrtFacts(base);
                result.exact = result.exact && argument(1).exact;
                return result;
            }

            if (!exponentExpression.asNumber().asReal().isInteger())
                return ValueFacts{NumericDomain::Complex, RealSign::Unknown, base.exact && argument(1).exact, false};

            const BigInt& exponent = exponentExpression.asNumber().asReal().asInteger();
            if (exponent.isZero())
                return ValueFacts{NumericDomain::Integer, RealSign::Positive, true, false};

            NumericDomain domain = base.domain;
            if (exponent.isNegative() && domain == NumericDomain::Integer)
                domain = NumericDomain::Rational;

            RealSign sign = RealSign::Unknown;
            if (base.isProvablyReal()) {
                if (base.sign == RealSign::Zero)
                    sign = exponent.isNegative() ? RealSign::Unknown : RealSign::Zero;
                else if (base.sign == RealSign::Positive)
                    sign = RealSign::Positive;
                else if (base.sign == RealSign::Negative) {
                    const bool even = exponent.abs().trailingZeroBits() != 0;
                    sign = even ? RealSign::Positive : RealSign::Negative;
                }
            }
            return ValueFacts{domain, sign, base.exact && argument(1).exact, false};
        }

    case BuiltinId::Root: {
        // canonical Rootはminimal polynomialへ縮約済み。degree>1ならQ上既約多項式の根なので
        // Rationalではあり得ず，したがってIntegerでもない。
        const bool irrational = call.algebraicValue && call.algebraicValue->polynomial().size() > 2;
        if (call.arguments.size() == 2) {
            ValueFacts rootFacts{NumericDomain::Real, RealSign::Unknown, true, false};
            rootFacts.provablyNonInteger = irrational;
            rootFacts.provablyNonRational = irrational;
            return rootFacts;
        }
        if (call.arguments.size() == 3 && call.arguments[2].isSymbol()
            && call.arguments[2].asSymbol().view() == "Complex") {
            ValueFacts rootFacts{NumericDomain::Complex, RealSign::Unknown, true, false};
            rootFacts.provablyNonInteger = irrational;
            rootFacts.provablyNonRational = irrational;
            return rootFacts;
        }
        return {};
    }

    case BuiltinId::Cbrt:
        if (call.arguments.size() != 1 || !argument(0).isProvablyReal())
            return {};
        return ValueFacts{
            NumericDomain::Real, argument(0).sign, argument(0).exact, false};

    case BuiltinId::Hypot:
        if (call.arguments.size() != 2
            || !argument(0).isProvablyReal() || !argument(1).isProvablyReal())
            return {};
        {
            RealSign sign = RealSign::NonNegative;
            if (argument(0).sign == RealSign::Zero && argument(1).sign == RealSign::Zero)
                sign = RealSign::Zero;
            else {
                const auto nonZero = [](RealSign value) noexcept {
                    return value == RealSign::Positive || value == RealSign::Negative
                        || value == RealSign::NonZero;
                };
                if (nonZero(argument(0).sign) || nonZero(argument(1).sign))
                    sign = RealSign::Positive;
            }
            return ValueFacts{
                NumericDomain::Real, sign, argument(0).exact && argument(1).exact, false};
        }

    case BuiltinId::Cis:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        return ValueFacts{
            NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Polar:
        if (call.arguments.size() != 2
            || !argument(0).isNumeric() || !argument(1).isNumeric())
            return {};
        return ValueFacts{
            NumericDomain::Complex, RealSign::Unknown,
            argument(0).exact && argument(1).exact, false};

    case BuiltinId::DegreeToRadian:
    case BuiltinId::DegreeToGradian:
    case BuiltinId::RadianToDegree:
    case BuiltinId::RadianToGradian:
    case BuiltinId::GradianToDegree:
    case BuiltinId::GradianToRadian:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        // 正の実定数倍なのでdomain/sign/exactnessをそのまま保てる。
        return argument(0);

    case BuiltinId::NextPow2:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        return ValueFacts{NumericDomain::Integer, RealSign::Unknown, true, false};

    case BuiltinId::IsPrime:
        return {}; // Boolean result.
    case BuiltinId::NextPrime:
    case BuiltinId::PreviousPrime:
    case BuiltinId::Totient:
        return ValueFacts{NumericDomain::Integer, RealSign::Positive, true, false};
    case BuiltinId::FactorInteger:
        return {}; // List result.

    case BuiltinId::Zeta:
    case BuiltinId::Digamma:
    case BuiltinId::Trigamma:
    case BuiltinId::IncompleteBeta:
        // 詳細なdomain/signは各函数の定義域・parameterに依存するため、現段階では推論しない。
        return {};

    case BuiltinId::Sqrt:
        if (call.arguments.size() != 1)
            return {};
        return sqrtFacts(argument(0));

    case BuiltinId::Abs:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        return ValueFacts{
            NumericDomain::Real,
            argument(0).sign == RealSign::Zero ? RealSign::Zero : RealSign::NonNegative,
            argument(0).exact,
            false};

    case BuiltinId::Sign:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (argument(0).isProvablyReal()) {
            RealSign sign = RealSign::Unknown;
            switch (argument(0).sign) {
            case RealSign::Negative: sign = RealSign::Negative; break;
            case RealSign::Zero: sign = RealSign::Zero; break;
            case RealSign::Positive: sign = RealSign::Positive; break;
            case RealSign::NonZero: sign = RealSign::NonZero; break;
            case RealSign::NonPositive:
            case RealSign::NonNegative:
            case RealSign::Unknown:
                break;
            }
            return ValueFacts{NumericDomain::Real, sign, argument(0).exact, false};
        }
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Re:
    case BuiltinId::Im:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        return ValueFacts{NumericDomain::Real, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Conj:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        return ValueFacts{
            argument(0).domain,
            argument(0).isProvablyReal() ? argument(0).sign : RealSign::Unknown,
            argument(0).exact,
            argument(0).provablyNonReal};

    case BuiltinId::Sin:
    case BuiltinId::Cos:
    case BuiltinId::Tan:
    case BuiltinId::Cot:
    case BuiltinId::Sec:
    case BuiltinId::Csc:
        if (call.arguments.size() != 1)
            return {};
        {
            const ValueFacts& input = argument(0);
            if (input.isProvablyReal())
                return ValueFacts{NumericDomain::Real, RealSign::Unknown, input.exact, false};
            if (input.domain == NumericDomain::Complex)
                return ValueFacts{NumericDomain::Complex, RealSign::Unknown, input.exact, false};
            return {};
        }

    case BuiltinId::Sinh:
    case BuiltinId::Tanh:
    case BuiltinId::Csch:
    case BuiltinId::Coth:
        if (call.arguments.size() != 1)
            return {};
        {
            const ValueFacts& input = argument(0);
            if (input.isProvablyReal())
                return ValueFacts{NumericDomain::Real, input.sign, input.exact, false};
            if (input.domain == NumericDomain::Complex)
                return ValueFacts{NumericDomain::Complex, RealSign::Unknown, input.exact, false};
            return {};
        }

    case BuiltinId::Cosh:
    case BuiltinId::Sech:
        if (call.arguments.size() != 1)
            return {};
        {
            const ValueFacts& input = argument(0);
            if (input.isProvablyReal())
                return ValueFacts{NumericDomain::Real, RealSign::Positive, input.exact, false};
            if (input.domain == NumericDomain::Complex)
                return ValueFacts{NumericDomain::Complex, RealSign::Unknown, input.exact, false};
            return {};
        }

    case BuiltinId::Atan:
    case BuiltinId::Asinh:
        if (call.arguments.size() != 1)
            return {};
        {
            const ValueFacts& input = argument(0);
            if (input.isProvablyReal())
                return ValueFacts{NumericDomain::Real, input.sign, input.exact, false};
            if (input.domain == NumericDomain::Complex)
                return ValueFacts{NumericDomain::Complex, RealSign::Unknown, input.exact, false};
            return {};
        }

    case BuiltinId::Asin:
    case BuiltinId::Acos:
    case BuiltinId::Acosh:
    case BuiltinId::Atanh:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (call.arguments[0].isNumber() && call.arguments[0].asNumber().isReal()) {
            const auto value = call.arguments[0].asNumber().asReal().toRational();
            const numeric::Rational minusOne{numeric::BigInt{-1}};
            const numeric::Rational one{numeric::BigInt{1}};
            const bool realPrincipal =
                ((builtin->id == BuiltinId::Asin || builtin->id == BuiltinId::Acos)
                    && minusOne <= value && value <= one)
                || (builtin->id == BuiltinId::Atanh && minusOne < value && value < one)
                || (builtin->id == BuiltinId::Acosh && value >= one);
            if (realPrincipal)
                return ValueFacts{NumericDomain::Real,
                    builtin->id == BuiltinId::Atanh ? argument(0).sign : RealSign::Unknown,
                    argument(0).exact, false};
        }
        // 実入力でもprincipal値が複素数になる領域を持つため、範囲条件を証明できない一般式ではComplexを保守的な上界とする。
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Atan2:
        if (call.arguments.size() != 2
            || !argument(0).isProvablyReal() || !argument(1).isProvablyReal())
            return {};
        return ValueFacts{
            NumericDomain::Real, RealSign::Unknown, argument(0).exact && argument(1).exact, false};

    case BuiltinId::Arg:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        // principal Argは非零複素数を実数 (-Pi, Pi] へ写す。
        // 0で未定義かどうかはEvaluatorが値を見て判定する。
        return ValueFacts{NumericDomain::Real, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Log:
        if (call.arguments.size() == 1 && argument(0).isNumeric()) {
            const ValueFacts& input = argument(0);
            if (input.isProvablyReal() && input.sign == RealSign::Positive) {
                RealSign logSign = RealSign::Unknown;
                if (call.arguments[0].isNumber() && call.arguments[0].asNumber().isReal()) {
                    const auto value = call.arguments[0].asNumber().asReal().toRational();
                    const numeric::Rational one{BigInt{1}};
                    logSign = value == one ? RealSign::Zero
                        : value > one ? RealSign::Positive : RealSign::Negative;
                }
                return ValueFacts{NumericDomain::Real, logSign, input.exact, false};
            }
            if (input.isProvablyNegativeReal())
                return ValueFacts{NumericDomain::Complex, RealSign::Unknown, input.exact, true};
            return ValueFacts{NumericDomain::Complex, RealSign::Unknown, input.exact, false};
        }
        if (call.arguments.size() == 2
            && argument(0).isNumeric() && argument(1).isNumeric()) {
            const ValueFacts& base = argument(0);
            const ValueFacts& value = argument(1);
            const bool exact = base.exact && value.exact;
            if (base.isProvablyReal() && base.sign == RealSign::Positive
                && value.isProvablyReal() && value.sign == RealSign::Positive)
                return ValueFacts{NumericDomain::Real, RealSign::Unknown, exact, false};
            return ValueFacts{NumericDomain::Complex, RealSign::Unknown, exact, false};
        }
        return {};

    case BuiltinId::Log2:
    case BuiltinId::Log10:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (argument(0).isProvablyReal() && argument(0).sign == RealSign::Positive)
            return ValueFacts{NumericDomain::Real, RealSign::Unknown, argument(0).exact, false};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Gamma:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (argument(0).isProvablyReal()) {
            const RealSign sign = argument(0).sign == RealSign::Positive
                ? RealSign::Positive : RealSign::Unknown;
            return ValueFacts{NumericDomain::Real, sign, argument(0).exact, false};
        }
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::LogGamma:
        if (call.arguments.size() == 1 && argument(0).isProvablyReal())
            return ValueFacts{NumericDomain::Real, RealSign::Unknown, argument(0).exact, false};
        return {};

    case BuiltinId::Erf:
    case BuiltinId::Erfc:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (argument(0).isProvablyReal())
            return ValueFacts{NumericDomain::Real, RealSign::Unknown, argument(0).exact, false};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (argument(0).isProvablyReal())
            return ValueFacts{NumericDomain::Real, RealSign::Unknown, argument(0).exact, false};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Hypergeometric1F1:
        if (call.arguments.size() != 3
            || !argument(0).isNumeric() || !argument(1).isNumeric() || !argument(2).isNumeric())
            return {};
        if (argument(0).isProvablyReal()
            && argument(1).isProvablyReal()
            && argument(2).isProvablyReal())
            return ValueFacts{NumericDomain::Real, RealSign::Unknown,
                argument(0).exact && argument(1).exact && argument(2).exact, false};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown,
            argument(0).exact && argument(1).exact && argument(2).exact, false};

    case BuiltinId::Hypergeometric2F1:
        if (call.arguments.size() != 4)
            return {};
        // principal branchでは実parameterでもz>1で複素値になり得る。
        // 現Facts層ではbranch cut判定を捏造せずComplex可能性を保持する。
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown,
            argument(0).exact && argument(1).exact && argument(2).exact && argument(3).exact, false};

    case BuiltinId::EllipticF:
    case BuiltinId::EllipticE:
        if (call.arguments.size() != 2)
            return {};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown,
            argument(0).exact && argument(1).exact, false};

    case BuiltinId::EllipticPi:
        if (call.arguments.size() != 3)
            return {};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown,
            argument(0).exact && argument(1).exact && argument(2).exact, false};

    case BuiltinId::ExponentialIntegralEi:
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
    case BuiltinId::LogarithmicIntegralLi:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (builtin->id == BuiltinId::SineIntegralSi && argument(0).isProvablyReal())
            return ValueFacts{NumericDomain::Real, RealSign::Unknown, argument(0).exact, false};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Polylog:
        if (call.arguments.size() != 2 || !argument(0).isNumeric() || !argument(1).isNumeric())
            return {};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown,
            argument(0).exact && argument(1).exact, false};

    case BuiltinId::Beta:
    case BuiltinId::BetaLog:
        if (call.arguments.size() != 2
            || !argument(0).isProvablyReal() || !argument(1).isProvablyReal())
            return {};
        return ValueFacts{NumericDomain::Real,
            builtin->id == BuiltinId::Beta ? RealSign::Positive : RealSign::Unknown,
            argument(0).exact && argument(1).exact, false};

    case BuiltinId::LambertW: {
        if (call.arguments.empty() || call.arguments.size() > 2)
            return {};
        const std::size_t valueIndex = call.arguments.size() - 1;
        const ValueFacts& value = argument(valueIndex);
        bool principal = call.arguments.size() == 1;
        if (call.arguments.size() == 2 && call.arguments[0].isNumber()
            && call.arguments[0].asNumber().isReal()
            && call.arguments[0].asNumber().asReal().isInteger())
            principal = call.arguments[0].asNumber().asReal().asInteger().isZero();
        // W_0(x) は x>=0 で実かつ非負。負側のreal branch判定には -1/e 境界の
        // 証明が必要なので、ここで「負の実数ならReal」とは推測しない。
        if (principal && value.isProvablyReal()
            && (value.sign == RealSign::Zero || value.sign == RealSign::Positive
                || value.sign == RealSign::NonNegative))
            return ValueFacts{NumericDomain::Real, value.sign, value.exact, false};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, value.exact, false};
    }

    case BuiltinId::Exp:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        {
            const ValueFacts& input = argument(0);
            if (input.isProvablyReal())
                return ValueFacts{NumericDomain::Real, RealSign::Positive, input.exact, false};
            return ValueFacts{NumericDomain::Complex, RealSign::Unknown, input.exact, false};
        }

    case BuiltinId::Expm1:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (argument(0).isProvablyReal())
            return ValueFacts{NumericDomain::Real, argument(0).sign, argument(0).exact, false};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Log1p:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (argument(0).isProvablyNonNegativeReal())
            return ValueFacts{NumericDomain::Real, argument(0).sign, argument(0).exact, false};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc:
    case BuiltinId::Sinhc:
    case BuiltinId::Tanhc:
    case BuiltinId::Expc:
        if (call.arguments.size() != 1 || !argument(0).isNumeric())
            return {};
        if (argument(0).isProvablyReal())
            return ValueFacts{NumericDomain::Real, RealSign::Unknown, argument(0).exact, false};
        return ValueFacts{NumericDomain::Complex, RealSign::Unknown, argument(0).exact, false};

    case BuiltinId::NumericalApproximation:
        if (call.arguments.empty())
            return {};
        {
            ValueFacts result = argument(0);
            result.exact = false;
            return result;
        }

    case BuiltinId::BitAnd:
    case BuiltinId::BitOr:
    case BuiltinId::BitXor:
    case BuiltinId::BitNot:
    case BuiltinId::BitShiftLeft:
    case BuiltinId::BitShiftRight:
    case BuiltinId::BitLength:
    case BuiltinId::BitCount:
    case BuiltinId::BitGet:
        return ValueFacts{NumericDomain::Integer, RealSign::Unknown, true, false};

    case BuiltinId::Fma:
        if (call.arguments.size() != 3 || !argument(0).isNumeric()
            || !argument(1).isNumeric() || !argument(2).isNumeric())
            return {};
        {
            const bool real = argument(0).isProvablyReal()
                && argument(1).isProvablyReal() && argument(2).isProvablyReal();
            return ValueFacts{real ? NumericDomain::Real : NumericDomain::Complex,
                RealSign::Unknown,
                argument(0).exact && argument(1).exact && argument(2).exact, false};
        }

    case BuiltinId::Clamp:
        if (call.arguments.size() != 3 || !argument(0).isProvablyReal()
            || !argument(1).isProvablyReal() || !argument(2).isProvablyReal())
            return {};
        return ValueFacts{NumericDomain::Real, RealSign::Unknown,
            argument(0).exact && argument(1).exact && argument(2).exact, false};

    case BuiltinId::Proj:
        if (call.arguments.size() != 1)
            return {};
        return argument(0);

    case BuiltinId::Factorial:
        if (call.arguments.size() == 1 && argument(0).domain == NumericDomain::Integer)
            return ValueFacts{NumericDomain::Integer, RealSign::NonNegative, true, false};
        return {};

    case BuiltinId::UnitApplied:
        if (call.arguments.empty())
            return {};
        return argument(0);

    case BuiltinId::Derivative:
    case BuiltinId::SymbolicIntegral:
    case BuiltinId::Limit:
    case BuiltinId::Floor:
    case BuiltinId::Ceil:
    case BuiltinId::Trunc:
    case BuiltinId::Round:
    case BuiltinId::Frac:
    case BuiltinId::Gcd:
    case BuiltinId::Lcm:
    case BuiltinId::Mod:
    case BuiltinId::Rem:
    case BuiltinId::Quotient:
    case BuiltinId::Permutation:
    case BuiltinId::Combination:
    case BuiltinId::Fibonacci:
    case BuiltinId::DiscreteFourierTransform:
    case BuiltinId::FastFourierTransform:
    case BuiltinId::InverseFourierTransform:
    case BuiltinId::Convolution:
    case BuiltinId::Transpose:
    case BuiltinId::MatrixAdd:
    case BuiltinId::MatrixMultiply:
    case BuiltinId::Determinant:
    case BuiltinId::Inverse:
    case BuiltinId::Rref:
    case BuiltinId::Rank:
    case BuiltinId::SolveLinear:
    case BuiltinId::NullSpace:
    case BuiltinId::LuDecomposition:
    case BuiltinId::QrDecomposition:
    case BuiltinId::SingularValueDecomposition:
    case BuiltinId::Eigenvalues:
    case BuiltinId::Eigenvectors:
    case BuiltinId::Eigensystem:
    case BuiltinId::ConjugateTranspose:
    case BuiltinId::Length:
    case BuiltinId::NumericDerivative:
    case BuiltinId::NumericIntegral:
    case BuiltinId::Precision:
    case BuiltinId::Accuracy:
    case BuiltinId::Explain:
    case BuiltinId::Rationalize:
    case BuiltinId::GeneralizedBinomial:
    case BuiltinId::FallingFactorial:
    case BuiltinId::RisingFactorial:
    case BuiltinId::RandSeed:
    case BuiltinId::Rand:
    case BuiltinId::RandInt:
    case BuiltinId::Choice:
    case BuiltinId::RandN:
    case BuiltinId::Sum:
    case BuiltinId::Product:
    case BuiltinId::Min:
    case BuiltinId::Max:
    case BuiltinId::Mean:
    case BuiltinId::Median:
    case BuiltinId::Mode:
    case BuiltinId::Quantile:
    case BuiltinId::Percentile:
    case BuiltinId::VariancePopulation:
    case BuiltinId::VarianceSample:
    case BuiltinId::StddevPopulation:
    case BuiltinId::StddevSample:
    case BuiltinId::GeometricMean:
    case BuiltinId::HarmonicMean:
    case BuiltinId::Rms:
    case BuiltinId::MedianAbsoluteDeviation:
    case BuiltinId::MeanAbsoluteDeviation:
    case BuiltinId::Skewness:
    case BuiltinId::KurtosisPopulation:
    case BuiltinId::KurtosisSample:
    case BuiltinId::CoefficientVariation:
    case BuiltinId::StandardError:
    case BuiltinId::ZScore:
    case BuiltinId::Iqr:
    case BuiltinId::TrimMean:
    case BuiltinId::WinsorMean:
    case BuiltinId::Winsorized:
    case BuiltinId::Covariance:
    case BuiltinId::Correlation:
    case BuiltinId::SpearmanCorrelation:
    case BuiltinId::PercentRank:
    case BuiltinId::Dimensions:
    case BuiltinId::ArrayRank:
    case BuiltinId::ArrayGet:
    case BuiltinId::Reshape:
    case BuiltinId::Identity:
    case BuiltinId::Zeros:
    case BuiltinId::MatrixGet:
    case BuiltinId::Trace:
    case BuiltinId::Rows:
    case BuiltinId::Cols:
    case BuiltinId::Diag:
    case BuiltinId::VectorAdd:
    case BuiltinId::VectorSubtract:
    case BuiltinId::VectorScale:
    case BuiltinId::VectorDot:
    case BuiltinId::VectorCross:
    case BuiltinId::VectorNorm:
    case BuiltinId::VectorManhattan:
    case BuiltinId::VectorEuclidean:
    case BuiltinId::VectorNormalize:
    case BuiltinId::VectorProject:
    case BuiltinId::VectorAngle:
    case BuiltinId::VectorReflect:
    case BuiltinId::VectorReflectAxis:
    case BuiltinId::VectorSum:
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify:
    case BuiltinId::Map:
    case BuiltinId::Range:
    case BuiltinId::Table:
    case BuiltinId::Expand:
    case BuiltinId::Factor:
    case BuiltinId::Collect:
    case BuiltinId::Solve:
    case BuiltinId::Set:
    case BuiltinId::SetDelayed:
    case BuiltinId::Less:
    case BuiltinId::LessEqual:
    case BuiltinId::Greater:
    case BuiltinId::GreaterEqual:
    case BuiltinId::Equal:
    case BuiltinId::NotEqual:
    case BuiltinId::LogicalAnd:
    case BuiltinId::Element:
    case BuiltinId::If:
    case BuiltinId::History:
    case BuiltinId::InputHistory:
    case BuiltinId::OutputHistory:
    case BuiltinId::Exit:
    case BuiltinId::Clear:
    case BuiltinId::Definitions:
    case BuiltinId::Undefine:
    case BuiltinId::AngleMode:
        return {};
    }

    static_cast<void>(mathematics);
    return {};
}

} // namespace

bool ValueFacts::isNumeric() const noexcept {
    return domain != NumericDomain::Unknown;
}

bool ValueFacts::isProvablyReal() const noexcept {
    return domain == NumericDomain::Integer
        || domain == NumericDomain::Rational
        || domain == NumericDomain::Real;
}

bool ValueFacts::isProvablyComplex() const noexcept {
    return isNumeric();
}

bool ValueFacts::isProvablyNonNegativeReal() const noexcept {
    return isProvablyReal()
        && (sign == RealSign::Zero
            || sign == RealSign::Positive
            || sign == RealSign::NonNegative);
}

bool ValueFacts::isProvablyNegativeReal() const noexcept {
    return isProvablyReal() && sign == RealSign::Negative;
}

namespace {

ValueFacts inferValueFactsImpl(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AssumptionSet* assumptions) {
    struct WorkItem final {
        Expr expression;
        bool expanded = false;
    };

    std::vector<WorkItem> stack;
    stack.push_back(WorkItem{expression, false});
    std::unordered_map<const void*, ValueFacts> facts;

    while (!stack.empty()) {
        WorkItem current = std::move(stack.back());
        stack.pop_back();

        if (facts.contains(current.expression.identity()))
            continue;

        if (current.expression.isCall() && !current.expanded) {
            stack.push_back(WorkItem{current.expression, true});
            const auto& arguments = current.expression.asCall().arguments;
            for (auto iterator = arguments.rbegin(); iterator != arguments.rend(); ++iterator) {
                if (!facts.contains(iterator->identity()))
                    stack.push_back(WorkItem{*iterator, false});
            }
            continue;
        }

        ValueFacts result;
        switch (current.expression.kind()) {
        case expression::ExprKind::Number:
            result = factsForNumber(current.expression.asNumber());
            break;
        case expression::ExprKind::DecimalApproximation:
            result = ValueFacts{NumericDomain::Real, RealSign::Unknown, false, false};
            break;
        case expression::ExprKind::ComplexDecimalApproximation:
            result = ValueFacts{NumericDomain::Complex, RealSign::Unknown, false, false};
            break;
        case expression::ExprKind::Symbol:
            result = factsForSymbol(current.expression.asSymbol(), mathematics);
            break;
        case expression::ExprKind::Call:
            result = inferCallFacts(current.expression.asCall(), facts, builtins, mathematics);
            break;
        case expression::ExprKind::Boolean:
        case expression::ExprKind::String:
        case expression::ExprKind::Array:
        case expression::ExprKind::List:
        case expression::ExprKind::SolutionSet:
            result = {};
            break;
        }

        result = applyAssumptions(current.expression, result, assumptions);
        facts.emplace(current.expression.identity(), result);
    }

    return facts.at(expression.identity());
}

} // namespace

ValueFacts inferValueFacts(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    return inferValueFactsImpl(expression, builtins, mathematics, nullptr);
}

ValueFacts inferValueFacts(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AssumptionSet& assumptions) {
    return inferValueFactsImpl(expression, builtins, mathematics, &assumptions);
}

} // namespace mmcal::mathematics
