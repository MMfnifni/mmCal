// 安全な標準式簡約
#include "simplifier.hpp"
#include "expression/exact_value.hpp"
#include "expression/array_utils.hpp"

#include "expression_ordering.hpp"

#include "error/error_message.hpp"
#include "mathematics/definedness.hpp"
#include "mathematics/exact_algebra.hpp"
#include "mathematics/exact_hyperbolic.hpp"
#include "mathematics/exact_roots.hpp"
#include "mathematics/exact_transcendental.hpp"
#include "mathematics/exact_trigonometry.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mmcal::simplification {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::RelationKind;
using mathematics::TruthValue;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;


[[nodiscard]] bool provablyDefined(
    const Expr& expression,
    const SimplificationContext& context) {
    if (context.assumeExpressionsDefined)
        return true;

    const auto conditions = mathematics::expressionDomainConditions(
        expression, context.builtins, context.mathematics);
    if (!conditions)
        return false;

    const mathematics::KnowledgeContext knowledge = context.knowledge();
    for (const mathematics::Predicate& condition : conditions->predicates())
        if (knowledge.prove(condition) != TruthValue::True)
            return false;
    return true;
}

[[nodiscard]] std::optional<mathematics::Predicate> conditionPredicate(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    std::optional<RelationKind> kind;
    switch (definition->id) {
    case BuiltinId::Equal: kind = RelationKind::Equal; break;
    case BuiltinId::NotEqual: kind = RelationKind::NotEqual; break;
    case BuiltinId::Less: kind = RelationKind::Less; break;
    case BuiltinId::LessEqual: kind = RelationKind::LessEqual; break;
    case BuiltinId::Greater: kind = RelationKind::Greater; break;
    case BuiltinId::GreaterEqual: kind = RelationKind::GreaterEqual; break;
    default: return std::nullopt;
    }
    return mathematics::relation(
        *kind,
        expression.asCall().arguments[0],
        expression.asCall().arguments[1]);
}

[[nodiscard]] Expr integerExpr(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return builtins.isCallTo(expression, id);
}

[[nodiscard]] bool isExactReal(const Expr& expression, std::int64_t value) {
    return expression.isNumber()
        && expression.asNumber().isReal()
        && expression.asNumber().asReal() == RealNumber{BigInt{value}};
}

[[nodiscard]] std::optional<BigInt> positiveExactInteger(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    const RealNumber& real = expression.asNumber().asReal();
    if (!real.isInteger() || real.isNegative() || real.isZero())
        return std::nullopt;
    return real.asInteger();
}

[[nodiscard]] bool isMathematicalConstant(
    const Expr& expression,
    const mathematics::MathRegistry& mathematics,
    mathematics::ConstantId id) {
    if (!expression.isSymbol())
        return false;
    const auto* definition = mathematics.findConstant(id);
    return definition && expression.asSymbol().sameIdentity(definition->symbol);
}

void sortCanonical(std::vector<Expr>& expressions) {
    // total-order keyは式ごとに一度だけ生成し、n log n回の再帰serializeを避ける。
    std::vector<std::pair<std::string, Expr>> keyed;
    keyed.reserve(expressions.size());
    for (Expr& expression : expressions)
        keyed.emplace_back(expressionOrderKey(expression), std::move(expression));

    std::sort(keyed.begin(), keyed.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });
    for (std::size_t i = 0; i < keyed.size(); ++i)
        expressions[i] = std::move(keyed[i].second);
}

[[nodiscard]] bool containsInfinity(const Expr& expression) {
    if (expression.isSymbol()) {
        const std::string_view name = expression.asSymbol().view();
        return name == "Infinity" || name == "ComplexInfinity"
            || name == "Indeterminate";
    }
    if (expression.isCall()) {
        for (const Expr& argument : expression.asCall().arguments)
            if (containsInfinity(argument))
                return true;
    }
    else if (expression.isArray()) {
        const auto& array = expression.asArray();
        for (std::size_t i = 0; i < array.size(); ++i)
            if (containsInfinity(array.element(i)))
                return true;
    }
    else if (expression.isList()) {
        for (const Expr& element : expression.asList().elements)
            if (containsInfinity(element))
                return true;
    }
    return false;
}

[[nodiscard]] std::optional<BuiltinId> principalInverseOf(BuiltinId function) noexcept {
    switch (function) {
    case BuiltinId::Sin: return BuiltinId::Asin;
    case BuiltinId::Cos: return BuiltinId::Acos;
    case BuiltinId::Tan: return BuiltinId::Atan;
    case BuiltinId::Sinh: return BuiltinId::Asinh;
    case BuiltinId::Cosh: return BuiltinId::Acosh;
    case BuiltinId::Tanh: return BuiltinId::Atanh;
    default: return std::nullopt;
    }
}

[[nodiscard]] std::optional<Expr> simplifyPrincipalInverseComposition(
    BuiltinId outer,
    const Expr& argument,
    const SimplificationContext& context) {
    const auto inverse = principalInverseOf(outer);
    if (!inverse || !isHead(argument, context.builtins, *inverse)
        || argument.asCall().arguments.size() != 1)
        return std::nullopt;

    const Expr& inner = argument.asCall().arguments.front();
    // f[f^-1[z]] の向きだけを縮約する。逆向きは周期性・branch cutのため一般には成立しない。
    // extended/special infinityはinverseの有限平面上の恒等式の対象外なので，従来の未評価を保つ。
    // atan/atanh のbranch point等で未定義な点を消さないよう，inverse call全体のdefinednessも要求する。
    if (containsInfinity(inner) || !provablyDefined(argument, context))
        return std::nullopt;
    return inner;
}

[[nodiscard]] std::optional<Expr> predefinedValue(
    const SimplificationContext& context,
    symbols::PredefinedSymbolId id) {
    if (!context.predefinedSymbols)
        return std::nullopt;
    const auto* definition = context.predefinedSymbols->find(id);
    if (!definition)
        return std::nullopt;
    return Expr{definition->symbol};
}

[[nodiscard]] bool isPositiveInfinity(const Expr& expression) noexcept {
    return expression.isSymbol() && expression.asSymbol().view() == "Infinity";
}

template <typename Visitor>
void visitSignedAddTerms(
    const Expr& expression,
    bool negative,
    const evaluation::BuiltinRegistry& builtins,
    Visitor&& visitor) {
    if (isHead(expression, builtins, BuiltinId::Add)) {
        for (const Expr& nested : expression.asCall().arguments)
            visitSignedAddTerms(nested, negative, builtins, visitor);
        return;
    }

    if (isHead(expression, builtins, BuiltinId::Subtract)
        && expression.asCall().arguments.size() == 2) {
        visitSignedAddTerms(expression.asCall().arguments[0], negative, builtins, visitor);
        visitSignedAddTerms(expression.asCall().arguments[1], !negative, builtins, visitor);
        return;
    }

    if (isHead(expression, builtins, BuiltinId::Negate)
        && expression.asCall().arguments.size() == 1) {
        visitSignedAddTerms(expression.asCall().arguments.front(), !negative, builtins, visitor);
        return;
    }

    visitor(expression, negative);
}

[[nodiscard]] Expr buildProduct(
    std::vector<Expr> factors,
    const evaluation::BuiltinRegistry& builtins) {
    if (factors.empty())
        return integerExpr(1);
    if (factors.size() == 1)
        return factors.front();
    sortCanonical(factors);
    return Expr::call(builtins.symbol(BuiltinId::Multiply), std::move(factors));
}

struct LinearTerm final {
    Rational coefficient;
    Expr atom;
};

[[nodiscard]] LinearTerm extractLinearTerm(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (isHead(expression, builtins, BuiltinId::Negate)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() == 1) {
            LinearTerm inner = extractLinearTerm(arguments.front(), builtins);
            inner.coefficient = -inner.coefficient;
            return inner;
        }
    }

    if (isHead(expression, builtins, BuiltinId::Divide)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() == 2) {
            const auto denominator = expression::exact::realRational(arguments[1]);
            if (denominator && !denominator->isZero()) {
                LinearTerm numerator = extractLinearTerm(arguments[0], builtins);
                numerator.coefficient /= *denominator;
                return numerator;
            }
        }
    }

    if (isHead(expression, builtins, BuiltinId::Multiply)) {
        Rational coefficient{BigInt{1}};
        std::vector<Expr> factors;
        for (const Expr& factor : expression.asCall().arguments) {
            if (const auto numeric = expression::exact::realRational(factor)) {
                coefficient *= *numeric;
                continue;
            }
            factors.push_back(factor);
        }

        if (!factors.empty())
            return LinearTerm{std::move(coefficient), buildProduct(std::move(factors), builtins)};
    }

    return LinearTerm{Rational{BigInt{1}}, expression};
}

struct TrigSquare final {
    BuiltinId function = BuiltinId::Sin;
    Expr argument;
};

[[nodiscard]] std::optional<TrigSquare> trigSquare(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Power))
        return std::nullopt;
    const auto& powerArguments = expression.asCall().arguments;
    if (powerArguments.size() != 2 || !isExactReal(powerArguments[1], 2)
        || !powerArguments[0].isCall()
        || powerArguments[0].asCall().arguments.size() != 1)
        return std::nullopt;

    if (isHead(powerArguments[0], builtins, BuiltinId::Sin))
        return TrigSquare{BuiltinId::Sin, powerArguments[0].asCall().arguments.front()};
    if (isHead(powerArguments[0], builtins, BuiltinId::Cos))
        return TrigSquare{BuiltinId::Cos, powerArguments[0].asCall().arguments.front()};
    return std::nullopt;
}

[[nodiscard]] Expr canonicalAdd(
    const std::vector<Expr>& arguments,
    const SimplificationContext& context) {
    const auto& builtins = context.builtins;
    Number numericSum{BigInt{0}};
    std::vector<Expr> terms;
    terms.reserve(arguments.size());
    for (const Expr& argument : arguments) {
        visitSignedAddTerms(argument, false, builtins, [&](const Expr& term, bool negative) {
            if (term.isNumber()) {
                numericSum += negative ? -term.asNumber() : term.asNumber();
                return;
            }
            if (negative)
                terms.push_back(Expr::call(builtins.symbol(BuiltinId::Negate), {term}));
            else
                terms.push_back(term);
        });
    }

    struct Group final {
        Expr atom;
        Rational coefficient;
    };
    std::vector<Group> groups;
    groups.reserve(terms.size());
    std::unordered_map<std::string, std::size_t> groupByKey;
    groupByKey.reserve(terms.size());

    for (const Expr& term : terms) {
        LinearTerm linear = extractLinearTerm(term, builtins);
        // Infinityは現在extended-real sentinelであり，通常symbolのように
        // c*x + d*xとして係数相殺すると Infinity-Infinity -> 0 を捏造する。
        // 拡張実数算術を完全実装するまではatomic termとして保持する。
        if (containsInfinity(linear.atom)) {
            groups.push_back(Group{std::move(linear.atom), std::move(linear.coefficient)});
            continue;
        }
        const std::string key = expressionOrderKey(linear.atom);
        auto found = groupByKey.find(key);
        if (found != groupByKey.end() && groups[found->second].atom == linear.atom) {
            groups[found->second].coefficient += linear.coefficient;
            continue;
        }

        // expressionOrderKeyはcanonical ordering用で数学的identityそのものではないため、万一文字列keyが衝突した場合はExpr equalityで安全にfallbackする。
        if (found != groupByKey.end()) {
            const auto iterator = std::find_if(
                groups.begin(), groups.end(),
                [&](const Group& group) { return group.atom == linear.atom; });
            if (iterator != groups.end()) {
                iterator->coefficient += linear.coefficient;
                continue;
            }
        }

        const std::size_t index = groups.size();
        groups.push_back(Group{std::move(linear.atom), std::move(linear.coefficient)});
        groupByKey.emplace(std::move(key), index);
    }

    // Pythagorean identity は実数角に限らず複素引数でも恒等的に成立する。
    // mmCalの裸の角度は現在のAngleSemantics（既定Radian）で解釈されるが、sin/cosの双方に同じ変換が適用されるため sin[u]^2 + cos[u]^2 = 1 は単位指定の有無に依存しない。
    // 同じ有理係数 c が掛かった対も c*(sin^2+cos^2) -> c としてまとめる。
    for (std::size_t i = 0; i < groups.size(); ++i) {
        if (groups[i].coefficient.isZero())
            continue;
        const auto lhs = trigSquare(groups[i].atom, builtins);
        if (!lhs || lhs->function != BuiltinId::Sin
            || !provablyDefined(lhs->argument, context))
            continue;
        for (std::size_t j = 0; j < groups.size(); ++j) {
            if (groups[j].coefficient != groups[i].coefficient)
                continue;
            const auto rhs = trigSquare(groups[j].atom, builtins);
            if (!rhs || rhs->function != BuiltinId::Cos || !(rhs->argument == lhs->argument))
                continue;
            numericSum += Number{groups[i].coefficient};
            groups[i].coefficient = Rational{BigInt{0}};
            groups[j].coefficient = Rational{BigInt{0}};
            break;
        }
    }

    // 加法の項順は係数ではなくatomで決める。1/3*x^3のように係数を
    // quotient normal formへ移してもx, x^2, x^3の順序が揺れない。
    std::sort(groups.begin(), groups.end(), [](const Group& lhs, const Group& rhs) {
        return ExpressionLess{}(lhs.atom, rhs.atom);
    });

    std::vector<Expr> result;
    if (!numericSum.isZero())
        result.emplace_back(std::move(numericSum));

    for (Group& group : groups) {
        if (group.coefficient.isZero()) {
            if (!provablyDefined(group.atom, context))
                result.push_back(Expr::call(
                    builtins.symbol(BuiltinId::Subtract), {group.atom, group.atom}));
            continue;
        }
        result.push_back(mathematics::scaleExactExpression(
            group.coefficient, group.atom, builtins));
    }

    if (result.empty())
        return integerExpr(0);
    if (result.size() == 1)
        return result.front();
    if (result.size() == 2
        && isHead(result[1], builtins, BuiltinId::Negate)
        && result[1].asCall().arguments.size() == 1)
        return Expr::call(
            builtins.symbol(BuiltinId::Subtract),
            {result[0], result[1].asCall().arguments.front()});
    return Expr::call(builtins.symbol(BuiltinId::Add), std::move(result));
}

[[nodiscard]] std::vector<Expr> groupRepeatedFactors(
    std::vector<Expr> factors,
    const evaluation::BuiltinRegistry& builtins) {
    sortCanonical(factors);
    std::vector<Expr> grouped;
    grouped.reserve(factors.size());
    for (std::size_t i = 0; i < factors.size();) {
        std::size_t j = i + 1;
        while (j < factors.size() && factors[j] == factors[i])
            ++j;
        const std::size_t count = j - i;
        if (count == 1)
            grouped.push_back(factors[i]);
        else
            grouped.push_back(Expr::call(
                builtins.symbol(BuiltinId::Power),
                {factors[i], Expr{Number{BigInt::parse(std::to_string(count))}}}));
        i = j;
    }
    return grouped;
}

struct ProductParts final {
    Number coefficient{BigInt{1}};
    std::vector<Expr> numerator;
    std::vector<Expr> denominator;
};

void collectProductParts(
    const Expr& expression,
    bool reciprocal,
    ProductParts& parts,
    const evaluation::BuiltinRegistry& builtins) {
    if (isHead(expression, builtins, BuiltinId::Negate)
        && expression.asCall().arguments.size() == 1) {
        parts.coefficient *= Number{BigInt{-1}};
        collectProductParts(
            expression.asCall().arguments.front(), reciprocal, parts, builtins);
        return;
    }

    if (isHead(expression, builtins, BuiltinId::Multiply)) {
        for (const Expr& factor : expression.asCall().arguments)
            collectProductParts(factor, reciprocal, parts, builtins);
        return;
    }

    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        // (a/b)*c は安全に ac/b へ平坦化できるが、a/(b/c) を ac/b へすると
        // c=0 という元式のholeを失う。逆数側で遭遇したDivideはatomic denominator
        // として保持し、definednessを変えない範囲だけをnormal form化する。
        if (reciprocal) {
            parts.denominator.push_back(expression);
            return;
        }
        const auto& arguments = expression.asCall().arguments;
        collectProductParts(arguments[0], false, parts, builtins);
        collectProductParts(arguments[1], true, parts, builtins);
        return;
    }

    if (expression.isNumber()) {
        if (reciprocal) {
            if (expression.asNumber().isZero())
                error::throwCalcError(error::CalcErrorType::Domain, "Division by zero");
            parts.coefficient /= expression.asNumber();
        }
        else
            parts.coefficient *= expression.asNumber();
        return;
    }

    (reciprocal ? parts.denominator : parts.numerator).push_back(expression);
}

[[nodiscard]] Expr productFromFactors(
    std::vector<Expr> factors,
    const evaluation::BuiltinRegistry& builtins) {
    factors = groupRepeatedFactors(std::move(factors), builtins);
    if (factors.empty())
        return integerExpr(1);
    if (factors.size() == 1)
        return factors.front();
    return Expr::call(builtins.symbol(BuiltinId::Multiply), std::move(factors));
}

[[nodiscard]] Expr buildProductNormalForm(
    ProductParts parts,
    const evaluation::BuiltinRegistry& builtins) {
    // denominatorを含まない0積は従来どおり0へ畳み込む。
    // 0*(1/x) のように明示的な分母を含む場合だけholeを失わない形を保持する。
    const bool hasInfinityFactor = std::any_of(
        parts.numerator.begin(), parts.numerator.end(), containsInfinity);
    if (parts.coefficient.isZero() && parts.denominator.empty() && !hasInfinityFactor)
        return integerExpr(0);

    // exact real係数は整数の分子・分母へ分解し、記号分母と同じDivideへ集約する。
    // これにより (a/2)*b と a*b/2 が同一構造になる。
    bool negative = false;
    if (parts.coefficient.isReal()) {
        Rational coefficient = parts.coefficient.asReal().toRational();
        negative = coefficient.numerator().isNegative();
        const BigInt numerator = coefficient.numerator().abs();
        const BigInt denominator = coefficient.denominator();
        if (!(numerator == BigInt{1}) || parts.numerator.empty())
            parts.numerator.emplace_back(Number{numerator});
        if (!(denominator == BigInt{1}))
            parts.denominator.emplace_back(Number{denominator});
    }
    else if (parts.coefficient.realPart().isZero()) {
        // 純虚数の有理係数 qI は q と I に分離し、real係数と同じ分母へ集約する。
        // これにより (I/2)*Pi と I*Pi/2 が同一構造になる。
        Rational coefficient = parts.coefficient.imaginaryPart().toRational();
        negative = coefficient.numerator().isNegative();
        const BigInt numerator = coefficient.numerator().abs();
        const BigInt denominator = coefficient.denominator();
        if (!(numerator == BigInt{1}))
            parts.numerator.emplace_back(Number{numerator});
        parts.numerator.emplace_back(Number::complex(RealNumber{}, RealNumber{BigInt{1}}));
        if (!(denominator == BigInt{1}))
            parts.denominator.emplace_back(Number{denominator});
    }
    else if (!(parts.coefficient == Number{BigInt{1}})) {
        parts.numerator.emplace_back(std::move(parts.coefficient));
    }

    Expr numerator = productFromFactors(std::move(parts.numerator), builtins);
    Expr result = std::move(numerator);
    if (!parts.denominator.empty()) {
        Expr denominator = productFromFactors(std::move(parts.denominator), builtins);
        result = Expr::call(
            builtins.symbol(BuiltinId::Divide),
            {std::move(result), std::move(denominator)});
    }

    // 符号はquotient全体の外へ出す。Addの線形項抽出が
    // -F + F を構造的に相殺でき、分子内Negateという別形を作らない。
    if (negative) {
        if (result.isNumber())
            return Expr{-result.asNumber()};
        return Expr::call(builtins.symbol(BuiltinId::Negate), {std::move(result)});
    }
    return result;
}

[[nodiscard]] Expr canonicalMultiply(
    const std::vector<Expr>& arguments,
    const evaluation::BuiltinRegistry& builtins) {
    ProductParts parts;
    for (const Expr& argument : arguments)
        collectProductParts(argument, false, parts, builtins);
    return buildProductNormalForm(std::move(parts), builtins);
}

[[nodiscard]] Expr canonicalDivide(
    const Expr& numerator,
    const Expr& denominator,
    const evaluation::BuiltinRegistry& builtins) {
    ProductParts parts;
    collectProductParts(numerator, false, parts, builtins);
    collectProductParts(denominator, true, parts, builtins);
    return buildProductNormalForm(std::move(parts), builtins);
}


[[nodiscard]] const std::vector<Expr>* hypergeometric1F1Arguments(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Hypergeometric1F1)
        || expression.asCall().arguments.size() != 3)
        return nullptr;
    return &expression.asCall().arguments;
}

[[nodiscard]] const std::vector<Expr>* hypergeometric2F1Arguments(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Hypergeometric2F1)
        || expression.asCall().arguments.size() != 4)
        return nullptr;
    return &expression.asCall().arguments;
}

struct PositiveIntegerPower final {
    Expr base;
    std::uint64_t exponent = 1;
};

[[nodiscard]] std::optional<PositiveIntegerPower> positiveIntegerPower(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Power))
        return PositiveIntegerPower{expression, 1};
    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() != 2)
        return std::nullopt;
    const auto exponent = expression::exact::realRational(arguments[1]);
    if (!exponent || !exponent->isInteger() || !exponent->numerator().isPositive())
        return std::nullopt;
    const auto count = numeric::tryToUint64(exponent->numerator());
    if (!count)
        return std::nullopt;
    return PositiveIntegerPower{arguments[0], *count};
}

[[nodiscard]] bool samePositiveMonomial(
    std::span<const Expr> factors,
    const Expr& expected,
    const evaluation::BuiltinRegistry& builtins) {
    const auto target = positiveIntegerPower(expected, builtins);
    if (!target)
        return false;
    std::uint64_t total = 0;
    for (const Expr& factor : factors) {
        const auto power = positiveIntegerPower(factor, builtins);
        if (!power || !(power->base == target->base)
            || power->exponent > std::numeric_limits<std::uint64_t>::max() - total)
            return false;
        total += power->exponent;
    }
    return total == target->exponent;
}

[[nodiscard]] bool productIsExactly(
    const Expr& expression,
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins) {
    if (expression == lhs && isExactReal(rhs, 1))
        return true;
    if (!isHead(expression, builtins, BuiltinId::Multiply))
        return false;

    const auto& factors = expression.asCall().arguments;
    std::vector<Expr> remaining;
    bool foundRhs = false;
    remaining.reserve(factors.size());
    for (const Expr& factor : factors) {
        if (!foundRhs && factor == rhs) {
            foundRhs = true;
            continue;
        }
        remaining.push_back(factor);
    }
    return foundRhs && !remaining.empty()
        && samePositiveMonomial(remaining, lhs, builtins);
}

[[nodiscard]] std::optional<Expr> simplifyHypergeometricContiguous(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Add)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const auto tryOrientation = [&](const Expr& first, const Expr& second) -> std::optional<Expr> {
        LinearTerm baseTerm = extractLinearTerm(first, builtins);
        if (baseTerm.coefficient != Rational{BigInt{1}})
            return std::nullopt;
        const auto* base = hypergeometric1F1Arguments(baseTerm.atom, builtins);
        if (!base)
            return std::nullopt;

        const auto a = expression::exact::realRational((*base)[0]);
        const auto b = expression::exact::realRational((*base)[1]);
        if (!a || !b || a->isZero() || *b != *a + Rational{BigInt{1}})
            return std::nullopt;
        // 積分器が生成するa=1/n>0の領域だけに限定し、parameter poleを跨ぐ一般変形にしない。
        if (a->numerator().isNegative())
            return std::nullopt;

        const Expr& z = (*base)[2];
        LinearTerm zTerm = extractLinearTerm(z, builtins);
        LinearTerm correction = extractLinearTerm(second, builtins);
        const Rational expectedCoefficient = zTerm.coefficient / *b;
        if (correction.coefficient != expectedCoefficient)
            return std::nullopt;

        const auto* shifted = hypergeometric1F1Arguments(
            isHead(correction.atom, builtins, BuiltinId::Multiply)
                ? [&]() -> const Expr& {
                    for (const Expr& factor : correction.atom.asCall().arguments)
                        if (hypergeometric1F1Arguments(factor, builtins))
                            return factor;
                    return correction.atom;
                }()
                : correction.atom,
            builtins);
        if (!shifted)
            return std::nullopt;

        const auto shiftedA = expression::exact::realRational((*shifted)[0]);
        const auto shiftedB = expression::exact::realRational((*shifted)[1]);
        if (!shiftedA || !shiftedB
            || *shiftedA != *a + Rational{BigInt{1}}
            || *shiftedB != *b + Rational{BigInt{1}}
            || !((*shifted)[2] == z))
            return std::nullopt;

        const Expr shiftedCall = Expr::call(
            builtins.symbol(BuiltinId::Hypergeometric1F1), *shifted);
        if (!productIsExactly(correction.atom, zTerm.atom, shiftedCall, builtins))
            return std::nullopt;

        // M(a,a+1,z) + z/(a+1) M(a+1,a+2,z) = exp(z)。
        return Expr::call(builtins.symbol(BuiltinId::Exp), {z});
    };

    const auto& terms = expression.asCall().arguments;
    if (auto result = tryOrientation(terms[0], terms[1]))
        return result;
    return tryOrientation(terms[1], terms[0]);
}


[[nodiscard]] std::optional<Expr> simplifyHypergeometric2F1Contiguous(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Add)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const auto tryOrientation = [&](const Expr& first, const Expr& second) -> std::optional<Expr> {
        LinearTerm baseTerm = extractLinearTerm(first, builtins);
        if (baseTerm.coefficient != Rational{BigInt{1}})
            return std::nullopt;
        const auto* base = hypergeometric2F1Arguments(baseTerm.atom, builtins);
        if (!base)
            return std::nullopt;
        const auto a = expression::exact::realRational((*base)[0]);
        const auto b = expression::exact::realRational((*base)[1]);
        const auto c = expression::exact::realRational((*base)[2]);
        if (!a || !b || !c || b->numerator().isNegative() || b->isZero()
            || *c != *b + Rational{BigInt{1}})
            return std::nullopt;

        const Expr& z = (*base)[3];
        LinearTerm zTerm = extractLinearTerm(z, builtins);
        LinearTerm correction = extractLinearTerm(second, builtins);
        const Rational expectedCoefficient = zTerm.coefficient
            * *a / (*b + Rational{BigInt{1}});
        if (correction.coefficient != expectedCoefficient)
            return std::nullopt;

        const Expr* shiftedAtom = &correction.atom;
        if (isHead(correction.atom, builtins, BuiltinId::Multiply)) {
            for (const Expr& factor : correction.atom.asCall().arguments)
                if (hypergeometric2F1Arguments(factor, builtins)) {
                    shiftedAtom = &factor;
                    break;
                }
        }
        const auto* shifted = hypergeometric2F1Arguments(*shiftedAtom, builtins);
        if (!shifted)
            return std::nullopt;
        const auto shiftedA = expression::exact::realRational((*shifted)[0]);
        const auto shiftedB = expression::exact::realRational((*shifted)[1]);
        const auto shiftedC = expression::exact::realRational((*shifted)[2]);
        if (!shiftedA || !shiftedB || !shiftedC
            || *shiftedA != *a + Rational{BigInt{1}}
            || *shiftedB != *b + Rational{BigInt{1}}
            || *shiftedC != *c + Rational{BigInt{1}}
            || !((*shifted)[3] == z))
            return std::nullopt;

        const Expr shiftedCall = Expr::call(
            builtins.symbol(BuiltinId::Hypergeometric2F1), *shifted);
        if (!productIsExactly(correction.atom, zTerm.atom, shiftedCall, builtins))
            return std::nullopt;

        // F(a,b;b+1;z)+a z/(b+1) F(a+1,b+1;b+2;z)=(1-z)^(-a)。
        // binomial-power積分が生成するb>0の限定形だけをproof Knowledgeとして使う。
        return Expr::call(builtins.symbol(BuiltinId::Power), {
            Expr::call(builtins.symbol(BuiltinId::Subtract), {integerExpr(1), z}),
            Expr{Number{RealNumber{-*a}}}});
    };

    const auto& terms = expression.asCall().arguments;
    if (auto result = tryOrientation(terms[0], terms[1]))
        return result;
    return tryOrientation(terms[1], terms[0]);
}

[[nodiscard]] TruthValue proveNonZero(
    const Expr& expression,
    const mathematics::KnowledgeContext& knowledge) {
    return knowledge.prove(mathematics::relation(
        RelationKind::NotEqual, expression, integerExpr(0)));
}

[[nodiscard]] Expr localRewrite(
    const Expr& expression,
    const SimplificationContext& context) {
    if (!expression.isCall())
        return expression;

    const auto* definition = context.builtins.find(expression.asCall().head);
    if (!definition)
        return expression;

    const auto& arguments = expression.asCall().arguments;
    const mathematics::KnowledgeContext knowledge = context.knowledge();

    switch (definition->id) {
    case BuiltinId::Add: {
        Expr canonical = canonicalAdd(arguments, context);
        if (auto hypergeometric = simplifyHypergeometricContiguous(canonical, context.builtins))
            return *hypergeometric;
        if (auto hypergeometric = simplifyHypergeometric2F1Contiguous(canonical, context.builtins))
            return *hypergeometric;
        return canonical;
    }

    case BuiltinId::Multiply:
        if (std::any_of(arguments.begin(), arguments.end(), [](const Expr& argument) {
                return argument.isNumber() && argument.asNumber().isZero();
            })) {
            for (const Expr& argument : arguments)
                if (!(argument.isNumber() && argument.asNumber().isZero())
                    && !provablyDefined(argument, context))
                    return expression;
        }
        return canonicalMultiply(arguments, context.builtins);

    case BuiltinId::Subtract: {
        if (arguments.size() != 2)
            return expression;
        if (arguments[0].isNumber() && arguments[1].isNumber())
            return Expr{arguments[0].asNumber() - arguments[1].asNumber()};
        if (isExactReal(arguments[1], 0))
            return arguments[0];
        if (arguments[0] == arguments[1] && !containsInfinity(arguments[0])
            && provablyDefined(arguments[0], context))
            return integerExpr(0);
        if (isExactReal(arguments[0], 0))
            return Expr::call(
                context.builtins.symbol(BuiltinId::Negate), {arguments[1]});
        if (arguments[1].isNumber() && arguments[1].asNumber().isReal()
            && arguments[1].asNumber().asReal().isNegative())
            return Expr::call(
                context.builtins.symbol(BuiltinId::Add),
                {arguments[0], Expr{-arguments[1].asNumber()}});
        // (a-b)-c でb,cがexact Rationalなら定数項だけを先に畳む。
        // subtraction自体のdomainは変わらないためsymbolicなaにも安全に適用できる。
        if (isHead(arguments[0], context.builtins, BuiltinId::Subtract)
            && arguments[0].asCall().arguments.size() == 2) {
            const auto& inner = arguments[0].asCall().arguments;
            const auto innerConstant = expression::exact::realRational(inner[1]);
            const auto outerConstant = expression::exact::realRational(arguments[1]);
            if (innerConstant && outerConstant)
                return Expr::call(context.builtins.symbol(BuiltinId::Subtract), {
                    inner[0], Expr{Number{*innerConstant + *outerConstant}}});
        }
        if (isHead(arguments[1], context.builtins, BuiltinId::Negate)
            && arguments[1].asCall().arguments.size() == 1)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Add),
                {arguments[0], arguments[1].asCall().arguments.front()});
        // builtin評価で(a+b)-aが後から現れる場合だけ，同一のdefined termを安全に相殺する。
        // 一般のsubtractionを加法normal formへ並べ替えないため，積分結果等の既存表示順は維持する。
        if (isHead(arguments[0], context.builtins, BuiltinId::Add)
            && provablyDefined(arguments[1], context)) {
            std::vector<Expr> remaining = arguments[0].asCall().arguments;
            const auto match = std::find(remaining.begin(), remaining.end(), arguments[1]);
            if (match != remaining.end()) {
                remaining.erase(match);
                if (remaining.empty())
                    return integerExpr(0);
                return canonicalAdd(remaining, context);
            }
        }

        return expression;
    }

    case BuiltinId::Divide:
        if (arguments.size() != 2)
            return expression;
        if (arguments[1].isNumber() && arguments[1].asNumber().isZero())
            error::throwCalcError(error::CalcErrorType::Domain, "Division by zero");
        if (arguments[0].isNumber() && arguments[1].isNumber())
            return Expr{arguments[0].asNumber() / arguments[1].asNumber()};
        if (isExactReal(arguments[1], 1))
            return arguments[0];
        if ((arguments[0] == arguments[1] || isExactReal(arguments[0], 0))
            && provablyDefined(arguments[1], context)
            && proveNonZero(arguments[1], knowledge) == TruthValue::True)
            return arguments[0] == arguments[1] ? integerExpr(1) : integerExpr(0);

        // 分子積に分母そのものが含まれ、その因子が非零と証明できる場合だけ約分する。
        // 旧canonicalDivideはdomain hole保持のためsymbolic因子を一切cancelしなかったが、
        // Pi*x/PiのようにKnowledgeが非零を証明できるケースまで残していた。
        if (isHead(arguments[0], context.builtins, BuiltinId::Multiply)
            && provablyDefined(arguments[1], context)
            && proveNonZero(arguments[1], knowledge) == TruthValue::True) {
            std::vector<Expr> factors = arguments[0].asCall().arguments;
            const auto match = std::find(factors.begin(), factors.end(), arguments[1]);
            if (match != factors.end()) {
                factors.erase(match);
                if (factors.empty())
                    return integerExpr(1);
                if (factors.size() == 1)
                    return factors.front();
                return Expr::call(context.builtins.symbol(BuiltinId::Multiply), std::move(factors));
            }
        }

        // 正のexact Rationalの平方根を分母に持つ場合だけ共役化する。
        // c/sqrt[r] = c sqrt[r]/r (r>0) はprincipal branchでも安全で、atan[1/sqrt[3]] 等を既存のexact inverse-trig知識へ正規化できる。
        if (const auto numerator = expression::exact::realRational(arguments[0]); numerator
            && isHead(arguments[1], context.builtins, BuiltinId::Sqrt)
            && arguments[1].asCall().arguments.size() == 1) {
            if (const auto radicand = expression::exact::realRational(
                    arguments[1].asCall().arguments[0]);
                radicand && *radicand > Rational{BigInt{0}}) {
                return mathematics::scaleExactExpression(
                    *numerator / *radicand, arguments[1], context.builtins);
            }
        }

        // exactな有理分母なら、分子に既に付いている有理係数と先に約分する。
        // 例: (4 sqrt[2])/2 -> 2 sqrt[2]。branchやdomainには触れず、Rational係数だけを整理するため常に安全。
        if (const auto denominator = expression::exact::realRational(arguments[1]);
            denominator && !denominator->isZero()) {
            LinearTerm numerator = extractLinearTerm(arguments[0], context.builtins);
            if (!(numerator.coefficient == Rational{BigInt{1}})
                || !(numerator.atom == arguments[0])) {
                return mathematics::scaleExactExpression(
                    numerator.coefficient / *denominator,
                    numerator.atom,
                    context.builtins);
            }
        }
        return canonicalDivide(arguments[0], arguments[1], context.builtins);

    case BuiltinId::Negate:
        if (arguments.size() != 1)
            return expression;
        if (arguments[0].isNumber())
            return Expr{-arguments[0].asNumber()};
        if (isHead(arguments[0], context.builtins, BuiltinId::Negate)
            && arguments[0].asCall().arguments.size() == 1)
            return arguments[0].asCall().arguments.front();
        if (isHead(arguments[0], context.builtins, BuiltinId::If)
            && arguments[0].asCall().arguments.size() == 3) {
            const auto& branch = arguments[0].asCall().arguments;
            return Expr::call(context.builtins.symbol(BuiltinId::If), {
                branch[0],
                Expr::call(context.builtins.symbol(BuiltinId::Negate), {branch[1]}),
                Expr::call(context.builtins.symbol(BuiltinId::Negate), {branch[2]})});
        }
        if (isHead(arguments[0], context.builtins, BuiltinId::Cases)) {
            std::vector<Expr> branches;
            branches.reserve(arguments[0].asCall().arguments.size());
            for (const Expr& branchExpression : arguments[0].asCall().arguments) {
                if (!isHead(branchExpression, context.builtins, BuiltinId::CaseBranch)
                    || branchExpression.asCall().arguments.empty()
                    || branchExpression.asCall().arguments.size() > 2)
                    return expression;
                const auto& branch = branchExpression.asCall().arguments;
                std::vector<Expr> branchArguments{
                    Expr::call(context.builtins.symbol(BuiltinId::Negate), {branch[0]})};
                if (branch.size() == 2)
                    branchArguments.push_back(branch[1]);
                branches.push_back(Expr::call(
                    context.builtins.symbol(BuiltinId::CaseBranch),
                    std::move(branchArguments)));
            }
            return Expr::call(context.builtins.symbol(BuiltinId::Cases), std::move(branches));
        }
        return expression;

    case BuiltinId::Power:
        if (arguments.size() != 2)
            return expression;
        if (arguments[1].isSymbol()
            && arguments[1].asSymbol().view() == "ComplexInfinity") {
            if (const auto value = predefinedValue(
                    context, symbols::PredefinedSymbolId::Indeterminate))
                return *value;
        }
        if (isPositiveInfinity(arguments[1])) {
            const Expr magnitude = Expr::call(
                context.builtins.symbol(BuiltinId::Abs), {arguments[0]});
            if (knowledge.prove(mathematics::relation(
                    RelationKind::Equal, magnitude, integerExpr(1))) == TruthValue::True) {
                if (const auto value = predefinedValue(
                        context, symbols::PredefinedSymbolId::Indeterminate))
                    return *value;
            }
        }
        // E^z = exp[z] はprincipal branchに依存しない全域恒等式。
        // 指数函数の知識をPower/Expで二重化せず、canonicalなExp headへ寄せる。
        if (isMathematicalConstant(
                arguments[0], context.mathematics, mathematics::ConstantId::E))
            return Expr::call(context.builtins.symbol(BuiltinId::Exp), {arguments[1]});
        if (arguments[0].isNumber() && arguments[1].isNumber()
            && arguments[1].asNumber().isReal()
            && arguments[1].asNumber().asReal().isInteger()) {
            const BigInt& exponent = arguments[1].asNumber().asReal().asInteger();
            if (arguments[0].asNumber().isZero()) {
                if (exponent.isZero()) {
                    if (const auto value = predefinedValue(
                            context, symbols::PredefinedSymbolId::Indeterminate))
                        return *value;
                    return expression;
                }
                if (exponent.isNegative())
                    error::throwCalcError(error::CalcErrorType::Domain,
                        "Zero cannot be raised to a negative power");
            }
            if (const auto magnitude = numeric::tryToUint64(exponent.abs())) {
                Number value = numeric::integerPower(arguments[0].asNumber(), *magnitude);
                if (exponent.isNegative())
                    value = Number{BigInt{1}} / value;
                return Expr{std::move(value)};
            }
        }
        // (a^m)^n=a^(mn) はm,nが正のexact integerなら複素数全体で安全。
        // 一般複素指数には拡張せず，Powerのbranch semanticsを保持する。
        if (const auto outerExponent = positiveExactInteger(arguments[1]); outerExponent
            && isHead(arguments[0], context.builtins, BuiltinId::Power)
            && arguments[0].asCall().arguments.size() == 2) {
            const auto& innerArguments = arguments[0].asCall().arguments;
            if (const auto innerExponent = positiveExactInteger(innerArguments[1])) {
                const BigInt combined = *innerExponent * *outerExponent;
                return Expr::call(context.builtins.symbol(BuiltinId::Power), {
                    innerArguments[0], Expr{Number{combined}}});
            }
        }
        // (a/b)^n = a^n/b^n は正のexact integer n ならprincipal branchに依存せず安全。
        // sqrt等を含む有理函数の微分後に (u/sqrt[c])^2 をu^2/cへ落とせるようにする。
        if (const auto exponent = positiveExactInteger(arguments[1]); exponent
            && isHead(arguments[0], context.builtins, BuiltinId::Divide)
            && arguments[0].asCall().arguments.size() == 2) {
            const auto magnitude = numeric::tryToUint64(*exponent);
            if (magnitude && *magnitude <= 64) {
                const auto& quotient = arguments[0].asCall().arguments;
                return Expr::call(context.builtins.symbol(BuiltinId::Divide), {
                    Expr::call(context.builtins.symbol(BuiltinId::Power),
                        {quotient[0], arguments[1]}),
                    Expr::call(context.builtins.symbol(BuiltinId::Power),
                        {quotient[1], arguments[1]})});
            }
        }

        // (c*a)^n でnが正のexact integerなら、exact numeric係数cだけを外へ出す。
        // (ab)^z=a^z b^z を一般複素指数へ拡張せず、式サイズも増やさない限定形。
        // 例: (861(1-x)^64)^3 -> 638277381((1-x)^64)^3 -> 638277381(1-x)^192。
        if (const auto outerExponent = positiveExactInteger(arguments[1]); outerExponent
            && isHead(arguments[0], context.builtins, BuiltinId::Multiply)) {
            const auto& factors = arguments[0].asCall().arguments;
            Number coefficient{BigInt{1}};
            std::vector<Expr> symbolicFactors;
            symbolicFactors.reserve(factors.size());
            for (const Expr& factor : factors) {
                if (factor.isNumber())
                    coefficient *= factor.asNumber();
                else
                    symbolicFactors.push_back(factor);
            }

            if (!(coefficient == Number{BigInt{1}}) && !symbolicFactors.empty()) {
                Expr symbolicBase = productFromFactors(
                    std::move(symbolicFactors), context.builtins);
                Expr coefficientPower = Expr::call(
                    context.builtins.symbol(BuiltinId::Power),
                    {Expr{coefficient}, arguments[1]});
                Expr symbolicPower = Expr::call(
                    context.builtins.symbol(BuiltinId::Power),
                    {std::move(symbolicBase), arguments[1]});
                return Expr::call(
                    context.builtins.symbol(BuiltinId::Multiply),
                    {std::move(coefficientPower), std::move(symbolicPower)});
            }
        }
        if (isExactReal(arguments[1], 1))
            return arguments[0];
        if (isExactReal(arguments[1], 0)) {
            if (provablyDefined(arguments[0], context)
                && proveNonZero(arguments[0], knowledge) == TruthValue::True)
                return integerExpr(1);
            return expression;
        }
        if (const auto exponent = expression::exact::realRational(arguments[1]); exponent) {
            if (*exponent == Rational{BigInt{1}, BigInt{2}})
                return Expr::call(context.builtins.symbol(BuiltinId::Sqrt), {arguments[0]});
            // principal sqrtは定義上 w^2=z を満たすため、(sqrt[z])^2=z はbranch cutを跨いでも安全。
            // sqrt[z^2]とは異なり、こちらは条件なしで簡約できる。
            if (*exponent == Rational{BigInt{2}}
                && isHead(arguments[0], context.builtins, BuiltinId::Sqrt)
                && arguments[0].asCall().arguments.size() == 1)
                return arguments[0].asCall().arguments.front();
            // real cbrtは実軸全体で一価なcubeの逆函数なので、引数がRealと証明済みなら
            // (cbrt[z])^3=z。未証明の複素zへ拡張するとcbrtのdomain holeを消すため行わない。
            if (*exponent == Rational{BigInt{3}}
                && isHead(arguments[0], context.builtins, BuiltinId::Cbrt)
                && arguments[0].asCall().arguments.size() == 1
                && knowledge.facts(arguments[0].asCall().arguments.front()).isProvablyReal())
                return arguments[0].asCall().arguments.front();
        }
        return expression;

    case BuiltinId::Cbrt: {
        if (arguments.size() != 1)
            return expression;
        const Expr& input = arguments[0];
        if (input.isNumber()) {
            if (!input.asNumber().isReal())
                return expression;
            const RealNumber& real = input.asNumber().asReal();
            if (const auto root = mathematics::exactRealCubeRoot(real))
                return Expr{Number{*root}};
            if (real.isNegative())
                return Expr::call(
                    context.builtins.symbol(BuiltinId::Negate),
                    {Expr::call(
                        context.builtins.symbol(BuiltinId::Cbrt),
                        {Expr{Number{real.abs()}}})});
        }
        if (isHead(input, context.builtins, BuiltinId::Negate)
            && input.asCall().arguments.size() == 1
            && knowledge.facts(input.asCall().arguments.front()).isProvablyReal())
            return Expr::call(
                context.builtins.symbol(BuiltinId::Negate),
                {Expr::call(
                    context.builtins.symbol(BuiltinId::Cbrt),
                    {input.asCall().arguments.front()})});
        return expression;
    }

    case BuiltinId::Hypot: {
        if (arguments.size() != 2)
            return expression;
        if (!knowledge.facts(arguments[0]).isProvablyReal()
            || !knowledge.facts(arguments[1]).isProvablyReal())
            return expression;
        const Expr x2 = Expr::call(
            context.builtins.symbol(BuiltinId::Power), {arguments[0], integerExpr(2)});
        const Expr y2 = Expr::call(
            context.builtins.symbol(BuiltinId::Power), {arguments[1], integerExpr(2)});
        return Expr::call(
            context.builtins.symbol(BuiltinId::Sqrt),
            {Expr::call(context.builtins.symbol(BuiltinId::Add), {x2, y2})});
    }

    case BuiltinId::Cis: {
        if (arguments.size() != 1)
            return expression;
        const auto cosine = mathematics::simplifyExactTrig(
            mathematics::FunctionId::Cos, arguments[0], context.builtins,
            context.mathematics, context.angleSemantics);
        const auto sine = mathematics::simplifyExactTrig(
            mathematics::FunctionId::Sin, arguments[0], context.builtins,
            context.mathematics, context.angleSemantics);
        if (!cosine || !sine)
            return expression;
        const Expr imaginaryUnit{Number::complex(RealNumber{}, RealNumber{BigInt{1}})};
        return Expr::call(
            context.builtins.symbol(BuiltinId::Add),
            {*cosine, Expr::call(
                context.builtins.symbol(BuiltinId::Multiply),
                {imaginaryUnit, *sine})});
    }

    case BuiltinId::Sqrt: {
        if (arguments.size() != 1)
            return expression;
        if (arguments[0].isNumber()) {
            if (const auto root = mathematics::exactPrincipalSquareRoot(arguments[0].asNumber()))
                return Expr{*root};

            // exactな有理実数では、平方因子をradicalの外へ出してcanonical化する。
            // sqrt[8] -> 2 sqrt[2], sqrt[32] -> 4 sqrt[2], sqrt[2/3] -> sqrt[6]/3。負数はprincipal branchなのでIを掛ける。
            if (arguments[0].asNumber().isReal()) {
                const RealNumber& real = arguments[0].asNumber().asReal();
                if (!real.isZero()) {
                    const bool negative = real.isNegative();
                    const Rational magnitude = real.abs().toRational();
                    const auto decomposition =
                        mathematics::decomposePositiveRationalSquareRoot(magnitude);

                    const bool unchangedPositiveInteger =
                        decomposition.coefficient == Rational{BigInt{1}}
                        && magnitude.denominator() == BigInt{1}
                        && decomposition.radicand == magnitude.numerator();

                    if (negative || !unchangedPositiveInteger) {
                        Expr radical = decomposition.radicand == BigInt{1}
                            ? integerExpr(1)
                            : Expr::call(
                                context.builtins.symbol(BuiltinId::Sqrt),
                                {Expr{Number{decomposition.radicand}}});
                        Expr normalized = mathematics::scaleExactExpression(
                            decomposition.coefficient, radical, context.builtins);
                        if (!negative)
                            return normalized;

                        const Expr imaginaryUnit{Number::complex(
                            RealNumber{}, RealNumber{BigInt{1}})};
                        return Expr::call(
                            context.builtins.symbol(BuiltinId::Multiply),
                            {imaginaryUnit, std::move(normalized)});
                    }
                }
            }
        }

        // 負実数であることが証明できる場合は、principal sqrtをexact Complexへ昇格する。
        // sqrt[-a] = I sqrt[a] (a > 0)。数値近似は一切使わない。
        const mathematics::ValueFacts inputFacts = knowledge.facts(arguments[0]);
        if (inputFacts.isProvablyNegativeReal()) {
            Expr magnitude = arguments[0];
            if (arguments[0].isNumber()) {
                magnitude = Expr{Number{arguments[0].asNumber().asReal().abs()}};
            }
            else if (isHead(arguments[0], context.builtins, BuiltinId::Negate)
                && arguments[0].asCall().arguments.size() == 1) {
                magnitude = arguments[0].asCall().arguments.front();
            }
            else {
                magnitude = Expr::call(
                    context.builtins.symbol(BuiltinId::Negate), {arguments[0]});
            }

            const Expr positiveRoot = Expr::call(
                context.builtins.symbol(BuiltinId::Sqrt), {std::move(magnitude)});
            const Expr imaginaryUnit{Number::complex(
                RealNumber{}, RealNumber{BigInt{1}})};
            return Expr::call(
                context.builtins.symbol(BuiltinId::Multiply),
                {imaginaryUnit, positiveRoot});
        }

        if (isHead(arguments[0], context.builtins, BuiltinId::Power)) {
            const auto& powerArguments = arguments[0].asCall().arguments;
            if (powerArguments.size() == 2 && isExactReal(powerArguments[1], 2)) {
                const Expr& base = powerArguments[0];
                if (knowledge.prove(mathematics::relation(
                    RelationKind::GreaterEqual, base, integerExpr(0))) == TruthValue::True)
                    return base;
                if (knowledge.prove(mathematics::relation(
                    RelationKind::LessEqual, base, integerExpr(0))) == TruthValue::True)
                    return Expr::call(context.builtins.symbol(BuiltinId::Negate), {base});
                // principal sqrtは一価函数なので、実数xについて sqrt[x^2] = |x|。
                // 方程式を逆に解く際の ±x (PlusMinus) とは意味が異なる。
                if (knowledge.facts(base).isProvablyReal())
                    return Expr::call(context.builtins.symbol(BuiltinId::Abs), {base});
            }
        }
        return expression;
    }

    case BuiltinId::Abs: {
        if (arguments.size() != 1)
            return expression;
        const Expr& input = arguments[0];
        if (input.isNumber()) {
            const Number& number = input.asNumber();
            if (number.isReal())
                return Expr{Number{number.asReal().abs()}};

            const auto& complex = number.asComplex();
            const RealNumber magnitudeSquared =
                complex.real * complex.real + complex.imaginary * complex.imaginary;
            if (const auto exact = mathematics::exactPrincipalSquareRoot(Number{magnitudeSquared}))
                return Expr{*exact};
            return Expr::call(
                context.builtins.symbol(BuiltinId::Sqrt),
                {Expr{Number{magnitudeSquared}}});
        }

        const mathematics::ValueFacts facts = knowledge.facts(input);
        if (facts.sign == mathematics::RealSign::Zero)
            return integerExpr(0);
        if (facts.isProvablyReal()) {
            if (facts.sign == mathematics::RealSign::Positive
                || facts.sign == mathematics::RealSign::NonNegative)
                return input;
            if (facts.sign == mathematics::RealSign::Negative
                || facts.sign == mathematics::RealSign::NonPositive)
                return Expr::call(context.builtins.symbol(BuiltinId::Negate), {input});
        }
        if (isHead(input, context.builtins, BuiltinId::Negate)
            && input.asCall().arguments.size() == 1)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Abs),
                {input.asCall().arguments.front()});
        if (isHead(input, context.builtins, BuiltinId::Abs)
            && input.asCall().arguments.size() == 1)
            return input;
        if (isHead(input, context.builtins, BuiltinId::Conj)
            && input.asCall().arguments.size() == 1)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Abs),
                {input.asCall().arguments.front()});
        return expression;
    }

    case BuiltinId::Sign: {
        if (arguments.size() != 1)
            return expression;
        const Expr& input = arguments[0];
        if (input.isNumber()) {
            const Number& number = input.asNumber();
            if (number.isZero())
                return integerExpr(0);
            if (number.isReal())
                return integerExpr(number.asReal().isNegative() ? -1 : 1);
            return Expr::call(
                context.builtins.symbol(BuiltinId::Divide),
                {input, Expr::call(context.builtins.symbol(BuiltinId::Abs), {input})});
        }

        const mathematics::ValueFacts facts = knowledge.facts(input);
        if (provablyDefined(input, context)) {
            if (facts.sign == mathematics::RealSign::Zero)
                return integerExpr(0);
            if (facts.sign == mathematics::RealSign::Positive)
                return integerExpr(1);
            if (facts.sign == mathematics::RealSign::Negative)
                return integerExpr(-1);
        }
        if (isHead(input, context.builtins, BuiltinId::Negate)
            && input.asCall().arguments.size() == 1)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Negate),
                {Expr::call(context.builtins.symbol(BuiltinId::Sign),
                    {input.asCall().arguments.front()})});
        return expression;
    }

    case BuiltinId::Re: {
        if (arguments.size() != 1)
            return expression;
        const Expr& input = arguments[0];
        if (input.isNumber())
            return Expr{Number{input.asNumber().realPart()}};
        if (knowledge.facts(input).isProvablyReal())
            return input;
        if (isHead(input, context.builtins, BuiltinId::Conj)
            && input.asCall().arguments.size() == 1)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Re),
                {input.asCall().arguments.front()});
        if (isHead(input, context.builtins, BuiltinId::Negate)
            && input.asCall().arguments.size() == 1)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Negate),
                {Expr::call(context.builtins.symbol(BuiltinId::Re),
                    {input.asCall().arguments.front()})});
        if (isHead(input, context.builtins, BuiltinId::Add)) {
            std::vector<Expr> parts;
            parts.reserve(input.asCall().arguments.size());
            for (const Expr& part : input.asCall().arguments)
                parts.push_back(Expr::call(context.builtins.symbol(BuiltinId::Re), {part}));
            return Expr::call(context.builtins.symbol(BuiltinId::Add), std::move(parts));
        }
        return expression;
    }

    case BuiltinId::Im: {
        if (arguments.size() != 1)
            return expression;
        const Expr& input = arguments[0];
        if (input.isNumber())
            return Expr{Number{input.asNumber().imaginaryPart()}};
        if (knowledge.facts(input).isProvablyReal()
            && provablyDefined(input, context))
            return integerExpr(0);
        if (isHead(input, context.builtins, BuiltinId::Conj)
            && input.asCall().arguments.size() == 1)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Negate),
                {Expr::call(context.builtins.symbol(BuiltinId::Im),
                    {input.asCall().arguments.front()})});
        if (isHead(input, context.builtins, BuiltinId::Negate)
            && input.asCall().arguments.size() == 1)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Negate),
                {Expr::call(context.builtins.symbol(BuiltinId::Im),
                    {input.asCall().arguments.front()})});
        if (isHead(input, context.builtins, BuiltinId::Add)) {
            std::vector<Expr> parts;
            parts.reserve(input.asCall().arguments.size());
            for (const Expr& part : input.asCall().arguments)
                parts.push_back(Expr::call(context.builtins.symbol(BuiltinId::Im), {part}));
            return Expr::call(context.builtins.symbol(BuiltinId::Add), std::move(parts));
        }
        return expression;
    }

    case BuiltinId::Conj: {
        if (arguments.size() != 1)
            return expression;
        const Expr& input = arguments[0];
        if (input.isNumber())
            return Expr{input.asNumber().conjugate()};
        if (knowledge.facts(input).isProvablyReal())
            return input;
        if (isHead(input, context.builtins, BuiltinId::Conj)
            && input.asCall().arguments.size() == 1)
            return input.asCall().arguments.front();
        if (isHead(input, context.builtins, BuiltinId::Negate)
            && input.asCall().arguments.size() == 1)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Negate),
                {Expr::call(context.builtins.symbol(BuiltinId::Conj),
                    {input.asCall().arguments.front()})});
        if (isHead(input, context.builtins, BuiltinId::Add)
            || isHead(input, context.builtins, BuiltinId::Multiply)) {
            const BuiltinId operation = isHead(input, context.builtins, BuiltinId::Add)
                ? BuiltinId::Add : BuiltinId::Multiply;
            std::vector<Expr> parts;
            parts.reserve(input.asCall().arguments.size());
            for (const Expr& part : input.asCall().arguments)
                parts.push_back(Expr::call(context.builtins.symbol(BuiltinId::Conj), {part}));
            return Expr::call(context.builtins.symbol(operation), std::move(parts));
        }
        if (isHead(input, context.builtins, BuiltinId::Divide)
            && input.asCall().arguments.size() == 2)
            return Expr::call(
                context.builtins.symbol(BuiltinId::Divide),
                {Expr::call(context.builtins.symbol(BuiltinId::Conj), {input.asCall().arguments[0]}),
                 Expr::call(context.builtins.symbol(BuiltinId::Conj), {input.asCall().arguments[1]})});
        return expression;
    }

    case BuiltinId::Sin:
    case BuiltinId::Cos:
    case BuiltinId::Tan:
    case BuiltinId::Cot:
    case BuiltinId::Sec:
    case BuiltinId::Csc:
    case BuiltinId::Asin:
    case BuiltinId::Acos:
    case BuiltinId::Atan: {
        if (arguments.size() != 1)
            return expression;

        const auto* mathFunction = context.mathematics.findFunction(expression.asCall().head);
        if (mathFunction
            && isHead(arguments[0], context.builtins, BuiltinId::Negate)
            && arguments[0].asCall().arguments.size() == 1) {
            const Expr inner = arguments[0].asCall().arguments.front();
            const Expr sameFunction = Expr::call(expression.asCall().head, {inner});
            if (mathFunction->parity == mathematics::FunctionParity::Even)
                return sameFunction;
            if (mathFunction->parity == mathematics::FunctionParity::Odd)
                return Expr::call(context.builtins.symbol(BuiltinId::Negate), {sameFunction});
        }

        if (!mathFunction)
            return expression;
        if (auto composed = simplifyPrincipalInverseComposition(
            definition->id, arguments[0], context))
            return std::move(*composed);
        if (definition->id == BuiltinId::Asin
            || definition->id == BuiltinId::Acos
            || definition->id == BuiltinId::Atan) {
            if (auto exact = mathematics::simplifyExactInverseTrig(
                mathFunction->id, arguments[0], context.builtins, context.mathematics,
                context.angleSemantics))
                return std::move(*exact);
            return expression;
        }

        if (auto exact = mathematics::simplifyExactTrig(
            mathFunction->id, arguments[0], context.builtins, context.mathematics,
            context.angleSemantics))
            return std::move(*exact);
        return expression;
    }

    case BuiltinId::Atan2:
        if (arguments.size() == 2) {
            if (auto exact = mathematics::simplifyExactAtan2(
                arguments[0], arguments[1], context.builtins, context.mathematics,
                context.angleSemantics))
                return std::move(*exact);
        }
        return expression;

    case BuiltinId::Sinh:
    case BuiltinId::Cosh:
    case BuiltinId::Tanh:
    case BuiltinId::Asinh:
    case BuiltinId::Acosh:
    case BuiltinId::Atanh:
    case BuiltinId::Csch:
    case BuiltinId::Sech:
    case BuiltinId::Coth: {
        if (arguments.size() != 1)
            return expression;

        const auto* mathFunction = context.mathematics.findFunction(expression.asCall().head);
        if (!mathFunction)
            return expression;
        if (isHead(arguments[0], context.builtins, BuiltinId::Negate)
            && arguments[0].asCall().arguments.size() == 1) {
            const Expr inner = arguments[0].asCall().arguments.front();
            const Expr sameFunction = Expr::call(expression.asCall().head, {inner});
            if (mathFunction->parity == mathematics::FunctionParity::Even)
                return sameFunction;
            if (mathFunction->parity == mathematics::FunctionParity::Odd)
                return Expr::call(context.builtins.symbol(BuiltinId::Negate), {sameFunction});
        }

        if (auto composed = simplifyPrincipalInverseComposition(
            definition->id, arguments[0], context))
            return std::move(*composed);

        if (auto exact = mathematics::simplifyExactHyperbolic(
            mathFunction->id, arguments[0], context.builtins, context.mathematics))
            return std::move(*exact);
        return expression;
    }

    case BuiltinId::Arg:
        if (arguments.size() == 1) {
            if (auto exact = mathematics::simplifyExactArg(
                arguments[0], context.builtins, context.mathematics))
                return std::move(*exact);
        }
        return expression;

    case BuiltinId::Log:
        if (arguments.size() == 1) {
            if (isHead(arguments[0], context.builtins, BuiltinId::Exp)
                && arguments[0].asCall().arguments.size() == 1) {
                const Expr& inner = arguments[0].asCall().arguments.front();
                // principal Log[Exp[z]] = z は一般複素数では成立しない。
                // zが実数と証明できる場合はIm[z]=0でprincipal strip内なので安全。
                if (knowledge.facts(inner).isProvablyReal())
                    return inner;
            }
            if (auto exact = mathematics::simplifyExactLog(
                arguments[0], context.builtins, context.mathematics))
                return std::move(*exact);
            return expression;
        }
        if (arguments.size() == 2) {
            if (auto exact = mathematics::simplifyExactLog(
                arguments[0], arguments[1], context.builtins, context.mathematics))
                return std::move(*exact);

            // symbolic baseでは base!=0,1 が証明できるときだけ恒等式を使う。
            const Expr baseMinusOne = Expr::call(
                context.builtins.symbol(BuiltinId::Subtract),
                {arguments[0], integerExpr(1)});
            const bool validBase = provablyDefined(arguments[0], context)
                && proveNonZero(arguments[0], knowledge) == TruthValue::True
                && proveNonZero(baseMinusOne, knowledge) == TruthValue::True;
            if (validBase && isExactReal(arguments[1], 1))
                return integerExpr(0);
            if (validBase && arguments[0] == arguments[1])
                return integerExpr(1);
            return expression;
        }
        return expression;

    case BuiltinId::Log2:
    case BuiltinId::Log10:
        if (arguments.size() == 1) {
            const std::int64_t base = definition->id == BuiltinId::Log2 ? 2 : 10;
            return Expr::call(
                context.builtins.symbol(BuiltinId::Log),
                {integerExpr(base), arguments[0]});
        }
        return expression;

    case BuiltinId::Erf:
        if (arguments.size() == 1 && isExactReal(arguments[0], 0))
            return integerExpr(0);
        if (arguments.size() == 1
            && isHead(arguments[0], context.builtins, BuiltinId::Negate)
            && arguments[0].asCall().arguments.size() == 1)
            return Expr::call(context.builtins.symbol(BuiltinId::Negate), {
                Expr::call(context.builtins.symbol(BuiltinId::Erf),
                    {arguments[0].asCall().arguments.front()})});
        return expression;

    case BuiltinId::Erfc:
        if (arguments.size() == 1 && isExactReal(arguments[0], 0))
            return integerExpr(1);
        if (arguments.size() == 1
            && isHead(arguments[0], context.builtins, BuiltinId::Negate)
            && arguments[0].asCall().arguments.size() == 1)
            return Expr::call(context.builtins.symbol(BuiltinId::Subtract), {
                integerExpr(2),
                Expr::call(context.builtins.symbol(BuiltinId::Erfc),
                    {arguments[0].asCall().arguments.front()})});
        return expression;

    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS:
        if (arguments.size() == 1
            && isHead(arguments[0], context.builtins, BuiltinId::Negate)
            && arguments[0].asCall().arguments.size() == 1)
            return Expr::call(context.builtins.symbol(BuiltinId::Negate), {
                Expr::call(context.builtins.symbol(definition->id),
                    {arguments[0].asCall().arguments.front()})});
        return expression;

    case BuiltinId::Hypergeometric1F1:
        if (arguments.size() == 3
            && isExactReal(arguments[0], 0)
            && provablyDefined(arguments[0], context)
            && provablyDefined(arguments[1], context)
            && provablyDefined(arguments[2], context))
            return integerExpr(1);
        return expression;

    case BuiltinId::Hypergeometric2F1:
        if (arguments.size() == 4
            && (isExactReal(arguments[0], 0) || isExactReal(arguments[1], 0))
            && provablyDefined(arguments[0], context)
            && provablyDefined(arguments[1], context)
            && provablyDefined(arguments[2], context)
            && provablyDefined(arguments[3], context))
            return integerExpr(1);
        return expression;

    case BuiltinId::EllipticF:
    case BuiltinId::EllipticE:
        if (arguments.size() == 2) {
            if (isExactReal(arguments[0], 0)
                && provablyDefined(arguments[1], context))
                return integerExpr(0);
            if (isExactReal(arguments[1], 0))
                return arguments[0];
        }
        return expression;

    case BuiltinId::EllipticPi:
        if (arguments.size() == 3) {
            if (isExactReal(arguments[1], 0)
                && provablyDefined(arguments[0], context)
                && provablyDefined(arguments[2], context))
                return integerExpr(0);
            if (isExactReal(arguments[0], 0))
                return Expr::call(context.builtins.symbol(BuiltinId::EllipticF),
                    {arguments[1], arguments[2]});
        }
        return expression;

    case BuiltinId::Exp:
        if (arguments.size() != 1)
            return expression;
        if (isHead(arguments[0], context.builtins, BuiltinId::Log)
            && arguments[0].asCall().arguments.size() == 1) {
            const Expr& inner = arguments[0].asCall().arguments.front();
            if (proveNonZero(inner, knowledge) == TruthValue::True)
                return inner;
        }
        if (auto exact = mathematics::simplifyExactExp(
            arguments[0], context.builtins, context.mathematics))
            return std::move(*exact);
        return expression;

    case BuiltinId::Polylog:
        if (arguments.size() == 2 && isExactReal(arguments[1], 0)
            && provablyDefined(arguments[0], context))
            return integerExpr(0);
        return expression;

    case BuiltinId::GeneralizedBinomial:
    case BuiltinId::FallingFactorial:
    case BuiltinId::RisingFactorial:
        if (arguments.size() == 2 && isExactReal(arguments[1], 0)
            && provablyDefined(arguments[0], context))
            return integerExpr(1);
        return expression;

    case BuiltinId::Polar:
    case BuiltinId::NextPow2:
    case BuiltinId::DegreeToRadian:
    case BuiltinId::DegreeToGradian:
    case BuiltinId::RadianToDegree:
    case BuiltinId::RadianToGradian:
    case BuiltinId::GradianToDegree:
    case BuiltinId::GradianToRadian:
    case BuiltinId::Expm1:
    case BuiltinId::Log1p:
    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc:
    case BuiltinId::Sinhc:
    case BuiltinId::Tanhc:
    case BuiltinId::Expc:
    case BuiltinId::Gamma:
    case BuiltinId::LogGamma:
    case BuiltinId::LambertW:
    case BuiltinId::ExponentialIntegralEi:
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
    case BuiltinId::LogarithmicIntegralLi:
    case BuiltinId::Beta:
    case BuiltinId::BetaLog:
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
    case BuiltinId::Trace:
    case BuiltinId::Rows:
    case BuiltinId::Cols:
    case BuiltinId::Diag:
    case BuiltinId::VectorAdd:
    case BuiltinId::VectorSubtract:
    case BuiltinId::VectorScale:
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
    case BuiltinId::VectorInner:
    case BuiltinId::VectorOuter:
    case BuiltinId::VectorRejection:
    case BuiltinId::OrthogonalQ:
    case BuiltinId::OrthonormalQ:
    case BuiltinId::LinearIndependentQ:
    case BuiltinId::GramSchmidt:
    case BuiltinId::Gradient:
    case BuiltinId::Divergence:
    case BuiltinId::Curl:
    case BuiltinId::Laplacian:
    case BuiltinId::Jacobian:
    case BuiltinId::Hessian:
    case BuiltinId::DirectionalDerivative:
    case BuiltinId::Factorial:
    case BuiltinId::Derivative:
    case BuiltinId::SymbolicIntegral:
    case BuiltinId::Limit:
    case BuiltinId::Floor:
    case BuiltinId::Ceil:
    case BuiltinId::Trunc:
    case BuiltinId::Round:
    case BuiltinId::Frac:
    case BuiltinId::BitAnd:
    case BuiltinId::BitOr:
    case BuiltinId::BitXor:
    case BuiltinId::BitNot:
    case BuiltinId::BitShiftLeft:
    case BuiltinId::BitShiftRight:
    case BuiltinId::BitLength:
    case BuiltinId::BitCount:
    case BuiltinId::BitGet:
    case BuiltinId::Fma:
    case BuiltinId::Clamp:
    case BuiltinId::Proj:
    case BuiltinId::Gcd:
    case BuiltinId::Lcm:
    case BuiltinId::Mod:
    case BuiltinId::Rem:
    case BuiltinId::Quotient:
    case BuiltinId::IsPrime:
    case BuiltinId::NextPrime:
    case BuiltinId::PreviousPrime:
    case BuiltinId::FactorInteger:
    case BuiltinId::Totient:
    case BuiltinId::Zeta:
    case BuiltinId::Digamma:
    case BuiltinId::Trigamma:
    case BuiltinId::IncompleteBeta:
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
    case BuiltinId::ConditionNumber:
    case BuiltinId::LeastSquares:
    case BuiltinId::PseudoInverse:
    case BuiltinId::Eigenvalues:
    case BuiltinId::Eigenvectors:
    case BuiltinId::Eigensystem:
    case BuiltinId::ConjugateTranspose:
    case BuiltinId::Length:
    case BuiltinId::NumericDerivative:
    case BuiltinId::NumericIntegral:
    case BuiltinId::NumericalApproximation:
    case BuiltinId::Precision:
    case BuiltinId::Accuracy:
    case BuiltinId::Explain:
    case BuiltinId::Map:
    case BuiltinId::Range:
    case BuiltinId::Table:
    case BuiltinId::Rationalize:
    case BuiltinId::Root:
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify:
    case BuiltinId::Expand:
    case BuiltinId::Factor:
    case BuiltinId::Collect:
    case BuiltinId::Cases: {
        std::vector<Expr> result;
        result.reserve(arguments.size());
        bool unresolvedBefore = false;
        for (const Expr& branchExpression : arguments) {
            if (!isHead(branchExpression, context.builtins, BuiltinId::CaseBranch)
                || branchExpression.asCall().arguments.empty()
                || branchExpression.asCall().arguments.size() > 2)
                return expression;
            const auto& branch = branchExpression.asCall().arguments;
            if (branch.size() == 1) {
                if (!unresolvedBefore)
                    return branch[0];
                result.push_back(branchExpression);
                break;
            }
            TruthValue branchTruth = TruthValue::Unknown;
            if (branch[1].isBoolean())
                branchTruth = branch[1].asBoolean() ? TruthValue::True : TruthValue::False;
            else if (const auto predicate = conditionPredicate(branch[1], context.builtins))
                branchTruth = knowledge.prove(*predicate);
            if (branchTruth != TruthValue::Unknown) {
                if (branchTruth == TruthValue::False)
                    continue;
                if (!unresolvedBefore)
                    return branch[0];
                result.push_back(Expr::call(
                    context.builtins.symbol(BuiltinId::CaseBranch), {branch[0]}));
                break;
            }
            unresolvedBefore = true;
            result.push_back(branchExpression);
        }
        if (result.empty()) {
            if (const auto value = predefinedValue(
                    context, symbols::PredefinedSymbolId::Indeterminate))
                return *value;
            return expression;
        }
        return Expr::call(context.builtins.symbol(BuiltinId::Cases), std::move(result));
    }
    case BuiltinId::CaseBranch:
        return expression;
    case BuiltinId::Solve:
    case BuiltinId::GroebnerBasis:
    case BuiltinId::PolynomialReduce:
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
    case BuiltinId::Series:
    case BuiltinId::SeriesData:
    case BuiltinId::Normal:
    case BuiltinId::ToNormal:
    case BuiltinId::UnitApplied:
        return expression;
    }

    return expression;
}

[[nodiscard]] Expr simplifyOnePass(
    const Expr& root,
    const SimplificationContext& context) {
    struct Frame final {
        Expr expression;
        bool childrenDone = false;
    };

    std::vector<Frame> stack{{root, false}};
    std::unordered_map<const void*, Expr> completed;
    completed.reserve(64);

    while (!stack.empty()) {
        Frame current = std::move(stack.back());
        stack.pop_back();

        if (completed.find(current.expression.identity()) != completed.end())
            continue;

        if (!current.childrenDone) {
            if ((!current.expression.isCall() && !current.expression.isArray()
                    && !current.expression.isList())
                || (current.expression.isArray()
                    && !current.expression.asArray().hasStoredExpressions())) {
                completed.emplace(current.expression.identity(), current.expression);
                continue;
            }

            stack.push_back(Frame{current.expression, true});
            if (current.expression.isCall()) {
                const auto& arguments = current.expression.asCall().arguments;
                for (auto iterator = arguments.rbegin(); iterator != arguments.rend(); ++iterator)
                    if (completed.find(iterator->identity()) == completed.end())
                        stack.push_back(Frame{*iterator, false});
            }
            else if (current.expression.isArray()) {
                const auto elements = current.expression.asArray().storedExpressions();
                for (auto iterator = elements.rbegin(); iterator != elements.rend(); ++iterator)
                    if (completed.find(iterator->identity()) == completed.end())
                        stack.push_back(Frame{*iterator, false});
            }
            else {
                const auto& elements = current.expression.asList().elements;
                for (auto iterator = elements.rbegin(); iterator != elements.rend(); ++iterator)
                    if (completed.find(iterator->identity()) == completed.end())
                        stack.push_back(Frame{*iterator, false});
            }
            continue;
        }

        Expr rebuilt = current.expression;
        if (current.expression.isCall()) {
            const auto& sourceCall = current.expression.asCall();
            std::vector<Expr> arguments;
            arguments.reserve(sourceCall.arguments.size());
            for (const Expr& child : sourceCall.arguments)
                arguments.push_back(completed.at(child.identity()));
            rebuilt = Expr::rebuildCall(sourceCall, std::move(arguments));
            rebuilt = localRewrite(rebuilt, context);
        }
        else if (current.expression.isArray()) {
            const auto& array = current.expression.asArray();
            const auto entries = array.expressionEntries();
            std::vector<std::size_t> indices;
            std::vector<Expr> elements;
            indices.reserve(entries.size());
            elements.reserve(entries.size());
            for (const auto& entry : entries) {
                indices.push_back(entry.index);
                elements.push_back(completed.at(entry.expression.identity()));
            }
            rebuilt = Expr::array(array.replacedExpressions(indices, std::move(elements)));
        }
        else if (current.expression.isList()) {
            std::vector<Expr> elements;
            elements.reserve(current.expression.asList().elements.size());
            for (const Expr& child : current.expression.asList().elements)
                elements.push_back(completed.at(child.identity()));
            rebuilt = expression::braceValue(std::move(elements));
        }

        completed.emplace(current.expression.identity(), std::move(rebuilt));
    }

    return completed.at(root.identity());
}

} // namespace

namespace {

struct ExplicitLinearTerm final {
    Rational coefficient{BigInt{0}};
    Expr atom{Number{BigInt{1}}};
};

[[nodiscard]] bool isAdditiveForExplicitLinearCombination(
    const Expr& expression,
    const SimplificationContext& context,
    BuiltinId* id = nullptr) {
    if (!expression.isCall()) return false;
    const auto* definition = context.builtins.find(expression.asCall().head);
    if (!definition) return false;
    if (id) *id = definition->id;
    return definition->id == BuiltinId::Add
        || definition->id == BuiltinId::Subtract
        || definition->id == BuiltinId::Negate;
}

void collectExplicitLinearTerms(
    const Expr& expression,
    Rational coefficient,
    std::vector<ExplicitLinearTerm>& terms,
    const SimplificationContext& context,
    std::size_t depth = 0) {
    if (depth > 256) {
        terms.push_back(ExplicitLinearTerm{std::move(coefficient), expression});
        return;
    }

    BuiltinId id = BuiltinId::Add;
    if (isAdditiveForExplicitLinearCombination(expression, context, &id)) {
        const auto& arguments = expression.asCall().arguments;
        if (id == BuiltinId::Add) {
            for (const Expr& argument : arguments)
                collectExplicitLinearTerms(
                    argument, coefficient, terms, context, depth + 1);
            return;
        }
        if (id == BuiltinId::Subtract && arguments.size() == 2) {
            collectExplicitLinearTerms(
                arguments[0], coefficient, terms, context, depth + 1);
            collectExplicitLinearTerms(
                arguments[1], -coefficient, terms, context, depth + 1);
            return;
        }
        if (id == BuiltinId::Negate && arguments.size() == 1) {
            collectExplicitLinearTerms(
                arguments[0], -coefficient, terms, context, depth + 1);
            return;
        }
    }

    if (expression.isCall()) {
        const auto* definition = context.builtins.find(expression.asCall().head);
        if (definition && definition->id == BuiltinId::Multiply) {
            Rational scalar{BigInt{1}};
            std::vector<Expr> nonNumeric;
            nonNumeric.reserve(expression.asCall().arguments.size());
            for (const Expr& factor : expression.asCall().arguments) {
                if (factor.isNumber() && factor.asNumber().isReal()) {
                    scalar *= factor.asNumber().asReal().toRational();
                    continue;
                }
                nonNumeric.push_back(factor);
            }
            coefficient *= scalar;
            if (nonNumeric.empty()) {
                terms.push_back(ExplicitLinearTerm{
                    std::move(coefficient), Expr{Number{BigInt{1}}}});
                return;
            }
            if (nonNumeric.size() == 1
                && isAdditiveForExplicitLinearCombination(nonNumeric.front(), context)) {
                collectExplicitLinearTerms(
                    nonNumeric.front(), coefficient, terms, context, depth + 1);
                return;
            }
            Expr atom = nonNumeric.size() == 1
                ? nonNumeric.front()
                : Expr::call(
                    context.builtins.symbol(BuiltinId::Multiply),
                    std::move(nonNumeric));
            terms.push_back(ExplicitLinearTerm{std::move(coefficient), std::move(atom)});
            return;
        }
    }

    terms.push_back(ExplicitLinearTerm{std::move(coefficient), expression});
}

[[nodiscard]] Expr buildExplicitLinearTerm(
    const Rational& coefficient,
    const Expr& atom,
    const SimplificationContext& context) {
    if (coefficient == Rational{BigInt{1}})
        return atom;
    if (coefficient == Rational{BigInt{-1}})
        return Expr::call(
            context.builtins.symbol(BuiltinId::Negate), {atom});
    if (atom.isNumber() && atom.asNumber().isReal()
        && atom.asNumber().asReal().toRational() == Rational{BigInt{1}})
        return Expr{Number{coefficient}};
    return Expr::call(
        context.builtins.symbol(BuiltinId::Multiply),
        {Expr{Number{coefficient}}, atom});
}

} // namespace

Expr simplifyExplicitLinearCombination(
    const Expr& expression,
    const SimplificationContext& context) {
    Expr simplified = Simplifier{}.simplify(expression, context);
    if (!isAdditiveForExplicitLinearCombination(simplified, context))
        return simplified;

    std::vector<ExplicitLinearTerm> terms;
    terms.reserve(32);
    collectExplicitLinearTerms(
        simplified, Rational{BigInt{1}}, terms, context);
    if (terms.size() < 2)
        return simplified;

    struct Group final {
        Expr atom;
        Rational sum{BigInt{0}};
    };
    std::vector<Group> groups;
    groups.reserve(terms.size());
    for (const auto& term : terms) {
        auto iterator = std::find_if(
            groups.begin(), groups.end(),
            [&](const Group& group) { return group.atom == term.atom; });
        if (iterator == groups.end())
            groups.push_back(Group{term.atom, term.coefficient});
        else
            iterator->sum += term.coefficient;
    }

    // 0への相殺は定義域を消し得る。明示前提込みでatomのdefinednessを証明できない
    // groupが一つでもあれば，通常simplifyの結果をそのまま返す。
    for (const Group& group : groups)
        if (group.sum.isZero() && !provablyDefined(group.atom, context))
            return simplified;

    std::vector<Expr> rebuilt;
    rebuilt.reserve(groups.size());
    for (const Group& group : groups) {
        if (group.sum.isZero())
            continue;
        rebuilt.push_back(buildExplicitLinearTerm(group.sum, group.atom, context));
    }
    if (rebuilt.empty())
        return Expr{Number{BigInt{0}}};
    Expr result = rebuilt.size() == 1
        ? std::move(rebuilt.front())
        : Expr::call(context.builtins.symbol(BuiltinId::Add), std::move(rebuilt));
    return Simplifier{}.simplify(result, context);
}

Expr Simplifier::simplify(
    const Expr& expression,
    const SimplificationContext& context) const {
    if (context.maximumPasses == 0)
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "Simplification pass limit must be positive");

    Expr current = expression;
    for (std::size_t pass = 0; pass < context.maximumPasses; ++pass) {
        if (context.budget)
            context.budget->consume(evaluation::EvaluationResource::SimplificationCandidate);
        else
            evaluation::consumeEvaluationBudget(
                evaluation::EvaluationResource::SimplificationCandidate);
        Expr next = simplifyOnePass(current, context);
        if (next == current)
            return next;
        current = std::move(next);
    }

    error::throwCalcError(
        error::CalcErrorType::Internal,
        "Simplification pass limit exceeded");
}

} // namespace mmcal::simplification
