// 候補探索型FullSimplify
#include "full_simplifier.hpp"
#include "expression/array_utils.hpp"

#include "expression_cost.hpp"
#include "mathematics/predicate.hpp"
#include "mathematics/trigonometric_polynomial.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "simplifier.hpp"
#include "symbolic/algebra_transforms.hpp"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <optional>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace mmcal::simplification {
namespace {

using expression::Expr;

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    evaluation::BuiltinId id) {
    return builtins.isCallTo(expression, id);
}

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{numeric::Number{numeric::BigInt{value}}};
}

[[nodiscard]] Expr product(
    std::vector<Expr> factors,
    const evaluation::BuiltinRegistry& builtins) {
    if (factors.empty())
        return integer(1);
    if (factors.size() == 1)
        return factors.front();
    return Expr::call(builtins.symbol(evaluation::BuiltinId::Multiply), std::move(factors));
}

[[nodiscard]] std::optional<Expr> cancelProvablyNonzeroCommonFactor(
    const Expr& expression,
    const SimplificationContext& context) {
    if (!isHead(expression, context.builtins, evaluation::BuiltinId::Divide))
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() != 2)
        return std::nullopt;

    std::vector<Expr> numeratorFactors =
        isHead(arguments[0], context.builtins, evaluation::BuiltinId::Multiply)
        ? arguments[0].asCall().arguments : std::vector<Expr>{arguments[0]};
    std::vector<Expr> denominatorFactors =
        isHead(arguments[1], context.builtins, evaluation::BuiltinId::Multiply)
        ? arguments[1].asCall().arguments : std::vector<Expr>{arguments[1]};

    const mathematics::KnowledgeContext knowledge = context.knowledge();
    for (std::size_t i = 0; i < numeratorFactors.size(); ++i) {
        for (std::size_t j = 0; j < denominatorFactors.size(); ++j) {
            if (!(numeratorFactors[i] == denominatorFactors[j]))
                continue;
            if (knowledge.prove(mathematics::relation(
                    mathematics::RelationKind::NotEqual,
                    numeratorFactors[i], integer(0))) != mathematics::TruthValue::True)
                continue;

            numeratorFactors.erase(numeratorFactors.begin() + static_cast<std::ptrdiff_t>(i));
            denominatorFactors.erase(denominatorFactors.begin() + static_cast<std::ptrdiff_t>(j));
            Expr numerator = product(std::move(numeratorFactors), context.builtins);
            Expr denominator = product(std::move(denominatorFactors), context.builtins);
            if (denominator.isNumber() && denominator.asNumber() == numeric::Number{numeric::BigInt{1}})
                return numerator;
            return Expr::call(
                context.builtins.symbol(evaluation::BuiltinId::Divide),
                {std::move(numerator), std::move(denominator)});
        }
    }
    return std::nullopt;
}


[[nodiscard]] std::optional<Expr> reciprocalTrigVariant(
    const Expr& expression,
    const SimplificationContext& context) {
    using evaluation::BuiltinId;
    const auto negativeInteger = [&](std::int64_t magnitude) {
        return Expr{numeric::Number{numeric::BigInt{-magnitude}}};
    };
    const Expr* base = &expression;
    std::optional<std::int64_t> exponent;
    if (isHead(expression, context.builtins, BuiltinId::Power)
        && expression.asCall().arguments.size() == 2) {
        base = &expression.asCall().arguments[0];
        const Expr& power = expression.asCall().arguments[1];
        if (!power.isNumber() || !power.asNumber().isReal()
            || !power.asNumber().asReal().isInteger())
            return std::nullopt;
        const auto converted = numeric::tryToUint64(power.asNumber().asReal().asInteger());
        if (!converted || *converted == 0 || *converted > static_cast<std::uint64_t>(INT64_MAX))
            return std::nullopt;
        exponent = static_cast<std::int64_t>(*converted);
    }

    if (!base->isCall() || base->asCall().arguments.size() != 1)
        return std::nullopt;
    const auto* definition = context.builtins.find(base->asCall().head);
    if (!definition)
        return std::nullopt;
    const Expr& u = base->asCall().arguments[0];

    BuiltinId direct = BuiltinId::Sin;
    switch (definition->id) {
    case BuiltinId::Csc: direct = BuiltinId::Sin; break;
    case BuiltinId::Sec: direct = BuiltinId::Cos; break;
    default: return std::nullopt;
    }
    const std::int64_t powerValue = exponent.value_or(1);
    return Expr::call(context.builtins.symbol(BuiltinId::Power), {
        Expr::call(context.builtins.symbol(direct), {u}), negativeInteger(powerValue)});
}

[[nodiscard]] std::optional<Expr> distributeProductSquareVariant(
    const Expr& expression,
    const SimplificationContext& context) {
    using evaluation::BuiltinId;
    if (!isHead(expression, context.builtins, BuiltinId::Power)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    const Expr& exponent = expression.asCall().arguments[1];
    if (!exponent.isNumber() || !exponent.asNumber().isReal()
        || !exponent.asNumber().asReal().isInteger()
        || exponent.asNumber().asReal().asInteger() != numeric::BigInt{2})
        return std::nullopt;
    const Expr& base = expression.asCall().arguments[0];
    if (!isHead(base, context.builtins, BuiltinId::Multiply))
        return std::nullopt;

    // 整数二乗 (ab)^2=a^2b^2 は複素数全体でbranch条件なしに成立する。
    // 通常Simplifierでは式膨張を避け、FullSimplifyのproof候補としてだけ使う。
    std::vector<Expr> factors;
    factors.reserve(base.asCall().arguments.size());
    for (const Expr& factor : base.asCall().arguments)
        factors.push_back(Expr::call(context.builtins.symbol(BuiltinId::Power), {factor, exponent}));
    return Expr::call(context.builtins.symbol(BuiltinId::Multiply), std::move(factors));
}

[[nodiscard]] std::optional<std::uint64_t> positiveIntegerExponent(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    const auto converted = numeric::tryToUint64(expression.asNumber().asReal().asInteger());
    if (!converted || *converted == 0)
        return std::nullopt;
    return *converted;
}

[[nodiscard]] Expr proofTrigPowerRewrite(
    const Expr& expression,
    const SimplificationContext& context,
    bool& changed) {
    using evaluation::BuiltinId;

    Expr current = expression;
    if (expression.isCall()) {
        const auto& call = expression.asCall();
        std::vector<Expr> arguments;
        arguments.reserve(call.arguments.size());
        for (const Expr& argument : call.arguments)
            arguments.push_back(proofTrigPowerRewrite(argument, context, changed));
        current = Expr::rebuildCall(call, std::move(arguments));
    }
    else if (expression.isArray()) {
        const auto& array = expression.asArray();
        const auto entries = array.expressionEntries();
        if (!entries.empty()) {
            std::vector<std::size_t> indices;
            std::vector<Expr> elements;
            indices.reserve(entries.size());
            elements.reserve(entries.size());
            for (const auto& entry : entries) {
                indices.push_back(entry.index);
                elements.push_back(proofTrigPowerRewrite(entry.expression, context, changed));
            }
            current = Expr::array(array.replacedExpressions(indices, std::move(elements)));
        }
    }
    else if (expression.isList()) {
        const auto& list = expression.asList();
        std::vector<Expr> elements;
        elements.reserve(list.elements.size());
        for (const Expr& element : list.elements)
            elements.push_back(proofTrigPowerRewrite(element, context, changed));
        current = expression::braceValue(std::move(elements));
    }

    // 明示角度単位が現在の既定角度と同じなら、同じsession意味論では完全に同値。
    // Fresnel微分は常にRadを明示するため、default Radのintegrandとのderivative-backで
    // UnitAppliedだけが差として残らないようproof候補内で正規化する。
    if (current.isCall() && current.asCall().arguments.size() == 1) {
        const auto* definition = context.builtins.find(current.asCall().head);
        if (definition) {
            const bool directTrig = definition->id == BuiltinId::Sin
                || definition->id == BuiltinId::Cos
                || definition->id == BuiltinId::Tan
                || definition->id == BuiltinId::Cot
                || definition->id == BuiltinId::Sec
                || definition->id == BuiltinId::Csc
                || definition->id == BuiltinId::Cis;
            const Expr& argument = current.asCall().arguments[0];
            if (directTrig && isHead(argument, context.builtins, BuiltinId::UnitApplied)
                && argument.asCall().arguments.size() == 2
                && argument.asCall().arguments[1].isString()) {
                const auto explicitUnit = mathematics::AngleSemantics::parseUnit(
                    argument.asCall().arguments[1].asString());
                if (explicitUnit && *explicitUnit == context.angleSemantics.defaultUnit()) {
                    changed = true;
                    return Expr::call(current.asCall().head, {argument.asCall().arguments[0]});
                }
            }
        }
    }

    if (!isHead(current, context.builtins, BuiltinId::Power)
        || current.asCall().arguments.size() != 2)
        return current;

    const Expr& base = current.asCall().arguments[0];
    const Expr& exponent = current.asCall().arguments[1];
    const auto outerPower = positiveIntegerExponent(exponent);
    if (!outerPower)
        return current;

    // (ab)^2=a^2b^2 は整数二乗なのでprincipal branchに依存しない。
    // FresnelのD[F]-f等では対象部分式が三角函数のさらに内側にあるため、
    // FullSimplifyのproof候補として深い位置にも一括適用する。
    if (*outerPower == 2 && isHead(base, context.builtins, BuiltinId::Multiply)) {
        std::vector<Expr> factors;
        factors.reserve(base.asCall().arguments.size());
        for (const Expr& factor : base.asCall().arguments)
            factors.push_back(Expr::call(
                context.builtins.symbol(BuiltinId::Power), {factor, exponent}));
        changed = true;
        return Expr::call(context.builtins.symbol(BuiltinId::Multiply), std::move(factors));
    }

    // (a^m)^n=a^(mn) はm,nが正整数なら0を含め複素数全体で安全。
    // 一般複素指数のPower則へは拡張せず、proof候補に限定する。
    if (isHead(base, context.builtins, BuiltinId::Power)
        && base.asCall().arguments.size() == 2) {
        const auto innerPower = positiveIntegerExponent(base.asCall().arguments[1]);
        if (innerPower && *innerPower <= UINT64_MAX / *outerPower) {
            changed = true;
            return Expr::call(context.builtins.symbol(BuiltinId::Power), {
                base.asCall().arguments[0],
                Expr{numeric::Number{numeric::BigInt::fromUnsigned(*innerPower * *outerPower)}}});
        }
    }

    if (!base.isCall() || base.asCall().arguments.size() != 1)
        return current;
    const auto* definition = context.builtins.find(base.asCall().head);
    if (!definition)
        return current;
    const Expr& u = base.asCall().arguments[0];

    // tan^2=cos^-2-1, cot^2=sin^-2-1 はmeromorphic identityとして同じpole集合で成立する。
    // sec/cscを一度挟まずdirect reciprocal powerへ落とすことで、derivative-back中の
    // sec^n/cos^-n混在を一候補で同じnormal formへ揃える。通常Simplifierには常設しない。
    if (*outerPower == 2
        && (definition->id == BuiltinId::Tan || definition->id == BuiltinId::Cot)) {
        const BuiltinId direct = definition->id == BuiltinId::Tan
            ? BuiltinId::Cos : BuiltinId::Sin;
        Expr reciprocalSquared = Expr::call(context.builtins.symbol(BuiltinId::Power), {
            Expr::call(context.builtins.symbol(direct), {u}), integer(-2)});
        changed = true;
        return Expr::call(context.builtins.symbol(BuiltinId::Subtract), {
            std::move(reciprocalSquared), integer(1)});
    }

    // sec^n=cos^-n, csc^n=sin^-n（n>0）も同じpole集合を持つ安全な整数冪恒等式。
    if (definition->id == BuiltinId::Sec || definition->id == BuiltinId::Csc) {
        const BuiltinId direct = definition->id == BuiltinId::Sec
            ? BuiltinId::Cos : BuiltinId::Sin;
        const numeric::BigInt negativePower = -numeric::BigInt::fromUnsigned(*outerPower);
        changed = true;
        return Expr::call(context.builtins.symbol(BuiltinId::Power), {
            Expr::call(context.builtins.symbol(direct), {u}),
            Expr{numeric::Number{negativePower}}});
    }

    return current;
}

[[nodiscard]] std::optional<Expr> deepProofTrigPowerVariant(
    const Expr& expression,
    const SimplificationContext& context) {
    bool changed = false;
    Expr candidate = proofTrigPowerRewrite(expression, context, changed);
    if (!changed)
        return std::nullopt;
    return candidate;
}

[[nodiscard]] bool lessCost(
    const ExpressionCost& lhs,
    const ExpressionCost& rhs) noexcept {
    return std::tie(lhs.nodes, lhs.depth, lhs.leaves)
        < std::tie(rhs.nodes, rhs.depth, rhs.leaves);
}

[[nodiscard]] bool contains(
    const std::vector<Expr>& expressions,
    const Expr& candidate) {
    return std::find(expressions.begin(), expressions.end(), candidate) != expressions.end();
}

[[nodiscard]] std::vector<expression::Symbol> collectVariables(
    const Expr& root,
    const mathematics::MathRegistry& mathematics) {
    std::vector<expression::Symbol> variables;
    std::vector<Expr> pending{root};
    while (!pending.empty()) {
        Expr current = std::move(pending.back());
        pending.pop_back();

        if (current.isSymbol()) {
            if (!mathematics.findConstant(current.asSymbol())) {
                const auto duplicate = std::find_if(
                    variables.begin(), variables.end(), [&](const expression::Symbol& symbol) {
                        return symbol.sameIdentity(current.asSymbol());
                    });
                if (duplicate == variables.end())
                    variables.push_back(current.asSymbol());
            }
            continue;
        }

        if (current.isCall()) {
            for (const Expr& argument : current.asCall().arguments)
                pending.push_back(argument);
        }
        else if (current.isArray()
            && current.asArray().storageKind() == expression::ArrayStorageKind::Generic) {
            for (const Expr& element : current.asArray().storedExpressions())
                pending.push_back(element);
        }
        else if (current.isList()) {
            for (const Expr& element : current.asList().elements)
                pending.push_back(element);
        }
    }
    return variables;
}

[[nodiscard]] std::vector<Expr> rootVariants(
    const Expr& expression,
    const SimplificationContext& context) {
    std::vector<Expr> variants;
    variants.push_back(symbolic::expandExpression(
        expression, context.builtins, context.mathematics, context.angleSemantics));
    variants.push_back(symbolic::factorExpression(
        expression, context.builtins, context.mathematics, context.angleSemantics));

    const auto variables = collectVariables(expression, context.mathematics);
    for (const expression::Symbol& variable : variables)
        variants.push_back(symbolic::collectExpression(
            expression, variable,
            context.builtins, context.mathematics, context.angleSemantics));
    if (const auto cancelled = cancelProvablyNonzeroCommonFactor(expression, context))
        variants.push_back(*cancelled);

    if (const auto reciprocal = reciprocalTrigVariant(expression, context))
        variants.push_back(*reciprocal);
    if (const auto distributedSquare = distributeProductSquareVariant(expression, context))
        variants.push_back(*distributedSquare);
    if (const auto deepProof = deepProofTrigPowerVariant(expression, context)) {
        variants.push_back(*deepProof);
        // tan^2 -> sec^2-1 のようなproof恒等式は積の中では展開後に相殺が見える。
        // candidate上限を探索順だけで浪費しないよう、この安全なproof候補だけ直接expandも試す。
        variants.push_back(symbolic::expandExpression(
            *deepProof, context.builtins, context.mathematics, context.angleSemantics));
    }

    // 高次三角冪を常時展開すると式が大きくなるためSimplifierの既定規則にはしない。
    // FullSimplifyのbounded candidateとしてだけ有限Fourier恒等式を試し、
    // D[F]-fのように展開後に大きく相殺できる場合だけcost比較で採用する。
    if (const auto reduced = mathematics::reduceTrigMonomial(expression, context.builtins))
        variants.push_back(*reduced);
    if (const auto reduced = mathematics::reduceTrigProduct(expression, context.builtins))
        variants.push_back(*reduced);
    return variants;
}

[[nodiscard]] std::vector<Expr> childVariants(
    const Expr& expression,
    const SimplificationContext& context) {
    std::vector<Expr> result;
    if (expression.isCall()) {
        const auto& call = expression.asCall();
        for (std::size_t i = 0; i < call.arguments.size(); ++i) {
            const std::vector<Expr> transformed = rootVariants(call.arguments[i], context);
            for (const Expr& replacement : transformed) {
                if (replacement == call.arguments[i])
                    continue;
                std::vector<Expr> arguments = call.arguments;
                arguments[i] = replacement;
                result.push_back(Expr::call(call.head, std::move(arguments)));
            }
        }
    }
    else if (expression.isArray()) {
        const auto& array = expression.asArray();
        for (const auto& entry : array.expressionEntries()) {
            const std::vector<Expr> transformed = rootVariants(entry.expression, context);
            for (const Expr& replacement : transformed) {
                if (replacement == entry.expression)
                    continue;
                const std::size_t index = entry.index;
                result.push_back(Expr::array(array.replacedExpressions(
                    std::span<const std::size_t>{&index, 1}, {replacement})));
            }
        }
    }
    else if (expression.isList()) {
        const auto& list = expression.asList();
        for (std::size_t i = 0; i < list.elements.size(); ++i) {
            const std::vector<Expr> transformed = rootVariants(list.elements[i], context);
            for (const Expr& replacement : transformed) {
                if (replacement == list.elements[i])
                    continue;
                std::vector<Expr> elements = list.elements;
                elements[i] = replacement;
                result.push_back(expression::braceValue(std::move(elements)));
            }
        }
    }
    return result;
}

} // namespace

Expr fullSimplify(
    const Expr& expression,
    const SimplificationContext& context,
    FullSimplificationOptions options) {
    if (options.maximumCandidates == 0)
        return Simplifier{}.simplify(expression, context);

    const Simplifier simplifier;
    Expr initial = simplifier.simplify(expression, context);
    Expr best = initial;
    ExpressionCost bestCost = measureExpressionCost(best);

    std::vector<Expr> seen;
    seen.reserve(options.maximumCandidates);
    seen.push_back(initial);
    std::deque<Expr> queue;
    queue.push_back(initial);

    auto consider = [&](Expr candidate) {
        if (seen.size() >= options.maximumCandidates)
            return;
        candidate = simplifier.simplify(candidate, context);
        if (contains(seen, candidate))
            return;

        const ExpressionCost cost = measureExpressionCost(candidate);
        if (lessCost(cost, bestCost)) {
            best = candidate;
            bestCost = cost;
        }
        seen.push_back(candidate);
        queue.push_back(std::move(candidate));
    };

    while (!queue.empty() && seen.size() < options.maximumCandidates) {
        Expr current = std::move(queue.front());
        queue.pop_front();
        for (Expr candidate : rootVariants(current, context))
            consider(std::move(candidate));
        for (Expr candidate : childVariants(current, context))
            consider(std::move(candidate));
    }

    return best;
}

} // namespace mmcal::simplification
