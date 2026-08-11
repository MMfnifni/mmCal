// 候補探索型FullSimplify
#include "full_simplifier.hpp"

#include "expression_cost.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "simplifier.hpp"
#include "symbolic/algebra_transforms.hpp"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <optional>
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
    return expression.isCall()
        && expression.asCall().head.sameIdentity(builtins.symbol(id));
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
        else if (current.isArray()) {
            for (const Expr& element : current.asArray().elements)
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
        for (std::size_t i = 0; i < array.elements.size(); ++i) {
            const std::vector<Expr> transformed = rootVariants(array.elements[i], context);
            for (const Expr& replacement : transformed) {
                if (replacement == array.elements[i])
                    continue;
                std::vector<Expr> elements = array.elements;
                elements[i] = replacement;
                result.push_back(Expr::array(array.shape, std::move(elements)));
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
