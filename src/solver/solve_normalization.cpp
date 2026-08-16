// Solve専用の副作用を持たないsemantic normalization
#include "solve_normalization.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/array_utils.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"

#include <vector>

namespace mmcal::solver {
namespace {

using expression::Expr;

[[nodiscard]] Expr canonicalizeBuiltinHeads(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (expression.isCall()) {
        const auto& call = expression.asCall();
        std::vector<Expr> arguments;
        arguments.reserve(call.arguments.size());
        bool changed = false;
        for (const Expr& argument : call.arguments) {
            Expr normalized = canonicalizeBuiltinHeads(argument, builtins);
            changed = changed || !(normalized == argument);
            arguments.push_back(std::move(normalized));
        }

        expression::Symbol head = call.head;
        if (const auto* definition = builtins.find(call.head)) {
            const expression::Symbol& canonical = builtins.symbol(definition->id);
            if (!head.sameIdentity(canonical)) {
                head = canonical;
                changed = true;
            }
        }

        if (!changed)
            return expression;
        if (head.sameIdentity(call.head))
            return Expr::rebuildCall(call, std::move(arguments));
        return Expr::call(std::move(head), std::move(arguments));
    }

    if (expression.isArray()) {
        const auto& array = expression.asArray();
        std::vector<Expr> elements;
        elements.reserve(array.size());
        bool changed = false;
        for (std::size_t i = 0; i < array.size(); ++i) {
            Expr normalized = canonicalizeBuiltinHeads(array.element(i), builtins);
            changed = changed || !(normalized == array.element(i));
            elements.push_back(std::move(normalized));
        }
        return changed ? Expr::array(array.shape, std::move(elements)) : expression;
    }

    if (expression.isList()) {
        const auto& list = expression.asList();
        std::vector<Expr> elements;
        elements.reserve(list.elements.size());
        bool changed = false;
        for (const Expr& element : list.elements) {
            Expr normalized = canonicalizeBuiltinHeads(element, builtins);
            changed = changed || !(normalized == element);
            elements.push_back(std::move(normalized));
        }
        return changed ? expression::braceValue(std::move(elements)) : expression;
    }

    return expression;
}

} // namespace

expression::Expr normalizeForSolve(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr normalized = canonicalizeBuiltinHeads(expression, builtins);
    return simplification::Simplifier{}.simplify(
        normalized,
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
}

} // namespace mmcal::solver
