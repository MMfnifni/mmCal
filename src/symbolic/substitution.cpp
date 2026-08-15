// 式中Symbolの置換
#include "substitution.hpp"
#include "expression/array_utils.hpp"

#include <utility>
#include <vector>

namespace mmcal::symbolic {

expression::Expr substituteSymbol(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const expression::Expr& value) {
    if (expression.isSymbol() && expression.asSymbol().sameIdentity(variable))
        return value;

    if (expression.isCall()) {
        std::vector<expression::Expr> arguments;
        arguments.reserve(expression.asCall().arguments.size());
        for (const expression::Expr& argument : expression.asCall().arguments)
            arguments.push_back(substituteSymbol(argument, variable, value));
        return expression::Expr::call(expression.asCall().head, std::move(arguments));
    }

    if (expression.isArray()) {
        const auto& array = expression.asArray();
        const auto entries = array.expressionEntries();
        if (entries.empty())
            return expression;
        std::vector<std::size_t> indices;
        std::vector<expression::Expr> elements;
        indices.reserve(entries.size());
        elements.reserve(entries.size());
        for (const auto& entry : entries) {
            indices.push_back(entry.index);
            elements.push_back(substituteSymbol(entry.expression, variable, value));
        }
        return expression::Expr::array(array.replacedExpressions(indices, std::move(elements)));
    }

    if (expression.isList()) {
        std::vector<expression::Expr> elements;
        elements.reserve(expression.asList().elements.size());
        for (const expression::Expr& element : expression.asList().elements)
            elements.push_back(substituteSymbol(element, variable, value));
        return expression::braceValue(std::move(elements));
    }

    return expression;
}

} // namespace mmcal::symbolic
