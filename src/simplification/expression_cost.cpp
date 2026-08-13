// 式の複雑度評価
#include "expression_cost.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace mmcal::simplification {

ExpressionCost measureExpressionCost(const expression::Expr& expression) {
    struct Entry final {
        expression::Expr expression;
        std::size_t depth = 1;
    };

    ExpressionCost result;
    std::vector<Entry> stack{{expression, 1}};

    while (!stack.empty()) {
        Entry current = std::move(stack.back());
        stack.pop_back();
        ++result.nodes;
        result.depth = std::max(result.depth, current.depth);

        if (current.expression.isCall()) {
            for (const auto& argument : current.expression.asCall().arguments)
                stack.push_back(Entry{argument, current.depth + 1});
            continue;
        }
        if (current.expression.isArray()) {
            for (const auto& element : current.expression.asArray().elements)
                stack.push_back(Entry{element, current.depth + 1});
            continue;
        }
        if (current.expression.isList()) {
            for (const auto& element : current.expression.asList().elements)
                stack.push_back(Entry{element, current.depth + 1});
            continue;
        }

        ++result.leaves;
    }

    return result;
}

} // namespace mmcal::simplification
