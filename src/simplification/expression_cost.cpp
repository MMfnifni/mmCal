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
            const auto& array = current.expression.asArray();
            const auto entries = array.expressionEntries();
            const std::size_t packedLeaves = array.size() - entries.size();
            result.nodes += packedLeaves;
            result.leaves += packedLeaves;
            if (!array.empty())
                result.depth = std::max(result.depth, current.depth + 1);
            for (const auto& entry : entries)
                stack.push_back(Entry{entry.expression, current.depth + 1});
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
