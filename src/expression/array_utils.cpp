// Array shapeと評価後flattenの共通処理
#include "array_utils.hpp"

#include <limits>
#include <stdexcept>
#include <utility>

namespace mmcal::expression {

std::size_t arrayElementCount(std::span<const std::size_t> shape) {
    if (shape.empty())
        throw std::invalid_argument("Array shape must contain at least one dimension");

    std::size_t count = 1;
    for (const std::size_t dimension : shape) {
        if (dimension == 0)
            return 0;
        if (count > std::numeric_limits<std::size_t>::max() / dimension)
            throw std::length_error("Array element count exceeds the size_t range");
        count *= dimension;
    }
    return count;
}

Expr rebuildEvaluatedArray(
    std::vector<std::size_t> outerShape,
    std::vector<Expr> evaluatedElements) {
    // Expr::arrayをArray生成の唯一の正規化境界とし、Evaluator固有の別規則を持たない。
    return Expr::array(std::move(outerShape), std::move(evaluatedElements));
}

Expr braceValue(std::vector<Expr> elements) {
    const std::size_t count = elements.size();
    if (elements.empty())
        return Expr::array({0}, {});

    // 例外を矩形性判定へ使わず，Arrayへ昇格できる場合だけ一度moveする。
    // 大きな{Q,R}等で不要なExpr vector copyを発生させない。
    const bool firstArray = elements.front().isArray();
    if (elements.front().isList())
        return Expr::list(std::move(elements));
    if (firstArray) {
        const auto& shape = elements.front().asArray().shape;
        for (std::size_t i = 1; i < elements.size(); ++i)
            if (!elements[i].isArray() || elements[i].asArray().shape != shape)
                return Expr::list(std::move(elements));
    }
    else {
        for (std::size_t i = 1; i < elements.size(); ++i)
            if (elements[i].isArray() || elements[i].isList())
                return Expr::list(std::move(elements));
    }
    return Expr::array({count}, std::move(elements));
}

std::vector<std::size_t> commonBraceDimensions(const Expr& expression) {
    if (expression.isArray())
        return expression.asArray().shape;
    if (!expression.isList())
        return {};

    const auto& list = expression.asList();
    std::vector<std::size_t> result{list.elements.size()};
    if (list.elements.empty())
        return result;

    std::vector<std::size_t> child = commonBraceDimensions(list.elements.front());
    if (child.empty())
        return result;
    for (std::size_t i = 1; i < list.elements.size(); ++i)
        if (commonBraceDimensions(list.elements[i]) != child)
            return result;
    result.insert(result.end(), child.begin(), child.end());
    return result;
}

bool braceLiteralPreservesShape(std::span<const std::size_t> shape) noexcept {
    for (std::size_t i = 0; i < shape.size(); ++i)
        if (shape[i] == 0)
            return i + 1 == shape.size();
    return true;
}

} // namespace mmcal::expression
