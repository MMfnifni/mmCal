#pragma once

#include "expression/expr.hpp"

#include <cstddef>

namespace mmcal::simplification {

// FullSimplifyが複数の同値候補を比較するときのための、意味を持たない構造コスト。
// では探索は行わず、将来の評価尺度だけ先に独立させる。
struct ExpressionCost final {
    std::size_t nodes = 0;
    std::size_t leaves = 0;
    std::size_t depth = 0;

    [[nodiscard]] bool operator==(const ExpressionCost&) const = default;
};

[[nodiscard]] ExpressionCost measureExpressionCost(const expression::Expr& expression);

} // namespace mmcal::simplification
