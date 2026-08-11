#pragma once

#include "expression/expr.hpp"

#include <string>

namespace mmcal::simplification {

// 可換演算のcanonical化専用の決定的な構造順序。
// 数学上の大小関係ではなく、同じExpr集合を常に同じAST順へ並べるために使う。
[[nodiscard]] std::string expressionOrderKey(const expression::Expr& expression);

struct ExpressionLess final {
    [[nodiscard]] bool operator()(
        const expression::Expr& lhs,
        const expression::Expr& rhs) const;
};

} // namespace mmcal::simplification
