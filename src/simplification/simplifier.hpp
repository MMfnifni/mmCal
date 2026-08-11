#pragma once

#include "expression/expr.hpp"
#include "simplification_context.hpp"

namespace mmcal::simplification {

// 副作用を持たない保守的な自動簡約器。
// 「成立条件を証明できない変形は行わない」を最優先の契約とする。
class Simplifier final {
public:
    [[nodiscard]] expression::Expr simplify(
        const expression::Expr& expression,
        const SimplificationContext& context) const;
};

} // namespace mmcal::simplification
