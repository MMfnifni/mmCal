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

// 利用者が明示的にsimplify/fullSimplifyを要求した場合だけ使う強めの後処理。
// exact有理係数の線形結合を平坦化し，同一項を結合する。自動簡約の契約には混ぜない。
[[nodiscard]] expression::Expr simplifyExplicitLinearCombination(
    const expression::Expr& expression,
    const SimplificationContext& context);

} // namespace mmcal::simplification
