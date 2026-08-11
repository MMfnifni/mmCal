#pragma once

#include "assumption_set.hpp"
#include "expression/expr.hpp"

#include <optional>

namespace mmcal::evaluation {
class BuiltinRegistry;
}

namespace mmcal::mathematics {
class MathRegistry;

// scalarな数学式が値を持つために必要な条件を収集する。
// 未対応headや現在のPredicate表現で安全に記述できないdomainはnulloptを返す。
[[nodiscard]] std::optional<AssumptionSet> expressionDomainConditions(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics);

} // namespace mmcal::mathematics
