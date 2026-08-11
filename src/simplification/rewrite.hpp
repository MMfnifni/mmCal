#pragma once

#include "expression/expr.hpp"

namespace mmcal::simplification {

// 自動簡約の適用順序。後段ほど数学的な前提を多く必要とする。
enum class RewritePhase {
    Canonical,
    Exact,
    Conditional
};

struct RewriteResult final {
    expression::Expr expression;
    bool changed = false;
    RewritePhase phase = RewritePhase::Canonical;
};

} // namespace mmcal::simplification
