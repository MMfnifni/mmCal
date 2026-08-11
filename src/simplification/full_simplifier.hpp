#pragma once

#include "expression/expr.hpp"
#include "simplification_context.hpp"

#include <cstddef>

namespace mmcal::simplification {

struct FullSimplificationOptions final {
    std::size_t maximumCandidates = 96;
};

// canonical Simplifierに加え、複数のexact同値形を探索して最小コストを選ぶ。
[[nodiscard]] expression::Expr fullSimplify(
    const expression::Expr& expression,
    const SimplificationContext& context,
    FullSimplificationOptions options = {});

} // namespace mmcal::simplification
