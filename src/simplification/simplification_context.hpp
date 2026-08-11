#pragma once

#include "evaluation/builtin_registry.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/math_registry.hpp"

#include <cstddef>

namespace mmcal::simplification {

// Simplifierへ渡す読み取り専用の数学的文脈。
// 将来 simplify[..., assumptions] / Refine / Solver が同じKnowledgeContextを共有できる。
struct SimplificationContext final {
    const evaluation::BuiltinRegistry& builtins;
    const mathematics::MathRegistry& mathematics;
    const mathematics::AngleSemantics& angleSemantics;
    mathematics::AssumptionSet assumptions{};
    std::size_t maximumPasses = 32;

    [[nodiscard]] mathematics::KnowledgeContext knowledge() const {
        return mathematics::KnowledgeContext{builtins, mathematics, assumptions};
    }
};

} // namespace mmcal::simplification
