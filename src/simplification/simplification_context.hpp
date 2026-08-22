#pragma once

#include "evaluation/builtin_registry.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/math_registry.hpp"
#include "symbols/symbol_registry.hpp"

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
    evaluation::EvaluationBudget* budget = nullptr;
    // Session固有のprotected exceptional valueを返す必要があるtop-level評価で指定する。
    // 純粋な内部変形で未指定なら，不定形を勝手に別SymbolTableの値へ置換しない。
    const symbols::SymbolRegistry* predefinedSymbols = nullptr;
    // 恒等式証明など，入力式が定義される点上だけで値の一致を調べる内部用途。
    // 通常のsimplify/fullSimplifyではfalseのままとし，定義域を広げる簡約を禁止する。
    bool assumeExpressionsDefined = false;

    [[nodiscard]] mathematics::KnowledgeContext knowledge() const {
        return mathematics::KnowledgeContext{builtins, mathematics, assumptions};
    }
};

} // namespace mmcal::simplification
