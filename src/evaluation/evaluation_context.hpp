#pragma once

#include "diagnostic.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"

#include <optional>
#include <span>
#include <vector>

namespace mmcal::evaluation {

// 1回の評価要求にだけ属する外部状態をまとめる。
// Evaluator本体の永続状態(Environment/registry等)と、履歴・diagnostic・session flagを分離する。
struct EvaluationContext final {
    // % / %% と Out[-n] 用。成功した出力だけを相対順に保持する。
    std::span<const expression::Expr> history;
    // 正の In[n] / Out[n] 用。index n-1 が画面上の絶対入力番号nに対応する。
    // 負のIn参照は現在入力slotを除いたinputsから相対参照する。
    std::span<const std::optional<expression::Expr>> inputs;
    std::span<const std::optional<expression::Expr>> outputs;
    std::vector<EvaluationDiagnostic>* diagnostics = nullptr;
    bool* exitRequested = nullptr;
    bool* clearRequested = nullptr;
    bool* definitionsChanged = nullptr;
    mathematics::AngleSemantics* angleSemantics = nullptr;
};

} // namespace mmcal::evaluation
