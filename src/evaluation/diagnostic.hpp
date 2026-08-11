// Warning・Info診断情報
#pragma once

#include "expression/expr.hpp"

#include <optional>
#include <string>

namespace mmcal::evaluation {

enum class DiagnosticSeverity {
    Info,
    Warning
};

// 評価は成功したが、Frontendへ補足通知すべき情報を保持する。Warningは未評価・未解決、Infoは再定義など正常な状態変更に使う。
struct EvaluationDiagnostic final {
    DiagnosticSeverity severity = DiagnosticSeverity::Warning;
    std::string code;
    std::string message;
    std::optional<expression::Expr> previousExpression;

    [[nodiscard]] bool operator==(const EvaluationDiagnostic&) const = default;
};

} // namespace mmcal::evaluation
