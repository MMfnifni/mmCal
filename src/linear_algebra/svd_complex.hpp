#pragma once

#include "approximation/certified_evaluator.hpp"
#include "expression/expr.hpp"

#include <cstddef>
#include <optional>

namespace mmcal::linear_algebra::detail {

// 一般複素行列のreduced numerical SVDを一つの作業precisionで試行する。
// 成功条件はA=U S V^Hの再構成とU/Vの直交性が要求桁より厳しい区間監査を通ること。
[[nodiscard]] std::optional<expression::Expr> approximateComplexSvdAtPrecision(
    const expression::ArrayExpr& source,
    std::size_t bits,
    std::size_t digits,
    const approximation::CertifiedEvaluator& certified);

} // namespace mmcal::linear_algebra::detail
