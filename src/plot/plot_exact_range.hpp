#pragma once

#include "plot_exact_geometry.hpp"
#include "plot_range_estimator.hpp"

namespace mmcal::plot {

// exact geometryからsampling密度に依存しないAutomatic rangeを作る。
// EllipticArcで一般の極値phaseをexactに決められない場合だけ，完全楕円の
// 解析的bboxへ安全に包絡する。その場合tight=falseで区別する。
struct ExactCurveRange2D final {
    CurveRangeEstimate1D x;
    CurveRangeEstimate1D y;
    bool tight = true;
};

[[nodiscard]] std::optional<ExactCurveRange2D> estimateExactCurveRange(
    const ExactCurveGeometry& geometry,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options = {});

} // namespace mmcal::plot
