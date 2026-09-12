#pragma once

#include "plot_adaptive_sampling.hpp"
#include "plot_analysis.hpp"
#include "plot_axis_layout.hpp"
#include "plot_graphics_lowering.hpp"
#include "plot_exact_geometry.hpp"
#include "plot_exact_range.hpp"
#include "plot_intersections.hpp"
#include "plot_periodicity.hpp"
#include "plot_range_estimator.hpp"
#include "plot_request.hpp"
#include "plot_scene.hpp"
#include "plot_view_transform.hpp"

#include "graphics/graphics_scene.hpp"
#include "mathematics/assumption_set.hpp"

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace mmcal::plot {

struct PlotPipelineOptions final {
    CoarseSamplingOptions coarseSampling;
    PlotRangeEstimatorOptions rangeEstimator;
    AdaptiveSamplingOptions adaptiveSampling;
    CurveIntersectionOptions intersections;
    PlotAxisLayoutOptions axisLayout;
    PlotGraphicsLoweringOptions graphicsLowering;
    PlotViewportMm viewport;
    // toNormal等，exact primitive自体ではなく最終sample座標列が必要なconsumer向け。
    // 通常描画ではfalseのままにし，Line/Bezier/Ellipse/Arcのexact fast pathを保持する。
    bool forceSampledGeometry = false;
};

enum class PlotPipelineStatus {
    Success,
    EmptyRequest,
    InvalidOptions,
    UnsafeDiscontinuity,
    CompileFailed,
    SamplingFailed,
    RangeEstimationFailed,
    ViewTransformFailed,
    AdaptiveSamplingFailed,
    IntersectionFailed,
    SceneBuildFailed,
    AxisLayoutFailed,
    GraphicsLoweringFailed,
    SvgRenderFailed,
    DomainIntersectionFailed
};

struct PlotPipelineOutput final {
    std::vector<PlotAnalysis> analyses;
    std::vector<PlotProgram> programs;
    std::vector<SampledCurve> curves;
    PlotRangeEstimate range;
    PlotViewTransform transform;
    PlotAxesLayout axes;
    PlotScene plotScene;
    graphics::GraphicsScene graphicsScene;
};

struct PlotPipelineResult final {
    PlotPipelineStatus status = PlotPipelineStatus::EmptyRequest;
    std::size_t failedCurveIndex = std::numeric_limits<std::size_t>::max();
    std::optional<PlotPipelineOutput> output;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PlotPipelineStatus::Success && output.has_value();
    }
};


struct ParametricPlotPipelineOutput final {
    std::vector<PlotAnalysis> xAnalyses;
    std::vector<PlotAnalysis> yAnalyses;
    std::vector<std::pair<PlotProgram, PlotProgram>> programs;
    std::vector<std::optional<ParametricPeriodReduction>> periodReductions;
    std::vector<std::optional<ExactCurveGeometry>> exactGeometries;
    std::vector<std::optional<ExactCurveRange2D>> exactRanges;
    std::vector<SampledCurve2D> curves;
    CurveViewRange2D range;
    PlotViewTransform transform;
    PlotAxesLayout axes;
    PlotScene plotScene;
    graphics::GraphicsScene graphicsScene;
};

struct ParametricPlotPipelineResult final {
    PlotPipelineStatus status = PlotPipelineStatus::EmptyRequest;
    std::size_t failedCurveIndex = std::numeric_limits<std::size_t>::max();
    std::optional<ParametricPlotPipelineOutput> output;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PlotPipelineStatus::Success && output.has_value();
    }
};

struct ParametricPlotSvgPipelineResult final {
    PlotPipelineStatus status = PlotPipelineStatus::EmptyRequest;
    std::size_t failedCurveIndex = std::numeric_limits<std::size_t>::max();
    std::optional<ParametricPlotPipelineOutput> output;
    std::optional<std::string> svg;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PlotPipelineStatus::Success && output.has_value() && svg.has_value();
    }
};

struct PlotSvgPipelineResult final {
    PlotPipelineStatus status = PlotPipelineStatus::EmptyRequest;
    std::size_t failedCurveIndex = std::numeric_limits<std::size_t>::max();
    std::optional<PlotPipelineOutput> output;
    std::optional<std::string> svg;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PlotPipelineStatus::Success && output.has_value() && svg.has_value();
    }
};

// Plotの内部工程を一つに束ねる。各段の数学的責務は既存層へ委譲し，
// この層自身ではdomain推測や追加simplificationを行わない。
[[nodiscard]] PlotPipelineResult buildPlotPipeline(
    const std::vector<PlotRequest>& requests,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

[[nodiscard]] PlotPipelineResult buildPlotPipeline(
    const PlotRequestSet& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

[[nodiscard]] PlotPipelineResult buildPlotPipeline(
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

// buildPlotPipelineの結果をそのままSVG backendへ送るCLI向け入口。
[[nodiscard]] PlotSvgPipelineResult renderPlotSvg(
    const std::vector<PlotRequest>& requests,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

[[nodiscard]] PlotSvgPipelineResult renderPlotSvg(
    const PlotRequestSet& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

[[nodiscard]] PlotSvgPipelineResult renderPlotSvg(
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});


// ParametricPlotはx/yのscalar analysisを独立に行い，domainの積集合上だけを
// 共有2D curve pipelineへ流す。未証明domainを跨いだ補間は行わない。
[[nodiscard]] ParametricPlotPipelineResult buildParametricPlotPipeline(
    const std::vector<ParametricCurveRequest>& requests,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

[[nodiscard]] ParametricPlotPipelineResult buildParametricPlotPipeline(
    const ParametricPlotRequestSet& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

[[nodiscard]] ParametricPlotPipelineResult buildParametricPlotPipeline(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

[[nodiscard]] ParametricPlotSvgPipelineResult renderParametricPlotSvg(
    const std::vector<ParametricCurveRequest>& requests,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

[[nodiscard]] ParametricPlotSvgPipelineResult renderParametricPlotSvg(
    const ParametricPlotRequestSet& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

[[nodiscard]] ParametricPlotSvgPipelineResult renderParametricPlotSvg(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {},
    const PlotPipelineOptions& options = {});

} // namespace mmcal::plot
