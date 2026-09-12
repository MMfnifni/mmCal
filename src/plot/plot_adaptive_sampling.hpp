#pragma once

#include "plot_program.hpp"
#include "plot_sampling.hpp"
#include "plot_view_transform.hpp"

#include "mathematics/angle.hpp"

#include <cstddef>
#include <optional>

namespace mmcal::plot {

struct AdaptiveSamplingOptions final {
    std::size_t maxRecursion = 10;
    std::size_t maxTotalSamples = 16384;
    double chordToleranceMm = 0.10;
    double minimumSpanMm = 0.10;
    // open domain端はendpoint自体を評価できないため，通常のcoarse sampleだけでは
    // tan[x]^(-1/2)の有限極限やpole近傍が端まで届かない。最終viewport確定後に
    // この物理距離まで内側へ専用sampleを一つ置く。
    double openBoundaryApproachMm = 0.005;
    std::size_t rootRefinementIterations = 32;
    std::size_t extremumRefinementIterations = 12;
    bool preserveAxisIntersections = true;
    bool preserveLocalExtrema = true;
};

enum class AdaptiveSamplingStatus {
    Success,
    InvalidOptions,
    ProgramInitializationFailed,
    ResourceLimit
};

struct AdaptiveSamplingStatistics final {
    std::size_t evaluations = 0;
    std::size_t insertedSamples = 0;
    std::size_t xAxisIntercepts = 0;
    std::size_t yAxisIntercepts = 0;
    std::size_t localExtrema = 0;
    std::size_t recursionLimitHits = 0;
};

struct AdaptiveSamplingResult final {
    AdaptiveSamplingStatus status = AdaptiveSamplingStatus::InvalidOptions;
    std::optional<SampledCurve> curve;
    AdaptiveSamplingStatistics statistics;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == AdaptiveSamplingStatus::Success && curve.has_value();
    }
};

// 粗sampleをplot-layout-space(mm)のchord errorで細分化する。
// 意味のあるplot pointはsample vertexへ強制挿入し，後のGraphics Pathへそのまま残せる。
[[nodiscard]] AdaptiveSamplingResult refinePlotSamples(
    const PlotProgram& program,
    const SampledCurve& coarse,
    const PlotViewTransform& transform,
    mathematics::AngleSemantics angles,
    const AdaptiveSamplingOptions& options = {});

[[nodiscard]] AdaptiveSamplingResult refineParametricPlotSamples(
    const PlotProgram& xProgram,
    const PlotProgram& yProgram,
    const SampledCurve2D& coarse,
    const PlotViewTransform& transform,
    mathematics::AngleSemantics angles,
    const AdaptiveSamplingOptions& options = {});

} // namespace mmcal::plot
