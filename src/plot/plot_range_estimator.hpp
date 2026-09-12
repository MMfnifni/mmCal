#pragma once

#include "plot_analysis.hpp"
#include "plot_sample_analysis.hpp"
#include "plot_sampling.hpp"

#include "numeric/big_float.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace mmcal::plot {

// 粗sampleから得るAutomatic PlotRange用の仮range。
// observedは全finite sample，dataは既知asymptote近傍等を除いたrange，viewはpadding後を表す。
struct CurveRangeEstimate1D final {
    numeric::BigFloat observedMinimum;
    numeric::BigFloat observedMaximum;
    numeric::BigFloat dataMinimum;
    numeric::BigFloat dataMaximum;
    numeric::BigFloat viewMinimum;
    numeric::BigFloat viewMaximum;
    std::size_t includedFiniteSamples = 0;
    std::size_t excludedFiniteSamples = 0;
    bool boundaryTrimmed = false;
    bool flatExpanded = false;
};

// 既存Plot APIの互換名。range推定自体は1次元量としてx/y共通で扱う。
using PlotRangeEstimate = CurveRangeEstimate1D;

enum class CurveCoordinate2D {
    X,
    Y
};

struct CurveViewRange2D final {
    CurveRangeEstimate1D x;
    CurveRangeEstimate1D y;
};

struct PlotRangeEstimatorOptions final {
    // proven pole/asymptoteのopen endpointから何sample分をrange推定だけで無視するか。
    // 描画sample自体は捨てず，後段adaptive refinement用に保持する。
    std::size_t asymptoteGuardSamples = 2;
    // 未解析domain内にnonfinite sampleが出た場合，その周囲をrange推定から除く。
    std::size_t nonFiniteGuardSamples = 1;
    // trimming後にこれ未満しか残らなければ，全finite sampleへ安全にfallbackする。
    std::size_t minimumIncludedFiniteSamples = 3;
    // proven asymptoteのtailが固定guardより長い場合だけ，内部sampleの代表scaleに
    // 対してこの倍率を超える連続tailを追加でrange推定から外す。
    // 一般の急増函数には適用せず，既知asymptote境界に限定する。
    std::size_t asymptoteTailScaleFactor = 8;
    // view rangeへ加える余白。既定はdata spanの5%。
    std::size_t paddingNumerator = 1;
    std::size_t paddingDenominator = 20;
};


// 既知の解析的extentから，通常のAutomatic rangeと同じpadding規則でrangeを作る。
[[nodiscard]] std::optional<CurveRangeEstimate1D> makeCurveRangeEstimate(
    const numeric::BigFloat& minimum,
    const numeric::BigFloat& maximum,
    std::size_t precisionBits,
    std::size_t contributingPoints = 2,
    const PlotRangeEstimatorOptions& options = {});

// 2D sampled curveの指定coordinateから汎用1D rangeを作る。
[[nodiscard]] std::optional<CurveRangeEstimate1D> estimateCurveCoordinateRange(
    const SampledCurve2D& curve,
    CurveCoordinate2D coordinate,
    const PlotRangeEstimatorOptions& options = {});

// coordinate自身のPlotAnalysisを併用し，proven pole/asymptote境界のtailだけを
// Automatic rangeから除外する。ParametricPlotのx/y rangeで使用する。
[[nodiscard]] std::optional<CurveRangeEstimate1D> estimateAnalyzedCurveCoordinateRange(
    const PlotAnalysis& analysis,
    const SampledCurve2D& curve,
    CurveCoordinate2D coordinate,
    const PlotRangeEstimatorOptions& options = {});

// 1D range群をdata extentで統合する共通処理。
[[nodiscard]] std::optional<CurveRangeEstimate1D> combineCurveRangeEstimates(
    const std::vector<CurveRangeEstimate1D>& estimates,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options = {});

// coarse inspectionを使って仮viewport rangeを決める。
// percentileだけで値を切らず，既知asymptote/nonfinite近傍だけを局所的に除外する。
[[nodiscard]] std::optional<PlotRangeEstimate> estimatePlotRange(
    const PlotAnalysis& analysis,
    const SampledCurve& curve,
    const PlotCoarseInspection& inspection,
    const PlotRangeEstimatorOptions& options = {});

// 複数curveのrangeをdata extentで統合し，paddingは統合後に一度だけ適用する。
// flat curveごとの局所paddingがmulti-curve viewportを不必要に膨らませるのを避ける。
[[nodiscard]] std::optional<PlotRangeEstimate> combinePlotRangeEstimates(
    const std::vector<PlotRangeEstimate>& estimates,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options = {});

} // namespace mmcal::plot
