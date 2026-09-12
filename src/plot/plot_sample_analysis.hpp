#pragma once

#include "plot_analysis.hpp"
#include "plot_sampling.hpp"

#include "numeric/big_float.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace mmcal::plot {

// 粗sampleから得られる有限値の素朴なrange。外れ値除去等は後段のRangeEstimatorで行う。
struct PlotFiniteRangeCandidate final {
    numeric::BigFloat minimum;
    numeric::BigFloat maximum;
    std::size_t finiteSampleCount = 0;
};

enum class PlotSuspicionKind {
    NonFiniteSample,
    RapidVariation,
    Oscillatory,
    UnresolvedBoundary
};

// 粗sample列のうち，後段で追加評価すべき範囲をsample indexで保持する。
// boundary由来の場合も，exact endpoint自体はsourceIntervalから復元できる。
struct PlotSuspiciousSpan final {
    std::size_t segmentIndex = 0;
    std::size_t firstSampleIndex = 0;
    std::size_t lastSampleIndex = 0;
    PlotSuspicionKind kind = PlotSuspicionKind::NonFiniteSample;
    PlotNumericStatus numericStatus = PlotNumericStatus::Finite;
};

struct PlotCoarseInspection final {
    std::optional<PlotFiniteRangeCandidate> finiteRange;
    std::vector<PlotSuspiciousSpan> suspiciousSpans;
};

struct PlotCoarseInspectionOptions final {
    // first differenceの向きがこの回数以上反転したsegmentをoscillatory candidateとする。
    std::size_t minimumOscillationTurns = 3;
    std::size_t maxSuspiciousSpans = 1024;
};

// 粗sampling結果を，range候補とadaptive refinement候補へ整理する。
// この層は数学的な特異点を新規に証明せず，PlotAnalysisの既知landmarkと数値列だけを使う。
[[nodiscard]] PlotCoarseInspection inspectCoarseSamples(
    const PlotAnalysis& analysis,
    const SampledCurve& curve,
    const PlotCoarseInspectionOptions& options = {});

} // namespace mmcal::plot
