#pragma once

#include "plot_program.hpp"
#include "plot_sampling.hpp"
#include "plot_view_transform.hpp"

#include "mathematics/angle.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace mmcal::plot {

struct CurveIntersectionOptions final {
    std::size_t refinementIterations = 36;
    std::size_t tangencyRefinementIterations = 32;
    std::size_t maxTangencyCandidates = 512;
    std::size_t maxIntersections = 256;
    double minimumRefinementSpanMm = 0.10;
    double maximumSnapDistanceMm = 0.10;
};

enum class CurveIntersectionStatus {
    Success,
    InvalidInput,
    ProgramInitializationFailed,
    ResourceLimit
};

struct CurveIntersectionStatistics final {
    std::size_t evaluations = 0;
    std::size_t intersections = 0;
    std::size_t tangentialIntersections = 0;
    std::size_t tangencyCandidates = 0;
    std::size_t rejectedTangencyCandidates = 0;
    std::size_t rejectedUnsafeSnaps = 0;
    std::size_t coincidentAffinePairs = 0;
    std::size_t truncatedPairs = 0;
};

struct CurveIntersectionResult final {
    CurveIntersectionStatus status = CurveIntersectionStatus::InvalidInput;
    std::optional<std::vector<SampledCurve>> curves;
    CurveIntersectionStatistics statistics;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == CurveIntersectionStatus::Success && curves.has_value();
    }
};

// 同一x軸上の複数函数について，sampleでbracketできる横切り交点に加え，
// |f-g|の局所極小から接触交点候補をbounded refinementする。
// 確定した交点は両曲線へ共通座標のsemantic vertexとして挿入する。
// anchorは描画補助情報なので，高密度交点では上限到達後にpairだけ打ち切り，曲線描画自体は失敗させない。
[[nodiscard]] CurveIntersectionResult refineCurveIntersections(
    const std::vector<PlotProgram>& programs,
    const std::vector<SampledCurve>& curves,
    mathematics::AngleSemantics angles,
    const CurveIntersectionOptions& options = {},
    const PlotViewTransform* transform = nullptr);

} // namespace mmcal::plot
