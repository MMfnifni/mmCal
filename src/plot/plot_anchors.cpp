// PlotScene手前でsampling vertexを共有semantic anchorへ昇格する。
#include "plot_anchors.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace mmcal::plot {
namespace {

[[nodiscard]] PlotAnchorKind anchorKinds(PlotPointTag tags) noexcept {
    PlotAnchorKind kinds = PlotAnchorKind::None;
    if (hasPlotPointTag(tags, PlotPointTag::XAxisIntercept))
        kinds |= PlotAnchorKind::XAxisIntercept;
    if (hasPlotPointTag(tags, PlotPointTag::YAxisIntercept))
        kinds |= PlotAnchorKind::YAxisIntercept;
    if (hasPlotPointTag(tags, PlotPointTag::LocalExtremum))
        kinds |= PlotAnchorKind::LocalExtremum;
    if (hasPlotPointTag(tags, PlotPointTag::CurveIntersection))
        kinds |= PlotAnchorKind::CurveIntersection;
    return kinds;
}

struct AnchorCandidate final {
    numeric::BigFloat x;
    numeric::BigFloat y;
    PlotAnchorKind kinds = PlotAnchorKind::None;
    PlotAnchorVertexRef vertex;
};

} // namespace

PlotAnchorSet collectPlotAnchors(const std::vector<SampledCurve>& curves) {
    std::vector<AnchorCandidate> candidates;
    for (std::size_t curveIndex = 0; curveIndex < curves.size(); ++curveIndex) {
        for (std::size_t segmentIndex = 0; segmentIndex < curves[curveIndex].segments.size(); ++segmentIndex) {
            const auto& segment = curves[curveIndex].segments[segmentIndex];
            for (std::size_t sampleIndex = 0; sampleIndex < segment.samples.size(); ++sampleIndex) {
                const auto& sample = segment.samples[sampleIndex];
                if (!sample.finite())
                    continue;
                const auto kinds = anchorKinds(sample.tags);
                if (kinds == PlotAnchorKind::None)
                    continue;
                candidates.push_back(AnchorCandidate{
                    sample.x, sample.y, kinds,
                    PlotAnchorVertexRef{curveIndex, segmentIndex, sampleIndex}});
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const AnchorCandidate& lhs, const AnchorCandidate& rhs) {
        if (lhs.x != rhs.x)
            return lhs.x < rhs.x;
        if (lhs.y != rhs.y)
            return lhs.y < rhs.y;
        if (lhs.vertex.curveIndex != rhs.vertex.curveIndex)
            return lhs.vertex.curveIndex < rhs.vertex.curveIndex;
        if (lhs.vertex.segmentIndex != rhs.vertex.segmentIndex)
            return lhs.vertex.segmentIndex < rhs.vertex.segmentIndex;
        return lhs.vertex.sampleIndex < rhs.vertex.sampleIndex;
    });

    PlotAnchorSet result;
    for (const auto& candidate : candidates) {
        if (result.anchors.empty()
            || result.anchors.back().x != candidate.x
            || result.anchors.back().y != candidate.y) {
            PlotAnchor anchor;
            anchor.id = static_cast<PlotAnchorId>(result.anchors.size() + 1);
            anchor.kinds = candidate.kinds;
            anchor.x = candidate.x;
            anchor.y = candidate.y;
            anchor.curveIndices.push_back(candidate.vertex.curveIndex);
            anchor.vertices.push_back(candidate.vertex);
            result.anchors.push_back(std::move(anchor));
            continue;
        }

        auto& anchor = result.anchors.back();
        anchor.kinds |= candidate.kinds;
        if (anchor.curveIndices.empty() || anchor.curveIndices.back() != candidate.vertex.curveIndex)
            anchor.curveIndices.push_back(candidate.vertex.curveIndex);
        anchor.vertices.push_back(candidate.vertex);
    }
    return result;
}

} // namespace mmcal::plot
