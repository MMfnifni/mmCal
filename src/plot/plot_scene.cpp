#include "plot_scene.hpp"

#include <utility>
#include <vector>

namespace mmcal::plot {
namespace {

struct SceneVertexLocation final {
    bool valid = false;
    std::size_t segmentIndex = 0;
    std::size_t vertexIndex = 0;
};

[[nodiscard]] PlotSceneCurveId curveId(std::size_t curveIndex) noexcept {
    return static_cast<PlotSceneCurveId>(curveIndex + 1);
}

} // namespace

PlotSceneBuildResult buildPlotScene(
    const std::vector<SampledCurve>& curves,
    std::optional<PlotAnchorSet> anchors) {
    PlotSceneBuildResult result;
    if (!anchors)
        anchors = collectPlotAnchors(curves);

    // numeric holeを消して両側を再接続しないよう，有限sampleの連続runごとにScene segmentを作る。
    std::vector<std::vector<std::vector<SceneVertexLocation>>> locations(curves.size());
    PlotScene scene;
    scene.curves.reserve(curves.size());
    for (std::size_t curveIndex = 0; curveIndex < curves.size(); ++curveIndex) {
        PlotSceneCurve curve;
        curve.id = curveId(curveIndex);
        locations[curveIndex].resize(curves[curveIndex].segments.size());
        for (std::size_t sourceSegmentIndex = 0;
             sourceSegmentIndex < curves[curveIndex].segments.size(); ++sourceSegmentIndex) {
            const auto& sourceSegment = curves[curveIndex].segments[sourceSegmentIndex];
            auto& segmentLocations = locations[curveIndex][sourceSegmentIndex];
            segmentLocations.resize(sourceSegment.samples.size());

            std::optional<std::size_t> activeSceneSegment;
            for (std::size_t sampleIndex = 0; sampleIndex < sourceSegment.samples.size(); ++sampleIndex) {
                const auto& sample = sourceSegment.samples[sampleIndex];
                if (!sample.finite()) {
                    activeSceneSegment.reset();
                    continue;
                }

                if (!activeSceneSegment) {
                    curve.segments.push_back(PlotSceneSegment{
                        sourceSegment.sourceInterval, sourceSegment.geometryKind, {},
                        sourceSegment.quadraticControlPoint,
                        sourceSegment.bezierControlPoints,
                        sourceSegment.ellipseGeometry,
                        sourceSegment.ellipticArcGeometry,
                        sourceSegment.markLowerEndpoint && sampleIndex == 0,
                        false,
                        sampleIndex == 0 ? sourceSegment.lowerEndpointPoint : std::nullopt,
                        std::nullopt});
                    activeSceneSegment = curve.segments.size() - 1;
                }

                auto& sceneSegment = curve.segments[*activeSceneSegment];
                const auto vertexIndex = sceneSegment.vertices.size();
                sceneSegment.vertices.push_back(PlotSceneVertex{sample.x, sample.y, std::nullopt});
                if (sampleIndex + 1 == sourceSegment.samples.size()) {
                    sceneSegment.markUpperEndpoint = sourceSegment.markUpperEndpoint;
                    sceneSegment.upperEndpointPoint = sourceSegment.upperEndpointPoint;
                }
                segmentLocations[sampleIndex] = SceneVertexLocation{
                    true, *activeSceneSegment, vertexIndex};
            }
        }
        scene.curves.push_back(std::move(curve));
    }

    scene.anchors.reserve(anchors->anchors.size());
    for (const auto& sourceAnchor : anchors->anchors) {
        PlotSceneAnchor anchor;
        anchor.id = sourceAnchor.id;
        anchor.kinds = sourceAnchor.kinds;
        anchor.x = sourceAnchor.x;
        anchor.y = sourceAnchor.y;
        anchor.curveIds.reserve(sourceAnchor.curveIndices.size());
        for (const auto curveIndex : sourceAnchor.curveIndices) {
            if (curveIndex >= curves.size())
                return result;
            anchor.curveIds.push_back(curveId(curveIndex));
        }

        anchor.vertices.reserve(sourceAnchor.vertices.size());
        for (const auto& sourceRef : sourceAnchor.vertices) {
            if (sourceRef.curveIndex >= curves.size()
                || sourceRef.segmentIndex >= curves[sourceRef.curveIndex].segments.size()
                || sourceRef.sampleIndex >= curves[sourceRef.curveIndex]
                    .segments[sourceRef.segmentIndex].samples.size())
                return result;

            const auto& sourceSample = curves[sourceRef.curveIndex]
                .segments[sourceRef.segmentIndex].samples[sourceRef.sampleIndex];
            if (!sourceSample.finite())
                return result;
            if (sourceSample.x != sourceAnchor.x || sourceSample.y != sourceAnchor.y) {
                result.status = PlotSceneBuildStatus::AnchorCoordinateMismatch;
                return result;
            }

            const auto location = locations[sourceRef.curveIndex]
                [sourceRef.segmentIndex][sourceRef.sampleIndex];
            if (!location.valid)
                return result;
            auto& vertex = scene.curves[sourceRef.curveIndex]
                .segments[location.segmentIndex].vertices[location.vertexIndex];
            if (vertex.anchorId && *vertex.anchorId != sourceAnchor.id)
                return result;
            vertex.anchorId = sourceAnchor.id;
            anchor.vertices.push_back(PlotSceneAnchorVertexRef{
                curveId(sourceRef.curveIndex), location.segmentIndex, location.vertexIndex});
        }
        scene.anchors.push_back(std::move(anchor));
    }

    result.status = PlotSceneBuildStatus::Success;
    result.scene = std::move(scene);
    return result;
}

} // namespace mmcal::plot
