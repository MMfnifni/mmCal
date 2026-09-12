#pragma once

#include "plot_anchors.hpp"
#include "plot_sampling.hpp"

#include "numeric/big_float.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace mmcal::plot {

using PlotSceneCurveId = std::uint64_t;

struct PlotSceneVertex final {
    numeric::BigFloat x;
    numeric::BigFloat y;
    std::optional<PlotAnchorId> anchorId;
};

struct PlotSceneSegment final {
    PlotInterval sourceInterval;
    PlotSegmentGeometryKind geometryKind = PlotSegmentGeometryKind::Polyline;
    std::vector<PlotSceneVertex> vertices;
    std::optional<PlotEndpointPoint> quadraticControlPoint;
    std::optional<PlotBezierControlPoints> bezierControlPoints;
    std::optional<PlotEllipseGeometry> ellipseGeometry;
    std::optional<PlotEllipticArcGeometry> ellipticArcGeometry;
    bool markLowerEndpoint = false;
    bool markUpperEndpoint = false;
    std::optional<PlotEndpointPoint> lowerEndpointPoint;
    std::optional<PlotEndpointPoint> upperEndpointPoint;
};

struct PlotSceneCurve final {
    PlotSceneCurveId id = 0;
    std::vector<PlotSceneSegment> segments;
};

struct PlotSceneAnchorVertexRef final {
    PlotSceneCurveId curveId = 0;
    std::size_t segmentIndex = 0;
    std::size_t vertexIndex = 0;
};

// backend非依存の共有semantic point。sampling上のindexではなくScene内vertexを参照する。
struct PlotSceneAnchor final {
    PlotAnchorId id = 0;
    PlotAnchorKind kinds = PlotAnchorKind::None;
    numeric::BigFloat x;
    numeric::BigFloat y;
    std::vector<PlotSceneCurveId> curveIds;
    std::vector<PlotSceneAnchorVertexRef> vertices;
};

// PlotSceneは数学曲線と意味点だけを保持する。軸・tick・style等は後段で拡張する。
struct PlotScene final {
    std::vector<PlotSceneCurve> curves;
    std::vector<PlotSceneAnchor> anchors;
};

enum class PlotSceneBuildStatus {
    Success,
    InvalidAnchorReference,
    AnchorCoordinateMismatch
};

struct PlotSceneBuildResult final {
    PlotSceneBuildStatus status = PlotSceneBuildStatus::InvalidAnchorReference;
    std::optional<PlotScene> scene;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PlotSceneBuildStatus::Success && scene.has_value();
    }
};

// SampledCurveから純sampling tagを除去してPlotSceneへloweringする。
// anchorsを省略した場合はsemantic tagから自動収集する。
[[nodiscard]] PlotSceneBuildResult buildPlotScene(
    const std::vector<SampledCurve>& curves,
    std::optional<PlotAnchorSet> anchors = std::nullopt);

} // namespace mmcal::plot
