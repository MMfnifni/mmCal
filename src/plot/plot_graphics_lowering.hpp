#pragma once

#include "plot_axis_layout.hpp"
#include "plot_scene.hpp"
#include "plot_view_transform.hpp"

#include "graphics/graphics_scene.hpp"

#include <cstddef>
#include <optional>

namespace mmcal::plot {

struct PlotGraphicsLoweringOptions final {
    double curveStrokeWidthMm = 0.8;
    double axisStrokeWidthMm = 0.5;
    double tickStrokeWidthMm = 0.3;
    double majorTickLengthMm = 1.8;
    double endpointMarkerDiameterMm = 1.8;
    double endpointMarkerStrokeWidthMm = 0.3;
    bool showMajorTicks = true;
    bool showMajorTickLabels = true;
    bool showEndpointMarkers = true;
    double tickLabelFontSizeMm = 3.2;
    double tickLabelGapMm = 1.0;
    double minimumOuterMarginMm = 2.0;
    std::size_t maxTickLabelCharacters = 16;
};

enum class PlotGraphicsLoweringStatus {
    Success,
    InvalidOptions,
    MappingFailed
};

struct PlotGraphicsLoweringResult final {
    PlotGraphicsLoweringStatus status = PlotGraphicsLoweringStatus::InvalidOptions;
    std::optional<graphics::GraphicsScene> scene;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PlotGraphicsLoweringStatus::Success && scene.has_value();
    }
};

// PlotSceneと軸layoutをmm単位の汎用GraphicsSceneへloweringする。
[[nodiscard]] PlotGraphicsLoweringResult lowerPlotToGraphics(
    const PlotScene& plotScene,
    const PlotAxesLayout& axes,
    const PlotViewTransform& transform,
    const PlotGraphicsLoweringOptions& options = {});

} // namespace mmcal::plot
