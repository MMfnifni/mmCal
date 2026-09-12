#pragma once

#include "plot_view_transform.hpp"

#include "numeric/big_float.hpp"
#include "numeric/rational.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace mmcal::plot {

enum class PlotAxisOrientation {
    Horizontal,
    Vertical
};

enum class PlotAxisPlacement {
    CrossZero,
    MinimumEdge,
    MaximumEdge
};

struct PlotTickLayout final {
    numeric::BigFloat value;
    double positionMm = 0.0;
    // 1-2-5 tick生成時のexact値。表示でBigFloat丸め誤差を露出させない。
    std::optional<numeric::Rational> exactValue;
};

struct PlotAxisLayout final {
    PlotAxisOrientation orientation = PlotAxisOrientation::Horizontal;
    PlotAxisPlacement placement = PlotAxisPlacement::CrossZero;
    double axisPositionMm = 0.0;
    std::vector<PlotTickLayout> majorTicks;
};

struct PlotAxesLayout final {
    PlotAxisLayout xAxis;
    PlotAxisLayout yAxis;
};

struct PlotAxisLayoutOptions final {
    // Ticks->Falseでは軸位置だけを計算し，tick列挙自体を行わない。
    bool generateMajorTicks = true;
    double targetMajorTickSpacingMm = 20.0;
    std::size_t minimumMajorTicks = 2;
    std::size_t maximumMajorTicks = 12;
    std::size_t maxDecimalExponentMagnitude = 1024;
    std::size_t maxTickIndexBits = 4096;
};

enum class PlotAxisLayoutStatus {
    Success,
    InvalidOptions,
    ResourceLimit,
    MappingFailed
};

struct PlotAxisLayoutResult final {
    PlotAxisLayoutStatus status = PlotAxisLayoutStatus::InvalidOptions;
    std::optional<PlotAxesLayout> layout;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PlotAxisLayoutStatus::Success && layout.has_value();
    }
};

// 1-2-5×10^nのmajor tickをmm間隔から決める。文字labelの衝突回避はTextMetrics層へ残す。
[[nodiscard]] PlotAxisLayoutResult layoutPlotAxes(
    const PlotViewTransform& transform,
    const PlotAxisLayoutOptions& options = {});

} // namespace mmcal::plot
