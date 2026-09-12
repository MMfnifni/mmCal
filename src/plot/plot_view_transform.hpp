#pragma once

#include "plot_program.hpp"
#include "plot_range_estimator.hpp"
#include "plot_request.hpp"

#include "evaluation/builtin_registry.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/big_float.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace mmcal::plot {

struct PlotViewportMm final {
    double widthMm = 150.0;
    double heightMm = 100.0;
};

struct PlotViewPointMm final {
    double xMm = 0.0;
    double yMm = 0.0;
};

// data座標からPlotArea内のmm座標へ写すtransform。
// Plot/Graphics内部では物理長さをmmで統一し，pixel変換はRaster/GUI backendへ残す。
// yは上向きを正とし，backend固有の上下反転もGraphics層の責務とする。
class PlotViewTransform final {
public:
    PlotViewTransform(
        numeric::BigFloat xMinimum,
        numeric::BigFloat xMaximum,
        numeric::BigFloat yMinimum,
        numeric::BigFloat yMaximum,
        std::size_t precisionBits,
        PlotViewportMm viewport);

    [[nodiscard]] std::optional<PlotViewPointMm> map(
        const numeric::BigFloat& x,
        const numeric::BigFloat& y) const;
    [[nodiscard]] std::optional<double> mapX(const numeric::BigFloat& x) const;
    [[nodiscard]] std::optional<double> mapY(const numeric::BigFloat& y) const;

    [[nodiscard]] const numeric::BigFloat& xMinimum() const noexcept { return xMinimum_; }
    [[nodiscard]] const numeric::BigFloat& xMaximum() const noexcept { return xMaximum_; }
    [[nodiscard]] const numeric::BigFloat& yMinimum() const noexcept { return yMinimum_; }
    [[nodiscard]] const numeric::BigFloat& yMaximum() const noexcept { return yMaximum_; }
    [[nodiscard]] std::size_t precisionBits() const noexcept { return precisionBits_; }
    [[nodiscard]] const PlotViewportMm& viewport() const noexcept { return viewport_; }

private:
    numeric::BigFloat xMinimum_;
    numeric::BigFloat xMaximum_;
    numeric::BigFloat yMinimum_;
    numeric::BigFloat yMaximum_;
    std::size_t precisionBits_ = 0;
    PlotViewportMm viewport_;
};


// 既に確定した2D view rangeからbackend-independentなtransformを作る。
// ParametricPlotではx/y双方をAutomatic推定してこの経路へ入れる。
[[nodiscard]] std::optional<PlotViewTransform> makePlotViewTransform(
    const CurveViewRange2D& range,
    std::size_t precisionBits,
    PlotViewportMm viewport = {});

[[nodiscard]] std::optional<PlotViewTransform> makePlotViewTransform(
    const std::vector<PlotRequest>& requests,
    const PlotRangeEstimate& range,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits,
    PlotViewportMm viewport = {});

[[nodiscard]] std::optional<PlotViewTransform> makePlotViewTransform(
    const PlotRequest& request,
    const PlotRangeEstimate& range,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits,
    PlotViewportMm viewport = {});

} // namespace mmcal::plot
