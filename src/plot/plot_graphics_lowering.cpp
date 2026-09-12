#include "plot_graphics_lowering.hpp"

#include "graphics/text_metrics.hpp"
#include "numeric/integer_algorithms.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <limits>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::plot {
namespace {

using graphics::GraphicsCircleNode;
using graphics::GraphicsEllipseNode;
using graphics::GraphicsEllipticArcNode;
using graphics::GraphicsCubicTo;
using graphics::GraphicsFillStyle;
using graphics::GraphicsLineTo;
using graphics::GraphicsMoveTo;
using graphics::GraphicsPathNode;
using graphics::GraphicsPointMm;
using graphics::GraphicsQuadraticTo;
using graphics::GraphicsSemanticKind;
using graphics::GraphicsSemanticRef;
using graphics::GraphicsStrokeStyle;
using graphics::GraphicsTextAnchor;
using graphics::GraphicsTextMetrics;
using graphics::GraphicsTextNode;
using numeric::BigInt;
using numeric::Rational;


[[nodiscard]] double rationalToDouble(const Rational& value) {
    try {
        const long double numerator = std::stold(value.numerator().toString());
        const long double denominator = std::stold(value.denominator().toString());
        return static_cast<double>(numerator / denominator);
    }
    catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

struct TickLabel final {
    std::string text;
    GraphicsTextMetrics metrics;
};

struct PlotMarginsMm final {
    double left = 0.0;
    double right = 0.0;
    double bottom = 0.0;
    double top = 0.0;
};

[[nodiscard]] GraphicsStrokeStyle stroke(double widthMm) {
    GraphicsStrokeStyle result;
    result.widthMm = widthMm;
    return result;
}

[[nodiscard]] GraphicsPathNode line(
    GraphicsPointMm from,
    GraphicsPointMm to,
    double widthMm,
    GraphicsSemanticRef semantic) {
    GraphicsPathNode path;
    path.commands.push_back(GraphicsMoveTo{from});
    path.commands.push_back(GraphicsLineTo{to});
    path.stroke = stroke(widthMm);
    path.semantic = semantic;
    return path;
}



struct PendingEndpointMarker final {
    GraphicsPointMm center;
    bool filled = false;
    GraphicsSemanticRef semantic;
};

[[nodiscard]] GraphicsFillStyle fill(const graphics::GraphicsColor& color) {
    return GraphicsFillStyle{color};
}

[[nodiscard]] graphics::GraphicsColor black() noexcept {
    return graphics::GraphicsColor{};
}

[[nodiscard]] graphics::GraphicsColor white() noexcept {
    return graphics::GraphicsColor{255, 255, 255, 255};
}

[[nodiscard]] bool samePoint(GraphicsPointMm lhs, GraphicsPointMm rhs) noexcept {
    constexpr double tolerance = 1e-9;
    return std::abs(lhs.xMm - rhs.xMm) <= tolerance && std::abs(lhs.yMm - rhs.yMm) <= tolerance;
}

[[nodiscard]] GraphicsPointMm interpolatePoint(
    GraphicsPointMm lhs,
    GraphicsPointMm rhs,
    double t) noexcept {
    return {
        lhs.xMm + (rhs.xMm - lhs.xMm) * t,
        lhs.yMm + (rhs.yMm - lhs.yMm) * t};
}


[[nodiscard]] bool clipLineSegment(
    GraphicsPointMm from,
    GraphicsPointMm to,
    const graphics::GraphicsRectMm& rect,
    GraphicsPointMm& clippedFrom,
    GraphicsPointMm& clippedTo) noexcept {
    const double minX = rect.xMm;
    const double maxX = rect.xMm + rect.widthMm;
    const double minY = rect.yMm;
    const double maxY = rect.yMm + rect.heightMm;
    const double dx = to.xMm - from.xMm;
    const double dy = to.yMm - from.yMm;
    double lower = 0.0;
    double upper = 1.0;

    const auto clipBoundary = [&](double p, double q) {
        constexpr double epsilon = 1e-15;
        if (std::abs(p) <= epsilon)
            return q >= 0.0;
        const double ratio = q / p;
        if (p < 0.0) {
            if (ratio > upper)
                return false;
            lower = std::max(lower, ratio);
        }
        else {
            if (ratio < lower)
                return false;
            upper = std::min(upper, ratio);
        }
        return lower <= upper;
    };

    if (!clipBoundary(-dx, from.xMm - minX)
        || !clipBoundary(dx, maxX - from.xMm)
        || !clipBoundary(-dy, from.yMm - minY)
        || !clipBoundary(dy, maxY - from.yMm))
        return false;

    clippedFrom = interpolatePoint(from, to, lower);
    clippedTo = interpolatePoint(from, to, upper);
    return true;
}

void appendMoveTo(
    GraphicsPathNode& path,
    std::optional<GraphicsPointMm>& current,
    GraphicsPointMm point) {
    if (!current || !samePoint(*current, point))
        path.commands.push_back(GraphicsMoveTo{point});
    current = point;
}

void appendClippedLine(
    GraphicsPathNode& path,
    std::optional<GraphicsPointMm>& current,
    GraphicsPointMm from,
    GraphicsPointMm to,
    const graphics::GraphicsRectMm& rect) {
    GraphicsPointMm clippedFrom;
    GraphicsPointMm clippedTo;
    if (!clipLineSegment(from, to, rect, clippedFrom, clippedTo)) {
        current.reset();
        return;
    }

    // 矩形への接触が一点だけなら描画primitiveを持たない。
    // MoveToだけを残すとbackendへ範囲外座標が露出し得るため，完全に捨てる。
    if (samePoint(clippedFrom, clippedTo)) {
        if (!current || !samePoint(*current, clippedFrom))
            current.reset();
        return;
    }

    appendMoveTo(path, current, clippedFrom);
    path.commands.push_back(GraphicsLineTo{clippedTo});
    current = clippedTo;
}

void appendClippedPolylineCommands(
    GraphicsPathNode& path,
    const std::vector<GraphicsPointMm>& points,
    const graphics::GraphicsRectMm& rect) {
    if (points.size() < 2)
        return;

    std::optional<GraphicsPointMm> current;
    for (std::size_t i = 1; i < points.size(); ++i)
        appendClippedLine(path, current, points[i - 1], points[i], rect);
}

struct QuadraticBezierMm final {
    GraphicsPointMm p0;
    GraphicsPointMm p1;
    GraphicsPointMm p2;
};

struct CubicBezierMm final {
    GraphicsPointMm p0;
    GraphicsPointMm p1;
    GraphicsPointMm p2;
    GraphicsPointMm p3;
};


struct BezierBoundsMm final {
    double minX = 0.0;
    double maxX = 0.0;
    double minY = 0.0;
    double maxY = 0.0;
};

[[nodiscard]] BezierBoundsMm bounds(const QuadraticBezierMm& curve) noexcept {
    return {
        std::min({curve.p0.xMm, curve.p1.xMm, curve.p2.xMm}),
        std::max({curve.p0.xMm, curve.p1.xMm, curve.p2.xMm}),
        std::min({curve.p0.yMm, curve.p1.yMm, curve.p2.yMm}),
        std::max({curve.p0.yMm, curve.p1.yMm, curve.p2.yMm})};
}

[[nodiscard]] BezierBoundsMm bounds(const CubicBezierMm& curve) noexcept {
    return {
        std::min({curve.p0.xMm, curve.p1.xMm, curve.p2.xMm, curve.p3.xMm}),
        std::max({curve.p0.xMm, curve.p1.xMm, curve.p2.xMm, curve.p3.xMm}),
        std::min({curve.p0.yMm, curve.p1.yMm, curve.p2.yMm, curve.p3.yMm}),
        std::max({curve.p0.yMm, curve.p1.yMm, curve.p2.yMm, curve.p3.yMm})};
}

[[nodiscard]] bool outside(
    const BezierBoundsMm& box,
    const graphics::GraphicsRectMm& rect) noexcept {
    return box.maxX < rect.xMm || box.minX > rect.xMm + rect.widthMm
        || box.maxY < rect.yMm || box.minY > rect.yMm + rect.heightMm;
}

[[nodiscard]] bool insideExpanded(
    const BezierBoundsMm& box,
    const graphics::GraphicsRectMm& rect,
    double marginMm) noexcept {
    return box.minX >= rect.xMm - marginMm
        && box.maxX <= rect.xMm + rect.widthMm + marginMm
        && box.minY >= rect.yMm - marginMm
        && box.maxY <= rect.yMm + rect.heightMm + marginMm;
}

[[nodiscard]] std::pair<QuadraticBezierMm, QuadraticBezierMm> splitBezier(
    const QuadraticBezierMm& curve,
    double t) noexcept {
    const auto p01 = interpolatePoint(curve.p0, curve.p1, t);
    const auto p12 = interpolatePoint(curve.p1, curve.p2, t);
    const auto point = interpolatePoint(p01, p12, t);
    return {
        QuadraticBezierMm{curve.p0, p01, point},
        QuadraticBezierMm{point, p12, curve.p2}};
}

[[nodiscard]] std::pair<CubicBezierMm, CubicBezierMm> splitBezier(
    const CubicBezierMm& curve,
    double t) noexcept {
    const auto p01 = interpolatePoint(curve.p0, curve.p1, t);
    const auto p12 = interpolatePoint(curve.p1, curve.p2, t);
    const auto p23 = interpolatePoint(curve.p2, curve.p3, t);
    const auto p012 = interpolatePoint(p01, p12, t);
    const auto p123 = interpolatePoint(p12, p23, t);
    const auto point = interpolatePoint(p012, p123, t);
    return {
        CubicBezierMm{curve.p0, p01, p012, point},
        CubicBezierMm{point, p123, p23, curve.p3}};
}


void appendClippedQuadraticBezier(
    GraphicsPathNode& path,
    std::optional<GraphicsPointMm>& current,
    const QuadraticBezierMm& curve,
    const graphics::GraphicsRectMm& rect,
    std::size_t depth = 0) {
    const auto box = bounds(curve);
    if (outside(box, rect)) {
        current.reset();
        return;
    }

    // clip境界近傍までde Casteljau分割し，巨大なoff-screen制御点をbackendへ渡さない。
    constexpr double boundaryMarginMm = 0.01;
    constexpr std::size_t maxDepth = 48;
    if (insideExpanded(box, rect, boundaryMarginMm)) {
        appendMoveTo(path, current, curve.p0);
        path.commands.push_back(GraphicsQuadraticTo{curve.p1, curve.p2});
        current = curve.p2;
        return;
    }
    if (depth >= maxDepth) {
        appendClippedLine(path, current, curve.p0, curve.p2, rect);
        return;
    }

    const auto [left, right] = splitBezier(curve, 0.5);
    appendClippedQuadraticBezier(path, current, left, rect, depth + 1);
    appendClippedQuadraticBezier(path, current, right, rect, depth + 1);
}

void appendClippedCubicBezier(
    GraphicsPathNode& path,
    std::optional<GraphicsPointMm>& current,
    const CubicBezierMm& curve,
    const graphics::GraphicsRectMm& rect,
    std::size_t depth = 0) {
    const auto box = bounds(curve);
    if (outside(box, rect)) {
        current.reset();
        return;
    }

    constexpr double boundaryMarginMm = 0.01;
    constexpr std::size_t maxDepth = 48;
    if (insideExpanded(box, rect, boundaryMarginMm)) {
        appendMoveTo(path, current, curve.p0);
        path.commands.push_back(GraphicsCubicTo{curve.p1, curve.p2, curve.p3});
        current = curve.p3;
        return;
    }
    if (depth >= maxDepth) {
        appendClippedLine(path, current, curve.p0, curve.p3, rect);
        return;
    }

    const auto [left, right] = splitBezier(curve, 0.5);
    appendClippedCubicBezier(path, current, left, rect, depth + 1);
    appendClippedCubicBezier(path, current, right, rect, depth + 1);
}

[[nodiscard]] std::vector<double> anchorParameters(
    const PlotSceneSegment& segment,
    const std::vector<GraphicsPointMm>& mappedVertices) {
    std::vector<double> result;
    if (segment.vertices.size() != mappedVertices.size() || mappedVertices.size() < 2)
        return result;
    const double x0 = mappedVertices.front().xMm;
    const double x1 = mappedVertices.back().xMm;
    const double width = x1 - x0;
    if (!(width > 0.0))
        return result;

    for (std::size_t i = 1; i + 1 < segment.vertices.size(); ++i) {
        if (!segment.vertices[i].anchorId)
            continue;
        const double t = (mappedVertices[i].xMm - x0) / width;
        if (t > 0.0 && t < 1.0)
            result.push_back(t);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end(), [](double lhs, double rhs) {
        return std::abs(lhs - rhs) <= 1e-12;
    }), result.end());
    return result;
}

void appendQuadraticBezierCommands(
    GraphicsPathNode& path,
    QuadraticBezierMm curve,
    const std::vector<double>& parameters,
    const graphics::GraphicsRectMm& clipRect) {
    double previous = 0.0;
    std::optional<GraphicsPointMm> current;
    for (const double globalT : parameters) {
        const double denominator = 1.0 - previous;
        if (!(denominator > 0.0))
            break;
        const double localT = std::clamp((globalT - previous) / denominator, 0.0, 1.0);
        const auto [left, right] = splitBezier(curve, localT);
        appendClippedQuadraticBezier(path, current, left, clipRect);
        curve = right;
        previous = globalT;
    }
    appendClippedQuadraticBezier(path, current, curve, clipRect);
}

void appendCubicBezierCommands(
    GraphicsPathNode& path,
    CubicBezierMm curve,
    const std::vector<double>& parameters,
    const graphics::GraphicsRectMm& clipRect) {
    double previous = 0.0;
    std::optional<GraphicsPointMm> current;
    for (const double globalT : parameters) {
        const double denominator = 1.0 - previous;
        if (!(denominator > 0.0))
            break;
        const double localT = std::clamp((globalT - previous) / denominator, 0.0, 1.0);
        const auto [left, right] = splitBezier(curve, localT);
        appendClippedCubicBezier(path, current, left, clipRect);
        curve = right;
        previous = globalT;
    }
    appendClippedCubicBezier(path, current, curve, clipRect);
}

void appendMarker(std::vector<PendingEndpointMarker>& markers, PendingEndpointMarker marker) {
    for (const auto& existing : markers) {
        if (existing.filled == marker.filled && samePoint(existing.center, marker.center))
            return;
    }
    markers.push_back(marker);
}

[[nodiscard]] std::optional<GraphicsPointMm> reconstructStraightLineEndpoint(
    const std::vector<GraphicsPointMm>& mappedVertices,
    PlotEndpointInclusion lowerInclusion,
    PlotEndpointInclusion upperInclusion,
    bool lowerEndpoint) {
    if (mappedVertices.size() < 2)
        return std::nullopt;

    const auto& first = mappedVertices.front();
    const auto& last = mappedVertices.back();
    if (lowerEndpoint && lowerInclusion == PlotEndpointInclusion::Closed)
        return first;
    if (!lowerEndpoint && upperInclusion == PlotEndpointInclusion::Closed)
        return last;

    const double denominator = lowerInclusion == PlotEndpointInclusion::Open
            && upperInclusion == PlotEndpointInclusion::Open
        ? 1022.0
        : 1023.0;
    const double dx = last.xMm - first.xMm;
    const double dy = last.yMm - first.yMm;
    if (lowerEndpoint) {
        return GraphicsPointMm{first.xMm - dx / denominator, first.yMm - dy / denominator};
    }
    return GraphicsPointMm{last.xMm + dx / denominator, last.yMm + dy / denominator};
}

[[nodiscard]] std::optional<GraphicsPointMm> reconstructPiecewiseConstantEndpoint(
    const std::vector<GraphicsPointMm>& mappedVertices,
    bool lowerEndpoint) {
    if (mappedVertices.size() < 2)
        return std::nullopt;

    // samplingは両端を1/1024だけ内側へ寄せているため，2点間は1022/1024幅。
    // yはpiece内部の定数値をそのまま使い，xだけ数学的境界へ戻す。
    constexpr double denominator = 1022.0;
    const auto& first = mappedVertices.front();
    const auto& last = mappedVertices.back();
    const double dx = last.xMm - first.xMm;
    return lowerEndpoint
        ? GraphicsPointMm{first.xMm - dx / denominator, first.yMm}
        : GraphicsPointMm{last.xMm + dx / denominator, last.yMm};
}

void appendSegmentEndpointMarkers(
    std::vector<PendingEndpointMarker>& markers,
    PlotSceneCurveId curveId,
    const PlotSceneSegment& segment,
    const std::vector<GraphicsPointMm>& mappedVertices,
    std::optional<GraphicsPointMm> mappedLowerEndpointPoint,
    std::optional<GraphicsPointMm> mappedUpperEndpointPoint) {
    if (mappedVertices.empty())
        return;
    if (mappedVertices.size() == 1) {
        if (segment.markLowerEndpoint || segment.markUpperEndpoint) {
            appendMarker(markers, PendingEndpointMarker{
                mappedVertices.front(), true,
                GraphicsSemanticRef{GraphicsSemanticKind::Point, curveId}});
        }
        return;
    }

    const auto endpointCenter = [&](bool lowerEndpoint) -> std::optional<GraphicsPointMm> {
        const auto explicitPoint = lowerEndpoint
            ? mappedLowerEndpointPoint : mappedUpperEndpointPoint;
        if (explicitPoint)
            return explicitPoint;
        if (segment.geometryKind == PlotSegmentGeometryKind::PiecewiseConstant)
            return reconstructPiecewiseConstantEndpoint(mappedVertices, lowerEndpoint);

        const auto inclusion = lowerEndpoint
            ? segment.sourceInterval.lowerInclusion
            : segment.sourceInterval.upperInclusion;
        if (inclusion == PlotEndpointInclusion::Closed)
            return lowerEndpoint ? mappedVertices.front() : mappedVertices.back();
        if (segment.geometryKind == PlotSegmentGeometryKind::StraightLine)
            return reconstructStraightLineEndpoint(
                mappedVertices,
                segment.sourceInterval.lowerInclusion,
                segment.sourceInterval.upperInclusion,
                lowerEndpoint);
        return std::nullopt;
    };

    if (segment.markLowerEndpoint) {
        if (const auto center = endpointCenter(true)) {
            appendMarker(markers, PendingEndpointMarker{
                *center,
                segment.sourceInterval.lowerInclusion == PlotEndpointInclusion::Closed,
                GraphicsSemanticRef{GraphicsSemanticKind::CurveEndpoint, curveId}});
        }
    }

    if (segment.markUpperEndpoint) {
        if (const auto center = endpointCenter(false)) {
            appendMarker(markers, PendingEndpointMarker{
                *center,
                segment.sourceInterval.upperInclusion == PlotEndpointInclusion::Closed,
                GraphicsSemanticRef{GraphicsSemanticKind::CurveEndpoint, curveId}});
        }
    }
}

[[nodiscard]] bool markerIntersectsClip(
    const PendingEndpointMarker& marker,
    const PlotGraphicsLoweringOptions& options,
    const graphics::GraphicsRectMm& clipRect) noexcept {
    const double radius = options.endpointMarkerDiameterMm / 2.0;
    return marker.center.xMm + radius >= clipRect.xMm
        && marker.center.xMm - radius <= clipRect.xMm + clipRect.widthMm
        && marker.center.yMm + radius >= clipRect.yMm
        && marker.center.yMm - radius <= clipRect.yMm + clipRect.heightMm;
}

[[nodiscard]] GraphicsCircleNode endpointCircle(
    const PendingEndpointMarker& marker,
    const PlotGraphicsLoweringOptions& options,
    graphics::GraphicsRectMm clipRect) {
    GraphicsCircleNode circle;
    circle.center = marker.center;
    circle.radiusMm = options.endpointMarkerDiameterMm / 2.0;
    circle.stroke = stroke(options.endpointMarkerStrokeWidthMm);
    circle.fill = fill(marker.filled ? black() : white());
    circle.semantic = marker.semantic;
    circle.clipRect = clipRect;
    return circle;
}

[[nodiscard]] std::pair<double, double> tickOffsets(
    PlotAxisPlacement placement,
    double lengthMm) noexcept {
    if (placement == PlotAxisPlacement::MinimumEdge)
        return {0.0, lengthMm};
    if (placement == PlotAxisPlacement::MaximumEdge)
        return {-lengthMm, 0.0};
    return {-lengthMm / 2.0, lengthMm / 2.0};
}

[[nodiscard]] bool divideFactor(BigInt& value, std::int64_t factor) {
    const BigInt divisor{factor};
    const auto divided = numeric::divmod(value, divisor);
    if (!divided.remainder.isZero())
        return false;
    value = divided.quotient;
    return true;
}

struct DecimalValue final {
    bool negative = false;
    std::string digits;
    std::size_t scale = 0;
};

// nice tickは有限小数になるため，exact Rationalから10進表記を組み立てる。
// BigFloatを文字列化すると0.2等でbinary丸め桁が露出するので表示には使わない。
[[nodiscard]] std::optional<DecimalValue> finiteDecimal(const Rational& value) {
    BigInt denominator = value.denominator();
    std::size_t twos = 0;
    std::size_t fives = 0;
    while (divideFactor(denominator, 2))
        ++twos;
    while (divideFactor(denominator, 5))
        ++fives;
    if (denominator != BigInt{1})
        return std::nullopt;

    const std::size_t scale = std::max(twos, fives);
    BigInt scaled = value.numerator().abs();
    if (scale > twos)
        scaled *= numeric::pow(BigInt{2}, static_cast<std::uint64_t>(scale - twos));
    if (scale > fives)
        scaled *= numeric::pow(BigInt{5}, static_cast<std::uint64_t>(scale - fives));

    return DecimalValue{value.numerator().isNegative(), scaled.toString(), scale};
}

[[nodiscard]] std::string fixedDecimal(const DecimalValue& value) {
    if (value.digits == "0")
        return "0";

    std::string digits = value.digits;
    if (digits.size() <= value.scale)
        digits.insert(0, value.scale + 1 - digits.size(), '0');

    const std::size_t point = digits.size() - value.scale;
    std::string result = digits.substr(0, point);
    if (value.scale != 0) {
        std::string fractional = digits.substr(point);
        while (!fractional.empty() && fractional.back() == '0')
            fractional.pop_back();
        if (!fractional.empty()) {
            result.push_back('.');
            result += fractional;
        }
    }
    if (value.negative)
        result.insert(result.begin(), '-');
    return result;
}

[[nodiscard]] std::string scientificDecimal(
    const DecimalValue& value,
    std::size_t significantDigits) {
    if (value.digits == "0")
        return "0";

    const auto exponent = static_cast<std::int64_t>(value.digits.size())
        - static_cast<std::int64_t>(value.scale) - 1;
    const std::size_t count = std::min(significantDigits, value.digits.size());
    std::string mantissa = value.digits.substr(0, count);
    while (mantissa.size() > 1 && mantissa.back() == '0')
        mantissa.pop_back();
    if (mantissa.size() > 1)
        mantissa.insert(mantissa.begin() + 1, '.');
    if (value.negative)
        mantissa.insert(mantissa.begin(), '-');
    return mantissa + "e" + (exponent >= 0 ? "+" : "") + std::to_string(exponent);
}

[[nodiscard]] std::string formatTickValue(
    const PlotTickLayout& tick,
    std::size_t maxCharacters) {
    if (!tick.exactValue)
        return tick.value.toRational().toString();
    const auto decimal = finiteDecimal(*tick.exactValue);
    if (!decimal)
        return tick.exactValue->toString();

    const std::string fixed = fixedDecimal(*decimal);
    if (fixed.size() <= maxCharacters)
        return fixed;
    return scientificDecimal(*decimal, std::max<std::size_t>(1, maxCharacters / 2));
}

[[nodiscard]] std::optional<TickLabel> makeTickLabel(
    const PlotTickLayout& tick,
    const PlotGraphicsLoweringOptions& options) {
    std::string text = formatTickValue(tick, options.maxTickLabelCharacters);
    const auto metrics = graphics::measureTextApproximate(text, options.tickLabelFontSizeMm);
    if (!metrics)
        return std::nullopt;
    return TickLabel{std::move(text), *metrics};
}

[[nodiscard]] PlotMarginsMm computeMargins(
    const PlotAxesLayout& axes,
    const std::vector<TickLabel>& xLabels,
    const std::vector<TickLabel>& yLabels,
    const PlotGraphicsLoweringOptions& options) {
    PlotMarginsMm margins{
        options.minimumOuterMarginMm,
        options.minimumOuterMarginMm,
        options.minimumOuterMarginMm,
        options.minimumOuterMarginMm};
    if (!options.showMajorTickLabels)
        return margins;

    double maxXWidth = 0.0;
    double maxXHeight = 0.0;
    for (const auto& label : xLabels) {
        maxXWidth = std::max(maxXWidth, label.metrics.advanceWidthMm);
        maxXHeight = std::max(maxXHeight, label.metrics.heightMm());
    }
    double maxYWidth = 0.0;
    double maxYHeight = 0.0;
    for (const auto& label : yLabels) {
        maxYWidth = std::max(maxYWidth, label.metrics.advanceWidthMm);
        maxYHeight = std::max(maxYHeight, label.metrics.heightMm());
    }

    // 端のlabelがplot areaから半分はみ出す分もcanvas側へ確保する。
    margins.left = std::max(margins.left, maxXWidth / 2.0 + options.tickLabelGapMm);
    margins.right = std::max(margins.right, maxXWidth / 2.0 + options.tickLabelGapMm);
    margins.bottom = std::max(margins.bottom, maxYHeight / 2.0 + options.tickLabelGapMm);
    margins.top = std::max(margins.top, maxYHeight / 2.0 + options.tickLabelGapMm);

    if (axes.xAxis.placement == PlotAxisPlacement::MinimumEdge)
        margins.bottom = std::max(
            margins.bottom, options.majorTickLengthMm + options.tickLabelGapMm + maxXHeight);
    else if (axes.xAxis.placement == PlotAxisPlacement::MaximumEdge)
        margins.top = std::max(
            margins.top, options.majorTickLengthMm + options.tickLabelGapMm + maxXHeight);

    if (axes.yAxis.placement == PlotAxisPlacement::MinimumEdge)
        margins.left = std::max(
            margins.left, options.majorTickLengthMm + options.tickLabelGapMm + maxYWidth);
    else if (axes.yAxis.placement == PlotAxisPlacement::MaximumEdge)
        margins.right = std::max(
            margins.right, options.majorTickLengthMm + options.tickLabelGapMm + maxYWidth);
    return margins;
}

struct PlotAreaMm final {
    double left = 0.0;
    double bottom = 0.0;
    double width = 0.0;
    double height = 0.0;
};

[[nodiscard]] std::optional<PlotAreaMm> plotArea(
    const PlotViewportMm& canvas,
    const PlotMarginsMm& margins) noexcept {
    const double width = canvas.widthMm - margins.left - margins.right;
    const double height = canvas.heightMm - margins.bottom - margins.top;
    if (!(width > 0.0) || !(height > 0.0))
        return std::nullopt;
    return PlotAreaMm{margins.left, margins.bottom, width, height};
}


[[nodiscard]] std::optional<double> mapCoordinateForGraphics(
    const numeric::BigFloat& value,
    const numeric::BigFloat& minimum,
    const numeric::BigFloat& maximum,
    std::optional<double> mapped,
    double viewportSizeMm) {
    if (mapped)
        return mapped;

    // BigFloat値がdouble換算範囲を越えても，viewport外であることをexact比較できれば
    // lowering用の有限sentinelへ置く。直後の幾何clipで捨てるため数学rangeは変更しない。
    constexpr double overflowGuard = 1.0e6;
    if (value < minimum)
        return -viewportSizeMm * overflowGuard;
    if (value > maximum)
        return viewportSizeMm * (1.0 + overflowGuard);
    return std::nullopt;
}

[[nodiscard]] std::optional<PlotViewPointMm> mapForGraphics(
    const PlotViewTransform& transform,
    const numeric::BigFloat& x,
    const numeric::BigFloat& y) {
    const auto xMm = mapCoordinateForGraphics(
        x, transform.xMinimum(), transform.xMaximum(), transform.mapX(x),
        transform.viewport().widthMm);
    const auto yMm = mapCoordinateForGraphics(
        y, transform.yMinimum(), transform.yMaximum(), transform.mapY(y),
        transform.viewport().heightMm);
    if (!xMm || !yMm)
        return std::nullopt;
    return PlotViewPointMm{*xMm, *yMm};
}

[[nodiscard]] GraphicsPointMm fitToPlotArea(
    GraphicsPointMm point,
    const PlotViewportMm& canvas,
    const PlotAreaMm& area) noexcept {
    return {
        area.left + point.xMm * area.width / canvas.widthMm,
        area.bottom + point.yMm * area.height / canvas.heightMm};
}

[[nodiscard]] bool skipOriginTick(
    const PlotAxesLayout& axes,
    const PlotTickLayout& tick,
    bool yAxis) noexcept {
    if (!tick.value.isZero())
        return false;
    return yAxis
        ? axes.xAxis.placement == PlotAxisPlacement::CrossZero
        : axes.yAxis.placement == PlotAxisPlacement::CrossZero;
}

} // namespace

PlotGraphicsLoweringResult lowerPlotToGraphics(
    const PlotScene& plotScene,
    const PlotAxesLayout& axes,
    const PlotViewTransform& transform,
    const PlotGraphicsLoweringOptions& options) {
    PlotGraphicsLoweringResult result;
    if (!(options.curveStrokeWidthMm > 0.0)
        || !(options.axisStrokeWidthMm > 0.0)
        || !(options.tickStrokeWidthMm > 0.0)
        || !(options.majorTickLengthMm > 0.0)
        || !(options.endpointMarkerDiameterMm > 0.0)
        || !(options.endpointMarkerStrokeWidthMm > 0.0)
        || !(options.tickLabelFontSizeMm > 0.0)
        || !(options.tickLabelGapMm >= 0.0)
        || !(options.minimumOuterMarginMm >= 0.0)
        || options.maxTickLabelCharacters < 4)
        return result;

    std::vector<TickLabel> xLabels;
    std::vector<TickLabel> yLabels;
    if (options.showMajorTickLabels) {
        xLabels.reserve(axes.xAxis.majorTicks.size());
        yLabels.reserve(axes.yAxis.majorTicks.size());
        for (const auto& tick : axes.xAxis.majorTicks) {
            const auto label = makeTickLabel(tick, options);
            if (!label)
                return result;
            xLabels.push_back(*label);
        }
        for (const auto& tick : axes.yAxis.majorTicks) {
            const auto label = makeTickLabel(tick, options);
            if (!label)
                return result;
            yLabels.push_back(*label);
        }
    }

    const PlotMarginsMm margins = computeMargins(axes, xLabels, yLabels, options);
    const auto area = plotArea(transform.viewport(), margins);
    if (!area)
        return result;

    graphics::GraphicsScene scene;
    std::vector<PendingEndpointMarker> endpointMarkers;
    // viewportは最終canvas寸法である。label marginはcanvas外へ足さず，内部plot areaを縮める。
    scene.extent = graphics::GraphicsExtentMm{
        transform.viewport().widthMm, transform.viewport().heightMm};

    for (const auto& curve : plotScene.curves) {
        for (const auto& segment : curve.segments) {
            if (segment.vertices.empty())
                continue;
            std::vector<GraphicsPointMm> mappedVertices;
            mappedVertices.reserve(segment.vertices.size());
            GraphicsPathNode path;
            path.stroke = stroke(options.curveStrokeWidthMm);
            path.semantic = GraphicsSemanticRef{GraphicsSemanticKind::Curve, curve.id};
            path.clipRect = graphics::GraphicsRectMm{
                area->left, area->bottom, area->width, area->height};
            for (const auto& vertex : segment.vertices) {
                const auto point = mapForGraphics(transform, vertex.x, vertex.y);
                if (!point) {
                    result.status = PlotGraphicsLoweringStatus::MappingFailed;
                    return result;
                }
                mappedVertices.push_back(fitToPlotArea(
                    {point->xMm, point->yMm}, transform.viewport(), *area));
            }

            if (segment.geometryKind == PlotSegmentGeometryKind::Ellipse
                && segment.ellipseGeometry) {
                const auto mapDataPoint = [&](const PlotEndpointPoint& point)
                    -> std::optional<GraphicsPointMm> {
                    const auto mapped = mapForGraphics(transform, point.x, point.y);
                    if (!mapped)
                        return std::nullopt;
                    return fitToPlotArea(
                        {mapped->xMm, mapped->yMm}, transform.viewport(), *area);
                };
                const auto center = mapDataPoint(segment.ellipseGeometry->center);
                const auto offsetEndpoint = [&](const PlotEndpointPoint& vector) {
                    return PlotEndpointPoint{
                        numeric::add(
                            segment.ellipseGeometry->center.x, vector.x,
                            transform.precisionBits(), numeric::RoundingMode::NearestEven),
                        numeric::add(
                            segment.ellipseGeometry->center.y, vector.y,
                            transform.precisionBits(), numeric::RoundingMode::NearestEven)};
                };
                const auto cosinePoint = mapDataPoint(
                    offsetEndpoint(segment.ellipseGeometry->cosineAxis));
                const auto sinePoint = mapDataPoint(
                    offsetEndpoint(segment.ellipseGeometry->sineAxis));
                if (!center || !cosinePoint || !sinePoint) {
                    result.status = PlotGraphicsLoweringStatus::MappingFailed;
                    return result;
                }
                GraphicsEllipseNode ellipse;
                ellipse.center = *center;
                ellipse.cosineAxis = {
                    cosinePoint->xMm - center->xMm, cosinePoint->yMm - center->yMm};
                ellipse.sineAxis = {
                    sinePoint->xMm - center->xMm, sinePoint->yMm - center->yMm};
                ellipse.stroke = stroke(options.curveStrokeWidthMm);
                ellipse.semantic = GraphicsSemanticRef{GraphicsSemanticKind::Curve, curve.id};
                ellipse.clipRect = graphics::GraphicsRectMm{
                    area->left, area->bottom, area->width, area->height};
                scene.nodes.push_back(std::move(ellipse));
                continue;
            }

            if (segment.geometryKind == PlotSegmentGeometryKind::EllipticArc
                && segment.ellipticArcGeometry) {
                const auto mapDataPoint = [&](const PlotEndpointPoint& point)
                    -> std::optional<GraphicsPointMm> {
                    const auto mapped = mapForGraphics(transform, point.x, point.y);
                    if (!mapped)
                        return std::nullopt;
                    return fitToPlotArea(
                        {mapped->xMm, mapped->yMm}, transform.viewport(), *area);
                };
                const auto center = mapDataPoint(segment.ellipticArcGeometry->center);
                const auto offsetEndpoint = [&](const PlotEndpointPoint& vector) {
                    return PlotEndpointPoint{
                        numeric::add(
                            segment.ellipticArcGeometry->center.x, vector.x,
                            transform.precisionBits(), numeric::RoundingMode::NearestEven),
                        numeric::add(
                            segment.ellipticArcGeometry->center.y, vector.y,
                            transform.precisionBits(), numeric::RoundingMode::NearestEven)};
                };
                const auto cosinePoint = mapDataPoint(
                    offsetEndpoint(segment.ellipticArcGeometry->cosineAxis));
                const auto sinePoint = mapDataPoint(
                    offsetEndpoint(segment.ellipticArcGeometry->sineAxis));
                const double startRadians = 2.0 * std::numbers::pi
                    * rationalToDouble(segment.ellipticArcGeometry->startTurns);
                const double sweepRadians = 2.0 * std::numbers::pi
                    * rationalToDouble(segment.ellipticArcGeometry->sweepTurns);
                if (!center || !cosinePoint || !sinePoint
                    || !std::isfinite(startRadians) || !std::isfinite(sweepRadians)
                    || sweepRadians == 0.0) {
                    result.status = PlotGraphicsLoweringStatus::MappingFailed;
                    return result;
                }
                GraphicsEllipticArcNode arc;
                arc.center = *center;
                arc.cosineAxis = {
                    cosinePoint->xMm - center->xMm, cosinePoint->yMm - center->yMm};
                arc.sineAxis = {
                    sinePoint->xMm - center->xMm, sinePoint->yMm - center->yMm};
                arc.startRadians = startRadians;
                arc.sweepRadians = sweepRadians;
                arc.stroke = stroke(options.curveStrokeWidthMm);
                arc.semantic = GraphicsSemanticRef{GraphicsSemanticKind::Curve, curve.id};
                arc.clipRect = graphics::GraphicsRectMm{
                    area->left, area->bottom, area->width, area->height};
                scene.nodes.push_back(std::move(arc));
                continue;
            }

            const auto mapBezierPoint = [&](const PlotEndpointPoint& point)
                -> std::optional<GraphicsPointMm> {
                const auto mapped = mapForGraphics(transform, point.x, point.y);
                if (!mapped)
                    return std::nullopt;
                return fitToPlotArea(
                    {mapped->xMm, mapped->yMm}, transform.viewport(), *area);
            };
            const auto splitParameters = anchorParameters(segment, mappedVertices);
            if (segment.geometryKind == PlotSegmentGeometryKind::QuadraticBezier
                && segment.quadraticControlPoint && mappedVertices.size() >= 2) {
                const auto control = mapBezierPoint(*segment.quadraticControlPoint);
                if (!control) {
                    result.status = PlotGraphicsLoweringStatus::MappingFailed;
                    return result;
                }
                appendQuadraticBezierCommands(
                    path,
                    QuadraticBezierMm{mappedVertices.front(), *control, mappedVertices.back()},
                    splitParameters, *path.clipRect);
            }
            else if (segment.geometryKind == PlotSegmentGeometryKind::CubicBezier
                && segment.bezierControlPoints && mappedVertices.size() >= 2) {
                const auto control1 = mapBezierPoint(segment.bezierControlPoints->control1);
                const auto control2 = mapBezierPoint(segment.bezierControlPoints->control2);
                if (!control1 || !control2) {
                    result.status = PlotGraphicsLoweringStatus::MappingFailed;
                    return result;
                }
                appendCubicBezierCommands(
                    path,
                    CubicBezierMm{
                        mappedVertices.front(), *control1, *control2, mappedVertices.back()},
                    splitParameters, *path.clipRect);
            }
            else {
                appendClippedPolylineCommands(path, mappedVertices, *path.clipRect);
            }
            if (!path.commands.empty())
                scene.nodes.push_back(std::move(path));
            if (options.showEndpointMarkers) {
                const auto mapEndpointPoint = [&](const std::optional<PlotEndpointPoint>& point)
                    -> std::optional<GraphicsPointMm> {
                    if (!point)
                        return std::nullopt;
                    const auto mapped = mapForGraphics(transform, point->x, point->y);
                    if (!mapped)
                        return std::nullopt;
                    return fitToPlotArea(
                        {mapped->xMm, mapped->yMm}, transform.viewport(), *area);
                };
                appendSegmentEndpointMarkers(
                    endpointMarkers, curve.id, segment, mappedVertices,
                    mapEndpointPoint(segment.lowerEndpointPoint),
                    mapEndpointPoint(segment.upperEndpointPoint));
            }
        }
    }

    const double xAxisY = area->bottom
        + axes.xAxis.axisPositionMm * area->height / transform.viewport().heightMm;
    const double yAxisX = area->left
        + axes.yAxis.axisPositionMm * area->width / transform.viewport().widthMm;
    scene.nodes.push_back(line(
        {area->left, xAxisY},
        {area->left + area->width, xAxisY},
        options.axisStrokeWidthMm,
        {GraphicsSemanticKind::XAxis, 1}));
    scene.nodes.push_back(line(
        {yAxisX, area->bottom},
        {yAxisX, area->bottom + area->height},
        options.axisStrokeWidthMm,
        {GraphicsSemanticKind::YAxis, 1}));

    if (options.showMajorTicks) {
        const auto [xBefore, xAfter] = tickOffsets(
            axes.xAxis.placement, options.majorTickLengthMm);
        for (std::size_t i = 0; i < axes.xAxis.majorTicks.size(); ++i) {
            if (skipOriginTick(axes, axes.xAxis.majorTicks[i], false))
                continue;
            const auto x = area->left
                + axes.xAxis.majorTicks[i].positionMm * area->width / transform.viewport().widthMm;
            scene.nodes.push_back(line(
                {x, xAxisY + xBefore},
                {x, xAxisY + xAfter},
                options.tickStrokeWidthMm,
                {GraphicsSemanticKind::XTick, static_cast<std::uint64_t>(i + 1)}));
        }

        const auto [yBefore, yAfter] = tickOffsets(
            axes.yAxis.placement, options.majorTickLengthMm);
        for (std::size_t i = 0; i < axes.yAxis.majorTicks.size(); ++i) {
            if (skipOriginTick(axes, axes.yAxis.majorTicks[i], true))
                continue;
            const auto y = area->bottom
                + axes.yAxis.majorTicks[i].positionMm * area->height / transform.viewport().heightMm;
            scene.nodes.push_back(line(
                {yAxisX + yBefore, y},
                {yAxisX + yAfter, y},
                options.tickStrokeWidthMm,
                {GraphicsSemanticKind::YTick, static_cast<std::uint64_t>(i + 1)}));
        }
    }

    if (options.showMajorTickLabels) {
        const bool xAbove = axes.xAxis.placement == PlotAxisPlacement::MaximumEdge;
        const double xTickExtent = axes.xAxis.placement == PlotAxisPlacement::CrossZero
            ? options.majorTickLengthMm / 2.0 : options.majorTickLengthMm;
        for (std::size_t i = 0; i < axes.xAxis.majorTicks.size(); ++i) {
            if (skipOriginTick(axes, axes.xAxis.majorTicks[i], false))
                continue;
            const auto& metrics = xLabels[i].metrics;
            GraphicsTextNode text;
            text.origin.xMm = area->left
                + axes.xAxis.majorTicks[i].positionMm * area->width / transform.viewport().widthMm;
            text.origin.yMm = xAbove
                ? xAxisY + xTickExtent + options.tickLabelGapMm + metrics.descentMm
                : xAxisY - xTickExtent - options.tickLabelGapMm - metrics.ascentMm;
            text.text = xLabels[i].text;
            text.fontSizeMm = options.tickLabelFontSizeMm;
            text.anchor = GraphicsTextAnchor::Middle;
            text.semantic = GraphicsSemanticRef{
                GraphicsSemanticKind::XTickLabel, static_cast<std::uint64_t>(i + 1)};
            scene.nodes.push_back(std::move(text));
        }

        const bool yRight = axes.yAxis.placement == PlotAxisPlacement::MaximumEdge;
        const double yTickExtent = axes.yAxis.placement == PlotAxisPlacement::CrossZero
            ? options.majorTickLengthMm / 2.0 : options.majorTickLengthMm;
        for (std::size_t i = 0; i < axes.yAxis.majorTicks.size(); ++i) {
            const auto& tick = axes.yAxis.majorTicks[i];
            if (skipOriginTick(axes, tick, true))
                continue;
            const auto& metrics = yLabels[i].metrics;
            GraphicsTextNode text;
            text.origin.xMm = yRight
                ? yAxisX + yTickExtent + options.tickLabelGapMm
                : yAxisX - yTickExtent - options.tickLabelGapMm;
            text.origin.yMm = area->bottom
                + tick.positionMm * area->height / transform.viewport().heightMm
                - (metrics.ascentMm - metrics.descentMm) / 2.0;
            text.text = yLabels[i].text;
            text.fontSizeMm = options.tickLabelFontSizeMm;
            text.anchor = yRight ? GraphicsTextAnchor::Start : GraphicsTextAnchor::End;
            text.semantic = GraphicsSemanticRef{
                GraphicsSemanticKind::YTickLabel, static_cast<std::uint64_t>(i + 1)};
            scene.nodes.push_back(std::move(text));
        }
    }

    const graphics::GraphicsRectMm curveClip{
        area->left, area->bottom, area->width, area->height};
    for (const auto& marker : endpointMarkers) {
        if (markerIntersectsClip(marker, options, curveClip))
            scene.nodes.push_back(endpointCircle(marker, options, curveClip));
    }

    result.status = PlotGraphicsLoweringStatus::Success;
    result.scene = std::move(scene);
    return result;
}

} // namespace mmcal::plot
