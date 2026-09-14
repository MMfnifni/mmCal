#include "raster_backend.hpp"

#include "text_metrics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcal::graphics {
namespace {

constexpr std::uint64_t maxInternalPixels = 128ull * 1024ull * 1024ull;
constexpr std::uint32_t maxOutputDimension = 32768;
constexpr double rasterFlatnessPx = 0.25;

struct PremulPixel final {
    std::uint8_t red = 0;
    std::uint8_t green = 0;
    std::uint8_t blue = 0;
    std::uint8_t alpha = 0;
};

struct PointPx final {
    double x = 0.0;
    double y = 0.0;
};

struct RectPx final {
    double left = 0.0;
    double top = 0.0;
    double right = 0.0;
    double bottom = 0.0;
};

struct Surface final {
    std::uint32_t width = 0;
    std::uint32_t height = 0;
    std::vector<PremulPixel> pixels;

    [[nodiscard]] PremulPixel& at(std::uint32_t x, std::uint32_t y) noexcept {
        return pixels[static_cast<std::size_t>(y) * width + x];
    }
};

struct RasterContext final {
    const GraphicsScene& scene;
    Surface surface;
    double scale = 1.0;
    double heightMm = 0.0;
};

[[nodiscard]] bool finite(GraphicsPointMm point) noexcept {
    return std::isfinite(point.xMm) && std::isfinite(point.yMm);
}

[[nodiscard]] bool validRect(const GraphicsRectMm& rect) noexcept {
    return std::isfinite(rect.xMm) && std::isfinite(rect.yMm)
        && std::isfinite(rect.widthMm) && std::isfinite(rect.heightMm)
        && rect.widthMm > 0.0 && rect.heightMm > 0.0;
}

[[nodiscard]] PointPx toPx(const RasterContext& context, GraphicsPointMm point) noexcept {
    return {point.xMm * context.scale, (context.heightMm - point.yMm) * context.scale};
}

[[nodiscard]] RectPx toPx(const RasterContext& context, const GraphicsRectMm& rect) noexcept {
    return {
        rect.xMm * context.scale,
        (context.heightMm - rect.yMm - rect.heightMm) * context.scale,
        (rect.xMm + rect.widthMm) * context.scale,
        (context.heightMm - rect.yMm) * context.scale};
}

[[nodiscard]] bool insideClip(double x, double y, const std::optional<RectPx>& clip) noexcept {
    if (!clip)
        return true;
    return x >= clip->left && x < clip->right && y >= clip->top && y < clip->bottom;
}

void blend(PremulPixel& destination, const GraphicsColor& source) noexcept {
    const std::uint32_t alpha = source.alpha;
    const std::uint32_t inverse = 255u - alpha;
    const auto sourcePremul = [alpha](std::uint8_t channel) {
        return (static_cast<std::uint32_t>(channel) * alpha + 127u) / 255u;
    };
    const auto blendChannel = [inverse](std::uint32_t src, std::uint8_t dst) {
        return static_cast<std::uint8_t>(std::min<std::uint32_t>(
            255u, src + (static_cast<std::uint32_t>(dst) * inverse + 127u) / 255u));
    };

    destination.red = blendChannel(sourcePremul(source.red), destination.red);
    destination.green = blendChannel(sourcePremul(source.green), destination.green);
    destination.blue = blendChannel(sourcePremul(source.blue), destination.blue);
    destination.alpha = static_cast<std::uint8_t>(std::min<std::uint32_t>(
        255u, alpha + (static_cast<std::uint32_t>(destination.alpha) * inverse + 127u) / 255u));
}

void blendPixel(
    RasterContext& context,
    int x,
    int y,
    const GraphicsColor& color,
    const std::optional<RectPx>& clip) noexcept {
    if (x < 0 || y < 0
        || x >= static_cast<int>(context.surface.width)
        || y >= static_cast<int>(context.surface.height))
        return;
    const double centerX = static_cast<double>(x) + 0.5;
    const double centerY = static_cast<double>(y) + 0.5;
    if (!insideClip(centerX, centerY, clip))
        return;
    blend(context.surface.at(static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(y)), color);
}

[[nodiscard]] double squaredDistance(PointPx lhs, PointPx rhs) noexcept {
    const double dx = lhs.x - rhs.x;
    const double dy = lhs.y - rhs.y;
    return dx * dx + dy * dy;
}

void drawDisk(
    RasterContext& context,
    PointPx center,
    double radius,
    const GraphicsColor& color,
    const std::optional<RectPx>& clip) {
    if (!(radius >= 0.0) || !std::isfinite(radius))
        return;
    const int minX = static_cast<int>(std::floor(center.x - radius - 1.0));
    const int maxX = static_cast<int>(std::ceil(center.x + radius + 1.0));
    const int minY = static_cast<int>(std::floor(center.y - radius - 1.0));
    const int maxY = static_cast<int>(std::ceil(center.y + radius + 1.0));
    const double radiusSquared = radius * radius;
    for (int y = minY; y <= maxY; ++y) {
        for (int x = minX; x <= maxX; ++x) {
            const PointPx sample{static_cast<double>(x) + 0.5, static_cast<double>(y) + 0.5};
            if (squaredDistance(sample, center) <= radiusSquared)
                blendPixel(context, x, y, color, clip);
        }
    }
}

void drawCapsule(
    RasterContext& context,
    PointPx from,
    PointPx to,
    double width,
    GraphicsLineCap cap,
    const GraphicsColor& color,
    const std::optional<RectPx>& clip) {
    if (!(width > 0.0) || !std::isfinite(width))
        return;
    const double radius = width / 2.0;
    const double dx = to.x - from.x;
    const double dy = to.y - from.y;
    const double lengthSquared = dx * dx + dy * dy;
    if (!(lengthSquared > 1e-24)) {
        if (cap == GraphicsLineCap::Round)
            drawDisk(context, from, radius, color, clip);
        return;
    }

    PointPx boundFrom = from;
    PointPx boundTo = to;
    if (cap == GraphicsLineCap::Square) {
        const double length = std::sqrt(lengthSquared);
        const double ux = dx / length;
        const double uy = dy / length;
        boundFrom = {from.x - ux * radius, from.y - uy * radius};
        boundTo = {to.x + ux * radius, to.y + uy * radius};
    }

    const int minX = static_cast<int>(std::floor(std::min(boundFrom.x, boundTo.x) - radius - 1.0));
    const int maxX = static_cast<int>(std::ceil(std::max(boundFrom.x, boundTo.x) + radius + 1.0));
    const int minY = static_cast<int>(std::floor(std::min(boundFrom.y, boundTo.y) - radius - 1.0));
    const int maxY = static_cast<int>(std::ceil(std::max(boundFrom.y, boundTo.y) + radius + 1.0));
    const double radiusSquared = radius * radius;

    for (int y = minY; y <= maxY; ++y) {
        for (int x = minX; x <= maxX; ++x) {
            const PointPx sample{static_cast<double>(x) + 0.5, static_cast<double>(y) + 0.5};
            double t = ((sample.x - from.x) * dx + (sample.y - from.y) * dy) / lengthSquared;
            if (cap == GraphicsLineCap::Butt && (t < 0.0 || t > 1.0))
                continue;
            if (cap == GraphicsLineCap::Square) {
                const double length = std::sqrt(lengthSquared);
                const double extension = radius / length;
                if (t < -extension || t > 1.0 + extension)
                    continue;
                t = std::clamp(t, -extension, 1.0 + extension);
            }
            else
                t = std::clamp(t, 0.0, 1.0);
            const PointPx closest{from.x + t * dx, from.y + t * dy};
            if (squaredDistance(sample, closest) <= radiusSquared)
                blendPixel(context, x, y, color, clip);
        }
    }
}

[[nodiscard]] double pointLineDistance(PointPx point, PointPx from, PointPx to) noexcept {
    const double dx = to.x - from.x;
    const double dy = to.y - from.y;
    const double lengthSquared = dx * dx + dy * dy;
    if (!(lengthSquared > 1e-24))
        return std::sqrt(squaredDistance(point, from));
    const double areaTwice = std::abs(dx * (from.y - point.y) - (from.x - point.x) * dy);
    return areaTwice / std::sqrt(lengthSquared);
}

[[nodiscard]] PointPx midpoint(PointPx lhs, PointPx rhs) noexcept {
    return {(lhs.x + rhs.x) * 0.5, (lhs.y + rhs.y) * 0.5};
}

void flattenQuadratic(
    PointPx p0,
    PointPx p1,
    PointPx p2,
    std::vector<PointPx>& output,
    std::size_t depth = 0) {
    constexpr std::size_t maxDepth = 18;
    if (depth >= maxDepth || pointLineDistance(p1, p0, p2) <= rasterFlatnessPx) {
        output.push_back(p2);
        return;
    }
    const PointPx p01 = midpoint(p0, p1);
    const PointPx p12 = midpoint(p1, p2);
    const PointPx p012 = midpoint(p01, p12);
    flattenQuadratic(p0, p01, p012, output, depth + 1);
    flattenQuadratic(p012, p12, p2, output, depth + 1);
}

void flattenCubic(
    PointPx p0,
    PointPx p1,
    PointPx p2,
    PointPx p3,
    std::vector<PointPx>& output,
    std::size_t depth = 0) {
    constexpr std::size_t maxDepth = 18;
    const double flatness = std::max(
        pointLineDistance(p1, p0, p3),
        pointLineDistance(p2, p0, p3));
    if (depth >= maxDepth || flatness <= rasterFlatnessPx) {
        output.push_back(p3);
        return;
    }
    const PointPx p01 = midpoint(p0, p1);
    const PointPx p12 = midpoint(p1, p2);
    const PointPx p23 = midpoint(p2, p3);
    const PointPx p012 = midpoint(p01, p12);
    const PointPx p123 = midpoint(p12, p23);
    const PointPx p0123 = midpoint(p012, p123);
    flattenCubic(p0, p01, p012, p0123, output, depth + 1);
    flattenCubic(p0123, p123, p23, p3, output, depth + 1);
}

struct FlattenedSubpath final {
    std::vector<PointPx> points;
    bool closed = false;
};

[[nodiscard]] bool flattenPath(
    const RasterContext& context,
    const GraphicsPathNode& path,
    std::vector<FlattenedSubpath>& subpaths) {
    std::optional<PointPx> current;
    std::optional<PointPx> subpathStart;
    FlattenedSubpath* subpath = nullptr;

    for (const auto& command : path.commands) {
        const bool ok = std::visit([&](const auto& value) -> bool {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, GraphicsMoveTo>) {
                if (!finite(value.point))
                    return false;
                subpaths.push_back(FlattenedSubpath{});
                subpath = &subpaths.back();
                current = toPx(context, value.point);
                subpathStart = current;
                subpath->points.push_back(*current);
            }
            else if constexpr (std::is_same_v<T, GraphicsLineTo>) {
                if (!current || !finite(value.point))
                    return false;
                const PointPx endpoint = toPx(context, value.point);
                subpath->points.push_back(endpoint);
                current = endpoint;
            }
            else if constexpr (std::is_same_v<T, GraphicsQuadraticTo>) {
                if (!current || !finite(value.control) || !finite(value.point))
                    return false;
                const PointPx control = toPx(context, value.control);
                const PointPx endpoint = toPx(context, value.point);
                flattenQuadratic(*current, control, endpoint, subpath->points);
                current = endpoint;
            }
            else if constexpr (std::is_same_v<T, GraphicsCubicTo>) {
                if (!current || !finite(value.control1) || !finite(value.control2)
                    || !finite(value.point))
                    return false;
                const PointPx control1 = toPx(context, value.control1);
                const PointPx control2 = toPx(context, value.control2);
                const PointPx endpoint = toPx(context, value.point);
                flattenCubic(*current, control1, control2, endpoint, subpath->points);
                current = endpoint;
            }
            else if constexpr (std::is_same_v<T, GraphicsClosePath>) {
                if (!current || !subpathStart || !subpath)
                    return false;
                if (squaredDistance(*current, *subpathStart) > 1e-18)
                    subpath->points.push_back(*subpathStart);
                subpath->closed = true;
                current = subpathStart;
            }
            return true;
        }, command);
        if (!ok)
            return false;
    }
    return true;
}

[[nodiscard]] bool pointInPolygon(PointPx point, const std::vector<PointPx>& polygon) noexcept {
    if (polygon.size() < 3)
        return false;
    bool inside = false;
    std::size_t j = polygon.size() - 1;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const PointPx a = polygon[i];
        const PointPx b = polygon[j];
        const bool crosses = (a.y > point.y) != (b.y > point.y);
        if (crosses) {
            const double x = (b.x - a.x) * (point.y - a.y) / (b.y - a.y) + a.x;
            if (point.x < x)
                inside = !inside;
        }
        j = i;
    }
    return inside;
}

void fillSubpaths(
    RasterContext& context,
    const std::vector<FlattenedSubpath>& subpaths,
    const GraphicsColor& color,
    const std::optional<RectPx>& clip) {
    bool any = false;
    double minX = std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    for (const auto& subpath : subpaths) {
        for (const auto point : subpath.points) {
            any = true;
            minX = std::min(minX, point.x);
            minY = std::min(minY, point.y);
            maxX = std::max(maxX, point.x);
            maxY = std::max(maxY, point.y);
        }
    }
    if (!any)
        return;

    const int x0 = static_cast<int>(std::floor(minX));
    const int x1 = static_cast<int>(std::ceil(maxX));
    const int y0 = static_cast<int>(std::floor(minY));
    const int y1 = static_cast<int>(std::ceil(maxY));
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            const PointPx sample{static_cast<double>(x) + 0.5, static_cast<double>(y) + 0.5};
            bool inside = false;
            for (const auto& subpath : subpaths)
                inside ^= pointInPolygon(sample, subpath.points);
            if (inside)
                blendPixel(context, x, y, color, clip);
        }
    }
}

void strokeSubpaths(
    RasterContext& context,
    const std::vector<FlattenedSubpath>& subpaths,
    const GraphicsStrokeStyle& stroke,
    const std::optional<RectPx>& clip) {
    const double width = stroke.widthMm * context.scale;
    if (!(width > 0.0) || !std::isfinite(width))
        return;
    for (const auto& subpath : subpaths) {
        if (subpath.points.size() < 2)
            continue;
        for (std::size_t i = 1; i < subpath.points.size(); ++i) {
            GraphicsLineCap cap = stroke.lineCap;
            if (stroke.lineJoin == GraphicsLineJoin::Round)
                cap = GraphicsLineCap::Round;
            drawCapsule(
                context, subpath.points[i - 1], subpath.points[i], width,
                cap, stroke.color, clip);
        }
        if (stroke.lineJoin == GraphicsLineJoin::Round && subpath.points.size() > 2) {
            const double radius = width / 2.0;
            for (std::size_t i = 1; i + 1 < subpath.points.size(); ++i)
                drawDisk(context, subpath.points[i], radius, stroke.color, clip);
        }
    }
}

[[nodiscard]] std::vector<PointPx> flattenEllipse(
    const RasterContext& context,
    GraphicsPointMm center,
    GraphicsPointMm cosineAxis,
    GraphicsPointMm sineAxis,
    double startRadians,
    double sweepRadians) {
    const PointPx c = toPx(context, center);
    const double ax = cosineAxis.xMm * context.scale;
    const double ay = -cosineAxis.yMm * context.scale;
    const double bx = sineAxis.xMm * context.scale;
    const double by = -sineAxis.yMm * context.scale;
    const double radiusBound = std::max(std::hypot(ax, ay), std::hypot(bx, by));
    const double sweep = std::abs(sweepRadians);
    std::size_t segments = 1;
    if (radiusBound > rasterFlatnessPx && sweep > 0.0) {
        const double ratio = std::clamp(1.0 - rasterFlatnessPx / radiusBound, -1.0, 1.0);
        const double maxStep = 2.0 * std::acos(ratio);
        if (maxStep > 1e-9)
            segments = static_cast<std::size_t>(std::ceil(sweep / maxStep));
    }
    segments = std::clamp<std::size_t>(segments, 1, 32768);

    std::vector<PointPx> points;
    points.reserve(segments + 1);
    for (std::size_t i = 0; i <= segments; ++i) {
        const double t = startRadians
            + sweepRadians * static_cast<double>(i) / static_cast<double>(segments);
        const double cosine = std::cos(t);
        const double sine = std::sin(t);
        points.push_back({
            c.x + ax * cosine + bx * sine,
            c.y + ay * cosine + by * sine});
    }
    return points;
}

[[nodiscard]] std::optional<RectPx> nodeClip(
    const RasterContext& context,
    const std::optional<GraphicsRectMm>& rect) {
    if (!rect)
        return std::nullopt;
    if (!validRect(*rect))
        return std::nullopt;
    return toPx(context, *rect);
}

[[nodiscard]] bool renderPath(RasterContext& context, const GraphicsPathNode& path) {
    if (path.clipRect && !validRect(*path.clipRect))
        return false;
    if (path.stroke && (!(path.stroke->widthMm > 0.0) || !std::isfinite(path.stroke->widthMm)))
        return false;
    std::vector<FlattenedSubpath> subpaths;
    if (!flattenPath(context, path, subpaths))
        return false;
    const auto clip = nodeClip(context, path.clipRect);
    if (path.fill)
        fillSubpaths(context, subpaths, path.fill->color, clip);
    if (path.stroke)
        strokeSubpaths(context, subpaths, *path.stroke, clip);
    return true;
}

[[nodiscard]] bool renderEllipse(
    RasterContext& context,
    const GraphicsEllipseNode& ellipse) {
    if (!finite(ellipse.center) || !finite(ellipse.cosineAxis) || !finite(ellipse.sineAxis)
        || (ellipse.clipRect && !validRect(*ellipse.clipRect)))
        return false;
    const auto points = flattenEllipse(
        context, ellipse.center, ellipse.cosineAxis, ellipse.sineAxis,
        0.0, 2.0 * std::numbers::pi);
    FlattenedSubpath subpath{points, true};
    if (!subpath.points.empty())
        subpath.points.push_back(subpath.points.front());
    const std::vector<FlattenedSubpath> subpaths{subpath};
    const auto clip = nodeClip(context, ellipse.clipRect);
    if (ellipse.fill)
        fillSubpaths(context, subpaths, ellipse.fill->color, clip);
    if (ellipse.stroke)
        strokeSubpaths(context, subpaths, *ellipse.stroke, clip);
    return true;
}

[[nodiscard]] bool renderEllipticArc(
    RasterContext& context,
    const GraphicsEllipticArcNode& arc) {
    if (!finite(arc.center) || !finite(arc.cosineAxis) || !finite(arc.sineAxis)
        || !std::isfinite(arc.startRadians) || !std::isfinite(arc.sweepRadians)
        || (arc.clipRect && !validRect(*arc.clipRect)))
        return false;
    FlattenedSubpath subpath{
        flattenEllipse(
            context, arc.center, arc.cosineAxis, arc.sineAxis,
            arc.startRadians, arc.sweepRadians),
        false};
    const std::vector<FlattenedSubpath> subpaths{subpath};
    const auto clip = nodeClip(context, arc.clipRect);
    if (arc.fill)
        fillSubpaths(context, subpaths, arc.fill->color, clip);
    if (arc.stroke)
        strokeSubpaths(context, subpaths, *arc.stroke, clip);
    return true;
}

[[nodiscard]] bool renderCircle(RasterContext& context, const GraphicsCircleNode& circle) {
    if (!finite(circle.center) || !(circle.radiusMm >= 0.0) || !std::isfinite(circle.radiusMm)
        || (circle.clipRect && !validRect(*circle.clipRect)))
        return false;
    const PointPx center = toPx(context, circle.center);
    const double radius = circle.radiusMm * context.scale;
    const auto clip = nodeClip(context, circle.clipRect);
    if (circle.fill)
        drawDisk(context, center, radius, circle.fill->color, clip);
    if (circle.stroke) {
        if (!(circle.stroke->widthMm > 0.0) || !std::isfinite(circle.stroke->widthMm))
            return false;
        const double width = circle.stroke->widthMm * context.scale;
        const double outer = radius + width / 2.0;
        const double inner = std::max(0.0, radius - width / 2.0);
        const int minX = static_cast<int>(std::floor(center.x - outer - 1.0));
        const int maxX = static_cast<int>(std::ceil(center.x + outer + 1.0));
        const int minY = static_cast<int>(std::floor(center.y - outer - 1.0));
        const int maxY = static_cast<int>(std::ceil(center.y + outer + 1.0));
        const double outerSquared = outer * outer;
        const double innerSquared = inner * inner;
        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                const PointPx sample{static_cast<double>(x) + 0.5, static_cast<double>(y) + 0.5};
                const double distanceSquared = squaredDistance(sample, center);
                if (distanceSquared <= outerSquared && distanceSquared >= innerSquared)
                    blendPixel(context, x, y, circle.stroke->color, clip);
            }
        }
    }
    return true;
}

struct StrokePoint final {
    double x = 0.0;
    double y = 0.0;
};

struct StrokeSegment final {
    StrokePoint from;
    StrokePoint to;
};

struct StrokeGlyph final {
    std::span<const StrokeSegment> segments;
    // 添付の字形比率を保つため，高さに対する幅で保持する。
    double widthInHeights = 0.0;
};

constexpr std::array<StrokeSegment, 1> glyphMinus{{
    // 線端が隣接数字へ食い込まないよう，字幅内に左右のside bearingを持たせる。
    {{0.064819, 0.500000}, {0.574819, 0.500000}},
}};

constexpr std::array<StrokeSegment, 8> glyph0{{
    {{0.477483, 0.000000}, {0.594595, 0.153155}},
    {{0.594595, 0.153155}, {0.594595, 0.855862}},
    {{0.594595, 0.855862}, {0.477483, 1.000000}},
    {{0.477483, 0.000000}, {0.117121, 0.000000}},
    {{0.117121, 0.000000}, {0.000000, 0.153155}},
    {{0.000000, 0.153155}, {0.000000, 0.855862}},
    {{0.000000, 0.855862}, {0.117121, 1.000000}},
    {{0.117121, 1.000000}, {0.477483, 1.000000}},
}};

constexpr std::array<StrokeSegment, 3> glyph1{{
    {{0.000000, 0.855862}, {0.081087, 1.000000}},
    {{0.081087, 1.000000}, {0.081087, 0.000000}},
    {{0.000000, 0.000000}, {0.162164, 0.000000}},
}};

constexpr std::array<StrokeSegment, 14> glyph2{{
    {{0.000000, 0.702707}, {0.000000, 0.855862}},
    {{0.000000, 0.855862}, {0.081077, 0.945948}},
    {{0.081077, 0.945948}, {0.162164, 1.000000}},
    {{0.162164, 1.000000}, {0.396397, 1.000000}},
    {{0.441440, 1.000000}, {0.522517, 0.945948}},
    {{0.522517, 0.945948}, {0.594586, 0.855862}},
    {{0.594586, 0.855862}, {0.594586, 0.603603}},
    {{0.594586, 0.603603}, {0.441440, 0.504500}},
    {{0.441440, 0.504500}, {0.234233, 0.450448}},
    {{0.234233, 0.450448}, {0.117112, 0.351344}},
    {{0.117112, 0.351344}, {0.045043, 0.207207}},
    {{0.045043, 0.207207}, {0.000000, 0.000000}},
    {{0.000000, 0.000000}, {0.594586, 0.000000}},
    {{0.396397, 1.000000}, {0.441440, 1.000000}},
}};

constexpr std::array<StrokeSegment, 17> glyph3{{
    {{0.396397, 0.549552}, {0.198198, 0.549552}},
    {{0.000000, 0.846845}, {0.072073, 0.945948}},
    {{0.072073, 0.945948}, {0.153151, 0.990982}},
    {{0.153151, 0.990982}, {0.396397, 0.990982}},
    {{0.396397, 0.990982}, {0.513513, 0.945948}},
    {{0.513513, 0.945948}, {0.594590, 0.846845}},
    {{0.594590, 0.846845}, {0.594590, 0.693689}},
    {{0.594590, 0.693689}, {0.513513, 0.594604}},
    {{0.513513, 0.594604}, {0.396397, 0.549552}},
    {{0.396397, 0.549552}, {0.513513, 0.495500}},
    {{0.513513, 0.495500}, {0.594590, 0.396397}},
    {{0.594590, 0.396397}, {0.594590, 0.153155}},
    {{0.594590, 0.153155}, {0.513513, 0.054052}},
    {{0.513513, 0.054052}, {0.396397, 0.000000}},
    {{0.396397, 0.000000}, {0.153151, 0.000000}},
    {{0.153151, 0.000000}, {0.036034, 0.054052}},
    {{0.036034, 0.054052}, {0.000000, 0.153155}},
}};

constexpr std::array<StrokeSegment, 3> glyph4{{
    {{0.477483, 0.990987}, {0.000000, 0.243241}},
    {{0.000000, 0.243241}, {0.594586, 0.243241}},
    {{0.477483, 0.000000}, {0.477483, 0.990987}},
}};

constexpr std::array<StrokeSegment, 13> glyph5{{
    {{0.594595, 0.504500}, {0.594595, 0.153155}},
    {{0.558560, 1.000000}, {0.000000, 1.000000}},
    {{0.000000, 1.000000}, {0.000000, 0.504500}},
    {{0.000000, 0.504500}, {0.081077, 0.603603}},
    {{0.081077, 0.603603}, {0.198198, 0.648656}},
    {{0.198198, 0.648656}, {0.396397, 0.648656}},
    {{0.594595, 0.153155}, {0.513508, 0.054052}},
    {{0.513508, 0.054052}, {0.396397, 0.000000}},
    {{0.396397, 0.000000}, {0.198198, 0.000000}},
    {{0.198198, 0.000000}, {0.081077, 0.054052}},
    {{0.081077, 0.054052}, {0.000000, 0.108103}},
    {{0.396397, 0.648656}, {0.513508, 0.603603}},
    {{0.513508, 0.603603}, {0.594595, 0.504500}},
}};

constexpr std::array<StrokeSegment, 17> glyph6{{
    {{0.396397, 0.549570}, {0.189190, 0.549570}},
    {{0.000000, 0.846863}, {0.072078, 0.945948}},
    {{0.072078, 0.945948}, {0.153155, 0.991001}},
    {{0.153155, 0.991001}, {0.396397, 0.991001}},
    {{0.396397, 0.991001}, {0.513517, 0.945948}},
    {{0.513517, 0.945948}, {0.594595, 0.846863}},
    {{0.396397, 0.549570}, {0.513517, 0.504518}},
    {{0.513517, 0.504518}, {0.594595, 0.396415}},
    {{0.594595, 0.396415}, {0.594595, 0.153155}},
    {{0.594595, 0.153155}, {0.513517, 0.054070}},
    {{0.513517, 0.054070}, {0.396397, 0.000000}},
    {{0.396397, 0.000000}, {0.189190, 0.000000}},
    {{0.189190, 0.000000}, {0.072078, 0.054070}},
    {{0.072078, 0.054070}, {0.000000, 0.153155}},
    {{0.000000, 0.153155}, {0.000000, 0.846863}},
    {{0.189190, 0.549570}, {0.072078, 0.504518}},
    {{0.072078, 0.504518}, {0.000000, 0.396415}},
}};

constexpr std::array<StrokeSegment, 6> glyph7{{
    {{0.000000, 1.000000}, {0.594586, 1.000000}},
    {{0.594586, 1.000000}, {0.594586, 0.891879}},
    {{0.594586, 0.891879}, {0.396387, 0.693689}},
    {{0.396387, 0.693689}, {0.315310, 0.549552}},
    {{0.315310, 0.549552}, {0.279276, 0.396397}},
    {{0.279276, 0.396397}, {0.279276, 0.000000}},
}};

constexpr std::array<StrokeSegment, 23> glyph8{{
    {{0.153155, 0.549534}, {0.036034, 0.594586}},
    {{0.036034, 0.594586}, {0.000000, 0.702689}},
    {{0.000000, 0.702689}, {0.000000, 0.846845}},
    {{0.000000, 0.846845}, {0.036034, 0.945930}},
    {{0.036034, 0.945930}, {0.153155, 0.999982}},
    {{0.153155, 0.999982}, {0.432431, 0.999982}},
    {{0.432431, 0.999982}, {0.558570, 0.945930}},
    {{0.558570, 0.945930}, {0.594604, 0.846845}},
    {{0.594604, 0.846845}, {0.594604, 0.702689}},
    {{0.594604, 0.702689}, {0.558570, 0.594586}},
    {{0.558570, 0.594586}, {0.432431, 0.549534}},
    {{0.432431, 0.549534}, {0.153155, 0.549534}},
    {{0.153155, 0.549534}, {0.036034, 0.450448}},
    {{0.036034, 0.450448}, {0.000000, 0.351344}},
    {{0.000000, 0.351344}, {0.000000, 0.153155}},
    {{0.000000, 0.153155}, {0.036034, 0.054052}},
    {{0.036034, 0.054052}, {0.153155, 0.000000}},
    {{0.153155, 0.000000}, {0.432431, 0.000000}},
    {{0.432431, 0.000000}, {0.558570, 0.054052}},
    {{0.558570, 0.054052}, {0.594604, 0.153155}},
    {{0.594604, 0.153155}, {0.594604, 0.351344}},
    {{0.594604, 0.351344}, {0.558570, 0.450448}},
    {{0.558570, 0.450448}, {0.432431, 0.549534}},
}};

constexpr std::array<StrokeSegment, 17> glyph9{{
    {{0.198198, 0.441430}, {0.405405, 0.441430}},
    {{0.594595, 0.144138}, {0.522517, 0.045052}},
    {{0.522517, 0.045052}, {0.441440, 0.000000}},
    {{0.441440, 0.000000}, {0.198198, 0.000000}},
    {{0.198198, 0.000000}, {0.081077, 0.045052}},
    {{0.081077, 0.045052}, {0.000000, 0.144138}},
    {{0.198198, 0.441430}, {0.081077, 0.486483}},
    {{0.081077, 0.486483}, {0.000000, 0.594586}},
    {{0.000000, 0.594586}, {0.000000, 0.837845}},
    {{0.000000, 0.837845}, {0.081077, 0.936931}},
    {{0.081077, 0.936931}, {0.198198, 0.991001}},
    {{0.198198, 0.991001}, {0.405405, 0.991001}},
    {{0.405405, 0.991001}, {0.522517, 0.936931}},
    {{0.522517, 0.936931}, {0.594595, 0.837845}},
    {{0.594595, 0.837845}, {0.594595, 0.144138}},
    {{0.405405, 0.441430}, {0.522517, 0.486483}},
    {{0.522517, 0.486483}, {0.594595, 0.594586}},
}};

constexpr std::array<StrokeSegment, 11> glyphPeriod{{
    {{0.045052, 0.000000}, {0.126139, 0.000000}},
    {{0.171172, 0.207189}, {0.171172, 0.054033}},
    {{0.126139, 0.000000}, {0.171172, 0.054033}},
    {{0.045052, 0.000000}, {0.000000, 0.054033}},
    {{0.000000, 0.054033}, {0.000000, 0.207189}},
    {{0.126139, 0.261258}, {0.045052, 0.261258}},
    {{0.045052, 0.261258}, {0.000000, 0.207189}},
    {{0.171172, 0.207189}, {0.126139, 0.261258}},
    {{0.000000, 0.054033}, {0.171172, 0.207189}},
    {{0.171172, 0.054033}, {0.000000, 0.207189}},
    {{0.085595, 0.261258}, {0.085595, 0.000000}},
}};


template <std::size_t N>
[[nodiscard]] constexpr StrokeGlyph makeStrokeGlyph(
    const std::array<StrokeSegment, N>& segments,
    double widthInHeights) noexcept {
    return StrokeGlyph{std::span<const StrokeSegment>{segments}, widthInHeights};
}

[[nodiscard]] std::optional<StrokeGlyph> glyph(char ch) noexcept {
    switch (ch) {
    case '0': return makeStrokeGlyph(glyph0, 0.594595);
    case '1': return makeStrokeGlyph(glyph1, 0.162164);
    case '2': return makeStrokeGlyph(glyph2, 0.594586);
    case '3': return makeStrokeGlyph(glyph3, 0.594590);
    case '4': return makeStrokeGlyph(glyph4, 0.594586);
    case '5': return makeStrokeGlyph(glyph5, 0.594595);
    case '6': return makeStrokeGlyph(glyph6, 0.594595);
    case '7': return makeStrokeGlyph(glyph7, 0.594586);
    case '8': return makeStrokeGlyph(glyph8, 0.594604);
    case '9': return makeStrokeGlyph(glyph9, 0.594595);
    case '-': return makeStrokeGlyph(glyphMinus, 0.639638);
    case '.': return makeStrokeGlyph(glyphPeriod, 0.171172);
    default: return std::nullopt;
    }
}

[[nodiscard]] double characterAdvanceMm(char ch, double fontSizeMm) noexcept {
    const unsigned char uch = static_cast<unsigned char>(ch);
    if (uch >= '0' && uch <= '9')
        return 0.56 * fontSizeMm;
    if (ch == '.')
        return 0.28 * fontSizeMm;
    if (ch == '-')
        return 0.52 * fontSizeMm;
    return 0.60 * fontSizeMm;
}

// 数字・符号・小数点の組合せで線端が詰まって見えないよう，glyph間だけ僅かなtrackingを入れる。
// 先頭・末尾には加えず，anchor位置と外形の余白を不必要に動かさない。
constexpr double strokeGlyphSpacingEm = 0.04;

[[nodiscard]] double strokeTextAdvanceMm(std::string_view text, double fontSizeMm) noexcept {
    double width = 0.0;
    for (std::size_t i = 0; i < text.size(); ++i) {
        width += characterAdvanceMm(text[i], fontSizeMm);
        if (i + 1 < text.size())
            width += strokeGlyphSpacingEm * fontSizeMm;
    }
    return width;
}

[[nodiscard]] bool renderText(RasterContext& context, const GraphicsTextNode& text) {
    if (!finite(text.origin) || !(text.fontSizeMm > 0.0) || !std::isfinite(text.fontSizeMm))
        return false;
    for (const char ch : text.text) {
        if (!glyph(ch))
            return false;
    }

    const double textAdvance = strokeTextAdvanceMm(text.text, text.fontSizeMm);
    double x = text.origin.xMm;
    if (text.anchor == GraphicsTextAnchor::Middle)
        x -= textAdvance / 2.0;
    else if (text.anchor == GraphicsTextAnchor::End)
        x -= textAdvance;

    // 数値tick専用の単線font。添付字形の縦横比を保ち，線幅もfont sizeに比例させる。
    const double glyphHeight = 0.78 * text.fontSizeMm;
    const double strokeWidthPx = 0.082 * text.fontSizeMm * context.scale;
    for (std::size_t i = 0; i < text.text.size(); ++i) {
        const char ch = text.text[i];
        const double advance = characterAdvanceMm(ch, text.fontSizeMm);
        const auto shape = glyph(ch);
        if (!shape)
            return false;
        const double drawingWidth = shape->widthInHeights * glyphHeight;
        const double left = x + (advance - drawingWidth) / 2.0;
        for (const auto& segment : shape->segments) {
            const GraphicsPointMm fromMm{
                left + segment.from.x * glyphHeight,
                text.origin.yMm + segment.from.y * glyphHeight};
            const GraphicsPointMm toMm{
                left + segment.to.x * glyphHeight,
                text.origin.yMm + segment.to.y * glyphHeight};
            drawCapsule(
                context, toPx(context, fromMm), toPx(context, toMm), strokeWidthPx,
                GraphicsLineCap::Round, text.color, std::nullopt);
        }
        x += advance;
        if (i + 1 < text.text.size())
            x += strokeGlyphSpacingEm * text.fontSizeMm;
    }
    return true;
}

[[nodiscard]] RasterRenderStatus renderNode(RasterContext& context, const GraphicsNode& node) {
    return std::visit([&](const auto& value) -> RasterRenderStatus {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, GraphicsPathNode>)
            return renderPath(context, value) ? RasterRenderStatus::Success : RasterRenderStatus::InvalidScene;
        else if constexpr (std::is_same_v<T, GraphicsEllipseNode>)
            return renderEllipse(context, value) ? RasterRenderStatus::Success : RasterRenderStatus::InvalidScene;
        else if constexpr (std::is_same_v<T, GraphicsEllipticArcNode>)
            return renderEllipticArc(context, value) ? RasterRenderStatus::Success : RasterRenderStatus::InvalidScene;
        else if constexpr (std::is_same_v<T, GraphicsCircleNode>)
            return renderCircle(context, value) ? RasterRenderStatus::Success : RasterRenderStatus::InvalidScene;
        else {
            if (renderText(context, value))
                return RasterRenderStatus::Success;
            for (const char ch : value.text) {
                if (!glyph(ch))
                    return RasterRenderStatus::UnsupportedText;
            }
            return RasterRenderStatus::InvalidScene;
        }
    }, node);
}

[[nodiscard]] std::vector<std::uint8_t> downsample(
    const Surface& source,
    std::uint32_t width,
    std::uint32_t height,
    std::uint32_t factor) {
    std::vector<std::uint8_t> rgba(
        static_cast<std::size_t>(width) * height * 4u, 0);
    const std::uint32_t sampleCount = factor * factor;
    for (std::uint32_t y = 0; y < height; ++y) {
        for (std::uint32_t x = 0; x < width; ++x) {
            std::uint32_t red = 0;
            std::uint32_t green = 0;
            std::uint32_t blue = 0;
            std::uint32_t alpha = 0;
            for (std::uint32_t sy = 0; sy < factor; ++sy) {
                for (std::uint32_t sx = 0; sx < factor; ++sx) {
                    const auto& pixel = source.pixels[
                        static_cast<std::size_t>(y * factor + sy) * source.width
                        + (x * factor + sx)];
                    red += pixel.red;
                    green += pixel.green;
                    blue += pixel.blue;
                    alpha += pixel.alpha;
                }
            }
            red = (red + sampleCount / 2u) / sampleCount;
            green = (green + sampleCount / 2u) / sampleCount;
            blue = (blue + sampleCount / 2u) / sampleCount;
            alpha = (alpha + sampleCount / 2u) / sampleCount;

            const std::size_t index = (static_cast<std::size_t>(y) * width + x) * 4u;
            rgba[index + 3] = static_cast<std::uint8_t>(alpha);
            if (alpha == 0)
                continue;
            rgba[index + 0] = static_cast<std::uint8_t>(std::min<std::uint32_t>(
                255u, (red * 255u + alpha / 2u) / alpha));
            rgba[index + 1] = static_cast<std::uint8_t>(std::min<std::uint32_t>(
                255u, (green * 255u + alpha / 2u) / alpha));
            rgba[index + 2] = static_cast<std::uint8_t>(std::min<std::uint32_t>(
                255u, (blue * 255u + alpha / 2u) / alpha));
        }
    }
    return rgba;
}

} // namespace

RasterRenderResult renderRaster(const GraphicsScene& scene, const RasterRenderOptions& options) {
    if (!(scene.extent.widthMm > 0.0) || !(scene.extent.heightMm > 0.0)
        || !std::isfinite(scene.extent.widthMm) || !std::isfinite(scene.extent.heightMm))
        return RasterRenderResult{RasterRenderStatus::InvalidScene, std::nullopt};
    if (!(options.dpi > 0.0) || !std::isfinite(options.dpi)
        || (options.antialiasing != 1 && options.antialiasing != 2 && options.antialiasing != 4)
        || ((options.widthPx == 0) != (options.heightPx == 0)))
        return RasterRenderResult{RasterRenderStatus::InvalidOptions, std::nullopt};

    std::uint32_t width = options.widthPx;
    std::uint32_t height = options.heightPx;
    double dpiX = options.dpi;
    double dpiY = options.dpi;
    if (width == 0) {
        const double widthValue = scene.extent.widthMm * options.dpi / 25.4;
        const double heightValue = scene.extent.heightMm * options.dpi / 25.4;
        if (!(widthValue >= 1.0) || !(heightValue >= 1.0)
            || widthValue > maxOutputDimension || heightValue > maxOutputDimension)
            return RasterRenderResult{RasterRenderStatus::ResourceLimit, std::nullopt};
        width = static_cast<std::uint32_t>(std::llround(widthValue));
        height = static_cast<std::uint32_t>(std::llround(heightValue));
    }
    else {
        if (width > maxOutputDimension || height > maxOutputDimension)
            return RasterRenderResult{RasterRenderStatus::ResourceLimit, std::nullopt};
        dpiX = static_cast<double>(width) * 25.4 / scene.extent.widthMm;
        dpiY = static_cast<double>(height) * 25.4 / scene.extent.heightMm;
        const double relativeDifference = std::abs(dpiX - dpiY) / std::max(dpiX, dpiY);
        if (relativeDifference > 0.002)
            return RasterRenderResult{RasterRenderStatus::InvalidOptions, std::nullopt};
    }
    width = std::max<std::uint32_t>(1, width);
    height = std::max<std::uint32_t>(1, height);

    const std::uint64_t internalWidth = static_cast<std::uint64_t>(width) * options.antialiasing;
    const std::uint64_t internalHeight = static_cast<std::uint64_t>(height) * options.antialiasing;
    if (internalWidth > maxOutputDimension * 4ull || internalHeight > maxOutputDimension * 4ull
        || internalWidth * internalHeight > maxInternalPixels)
        return RasterRenderResult{RasterRenderStatus::ResourceLimit, std::nullopt};

    try {
        const PremulPixel background = options.background == RasterBackground::White
            ? PremulPixel{255, 255, 255, 255}
            : PremulPixel{0, 0, 0, 0};
        RasterContext context{
            scene,
            Surface{
                static_cast<std::uint32_t>(internalWidth),
                static_cast<std::uint32_t>(internalHeight),
                std::vector<PremulPixel>(
                    static_cast<std::size_t>(internalWidth * internalHeight),
                    background)},
            (static_cast<double>(width) / scene.extent.widthMm) * options.antialiasing,
            scene.extent.heightMm};

        for (const auto& node : scene.nodes) {
            const RasterRenderStatus status = renderNode(context, node);
            if (status != RasterRenderStatus::Success)
                return RasterRenderResult{status, std::nullopt};
        }

        RasterImage image;
        image.widthPx = width;
        image.heightPx = height;
        image.dpiX = dpiX;
        image.dpiY = dpiY;
        image.rgba = downsample(context.surface, width, height, options.antialiasing);
        return RasterRenderResult{RasterRenderStatus::Success, std::move(image)};
    }
    catch (const std::bad_alloc&) {
        return RasterRenderResult{RasterRenderStatus::ResourceLimit, std::nullopt};
    }
    catch (...) {
        return RasterRenderResult{RasterRenderStatus::RenderFailed, std::nullopt};
    }
}

} // namespace mmcal::graphics
