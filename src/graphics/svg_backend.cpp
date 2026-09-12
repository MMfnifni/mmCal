#include "svg_backend.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>

namespace mmcal::graphics {
namespace {

[[nodiscard]] bool finite(GraphicsPointMm point) noexcept {
    return std::isfinite(point.xMm) && std::isfinite(point.yMm);
}

[[nodiscard]] std::string number(double value) {
    if (value == 0.0)
        value = 0.0;
    std::ostringstream out;
    out << std::setprecision(15) << value;
    return out.str();
}

[[nodiscard]] std::string color(const GraphicsColor& value) {
    std::ostringstream out;
    out << '#' << std::hex << std::setfill('0')
        << std::setw(2) << static_cast<unsigned>(value.red)
        << std::setw(2) << static_cast<unsigned>(value.green)
        << std::setw(2) << static_cast<unsigned>(value.blue);
    return out.str();
}

[[nodiscard]] std::string escapeXml(std::string_view text) {
    std::string out;
    out.reserve(text.size());
    for (const char ch : text) {
        switch (ch) {
        case '&': out += "&amp;"; break;
        case '<': out += "&lt;"; break;
        case '>': out += "&gt;"; break;
        case '"': out += "&quot;"; break;
        case '\'': out += "&apos;"; break;
        default: out += ch; break;
        }
    }
    return out;
}

[[nodiscard]] const char* lineCap(GraphicsLineCap value) noexcept {
    switch (value) {
    case GraphicsLineCap::Butt: return "butt";
    case GraphicsLineCap::Round: return "round";
    case GraphicsLineCap::Square: return "square";
    }
    return "butt";
}

[[nodiscard]] const char* lineJoin(GraphicsLineJoin value) noexcept {
    switch (value) {
    case GraphicsLineJoin::Miter: return "miter";
    case GraphicsLineJoin::Round: return "round";
    case GraphicsLineJoin::Bevel: return "bevel";
    }
    return "miter";
}

[[nodiscard]] const char* textAnchor(GraphicsTextAnchor value) noexcept {
    switch (value) {
    case GraphicsTextAnchor::Start: return "start";
    case GraphicsTextAnchor::Middle: return "middle";
    case GraphicsTextAnchor::End: return "end";
    }
    return "start";
}

[[nodiscard]] const char* semanticKind(GraphicsSemanticKind value) noexcept {
    switch (value) {
    case GraphicsSemanticKind::Curve: return "curve";
    case GraphicsSemanticKind::CurveEndpoint: return "curve-endpoint";
    case GraphicsSemanticKind::Point: return "point";
    case GraphicsSemanticKind::XAxis: return "x-axis";
    case GraphicsSemanticKind::YAxis: return "y-axis";
    case GraphicsSemanticKind::XTick: return "x-tick";
    case GraphicsSemanticKind::YTick: return "y-tick";
    case GraphicsSemanticKind::XTickLabel: return "x-tick-label";
    case GraphicsSemanticKind::YTickLabel: return "y-tick-label";
    case GraphicsSemanticKind::Anchor: return "anchor";
    case GraphicsSemanticKind::Label: return "label";
    }
    return "unknown";
}

void writeSemantic(std::ostringstream& out, const std::optional<GraphicsSemanticRef>& semantic) {
    if (!semantic)
        return;
    out << " data-mmcal-kind=\"" << semanticKind(semantic->kind)
        << "\" data-mmcal-id=\"" << semantic->id << '"';
}

[[nodiscard]] bool writePath(
    std::ostringstream& out,
    const GraphicsPathNode& path,
    double heightMm,
    std::optional<std::size_t> clipId,
    std::string_view indent = "  ") {
    out << indent << "<path d=\"";
    for (const auto& command : path.commands) {
        bool commandFinite = true;
        std::visit([&](const auto& value) {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, GraphicsMoveTo>) {
                commandFinite = finite(value.point);
                if (commandFinite)
                    out << 'M' << number(value.point.xMm) << ' '
                        << number(heightMm - value.point.yMm) << ' ';
            }
            else if constexpr (std::is_same_v<T, GraphicsLineTo>) {
                commandFinite = finite(value.point);
                if (commandFinite)
                    out << 'L' << number(value.point.xMm) << ' '
                        << number(heightMm - value.point.yMm) << ' ';
            }
            else if constexpr (std::is_same_v<T, GraphicsQuadraticTo>) {
                commandFinite = finite(value.control) && finite(value.point);
                if (commandFinite)
                    out << 'Q' << number(value.control.xMm) << ' '
                        << number(heightMm - value.control.yMm) << ' '
                        << number(value.point.xMm) << ' '
                        << number(heightMm - value.point.yMm) << ' ';
            }
            else if constexpr (std::is_same_v<T, GraphicsCubicTo>) {
                commandFinite = finite(value.control1) && finite(value.control2) && finite(value.point);
                if (commandFinite)
                    out << 'C' << number(value.control1.xMm) << ' '
                        << number(heightMm - value.control1.yMm) << ' '
                        << number(value.control2.xMm) << ' '
                        << number(heightMm - value.control2.yMm) << ' '
                        << number(value.point.xMm) << ' '
                        << number(heightMm - value.point.yMm) << ' ';
            }
            else if constexpr (std::is_same_v<T, GraphicsClosePath>)
                out << 'Z';
        }, command);
        if (!commandFinite)
            return false;
    }
    out << '"';
    if (path.stroke) {
        out << " stroke=\"" << color(path.stroke->color) << '"'
            << " stroke-width=\"" << number(path.stroke->widthMm) << '"'
            << " stroke-linecap=\"" << lineCap(path.stroke->lineCap) << '"'
            << " stroke-linejoin=\"" << lineJoin(path.stroke->lineJoin) << '"';
        if (path.stroke->color.alpha != 255)
            out << " stroke-opacity=\""
                << number(static_cast<double>(path.stroke->color.alpha) / 255.0) << '"';
    }
    else
        out << " stroke=\"none\"";
    if (path.fill) {
        out << " fill=\"" << color(path.fill->color) << '"';
        if (path.fill->color.alpha != 255)
            out << " fill-opacity=\""
                << number(static_cast<double>(path.fill->color.alpha) / 255.0) << '"';
    }
    else
        out << " fill=\"none\"";
    if (clipId)
        out << " clip-path=\"url(#mmcal-clip-" << *clipId << ")\"";
    writeSemantic(out, path.semantic);
    out << "/>\n";
    return true;
}


struct SvgEllipseGeometry final {
    double radiusX = 0.0;
    double radiusY = 0.0;
    double rotationRadians = 0.0;
    double determinant = 0.0;
};

[[nodiscard]] std::optional<SvgEllipseGeometry> canonicalSvgEllipse(
    GraphicsPointMm cosineAxis,
    GraphicsPointMm sineAxis) noexcept {
    // GraphicsSceneはy上向き，SVGはy下向きなので，まずSVG座標系の2x2 affine部分へ変換する。
    const double a = cosineAxis.xMm;
    const double b = sineAxis.xMm;
    const double c = -cosineAxis.yMm;
    const double d = -sineAxis.yMm;
    const double determinant = a * d - b * c;
    if (!std::isfinite(determinant) || std::abs(determinant) <= 1.0e-15)
        return std::nullopt;

    // A*A^Tの固有値・固有vectorから，scaleを持たないcanonical ellipseへ分解する。
    // unit circle側の直交変換はfull ellipseの集合を変えないのでSVGへ保持する必要はない。
    const double xx = a * a + b * b;
    const double xy = a * c + b * d;
    const double yy = c * c + d * d;
    const double trace = xx + yy;
    const double discriminant = std::hypot(xx - yy, 2.0 * xy);
    const double lambdaX = 0.5 * (trace + discriminant);
    if (!(lambdaX > 0.0) || !std::isfinite(lambdaX))
        return std::nullopt;
    // 小さい固有値はtrace-discriminantだと扁平な楕円で桁落ちするためdet(A)^2/lambdaXから得る。
    const double lambdaY = (determinant * determinant) / lambdaX;
    if (!(lambdaY > 0.0) || !std::isfinite(lambdaY))
        return std::nullopt;

    const double radiusX = std::sqrt(lambdaX);
    const double radiusY = std::sqrt(lambdaY);
    if (!std::isfinite(radiusX) || !std::isfinite(radiusY))
        return std::nullopt;

    double rotationRadians = 0.0;
    const double scale = std::max({xx, yy, 1.0});
    if (discriminant > scale * 1.0e-14)
        rotationRadians = 0.5 * std::atan2(2.0 * xy, xx - yy);

    return SvgEllipseGeometry{radiusX, radiusY, rotationRadians, determinant};
}

[[nodiscard]] GraphicsPointMm svgEllipsePoint(
    GraphicsPointMm center,
    GraphicsPointMm cosineAxis,
    GraphicsPointMm sineAxis,
    double angle,
    double heightMm) noexcept {
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    return GraphicsPointMm{
        center.xMm + cosine * cosineAxis.xMm + sine * sineAxis.xMm,
        heightMm - (center.yMm + cosine * cosineAxis.yMm + sine * sineAxis.yMm)};
}

[[nodiscard]] bool writeEllipse(
    std::ostringstream& out,
    const GraphicsEllipseNode& ellipse,
    double heightMm,
    std::optional<std::size_t> clipId,
    std::string_view indent = "  ") {
    if (!finite(ellipse.center) || !finite(ellipse.cosineAxis) || !finite(ellipse.sineAxis))
        return false;
    const auto geometry = canonicalSvgEllipse(ellipse.cosineAxis, ellipse.sineAxis);
    if (!geometry)
        return false;

    const double centerX = ellipse.center.xMm;
    const double centerY = heightMm - ellipse.center.yMm;
    const double degrees = geometry->rotationRadians * 180.0 / std::acos(-1.0);

    // scale transform + non-scaling-strokeはIllustrator等で展開時に線幅がscaleされることがある。
    // 半径へscaleを焼き込み，必要な回転だけを残してstroke-widthを通常のmm値として保持する。
    out << indent << "<ellipse cx=\"" << number(centerX)
        << "\" cy=\"" << number(centerY)
        << "\" rx=\"" << number(geometry->radiusX)
        << "\" ry=\"" << number(geometry->radiusY) << '"';
    if (std::abs(geometry->rotationRadians) > 1.0e-14)
        out << " transform=\"rotate(" << number(degrees) << ' '
            << number(centerX) << ' ' << number(centerY) << ")\"";
    if (ellipse.stroke) {
        out << " stroke=\"" << color(ellipse.stroke->color) << '"'
            << " stroke-width=\"" << number(ellipse.stroke->widthMm) << '"'
            << " stroke-linecap=\"" << lineCap(ellipse.stroke->lineCap) << '"'
            << " stroke-linejoin=\"" << lineJoin(ellipse.stroke->lineJoin) << '"';
        if (ellipse.stroke->color.alpha != 255)
            out << " stroke-opacity=\""
                << number(static_cast<double>(ellipse.stroke->color.alpha) / 255.0) << '"';
    }
    else
        out << " stroke=\"none\"";
    if (ellipse.fill) {
        out << " fill=\"" << color(ellipse.fill->color) << '"';
        if (ellipse.fill->color.alpha != 255)
            out << " fill-opacity=\""
                << number(static_cast<double>(ellipse.fill->color.alpha) / 255.0) << '"';
    }
    else
        out << " fill=\"none\"";
    if (clipId)
        out << " clip-path=\"url(#mmcal-clip-" << *clipId << ")\"";
    writeSemantic(out, ellipse.semantic);
    out << "/>\n";
    return true;
}

[[nodiscard]] bool writeEllipticArc(
    std::ostringstream& out,
    const GraphicsEllipticArcNode& arc,
    double heightMm,
    std::optional<std::size_t> clipId,
    std::string_view indent = "  ") {
    if (!finite(arc.center) || !finite(arc.cosineAxis) || !finite(arc.sineAxis)
        || !std::isfinite(arc.startRadians) || !std::isfinite(arc.sweepRadians)
        || arc.sweepRadians == 0.0)
        return false;
    const auto geometry = canonicalSvgEllipse(arc.cosineAxis, arc.sineAxis);
    if (!geometry)
        return false;

    const double halfPi = std::acos(-1.0) / 2.0;
    const std::size_t pieces = std::max<std::size_t>(1,
        static_cast<std::size_t>(std::ceil(std::abs(arc.sweepRadians) / halfPi)));
    const double delta = arc.sweepRadians / static_cast<double>(pieces);
    const double degrees = geometry->rotationRadians * 180.0 / std::acos(-1.0);
    const int sweepFlag = delta * geometry->determinant > 0.0 ? 1 : 0;

    const auto start = svgEllipsePoint(
        arc.center, arc.cosineAxis, arc.sineAxis, arc.startRadians, heightMm);
    out << indent << "<path d=\"M" << number(start.xMm) << ' ' << number(start.yMm);
    for (std::size_t i = 1; i <= pieces; ++i) {
        const auto point = svgEllipsePoint(
            arc.center,
            arc.cosineAxis,
            arc.sineAxis,
            arc.startRadians + delta * static_cast<double>(i),
            heightMm);
        out << " A" << number(geometry->radiusX) << ' ' << number(geometry->radiusY)
            << ' ' << number(degrees) << " 0 " << sweepFlag << ' '
            << number(point.xMm) << ' ' << number(point.yMm);
    }
    out << '"';
    if (arc.stroke) {
        out << " stroke=\"" << color(arc.stroke->color) << '\"'
            << " stroke-width=\"" << number(arc.stroke->widthMm) << '\"'
            << " stroke-linecap=\"" << lineCap(arc.stroke->lineCap) << '\"'
            << " stroke-linejoin=\"" << lineJoin(arc.stroke->lineJoin) << '\"';
        if (arc.stroke->color.alpha != 255)
            out << " stroke-opacity=\""
                << number(static_cast<double>(arc.stroke->color.alpha) / 255.0) << '\"';
    }
    else
        out << " stroke=\"none\"";
    if (arc.fill) {
        out << " fill=\"" << color(arc.fill->color) << '\"';
        if (arc.fill->color.alpha != 255)
            out << " fill-opacity=\""
                << number(static_cast<double>(arc.fill->color.alpha) / 255.0) << '\"';
    }
    else
        out << " fill=\"none\"";
    if (clipId)
        out << " clip-path=\"url(#mmcal-clip-" << *clipId << ")\"";
    writeSemantic(out, arc.semantic);
    out << "/>\n";
    return true;
}

[[nodiscard]] bool writeCircle(
    std::ostringstream& out,
    const GraphicsCircleNode& circle,
    double heightMm,
    std::optional<std::size_t> clipId,
    std::string_view indent = "  ") {
    if (!finite(circle.center) || !(circle.radiusMm > 0.0) || !std::isfinite(circle.radiusMm))
        return false;
    out << indent << "<circle cx=\"" << number(circle.center.xMm)
        << "\" cy=\"" << number(heightMm - circle.center.yMm)
        << "\" r=\"" << number(circle.radiusMm) << '"';
    if (circle.stroke) {
        out << " stroke=\"" << color(circle.stroke->color) << '"'
            << " stroke-width=\"" << number(circle.stroke->widthMm) << '"'
            << " stroke-linecap=\"" << lineCap(circle.stroke->lineCap) << '"'
            << " stroke-linejoin=\"" << lineJoin(circle.stroke->lineJoin) << '"';
        if (circle.stroke->color.alpha != 255)
            out << " stroke-opacity=\""
                << number(static_cast<double>(circle.stroke->color.alpha) / 255.0) << '"';
    }
    else
        out << " stroke=\"none\"";
    if (circle.fill) {
        out << " fill=\"" << color(circle.fill->color) << '"';
        if (circle.fill->color.alpha != 255)
            out << " fill-opacity=\""
                << number(static_cast<double>(circle.fill->color.alpha) / 255.0) << '"';
    }
    else
        out << " fill=\"none\"";
    if (clipId)
        out << " clip-path=\"url(#mmcal-clip-" << *clipId << ")\"";
    writeSemantic(out, circle.semantic);
    out << "/>\n";
    return true;
}

[[nodiscard]] bool writeText(
    std::ostringstream& out,
    const GraphicsTextNode& text,
    double heightMm,
    std::string_view indent = "  ") {
    if (!finite(text.origin) || !(text.fontSizeMm > 0.0) || !std::isfinite(text.fontSizeMm))
        return false;
    out << indent << "<text x=\"" << number(text.origin.xMm)
        << "\" y=\"" << number(heightMm - text.origin.yMm)
        << "\" font-size=\"" << number(text.fontSizeMm)
        << "\" font-family=\"" << escapeXml(text.fontFamily)
        << "\" text-anchor=\"" << textAnchor(text.anchor)
        << "\" fill=\"" << color(text.color) << '"';
    if (text.color.alpha != 255)
        out << " fill-opacity=\""
            << number(static_cast<double>(text.color.alpha) / 255.0) << '"';
    writeSemantic(out, text.semantic);
    out << '>' << escapeXml(text.text) << "</text>\n";
    return true;
}

[[nodiscard]] bool hasSemanticKind(
    const GraphicsNode& node,
    GraphicsSemanticKind kind) noexcept {
    return std::visit([kind](const auto& value) {
        return value.semantic && value.semantic->kind == kind;
    }, node);
}

[[nodiscard]] bool validClipRect(const GraphicsRectMm& rect) noexcept {
    return std::isfinite(rect.xMm) && std::isfinite(rect.yMm)
        && std::isfinite(rect.widthMm) && std::isfinite(rect.heightMm)
        && rect.widthMm > 0.0 && rect.heightMm > 0.0;
}

void writeClipDefinition(
    std::ostringstream& out,
    const GraphicsRectMm& rect,
    double heightMm,
    std::size_t clipId,
    std::string_view indent) {
    out << indent << "<defs><clipPath id=\"mmcal-clip-" << clipId
        << "\"><rect x=\"" << number(rect.xMm)
        << "\" y=\"" << number(heightMm - rect.yMm - rect.heightMm)
        << "\" width=\"" << number(rect.widthMm)
        << "\" height=\"" << number(rect.heightMm)
        << "\"/></clipPath></defs>\n";
}

[[nodiscard]] bool writeNode(
    std::ostringstream& out,
    const GraphicsNode& node,
    double heightMm,
    std::size_t& clipCounter,
    std::string_view indent = "  ") {
    return std::visit([&](const auto& value) {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, GraphicsPathNode>
            || std::is_same_v<T, GraphicsEllipseNode>
            || std::is_same_v<T, GraphicsEllipticArcNode>
            || std::is_same_v<T, GraphicsCircleNode>) {
            std::optional<std::size_t> clipId;
            if (value.clipRect) {
                if (!validClipRect(*value.clipRect))
                    return false;
                clipId = clipCounter++;
                writeClipDefinition(out, *value.clipRect, heightMm, *clipId, indent);
            }
            if constexpr (std::is_same_v<T, GraphicsPathNode>)
                return writePath(out, value, heightMm, clipId, indent);
            else if constexpr (std::is_same_v<T, GraphicsEllipseNode>)
                return writeEllipse(out, value, heightMm, clipId, indent);
            else if constexpr (std::is_same_v<T, GraphicsEllipticArcNode>)
                return writeEllipticArc(out, value, heightMm, clipId, indent);
            else
                return writeCircle(out, value, heightMm, clipId, indent);
        }
        else
            return writeText(out, value, heightMm, indent);
    }, node);
}

[[nodiscard]] std::optional<std::size_t> writeGroupedNodes(
    std::ostringstream& out,
    const GraphicsScene& scene,
    std::size_t startIndex,
    GraphicsSemanticKind semanticKind,
    std::string_view groupKind,
    std::size_t& clipCounter) {
    std::size_t i = startIndex;
    if (i >= scene.nodes.size() || !hasSemanticKind(scene.nodes[i], semanticKind))
        return startIndex;

    out << "  <g data-mmcal-kind=\"" << groupKind << "\">\n";
    while (i < scene.nodes.size() && hasSemanticKind(scene.nodes[i], semanticKind)) {
        if (!writeNode(out, scene.nodes[i], scene.extent.heightMm, clipCounter, "    "))
            return std::nullopt;
        ++i;
    }
    out << "  </g>\n";
    return i;
}

} // namespace

SvgRenderResult renderSvg(const GraphicsScene& scene) {
    SvgRenderResult result;
    if (!(scene.extent.widthMm > 0.0) || !(scene.extent.heightMm > 0.0)
        || !std::isfinite(scene.extent.widthMm) || !std::isfinite(scene.extent.heightMm))
        return result;

    std::ostringstream out;
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
        << number(scene.extent.widthMm) << "mm\" height=\""
        << number(scene.extent.heightMm) << "mm\" viewBox=\"0 0 "
        << number(scene.extent.widthMm) << ' ' << number(scene.extent.heightMm) << "\">\n";

    std::size_t clipCounter = 0;
    for (std::size_t i = 0; i < scene.nodes.size();) {
        const auto xTickGroup = writeGroupedNodes(
            out, scene, i, GraphicsSemanticKind::XTick, "x-ticks", clipCounter);
        if (!xTickGroup) {
            result.status = SvgRenderStatus::NonFiniteGeometry;
            return result;
        }
        if (*xTickGroup != i) {
            i = *xTickGroup;
            continue;
        }

        const auto yTickGroup = writeGroupedNodes(
            out, scene, i, GraphicsSemanticKind::YTick, "y-ticks", clipCounter);
        if (!yTickGroup) {
            result.status = SvgRenderStatus::NonFiniteGeometry;
            return result;
        }
        if (*yTickGroup != i) {
            i = *yTickGroup;
            continue;
        }

        const auto xLabelGroup = writeGroupedNodes(
            out, scene, i, GraphicsSemanticKind::XTickLabel, "x-tick-labels", clipCounter);
        if (!xLabelGroup) {
            result.status = SvgRenderStatus::NonFiniteGeometry;
            return result;
        }
        if (*xLabelGroup != i) {
            i = *xLabelGroup;
            continue;
        }

        const auto yLabelGroup = writeGroupedNodes(
            out, scene, i, GraphicsSemanticKind::YTickLabel, "y-tick-labels", clipCounter);
        if (!yLabelGroup) {
            result.status = SvgRenderStatus::NonFiniteGeometry;
            return result;
        }
        if (*yLabelGroup != i) {
            i = *yLabelGroup;
            continue;
        }

        const auto& node = scene.nodes[i];
        if (!writeNode(out, node, scene.extent.heightMm, clipCounter)) {
            result.status = SvgRenderStatus::NonFiniteGeometry;
            return result;
        }
        ++i;
    }
    out << "</svg>\n";

    result.status = SvgRenderStatus::Success;
    result.svg = out.str();
    return result;
}

} // namespace mmcal::graphics
