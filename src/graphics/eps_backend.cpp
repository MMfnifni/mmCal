#include "eps_backend.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>

namespace mmcal::graphics {
namespace {

constexpr double pointsPerMillimeter = 72.0 / 25.4;

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

[[nodiscard]] bool opaque(const GraphicsColor& color) noexcept {
    return color.alpha == 255;
}

[[nodiscard]] std::string rgb(const GraphicsColor& color) {
    std::ostringstream out;
    out << std::setprecision(9)
        << static_cast<double>(color.red) / 255.0 << ' '
        << static_cast<double>(color.green) / 255.0 << ' '
        << static_cast<double>(color.blue) / 255.0;
    return out.str();
}

[[nodiscard]] int lineCap(GraphicsLineCap value) noexcept {
    switch (value) {
    case GraphicsLineCap::Butt: return 0;
    case GraphicsLineCap::Round: return 1;
    case GraphicsLineCap::Square: return 2;
    }
    return 0;
}

[[nodiscard]] int lineJoin(GraphicsLineJoin value) noexcept {
    switch (value) {
    case GraphicsLineJoin::Miter: return 0;
    case GraphicsLineJoin::Round: return 1;
    case GraphicsLineJoin::Bevel: return 2;
    }
    return 0;
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

void writeSemantic(
    std::ostringstream& out,
    const std::optional<GraphicsSemanticRef>& semantic) {
    if (semantic)
        out << "% mmCal-semantic: " << semanticKind(semantic->kind)
            << ' ' << semantic->id << '\n';
}

[[nodiscard]] bool validStroke(const GraphicsStrokeStyle& stroke) noexcept {
    return stroke.widthMm > 0.0 && std::isfinite(stroke.widthMm) && opaque(stroke.color);
}

[[nodiscard]] bool validFill(const GraphicsFillStyle& fill) noexcept {
    return opaque(fill.color);
}

void writeStrokeStyle(std::ostringstream& out, const GraphicsStrokeStyle& stroke) {
    out << number(stroke.widthMm) << " setlinewidth\n"
        << lineCap(stroke.lineCap) << " setlinecap\n"
        << lineJoin(stroke.lineJoin) << " setlinejoin\n"
        << rgb(stroke.color) << " setrgbcolor\n";
}

void writePaint(
    std::ostringstream& out,
    const std::optional<GraphicsStrokeStyle>& stroke,
    const std::optional<GraphicsFillStyle>& fill) {
    if (fill && stroke) {
        out << "gsave\n" << rgb(fill->color) << " setrgbcolor\nfill\ngrestore\n";
        writeStrokeStyle(out, *stroke);
        out << "stroke\n";
    }
    else if (fill)
        out << rgb(fill->color) << " setrgbcolor\nfill\n";
    else if (stroke) {
        writeStrokeStyle(out, *stroke);
        out << "stroke\n";
    }
    else
        out << "newpath\n";
}

[[nodiscard]] GraphicsPointMm interpolate(
    GraphicsPointMm lhs,
    GraphicsPointMm rhs,
    double t) noexcept {
    return {
        lhs.xMm + (rhs.xMm - lhs.xMm) * t,
        lhs.yMm + (rhs.yMm - lhs.yMm) * t};
}

[[nodiscard]] bool validClipRect(const GraphicsRectMm& rect) noexcept {
    return std::isfinite(rect.xMm) && std::isfinite(rect.yMm)
        && std::isfinite(rect.widthMm) && std::isfinite(rect.heightMm)
        && rect.widthMm > 0.0 && rect.heightMm > 0.0;
}

void beginClip(std::ostringstream& out, const GraphicsRectMm& rect) {
    out << "gsave\nnewpath\n"
        << number(rect.xMm) << ' ' << number(rect.yMm) << " moveto\n"
        << number(rect.xMm + rect.widthMm) << ' ' << number(rect.yMm) << " lineto\n"
        << number(rect.xMm + rect.widthMm) << ' ' << number(rect.yMm + rect.heightMm) << " lineto\n"
        << number(rect.xMm) << ' ' << number(rect.yMm + rect.heightMm) << " lineto\n"
        << "closepath\nclip\nnewpath\n";
}

[[nodiscard]] bool writePath(std::ostringstream& out, const GraphicsPathNode& path) {
    if (path.stroke && !validStroke(*path.stroke))
        return false;
    if (path.fill && !validFill(*path.fill))
        return false;
    if (path.clipRect) {
        if (!validClipRect(*path.clipRect))
            return false;
        beginClip(out, *path.clipRect);
    }

    writeSemantic(out, path.semantic);
    out << "newpath\n";
    std::optional<GraphicsPointMm> current;
    std::optional<GraphicsPointMm> subpathStart;

    for (const auto& command : path.commands) {
        bool valid = true;
        std::visit([&](const auto& value) {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, GraphicsMoveTo>) {
                valid = finite(value.point);
                if (valid) {
                    out << number(value.point.xMm) << ' ' << number(value.point.yMm) << " moveto\n";
                    current = value.point;
                    subpathStart = value.point;
                }
            }
            else if constexpr (std::is_same_v<T, GraphicsLineTo>) {
                valid = current.has_value() && finite(value.point);
                if (valid) {
                    out << number(value.point.xMm) << ' ' << number(value.point.yMm) << " lineto\n";
                    current = value.point;
                }
            }
            else if constexpr (std::is_same_v<T, GraphicsQuadraticTo>) {
                valid = current.has_value() && finite(value.control) && finite(value.point);
                if (valid) {
                    // PostScriptにはquadratic operatorが無いので，形状を変えずcubicへdegree elevationする。
                    const auto c1 = interpolate(*current, value.control, 2.0 / 3.0);
                    const auto c2 = interpolate(value.point, value.control, 2.0 / 3.0);
                    out << number(c1.xMm) << ' ' << number(c1.yMm) << ' '
                        << number(c2.xMm) << ' ' << number(c2.yMm) << ' '
                        << number(value.point.xMm) << ' ' << number(value.point.yMm)
                        << " curveto\n";
                    current = value.point;
                }
            }
            else if constexpr (std::is_same_v<T, GraphicsCubicTo>) {
                valid = current.has_value()
                    && finite(value.control1) && finite(value.control2) && finite(value.point);
                if (valid) {
                    out << number(value.control1.xMm) << ' ' << number(value.control1.yMm) << ' '
                        << number(value.control2.xMm) << ' ' << number(value.control2.yMm) << ' '
                        << number(value.point.xMm) << ' ' << number(value.point.yMm)
                        << " curveto\n";
                    current = value.point;
                }
            }
            else if constexpr (std::is_same_v<T, GraphicsClosePath>) {
                valid = subpathStart.has_value();
                if (valid) {
                    out << "closepath\n";
                    current = subpathStart;
                }
            }
        }, command);
        if (!valid)
            return false;
    }

    writePaint(out, path.stroke, path.fill);
    if (path.clipRect)
        out << "grestore\n";
    return true;
}

[[nodiscard]] bool writeEllipse(std::ostringstream& out, const GraphicsEllipseNode& ellipse) {
    if (!finite(ellipse.center) || !finite(ellipse.cosineAxis) || !finite(ellipse.sineAxis))
        return false;
    const double determinant = ellipse.cosineAxis.xMm * ellipse.sineAxis.yMm
        - ellipse.cosineAxis.yMm * ellipse.sineAxis.xMm;
    if (!std::isfinite(determinant) || std::abs(determinant) <= 1.0e-15)
        return false;
    if (ellipse.stroke && !validStroke(*ellipse.stroke))
        return false;
    if (ellipse.fill && !validFill(*ellipse.fill))
        return false;
    if (ellipse.clipRect) {
        if (!validClipRect(*ellipse.clipRect))
            return false;
        beginClip(out, *ellipse.clipRect);
    }

    writeSemantic(out, ellipse.semantic);
    // PostScriptのarcへunit circleを渡し，CTMだけ楕円のaffine mapへ一時変更する。
    // path構築後にCTMを戻すため，stroke幅はGraphicsSceneのmm指定を維持する。
    out << "matrix currentmatrix\n["
        << number(ellipse.cosineAxis.xMm) << ' ' << number(ellipse.cosineAxis.yMm) << ' '
        << number(ellipse.sineAxis.xMm) << ' ' << number(ellipse.sineAxis.yMm) << ' '
        << number(ellipse.center.xMm) << ' ' << number(ellipse.center.yMm)
        << "] concat\nnewpath\n0 0 1 0 360 arc\nclosepath\nsetmatrix\n";
    writePaint(out, ellipse.stroke, ellipse.fill);
    if (ellipse.clipRect)
        out << "grestore\n";
    return true;
}

[[nodiscard]] bool writeEllipticArc(std::ostringstream& out, const GraphicsEllipticArcNode& arc) {
    if (!finite(arc.center) || !finite(arc.cosineAxis) || !finite(arc.sineAxis)
        || !std::isfinite(arc.startRadians) || !std::isfinite(arc.sweepRadians)
        || arc.sweepRadians == 0.0)
        return false;
    const double determinant = arc.cosineAxis.xMm * arc.sineAxis.yMm
        - arc.cosineAxis.yMm * arc.sineAxis.xMm;
    if (!std::isfinite(determinant) || std::abs(determinant) <= 1.0e-15)
        return false;
    if (arc.stroke && !validStroke(*arc.stroke))
        return false;
    if (arc.fill && !validFill(*arc.fill))
        return false;
    if (arc.clipRect) {
        if (!validClipRect(*arc.clipRect))
            return false;
        beginClip(out, *arc.clipRect);
    }

    const double halfPi = std::acos(-1.0) / 2.0;
    const std::size_t pieces = std::max<std::size_t>(1,
        static_cast<std::size_t>(std::ceil(std::abs(arc.sweepRadians) / halfPi)));
    const double delta = arc.sweepRadians / static_cast<double>(pieces);

    writeSemantic(out, arc.semantic);
    out << "matrix currentmatrix\n["
        << number(arc.cosineAxis.xMm) << ' ' << number(arc.cosineAxis.yMm) << ' '
        << number(arc.sineAxis.xMm) << ' ' << number(arc.sineAxis.yMm) << ' '
        << number(arc.center.xMm) << ' ' << number(arc.center.yMm)
        << "] concat\nnewpath\n";
    double startDegrees = arc.startRadians * 180.0 / std::acos(-1.0);
    const double deltaDegrees = delta * 180.0 / std::acos(-1.0);
    for (std::size_t i = 0; i < pieces; ++i) {
        const double endDegrees = startDegrees + deltaDegrees;
        out << "0 0 1 " << number(startDegrees) << ' ' << number(endDegrees)
            << (delta > 0.0 ? " arc\n" : " arcn\n");
        startDegrees = endDegrees;
    }
    out << "setmatrix\n";
    writePaint(out, arc.stroke, arc.fill);
    if (arc.clipRect)
        out << "grestore\n";
    return true;
}

[[nodiscard]] bool writeCircle(std::ostringstream& out, const GraphicsCircleNode& circle) {
    if (!finite(circle.center) || !(circle.radiusMm > 0.0) || !std::isfinite(circle.radiusMm))
        return false;
    if (circle.stroke && !validStroke(*circle.stroke))
        return false;
    if (circle.fill && !validFill(*circle.fill))
        return false;
    if (circle.clipRect) {
        if (!validClipRect(*circle.clipRect))
            return false;
        beginClip(out, *circle.clipRect);
    }

    writeSemantic(out, circle.semantic);
    out << "newpath\n"
        << number(circle.center.xMm) << ' ' << number(circle.center.yMm) << ' '
        << number(circle.radiusMm) << " 0 360 arc\nclosepath\n";
    writePaint(out, circle.stroke, circle.fill);
    if (circle.clipRect)
        out << "grestore\n";
    return true;
}

[[nodiscard]] std::optional<std::string> postScriptString(std::string_view text) {
    std::string result;
    result.reserve(text.size() + 8);
    for (const unsigned char ch : text) {
        if (ch < 0x20 || ch > 0x7e)
            return std::nullopt;
        if (ch == '(' || ch == ')' || ch == '\\')
            result.push_back('\\');
        result.push_back(static_cast<char>(ch));
    }
    return result;
}

[[nodiscard]] const char* postScriptFont(std::string_view family) noexcept {
    if (family == "serif")
        return "Times-Roman";
    if (family == "monospace")
        return "Courier";
    return "Helvetica";
}

[[nodiscard]] bool writeText(std::ostringstream& out, const GraphicsTextNode& text) {
    if (!finite(text.origin) || !(text.fontSizeMm > 0.0) || !std::isfinite(text.fontSizeMm)
        || !opaque(text.color))
        return false;
    const auto escaped = postScriptString(text.text);
    if (!escaped)
        return false;

    writeSemantic(out, text.semantic);
    out << '/' << postScriptFont(text.fontFamily) << " findfont "
        << number(text.fontSizeMm) << " scalefont setfont\n"
        << rgb(text.color) << " setrgbcolor\n"
        << number(text.origin.xMm) << ' ' << number(text.origin.yMm) << " moveto\n";
    if (text.anchor == GraphicsTextAnchor::Middle)
        out << '(' << *escaped << ") dup stringwidth pop -0.5 mul 0 rmoveto show\n";
    else if (text.anchor == GraphicsTextAnchor::End)
        out << '(' << *escaped << ") dup stringwidth pop neg 0 rmoveto show\n";
    else
        out << '(' << *escaped << ") show\n";
    return true;
}

[[nodiscard]] bool hasUnsupportedTransparency(const GraphicsNode& node) noexcept {
    return std::visit([](const auto& value) {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, GraphicsPathNode>
            || std::is_same_v<T, GraphicsEllipseNode>
            || std::is_same_v<T, GraphicsEllipticArcNode>
            || std::is_same_v<T, GraphicsCircleNode>)
            return (value.stroke && !opaque(value.stroke->color))
                || (value.fill && !opaque(value.fill->color));
        else
            return !opaque(value.color);
    }, node);
}

[[nodiscard]] bool hasUnsupportedText(const GraphicsNode& node) {
    const auto* text = std::get_if<GraphicsTextNode>(&node);
    return text && !postScriptString(text->text).has_value();
}

} // namespace

EpsRenderResult renderEps(const GraphicsScene& scene) {
    EpsRenderResult result;
    if (!(scene.extent.widthMm > 0.0) || !(scene.extent.heightMm > 0.0)
        || !std::isfinite(scene.extent.widthMm) || !std::isfinite(scene.extent.heightMm))
        return result;

    for (const auto& node : scene.nodes) {
        if (hasUnsupportedTransparency(node)) {
            result.status = EpsRenderStatus::UnsupportedTransparency;
            return result;
        }
        if (hasUnsupportedText(node)) {
            result.status = EpsRenderStatus::UnsupportedText;
            return result;
        }
    }

    const double widthPt = scene.extent.widthMm * pointsPerMillimeter;
    const double heightPt = scene.extent.heightMm * pointsPerMillimeter;
    std::ostringstream out;
    out << "%!PS-Adobe-3.0 EPSF-3.0\n"
        << "%%BoundingBox: 0 0 " << static_cast<long long>(std::ceil(widthPt))
        << ' ' << static_cast<long long>(std::ceil(heightPt)) << "\n"
        << "%%HiResBoundingBox: 0 0 " << number(widthPt) << ' ' << number(heightPt) << "\n"
        << "%%Creator: mmCal\n"
        << "%%LanguageLevel: 2\n"
        << "%%Pages: 1\n"
        << "%%EndComments\n"
        << "gsave\n"
        // user spaceをmmへする。BoundingBoxだけはDSC規約どおりPostScript pointで記述する。
        << number(pointsPerMillimeter) << ' ' << number(pointsPerMillimeter) << " scale\n"
        // EPSは埋め込み用途なので，canvas外へ伸びる漸近線等を明示的にviewportでclipする。
        << "newpath\n0 0 moveto\n"
        << number(scene.extent.widthMm) << " 0 lineto\n"
        << number(scene.extent.widthMm) << ' ' << number(scene.extent.heightMm) << " lineto\n"
        << "0 " << number(scene.extent.heightMm) << " lineto\nclosepath\nclip\nnewpath\n";

    for (const auto& node : scene.nodes) {
        const bool written = std::visit([&](const auto& value) {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, GraphicsPathNode>)
                return writePath(out, value);
            else if constexpr (std::is_same_v<T, GraphicsEllipseNode>)
                return writeEllipse(out, value);
            else if constexpr (std::is_same_v<T, GraphicsEllipticArcNode>)
                return writeEllipticArc(out, value);
            else if constexpr (std::is_same_v<T, GraphicsCircleNode>)
                return writeCircle(out, value);
            else
                return writeText(out, value);
        }, node);
        if (!written) {
            result.status = EpsRenderStatus::NonFiniteGeometry;
            return result;
        }
    }

    out << "grestore\n%%EOF\n";
    result.status = EpsRenderStatus::Success;
    result.eps = out.str();
    return result;
}

} // namespace mmcal::graphics
