#include "pdf_backend.hpp"

#include "text_metrics.hpp"
#include "../../version.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#else
#include <pwd.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace mmcal::graphics {
namespace {

constexpr double pointsPerMillimeter = 72.0 / 25.4;
constexpr double circleKappa = 0.5522847498307936;
constexpr std::string_view xmpNamespace = "urn:mmcal:metadata:1.0";

struct PdfTimestamp final {
    std::string pdf;
    std::string xmp;
};

struct PdfId final {
    std::array<std::uint8_t, 16> bytes{};
};

[[nodiscard]] bool finite(GraphicsPointMm point) noexcept {
    return std::isfinite(point.xMm) && std::isfinite(point.yMm);
}

[[nodiscard]] bool opaque(const GraphicsColor& color) noexcept {
    return color.alpha == 255;
}

[[nodiscard]] std::string number(double value) {
    if (value == 0.0)
        value = 0.0;
    std::ostringstream out;
    out << std::setprecision(15) << value;
    return out.str();
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

[[nodiscard]] std::optional<std::string> pdfAsciiString(std::string_view text) {
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

[[nodiscard]] std::string xmlEscape(std::string_view text) {
    std::string result;
    result.reserve(text.size() + 16);
    for (const char ch : text) {
        switch (ch) {
        case '&': result += "&amp;"; break;
        case '<': result += "&lt;"; break;
        case '>': result += "&gt;"; break;
        case '\"': result += "&quot;"; break;
        case '\'': result += "&apos;"; break;
        default: result.push_back(ch); break;
        }
    }
    return result;
}

[[nodiscard]] const char* pdfFontResource(std::string_view family) noexcept {
    if (family == "serif")
        return "F2";
    if (family == "monospace")
        return "F3";
    return "F1";
}

[[nodiscard]] bool validStroke(const GraphicsStrokeStyle& stroke) noexcept {
    return stroke.widthMm > 0.0 && std::isfinite(stroke.widthMm) && opaque(stroke.color);
}

[[nodiscard]] bool validFill(const GraphicsFillStyle& fill) noexcept {
    return opaque(fill.color);
}

void writeStrokeStyle(std::ostringstream& out, const GraphicsStrokeStyle& stroke) {
    out << number(stroke.widthMm) << " w\n"
        << lineCap(stroke.lineCap) << " J\n"
        << lineJoin(stroke.lineJoin) << " j\n"
        << rgb(stroke.color) << " RG\n";
}

void writePaint(
    std::ostringstream& out,
    const std::optional<GraphicsStrokeStyle>& stroke,
    const std::optional<GraphicsFillStyle>& fill) {
    if (fill)
        out << rgb(fill->color) << " rg\n";
    if (stroke)
        writeStrokeStyle(out, *stroke);
    if (fill && stroke)
        out << "B\n";
    else if (fill)
        out << "f\n";
    else if (stroke)
        out << "S\n";
    else
        out << "n\n";
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
    out << "q\n"
        << number(rect.xMm) << ' ' << number(rect.yMm) << ' '
        << number(rect.widthMm) << ' ' << number(rect.heightMm) << " re W n\n";
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
    std::optional<GraphicsPointMm> current;
    std::optional<GraphicsPointMm> subpathStart;
    for (const auto& command : path.commands) {
        bool valid = true;
        std::visit([&](const auto& value) {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, GraphicsMoveTo>) {
                valid = finite(value.point);
                if (valid) {
                    out << number(value.point.xMm) << ' ' << number(value.point.yMm) << " m\n";
                    current = value.point;
                    subpathStart = value.point;
                }
            }
            else if constexpr (std::is_same_v<T, GraphicsLineTo>) {
                valid = current.has_value() && finite(value.point);
                if (valid) {
                    out << number(value.point.xMm) << ' ' << number(value.point.yMm) << " l\n";
                    current = value.point;
                }
            }
            else if constexpr (std::is_same_v<T, GraphicsQuadraticTo>) {
                valid = current.has_value() && finite(value.control) && finite(value.point);
                if (valid) {
                    // PDFにもquadratic operatorは無いので，exact degree elevationしてcubicにする。
                    const auto c1 = interpolate(*current, value.control, 2.0 / 3.0);
                    const auto c2 = interpolate(value.point, value.control, 2.0 / 3.0);
                    out << number(c1.xMm) << ' ' << number(c1.yMm) << ' '
                        << number(c2.xMm) << ' ' << number(c2.yMm) << ' '
                        << number(value.point.xMm) << ' ' << number(value.point.yMm)
                        << " c\n";
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
                        << " c\n";
                    current = value.point;
                }
            }
            else if constexpr (std::is_same_v<T, GraphicsClosePath>) {
                valid = subpathStart.has_value();
                if (valid) {
                    out << "h\n";
                    current = subpathStart;
                }
            }
        }, command);
        if (!valid)
            return false;
    }

    writePaint(out, path.stroke, path.fill);
    if (path.clipRect)
        out << "Q\n";
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

    const auto point = [&](double ca, double sa) {
        return GraphicsPointMm{
            ellipse.center.xMm + ca * ellipse.cosineAxis.xMm + sa * ellipse.sineAxis.xMm,
            ellipse.center.yMm + ca * ellipse.cosineAxis.yMm + sa * ellipse.sineAxis.yMm};
    };
    const auto offset = [&](GraphicsPointMm p, double ca, double sa) {
        return GraphicsPointMm{
            p.xMm + ca * ellipse.cosineAxis.xMm + sa * ellipse.sineAxis.xMm,
            p.yMm + ca * ellipse.cosineAxis.yMm + sa * ellipse.sineAxis.yMm};
    };
    const auto p0 = point(1.0, 0.0);
    const auto p1 = point(0.0, 1.0);
    const auto p2 = point(-1.0, 0.0);
    const auto p3 = point(0.0, -1.0);
    const double k = circleKappa;

    writeSemantic(out, ellipse.semantic);
    out << number(p0.xMm) << ' ' << number(p0.yMm) << " m\n";
    auto cubic = [&](GraphicsPointMm c1, GraphicsPointMm c2, GraphicsPointMm p) {
        out << number(c1.xMm) << ' ' << number(c1.yMm) << ' '
            << number(c2.xMm) << ' ' << number(c2.yMm) << ' '
            << number(p.xMm) << ' ' << number(p.yMm) << " c\n";
    };
    cubic(offset(p0, 0.0, k), offset(p1, k, 0.0), p1);
    cubic(offset(p1, -k, 0.0), offset(p2, 0.0, k), p2);
    cubic(offset(p2, 0.0, -k), offset(p3, -k, 0.0), p3);
    cubic(offset(p3, k, 0.0), offset(p0, 0.0, -k), p0);
    out << "h\n";
    writePaint(out, ellipse.stroke, ellipse.fill);
    if (ellipse.clipRect)
        out << "Q\n";
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
    const auto point = [&](double angle) {
        return GraphicsPointMm{
            arc.center.xMm + std::cos(angle) * arc.cosineAxis.xMm
                + std::sin(angle) * arc.sineAxis.xMm,
            arc.center.yMm + std::cos(angle) * arc.cosineAxis.yMm
                + std::sin(angle) * arc.sineAxis.yMm};
    };
    const auto tangent = [&](double angle) {
        return GraphicsPointMm{
            -std::sin(angle) * arc.cosineAxis.xMm
                + std::cos(angle) * arc.sineAxis.xMm,
            -std::sin(angle) * arc.cosineAxis.yMm
                + std::cos(angle) * arc.sineAxis.yMm};
    };

    writeSemantic(out, arc.semantic);
    GraphicsPointMm p0 = point(arc.startRadians);
    out << number(p0.xMm) << ' ' << number(p0.yMm) << " m\n";
    double theta0 = arc.startRadians;
    for (std::size_t i = 0; i < pieces; ++i) {
        const double theta1 = theta0 + delta;
        const GraphicsPointMm p1 = point(theta1);
        const GraphicsPointMm d0 = tangent(theta0);
        const GraphicsPointMm d1 = tangent(theta1);
        const double k = (4.0 / 3.0) * std::tan(delta / 4.0);
        const GraphicsPointMm c1{p0.xMm + k * d0.xMm, p0.yMm + k * d0.yMm};
        const GraphicsPointMm c2{p1.xMm - k * d1.xMm, p1.yMm - k * d1.yMm};
        out << number(c1.xMm) << ' ' << number(c1.yMm) << ' '
            << number(c2.xMm) << ' ' << number(c2.yMm) << ' '
            << number(p1.xMm) << ' ' << number(p1.yMm) << " c\n";
        theta0 = theta1;
        p0 = p1;
    }
    writePaint(out, arc.stroke, arc.fill);
    if (arc.clipRect)
        out << "Q\n";
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
    const double x = circle.center.xMm;
    const double y = circle.center.yMm;
    const double r = circle.radiusMm;
    const double k = r * circleKappa;
    out << number(x + r) << ' ' << number(y) << " m\n"
        << number(x + r) << ' ' << number(y + k) << ' '
        << number(x + k) << ' ' << number(y + r) << ' '
        << number(x) << ' ' << number(y + r) << " c\n"
        << number(x - k) << ' ' << number(y + r) << ' '
        << number(x - r) << ' ' << number(y + k) << ' '
        << number(x - r) << ' ' << number(y) << " c\n"
        << number(x - r) << ' ' << number(y - k) << ' '
        << number(x - k) << ' ' << number(y - r) << ' '
        << number(x) << ' ' << number(y - r) << " c\n"
        << number(x + k) << ' ' << number(y - r) << ' '
        << number(x + r) << ' ' << number(y - k) << ' '
        << number(x + r) << ' ' << number(y) << " c\nh\n";
    writePaint(out, circle.stroke, circle.fill);
    if (circle.clipRect)
        out << "Q\n";
    return true;
}

[[nodiscard]] bool writeText(std::ostringstream& out, const GraphicsTextNode& text) {
    if (!finite(text.origin) || !(text.fontSizeMm > 0.0) || !std::isfinite(text.fontSizeMm)
        || !opaque(text.color))
        return false;
    const auto escaped = pdfAsciiString(text.text);
    const auto metrics = measureTextApproximate(text.text, text.fontSizeMm);
    if (!escaped || !metrics)
        return false;

    double x = text.origin.xMm;
    if (text.anchor == GraphicsTextAnchor::Middle)
        x -= metrics->advanceWidthMm / 2.0;
    else if (text.anchor == GraphicsTextAnchor::End)
        x -= metrics->advanceWidthMm;

    writeSemantic(out, text.semantic);
    out << "BT\n/" << pdfFontResource(text.fontFamily) << ' '
        << number(text.fontSizeMm) << " Tf\n"
        << rgb(text.color) << " rg\n"
        << "1 0 0 1 " << number(x) << ' ' << number(text.origin.yMm) << " Tm\n"
        << '(' << *escaped << ") Tj\nET\n";
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
    return text && !pdfAsciiString(text->text).has_value();
}

[[nodiscard]] std::tm utcTm(std::time_t time) {
    std::tm result{};
#if defined(_WIN32)
    gmtime_s(&result, &time);
#else
    gmtime_r(&time, &result);
#endif
    return result;
}

[[nodiscard]] PdfTimestamp makeTimestamp(bool deterministic) {
    const std::time_t time = deterministic
        ? std::time_t{946684800} // 2000-01-01T00:00:00Z
        : std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    const std::tm tm = utcTm(time);

    std::ostringstream pdf;
    pdf << "D:" << std::put_time(&tm, "%Y%m%d%H%M%SZ");
    std::ostringstream xmp;
    xmp << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return PdfTimestamp{pdf.str(), xmp.str()};
}

[[nodiscard]] PdfId makeId(bool deterministic, std::uint8_t discriminator) {
    PdfId result;
    if (deterministic) {
        constexpr std::array<std::uint8_t, 16> base{
            0x6d, 0x6d, 0x43, 0x61, 0x6c, 0x50, 0x44, 0x46,
            0x80, 0x00, 0x14, 0x00, 0x00, 0x00, 0x00, 0x00};
        result.bytes = base;
        result.bytes.back() = discriminator;
        return result;
    }

    std::random_device random;
    for (auto& byte : result.bytes)
        byte = static_cast<std::uint8_t>(random());
    // UUID v4/variant bits。PDF trailer自体はUUIDを要求しないがXMPとの共有IDに使う。
    result.bytes[6] = static_cast<std::uint8_t>((result.bytes[6] & 0x0fU) | 0x40U);
    result.bytes[8] = static_cast<std::uint8_t>((result.bytes[8] & 0x3fU) | 0x80U);
    return result;
}

[[nodiscard]] std::string idHex(const PdfId& id) {
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (const auto byte : id.bytes)
        out << std::setw(2) << static_cast<unsigned>(byte);
    return out.str();
}

[[nodiscard]] std::string idUuid(const PdfId& id) {
    const std::string hex = idHex(id);
    return "uuid:" + hex.substr(0, 8) + '-' + hex.substr(8, 4) + '-'
        + hex.substr(12, 4) + '-' + hex.substr(16, 4) + '-' + hex.substr(20, 12);
}

#if defined(_WIN32)
[[nodiscard]] std::string wideToUtf8(std::wstring_view value) {
    if (value.empty())
        return {};
    const int size = WideCharToMultiByte(
        CP_UTF8, 0, value.data(), static_cast<int>(value.size()), nullptr, 0, nullptr, nullptr);
    if (size <= 0)
        return {};
    std::string result(static_cast<std::size_t>(size), '\0');
    WideCharToMultiByte(
        CP_UTF8, 0, value.data(), static_cast<int>(value.size()), result.data(), size, nullptr, nullptr);
    return result;
}
#endif

[[nodiscard]] std::optional<std::string> currentUserName() {
#if defined(_WIN32)
    const DWORD size = GetEnvironmentVariableW(L"USERNAME", nullptr, 0);
    if (size == 0)
        return std::nullopt;
    std::wstring value(size, L'\0');
    const DWORD written = GetEnvironmentVariableW(L"USERNAME", value.data(), size);
    if (written == 0 || written >= size)
        return std::nullopt;
    value.resize(written);
    std::string utf8 = wideToUtf8(value);
    return utf8.empty() ? std::nullopt : std::optional<std::string>{std::move(utf8)};
#else
    const uid_t uid = geteuid();
    long bufferSize = sysconf(_SC_GETPW_R_SIZE_MAX);
    if (bufferSize < 1024)
        bufferSize = 16384;
    std::vector<char> buffer(static_cast<std::size_t>(bufferSize));
    passwd record{};
    passwd* found = nullptr;
    if (getpwuid_r(uid, &record, buffer.data(), buffer.size(), &found) != 0
        || !found || !record.pw_name || *record.pw_name == '\0')
        return std::nullopt;
    return std::string{record.pw_name};
#endif
}

[[nodiscard]] std::string xmpPacket(
    const PdfTimestamp& timestamp,
    const PdfId& documentId,
    const PdfId& instanceId,
    const std::optional<std::string>& userName) {
    const std::string version = MMCAL_VERSION_STRING;
    const std::string title = "mmCal " + version + " Plot";
    const std::string creatorTool = "mmCal " + version;
    const std::string producer = "mmCal " + version + " PDF Plotter";

    std::ostringstream out;
    out << "<?xpacket begin=\"\" id=\"W5M0MpCehiHzreSzNTczkc9d\"?>\n"
        << "<x:xmpmeta xmlns:x=\"adobe:ns:meta/\">\n"
        << " <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\">\n"
        << "  <rdf:Description rdf:about=\"\"\n"
        << "   xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n"
        << "   xmlns:pdf=\"http://ns.adobe.com/pdf/1.3/\"\n"
        << "   xmlns:xmp=\"http://ns.adobe.com/xap/1.0/\"\n"
        << "   xmlns:xmpMM=\"http://ns.adobe.com/xap/1.0/mm/\"\n"
        << "   xmlns:mmcal=\"" << xmpNamespace << "\">\n"
        << "   <dc:title><rdf:Alt><rdf:li xml:lang=\"x-default\">"
        << xmlEscape(title) << "</rdf:li></rdf:Alt></dc:title>\n"
        << "   <dc:language><rdf:Bag><rdf:li>ja-JP</rdf:li></rdf:Bag></dc:language>\n"
        << "   <pdf:Producer>" << xmlEscape(producer) << "</pdf:Producer>\n"
        << "   <xmp:CreatorTool>" << xmlEscape(creatorTool) << "</xmp:CreatorTool>\n"
        << "   <xmp:CreateDate>" << timestamp.xmp << "</xmp:CreateDate>\n"
        << "   <xmp:ModifyDate>" << timestamp.xmp << "</xmp:ModifyDate>\n"
        << "   <xmp:MetadataDate>" << timestamp.xmp << "</xmp:MetadataDate>\n"
        << "   <xmpMM:DocumentID>" << idUuid(documentId) << "</xmpMM:DocumentID>\n"
        << "   <xmpMM:InstanceID>" << idUuid(instanceId) << "</xmpMM:InstanceID>\n"
        << "   <mmcal:Version>" << xmlEscape(version) << "</mmcal:Version>\n";
    if (userName)
        out << "   <mmcal:UserName>" << xmlEscape(*userName) << "</mmcal:UserName>\n";
    out << "  </rdf:Description>\n"
        << " </rdf:RDF>\n"
        << "</x:xmpmeta>\n"
        << "<?xpacket end=\"w\"?>\n";
    return out.str();
}

[[nodiscard]] std::string streamObject(std::string_view data, std::string_view extraDictionary = {}) {
    std::ostringstream out;
    out << "<< /Length " << data.size();
    if (!extraDictionary.empty())
        out << ' ' << extraDictionary;
    out << " >>\nstream\n";
    std::string result = out.str();
    result.append(data);
    if (result.empty() || result.back() != '\n')
        result.push_back('\n');
    result += "endstream";
    return result;
}

[[nodiscard]] std::string buildPdf(
    const GraphicsScene& scene,
    std::string content,
    const PdfTimestamp& timestamp,
    const PdfId& documentId,
    const PdfId& instanceId,
    const std::optional<std::string>& userName) {
    const std::string version = MMCAL_VERSION_STRING;
    const std::string title = "mmCal " + version + " Plot";
    const std::string creator = "mmCal " + version;
    const std::string producer = "mmCal " + version + " PDF Plotter";
    const auto titlePdf = pdfAsciiString(title).value();
    const auto creatorPdf = pdfAsciiString(creator).value();
    const auto producerPdf = pdfAsciiString(producer).value();

    const std::string xmp = xmpPacket(timestamp, documentId, instanceId, userName);
    const double widthPt = scene.extent.widthMm * pointsPerMillimeter;
    const double heightPt = scene.extent.heightMm * pointsPerMillimeter;

    std::vector<std::string> objects(9);
    objects[0] = "<< /Type /Catalog /Pages 2 0 R /Metadata 8 0 R /Lang (ja-JP) >>";
    objects[1] = "<< /Type /Pages /Kids [3 0 R] /Count 1 >>";
    {
        std::ostringstream page;
        page << "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 "
            << number(widthPt) << ' ' << number(heightPt) << "] "
            << "/Resources << /ProcSet [/PDF /Text] /Font << "
            << "/F1 5 0 R /F2 6 0 R /F3 7 0 R >> >> "
            << "/Contents 4 0 R >>";
        objects[2] = page.str();
    }
    objects[3] = streamObject(content);
    objects[4] = "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>";
    objects[5] = "<< /Type /Font /Subtype /Type1 /BaseFont /Times-Roman /Encoding /WinAnsiEncoding >>";
    objects[6] = "<< /Type /Font /Subtype /Type1 /BaseFont /Courier /Encoding /WinAnsiEncoding >>";
    objects[7] = streamObject(xmp, "/Type /Metadata /Subtype /XML");
    {
        std::ostringstream info;
        info << "<< /Title (" << titlePdf << ")"
            << " /Creator (" << creatorPdf << ")"
            << " /Producer (" << producerPdf << ")"
            << " /CreationDate (" << timestamp.pdf << ")"
            << " /ModDate (" << timestamp.pdf << ") >>";
        objects[8] = info.str();
    }

    std::string result = "%PDF-1.4\n%\xE2\xE3\xCF\xD3\n";
    std::vector<std::size_t> offsets(objects.size() + 1, 0);
    for (std::size_t i = 0; i < objects.size(); ++i) {
        offsets[i + 1] = result.size();
        result += std::to_string(i + 1) + " 0 obj\n";
        result += objects[i];
        result += "\nendobj\n";
    }

    const std::size_t xrefOffset = result.size();
    std::ostringstream xref;
    xref << "xref\n0 " << (objects.size() + 1) << "\n"
        << "0000000000 65535 f \n";
    for (std::size_t i = 1; i < offsets.size(); ++i)
        xref << std::setw(10) << std::setfill('0') << offsets[i] << " 00000 n \n";
    xref << "trailer\n<< /Size " << (objects.size() + 1)
        << " /Root 1 0 R /Info 9 0 R /ID [<" << idHex(documentId)
        << "><" << idHex(instanceId) << ">] >>\n"
        << "startxref\n" << xrefOffset << "\n%%EOF\n";
    result += xref.str();
    return result;
}

} // namespace

PdfRenderResult renderPdf(
    const GraphicsScene& scene,
    const PdfRenderOptions& options) {
    PdfRenderResult result;
    if (!(scene.extent.widthMm > 0.0) || !(scene.extent.heightMm > 0.0)
        || !std::isfinite(scene.extent.widthMm) || !std::isfinite(scene.extent.heightMm))
        return result;

    for (const auto& node : scene.nodes) {
        if (hasUnsupportedTransparency(node)) {
            result.status = PdfRenderStatus::UnsupportedTransparency;
            return result;
        }
        if (hasUnsupportedText(node)) {
            result.status = PdfRenderStatus::UnsupportedText;
            return result;
        }
    }

    std::ostringstream content;
    content << "% mmCal PDF 1.4 content stream; user space below is millimeters.\n"
        << "q\n"
        << number(pointsPerMillimeter) << " 0 0 " << number(pointsPerMillimeter) << " 0 0 cm\n"
        // 漸近線などcanvas外geometryを埋め込み先へ漏らさない。
        << "0 0 " << number(scene.extent.widthMm) << ' '
        << number(scene.extent.heightMm) << " re W n\n";

    for (const auto& node : scene.nodes) {
        const bool written = std::visit([&](const auto& value) {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, GraphicsPathNode>)
                return writePath(content, value);
            else if constexpr (std::is_same_v<T, GraphicsEllipseNode>)
                return writeEllipse(content, value);
            else if constexpr (std::is_same_v<T, GraphicsEllipticArcNode>)
                return writeEllipticArc(content, value);
            else if constexpr (std::is_same_v<T, GraphicsCircleNode>)
                return writeCircle(content, value);
            else
                return writeText(content, value);
        }, node);
        if (!written) {
            result.status = PdfRenderStatus::NonFiniteGeometry;
            return result;
        }
    }
    content << "Q\n";

    const PdfTimestamp timestamp = makeTimestamp(options.deterministic);
    const PdfId documentId = makeId(options.deterministic, 1);
    const PdfId instanceId = makeId(options.deterministic, 2);
    const std::optional<std::string> userName = options.userName
        ? options.userName : currentUserName();

    try {
        result.pdf = buildPdf(
            scene, content.str(), timestamp, documentId, instanceId, userName);
        result.status = PdfRenderStatus::Success;
    }
    catch (...) {
        result.status = PdfRenderStatus::RenderFailed;
        result.pdf.reset();
    }
    return result;
}

} // namespace mmcal::graphics
