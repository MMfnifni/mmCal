#include "graphics_backend.hpp"

#include "eps_backend.hpp"
#include "pdf_backend.hpp"
#include "svg_backend.hpp"

#include <algorithm>
#include <cctype>
#include <string>
#include <utility>

namespace mmcal::graphics {
namespace {

[[nodiscard]] std::string asciiLower(std::string_view value) {
    std::string result{value};
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return result;
}

} // namespace

std::optional<GraphicsFormat> parseGraphicsFormat(std::string_view format) {
    const std::string value = asciiLower(format);
    if (value == "svg")
        return GraphicsFormat::Svg;
    if (value == "eps")
        return GraphicsFormat::Eps;
    if (value == "pdf")
        return GraphicsFormat::Pdf;
    return std::nullopt;
}

std::optional<GraphicsFormat> graphicsFormatFromExtension(std::string_view extension) {
    std::string value = asciiLower(extension);
    if (!value.empty() && value.front() == '.')
        value.erase(value.begin());
    return parseGraphicsFormat(value);
}

const char* graphicsFormatName(GraphicsFormat format) noexcept {
    switch (format) {
    case GraphicsFormat::Svg: return "SVG";
    case GraphicsFormat::Eps: return "EPS";
    case GraphicsFormat::Pdf: return "PDF";
    }
    return "Unknown";
}

const char* graphicsFormatExtension(GraphicsFormat format) noexcept {
    switch (format) {
    case GraphicsFormat::Svg: return ".svg";
    case GraphicsFormat::Eps: return ".eps";
    case GraphicsFormat::Pdf: return ".pdf";
    }
    return "";
}

GraphicsRenderResult renderGraphics(const GraphicsScene& scene, GraphicsFormat format) {
    switch (format) {
    case GraphicsFormat::Svg: {
        auto rendered = renderSvg(scene);
        if (rendered)
            return GraphicsRenderResult{GraphicsRenderStatus::Success, std::move(rendered.svg)};
        return GraphicsRenderResult{
            rendered.status == SvgRenderStatus::InvalidScene
                ? GraphicsRenderStatus::InvalidScene
                : GraphicsRenderStatus::RenderFailed,
            std::nullopt};
    }
    case GraphicsFormat::Eps: {
        auto rendered = renderEps(scene);
        if (rendered)
            return GraphicsRenderResult{GraphicsRenderStatus::Success, std::move(rendered.eps)};
        if (rendered.status == EpsRenderStatus::InvalidScene)
            return GraphicsRenderResult{GraphicsRenderStatus::InvalidScene, std::nullopt};
        if (rendered.status == EpsRenderStatus::UnsupportedTransparency
            || rendered.status == EpsRenderStatus::UnsupportedText)
            return GraphicsRenderResult{GraphicsRenderStatus::UnsupportedFeature, std::nullopt};
        return GraphicsRenderResult{GraphicsRenderStatus::RenderFailed, std::nullopt};
    }
    case GraphicsFormat::Pdf: {
        auto rendered = renderPdf(scene);
        if (rendered)
            return GraphicsRenderResult{GraphicsRenderStatus::Success, std::move(rendered.pdf)};
        if (rendered.status == PdfRenderStatus::InvalidScene)
            return GraphicsRenderResult{GraphicsRenderStatus::InvalidScene, std::nullopt};
        if (rendered.status == PdfRenderStatus::UnsupportedTransparency
            || rendered.status == PdfRenderStatus::UnsupportedText)
            return GraphicsRenderResult{GraphicsRenderStatus::UnsupportedFeature, std::nullopt};
        return GraphicsRenderResult{GraphicsRenderStatus::RenderFailed, std::nullopt};
    }
    }
    return GraphicsRenderResult{GraphicsRenderStatus::RenderFailed, std::nullopt};
}

} // namespace mmcal::graphics
