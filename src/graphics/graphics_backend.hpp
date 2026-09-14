#pragma once

#include "graphics_scene.hpp"
#include "raster_backend.hpp"

#include <optional>
#include <string>
#include <string_view>

namespace mmcal::graphics {

enum class GraphicsFormat {
    Svg,
    Eps,
    Pdf,
    Png,
    Webp
};

enum class GraphicsRenderStatus {
    Success,
    InvalidScene,
    InvalidOptions,
    UnsupportedFeature,
    RenderFailed
};

struct GraphicsRenderOptions final {
    RasterRenderOptions raster;
};

struct GraphicsRenderResult final {
    GraphicsRenderStatus status = GraphicsRenderStatus::InvalidScene;
    std::optional<std::string> data;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == GraphicsRenderStatus::Success && data.has_value();
    }
};

[[nodiscard]] std::optional<GraphicsFormat> parseGraphicsFormat(std::string_view format);
[[nodiscard]] std::optional<GraphicsFormat> graphicsFormatFromExtension(std::string_view extension);
[[nodiscard]] const char* graphicsFormatName(GraphicsFormat format) noexcept;
[[nodiscard]] const char* graphicsFormatExtension(GraphicsFormat format) noexcept;

// GraphicsSceneをbackendへ送る共通入口。Plot側は出力形式を意識しない。
[[nodiscard]] GraphicsRenderResult renderGraphics(
    const GraphicsScene& scene,
    GraphicsFormat format,
    const GraphicsRenderOptions& options = {});

} // namespace mmcal::graphics
