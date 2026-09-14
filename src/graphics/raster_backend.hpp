#pragma once

#include "graphics_scene.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace mmcal::graphics {

enum class RasterBackground {
    White,
    None
};

struct RasterRenderOptions final {
    // 0ならdpiとsceneの物理寸法から決定する。片方だけの指定は許可しない。
    std::uint32_t widthPx = 0;
    std::uint32_t heightPx = 0;
    double dpi = 254.0;
    std::uint32_t antialiasing = 2;
    RasterBackground background = RasterBackground::White;
};

struct RasterImage final {
    std::uint32_t widthPx = 0;
    std::uint32_t heightPx = 0;
    // PNG pHYs等へそのまま渡せる実効解像度。ImageSize指定ではx/yが異なり得る。
    double dpiX = 0.0;
    double dpiY = 0.0;
    std::vector<std::uint8_t> rgba;
};

enum class RasterRenderStatus {
    Success,
    InvalidScene,
    InvalidOptions,
    UnsupportedText,
    ResourceLimit,
    RenderFailed
};

struct RasterRenderResult final {
    RasterRenderStatus status = RasterRenderStatus::InvalidScene;
    std::optional<RasterImage> image;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == RasterRenderStatus::Success && image.has_value();
    }
};

// GraphicsSceneをsquare-pixel rasterへ落とす共通backend。
// AAはsupersampling + premultiplied-alpha box downsampleで行う。
[[nodiscard]] RasterRenderResult renderRaster(
    const GraphicsScene& scene,
    const RasterRenderOptions& options = {});

} // namespace mmcal::graphics
