#pragma once

#include "raster_backend.hpp"

#include <optional>
#include <string>

namespace mmcal::graphics {

enum class PngRenderStatus {
    Success,
    InvalidScene,
    InvalidOptions,
    UnsupportedText,
    ResourceLimit,
    EncodeFailed
};

struct PngRenderResult final {
    PngRenderStatus status = PngRenderStatus::InvalidScene;
    std::optional<std::string> png;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PngRenderStatus::Success && png.has_value();
    }
};

// PNGはRGBA8，adaptive scanline filter，zlib/DEFLATEを自前生成し，dynamic/fixed Huffmanから小さい方を採用する。
[[nodiscard]] PngRenderResult renderPng(
    const GraphicsScene& scene,
    const RasterRenderOptions& options = {});

} // namespace mmcal::graphics
