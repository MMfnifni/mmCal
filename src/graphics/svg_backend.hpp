#pragma once

#include "graphics_scene.hpp"

#include <optional>
#include <string>

namespace mmcal::graphics {

enum class SvgRenderStatus {
    Success,
    InvalidScene,
    NonFiniteGeometry
};

struct SvgRenderResult final {
    SvgRenderStatus status = SvgRenderStatus::InvalidScene;
    std::optional<std::string> svg;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == SvgRenderStatus::Success && svg.has_value();
    }
};

// GraphicsSceneのmm座標をそのままSVG viewBox unitへ対応させる。
// Graphicsのy上向きからSVGのy下向きへの反転はbackend内だけで行う。
[[nodiscard]] SvgRenderResult renderSvg(const GraphicsScene& scene);

} // namespace mmcal::graphics
