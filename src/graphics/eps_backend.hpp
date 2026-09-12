#pragma once

#include "graphics_scene.hpp"

#include <optional>
#include <string>

namespace mmcal::graphics {

enum class EpsRenderStatus {
    Success,
    InvalidScene,
    NonFiniteGeometry,
    UnsupportedTransparency,
    UnsupportedText
};

struct EpsRenderResult final {
    EpsRenderStatus status = EpsRenderStatus::InvalidScene;
    std::optional<std::string> eps;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == EpsRenderStatus::Success && eps.has_value();
    }
};

// EPSはPostScriptのy上向き座標をそのまま使えるため，GraphicsSceneのmm座標を
// 冒頭のCTMでmm->ptへ変換して描画する。QuadraticはCubicへ厳密にdegree elevationする。
[[nodiscard]] EpsRenderResult renderEps(const GraphicsScene& scene);

} // namespace mmcal::graphics
