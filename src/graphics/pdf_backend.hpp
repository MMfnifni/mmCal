#pragma once

#include "graphics_scene.hpp"

#include <optional>
#include <string>

namespace mmcal::graphics {

enum class PdfRenderStatus {
    Success,
    InvalidScene,
    NonFiniteGeometry,
    UnsupportedTransparency,
    UnsupportedText,
    RenderFailed
};

struct PdfRenderOptions final {
    // byte-for-byte回帰用。日時・DocumentID・InstanceIDを固定する。
    bool deterministic = false;
    // テストや埋め込み環境ではOS取得を避けて明示値を渡せる。
    std::optional<std::string> userName;
};

struct PdfRenderResult final {
    PdfRenderStatus status = PdfRenderStatus::InvalidScene;
    std::optional<std::string> pdf;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PdfRenderStatus::Success && pdf.has_value();
    }
};

// PDF 1.4 / classic xref / 非圧縮object・streamで出力する。
// GraphicsSceneはmm・y上向きのまま扱い，backend先頭のCTMでpointへ変換する。
[[nodiscard]] PdfRenderResult renderPdf(
    const GraphicsScene& scene,
    const PdfRenderOptions& options = {});

} // namespace mmcal::graphics
