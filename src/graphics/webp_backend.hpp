#pragma once

#include "raster_backend.hpp"

#include <optional>
#include <string>

namespace mmcal::graphics {

enum class WebpRenderStatus {
    Success,
    InvalidScene,
    InvalidOptions,
    UnsupportedText,
    ResourceLimit,
    EncodeFailed
};

struct WebpRenderResult final {
    WebpRenderStatus status = WebpRenderStatus::InvalidScene;
    std::optional<std::string> webp;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == WebpRenderStatus::Success && webp.has_value();
    }
};

// WebP Lossless (VP8L) のcorrectness-first encoder。
// literal・color cache・bounded LZ77 backward referenceに加え，Predictor / Subtract Green Transformを扱う。
// Predictorは16x16 blockごとに14 modeを評価し，Subtract GreenはR/BからGをmod-256で減算する。
// None / Predictor / Subtract Green / 両者併用を実encodeし，最小VP8L payloadだけ採用する。
// 短いmatchでは1 pixel lazy matchingを行い，greedy継続と次tokenまで比較する。
// backward distanceはVP8Lの近傍2D mappingを優先し，一般scan-line mappingへfallbackする。
// color cacheは16 entryの仕様準拠hash cacheを用いる。
// G/length/cache, R, B, A, distanceの5 prefix treeは実token頻度からcanonical Huffman codeを構築し，
// code-length列も16/17/18 repeat codeで圧縮してVP8Lへ格納する。
[[nodiscard]] WebpRenderResult renderWebp(
    const GraphicsScene& scene,
    const RasterRenderOptions& options = {});

} // namespace mmcal::graphics
