#pragma once

#include <cctype>
#include <cmath>
#include <optional>
#include <string_view>

namespace mmcal::graphics {

struct GraphicsTextMetrics final {
    double advanceWidthMm = 0.0;
    double ascentMm = 0.0;
    double descentMm = 0.0;

    [[nodiscard]] double heightMm() const noexcept { return ascentMm + descentMm; }
};

// Font backend導入前のlayout用近似値。tick labelはASCII数値中心なので，
// 一律0.6emより文字種別幅を持たせた方がmargin推定の過不足が小さい。
[[nodiscard]] inline std::optional<GraphicsTextMetrics> measureTextApproximate(
    std::string_view text,
    double fontSizeMm) {
    if (!(fontSizeMm > 0.0) || !std::isfinite(fontSizeMm))
        return std::nullopt;

    double em = 0.0;
    for (const unsigned char ch : text) {
        if (std::isdigit(ch))
            em += 0.56;
        else if (ch == '.' || ch == ',')
            em += 0.28;
        else if (ch == '-')
            em += 0.52;
        else if (ch == '+')
            em += 0.36;
        else if (ch == 'e' || ch == 'E')
            em += 0.52;
        else if (ch == ' ')
            em += 0.28;
        else
            em += 0.60;
    }
    return GraphicsTextMetrics{em * fontSizeMm, 0.78 * fontSizeMm, 0.22 * fontSizeMm};
}

} // namespace mmcal::graphics
