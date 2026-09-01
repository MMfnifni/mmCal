#pragma once

#include <string_view>

namespace mmcal::cli {

enum class OutputLayout {
    Auto,
    Single,
    Multi
};

[[nodiscard]] inline std::string_view outputLayoutName(OutputLayout layout) noexcept {
    switch (layout) {
    case OutputLayout::Auto: return "Auto";
    case OutputLayout::Single: return "Single";
    case OutputLayout::Multi: return "Multi";
    }
    return "Auto";
}

} // namespace mmcal::cli
