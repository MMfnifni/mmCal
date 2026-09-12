#pragma once

#include "plot_sampling.hpp"

#include "numeric/big_float.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace mmcal::plot {

using PlotAnchorId = std::uint64_t;

enum class PlotAnchorKind : std::uint32_t {
    None = 0,
    XAxisIntercept = 1u << 0,
    YAxisIntercept = 1u << 1,
    LocalExtremum = 1u << 2,
    CurveIntersection = 1u << 3
};

[[nodiscard]] constexpr PlotAnchorKind operator|(PlotAnchorKind lhs, PlotAnchorKind rhs) noexcept {
    return static_cast<PlotAnchorKind>(
        static_cast<std::uint32_t>(lhs) | static_cast<std::uint32_t>(rhs));
}

inline PlotAnchorKind& operator|=(PlotAnchorKind& lhs, PlotAnchorKind rhs) noexcept {
    lhs = lhs | rhs;
    return lhs;
}

[[nodiscard]] constexpr bool hasPlotAnchorKind(PlotAnchorKind value, PlotAnchorKind kind) noexcept {
    return (static_cast<std::uint32_t>(value) & static_cast<std::uint32_t>(kind)) != 0;
}

struct PlotAnchorVertexRef final {
    std::size_t curveIndex = 0;
    std::size_t segmentIndex = 0;
    std::size_t sampleIndex = 0;
};

// PlotSceneへ渡す共有semantic point。座標が一致する複数curveの意味点を一つへ束ねる。
struct PlotAnchor final {
    PlotAnchorId id = 0;
    PlotAnchorKind kinds = PlotAnchorKind::None;
    numeric::BigFloat x;
    numeric::BigFloat y;
    std::vector<std::size_t> curveIndices;
    std::vector<PlotAnchorVertexRef> vertices;
};

struct PlotAnchorSet final {
    std::vector<PlotAnchor> anchors;
};

// Coarse/Adaptive等の純sampling tagは捨て，意味を持つvertexだけを共有anchorへ昇格する。
[[nodiscard]] PlotAnchorSet collectPlotAnchors(const std::vector<SampledCurve>& curves);

} // namespace mmcal::plot
