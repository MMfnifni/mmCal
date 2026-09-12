#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace mmcal::graphics {

struct GraphicsPointMm final {
    double xMm = 0.0;
    double yMm = 0.0;
};

struct GraphicsExtentMm final {
    double widthMm = 0.0;
    double heightMm = 0.0;
};

struct GraphicsRectMm final {
    double xMm = 0.0;
    double yMm = 0.0;
    double widthMm = 0.0;
    double heightMm = 0.0;
};

struct GraphicsColor final {
    std::uint8_t red = 0;
    std::uint8_t green = 0;
    std::uint8_t blue = 0;
    std::uint8_t alpha = 255;
};

enum class GraphicsLineCap {
    Butt,
    Round,
    Square
};

enum class GraphicsLineJoin {
    Miter,
    Round,
    Bevel
};

struct GraphicsStrokeStyle final {
    GraphicsColor color;
    double widthMm = 0.25;
    GraphicsLineCap lineCap = GraphicsLineCap::Round;
    GraphicsLineJoin lineJoin = GraphicsLineJoin::Round;
};

struct GraphicsFillStyle final {
    GraphicsColor color;
};

struct GraphicsMoveTo final { GraphicsPointMm point; };
struct GraphicsLineTo final { GraphicsPointMm point; };
struct GraphicsQuadraticTo final {
    GraphicsPointMm control;
    GraphicsPointMm point;
};
struct GraphicsCubicTo final {
    GraphicsPointMm control1;
    GraphicsPointMm control2;
    GraphicsPointMm point;
};
struct GraphicsClosePath final {};
using GraphicsPathCommand = std::variant<
    GraphicsMoveTo, GraphicsLineTo, GraphicsQuadraticTo, GraphicsCubicTo, GraphicsClosePath>;

enum class GraphicsSemanticKind {
    Curve,
    CurveEndpoint,
    Point,
    XAxis,
    YAxis,
    XTick,
    YTick,
    XTickLabel,
    YTickLabel,
    Anchor,
    Label
};

struct GraphicsSemanticRef final {
    GraphicsSemanticKind kind = GraphicsSemanticKind::Curve;
    std::uint64_t id = 0;
};

struct GraphicsPathNode final {
    std::vector<GraphicsPathCommand> commands;
    std::optional<GraphicsStrokeStyle> stroke;
    std::optional<GraphicsFillStyle> fill;
    std::optional<GraphicsSemanticRef> semantic;
    std::optional<GraphicsRectMm> clipRect;
};


struct GraphicsEllipseNode final {
    GraphicsPointMm center;
    // unit circle (cos t, sin t) をmm空間へ写す2本の列vector。
    // 直交を要求しないため，回転楕円を含む任意のaffine imageを保持できる。
    GraphicsPointMm cosineAxis;
    GraphicsPointMm sineAxis;
    std::optional<GraphicsStrokeStyle> stroke;
    std::optional<GraphicsFillStyle> fill;
    std::optional<GraphicsSemanticRef> semantic;
    std::optional<GraphicsRectMm> clipRect;
};

struct GraphicsEllipticArcNode final {
    GraphicsPointMm center;
    GraphicsPointMm cosineAxis;
    GraphicsPointMm sineAxis;
    double startRadians = 0.0;
    double sweepRadians = 0.0;
    std::optional<GraphicsStrokeStyle> stroke;
    std::optional<GraphicsFillStyle> fill;
    std::optional<GraphicsSemanticRef> semantic;
    std::optional<GraphicsRectMm> clipRect;
};

struct GraphicsCircleNode final {
    GraphicsPointMm center;
    double radiusMm = 0.0;
    std::optional<GraphicsStrokeStyle> stroke;
    std::optional<GraphicsFillStyle> fill;
    std::optional<GraphicsSemanticRef> semantic;
    std::optional<GraphicsRectMm> clipRect;
};

enum class GraphicsTextAnchor {
    Start,
    Middle,
    End
};

struct GraphicsTextNode final {
    GraphicsPointMm origin;
    std::string text;
    double fontSizeMm = 3.5;
    std::string fontFamily = "sans-serif";
    GraphicsTextAnchor anchor = GraphicsTextAnchor::Start;
    GraphicsColor color;
    std::optional<GraphicsSemanticRef> semantic;
};

using GraphicsNode = std::variant<
    GraphicsPathNode, GraphicsEllipseNode, GraphicsEllipticArcNode,
    GraphicsCircleNode, GraphicsTextNode>;

// backend非依存の物理描画scene。座標・線幅・文字寸法はすべてmm。
struct GraphicsScene final {
    GraphicsExtentMm extent;
    std::vector<GraphicsNode> nodes;
};

} // namespace mmcal::graphics
