#pragma once

#include "plot_analysis.hpp"
#include "plot_exact_geometry.hpp"
#include "plot_program.hpp"
#include "plot_request.hpp"

#include "evaluation/builtin_registry.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/big_float.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::plot {

enum class PlotPointTag : std::uint32_t {
    None = 0,
    Coarse = 1u << 0,
    Adaptive = 1u << 1,
    XAxisIntercept = 1u << 2,
    YAxisIntercept = 1u << 3,
    LocalExtremum = 1u << 4,
    CurveIntersection = 1u << 5
};

[[nodiscard]] constexpr PlotPointTag operator|(PlotPointTag lhs, PlotPointTag rhs) noexcept {
    return static_cast<PlotPointTag>(
        static_cast<std::uint32_t>(lhs) | static_cast<std::uint32_t>(rhs));
}

inline PlotPointTag& operator|=(PlotPointTag& lhs, PlotPointTag rhs) noexcept {
    lhs = lhs | rhs;
    return lhs;
}

[[nodiscard]] constexpr bool hasPlotPointTag(PlotPointTag value, PlotPointTag tag) noexcept {
    return (static_cast<std::uint32_t>(value) & static_cast<std::uint32_t>(tag)) != 0;
}

// 1次元parameterから2次元曲線へ写したsampling点。
// 通常のPlotでは parameter == x だが，ParametricPlotでは両者を分離する。
// 非有限値でもparameter/xは保持し，geometryの有効性はstatusで表す。
struct CurveSample2D final {
    numeric::BigFloat parameter;
    numeric::BigFloat x;
    numeric::BigFloat y;
    PlotNumericStatus status = PlotNumericStatus::Undefined;
    PlotPointTag tags = PlotPointTag::None;

    CurveSample2D(
        numeric::BigFloat parameterValue,
        numeric::BigFloat xValue,
        numeric::BigFloat yValue,
        PlotNumericStatus sampleStatus,
        PlotPointTag sampleTags = PlotPointTag::None)
        : parameter(std::move(parameterValue)),
          x(std::move(xValue)),
          y(std::move(yValue)),
          status(sampleStatus),
          tags(sampleTags) {}

    // graph y=f(x)用の互換constructor。parameterとgeometry xを同じ値にする。
    CurveSample2D(
        numeric::BigFloat xValue,
        numeric::BigFloat yValue,
        PlotNumericStatus sampleStatus,
        PlotPointTag sampleTags = PlotPointTag::None)
        : parameter(xValue),
          x(std::move(xValue)),
          y(std::move(yValue)),
          status(sampleStatus),
          tags(sampleTags) {}

    [[nodiscard]] bool finite() const noexcept {
        return status == PlotNumericStatus::Finite;
    }
};

// 既存Plot内部APIの互換名。新規共通処理はCurveSample2Dを使う。
using PlotSample = CurveSample2D;

enum class PlotSegmentGeometryKind {
    Polyline,
    StraightLine,
    QuadraticBezier,
    CubicBezier,
    Ellipse,
    EllipticArc,
    PiecewiseConstant
};

// open domain境界の有限な片側極限を，曲線sampleとは分離して保持する。
// 数学的には未定義点なのでsamplesへ混ぜず，endpoint markerの位置だけに使う。
struct PlotEndpointPoint final {
    numeric::BigFloat x;
    numeric::BigFloat y;
};

struct PlotBezierControlPoints final {
    PlotEndpointPoint control1;
    PlotEndpointPoint control2;
};

struct PlotEllipseGeometry final {
    PlotEndpointPoint center;
    PlotEndpointPoint cosineAxis;
    PlotEndpointPoint sineAxis;
};

struct PlotEllipticArcGeometry final {
    PlotEndpointPoint center;
    PlotEndpointPoint cosineAxis;
    PlotEndpointPoint sineAxis;
    numeric::Rational startTurns;
    numeric::Rational sweepTurns;
};

// PlotAnalysisで分割された一つの数学区間に対応するsampling結果。
// affine/constant函数はStraightLine，step函数はPiecewiseConstantとして明示し，
// 後段で不要な細分化や不連続点を跨ぐ補間を避ける。
struct SampledCurveSegment2D final {
    PlotInterval sourceInterval;
    std::vector<CurveSample2D> samples;
    PlotSegmentGeometryKind geometryKind = PlotSegmentGeometryKind::Polyline;
    std::optional<PlotEndpointPoint> quadraticControlPoint;
    std::optional<PlotBezierControlPoints> bezierControlPoints;
    std::optional<PlotEllipseGeometry> ellipseGeometry;
    std::optional<PlotEllipticArcGeometry> ellipticArcGeometry;
    // requestそのものの描画端にはendpoint markerを出さない。
    // domain分割やstep函数等で生じた内部境界だけを可視化するための情報を保持する。
    bool markLowerEndpoint = false;
    bool markUpperEndpoint = false;
    std::optional<PlotEndpointPoint> lowerEndpointPoint;
    std::optional<PlotEndpointPoint> upperEndpointPoint;
    // adaptive samplingでopen boundaryへ物理解像度基準で近づくための数値parameter。
    // endpoint自体の函数値は評価せず，parameterだけをsampling時に確定して保持する。
    std::optional<numeric::BigFloat> lowerBoundaryParameter;
    std::optional<numeric::BigFloat> upperBoundaryParameter;
    // asymptoteまたは有限極限が証明できたopen boundaryだけを専用refinement対象にする。
    // sin[1/x]のようなoscillatory holeへ無意味に近づいて性能を落とさない。
    bool approachLowerBoundary = false;
    bool approachUpperBoundary = false;
};

struct SampledCurve2D final {
    std::vector<SampledCurveSegment2D> segments;
    std::size_t precisionBits = 0;
};

// 既存Plot内部APIの互換名。ParametricPlot導入時はSampledCurve2Dを直接共有する。
using SampledCurveSegment = SampledCurveSegment2D;
using SampledCurve = SampledCurve2D;

struct CoarseSamplingOptions final {
    std::size_t precisionBits = 64;
    std::size_t samplesPerInterval = 49;
    std::size_t maxTotalSamples = 4096;
    // exact geometryを座標列へmaterializeする内部consumer専用。
    // PiecewiseConstantは不連続意味論を保つためこの指定でも専用geometryを維持する。
    bool forcePolylineGeometry = false;
};

enum class PlotSamplingStatus {
    Success,
    InvalidOptions,
    EndpointCompilationFailed,
    EndpointEvaluationFailed,
    ProgramInitializationFailed,
    InvalidInterval,
    PrecisionInsufficient,
    ResourceLimit
};

struct PlotSamplingResult final {
    PlotSamplingStatus status = PlotSamplingStatus::InvalidOptions;
    std::optional<SampledCurve> curve;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PlotSamplingStatus::Success && curve.has_value();
    }
};

// PlotAnalysisの区間を跨がず，各区間へ低密度の初期sampleを配置する。
// open endpointは評価しないため，証明済みpole/branch boundaryを誤ってsampleしない。
[[nodiscard]] PlotSamplingResult coarseSamplePlot(
    const PlotRequest& request,
    const PlotAnalysis& analysis,
    const PlotProgram& program,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const CoarseSamplingOptions& options = {});

// ParametricPlotの共通domain上で，同一parameterをx/yの2本のPlotProgramへ渡して
// CurveSample2Dを構築する。domain分割は両coordinateの解析結果を交差したものを受け取る。
[[nodiscard]] PlotSamplingResult coarseSampleParametricPlot(
    const ParametricCurveRequest& request,
    const PlotDomain& domain,
    const PlotProgram& xProgram,
    const PlotProgram& yProgram,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const CoarseSamplingOptions& options = {},
    const ExactCurveGeometry* exactGeometry = nullptr);

} // namespace mmcal::plot
