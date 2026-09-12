// Plotの初期samplingとSampledCurve IR
#include "plot_sampling.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "numeric/rounding_mode.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <optional>
#include <utility>

namespace mmcal::plot {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;
using numeric::RoundingMode;

[[nodiscard]] std::optional<BigFloat> evaluateEndpoint(
    const expression::Expr& expression,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits,
    PlotSamplingStatus& failure) {
    const PlotCompileResult compiled = compilePlotProgram(
        expression, request.variable, builtins, mathematics);
    if (!compiled) {
        failure = PlotSamplingStatus::EndpointCompilationFailed;
        return std::nullopt;
    }
    if (compiled.program->variableDependent[compiled.program->resultRegister]) {
        failure = PlotSamplingStatus::EndpointCompilationFailed;
        return std::nullopt;
    }

    try {
        BigFloatPlotExecutor executor{*compiled.program, precisionBits, angles};
        const BigFloat zero = BigFloat::fromBigInt(
            BigInt{0}, precisionBits, RoundingMode::NearestEven);
        const PlotNumericResult result = executor.evaluate(zero);
        if (!result.finite()) {
            failure = PlotSamplingStatus::EndpointEvaluationFailed;
            return std::nullopt;
        }
        return result.value;
    }
    catch (...) {
        failure = PlotSamplingStatus::EndpointEvaluationFailed;
        return std::nullopt;
    }
}

struct SampleFraction final {
    std::size_t numerator = 0;
    std::size_t denominator = 1;
};

[[nodiscard]] SampleFraction sampleFraction(
    std::size_t index,
    std::size_t count,
    PlotEndpointInclusion lower,
    PlotEndpointInclusion upper) {
    const bool lowerClosed = lower == PlotEndpointInclusion::Closed;
    const bool upperClosed = upper == PlotEndpointInclusion::Closed;

    if (lowerClosed && upperClosed)
        return SampleFraction{index, count - 1};
    if (lowerClosed)
        return SampleFraction{index, count};
    if (upperClosed)
        return SampleFraction{index + 1, count};
    return SampleFraction{index + 1, count + 1};
}

// StraightLineは2点だけで十分だが，generic samplingと同じfraction規則を使うと
// open endpoint側が区間の半分までしか伸びない。未定義端点そのものは評価せず，
// 描画上ほぼ端点まで届くexact rational insetを使う。
[[nodiscard]] SampleFraction straightLineSampleFraction(
    std::size_t index,
    PlotEndpointInclusion lower,
    PlotEndpointInclusion upper) {
    constexpr std::size_t insetDenominator = 1024;
    if (index == 0)
        return lower == PlotEndpointInclusion::Closed
            ? SampleFraction{0, 1}
            : SampleFraction{1, insetDenominator};
    return upper == PlotEndpointInclusion::Closed
        ? SampleFraction{1, 1}
        : SampleFraction{insetDenominator - 1, insetDenominator};
}


// step函数の境界値はsin[Pi]のように数値評価の丸めだけで隣の整数段へ落ち得る。
// 非singleton pieceは包含状態に関係なく両端を僅かに内側から評価し，区間内部の
// 定数値だけで水平線を構成する。開閉markerのx位置はlowering側で復元する。
[[nodiscard]] SampleFraction piecewiseConstantSampleFraction(std::size_t index) {
    constexpr std::size_t insetDenominator = 1024;
    return index == 0
        ? SampleFraction{1, insetDenominator}
        : SampleFraction{insetDenominator - 1, insetDenominator};
}

[[nodiscard]] BigFloat interpolate(
    const BigFloat& lower,
    const BigFloat& upper,
    SampleFraction fraction,
    std::size_t precisionBits) {
    if (fraction.numerator == 0)
        return lower;
    if (fraction.numerator == fraction.denominator)
        return upper;

    const BigFloat width = numeric::subtract(
        upper, lower, precisionBits, RoundingMode::NearestEven);
    const Rational exactFraction{
        BigInt{static_cast<std::int64_t>(fraction.numerator)},
        BigInt{static_cast<std::int64_t>(fraction.denominator)}};
    const BigFloat t = BigFloat::fromRational(
        exactFraction, precisionBits, RoundingMode::NearestEven);
    const BigFloat offset = numeric::multiply(
        width, t, precisionBits, RoundingMode::NearestEven);
    return numeric::add(
        lower, offset, precisionBits, RoundingMode::NearestEven);
}

[[nodiscard]] PlotSample evaluateSample(
    BigFloatPlotExecutor& executor,
    BigFloat x) {
    PlotNumericResult result = executor.evaluate(x);
    return PlotSample{std::move(x), std::move(result.value), result.status, PlotPointTag::Coarse};
}

[[nodiscard]] std::optional<PlotEndpointPoint> quadraticBezierControlPoint(
    const PlotSample& p0,
    const PlotSample& midpoint,
    const PlotSample& p2,
    std::size_t precisionBits) {
    if (!p0.finite() || !midpoint.finite() || !p2.finite())
        return std::nullopt;
    const BigFloat two = BigFloat::fromBigInt(
        BigInt{2}, precisionBits, RoundingMode::NearestEven);
    const BigFloat half = BigFloat::fromRational(
        Rational{BigInt{1}, BigInt{2}}, precisionBits, RoundingMode::NearestEven);
    const BigFloat controlX = numeric::multiply(
        numeric::add(p0.x, p2.x, precisionBits, RoundingMode::NearestEven),
        half, precisionBits, RoundingMode::NearestEven);
    const BigFloat endpointAverage = numeric::multiply(
        numeric::add(p0.y, p2.y, precisionBits, RoundingMode::NearestEven),
        half, precisionBits, RoundingMode::NearestEven);
    const BigFloat controlY = numeric::subtract(
        numeric::multiply(midpoint.y, two, precisionBits, RoundingMode::NearestEven),
        endpointAverage, precisionBits, RoundingMode::NearestEven);
    return PlotEndpointPoint{controlX, controlY};
}

[[nodiscard]] std::optional<PlotBezierControlPoints> cubicBezierControlPoints(
    const PlotSample& p0,
    const PlotSample& oneThird,
    const PlotSample& twoThirds,
    const PlotSample& p3,
    std::size_t precisionBits) {
    if (!p0.finite() || !oneThird.finite() || !twoThirds.finite() || !p3.finite())
        return std::nullopt;

    const BigFloat eighteen = BigFloat::fromBigInt(
        BigInt{18}, precisionBits, RoundingMode::NearestEven);
    const auto scaled = [&](const BigFloat& value, std::int64_t factor) {
        return numeric::multiply(
            value,
            BigFloat::fromBigInt(BigInt{factor}, precisionBits, RoundingMode::NearestEven),
            precisionBits, RoundingMode::NearestEven);
    };

    const BigFloat a = numeric::subtract(
        numeric::subtract(
            scaled(oneThird.y, 27), scaled(p0.y, 8),
            precisionBits, RoundingMode::NearestEven),
        p3.y, precisionBits, RoundingMode::NearestEven);
    const BigFloat b = numeric::subtract(
        numeric::subtract(
            scaled(twoThirds.y, 27), p0.y,
            precisionBits, RoundingMode::NearestEven),
        scaled(p3.y, 8), precisionBits, RoundingMode::NearestEven);
    const BigFloat c1y = numeric::divide(
        numeric::subtract(scaled(a, 2), b, precisionBits, RoundingMode::NearestEven),
        eighteen, precisionBits, RoundingMode::NearestEven);
    const BigFloat c2y = numeric::divide(
        numeric::subtract(scaled(b, 2), a, precisionBits, RoundingMode::NearestEven),
        eighteen, precisionBits, RoundingMode::NearestEven);

    const BigFloat three = BigFloat::fromBigInt(
        BigInt{3}, precisionBits, RoundingMode::NearestEven);
    const BigFloat c1x = numeric::divide(
        numeric::add(scaled(p0.x, 2), p3.x, precisionBits, RoundingMode::NearestEven),
        three, precisionBits, RoundingMode::NearestEven);
    const BigFloat c2x = numeric::divide(
        numeric::add(p0.x, scaled(p3.x, 2), precisionBits, RoundingMode::NearestEven),
        three, precisionBits, RoundingMode::NearestEven);
    return PlotBezierControlPoints{{c1x, c1y}, {c2x, c2y}};
}

[[nodiscard]] std::optional<BigFloat> evaluateParametricEndpoint(
    const expression::Expr& expression,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits,
    PlotSamplingStatus& failure) {
    const PlotCompileResult compiled = compilePlotProgram(
        expression, request.parameter, builtins, mathematics);
    if (!compiled || compiled.program->variableDependent[compiled.program->resultRegister]) {
        failure = PlotSamplingStatus::EndpointCompilationFailed;
        return std::nullopt;
    }
    try {
        BigFloatPlotExecutor executor{*compiled.program, precisionBits, angles};
        const BigFloat zero = BigFloat::fromBigInt(
            BigInt{0}, precisionBits, RoundingMode::NearestEven);
        const PlotNumericResult result = executor.evaluate(zero);
        if (!result.finite()) {
            failure = PlotSamplingStatus::EndpointEvaluationFailed;
            return std::nullopt;
        }
        return result.value;
    }
    catch (...) {
        failure = PlotSamplingStatus::EndpointEvaluationFailed;
        return std::nullopt;
    }
}


[[nodiscard]] std::optional<PlotEndpointPoint> evaluateSymbolicPoint(
    const SymbolicPoint2D& point,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits,
    PlotSamplingStatus& failure) {
    const auto x = evaluateParametricEndpoint(
        point.x, request, builtins, mathematics, angles, precisionBits, failure);
    if (!x)
        return std::nullopt;
    const auto y = evaluateParametricEndpoint(
        point.y, request, builtins, mathematics, angles, precisionBits, failure);
    if (!y)
        return std::nullopt;
    return PlotEndpointPoint{*x, *y};
}

[[nodiscard]] PlotSegmentGeometryKind exactGeometryKind(const ExactCurveGeometry& geometry) {
    return std::visit([](const auto& value) -> PlotSegmentGeometryKind {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, ExactPointGeometry>
            || std::is_same_v<T, ExactLineGeometry>)
            return PlotSegmentGeometryKind::StraightLine;
        else if constexpr (std::is_same_v<T, ExactQuadraticBezierGeometry>)
            return PlotSegmentGeometryKind::QuadraticBezier;
        else if constexpr (std::is_same_v<T, ExactCubicBezierGeometry>)
            return PlotSegmentGeometryKind::CubicBezier;
        else if constexpr (std::is_same_v<T, ExactEllipticArcGeometry>)
            return PlotSegmentGeometryKind::EllipticArc;
        else {
            static_assert(std::is_same_v<T, ExactCircleGeometry>
                || std::is_same_v<T, ExactEllipseGeometry>);
            return PlotSegmentGeometryKind::Ellipse;
        }
    }, geometry.value);
}

[[nodiscard]] PlotNumericStatus combinedStatus(
    const PlotNumericResult& x,
    const PlotNumericResult& y) noexcept {
    if (!x.finite())
        return x.status;
    return y.status;
}

[[nodiscard]] CurveSample2D evaluateParametricSample(
    BigFloatPlotExecutor& xExecutor,
    BigFloatPlotExecutor& yExecutor,
    BigFloat parameter,
    PlotPointTag tags = PlotPointTag::Coarse) {
    auto x = xExecutor.evaluate(parameter);
    auto y = yExecutor.evaluate(parameter);
    return CurveSample2D{
        std::move(parameter), std::move(x.value), std::move(y.value),
        combinedStatus(x, y), tags};
}

[[nodiscard]] int polynomialDegreeBound(PlotCurveGeometryKind kind) noexcept {
    switch (kind) {
    case PlotCurveGeometryKind::Constant:
    case PlotCurveGeometryKind::PiecewiseConstant:
        return 0;
    case PlotCurveGeometryKind::Affine:
        return 1;
    case PlotCurveGeometryKind::QuadraticPolynomial:
        return 2;
    case PlotCurveGeometryKind::CubicPolynomial:
        return 3;
    case PlotCurveGeometryKind::Generic:
        return 4;
    }
    return 4;
}

[[nodiscard]] PlotSegmentGeometryKind parametricGeometryKind(
    const PlotProgram& xProgram,
    const PlotProgram& yProgram,
    const PlotInterval& interval) noexcept {
    const int degree = std::max(
        polynomialDegreeBound(xProgram.geometryKind),
        polynomialDegreeBound(yProgram.geometryKind));
    if (degree <= 1)
        return PlotSegmentGeometryKind::StraightLine;
    if (interval.lowerInclusion == PlotEndpointInclusion::Closed
        && interval.upperInclusion == PlotEndpointInclusion::Closed) {
        if (degree == 2)
            return PlotSegmentGeometryKind::QuadraticBezier;
        if (degree == 3)
            return PlotSegmentGeometryKind::CubicBezier;
    }
    return PlotSegmentGeometryKind::Polyline;
}

[[nodiscard]] BigFloat quadraticBezierCoordinate(
    const BigFloat& p0,
    const BigFloat& midpoint,
    const BigFloat& p2,
    std::size_t precisionBits) {
    const BigFloat two = BigFloat::fromBigInt(
        BigInt{2}, precisionBits, RoundingMode::NearestEven);
    const BigFloat half = BigFloat::fromRational(
        Rational{BigInt{1}, BigInt{2}}, precisionBits, RoundingMode::NearestEven);
    const BigFloat endpointAverage = numeric::multiply(
        numeric::add(p0, p2, precisionBits, RoundingMode::NearestEven),
        half, precisionBits, RoundingMode::NearestEven);
    return numeric::subtract(
        numeric::multiply(midpoint, two, precisionBits, RoundingMode::NearestEven),
        endpointAverage, precisionBits, RoundingMode::NearestEven);
}

[[nodiscard]] std::optional<PlotEndpointPoint> parametricQuadraticBezierControlPoint(
    const CurveSample2D& p0,
    const CurveSample2D& midpoint,
    const CurveSample2D& p2,
    std::size_t precisionBits) {
    if (!p0.finite() || !midpoint.finite() || !p2.finite())
        return std::nullopt;
    return PlotEndpointPoint{
        quadraticBezierCoordinate(p0.x, midpoint.x, p2.x, precisionBits),
        quadraticBezierCoordinate(p0.y, midpoint.y, p2.y, precisionBits)};
}

[[nodiscard]] std::pair<BigFloat, BigFloat> cubicBezierCoordinateControls(
    const BigFloat& p0,
    const BigFloat& oneThird,
    const BigFloat& twoThirds,
    const BigFloat& p3,
    std::size_t precisionBits) {
    const BigFloat eighteen = BigFloat::fromBigInt(
        BigInt{18}, precisionBits, RoundingMode::NearestEven);
    const auto scaled = [&](const BigFloat& value, std::int64_t factor) {
        return numeric::multiply(
            value,
            BigFloat::fromBigInt(BigInt{factor}, precisionBits, RoundingMode::NearestEven),
            precisionBits, RoundingMode::NearestEven);
    };
    const BigFloat a = numeric::subtract(
        numeric::subtract(
            scaled(oneThird, 27), scaled(p0, 8),
            precisionBits, RoundingMode::NearestEven),
        p3, precisionBits, RoundingMode::NearestEven);
    const BigFloat b = numeric::subtract(
        numeric::subtract(
            scaled(twoThirds, 27), p0, precisionBits, RoundingMode::NearestEven),
        scaled(p3, 8), precisionBits, RoundingMode::NearestEven);
    return {
        numeric::divide(
            numeric::subtract(scaled(a, 2), b, precisionBits, RoundingMode::NearestEven),
            eighteen, precisionBits, RoundingMode::NearestEven),
        numeric::divide(
            numeric::subtract(scaled(b, 2), a, precisionBits, RoundingMode::NearestEven),
            eighteen, precisionBits, RoundingMode::NearestEven)};
}

[[nodiscard]] std::optional<PlotBezierControlPoints> parametricCubicBezierControlPoints(
    const CurveSample2D& p0,
    const CurveSample2D& oneThird,
    const CurveSample2D& twoThirds,
    const CurveSample2D& p3,
    std::size_t precisionBits) {
    if (!p0.finite() || !oneThird.finite() || !twoThirds.finite() || !p3.finite())
        return std::nullopt;
    const auto x = cubicBezierCoordinateControls(
        p0.x, oneThird.x, twoThirds.x, p3.x, precisionBits);
    const auto y = cubicBezierCoordinateControls(
        p0.y, oneThird.y, twoThirds.y, p3.y, precisionBits);
    return PlotBezierControlPoints{{x.first, y.first}, {x.second, y.second}};
}

[[nodiscard]] bool needsOpenBoundaryApproach(
    const PlotAnalysis& analysis,
    const expression::Expr& position) {
    return std::any_of(
        analysis.landmarks.begin(), analysis.landmarks.end(),
        [&](const PlotLandmark& landmark) {
            if (landmark.position != position)
                return false;
            return landmark.kind == PlotLandmarkKind::Pole
                || landmark.kind == PlotLandmarkKind::VerticalAsymptote
                || landmark.kind == PlotLandmarkKind::RemovableSingularity;
        });
}

} // namespace

PlotSamplingResult coarseSamplePlot(
    const PlotRequest& request,
    const PlotAnalysis& analysis,
    const PlotProgram& program,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const CoarseSamplingOptions& options) {
    if (options.precisionBits < 8
        || options.samplesPerInterval < 2
        || options.maxTotalSamples == 0)
        return PlotSamplingResult{PlotSamplingStatus::InvalidOptions, std::nullopt};

    const bool forcePolyline = options.forcePolylineGeometry
        && program.geometryKind != PlotCurveGeometryKind::PiecewiseConstant;
    const std::size_t plannedSamplesPerInterval = forcePolyline
        || program.geometryKind == PlotCurveGeometryKind::Generic
        || program.geometryKind == PlotCurveGeometryKind::QuadraticPolynomial
        || program.geometryKind == PlotCurveGeometryKind::CubicPolynomial
        ? options.samplesPerInterval
        : std::size_t{2};
    if (!analysis.domain.intervals.empty()
        && plannedSamplesPerInterval > options.maxTotalSamples / analysis.domain.intervals.size())
        return PlotSamplingResult{PlotSamplingStatus::ResourceLimit, std::nullopt};

    std::optional<BigFloatPlotExecutor> executor;
    try {
        executor.emplace(program, options.precisionBits, angles);
    }
    catch (...) {
        return PlotSamplingResult{
            PlotSamplingStatus::ProgramInitializationFailed, std::nullopt};
    }

    SampledCurve curve;
    curve.precisionBits = options.precisionBits;
    curve.segments.reserve(analysis.domain.intervals.size());

    PlotSamplingStatus requestEndpointFailure = PlotSamplingStatus::Success;
    const auto requestedLower = evaluateEndpoint(
        request.lower, request, builtins, mathematics, angles,
        options.precisionBits, requestEndpointFailure);
    if (!requestedLower)
        return PlotSamplingResult{requestEndpointFailure, std::nullopt};
    const auto requestedUpper = evaluateEndpoint(
        request.upper, request, builtins, mathematics, angles,
        options.precisionBits, requestEndpointFailure);
    if (!requestedUpper)
        return PlotSamplingResult{requestEndpointFailure, std::nullopt};

    std::size_t totalSamples = 0;
    for (const PlotInterval& interval : analysis.domain.intervals) {
        PlotSamplingStatus endpointFailure = PlotSamplingStatus::Success;
        const auto lower = evaluateEndpoint(
            interval.lower, request, builtins, mathematics, angles,
            options.precisionBits, endpointFailure);
        if (!lower)
            return PlotSamplingResult{endpointFailure, std::nullopt};
        const auto upper = evaluateEndpoint(
            interval.upper, request, builtins, mathematics, angles,
            options.precisionBits, endpointFailure);
        if (!upper)
            return PlotSamplingResult{endpointFailure, std::nullopt};

        if (*lower > *upper)
            return PlotSamplingResult{PlotSamplingStatus::InvalidInterval, std::nullopt};

        PlotSegmentGeometryKind segmentGeometry = PlotSegmentGeometryKind::Polyline;
        if (!forcePolyline && (program.geometryKind == PlotCurveGeometryKind::Constant
            || program.geometryKind == PlotCurveGeometryKind::Affine))
            segmentGeometry = PlotSegmentGeometryKind::StraightLine;
        else if (!forcePolyline && program.geometryKind == PlotCurveGeometryKind::QuadraticPolynomial
            && interval.lowerInclusion == PlotEndpointInclusion::Closed
            && interval.upperInclusion == PlotEndpointInclusion::Closed)
            segmentGeometry = PlotSegmentGeometryKind::QuadraticBezier;
        else if (!forcePolyline && program.geometryKind == PlotCurveGeometryKind::CubicPolynomial
            && interval.lowerInclusion == PlotEndpointInclusion::Closed
            && interval.upperInclusion == PlotEndpointInclusion::Closed)
            segmentGeometry = PlotSegmentGeometryKind::CubicBezier;
        else if (program.geometryKind == PlotCurveGeometryKind::PiecewiseConstant)
            segmentGeometry = PlotSegmentGeometryKind::PiecewiseConstant;

        // 内部の特異点だけでなく，request端そのものがopen domain境界である場合も
        // 曲線を物理解像度まで境界へ近づける。x^x on {0,3} のような右極限holeは
        // これが無いと最初の粗sample(≈1/49)から始まり，不自然に軸から離れて見える。
        const bool lowerAtRequestedBoundary = *lower == *requestedLower;
        const bool upperAtRequestedBoundary = *upper == *requestedUpper;
        SampledCurveSegment segment{
            interval, {}, segmentGeometry, std::nullopt, std::nullopt, std::nullopt,
            std::nullopt,
            *lower != *requestedLower,
            *upper != *requestedUpper,
            std::nullopt, std::nullopt,
            *lower, *upper,
            interval.lowerInclusion == PlotEndpointInclusion::Open
                && (lowerAtRequestedBoundary || needsOpenBoundaryApproach(analysis, interval.lower)),
            interval.upperInclusion == PlotEndpointInclusion::Open
                && (upperAtRequestedBoundary || needsOpenBoundaryApproach(analysis, interval.upper))};

        if (*lower == *upper) {
            if (interval.lower != interval.upper)
                return PlotSamplingResult{
                    PlotSamplingStatus::PrecisionInsufficient, std::nullopt};
            if (interval.lowerInclusion != PlotEndpointInclusion::Closed
                || interval.upperInclusion != PlotEndpointInclusion::Closed)
                return PlotSamplingResult{PlotSamplingStatus::InvalidInterval, std::nullopt};
            if (++totalSamples > options.maxTotalSamples)
                return PlotSamplingResult{PlotSamplingStatus::ResourceLimit, std::nullopt};
            segment.samples.push_back(evaluateSample(*executor, *lower));
            curve.segments.push_back(std::move(segment));
            continue;
        }

        const std::size_t count = [&] {
            if (segmentGeometry == PlotSegmentGeometryKind::Polyline
                || segmentGeometry == PlotSegmentGeometryKind::QuadraticBezier
                || segmentGeometry == PlotSegmentGeometryKind::CubicBezier)
                return options.samplesPerInterval;
            return std::size_t{2};
        }();
        if (totalSamples > options.maxTotalSamples - count)
            return PlotSamplingResult{PlotSamplingStatus::ResourceLimit, std::nullopt};
        totalSamples += count;
        segment.samples.reserve(count);

        for (std::size_t i = 0; i < count; ++i) {
            SampleFraction fraction;
            if (segmentGeometry == PlotSegmentGeometryKind::Polyline
                || segmentGeometry == PlotSegmentGeometryKind::QuadraticBezier
                || segmentGeometry == PlotSegmentGeometryKind::CubicBezier)
                fraction = sampleFraction(
                    i, count, interval.lowerInclusion, interval.upperInclusion);
            else if (segmentGeometry == PlotSegmentGeometryKind::PiecewiseConstant)
                fraction = piecewiseConstantSampleFraction(i);
            else
                fraction = straightLineSampleFraction(
                    i, interval.lowerInclusion, interval.upperInclusion);
            BigFloat x = interpolate(*lower, *upper, fraction, options.precisionBits);
            segment.samples.push_back(evaluateSample(*executor, std::move(x)));
        }
        if (segmentGeometry == PlotSegmentGeometryKind::QuadraticBezier) {
            PlotSample middle = evaluateSample(
                *executor, interpolate(*lower, *upper, SampleFraction{1, 2}, options.precisionBits));
            segment.quadraticControlPoint = quadraticBezierControlPoint(
                segment.samples.front(), middle, segment.samples.back(), options.precisionBits);
            if (!segment.quadraticControlPoint)
                segment.geometryKind = PlotSegmentGeometryKind::Polyline;
        }
        else if (segmentGeometry == PlotSegmentGeometryKind::CubicBezier) {
            PlotSample oneThird = evaluateSample(
                *executor, interpolate(*lower, *upper, SampleFraction{1, 3}, options.precisionBits));
            PlotSample twoThirds = evaluateSample(
                *executor, interpolate(*lower, *upper, SampleFraction{2, 3}, options.precisionBits));
            segment.bezierControlPoints = cubicBezierControlPoints(
                segment.samples.front(), oneThird, twoThirds,
                segment.samples.back(), options.precisionBits);
            if (!segment.bezierControlPoints)
                segment.geometryKind = PlotSegmentGeometryKind::Polyline;
        }
        curve.segments.push_back(std::move(segment));
    }

    return PlotSamplingResult{PlotSamplingStatus::Success, std::move(curve)};
}

PlotSamplingResult coarseSampleParametricPlot(
    const ParametricCurveRequest& request,
    const PlotDomain& domain,
    const PlotProgram& xProgram,
    const PlotProgram& yProgram,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const CoarseSamplingOptions& options,
    const ExactCurveGeometry* exactGeometry) {
    if (options.precisionBits < 8 || options.samplesPerInterval < 2
        || options.maxTotalSamples == 0)
        return PlotSamplingResult{PlotSamplingStatus::InvalidOptions, std::nullopt};
    if (!domain.intervals.empty()
        && options.samplesPerInterval > options.maxTotalSamples / domain.intervals.size())
        return PlotSamplingResult{PlotSamplingStatus::ResourceLimit, std::nullopt};

    std::optional<BigFloatPlotExecutor> xExecutor;
    std::optional<BigFloatPlotExecutor> yExecutor;
    try {
        xExecutor.emplace(xProgram, options.precisionBits, angles);
        yExecutor.emplace(yProgram, options.precisionBits, angles);
    }
    catch (...) {
        return PlotSamplingResult{
            PlotSamplingStatus::ProgramInitializationFailed, std::nullopt};
    }

    PlotSamplingStatus endpointFailure = PlotSamplingStatus::Success;
    const auto requestedLower = evaluateParametricEndpoint(
        request.lower, request, builtins, mathematics, angles,
        options.precisionBits, endpointFailure);
    if (!requestedLower)
        return PlotSamplingResult{endpointFailure, std::nullopt};
    const auto requestedUpper = evaluateParametricEndpoint(
        request.upper, request, builtins, mathematics, angles,
        options.precisionBits, endpointFailure);
    if (!requestedUpper)
        return PlotSamplingResult{endpointFailure, std::nullopt};

    SampledCurve2D curve;
    curve.precisionBits = options.precisionBits;
    curve.segments.reserve(domain.intervals.size());
    std::size_t totalSamples = 0;

    for (const PlotInterval& interval : domain.intervals) {
        endpointFailure = PlotSamplingStatus::Success;
        const auto lower = evaluateParametricEndpoint(
            interval.lower, request, builtins, mathematics, angles,
            options.precisionBits, endpointFailure);
        if (!lower)
            return PlotSamplingResult{endpointFailure, std::nullopt};
        const auto upper = evaluateParametricEndpoint(
            interval.upper, request, builtins, mathematics, angles,
            options.precisionBits, endpointFailure);
        if (!upper)
            return PlotSamplingResult{endpointFailure, std::nullopt};
        if (*lower > *upper)
            return PlotSamplingResult{PlotSamplingStatus::InvalidInterval, std::nullopt};

        const PlotSegmentGeometryKind geometry = options.forcePolylineGeometry
            ? PlotSegmentGeometryKind::Polyline
            : exactGeometry && domain.intervals.size() == 1
                ? exactGeometryKind(*exactGeometry)
                : parametricGeometryKind(xProgram, yProgram, interval);
        SampledCurveSegment2D segment{
            interval, {}, geometry, std::nullopt, std::nullopt, std::nullopt,
            std::nullopt,
            *lower != *requestedLower, *upper != *requestedUpper,
            std::nullopt, std::nullopt, *lower, *upper, false, false};

        if (*lower == *upper) {
            if (interval.lower != interval.upper)
                return PlotSamplingResult{
                    PlotSamplingStatus::PrecisionInsufficient, std::nullopt};
            if (interval.lowerInclusion != PlotEndpointInclusion::Closed
                || interval.upperInclusion != PlotEndpointInclusion::Closed)
                return PlotSamplingResult{PlotSamplingStatus::InvalidInterval, std::nullopt};
            if (++totalSamples > options.maxTotalSamples)
                return PlotSamplingResult{PlotSamplingStatus::ResourceLimit, std::nullopt};
            segment.samples.push_back(evaluateParametricSample(
                *xExecutor, *yExecutor, *lower));
            curve.segments.push_back(std::move(segment));
            continue;
        }

        const std::size_t count = geometry == PlotSegmentGeometryKind::StraightLine
            ? std::size_t{2} : options.samplesPerInterval;
        if (totalSamples > options.maxTotalSamples - count)
            return PlotSamplingResult{PlotSamplingStatus::ResourceLimit, std::nullopt};
        totalSamples += count;
        segment.samples.reserve(count);

        for (std::size_t i = 0; i < count; ++i) {
            const SampleFraction fraction = geometry == PlotSegmentGeometryKind::StraightLine
                ? straightLineSampleFraction(i, interval.lowerInclusion, interval.upperInclusion)
                : sampleFraction(i, count, interval.lowerInclusion, interval.upperInclusion);
            BigFloat parameter = interpolate(
                *lower, *upper, fraction, options.precisionBits);
            segment.samples.push_back(evaluateParametricSample(
                *xExecutor, *yExecutor, std::move(parameter)));
        }

        if (geometry == PlotSegmentGeometryKind::QuadraticBezier) {
            if (exactGeometry) {
                if (const auto* quadratic = std::get_if<ExactQuadraticBezierGeometry>(
                        &exactGeometry->value)) {
                    endpointFailure = PlotSamplingStatus::Success;
                    segment.quadraticControlPoint = evaluateSymbolicPoint(
                        quadratic->control, request, builtins, mathematics, angles,
                        options.precisionBits, endpointFailure);
                }
            }
            if (!segment.quadraticControlPoint) {
                auto midpointSample = evaluateParametricSample(
                    *xExecutor, *yExecutor,
                    interpolate(*lower, *upper, SampleFraction{1, 2}, options.precisionBits));
                segment.quadraticControlPoint = parametricQuadraticBezierControlPoint(
                    segment.samples.front(), midpointSample, segment.samples.back(),
                    options.precisionBits);
            }
            if (!segment.quadraticControlPoint)
                segment.geometryKind = PlotSegmentGeometryKind::Polyline;
        }
        else if (geometry == PlotSegmentGeometryKind::CubicBezier) {
            if (exactGeometry) {
                if (const auto* cubic = std::get_if<ExactCubicBezierGeometry>(
                        &exactGeometry->value)) {
                    endpointFailure = PlotSamplingStatus::Success;
                    const auto control1 = evaluateSymbolicPoint(
                        cubic->control1, request, builtins, mathematics, angles,
                        options.precisionBits, endpointFailure);
                    const auto control2 = control1 ? evaluateSymbolicPoint(
                        cubic->control2, request, builtins, mathematics, angles,
                        options.precisionBits, endpointFailure) : std::nullopt;
                    if (control1 && control2)
                        segment.bezierControlPoints = PlotBezierControlPoints{*control1, *control2};
                }
            }
            if (!segment.bezierControlPoints) {
                auto oneThird = evaluateParametricSample(
                    *xExecutor, *yExecutor,
                    interpolate(*lower, *upper, SampleFraction{1, 3}, options.precisionBits));
                auto twoThirds = evaluateParametricSample(
                    *xExecutor, *yExecutor,
                    interpolate(*lower, *upper, SampleFraction{2, 3}, options.precisionBits));
                segment.bezierControlPoints = parametricCubicBezierControlPoints(
                    segment.samples.front(), oneThird, twoThirds, segment.samples.back(),
                    options.precisionBits);
            }
            if (!segment.bezierControlPoints)
                segment.geometryKind = PlotSegmentGeometryKind::Polyline;
        }
        else if (geometry == PlotSegmentGeometryKind::Ellipse) {
            if (exactGeometry) {
                const SymbolicPoint2D* center = nullptr;
                const SymbolicPoint2D* cosineAxis = nullptr;
                const SymbolicPoint2D* sineAxis = nullptr;
                if (const auto* circle = std::get_if<ExactCircleGeometry>(&exactGeometry->value)) {
                    center = &circle->center;
                    cosineAxis = &circle->cosineAxis;
                    sineAxis = &circle->sineAxis;
                }
                else if (const auto* ellipse = std::get_if<ExactEllipseGeometry>(&exactGeometry->value)) {
                    center = &ellipse->center;
                    cosineAxis = &ellipse->cosineAxis;
                    sineAxis = &ellipse->sineAxis;
                }
                if (center && cosineAxis && sineAxis) {
                    endpointFailure = PlotSamplingStatus::Success;
                    const auto c = evaluateSymbolicPoint(
                        *center, request, builtins, mathematics, angles,
                        options.precisionBits, endpointFailure);
                    const auto a = c ? evaluateSymbolicPoint(
                        *cosineAxis, request, builtins, mathematics, angles,
                        options.precisionBits, endpointFailure) : std::nullopt;
                    const auto b = a ? evaluateSymbolicPoint(
                        *sineAxis, request, builtins, mathematics, angles,
                        options.precisionBits, endpointFailure) : std::nullopt;
                    if (c && a && b)
                        segment.ellipseGeometry = PlotEllipseGeometry{*c, *a, *b};
                }
            }
            if (!segment.ellipseGeometry)
                segment.geometryKind = PlotSegmentGeometryKind::Polyline;
        }
        else if (geometry == PlotSegmentGeometryKind::EllipticArc) {
            if (exactGeometry) {
                if (const auto* arc = std::get_if<ExactEllipticArcGeometry>(&exactGeometry->value)) {
                    endpointFailure = PlotSamplingStatus::Success;
                    const auto c = evaluateSymbolicPoint(
                        arc->center, request, builtins, mathematics, angles,
                        options.precisionBits, endpointFailure);
                    const auto a = c ? evaluateSymbolicPoint(
                        arc->cosineAxis, request, builtins, mathematics, angles,
                        options.precisionBits, endpointFailure) : std::nullopt;
                    const auto b = a ? evaluateSymbolicPoint(
                        arc->sineAxis, request, builtins, mathematics, angles,
                        options.precisionBits, endpointFailure) : std::nullopt;
                    if (c && a && b) {
                        segment.ellipticArcGeometry = PlotEllipticArcGeometry{
                            *c, *a, *b, arc->startTurns, arc->sweepTurns};
                    }
                }
            }
            if (!segment.ellipticArcGeometry)
                segment.geometryKind = PlotSegmentGeometryKind::Polyline;
        }
        curve.segments.push_back(std::move(segment));
    }

    return PlotSamplingResult{PlotSamplingStatus::Success, std::move(curve)};
}

} // namespace mmcal::plot
