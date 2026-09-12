#include "plot_pipeline.hpp"
#include "plot_periodicity.hpp"
#include "plot_exact_geometry.hpp"
#include "plot_exact_range.hpp"

#include "graphics/graphics_backend.hpp"
#include "plot_program.hpp"
#include "plot_sample_analysis.hpp"
#include "plot_sampling.hpp"
#include "simplification/simplifier.hpp"
#include "simplification/simplification_context.hpp"
#include "symbolic/limit.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rounding_mode.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>

namespace mmcal::plot {
namespace {

[[nodiscard]] PlotPipelineResult failure(
    PlotPipelineStatus status,
    std::size_t curveIndex = std::numeric_limits<std::size_t>::max()) {
    PlotPipelineResult result;
    result.status = status;
    result.failedCurveIndex = curveIndex;
    return result;
}


[[nodiscard]] std::optional<PlotRequestOptions> mergedPlotRequestOptions(
    const std::vector<PlotRequest>& requests) {
    PlotRequestOptions merged;
    for (const auto& request : requests) {
        if (request.options.plotRange) {
            if (merged.plotRange && *merged.plotRange != *request.options.plotRange)
                return std::nullopt;
            merged.plotRange = request.options.plotRange;
        }
        if (request.options.aspectRatio) {
            if (merged.aspectRatio && *merged.aspectRatio != *request.options.aspectRatio)
                return std::nullopt;
            merged.aspectRatio = request.options.aspectRatio;
        }
        if (request.options.ticks) {
            if (merged.ticks && *merged.ticks != *request.options.ticks)
                return std::nullopt;
            merged.ticks = request.options.ticks;
        }
    }
    return merged;
}


[[nodiscard]] std::size_t scaledSamplingCount(
    std::size_t base, std::size_t plotPoints, std::size_t minimum) {
    constexpr std::size_t defaultPlotPoints = 100;
    if (plotPoints == defaultPlotPoints)
        return std::max(base, minimum);

    const std::size_t max = std::numeric_limits<std::size_t>::max();
    if (base > max / plotPoints)
        return max;
    const std::size_t product = base * plotPoints;
    const std::size_t scaled = product > max - (defaultPlotPoints - 1)
        ? max : (product + defaultPlotPoints - 1) / defaultPlotPoints;
    return std::max(scaled, minimum);
}

[[nodiscard]] CoarseSamplingOptions samplingDensityOptions(
    const CoarseSamplingOptions& base, const PlotRequestOptions& requestOptions,
    bool densityApplies) {
    if (!densityApplies || !requestOptions.plotPoints)
        return base;
    CoarseSamplingOptions result = base;
    result.samplesPerInterval = scaledSamplingCount(
        base.samplesPerInterval, *requestOptions.plotPoints, 2);
    result.maxTotalSamples = scaledSamplingCount(
        base.maxTotalSamples, *requestOptions.plotPoints, 1);
    return result;
}

[[nodiscard]] AdaptiveSamplingOptions adaptiveDensityOptions(
    const AdaptiveSamplingOptions& base, const PlotRequestOptions& requestOptions,
    bool densityApplies) {
    if (!densityApplies || !requestOptions.plotPoints)
        return base;
    AdaptiveSamplingOptions result = base;
    // PlotPointsは誤差許容値や再帰深度を変えない。初期mesh増加で既定resource capだけを
    // 同率に広げ，指定密度そのものがResourceLimitになることを避ける。
    result.maxTotalSamples = scaledSamplingCount(
        base.maxTotalSamples, *requestOptions.plotPoints, 1);
    return result;
}

[[nodiscard]] double bigFloatToDouble(const numeric::BigFloat& value) {
    if (value.isZero())
        return 0.0;
    long double significand = 0.0L;
    try {
        significand = std::stold(value.significand().toString());
    }
    catch (...) {
        return value.isNegative() ? -std::numeric_limits<double>::infinity()
                                  : std::numeric_limits<double>::infinity();
    }
    const auto exponent = value.exponent();
    if (exponent < static_cast<numeric::BigFloat::exponent_type>(
            std::numeric_limits<int>::min()))
        return 0.0;
    if (exponent > static_cast<numeric::BigFloat::exponent_type>(
            std::numeric_limits<int>::max()))
        return value.isNegative() ? -std::numeric_limits<double>::infinity()
                                  : std::numeric_limits<double>::infinity();
    return static_cast<double>(std::ldexp(significand, static_cast<int>(exponent)));
}
[[nodiscard]] std::optional<numeric::BigFloat> finiteConstantValue(
    const expression::Expr& expression,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits) {
    const auto compiled = compilePlotProgram(
        expression, request.variable, builtins, mathematics);
    if (!compiled
        || compiled.program->variableDependent[compiled.program->resultRegister])
        return std::nullopt;

    try {
        BigFloatPlotExecutor executor{*compiled.program, precisionBits, angles};
        const auto zero = numeric::BigFloat::fromBigInt(
            numeric::BigInt{0}, precisionBits, numeric::RoundingMode::NearestEven);
        const auto value = executor.evaluate(zero);
        if (!value.finite())
            return std::nullopt;
        return value.value;
    }
    catch (...) {
        return std::nullopt;
    }
}

void annotateFiniteOpenEndpointLimits(
    const PlotRequest& request,
    const PlotAnalysis& analysis,
    SampledCurve& curve,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions) {
    // generic曲線のopen boundaryはpath sampleへ未定義点を混ぜない。その代わり，
    // 元式の片側極限が有限と証明できた場合だけmarker用の独立点を保持する。
    // poleやsin[1/x]のように極限が発散・不定・未解決なら何も追加しない。
    for (auto& segment : curve.segments) {
        if (segment.geometryKind != PlotSegmentGeometryKind::Polyline)
            continue;

        const auto annotate = [&](bool lowerEndpoint) {
            const auto inclusion = lowerEndpoint
                ? segment.sourceInterval.lowerInclusion
                : segment.sourceInterval.upperInclusion;
            if (inclusion != PlotEndpointInclusion::Open)
                return;

            const expression::Expr& point = lowerEndpoint
                ? segment.sourceInterval.lower : segment.sourceInterval.upper;
            const symbolic::LimitDirection direction = lowerEndpoint
                ? symbolic::LimitDirection::Right : symbolic::LimitDirection::Left;

            std::optional<expression::Expr> limit;
            bool provenDivergence = false;
            for (const auto& landmark : analysis.landmarks) {
                if (landmark.position != point)
                    continue;
                if (landmark.finiteLimit)
                    limit = *landmark.finiteLimit;
                if (landmark.kind == PlotLandmarkKind::Pole
                    || landmark.kind == PlotLandmarkKind::VerticalAsymptote)
                    provenDivergence = true;
            }
            // exact prepassが発散を証明済みなら，一般Limitへ戻して巨大な有限近似値を
            // removable holeと誤認しない。有限極限metadataが明示されている場合だけ優先する。
            if (!limit && provenDivergence)
                return;
            if (!limit) {
                try {
                    limit = symbolic::limitExpression(
                        request.expression, request.variable, point, direction,
                        builtins, mathematics, angles, infinitySymbol, assumptions);
                }
                catch (...) {
                    return;
                }
            }
            const auto y = finiteConstantValue(
                *limit, request, builtins, mathematics, angles, curve.precisionBits);
            if (!y)
                return;

            // endpoint自身のxは既にsampling時と同じcertified constant evaluatorで扱える。
            const auto x = finiteConstantValue(
                point, request, builtins, mathematics, angles, curve.precisionBits);
            if (!x)
                return;
            PlotEndpointPoint marker{*x, *y};
            if (lowerEndpoint) {
                segment.lowerEndpointPoint = std::move(marker);
                segment.approachLowerBoundary = true;
                // request端でも，函数側のopen endpointに有限極限holeがあるなら
                // 白丸で示す。単なる描画打ち切り端とは区別する。
                segment.markLowerEndpoint = true;
            }
            else {
                segment.upperEndpointPoint = std::move(marker);
                segment.approachUpperBoundary = true;
                segment.markUpperEndpoint = true;
            }
        };

        annotate(true);
        annotate(false);
    }
}

} // namespace

PlotPipelineResult buildPlotPipeline(
    const std::vector<PlotRequest>& requests,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    if (requests.empty())
        return failure(PlotPipelineStatus::EmptyRequest);

    const auto requestOptions = mergedPlotRequestOptions(requests);
    if (!requestOptions)
        return failure(PlotPipelineStatus::InvalidOptions);

    std::vector<PlotAnalysis> analyses;
    std::vector<PlotProgram> programs;
    std::vector<SampledCurve> coarseCurves;
    std::vector<PlotRangeEstimate> ranges;
    analyses.reserve(requests.size());
    programs.reserve(requests.size());
    coarseCurves.reserve(requests.size());
    ranges.reserve(requests.size());

    for (std::size_t curveIndex = 0; curveIndex < requests.size(); ++curveIndex) {
        const auto& curveRequest = requests[curveIndex];
        analyses.push_back(analyzePlotRequest(
            curveRequest, builtins, mathematics, angles, infinitySymbol, assumptions));
        if (analyses.back().samplingSafety != PlotSamplingSafety::Safe)
            return failure(PlotPipelineStatus::UnsafeDiscontinuity, curveIndex);

        expression::Expr evaluationExpression = curveRequest.expression;
        if (analyses.back().domain.coverage == PlotDomainCoverage::Complete) {
            simplification::SimplificationContext simplificationContext{
                builtins, mathematics, angles, assumptions};
            // PlotAnalysisが元式の実定義域を完全に証明した後だけ，
            // その各component上で成立する恒等式を数値評価用式へ適用する。
            // 元式のdomain interval自体は保持するため，sin[1/x]^2+cos[1/x]^2 -> 1
            // としてもx=0のholeを消さない。
            simplificationContext.assumeExpressionsDefined = true;
            evaluationExpression = simplification::Simplifier{}.simplify(
                evaluationExpression, simplificationContext);
        }

        const auto compiled = compilePlotProgram(
            evaluationExpression, curveRequest.variable, builtins, mathematics);
        if (!compiled)
            return failure(PlotPipelineStatus::CompileFailed, curveIndex);
        programs.push_back(*compiled.program);

        const bool densityApplies = options.forceSampledGeometry
            || programs.back().geometryKind == PlotCurveGeometryKind::Generic;
        auto coarseOptions = samplingDensityOptions(
            options.coarseSampling, curveRequest.options, densityApplies);
        coarseOptions.forcePolylineGeometry = options.forceSampledGeometry;
        const auto sampled = coarseSamplePlot(
            curveRequest, analyses.back(), programs.back(), builtins, mathematics,
            angles, coarseOptions);
        if (!sampled)
            return failure(PlotPipelineStatus::SamplingFailed, curveIndex);
        coarseCurves.push_back(*sampled.curve);
        annotateFiniteOpenEndpointLimits(
            curveRequest, analyses.back(), coarseCurves.back(), builtins, mathematics, angles,
            infinitySymbol, assumptions);

        const auto inspection = inspectCoarseSamples(analyses.back(), coarseCurves.back());
        const auto range = estimatePlotRange(
            analyses.back(), coarseCurves.back(), inspection, options.rangeEstimator);
        if (!range)
            return failure(PlotPipelineStatus::RangeEstimationFailed, curveIndex);
        ranges.push_back(*range);
    }

    const std::size_t precisionBits = options.coarseSampling.precisionBits;
    const auto combinedRange = combinePlotRangeEstimates(
        ranges, precisionBits, options.rangeEstimator);
    if (!combinedRange)
        return failure(PlotPipelineStatus::RangeEstimationFailed);

    PlotRangeEstimate effectiveRange = *combinedRange;
    if (requestOptions->plotRange) {
        const auto minimum = finiteConstantValue(
            requestOptions->plotRange->minimum, requests.front(), builtins, mathematics, angles,
            precisionBits);
        const auto maximum = finiteConstantValue(
            requestOptions->plotRange->maximum, requests.front(), builtins, mathematics, angles,
            precisionBits);
        if (!minimum || !maximum || !(*minimum < *maximum))
            return failure(PlotPipelineStatus::InvalidOptions);
        // 明示PlotRangeにはAutomaticのpaddingを加えない。observed/data rangeは監査用に保持し，
        // viewportだけを指定された数学座標へ固定する。
        effectiveRange.viewMinimum = *minimum;
        effectiveRange.viewMaximum = *maximum;
    }

    PlotViewportMm viewport = options.viewport;
    if (requestOptions->aspectRatio) {
        const auto ratioValue = finiteConstantValue(
            *requestOptions->aspectRatio, requests.front(), builtins, mathematics, angles,
            precisionBits);
        if (!ratioValue)
            return failure(PlotPipelineStatus::InvalidOptions);
        const double ratio = bigFloatToDouble(*ratioValue);
        if (!(ratio > 0.0) || !std::isfinite(ratio))
            return failure(PlotPipelineStatus::InvalidOptions);
        // AspectRatioはheight/width。既定幅150 mmを保ち，高さだけを比率に合わせる。
        // 1 -> 1:1，0.5 -> 1:2。未指定時だけ従来の150x100 mm (2:3) を保つ。
        viewport.heightMm = viewport.widthMm * ratio;
        if (!(viewport.heightMm > 0.0) || !std::isfinite(viewport.heightMm))
            return failure(PlotPipelineStatus::InvalidOptions);
    }

    const auto transform = makePlotViewTransform(
        requests, effectiveRange, builtins, mathematics, angles,
        precisionBits, viewport);
    if (!transform)
        return failure(PlotPipelineStatus::ViewTransformFailed);

    std::vector<SampledCurve> refinedCurves;
    refinedCurves.reserve(coarseCurves.size());
    for (std::size_t curveIndex = 0; curveIndex < coarseCurves.size(); ++curveIndex) {
        const bool densityApplies = options.forceSampledGeometry
            || programs[curveIndex].geometryKind == PlotCurveGeometryKind::Generic;
        const auto adaptiveOptions = adaptiveDensityOptions(
            options.adaptiveSampling, requests[curveIndex].options, densityApplies);
        const auto refined = refinePlotSamples(
            programs[curveIndex], coarseCurves[curveIndex], *transform,
            angles, adaptiveOptions);
        if (!refined)
            return failure(PlotPipelineStatus::AdaptiveSamplingFailed, curveIndex);
        refinedCurves.push_back(*refined.curve);
    }

    if (refinedCurves.size() > 1) {
        const auto intersections = refineCurveIntersections(
            programs, refinedCurves, angles, options.intersections, &*transform);
        if (!intersections)
            return failure(PlotPipelineStatus::IntersectionFailed);
        refinedCurves = *intersections.curves;
    }

    const auto plotScene = buildPlotScene(refinedCurves);
    if (!plotScene)
        return failure(PlotPipelineStatus::SceneBuildFailed);

    PlotAxisLayoutOptions axisLayoutOptions = options.axisLayout;
    PlotGraphicsLoweringOptions graphicsLoweringOptions = options.graphicsLowering;
    if (requestOptions->ticks && !*requestOptions->ticks) {
        axisLayoutOptions.generateMajorTicks = false;
        graphicsLoweringOptions.showMajorTicks = false;
        graphicsLoweringOptions.showMajorTickLabels = false;
    }

    const auto axes = layoutPlotAxes(*transform, axisLayoutOptions);
    if (!axes)
        return failure(PlotPipelineStatus::AxisLayoutFailed);

    const auto graphicsScene = lowerPlotToGraphics(
        *plotScene.scene, *axes.layout, *transform, graphicsLoweringOptions);
    if (!graphicsScene)
        return failure(PlotPipelineStatus::GraphicsLoweringFailed);

    PlotPipelineResult result;
    result.status = PlotPipelineStatus::Success;
    result.output = PlotPipelineOutput{
        std::move(analyses),
        std::move(programs),
        std::move(refinedCurves),
        std::move(effectiveRange),
        *transform,
        *axes.layout,
        *plotScene.scene,
        *graphicsScene.scene};
    return result;
}

PlotPipelineResult buildPlotPipeline(
    const PlotRequestSet& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    return buildPlotPipeline(
        splitPlotRequests(request), builtins, mathematics, angles,
        infinitySymbol, assumptions, options);
}

PlotPipelineResult buildPlotPipeline(
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    return buildPlotPipeline(
        std::vector<PlotRequest>{request}, builtins, mathematics, angles,
        infinitySymbol, assumptions, options);
}

PlotSvgPipelineResult renderPlotSvg(
    const std::vector<PlotRequest>& requests,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    auto built = buildPlotPipeline(
        requests, builtins, mathematics, angles, infinitySymbol, assumptions, options);

    PlotSvgPipelineResult result;
    result.status = built.status;
    result.failedCurveIndex = built.failedCurveIndex;
    if (!built)
        return result;

    auto svg = graphics::renderGraphics(
        built.output->graphicsScene, graphics::GraphicsFormat::Svg);
    if (!svg || !svg.data) {
        result.status = PlotPipelineStatus::SvgRenderFailed;
        return result;
    }

    result.status = PlotPipelineStatus::Success;
    result.svg = std::move(svg.data);
    result.output = std::move(built.output);
    return result;
}

PlotSvgPipelineResult renderPlotSvg(
    const PlotRequestSet& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    return renderPlotSvg(
        splitPlotRequests(request), builtins, mathematics, angles,
        infinitySymbol, assumptions, options);
}

PlotSvgPipelineResult renderPlotSvg(
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    return renderPlotSvg(
        std::vector<PlotRequest>{request}, builtins, mathematics, angles,
        infinitySymbol, assumptions, options);
}


namespace {

[[nodiscard]] ParametricPlotPipelineResult parametricFailure(
    PlotPipelineStatus status,
    std::size_t curveIndex = std::numeric_limits<std::size_t>::max()) {
    ParametricPlotPipelineResult result;
    result.status = status;
    result.failedCurveIndex = curveIndex;
    return result;
}

[[nodiscard]] std::optional<PlotDomain> intersectParametricDomains(
    const PlotDomain& left,
    const PlotDomain& right,
    const PlotRequest& endpointRequest,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits) {
    PlotDomain result;
    result.coverage = left.coverage == PlotDomainCoverage::Complete
            && right.coverage == PlotDomainCoverage::Complete
        ? PlotDomainCoverage::Complete : PlotDomainCoverage::Unknown;

    const auto endpointValue = [&](const expression::Expr& endpoint) {
        return finiteConstantValue(
            endpoint, endpointRequest, builtins, mathematics, angles, precisionBits);
    };

    for (const auto& a : left.intervals) {
        const auto al = endpointValue(a.lower);
        const auto au = endpointValue(a.upper);
        if (!al || !au)
            return std::nullopt;
        for (const auto& b : right.intervals) {
            const auto bl = endpointValue(b.lower);
            const auto bu = endpointValue(b.upper);
            if (!bl || !bu)
                return std::nullopt;

            const bool lowerFromA = !(*al < *bl);
            const bool upperFromA = !(*bu < *au);
            const auto& lowerExpr = lowerFromA ? a.lower : b.lower;
            const auto& upperExpr = upperFromA ? a.upper : b.upper;
            const auto lowerValue = lowerFromA ? *al : *bl;
            const auto upperValue = upperFromA ? *au : *bu;
            if (upperValue < lowerValue)
                continue;

            PlotEndpointInclusion lowerInclusion;
            if (*al == *bl)
                lowerInclusion = a.lowerInclusion == PlotEndpointInclusion::Closed
                        && b.lowerInclusion == PlotEndpointInclusion::Closed
                    ? PlotEndpointInclusion::Closed : PlotEndpointInclusion::Open;
            else
                lowerInclusion = lowerFromA ? a.lowerInclusion : b.lowerInclusion;

            PlotEndpointInclusion upperInclusion;
            if (*au == *bu)
                upperInclusion = a.upperInclusion == PlotEndpointInclusion::Closed
                        && b.upperInclusion == PlotEndpointInclusion::Closed
                    ? PlotEndpointInclusion::Closed : PlotEndpointInclusion::Open;
            else
                upperInclusion = upperFromA ? a.upperInclusion : b.upperInclusion;

            if (lowerValue == upperValue
                && (lowerInclusion != PlotEndpointInclusion::Closed
                    || upperInclusion != PlotEndpointInclusion::Closed))
                continue;
            result.intervals.push_back(PlotInterval{
                lowerExpr, lowerInclusion, upperExpr, upperInclusion});
        }
    }
    return result;
}

} // namespace

ParametricPlotPipelineResult buildParametricPlotPipeline(
    const std::vector<ParametricCurveRequest>& requests,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    if (requests.empty())
        return parametricFailure(PlotPipelineStatus::EmptyRequest);

    std::vector<PlotRequest> optionRequests;
    optionRequests.reserve(requests.size());
    for (const auto& request : requests)
        optionRequests.push_back(PlotRequest{
            request.xExpression, request.parameter, request.lower, request.upper, request.options});
    const auto requestOptions = mergedPlotRequestOptions(optionRequests);
    if (!requestOptions || requestOptions->plotRange)
        return parametricFailure(PlotPipelineStatus::InvalidOptions);

    std::vector<ParametricCurveRequest> effectiveRequests = requests;
    std::vector<std::optional<ParametricPeriodReduction>> periodReductions;
    periodReductions.reserve(effectiveRequests.size());
    for (auto& request : effectiveRequests)
        periodReductions.push_back(reduceParametricCurvePeriod(
            request, builtins, mathematics, angles, assumptions));

    const std::size_t precisionBits = options.coarseSampling.precisionBits;
    std::vector<PlotAnalysis> xAnalyses;
    std::vector<PlotAnalysis> yAnalyses;
    std::vector<std::pair<PlotProgram, PlotProgram>> programs;
    std::vector<std::optional<ExactCurveGeometry>> exactGeometries;
    std::vector<std::optional<ExactCurveRange2D>> exactRanges;
    std::vector<SampledCurve2D> coarseCurves;
    std::vector<CurveRangeEstimate1D> xRanges;
    std::vector<CurveRangeEstimate1D> yRanges;
    xAnalyses.reserve(requests.size());
    yAnalyses.reserve(requests.size());
    programs.reserve(requests.size());
    exactGeometries.reserve(requests.size());
    exactRanges.reserve(requests.size());
    coarseCurves.reserve(requests.size());

    for (std::size_t i = 0; i < effectiveRequests.size(); ++i) {
        const auto& request = effectiveRequests[i];
        PlotRequest xr{request.xExpression, request.parameter, request.lower, request.upper, request.options};
        PlotRequest yr{request.yExpression, request.parameter, request.lower, request.upper, request.options};
        xAnalyses.push_back(analyzePlotRequest(xr, builtins, mathematics, angles, infinitySymbol, assumptions));
        yAnalyses.push_back(analyzePlotRequest(yr, builtins, mathematics, angles, infinitySymbol, assumptions));
        if (xAnalyses.back().samplingSafety != PlotSamplingSafety::Safe
            || yAnalyses.back().samplingSafety != PlotSamplingSafety::Safe)
            return parametricFailure(PlotPipelineStatus::UnsafeDiscontinuity, i);

        auto compileCoordinate = [&](expression::Expr expression, const PlotAnalysis& analysis) -> std::optional<PlotProgram> {
            if (analysis.domain.coverage == PlotDomainCoverage::Complete) {
                simplification::SimplificationContext context{builtins, mathematics, angles, assumptions};
                context.assumeExpressionsDefined = true;
                expression = simplification::Simplifier{}.simplify(expression, context);
            }
            const auto compiled = compilePlotProgram(expression, request.parameter, builtins, mathematics);
            if (!compiled)
                return std::nullopt;
            return *compiled.program;
        };
        auto xp = compileCoordinate(request.xExpression, xAnalyses.back());
        auto yp = compileCoordinate(request.yExpression, yAnalyses.back());
        if (!xp || !yp)
            return parametricFailure(PlotPipelineStatus::CompileFailed, i);
        programs.emplace_back(std::move(*xp), std::move(*yp));

        const auto domain = intersectParametricDomains(
            xAnalyses.back().domain, yAnalyses.back().domain, xr,
            builtins, mathematics, angles, precisionBits);
        if (!domain)
            return parametricFailure(PlotPipelineStatus::DomainIntersectionFailed, i);
        auto exactGeometry = recognizeExactParametricGeometry(
            request, builtins, mathematics, angles, assumptions);
        const bool exactGeometryApplicable = exactGeometry
            && domain->intervals.size() == 1
            && domain->intervals.front().lowerInclusion == PlotEndpointInclusion::Closed
            && domain->intervals.front().upperInclusion == PlotEndpointInclusion::Closed;
        const ExactCurveGeometry* exactGeometryForSampling =
            exactGeometryApplicable && !options.forceSampledGeometry
                ? &*exactGeometry : nullptr;
        auto coarseOptions = samplingDensityOptions(
            options.coarseSampling, request.options,
            options.forceSampledGeometry || !exactGeometryApplicable);
        coarseOptions.forcePolylineGeometry = options.forceSampledGeometry;
        const auto sampled = coarseSampleParametricPlot(
            request, *domain, programs.back().first, programs.back().second,
            builtins, mathematics, angles, coarseOptions, exactGeometryForSampling);
        std::optional<ExactCurveRange2D> exactRange;
        if (exactGeometryApplicable) {
            exactRange = estimateExactCurveRange(
                *exactGeometry, request, builtins, mathematics, angles, assumptions,
                precisionBits, options.rangeEstimator);
        }
        exactGeometries.push_back(std::move(exactGeometry));
        exactRanges.push_back(exactRange);
        if (!sampled)
            return parametricFailure(PlotPipelineStatus::SamplingFailed, i);
        coarseCurves.push_back(*sampled.curve);
        if (exactRange) {
            xRanges.push_back(exactRange->x);
            yRanges.push_back(exactRange->y);
        }
        else {
            const auto xrng = estimateAnalyzedCurveCoordinateRange(
                xAnalyses.back(), coarseCurves.back(), CurveCoordinate2D::X,
                options.rangeEstimator);
            const auto yrng = estimateAnalyzedCurveCoordinateRange(
                yAnalyses.back(), coarseCurves.back(), CurveCoordinate2D::Y,
                options.rangeEstimator);
            if (!xrng || !yrng)
                return parametricFailure(PlotPipelineStatus::RangeEstimationFailed, i);
            xRanges.push_back(*xrng);
            yRanges.push_back(*yrng);
        }
    }

    const auto xRange = combineCurveRangeEstimates(xRanges, precisionBits, options.rangeEstimator);
    const auto yRange = combineCurveRangeEstimates(yRanges, precisionBits, options.rangeEstimator);
    if (!xRange || !yRange)
        return parametricFailure(PlotPipelineStatus::RangeEstimationFailed);
    CurveViewRange2D range{*xRange, *yRange};

    PlotViewportMm viewport = options.viewport;
    if (requestOptions->aspectRatio) {
        const auto ratioValue = finiteConstantValue(
            *requestOptions->aspectRatio, optionRequests.front(), builtins, mathematics, angles,
            precisionBits);
        if (!ratioValue)
            return parametricFailure(PlotPipelineStatus::InvalidOptions);
        const double ratio = bigFloatToDouble(*ratioValue);
        if (!(ratio > 0.0) || !std::isfinite(ratio))
            return parametricFailure(PlotPipelineStatus::InvalidOptions);
        viewport.heightMm = viewport.widthMm * ratio;
    }

    const auto transform = makePlotViewTransform(range, precisionBits, viewport);
    if (!transform)
        return parametricFailure(PlotPipelineStatus::ViewTransformFailed);

    std::vector<SampledCurve2D> refinedCurves;
    refinedCurves.reserve(coarseCurves.size());
    for (std::size_t i = 0; i < coarseCurves.size(); ++i) {
        const bool densityApplies = std::any_of(
            coarseCurves[i].segments.begin(), coarseCurves[i].segments.end(),
            [](const auto& segment) {
                return segment.geometryKind == PlotSegmentGeometryKind::Polyline;
            });
        const auto adaptiveOptions = adaptiveDensityOptions(
            options.adaptiveSampling, effectiveRequests[i].options, densityApplies);
        const auto refined = refineParametricPlotSamples(
            programs[i].first, programs[i].second, coarseCurves[i], *transform,
            angles, adaptiveOptions);
        if (!refined)
            return parametricFailure(PlotPipelineStatus::AdaptiveSamplingFailed, i);
        refinedCurves.push_back(*refined.curve);
    }

    const auto plotScene = buildPlotScene(refinedCurves);
    if (!plotScene)
        return parametricFailure(PlotPipelineStatus::SceneBuildFailed);

    PlotAxisLayoutOptions axisOptions = options.axisLayout;
    PlotGraphicsLoweringOptions loweringOptions = options.graphicsLowering;
    if (requestOptions->ticks && !*requestOptions->ticks) {
        axisOptions.generateMajorTicks = false;
        loweringOptions.showMajorTicks = false;
        loweringOptions.showMajorTickLabels = false;
    }
    const auto axes = layoutPlotAxes(*transform, axisOptions);
    if (!axes)
        return parametricFailure(PlotPipelineStatus::AxisLayoutFailed);
    const auto graphicsScene = lowerPlotToGraphics(
        *plotScene.scene, *axes.layout, *transform, loweringOptions);
    if (!graphicsScene)
        return parametricFailure(PlotPipelineStatus::GraphicsLoweringFailed);

    ParametricPlotPipelineResult result;
    result.status = PlotPipelineStatus::Success;
    result.output = ParametricPlotPipelineOutput{
        std::move(xAnalyses), std::move(yAnalyses), std::move(programs),
        std::move(periodReductions), std::move(exactGeometries), std::move(exactRanges),
        std::move(refinedCurves), std::move(range),
        *transform, *axes.layout, *plotScene.scene, *graphicsScene.scene};
    return result;
}

ParametricPlotPipelineResult buildParametricPlotPipeline(
    const ParametricPlotRequestSet& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    return buildParametricPlotPipeline(
        splitParametricPlotRequests(request), builtins, mathematics, angles,
        infinitySymbol, assumptions, options);
}

ParametricPlotPipelineResult buildParametricPlotPipeline(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    return buildParametricPlotPipeline(
        std::vector<ParametricCurveRequest>{request}, builtins, mathematics, angles,
        infinitySymbol, assumptions, options);
}

ParametricPlotSvgPipelineResult renderParametricPlotSvg(
    const std::vector<ParametricCurveRequest>& requests,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    auto built = buildParametricPlotPipeline(
        requests, builtins, mathematics, angles, infinitySymbol, assumptions, options);
    ParametricPlotSvgPipelineResult result;
    result.status = built.status;
    result.failedCurveIndex = built.failedCurveIndex;
    if (!built)
        return result;
    auto svg = graphics::renderGraphics(built.output->graphicsScene, graphics::GraphicsFormat::Svg);
    if (!svg || !svg.data) {
        result.status = PlotPipelineStatus::SvgRenderFailed;
        return result;
    }
    result.status = PlotPipelineStatus::Success;
    result.svg = std::move(svg.data);
    result.output = std::move(built.output);
    return result;
}

ParametricPlotSvgPipelineResult renderParametricPlotSvg(
    const ParametricPlotRequestSet& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    return renderParametricPlotSvg(
        splitParametricPlotRequests(request), builtins, mathematics, angles,
        infinitySymbol, assumptions, options);
}

ParametricPlotSvgPipelineResult renderParametricPlotSvg(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const PlotPipelineOptions& options) {
    return renderParametricPlotSvg(
        std::vector<ParametricCurveRequest>{request}, builtins, mathematics, angles,
        infinitySymbol, assumptions, options);
}

} // namespace mmcal::plot
