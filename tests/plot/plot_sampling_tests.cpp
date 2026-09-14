// Plot粗samplingの区間・status契約
#include "plot_sampling_tests.hpp"
#include "../../version.h"

#include "kernel/kernel_session.hpp"
#include "mathematics/angle.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "numeric/rounding_mode.hpp"
#include "plot/plot_adaptive_sampling.hpp"
#include "plot/plot_axis_layout.hpp"
#include "plot/plot_anchors.hpp"
#include "plot/plot_intersections.hpp"
#include "plot/plot_analysis.hpp"
#include "plot/plot_program.hpp"
#include "plot/plot_graphics_lowering.hpp"
#include "plot/plot_exact_geometry.hpp"
#include "plot/plot_periodicity.hpp"
#include "plot/plot_pipeline.hpp"
#include "graphics/eps_backend.hpp"
#include "graphics/graphics_backend.hpp"
#include "graphics/pdf_backend.hpp"
#include "graphics/png_backend.hpp"
#include "graphics/svg_backend.hpp"
#include "graphics/webp_backend.hpp"
#include "plot/plot_sample_analysis.hpp"
#include "plot/plot_request.hpp"
#include "plot/plot_range_estimator.hpp"
#include "plot/plot_sampling.hpp"
#include "plot/plot_scene.hpp"
#include "plot/plot_view_transform.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string_view>
#include <type_traits>
#include <variant>

namespace mmcal::tests {
namespace {

using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] plot::PlotRequest request(
    kernel::KernelSession& session,
    std::string_view expression,
    std::string_view lower,
    std::string_view upper) {
    const auto x = session.evaluate("x").asSymbol();
    return plot::PlotRequest{
        session.evaluate(expression), x,
        session.evaluate(lower), session.evaluate(upper)};
}

[[nodiscard]] plot::PlotCompileResult compile(
    kernel::KernelSession& session,
    const plot::PlotRequest& request) {
    return plot::compilePlotProgram(
        request.expression, request.variable,
        session.builtinRegistry(), session.mathRegistry());
}

[[nodiscard]] plot::PlotAnalysis analyze(
    kernel::KernelSession& session,
    const plot::PlotRequest& request) {
    const auto* infinity = session.symbolRegistry().find("Infinity");
    if (!infinity)
        return plot::makeInitialPlotAnalysis(request);
    return plot::analyzePlotRequest(
        request,
        session.builtinRegistry(), session.mathRegistry(),
        mathematics::defaultAngleSemantics(), infinity->symbol);
}

[[nodiscard]] bool allFinite(const plot::SampledCurve& curve) {
    for (const auto& segment : curve.segments)
        for (const auto& sample : segment.samples)
            if (!sample.finite())
                return false;
    return true;
}




[[nodiscard]] bool hasSuspicion(
    const plot::PlotCoarseInspection& inspection,
    plot::PlotSuspicionKind kind) {
    for (const auto& span : inspection.suspiciousSpans)
        if (span.kind == kind)
            return true;
    return false;
}

[[nodiscard]] std::size_t countTagged(
    const plot::SampledCurve& curve,
    plot::PlotPointTag tag) {
    std::size_t count = 0;
    for (const auto& segment : curve.segments)
        for (const auto& sample : segment.samples)
            if (plot::hasPlotPointTag(sample.tags, tag))
                ++count;
    return count;
}


[[nodiscard]] bool curvePathCoordinatesStayNearClip(
    const graphics::GraphicsScene& scene,
    double toleranceMm = 0.011) {
    for (const auto& node : scene.nodes) {
        const auto* path = std::get_if<graphics::GraphicsPathNode>(&node);
        if (!path || !path->semantic || !path->clipRect
            || path->semantic->kind != graphics::GraphicsSemanticKind::Curve)
            continue;
        const auto& rect = *path->clipRect;
        const auto inside = [&](const graphics::GraphicsPointMm& point) {
            return point.xMm >= rect.xMm - toleranceMm
                && point.xMm <= rect.xMm + rect.widthMm + toleranceMm
                && point.yMm >= rect.yMm - toleranceMm
                && point.yMm <= rect.yMm + rect.heightMm + toleranceMm;
        };
        for (const auto& command : path->commands) {
            const bool valid = std::visit([&](const auto& value) {
                using T = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<T, graphics::GraphicsMoveTo>
                    || std::is_same_v<T, graphics::GraphicsLineTo>)
                    return inside(value.point);
                else if constexpr (std::is_same_v<T, graphics::GraphicsQuadraticTo>)
                    return inside(value.control) && inside(value.point);
                else if constexpr (std::is_same_v<T, graphics::GraphicsCubicTo>)
                    return inside(value.control1) && inside(value.control2)
                        && inside(value.point);
                else
                    return true;
            }, command);
            if (!valid)
                return false;
        }
    }
    return true;
}


[[nodiscard]] bool curvePathsHaveNoDanglingMoveTo(
    const graphics::GraphicsScene& scene) {
    for (const auto& node : scene.nodes) {
        const auto* path = std::get_if<graphics::GraphicsPathNode>(&node);
        if (!path || !path->semantic
            || path->semantic->kind != graphics::GraphicsSemanticKind::Curve)
            continue;
        bool pendingMove = false;
        for (const auto& command : path->commands) {
            if (std::holds_alternative<graphics::GraphicsMoveTo>(command)) {
                if (pendingMove)
                    return false;
                pendingMove = true;
            }
            else if (std::holds_alternative<graphics::GraphicsLineTo>(command)
                || std::holds_alternative<graphics::GraphicsQuadraticTo>(command)
                || std::holds_alternative<graphics::GraphicsCubicTo>(command)) {
                pendingMove = false;
            }
            else if (std::holds_alternative<graphics::GraphicsClosePath>(command)) {
                if (pendingMove)
                    return false;
                pendingMove = false;
            }
        }
        if (pendingMove)
            return false;
    }
    return true;
}


[[nodiscard]] bool clippedEndpointMarkersIntersectViewport(
    const graphics::GraphicsScene& scene) {
    for (const auto& node : scene.nodes) {
        const auto* circle = std::get_if<graphics::GraphicsCircleNode>(&node);
        if (!circle || !circle->semantic || !circle->clipRect
            || circle->semantic->kind != graphics::GraphicsSemanticKind::CurveEndpoint)
            continue;
        const auto& rect = *circle->clipRect;
        if (circle->center.xMm + circle->radiusMm < rect.xMm
            || circle->center.xMm - circle->radiusMm > rect.xMm + rect.widthMm
            || circle->center.yMm + circle->radiusMm < rect.yMm
            || circle->center.yMm - circle->radiusMm > rect.yMm + rect.heightMm)
            return false;
    }
    return true;
}

} // namespace

void runPlotSamplingTests(TestRunner& tests) {
    kernel::KernelSession session;
    const plot::CoarseSamplingOptions options{96, 5, 64};

    {
        const auto one = numeric::BigFloat::fromBigInt(
            BigInt{1}, 96, numeric::RoundingMode::NearestEven);
        const auto two = numeric::BigFloat::fromBigInt(
            BigInt{2}, 96, numeric::RoundingMode::NearestEven);
        const auto three = numeric::BigFloat::fromBigInt(
            BigInt{3}, 96, numeric::RoundingMode::NearestEven);
        const plot::CurveSample2D sample{
            one, two, three, plot::PlotNumericStatus::Finite};
        tests.expect(
            sample.parameter.toRational() == Rational{BigInt{1}}
                && sample.x.toRational() == Rational{BigInt{2}}
                && sample.y.toRational() == Rational{BigInt{3}},
            "CurveSample2D: parameter is independent from geometry x/y");
    }


    {
        const auto* infinity = session.symbolRegistry().find("Infinity");
        if (infinity) {
            plot::PlotPipelineOptions densityOptions;
            densityOptions.adaptiveSampling.chordToleranceMm = 1.0e6;
            densityOptions.adaptiveSampling.preserveAxisIntersections = false;
            densityOptions.adaptiveSampling.preserveLocalExtrema = false;

            auto defaultDensity = request(session, "sin[x]", "-Pi", "Pi");
            defaultDensity.options.plotPoints = 100;
            auto doubleDensity = defaultDensity;
            doubleDensity.options.plotPoints = 200;
            const auto defaultBuilt = plot::buildPlotPipeline(
                defaultDensity, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol, {}, densityOptions);
            const auto doubleBuilt = plot::buildPlotPipeline(
                doubleDensity, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol, {}, densityOptions);
            const auto sampleCount = [](const auto& curve) {
                std::size_t count = 0;
                for (const auto& segment : curve.segments)
                    count += segment.samples.size();
                return count;
            };
            tests.expect(
                defaultBuilt && doubleBuilt
                    && sampleCount(doubleBuilt.output->curves.front()) + 2
                        >= 2 * sampleCount(defaultBuilt.output->curves.front())
                    && sampleCount(doubleBuilt.output->curves.front())
                        <= 2 * sampleCount(defaultBuilt.output->curves.front()) + 2,
                "PlotPoints: 100 preserves the established coarse mesh and 200 doubles generic Plot density");

            const auto u = session.evaluate("u").asSymbol();
            plot::ParametricCurveRequest parametricDefault{
                session.evaluate("cos[7u]"), session.evaluate("sin[5u]"),
                u, session.evaluate("0"), session.evaluate("2Pi")};
            parametricDefault.options.plotPoints = 100;
            auto parametricDouble = parametricDefault;
            parametricDouble.options.plotPoints = 200;
            const auto parametricDefaultBuilt = plot::buildParametricPlotPipeline(
                parametricDefault, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol, {}, densityOptions);
            const auto parametricDoubleBuilt = plot::buildParametricPlotPipeline(
                parametricDouble, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol, {}, densityOptions);
            tests.expect(
                parametricDefaultBuilt && parametricDoubleBuilt
                    && sampleCount(parametricDoubleBuilt.output->curves.front()) + 2
                        >= 2 * sampleCount(parametricDefaultBuilt.output->curves.front())
                    && sampleCount(parametricDoubleBuilt.output->curves.front())
                        <= 2 * sampleCount(parametricDefaultBuilt.output->curves.front()) + 2,
                "PlotPoints: generic ParametricPlot uses the same relative coarse-density semantics");

            plot::ParametricCurveRequest exactCircle{
                session.evaluate("cos[u]"), session.evaluate("sin[u]"),
                u, session.evaluate("0"), session.evaluate("2Pi")};
            exactCircle.options.plotPoints = 100;
            auto exactCircleMaximum = exactCircle;
            exactCircleMaximum.options.plotPoints = 1024;
            const auto circle100 = plot::buildParametricPlotPipeline(
                exactCircle, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol, {}, densityOptions);
            const auto circle1024 = plot::buildParametricPlotPipeline(
                exactCircleMaximum, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol, {}, densityOptions);
            tests.expect(
                circle100 && circle1024
                    && sampleCount(circle100.output->curves.front())
                        == sampleCount(circle1024.output->curves.front())
                    && circle100.output->exactGeometries.front()
                    && circle1024.output->exactGeometries.front()
                    && std::holds_alternative<plot::ExactCircleGeometry>(
                        circle100.output->exactGeometries.front()->value)
                    && std::holds_alternative<plot::ExactCircleGeometry>(
                        circle1024.output->exactGeometries.front()->value),
                "PlotPoints: exact circle geometry ignores density increases and remains exact");
        }
    }

    {
        const auto u = session.evaluate("u").asSymbol();
        const auto spectrum = plot::analyzeHarmonicSpectrum(
            session.evaluate("sin[u]^2"), u,
            session.builtinRegistry(), session.mathRegistry(),
            mathematics::defaultAngleSemantics());
        tests.expect(
            spectrum && spectrum->frequencies
                == std::vector<Rational>{Rational{BigInt{0}}, Rational{BigInt{2}}},
            "HarmonicSpectrum: integer trig powers expose reduced finite frequency support");

        const plot::ParametricCurveRequest harmonicRequest{
            session.evaluate("cos[u]+cos[2u]/2"),
            session.evaluate("sin[u]+sin[4u]/4"),
            u, session.evaluate("0"), session.evaluate("4Pi")};
        const auto certificate = plot::analyzeParametricCurvePeriod(
            harmonicRequest,
            session.builtinRegistry(), session.mathRegistry(),
            mathematics::defaultAngleSemantics());
        tests.expect(
            certificate
                && certificate->turns == Rational{BigInt{1}}
                && certificate->proofKind == plot::PeriodProofKind::HarmonicSpectrum
                && certificate->spectrum.has_value(),
            "PeriodCertificate: mixed harmonics prove a 2Pi common period through spectrum analysis");

        const plot::ParametricCurveRequest squareRequest{
            session.evaluate("sin[u]^2"), session.evaluate("cos[u]^2"),
            u, session.evaluate("0"), session.evaluate("2Pi")};
        const auto squareCertificate = plot::analyzeParametricCurvePeriod(
            squareRequest,
            session.builtinRegistry(), session.mathRegistry(),
            mathematics::defaultAngleSemantics());
        tests.expect(
            squareCertificate
                && squareCertificate->turns == Rational{BigInt{1}, BigInt{2}}
                && squareCertificate->proofKind == plot::PeriodProofKind::HarmonicSpectrum,
            "PeriodCertificate: sin^2/cos^2 prove the shorter Pi period exactly");

        const plot::ParametricCurveRequest quadraticRequest{
            session.evaluate("u"), session.evaluate("u^2"),
            u, session.evaluate("-2"), session.evaluate("2")};
        const auto quadratic = plot::recognizeExactParametricGeometry(
            quadraticRequest,
            session.builtinRegistry(), session.mathRegistry(),
            mathematics::defaultAngleSemantics());
        const auto* q = quadratic
            ? std::get_if<plot::ExactQuadraticBezierGeometry>(&quadratic->value)
            : nullptr;
        tests.expect(
            q && q->start == plot::SymbolicPoint2D{session.evaluate("-2"), session.evaluate("4")}
                && q->control == plot::SymbolicPoint2D{session.evaluate("0"), session.evaluate("-4")}
                && q->end == plot::SymbolicPoint2D{session.evaluate("2"), session.evaluate("4")},
            "ExactCurveGeometry: quadratic parametric polynomials retain exact symbolic Bezier controls");

        const plot::ParametricCurveRequest cubicRequest{
            session.evaluate("u^3-3u"), session.evaluate("u^2-1"),
            u, session.evaluate("-2"), session.evaluate("2")};
        const auto cubic = plot::recognizeExactParametricGeometry(
            cubicRequest,
            session.builtinRegistry(), session.mathRegistry(),
            mathematics::defaultAngleSemantics());
        const auto* c = cubic
            ? std::get_if<plot::ExactCubicBezierGeometry>(&cubic->value)
            : nullptr;
        tests.expect(
            c && c->start == plot::SymbolicPoint2D{session.evaluate("-2"), session.evaluate("3")}
                && c->control1 == plot::SymbolicPoint2D{session.evaluate("10"), session.evaluate("-7/3")}
                && c->control2 == plot::SymbolicPoint2D{session.evaluate("-10"), session.evaluate("-7/3")}
                && c->end == plot::SymbolicPoint2D{session.evaluate("2"), session.evaluate("3")},
            "ExactCurveGeometry: cubic parametric polynomials retain exact symbolic Bezier controls");

        const plot::ParametricCurveRequest circleRequest{
            session.evaluate("cos[u]"), session.evaluate("sin[u]"),
            u, session.evaluate("0"), session.evaluate("2Pi")};
        const auto circle = plot::recognizeExactParametricGeometry(
            circleRequest,
            session.builtinRegistry(), session.mathRegistry(),
            mathematics::defaultAngleSemantics());
        tests.expect(
            circle && std::holds_alternative<plot::ExactCircleGeometry>(circle->value),
            "ExactCurveGeometry: a complete unit circle is classified symbolically as Circle");

        const plot::ParametricCurveRequest ellipseRequest{
            session.evaluate("3+2cos[u]+sin[u]"),
            session.evaluate("-1+4cos[u]-2sin[u]"),
            u, session.evaluate("0"), session.evaluate("2Pi")};
        const auto ellipse = plot::recognizeExactParametricGeometry(
            ellipseRequest,
            session.builtinRegistry(), session.mathRegistry(),
            mathematics::defaultAngleSemantics());
        tests.expect(
            ellipse && std::holds_alternative<plot::ExactEllipseGeometry>(ellipse->value),
            "ExactCurveGeometry: a non-orthogonal affine sin/cos image remains an exact Ellipse");

        const plot::ParametricCurveRequest arcRequest{
            session.evaluate("cos[u]"), session.evaluate("sin[u]"),
            u, session.evaluate("0"), session.evaluate("Pi")};
        const auto arc = plot::recognizeExactParametricGeometry(
            arcRequest,
            session.builtinRegistry(), session.mathRegistry(),
            mathematics::defaultAngleSemantics());
        const auto* ellipticArc = arc
            ? std::get_if<plot::ExactEllipticArcGeometry>(&arc->value)
            : nullptr;
        tests.expect(
            ellipticArc
                && ellipticArc->startTurns == Rational{BigInt{0}}
                && ellipticArc->sweepTurns == Rational{BigInt{1}, BigInt{2}},
            "ExactCurveGeometry: a semicircle retains exact EllipticArc phase in turns");

        const auto* infinity = session.symbolRegistry().find("Infinity");
        if (infinity) {
            const auto built = plot::buildParametricPlotPipeline(
                quadraticRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            const bool hasSymbolicQuadratic = built
                && built.output->exactGeometries.size() == 1
                && built.output->exactGeometries.front()
                && std::holds_alternative<plot::ExactQuadraticBezierGeometry>(
                    built.output->exactGeometries.front()->value);
            const bool symbolicControlDrivesSampling = built
                && built.output->curves.size() == 1
                && built.output->curves.front().segments.size() == 1
                && built.output->curves.front().segments.front().quadraticControlPoint
                && built.output->curves.front().segments.front().quadraticControlPoint->x.toRational()
                    == Rational{BigInt{0}}
                && built.output->curves.front().segments.front().quadraticControlPoint->y.toRational()
                    == Rational{BigInt{-4}};
            tests.expect(
                hasSymbolicQuadratic && symbolicControlDrivesSampling,
                "ParametricPlot pipeline: symbolic exact geometry drives the numeric Bezier lowering path");

            const auto arcBuilt = plot::buildParametricPlotPipeline(
                arcRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            const bool arcPrimitive = arcBuilt
                && arcBuilt.output->curves.size() == 1
                && arcBuilt.output->curves.front().segments.size() == 1
                && arcBuilt.output->curves.front().segments.front().geometryKind
                    == plot::PlotSegmentGeometryKind::EllipticArc
                && arcBuilt.output->curves.front().segments.front().ellipticArcGeometry.has_value()
                && std::any_of(
                    arcBuilt.output->graphicsScene.nodes.begin(),
                    arcBuilt.output->graphicsScene.nodes.end(),
                    [](const auto& node) {
                        return std::holds_alternative<graphics::GraphicsEllipticArcNode>(node);
                    });
            const bool arcRange = arcBuilt
                && arcBuilt.output->exactRanges.size() == 1
                && arcBuilt.output->exactRanges.front()
                && arcBuilt.output->exactRanges.front()->tight
                && arcBuilt.output->exactRanges.front()->x.dataMinimum.toRational()
                    == Rational{BigInt{-1}}
                && arcBuilt.output->exactRanges.front()->x.dataMaximum.toRational()
                    == Rational{BigInt{1}}
                && arcBuilt.output->exactRanges.front()->y.dataMinimum.toRational()
                    == Rational{BigInt{0}}
                && arcBuilt.output->exactRanges.front()->y.dataMaximum.toRational()
                    == Rational{BigInt{1}};
            tests.expect(
                arcPrimitive && arcRange,
                "ParametricPlot pipeline: exact EllipticArc drives primitive lowering and tight analytic range");

            const bool quadraticRange = built
                && built.output->exactRanges.size() == 1
                && built.output->exactRanges.front()
                && built.output->exactRanges.front()->tight
                && built.output->exactRanges.front()->x.dataMinimum.toRational()
                    == Rational{BigInt{-2}}
                && built.output->exactRanges.front()->x.dataMaximum.toRational()
                    == Rational{BigInt{2}}
                && built.output->exactRanges.front()->y.dataMinimum.toRational()
                    == Rational{BigInt{0}}
                && built.output->exactRanges.front()->y.dataMaximum.toRational()
                    == Rational{BigInt{4}};
            tests.expect(
                quadraticRange,
                "ParametricPlot pipeline: quadratic Bezier Automatic range uses analytic derivative extrema");

            const auto cubicBuilt = plot::buildParametricPlotPipeline(
                cubicRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            const bool cubicRange = cubicBuilt
                && cubicBuilt.output->exactRanges.size() == 1
                && cubicBuilt.output->exactRanges.front()
                && cubicBuilt.output->exactRanges.front()->tight
                && cubicBuilt.output->exactRanges.front()->x.dataMinimum.toRational()
                    == Rational{BigInt{-2}}
                && cubicBuilt.output->exactRanges.front()->x.dataMaximum.toRational()
                    == Rational{BigInt{2}}
                && cubicBuilt.output->exactRanges.front()->y.dataMinimum.toRational()
                    == Rational{BigInt{-1}}
                && cubicBuilt.output->exactRanges.front()->y.dataMaximum.toRational()
                    == Rational{BigInt{3}};
            tests.expect(
                cubicRange,
                "ParametricPlot pipeline: cubic Bezier Automatic range solves analytic derivative roots");

            const plot::ParametricCurveRequest obliqueArcRequest{
                session.evaluate("2cos[u]+sin[u]"),
                session.evaluate("4cos[u]-2sin[u]"),
                u, session.evaluate("0"), session.evaluate("Pi/3")};
            const auto obliqueArcBuilt = plot::buildParametricPlotPipeline(
                obliqueArcRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(
                obliqueArcBuilt
                    && obliqueArcBuilt.output->exactRanges.size() == 1
                    && obliqueArcBuilt.output->exactRanges.front()
                    && !obliqueArcBuilt.output->exactRanges.front()->tight,
                "ParametricPlot pipeline: uncertifiable oblique arc extrema use the safe analytic ellipse envelope");


            const plot::ParametricCurveRequest reciprocalCurveRequest{
                session.evaluate("u"), session.evaluate("1/u"),
                u, session.evaluate("-4"), session.evaluate("4")};
            const auto reciprocalCurve = plot::buildParametricPlotPipeline(
                reciprocalCurveRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(
                reciprocalCurve
                    && !reciprocalCurve.output->range.x.boundaryTrimmed
                    && reciprocalCurve.output->range.y.boundaryTrimmed
                    && reciprocalCurve.output->range.y.dataMinimum
                        > reciprocalCurve.output->range.y.observedMinimum
                    && reciprocalCurve.output->range.y.dataMaximum
                        < reciprocalCurve.output->range.y.observedMaximum,
                "ParametricPlot range: a pole in one coordinate trims only that coordinate tail");

            const plot::ParametricCurveRequest tangentCurveRequest{
                session.evaluate("tan[u]"), session.evaluate("u"),
                u, session.evaluate("-Pi"), session.evaluate("Pi")};
            const auto tangentCurve = plot::buildParametricPlotPipeline(
                tangentCurveRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(
                tangentCurve
                    && tangentCurve.output->range.x.boundaryTrimmed
                    && !tangentCurve.output->range.y.boundaryTrimmed
                    && tangentCurve.output->range.x.dataMinimum
                        > tangentCurve.output->range.x.observedMinimum
                    && tangentCurve.output->range.x.dataMaximum
                        < tangentCurve.output->range.x.observedMaximum,
                "ParametricPlot range: coordinate-wise asymptote trimming prevents unrelated axes from being flattened");
        }
    }

    const auto runCoreSampling = [&] {
        const auto polynomialRequest = request(session, "x^2", "-1", "1");
        const auto polynomialProgram = compile(session, polynomialRequest);
        tests.expect(static_cast<bool>(polynomialProgram),
            "PlotSampling: polynomial prerequisite compiles");
        if (polynomialProgram) {
            const auto sampled = plot::coarseSamplePlot(
                polynomialRequest, plot::makeInitialPlotAnalysis(polynomialRequest),
                *polynomialProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            tests.expect(sampled && sampled.curve->segments.size() == 1
                    && sampled.curve->segments[0].samples.size() == 5,
                "PlotSampling: closed interval receives the requested coarse samples");
            if (sampled && sampled.curve->segments[0].samples.size() == 5) {
                const auto& samples = sampled.curve->segments[0].samples;
                tests.expect(samples.front().x.toRational() == Rational{BigInt{-1}}
                        && samples.front().y.toRational() == Rational{BigInt{1}}
                        && samples[2].x.toRational().isZero()
                        && samples[2].y.toRational().isZero()
                        && samples.back().x.toRational() == Rational{BigInt{1}},
                    "PlotSampling: closed endpoints and midpoint are sampled exactly");
                tests.expect(std::all_of(
                        samples.begin(), samples.end(),
                        [](const plot::CurveSample2D& sample) {
                            return sample.parameter == sample.x;
                        }),
                    "PlotSampling: graph Plot lowers to the shared 2D curve IR as t -> (t,f(t))");

                const auto xRange = plot::estimateCurveCoordinateRange(
                    *sampled.curve, plot::CurveCoordinate2D::X);
                const auto yRange = plot::estimateCurveCoordinateRange(
                    *sampled.curve, plot::CurveCoordinate2D::Y);
                tests.expect(xRange && yRange
                        && xRange->dataMinimum.toRational() == Rational{BigInt{-1}}
                        && xRange->dataMaximum.toRational() == Rational{BigInt{1}}
                        && yRange->dataMinimum.toRational().isZero()
                        && yRange->dataMaximum.toRational() == Rational{BigInt{1}},
                    "PlotRange: shared 1D estimator handles x/y coordinates independently");
                const auto transform = xRange && yRange
                    ? plot::makePlotViewTransform(
                        plot::CurveViewRange2D{*xRange, *yRange}, 96)
                    : std::nullopt;
                tests.expect(static_cast<bool>(transform),
                    "PlotViewTransform: shared 2D view range constructs a transform directly");
            }
        }
    
        const auto affineRequest = request(session, "2x+3", "-1", "1");
        const auto affineProgram = compile(session, affineRequest);
        if (affineProgram) {
            const auto sampled = plot::coarseSamplePlot(
                affineRequest, plot::makeInitialPlotAnalysis(affineRequest),
                *affineProgram.program, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            tests.expect(sampled && sampled.curve->segments.size() == 1
                    && sampled.curve->segments[0].samples.size() == 2
                    && sampled.curve->segments[0].geometryKind
                        == plot::PlotSegmentGeometryKind::StraightLine,
                "PlotSampling: affine functions lower directly to a two-endpoint straight segment");
            if (sampled) {
                const auto scene = plot::buildPlotScene(std::vector<plot::SampledCurve>{*sampled.curve});
                tests.expect(scene && scene.scene->curves.size() == 1
                        && scene.scene->curves[0].segments.size() == 1
                        && scene.scene->curves[0].segments[0].geometryKind
                            == plot::PlotSegmentGeometryKind::StraightLine
                        && scene.scene->curves[0].segments[0].vertices.size() == 2,
                    "PlotScene: affine geometry specialization survives semantic lowering");
            }
        }

        const auto cubicRequest = request(session, "x^3-2x", "-2", "2");
        const auto cubicProgram = compile(session, cubicRequest);
        if (cubicProgram) {
            const auto sampled = plot::coarseSamplePlot(
                cubicRequest, plot::makeInitialPlotAnalysis(cubicRequest),
                *cubicProgram.program, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            tests.expect(sampled && sampled.curve->segments.size() == 1
                    && sampled.curve->segments[0].geometryKind
                        == plot::PlotSegmentGeometryKind::CubicBezier
                    && sampled.curve->segments[0].samples.size() == 5
                    && sampled.curve->segments[0].bezierControlPoints.has_value(),
                "PlotSampling: cubic polynomials retain coarse analysis samples and exact cubic-Bezier control points");
            if (sampled) {
                const auto scene = plot::buildPlotScene(std::vector<plot::SampledCurve>{*sampled.curve});
                const auto inspection = plot::inspectCoarseSamples(
                    plot::makeInitialPlotAnalysis(cubicRequest), *sampled.curve);
                const auto range = plot::estimatePlotRange(
                    plot::makeInitialPlotAnalysis(cubicRequest), *sampled.curve, inspection);
                const auto transform = range ? plot::makePlotViewTransform(
                    cubicRequest, *range,
                    session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), 96) : std::nullopt;
                const auto axes = transform
                    ? plot::layoutPlotAxes(*transform)
                    : plot::PlotAxisLayoutResult{};
                const auto graphics = scene && transform && axes
                    ? plot::lowerPlotToGraphics(*scene.scene, *axes.layout, *transform)
                    : plot::PlotGraphicsLoweringResult{};
                bool hasCubicCommand = false;
                if (graphics) {
                    for (const auto& node : graphics.scene->nodes) {
                        const auto* path = std::get_if<graphics::GraphicsPathNode>(&node);
                        if (!path || !path->semantic
                            || path->semantic->kind != graphics::GraphicsSemanticKind::Curve)
                            continue;
                        for (const auto& command : path->commands)
                            hasCubicCommand = hasCubicCommand
                                || std::holds_alternative<graphics::GraphicsCubicTo>(command);
                    }
                }
                tests.expect(scene
                        && scene.scene->curves[0].segments[0].geometryKind
                            == plot::PlotSegmentGeometryKind::CubicBezier
                        && scene.scene->curves[0].segments[0].bezierControlPoints.has_value()
                        && hasCubicCommand,
                    "PlotGraphics: cubic-polynomial segments lower to SVG-ready cubic path commands");
            }
        }

        const auto quadraticRequest = request(session, "x^2-x-1", "-2", "2");
        const auto quadraticProgram = compile(session, quadraticRequest);
        if (quadraticProgram) {
            const auto sampled = plot::coarseSamplePlot(
                quadraticRequest, plot::makeInitialPlotAnalysis(quadraticRequest),
                *quadraticProgram.program, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            bool hasQuadraticCommand = false;
            if (sampled) {
                const auto scene = plot::buildPlotScene(std::vector<plot::SampledCurve>{*sampled.curve});
                const auto inspection = plot::inspectCoarseSamples(
                    plot::makeInitialPlotAnalysis(quadraticRequest), *sampled.curve);
                const auto range = plot::estimatePlotRange(
                    plot::makeInitialPlotAnalysis(quadraticRequest), *sampled.curve, inspection);
                const auto transform = range ? plot::makePlotViewTransform(
                    quadraticRequest, *range,
                    session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), 96) : std::nullopt;
                const auto axes = transform
                    ? plot::layoutPlotAxes(*transform)
                    : plot::PlotAxisLayoutResult{};
                const auto graphics = scene && transform && axes
                    ? plot::lowerPlotToGraphics(*scene.scene, *axes.layout, *transform)
                    : plot::PlotGraphicsLoweringResult{};
                if (graphics) {
                    for (const auto& node : graphics.scene->nodes) {
                        const auto* path = std::get_if<graphics::GraphicsPathNode>(&node);
                        if (!path || !path->semantic
                            || path->semantic->kind != graphics::GraphicsSemanticKind::Curve)
                            continue;
                        for (const auto& command : path->commands)
                            hasQuadraticCommand = hasQuadraticCommand
                                || std::holds_alternative<graphics::GraphicsQuadraticTo>(command);
                    }
                }
            }
            tests.expect(sampled
                    && sampled.curve->segments[0].geometryKind
                        == plot::PlotSegmentGeometryKind::QuadraticBezier
                    && sampled.curve->segments[0].quadraticControlPoint.has_value()
                    && hasQuadraticCommand,
                "PlotGraphics: quadratic polynomials lower exactly to quadratic Bezier path commands");
        }
    
        const auto periodicStepRequest = request(session, "floor[sin[x]]", "-Pi", "Pi");
        const auto periodicStepProgram = compile(session, periodicStepRequest);
        if (periodicStepProgram) {
            const auto sampled = plot::coarseSamplePlot(
                periodicStepRequest, analyze(session, periodicStepRequest),
                *periodicStepProgram.program, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            bool horizontalPieces = sampled && sampled.curve->segments.size() == 5;
            if (horizontalPieces) {
                for (const auto& segment : sampled.curve->segments) {
                    if (segment.samples.size() == 1)
                        continue;
                    if (segment.geometryKind != plot::PlotSegmentGeometryKind::PiecewiseConstant
                        || segment.samples.size() != 2
                        || !segment.samples[0].finite() || !segment.samples[1].finite()
                        || segment.samples[0].y != segment.samples[1].y) {
                        horizontalPieces = false;
                        break;
                    }
                }
            }
            tests.expect(horizontalPieces,
                "PlotSampling: periodic step pieces sample strictly inside boundaries and remain horizontal");
        }

        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto rationalPower = plot::buildPlotPipeline(
                request(session, "x^(3/2)", "-3", "3"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(rationalPower
                    && rationalPower.output->analyses.front().domain.coverage
                        == plot::PlotDomainCoverage::Complete
                    && rationalPower.output->analyses.front().domain.intervals.size() == 1,
                "PlotPipeline: noninteger rational Power reaches sampling only on its proven real branch");

            const auto positiveSelfPower = plot::buildPlotPipeline(
                request(session, "x^x", "0", "3"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            const auto mixedSelfPower = plot::buildPlotPipeline(
                request(session, "x^x", "-3", "3"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(positiveSelfPower
                    && mixedSelfPower.status == plot::PlotPipelineStatus::UnsafeDiscontinuity,
                "PlotPipeline: x^x uses request-bound proof on the positive branch but rejects the mixed real/complex branch safely");

            bool positiveSelfPowerTouchesYAxisHole = positiveSelfPower
                && positiveSelfPower.output->curves.size() == 1
                && !positiveSelfPower.output->curves.front().segments.empty();
            std::size_t positiveSelfPowerOpenMarkers = 0;
            if (positiveSelfPowerTouchesYAxisHole) {
                const auto& output = *positiveSelfPower.output;
                constexpr double toleranceMm = 0.0051;
                for (const auto& segment : output.curves.front().segments) {
                    if (segment.sourceInterval.lowerInclusion == plot::PlotEndpointInclusion::Open
                        && segment.lowerBoundaryParameter && !segment.samples.empty()) {
                        const auto boundary = output.transform.mapX(*segment.lowerBoundaryParameter);
                        const auto sample = output.transform.mapX(segment.samples.front().x);
                        positiveSelfPowerTouchesYAxisHole = boundary && sample
                            && std::abs(*sample - *boundary) <= toleranceMm;
                    }
                    if (segment.lowerEndpointPoint)
                        ++positiveSelfPowerOpenMarkers;
                }
            }
            tests.expect(positiveSelfPowerTouchesYAxisHole && positiveSelfPowerOpenMarkers == 1,
                "PlotPipeline: x^x on a nonnegative request approaches the x=0 hole physically and preserves its finite limit marker");

            const auto dynamicIntegerPower = plot::buildPlotPipeline(
                request(session, "(-2)^floor[x]", "-3", "3"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            bool dynamicStepsDisconnected = dynamicIntegerPower
                && dynamicIntegerPower.output->curves.size() == 1
                && dynamicIntegerPower.output->curves.front().segments.size() == 7;
            if (dynamicStepsDisconnected) {
                for (const auto& segment : dynamicIntegerPower.output->curves.front().segments) {
                    if (segment.samples.size() <= 1)
                        continue;
                    const auto y = segment.samples.front().y;
                    for (const auto& sample : segment.samples) {
                        if (!sample.finite() || sample.y != y) {
                            dynamicStepsDisconnected = false;
                            break;
                        }
                    }
                    if (!dynamicStepsDisconnected)
                        break;
                }
            }
            tests.expect(dynamicStepsDisconnected,
                "PlotPipeline: integer-valued Power exponents preserve nested step jumps without vertical connector artifacts");

            const auto reciprocalRootTan = plot::buildPlotPipeline(
                request(session, "tan[x]^(-1/2)", "-Pi", "Pi"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            bool reciprocalRootTouchesFinitePoleLimits = reciprocalRootTan
                && reciprocalRootTan.output->curves.size() == 1;
            std::size_t zeroLimitEndpoints = 0;
            bool reciprocalRootRangeContainsZero = static_cast<bool>(reciprocalRootTan);
            if (reciprocalRootTouchesFinitePoleLimits) {
                const auto& output = *reciprocalRootTan.output;
                const auto zero = numeric::BigFloat::fromBigInt(
                    BigInt{0}, output.transform.precisionBits(), numeric::RoundingMode::NearestEven);
                reciprocalRootRangeContainsZero = output.range.dataMinimum <= zero
                    && zero <= output.range.dataMaximum
                    && output.range.viewMinimum <= zero
                    && zero <= output.range.viewMaximum;
                constexpr double toleranceMm = 0.0051;
                for (const auto& segment : output.curves.front().segments) {
                    if (segment.sourceInterval.lowerInclusion == plot::PlotEndpointInclusion::Open
                        && segment.lowerBoundaryParameter && !segment.samples.empty()) {
                        const auto boundary = output.transform.mapX(*segment.lowerBoundaryParameter);
                        const auto sample = output.transform.mapX(segment.samples.front().x);
                        reciprocalRootTouchesFinitePoleLimits = boundary && sample
                            && std::abs(*sample - *boundary) <= toleranceMm;
                    }
                    if (segment.sourceInterval.upperInclusion == plot::PlotEndpointInclusion::Open
                        && segment.upperBoundaryParameter && !segment.samples.empty()) {
                        const auto boundary = output.transform.mapX(*segment.upperBoundaryParameter);
                        const auto sample = output.transform.mapX(segment.samples.back().x);
                        reciprocalRootTouchesFinitePoleLimits = reciprocalRootTouchesFinitePoleLimits
                            && boundary && sample
                            && std::abs(*sample - *boundary) <= toleranceMm;
                    }
                    const auto countZero = [&](const std::optional<plot::PlotEndpointPoint>& point) {
                        if (point && point->y.isZero())
                            ++zeroLimitEndpoints;
                    };
                    countZero(segment.lowerEndpointPoint);
                    countZero(segment.upperEndpointPoint);
                }
            }
            tests.expect(reciprocalRootTouchesFinitePoleLimits
                    && zeroLimitEndpoints == 2
                    && reciprocalRootRangeContainsZero,
                "PlotPipeline: tan[x]^(-1/2) approaches every open boundary physically and keeps zero-limit holes inside the automatic range");

            const auto absoluteTanPower = plot::buildPlotPipeline(
                request(session, "abs[tan[x]]^(3/2)", "-Pi", "Pi"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            bool absoluteTanPolesAreSymmetricAsymptotes = absoluteTanPower
                && absoluteTanPower.output->range.boundaryTrimmed
                && absoluteTanPower.output->curves.size() == 1;
            std::size_t poleApproachesBeyondView = 0;
            if (absoluteTanPolesAreSymmetricAsymptotes) {
                const auto& output = *absoluteTanPower.output;
                for (const auto& segment : output.curves.front().segments) {
                    const auto inspect = [&](bool lowerEndpoint) {
                        const auto inclusion = lowerEndpoint
                            ? segment.sourceInterval.lowerInclusion
                            : segment.sourceInterval.upperInclusion;
                        const auto& boundaryX = lowerEndpoint
                            ? segment.lowerBoundaryParameter : segment.upperBoundaryParameter;
                        if (inclusion != plot::PlotEndpointInclusion::Open
                            || !boundaryX || segment.samples.empty())
                            return;
                        const auto& sample = lowerEndpoint
                            ? segment.samples.front() : segment.samples.back();
                        const auto boundary = output.transform.mapX(*boundaryX);
                        const auto sampleX = output.transform.mapX(sample.x);
                        const auto sampleY = output.transform.mapY(sample.y);
                        if (boundary && sampleX && sampleY
                            && std::abs(*sampleX - *boundary) <= 0.0051
                            && *sampleY > output.transform.viewport().heightMm)
                            ++poleApproachesBeyondView;
                    };
                    inspect(true);
                    inspect(false);
                }
            }
            tests.expect(absoluteTanPolesAreSymmetricAsymptotes
                    && poleApproachesBeyondView == 4,
                "PlotPipeline: abs[tan[x]]^(3/2) classifies both tan poles as asymptotes and approaches all four sides with the same physical inset");


            const auto tangentClip = plot::buildPlotPipeline(
                request(session, "tan[x]", "-Pi", "Pi"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            const auto reciprocalClip = plot::buildPlotPipeline(
                request(session, "1/x", "-4", "4"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(
                tangentClip && reciprocalClip && absoluteTanPower
                    && curvePathCoordinatesStayNearClip(tangentClip.output->graphicsScene)
                    && curvePathCoordinatesStayNearClip(reciprocalClip.output->graphicsScene)
                    && curvePathCoordinatesStayNearClip(absoluteTanPower.output->graphicsScene),
                "PlotGraphics: pole/asymptote paths are geometrically clipped before backend serialization");

            auto explicitCubicRequest = request(session, "x^3", "-10", "10");
            explicitCubicRequest.options.plotRange = plot::PlotRangeOption{
                session.evaluate("-1"), session.evaluate("1")};
            const auto explicitCubic = plot::buildPlotPipeline(
                explicitCubicRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(
                explicitCubic
                    && curvePathCoordinatesStayNearClip(explicitCubic.output->graphicsScene),
                "PlotGraphics: exact Bezier paths with explicit PlotRange discard huge off-screen control geometry");


            auto farHoleRequest = request(
                session, "1000*((x^2-1)/(x-1))", "0", "2");
            farHoleRequest.options.plotRange = plot::PlotRangeOption{
                session.evaluate("-1"), session.evaluate("1")};
            const auto farHole = plot::buildPlotPipeline(
                farHoleRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(
                farHole
                    && clippedEndpointMarkersIntersectViewport(farHole.output->graphicsScene)
                    && std::none_of(
                        farHole.output->graphicsScene.nodes.begin(),
                        farHole.output->graphicsScene.nodes.end(),
                        [](const auto& node) {
                            const auto* circle = std::get_if<graphics::GraphicsCircleNode>(&node);
                            return circle && circle->semantic
                                && circle->semantic->kind
                                    == graphics::GraphicsSemanticKind::CurveEndpoint;
                        }),
                "PlotGraphics: fully off-screen hole markers are culled before serialization");

            const auto mixedClip = plot::buildPlotPipeline(
                std::vector<plot::PlotRequest>{
                    request(session, "tan[x]", "-Pi", "Pi"),
                    request(session, "1/x", "-Pi", "Pi"),
                    request(session, "sin[x]", "-Pi", "Pi")},
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(
                mixedClip && curvePathCoordinatesStayNearClip(mixedClip.output->graphicsScene),
                "PlotGraphics: multi-curve clipping remains bounded across poles and regular curves");

            const auto extremeReciprocalPower = plot::buildPlotPipeline(
                request(session, "(1/x)^100", "-1", "1"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(
                extremeReciprocalPower
                    && curvePathCoordinatesStayNearClip(
                        extremeReciprocalPower.output->graphicsScene)
                    && curvePathsHaveNoDanglingMoveTo(
                        extremeReciprocalPower.output->graphicsScene),
                "PlotGraphics: extreme pole overflow is clipped without dangling MoveTo commands");

            const auto cubicWithLine = plot::buildPlotPipeline(
                std::vector<plot::PlotRequest>{
                    request(session, "x^3-x", "-2", "2"),
                    request(session, "0", "-2", "2")},
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            bool cubicAnchorsSplitBezier = cubicWithLine
                && cubicWithLine.output->plotScene.curves.size() == 2
                && !cubicWithLine.output->plotScene.curves[0].segments.empty();
            std::size_t intersectionAnchors = 0;
            std::size_t extremumAnchors = 0;
            std::size_t internalAnchoredVertices = 0;
            std::size_t cubicCommands = 0;
            if (cubicAnchorsSplitBezier) {
                const auto& output = *cubicWithLine.output;
                for (const auto& anchor : output.plotScene.anchors) {
                    if (plot::hasPlotAnchorKind(anchor.kinds, plot::PlotAnchorKind::CurveIntersection))
                        ++intersectionAnchors;
                    if (plot::hasPlotAnchorKind(anchor.kinds, plot::PlotAnchorKind::LocalExtremum))
                        ++extremumAnchors;
                }
                const auto& segment = output.plotScene.curves[0].segments[0];
                for (std::size_t i = 1; i + 1 < segment.vertices.size(); ++i)
                    internalAnchoredVertices += segment.vertices[i].anchorId ? 1U : 0U;
                for (const auto& node : output.graphicsScene.nodes) {
                    const auto* path = std::get_if<graphics::GraphicsPathNode>(&node);
                    if (!path || !path->semantic
                        || path->semantic->kind != graphics::GraphicsSemanticKind::Curve
                        || path->semantic->id != output.plotScene.curves[0].id)
                        continue;
                    for (const auto& command : path->commands)
                        cubicCommands += std::holds_alternative<graphics::GraphicsCubicTo>(command)
                            ? 1U : 0U;
                }
            }
            tests.expect(cubicAnchorsSplitBezier
                    && intersectionAnchors >= 3
                    && extremumAnchors >= 2
                    && internalAnchoredVertices > 0
                    && cubicCommands == internalAnchoredVertices + 1,
                "PlotGraphics: cubic Beziers are de Casteljau-split at intersection/extremum semantic anchors without reverting to polylines");
        }

        const auto rationalRequest = request(session, "1/(x^2-1)", "-3", "3");
        const auto rationalProgram = compile(session, rationalRequest);
        if (rationalProgram) {
            const auto sampled = plot::coarseSamplePlot(
                rationalRequest, analyze(session, rationalRequest), *rationalProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            tests.expect(sampled && sampled.curve->segments.size() == 3 && allFinite(*sampled.curve),
                "PlotSampling: symbolic pole partition prevents coarse samples from hitting poles");
            if (sampled) {
                bool touchedPole = false;
                for (const auto& segment : sampled.curve->segments)
                    for (const auto& sample : segment.samples)
                        touchedPole = touchedPole
                            || sample.x.toRational() == Rational{BigInt{-1}}
                            || sample.x.toRational() == Rational{BigInt{1}};
                tests.expect(!touchedPole,
                    "PlotSampling: open pole endpoints are never evaluated");
            }
        }
    
        const auto logarithmRequest = request(session, "log[x]", "-3", "3");
        const auto logarithmProgram = compile(session, logarithmRequest);
        if (logarithmProgram) {
            const auto sampled = plot::coarseSamplePlot(
                logarithmRequest, analyze(session, logarithmRequest), *logarithmProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            tests.expect(sampled && sampled.curve->segments.size() == 1
                    && allFinite(*sampled.curve)
                    && sampled.curve->segments[0].samples.front().x.isPositive(),
                "PlotSampling: open log branch boundary is skipped before numeric evaluation");
        }

        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto si = plot::buildPlotPipeline(
                request(session, "Si[x]", "-3", "3"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            tests.expect(si
                    && si.output->analyses.front().domain.coverage == plot::PlotDomainCoverage::Complete
                    && si.output->analyses.front().domain.intervals.size() == 1,
                "PlotPipeline: Si uses its entire real domain without artificial cuts");

            const auto ei = plot::buildPlotPipeline(
                request(session, "Ei[x]", "-3", "3"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            bool eiZeroAsymptote = false;
            if (ei) {
                for (const auto& landmark : ei.output->analyses.front().landmarks) {
                    if (landmark.kind == plot::PlotLandmarkKind::VerticalAsymptote
                        && landmark.position == session.evaluate("0")) {
                        eiZeroAsymptote = true;
                        break;
                    }
                }
            }
            tests.expect(ei
                    && ei.output->analyses.front().domain.coverage == plot::PlotDomainCoverage::Complete
                    && ei.output->analyses.front().domain.intervals.size() == 2
                    && eiZeroAsymptote,
                "PlotPipeline: Ei keeps both real half-axes and treats zero as a vertical asymptote");

            const auto ci = plot::buildPlotPipeline(
                request(session, "Ci[x]", "-3", "3"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            bool ciZeroAsymptote = false;
            if (ci) {
                for (const auto& landmark : ci.output->analyses.front().landmarks) {
                    if (landmark.kind == plot::PlotLandmarkKind::VerticalAsymptote
                        && landmark.position == session.evaluate("0")) {
                        ciZeroAsymptote = true;
                        break;
                    }
                }
            }
            tests.expect(ci
                    && ci.output->analyses.front().domain.coverage == plot::PlotDomainCoverage::Complete
                    && ci.output->analyses.front().domain.intervals.size() == 1
                    && ci.output->analyses.front().domain.intervals.front().lower
                        == session.evaluate("0")
                    && ci.output->analyses.front().domain.intervals.front().lowerInclusion
                        == plot::PlotEndpointInclusion::Open
                    && ciZeroAsymptote,
                "PlotPipeline: principal Ci restricts real Plot to x>0 and marks the zero asymptote");
        }
    
        const auto unknownRequest = request(session, "1/x", "-1", "1");
        const auto unknownProgram = compile(session, unknownRequest);
        if (unknownProgram) {
            const auto sampled = plot::coarseSamplePlot(
                unknownRequest, plot::makeInitialPlotAnalysis(unknownRequest),
                *unknownProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            tests.expect(sampled && sampled.curve->segments.size() == 1
                    && sampled.curve->segments[0].samples[2].status
                        == plot::PlotNumericStatus::DivisionByZero,
                "PlotSampling: unknown-domain coarse pass preserves nonfinite sample status");
            if (sampled) {
                const auto scene = plot::buildPlotScene(std::vector<plot::SampledCurve>{*sampled.curve});
                tests.expect(scene && scene.scene->curves.size() == 1
                        && scene.scene->curves[0].segments.size() == 2
                        && scene.scene->curves[0].segments[0].vertices.size() == 2
                        && scene.scene->curves[0].segments[1].vertices.size() == 2,
                    "PlotScene: nonfinite numeric holes split paths instead of reconnecting across them");
            }
        }
    
        const auto tangentRequest = request(session, "tan[x]", "-Pi", "Pi");
        const auto tangentProgram = compile(session, tangentRequest);
        if (tangentProgram) {
            const auto sampled = plot::coarseSamplePlot(
                tangentRequest, analyze(session, tangentRequest), *tangentProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            tests.expect(sampled && sampled.curve->segments.size() == 3 && allFinite(*sampled.curve),
                "PlotSampling: exact Pi endpoints lower successfully and tan poles stay excluded");
        }
    
        const auto oscillatoryRequest = request(session, "sin[1/x]", "-1", "1");
        const auto oscillatoryProgram = compile(session, oscillatoryRequest);
        if (oscillatoryProgram) {
            const auto oscillatoryAnalysis = analyze(session, oscillatoryRequest);
            plot::CoarseSamplingOptions oscillatoryOptions{96, 33, 128};
            const auto sampled = plot::coarseSamplePlot(
                oscillatoryRequest, oscillatoryAnalysis, *oscillatoryProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), oscillatoryOptions);
            tests.expect(sampled && sampled.curve->segments.size() == 2 && allFinite(*sampled.curve),
                "PlotSampling: sin(1/x) never samples its excluded origin");
            if (sampled) {
                const auto inspection = plot::inspectCoarseSamples(
                    oscillatoryAnalysis, *sampled.curve);
                tests.expect(inspection.finiteRange.has_value()
                        && inspection.finiteRange->minimum >= numeric::BigFloat::fromBigInt(
                            BigInt{-1}, 96, numeric::RoundingMode::NearestEven)
                        && inspection.finiteRange->maximum <= numeric::BigFloat::fromBigInt(
                            BigInt{1}, 96, numeric::RoundingMode::NearestEven),
                    "PlotSampling: sin(1/x) coarse finite range remains inside [-1,1]");
                tests.expect(hasSuspicion(inspection, plot::PlotSuspicionKind::Oscillatory),
                    "PlotSampling: sin(1/x) coarse pass marks repeated turning as oscillatory");
                tests.expect(hasSuspicion(inspection, plot::PlotSuspicionKind::UnresolvedBoundary),
                    "PlotSampling: sin(1/x) punctured origin remains an unresolved boundary");
                const auto range = plot::estimatePlotRange(
                    oscillatoryAnalysis, *sampled.curve, inspection);
                tests.expect(range.has_value() && !range->boundaryTrimmed
                        && range->dataMinimum == inspection.finiteRange->minimum
                        && range->dataMaximum == inspection.finiteRange->maximum,
                    "PlotRange: unresolved oscillatory boundary does not discard bounded data");
            }
        }
    
        const auto constantRequest = request(session, "2", "-1", "1");
        const auto constantProgram = compile(session, constantRequest);
        if (constantProgram) {
            const auto sampled = plot::coarseSamplePlot(
                constantRequest, plot::makeInitialPlotAnalysis(constantRequest),
                *constantProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            if (sampled) {
                const auto inspection = plot::inspectCoarseSamples(
                    plot::makeInitialPlotAnalysis(constantRequest), *sampled.curve);
                const auto range = plot::estimatePlotRange(
                    plot::makeInitialPlotAnalysis(constantRequest), *sampled.curve, inspection);
                tests.expect(range.has_value() && range->flatExpanded
                        && range->dataMinimum == range->dataMaximum
                        && range->viewMinimum < range->dataMinimum
                        && range->viewMaximum > range->dataMaximum,
                    "PlotRange: flat data receives a finite symmetric viewport margin");
            }
        }
    
        const auto expRequest = request(session, "exp[5x]", "-1", "1");
        const auto expProgram = compile(session, expRequest);
        if (expProgram) {
            const auto expAnalysis = plot::makeInitialPlotAnalysis(expRequest);
            plot::CoarseSamplingOptions denseOptions{96, 33, 128};
            const auto sampled = plot::coarseSamplePlot(
                expRequest, expAnalysis, *expProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), denseOptions);
            if (sampled) {
                const auto inspection = plot::inspectCoarseSamples(expAnalysis, *sampled.curve);
                const auto range = plot::estimatePlotRange(expAnalysis, *sampled.curve, inspection);
                tests.expect(range.has_value() && !range->boundaryTrimmed
                        && range->dataMinimum == inspection.finiteRange->minimum
                        && range->dataMaximum == inspection.finiteRange->maximum,
                    "PlotRange: legitimate exponential growth is not percentile-trimmed");
            }
        }
    
        const auto poleRangeRequest = request(session, "1/(x^2-1)", "-3", "3");
        const auto poleRangeProgram = compile(session, poleRangeRequest);
        if (poleRangeProgram) {
            const auto poleAnalysis = analyze(session, poleRangeRequest);
            plot::CoarseSamplingOptions denseOptions{96, 33, 256};
            const auto sampled = plot::coarseSamplePlot(
                poleRangeRequest, poleAnalysis, *poleRangeProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), denseOptions);
            if (sampled) {
                const auto inspection = plot::inspectCoarseSamples(poleAnalysis, *sampled.curve);
                const auto range = plot::estimatePlotRange(poleAnalysis, *sampled.curve, inspection);
                tests.expect(range.has_value() && range->boundaryTrimmed
                        && range->excludedFiniteSamples > 0
                        && (range->dataMinimum > range->observedMinimum
                            || range->dataMaximum < range->observedMaximum),
                    "PlotRange: proven asymptote-adjacent samples do not dominate Automatic range");
            }
        }
    
    
        if (polynomialProgram) {
            const auto analysis = plot::makeInitialPlotAnalysis(polynomialRequest);
            plot::CoarseSamplingOptions coarseOptions{96, 5, 128};
            const auto sampled = plot::coarseSamplePlot(
                polynomialRequest, analysis, *polynomialProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), coarseOptions);
            if (sampled) {
                const auto inspection = plot::inspectCoarseSamples(analysis, *sampled.curve);
                const auto range = plot::estimatePlotRange(analysis, *sampled.curve, inspection);
                const auto transform = range ? plot::makePlotViewTransform(
                    polynomialRequest, *range,
                    session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), 96,
                    plot::PlotViewportMm{160.0, 100.0}) : std::nullopt;
                tests.expect(transform.has_value(),
                    "PlotAdaptive: request/range lower to a finite mm view transform");
                if (transform)
                    tests.expect(transform->viewport().widthMm == 160.0
                            && transform->viewport().heightMm == 100.0,
                        "PlotAdaptive: internal layout extent is expressed in millimeters");
                if (transform) {
                    const auto axes = plot::layoutPlotAxes(*transform);
                    tests.expect(axes
                            && axes.layout->xAxis.placement == plot::PlotAxisPlacement::CrossZero
                            && axes.layout->yAxis.placement == plot::PlotAxisPlacement::CrossZero
                            && !axes.layout->xAxis.majorTicks.empty()
                            && !axes.layout->yAxis.majorTicks.empty(),
                        "PlotAxisLayout: Automatic range padding keeps visible zero-crossing axes");
                }
                if (transform) {
                    plot::AdaptiveSamplingOptions adaptiveOptions;
                    adaptiveOptions.maxRecursion = 8;
                    adaptiveOptions.maxTotalSamples = 1024;
                    adaptiveOptions.chordToleranceMm = 0.15;
                    const auto refined = plot::refinePlotSamples(
                        *polynomialProgram.program, *sampled.curve, *transform,
                        mathematics::defaultAngleSemantics(), adaptiveOptions);
                    tests.expect(refined
                            && refined.curve->segments[0].geometryKind
                                == plot::PlotSegmentGeometryKind::QuadraticBezier,
                        "PlotAdaptive: exact quadratic Bezier geometry bypasses chord refinement");
                    if (refined) {
                        tests.expect(countTagged(*refined.curve, plot::PlotPointTag::YAxisIntercept) == 1,
                            "PlotAdaptive: y-axis intersection is preserved as a path vertex");
                        tests.expect(countTagged(*refined.curve, plot::PlotPointTag::LocalExtremum) >= 1,
                            "PlotAdaptive: local extremum candidate is preserved as a path vertex");
                    }
                }
            }
        }
    
    
        {
            const auto bf = [](std::int64_t value) {
                return numeric::BigFloat::fromBigInt(
                    BigInt{value}, 96, numeric::RoundingMode::NearestEven);
            };
            const plot::PlotViewTransform transform{
                bf(-10), bf(10), bf(-5), bf(5), 96, plot::PlotViewportMm{160.0, 100.0}};
            const auto axes = plot::layoutPlotAxes(transform);
            bool hasZeroX = false;
            bool hasZeroY = false;
            if (axes) {
                for (const auto& tick : axes.layout->xAxis.majorTicks)
                    hasZeroX = hasZeroX || tick.value.isZero();
                for (const auto& tick : axes.layout->yAxis.majorTicks)
                    hasZeroY = hasZeroY || tick.value.isZero();
            }
            tests.expect(axes
                    && axes.layout->xAxis.placement == plot::PlotAxisPlacement::CrossZero
                    && axes.layout->yAxis.placement == plot::PlotAxisPlacement::CrossZero
                    && std::abs(axes.layout->xAxis.axisPositionMm - 50.0) < 1e-12
                    && std::abs(axes.layout->yAxis.axisPositionMm - 80.0) < 1e-12
                    && hasZeroX && hasZeroY,
                "PlotAxisLayout: symmetric ranges place crossing axes and 1-2-5 ticks in mm space");
        }
        {
            const auto bf = [](std::int64_t value) {
                return numeric::BigFloat::fromBigInt(
                    BigInt{value}, 96, numeric::RoundingMode::NearestEven);
            };
            const plot::PlotViewTransform transform{
                bf(2), bf(12), bf(3), bf(8), 96, plot::PlotViewportMm{160.0, 100.0}};
            const auto axes = plot::layoutPlotAxes(transform);
            tests.expect(axes
                    && axes.layout->xAxis.placement == plot::PlotAxisPlacement::MinimumEdge
                    && axes.layout->yAxis.placement == plot::PlotAxisPlacement::MinimumEdge
                    && std::abs(axes.layout->xAxis.axisPositionMm) < 1e-12
                    && std::abs(axes.layout->yAxis.axisPositionMm) < 1e-12,
                "PlotAxisLayout: ranges wholly above/right of zero place axes on minimum edges");
        }
    
        {
            graphics::GraphicsScene scene;
            scene.extent = {20.0, 10.0};
            graphics::GraphicsTextNode text;
            text.origin = {2.0, 3.0};
            text.text = "a<b&c";
            text.semantic = graphics::GraphicsSemanticRef{graphics::GraphicsSemanticKind::Label, 7};
            scene.nodes.push_back(text);
            const auto svg = graphics::renderSvg(scene);
            tests.expect(svg
                    && svg.svg->find("y=\"7\"") != std::string::npos
                    && svg.svg->find("a&lt;b&amp;c") != std::string::npos
                    && svg.svg->find("data-mmcal-id=\"7\"") != std::string::npos,
                "SVG: backend flips y coordinates and XML-escapes text without changing GraphicsScene");
        }


        {
            graphics::GraphicsScene scene;
            scene.extent = {20.0, 10.0};
            graphics::GraphicsPathNode path;
            path.commands.push_back(graphics::GraphicsMoveTo{{2.0, 3.0}});
            path.commands.push_back(graphics::GraphicsQuadraticTo{{6.0, 8.0}, {10.0, 3.0}});
            path.stroke = graphics::GraphicsStrokeStyle{};
            path.semantic = graphics::GraphicsSemanticRef{graphics::GraphicsSemanticKind::Curve, 9};
            scene.nodes.push_back(path);
            const auto eps = graphics::renderEps(scene);
            tests.expect(eps
                    && eps.eps->find("%!PS-Adobe-3.0 EPSF-3.0") == 0
                    && eps.eps->find("%%BoundingBox: 0 0 57 29") != std::string::npos
                    && eps.eps->find("2 3 moveto") != std::string::npos
                    && eps.eps->find("curveto") != std::string::npos
                    && eps.eps->find("clip\n") != std::string::npos
                    && eps.eps->find("% mmCal-semantic: curve 9") != std::string::npos,
                "EPS: backend keeps y-up mm geometry, clips to the canvas, and degree-elevates quadratic paths exactly");

            const auto inferredSvg = graphics::graphicsFormatFromExtension(".SVG");
            const auto inferredEps = graphics::graphicsFormatFromExtension("eps");
            const auto inferredPdf = graphics::graphicsFormatFromExtension(".PDF");
            const auto inferredPng = graphics::graphicsFormatFromExtension(".PNG");
            const auto inferredWebp = graphics::graphicsFormatFromExtension(".WEBP");
            const auto rendered = inferredEps
                ? graphics::renderGraphics(scene, *inferredEps)
                : graphics::GraphicsRenderResult{};
            tests.expect(inferredSvg == graphics::GraphicsFormat::Svg
                    && inferredEps == graphics::GraphicsFormat::Eps
                    && inferredPdf == graphics::GraphicsFormat::Pdf
                    && inferredPng == graphics::GraphicsFormat::Png
                    && inferredWebp == graphics::GraphicsFormat::Webp
                    && rendered && rendered.data
                    && rendered.data->find("%!PS-Adobe-3.0 EPSF-3.0") == 0,
                "GraphicsBackend: case-insensitive format parsing dispatches SVG/EPS/PDF/PNG/WEBP behind one scene renderer");

            graphics::PdfRenderOptions deterministicPdf;
            deterministicPdf.deterministic = true;
            deterministicPdf.userName = "mmcal-test-user";
            const auto pdfA = graphics::renderPdf(scene, deterministicPdf);
            const auto pdfB = graphics::renderPdf(scene, deterministicPdf);
            tests.expect(pdfA && pdfB && pdfA.pdf == pdfB.pdf
                    && pdfA.pdf->find("%PDF-1.4") == 0
                    && pdfA.pdf->find("/Title (mmCal " MMCAL_VERSION_STRING " Plot)") != std::string::npos
                    && pdfA.pdf->find("/Creator (mmCal " MMCAL_VERSION_STRING ")") != std::string::npos
                    && pdfA.pdf->find("/Producer (mmCal " MMCAL_VERSION_STRING " PDF Plotter)") != std::string::npos
                    && pdfA.pdf->find("/Lang (ja-JP)") != std::string::npos
                    && pdfA.pdf->find("<mmcal:UserName>mmcal-test-user</mmcal:UserName>") != std::string::npos
                    && pdfA.pdf->find("<mmcal:Version>" MMCAL_VERSION_STRING "</mmcal:Version>") != std::string::npos
                    && pdfA.pdf->find("D:20000101000000Z") != std::string::npos
                    && pdfA.pdf->find("xref\n0 ") != std::string::npos
                    && pdfA.pdf->find("/ObjStm") == std::string::npos
                    && pdfA.pdf->find("/XRef") == std::string::npos
                    && pdfA.pdf->find("/Filter") == std::string::npos
                    && pdfA.pdf->find(" re W n\n") != std::string::npos
                    && pdfA.pdf->find(" c\n") != std::string::npos,
                "PDF: deterministic PDF 1.4 keeps metadata, classic xref, clipping, and uncompressed readable vector content");

            graphics::GraphicsScene rasterScene;
            rasterScene.extent = {20.0, 10.0};
            graphics::GraphicsPathNode rasterPath;
            rasterPath.commands.push_back(graphics::GraphicsMoveTo{{1.0, 1.0}});
            rasterPath.commands.push_back(graphics::GraphicsCubicTo{
                {5.0, 9.0}, {15.0, 1.0}, {19.0, 9.0}});
            rasterPath.stroke = graphics::GraphicsStrokeStyle{};
            rasterScene.nodes.push_back(rasterPath);
            graphics::GraphicsTextNode rasterText;
            rasterText.origin = {10.0, 5.0};
            rasterText.text = "-0.123456789";
            rasterText.fontSizeMm = 2.0;
            rasterText.anchor = graphics::GraphicsTextAnchor::Middle;
            rasterScene.nodes.push_back(rasterText);
            graphics::RasterRenderOptions pngOptions;
            pngOptions.dpi = 254.0;
            pngOptions.antialiasing = 2;
            const auto png = graphics::renderPng(rasterScene, pngOptions);
            const bool pngHeader = png && png.png && png.png->size() > 64
                && png.png->compare(0, 8, "\x89PNG\r\n\x1a\n", 8) == 0;
            const auto byte = [&](std::size_t offset) {
                return static_cast<unsigned char>((*png.png)[offset]);
            };
            const std::uint32_t pngWidth = pngHeader
                ? (static_cast<std::uint32_t>(byte(16)) << 24u)
                    | (static_cast<std::uint32_t>(byte(17)) << 16u)
                    | (static_cast<std::uint32_t>(byte(18)) << 8u)
                    | static_cast<std::uint32_t>(byte(19))
                : 0u;
            const std::uint32_t pngHeight = pngHeader
                ? (static_cast<std::uint32_t>(byte(20)) << 24u)
                    | (static_cast<std::uint32_t>(byte(21)) << 16u)
                    | (static_cast<std::uint32_t>(byte(22)) << 8u)
                    | static_cast<std::uint32_t>(byte(23))
                : 0u;
            const std::size_t pngIdat = pngHeader ? png.png->find("IDAT") : std::string::npos;
            const bool pngDynamicDeflate = pngIdat != std::string::npos
                && pngIdat + 6u < png.png->size()
                && ((static_cast<unsigned char>((*png.png)[pngIdat + 6u]) >> 1u) & 0x03u) == 2u;
            tests.expect(pngHeader && pngWidth == 200 && pngHeight == 100
                    && png.png->find("pHYs") != std::string::npos
                    && pngIdat != std::string::npos
                    && pngDynamicDeflate
                    && png.png->size() < 20u * 10u * 100u,
                "PNG: raster backend writes 254-dpi dimensions, pHYs, filtered IDAT, and dynamic-Huffman DEFLATE");

            // 長い258-byte matchが頻発する白canvasで，専用length code 285を使う回帰。
            // code 284 + extra=31へ落ちると同じ画像でも大幅に肥大化する。
            graphics::GraphicsScene pngLongRunScene;
            pngLongRunScene.extent = {40.0, 30.0};
            graphics::RasterRenderOptions pngLongRunOptions;
            pngLongRunOptions.dpi = 254.0;
            pngLongRunOptions.antialiasing = 1;
            const auto pngLongRun = graphics::renderPng(pngLongRunScene, pngLongRunOptions);
            tests.expect(pngLongRun && pngLongRun.png && pngLongRun.png->size() < 1600u,
                "PNG DEFLATE: 258-byte LZ77 matches use the dedicated RFC 1951 length code 285");

            const auto webp = graphics::renderWebp(rasterScene, pngOptions);
            const bool webpHeader = webp && webp.webp && webp.webp->size() > 25
                && webp.webp->compare(0, 4, "RIFF", 4) == 0
                && webp.webp->compare(8, 8, "WEBPVP8L", 8) == 0
                && static_cast<unsigned char>((*webp.webp)[20]) == 0x2f;
            // Transform選択はentropy tree最適化で変わり得るため，固定offsetのcache bitには依存しない。
            // main imageのencoder経路自体は常に16-entry color cacheとadaptive prefix treeを使う。
            tests.expect(webpHeader && webp.webp->size() < 12000u,
                "WEBP: lossless backend writes bounded LZ77/color-cache data with adaptive prefix trees");

            graphics::GraphicsScene predictorScene;
            predictorScene.extent = {25.6, 6.4};
            for (unsigned shade = 0; shade < 256u; ++shade) {
                graphics::GraphicsPathNode strip;
                strip.commands.push_back(graphics::GraphicsMoveTo{{shade * 0.1, 0.0}});
                strip.commands.push_back(graphics::GraphicsLineTo{{(shade + 1u) * 0.1, 0.0}});
                strip.commands.push_back(graphics::GraphicsLineTo{{(shade + 1u) * 0.1, 6.4}});
                strip.commands.push_back(graphics::GraphicsLineTo{{shade * 0.1, 6.4}});
                strip.commands.push_back(graphics::GraphicsClosePath{});
                const auto component = static_cast<std::uint8_t>(shade);
                strip.fill = graphics::GraphicsFillStyle{
                    graphics::GraphicsColor{component, component, component, 255u}};
                predictorScene.nodes.push_back(std::move(strip));
            }
            graphics::RasterRenderOptions predictorOptions;
            predictorOptions.widthPx = 256;
            predictorOptions.heightPx = 64;
            predictorOptions.antialiasing = 1;
            const auto predictorWebp = graphics::renderWebp(predictorScene, predictorOptions);
            const bool predictorHeader = predictorWebp && predictorWebp.webp
                && predictorWebp.webp->size() > 26u
                && predictorWebp.webp->compare(0, 4, "RIFF", 4) == 0
                && predictorWebp.webp->compare(8, 8, "WEBPVP8L", 8) == 0;
            // VP8L header直後の最初のtransform bit=1，type=00ならPredictor Transform。
            const bool predictorTransform = predictorHeader
                && (static_cast<unsigned char>((*predictorWebp.webp)[25]) & 0x07u) == 0x01u;
            tests.expect(predictorTransform && predictorWebp.webp->size() < 1000u,
                "WEBP predictor: adaptive VP8L predictor transform is selected for a smooth grayscale ramp");

            // Gを擬似的に並べ替えつつR=G+20, B=G+40とする。
            // 原画像ではRGB各channelが広いalphabetを持つがSubtract Green後はR/Bが定数となるため，
            // adaptive prefix tree込みでもSubtract Greenがdeterministicに有利になる。
            graphics::GraphicsScene subtractGreenScene;
            subtractGreenScene.extent = {256.0, 1.0};
            for (unsigned x = 0u; x < 256u; ++x) {
                const auto green = static_cast<std::uint8_t>((x * 73u) & 0xffu);
                graphics::GraphicsPathNode pixel;
                pixel.commands.push_back(graphics::GraphicsMoveTo{{static_cast<double>(x), 0.0}});
                pixel.commands.push_back(graphics::GraphicsLineTo{{static_cast<double>(x + 1u), 0.0}});
                pixel.commands.push_back(graphics::GraphicsLineTo{{static_cast<double>(x + 1u), 1.0}});
                pixel.commands.push_back(graphics::GraphicsLineTo{{static_cast<double>(x), 1.0}});
                pixel.commands.push_back(graphics::GraphicsClosePath{});
                pixel.fill = graphics::GraphicsFillStyle{graphics::GraphicsColor{
                    static_cast<std::uint8_t>(green + 20u),
                    green,
                    static_cast<std::uint8_t>(green + 40u),
                    255u}};
                subtractGreenScene.nodes.push_back(std::move(pixel));
            }

            graphics::RasterRenderOptions subtractGreenOptions;
            subtractGreenOptions.widthPx = 256u;
            subtractGreenOptions.heightPx = 1u;
            subtractGreenOptions.antialiasing = 1u;
            const auto subtractGreenWebp =
                graphics::renderWebp(subtractGreenScene, subtractGreenOptions);
            // VP8L header直後: transform-present=1, type=10(Subtract Green)なのでlow 3 bitsは101。
            const bool subtractGreenTransform = subtractGreenWebp && subtractGreenWebp.webp
                && subtractGreenWebp.webp->size() > 26u
                && (static_cast<unsigned char>((*subtractGreenWebp.webp)[25]) & 0x07u) == 0x05u;
            tests.expect(subtractGreenTransform && subtractGreenWebp.webp->size() < 256u,
                "WEBP subtract-green: correlated RGB channels select Subtract Green with adaptive prefix trees");

            graphics::GraphicsScene backgroundScene;
            backgroundScene.extent = {2.0, 1.0};
            graphics::RasterRenderOptions whiteBackground;
            whiteBackground.widthPx = 2;
            whiteBackground.heightPx = 1;
            whiteBackground.antialiasing = 1;
            const auto whiteRaster = graphics::renderRaster(backgroundScene, whiteBackground);
            graphics::RasterRenderOptions transparentBackground = whiteBackground;
            transparentBackground.background = graphics::RasterBackground::None;
            const auto transparentRaster = graphics::renderRaster(backgroundScene, transparentBackground);
            tests.expect(whiteRaster && transparentRaster
                    && whiteRaster.image->rgba == std::vector<std::uint8_t>{
                        255, 255, 255, 255, 255, 255, 255, 255}
                    && transparentRaster.image->rgba == std::vector<std::uint8_t>{
                        0, 0, 0, 0, 0, 0, 0, 0},
                "PNG background: White is default and None preserves a fully transparent RGBA canvas");

            graphics::GraphicsScene unsupportedTextScene = rasterScene;
            std::get<graphics::GraphicsTextNode>(unsupportedTextScene.nodes.back()).text = "1e3";
            const auto unsupportedText = graphics::renderPng(unsupportedTextScene, pngOptions);
            tests.expect(unsupportedText.status == graphics::PngRenderStatus::UnsupportedText,
                "PNG text: built-in stroke font is deliberately limited to digits, minus, and decimal point");

            graphics::GraphicsScene transparentScene;
            transparentScene.extent = {10.0, 10.0};
            graphics::GraphicsCircleNode circle;
            circle.center = {5.0, 5.0};
            circle.radiusMm = 1.0;
            circle.fill = graphics::GraphicsFillStyle{
                graphics::GraphicsColor{0, 0, 0, 128}};
            transparentScene.nodes.push_back(circle);
            const auto unsupported = graphics::renderGraphics(
                transparentScene, graphics::GraphicsFormat::Eps);
            tests.expect(unsupported.status == graphics::GraphicsRenderStatus::UnsupportedFeature,
                "EPS: unsupported transparency fails explicitly instead of silently changing appearance");
        }
    
        const auto extremumRequest = request(session, "(x-1/3)^2+1", "-1", "1");
        const auto extremumProgram = compile(session, extremumRequest);
        if (extremumProgram) {
            const auto extremumAnalysis = plot::makeInitialPlotAnalysis(extremumRequest);
            const auto sampled = plot::coarseSamplePlot(
                extremumRequest, extremumAnalysis, *extremumProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), plot::CoarseSamplingOptions{96, 5, 128});
            if (sampled) {
                const auto inspection = plot::inspectCoarseSamples(extremumAnalysis, *sampled.curve);
                const auto range = plot::estimatePlotRange(extremumAnalysis, *sampled.curve, inspection);
                const auto transform = range ? plot::makePlotViewTransform(
                    extremumRequest, *range,
                    session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), 96) : std::nullopt;
                if (transform) {
                    plot::AdaptiveSamplingOptions adaptiveOptions;
                    adaptiveOptions.maxTotalSamples = 2048;
                    adaptiveOptions.extremumRefinementIterations = 10;
                    const auto refined = plot::refinePlotSamples(
                        *extremumProgram.program, *sampled.curve, *transform,
                        mathematics::defaultAngleSemantics(), adaptiveOptions);
                    tests.expect(refined
                            && countTagged(*refined.curve, plot::PlotPointTag::LocalExtremum) >= 1
                            && refined.statistics.localExtrema >= 1,
                        "PlotAdaptive: off-grid local extrema are numerically refined into path vertices");
                }
            }
        }
    
        const auto axisRootsRequest = request(session, "x^2-1/5", "-1", "1");
        const auto rootsProgram = compile(session, axisRootsRequest);
        if (rootsProgram) {
            const auto rootsAnalysis = plot::makeInitialPlotAnalysis(axisRootsRequest);
            plot::CoarseSamplingOptions rootsCoarse{96, 5, 128};
            const auto sampled = plot::coarseSamplePlot(
                axisRootsRequest, rootsAnalysis, *rootsProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), rootsCoarse);
            if (sampled) {
                const auto inspection = plot::inspectCoarseSamples(rootsAnalysis, *sampled.curve);
                const auto range = plot::estimatePlotRange(rootsAnalysis, *sampled.curve, inspection);
                const auto transform = range ? plot::makePlotViewTransform(
                    axisRootsRequest, *range,
                    session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), 96) : std::nullopt;
                if (transform) {
                    plot::AdaptiveSamplingOptions adaptiveOptions;
                    adaptiveOptions.maxTotalSamples = 1024;
                    const auto refined = plot::refinePlotSamples(
                        *rootsProgram.program, *sampled.curve, *transform,
                        mathematics::defaultAngleSemantics(), adaptiveOptions);
                    tests.expect(refined
                            && countTagged(*refined.curve, plot::PlotPointTag::XAxisIntercept) == 2,
                        "PlotAdaptive: bracketed x-axis intersections become explicit path vertices");
                }
            }
        }
    
        if (oscillatoryProgram) {
            const auto oscillatoryAnalysis = analyze(session, oscillatoryRequest);
            plot::CoarseSamplingOptions oscillatoryCoarse{96, 17, 256};
            const auto sampled = plot::coarseSamplePlot(
                oscillatoryRequest, oscillatoryAnalysis, *oscillatoryProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), oscillatoryCoarse);
            if (sampled) {
                const auto inspection = plot::inspectCoarseSamples(oscillatoryAnalysis, *sampled.curve);
                const auto range = plot::estimatePlotRange(oscillatoryAnalysis, *sampled.curve, inspection);
                const auto transform = range ? plot::makePlotViewTransform(
                    oscillatoryRequest, *range,
                    session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), 96,
                    plot::PlotViewportMm{160.0, 100.0}) : std::nullopt;
                if (transform) {
                    plot::AdaptiveSamplingOptions adaptiveOptions;
                    adaptiveOptions.maxRecursion = 12;
                    adaptiveOptions.maxTotalSamples = 4096;
                    adaptiveOptions.minimumSpanMm = 0.15;
                    const auto refined = plot::refinePlotSamples(
                        *oscillatoryProgram.program, *sampled.curve, *transform,
                        mathematics::defaultAngleSemantics(), adaptiveOptions);
                    tests.expect(refined && refined.curve->segments.size() == 2
                            && refined.statistics.insertedSamples > 0,
                        "PlotAdaptive: sin(1/x) refines without reconnecting across the puncture");
                    tests.expect(!refined || refined.curve->segments[0].samples.size()
                            + refined.curve->segments[1].samples.size() <= adaptiveOptions.maxTotalSamples,
                        "PlotAdaptive: infinitely rapid oscillation remains bounded by mm-layout/resource limits");
                }
            }
        }
    
        if (affineProgram) {
            const auto sampled = plot::coarseSamplePlot(
                affineRequest, plot::makeInitialPlotAnalysis(affineRequest),
                *affineProgram.program, session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            if (sampled) {
                const auto inspection = plot::inspectCoarseSamples(
                    plot::makeInitialPlotAnalysis(affineRequest), *sampled.curve);
                const auto range = plot::estimatePlotRange(
                    plot::makeInitialPlotAnalysis(affineRequest), *sampled.curve, inspection);
                const auto transform = range ? plot::makePlotViewTransform(
                    affineRequest, *range, session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), 96) : std::nullopt;
                if (transform) {
                    plot::AdaptiveSamplingOptions noAnchors;
                    noAnchors.preserveAxisIntersections = false;
                    noAnchors.preserveLocalExtrema = false;
                    const auto refined = plot::refinePlotSamples(
                        *affineProgram.program, *sampled.curve, *transform,
                        mathematics::defaultAngleSemantics(), noAnchors);
                    tests.expect(refined && refined.curve->segments[0].samples.size() == 2
                            && refined.statistics.insertedSamples == 0,
                        "PlotAdaptive: affine straight segments bypass chord refinement");
                    if (refined) {
                        const auto plotScene = plot::buildPlotScene(
                            std::vector<plot::SampledCurve>{*refined.curve});
                        const auto axes = plot::layoutPlotAxes(*transform);
                        const auto graphics = plotScene && axes
                            ? plot::lowerPlotToGraphics(*plotScene.scene, *axes.layout, *transform)
                            : plot::PlotGraphicsLoweringResult{};
                        bool hasCurvePath = false;
                        bool hasXAxis = false;
                        bool hasYAxis = false;
                        bool hasXTick = false;
                        bool hasYTick = false;
                        bool defaultStrokeWidthsMatch = true;
                        bool hasXTickLabel = false;
                        bool hasYTickLabel = false;
                        if (graphics) {
                            for (const auto& node : graphics.scene->nodes) {
                                if (const auto* text = std::get_if<graphics::GraphicsTextNode>(&node);
                                    text && text->semantic) {
                                    hasXTickLabel = hasXTickLabel
                                        || text->semantic->kind == graphics::GraphicsSemanticKind::XTickLabel;
                                    hasYTickLabel = hasYTickLabel
                                        || text->semantic->kind == graphics::GraphicsSemanticKind::YTickLabel;
                                }
                                const auto* path = std::get_if<graphics::GraphicsPathNode>(&node);
                                if (!path || !path->semantic)
                                    continue;
                                const auto kind = path->semantic->kind;
                                hasCurvePath = hasCurvePath || kind == graphics::GraphicsSemanticKind::Curve;
                                hasXAxis = hasXAxis || kind == graphics::GraphicsSemanticKind::XAxis;
                                hasYAxis = hasYAxis || kind == graphics::GraphicsSemanticKind::YAxis;
                                hasXTick = hasXTick || kind == graphics::GraphicsSemanticKind::XTick;
                                hasYTick = hasYTick || kind == graphics::GraphicsSemanticKind::YTick;
                                if (!path->stroke)
                                    defaultStrokeWidthsMatch = false;
                                else if (kind == graphics::GraphicsSemanticKind::Curve)
                                    defaultStrokeWidthsMatch = defaultStrokeWidthsMatch
                                        && path->stroke->widthMm == 0.8;
                                else if (kind == graphics::GraphicsSemanticKind::XAxis
                                    || kind == graphics::GraphicsSemanticKind::YAxis)
                                    defaultStrokeWidthsMatch = defaultStrokeWidthsMatch
                                        && path->stroke->widthMm == 0.5;
                                else if (kind == graphics::GraphicsSemanticKind::XTick
                                    || kind == graphics::GraphicsSemanticKind::YTick)
                                    defaultStrokeWidthsMatch = defaultStrokeWidthsMatch
                                        && path->stroke->widthMm == 0.3;
                            }
                        }
                        tests.expect(graphics
                                && graphics.scene->extent.widthMm == transform->viewport().widthMm
                                && graphics.scene->extent.heightMm == transform->viewport().heightMm
                                && hasCurvePath && hasXAxis && hasYAxis && hasXTick && hasYTick
                                && hasXTickLabel && hasYTickLabel && defaultStrokeWidthsMatch,
                            "PlotGraphics: label margins remain inside the fixed physical canvas");
                        if (graphics) {
                            const auto svg = graphics::renderSvg(*graphics.scene);
                            tests.expect(svg
                                    && svg.svg->find("data-mmcal-kind=\"curve\"") != std::string::npos
                                    && svg.svg->find("data-mmcal-kind=\"x-axis\"") != std::string::npos
                                    && svg.svg->find("<g data-mmcal-kind=\"x-ticks\">") != std::string::npos
                                    && svg.svg->find("<g data-mmcal-kind=\"y-ticks\">") != std::string::npos
                                    && svg.svg->find("<g data-mmcal-kind=\"x-tick-labels\">") != std::string::npos
                                    && svg.svg->find("<g data-mmcal-kind=\"y-tick-labels\">") != std::string::npos
                                    && svg.svg->find("font-family=\"sans-serif\"") != std::string::npos
                                    && svg.svg->find("marker-start=") == std::string::npos
                                    && svg.svg->find("marker-end=") == std::string::npos,
                                "SVG: tick paths and labels remain grouped without axis arrow markers");
                        }
                    }
                }
            }
        }
    
    
        if (polynomialProgram) {
            plot::CoarseSamplingOptions smallBudget = options;
            smallBudget.maxTotalSamples = 4;
            const auto sampled = plot::coarseSamplePlot(
                polynomialRequest, plot::makeInitialPlotAnalysis(polynomialRequest),
                *polynomialProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), smallBudget);
            tests.expect(sampled.status == plot::PlotSamplingStatus::ResourceLimit,
                "PlotSampling: total coarse sample count has an explicit resource bound");
        }
    };
    runCoreSampling();

    const auto runPipelineAndGraphics = [&] {
        // orchestratorは単一curveの内部工程をSVGまで一度に通す。
        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto sineRequest = request(session, "sin[x]", "-Pi", "Pi");
            plot::PlotPipelineOptions pipelineOptions;
            pipelineOptions.coarseSampling.precisionBits = 96;
            const auto rendered = plot::renderPlotSvg(
                sineRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol,
                {}, pipelineOptions);
            std::size_t sineVertexCount = 0;
            std::size_t sineEndpointMarkerCount = 0;
            if (rendered && !rendered.output->curves.empty()
                && !rendered.output->curves.front().segments.empty())
                sineVertexCount = rendered.output->curves.front().segments.front().samples.size();
            if (rendered) {
                for (const auto& node : rendered.output->graphicsScene.nodes)
                    if (std::get_if<graphics::GraphicsCircleNode>(&node))
                        ++sineEndpointMarkerCount;
            }
            tests.expect(rendered
                    && rendered.output->curves.size() == 1
                    && rendered.output->plotScene.curves.size() == 1
                    && rendered.output->graphicsScene.extent.widthMm == 150.0
                    && rendered.output->graphicsScene.extent.heightMm == 100.0
                    && sineVertexCount >= 90 && sineVertexCount <= 110
                    && rendered.svg->find("stroke-width=\"0.8\"") != std::string::npos
                    && rendered.svg->find("stroke-width=\"0.5\"") != std::string::npos
                    && rendered.svg->find("stroke-width=\"0.3\"") != std::string::npos
                    && rendered.svg->find("<g data-mmcal-kind=\"x-ticks\">") != std::string::npos
                    && rendered.svg->find("<g data-mmcal-kind=\"y-ticks\">") != std::string::npos
                    && rendered.svg->find("<g data-mmcal-kind=\"x-tick-labels\">") != std::string::npos
                    && rendered.svg->find("<g data-mmcal-kind=\"y-tick-labels\">") != std::string::npos
                    && rendered.svg->find(">-0.5</text>") != std::string::npos
                    && rendered.svg->find(">0</text>") == std::string::npos
                    && sineEndpointMarkerCount == 0
                    && rendered.svg->find("width=\"150mm\" height=\"100mm\"") != std::string::npos,
                "PlotPipeline: default sin[x] uses about 100 samples on a fixed 150x100 mm canvas without origin ticks or graph-edge markers");
        }
    
        // Plot用簡約は元式のcomplete domainを保持したまま評価式だけへ適用する。
        // sin[1/x]^2+cos[1/x]^2はx!=0で1だが，x=0のholeを消してはいけない。
        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto identityRequest = request(
                session, "sin[1/x]^2+cos[1/x]^2", "-Pi", "Pi");
            const auto built = plot::buildPlotPipeline(
                identityRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            bool domainPreserved = false;
            bool constantGeometry = false;
            if (built && built.output->analyses.size() == 1
                && built.output->programs.size() == 1
                && built.output->curves.size() == 1) {
                const auto& analysis = built.output->analyses.front();
                const auto& curve = built.output->curves.front();
                domainPreserved = analysis.domain.coverage == plot::PlotDomainCoverage::Complete
                    && analysis.domain.intervals.size() == 2
                    && curve.segments.size() == 2;
                constantGeometry = built.output->programs.front().geometryKind
                        == plot::PlotCurveGeometryKind::Constant
                    && std::all_of(
                        curve.segments.begin(), curve.segments.end(),
                        [](const plot::SampledCurveSegment& segment) {
                            return segment.geometryKind == plot::PlotSegmentGeometryKind::StraightLine
                                && segment.samples.size() == 2
                                && segment.samples[0].finite() && segment.samples[1].finite()
                                && segment.samples[0].y == segment.samples[1].y;
                        });
                if (domainPreserved && constantGeometry) {
                    const auto leftNearHole = built.output->transform.mapX(
                        curve.segments[0].samples.back().x);
                    const auto rightNearHole = built.output->transform.mapX(
                        curve.segments[1].samples.front().x);
                    const auto zero = built.output->transform.mapX(numeric::BigFloat::fromBigInt(
                        BigInt{0}, curve.precisionBits, numeric::RoundingMode::NearestEven));
                    constantGeometry = leftNearHole && rightNearHole && zero
                        && *leftNearHole < *zero && *zero < *rightNearHole
                        && *rightNearHole - *leftNearHole < 0.5;
                }
            }
            tests.expect(built && domainPreserved && constantGeometry,
                "PlotPipeline: domain-preserving simplification keeps the x=0 hole while reducing the trigonometric identity to constant lines");
        }
    
        // endpoint markerは開区間端を白抜き，閉区間端とsingleton pointを塗りつぶしで描く。
        // sin[1/x]^2+cos[1/x]^2はx=0のholeを1個の開丸へ集約し，asin[sec[x]]は3個の点だけを示す。
        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto identityRequest = request(session, "sin[1/x]^2+cos[1/x]^2", "-Pi", "Pi");
            const auto identityRendered = plot::renderPlotSvg(
                identityRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            const auto secantRequest = request(session, "asin[sec[x]]", "-Pi", "Pi");
            const auto secantRendered = plot::renderPlotSvg(
                secantRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            const auto floorRequest = request(session, "floor[x]", "-2", "2");
            const auto floorRendered = plot::renderPlotSvg(
                floorRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            bool identityHasSingleOpenHole = false;
            bool secantHasThreePoints = false;
            bool floorHasMixedOpenClosed = false;
            if (identityRendered) {
                std::size_t openMarkers = 0;
                for (const auto& node : identityRendered.output->graphicsScene.nodes) {
                    const auto* circle = std::get_if<graphics::GraphicsCircleNode>(&node);
                    if (!circle || !circle->fill)
                        continue;
                    if (circle->fill->color.red == 255
                        && circle->fill->color.green == 255
                        && circle->fill->color.blue == 255)
                        ++openMarkers;
                }
                identityHasSingleOpenHole = openMarkers == 1;
            }
            if (secantRendered) {
                std::size_t points = 0;
                for (const auto& node : secantRendered.output->graphicsScene.nodes) {
                    const auto* circle = std::get_if<graphics::GraphicsCircleNode>(&node);
                    if (circle && circle->semantic
                        && circle->semantic->kind == graphics::GraphicsSemanticKind::Point)
                        ++points;
                }
                secantHasThreePoints = points == 3;
            }
            if (floorRendered) {
                bool hasOpen = false;
                bool hasClosed = false;
                for (const auto& node : floorRendered.output->graphicsScene.nodes) {
                    const auto* circle = std::get_if<graphics::GraphicsCircleNode>(&node);
                    if (!circle || !circle->fill)
                        continue;
                    const auto& c = circle->fill->color;
                    hasOpen = hasOpen || (c.red == 255 && c.green == 255 && c.blue == 255);
                    hasClosed = hasClosed || (c.red == 0 && c.green == 0 && c.blue == 0);
                }
                floorHasMixedOpenClosed = hasOpen && hasClosed;
            }
            tests.expect(identityRendered && secantRendered && floorRendered
                    && identityHasSingleOpenHole && secantHasThreePoints && floorHasMixedOpenClosed,
                "PlotGraphics: endpoint markers use one open hole marker, filled singleton points, and mixed open/closed step endpoints");
        }

        // generic polylineでも，元式のopen boundaryに有限な片側極限があれば
        // 未定義点をpath sampleへ混ぜずに白丸だけを数学的境界へ置く。
        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto countOpenMarkers = [&](std::string_view source) {
                const auto rendered = plot::renderPlotSvg(
                    request(session, source, "-3", "3"),
                    session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), infinity->symbol);
                if (!rendered)
                    return std::optional<std::size_t>{};
                std::size_t openMarkers = 0;
                for (const auto& node : rendered.output->graphicsScene.nodes) {
                    const auto* circle = std::get_if<graphics::GraphicsCircleNode>(&node);
                    if (!circle || !circle->fill)
                        continue;
                    const auto& c = circle->fill->color;
                    if (c.red == 255 && c.green == 255 && c.blue == 255)
                        ++openMarkers;
                }
                return std::optional<std::size_t>{openMarkers};
            };

            const auto xOverX = countOpenMarkers("x/x");
            const auto factorHole = countOpenMarkers("(x^2-1)/(x-1)");
            const auto reciprocalHole = countOpenMarkers("1/(1/x)");
            const auto zeroPowerHole = countOpenMarkers("x^0");
            const auto oscillatory = countOpenMarkers("sin[1/x]");
            const auto selfPowerRendered = plot::renderPlotSvg(
                request(session, "x^x", "0", "3"),
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            std::optional<std::size_t> selfPower;
            if (selfPowerRendered) {
                std::size_t openMarkers = 0;
                for (const auto& node : selfPowerRendered.output->graphicsScene.nodes) {
                    const auto* circle = std::get_if<graphics::GraphicsCircleNode>(&node);
                    if (!circle || !circle->fill)
                        continue;
                    const auto& c = circle->fill->color;
                    if (c.red == 255 && c.green == 255 && c.blue == 255)
                        ++openMarkers;
                }
                selfPower = openMarkers;
            }
            tests.expect(xOverX && *xOverX == 1
                    && factorHole && *factorHole == 1
                    && reciprocalHole && *reciprocalHole == 1
                    && zeroPowerHole && *zeroPowerHole == 1
                    && oscillatory && *oscillatory == 0
                    && selfPower && *selfPower == 1,
                "PlotGraphics: finite one-sided limits create generic removable-hole markers, including request-edge holes such as x^x at x=0, without inventing one for sin[1/x]");
        }

        // 高密度交点をmm解像度で打ち切る場合も，近似交点へ既存curveを強制snapして
        // smooth側のgeometryを壊してはいけない。sin[x]はsin[1/x]との併記でも本来の形を保つ。
        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const std::vector<plot::PlotRequest> requests{
                request(session, "sin[x]", "-Pi", "Pi"),
                request(session, "sin[1/x]", "-Pi", "Pi")};
            const auto built = plot::buildPlotPipeline(
                requests,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            bool smoothCurvePreserved = false;
            if (built && built.output->curves.size() == 2 && built.output->programs.size() == 2) {
                try {
                    plot::BigFloatPlotExecutor executor{
                        built.output->programs[0], built.output->curves[0].precisionBits,
                        mathematics::defaultAngleSemantics()};
                    smoothCurvePreserved = true;
                    for (const auto& segment : built.output->curves[0].segments) {
                        for (const auto& sample : segment.samples) {
                            const auto exact = executor.evaluate(sample.x);
                            const auto sampledY = built.output->transform.mapY(sample.y);
                            const auto exactY = exact.finite()
                                ? built.output->transform.mapY(exact.value) : std::nullopt;
                            if (!sampledY || !exactY || std::abs(*sampledY - *exactY) > 0.051) {
                                smoothCurvePreserved = false;
                                break;
                            }
                        }
                        if (!smoothCurvePreserved)
                            break;
                    }
                }
                catch (...) {
                    smoothCurvePreserved = false;
                }
            }
            tests.expect(built && smoothCurvePreserved,
                "PlotIntersections: dense sin[x]/sin[1/x] crossings never visibly distort the smooth curve while creating semantic anchors");
        }
    
        // 明示角度単位のdomain endpointはactive angle unitへ変換してからsamplingする。
        // Radian既定では±180Degが±Piへ落ちるため，一周期を巨大な180-radian区間として扱わない。
        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto degreeRequest = request(session, "sin[x]", "-180Deg", "180Deg");
            const auto rendered = plot::renderPlotSvg(
                degreeRequest,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol);
            std::size_t samples = 0;
            if (rendered && !rendered.output->curves.empty()
                && !rendered.output->curves.front().segments.empty())
                samples = rendered.output->curves.front().segments.front().samples.size();
            tests.expect(rendered && samples >= 90 && samples <= 110,
                "PlotPipeline: explicit degree endpoints render one sine period with normal sample density");
        }
    
        // multi-curve rangeは各flat curveの局所paddingではなく，統合data spanへ一度だけpaddingする。
        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto x = session.evaluate("x").asSymbol();
            const plot::PlotRequestSet flatSet{
                {session.evaluate("100"), session.evaluate("101")},
                x, session.evaluate("-1"), session.evaluate("1")};
            plot::PlotPipelineOptions pipelineOptions;
            pipelineOptions.coarseSampling.precisionBits = 96;
            const auto built = plot::buildPlotPipeline(
                flatSet,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol,
                {}, pipelineOptions);
            bool tightCombinedRange = false;
            if (built) {
                const auto y100 = built.output->transform.mapY(numeric::BigFloat::fromBigInt(
                    BigInt{100}, 96, numeric::RoundingMode::NearestEven));
                const auto y101 = built.output->transform.mapY(numeric::BigFloat::fromBigInt(
                    BigInt{101}, 96, numeric::RoundingMode::NearestEven));
                tightCombinedRange = y100 && y101 && *y100 > 4.0 && *y100 < 6.0
                    && *y101 > 94.0 && *y101 < 96.0;
            }
            tests.expect(built && tightCombinedRange,
                "PlotPipeline: multi-curve Automatic range pads the combined data extent only once");
        }
    
        // 異なるPlot区間を合成するとviewportは大域包絡へ統一するが，各曲線のsampling区間は延長しない。
        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto squareRequest = request(session, "x^2", "-2", "2");
            const auto tangentRequest = request(session, "tan[x]", "-Pi", "Pi");
            plot::PlotPipelineOptions pipelineOptions;
            pipelineOptions.coarseSampling.precisionBits = 96;
            const auto built = plot::buildPlotPipeline(
                std::vector<plot::PlotRequest>{squareRequest, tangentRequest},
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol,
                {}, pipelineOptions);
            bool preservedLocalDomain = false;
            bool globalViewport = false;
            if (built && built.output->curves.size() == 2
                && !built.output->curves[0].segments.empty()) {
                const auto& squareSegments = built.output->curves[0].segments;
                const auto& first = squareSegments.front().samples.front();
                const auto& last = squareSegments.back().samples.back();
                preservedLocalDomain = first.x.toRational() == numeric::Rational{BigInt{-2}}
                    && last.x.toRational() == numeric::Rational{BigInt{2}};
    
                const auto left = built.output->transform.mapX(first.x);
                const auto right = built.output->transform.mapX(last.x);
                globalViewport = left && right && *left > 0.0
                    && *right < built.output->transform.viewport().widthMm;
            }
            tests.expect(built && preservedLocalDomain && globalViewport,
                "PlotPipeline: mixed plot domains use a global x viewport without extrapolating local curves");
        }
    
        // multi-curve orchestratorはintersectionをshared PlotAnchorとしてSceneまで保持する。
        if (const auto* infinity = session.symbolRegistry().find("Infinity")) {
            const auto x = session.evaluate("x").asSymbol();
            const plot::PlotRequestSet crossingSet{
                {session.evaluate("x"), session.evaluate("-x")},
                x, session.evaluate("-1"), session.evaluate("1")};
            plot::PlotPipelineOptions pipelineOptions;
            pipelineOptions.coarseSampling.precisionBits = 96;
            const auto built = plot::buildPlotPipeline(
                crossingSet,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), infinity->symbol,
                {}, pipelineOptions);
            bool sharedIntersection = false;
            if (built) {
                for (const auto& anchor : built.output->plotScene.anchors) {
                    if (plot::hasPlotAnchorKind(anchor.kinds, plot::PlotAnchorKind::CurveIntersection)
                        && anchor.curveIds.size() == 2 && anchor.vertices.size() == 2) {
                        sharedIntersection = true;
                        break;
                    }
                }
            }
            tests.expect(built && sharedIntersection,
                "PlotPipeline: multi-curve intersections survive as shared semantic anchors");
        }
    
    };
    runPipelineAndGraphics();

    const auto runIntersections = [&] {
        const auto xProgram = compile(session, request(session, "x", "-1", "1"));
        const auto minusXProgram = compile(session, request(session, "-x", "-1", "1"));
        if (xProgram && minusXProgram) {
            const auto xRequest = request(session, "x", "-1", "1");
            const auto minusXRequest = request(session, "-x", "-1", "1");
            const auto xSamples = plot::coarseSamplePlot(
                xRequest, plot::makeInitialPlotAnalysis(xRequest), *xProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            const auto minusXSamples = plot::coarseSamplePlot(
                minusXRequest, plot::makeInitialPlotAnalysis(minusXRequest), *minusXProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            if (xSamples && minusXSamples) {
                std::vector<plot::PlotProgram> programs{*xProgram.program, *minusXProgram.program};
                std::vector<plot::SampledCurve> curves{*xSamples.curve, *minusXSamples.curve};
                const auto intersections = plot::refineCurveIntersections(
                    programs, curves, mathematics::defaultAngleSemantics());
                tests.expect(intersections && intersections.statistics.intersections == 1
                        && countTagged((*intersections.curves)[0], plot::PlotPointTag::CurveIntersection) == 1
                        && countTagged((*intersections.curves)[1], plot::PlotPointTag::CurveIntersection) == 1,
                    "PlotIntersections: crossing curves share an explicit editable intersection vertex");
                if (intersections) {
                    bool atOrigin = false;
                    for (const auto& sample : (*intersections.curves)[0].segments[0].samples)
                        if (plot::hasPlotPointTag(sample.tags, plot::PlotPointTag::CurveIntersection)
                            && sample.x.toRational().isZero() && sample.y.toRational().isZero())
                            atOrigin = true;
                    tests.expect(atOrigin,
                        "PlotIntersections: affine x and -x intersection is retained at the exact origin");
                }
            }
        }
    
    
        // 交点は2曲線だけでなくpairwiseに3本以上を処理し，共通座標へ正規化する。
        const auto zeroProgram = compile(session, request(session, "0", "-2", "2"));
        const auto shiftedSquareProgram = compile(session, request(session, "(x-1/3)^2", "-2", "2"));
        if (xProgram && minusXProgram && zeroProgram) {
            const std::vector<plot::PlotRequest> requests{
                request(session, "x", "-1", "1"),
                request(session, "-x", "-1", "1"),
                request(session, "0", "-1", "1")};
            std::vector<plot::PlotProgram> programs{
                *xProgram.program, *minusXProgram.program, *zeroProgram.program};
            std::vector<plot::SampledCurve> curves;
            for (std::size_t i = 0; i < requests.size(); ++i) {
                const auto sampled = plot::coarseSamplePlot(
                    requests[i], plot::makeInitialPlotAnalysis(requests[i]), programs[i],
                    session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), options);
                if (sampled)
                    curves.push_back(*sampled.curve);
            }
            if (curves.size() == 3) {
                const auto intersections = plot::refineCurveIntersections(
                    programs, curves, mathematics::defaultAngleSemantics());
                const auto anchors = intersections
                    ? plot::collectPlotAnchors(*intersections.curves) : plot::PlotAnchorSet{};
                bool sharedOrigin = false;
                for (const auto& anchor : anchors.anchors)
                    if (plot::hasPlotAnchorKind(anchor.kinds, plot::PlotAnchorKind::CurveIntersection)
                        && anchor.x.toRational().isZero() && anchor.y.toRational().isZero()
                        && anchor.curveIndices.size() == 3 && anchor.vertices.size() == 3)
                        sharedOrigin = true;
                tests.expect(intersections && intersections.statistics.intersections == 3 && sharedOrigin,
                    "PlotIntersections: three curves collapse pairwise origin hits into one shared anchor");
                if (intersections) {
                    const auto scene = plot::buildPlotScene(*intersections.curves, anchors);
                    bool linked = false;
                    if (scene && scene.scene->anchors.size() == 1) {
                        const auto& anchor = scene.scene->anchors[0];
                        linked = anchor.curveIds.size() == 3 && anchor.vertices.size() == 3;
                        for (const auto& ref : anchor.vertices) {
                            if (ref.curveId == 0 || ref.curveId > scene.scene->curves.size()) {
                                linked = false;
                                break;
                            }
                            const auto& vertex = scene.scene->curves[ref.curveId - 1]
                                .segments[ref.segmentIndex].vertices[ref.vertexIndex];
                            linked = linked && vertex.anchorId == anchor.id
                                && vertex.x == anchor.x && vertex.y == anchor.y;
                        }
                    }
                    tests.expect(linked,
                        "PlotScene: shared anchors and curve vertices keep bidirectional identity");
                }
            }
        }
    
        // off-gridの偶数重根は符号反転しないため，|f-g|の局所極小をbounded refinementする。
        if (shiftedSquareProgram && zeroProgram) {
            const auto squareRequest = request(session, "(x-1/3)^2", "-1", "1");
            const auto zeroRequest = request(session, "0", "-1", "1");
            plot::CoarseSamplingOptions tangencyOptions{96, 6, 64};
            const auto squareSamples = plot::coarseSamplePlot(
                squareRequest, plot::makeInitialPlotAnalysis(squareRequest), *shiftedSquareProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), tangencyOptions);
            const auto zeroSamples = plot::coarseSamplePlot(
                zeroRequest, plot::makeInitialPlotAnalysis(zeroRequest), *zeroProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), tangencyOptions);
            if (squareSamples && zeroSamples) {
                const auto intersections = plot::refineCurveIntersections(
                    std::vector<plot::PlotProgram>{*shiftedSquareProgram.program, *zeroProgram.program},
                    std::vector<plot::SampledCurve>{*squareSamples.curve, *zeroSamples.curve},
                    mathematics::defaultAngleSemantics());
                tests.expect(intersections && intersections.statistics.intersections == 1
                        && intersections.statistics.tangentialIntersections == 1
                        && countTagged((*intersections.curves)[0], plot::PlotPointTag::CurveIntersection) == 1
                        && countTagged((*intersections.curves)[1], plot::PlotPointTag::CurveIntersection) == 1,
                    "PlotIntersections: off-grid tangential roots become explicit shared vertices");
                if (intersections) {
                    const auto anchors = plot::collectPlotAnchors(*intersections.curves);
                    tests.expect(anchors.anchors.size() == 1
                            && plot::hasPlotAnchorKind(
                                anchors.anchors[0].kinds, plot::PlotAnchorKind::CurveIntersection)
                            && anchors.anchors[0].curveIndices.size() == 2
                            && anchors.anchors[0].vertices.size() == 2,
                        "PlotAnchors: tangential intersection is represented by one shared semantic object");
                }
            }
        }
    
        // 小さいが非零の局所minimumは接触交点へ昇格させない。
        const auto nearMissProgram = compile(session, request(session, "(x-1/3)^2+1/1000", "-1", "1"));
        if (nearMissProgram && zeroProgram) {
            const auto nearMissRequest = request(session, "(x-1/3)^2+1/1000", "-1", "1");
            const auto zeroRequest = request(session, "0", "-1", "1");
            plot::CoarseSamplingOptions tangencyOptions{96, 6, 64};
            const auto nearMissSamples = plot::coarseSamplePlot(
                nearMissRequest, plot::makeInitialPlotAnalysis(nearMissRequest), *nearMissProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), tangencyOptions);
            const auto zeroSamples = plot::coarseSamplePlot(
                zeroRequest, plot::makeInitialPlotAnalysis(zeroRequest), *zeroProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), tangencyOptions);
            if (nearMissSamples && zeroSamples) {
                const auto intersections = plot::refineCurveIntersections(
                    std::vector<plot::PlotProgram>{*nearMissProgram.program, *zeroProgram.program},
                    std::vector<plot::SampledCurve>{*nearMissSamples.curve, *zeroSamples.curve},
                    mathematics::defaultAngleSemantics());
                tests.expect(intersections && intersections.statistics.intersections == 0
                        && intersections.statistics.tangencyCandidates >= 1
                        && intersections.statistics.rejectedTangencyCandidates >= 1,
                    "PlotIntersections: nonzero near-miss minima are not misclassified as tangencies");
            }
        }
    
        // 一組に複数の横切り交点があっても全て保持する。
        const auto cubicProgram = compile(session, request(session, "x^3-x", "-2", "2"));
        if (cubicProgram && zeroProgram) {
            const auto cubicRequest = request(session, "x^3-x", "-2", "2");
            const auto zeroRequest = request(session, "0", "-2", "2");
            plot::CoarseSamplingOptions multiRootOptions{96, 9, 64};
            const auto cubicSamples = plot::coarseSamplePlot(
                cubicRequest, plot::makeInitialPlotAnalysis(cubicRequest), *cubicProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), multiRootOptions);
            const auto zeroSamples = plot::coarseSamplePlot(
                zeroRequest, plot::makeInitialPlotAnalysis(zeroRequest), *zeroProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), multiRootOptions);
            if (cubicSamples && zeroSamples) {
                const auto intersections = plot::refineCurveIntersections(
                    std::vector<plot::PlotProgram>{*cubicProgram.program, *zeroProgram.program},
                    std::vector<plot::SampledCurve>{*cubicSamples.curve, *zeroSamples.curve},
                    mathematics::defaultAngleSemantics());
                tests.expect(intersections && intersections.statistics.intersections == 3,
                    "PlotIntersections: one curve pair may retain multiple intersections");
            }
        }
    
        // affine同一線は無限個の交点を列挙しない。
        if (xProgram) {
            const auto xRequestA = request(session, "x", "-1", "1");
            const auto xRequestB = request(session, "x", "-1", "1");
            const auto samplesA = plot::coarseSamplePlot(
                xRequestA, plot::makeInitialPlotAnalysis(xRequestA), *xProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            const auto samplesB = plot::coarseSamplePlot(
                xRequestB, plot::makeInitialPlotAnalysis(xRequestB), *xProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            if (samplesA && samplesB) {
                const auto intersections = plot::refineCurveIntersections(
                    std::vector<plot::PlotProgram>{*xProgram.program, *xProgram.program},
                    std::vector<plot::SampledCurve>{*samplesA.curve, *samplesB.curve},
                    mathematics::defaultAngleSemantics());
                tests.expect(intersections && intersections.statistics.intersections == 0
                        && intersections.statistics.coincidentAffinePairs == 1,
                    "PlotIntersections: coincident affine lines are recorded without enumerating points");
            }
        }
    
        // domain分割された曲線同士ではpoleを跨いだ偽交点を作らない。
        const auto reciprocalProgram = compile(session, request(session, "1/(x-1)", "0", "2"));
        if (reciprocalProgram && zeroProgram) {
            const auto reciprocalRequest = request(session, "1/(x-1)", "0", "2");
            const auto zeroRequest = request(session, "0", "0", "2");
            const auto reciprocalSamples = plot::coarseSamplePlot(
                reciprocalRequest, analyze(session, reciprocalRequest), *reciprocalProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            const auto zeroSamples = plot::coarseSamplePlot(
                zeroRequest, plot::makeInitialPlotAnalysis(zeroRequest), *zeroProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), options);
            if (reciprocalSamples && zeroSamples) {
                const auto intersections = plot::refineCurveIntersections(
                    std::vector<plot::PlotProgram>{*reciprocalProgram.program, *zeroProgram.program},
                    std::vector<plot::SampledCurve>{*reciprocalSamples.curve, *zeroSamples.curve},
                    mathematics::defaultAngleSemantics());
                tests.expect(intersections && reciprocalSamples.curve->segments.size() == 2
                        && intersections.statistics.intersections == 0,
                    "PlotIntersections: domain partitions prevent false roots across a pole");
            }
        }
    
        // 近接した二つの横切り根も別anchorとして残す。
        const auto closeRootsProgram = compile(session, request(session, "x^2-1/10000", "-1/10", "1/10"));
        if (closeRootsProgram && zeroProgram) {
            const auto closeRootsRequest = request(session, "x^2-1/10000", "-1/10", "1/10");
            const auto zeroRequest = request(session, "0", "-1/10", "1/10");
            plot::CoarseSamplingOptions closeRootOptions{96, 9, 64};
            const auto rootsSamples = plot::coarseSamplePlot(
                closeRootsRequest, plot::makeInitialPlotAnalysis(closeRootsRequest), *closeRootsProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), closeRootOptions);
            const auto zeroSamples = plot::coarseSamplePlot(
                zeroRequest, plot::makeInitialPlotAnalysis(zeroRequest), *zeroProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), closeRootOptions);
            if (rootsSamples && zeroSamples) {
                const auto intersections = plot::refineCurveIntersections(
                    std::vector<plot::PlotProgram>{*closeRootsProgram.program, *zeroProgram.program},
                    std::vector<plot::SampledCurve>{*rootsSamples.curve, *zeroSamples.curve},
                    mathematics::defaultAngleSemantics());
                tests.expect(intersections && intersections.statistics.intersections == 2
                        && plot::collectPlotAnchors(*intersections.curves).anchors.size() == 2,
                    "PlotIntersections: nearby roots remain distinct shared anchors");
                plot::CurveIntersectionOptions capped;
                capped.maxIntersections = 1;
                const auto truncated = plot::refineCurveIntersections(
                    std::vector<plot::PlotProgram>{*closeRootsProgram.program, *zeroProgram.program},
                    std::vector<plot::SampledCurve>{*rootsSamples.curve, *zeroSamples.curve},
                    mathematics::defaultAngleSemantics(), capped);
                tests.expect(truncated && truncated.statistics.intersections == 1
                        && truncated.statistics.truncatedPairs == 1,
                    "PlotIntersections: semantic anchor saturation truncates the pair without failing rendering");
            }
        }
    
        // 接触候補の探索自体にも独立したresource boundを持たせる。
        if (shiftedSquareProgram && zeroProgram) {
            const auto squareRequest = request(session, "(x-1/3)^2", "-1", "1");
            const auto zeroRequest = request(session, "0", "-1", "1");
            plot::CoarseSamplingOptions tangencyOptions{96, 6, 64};
            const auto squareSamples = plot::coarseSamplePlot(
                squareRequest, plot::makeInitialPlotAnalysis(squareRequest), *shiftedSquareProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), tangencyOptions);
            const auto zeroSamples = plot::coarseSamplePlot(
                zeroRequest, plot::makeInitialPlotAnalysis(zeroRequest), *zeroProgram.program,
                session.builtinRegistry(), session.mathRegistry(),
                mathematics::defaultAngleSemantics(), tangencyOptions);
            if (squareSamples && zeroSamples) {
                plot::CurveIntersectionOptions bounded;
                bounded.maxTangencyCandidates = 1;
                bounded.maxIntersections = 1;
                const auto intersections = plot::refineCurveIntersections(
                    std::vector<plot::PlotProgram>{*shiftedSquareProgram.program, *zeroProgram.program},
                    std::vector<plot::SampledCurve>{*squareSamples.curve, *zeroSamples.curve},
                    mathematics::defaultAngleSemantics(), bounded);
                tests.expect(intersections && intersections.statistics.intersections == 1,
                    "PlotIntersections: bounded tangency refinement still accepts one valid candidate");
            }
        }
    };
    runIntersections();
}

} // namespace mmcal::tests
