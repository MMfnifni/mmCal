// PlotRequest / PlotAnalysisの基礎契約
#include "plot_analysis_tests.hpp"

#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "plot/plot_analysis.hpp"
#include "test_framework.hpp"

namespace mmcal::tests {
namespace {

[[nodiscard]] expression::Expr integerExpr(long long value) {
    return expression::Expr{numeric::Number{numeric::BigInt{value}}};
}

[[nodiscard]] std::string intervalText(const plot::PlotInterval& interval) {
    std::string result = interval.lowerInclusion == plot::PlotEndpointInclusion::Closed ? "[" : "(";
    result += formatting::formatExpr(interval.lower);
    result += ", ";
    result += formatting::formatExpr(interval.upper);
    result += interval.upperInclusion == plot::PlotEndpointInclusion::Closed ? "]" : ")";
    return result;
}

[[nodiscard]] bool hasLandmark(
    const plot::PlotAnalysis& analysis,
    std::string_view position,
    plot::PlotLandmarkKind kind) {
    for (const auto& landmark : analysis.landmarks) {
        if (landmark.kind == kind
            && formatting::formatExpr(landmark.position) == position)
            return true;
    }
    return false;
}

} // namespace

void runPlotAnalysisTests(TestRunner& tests) {
    using expression::Expr;
    using expression::Symbol;
    using plot::PlotDomainCoverage;
    using plot::PlotEndpointInclusion;
    using plot::PlotLandmark;
    using plot::PlotLandmarkConfidence;
    using plot::PlotLandmarkKind;
    using plot::PlotRequest;

    const Symbol x{"x"};
    const PlotRequest request{
        Expr{x},
        x,
        integerExpr(-3),
        integerExpr(3)};

    const auto initial = plot::makeInitialPlotAnalysis(request);
    tests.expect(initial.domain.coverage == PlotDomainCoverage::Unknown,
        "PlotAnalysis: initial domain remains unproven");
    tests.expectEqual(initial.domain.intervals.size(), std::size_t{1},
        "PlotAnalysis: initial domain contains one request interval");
    tests.expect(initial.domain.intervals[0].lower == request.lower
            && initial.domain.intervals[0].upper == request.upper,
        "PlotAnalysis: initial interval preserves exact request endpoints");
    tests.expect(initial.domain.intervals[0].lowerInclusion == PlotEndpointInclusion::Closed
            && initial.domain.intervals[0].upperInclusion == PlotEndpointInclusion::Closed,
        "PlotAnalysis: initial request endpoints are closed");
    tests.expect(initial.landmarks.empty(),
        "PlotAnalysis: initial analysis has no inferred landmarks");

    const PlotLandmark pole{
        integerExpr(1),
        PlotLandmarkKind::Pole,
        PlotLandmarkConfidence::Proven,
        std::nullopt};
    const PlotLandmark candidate{
        integerExpr(2),
        PlotLandmarkKind::VerticalAsymptote,
        PlotLandmarkConfidence::Candidate,
        std::nullopt};
    tests.expect(pole.position == integerExpr(1)
            && pole.kind == PlotLandmarkKind::Pole
            && pole.confidence == PlotLandmarkConfidence::Proven,
        "PlotLandmark: preserves exact proven landmark facts");
    tests.expect(candidate.confidence == PlotLandmarkConfidence::Candidate,
        "PlotLandmark: candidate confidence remains distinct from proof");


    kernel::KernelSession session;
    const Expr xExpression = session.evaluate("x");
    const auto* infinity = session.symbolRegistry().find("Infinity");
    tests.expect(xExpression.isSymbol() && infinity != nullptr,
        "PlotAnalysis: symbolic prepass prerequisites are registered");
    if (!xExpression.isSymbol() || !infinity)
        return;

    const Symbol plotX = xExpression.asSymbol();
    const auto analyze = [&](std::string_view source, std::string_view lower, std::string_view upper) {
        const PlotRequest plotRequest{
            session.evaluate(source),
            plotX,
            session.evaluate(lower),
            session.evaluate(upper)};
        return plot::analyzePlotRequest(
            plotRequest, session.builtinRegistry(), session.mathRegistry(),
            mathematics::defaultAngleSemantics(), infinity->symbol);
    };

    const auto rational = analyze("1/(x^2-1)", "-3", "3");
    tests.expect(rational.domain.coverage == PlotDomainCoverage::Complete
            && rational.domain.intervals.size() == 3,
        "PlotAnalysis: rational poles split the real sampling domain completely");
    if (rational.domain.intervals.size() == 3) {
        tests.expectEqual(intervalText(rational.domain.intervals[0]), std::string{"[-3, -1)"},
            "PlotAnalysis: rational left segment stops before the first pole");
        tests.expectEqual(intervalText(rational.domain.intervals[1]), std::string{"(-1, 1)"},
            "PlotAnalysis: rational middle segment stays between poles");
        tests.expectEqual(intervalText(rational.domain.intervals[2]), std::string{"(1, 3]"},
            "PlotAnalysis: rational right segment starts after the second pole");
    }
    for (std::string_view point : {"-1", "1"}) {
        tests.expect(hasLandmark(rational, point, PlotLandmarkKind::Pole)
                && hasLandmark(rational, point, PlotLandmarkKind::VerticalAsymptote)
                && hasLandmark(rational, point, PlotLandmarkKind::UndefinedPoint),
            "PlotAnalysis: rational denominator zero is a proven pole and vertical asymptote");
    }

    const auto zeroPower = analyze("x^0", "-3", "3");
    tests.expect(zeroPower.domain.coverage == PlotDomainCoverage::Complete
            && zeroPower.domain.intervals.size() == 2,
        "PlotAnalysis: zero exponent preserves the base-zero hole from 0^0 semantics");
    if (zeroPower.domain.intervals.size() == 2) {
        tests.expectEqual(intervalText(zeroPower.domain.intervals[0]), std::string{"[-3, 0)"},
            "PlotAnalysis: x^0 left component stops before the undefined origin");
        tests.expectEqual(intervalText(zeroPower.domain.intervals[1]), std::string{"(0, 3]"},
            "PlotAnalysis: x^0 right component starts after the undefined origin");
    }

    const auto rationalPositivePower = analyze("x^(3/2)", "-3", "3");
    tests.expect(rationalPositivePower.domain.coverage == PlotDomainCoverage::Complete
            && rationalPositivePower.domain.intervals.size() == 1,
        "PlotAnalysis: positive noninteger principal Power keeps only the nonnegative base branch");
    if (rationalPositivePower.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(rationalPositivePower.domain.intervals[0]), std::string{"[0, 3]"},
            "PlotAnalysis: x^(3/2) includes the defined zero endpoint");

    const auto rationalNegativePower = analyze("x^(-1/2)", "-3", "3");
    tests.expect(rationalNegativePower.domain.coverage == PlotDomainCoverage::Complete
            && rationalNegativePower.domain.intervals.size() == 1,
        "PlotAnalysis: negative noninteger principal Power requires a strictly positive base");
    if (rationalNegativePower.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(rationalNegativePower.domain.intervals[0]), std::string{"(0, 3]"},
            "PlotAnalysis: x^(-1/2) excludes zero exactly");

    const auto irrationalPower = analyze("x^Pi", "-3", "3");
    tests.expect(irrationalPower.domain.coverage == PlotDomainCoverage::Complete
            && irrationalPower.domain.intervals.size() == 1,
        "PlotAnalysis: known irrational positive exponents use principal-Power branch semantics");
    if (irrationalPower.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(irrationalPower.domain.intervals[0]), std::string{"[0, 3]"},
            "PlotAnalysis: x^Pi is real on the nonnegative real base only");

    const auto algebraicIrrationalPower = analyze("x^sqrt[2]", "-3", "3");
    tests.expect(algebraicIrrationalPower.domain.coverage == PlotDomainCoverage::Complete
            && algebraicIrrationalPower.domain.intervals.size() == 1,
        "PlotAnalysis: exact algebraic noninteger exponents reuse Integer-domain proofs");

    const auto polynomialRationalPower = analyze("(x^2-1)^(3/2)", "-3", "3");
    tests.expect(polynomialRationalPower.domain.coverage == PlotDomainCoverage::Complete
            && polynomialRationalPower.domain.intervals.size() == 2,
        "PlotAnalysis: polynomial bases of noninteger Power solve the exact branch inequality");
    if (polynomialRationalPower.domain.intervals.size() == 2) {
        tests.expectEqual(intervalText(polynomialRationalPower.domain.intervals[0]), std::string{"[-3, -1]"},
            "PlotAnalysis: positive rational Power keeps the left nonnegative polynomial branch");
        tests.expectEqual(intervalText(polynomialRationalPower.domain.intervals[1]), std::string{"[1, 3]"},
            "PlotAnalysis: positive rational Power keeps the right nonnegative polynomial branch");
    }

    const auto positiveBaseVariablePower = analyze("(x^2+1)^sin[x]", "-3", "3");
    tests.expect(positiveBaseVariablePower.domain.coverage == PlotDomainCoverage::Complete
            && positiveBaseVariablePower.domain.intervals.size() == 1,
        "PlotAnalysis: provably positive bases allow arbitrary real variable exponents");

    const auto zeroBaseVariablePower = analyze("0^x", "-3", "3");
    tests.expect(zeroBaseVariablePower.domain.coverage == PlotDomainCoverage::Complete
            && zeroBaseVariablePower.domain.intervals.size() == 1,
        "PlotAnalysis: zero-base variable Power reduces exactly to exponent > 0");
    if (zeroBaseVariablePower.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(zeroBaseVariablePower.domain.intervals[0]), std::string{"(0, 3]"},
            "PlotAnalysis: 0^x excludes zero and negative exponents");

    const auto selfPowerPositiveRequest = analyze("x^x", "0", "3");
    tests.expect(selfPowerPositiveRequest.domain.coverage == PlotDomainCoverage::Complete
            && selfPowerPositiveRequest.domain.intervals.size() == 1,
        "PlotAnalysis: request bounds prove the positive real branch of x^x exactly");
    if (selfPowerPositiveRequest.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(selfPowerPositiveRequest.domain.intervals[0]), std::string{"(0, 3]"},
            "PlotAnalysis: x^x preserves the 0^0 hole on a nonnegative request");

    const auto unresolvedSelfPower = analyze("x^x", "-3", "3");
    tests.expect(unresolvedSelfPower.samplingSafety == plot::PlotSamplingSafety::UnsupportedDiscontinuity,
        "PlotAnalysis: x^x with negative bases fails safely instead of dropping isolated integer-exponent points");

    const auto squareRoot = analyze("sqrt[x]", "-3", "3");
    tests.expect(squareRoot.domain.coverage == PlotDomainCoverage::Complete
            && squareRoot.domain.intervals.size() == 1,
        "PlotAnalysis: sqrt obtains a complete clipped real domain");
    if (squareRoot.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(squareRoot.domain.intervals[0]), std::string{"[0, 3]"},
            "PlotAnalysis: sqrt keeps its exact closed branch boundary");
    tests.expect(hasLandmark(squareRoot, "0", PlotLandmarkKind::DomainBoundary)
            && hasLandmark(squareRoot, "0", PlotLandmarkKind::BranchPoint)
            && !hasLandmark(squareRoot, "0", PlotLandmarkKind::UndefinedPoint),
        "PlotAnalysis: sqrt zero is a defined branch/domain boundary");

    const auto logarithm = analyze("log[x]", "-3", "3");
    tests.expect(logarithm.domain.coverage == PlotDomainCoverage::Complete
            && logarithm.domain.intervals.size() == 1,
        "PlotAnalysis: log obtains a complete clipped real domain");
    if (logarithm.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(logarithm.domain.intervals[0]), std::string{"(0, 3]"},
            "PlotAnalysis: log excludes its exact branch boundary");
    tests.expect(hasLandmark(logarithm, "0", PlotLandmarkKind::BranchPoint)
            && hasLandmark(logarithm, "0", PlotLandmarkKind::UndefinedPoint)
            && hasLandmark(logarithm, "0", PlotLandmarkKind::VerticalAsymptote)
            && !hasLandmark(logarithm, "0", PlotLandmarkKind::Pole),
        "PlotAnalysis: log zero is a branch asymptote but not a meromorphic pole");

    const auto tangent = analyze("tan[x]", "-Pi", "Pi");
    tests.expect(tangent.domain.coverage == PlotDomainCoverage::Complete
            && tangent.domain.intervals.size() == 3,
        "PlotAnalysis: tan enumerates exact periodic poles inside an exact-angle request");
    if (tangent.domain.intervals.size() == 3) {
        tests.expectEqual(intervalText(tangent.domain.intervals[0]), std::string{"[-Pi, -Pi/2)"},
            "PlotAnalysis: tan left periodic segment is exact");
        tests.expectEqual(intervalText(tangent.domain.intervals[1]), std::string{"(-Pi/2, Pi/2)"},
            "PlotAnalysis: tan middle periodic segment is exact");
        tests.expectEqual(intervalText(tangent.domain.intervals[2]), std::string{"(Pi/2, Pi]"},
            "PlotAnalysis: tan right periodic segment is exact");
    }
    for (std::string_view point : {"-Pi/2", "Pi/2"}) {
        tests.expect(hasLandmark(tangent, point, PlotLandmarkKind::Pole)
                && hasLandmark(tangent, point, PlotLandmarkKind::VerticalAsymptote),
            "PlotAnalysis: tan periodic singularities are proven poles");
    }


    const auto affineTangent = analyze("tan[2x+Pi/4]", "-Pi", "Pi");
    tests.expect(affineTangent.domain.coverage == PlotDomainCoverage::Complete
            && affineTangent.domain.intervals.size() == 5,
        "PlotAnalysis: affine tan arguments enumerate every periodic pole in the request");
    for (std::string_view point : {"-7Pi/8", "-3Pi/8", "Pi/8", "5Pi/8"})
        tests.expect(hasLandmark(affineTangent, point, PlotLandmarkKind::Pole),
            "PlotAnalysis: affine tan pole positions remain exact");

    const auto cotangent = analyze("cot[x]", "-Pi", "Pi");
    tests.expect(cotangent.domain.coverage == PlotDomainCoverage::Complete
            && cotangent.domain.intervals.size() == 2,
        "PlotAnalysis: cot poles at request endpoints and zero split the domain safely");
    if (cotangent.domain.intervals.size() == 2) {
        tests.expectEqual(intervalText(cotangent.domain.intervals[0]), std::string{"(-Pi, 0)"},
            "PlotAnalysis: cot left segment excludes both poles");
        tests.expectEqual(intervalText(cotangent.domain.intervals[1]), std::string{"(0, Pi)"},
            "PlotAnalysis: cot right segment excludes both poles");
    }

    const auto sqrtTangent = analyze("sqrt[tan[x]]", "-Pi", "Pi");
    tests.expect(sqrtTangent.domain.coverage == PlotDomainCoverage::Complete
            && sqrtTangent.domain.intervals.size() == 3,
        "PlotAnalysis: real sqrt[tan[x]] keeps only non-negative tan branches");
    if (sqrtTangent.domain.intervals.size() == 3) {
        tests.expectEqual(intervalText(sqrtTangent.domain.intervals[0]), std::string{"[-Pi, -Pi/2)"},
            "PlotAnalysis: sqrt[tan] keeps the left positive branch");
        tests.expectEqual(intervalText(sqrtTangent.domain.intervals[1]), std::string{"[0, Pi/2)"},
            "PlotAnalysis: sqrt[tan] keeps the central positive branch");
        tests.expectEqual(intervalText(sqrtTangent.domain.intervals[2]), std::string{"[Pi, Pi]"},
            "PlotAnalysis: sqrt[tan] preserves an isolated endpoint zero");
    }

    const auto logTangent = analyze("log[tan[x]]", "-Pi", "Pi");
    tests.expect(logTangent.domain.coverage == PlotDomainCoverage::Complete
            && logTangent.domain.intervals.size() == 2,
        "PlotAnalysis: real log[tan[x]] keeps only strictly positive tan branches");
    if (logTangent.domain.intervals.size() == 2) {
        tests.expectEqual(intervalText(logTangent.domain.intervals[0]), std::string{"(-Pi, -Pi/2)"},
            "PlotAnalysis: log[tan] excludes the zero and pole on the left branch");
        tests.expectEqual(intervalText(logTangent.domain.intervals[1]), std::string{"(0, Pi/2)"},
            "PlotAnalysis: log[tan] excludes the zero and pole on the central branch");
    }

    const auto shiftedTangent = analyze("tan[x]+1", "-Pi", "Pi");
    tests.expect(shiftedTangent.samplingSafety == plot::PlotSamplingSafety::Safe
            && shiftedTangent.domain.intervals.size() == 3,
        "PlotAnalysis: periodic holes propagate through arithmetic wrappers");

    const auto wrappedTangent = analyze("atan[tan[x]]", "-Pi", "Pi");
    tests.expect(wrappedTangent.samplingSafety == plot::PlotSamplingSafety::Safe
            && wrappedTangent.domain.intervals.size() == 3,
        "PlotAnalysis: periodic holes propagate through real-safe outer functions");

    const auto unsupportedShiftedSqrt = analyze("sqrt[tan[x]+1]", "-Pi", "Pi");
    tests.expect(unsupportedShiftedSqrt.samplingSafety
            == plot::PlotSamplingSafety::UnsupportedDiscontinuity,
        "PlotAnalysis: unproved periodic branch inequalities fail safely");

    const auto reciprocalTangent = analyze("1/tan[x]", "-Pi", "Pi");
    tests.expect(reciprocalTangent.domain.coverage == PlotDomainCoverage::Complete
            && reciprocalTangent.domain.intervals.size() == 4,
        "PlotAnalysis: reciprocal tan excludes both tan zeros and tan poles exactly");
    if (reciprocalTangent.domain.intervals.size() == 4) {
        tests.expectEqual(intervalText(reciprocalTangent.domain.intervals[0]), std::string{"(-Pi, -Pi/2)"},
            "PlotAnalysis: reciprocal tan first component is exact");
        tests.expectEqual(intervalText(reciprocalTangent.domain.intervals[1]), std::string{"(-Pi/2, 0)"},
            "PlotAnalysis: reciprocal tan second component is exact");
        tests.expectEqual(intervalText(reciprocalTangent.domain.intervals[2]), std::string{"(0, Pi/2)"},
            "PlotAnalysis: reciprocal tan third component is exact");
        tests.expectEqual(intervalText(reciprocalTangent.domain.intervals[3]), std::string{"(Pi/2, Pi)"},
            "PlotAnalysis: reciprocal tan fourth component is exact");
    }

    const auto sqrtAbsTangent = analyze("sqrt[abs[tan[x]]]", "-Pi", "Pi");
    tests.expect(sqrtAbsTangent.domain.coverage == PlotDomainCoverage::Complete
            && sqrtAbsTangent.domain.intervals.size() == 3,
        "PlotAnalysis: sqrt[abs[tan]] keeps every real tan branch while preserving poles");

    const auto logAbsTangent = analyze("log[abs[tan[x]]]", "-Pi", "Pi");
    tests.expect(logAbsTangent.domain.coverage == PlotDomainCoverage::Complete
            && logAbsTangent.domain.intervals.size() == 4,
        "PlotAnalysis: log[abs[tan]] excludes both tan zeros and poles exactly");

    const auto sqrtSquaredTangent = analyze("sqrt[tan[x]^2]", "-Pi", "Pi");
    tests.expect(sqrtSquaredTangent.domain.coverage == PlotDomainCoverage::Complete
            && sqrtSquaredTangent.domain.intervals.size() == 3,
        "PlotAnalysis: sqrt of an even trig power inherits only the inner definedness holes");

    const auto logSquaredSine = analyze("log[sin[x]^2]", "-Pi", "Pi");
    tests.expect(logSquaredSine.domain.coverage == PlotDomainCoverage::Complete
            && logSquaredSine.domain.intervals.size() == 2,
        "PlotAnalysis: log of an even sine power excludes exact sine zeros");

    const auto shiftedAtanhSine = analyze("atanh[sin[x+Pi/7]]", "-Pi", "Pi");
    tests.expect(shiftedAtanhSine.domain.coverage == PlotDomainCoverage::Complete
            && shiftedAtanhSine.domain.intervals.size() == 3,
        "PlotAnalysis: shifted atanh[sin] excludes exact periodic extrema instead of relying on the sample grid");
    tests.expect(hasLandmark(shiftedAtanhSine, "-9Pi/14", PlotLandmarkKind::UndefinedPoint)
            && hasLandmark(shiftedAtanhSine, "5Pi/14", PlotLandmarkKind::UndefinedPoint),
        "PlotAnalysis: shifted atanh[sin] extrema remain exact undefined boundaries");

    const auto asinTangent = analyze("asin[tan[x]]", "-Pi", "Pi");
    tests.expect(asinTangent.domain.coverage == PlotDomainCoverage::Complete
            && asinTangent.domain.intervals.size() == 3,
        "PlotAnalysis: asin[tan] keeps exactly the |tan| <= 1 arcs");
    if (asinTangent.domain.intervals.size() == 3) {
        tests.expectEqual(intervalText(asinTangent.domain.intervals[0]), std::string{"[-Pi, -3Pi/4]"},
            "PlotAnalysis: asin[tan] left unit interval is exact");
        tests.expectEqual(intervalText(asinTangent.domain.intervals[1]), std::string{"[-Pi/4, Pi/4]"},
            "PlotAnalysis: asin[tan] central unit interval is exact");
        tests.expectEqual(intervalText(asinTangent.domain.intervals[2]), std::string{"[3Pi/4, Pi]"},
            "PlotAnalysis: asin[tan] right unit interval is exact");
    }

    const auto atanhTangent = analyze("atanh[tan[x]]", "-Pi", "Pi");
    tests.expect(atanhTangent.domain.coverage == PlotDomainCoverage::Complete
            && atanhTangent.domain.intervals.size() == 3,
        "PlotAnalysis: atanh[tan] keeps exactly the strict |tan| < 1 arcs");

    const auto asinSecant = analyze("asin[sec[x]]", "-Pi", "Pi");
    tests.expect(asinSecant.domain.coverage == PlotDomainCoverage::Complete
            && asinSecant.domain.intervals.size() == 3,
        "PlotAnalysis: asin[sec] preserves only the isolated |sec| = 1 points");
    if (asinSecant.domain.intervals.size() == 3) {
        tests.expectEqual(intervalText(asinSecant.domain.intervals[0]), std::string{"[-Pi, -Pi]"},
            "PlotAnalysis: asin[sec] preserves the left isolated endpoint");
        tests.expectEqual(intervalText(asinSecant.domain.intervals[1]), std::string{"[0, 0]"},
            "PlotAnalysis: asin[sec] preserves the central isolated point");
        tests.expectEqual(intervalText(asinSecant.domain.intervals[2]), std::string{"[Pi, Pi]"},
            "PlotAnalysis: asin[sec] preserves the right isolated endpoint");
    }

    const auto sqrtSine = analyze("sqrt[sin[x]]", "-Pi", "Pi");
    tests.expect(sqrtSine.domain.coverage == PlotDomainCoverage::Complete
            && sqrtSine.domain.intervals.size() == 2,
        "PlotAnalysis: real sqrt[sin[x]] retains positive arcs and isolated endpoint zeros");
    if (sqrtSine.domain.intervals.size() == 2) {
        tests.expectEqual(intervalText(sqrtSine.domain.intervals[0]), std::string{"[-Pi, -Pi]"},
            "PlotAnalysis: sqrt[sin] preserves the isolated left endpoint zero");
        tests.expectEqual(intervalText(sqrtSine.domain.intervals[1]), std::string{"[0, Pi]"},
            "PlotAnalysis: sqrt[sin] keeps the non-negative half-period");
    }

    const auto logSine = analyze("log[sin[x]]", "-Pi", "Pi");
    tests.expect(logSine.domain.coverage == PlotDomainCoverage::Complete
            && logSine.domain.intervals.size() == 1
            && intervalText(logSine.domain.intervals[0]) == "(0, Pi)",
        "PlotAnalysis: real log[sin[x]] keeps only the positive half-period");

    const auto reciprocalSquareRoot = analyze("1/sqrt[x]", "-1", "2");
    tests.expect(reciprocalSquareRoot.domain.coverage == PlotDomainCoverage::Complete
            && reciprocalSquareRoot.domain.intervals.size() == 1,
        "PlotAnalysis: reciprocal sqrt obtains a complete positive real domain");
    if (reciprocalSquareRoot.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(reciprocalSquareRoot.domain.intervals[0]), std::string{"(0, 2]"},
            "PlotAnalysis: reciprocal sqrt excludes its zero denominator exactly");

    const auto cbrtNested = analyze("asin[cbrt[x]]", "-2", "2");
    tests.expect(cbrtNested.domain.coverage == PlotDomainCoverage::Complete
            && cbrtNested.domain.intervals.size() == 1,
        "PlotAnalysis: cbrt nested under asin propagates both exact branch bounds");
    if (cbrtNested.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(cbrtNested.domain.intervals[0]), std::string{"[-1, 1]"},
            "PlotAnalysis: cbrt inverse maps asin domain bounds back exactly");

    const Expr asinhMinusOne = session.evaluate("asinh[-1]");
    const Expr asinhOne = session.evaluate("asinh[1]");
    const auto sinhNested = analyze("asin[sinh[x]]", "-2", "2");
    tests.expect(sinhNested.domain.coverage == PlotDomainCoverage::Complete
            && sinhNested.domain.intervals.size() == 1,
        "PlotAnalysis: globally injective monotone inner functions propagate branch bounds");
    if (sinhNested.domain.intervals.size() == 1)
        tests.expect(sinhNested.domain.intervals[0].lower == asinhMinusOne
                && sinhNested.domain.intervals[0].upper == asinhOne
                && sinhNested.domain.intervals[0].lowerInclusion == PlotEndpointInclusion::Closed
                && sinhNested.domain.intervals[0].upperInclusion == PlotEndpointInclusion::Closed,
            "PlotAnalysis: sinh bounds are inverted through MathRegistry knowledge");

    const auto strictSinhNested = analyze("atanh[sinh[x]]", "-2", "2");
    tests.expect(strictSinhNested.domain.coverage == PlotDomainCoverage::Complete
            && strictSinhNested.domain.intervals.size() == 1,
        "PlotAnalysis: strict inverse-hyperbolic branch bounds remain complete");
    if (strictSinhNested.domain.intervals.size() == 1)
        tests.expect(strictSinhNested.domain.intervals[0].lower == asinhMinusOne
                && strictSinhNested.domain.intervals[0].upper == asinhOne
                && strictSinhNested.domain.intervals[0].lowerInclusion == PlotEndpointInclusion::Open
                && strictSinhNested.domain.intervals[0].upperInclusion == PlotEndpointInclusion::Open,
            "PlotAnalysis: atanh strict bounds survive monotone inversion");

    const auto shiftedCbrt = analyze("sqrt[cbrt[x]-2]", "0", "10");
    tests.expect(shiftedCbrt.domain.coverage == PlotDomainCoverage::Complete
            && shiftedCbrt.domain.intervals.size() == 1,
        "PlotAnalysis: cbrt arbitrary rational thresholds reduce exactly");
    if (shiftedCbrt.domain.intervals.size() == 1)
        tests.expectEqual(intervalText(shiftedCbrt.domain.intervals[0]), std::string{"[8, 10]"},
            "PlotAnalysis: cbrt threshold two maps to x >= 8");

    const auto floorStep = analyze("floor[x]", "-2", "2");
    tests.expect(floorStep.samplingSafety == plot::PlotSamplingSafety::Safe
            && floorStep.domain.coverage == PlotDomainCoverage::Complete
            && floorStep.domain.intervals.size() == 5,
        "PlotAnalysis: floor affine jumps split every integer boundary including a discontinuous request endpoint");
    if (floorStep.domain.intervals.size() == 5) {
        tests.expectEqual(intervalText(floorStep.domain.intervals[0]), std::string{"[-2, -1)"},
            "PlotAnalysis: floor keeps the jump point on the right branch");
        tests.expectEqual(intervalText(floorStep.domain.intervals[1]), std::string{"[-1, 0)"},
            "PlotAnalysis: floor middle step remains disconnected");
        tests.expectEqual(intervalText(floorStep.domain.intervals[2]), std::string{"[0, 1)"},
            "PlotAnalysis: floor zero step remains disconnected");
        tests.expectEqual(intervalText(floorStep.domain.intervals[3]), std::string{"[1, 2)"},
            "PlotAnalysis: floor final interior step stays open at the discontinuous upper request endpoint");
        tests.expectEqual(intervalText(floorStep.domain.intervals[4]), std::string{"[2, 2]"},
            "PlotAnalysis: floor preserves the exact upper-endpoint value as a singleton");
    }
    for (std::string_view point : {"-1", "0", "1"})
        tests.expect(hasLandmark(floorStep, point, PlotLandmarkKind::JumpDiscontinuity),
            "PlotAnalysis: floor integer crossings are proven jump landmarks");

    const auto ceilStep = analyze("ceil[x]", "-2", "2");
    tests.expect(ceilStep.samplingSafety == plot::PlotSamplingSafety::Safe
            && ceilStep.domain.intervals.size() == 5,
        "PlotAnalysis: ceil affine jumps split every integer boundary including a discontinuous request endpoint");
    if (ceilStep.domain.intervals.size() == 5) {
        tests.expectEqual(intervalText(ceilStep.domain.intervals[0]), std::string{"[-2, -2]"},
            "PlotAnalysis: ceil preserves the exact lower-endpoint value as a singleton");
        tests.expectEqual(intervalText(ceilStep.domain.intervals[1]), std::string{"(-2, -1]"},
            "PlotAnalysis: ceil first visible step starts open after the discontinuous request endpoint");
        tests.expectEqual(intervalText(ceilStep.domain.intervals[2]), std::string{"(-1, 0]"},
            "PlotAnalysis: ceil middle step remains disconnected");
        tests.expectEqual(intervalText(ceilStep.domain.intervals[3]), std::string{"(0, 1]"},
            "PlotAnalysis: ceil zero step remains disconnected");
        tests.expectEqual(intervalText(ceilStep.domain.intervals[4]), std::string{"(1, 2]"},
            "PlotAnalysis: ceil final step preserves the request endpoint");
    }

    const auto signStep = analyze("sign[x]", "-2", "2");
    tests.expect(signStep.samplingSafety == plot::PlotSamplingSafety::Safe
            && signStep.domain.intervals.size() == 3,
        "PlotAnalysis: sign isolates its defined zero from both one-sided jumps");
    if (signStep.domain.intervals.size() == 3) {
        tests.expectEqual(intervalText(signStep.domain.intervals[0]), std::string{"[-2, 0)"},
            "PlotAnalysis: sign left branch stops before zero");
        tests.expectEqual(intervalText(signStep.domain.intervals[1]), std::string{"[0, 0]"},
            "PlotAnalysis: sign preserves the exact defined zero as a singleton");
        tests.expectEqual(intervalText(signStep.domain.intervals[2]), std::string{"(0, 2]"},
            "PlotAnalysis: sign right branch starts after zero");
    }

    const auto nonlinearFloor = analyze("floor[x^2]", "-2", "2");
    tests.expect(nonlinearFloor.samplingSafety == plot::PlotSamplingSafety::Safe
            && nonlinearFloor.domain.intervals.size() == 11,
        "PlotAnalysis: polynomial floor enumerates integer level preimages exactly");
    tests.expect(hasLandmark(nonlinearFloor, "0", PlotLandmarkKind::JumpDiscontinuity),
        "PlotAnalysis: tangential integer levels are conservatively isolated as discontinuity candidates");

    const auto polynomialSign = analyze("sign[x^2]", "-2", "2");
    tests.expect(polynomialSign.samplingSafety == plot::PlotSamplingSafety::Safe
            && polynomialSign.domain.intervals.size() == 3,
        "PlotAnalysis: polynomial sign isolates exact real zero sets");

    const auto nestedPolynomialStep = analyze("sin[floor[x^2]]", "-2", "2");
    tests.expect(nestedPolynomialStep.samplingSafety == plot::PlotSamplingSafety::Safe
            && nestedPolynomialStep.domain.intervals.size() == 11,
        "PlotAnalysis: polynomial floor cuts survive inside continuous outer functions");

    const auto floorPiRange = analyze("floor[x]", "-Pi", "Pi");
    tests.expect(floorPiRange.samplingSafety == plot::PlotSamplingSafety::Safe
            && floorPiRange.domain.intervals.size() == 8,
        "PlotAnalysis: certified endpoint bounds enumerate affine floor jumps across Pi endpoints");

    const auto periodicFloor = analyze("floor[sin[x]]", "-Pi", "Pi");
    tests.expect(periodicFloor.samplingSafety == plot::PlotSamplingSafety::Safe
            && periodicFloor.domain.intervals.size() == 5,
        "PlotAnalysis: floor[sin[x]] isolates zero jumps and the +1 extremum singleton exactly");
    if (periodicFloor.domain.intervals.size() == 5) {
        tests.expectEqual(intervalText(periodicFloor.domain.intervals[0]), std::string{"[-Pi, -Pi]"},
            "PlotAnalysis: floor[sin[x]] keeps the exact endpoint zero as a singleton");
        tests.expectEqual(intervalText(periodicFloor.domain.intervals[1]), std::string{"(-Pi, 0)"},
            "PlotAnalysis: floor[sin[x]] negative branch starts open after the endpoint jump");
        tests.expectEqual(intervalText(periodicFloor.domain.intervals[2]), std::string{"[0, Pi/2)"},
            "PlotAnalysis: floor[sin[x]] positive branch stops before the isolated maximum");
        tests.expectEqual(intervalText(periodicFloor.domain.intervals[3]), std::string{"[Pi/2, Pi/2]"},
            "PlotAnalysis: floor[sin[x]] preserves the isolated value one at Pi/2");
        tests.expectEqual(intervalText(periodicFloor.domain.intervals[4]), std::string{"(Pi/2, Pi]"},
            "PlotAnalysis: floor[sin[x]] resumes after the isolated maximum");
    }

    const auto periodicCeil = analyze("ceil[sin[x]]", "-Pi", "Pi");
    tests.expect(periodicCeil.samplingSafety == plot::PlotSamplingSafety::Safe
            && periodicCeil.domain.intervals.size() == 5,
        "PlotAnalysis: ceil[sin[x]] isolates the -1 extremum and zero jump exactly");

    const auto periodicSign = analyze("sign[sin[x]]", "-Pi", "Pi");
    tests.expect(periodicSign.samplingSafety == plot::PlotSamplingSafety::Safe
            && periodicSign.domain.intervals.size() == 5,
        "PlotAnalysis: sign[sin[x]] isolates each interior zero as a singleton");

    const auto periodicCosFloor = analyze("floor[cos[x]]", "-Pi", "Pi");
    tests.expect(periodicCosFloor.samplingSafety == plot::PlotSamplingSafety::Safe
            && periodicCosFloor.domain.intervals.size() == 5,
        "PlotAnalysis: floor[cos[x]] handles zero crossings and the isolated +1 maximum");

    const auto unsupportedPeriodicStep = analyze("floor[sin[x^2]]", "-Pi", "Pi");
    tests.expect(unsupportedPeriodicStep.samplingSafety
            == plot::PlotSamplingSafety::UnsupportedDiscontinuity,
        "PlotAnalysis: non-affine periodic step phases still fail safely instead of drawing false connections");

    const auto oscillatory = analyze("sin[1/x]", "-1", "1");
    tests.expect(oscillatory.domain.coverage == PlotDomainCoverage::Complete
            && oscillatory.domain.intervals.size() == 2,
        "PlotAnalysis: sin(1/x) splits exactly at its excluded origin");
    if (oscillatory.domain.intervals.size() == 2) {
        tests.expectEqual(intervalText(oscillatory.domain.intervals[0]), std::string{"[-1, 0)"},
            "PlotAnalysis: sin(1/x) preserves the left punctured interval");
        tests.expectEqual(intervalText(oscillatory.domain.intervals[1]), std::string{"(0, 1]"},
            "PlotAnalysis: sin(1/x) preserves the right punctured interval");
    }
    tests.expect(hasLandmark(oscillatory, "0", PlotLandmarkKind::DomainBoundary)
            && hasLandmark(oscillatory, "0", PlotLandmarkKind::UndefinedPoint)
            && !hasLandmark(oscillatory, "0", PlotLandmarkKind::Pole)
            && !hasLandmark(oscillatory, "0", PlotLandmarkKind::VerticalAsymptote),
        "PlotAnalysis: sin(1/x) origin is undefined but not a pole/asymptote");
}

} // namespace mmcal::tests
