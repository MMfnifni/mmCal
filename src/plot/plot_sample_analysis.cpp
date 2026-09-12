// Plot粗sampleのrange候補とrefinement候補抽出
#include "plot_sample_analysis.hpp"

#include "numeric/rounding_mode.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>

namespace mmcal::plot {
namespace {

using numeric::BigFloat;
using numeric::RoundingMode;

[[nodiscard]] BigFloat absolute(BigFloat value) {
    return value.isNegative() ? -value : value;
}

[[nodiscard]] int differenceSign(
    const BigFloat& lhs,
    const BigFloat& rhs,
    std::size_t precisionBits) {
    const BigFloat difference = numeric::subtract(
        rhs, lhs, precisionBits, RoundingMode::NearestEven);
    if (difference.isPositive())
        return 1;
    if (difference.isNegative())
        return -1;
    return 0;
}

[[nodiscard]] bool hasLandmarkAt(
    const PlotAnalysis& analysis,
    const expression::Expr& position,
    PlotLandmarkKind kind) {
    return std::any_of(
        analysis.landmarks.begin(), analysis.landmarks.end(),
        [&](const PlotLandmark& landmark) {
            return landmark.position == position && landmark.kind == kind;
        });
}

[[nodiscard]] bool isUnresolvedOpenBoundary(
    const PlotAnalysis& analysis,
    const expression::Expr& position) {
    if (!hasLandmarkAt(analysis, position, PlotLandmarkKind::UndefinedPoint))
        return false;
    // pole/asymptoteは既に解析的な扱いが決まっているので，ここでは未知境界に含めない。
    return !hasLandmarkAt(analysis, position, PlotLandmarkKind::Pole)
        && !hasLandmarkAt(analysis, position, PlotLandmarkKind::VerticalAsymptote);
}

void appendSuspicion(
    PlotCoarseInspection& result,
    const PlotCoarseInspectionOptions& options,
    PlotSuspiciousSpan span) {
    if (result.suspiciousSpans.size() >= options.maxSuspiciousSpans)
        return;
    if (span.firstSampleIndex > span.lastSampleIndex)
        std::swap(span.firstSampleIndex, span.lastSampleIndex);
    result.suspiciousSpans.push_back(std::move(span));
}

[[nodiscard]] std::optional<PlotFiniteRangeCandidate> finiteRange(
    const SampledCurve& curve) {
    std::optional<PlotFiniteRangeCandidate> result;
    for (const SampledCurveSegment& segment : curve.segments) {
        for (const PlotSample& sample : segment.samples) {
            if (!sample.finite())
                continue;
            if (!result) {
                result = PlotFiniteRangeCandidate{sample.y, sample.y, 1};
                continue;
            }
            if (sample.y < result->minimum)
                result->minimum = sample.y;
            if (sample.y > result->maximum)
                result->maximum = sample.y;
            ++result->finiteSampleCount;
        }
    }
    return result;
}

void inspectNonFinite(
    PlotCoarseInspection& result,
    const SampledCurveSegment& segment,
    std::size_t segmentIndex,
    const PlotCoarseInspectionOptions& options) {
    for (std::size_t i = 0; i < segment.samples.size(); ++i) {
        const PlotSample& sample = segment.samples[i];
        if (sample.finite())
            continue;
        const std::size_t first = i == 0 ? 0 : i - 1;
        const std::size_t last = std::min(i + 1, segment.samples.size() - 1);
        appendSuspicion(result, options, PlotSuspiciousSpan{
            segmentIndex, first, last,
            PlotSuspicionKind::NonFiniteSample, sample.status});
    }
}

void inspectOscillation(
    PlotCoarseInspection& result,
    const SampledCurveSegment& segment,
    std::size_t segmentIndex,
    std::size_t precisionBits,
    const PlotCoarseInspectionOptions& options) {
    if (options.minimumOscillationTurns == 0 || segment.samples.size() < 4)
        return;

    std::size_t turns = 0;
    std::optional<std::size_t> firstTurn;
    std::size_t lastTurn = 0;
    int previousSign = 0;

    for (std::size_t i = 0; i + 1 < segment.samples.size(); ++i) {
        const PlotSample& lhs = segment.samples[i];
        const PlotSample& rhs = segment.samples[i + 1];
        if (!lhs.finite() || !rhs.finite()) {
            previousSign = 0;
            continue;
        }
        const int sign = differenceSign(lhs.y, rhs.y, precisionBits);
        if (sign == 0)
            continue;
        if (previousSign != 0 && sign != previousSign) {
            ++turns;
            if (!firstTurn)
                firstTurn = i;
            lastTurn = i + 1;
        }
        previousSign = sign;
    }

    if (turns < options.minimumOscillationTurns || !firstTurn)
        return;
    const std::size_t first = *firstTurn == 0 ? 0 : *firstTurn - 1;
    const std::size_t last = std::min(lastTurn + 1, segment.samples.size() - 1);
    appendSuspicion(result, options, PlotSuspiciousSpan{
        segmentIndex, first, last,
        PlotSuspicionKind::Oscillatory, PlotNumericStatus::Finite});
}

void inspectRapidVariation(
    PlotCoarseInspection& result,
    const SampledCurveSegment& segment,
    std::size_t segmentIndex,
    std::size_t precisionBits,
    const PlotFiniteRangeCandidate& range,
    const PlotCoarseInspectionOptions& options) {
    if (segment.samples.size() < 2 || range.minimum == range.maximum)
        return;

    const BigFloat totalRange = numeric::subtract(
        range.maximum, range.minimum, precisionBits, RoundingMode::NearestEven);
    const BigFloat twiceRange = numeric::add(
        totalRange, totalRange, precisionBits, RoundingMode::NearestEven);

    for (std::size_t i = 0; i + 1 < segment.samples.size(); ++i) {
        const PlotSample& lhs = segment.samples[i];
        const PlotSample& rhs = segment.samples[i + 1];
        if (!lhs.finite() || !rhs.finite())
            continue;
        const BigFloat delta = absolute(numeric::subtract(
            rhs.y, lhs.y, precisionBits, RoundingMode::NearestEven));
        const BigFloat fourDelta = numeric::add(
            numeric::add(delta, delta, precisionBits, RoundingMode::NearestEven),
            numeric::add(delta, delta, precisionBits, RoundingMode::NearestEven),
            precisionBits, RoundingMode::NearestEven);
        // 1区間だけで全finite rangeの半分超を飛ぶ場合だけ候補化する。
        if (fourDelta > twiceRange)
            appendSuspicion(result, options, PlotSuspiciousSpan{
                segmentIndex, i, i + 1,
                PlotSuspicionKind::RapidVariation, PlotNumericStatus::Finite});
    }
}

void inspectOpenBoundaries(
    PlotCoarseInspection& result,
    const PlotAnalysis& analysis,
    const SampledCurveSegment& segment,
    std::size_t segmentIndex,
    const PlotCoarseInspectionOptions& options) {
    if (segment.samples.empty())
        return;

    if (segment.sourceInterval.lowerInclusion == PlotEndpointInclusion::Open
        && isUnresolvedOpenBoundary(analysis, segment.sourceInterval.lower)) {
        appendSuspicion(result, options, PlotSuspiciousSpan{
            segmentIndex, 0, std::min<std::size_t>(1, segment.samples.size() - 1),
            PlotSuspicionKind::UnresolvedBoundary, PlotNumericStatus::Finite});
    }
    if (segment.sourceInterval.upperInclusion == PlotEndpointInclusion::Open
        && isUnresolvedOpenBoundary(analysis, segment.sourceInterval.upper)) {
        const std::size_t last = segment.samples.size() - 1;
        appendSuspicion(result, options, PlotSuspiciousSpan{
            segmentIndex, last == 0 ? 0 : last - 1, last,
            PlotSuspicionKind::UnresolvedBoundary, PlotNumericStatus::Finite});
    }
}

} // namespace

PlotCoarseInspection inspectCoarseSamples(
    const PlotAnalysis& analysis,
    const SampledCurve& curve,
    const PlotCoarseInspectionOptions& options) {
    PlotCoarseInspection result;
    result.finiteRange = finiteRange(curve);

    for (std::size_t i = 0; i < curve.segments.size(); ++i) {
        const SampledCurveSegment& segment = curve.segments[i];
        if (segment.samples.empty())
            continue;
        inspectNonFinite(result, segment, i, options);
        inspectOpenBoundaries(result, analysis, segment, i, options);
        inspectOscillation(result, segment, i, curve.precisionBits, options);
        if (result.finiteRange)
            inspectRapidVariation(
                result, segment, i, curve.precisionBits, *result.finiteRange, options);
    }
    return result;
}

} // namespace mmcal::plot
