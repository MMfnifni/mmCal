// Plot粗sampleから仮viewport rangeを推定する
#include "plot_range_estimator.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rounding_mode.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace mmcal::plot {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::RoundingMode;

[[nodiscard]] BigFloat absolute(BigFloat value) {
    return value.isNegative() ? -value : value;
}

[[nodiscard]] const BigFloat& coordinateValue(
    const CurveSample2D& sample,
    CurveCoordinate2D coordinate) noexcept {
    return coordinate == CurveCoordinate2D::X ? sample.x : sample.y;
}

[[nodiscard]] const BigFloat& coordinateValue(
    const PlotEndpointPoint& point,
    CurveCoordinate2D coordinate) noexcept {
    return coordinate == CurveCoordinate2D::X ? point.x : point.y;
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

[[nodiscard]] bool isKnownAsymptoticBoundary(
    const PlotAnalysis& analysis,
    const expression::Expr& position) {
    return hasLandmarkAt(analysis, position, PlotLandmarkKind::Pole)
        || hasLandmarkAt(analysis, position, PlotLandmarkKind::VerticalAsymptote);
}

void excludeRange(
    std::vector<bool>& excluded,
    std::size_t first,
    std::size_t last) {
    if (excluded.empty())
        return;
    first = std::min(first, excluded.size() - 1);
    last = std::min(last, excluded.size() - 1);
    if (first > last)
        std::swap(first, last);
    for (std::size_t i = first; i <= last; ++i)
        excluded[i] = true;
}

[[nodiscard]] std::vector<std::vector<bool>> excludedSamples(
    const PlotAnalysis& analysis,
    const SampledCurve& curve,
    CurveCoordinate2D coordinate,
    const PlotCoarseInspection* inspection,
    const PlotRangeEstimatorOptions& options) {
    std::vector<std::vector<bool>> result;
    result.reserve(curve.segments.size());
    for (const SampledCurveSegment& segment : curve.segments)
        result.emplace_back(segment.samples.size(), false);

    for (std::size_t segmentIndex = 0; segmentIndex < curve.segments.size(); ++segmentIndex) {
        const SampledCurveSegment& segment = curve.segments[segmentIndex];
        if (segment.samples.empty())
            continue;

        if (options.asymptoteGuardSamples > 0
            && segment.sourceInterval.lowerInclusion == PlotEndpointInclusion::Open
            && isKnownAsymptoticBoundary(analysis, segment.sourceInterval.lower)) {
            excludeRange(
                result[segmentIndex], 0,
                std::min(options.asymptoteGuardSamples, segment.samples.size()) - 1);
        }
        if (options.asymptoteGuardSamples > 0
            && segment.sourceInterval.upperInclusion == PlotEndpointInclusion::Open
            && isKnownAsymptoticBoundary(analysis, segment.sourceInterval.upper)) {
            const std::size_t count = std::min(options.asymptoteGuardSamples, segment.samples.size());
            excludeRange(
                result[segmentIndex], segment.samples.size() - count,
                segment.samples.size() - 1);
        }
    }

    // fixed guardだけでは1/x^4やtan(x)^nのような高次poleで，3点目以降も
    // Automatic rangeを支配し得る。既知asymptote境界に限って，残った内部sampleの
    // median absolute valueをscaleとし，その一定倍を超える連続tailを追加除外する。
    // exp等の正当な急増にはasymptote landmarkが無いため影響しない。
    if (options.asymptoteTailScaleFactor > 0) {
        const std::size_t precisionBits = std::max<std::size_t>(2, curve.precisionBits);
        const BigFloat one = BigFloat::fromBigInt(
            BigInt{1}, precisionBits, RoundingMode::NearestEven);
        const BigFloat factor = BigFloat::fromBigInt(
            BigInt{static_cast<std::int64_t>(options.asymptoteTailScaleFactor)},
            precisionBits, RoundingMode::NearestEven);

        for (std::size_t segmentIndex = 0; segmentIndex < curve.segments.size(); ++segmentIndex) {
            const SampledCurveSegment& segment = curve.segments[segmentIndex];
            if (segment.samples.empty())
                continue;
            const bool lowerAsymptote = segment.sourceInterval.lowerInclusion
                    == PlotEndpointInclusion::Open
                && isKnownAsymptoticBoundary(analysis, segment.sourceInterval.lower);
            const bool upperAsymptote = segment.sourceInterval.upperInclusion
                    == PlotEndpointInclusion::Open
                && isKnownAsymptoticBoundary(analysis, segment.sourceInterval.upper);
            if (!lowerAsymptote && !upperAsymptote)
                continue;

            std::vector<BigFloat> magnitudes;
            magnitudes.reserve(segment.samples.size());
            std::size_t includedFinite = 0;
            for (std::size_t sampleIndex = 0; sampleIndex < segment.samples.size(); ++sampleIndex) {
                const PlotSample& sample = segment.samples[sampleIndex];
                if (!sample.finite() || result[segmentIndex][sampleIndex])
                    continue;
                magnitudes.push_back(absolute(coordinateValue(sample, coordinate)));
                ++includedFinite;
            }
            if (magnitudes.size() < options.minimumIncludedFiniteSamples)
                continue;
            std::sort(magnitudes.begin(), magnitudes.end());
            BigFloat scale = magnitudes[magnitudes.size() / 2];
            if (scale < one)
                scale = one;
            const BigFloat threshold = numeric::multiply(
                scale, factor, precisionBits, RoundingMode::NearestEven);

            auto trimFromBoundary = [&](bool lower) {
                if (includedFinite <= options.minimumIncludedFiniteSamples)
                    return;
                if (lower) {
                    for (std::size_t i = 0; i < segment.samples.size(); ++i) {
                        if (result[segmentIndex][i] || !segment.samples[i].finite())
                            continue;
                        if (!(absolute(coordinateValue(segment.samples[i], coordinate)) > threshold))
                            break;
                        if (includedFinite <= options.minimumIncludedFiniteSamples)
                            break;
                        result[segmentIndex][i] = true;
                        --includedFinite;
                    }
                }
                else {
                    for (std::size_t offset = 0; offset < segment.samples.size(); ++offset) {
                        const std::size_t i = segment.samples.size() - 1 - offset;
                        if (result[segmentIndex][i] || !segment.samples[i].finite())
                            continue;
                        if (!(absolute(coordinateValue(segment.samples[i], coordinate)) > threshold))
                            break;
                        if (includedFinite <= options.minimumIncludedFiniteSamples)
                            break;
                        result[segmentIndex][i] = true;
                        --includedFinite;
                    }
                }
            };
            if (lowerAsymptote)
                trimFromBoundary(true);
            if (upperAsymptote)
                trimFromBoundary(false);
        }
    }

    // NonFiniteSampleだけは近傍もrangeから外す。RapidVariation/Oscillatoryだけでは
    // exp等の正当な急増まで切ってしまうため，ここでは除外根拠にしない。
    if (inspection) {
        for (const PlotSuspiciousSpan& span : inspection->suspiciousSpans) {
            if (span.kind != PlotSuspicionKind::NonFiniteSample
                || span.segmentIndex >= result.size()
                || result[span.segmentIndex].empty())
                continue;
            const std::size_t first = span.firstSampleIndex > options.nonFiniteGuardSamples
                ? span.firstSampleIndex - options.nonFiniteGuardSamples
                : 0;
            const std::size_t last = std::min(
                span.lastSampleIndex + options.nonFiniteGuardSamples,
                result[span.segmentIndex].size() - 1);
            excludeRange(result[span.segmentIndex], first, last);
        }
    }
    return result;
}

struct FiniteExtent final {
    BigFloat minimum;
    BigFloat maximum;
    std::size_t count = 0;
    std::size_t excludedFiniteCount = 0;
};

void includePoint(
    std::optional<FiniteExtent>& extent,
    const BigFloat& value) {
    if (!extent) {
        extent = FiniteExtent{value, value, 1, 0};
        return;
    }
    if (value < extent->minimum)
        extent->minimum = value;
    if (value > extent->maximum)
        extent->maximum = value;
    ++extent->count;
}

[[nodiscard]] std::optional<FiniteExtent> finiteExtent(
    const SampledCurve& curve,
    CurveCoordinate2D coordinate,
    const std::vector<std::vector<bool>>& excluded,
    bool honorExclusions) {
    std::optional<FiniteExtent> result;
    for (std::size_t segmentIndex = 0; segmentIndex < curve.segments.size(); ++segmentIndex) {
        const SampledCurveSegment& segment = curve.segments[segmentIndex];
        for (std::size_t sampleIndex = 0; sampleIndex < segment.samples.size(); ++sampleIndex) {
            const PlotSample& sample = segment.samples[sampleIndex];
            if (!sample.finite())
                continue;
            const bool isExcluded = honorExclusions
                && segmentIndex < excluded.size()
                && sampleIndex < excluded[segmentIndex].size()
                && excluded[segmentIndex][sampleIndex];
            if (isExcluded) {
                if (result)
                    ++result->excludedFiniteCount;
                continue;
            }
            includePoint(result, coordinateValue(sample, coordinate));
        }

        if (segment.lowerEndpointPoint)
            includePoint(result, coordinateValue(*segment.lowerEndpointPoint, coordinate));
        if (segment.upperEndpointPoint)
            includePoint(result, coordinateValue(*segment.upperEndpointPoint, coordinate));
    }

    if (result && honorExclusions) {
        // 最初のincluded sampleより前に除外されたfinite sampleも数える。
        std::size_t excludedCount = 0;
        for (std::size_t segmentIndex = 0; segmentIndex < curve.segments.size(); ++segmentIndex) {
            const SampledCurveSegment& segment = curve.segments[segmentIndex];
            for (std::size_t sampleIndex = 0; sampleIndex < segment.samples.size(); ++sampleIndex) {
                if (segment.samples[sampleIndex].finite()
                    && segmentIndex < excluded.size()
                    && sampleIndex < excluded[segmentIndex].size()
                    && excluded[segmentIndex][sampleIndex])
                    ++excludedCount;
            }
        }
        result->excludedFiniteCount = excludedCount;
    }
    return result;
}

[[nodiscard]] BigFloat paddingFor(
    const BigFloat& minimum,
    const BigFloat& maximum,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options,
    bool& flatExpanded) {
    const BigFloat numerator = BigFloat::fromBigInt(
        BigInt{static_cast<std::int64_t>(options.paddingNumerator)},
        precisionBits, RoundingMode::NearestEven);
    const BigFloat denominator = BigFloat::fromBigInt(
        BigInt{static_cast<std::int64_t>(std::max<std::size_t>(1, options.paddingDenominator))},
        precisionBits, RoundingMode::NearestEven);

    BigFloat span = numeric::subtract(
        maximum, minimum, precisionBits, RoundingMode::NearestEven);
    if (span.isZero()) {
        flatExpanded = true;
        const BigFloat one = BigFloat::fromBigInt(
            BigInt{1}, precisionBits, RoundingMode::NearestEven);
        BigFloat scale = absolute(minimum);
        if (scale < one)
            scale = one;
        span = scale;
    }
    const BigFloat fraction = numeric::divide(
        numerator, denominator, precisionBits, RoundingMode::NearestEven);
    return numeric::multiply(span, fraction, precisionBits, RoundingMode::NearestEven);
}

} // namespace

std::optional<CurveRangeEstimate1D> makeCurveRangeEstimate(
    const BigFloat& minimum,
    const BigFloat& maximum,
    std::size_t precisionBits,
    std::size_t contributingPoints,
    const PlotRangeEstimatorOptions& options) {
    if (precisionBits < 2 || maximum < minimum)
        return std::nullopt;
    bool flatExpanded = false;
    const BigFloat padding = paddingFor(
        minimum, maximum, precisionBits, options, flatExpanded);
    return CurveRangeEstimate1D{
        minimum, maximum, minimum, maximum,
        numeric::subtract(minimum, padding, precisionBits, RoundingMode::NearestEven),
        numeric::add(maximum, padding, precisionBits, RoundingMode::NearestEven),
        contributingPoints, 0, false, flatExpanded};
}

std::optional<CurveRangeEstimate1D> estimateCurveCoordinateRange(
    const SampledCurve2D& curve,
    CurveCoordinate2D coordinate,
    const PlotRangeEstimatorOptions& options) {
    if (curve.precisionBits < 2)
        return std::nullopt;

    std::optional<FiniteExtent> extent;
    const auto includeCoordinate = [&](const PlotEndpointPoint& point) {
        includePoint(extent, coordinate == CurveCoordinate2D::X ? point.x : point.y);
    };
    for (const SampledCurveSegment2D& segment : curve.segments) {
        for (const CurveSample2D& sample : segment.samples) {
            if (!sample.finite())
                continue;
            includePoint(extent, coordinate == CurveCoordinate2D::X ? sample.x : sample.y);
        }
        if (segment.lowerEndpointPoint)
            includeCoordinate(*segment.lowerEndpointPoint);
        if (segment.upperEndpointPoint)
            includeCoordinate(*segment.upperEndpointPoint);
    }
    if (!extent)
        return std::nullopt;

    const std::size_t precisionBits = std::max<std::size_t>(2, curve.precisionBits);
    bool flatExpanded = false;
    const BigFloat padding = paddingFor(
        extent->minimum, extent->maximum, precisionBits, options, flatExpanded);
    return CurveRangeEstimate1D{
        extent->minimum, extent->maximum, extent->minimum, extent->maximum,
        numeric::subtract(extent->minimum, padding, precisionBits, RoundingMode::NearestEven),
        numeric::add(extent->maximum, padding, precisionBits, RoundingMode::NearestEven),
        extent->count, 0, false, flatExpanded};
}


std::optional<CurveRangeEstimate1D> estimateAnalyzedCurveCoordinateRange(
    const PlotAnalysis& analysis,
    const SampledCurve2D& curve,
    CurveCoordinate2D coordinate,
    const PlotRangeEstimatorOptions& options) {
    if (curve.precisionBits < 2)
        return std::nullopt;

    std::vector<std::vector<bool>> none;
    none.reserve(curve.segments.size());
    for (const auto& segment : curve.segments)
        none.emplace_back(segment.samples.size(), false);
    const auto observed = finiteExtent(curve, coordinate, none, false);
    if (!observed)
        return std::nullopt;

    const auto excluded = excludedSamples(analysis, curve, coordinate, nullptr, options);
    auto extent = finiteExtent(curve, coordinate, excluded, true);
    bool boundaryTrimmed = extent && extent->excludedFiniteCount > 0;
    if (!extent || extent->count < options.minimumIncludedFiniteSamples) {
        extent = finiteExtent(curve, coordinate, excluded, false);
        boundaryTrimmed = false;
    }
    if (!extent)
        return std::nullopt;

    const std::size_t precisionBits = std::max<std::size_t>(2, curve.precisionBits);
    bool flatExpanded = false;
    const BigFloat padding = paddingFor(
        extent->minimum, extent->maximum, precisionBits, options, flatExpanded);
    return CurveRangeEstimate1D{
        observed->minimum, observed->maximum,
        extent->minimum, extent->maximum,
        numeric::subtract(extent->minimum, padding, precisionBits, RoundingMode::NearestEven),
        numeric::add(extent->maximum, padding, precisionBits, RoundingMode::NearestEven),
        extent->count,
        boundaryTrimmed ? extent->excludedFiniteCount : 0,
        boundaryTrimmed,
        flatExpanded};
}

std::optional<CurveRangeEstimate1D> combineCurveRangeEstimates(
    const std::vector<CurveRangeEstimate1D>& estimates,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options) {
    if (estimates.empty() || precisionBits < 2)
        return std::nullopt;

    CurveRangeEstimate1D result = estimates.front();
    result.includedFiniteSamples = 0;
    result.excludedFiniteSamples = 0;
    result.boundaryTrimmed = false;
    result.flatExpanded = false;

    for (const auto& estimate : estimates) {
        if (estimate.observedMinimum < result.observedMinimum)
            result.observedMinimum = estimate.observedMinimum;
        if (estimate.observedMaximum > result.observedMaximum)
            result.observedMaximum = estimate.observedMaximum;
        if (estimate.dataMinimum < result.dataMinimum)
            result.dataMinimum = estimate.dataMinimum;
        if (estimate.dataMaximum > result.dataMaximum)
            result.dataMaximum = estimate.dataMaximum;
        result.includedFiniteSamples += estimate.includedFiniteSamples;
        result.excludedFiniteSamples += estimate.excludedFiniteSamples;
        result.boundaryTrimmed = result.boundaryTrimmed || estimate.boundaryTrimmed;
    }

    bool flatExpanded = false;
    const BigFloat padding = paddingFor(
        result.dataMinimum, result.dataMaximum, precisionBits, options, flatExpanded);
    result.viewMinimum = numeric::subtract(
        result.dataMinimum, padding, precisionBits, RoundingMode::NearestEven);
    result.viewMaximum = numeric::add(
        result.dataMaximum, padding, precisionBits, RoundingMode::NearestEven);
    result.flatExpanded = flatExpanded;
    return result;
}

std::optional<PlotRangeEstimate> estimatePlotRange(
    const PlotAnalysis& analysis,
    const SampledCurve& curve,
    const PlotCoarseInspection& inspection,
    const PlotRangeEstimatorOptions& options) {
    if (!inspection.finiteRange)
        return std::nullopt;

    const auto excluded = excludedSamples(
        analysis, curve, CurveCoordinate2D::Y, &inspection, options);
    auto extent = finiteExtent(curve, CurveCoordinate2D::Y, excluded, true);
    bool boundaryTrimmed = extent && extent->excludedFiniteCount > 0;

    if (!extent || extent->count < options.minimumIncludedFiniteSamples) {
        extent = finiteExtent(curve, CurveCoordinate2D::Y, excluded, false);
        boundaryTrimmed = false;
    }
    if (!extent)
        return std::nullopt;

    const std::size_t precisionBits = std::max<std::size_t>(2, curve.precisionBits);
    bool flatExpanded = false;
    const BigFloat padding = paddingFor(
        extent->minimum, extent->maximum, precisionBits, options, flatExpanded);

    PlotRangeEstimate result{
        inspection.finiteRange->minimum,
        inspection.finiteRange->maximum,
        extent->minimum,
        extent->maximum,
        numeric::subtract(extent->minimum, padding, precisionBits, RoundingMode::NearestEven),
        numeric::add(extent->maximum, padding, precisionBits, RoundingMode::NearestEven),
        extent->count,
        boundaryTrimmed ? extent->excludedFiniteCount : 0,
        boundaryTrimmed,
        flatExpanded};
    return result;
}


std::optional<PlotRangeEstimate> combinePlotRangeEstimates(
    const std::vector<PlotRangeEstimate>& estimates,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options) {
    return combineCurveRangeEstimates(estimates, precisionBits, options);
}

} // namespace mmcal::plot
