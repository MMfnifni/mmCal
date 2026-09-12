#include "plot_adaptive_sampling.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rounding_mode.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::plot {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::RoundingMode;

[[nodiscard]] BigFloat midpoint(const BigFloat& lhs, const BigFloat& rhs, std::size_t precisionBits) {
    const auto sum = numeric::add(lhs, rhs, precisionBits, RoundingMode::NearestEven);
    const auto two = BigFloat::fromBigInt(BigInt{2}, precisionBits, RoundingMode::NearestEven);
    return numeric::divide(sum, two, precisionBits, RoundingMode::NearestEven);
}

// layout側のdouble比率をBigFloat演算へ持ち込むため，doubleの有効桁内でdyadic化する。
// Plot layoutは元々double mmなので，ここで数学的domain判定のexactnessを緩めるものではない。
[[nodiscard]] BigFloat layoutRatio(double value, std::size_t precisionBits) {
    int exponent = 0;
    const double fraction = std::frexp(value, &exponent);
    constexpr int mantissaBits = 52;
    const auto significand = static_cast<std::int64_t>(
        std::llround(std::ldexp(fraction, mantissaBits)));
    return BigFloat::fromDyadic(
        BigInt{significand}, exponent - mantissaBits,
        precisionBits, RoundingMode::NearestEven);
}

// adaptive chord coreはparameterとgeometry x/yを分離して扱う。
// 現行Plotは t -> (t, f(t)) を供給し，将来のParametricPlotは
// t -> (x(t), y(t)) evaluatorを供給すれば同じrefinementを使える。
class Curve2DEvaluator {
public:
    virtual ~Curve2DEvaluator() = default;

    [[nodiscard]] virtual CurveSample2D evaluate(
        BigFloat parameter, PlotPointTag tags) = 0;

    // 再帰停止用のparameter spanをlayout mmで返す。通常Plotはx方向距離を
    // そのまま使い，ParametricPlotでは2D chord等の別metricを選べる。
    [[nodiscard]] virtual std::optional<double> refinementSpanMm(
        const CurveSample2D& left,
        const CurveSample2D& right,
        const PlotViewTransform& transform) const = 0;
};

class PlotGraphEvaluator final : public Curve2DEvaluator {
public:
    PlotGraphEvaluator(
        BigFloatPlotExecutor& executor,
        AdaptiveSamplingStatistics& statistics)
        : executor_(executor), statistics_(statistics) {}

    [[nodiscard]] CurveSample2D evaluate(
        BigFloat parameter, PlotPointTag tags) override {
        ++statistics_.evaluations;
        const BigFloat x = parameter;
        auto result = executor_.evaluate(parameter);
        return CurveSample2D{
            std::move(parameter), x, std::move(result.value), result.status, tags};
    }

    [[nodiscard]] std::optional<double> refinementSpanMm(
        const CurveSample2D& left,
        const CurveSample2D& right,
        const PlotViewTransform& transform) const override {
        const auto leftMm = transform.mapX(left.x);
        const auto rightMm = transform.mapX(right.x);
        if (!leftMm || !rightMm)
            return std::nullopt;
        return std::abs(*rightMm - *leftMm);
    }

private:
    BigFloatPlotExecutor& executor_;
    AdaptiveSamplingStatistics& statistics_;
};

class ParametricCurveEvaluator final : public Curve2DEvaluator {
public:
    ParametricCurveEvaluator(
        BigFloatPlotExecutor& xExecutor,
        BigFloatPlotExecutor& yExecutor,
        AdaptiveSamplingStatistics& statistics)
        : xExecutor_(xExecutor), yExecutor_(yExecutor), statistics_(statistics) {}

    [[nodiscard]] CurveSample2D evaluate(
        BigFloat parameter, PlotPointTag tags) override {
        statistics_.evaluations += 2;
        auto x = xExecutor_.evaluate(parameter);
        auto y = yExecutor_.evaluate(parameter);
        const PlotNumericStatus status = !x.finite() ? x.status : y.status;
        return CurveSample2D{
            std::move(parameter), std::move(x.value), std::move(y.value), status, tags};
    }

    [[nodiscard]] std::optional<double> refinementSpanMm(
        const CurveSample2D&,
        const CurveSample2D&,
        const PlotViewTransform&) const override {
        // Parametric curveではendpoint chordが0でも途中でloopし得る。
        // midpointを見る前の距離だけで打ち切らず，chord error判定を必ず一度通す。
        return std::nullopt;
    }

private:
    BigFloatPlotExecutor& xExecutor_;
    BigFloatPlotExecutor& yExecutor_;
    AdaptiveSamplingStatistics& statistics_;
};

[[nodiscard]] double pointSegmentDistance(
    const PlotViewPointMm& point,
    const PlotViewPointMm& a,
    const PlotViewPointMm& b) {
    const double dx = b.xMm - a.xMm;
    const double dy = b.yMm - a.yMm;
    const double length2 = dx * dx + dy * dy;
    if (length2 == 0.0)
        return std::hypot(point.xMm - a.xMm, point.yMm - a.yMm);
    const double t = std::clamp(
        ((point.xMm - a.xMm) * dx + (point.yMm - a.yMm) * dy) / length2, 0.0, 1.0);
    const double px = a.xMm + t * dx;
    const double py = a.yMm + t * dy;
    return std::hypot(point.xMm - px, point.yMm - py);
}

struct RefineContext final {
    Curve2DEvaluator& evaluator;
    const PlotViewTransform& transform;
    const AdaptiveSamplingOptions& options;
    AdaptiveSamplingStatistics& statistics;
    std::size_t precisionBits = 0;
    std::size_t totalSamples = 0;
    bool resourceLimit = false;
};

void insertOpenBoundaryApproachSample(
    std::vector<PlotSample>& samples,
    bool lowerEndpoint,
    const std::optional<BigFloat>& boundaryX,
    PlotGraphEvaluator& evaluator,
    const PlotViewTransform& transform,
    const AdaptiveSamplingOptions& options,
    AdaptiveSamplingStatistics& statistics,
    std::size_t& totalSamples,
    bool& resourceLimit) {
    if (!boundaryX || samples.empty() || !(options.openBoundaryApproachMm > 0.0))
        return;

    const PlotSample& nearest = lowerEndpoint ? samples.front() : samples.back();
    const auto boundaryMm = transform.mapX(*boundaryX);
    const auto nearestMm = transform.mapX(nearest.x);
    if (!boundaryMm || !nearestMm)
        return;
    const double distanceMm = std::abs(*nearestMm - *boundaryMm);
    if (!(distanceMm > options.openBoundaryApproachMm))
        return;
    if (totalSamples >= options.maxTotalSamples) {
        resourceLimit = true;
        return;
    }

    const double ratio = options.openBoundaryApproachMm / distanceMm;
    if (!(ratio > 0.0 && ratio < 1.0))
        return;
    const BigFloat ratioValue = layoutRatio(ratio, transform.precisionBits());
    const BigFloat delta = numeric::subtract(
        nearest.parameter, *boundaryX, transform.precisionBits(), RoundingMode::NearestEven);
    const BigFloat offset = numeric::multiply(
        delta, ratioValue, transform.precisionBits(), RoundingMode::NearestEven);
    BigFloat x = numeric::add(
        *boundaryX, offset, transform.precisionBits(), RoundingMode::NearestEven);

    if (lowerEndpoint) {
        if (!(*boundaryX < x && x < nearest.parameter))
            return;
    }
    else if (!(nearest.parameter < x && x < *boundaryX))
        return;

    CurveSample2D sample = evaluator.evaluate(
        std::move(x), PlotPointTag::Adaptive);
    if (!sample.finite())
        return;
    if (lowerEndpoint)
        samples.insert(samples.begin(), std::move(sample));
    else
        samples.push_back(std::move(sample));
    ++totalSamples;
    ++statistics.insertedSamples;
}

void refinePair(
    const CurveSample2D& left,
    const CurveSample2D& right,
    std::size_t depth,
    RefineContext& context,
    std::vector<CurveSample2D>& output) {
    if (context.resourceLimit)
        return;
    if (depth >= context.options.maxRecursion) {
        ++context.statistics.recursionLimitHits;
        return;
    }
    if (context.totalSamples >= context.options.maxTotalSamples) {
        context.resourceLimit = true;
        return;
    }

    std::optional<PlotViewPointMm> leftView;
    std::optional<PlotViewPointMm> rightView;
    if (left.finite())
        leftView = context.transform.map(left.x, left.y);
    if (right.finite())
        rightView = context.transform.map(right.x, right.y);
    const auto refinementSpan = context.evaluator.refinementSpanMm(
        left, right, context.transform);
    if (refinementSpan && *refinementSpan <= context.options.minimumSpanMm)
        return;

    BigFloat parameter = midpoint(left.parameter, right.parameter, context.precisionBits);
    if (parameter == left.parameter || parameter == right.parameter)
        return;
    CurveSample2D middle = context.evaluator.evaluate(
        std::move(parameter), PlotPointTag::Adaptive);

    bool needsRefinement = !middle.finite() || !left.finite() || !right.finite();
    if (!needsRefinement && leftView && rightView) {
        const auto middleView = context.transform.map(middle.x, middle.y);
        needsRefinement = middleView
            && pointSegmentDistance(*middleView, *leftView, *rightView)
                > context.options.chordToleranceMm;
    }

    // 中点はchord判定のため既に評価済み。誤差内でも捨てずにvertexとして保持すれば，
    // coarse密度を半分へ下げても最終表示密度を維持でき，高価な特殊函数の重複評価を減らせる。
    ++context.totalSamples;
    ++context.statistics.insertedSamples;
    if (!needsRefinement) {
        output.push_back(std::move(middle));
        return;
    }

    refinePair(left, middle, depth + 1, context, output);
    output.push_back(middle);
    refinePair(output.back(), right, depth + 1, context, output);
}

[[nodiscard]] int signOf(const BigFloat& value) noexcept {
    if (value.isZero())
        return 0;
    return value.isNegative() ? -1 : 1;
}

[[nodiscard]] const BigFloat& coordinateValue(
    const CurveSample2D& sample,
    CurveCoordinate2D coordinate) noexcept {
    return coordinate == CurveCoordinate2D::X ? sample.x : sample.y;
}

void setCoordinateZero(
    CurveSample2D& sample,
    CurveCoordinate2D coordinate,
    std::size_t precisionBits) {
    BigFloat zero = BigFloat::fromBigInt(
        BigInt{0}, precisionBits, RoundingMode::NearestEven);
    if (coordinate == CurveCoordinate2D::X)
        sample.x = std::move(zero);
    else
        sample.y = std::move(zero);
}

[[nodiscard]] PlotPointTag axisTagForCoordinate(CurveCoordinate2D coordinate) noexcept {
    return coordinate == CurveCoordinate2D::X
        ? PlotPointTag::YAxisIntercept
        : PlotPointTag::XAxisIntercept;
}

void insertParametricAxisRoots(
    std::vector<CurveSample2D>& samples,
    ParametricCurveEvaluator& evaluator,
    CurveCoordinate2D coordinate,
    std::size_t precisionBits,
    const AdaptiveSamplingOptions& options,
    AdaptiveSamplingStatistics& statistics,
    std::size_t& totalSamples,
    bool& resourceLimit) {
    if (samples.empty())
        return;

    const PlotPointTag tag = axisTagForCoordinate(coordinate);
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (!samples[i].finite() || !coordinateValue(samples[i], coordinate).isZero())
            continue;
        const bool leftZero = i > 0 && samples[i - 1].finite()
            && coordinateValue(samples[i - 1], coordinate).isZero();
        const bool rightZero = i + 1 < samples.size() && samples[i + 1].finite()
            && coordinateValue(samples[i + 1], coordinate).isZero();
        // 軸上に連続して乗る区間は離散intersectionではないのでanchor化しない。
        if (leftZero || rightZero)
            continue;
        if (!hasPlotPointTag(samples[i].tags, tag)) {
            samples[i].tags |= tag;
            if (coordinate == CurveCoordinate2D::X)
                ++statistics.yAxisIntercepts;
            else
                ++statistics.xAxisIntercepts;
        }
    }

    std::size_t i = 1;
    while (i < samples.size()) {
        if (!samples[i - 1].finite() || !samples[i].finite()) {
            ++i;
            continue;
        }
        const int leftSign = signOf(coordinateValue(samples[i - 1], coordinate));
        const int rightSign = signOf(coordinateValue(samples[i], coordinate));
        if (leftSign == 0 || rightSign == 0 || leftSign == rightSign) {
            ++i;
            continue;
        }
        if (totalSamples >= options.maxTotalSamples) {
            resourceLimit = true;
            return;
        }

        CurveSample2D left = samples[i - 1];
        CurveSample2D right = samples[i];
        CurveSample2D best = left;
        for (std::size_t iter = 0; iter < options.rootRefinementIterations; ++iter) {
            BigFloat parameter = midpoint(left.parameter, right.parameter, precisionBits);
            if (parameter == left.parameter || parameter == right.parameter)
                break;
            auto middle = evaluator.evaluate(std::move(parameter), PlotPointTag::Adaptive);
            if (!middle.finite())
                break;
            best = middle;
            const int middleSign = signOf(coordinateValue(middle, coordinate));
            if (middleSign == 0)
                break;
            if (middleSign == signOf(coordinateValue(left, coordinate)))
                left = std::move(middle);
            else
                right = std::move(middle);
        }
        if (!best.finite()) {
            ++i;
            continue;
        }
        setCoordinateZero(best, coordinate, precisionBits);
        best.tags |= tag;
        samples.insert(samples.begin() + static_cast<std::ptrdiff_t>(i), std::move(best));
        ++totalSamples;
        ++statistics.insertedSamples;
        if (coordinate == CurveCoordinate2D::X)
            ++statistics.yAxisIntercepts;
        else
            ++statistics.xAxisIntercepts;
        i += 2;
    }
}

void refineParametricCoordinateExtrema(
    std::vector<CurveSample2D>& samples,
    ParametricCurveEvaluator& evaluator,
    CurveCoordinate2D coordinate,
    std::size_t precisionBits,
    const AdaptiveSamplingOptions& options,
    AdaptiveSamplingStatistics& statistics,
    std::size_t& totalSamples,
    bool& resourceLimit) {
    if (samples.size() < 3)
        return;

    struct Candidate final {
        std::size_t centerIndex = 0;
        bool maximum = false;
    };
    std::vector<Candidate> candidates;
    for (std::size_t i = 1; i + 1 < samples.size(); ++i) {
        const auto& a = samples[i - 1];
        const auto& b = samples[i];
        const auto& c = samples[i + 1];
        if (!a.finite() || !b.finite() || !c.finite())
            continue;
        const BigFloat leftDelta = numeric::subtract(
            coordinateValue(b, coordinate), coordinateValue(a, coordinate),
            precisionBits, RoundingMode::NearestEven);
        const BigFloat rightDelta = numeric::subtract(
            coordinateValue(c, coordinate), coordinateValue(b, coordinate),
            precisionBits, RoundingMode::NearestEven);
        if (leftDelta.isPositive() && rightDelta.isNegative())
            candidates.push_back(Candidate{i, true});
        else if (leftDelta.isNegative() && rightDelta.isPositive())
            candidates.push_back(Candidate{i, false});
    }

    std::vector<CurveSample2D> anchors;
    anchors.reserve(candidates.size());
    for (const Candidate candidate : candidates) {
        CurveSample2D left = samples[candidate.centerIndex - 1];
        CurveSample2D center = samples[candidate.centerIndex];
        CurveSample2D right = samples[candidate.centerIndex + 1];
        const auto better = [coordinate, maximum = candidate.maximum](
            const CurveSample2D& lhs, const CurveSample2D& rhs) {
            const auto& l = coordinateValue(lhs, coordinate);
            const auto& r = coordinateValue(rhs, coordinate);
            return maximum ? l > r : l < r;
        };

        for (std::size_t iter = 0; iter < options.extremumRefinementIterations; ++iter) {
            if (totalSamples + 2 > options.maxTotalSamples) {
                resourceLimit = true;
                return;
            }
            BigFloat tl = midpoint(left.parameter, center.parameter, precisionBits);
            BigFloat tr = midpoint(center.parameter, right.parameter, precisionBits);
            if (tl == left.parameter || tl == center.parameter
                || tr == center.parameter || tr == right.parameter)
                break;
            auto midLeft = evaluator.evaluate(std::move(tl), PlotPointTag::Adaptive);
            auto midRight = evaluator.evaluate(std::move(tr), PlotPointTag::Adaptive);
            totalSamples += 2;
            if (!midLeft.finite() || !midRight.finite())
                break;

            if (better(midLeft, center) && !better(midRight, midLeft)) {
                right = center;
                center = std::move(midLeft);
            }
            else if (better(midRight, center) && !better(midLeft, midRight)) {
                left = center;
                center = std::move(midRight);
            }
            else {
                left = std::move(midLeft);
                right = std::move(midRight);
            }
        }
        center.tags |= PlotPointTag::LocalExtremum;
        anchors.push_back(std::move(center));
        ++statistics.localExtrema;
    }

    for (auto& anchor : anchors) {
        bool merged = false;
        for (auto& sample : samples) {
            if (sample.parameter == anchor.parameter) {
                sample.tags |= PlotPointTag::LocalExtremum;
                merged = true;
                break;
            }
        }
        if (!merged) {
            if (totalSamples >= options.maxTotalSamples) {
                resourceLimit = true;
                return;
            }
            samples.push_back(std::move(anchor));
            ++totalSamples;
            ++statistics.insertedSamples;
        }
    }
    std::sort(samples.begin(), samples.end(), [](const CurveSample2D& lhs, const CurveSample2D& rhs) {
        return lhs.parameter < rhs.parameter;
    });
}

void markOrInsertYAxis(
    std::vector<CurveSample2D>& samples,
    PlotGraphEvaluator& evaluator,
    std::size_t precisionBits,
    AdaptiveSamplingStatistics& statistics,
    std::size_t maxTotalSamples,
    std::size_t& totalSamples,
    bool& resourceLimit) {
    if (samples.empty())
        return;
    const BigFloat zero = BigFloat::fromBigInt(BigInt{0}, precisionBits, RoundingMode::NearestEven);
    for (auto& sample : samples) {
        if (sample.x == zero) {
            sample.tags |= PlotPointTag::YAxisIntercept;
            ++statistics.yAxisIntercepts;
            return;
        }
    }
    for (std::size_t i = 1; i < samples.size(); ++i) {
        if (samples[i - 1].x < zero && zero < samples[i].x) {
            if (totalSamples >= maxTotalSamples) {
                resourceLimit = true;
                return;
            }
            auto sample = evaluator.evaluate(zero, PlotPointTag::YAxisIntercept);
            samples.insert(samples.begin() + static_cast<std::ptrdiff_t>(i), std::move(sample));
            ++totalSamples;
            ++statistics.insertedSamples;
            ++statistics.yAxisIntercepts;
            return;
        }
    }
}

void insertXAxisRoots(
    std::vector<CurveSample2D>& samples,
    PlotGraphEvaluator& evaluator,
    std::size_t precisionBits,
    const AdaptiveSamplingOptions& options,
    AdaptiveSamplingStatistics& statistics,
    std::size_t& totalSamples,
    bool& resourceLimit) {
    for (auto& sample : samples) {
        if (sample.finite() && sample.y.isZero()) {
            if (!hasPlotPointTag(sample.tags, PlotPointTag::XAxisIntercept)) {
                sample.tags |= PlotPointTag::XAxisIntercept;
                ++statistics.xAxisIntercepts;
            }
        }
    }

    std::size_t i = 1;
    while (i < samples.size()) {
        if (!samples[i - 1].finite() || !samples[i].finite()
            || signOf(samples[i - 1].y) == 0 || signOf(samples[i].y) == 0
            || signOf(samples[i - 1].y) == signOf(samples[i].y)) {
            ++i;
            continue;
        }
        if (totalSamples >= options.maxTotalSamples) {
            resourceLimit = true;
            return;
        }

        CurveSample2D left = samples[i - 1];
        CurveSample2D right = samples[i];
        CurveSample2D best = left;
        for (std::size_t iter = 0; iter < options.rootRefinementIterations; ++iter) {
            BigFloat parameter = midpoint(left.parameter, right.parameter, precisionBits);
            if (parameter == left.parameter || parameter == right.parameter)
                break;
            auto middle = evaluator.evaluate(std::move(parameter), PlotPointTag::Adaptive);
            if (!middle.finite())
                break;
            best = middle;
            const int sm = signOf(middle.y);
            if (sm == 0)
                break;
            if (sm == signOf(left.y))
                left = std::move(middle);
            else
                right = std::move(middle);
        }
        if (!best.finite()) {
            ++i;
            continue;
        }
        // 符号bracketで採用したsemantic rootはPath上でも厳密にx軸へ載せる。
        best.y = BigFloat::fromBigInt(BigInt{0}, precisionBits, RoundingMode::NearestEven);
        best.tags |= PlotPointTag::XAxisIntercept;
        samples.insert(samples.begin() + static_cast<std::ptrdiff_t>(i), std::move(best));
        ++totalSamples;
        ++statistics.insertedSamples;
        ++statistics.xAxisIntercepts;
        i += 2;
    }
}

void refineLocalExtrema(
    std::vector<CurveSample2D>& samples,
    PlotGraphEvaluator& evaluator,
    std::size_t precisionBits,
    const AdaptiveSamplingOptions& options,
    AdaptiveSamplingStatistics& statistics,
    std::size_t& totalSamples,
    bool& resourceLimit) {
    if (samples.size() < 3)
        return;

    struct Candidate final {
        std::size_t centerIndex = 0;
        bool maximum = false;
    };
    std::vector<Candidate> candidates;
    for (std::size_t i = 1; i + 1 < samples.size(); ++i) {
        const auto& a = samples[i - 1];
        const auto& b = samples[i];
        const auto& c = samples[i + 1];
        if (!a.finite() || !b.finite() || !c.finite())
            continue;
        const auto leftDelta = numeric::subtract(
            b.y, a.y, precisionBits, RoundingMode::NearestEven);
        const auto rightDelta = numeric::subtract(
            c.y, b.y, precisionBits, RoundingMode::NearestEven);
        if (leftDelta.isPositive() && rightDelta.isNegative())
            candidates.push_back(Candidate{i, true});
        else if (leftDelta.isNegative() && rightDelta.isPositive())
            candidates.push_back(Candidate{i, false});
    }

    std::vector<CurveSample2D> anchors;
    anchors.reserve(candidates.size());
    for (const Candidate candidate : candidates) {
        CurveSample2D left = samples[candidate.centerIndex - 1];
        CurveSample2D center = samples[candidate.centerIndex];
        CurveSample2D right = samples[candidate.centerIndex + 1];

        auto better = [maximum = candidate.maximum](
            const CurveSample2D& lhs, const CurveSample2D& rhs) {
            return maximum ? lhs.y > rhs.y : lhs.y < rhs.y;
        };

        for (std::size_t iter = 0; iter < options.extremumRefinementIterations; ++iter) {
            if (totalSamples + 2 > options.maxTotalSamples) {
                resourceLimit = true;
                return;
            }
            BigFloat tl = midpoint(left.parameter, center.parameter, precisionBits);
            BigFloat tr = midpoint(center.parameter, right.parameter, precisionBits);
            if (tl == left.parameter || tl == center.parameter
                || tr == center.parameter || tr == right.parameter)
                break;
            auto midLeft = evaluator.evaluate(
                std::move(tl), PlotPointTag::Adaptive);
            auto midRight = evaluator.evaluate(
                std::move(tr), PlotPointTag::Adaptive);
            totalSamples += 2;
            if (!midLeft.finite() || !midRight.finite())
                break;

            if (better(midLeft, center) && !better(midRight, midLeft)) {
                right = center;
                center = std::move(midLeft);
            }
            else if (better(midRight, center) && !better(midLeft, midRight)) {
                left = center;
                center = std::move(midRight);
            }
            else {
                left = std::move(midLeft);
                right = std::move(midRight);
            }
        }
        center.tags |= PlotPointTag::LocalExtremum;
        anchors.push_back(std::move(center));
        ++statistics.localExtrema;
    }

    if (anchors.empty())
        return;
    for (auto& anchor : anchors) {
        bool merged = false;
        for (auto& sample : samples) {
            if (sample.parameter == anchor.parameter) {
                sample.tags |= PlotPointTag::LocalExtremum;
                merged = true;
                break;
            }
        }
        if (!merged) {
            if (samples.size() >= options.maxTotalSamples) {
                resourceLimit = true;
                return;
            }
            samples.push_back(std::move(anchor));
            ++statistics.insertedSamples;
        }
    }
    std::sort(samples.begin(), samples.end(), [](const CurveSample2D& lhs, const CurveSample2D& rhs) {
        return lhs.parameter < rhs.parameter;
    });
}

} // namespace

AdaptiveSamplingResult refinePlotSamples(
    const PlotProgram& program,
    const SampledCurve& coarse,
    const PlotViewTransform& transform,
    mathematics::AngleSemantics angles,
    const AdaptiveSamplingOptions& options) {
    if (coarse.precisionBits < 8 || options.maxRecursion == 0
        || options.maxTotalSamples == 0 || !(options.chordToleranceMm > 0.0)
        || !(options.minimumSpanMm > 0.0)
        || !(options.openBoundaryApproachMm > 0.0))
        return {};

    BigFloatPlotExecutor executor{program, coarse.precisionBits, angles};
    AdaptiveSamplingStatistics statistics;
    PlotGraphEvaluator evaluator{executor, statistics};
    SampledCurve refined;
    refined.precisionBits = coarse.precisionBits;
    refined.segments.reserve(coarse.segments.size());

    std::size_t totalSamples = 0;
    for (const auto& segment : coarse.segments)
        totalSamples += segment.samples.size();
    if (totalSamples > options.maxTotalSamples)
        return AdaptiveSamplingResult{AdaptiveSamplingStatus::ResourceLimit, std::nullopt, statistics};

    bool resourceLimit = false;
    for (const auto& source : coarse.segments) {
        SampledCurveSegment target = source;
        target.samples.clear();
        if (source.samples.empty()) {
            refined.segments.push_back(std::move(target));
            continue;
        }

        std::vector<CurveSample2D> seeds = source.samples;
        if (source.geometryKind == PlotSegmentGeometryKind::Polyline) {
            if (source.sourceInterval.lowerInclusion == PlotEndpointInclusion::Open
                && source.approachLowerBoundary)
                insertOpenBoundaryApproachSample(
                    seeds, true, source.lowerBoundaryParameter, evaluator, transform, options,
                    statistics, totalSamples, resourceLimit);
            if (!resourceLimit
                && source.sourceInterval.upperInclusion == PlotEndpointInclusion::Open
                && source.approachUpperBoundary)
                insertOpenBoundaryApproachSample(
                    seeds, false, source.upperBoundaryParameter, evaluator, transform, options,
                    statistics, totalSamples, resourceLimit);
        }
        if (resourceLimit) {
            refined.segments.push_back(std::move(target));
            break;
        }

        target.samples.reserve(seeds.size() * 2);
        for (std::size_t i = 0; i + 1 < seeds.size(); ++i) {
            target.samples.push_back(seeds[i]);
            if (source.geometryKind == PlotSegmentGeometryKind::StraightLine
                || source.geometryKind == PlotSegmentGeometryKind::QuadraticBezier
                || source.geometryKind == PlotSegmentGeometryKind::CubicBezier
                || source.geometryKind == PlotSegmentGeometryKind::Ellipse
                || source.geometryKind == PlotSegmentGeometryKind::EllipticArc
                || source.geometryKind == PlotSegmentGeometryKind::PiecewiseConstant)
                continue;
            RefineContext context{
                evaluator, transform, options, statistics,
                coarse.precisionBits, totalSamples, false};
            refinePair(seeds[i], seeds[i + 1], 0, context, target.samples);
            totalSamples = context.totalSamples;
            if (context.resourceLimit) {
                resourceLimit = true;
                break;
            }
        }
        if (!resourceLimit)
            target.samples.push_back(seeds.back());

        if (!resourceLimit && options.preserveAxisIntersections) {
            markOrInsertYAxis(
                target.samples, evaluator, coarse.precisionBits, statistics,
                options.maxTotalSamples, totalSamples, resourceLimit);
            if (!resourceLimit)
                insertXAxisRoots(
                    target.samples, evaluator, coarse.precisionBits, options,
                    statistics, totalSamples, resourceLimit);
        }
        if (!resourceLimit && options.preserveLocalExtrema)
            refineLocalExtrema(
                target.samples, evaluator, coarse.precisionBits, options,
                statistics, totalSamples, resourceLimit);

        refined.segments.push_back(std::move(target));
        if (resourceLimit)
            break;
    }

    if (resourceLimit)
        return AdaptiveSamplingResult{
            AdaptiveSamplingStatus::ResourceLimit, std::nullopt, statistics};
    return AdaptiveSamplingResult{
        AdaptiveSamplingStatus::Success, std::move(refined), statistics};
}

AdaptiveSamplingResult refineParametricPlotSamples(
    const PlotProgram& xProgram,
    const PlotProgram& yProgram,
    const SampledCurve2D& coarse,
    const PlotViewTransform& transform,
    mathematics::AngleSemantics angles,
    const AdaptiveSamplingOptions& options) {
    if (coarse.precisionBits < 8 || options.maxRecursion == 0
        || options.maxTotalSamples == 0 || !(options.chordToleranceMm > 0.0)
        || !(options.minimumSpanMm > 0.0)
        || !(options.openBoundaryApproachMm > 0.0))
        return {};

    std::optional<BigFloatPlotExecutor> xExecutor;
    std::optional<BigFloatPlotExecutor> yExecutor;
    try {
        xExecutor.emplace(xProgram, coarse.precisionBits, angles);
        yExecutor.emplace(yProgram, coarse.precisionBits, angles);
    }
    catch (...) {
        return AdaptiveSamplingResult{
            AdaptiveSamplingStatus::ProgramInitializationFailed, std::nullopt, {}};
    }

    AdaptiveSamplingStatistics statistics;
    ParametricCurveEvaluator evaluator{*xExecutor, *yExecutor, statistics};
    SampledCurve2D refined;
    refined.precisionBits = coarse.precisionBits;
    refined.segments.reserve(coarse.segments.size());

    std::size_t totalSamples = 0;
    for (const auto& segment : coarse.segments)
        totalSamples += segment.samples.size();
    if (totalSamples > options.maxTotalSamples)
        return AdaptiveSamplingResult{
            AdaptiveSamplingStatus::ResourceLimit, std::nullopt, statistics};

    bool resourceLimit = false;
    for (const auto& source : coarse.segments) {
        SampledCurveSegment2D target = source;
        target.samples.clear();
        if (source.samples.empty()) {
            refined.segments.push_back(std::move(target));
            continue;
        }

        target.samples.reserve(source.samples.size() * 2);
        for (std::size_t i = 0; i + 1 < source.samples.size(); ++i) {
            target.samples.push_back(source.samples[i]);
            if (source.geometryKind != PlotSegmentGeometryKind::Polyline)
                continue;
            RefineContext context{
                evaluator, transform, options, statistics,
                coarse.precisionBits, totalSamples, false};
            refinePair(source.samples[i], source.samples[i + 1], 0, context, target.samples);
            totalSamples = context.totalSamples;
            if (context.resourceLimit) {
                resourceLimit = true;
                break;
            }
        }
        if (!resourceLimit)
            target.samples.push_back(source.samples.back());

        if (!resourceLimit && options.preserveAxisIntersections) {
            insertParametricAxisRoots(
                target.samples, evaluator, CurveCoordinate2D::Y,
                coarse.precisionBits, options, statistics, totalSamples, resourceLimit);
            if (!resourceLimit)
                insertParametricAxisRoots(
                    target.samples, evaluator, CurveCoordinate2D::X,
                    coarse.precisionBits, options, statistics, totalSamples, resourceLimit);
        }
        if (!resourceLimit && options.preserveLocalExtrema) {
            refineParametricCoordinateExtrema(
                target.samples, evaluator, CurveCoordinate2D::X,
                coarse.precisionBits, options, statistics, totalSamples, resourceLimit);
            if (!resourceLimit)
                refineParametricCoordinateExtrema(
                    target.samples, evaluator, CurveCoordinate2D::Y,
                    coarse.precisionBits, options, statistics, totalSamples, resourceLimit);
        }

        refined.segments.push_back(std::move(target));
        if (resourceLimit)
            break;
    }

    if (resourceLimit)
        return AdaptiveSamplingResult{
            AdaptiveSamplingStatus::ResourceLimit, std::nullopt, statistics};
    return AdaptiveSamplingResult{
        AdaptiveSamplingStatus::Success, std::move(refined), statistics};
}

} // namespace mmcal::plot
