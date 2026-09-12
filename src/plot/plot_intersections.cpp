// 複数Plot曲線の交点anchor refinement
#include "plot_intersections.hpp"

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

[[nodiscard]] int signOf(const BigFloat& value) noexcept {
    if (value.isZero())
        return 0;
    return value.isNegative() ? -1 : 1;
}

[[nodiscard]] BigFloat absolute(BigFloat value) {
    return value.isNegative() ? -value : value;
}

[[nodiscard]] BigFloat midpoint(const BigFloat& lhs, const BigFloat& rhs, std::size_t bits) {
    const auto sum = numeric::add(lhs, rhs, bits, RoundingMode::NearestEven);
    const auto two = BigFloat::fromBigInt(BigInt{2}, bits, RoundingMode::NearestEven);
    return numeric::divide(sum, two, bits, RoundingMode::NearestEven);
}

struct DifferenceValue final {
    BigFloat x;
    BigFloat lhs;
    BigFloat rhs;
    BigFloat difference;
};

[[nodiscard]] std::optional<DifferenceValue> evaluateDifference(
    BigFloatPlotExecutor& lhsExecutor,
    BigFloatPlotExecutor& rhsExecutor,
    const BigFloat& x,
    std::size_t bits,
    CurveIntersectionStatistics& statistics) {
    ++statistics.evaluations;
    const auto lhs = lhsExecutor.evaluate(x);
    ++statistics.evaluations;
    const auto rhs = rhsExecutor.evaluate(x);
    if (!lhs.finite() || !rhs.finite())
        return std::nullopt;
    return DifferenceValue{
        x, lhs.value, rhs.value,
        numeric::subtract(lhs.value, rhs.value, bits, RoundingMode::NearestEven)};
}

[[nodiscard]] std::optional<BigFloat> sampledValueAt(
    const SampledCurveSegment& segment,
    const BigFloat& x) {
    const auto it = std::lower_bound(
        segment.samples.begin(), segment.samples.end(), x,
        [](const PlotSample& sample, const BigFloat& value) { return sample.x < value; });
    if (it == segment.samples.end() || it->x != x || !it->finite())
        return std::nullopt;
    return it->y;
}

[[nodiscard]] std::optional<DifferenceValue> evaluateDifferenceCached(
    const SampledCurveSegment& lhsSegment,
    const SampledCurveSegment& rhsSegment,
    BigFloatPlotExecutor& lhsExecutor,
    BigFloatPlotExecutor& rhsExecutor,
    const BigFloat& x,
    std::size_t bits,
    CurveIntersectionStatistics& statistics) {
    std::optional<BigFloat> lhsValue = sampledValueAt(lhsSegment, x);
    if (!lhsValue) {
        ++statistics.evaluations;
        const auto lhs = lhsExecutor.evaluate(x);
        if (!lhs.finite())
            return std::nullopt;
        lhsValue = lhs.value;
    }

    std::optional<BigFloat> rhsValue = sampledValueAt(rhsSegment, x);
    if (!rhsValue) {
        ++statistics.evaluations;
        const auto rhs = rhsExecutor.evaluate(x);
        if (!rhs.finite())
            return std::nullopt;
        rhsValue = rhs.value;
    }

    return DifferenceValue{
        x, *lhsValue, *rhsValue,
        numeric::subtract(*lhsValue, *rhsValue, bits, RoundingMode::NearestEven)};
}

[[nodiscard]] bool atLayoutResolution(
    const BigFloat& lhs,
    const BigFloat& rhs,
    const PlotViewTransform* transform,
    double minimumSpanMm) {
    if (!transform)
        return false;
    const auto left = transform->mapX(lhs);
    const auto right = transform->mapX(rhs);
    return left && right && std::abs(*right - *left) <= minimumSpanMm;
}

[[nodiscard]] bool safeToSnapIntersection(
    const DifferenceValue& root,
    const PlotViewTransform* transform,
    double maximumDistanceMm) {
    if (!transform)
        return true;
    const auto lhsY = transform->mapY(root.lhs);
    const auto rhsY = transform->mapY(root.rhs);
    return lhsY && rhsY && std::abs(*lhsY - *rhsY) <= maximumDistanceMm;
}

void insertAnchor(
    SampledCurveSegment& segment,
    const BigFloat& x,
    const BigFloat& y) {
    for (auto& sample : segment.samples) {
        if (sample.x == x) {
            // 交点は共有semantic objectへ昇格する前提なので，両Path上の座標も一致させる。
            sample.y = y;
            sample.status = PlotNumericStatus::Finite;
            sample.tags |= PlotPointTag::CurveIntersection;
            return;
        }
    }
    PlotSample sample{x, y, PlotNumericStatus::Finite, PlotPointTag::CurveIntersection};
    const auto position = std::lower_bound(
        segment.samples.begin(), segment.samples.end(), sample.x,
        [](const PlotSample& existing, const BigFloat& value) { return existing.x < value; });
    segment.samples.insert(position, std::move(sample));
}

[[nodiscard]] std::vector<BigFloat> mergedGrid(
    const SampledCurveSegment& lhs,
    const SampledCurveSegment& rhs) {
    std::vector<BigFloat> xs;
    if (lhs.samples.empty() || rhs.samples.empty())
        return xs;
    const BigFloat lower = lhs.samples.front().x < rhs.samples.front().x
        ? rhs.samples.front().x : lhs.samples.front().x;
    const BigFloat upper = lhs.samples.back().x < rhs.samples.back().x
        ? lhs.samples.back().x : rhs.samples.back().x;
    if (upper < lower)
        return xs;
    const auto append = [&](const SampledCurveSegment& segment) {
        for (const auto& sample : segment.samples)
            if (!(sample.x < lower) && !(upper < sample.x))
                xs.push_back(sample.x);
    };
    append(lhs);
    append(rhs);
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    return xs;
}

[[nodiscard]] DifferenceValue refineBracket(
    DifferenceValue left,
    DifferenceValue right,
    BigFloatPlotExecutor& lhsExecutor,
    BigFloatPlotExecutor& rhsExecutor,
    std::size_t bits,
    const PlotViewTransform* transform,
    const CurveIntersectionOptions& options,
    CurveIntersectionStatistics& statistics) {
    DifferenceValue best = left;
    for (std::size_t i = 0; i < options.refinementIterations; ++i) {
        if (atLayoutResolution(
                left.x, right.x, transform, options.minimumRefinementSpanMm))
            break;
        BigFloat x = midpoint(left.x, right.x, bits);
        if (x == left.x || x == right.x)
            break;
        auto middle = evaluateDifference(lhsExecutor, rhsExecutor, x, bits, statistics);
        if (!middle)
            break;
        best = *middle;
        const int sign = signOf(middle->difference);
        if (sign == 0)
            return *middle;
        if (sign == signOf(left.difference))
            left = std::move(*middle);
        else
            right = std::move(*middle);
    }
    BigFloat x = midpoint(left.x, right.x, bits);
    if (x != left.x && x != right.x) {
        if (auto final = evaluateDifference(lhsExecutor, rhsExecutor, x, bits, statistics))
            return *final;
    }
    return best;
}

[[nodiscard]] bool betterTangencySample(
    const DifferenceValue& lhs,
    const DifferenceValue& rhs) {
    return absolute(lhs.difference) < absolute(rhs.difference);
}

[[nodiscard]] std::optional<DifferenceValue> refineTangencyCandidate(
    DifferenceValue left,
    DifferenceValue center,
    DifferenceValue right,
    BigFloatPlotExecutor& lhsExecutor,
    BigFloatPlotExecutor& rhsExecutor,
    std::size_t bits,
    const PlotViewTransform* transform,
    const CurveIntersectionOptions& options,
    CurveIntersectionStatistics& statistics) {
    BigFloat scale = absolute(left.difference);
    if (scale < absolute(center.difference))
        scale = absolute(center.difference);
    if (scale < absolute(right.difference))
        scale = absolute(right.difference);
    if (scale.isZero())
        return center;

    for (std::size_t i = 0; i < options.tangencyRefinementIterations; ++i) {
        if (atLayoutResolution(
                left.x, right.x, transform, options.minimumRefinementSpanMm))
            break;
        BigFloat xl = midpoint(left.x, center.x, bits);
        BigFloat xr = midpoint(center.x, right.x, bits);
        if (xl == left.x || xl == center.x || xr == center.x || xr == right.x)
            break;
        auto midLeft = evaluateDifference(lhsExecutor, rhsExecutor, xl, bits, statistics);
        auto midRight = evaluateDifference(lhsExecutor, rhsExecutor, xr, bits, statistics);
        if (!midLeft || !midRight)
            return std::nullopt;

        if (midLeft->difference.isZero())
            return midLeft;
        if (midRight->difference.isZero())
            return midRight;

        if (betterTangencySample(*midLeft, center)
            && !betterTangencySample(*midRight, *midLeft)) {
            right = std::move(center);
            center = std::move(*midLeft);
        }
        else if (betterTangencySample(*midRight, center)
            && !betterTangencySample(*midLeft, *midRight)) {
            left = std::move(center);
            center = std::move(*midRight);
        }
        else {
            left = std::move(*midLeft);
            right = std::move(*midRight);
        }
    }

    // 接触根なら局所最小値は区間縮小に対して急速に0へ落ちる。
    // 非零のnear-missを安易に交点化しないよう，初期局所scaleに対して
    // 有効precisionの半分以上を失うほど小さくなった場合だけ受理する。
    const std::size_t verificationBits = std::max<std::size_t>(8, bits / 2);
    const auto relativeTolerance = BigFloat::fromDyadic(
        BigInt{1}, -static_cast<BigFloat::exponent_type>(verificationBits),
        bits, RoundingMode::NearestEven);
    const auto tolerance = numeric::multiply(
        scale, relativeTolerance, bits, RoundingMode::NearestEven);
    if (absolute(center.difference) <= tolerance)
        return center;
    return std::nullopt;
}

[[nodiscard]] bool isTangencyCandidate(
    const DifferenceValue& left,
    const DifferenceValue& center,
    const DifferenceValue& right) {
    const int sign = signOf(center.difference);
    if (sign == 0 || signOf(left.difference) != sign || signOf(right.difference) != sign)
        return false;
    return absolute(center.difference) < absolute(left.difference)
        && absolute(center.difference) < absolute(right.difference);
}

} // namespace

CurveIntersectionResult refineCurveIntersections(
    const std::vector<PlotProgram>& programs,
    const std::vector<SampledCurve>& curves,
    mathematics::AngleSemantics angles,
    const CurveIntersectionOptions& options,
    const PlotViewTransform* transform) {
    if (programs.size() != curves.size() || programs.size() < 2
        || options.refinementIterations == 0 || options.tangencyRefinementIterations == 0
        || options.maxTangencyCandidates == 0 || options.maxIntersections == 0
        || !(options.minimumRefinementSpanMm > 0.0)
        || !(options.maximumSnapDistanceMm > 0.0))
        return {};
    const std::size_t bits = curves.front().precisionBits;
    if (bits < 8)
        return {};
    for (const auto& curve : curves)
        if (curve.precisionBits != bits)
            return {};

    std::vector<SampledCurve> result = curves;
    CurveIntersectionStatistics statistics;

    try {
        for (std::size_t i = 0; i + 1 < programs.size(); ++i) {
            for (std::size_t j = i + 1; j < programs.size(); ++j) {
                BigFloatPlotExecutor lhsExecutor{programs[i], bits, angles};
                BigFloatPlotExecutor rhsExecutor{programs[j], bits, angles};
                std::vector<BigFloat> pairRoots;
                std::size_t pairTangencyCandidates = 0;
                bool truncatePair = false;
                bool tangencyBudgetExhausted = false;

                for (std::size_t si = 0; si < result[i].segments.size() && !truncatePair; ++si) {
                    for (std::size_t sj = 0; sj < result[j].segments.size() && !truncatePair; ++sj) {
                        auto xs = mergedGrid(result[i].segments[si], result[j].segments[sj]);
                        if (xs.empty())
                            continue;

                        std::vector<std::optional<DifferenceValue>> values;
                        values.reserve(xs.size());
                        std::size_t exactZeroCount = 0;
                        for (const auto& x : xs) {
                            auto value = evaluateDifferenceCached(
                                result[i].segments[si], result[j].segments[sj],
                                lhsExecutor, rhsExecutor, x, bits, statistics);
                            if (value && value->difference.isZero())
                                ++exactZeroCount;
                            values.push_back(std::move(value));
                        }
                        const auto affineLike = [](PlotCurveGeometryKind kind) {
                            return kind == PlotCurveGeometryKind::Constant
                                || kind == PlotCurveGeometryKind::Affine;
                        };
                        if (affineLike(programs[i].geometryKind)
                            && affineLike(programs[j].geometryKind)
                            && exactZeroCount >= 2) {
                            ++statistics.coincidentAffinePairs;
                            continue;
                        }

                        const auto commit = [&](DifferenceValue root, bool tangential) -> bool {
                            if (std::find(pairRoots.begin(), pairRoots.end(), root.x) != pairRoots.end())
                                return true;
                            // mm解像度でrefinementを止めても，両函数値がまだ離れている場合がある。
                            // その状態で両Pathを平均yへ強制すると曲線そのものを歪めるため，
                            // 可視上同一点とみなせる場合だけ共有vertexへsnapする。
                            if (!safeToSnapIntersection(
                                    root, transform, options.maximumSnapDistanceMm)) {
                                ++statistics.rejectedUnsafeSnaps;
                                return true;
                            }
                            if (pairRoots.size() >= options.maxIntersections) {
                                truncatePair = true;
                                return false;
                            }
                            pairRoots.push_back(root.x);
                            const auto y = midpoint(root.lhs, root.rhs, bits);
                            insertAnchor(result[i].segments[si], root.x, y);
                            insertAnchor(result[j].segments[sj], root.x, y);
                            ++statistics.intersections;
                            if (tangential)
                                ++statistics.tangentialIntersections;
                            return true;
                        };

                        for (std::size_t k = 0; k < values.size(); ++k) {
                            if (values[k] && values[k]->difference.isZero()) {
                                if (!commit(*values[k], false))
                                    break;
                            }
                            if (k == 0 || !values[k - 1] || !values[k]
                                || values[k - 1]->difference.isZero()
                                || values[k]->difference.isZero()
                                || signOf(values[k - 1]->difference)
                                    == signOf(values[k]->difference))
                                continue;
                            auto root = refineBracket(
                                *values[k - 1], *values[k], lhsExecutor, rhsExecutor,
                                bits, transform, options, statistics);
                            if (!commit(std::move(root), false))
                                break;
                        }

                        if (truncatePair)
                            continue;

                        for (std::size_t k = 1; k + 1 < values.size(); ++k) {
                            if (tangencyBudgetExhausted)
                                break;
                            if (!values[k - 1] || !values[k] || !values[k + 1]
                                || !isTangencyCandidate(*values[k - 1], *values[k], *values[k + 1]))
                                continue;
                            if (pairTangencyCandidates >= options.maxTangencyCandidates) {
                                tangencyBudgetExhausted = true;
                                break;
                            }
                            ++pairTangencyCandidates;
                            ++statistics.tangencyCandidates;
                            auto root = refineTangencyCandidate(
                                *values[k - 1], *values[k], *values[k + 1],
                                lhsExecutor, rhsExecutor, bits, transform, options, statistics);
                            if (!root) {
                                ++statistics.rejectedTangencyCandidates;
                                continue;
                            }
                            if (!commit(std::move(*root), true))
                                break;
                        }
                    }
                }
                if (truncatePair || tangencyBudgetExhausted)
                    ++statistics.truncatedPairs;
            }
        }
    }
    catch (...) {
        return CurveIntersectionResult{
            CurveIntersectionStatus::ProgramInitializationFailed,
            std::nullopt, statistics};
    }

    return CurveIntersectionResult{
        CurveIntersectionStatus::Success, std::move(result), statistics};
}

} // namespace mmcal::plot
