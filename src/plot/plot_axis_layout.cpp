#include "plot_axis_layout.hpp"

#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"
#include "numeric/rounding_mode.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace mmcal::plot {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;
using numeric::RoundingMode;

[[nodiscard]] BigFloat zero(std::size_t precisionBits) {
    return BigFloat::fromBigInt(BigInt{0}, precisionBits, RoundingMode::NearestEven);
}

[[nodiscard]] std::optional<double> log10Abs(const BigFloat& value) {
    if (value.isZero())
        return std::nullopt;
    std::string digits = value.significand().abs().toString();
    const auto prefixLength = std::min<std::size_t>(digits.size(), 16);
    const double prefix = std::stod(digits.substr(0, prefixLength));
    const double mantissa = prefix / std::pow(10.0, static_cast<double>(prefixLength - 1));
    return static_cast<double>(digits.size() - 1) + std::log10(mantissa)
        + static_cast<double>(value.exponent()) * std::log10(2.0);
}

[[nodiscard]] std::optional<Rational> decimalStep(
    const BigFloat& minimum,
    const BigFloat& maximum,
    std::size_t targetTickCount,
    const PlotAxisLayoutOptions& options,
    std::size_t precisionBits) {
    const auto width = numeric::subtract(maximum, minimum, precisionBits, RoundingMode::NearestEven);
    const auto widthLog10 = log10Abs(width);
    if (!widthLog10 || targetTickCount == 0)
        return std::nullopt;

    const double rawLog10 = *widthLog10 - std::log10(static_cast<double>(targetTickCount));
    if (!std::isfinite(rawLog10)
        || std::abs(rawLog10) > static_cast<double>(options.maxDecimalExponentMagnitude) + 2.0)
        return std::nullopt;
    const auto exponent = static_cast<std::int64_t>(std::floor(rawLog10));
    if (exponent > static_cast<std::int64_t>(options.maxDecimalExponentMagnitude)
        || exponent < -static_cast<std::int64_t>(options.maxDecimalExponentMagnitude))
        return std::nullopt;

    const double mantissa = std::pow(10.0, rawLog10 - static_cast<double>(exponent));
    std::int64_t nice = mantissa <= 1.0 ? 1 : mantissa <= 2.0 ? 2 : mantissa <= 5.0 ? 5 : 10;
    std::int64_t stepExponent = exponent;
    if (nice == 10) {
        nice = 1;
        ++stepExponent;
    }
    if (stepExponent > static_cast<std::int64_t>(options.maxDecimalExponentMagnitude)
        || stepExponent < -static_cast<std::int64_t>(options.maxDecimalExponentMagnitude))
        return std::nullopt;
    if (stepExponent >= 0)
        return Rational{BigInt{nice} * numeric::pow(
            BigInt{10}, static_cast<std::uint64_t>(stepExponent))};
    return Rational{
        BigInt{nice}, numeric::pow(BigInt{10}, static_cast<std::uint64_t>(-stepExponent))};
}

[[nodiscard]] bool rationalConversionWithinBound(
    const BigFloat& value,
    std::size_t maxBits) noexcept {
    const auto exponent = value.exponent();
    const auto extra = exponent < 0
        ? static_cast<std::uint64_t>(-(exponent + 1)) + 1
        : static_cast<std::uint64_t>(exponent);
    return extra <= maxBits && value.significand().bitLength() <= maxBits - static_cast<std::size_t>(extra);
}

[[nodiscard]] BigInt ceilRational(const Rational& value) {
    const auto division = numeric::divmod(value.numerator(), value.denominator());
    if (division.remainder.isZero() || value.numerator().isNegative())
        return division.quotient;
    return division.quotient + BigInt{1};
}

struct TickGenerationResult final {
    bool resourceLimited = false;
    std::vector<Rational> values;
};

[[nodiscard]] TickGenerationResult generateTicks(
    const BigFloat& minimum,
    const BigFloat& maximum,
    double spanMm,
    const PlotAxisLayoutOptions& options,
    std::size_t precisionBits) {
    TickGenerationResult result;
    const auto targetCount = std::clamp<std::size_t>(
        static_cast<std::size_t>(std::max(1.0, std::round(spanMm / options.targetMajorTickSpacingMm))),
        options.minimumMajorTicks, options.maximumMajorTicks);
    const auto step = decimalStep(minimum, maximum, targetCount, options, precisionBits);
    if (!step) {
        result.resourceLimited = true;
        return result;
    }
    if (!rationalConversionWithinBound(minimum, options.maxTickIndexBits)
        || !rationalConversionWithinBound(maximum, options.maxTickIndexBits)) {
        result.resourceLimited = true;
        return result;
    }

    const Rational minExact = minimum.toRational();
    const Rational maxExact = maximum.toRational();
    BigInt index = ceilRational(minExact / *step);
    const std::size_t hardLimit = options.maximumMajorTicks + 2;
    while (result.values.size() < hardLimit) {
        const Rational tickExact = Rational{index} * *step;
        if (tickExact > maxExact)
            break;
        result.values.push_back(tickExact);
        index += BigInt{1};
    }
    if (result.values.size() == hardLimit) {
        result.resourceLimited = true;
        result.values.clear();
    }
    return result;
}

[[nodiscard]] PlotAxisPlacement axisPlacement(
    const BigFloat& minimum,
    const BigFloat& maximum,
    const BigFloat& zeroValue) noexcept {
    if (minimum <= zeroValue && zeroValue <= maximum)
        return PlotAxisPlacement::CrossZero;
    return zeroValue < minimum ? PlotAxisPlacement::MinimumEdge : PlotAxisPlacement::MaximumEdge;
}

} // namespace

PlotAxisLayoutResult layoutPlotAxes(
    const PlotViewTransform& transform,
    const PlotAxisLayoutOptions& options) {
    PlotAxisLayoutResult result;
    if ((options.generateMajorTicks && !(options.targetMajorTickSpacingMm > 0.0))
        || options.minimumMajorTicks == 0
        || options.minimumMajorTicks > options.maximumMajorTicks
        || options.maximumMajorTicks > 128
        || options.maxDecimalExponentMagnitude == 0
        || options.maxTickIndexBits < 64)
        return result;

    const auto precisionBits = transform.precisionBits();
    const auto zeroValue = zero(precisionBits);
    const auto xPlacement = axisPlacement(transform.yMinimum(), transform.yMaximum(), zeroValue);
    const auto yPlacement = axisPlacement(transform.xMinimum(), transform.xMaximum(), zeroValue);

    double xAxisPositionMm = 0.0;
    if (xPlacement == PlotAxisPlacement::CrossZero) {
        const auto mapped = transform.mapY(zeroValue);
        if (!mapped) {
            result.status = PlotAxisLayoutStatus::MappingFailed;
            return result;
        }
        xAxisPositionMm = *mapped;
    }
    else if (xPlacement == PlotAxisPlacement::MaximumEdge)
        xAxisPositionMm = transform.viewport().heightMm;

    double yAxisPositionMm = 0.0;
    if (yPlacement == PlotAxisPlacement::CrossZero) {
        const auto mapped = transform.mapX(zeroValue);
        if (!mapped) {
            result.status = PlotAxisLayoutStatus::MappingFailed;
            return result;
        }
        yAxisPositionMm = *mapped;
    }
    else if (yPlacement == PlotAxisPlacement::MaximumEdge)
        yAxisPositionMm = transform.viewport().widthMm;

    PlotAxesLayout layout{
        PlotAxisLayout{PlotAxisOrientation::Horizontal, xPlacement, xAxisPositionMm, {}},
        PlotAxisLayout{PlotAxisOrientation::Vertical, yPlacement, yAxisPositionMm, {}}};
    if (!options.generateMajorTicks) {
        result.status = PlotAxisLayoutStatus::Success;
        result.layout = std::move(layout);
        return result;
    }

    const auto xValues = generateTicks(
        transform.xMinimum(), transform.xMaximum(), transform.viewport().widthMm,
        options, precisionBits);
    const auto yValues = generateTicks(
        transform.yMinimum(), transform.yMaximum(), transform.viewport().heightMm,
        options, precisionBits);
    if (xValues.resourceLimited || yValues.resourceLimited) {
        result.status = PlotAxisLayoutStatus::ResourceLimit;
        return result;
    }

    layout.xAxis.majorTicks.reserve(xValues.values.size());
    for (const auto& exactValue : xValues.values) {
        const auto value = BigFloat::fromRational(
            exactValue, precisionBits, RoundingMode::NearestEven);
        const auto position = transform.mapX(value);
        if (!position) {
            result.status = PlotAxisLayoutStatus::MappingFailed;
            return result;
        }
        layout.xAxis.majorTicks.push_back(PlotTickLayout{value, *position, exactValue});
    }
    layout.yAxis.majorTicks.reserve(yValues.values.size());
    for (const auto& exactValue : yValues.values) {
        const auto value = BigFloat::fromRational(
            exactValue, precisionBits, RoundingMode::NearestEven);
        const auto position = transform.mapY(value);
        if (!position) {
            result.status = PlotAxisLayoutStatus::MappingFailed;
            return result;
        }
        layout.yAxis.majorTicks.push_back(PlotTickLayout{value, *position, exactValue});
    }

    result.status = PlotAxisLayoutStatus::Success;
    result.layout = std::move(layout);
    return result;
}

} // namespace mmcal::plot
