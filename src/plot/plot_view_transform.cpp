#include "plot_view_transform.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rounding_mode.hpp"

#include <cmath>
#include <limits>
#include <optional>

namespace mmcal::plot {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::RoundingMode;

[[nodiscard]] std::optional<BigFloat> evaluateConstant(
    const expression::Expr& expression,
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits) {
    const auto compiled = compilePlotProgram(expression, request.variable, builtins, mathematics);
    if (!compiled || compiled.program->variableDependent[compiled.program->resultRegister])
        return std::nullopt;
    try {
        BigFloatPlotExecutor executor{*compiled.program, precisionBits, angles};
        const BigFloat zero = BigFloat::fromBigInt(BigInt{0}, precisionBits, RoundingMode::NearestEven);
        const auto result = executor.evaluate(zero);
        if (!result.finite())
            return std::nullopt;
        return result.value;
    }
    catch (...) {
        return std::nullopt;
    }
}

[[nodiscard]] double unitToDouble(const BigFloat& value) {
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
    if (exponent < static_cast<BigFloat::exponent_type>(std::numeric_limits<int>::min()))
        return 0.0;
    if (exponent > static_cast<BigFloat::exponent_type>(std::numeric_limits<int>::max()))
        return value.isNegative() ? -std::numeric_limits<double>::infinity()
                                  : std::numeric_limits<double>::infinity();
    return static_cast<double>(std::ldexp(significand, static_cast<int>(exponent)));
}

[[nodiscard]] std::optional<double> normalized(
    const BigFloat& value,
    const BigFloat& minimum,
    const BigFloat& maximum,
    std::size_t precisionBits) {
    if (!(minimum < maximum))
        return std::nullopt;
    const auto width = numeric::subtract(maximum, minimum, precisionBits, RoundingMode::NearestEven);
    const auto offset = numeric::subtract(value, minimum, precisionBits, RoundingMode::NearestEven);
    const auto ratio = numeric::divide(offset, width, precisionBits, RoundingMode::NearestEven);
    const double result = unitToDouble(ratio);
    if (!std::isfinite(result))
        return std::nullopt;
    return result;
}

} // namespace

PlotViewTransform::PlotViewTransform(
    BigFloat xMinimum,
    BigFloat xMaximum,
    BigFloat yMinimum,
    BigFloat yMaximum,
    std::size_t precisionBits,
    PlotViewportMm viewport)
    : xMinimum_(std::move(xMinimum)),
      xMaximum_(std::move(xMaximum)),
      yMinimum_(std::move(yMinimum)),
      yMaximum_(std::move(yMaximum)),
      precisionBits_(precisionBits),
      viewport_(viewport) {}

std::optional<PlotViewPointMm> PlotViewTransform::map(
    const BigFloat& x,
    const BigFloat& y) const {
    const auto xMm = mapX(x);
    const auto yMm = mapY(y);
    if (!xMm || !yMm)
        return std::nullopt;
    return PlotViewPointMm{*xMm, *yMm};
}

std::optional<double> PlotViewTransform::mapX(const BigFloat& x) const {
    const auto nx = normalized(x, xMinimum_, xMaximum_, precisionBits_);
    if (!nx)
        return std::nullopt;
    return *nx * viewport_.widthMm;
}

std::optional<double> PlotViewTransform::mapY(const BigFloat& y) const {
    const auto ny = normalized(y, yMinimum_, yMaximum_, precisionBits_);
    if (!ny)
        return std::nullopt;
    return *ny * viewport_.heightMm;
}

std::optional<PlotViewTransform> makePlotViewTransform(
    const CurveViewRange2D& range,
    std::size_t precisionBits,
    PlotViewportMm viewport) {
    if (precisionBits < 8 || !(viewport.widthMm > 0.0) || !(viewport.heightMm > 0.0)
        || !(range.x.viewMinimum < range.x.viewMaximum)
        || !(range.y.viewMinimum < range.y.viewMaximum))
        return std::nullopt;
    return PlotViewTransform{
        range.x.viewMinimum, range.x.viewMaximum,
        range.y.viewMinimum, range.y.viewMaximum, precisionBits, viewport};
}

std::optional<PlotViewTransform> makePlotViewTransform(
    const std::vector<PlotRequest>& requests,
    const PlotRangeEstimate& range,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits,
    PlotViewportMm viewport) {
    if (requests.empty() || precisionBits < 8
        || !(viewport.widthMm > 0.0) || !(viewport.heightMm > 0.0))
        return std::nullopt;

    std::optional<BigFloat> xMinimum;
    std::optional<BigFloat> xMaximum;
    for (const auto& request : requests) {
        const auto lower = evaluateConstant(
            request.lower, request, builtins, mathematics, angles, precisionBits);
        const auto upper = evaluateConstant(
            request.upper, request, builtins, mathematics, angles, precisionBits);
        if (!lower || !upper || !(*lower < *upper))
            return std::nullopt;
        if (!xMinimum || *lower < *xMinimum)
            xMinimum = *lower;
        if (!xMaximum || *xMaximum < *upper)
            xMaximum = *upper;
    }

    if (!xMinimum || !xMaximum || !(*xMinimum < *xMaximum)
        || !(range.viewMinimum < range.viewMaximum))
        return std::nullopt;

    const CurveRangeEstimate1D xRange{
        *xMinimum, *xMaximum, *xMinimum, *xMaximum, *xMinimum, *xMaximum,
        0, 0, false, false};
    return makePlotViewTransform(
        CurveViewRange2D{xRange, range}, precisionBits, viewport);
}

std::optional<PlotViewTransform> makePlotViewTransform(
    const PlotRequest& request,
    const PlotRangeEstimate& range,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits,
    PlotViewportMm viewport) {
    return makePlotViewTransform(
        std::vector<PlotRequest>{request}, range, builtins, mathematics, angles,
        precisionBits, viewport);
}

} // namespace mmcal::plot
