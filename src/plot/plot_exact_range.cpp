#include "plot_exact_range.hpp"

#include "mathematics/exact_trigonometry.hpp"
#include "plot_program.hpp"
#include "simplification/simplifier.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::plot {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;
using numeric::RoundingMode;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{numeric::Number{BigInt{value}}};
}

[[nodiscard]] Expr simplify(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    simplification::SimplificationContext context{builtins, mathematics, angles, assumptions};
    return simplification::Simplifier{}.simplify(std::move(expression), context);
}

[[nodiscard]] Expr binary(
    BuiltinId id,
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplify(
        Expr::call(builtins.symbol(id), {std::move(lhs), std::move(rhs)}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr square(
    Expr value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return binary(BuiltinId::Power, std::move(value), integer(2),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr sqrtExpr(
    Expr value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplify(
        Expr::call(builtins.symbol(BuiltinId::Sqrt), {std::move(value)}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] bool exactZero(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] std::optional<BigFloat> constantValue(
    const Expr& expression,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits) {
    const auto compiled = compilePlotProgram(
        expression, request.parameter, builtins, mathematics);
    if (!compiled || compiled.program->variableDependent[compiled.program->resultRegister])
        return std::nullopt;
    try {
        BigFloatPlotExecutor executor{*compiled.program, precisionBits, angles};
        const BigFloat zero = BigFloat::fromBigInt(
            BigInt{0}, precisionBits, RoundingMode::NearestEven);
        const auto value = executor.evaluate(zero);
        if (!value.finite())
            return std::nullopt;
        return value.value;
    }
    catch (...) {
        return std::nullopt;
    }
}

struct ScalarExtent final {
    BigFloat minimum;
    BigFloat maximum;
    std::size_t count = 0;
};

void include(ScalarExtent& extent, const BigFloat& value) {
    if (extent.count == 0) {
        extent.minimum = value;
        extent.maximum = value;
    }
    else {
        if (value < extent.minimum)
            extent.minimum = value;
        if (value > extent.maximum)
            extent.maximum = value;
    }
    ++extent.count;
}

[[nodiscard]] std::optional<CurveRangeEstimate1D> finish(
    const ScalarExtent& extent,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options) {
    if (extent.count == 0)
        return std::nullopt;
    return makeCurveRangeEstimate(
        extent.minimum, extent.maximum, precisionBits, extent.count, options);
}

[[nodiscard]] Expr quadraticValue(
    const Expr& p0,
    const Expr& p1,
    const Expr& p2,
    const Expr& t,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const Expr oneMinus = binary(BuiltinId::Subtract, integer(1), t,
        builtins, mathematics, angles, assumptions);
    const Expr first = binary(BuiltinId::Multiply,
        square(oneMinus, builtins, mathematics, angles, assumptions), p0,
        builtins, mathematics, angles, assumptions);
    const Expr middle = binary(BuiltinId::Multiply,
        binary(BuiltinId::Multiply, integer(2),
            binary(BuiltinId::Multiply, oneMinus, t,
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        p1, builtins, mathematics, angles, assumptions);
    const Expr last = binary(BuiltinId::Multiply,
        square(t, builtins, mathematics, angles, assumptions), p2,
        builtins, mathematics, angles, assumptions);
    return binary(BuiltinId::Add,
        binary(BuiltinId::Add, first, middle, builtins, mathematics, angles, assumptions),
        last, builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr cubicValue(
    const Expr& p0,
    const Expr& p1,
    const Expr& p2,
    const Expr& p3,
    const Expr& t,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const Expr oneMinus = binary(BuiltinId::Subtract, integer(1), t,
        builtins, mathematics, angles, assumptions);
    const Expr omt2 = square(oneMinus, builtins, mathematics, angles, assumptions);
    const Expr t2 = square(t, builtins, mathematics, angles, assumptions);
    const Expr omt3 = binary(BuiltinId::Multiply, omt2, oneMinus,
        builtins, mathematics, angles, assumptions);
    const Expr t3 = binary(BuiltinId::Multiply, t2, t,
        builtins, mathematics, angles, assumptions);
    const Expr first = binary(BuiltinId::Multiply, omt3, p0,
        builtins, mathematics, angles, assumptions);
    const Expr second = binary(BuiltinId::Multiply,
        binary(BuiltinId::Multiply, integer(3),
            binary(BuiltinId::Multiply, omt2, t,
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        p1, builtins, mathematics, angles, assumptions);
    const Expr third = binary(BuiltinId::Multiply,
        binary(BuiltinId::Multiply, integer(3),
            binary(BuiltinId::Multiply, oneMinus, t2,
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        p2, builtins, mathematics, angles, assumptions);
    const Expr fourth = binary(BuiltinId::Multiply, t3, p3,
        builtins, mathematics, angles, assumptions);
    return binary(BuiltinId::Add,
        binary(BuiltinId::Add,
            binary(BuiltinId::Add, first, second, builtins, mathematics, angles, assumptions),
            third, builtins, mathematics, angles, assumptions),
        fourth, builtins, mathematics, angles, assumptions);
}

[[nodiscard]] bool interiorUnitParameter(
    const Expr& t,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits) {
    const auto value = constantValue(t, request, builtins, mathematics, angles, precisionBits);
    if (!value)
        return false;
    const BigFloat zero = BigFloat::fromBigInt(BigInt{0}, precisionBits, RoundingMode::NearestEven);
    const BigFloat one = BigFloat::fromBigInt(BigInt{1}, precisionBits, RoundingMode::NearestEven);
    return *value > zero && *value < one;
}

[[nodiscard]] std::optional<CurveRangeEstimate1D> quadraticRange(
    const Expr& p0,
    const Expr& p1,
    const Expr& p2,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options) {
    ScalarExtent extent;
    const auto v0 = constantValue(p0, request, builtins, mathematics, angles, precisionBits);
    const auto v2 = constantValue(p2, request, builtins, mathematics, angles, precisionBits);
    if (!v0 || !v2)
        return std::nullopt;
    include(extent, *v0);
    include(extent, *v2);

    const Expr denominator = binary(BuiltinId::Add,
        binary(BuiltinId::Subtract, p0,
            binary(BuiltinId::Multiply, integer(2), p1,
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        p2, builtins, mathematics, angles, assumptions);
    if (!exactZero(denominator)) {
        const Expr t = binary(BuiltinId::Divide,
            binary(BuiltinId::Subtract, p0, p1, builtins, mathematics, angles, assumptions),
            denominator, builtins, mathematics, angles, assumptions);
        if (interiorUnitParameter(t, request, builtins, mathematics, angles, precisionBits)) {
            const auto value = constantValue(
                quadraticValue(p0, p1, p2, t, builtins, mathematics, angles, assumptions),
                request, builtins, mathematics, angles, precisionBits);
            if (!value)
                return std::nullopt;
            include(extent, *value);
        }
    }
    return finish(extent, precisionBits, options);
}

[[nodiscard]] std::optional<CurveRangeEstimate1D> cubicRange(
    const Expr& p0,
    const Expr& p1,
    const Expr& p2,
    const Expr& p3,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options) {
    ScalarExtent extent;
    const auto v0 = constantValue(p0, request, builtins, mathematics, angles, precisionBits);
    const auto v3 = constantValue(p3, request, builtins, mathematics, angles, precisionBits);
    if (!v0 || !v3)
        return std::nullopt;
    include(extent, *v0);
    include(extent, *v3);

    const Expr a = binary(BuiltinId::Add,
        binary(BuiltinId::Subtract,
            binary(BuiltinId::Add,
                binary(BuiltinId::Multiply, integer(-1), p0, builtins, mathematics, angles, assumptions),
                binary(BuiltinId::Multiply, integer(3), p1, builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions),
            binary(BuiltinId::Multiply, integer(3), p2, builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        p3, builtins, mathematics, angles, assumptions);
    const Expr b = binary(BuiltinId::Add,
        binary(BuiltinId::Subtract,
            binary(BuiltinId::Multiply, integer(3), p0, builtins, mathematics, angles, assumptions),
            binary(BuiltinId::Multiply, integer(6), p1, builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        binary(BuiltinId::Multiply, integer(3), p2, builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
    const Expr c = binary(BuiltinId::Add,
        binary(BuiltinId::Multiply, integer(-3), p0, builtins, mathematics, angles, assumptions),
        binary(BuiltinId::Multiply, integer(3), p1, builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);

    std::vector<Expr> roots;
    if (exactZero(a)) {
        if (!exactZero(b)) {
            roots.push_back(binary(BuiltinId::Divide,
                binary(BuiltinId::Multiply, integer(-1), c,
                    builtins, mathematics, angles, assumptions),
                binary(BuiltinId::Multiply, integer(2), b,
                    builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions));
        }
    }
    else {
        const Expr discriminant = binary(BuiltinId::Subtract,
            square(b, builtins, mathematics, angles, assumptions),
            binary(BuiltinId::Multiply, integer(3),
                binary(BuiltinId::Multiply, a, c,
                    builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions);
        const Expr radical = sqrtExpr(discriminant, builtins, mathematics, angles, assumptions);
        const Expr denominator = binary(BuiltinId::Multiply, integer(3), a,
            builtins, mathematics, angles, assumptions);
        for (BuiltinId operation : {BuiltinId::Add, BuiltinId::Subtract}) {
            roots.push_back(binary(BuiltinId::Divide,
                binary(operation,
                    binary(BuiltinId::Multiply, integer(-1), b,
                        builtins, mathematics, angles, assumptions),
                    radical, builtins, mathematics, angles, assumptions),
                denominator, builtins, mathematics, angles, assumptions));
        }
    }

    for (const Expr& t : roots) {
        if (!interiorUnitParameter(t, request, builtins, mathematics, angles, precisionBits))
            continue;
        const auto value = constantValue(
            cubicValue(p0, p1, p2, p3, t, builtins, mathematics, angles, assumptions),
            request, builtins, mathematics, angles, precisionBits);
        if (!value)
            return std::nullopt;
        include(extent, *value);
    }
    return finish(extent, precisionBits, options);
}

[[nodiscard]] std::optional<Expr> angleExpressionFromTurns(
    const Rational& turns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Degree:
        return Expr{numeric::Number{turns * Rational{BigInt{360}}}};
    case mathematics::AngleUnit::Gradian:
        return Expr{numeric::Number{turns * Rational{BigInt{400}}}};
    case mathematics::AngleUnit::Radian: {
        const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
        if (!pi)
            return std::nullopt;
        const Rational coefficient = turns * Rational{BigInt{2}};
        if (coefficient.isZero())
            return integer(0);
        if (coefficient == Rational{BigInt{1}})
            return Expr{pi->symbol};
        return Expr::call(
            builtins.symbol(BuiltinId::Multiply),
            {Expr{numeric::Number{coefficient}}, Expr{pi->symbol}});
    }
    }
    return std::nullopt;
}

[[nodiscard]] Expr ellipseCoordinateAt(
    const Expr& center,
    const Expr& cosineAxis,
    const Expr& sineAxis,
    const Expr& angle,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const Expr cosine = simplify(
        Expr::call(builtins.symbol(BuiltinId::Cos), {angle}),
        builtins, mathematics, angles, assumptions);
    const Expr sine = simplify(
        Expr::call(builtins.symbol(BuiltinId::Sin), {angle}),
        builtins, mathematics, angles, assumptions);
    return binary(BuiltinId::Add,
        binary(BuiltinId::Add, center,
            binary(BuiltinId::Multiply, cosineAxis, cosine,
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        binary(BuiltinId::Multiply, sineAxis, sine,
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Rational absolute(Rational value) {
    return value.numerator().isNegative() ? -value : value;
}

[[nodiscard]] Rational positiveModuloOne(const Rational& value) {
    BigInt remainder = value.numerator() % value.denominator();
    if (remainder.isNegative())
        remainder += value.denominator();
    return Rational{std::move(remainder), value.denominator()};
}

[[nodiscard]] bool turnOccursOnArc(
    const Rational& candidate,
    const Rational& start,
    const Rational& sweep) {
    const Rational span = absolute(sweep);
    if (span >= Rational{BigInt{1}})
        return true;
    const Rational delta = sweep.numerator().isNegative()
        ? positiveModuloOne(start - candidate)
        : positiveModuloOne(candidate - start);
    return delta <= span;
}

struct EllipseCoordinateRange final {
    CurveRangeEstimate1D range;
    bool tight = true;
};

[[nodiscard]] std::optional<EllipseCoordinateRange> ellipseCoordinateRange(
    const Expr& center,
    const Expr& cosineAxis,
    const Expr& sineAxis,
    std::optional<std::pair<Rational, Rational>> arc,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options) {
    const Expr radius = sqrtExpr(
        binary(BuiltinId::Add,
            square(cosineAxis, builtins, mathematics, angles, assumptions),
            square(sineAxis, builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
    const Expr fullMinimum = binary(BuiltinId::Subtract, center, radius,
        builtins, mathematics, angles, assumptions);
    const Expr fullMaximum = binary(BuiltinId::Add, center, radius,
        builtins, mathematics, angles, assumptions);
    const auto fullMinValue = constantValue(
        fullMinimum, request, builtins, mathematics, angles, precisionBits);
    const auto fullMaxValue = constantValue(
        fullMaximum, request, builtins, mathematics, angles, precisionBits);
    if (!fullMinValue || !fullMaxValue)
        return std::nullopt;
    if (!arc) {
        const auto range = makeCurveRangeEstimate(
            *fullMinValue, *fullMaxValue, precisionBits, 2, options);
        return range ? std::optional<EllipseCoordinateRange>{{*range, true}} : std::nullopt;
    }

    ScalarExtent extent;
    const Rational start = arc->first;
    const Rational sweep = arc->second;
    const auto startAngle = angleExpressionFromTurns(start, builtins, mathematics, angles);
    const auto endAngle = angleExpressionFromTurns(start + sweep, builtins, mathematics, angles);
    if (!startAngle || !endAngle)
        return std::nullopt;
    const auto startValue = constantValue(
        ellipseCoordinateAt(center, cosineAxis, sineAxis, *startAngle,
            builtins, mathematics, angles, assumptions),
        request, builtins, mathematics, angles, precisionBits);
    const auto endValue = constantValue(
        ellipseCoordinateAt(center, cosineAxis, sineAxis, *endAngle,
            builtins, mathematics, angles, assumptions),
        request, builtins, mathematics, angles, precisionBits);
    if (!startValue || !endValue)
        return std::nullopt;
    include(extent, *startValue);
    include(extent, *endValue);

    if (exactZero(cosineAxis) && exactZero(sineAxis)) {
        const auto range = finish(extent, precisionBits, options);
        return range ? std::optional<EllipseCoordinateRange>{{*range, true}} : std::nullopt;
    }

    std::optional<Rational> maximumTurn;
    try {
        if (const auto angle = mathematics::simplifyExactAtan2(
                sineAxis, cosineAxis, builtins, mathematics, angles)) {
            if (const auto exact = mathematics::extractExactAngle(
                    *angle, builtins, mathematics, angles))
                maximumTurn = exact->turns;
        }
    }
    catch (...) {
        maximumTurn.reset();
    }

    if (!maximumTurn) {
        // 一般のatan2を近似してarc内外を決めない。完全楕円bboxなら常に安全。
        const auto range = makeCurveRangeEstimate(
            *fullMinValue, *fullMaxValue, precisionBits, 2, options);
        return range ? std::optional<EllipseCoordinateRange>{{*range, false}} : std::nullopt;
    }

    if (turnOccursOnArc(*maximumTurn, start, sweep))
        include(extent, *fullMaxValue);
    const Rational minimumTurn = *maximumTurn + Rational{BigInt{1}, BigInt{2}};
    if (turnOccursOnArc(minimumTurn, start, sweep))
        include(extent, *fullMinValue);
    const auto range = finish(extent, precisionBits, options);
    return range ? std::optional<EllipseCoordinateRange>{{*range, true}} : std::nullopt;
}

[[nodiscard]] std::optional<ExactCurveRange2D> pointOrLineRange(
    const std::vector<SymbolicPoint2D>& points,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options) {
    ScalarExtent x;
    ScalarExtent y;
    for (const auto& point : points) {
        const auto px = constantValue(point.x, request, builtins, mathematics, angles, precisionBits);
        const auto py = constantValue(point.y, request, builtins, mathematics, angles, precisionBits);
        if (!px || !py)
            return std::nullopt;
        include(x, *px);
        include(y, *py);
    }
    const auto xr = finish(x, precisionBits, options);
    const auto yr = finish(y, precisionBits, options);
    if (!xr || !yr)
        return std::nullopt;
    return ExactCurveRange2D{*xr, *yr, true};
}

} // namespace

std::optional<ExactCurveRange2D> estimateExactCurveRange(
    const ExactCurveGeometry& geometry,
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::size_t precisionBits,
    const PlotRangeEstimatorOptions& options) {
    return std::visit([&](const auto& value) -> std::optional<ExactCurveRange2D> {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, ExactPointGeometry>)
            return pointOrLineRange(
                {value.point}, request, builtins, mathematics, angles, precisionBits, options);
        else if constexpr (std::is_same_v<T, ExactLineGeometry>)
            return pointOrLineRange(
                {value.start, value.end}, request, builtins, mathematics, angles, precisionBits, options);
        else if constexpr (std::is_same_v<T, ExactQuadraticBezierGeometry>) {
            const auto x = quadraticRange(
                value.start.x, value.control.x, value.end.x,
                request, builtins, mathematics, angles, assumptions, precisionBits, options);
            const auto y = quadraticRange(
                value.start.y, value.control.y, value.end.y,
                request, builtins, mathematics, angles, assumptions, precisionBits, options);
            if (!x || !y)
                return std::nullopt;
            return ExactCurveRange2D{*x, *y, true};
        }
        else if constexpr (std::is_same_v<T, ExactCubicBezierGeometry>) {
            const auto x = cubicRange(
                value.start.x, value.control1.x, value.control2.x, value.end.x,
                request, builtins, mathematics, angles, assumptions, precisionBits, options);
            const auto y = cubicRange(
                value.start.y, value.control1.y, value.control2.y, value.end.y,
                request, builtins, mathematics, angles, assumptions, precisionBits, options);
            if (!x || !y)
                return std::nullopt;
            return ExactCurveRange2D{*x, *y, true};
        }
        else if constexpr (std::is_same_v<T, ExactCircleGeometry>
            || std::is_same_v<T, ExactEllipseGeometry>) {
            const auto x = ellipseCoordinateRange(
                value.center.x, value.cosineAxis.x, value.sineAxis.x, std::nullopt,
                request, builtins, mathematics, angles, assumptions, precisionBits, options);
            const auto y = ellipseCoordinateRange(
                value.center.y, value.cosineAxis.y, value.sineAxis.y, std::nullopt,
                request, builtins, mathematics, angles, assumptions, precisionBits, options);
            if (!x || !y)
                return std::nullopt;
            return ExactCurveRange2D{x->range, y->range, true};
        }
        else {
            const auto x = ellipseCoordinateRange(
                value.center.x, value.cosineAxis.x, value.sineAxis.x,
                std::pair<Rational, Rational>{value.startTurns, value.sweepTurns},
                request, builtins, mathematics, angles, assumptions, precisionBits, options);
            const auto y = ellipseCoordinateRange(
                value.center.y, value.cosineAxis.y, value.sineAxis.y,
                std::pair<Rational, Rational>{value.startTurns, value.sweepTurns},
                request, builtins, mathematics, angles, assumptions, precisionBits, options);
            if (!x || !y)
                return std::nullopt;
            return ExactCurveRange2D{x->range, y->range, x->tight && y->tight};
        }
    }, geometry.value);
}

} // namespace mmcal::plot
