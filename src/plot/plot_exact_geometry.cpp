#include "plot_exact_geometry.hpp"

#include "plot_periodicity.hpp"
#include "simplification/simplifier.hpp"
#include "mathematics/exact_trigonometry.hpp"
#include "symbolic/polynomial.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <variant>

namespace mmcal::plot {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Rational;

struct HarmonicLinearForm final {
    Expr constant{numeric::Number{BigInt{0}}};
    Expr cosine{numeric::Number{BigInt{0}}};
    Expr sine{numeric::Number{BigInt{0}}};
    std::optional<Expr> argument;
};

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

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{numeric::Number{BigInt{value}}};
}

[[nodiscard]] bool exactZero(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] bool exactNonZeroRealNumber(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isReal()
        && !expression.asNumber().isZero();
}

[[nodiscard]] Expr polynomialValue(
    const symbolic::ExpressionPolynomial& polynomial,
    const Expr& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr value = polynomial.coefficient(polynomial.degree());
    for (std::size_t exponent = polynomial.degree(); exponent-- > 0;) {
        value = binary(BuiltinId::Add,
            binary(BuiltinId::Multiply, std::move(value), parameter,
                builtins, mathematics, angles, assumptions),
            polynomial.coefficient(exponent),
            builtins, mathematics, angles, assumptions);
    }
    return value;
}

[[nodiscard]] Expr polynomialDerivativeValue(
    const symbolic::ExpressionPolynomial& polynomial,
    const Expr& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (polynomial.degree() == 0)
        return integer(0);
    Expr value = binary(BuiltinId::Multiply,
        polynomial.coefficient(polynomial.degree()),
        integer(static_cast<std::int64_t>(polynomial.degree())),
        builtins, mathematics, angles, assumptions);
    for (std::size_t exponent = polynomial.degree() - 1; exponent > 0; --exponent) {
        value = binary(BuiltinId::Add,
            binary(BuiltinId::Multiply, std::move(value), parameter,
                builtins, mathematics, angles, assumptions),
            binary(BuiltinId::Multiply,
                polynomial.coefficient(exponent),
                integer(static_cast<std::int64_t>(exponent)),
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions);
    }
    return value;
}

[[nodiscard]] SymbolicPoint2D polynomialPoint(
    const symbolic::ExpressionPolynomial& x,
    const symbolic::ExpressionPolynomial& y,
    const Expr& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return SymbolicPoint2D{
        polynomialValue(x, parameter, builtins, mathematics, angles, assumptions),
        polynomialValue(y, parameter, builtins, mathematics, angles, assumptions)};
}

[[nodiscard]] SymbolicPoint2D polynomialDerivativePoint(
    const symbolic::ExpressionPolynomial& x,
    const symbolic::ExpressionPolynomial& y,
    const Expr& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return SymbolicPoint2D{
        polynomialDerivativeValue(x, parameter, builtins, mathematics, angles, assumptions),
        polynomialDerivativeValue(y, parameter, builtins, mathematics, angles, assumptions)};
}

[[nodiscard]] SymbolicPoint2D addScaled(
    const SymbolicPoint2D& point,
    const SymbolicPoint2D& derivative,
    const Expr& scale,
    BuiltinId operation,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return SymbolicPoint2D{
        binary(operation, point.x,
            binary(BuiltinId::Multiply, scale, derivative.x,
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        binary(operation, point.y,
            binary(BuiltinId::Multiply, scale, derivative.y,
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions)};
}

[[nodiscard]] std::optional<ExactCurveGeometry> recognizePolynomialGeometry(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto x = symbolic::toExpressionPolynomial(
        request.xExpression, request.parameter, builtins, mathematics, angles);
    const auto y = symbolic::toExpressionPolynomial(
        request.yExpression, request.parameter, builtins, mathematics, angles);
    if (!x || !y)
        return std::nullopt;
    const std::size_t degree = std::max(x->degree(), y->degree());
    if (degree > 3)
        return std::nullopt;

    const SymbolicPoint2D start = polynomialPoint(
        *x, *y, request.lower, builtins, mathematics, angles, assumptions);
    if (degree == 0)
        return ExactCurveGeometry{ExactPointGeometry{start}};

    const SymbolicPoint2D end = polynomialPoint(
        *x, *y, request.upper, builtins, mathematics, angles, assumptions);
    if (degree == 1)
        return ExactCurveGeometry{ExactLineGeometry{start, end}};

    const Expr span = binary(BuiltinId::Subtract, request.upper, request.lower,
        builtins, mathematics, angles, assumptions);
    const SymbolicPoint2D startDerivative = polynomialDerivativePoint(
        *x, *y, request.lower, builtins, mathematics, angles, assumptions);
    if (degree == 2) {
        const Expr halfSpan = binary(BuiltinId::Divide, span, integer(2),
            builtins, mathematics, angles, assumptions);
        return ExactCurveGeometry{ExactQuadraticBezierGeometry{
            start,
            addScaled(start, startDerivative, halfSpan, BuiltinId::Add,
                builtins, mathematics, angles, assumptions),
            end}};
    }

    const Expr thirdSpan = binary(BuiltinId::Divide, span, integer(3),
        builtins, mathematics, angles, assumptions);
    const SymbolicPoint2D endDerivative = polynomialDerivativePoint(
        *x, *y, request.upper, builtins, mathematics, angles, assumptions);
    return ExactCurveGeometry{ExactCubicBezierGeometry{
        start,
        addScaled(start, startDerivative, thirdSpan, BuiltinId::Add,
            builtins, mathematics, angles, assumptions),
        addScaled(end, endDerivative, thirdSpan, BuiltinId::Subtract,
            builtins, mathematics, angles, assumptions),
        end}};
}

[[nodiscard]] std::optional<HarmonicLinearForm> merge(
    HarmonicLinearForm lhs,
    HarmonicLinearForm rhs,
    BuiltinId operation,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (lhs.argument && rhs.argument && *lhs.argument != *rhs.argument)
        return std::nullopt;
    if (!lhs.argument)
        lhs.argument = rhs.argument;
    lhs.constant = binary(operation, lhs.constant, rhs.constant,
        builtins, mathematics, angles, assumptions);
    lhs.cosine = binary(operation, lhs.cosine, rhs.cosine,
        builtins, mathematics, angles, assumptions);
    lhs.sine = binary(operation, lhs.sine, rhs.sine,
        builtins, mathematics, angles, assumptions);
    return lhs;
}

[[nodiscard]] HarmonicLinearForm scale(
    HarmonicLinearForm form,
    const Expr& factor,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    form.constant = binary(BuiltinId::Multiply, form.constant, factor,
        builtins, mathematics, angles, assumptions);
    form.cosine = binary(BuiltinId::Multiply, form.cosine, factor,
        builtins, mathematics, angles, assumptions);
    form.sine = binary(BuiltinId::Multiply, form.sine, factor,
        builtins, mathematics, angles, assumptions);
    return form;
}

[[nodiscard]] std::optional<HarmonicLinearForm> decompose(
    const Expr& expression,
    const expression::Symbol& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!symbolic::containsSymbol(expression, parameter))
        return HarmonicLinearForm{expression, integer(0), integer(0), std::nullopt};
    if (!expression.isCall())
        return std::nullopt;

    const auto& call = expression.asCall();
    const auto* builtin = builtins.find(call.head);
    if (!builtin)
        return std::nullopt;

    if ((builtin->id == BuiltinId::Sin || builtin->id == BuiltinId::Cos)
        && call.arguments.size() == 1) {
        HarmonicLinearForm form;
        form.argument = simplify(call.arguments[0], builtins, mathematics, angles, assumptions);
        if (builtin->id == BuiltinId::Cos)
            form.cosine = integer(1);
        else
            form.sine = integer(1);
        return form;
    }

    if (builtin->id == BuiltinId::Negate && call.arguments.size() == 1) {
        auto child = decompose(call.arguments[0], parameter, builtins, mathematics, angles, assumptions);
        if (!child)
            return std::nullopt;
        return scale(std::move(*child), integer(-1), builtins, mathematics, angles, assumptions);
    }

    if ((builtin->id == BuiltinId::Add || builtin->id == BuiltinId::Subtract)
        && !call.arguments.empty()) {
        auto result = decompose(call.arguments.front(), parameter, builtins, mathematics, angles, assumptions);
        if (!result)
            return std::nullopt;
        for (std::size_t i = 1; i < call.arguments.size(); ++i) {
            auto child = decompose(call.arguments[i], parameter, builtins, mathematics, angles, assumptions);
            if (!child)
                return std::nullopt;
            auto merged = merge(std::move(*result), std::move(*child),
                builtin->id == BuiltinId::Add ? BuiltinId::Add : BuiltinId::Subtract,
                builtins, mathematics, angles, assumptions);
            if (!merged)
                return std::nullopt;
            result = std::move(merged);
        }
        return result;
    }

    if (builtin->id == BuiltinId::Multiply && !call.arguments.empty()) {
        std::optional<std::size_t> dependentIndex;
        Expr factor = integer(1);
        for (std::size_t i = 0; i < call.arguments.size(); ++i) {
            if (symbolic::containsSymbol(call.arguments[i], parameter)) {
                if (dependentIndex)
                    return std::nullopt;
                dependentIndex = i;
            }
            else {
                factor = binary(BuiltinId::Multiply, factor, call.arguments[i],
                    builtins, mathematics, angles, assumptions);
            }
        }
        if (!dependentIndex)
            return HarmonicLinearForm{expression, integer(0), integer(0), std::nullopt};
        auto child = decompose(call.arguments[*dependentIndex], parameter,
            builtins, mathematics, angles, assumptions);
        if (!child)
            return std::nullopt;
        return scale(std::move(*child), factor, builtins, mathematics, angles, assumptions);
    }

    if (builtin->id == BuiltinId::Divide && call.arguments.size() == 2
        && !symbolic::containsSymbol(call.arguments[1], parameter)) {
        auto child = decompose(call.arguments[0], parameter,
            builtins, mathematics, angles, assumptions);
        if (!child)
            return std::nullopt;
        const Expr reciprocal = binary(BuiltinId::Divide, integer(1), call.arguments[1],
            builtins, mathematics, angles, assumptions);
        return scale(std::move(*child), reciprocal,
            builtins, mathematics, angles, assumptions);
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<Rational> affineRationalCoefficient(
    const Expr& argument,
    const expression::Symbol& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto polynomial = symbolic::toExpressionPolynomial(
        argument, parameter, builtins, mathematics, angles);
    if (!polynomial || polynomial->degree() != 1)
        return std::nullopt;
    const Expr& coefficient = polynomial->coefficient(1);
    if (!coefficient.isNumber() || !coefficient.asNumber().isReal())
        return std::nullopt;
    Rational value = coefficient.asNumber().asReal().toRational();
    if (value.isZero())
        return std::nullopt;
    if (value < Rational{BigInt{0}})
        value = -value;
    return value;
}

[[nodiscard]] std::optional<ExactCurveGeometry> recognizeEllipseGeometry(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    auto x = decompose(request.xExpression, request.parameter,
        builtins, mathematics, angles, assumptions);
    auto y = decompose(request.yExpression, request.parameter,
        builtins, mathematics, angles, assumptions);
    if (!x || !y || !x->argument || !y->argument || *x->argument != *y->argument)
        return std::nullopt;

    const auto coefficient = affineRationalCoefficient(
        *x->argument, request.parameter, builtins, mathematics, angles);
    if (!coefficient)
        return std::nullopt;
    const Rational periodTurns = Rational{BigInt{1}} / *coefficient;
    const auto repetitions = exactParametricPeriodRepetitions(
        request, periodTurns, builtins, mathematics, angles, assumptions);

    const SymbolicPoint2D center{x->constant, y->constant};
    const SymbolicPoint2D cosineAxis{x->cosine, y->cosine};
    const SymbolicPoint2D sineAxis{x->sine, y->sine};
    const Expr determinant = binary(BuiltinId::Subtract,
        binary(BuiltinId::Multiply, cosineAxis.x, sineAxis.y,
            builtins, mathematics, angles, assumptions),
        binary(BuiltinId::Multiply, cosineAxis.y, sineAxis.x,
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
    if (!exactNonZeroRealNumber(determinant))
        return std::nullopt;

    const Expr dot = binary(BuiltinId::Add,
        binary(BuiltinId::Multiply, cosineAxis.x, sineAxis.x,
            builtins, mathematics, angles, assumptions),
        binary(BuiltinId::Multiply, cosineAxis.y, sineAxis.y,
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
    const Expr cosineNorm2 = binary(BuiltinId::Add,
        binary(BuiltinId::Power, cosineAxis.x, integer(2),
            builtins, mathematics, angles, assumptions),
        binary(BuiltinId::Power, cosineAxis.y, integer(2),
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
    const Expr sineNorm2 = binary(BuiltinId::Add,
        binary(BuiltinId::Power, sineAxis.x, integer(2),
            builtins, mathematics, angles, assumptions),
        binary(BuiltinId::Power, sineAxis.y, integer(2),
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
    const Expr normDifference = binary(BuiltinId::Subtract, cosineNorm2, sineNorm2,
        builtins, mathematics, angles, assumptions);

    const bool circle = exactZero(dot) && exactZero(normDifference)
        && exactNonZeroRealNumber(cosineNorm2);
    if (repetitions == std::optional<std::size_t>{1}) {
        if (circle) {
            const Expr radius = simplify(
                Expr::call(builtins.symbol(BuiltinId::Sqrt), {cosineNorm2}),
                builtins, mathematics, angles, assumptions);
            return ExactCurveGeometry{ExactCircleGeometry{
                center, radius, cosineAxis, sineAxis, periodTurns}};
        }
        return ExactCurveGeometry{ExactEllipseGeometry{
            center, cosineAxis, sineAxis, periodTurns}};
    }

    // 部分弧は開始phaseとsweepまでexactに確定できる場合だけprimitive化する。
    // endpointが近似値等でturnへexact変換できなければ，samplingへ安全に戻す。
    const auto argumentPolynomial = symbolic::toExpressionPolynomial(
        *x->argument, request.parameter, builtins, mathematics, angles);
    if (!argumentPolynomial || argumentPolynomial->degree() != 1)
        return std::nullopt;
    const Expr startArgument = polynomialValue(
        *argumentPolynomial, request.lower, builtins, mathematics, angles, assumptions);
    const Expr endArgument = polynomialValue(
        *argumentPolynomial, request.upper, builtins, mathematics, angles, assumptions);
    const auto startAngle = mathematics::extractExactAngle(
        startArgument, builtins, mathematics, angles);
    const auto endAngle = mathematics::extractExactAngle(
        endArgument, builtins, mathematics, angles);
    if (!startAngle || !endAngle || startAngle->turns == endAngle->turns)
        return std::nullopt;

    return ExactCurveGeometry{ExactEllipticArcGeometry{
        center, cosineAxis, sineAxis,
        startAngle->turns, endAngle->turns - startAngle->turns, periodTurns}};
}

} // namespace

std::optional<ExactCurveGeometry> recognizeExactParametricGeometry(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (auto polynomial = recognizePolynomialGeometry(
            request, builtins, mathematics, angles, assumptions))
        return polynomial;
    return recognizeEllipseGeometry(request, builtins, mathematics, angles, assumptions);
}

std::optional<ParametricEllipseRecognition> recognizeParametricEllipse(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto geometry = recognizeEllipseGeometry(
        request, builtins, mathematics, angles, assumptions);
    if (!geometry)
        return std::nullopt;
    if (const auto* circle = std::get_if<ExactCircleGeometry>(&geometry->value))
        return ParametricEllipseRecognition{circle->periodTurns};
    if (const auto* ellipse = std::get_if<ExactEllipseGeometry>(&geometry->value))
        return ParametricEllipseRecognition{ellipse->periodTurns};
    return std::nullopt;
}

} // namespace mmcal::plot
