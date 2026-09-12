#pragma once

#include "plot_request.hpp"

#include "evaluation/builtin_registry.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/rational.hpp"

#include <optional>
#include <variant>

namespace mmcal::plot {

// sampling/layout精度から独立した，式そのものに由来するexact 2D geometry。
// 座標はExprのまま保持し，BigFloat化は後段のview/lowering境界まで遅延する。
struct SymbolicPoint2D final {
    expression::Expr x;
    expression::Expr y;

    [[nodiscard]] bool operator==(const SymbolicPoint2D&) const = default;
};

struct ExactPointGeometry final {
    SymbolicPoint2D point;
    [[nodiscard]] bool operator==(const ExactPointGeometry&) const = default;
};

struct ExactLineGeometry final {
    SymbolicPoint2D start;
    SymbolicPoint2D end;
    [[nodiscard]] bool operator==(const ExactLineGeometry&) const = default;
};

struct ExactQuadraticBezierGeometry final {
    SymbolicPoint2D start;
    SymbolicPoint2D control;
    SymbolicPoint2D end;
    [[nodiscard]] bool operator==(const ExactQuadraticBezierGeometry&) const = default;
};

struct ExactCubicBezierGeometry final {
    SymbolicPoint2D start;
    SymbolicPoint2D control1;
    SymbolicPoint2D control2;
    SymbolicPoint2D end;
    [[nodiscard]] bool operator==(const ExactCubicBezierGeometry&) const = default;
};

struct ExactCircleGeometry final {
    SymbolicPoint2D center;
    expression::Expr radius;
    // 元parameterizationのphase基底。後段で向き/開始phaseを必要とするときに再利用する。
    SymbolicPoint2D cosineAxis;
    SymbolicPoint2D sineAxis;
    numeric::Rational periodTurns;
    [[nodiscard]] bool operator==(const ExactCircleGeometry&) const = default;
};

struct ExactEllipseGeometry final {
    SymbolicPoint2D center;
    SymbolicPoint2D cosineAxis;
    SymbolicPoint2D sineAxis;
    numeric::Rational periodTurns;
    [[nodiscard]] bool operator==(const ExactEllipseGeometry&) const = default;
};

// 楕円の部分弧または複数周回を，元trig phaseのturnでexactに保持する。
// 1 turn = 360 Deg = 2 Pi Rad = 400 Grad。負sweepは逆向きを表す。
struct ExactEllipticArcGeometry final {
    SymbolicPoint2D center;
    SymbolicPoint2D cosineAxis;
    SymbolicPoint2D sineAxis;
    numeric::Rational startTurns;
    numeric::Rational sweepTurns;
    numeric::Rational periodTurns;
    [[nodiscard]] bool operator==(const ExactEllipticArcGeometry&) const = default;
};

using ExactCurveGeometryValue = std::variant<
    ExactPointGeometry,
    ExactLineGeometry,
    ExactQuadraticBezierGeometry,
    ExactCubicBezierGeometry,
    ExactCircleGeometry,
    ExactEllipseGeometry,
    ExactEllipticArcGeometry>;

struct ExactCurveGeometry final {
    ExactCurveGeometryValue value;
};

// parameterについて3次以下の多項式，または
// C + A cos(qt+b) + B sin(qt+b) をsymbolic geometryへ落とす。
// trig affine像は完全1周期ならCircle/Ellipse，それ以外はEllipticArcとして保持する。
[[nodiscard]] std::optional<ExactCurveGeometry> recognizeExactParametricGeometry(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// 旧pipelineとの互換用。新規処理はExactCurveGeometryを利用する。
struct ParametricEllipseRecognition final {
    numeric::Rational periodTurns;
};

[[nodiscard]] std::optional<ParametricEllipseRecognition> recognizeParametricEllipse(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::plot
