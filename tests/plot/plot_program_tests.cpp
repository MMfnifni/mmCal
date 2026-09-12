// PlotProgram / BigFloat executorの基礎契約
#include "plot_program_tests.hpp"

#include "kernel/kernel_session.hpp"
#include "mathematics/angle.hpp"
#include "numeric/big_float.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "numeric/rounding_mode.hpp"
#include "plot/plot_program.hpp"
#include "test_framework.hpp"

#include <cstddef>
#include <string_view>

namespace mmcal::tests {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;
using numeric::RoundingMode;

[[nodiscard]] BigFloat exact(std::int64_t value, std::size_t bits = 96) {
    return BigFloat::fromBigInt(BigInt{value}, bits, RoundingMode::NearestEven);
}

[[nodiscard]] bool near(
    const BigFloat& value,
    const Rational& expected,
    const Rational& tolerance) {
    const Rational actual = value.toRational();
    const Rational difference = actual >= expected ? actual - expected : expected - actual;
    return difference <= tolerance;
}

} // namespace

void runPlotProgramTests(TestRunner& tests) {
    kernel::KernelSession session;
    const auto xExpression = session.evaluate("x");
    tests.expect(xExpression.isSymbol(),
        "PlotProgram: symbolic variable prerequisite is registered");
    if (!xExpression.isSymbol())
        return;
    const auto x = xExpression.asSymbol();

    const auto compile = [&](std::string_view source) {
        return plot::compilePlotProgram(
            session.evaluate(source), x,
            session.builtinRegistry(), session.mathRegistry());
    };

    const auto line = compile("2x+3");
    tests.expect(line && line.program->geometryKind == plot::PlotCurveGeometryKind::Affine,
        "PlotProgram: recognizes affine functions for direct line lowering");
    const auto constantLine = compile("7");
    tests.expect(constantLine
            && constantLine.program->geometryKind == plot::PlotCurveGeometryKind::Constant,
        "PlotProgram: recognizes constant functions as straight geometry");

    const auto floorStep = compile("floor[x]");
    tests.expect(floorStep
            && floorStep.program->geometryKind == plot::PlotCurveGeometryKind::PiecewiseConstant,
        "PlotProgram: recognizes direct floor/ceil/sign outputs as piecewise-constant geometry");

    const auto polynomial = compile("x^2+2x+1");
    tests.expect(static_cast<bool>(polynomial)
            && polynomial.program->geometryKind == plot::PlotCurveGeometryKind::QuadraticPolynomial,
        "PlotProgram: compiles polynomial arithmetic and recognizes exact quadratic geometry");
    if (polynomial) {
        plot::BigFloatPlotExecutor executor{*polynomial.program, 96};
        const auto result = executor.evaluate(exact(2));
        tests.expect(result.status == plot::PlotNumericStatus::Finite
                && result.value.toRational() == Rational{BigInt{9}},
            "PlotProgram: evaluates polynomial without returning to Expr");
    }

    const auto cubicPolynomial = compile("3x^3-2x^2+x-5");
    tests.expect(cubicPolynomial
            && cubicPolynomial.program->geometryKind == plot::PlotCurveGeometryKind::CubicPolynomial,
        "PlotProgram: recognizes cubic polynomials for exact cubic-Bezier lowering");

    const auto quadraticPolynomial = compile("3x^2-2x+5");
    tests.expect(quadraticPolynomial
            && quadraticPolynomial.program->geometryKind
                == plot::PlotCurveGeometryKind::QuadraticPolynomial,
        "PlotProgram: recognizes quadratic polynomials for exact quadratic-Bezier lowering");

    const auto elementary = compile("sin[x]*exp[-x^2]+1/(x^2+1)");
    tests.expect(static_cast<bool>(elementary),
        "PlotProgram: compiles initial elementary-function set");
    if (elementary) {
        plot::BigFloatPlotExecutor executor{*elementary.program, 96};
        const auto zero = executor.evaluate(exact(0));
        tests.expect(zero.status == plot::PlotNumericStatus::Finite
                && zero.value.toRational() == Rational{BigInt{1}},
            "PlotProgram: elementary pipeline evaluates an exact easy point");
        const auto one = executor.evaluate(exact(1));
        tests.expect(one.status == plot::PlotNumericStatus::Finite,
            "PlotProgram: elementary pipeline evaluates a generic finite point");
    }

    const auto invariant = compile("sqrt[2]+x");
    tests.expect(static_cast<bool>(invariant),
        "PlotProgram: compiles invariant transcendental subexpressions");
    if (invariant) {
        std::size_t invariantRegisters = 0;
        for (bool dependent : invariant.program->variableDependent)
            invariantRegisters += dependent ? 0U : 1U;
        tests.expect(invariantRegisters >= 2,
            "PlotProgram: marks constant/sqrt registers for one-time evaluation");
        plot::BigFloatPlotExecutor executor{*invariant.program, 96};
        tests.expect(executor.evaluate(exact(0)).status == plot::PlotNumericStatus::Finite
                && executor.evaluate(exact(3)).status == plot::PlotNumericStatus::Finite,
            "PlotProgram: reuses invariant registers across samples");
    }

    const auto squareRoot = compile("sqrt[x]");
    if (squareRoot) {
        plot::BigFloatPlotExecutor executor{*squareRoot.program, 96};
        tests.expect(executor.evaluate(exact(-1)).status == plot::PlotNumericStatus::NonReal,
            "PlotProgram: reports real-plot sqrt domain exits explicitly");
    }

    const auto reciprocal = compile("1/x");
    if (reciprocal) {
        plot::BigFloatPlotExecutor executor{*reciprocal.program, 96};
        tests.expect(executor.evaluate(exact(0)).status == plot::PlotNumericStatus::DivisionByZero,
            "PlotProgram: keeps division-by-zero distinct from generic undefined values");
    }

    const auto logarithm = compile("log[x]");
    if (logarithm) {
        plot::BigFloatPlotExecutor executor{*logarithm.program, 96};
        tests.expect(executor.evaluate(exact(0)).status == plot::PlotNumericStatus::Undefined
                && executor.evaluate(exact(-1)).status == plot::PlotNumericStatus::NonReal,
            "PlotProgram: distinguishes log zero from negative-real branch exit");
    }

    const auto explicitDegree = compile("180Deg");
    tests.expect(static_cast<bool>(explicitDegree),
        "PlotProgram: explicit angle units compile as unit conversion instructions");
    if (explicitDegree) {
        plot::BigFloatPlotExecutor executor{*explicitDegree.program, 96};
        const auto result = executor.evaluate(exact(0));
        tests.expect(result.status == plot::PlotNumericStatus::Finite
                && near(
                    result.value,
                    Rational{BigInt{314159}, BigInt{100000}},
                    Rational{BigInt{1}, BigInt{100000}}),
            "PlotProgram: 180Deg converts to Pi in the default radian coordinate");
    }

    const auto explicitDegreeVariable = compile("x Deg");
    tests.expect(explicitDegreeVariable
            && explicitDegreeVariable.program->geometryKind == plot::PlotCurveGeometryKind::Affine,
        "PlotProgram: variable angle-unit conversion preserves affine geometry");
    if (explicitDegreeVariable) {
        plot::BigFloatPlotExecutor executor{*explicitDegreeVariable.program, 96};
        const auto result = executor.evaluate(exact(180));
        tests.expect(result.status == plot::PlotNumericStatus::Finite
                && near(
                    result.value,
                    Rational{BigInt{314159}, BigInt{100000}},
                    Rational{BigInt{1}, BigInt{100000}}),
            "PlotProgram: x Deg converts variable magnitudes into the active angle unit");
    }

    const auto degreeSine = compile("sin[x]");
    if (degreeSine) {
        plot::BigFloatPlotExecutor executor{
            *degreeSine.program, 96,
            mathematics::AngleSemantics{mathematics::AngleUnit::Degree}};
        const auto result = executor.evaluate(exact(30));
        tests.expect(result.status == plot::PlotNumericStatus::Finite
                && near(
                    result.value,
                    Rational{BigInt{1}, BigInt{2}},
                    Rational{BigInt{1}, BigInt{1000000}}),
            "PlotProgram: trigonometric executor respects angle semantics");
    }


    const auto extendedElementary = compile(
        "cbrt[x]+expm1[x]+log1p[x]+asin[x]+acos[x]+atan[x]+asinh[x]+acosh[x+2]+atanh[x/2]+erf[x]+erfc[x]");
    tests.expect(static_cast<bool>(extendedElementary),
        "PlotProgram: compiles the certified extended elementary-function set");
    if (extendedElementary) {
        plot::BigFloatPlotExecutor executor{*extendedElementary.program, 96};
        tests.expect(executor.evaluate(exact(0)).status == plot::PlotNumericStatus::Finite,
            "PlotProgram: extended elementary functions evaluate through the dedicated IR");
    }

    const auto logarithmBases = compile("log2[x]+log10[x]+log[2,x]");
    tests.expect(static_cast<bool>(logarithmBases),
        "PlotProgram: base-2, base-10, and explicit-base logarithms compile");
    if (logarithmBases) {
        plot::BigFloatPlotExecutor executor{*logarithmBases.program, 96};
        const auto value = executor.evaluate(exact(100));
        tests.expect(value.status == plot::PlotNumericStatus::Finite,
            "PlotProgram: logarithm base variants evaluate at positive real samples");
    }

    const auto reciprocalTrig = compile("cot[x]+sec[x]+csc[x]");
    tests.expect(static_cast<bool>(reciprocalTrig),
        "PlotProgram: reciprocal trigonometric functions compile without generic evaluation");

    const auto specialReal = compile("fresnelc[x]+fresnels[x]");
    tests.expect(static_cast<bool>(specialReal),
        "PlotProgram: real Fresnel functions reuse certified point backends");

    const auto normalizedElementary = compile(
        "sinc[x]+cosc[x]+tanc[x]+sinhc[x]+tanhc[x]+expc[x]");
    tests.expect(static_cast<bool>(normalizedElementary),
        "PlotProgram: normalized trigonometric/hyperbolic/exponential functions compile");
    if (normalizedElementary) {
        plot::BigFloatPlotExecutor executor{*normalizedElementary.program, 96};
        const auto zero = executor.evaluate(exact(0));
        tests.expect(zero.status == plot::PlotNumericStatus::Finite
                && zero.value.toRational() == Rational{BigInt{5}},
            "PlotProgram: normalized functions use their removable values at zero without numeric division");
    }

    const auto sineIntegral = compile("Si[x]");
    tests.expect(static_cast<bool>(sineIntegral),
        "PlotProgram: SineIntegral Si compiles through the certified real backend");
    if (sineIntegral) {
        plot::BigFloatPlotExecutor executor{*sineIntegral.program, 96};
        const auto zero = executor.evaluate(exact(0));
        tests.expect(zero.status == plot::PlotNumericStatus::Finite && zero.value.isZero()
                && executor.evaluate(exact(-1)).status == plot::PlotNumericStatus::Finite
                && executor.evaluate(exact(1)).status == plot::PlotNumericStatus::Finite,
            "PlotProgram: Si remains real on the whole real axis and keeps Si(0)=0");
    }

    const auto exponentialIntegral = compile("Ei[x]");
    tests.expect(static_cast<bool>(exponentialIntegral),
        "PlotProgram: ExponentialIntegral Ei compiles through the certified real backend");
    if (exponentialIntegral) {
        plot::BigFloatPlotExecutor executor{*exponentialIntegral.program, 96};
        tests.expect(executor.evaluate(exact(-1)).status == plot::PlotNumericStatus::Finite
                && executor.evaluate(exact(0)).status == plot::PlotNumericStatus::Undefined
                && executor.evaluate(exact(1)).status == plot::PlotNumericStatus::Finite,
            "PlotProgram: real Ei is finite on both half-axes and undefined at zero");
    }

    const auto cosineIntegral = compile("Ci[x]");
    tests.expect(static_cast<bool>(cosineIntegral),
        "PlotProgram: CosineIntegral Ci compiles through the positive-real certified backend");
    if (cosineIntegral) {
        plot::BigFloatPlotExecutor executor{*cosineIntegral.program, 96};
        tests.expect(executor.evaluate(exact(-1)).status == plot::PlotNumericStatus::NonReal
                && executor.evaluate(exact(0)).status == plot::PlotNumericStatus::Undefined
                && executor.evaluate(exact(1)).status == plot::PlotNumericStatus::Finite,
            "PlotProgram: principal Ci rejects the negative real branch and zero in real Plot");
    }

    const auto degreeInverse = compile("asin[x]");
    if (degreeInverse) {
        plot::BigFloatPlotExecutor executor{
            *degreeInverse.program, 96,
            mathematics::AngleSemantics{mathematics::AngleUnit::Degree}};
        const auto result = executor.evaluate(exact(1));
        tests.expect(result.status == plot::PlotNumericStatus::Finite
                && near(result.value, Rational{BigInt{90}}, Rational{BigInt{1}, BigInt{1000000}}),
            "PlotProgram: inverse trigonometric outputs respect degree semantics");
    }

    const auto discreteReal = compile("floor[x]+ceil[x]+sign[x]");
    tests.expect(static_cast<bool>(discreteReal),
        "PlotProgram: floor, ceil, and real sign compile into the dedicated IR");
    if (discreteReal) {
        plot::BigFloatPlotExecutor executor{*discreteReal.program, 96};
        const auto negative = executor.evaluate(BigFloat::fromRational(
            Rational{BigInt{-3}, BigInt{2}}, 96, RoundingMode::NearestEven));
        tests.expect(negative.status == plot::PlotNumericStatus::Finite
                && negative.value.toRational() == Rational{BigInt{-4}},
            "PlotProgram: floor/ceil/sign preserve exact real step semantics");
        const auto zero = executor.evaluate(exact(0));
        tests.expect(zero.status == plot::PlotNumericStatus::Finite
                && zero.value.isZero(),
            "PlotProgram: sign zero remains exactly zero");
    }


    const auto lowCostDiscrete = compile("trunc[x]+round[x]+frac[x]");
    tests.expect(static_cast<bool>(lowCostDiscrete),
        "PlotProgram: trunc/round/frac compile into the dedicated real Plot IR");
    if (lowCostDiscrete) {
        plot::BigFloatPlotExecutor executor{*lowCostDiscrete.program, 96};
        const auto positive = executor.evaluate(BigFloat::fromRational(
            Rational{BigInt{3}, BigInt{2}}, 96, RoundingMode::NearestEven));
        const auto negative = executor.evaluate(BigFloat::fromRational(
            Rational{BigInt{-3}, BigInt{2}}, 96, RoundingMode::NearestEven));
        tests.expect(positive.status == plot::PlotNumericStatus::Finite
                && positive.value.toRational() == Rational{BigInt{7}, BigInt{2}}
                && negative.status == plot::PlotNumericStatus::Finite
                && negative.value.toRational() == Rational{BigInt{-5}, BigInt{2}},
            "PlotProgram: trunc/nearest-even round/frac preserve exact step and fractional semantics");
    }

    const auto piConstant = compile("Pi+x");
    tests.expect(static_cast<bool>(piConstant),
        "PlotProgram: compiles registered mathematical constants");
    if (piConstant) {
        plot::BigFloatPlotExecutor executor{*piConstant.program, 96};
        const auto result = executor.evaluate(exact(0));
        tests.expect(result.status == plot::PlotNumericStatus::Finite
                && near(
                    result.value,
                    Rational{BigInt{314159}, BigInt{100000}},
                    Rational{BigInt{1}, BigInt{100000}}),
            "PlotProgram: numerical constant pool is available to the executor");
    }

    const auto unsupportedFunction = compile("gamma[x]");
    tests.expect(unsupportedFunction.status == plot::PlotCompileStatus::UnsupportedExpression,
        "PlotProgram: unsupported functions fail compilation instead of re-entering evaluator");

    const auto realPower = compile("2^x");
    tests.expect(static_cast<bool>(realPower),
        "PlotProgram: positive-base variable exponents compile to principal real Power");
    if (realPower) {
        plot::BigFloatPlotExecutor executor{*realPower.program, 96};
        const auto value = executor.evaluate(exact(3));
        tests.expect(value.status == plot::PlotNumericStatus::Finite
                && near(value.value, Rational{BigInt{8}}, Rational{BigInt{1}, BigInt{1000000}}),
            "PlotProgram: positive-base variable Power evaluates through certified log/exp");
    }

    const auto rationalPower = compile("x^(3/2)");
    tests.expect(static_cast<bool>(rationalPower),
        "PlotProgram: exact noninteger rational exponents compile");
    if (rationalPower) {
        plot::BigFloatPlotExecutor executor{*rationalPower.program, 96};
        const auto positive = executor.evaluate(exact(4));
        tests.expect(positive.status == plot::PlotNumericStatus::Finite
                && near(positive.value, Rational{BigInt{8}}, Rational{BigInt{1}, BigInt{1000000}})
                && executor.evaluate(exact(-1)).status == plot::PlotNumericStatus::NonReal
                && executor.evaluate(exact(0)).status == plot::PlotNumericStatus::Finite,
            "PlotProgram: principal rational Power keeps negative-base branch exits explicit");
    }

    const auto irrationalPower = compile("x^Pi");
    tests.expect(static_cast<bool>(irrationalPower),
        "PlotProgram: invariant irrational exponents compile instead of requiring integer Power");
    if (irrationalPower) {
        plot::BigFloatPlotExecutor executor{*irrationalPower.program, 96};
        tests.expect(executor.evaluate(exact(1)).status == plot::PlotNumericStatus::Finite
                && executor.evaluate(exact(-1)).status == plot::PlotNumericStatus::NonReal,
            "PlotProgram: irrational principal Power is real only on the proven real branch");
    }

    const auto integerValuedPower = compile("(-2)^floor[x]");
    tests.expect(static_cast<bool>(integerValuedPower),
        "PlotProgram: integer-valued dynamic exponents compile");
    if (integerValuedPower) {
        plot::BigFloatPlotExecutor executor{*integerValuedPower.program, 96};
        const auto half = BigFloat::fromRational(
            Rational{BigInt{3}, BigInt{2}}, 96, RoundingMode::NearestEven);
        const auto value = executor.evaluate(half);
        tests.expect(value.status == plot::PlotNumericStatus::Finite
                && value.value.toRational() == Rational{BigInt{-2}},
            "PlotProgram: negative bases accept runtime exponents only when they are exact integers");
    }
}

} // namespace mmcal::tests
