// 三角函数のexact特殊値の回帰テスト
#include "exact_trigonometry_tests.hpp"

#include "builtins/names.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "formatting/expr_formatter.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/exact_trigonometry.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

namespace mmcal::tests {
namespace {

expression::Expr integer(std::int64_t value) {
    return expression::Expr{numeric::Number{numeric::BigInt{value}}};
}

expression::Expr rational(std::int64_t numerator, std::int64_t denominator) {
    return expression::Expr{numeric::Number{numeric::Rational{
        numeric::BigInt{numerator}, numeric::BigInt{denominator}}}};
}

} // namespace

void runExactTrigonometryTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const evaluation::BuiltinRegistry builtins = evaluation::BuiltinRegistry::defaults(symbols);
    const mathematics::MathRegistry mathematics = mathematics::MathRegistry::defaults(symbols, builtins);
    const mathematics::AngleSemantics defaultAngles;
    const mathematics::AngleSemantics angles{mathematics::AngleUnit::Degree};

    const expression::Expr pi{symbols.intern("Pi")};
    const expression::Expr twoPi = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Multiply), {integer(2), pi});
    const expression::Expr piOverSix = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Divide), {pi, integer(6)});

    tests.expect(mathematics::extractRationalPiMultiple(pi, builtins, mathematics)->toString() == "1",
        "ExactTrig: recognizes Pi coefficient");
    tests.expect(mathematics::extractRationalPiMultiple(twoPi, builtins, mathematics)->toString() == "2",
        "ExactTrig: recognizes numeric Pi multiple");
    tests.expect(mathematics::extractRationalPiMultiple(piOverSix, builtins, mathematics)->toString() == "1/6",
        "ExactTrig: recognizes divided Pi multiple");

    const auto sin30 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Sin, integer(30), builtins, mathematics, angles);
    tests.expect(sin30 && formatting::formatExpr(*sin30) == "1/2",
        "ExactTrig: explicit Degree semantics keeps 30 degrees exact");

    const auto cos180 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Cos, integer(180), builtins, mathematics, angles);
    tests.expect(cos180 && formatting::formatExpr(*cos180) == "-1",
        "ExactTrig: cos 180 degrees is exact -1");


    const auto sin15 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Sin, integer(15), builtins, mathematics, angles);
    tests.expect(sin15 && formatting::formatExpr(*sin15) == "(sqrt[6]-sqrt[2])/4",
        "ExactTrig: sin 15 degrees is represented by exact radicals");

    const auto cos75 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Cos, integer(75), builtins, mathematics, angles);
    tests.expect(cos75 && formatting::formatExpr(*cos75) == "(sqrt[6]-sqrt[2])/4",
        "ExactTrig: cos 75 degrees shares the same exact radical value");

    const auto cos15 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Cos, integer(15), builtins, mathematics, angles);
    tests.expect(cos15 && formatting::formatExpr(*cos15) == "(sqrt[6]+sqrt[2])/4",
        "ExactTrig: cos 15 degrees is represented by exact radicals");

    const auto tan45 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Tan, integer(45), builtins, mathematics, angles);
    tests.expect(tan45 && formatting::formatExpr(*tan45) == "1",
        "ExactTrig: tan 45 degrees is exact one");

    const auto tan30 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Tan, integer(30), builtins, mathematics, angles);
    tests.expect(tan30 && formatting::formatExpr(*tan30) == "sqrt[3]/3",
        "ExactTrig: tan 30 degrees uses an exact radical");

    const auto tan135 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Tan, integer(135), builtins, mathematics, angles);
    tests.expect(tan135 && formatting::formatExpr(*tan135) == "-1",
        "ExactTrig: tan reduction preserves quadrant sign and half-turn periodicity");


    const auto cot90 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Cot, integer(90), builtins, mathematics, angles);
    tests.expect(cot90 && formatting::formatExpr(*cot90) == "0",
        "ExactTrig: cot 90 degrees is exact zero");

    const auto sec60 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Sec, integer(60), builtins, mathematics, angles);
    tests.expect(sec60 && formatting::formatExpr(*sec60) == "2",
        "ExactTrig: sec 60 degrees is exact two");

    const auto sec45 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Sec, integer(45), builtins, mathematics, angles);
    tests.expect(sec45 && formatting::formatExpr(*sec45) == "sqrt[2]",
        "ExactTrig: sec 45 degrees is rationalized to sqrt two");

    const auto cot60 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Cot, integer(60), builtins, mathematics, angles);
    tests.expect(cot60 && formatting::formatExpr(*cot60) == "sqrt[3]/3",
        "ExactTrig: cot 60 degrees uses the canonical rationalized radical");

    const auto csc30 = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Csc, integer(30), builtins, mathematics, angles);
    tests.expect(csc30 && formatting::formatExpr(*csc30) == "2",
        "ExactTrig: csc 30 degrees is exact two");

    const auto asinHalf = mathematics::simplifyExactInverseTrig(
        mathematics::FunctionId::Asin, rational(1, 2), builtins, mathematics, angles);
    tests.expect(asinHalf && formatting::formatExpr(*asinHalf) == "30",
        "ExactTrig: asin one half returns the default Degree angle");

    const auto sinPiDefault = mathematics::simplifyExactTrig(
        mathematics::FunctionId::Sin, pi, builtins, mathematics, defaultAngles);
    tests.expect(sinPiDefault && formatting::formatExpr(*sinPiDefault) == "0",
        "ExactTrig: default angle semantics is Radian");

    const mathematics::AngleSemantics radians{mathematics::AngleUnit::Radian};
    const auto asinHalfRadians = mathematics::simplifyExactInverseTrig(
        mathematics::FunctionId::Asin, rational(1, 2), builtins, mathematics, radians);
    tests.expect(asinHalfRadians && formatting::formatExpr(*asinHalfRadians) == "Pi/6",
        "ExactTrig: inverse trig exact output follows Radian angle semantics");

    const mathematics::AngleSemantics gradians{mathematics::AngleUnit::Gradian};
    const auto asinHalfGradians = mathematics::simplifyExactInverseTrig(
        mathematics::FunctionId::Asin, rational(1, 2), builtins, mathematics, gradians);
    tests.expect(asinHalfGradians && formatting::formatExpr(*asinHalfGradians) == "100/3",
        "ExactTrig: inverse trig exact output follows Gradian angle semantics");

    const auto atan2NorthWest = mathematics::simplifyExactAtan2(
        integer(1), integer(-1), builtins, mathematics, angles);
    tests.expect(atan2NorthWest && formatting::formatExpr(*atan2NorthWest) == "135",
        "ExactTrig: atan2 preserves the quadrant in Degree mode");
}

} // namespace mmcal::tests
