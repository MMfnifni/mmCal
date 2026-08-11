// 記号微分Dの回帰テスト
#include "differentiation_tests.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "formatting/expr_formatter.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "symbolic/differentiation.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

#include <string>

namespace mmcal::tests {

void runDifferentiationTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const auto builtins = evaluation::BuiltinRegistry::defaults(symbols);
    const auto mathematics = mathematics::MathRegistry::defaults(symbols, builtins);
    const expression::Symbol x = symbols.intern("x");

    const expression::Expr bareSin = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sin), {expression::Expr{x}});

    const expression::Expr degreeDerivative = symbolic::differentiateExpression(
        bareSin, x, builtins, mathematics,
        mathematics::AngleSemantics{mathematics::AngleUnit::Degree});
    tests.expectEqual(
        formatting::formatExpr(degreeDerivative),
        std::string{"Pi/180cos[x]"},
        "Differentiation: bare trig uses Degree scale in Degree mode");

    const expression::Expr radianDerivative = symbolic::differentiateExpression(
        bareSin, x, builtins, mathematics,
        mathematics::AngleSemantics{mathematics::AngleUnit::Radian});
    tests.expectEqual(
        formatting::formatExpr(radianDerivative),
        std::string{"cos[x]"},
        "Differentiation: bare trig has unit scale in Radian mode");

    const expression::Expr explicitRadian = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::UnitApplied),
        {expression::Expr{x}, expression::Expr{std::string{"Rad"}}});
    const expression::Expr explicitSin = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sin), {explicitRadian});
    const expression::Expr explicitDerivative = symbolic::differentiateExpression(
        explicitSin, x, builtins, mathematics,
        mathematics::AngleSemantics{mathematics::AngleUnit::Degree});
    tests.expectEqual(
        formatting::formatExpr(explicitDerivative),
        std::string{"cos[x Rad]"},
        "Differentiation: explicit Radian overrides the Degree default");
}

} // namespace mmcal::tests
