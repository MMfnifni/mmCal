// exact代数計算の回帰テスト
#include "exact_algebra_tests.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "formatting/expr_formatter.hpp"
#include "mathematics/exact_algebra.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

#include <vector>

namespace mmcal::tests {
namespace {

expression::Expr integer(std::int64_t value) {
    return expression::Expr{numeric::Number{numeric::BigInt{value}}};
}

expression::Expr sqrtExpr(
    std::int64_t value,
    const evaluation::BuiltinRegistry& builtins) {
    return expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sqrt),
        {integer(value)});
}

} // namespace

void runExactAlgebraTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const evaluation::BuiltinRegistry builtins = evaluation::BuiltinRegistry::defaults(symbols);

    const expression::Expr difference = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Subtract),
        {sqrtExpr(6, builtins), sqrtExpr(2, builtins)});
    const expression::Expr quarter = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Divide),
        {difference, integer(4)});

    const expression::Expr doubled = mathematics::scaleExactExpression(
        numeric::Rational{numeric::BigInt{2}},
        quarter,
        builtins);
    tests.expectEqual(
        formatting::formatExpr(doubled),
        std::string{"(sqrt[6]-sqrt[2])/2"},
        "ExactAlgebra: rational scaling reduces an existing exact denominator");

    const std::vector<expression::Expr> duplicateTerms{quarter, quarter};
    const std::vector<expression::Expr> combined =
        mathematics::combineStructurallyIdenticalTerms(duplicateTerms, builtins);
    tests.expect(combined.size() == 1,
        "ExactAlgebra: structurally identical terms are collected");
    tests.expectEqual(
        formatting::formatExpr(combined.front()),
        std::string{"(sqrt[6]-sqrt[2])/2"},
        "ExactAlgebra: collected identical radical terms keep exact form");

    const expression::Expr sqrt2 = sqrtExpr(2, builtins);
    const expression::Expr sqrt3 = sqrtExpr(3, builtins);
    const std::vector<expression::Expr> distinctTerms{sqrt2, sqrt3};
    const std::vector<expression::Expr> distinct =
        mathematics::combineStructurallyIdenticalTerms(distinctTerms, builtins);
    tests.expect(distinct.size() == 2,
        "ExactAlgebra: mathematically unrelated structures are not guessed equal");
}

} // namespace mmcal::tests
