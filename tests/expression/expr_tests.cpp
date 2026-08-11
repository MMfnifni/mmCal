// 式ASTと値表現の回帰テスト
#include "expr_tests.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "formatting/expr_formatter.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <vector>

namespace mmcal::tests {

void runExprTests(TestRunner& tests) {
    using expression::Expr;
    using expression::ExprKind;
    using expression::Symbol;
    using formatting::formatExpr;
    using numeric::BigInt;
    using numeric::Number;
    using numeric::RealNumber;

    const Expr integer{Number{BigInt{42}}};
    tests.expect(integer.kind() == ExprKind::Number, "Expr: number node kind");
    tests.expectEqual(formatExpr(integer), "42", "Expr: formats number node");

    const Expr truth{true};
    tests.expect(truth.kind() == ExprKind::Boolean, "Expr: boolean node kind");
    tests.expect(truth.asBoolean(), "Expr: reads boolean value");
    tests.expectEqual(formatExpr(truth), "True", "Expr: formats boolean node");

    const Expr x{Symbol{"x"}};
    tests.expect(x.kind() == ExprKind::Symbol, "Expr: symbol node kind");
    tests.expectEqual(formatExpr(x), "x", "Expr: formats symbol");

    const Expr vector = Expr::array(
        {3},
        {Expr{Number{BigInt{1}}}, x, Expr{Number{BigInt{3}}}});
    tests.expect(vector.isArray(), "Expr: one-dimensional array");
    tests.expectEqual(vector.asArray().rank(), std::size_t{1}, "Expr: vector rank");
    tests.expectEqual(vector.asArray().size(), std::size_t{3}, "Expr: vector element count");
    tests.expectEqual(formatExpr(vector), "{1, x, 3}", "Expr: formats vector");

    const Expr matrix = Expr::array(
        {2, 2},
        {
            Expr{Number{BigInt{1}}}, Expr{Number{BigInt{2}}},
            Expr{Number{BigInt{3}}}, Expr{Number{BigInt{4}}}
        });
    tests.expectEqual(formatExpr(matrix), "{{1, 2}, {3, 4}}", "Expr: formats matrix");

    const Expr emptyMatrix = Expr::array({2, 0}, {});
    tests.expectEqual(formatExpr(emptyMatrix), "{{}, {}}", "Expr: formats zero-column matrix");

    const Expr sqrtTwo = Expr::call(Symbol{"sqrt"}, {Expr{Number{BigInt{2}}}});
    tests.expect(sqrtTwo.isCall(), "Expr: function call node");
    tests.expectEqual(formatExpr(sqrtTwo), "sqrt[2]", "Expr: formats function call");

    const Expr complex{Number::complex(RealNumber{BigInt{1}}, RealNumber{BigInt{2}})};
    tests.expectEqual(formatExpr(complex), "1+2I", "Expr: formats complex number node");

    tests.expectEqual(
        formatting::trimRedundantFractionalZeros("3.10000"),
        "3.1",
        "Formatter: trims redundant fixed fractional zeros");
    tests.expectEqual(
        formatting::trimRedundantFractionalZeros("{3.0000, -2.5000, 0.33300}"),
        "{3, -2.5, 0.333}",
        "Formatter: trims fixed fractional zeros inside formatted expressions");
    tests.expectEqual(
        formatting::trimRedundantFractionalZeros("{\"3.1000\", 3.1000}"),
        "{\"3.1000\", 3.1}",
        "Formatter: preserves decimal-looking text inside string literals");

    const auto power = Expr::call(
        evaluation::defaultBuiltinRegistry().symbol(evaluation::BuiltinId::Power),
        {Expr{numeric::Number{numeric::BigInt{-8}}},
         Expr{numeric::Number{numeric::Rational{numeric::BigInt{1}, numeric::BigInt{3}}}}});
    tests.expectEqual(
        formatExpr(power),
        "(-8)^(1/3)",
        "Formatter: parenthesizes signed and rational Power operands for round-trip safety");

    tests.expect(x == Expr{Symbol{"x"}}, "Expr: compares equal expression nodes");
    tests.expect(!(x == Expr{Symbol{"y"}}), "Expr: distinguishes different symbols");
    tests.expect(matrix == Expr::array(
        {2, 2},
        {
            Expr{Number{BigInt{1}}}, Expr{Number{BigInt{2}}},
            Expr{Number{BigInt{3}}}, Expr{Number{BigInt{4}}}
        }), "Expr: compares equal array expressions");

    tests.expectThrows<std::invalid_argument>(
        [] { static_cast<void>(Symbol{""}); },
        "Symbol: rejects empty name");
    tests.expectThrows<std::invalid_argument>(
        [] { static_cast<void>(Expr::array({}, {})); },
        "Expr: rejects rank-zero array");
    tests.expectThrows<std::invalid_argument>(
        [] {
            static_cast<void>(Expr::array(
                {2, 2},
                {Expr{Number{BigInt{1}}}}));
        },
        "Expr: rejects shape and element mismatch");
    tests.expectThrows<std::logic_error>(
        [&] { static_cast<void>(x.asNumber()); },
        "Expr: rejects symbol number access");
}

} // namespace mmcal::tests
