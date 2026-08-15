// 式ASTと値表現の回帰テスト
#include "expr_tests.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "formatting/expr_formatter.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "test_framework.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace mmcal::tests {

void runExprTests(TestRunner& tests) {
    using expression::ArrayBuilder;
    using expression::ArrayStorageKind;
    using expression::Expr;
    using expression::ExprKind;
    using expression::Symbol;
    using formatting::formatExpr;
    using numeric::BigInt;
    using numeric::Number;
    using numeric::Rational;
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
    tests.expect(matrix.asArray().storageKind() == ArrayStorageKind::Integer,
        "Expr: exact integer array uses packed integer storage");
    tests.expectEqual(
        matrix.asArray().size(), std::size_t{4},
        "Expr: packed integer storage exposes all elements");
    tests.expectEqual(
        matrix.asArray().integerAt(3).toString(), std::string{"4"},
        "Expr: packed integer page exposes stored values");
    tests.expectEqual(
        formatExpr(matrix.asArray().element(2)), "3",
        "Expr: packed integer element materializes with unchanged value");

    const Expr rationalVector = Expr::array(
        {2},
        {Expr{Number{BigInt{1}}},
         Expr{Number{numeric::Rational{BigInt{1}, BigInt{2}}}}});
    tests.expect(rationalVector.asArray().storageKind() == ArrayStorageKind::Rational,
        "Expr: mixed integer and rational array promotes to packed rational storage");

    const Expr complexVector = Expr::array(
        {2},
        {Expr{Number{BigInt{1}}},
         Expr{Number::complex(RealNumber{BigInt{0}}, RealNumber{BigInt{1}})}});
    tests.expect(complexVector.asArray().storageKind() == ArrayStorageKind::Number,
        "Expr: exact complex array uses packed Number storage");

    const Expr genericVector = Expr::array(
        {2}, {Expr{Number{BigInt{1}}}, x});
    tests.expect(genericVector.asArray().storageKind() == ArrayStorageKind::Generic,
        "Expr: symbolic array keeps generic Expr storage");

    ArrayBuilder lateSymbolBuilder;
    lateSymbolBuilder.reserve(2051);
    for (std::size_t i = 0; i < 2050; ++i)
        lateSymbolBuilder.append(BigInt{static_cast<std::int64_t>(i + 1)});
    lateSymbolBuilder.append(Expr{Symbol{"late"}});
    const auto lateSymbolArray = lateSymbolBuilder.finish({2051});
    tests.expect(lateSymbolArray.storageKind() == ArrayStorageKind::Generic,
        "Expr: late symbol widens the logical Array domain to Generic");
    tests.expect(lateSymbolArray.storedKindAt(0) == ArrayStorageKind::Integer
            && lateSymbolArray.storedKindAt(2047) == ArrayStorageKind::Integer,
        "Expr: late symbol does not rebuild completed integer pages");
    tests.expect(lateSymbolArray.storedKindAt(2048) == ArrayStorageKind::Generic
            && lateSymbolArray.storedKindAt(2050) == ArrayStorageKind::Generic,
        "Expr: late symbol promotes only the unfinished page to Generic");

    ArrayBuilder lateRationalBuilder;
    lateRationalBuilder.reserve(2049);
    for (std::size_t i = 0; i < 2048; ++i)
        lateRationalBuilder.append(BigInt{static_cast<std::int64_t>(i + 1)});
    lateRationalBuilder.append(Rational{BigInt{1}, BigInt{2}});
    const auto lateRationalArray = lateRationalBuilder.finish({2049});
    tests.expect(lateRationalArray.storageKind() == ArrayStorageKind::Rational,
        "Expr: exact pages expose their narrowest common Rational domain");
    tests.expect(lateRationalArray.storedKindAt(1024) == ArrayStorageKind::Integer
            && lateRationalArray.storedKindAt(2048) == ArrayStorageKind::Rational,
        "Expr: exact promotion across a page boundary preserves completed integer pages");

    const Expr integerRow = Expr::array(
        {2}, {Expr{Number{BigInt{1}}}, Expr{Number{BigInt{2}}}});
    const Expr rationalRow = Expr::array(
        {2},
        {Expr{Number{numeric::Rational{BigInt{1}, BigInt{2}}}},
         Expr{Number{BigInt{3}}}});
    const Expr promotedMatrix = Expr::array({2}, {integerRow, rationalRow});
    tests.expect(promotedMatrix.asArray().shape == std::vector<std::size_t>({2, 2}),
        "Expr: nested packed rows preserve rectangular shape");
    tests.expect(promotedMatrix.asArray().storageKind() == ArrayStorageKind::Rational,
        "Expr: nested packed rows promote exact storage without materializing Expr nodes");
    tests.expectEqual(formatExpr(promotedMatrix), "{{1, 2}, {1/2, 3}}",
        "Expr: nested packed promotion preserves formatted value");

    const auto transposed = matrix.asArray().transposed();
    tests.expect(!transposed.isContiguous(),
        "Expr: transpose represents layout with shared non-contiguous strides");
    tests.expectEqual(formatExpr(Expr::array(transposed)), "{{1, 3}, {2, 4}}",
        "Expr: transposed packed view preserves logical row-major order");
    tests.expect(transposed.transposed() == matrix.asArray(),
        "Expr: double transpose restores the original Array value");
    const auto transposedRow = transposed.sliced({2}, 0, 2);
    tests.expectEqual(formatExpr(Expr::array(transposedRow)), "{1, 3}",
        "Expr: packed slice follows non-contiguous transpose strides");
    const auto transposedReshaped = transposed.reshaped({4});
    tests.expectEqual(formatExpr(Expr::array(transposedReshaped)), "{1, 3, 2, 4}",
        "Expr: reshape materializes non-contiguous packed views in logical order");

    const auto reshaped = matrix.asArray().reshaped({4});
    tests.expect(reshaped.storageKind() == ArrayStorageKind::Integer,
        "Expr: reshape preserves packed storage");
    const auto sliced = matrix.asArray().sliced({2}, 1, 2);
    tests.expect(sliced.storageKind() == ArrayStorageKind::Integer,
        "Expr: slice preserves packed storage");
    tests.expectEqual(formatExpr(Expr::array(sliced)), "{2, 3}",
        "Expr: packed slice preserves values");

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

    const Expr xCopy = x;
    const Expr xIndependent{Symbol{"x"}};
    tests.expect(xCopy.identity() == x.identity(),
        "Expr: copied handles preserve node identity");
    tests.expect(xIndependent.identity() != x.identity(),
        "Expr: independently constructed equal expressions keep distinct node identity");
    tests.expect(x == xIndependent, "Expr: compares equal expression nodes structurally");
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
