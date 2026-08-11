// 任意精度有理数Rationalの回帰テスト
#include "rational_tests.hpp"

#include "numeric/rational.hpp"
#include "test_framework.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>

namespace mmcal::tests {
namespace {

using mmcal::numeric::BigInt;
using mmcal::numeric::Rational;

void testNormalization(TestRunner& tests) {
    tests.expectEqual(Rational{}.toString(), std::string{"0"}, "Rational default is zero");
    tests.expectEqual(Rational{BigInt{2}, BigInt{4}}.toString(), std::string{"1/2"}, "Rational reduces common factors");
    tests.expectEqual(Rational{BigInt{-2}, BigInt{-4}}.toString(), std::string{"1/2"}, "Rational normalizes double negative");
    tests.expectEqual(Rational{BigInt{2}, BigInt{-4}}.toString(), std::string{"-1/2"}, "Rational keeps denominator positive");
    tests.expectEqual(Rational{BigInt{}, BigInt{-999}}.toString(), std::string{"0"}, "Rational canonical zero is zero over one");
    tests.expect(Rational{BigInt{8}, BigInt{4}}.isInteger(), "Rational recognizes normalized integers");

    tests.expectThrows<std::domain_error>(
        [] { (void)Rational{BigInt{1}, BigInt{0}}; },
        "Rational rejects a zero denominator");
}

void testParsing(TestRunner& tests) {
    tests.expectEqual(Rational::parse("0.1").toString(), std::string{"1/10"}, "Rational parses decimal exactly");
    tests.expectEqual(Rational::parse("-12.500").toString(), std::string{"-25/2"}, "Rational parses signed decimal exactly");
    tests.expectEqual(Rational::parse(".125").toString(), std::string{"1/8"}, "Rational accepts an omitted whole part");
    tests.expectEqual(Rational::parse("12.").toString(), std::string{"12"}, "Rational accepts an omitted fractional part");
    tests.expectEqual(Rational::parse("101.01", 2).toString(), std::string{"21/4"}, "Rational parses binary fractions exactly");
    tests.expectEqual(Rational::parse("A.F", 16).toString(16), std::string{"AF/10"}, "Rational parses and formats hexadecimal fraction values");

    tests.expectThrows<std::invalid_argument>(
        [] { (void)Rational::parse("."); },
        "Rational rejects a bare radix point");
    tests.expectThrows<std::invalid_argument>(
        [] { (void)Rational::parse("1.2.3"); },
        "Rational rejects multiple radix points");
}

void testArithmetic(TestRunner& tests) {
    const Rational oneThird{BigInt{1}, BigInt{3}};
    const Rational oneSixth{BigInt{1}, BigInt{6}};

    tests.expectEqual((oneThird + oneSixth).toString(), std::string{"1/2"}, "Rational addition");
    tests.expectEqual((oneThird - oneSixth).toString(), std::string{"1/6"}, "Rational subtraction");
    tests.expectEqual((oneThird * oneSixth).toString(), std::string{"1/18"}, "Rational multiplication");
    tests.expectEqual((oneThird / oneSixth).toString(), std::string{"2"}, "Rational division");
    tests.expectEqual((-oneThird).toString(), std::string{"-1/3"}, "Rational unary negation");

    auto selfAdd = oneThird;
    selfAdd += selfAdd;
    tests.expectEqual(selfAdd.toString(), std::string{"2/3"}, "Rational self addition");

    auto selfMultiply = Rational{BigInt{2}, BigInt{3}};
    selfMultiply *= selfMultiply;
    tests.expectEqual(selfMultiply.toString(), std::string{"4/9"}, "Rational self multiplication");

    auto selfDivide = Rational{BigInt{-7}, BigInt{11}};
    selfDivide /= selfDivide;
    tests.expectEqual(selfDivide.toString(), std::string{"1"}, "Rational self division");

    const auto largeA = Rational{
        BigInt::parse("123456789012345678901234567890"),
        BigInt::parse("98765432109876543210987654321")};
    const auto largeB = Rational{
        BigInt::parse("111111111111111111111111111111"),
        BigInt::parse("222222222222222222222222222222")};
    const auto result = largeA * largeB;
    tests.expect(
        result * Rational{BigInt{2}} == largeA,
        "Rational cross-cancellation preserves a large product");

    tests.expectThrows<std::domain_error>(
        [&] { (void)(oneThird / Rational{}); },
        "Rational rejects division by zero");


    // 演算子内の交差約分最適化が、constructorによる完全normalizeと一致することを
    // 符号・互いに素でない分母を含む決定的な入力列で照合する。
    bool optimizedArithmeticMatchesReference = true;
    for (std::int64_t i = 1; i <= 160; ++i) {
        const BigInt a{i * 37 - 2500};
        const BigInt b{i * 11 + 1};
        const BigInt c{i * 53 - 1700};
        const BigInt d{i * 7 + 3};
        const Rational lhs{a, b};
        const Rational rhs{c, d};

        const Rational addReference{
            lhs.numerator() * rhs.denominator() + rhs.numerator() * lhs.denominator(),
            lhs.denominator() * rhs.denominator()};
        const Rational mulReference{
            lhs.numerator() * rhs.numerator(),
            lhs.denominator() * rhs.denominator()};
        if (!(lhs + rhs == addReference) || !(lhs * rhs == mulReference))
            optimizedArithmeticMatchesReference = false;

        if (!rhs.isZero()) {
            const Rational divReference{
                lhs.numerator() * rhs.denominator(),
                lhs.denominator() * rhs.numerator()};
            if (!(lhs / rhs == divReference))
                optimizedArithmeticMatchesReference = false;
        }
    }
    tests.expect(optimizedArithmeticMatchesReference,
        "Rational optimized arithmetic matches full-normalization references");
}

void testComparisonAndExactDecimals(TestRunner& tests) {
    tests.expect(Rational{BigInt{1}, BigInt{3}} < Rational{BigInt{1}, BigInt{2}}, "Rational exact comparison");
    tests.expect(Rational{BigInt{-2}, BigInt{3}} < Rational{BigInt{-1}, BigInt{2}}, "Rational negative comparison");

    const auto decimalSum = Rational::parse("0.1") + Rational::parse("0.2");
    tests.expectEqual(decimalSum.toString(), std::string{"3/10"}, "decimal addition has no binary floating error");
    tests.expect(decimalSum == Rational::parse("0.3"), "equivalent decimal literals compare equal");
}

} // namespace

void runRationalTests(TestRunner& tests) {
    testNormalization(tests);
    testParsing(tests);
    testArithmetic(tests);
    testComparisonAndExactDecimals(tests);
}

} // namespace mmcal::tests
