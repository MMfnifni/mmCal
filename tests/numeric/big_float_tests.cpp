// 任意精度二進浮動小数の回帰テスト
#include "big_float_tests.hpp"

#include "numeric/big_float.hpp"
#include "numeric/rounding_mode.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>

namespace mmcal::tests {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;
using numeric::RoundingMode;

Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

void expectRational(
    TestRunner& tests,
    const BigFloat& value,
    const Rational& expected,
    std::string_view name) {
    tests.expectEqual(value.toRational().toString(), expected.toString(), name);
}

void testConstructionAndNormalization(TestRunner& tests) {
    const BigFloat zero;
    tests.expect(zero.isZero(), "BigFloat: default value is zero");
    tests.expectEqual(zero.exponent(), BigFloat::exponent_type{0},
        "BigFloat: zero has canonical exponent");

    const auto normalized = BigFloat::fromDyadic(BigInt{12}, -3, 16);
    tests.expectEqual(normalized.significand().toString(), std::string{"3"},
        "BigFloat: removes powers of two from significand");
    tests.expectEqual(normalized.exponent(), BigFloat::exponent_type{-1},
        "BigFloat: normalization transfers powers of two to exponent");
    expectRational(tests, normalized, rational(3, 2),
        "BigFloat: normalized dyadic keeps exact value");

    const auto half = BigFloat::fromRational(rational(1, 2), 8);
    expectRational(tests, half, rational(1, 2),
        "BigFloat: dyadic rational converts exactly");
    tests.expectEqual(half.precisionBits(), std::size_t{8},
        "BigFloat: stores requested binary precision");

    tests.expectThrows<std::invalid_argument>([] {
        static_cast<void>(BigFloat::fromBigInt(BigInt{1}, 0));
    }, "BigFloat: rejects zero-bit precision");
}

void testRoundingModes(TestRunner& tests) {
    // 1/3 = 1.010101... * 2^-2。
    // 4bitでは下側が 5/16、上側が 11/32 になる。
    expectRational(tests,
        BigFloat::fromRational(rational(1, 3), 4, RoundingMode::TowardNegative),
        rational(5, 16),
        "BigFloat: positive directed rounding downward");
    expectRational(tests,
        BigFloat::fromRational(rational(1, 3), 4, RoundingMode::TowardPositive),
        rational(11, 32),
        "BigFloat: positive directed rounding upward");
    expectRational(tests,
        BigFloat::fromRational(rational(1, 3), 4, RoundingMode::NearestEven),
        rational(11, 32),
        "BigFloat: nearest rounding for one third");

    expectRational(tests,
        BigFloat::fromRational(rational(-1, 3), 4, RoundingMode::TowardPositive),
        rational(-5, 16),
        "BigFloat: negative upward rounding moves toward positive infinity");
    expectRational(tests,
        BigFloat::fromRational(rational(-1, 3), 4, RoundingMode::TowardNegative),
        rational(-11, 32),
        "BigFloat: negative downward rounding moves toward negative infinity");

    // 9/8 は precision=3 で 1 と 5/4 のちょうど中点。
    // 1 の保持仮数100bは偶数なので ties-to-even では1を選ぶ。
    expectRational(tests,
        BigFloat::fromRational(rational(9, 8), 3, RoundingMode::NearestEven),
        rational(1),
        "BigFloat: nearest-even chooses even lower midpoint");

    // 11/8 は 5/4(101b) と 3/2(110b) の中点なので、偶数側3/2を選ぶ。
    expectRational(tests,
        BigFloat::fromRational(rational(11, 8), 3, RoundingMode::NearestEven),
        rational(3, 2),
        "BigFloat: nearest-even chooses even upper midpoint");
}

void testArithmetic(TestRunner& tests) {
    const auto threeHalves = BigFloat::fromRational(rational(3, 2), 16);
    const auto oneQuarter = BigFloat::fromRational(rational(1, 4), 16);

    expectRational(tests,
        add(threeHalves, oneQuarter, 16, RoundingMode::NearestEven),
        rational(7, 4),
        "BigFloat: exact addition");
    expectRational(tests,
        subtract(threeHalves, oneQuarter, 16, RoundingMode::NearestEven),
        rational(5, 4),
        "BigFloat: exact subtraction");
    expectRational(tests,
        multiply(threeHalves, oneQuarter, 16, RoundingMode::NearestEven),
        rational(3, 8),
        "BigFloat: exact multiplication");
    expectRational(tests,
        divide(threeHalves, oneQuarter, 16, RoundingMode::NearestEven),
        rational(6),
        "BigFloat: exact division");

    const auto one = BigFloat::fromBigInt(BigInt{1}, 16);
    const auto three = BigFloat::fromBigInt(BigInt{3}, 16);
    const auto lower = divide(one, three, 10, RoundingMode::TowardNegative);
    const auto upper = divide(one, three, 10, RoundingMode::TowardPositive);
    tests.expect(lower.toRational() < rational(1, 3),
        "BigFloat: directed division lower bound is below exact value");
    tests.expect(upper.toRational() > rational(1, 3),
        "BigFloat: directed division upper bound is above exact value");

    tests.expectThrows<std::domain_error>([&] {
        static_cast<void>(divide(one, BigFloat{}, 16, RoundingMode::NearestEven));
    }, "BigFloat: division by zero throws");
}

void testComparison(TestRunner& tests) {
    const auto a = BigFloat::fromDyadic(BigInt{3}, -1, 16);   // 1.5
    const auto b = BigFloat::fromDyadic(BigInt{6}, -2, 16);   // same value, normalized differently at input
    const auto c = BigFloat::fromDyadic(BigInt{7}, -2, 16);   // 1.75
    const auto negative = -a;

    tests.expect(a == b, "BigFloat: equality compares mathematical dyadic value");
    tests.expect(a < c, "BigFloat: exact positive comparison");
    tests.expect(negative < a, "BigFloat: signed comparison");
    tests.expect(BigFloat{} < a, "BigFloat: zero compares below a positive value");
    tests.expect(negative < BigFloat{}, "BigFloat: negative value compares below zero");
}

void testDirectedRoundingProperties(TestRunner& tests) {
    bool allBoundsContainExactValue = true;
    bool allNearestValuesStayInsideDirectedBounds = true;

    for (std::int64_t numerator = -31; numerator <= 31; ++numerator) {
        for (std::int64_t denominator = 1; denominator <= 23; ++denominator) {
            const Rational exact{BigInt{numerator}, BigInt{denominator}};

            for (std::size_t precision = 1; precision <= 18; ++precision) {
                const auto lower = BigFloat::fromRational(
                    exact, precision, RoundingMode::TowardNegative);
                const auto upper = BigFloat::fromRational(
                    exact, precision, RoundingMode::TowardPositive);
                const auto nearest = BigFloat::fromRational(
                    exact, precision, RoundingMode::NearestEven);

                if (lower.toRational() > exact || upper.toRational() < exact)
                    allBoundsContainExactValue = false;
                if (nearest < lower || nearest > upper)
                    allNearestValuesStayInsideDirectedBounds = false;
            }
        }
    }

    tests.expect(allBoundsContainExactValue,
        "BigFloat: directed rational conversion encloses exact values across a grid");
    tests.expect(allNearestValuesStayInsideDirectedBounds,
        "BigFloat: nearest values remain inside directed-rounding bounds");
}


void testDirectedArithmeticProperties(TestRunner& tests) {
    bool addBoundsCorrect = true;
    bool multiplyBoundsCorrect = true;
    bool divideBoundsCorrect = true;

    for (std::int64_t a = -9; a <= 9; ++a) {
        for (std::int64_t b = -9; b <= 9; ++b) {
            const auto lhs = BigFloat::fromDyadic(BigInt{a}, -3, 32);
            const auto rhs = BigFloat::fromDyadic(BigInt{b}, -2, 32);
            const Rational lhsExact = lhs.toRational();
            const Rational rhsExact = rhs.toRational();

            for (std::size_t precision = 1; precision <= 12; ++precision) {
                const auto addLower = add(
                    lhs, rhs, precision, RoundingMode::TowardNegative);
                const auto addUpper = add(
                    lhs, rhs, precision, RoundingMode::TowardPositive);
                const Rational exactSum = lhsExact + rhsExact;
                if (addLower.toRational() > exactSum || addUpper.toRational() < exactSum)
                    addBoundsCorrect = false;

                const auto mulLower = multiply(
                    lhs, rhs, precision, RoundingMode::TowardNegative);
                const auto mulUpper = multiply(
                    lhs, rhs, precision, RoundingMode::TowardPositive);
                const Rational exactProduct = lhsExact * rhsExact;
                if (mulLower.toRational() > exactProduct || mulUpper.toRational() < exactProduct)
                    multiplyBoundsCorrect = false;

                if (!rhs.isZero()) {
                    const auto divLower = divide(
                        lhs, rhs, precision, RoundingMode::TowardNegative);
                    const auto divUpper = divide(
                        lhs, rhs, precision, RoundingMode::TowardPositive);
                    const Rational exactQuotient = lhsExact / rhsExact;
                    if (divLower.toRational() > exactQuotient
                        || divUpper.toRational() < exactQuotient)
                        divideBoundsCorrect = false;
                }
            }
        }
    }

    tests.expect(addBoundsCorrect,
        "BigFloat: directed addition encloses exact dyadic sums across a grid");
    tests.expect(multiplyBoundsCorrect,
        "BigFloat: directed multiplication encloses exact dyadic products across a grid");
    tests.expect(divideBoundsCorrect,
        "BigFloat: directed division encloses exact dyadic quotients across a grid");
}

void testIntegerRounding(TestRunner& tests) {
    // 13 = 1101b。3bitへ丸めると12と14の中点で、110bが偶数なので12。
    expectRational(tests,
        BigFloat::fromBigInt(BigInt{13}, 3, RoundingMode::NearestEven),
        rational(12),
        "BigFloat: integer midpoint uses ties-to-even");
    expectRational(tests,
        BigFloat::fromBigInt(BigInt{13}, 3, RoundingMode::TowardPositive),
        rational(14),
        "BigFloat: integer upward rounding");
    expectRational(tests,
        BigFloat::fromBigInt(BigInt{-13}, 3, RoundingMode::TowardZero),
        rational(-12),
        "BigFloat: integer rounding toward zero");
}

} // namespace

void runBigFloatTests(TestRunner& tests) {
    testConstructionAndNormalization(tests);
    testRoundingModes(tests);
    testArithmetic(tests);
    testComparison(tests);
    testDirectedRoundingProperties(tests);
    testDirectedArithmeticProperties(tests);
    testIntegerRounding(tests);
}

} // namespace mmcal::tests
