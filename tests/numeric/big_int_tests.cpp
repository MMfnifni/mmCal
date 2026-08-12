// 任意精度整数BigIntの回帰テスト
#include "big_int_tests.hpp"

#include "numeric/big_int.hpp"
#include "test_framework.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace mmcal::tests {
namespace {

using mmcal::numeric::BigInt;
using mmcal::numeric::divmod;

void testConstructionAndFormatting(TestRunner& tests) {
    tests.expectEqual(BigInt{}.toString(), std::string{"0"}, "BigInt default is zero");
    tests.expectEqual(BigInt{42}.toString(), std::string{"42"}, "BigInt positive construction");
    tests.expectEqual(BigInt{-42}.toString(), std::string{"-42"}, "BigInt negative construction");
    tests.expectEqual(
        BigInt{std::numeric_limits<std::int64_t>::min()}.toString(),
        std::string{"-9223372036854775808"},
        "BigInt constructs INT64_MIN without overflow");

    tests.expectEqual(
        BigInt::parse("-123456789012345678901234567890").toString(),
        std::string{"-123456789012345678901234567890"},
        "BigInt parses arbitrary negative values");
    tests.expectEqual(BigInt::parse("-0").toString(), std::string{"0"}, "BigInt normalizes negative zero");
    tests.expectEqual(
        BigInt::fromUnsigned(UINT64_MAX).toString(),
        std::string{"18446744073709551615"},
        "BigInt constructs the full uint64 range without decimal parsing");
    tests.expectEqual(BigInt::parse("-FF", 16).toString(16), std::string{"-FF"}, "BigInt signed radix conversion");

    tests.expectThrows<std::invalid_argument>(
        [] { (void)BigInt::parse("-"); },
        "BigInt rejects a sign without digits");
}

void testComparisonAndUnarySign(TestRunner& tests) {
    const auto large = BigInt::parse("123456789012345678901234567890");
    const auto negativeLarge = -large;

    tests.expect(negativeLarge < BigInt{-1}, "BigInt compares negative magnitudes in reverse");
    tests.expect(BigInt{-1} < BigInt{}, "BigInt negative is less than zero");
    tests.expect(BigInt{} < BigInt{1}, "BigInt zero is less than positive");
    tests.expect(-negativeLarge == large, "BigInt double negation");
    tests.expect((-BigInt{}).isZero(), "BigInt negated zero remains zero");
    tests.expect(large.abs() == large, "BigInt positive absolute value");
    tests.expect(negativeLarge.abs() == large, "BigInt negative absolute value");
}

void testArithmetic(TestRunner& tests) {
    const auto a = BigInt::parse("123456789012345678901234567890");
    const auto b = BigInt::parse("987654321098765432109876543210");

    tests.expectEqual((a + b).toString(), std::string{"1111111110111111111011111111100"}, "BigInt large addition");
    tests.expectEqual((a - b).toString(), std::string{"-864197532086419753208641975320"}, "BigInt positive subtraction crossing zero");
    tests.expectEqual((-a + b).toString(), std::string{"864197532086419753208641975320"}, "BigInt mixed-sign addition");
    tests.expectEqual((-a - b).toString(), std::string{"-1111111110111111111011111111100"}, "BigInt negative subtraction");
    tests.expectEqual((BigInt{-7} * BigInt{-9}).toString(), std::string{"63"}, "BigInt negative multiplication");
    tests.expectEqual((BigInt{-7} * BigInt{9}).toString(), std::string{"-63"}, "BigInt mixed-sign multiplication");

    auto self = BigInt::parse("-12345678901234567890");
    self += self;
    tests.expectEqual(self.toString(), std::string{"-24691357802469135780"}, "BigInt self addition");
    self -= self;
    tests.expect(self.isZero() && !self.isNegative(), "BigInt self subtraction normalizes zero sign");
}

void testLargeMultiplication(TestRunner& tests) {
    constexpr std::size_t balancedBits = 32 * 96;
    BigInt power{1};
    power <<= balancedBits;

    const BigInt below = power - BigInt{1};
    const BigInt above = power + BigInt{1};
    BigInt expected{1};
    expected <<= balancedBits * 2;
    expected -= BigInt{1};

    tests.expect(
        below * above == expected,
        "BigInt Karatsuba multiplication preserves a large difference-of-squares identity");

    const BigInt negative = -below;
    tests.expect(
        negative * above == -expected,
        "BigInt Karatsuba multiplication preserves signed results");

    BigInt self = above;
    self *= self;
    tests.expect(
        self == above * above,
        "BigInt Karatsuba multiplication supports self multiplication");

    constexpr std::size_t largeBits = 32 * 200;
    constexpr std::size_t smallBits = 32 * 20;
    BigInt large{1};
    large <<= largeBits;
    large += BigInt{1};
    BigInt small{1};
    small <<= smallBits;
    small += BigInt{1};

    BigInt unbalancedExpected{1};
    unbalancedExpected <<= largeBits + smallBits;
    BigInt largeTerm{1};
    largeTerm <<= largeBits;
    BigInt smallTerm{1};
    smallTerm <<= smallBits;
    unbalancedExpected += largeTerm;
    unbalancedExpected += smallTerm;
    unbalancedExpected += BigInt{1};

    tests.expect(
        large * small == unbalancedExpected,
        "BigInt unbalanced multiplication keeps the exact schoolbook fallback result");

    constexpr std::size_t toomBits = 32 * 1800;
    BigInt toomPower{1};
    toomPower <<= toomBits;
    const BigInt toomBelow = toomPower - BigInt{1};
    const BigInt toomAbove = toomPower + BigInt{1};
    BigInt toomExpected{1};
    toomExpected <<= toomBits * 2;
    toomExpected -= BigInt{1};
    tests.expect(
        toomBelow * toomAbove == toomExpected,
        "BigInt Toom-3 multiplication preserves a huge difference-of-squares identity");

    // 専用square経路は一般乗算と独立した対称算法なので、巨大値で恒等式を直接固定する。
    constexpr std::size_t squareBits = 32 * 512;
    BigInt squareBase{1};
    squareBase <<= squareBits;
    squareBase += BigInt{1234567};
    BigInt squareExpected{1};
    squareExpected <<= squareBits * 2;
    BigInt squareCross{2469134};
    squareCross <<= squareBits;
    squareExpected += squareCross;
    squareExpected += BigInt{1524155677489LL};
    tests.expect(
        squareBase * squareBase == squareExpected,
        "BigInt dedicated square preserves a 512-limb exact identity");
}

void testDivision(TestRunner& tests) {
    struct Case final {
        std::int64_t dividend;
        std::int64_t divisor;
        std::int64_t quotient;
        std::int64_t remainder;
    };

    constexpr Case cases[] = {
        {7, 3, 2, 1},
        {-7, 3, -2, -1},
        {7, -3, -2, 1},
        {-7, -3, 2, -1},
        {1, 2, 0, 1},
        {-1, 2, 0, -1},
    };

    bool allCasesMatch = true;
    for (const auto& item : cases) {
        const auto result = divmod(BigInt{item.dividend}, BigInt{item.divisor});
        allCasesMatch = allCasesMatch
            && result.quotient == BigInt{item.quotient}
            && result.remainder == BigInt{item.remainder};
    }
    tests.expect(allCasesMatch, "BigInt division truncates toward zero and keeps dividend remainder sign");

    const auto dividend = BigInt::parse("-1234567890123456789012345678901234567890");
    const auto divisor = BigInt::parse("98765432109876543210");
    const auto result = divmod(dividend, divisor);

    tests.expectEqual(result.quotient.toString(), std::string{"-12499999886093750001"}, "BigInt large signed quotient");
    tests.expectEqual(result.remainder.toString(), std::string{"-54205246805420524680"}, "BigInt large signed remainder");
    tests.expect(
        result.quotient * divisor + result.remainder == dividend,
        "BigInt signed division reconstruction");
    tests.expect(result.remainder.abs() < divisor.abs(), "BigInt signed remainder magnitude bound");

    BigInt powerOfTwo{1};
    powerOfTwo <<= 4097;
    BigInt powerDividend = powerOfTwo * BigInt{123456789} + BigInt{987654321};
    const auto powerResult = divmod(powerDividend, powerOfTwo);
    tests.expectEqual(powerResult.quotient.toString(), std::string{"123456789"},
        "BigInt power-of-two division uses the exact shift quotient");
    tests.expectEqual(powerResult.remainder.toString(), std::string{"987654321"},
        "BigInt power-of-two division preserves low-bit remainder");

    tests.expectThrows<std::domain_error>(
        [] { (void)(BigInt{1} / BigInt{}); },
        "BigInt rejects division by zero");
}

void testAgainstInt64(TestRunner& tests) {
    bool matches = true;

    for (std::int64_t lhs = -125; lhs <= 125 && matches; ++lhs) {
        for (std::int64_t rhs = -31; rhs <= 31; ++rhs) {
            const BigInt a{lhs};
            const BigInt b{rhs};

            matches = matches
                && (a + b).toString() == std::to_string(lhs + rhs)
                && (a - b).toString() == std::to_string(lhs - rhs)
                && (a * b).toString() == std::to_string(lhs * rhs);

            if (rhs != 0)
                matches = matches
                    && (a / b).toString() == std::to_string(lhs / rhs)
                    && (a % b).toString() == std::to_string(lhs % rhs);

            if (!matches)
                break;
        }
    }

    tests.expect(matches, "BigInt arithmetic agrees with int64_t over a signed test grid");
}

} // namespace

void runBigIntTests(TestRunner& tests) {
    testConstructionAndFormatting(tests);
    testComparisonAndUnarySign(tests);
    testArithmetic(tests);
    testLargeMultiplication(tests);
    testDivision(tests);
    testAgainstInt64(tests);
}

} // namespace mmcal::tests
