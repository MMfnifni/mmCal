// 整数平方根・階乗などの整数算法の回帰テスト
#include "integer_algorithms_tests.hpp"

#include "numeric/integer_algorithms.hpp"
#include "test_framework.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>

namespace mmcal::tests {
namespace {

using mmcal::numeric::BigInt;
using mmcal::numeric::gcd;
using mmcal::numeric::factorial;
using mmcal::numeric::integerCubeRoot;
using mmcal::numeric::integerSqrt;
using mmcal::numeric::isPerfectSquare;
using mmcal::numeric::lcm;
using mmcal::numeric::pow;

void testGcdAndLcm(TestRunner& tests) {
    tests.expectEqual(gcd(BigInt{}, BigInt{}).toString(), std::string{"0"}, "gcd zero zero is zero");
    tests.expectEqual(gcd(BigInt{-84}, BigInt{30}).toString(), std::string{"6"}, "gcd ignores signs");
    tests.expectEqual(lcm(BigInt{-21}, BigInt{6}).toString(), std::string{"42"}, "lcm is nonnegative");
    tests.expectEqual(lcm(BigInt{}, BigInt{99}).toString(), std::string{"0"}, "lcm with zero is zero");

    const auto a = BigInt::parse("123456789012345678901234567890");
    const auto b = BigInt::parse("98765432109876543210");
    const auto factor = gcd(a, b);
    tests.expect((a % factor).isZero() && (b % factor).isZero(), "large gcd divides both operands");
}

void testPower(TestRunner& tests) {
    tests.expectEqual(pow(BigInt{0}, 0).toString(), std::string{"1"}, "integer pow defines zero to zero as one");
    tests.expectEqual(pow(BigInt{-2}, 9).toString(), std::string{"-512"}, "integer pow odd negative exponent parity");
    tests.expectEqual(pow(BigInt{-2}, 10).toString(), std::string{"1024"}, "integer pow even negative exponent parity");
    tests.expectEqual(
        pow(BigInt{2}, 256).toString(),
        std::string{"115792089237316195423570985008687907853269984665640564039457584007913129639936"},
        "integer pow handles arbitrary-size results");
}

void testFactorial(TestRunner& tests) {
    tests.expectEqual(factorial(0).toString(), std::string{"1"},
        "integer factorial defines zero factorial exactly");
    tests.expectEqual(factorial(10).toString(), std::string{"3628800"},
        "integer factorial product tree returns the expected exact value");
    tests.expectEqual(factorial(100).toString(),
        std::string{"933262154439441526816992388562667004907159682643816214685929"
                    "638952175999932299156089414639761565182862536979208272237582"
                    "51185210916864000000000000000000000000"},
        "integer factorial supports arbitrary-precision results");
}

void testUint64Conversion(TestRunner& tests) {
    const auto maximum = tryToUint64(BigInt::parse("18446744073709551615"));
    tests.expect(maximum && *maximum == UINT64_MAX,
        "BigInt uint64 conversion accepts the exact maximum");
    tests.expect(!tryToUint64(BigInt::parse("18446744073709551616")),
        "BigInt uint64 conversion rejects overflow");
    tests.expect(!tryToUint64(BigInt{-1}),
        "BigInt uint64 conversion rejects negative values");
}

void testIntegerSquareRoot(TestRunner& tests) {
    const auto exact = integerSqrt(BigInt::parse("15241578750190521"));
    tests.expectEqual(exact.root.toString(), std::string{"123456789"}, "integer sqrt exact root");
    tests.expect(exact.remainder.isZero(), "integer sqrt exact remainder");

    const auto inexact = integerSqrt(BigInt{27});
    tests.expectEqual(inexact.root.toString(), std::string{"5"}, "integer sqrt floors inexact root");
    tests.expectEqual(inexact.remainder.toString(), std::string{"2"}, "integer sqrt reports residual");

    const auto hugeRoot = BigInt::parse("123456789012345678901234567890");
    const auto hugeValue = hugeRoot * hugeRoot + BigInt{12345};
    const auto hugeResult = integerSqrt(hugeValue);
    tests.expect(hugeResult.root == hugeRoot, "integer sqrt large root");
    tests.expectEqual(hugeResult.remainder.toString(), std::string{"12345"}, "integer sqrt large residual");

    bool invariantHolds = true;
    for (std::int64_t value = 0; value <= 5000; ++value) {
        const auto result = integerSqrt(BigInt{value});
        const auto next = result.root + BigInt{1};
        invariantHolds = invariantHolds
            && result.root * result.root + result.remainder == BigInt{value}
            && result.remainder >= BigInt{}
            && next * next > BigInt{value};

        if (!invariantHolds)
            break;
    }
    tests.expect(invariantHolds, "integer sqrt invariants hold over a test range");

    tests.expect(isPerfectSquare(BigInt{144}), "perfect-square detection accepts a square");
    tests.expect(!isPerfectSquare(BigInt{145}), "perfect-square detection rejects a nonsquare");
    tests.expect(!isPerfectSquare(BigInt{-1}), "perfect-square detection rejects negatives");

    tests.expectThrows<std::domain_error>(
        [] { (void)integerSqrt(BigInt{-1}); },
        "integer sqrt rejects negatives");
}


void testIntegerCubeRoot(TestRunner& tests) {
    const auto exact = integerCubeRoot(BigInt::parse("1881676371789154860897069"));
    tests.expectEqual(exact.root.toString(), std::string{"123456789"},
        "integer cube root exact root");
    tests.expect(exact.remainder.isZero(), "integer cube root exact remainder");

    const auto inexact = integerCubeRoot(BigInt{30});
    tests.expectEqual(inexact.root.toString(), std::string{"3"},
        "integer cube root floors inexact root");
    tests.expectEqual(inexact.remainder.toString(), std::string{"3"},
        "integer cube root reports residual");

    bool invariantHolds = true;
    for (std::int64_t value = 0; value <= 5000; ++value) {
        const auto result = integerCubeRoot(BigInt{value});
        const auto next = result.root + BigInt{1};
        invariantHolds = invariantHolds
            && result.root * result.root * result.root + result.remainder == BigInt{value}
            && result.remainder >= BigInt{}
            && next * next * next > BigInt{value};
        if (!invariantHolds)
            break;
    }
    tests.expect(invariantHolds, "integer cube root invariants hold over a test range");

    tests.expectThrows<std::domain_error>(
        [] { (void)integerCubeRoot(BigInt{-1}); },
        "integer cube root primitive rejects negatives");
}

} // namespace

void runIntegerAlgorithmTests(TestRunner& tests) {
    testGcdAndLcm(tests);
    testPower(tests);
    testFactorial(tests);
    testUint64Conversion(tests);
    testIntegerSquareRoot(tests);
    testIntegerCubeRoot(tests);
}

} // namespace mmcal::tests
