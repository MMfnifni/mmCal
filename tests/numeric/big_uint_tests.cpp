// BigInt内部の符号なし多倍長整数の回帰テスト
#include "numeric/big_uint_tests.hpp"

#include "numeric/detail/big_uint.hpp"
#include "test_framework.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace mmcal::tests {
namespace {

using mmcal::numeric::detail::BigUInt;
using mmcal::numeric::detail::BigUIntDivModResult;

BigUIntDivModResult binaryDivmodReference(
    const BigUInt& dividend,
    const BigUInt& divisor) {
    if (divisor.isZero())
        throw std::domain_error("BigUInt division by zero");

    if (dividend < divisor)
        return {BigUInt{}, dividend};

    BigUInt quotient;
    BigUInt remainder = dividend;
    std::size_t shift = remainder.bitLength() - divisor.bitLength();
    BigUInt shiftedDivisor = divisor << shift;
    BigUInt quotientBit = BigUInt{1} << shift;

    for (;;) {
        if (remainder >= shiftedDivisor) {
            remainder -= shiftedDivisor;
            quotient += quotientBit;
        }

        if (shift == 0)
            break;

        shiftedDivisor >>= 1;
        quotientBit >>= 1;
        --shift;
    }

    return {std::move(quotient), std::move(remainder)};
}

std::string randomHex(
    std::mt19937_64& generator,
    std::size_t limbCount) {
    constexpr std::string_view digits = "0123456789ABCDEF";
    std::string text(limbCount * 8, '0');

    for (char& digit : text)
        digit = digits[generator() & 0xF];

    text.front() = digits[1 + generator() % 15];
    return text;
}

void testConstruction(TestRunner& tests) {
    const BigUInt zero;
    tests.expect(zero.isZero(), "default value is zero");
    tests.expectEqual(zero.limbCount(), std::size_t{0}, "zero has no limbs");
    tests.expectEqual(zero.bitLength(), std::size_t{0}, "zero bit length");

    const BigUInt one{1};
    tests.expect(!one.isZero(), "one is nonzero");
    tests.expectEqual(one.bitLength(), std::size_t{1}, "one bit length");

    const BigUInt max64{std::numeric_limits<std::uint64_t>::max()};
    tests.expectEqual(max64.limbCount(), std::size_t{2}, "uint64 max uses two limbs");
    tests.expectEqual(max64.bitLength(), std::size_t{64}, "uint64 max bit length");
}

void testDecimalConversion(TestRunner& tests) {
    const auto max64 = BigUInt::parse("18446744073709551615");
    tests.expectEqual(
        max64.toString(),
        std::string{"18446744073709551615"},
        "parse uint64 max as decimal");

    const auto beyond64 = BigUInt::parse("18446744073709551616");
    tests.expectEqual(
        beyond64.toString(),
        std::string{"18446744073709551616"},
        "parse value beyond uint64");

    const std::string huge = "1234567890123456789012345678901234567890";
    tests.expectEqual(
        BigUInt::parse(huge).toString(),
        huge,
        "large decimal round trip");

    // 10^9 chunk境界をまたぐ長さでも先頭chunkと0埋めを崩さないことを確認する。
    std::string chunked = "7";
    for (std::size_t i = 0; i < 120; ++i)
        chunked += "123456789";
    tests.expectEqual(
        BigUInt::parse(chunked).toString(),
        chunked,
        "chunked decimal conversion round trip");

    // 128 limbsを超えてdivide-and-conquer 10進変換へ入る長さでも、先頭0やblock幅を作らない。
    std::string divideAndConquer = "9";
    for (std::size_t i = 0; i < 360; ++i)
        divideAndConquer += "314159265";
    tests.expectEqual(
        BigUInt::parse(divideAndConquer).toString(),
        divideAndConquer,
        "divide-and-conquer decimal conversion round trip");
}

void testRadixConversion(TestRunner& tests) {
    const auto value = BigUInt::parse("FFFFFFFFFFFFFFFFFFFFFFFF", 16);
    tests.expectEqual(
        value.toString(16),
        std::string{"FFFFFFFFFFFFFFFFFFFFFFFF"},
        "hexadecimal round trip");

    tests.expectEqual(
        BigUInt::parse("101010", 2).toString(),
        std::string{"42"},
        "binary to decimal");

    tests.expectEqual(
        BigUInt::parse("Z", 36).toString(),
        std::string{"35"},
        "base 36 digit");

    const auto decimal = BigUInt::parse("123456789012345678901234567890");
    for (unsigned radix = 2; radix <= 36; ++radix) {
        const auto encoded = decimal.toString(radix);
        const auto decoded = BigUInt::parse(encoded, radix);
        tests.expect(decoded == decimal, "radix round trip");
    }
}

void testComparison(TestRunner& tests) {
    const auto a = BigUInt::parse("4294967295");
    const auto b = BigUInt::parse("4294967296");
    const auto c = BigUInt::parse("4294967296");

    tests.expect(a < b, "comparison across limb boundary");
    tests.expect(b == c, "equal parsed values");
    tests.expect(c > a, "greater comparison");
}

void testAddition(TestRunner& tests) {
    const BigUInt zero;
    const BigUInt one{1};

    tests.expectEqual((zero + zero).toString(), std::string{"0"}, "zero plus zero");
    tests.expectEqual((zero + one).toString(), std::string{"1"}, "zero plus one");
    tests.expectEqual((one + zero).toString(), std::string{"1"}, "one plus zero");

    const auto limbMax = BigUInt::parse("FFFFFFFF", 16);
    tests.expectEqual(
        (limbMax + one).toString(16),
        std::string{"100000000"},
        "addition carries into a new limb");

    const auto threeLimbMax = BigUInt::parse("FFFFFFFFFFFFFFFFFFFFFFFF", 16);
    tests.expectEqual(
        (threeLimbMax + one).toString(16),
        std::string{"1000000000000000000000000"},
        "addition propagates carry across multiple limbs");

    const auto a = BigUInt::parse("123456789012345678901234567890");
    const auto b = BigUInt::parse("987654321098765432109876543210");
    tests.expectEqual(
        (a + b).toString(),
        std::string{"1111111110111111111011111111100"},
        "large arbitrary addition");

    auto self = BigUInt::parse("18446744073709551616");
    self += self;
    tests.expectEqual(
        self.toString(),
        std::string{"36893488147419103232"},
        "self addition is safe");

    const auto sum = a + b;
    tests.expect(sum == b + a, "addition is commutative");
}

void testSubtraction(TestRunner& tests) {
    const BigUInt zero;
    const BigUInt one{1};

    tests.expectEqual((one - zero).toString(), std::string{"1"}, "subtract zero");
    tests.expectEqual((one - one).toString(), std::string{"0"}, "subtract equal values");

    const auto oneLimbBoundary = BigUInt::parse("100000000", 16);
    tests.expectEqual(
        (oneLimbBoundary - one).toString(16),
        std::string{"FFFFFFFF"},
        "subtraction borrows across a limb boundary");

    const auto threeLimbBoundary = BigUInt::parse("1000000000000000000000000", 16);
    tests.expectEqual(
        (threeLimbBoundary - one).toString(16),
        std::string{"FFFFFFFFFFFFFFFFFFFFFFFF"},
        "subtraction propagates borrow across multiple limbs");

    const auto larger = BigUInt::parse("987654321098765432109876543210");
    const auto smaller = BigUInt::parse("123456789012345678901234567890");
    tests.expectEqual(
        (larger - smaller).toString(),
        std::string{"864197532086419753208641975320"},
        "large arbitrary subtraction");

    auto self = BigUInt::parse("12345678901234567890");
    self -= self;
    tests.expect(self.isZero(), "self subtraction produces normalized zero");
    tests.expectEqual(self.limbCount(), std::size_t{0}, "self subtraction removes zero limbs");

    const auto sum = larger + smaller;
    tests.expect(sum - smaller == larger, "subtraction reverses addition");

    tests.expectThrows<std::underflow_error>(
        [] {
            BigUInt value{1};
            value -= BigUInt{2};
        },
        "negative unsigned subtraction is rejected");
}


void testMultiplication(TestRunner& tests) {
    const BigUInt zero;
    const BigUInt one{1};
    const BigUInt two{2};

    tests.expectEqual((zero * zero).toString(), std::string{"0"}, "zero times zero");
    tests.expectEqual((zero * two).toString(), std::string{"0"}, "zero times nonzero");
    tests.expectEqual((two * zero).toString(), std::string{"0"}, "nonzero times zero");
    tests.expectEqual((one * two).toString(), std::string{"2"}, "one is multiplicative identity");

    const auto limbMax = BigUInt::parse("FFFFFFFF", 16);
    tests.expectEqual(
        (limbMax * limbMax).toString(16),
        std::string{"FFFFFFFE00000001"},
        "maximum limbs multiply without overflowing the accumulator");

    const auto threeLimbMax = BigUInt::parse("FFFFFFFFFFFFFFFFFFFFFFFF", 16);
    tests.expectEqual(
        (threeLimbMax * threeLimbMax).toString(16),
        std::string{"FFFFFFFFFFFFFFFFFFFFFFFE000000000000000000000001"},
        "multiplication carries across multiple limbs");

    const auto a = BigUInt::parse("1234567890123456789012345678901234567890");
    const auto b = BigUInt::parse("9876543210987654321098765432109876543210");
    tests.expectEqual(
        (a * b).toString(),
        std::string{
            "12193263113702179522618503273386678859448712086533622923332237463801111263526900"},
        "large arbitrary multiplication");

    tests.expect(a * b == b * a, "multiplication is commutative");

    auto self = BigUInt::parse("18446744073709551616");
    self *= self;
    tests.expectEqual(
        self.toString(),
        std::string{"340282366920938463463374607431768211456"},
        "self multiplication is safe");

    const auto x = BigUInt::parse("12345678901234567890");
    const auto y = BigUInt::parse("11111111111111111111");
    const auto z = BigUInt::parse("22222222222222222222");
    tests.expect(
        x * (y + z) == x * y + x * z,
        "multiplication distributes over addition");
}

void testBitShifts(TestRunner& tests) {
    const BigUInt zero;
    const BigUInt one{1};

    tests.expectEqual((zero << 1000).toString(), std::string{"0"}, "left shift keeps zero normalized");
    tests.expectEqual((zero >> 1000).toString(), std::string{"0"}, "right shift keeps zero normalized");
    tests.expectEqual((one << 0).toString(), std::string{"1"}, "zero-bit left shift is identity");
    tests.expectEqual((one >> 0).toString(), std::string{"1"}, "zero-bit right shift is identity");

    tests.expectEqual((one << 1).toString(), std::string{"2"}, "left shift by one bit");
    tests.expectEqual((one << 31).toString(16), std::string{"80000000"}, "left shift within one limb");
    tests.expectEqual((one << 32).toString(16), std::string{"100000000"}, "left shift by one whole limb");
    tests.expectEqual((one << 65).toString(16), std::string{"20000000000000000"}, "left shift by limbs and bits");

    const auto limbMax = BigUInt::parse("FFFFFFFF", 16);
    tests.expectEqual(
        (limbMax << 1).toString(16),
        std::string{"1FFFFFFFE"},
        "left shift carries into a new limb");

    const auto mixed = BigUInt::parse("123456789ABCDEF0123456789ABCDEF", 16);
    tests.expectEqual(
        (mixed << 36).toString(16),
        std::string{"123456789ABCDEF0123456789ABCDEF000000000"},
        "large value left shift");

    tests.expectEqual(
        (BigUInt::parse("100000000", 16) >> 32).toString(),
        std::string{"1"},
        "right shift by one whole limb");
    tests.expectEqual(
        (BigUInt::parse("1FFFFFFFE", 16) >> 1).toString(16),
        std::string{"FFFFFFFF"},
        "right shift combines adjacent limbs");
    tests.expectEqual(
        (BigUInt::parse("123456789ABCDEF000000000", 16) >> 36).toString(16),
        std::string{"123456789ABCDEF"},
        "right shift by limbs and bits");

    auto mutableValue = BigUInt::parse("FEDCBA9876543210", 16);
    mutableValue <<= 37;
    mutableValue >>= 37;
    tests.expectEqual(
        mutableValue.toString(16),
        std::string{"FEDCBA9876543210"},
        "matching shifts restore a value when no low bits are discarded");

    const auto arbitrary = BigUInt::parse("123456789012345678901234567890");
    constexpr std::array<std::size_t, 8> shifts{1, 7, 31, 32, 33, 63, 64, 95};
    for (const std::size_t shift : shifts) {
        tests.expect(
            ((arbitrary << shift) >> shift) == arbitrary,
            "left then right shift round trip");
    }

    auto discarded = BigUInt::parse("FFFFFFFFFFFFFFFF", 16);
    discarded >>= 1000;
    tests.expect(discarded.isZero(), "oversized right shift produces zero");
    tests.expectEqual(discarded.limbCount(), std::size_t{0}, "oversized shift keeps canonical zero");
}


void testDivision(TestRunner& tests) {
    const BigUInt zero;
    const BigUInt one{1};
    const BigUInt two{2};

    const auto zeroByOne = divmod(zero, one);
    tests.expect(zeroByOne.quotient.isZero(), "zero divided by one has zero quotient");
    tests.expect(zeroByOne.remainder.isZero(), "zero divided by one has zero remainder");

    const auto oneByOne = divmod(one, one);
    tests.expectEqual(oneByOne.quotient.toString(), std::string{"1"}, "one divided by one");
    tests.expect(oneByOne.remainder.isZero(), "one divided by one is exact");

    const auto smaller = divmod(BigUInt{7}, BigUInt{11});
    tests.expect(smaller.quotient.isZero(), "smaller dividend produces zero quotient");
    tests.expectEqual(smaller.remainder.toString(), std::string{"7"}, "smaller dividend becomes remainder");

    const auto simple = divmod(BigUInt{100}, BigUInt{7});
    tests.expectEqual(simple.quotient.toString(), std::string{"14"}, "simple quotient");
    tests.expectEqual(simple.remainder.toString(), std::string{"2"}, "simple remainder");

    const auto limbExact = divmod(
        BigUInt::parse("FFFFFFFFFFFFFFFF", 16),
        BigUInt::parse("FFFFFFFF", 16));
    tests.expectEqual(
        limbExact.quotient.toString(16),
        std::string{"100000001"},
        "division across limb boundary");
    tests.expect(limbExact.remainder.isZero(), "limb-boundary division is exact");

    const auto dividend = BigUInt::parse("1234567890123456789012345678901234567890");
    const auto divisor = BigUInt::parse("98765432109876543210");
    const auto large = divmod(dividend, divisor);
    tests.expectEqual(
        large.quotient.toString(),
        std::string{"12499999886093750001"},
        "large arbitrary quotient");
    tests.expectEqual(
        large.remainder.toString(),
        std::string{"54205246805420524680"},
        "large arbitrary remainder");

    tests.expect(
        large.quotient * divisor + large.remainder == dividend,
        "division reconstructs the dividend");
    tests.expect(large.remainder < divisor, "remainder is smaller than divisor");

    const auto powerDividend = one << 257;
    const auto powerDivisor = one << 129;
    const auto powers = divmod(powerDividend, powerDivisor);
    tests.expect(powers.quotient == (one << 128), "power-of-two quotient");
    tests.expect(powers.remainder.isZero(), "power-of-two division is exact");

    tests.expectEqual((BigUInt{100} / BigUInt{7}).toString(), std::string{"14"}, "operator slash");
    tests.expectEqual((BigUInt{100} % BigUInt{7}).toString(), std::string{"2"}, "operator percent");

    auto divideAssign = dividend;
    divideAssign /= divisor;
    tests.expect(divideAssign == large.quotient, "division assignment");

    auto remainderAssign = dividend;
    remainderAssign %= divisor;
    tests.expect(remainderAssign == large.remainder, "remainder assignment");

    auto selfDivision = BigUInt::parse("999999999999999999999999999999");
    selfDivision /= selfDivision;
    tests.expectEqual(selfDivision.toString(), std::string{"1"}, "self division is safe");

    auto selfRemainder = BigUInt::parse("999999999999999999999999999999");
    selfRemainder %= selfRemainder;
    tests.expect(selfRemainder.isZero(), "self remainder is zero");

    constexpr std::array<std::pair<std::uint64_t, std::uint64_t>, 8> cases{{
        {1, 1},
        {10, 3},
        {std::numeric_limits<std::uint32_t>::max(), 97},
        {std::uint64_t{1} << 32, 65537},
        {std::numeric_limits<std::uint64_t>::max(), 2},
        {std::numeric_limits<std::uint64_t>::max(), 3},
        {std::numeric_limits<std::uint64_t>::max(), std::numeric_limits<std::uint32_t>::max()},
        {9876543210123456789ULL, 123456789ULL}
    }};

    for (const auto& [nativeDividend, nativeDivisor] : cases) {
        const auto result = divmod(BigUInt{nativeDividend}, BigUInt{nativeDivisor});
        tests.expectEqual(
            result.quotient.toString(),
            std::to_string(nativeDividend / nativeDivisor),
            "native quotient cross-check");
        tests.expectEqual(
            result.remainder.toString(),
            std::to_string(nativeDividend % nativeDivisor),
            "native remainder cross-check");
    }

    tests.expectThrows<std::domain_error>(
        [] { static_cast<void>(divmod(BigUInt{1}, BigUInt{})); },
        "divmod rejects division by zero");
    tests.expectThrows<std::domain_error>(
        [] {
            auto value = BigUInt{1};
            value /= BigUInt{};
        },
        "operator slash assignment rejects division by zero");
    tests.expectThrows<std::domain_error>(
        [] {
            auto value = BigUInt{1};
            value %= BigUInt{};
        },
        "operator percent assignment rejects division by zero");

    bool normalizationCasesPass = true;
    for (unsigned topBit = 0; topBit < 32; ++topBit) {
        const BigUInt generatedDivisor =
            (BigUInt{1} << (32 + topBit)) + BigUInt{0xA5A5A5A5};
        const BigUInt generatedQuotient =
            BigUInt::parse("123456789ABCDEF012345", 16);
        const BigUInt generatedRemainder =
            BigUInt{0x13579BDF} % generatedDivisor;
        const BigUInt generatedDividend =
            generatedDivisor * generatedQuotient + generatedRemainder;
        const auto generated = divmod(generatedDividend, generatedDivisor);

        normalizationCasesPass = normalizationCasesPass
            && generated.quotient == generatedQuotient
            && generated.remainder == generatedRemainder;
    }
    tests.expect(
        normalizationCasesPass,
        "normalized division handles every top-limb bit position");

    std::mt19937_64 generator{0x5A17D1A1D3ULL};
    bool referenceCasesPass = true;

    for (std::size_t caseIndex = 0; caseIndex < 128; ++caseIndex) {
        const std::size_t dividendLimbs = 1 + generator() % 12;
        const std::size_t divisorLimbs = 2 + generator() % 8;
        const BigUInt generatedDividend =
            BigUInt::parse(randomHex(generator, dividendLimbs), 16);
        const BigUInt generatedDivisor =
            BigUInt::parse(randomHex(generator, divisorLimbs), 16);

        const auto actual = divmod(generatedDividend, generatedDivisor);
        const auto reference =
            binaryDivmodReference(generatedDividend, generatedDivisor);

        referenceCasesPass = referenceCasesPass
            && actual.quotient == reference.quotient
            && actual.remainder == reference.remainder
            && actual.quotient * generatedDivisor + actual.remainder
                == generatedDividend
            && actual.remainder < generatedDivisor;

        if (!referenceCasesPass)
            break;
    }

    tests.expect(
        referenceCasesPass,
        "normalized division matches the binary reference implementation");
}

void testErrors(TestRunner& tests) {
    tests.expectThrows<std::invalid_argument>(
        [] { static_cast<void>(BigUInt::parse("2", 2)); },
        "invalid digit is rejected");

    tests.expectThrows<std::invalid_argument>(
        [] { static_cast<void>(BigUInt::parse("-1")); },
        "negative value is rejected");
}

} // namespace

void runBigUIntTests(TestRunner& tests) {
    testConstruction(tests);
    testDecimalConversion(tests);
    testRadixConversion(tests);
    testComparison(tests);
    testAddition(tests);
    testSubtraction(tests);
    testMultiplication(tests);
    testBitShifts(tests);
    testDivision(tests);
    testErrors(tests);
}

} // namespace mmcal::tests
