// 整数・有理数を統合するRealNumberの回帰テスト
#include "real_number_tests.hpp"

#include "numeric/real_number.hpp"
#include "test_framework.hpp"

#include <stdexcept>

namespace mmcal::tests {

void runRealNumberTests(TestRunner& tests) {
    using numeric::BigInt;
    using numeric::Rational;
    using numeric::RealNumber;

    const RealNumber integer{BigInt::parse("123456789012345678901234567890")};
    tests.expect(integer.isInteger(), "RealNumber: stores BigInt as integer");
    tests.expectEqual(integer.toString(), "123456789012345678901234567890", "RealNumber: formats a large integer");

    const RealNumber reduced{Rational{BigInt{6}, BigInt{3}}};
    tests.expect(reduced.isInteger(), "RealNumber: reduces denominator-one rational to integer");
    tests.expectEqual(reduced.toString(), "2", "RealNumber: formats reduced integer");

    const RealNumber fraction{Rational{BigInt{7}, BigInt{3}}};
    tests.expect(fraction.isRational(), "RealNumber: stores non-Integer rational");
    tests.expectEqual(fraction.toString(), "7/3", "RealNumber: formats rational");

    tests.expectEqual(
        (RealNumber{BigInt{2}} + fraction).toString(),
        "13/3",
        "RealNumber: adds integer and rational");
    tests.expectEqual(
        (fraction - RealNumber{BigInt{3}}).toString(),
        "-2/3",
        "RealNumber: subtracts integer and rational");
    tests.expectEqual(
        (fraction * RealNumber{Rational{BigInt{9}, BigInt{14}}}).toString(),
        "3/2",
        "RealNumber: multiplies rationals");
    tests.expectEqual(
        (RealNumber{BigInt{7}} / RealNumber{BigInt{2}}).toString(),
        "7/2",
        "RealNumber: keeps integer division exact");
    tests.expectEqual(
        (RealNumber{BigInt{8}} / RealNumber{BigInt{2}}).toString(),
        "4",
        "RealNumber: reduces exact division to integer");

    tests.expect(
        RealNumber{BigInt{1}} == RealNumber{Rational{BigInt{2}, BigInt{2}}},
        "RealNumber: compares equivalent internal representations");
    tests.expect(
        RealNumber{Rational{BigInt{-1}, BigInt{2}}} < RealNumber{BigInt{0}},
        "RealNumber: orders integer and rational");
    tests.expectEqual(
        RealNumber{Rational{BigInt{-3}, BigInt{4}}}.abs().toString(),
        "3/4",
        "RealNumber: absolute value");

    tests.expectThrows<std::logic_error>(
        [&] { static_cast<void>(fraction.asInteger()); },
        "RealNumber: rejects rational integer access");
    tests.expectThrows<std::domain_error>(
        [] { static_cast<void>(RealNumber{BigInt{1}} / RealNumber{}); },
        "RealNumber: rejects division by zero");
}

} // namespace mmcal::tests
