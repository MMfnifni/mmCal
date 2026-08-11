// 十進近似値metadataの回帰テスト
#include "decimal_approximation_tests.hpp"

#include "numeric/decimal_approximation.hpp"
#include "numeric/rational.hpp"
#include "numeric/real_number.hpp"
#include "test_framework.hpp"

#include <string>

namespace mmcal::tests {

void runDecimalApproximationTests(TestRunner& tests) {
    using numeric::BigInt;
    using numeric::DecimalApproximation;
    using numeric::Rational;
    using numeric::RealNumber;

    const auto half = DecimalApproximation::fromReal(
        RealNumber{Rational{BigInt{1}, BigInt{2}}});
    tests.expectEqual(std::string{half.text()}, std::string{"0.5"},
        "DecimalApproximation: preserves a terminating decimal");
    tests.expect(!half.isRounded(),
        "DecimalApproximation: marks a terminating decimal as exact");

    const auto oneThird = DecimalApproximation::fromReal(
        RealNumber{Rational{BigInt{1}, BigInt{3}}});
    tests.expectEqual(std::string{oneThird.text()},
        std::string{"0.3333333333333333"},
        "DecimalApproximation: uses sixteen fractional digits by default");
    tests.expect(oneThird.isRounded(),
        "DecimalApproximation: marks a repeating decimal as rounded");

    const auto twoThirds = DecimalApproximation::fromReal(
        RealNumber{Rational{BigInt{2}, BigInt{3}}});
    tests.expectEqual(std::string{twoThirds.text()},
        std::string{"0.6666666666666667"},
        "DecimalApproximation: rounds the last fractional digit");

    const auto preciseThird = DecimalApproximation::fromReal(
        RealNumber{Rational{BigInt{1}, BigInt{3}}},
        20);
    tests.expectEqual(std::string{preciseThird.text()},
        std::string{"0.33333333333333333333"},
        "DecimalApproximation: accepts an explicit fractional digit count");

    const auto negative = DecimalApproximation::fromReal(
        RealNumber{Rational{BigInt{-1}, BigInt{8}}});
    tests.expectEqual(std::string{negative.text()}, std::string{"-0.125"},
        "DecimalApproximation: preserves the sign");

    const auto roundedCarry = DecimalApproximation::fromReal(
        RealNumber{Rational{BigInt{9999}, BigInt{10001}}},
        3);
    tests.expectEqual(std::string{roundedCarry.text()}, std::string{"1.000"},
        "DecimalApproximation: carries rounding into the integer part");


    const auto fixedHalf = DecimalApproximation::fromRealFixed(
        RealNumber{Rational{BigInt{1}, BigInt{2}}},
        4);
    tests.expectEqual(std::string{fixedHalf.text()}, std::string{"0.5000"},
        "DecimalApproximation: fixed form pads exact terminating decimals");

    const auto fixedZeroDigits = DecimalApproximation::fromRealFixed(
        RealNumber{Rational{BigInt{5}, BigInt{2}}},
        0);
    tests.expectEqual(std::string{fixedZeroDigits.text()}, std::string{"2"},
        "DecimalApproximation: zero fixed digits uses nearest-even integer rounding");

    const auto certifiedStable = DecimalApproximation::fromCertifiedInterval(
        Rational{BigInt{33330}, BigInt{100000}},
        Rational{BigInt{33331}, BigInt{100000}},
        3);
    tests.expect(certifiedStable.has_value()
        && certifiedStable->text() == std::string_view{"0.333"},
        "DecimalApproximation: certifies a whole interval with one rounded result");

    const auto certifiedUnstable = DecimalApproximation::fromCertifiedInterval(
        Rational{BigInt{3334}, BigInt{10000}},
        Rational{BigInt{3336}, BigInt{10000}},
        3);
    tests.expect(!certifiedUnstable.has_value(),
        "DecimalApproximation: refuses an interval that crosses a rounding boundary");

    tests.expectThrows<std::invalid_argument>([] {
        static_cast<void>(DecimalApproximation::fromReal(
            RealNumber{Rational{BigInt{1}, BigInt{3}}},
            0));
    }, "DecimalApproximation: rejects zero precision");
}

} // namespace mmcal::tests
