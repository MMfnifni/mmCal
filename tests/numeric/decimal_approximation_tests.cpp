// 十進近似値metadataの回帰テスト
#include "decimal_approximation_tests.hpp"

#include "numeric/approximation_quality.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/integer_algorithms.hpp"
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

    const auto certifiedTrailingZeros = DecimalApproximation::fromCertifiedInterval(
        Rational{BigInt{149999999}, BigInt{100000000}},
        Rational{BigInt{150000001}, BigInt{100000000}},
        6);
    tests.expect(certifiedTrailingZeros.has_value()
        && certifiedTrailingZeros->text() == std::string_view{"1.50"}
        && certifiedTrailingZeros->fractionalDigits() == 2
        && certifiedTrailingZeros->requestedFractionalDigits() == 6,
        "DecimalApproximation: certified display compacts redundant trailing zeros but retains one lower-order zero");

    const auto certifiedSignificantTrailingZero = DecimalApproximation::fromCertifiedInterval(
        Rational{BigInt{122999999}, BigInt{100000000}},
        Rational{BigInt{123000001}, BigInt{100000000}},
        6);
    tests.expect(certifiedSignificantTrailingZero.has_value()
        && certifiedSignificantTrailingZero->text() == std::string_view{"1.230"}
        && certifiedSignificantTrailingZero->fractionalDigits() == 3
        && certifiedSignificantTrailingZero->requestedFractionalDigits() == 6,
        "DecimalApproximation: certified display retains one zero below the last observed nonzero digit");

    const auto certifiedInteger = DecimalApproximation::fromCertifiedInterval(
        Rational{BigInt{99999999}, BigInt{100000000}},
        Rational{BigInt{100000001}, BigInt{100000000}},
        6);
    tests.expect(certifiedInteger.has_value()
        && certifiedInteger->text() == std::string_view{"1.0"}
        && certifiedInteger->requestedFractionalDigits() == 6,
        "DecimalApproximation: certified integer-like output keeps one fractional zero");

    const auto certifiedUnstable = DecimalApproximation::fromCertifiedInterval(
        Rational{BigInt{3334}, BigInt{10000}},
        Rational{BigInt{3336}, BigInt{10000}},
        3);
    tests.expect(!certifiedUnstable.has_value(),
        "DecimalApproximation: refuses an interval that crosses a rounding boundary");

    const auto withInformation = DecimalApproximation::fromCertifiedIntervalWithInformation(
        Rational{BigInt{33330}, BigInt{100000}},
        Rational{BigInt{33331}, BigInt{100000}},
        Rational{BigInt{3332}, BigInt{10000}},
        Rational{BigInt{3334}, BigInt{10000}},
        3);
    tests.expect(withInformation.has_value()
        && withInformation->certifiedLower() == Rational{BigInt{33330}, BigInt{100000}}
        && withInformation->certifiedUpper() == Rational{BigInt{33331}, BigInt{100000}}
        && withInformation->informationLower() == Rational{BigInt{665}, BigInt{2000}}
        && withInformation->informationUpper() == Rational{BigInt{667}, BigInt{2000}},
        "DecimalApproximation: preserves propagated information and the output rounding quantum");

    const auto exactSignificantTerminating = DecimalApproximation::fromRealSignificant(
        RealNumber{Rational{BigInt{617}, BigInt{500}}},
        20);
    tests.expectEqual(std::string{exactSignificantTerminating.text()}, std::string{"1.2340"},
        "DecimalApproximation: significant approximation keeps one provenance zero for a terminating value");

    const auto significantTrailingZeros = DecimalApproximation::fromCertifiedIntervalSignificant(
        Rational{BigInt::parse("1233999999999999999999"), numeric::pow(BigInt{10}, 21)},
        Rational{BigInt::parse("1234000000000000000001"), numeric::pow(BigInt{10}, 21)},
        20);
    tests.expect(significantTrailingZeros.has_value()
        && significantTrailingZeros->text() == std::string_view{"1.2340"},
        "DecimalApproximation: non-point significant display retains one trailing zero");

    const auto zeroCentered = DecimalApproximation::fromCertifiedIntervalSignificant(
        Rational{BigInt{-1}, numeric::pow(BigInt{10}, 40)},
        Rational{BigInt{1}, numeric::pow(BigInt{10}, 40)},
        20);
    tests.expect(zeroCentered.has_value()
        && zeroCentered->text() == std::string_view{"0.0"}
        && zeroCentered->fractionalDigits() == 1
        && zeroCentered->requestedFractionalDigits() > 20,
        "DecimalApproximation: zero-centered display stays compact while metadata preserves accuracy");

    tests.expect(preciseThird.informationLower() <= preciseThird.certifiedLower()
        && preciseThird.certifiedUpper() <= preciseThird.informationUpper(),
        "DecimalApproximation: information enclosure always contains certified truth enclosure");


    const auto defaultPointInformation = DecimalApproximation::fromCertifiedIntervalWithInformationSignificant(
        Rational{}, Rational{}, Rational{}, Rational{}, 5);
    tests.expect(defaultPointInformation.has_value()
        && defaultPointInformation->certifiedExactlyZero()
        && !defaultPointInformation->informationExactlyZero(),
        "DecimalApproximation: default explicit-information construction keeps display quantization around a certified zero point");

    const auto preservedPointInformation = DecimalApproximation::fromCertifiedIntervalWithInformationSignificant(
        Rational{}, Rational{}, Rational{}, Rational{}, 5,
        numeric::InformationQuantization::PreserveExactPoint);
    tests.expect(preservedPointInformation.has_value()
        && preservedPointInformation->certifiedExactlyZero()
        && preservedPointInformation->informationExactlyZero()
        && preservedPointInformation->informationLower().isZero()
        && preservedPointInformation->informationUpper().isZero(),
        "DecimalApproximation: PreserveExactPoint keeps proven point information without reintroducing display quantum");

    const auto finiteZero = DecimalApproximation::fromRealSignificant(RealNumber{Rational{}}, 5);
    tests.expect(finiteZero.certifiedExactlyZero()
        && !finiteZero.informationExactlyZero(),
        "DecimalApproximation: N-style finite precision zero distinguishes certified truth from reusable information");

    const auto exactImaginary = DecimalApproximation::fromRealSignificant(RealNumber{Rational{BigInt{1}}}, 20);
    const auto pureImaginary = numeric::ComplexDecimalApproximation::fromComponents(
        *defaultPointInformation, exactImaginary, true, false);
    tests.expect(pureImaginary.realCertifiedExactlyZero()
        && pureImaginary.realInformationExactlyZero()
        && numeric::precisionDigits(pureImaginary) >= 19,
        "ComplexDecimalApproximation: information-exact zero component does not cap whole-complex precision");

    tests.expectThrows<std::invalid_argument>([] {
        static_cast<void>(DecimalApproximation::fromCertifiedIntervalWithInformation(
            Rational{BigInt{1}, BigInt{3}},
            Rational{BigInt{1}, BigInt{3}},
            Rational{BigInt{34}, BigInt{100}},
            Rational{BigInt{35}, BigInt{100}},
            3));
    }, "DecimalApproximation: rejects information enclosure that excludes certified truth");

    tests.expectThrows<std::invalid_argument>([] {
        static_cast<void>(DecimalApproximation::fromReal(
            RealNumber{Rational{BigInt{1}, BigInt{3}}},
            0));
    }, "DecimalApproximation: rejects zero precision");
}

} // namespace mmcal::tests
