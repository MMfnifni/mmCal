// certified・transcendentalの回帰テスト
#include "certified_transcendental_tests.hpp"

#include "approximation/certified_exponential.hpp"
#include "approximation/certified_logarithm.hpp"
#include "approximation/certified_special_functions.hpp"
#include "approximation/certified_trigonometry.hpp"
#include "approximation/real_interval.hpp"
#include "numeric/big_int.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/rational.hpp"
#include "test_framework.hpp"

namespace mmcal::tests {

void runCertifiedTranscendentalTests(TestRunner& tests) {
    const numeric::Rational one{numeric::BigInt{1}};
    const numeric::Rational two{numeric::BigInt{2}};

    const auto expOne = approximation::encloseExp(
        approximation::RealInterval::fromRational(one, 192), 192);
    const numeric::Rational eProbe = numeric::Rational::parse(
        "2.71828182845904523536028747135266249775724709369995957496696762772");
    tests.expect(expOne.interval.contains(eProbe),
        "Certified Exp: enclosure contains a high-precision E probe");
    tests.expect(expOne.termsUsed > 0,
        "Certified Exp: reports actual binary-splitting Taylor work");

    const auto w0One = approximation::encloseLambertWReal(
        approximation::RealInterval::fromRational(one, 192), 0, 160);
    tests.expect(
        w0One.lower().toRational() > numeric::Rational::parse("0.5671432904")
            && w0One.upper().toRational() < numeric::Rational::parse("0.5671432905"),
        "Certified Lambert W: principal real branch is enclosed by an independent decimal bracket");

    const numeric::Rational minusOneTenth{numeric::BigInt{-1}, numeric::BigInt{10}};
    const auto wMinusOne = approximation::encloseLambertWReal(
        approximation::RealInterval::fromRational(minusOneTenth, 192), -1, 160);
    tests.expect(
        wMinusOne.lower().toRational() > numeric::Rational::parse("-3.5771520640")
            && wMinusOne.upper().toRational() < numeric::Rational::parse("-3.5771520639"),
        "Certified Lambert W: lower real branch is enclosed by an independent decimal bracket");

    const auto logTwo = approximation::encloseLogPositive(
        approximation::RealInterval::fromRational(two, 192), 192);
    const numeric::Rational logTwoProbe = numeric::Rational::parse(
        "0.69314718055994530941723212145817656807550013436025525412068000949");
    tests.expect(logTwo.interval.contains(logTwoProbe),
        "Certified Log: enclosure contains a high-precision log(2) probe");
    tests.expect(logTwo.termsUsed > 0,
        "Certified Log: reports actual binary-splitting atanh-series work");

    const auto expHigh = approximation::encloseExp(
        approximation::RealInterval::fromRational(one, 400), 400);
    const auto expHighDecimal = numeric::DecimalApproximation::fromCertifiedInterval(
        expHigh.interval.lower().toRational(), expHigh.interval.upper().toRational(), 100);
    tests.expect(expHighDecimal.has_value()
            && expHighDecimal->text()
                == "2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274",
        "Certified Exp: binary-splitting series certifies 100 fractional digits");

    const auto logHigh = approximation::encloseLogPositive(
        approximation::RealInterval::fromRational(two, 400), 400);
    const auto logHighDecimal = numeric::DecimalApproximation::fromCertifiedInterval(
        logHigh.interval.lower().toRational(), logHigh.interval.upper().toRational(), 100);
    tests.expect(logHighDecimal.has_value()
            && logHighDecimal->text()
                == "0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875",
        "Certified Log: binary-splitting series certifies 100 fractional digits");


    // 極端な大きさでも、単発のdecimal文字列ではなく包含不変量を検証する。
    const numeric::Rational hugePositive{numeric::BigInt{600}};
    const numeric::Rational hugeNegative{numeric::BigInt{-600}};
    const auto expPositive = approximation::encloseExp(
        approximation::RealInterval::fromRational(hugePositive, 320), 320);
    const auto expNegative = approximation::encloseExp(
        approximation::RealInterval::fromRational(hugeNegative, 320), 320);
    const auto expProduct = approximation::multiply(
        expPositive.interval, expNegative.interval, 300);
    tests.expect(expProduct.contains(one),
        "Certified Exp extreme: exp(600)*exp(-600) enclosure contains 1");
    tests.expect(expPositive.squarings > 0 && expNegative.squarings > 0,
        "Certified Exp extreme: large arguments use range reduction without overflow");

    const numeric::BigInt powerOfTwo = numeric::BigInt{1} << 768;
    const numeric::Rational largePower{powerOfTwo};
    const numeric::Rational reciprocalPower{numeric::BigInt{1}, powerOfTwo};
    const auto logLarge = approximation::encloseLogPositive(
        approximation::RealInterval::fromRational(largePower, 320), 320);
    const auto logSmall = approximation::encloseLogPositive(
        approximation::RealInterval::fromRational(reciprocalPower, 320), 320);
    const auto logSymmetry = approximation::add(logLarge.interval, logSmall.interval, 300);
    tests.expect(logSymmetry.contains(numeric::Rational{numeric::BigInt{0}}),
        "Certified Log extreme: log(2^768)+log(2^-768) enclosure contains 0");

    const numeric::BigInt hugeTurnsBase = (numeric::BigInt{1} << 1024) + numeric::BigInt{123};
    const numeric::Rational hugeTurns{
        hugeTurnsBase * numeric::BigInt{6} + numeric::BigInt{1}, numeric::BigInt{6}};
    const numeric::Rational oneSixth{numeric::BigInt{1}, numeric::BigInt{6}};
    tests.expectEqual(
        approximation::approximateSinTurns(hugeTurns, 50).value.text(),
        approximation::approximateSinTurns(oneSixth, 50).value.text(),
        "Certified Trig extreme: huge whole-turn offsets reduce exactly before approximation");
    tests.expectEqual(
        approximation::approximateCosTurns(hugeTurns, 50).value.text(),
        approximation::approximateCosTurns(oneSixth, 50).value.text(),
        "Certified Trig extreme: huge angle reduction preserves cosine rounding");


    const auto gammaOne = approximation::encloseGammaReal(
        approximation::RealInterval::fromRational(one, 192), 192);
    tests.expect(gammaOne.isPoint() && gammaOne.lower().toRational() == one,
        "Certified Gamma: Gamma(1) uses the exact fast path");

    const numeric::Rational oneThird{numeric::BigInt{1}, numeric::BigInt{3}};
    const auto gammaOneThird = approximation::encloseGammaReal(
        approximation::RealInterval::fromRational(oneThird, 384), 320);
    const auto gammaOneThirdDecimal = numeric::DecimalApproximation::fromCertifiedInterval(
        gammaOneThird.lower().toRational(), gammaOneThird.upper().toRational(), 49);
    tests.expect(gammaOneThirdDecimal.has_value()
            && gammaOneThirdDecimal->text()
                == "2.6789385347077476336556929409746776441286893779573",
        "Certified Gamma: interval Stirling path keeps a high-precision enclosure");

    const auto gammaOneThirdExact = approximation::encloseGammaRational(oneThird, 1024);
    const auto gammaOneThirdExactDecimal = numeric::DecimalApproximation::fromCertifiedInterval(
        gammaOneThirdExact.lower().toRational(), gammaOneThirdExact.upper().toRational(), 49);
    tests.expect(gammaOneThirdExactDecimal.has_value()
            && gammaOneThirdExactDecimal->text()
                == "2.6789385347077476336556929409746776441286893779573",
        "Certified Gamma: exact Rational dispatch and high-precision Stirling planner preserve the enclosure");

    const numeric::Rational twoThirds{numeric::BigInt{2}, numeric::BigInt{3}};
    const numeric::Rational oneQuarter{numeric::BigInt{1}, numeric::BigInt{4}};
    const auto ibetaPoint = approximation::encloseIncompleteBetaRegularized(
        oneThird, twoThirds,
        approximation::RealInterval::fromRational(oneQuarter, 224), 160);
    const auto ibetaDecimal = numeric::DecimalApproximation::fromCertifiedInterval(
        ibetaPoint.lower().toRational(), ibetaPoint.upper().toRational(), 20);
    tests.expect(ibetaDecimal.has_value() && ibetaDecimal->text() == "0.53302858123542523627",
        "Certified ibeta: shared point normalization preserves the certified enclosure");
}

} // namespace mmcal::tests
