// certified・transcendentalの回帰テスト
#include "certified_transcendental_tests.hpp"

#include "approximation/certified_exponential.hpp"
#include "approximation/certified_logarithm.hpp"
#include "approximation/certified_trigonometry.hpp"
#include "approximation/real_interval.hpp"
#include "numeric/big_int.hpp"
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
        "Certified Exp: reports actual Taylor work");

    const auto logTwo = approximation::encloseLogPositive(
        approximation::RealInterval::fromRational(two, 192), 192);
    const numeric::Rational logTwoProbe = numeric::Rational::parse(
        "0.69314718055994530941723212145817656807550013436025525412068000949");
    tests.expect(logTwo.interval.contains(logTwoProbe),
        "Certified Log: enclosure contains a high-precision log(2) probe");
    tests.expect(logTwo.termsUsed > 0,
        "Certified Log: reports actual atanh-series work");


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
}

} // namespace mmcal::tests
