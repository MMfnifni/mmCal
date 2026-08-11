// certified・transcendentalの回帰テスト
#include "certified_transcendental_tests.hpp"

#include "approximation/certified_exponential.hpp"
#include "approximation/certified_logarithm.hpp"
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
}

} // namespace mmcal::tests
