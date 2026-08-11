// Piなど定数の保証付き評価の回帰テスト
#include "certified_constants_tests.hpp"

#include "approximation/certified_constants.hpp"
#include "mathematics/math_ids.hpp"
#include "numeric/rational.hpp"
#include "test_framework.hpp"

#include <string>

namespace mmcal::tests {

void runCertifiedConstantTests(TestRunner& tests) {
    const auto enclosure = approximation::enclosePi(160);
    const auto piProbe = numeric::Rational::parse(
        "3.14159265358979323846264338327950288419716939937510");
    tests.expect(enclosure.interval.contains(piProbe),
        "CertifiedPi: 160-bit enclosure contains a high-precision Pi probe");
    tests.expect(enclosure.termsUsed > 0,
        "CertifiedPi: reports actual series work");

    const auto pi16 = approximation::approximatePi(16);
    tests.expectEqual(std::string{pi16.text()},
        std::string{"3.1415926535897932"},
        "CertifiedPi: rounds Pi to arbitrary requested decimal precision");

    const auto pi50 = approximation::approximatePi(50);
    tests.expectEqual(std::string{pi50.text()},
        std::string{"3.14159265358979323846264338327950288419716939937511"},
        "CertifiedPi: certifies 50 fractional digits");

    const auto pi100 = approximation::approximatePi(100);
    tests.expectEqual(std::string{pi100.text()},
        std::string{"3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170680"},
        "CertifiedPi: certifies 100 fractional digits without a fixed iteration count");

    const auto e30 = approximation::approximateConstant(mathematics::ConstantId::E, 30);
    tests.expect(e30.has_value(),
        "Certified constants: E is generated from the certified Exp[1] definition");
    tests.expect(e30 && std::string{e30->text()} == "2.718281828459045235360287471353",
        "Certified constants: E rounds correctly to 30 fractional digits");
}

} // namespace mmcal::tests
