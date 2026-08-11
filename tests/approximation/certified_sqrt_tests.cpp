// 実平方根の保証付き評価の回帰テスト
#include "certified_sqrt_tests.hpp"

#include "approximation/certified_sqrt.hpp"
#include "approximation/certified_complex_sqrt.hpp"
#include "approximation/complex_interval.hpp"
#include "approximation/real_interval.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "numeric/real_number.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>

namespace mmcal::tests {
namespace {

using numeric::BigInt;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

} // namespace

void runCertifiedSqrtTests(TestRunner& tests) {
    const auto sqrt2 = approximation::encloseSqrt(rational(2), 128);
    const Rational lower = sqrt2.interval.lower().toRational();
    const Rational upper = sqrt2.interval.upper().toRational();
    tests.expect(lower * lower <= rational(2) && upper * upper >= rational(2),
        "CertifiedSqrt: interval squares bracket exact value");
    tests.expect(lower >= rational(1) && upper <= rational(2),
        "CertifiedSqrt: sqrt(2) enclosure stays in obvious mathematical bounds");

    const auto exact = approximation::encloseSqrt(rational(9, 16), 64);
    tests.expect(exact.interval.contains(rational(3, 4)),
        "CertifiedSqrt: exact dyadic square root is enclosed");

    tests.expectEqual(
        std::string{approximation::approximateSqrt(RealNumber{rational(2)}, 50).text()},
        std::string{"1.41421356237309504880168872420969807856967187537695"},
        "CertifiedSqrt: sqrt(2) rounds correctly to 50 fractional digits");

    tests.expectThrows<std::domain_error>([] {
        static_cast<void>(approximation::encloseSqrt(rational(-2), 64));
    }, "CertifiedSqrt: real enclosure rejects negative input");

    const approximation::ComplexInterval threePlusFourI{
        approximation::RealInterval::fromRational(rational(3), 96),
        approximation::RealInterval::fromRational(rational(4), 96)
    };
    const auto complexRoot = approximation::enclosePrincipalComplexSqrt(
        threePlusFourI, 96);
    tests.expect(complexRoot.real().contains(rational(2))
        && complexRoot.imaginary().contains(rational(1)),
        "CertifiedComplexSqrt: principal sqrt of 3+4I encloses 2+I");

    const approximation::ComplexInterval minusThreeMinusFourI{
        approximation::RealInterval::fromRational(rational(-3), 96),
        approximation::RealInterval::fromRational(rational(-4), 96)
    };
    const auto lowerHalfRoot = approximation::enclosePrincipalComplexSqrt(
        minusThreeMinusFourI, 96);
    tests.expect(lowerHalfRoot.real().contains(rational(1))
        && lowerHalfRoot.imaginary().contains(rational(-2)),
        "CertifiedComplexSqrt: lower-half-plane input selects negative imaginary root");

    const approximation::ComplexInterval negativeRealAxis{
        approximation::RealInterval::fromRational(rational(-4), 96),
        approximation::RealInterval::fromRational(rational(0), 96)
    };
    const auto cutValue = approximation::enclosePrincipalComplexSqrt(
        negativeRealAxis, 96);
    tests.expect(cutValue.real().contains(rational(0))
        && cutValue.imaginary().contains(rational(2))
        && cutValue.imaginary().lower() >= numeric::BigFloat{},
        "CertifiedComplexSqrt: negative real axis uses the positive-imaginary principal value");
}

} // namespace mmcal::tests
