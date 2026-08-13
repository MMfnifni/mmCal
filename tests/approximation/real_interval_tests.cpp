// 実区間の回帰テスト
#include "real_interval_tests.hpp"

#include "approximation/real_interval.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "test_framework.hpp"

#include <stdexcept>

namespace mmcal::tests {
namespace {

using approximation::RealInterval;
using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

} // namespace

void runRealIntervalTests(TestRunner& tests) {
    const auto third = RealInterval::fromRational(rational(1, 3), 12);
    tests.expect(third.contains(rational(1, 3)),
        "RealInterval: directed conversion contains exact rational");
    tests.expect(!third.isPoint(),
        "RealInterval: non-dyadic rational is not falsely represented as a point");

    const auto half = RealInterval::fromRational(rational(1, 2), 12);
    tests.expect(half.isPoint(),
        "RealInterval: exact dyadic rational becomes a point interval");

    const auto negativePoint = RealInterval::point(
        numeric::BigFloat::fromRational(rational(-3, 8), 16));
    tests.expect(negativePoint.isPoint() && negativePoint.contains(rational(-3, 8)),
        "RealInterval: negative point construction is copy/move order independent");

    const auto sum = approximation::add(third, half, 10);
    tests.expect(sum.contains(rational(5, 6)),
        "RealInterval: outward addition encloses exact result");

    const auto difference = approximation::subtract(third, half, 10);
    tests.expect(difference.contains(rational(-1, 6)),
        "RealInterval: outward subtraction encloses exact result");

    const auto product = approximation::multiply(third, half, 10);
    tests.expect(product.contains(rational(1, 6)),
        "RealInterval: outward multiplication encloses exact result");

    const auto quotient = approximation::divide(third, half, 10);
    tests.expect(quotient.contains(rational(2, 3)),
        "RealInterval: outward division encloses exact result");

    const auto negativeThird = approximation::negate(third);
    tests.expect(negativeThird.contains(rational(-1, 3)),
        "RealInterval: negation reverses and negates endpoints");

    const auto joined = approximation::hull(third, negativeThird);
    tests.expect(joined.containsZero() && joined.contains(rational(1, 3))
        && joined.contains(rational(-1, 3)),
        "RealInterval: hull contains both source intervals");

    const auto aroundZero = RealInterval{
        numeric::BigFloat::fromBigInt(BigInt{-1}, 16),
        numeric::BigFloat::fromBigInt(BigInt{1}, 16)};
    tests.expectThrows<std::domain_error>([&] {
        static_cast<void>(approximation::divide(half, aroundZero, 16));
    }, "RealInterval: division rejects denominator interval containing zero");
}

} // namespace mmcal::tests
