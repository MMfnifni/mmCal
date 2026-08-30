// 複素超越函数の保証付き評価の回帰テスト
#include "certified_complex_transcendental_tests.hpp"

#include "approximation/certified_complex_transcendental.hpp"
#include "approximation/complex_interval.hpp"
#include "approximation/real_interval.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "test_framework.hpp"

namespace mmcal::tests {
namespace {

[[nodiscard]] numeric::Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return numeric::Rational{numeric::BigInt{numerator}, numeric::BigInt{denominator}};
}

[[nodiscard]] approximation::RealInterval point(
    const numeric::Rational& value,
    std::size_t bits) {
    return approximation::RealInterval::fromRational(value, bits);
}

[[nodiscard]] approximation::ComplexInterval complexPoint(
    const numeric::Rational& real,
    const numeric::Rational& imaginary,
    std::size_t bits) {
    return approximation::ComplexInterval{point(real, bits), point(imaginary, bits)};
}

} // namespace

void runCertifiedComplexTranscendentalTests(TestRunner& tests) {
    constexpr std::size_t bits = 224;

    const auto arg = approximation::enclosePrincipalArgument(
        complexPoint(rational(2), rational(3), bits), bits);
    tests.expect(
        arg.interval.lower().toRational() > numeric::Rational::parse(
            "0.98279372324732906798571061101466601449687745363162855")
        && arg.interval.upper().toRational() < numeric::Rational::parse(
            "0.98279372324732906798571061101466601449687745363162857"),
        "Certified Arg: atan2 enclosure lies inside an independent decimal bracket");

    const auto log = approximation::enclosePrincipalComplexLog(
        complexPoint(rational(1), rational(1), bits), bits);
    tests.expect(
        log.interval.real().lower().toRational() > numeric::Rational::parse(
            "0.34657359027997265470861606072908828403775006718012762")
        && log.interval.real().upper().toRational() < numeric::Rational::parse(
            "0.34657359027997265470861606072908828403775006718012764"),
        "Certified complex Log: real part lies inside a log(sqrt(2)) bracket");
    tests.expect(
        log.interval.imaginary().lower().toRational() > numeric::Rational::parse(
            "0.78539816339744830961566084581987572104929234984377645")
        && log.interval.imaginary().upper().toRational() < numeric::Rational::parse(
            "0.78539816339744830961566084581987572104929234984377646"),
        "Certified complex Log: imaginary part lies inside a Pi/4 bracket");

    const auto nonDyadicLog = approximation::enclosePrincipalComplexLog(
        complexPoint(rational(-23, 5), rational(-2), bits), bits);
    tests.expect(
        nonDyadicLog.interval.real().lower().toRational() > numeric::Rational::parse(
            "1.61262771591611988770805581122609921021853894560817275")
        && nonDyadicLog.interval.real().upper().toRational() < numeric::Rational::parse(
            "1.61262771591611988770805581122609921021853894560817276"),
        "Certified complex Log: non-dyadic real component keeps a tight real enclosure");
    tests.expect(
        nonDyadicLog.interval.imaginary().lower().toRational() > numeric::Rational::parse(
            "-2.73146531304830224601143261242428361928463864052889394")
        && nonDyadicLog.interval.imaginary().upper().toRational() < numeric::Rational::parse(
            "-2.73146531304830224601143261242428361928463864052889393"),
        "Certified complex Log: non-dyadic Arg stays certified without Rational-series growth");

    const auto exponential = approximation::encloseComplexExp(
        complexPoint(rational(1), rational(1), bits), bits);
    tests.expect(
        exponential.interval.real().lower().toRational() > numeric::Rational::parse(
            "1.46869393991588515713896759732660426132695673662900872")
        && exponential.interval.real().upper().toRational() < numeric::Rational::parse(
            "1.46869393991588515713896759732660426132695673662900873"),
        "Certified complex Exp: real part lies inside an independent decimal bracket");
    tests.expect(
        exponential.interval.imaginary().lower().toRational() > numeric::Rational::parse(
            "2.28735528717884239120817190670050180895558625666835568")
        && exponential.interval.imaginary().upper().toRational() < numeric::Rational::parse(
            "2.28735528717884239120817190670050180895558625666835569"),
        "Certified complex Exp: imaginary part lies inside an independent decimal bracket");

    const auto principalPower = approximation::enclosePrincipalPower(
        complexPoint(rational(-8), rational(0), bits),
        complexPoint(rational(1, 3), rational(0), bits),
        bits);
    tests.expect(principalPower.interval.real().contains(rational(1)),
        "Certified Power: principal cube root of -8 contains real part 1");
    tests.expect(
        principalPower.interval.imaginary().lower().toRational() > numeric::Rational::parse(
            "1.73205080756887729352744634150587236694280525381038062")
        && principalPower.interval.imaginary().upper().toRational() < numeric::Rational::parse(
            "1.73205080756887729352744634150587236694280525381038064"),
        "Certified Power: principal cube root of -8 lies inside a sqrt(3) bracket");
}

} // namespace mmcal::tests
