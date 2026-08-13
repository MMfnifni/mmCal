// mmCalの性能測定・ランダム不変量試験を本体testから分離して実行する開発用runner
#include "approximation/certified_exponential.hpp"
#include "approximation/certified_logarithm.hpp"
#include "approximation/certified_special_functions.hpp"
#include "approximation/precision.hpp"
#include "approximation/real_interval.hpp"
#include "builtins/signal_processing.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "formatting/expr_formatter.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"
#include "numeric/number.hpp"
#include "symbols/symbol_table.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string_view>

namespace {

using Clock = std::chrono::steady_clock;
using mmcal::numeric::BigInt;
using mmcal::numeric::Rational;

[[nodiscard]] BigInt randomPositiveBigInt(
    std::mt19937_64& rng,
    std::size_t bits) {
    if (bits == 0)
        return BigInt{};

    BigInt value;
    std::size_t remaining = bits;
    while (remaining != 0) {
        const std::size_t chunkBits = std::min<std::size_t>(64, remaining);
        std::uint64_t chunk = rng();
        if (chunkBits < 64)
            chunk &= (std::uint64_t{1} << chunkBits) - 1;
        value <<= chunkBits;
        value += BigInt::fromUnsigned(chunk);
        remaining -= chunkBits;
    }

    value += BigInt{1} << (bits - 1);
    return value;
}

[[nodiscard]] BigInt makeOperand(std::size_t limbs, unsigned salt) {
    BigInt value{1};
    value <<= limbs * 32 - 1;
    value += BigInt{static_cast<std::int64_t>(0x1234567u + salt)};

    BigInt middle{1};
    middle <<= std::max<std::size_t>(1, limbs * 13);
    value += middle;
    return value;
}

[[nodiscard]] double benchmarkMultiply(std::size_t limbs, int iterations) {
    const BigInt lhs = makeOperand(limbs, 17);
    const BigInt rhs = makeOperand(limbs, 93);
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i)
        checksum += (lhs * rhs).bitLength();
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::micro>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkSquare(std::size_t limbs, int iterations) {
    const BigInt value = makeOperand(limbs, 41);
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i)
        checksum += (value * value).bitLength();
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::micro>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkDivision(std::size_t limbs, int iterations) {
    const BigInt divisor = makeOperand(limbs, 53);
    const BigInt quotient = makeOperand(limbs, 71);
    const BigInt dividend = divisor * quotient + BigInt{123};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i)
        checksum += (dividend / divisor).bitLength();
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::micro>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkFactorial(std::uint64_t n, int iterations) {
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i)
        checksum += mmcal::numeric::factorial(n).bitLength();
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

void benchmarkDecimalConversion(std::uint64_t n) {
    const BigInt value = mmcal::numeric::factorial(n);
    const auto formatStart = Clock::now();
    const std::string text = value.toString();
    const auto formatEnd = Clock::now();
    const BigInt parsed = BigInt::parse(text);
    const auto parseEnd = Clock::now();
    if (parsed != value)
        std::abort();
    std::cout << std::setw(6) << n << "!  " << std::setw(8) << text.size()
              << " digits  toString="
              << std::chrono::duration<double, std::milli>(formatEnd - formatStart).count()
              << " ms  parse="
              << std::chrono::duration<double, std::milli>(parseEnd - formatEnd).count()
              << " ms\n";
}

[[nodiscard]] double benchmarkExp(std::size_t digits, int iterations) {
    const std::size_t bits = mmcal::approximation::decimalDigitsToBinaryBits(digits + 12);
    const auto one = mmcal::approximation::RealInterval::fromRational(
        Rational{BigInt{1}}, bits);
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i)
        checksum += mmcal::approximation::encloseExp(one, bits).termsUsed;
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkLog2(std::size_t digits, int iterations) {
    const std::size_t bits = mmcal::approximation::decimalDigitsToBinaryBits(digits + 12);
    const auto two = mmcal::approximation::RealInterval::fromRational(
        Rational{BigInt{2}}, bits);
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i)
        checksum += mmcal::approximation::encloseLogPositive(two, bits).termsUsed;
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkGenericExp(std::size_t digits) {
    const std::size_t bits = mmcal::approximation::decimalDigitsToBinaryBits(digits + 12);
    const Rational value{BigInt{123456789}, BigInt{987654321}};
    const auto input = mmcal::approximation::RealInterval::fromRational(value, bits);
    const auto start = Clock::now();
    const auto result = mmcal::approximation::encloseExp(input, bits);
    const auto end = Clock::now();
    if (result.termsUsed == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

[[nodiscard]] double benchmarkGenericLog(std::size_t digits) {
    const std::size_t bits = mmcal::approximation::decimalDigitsToBinaryBits(digits + 12);
    const Rational value{BigInt{123456789}, BigInt{987654321}};
    const auto input = mmcal::approximation::RealInterval::fromRational(value, bits);
    const auto start = Clock::now();
    const auto result = mmcal::approximation::encloseLogPositive(input, bits);
    const auto end = Clock::now();
    if (result.termsUsed == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

struct FourierFixture final {
    mmcal::symbols::SymbolTable symbols;
    mmcal::evaluation::BuiltinRegistry registry;
    mmcal::mathematics::MathRegistry mathematics;
    mmcal::mathematics::AngleSemantics angles;
    mmcal::builtins::FourierTransformCache cache;

    FourierFixture()
        : registry(mmcal::evaluation::BuiltinRegistry::defaults(symbols)),
          mathematics(mmcal::mathematics::MathRegistry::defaults(symbols, registry)),
          angles(mmcal::mathematics::defaultAngleSemantics()) {}
};

[[nodiscard]] mmcal::expression::Expr fourierInput(std::size_t size) {
    std::vector<mmcal::expression::Expr> values;
    values.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        const std::int64_t value = static_cast<std::int64_t>((i * 37 + 11) % 101) - 50;
        values.emplace_back(mmcal::numeric::Number{BigInt{value}});
    }
    return mmcal::expression::Expr::array({size}, std::move(values));
}

[[nodiscard]] double benchmarkExactFft(std::size_t size, int iterations) {
    FourierFixture fixture;
    const mmcal::expression::Expr input = fourierInput(size);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    static_cast<void>(mmcal::builtins::evaluateFft(
        arguments, fixture.registry, fixture.mathematics, fixture.angles, fixture.cache));

    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateFft(
            arguments, fixture.registry, fixture.mathematics, fixture.angles, fixture.cache);
        checksum += result.asArray().elements.size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkApproximateDft(
    std::size_t size,
    std::size_t digits,
    int iterations) {
    FourierFixture fixture;
    const mmcal::expression::Expr input = fourierInput(size);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateDft(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result)
            std::abort();
        checksum += result->asArray().elements.size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkApproximateFft(
    std::size_t size,
    std::size_t digits,
    int iterations) {
    FourierFixture fixture;
    const mmcal::expression::Expr input = fourierInput(size);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateFft(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result)
            std::abort();
        checksum += result->asArray().elements.size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] bool approximateContainsInteger(
    const mmcal::expression::Expr& value,
    const Rational& expected) {
    if (value.isDecimalApproximation()) {
        const auto& decimal = value.asDecimalApproximation();
        return decimal.certifiedLower() <= expected && expected <= decimal.certifiedUpper();
    }
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        const Rational zero{};
        return complex.real().certifiedLower() <= expected
            && expected <= complex.real().certifiedUpper()
            && complex.imaginary().certifiedLower() <= zero
            && zero <= complex.imaginary().certifiedUpper();
    }
    if (value.isNumber() && value.asNumber().isReal())
        return value.asNumber().asReal().toRational() == expected;
    return false;
}

[[nodiscard]] bool runRandomFourierChecks(std::size_t count) {
    FourierFixture fixture;
    std::mt19937_64 rng{0x46465443455254ULL};
    constexpr std::array<std::size_t, 10> sizes{3, 5, 7, 8, 15, 16, 17, 31, 33, 65};

    for (std::size_t caseIndex = 0; caseIndex < count; ++caseIndex) {
        const std::size_t n = sizes[rng() % sizes.size()];
        std::vector<mmcal::expression::Expr> values;
        values.reserve(n);
        std::vector<Rational> expected;
        expected.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            const std::int64_t integer = static_cast<std::int64_t>(rng() % 31) - 15;
            values.emplace_back(mmcal::numeric::Number{BigInt{integer}});
            expected.emplace_back(BigInt{integer});
        }
        const auto input = mmcal::expression::Expr::array({n}, std::move(values));
        const std::array<mmcal::expression::Expr, 1> arguments{input};
        const auto transformed = mmcal::builtins::evaluateApproximateFft(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!transformed)
            return false;

        const std::array<mmcal::expression::Expr, 1> inverseArguments{*transformed};
        const auto roundTrip = mmcal::builtins::evaluateApproximateIfft(
            inverseArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!roundTrip || roundTrip->asArray().elements.size() != n)
            return false;

        for (std::size_t i = 0; i < n; ++i)
            if (!approximateContainsInteger(roundTrip->asArray().elements[i], expected[i]))
                return false;
    }
    return true;
}

[[nodiscard]] bool runRandomBigIntChecks(std::size_t count) {
    std::mt19937_64 rng{0x4D4D43414CULL};
    for (std::size_t i = 0; i < count; ++i) {
        const std::size_t divisorBits = 96 + (rng() % 8096);
        const std::size_t quotientBits = 64 + (rng() % 8096);
        const BigInt divisor = randomPositiveBigInt(rng, divisorBits);
        const BigInt expectedQuotient = randomPositiveBigInt(rng, quotientBits);
        const BigInt expectedRemainder = BigInt::fromUnsigned(rng() & 0xffffu);
        const BigInt dividend = divisor * expectedQuotient + expectedRemainder;

        const auto dm = mmcal::numeric::divmod(dividend, divisor);
        if (dm.quotient != expectedQuotient || dm.remainder != expectedRemainder)
            return false;
        if (dm.quotient * divisor + dm.remainder != dividend)
            return false;

        const BigInt product = divisor * expectedQuotient;
        if (product / divisor != expectedQuotient || product % divisor != BigInt{})
            return false;

        const std::string decimal = dividend.toString();
        if (BigInt::parse(decimal) != dividend)
            return false;
    }
    return true;
}

[[nodiscard]] bool runRandomCertifiedChecks(std::size_t count) {
    std::mt19937_64 rng{0x4558504C4F47ULL};
    constexpr std::size_t bits = 192;
    const Rational one{BigInt{1}};
    const Rational zero{BigInt{0}};

    for (std::size_t i = 0; i < count; ++i) {
        const std::int64_t numerator = static_cast<std::int64_t>(rng() % 161) - 80;
        const std::int64_t denominator = 1 + static_cast<std::int64_t>(rng() % 31);
        const Rational x{BigInt{numerator}, BigInt{denominator}};
        const auto xInterval = mmcal::approximation::RealInterval::fromRational(x, bits);
        const auto minusXInterval = mmcal::approximation::RealInterval::fromRational(-x, bits);
        const auto expX = mmcal::approximation::encloseExp(xInterval, bits);
        const auto expMinusX = mmcal::approximation::encloseExp(minusXInterval, bits);
        const auto expProduct = mmcal::approximation::multiply(
            expX.interval, expMinusX.interval, bits - 16);
        if (!expProduct.contains(one))
            return false;

        const Rational positive{
            BigInt{1 + static_cast<std::int64_t>(rng() % 160)},
            BigInt{1 + static_cast<std::int64_t>(rng() % 31)}};
        const Rational reciprocal = one / positive;
        const auto logPositive = mmcal::approximation::encloseLogPositive(
            mmcal::approximation::RealInterval::fromRational(positive, bits), bits);
        const auto logReciprocal = mmcal::approximation::encloseLogPositive(
            mmcal::approximation::RealInterval::fromRational(reciprocal, bits), bits);
        const auto logSum = mmcal::approximation::add(
            logPositive.interval, logReciprocal.interval, bits - 16);
        if (!logSum.contains(zero))
            return false;
    }
    return true;
}

[[nodiscard]] bool runRandomSpecialFunctionChecks(std::size_t count) {
    std::mt19937_64 rng{0x324631454C4C4950ULL};
    constexpr std::size_t bits = 192;
    const Rational one{BigInt{1}};

    for (std::size_t i = 0; i < count; ++i) {
        // 2F1(a,b;b,z)=(1-z)^(-a) のうち a=1 を使い、
        // randomなsafe-domain pointでcertified enclosureがexact値を含むことを監視する。
        const Rational b{BigInt{1 + static_cast<std::int64_t>(rng() % 9)},
            BigInt{1 + static_cast<std::int64_t>(rng() % 7)}};
        const Rational z{BigInt{static_cast<std::int64_t>(rng() % 15) - 7}, BigInt{16}};
        const Rational expected2F1 = one / (one - z);
        const auto hyper = mmcal::approximation::encloseHypergeometric2F1Real(
            one, b, b, z, bits);
        if (!hyper.contains(expected2F1))
            return false;

        // m=0では三種の不完全楕円積分が振幅phiへexactに退化する。
        // symbolic簡約を経由せず数値backend自体の包含性をrandom pointで確認する。
        const Rational phi{BigInt{static_cast<std::int64_t>(rng() % 31) - 15}, BigInt{8}};
        if (!mmcal::approximation::encloseEllipticFReal(phi, Rational{BigInt{0}}, bits).contains(phi))
            return false;
        if (!mmcal::approximation::encloseEllipticEReal(phi, Rational{BigInt{0}}, bits).contains(phi))
            return false;
        if (!mmcal::approximation::encloseEllipticPiReal(
                Rational{BigInt{0}}, phi, Rational{BigInt{0}}, bits).contains(phi))
            return false;
    }
    return true;
}

void runBenchmarks(bool full) {
    std::cout << "BigInt multiply/division\n";
    for (const std::size_t limbs : full
        ? std::initializer_list<std::size_t>{64, 128, 256, 512, 1024, 2048, 4096}
        : std::initializer_list<std::size_t>{64, 128, 256, 512, 1024}) {
        const int iterations = limbs <= 256 ? 200 : limbs <= 1024 ? 40 : 8;
        std::cout << std::setw(5) << limbs << " limbs  mul="
                  << std::fixed << std::setprecision(3)
                  << benchmarkMultiply(limbs, iterations) << " us  square="
                  << benchmarkSquare(limbs, iterations) << " us  div="
                  << benchmarkDivision(limbs, iterations) << " us\n";
    }

    std::cout << "factorial\n";
    for (const std::uint64_t n : full
        ? std::initializer_list<std::uint64_t>{10000, 20000, 40000, 80000}
        : std::initializer_list<std::uint64_t>{10000, 20000, 40000}) {
        std::cout << std::setw(6) << n << "!  "
                  << benchmarkFactorial(n, n <= 10000 ? 3 : 1) << " ms\n";
    }

    std::cout << "decimal conversion\n";
    for (const std::uint64_t n : full
        ? std::initializer_list<std::uint64_t>{10000, 20000, 40000}
        : std::initializer_list<std::uint64_t>{10000, 20000})
        benchmarkDecimalConversion(n);

    std::cout << "exact/certified FFT\n";
    for (const std::size_t size : full
        ? std::initializer_list<std::size_t>{32, 64, 128, 256, 512}
        : std::initializer_list<std::size_t>{32, 64, 128}) {
        const int exactIterations = size <= 64 ? 3 : 1;
        const int approximateIterations = size <= 128 ? 5 : 2;
        std::cout << std::setw(5) << size << " points  exact="
                  << benchmarkExactFft(size, exactIterations) << " ms  N[...,16]="
                  << benchmarkApproximateFft(size, 16, approximateIterations) << " ms\n";
    }
    for (const std::size_t size : full
        ? std::initializer_list<std::size_t>{65, 127, 257, 509}
        : std::initializer_list<std::size_t>{65, 127})
        std::cout << std::setw(5) << size << " points  direct="
                  << benchmarkApproximateDft(size, 16, 2) << " ms  FFT="
                  << benchmarkApproximateFft(size, 16, 2) << " ms\n";

    std::cout << "certified exp/log\n";
    for (const std::size_t digits : full
        ? std::initializer_list<std::size_t>{100, 500, 1000, 2000, 5000, 10000}
        : std::initializer_list<std::size_t>{100, 500, 1000, 2000}) {
        const int iterations = digits <= 500 ? 5 : digits <= 2000 ? 2 : 1;
        std::cout << std::setw(5) << digits << " digits  exp(1)="
                  << benchmarkExp(digits, iterations) << " ms  log(2)="
                  << benchmarkLog2(digits, iterations) << " ms\n";
    }

    if (full) {
        std::cout << "generic certified rational (123456789/987654321)\n";
        for (const std::size_t digits : {1000u, 2000u, 5000u})
            std::cout << std::setw(5) << digits << " digits  exp="
                      << benchmarkGenericExp(digits) << " ms  log="
                      << benchmarkGenericLog(digits) << " ms\n";
    }
}

void printUsage() {
    std::cout
        << "mmCal.Benchmarks [--full] [--random-only] [--benchmark-only]\n"
        << "  default          quick random checks + quick benchmarks\n"
        << "  --full           larger random set and high-precision benchmarks\n"
        << "  --random-only    run invariant checks only\n"
        << "  --benchmark-only run timings only\n";
}

} // namespace

int main(int argc, char** argv) {
    bool full = false;
    bool randomOnly = false;
    bool benchmarkOnly = false;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--full")
            full = true;
        else if (arg == "--random-only")
            randomOnly = true;
        else if (arg == "--benchmark-only")
            benchmarkOnly = true;
        else if (arg == "--help" || arg == "-h") {
            printUsage();
            return 0;
        }
        else {
            std::cerr << "Unknown option: " << arg << '\n';
            printUsage();
            return 2;
        }
    }

    if (randomOnly && benchmarkOnly) {
        std::cerr << "--random-only and --benchmark-only cannot be combined\n";
        return 2;
    }

    if (!benchmarkOnly) {
        const std::size_t bigIntCases = full ? 2000 : 300;
        const std::size_t certifiedCases = full ? 200 : 40;
        std::cout << "Random BigInt checks: " << bigIntCases << " cases\n";
        if (!runRandomBigIntChecks(bigIntCases)) {
            std::cerr << "Random BigInt check failed\n";
            return 1;
        }
        std::cout << "Random certified exp/log checks: " << certifiedCases << " cases\n";
        if (!runRandomCertifiedChecks(certifiedCases)) {
            std::cerr << "Random certified check failed\n";
            return 1;
        }
        std::cout << "Random 2F1/elliptic checks: " << certifiedCases << " cases\n";
        if (!runRandomSpecialFunctionChecks(certifiedCases)) {
            std::cerr << "Random special-function check failed\n";
            return 1;
        }
        std::cout << "Random certified FFT checks: " << certifiedCases << " cases\n";
        if (!runRandomFourierChecks(certifiedCases)) {
            std::cerr << "Random certified FFT check failed\n";
            return 1;
        }
        std::cout << "Random checks: PASS\n";
    }

    if (!randomOnly)
        runBenchmarks(full);
    return 0;
}
