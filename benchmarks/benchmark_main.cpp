// mmCalの性能測定・ランダム不変量試験を本体testから分離して実行する開発用runner
#include "approximation/certified_exponential.hpp"
#include "approximation/certified_logarithm.hpp"
#include "approximation/precision.hpp"
#include "approximation/real_interval.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
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
        std::cout << "Random checks: PASS\n";
    }

    if (!randomOnly)
        runBenchmarks(full);
    return 0;
}
