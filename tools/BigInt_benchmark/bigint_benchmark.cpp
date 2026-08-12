// BigInt乗算のKaratsuba crossoverとfactorial性能を測る簡易benchmark
// GCC/Clang例:
// g++ -O3 -DNDEBUG -std=c++20 -DMMCAL_KARATSUBA_THRESHOLD_LIMBS=48 -Isrc \
//   tools/BigInt_benchmark/bigint_benchmark.cpp src/numeric/detail/big_uint.cpp \
//   src/numeric/big_int.cpp src/numeric/integer_algorithms.cpp -o bigint_benchmark
// MSVCでは /O2 /DNDEBUG /DMMCAL_KARATSUBA_THRESHOLD_LIMBS=48 を指定する。

#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {

using mmcal::numeric::BigInt;
using Clock = std::chrono::steady_clock;

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
    for (int iteration = 0; iteration < iterations; ++iteration) {
        const BigInt result = lhs * rhs;
        checksum += result.bitLength();
    }
    const auto end = Clock::now();

    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::micro>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkFactorial(std::uint64_t n, int iterations) {
    std::size_t checksum = 0;

    const auto start = Clock::now();
    for (int iteration = 0; iteration < iterations; ++iteration) {
        const BigInt result = mmcal::numeric::factorial(n);
        checksum += result.bitLength();
    }
    const auto end = Clock::now();

    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

} // namespace

int main() {
    std::cout << "multiply_us_per_op\n";
    for (const std::size_t limbs : {16u, 24u, 32u, 40u, 48u, 56u, 64u, 80u,
             96u, 128u, 192u, 256u, 384u, 512u}) {
        const int iterations = limbs <= 64 ? 5000 : limbs <= 192 ? 1500 : 300;
        std::cout << std::setw(4) << limbs << " limbs  "
                  << std::fixed << std::setprecision(3)
                  << benchmarkMultiply(limbs, iterations) << " us\n";
    }

    std::cout << "factorial_ms_per_op\n";
    for (const std::uint64_t n : {1000u, 2000u, 5000u, 10000u, 20000u, 40000u}) {
        const int iterations = n <= 10000 ? 4 : 2;
        std::cout << std::setw(5) << n << "!  "
                  << std::fixed << std::setprecision(3)
                  << benchmarkFactorial(n, iterations) << " ms\n";
    }
}
