// BigInt乗算算法のcrossover、factorial、10進変換を測る簡易benchmark
// GCC/Clang例:
// g++ -O3 -DNDEBUG -std=c++20 -DMMCAL_TOOM3_THRESHOLD_LIMBS=1280 \
//   -DMMCAL_TOOM3_RECURSIVE_THRESHOLD_LIMBS=448 -Isrc \
//   tools/BigInt_benchmark/bigint_benchmark.cpp src/numeric/detail/big_uint.cpp \
//   src/numeric/big_int.cpp src/numeric/integer_algorithms.cpp -o bigint_benchmark
// Prime-Swing比較時は -DMMCAL_USE_PRIME_SWING_FACTORIAL を追加する。
// MSVCでは /O2 /DNDEBUG と同じ /D 定義を指定する。
// Karatsuba scratch/workspace化はpool型・再帰深度型ともGCC実測で遅くなったため本番採用していない。
// allocator特性が異なるMSVCでは、このbenchmarkで再評価してから導入する。

#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

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

void benchmarkDecimalConversion(std::uint64_t n) {
    const BigInt value = mmcal::numeric::factorial(n);

    const auto formatStart = Clock::now();
    const std::string text = value.toString();
    const auto formatEnd = Clock::now();
    const BigInt parsed = BigInt::parse(text);
    const auto parseEnd = Clock::now();

    if (parsed != value)
        std::abort();

    std::cout << std::setw(6) << n << "!  "
              << std::setw(8) << text.size() << " digits  toString="
              << std::fixed << std::setprecision(3)
              << std::chrono::duration<double, std::milli>(formatEnd - formatStart).count()
              << " ms  parse="
              << std::chrono::duration<double, std::milli>(parseEnd - formatEnd).count()
              << " ms\n";
}

} // namespace

int main() {
    std::cout << "multiply_us_per_op\n";
    for (const std::size_t limbs : {16u, 24u, 32u, 40u, 48u, 56u, 64u, 96u,
             128u, 256u, 512u, 768u, 1024u, 1152u, 1536u, 2048u, 4096u, 8192u}) {
        const int iterations = limbs <= 128 ? 2000 : limbs <= 1024 ? 300 : limbs <= 4096 ? 60 : 15;
        std::cout << std::setw(5) << limbs << " limbs  "
                  << std::fixed << std::setprecision(3)
                  << benchmarkMultiply(limbs, iterations) << " us\n";
    }

    std::cout << "factorial_ms_per_op\n";
    for (const std::uint64_t n : {1000u, 5000u, 10000u, 20000u, 40000u, 80000u}) {
        const int iterations = n <= 10000 ? 4 : n <= 40000 ? 2 : 1;
        std::cout << std::setw(6) << n << "!  "
                  << std::fixed << std::setprecision(3)
                  << benchmarkFactorial(n, iterations) << " ms\n";
    }

    std::cout << "decimal_conversion\n";
    for (const std::uint64_t n : {10000u, 20000u, 40000u})
        benchmarkDecimalConversion(n);
}
