// mmCalの性能測定・ランダム不変量試験を本体testから分離して実行する開発用runner
#include "approximation/certified_exponential.hpp"
#include "approximation/certified_logarithm.hpp"
#include "approximation/certified_special_functions.hpp"
#include "approximation/precision.hpp"
#include "approximation/real_interval.hpp"
#include "builtins/signal_processing.hpp"
#include "builtins/linear_algebra.hpp"
#include "linear_algebra/decomposition.hpp"
#include "linear_algebra/fraction_free_elimination.hpp"
#include "linear_algebra/modular_linear_algebra.hpp"
#include "kernel/kernel_session.hpp"
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
#include "symbolic/algebraic_number.hpp"
#include "symbolic/number_field.hpp"
#include "certification_boundary_fuzzer.hpp"
#include "performance_cliff_audit.hpp"
#include "random_expression_fuzzer.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <span>
#include <string>
#include <string_view>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using mmcal::numeric::BigInt;
using mmcal::numeric::Rational;
using mmcal::numeric::RealNumber;

class BenchmarkFailure final : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

[[noreturn]] void reportUnavailableBenchmarkResult(
    std::string_view operation,
    std::size_t size,
    std::size_t digits) {
    throw BenchmarkFailure{
        std::string{operation} + " returned no certified result (size="
        + std::to_string(size) + ", digits=" + std::to_string(digits) + ")"};
}

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
        checksum += result.asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

struct ExactFftRoundTripTiming final {
    double firstMilliseconds = 0.0;
    double warmMilliseconds = 0.0;
};

[[nodiscard]] ExactFftRoundTripTiming benchmarkExactFftRoundTrip(
    std::size_t size,
    std::size_t iterations) {
    FourierFixture fixture;
    const mmcal::expression::Expr input = fourierInput(size);
    const std::array<mmcal::expression::Expr, 1> forwardArguments{input};

    const auto firstStart = Clock::now();
    const auto firstForward = mmcal::builtins::evaluateFft(
        forwardArguments, fixture.registry, fixture.mathematics, fixture.angles, fixture.cache);
    const std::array<mmcal::expression::Expr, 1> firstInverseArguments{firstForward};
    const auto firstRoundTrip = mmcal::builtins::evaluateIfft(
        firstInverseArguments, fixture.registry, fixture.mathematics, fixture.angles, fixture.cache);
    const auto firstEnd = Clock::now();
    if (firstRoundTrip != input)
        std::abort();

    double warmMilliseconds = 0.0;
    for (std::size_t i = 0; i < iterations; ++i) {
        const auto start = Clock::now();
        const auto forward = mmcal::builtins::evaluateFft(
            forwardArguments, fixture.registry, fixture.mathematics, fixture.angles, fixture.cache);
        const std::array<mmcal::expression::Expr, 1> inverseArguments{forward};
        const auto roundTrip = mmcal::builtins::evaluateIfft(
            inverseArguments, fixture.registry, fixture.mathematics, fixture.angles, fixture.cache);
        const auto end = Clock::now();
        if (roundTrip != input)
            std::abort();
        warmMilliseconds += std::chrono::duration<double, std::milli>(end - start).count();
    }

    return ExactFftRoundTripTiming{
        std::chrono::duration<double, std::milli>(firstEnd - firstStart).count(),
        warmMilliseconds / static_cast<double>(iterations)};
}

void runExactCyclotomicFftBenchmark(std::size_t iterations) {
    std::cout << "exact cyclotomic FFT round-trip\n";
    for (const std::size_t size : std::initializer_list<std::size_t>{5, 7, 10, 12, 15, 21}) {
        const auto timing = benchmarkExactFftRoundTrip(size, iterations);
        std::cout << "  " << std::setw(3) << size << " points"
                  << "  first=" << timing.firstMilliseconds << " ms"
                  << "  warm=" << timing.warmMilliseconds << " ms\n";
    }
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
        checksum += result->asArray().size();
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
        checksum += result->asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

struct ApproximateFftCrossoverTiming final {
    double directMilliseconds = 0.0;
    double bluesteinMilliseconds = 0.0;
};

[[nodiscard]] ApproximateFftCrossoverTiming benchmarkApproximateFftCrossover(
    std::size_t size,
    std::size_t digits,
    std::size_t iterations) {
    FourierFixture fixture;
    const mmcal::expression::Expr input = fourierInput(size);
    const std::array<mmcal::expression::Expr, 1> arguments{input};

    std::optional<mmcal::expression::Expr> directResult;
    const auto directStart = Clock::now();
    for (std::size_t i = 0; i < iterations; ++i) {
        directResult = mmcal::builtins::evaluateApproximateDft(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!directResult)
            std::abort();
    }
    const auto directEnd = Clock::now();

    std::optional<mmcal::expression::Expr> bluesteinResult;
    const auto bluesteinStart = Clock::now();
    for (std::size_t i = 0; i < iterations; ++i) {
        bluesteinResult = mmcal::builtins::evaluateApproximateBluesteinFftForBenchmark(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!bluesteinResult)
            std::abort();
    }
    const auto bluesteinEnd = Clock::now();

    const auto displayedComponents = [](const mmcal::expression::Expr& value)
        -> std::optional<std::pair<Rational, Rational>> {
        if (value.isDecimalApproximation())
            return std::pair{value.asDecimalApproximation().displayedValue(), Rational{}};
        if (value.isComplexDecimalApproximation()) {
            const auto& complex = value.asComplexDecimalApproximation();
            return std::pair{
                complex.real().displayedValue(), complex.imaginary().displayedValue()};
        }
        return std::nullopt;
    };
    const auto sameDisplayedArray = [&] {
        const auto& directArray = directResult->asArray();
        const auto& bluesteinArray = bluesteinResult->asArray();
        if (directArray.shape != bluesteinArray.shape)
            return false;
        for (std::size_t index = 0; index < directArray.size(); ++index) {
            const auto direct = displayedComponents(directArray.element(index));
            const auto bluestein = displayedComponents(bluesteinArray.element(index));
            if (!direct || !bluestein || *direct != *bluestein)
                return false;
        }
        return true;
    };

    if (!sameDisplayedArray()) {
        const auto& directArray = directResult->asArray();
        const auto& bluesteinArray = bluesteinResult->asArray();
        for (std::size_t index = 0;
             index < std::min(directArray.size(), bluesteinArray.size()); ++index) {
            const auto direct = displayedComponents(directArray.element(index));
            const auto bluestein = displayedComponents(bluesteinArray.element(index));
            if (!direct || !bluestein || *direct != *bluestein) {
                throw BenchmarkFailure{
                    "Direct DFT and forced Bluestein differ at size "
                    + std::to_string(size) + ", index " + std::to_string(index)
                    + ": direct=" + mmcal::formatting::formatExpr(directArray.element(index))
                    + ", Bluestein="
                    + mmcal::formatting::formatExpr(bluesteinArray.element(index))};
            }
        }
        throw BenchmarkFailure{
            "Direct DFT and forced Bluestein result shapes differ at size "
            + std::to_string(size)};
    }

    return ApproximateFftCrossoverTiming{
        std::chrono::duration<double, std::milli>(directEnd - directStart).count()
            / static_cast<double>(iterations),
        std::chrono::duration<double, std::milli>(bluesteinEnd - bluesteinStart).count()
            / static_cast<double>(iterations)};
}

void runApproximateFftThresholdBenchmark(std::size_t iterations) {
    if (iterations == 0)
        throw std::invalid_argument("FFT threshold benchmark iterations must be positive");

    constexpr std::array<std::size_t, 15> sizes{
        65, 95, 127, 191, 255, 257, 319, 335, 351, 367, 383, 384, 385, 447, 509};
    std::cout << "certified FFT direct/Bluestein threshold sweep"
              << " (16 digits, " << iterations << " iteration"
              << (iterations == 1 ? "" : "s") << ")\n"
              << "current policy: direct below "
              << mmcal::builtins::approximateFftBluesteinThreshold
              << " points, Bluestein at or above it\n";
    for (const std::size_t size : sizes) {
        const ApproximateFftCrossoverTiming timing =
            benchmarkApproximateFftCrossover(size, 16, iterations);
        std::cout << std::setw(5) << size << " points  direct="
                  << std::fixed << std::setprecision(3) << timing.directMilliseconds
                  << " ms  Bluestein=" << timing.bluesteinMilliseconds
                  << " ms  faster="
                  << (timing.directMilliseconds <= timing.bluesteinMilliseconds
                      ? "direct" : "Bluestein") << '\n';
    }
}

struct MatrixFixture final {
    mmcal::symbols::SymbolTable symbols;
    mmcal::evaluation::BuiltinRegistry registry;
    mmcal::mathematics::MathRegistry mathematics;
    mmcal::mathematics::AngleSemantics angles;

    MatrixFixture()
        : registry(mmcal::evaluation::BuiltinRegistry::defaults(symbols)),
          mathematics(mmcal::mathematics::MathRegistry::defaults(symbols, registry)),
          angles(mmcal::mathematics::defaultAngleSemantics()) {}
};

[[nodiscard]] mmcal::expression::Expr matrixInput(std::size_t size, unsigned salt) {
    mmcal::expression::ArrayBuilder builder;
    builder.reserve(size * size);
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            const std::int64_t value = row == column
                ? static_cast<std::int64_t>(4 * size + 1 + salt % 3)
                : static_cast<std::int64_t>((row * 17 + column * 29 + salt) % 5) - 2;
            builder.append(BigInt{value});
        }
    }
    return mmcal::expression::Expr::array(builder.finish({size, size}));
}

// Eigen / eigensystemのtimingにはsimple spectrumを構成的に保証した密実対称行列を使う。
// 非対角成分の絶対値は高々2，隣接する対角値の間隔は4n+1なので，各Gershgorin円板は
// 互いに素になる。したがって各円板は固有値をちょうど1個含み，重複固有値は生じない。
[[nodiscard]] mmcal::expression::Expr simpleSpectrumMatrixInput(
    std::size_t size,
    unsigned salt) {
    mmcal::expression::ArrayBuilder builder;
    builder.reserve(size * size);
    const std::int64_t spacing = static_cast<std::int64_t>(4 * size + 1);
    const std::int64_t diagonalBase = spacing + static_cast<std::int64_t>(salt % 3);
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            std::int64_t value = 0;
            if (row == column) {
                value = diagonalBase + spacing * static_cast<std::int64_t>(row);
            }
            else {
                const std::size_t low = std::min(row, column);
                const std::size_t high = std::max(row, column);
                value = static_cast<std::int64_t>((low * 17 + high * 29 + salt) % 5) - 2;
            }
            builder.append(BigInt{value});
        }
    }
    return mmcal::expression::Expr::array(builder.finish({size, size}));
}

[[nodiscard]] mmcal::expression::Expr randomDecimalMatrix(
    std::size_t size,
    std::uint64_t seed = 1234) {
    // 添付fft_gen.pyと同じ [-1,1] / 小数10桁という負荷特性を，
    // parser時間と算法時間を分離するためexact Rationalとして直接構築する。
    constexpr std::int64_t scale = 10000000000LL;
    std::mt19937_64 rng{seed};
    std::uniform_int_distribution<std::int64_t> distribution{-scale, scale};
    mmcal::expression::ArrayBuilder builder;
    builder.reserve(size * size);
    for (std::size_t i = 0; i < size * size; ++i)
        builder.append(Rational{BigInt{distribution(rng)}, BigInt{scale}});
    return mmcal::expression::Expr::array(builder.finish({size, size}));
}


[[nodiscard]] mmcal::expression::Expr rowSumVector(
    const mmcal::expression::Expr& matrix) {
    const auto& array = matrix.asArray();
    const std::size_t size = array.shape[0];
    std::vector<mmcal::expression::Expr> values;
    values.reserve(size);
    for (std::size_t row = 0; row < size; ++row) {
        Rational sum{BigInt{0}};
        for (std::size_t column = 0; column < size; ++column)
            sum += array.exactNumber(row * size + column).realPart().toRational();
        values.emplace_back(mmcal::numeric::Number{std::move(sum)});
    }
    return mmcal::expression::Expr::array({size}, std::move(values));
}

[[nodiscard]] double benchmarkExactDot(std::size_t size, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr lhs = matrixInput(size, 3);
    const mmcal::expression::Expr rhs = matrixInput(size, 11);
    const std::array<mmcal::expression::Expr, 2> arguments{lhs, rhs};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateDot(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum += result.asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkExactDeterminant(std::size_t size, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 5);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateDeterminant(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum += result.asNumber().asReal().toRational().numerator().bitLength();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkExactRref(std::size_t size, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 7);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateRref(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum += result.asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] mmcal::expression::Expr integerSolutionVector(std::size_t size) {
    std::vector<mmcal::expression::Expr> values;
    values.reserve(size);
    for (std::size_t i = 0; i < size; ++i)
        values.emplace_back(mmcal::numeric::Number{BigInt{static_cast<std::int64_t>(i + 1)}});
    return mmcal::expression::Expr::array({size}, std::move(values));
}

[[nodiscard]] double benchmarkExactSolveLinear(std::size_t size, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr matrix = matrixInput(size, 9);
    const mmcal::expression::Expr expected = integerSolutionVector(size);
    const std::array<mmcal::expression::Expr, 2> dotArguments{matrix, expected};
    const mmcal::expression::Expr rhs = mmcal::builtins::evaluateDot(
        dotArguments, fixture.registry, fixture.mathematics, fixture.angles);
    const std::array<mmcal::expression::Expr, 2> arguments{matrix, rhs};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateSolveLinear(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum += result.asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] mmcal::expression::Expr matrixWithDependentColumn(
    const mmcal::expression::Expr& matrix) {
    const auto& array = matrix.asArray();
    const std::size_t rows = array.shape[0];
    const std::size_t columns = array.shape[1];
    std::vector<mmcal::expression::Expr> values;
    values.reserve(rows * (columns + 1));

    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t column = 0; column < columns; ++column)
            values.push_back(array.element(row * columns + column));

        const auto first = array.exactNumber(row * columns);
        const auto second = columns > 1
            ? array.exactNumber(row * columns + 1)
            : mmcal::numeric::Number{};
        values.emplace_back(first + mmcal::numeric::Number{BigInt{2}} * second);
    }
    return mmcal::expression::Expr::array({rows, columns + 1}, std::move(values));
}

[[nodiscard]] double benchmarkExactNullSpace(std::size_t size, int iterations) {
    MatrixFixture fixture;
    const auto matrix = matrixWithDependentColumn(matrixInput(size, 13));
    const std::array<mmcal::expression::Expr, 1> arguments{matrix};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateNullSpace(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum += result.asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkExactLu(std::size_t size, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 19);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateLuDecomposition(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum += result.asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkExactQr(std::size_t size, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 17);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateQrDecomposition(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum += result.asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkApproximateLu(
    std::size_t size, std::size_t digits, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 19);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateLuDecomposition(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result)
            std::abort();
        checksum += result->asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkApproximateQr(
    std::size_t size, std::size_t digits, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 17);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateQrDecomposition(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result)
            std::abort();
        checksum += result->asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkApproximateSvd(
    std::size_t size, std::size_t digits, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 23);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateSingularValueDecomposition(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result)
            std::abort();
        checksum += result->isArray()
            ? result->asArray().size()
            : result->asList().elements.size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkApproximateEigenvalues(
    std::size_t size, std::size_t digits, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = simpleSpectrumMatrixInput(size, 29);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateEigenvalues(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result)
            reportUnavailableBenchmarkResult("N[eigenvalues]", size, digits);
        checksum += result->asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkApproximateEigensystem(
    std::size_t size, std::size_t digits, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = simpleSpectrumMatrixInput(size, 29);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateEigensystem(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result)
            reportUnavailableBenchmarkResult("N[eigensystem]", size, digits);
        checksum += result->isArray() ? result->asArray().size() : result->asList().elements.size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkApproximateQrBlock(
    std::size_t size, std::size_t digits, int iterations, std::size_t blockColumns) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 17);
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::linear_algebra::approximateQrDecompositionWithBlockSize(
            input.asArray(), fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits}, blockColumns);
        if (!result)
            std::abort();
        checksum += result->asArray().size();
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count() / iterations;
}

[[nodiscard]] double benchmarkApproximateDeterminant(
    std::size_t size, std::size_t digits, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 5);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateDeterminant(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result)
            std::abort();
        checksum += result->isDecimalApproximation() ? 1 : 2;
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

[[nodiscard]] mmcal::expression::Expr randomComplexMatrix(
    std::mt19937_64& rng,
    std::size_t size) {
    std::vector<mmcal::expression::Expr> values;
    values.reserve(size * size);
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            std::int64_t real = static_cast<std::int64_t>(rng() % 7) - 3;
            const std::int64_t imaginary = static_cast<std::int64_t>(rng() % 7) - 3;
            if (row == column)
                real += static_cast<std::int64_t>(3 * size + 1);
            values.emplace_back(mmcal::numeric::Number::complex(
                RealNumber{BigInt{real}}, RealNumber{BigInt{imaginary}}));
        }
    }
    return mmcal::expression::Expr::array({size, size}, std::move(values));
}

[[nodiscard]] mmcal::expression::Expr randomInvertibleMatrix(
    std::mt19937_64& rng,
    std::size_t size) {
    std::vector<mmcal::expression::Expr> values;
    values.reserve(size * size);
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            const std::int64_t value = row == column
                ? static_cast<std::int64_t>(4 * size + 1 + rng() % 5)
                : static_cast<std::int64_t>(rng() % 7) - 3;
            values.emplace_back(mmcal::numeric::Number{BigInt{value}});
        }
    }
    return mmcal::expression::Expr::array({size, size}, std::move(values));
}

// eigensystemには完全な固有vector基底が存在する入力を渡す。可逆性だけでは
// defective Jordan blockを排除できないため，互いに異なる対角値を持つ上三角行列を
// 構成し，distinct eigenvalueによる対角化可能性を保証する。
[[nodiscard]] mmcal::expression::Expr randomSimpleSpectrumMatrix(
    std::mt19937_64& rng,
    std::size_t size) {
    std::vector<mmcal::expression::Expr> values;
    values.reserve(size * size);
    const std::int64_t diagonalBase = static_cast<std::int64_t>(5 * size + 3);
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            std::int64_t value = 0;
            if (row == column)
                value = diagonalBase + static_cast<std::int64_t>(7 * row);
            else if (row < column)
                value = static_cast<std::int64_t>(rng() % 7) - 3;
            values.emplace_back(mmcal::numeric::Number{BigInt{value}});
        }
    }
    return mmcal::expression::Expr::array({size, size}, std::move(values));
}

[[nodiscard]] mmcal::expression::Expr rationallyScaledMatrix(
    const mmcal::expression::Expr& matrix) {
    const auto& array = matrix.asArray();
    const std::size_t rows = array.shape[0];
    const std::size_t columns = array.shape[1];
    std::vector<mmcal::expression::Expr> values;
    values.reserve(array.size());

    for (std::size_t row = 0; row < rows; ++row)
        for (std::size_t column = 0; column < columns; ++column) {
            const auto rational = array.exactNumber(row * columns + column)
                .asReal().toRational();
            const BigInt denominator = BigInt{static_cast<std::int64_t>((row + 2) * (column + 3))};
            values.emplace_back(mmcal::numeric::Number{Rational{
                rational.numerator(), rational.denominator() * denominator}});
        }
    return mmcal::expression::Expr::array({rows, columns}, std::move(values));
}

[[nodiscard]] bool isExactIdentity(const mmcal::expression::Expr& matrix) {
    if (!matrix.isArray() || !matrix.asArray().isMatrix()
        || matrix.asArray().shape[0] != matrix.asArray().shape[1])
        return false;
    const std::size_t size = matrix.asArray().shape[0];
    for (std::size_t row = 0; row < size; ++row)
        for (std::size_t column = 0; column < size; ++column) {
            const auto& array = matrix.asArray();
            if (!array.hasExactNumberStorage())
                return false;
            const mmcal::numeric::Number expected{BigInt{row == column ? 1 : 0}};
            if (!(array.exactNumber(row * size + column) == expected))
                return false;
        }
    return true;
}

[[nodiscard]] bool isExactZeroVector(const mmcal::expression::Expr& vector) {
    if (!vector.isArray() || !vector.asArray().isVector())
        return false;
    const auto& array = vector.asArray();
    if (!array.hasExactNumberStorage())
        return false;
    for (std::size_t i = 0; i < array.size(); ++i)
        if (!array.exactNumber(i).isZero())
            return false;
    return true;
}

[[nodiscard]] bool validatesNullSpace(
    const mmcal::expression::Expr& matrix,
    const mmcal::expression::Expr& basis,
    MatrixFixture& fixture) {
    if (!basis.isArray() || !basis.asArray().isMatrix())
        return false;
    const auto& basisArray = basis.asArray();
    const std::size_t variables = matrix.asArray().shape[1];
    if (basisArray.shape[1] != variables)
        return false;

    for (std::size_t row = 0; row < basisArray.shape[0]; ++row) {
        std::vector<mmcal::expression::Expr> vectorElements;
        vectorElements.reserve(variables);
        for (std::size_t column = 0; column < variables; ++column)
            vectorElements.push_back(basisArray.element(row * variables + column));
        const auto vector = mmcal::expression::Expr::array({variables}, std::move(vectorElements));
        const std::array<mmcal::expression::Expr, 2> arguments{matrix, vector};
        if (!isExactZeroVector(mmcal::builtins::evaluateDot(
                arguments, fixture.registry, fixture.mathematics, fixture.angles)))
            return false;
    }
    return true;
}


[[nodiscard]] mmcal::expression::Expr packedMatrixFactor(
    const mmcal::expression::Expr& packed,
    std::size_t factor) {
    const auto& array = packed.asArray();
    const std::size_t rows = array.shape[1];
    const std::size_t columns = array.shape[2];
    const std::size_t factorSize = rows * columns;
    return mmcal::expression::Expr::array(
        array.sliced({rows, columns}, factor * factorSize, factorSize));
}

[[nodiscard]] bool approximateContainsNumber(
    const mmcal::expression::Expr& approximate,
    const mmcal::numeric::Number& exact) {
    const Rational real = exact.realPart().toRational();
    const Rational imaginary = exact.imaginaryPart().toRational();

    if (approximate.isDecimalApproximation()) {
        if (!imaginary.isZero())
            return false;
        const auto& decimal = approximate.asDecimalApproximation();
        return decimal.certifiedLower() <= real && real <= decimal.certifiedUpper();
    }
    if (approximate.isComplexDecimalApproximation()) {
        const auto& complex = approximate.asComplexDecimalApproximation();
        return complex.real().certifiedLower() <= real
            && real <= complex.real().certifiedUpper()
            && complex.imaginary().certifiedLower() <= imaginary
            && imaginary <= complex.imaginary().certifiedUpper();
    }
    if (approximate.isNumber())
        return approximate.asNumber() == exact;
    return false;
}

[[nodiscard]] bool approximateMatrixContainsExact(
    const mmcal::expression::Expr& approximate,
    const mmcal::expression::Expr& exact) {
    if (!approximate.isArray() || !exact.isArray()
        || approximate.asArray().shape != exact.asArray().shape)
        return false;
    const auto& lhs = approximate.asArray();
    const auto& rhs = exact.asArray();
    if (!rhs.hasExactNumberStorage())
        return false;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (!approximateContainsNumber(lhs.element(i), rhs.exactNumber(i)))
            return false;
    }
    return true;
}

struct DisplayedComplex final {
    Rational real;
    Rational imaginary;
};

[[nodiscard]] std::optional<DisplayedComplex> displayedComplex(
    const mmcal::expression::Expr& value) {
    if (value.isNumber())
        return DisplayedComplex{
            value.asNumber().realPart().toRational(),
            value.asNumber().imaginaryPart().toRational()};
    if (value.isDecimalApproximation())
        return DisplayedComplex{
            value.asDecimalApproximation().displayedValue(), Rational{BigInt{0}}};
    if (value.isComplexDecimalApproximation())
        return DisplayedComplex{
            value.asComplexDecimalApproximation().real().displayedValue(),
            value.asComplexDecimalApproximation().imaginary().displayedValue()};
    return std::nullopt;
}

[[nodiscard]] DisplayedComplex multiplyDisplayed(
    const DisplayedComplex& lhs,
    const DisplayedComplex& rhs) {
    return DisplayedComplex{
        lhs.real * rhs.real - lhs.imaginary * rhs.imaginary,
        lhs.real * rhs.imaginary + lhs.imaginary * rhs.real};
}

[[nodiscard]] bool verifyDisplayedEigenRelation(
    const mmcal::expression::Expr& exactMatrix,
    const mmcal::expression::Expr& values,
    const mmcal::expression::Expr& vectors) {
    const auto& a = exactMatrix.asArray();
    const auto& lambdas = values.asArray();
    const auto& v = vectors.asArray();
    const std::size_t size = a.shape[0];
    const Rational tolerance{BigInt{1}, mmcal::numeric::pow(BigInt{10}, 10)};
    auto absoluteRational = [](const Rational& x) { return x < Rational{BigInt{0}} ? -x : x; };

    for (std::size_t column = 0; column < size; ++column) {
        const auto lambda = displayedComplex(lambdas.element(column));
        if (!lambda)
            return false;
        bool nonzeroVector = false;
        for (std::size_t row = 0; row < size; ++row) {
            DisplayedComplex av{Rational{BigInt{0}}, Rational{BigInt{0}}};
            for (std::size_t k = 0; k < size; ++k) {
                const auto aik = displayedComplex(a.element(row * size + k));
                const auto vk = displayedComplex(v.element(k * size + column));
                if (!aik || !vk)
                    return false;
                const DisplayedComplex product = multiplyDisplayed(*aik, *vk);
                av.real += product.real;
                av.imaginary += product.imaginary;
            }
            const auto vr = displayedComplex(v.element(row * size + column));
            if (!vr)
                return false;
            nonzeroVector = nonzeroVector
                || absoluteRational(vr->real) > tolerance
                || absoluteRational(vr->imaginary) > tolerance;
            const DisplayedComplex lv = multiplyDisplayed(*lambda, *vr);
            if (absoluteRational(av.real - lv.real) > tolerance
                || absoluteRational(av.imaginary - lv.imaginary) > tolerance)
                return false;
        }
        if (!nonzeroVector)
            return false;
    }
    return true;
}

[[nodiscard]] bool runRandomMatrixChecks(std::size_t count) {
    MatrixFixture fixture;
    std::mt19937_64 rng{0x4D41545249584345ULL};
    std::mt19937_64 eigenRng{0x454947454E434552ULL};

    // timingで使う全サイズをrandom-onlyでも直接監視し，測定開始後のabortを防ぐ。
    constexpr std::array<std::size_t, 4> eigenTimingSizes{4, 8, 12, 16};
    for (const std::size_t size : eigenTimingSizes) {
        const auto matrix = simpleSpectrumMatrixInput(size, 29);
        const std::array<mmcal::expression::Expr, 1> arguments{matrix};
        const auto system = mmcal::builtins::evaluateApproximateEigensystem(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{16});
        if (!system || !system->isList() || system->asList().elements.size() != 2) {
            std::cerr << "Fixed Eigen timing input failure: phase=shape size=" << size
                << " matrix=" << mmcal::formatting::formatExpr(matrix) << '\n';
            return false;
        }
        const auto& values = system->asList().elements[0];
        const auto& vectors = system->asList().elements[1];
        if (!values.isArray() || values.asArray().shape != std::vector<std::size_t>{size}
            || !vectors.isArray()
            || vectors.asArray().shape != std::vector<std::size_t>{size, size}
            || !verifyDisplayedEigenRelation(matrix, values, vectors)) {
            std::cerr << "Fixed Eigen timing input failure: phase=relation size=" << size
                << " matrix=" << mmcal::formatting::formatExpr(matrix) << '\n';
            return false;
        }
    }

    for (std::size_t caseIndex = 0; caseIndex < count; ++caseIndex) {
        const std::size_t size = 1 + rng() % 5;
        const auto matrix = randomInvertibleMatrix(rng, size);
        const std::array<mmcal::expression::Expr, 1> unary{matrix};
        const auto reportFailure = [&](std::string_view stage,
            const mmcal::expression::Expr& subject) {
            std::cerr << "Random certified Matrix failure: phase=" << stage
                << " case=" << (caseIndex + 1)
                << " case-index=" << caseIndex
                << " size=" << subject.asArray().shape[0]
                << " matrix=" << mmcal::formatting::formatExpr(subject) << '\n';
            return false;
        };
        const auto fail = [&](std::string_view stage) {
            return reportFailure(stage, matrix);
        };

        const auto determinant = mmcal::builtins::evaluateDeterminant(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        if (!determinant.isNumber() || !determinant.asNumber().isReal()
            || determinant.asNumber().isZero())
            return fail("exact-determinant");

        const auto transposed = mmcal::builtins::evaluateTranspose(unary, fixture.registry);
        const std::array<mmcal::expression::Expr, 1> transposedArguments{transposed};
        const auto transposedDeterminant = mmcal::builtins::evaluateDeterminant(
            transposedArguments, fixture.registry, fixture.mathematics, fixture.angles);
        if (!(transposedDeterminant == determinant))
            return fail("transpose-determinant");

        const auto inverse = mmcal::builtins::evaluateMatrixInverse(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> productArguments{matrix, inverse};
        if (!isExactIdentity(mmcal::builtins::evaluateDot(
                productArguments, fixture.registry, fixture.mathematics, fixture.angles)))
            return fail("exact-inverse");

        const auto lu = mmcal::builtins::evaluateLuDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        if (!lu.isArray() || lu.asArray().shape != std::vector<std::size_t>{3, size, size})
            return fail("exact-lu-shape");
        const auto permutation = packedMatrixFactor(lu, 0);
        const auto lower = packedMatrixFactor(lu, 1);
        const auto upper = packedMatrixFactor(lu, 2);
        const std::array<mmcal::expression::Expr, 2> paArguments{permutation, matrix};
        const std::array<mmcal::expression::Expr, 2> luArguments{lower, upper};
        if (!(mmcal::builtins::evaluateDot(
                paArguments, fixture.registry, fixture.mathematics, fixture.angles)
            == mmcal::builtins::evaluateDot(
                luArguments, fixture.registry, fixture.mathematics, fixture.angles)))
            return fail("exact-lu-reconstruction");

        const auto approximateQr = mmcal::builtins::evaluateApproximateQrDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!approximateQr || !approximateQr->isArray()
            || approximateQr->asArray().shape != std::vector<std::size_t>{2, size, size})
            return fail("certified-qr-shape");
        const auto q = packedMatrixFactor(*approximateQr, 0);
        const auto r = packedMatrixFactor(*approximateQr, 1);
        const std::array<mmcal::expression::Expr, 2> qrArguments{q, r};
        const auto reconstructed = mmcal::builtins::evaluateApproximateDot(
            qrArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{18});
        if (!reconstructed || !approximateMatrixContainsExact(*reconstructed, matrix))
            return fail("certified-qr-reconstruction");

        const std::array<mmcal::expression::Expr, 1> qUnary{q};
        const auto qt = mmcal::builtins::evaluateTranspose(qUnary, fixture.registry);
        const std::array<mmcal::expression::Expr, 2> orthogonalArguments{qt, q};
        const auto qtq = mmcal::builtins::evaluateApproximateDot(
            orthogonalArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{18});
        std::vector<mmcal::expression::Expr> identityElements;
        identityElements.reserve(size * size);
        for (std::size_t row = 0; row < size; ++row)
            for (std::size_t column = 0; column < size; ++column)
                identityElements.emplace_back(mmcal::numeric::Number{BigInt{row == column ? 1 : 0}});
        const auto identity = mmcal::expression::Expr::array({size, size}, std::move(identityElements));
        if (!qtq || !approximateMatrixContainsExact(*qtq, identity))
            return fail("certified-qr-orthogonality");

        const auto approximateSvd = mmcal::builtins::evaluateApproximateSingularValueDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!approximateSvd || !approximateSvd->isArray()
            || approximateSvd->asArray().shape != std::vector<std::size_t>{3, size, size})
            return fail("certified-svd-shape");
        const auto svdU = packedMatrixFactor(*approximateSvd, 0);
        const auto svdS = packedMatrixFactor(*approximateSvd, 1);
        const auto svdV = packedMatrixFactor(*approximateSvd, 2);
        const std::array<mmcal::expression::Expr, 2> usArguments{svdU, svdS};
        const auto us = mmcal::builtins::evaluateApproximateDot(
            usArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{18});
        if (!us)
            return fail("certified-svd-us-product");
        const std::array<mmcal::expression::Expr, 1> svdVUnary{svdV};
        const auto svdVt = mmcal::builtins::evaluateTranspose(svdVUnary, fixture.registry);
        const std::array<mmcal::expression::Expr, 2> svdReconstructArguments{*us, svdVt};
        const auto svdReconstructed = mmcal::builtins::evaluateApproximateDot(
            svdReconstructArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{18});
        if (!svdReconstructed || !approximateMatrixContainsExact(*svdReconstructed, matrix))
            return fail("certified-svd-reconstruction");

        // SVD backend自体は20桁でreconstruction/直交性をcertifyする。
        // ここでは20桁へ丸めて公開されたfactorを12桁interval演算へ再投入し，
        // 表示materialization後にも分解関係が十分保たれることを別に監査する。
        const std::size_t complexSize = 1 + rng() % 4;
        const auto complexMatrix = randomComplexMatrix(rng, complexSize);
        const std::array<mmcal::expression::Expr, 1> complexUnary{complexMatrix};
        const auto complexSvd = mmcal::builtins::evaluateApproximateSingularValueDecomposition(
            complexUnary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!complexSvd || !complexSvd->isArray()
            || complexSvd->asArray().shape != std::vector<std::size_t>{3, complexSize, complexSize})
            return reportFailure("certified-complex-svd-shape", complexMatrix);
        const auto complexU = packedMatrixFactor(*complexSvd, 0);
        const auto complexS = packedMatrixFactor(*complexSvd, 1);
        const auto complexV = packedMatrixFactor(*complexSvd, 2);
        const std::array<mmcal::expression::Expr, 2> complexUsArguments{complexU, complexS};
        const auto complexUs = mmcal::builtins::evaluateApproximateDot(
            complexUsArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{12});
        if (!complexUs)
            return reportFailure("certified-complex-svd-us-product", complexMatrix);
        const std::array<mmcal::expression::Expr, 1> complexVUnary{complexV};
        const auto complexVh = mmcal::builtins::evaluateConjugateTranspose(
            complexVUnary, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> complexReconstructArguments{*complexUs, complexVh};
        const auto complexReconstructed = mmcal::builtins::evaluateApproximateDot(
            complexReconstructArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{12});
        if (!complexReconstructed || !approximateMatrixContainsExact(*complexReconstructed, complexMatrix))
            return reportFailure("certified-complex-svd-reconstruction", complexMatrix);

        std::vector<mmcal::expression::Expr> complexIdentityElements;
        complexIdentityElements.reserve(complexSize * complexSize);
        for (std::size_t row = 0; row < complexSize; ++row)
            for (std::size_t column = 0; column < complexSize; ++column)
                complexIdentityElements.emplace_back(
                    mmcal::numeric::Number{BigInt{row == column ? 1 : 0}});
        const auto complexIdentity = mmcal::expression::Expr::array(
            {complexSize, complexSize}, std::move(complexIdentityElements));
        const std::array<mmcal::expression::Expr, 1> complexUUnary{complexU};
        const auto complexUh = mmcal::builtins::evaluateConjugateTranspose(
            complexUUnary, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> complexUOrthogonalArguments{complexUh, complexU};
        const auto complexUhu = mmcal::builtins::evaluateApproximateDot(
            complexUOrthogonalArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{12});
        if (!complexUhu || !approximateMatrixContainsExact(*complexUhu, complexIdentity))
            return reportFailure("certified-complex-svd-u-orthogonality", complexMatrix);
        const std::array<mmcal::expression::Expr, 2> complexVOrthogonalArguments{complexVh, complexV};
        const auto complexVhv = mmcal::builtins::evaluateApproximateDot(
            complexVOrthogonalArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{12});
        if (!complexVhv || !approximateMatrixContainsExact(*complexVhv, complexIdentity))
            return reportFailure("certified-complex-svd-v-orthogonality", complexMatrix);

        const auto eigenMatrix = randomSimpleSpectrumMatrix(eigenRng, size);
        const std::array<mmcal::expression::Expr, 1> eigenUnary{eigenMatrix};
        const auto eigenSystem = mmcal::builtins::evaluateApproximateEigensystem(
            eigenUnary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!eigenSystem || !eigenSystem->isList() || eigenSystem->asList().elements.size() != 2)
            return reportFailure("certified-eigensystem-shape", eigenMatrix);
        const auto& eigenValues = eigenSystem->asList().elements[0];
        const auto& eigenVectors = eigenSystem->asList().elements[1];
        if (!eigenValues.isArray() || eigenValues.asArray().shape != std::vector<std::size_t>{size}
            || !eigenVectors.isArray()
            || eigenVectors.asArray().shape != std::vector<std::size_t>{size, size})
            return reportFailure("certified-eigensystem-component-shape", eigenMatrix);
        if (!verifyDisplayedEigenRelation(eigenMatrix, eigenValues, eigenVectors))
            return reportFailure("certified-eigen-relation", eigenMatrix);

        if (!isExactIdentity(mmcal::builtins::evaluateRref(
                unary, fixture.registry, fixture.mathematics, fixture.angles)))
            return fail("exact-rref");

        const auto rationalMatrix = rationallyScaledMatrix(matrix);
        const std::array<mmcal::expression::Expr, 1> rationalUnary{rationalMatrix};
        const auto rationalInverse = mmcal::builtins::evaluateMatrixInverse(
            rationalUnary, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> rationalProductArguments{
            rationalMatrix, rationalInverse};
        if (!isExactIdentity(mmcal::builtins::evaluateDot(
                rationalProductArguments, fixture.registry, fixture.mathematics, fixture.angles)))
            return fail("rational-inverse");
        if (!isExactIdentity(mmcal::builtins::evaluateRref(
                rationalUnary, fixture.registry, fixture.mathematics, fixture.angles)))
            return fail("rational-rref");

        const auto expectedSolution = integerSolutionVector(size);
        const std::array<mmcal::expression::Expr, 2> rhsArguments{matrix, expectedSolution};
        const auto rhs = mmcal::builtins::evaluateDot(
            rhsArguments, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> solveArguments{matrix, rhs};
        const auto solution = mmcal::builtins::evaluateSolveLinear(
            solveArguments, fixture.registry, fixture.mathematics, fixture.angles);
        if (!(solution == expectedSolution))
            return fail("exact-linear-solve");

        const std::array<mmcal::expression::Expr, 2> rationalRhsArguments{
            rationalMatrix, expectedSolution};
        const auto rationalRhs = mmcal::builtins::evaluateDot(
            rationalRhsArguments, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> rationalSolveArguments{
            rationalMatrix, rationalRhs};
        if (!(mmcal::builtins::evaluateSolveLinear(
                rationalSolveArguments, fixture.registry, fixture.mathematics, fixture.angles)
                == expectedSolution))
            return fail("rational-linear-solve");

        const auto dependentMatrix = matrixWithDependentColumn(matrix);
        const std::array<mmcal::expression::Expr, 1> dependentUnary{dependentMatrix};
        const auto nullSpace = mmcal::builtins::evaluateNullSpace(
            dependentUnary, fixture.registry, fixture.mathematics, fixture.angles);
        if (!nullSpace.isArray() || nullSpace.asArray().shape != std::vector<std::size_t>{1, size + 1}
            || !validatesNullSpace(dependentMatrix, nullSpace, fixture))
            return fail("exact-null-space");

        const auto rationalDependentMatrix = rationallyScaledMatrix(dependentMatrix);
        const std::array<mmcal::expression::Expr, 1> rationalDependentUnary{rationalDependentMatrix};
        const auto rationalNullSpace = mmcal::builtins::evaluateNullSpace(
            rationalDependentUnary, fixture.registry, fixture.mathematics, fixture.angles);
        if (!rationalNullSpace.isArray()
            || rationalNullSpace.asArray().shape != std::vector<std::size_t>{1, size + 1}
            || !validatesNullSpace(rationalDependentMatrix, rationalNullSpace, fixture))
            return fail("rational-null-space");

        const auto approximateSolution = mmcal::builtins::evaluateApproximateSolveLinear(
            solveArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!approximateSolution || !approximateSolution->isArray()
            || approximateSolution->asArray().size() != size)
            return fail("certified-linear-solve-shape");
        for (std::size_t i = 0; i < size; ++i)
            if (!approximateContainsInteger(approximateSolution->asArray().element(i),
                    Rational{BigInt{static_cast<std::int64_t>(i + 1)}}))
                return fail("certified-linear-solve-containment");

        const auto approximateDeterminant = mmcal::builtins::evaluateApproximateDeterminant(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!approximateDeterminant || !approximateContainsInteger(
                *approximateDeterminant, determinant.asNumber().asReal().toRational()))
            return fail("certified-determinant-containment");
    }
    return true;
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
        if (!roundTrip || roundTrip->asArray().size() != n)
            return false;

        for (std::size_t i = 0; i < n; ++i)
            if (!approximateContainsInteger(roundTrip->asArray().element(i), expected[i]))
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

[[nodiscard]] double benchmarkLargeMatrixOperation(
    std::string_view operation,
    std::size_t size,
    std::size_t digits) {
    MatrixFixture fixture;
    const mmcal::expression::Expr matrix = randomDecimalMatrix(size);
    const std::array<mmcal::expression::Expr, 1> unary{matrix};
    const std::array<mmcal::expression::Expr, 2> binary{matrix, matrix};
    const mmcal::expression::Expr rhs = rowSumVector(matrix); // x={1,...,1} の既知解。
    const std::array<mmcal::expression::Expr, 2> linearSystem{matrix, rhs};
    const auto start = Clock::now();
    std::size_t checksum = 0;

    if (operation == "transpose") {
        const auto result = mmcal::builtins::evaluateTranspose(unary, fixture.registry);
        checksum = result.asArray().size();
    }
    else if (operation == "trace") {
        const auto result = mmcal::builtins::evaluateTrace(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.isNumber() ? 1 : 2;
    }
    else if (operation == "ndot") {
        const auto result = mmcal::builtins::evaluateApproximateDot(
            binary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->asArray().size();
    }
    else if (operation == "det") {
        const auto result = mmcal::builtins::evaluateDeterminant(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.isNumber() ? result.asNumber().realPart().toRational().numerator().bitLength() : 1;
    }
    else if (operation == "ndet") {
        const auto result = mmcal::builtins::evaluateApproximateDeterminant(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = 1;
    }
    else if (operation == "ninv") {
        const auto result = mmcal::builtins::evaluateApproximateInverse(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->asArray().size();
    }
    else if (operation == "nrank") {
        const auto result = mmcal::builtins::evaluateApproximateMatrixRank(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = 1;
    }
    else if (operation == "nsolve") {
        const auto result = mmcal::builtins::evaluateApproximateSolveLinear(
            linearSystem, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->asArray().size();
    }
    else if (operation == "nnull") {
        const auto result = mmcal::builtins::evaluateApproximateNullSpace(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->asArray().size() + 1;
    }
    else if (operation == "rref") {
        const auto result = mmcal::builtins::evaluateRref(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.asArray().size();
    }
    else if (operation == "rank") {
        const auto result = mmcal::builtins::evaluateMatrixRank(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.isNumber() ? 1 : 2;
    }
    else if (operation == "lu") {
        const auto result = mmcal::builtins::evaluateLuDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.isArray() ? result.asArray().size() : result.asList().elements.size();
    }
    else if (operation == "nlu") {
        const auto result = mmcal::builtins::evaluateApproximateLuDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->isArray() ? result->asArray().size() : result->asList().elements.size();
    }
    else if (operation == "nqr") {
        const auto result = mmcal::builtins::evaluateApproximateQrDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->isArray() ? result->asArray().size() : result->asList().elements.size();
    }
    else if (operation == "nsvd") {
        const auto result = mmcal::builtins::evaluateApproximateSingularValueDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->isArray() ? result->asArray().size() : result->asList().elements.size();
    }
    else if (operation == "neigen") {
        const auto result = mmcal::builtins::evaluateApproximateEigenvalues(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->asArray().size();
    }
    else if (operation == "neigensystem") {
        const auto result = mmcal::builtins::evaluateApproximateEigensystem(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->isArray() ? result->asArray().size() : result->asList().elements.size();
    }
    else {
        throw std::invalid_argument("Unknown large matrix operation");
    }

    const auto end = Clock::now();
    if (checksum == 0) std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

void runLargeMatrixBenchmark(std::string_view operation, std::size_t size, std::size_t digits) {
    std::cout << "large random decimal matrix: op=" << operation
              << " size=" << size << " digits=" << digits << '\n';
    std::cout << benchmarkLargeMatrixOperation(operation, size, digits) << " ms\n";
}

void runBenchmarks(bool full) {
    std::cout << "BigInt multiply/division\n";
    constexpr std::array<std::size_t, 7> limbCases{64, 128, 256, 512, 1024, 2048, 4096};
    for (const std::size_t limbs : std::span{limbCases}.first(full ? limbCases.size() : 5)) {
        const int iterations = limbs <= 256 ? 200 : limbs <= 1024 ? 40 : 8;
        std::cout << std::setw(5) << limbs << " limbs  mul="
                  << std::fixed << std::setprecision(3)
                  << benchmarkMultiply(limbs, iterations) << " us  square="
                  << benchmarkSquare(limbs, iterations) << " us  div="
                  << benchmarkDivision(limbs, iterations) << " us\n";
    }

    std::cout << "factorial\n";
    constexpr std::array<std::uint64_t, 4> factorialCases{10000, 20000, 40000, 80000};
    for (const std::uint64_t n : std::span{factorialCases}.first(full ? factorialCases.size() : 3)) {
        std::cout << std::setw(6) << n << "!  "
                  << benchmarkFactorial(n, n <= 10000 ? 3 : 1) << " ms\n";
    }

    std::cout << "decimal conversion\n";
    constexpr std::array<std::uint64_t, 3> decimalCases{10000, 20000, 40000};
    for (const std::uint64_t n : std::span{decimalCases}.first(full ? decimalCases.size() : 2))
        benchmarkDecimalConversion(n);

    std::cout << "exact/certified Matrix\n";
    constexpr std::array<std::size_t, 4> exactDotCases{16, 32, 64, 96};
    for (const std::size_t size : std::span{exactDotCases}.first(full ? exactDotCases.size() : 3)) {
        const int iterations = size <= 32 ? 5 : 2;
        std::cout << std::setw(5) << size << "x" << size << " dot="
                  << benchmarkExactDot(size, iterations) << " ms\n";
    }
    constexpr std::array<std::size_t, 4> exactMatrixCases{8, 12, 16, 20};
    for (const std::size_t size : std::span{exactMatrixCases}.first(full ? exactMatrixCases.size() : 3)) {
        const int iterations = size <= 12 ? 3 : 1;
        std::cout << std::setw(5) << size << "x" << size << " det="
                  << benchmarkExactDeterminant(size, iterations) << " ms  rref="
                  << benchmarkExactRref(size, iterations) << " ms  solve="
                  << benchmarkExactSolveLinear(size, iterations) << " ms  null="
                  << benchmarkExactNullSpace(size, iterations) << " ms  LU="
                  << benchmarkExactLu(size, iterations) << " ms  N[det,16]="
                  << benchmarkApproximateDeterminant(size, 16, iterations) << " ms\n";
        std::cout << "      certified decomposition  N[LU,16]="
                  << benchmarkApproximateLu(size, 16, iterations) << " ms  N[QR,16]="
                  << benchmarkApproximateQr(size, 16, iterations) << " ms\n";
    }

    std::cout << "exact fraction-free QR\n";
    for (const std::size_t size : std::initializer_list<std::size_t>{2, 4, 8, 16})
        std::cout << std::setw(5) << size << "x" << size << " exact QR="
                  << benchmarkExactQr(size, 1) << " ms\n";

    std::cout << "certified Householder QR block sweep\n";
    constexpr std::array<std::size_t, 4> qrBlockCases{8, 16, 24, 32};
    for (const std::size_t size : std::span{qrBlockCases}.first(full ? qrBlockCases.size() : 3)) {
        const int iterations = size <= 16 ? 2 : 1;
        std::cout << std::setw(5) << size << "x" << size
                  << " b1=" << benchmarkApproximateQrBlock(size, 16, iterations, 1) << " ms"
                  << " b8=" << benchmarkApproximateQrBlock(size, 16, iterations, 8) << " ms"
                  << " b16=" << benchmarkApproximateQrBlock(size, 16, iterations, 16) << " ms"
                  << " b32=" << benchmarkApproximateQrBlock(size, 16, iterations, 32) << " ms\n";
    }

    std::cout << "certified reduced SVD\n";
    for (const std::size_t size : std::initializer_list<std::size_t>{4, 8, 12, 16}) {
        const int iterations = size <= 8 ? 2 : 1;
        std::cout << std::setw(5) << size << "x" << size << " N[SVD,16]="
                  << benchmarkApproximateSvd(size, 16, iterations) << " ms\n";
    }

    std::cout << "certified eigen / eigensystem\n";
    for (const std::size_t size : std::initializer_list<std::size_t>{4, 8, 12, 16}) {
        const int iterations = size <= 8 ? 2 : 1;
        const double eigenvaluesMilliseconds =
            benchmarkApproximateEigenvalues(size, 16, iterations);
        const double eigensystemMilliseconds =
            benchmarkApproximateEigensystem(size, 16, iterations);
        std::cout << std::setw(5) << size << "x" << size << " N[eigenvalues,16]="
                  << eigenvaluesMilliseconds << " ms  N[eigensystem,16]="
                  << eigensystemMilliseconds << " ms\n";
    }

    std::cout << "exact/certified FFT\n";
    constexpr std::array<std::size_t, 5> radixTwoFftCases{32, 64, 128, 256, 512};
    for (const std::size_t size : std::span{radixTwoFftCases}.first(
            full ? radixTwoFftCases.size() : 3)) {
        const int exactIterations = size <= 64 ? 3 : 1;
        const int approximateIterations = size <= 128 ? 5 : 2;
        std::cout << std::setw(5) << size << " points  exact="
                  << benchmarkExactFft(size, exactIterations) << " ms  N[...,16]="
                  << benchmarkApproximateFft(size, 16, approximateIterations) << " ms\n";
    }
    constexpr std::array<std::size_t, 4> directFftCases{65, 127, 257, 509};
    for (const std::size_t size : std::span{directFftCases}.first(
            full ? directFftCases.size() : 2))
        std::cout << std::setw(5) << size << " points  direct="
                  << benchmarkApproximateDft(size, 16, 2) << " ms  FFT="
                  << benchmarkApproximateFft(size, 16, 2) << " ms\n";

    std::cout << "certified exp/log\n";
    constexpr std::array<std::size_t, 6> precisionCases{100, 500, 1000, 2000, 5000, 10000};
    for (const std::size_t digits : std::span{precisionCases}.first(
            full ? precisionCases.size() : 4)) {
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

void runAlgebraicFieldBenchmark(std::size_t iterations) {
    if (iterations == 0)
        throw std::invalid_argument("Algebraic field benchmark iterations must be positive");

    constexpr std::string_view expression =
        "(root[{-2,0,1},2]+root[{-3,0,0,1},1])*(root[{-2,0,1},2]-root[{-3,0,0,1},1])";
    constexpr std::string_view expected = "root[{1, 12, -6, 1}, 1]";

    mmcal::kernel::KernelSession session;
    double firstMilliseconds = 0.0;
    double warmMilliseconds = 0.0;
    for (std::size_t i = 0; i < iterations; ++i) {
        const auto start = Clock::now();
        const auto result = session.evaluate(expression);
        const auto end = Clock::now();
        if (mmcal::formatting::formatExpr(result) != expected)
            std::abort();

        const double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
        if (i == 0)
            firstMilliseconds = elapsed;
        else
            warmMilliseconds += elapsed;
    }

    std::cout << "algebraic compositum reuse\n"
              << "  expression: " << expression << '\n'
              << "  first: " << firstMilliseconds << " ms\n";
    if (iterations > 1)
        std::cout << "  warm avg: " << warmMilliseconds / static_cast<double>(iterations - 1)
                  << " ms (" << (iterations - 1) << " runs)\n";

    // x^12-2はEisensteinで既約。高めのdegreeでsame-field inverse reuseを単独測定する。
    std::vector<Rational> polynomial(13);
    polynomial[0] = Rational{BigInt{-2}};
    polynomial[12] = Rational{BigInt{1}};
    const auto generator = mmcal::symbolic::AlgebraicNumber::create(
        polynomial, 2, mmcal::symbolic::AlgebraicRootDomain::Real);
    auto field = generator ? mmcal::symbolic::NumberFieldContext::create(*generator) : nullptr;
    if (!field)
        std::abort();

    std::vector<Rational> lhsCoordinates(field->degree());
    std::vector<Rational> denominatorCoordinates(field->degree());
    for (std::size_t i = 0; i < field->degree(); ++i) {
        lhsCoordinates[i] = Rational{BigInt::fromUnsigned(i + 1)};
        denominatorCoordinates[i] = Rational{
            BigInt{static_cast<std::int64_t>(i % 5) - 2}};
    }
    denominatorCoordinates[0] = Rational{BigInt{2}};
    auto lhs = mmcal::symbolic::AlgebraicElement::create(field, std::move(lhsCoordinates));
    auto denominator = mmcal::symbolic::AlgebraicElement::create(
        field, std::move(denominatorCoordinates));
    if (!lhs || !denominator)
        std::abort();

    const auto reciprocalStart = Clock::now();
    const auto reciprocal = field->reciprocal(denominator->coefficients());
    const auto reciprocalEnd = Clock::now();
    if (!reciprocal)
        std::abort();
    const auto product = field->multiply(denominator->coefficients(), *reciprocal);
    if (product.empty() || product[0] != Rational{BigInt{1}}
        || !std::all_of(product.begin() + 1, product.end(),
            [](const Rational& coefficient) { return coefficient.isZero(); }))
        std::abort();

    const std::size_t microIterations = std::max<std::size_t>(256, iterations * 64);
    std::size_t checksum = 0;
    const auto warmReciprocalStart = Clock::now();
    for (std::size_t i = 0; i < microIterations; ++i) {
        const auto value = field->reciprocal(denominator->coefficients());
        if (!value)
            std::abort();
        checksum += value->size();
    }
    const auto warmReciprocalEnd = Clock::now();

    const auto warmDivideStart = Clock::now();
    for (std::size_t i = 0; i < microIterations; ++i) {
        const auto value = lhs->divide(*denominator);
        if (!value)
            std::abort();
        checksum += value->coefficients().size();
    }
    const auto warmDivideEnd = Clock::now();
    if (checksum == 0)
        std::abort();

    const auto minimalPolynomialStart = Clock::now();
    const auto minimalPolynomial = denominator->minimalPolynomial();
    const auto minimalPolynomialEnd = Clock::now();
    if (!minimalPolynomial || minimalPolynomial->empty())
        std::abort();

    const std::size_t minimalPolynomialIterations = std::max<std::size_t>(256, iterations * 64);
    std::size_t minimalPolynomialChecksum = 0;
    const auto warmMinimalPolynomialStart = Clock::now();
    for (std::size_t i = 0; i < minimalPolynomialIterations; ++i) {
        const auto minimal = denominator->minimalPolynomial();
        if (!minimal || *minimal != *minimalPolynomial)
            std::abort();
        minimalPolynomialChecksum += minimal->size();
    }
    const auto warmMinimalPolynomialEnd = Clock::now();
    if (minimalPolynomialChecksum == 0)
        std::abort();

    std::cout << "algebraic reciprocal reuse (degree 12)\n"
              << "  first reciprocal: "
              << std::chrono::duration<double, std::micro>(
                     reciprocalEnd - reciprocalStart).count()
              << " us\n"
              << "  warm reciprocal avg: "
              << std::chrono::duration<double, std::micro>(
                     warmReciprocalEnd - warmReciprocalStart).count()
                     / static_cast<double>(microIterations)
              << " us (" << microIterations << " runs)\n"
              << "  warm divide avg: "
              << std::chrono::duration<double, std::micro>(
                     warmDivideEnd - warmDivideStart).count()
                     / static_cast<double>(microIterations)
              << " us (" << microIterations << " runs)\n"
              << "  first minimal polynomial: "
              << std::chrono::duration<double, std::micro>(
                     minimalPolynomialEnd - minimalPolynomialStart).count()
              << " us\n"
              << "  warm minimal polynomial avg: "
              << std::chrono::duration<double, std::micro>(
                     warmMinimalPolynomialEnd - warmMinimalPolynomialStart).count()
                     / static_cast<double>(minimalPolynomialIterations)
              << " us (" << minimalPolynomialIterations << " runs)\n";
}

void runCertifiedSpecialFunctionBenchmark(std::size_t iterations) {
    const Rational oneThird{BigInt{1}, BigInt{3}};
    const Rational twoThirds{BigInt{2}, BigInt{3}};
    const Rational oneQuarter{BigInt{1}, BigInt{4}};

    std::cout << "certified special functions\n";
    for (const std::size_t bits : {80U, 160U, 320U, 640U, 1280U}) {
        const auto betaInput = mmcal::approximation::RealInterval::fromRational(
            oneQuarter, bits + 64);

        double gammaMilliseconds = 0.0;
        for (std::size_t i = 0; i < iterations; ++i) {
            const auto start = Clock::now();
            const auto value = mmcal::approximation::encloseGammaRational(oneThird, bits);
            const auto end = Clock::now();
            if (value.lower().toRational() <= Rational{BigInt{0}})
                std::abort();
            gammaMilliseconds += std::chrono::duration<double, std::milli>(end - start).count();
        }

        const std::size_t ibetaIterations = bits >= 640 ? 1 : iterations;
        double ibetaMilliseconds = 0.0;
        for (std::size_t i = 0; i < ibetaIterations; ++i) {
            const auto start = Clock::now();
            const auto value = mmcal::approximation::encloseIncompleteBetaRegularized(
                oneThird, twoThirds, betaInput, bits);
            const auto end = Clock::now();
            if (value.lower().toRational() < Rational{BigInt{0}}
                || value.upper().toRational() > Rational{BigInt{1}})
                std::abort();
            ibetaMilliseconds += std::chrono::duration<double, std::milli>(end - start).count();
        }

        std::cout << "  " << std::setw(4) << bits << " bit"
                  << "  gamma[1/3]=" << gammaMilliseconds / static_cast<double>(iterations) << " ms"
                  << "  ibeta[1/3,2/3,1/4]="
                  << ibetaMilliseconds / static_cast<double>(ibetaIterations) << " ms\n";
    }
}


[[nodiscard]] mmcal::linear_algebra::IntegerMatrixBuffer exactBackendMatrix(
    std::size_t size,
    std::size_t coefficientBits,
    unsigned salt) {
    mmcal::linear_algebra::IntegerMatrixBuffer matrix{size, size};
    BigInt diagonalBase{1};
    if (coefficientBits > 1)
        diagonalBase <<= coefficientBits - 1;

    const BigInt diagonalScale{static_cast<std::int64_t>(2 * size + 3)};
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            const std::int64_t perturbation = 1 + static_cast<std::int64_t>(
                (row * 37 + column * 53 + salt) % 127);
            if (row == column) {
                matrix(row, column) = diagonalBase * diagonalScale + BigInt{perturbation};
                continue;
            }

            BigInt value = diagonalBase + BigInt{perturbation};
            if (((row * 11 + column * 7 + salt) & 1U) != 0)
                value = -value;
            matrix(row, column) = std::move(value);
        }
    }
    return matrix;
}

[[nodiscard]] mmcal::expression::Expr exactBackendExpr(
    const mmcal::linear_algebra::IntegerMatrixBuffer& matrix) {
    std::vector<mmcal::expression::Expr> elements;
    elements.reserve(matrix.size());
    for (const BigInt& value : matrix.elements())
        elements.emplace_back(mmcal::numeric::Number{value});
    return mmcal::expression::Expr::array(
        {matrix.rows(), matrix.columns()}, std::move(elements));
}

[[nodiscard]] double benchmarkExactPublicPath(
    std::string_view operation,
    const mmcal::expression::Expr& input) {
    MatrixFixture fixture;
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    if (operation == "inverse") {
        const auto result = mmcal::builtins::evaluateMatrixInverse(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.asArray().size();
    } else if (operation == "rref") {
        const auto result = mmcal::builtins::evaluateRref(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.asArray().size();
    } else if (operation == "rank") {
        const auto result = mmcal::builtins::evaluateMatrixRank(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.asNumber().asReal().toRational().numerator().bitLength();
    } else if (operation == "null") {
        const auto result = mmcal::builtins::evaluateNullSpace(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.asArray().size();
    } else if (operation == "qr") {
        const auto result = mmcal::builtins::evaluateQrDecomposition(
            arguments, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.asArray().size();
    } else {
        throw std::invalid_argument("Unknown exact public-path benchmark operation");
    }
    const auto end = Clock::now();
    if (checksum == 0 && operation != "null")
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

[[nodiscard]] mmcal::linear_algebra::IntegerMatrixBuffer exactBackendAugmented(
    const mmcal::linear_algebra::IntegerMatrixBuffer& matrix,
    bool inverseRhs) {
    const std::size_t size = matrix.rows();
    const std::size_t rhsColumns = inverseRhs ? size : 1;
    mmcal::linear_algebra::IntegerMatrixBuffer augmented{size, size + rhsColumns};
    for (std::size_t row = 0; row < size; ++row) {
        BigInt rowSum;
        for (std::size_t column = 0; column < size; ++column) {
            augmented(row, column) = matrix(row, column);
            rowSum += matrix(row, column);
        }
        if (inverseRhs)
            augmented(row, size + row) = BigInt{1};
        else
            augmented(row, size) = std::move(rowSum);
    }
    return augmented;
}

[[nodiscard]] double benchmarkBareissDeterminantBackend(
    const mmcal::linear_algebra::IntegerMatrixBuffer& matrix,
    std::size_t iterations) {
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (std::size_t i = 0; i < iterations; ++i)
        checksum += mmcal::linear_algebra::bareissDeterminant(matrix).bitLength();
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count()
        / static_cast<double>(iterations);
}

[[nodiscard]] double benchmarkModularDeterminantBackend(
    const mmcal::linear_algebra::IntegerMatrixBuffer& matrix,
    std::size_t iterations,
    std::size_t& primes) {
    std::size_t checksum = 0;
    primes = 0;
    const auto start = Clock::now();
    for (std::size_t i = 0; i < iterations; ++i) {
        mmcal::linear_algebra::ModularLinearAlgebraStats stats;
        checksum += mmcal::linear_algebra::modularDeterminant(matrix, &stats).bitLength();
        primes += stats.primesAccepted;
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    primes /= iterations;
    return std::chrono::duration<double, std::milli>(end - start).count()
        / static_cast<double>(iterations);
}

[[nodiscard]] double benchmarkBareissAugmentedCore(
    const mmcal::linear_algebra::IntegerMatrixBuffer& augmented,
    std::size_t variables,
    std::size_t iterations) {
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (std::size_t i = 0; i < iterations; ++i)
        checksum += mmcal::linear_algebra::bareissEchelon(augmented, variables).pivotColumns.size();
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    return std::chrono::duration<double, std::milli>(end - start).count()
        / static_cast<double>(iterations);
}

[[nodiscard]] double benchmarkModularSolveBackend(
    const mmcal::linear_algebra::IntegerMatrixBuffer& augmented,
    std::size_t variables,
    std::size_t iterations,
    std::size_t& primes) {
    std::size_t checksum = 0;
    primes = 0;
    const auto start = Clock::now();
    for (std::size_t i = 0; i < iterations; ++i) {
        mmcal::linear_algebra::ModularLinearAlgebraStats stats;
        const auto solution = mmcal::linear_algebra::modularSolve(augmented, variables, &stats);
        if (!solution)
            std::abort();
        checksum += solution->size();
        primes += stats.primesAccepted;
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    primes /= iterations;
    return std::chrono::duration<double, std::milli>(end - start).count()
        / static_cast<double>(iterations);
}

[[nodiscard]] double benchmarkModularInverseBackend(
    const mmcal::linear_algebra::IntegerMatrixBuffer& matrix,
    std::size_t iterations,
    std::size_t& primes) {
    std::size_t checksum = 0;
    primes = 0;
    const auto start = Clock::now();
    for (std::size_t i = 0; i < iterations; ++i) {
        mmcal::linear_algebra::ModularLinearAlgebraStats stats;
        const auto inverse = mmcal::linear_algebra::modularInverse(matrix, &stats);
        if (!inverse)
            std::abort();
        checksum += inverse->size();
        primes += stats.primesAccepted;
    }
    const auto end = Clock::now();
    if (checksum == 0)
        std::abort();
    primes /= iterations;
    return std::chrono::duration<double, std::milli>(end - start).count()
        / static_cast<double>(iterations);
}

void runExactLinearAlgebraBackendBenchmark(std::size_t iterations) {
    if (iterations == 0)
        throw std::invalid_argument("Exact linear algebra benchmark iterations must be positive");

    std::cout << "exact integer linear algebra backend crossover\n";
    std::cout << "determinant / solveLinear\n";
    for (const std::size_t coefficientBits : {16U, 96U, 256U, 512U}) {
        std::cout << "coefficient height: " << coefficientBits << " bits\n";
        for (const std::size_t size : {4U, 6U, 8U, 10U, 12U, 16U, 20U, 24U, 32U}) {
            const auto matrix = exactBackendMatrix(size, coefficientBits, 31);
            const auto solveAugmented = exactBackendAugmented(matrix, false);
            const std::size_t localIterations = size <= 8 ? iterations : 1;

            std::size_t detPrimes = 0;
            std::size_t solvePrimes = 0;
            const double bareissDet = benchmarkBareissDeterminantBackend(
                matrix, localIterations);
            const double modularDet = benchmarkModularDeterminantBackend(
                matrix, localIterations, detPrimes);
            const double bareissSolve = benchmarkBareissAugmentedCore(
                solveAugmented, size, localIterations);
            const double modularSolve = benchmarkModularSolveBackend(
                solveAugmented, size, localIterations, solvePrimes);

            std::cout << std::setw(4) << size << "x" << size
                      << " det[B/M]=" << bareissDet << '/' << modularDet << " ms"
                      << " p=" << detPrimes
                      << " auto=" << (mmcal::linear_algebra::preferModularDeterminant(matrix)
                          ? "M" : "B")
                      << "  solve-core[B/M]=" << bareissSolve << '/' << modularSolve << " ms"
                      << " p=" << solvePrimes
                      << " auto=" << (mmcal::linear_algebra::preferModularSolve(solveAugmented, size)
                          ? "M" : "B") << '\n';
        }
    }

    // inverseは32x32・256-bitまでcrossoverが観測されていないため別枠で測る。
    // 512-bit sweepまで毎回含めると開発用benchmark自体が過度に重くなる。
    std::cout << "inverse\n";
    for (const std::size_t coefficientBits : {16U, 96U, 256U}) {
        std::cout << "coefficient height: " << coefficientBits << " bits\n";
        for (const std::size_t size : {4U, 6U, 8U, 10U, 12U, 16U, 24U, 32U}) {
            const auto matrix = exactBackendMatrix(size, coefficientBits, 31);
            const auto inverseAugmented = exactBackendAugmented(matrix, true);
            const std::size_t localIterations = size <= 8 ? iterations : 1;

            std::size_t inversePrimes = 0;
            const double bareissInverse = benchmarkBareissAugmentedCore(
                inverseAugmented, size, localIterations);
            const double modularInverse = benchmarkModularInverseBackend(
                matrix, localIterations, inversePrimes);
            std::cout << std::setw(4) << size << "x" << size
                      << " inverse-core[B/M]=" << bareissInverse << '/' << modularInverse << " ms"
                      << " p=" << inversePrimes
                      << " auto=B" << '\n';
        }
    }

    std::cout << "public exact paths (including final Rational/radical materialization)\n";
    for (const std::size_t coefficientBits : {16U, 96U, 256U}) {
        std::cout << "coefficient height: " << coefficientBits << " bits\n";
        for (const std::size_t size : {8U, 16U, 24U, 32U}) {
            const auto matrix = exactBackendExpr(exactBackendMatrix(size, coefficientBits, 47));
            const auto dependent = matrixWithDependentColumn(matrix);
            std::cout << std::setw(4) << size << "x" << size
                      << " inverse=" << benchmarkExactPublicPath("inverse", matrix) << " ms"
                      << " rref=" << benchmarkExactPublicPath("rref", matrix) << " ms"
                      << " rank=" << benchmarkExactPublicPath("rank", matrix) << " ms"
                      << " null(dep)=" << benchmarkExactPublicPath("null", dependent) << " ms"
                      << " qr=" << benchmarkExactPublicPath("qr", matrix) << " ms\n";
        }
    }

    std::cout << "rank-deficient structural paths\n";
    for (const std::size_t coefficientBits : {16U, 96U, 256U}) {
        std::cout << "coefficient height: " << coefficientBits << " bits\n";
        for (const std::size_t size : {48U, 64U}) {
            if (coefficientBits >= 256 && size > 48)
                continue;
            const auto matrix = exactBackendExpr(exactBackendMatrix(size, coefficientBits, 59));
            const auto dependent = matrixWithDependentColumn(matrix);
            std::cout << std::setw(4) << size << "x" << (size + 1)
                      << " rref=" << benchmarkExactPublicPath("rref", dependent) << " ms"
                      << " rank=" << benchmarkExactPublicPath("rank", dependent) << " ms"
                      << " null=" << benchmarkExactPublicPath("null", dependent) << " ms\n";
        }
    }
}

[[nodiscard]] std::string denseShiftedIdentityExpression(std::size_t order) {
    std::string result{"{"};
    for (std::size_t row = 0; row < order; ++row) {
        if (row != 0)
            result += ',';
        result += '{';
        for (std::size_t column = 0; column < order; ++column) {
            if (column != 0)
                result += ',';
            result += row == column ? "2" : "1";
        }
        result += '}';
    }
    result += '}';
    return result;
}

[[nodiscard]] std::string repeatedIntegerVectorExpression(
    std::size_t count,
    std::size_t value) {
    std::string result{"{"};
    for (std::size_t index = 0; index < count; ++index) {
        if (index != 0)
            result += ',';
        result += std::to_string(value);
    }
    result += '}';
    return result;
}

void printEvaluationUsage(
    std::string_view label,
    std::string_view expression,
    mmcal::kernel::KernelSession& session) {
    const auto start = Clock::now();
    static_cast<void>(session.evaluate(expression));
    const auto end = Clock::now();
    const auto& usage = session.lastEvaluationUsage();
    std::cout << label << "  "
              << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n"
              << "  input=" << usage.inputBytes
              << " work=" << usage.evaluationSteps
              << " depth=" << usage.maximumDepth
              << " nodes=" << usage.generatedNodes
              << " simplify=" << usage.simplificationCandidates
              << " solver=" << usage.solverBranches
              << " integrate=" << usage.integrationCandidates
              << " certified=" << usage.certifiedRefinements << '\n'
              << "  dense=" << usage.denseArrayElements
              << " matrix-temp=" << usage.temporaryMatrixElements
              << " bigint-max=" << usage.maximumBigIntegerBits
              << " precision-max=" << usage.maximumRequestedPrecisionDigits
              << " algebraic-degree-max=" << usage.maximumAlgebraicDegree
              << " algebraic-refine=" << usage.algebraicRefinements
              << " modular-primes=" << usage.modularPrimes << '\n';
}

void runEvaluationBudgetTelemetry() {
    mmcal::kernel::KernelSession session;
    std::cout << "EvaluationBudget representative workload telemetry\n";
    printEvaluationUsage("algebra", "factor[expand[(x+1)^12]]", session);
    const std::string dense48 = denseShiftedIdentityExpression(48);
    printEvaluationUsage("exact-det-modular", "det[" + dense48 + "]", session);
    const std::string dense32 = denseShiftedIdentityExpression(32);
    printEvaluationUsage("exact-solve-modular",
        "solveLinear[" + dense32 + "," + repeatedIntegerVectorExpression(32, 33) + "]",
        session);
    printEvaluationUsage("certified-exp", "N[exp[1],100]", session);
    printEvaluationUsage("integration", "integrate[sin[x]^2,x]", session);
}


void printUsage() {
    std::cout
        << "mmCal.Benchmarks [--full] [--random-only] [--benchmark-only]\n"
        << "  default          quick random checks + quick benchmarks\n"
        << "  --full           larger random set and high-precision benchmarks\n"
        << "  --random-only    run invariant checks only\n"
        << "  --benchmark-only run timings only\n"
        << "  --matrix-large <op> <size> [digits]\n"
        << "  --algebraic-field [iterations]  benchmark persistent compositum/embedding reuse\n"
        << "  --special-functions [iterations] benchmark certified gamma/ibeta scaling\n"
        << "  --exact-cyclotomic-fft [iterations] benchmark exact non-power-of-two FFT round trips\n"
        << "  --fft-threshold [iterations] compare certified direct DFT with forced Bluestein\n"
        << "  --exact-linear-algebra [iterations] compare Bareiss and modular exact backends\n"
        << "  --budget-telemetry collect representative EvaluationBudget usage\n"
        << "  --performance-cliffs [iterations] audit special-function bounded-work/performance cliffs\n"
        << "  --algebraic-root-cliffs [iterations] audit Complex Root isolation/refinement cliffs\n"
        << "    op: transpose trace ndot det ndet ninv rref rank nrank nsolve nnull lu nlu nqr nsvd neigen neigensystem\n"
        << "  --random-expressions [--loop|--nostop-loop] [--threads N] [--seed N] [--case N] [--cases N] [--max-depth N] [--report-every N]\n"
        << "    grammar-aware semantic fuzzer; --loop stops on the first FAIL, --nostop-loop reports FAILs and continues\n"
        << "  --certification-boundaries [--loop|--nostop-loop] [--threads N] [--seed N] [--case N] [--cases N] [--report-every N] [--timeout-ms N]\n"
        << "    certification fuzzer; boundary classification plus N/provenance metamorphic invariants\n";
}

} // namespace

int runBenchmarkMain(int argc, char** argv) {
    bool full = false;
    bool randomOnly = false;
    bool benchmarkOnly = false;
    std::optional<std::string> largeMatrixOperation;
    std::size_t largeMatrixSize = 0;
    std::size_t largeMatrixDigits = 16;
    bool algebraicFieldBenchmark = false;
    std::size_t algebraicFieldIterations = 8;
    bool specialFunctionBenchmark = false;
    std::size_t specialFunctionIterations = 1;
    bool performanceCliffAudit = false;
    std::size_t performanceCliffIterations = 1;
    bool algebraicRootCliffAudit = false;
    std::size_t algebraicRootCliffIterations = 1;
    bool exactCyclotomicFftBenchmark = false;
    std::size_t exactCyclotomicFftIterations = 3;
    bool approximateFftThresholdBenchmark = false;
    std::size_t approximateFftThresholdIterations = 1;
    bool exactLinearAlgebraBenchmark = false;
    std::size_t exactLinearAlgebraIterations = 3;
    bool budgetTelemetry = false;
    bool randomExpressions = false;
    bool certificationBoundaries = false;
    mmcal::benchmarks::RandomExpressionFuzzerOptions expressionOptions;
    mmcal::benchmarks::CertificationBoundaryFuzzerOptions boundaryOptions;
    const std::uint64_t defaultFuzzSeed = static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count())
        ^ (static_cast<std::uint64_t>(std::random_device{}()) << 32)
        ^ static_cast<std::uint64_t>(std::random_device{}());
    expressionOptions.seed = defaultFuzzSeed;
    boundaryOptions.seed = defaultFuzzSeed;

    // 共通の--seed/--case等を引数順に依存せずrouteするため，modeだけ先に判定する。
    for (int i = 1; i < argc; ++i) {
        if (std::string_view{argv[i]} == "--certification-boundaries")
            certificationBoundaries = true;
    }

    for (int i = 1; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--full")
            full = true;
        else if (arg == "--random-only")
            randomOnly = true;
        else if (arg == "--benchmark-only")
            benchmarkOnly = true;
        else if (arg == "--matrix-large") {
            if (i + 2 >= argc) {
                std::cerr << "--matrix-large requires <op> <size> [digits]\n";
                return 2;
            }
            largeMatrixOperation = std::string{argv[++i]};
            largeMatrixSize = static_cast<std::size_t>(std::stoull(argv[++i]));
            if (i + 1 < argc && argv[i + 1][0] != '-')
                largeMatrixDigits = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--algebraic-field") {
            algebraicFieldBenchmark = true;
            if (i + 1 < argc && argv[i + 1][0] != '-')
                algebraicFieldIterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--special-functions") {
            specialFunctionBenchmark = true;
            if (i + 1 < argc && argv[i + 1][0] != '-')
                specialFunctionIterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--performance-cliffs") {
            performanceCliffAudit = true;
            if (i + 1 < argc && argv[i + 1][0] != '-')
                performanceCliffIterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--algebraic-root-cliffs") {
            algebraicRootCliffAudit = true;
            if (i + 1 < argc && argv[i + 1][0] != '-')
                algebraicRootCliffIterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--exact-cyclotomic-fft") {
            exactCyclotomicFftBenchmark = true;
            if (i + 1 < argc && argv[i + 1][0] != '-')
                exactCyclotomicFftIterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--fft-threshold") {
            approximateFftThresholdBenchmark = true;
            if (i + 1 < argc && argv[i + 1][0] != '-')
                approximateFftThresholdIterations =
                    static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--exact-linear-algebra") {
            exactLinearAlgebraBenchmark = true;
            if (i + 1 < argc && argv[i + 1][0] != '-')
                exactLinearAlgebraIterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--budget-telemetry")
            budgetTelemetry = true;
        else if (arg == "--random-expressions")
            randomExpressions = true;
        else if (arg == "--certification-boundaries")
            certificationBoundaries = true;
        else if (arg == "--loop") {
            if (certificationBoundaries)
                boundaryOptions.loop = true;
            else {
                expressionOptions.loop = true;
                randomExpressions = true;
            }
        }
        else if (arg == "--nostop-loop") {
            if (certificationBoundaries) {
                boundaryOptions.loop = true;
                boundaryOptions.noStopLoop = true;
            }
            else {
                expressionOptions.loop = true;
                expressionOptions.noStopLoop = true;
                randomExpressions = true;
            }
        }
        else if (arg == "--threads") {
            if (i + 1 >= argc) {
                std::cerr << "--threads requires N\n";
                return 2;
            }
            const auto value = static_cast<std::size_t>(std::stoull(argv[++i]));
            if (certificationBoundaries)
                boundaryOptions.threads = value;
            else {
                expressionOptions.threads = value;
                randomExpressions = true;
            }
        }
        else if (arg == "--seed") {
            if (i + 1 >= argc) {
                std::cerr << "--seed requires N\n";
                return 2;
            }
            const auto value = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            if (certificationBoundaries)
                boundaryOptions.seed = value;
            else {
                expressionOptions.seed = value;
                randomExpressions = true;
            }
        }
        else if (arg == "--case") {
            if (i + 1 >= argc) {
                std::cerr << "--case requires N\n";
                return 2;
            }
            const auto value = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            if (certificationBoundaries)
                boundaryOptions.singleCase = value;
            else {
                expressionOptions.singleCase = value;
                randomExpressions = true;
            }
        }
        else if (arg == "--cases") {
            if (i + 1 >= argc) {
                std::cerr << "--cases requires N\n";
                return 2;
            }
            const auto value = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            if (certificationBoundaries)
                boundaryOptions.cases = value;
            else {
                expressionOptions.cases = value;
                randomExpressions = true;
            }
        }
        else if (arg == "--timeout-ms") {
            if (!certificationBoundaries) {
                std::cerr << "--timeout-ms is only valid with --certification-boundaries\n";
                return 2;
            }
            if (i + 1 >= argc) {
                std::cerr << "--timeout-ms requires N\n";
                return 2;
            }
            boundaryOptions.caseTimeoutMilliseconds =
                static_cast<std::uint64_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--max-depth") {
            if (certificationBoundaries) {
                std::cerr << "--max-depth is only valid with --random-expressions\n";
                return 2;
            }
            randomExpressions = true;
            if (i + 1 >= argc) {
                std::cerr << "--max-depth requires N\n";
                return 2;
            }
            expressionOptions.maxDepth = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--report-every") {
            if (i + 1 >= argc) {
                std::cerr << "--report-every requires N\n";
                return 2;
            }
            const auto value = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            if (certificationBoundaries)
                boundaryOptions.reportEvery = value;
            else {
                expressionOptions.reportEvery = value;
                randomExpressions = true;
            }
        }
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
    if (randomExpressions && certificationBoundaries) {
        std::cerr << "--random-expressions and --certification-boundaries cannot be combined\n";
        return 2;
    }

    if (largeMatrixOperation) {
        runLargeMatrixBenchmark(*largeMatrixOperation, largeMatrixSize, largeMatrixDigits);
        return 0;
    }

    if (algebraicFieldBenchmark) {
        if (algebraicFieldIterations == 0) {
            std::cerr << "--algebraic-field iterations must be at least 1\n";
            return 2;
        }
        runAlgebraicFieldBenchmark(algebraicFieldIterations);
        return 0;
    }

    if (specialFunctionBenchmark) {
        if (specialFunctionIterations == 0) {
            std::cerr << "--special-functions iterations must be at least 1\n";
            return 2;
        }
        runCertifiedSpecialFunctionBenchmark(specialFunctionIterations);
        return 0;
    }

    if (performanceCliffAudit) {
        if (performanceCliffIterations == 0) {
            std::cerr << "--performance-cliffs iterations must be at least 1\n";
            return 2;
        }
        mmcal::benchmarks::runSpecialFunctionPerformanceCliffAudit(
            {.iterations = performanceCliffIterations});
        return 0;
    }

    if (algebraicRootCliffAudit) {
        if (algebraicRootCliffIterations == 0) {
            std::cerr << "--algebraic-root-cliffs iterations must be at least 1\n";
            return 2;
        }
        mmcal::benchmarks::runAlgebraicRootPerformanceCliffAudit(
            {.iterations = algebraicRootCliffIterations});
        return 0;
    }

    if (exactCyclotomicFftBenchmark) {
        if (exactCyclotomicFftIterations == 0) {
            std::cerr << "--exact-cyclotomic-fft iterations must be at least 1\n";
            return 2;
        }
        runExactCyclotomicFftBenchmark(exactCyclotomicFftIterations);
        return 0;
    }

    if (approximateFftThresholdBenchmark) {
        if (approximateFftThresholdIterations == 0) {
            std::cerr << "--fft-threshold iterations must be at least 1\n";
            return 2;
        }
        runApproximateFftThresholdBenchmark(approximateFftThresholdIterations);
        return 0;
    }

    if (exactLinearAlgebraBenchmark) {
        if (exactLinearAlgebraIterations == 0) {
            std::cerr << "--exact-linear-algebra iterations must be at least 1\n";
            return 2;
        }
        runExactLinearAlgebraBackendBenchmark(exactLinearAlgebraIterations);
        return 0;
    }

    if (budgetTelemetry) {
        runEvaluationBudgetTelemetry();
        return 0;
    }

    if (certificationBoundaries) {
        if (boundaryOptions.caseTimeoutMilliseconds == 0) {
            std::cerr << "--timeout-ms must be at least 1\n";
            return 2;
        }
        if (boundaryOptions.threads == 0) {
            std::cerr << "--threads must be at least 1\n";
            return 2;
        }
        if (boundaryOptions.noStopLoop && boundaryOptions.singleCase) {
            std::cerr << "--nostop-loop cannot be combined with --case\n";
            return 2;
        }
        if (boundaryOptions.singleCase && *boundaryOptions.singleCase == 0) {
            std::cerr << "--case is 1-based and must be at least 1\n";
            return 2;
        }
        if (!boundaryOptions.loop && !boundaryOptions.singleCase && boundaryOptions.cases == 0) {
            std::cerr << "--cases must be at least 1\n";
            return 2;
        }
        return mmcal::benchmarks::runCertificationBoundaryFuzzer(boundaryOptions) ? 0 : 1;
    }

    if (randomExpressions) {
        if (expressionOptions.maxDepth == 0) {
            std::cerr << "--max-depth must be at least 1\n";
            return 2;
        }
        if (expressionOptions.threads == 0) {
            std::cerr << "--threads must be at least 1\n";
            return 2;
        }
        if (expressionOptions.noStopLoop && expressionOptions.singleCase) {
            std::cerr << "--nostop-loop cannot be combined with --case\n";
            return 2;
        }
        if (expressionOptions.singleCase && *expressionOptions.singleCase == 0) {
            std::cerr << "--case is 1-based and must be at least 1\n";
            return 2;
        }
        if (!expressionOptions.loop && !expressionOptions.singleCase && expressionOptions.cases == 0) {
            std::cerr << "--cases must be at least 1\n";
            return 2;
        }
        return mmcal::benchmarks::runRandomExpressionFuzzer(expressionOptions) ? 0 : 1;
    }

    if (!benchmarkOnly) {
        const std::size_t bigIntCases = full ? 2000 : 300;
        const std::size_t certifiedCases = full ? 200 : 40;
        const std::size_t matrixCases = full ? 200 : 40;
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
        std::cout << "Random certified Matrix checks: " << matrixCases << " cases\n";
        if (!runRandomMatrixChecks(matrixCases)) {
            std::cerr << "Random certified Matrix check failed\n";
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

int main(int argc, char** argv) {
    try {
        return runBenchmarkMain(argc, argv);
    }
    catch (const BenchmarkFailure& error) {
        std::cerr << "Benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
