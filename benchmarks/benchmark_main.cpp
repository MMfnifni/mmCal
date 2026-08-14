// mmCalの性能測定・ランダム不変量試験を本体testから分離して実行する開発用runner
#include "approximation/certified_exponential.hpp"
#include "approximation/certified_logarithm.hpp"
#include "approximation/certified_special_functions.hpp"
#include "approximation/precision.hpp"
#include "approximation/real_interval.hpp"
#include "builtins/signal_processing.hpp"
#include "builtins/linear_algebra.hpp"
#include "linear_algebra/decomposition.hpp"
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
#include <string>
#include <string_view>
#include <stdexcept>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using mmcal::numeric::BigInt;
using mmcal::numeric::Rational;
using mmcal::numeric::RealNumber;

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
    std::vector<mmcal::expression::Expr> values;
    values.reserve(size * size);
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            const std::int64_t value = row == column
                ? static_cast<std::int64_t>(4 * size + 1 + salt % 3)
                : static_cast<std::int64_t>((row * 17 + column * 29 + salt) % 5) - 2;
            values.emplace_back(mmcal::numeric::Number{BigInt{value}});
        }
    }
    return mmcal::expression::Expr::array({size, size}, std::move(values));
}

[[nodiscard]] mmcal::expression::Expr randomDecimalMatrix(
    std::size_t size,
    std::uint64_t seed = 1234) {
    // 添付fft_gen.pyと同じ [-1,1] / 小数10桁という負荷特性を，
    // parser時間と算法時間を分離するためexact Rationalとして直接構築する。
    constexpr std::int64_t scale = 10000000000LL;
    std::mt19937_64 rng{seed};
    std::uniform_int_distribution<std::int64_t> distribution{-scale, scale};
    std::vector<mmcal::expression::Expr> values;
    values.reserve(size * size);
    for (std::size_t i = 0; i < size * size; ++i)
        values.emplace_back(mmcal::numeric::Number{Rational{
            BigInt{distribution(rng)}, BigInt{scale}}});
    return mmcal::expression::Expr::array({size, size}, std::move(values));
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
            sum += array.elements[row * size + column].asNumber().realPart().toRational();
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
        checksum += result.asArray().elements.size();
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
        checksum += result.asArray().elements.size();
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
        checksum += result.asArray().elements.size();
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
            values.push_back(array.elements[row * columns + column]);

        const auto first = array.elements[row * columns].asNumber();
        const auto second = columns > 1
            ? array.elements[row * columns + 1].asNumber()
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
        checksum += result.asArray().elements.size();
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
        checksum += result.asArray().elements.size();
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
        checksum += result.asArray().elements.size();
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
        checksum += result->asArray().elements.size();
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
        checksum += result->asArray().elements.size();
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
            ? result->asArray().elements.size()
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
    const mmcal::expression::Expr input = matrixInput(size, 29);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateEigenvalues(
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

[[nodiscard]] double benchmarkApproximateEigensystem(
    std::size_t size, std::size_t digits, int iterations) {
    MatrixFixture fixture;
    const mmcal::expression::Expr input = matrixInput(size, 29);
    const std::array<mmcal::expression::Expr, 1> arguments{input};
    std::size_t checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto result = mmcal::builtins::evaluateApproximateEigensystem(
            arguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result)
            std::abort();
        checksum += result->isArray() ? result->asArray().elements.size() : result->asList().elements.size();
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
        checksum += result->asArray().elements.size();
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

[[nodiscard]] mmcal::expression::Expr rationallyScaledMatrix(
    const mmcal::expression::Expr& matrix) {
    const auto& array = matrix.asArray();
    const std::size_t rows = array.shape[0];
    const std::size_t columns = array.shape[1];
    std::vector<mmcal::expression::Expr> values;
    values.reserve(array.elements.size());

    for (std::size_t row = 0; row < rows; ++row)
        for (std::size_t column = 0; column < columns; ++column) {
            const auto rational = array.elements[row * columns + column]
                .asNumber().asReal().toRational();
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
            const auto& value = matrix.asArray().elements[row * size + column];
            if (!value.isNumber())
                return false;
            const mmcal::numeric::Number expected{BigInt{row == column ? 1 : 0}};
            if (!(value.asNumber() == expected))
                return false;
        }
    return true;
}

[[nodiscard]] bool isExactZeroVector(const mmcal::expression::Expr& vector) {
    if (!vector.isArray() || !vector.asArray().isVector())
        return false;
    for (const auto& value : vector.asArray().elements)
        if (!value.isNumber() || !value.asNumber().isZero())
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
            vectorElements.push_back(basisArray.elements[row * variables + column]);
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
    const auto begin = array.elements.begin() + static_cast<std::ptrdiff_t>(factor * factorSize);
    std::vector<mmcal::expression::Expr> values(begin, begin + static_cast<std::ptrdiff_t>(factorSize));
    return mmcal::expression::Expr::array({rows, columns}, std::move(values));
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
    const auto& lhs = approximate.asArray().elements;
    const auto& rhs = exact.asArray().elements;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (!rhs[i].isNumber() || !approximateContainsNumber(lhs[i], rhs[i].asNumber()))
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
        const auto lambda = displayedComplex(lambdas.elements[column]);
        if (!lambda)
            return false;
        for (std::size_t row = 0; row < size; ++row) {
            DisplayedComplex av{Rational{BigInt{0}}, Rational{BigInt{0}}};
            for (std::size_t k = 0; k < size; ++k) {
                const auto aik = displayedComplex(a.elements[row * size + k]);
                const auto vk = displayedComplex(v.elements[k * size + column]);
                if (!aik || !vk)
                    return false;
                const DisplayedComplex product = multiplyDisplayed(*aik, *vk);
                av.real += product.real;
                av.imaginary += product.imaginary;
            }
            const auto vr = displayedComplex(v.elements[row * size + column]);
            if (!vr)
                return false;
            const DisplayedComplex lv = multiplyDisplayed(*lambda, *vr);
            if (absoluteRational(av.real - lv.real) > tolerance
                || absoluteRational(av.imaginary - lv.imaginary) > tolerance)
                return false;
        }
    }
    return true;
}

[[nodiscard]] bool runRandomMatrixChecks(std::size_t count) {
    MatrixFixture fixture;
    std::mt19937_64 rng{0x4D41545249584345ULL};

    for (std::size_t caseIndex = 0; caseIndex < count; ++caseIndex) {
        const std::size_t size = 1 + rng() % 5;
        const auto matrix = randomInvertibleMatrix(rng, size);
        const std::array<mmcal::expression::Expr, 1> unary{matrix};

        const auto determinant = mmcal::builtins::evaluateDeterminant(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        if (!determinant.isNumber() || !determinant.asNumber().isReal()
            || determinant.asNumber().isZero())
            return false;

        const auto transposed = mmcal::builtins::evaluateTranspose(unary, fixture.registry);
        const std::array<mmcal::expression::Expr, 1> transposedArguments{transposed};
        const auto transposedDeterminant = mmcal::builtins::evaluateDeterminant(
            transposedArguments, fixture.registry, fixture.mathematics, fixture.angles);
        if (!(transposedDeterminant == determinant))
            return false;

        const auto inverse = mmcal::builtins::evaluateMatrixInverse(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> productArguments{matrix, inverse};
        if (!isExactIdentity(mmcal::builtins::evaluateDot(
                productArguments, fixture.registry, fixture.mathematics, fixture.angles)))
            return false;

        const auto lu = mmcal::builtins::evaluateLuDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        if (!lu.isArray() || lu.asArray().shape != std::vector<std::size_t>{3, size, size})
            return false;
        const auto permutation = packedMatrixFactor(lu, 0);
        const auto lower = packedMatrixFactor(lu, 1);
        const auto upper = packedMatrixFactor(lu, 2);
        const std::array<mmcal::expression::Expr, 2> paArguments{permutation, matrix};
        const std::array<mmcal::expression::Expr, 2> luArguments{lower, upper};
        if (!(mmcal::builtins::evaluateDot(
                paArguments, fixture.registry, fixture.mathematics, fixture.angles)
            == mmcal::builtins::evaluateDot(
                luArguments, fixture.registry, fixture.mathematics, fixture.angles)))
            return false;

        const auto approximateQr = mmcal::builtins::evaluateApproximateQrDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!approximateQr || !approximateQr->isArray()
            || approximateQr->asArray().shape != std::vector<std::size_t>{2, size, size})
            return false;
        const auto q = packedMatrixFactor(*approximateQr, 0);
        const auto r = packedMatrixFactor(*approximateQr, 1);
        const std::array<mmcal::expression::Expr, 2> qrArguments{q, r};
        const auto reconstructed = mmcal::builtins::evaluateApproximateDot(
            qrArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{18});
        if (!reconstructed || !approximateMatrixContainsExact(*reconstructed, matrix))
            return false;

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
            return false;

        const auto approximateSvd = mmcal::builtins::evaluateApproximateSingularValueDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!approximateSvd || !approximateSvd->isArray()
            || approximateSvd->asArray().shape != std::vector<std::size_t>{3, size, size})
            return false;
        const auto svdU = packedMatrixFactor(*approximateSvd, 0);
        const auto svdS = packedMatrixFactor(*approximateSvd, 1);
        const auto svdV = packedMatrixFactor(*approximateSvd, 2);
        const std::array<mmcal::expression::Expr, 2> usArguments{svdU, svdS};
        const auto us = mmcal::builtins::evaluateApproximateDot(
            usArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{18});
        if (!us)
            return false;
        const std::array<mmcal::expression::Expr, 1> svdVUnary{svdV};
        const auto svdVt = mmcal::builtins::evaluateTranspose(svdVUnary, fixture.registry);
        const std::array<mmcal::expression::Expr, 2> svdReconstructArguments{*us, svdVt};
        const auto svdReconstructed = mmcal::builtins::evaluateApproximateDot(
            svdReconstructArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{18});
        if (!svdReconstructed || !approximateMatrixContainsExact(*svdReconstructed, matrix))
            return false;

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
            || complexSvd->asArray().shape != std::vector<std::size_t>{3, complexSize, complexSize}) {
            std::cerr << "complex SVD factor failure case=" << caseIndex << " size=" << complexSize << "\n";
            return false;
        }
        const auto complexU = packedMatrixFactor(*complexSvd, 0);
        const auto complexS = packedMatrixFactor(*complexSvd, 1);
        const auto complexV = packedMatrixFactor(*complexSvd, 2);
        const std::array<mmcal::expression::Expr, 2> complexUsArguments{complexU, complexS};
        const auto complexUs = mmcal::builtins::evaluateApproximateDot(
            complexUsArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{12});
        if (!complexUs)
            return false;
        const std::array<mmcal::expression::Expr, 1> complexVUnary{complexV};
        const auto complexVh = mmcal::builtins::evaluateConjugateTranspose(
            complexVUnary, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> complexReconstructArguments{*complexUs, complexVh};
        const auto complexReconstructed = mmcal::builtins::evaluateApproximateDot(
            complexReconstructArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{12});
        if (!complexReconstructed || !approximateMatrixContainsExact(*complexReconstructed, complexMatrix)) {
            std::cerr << "complex SVD reconstruction failure case=" << caseIndex << " size=" << complexSize << "\n";
            return false;
        }

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
        if (!complexUhu || !approximateMatrixContainsExact(*complexUhu, complexIdentity)) {
            std::cerr << "complex SVD U orthogonality failure case=" << caseIndex << " size=" << complexSize << "\n";
            return false;
        }
        const std::array<mmcal::expression::Expr, 2> complexVOrthogonalArguments{complexVh, complexV};
        const auto complexVhv = mmcal::builtins::evaluateApproximateDot(
            complexVOrthogonalArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{12});
        if (!complexVhv || !approximateMatrixContainsExact(*complexVhv, complexIdentity)) {
            std::cerr << "complex SVD V orthogonality failure case=" << caseIndex << " size=" << complexSize << "\n";
            return false;
        }

        const auto eigenSystem = mmcal::builtins::evaluateApproximateEigensystem(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!eigenSystem || !eigenSystem->isList() || eigenSystem->asList().elements.size() != 2)
            return false;
        const auto& eigenValues = eigenSystem->asList().elements[0];
        const auto& eigenVectors = eigenSystem->asList().elements[1];
        if (!eigenValues.isArray() || eigenValues.asArray().shape != std::vector<std::size_t>{size}
            || !eigenVectors.isArray()
            || eigenVectors.asArray().shape != std::vector<std::size_t>{size, size})
            return false;
        if (!verifyDisplayedEigenRelation(matrix, eigenValues, eigenVectors))
            return false;

        if (!isExactIdentity(mmcal::builtins::evaluateRref(
                unary, fixture.registry, fixture.mathematics, fixture.angles)))
            return false;

        const auto rationalMatrix = rationallyScaledMatrix(matrix);
        const std::array<mmcal::expression::Expr, 1> rationalUnary{rationalMatrix};
        const auto rationalInverse = mmcal::builtins::evaluateMatrixInverse(
            rationalUnary, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> rationalProductArguments{
            rationalMatrix, rationalInverse};
        if (!isExactIdentity(mmcal::builtins::evaluateDot(
                rationalProductArguments, fixture.registry, fixture.mathematics, fixture.angles)))
            return false;
        if (!isExactIdentity(mmcal::builtins::evaluateRref(
                rationalUnary, fixture.registry, fixture.mathematics, fixture.angles)))
            return false;

        const auto expectedSolution = integerSolutionVector(size);
        const std::array<mmcal::expression::Expr, 2> rhsArguments{matrix, expectedSolution};
        const auto rhs = mmcal::builtins::evaluateDot(
            rhsArguments, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> solveArguments{matrix, rhs};
        const auto solution = mmcal::builtins::evaluateSolveLinear(
            solveArguments, fixture.registry, fixture.mathematics, fixture.angles);
        if (!(solution == expectedSolution))
            return false;

        const std::array<mmcal::expression::Expr, 2> rationalRhsArguments{
            rationalMatrix, expectedSolution};
        const auto rationalRhs = mmcal::builtins::evaluateDot(
            rationalRhsArguments, fixture.registry, fixture.mathematics, fixture.angles);
        const std::array<mmcal::expression::Expr, 2> rationalSolveArguments{
            rationalMatrix, rationalRhs};
        if (!(mmcal::builtins::evaluateSolveLinear(
                rationalSolveArguments, fixture.registry, fixture.mathematics, fixture.angles)
                == expectedSolution))
            return false;

        const auto dependentMatrix = matrixWithDependentColumn(matrix);
        const std::array<mmcal::expression::Expr, 1> dependentUnary{dependentMatrix};
        const auto nullSpace = mmcal::builtins::evaluateNullSpace(
            dependentUnary, fixture.registry, fixture.mathematics, fixture.angles);
        if (!nullSpace.isArray() || nullSpace.asArray().shape != std::vector<std::size_t>{1, size + 1}
            || !validatesNullSpace(dependentMatrix, nullSpace, fixture))
            return false;

        const auto rationalDependentMatrix = rationallyScaledMatrix(dependentMatrix);
        const std::array<mmcal::expression::Expr, 1> rationalDependentUnary{rationalDependentMatrix};
        const auto rationalNullSpace = mmcal::builtins::evaluateNullSpace(
            rationalDependentUnary, fixture.registry, fixture.mathematics, fixture.angles);
        if (!rationalNullSpace.isArray()
            || rationalNullSpace.asArray().shape != std::vector<std::size_t>{1, size + 1}
            || !validatesNullSpace(rationalDependentMatrix, rationalNullSpace, fixture))
            return false;

        const auto approximateSolution = mmcal::builtins::evaluateApproximateSolveLinear(
            solveArguments, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!approximateSolution || !approximateSolution->isArray()
            || approximateSolution->asArray().elements.size() != size)
            return false;
        for (std::size_t i = 0; i < size; ++i)
            if (!approximateContainsInteger(approximateSolution->asArray().elements[i],
                    Rational{BigInt{static_cast<std::int64_t>(i + 1)}}))
                return false;

        const auto approximateDeterminant = mmcal::builtins::evaluateApproximateDeterminant(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{20});
        if (!approximateDeterminant || !approximateContainsInteger(
                *approximateDeterminant, determinant.asNumber().asReal().toRational()))
            return false;
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
        checksum = result.asArray().elements.size();
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
        checksum = result->asArray().elements.size();
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
        checksum = result->asArray().elements.size();
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
        checksum = result->asArray().elements.size();
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
        checksum = result.asArray().elements.size();
    }
    else if (operation == "rank") {
        const auto result = mmcal::builtins::evaluateMatrixRank(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.isNumber() ? 1 : 2;
    }
    else if (operation == "lu") {
        const auto result = mmcal::builtins::evaluateLuDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles);
        checksum = result.isArray() ? result.asArray().elements.size() : result.asList().elements.size();
    }
    else if (operation == "nlu") {
        const auto result = mmcal::builtins::evaluateApproximateLuDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->isArray() ? result->asArray().elements.size() : result->asList().elements.size();
    }
    else if (operation == "nqr") {
        const auto result = mmcal::builtins::evaluateApproximateQrDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->isArray() ? result->asArray().elements.size() : result->asList().elements.size();
    }
    else if (operation == "nsvd") {
        const auto result = mmcal::builtins::evaluateApproximateSingularValueDecomposition(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->isArray() ? result->asArray().elements.size() : result->asList().elements.size();
    }
    else if (operation == "neigen") {
        const auto result = mmcal::builtins::evaluateApproximateEigenvalues(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->asArray().elements.size();
    }
    else if (operation == "neigensystem") {
        const auto result = mmcal::builtins::evaluateApproximateEigensystem(
            unary, fixture.registry, fixture.mathematics, fixture.angles,
            mmcal::approximation::ApproximationContext{digits});
        if (!result) std::abort();
        checksum = result->isArray() ? result->asArray().elements.size() : result->asList().elements.size();
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

    std::cout << "exact/certified Matrix\n";
    for (const std::size_t size : full
        ? std::initializer_list<std::size_t>{16, 32, 64, 96}
        : std::initializer_list<std::size_t>{16, 32, 64}) {
        const int iterations = size <= 32 ? 5 : 2;
        std::cout << std::setw(5) << size << "x" << size << " dot="
                  << benchmarkExactDot(size, iterations) << " ms\n";
    }
    for (const std::size_t size : full
        ? std::initializer_list<std::size_t>{8, 12, 16, 20}
        : std::initializer_list<std::size_t>{8, 12, 16}) {
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

    std::cout << "exact Householder QR (small-order symbolic path)\n";
    for (const std::size_t size : std::initializer_list<std::size_t>{2, 3})
        std::cout << std::setw(5) << size << "x" << size << " exact QR="
                  << benchmarkExactQr(size, 1) << " ms\n";

    std::cout << "certified Householder QR block sweep\n";
    for (const std::size_t size : full
        ? std::initializer_list<std::size_t>{8, 16, 24, 32}
        : std::initializer_list<std::size_t>{8, 16, 24}) {
        const int iterations = size <= 16 ? 2 : 1;
        std::cout << std::setw(5) << size << "x" << size
                  << " b1=" << benchmarkApproximateQrBlock(size, 16, iterations, 1) << " ms"
                  << " b8=" << benchmarkApproximateQrBlock(size, 16, iterations, 8) << " ms"
                  << " b16=" << benchmarkApproximateQrBlock(size, 16, iterations, 16) << " ms"
                  << " b32=" << benchmarkApproximateQrBlock(size, 16, iterations, 32) << " ms\n";
    }

    std::cout << "certified reduced SVD\n";
    for (const std::size_t size : full
        ? std::initializer_list<std::size_t>{4, 8, 12, 16}
        : std::initializer_list<std::size_t>{4, 8, 12, 16}) {
        const int iterations = size <= 8 ? 2 : 1;
        std::cout << std::setw(5) << size << "x" << size << " N[SVD,16]="
                  << benchmarkApproximateSvd(size, 16, iterations) << " ms\n";
    }

    std::cout << "certified eigen / eigensystem\n";
    for (const std::size_t size : std::initializer_list<std::size_t>{4, 8, 12, 16}) {
        const int iterations = size <= 8 ? 2 : 1;
        std::cout << std::setw(5) << size << "x" << size << " N[eigenvalues,16]="
                  << benchmarkApproximateEigenvalues(size, 16, iterations) << " ms  N[eigensystem,16]="
                  << benchmarkApproximateEigensystem(size, 16, iterations) << " ms\n";
    }

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
        << "  --benchmark-only run timings only\n"
        << "  --matrix-large <op> <size> [digits]\n"
        << "    op: transpose trace ndot det ndet ninv rref rank nrank nsolve nnull lu nlu nqr nsvd neigen neigensystem\n"
        << "  --random-expressions [--loop|--nostop-loop] [--threads N] [--seed N] [--case N] [--cases N] [--max-depth N] [--report-every N]\n"
        << "    grammar-aware semantic fuzzer; --loop stops on the first FAIL, --nostop-loop reports FAILs and continues\n";
}

} // namespace

int main(int argc, char** argv) {
    bool full = false;
    bool randomOnly = false;
    bool benchmarkOnly = false;
    std::optional<std::string> largeMatrixOperation;
    std::size_t largeMatrixSize = 0;
    std::size_t largeMatrixDigits = 16;
    bool randomExpressions = false;
    mmcal::benchmarks::RandomExpressionFuzzerOptions expressionOptions;
    expressionOptions.seed = static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count())
        ^ (static_cast<std::uint64_t>(std::random_device{}()) << 32)
        ^ static_cast<std::uint64_t>(std::random_device{}());
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
        else if (arg == "--random-expressions")
            randomExpressions = true;
        else if (arg == "--loop") {
            expressionOptions.loop = true;
            randomExpressions = true;
        }
        else if (arg == "--nostop-loop") {
            expressionOptions.loop = true;
            expressionOptions.noStopLoop = true;
            randomExpressions = true;
        }
        else if (arg == "--threads") {
            randomExpressions = true;
            if (i + 1 >= argc) {
                std::cerr << "--threads requires N\n";
                return 2;
            }
            expressionOptions.threads = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--seed") {
            randomExpressions = true;
            if (i + 1 >= argc) {
                std::cerr << "--seed requires N\n";
                return 2;
            }
            expressionOptions.seed = static_cast<std::uint64_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--case") {
            if (i + 1 >= argc) {
                std::cerr << "--case requires N\n";
                return 2;
            }
            expressionOptions.singleCase = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            randomExpressions = true;
        }
        else if (arg == "--cases") {
            randomExpressions = true;
            if (i + 1 >= argc) {
                std::cerr << "--cases requires N\n";
                return 2;
            }
            expressionOptions.cases = static_cast<std::uint64_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--max-depth") {
            randomExpressions = true;
            if (i + 1 >= argc) {
                std::cerr << "--max-depth requires N\n";
                return 2;
            }
            expressionOptions.maxDepth = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--report-every") {
            randomExpressions = true;
            if (i + 1 >= argc) {
                std::cerr << "--report-every requires N\n";
                return 2;
            }
            expressionOptions.reportEvery = static_cast<std::uint64_t>(std::stoull(argv[++i]));
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

    if (largeMatrixOperation) {
        runLargeMatrixBenchmark(*largeMatrixOperation, largeMatrixSize, largeMatrixDigits);
        return 0;
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
