// combinatorics・signalの回帰テスト
#include "combinatorics_signal_tests.hpp"

#include "builtins/signal_processing.hpp"
#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "evaluation/builtin_registry.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

#include <array>
#include <stdexcept>
#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

std::string eval(kernel::KernelSession& session, std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

error::CalcError evalError(kernel::KernelSession& session, std::string_view source) {
    try { static_cast<void>(session.evaluate(source)); }
    catch (const error::CalcError& e) { return e; }
    throw std::logic_error("Expected CalcError");
}

} // namespace

void runCombinatoricsSignalTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(eval(session, "perm[10,3]"), std::string{"720"},
        "permutation is exact");
    tests.expectEqual(eval(session, "perm[5,7]"), std::string{"0"},
        "permutation is zero when r exceeds n");
    tests.expectEqual(eval(session, "comb[10,3]"), std::string{"120"},
        "combination is exact");
    tests.expectEqual(eval(session, "comb[100,50]"),
        std::string{"100891344545564193334812497256"},
        "combination uses arbitrary-precision integers");
    tests.expectEqual(eval(session, "fib[0]"), std::string{"0"},
        "Fibonacci zero is exact");
    tests.expectEqual(eval(session, "fib[25]"), std::string{"75025"},
        "Fibonacci fast doubling returns the expected value");
    tests.expectEqual(eval(session, "fib[100]"), std::string{"354224848179261915075"},
        "Fibonacci is not limited by machine integer width");
    tests.expectEqual(eval(session, "perm[n,2]"), std::string{"perm[n, 2]"},
        "symbolic combinatorics remain unevaluated");
    tests.expect(evalError(session, "comb[-1,2]").type() == error::CalcErrorType::Domain,
        "combination rejects negative n");
    tests.expect(evalError(session, "fib[-1]").type() == error::CalcErrorType::Domain,
        "Fibonacci rejects negative indices until a signed extension is specified");

    tests.expectEqual(eval(session, "isprime[97]"), std::string{"True"},
        "isprime proves a uint64 prime deterministically");
    tests.expectEqual(eval(session, "isprime[561]"), std::string{"False"},
        "isprime rejects a Carmichael composite");
    tests.expectEqual(eval(session, "isprime[18446744073709551557]"), std::string{"True"},
        "isprime covers the deterministic uint64 boundary region");
    tests.expectEqual(eval(session, "nextprime[14]"), std::string{"17"},
        "nextprime returns the next exact prime");
    tests.expectEqual(eval(session, "prevprime[14]"), std::string{"13"},
        "prevprime returns the previous exact prime");
    tests.expect(evalError(session, "prevprime[2]").type() == error::CalcErrorType::Domain,
        "prevprime rejects values without a positive predecessor prime");
    tests.expectEqual(eval(session, "factorint[360]"), std::string{"{2, 2, 2, 3, 3, 5}"},
        "factorint returns flattened sorted prime factors");
    tests.expectEqual(eval(session, "factorint[-84]"), std::string{"{-1, 2, 2, 3, 7}"},
        "factorint preserves the sign as a leading minus-one factor");
    tests.expectEqual(eval(session, "factorint[1000036000099]"),
        std::string{"{1000003, 1000033}"},
        "factorint splits a nontrivial uint64 semiprime");
    tests.expectEqual(eval(session, "totient[1]"), std::string{"1"},
        "totient one is exact");
    tests.expectEqual(eval(session, "totient[9]"), std::string{"6"},
        "totient uses exact prime factorization");
    tests.expect(evalError(session, "factorint[0]").type() == error::CalcErrorType::Domain,
        "factorint rejects zero");
    tests.expect(evalError(session, "totient[0]").type() == error::CalcErrorType::Domain,
        "totient requires a positive integer");

    tests.expectEqual(eval(session, "dft[{1,2,3,4}]"),
        std::string{"{10, -2+2I, -2, -2-2I}"},
        "DFT remains exact for fourth roots of unity");
    tests.expectEqual(eval(session, "fft[{1,2,3,4}]"),
        std::string{"{10, -2+2I, -2, -2-2I}"},
        "radix-2 FFT agrees with exact DFT");
    tests.expectEqual(eval(session, "ifft[fft[{1+I,2-I,3+2I,4-3I}]]"),
        std::string{"{1+I, 2-I, 3+2I, 4-3I}"},
        "exact radix-2 FFT round trip is lossless");
    tests.expectEqual(eval(session, "fft[{1,2,3}]"), eval(session, "dft[{1,2,3}]"),
        "non-power-of-two FFT falls back to exact DFT");
    tests.expectEqual(eval(session, "convolve[{1,2},{3,4}]"), std::string{"{3, 10, 8}"},
        "convolution is exact");
    tests.expectEqual(eval(session, "convolve[{a,b},{c,d}]"),
        std::string{"{a c, a d+b c, b d}"},
        "convolution supports symbolic coefficients");
    tests.expectEqual(eval(session, "N[dft[{1,2,3}],12]"),
        std::string{"{6, -1.5+0.866025403784I, -1.5-0.866025403784I}"},
        "N recursively approximates exact transform arrays");
    tests.expectEqual(eval(session, "N[fft[{1,2,3}],12]"),
        std::string{"{6, -1.5+0.866025403784I, -1.5-0.866025403784I}"},
        "N pushes requested precision into the FFT backend before exact expansion");
    tests.expectEqual(eval(session, "fft[N[{1,2,3,4},12]]"),
        std::string{"{10, -2+2I, -2, -2-2I}"},
        "FFT dispatches approximate inputs directly to the certified backend");

    // Fourier位相は数学上rad固定であり、ユーザーの三角函数既定単位から独立する。
    session.setDefaultAngleUnit(mathematics::AngleUnit::Degree);
    tests.expectEqual(eval(session, "fft[{1,2,3,4}]"),
        std::string{"{10, -2+2I, -2, -2-2I}"},
        "FFT definition is invariant under the session angle mode");

    tests.expect(evalError(session, "dft[{{1,2},{3,4}}]").type() == error::CalcErrorType::Type,
        "DFT requires a rank-1 array");


    symbols::SymbolTable cacheSymbols;
    const auto cacheRegistry = evaluation::BuiltinRegistry::defaults(cacheSymbols);
    const auto cacheMath = mathematics::MathRegistry::defaults(cacheSymbols, cacheRegistry);
    const auto cacheAngles = mathematics::defaultAngleSemantics();
    builtins::FourierTransformCache cache;
    const auto integerExpr = [](std::int64_t value) {
        return expression::Expr{numeric::Number{numeric::BigInt{value}}};
    };
    const expression::Expr four = expression::Expr::array(
        {4}, {integerExpr(1), integerExpr(2), integerExpr(3), integerExpr(4)});
    const std::array<expression::Expr, 1> fftArguments{four};
    const expression::Expr firstTransform = builtins::evaluateFft(
        fftArguments, cacheRegistry, cacheMath, cacheAngles, cache);
    tests.expectEqual(cache.planCount(), std::size_t{1},
        "FFT cache: first radix-2 transform creates one plan");

    const std::array<expression::Expr, 1> ifftArguments{firstTransform};
    static_cast<void>(builtins::evaluateIfft(
        ifftArguments, cacheRegistry, cacheMath, cacheAngles, cache));
    tests.expectEqual(cache.planCount(), std::size_t{1},
        "FFT cache: IFFT reuses the same-size transform plan");

    const expression::Expr eight = expression::Expr::array(
        {8}, {integerExpr(1), integerExpr(2), integerExpr(3), integerExpr(4),
            integerExpr(5), integerExpr(6), integerExpr(7), integerExpr(8)});
    const std::array<expression::Expr, 1> fftEightArguments{eight};
    static_cast<void>(builtins::evaluateFft(
        fftEightArguments, cacheRegistry, cacheMath, cacheAngles, cache));
    tests.expectEqual(cache.planCount(), std::size_t{2},
        "FFT cache: a different radix-2 size creates a separate plan");
    cache.clear();
    tests.expectEqual(cache.planCount(), std::size_t{0},
        "FFT cache: plans can be cleared with evaluator lifetime semantics");
}

} // namespace mmcal::tests
