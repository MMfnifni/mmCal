// combinatorics・signalの回帰テスト
#include "combinatorics_signal_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "mathematics/angle.hpp"
#include "test_framework.hpp"

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
        std::string{"{6, -1.500000000000+0.866025403784I, -1.500000000000-0.866025403784I}"},
        "N recursively approximates exact transform arrays");

    // Fourier位相は数学上rad固定であり、ユーザーの三角函数既定単位から独立する。
    session.setDefaultAngleUnit(mathematics::AngleUnit::Degree);
    tests.expectEqual(eval(session, "fft[{1,2,3,4}]"),
        std::string{"{10, -2+2I, -2, -2-2I}"},
        "FFT definition is invariant under the session angle mode");

    tests.expect(evalError(session, "dft[{{1,2},{3,4}}]").type() == error::CalcErrorType::Type,
        "DFT requires a rank-1 array");
}

} // namespace mmcal::tests
