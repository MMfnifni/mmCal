// random・functionの回帰テスト
#include "random_function_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "numeric/big_int.hpp"
#include "numeric/real_number.hpp"
#include "test_framework.hpp"

#include <array>
#include <stdexcept>
#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

[[nodiscard]] std::string eval(kernel::KernelSession& session, std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

[[nodiscard]] error::CalcError evalError(kernel::KernelSession& session, std::string_view source) {
    try { static_cast<void>(session.evaluate(source)); }
    catch (const error::CalcError& exception) { return exception; }
    throw std::logic_error("Expected CalcError");
}

} // namespace

void runRandomFunctionTests(TestRunner& tests) {
    kernel::KernelSession first;
    kernel::KernelSession second;

    tests.expectEqual(eval(first, "randSeed[42]"), std::string{"42"},
        "randSeed returns the explicit seed");
    tests.expectEqual(eval(second, "randSeed[42]"), std::string{"42"},
        "equal seeds initialize independent sessions identically");

    for (int i = 0; i < 6; ++i)
        tests.expectEqual(eval(first, "rand[]"), eval(second, "rand[]"),
            "seeded rand sequence is reproducible across sessions");

    static_cast<void>(eval(first, "randSeed[123456789]"));
    const std::string firstSample = eval(first, "rand[]");
    static_cast<void>(eval(first, "rand[]"));
    static_cast<void>(eval(first, "randSeed[123456789]"));
    tests.expectEqual(eval(first, "rand[]"), firstSample,
        "reseeding rewinds the deterministic stream");

    static_cast<void>(eval(first, "randSeed[7]"));
    for (int i = 0; i < 32; ++i) {
        const auto value = first.evaluate("rand[]");
        tests.expect(value.isNumber() && value.asNumber().isReal()
            && value.asNumber().asReal() >= numeric::RealNumber{numeric::BigInt{}}
            && value.asNumber().asReal() < numeric::RealNumber{numeric::BigInt{1}},
            "rand returns an exact value in [0,1)");
    }

    static_cast<void>(eval(first, "randSeed[7]"));
    for (int i = 0; i < 32; ++i) {
        const auto value = first.evaluate("rand[-3,5]");
        tests.expect(value.isNumber() && value.asNumber().isReal()
            && value.asNumber().asReal() >= numeric::RealNumber{numeric::BigInt{-3}}
            && value.asNumber().asReal() < numeric::RealNumber{numeric::BigInt{5}},
            "bounded rand stays inside [lo,hi)");
    }
    tests.expectEqual(eval(first, "rand[0]"), std::string{"0"},
        "rand zero upper bound is exactly zero");
    tests.expect(evalError(first, "rand[-1]").type() == error::CalcErrorType::Domain,
        "one-argument rand rejects a negative upper bound");
    tests.expect(evalError(first, "rand[2,1]").type() == error::CalcErrorType::Domain,
        "two-argument rand rejects reversed bounds");

    static_cast<void>(eval(first, "randSeed[91]"));
    for (int i = 0; i < 64; ++i) {
        const auto value = first.evaluate("randint[-2,2]");
        tests.expect(value.isNumber() && value.asNumber().isReal()
            && value.asNumber().asReal().isInteger()
            && value.asNumber().asReal() >= numeric::RealNumber{numeric::BigInt{-2}}
            && value.asNumber().asReal() <= numeric::RealNumber{numeric::BigInt{2}},
            "randint returns an inclusive exact integer sample");
    }
    tests.expectEqual(eval(first, "randint[0]"), std::string{"0"},
        "randint zero bound is exactly zero");
    tests.expect(evalError(first, "randint[3,2]").type() == error::CalcErrorType::Domain,
        "randint rejects reversed bounds");

    static_cast<void>(eval(first, "randSeed[20260809]"));
    const numeric::BigInt hugeLower = numeric::BigInt::parse("1000000000000000000000000000000");
    const numeric::BigInt hugeUpper = numeric::BigInt::parse("1000000000000000000000000000010");
    for (int i = 0; i < 16; ++i) {
        const auto value = first.evaluate(
            "randint[1000000000000000000000000000000,1000000000000000000000000000010]");
        tests.expect(value.isNumber() && value.asNumber().isReal()
            && value.asNumber().asReal().isInteger()
            && value.asNumber().asReal().asInteger() >= hugeLower
            && value.asNumber().asReal().asInteger() <= hugeUpper,
            "randint supports inclusive ranges wider than uint64");
    }

    static_cast<void>(eval(first, "randSeed[314159]"));
    const std::string selectedA = eval(first, "choice[{2,3,5,7,11}] ");
    static_cast<void>(eval(first, "randSeed[314159]"));
    tests.expectEqual(eval(first, "choice[2,3,5,7,11]"), selectedA,
        "array and variadic choice share the same selection rule");
    tests.expect(evalError(first, "choice[{}]").type() == error::CalcErrorType::Domain,
        "choice rejects an empty array");

    static_cast<void>(eval(first, "randSeed[271828]"));
    const std::string normalA = eval(first, "randn[]");
    static_cast<void>(eval(first, "randSeed[271828]"));
    tests.expectEqual(eval(first, "randn[]"), normalA,
        "randn is reproducible as an exact Box-Muller expression");
    static_cast<void>(eval(first, "randSeed[271828]"));
    const std::string approximateNormal = eval(first, "N[randn[],12]");
    tests.expect(!approximateNormal.empty() && approximateNormal.find('[') == std::string::npos,
        "randn can be certified to a numeric decimal with N");
    tests.expectEqual(eval(first, "randn[5,0]"), std::string{"5"},
        "zero normal standard deviation returns the mean exactly");
    tests.expect(evalError(first, "randn[0,-1]").type() == error::CalcErrorType::Domain,
        "randn rejects a negative sigma");

    kernel::KernelSession entropySession;
    const auto entropySeed = entropySession.evaluate("randSeed[]");
    tests.expect(entropySeed.isNumber() && entropySeed.asNumber().isReal()
        && entropySeed.asNumber().asReal().isInteger(),
        "randSeed without arguments returns a replayable integer seed");
}

} // namespace mmcal::tests
