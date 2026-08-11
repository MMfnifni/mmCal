// 指数・対数・特殊値のexact計算の回帰テスト
#include "exact_transcendental_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

[[nodiscard]] std::string evaluateAndFormat(
    kernel::KernelSession& session,
    std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

[[nodiscard]] error::CalcError evaluateError(
    kernel::KernelSession& session,
    std::string_view source) {
    try {
        static_cast<void>(session.evaluate(source));
    }
    catch (const error::CalcError& exception) {
        return exception;
    }
    throw std::logic_error("Expected CalcError was not thrown");
}

} // namespace

void runExactTranscendentalTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(evaluateAndFormat(session, "arg[1]"), std::string{"0 Rad"},
        "Exact transcendental: principal Arg of positive real is zero radians");
    tests.expectEqual(evaluateAndFormat(session, "arg[-1]"), std::string{"Pi Rad"},
        "Exact transcendental: principal Arg uses +Pi radians on negative real axis");
    tests.expectEqual(evaluateAndFormat(session, "arg[I]"), std::string{"Pi/2 Rad"},
        "Exact transcendental: Arg[I] is Pi/2 radians");
    tests.expectEqual(evaluateAndFormat(session, "arg[-I]"), std::string{"-Pi/2 Rad"},
        "Exact transcendental: Arg[-I] is -Pi/2 radians");
    tests.expectEqual(evaluateAndFormat(session, "arg[1 + I]"), std::string{"Pi/4 Rad"},
        "Exact transcendental: diagonal complex Arg is exact radians");
    tests.expectEqual(evaluateAndFormat(session, "sin[arg[I]]"), std::string{"1"},
        "Exact transcendental: Arg carries explicit radians into explicit-radian trig");

    const error::CalcError argZero = evaluateError(session, "arg[0]");
    tests.expect(argZero.type() == error::CalcErrorType::Domain,
        "Exact transcendental: Arg[0] is undefined");

    tests.expectEqual(evaluateAndFormat(session, "log[1]"), std::string{"0"},
        "Exact transcendental: Log[1] is zero");
    tests.expectEqual(evaluateAndFormat(session, "log[E]"), std::string{"1"},
        "Exact transcendental: E is tied to natural logarithm exactly");
    tests.expectEqual(evaluateAndFormat(session, "log[-1]"), std::string{"I Pi"},
        "Exact transcendental: principal Log[-1] uses +I Pi");
    tests.expectEqual(evaluateAndFormat(session, "log[I]"), std::string{"I Pi/2"},
        "Exact transcendental: principal Log[I] is I Pi/2");
    tests.expectEqual(evaluateAndFormat(session, "log[-E]"), std::string{"1+I Pi"},
        "Exact transcendental: principal Log[-E] keeps exact real and branch parts");

    const error::CalcError logZero = evaluateError(session, "log[0]");
    tests.expect(logZero.type() == error::CalcErrorType::Domain,
        "Exact transcendental: Log[0] is a domain error, not Infinity");

    tests.expectEqual(evaluateAndFormat(session, "log[10,1000]"), std::string{"3"},
        "Exact transcendental: arbitrary-base log recognizes an integer power");
    tests.expectEqual(evaluateAndFormat(session, "log[2,1/8]"), std::string{"-3"},
        "Exact transcendental: arbitrary-base log recognizes a negative integer exponent");
    tests.expectEqual(evaluateAndFormat(session, "log[4,2]"), std::string{"1/2"},
        "Exact transcendental: arbitrary-base log recognizes a reciprocal integer exponent");
    tests.expectEqual(evaluateAndFormat(session, "log[1/4,2]"), std::string{"-1/2"},
        "Exact transcendental: reciprocal bases preserve the exponent sign");
    tests.expectEqual(evaluateAndFormat(session, "log[E,E]"), std::string{"1"},
        "Exact transcendental: base E agrees with natural logarithm");

    const error::CalcError logBaseZero = evaluateError(session, "log[0,10]");
    tests.expect(logBaseZero.type() == error::CalcErrorType::Domain,
        "Exact transcendental: logarithm base zero is rejected");
    const error::CalcError logBaseOne = evaluateError(session, "log[1,10]");
    tests.expect(logBaseOne.type() == error::CalcErrorType::Domain,
        "Exact transcendental: logarithm base one is rejected");
    const error::CalcError logValueZero = evaluateError(session, "log[10,0]");
    tests.expect(logValueZero.type() == error::CalcErrorType::Domain,
        "Exact transcendental: arbitrary-base logarithm rejects value zero");

    tests.expectEqual(evaluateAndFormat(session, "exp[0]"), std::string{"1"},
        "Exact transcendental: Exp[0] is one");
    tests.expectEqual(evaluateAndFormat(session, "exp[1]"), std::string{"E"},
        "Exact transcendental: E is Exp[1]");
    tests.expectEqual(evaluateAndFormat(session, "exp[I Pi]"), std::string{"-1"},
        "Exact transcendental: Euler identity is exact");
    tests.expectEqual(evaluateAndFormat(session, "exp[I Pi / 2]"), std::string{"I"},
        "Exact transcendental: quarter-turn exponential is exact");

    tests.expectEqual(evaluateAndFormat(session, "Infinity"), std::string{"Infinity"},
        "Exact transcendental: Infinity is not predefined; temporary free-symbol mode preserves the name");
}

} // namespace mmcal::tests
