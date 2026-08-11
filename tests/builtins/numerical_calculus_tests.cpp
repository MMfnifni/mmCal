// 数値微分・保証付き数値積分の回帰テスト
#include "numerical_calculus_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
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

void runNumericalCalculusTests(TestRunner& tests) {
    kernel::KernelSession session;
    tests.expectEqual(eval(session, "diff[x^2,x,3]"), std::string{"6.0000000000000000"},
        "Numerical calculus: diff reuses symbolic D and certified point evaluation");
    tests.expectEqual(eval(session, "diff[sin[x],x,0,20]"),
        std::string{"1.00000000000000000000"},
        "Numerical calculus: diff respects default radian semantics");
    tests.expectEqual(eval(session, "diff[asin[x],x,0,20]"),
        std::string{"1.00000000000000000000"},
        "Numerical calculus: inverse trig derivative is radian-scaled by default");
    tests.expectEqual(eval(session, "nintegrate[x^2,{x,0,1},12]"),
        std::string{"0.333333333333"},
        "Numerical calculus: polynomial integral is certified");
    tests.expectEqual(eval(session, "nintegrate[sin[x],{x,0,Pi},12]"),
        std::string{"2.000000000000"},
        "Numerical calculus: integration uses radian trigonometry by default");
    tests.expectEqual(eval(session, "nintegrate[cos[x],{x,0,Pi/2},12]"),
        std::string{"1.000000000000"},
        "Numerical calculus: cosine radian integral is certified");
    tests.expectEqual(eval(session, "nintegrate[sin[x Deg],{x,0,180},12]"),
        std::string{"114.591559026165"},
        "Numerical calculus: a unit suffix can be applied to the iterator expression");
    tests.expectEqual(eval(session, "nintegrate[exp[x],{x,0,1},12]"),
        std::string{"1.718281828459"},
        "Numerical calculus: transcendental integral is certified");
    tests.expectEqual(eval(session, "nintegrate[x^2,{x,1,0},12]"),
        std::string{"-0.333333333333"},
        "Numerical calculus: reversed limits negate the result");
    tests.expectEqual(eval(session, "nintegrate[x^2,{x,1,1}]"), std::string{"0"},
        "Numerical calculus: zero-width interval is exact zero");
    const auto singular = evalError(session, "nintegrate[1/x,{x,-1,1},8]");
    tests.expect(singular.type() == error::CalcErrorType::Domain,
        "Numerical calculus: possible singularity is rejected rather than sampled through");
    tests.expectEqual(eval(session, "nintegrate[x,{x,0,Pi},12]"),
        std::string{"4.934802200545"},
        "Numerical calculus: certified transcendental bounds are normalized to a rational parameter interval");
    tests.expectEqual(eval(session, "nintegrate[2*x,{x,0,sqrt[2]},12]"),
        std::string{"2.000000000000"},
        "Numerical calculus: certified radical bounds are supported");

    static_cast<void>(session.evaluate("x := 100"));
    tests.expectEqual(eval(session, "nintegrate[x,{x,0,1},12]"),
        std::string{"0.500000000000"},
        "Numerical calculus: iterator variable is held even when it has a global value");

    static_cast<void>(session.evaluate("lo := 0"));
    static_cast<void>(session.evaluate("hi := 1"));
    tests.expectEqual(eval(session, "nintegrate[x,{x,lo,hi},12]"),
        std::string{"0.500000000000"},
        "Numerical calculus: iterator bounds are normally evaluated");

    tests.expectEqual(eval(session, "integrate[x^2,{x,0,1}]"),
        std::string{"1/3"},
        "Numerical calculus: integrate is now the exact symbolic definite integral");

    const auto malformedIterator = evalError(session, "nintegrate[x,{x,0},8]");
    tests.expect(malformedIterator.type() == error::CalcErrorType::Type,
        "Numerical calculus: nintegrate rejects malformed iterator specifications");
}

} // namespace mmcal::tests
