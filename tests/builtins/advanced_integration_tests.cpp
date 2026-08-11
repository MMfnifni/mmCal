// advanced・integrationの回帰テスト
#include "advanced_integration_tests.hpp"

#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "test_framework.hpp"

#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

[[nodiscard]] std::string eval(kernel::KernelSession& session, std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

[[nodiscard]] const evaluation::EvaluationDiagnostic* findDiagnostic(
    const kernel::KernelSession& session,
    std::string_view code) {
    for (const auto& diagnostic : session.diagnostics())
        if (diagnostic.code == code)
            return &diagnostic;
    return nullptr;
}

} // namespace

void runAdvancedIntegrationTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(eval(session, "integrate[x^2+sin[x],x]"),
        std::string{"x^3/3-cos[x]"},
        "integration is linear across polynomial and trigonometric terms");
    tests.expectEqual(eval(session, "D[sin[x],x]"), std::string{"cos[x]"},
        "canonical trigonometric names share derivative knowledge");

    tests.expectEqual(eval(session, "integrate[E^x*cos[x],x]"),
        std::string{"(cos[x]+sin[x])/2*exp[x]"},
        "E^x canonicalizes to exp and exponential-trigonometric products integrate exactly");
    tests.expectEqual(eval(session, "D[E^x,x]"), std::string{"exp[x]"},
        "D reuses the E^x to exp canonical identity");
    tests.expectEqual(eval(session, "(cos[x]+sin[x])/2exp[x]"),
        std::string{"(cos[x]+sin[x])/2*exp[x]"},
        "formatted exponential products can be parsed back through implicit multiplication");
    tests.expectEqual(eval(session,
        "integrate[E^x*cos[x],x]==(cos[x]+sin[x])/2exp[x]"),
        std::string{"True"},
        "the exact exponential-trigonometric antiderivative compares equal after round-trip parsing");


    tests.expectEqual(eval(session, "D[integrate[1/(x^3+1),x],x]"),
        std::string{"1/(1+x^3)"},
        "exact partial fractions differentiate back to the cubic rational integrand");
    tests.expectEqual(eval(session, "integrate[1/(x^3+1),{x,0,1}]"),
        std::string{"log[2]/3+Pi/6*2sqrt[3]/3"},
        "cubic rational definite integral reduces to exact log and Pi terms");
    tests.expectEqual(eval(session, "integrate[1/((x+1)^2*(x-1)),x]"),
        std::string{"log[-1+x]/4-log[1+x]/4+(1+x)^(-1)/2"},
        "partial fractions support repeated linear factors");

    tests.expectEqual(eval(session, "integrate[1/sqrt[x^2-1],x]"),
        std::string{"log[sqrt[x^2-1]+x]"},
        "local antiderivative is accepted without installing a global sqrt factorization identity");
    tests.expectEqual(eval(session, "D[integrate[1/sqrt[x^2-1],x],x]"),
        std::string{"1/sqrt[x^2-1]"},
        "D applies the fundamental theorem to the local quadratic-root primitive request");

    const std::string nested = eval(session, "integrate[sqrt[x+sqrt[x]],x]");
    tests.expect(nested.find("integrate[") == std::string::npos
            && nested.find("log[") != std::string::npos,
        "sqrt substitution handles nested quadratic radicals symbolically");
    tests.expect(findDiagnostic(session, "integrate::unevaluated") == nullptr,
        "supported nested radical does not emit an unevaluated warning");
    tests.expect(nested.find("+-") == std::string::npos
            && nested.find(" + ") == std::string::npos
            && nested.find(" - ") == std::string::npos
            && nested.find(" / ") == std::string::npos
            && nested.find(" ^ ") == std::string::npos,
        "formatter removes redundant operator spacing and plus-negative forms");
    tests.expectEqual(eval(session, nested), nested,
        "compact nested-radical formatting round-trips through the parser");

    const std::string partial = eval(session, "integrate[x^2+gamma[x],x]");
    tests.expect(partial.find("x^3/3") != std::string::npos
            && partial.find("integrate[gamma[x], x]") != std::string::npos,
        "linear integration preserves solved terms when one term remains unresolved");
    tests.expect(findDiagnostic(session, "integrate::unevaluated") != nullptr,
        "partially evaluated integral still emits an explicit warning");

    tests.expectEqual(eval(session, "D[sin[x],{x,4}]"), std::string{"sin[x]"},
        "D supports exact higher derivative order specifications");
    tests.expectEqual(eval(session, "D[x^2*y^3,x,y]"), std::string{"6x y^2"},
        "D supports sequential mixed derivative specifications");
    tests.expectEqual(eval(session, "D[x^5,{x,0}]"), std::string{"x^5"},
        "zeroth derivative order returns the expression unchanged");
}

} // namespace mmcal::tests
