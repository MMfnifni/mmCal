// 記号積分integrateの回帰テスト
#include "integration_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "test_framework.hpp"

#include <stdexcept>
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

void runIntegrationTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(eval(session, "integrate[5,x]"), std::string{"5x"},
        "integrates constants without an integration constant");
    tests.expectEqual(eval(session, "integrate[x,x]"), std::string{"x^2/2"},
        "integrates x");
    tests.expectEqual(eval(session, "integrate[x^2+2*x+1,x]"),
        std::string{"x+x^2+x^3/3"},
        "integrates exact polynomials termwise");
    tests.expectEqual(eval(session, "integrate[1/x,x]"), std::string{"log[x]"},
        "integrates reciprocal by principal Log");
    tests.expectEqual(eval(session, "integrate[exp[x],x]"), std::string{"exp[x]"},
        "integrates exp");
    tests.expectEqual(eval(session, "integrate[sin[x],x]"), std::string{"-cos[x]"},
        "integrates sine in the default radian semantics");
    tests.expectEqual(eval(session, "integrate[cos[x],x]"), std::string{"sin[x]"},
        "integrates cosine");
    tests.expectEqual(eval(session, "integrate[2*x*cos[x^2],x]"), std::string{"sin[x^2]"},
        "reverse chain rule is discovered and verified through D");
    tests.expectEqual(eval(session, "integrate[2*x/(x^2+1),x]"), std::string{"log[1+x^2]"},
        "logarithmic substitution is discovered and verified");
    tests.expectEqual(eval(session, "integrate[1/(1+x^2),x]"), std::string{"atan[x]"},
        "inverse tangent derivative knowledge is reusable by integration");
    tests.expectEqual(eval(session, "integrate[log[x],x]"), std::string{"x log[x]-x"},
        "integrates Log by the standard primitive");
    tests.expectEqual(eval(session, "integrate[x*exp[x],x]"), std::string{"x exp[x]-exp[x]"},
        "integration by parts handles polynomial times exp");
    tests.expectEqual(eval(session, "integrate[erf[x],x]"),
        std::string{"x erf[x]+exp[-x^2]/sqrt[Pi]"},
        "integrates erf using shared special-function derivative knowledge");
    tests.expectEqual(eval(session, "integrate[sqrt[x],x]"),
        std::string{"2/3x sqrt[x]"},
        "integrates principal square root when D verifies the result");

    tests.expectEqual(eval(session, "integrate[(2*x+3)^5,x]"),
        std::string{"(3+2x)^6/12"},
        "affine power rule keeps a compact exact primitive");
    tests.expectEqual(eval(session, "integrate[1/(2*x+3),x]"),
        std::string{"log[3+2x]/2"},
        "affine reciprocal integrates to principal Log");
    tests.expectEqual(eval(session, "integrate[exp[2*x+1],x]"),
        std::string{"exp[1+2x]/2"},
        "affine exponential chain rule");
    tests.expectEqual(eval(session, "integrate[tan[x],x]"),
        std::string{"-log[cos[x]]"},
        "tangent primitive uses principal Log");
    tests.expectEqual(eval(session, "integrate[sech[x],x]"),
        std::string{"atan[sinh[x]]"},
        "hyperbolic secant primitive reuses inverse-trig knowledge");
    tests.expectEqual(eval(session, "integrate[1/(x^2-1),x]"),
        std::string{"-atanh[x]"},
        "polynomial proportionality proves a reverse atanh chain");
    tests.expectEqual(eval(session, "integrate[x*log[x],x]"),
        std::string{"x^2/2log[x]-x^2/4"},
        "monomial times Log uses the exact integration-by-parts formula");
    tests.expectEqual(eval(session, "integrate[acosh[x],x]"),
        std::string{"x acosh[x]-sqrt[1+x]sqrt[x-1]"},
        "principal acosh primitive uses the same branch convention as D");
    tests.expectEqual(eval(session, "integrate[sin[x Deg],x]"),
        std::string{"-(180/Pi cos[x Deg])"},
        "explicit Degree suffix is integrated with the correct radian scale");

    tests.expectEqual(eval(session, "integrate[sin[x]^2,x]"),
        std::string{"(x-sin[2x]/2)/2"},
        "trigonometric square uses an exact double-angle identity");
    tests.expectEqual(eval(session, "integrate[tan[x]^2,x]"),
        std::string{"tan[x]-x"},
        "tangent square reduces through secant squared");
    tests.expectEqual(eval(session, "integrate[sin[x]*cos[x],x]"),
        std::string{"sin[x]^2/2"},
        "reverse-chain candidate set recognizes trig products");
    tests.expectEqual(eval(session, "integrate[2*x/sqrt[1+x^2],x]"),
        std::string{"2sqrt[1+x^2]"},
        "reverse-chain candidate set recognizes square-root derivatives");
    tests.expectEqual(eval(session, "integrate[1/(x^2+4),x]"),
        std::string{"atan[x/2]/2"},
        "exact quadratic reciprocal completes the square to atan");
    tests.expectEqual(eval(session, "integrate[1/(4-x^2),x]"),
        std::string{"atanh[x/2]/2"},
        "exact quadratic reciprocal completes the square to atanh");
    tests.expectEqual(eval(session, "integrate[1/(x+1)^2,x]"),
        std::string{"-1/(1+x)"},
        "repeated linear quadratic uses the rational primitive");
    tests.expectEqual(eval(session, "integrate[1/sqrt[4-x^2],x]"),
        std::string{"asin[x/2]"},
        "negative-leading quadratic root uses asin when the positive scale is provable");
    tests.expectEqual(eval(session, "integrate[1/sqrt[x^2+4],x]"),
        std::string{"asinh[x/2]"},
        "positive-definite quadratic root uses asinh");
    tests.expectEqual(eval(session, "integrate[log[2,x],x]"),
        std::string{"(x log[x]-x)/log[2]"},
        "arbitrary-base logarithm reuses the natural-log primitive");

    tests.expectEqual(eval(session, "integrate[(x+1)/(x+2),x]"),
        std::string{"x-log[2+x]"},
        "rational function with a linear denominator uses exact polynomial division");
    tests.expectEqual(eval(session, "integrate[(x+1)/(x^2+4),x]"),
        std::string{"log[4+x^2]/2+atan[x/2]/2"},
        "rational function over a quadratic splits into log and reciprocal parts");

    kernel::KernelSession degreeIntegrals;
    degreeIntegrals.setDefaultAngleUnit(mathematics::AngleUnit::Degree);
    tests.expectEqual(eval(degreeIntegrals, "integrate[1/(1+x^2),x]"),
        std::string{"Pi/180atan[x]"},
        "inverse-trig primitives respect a Degree session without changing the integrand");

    tests.expectEqual(eval(session, "integrate[x^2,{x,0,1}]"), std::string{"1/3"},
        "exact definite polynomial integral");
    tests.expectEqual(eval(session, "integrate[sin[x],{x,0,Pi}]"), std::string{"2"},
        "exact definite trigonometric integral");
    tests.expectEqual(eval(session, "integrate[1/x,{x,1,2}]"), std::string{"log[2]"},
        "definite integral with a certified nonzero denominator");

    tests.expectEqual(eval(session, "integrate[log[x],{x,1,Pi}]"),
        std::string{"1+Pi log[Pi]-Pi"},
        "certified transcendental bounds permit safe exact definite evaluation");

    const std::string pole = eval(session, "integrate[tan[x],{x,0,2}]");
    tests.expect(pole.find("integrate[") == 0,
        "definite integration refuses an interval crossing a trig pole");
    tests.expect(findDiagnostic(session, "integrate::unevaluated") != nullptr,
        "unsafe trig definite integral emits a warning");

    const std::string singular = eval(session, "integrate[1/x,{x,-1,1}]");
    tests.expect(singular.find("integrate[") == 0,
        "integral crossing a pole stays unevaluated");
    tests.expect(findDiagnostic(session, "integrate::unevaluated") != nullptr,
        "unresolved definite integral emits a warning");

    tests.expectEqual(eval(session, "integrate[1/sqrt[x^2-1],x]"),
        std::string{"log[sqrt[x^2-1]+x]"},
        "branch-sensitive quadratic root uses a verified local primitive without a global sqrt identity");

    const std::string unsupported = eval(session, "integrate[gamma[x],x]");
    tests.expect(unsupported.find("integrate[") == 0,
        "unsupported special-function antiderivative stays symbolic");
    tests.expect(findDiagnostic(session, "integrate::unevaluated") != nullptr,
        "unsupported indefinite integral emits a warning");

    tests.expectEqual(eval(session, "D[integrate[f[x],x],x]"), std::string{"f[x]"},
        "D applies the fundamental theorem to an unresolved indefinite integral");
    tests.expectEqual(eval(session, "D[integrate[t^2,{t,0,x}],x]"), std::string{"x^2"},
        "D applies the Leibniz rule to a variable upper bound");

    kernel::KernelSession boundVariable;
    static_cast<void>(eval(boundVariable, "x:=5"));
    tests.expectEqual(eval(boundVariable, "integrate[x,x]"), std::string{"x^2/2"},
        "integration variable remains held even when a global value exists");
}

} // namespace mmcal::tests
