// calculus・knowledgeの回帰テスト
#include "calculus_knowledge_tests.hpp"

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

void runCalculusKnowledgeTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(eval(session, "limit[sin[x]/x,x,0]"), std::string{"1"},
        "limit uses D-backed l'Hopital for sin(x)/x");
    tests.expectEqual(eval(session, "limit[(1-cos[x])/x^2,x,0]"), std::string{"1/2"},
        "repeated l'Hopital handles second-order 0/0 limits");
    tests.expectEqual(eval(session, "limit[1/x,x,0,1]"), std::string{"Infinity"},
        "right one-sided rational pole returns positive Infinity");
    tests.expectEqual(eval(session, "limit[1/x,x,0,-1]"), std::string{"-Infinity"},
        "left one-sided rational pole returns negative Infinity");
    tests.expectEqual(eval(session, "limit[abs[x]/x,x,0,1]"), std::string{"1"},
        "right-limit direction is injected into local assumptions");
    tests.expectEqual(eval(session, "limit[abs[x]/x,x,0,-1]"), std::string{"-1"},
        "left-limit direction is injected into local assumptions");
    tests.expectEqual(eval(session, "limit[atan[x],x,Infinity]"), std::string{"Pi/2"},
        "inverse-trigonometric infinity limit respects Radian semantics");
    tests.expectEqual(eval(session, "limit[exp[-x],x,Infinity]"), std::string{"0"},
        "affine exponential asymptotics are known");

    tests.expectEqual(eval(session, "limit[1/x,x,0]"), std::string{"limit[1/x, x, 0]"},
        "non-existent two-sided pole is not collapsed to a principal value");
    tests.expect(findDiagnostic(session, "limit::unevaluated") != nullptr,
        "unresolved two-sided limit emits a warning");

    tests.expectEqual(eval(session, "integrate[exp[-x],{x,0,Infinity}]"), std::string{"1"},
        "convergent semi-infinite exponential integral uses endpoint limits");
    tests.expectEqual(eval(session, "integrate[1/x^2,{x,1,Infinity}]"), std::string{"1"},
        "convergent rational improper integral is exact");
    tests.expectEqual(eval(session, "integrate[1/(1+x^2),{x,-Infinity,Infinity}]"), std::string{"Pi"},
        "two-sided improper rational integral uses atan infinity limits");
    tests.expectEqual(eval(session, "integrate[1/sqrt[x],{x,0,1}]"), std::string{"2"},
        "integrable endpoint square-root singularity is accepted");
    tests.expectEqual(eval(session, "integrate[log[x],{x,0,1}]"), std::string{"-1"},
        "logarithmic endpoint singularity is handled by a one-sided limit");

    const std::string internalPole = eval(session, "integrate[1/(x-2),{x,1,Infinity}]");
    tests.expect(internalPole.find("integrate[") == 0,
        "improper integration refuses an internal rational pole");
    tests.expect(findDiagnostic(session, "integrate::conditionsRequired") != nullptr,
        "unsafe improper integral emits an explicit warning");

    tests.expectEqual(eval(session, "integrate[abs[x],x,x>=0]"), std::string{"x^2/2"},
        "integrate accepts assumptions and reuses assumption-aware simplification");
    tests.expectEqual(eval(session, "integrate[abs[x],x,x<=0]"), std::string{"-x^2/2"},
        "integrate assumptions select the opposite absolute-value branch safely");

    tests.expectEqual(eval(session, "solve[exp[x]==2,x,Real]"), std::string{"{x==log[2]}"},
        "Solve consumes global-real inverse knowledge for exp/log");
    tests.expectEqual(eval(session, "solve[log[x]==2,x,Real]"), std::string{"{x==exp[2]}"},
        "Solve reverses real Log through Exp without losing its natural domain");
    tests.expectEqual(eval(session, "solve[sinh[3x]==2,x,Real]"), std::string{"{x==asinh[2]/3}"},
        "Solve combines inverse-function knowledge with exact polynomial solving");
    tests.expectEqual(eval(session, "solve[tanh[x]==1/2,x,Real]"), std::string{"{x==atanh[1/2]}"},
        "Tanh inversion checks its open (-1,1) real range");
    tests.expectEqual(eval(session, "solve[tanh[x]==2,x,Real]"), std::string{"{}"},
        "Solve rejects values outside Tanh's real range");
    tests.expectEqual(eval(session, "solve[exp[x]==a,x,Real]"),
        std::string{"{x==log[a] if a in Real&&a>0}"},
        "symbolic inverse solution retains the real range condition");

    tests.expectEqual(eval(session, "solve[ellipticF[x,0]==2,x]"), std::string{"{x==2}"},
        "Solve consumes exact ellipticF degeneration before polynomial solving");
    tests.expectEqual(eval(session, "solve[ellipticE[x,0]==3,x]"), std::string{"{x==3}"},
        "Solve consumes exact ellipticE degeneration before polynomial solving");
    tests.expectEqual(eval(session, "solve[ellipticPi[0,x,0]==4,x]"), std::string{"{x==4}"},
        "Solve consumes exact ellipticPi degeneration before polynomial solving");
    tests.expectEqual(eval(session, "solve[hypergeometric2F1[0,2,3,x]==1,x]"), std::string{"All"},
        "Solve recognizes a globally constant terminating 2F1 equation");

    const std::string ellipticUnresolved = eval(session, "solve[ellipticF[x,1/3]==2,x]");
    tests.expect(ellipticUnresolved == "UnresolvedSolutionSet[x]"
            && findDiagnostic(session, "solve::unresolved") != nullptr,
        "Solve does not invent a principal inverse for a general elliptic function");
    const std::string hypergeometricUnresolved = eval(
        session, "solve[hypergeometric2F1[1,1,2,x]==2,x]");
    tests.expect(hypergeometricUnresolved == "UnresolvedSolutionSet[x]"
            && findDiagnostic(session, "solve::unresolved") != nullptr,
        "Solve does not invent a global inverse for a general Gauss hypergeometric function");

    const std::string periodic = eval(session, "solve[sin[x]==0,x,Real]");
    tests.expect(periodic == "UnresolvedSolutionSet[x]"
            && findDiagnostic(session, "solve::unresolved") != nullptr,
        "periodic equations stay unresolved until integer-parameter solution families exist");
}

} // namespace mmcal::tests
