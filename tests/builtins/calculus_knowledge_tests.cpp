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
    tests.expectEqual(eval(session, "integrate[exp[-x^2],{x,0,Infinity}]"),
        std::string{"sqrt[Pi]/2"},
        "Gaussian semi-infinite integral closes through the canonical erf endpoint limit");
    tests.expectEqual(eval(session, "integrate[exp[-x^2],{x,-Infinity,Infinity}]"),
        std::string{"sqrt[Pi]"},
        "two-sided Gaussian integral closes exactly through erf infinity limits");

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

    tests.expectEqual(eval(session, "root[{-4,0,2},2]"),
        std::string{"root[{-2, 0, 1}, 2]"},
        "root canonicalizes its defining polynomial to monic form");
    tests.expectEqual(eval(session, "root[{4,0,-4,0,1},2]"),
        std::string{"root[{-2, 0, 1}, 2]"},
        "root removes repeated polynomial factors because indices count distinct real roots");
    tests.expectEqual(eval(session, "N[root[{-2,0,1},2],30]"),
        std::string{"1.41421356237309504880168872421"},
        "root has certified arbitrary-precision numerical refinement");
    tests.expectEqual(eval(session, "N[root[{0,-2,0,1},2],30]"),
        std::string{"0"},
        "root refinement recognizes an exact zero root without relative-precision stalling");
    tests.expectEqual(eval(session, "solve[x^5-x+1==0,x,Real]"),
        std::string{"{x==root[{1, -1, 0, 0, 0, 1}, 1]}"},
        "Real Solve falls back to an exact Root representation for unresolved rational polynomials");
    tests.expectEqual(eval(session, "solve[(x^2-2)^2==0,x,Real]"),
        std::string{"{x==root[{-2, 0, 1}, 1], x==root[{-2, 0, 1}, 2]}"},
        "algebraic Root fallback preserves distinct repeated-polynomial roots canonically");
    tests.expectEqual(eval(session, "solve[x^4+1==0,x,Real]"),
        std::string{"{}"},
        "Real Root isolation proves that a rational polynomial has no real roots");
    const std::string complexAlgebraicUnresolved = eval(session, "solve[x^5-x+1==0,x]");
    tests.expect(complexAlgebraicUnresolved == "UnresolvedSolutionSet[x]"
            && findDiagnostic(session, "solve::unresolved") != nullptr,
        "general complex algebraic Root isolation remains intentionally unresolved");
    const std::string rootExplanation = eval(session, "explain[root[{-2,0,1},2]]");
    tests.expect(rootExplanation.find("{\"Kind\", \"AlgebraicNumber\"}") != std::string::npos
            && rootExplanation.find("{\"PolynomialDegree\", 2}") != std::string::npos
            && rootExplanation.find("{\"RootIndex\", 2}") != std::string::npos,
        "explain exposes Root as an exact real algebraic number without evaluating it numerically");

    tests.expectEqual(eval(session, "solve[sin[x]==0,x,Real]"),
        std::string{"{x==Pi k where k in Integer}"},
        "Solve represents real sine zeros as integer-parameter solution families");
    tests.expectEqual(eval(session, "solve[cos[x]==0,x,Real]"),
        std::string{"{x==Pi/2+Pi k where k in Integer}"},
        "Solve represents real cosine zeros as periodic families");
    tests.expectEqual(eval(session, "solve[tan[x]==1,x,Real]"),
        std::string{"{x==Pi/4+Pi k where k in Integer}"},
        "Solve uses the half-turn period of tangent");
    tests.expectEqual(eval(session, "solve[sin[2x+1]==0,x,Real]"),
        std::string{"{x==-(1-Pi k)/2 where k in Integer}"},
        "Solve propagates periodic targets through an exact affine argument");
    tests.expectEqual(eval(session, "solve[sin[x]==1,x,Real]"),
        std::string{"{x==Pi/2+2Pi k where k in Integer}"},
        "Sine endpoint targets avoid duplicate periodic branches");
    tests.expectEqual(eval(session, "solve[cos[x]==-1,x,Real]"),
        std::string{"{x==Pi+2Pi k where k in Integer}"},
        "Cosine endpoint targets avoid duplicate periodic branches");
    tests.expectEqual(eval(session, "solve[sin[x]==2,x,Real]"), std::string{"{}"},
        "periodic Solve rejects exact targets outside the real range");
    tests.expectEqual(eval(session, "solve[sin[x^2]==0,x,Real]"),
        std::string{"UnresolvedSolutionSet[x]"},
        "first periodic Solve implementation remains conservative for nonlinear arguments");
    tests.expectEqual(eval(session, "solve[sin[x]==k,x,Real]").find("where k1 in Integer") != std::string::npos, true,
        "periodic Solve chooses a fresh formal parameter when k already occurs in the relation");

    kernel::KernelSession degreePeriodic;
    static_cast<void>(eval(degreePeriodic, "angleMode[Deg]"));
    tests.expectEqual(eval(degreePeriodic, "solve[sin[x]==0,x,Real]"),
        std::string{"{x==180k where k in Integer}"},
        "periodic Solve respects Degree session semantics");
    tests.expectEqual(eval(degreePeriodic, "solve[tan[x]==1,x,Real]"),
        std::string{"{x==45+180k where k in Integer}"},
        "tangent periodic families respect Degree semantics");
}

} // namespace mmcal::tests
