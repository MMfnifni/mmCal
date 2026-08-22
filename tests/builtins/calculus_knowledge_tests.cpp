// calculus・knowledgeの回帰テスト
#include "calculus_knowledge_tests.hpp"

#include "error/error_message.hpp"
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
    tests.expectEqual(eval(session, "limit[sin[1/x],x,0]"), std::string{"Indeterminate"},
        "two-sided essential oscillation is reported as a proved non-existent limit");
    tests.expectEqual(eval(session, "limit[sin[1/x],x,0,1]"), std::string{"Indeterminate"},
        "right-sided essential oscillation is reported as Indeterminate");
    tests.expectEqual(eval(session, "limit[sin[1/x],x,0,-1]"), std::string{"Indeterminate"},
        "left-sided essential oscillation is reported as Indeterminate");
    tests.expectEqual(eval(session, "limit[cos[x],x,Infinity]"), std::string{"Indeterminate"},
        "periodic cosine has no single real limit at positive infinity");
    tests.expectEqual(eval(session, "limit[tan[x],x,Infinity]"), std::string{"Indeterminate"},
        "periodic tangent has no single real limit at positive infinity");
    tests.expectEqual(eval(session, "limit[x sin[1/x],x,0]"), std::string{"0"},
        "squeeze theorem closes a vanishing factor times bounded sine");
    tests.expectEqual(eval(session, "limit[Ei[x],x,0,-1]"), std::string{"-Infinity"},
        "Ei approaches negative infinity from the real left at zero");
    tests.expectEqual(eval(session, "limit[Ei[x],x,0,1]"), std::string{"-Infinity"},
        "Ei approaches negative infinity from the real right at zero");
    tests.expectEqual(eval(session, "limit[Ei[x],x,0]"), std::string{"-Infinity"},
        "Ei has the same two-sided real directed limit at zero");
    tests.expectEqual(eval(session, "limit[Ei[x],x,Infinity]"), std::string{"Infinity"},
        "Ei grows exponentially on the positive real axis");
    tests.expectEqual(eval(session, "limit[Ei[x],x,-Infinity]"), std::string{"0"},
        "Ei tends to zero along the negative real axis");
    tests.expectEqual(eval(session, "limit[Ci[x],x,0,-1]"), std::string{"-Infinity"},
        "principal Ci has negative-real directed infinity at zero from the left");
    tests.expectEqual(eval(session, "limit[Ci[x],x,0,1]"), std::string{"-Infinity"},
        "Ci approaches negative infinity at zero from the positive real side");
    tests.expectEqual(eval(session, "limit[Ci[x],x,0]"), std::string{"-Infinity"},
        "Ci has the same two-sided real directed infinity at zero");
    tests.expectEqual(eval(session, "limit[Ci[x],x,Infinity]"), std::string{"0"},
        "Ci tends to zero along the positive real axis");
    tests.expectEqual(eval(session, "limit[Ci[x],x,-Infinity]"), std::string{"I Pi"},
        "principal Ci retains its branch-cut offset I Pi at negative infinity");
    tests.expectEqual(eval(session, "limit[li[x],x,0]"), std::string{"0"},
        "principal li tends to zero at its origin branch point");
    tests.expectEqual(eval(session, "limit[li[x],x,1,-1]"), std::string{"-Infinity"},
        "li approaches negative infinity from the left at its logarithmic singularity");
    tests.expectEqual(eval(session, "limit[li[x],x,1,1]"), std::string{"-Infinity"},
        "li approaches negative infinity from the right at its logarithmic singularity");
    tests.expectEqual(eval(session, "limit[li[x],x,1]"), std::string{"-Infinity"},
        "li has the same two-sided real limit at one");
    tests.expectEqual(eval(session, "limit[li[x],{x,1,1}]"), std::string{"-Infinity"},
        "limit accepts the compact {variable, point, direction} form");
    tests.expectEqual(eval(session, "limit[li[x],x,Infinity]"), std::string{"Infinity"},
        "li grows without bound on the positive real axis");
    tests.expectEqual(eval(session, "limit[li[x],x,-Infinity]"),
        std::string{"ComplexInfinity"},
        "principal li has unbounded complex magnitude along the negative real axis");

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

    tests.expectEqual(eval(session, "solve[exp[x]==2,x,Real]"), std::string{"{x == log[2]}"},
        "Solve consumes global-real inverse knowledge for exp/log");
    tests.expectEqual(eval(session, "solve[E^x==8,x,Real]"), std::string{"{x == log[8]}"},
        "Solve-safe normalization canonicalizes E^x before transcendental classification");
    tests.expectEqual(eval(session, "solve[2^x==8,x,Real]"), std::string{"{x == 3}"},
        "Solve closes positive constant-base exponential equations exactly");
    tests.expectEqual(eval(session, "solve[2^(2x+1)==8,x,Real]"), std::string{"{x == 1}"},
        "Solve combines constant-base exponential inversion with affine polynomial solving");
    tests.expectEqual(eval(session, "solve[log[x]==2,x,Real]"), std::string{"{x == exp[2]}"},
        "Solve reverses real Log through Exp without losing its natural domain");
    tests.expectEqual(eval(session, "solve[ln[x]==2,x,Real]"), std::string{"{x == exp[2]}"},
        "Solve-safe normalization makes ln share canonical log semantics");
    tests.expectEqual(eval(session, "solve[log2[x]==3,x,Real]"), std::string{"{x == 8}"},
        "Solve inverts canonicalized base-2 logarithms exactly");
    tests.expectEqual(eval(session, "solve[log10[x]==2,x,Real]"), std::string{"{x == 100}"},
        "Solve inverts canonicalized base-10 logarithms exactly");
    tests.expectEqual(eval(session, "solve[2^x==-1,x,Real]"), std::string{"{}"},
        "positive real exponentials reject non-positive real right-hand sides");
    tests.expectEqual(eval(session, "solve[(1/2)^x==4,x,Real]"), std::string{"{x == -2}"},
        "constant-base exponential inversion also handles positive bases below one");
    tests.expectEqual(eval(session, "solve[1^x==1,x,Real]"), std::string{"All"},
        "the constant base-one exponential recognizes the universal equation");
    tests.expectEqual(eval(session, "solve[1^x==2,x,Real]"), std::string{"{}"},
        "the constant base-one exponential rejects a different constant target");
    tests.expectEqual(eval(session, "solve[sinh[3x]==2,x,Real]"), std::string{"{x == asinh[2]/3}"},
        "Solve combines inverse-function knowledge with exact polynomial solving");
    tests.expectEqual(eval(session, "solve[tanh[x]==1/2,x,Real]"), std::string{"{x == atanh[1/2]}"},
        "Tanh inversion checks its open (-1,1) real range");
    tests.expectEqual(eval(session, "solve[tanh[x]==2,x,Real]"), std::string{"{}"},
        "Solve rejects values outside Tanh's real range");
    tests.expectEqual(eval(session, "solve[exp[x]==a,x,Real]"),
        std::string{"{x == log[a] if a in Real && a > 0}"},
        "symbolic inverse solution retains the real range condition");

    tests.expectEqual(eval(session, "solve[x^2==4,Real]"), std::string{"{x == 2, x == -2}"},
        "solve[equation,Real] infers a unique unknown without treating Real as the variable");
    tests.expectEqual(eval(session, "solve[1.1^x==0,Real]"), std::string{"{}"},
        "positive real exponential is certified nonzero in the Real-domain solve shorthand");
    tests.expectEqual(eval(session, "solve[1.1^x==0,x,Real]"), std::string{"{}"},
        "positive real exponential is certified nonzero with an explicit solve variable");
    tests.expectEqual(eval(session, "solve[x==1/2,x,Integer]"), std::string{"{}"},
        "Solve uses negative domain knowledge to reject a non-integral rational candidate");
    tests.expectEqual(eval(session, "solve[x==Pi,x,Integer]"), std::string{"{}"},
        "Solve uses transcendence metadata to reject an Integer candidate");
    tests.expectEqual(eval(session, "solve[x==Phi,x,Rational]"), std::string{"{}"},
        "Solve uses irrationality metadata to reject a Rational candidate");
    tests.expectEqual(eval(session, "solve[1.1^x==x^2,x,Real]"),
        std::string{"{x == -2lambertw[log[11/10]/2]/log[11/10], x == -2lambertw[-log[11/10]/2]/log[11/10], x == -2lambertw[-1, -log[11/10]/2]/log[11/10]}"},
        "Real exponential-square equations close exactly through the two real Lambert W branches");
    tests.expectEqual(eval(session, "solve[1.1^x==x^2,Real]"),
        std::string{"{x == -2lambertw[log[11/10]/2]/log[11/10], x == -2lambertw[-log[11/10]/2]/log[11/10], x == -2lambertw[-1, -log[11/10]/2]/log[11/10]}"},
        "Real solve shorthand reuses the Lambert W exponential classifier");
    tests.expectEqual(eval(session, "solve[3^x==x^2,Real]"),
        std::string{"{x == -2lambertw[log[3]/2]/log[3]}"},
        "Lambert W classifier omits non-real W branches when the real threshold is exceeded");
    tests.expectEqual(eval(session, "solve[(1/2)^x==x^2,Real]"),
        std::string{"{x == -2lambertw[-log[1/2]/2]/log[1/2], x == -2lambertw[log[1/2]/2]/log[1/2], x == -2lambertw[-1, log[1/2]/2]/log[1/2]}"},
        "Lambert W classifier handles positive bases below one without double-negation artifacts");

    tests.expectThrows<error::CalcError>([&] { (void)session.evaluate("solve[a*x==1,Real]"); },
        "solve[equation,domain] rejects ambiguous unknown inference");
    tests.expectThrows<error::CalcError>([&] { (void)session.evaluate("solve[x==1,Pi]"); },
        "solve rejects protected constants as explicit solver variables");

    tests.expectEqual(eval(session, "lambertw[0]"), std::string{"0"},
        "Lambert W principal branch knows W(0)=0 exactly");
    tests.expectEqual(eval(session, "lambertw[E]"), std::string{"1"},
        "Lambert W principal branch knows W(E)=1 exactly");
    tests.expectEqual(eval(session, "lambertw[-1/E]"), std::string{"-1"},
        "Lambert W principal branch knows the exact real branch point");
    tests.expectEqual(eval(session, "lambertw[-1,-1/E]"), std::string{"-1"},
        "Lambert W lower branch shares the exact real branch point");
    tests.expectEqual(eval(session, "N[lambertw[1],20]"),
        std::string{"0.5671432904097838730"},
        "Lambert W principal real branch has certified numerical evaluation");
    tests.expectEqual(eval(session, "N[lambertw[-1,-1/10],20]"),
        std::string{"-3.5771520639572972184"},
        "Lambert W lower real branch has certified numerical evaluation");
    tests.expectEqual(eval(session, "N[solve[1.1^x==x^2,x,Real],20]"),
        std::string{"{x == -0.95548727594562198165, x == 1.0513800237472769374, x == 95.716830168405222740}"},
        "N approximates numerically closed SolutionSet binding values while preserving variables");
    tests.expectEqual(eval(session, "D[lambertw[x],x]"),
        std::string{"cases[lambertw[x]/(x*(1+lambertw[x])) if x != 0; 1 if x == 0]"},
        "Lambert W principal derivative preserves its removable value at zero");

    tests.expectEqual(eval(session, "solve[ellipticF[x,0]==2,x]"), std::string{"{x == 2}"},
        "Solve consumes exact ellipticF degeneration before polynomial solving");
    tests.expectEqual(eval(session, "solve[ellipticE[x,0]==3,x]"), std::string{"{x == 3}"},
        "Solve consumes exact ellipticE degeneration before polynomial solving");
    tests.expectEqual(eval(session, "solve[ellipticPi[0,x,0]==4,x]"), std::string{"{x == 4}"},
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


    tests.expectEqual(eval(session, "root[{-2,0,1},2]==sqrt[2]"), std::string{"True"},
        "Root and radical representations compare by exact algebraic identity");
    tests.expectEqual(eval(session, "root[{-2,0,1},1]==-sqrt[2]"), std::string{"True"},
        "negative conjugate Root compares equal to the matching radical expression");
    tests.expectEqual(eval(session, "root[{-2,0,1},2]<sqrt[3]"), std::string{"True"},
        "Root-versus-radical real ordering is exact");
    tests.expectEqual(eval(session, "root[{-2,0,0,1},1]==cbrt[2]"), std::string{"True"},
        "real cubic Root compares equal to the matching positive cbrt expression");
    tests.expectEqual(eval(session, "root[{2,0,0,1},1]==-cbrt[2]"), std::string{"True"},
        "real cubic Root compares equal to the matching negative cbrt expression");
    tests.expectEqual(eval(session, "Phi==root[{-1,-1,1},2]"), std::string{"True"},
        "algebraic constants participate in the same exact comparison bridge");
    tests.expectEqual(eval(session, "solve[x==root[{-2,0,1},2],x,Rational]"), std::string{"{}"},
        "direct algebraic Root bindings are filtered by Rational domain knowledge");
    tests.expectEqual(eval(session, "solve[x==root[{-2,0,1},2],x,Real]"),
        std::string{"{x == root[{-2, 0, 1}, 2]}"},
        "direct algebraic Root bindings survive compatible Real constraints");
    tests.expectEqual(eval(session, "solve[x==sqrt[2],x,Rational]"), std::string{"{}"},
        "radical bindings reuse algebraic irrationality in Solve constraints");
    tests.expectEqual(eval(session, "solve[x==Phi,x,Rational]"), std::string{"{}"},
        "algebraic constants reuse exact non-Rational knowledge in Solve constraints");

    tests.expectEqual(eval(session, "root[{-4,0,2},2]"),
        std::string{"root[{-2, 0, 1}, 2]"},
        "root canonicalizes its defining polynomial to monic form");
    tests.expectEqual(eval(session, "root[{4,0,-4,0,1},2]"),
        std::string{"root[{-2, 0, 1}, 2]"},
        "root removes repeated polynomial factors because indices count distinct real roots");
    tests.expectEqual(eval(session, "root[{6,0,-5,0,1},1]"),
        std::string{"root[{-3, 0, 1}, 1]"},
        "root reduces a square-free reducible polynomial to the selected real minimal factor");
    tests.expectEqual(eval(session, "root[{6,0,-5,0,1},3]"),
        std::string{"root[{-2, 0, 1}, 2]"},
        "real minimal-polynomial reduction preserves the selected conjugate root");
    tests.expectEqual(eval(session, "root[{2,0,3,0,1},1,Complex]"),
        std::string{"root[{2, 0, 1}, 1, Complex]"},
        "complex Root reduces a reducible defining polynomial to the certified selected factor");
    tests.expectEqual(eval(session, "N[root[{-2,0,1},2],30]"),
        std::string{"1.41421356237309504880168872421"},
        "root has certified arbitrary-precision numerical refinement");
    tests.expectEqual(eval(session, "N[root[{0,-2,0,1},2],30]"),
        std::string{"0.0"},
        "root refinement recognizes an exact zero root without relative-precision stalling");
    tests.expectEqual(eval(session, "solve[x^5-x+1==0,x,Real]"),
        std::string{"{x == root[{1, -1, 0, 0, 0, 1}, 1]}"},
        "Real Solve falls back to an exact Root representation for unresolved rational polynomials");
    tests.expectEqual(eval(session, "solve[(x^2-2)^2==0,x,Real]"),
        std::string{"{x == root[{-2, 0, 1}, 1], x == root[{-2, 0, 1}, 2]}"},
        "algebraic Root fallback preserves distinct repeated-polynomial roots canonically");
    tests.expectEqual(eval(session, "solve[x^4+1==0,x,Real]"),
        std::string{"{}"},
        "Real Root isolation proves that a rational polynomial has no real roots");
    tests.expectEqual(eval(session, "solve[(x^2-2)*(x^2-3)==0,x,Real]"),
        std::string{"{x == root[{-3, 0, 1}, 1], x == root[{-2, 0, 1}, 1], x == root[{-2, 0, 1}, 2], x == root[{-3, 0, 1}, 2]}"},
        "Real Solve canonicalizes reducible polynomial roots to their minimal factors");
    tests.expectEqual(eval(session, "solve[x^5-x+1==0,x]"),
        std::string{"{x == root[{1, -1, 0, 0, 0, 1}, 1, Complex], x == root[{1, -1, 0, 0, 0, 1}, 2, Complex], x == root[{1, -1, 0, 0, 0, 1}, 3, Complex], x == root[{1, -1, 0, 0, 0, 1}, 4, Complex], x == root[{1, -1, 0, 0, 0, 1}, 5, Complex]}"},
        "Complex Solve falls back to certified complex Root isolation for unresolved rational polynomials");
    tests.expectEqual(eval(session, "N[root[{1,0,1},2,Complex],30]"),
        std::string{"1.0I"},
        "complex Root refinement certifies an exact imaginary algebraic root");
    tests.expectEqual(eval(session, "root[{-2,0,1},2]*root[{-2,0,1},2]"),
        std::string{"2"},
        "bounded AlgebraicNumber arithmetic re-identifies an exact rational product");
    tests.expectEqual(eval(session,
        "(root[{-2,0,1},2]+root[{-3,0,1},2])^2"),
        std::string{"root[{1, -10, 1}, 2]"},
        "minimal-polynomial reduction keeps the correct positive conjugate after algebraic squaring");
    tests.expectEqual(eval(session,
        "(root[{-2,0,1},2]+root[{-3,0,1},2])^3"),
        std::string{"root[{1, 0, -970, 0, 1}, 4]"},
        "chained algebraic powers retain the certified result root while reducing the polynomial");
    tests.expectEqual(eval(session,
        "root[{-2,0,1},2]+root[{-3,0,0,1},1]"),
        std::string{"root[{1, -36, 12, -6, -6, 0, 1}, 2]"},
        "primitive-element reduction computes a degree-six minimal polynomial in a linearly disjoint compositum");

    tests.expectEqual(eval(session,
        "(root[{-2,0,1},2]+root[{-3,0,0,1},1])^2"),
        std::string{"root[{1, -1272, -300, -178, 60, -12, 1}, 2]"},
        "persistent number-field coordinates reuse the degree-six compositum for squaring");
    tests.expectEqual(eval(session,
        "(root[{-2,0,1},2]+root[{-3,0,0,1},1])^3"),
        std::string{"root[{1, -45378, -29049, -4140, 111, -18, 1}, 2]"},
        "persistent number-field coordinates survive rational identity factors in binary powering");
    tests.expectEqual(eval(session,
        "(root[{-2,0,1},2]+root[{-3,0,0,1},1])^-1"),
        std::string{"root[{1, 0, -6, -6, 12, -36, 1}, 1]"},
        "persistent number-field division uses an exact inverse in Q(theta)");
    (void)eval(session, "nfStage2A:=root[{-2,0,1},2]+root[{-3,0,0,1},1]");
    (void)eval(session, "nfStage2B:=simplify[nfStage2A]");
    tests.expectEqual(eval(session, "nfStage2B^2"),
        std::string{"root[{1, -1272, -300, -178, 60, -12, 1}, 2]"},
        "persistent number-field cache survives unchanged Simplifier Call reconstruction");

    tests.expectEqual(eval(session, "root[{1,0,1},2,Complex]+I"),
        std::string{"2I"},
        "AlgebraicNumber arithmetic mixes complex Root values with exact complex rationals");
    tests.expectEqual(eval(session, "root[{1,-1,0,0,0,1},1,Complex]+1"),
        std::string{"root[{1, 4, -10, 10, -5, 1}, 1, Complex]"},
        "high-degree algebraic roots support exact rational translation without numerical fallback");

    tests.expectEqual(eval(session, "root[{1,-1,0,0,0,1},1,Complex]^2"),
        std::string{"root[{-1, 1, 0, -2, 0, 1}, 3, Complex]"},
        "persistent generator fields keep same-Root powers inside the original degree-five extension");

    tests.expectEqual(eval(session, "root[{-2,0,1},2]>1"), std::string{"True"},
        "algebraic Stage3 orders a real Root against an exact rational without approximation");
    tests.expectEqual(eval(session, "root[{-2,0,1},1]<-1"), std::string{"True"},
        "algebraic Stage3 proves the sign of a negative real Root exactly");
    tests.expectEqual(eval(session, "root[{-2,0,1},2]==1"), std::string{"False"},
        "algebraic Stage3 proves inequality between distinct exact algebraic values");
    tests.expectEqual(eval(session,
        "root[{-2,0,1},2]!=root[{-3,0,1},2]"), std::string{"True"},
        "algebraic Stage3 rejects roots with distinct proven minimal polynomials");
    tests.expectEqual(eval(session,
        "root[{-2,0,1},1]<root[{-2,0,1},2]"), std::string{"True"},
        "algebraic Stage3 orders conjugate real roots by certified isolation");
    tests.expectEqual(eval(session,
        "root[{1,0,1},1,Complex]!=root[{1,0,1},2,Complex]"), std::string{"True"},
        "algebraic Stage3 equality distinguishes complex roots without inventing an order");
    tests.expectEqual(eval(session,
        "root[{1,0,1},1,Complex]<root[{1,0,1},2,Complex]"),
        std::string{"root[{1, 0, 1}, 1, Complex] < root[{1, 0, 1}, 2, Complex]"},
        "algebraic Stage3 leaves mathematical ordering of complex roots undefined");
    const std::string rootExplanation = eval(session, "explain[root[{-2,0,1},2]]");
    tests.expect(rootExplanation.find("{\"Kind\", \"AlgebraicNumber\"}") != std::string::npos
            && rootExplanation.find("{\"PolynomialDegree\", 2}") != std::string::npos
            && rootExplanation.find("{\"RootIndex\", 2}") != std::string::npos,
        "explain exposes Root as an exact real algebraic number without evaluating it numerically");
    const std::string complexRootExplanation = eval(session, "explain[root[{1,0,1},1,Complex]]");
    tests.expect(complexRootExplanation.find("{\"Domain\", \"Complex\"}") != std::string::npos
            && complexRootExplanation.find("{\"RootIndex\", 1}") != std::string::npos,
        "explain exposes complex Root values as exact complex algebraic numbers");

    tests.expectEqual(eval(session, "solve[sin[x]==0,x,Real]"),
        std::string{"{x == Pi k where k in Integer}"},
        "Solve represents real sine zeros as integer-parameter solution families");
    tests.expectEqual(eval(session, "solve[cos[x]==0,x,Real]"),
        std::string{"{x == Pi/2+Pi k where k in Integer}"},
        "Solve represents real cosine zeros as periodic families");
    tests.expectEqual(eval(session, "solve[tan[x]==1,x,Real]"),
        std::string{"{x == Pi/4+Pi k where k in Integer}"},
        "Solve uses the half-turn period of tangent");
    tests.expectEqual(eval(session, "solve[sin[2x+1]==0,x,Real]"),
        std::string{"{x == -(1-Pi k)/2 where k in Integer}"},
        "Solve propagates periodic targets through an exact affine argument");
    tests.expectEqual(eval(session, "solve[sin[x]==1,x,Real]"),
        std::string{"{x == Pi/2+2Pi k where k in Integer}"},
        "Sine endpoint targets avoid duplicate periodic branches");
    tests.expectEqual(eval(session, "solve[cos[x]==-1,x,Real]"),
        std::string{"{x == Pi+2Pi k where k in Integer}"},
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
        std::string{"{x == 180k where k in Integer}"},
        "periodic Solve respects Degree session semantics");
    tests.expectEqual(eval(degreePeriodic, "solve[tan[x]==1,x,Real]"),
        std::string{"{x == 180k+45 where k in Integer}"},
        "tangent periodic families respect Degree semantics");
}

} // namespace mmcal::tests
