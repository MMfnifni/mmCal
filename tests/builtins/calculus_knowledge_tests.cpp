// calculus・knowledgeの回帰テスト
#include "calculus_knowledge_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "mathematics/angle.hpp"
#include "solver/real_function_analysis.hpp"
#include "test_framework.hpp"

#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

[[nodiscard]] std::string eval(kernel::KernelSession& session, std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}


[[nodiscard]] std::string intervalText(const solver::RealDomainInterval& interval) {
    std::string result = interval.lowerInclusive ? "[" : "(";
    result += interval.lower ? formatting::formatExpr(*interval.lower) : "-Infinity";
    result += ", ";
    result += interval.upper ? formatting::formatExpr(*interval.upper) : "Infinity";
    result += interval.upperInclusive ? "]" : ")";
    return result;
}

[[nodiscard]] std::string rangeText(const solver::RealValueRange& range) {
    std::string result = range.lowerInclusive ? "[" : "(";
    result += formatting::formatExpr(range.lower);
    result += ", ";
    result += formatting::formatExpr(range.upper);
    result += range.upperInclusive ? "]" : ")";
    return result;
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
    tests.expectEqual(eval(session, "limit[D[abs[x],x],x,0,1]"), std::string{"1"},
        "right-limit assumptions re-materialize a held derivative after abs simplifies");
    tests.expectEqual(eval(session, "limit[D[abs[x],x],x,0,-1]"), std::string{"-1"},
        "left-limit assumptions re-materialize the opposite held derivative branch");
    tests.expectEqual(eval(session, "limit[D[abs[x],{x,2}],x,0,1]"), std::string{"0"},
        "held repeated derivatives are materialized after directional simplification");
    tests.expectEqual(eval(session, "limit[cases[x if x>=0;-x if x<0],x,0]"),
        std::string{"0"},
        "two-sided limits of variable-dependent cases compare separately simplified side branches");
    tests.expectEqual(eval(session, "limit[cases[1 if x>0;2 if x<0],x,0]"),
        std::string{"Indeterminate"},
        "two-sided cases limits report disagreement between proven side limits");
    tests.expectEqual(eval(session, "limit[solve[x==1,x],x,0]"),
        std::string{"solve[x == 1, x]"},
        "limit does not capture a solve-bound variable during direct substitution");
    tests.expectEqual(eval(session, "limit[solve[y==x,y],x,0]"),
        std::string{"solve[y == 0, y]"},
        "limit still substitutes a free parameter inside solve relations");
    tests.expectEqual(eval(session, "limit[solve[{x+y==1,x>0},{x,y}],x,0]"),
        std::string{"solve[{x+y == 1, x > 0}, {x, y}]"},
        "limit protects solve variable arrays and their constraints from capture");
    tests.expectEqual(eval(session, "limit[grad[x^2,{x}],x,0]"),
        std::string{"{0}"},
        "limit materializes gradient before substituting its coordinate variable");
    tests.expectEqual(eval(session, "limit[divergence[{x^2,y},{x,y}],x,0]"),
        std::string{"1"},
        "limit materializes divergence before substituting its coordinate variables");
    tests.expectEqual(eval(session, "limit[curl[{-y,x},{x,y}],x,0]"),
        std::string{"2"},
        "limit materializes curl before substituting its coordinate variables");
    tests.expectEqual(eval(session, "limit[laplacian[x^2+y^2,{x,y}],x,0]"),
        std::string{"4"},
        "limit materializes laplacian before substituting its coordinate variables");
    tests.expectEqual(eval(session, "limit[jacobian[{x^2,x*y},{x,y}],x,0]"),
        std::string{"{{0, 0}, {y, 0}}"},
        "limit materializes jacobian before substituting its coordinate variables");
    tests.expectEqual(eval(session, "limit[hessian[x^2+x*y+y^2,{x,y}],x,0]"),
        std::string{"{{2, 1}, {1, 2}}"},
        "limit materializes hessian before substituting its coordinate variables");
    tests.expectEqual(eval(session, "limit[directionalDerivative[x^2+y^2,{1,0},{x,y}],x,0]"),
        std::string{"0"},
        "limit materializes directional derivatives before substituting coordinates");
    tests.expectEqual(eval(session, "limit[diff[x^2,x,1],x,0]"),
        std::string{"diff[x^2, x, 1]"},
        "limit does not capture numeric-derivative control variables");
    tests.expectEqual(eval(session, "limit[diff[x^2,x,x],x,0]"),
        std::string{"diff[x^2, x, 0]"},
        "limit substitutes free numeric-derivative evaluation points without touching the binder");
    tests.expectEqual(eval(session, "limit[collect[x^2+x,x],x,0]"),
        std::string{"0"},
        "limit materializes collect before substituting its control variable");
    tests.expectEqual(eval(session, "limit[groebnerBasis[{x^2+y,x*y-1},{x,y}],x,0]"),
        std::string{"{y, -1, y^2}"},
        "limit materializes Groebner basis before substituting polynomial-ring variables");
    tests.expectEqual(eval(session, "limit[polynomialReduce[x^2+y,{x+y},{x,y}],x,0]"),
        std::string{"{{-y}, y^2+y}"},
        "limit materializes polynomialReduce before substituting polynomial-ring variables");
    tests.expectEqual(eval(session, "limit[expand[(x+1)^2],x,0]"),
        std::string{"1"},
        "limit materializes expand before direct point substitution");
    tests.expectEqual(eval(session, "limit[factor[x^2-1],x,1]"),
        std::string{"0"},
        "limit materializes factor before direct point substitution");
    tests.expectEqual(eval(session, "limit[normal[series[exp[x],{x,0,2}]],x,0]"),
        std::string{"1"},
        "limit connects held Series through Normal before taking the finite-point limit");
    tests.expectEqual(eval(session, "limit[toNormal[series[exp[x],{x,0,2}]],x,0]"),
        std::string{"1"},
        "limit connects held Series through toNormal before taking the finite-point limit");
    tests.expectEqual(eval(session, "limit[{{x-y},y^2+y},x,0]"),
        std::string{"{{-y}, y^2+y}"},
        "limit evaluates ragged brace values componentwise");

    tests.expectEqual(eval(session, "integrate[D[x^2,x],x]"), std::string{"x^2"},
        "integrate materializes a held derivative before constructing the primitive");
    tests.expectEqual(eval(session, "integrate[normal[series[exp[x],{x,0,2}]],x]"),
        std::string{"x+x^2/2+x^3/6"},
        "integrate materializes Normal[Series] before integration");
    tests.expectEqual(eval(session, "normal[series[D[exp[x],x],{x,0,5}]]"),
        std::string{"1+x+x^2/2+x^3/6+x^4/24+x^5/120"},
        "series materializes a held derivative before coefficient extraction");
    tests.expectEqual(eval(session,
        "normal[series[normal[series[exp[x],{x,0,5}]],{x,0,3}]]"),
        std::string{"1+x+x^2/2+x^3/6"},
        "nested Normal[Series] pipelines compose without leaving held frontends");
    tests.expectEqual(eval(session, "solve[D[x^2,x]==2*x,x,Real]"), std::string{"All"},
        "solve materializes held derivatives before relation normalization");
    tests.expectEqual(eval(session,
        "solve[normal[series[exp[x],{x,0,2}]]==1+x+x^2/2,x,Real]"),
        std::string{"All"},
        "solve materializes Normal[Series] before proving an identity");
    tests.expectEqual(eval(session, "integrate[1+D[x^2,x],x]"),
        std::string{"x^2+x"},
        "integrate materializes a nested held derivative inside ordinary arithmetic");
    tests.expectEqual(eval(session, "integrate[sin[D[x^2,x]],x]"),
        std::string{"-cos[2x]/2"},
        "integrate materializes a held derivative nested inside an elementary function");
    tests.expectEqual(eval(session,
        "integrate[1+normal[series[exp[x],{x,0,2}]],x]"),
        std::string{"2x+x^2/2+x^3/6"},
        "integrate materializes nested Normal[Series] pipelines");
    tests.expectEqual(eval(session,
        "normal[series[1+D[exp[x],x],{x,0,3}]]"),
        std::string{"2+x+x^2/2+x^3/6"},
        "series materializes a nested held derivative inside ordinary arithmetic");
    tests.expectEqual(eval(session,
        "normal[series[1+normal[series[exp[x],{x,0,2}]],{x,0,2}]]"),
        std::string{"2+x+x^2/2"},
        "series materializes nested Normal[Series] inside a larger expression");
    tests.expectEqual(eval(session,
        "D[1+normal[series[exp[x],{x,0,2}]],x]"),
        std::string{"x+1"},
        "D materializes nested Normal[Series] inside ordinary arithmetic");
    tests.expectEqual(eval(session,
        "D[sin[normal[series[x,{x,0,2}]]],x]"),
        std::string{"cos[x]"},
        "D materializes a Series pipeline nested inside an elementary function");
    tests.expectEqual(eval(session, "grad[D[x^2,x],{x}]"),
        std::string{"{2}"},
        "vector calculus materializes a held derivative before differentiating the field");
    tests.expectEqual(eval(session,
        "jacobian[{D[x^2,x],normal[series[exp[y],{y,0,2}]]},{x,y}]"),
        std::string{"{{2, 0}, {0, y+1}}"},
        "vector calculus materializes mixed held derivative and Series pipelines");
    tests.expectEqual(eval(session,
        "laplacian[normal[series[exp[x],{x,0,3}]],{x}]"),
        std::string{"x+1"},
        "laplacian materializes Normal[Series] before repeated differentiation");
    tests.expectEqual(eval(session, "D[expand[(x+1)^2],x]"),
        std::string{"2x+2"},
        "D materializes an explicit expand frontend without evaluating the control variable");
    tests.expectEqual(eval(session, "integrate[factor[x^2-1],x]"),
        std::string{"-x+x^3/3"},
        "integrate materializes an explicit factor frontend before integration");
    tests.expectEqual(eval(session,
        "normal[series[expand[(x+1)^3],{x,0,2}]]"),
        std::string{"3x^2+3x+1"},
        "series materializes an explicit expand frontend before coefficient extraction");
    tests.expectEqual(eval(session,
        "solve[expand[(x+1)^2]==x^2+2*x+1,x,Real]"),
        std::string{"All"},
        "solve materializes explicit algebra transforms before relation normalization");

    kernel::KernelSession heldTransformSession;
    static_cast<void>(eval(heldTransformSession, "x:=5"));
    tests.expectEqual(eval(heldTransformSession, "D[expand[(x+1)^2],x]"),
        std::string{"2x+2"},
        "held algebra-transform materialization does not resolve a session definition of the D variable");
    tests.expectEqual(eval(heldTransformSession, "integrate[expand[(x+1)^2],x]"),
        std::string{"x+x^2+x^3/3"},
        "held algebra-transform materialization does not resolve a session definition of the integration variable");

    tests.expectEqual(eval(session, "D[limit[x*y,x,0],y]"), std::string{"0"},
        "D materializes an inner limit before differentiating");
    tests.expectEqual(eval(session, "D[integrate[x*y,x],y]"), std::string{"x^2/2"},
        "D materializes an inner symbolic integral before differentiating");
    tests.expectEqual(eval(session, "integrate[limit[x*y,x,0],y]"), std::string{"0"},
        "integrate materializes an inner limit before constructing the outer primitive");
    tests.expectEqual(eval(session,
        "normal[series[limit[exp[x*y],x,0],{y,0,3}]]"), std::string{"1"},
        "series materializes an inner limit before coefficient extraction");
    tests.expectEqual(eval(session, "solve[limit[x*y,x,0]==0,y,Real]"), std::string{"All"},
        "solve materializes an inner limit before relation normalization");
    tests.expectEqual(eval(session, "grad[limit[x*y,x,y],{y}]"), std::string{"{2y}"},
        "vector calculus materializes an inner limit before differentiation");

    kernel::KernelSession nestedBinderSession;
    static_cast<void>(eval(nestedBinderSession, "x:=7"));
    static_cast<void>(eval(nestedBinderSession, "y:=5"));
    tests.expectEqual(eval(nestedBinderSession, "D[limit[x*y,x,y],y]"),
        std::string{"2y"},
        "inner limit materialization protects the outer D variable from session bindings");
    tests.expectEqual(eval(nestedBinderSession, "D[integrate[x*y,{x,0,y}],y]"),
        std::string{"3y^2/2"},
        "inner definite integration protects the outer D variable in endpoint expressions");
    tests.expectEqual(eval(session, "limit[integrate[x*y,x],y,0]"), std::string{"0"},
        "limit closes a residual symbolic integral after point substitution");
    tests.expectEqual(eval(session, "limit[1+integrate[x*y,x],y,0]"), std::string{"1"},
        "limit closes residual symbolic integrals inside ordinary arithmetic");
    tests.expectEqual(eval(session, "limit[limit[x*y,x,y],y,0]"), std::string{"0"},
        "limit closes an inner limit whose binder differs from the outer variable");
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
    tests.expectEqual(eval(session, "limit[x*Ei[x],x,0,1]"), std::string{"0"},
        "local Series proves a vanishing factor times logarithmic Ei at zero");
    tests.expectEqual(eval(session, "limit[sqrt[x]*log[x],x,0,1]"), std::string{"0"},
        "local Series proves a positive Puiseux power dominates log at zero");
    tests.expectEqual(eval(session, "limit[Ei[x]-log[x],x,0,1]"), std::string{"-digamma[1]"},
        "local Series resolves finite cancellation between Ei and log singularities");
    tests.expectEqual(eval(session, "limit[Ci[x]-log[x],x,0,1]"), std::string{"-digamma[1]"},
        "local Series resolves finite cancellation between Ci and log singularities");
    tests.expectEqual(eval(session, "limit[log[2*x]-log[x],x,0,1]"), std::string{"log[2]"},
        "local Series resolves logarithmic cancellation with a positive scale factor");
    tests.expectEqual(eval(session, "limit[log[x]^2,x,0,1]"),
        std::string{"limit[log[x]^2, x, 0, 1]"},
        "Series fallback remains conservative when a constant-order logarithmic divergence remains");
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
    tests.expectEqual(eval(session, "limit[Si[x],x,Infinity]"), std::string{"Pi/2"},
        "Si positive-infinity limit uses the exact DLMF constant");
    tests.expectEqual(eval(session, "limit[fresnelc[x],x,Infinity]"), std::string{"1/2"},
        "Fresnel C positive-infinity limit is exact");
    tests.expectEqual(eval(session, "limit[fresnels[x],x,Infinity]"), std::string{"1/2"},
        "Fresnel S positive-infinity limit is exact");
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
    tests.expectEqual(eval(session, "integrate[sin[x]/x,{x,0,Infinity}]"),
        std::string{"Pi/2"},
        "Dirichlet integral closes through the exact Si infinity limit");
    tests.expectEqual(eval(session, "integrate[cos[Pi*x^2/2],{x,0,Infinity}]"),
        std::string{"1/2"},
        "Fresnel cosine integral closes through the exact endpoint limit");
    tests.expectEqual(eval(session, "integrate[sin[Pi*x^2/2],{x,0,Infinity}]"),
        std::string{"1/2"},
        "Fresnel sine integral closes through the exact endpoint limit");
    tests.expectEqual(eval(session, "integrate[exp[-a*x],{x,0,Infinity},a>0]"),
        std::string{"1/a"},
        "assumption-aware Gamma kernels prove the positive exponential scale before reducing");
    tests.expectEqual(eval(session, "integrate[exp[-a*x^2],{x,-Infinity,Infinity},a>0]"),
        std::string{"sqrt[Pi]/sqrt[a]"},
        "assumption-aware Gaussian whole-line integral preserves the positive scale condition");
    tests.expectEqual(eval(session, "integrate[x^(s-1)*exp[-x],{x,0,Infinity},s>0]"),
        std::string{"gamma[s]"},
        "DLMF Gamma integral family closes under a proved positive parameter");
    tests.expectEqual(eval(session, "integrate[x^(a-1)*(1-x)^(b-1),{x,0,1},{a>0,b>0}]"),
        std::string{"beta[a, b]"},
        "Euler Beta integral uses both endpoint convergence assumptions");
    tests.expectEqual(eval(session, "integrate[1/(1+x^4),{x,0,Infinity}]"),
        std::string{"Pi/(2sqrt[2])"},
        "Beta reflection after t=x^q evaluates reciprocal quartic improper integrals");
    tests.expectEqual(eval(session, "integrate[x/(1+x^4),{x,0,Infinity}]"),
        std::string{"Pi/4"},
        "generalized Beta reflection handles monomial numerators");
    tests.expectEqual(eval(session, "integrate[log[x]^2,{x,0,1}]"),
        std::string{"2"},
        "log moments on the unit interval use exact parameter differentiation");
    tests.expectEqual(eval(session, "integrate[x^2*log[x],{x,0,1}]"),
        std::string{"-1/9"},
        "weighted log moments preserve exact endpoint convergence");
    tests.expectEqual(eval(session, "integrate[1/sqrt[1-x^4],{x,0,1}]"),
        std::string{"beta[1/4, 1/2]/4"},
        "endpoint algebraic singularities reduce to Euler Beta after t=x^q");
    tests.expectEqual(eval(session, "integrate[x^n,{x,0,a},{a>0,n>-1}]"),
        std::string{"a^(n+1)/(n+1)"},
        "parameterized power definite integration consumes the exact n>-1 convergence condition");
    tests.expectEqual(eval(session, "integrate[1/(x-a),{x,0,1},{a<0}]"),
        std::string{"log[-(1-a)/a]"},
        "symbolic poles proved below a finite interval are integrated without branch ambiguity");
    tests.expectEqual(eval(session, "integrate[1/(x-a),{x,0,1},{a>1}]"),
        std::string{"log[-(1-a)/a]"},
        "symbolic poles proved above a finite interval use the same positive endpoint ratio");

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
    tests.expectEqual(eval(session, "solve[log[x]==log[x],x]"), std::string{"All if x != 0"},
        "Solve retains the finite domain condition of an identical logarithm equation");
    tests.expectEqual(eval(session, "solve[zeta[x]==zeta[x],x]"), std::string{"All if x != 1"},
        "Solve retains the finite pole exclusion of an identical zeta equation");
    tests.expectEqual(eval(session, "solve[exp[log[x]]==x,x]"), std::string{"All if x != 0"},
        "Solve proves inverse-composition identities under their complete definedness conditions");
    tests.expectEqual(eval(session, "solve[log[exp[x]]==x,x,Real]"), std::string{"All"},
        "Real Solve uses its ambient domain while normalizing principal log-exp composition");
    tests.expectEqual(eval(session, "solve[sqrt[x^2]==x,x,Real]"),
        std::string{"{x in Real if x >= 0}"},
        "Real Solve normalizes sqrt of a square through abs and preserves the positive branch");
    tests.expectEqual(eval(session, "solve[sqrt[x^2]==-x,x,Real]"),
        std::string{"{x in Real if x <= 0}"},
        "Real Solve preserves the negative branch of sqrt[x^2]==-x");
    tests.expectEqual(eval(session, "solve[abs[x]==x,x,Real]"),
        std::string{"{x in Real if x >= 0}"},
        "absolute-value self equality reduces to its exact nonnegative domain");
    tests.expectEqual(eval(session, "solve[abs[x]==-x,x,Real]"),
        std::string{"{x in Real if x <= 0}"},
        "absolute-value negated self equality reduces to its exact nonpositive domain");
    tests.expectEqual(eval(session, "solve[exp[x]==x,x,Real]"), std::string{"{}"},
        "real Exp has no fixed point because exp[u] is strictly above u");
    tests.expectEqual(eval(session, "solve[exp[x+1]==x+1,x,Real]"), std::string{"{}"},
        "real Exp fixed-point exclusion applies to a proved-real composite argument");
    tests.expectEqual(eval(session, "solve[log[x]==x,x,Real]"), std::string{"{}"},
        "principal real Log has no fixed point");
    tests.expectEqual(eval(session, "solve[x+log[x]==0,x,Real]"),
        std::string{"{x == lambertw[1]}"},
        "Log plus its real argument normalizes to the Lambert W normal form");
    tests.expectEqual(eval(session, "solve[exp[x]+x==0,x,Real]"),
        std::string{"{x == -lambertw[1]}"},
        "Exp plus its real argument normalizes through a negated Lambert W target");
    tests.expectEqual(eval(session, "solve[x+1+log[x+1]==0,x,Real]"),
        std::string{"{x == lambertw[1]-1}"},
        "Lambert self-relation normalization composes with affine target solving");
    tests.expectEqual(eval(session, "solve[exp[x+1]+x+1==0,x,Real]"),
        std::string{"{x == -(1+lambertw[1])}"},
        "negated Lambert self-relation normalization composes with affine targets");
    tests.expectEqual(eval(session, "solve[ln[x]==2,x,Real]"), std::string{"{x == exp[2]}"},
        "Solve-safe normalization makes ln share canonical log semantics");
    tests.expectEqual(eval(session, "solve[log2[x]==3,x,Real]"), std::string{"{x == 8}"},
        "Solve inverts canonicalized base-2 logarithms exactly");
    tests.expectEqual(eval(session, "solve[log10[x]==2,x,Real]"), std::string{"{x == 100}"},
        "Solve inverts canonicalized base-10 logarithms exactly");
    tests.expectEqual(eval(session, "solve[log2[x]==y,x,Real]"),
        std::string{"{x == 2^y if y in Real}"},
        "symbolic base-2 logarithm inversion retains the real target condition");
    tests.expectEqual(eval(session, "solve[log10[x]==y,x,Real]"),
        std::string{"{x == 10^y if y in Real}"},
        "symbolic base-10 logarithm inversion retains the real target condition");
    tests.expectEqual(eval(session, "solve[asin[x]==Pi/2,x,Real]"), std::string{"{x == 1}"},
        "principal asin inversion includes its closed upper endpoint");
    tests.expectEqual(eval(session, "solve[asin[x]==-Pi/2,x,Real]"), std::string{"{x == -1}"},
        "principal asin inversion includes its closed lower endpoint without a tautological condition");
    tests.expectEqual(eval(session, "solve[asin[x]==2,x,Real]"), std::string{"{}"},
        "principal asin inversion rejects targets above its active-angle range");
    tests.expectEqual(eval(session, "solve[acos[x]==Pi,x,Real]"), std::string{"{x == -1}"},
        "principal acos inversion includes its closed upper endpoint");
    tests.expectEqual(eval(session, "solve[atan[x]==Pi/2,x,Real]"), std::string{"{}"},
        "principal atan inversion rejects its open endpoint before evaluating tan");
    tests.expectEqual(eval(session, "solve[atan[x]==Pi/4,x,Real]"), std::string{"{x == 1}"},
        "principal atan inversion accepts an interior exact angle");
    tests.expectEqual(eval(session, "solve[acosh[x]==y,x,Real]"),
        std::string{"{x == cosh[y] if y in Real && y >= 0}"},
        "principal acosh inversion retains its nonnegative real range condition");
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
    tests.expectEqual(eval(session, "solve[x*exp[x]==1,x,Real]"),
        std::string{"{x == lambertw[1]}"},
        "canonical u Exp[u] equations normalize directly to Lambert W");
    tests.expectEqual(eval(session, "solve[x*exp[x]==a,x,Real]"),
        std::string{"{x == lambertw[a] if a in Real && a >= -1/E, x == lambertw[-1, a] if a in Real && a > -1/E && a < 0}"},
        "symbolic Lambert normal form preserves both exact real branch conditions");
    tests.expectEqual(eval(session, "solve[exp[-x]==x,x,Real]"),
        std::string{"{x == lambertw[1]}"},
        "negative exponential fixed points normalize through u Exp[u]==1");
    tests.expectEqual(eval(session, "solve[lambertw[x]==1,x]"),
        std::string{"{x == E}"},
        "principal Lambert W inversion uses w Exp[w] on its proved real range");
    tests.expectEqual(eval(session, "solve[cosh[x]==2,x,Real]"),
        std::string{"{x == acosh[2], x == -acosh[2]}"},
        "real cosh inversion returns both even branches");
    tests.expectEqual(eval(session, "solve[cosh[x]==1,x,Real]"),
        std::string{"{x == 0}"},
        "real cosh inversion coalesces the two branches at its minimum");
    tests.expectEqual(eval(session, "solve[{cosh[x]==2,x>0},x,Real]"),
        std::string{"{x == acosh[2]}"},
        "constraint filtering uses exact acosh sign knowledge to discard the negative branch");
    tests.expectEqual(eval(session, "solve[{cosh[x]==2,x<0},x,Real]"),
        std::string{"{x == -acosh[2]}"},
        "constraint filtering keeps only the negative cosh branch below zero");
    tests.expectEqual(eval(session, "solve[erf[x]==0,x,Real]"),
        std::string{"{x == 0}"},
        "real erf uses strict monotonicity to certify its unique zero");
    tests.expectEqual(eval(session, "solve[erf[x^2-1]==0,x,Real]"),
        std::string{"{x == 1, x == -1}"},
        "real erf zero knowledge composes with exact polynomial solving");
    tests.expectEqual(eval(session, "solve[erf[x^2+1]==0,x,Real]"),
        std::string{"{}"},
        "real erf zero knowledge can certify that a composite equation has no solution");
    tests.expectEqual(eval(session, "solve[{sqrt[x]==2,x<10},x,Real]"),
        std::string{"{x == 4}"},
        "relation-array Solve reuses the scalar principal-square-root dispatcher");
    tests.expectEqual(eval(session, "solve[{sin[x]==0,x>0},x,Real]"),
        std::string{"{x == Pi k where k in Integer if Pi k > 0}"},
        "relation-array Solve preserves periodic scalar solutions and applies constraints");
    tests.expectEqual(eval(session, "solve[{exp[x]==1,x>=0},x,Real]"),
        std::string{"{x == 0}"},
        "relation-array Solve reuses scalar transcendental dispatch before filtering constraints");
    {
        const expression::Expr xExpression = session.evaluate("x");
        const auto* infinity = session.symbolRegistry().find("Infinity");
        tests.expect(infinity != nullptr,
            "Real function analysis requires the registered Infinity sentinel");
        if (infinity && xExpression.isSymbol()) {
            const expression::Symbol x = xExpression.asSymbol();
            const auto analyze = [&](std::string_view source) {
                return solver::analyzeRealFunction(
                    session.evaluate(source), x,
                    session.builtinRegistry(), session.mathRegistry(),
                    mathematics::defaultAngleSemantics(), infinity->symbol);
            };

            const auto logarithm = analyze("log[x]");
            tests.expect(logarithm.domainComplete && logarithm.domain.size() == 1,
                "Real function analysis finds the complete positive domain of principal Log");
            if (logarithm.domainComplete && logarithm.domain.size() == 1)
                tests.expectEqual(intervalText(logarithm.domain.front()), std::string{"(0, Infinity)"},
                    "principal Log real domain is the positive half-line");
            tests.expect(logarithm.pieces.size() == 1
                    && logarithm.pieces.front().monotonicity
                        == solver::RealIntervalMonotonicity::Increasing,
                "principal Log is certified increasing on its real domain");
            if (logarithm.pieces.size() == 1 && logarithm.pieces.front().range)
                tests.expectEqual(rangeText(*logarithm.pieces.front().range),
                    std::string{"(-Infinity, Infinity)"},
                    "principal Log endpoint limits certify its full real range");

            const auto squareRoot = analyze("sqrt[x]");
            tests.expect(squareRoot.domainComplete && squareRoot.domain.size() == 1,
                "Real function analysis finds the complete principal square-root domain");
            if (squareRoot.domainComplete && squareRoot.domain.size() == 1)
                tests.expectEqual(intervalText(squareRoot.domain.front()), std::string{"[0, Infinity)"},
                    "principal square-root real domain includes zero");
            tests.expect(squareRoot.pieces.size() == 1
                    && squareRoot.pieces.front().range.has_value(),
                "principal square-root range is certified from endpoint limits and monotonicity");
            if (squareRoot.pieces.size() == 1 && squareRoot.pieces.front().range)
                tests.expectEqual(rangeText(*squareRoot.pieces.front().range),
                    std::string{"[0, Infinity)"},
                    "principal square-root range follows from endpoint limits and monotonicity");

            const auto inverseHyperbolic = analyze("atanh[x]");
            tests.expect(inverseHyperbolic.domainComplete && inverseHyperbolic.domain.size() == 1,
                "Real function analysis respects the principal atanh interval");
            if (inverseHyperbolic.domainComplete && inverseHyperbolic.domain.size() == 1)
                tests.expectEqual(intervalText(inverseHyperbolic.domain.front()),
                    std::string{"(-1, 1)"},
                    "principal atanh real domain is open at both branch points");
            tests.expect(inverseHyperbolic.pieces.size() == 1
                    && inverseHyperbolic.pieces.front().range.has_value(),
                "principal atanh range is certified from branch-point limits");
            if (inverseHyperbolic.pieces.size() == 1 && inverseHyperbolic.pieces.front().range)
                tests.expectEqual(rangeText(*inverseHyperbolic.pieces.front().range),
                    std::string{"(-Infinity, Infinity)"},
                    "principal atanh endpoint limits certify its full real range");

            const auto rational = analyze("1/(x^2-1)");
            tests.expect(rational.domainComplete && rational.domain.size() == 3,
                "rational real domain is split at every exact pole");
            if (rational.domainComplete && rational.domain.size() == 3) {
                tests.expectEqual(intervalText(rational.domain[0]),
                    std::string{"(-Infinity, -1)"},
                    "rational domain keeps the left connected component");
                tests.expectEqual(intervalText(rational.domain[1]),
                    std::string{"(-1, 1)"},
                    "rational domain keeps the middle connected component");
                tests.expectEqual(intervalText(rational.domain[2]),
                    std::string{"(1, Infinity)"},
                    "rational domain keeps the right connected component");
            }

            const auto logarithmicResidual = analyze("log[x]-x+1");
            tests.expect(logarithmicResidual.domainComplete
                    && logarithmicResidual.pieces.size() == 2,
                "critical-point partition splits a half-line real domain exactly");
            if (logarithmicResidual.pieces.size() == 2) {
                tests.expectEqual(intervalText(logarithmicResidual.pieces[0].domain),
                    std::string{"(0, 1]"},
                    "logarithmic residual first monotone piece ends at its exact critical point");
                tests.expectEqual(intervalText(logarithmicResidual.pieces[1].domain),
                    std::string{"[1, Infinity)"},
                    "logarithmic residual second monotone piece starts at its exact critical point");
                tests.expect(logarithmicResidual.pieces[0].monotonicity
                        == solver::RealIntervalMonotonicity::Increasing
                        && logarithmicResidual.pieces[1].monotonicity
                            == solver::RealIntervalMonotonicity::Decreasing,
                    "rational derivative sign charts certify monotonicity on critical-point pieces");
            }
        }
    }

    tests.expectEqual(eval(session, "solve[log[x]-x-1==0,x,Real]"),
        std::string{"{}"},
        "interval endpoint bounds prove nonexistence on a half-line real domain");
    tests.expectEqual(eval(session, "solve[log[x]-x+1==0,x,Real]"),
        std::string{"{x == 1}"},
        "domain decomposition plus a unique exact extremum closes a principal-log equation");
    tests.expectEqual(eval(session, "solve[log[x]-2x+2==0,x,Real]"),
        std::string{"UnresolvedSolutionSet[x]"},
        "interval analysis stays unresolved when another exact root cannot be represented");

    tests.expectEqual(eval(session, "solve[sin[x]==x,x,Real]"),
        std::string{"{x == 0}"},
        "non-strict derivative plus a countable zero set proves strict monotonicity");
    tests.expectEqual(eval(session, "solve[x==sin[x],x,Real]"),
        std::string{"{x == 0}"},
        "non-strict strictness proof is invariant under relation orientation");

    tests.expectEqual(eval(session, "solve[abs[x]<2,x]"),
        std::string{"{x in Real if x > -2 && x < 2}"},
        "absolute-value strict inequality reduces to an exact real polynomial sign chart");
    tests.expectEqual(eval(session, "solve[abs[x-1]>=3,x]"),
        std::string{"{x in Real if x <= -2, x in Real if x >= 4}"},
        "absolute-value exterior inequality preserves both connected components");
    tests.expectEqual(eval(session, "solve[abs[2x-1]<=3,x]"),
        std::string{"{x in Real if x >= -1 && x <= 2}"},
        "affine absolute-value inequality reuses polynomial interval solving");
    tests.expectEqual(eval(session, "solve[abs[x]==2,x,Real]"),
        std::string{"{x == 2, x == -2}"},
        "real absolute-value equality splits into its two exact target equations");
    tests.expectEqual(eval(session, "solve[abs[x]==-2,x,Real]"),
        std::string{"{}"},
        "negative absolute-value target is proved impossible over Real");
    tests.expectEqual(eval(session, "solve[abs[x]!=2,x,Real]"),
        std::string{"{x in Real if x != 2 && x != -2}"},
        "real absolute-value disequality excludes both exact target points");
    tests.expectEqual(eval(session, "solve[abs[x]==2,x]"),
        std::string{"UnresolvedSolutionSet[x]"},
        "absolute-value equality does not collapse the default Complex locus to two real points");
    tests.expectEqual(eval(session, "solve[abs[1/x]>=0,x]"),
        std::string{"All if x != 0"},
        "absolute-value range proof preserves a rational argument domain hole");
    tests.expectEqual(eval(session, "solve[abs[log[x]]>=0,x]"),
        std::string{"All if x != 0"},
        "absolute-value range proof preserves principal logarithm definedness");

    tests.expectEqual(eval(session, "solve[exp[x]+x^2+1==0,x,Real]"),
        std::string{"{}"},
        "Real equation prover excludes a globally positive residual without root search");
    tests.expectEqual(eval(session, "solve[cosh[x]+x^4==0,x,Real]"),
        std::string{"{}"},
        "Real equation prover combines nonnegative even powers with a positive function range");
    tests.expectEqual(eval(session, "solve[x+exp[x]-1==0,x,Real]"),
        std::string{"{x == 0}"},
        "strict monotonicity plus an exact anchor proves a unique real root");
    tests.expectEqual(eval(session, "solve[x+erf[x]==0,x,Real]"),
        std::string{"{x == 0}"},
        "MathRegistry monotonicity facts participate in exact uniqueness proofs");
    tests.expectEqual(eval(session, "solve[exp[x]==x+1,x,Real]"),
        std::string{"{x == 0}"},
        "strict convexity and an exact global minimum prove a tangent unique root");
    tests.expectEqual(eval(session, "solve[exp[x]-x+5==0,x,Real]"),
        std::string{"{}"},
        "strict convexity and a positive exact minimum prove real nonexistence");
    tests.expectEqual(eval(session, "solve[erf[x]==1,x,Real]"),
        std::string{"{}"},
        "open real function range proves that an endpoint target has no finite solution");
    tests.expectEqual(eval(session, "solve[erfc[x]==1,x,Real]"),
        std::string{"{x == 0}"},
        "global injectivity plus an exact anchor closes an inverse-free real equation");
    tests.expectEqual(eval(session, "solve[erfc[x]==0,x,Real]"),
        std::string{"{}"},
        "open complementary-error-function range excludes its limiting endpoint");
    tests.expectEqual(eval(session, "solve[erfc[x]==2,x,Real]"),
        std::string{"{}"},
        "open complementary-error-function range excludes its opposite limiting endpoint");
    tests.expectEqual(eval(session, "solve[exp[x]==x+2,x,Real]"),
        std::string{"{x == -lambertw[-exp[-2]]-2, x == -lambertw[-1, -exp[-2]]-2}"},
        "affine exponential equations close through both real Lambert W branches");
    tests.expectEqual(eval(session, "solve[2^x==x,x,Real]"),
        std::string{"{}"},
        "constant-base affine exponential equations prove branch-point nonexistence exactly");
    tests.expectEqual(eval(session, "solve[(4/3)^x==x,x,Real]"),
        std::string{"{x == -lambertw[-log[4/3]]/log[4/3], x == -lambertw[-1, -log[4/3]]/log[4/3]}"},
        "constant-base affine exponential equations preserve both real Lambert branches");
    tests.expectEqual(eval(session, "solve[(1/2)^x==x,x,Real]"),
        std::string{"{x == -lambertw[-log[1/2]]/log[1/2]}"},
        "constant bases below one produce the unique positive real Lambert branch");
    tests.expectEqual(eval(session, "solve[x^x==1,x,Real]"),
        std::string{"{x == 1}"},
        "principal self-power equation closes exactly at the unit target");
    tests.expectEqual(eval(session, "solve[x^x==2,x,Real]"),
        std::string{"{x == exp[lambertw[log[2]]]}"},
        "principal self-power target above one closes through Lambert W");
    tests.expectEqual(eval(session, "solve[x^x==-1,x,Real]"),
        std::string{"{x == -1}"},
        "principal self-power keeps the unique negative unit target");
    tests.expectEqual(eval(session, "solve[x^x==1/4,x,Real]"),
        std::string{"UnresolvedSolutionSet[x]"},
        "self-power targets inside the unit interval remain unresolved when negative integer roots may occur");
    tests.expectEqual(eval(session, "solve[erf[x]==1/2,x,Real]"),
        std::string{"UnresolvedSolutionSet[x]"},
        "uniqueness without an exact representable root remains unresolved");
    tests.expectEqual(eval(session, "solve[exp[x]-x+5==0,x,Complex]"),
        std::string{"UnresolvedSolutionSet[x]"},
        "real nonexistence proofs do not leak into the complex solve domain");
    tests.expectEqual(eval(session, "solve[sqrt[x]==2,x]"),
        std::string{"{x == 4}"},
        "principal square-root equations invert through an exact squared candidate");
    tests.expectEqual(eval(session, "solve[sqrt[x]==-2,x]"),
        std::string{"{}"},
        "principal square-root inversion rejects the negative real image");
    tests.expectEqual(eval(session, "solve[sqrt[x]==I,x]"),
        std::string{"{x == -1}"},
        "principal square-root inversion accepts the positive imaginary boundary image");
    tests.expectEqual(eval(session, "solve[sqrt[x]==y,x]"),
        std::string{"{x == y^2 if re[y] > 0, x == y^2 if re[y] == 0 && im[y] >= 0}"},
        "symbolic square-root targets retain the complete principal-range condition");
    tests.expectEqual(eval(session, "solve[sqrt[-x]==y,x]"),
        std::string{"{x == -y^2 if re[y] > 0, x == -y^2 if re[y] == 0 && im[y] >= 0}"},
        "symbolic square-root targets preserve the same principal-range condition after inversion");
    tests.expectEqual(eval(session, "solve[sqrt[x]==-I,x]"),
        std::string{"{}"},
        "principal square-root inversion rejects the negative imaginary boundary image");
    tests.expectEqual(eval(session, "solve[sqrt[x]==-2+3I,x]"),
        std::string{"{}"},
        "principal square-root inversion reduces exact complex targets before range testing");
    tests.expectEqual(eval(session, "solve[sqrt[x]==2-3I,x]"),
        std::string{"{x == (2-3I)^2}"},
        "principal square-root inversion accepts exact complex targets in the right half-plane");
    tests.expectEqual(eval(session, "solve[sqrt[x+1]==x-1,x]"),
        std::string{"{x == 3}"},
        "square-root candidate filtering removes roots introduced by squaring");
    tests.expectEqual(eval(session, "solve[cbrt[x]==2,x]"),
        std::string{"{x == 8}"},
        "real cube-root equations invert through an exact cubed candidate");
    tests.expectEqual(eval(session, "solve[cbrt[x+1]==x-1,x]"),
        std::string{"{x == root[{-2, 2, -3, 1}, 1]}"},
        "real cube-root equations keep only the exact real transformed root");

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
    tests.expectEqual(eval(session, "N[lambertw[1+I],20]"),
        std::string{"0.65696606923043640587+0.32545033941341502999I"},
        "Lambert W principal branch has a certified complex backend");
    tests.expectEqual(eval(session, "N[lambertw[2,1],20]"),
        std::string{"-2.4015851048680028842+10.776299516115070898I"},
        "Lambert W accepts certified positive complex branch indices");
    tests.expectEqual(eval(session, "N[lambertw[-3,1],20]"),
        std::string{"-2.8535817554090378072-17.113535539412145913I"},
        "Lambert W accepts certified negative complex branch indices");
    tests.expectEqual(eval(session, "N[lambertw[2,1+I],20]"),
        std::string{"-2.1208839379437137158+11.600137110774577828I"},
        "Lambert W preserves explicit branch selection for complex input");
    tests.expectEqual(eval(session, "N[lambertw[-1/E+I/10^8],20]"),
        std::string{"-0.99983512787429915685+0.00016485400656056139308I"},
        "Lambert W principal branch certifies the upper side of the -1/e branch point");
    tests.expectEqual(eval(session, "N[lambertw[-1/E-I/10^8],20]"),
        std::string{"-0.99983512787429915685-0.00016485400656056139308I"},
        "Lambert W principal branch certifies the lower side of the -1/e branch point");
    tests.expectEqual(eval(session, "N[lambertw[-1,-1/E+I/10^8],20]"),
        std::string{"-1.0001648721257003724-0.00016489025031827418034I"},
        "Lambert W branch -1 follows the local branch above the -1/e cut");
    tests.expectEqual(eval(session, "N[lambertw[1,-1/E-I/10^8],20]"),
        std::string{"-1.0001648721257003724+0.00016489025031827418034I"},
        "Lambert W branch +1 follows the symmetric local branch below the -1/e cut");
    tests.expectEqual(eval(session, "N[lambertw[-1,-1/E-I/10^8],20]"),
        std::string{"-3.0888430122347265671-7.4614892575256769112I"},
        "Lambert W branch -1 does not misconnect to the local branch below the cut");
    tests.expectEqual(eval(session, "N[lambertw[-1/E+1/10^12],20]"),
        std::string{"-0.99999766835783058882"},
        "Lambert W principal real branch uses the local branch-point backend to the right of -1/e");
    tests.expectEqual(eval(session, "N[lambertw[-1,-1/E+1/10^12],20]"),
        std::string{"-1.0000023316457937869"},
        "Lambert W lower real branch uses the local branch-point backend to the right of -1/e");
    tests.expectEqual(eval(session, "N[solve[1.1^x==x^2,x,Real],20]"),
        std::string{"{x == -0.95548727594562198165, x == 1.0513800237472769374, x == 95.716830168405222740}"},
        "N approximates numerically closed SolutionSet binding values while preserving variables");
    tests.expectEqual(eval(session, "D[lambertw[x],x]"),
        std::string{"exp[-lambertw[x]]/(1+lambertw[x])"},
        "Lambert W derivative uses the DLMF form regular at the principal zero");
    tests.expectEqual(eval(session, "D[lambertw[x],{x,2}]"),
        std::string{"(-2-lambertw[x])exp[-2lambertw[x]]/(1+lambertw[x])^3"},
        "Lambert W repeated derivatives use the exact DLMF polynomial recurrence");
    tests.expectEqual(eval(session, "D[polylog[3,x],{x,3}]"),
        std::string{"cases[(x/(1-x)+3log[1-x]+2polylog[2, x])/x^3 if x != 0; 2/9 if x == 0]"},
        "polylog repeated derivatives close directly while preserving the removable value at zero");
    tests.expectEqual(eval(session, "D[exp[-x^2],{x,4}]"),
        std::string{"(16x^4-48x^2+12)exp[-x^2]"},
        "repeated derivatives of quadratic exponentials retain a collected polynomial factor");
    tests.expectEqual(eval(session, "D[exp[x]*cos[x],{x,12}]"),
        std::string{"-64cos[x]exp[x]"},
        "repeated derivatives flatten exact linear combinations instead of accumulating nested product-rule terms");

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
    tests.expectEqual(eval(session,
        "solve[x^6-3x^5-x^4+2x^3+2x^2-2x-1==0,x]"),
        std::string{"{x == root[{-1, -2, 2, 2, -1, -3, 1}, 1, Complex], x == root[{-1, -2, 2, 2, -1, -3, 1}, 2, Complex], x == root[{-1, -2, 2, 2, -1, -3, 1}, 3, Complex], x == root[{-1, -2, 2, 2, -1, -3, 1}, 4, Complex], x == root[{-1, -2, 2, 2, -1, -3, 1}, 5, Complex], x == root[{-1, -2, 2, 2, -1, -3, 1}, 6, Complex]}"},
        "Complex Solve proves irreducibility from intersected modular factor-degree constraints");
    tests.expectEqual(eval(session,
        "solve[(x^3-x-1)*(x^3+x+1)==0,x]"),
        std::string{"{x == root[{1, 1, 0, 1}, 1, Complex], x == root[{-1, -1, 0, 1}, 1, Complex], x == root[{1, 1, 0, 1}, 2, Complex], x == root[{-1, -1, 0, 1}, 2, Complex], x == root[{-1, -1, 0, 1}, 3, Complex], x == root[{1, 1, 0, 1}, 3, Complex]}"},
        "Complex Solve batch-canonicalizes reducible all-root isolation to minimal factors");
    tests.expectEqual(eval(session, "N[root[{1,0,1},2,Complex],30]"),
        std::string{"1.0I"},
        "complex Root refinement certifies an exact imaginary algebraic root");
    (void)eval(session, "root[{-2,0,0,0,0,0,0,0,1},1,Complex]");
    tests.expectEqual(eval(session, "N[Out[-1],30]"),
        std::string{"0.0-1.09050773266525765920701065576I"},
        "complex Root numerical evaluation reuses the cached isolating disk from exact evaluation");
    tests.expectEqual(eval(session, "N[root[{-2,0,0,0,0,0,0,0,0,0,1},1,Complex],30]"),
        std::string{"-0.331196214043795628507030588621-1.01931713553736126627822937193I"},
        "complex Root isolation remains stable for sparse symmetric degree-ten polynomials");
    tests.expectEqual(eval(session, "root[{-2,0,1},2]*root[{-2,0,1},2]"),
        std::string{"2"},
        "bounded AlgebraicNumber arithmetic re-identifies an exact rational product");
    tests.expectEqual(eval(session, "root[{-2,0,1},2]+root[{-3,0,1},2]"),
        std::string{"root[{1, 0, -10, 0, 1}, 4]"},
        "pure quadratic addition keeps the exact positive conjugate without general primitive-element construction");
    tests.expectEqual(eval(session, "root[{-2,0,1},2]-root[{-3,0,1},2]"),
        std::string{"root[{1, 0, -10, 0, 1}, 2]"},
        "pure quadratic subtraction keeps the exact negative conjugate without general primitive-element construction");
    tests.expectEqual(eval(session, "root[{-2,0,1},2]*root[{-3,0,1},2]"),
        std::string{"root[{-6, 0, 1}, 2]"},
        "pure quadratic multiplication avoids a general resultant");
    tests.expectEqual(eval(session, "root[{-2,0,1},2]/root[{-3,0,1},2]"),
        std::string{"root[{-2/3, 0, 1}, 2]"},
        "pure quadratic division avoids a general resultant");
    tests.expectEqual(eval(session, "root[{-2,0,1},1]+root[{-3,0,1},2]"),
        std::string{"root[{1, 0, -10, 0, 1}, 3]"},
        "pure quadratic addition preserves the mixed-sign conjugate index");
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
