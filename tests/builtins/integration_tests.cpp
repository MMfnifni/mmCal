// 記号積分integrateの回帰テスト
#include "integration_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "mathematics/definedness.hpp"
#include "simplification/full_simplifier.hpp"
#include "simplification/simplification_context.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace mmcal::tests {
namespace {

[[nodiscard]] std::string eval(kernel::KernelSession& session, std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

// primitiveとintegrandは，双方が定義される点上で一致すればよい。
// 公開fullSimplifyのdomain-hole保持とは分離して，test proofだけ定義条件を仮定する。
[[nodiscard]] std::string derivativeBackProof(
    kernel::KernelSession& session,
    const std::string& primitive,
    std::string_view integrand) {
    expression::Expr residual = session.evaluate(
        std::string{"D[("} + primitive + "),x]-(" + std::string{integrand} + ")");

    mathematics::AssumptionSet assumptions;
    if (auto conditions = mathematics::expressionDomainConditions(
            residual, session.builtinRegistry(), session.mathRegistry()))
        assumptions = std::move(*conditions);

    const mathematics::AngleSemantics angles{session.defaultAngleUnit()};
    simplification::SimplificationContext context{
        session.builtinRegistry(), session.mathRegistry(), angles, std::move(assumptions)};
    context.assumeExpressionsDefined = true;
    return formatting::formatExpr(simplification::fullSimplify(residual, context));
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
    tests.expectEqual(eval(session, "integrate[2/log[3*x+1],x]"), std::string{"2li[3x+1]/3"},
        "affine reciprocal Log integrates directly through li");
    tests.expectEqual(eval(session, "integrate[exp[x],x]"), std::string{"exp[x]"},
        "integrates exp");
    tests.expectEqual(eval(session, "integrate[sin[x],x]"), std::string{"-cos[x]"},
        "integrates sine in the default radian semantics");
    tests.expectEqual(eval(session, "integrate[cos[x],x]"), std::string{"sin[x]"},
        "integrates cosine");
    tests.expectEqual(eval(session, "integrate[2*x*cos[x^2],x]"), std::string{"sin[x^2]"},
        "reverse chain rule is discovered and verified through D");
    tests.expectEqual(eval(session, "integrate[2*x/(x^2+1),x]"), std::string{"log[x^2+1]"},
        "logarithmic substitution is discovered and verified");
    tests.expectEqual(eval(session, "integrate[1/(1+x^2),x]"), std::string{"atan[x]"},
        "inverse tangent derivative knowledge is reusable by integration");
    tests.expectEqual(eval(session, "integrate[log[x],x]"), std::string{"x log[x]-x"},
        "integrates Log by the standard primitive");
    tests.expectEqual(eval(session, "integrate[x*exp[x],x]"), std::string{"x exp[x]-exp[x]"},
        "integration by parts handles polynomial times exp");
    tests.expectEqual(eval(session, "integrate[erf[x],x]"),
        std::string{"exp[-x^2]/sqrt[Pi]+x erf[x]"},
        "integrates erf using shared special-function derivative knowledge");
    tests.expectEqual(eval(session, "integrate[sqrt[x],x]"),
        std::string{"2x sqrt[x]/3"},
        "integrates principal square root when D verifies the result");

    tests.expectEqual(eval(session, "integrate[(2*x+3)^5,x]"),
        std::string{"(2x+3)^6/12"},
        "affine power rule keeps a compact exact primitive");
    tests.expectEqual(eval(session, "integrate[1/(2*x+3),x]"),
        std::string{"log[2x+3]/2"},
        "affine reciprocal integrates to principal Log");
    tests.expectEqual(eval(session, "integrate[exp[2*x+1],x]"),
        std::string{"exp[2x+1]/2"},
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
        std::string{"x^2log[x]/2-x^2/4"},
        "monomial times Log uses the exact integration-by-parts formula");
    tests.expectEqual(eval(session, "integrate[acosh[x],x]"),
        std::string{"x acosh[x]-sqrt[x+1]sqrt[x-1]"},
        "principal acosh primitive uses the same branch convention as D");
    tests.expectEqual(eval(session, "integrate[sin[x Deg],x]"),
        std::string{"-180cos[x Deg]/Pi"},
        "explicit Degree suffix is integrated with the correct radian scale");

    tests.expectEqual(eval(session, "integrate[sin[x]^2,x]"),
        std::string{"(x-sin[2x]/2)/2"},
        "trigonometric square uses an exact double-angle identity");
    tests.expectEqual(eval(session, "integrate[sin[2x]^6,x]"),
        std::string{"5x/16-sin[12x]/384-15sin[4x]/128+3sin[8x]/128"},
        "even trigonometric powers reduce to a finite Fourier polynomial");
    tests.expectEqual(eval(session, "integrate[sin[2x]^(-2),x]"),
        std::string{"-cot[2x]/2"},
        "negative sine powers reuse the reciprocal-trigonometric reduction knowledge");
    tests.expectEqual(eval(session, "integrate[cos[3x]^(-2),x]"),
        std::string{"tan[3x]/3"},
        "negative cosine powers reuse the secant recurrence");
    tests.expectEqual(eval(session, "integrate[sin[x]^5*cos[x]^4,x]"),
        std::string{"-3cos[x]/128-cos[3x]/192+cos[5x]/320+cos[7x]/1792-cos[9x]/2304"},
        "mixed integer sine/cosine powers share the same finite Fourier reduction");
    tests.expectEqual(eval(session, "integrate[sin[2x]*cos[3x],x]"),
        std::string{"(cos[x]-cos[5x]/5)/2"},
        "different trigonometric frequencies fall back to product-to-sum");
    tests.expectEqual(eval(session, "integrate[tan[x]^2,x]"),
        std::string{"tan[x]-x"},
        "tangent square reduces through secant squared");
    tests.expectEqual(eval(session, "integrate[sin[x]*cos[x],x]"),
        std::string{"-cos[2x]/4"},
        "trigonometric product-to-sum keeps the primitive flat");
    tests.expectEqual(eval(session, "integrate[2*x/sqrt[1+x^2],x]"),
        std::string{"2sqrt[x^2+1]"},
        "reverse-chain candidate set recognizes square-root derivatives");
    tests.expectEqual(eval(session, "integrate[2*x*(1+x^2)^5,x]"),
        std::string{"(x^2+1)^6/6"},
        "power-chain Knowledge keeps f-prime times f-to-a-power compact");
    tests.expectEqual(eval(session, "integrate[3*x^2*sqrt[1+x^3],x]"),
        std::string{"2(x^3+1)^(3/2)/3"},
        "power-chain Knowledge recognizes square-root compositions without relying on D expansion");
    tests.expectEqual(eval(session, "integrate[x*sqrt[1+x^2],x]"),
        std::string{"(x^2+1)^(3/2)/3"},
        "power-chain Knowledge determines the proportional derivative factor exactly");
    tests.expectEqual(eval(session, "integrate[sqrt[x]/(1+x),x]"),
        std::string{"-2atan[sqrt[x]]+2sqrt[x]"},
        "principal square-root substitution handles algebraic rational kernels");
    tests.expectEqual(eval(session, "integrate[x/(1+x^4),x]"),
        std::string{"atan[x^2]/2"},
        "generated monomial substitution candidates discover u=x squared");
    tests.expectEqual(eval(session, "integrate[x/sqrt[1-x^4],x]"),
        std::string{"asin[x^2]/2"},
        "generated monomial substitution candidates compose with inverse-trig Knowledge");
    tests.expectEqual(eval(session, "integrate[1/(x^2+4),x]"),
        std::string{"atan[x/2]/2"},
        "exact quadratic reciprocal completes the square to atan");
    tests.expectEqual(eval(session, "integrate[1/(4-x^2),x]"),
        std::string{"atanh[x/2]/2"},
        "exact quadratic reciprocal completes the square to atanh");
    tests.expectEqual(eval(session, "integrate[1/(x+1)^2,x]"),
        std::string{"-1/(x+1)"},
        "repeated linear quadratic uses the rational primitive");
    tests.expectEqual(eval(session, "integrate[1/sqrt[4-x^2],x]"),
        std::string{"asin[x/2]"},
        "negative-leading quadratic root uses asin when the positive scale is provable");
    tests.expectEqual(eval(session, "integrate[1/sqrt[x^2+4],x]"),
        std::string{"asinh[x/2]"},
        "positive-definite quadratic root uses asinh");
    tests.expectEqual(eval(session, "integrate[sqrt[4-x^2],x]"),
        std::string{"x sqrt[4-x^2]/2+2asin[x/2]"},
        "quadratic square-root primitive handles circular completion");
    tests.expectEqual(eval(session, "fullSimplify[D[integrate[sqrt[x^2+4],x],x]-sqrt[x^2+4]]"),
        std::string{"0"},
        "quadratic square-root hyperbolic primitive differentiates back exactly");
    tests.expectEqual(eval(session, "fullSimplify[D[integrate[sqrt[x^2-4],x],x]-sqrt[x^2-4]]"),
        std::string{"0"},
        "quadratic square-root logarithmic primitive differentiates back exactly");
    tests.expectEqual(eval(session, "integrate[log[2,x],x]"),
        std::string{"(x log[x]-x)/log[2]"},
        "arbitrary-base logarithm reuses the natural-log primitive");
    tests.expectEqual(eval(session, "integrate[cos[4x^2],x]"),
        std::string{"fresnelc[x sqrt[8/Pi]]/sqrt[8/Pi]"},
        "quadratic cosine phase closes through Fresnel C");
    tests.expectEqual(eval(session, "integrate[sin[8x^2],x]"),
        std::string{"fresnels[x sqrt[16/Pi]]/sqrt[16/Pi]"},
        "quadratic sine phase closes through Fresnel S");
    tests.expectEqual(eval(session, "integrate[cos[Pi*x^2/2],x]"),
        std::string{"fresnelc[x]"},
        "Fresnel C defining kernel recognizes a symbolic Pi scale");
    tests.expectEqual(eval(session, "integrate[sin[Pi*x^2/2],x]"),
        std::string{"fresnels[x]"},
        "Fresnel S defining kernel recognizes a symbolic Pi scale");
    tests.expectEqual(eval(session, "integrate[sin[2x^2]^4,x]"),
        std::string{"3x/8+fresnelc[x sqrt[16/Pi]]/(8sqrt[16/Pi])-fresnelc[x sqrt[8/Pi]]/sqrt[8/Pi]/2"},
        "trigonometric power reduction composes with Fresnel integration");
    tests.expectEqual(eval(session, "integrate[fresnelc[x],x]"),
        std::string{"x fresnelc[x]-sin[Pi x^2/2 Rad]/Pi"},
        "Fresnel C itself has an exact primitive from shared derivative knowledge");

    tests.expectEqual(eval(session, "integrate[(x+1)/(x+2),x]"),
        std::string{"x-log[x+2]"},
        "rational function with a linear denominator uses exact polynomial division");
    tests.expectEqual(eval(session, "integrate[(x+1)/(x^2+4),x]"),
        std::string{"atan[x/2]/2+log[x^2+4]/2"},
        "rational function over a quadratic splits into log and reciprocal parts");
    tests.expectEqual(derivativeBackProof(
            session, eval(session, "integrate[1/(x*(x+1)),x]"), "1/(x*(x+1))"),
        std::string{"0"},
        "explicit multiplication syntax exercises exact partial fractions without confusing x(...) with a function call");

    tests.expectEqual(eval(session, "integrate[log[1+x]/x,x]"),
        std::string{"-polylog[2, -x]"},
        "dilogarithm Knowledge connects the plus-sign logarithmic kernel");
    tests.expectEqual(eval(session, "integrate[log[x^2+1]/x,x]"),
        std::string{"-polylog[2, -x^2]/2"},
        "dilogarithm Knowledge generalizes to monomial arguments");
    tests.expectEqual(eval(session, "integrate[exp[-x^4],x]"),
        std::string{"x hypergeometric1F1[1/4, 5/4, -x^4]"},
        "1F1 monomial-exponential reduction accepts a negative coefficient");
    tests.expectEqual(eval(session,
        "fullSimplify[D[integrate[exp[x]*(sin[x]+cos[x]),x],x]-exp[x]*(sin[x]+cos[x])]"),
        std::string{"0"},
        "linearity distributes a product over a short sum before exp-trig integration");

    tests.expectEqual(eval(session, "integrate[1/(1+cos[x]),x]"),
        std::string{"tan[x/2]"},
        "Weierstrass substitution reduces a rational cosine kernel");
    tests.expectEqual(eval(session, "integrate[1/(1-cos[x]),x]"),
        std::string{"-1/tan[x/2]"},
        "Weierstrass substitution reduces the complementary cosine kernel");
    tests.expectEqual(eval(session,
        "fullSimplify[D[integrate[1/(1+sin[x]),x],x]-1/(1+sin[x]),"
        "{cos[x/2]!=0,1+tan[x/2]!=0,1+sin[x]!=0}]"),
        std::string{"0"},
        "Weierstrass substitution differentiates back on the common defined domain for the positive sine kernel");
    tests.expectEqual(eval(session,
        "fullSimplify[D[integrate[1/(1-sin[x]),x],x]-1/(1-sin[x]),"
        "{cos[x/2]!=0,-1+tan[x/2]!=0,1-sin[x]!=0}]"),
        std::string{"0"},
        "Weierstrass substitution differentiates back on the common defined domain for the negative sine kernel");
    tests.expect(eval(session, "integrate[1/(2+cos[x]),x]").find("integrate[") == std::string::npos,
        "Weierstrass substitution is a rational transformation rather than four hard-coded identities");

    kernel::KernelSession degreeIntegrals;
    degreeIntegrals.setDefaultAngleUnit(mathematics::AngleUnit::Degree);
    tests.expectEqual(eval(degreeIntegrals, "integrate[1/(1+x^2),x]"),
        std::string{"Pi atan[x]/180"},
        "inverse-trig primitives respect a Degree session without changing the integrand");

    tests.expectEqual(eval(session, "integrate[x^2,{x,0,1}]"), std::string{"1/3"},
        "exact definite polynomial integral");
    tests.expectEqual(eval(session, "integrate[sin[x],{x,0,Pi}]"), std::string{"2"},
        "exact definite trigonometric integral");
    tests.expectEqual(eval(session, "integrate[1/x,{x,1,2}]"), std::string{"log[2]"},
        "definite integral with a certified nonzero denominator");

    tests.expectEqual(eval(session, "integrate[log[x],{x,1,Pi}]"),
        std::string{"1-Pi+Pi log[Pi]"},
        "certified transcendental bounds permit safe exact definite evaluation");

    const std::string pole = eval(session, "integrate[tan[x],{x,0,2}]");
    tests.expect(pole.find("integrate[") == 0,
        "definite integration refuses an interval crossing a trig pole");
    tests.expect(findDiagnostic(session, "integrate::conditionsRequired") != nullptr,
        "unsafe trig definite integral emits a warning");

    const std::string singular = eval(session, "integrate[1/x,{x,-1,1}]");
    tests.expect(singular.find("integrate[") == 0,
        "integral crossing a pole stays unevaluated");
    tests.expect(findDiagnostic(session, "integrate::conditionsRequired") != nullptr,
        "unresolved definite integral emits a warning");

    tests.expectEqual(eval(session, "integrate[1/sqrt[x^2-1],x]"),
        std::string{"log[x+sqrt[x^2-1]]"},
        "branch-sensitive quadratic root uses a verified local primitive without a global sqrt identity");

    // 上側不完全Gammaによる局所表示はbranchとx=0のremovable holeを持つため採用しない。
    // 1F1基盤導入後は原点でentireな形を安全な標準形として使う。
    tests.expectEqual(eval(session, "integrate[exp[x^6],x]"),
        std::string{"x hypergeometric1F1[1/6, 7/6, x^6]"},
        "exponential monomial uses the branch-safe confluent hypergeometric primitive");

    const std::string unsupported = eval(session, "integrate[gamma[x],x]");
    tests.expect(unsupported.find("integrate[") == 0,
        "unsupported special-function antiderivative stays symbolic");
    tests.expect(findDiagnostic(session, "integrate::unsupported") != nullptr,
        "unsupported indefinite integral emits a warning");

    const std::string noClosedForm = eval(session, "integrate[x^x,x]");
    tests.expect(noClosedForm.find("integrate[") == 0,
        "recognized non-closed-form family stays symbolic");
    tests.expect(findDiagnostic(session, "integrate::noKnownClosedForm") != nullptr,
        "recognized non-closed-form family is distinguished from an implementation gap");

    kernel::KernelSession domainDiagnostic;
    const std::string absWithoutDomain = eval(domainDiagnostic, "integrate[abs[x],x]");
    tests.expect(absWithoutDomain.find("integrate[") == 0,
        "absolute-value integral stays symbolic without a real-domain assumption");
    tests.expect(findDiagnostic(domainDiagnostic, "integrate::conditionsRequired") != nullptr,
        "domain-dependent integral reports that assumptions are required");

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
