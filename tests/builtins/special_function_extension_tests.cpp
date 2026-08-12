// special・function・extensionの回帰テスト
#include "special_function_extension_tests.hpp"

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

void runSpecialFunctionExtensionTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(eval(session, "log2[8]"), std::string{"3"},
        "log2 reuses exact arbitrary-base logarithm semantics");
    tests.expectEqual(eval(session, "log10[1000]"), std::string{"3"},
        "log10 reuses exact arbitrary-base logarithm semantics");
    tests.expectEqual(eval(session, "D[log2[x],x]"), std::string{"1/(x log[2])"},
        "log2 derivative uses the canonical natural-log denominator");

    tests.expectEqual(eval(session, "gamma[5]"), std::string{"24"},
        "Gamma at a positive integer is exact factorial");
    tests.expectEqual(eval(session, "gamma[1/2]"), std::string{"sqrt[Pi]"},
        "Gamma one-half is exact");
    tests.expectEqual(eval(session, "gamma[-1/2]"), std::string{"-2sqrt[Pi]"},
        "Gamma negative half-integer is exact");
    tests.expect(evalError(session, "gamma[0]").type() == error::CalcErrorType::Domain,
        "Gamma rejects non-positive integer poles");
    tests.expectEqual(eval(session, "N[gamma[1/3],20]"),
        std::string{"2.67893853470774763366"},
        "certified Gamma encloses a general positive real value");
    tests.expectEqual(eval(session, "N[gamma[-1/3],20]"),
        std::string{"-4.06235381827920125084"},
        "certified Gamma uses reflection on the negative real axis");
    tests.expectEqual(eval(session, "N[lgamma[1/3],20]"),
        std::string{"0.98542064692776706919"},
        "lgamma is certified on the real axis");

    tests.expectEqual(eval(session, "erf[0]"), std::string{"0"},
        "erf zero is exact");
    tests.expectEqual(eval(session, "erfc[0]"), std::string{"1"},
        "erfc zero is exact");
    tests.expectEqual(eval(session, "N[erf[1],20]"),
        std::string{"0.84270079294971486934"},
        "erf certified Maclaurin evaluation is accurate");
    tests.expectEqual(eval(session, "N[erfc[1],20]"),
        std::string{"0.15729920705028513066"},
        "erfc certified evaluation is accurate");
    tests.expectEqual(eval(session, "D[erf[x],x]"),
        std::string{"2*exp[-x^2]/sqrt[Pi]"},
        "erf derivative is symbolic and exact");
    tests.expectEqual(eval(session, "D[erfc[x],x]"),
        std::string{"-2*exp[-x^2]/sqrt[Pi]"},
        "erfc derivative is symbolic and exact");

    tests.expectEqual(eval(session, "fresnelc[0]"), std::string{"0"},
        "Fresnel C is exact at zero");
    tests.expectEqual(eval(session, "fresnels[0]"), std::string{"0"},
        "Fresnel S is exact at zero");
    tests.expectEqual(eval(session, "fresnelc[-1]"), std::string{"-fresnelc[1]"},
        "Fresnel C uses exact odd parity");
    tests.expectEqual(eval(session, "N[fresnelc[1],20]"),
        std::string{"0.77989340037682282947"},
        "Fresnel C has a certified real numerical backend");
    tests.expectEqual(eval(session, "N[fresnels[1],20]"),
        std::string{"0.43825914739035476608"},
        "Fresnel S has a certified real numerical backend");
    tests.expectEqual(eval(session, "N[fresnelc[10],20]"),
        std::string{"0.49989869420551572361"},
        "Fresnel C switches safely to the large-argument backend");
    tests.expectEqual(eval(session, "D[fresnelc[x],x]"),
        std::string{"cos[Pi x^2/2 Rad]"},
        "Fresnel C derivative is angle-mode independent and explicitly radian");
    tests.expectEqual(eval(session, "D[fresnels[x],x]"),
        std::string{"sin[Pi x^2/2 Rad]"},
        "Fresnel S derivative is angle-mode independent and explicitly radian");

    tests.expectEqual(eval(session, "hypergeometric1F1[0,3,2]"), std::string{"1"},
        "1F1 with a=0 terminates to one exactly");
    tests.expectEqual(eval(session, "hypergeometric1F1[-2,3,2]"), std::string{"0"},
        "terminating 1F1 polynomial is evaluated exactly");
    tests.expectEqual(eval(session, "hypergeometric1F1[2,2,1]"), std::string{"E"},
        "1F1(a;a;z) reuses exp when the parameter is away from its poles");
    tests.expectEqual(eval(session, "N[hypergeometric1F1[1/6,7/6,1],20]"),
        std::string{"1.19206880798188830082"},
        "1F1 has a certified exact-Rational real series backend");
    tests.expectEqual(eval(session, "D[hypergeometric1F1[1/6,7/6,x],x]"),
        std::string{"hypergeometric1F1[7/6, 13/6, x]/7"},
        "1F1 derivative uses the exact contiguous derivative identity");

    tests.expectEqual(eval(session, "hypergeometric2F1[-2,1,3,1/2]"), std::string{"17/24"},
        "terminating 2F1 polynomial is evaluated exactly");
    tests.expectEqual(eval(session, "hypergeometric2F1[0,2,3,x]"), std::string{"1"},
        "2F1 with a zero numerator parameter reduces before numerical evaluation");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/2,3/2,1/4],20]"),
        std::string{"1.04719755119659774615"},
        "2F1 has a certified real backend inside its principal Gauss-series disk");
    tests.expectEqual(eval(session, "D[hypergeometric2F1[1/2,1/3,5/4,x],x]"),
        std::string{"2hypergeometric2F1[3/2, 4/3, 9/4, x]/15"},
        "2F1 differentiates with respect to its argument when parameters are constant");
    tests.expectEqual(eval(session, "D[hypergeometric2F1[x,1,2,x],x]"),
        std::string{"D[hypergeometric2F1[x, 1, 2, x], x]"},
        "2F1 parameter derivatives remain unresolved instead of applying an incomplete z-only rule");

    tests.expectEqual(eval(session, "ellipticF[x,0]"), std::string{"x"},
        "ellipticF degenerates exactly at m=0");
    tests.expectEqual(eval(session, "ellipticE[x,0]"), std::string{"x"},
        "ellipticE degenerates exactly at m=0");
    tests.expectEqual(eval(session, "ellipticPi[0,x,0]"), std::string{"x"},
        "ellipticPi reduces through ellipticF when n=0");
    tests.expectEqual(eval(session, "N[ellipticF[1/2,1/3],20]"),
        std::string{"0.50684775626543110920"},
        "ellipticF has a certified real-amplitude backend");
    tests.expectEqual(eval(session, "N[ellipticE[1/2,1/3],20]"),
        std::string{"0.49331536201475850521"},
        "ellipticE has a certified real-amplitude backend");
    tests.expectEqual(eval(session, "N[ellipticPi[1/5,1/2,1/3],20]"),
        std::string{"0.51520338216141386085"},
        "ellipticPi has a certified real-amplitude backend away from poles");
    tests.expectEqual(eval(session, "D[ellipticF[x,1/3],x]"),
        std::string{"1/sqrt[1-sin[x Rad]^2/3]"},
        "ellipticF amplitude derivative is exact and explicitly radian");
    tests.expectEqual(eval(session, "D[ellipticE[x,1/3],x]"),
        std::string{"sqrt[1-sin[x Rad]^2/3]"},
        "ellipticE amplitude derivative is exact and explicitly radian");
    tests.expectEqual(eval(session, "D[ellipticPi[1/5,x,1/3],x]"),
        std::string{"1/((1-sin[x Rad]^2/5)sqrt[1-sin[x Rad]^2/3])"},
        "ellipticPi amplitude derivative is exact and explicitly radian");
    tests.expectEqual(eval(session, "D[ellipticF[x,x],x]"),
        std::string{"D[ellipticF[x, x], x]"},
        "elliptic parameter derivatives remain unresolved until their complete formulas are implemented");

    tests.expect(evalError(session, "Ei[0]").type() == error::CalcErrorType::Domain,
        "Ei rejects its logarithmic singularity at zero");
    tests.expectEqual(eval(session, "Si[0]"), std::string{"0"},
        "Si is exact at zero");
    tests.expectEqual(eval(session, "Si[-1]"), std::string{"-Si[1]"},
        "Si uses exact odd parity");
    tests.expect(evalError(session, "Ci[0]").type() == error::CalcErrorType::Domain,
        "Ci rejects its logarithmic singularity at zero");
    tests.expect(evalError(session, "li[1]").type() == error::CalcErrorType::Domain,
        "li rejects its logarithmic singularity at one");
    tests.expectEqual(eval(session, "N[Ei[1],20]"),
        std::string{"1.89511781635593675547"},
        "Ei has a certified real backend near the origin");
    tests.expectEqual(eval(session, "N[Si[1],20]"),
        std::string{"0.94608307036718301494"},
        "Si has a certified real backend");
    tests.expectEqual(eval(session, "N[Ci[1],20]"),
        std::string{"0.33740392290096813466"},
        "Ci has a certified positive-real backend");
    tests.expectEqual(eval(session, "N[li[2],20]"),
        std::string{"1.04516378011749278484"},
        "li reuses the certified log and Ei backends");
    tests.expectEqual(eval(session, "D[Ei[x],x]"), std::string{"exp[x]/x"},
        "Ei derivative is exact");
    tests.expectEqual(eval(session, "D[Si[x],x]"), std::string{"sin[x]/x"},
        "Si derivative uses the default Radian syntax without a redundant unit wrapper");
    tests.expectEqual(eval(session, "D[Ci[x],x]"), std::string{"cos[x]/x"},
        "Ci derivative uses the default Radian syntax without a redundant unit wrapper");
    tests.expectEqual(eval(session, "D[li[x],x]"), std::string{"1/log[x]"},
        "li derivative is exact");

    tests.expectEqual(eval(session, "polylog[0,x]"), std::string{"x/(1-x)"},
        "polylog order zero reduces to a rational function");
    tests.expectEqual(eval(session, "polylog[1,x]"), std::string{"-log[1-x]"},
        "polylog order one reduces to principal Log");
    tests.expectEqual(eval(session, "polylog[2,1]"), std::string{"Pi^2/6"},
        "dilogarithm at one is exact");
    tests.expectEqual(eval(session, "polylog[2,-1]"), std::string{"-Pi^2/12"},
        "dilogarithm at minus one is exact");
    tests.expectEqual(eval(session, "N[polylog[2,1/2],20]"),
        std::string{"0.58224052646501250590"},
        "polylog has a certified |z|<1 real series backend for positive integer order");
    tests.expectEqual(eval(session, "D[polylog[2,x],x]"),
        std::string{"-log[1-x]/x"},
        "dilogarithm derivative reduces Li_1 directly to Log");
    tests.expectEqual(eval(session, "D[polylog[3,x],x]"),
        std::string{"polylog[2, x]/x"},
        "general positive-order polylog derivative lowers the order by one");

    tests.expectEqual(eval(session, "beta[2,3]"), std::string{"1/12"},
        "Beta at positive integers is exact");
    tests.expectEqual(eval(session, "beta[1/2,1/2]"), std::string{"Pi"},
        "exact Gamma values simplify Beta one-half exactly");
    tests.expectEqual(eval(session, "beta[1/2,1]"), std::string{"2"},
        "Beta(x,1)=1/x is exact on the positive-real domain");
    tests.expectEqual(eval(session, "betaln[1/2,1/2]"), std::string{"log[Pi]"},
        "betaln preserves an exact logarithmic result");
    tests.expectEqual(eval(session, "N[beta[1/3,2/3],20]"),
        std::string{"3.62759872846843570119"},
        "Beta has a certified positive-real numerical backend");
    tests.expect(evalError(session, "beta[-1/2,2]").type() == error::CalcErrorType::Domain,
        "current Beta contract rejects non-positive real arguments");

    tests.expectEqual(eval(session, "binom[1/2,2]"), std::string{"-1/8"},
        "generalized binomial with nonnegative integer order is exact");
    tests.expectEqual(eval(session, "fallingfact[5,3]"), std::string{"60"},
        "falling factorial is exact");
    tests.expectEqual(eval(session, "risingfact[5,3]"), std::string{"210"},
        "rising factorial is exact");

    tests.expectEqual(eval(session, "sqrt[Pi]^2"), std::string{"Pi"},
        "squaring a principal square root is branch-safe and exact");
}

} // namespace mmcal::tests
