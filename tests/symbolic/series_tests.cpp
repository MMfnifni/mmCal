// SeriesData / series / normal の回帰テスト
#include "series_tests.hpp"

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

} // namespace

void runSeriesTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(
        eval(session, "series[x,{x,0,3}]"),
        std::string{"seriesData[x, 0, {1, 0, 0}, 1, 4, 1]"},
        "Series: minimumExponent records the first nonzero Taylor exponent");

    tests.expectEqual(
        eval(session, "series[(1+x)^3,{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 3, 3, 1, 0, 0}, 0, 6, 1]"},
        "Series: TPSA expands a polynomial power without repeated differentiation");

    tests.expectEqual(
        eval(session, "normal[series[(1+x)^3,{x,0,5}]]"),
        std::string{"x^3+3x^2+3x+1"},
        "Series: normal removes the order term and returns the truncated exact expression");

    tests.expectEqual(
        eval(session, "series[(x+1)*(x-2),{x,1,3}]"),
        std::string{"seriesData[x, 1, {-2, 1, 1, 0}, 0, 4, 1]"},
        "Series: nonzero expansion centers are represented exactly");

    tests.expectEqual(
        eval(session, "normal[series[(x+1)*(x-2),{x,1,3}]]"),
        std::string{"-3+x+(x-1)^2"},
        "Series: normal reconstructs powers of x-center exactly");

    tests.expectEqual(
        eval(session, "series[1/(1-x),{x,0,6}]"),
        std::string{"seriesData[x, 0, {1, 1, 1, 1, 1, 1, 1}, 0, 7, 1]"},
        "Series: formal inversion expands a regular reciprocal");

    tests.expectEqual(
        eval(session, "series[1/x,{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, 0, 0, 0, 0, 0}, -1, 5, 1]"},
        "Series: Laurent data represents a simple pole exactly");

    tests.expectEqual(
        eval(session, "normal[series[1/x,{x,0,4}]]"),
        std::string{"x^(-1)"},
        "Series: normal reconstructs negative integer powers");

    tests.expectEqual(
        eval(session, "series[(1+x)^-2,{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, -2, 3, -4, 5, -6}, 0, 6, 1]"},
        "Series: negative integer powers use formal inversion rather than differentiation");

    tests.expectEqual(
        eval(session, "series[(1+x)^-4097,{x,0,3}]"),
        std::string{"seriesData[x, 0, {1, -4097, 8394753, -11470030849}, 0, 4, 1]"},
        "Series: analytic negative integer powers are not cut off at the former magnitude 4096 boundary");
    tests.expectEqual(
        eval(session, "series[(1+x)^-1000000,{x,0,2}]"),
        std::string{"seriesData[x, 0, {1, -1000000, 500000500000}, 0, 3, 1]"},
        "Series: truncated integer powers scale logarithmically in the exponent at fixed order");

    tests.expectEqual(
        eval(session, "series[1/(x^2+x^3),{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, -1, 1, -1, 1, -1, 1}, -2, 5, 1]"},
        "Series: Laurent inversion propagates a higher-order pole");

    tests.expectEqual(
        eval(session, "series[(x-1)^-2,{x,1,3}]"),
        std::string{"seriesData[x, 1, {1, 0, 0, 0, 0, 0}, -2, 4, 1]"},
        "Series: shifted Laurent poles use the expansion center exactly");

    tests.expectEqual(
        eval(session, "series[1/(a+x),{x,0,3},a!=0]"),
        std::string{"seriesData[x, 0, {1/a, -1/a^2, 1/(a a^2), -1/a^4}, 0, 4, 1]"},
        "Series: symbolic leading coefficients are inverted only under a nonzero assumption");

    tests.expectEqual(
        eval(session, "series[1/(a+x),{x,0,3}]"),
        std::string{"series[1/(a+x), {x, 0, 3}]"},
        "Series: an unproved symbolic leading coefficient remains unevaluated");
    tests.expect(!session.diagnostics().empty()
            && session.diagnostics().back().code == "series::unsupported",
        "Series: uncertain symbolic inversion emits an explicit diagnostic");

    tests.expectEqual(
        eval(session, "series[exp[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120}, 0, 6, 1]"},
        "Series: exp uses a TPSA coefficient recurrence");

    tests.expectEqual(
        eval(session, "series[log[1+x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, -1/2, 1/3, -1/4, 1/5}, 1, 6, 1]"},
        "Series: log uses A'/A without repeated differentiation");

    tests.expectEqual(
        eval(session, "series[sin[x],{x,0,7}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/6, 0, 1/120, 0, -1/5040}, 1, 8, 1]"},
        "Series: sin and cos use a coupled TPSA recurrence");

    tests.expectEqual(
        eval(session, "series[cos[x],{x,Pi/2,4}]"),
        std::string{"seriesData[x, Pi/2, {-1, 0, 1/6, 0}, 1, 5, 1]"},
        "Series: trigonometric recurrence supports exact nonzero centers");

    tests.expectEqual(
        eval(session, "series[sinh[x],{x,0,7}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/6, 0, 1/120, 0, 1/5040}, 1, 8, 1]"},
        "Series: sinh and cosh use a coupled TPSA recurrence");

    tests.expectEqual(
        eval(session, "series[exp[sin[x]],{x,0,6}]"),
        std::string{"seriesData[x, 0, {1, 1, 1/2, 0, -1/8, -1/15, -1/240}, 0, 7, 1]"},
        "Series: analytic TPSA expansions compose without higher derivatives");

    tests.expectEqual(
        eval(session, "series[1/sin[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/6, 0, 7/360, 0, 31/15120}, -1, 6, 1]"},
        "Series: analytic valuation composes with Laurent inversion");

    tests.expectEqual(
        eval(session, "series[1/log[1+x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 1/2, -1/12, 1/24, -19/720, 3/160, -863/60480}, -1, 6, 1]"},
        "Series: logarithmic valuation composes with Laurent inversion");

    tests.expectEqual(
        eval(session, "series[sqrt[1+x],{x,0,6}]"),
        std::string{"seriesData[x, 0, {1, 1/2, -1/8, 1/16, -5/128, 7/256, -21/1024}, 0, 7, 1]"},
        "Series: principal sqrt uses the rational-power TPSA recurrence");

    tests.expectEqual(
        eval(session, "series[(1+x)^(3/2),{x,0,6}]"),
        std::string{"seriesData[x, 0, {1, 3/2, 3/8, -1/16, 3/128, -3/256, 7/1024}, 0, 7, 1]"},
        "Series: exact rational powers expand directly without exp-log rewriting");

    tests.expectEqual(
        eval(session, "series[(2+x)^(-1/3),{x,0,5}]"),
        std::string{"seriesData[x, 0, {2^(-1/3), -2^(-1/3)/6, 2^(-1/3)/18, -7*2^(-1/3)/324, 35*2^(-1/3)/3888, -91*2^(-1/3)/23328}, 0, 6, 1]"},
        "Series: negative rational powers retain an exact principal constant term");

    tests.expectEqual(
        eval(session, "series[1/sqrt[1+x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, -1/2, 3/8, -5/16, 35/128, -63/256}, 0, 6, 1]"},
        "Series: principal sqrt composes with Laurent inversion");

    tests.expectEqual(
        eval(session, "series[sqrt[x^2],{x,0,5}]"),
        std::string{"series[sqrt[x^2], {x, 0, 5}]"},
        "Series: a higher-multiplicity principal branch point is not collapsed to one local branch");

    tests.expectEqual(
        eval(session, "series[sqrt[a+x],{x,0,3},a>0]"),
        std::string{"seriesData[x, 0, {sqrt[a], sqrt[a]/(2a), -sqrt[a]/(8a^2), sqrt[a]/(16a a^2)}, 0, 4, 1]"},
        "Series: assumptions certify a symbolic principal-power expansion center");

    tests.expectEqual(
        eval(session, "series[log[x],{x,-1,3}]"),
        std::string{"series[log[x], {x, -1, 3}]"},
        "Series: principal log refuses a local Taylor series on the negative-real branch cut");

    tests.expectEqual(
        eval(session, "series[log[x],{x,I,3}]"),
        std::string{"seriesData[x, I, {I Pi/2, -I, 1/2, I/3}, 0, 4, 1]"},
        "Series: principal log expands at a certified nonreal center");

    kernel::KernelSession degreeSession;
    static_cast<void>(degreeSession.evaluate("angleMode[Deg]"));
    tests.expectEqual(
        eval(degreeSession, "series[sin[x],{x,0,1}]"),
        std::string{"seriesData[x, 0, {Pi/180}, 1, 2, 1]"},
        "Series: direct trigonometric recurrence honors the session angle mode");
    tests.expectEqual(
        eval(degreeSession, "series[sin[x Rad],{x,0,3}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/6}, 1, 4, 1]"},
        "Series: an explicit angle unit overrides the session angle mode");

    tests.expectEqual(
        eval(session, "series[sqrt[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 1, 11, 2]"},
        "Series: branch-point sqrt uses a half-integer Puiseux grid");

    tests.expectEqual(
        eval(session, "series[x^(1/3),{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 1, 13, 3]"},
        "Series: exact rational powers introduce the required Puiseux denominator");

    tests.expectEqual(
        eval(session, "series[sqrt[x]*(1+x),{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, 0, 1, 0, 0, 0, 0, 0}, 1, 9, 2]"},
        "Series: Puiseux grids compose with ordinary Taylor factors");

    tests.expectEqual(
        eval(session, "series[1/sqrt[x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, -1, 9, 2]"},
        "Series: Puiseux inversion represents fractional poles exactly");

    tests.expectEqual(
        eval(session, "series[sqrt[x-1],{x,1,4}]"),
        std::string{"seriesData[x, 1, {1, 0, 0, 0, 0, 0, 0, 0}, 1, 9, 2]"},
        "Series: shifted simple branch points retain the exact expansion center");

    tests.expectEqual(
        eval(session, "series[exp[sqrt[x]],{x,0,3}]"),
        std::string{"seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120, 1/720}, 0, 7, 2]"},
        "Series: analytic TPSA recurrences compose over a Puiseux grid");

    tests.expectEqual(
        eval(session, "normal[series[sqrt[x]*(1+x),{x,0,4}]]"),
        std::string{"x^(3/2)+sqrt[x]"},
        "Series: normal reconstructs fractional powers from SeriesData");

    tests.expectEqual(
        eval(session, "series[(-x)^(1/2),{x,0,4}]"),
        std::string{"series[sqrt[-x], {x, 0, 4}]"},
        "Series: a negative leading direction is not assigned an uncertified principal Puiseux branch");

    tests.expectEqual(
        eval(session, "D[series[exp[x],{x,0,5}],x]"),
        std::string{"seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24}, 0, 5, 1]"},
        "Series: D differentiates SeriesData coefficients directly");

    tests.expectEqual(
        eval(session, "D[series[x^(1/3),{x,0,4}],x]"),
        std::string{"seriesData[x, 0, {1/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, -2, 10, 3]"},
        "Series: D preserves the Puiseux exponent grid");

    tests.expectEqual(
        eval(session, "integrate[series[exp[x],{x,0,5}],x]"),
        std::string{"seriesData[x, 0, {1, 1/2, 1/6, 1/24, 1/120, 1/720}, 1, 7, 1]"},
        "Series: integrate integrates SeriesData coefficients directly");

    tests.expectEqual(
        eval(session, "integrate[series[sqrt[x],{x,0,4}],x]"),
        std::string{"seriesData[x, 0, {2/3, 0, 0, 0, 0, 0, 0, 0}, 3, 11, 2]"},
        "Series: integrate preserves fractional exponents");

    tests.expectEqual(
        eval(session, "integrate[series[1/sqrt[x],{x,0,3}],x]"),
        std::string{"seriesData[x, 0, {2, 0, 0, 0, 0, 0, 0, 0}, 1, 9, 2]"},
        "Series: fractional poles other than x^-1 integrate within Puiseux data");

    tests.expectEqual(
        eval(session, "integrate[series[1/x,{x,0,3}],x]"),
        std::string{"seriesData[x, 0, {0, 0, 0, 0, 0}, 0, 5, 1, {{1, 0, 0, 0, 0}}]"},
        "Series: integrating an x^-1 term closes into logarithmic SeriesData");

    tests.expectEqual(
        eval(session, "series[erf[x],{x,0,7}]"),
        std::string{"seriesData[x, 0, {2/sqrt[Pi], 0, -2/sqrt[Pi]/3, 0, 1/(5sqrt[Pi]), 0, -1/(3sqrt[Pi])/7}, 1, 8, 1]"},
        "Series: erf uses its entire derivative kernel over TPSA");

    tests.expectEqual(
        eval(session, "series[Si[x],{x,0,7}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/18, 0, 1/600, 0, -1/35280}, 1, 8, 1]"},
        "Series: Si fills its removable derivative singularity at zero exactly");

    tests.expectEqual(
        eval(session, "series[Ei[1+x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {Ei[1], E, 0, E/6, -E/12}, 0, 5, 1]"},
        "Series: principal Ei expands at a positive regular center");

    tests.expectEqual(
        eval(session, "series[Ci[1+x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {Ci[1], cos[1 Rad], (-cos[1 Rad]-sin[1 Rad])/2, (cos[1 Rad]/2+sin[1 Rad])/3, (-cos[1 Rad]/2-5sin[1 Rad]/6)/4}, 0, 5, 1]"},
        "Series: principal Ci expands with a radian cosine kernel");

    tests.expectEqual(
        eval(session, "series[li[2+x],{x,0,2}]"),
        std::string{"seriesData[x, 0, {li[2], 1/log[2], -1/(4log[2]^2)}, 0, 3, 1]"},
        "Series: principal li expands at a regular center above one through 1/log[z]");

    tests.expectEqual(
        eval(session, "series[li[a+x],{x,0,2},a>1]"),
        std::string{"seriesData[x, 0, {li[a], 1/log[a], -1/(a log[a]^2)/2}, 0, 3, 1]"},
        "Series: li accepts a symbolic center when assumptions prove the regular principal branch");

    tests.expectEqual(
        eval(session, "series[li[2+sqrt[x]],{x,0,2}]"),
        std::string{"seriesData[x, 0, {li[2], 1/log[2], -1/(4log[2]^2), (1/(4log[2]^3)+1/(8log[2]^2))/3, (-1/(24log[2]^2)-1/(8log[2]^3)-1/(8log[2]^4))/4}, 0, 5, 2]"},
        "Series: li composes over a Puiseux argument away from its cuts and singularity");

    tests.expectEqual(
        eval(session, "series[erf[sqrt[x]],{x,0,4}]"),
        std::string{"seriesData[x, 0, {2/sqrt[Pi], 0, -2/sqrt[Pi]/3, 0, 1/(5sqrt[Pi]), 0, -1/(3sqrt[Pi])/7, 0}, 1, 9, 2]"},
        "Series: special-function primitive providers compose over Puiseux grids");

    tests.expectEqual(
        eval(session, "series[erfc[x],{x,0,7}]"),
        std::string{"seriesData[x, 0, {1, -2/sqrt[Pi], 0, 2/(3sqrt[Pi]), 0, -1/sqrt[Pi]/5, 0, 1/(21sqrt[Pi])}, 0, 8, 1]"},
        "Series: erfc reuses the entire error-function derivative kernel with the complementary sign");

    tests.expectEqual(
        eval(session, "series[fresnelc[x],{x,0,9}]"),
        std::string{"seriesData[x, 0, {1, 0, 0, 0, -Pi^2/40, 0, 0, 0, Pi^4/3456}, 1, 10, 1]"},
        "Series: Fresnel C matches the DLMF Maclaurin coefficients through the shared trigonometric TPSA kernel");

    tests.expectEqual(
        eval(session, "series[fresnels[x],{x,0,11}]"),
        std::string{"seriesData[x, 0, {Pi/6, 0, 0, 0, -Pi Pi^2/336, 0, 0, 0, Pi Pi^4/42240}, 3, 12, 1]"},
        "Series: Fresnel S starts at cubic order and matches the DLMF Maclaurin coefficients");

    tests.expectEqual(
        eval(session, "series[fresnelc[sqrt[x]],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, 0, 0, 0, -Pi^2/40, 0, 0, 0}, 1, 9, 2]"},
        "Series: Fresnel C composes over a half-integer Puiseux grid");

    tests.expectEqual(
        eval(session, "series[fresnels[sqrt[x]],{x,0,4}]"),
        std::string{"seriesData[x, 0, {Pi/6, 0, 0, 0, -Pi Pi^2/336, 0}, 3, 9, 2]"},
        "Series: Fresnel S preserves its cubic zero under Puiseux composition");

    tests.expectEqual(
        eval(session, "series[(1+x)*Si[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 1, -1/18, -1/18, 1/600}, 1, 6, 1]"},
        "Series: special-function valuation composes Si with ordinary products");

    tests.expectEqual(
        eval(session, "series[1/Si[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/18, 0, 23/16200, 0, 209/14288400}, -1, 6, 1]"},
        "Series: special-function valuation exposes the removable Si zero to Laurent inversion");

    tests.expectEqual(
        eval(session, "series[(1+x)*fresnels[x],{x,0,6}]"),
        std::string{"seriesData[x, 0, {Pi/6, Pi/6, 0, 0}, 3, 7, 1]"},
        "Series: special-function valuation preserves the cubic zero of Fresnel S in products");

    tests.expectEqual(
        eval(session, "series[Ei[x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {-digamma[1], 1, 1/4, 1/18, 1/96}, 0, 5, 1, {{1, 0, 0, 0, 0}}]"},
        "Series: Ei at its logarithmic singularity uses logarithmic SeriesData");

    tests.expectEqual(
        eval(session, "series[Ci[-1+x],{x,0,3}]"),
        std::string{"series[Ci[x-1], {x, 0, 3}]"},
        "Series: Ci does not construct a complex neighborhood centered on the principal cut");

    tests.expectEqual(
        eval(session, "series[li[1+x],{x,0,3}]"),
        std::string{"series[li[x+1], {x, 0, 3}]"},
        "Series: li does not cross its singular point at one");

    tests.expectEqual(
        eval(session, "series[li[x],{x,0,3}]"),
        std::string{"series[li[x], {x, 0, 3}]"},
        "Series: li at zero is not forced into pure power or Puiseux data");

    tests.expectEqual(
        eval(session, "series[li[1/2+x],{x,0,3}]"),
        std::string{"series[li[x+1/2], {x, 0, 3}]"},
        "Series: li does not invent a two-sided principal neighborhood on the inherited Ei cut");

    tests.expectEqual(
        eval(session, "series[asin[x],{x,0,7}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/6, 0, 3/40, 0, 5/112}, 1, 8, 1]"},
        "Series: asin matches the DLMF Maclaurin coefficients through its derivative kernel");

    tests.expectEqual(
        eval(session, "series[acos[x],{x,0,6}]"),
        std::string{"seriesData[x, 0, {Pi/2, -1, 0, -1/6, 0, -3/40, 0}, 0, 7, 1]"},
        "Series: acos preserves its principal constant and the negative asin derivative kernel");

    tests.expectEqual(
        eval(session, "series[atan[x],{x,0,7}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/3, 0, 1/5, 0, -1/7}, 1, 8, 1]"},
        "Series: atan matches the alternating DLMF Maclaurin series");

    tests.expectEqual(
        eval(session, "series[asin[a+x],{x,0,2},-1<a<1]"),
        std::string{"seriesData[x, 0, {asin[a], (1-a^2)^(-1/2), a*(1-a^2)^(-1/2)/(2(1-a^2))}, 0, 3, 1]"},
        "Series: assumptions certify a symbolic regular center for principal asin");

    tests.expectEqual(
        eval(session, "series[atan[a+x],{x,0,2},a>0]"),
        std::string{"seriesData[x, 0, {atan[a], 1/(a^2+1), -2a/(a^2+1)^2/2}, 0, 3, 1]"},
        "Series: real assumptions certify regular atan centers");

    tests.expectEqual(
        eval(session, "series[asin[I+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {asin[I], 2^(-1/2), I*2^(-1/2)/4, -2^(-1/2)/24}, 0, 4, 1]"},
        "Series: asin expands at a nonreal center away from its real-axis cuts");

    tests.expectEqual(
        eval(session, "series[atan[1+I+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {atan[1+I], 1/5-2I/5, -1/25+7I/25, -1/375-68I/375}, 0, 4, 1]"},
        "Series: atan accepts complex centers with a provably nonzero real part");

    tests.expectEqual(
        eval(session, "series[atan[I/2+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {atan[I/2], 4/3, -8I/9, -112/81}, 0, 4, 1]"},
        "Series: atan accepts the imaginary-axis segment between its branch points");

    tests.expectEqual(
        eval(session, "series[asin[sqrt[x]],{x,0,3}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/6, 0, 3/40, 0}, 1, 7, 2]"},
        "Series: inverse trigonometric providers compose over an existing Puiseux grid");

    tests.expectEqual(
        eval(session, "series[1/asin[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/6, 0, -17/360, 0, -367/15120}, -1, 6, 1]"},
        "Series: inverse-trigonometric valuation exposes the simple asin zero to Laurent inversion");

    tests.expectEqual(
        eval(session, "series[asin[1+x],{x,0,3}]"),
        std::string{"series[asin[x+1], {x, 0, 3}]"},
        "Series: asin does not force a Taylor series at its branch point");
    tests.expect(!session.diagnostics().empty()
            && session.diagnostics().back().code == "series::unsupported",
        "Series: inverse-trigonometric branch-point rejection emits an explicit diagnostic");

    tests.expectEqual(
        eval(session, "series[atan[I+x],{x,0,3}]"),
        std::string{"series[atan[x+I], {x, 0, 3}]"},
        "Series: atan does not cross its principal branch point at I");

    tests.expectEqual(
        eval(session, "series[lambertw[x],{x,0,7}]"),
        std::string{"seriesData[x, 0, {1, -1, 3/2, -8/3, 125/24, -54/5, 16807/720}, 1, 8, 1]"},
        "Series: principal Lambert W matches the DLMF Maclaurin coefficients");

    tests.expectEqual(
        eval(session, "series[lambertw[E+x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, exp[-1]/2, -3*exp[-2]/16, 19*exp[-3]/192, -185*exp[-4]/3072}, 0, 5, 1]"},
        "Series: Lambert W uses exact special values at a regular center");

    tests.expectEqual(
        eval(session, "series[lambertw[a+x],{x,0,2},a>-1/E]"),
        std::string{"seriesData[x, 0, {lambertw[a], exp[-lambertw[a]]/(1+lambertw[a]), (-2-lambertw[a])exp[-2lambertw[a]]/(2(1+lambertw[a])^3)}, 0, 3, 1]"},
        "Series: assumptions certify symbolic regular centers for the principal Lambert W branch");

    tests.expectEqual(
        eval(session, "series[lambertw[-1,a+x],{x,0,2},a>0]"),
        std::string{"seriesData[x, 0, {lambertw[-1, a], exp[-lambertw[-1, a]]/(1+lambertw[-1, a]), (-2-lambertw[-1, a])exp[-2lambertw[-1, a]]/(2(1+lambertw[-1, a])^3)}, 0, 3, 1]"},
        "Series: explicit nonprincipal Lambert W branches expand away from their cut");

    tests.expectEqual(
        eval(session, "series[lambertw[sqrt[x]],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, -1, 3/2, -8/3, 125/24, -54/5, 16807/720, -16384/315}, 1, 9, 2]"},
        "Series: Lambert W composes over an existing Puiseux grid");

    tests.expectEqual(
        eval(session, "series[1/lambertw[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 1, -1/2, 2/3, -9/8, 32/15, -625/144}, -1, 6, 1]"},
        "Series: Lambert W valuation exposes the simple principal zero to Laurent inversion");

    tests.expectEqual(
        eval(session, "series[lambertw[-1/E+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {-1, sqrt[2]sqrt[E], -2*E/3, 11*E sqrt[2]sqrt[E]/36, -43*exp[2]/135, 769*exp[2]sqrt[2]sqrt[E]/4320, -1768*E exp[2]/8505}, 0, 7, 2]"},
        "Series: principal Lambert W uses the DLMF square-root Puiseux expansion at -1/E");

    tests.expectEqual(
        eval(session, "series[lambertw[-1,-1/E+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {-1, -sqrt[2]sqrt[E], -2*E/3, -11*E sqrt[2]sqrt[E]/36, -43*exp[2]/135, -769*exp[2]sqrt[2]sqrt[E]/4320, -1768*E exp[2]/8505}, 0, 7, 2]"},
        "Series: the -1 Lambert W branch selects the opposite square-root sheet at the branch point");

    tests.expectEqual(
        eval(session, "series[1/lambertw[-1/E+x],{x,0,2}]"),
        std::string{"seriesData[x, 0, {-1, -sqrt[2]sqrt[E], -4*E/3, -35*E sqrt[2]sqrt[E]/36, -182*exp[2]/135}, 0, 5, 2]"},
        "Series: branch-point Lambert W valuation remains available to Laurent-style inversion");

    tests.expectEqual(
        eval(session, "series[lambertw[-1/E+sqrt[x]],{x,0,2}]"),
        std::string{"seriesData[x, 0, {-1, sqrt[2]sqrt[E], -2*E/3, 11*E sqrt[2]sqrt[E]/36, -43*exp[2]/135, 769*exp[2]sqrt[2]sqrt[E]/4320, -1768*E exp[2]/8505, 680863*E exp[2]sqrt[2]sqrt[E]/5443200, -3926*exp[2]^2/25515}, 0, 9, 4]"},
        "Series: Lambert W branch-point expansion composes over an existing Puiseux grid");

    tests.expectEqual(
        eval(session, "series[lambertw[-1/E-x],{x,0,2}]"),
        std::string{"series[lambertw[-1/E-x], {x, 0, 2}]"},
        "Series: Lambert W branch-point expansion does not choose a principal square root across the cut");
    tests.expect(!session.diagnostics().empty()
            && session.diagnostics().back().code == "series::unsupported",
        "Series: rejected Lambert W branch direction emits an explicit diagnostic");

    tests.expectEqual(
        eval(session, "series[lambertw[-1/E+x^2],{x,0,2}]"),
        std::string{"series[lambertw[-1/E+x^2], {x, 0, 2}]"},
        "Series: Lambert W does not collapse the ambiguous sqrt[x^2] branch at a multiple contact");

    tests.expectEqual(
        eval(session, "series[lambertw[1,-1/E+x],{x,0,2}]"),
        std::string{"series[lambertw[1, x-1/E], {x, 0, 2}]"},
        "Series: only the principal and -1 sheets use the implemented -1/E Puiseux provider");

    tests.expectEqual(
        eval(session, "series[lambertw[-1,-1/10+x],{x,0,3}]"),
        std::string{"series[lambertw[-1, x-1/10], {x, 0, 3}]"},
        "Series: nonprincipal Lambert W does not invent a two-sided neighborhood on its cut");

    tests.expectEqual(
        eval(session, "series[lambertw[k,1+x],{x,0,2}]"),
        std::string{"series[lambertw[k, x+1], {x, 0, 2}]"},
        "Series: symbolic Lambert W branch indices remain unevaluated");

    tests.expectEqual(
        eval(session, "series[gamma[1+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {1, digamma[1], (Pi^2/6+digamma[1]^2)/2, ((Pi^2/6+digamma[1]^2)digamma[1]/2+Pi^2digamma[1]/6-zeta[3])/3}, 0, 4, 1]"},
        "Series: gamma at one follows the DLMF log-Gamma coefficients and TPSA exponentiation");

    tests.expectEqual(
        eval(session, "series[lgamma[1+x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {digamma[1], Pi^2/12, -zeta[3]/3, Pi^4/360, -zeta[5]/5}, 1, 6, 1]"},
        "Series: lgamma at one matches the DLMF Taylor series without a cancellation-heavy Gamma round trip");

    tests.expectEqual(
        eval(session, "series[gamma[1/2+x],{x,0,2}]"),
        std::string{"seriesData[x, 0, {sqrt[Pi], (digamma[1]-2log[2])sqrt[Pi], (Pi^2/2+(digamma[1]-2log[2])^2)sqrt[Pi]/2}, 0, 3, 1]"},
        "Series: gamma uses the exact half-integer base expansion");

    tests.expectEqual(
        eval(session, "series[gamma[-1/2+x],{x,0,2}]"),
        std::string{"seriesData[x, 0, {-2sqrt[Pi], -2(2+digamma[1]-2log[2])sqrt[Pi], -(2(2+Pi^2/4)+(2+digamma[1]-2log[2])^2)sqrt[Pi]}, 0, 3, 1]"},
        "Series: gamma shifts the half-integer base across a regular negative center");

    tests.expectEqual(
        eval(session, "series[lgamma[-1/2+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {log[2sqrt[Pi]], 2+digamma[1]-2log[2], 2+Pi^2/4, 8/3-7zeta[3]/3}, 0, 4, 1]"},
        "Series: real lgamma preserves log-absolute-value semantics at a negative half-integer center");

    tests.expectEqual(
        eval(session, "series[gamma[1+sqrt[x]],{x,0,2}]"),
        std::string{"seriesData[x, 0, {1, digamma[1], (Pi^2/6+digamma[1]^2)/2, ((Pi^2/6+digamma[1]^2)digamma[1]/2+Pi^2digamma[1]/6-zeta[3])/3, (((Pi^2/6+digamma[1]^2)digamma[1]/2+Pi^2digamma[1]/6-zeta[3])digamma[1]/3+(Pi^2/6+digamma[1]^2)Pi^2/12-digamma[1]zeta[3]+Pi^4/90)/4}, 0, 5, 2]"},
        "Series: gamma composes over an existing Puiseux grid");

    tests.expectEqual(
        eval(session, "series[gamma[1+I*x],{x,0,2}]"),
        std::string{"seriesData[x, 0, {1, I digamma[1], (-Pi^2/6-digamma[1]^2)/2}, 0, 3, 1]"},
        "Series: meromorphic gamma permits complex local directions away from its poles");

    tests.expectEqual(
        eval(session, "series[(1+x)*lgamma[1+x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {digamma[1], Pi^2/12+digamma[1], Pi^2/12-zeta[3]/3, Pi^4/360-zeta[3]/3}, 1, 5, 1]"},
        "Series: lgamma valuation preserves its simple zero at one in products");

    tests.expectEqual(
        eval(session, "series[gamma[x],{x,0,3}]"),
        std::string{"series[gamma[x], {x, 0, 3}]"},
        "Series: gamma poles are not forced into the regular-center provider");
    tests.expect(!session.diagnostics().empty()
            && session.diagnostics().back().code == "series::unsupported",
        "Series: unsupported gamma pole expansions emit an explicit diagnostic");

    tests.expectEqual(
        eval(session, "series[gamma[1/3+x],{x,0,3}]"),
        std::string{"series[gamma[x+1/3], {x, 0, 3}]"},
        "Series: gamma centers without an exact local coefficient basis remain unevaluated");

    tests.expectEqual(
        eval(session, "series[lgamma[1+I*x],{x,0,3}]"),
        std::string{"series[lgamma[I x+1], {x, 0, 3}]"},
        "Series: lgamma does not invent a complex analytic continuation beyond its real log-absolute-value contract");

    tests.expectEqual(
        eval(session, "series[digamma[1+x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {digamma[1], Pi^2/6, -zeta[3], Pi^4/90, -zeta[5], zeta[6]}, 0, 6, 1]"},
        "Series: digamma at one follows the derivative of the DLMF log-Gamma expansion");

    tests.expectEqual(
        eval(session, "series[trigamma[1+x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {Pi^2/6, -2zeta[3], Pi^4/30, -4zeta[5], 5zeta[6]}, 0, 5, 1]"},
        "Series: trigamma at one differentiates the digamma coefficient family exactly");

    tests.expectEqual(
        eval(session, "series[digamma[2+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {1+digamma[1], -1+Pi^2/6, 1-zeta[3], -1+Pi^4/90}, 0, 4, 1]"},
        "Series: digamma recurrence shifts the base expansion to a positive integer center");

    tests.expectEqual(
        eval(session, "series[trigamma[2+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {Pi^2/6-1, 2-2zeta[3], Pi^4/30-3, 4-4zeta[5]}, 0, 4, 1]"},
        "Series: trigamma recurrence shifts the base expansion to a positive integer center");

    tests.expectEqual(
        eval(session, "series[digamma[1/2+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {digamma[1]-2log[2], Pi^2/2, -7zeta[3], Pi^4/6}, 0, 4, 1]"},
        "Series: digamma uses the exact half-integer local coefficient basis");

    tests.expectEqual(
        eval(session, "series[trigamma[-1/2+x],{x,0,2}]"),
        std::string{"seriesData[x, 0, {4+Pi^2/2, 16-14zeta[3], 48+Pi^4/2}, 0, 3, 1]"},
        "Series: trigamma recurrence crosses to a regular negative half-integer center");

    tests.expectEqual(
        eval(session, "series[digamma[1+sqrt[x]],{x,0,2}]"),
        std::string{"seriesData[x, 0, {digamma[1], Pi^2/6, -zeta[3], Pi^4/90, -zeta[5]}, 0, 5, 2]"},
        "Series: digamma composes over an existing Puiseux grid");

    tests.expectEqual(
        eval(session, "series[trigamma[1+sqrt[x]],{x,0,2}]"),
        std::string{"seriesData[x, 0, {Pi^2/6, -2zeta[3], Pi^4/30, -4zeta[5], 5zeta[6]}, 0, 5, 2]"},
        "Series: trigamma composes over an existing Puiseux grid");

    tests.expectEqual(
        eval(session, "series[trigamma[1+I*x],{x,0,2}]"),
        std::string{"seriesData[x, 0, {Pi^2/6, -2I zeta[3], -Pi^4/30}, 0, 3, 1]"},
        "Series: trigamma permits complex local directions away from its poles");

    tests.expectEqual(
        eval(session, "D[series[digamma[1+x],{x,0,5}],x]"),
        std::string{"seriesData[x, 0, {Pi^2/6, -2zeta[3], Pi^4/30, -4zeta[5], 5zeta[6]}, 0, 5, 1]"},
        "Series: differentiating a digamma SeriesData agrees with the trigamma provider");

    tests.expectEqual(
        eval(session, "series[(1+x)*digamma[1+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {digamma[1], Pi^2/6+digamma[1], Pi^2/6-zeta[3], Pi^4/90-zeta[3]}, 0, 4, 1]"},
        "Series: digamma valuation remains finite in products at supported regular centers");

    tests.expectEqual(
        eval(session, "series[digamma[x],{x,0,3}]"),
        std::string{"series[digamma[x], {x, 0, 3}]"},
        "Series: digamma poles are not forced into the regular-center provider");
    tests.expect(!session.diagnostics().empty()
            && session.diagnostics().back().code == "series::unsupported",
        "Series: unsupported digamma pole expansions emit an explicit diagnostic");

    tests.expectEqual(
        eval(session, "series[trigamma[x],{x,0,3}]"),
        std::string{"series[trigamma[x], {x, 0, 3}]"},
        "Series: trigamma poles are not forced into the regular-center provider");

    tests.expectEqual(
        eval(session, "series[digamma[1/3+x],{x,0,3}]"),
        std::string{"series[digamma[x+1/3], {x, 0, 3}]"},
        "Series: digamma centers without an exact local coefficient basis remain unevaluated");

    tests.expectEqual(
        eval(session, "series[polylog[2,x],{x,0,6}]"),
        std::string{"seriesData[x, 0, {1, 1/4, 1/9, 1/16, 1/25, 1/36}, 1, 7, 1]"},
        "Series: polylog at the origin follows the DLMF defining power series");

    tests.expectEqual(
        eval(session, "series[polylog[a,x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, 2^(-a), 3^(-a), 4^(-a)}, 1, 5, 1]"},
        "Series: polylog keeps a variable-independent symbolic order exact at the origin");

    tests.expectEqual(
        eval(session, "series[polylog[1/2,x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, 2^(-1/2), 3^(-1/2), 4^(-1/2)}, 1, 5, 1]"},
        "Series: polylog supports exact noninteger order coefficients at the origin");

    tests.expectEqual(
        eval(session, "series[polylog[2,x+x^2],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 5/4, 11/18, 31/48, 187/300}, 1, 6, 1]"},
        "Series: polylog origin coefficients compose through an ordinary TPSA argument");

    tests.expectEqual(
        eval(session, "series[polylog[2,sqrt[x]],{x,0,3}]"),
        std::string{"seriesData[x, 0, {1, 1/4, 1/9, 1/16, 1/25, 1/36}, 1, 7, 2]"},
        "Series: polylog origin coefficients compose over a square-root Puiseux grid");

    tests.expectEqual(
        eval(session, "series[polylog[2,x^(1/3)],{x,0,2}]"),
        std::string{"seriesData[x, 0, {1, 1/4, 1/9, 1/16, 1/25, 1/36}, 1, 7, 3]"},
        "Series: polylog origin coefficients preserve a cubic Puiseux grid");

    tests.expectEqual(
        eval(session, "series[polylog[2,-sqrt[x]],{x,0,2}]"),
        std::string{"seriesData[x, 0, {-1, 1/4, -1/9, 1/16}, 1, 5, 2]"},
        "Series: polylog origin composition preserves the sign of the local Puiseux parameter");

    tests.expectEqual(
        eval(session, "series[1/polylog[2,x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, -1/4, -7/144, -13/576, -6911/518400, -6151/691200}, -1, 5, 1]"},
        "Series: polylog valuation exposes its simple zero at the origin to Laurent inversion");

    tests.expectEqual(
        eval(session, "D[series[polylog[3,x],{x,0,5}],x]"),
        std::string{"seriesData[x, 0, {1, 1/4, 1/9, 1/16, 1/25}, 0, 5, 1]"},
        "Series: differentiating the Li_3 origin SeriesData matches Li_2 divided by x");

    tests.expectEqual(
        eval(session, "series[polylog[2,x]/x,{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, 1/4, 1/9, 1/16, 1/25}, 0, 5, 1]"},
        "Series: polylog derivative kernel composes through Laurent cancellation at the origin");

    tests.expectEqual(
        eval(session, "series[(1+x)*polylog[2,x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1, 5/4, 13/36, 25/144}, 1, 5, 1]"},
        "Series: polylog origin valuation composes in products");

    tests.expectEqual(
        eval(session, "series[polylog[1,1/2+x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {-log[1/2], 2, 2, 8/3, 4}, 0, 5, 1]"},
        "Series: Li_1 at a regular nonzero center reduces to the exact principal logarithm coefficients");

    tests.expectEqual(
        eval(session, "series[polylog[2,1/2+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {polylog[2, 1/2], -2log[1/2], 2(1+log[1/2]), 4(-1-2log[1/2])/3}, 0, 4, 1]"},
        "Series: positive-integer polylog expands at an exact regular center below the principal cut");

    tests.expectEqual(
        eval(session, "series[polylog[3,-1+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {-3zeta[3]/4, Pi^2/12, (Pi^2/12-log[2])/2, -(-1/2-Pi^2/6+3log[2])/6}, 0, 4, 1]"},
        "Series: polylog regular-center recursion reaches negative real centers off the principal cut");

    tests.expectEqual(
        eval(session, "series[polylog[2,I+x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {polylog[2, I], I log[1-I], -(-1/2+I/2+log[1-I])/2, I*(1-3I/2-2log[1-I])/6}, 0, 4, 1]"},
        "Series: polylog regular-center recursion supports nonreal centers");

    tests.expectEqual(
        eval(session, "series[polylog[2,a+x],{x,0,2},{a>0,a<1}]"),
        std::string{"seriesData[x, 0, {polylog[2, a], -log[1-a]/a, (a/(1-a)+log[1-a])/(2a^2)}, 0, 3, 1]"},
        "Series: assumptions can certify a symbolic polylog center inside the real regular interval");

    tests.expectEqual(
        eval(session, "series[polylog[2,1/2+sqrt[x]],{x,0,2}]"),
        std::string{"seriesData[x, 0, {polylog[2, 1/2], -2log[1/2], 2(1+log[1/2]), 4(-1-2log[1/2])/3, 2(5+6log[1/2])/3}, 0, 5, 2]"},
        "Series: regular-center polylog composes over the existing Puiseux grid");

    tests.expectEqual(
        eval(session, "series[1/polylog[2,1/2+x],{x,0,2}]"),
        std::string{"seriesData[x, 0, {1/polylog[2, 1/2], 2log[1/2]/polylog[2, 1/2]^2, -(2(1+log[1/2])/polylog[2, 1/2]-4log[1/2]^2/polylog[2, 1/2]^2)/polylog[2, 1/2]}, 0, 3, 1]"},
        "Series: positive polylog values on zero-to-one propagate nonzero evidence into reciprocal inversion");

    tests.expectEqual(
        eval(session, "series[polylog[2,1/2+x]^-1,{x,0,2}]"),
        std::string{"seriesData[x, 0, {1/polylog[2, 1/2], 2log[1/2]/polylog[2, 1/2]^2, -(2(1+log[1/2])/polylog[2, 1/2]-4log[1/2]^2/polylog[2, 1/2]^2)/polylog[2, 1/2]}, 0, 3, 1]"},
        "Series: negative integer powers receive the same certified polylog leading-coefficient evidence");

    tests.expectEqual(
        eval(session, "series[polylog[2,2+x],{x,0,2}]"),
        std::string{"series[polylog[2, x+2], {x, 0, 2}]"},
        "Series: polylog centers on the principal cut remain unevaluated");
    tests.expect(!session.diagnostics().empty()
            && session.diagnostics().back().code == "series::unsupported",
        "Series: unsupported polylog cut centers emit an explicit diagnostic");

    tests.expectEqual(
        eval(session, "series[polylog[1/2,1/2+x],{x,0,2}]"),
        std::string{"series[polylog[1/2, x+1/2], {x, 0, 2}]"},
        "Series: noninteger-order polylog regular centers remain outside the exact provider");

    tests.expectEqual(
        eval(session, "D[series[polylog[3,1/2+x],{x,0,3}],x]"),
        std::string{"seriesData[x, 0, {2polylog[2, 1/2], 4(-log[1/2]-polylog[2, 1/2]), 4(1+3log[1/2]+2polylog[2, 1/2])}, 0, 3, 1]"},
        "Series: differentiated regular-center Li_3 agrees with the Li_2 over z recurrence");

    tests.expectEqual(
        eval(session, "series[polylog[x,x],{x,0,3}]"),
        std::string{"series[polylog[x, x], {x, 0, 3}]"},
        "Series: polylog order parameters depending on the expansion variable remain unevaluated");

    tests.expectEqual(
        eval(session, "series[log[x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {0, 0, 0, 0, 0}, 0, 5, 1, {{1, 0, 0, 0, 0}}]"},
        "Series: logarithmic SeriesData represents the principal log singularity at the origin");

    tests.expectEqual(
        eval(session, "normal[series[log[x],{x,0,4}]]"),
        std::string{"log[x]"},
        "Series: normal reconstructs logarithmic coefficient layers");

    tests.expectEqual(
        eval(session, "D[series[log[x],{x,0,4}],x]"),
        std::string{"seriesData[x, 0, {1, 0, 0, 0, 0}, -1, 4, 1]"},
        "Series: differentiating a logarithmic layer produces the expected Laurent pole");

    tests.expectEqual(
        eval(session, "integrate[series[1/x,{x,0,4}],x]"),
        std::string{"seriesData[x, 0, {0, 0, 0, 0, 0, 0}, 0, 6, 1, {{1, 0, 0, 0, 0, 0}}]"},
        "Series: integrating a simple Laurent pole closes into a logarithmic layer");

    tests.expectEqual(
        eval(session, "normal[seriesData[x,0,{0},0,1,1,{{0},{1}}]]"),
        std::string{"log[x]^2"},
        "Series: logarithmic SeriesData supports higher powers of the local logarithm");

    tests.expectEqual(
        eval(session, "D[seriesData[x,0,{0},0,1,1,{{0},{1}}],x]"),
        std::string{"seriesData[x, 0, {0}, -1, 0, 1, {{2}}]"},
        "Series: differentiating log squared lowers the logarithmic degree exactly");

    tests.expectEqual(
        eval(session, "integrate[seriesData[x,0,{0},-1,0,1,{{2}}],x]"),
        std::string{"seriesData[x, 0, {0}, 0, 1, 1, {{0}, {1}}]"},
        "Series: integrating x^-1 times log raises the logarithmic degree exactly");

    tests.expectEqual(
        eval(session, "series[Ei[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {-digamma[1], 1, 1/4, 1/18, 1/96, 1/600}, 0, 6, 1, {{1, 0, 0, 0, 0, 0}}]"},
        "Series: Ei at the origin follows the DLMF logarithmic power series");

    tests.expectEqual(
        eval(session, "series[Ci[x],{x,0,6}]"),
        std::string{"seriesData[x, 0, {-digamma[1], 0, -1/4, 0, 1/96, 0, -1/4320}, 0, 7, 1, {{1, 0, 0, 0, 0, 0, 0}}]"},
        "Series: Ci at the origin follows the DLMF logarithmic even-power series");

    tests.expectEqual(
        eval(session, "D[series[Ei[x],{x,0,5}],x]"),
        std::string{"seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120}, -1, 5, 1]"},
        "Series: differentiated Ei logarithmic data matches exp[x]/x locally");

    tests.expectEqual(
        eval(session, "series[x*log[x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {0, 0, 0, 0}, 1, 5, 1, {{1, 0, 0, 0}}]"},
        "Series: logarithmic layers multiply by ordinary Taylor monomials");

    tests.expectEqual(
        eval(session, "series[log[x]^2,{x,0,4}]"),
        std::string{"seriesData[x, 0, {0, 0, 0, 0, 0}, 0, 5, 1, {{0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}}]"},
        "Series: nonnegative integer powers of a logarithmic layer remain closed");

    tests.expectEqual(
        eval(session, "series[Ei[x]+Ci[x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {-2digamma[1], 1, 0, 1/18, 1/48}, 0, 5, 1, {{2, 0, 0, 0, 0}}]"},
        "Series: logarithmic special-function layers add coefficient-wise");

    tests.expectEqual(
        eval(session, "series[sqrt[x]*log[x],{x,0,3}]"),
        std::string{"seriesData[x, 0, {0, 0, 0, 0, 0, 0}, 1, 7, 2, {{1, 0, 0, 0, 0, 0}}]"},
        "Series: logarithmic layers align with an existing Puiseux exponent grid");

    tests.expectEqual(
        eval(session, "D[series[x*log[x],{x,0,4}],x]"),
        std::string{"seriesData[x, 0, {1, 0, 0, 0}, 0, 4, 1, {{1, 0, 0, 0}}]"},
        "Series: differentiated x log x retains both power and logarithmic contributions");

    tests.expectEqual(
        eval(session, "integrate[series[log[x]^2,{x,0,3}],x]"),
        std::string{"seriesData[x, 0, {2, 0, 0, 0}, 1, 5, 1, {{-2, 0, 0, 0}, {1, 0, 0, 0}}]"},
        "Series: integration lowers logarithmic degree terms by exact recurrence");

    tests.expectEqual(
        eval(session, "series[log[2*x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {log[2], 0, 0, 0, 0}, 0, 5, 1, {{1, 0, 0, 0, 0}}]"},
        "Series: principal logarithmic composition extracts a positive leading constant");

    tests.expectEqual(
        eval(session, "series[Ei[2*x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {-digamma[1]+log[2], 2, 1, 4/9, 1/6}, 0, 5, 1, {{1, 0, 0, 0, 0}}]"},
        "Series: Ei logarithmic composition follows a positive simple zero");

    tests.expectEqual(
        eval(session, "series[Ci[2*x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {-digamma[1]+log[2], 0, -1, 0, 1/6}, 0, 5, 1, {{1, 0, 0, 0, 0}}]"},
        "Series: Ci logarithmic composition follows a positive simple zero");

    tests.expectEqual(
        eval(session, "series[log[sqrt[x]],{x,0,4}]"),
        std::string{"seriesData[x, 0, {0, 0, 0, 0, 0, 0, 0, 0, 0}, 0, 9, 2, {{1/2, 0, 0, 0, 0, 0, 0, 0, 0}}]"},
        "Series: principal logarithmic composition preserves a certified Puiseux zero");

    tests.expectEqual(
        eval(session, "series[Ei[sqrt[x]],{x,0,3}]"),
        std::string{"seriesData[x, 0, {-digamma[1], 1, 1/4, 1/18, 1/96, 1/600, 1/4320}, 0, 7, 2, {{1/2, 0, 0, 0, 0, 0, 0}}]"},
        "Series: Ei logarithmic composition preserves a certified Puiseux zero");

    tests.expectEqual(
        eval(session, "series[log[x+x^2],{x,0,4}]"),
        std::string{"seriesData[x, 0, {0, 1, -1/2, 1/3, -1/4}, 0, 5, 1, {{1, 0, 0, 0, 0}}]"},
        "Series: logarithmic composition separates a simple zero from its regular unit");

    tests.expectEqual(
        eval(session, "series[log[x^2],{x,0,4}]"),
        std::string{"series[log[x^2], {x, 0, 4}]"},
        "Series: higher winding logarithmic zeros remain unevaluated rather than changing the principal branch");

    tests.expectEqual(
        eval(session, "series[log[-x],{x,0,4}]"),
        std::string{"series[log[-x], {x, 0, 4}]"},
        "Series: negative leading logarithmic directions remain unevaluated");

    kernel::KernelSession specialDegreeSession;
    static_cast<void>(specialDegreeSession.evaluate("angleMode[Deg]"));
    tests.expectEqual(
        eval(specialDegreeSession, "series[asin[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {180/Pi, 0, 30/Pi, 0, 27/(2Pi)}, 1, 6, 1]"},
        "Series: inverse trigonometric coefficients use the session output angle unit");
    tests.expectEqual(
        eval(specialDegreeSession, "series[Si[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/18, 0, 1/600}, 1, 6, 1]"},
        "Series: Si uses its intrinsic radian kernel independently of session angle mode");

    tests.expectEqual(
        eval(specialDegreeSession, "series[fresnelc[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 0, 0, -Pi^2/40}, 1, 6, 1]"},
        "Series: Fresnel definitions use intrinsic radians independently of session angle mode");


    tests.expectEqual(
        eval(specialDegreeSession, "series[fresnels[x],{x,0,7}]"),
        std::string{"seriesData[x, 0, {Pi/6, 0, 0, 0, -Pi Pi^2/336}, 3, 8, 1]"},
        "Series: Fresnel S also keeps its intrinsic radian kernel in degree mode");

    tests.expectEqual(
        eval(session, "toNormal[series[(1+x)^2,{x,0,3}]]"),
        std::string{"x^2+2x+1"},
        "Series: toNormal converts a top-level SeriesData object");

    tests.expectEqual(
        eval(session, "toNormal[{series[(1+x)^2,{x,0,3}],series[log[x],{x,0,2}]}]"),
        std::string{"{x^2+2x+1, log[x]}"},
        "Series: toNormal recursively converts SeriesData inside lists");

    tests.expectEqual(
        eval(session, "toNormal[{{series[(1+x)^2,{x,0,3}]},{series[log[x],{x,0,2}]}}]"),
        std::string{"{{x^2+2x+1}, {log[x]}}"},
        "Series: toNormal recursively converts nested containers");

    const std::string toNormalOnce =
        eval(session, "toNormal[series[exp[x],{x,0,3}]]");
    tests.expectEqual(
        eval(session, "toNormal[toNormal[series[exp[x],{x,0,3}]]]"),
        toNormalOnce,
        "Series: toNormal is idempotent for converted SeriesData");

    tests.expectEqual(
        eval(session, "toNormal[42]"), std::string{"42"},
        "Series: toNormal leaves ordinary expressions unchanged");

    tests.expectEqual(
        eval(session, "normal[{series[(1+x)^2,{x,0,3}]}]"),
        std::string{"{seriesData[x, 0, {1, 2, 1, 0}, 0, 4, 1]}"},
        "Series: normal retains its top-level-only compatibility semantics");

    tests.expectEqual(
        eval(session, "toNormal[cases[seriesData[x,0,{1,2,1},0,3,1] if a>0;1 if a<=0]]"),
        std::string{"cases[x^2+2x+1 if a > 0; 1 if a <= 0]"},
        "Series: toNormal recursively converts SeriesData retained inside calls");

    tests.expectEqual(
        eval(session, "series[x,{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {1, 0, 0, 0, 0}, -1, 4, 1]"},
        "Series: Infinity center uses the reciprocal local variable");

    tests.expectEqual(
        eval(session, "series[1/(x+1),{x,Infinity,4}]"),
        std::string{"seriesData[x, Infinity, {1, -1, 1, -1}, 1, 5, 1]"},
        "Series: rational functions expand at positive infinity");

    tests.expectEqual(
        eval(session, "series[exp[1/x],{x,Infinity,4}]"),
        std::string{"seriesData[x, Infinity, {1, 1, 1/2, 1/6, 1/24}, 0, 5, 1]"},
        "Series: analytic functions of the reciprocal reuse TPSA at infinity");

    tests.expectEqual(
        eval(session, "series[sqrt[x],{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {1, 0, 0, 0, 0, 0, 0, 0}, -1, 7, 2]"},
        "Series: Puiseux exponents are retained at positive infinity");

    tests.expectEqual(
        eval(session, "series[log[1/x],{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {0, 0, 0, 0}, 0, 4, 1, {{1, 0, 0, 0}}]"},
        "Series: logarithmic layers use log of the reciprocal at infinity");

    tests.expectEqual(
        eval(session, "normal[series[x^2+1/x,{x,Infinity,3}]]"),
        std::string{"x^(-1)+x^2"},
        "Series: normal emits ordinary powers of x for Infinity SeriesData");

    tests.expectEqual(
        eval(session, "toNormal[{series[1/(x+1),{x,Infinity,3}]}]"),
        std::string{"{x^(-1)-x^(-2)+x^(-3)}"},
        "Series: toNormal recursively normalizes Infinity SeriesData");

    tests.expectEqual(
        eval(session, "D[series[1/(x+1),{x,Infinity,4}],x]"),
        std::string{"seriesData[x, Infinity, {-1, 2, -3, 4}, 2, 6, 1]"},
        "Series: differentiation uses the reciprocal-variable chain rule at infinity");

    tests.expectEqual(
        eval(session, "integrate[series[1/x^2,{x,Infinity,3}],x]"),
        std::string{"seriesData[x, Infinity, {-1, 0}, 1, 3, 1]"},
        "Series: integration shifts reciprocal exponents at infinity");

    tests.expectEqual(
        eval(session, "integrate[series[1/x,{x,Infinity,3}],x]"),
        std::string{"seriesData[x, Infinity, {0, 0, 0}, 0, 3, 1, {{-1, 0, 0}}]"},
        "Series: integration of 1/x closes into the logarithmic layer at infinity");

    tests.expectEqual(
        eval(session, "integrate[series[1/x,{x,Infinity,0}],x]"),
        std::string{"integrate[seriesData[x, Infinity, {0}, 0, 1, 1], x]"},
        "Series: integration rejects an O(1/x) remainder that can generate an unknown logarithm");

    tests.expectEqual(
        eval(session, "series[sin[x],{x,Infinity,3}]"),
        std::string{"series[sin[x], {x, Infinity, 3}]"},
        "Series: oscillatory infinity expansions remain unevaluated");

    tests.expectEqual(
        eval(session, "series[exp[x],{x,Infinity,3}]"),
        std::string{"series[exp[x], {x, Infinity, 3}]"},
        "Series: essential exponential growth at infinity remains unevaluated");

    tests.expectEqual(
        eval(session, "series[log[x],{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {0, 0, 0, 0}, 0, 4, 1, {{-1, 0, 0, 0}}]"},
        "Series: direct logarithmic growth at positive infinity uses the logarithmic layer");

    tests.expectEqual(
        eval(session, "normal[series[log[x],{x,Infinity,3}]]"),
        std::string{"log[x]"},
        "Series: normal emits log[x] for positive-infinity logarithmic layers");

    tests.expectEqual(
        eval(session, "series[log[2*x],{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {log[2], 0, 0, 0}, 0, 4, 1, {{-1, 0, 0, 0}}]"},
        "Series: positive leading constants are separated from logarithmic growth at infinity");

    tests.expectEqual(
        eval(session, "series[log[a*x],{x,Infinity,3},a>0]"),
        std::string{"seriesData[x, Infinity, {log[a], 0, 0, 0}, 0, 4, 1, {{-1, 0, 0, 0}}]"},
        "Series: positive assumptions certify symbolic leading constants at infinity");

    tests.expectEqual(
        eval(session, "series[log[a*x],{x,Infinity,3}]"),
        std::string{"series[log[a x], {x, Infinity, 3}]"},
        "Series: symbolic logarithmic leading constants require a positive proof at infinity");

    tests.expectEqual(
        eval(session, "series[log[(x+1)/(x+2)],{x,Infinity,4}]"),
        std::string{"seriesData[x, Infinity, {0, -1, 3/2, -7/3, 15/4}, 0, 5, 1]"},
        "Series: logarithmic infinity provider handles regular ratios with zero leading exponent");

    tests.expectEqual(
        eval(session, "series[log[x+1],{x,Infinity,4}]"),
        std::string{"seriesData[x, Infinity, {0, 1, -1/2, 1/3, -1/4}, 0, 5, 1, {{-1, 0, 0, 0, 0}}]"},
        "Series: logarithmic asymptotic composition expands regular reciprocal corrections");

    tests.expectEqual(
        eval(session, "series[log[x^2+1],{x,Infinity,4}]"),
        std::string{"seriesData[x, Infinity, {0, 0, 1, 0, -1/2}, 0, 5, 1, {{-2, 0, 0, 0, 0}}]"},
        "Series: logarithmic asymptotics retain the leading power at infinity");

    tests.expectEqual(
        eval(session, "series[log[sqrt[x]+1],{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {0, 1, -1/2, 1/3, -1/4, 1/5, -1/6}, 0, 7, 2, {{-1/2, 0, 0, 0, 0, 0, 0}}]"},
        "Series: logarithmic asymptotics compose with Puiseux powers at infinity");

    tests.expectEqual(
        eval(session, "series[log[x]^2,{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {0, 0, 0, 0}, 0, 4, 1, {{0, 0, 0, 0}, {1, 0, 0, 0}}]"},
        "Series: logarithmic layers multiply at positive infinity");

    tests.expectEqual(
        eval(session, "series[log[x]/x,{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {0, 0, 0}, 1, 4, 1, {{-1, 0, 0}}]"},
        "Series: logarithmic layers multiply reciprocal powers at infinity");

    tests.expectEqual(
        eval(session, "series[1/log[x],{x,Infinity,3}]"),
        std::string{"series[1/log[x], {x, Infinity, 3}]"},
        "Series: reciprocal logarithmic transseries remain unevaluated");

    tests.expectEqual(
        eval(session, "series[log[-x],{x,Infinity,3}]"),
        std::string{"series[log[-x], {x, Infinity, 3}]"},
        "Series: logarithmic infinity provider does not guess a negative-axis branch");

    tests.expectEqual(
        eval(session, "series[log[I*x],{x,Infinity,3}]"),
        std::string{"series[log[I x], {x, Infinity, 3}]"},
        "Series: logarithmic infinity provider does not guess a complex branch");

    kernel::KernelSession lowCostSeriesSession;
    static_cast<void>(lowCostSeriesSession.evaluate("angleMode[Rad]"));

    tests.expectEqual(
        eval(lowCostSeriesSession, "series[tan[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/3, 0, 2/15}, 1, 6, 1]"},
        "Series: tan reuses sin/cos TPSA without high-order differentiation");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[cot[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/3, 0, -1/45, 0, -2/945}, -1, 6, 1]"},
        "Series: cot reuses Laurent inversion at its simple pole");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[sec[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/2, 0, 5/24, 0}, 0, 6, 1]"},
        "Series: sec reuses cosine inversion");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[csc[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/6, 0, 7/360, 0, 31/15120}, -1, 6, 1]"},
        "Series: csc reuses sine Laurent inversion");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[tanh[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/3, 0, 2/15}, 1, 6, 1]"},
        "Series: tanh reuses sinh/cosh TPSA");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[coth[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/3, 0, -1/45, 0, 2/945}, -1, 6, 1]"},
        "Series: coth reuses hyperbolic Laurent inversion");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[expm1[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 1/2, 1/6, 1/24, 1/120}, 1, 6, 1]"},
        "Series: expm1 reuses the exponential TPSA kernel");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[log1p[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, -1/2, 1/3, -1/4, 1/5}, 1, 6, 1]"},
        "Series: log1p reuses the branch-aware logarithm kernel");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[sinc[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/6, 0, 1/120, 0}, 0, 6, 1]"},
        "Series: sinc removes its origin singularity through the trig TPSA kernel");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[cosc[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1/2, 0, -1/24, 0, 1/720}, 1, 6, 1]"},
        "Series: cosc removes its origin singularity through a half-angle identity");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[tanc[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/3, 0, 2/15, 0}, 0, 6, 1]"},
        "Series: tanc composes tangent TPSA with its removable origin singularity");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[sinhc[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, 1/6, 0, 1/120, 0}, 0, 6, 1]"},
        "Series: sinhc reuses the hyperbolic TPSA kernel");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[tanhc[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 0, -1/3, 0, 2/15, 0}, 0, 6, 1]"},
        "Series: tanhc reuses the hyperbolic quotient kernel");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[expc[x],{x,0,5}]"),
        std::string{"seriesData[x, 0, {1, 1/2, 1/6, 1/24, 1/120, 1/720}, 0, 6, 1]"},
        "Series: expc removes the expm1 origin singularity");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[log2[1+x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1/log[2], -1/(2log[2]), 1/(3log[2]), -1/(4log[2])}, 1, 5, 1]"},
        "Series: log2 uses the general two-argument logarithm rewrite");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[log10[1+x],{x,0,4}]"),
        std::string{"seriesData[x, 0, {1/log[10], -1/(2log[10]), 1/(3log[10]), -1/(4log[10])}, 1, 5, 1]"},
        "Series: log10 uses the general two-argument logarithm rewrite");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[log2[x],{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {0, 0, 0, 0}, 0, 4, 1, {{-1/log[2], 0, 0, 0}}]"},
        "Series: log2 keeps the logarithmic layer at positive infinity");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[log10[x],{x,Infinity,3}]"),
        std::string{"seriesData[x, Infinity, {0, 0, 0, 0}, 0, 4, 1, {{-1/log[10], 0, 0, 0}}]"},
        "Series: log10 keeps the logarithmic layer at positive infinity");
    tests.expectEqual(
        eval(lowCostSeriesSession, "series[1/log2[x],{x,Infinity,3}]"),
        std::string{"series[1/log[2, x], {x, Infinity, 3}]"},
        "Series: reciprocal logarithmic transseries remain unevaluated after log2 canonicalization");

    kernel::KernelSession lowCostDegreeSession;
    static_cast<void>(lowCostDegreeSession.evaluate("angleMode[Deg]"));
    tests.expectEqual(
        eval(lowCostDegreeSession, "series[tan[x],{x,0,1}]"),
        std::string{"seriesData[x, 0, {Pi/180}, 1, 2, 1]"},
        "Series: low-cost trig rewrites apply Degree-to-Radian scaling exactly once");

    tests.expectEqual(
        eval(session, "N[series[exp[x],{x,0,3}],30]"),
        std::string{"seriesData[x, 0.0, {1.0, 1.0, 0.50, 0.166666666666666666666666666667}, 0, 4, 1]"},
        "Series: N approximates coefficients while preserving structural exponent metadata");
    tests.expectEqual(
        eval(session, "D[N[series[exp[x],{x,0,3}],30],x]"),
        std::string{"seriesData[x, 0.0, {1.0, 1.0, 0.50}, 0, 3, 1]"},
        "Series: differentiated numerical SeriesData remains structurally parseable");
    tests.expectEqual(
        eval(session, "series[tan[1/x],{x,Infinity,6}]"),
        std::string{"seriesData[x, Infinity, {1, 0, 1/3, 0, 2/15, 0}, 1, 7, 1]"},
        "Series: positive-infinity reciprocal composition expands tangent without a false zero division");
    tests.expectEqual(
        eval(session, "series[sec[1/x],{x,Infinity,6}]"),
        std::string{"seriesData[x, Infinity, {1, 0, 1/2, 0, 5/24, 0, 61/720}, 0, 7, 1]"},
        "Series: positive-infinity reciprocal composition expands secant");
    tests.expectEqual(
        eval(session, "series[cot[1/x],{x,Infinity,6}]"),
        std::string{"seriesData[x, Infinity, {1, 0, -1/3, 0, -1/45, 0, -2/945, 0}, -1, 7, 1]"},
        "Series: positive-infinity reciprocal composition preserves cotangent Laurent order");
    tests.expectEqual(
        eval(session, "series[csc[1/x],{x,Infinity,6}]"),
        std::string{"seriesData[x, Infinity, {1, 0, 1/6, 0, 7/360, 0, 31/15120, 0}, -1, 7, 1]"},
        "Series: positive-infinity reciprocal composition preserves cosecant Laurent order");
    tests.expectEqual(
        eval(session, "series[sin[1/x]/cos[1/x],{x,Infinity,6}]"),
        std::string{"seriesData[x, Infinity, {1, 0, 1/3, 0, 2/15, 0}, 1, 7, 1]"},
        "Series: explicit reciprocal trig quotients share the positive-infinity normalization path");

    tests.expectEqual(
        eval(session, "normal[42]"), std::string{"42"},
        "Series: normal leaves non-series expressions unchanged");
}

} // namespace mmcal::tests
