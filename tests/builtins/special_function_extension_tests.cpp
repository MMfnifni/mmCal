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
        std::string{"2.6789385347077476337"},
        "certified Gamma encloses a general positive real value");
    tests.expectEqual(eval(session, "N[gamma[-1/3],20]"),
        std::string{"-4.0623538182792012508"},
        "certified Gamma uses reflection on the negative real axis");
    tests.expectEqual(eval(session, "N[gamma[1+I],30]"),
        std::string{"0.498015668118356042713691117462-0.154949828301810685124955130484I"},
        "Gamma has a certified principal complex backend");
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
    tests.expectEqual(eval(session, "N[erf[1+I],30]"),
        std::string{"1.31615128169794764488027108024+0.190453469237834686284108861969I"},
        "erf has a certified entire complex series backend");
    tests.expectEqual(eval(session, "N[erfc[1+I],30]"),
        std::string{"-0.316151281697947644880271080244-0.190453469237834686284108861969I"},
        "erfc reuses the certified complex erf enclosure");
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
    tests.expectEqual(eval(session, "N[fresnelc[1+I],30]"),
        std::string{"2.55579377810243902463452238835+2.55579377810243902463452238835I"},
        "Fresnel C has a certified entire complex series backend");
    tests.expectEqual(eval(session, "N[fresnels[1+I],30]"),
        std::string{"-2.06188821919484046808071653669+2.06188821919484046808071653669I"},
        "Fresnel S has a certified entire complex series backend");
    tests.expectEqual(eval(session, "N[fresnelc[7+I],20]"),
        std::string{"10840706.412719365363+79376633.864722432483I"},
        "Fresnel C remains certified in the former complex series-boundary region");
    tests.expectEqual(eval(session, "N[fresnels[7+I],20]"),
        std::string{"-79376633.364722432495+10840705.912719365365I"},
        "Fresnel S shares the optimized certified complex series backend");
    tests.expectEqual(eval(session, "N[fresnelc[8+I],20]"),
        std::string{"-1613473902.7057768968+193865931.68895630904I"},
        "complex Fresnel no longer has the historical |z|=8 work boundary");
    tests.expectEqual(eval(session, "N[fresnelc[32+I],20]"),
        std::string{"-227144419018168036370000000000000000000000.0+7027787067654522305300000000000000000000.0I"},
        "complex Fresnel uses the certified asymptotic wedge for large near-axis arguments");
    tests.expectEqual(eval(session, "N[fresnels[1+32I],20]"),
        std::string{"227144419018168036370000000000000000000000.0+7027787067654522305300000000000000000000.0I"},
        "complex Fresnel quarter-turn symmetry shares the asymptotic wedge across axes");
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
        std::string{"1.1920688079818883008"},
        "1F1 has a certified exact-Rational real series backend");
    tests.expectEqual(eval(session, "D[hypergeometric1F1[1/6,7/6,x],x]"),
        std::string{"hypergeometric1F1[7/6, 13/6, x]/7"},
        "1F1 derivative uses the exact contiguous derivative identity");
    tests.expectEqual(eval(session, "N[hypergeometric1F1[1/2,5/4,1+I],30]"),
        std::string{"1.29886920866742010665179517472+0.728624123036746834338241315805I"},
        "1F1 has a certified entire complex series backend");
    tests.expectEqual(eval(session, "N[hypergeometric1F1[1,2,N[1/2,5]],20]"),
        std::string{"1.2974"},
        "1F1 projects a finite-precision all-real complex enclosure back to a real result");
    tests.expectEqual(eval(session, "hypergeometric1F1[0,0,2+I]"), std::string{"1"},
        "terminating 1F1 defined before a denominator pole is independent of a complex z");
    tests.expectEqual(eval(session, "hypergeometric1F1[0,3,1/x]"),
        std::string{"hypergeometric1F1[0, 3, 1/x]"},
        "1F1 zero-parameter degeneration preserves an undefined disappearing argument");
    tests.expectEqual(eval(session, "simplify[hypergeometric1F1[0,3,1/x],x!=0]"),
        std::string{"1"},
        "1F1 zero-parameter degeneration resumes when the disappearing argument is proven defined");
    tests.expectEqual(eval(session, "hypergeometric1F1[-1,3,2+I]"),
        std::string{"1/3-I/3"},
        "terminating 1F1 polynomials evaluate exact complex-number arguments without a real-only dispatch");

    tests.expectEqual(eval(session, "hypergeometric2F1[-2,1,3,1/2]"), std::string{"17/24"},
        "terminating 2F1 polynomial is evaluated exactly");
    tests.expectEqual(eval(session, "hypergeometric2F1[0,2,3,x]"), std::string{"1"},
        "2F1 with a zero numerator parameter reduces before numerical evaluation");
    tests.expectEqual(eval(session, "hypergeometric2F1[0,2,3,1/x]"),
        std::string{"hypergeometric2F1[0, 2, 3, 1/x]"},
        "2F1 zero-parameter degeneration preserves an undefined disappearing argument");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/2,3/2,1/4],20]"),
        std::string{"1.0471975511965977462"},
        "2F1 has a certified real backend inside its principal Gauss-series disk");
    tests.expectEqual(eval(session, "D[hypergeometric2F1[1/2,1/3,5/4,x],x]"),
        std::string{"2hypergeometric2F1[3/2, 4/3, 9/4, x]/15"},
        "2F1 differentiates with respect to its argument when parameters are constant");
    tests.expectEqual(eval(session, "D[hypergeometric2F1[x,1,2,x],x]"),
        std::string{"D[hypergeometric2F1[x, 1, 2, x], x]"},
        "2F1 parameter derivatives remain unresolved instead of applying an incomplete z-only rule");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/3,5/4,1/2+I/4],30]"),
        std::string{"1.07686246822300103290787411965+0.0572588164342812322158675308896I"},
        "2F1 has a certified complex Gauss-series backend inside the unit disk");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1,1,2,N[1/2,5]],20]"),
        std::string{"1.3863"},
        "2F1 keeps finite-precision real inputs real below the principal branch cut");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[3.4,5.6,4+I,4.6+2I],20]"),
        std::string{"0.0046136876612922014955+0.0019659119401108294965I"},
        "2F1 uses a certified principal 1/z connection outside the unit disk when nondegenerate");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[3.4,5.6,4+I,4.6+2I],40]"),
        std::string{"0.004613687661292201495472361244961975817392+0.001965911940110829496460833819409636768435I"},
        "2F1 principal 1/z continuation remains certified at higher precision");
    tests.expectEqual(eval(session, "hypergeometric2F1[-1,2,3,2+I]"),
        std::string{"-1/3-2I/3"},
        "terminating 2F1 polynomials bypass continuation and evaluate exact complex-number arguments");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/3,5/4,2],20]"),
        std::string{"1.0930912348251107981-0.47781127033220940796I"},
        "2F1 exact real z on the principal cut uses the same principal continuation convention as complex z");
    tests.expectEqual(eval(session, "N[hypergeometric1F1[1/2,5/4,160],20]"),
        std::string{"34978276269386186023000000000000000000000000000000000000000000000000.0"},
        "1F1 remains certified at the bounded-work argument threshold");
    tests.expectEqual(eval(session, "N[hypergeometric1F1[1/2,5/4,161],20]"),
        std::string{"94636147547763363460000000000000000000000000000000000000000000000000.0"},
        "1F1 no longer imposes the former |z| <= 160 work boundary");
    tests.expectEqual(eval(session, "N[hypergeometric1F1[1/2,5/4,-512],20]"),
        std::string{"0.032697046018868862934"},
        "1F1 keeps large negative real arguments certified without exact-Rational blow-up");
    tests.expectEqual(eval(session, "precision[N[hypergeometric1F1[1/2,5/4,256+I],20]]"),
        std::string{"19"},
        "complex 1F1 remains precision-carrying well beyond the former magnitude boundary");
    tests.expectEqual(eval(session, "precision[N[hypergeometric1F1[1/2,5/4,N[161,5]],20]]"),
        std::string{"2"},
        "1F1 does not recover hidden precision when a finite-precision argument exceeds the former boundary");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/3,5/4,9/10],20]"),
        std::string{"1.2540597304760187485"},
        "2F1 remains certified at the bounded-work Gauss-series threshold");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/3,5/4,19/20],20]"),
        std::string{"1.3068052887857217891"},
        "2F1 certifies the full open unit disk instead of imposing the former 9/10 work boundary");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/3,5/4,99/100],20]"),
        std::string{"1.3918605978010497470"},
        "2F1 remains certified close to the unit-circle convergence boundary");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/3,5/4,1],20]"),
        std::string{"1.4908745724194916165"},
        "2F1 uses Gauss summation at z=1 when Re(c-a-b)>0");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/3,1/2,1],20]"),
        std::string{"hypergeometric2F1[1/2, 1/3, 1/2, 1]"},
        "2F1 keeps unsupported unit-circle boundary cases unevaluated");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/3,1/10^20,9/10],20]"),
        std::string{"131575556603428593630.0"},
        "2F1 certifies a large real value near a denominator-parameter pole without Rational blow-up");
    tests.expectEqual(eval(session, "N[hypergeometric2F1[1/2,1/3,1/10^20,1/2+I/4],20]"),
        std::string{"10157290331004504305.0+13137475151258192356.0I"},
        "complex 2F1 uses a relative tail bound for large values near a denominator-parameter pole");

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
    tests.expectEqual(eval(session, "N[ellipticF[N[1/2,5],1/3],20]"),
        std::string{"0.50685"},
        "ellipticF accepts finite-precision real interval inputs without machine fallback");
    tests.expectEqual(eval(session, "N[ellipticE[N[1/2,5],1/3],20]"),
        std::string{"0.49332"},
        "ellipticE accepts finite-precision real interval inputs without machine fallback");
    tests.expectEqual(eval(session, "N[ellipticPi[1/4,N[1/2,5],1/3],20]"),
        std::string{"0.51737"},
        "ellipticPi accepts finite-precision real interval inputs without machine fallback");
    tests.expectEqual(eval(session, "N[ellipticF[1/2,9/10],20]"),
        std::string{"0.51976394243782869496"},
        "ellipticF remains certified at the bounded-work parameter threshold");
    tests.expectEqual(eval(session, "N[ellipticE[1/2,9/10],20]"),
        std::string{"0.48155710744262109430"},
        "ellipticE remains certified at the bounded-work parameter threshold");
    tests.expectEqual(eval(session, "N[ellipticPi[9/10,1/2,1/3],20]"),
        std::string{"0.54883454168831276122"},
        "ellipticPi remains certified at the characteristic work threshold");
    tests.expectEqual(eval(session, "N[ellipticF[3/2,9/10],20]"),
        std::string{"2.3558627383594486071"},
        "ellipticF amplitude-aware tail proof remains certified near Pi/2");
    tests.expectEqual(eval(session, "N[ellipticE[3/2,9/10],20]"),
        std::string{"1.0822199401249235350"},
        "ellipticE amplitude-aware tail proof remains certified near Pi/2");
    tests.expectEqual(eval(session, "N[ellipticPi[9/10,3/2,1/3],20]"),
        std::string{"4.9290739403905673975"},
        "ellipticPi amplitude-aware tail proof remains certified near Pi/2");
    tests.expectEqual(eval(session, "N[ellipticF[2,9/10],20]"),
        std::string{"3.7114432264647167179"},
        "ellipticF keeps the conservative tail proof beyond the monotone amplitude region");
    tests.expectEqual(eval(session, "N[ellipticF[1/2,19/20],20]"),
        std::string{"0.52099294480389383536"},
        "ellipticF crosses the former 9/10 work boundary through the Carlson backend");
    tests.expectEqual(eval(session, "N[ellipticE[1/2,19/20],20]"),
        std::string{"0.48049357622814572968"},
        "ellipticE crosses the former 9/10 work boundary through the Carlson backend");
    tests.expectEqual(eval(session, "N[ellipticPi[19/20,1/2,1/3],20]"),
        std::string{"0.55154859568790941306"},
        "ellipticPi crosses the former characteristic work boundary through RJ");
    tests.expectEqual(eval(session, "N[ellipticF[1/2,99/100],20]"),
        std::string{"0.52198775871658283077"},
        "ellipticF remains certified close to the complete-parameter boundary");
    tests.expectEqual(eval(session, "N[ellipticE[1/2,99/100],20]"),
        std::string{"0.47963950999048385681"},
        "ellipticE remains certified close to the complete-parameter boundary");
    tests.expectEqual(eval(session, "N[ellipticPi[99/100,1/2,99/100],20]"),
        std::string{"0.57149053861428088394"},
        "ellipticPi certifies simultaneous near-one m and n through RF/RJ");
    tests.expectEqual(eval(session, "N[ellipticF[2,19/20],20]"),
        std::string{"4.3357836260680930745"},
        "ellipticF Carlson backend preserves real-period reduction");
    tests.expectEqual(eval(session, "N[ellipticPi[19/20,2,1/3],20]"),
        std::string{"14.335569012614443776"},
        "ellipticPi RJ path preserves real-period reduction without RC cancellation loss");
    tests.expectEqual(eval(session, "N[ellipticF[Pi/2,19/20],20]"),
        std::string{"2.9083372484445521001"},
        "ellipticF reduces exact rational Pi multiples without interval quotient ambiguity");
    tests.expectEqual(eval(session, "N[ellipticE[Pi/2,19/20],20]"),
        std::string{"1.0604737277662782427"},
        "ellipticE reduces exact rational Pi multiples without interval quotient ambiguity");
    tests.expectEqual(eval(session, "N[ellipticE[Pi/2,1],20]"),
        std::string{"1.0"},
        "ellipticE keeps the finite exact m=1 complete-boundary degeneration");
    tests.expectEqual(eval(session, "N[ellipticE[3Pi/2,1],20]"),
        std::string{"3.0"},
        "ellipticE m=1 degeneration preserves real period accumulation");
    tests.expectEqual(eval(session, "N[ellipticPi[19/20,Pi/2,1/3],20]"),
        std::string{"8.2772815079796133126"},
        "ellipticPi reduces exact rational Pi multiples without interval quotient ambiguity");
    tests.expectEqual(eval(session, "N[ellipticPi[19/20,3Pi/2,1/3],20]"),
        std::string{"24.831844523938839938"},
        "ellipticPi reuses the complete value for exact half-integer pi amplitudes");
    tests.expectEqual(eval(session, "N[ellipticF[1/2,2],20]"),
        std::string{"0.55135887907967981413"},
        "ellipticF accepts m greater than one while the integration path stays on the real branch");
    tests.expectEqual(eval(session, "N[ellipticPi[2,1/2,1/3],20]"),
        std::string{"0.62289042574026551304"},
        "ellipticPi accepts n greater than one before the first pole");
    tests.expectEqual(eval(session, "N[ellipticF[4/5,2],20]"),
        std::string{"ellipticF[4/5, 2]"},
        "ellipticF does not force a real value after the principal path reaches a complex branch");
    tests.expectEqual(eval(session, "N[ellipticPi[2,4/5,1/3],20]"),
        std::string{"ellipticPi[2, 4/5, 1/3]"},
        "ellipticPi does not cross a real pole without a principal-value backend");
    tests.expectEqual(eval(session, "N[ellipticF[1/2,N[19/20,5]],20]"),
        std::string{"0.52099"},
        "ellipticF Carlson evaluation preserves the finite-precision parameter information floor");
    tests.expectEqual(eval(session, "N[ellipticPi[N[19/20,5],1/2,1/3],20]"),
        std::string{"0.55155"},
        "ellipticPi Carlson evaluation preserves the finite-precision characteristic information floor");
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
    tests.expectEqual(eval(session, "li[0]"), std::string{"0"},
        "li has the exact principal value zero at the origin");
    tests.expectEqual(eval(session, "N[Ei[1],20]"),
        std::string{"1.8951178163559367555"},
        "Ei has a certified real backend near the origin");
    tests.expectEqual(eval(session, "N[Si[1],20]"),
        std::string{"0.94608307036718301494"},
        "Si has a certified real backend");
    tests.expectEqual(eval(session, "N[Ci[1],20]"),
        std::string{"0.33740392290096813466"},
        "Ci has a certified positive-real backend");
    tests.expectEqual(eval(session, "N[li[2],20]"),
        std::string{"1.0451637801174927848"},
        "li reuses the certified log and Ei backends");
    tests.expectEqual(eval(session, "N[li[-2],20]"),
        std::string{"0.035532275913560467787+3.7351504825512546398I"},
        "negative-real li uses the principal complex Log/Ei path instead of the positive-real backend");
    tests.expectEqual(eval(session, "N[li[-2+I],20]"),
        std::string{"0.35983768206786498306+3.8175912081992291774I"},
        "li supports principal complex certified evaluation away from its Log branch cut");
    tests.expectEqual(eval(session, "N[li[-2+I],40]"),
        std::string{"0.3598376820678649830641857502620169916140+3.817591208199229177430772175717100104837I"},
        "complex li remains certified when the Log/Arg path is evaluated at higher precision");
    tests.expectEqual(eval(session, "N[Ei[1+I],30]"),
        std::string{"1.76462598556385406842673816135+2.38776985151052241926279208910I"},
        "Ei has a certified principal complex series backend away from its cut");
    tests.expectEqual(eval(session, "N[Si[1+I],30]"),
        std::string{"1.10422265823558173955875396985+0.882453805007917743376124044695I"},
        "Si has a certified entire complex series backend");
    tests.expectEqual(eval(session, "N[Ci[1+I],30]"),
        std::string{"0.882172180555936325050614116656+0.287249133519955939527283572386I"},
        "Ci has a certified principal complex series backend away from its cut");
    tests.expectEqual(eval(session, "N[Ci[120+I],20]"),
        std::string{"0.0074445393237725791671+0.0079585143427773070995I"},
        "complex Ci requires whole-value relative precision instead of per-component zero rounding");
    tests.expectEqual(eval(session, "N[Ci[140+I],20]"),
        std::string{"0.01080710226306548880-0.0016788715636035789868I"},
        "complex Ci crosses the former |z| <= 128 boundary through the certified asymptotic backend");
    tests.expectEqual(eval(session, "N[Ei[513I],20]"),
        std::string{"-0.0015490370545520353789+3.1427759880886742510I"},
        "complex Ei crosses the former |z| <= 512 boundary without reverting to the cancellation-heavy series");
    tests.expectEqual(eval(session, "N[Ei[1000+I]/exp[1000+I],20]"),
        std::string{"0.0010010010030130653917-0.0000010020050201016102752I"},
        "complex Ei large-argument E1 continuation preserves relative precision");
    tests.expectEqual(eval(session, "N[Ci[1000+I],20]"),
        std::string{"0.0012757330374975137106+0.00066060412280668393821I"},
        "complex Ci large-argument E1 relation certifies the right half-plane");
    tests.expectEqual(eval(session, "N[Ci[-1000+I],20]"),
        std::string{"0.0012757330374975137106+3.1409320494669865545I"},
        "complex Ci left-half-plane continuation preserves the principal Log branch offset");
    tests.expectEqual(eval(session, "N[Ei[97],20]"),
        std::string{"13942532424263388810000000000000000000000.0"},
        "Ei no longer imposes the former real |x| <= 96 work boundary");
    tests.expectEqual(eval(session, "N[Si[97],20]"),
        std::string{"1.5802915860056074271"},
        "Si no longer imposes the former real |x| <= 96 work boundary");
    tests.expectEqual(eval(session, "N[Ci[97],20]"),
        std::string{"0.0040109142844998040983"},
        "Ci no longer imposes the former real x <= 96 work boundary");
    tests.expectEqual(eval(session, "N[Ei[-64],20]"),
        std::string{"-0.0000000000000000000000000000024679685594526945427"},
        "negative Ei preserves significant digits instead of collapsing to 0.0");
    tests.expectEqual(eval(session, "N[Si[512],20]"),
        std::string{"1.5727429488260625276"},
        "Si remains certified well beyond the former real-series boundary");
    tests.expectEqual(eval(session, "N[Ci[512],20]"),
        std::string{"0.00015911090433721849048"},
        "Ci remains certified well beyond the former real-series boundary");
    tests.expectEqual(eval(session, "N[Si[10000],20]"),
        std::string{"1.5708915453859619157"},
        "Si large-positive asymptotic backend certifies values far beyond the former boundary");
    tests.expectEqual(eval(session, "N[Ci[10000],20]"),
        std::string{"-0.000030551916724485212665"},
        "Ci large-positive asymptotic backend avoids EulerGamma cancellation at large arguments");
    tests.expectEqual(eval(session, "N[Ei[1000]/exp[1000],20]"),
        std::string{"0.0010010020060241207251"},
        "positive Ei asymptotic backend preserves relative precision at large arguments");
    tests.expectEqual(eval(session, "N[Ei[10000]/exp[10000],20]"),
        std::string{"0.00010001000200060024012"},
        "positive Ei asymptotic backend scales without a fixed magnitude threshold");
    tests.expectEqual(eval(session, "N[Ei[-1+I],20]"),
        std::string{"-0.00028162445198141832551+2.9622681185504342983I"},
        "complex Ei remains stable in the left half-plane away from the principal cut");
    tests.expectEqual(eval(session, "N[Ci[-1+I],20]"),
        std::string{"0.88217218055593632505+2.8543435200698372989I"},
        "complex Ci remains stable in the left half-plane away from the principal cut");
    tests.expectEqual(eval(session, "D[Ei[x],x]"), std::string{"exp[x]/x"},
        "Ei derivative is exact");
    tests.expectEqual(eval(session, "D[Si[x],x]"), std::string{"sinc[x Rad]"},
        "Si derivative preserves the removable singularity and explicit radian semantics");
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
    tests.expectEqual(eval(session, "polylog[3,1]"), std::string{"zeta[3]"},
        "positive integer polylog at one reduces exactly to zeta");
    tests.expectEqual(eval(session, "polylog[3,-1]"), std::string{"-3zeta[3]/4"},
        "positive integer polylog at minus one reduces exactly to eta and zeta");
    tests.expectEqual(eval(session, "N[polylog[2,1/2],20]"),
        std::string{"0.58224052646501250590"},
        "polylog has a certified |z|<1 real series backend for positive integer order");
    tests.expectEqual(eval(session, "D[polylog[2,x],x]"),
        std::string{"cases[-log[1-x]/x if x != 0; 1 if x == 0]"},
        "dilogarithm derivative preserves its removable singularity at zero");
    tests.expectEqual(eval(session, "D[polylog[3,x],x]"),
        std::string{"cases[polylog[2, x]/x if x != 0; 1 if x == 0]"},
        "general polylog derivative lowers the order while preserving the value at zero");
    tests.expectEqual(eval(session, "N[polylog[2,1/2+I/4],30]"),
        std::string{"0.545867504964079626756361225341+0.339137699239769040815315138383I"},
        "polylog has a certified complex series backend inside the unit disk");
    tests.expectEqual(eval(session, "N[polylog[2,N[1/2,5]],20]"),
        std::string{"0.58224"},
        "polylog reuses the complex interval series for finite-precision real inputs and returns a real enclosure");
    tests.expectEqual(eval(session, "N[polylog[2,2],20]"), std::string{"polylog[2, 2]"},
        "dilogarithm keeps an exact point on the principal positive-real cut unevaluated");
    tests.expectEqual(eval(session, "N[polylog[2,49/50],20]"),
        std::string{"1.5457997120314656097"},
        "polylog remains certified across the former bounded-work threshold");
    tests.expectEqual(eval(session, "N[polylog[2,99/100],20]"),
        std::string{"1.5886254480763753270"},
        "dilogarithm uses a certified reflection backend near one");
    tests.expectEqual(eval(session, "N[polylog[2,999/1000],20]"),
        std::string{"1.6370226052761177427"},
        "dilogarithm remains fast and certified close to the unit point");
    tests.expectEqual(eval(session, "N[polylog[3,999/1000],20]"),
        std::string{"1.2004153539954643452"},
        "higher positive integer polylog uses the certified near-one logarithmic expansion");
    tests.expectEqual(eval(session, "N[polylog[3,N[999/1000,8]],20]"),
        std::string{"1.2004154"},
        "near-one polylog preserves finite input information instead of recovering hidden precision");
    tests.expectEqual(eval(session, "N[polylog[2,-2],20]"),
        std::string{"-1.4367463668836809464"},
        "dilogarithm maps the negative real axis through the DLMF connection formula");
    tests.expectEqual(eval(session, "N[polylog[2,I],20]"),
        std::string{"-0.20561675835602830456+0.91596559417721901505I"},
        "dilogarithm certifies a representative unit-circle point by continuation");
    tests.expectEqual(eval(session, "N[polylog[2,2+I],20]"),
        std::string{"1.1866885370000578311+2.4077407693457720017I"},
        "dilogarithm inversion certifies complex values outside the unit disk away from the cut");

    tests.expectEqual(eval(session, "beta[2,3]"), std::string{"1/12"},
        "Beta at positive integers is exact");
    tests.expectEqual(eval(session, "beta[1/2,1/2]"), std::string{"Pi"},
        "exact Gamma values simplify Beta one-half exactly");
    tests.expectEqual(eval(session, "beta[1/2,1]"), std::string{"2"},
        "Beta(x,1)=1/x is exact on the positive-real domain");
    tests.expectEqual(eval(session, "betaln[1/2,1/2]"), std::string{"log[Pi]"},
        "betaln preserves an exact logarithmic result");
    tests.expectEqual(eval(session, "N[beta[1/3,2/3],20]"),
        std::string{"3.6275987284684357012"},
        "Beta has a certified positive-real numerical backend");
    tests.expect(evalError(session, "beta[-1/2,2]").type() == error::CalcErrorType::Domain,
        "current Beta contract rejects non-positive real arguments");
    tests.expectEqual(eval(session, "D[beta[x,3],x]"),
        std::string{"(digamma[x]-digamma[x+3])beta[x, 3]"},
        "Beta differentiates through its logarithmic Gamma derivative");
    tests.expectEqual(eval(session, "D[betaln[x,3],x]"),
        std::string{"digamma[x]-digamma[x+3]"},
        "betaln differentiates through digamma without expanding Gamma quotients");

    tests.expectEqual(eval(session, "zeta[0]"), std::string{"-1/2"},
        "zeta zero is exact");
    tests.expectEqual(eval(session, "zeta[-2]"), std::string{"0"},
        "zeta negative even integers are exact trivial zeros");
    tests.expectEqual(eval(session, "zeta[2]"), std::string{"Pi^2/6"},
        "zeta two uses the exact Basel value");
    tests.expectEqual(eval(session, "N[zeta[3],20]"),
        std::string{"1.2020569031595942854"},
        "zeta has a certified real backend for s greater than one");
    tests.expectEqual(eval(session, "N[zeta[3/2],20]"),
        std::string{"2.6123753486854883433"},
        "zeta keeps non-integer real Euler-Maclaurin evaluation certified");
    tests.expectEqual(eval(session, "N[zeta[12/5],20]"),
        std::string{"1.3833428588407357282"},
        "zeta keeps general rational Euler-Maclaurin evaluation certified");
    tests.expectEqual(eval(session, "N[zeta[2+I],30]"),
        std::string{"1.15035570325490267174284993474-0.437530865919607881117527898593I"},
        "zeta has a certified complex Euler-Maclaurin backend for Re(s)>1");
    tests.expectEqual(eval(session, "N[zeta[-1/2],20]"),
        std::string{"-0.20788622497735456602"},
        "zeta uses the exact functional equation before certified evaluation on negative rationals");
    tests.expectEqual(eval(session, "N[zeta[1/2],20]"),
        std::string{"-1.4603545088095868129"},
        "zeta certifies the real critical strip by generalized Euler-Maclaurin evaluation");
    tests.expectEqual(eval(session, "N[zeta[1/2+I],20]"),
        std::string{"0.14393642707718906032-0.72209974353167308913I"},
        "zeta certifies complex values in the critical strip");
    tests.expectEqual(eval(session, "N[zeta[-3+I],20]"),
        std::string{"0.014382512185224971007+0.010349659644311743514I"},
        "zeta uses the functional equation for the complex left half-plane");
    tests.expect(evalError(session, "zeta[1]").type() == error::CalcErrorType::Domain,
        "zeta rejects its pole at one");

    tests.expectEqual(eval(session, "N[digamma[1],20]"),
        std::string{"-0.57721566490153286061"},
        "digamma has a certified positive-real backend");
    tests.expectEqual(eval(session, "N[trigamma[1],20]"),
        std::string{"1.6449340668482264365"},
        "trigamma has a certified positive-real backend");
    tests.expectEqual(eval(session, "N[digamma[1/3],20]"),
        std::string{"-3.1320337800208063230"},
        "digamma preserves an exact rational argument in the certified real backend");
    tests.expectEqual(eval(session, "N[trigamma[1/3],20]"),
        std::string{"10.095597125427094082"},
        "trigamma preserves an exact rational argument in the certified real backend");
    tests.expectEqual(eval(session, "N[digamma[1+I],30]"),
        std::string{"0.0946503206224769772718784827219+1.07667404746858117413405079475I"},
        "digamma has a certified complex recurrence plus asymptotic backend");
    tests.expectEqual(eval(session, "N[trigamma[1+I],30]"),
        std::string{"0.463000096622763786298326518184-0.794233542759318865583013617157I"},
        "trigamma has a certified complex Euler-Maclaurin backend");
    tests.expectEqual(eval(session, "N[digamma[2+I]-digamma[1+I]-1/(1+I),20]"),
        std::string{"0"},
        "exact-complex digamma recurrence cancels to an exact-source zero display");
    tests.expectEqual(eval(session, "N[trigamma[2+I]-trigamma[1+I]+1/(1+I)^2,20]"),
        std::string{"0"},
        "exact-complex trigamma recurrence cancels to an exact-source zero display");
    tests.expectEqual(eval(session, "N[digamma[-1/2],20]"),
        std::string{"0.036489973978576520559"},
        "digamma shifts negative noninteger rationals into the certified positive domain exactly");
    tests.expectEqual(eval(session, "N[trigamma[-1/2],20]"),
        std::string{"8.9348022005446793094"},
        "trigamma shifts negative noninteger rationals into the certified positive domain exactly");
    tests.expectEqual(eval(session, "N[digamma[N[-1/2,5]],20]"),
        std::string{"0.0365"},
        "finite-precision negative real digamma uses the real recurrence backend safely");
    tests.expectEqual(eval(session, "N[trigamma[N[-1/2,5]],20]"),
        std::string{"8.9348"},
        "finite-precision negative real trigamma uses the complex recurrence backend safely");
    tests.expectEqual(eval(session, "N[Ci[N[-1/2,5]],20]"),
        std::string{"-0.17778+3.1416I"},
        "finite-precision negative real Ci returns the principal complex value");
    tests.expectEqual(eval(session, "trigamma[2]"), std::string{"Pi^2/6-1"},
        "trigamma at positive integers reduces exactly to a harmonic correction");
    tests.expect(evalError(session, "digamma[0]").type() == error::CalcErrorType::Domain,
        "digamma rejects non-positive integer poles");
    tests.expect(evalError(session, "trigamma[-1]").type() == error::CalcErrorType::Domain,
        "trigamma rejects non-positive integer poles");
    tests.expectEqual(eval(session, "D[gamma[x],x]"), std::string{"digamma[x]gamma[x]"},
        "Gamma derivative uses digamma");
    tests.expectEqual(eval(session, "D[lgamma[x],x]"), std::string{"digamma[x]"},
        "LogGamma derivative uses digamma");
    tests.expectEqual(eval(session, "D[digamma[x],x]"), std::string{"trigamma[x]"},
        "digamma derivative uses trigamma");

    tests.expectEqual(eval(session, "ibeta[1,1,1/4]"), std::string{"1/4"},
        "regularized incomplete Beta reduces exactly for unit parameters");
    tests.expectEqual(eval(session, "ibeta[2,3,1/2]"), std::string{"11/16"},
        "regularized incomplete Beta is exact for positive integer parameters");
    tests.expectEqual(eval(session, "N[ibeta[1/3,2/3,1/4],20]"),
        std::string{"0.53302858123542523627"},
        "regularized incomplete Beta has a certified real backend");
    tests.expectEqual(eval(session, "N[ibeta[N[1/3,8],N[2/3,8],1/4],8]"),
        std::string{"0.53302858"},
        "ibeta propagates finite-precision positive parameter enclosures without recovering hidden exact values");
    tests.expect(evalError(session, "ibeta[1,1,2]").type() == error::CalcErrorType::Domain,
        "ibeta rejects x outside the real unit interval");
    tests.expectEqual(eval(session, "D[ibeta[2,3,x],x]"),
        std::string{"x*(1-x)^2/beta[2, 3]"},
        "ibeta differentiates with respect to x when its parameters are constant");

    tests.expectEqual(eval(session, "binom[1/2,2]"), std::string{"-1/8"},
        "generalized binomial with nonnegative integer order is exact");
    tests.expectEqual(eval(session, "fallingfact[5,3]"), std::string{"60"},
        "falling factorial is exact");
    tests.expectEqual(eval(session, "risingfact[5,3]"), std::string{"210"},
        "rising factorial is exact");
    tests.expectEqual(eval(session, "D[binom[x,3],x]"), std::string{"1/3-x+x^2/2"},
        "finite generalized binomial order differentiates through its exact polynomial");
    tests.expectEqual(eval(session, "D[fallingfact[x,3],x]"), std::string{"3x^2-6x+2"},
        "finite falling factorial differentiates through its exact polynomial");
    tests.expectEqual(eval(session, "D[risingfact[x,3],x]"), std::string{"3x^2+6x+2"},
        "finite rising factorial differentiates through its exact polynomial");

    tests.expectEqual(eval(session, "sqrt[Pi]^2"), std::string{"Pi"},
        "squaring a principal square root is branch-safe and exact");
}

} // namespace mmcal::tests
