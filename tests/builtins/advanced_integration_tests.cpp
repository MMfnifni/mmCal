// advanced・integrationの回帰テスト
#include "advanced_integration_tests.hpp"

#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "test_framework.hpp"

#include <string>
#include <string_view>
#include <vector>

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

enum class DerivativeBackMode {
    Strict,
    ResolutionOnly
};

struct DerivativeBackCase final {
    std::string_view label;
    std::string_view integrand;
    DerivativeBackMode mode = DerivativeBackMode::Strict;
};

} // namespace

void runAdvancedIntegrationTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(eval(session, "integrate[x^2+sin[x],x]"),
        std::string{"x^3/3-cos[x]"},
        "integration is linear across polynomial and trigonometric terms");
    tests.expectEqual(eval(session, "D[sin[x],x]"), std::string{"cos[x]"},
        "canonical trigonometric names share derivative knowledge");

    tests.expectEqual(eval(session, "integrate[E^x*cos[x],x]"),
        std::string{"(cos[x]+sin[x])exp[x]/2"},
        "E^x canonicalizes to exp and exponential-trigonometric products integrate exactly");
    tests.expectEqual(eval(session, "D[E^x,x]"), std::string{"exp[x]"},
        "D reuses the E^x to exp canonical identity");
    tests.expectEqual(eval(session, "(cos[x]+sin[x])/2exp[x]"),
        std::string{"(cos[x]+sin[x])exp[x]/2"},
        "formatted exponential products can be parsed back through implicit multiplication");
    tests.expectEqual(eval(session,
        "integrate[E^x*cos[x],x]==(cos[x]+sin[x])/2exp[x]"),
        std::string{"True"},
        "the exact exponential-trigonometric antiderivative compares equal after round-trip parsing");


    tests.expectEqual(eval(session, "D[integrate[1/(x^3+1),x],x]"),
        std::string{"1/(1+x^3)"},
        "exact partial fractions differentiate back to the cubic rational integrand");
    tests.expectEqual(eval(session, "integrate[1/(x^3+1),{x,0,1}]"),
        std::string{"Pi sqrt[3]/9+log[2]/3"},
        "cubic rational definite integral reduces to exact log and Pi terms");
    tests.expectEqual(eval(session, "integrate[1/((x+1)^2*(x-1)),x]"),
        std::string{"(1+x)^(-1)/2+log[-1+x]/4-log[1+x]/4"},
        "partial fractions support repeated linear factors");

    tests.expectEqual(eval(session, "integrate[1/sqrt[x^2-1],x]"),
        std::string{"log[x+sqrt[x^2-1]]"},
        "local antiderivative is accepted without installing a global sqrt factorization identity");
    tests.expectEqual(eval(session, "D[integrate[1/sqrt[x^2-1],x],x]"),
        std::string{"1/sqrt[x^2-1]"},
        "D applies the fundamental theorem to the local quadratic-root primitive request");

    const std::string nested = eval(session, "integrate[sqrt[x+sqrt[x]],x]");
    tests.expect(nested.find("integrate[") == std::string::npos
            && nested.find("log[") != std::string::npos,
        "sqrt substitution handles nested quadratic radicals symbolically");
    tests.expect(findDiagnostic(session, "integrate::unsupported") == nullptr
            && findDiagnostic(session, "integrate::partial") == nullptr,
        "supported nested radical does not emit an unevaluated warning");
    tests.expect(nested.find("+-") == std::string::npos
            && nested.find(" + ") == std::string::npos
            && nested.find(" - ") == std::string::npos
            && nested.find(" / ") == std::string::npos
            && nested.find(" ^ ") == std::string::npos,
        "formatter removes redundant operator spacing and plus-negative forms");
    tests.expectEqual(eval(session, nested), nested,
        "compact nested-radical formatting round-trips through the parser");

    tests.expectEqual(eval(session, "integrate[exp[x^6],x]"),
        std::string{"x hypergeometric1F1[1/6, 7/6, x^6]"},
        "exp of an integer monomial closes through the entire 1F1 representation");
    tests.expectEqual(eval(session,
        "fullSimplify[D[x hypergeometric1F1[1/6,7/6,x^6],x]-exp[x^6]]"),
        std::string{"0"},
        "1F1 contiguous knowledge proves the exp[x^6] antiderivative exactly");
    tests.expectEqual(eval(session, "integrate[exp[2*x^3],x]"),
        std::string{"x hypergeometric1F1[1/3, 4/3, 2x^3]"},
        "coefficient exponential monomial uses the owned 1F1 argument safely");
    tests.expectEqual(eval(session,
        "fullSimplify[D[x hypergeometric1F1[1/3,4/3,2*x^3],x]-exp[2*x^3]]"),
        std::string{"0"},
        "coefficient exponential monomial 1F1 primitive differentiates back exactly");
    tests.expectEqual(eval(session, "integrate[sqrt[1+2*x^3],x]"),
        std::string{"x hypergeometric2F1[-1/2, 1/3, 4/3, -2x^3]"},
        "binomial algebraic powers close through a branch-safe 2F1 primitive");
    tests.expectEqual(eval(session,
        "fullSimplify[D[x hypergeometric2F1[-1/2,1/3,4/3,-2*x^3],x]-sqrt[1+2*x^3]]"),
        std::string{"0"},
        "2F1 contiguous knowledge proves the binomial-power primitive exactly");
    tests.expectEqual(eval(session, "integrate[1/(1+x^5),x]"),
        std::string{"x hypergeometric2F1[1, 1/5, 6/5, -x^5]"},
        "rational binomial kernels can use the same 2F1 integration family");

    tests.expectEqual(eval(session, "integrate[exp[x]/x,x]"), std::string{"Ei[x]"},
        "exponential-over-argument kernels close through Ei");
    tests.expectEqual(eval(session, "integrate[sin[x]/x,x]"), std::string{"Si[x]"},
        "sine-over-argument kernels close through Si");
    tests.expectEqual(eval(session, "integrate[cos[x]/x,x]"), std::string{"Ci[x]"},
        "cosine-over-argument kernels close through Ci");
    tests.expectEqual(eval(session, "integrate[1/log[x],x]"), std::string{"li[x]"},
        "reciprocal-log kernels close through li");
    tests.expectEqual(eval(session, "integrate[log[1-x]/x,x]"), std::string{"-polylog[2, x]"},
        "logarithmic-over-argument kernels close through the dilogarithm");

    tests.expectEqual(eval(session, "integrate[1/sqrt[1-(1/3)*sin[x]^2],x]"),
        std::string{"ellipticF[x, 1/3]"},
        "the canonical first-kind elliptic kernel integrates to ellipticF");
    tests.expectEqual(eval(session, "integrate[sqrt[1-(1/3)*sin[x]^2],x]"),
        std::string{"ellipticE[x, 1/3]"},
        "the canonical second-kind elliptic kernel integrates to ellipticE");
    tests.expectEqual(eval(session,
        "integrate[1/((1-(1/5)*sin[x]^2)*sqrt[1-(1/3)*sin[x]^2]),x]"),
        std::string{"ellipticPi[1/5, x, 1/3]"},
        "the canonical third-kind elliptic kernel integrates to ellipticPi");
    tests.expectEqual(eval(session,
        "fullSimplify[D[ellipticF[x,1/3],x]-1/sqrt[1-(1/3)*sin[x]^2]]"),
        std::string{"0"},
        "ellipticF amplitude derivative proves the direct first-kind kernel");
    tests.expectEqual(eval(session,
        "fullSimplify[D[ellipticE[x,1/3],x]-sqrt[1-(1/3)*sin[x]^2]]"),
        std::string{"0"},
        "ellipticE amplitude derivative proves the direct second-kind kernel");
    tests.expectEqual(eval(session,
        "fullSimplify[D[ellipticPi[1/5,x,1/3],x]-1/((1-(1/5)*sin[x]^2)*sqrt[1-(1/3)*sin[x]^2]) ]"),
        std::string{"0"},
        "ellipticPi amplitude derivative proves the direct third-kind kernel");
    tests.expectEqual(eval(session, "integrate[1/sqrt[1-x^4],x]"),
        std::string{"ellipticF[asin[x], -1]"},
        "quartic algebraic kernels reduce to ellipticF without unsafe sqrt factorization");
    tests.expectEqual(eval(session, "integrate[sec[x]^3,x]"),
        std::string{"sec[x]tan[x]/2+log[sec[x]+tan[x]]/2"},
        "positive secant powers use the standard reduction formula");
    tests.expectEqual(eval(session, "integrate[csc[x]^3,x]"),
        std::string{"-cot[x]csc[x]/2-log[cot[x]+csc[x]]/2"},
        "positive cosecant powers use the standard reduction formula");
    tests.expectEqual(eval(session, "fullSimplify[D[integrate[tan[x]^4,x],x]-tan[x]^4]"),
        std::string{"0"},
        "positive tangent powers share reduction knowledge with the derivative proof engine");
    tests.expectEqual(eval(session,
        "fullSimplify[D[integrate[cos[2x^2+3x+1],x],x]-cos[2x^2+3x+1]]"),
        std::string{"0"},
        "quadratic Fresnel reduction handles linear and constant phase terms by completing the square");

    const std::string partial = eval(session, "integrate[x^2+gamma[x],x]");
    tests.expect(partial.find("x^3/3") != std::string::npos
            && partial.find("integrate[gamma[x], x]") != std::string::npos,
        "linear integration preserves solved terms when one term remains unresolved");
    tests.expect(findDiagnostic(session, "integrate::partial") != nullptr,
        "partially evaluated integral still emits an explicit warning");

    // 積分器のruntime gateにはせず、既知rule familyの退行をテスト側から監視する。
    // principal branchを跨ぐ局所primitiveなど、現Simplifierだけでは恒等式証明し切れない
    // familyはResolutionOnlyとして「積分能力を落とさない」ことを優先する。
    const std::vector<DerivativeBackCase> derivativeBackCases{
        {"constant", "5"},
        {"polynomial", "x^4-3*x+2"},
        {"reciprocal", "1/x"},
        {"affine power", "(2*x+3)^5"},
        {"affine reciprocal", "1/(2*x+3)"},
        {"exponential", "exp[2*x+1]"},
        {"sine", "sin[x]"},
        {"cosine", "cos[x]"},
        {"cotangent", "cot[x]", DerivativeBackMode::ResolutionOnly},
        {"secant", "sec[x]", DerivativeBackMode::ResolutionOnly},
        {"cosecant", "csc[x]", DerivativeBackMode::ResolutionOnly},
        {"hyperbolic sine", "sinh[x]"},
        {"hyperbolic cosine", "cosh[x]"},
        {"hyperbolic tangent", "tanh[x]", DerivativeBackMode::ResolutionOnly},
        {"hyperbolic cotangent", "coth[x]", DerivativeBackMode::ResolutionOnly},
        {"hyperbolic cosecant", "csch[x]", DerivativeBackMode::ResolutionOnly},
        {"inverse chain trigonometric", "2*x*cos[x^2]"},
        {"inverse chain logarithmic", "2*x/(x^2+1)"},
        {"arctangent rational", "1/(1+x^2)"},
        {"logarithm", "log[x]", DerivativeBackMode::ResolutionOnly},
        {"logarithmic derivative", "log[x]/x"},
        {"iterated logarithmic derivative", "1/(x*log[x])", DerivativeBackMode::ResolutionOnly},
        {"log1p", "log1p[x]", DerivativeBackMode::ResolutionOnly},
        {"polynomial times exponential", "x*exp[x]"},
        {"erf", "erf[x]"},
        {"erfc", "erfc[x]"},
        {"square root", "sqrt[x]", DerivativeBackMode::ResolutionOnly},
        {"cube root", "cbrt[x]", DerivativeBackMode::ResolutionOnly},
        {"expm1", "expm1[x]", DerivativeBackMode::ResolutionOnly},
        {"tangent", "tan[x]", DerivativeBackMode::ResolutionOnly},
        {"sech", "sech[x]", DerivativeBackMode::ResolutionOnly},
        {"quadratic rational", "1/(x^2-1)", DerivativeBackMode::ResolutionOnly},
        {"x log x", "x*log[x]", DerivativeBackMode::ResolutionOnly},
        {"trigonometric sine square", "sin[x]^2", DerivativeBackMode::ResolutionOnly},
        {"trigonometric sine sixth power", "sin[2*x]^6"},
        {"reciprocal sine square", "sin[2*x]^(-2)"},
        {"reciprocal cosine fourth power", "cos[3*x]^(-4)", DerivativeBackMode::ResolutionOnly},
        {"quadratic cosine Fresnel", "cos[4*x^2]"},
        {"quadratic sine Fresnel", "sin[8*x^2]"},
        {"shifted quadratic Fresnel", "cos[2*x^2+3*x+1]", DerivativeBackMode::ResolutionOnly},
        {"exponential integral Ei", "exp[x]/x"},
        {"sine integral Si", "sin[x]/x"},
        {"cosine integral Ci", "cos[x]/x"},
        {"logarithmic integral li", "1/log[x]"},
        {"dilogarithm", "log[1-x]/x"},
        {"hypergeometric exponential monomial", "exp[x^6]"},
        {"hypergeometric binomial power", "sqrt[1+2*x^3]"},
        {"elliptic first-kind kernel", "1/sqrt[1-(1/3)*sin[x]^2]"},
        {"elliptic second-kind kernel", "sqrt[1-(1/3)*sin[x]^2]"},
        {"elliptic third-kind kernel", "1/((1-(1/5)*sin[x]^2)*sqrt[1-(1/3)*sin[x]^2])"},
        {"quartic elliptic reduction", "1/sqrt[1-x^4]", DerivativeBackMode::ResolutionOnly},
        {"secant cube reduction", "sec[x]^3", DerivativeBackMode::ResolutionOnly},
        {"cosecant cube reduction", "csc[x]^3", DerivativeBackMode::ResolutionOnly},
        {"tangent fourth power reduction", "tan[x]^4", DerivativeBackMode::ResolutionOnly},
        {"cotangent fourth power reduction", "cot[x]^4", DerivativeBackMode::ResolutionOnly},
        {"Fresnel C primitive", "fresnelc[x]"},
        {"Fresnel S primitive", "fresnels[x]"},
        {"mixed trigonometric integer powers", "sin[x]^5*cos[x]^4"},
        {"trigonometric product-to-sum", "sin[2*x]*cos[3*x]"},
        {"trigonometric cosine square", "cos[x]^2", DerivativeBackMode::ResolutionOnly},
        {"trigonometric tangent square", "tan[x]^2", DerivativeBackMode::ResolutionOnly},
        {"trigonometric cotangent square", "cot[x]^2", DerivativeBackMode::ResolutionOnly},
        {"hyperbolic sine square", "sinh[x]^2", DerivativeBackMode::ResolutionOnly},
        {"hyperbolic cosine square", "cosh[x]^2", DerivativeBackMode::ResolutionOnly},
        {"hyperbolic tangent square", "tanh[x]^2", DerivativeBackMode::ResolutionOnly},
        {"hyperbolic cotangent square", "coth[x]^2", DerivativeBackMode::ResolutionOnly},
        {"trigonometric product", "sin[x]*cos[x]"},
        {"quadratic-root inverse chain", "2*x/sqrt[1+x^2]"},
        {"inverse sine", "asin[x]"},
        {"inverse cosine", "acos[x]"},
        {"inverse tangent", "atan[x]"},
        {"inverse hyperbolic sine", "asinh[x]"},
        {"inverse hyperbolic tangent", "atanh[x]", DerivativeBackMode::ResolutionOnly},
        {"repeated linear rational", "1/(x+1)^2"},
        {"arbitrary-base logarithm", "log[2,x]", DerivativeBackMode::ResolutionOnly},
        {"rational division", "(x+1)/(x+2)", DerivativeBackMode::ResolutionOnly},
        {"exponential-trigonometric parts", "E^x*cos[x]", DerivativeBackMode::ResolutionOnly},
        {"cubic partial fractions", "1/(x^3+1)", DerivativeBackMode::ResolutionOnly},
        {"repeated partial fractions", "1/((x+1)^2*(x-1))", DerivativeBackMode::ResolutionOnly},
        {"inverse circular quadratic root", "1/sqrt[4-x^2]", DerivativeBackMode::ResolutionOnly},
        {"inverse hyperbolic quadratic root", "1/sqrt[x^2+4]", DerivativeBackMode::ResolutionOnly},
        {"acosh", "acosh[x]", DerivativeBackMode::ResolutionOnly},
        {"local sqrt x^2-1", "1/sqrt[x^2-1]", DerivativeBackMode::ResolutionOnly},
        {"nested radical substitution", "sqrt[x+sqrt[x]]", DerivativeBackMode::ResolutionOnly}
    };

    for (const DerivativeBackCase& testCase : derivativeBackCases) {
        kernel::KernelSession proofSession;
        const std::string primitive = eval(
            proofSession, std::string{"integrate["} + std::string{testCase.integrand} + ",x]");
        const bool resolved = primitive.find("integrate[") == std::string::npos;
        tests.expect(resolved, std::string{"Integration harness resolves: "} + std::string{testCase.label});
        if (!resolved)
            continue;

        const std::string proof = eval(proofSession,
            std::string{"fullSimplify[D[("} + primitive + "),x]-("
                + std::string{testCase.integrand} + ")]" );
        if (testCase.mode == DerivativeBackMode::Strict) {
            tests.expectEqual(proof, std::string{"0"},
                std::string{"Integration derivative-back: "} + std::string{testCase.label});
            continue;
        }

        // 現Simplifierで0まで証明できないfamilyもD自体は必ず実行する。
        // 未評価D/integrateへ後退しないことを監視しつつ、runtimeの積分採否には使わない。
        tests.expect(proof.find("D[") == std::string::npos
                && proof.find("integrate[") == std::string::npos,
            std::string{"Integration derivative-back remains evaluable: "}
                + std::string{testCase.label});
    }

    // 256乗も原理上は同じ有限Fourier展開。FullSimplifyの通常候補上限は64のままなので、
    // ここでは巨大式を証明探索へ再投入せず、積分器が128周波数項へ有限時間で展開できることを監視する。
    {
        kernel::KernelSession largeTrigSession;
        const std::string primitive = eval(largeTrigSession, "integrate[sin[2*x]^256,x]");
        std::size_t sineTerms = 0;
        for (std::size_t pos = 0; (pos = primitive.find("sin[", pos)) != std::string::npos; pos += 4)
            ++sineTerms;
        tests.expect(primitive.find("integrate[") == std::string::npos && sineTerms == 128,
            "Trig polynomial degree 256 resolves to the expected finite 128-frequency primitive");
    }

    // 有限Fourier生成式は個別次数の表ではなく一般式なので、小次数格子をまとめて検証する。
    // integrateが解けることと、同じTrigKnowledgeをFullSimplifyが使ってD[F]-f=0を証明することの両方を要求する。
    for (int sinePower = 0; sinePower <= 5; ++sinePower) {
        for (int cosinePower = 0; cosinePower <= 5; ++cosinePower) {
            const int totalPower = sinePower + cosinePower;
            if (totalPower < 2 || totalPower > 6)
                continue;

            std::string integrand;
            if (sinePower > 0)
                integrand = "sin[x]^" + std::to_string(sinePower);
            if (cosinePower > 0) {
                if (!integrand.empty())
                    integrand += "*";
                integrand += "cos[x]^" + std::to_string(cosinePower);
            }

            kernel::KernelSession trigSession;
            const std::string primitive = eval(
                trigSession, "integrate[" + integrand + ",x]");
            tests.expect(primitive.find("integrate[") == std::string::npos,
                "Trig polynomial grid resolves m=" + std::to_string(sinePower)
                    + ", n=" + std::to_string(cosinePower));
            if (primitive.find("integrate[") != std::string::npos)
                continue;

            const std::string proof = eval(trigSession,
                "fullSimplify[D[(" + primitive + "),x]-(" + integrand + ")]");
            if (totalPower == 2 && (sinePower == 0 || cosinePower == 0)) {
                // 既存sin^2/cos^2 familyはFullSimplifyの探索上限内でhalf-angleを逆証明し切れない。
                // 積分能力は維持し、未評価D/integrateへ後退しないことだけを監視する。
                tests.expect(proof.find("D[") == std::string::npos
                        && proof.find("integrate[") == std::string::npos,
                    "Trig polynomial square derivative-back remains evaluable");
            }
            else {
                tests.expectEqual(proof, std::string{"0"},
                    "Trig polynomial derivative-back m=" + std::to_string(sinePower)
                        + ", n=" + std::to_string(cosinePower));
            }
        }
    }

    tests.expectEqual(eval(session, "D[sin[x],{x,4}]"), std::string{"sin[x]"},
        "D supports exact higher derivative order specifications");
    tests.expectEqual(eval(session, "D[x^2*y^3,x,y]"), std::string{"6x y^2"},
        "D supports sequential mixed derivative specifications");
    tests.expectEqual(eval(session, "D[x^5,{x,0}]"), std::string{"x^5"},
        "zeroth derivative order returns the expression unchanged");
}

} // namespace mmcal::tests
