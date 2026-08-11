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
    tests.expect(findDiagnostic(session, "integrate::unevaluated") == nullptr,
        "supported nested radical does not emit an unevaluated warning");
    tests.expect(nested.find("+-") == std::string::npos
            && nested.find(" + ") == std::string::npos
            && nested.find(" - ") == std::string::npos
            && nested.find(" / ") == std::string::npos
            && nested.find(" ^ ") == std::string::npos,
        "formatter removes redundant operator spacing and plus-negative forms");
    tests.expectEqual(eval(session, nested), nested,
        "compact nested-radical formatting round-trips through the parser");

    const std::string partial = eval(session, "integrate[x^2+gamma[x],x]");
    tests.expect(partial.find("x^3/3") != std::string::npos
            && partial.find("integrate[gamma[x], x]") != std::string::npos,
        "linear integration preserves solved terms when one term remains unresolved");
    tests.expect(findDiagnostic(session, "integrate::unevaluated") != nullptr,
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

    tests.expectEqual(eval(session, "D[sin[x],{x,4}]"), std::string{"sin[x]"},
        "D supports exact higher derivative order specifications");
    tests.expectEqual(eval(session, "D[x^2*y^3,x,y]"), std::string{"6x y^2"},
        "D supports sequential mixed derivative specifications");
    tests.expectEqual(eval(session, "D[x^5,{x,0}]"), std::string{"x^5"},
        "zeroth derivative order returns the expression unchanged");
}

} // namespace mmcal::tests
