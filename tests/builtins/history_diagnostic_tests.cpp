// history・diagnosticの回帰テスト
#include "history_diagnostic_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "numeric/decimal_approximation.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

[[nodiscard]] std::string eval(kernel::KernelSession& session, std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

[[nodiscard]] error::CalcError evalError(kernel::KernelSession& session, std::string_view source) {
    try { static_cast<void>(session.evaluate(source)); }
    catch (const error::CalcError& exception) { return exception; }
    throw std::logic_error("Expected CalcError");
}

[[nodiscard]] bool hasWarning(const kernel::KernelSession& session, std::string_view code) {
    for (const auto& diagnostic : session.diagnostics())
        if (diagnostic.severity == evaluation::DiagnosticSeverity::Warning
            && diagnostic.code == code)
            return true;
    return false;
}

[[nodiscard]] std::size_t warningCount(
    const kernel::KernelSession& session,
    std::string_view code) {
    std::size_t count = 0;
    for (const auto& diagnostic : session.diagnostics())
        if (diagnostic.severity == evaluation::DiagnosticSeverity::Warning
            && diagnostic.code == code)
            ++count;
    return count;
}

} // namespace

void runHistoryDiagnosticTests(TestRunner& tests) {
    kernel::KernelSession history;
    tests.expectEqual(eval(history, "2+3"), std::string{"5"},
        "first output is stored");
    tests.expectEqual(eval(history, "10*10"), std::string{"100"},
        "second output is stored");
    tests.expectEqual(eval(history, "In[1]"), std::string{"5"},
        "In[n] re-evaluates the stored input expression");
    tests.expectEqual(eval(history, "Out[2]"), std::string{"100"},
        "Out[n] returns the absolute output entry");
    tests.expectEqual(eval(history, "%"), std::string{"100"},
        "percent history remains relative to successful outputs");

    tests.expectEqual(eval(history, "In[-1]"), std::string{"100"},
        "In[-1] re-evaluates the immediately previous input");
    tests.expectEqual(eval(history, "Out[-1]"), std::string{"100"},
        "Out[-1] returns the most recent successful output");
    tests.expectEqual(eval(history, "@"), std::string{"100"},
        "at shorthand re-evaluates the immediately previous input");

    kernel::KernelSession symbolicHistory;
    tests.expectEqual(eval(symbolicHistory, "D[E^x cos[x],x]"),
        std::string{"cos[x]exp[x]-exp[x]sin[x]"},
        "symbolic history seed stores an exact derivative");
    tests.expectEqual(eval(symbolicHistory, "integrate[Out[1],x]"),
        std::string{"(cos[x]+sin[x])exp[x]/2-(sin[x]-cos[x])exp[x]/2"},
        "integrate resolves an absolute Out snapshot inside its held integrand");
    tests.expectEqual(eval(symbolicHistory, "integrate[In[1],x]"),
        std::string{"(cos[x]+sin[x])exp[x]/2-(sin[x]-cos[x])exp[x]/2"},
        "integrate re-evaluates an absolute In snapshot inside its held integrand");

    kernel::KernelSession nestedSymbolicHistory;
    static_cast<void>(eval(nestedSymbolicHistory, "D[E^x cos[x],x]"));
    tests.expectEqual(eval(nestedSymbolicHistory, "integrate[2*Out[1],x]"),
        std::string{"2((cos[x]+sin[x])exp[x]/2-(sin[x]-cos[x])exp[x]/2)"},
        "held symbolic operators resolve nested output-history references without evaluating unrelated terms");

    static_cast<void>(evalError(history, "log[0]"));
    const std::size_t failedIndex = history.inputCount();
    tests.expect(history.inputHistory(failedIndex) != nullptr,
        "lowered input remains available after evaluation failure");
    tests.expect(history.outputHistory(failedIndex) == nullptr,
        "failed evaluation has no Out entry");
    tests.expect(evalError(history, "In[-1]").type() == error::CalcErrorType::Domain,
        "negative In can re-evaluate the previous lowered input even when its evaluation failed");
    tests.expectEqual(eval(history, "Out[-1]"), std::string{"100"},
        "negative Out skips failed inputs like percent history");
    tests.expect(evalError(history, "In[0]").type() == error::CalcErrorType::Type,
        "In[0] is invalid");
    tests.expect(evalError(history, "Out[0]").type() == error::CalcErrorType::Type,
        "Out[0] is invalid");
    tests.expect(evalError(history, "Out[-100]").type() == error::CalcErrorType::Evaluation,
        "unavailable negative Out entry reports an evaluation error");
    tests.expect(evalError(history, "Out[1000]").type() == error::CalcErrorType::Evaluation,
        "unavailable Out entry reports an evaluation error");

    const std::size_t beforeClear = history.inputCount();
    history.clearHistory();
    tests.expectEqual(eval(history, "7+8"), std::string{"15"},
        "absolute history remains aligned after clearHistory");
    tests.expect(history.outputHistory(beforeClear + 1) != nullptr,
        "Out absolute index follows the preserved input counter after clearHistory");

    kernel::KernelSession explainHistory;
    static_cast<void>(eval(explainHistory, "2+3"));
    tests.expectEqual(eval(explainHistory, "explain[Out[1]]"),
        std::string{"{{\"Kind\", \"Number\"}, {\"Domain\", \"Integer\"}, {\"Exactness\", \"Exact\"}, {\"Zero\", False}, {\"Sign\", \"Positive\"}, {\"BitLength\", 3}}"},
        "explain inspects the already evaluated Out value");

    kernel::KernelSession warnings;
    static_cast<void>(warnings.evaluate("D[abs[x],x]"));
    tests.expect(hasWarning(warnings, "D::unevaluated"),
        "unresolved D emits a warning");
    static_cast<void>(warnings.evaluate("D[x^2,x]"));
    tests.expect(!hasWarning(warnings, "D::unevaluated"),
        "completed D does not emit an unevaluated warning");

    static_cast<void>(warnings.evaluate("solve[sin[x]==0,x]"));
    tests.expect(hasWarning(warnings, "solve::unresolved"),
        "unresolved solve emits a warning");

    const std::size_t genericWarningsBefore = warningCount(warnings, "N::unevaluated");
    static_cast<void>(warnings.evaluate("N[x,20]"));
    tests.expect(warningCount(warnings, "N::unevaluated") == genericWarningsBefore,
        "N keeps a free symbolic value without emitting a generic warning");

    tests.expectEqual(eval(warnings, "N[lambertw[2,1],20]"),
        std::string{"-2.4015851048680028842+10.776299516115070898I"},
        "N certifies an explicit non-principal Lambert W branch");
    tests.expect(!hasWarning(warnings, "N::unsupported"),
        "implemented complex Lambert W branches do not retain the old unsupported diagnostic");

    kernel::KernelSession lambertCut;
    tests.expectEqual(eval(lambertCut, "N[lambertw[-1,-1/2],30]"),
        std::string{"lambertw[-1, -1/2]"},
        "N leaves the exact W_-1 negative-real cut unevaluated without logarithmic contraction");
    tests.expect(hasWarning(lambertCut, "N::unsupported"),
        "W_-1 negative-real cut is classified as a backend branch-convention gap");

    kernel::KernelSession boundedPrecision;
    auto boundedLimits = boundedPrecision.evaluationLimits();
    boundedLimits.maxCertifiedRefinements = 5'000;
    boundedPrecision.setEvaluationLimits(boundedLimits);
    static_cast<void>(boundedPrecision.evaluate(
        "N[sqrt[-1+I*sin[N[Pi,5]]],20]"));
    tests.expect(hasWarning(boundedPrecision, "N::precision")
            && !hasWarning(boundedPrecision, "N::unsupported")
            && boundedPrecision.lastEvaluationUsage().certifiedRefinements < 5'000,
        "N stops persistent branch ambiguity locally instead of exhausting the global refinement budget");

    kernel::KernelSession refinableBranch;
    const auto refinableValue = refinableBranch.evaluate(
        "N[log[-1+I*sin[Pi+1/10^100]],20]");
    tests.expect(refinableValue.isComplexDecimalApproximation()
            && !hasWarning(refinableBranch, "N::precision"),
        "N bounded refinement still resolves an exact branch side that needs extra guard precision");

    kernel::KernelSession approximateParameterPole;
    static_cast<void>(approximateParameterPole.evaluate(
        "N[hypergeometric2F1[1,2,N[0,5],2],20]"));
    tests.expect(hasWarning(approximateParameterPole, "N::precision")
            && !hasWarning(approximateParameterPole, "N::unsupported"),
        "N keeps finite 2F1 denominator-pole ambiguity as PrecisionInsufficient");

    kernel::KernelSession approximate1F1Pole;
    static_cast<void>(approximate1F1Pole.evaluate(
        "N[hypergeometric1F1[1,N[0,5],2],20]"));
    tests.expect(hasWarning(approximate1F1Pole, "N::precision")
            && !hasWarning(approximate1F1Pole, "N::unsupported"),
        "N keeps finite 1F1 denominator-pole ambiguity as PrecisionInsufficient");

    kernel::KernelSession unsupportedConsumers;
    tests.expect(evalError(unsupportedConsumers, "round[lambertw[-1]]").type()
            == error::CalcErrorType::Type,
        "rounding a certified non-real Lambert W value reports TypeError");
    const auto lambertDerivative = unsupportedConsumers.evaluate(
        "diff[lambertw[x],x,-1,20]");
    tests.expect(lambertDerivative.isComplexDecimalApproximation(),
        "numeric differentiation consumes the certified complex Lambert W derivative");
    tests.expect(evalError(
            unsupportedConsumers,
            "nintegrate[lambertw[x],{x,-2,-1},20]").type()
            == error::CalcErrorType::Evaluation,
        "numeric integration reports a remaining Lambert W backend gap as EvaluationError");

    kernel::KernelSession nestedWarnings;
    static_cast<void>(nestedWarnings.evaluate("N[D[abs[x],x],20]"));
    tests.expect(warningCount(nestedWarnings, "D::unevaluated") == 1
            && warningCount(nestedWarnings, "N::unevaluated") == 0,
        "outer N does not duplicate a more specific warning emitted by the inner operation");
    static_cast<void>(nestedWarnings.evaluate("N[True,20]"));
    static_cast<void>(nestedWarnings.evaluate("N[Infinity,20]"));
    static_cast<void>(nestedWarnings.evaluate("N[ComplexInfinity,20]"));
    static_cast<void>(nestedWarnings.evaluate("N[Indeterminate,20]"));
    tests.expect(warningCount(nestedWarnings, "N::unevaluated") == 0,
        "N treats inert Boolean and exceptional numeric atoms as intentional exact values without warning");

    kernel::KernelSession closedUnsupported;
    static_cast<void>(closedUnsupported.evaluate("N[0*Infinity,20]"));
    tests.expect(warningCount(closedUnsupported, "N::unevaluated") == 0,
        "N preserves a canonical Indeterminate result without a generic warning");

    const auto exactApprox = warnings.evaluate("N[1/3,20]");
    tests.expect(exactApprox.isDecimalApproximation(),
        "N exact rational produces a decimal approximation");
    if (exactApprox.isDecimalApproximation()) {
        const auto& value = exactApprox.asDecimalApproximation();
        tests.expect(value.origin() == numeric::ApproximationOrigin::ExactValue
            && value.certifiedEnclosureIsPoint()
            && value.displayedValue() != value.certifiedLower()
            && value.requestedSignificantDigits() == 20,
            "exact approximation retains point enclosure and requested digits");
    }

    tests.expectEqual(eval(warnings, "explain[{{1,2},{3,4}}]"),
        std::string{"{{\"Kind\", \"Array\"}, {\"Domain\", \"Integer\"}, {\"Exactness\", \"Exact\"}, {\"ArrayRank\", 2}, {\"Dimensions\", {2, 2}}, {\"ElementCount\", 4}, {\"Rectangular\", True}, {\"Empty\", False}, {\"Matrix\", True}, {\"Square\", True}, {\"Order\", 2}}"},
        "explain reports zero-cost Array shape metadata");
    tests.expectEqual(eval(warnings, "explain[{{1,2,3},{4,5,6}},\"internal\"]"),
        std::string{"{{\"Kind\", \"Array\"}, {\"Domain\", \"Integer\"}, {\"Exactness\", \"Exact\"}, {\"ArrayRank\", 2}, {\"Dimensions\", {2, 3}}, {\"ElementCount\", 6}, {\"Rectangular\", True}, {\"Empty\", False}, {\"Matrix\", True}, {\"Square\", False}, {\"Representation\", \"ArrayExpr\"}, {\"Storage\", \"Integer\"}, {\"Contiguous\", True}, {\"StoredExpressions\", False}}"},
        "explain internal exposes representation metadata without scanning the Array");

    const auto piApprox = warnings.evaluate("N[Pi,20]");
    tests.expect(piApprox.isDecimalApproximation(),
        "N Pi produces a certified decimal approximation");
    if (piApprox.isDecimalApproximation()) {
        const auto& value = piApprox.asDecimalApproximation();
        tests.expect(value.origin() == numeric::ApproximationOrigin::CertifiedInterval
            && value.certifiedLower() <= value.certifiedUpper()
            && value.requestedSignificantDigits() == 20,
            "certified approximation retains its enclosure metadata");
    }

    tests.expectEqual(eval(warnings, "explain[Pi]"),
        std::string{"{{\"Kind\", \"Constant\"}, {\"Domain\", \"Real\"}, {\"Exactness\", \"Exact\"}, {\"Name\", \"Pi\"}, {\"Real\", True}, {\"Positive\", True}, {\"Irrational\", True}, {\"ArithmeticClass\", \"Transcendental\"}}"},
        "explain uses MathRegistry metadata for Pi");
    tests.expectEqual(eval(warnings, "explain[Infinity]"),
        std::string{"{{\"Kind\", \"Constant\"}, {\"Domain\", \"ExtendedReal\"}, {\"Exactness\", \"Exact\"}, {\"Name\", \"Infinity\"}, {\"Infinite\", True}, {\"Finite\", False}, {\"Sign\", \"Positive\"}}"},
        "explain uses predefined-symbol metadata for Infinity");
    tests.expectEqual(eval(warnings, "explain[ComplexInfinity]"),
        std::string{"{{\"Kind\", \"Constant\"}, {\"Domain\", \"ExtendedComplex\"}, {\"Exactness\", \"Exact\"}, {\"Name\", \"ComplexInfinity\"}, {\"Infinite\", True}, {\"Finite\", False}, {\"Direction\", \"Undetermined\"}}"},
        "explain distinguishes directionless infinity from positive Infinity");
    tests.expectEqual(eval(warnings, "explain[Indeterminate]"),
        std::string{"{{\"Kind\", \"Indeterminate\"}, {\"Domain\", \"Undefined\"}, {\"Exactness\", \"Indeterminate\"}, {\"Name\", \"Indeterminate\"}, {\"Numeric\", False}, {\"Defined\", False}}"},
        "explain reports Indeterminate as a protected nonnumeric exceptional value");

    tests.expectEqual(eval(warnings, "explain[sin]"),
        std::string{"{{\"Kind\", \"BuiltinFunction\"}, {\"Domain\", \"Function\"}, {\"Exactness\", \"Exact\"}, {\"Name\", \"sin\"}, {\"Arity\", 1}, {\"ArgumentEvaluation\", \"All\"}, {\"FunctionDomain\", \"ComplexToComplexRealPreserving\"}, {\"Parity\", \"Odd\"}, {\"Branch\", \"SingleValued\"}, {\"PeriodTurns\", 1}, {\"PrincipalInverse\", \"asin\"}, {\"RealGloballyInjective\", False}, {\"RealRange\", \"ClosedMinusOneToOne\"}}"},
        "explain exposes zero-cost BuiltinRegistry and MathRegistry metadata");
    tests.expectEqual(eval(warnings, "explain[table]"),
        std::string{"{{\"Kind\", \"BuiltinFunction\"}, {\"Domain\", \"Function\"}, {\"Exactness\", \"Exact\"}, {\"Name\", \"table\"}, {\"Arity\", 2}, {\"ArgumentEvaluation\", \"HoldFirstAndTableIteratorSpec\"}}"},
        "explain reports held-argument evaluation semantics for table");

    const std::string piExplain = eval(warnings, "explain[N[Pi,20]]");
    tests.expect(piExplain.find("{\"Exactness\", \"CertifiedApproximation\"}") != std::string::npos
        && piExplain.find("{\"CertifiedEnclosure\", {") != std::string::npos
        && piExplain.find("{\"InformationEnclosure\", {") != std::string::npos,
        "explain exposes certified and information enclosures separately");
    tests.expect(evalError(warnings, "explain[1,\"full\"]").type() == error::CalcErrorType::Domain,
        "explain rejects unknown modes instead of silently changing cost semantics");
}

} // namespace mmcal::tests
