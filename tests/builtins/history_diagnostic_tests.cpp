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

    static_cast<void>(evalError(history, "1/0"));
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

    static_cast<void>(warnings.evaluate("N[x,20]"));
    tests.expect(hasWarning(warnings, "N::unevaluated"),
        "unsupported N emits a warning");

    const auto exactApprox = warnings.evaluate("N[1/3,20]");
    tests.expect(exactApprox.isDecimalApproximation(),
        "N exact rational produces a decimal approximation");
    if (exactApprox.isDecimalApproximation()) {
        const auto& value = exactApprox.asDecimalApproximation();
        tests.expect(value.origin() == numeric::ApproximationOrigin::ExactValue
            && value.enclosureIsPoint()
            && value.displayedValue() != value.certifiedLower()
            && value.requestedFractionalDigits() == 20,
            "exact approximation retains point enclosure and requested digits");
    }

    const auto piApprox = warnings.evaluate("N[Pi,20]");
    tests.expect(piApprox.isDecimalApproximation(),
        "N Pi produces a certified decimal approximation");
    if (piApprox.isDecimalApproximation()) {
        const auto& value = piApprox.asDecimalApproximation();
        tests.expect(value.origin() == numeric::ApproximationOrigin::CertifiedInterval
            && value.certifiedLower() <= value.certifiedUpper()
            && value.requestedFractionalDigits() == 20,
            "certified approximation retains its enclosure metadata");
    }
}

} // namespace mmcal::tests
