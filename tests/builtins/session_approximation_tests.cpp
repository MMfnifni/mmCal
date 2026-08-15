// session・approximation
#include "session_approximation_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
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

[[nodiscard]] const evaluation::EvaluationDiagnostic* findDiagnostic(
    const kernel::KernelSession& session,
    std::string_view code) {
    for (const auto& diagnostic : session.diagnostics())
        if (diagnostic.code == code)
            return &diagnostic;
    return nullptr;
}

} // namespace

void runSessionApproximationTests(TestRunner& tests) {
    kernel::KernelSession history;
    tests.expectEqual(eval(history, "1+1"), std::string{"2"},
        "history setup first result");
    tests.expectEqual(eval(history, "2+2"), std::string{"4"},
        "history setup second result");
    tests.expectEqual(eval(history, "In[1]+In[2]"), std::string{"6"},
        "In expressions participate in normal evaluation");
    tests.expectEqual(eval(history, "In[3]"), std::string{"6"},
        "nested In references are recursively evaluated");
    tests.expectEqual(eval(history, "In[1]"), std::string{"2"},
        "In[1] re-evaluates 1+1");
    tests.expect(evalError(history, "In[6]").type() == error::CalcErrorType::Evaluation,
        "current In entry cannot self-reference");

    kernel::KernelSession definitions;
    tests.expectEqual(eval(definitions, "x:=2"), std::string{"2"},
        "variable definition succeeds");
    tests.expectEqual(eval(definitions, "f[t]:=t^2"), std::string{"f[t]:=t^2"},
        "function definition succeeds");
    const std::string defs = eval(definitions, "Defs[]");
    tests.expect(defs.find("x:=2") != std::string::npos
        && defs.find("f[t]:=t^2") != std::string::npos,
        "Defs lists variables and user functions");

    tests.expectEqual(eval(definitions, "x:=4"), std::string{"4"},
        "variable redefinition succeeds");
    const auto* redefined = findDiagnostic(definitions, "definition::redefined");
    tests.expect(redefined
        && redefined->severity == evaluation::DiagnosticSeverity::Info
        && redefined->previousExpression
        && formatting::formatExpr(*redefined->previousExpression) == "2",
        "variable redefinition reports the previous value");

    tests.expectEqual(eval(definitions, "UnDef[x,f]"), std::string{"2"},
        "UnDef removes multiple names and returns changed-name count");
    tests.expectEqual(eval(definitions, "Defs[]"), std::string{"{}"},
        "UnDef removes variable and function definitions");

    const error::CalcError undefineImaginary = evalError(definitions, "UnDef[I]");
    tests.expect(undefineImaginary.type() == error::CalcErrorType::Name
        && std::string{undefineImaginary.what()} == "Cannot undefine protected symbol: I",
        "UnDef rejects the protected imaginary-unit name before value lowering");

    static_cast<void>(eval(definitions, "y:=9"));
    static_cast<void>(eval(definitions, "Clear[]"));
    tests.expectEqual(definitions.nextInputNumber(), std::size_t{1},
        "Clear resets absolute input numbering");
    tests.expect(definitions.environment().size() == 0 && definitions.userFunctions().size() == 0,
        "Clear removes user definitions");

    kernel::KernelSession exit;
    static_cast<void>(eval(exit, "Exit[]"));
    tests.expect(exit.exitRequested(), "Exit requests frontend termination");

    kernel::KernelSession approximation;
    tests.expectEqual(eval(approximation, "precision[1/3]"), std::string{"Infinity"},
        "exact numbers have infinite precision");
    tests.expectEqual(eval(approximation, "accuracy[1/3]"), std::string{"Infinity"},
        "exact numbers have infinite accuracy");
    tests.expectEqual(eval(approximation, "accuracy[N[1/3,20]]"), std::string{"20"},
        "approximation reports guaranteed absolute decimal digits");
    tests.expectEqual(eval(approximation, "precision[N[1/3,20]]"), std::string{"19"},
        "approximation reports conservative relative decimal digits");
    tests.expectEqual(eval(approximation, "rationalize[N[1/3,20]]"), std::string{"1/3"},
        "rationalize recovers a rational from a certified point enclosure");
    tests.expectEqual(eval(approximation, "N[N[1/3,20],10]"),
        std::string{"0.3333333333"},
        "outer N can safely reduce the requested digits of an existing approximation");
    tests.expectEqual(eval(approximation, "N[N[1/3,10],20]"),
        std::string{"0.3333333333"},
        "outer N never invents precision beyond an existing approximation");
    tests.expectEqual(eval(approximation, "N[N[Pi,20],100]"),
        std::string{"3.14159265358979323846"},
        "outer N preserves the available guarantee of a certified approximation");
    tests.expectEqual(eval(approximation, "N[Pi,20]+1/3"),
        std::string{"3.4749259869231265718"},
        "certified approximations can add exact Rational operands");
    tests.expectEqual(eval(approximation,
        "N[226375608064910089/72057594037927936,1000]-"
        "N[905502432259640355/288230376151711744,1000]"),
        std::string{"0.0000000000000000034694469519536141888238489627838134765625"},
        "certified approximations participate directly in subtraction");
    tests.expectEqual(eval(approximation,
        "accuracy[N[Pi,20]*10000000000]"),
        std::string{"10"},
        "approximation arithmetic reduces reported accuracy when scale amplifies uncertainty");
    tests.expectEqual(eval(approximation,
        "accuracy[N[Pi,100]*10^50]"),
        std::string{"50"},
        "approximation arithmetic never recovers hidden guard digits beyond the input guarantee");
    tests.expectEqual(eval(approximation, "N[Pi+I,20]*N[E-I,20]"),
        std::string{"9.5397342226735670655-0.4233108251307480031I"},
        "complex certified approximations participate in arithmetic");
    tests.expect(evalError(approximation, "1/N[0,20]").type() == error::CalcErrorType::Domain,
        "division detects an approximation whose certified enclosure is exactly zero");

    tests.expectEqual(eval(approximation, "rationalize[N[1/3,20],0]"),
        std::string{"33333333333333333333/100000000000000000000"},
        "zero-tolerance rationalize preserves the displayed decimal exactly");
    tests.expectEqual(eval(approximation, "N[arg[-1],20]"),
        std::string{"3.14159265358979323846 Rad"},
        "N approximates the numeric value inside an explicit angle unit");
    tests.expectEqual(eval(approximation, "N[Phi,20]"),
        std::string{"1.61803398874989484820"},
        "N certifies Phi from its exact algebraic definition");
}

} // namespace mmcal::tests
