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
        "approximation reports scale-dependent guaranteed absolute decimal digits");
    tests.expectEqual(eval(approximation, "precision[N[1/3,20]]"), std::string{"19"},
        "approximation reports conservative relative decimal digits");
    tests.expectEqual(eval(approximation, "rationalize[N[1/3,20]]"), std::string{"1/3"},
        "rationalize recovers a rational consistent with the information enclosure");
    tests.expectEqual(eval(approximation, "rationalize[N[1001/999,2]]"), std::string{"1"},
        "rationalize does not recover hidden exact point information beyond declared approximation quality");
    tests.expectEqual(eval(approximation, "N[N[1/3,20],10]"),
        std::string{"0.3333333333"},
        "outer N can safely reduce the requested digits of an existing approximation");
    tests.expectEqual(eval(approximation, "N[N[1/3,10],20]"),
        std::string{"0.3333333333"},
        "outer N never invents precision beyond an existing approximation");
    tests.expectEqual(eval(approximation, "N[N[Pi,20],100]"),
        std::string{"3.1415926535897932385"},
        "outer N preserves the available guarantee of a certified approximation");
    tests.expectEqual(eval(approximation, "accuracy[N[N[Pi,20],100]]"), std::string{"19"},
        "outer N cannot narrow an existing information enclosure beyond its declared quality");
    tests.expectEqual(eval(approximation, "N[Pi,20]+1/3"),
        std::string{"3.474925986923126572"},
        "certified approximations can add exact Rational operands");
    tests.expectEqual(eval(approximation,
        "N[226375608064910089/72057594037927936,1000]-"
        "N[905502432259640355/288230376151711744,1000]"),
        std::string{"0.0000000000000000034694469519536141888238489627838134765625"},
        "certified approximations participate directly in subtraction");
    tests.expectEqual(eval(approximation,
        "accuracy[N[Pi,20]*10000000000]"),
        std::string{"8"},
        "approximation arithmetic reduces reported accuracy when scale amplifies uncertainty");
    tests.expectEqual(eval(approximation,
        "accuracy[N[Pi,100]*10^50]"),
        std::string{"48"},
        "approximation arithmetic never recovers hidden guard digits beyond the input guarantee");
    tests.expectEqual(eval(approximation,
        "accuracy[(N[226375608064910089/72057594037927936,1000]-N[905502432259640355/288230376151711744,1000])]"),
        std::string{"998"},
        "near cancellation preserves high absolute accuracy through the information enclosure");
    tests.expectEqual(eval(approximation,
        "precision[(N[226375608064910089/72057594037927936,1000]-N[905502432259640355/288230376151711744,1000])]"),
        std::string{"980"},
        "near cancellation loses relative precision while retaining absolute accuracy");
    tests.expectEqual(eval(approximation, "N[Pi+I,20]*N[E-I,20]"),
        std::string{"9.53973422267356707-0.423310825130748003I"},
        "complex certified approximations participate in arithmetic");
    tests.expect(evalError(approximation, "1/N[0,20]").type() == error::CalcErrorType::Domain,
        "division detects an approximation whose certified enclosure is exactly zero");

    tests.expectEqual(eval(approximation, "rationalize[N[1/3,20],0]"),
        std::string{"33333333333333333333/100000000000000000000"},
        "zero-tolerance rationalize preserves the displayed decimal exactly");
    tests.expectEqual(eval(approximation, "N[arg[-1],20]"),
        std::string{"3.1415926535897932385 Rad"},
        "N approximates the numeric value inside an explicit angle unit");
    tests.expectEqual(eval(approximation, "N[Phi,20]"),
        std::string{"1.6180339887498948482"},
        "N certifies Phi from its exact algebraic definition");

    tests.expectEqual(eval(approximation, "N[Pi*10^20,20]"),
        std::string{"314159265358979323850"},
        "N interprets precision as significant digits at large scale");
    tests.expectEqual(eval(approximation, "precision[N[Pi*10^20,20]]"), std::string{"19"},
        "significant-digit N preserves relative precision at large scale");
    tests.expectEqual(eval(approximation, "accuracy[N[Pi/10^20,20]]"), std::string{"39"},
        "significant-digit N increases absolute accuracy at small scale");
    tests.expectEqual(eval(approximation, "precision[N[Pi/10^100,50]]"), std::string{"49"},
        "significant-digit N does not collapse tiny nonzero values to zero");
    tests.expectEqual(eval(approximation, "accuracy[N[Pi/10^100,50]]"), std::string{"149"},
        "tiny significant-digit approximations retain their scale-dependent absolute accuracy");

    tests.expectEqual(eval(approximation, "sin[N[Pi,20]]"),
        std::string{"0.000000000000000000"},
        "approximation-valued transcendental calls can return a zero-centered certified result");
    tests.expectEqual(eval(approximation, "accuracy[sin[N[Pi,20]]]"), std::string{"18"},
        "zero-centered approximation preserves absolute accuracy");
    tests.expectEqual(eval(approximation, "precision[sin[N[Pi,20]]]"), std::string{"0"},
        "zero-centered approximation correctly reports no relative precision");
    tests.expectEqual(eval(approximation, "exp[N[1,20]]"),
        std::string{"2.718281828459045235"},
        "DecimalApproximation is a first-class input to certified transcendental functions");
    tests.expectEqual(eval(approximation, "log10[N[2,20]]"),
        std::string{"0.3010299956639811952"},
        "approximation input survives primitive rewrite before certified evaluation");
    tests.expectEqual(eval(approximation, "fract[N[Pi,20]]"),
        std::string{"0.141592653589793238"},
        "fractional-part rewrite preserves certified approximation arithmetic");
    tests.expectEqual(eval(approximation, "exp[N[1+I,20]]"),
        std::string{"1.46869393991588516+2.28735528717884239I"},
        "ComplexDecimalApproximation is a first-class input to certified transcendental functions");
    tests.expectEqual(eval(approximation, "sqrt[N[-2,20]]"),
        std::string{"1.414213562373095049I"},
        "approximate real input may promote to a certified complex function result");
    tests.expectEqual(eval(approximation, "log[N[-2,20]]"),
        std::string{"0.6931471805599453094+3.141592653589793238I"},
        "principal-branch certified functions accept approximate negative real input");

    tests.expectEqual(eval(approximation, "Infinity-Infinity"), std::string{"Infinity-Infinity"},
        "Infinity does not use finite-symbol self-cancellation");
    tests.expectEqual(eval(approximation, "0*Infinity"), std::string{"0*Infinity"},
        "zero times Infinity remains conservatively unevaluated");
    tests.expectEqual(eval(approximation, "Infinity/Infinity"), std::string{"Infinity/Infinity"},
        "indeterminate Infinity division remains conservatively unevaluated");

    tests.expectEqual(eval(approximation, "N[Pi,20]>3"), std::string{"True"},
        "ordered comparison can use an InformationEnclosure when the result is provable");
    tests.expectEqual(eval(approximation, "N[Pi,2]>3.1"), std::string{"3.1>31/10"},
        "ordered comparison remains unresolved when InformationEnclosures overlap the boundary");
    tests.expectEqual(eval(approximation, "N[Pi,20]==3"), std::string{"False"},
        "approximate equality proves inequality only from disjoint InformationEnclosures");
    tests.expectEqual(eval(approximation, "N[Pi,2]==3.1"), std::string{"3.1==31/10"},
        "approximate equality remains unresolved when InformationEnclosures overlap");
    tests.expectEqual(eval(approximation, "N[Pi,2]!=3.1"), std::string{"3.1!=31/10"},
        "approximate inequality remains unresolved when InformationEnclosures overlap");
    tests.expectEqual(eval(approximation, "min[N[Pi,20],3]"), std::string{"3"},
        "min accepts certified approximations when InformationEnclosures order the inputs");
    tests.expectEqual(eval(approximation, "max[N[Pi,20],3]"),
        std::string{"3.1415926535897932385"},
        "max preserves the selected approximation without rerounding it");
}

} // namespace mmcal::tests
