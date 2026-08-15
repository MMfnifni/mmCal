// 定義・履歴・角度・診断を持つセッションの回帰テスト
#include "kernel_session_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

[[nodiscard]] std::string evaluateAndFormat(
    kernel::KernelSession& session,
    std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

[[nodiscard]] error::CalcError evaluateError(
    kernel::KernelSession& session,
    std::string_view source) {
    try {
        static_cast<void>(session.evaluate(source));
    }
    catch (const error::CalcError& exception) {
        return exception;
    }

    throw std::logic_error("Expected CalcError was not thrown");
}

} // namespace

void runKernelSessionTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(evaluateAndFormat(session, "1 + 2 * 3"), std::string{"7"},
        "KernelSession: evaluates a complete input");
    tests.expectEqual(session.historySize(), std::size_t{1},
        "KernelSession: stores successful result history");
    tests.expectEqual(session.inputCount(), std::size_t{1},
        "KernelSession: counts successful inputs");
    tests.expectEqual(session.nextInputNumber(), std::size_t{2},
        "KernelSession: exposes the next input number");

    tests.expectEqual(evaluateAndFormat(session, "0.1"), std::string{"1/10"},
        "KernelSession: decimal literals are exact rationals");
    tests.expectEqual(evaluateAndFormat(session, "0.1 + 0.2 == 0.3"), std::string{"True"},
        "KernelSession: decimal arithmetic remains exact");

    tests.expectEqual(evaluateAndFormat(session, "Pi"), std::string{"Pi"},
        "KernelSession: preserves Pi as an exact symbolic constant");
    tests.expectEqual(evaluateAndFormat(session, "E"), std::string{"E"},
        "KernelSession: preserves E as an exact symbolic constant");
    tests.expectEqual(evaluateAndFormat(session, "I ^ 2"), std::string{"-1"},
        "KernelSession: lowers I to the exact imaginary-unit value");

    kernel::KernelSession angleSession;
    tests.expectEqual(evaluateAndFormat(angleSession, "angleMode[]"), std::string{"Rad"},
        "KernelSession: angleMode reports the default Radian session mode");
    tests.expectEqual(evaluateAndFormat(angleSession, "angleMode[Deg]"), std::string{"Deg"},
        "KernelSession: angleMode switches the session to Degree");
    tests.expectEqual(evaluateAndFormat(angleSession, "sin[30]"), std::string{"1/2"},
        "KernelSession: angleMode changes implicit trigonometric angle semantics");
    tests.expectEqual(evaluateAndFormat(angleSession, "angleMode[Grad]"), std::string{"Grad"},
        "KernelSession: angleMode accepts Gradian");
    tests.expectEqual(evaluateAndFormat(angleSession, "sin[100]"), std::string{"1"},
        "KernelSession: Gradian mode is used by implicit trigonometric arguments");
    tests.expect(evaluateError(angleSession, "angleMode[Pi]").type() == error::CalcErrorType::Domain,
        "KernelSession: angleMode rejects non-angle-mode symbols");
    tests.expect(evaluateError(angleSession, "Rad:=2").type() == error::CalcErrorType::Syntax,
        "KernelSession: angle-mode enumeration symbols are protected");

    tests.expectEqual(evaluateAndFormat(session, "Tau"), std::string{"Tau"},
        "KernelSession: removed Tau is no longer predefined and currently passes as a free symbol");
    tests.expectEqual(evaluateAndFormat(session, "NA"), std::string{"NA"},
        "KernelSession: removed NA is no longer predefined and currently passes as a free symbol");
    tests.expectEqual(evaluateAndFormat(session, "ESP"), std::string{"ESP"},
        "KernelSession: removed ESP is no longer predefined and currently passes as a free symbol");
    tests.expectEqual(evaluateAndFormat(session, "x + 1"), std::string{"1+x"},
        "KernelSession: an unbound symbol is temporarily preserved for symbolic work");
    tests.expectEqual(evaluateAndFormat(session, "simplify[sin[1]^2 + cos[1]^2]"), std::string{"1"},
        "KernelSession: Simplify knows the Pythagorean trigonometric identity");
    tests.expectEqual(evaluateAndFormat(session, "simplify[sin[x]^2 + cos[x]^2]"), std::string{"1"},
        "KernelSession: Pythagorean identity works for a free symbolic angle");
    kernel::KernelSession nestedPowerSession;
    tests.expectEqual(evaluateAndFormat(
        nestedPowerSession, "simplify[(((((x+5))^3)^4)^3)^4]"),
        std::string{"(5+x)^144"},
        "KernelSession: simplify flattens nested positive exact integer powers");
    tests.expectEqual(evaluateAndFormat(
        nestedPowerSession,
        "expand[(((((x+5))^3)^4)^3)^4]==expand[(x+5)^144]"),
        std::string{"True"},
        "KernelSession: expand sees the canonical flattened nested Power");

    tests.expectEqual(evaluateAndFormat(
        nestedPowerSession,
        "simplify[((((7*((5)^3-2))+(2+-2))*((((1-x))^4)^4)^4))^3]"),
        std::string{"638277381(1-x)^192"},
        "KernelSession: positive integer Power extracts and evaluates an exact product coefficient");
    tests.expectEqual(evaluateAndFormat(
        nestedPowerSession,
        "simplify[(2*x)^3]"),
        std::string{"8x^3"},
        "KernelSession: exact numeric product coefficient is powered without branch assumptions");


    kernel::KernelSession elementarySession;
    tests.expectEqual(evaluateAndFormat(elementarySession, "abs[-3]"), std::string{"3"},
        "KernelSession: abs evaluates an exact negative real");
    tests.expectEqual(evaluateAndFormat(elementarySession, "abs[3 + 4I]"), std::string{"5"},
        "KernelSession: abs computes exact complex magnitude");
    tests.expectEqual(evaluateAndFormat(elementarySession, "abs[1 + I]"), std::string{"sqrt[2]"},
        "KernelSession: abs preserves an irrational exact magnitude");
    tests.expectEqual(evaluateAndFormat(elementarySession, "sign[0]"), std::string{"0"},
        "KernelSession: sign of zero is exact zero");
    tests.expectEqual(evaluateAndFormat(elementarySession, "sign[3 + 4I]"),
        std::string{"3/5+4/5I"},
        "KernelSession: complex sign is z divided by its magnitude");
    tests.expectEqual(evaluateAndFormat(elementarySession, "re[3 + 4I]"), std::string{"3"},
        "KernelSession: re extracts exact real part");
    tests.expectEqual(evaluateAndFormat(elementarySession, "im[3 + 4I]"), std::string{"4"},
        "KernelSession: im extracts exact imaginary part");
    tests.expectEqual(evaluateAndFormat(elementarySession, "conj[3 + 4I]"), std::string{"3-4I"},
        "KernelSession: conj evaluates exact complex conjugation");
    tests.expectEqual(evaluateAndFormat(elementarySession, "abs[x]"), std::string{"abs[x]"},
        "KernelSession: abs of an unconstrained symbol remains symbolic");
    tests.expectEqual(evaluateAndFormat(elementarySession, "simplify[sqrt[x^2]]"),
        std::string{"sqrt[x^2]"},
        "KernelSession: principal sqrt does not invent a PlusMinus result for unknown x");
    tests.expectEqual(evaluateAndFormat(elementarySession,
        "simplify[sqrt[x^2], element[x, Real]]"), std::string{"abs[x]"},
        "KernelSession: real assumption gives sqrt[x^2] = abs[x]");
    tests.expectEqual(evaluateAndFormat(elementarySession,
        "simplify[sqrt[x^2], x >= 0]"), std::string{"x"},
        "KernelSession: nonnegative assumption removes abs from principal sqrt");
    tests.expectEqual(evaluateAndFormat(elementarySession,
        "simplify[sqrt[x^2], x <= 0]"), std::string{"-x"},
        "KernelSession: nonpositive assumption chooses the nonnegative principal root");
    tests.expectEqual(evaluateAndFormat(elementarySession,
        "simplify[abs[x], x >= 0]"), std::string{"x"},
        "KernelSession: abs uses positive sign assumptions");
    tests.expectEqual(evaluateAndFormat(elementarySession,
        "simplify[abs[x], {element[x, Real], x >= 0}]"), std::string{"x"},
        "KernelSession: simplify accepts an assumption array");
    tests.expectEqual(evaluateAndFormat(elementarySession,
        "simplify[conj[x], element[x, Real]]"), std::string{"x"},
        "KernelSession: conjugation uses real-domain assumptions");
    tests.expectEqual(evaluateAndFormat(elementarySession, "element[Pi, Real]"), std::string{"True"},
        "KernelSession: element proves a permanent real-domain fact");
    tests.expectEqual(evaluateAndFormat(elementarySession, "element[I, Real]"), std::string{"False"},
        "KernelSession: element rejects a provably non-real value");
    tests.expectEqual(evaluateAndFormat(elementarySession, "element[x, Real]"),
        std::string{"element[x, Real]"},
        "KernelSession: unresolved element predicate remains symbolic");
    tests.expectEqual(evaluateAndFormat(elementarySession, "N[abs[Pi + I], 20]"),
        std::string{"3.2969083094756151588"},
        "KernelSession: abs uses certified complex magnitude evaluation");
    tests.expectEqual(evaluateAndFormat(elementarySession, "N[sign[Pi + I], 20]"),
        std::string{"0.95289051398868735278+0.30331447105335286402I"},
        "KernelSession: sign uses certified complex normalization");
    const error::CalcError contradictoryAssumptions = evaluateError(
        elementarySession, "simplify[abs[x], {x < 0, x > 0}]");
    tests.expect(contradictoryAssumptions.type() == error::CalcErrorType::Domain,
        "KernelSession: contradictory sign assumptions are rejected");
    tests.expectEqual(evaluateAndFormat(session, "expand[(x + 1)^3]"), std::string{"x^3+3x^2+3x+1"},
        "KernelSession: Expand distributes a polynomial expression");
    tests.expectEqual(evaluateAndFormat(session, "factor[x^2 - 1]"), std::string{"(x-1)(x+1)"},
        "KernelSession: Factor handles a univariate rational quadratic");
    tests.expectEqual(evaluateAndFormat(session, "collect[(x + 1)^3, x]"), std::string{"x^3+3x^2+3x+1"},
        "KernelSession: Collect uses the shared polynomial representation");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 == 1, x]"), std::string{"{x==1, x==-1}"},
        "KernelSession: Solve returns an internal SolutionSet for a quadratic");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^3 == -8, x]"),
        std::string{"{x==-2, x==1-I sqrt[3], x==1+I sqrt[3]}"},
        "KernelSession: Solve distinguishes all cube roots from principal Power");

    tests.expectEqual(evaluateAndFormat(session, "sqrt[8]"), std::string{"2sqrt[2]"},
        "KernelSession: radical normalization extracts square factors");
    tests.expectEqual(evaluateAndFormat(session, "sqrt[32]"), std::string{"4sqrt[2]"},
        "KernelSession: radical normalization avoids opaque sqrt[32] output");
    tests.expectEqual(evaluateAndFormat(session, "sqrt[2/3]"), std::string{"sqrt[6]/3"},
        "KernelSession: rational radical normalization rationalizes the denominator exactly");
    tests.expectEqual(evaluateAndFormat(session, "sqrt[-8]"), std::string{"2I sqrt[2]"},
        "KernelSession: negative rational radicals normalize on the principal complex branch");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 - 8 == 0, x]"),
        std::string{"{x==2sqrt[2], x==-2sqrt[2]}"},
        "KernelSession: solver roots pass through radical normalization");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 - 32 == 0, x]"),
        std::string{"{x==4sqrt[2], x==-4sqrt[2]}"},
        "KernelSession: solver never leaves a reducible sqrt[32] root opaque");

    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 < 4, x]"),
        std::string{"{x in Real if x>-2&&x<2}"},
        "KernelSession: strict quadratic inequality returns an exact real interval branch");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 <= 4, x]"),
        std::string{"{x in Real if x>=-2&&x<=2}"},
        "KernelSession: non-strict quadratic inequality includes exact endpoints");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 > 4, x]"),
        std::string{"{x in Real if x<-2, x in Real if x>2}"},
        "KernelSession: quadratic inequality can return a union of real branches");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 - 8 < 0, x]"),
        std::string{"{x in Real if x>-2sqrt[2]&&x<2sqrt[2]}"},
        "KernelSession: inequality endpoints use canonical exact radicals");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^3 - x > 0, x]"),
        std::string{"{x in Real if x>-1&&x<0, x in Real if x>1}"},
        "KernelSession: higher-degree fully split polynomials use an exact sign chart");
    tests.expectEqual(evaluateAndFormat(session, "solve[(x^2-1)*(x^2-4) >= 0, x]"),
        std::string{"{x in Real if x<=-2, x in Real if x>=-1&&x<=1, x in Real if x>=2}"},
        "KernelSession: high-degree non-strict sign charts merge roots into adjacent intervals");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 < 4, x, Integer]"),
        std::string{"{x in Integer if x>-2&&x<2}"},
        "KernelSession: inequality solution regions can be restricted to an ordered subdomain");
    const error::CalcError complexInequality = evaluateError(session, "solve[x^2 < 4, x, Complex]");
    tests.expect(complexInequality.type() == error::CalcErrorType::Domain,
        "KernelSession: ordered inequalities reject an explicit Complex search domain");

    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 != 1, x]"),
        std::string{"{x in Complex if x!=1&&x!=-1}"},
        "KernelSession: not-equal polynomial relations return an exact complement branch");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^2 + 1 != 0, x, Real]"),
        std::string{"All"},
        "KernelSession: Real-domain knowledge removes exclusions that are provably non-real");
    tests.expectEqual(evaluateAndFormat(session, "solve[{x > 0, x < 2}, x]"),
        std::string{"{x in Real if x>0&&x<2}"},
        "KernelSession: one-variable relation arrays are treated as conjunctions");
    tests.expectEqual(evaluateAndFormat(session, "solve[{x^2 == 1, x > 0}, x]"),
        std::string{"{x==1}"},
        "KernelSession: equation candidates are filtered by inequality constraints exactly");
    tests.expectEqual(evaluateAndFormat(session, "solve[(x-1)*(x^2-2) > 0, x]"),
        std::string{"{x in Real if x>-sqrt[2]&&x<1, x in Real if x>sqrt[2]}"},
        "KernelSession: higher-degree sign charts retain an irrational quadratic residual");
    tests.expectEqual(evaluateAndFormat(session, "solve[(x-1)/(x+1) > 0, x]"),
        std::string{"{x in Real if x<-1, x in Real if x>1}"},
        "KernelSession: rational-function inequalities preserve poles in the sign chart");
    tests.expectEqual(evaluateAndFormat(session, "solve[(x-1)/(x+1) >= 0, x]"),
        std::string{"{x in Real if x<-1, x in Real if x>=1}"},
        "KernelSession: rational-function zeros may be closed while poles remain excluded");
    tests.expectEqual(evaluateAndFormat(session, "solve[1/(x+1) < 2, x]"),
        std::string{"{x in Real if x<-1, x in Real if x>-1/2}"},
        "KernelSession: rational-function conversion combines a nonzero right-hand side safely");
    tests.expectEqual(evaluateAndFormat(session, "solve[(x-1)/(x+1) == 0, x]"),
        std::string{"{x==1}"},
        "KernelSession: rational equations solve the numerator while preserving denominator definedness");
    tests.expectEqual(evaluateAndFormat(session, "solve[(x-1)/(x-1) == 0, x]"),
        std::string{"{}"},
        "KernelSession: rational equation cancellation does not resurrect an undefined root");

    tests.expectEqual(evaluateAndFormat(session, "expand[(x + y)^2]"),
        std::string{"x^2+2x y+y^2"},
        "KernelSession: polynomial output uses mathematical implicit multiplication");
    tests.expectEqual(evaluateAndFormat(session, "factor[x^2 - y^2]"),
        std::string{"(x-y)(x+y)"},
        "KernelSession: factor handles a multivariate difference of squares");
    tests.expectEqual(evaluateAndFormat(session, "factor[x^3 - y^3]"),
        std::string{"(x-y)(x^2+x y+y^2)"},
        "KernelSession: factor handles a multivariate difference of cubes");
    tests.expectEqual(evaluateAndFormat(session, "factor[x^6 - 1]"),
        std::string{"(x-1)(x^2+x+1)(x+1)(x^2-x+1)"},
        "KernelSession: factor recursively splits higher-degree univariate factors");
    tests.expectEqual(evaluateAndFormat(session, "factor[2*x^3 - 3*x^2 - 8*x + 12]"),
        std::string{"(x-2)(x+2)(2x-3)"},
        "KernelSession: factor uses primitive rational-root linear factors without recursion cycles");
    tests.expectEqual(evaluateAndFormat(session, "collect[x*y + x*z + y, x]"),
        std::string{"(y+z)x+y"},
        "KernelSession: collect permits polynomial coefficients containing other variables");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^4 == 1, x]"),
        std::string{"{x==1, x==-1, x==I, x==-I}"},
        "KernelSession: solve handles higher-degree binomials exactly");
    tests.expectEqual(evaluateAndFormat(session, "solve[x^5 - x == 0, x]"),
        std::string{"{x==0, x==1, x==-1, x==I, x==-I}"},
        "KernelSession: solve deflates exact rational roots before solving the residual factor");
    tests.expectEqual(evaluateAndFormat(session, "solve[{x + y == 3, x - y == 1}, {x, y}]"),
        std::string{"{{x==2, y==1}}"},
        "KernelSession: solve handles an exact multivariable linear system");

    tests.expectEqual(
        evaluateAndFormat(session, "collect[(sin[y]+1)*x^2 + (log[y]+E)*x + Phi, x]"),
        std::string{"(1+sin[y])x^2+(E+log[y])x+Phi"},
        "KernelSession: collect treats arbitrary x-free exact expressions as coefficients");
    tests.expectEqual(
        evaluateAndFormat(session, "collect[sin[z]*x*y + log[z]*x + E, {x,y}]"),
        std::string{"(sin[z]y+log[z])x+E"},
        "KernelSession: multivariable collect recursively preserves expression-valued coefficients");
    tests.expectEqual(
        evaluateAndFormat(session, "factor[sin[y]*x^2 + 2*sin[y]*x + sin[y]]"),
        std::string{"(1+x)^2sin[y]"},
        "KernelSession: factor extracts a structural function coefficient before polynomial factoring");
    tests.expectEqual(
        evaluateAndFormat(session, "factor[(a+b)*x^2 + 2*(a+b)*x + (a+b)]"),
        std::string{"(a+b)(1+x)^2"},
        "KernelSession: factor recovers an expression coefficient even after evaluation distributed it");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[a*x + b == 0, x]"),
        std::string{"cases[{x==-b/a} if a!=0; All if a==0&&b==0; {} if a==0&&b!=0]"},
        "KernelSession: symbolic linear solve preserves all degenerate coefficient cases");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[x^2 + b*x + c == 0, x]"),
        std::string{"cases[{x==(-b+sqrt[b^2-4c])/2, x==(-b-sqrt[b^2-4c])/2} if b^2-4c!=0; {x==-b/2 (multiplicity 2)} if b^2-4c==0]"},
        "KernelSession: symbolic quadratic solve distinguishes distinct and repeated roots exactly");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[a*x^2 + b*x + c == 0, x]"),
        std::string{"cases[{x==(-b+sqrt[b^2-4a c])/(2a), x==(-b-sqrt[b^2-4a c])/(2a)} if a!=0&&b^2-4a c!=0; {x==-b/(2a) (multiplicity 2)} if a!=0&&b^2-4a c==0; {x==-c/b} if a==0&&b!=0; All if a==0&&b==0&&c==0; {} if a==0&&b==0&&c!=0]"},
        "KernelSession: symbolic quadratic solve also preserves linear and constant degeneracies");

    tests.expectEqual(
        evaluateAndFormat(session, "solve[sin[y]*x + 1 == 0, x]"),
        std::string{"cases[{x==-1/sin[y]} if sin[y]!=0; {} if sin[y]==0]"},
        "KernelSession: solver accepts parameter coefficients from entire functions and keeps zero cases");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[tan[y]*x + 1 == 0, x]"),
        std::string{"cases[{x==-1/tan[y]} if tan[y]!=0; {} if tan[y]==0] if cos[y]!=0"},
        "KernelSession: solver preserves tangent pole conditions for symbolic coefficients");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[sec[y]*x + 1 == 0, x]"),
        std::string{"{x==-1/sec[y]} if cos[y]!=0"},
        "KernelSession: solver preserves secant pole conditions for symbolic coefficients");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[coth[y]*x + 1 == 0, x]"),
        std::string{"cases[{x==-1/coth[y]} if coth[y]!=0; {} if coth[y]==0] if sinh[y]!=0"},
        "KernelSession: solver preserves hyperbolic cotangent pole conditions");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[atanh[y]*x + 1 == 0, x]"),
        std::string{"cases[{x==-1/atanh[y]} if atanh[y]!=0; {} if atanh[y]==0] if 1-y^2!=0"},
        "KernelSession: solver preserves inverse hyperbolic tangent singularities");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[(a/b)*x + 1 == 0, x]"),
        std::string{"cases[{x==-1/(a/b)} if a/b!=0; {} if a/b==0] if b!=0"},
        "KernelSession: solver carries a symbolic denominator domain condition into the solution set");

    tests.expectEqual(
        evaluateAndFormat(session, "solve[x^2 + 1 == 0, x, Real]"),
        std::string{"{}"},
        "KernelSession: solve can restrict the ambient solution domain to Real");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[x^2 + 1 == 0, x, Complex]"),
        std::string{"{x==I, x==-I}"},
        "KernelSession: explicit Complex domain matches the default complex solve semantics");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[x*(x-1) == 0, x, x != 0]"),
        std::string{"{x==1}"},
        "KernelSession: solve filters exact roots using a nonzero constraint");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[x^2 == 1, x, {Real, x > 0}]"),
        std::string{"{x==1}"},
        "KernelSession: solve combines a domain and a relation constraint");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[0 == 0, x, x != 0]"),
        std::string{"All if x!=0"},
        "KernelSession: universal solutions preserve a nonzero restriction");

    tests.expectEqual(
        evaluateAndFormat(session, "solve[{x + y == 3}, {x, y}]"),
        std::string{"{x==3-y where y in Complex}"},
        "KernelSession: underdetermined rational linear systems return a parametric solution");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[{x + y == 3}, {x, y}, Real]"),
        std::string{"{x==3-y where y in Real}"},
        "KernelSession: parametric free variables inherit the requested Real domain");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[{x + y == 3}, {x, y}, y != 0]"),
        std::string{"{x==3-y where y in Complex if y!=0}"},
        "KernelSession: constraints on a free parameter remain attached to the parametric branch");

    tests.expectEqual(
        evaluateAndFormat(session, "solve[{a*x + y == 1, x + y == 2}, {x, y}]"),
        std::string{"cases[{{x==-1/(a-1), y==(2a-1)/(a-1)}} if a!=1; Unresolved if a==1]"},
        "KernelSession: symbolic-coefficient square linear systems use an exact determinant condition");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[{a*x + y == 1, x + y == 2}, {x, y}, a != 1]"),
        std::string{"{{x==-1/(a-1), y==(2a-1)/(a-1)} if a!=1}"},
        "KernelSession: a user constraint can select the nonsingular symbolic linear-system branch");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[log[y]*x + 1 == 0, x]"),
        std::string{"cases[{x==-1/log[y]} if log[y]!=0; {} if log[y]==0] if y!=0"},
        "KernelSession: coefficient domain conditions include the principal logarithm domain");
    tests.expectEqual(
        evaluateAndFormat(session, "solve[log[b,y]*x + 1 == 0, x]"),
        std::string{"cases[{x==-1/log[b, y]} if log[b, y]!=0; {} if log[b, y]==0] if b!=0&&b!=1&&y!=0"},
        "KernelSession: arbitrary-base logarithm propagates base/value definedness into Solver");

    tests.expectEqual(evaluateAndFormat(session, "solve[x == 1, x]"),
        std::string{"{x==1}"},
        "KernelSession: canonical solve uses the solver implementation");

    const error::CalcError protectedConstant = evaluateError(session, "Pi := 3");
    tests.expect(protectedConstant.type() == error::CalcErrorType::Syntax,
        "KernelSession: predefined symbolic constants are protected from assignment");
    const error::CalcError protectedBoolean = evaluateError(session, "True := 0");
    tests.expect(protectedBoolean.type() == error::CalcErrorType::Syntax,
        "KernelSession: predefined literal names are protected from assignment");
    const error::CalcError protectedBuiltin = evaluateError(session, "sqrt := 3");
    tests.expect(protectedBuiltin.type() == error::CalcErrorType::Syntax,
        "KernelSession: source-callable builtin names are protected from assignment");
    const error::CalcError parenthesizedBuiltinCall = evaluateError(session, "sqrt(4)");
    tests.expect(parenthesizedBuiltinCall.type() == error::CalcErrorType::Syntax,
        "KernelSession: built-in function calls require square brackets");

    tests.expectEqual(evaluateAndFormat(session, "x := 3"), std::string{"3"},
        "KernelSession: evaluates assignment");
    tests.expectEqual(evaluateAndFormat(session, "x + 2"), std::string{"5"},
        "KernelSession: preserves definitions between inputs");
    tests.expectEqual(evaluateAndFormat(session, "x(4)"), std::string{"12"},
        "KernelSession: remembers variables for parenthesized multiplication");

    kernel::KernelSession historySession;
    static_cast<void>(historySession.evaluate("10"));
    static_cast<void>(historySession.evaluate("20"));
    tests.expectEqual(evaluateAndFormat(historySession, "%%"), std::string{"10"},
        "KernelSession: resolves multi-depth history");
    tests.expectEqual(evaluateAndFormat(historySession, "% + 5"), std::string{"15"},
        "KernelSession: uses history inside expressions");
    tests.expectEqual(evaluateAndFormat(historySession, "Out[-2]"), std::string{"10"},
        "KernelSession: resolves negative output history");

    kernel::KernelSession relativeInputSession;
    static_cast<void>(relativeInputSession.evaluate("7*8"));
    tests.expectEqual(evaluateAndFormat(relativeInputSession, "In[-1]"), std::string{"56"},
        "KernelSession: negative input history re-evaluates previous input");

    kernel::KernelSession atHistorySession;
    static_cast<void>(atHistorySession.evaluate("7*8"));
    static_cast<void>(atHistorySession.evaluate("3+4"));
    static_cast<void>(atHistorySession.evaluate("2^5"));
    tests.expectEqual(evaluateAndFormat(atHistorySession, "@"), std::string{"32"},
        "KernelSession: at shorthand resolves previous input");

    kernel::KernelSession repeatedAtHistorySession;
    static_cast<void>(repeatedAtHistorySession.evaluate("7*8"));
    static_cast<void>(repeatedAtHistorySession.evaluate("3+4"));
    static_cast<void>(repeatedAtHistorySession.evaluate("2^5"));
    tests.expectEqual(evaluateAndFormat(repeatedAtHistorySession, "@@"), std::string{"7"},
        "KernelSession: repeated at shorthand resolves older input");

    kernel::KernelSession deepAtHistorySession;
    static_cast<void>(deepAtHistorySession.evaluate("7*8"));
    static_cast<void>(deepAtHistorySession.evaluate("3+4"));
    static_cast<void>(deepAtHistorySession.evaluate("2^5"));
    tests.expectEqual(evaluateAndFormat(deepAtHistorySession, "@@@"), std::string{"56"},
        "KernelSession: repeated at shorthand depth matches In[-n]");

    kernel::KernelSession deepPercentHistorySession;
    static_cast<void>(deepPercentHistorySession.evaluate("7*8"));
    static_cast<void>(deepPercentHistorySession.evaluate("3+4"));
    static_cast<void>(deepPercentHistorySession.evaluate("2^5"));
    tests.expectEqual(evaluateAndFormat(deepPercentHistorySession, "%%%"), std::string{"56"},
        "KernelSession: repeated percent shorthand resolves older output");

    kernel::KernelSession emptyHistory;
    const error::CalcError historyError = evaluateError(emptyHistory, "%");
    tests.expect(historyError.type() == error::CalcErrorType::Evaluation,
        "KernelSession: missing history is an evaluation error");
    tests.expect(historyError.span().has_value(),
        "KernelSession: missing history keeps source position");
    tests.expectEqual(emptyHistory.historySize(), std::size_t{0},
        "KernelSession: failed input is not added to history");
    tests.expectEqual(emptyHistory.inputCount(), std::size_t{1},
        "KernelSession: counts failed inputs");


    kernel::KernelSession numericalSession;
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[1/2]"),
        std::string{"0.5"},
        "KernelSession: N preserves terminating decimal length");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[1/3]"),
        std::string{"0.3333333333333333"},
        "KernelSession: N defaults to sixteen significant digits");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[1/3, 20]"),
        std::string{"0.33333333333333333333"},
        "KernelSession: N accepts explicit significant digits");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[Pi]"),
        std::string{"3.141592653589793"},
        "KernelSession: N invokes the certified Pi provider");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[Pi, 100]"),
        std::string{"3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068"},
        "KernelSession: Pi precision is arbitrary and certified");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[E, 30]"),
        std::string{"2.71828182845904523536028747135"},
        "KernelSession: E is certified from Exp[1], not stored machine precision");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[exp[1/3], 30]"),
        std::string{"1.3956124250860895286281253196"},
        "KernelSession: real Exp uses arbitrary-precision certified evaluation");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[log[2], 30]"),
        std::string{"0.693147180559945309417232121458"},
        "KernelSession: positive real Log uses certified arbitrary precision");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[log[10,2], 30]"),
        std::string{"0.301029995663981195213738894724"},
        "KernelSession: arbitrary-base Log uses certified change-of-base evaluation");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[arg[2 + 3 I], 30]"),
        std::string{"0.982793723247329067985710611015"},
        "KernelSession: certified atan2 evaluates a general principal Arg");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[log[1 + I], 30]"),
        std::string{"0.346573590279972654708616060729+0.78539816339744830961566084582I"},
        "KernelSession: N evaluates the principal complex logarithm");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[exp[1 + I], 30]"),
        std::string{"1.46869393991588515713896759733+2.2873552871788423912081719067I"},
        "KernelSession: N evaluates general complex Exp without machine floating point");
    tests.expectEqual(evaluateAndFormat(numericalSession, "(-8)^(1/3)"),
        std::string{"(-8)^(1/3)"},
        "KernelSession: exact principal Power remains symbolic when no simple exact form is implemented");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[(-8)^(1/3), 30]"),
        std::string{"1+1.73205080756887729352744634151I"},
        "KernelSession: non-integer Power uses Exp[w principal Log[z]]");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sqrt[2]]"),
        std::string{"1.414213562373095"},
        "KernelSession: N uses certified arbitrary-precision real sqrt");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sqrt[2], 50]"),
        std::string{"1.4142135623730950488016887242096980785696718753769"},
        "KernelSession: sqrt precision is not tied to machine floating point");
    tests.expectEqual(evaluateAndFormat(numericalSession, "(-2)^0.5"),
        std::string{"I sqrt[2]"},
        "KernelSession: decimal one-half exponent uses principal sqrt exactly");
    tests.expectEqual(evaluateAndFormat(numericalSession, "-2^0.5"),
        std::string{"-sqrt[2]"},
        "KernelSession: unary minus remains outside power without parentheses");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[(-2)^0.5, 30]"),
        std::string{"1.41421356237309504880168872421I"},
        "KernelSession: N evaluates exact complex half power through ComplexInterval");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sqrt[Pi], 30]"),
        std::string{"1.77245385090551602729816748334"},
        "KernelSession: certified evaluator composes Pi and sqrt");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[Pi + sqrt[2], 30]"),
        std::string{"4.55580621596288828726433210749"},
        "KernelSession: N evaluates a composed exact expression instead of per-function dispatch");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sqrt[-2]"),
        std::string{"I sqrt[2]"},
        "KernelSession: negative non-square real sqrt promotes exactly to complex");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sqrt[-Pi]"),
        std::string{"I sqrt[Pi]"},
        "KernelSession: domain knowledge promotes negative symbolic real sqrt to complex");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sqrt[3 + 4 I]"),
        std::string{"2+I"},
        "KernelSession: exact rational complex square root uses the principal branch");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sqrt[-3 - 4 I]"),
        std::string{"1-2I"},
        "KernelSession: principal complex square root keeps non-negative real part");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sqrt[1 + I], 30]"),
        std::string{"1.09868411346780996603980119524+0.455089860562227341304357757822I"},
        "KernelSession: N certifies a general upper-half-plane complex square root");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sqrt[-1 - I], 30]"),
        std::string{"0.455089860562227341304357757822-1.09868411346780996603980119524I"},
        "KernelSession: N certifies principal sqrt below the branch cut");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[1]"),
        std::string{"sin[1]"},
        "KernelSession: sin remains exact and symbolic without N");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cos[1]"),
        std::string{"cos[1]"},
        "KernelSession: cos remains exact and symbolic without N");
    tests.expectEqual(evaluateAndFormat(numericalSession, "tan[1]"),
        std::string{"tan[1]"},
        "KernelSession: tan remains exact and symbolic without N");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[180 Deg]"),
        std::string{"0"},
        "KernelSession: explicit Degree special angle simplifies");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cos[180 Deg]"),
        std::string{"-1"},
        "KernelSession: exact degree special angle simplifies");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[30 Deg]"),
        std::string{"1/2"},
        "KernelSession: sin 30 degrees is exact one half");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cos[60 Deg]"),
        std::string{"1/2"},
        "KernelSession: cos 60 degrees is exact one half");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[45 Deg]"),
        std::string{"sqrt[2]/2"},
        "KernelSession: exact special angles may produce symbolic radicals");
    tests.expectEqual(evaluateAndFormat(numericalSession, "tan[45 Deg]"),
        std::string{"1"},
        "KernelSession: tan 45 degrees is exact one");
    tests.expectEqual(evaluateAndFormat(numericalSession, "tan[30 Deg]"),
        std::string{"sqrt[3]/3"},
        "KernelSession: tan 30 degrees is an exact radical");
    tests.expectEqual(evaluateAndFormat(numericalSession, "tan[135 Deg]"),
        std::string{"-1"},
        "KernelSession: tangent uses a half-turn period and correct quadrant sign");
    const error::CalcError tangentPole = evaluateError(numericalSession, "tan[90 Deg]");
    tests.expect(tangentPole.type() == error::CalcErrorType::Domain,
        "KernelSession: exact tangent pole is a domain error");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[Pi]"),
        std::string{"0"},
        "KernelSession: bare trig arguments use Radian by default");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[Pi Rad]"),
        std::string{"0"},
        "KernelSession: Pi Rad is recognized exactly without approximating Pi");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cos[Pi Rad]"),
        std::string{"-1"},
        "KernelSession: cosine understands exact Pi radians");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[(Pi Rad) / 6]"),
        std::string{"1/2"},
        "KernelSession: explicit angle arithmetic remains exact in turns");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cos[2 * (Pi Rad)]"),
        std::string{"1"},
        "KernelSession: explicit radian angle arithmetic recognizes full turns");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[100 Grad]"),
        std::string{"1"},
        "KernelSession: explicit Gradian is converted exactly to turns");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[90 deg]"),
        std::string{"1"},
        "KernelSession: legacy lowercase angle unit remains accepted");

    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[15 Deg]"),
        std::string{"(sqrt[6]-sqrt[2])/4"},
        "KernelSession: sin 15 degrees keeps an exact radical value");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cos[75 Deg]"),
        std::string{"(sqrt[6]-sqrt[2])/4"},
        "KernelSession: cos 75 degrees keeps the matching exact radical value");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[15 Deg] + cos[75 Deg]"),
        std::string{"(sqrt[6]-sqrt[2])/2"},
        "KernelSession: identical exact radical terms are collected algebraically");
    tests.expectEqual(evaluateAndFormat(numericalSession, "2 * sin[15 Deg]"),
        std::string{"(sqrt[6]-sqrt[2])/2"},
        "KernelSession: rational scaling reduces exact symbolic denominators");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sin[75 Deg] + cos[15 Deg]"),
        std::string{"(sqrt[2]+sqrt[6])/2"},
        "KernelSession: exact trig values combine without numerical approximation");

    // 既定はRadian。明示Deg/GradだけがPiを介した角度変換を行う。
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sin[30 Deg], 100]"),
        std::string{"0.5"},
        "KernelSession: N sees exact special-angle simplification before approximation");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sin[1], 100]"),
        std::string{"0.8414709848078965066525023216302989996225630607983710656727517099919104043912396689486397435430526959"},
        "KernelSession: bare sin argument is certified in default Radian mode");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[cos[1], 100]"),
        std::string{"0.5403023058681397174009366074429766037323104206179222276700972553811003947744717645179518560871830893"},
        "KernelSession: bare cos argument is certified in default Radian mode");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sin[1 Grad], 50]"),
        std::string{"0.015707317311820675753295353309906770086948450733779"},
        "KernelSession: explicit Grad is converted through certified Pi");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sin[1 Rad], 100]"),
        std::string{"0.8414709848078965066525023216302989996225630607983710656727517099919104043912396689486397435430526959"},
        "KernelSession: explicit Rad keeps certified arbitrary precision sin backend");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[cos[1 Rad], 100]"),
        std::string{"0.5403023058681397174009366074429766037323104206179222276700972553811003947744717645179518560871830893"},
        "KernelSession: explicit Rad keeps certified arbitrary precision cos backend");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[tan[1 Rad], 100]"),
        std::string{"1.557407724654902230506974807458360173087250772381520038383946605698861397151727289555099965202242984"},
        "KernelSession: explicit Rad uses certified arbitrary precision tan backend");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[tan[1], 100]"),
        std::string{"1.557407724654902230506974807458360173087250772381520038383946605698861397151727289555099965202242984"},
        "KernelSession: bare tan argument is certified in default Radian mode");


    // 逆数三角函数も明示角度単位のexact angle reductionを共有する。
    tests.expectEqual(evaluateAndFormat(numericalSession, "cot[45 Deg]"), std::string{"1"},
        "KernelSession: cot 45 degrees is exact one");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cot[90 Deg]"), std::string{"0"},
        "KernelSession: cot 90 degrees is zero instead of inheriting tan's pole");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sec[60 Deg]"), std::string{"2"},
        "KernelSession: sec 60 degrees is exact two");
    tests.expectEqual(evaluateAndFormat(numericalSession, "csc[30 Deg]"), std::string{"2"},
        "KernelSession: csc 30 degrees is exact two");
    tests.expect(evaluateError(numericalSession, "cot[0]").type() == error::CalcErrorType::Domain,
        "KernelSession: cot rejects a sine zero");
    tests.expect(evaluateError(numericalSession, "sec[90 Deg]").type() == error::CalcErrorType::Domain,
        "KernelSession: sec rejects a cosine zero");
    tests.expect(evaluateError(numericalSession, "csc[0]").type() == error::CalcErrorType::Domain,
        "KernelSession: csc rejects a sine zero");

    // 逆三角函数の公開値は現在のAngleSemanticsに従い、既定ではRadian。
    tests.expectEqual(evaluateAndFormat(numericalSession, "asin[1/2]"), std::string{"Pi/6"},
        "KernelSession: asin one half returns Pi/6 radians exactly");
    tests.expectEqual(evaluateAndFormat(numericalSession, "acos[1/2]"), std::string{"Pi/3"},
        "KernelSession: acos one half returns Pi/3 radians exactly");
    tests.expectEqual(evaluateAndFormat(numericalSession, "atan[1]"), std::string{"Pi/4"},
        "KernelSession: atan one returns Pi/4 radians exactly");
    tests.expectEqual(evaluateAndFormat(numericalSession, "asin[sqrt[2]/2]"), std::string{"Pi/4"},
        "KernelSession: asin recognizes an exact radical special value");
    tests.expectEqual(evaluateAndFormat(numericalSession, "acos[-1/2]"), std::string{"2Pi/3"},
        "KernelSession: acos uses its principal range from zero through 180 degrees");
    tests.expectEqual(evaluateAndFormat(numericalSession, "atan[-sqrt[3]]"), std::string{"-Pi/3"},
        "KernelSession: atan preserves principal odd symmetry for exact radicals");
    tests.expectEqual(evaluateAndFormat(numericalSession, "atan2[1,1]"), std::string{"Pi/4"},
        "KernelSession: atan2 first quadrant is exact");
    tests.expectEqual(evaluateAndFormat(numericalSession, "atan2[1,-1]"), std::string{"3Pi/4"},
        "KernelSession: atan2 distinguishes the second quadrant");
    tests.expectEqual(evaluateAndFormat(numericalSession, "atan2[-1,-1]"), std::string{"-3Pi/4"},
        "KernelSession: atan2 distinguishes the third quadrant");
    tests.expectEqual(evaluateAndFormat(numericalSession, "atan2[0,-1]"), std::string{"Pi"},
        "KernelSession: atan2 uses the Arg-compatible negative-axis endpoint");
    tests.expect(evaluateError(numericalSession, "atan2[0,0]").type() == error::CalcErrorType::Domain,
        "KernelSession: atan2 is undefined at the origin");
    tests.expect(evaluateError(numericalSession, "atan2[I,1]").type() == error::CalcErrorType::Domain,
        "KernelSession: atan2 rejects exact complex coordinates without waiting for N");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[asin[1/3], 20]"),
        std::string{"0.3398369094541219371"},
        "KernelSession: asin arbitrary precision output uses radians by default");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[atan2[2,3], 20]"),
        std::string{"0.58800260354756755125"},
        "KernelSession: atan2 arbitrary precision keeps quadrant-aware radian output");

    // 複素逆三角函数はprincipal Log/Sqrtの枝を使う。負実軸のArgは+Pi側。
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[asin[2], 20]"),
        std::string{"1.5707963267948966192-1.3169578969248167086I"},
        "KernelSession: complex asin follows the principal branch");
    tests.expect(evaluateError(numericalSession, "N[atan[I],20]").type() == error::CalcErrorType::Domain,
        "KernelSession: atan detects its logarithmic singularity at I");

    // 双曲線函数は角度単位を通さず、無次元入力をそのまま評価する。
    tests.expectEqual(evaluateAndFormat(numericalSession, "sinh[0]"), std::string{"0"},
        "KernelSession: sinh zero is exact zero");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cosh[0]"), std::string{"1"},
        "KernelSession: cosh zero is exact one");
    tests.expectEqual(evaluateAndFormat(numericalSession, "tanh[0]"), std::string{"0"},
        "KernelSession: tanh zero is exact zero");
    tests.expectEqual(evaluateAndFormat(numericalSession, "asinh[0]"), std::string{"0"},
        "KernelSession: asinh zero is exact zero");
    tests.expectEqual(evaluateAndFormat(numericalSession, "acosh[1]"), std::string{"0"},
        "KernelSession: acosh one is exact zero");
    tests.expectEqual(evaluateAndFormat(numericalSession, "atanh[0]"), std::string{"0"},
        "KernelSession: atanh zero is exact zero");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sech[0]"), std::string{"1"},
        "KernelSession: sech zero is exact one");
    tests.expect(evaluateError(numericalSession, "csch[0]").type() == error::CalcErrorType::Domain,
        "KernelSession: csch rejects a sinh zero");
    tests.expect(evaluateError(numericalSession, "coth[0]").type() == error::CalcErrorType::Domain,
        "KernelSession: coth rejects a sinh zero");
    tests.expect(evaluateError(numericalSession, "atanh[1]").type() == error::CalcErrorType::Domain,
        "KernelSession: atanh rejects its positive logarithmic pole");
    tests.expect(evaluateError(numericalSession, "atanh[-1]").type() == error::CalcErrorType::Domain,
        "KernelSession: atanh rejects its negative logarithmic pole");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cosh[I Pi]"), std::string{"-1"},
        "KernelSession: cosh of I Pi reduces exactly through cosine");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sinh[I Pi/2]"), std::string{"I"},
        "KernelSession: sinh of I Pi over two reduces exactly through sine");
    tests.expect(evaluateError(numericalSession, "tanh[I Pi/2]").type() == error::CalcErrorType::Domain,
        "KernelSession: tanh detects an exact complex cosh zero");

    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sinh[1], 20]"),
        std::string{"1.1752011936438014569"},
        "KernelSession: sinh uses certified arbitrary precision evaluation");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[atanh[1/2], 20]"),
        std::string{"0.5493061443340548457"},
        "KernelSession: atanh uses certified arbitrary precision evaluation");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[csch[1], 20]"),
        std::string{"0.85091812823932154513"},
        "KernelSession: csch uses certified reciprocal sinh evaluation");

    tests.expectEqual(evaluateAndFormat(numericalSession, "N[acosh[-2], 20]"),
        std::string{"1.3169578969248167086+3.1415926535897932385I"},
        "KernelSession: complex acosh uses the Arg-compatible upper cut value");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[atanh[2], 20]"),
        std::string{"0.5493061443340548457-1.5707963267948966192I"},
        "KernelSession: complex atanh uses the principal Log boundary convention");

    // parity metadata is consumed by the Simplifier instead of being duplicated per function.
    tests.expectEqual(evaluateAndFormat(numericalSession, "cot[-x]"), std::string{"-cot[x]"},
        "KernelSession: cot oddness is simplified symbolically");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sec[-x]"), std::string{"sec[x]"},
        "KernelSession: sec evenness is simplified symbolically");
    tests.expectEqual(evaluateAndFormat(numericalSession, "sinh[-x]"), std::string{"-sinh[x]"},
        "KernelSession: sinh oddness is simplified symbolically");
    tests.expectEqual(evaluateAndFormat(numericalSession, "cosh[-x]"), std::string{"cosh[x]"},
        "KernelSession: cosh evenness is simplified symbolically");
    tests.expectEqual(evaluateAndFormat(numericalSession, "N[sin[1 Rad], 37]"),
        std::string{"0.8414709848078965066525023216302989996"},
        "KernelSession: N precision is not hard-coded to 100 digits");

    const error::CalcError precisionError = evaluateError(numericalSession, "N[1/3, 0]");
    tests.expect(precisionError.type() == error::CalcErrorType::Type,
        "KernelSession: N rejects non-positive precision");
    tests.expect(precisionError.span().has_value()
        && precisionError.span()->begin.column == 8,
        "KernelSession: N precision error points to the second argument");

    kernel::KernelSession functionSession;
    tests.expectEqual(evaluateAndFormat(functionSession, "f[x]:=x+1"),
        std::string{"f[x]:=x+1"},
        "KernelSession: stores delayed user-function definition");
    tests.expectEqual(functionSession.userFunctions().size(), std::size_t{1},
        "KernelSession: exposes user-function definitions");
    tests.expectEqual(evaluateAndFormat(functionSession, "f[3]"), std::string{"4"},
        "KernelSession: evaluates user function with bracket call");
    const error::CalcError parenthesizedUserCall = evaluateError(functionSession, "f(4)");
    tests.expect(parenthesizedUserCall.type() == error::CalcErrorType::Syntax,
        "KernelSession: rejects parenthesized user-function calls");

    static_cast<void>(functionSession.evaluate("x := 100"));
    tests.expectEqual(evaluateAndFormat(functionSession, "f[2]"), std::string{"3"},
        "KernelSession: function parameter shadows global variable");
    tests.expectEqual(evaluateAndFormat(functionSession, "x"), std::string{"100"},
        "KernelSession: function call preserves shadowed global variable");

    static_cast<void>(functionSession.evaluate("a := 10"));
    static_cast<void>(functionSession.evaluate("g[t] := t + a"));
    tests.expectEqual(evaluateAndFormat(functionSession, "g[2]"), std::string{"12"},
        "KernelSession: function body can read global definitions");
    static_cast<void>(functionSession.evaluate("a := 20"));
    tests.expectEqual(evaluateAndFormat(functionSession, "g[2]"), std::string{"22"},
        "KernelSession: delayed function body observes current global values");

    static_cast<void>(functionSession.evaluate("f[x, y] := x + y"));
    tests.expectEqual(evaluateAndFormat(functionSession, "f[2, 3]"), std::string{"5"},
        "KernelSession: supports user-function overloads by arity");
    tests.expectEqual(evaluateAndFormat(functionSession, "f[3]"), std::string{"4"},
        "KernelSession: arity overload preserves existing definition");

    static_cast<void>(functionSession.evaluate("f[x] := x * 2"));
    tests.expectEqual(evaluateAndFormat(functionSession, "f[3]"), std::string{"6"},
        "KernelSession: same-arity definition replaces previous definition");

    const error::CalcError arityError = evaluateError(functionSession, "f[1, 2, 3]");
    tests.expect(arityError.type() == error::CalcErrorType::Type,
        "KernelSession: undefined user-function arity is a type error");
    tests.expect(arityError.span().has_value(),
        "KernelSession: user-function arity error keeps call position");

    const error::CalcError builtinRedefinition = evaluateError(
        functionSession,
        "sqrt[x] := x");
    tests.expect(builtinRedefinition.type() == error::CalcErrorType::Syntax,
        "KernelSession: builtin functions cannot be redefined");

    const error::CalcError builtinParameter = evaluateError(
        functionSession,
        "h[sqrt] := sqrt");
    tests.expect(builtinParameter.type() == error::CalcErrorType::Syntax,
        "KernelSession: builtin function names cannot be function parameters");

    const error::CalcError constantParameter = evaluateError(
        functionSession,
        "h[Pi] := Pi");
    tests.expect(constantParameter.type() == error::CalcErrorType::Syntax,
        "KernelSession: constants cannot be function parameters");

    tests.expectEqual(
        evaluateAndFormat(functionSession, "if[True, 1, 1 / 0]"),
        std::string{"1"},
        "KernelSession: If does not evaluate the unselected false branch");
    tests.expectEqual(
        evaluateAndFormat(functionSession, "if[False, 1 / 0, 2]"),
        std::string{"2"},
        "KernelSession: If does not evaluate the unselected true branch");

    const error::CalcError ifConditionError = evaluateError(functionSession, "if[1, 2, 3]");
    tests.expect(ifConditionError.type() == error::CalcErrorType::Type,
        "KernelSession: If requires a Boolean condition");
    tests.expect(ifConditionError.span().has_value()
        && ifConditionError.span()->begin.column == 4,
        "KernelSession: If condition error points to the condition");

    static_cast<void>(functionSession.evaluate(
        "recursiveFact[n] := if[n == 0, 1, n * recursiveFact[n - 1]]"));
    tests.expectEqual(evaluateAndFormat(functionSession, "recursiveFact[0]"), std::string{"1"},
        "KernelSession: recursive function reaches its base case");
    tests.expectEqual(evaluateAndFormat(functionSession, "recursiveFact[20]"),
        std::string{"2432902008176640000"},
        "KernelSession: recursive function evaluates through local scopes");

    kernel::KernelSession provenanceSession;
    static_cast<void>(provenanceSession.evaluate("inner[x] := 1 / (x - x)"));
    static_cast<void>(provenanceSession.evaluate("outer[y] := inner[y]"));
    const error::CalcError provenanceError = evaluateError(provenanceSession, "outer[4]");
    tests.expect(provenanceError.document()
        && provenanceError.document()->inputNumber() == 1,
        "KernelSession: function-body error keeps definition input document");
    tests.expectEqual(provenanceError.sourceLabel(), std::string{"Defined at"},
        "KernelSession: function-body error labels its definition source");
    tests.expectEqual(provenanceError.trace().size(), std::size_t{2},
        "KernelSession: nested function error keeps call chain");
    tests.expect(provenanceError.trace().size() == 2
        && provenanceError.trace()[0].source.document->inputNumber() == 2
        && provenanceError.trace()[1].source.document->inputNumber() == 3,
        "KernelSession: call chain preserves inner-to-outer input order");

    const std::string provenanceMessage = error::errorMessage(provenanceError);
    tests.expect(provenanceMessage.find("Defined at In [1]") != std::string::npos,
        "KernelSession: formatted error shows definition input");
    tests.expect(provenanceMessage.find("Called from In [2]") != std::string::npos
        && provenanceMessage.find("Called from In [3]") != std::string::npos,
        "KernelSession: formatted error shows nested call sites");

    kernel::KernelSession syntaxDocumentSession;
    const error::CalcError syntaxDocumentError = evaluateError(syntaxDocumentSession, "1 + )");
    tests.expect(syntaxDocumentError.document()
        && syntaxDocumentError.document()->inputNumber() == 1,
        "KernelSession: parser errors are attached to the pending input document");
    tests.expectEqual(syntaxDocumentSession.inputCount(), std::size_t{0},
        "KernelSession: parser errors do not consume an input number");
    tests.expectEqual(syntaxDocumentSession.nextInputNumber(), std::size_t{1},
        "KernelSession: parser errors leave the next prompt number unchanged");

    const error::CalcError malformedLiteralError = evaluateError(syntaxDocumentSession, "2..2+a");
    tests.expect(malformedLiteralError.type() == error::CalcErrorType::Syntax
        && malformedLiteralError.document()
        && malformedLiteralError.document()->inputNumber() == 1,
        "KernelSession: malformed literals report the still-pending input number");
    tests.expectEqual(syntaxDocumentSession.inputCount(), std::size_t{0},
        "KernelSession: malformed literals do not consume an input number");
    tests.expectEqual(evaluateAndFormat(syntaxDocumentSession, "2+3"), std::string{"5"},
        "KernelSession: valid input after a frontend error reuses the pending number");
    tests.expectEqual(syntaxDocumentSession.inputCount(), std::size_t{1},
        "KernelSession: successfully lowered input commits the input number");

    static_cast<void>(functionSession.evaluate("loop[x] := loop[x]"));
    const error::CalcError recursionError = evaluateError(functionSession, "loop[1]");
    tests.expect(recursionError.type() == error::CalcErrorType::Evaluation,
        "KernelSession: identical recursive calls are detected as a cycle");
    tests.expect(std::string{recursionError.what()}.find("Cyclic function call") != std::string::npos,
        "KernelSession: cyclic function call reports a dedicated message");

    // 引数が変化し続ける再帰は同値呼出し検出では止まらないため、論理深度制限で停止する。
    static_cast<void>(functionSession.evaluate("grow[x] := grow[x + 1]"));
    functionSession.setEvaluationDepthLimit(8);
    const error::CalcError depthError = evaluateError(functionSession, "grow[1]");
    tests.expect(depthError.type() == error::CalcErrorType::Evaluation,
        "KernelSession: changing non-terminating recursion is stopped by the evaluation depth limit");
    functionSession.setEvaluationDepthLimit(1024);

    // 明示的評価スタックなら、深いユーザー再帰でもC++のcall stackを消費しない。
    kernel::KernelSession deepRecursionSession;
    static_cast<void>(deepRecursionSession.evaluate(
        "countdown[n] := if[n == 0, 0, countdown[n - 1]]"));
    deepRecursionSession.setEvaluationDepthLimit(20000);
    tests.expectEqual(evaluateAndFormat(deepRecursionSession, "countdown[2000]"), std::string{"0"},
        "KernelSession: deep recursion uses the explicit evaluation stack");

    kernel::KernelSession calculusSession;
    tests.expectEqual(evaluateAndFormat(calculusSession, "D[x^2,x]"), std::string{"2x"},
        "KernelSession: D differentiates an exact polynomial");
    tests.expectEqual(evaluateAndFormat(calculusSession, "D[exp[x^2],x]"),
        std::string{"2x exp[x^2]"},
        "KernelSession: D applies the chain rule to Exp");
    tests.expectEqual(evaluateAndFormat(calculusSession, "D[log[10,x],x]"),
        std::string{"1/(x log[10])"},
        "KernelSession: D differentiates arbitrary-base log without redundant cancellation");
    tests.expectEqual(evaluateAndFormat(calculusSession, "D[sin[x],x]"),
        std::string{"cos[x]"},
        "KernelSession: D respects the default Radian angle semantics");
    tests.expectEqual(evaluateAndFormat(calculusSession, "D[asin[x],x]"),
        std::string{"1/sqrt[1-x^2]"},
        "KernelSession: inverse trig derivatives use Radian output by default");
    tests.expectEqual(evaluateAndFormat(calculusSession, "D[abs[x],x]"),
        std::string{"D[abs[x], x]"},
        "KernelSession: D does not invent a complex derivative for abs");
    static_cast<void>(calculusSession.evaluate("x := 7"));
    tests.expectEqual(evaluateAndFormat(calculusSession, "D[x^2,x]"), std::string{"2x"},
        "KernelSession: D holds its differentiation variable even when it has a value");

    kernel::KernelSession discreteSession;
    tests.expectEqual(evaluateAndFormat(discreteSession, "floor[-3/2]"), std::string{"-2"},
        "KernelSession: floor rounds exact rationals toward negative infinity");
    tests.expectEqual(evaluateAndFormat(discreteSession, "ceil[-3/2]"), std::string{"-1"},
        "KernelSession: ceil rounds exact rationals toward positive infinity");
    tests.expectEqual(evaluateAndFormat(discreteSession, "trunc[-3/2]"), std::string{"-1"},
        "KernelSession: trunc rounds exact rationals toward zero");
    tests.expectEqual(evaluateAndFormat(discreteSession, "round[5/2]"), std::string{"2"},
        "KernelSession: round uses nearest-even for exact half ties");
    tests.expectEqual(evaluateAndFormat(discreteSession, "round[3/2]"), std::string{"2"},
        "KernelSession: nearest-even rounds the other half tie to the even integer");
    tests.expectEqual(evaluateAndFormat(discreteSession, "frac[-3/2]"), std::string{"1/2"},
        "KernelSession: frac is the mathematical fractional part in [0,1)");
    tests.expectEqual(evaluateAndFormat(discreteSession, "floor[Pi]"), std::string{"3"},
        "KernelSession: floor proves a transcendental integer part through certified intervals");
    tests.expectEqual(evaluateAndFormat(discreteSession, "frac[Pi]"), std::string{"Pi-3"},
        "KernelSession: frac keeps the exact transcendental remainder after certification");
    tests.expectEqual(evaluateAndFormat(discreteSession, "gcd[84,126,210]"), std::string{"42"},
        "KernelSession: gcd accepts multiple exact integers");
    tests.expectEqual(evaluateAndFormat(discreteSession, "lcm[6,8,9]"), std::string{"72"},
        "KernelSession: lcm accepts multiple exact integers");
    tests.expectEqual(evaluateAndFormat(discreteSession, "mod[-5,3]"), std::string{"1"},
        "KernelSession: mod uses floor-division semantics");
    tests.expectEqual(evaluateAndFormat(discreteSession, "rem[-5,3]"), std::string{"-2"},
        "KernelSession: rem uses truncation-quotient semantics");
    tests.expectEqual(evaluateAndFormat(discreteSession, "quotient[-5,3]"), std::string{"-1"},
        "KernelSession: quotient truncates integer division toward zero");

    kernel::KernelSession fullSimplifySession;
    tests.expectEqual(
        evaluateAndFormat(fullSimplifySession, "fullSimplify[x^2+2*x+1]"),
        std::string{"(1+x)^2"},
        "KernelSession: FullSimplify searches a factored equivalent with lower cost");
    tests.expectEqual(
        evaluateAndFormat(fullSimplifySession, "fullSimplify[(x^2-1)/(x-1)]"),
        std::string{"(x^2-1)/(x-1)"},
        "KernelSession: FullSimplify does not cancel a possibly zero factor");
    tests.expectEqual(
        evaluateAndFormat(fullSimplifySession, "fullSimplify[(x^2-1)/(x-1),x!=1]"),
        std::string{"1+x"},
        "KernelSession: FullSimplify uses assumptions to justify candidate cancellation");

    kernel::KernelSession compatibilitySession;
    tests.expectEqual(evaluateAndFormat(compatibilitySession, "ln[E]"), std::string{"1"},
        "KernelSession: ln is a canonical Log alias");
    tests.expectEqual(evaluateAndFormat(compatibilitySession, "pow[2,10]"), std::string{"1024"},
        "KernelSession: pow aliases the exact Power builtin");
    tests.expectEqual(evaluateAndFormat(compatibilitySession, "fact[5]"), std::string{"120"},
        "KernelSession: fact aliases factorial evaluation");
    tests.expectEqual(evaluateAndFormat(compatibilitySession, "fract[7/3]"), std::string{"1/3"},
        "KernelSession: fract aliases the mathematical fractional-part function");
    tests.expectEqual(evaluateAndFormat(compatibilitySession, "real[3+4I]"), std::string{"3"},
        "KernelSession: real aliases re");
    tests.expectEqual(evaluateAndFormat(compatibilitySession, "imag[3+4I]"), std::string{"4"},
        "KernelSession: imag aliases im");
    tests.expectEqual(evaluateAndFormat(compatibilitySession, "mag[3+4I]"), std::string{"5"},
        "KernelSession: mag aliases abs");
    tests.expectEqual(evaluateAndFormat(compatibilitySession, "unit[3+4I]"),
        std::string{"3/5+4/5I"},
        "KernelSession: unit aliases complex sign without a duplicate implementation");

    kernel::KernelSession utilitySession;
    tests.expectEqual(evaluateAndFormat(utilitySession, "cbrt[-8]"), std::string{"-2"},
        "KernelSession: cbrt selects the real cube root for a negative real");
    tests.expectEqual(evaluateAndFormat(utilitySession, "cbrt[1/8]"), std::string{"1/2"},
        "KernelSession: cbrt preserves exact rational cube roots");
    tests.expectEqual(evaluateAndFormat(utilitySession, "cbrt[-2]"), std::string{"-cbrt[2]"},
        "KernelSession: cbrt canonicalizes odd symmetry without using principal Power");
    tests.expectEqual(evaluateAndFormat(utilitySession, "N[cbrt[2],20]"),
        std::string{"1.2599210498948731648"},
        "KernelSession: cbrt has certified numerical evaluation");
    tests.expectEqual(evaluateAndFormat(utilitySession, "hypot[3,4]"), std::string{"5"},
        "KernelSession: hypot evaluates an exact Pythagorean triple");
    tests.expectEqual(evaluateAndFormat(utilitySession, "hypot[1,1]"), std::string{"sqrt[2]"},
        "KernelSession: hypot keeps an irrational exact result as a normalized radical");
    tests.expectEqual(evaluateAndFormat(utilitySession, "cis[60 Deg]"),
        std::string{"1/2+I sqrt[3]/2"},
        "KernelSession: cis accepts an explicit Degree angle");
    tests.expectEqual(evaluateAndFormat(utilitySession, "DtoR[180]"), std::string{"Pi"},
        "KernelSession: DtoR preserves exact Pi");
    tests.expectEqual(evaluateAndFormat(utilitySession, "DtoG[90]"), std::string{"100"},
        "KernelSession: DtoG is exact rational arithmetic");
    tests.expectEqual(evaluateAndFormat(utilitySession, "RtoD[Pi]"), std::string{"180"},
        "KernelSession: RtoD recognizes exact rational Pi multiples");
    tests.expectEqual(evaluateAndFormat(utilitySession, "RtoG[Pi/2]"), std::string{"100"},
        "KernelSession: RtoG recognizes exact rational Pi multiples");
    tests.expectEqual(evaluateAndFormat(utilitySession, "GtoD[200]"), std::string{"180"},
        "KernelSession: GtoD is exact rational arithmetic");
    tests.expectEqual(evaluateAndFormat(utilitySession, "GtoR[200]"), std::string{"Pi"},
        "KernelSession: GtoR preserves exact Pi");
    tests.expectEqual(evaluateAndFormat(utilitySession, "nextpow2[9]"), std::string{"4"},
        "KernelSession: nextpow2 returns the smallest exponent whose power covers x");
    tests.expectEqual(evaluateAndFormat(utilitySession, "nextpow2[1/10]"), std::string{"-3"},
        "KernelSession: nextpow2 extends exactly to positive rationals below one");
    tests.expectEqual(evaluateAndFormat(utilitySession, "nextpow2[Pi]"), std::string{"2"},
        "KernelSession: nextpow2 can certify symbolic exact constants");
    const error::CalcError nextPow2Domain = evaluateError(utilitySession, "nextpow2[0]");
    tests.expect(nextPow2Domain.type() == error::CalcErrorType::Domain,
        "KernelSession: nextpow2 rejects nonpositive arguments");
    tests.expectEqual(evaluateAndFormat(utilitySession, "D[cbrt[x],x]"),
        std::string{"1/(3cbrt[x]^2)"},
        "KernelSession: D differentiates real cube root without rewriting it as principal Power");
    tests.expectEqual(evaluateAndFormat(utilitySession, "D[cis[x],x]"),
        std::string{"I cis[x]"},
        "KernelSession: D applies the default Radian scale to cis");

    tests.expectEqual(evaluateAndFormat(utilitySession, "solve[cbrt[y]*x+1==0,x]"),
        std::string{"cases[{x==-1/cbrt[y]} if cbrt[y]!=0; {} if cbrt[y]==0] if y in Real"},
        "KernelSession: Solver inherits cbrt real-domain definedness from MathRegistry");

    kernel::KernelSession factorialSession;
    tests.expectEqual(evaluateAndFormat(factorialSession, "-5!"), std::string{"-120"},
        "KernelSession: factorial binds before unary minus");
    const error::CalcError negativeFactorial = evaluateError(factorialSession, "(-5)!");
    tests.expect(negativeFactorial.type() == error::CalcErrorType::Domain,
        "KernelSession: parenthesized negative factorial is undefined");

    session.clearHistory();
    tests.expectEqual(session.historySize(), std::size_t{0},
        "KernelSession: clears history independently");
    tests.expectEqual(evaluateAndFormat(session, "x"), std::string{"3"},
        "KernelSession: clearing history preserves definitions");

    session.clearDefinitions();
    tests.expectEqual(evaluateAndFormat(session, "x"), std::string{"x"},
        "KernelSession: clearing definitions leaves x as an unbound free symbol");

    functionSession.clearDefinitions();
    tests.expectEqual(functionSession.userFunctions().size(), std::size_t{0},
        "KernelSession: clearDefinitions removes user functions");
    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(functionSession.evaluate("f[1]"));
    }, "KernelSession: cleared user function is no longer callable");

    session.setEvaluationDepthLimit(2048);
    tests.expectEqual(session.evaluationDepthLimit(), std::size_t{2048},
        "KernelSession: forwards evaluation depth limit");

    kernel::KernelSession independentReset;
    kernel::KernelSession randomControl;
    tests.expectEqual(evaluateAndFormat(independentReset, "randSeed[424242]"), std::string{"424242"},
        "KernelSession: independent reset test seeds the random stream");
    tests.expectEqual(evaluateAndFormat(randomControl, "randSeed[424242]"), std::string{"424242"},
        "KernelSession: independent reset control uses the same random seed");
    tests.expectEqual(evaluateAndFormat(independentReset, "rand[]"), evaluateAndFormat(randomControl, "rand[]"),
        "KernelSession: seeded streams agree before independent reset");
    tests.expectEqual(evaluateAndFormat(independentReset, "temporary:=17"), std::string{"17"},
        "KernelSession: independent reset test creates a temporary definition");
    independentReset.resetForIndependentEvaluation();
    tests.expectEqual(independentReset.historySize(), std::size_t{0},
        "KernelSession: independent reset clears history");
    tests.expectEqual(independentReset.inputCount(), std::size_t{0},
        "KernelSession: independent reset clears input numbering");
    tests.expectEqual(evaluateAndFormat(independentReset, "temporary"), std::string{"temporary"},
        "KernelSession: independent reset clears definitions");
    tests.expectEqual(evaluateAndFormat(independentReset, "rand[]"), evaluateAndFormat(randomControl, "rand[]"),
        "KernelSession: independent reset preserves the deterministic random stream");

    session.reset();
    tests.expectEqual(session.historySize(), std::size_t{0},
        "KernelSession: reset clears history");
    tests.expectEqual(session.inputCount(), std::size_t{0},
        "KernelSession: reset clears input numbering");
}

} // namespace mmcal::tests
