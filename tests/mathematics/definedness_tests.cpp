// 式の定義条件生成の回帰テスト
#include "definedness_tests.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/definedness.hpp"
#include "mathematics/math_registry.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

#include <cstdint>
#include <utility>
#include <variant>
#include <vector>

namespace mmcal::tests {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::RelationKind;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr call(
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id,
    std::vector<Expr> arguments) {
    return Expr::call(builtins.symbol(id), std::move(arguments));
}

} // namespace

void runDefinednessTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const auto builtins = evaluation::BuiltinRegistry::defaults(symbols);
    const auto mathematics = mathematics::MathRegistry::defaults(symbols, builtins);
    const Expr y{symbols.intern("y")};
    const Expr x{symbols.intern("x")};

    const Expr tanY = call(builtins, BuiltinId::Tan, {y});
    const auto tanConditions = mathematics::expressionDomainConditions(
        tanY, builtins, mathematics);
    const Expr cosY = call(builtins, BuiltinId::Cos, {y});
    const auto expectedTan = mathematics::relation(RelationKind::NotEqual, cosY, integer(0));
    tests.expect(tanConditions && tanConditions->size() == 1
        && tanConditions->contains(expectedTan),
        "Definedness: tan requires cos(argument) != 0 from MathRegistry metadata");

    const Expr atanY = call(builtins, BuiltinId::Atan, {y});
    const Expr tanAtanY = call(builtins, BuiltinId::Tan, {atanY});
    const auto tanAtanConditions = mathematics::expressionDomainConditions(
        tanAtanY, builtins, mathematics);
    const Expr ySquarePlusOne = call(builtins, BuiltinId::Add, {
        integer(1), call(builtins, BuiltinId::Power, {y, integer(2)})});
    const auto expectedAtan = mathematics::relation(
        RelationKind::NotEqual, ySquarePlusOne, integer(0));
    tests.expect(tanAtanConditions && tanAtanConditions->size() == 1
        && tanAtanConditions->contains(expectedAtan),
        "Definedness: tan[atan[y]] inherits only atan branch-point exclusions");

    const Expr atanhY = call(builtins, BuiltinId::Atanh, {y});
    const Expr tanhAtanhY = call(builtins, BuiltinId::Tanh, {atanhY});
    const auto tanhAtanhConditions = mathematics::expressionDomainConditions(
        tanhAtanhY, builtins, mathematics);
    const Expr oneMinusYSquare = call(builtins, BuiltinId::Subtract, {
        integer(1), call(builtins, BuiltinId::Power, {y, integer(2)})});
    const auto expectedAtanh = mathematics::relation(
        RelationKind::NotEqual, oneMinusYSquare, integer(0));
    tests.expect(tanhAtanhConditions && tanhAtanhConditions->size() == 1
        && tanhAtanhConditions->contains(expectedAtanh),
        "Definedness: tanh[atanh[y]] inherits only atanh branch-point exclusions");

    const Expr logY = call(builtins, BuiltinId::Log, {y});
    const auto logConditions = mathematics::expressionDomainConditions(
        logY, builtins, mathematics);
    const auto expectedLog = mathematics::relation(RelationKind::NotEqual, y, integer(0));
    tests.expect(logConditions && logConditions->size() == 1
        && logConditions->contains(expectedLog),
        "Definedness: principal log requires a nonzero argument");

    const Expr b{symbols.intern("b")};
    const Expr logBase = call(builtins, BuiltinId::Log, {b, y});
    const auto logBaseConditions = mathematics::expressionDomainConditions(
        logBase, builtins, mathematics);
    const auto expectedBaseNonZero = mathematics::relation(
        RelationKind::NotEqual, b, integer(0));
    const Expr baseMinusOne = call(builtins, BuiltinId::Subtract, {b, integer(1)});
    const auto expectedBaseNotOne = mathematics::relation(
        RelationKind::NotEqual, baseMinusOne, integer(0));
    tests.expect(logBaseConditions && logBaseConditions->size() == 3
        && logBaseConditions->contains(expectedBaseNonZero)
        && logBaseConditions->contains(expectedBaseNotOne)
        && logBaseConditions->contains(expectedLog),
        "Definedness: arbitrary-base log requires base != 0,1 and value != 0");

    const Expr reciprocal = call(builtins, BuiltinId::Power, {y, integer(-1)});
    const auto reciprocalConditions = mathematics::expressionDomainConditions(
        reciprocal, builtins, mathematics);
    tests.expect(reciprocalConditions && reciprocalConditions->size() == 1
        && reciprocalConditions->contains(expectedLog),
        "Definedness: negative integer powers require a nonzero base");

    const Expr productDenominator = call(builtins, BuiltinId::Multiply, {
        x, call(builtins, BuiltinId::Add, {y, integer(1)})});
    const Expr reciprocalProduct = call(builtins, BuiltinId::Divide, {integer(1), productDenominator});
    const auto productConditions = mathematics::expressionDomainConditions(
        reciprocalProduct, builtins, mathematics);
    const auto expectedXNonZero = mathematics::relation(RelationKind::NotEqual, x, integer(0));
    const Expr yPlusOne = call(builtins, BuiltinId::Add, {y, integer(1)});
    const auto expectedYPlusOneNonZero = mathematics::relation(
        RelationKind::NotEqual, yPlusOne, integer(0));
    tests.expect(productConditions && productConditions->size() == 2
        && productConditions->contains(expectedXNonZero)
        && productConditions->contains(expectedYPlusOneNonZero),
        "Definedness: nonzero products decompose into reusable factor conditions");

    const Expr squareDenominator = call(builtins, BuiltinId::Power, {yPlusOne, integer(2)});
    const Expr reciprocalSquare = call(builtins, BuiltinId::Divide, {integer(1), squareDenominator});
    const auto squareConditions = mathematics::expressionDomainConditions(
        reciprocalSquare, builtins, mathematics);
    tests.expect(squareConditions && squareConditions->size() == 1
        && squareConditions->contains(expectedYPlusOneNonZero),
        "Definedness: nonzero integer powers reduce to a nonzero-base condition");

    const Expr zeroPower = call(builtins, BuiltinId::Power, {y, integer(0)});
    const auto zeroPowerConditions = mathematics::expressionDomainConditions(
        zeroPower, builtins, mathematics);
    tests.expect(zeroPowerConditions && zeroPowerConditions->size() == 1
        && zeroPowerConditions->contains(expectedLog),
        "Definedness: exponent zero retains the 0^0 exclusion");

    const Expr positivePower = call(builtins, BuiltinId::Power, {y, integer(2)});
    const auto positivePowerConditions = mathematics::expressionDomainConditions(
        positivePower, builtins, mathematics);
    tests.expect(positivePowerConditions && positivePowerConditions->empty(),
        "Definedness: positive integer powers do not require a nonzero base");

    const Expr positiveRationalPower = call(builtins, BuiltinId::Power, {
        y, Expr{Number{numeric::Rational{BigInt{3}, BigInt{2}}}}});
    const auto positiveRationalPowerConditions = mathematics::expressionDomainConditions(
        positiveRationalPower, builtins, mathematics);
    tests.expect(positiveRationalPowerConditions && positiveRationalPowerConditions->empty(),
        "Definedness: positive exact rational powers are defined at a zero base");

    const Expr negativeRationalPower = call(builtins, BuiltinId::Power, {
        y, Expr{Number{numeric::Rational{BigInt{-3}, BigInt{2}}}}});
    const auto negativeRationalPowerConditions = mathematics::expressionDomainConditions(
        negativeRationalPower, builtins, mathematics);
    tests.expect(negativeRationalPowerConditions && negativeRationalPowerConditions->size() == 1
        && negativeRationalPowerConditions->contains(expectedLog),
        "Definedness: negative exact rational powers require a nonzero base");

    const Expr cbrtY = call(builtins, BuiltinId::Cbrt, {y});
    const auto cbrtConditions = mathematics::expressionDomainConditions(
        cbrtY, builtins, mathematics);
    const auto expectedRealY = mathematics::elementOf(y, mathematics::NumericDomain::Real);
    tests.expect(cbrtConditions && cbrtConditions->size() == 1
        && cbrtConditions->contains(expectedRealY),
        "Definedness: real cbrt requires a real argument");

    const Expr hypotXY = call(builtins, BuiltinId::Hypot, {x, y});
    const auto hypotConditions = mathematics::expressionDomainConditions(
        hypotXY, builtins, mathematics);
    const auto expectedRealX = mathematics::elementOf(x, mathematics::NumericDomain::Real);
    tests.expect(hypotConditions && hypotConditions->size() == 2
        && hypotConditions->contains(expectedRealX)
        && hypotConditions->contains(expectedRealY),
        "Definedness: hypot requires both arguments to be real");

    const Expr expY = call(builtins, BuiltinId::Exp, {y});
    const auto expConditions = mathematics::expressionDomainConditions(
        expY, builtins, mathematics);
    tests.expect(expConditions && expConditions->empty(),
        "Definedness: entire functions add no domain conditions");

    const Expr reciprocalX = call(builtins, BuiltinId::Divide, {integer(1), x});
    const Expr arrayValue = Expr::array({2}, {reciprocalX, logY});
    const auto arrayConditions = mathematics::expressionDomainConditions(
        arrayValue, builtins, mathematics);
    tests.expect(arrayConditions && arrayConditions->size() == 2
        && arrayConditions->contains(expectedXNonZero)
        && arrayConditions->contains(expectedLog),
        "Definedness: Array values collect element domain conditions");

    const Expr raggedValue = Expr::list({
        Expr::array({1}, {reciprocalX}), logY});
    const auto raggedConditions = mathematics::expressionDomainConditions(
        raggedValue, builtins, mathematics);
    tests.expect(raggedConditions && raggedConditions->size() == 2
        && raggedConditions->contains(expectedXNonZero)
        && raggedConditions->contains(expectedLog),
        "Definedness: ragged brace Lists collect nested element domain conditions");

    const Expr gammaY = call(builtins, BuiltinId::Gamma, {y});
    tests.expect(!mathematics::expressionDomainConditions(gammaY, builtins, mathematics),
        "Definedness: Gamma pole complement stays unresolved rather than dropping infinitely many exclusions");

    const Expr digammaThree = call(builtins, BuiltinId::Digamma, {integer(3)});
    const auto digammaThreeConditions = mathematics::expressionDomainConditions(
        digammaThree, builtins, mathematics);
    tests.expect(digammaThreeConditions && digammaThreeConditions->empty(),
        "Definedness: exact positive Gamma-family arguments are known to avoid poles");

    const Expr digammaHalf = call(builtins, BuiltinId::Digamma, {
        Expr{Number{numeric::Rational{BigInt{1}, BigInt{2}}}}});
    const auto digammaHalfConditions = mathematics::expressionDomainConditions(
        digammaHalf, builtins, mathematics);
    tests.expect(digammaHalfConditions && digammaHalfConditions->empty(),
        "Definedness: exact noninteger Gamma-family arguments are known to avoid poles");

    const Expr digammaPole = call(builtins, BuiltinId::Digamma, {integer(0)});
    tests.expect(!mathematics::expressionDomainConditions(digammaPole, builtins, mathematics),
        "Definedness: exact Gamma-family poles are not certified as defined");

    const Expr zetaY = call(builtins, BuiltinId::Zeta, {y});
    const auto zetaConditions = mathematics::expressionDomainConditions(
        zetaY, builtins, mathematics);
    const auto expectedZeta = mathematics::relation(RelationKind::NotEqual, y, integer(1));
    tests.expect(zetaConditions && zetaConditions->size() == 1
        && zetaConditions->contains(expectedZeta),
        "Definedness: zeta has the exact finite-plane condition s != 1");

    const Expr betaXY = call(builtins, BuiltinId::Beta, {x, y});
    const auto betaConditions = mathematics::expressionDomainConditions(
        betaXY, builtins, mathematics);
    const auto expectedPositiveX = mathematics::relation(RelationKind::Greater, x, integer(0));
    const auto expectedPositiveY = mathematics::relation(RelationKind::Greater, y, integer(0));
    tests.expect(betaConditions && betaConditions->size() == 4
        && betaConditions->contains(expectedRealX)
        && betaConditions->contains(expectedRealY)
        && betaConditions->contains(expectedPositiveX)
        && betaConditions->contains(expectedPositiveY),
        "Definedness: current Beta contract requires both arguments to be positive reals");

    const Expr atan2 = call(builtins, BuiltinId::Atan2, {y, integer(1)});
    tests.expect(!mathematics::expressionDomainConditions(atan2, builtins, mathematics),
        "Definedness: unsupported compound predicates remain unresolved rather than weakened");
}

} // namespace mmcal::tests
