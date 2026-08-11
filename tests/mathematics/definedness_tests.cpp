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

    const Expr tanY = call(builtins, BuiltinId::Tan, {y});
    const auto tanConditions = mathematics::expressionDomainConditions(
        tanY, builtins, mathematics);
    const Expr cosY = call(builtins, BuiltinId::Cos, {y});
    const auto expectedTan = mathematics::relation(RelationKind::NotEqual, cosY, integer(0));
    tests.expect(tanConditions && tanConditions->size() == 1
        && tanConditions->contains(expectedTan),
        "Definedness: tan requires cos(argument) != 0 from MathRegistry metadata");

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

    const Expr cbrtY = call(builtins, BuiltinId::Cbrt, {y});
    const auto cbrtConditions = mathematics::expressionDomainConditions(
        cbrtY, builtins, mathematics);
    const auto expectedRealY = mathematics::elementOf(y, mathematics::NumericDomain::Real);
    tests.expect(cbrtConditions && cbrtConditions->size() == 1
        && cbrtConditions->contains(expectedRealY),
        "Definedness: real cbrt requires a real argument");

    const Expr x{symbols.intern("x")};
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

    const Expr gammaY = call(builtins, BuiltinId::Gamma, {y});
    tests.expect(!mathematics::expressionDomainConditions(gammaY, builtins, mathematics),
        "Definedness: Gamma pole complement stays unresolved rather than dropping infinitely many exclusions");

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
