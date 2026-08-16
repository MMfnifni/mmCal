// 局所仮定を含む数学知識コンテキストの回帰テスト
#include "knowledge_context_tests.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

namespace mmcal::tests {
namespace {

using expression::Expr;
using mathematics::NumericDomain;
using mathematics::Predicate;
using mathematics::RealSign;
using mathematics::RelationKind;
using mathematics::TruthValue;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

} // namespace

void runKnowledgeContextTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const auto builtins = evaluation::BuiltinRegistry::defaults(symbols);
    const auto math = mathematics::MathRegistry::defaults(symbols, builtins);
    mathematics::AssumptionSet none;
    const mathematics::KnowledgeContext baseKnowledge{builtins, math, none};

    const Expr pi{symbols.intern("Pi")};
    tests.expect(
        baseKnowledge.prove(mathematics::elementOf(pi, NumericDomain::Real))
            == TruthValue::True,
        "KnowledgeContext: Pi is provably real from permanent mathematical knowledge");
    tests.expect(
        baseKnowledge.prove(mathematics::elementOf(pi, NumericDomain::Rational))
            == TruthValue::False,
        "KnowledgeContext: Pi transcendence refutes Rational membership");
    tests.expect(
        baseKnowledge.prove(mathematics::elementOf(pi, NumericDomain::Integer))
            == TruthValue::False,
        "KnowledgeContext: Pi transcendence refutes Integer membership");

    const Expr half{numeric::Rational{BigInt{1}, BigInt{2}}};
    tests.expect(
        baseKnowledge.prove(mathematics::elementOf(half, NumericDomain::Integer))
            == TruthValue::False,
        "KnowledgeContext: a non-integral exact rational is provably not an integer");


    const Expr sqrtTwo = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sqrt), {integer(2)});
    const Expr sqrtThree = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sqrt), {integer(3)});
    const Expr rootSqrtTwo = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Root), {
            Expr::rationalArray({3}, {
                numeric::Rational{BigInt{-2}},
                numeric::Rational{},
                numeric::Rational{BigInt{1}}}),
            integer(2)});
    tests.expect(
        baseKnowledge.prove(mathematics::relation(
            RelationKind::Equal, rootSqrtTwo, sqrtTwo)) == TruthValue::True,
        "KnowledgeContext: canonical Root and equivalent radical share exact algebraic equality");
    tests.expect(
        baseKnowledge.prove(mathematics::relation(
            RelationKind::Less, sqrtTwo, sqrtThree)) == TruthValue::True,
        "KnowledgeContext: exact radical ordering reuses algebraic-number comparison");
    tests.expect(
        baseKnowledge.prove(mathematics::elementOf(sqrtTwo, NumericDomain::Rational))
            == TruthValue::False,
        "KnowledgeContext: irrational radicals are rejected from Rational by the algebraic bridge");

    const Expr cbrtTwo = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Cbrt), {integer(2)});
    const Expr rootCbrtTwo = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Root), {
            Expr::rationalArray({4}, {
                numeric::Rational{BigInt{-2}},
                numeric::Rational{},
                numeric::Rational{},
                numeric::Rational{BigInt{1}}}),
            integer(1)});
    tests.expect(
        baseKnowledge.prove(mathematics::relation(
            RelationKind::Equal, rootCbrtTwo, cbrtTwo)) == TruthValue::True,
        "KnowledgeContext: real cubic Root and cbrt share exact algebraic identity");

    const Expr imaginary{Number::complex(BigInt{0}, BigInt{1})};
    tests.expect(
        baseKnowledge.prove(mathematics::elementOf(imaginary, NumericDomain::Real))
            == TruthValue::False,
        "KnowledgeContext: exact non-real complex is provably not real");

    mathematics::AssumptionSet contradictory;
    contradictory.add(mathematics::elementOf(imaginary, NumericDomain::Real));
    const mathematics::KnowledgeContext contradictoryKnowledge{builtins, math, contradictory};
    tests.expect(
        contradictoryKnowledge.prove(mathematics::elementOf(imaginary, NumericDomain::Real))
            == TruthValue::False,
        "KnowledgeContext: an assumption cannot override a proven mathematical contradiction");

    const Expr x{symbols.intern("x")};
    tests.expect(
        baseKnowledge.prove(mathematics::elementOf(x, NumericDomain::Real))
            == TruthValue::Unknown,
        "KnowledgeContext: an unconstrained symbol is not silently assumed real");

    mathematics::AssumptionSet assumptions;
    assumptions.add(mathematics::elementOf(x, NumericDomain::Real));
    assumptions.add(mathematics::relation(RelationKind::Greater, x, integer(0)));
    const mathematics::KnowledgeContext assumedKnowledge{builtins, math, assumptions};

    const auto xFacts = assumedKnowledge.facts(x);
    tests.expect(
        xFacts.domain == NumericDomain::Real
            && xFacts.sign == RealSign::Positive
            && xFacts.exact,
        "KnowledgeContext: explicit x in Real and x > 0 assumptions become value facts");

    const Expr sqrtX = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sqrt),
        {x});
    const auto sqrtFacts = assumedKnowledge.facts(sqrtX);
    tests.expect(
        sqrtFacts.domain == NumericDomain::Real
            && sqrtFacts.sign == RealSign::Positive,
        "KnowledgeContext: assumptions propagate through principal sqrt");

    const Expr minusTwo = integer(-2);
    const Expr halfPower = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Power),
        {minusTwo, half});
    const auto halfPowerFacts = baseKnowledge.facts(halfPower);
    tests.expect(
        halfPowerFacts.domain == NumericDomain::Complex
            && halfPowerFacts.provablyNonReal,
        "KnowledgeContext: Power[x,1/2] shares principal sqrt domain semantics");

    tests.expect(
        assumedKnowledge.prove(
            mathematics::relation(RelationKind::NotEqual, x, integer(0)))
            == TruthValue::True,
        "KnowledgeContext: x > 0 proves x != 0 without numerical approximation");

    mathematics::AssumptionSet complexNonZero;
    complexNonZero.add(mathematics::elementOf(x, NumericDomain::Complex));
    complexNonZero.add(mathematics::relation(RelationKind::NotEqual, x, integer(0)));
    const mathematics::KnowledgeContext complexNonZeroKnowledge{builtins, math, complexNonZero};
    tests.expect(
        complexNonZeroKnowledge.facts(x).domain == NumericDomain::Complex,
        "KnowledgeContext: x != 0 does not incorrectly imply that a complex variable is real");
    tests.expect(
        complexNonZeroKnowledge.prove(mathematics::relation(RelationKind::Equal, x, integer(0)))
            == TruthValue::False,
        "KnowledgeContext: explicit != assumptions refute the complementary equality");

    tests.expect(
        assumedKnowledge.prove(mathematics::relation(RelationKind::Less, integer(0), x))
            == TruthValue::True,
        "KnowledgeContext: relation assumptions are recognized in reversed form");

    const Expr expX = Expr::call(builtins.symbol(evaluation::BuiltinId::Exp), {x});
    tests.expect(
        baseKnowledge.prove(mathematics::relation(RelationKind::NotEqual, expX, integer(0)))
            == TruthValue::True,
        "KnowledgeContext: globally defined zero-free functions prove nonzero");

    const Expr gammaX = Expr::call(builtins.symbol(evaluation::BuiltinId::Gamma), {x});
    tests.expect(
        baseKnowledge.prove(mathematics::relation(RelationKind::NotEqual, gammaX, integer(0)))
            == TruthValue::Unknown,
        "KnowledgeContext: zero-free knowledge does not erase unresolved function poles");
}

} // namespace mmcal::tests
