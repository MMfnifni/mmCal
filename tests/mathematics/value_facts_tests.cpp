// 式の符号・実数性などの事実推論の回帰テスト
#include "value_facts_tests.hpp"

#include "kernel/kernel_session.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/value_facts.hpp"
#include "test_framework.hpp"

namespace mmcal::tests {

void runValueFactsTests(TestRunner& tests) {
    kernel::KernelSession session;

    const auto pi = session.evaluate("Pi");
    const auto piFacts = mathematics::inferValueFacts(
        pi, session.builtinRegistry(), session.mathRegistry());
    tests.expect(piFacts.domain == mathematics::NumericDomain::Real
        && piFacts.sign == mathematics::RealSign::Positive
        && piFacts.exact,
        "ValueFacts: Pi is known exact positive real");
    tests.expect(piFacts.provablyNonInteger && piFacts.provablyNonRational,
        "ValueFacts: Pi transcendence implies non-integer and non-rational knowledge");

    const auto half = session.evaluate("1/2");
    const auto halfFacts = mathematics::inferValueFacts(
        half, session.builtinRegistry(), session.mathRegistry());
    tests.expect(halfFacts.domain == mathematics::NumericDomain::Rational
        && halfFacts.provablyNonInteger && !halfFacts.provablyNonRational,
        "ValueFacts: a non-integral exact rational is provably non-integer but remains rational");

    const auto phi = session.evaluate("Phi");
    const auto phiFacts = mathematics::inferValueFacts(
        phi, session.builtinRegistry(), session.mathRegistry());
    tests.expect(phiFacts.provablyNonInteger && phiFacts.provablyNonRational,
        "ValueFacts: known irrational algebraic constants refute Rational and Integer membership");

    const auto i = session.evaluate("I");
    const auto iFacts = mathematics::inferValueFacts(
        i, session.builtinRegistry(), session.mathRegistry());
    tests.expect(iFacts.domain == mathematics::NumericDomain::Complex
        && iFacts.provablyNonReal,
        "ValueFacts: I is provably non-real complex");

    const auto sqrtPi = session.evaluate("sqrt[Pi]");
    const auto sqrtPiFacts = mathematics::inferValueFacts(
        sqrtPi, session.builtinRegistry(), session.mathRegistry());
    tests.expect(sqrtPiFacts.domain == mathematics::NumericDomain::Real
        && sqrtPiFacts.sign == mathematics::RealSign::Positive,
        "ValueFacts: sqrt of positive real constant stays positive real");

    const auto sqrtNegative = session.evaluate("sqrt[-2]");
    const auto sqrtNegativeFacts = mathematics::inferValueFacts(
        sqrtNegative, session.builtinRegistry(), session.mathRegistry());
    tests.expect(sqrtNegativeFacts.domain == mathematics::NumericDomain::Complex
        && sqrtNegativeFacts.provablyNonReal,
        "ValueFacts: sqrt of negative real promotes to non-real complex");

    const auto sinPi = session.evaluate("sin[Pi]");
    const auto sinFacts = mathematics::inferValueFacts(
        sinPi, session.builtinRegistry(), session.mathRegistry());
    tests.expect(sinFacts.isProvablyReal(),
        "ValueFacts: sin maps a known real argument to real");

    const auto acoshTwoFacts = mathematics::inferValueFacts(
        session.evaluate("acosh[2]"), session.builtinRegistry(), session.mathRegistry());
    tests.expect(acoshTwoFacts.domain == mathematics::NumericDomain::Real
        && acoshTwoFacts.sign == mathematics::RealSign::Positive,
        "ValueFacts: principal acosh is positive for exact real arguments above one");
    const auto acosOneFacts = mathematics::inferValueFacts(
        session.evaluate("acos[1]"), session.builtinRegistry(), session.mathRegistry());
    tests.expect(acosOneFacts.sign == mathematics::RealSign::Zero,
        "ValueFacts: principal acos is zero at its exact upper endpoint");
    const auto asinNegativeFacts = mathematics::inferValueFacts(
        session.evaluate("asin[-1/2]"), session.builtinRegistry(), session.mathRegistry());
    tests.expect(asinNegativeFacts.sign == mathematics::RealSign::Negative,
        "ValueFacts: principal asin preserves sign on its exact real interval");

    const auto x = session.evaluate("x");
    mathematics::AssumptionSet realX;
    realX.add(mathematics::elementOf(x, mathematics::NumericDomain::Real));
    const auto squareFacts = mathematics::inferValueFacts(
        session.evaluate("x^2"), session.builtinRegistry(), session.mathRegistry(), realX);
    tests.expect(squareFacts.isProvablyReal()
        && squareFacts.sign == mathematics::RealSign::NonNegative,
        "ValueFacts: an even power of an unconstrained real base is nonnegative");
    const auto reciprocalSquareFacts = mathematics::inferValueFacts(
        session.evaluate("x^(-2)"), session.builtinRegistry(), session.mathRegistry(), realX);
    tests.expect(reciprocalSquareFacts.isProvablyReal()
        && reciprocalSquareFacts.sign == mathematics::RealSign::Positive,
        "ValueFacts: a negative even power is positive wherever it is defined");
}

} // namespace mmcal::tests
