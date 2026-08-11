// 式の符号・実数性などの事実推論の回帰テスト
#include "value_facts_tests.hpp"

#include "kernel/kernel_session.hpp"
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
}

} // namespace mmcal::tests
