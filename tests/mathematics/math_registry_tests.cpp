// 函数のdomain・逆函数・周期などの数学metadataの回帰テスト
#include "math_registry_tests.hpp"

#include "evaluation/builtin_registry.hpp"
#include "mathematics/math_registry.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

namespace mmcal::tests {

void runMathRegistryTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const evaluation::BuiltinRegistry builtins = evaluation::BuiltinRegistry::defaults(symbols);
    const mathematics::MathRegistry mathematics = mathematics::MathRegistry::defaults(symbols, builtins);

    const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
    tests.expect(pi != nullptr, "MathRegistry: Pi is registered");
    tests.expect(pi && pi->properties.exact && pi->properties.real
        && pi->properties.positive && pi->properties.irrational,
        "MathRegistry: Pi carries exact real positive irrational properties");
    tests.expect(pi && pi->properties.arithmeticClass == mathematics::ArithmeticClass::Transcendental,
        "MathRegistry: Pi is classified as transcendental");

    const auto* e = mathematics.findConstant(mathematics::ConstantId::E);
    tests.expect(e && e->properties.exact && e->properties.real && e->properties.positive
        && e->properties.arithmeticClass == mathematics::ArithmeticClass::Transcendental,
        "MathRegistry: E is an exact positive transcendental symbolic constant");

    const auto* phi = mathematics.findConstant(mathematics::ConstantId::Phi);
    tests.expect(phi && phi->properties.arithmeticClass == mathematics::ArithmeticClass::Algebraic,
        "MathRegistry: Phi is classified as algebraic");

    const auto* cbrt = mathematics.findFunction(mathematics::FunctionId::Cbrt);
    const auto* hypot = mathematics.findFunction(mathematics::FunctionId::Hypot);
    const auto* cis = mathematics.findFunction(mathematics::FunctionId::Cis);
    const auto* dtor = mathematics.findFunction(mathematics::FunctionId::DegreeToRadian);
    const auto* rtod = mathematics.findFunction(mathematics::FunctionId::RadianToDegree);
    const auto* sqrt = mathematics.findFunction(mathematics::FunctionId::Sqrt);
    const auto* abs = mathematics.findFunction(mathematics::FunctionId::Abs);
    const auto* sign = mathematics.findFunction(mathematics::FunctionId::Sign);
    const auto* re = mathematics.findFunction(mathematics::FunctionId::Re);
    const auto* im = mathematics.findFunction(mathematics::FunctionId::Im);
    const auto* conj = mathematics.findFunction(mathematics::FunctionId::Conj);
    const auto* sin = mathematics.findFunction(mathematics::FunctionId::Sin);
    const auto* cos = mathematics.findFunction(mathematics::FunctionId::Cos);
    const auto* tan = mathematics.findFunction(mathematics::FunctionId::Tan);
    tests.expect(cbrt
        && cbrt->parity == mathematics::FunctionParity::Odd
        && cbrt->domainRule == mathematics::FunctionDomainRule::RealToReal
        && cbrt->definednessRule == mathematics::FunctionDefinednessRule::ArgumentReal,
        "MathRegistry: cbrt is the single-valued real cube-root function");
    tests.expect(hypot && hypot->arity == 2
        && hypot->domainRule == mathematics::FunctionDomainRule::RealPairToReal
        && hypot->definednessRule == mathematics::FunctionDefinednessRule::ArgumentsReal,
        "MathRegistry: hypot requires two real arguments");
    tests.expect(cis && cis->periodTurns && cis->periodTurns->toString() == "1"
        && cis->branchRule == mathematics::FunctionBranchRule::SingleValued,
        "MathRegistry: cis is single-valued with a one-turn period");
    tests.expect(dtor && dtor->parity == mathematics::FunctionParity::Odd
        && rtod && rtod->parity == mathematics::FunctionParity::Odd,
        "MathRegistry: angle conversion functions are odd linear maps");
    tests.expect(sqrt
        && sqrt->domainRule == mathematics::FunctionDomainRule::ComplexToComplex
        && sqrt->branchRule == mathematics::FunctionBranchRule::PrincipalSquareRoot,
        "MathRegistry: sqrt separates complex domain from principal branch semantics");
    tests.expect(abs && abs->parity == mathematics::FunctionParity::Even
        && abs->domainRule == mathematics::FunctionDomainRule::ComplexToReal,
        "MathRegistry: abs is an even complex-to-real magnitude function");
    tests.expect(sign && sign->parity == mathematics::FunctionParity::Odd
        && re && re->domainRule == mathematics::FunctionDomainRule::ComplexToReal
        && im && im->domainRule == mathematics::FunctionDomainRule::ComplexToReal
        && conj && conj->parity == mathematics::FunctionParity::Odd,
        "MathRegistry: sign/re/im/conj mathematical properties are registered");
    tests.expect(sin && sin->parity == mathematics::FunctionParity::Odd
        && sin->domainRule == mathematics::FunctionDomainRule::ComplexToComplexRealPreserving
        && sin->branchRule == mathematics::FunctionBranchRule::SingleValued,
        "MathRegistry: sin is single-valued, complex-capable, and real-preserving");
    tests.expect(cos && cos->parity == mathematics::FunctionParity::Even,
        "MathRegistry: cos is even");
    tests.expect(sin && sin->periodTurns && sin->periodTurns->toString() == "1",
        "MathRegistry: sin period is one full turn");
    tests.expect(sin && sin->inverseFunction == mathematics::FunctionId::Asin
        && !sin->realGloballyInjective
        && sin->realRangeRule == mathematics::RealRangeRule::ClosedMinusOneToOne,
        "MathRegistry: sin records principal inverse/range knowledge without pretending global injectivity");
    tests.expect(tan && tan->parity == mathematics::FunctionParity::Odd
        && tan->periodTurns && tan->periodTurns->toString() == "1/2"
        && tan->definednessRule == mathematics::FunctionDefinednessRule::CosNonZero,
        "MathRegistry: tan records both periodicity and its cosine pole condition");


    const auto* cot = mathematics.findFunction(mathematics::FunctionId::Cot);
    const auto* sec = mathematics.findFunction(mathematics::FunctionId::Sec);
    const auto* csc = mathematics.findFunction(mathematics::FunctionId::Csc);
    const auto* asin = mathematics.findFunction(mathematics::FunctionId::Asin);
    const auto* acos = mathematics.findFunction(mathematics::FunctionId::Acos);
    const auto* atan = mathematics.findFunction(mathematics::FunctionId::Atan);
    const auto* atan2 = mathematics.findFunction(mathematics::FunctionId::Atan2);
    tests.expect(cot && cot->parity == mathematics::FunctionParity::Odd
        && cot->periodTurns && cot->periodTurns->toString() == "1/2",
        "MathRegistry: cot is odd with a half-turn period");
    tests.expect(sec && sec->parity == mathematics::FunctionParity::Even
        && sec->periodTurns && sec->periodTurns->toString() == "1",
        "MathRegistry: sec is even with a full-turn period");
    tests.expect(csc && csc->parity == mathematics::FunctionParity::Odd
        && csc->periodTurns && csc->periodTurns->toString() == "1",
        "MathRegistry: csc is odd with a full-turn period");
    tests.expect(asin && asin->branchRule == mathematics::FunctionBranchRule::PrincipalArcSine,
        "MathRegistry: asin records its principal branch");
    tests.expect(acos && acos->branchRule == mathematics::FunctionBranchRule::PrincipalArcCosine,
        "MathRegistry: acos records its principal branch");
    tests.expect(atan && atan->branchRule == mathematics::FunctionBranchRule::PrincipalArcTangent,
        "MathRegistry: atan records its principal branch");
    tests.expect(atan2 && atan2->arity == 2
        && atan2->domainRule == mathematics::FunctionDomainRule::RealPairToReal
        && atan2->branchRule == mathematics::FunctionBranchRule::PrincipalAtan2,
        "MathRegistry: atan2 is a real two-argument principal-angle function");

    const auto* sinh = mathematics.findFunction(mathematics::FunctionId::Sinh);
    const auto* cosh = mathematics.findFunction(mathematics::FunctionId::Cosh);
    const auto* tanh = mathematics.findFunction(mathematics::FunctionId::Tanh);
    const auto* asinh = mathematics.findFunction(mathematics::FunctionId::Asinh);
    const auto* acosh = mathematics.findFunction(mathematics::FunctionId::Acosh);
    const auto* atanh = mathematics.findFunction(mathematics::FunctionId::Atanh);
    const auto* csch = mathematics.findFunction(mathematics::FunctionId::Csch);
    const auto* sech = mathematics.findFunction(mathematics::FunctionId::Sech);
    const auto* coth = mathematics.findFunction(mathematics::FunctionId::Coth);
    tests.expect(sinh && sinh->parity == mathematics::FunctionParity::Odd
        && sinh->domainRule == mathematics::FunctionDomainRule::ComplexToComplexRealPreserving,
        "MathRegistry: sinh is odd and real-preserving");
    tests.expect(cosh && cosh->parity == mathematics::FunctionParity::Even,
        "MathRegistry: cosh is even");
    tests.expect(tanh && tanh->parity == mathematics::FunctionParity::Odd,
        "MathRegistry: tanh is odd");
    tests.expect(tanh && tanh->inverseFunction == mathematics::FunctionId::Atanh
        && tanh->realGloballyInjective
        && tanh->realRangeRule == mathematics::RealRangeRule::OpenMinusOneToOne,
        "MathRegistry: tanh records its global real inverse and open unit range");
    tests.expect(asinh && asinh->branchRule == mathematics::FunctionBranchRule::PrincipalAreaHyperbolicSine,
        "MathRegistry: asinh records its principal branch");
    tests.expect(acosh && acosh->branchRule == mathematics::FunctionBranchRule::PrincipalAreaHyperbolicCosine,
        "MathRegistry: acosh records its principal branch");
    tests.expect(atanh && atanh->branchRule == mathematics::FunctionBranchRule::PrincipalAreaHyperbolicTangent
        && atanh->definednessRule == mathematics::FunctionDefinednessRule::OneMinusSquareNonZero,
        "MathRegistry: atanh records its principal branch and pole condition");
    tests.expect(csch && csch->parity == mathematics::FunctionParity::Odd
        && sech && sech->parity == mathematics::FunctionParity::Even
        && coth && coth->parity == mathematics::FunctionParity::Odd,
        "MathRegistry: reciprocal hyperbolic parity is registered");

    const auto* arg = mathematics.findFunction(mathematics::FunctionId::Arg);
    const auto* log = mathematics.findFunction(mathematics::FunctionId::Log);
    const auto* exp = mathematics.findFunction(mathematics::FunctionId::Exp);
    tests.expect(arg
        && arg->domainRule == mathematics::FunctionDomainRule::ComplexToReal
        && arg->branchRule == mathematics::FunctionBranchRule::PrincipalArgument,
        "MathRegistry: Arg records its principal range/branch semantics");
    tests.expect(log
        && log->domainRule == mathematics::FunctionDomainRule::ComplexToComplex
        && log->branchRule == mathematics::FunctionBranchRule::PrincipalLogarithm
        && log->definednessRule == mathematics::FunctionDefinednessRule::Logarithm
        && log->arity == 1 && log->maximumArity == 2
        && log->acceptsArity(1) && log->acceptsArity(2) && !log->acceptsArity(3),
        "MathRegistry: Log records principal branch and unary/binary arity semantics");
    tests.expect(exp
        && exp->domainRule == mathematics::FunctionDomainRule::ComplexToComplexRealPreserving
        && exp->branchRule == mathematics::FunctionBranchRule::SingleValued
        && exp->definednessRule == mathematics::FunctionDefinednessRule::Everywhere,
        "MathRegistry: Exp is single-valued, real-preserving, and entire");
    tests.expect(exp && exp->inverseFunction == mathematics::FunctionId::Log
        && exp->realGloballyInjective
        && exp->realMonotonicity == mathematics::RealMonotonicity::Increasing
        && exp->realRangeRule == mathematics::RealRangeRule::Positive,
        "MathRegistry: Exp records its global real inverse and positive range");
    tests.expect(log && log->inverseFunction == mathematics::FunctionId::Exp
        && log->realGloballyInjective
        && log->realRangeRule == mathematics::RealRangeRule::AllReal,
        "MathRegistry: Log records the real inverse relation without changing its principal complex branch");

    const auto* gamma = mathematics.findFunction(mathematics::FunctionId::Gamma);
    const auto* logGamma = mathematics.findFunction(mathematics::FunctionId::LogGamma);
    const auto* erf = mathematics.findFunction(mathematics::FunctionId::Erf);
    const auto* erfc = mathematics.findFunction(mathematics::FunctionId::Erfc);
    const auto* beta = mathematics.findFunction(mathematics::FunctionId::Beta);
    const auto* betaLog = mathematics.findFunction(mathematics::FunctionId::BetaLog);
    tests.expect(gamma
        && gamma->domainRule == mathematics::FunctionDomainRule::ComplexToComplex
        && gamma->definednessRule == mathematics::FunctionDefinednessRule::GammaPoles,
        "MathRegistry: Gamma records its infinite discrete pole set without weakening it");
    tests.expect(logGamma
        && logGamma->domainRule == mathematics::FunctionDomainRule::RealToReal
        && logGamma->definednessRule == mathematics::FunctionDefinednessRule::GammaPoles,
        "MathRegistry: lgamma is currently the real log-absolute-Gamma function");
    tests.expect(erf && erf->parity == mathematics::FunctionParity::Odd
        && erf->domainRule == mathematics::FunctionDomainRule::ComplexToComplexRealPreserving
        && erfc && erfc->domainRule == mathematics::FunctionDomainRule::ComplexToComplexRealPreserving,
        "MathRegistry: erf/erfc are entire complex functions with real-axis preservation");
    tests.expect(beta && beta->arity == 2 && beta->maximumArity == 2
        && beta->domainRule == mathematics::FunctionDomainRule::RealPairToReal
        && beta->definednessRule == mathematics::FunctionDefinednessRule::ArgumentsPositiveReal
        && betaLog && betaLog->definednessRule == mathematics::FunctionDefinednessRule::ArgumentsPositiveReal,
        "MathRegistry: current beta/betaln contract is positive-real and binary");

    const auto* power = mathematics.findFunction(mathematics::FunctionId::Power);
    tests.expect(power
        && power->arity == 2
        && power->domainRule == mathematics::FunctionDomainRule::ComplexToComplex
        && power->branchRule == mathematics::FunctionBranchRule::PrincipalPower,
        "MathRegistry: Power records principal Exp[w Log[z]] branch semantics");
}

} // namespace mmcal::tests
