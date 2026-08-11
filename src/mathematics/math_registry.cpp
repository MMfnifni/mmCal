// 函数のdomain・逆函数・周期などの数学metadata
#include "math_registry.hpp"

#include "numeric/big_int.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::mathematics {
namespace {

[[nodiscard]] numeric::Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return numeric::Rational{numeric::BigInt{numerator}, numeric::BigInt{denominator}};
}

} // namespace

MathRegistry::MathRegistry(
    symbols::SymbolTable& symbolTable,
    const evaluation::BuiltinRegistry& builtins)
    : symbolTable_(symbolTable), builtins_(builtins) {}

MathRegistry MathRegistry::defaults(
    symbols::SymbolTable& symbolTable,
    const evaluation::BuiltinRegistry& builtins) {
    MathRegistry registry{symbolTable, builtins};

    registry.addConstant(
        "Pi",
        ConstantId::Pi,
        ConstantProperties{true, true, true, true, ArithmeticClass::Transcendental});
    registry.addConstant(
        "E",
        ConstantId::E,
        ConstantProperties{true, true, true, true, ArithmeticClass::Transcendental});
    registry.addConstant(
        "Phi",
        ConstantId::Phi,
        ConstantProperties{true, true, true, true, ArithmeticClass::Algebraic});

    registry.addFunction(
        evaluation::BuiltinId::Cbrt,
        FunctionId::Cbrt,
        FunctionParity::Odd,
        FunctionDomainRule::RealToReal,
        FunctionBranchRule::SingleValued,
        std::nullopt,
        FunctionDefinednessRule::ArgumentReal);
    registry.addFunction(
        evaluation::BuiltinId::Hypot,
        FunctionId::Hypot,
        FunctionParity::Neither,
        FunctionDomainRule::RealPairToReal,
        FunctionBranchRule::SingleValued,
        std::nullopt,
        FunctionDefinednessRule::ArgumentsReal,
        2);
    registry.addFunction(
        evaluation::BuiltinId::Cis,
        FunctionId::Cis,
        FunctionParity::Neither,
        FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::SingleValued,
        rational(1));
    registry.addFunction(
        evaluation::BuiltinId::Polar,
        FunctionId::Polar,
        FunctionParity::Neither,
        FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::SingleValued,
        std::nullopt,
        FunctionDefinednessRule::Everywhere,
        2);

    const struct {
        evaluation::BuiltinId builtin;
        FunctionId function;
    } angleConversions[] = {
        {evaluation::BuiltinId::DegreeToRadian, FunctionId::DegreeToRadian},
        {evaluation::BuiltinId::DegreeToGradian, FunctionId::DegreeToGradian},
        {evaluation::BuiltinId::RadianToDegree, FunctionId::RadianToDegree},
        {evaluation::BuiltinId::RadianToGradian, FunctionId::RadianToGradian},
        {evaluation::BuiltinId::GradianToDegree, FunctionId::GradianToDegree},
        {evaluation::BuiltinId::GradianToRadian, FunctionId::GradianToRadian}
    };
    for (const auto& conversion : angleConversions)
        registry.addFunction(
            conversion.builtin, conversion.function, FunctionParity::Odd,
            FunctionDomainRule::ComplexToComplexRealPreserving,
            FunctionBranchRule::SingleValued, std::nullopt);

    registry.addFunction(evaluation::BuiltinId::Expm1, FunctionId::Expm1,
        FunctionParity::Neither, FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued, std::nullopt);
    registry.addFunction(evaluation::BuiltinId::Log1p, FunctionId::Log1p,
        FunctionParity::Neither, FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::PrincipalLogarithm, std::nullopt,
        FunctionDefinednessRule::OnePlusArgumentNonZero);
    registry.addFunction(evaluation::BuiltinId::Sinc, FunctionId::Sinc,
        FunctionParity::Even, FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued, std::nullopt);
    registry.addFunction(evaluation::BuiltinId::Cosc, FunctionId::Cosc,
        FunctionParity::Odd, FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued, std::nullopt);
    registry.addFunction(evaluation::BuiltinId::Tanc, FunctionId::Tanc,
        FunctionParity::Even, FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued, std::nullopt, FunctionDefinednessRule::CosNonZero);
    registry.addFunction(evaluation::BuiltinId::Sinhc, FunctionId::Sinhc,
        FunctionParity::Even, FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued, std::nullopt);
    registry.addFunction(evaluation::BuiltinId::Tanhc, FunctionId::Tanhc,
        FunctionParity::Even, FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued, std::nullopt, FunctionDefinednessRule::CoshNonZero);
    registry.addFunction(evaluation::BuiltinId::Expc, FunctionId::Expc,
        FunctionParity::Neither, FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued, std::nullopt);

    // 特殊函数。Gammaのpole集合 {0,-1,-2,...} は現Predicateでは1個の有限条件へ落とせないため、
    // 専用definedness ruleとして保持してSolver側で不完全な条件を捏造しない。
    registry.addFunction(evaluation::BuiltinId::Gamma, FunctionId::Gamma,
        FunctionParity::Neither, FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::SingleValued, std::nullopt, FunctionDefinednessRule::GammaPoles);
    registry.addFunction(evaluation::BuiltinId::LogGamma, FunctionId::LogGamma,
        FunctionParity::Neither, FunctionDomainRule::RealToReal,
        FunctionBranchRule::SingleValued, std::nullopt, FunctionDefinednessRule::GammaPoles);
    registry.addFunction(evaluation::BuiltinId::Erf, FunctionId::Erf,
        FunctionParity::Odd, FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued, std::nullopt);
    registry.addFunction(evaluation::BuiltinId::Erfc, FunctionId::Erfc,
        FunctionParity::Neither, FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued, std::nullopt);
    registry.addFunction(evaluation::BuiltinId::Beta, FunctionId::Beta,
        FunctionParity::Neither, FunctionDomainRule::RealPairToReal,
        FunctionBranchRule::SingleValued, std::nullopt,
        FunctionDefinednessRule::ArgumentsPositiveReal, 2);
    registry.addFunction(evaluation::BuiltinId::BetaLog, FunctionId::BetaLog,
        FunctionParity::Neither, FunctionDomainRule::RealPairToReal,
        FunctionBranchRule::SingleValued, std::nullopt,
        FunctionDefinednessRule::ArgumentsPositiveReal, 2);

    registry.addFunction(
        evaluation::BuiltinId::Sqrt,
        FunctionId::Sqrt,
        FunctionParity::Neither,
        FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::PrincipalSquareRoot,
        std::nullopt);

    registry.addFunction(
        evaluation::BuiltinId::Abs,
        FunctionId::Abs,
        FunctionParity::Even,
        FunctionDomainRule::ComplexToReal,
        FunctionBranchRule::SingleValued,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Sign,
        FunctionId::Sign,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Re,
        FunctionId::Re,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToReal,
        FunctionBranchRule::SingleValued,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Im,
        FunctionId::Im,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToReal,
        FunctionBranchRule::SingleValued,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Conj,
        FunctionId::Conj,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        std::nullopt);

    const numeric::Rational fullTurn = rational(1);
    registry.addFunction(
        evaluation::BuiltinId::Sin,
        FunctionId::Sin,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        fullTurn);
    registry.addFunction(
        evaluation::BuiltinId::Cos,
        FunctionId::Cos,
        FunctionParity::Even,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        fullTurn);
    registry.addFunction(
        evaluation::BuiltinId::Tan,
        FunctionId::Tan,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        rational(1, 2),
        FunctionDefinednessRule::CosNonZero);

    registry.addFunction(
        evaluation::BuiltinId::Cot,
        FunctionId::Cot,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        rational(1, 2),
        FunctionDefinednessRule::SinNonZero);
    registry.addFunction(
        evaluation::BuiltinId::Sec,
        FunctionId::Sec,
        FunctionParity::Even,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        fullTurn,
        FunctionDefinednessRule::CosNonZero);
    registry.addFunction(
        evaluation::BuiltinId::Csc,
        FunctionId::Csc,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        fullTurn,
        FunctionDefinednessRule::SinNonZero);

    registry.addFunction(
        evaluation::BuiltinId::Asin,
        FunctionId::Asin,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::PrincipalArcSine,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Acos,
        FunctionId::Acos,
        FunctionParity::Neither,
        FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::PrincipalArcCosine,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Atan,
        FunctionId::Atan,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::PrincipalArcTangent,
        std::nullopt,
        FunctionDefinednessRule::OnePlusSquareNonZero);
    registry.addFunction(
        evaluation::BuiltinId::Atan2,
        FunctionId::Atan2,
        FunctionParity::Neither,
        FunctionDomainRule::RealPairToReal,
        FunctionBranchRule::PrincipalAtan2,
        std::nullopt,
        FunctionDefinednessRule::RealPairNotBothZero,
        2);

    registry.addFunction(
        evaluation::BuiltinId::Sinh,
        FunctionId::Sinh,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Cosh,
        FunctionId::Cosh,
        FunctionParity::Even,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Tanh,
        FunctionId::Tanh,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        std::nullopt,
        FunctionDefinednessRule::CoshNonZero);
    registry.addFunction(
        evaluation::BuiltinId::Asinh,
        FunctionId::Asinh,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::PrincipalAreaHyperbolicSine,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Acosh,
        FunctionId::Acosh,
        FunctionParity::Neither,
        FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::PrincipalAreaHyperbolicCosine,
        std::nullopt);
    registry.addFunction(
        evaluation::BuiltinId::Atanh,
        FunctionId::Atanh,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::PrincipalAreaHyperbolicTangent,
        std::nullopt,
        FunctionDefinednessRule::OneMinusSquareNonZero);
    registry.addFunction(
        evaluation::BuiltinId::Csch,
        FunctionId::Csch,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        std::nullopt,
        FunctionDefinednessRule::SinhNonZero);
    registry.addFunction(
        evaluation::BuiltinId::Sech,
        FunctionId::Sech,
        FunctionParity::Even,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        std::nullopt,
        FunctionDefinednessRule::CoshNonZero);
    registry.addFunction(
        evaluation::BuiltinId::Coth,
        FunctionId::Coth,
        FunctionParity::Odd,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        std::nullopt,
        FunctionDefinednessRule::SinhNonZero);

    // principal Arg は非零複素数を (-Pi, Pi] の実数へ写す。
    // 負実軸上では +Pi を採り、0では未定義。
    registry.addFunction(
        evaluation::BuiltinId::Arg,
        FunctionId::Arg,
        FunctionParity::Neither,
        FunctionDomainRule::ComplexToReal,
        FunctionBranchRule::PrincipalArgument,
        std::nullopt,
        FunctionDefinednessRule::ArgumentNonZero);

    // principal Log(z) = ln|z| + I Arg(z)。0では未定義で、branch cut は負実軸。Argと同じbranch規約を共有する。
    registry.addFunction(
        evaluation::BuiltinId::Log,
        FunctionId::Log,
        FunctionParity::Neither,
        FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::PrincipalLogarithm,
        std::nullopt,
        FunctionDefinednessRule::Logarithm,
        1,
        2);

    // Exp は複素平面全体で一価・正則。実数入力なら正の実数へ写す。
    registry.addFunction(
        evaluation::BuiltinId::Exp,
        FunctionId::Exp,
        FunctionParity::Neither,
        FunctionDomainRule::ComplexToComplexRealPreserving,
        FunctionBranchRule::SingleValued,
        std::nullopt);

    // principal Power(z,w) は非零baseについて
    //   Exp(w * principal Log(z))
    // と定義する。方程式 z^n=a の全解集合とは別概念であり、Power自身は一価。
    registry.addFunction(
        evaluation::BuiltinId::Power,
        FunctionId::Power,
        FunctionParity::Neither,
        FunctionDomainRule::ComplexToComplex,
        FunctionBranchRule::PrincipalPower,
        std::nullopt,
        FunctionDefinednessRule::PrincipalPower,
        2);

    // principal inverseと実軸上の全域単射性を分離して保持する。
    // sin/cos/tanはinverseを持つが周期函数なのでglobal injectiveにはしない。
    registry.setRealInverseKnowledge(FunctionId::Exp, FunctionId::Log, true,
        RealMonotonicity::Increasing, RealRangeRule::Positive);
    registry.setRealInverseKnowledge(FunctionId::Log, FunctionId::Exp, true,
        RealMonotonicity::Increasing, RealRangeRule::AllReal);
    registry.setRealInverseKnowledge(FunctionId::Sinh, FunctionId::Asinh, true,
        RealMonotonicity::Increasing, RealRangeRule::AllReal);
    registry.setRealInverseKnowledge(FunctionId::Asinh, FunctionId::Sinh, true,
        RealMonotonicity::Increasing, RealRangeRule::AllReal);
    registry.setRealInverseKnowledge(FunctionId::Tanh, FunctionId::Atanh, true,
        RealMonotonicity::Increasing, RealRangeRule::OpenMinusOneToOne);
    registry.setRealInverseKnowledge(FunctionId::Atanh, FunctionId::Tanh, true,
        RealMonotonicity::Increasing, RealRangeRule::AllReal);
    registry.setRealInverseKnowledge(FunctionId::Sin, FunctionId::Asin, false,
        RealMonotonicity::Unknown, RealRangeRule::ClosedMinusOneToOne);
    registry.setRealInverseKnowledge(FunctionId::Cos, FunctionId::Acos, false,
        RealMonotonicity::Unknown, RealRangeRule::ClosedMinusOneToOne);
    registry.setRealInverseKnowledge(FunctionId::Tan, FunctionId::Atan, false,
        RealMonotonicity::Unknown, RealRangeRule::AllReal);

    return registry;
}

void MathRegistry::addConstant(
    std::string_view name,
    ConstantId id,
    ConstantProperties properties) {
    const expression::Symbol symbol = symbolTable_.intern(name);
    const auto [iterator, inserted] = constants_.emplace(
        symbol.id(),
        ConstantDefinition{id, symbol, properties});
    static_cast<void>(iterator);
    if (!inserted)
        throw std::invalid_argument("Mathematical constant is already registered: " + symbol.name());

    const auto [idIterator, idInserted] = constantSymbols_.emplace(id, symbol);
    static_cast<void>(idIterator);
    if (!idInserted)
        throw std::invalid_argument("Mathematical constant ID is already registered");
}

void MathRegistry::addFunction(
    evaluation::BuiltinId builtinId,
    FunctionId id,
    FunctionParity parity,
    FunctionDomainRule domainRule,
    FunctionBranchRule branchRule,
    std::optional<numeric::Rational> periodTurns,
    FunctionDefinednessRule definednessRule,
    std::size_t arity,
    std::size_t maximumArity,
    std::optional<FunctionId> inverseFunction,
    bool realGloballyInjective,
    RealMonotonicity realMonotonicity,
    RealRangeRule realRangeRule) {
    const expression::Symbol symbol = builtins_.symbol(builtinId);
    if (maximumArity == 0)
        maximumArity = arity;
    if (maximumArity < arity)
        throw std::invalid_argument("Mathematical function maximum arity is smaller than minimum arity");
    const auto [iterator, inserted] = functions_.emplace(
        symbol.id(),
        FunctionDefinition{
            id, symbol, arity, maximumArity, parity, domainRule, branchRule, definednessRule,
            std::move(periodTurns), inverseFunction, realGloballyInjective,
            realMonotonicity, realRangeRule});
    static_cast<void>(iterator);
    if (!inserted)
        throw std::invalid_argument("Mathematical function is already registered: " + symbol.name());

    const auto [idIterator, idInserted] = functionSymbols_.emplace(id, symbol);
    static_cast<void>(idIterator);
    if (!idInserted)
        throw std::invalid_argument("Mathematical function ID is already registered");
}


void MathRegistry::setRealInverseKnowledge(
    FunctionId id,
    FunctionId inverse,
    bool globallyInjective,
    RealMonotonicity monotonicity,
    RealRangeRule rangeRule) {
    const auto symbolIterator = functionSymbols_.find(id);
    if (symbolIterator == functionSymbols_.end())
        throw std::invalid_argument("Mathematical function ID is not registered");
    auto iterator = functions_.find(symbolIterator->second.id());
    if (iterator == functions_.end())
        throw std::invalid_argument("Mathematical function symbol is not registered");
    iterator->second.inverseFunction = inverse;
    iterator->second.realGloballyInjective = globallyInjective;
    iterator->second.realMonotonicity = monotonicity;
    iterator->second.realRangeRule = rangeRule;
}

const ConstantDefinition* MathRegistry::findConstant(
    const expression::Symbol& symbol) const noexcept {
    const auto iterator = constants_.find(symbol.id());
    return iterator == constants_.end() ? nullptr : &iterator->second;
}

const ConstantDefinition* MathRegistry::findConstant(ConstantId id) const noexcept {
    const auto symbolIterator = constantSymbols_.find(id);
    if (symbolIterator == constantSymbols_.end())
        return nullptr;
    return findConstant(symbolIterator->second);
}

const FunctionDefinition* MathRegistry::findFunction(
    const expression::Symbol& symbol) const noexcept {
    const auto iterator = functions_.find(symbol.id());
    return iterator == functions_.end() ? nullptr : &iterator->second;
}

const FunctionDefinition* MathRegistry::findFunction(FunctionId id) const noexcept {
    const auto symbolIterator = functionSymbols_.find(id);
    if (symbolIterator == functionSymbols_.end())
        return nullptr;
    return findFunction(symbolIterator->second);
}

std::size_t MathRegistry::constantCount() const noexcept {
    return constants_.size();
}

std::size_t MathRegistry::functionCount() const noexcept {
    return functions_.size();
}

const MathRegistry& defaultMathRegistry() {
    static const MathRegistry registry = MathRegistry::defaults(
        symbols::defaultSymbolTable(),
        evaluation::defaultBuiltinRegistry());
    return registry;
}

} // namespace mmcal::mathematics
