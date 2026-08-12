#pragma once

#include "evaluation/builtin_registry.hpp"
#include "math_ids.hpp"
#include "expression/symbol.hpp"
#include "numeric/rational.hpp"
#include "symbols/symbol_id.hpp"
#include "symbols/symbol_table.hpp"

#include <cstddef>
#include <optional>
#include <string_view>
#include <unordered_map>

namespace mmcal::mathematics {

// algebraic / transcendental は相互排他的なのでboolを二つ持たない。
enum class ArithmeticClass {
    Unknown,
    Algebraic,
    Transcendental
};

struct ConstantProperties final {
    bool exact = true;
    bool real = false;
    bool positive = false;
    bool irrational = false;
    ArithmeticClass arithmeticClass = ArithmeticClass::Unknown;
};

struct ConstantDefinition final {
    ConstantId id = ConstantId::Pi;
    expression::Symbol symbol;
    ConstantProperties properties;
};

enum class FunctionParity {
    Neither,
    Odd,
    Even
};

enum class RealMonotonicity {
    Unknown,
    Increasing,
    Decreasing
};

enum class RealRangeRule {
    Unknown,
    AllReal,
    Positive,
    NonNegative,
    OpenMinusOneToOne,
    ClosedMinusOneToOne,
    OneToInfinity
};

// 定義域内で函数値が0を取り得るかという知識。
// definednessとは別であり、NeverZeroは「値が存在する点では0にならない」ことだけを表す。
enum class FunctionZeroRule {
    Unknown,
    NeverZero
};

// 函数が実数/複素数domainをどう写すかという数学的性質。
// branchの選び方とは独立にする。sin/cosは複素全体で定義され実数入力を実数へ写す一方、sqrtは複素全体へ延長できるが負実数入力を非実数へ写し得る。
enum class FunctionDomainRule {
    Unknown,
    RealToReal,
    ComplexToComplex,
    ComplexToComplexRealPreserving,
    ComplexToReal,
    RealPairToReal
};

// 多価函数を一価函数として扱う際のbranch規約。
// sqrtのprincipal branchは Re(sqrt(z)) >= 0、Re=0ならIm>=0 を採る。
// 将来Log/一般Powerを追加するときも、この「domain」と「branch」を混ぜない。
enum class FunctionBranchRule {
    SingleValued,
    PrincipalSquareRoot,
    PrincipalArgument,
    PrincipalLogarithm,
    PrincipalPower,
    PrincipalArcSine,
    PrincipalArcCosine,
    PrincipalArcTangent,
    PrincipalAtan2,
    PrincipalAreaHyperbolicSine,
    PrincipalAreaHyperbolicCosine,
    PrincipalAreaHyperbolicTangent,
    PrincipalHypergeometric2F1,
    PrincipalElliptic
};

// 函数が有限入力で値を持つために必要な追加条件。
// branch規約や値域とは分離し、Solver/Simplifierが同じdefinedness知識を共有する。
enum class FunctionDefinednessRule {
    Everywhere,
    ArgumentReal,
    ArgumentsReal,
    ArgumentNonZero,
    Logarithm,
    OnePlusArgumentNonZero,
    SinNonZero,
    CosNonZero,
    SinhNonZero,
    CoshNonZero,
    OnePlusSquareNonZero,
    OneMinusSquareNonZero,
    PrincipalPower,
    RealPairNotBothZero,
    ArgumentsPositiveReal,
    GammaPoles,
    Hypergeometric1F1Poles,
    Hypergeometric2F1Poles,
    EllipticPrincipal
};

struct FunctionDefinition final {
    FunctionId id = FunctionId::Sin;
    expression::Symbol symbol;
    // arityは最小引数数。maximumArity==arityなら固定arity。
    std::size_t arity = 1;
    std::size_t maximumArity = 1;
    FunctionParity parity = FunctionParity::Neither;
    FunctionDomainRule domainRule = FunctionDomainRule::Unknown;
    FunctionBranchRule branchRule = FunctionBranchRule::SingleValued;
    FunctionDefinednessRule definednessRule = FunctionDefinednessRule::Everywhere;

    // 周期は単位依存の「360」「2 Pi」ではなく、1回転=1という無次元量で記録する。
    // これによりDegree/Radian/Gradianのどれを入力しても同じ数学知識を再利用できる。
    std::optional<numeric::Rational> periodTurns;

    // 実軸上での逆函数知識。principal inverseと全解集合を混同しないため、周期函数ではrealGloballyInjective=falseのままinverseFunctionだけを保持する。
    std::optional<FunctionId> inverseFunction;
    bool realGloballyInjective = false;
    RealMonotonicity realMonotonicity = RealMonotonicity::Unknown;
    RealRangeRule realRangeRule = RealRangeRule::Unknown;
    FunctionZeroRule zeroRule = FunctionZeroRule::Unknown;

    [[nodiscard]] bool acceptsArity(std::size_t count) const noexcept {
        return count >= arity && count <= maximumArity;
    }
};

// SymbolRegistry/BuiltinRegistryが「言語上の名前」を管理するのに対し、MathRegistryはPiやsinが数学的に何者であるかを管理する。
class MathRegistry final {
public:
    MathRegistry(
        symbols::SymbolTable& symbolTable,
        const evaluation::BuiltinRegistry& builtins);

    [[nodiscard]] static MathRegistry defaults(
        symbols::SymbolTable& symbolTable,
        const evaluation::BuiltinRegistry& builtins);

    [[nodiscard]] const ConstantDefinition* findConstant(
        const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] const ConstantDefinition* findConstant(ConstantId id) const noexcept;
    [[nodiscard]] const FunctionDefinition* findFunction(
        const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] const FunctionDefinition* findFunction(FunctionId id) const noexcept;

    [[nodiscard]] std::size_t constantCount() const noexcept;
    [[nodiscard]] std::size_t functionCount() const noexcept;

private:
    symbols::SymbolTable& symbolTable_;
    const evaluation::BuiltinRegistry& builtins_;
    std::unordered_map<symbols::SymbolId, ConstantDefinition, symbols::SymbolIdHash> constants_;
    std::unordered_map<ConstantId, expression::Symbol> constantSymbols_;
    std::unordered_map<symbols::SymbolId, FunctionDefinition, symbols::SymbolIdHash> functions_;
    std::unordered_map<FunctionId, expression::Symbol> functionSymbols_;

    void addConstant(
        std::string_view name,
        ConstantId id,
        ConstantProperties properties);
    void addFunction(
        evaluation::BuiltinId builtinId,
        FunctionId id,
        FunctionParity parity,
        FunctionDomainRule domainRule,
        FunctionBranchRule branchRule,
        std::optional<numeric::Rational> periodTurns,
        FunctionDefinednessRule definednessRule = FunctionDefinednessRule::Everywhere,
        std::size_t arity = 1,
        std::size_t maximumArity = 0,
        std::optional<FunctionId> inverseFunction = std::nullopt,
        bool realGloballyInjective = false,
        RealMonotonicity realMonotonicity = RealMonotonicity::Unknown,
        RealRangeRule realRangeRule = RealRangeRule::Unknown);
    void setRealInverseKnowledge(
        FunctionId id,
        FunctionId inverse,
        bool globallyInjective,
        RealMonotonicity monotonicity,
        RealRangeRule rangeRule);
    void setZeroKnowledge(FunctionId id, FunctionZeroRule zeroRule);
};

[[nodiscard]] const MathRegistry& defaultMathRegistry();

} // namespace mmcal::mathematics
