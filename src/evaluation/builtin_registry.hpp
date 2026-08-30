#pragma once

#include "expression/symbol.hpp"
#include "symbols/symbol_id.hpp"
#include "symbols/symbol_table.hpp"

#include <cstddef>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace mmcal::expression {
class Expr;
}

namespace mmcal::evaluation {

// 評価器内の分岐を文字列比較から分離するための組み込み函数識別子。
enum class BuiltinId {
    Add,
    Subtract,
    Multiply,
    Divide,
    Power,
    Negate,
    Factorial,
    Derivative,
    SymbolicIntegral,
    Limit,
    Floor,
    Ceil,
    Trunc,
    Round,
    Frac,
    BitAnd,
    BitOr,
    BitXor,
    BitNot,
    BitShiftLeft,
    BitShiftRight,
    BitLength,
    BitCount,
    BitGet,
    Gcd,
    Lcm,
    Mod,
    Rem,
    Quotient,
    IsPrime,
    NextPrime,
    PreviousPrime,
    FactorInteger,
    Totient,
    Permutation,
    Combination,
    Fibonacci,
    DiscreteFourierTransform,
    FastFourierTransform,
    InverseFourierTransform,
    Convolution,
    Transpose,
    ConjugateTranspose,
    MatrixAdd,
    MatrixMultiply,
    Determinant,
    Inverse,
    Rref,
    Rank,
    SolveLinear,
    NullSpace,
    LuDecomposition,
    QrDecomposition,
    SingularValueDecomposition,
    ConditionNumber,
    LeastSquares,
    PseudoInverse,
    Eigenvalues,
    Eigenvectors,
    Eigensystem,
    NumericDerivative,
    NumericIntegral,
    Cbrt,
    Hypot,
    Fma,
    Clamp,
    Proj,
    Cis,
    Polar,
    NextPow2,
    DegreeToRadian,
    DegreeToGradian,
    RadianToDegree,
    RadianToGradian,
    GradianToDegree,
    GradianToRadian,
    Sum,
    Product,
    Map,
    Range,
    Table,
    Min,
    Max,
    Mean,
    Median,
    Mode,
    Quantile,
    Percentile,
    VariancePopulation,
    VarianceSample,
    StddevPopulation,
    StddevSample,
    GeometricMean,
    HarmonicMean,
    Rms,
    MedianAbsoluteDeviation,
    MeanAbsoluteDeviation,
    Skewness,
    KurtosisPopulation,
    KurtosisSample,
    CoefficientVariation,
    StandardError,
    ZScore,
    Iqr,
    TrimMean,
    WinsorMean,
    Winsorized,
    Covariance,
    Correlation,
    SpearmanCorrelation,
    PercentRank,
    Dimensions,
    ArrayRank,
    Length,
    ArrayGet,
    Reshape,
    Identity,
    Zeros,
    Trace,
    Rows,
    Cols,
    Diag,
    VectorAdd,
    VectorSubtract,
    VectorScale,
    VectorCross,
    VectorNorm,
    VectorManhattan,
    VectorEuclidean,
    VectorNormalize,
    VectorProject,
    VectorAngle,
    VectorReflect,
    VectorReflectAxis,
    VectorSum,
    Expm1,
    Log1p,
    Sinc,
    Cosc,
    Tanc,
    Sinhc,
    Tanhc,
    Expc,
    Log2,
    Log10,
    Gamma,
    LogGamma,
    LambertW,
    Erf,
    Erfc,
    FresnelC,
    FresnelS,
    Hypergeometric1F1,
    Hypergeometric2F1,
    EllipticF,
    EllipticE,
    EllipticPi,
    ExponentialIntegralEi,
    SineIntegralSi,
    CosineIntegralCi,
    LogarithmicIntegralLi,
    Polylog,
    Beta,
    BetaLog,
    Zeta,
    Digamma,
    Trigamma,
    IncompleteBeta,
    GeneralizedBinomial,
    FallingFactorial,
    RisingFactorial,
    RandSeed,
    Rand,
    RandInt,
    Choice,
    RandN,
    Sqrt,
    Abs,
    Sign,
    Re,
    Im,
    Conj,
    Sin,
    Cos,
    Tan,
    Cot,
    Sec,
    Csc,
    Asin,
    Acos,
    Atan,
    Atan2,
    Sinh,
    Cosh,
    Tanh,
    Asinh,
    Acosh,
    Atanh,
    Csch,
    Sech,
    Coth,
    Arg,
    Log,
    Exp,
    NumericalApproximation,
    Precision,
    Accuracy,
    Explain,
    Rationalize,
    Root,
    Simplify,
    FullSimplify,
    Expand,
    Factor,
    Collect,
    Solve,
    GroebnerBasis,
    PolynomialReduce,
    Cases,
    CaseBranch,
    Set,
    SetDelayed,
    Less,
    LessEqual,
    Greater,
    GreaterEqual,
    Equal,
    NotEqual,
    LogicalAnd,
    Element,
    If,
    History,
    InputHistory,
    OutputHistory,
    Exit,
    Clear,
    Definitions,
    Undefine,
    AngleMode,
    UnitApplied
};

// 特殊形式を含め、各引数を評価する時期を定義する。
enum class ArgumentEvaluation {
    All,
    HoldFirst,
    HoldFirstTwo,
    HoldFirstAndIteratorSpec,
    HoldFirstAndTableIteratorSpec,
    HoldAll
};

struct BuiltinDefinition final {
    static constexpr std::size_t unlimited = std::numeric_limits<std::size_t>::max();

    expression::Symbol symbol;
    BuiltinId id = BuiltinId::Add;
    std::size_t minimumArguments = 0;
    std::size_t maximumArguments = unlimited;
    ArgumentEvaluation argumentEvaluation = ArgumentEvaluation::All;
    bool sourceCallable = false;

    [[nodiscard]] bool acceptsArity(std::size_t count) const noexcept;
    [[nodiscard]] std::string_view name() const noexcept { return symbol.view(); }
};

// 組み込み函数のSymbol、引数個数、評価属性を一元管理する登録表。
class BuiltinRegistry final {
public:
    explicit BuiltinRegistry(symbols::SymbolTable& symbolTable);

    [[nodiscard]] static BuiltinRegistry defaults(symbols::SymbolTable& symbolTable);

    void add(
        std::string_view name,
        BuiltinId id,
        std::size_t minimumArguments,
        std::size_t maximumArguments,
        ArgumentEvaluation argumentEvaluation = ArgumentEvaluation::All,
        bool sourceCallable = false);
    // 既存BuiltinIdへ別名を割り当てる。arity/evaluationはcanonical定義を共有する。
    void addAlias(
        std::string_view name,
        BuiltinId targetId,
        bool sourceCallable = true);
    [[nodiscard]] const BuiltinDefinition* find(const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] const BuiltinDefinition* find(std::string_view name) const noexcept;
    [[nodiscard]] bool contains(const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] bool contains(std::string_view name) const noexcept;
    // AST側でBuiltinIdを直接比較する共通query。各subsystemのisHead重複を避ける。
    [[nodiscard]] bool isCallTo(const expression::Expr& expression, BuiltinId id) const noexcept;
    [[nodiscard]] const expression::Symbol& symbol(BuiltinId id) const;
    [[nodiscard]] std::unordered_set<std::string> sourceFunctionNames() const;
    [[nodiscard]] std::size_t size() const noexcept;

private:
    symbols::SymbolTable& symbolTable_;
    std::unordered_map<symbols::SymbolId, BuiltinDefinition, symbols::SymbolIdHash> definitions_;
    std::unordered_map<BuiltinId, expression::Symbol> symbolsById_;
};

// 単体利用向けの既定表。通常のKernelSessionは専用表を所有する。
[[nodiscard]] const BuiltinRegistry& defaultBuiltinRegistry();

} // namespace mmcal::evaluation
