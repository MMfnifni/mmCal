// 組込み函数登録と属性
#include "builtin_registry.hpp"

#include "builtins/names.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::evaluation {
namespace {

struct BuiltinSpec final {
    std::string_view name;
    BuiltinId id;
    std::size_t minimumArguments;
    std::size_t maximumArguments;
    ArgumentEvaluation argumentEvaluation = ArgumentEvaluation::All;
    bool sourceCallable = false;
};

struct BuiltinAliasSpec final {
    std::string_view name;
    BuiltinId targetId;
    bool sourceCallable = true;
};

[[nodiscard]] constexpr BuiltinSpec fixed(
    std::string_view name,
    BuiltinId id,
    std::size_t arity,
    ArgumentEvaluation evaluation = ArgumentEvaluation::All,
    bool sourceCallable = false) noexcept {
    return BuiltinSpec{name, id, arity, arity, evaluation, sourceCallable};
}

[[nodiscard]] constexpr BuiltinSpec range(
    std::string_view name,
    BuiltinId id,
    std::size_t minimumArguments,
    std::size_t maximumArguments,
    ArgumentEvaluation evaluation = ArgumentEvaluation::All,
    bool sourceCallable = false) noexcept {
    return BuiltinSpec{
        name, id, minimumArguments, maximumArguments, evaluation, sourceCallable};
}

[[nodiscard]] constexpr BuiltinSpec variadic(
    std::string_view name,
    BuiltinId id,
    std::size_t minimumArguments,
    ArgumentEvaluation evaluation = ArgumentEvaluation::All,
    bool sourceCallable = false) noexcept {
    return range(
        name, id, minimumArguments, BuiltinDefinition::unlimited, evaluation, sourceCallable);
}

constexpr BuiltinSpec defaultBuiltins[] = {
    variadic(builtins::names::add, BuiltinId::Add, 0),
    fixed(builtins::names::subtract, BuiltinId::Subtract, 2),
    variadic(builtins::names::multiply, BuiltinId::Multiply, 0),
    fixed(builtins::names::divide, BuiltinId::Divide, 2),
    fixed(builtins::names::power, BuiltinId::Power, 2),
    fixed(builtins::names::negate, BuiltinId::Negate, 1),
    fixed(builtins::names::factorial, BuiltinId::Factorial, 1),

    // 記号演算・離散数学。
    variadic(builtins::names::derivative, BuiltinId::Derivative, 2, ArgumentEvaluation::HoldAll, true),
    range(builtins::names::symbolicIntegral, BuiltinId::SymbolicIntegral, 2, 3,
        ArgumentEvaluation::HoldFirstAndIteratorSpec, true),
    range(builtins::names::limit, BuiltinId::Limit, 3, 4,
        ArgumentEvaluation::HoldFirstTwo, true),
    fixed(builtins::names::floor, BuiltinId::Floor, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::ceil, BuiltinId::Ceil, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::trunc, BuiltinId::Trunc, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::round, BuiltinId::Round, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::frac, BuiltinId::Frac, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::gcd, BuiltinId::Gcd, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::lcm, BuiltinId::Lcm, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::mod, BuiltinId::Mod, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::rem, BuiltinId::Rem, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::quotient, BuiltinId::Quotient, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::permutation, BuiltinId::Permutation, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::combination, BuiltinId::Combination, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::fibonacci, BuiltinId::Fibonacci, 1, ArgumentEvaluation::All, true),

    // exact/symbolic signal processing. FFT is radix-2 with exact DFT fallback.
    fixed(builtins::names::dft, BuiltinId::DiscreteFourierTransform, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::fft, BuiltinId::FastFourierTransform, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::ifft, BuiltinId::InverseFourierTransform, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::convolve, BuiltinId::Convolution, 2, ArgumentEvaluation::All, true),

    // exact-first線形代数。
    fixed(builtins::names::transpose, BuiltinId::Transpose, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::matrixAdd, BuiltinId::MatrixAdd, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::matrixMultiply, BuiltinId::MatrixMultiply, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::determinant, BuiltinId::Determinant, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::inverse, BuiltinId::Inverse, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::rref, BuiltinId::Rref, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::rank, BuiltinId::Rank, 1, ArgumentEvaluation::All, true),

    // 数値解析。binder型は被積分式とiterator変数だけを保持する。
    range(builtins::names::numericDerivative, BuiltinId::NumericDerivative, 3, 4,
        ArgumentEvaluation::HoldFirstTwo, true),
    range(builtins::names::numericIntegral, BuiltinId::NumericIntegral, 2, 3,
        ArgumentEvaluation::HoldFirstAndIteratorSpec, true),

    // 基本数学・複素補助。
    fixed(builtins::names::cbrt, BuiltinId::Cbrt, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::hypot, BuiltinId::Hypot, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::cis, BuiltinId::Cis, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::polar, BuiltinId::Polar, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::nextPow2, BuiltinId::NextPow2, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::degreeToRadian, BuiltinId::DegreeToRadian, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::degreeToGradian, BuiltinId::DegreeToGradian, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::radianToDegree, BuiltinId::RadianToDegree, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::radianToGradian, BuiltinId::RadianToGradian, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::gradianToDegree, BuiltinId::GradianToDegree, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::gradianToRadian, BuiltinId::GradianToRadian, 1, ArgumentEvaluation::All, true),

    // 集約函数。単一配列または可変長scalarを受ける。
    variadic(builtins::names::sum, BuiltinId::Sum, 0, ArgumentEvaluation::All, true),
    variadic(builtins::names::product, BuiltinId::Product, 0, ArgumentEvaluation::All, true),
    variadic(builtins::names::min, BuiltinId::Min, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::max, BuiltinId::Max, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::mean, BuiltinId::Mean, 1, ArgumentEvaluation::All, true),

    // 記述統計。単一配列またはscalar列を基本形とする。
    variadic(builtins::names::median, BuiltinId::Median, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::mode, BuiltinId::Mode, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::quantile, BuiltinId::Quantile, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::percentile, BuiltinId::Percentile, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::variancePopulation, BuiltinId::VariancePopulation, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::varianceSample, BuiltinId::VarianceSample, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::stddevPopulation, BuiltinId::StddevPopulation, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::stddevSample, BuiltinId::StddevSample, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::geometricMean, BuiltinId::GeometricMean, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::harmonicMean, BuiltinId::HarmonicMean, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::rms, BuiltinId::Rms, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::medianAbsoluteDeviation, BuiltinId::MedianAbsoluteDeviation, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::meanAbsoluteDeviation, BuiltinId::MeanAbsoluteDeviation, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::skewness, BuiltinId::Skewness, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::kurtosisPopulation, BuiltinId::KurtosisPopulation, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::kurtosisSample, BuiltinId::KurtosisSample, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::coefficientVariation, BuiltinId::CoefficientVariation, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::standardError, BuiltinId::StandardError, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::zScore, BuiltinId::ZScore, 3, ArgumentEvaluation::All, true),
    variadic(builtins::names::iqr, BuiltinId::Iqr, 1, ArgumentEvaluation::All, true),
    variadic(builtins::names::trimMean, BuiltinId::TrimMean, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::winsorMean, BuiltinId::WinsorMean, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::winsorized, BuiltinId::Winsorized, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::covariance, BuiltinId::Covariance, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::correlation, BuiltinId::Correlation, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::spearmanCorrelation, BuiltinId::SpearmanCorrelation, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::percentRank, BuiltinId::PercentRank, 2, ArgumentEvaluation::All, true),

    // 配列・ベクトル・行列utility。
    fixed(builtins::names::identity, BuiltinId::Identity, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::zeros, BuiltinId::Zeros, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::matrixGet, BuiltinId::MatrixGet, 3, ArgumentEvaluation::All, true),
    fixed(builtins::names::trace, BuiltinId::Trace, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::rows, BuiltinId::Rows, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::cols, BuiltinId::Cols, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::diag, BuiltinId::Diag, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorAdd, BuiltinId::VectorAdd, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorSubtract, BuiltinId::VectorSubtract, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorScale, BuiltinId::VectorScale, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorDot, BuiltinId::VectorDot, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorCross, BuiltinId::VectorCross, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorNorm, BuiltinId::VectorNorm, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorManhattan, BuiltinId::VectorManhattan, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorEuclidean, BuiltinId::VectorEuclidean, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorNormalize, BuiltinId::VectorNormalize, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorProject, BuiltinId::VectorProject, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorAngle, BuiltinId::VectorAngle, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorReflect, BuiltinId::VectorReflect, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorReflectAxis, BuiltinId::VectorReflectAxis, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::vectorSum, BuiltinId::VectorSum, 1, ArgumentEvaluation::All, true),

    // 0近傍で直接比を作ると桁落ちしやすい安定初等函数。
    fixed(builtins::names::expm1, BuiltinId::Expm1, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::log1p, BuiltinId::Log1p, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::sinc, BuiltinId::Sinc, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::cosc, BuiltinId::Cosc, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::tanc, BuiltinId::Tanc, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::sinhc, BuiltinId::Sinhc, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::tanhc, BuiltinId::Tanhc, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::expc, BuiltinId::Expc, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::log2, BuiltinId::Log2, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::log10, BuiltinId::Log10, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::gamma, BuiltinId::Gamma, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::logGamma, BuiltinId::LogGamma, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::erf, BuiltinId::Erf, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::erfc, BuiltinId::Erfc, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::fresnelC, BuiltinId::FresnelC, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::fresnelS, BuiltinId::FresnelS, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::hypergeometric1F1, BuiltinId::Hypergeometric1F1, 3, ArgumentEvaluation::All, true),
    fixed(builtins::names::hypergeometric2F1, BuiltinId::Hypergeometric2F1, 4, ArgumentEvaluation::All, true),
    fixed(builtins::names::ellipticF, BuiltinId::EllipticF, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::ellipticE, BuiltinId::EllipticE, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::ellipticPi, BuiltinId::EllipticPi, 3, ArgumentEvaluation::All, true),
    fixed(builtins::names::beta, BuiltinId::Beta, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::betaLog, BuiltinId::BetaLog, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::generalizedBinomial, BuiltinId::GeneralizedBinomial, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::fallingFactorial, BuiltinId::FallingFactorial, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::risingFactorial, BuiltinId::RisingFactorial, 2, ArgumentEvaluation::All, true),

    // 乱数。KernelSessionごとのstateを消費するため、数学的なpure functionとして扱わない。
    range(builtins::names::randSeed, BuiltinId::RandSeed, 0, 1, ArgumentEvaluation::All, true),
    range(builtins::names::rand, BuiltinId::Rand, 0, 2, ArgumentEvaluation::All, true),
    range(builtins::names::randInt, BuiltinId::RandInt, 0, 2, ArgumentEvaluation::All, true),
    variadic(builtins::names::choice, BuiltinId::Choice, 1, ArgumentEvaluation::All, true),
    range(builtins::names::randN, BuiltinId::RandN, 0, 2, ArgumentEvaluation::All, true),

    // 初等函数。
    fixed(builtins::names::sqrt, BuiltinId::Sqrt, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::abs, BuiltinId::Abs, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::sign, BuiltinId::Sign, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::re, BuiltinId::Re, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::im, BuiltinId::Im, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::conj, BuiltinId::Conj, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::sin, BuiltinId::Sin, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::cos, BuiltinId::Cos, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::tan, BuiltinId::Tan, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::cot, BuiltinId::Cot, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::sec, BuiltinId::Sec, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::csc, BuiltinId::Csc, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::asin, BuiltinId::Asin, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::acos, BuiltinId::Acos, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::atan, BuiltinId::Atan, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::atan2, BuiltinId::Atan2, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::sinh, BuiltinId::Sinh, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::cosh, BuiltinId::Cosh, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::tanh, BuiltinId::Tanh, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::asinh, BuiltinId::Asinh, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::acosh, BuiltinId::Acosh, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::atanh, BuiltinId::Atanh, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::csch, BuiltinId::Csch, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::sech, BuiltinId::Sech, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::coth, BuiltinId::Coth, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::arg, BuiltinId::Arg, 1, ArgumentEvaluation::All, true),
    range(builtins::names::log, BuiltinId::Log, 1, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::exp, BuiltinId::Exp, 1, ArgumentEvaluation::All, true),

    // 近似・式変形・Solver。
    range(builtins::names::numericalApproximation, BuiltinId::NumericalApproximation, 1, 2,
        ArgumentEvaluation::All, true),
    fixed(builtins::names::precision, BuiltinId::Precision, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::accuracy, BuiltinId::Accuracy, 1, ArgumentEvaluation::All, true),
    range(builtins::names::rationalize, BuiltinId::Rationalize, 1, 2, ArgumentEvaluation::All, true),
    range(builtins::names::simplify, BuiltinId::Simplify, 1, 2, ArgumentEvaluation::All, true),
    range(builtins::names::fullSimplify, BuiltinId::FullSimplify, 1, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::expand, BuiltinId::Expand, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::factor, BuiltinId::Factor, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::collect, BuiltinId::Collect, 2, ArgumentEvaluation::All, true),
    range(builtins::names::solve, BuiltinId::Solve, 2, 3, ArgumentEvaluation::HoldAll, true),

    // 言語・比較用head。
    fixed(builtins::names::set, BuiltinId::Set, 2, ArgumentEvaluation::HoldFirst),
    fixed(builtins::names::setDelayed, BuiltinId::SetDelayed, 2, ArgumentEvaluation::HoldAll),
    fixed(builtins::names::less, BuiltinId::Less, 2),
    fixed(builtins::names::lessEqual, BuiltinId::LessEqual, 2),
    fixed(builtins::names::greater, BuiltinId::Greater, 2),
    fixed(builtins::names::greaterEqual, BuiltinId::GreaterEqual, 2),
    fixed(builtins::names::equal, BuiltinId::Equal, 2),
    fixed(builtins::names::notEqual, BuiltinId::NotEqual, 2),
    variadic(builtins::names::logicalAnd, BuiltinId::LogicalAnd, 0),
    fixed(builtins::names::element, BuiltinId::Element, 2, ArgumentEvaluation::All, true),
    fixed(builtins::names::ifThenElse, BuiltinId::If, 3, ArgumentEvaluation::HoldAll, true),
    fixed(builtins::names::history, BuiltinId::History, 1),
    fixed(builtins::names::inputHistory, BuiltinId::InputHistory, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::outputHistory, BuiltinId::OutputHistory, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::exit, BuiltinId::Exit, 0, ArgumentEvaluation::HoldAll, true),
    fixed(builtins::names::clear, BuiltinId::Clear, 0, ArgumentEvaluation::HoldAll, true),
    fixed(builtins::names::definitions, BuiltinId::Definitions, 0, ArgumentEvaluation::HoldAll, true),
    variadic(builtins::names::undefine, BuiltinId::Undefine, 1, ArgumentEvaluation::HoldAll, true),
    range(builtins::names::angleMode, BuiltinId::AngleMode, 0, 1, ArgumentEvaluation::All, true),
    fixed(builtins::names::unitApplied, BuiltinId::UnitApplied, 2)
};

constexpr BuiltinAliasSpec defaultAliases[] = {
    {"pow", BuiltinId::Power},
    {"fact", BuiltinId::Factorial},
    {"fract", BuiltinId::Frac},
    {"ln", BuiltinId::Log},
    {"real", BuiltinId::Re},
    {"imag", BuiltinId::Im},
    {"mag", BuiltinId::Abs},
    {"unit", BuiltinId::Sign},
    {"csgn", BuiltinId::Sign},
    {"rect", BuiltinId::Polar},
    {"mmul", BuiltinId::MatrixMultiply},
    {"mtranspose", BuiltinId::Transpose},
    {"mdet", BuiltinId::Determinant},
    {"minverse", BuiltinId::Inverse},
    {"mrank", BuiltinId::Rank},
    {"ave", BuiltinId::Mean},
    {"mtrace", BuiltinId::Trace},
    {"mrows", BuiltinId::Rows},
    {"mcols", BuiltinId::Cols},
    {"mdiag", BuiltinId::Diag},
    {"vlength", BuiltinId::VectorNorm},
    {"vdistance", BuiltinId::VectorEuclidean},
    {"vunit", BuiltinId::VectorNormalize},


};

} // namespace

bool BuiltinDefinition::acceptsArity(std::size_t count) const noexcept {
    return count >= minimumArguments && count <= maximumArguments;
}

BuiltinRegistry::BuiltinRegistry(symbols::SymbolTable& symbolTable)
    : symbolTable_(symbolTable) {}

BuiltinRegistry BuiltinRegistry::defaults(symbols::SymbolTable& symbolTable) {
    BuiltinRegistry registry{symbolTable};
    for (const BuiltinSpec& spec : defaultBuiltins)
        registry.add(
            spec.name,
            spec.id,
            spec.minimumArguments,
            spec.maximumArguments,
            spec.argumentEvaluation,
            spec.sourceCallable);
    for (const BuiltinAliasSpec& alias : defaultAliases)
        registry.addAlias(alias.name, alias.targetId, alias.sourceCallable);
    return registry;
}

void BuiltinRegistry::add(
    std::string_view name,
    BuiltinId id,
    std::size_t minimumArguments,
    std::size_t maximumArguments,
    ArgumentEvaluation argumentEvaluation,
    bool sourceCallable) {
    if (minimumArguments > maximumArguments)
        throw std::invalid_argument("Builtin argument range is invalid");

    const expression::Symbol interned = symbolTable_.intern(name);
    const auto [iterator, inserted] = definitions_.emplace(
        interned.id(),
        BuiltinDefinition{
            interned,
            id,
            minimumArguments,
            maximumArguments,
            argumentEvaluation,
            sourceCallable
        });
    static_cast<void>(iterator);
    if (!inserted)
        throw std::invalid_argument("Builtin is already registered: " + interned.name());

    const auto [symbolIterator, symbolInserted] = symbolsById_.emplace(id, interned);
    static_cast<void>(symbolIterator);
    if (!symbolInserted)
        throw std::invalid_argument("Builtin ID is already registered");
}

void BuiltinRegistry::addAlias(
    std::string_view name,
    BuiltinId targetId,
    bool sourceCallable) {
    const auto canonical = symbolsById_.find(targetId);
    if (canonical == symbolsById_.end())
        throw std::invalid_argument("Builtin alias target is not registered");

    const BuiltinDefinition* target = find(canonical->second);
    if (!target)
        throw std::logic_error("Builtin canonical definition is missing");

    const expression::Symbol alias = symbolTable_.intern(name);
    const auto [iterator, inserted] = definitions_.emplace(
        alias.id(),
        BuiltinDefinition{
            alias,
            targetId,
            target->minimumArguments,
            target->maximumArguments,
            target->argumentEvaluation,
            sourceCallable
        });
    static_cast<void>(iterator);
    if (!inserted)
        throw std::invalid_argument("Builtin is already registered: " + alias.name());
}

const BuiltinDefinition* BuiltinRegistry::find(const expression::Symbol& symbol) const noexcept {
    const auto iterator = definitions_.find(symbol.id());
    return iterator == definitions_.end() ? nullptr : &iterator->second;
}

const BuiltinDefinition* BuiltinRegistry::find(std::string_view name) const noexcept {
    const expression::Symbol interned = symbolTable_.find(name);
    return interned.valid() ? find(interned) : nullptr;
}

bool BuiltinRegistry::contains(const expression::Symbol& symbol) const noexcept {
    return find(symbol) != nullptr;
}

bool BuiltinRegistry::contains(std::string_view name) const noexcept {
    return find(name) != nullptr;
}

const expression::Symbol& BuiltinRegistry::symbol(BuiltinId id) const {
    const auto iterator = symbolsById_.find(id);
    if (iterator == symbolsById_.end())
        throw std::logic_error("Builtin ID is not registered");
    return iterator->second;
}

std::unordered_set<std::string> BuiltinRegistry::sourceFunctionNames() const {
    std::unordered_set<std::string> result;
    for (const auto& [id, definition] : definitions_) {
        static_cast<void>(id);
        if (definition.sourceCallable)
            result.insert(definition.symbol.name());
    }
    return result;
}

std::size_t BuiltinRegistry::size() const noexcept {
    return definitions_.size();
}

const BuiltinRegistry& defaultBuiltinRegistry() {
    static const BuiltinRegistry registry = BuiltinRegistry::defaults(symbols::defaultSymbolTable());
    return registry;
}

} // namespace mmcal::evaluation
