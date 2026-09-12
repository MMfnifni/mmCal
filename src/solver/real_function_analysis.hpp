#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"

#include <optional>
#include <vector>

namespace mmcal::solver {

enum class RealIntervalMonotonicity {
    Unknown,
    Increasing,
    Decreasing,
    Constant
};

// lower==nullopt は -Infinity，upper==nullopt は +Infinity を表す。
struct RealDomainInterval final {
    std::optional<expression::Expr> lower;
    bool lowerInclusive = false;
    std::optional<expression::Expr> upper;
    bool upperInclusive = false;

    [[nodiscard]] bool operator==(const RealDomainInterval&) const = default;
};

// range境界はInfinity sentinelを含むExprとして保持する。
struct RealValueRange final {
    expression::Expr lower;
    bool lowerInclusive = false;
    expression::Expr upper;
    bool upperInclusive = false;

    [[nodiscard]] bool operator==(const RealValueRange&) const = default;
};

struct RealIntervalFunctionAnalysis final {
    RealDomainInterval domain;
    RealIntervalMonotonicity monotonicity = RealIntervalMonotonicity::Unknown;
    std::optional<expression::Expr> lowerLimit;
    std::optional<expression::Expr> upperLimit;
    std::optional<RealValueRange> range;

    [[nodiscard]] bool operator==(const RealIntervalFunctionAnalysis&) const = default;
};

struct RealDomainAnalysis final {
    // trueならintervalsのunionが実数入力に対する完全な実定義域である。
    bool complete = false;
    std::vector<RealDomainInterval> intervals;
};

struct RealFunctionAnalysis final {
    // trueならpiecesのunionが「実数入力に対して式が定義され実数値を取る集合」全体である。
    bool domainComplete = false;
    std::vector<RealDomainInterval> domain;
    std::optional<expression::Expr> derivative;
    // 定義域をexactな臨界点でさらに分割した解析区間。
    std::vector<RealIntervalFunctionAnalysis> pieces;
};

// Plot等が微分・値域解析を伴わずに利用する軽量な実定義域解析。
// complete=falseならintervalsを完全なdomainとして扱ってはならない。
[[nodiscard]] RealDomainAnalysis analyzeRealDomain(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// 一変数実函数をexact knowledgeだけで解析する内部基盤。
// 数値samplingは定義域・値域・単調性の証明には使わない。
// 初版はalgebraicなdomain条件と，exactに解けるcritical pointへbounded-workで限定する。
[[nodiscard]] RealFunctionAnalysis analyzeRealFunction(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::solver
