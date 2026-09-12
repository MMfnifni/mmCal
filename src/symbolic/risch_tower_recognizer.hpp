#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "symbolic/differential_tower.hpp"
#include "symbols/symbol_table.hpp"

#include <cstddef>

namespace mmcal::symbolic::risch {

struct DifferentialTowerRecognitionOptions final {
    std::size_t maximumTowerDepth = 16;
};

enum class DifferentialTowerRecognitionStatus {
    Complete,
    Partial,
    ResourceLimit,
    InvalidBaseVariable
};

struct DifferentialTowerRecognitionResult final {
    DifferentialTower tower;
    expression::Expr rewrittenExpression;
    DifferentialTowerRecognitionStatus status =
        DifferentialTowerRecognitionStatus::Partial;

    DifferentialTowerRecognitionResult(
        DifferentialTower recognizedTower,
        expression::Expr rewritten,
        DifferentialTowerRecognitionStatus recognitionStatus);

    [[nodiscard]] bool complete() const noexcept;
};

// Log[u] を primitive，Exp[u] を exponential generator としてbottom-upに認識する。
// differential coefficient が既に構築済みの下位有理函数体へ属す場合だけtowerへ加える。
// これは構文的なtower候補の構築であり，generatorの超越独立性は証明しない。
// したがって認識結果だけを非初等性certificateとして使ってはならない。
[[nodiscard]] DifferentialTowerRecognitionResult recognizeDifferentialTower(
    const expression::Expr& expression,
    const expression::Symbol& baseVariable,
    symbols::SymbolTable& symbols,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    DifferentialTowerRecognitionOptions options = {});

// 四則演算と整数冪だけを許す，現在のtowerに対する保守的なmembership判定。
[[nodiscard]] bool isRationalInDifferentialTower(
    const expression::Expr& expression,
    const DifferentialTower& tower,
    const evaluation::BuiltinRegistry& builtins);

} // namespace mmcal::symbolic::risch
