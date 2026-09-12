#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "symbolic/risch_core.hpp"

namespace mmcal::symbolic::risch {

// compactなLRT residue certificateを，有限個のRootとlogからなる公開Exprへ変換する。
// certificateを再確認できない入力はmaterializeしない。
[[nodiscard]] RischStageResult<expression::Expr> materializeLrtLogarithms(
    const RationalFunction& properSquareFreePart,
    const LrtResult& result,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const RischOptions& options = {});

[[nodiscard]] RischStageResult<expression::Expr> materializeRationalRischResult(
    const RischResult& result,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const RischOptions& options = {});

} // namespace mmcal::symbolic::risch
