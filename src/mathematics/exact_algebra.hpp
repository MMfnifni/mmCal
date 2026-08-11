#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "numeric/rational.hpp"

#include <span>
#include <vector>

namespace mmcal::mathematics {

// 厳密な有理数係数を、記号式そのものを近似せずに掛ける。
// 特に q * (expr / r) のような形では q/r をRationalとして先に約分し、
// 2 * ((a-b)/4) -> (a-b)/2 のような自明な係数整理を行う。
[[nodiscard]] expression::Expr scaleExactExpression(
    const numeric::Rational& coefficient,
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins);

// 構造的に完全一致する項だけを安全にまとめる。
// x+x -> 2*x は行うが、sqrt[8] と 2*sqrt[2] のように別の数学知識が必要な
// 同値判定はここでは行わない。後者はradical normalizerの責務に分離する。
[[nodiscard]] std::vector<expression::Expr> combineStructurallyIdenticalTerms(
    std::span<const expression::Expr> terms,
    const evaluation::BuiltinRegistry& builtins);

} // namespace mmcal::mathematics
