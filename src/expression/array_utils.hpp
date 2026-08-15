#pragma once

#include "expr.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace mmcal::expression {

// shapeから要素数を計算する。shapeは1次元以上、積はsize_t範囲内でなければならない。
[[nodiscard]] std::size_t arrayElementCount(std::span<const std::size_t> shape);

// Arrayの各leafを評価した後に、同一shapeのArrayが返った場合だけ末尾次元へ吸収する。
// scalar/Array混在やchild shape不一致はrectangular Arrayを壊すため拒否する。
[[nodiscard]] Expr rebuildEvaluatedArray(
    std::vector<std::size_t> outerShape,
    std::vector<Expr> evaluatedElements);


// packed ArrayのExpr pageだけを評価した結果で再構築する。scalar結果だけなら未変更pageを共有し，
// child Array/Listが返った場合だけ従来のflatten規則へfallbackする。
[[nodiscard]] Expr rebuildEvaluatedArray(
    const ArrayExpr& source,
    std::span<const std::size_t> replacedIndices,
    std::vector<Expr> evaluatedElements);

// {}の一般構築。矩形かつ同shapeならdense Arrayへ正規化し，そうでなければListとして保持する。
[[nodiscard]] Expr braceValue(std::vector<Expr> elements);

// 一般brace値の先頭から共通して保証できるdimensionsを返す。
// 例: {{1,2},{3}} -> {2}, {Q(3x2),R(2x2)} -> {2}。
[[nodiscard]] std::vector<std::size_t> commonBraceDimensions(const Expr& expression);

// brace literalだけでshapeを再入力可能か。0次元の後ろに次元が残るshapeは{}表記だけでは情報を失う。
[[nodiscard]] bool braceLiteralPreservesShape(std::span<const std::size_t> shape) noexcept;

} // namespace mmcal::expression
