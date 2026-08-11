#pragma once

#include "expression/expr.hpp"

#include <string>
#include <string_view>

namespace mmcal::formatting {

[[nodiscard]] std::string formatExpr(
    const expression::Expr& expression,
    unsigned radix = 10);

// 表示済み数式中の10進小数について、小数部末尾の不要な0だけを除去する。
// 文字列リテラル内は変更しないため、CLIのpresentation処理にも安全に利用できる。
[[nodiscard]] std::string trimRedundantFractionalZeros(std::string_view text);

} // namespace mmcal::formatting
