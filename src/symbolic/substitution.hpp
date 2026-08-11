#pragma once

#include "expression/expr.hpp"
#include "expression/symbol.hpp"

namespace mmcal::symbolic {

// 式中の自由な同一Symbolを値へ構造的に置換する。
// builtin headは式引数ではないため置換対象にしない。
[[nodiscard]] expression::Expr substituteSymbol(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const expression::Expr& value);

} // namespace mmcal::symbolic
