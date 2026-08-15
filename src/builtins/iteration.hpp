#pragma once

#include "expression/expr.hpp"

#include <span>
#include <string_view>
#include <vector>

namespace mmcal::builtins {

// range/tableが共有するexact有限列生成。1引数は1..end，2引数はstart..end，3引数はstep指定。
[[nodiscard]] std::vector<expression::Expr> exactRangeValues(
    std::span<const expression::Expr> arguments,
    std::string_view caller);

[[nodiscard]] expression::Expr evaluateRange(
    std::span<const expression::Expr> arguments);

} // namespace mmcal::builtins
