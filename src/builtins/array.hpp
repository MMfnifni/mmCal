#pragma once

#include "expression/expr.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluateDimensions(
    std::span<const expression::Expr> arguments);
[[nodiscard]] expression::Expr evaluateArrayRank(
    std::span<const expression::Expr> arguments);
[[nodiscard]] expression::Expr evaluateLength(
    std::span<const expression::Expr> arguments);
[[nodiscard]] expression::Expr evaluateArrayGet(
    std::span<const expression::Expr> arguments);
[[nodiscard]] expression::Expr evaluateReshape(
    std::span<const expression::Expr> arguments);

} // namespace mmcal::builtins
