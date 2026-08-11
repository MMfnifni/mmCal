#pragma once

#include "expression/expr.hpp"
#include "expression/symbol.hpp"

#include <optional>
#include <span>

namespace mmcal::builtins {

[[nodiscard]] std::optional<expression::Expr> evaluatePrecision(
    std::span<const expression::Expr> arguments,
    const expression::Symbol& infinity);
[[nodiscard]] std::optional<expression::Expr> evaluateAccuracy(
    std::span<const expression::Expr> arguments,
    const expression::Symbol& infinity);
[[nodiscard]] std::optional<expression::Expr> evaluateRationalize(
    std::span<const expression::Expr> arguments);

} // namespace mmcal::builtins
