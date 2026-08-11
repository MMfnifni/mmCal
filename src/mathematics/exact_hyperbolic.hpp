#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "math_registry.hpp"

#include <optional>

namespace mmcal::mathematics {

// 双曲線函数・逆双曲線函数について、branchに依存せずexactに確定できる値だけ簡約する。
// 一般式をExp/Logへ常時展開すると式が肥大化するため、記号入力はそのまま保持する。
[[nodiscard]] std::optional<expression::Expr> simplifyExactHyperbolic(
    FunctionId function,
    const expression::Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics);

} // namespace mmcal::mathematics
