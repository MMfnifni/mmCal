#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "math_registry.hpp"

#include <optional>

namespace mmcal::mathematics {

// principal Arg。値域は (-Pi, Pi]。0では未定義。
[[nodiscard]] std::optional<expression::Expr> simplifyExactArg(
    const expression::Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics);

// principal Log。Log[z] = ln|z| + I Arg[z]、Argは (-Pi, Pi]。
// exactに安全な変形だけを行い、一般値はnulloptで記号式を維持する。
[[nodiscard]] std::optional<expression::Expr> simplifyExactLog(
    const expression::Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics);

// 任意底対数 log[base, value] = principal Log[value] / principal Log[base]。
// base=0,1 または value=0 は未定義。exactに確定できる冪関係だけ畳み込む。
[[nodiscard]] std::optional<expression::Expr> simplifyExactLog(
    const expression::Expr& base,
    const expression::Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics);

// 複素指数函数Exp。全複素平面で一価。
// Exp[0]=1, Exp[1]=E, Eulerの特殊角などexactに確定できる場合だけ簡約する。
[[nodiscard]] std::optional<expression::Expr> simplifyExactExp(
    const expression::Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics);

} // namespace mmcal::mathematics
