#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <span>

namespace mmcal::builtins {

// の特殊函数群。exactに確定できる値だけここで簡約し、一般値は記号式として保持して N[...] のcertified backendへ渡す。
[[nodiscard]] expression::Expr evaluateSpecialFunction(
    evaluation::BuiltinId id,
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::builtins
