#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "random/random_engine.hpp"

#include <span>

namespace mmcal::builtins {

[[nodiscard]] expression::Expr evaluateRandSeed(
    std::span<const expression::Expr> arguments,
    random::RandomEngine& engine);
[[nodiscard]] expression::Expr evaluateRand(
    std::span<const expression::Expr> arguments,
    random::RandomEngine& engine);
[[nodiscard]] expression::Expr evaluateRandInt(
    std::span<const expression::Expr> arguments,
    random::RandomEngine& engine);
[[nodiscard]] expression::Expr evaluateChoice(
    std::span<const expression::Expr> arguments,
    random::RandomEngine& engine);
[[nodiscard]] expression::Expr evaluateRandN(
    std::span<const expression::Expr> arguments,
    random::RandomEngine& engine,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::builtins
