#pragma once

#include "error/error_message.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace mmcal::builtins {

// builtin共通のarity検査。CLIへ出る診断文を各実装で重複させない。
inline void requireArity(
    std::span<const expression::Expr> arguments,
    std::size_t expected,
    std::string_view name) {
    if (arguments.size() != expected)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects " + std::to_string(expected) + " argument(s)");
}

// 評価せず同じbuiltin callを保持する共通経路。
[[nodiscard]] inline expression::Expr holdBuiltin(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    evaluation::BuiltinId id,
    std::string_view name,
    std::size_t arity) {
    requireArity(arguments, arity, name);
    return expression::Expr::call(
        registry.symbol(id),
        std::vector<expression::Expr>{arguments.begin(), arguments.end()});
}

[[nodiscard]] inline expression::Expr holdUnaryBuiltin(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    evaluation::BuiltinId id,
    std::string_view name) {
    return holdBuiltin(arguments, registry, id, name, 1);
}

} // namespace mmcal::builtins
