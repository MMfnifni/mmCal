#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace mmcal::symbolic::detail {

// Cases/CaseBranchのAST形は微分・積分・極限で同一である。
// branchのcondition省略規則を一箇所へ固定し，subsystem間で構造がずれないようにする。
[[nodiscard]] inline expression::Expr makeCaseBranch(
    const evaluation::BuiltinRegistry& builtins,
    expression::Expr value,
    std::optional<expression::Expr> condition = std::nullopt) {
    std::vector<expression::Expr> arguments{std::move(value)};
    if (condition)
        arguments.push_back(std::move(*condition));
    return expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::CaseBranch), std::move(arguments));
}

[[nodiscard]] inline expression::Expr makeCases(
    const evaluation::BuiltinRegistry& builtins,
    std::vector<expression::Expr> branches) {
    return expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Cases), std::move(branches));
}

} // namespace mmcal::symbolic::detail
