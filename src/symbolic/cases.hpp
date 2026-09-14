#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace mmcal::symbolic::detail {

// Cases/CaseBranchのAST形は微分・積分・極限・at/toNormalで共通である。
// branchのcondition省略規則と構造検査を一箇所へ固定し，subsystem間で解釈がずれないようにする。
struct CaseBranchView final {
    const expression::Expr* value = nullptr;
    const expression::Expr* condition = nullptr;

    [[nodiscard]] bool hasCondition() const noexcept { return condition != nullptr; }
};

[[nodiscard]] inline std::optional<CaseBranchView> caseBranchView(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!builtins.isCallTo(expression, evaluation::BuiltinId::CaseBranch))
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.empty() || arguments.size() > 2)
        return std::nullopt;
    return CaseBranchView{
        &arguments[0],
        arguments.size() == 2 ? &arguments[1] : nullptr};
}

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
