#pragma once

#include "expression/expr.hpp"
#include "expression/origin_map.hpp"
#include "expression/symbol.hpp"
#include "symbols/symbol_id.hpp"

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace mmcal::evaluation {

struct UserFunctionDefinition final {
    UserFunctionDefinition(
        expression::Symbol functionName,
        std::vector<expression::Symbol> functionParameters,
        expression::Expr functionBody,
        expression::OriginMap functionOrigins = {});

    expression::Symbol name;
    std::vector<expression::Symbol> parameters;
    expression::Expr body;
    expression::OriginMap origins;

    [[nodiscard]] std::size_t arity() const noexcept;
};

// ユーザー定義函数をintern済みSymbolIdと引数個数ごとに保持する。
class UserFunctionRegistry final {
public:
    void define(UserFunctionDefinition definition);

    [[nodiscard]] const UserFunctionDefinition* find(
        const expression::Symbol& name,
        std::size_t arity) const noexcept;
    [[nodiscard]] bool contains(const expression::Symbol& name) const noexcept;
    [[nodiscard]] std::vector<std::size_t> arities(const expression::Symbol& name) const;
    [[nodiscard]] std::set<std::string> names() const;
    [[nodiscard]] std::vector<UserFunctionDefinition> definitions() const;

    [[nodiscard]] bool erase(const expression::Symbol& name);
    void clear() noexcept;
    [[nodiscard]] std::size_t size() const noexcept;

private:
    struct FunctionEntry final {
        expression::Symbol name;
        std::map<std::size_t, UserFunctionDefinition> overloads;
    };

    std::map<symbols::SymbolId, FunctionEntry> definitions_;
};

} // namespace mmcal::evaluation
