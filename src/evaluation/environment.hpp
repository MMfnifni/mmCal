#pragma once

#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "symbols/symbol_id.hpp"

#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mmcal::evaluation {

// シンボルと式の束縛を保持する評価環境。
// intern済みSymbolIdをキーにし、文字列比較を評価ホットパスから除く。
class Environment final {
public:
    void set(expression::Symbol symbol, expression::Expr value);
    void assign(expression::Symbol symbol, expression::Expr value);
    void setLocal(expression::Symbol symbol, expression::Expr value);

    [[nodiscard]] const expression::Expr* find(const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] bool contains(const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] bool containsLocal(const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] bool erase(const expression::Symbol& symbol);
    [[nodiscard]] std::vector<std::pair<expression::Symbol, expression::Expr>> definitions() const;

    void pushScope();
    void popScope();
    [[nodiscard]] std::size_t localDepth() const noexcept;

    void clear() noexcept;
    [[nodiscard]] std::size_t size() const noexcept;

private:
    struct Binding final {
        expression::Symbol symbol;
        expression::Expr value;
    };

    using Bindings = std::unordered_map<symbols::SymbolId, Binding, symbols::SymbolIdHash>;

    Bindings bindings_;
    std::vector<Bindings> localScopes_;
};

} // namespace mmcal::evaluation
