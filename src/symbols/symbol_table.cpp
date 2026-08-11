// 文字列とSymbolIDのintern管理
#include "symbol_table.hpp"

#include "expression/symbol.hpp"

#include <atomic>
#include <stdexcept>
#include <utility>

namespace mmcal::symbols {
namespace {

// SymbolIdはプロセス内で重複しない。SymbolTable自体の状態はグローバル共有しない。
std::atomic<std::uint64_t> nextSymbolId{1};

[[nodiscard]] SymbolId allocateSymbolId() {
    const std::uint64_t value = nextSymbolId.fetch_add(1, std::memory_order_relaxed);
    if (value == 0)
        throw std::overflow_error("Symbol ID space is exhausted");
    return SymbolId{value};
}

} // namespace

expression::Symbol SymbolTable::intern(std::string_view name) {
    if (name.empty())
        throw std::invalid_argument("Symbol name cannot be empty");

    if (const auto iterator = entries_.find(name); iterator != entries_.end())
        return expression::Symbol{iterator->second};

    auto entry = std::make_shared<const SymbolEntry>(SymbolEntry{
        allocateSymbolId(),
        std::string{name}
    });
    const auto [iterator, inserted] = entries_.emplace(entry->name, entry);
    static_cast<void>(inserted);
    return expression::Symbol{iterator->second};
}

expression::Symbol SymbolTable::find(std::string_view name) const {
    const auto iterator = entries_.find(name);
    return iterator == entries_.end() ? expression::Symbol{} : expression::Symbol{iterator->second};
}

bool SymbolTable::contains(std::string_view name) const noexcept {
    return entries_.contains(name);
}

std::size_t SymbolTable::size() const noexcept {
    return entries_.size();
}

SymbolTable& defaultSymbolTable() {
    static SymbolTable table;
    return table;
}

} // namespace mmcal::symbols
