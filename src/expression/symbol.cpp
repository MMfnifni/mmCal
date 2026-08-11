// 式中のSymbol表現
#include "symbol.hpp"

#include "symbols/symbol_table.hpp"

#include <utility>

namespace mmcal::expression {

Symbol::Symbol(std::string_view name)
    : Symbol(symbols::defaultSymbolTable().intern(name)) {}

Symbol::Symbol(std::shared_ptr<const symbols::SymbolEntry> entry) noexcept
    : entry_(std::move(entry)) {}

symbols::SymbolId Symbol::id() const noexcept {
    return entry_ ? entry_->id : symbols::SymbolId{};
}

const std::string& Symbol::name() const noexcept {
    static const std::string empty;
    return entry_ ? entry_->name : empty;
}

std::string_view Symbol::view() const noexcept {
    return name();
}

bool Symbol::valid() const noexcept {
    return static_cast<bool>(entry_);
}

bool Symbol::sameIdentity(const Symbol& rhs) const noexcept {
    return entry_ == rhs.entry_;
}

std::strong_ordering Symbol::operator<=>(const Symbol& rhs) const noexcept {
    return view() <=> rhs.view();
}

bool Symbol::operator==(const Symbol& rhs) const noexcept {
    return view() == rhs.view();
}

} // namespace mmcal::expression
