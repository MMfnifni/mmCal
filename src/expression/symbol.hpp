#pragma once

#include "symbols/symbol_id.hpp"

#include <compare>
#include <memory>
#include <string>
#include <string_view>

namespace mmcal::symbols {
struct SymbolEntry;
class SymbolTable;
}

namespace mmcal::expression {

// intern済みシンボルへの軽量ハンドル。名前の実体はSymbolTable側で共有する。
class Symbol final {
public:
    Symbol() noexcept = default;

    // 単体利用向けの互換入口。通常のKernel処理ではSession所有SymbolTableからinternする。
    explicit Symbol(std::string_view name);

    [[nodiscard]] symbols::SymbolId id() const noexcept;
    [[nodiscard]] const std::string& name() const noexcept;
    [[nodiscard]] std::string_view view() const noexcept;
    [[nodiscard]] bool valid() const noexcept;
    [[nodiscard]] bool sameIdentity(const Symbol& rhs) const noexcept;

    // 異なるSymbolTable由来でも同名シンボルは式として同値とみなす。
    [[nodiscard]] std::strong_ordering operator<=>(const Symbol& rhs) const noexcept;
    [[nodiscard]] bool operator==(const Symbol& rhs) const noexcept;

private:
    std::shared_ptr<const symbols::SymbolEntry> entry_;

    explicit Symbol(std::shared_ptr<const symbols::SymbolEntry> entry) noexcept;
    friend class symbols::SymbolTable;
};

} // namespace mmcal::expression
