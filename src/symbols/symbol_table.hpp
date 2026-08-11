#pragma once

#include "symbol_id.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>

namespace mmcal::expression {
class Symbol;
}

namespace mmcal::symbols {

// Symbolが共有するinternエントリ。文字列本体を式ノードごとに複製しない。
struct SymbolEntry final {
    SymbolId id;
    std::string name;
};

// 名前を一意なSymbolEntryへinternする表。KernelSessionが所有することを基本とする。
class SymbolTable final {
public:
    SymbolTable() = default;

    [[nodiscard]] expression::Symbol intern(std::string_view name);
    [[nodiscard]] expression::Symbol find(std::string_view name) const;
    [[nodiscard]] bool contains(std::string_view name) const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;

private:
    struct TransparentStringHash final {
        using is_transparent = void;

        [[nodiscard]] std::size_t operator()(std::string_view value) const noexcept {
            return std::hash<std::string_view>{}(value);
        }
    };

    struct TransparentStringEqual final {
        using is_transparent = void;

        [[nodiscard]] bool operator()(std::string_view lhs, std::string_view rhs) const noexcept {
            return lhs == rhs;
        }
    };

    std::unordered_map<
        std::string,
        std::shared_ptr<const SymbolEntry>,
        TransparentStringHash,
        TransparentStringEqual> entries_;
};

// 単体テストや独立したExpr構築用。通常のKernelSessionは専用SymbolTableを所有する。
[[nodiscard]] SymbolTable& defaultSymbolTable();

} // namespace mmcal::symbols
