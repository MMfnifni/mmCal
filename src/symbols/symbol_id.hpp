// Symbolの安定ID
#pragma once

#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>

namespace mmcal::symbols {

// intern済みシンボルを軽量に識別するID。0は無効値として予約する。
class SymbolId final {
public:
    constexpr SymbolId() noexcept = default;
    explicit constexpr SymbolId(std::uint64_t value) noexcept : value_(value) {}

    [[nodiscard]] constexpr std::uint64_t value() const noexcept { return value_; }
    [[nodiscard]] constexpr bool valid() const noexcept { return value_ != 0; }
    explicit constexpr operator bool() const noexcept { return valid(); }

    [[nodiscard]] constexpr std::strong_ordering operator<=>(const SymbolId&) const noexcept = default;
    [[nodiscard]] constexpr bool operator==(const SymbolId&) const noexcept = default;

private:
    std::uint64_t value_ = 0;
};

struct SymbolIdHash final {
    [[nodiscard]] std::size_t operator()(SymbolId id) const noexcept {
        return std::hash<std::uint64_t>{}(id.value());
    }
};

} // namespace mmcal::symbols
