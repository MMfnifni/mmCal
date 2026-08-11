// 入力範囲
#pragma once

#include <compare>
#include <cstddef>

namespace mmcal::source {

struct SourcePosition final {
    std::size_t offset = 0;
    std::size_t line = 1;
    std::size_t column = 1;

    [[nodiscard]] auto operator<=>(const SourcePosition&) const = default;
};

struct SourceSpan final {
    SourcePosition begin;
    SourcePosition end;

    [[nodiscard]] bool empty() const noexcept {
        return begin.offset == end.offset;
    }

    [[nodiscard]] auto operator<=>(const SourceSpan&) const = default;
};

[[nodiscard]] inline SourceSpan combine(SourceSpan lhs, SourceSpan rhs) noexcept {
    return SourceSpan{lhs.begin, rhs.end};
}

} // namespace mmcal::source
