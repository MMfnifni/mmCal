#pragma once

#include "source_span.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>

namespace mmcal::source {

// 1回の入力文字列と入力番号を共有所有し、後からエラー位置を再現できるようにする。
class SourceDocument final {
public:
    SourceDocument(
        std::size_t inputNumber,
        std::shared_ptr<const std::string> text);

    [[nodiscard]] std::size_t inputNumber() const noexcept;
    [[nodiscard]] std::string_view text() const noexcept;
    [[nodiscard]] const std::shared_ptr<const std::string>& textStorage() const noexcept;

private:
    std::size_t inputNumber_ = 0;
    std::shared_ptr<const std::string> text_;
};

struct SourceReference final {
    std::shared_ptr<const SourceDocument> document;
    SourceSpan span;

    [[nodiscard]] bool valid() const noexcept {
        return static_cast<bool>(document);
    }
};

} // namespace mmcal::source
