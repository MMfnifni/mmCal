// 入力テキスト文書
#include "source_document.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::source {

SourceDocument::SourceDocument(
    std::size_t inputNumber,
    std::shared_ptr<const std::string> text)
    : inputNumber_(inputNumber), text_(std::move(text)) {
    if (!text_)
        throw std::invalid_argument("Source document text cannot be null");
}

std::size_t SourceDocument::inputNumber() const noexcept {
    return inputNumber_;
}

std::string_view SourceDocument::text() const noexcept {
    return *text_;
}

const std::shared_ptr<const std::string>& SourceDocument::textStorage() const noexcept {
    return text_;
}

} // namespace mmcal::source
