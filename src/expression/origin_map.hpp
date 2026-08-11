#pragma once

#include "expr.hpp"
#include "source/source_document.hpp"

#include <memory>
#include <optional>
#include <unordered_map>

namespace mmcal::expression {

// Kernel式本体を汚さず、式ノードと入力文書上の位置を関連付ける補助表。
class OriginMap final {
public:
    OriginMap() = default;
    explicit OriginMap(std::shared_ptr<const source::SourceDocument> document)
        : document_(std::move(document)) {}

    void record(const Expr& expression, source::SourceSpan span) {
        origins_.try_emplace(
            expression.identity(),
            source::SourceReference{document_, span});
    }

    [[nodiscard]] std::optional<source::SourceReference> find(
        const Expr& expression) const {
        const auto iterator = origins_.find(expression.identity());
        if (iterator == origins_.end())
            return std::nullopt;

        return iterator->second;
    }

    [[nodiscard]] bool contains(const Expr& expression) const {
        return origins_.contains(expression.identity());
    }

    [[nodiscard]] const std::shared_ptr<const source::SourceDocument>& document() const noexcept {
        return document_;
    }

private:
    std::shared_ptr<const source::SourceDocument> document_;
    std::unordered_map<const void*, source::SourceReference> origins_;
};

} // namespace mmcal::expression
