// 構文木AST
#include "syntax_tree.hpp"

#include <utility>

namespace mmcal::syntax {

SyntaxTree::SyntaxTree(
    std::shared_ptr<const std::string> sourceText,
    std::vector<Token> tokens,
    SyntaxNodePtr root)
    : sourceText_(std::move(sourceText)),
      tokens_(std::move(tokens)),
      root_(std::move(root)) {}

std::string_view SyntaxTree::sourceText() const noexcept {
    return *sourceText_;
}

const std::vector<Token>& SyntaxTree::tokens() const noexcept {
    return tokens_;
}

const SyntaxNode& SyntaxTree::root() const noexcept {
    return *root_;
}

const SyntaxNodePtr& SyntaxTree::rootPtr() const noexcept {
    return root_;
}

std::string_view SyntaxTree::text(source::SourceSpan span) const noexcept {
    return sourceText().substr(
        span.begin.offset,
        span.end.offset - span.begin.offset);
}

} // namespace mmcal::syntax
