#pragma once

#include "syntax_node.hpp"
#include "token.hpp"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace mmcal::syntax {

// 元の入力、トークン列、位置情報付きASTをまとめて保持する。
class SyntaxTree final {
public:
    SyntaxTree(
        std::shared_ptr<const std::string> sourceText,
        std::vector<Token> tokens,
        SyntaxNodePtr root);

    [[nodiscard]] std::string_view sourceText() const noexcept;
    [[nodiscard]] const std::vector<Token>& tokens() const noexcept;
    [[nodiscard]] const SyntaxNode& root() const noexcept;
    [[nodiscard]] const SyntaxNodePtr& rootPtr() const noexcept;
    [[nodiscard]] std::string_view text(source::SourceSpan span) const noexcept;

private:
    std::shared_ptr<const std::string> sourceText_;
    std::vector<Token> tokens_;
    SyntaxNodePtr root_;
};

} // namespace mmcal::syntax
