#pragma once

#include "source/source_span.hpp"

#include <string_view>

namespace mmcal::syntax {

enum class TokenKind {
    End,
    Number,
    String,
    Identifier,
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    Bang,
    Assign,
    Less,
    LessEqual,
    Greater,
    GreaterEqual,
    EqualEqual,
    BangEqual,
    LParen,
    RParen,
    LBracket,
    RBracket,
    LBrace,
    RBrace,
    Comma,
    Percent,
    At
};

struct Token final {
    TokenKind kind = TokenKind::End;
    source::SourceSpan span;

    [[nodiscard]] std::string_view text(std::string_view sourceText) const noexcept {
        return sourceText.substr(span.begin.offset, span.end.offset - span.begin.offset);
    }
};

[[nodiscard]] std::string_view tokenKindName(TokenKind kind) noexcept;

} // namespace mmcal::syntax
