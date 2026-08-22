// token種別と位置情報
#include "token.hpp"

namespace mmcal::syntax {

std::string_view tokenKindName(TokenKind kind) noexcept {
    switch (kind) {
    case TokenKind::End:
        return "end of input";
    case TokenKind::Number:
        return "number";
    case TokenKind::String:
        return "string";
    case TokenKind::Identifier:
        return "identifier";
    case TokenKind::Plus:
        return "+";
    case TokenKind::Minus:
        return "-";
    case TokenKind::Star:
        return "*";
    case TokenKind::Slash:
        return "/";
    case TokenKind::Caret:
        return "^";
    case TokenKind::Bang:
        return "!";
    case TokenKind::Assign:
        return ":=";
    case TokenKind::Less:
        return "<";
    case TokenKind::LessEqual:
        return "<=";
    case TokenKind::Greater:
        return ">";
    case TokenKind::GreaterEqual:
        return ">=";
    case TokenKind::EqualEqual:
        return "==";
    case TokenKind::BangEqual:
        return "!=";
    case TokenKind::LParen:
        return "(";
    case TokenKind::RParen:
        return ")";
    case TokenKind::LBracket:
        return "[";
    case TokenKind::RBracket:
        return "]";
    case TokenKind::LBrace:
        return "{";
    case TokenKind::RBrace:
        return "}";
    case TokenKind::Comma:
        return ",";
    case TokenKind::Semicolon:
        return ";";
    case TokenKind::Percent:
        return "%";
    case TokenKind::At:
        return "@";
    }

    return "token";
}

} // namespace mmcal::syntax
