// 構文木node
#include "syntax_node.hpp"

namespace mmcal::syntax {

std::string_view syntaxKindName(const SyntaxNode& node) noexcept {
    switch (node.data.index()) {
    case 0:
        return "number literal";
    case 1:
        return "string literal";
    case 2:
        return "identifier";
    case 3:
        return "history reference";
    case 4:
        return "array literal";
    case 5:
        return "function call";
    case 6:
        return "cases expression";
    case 7:
        return "group";
    case 8:
        return "unary expression";
    case 9:
        return "binary expression";
    case 10:
        return "postfix expression";
    case 11:
        return "comparison";
    case 12:
        return "assignment";
    case 13:
        return "function signature";
    case 14:
        return "unit application";
    default:
        return "syntax node";
    }
}

} // namespace mmcal::syntax
