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
        return "group";
    case 7:
        return "unary expression";
    case 8:
        return "binary expression";
    case 9:
        return "postfix expression";
    case 10:
        return "comparison";
    case 11:
        return "assignment";
    case 12:
        return "function signature";
    case 13:
        return "unit application";
    default:
        return "syntax node";
    }
}

} // namespace mmcal::syntax
