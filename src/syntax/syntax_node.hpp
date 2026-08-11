#pragma once

#include "source/source_span.hpp"
#include "token.hpp"

#include <memory>
#include <string>
#include <variant>
#include <vector>

namespace mmcal::syntax {

struct SyntaxNode;
using SyntaxNodePtr = std::shared_ptr<const SyntaxNode>;

enum class UnaryOperator {
    Plus,
    Minus
};

enum class BinaryOperator {
    Add,
    Subtract,
    Multiply,
    Divide,
    Power,
    ImplicitMultiply
};

enum class PostfixOperator {
    Factorial
};

enum class ComparisonOperator {
    Less,
    LessEqual,
    Greater,
    GreaterEqual,
    Equal,
    NotEqual
};

enum class CallDelimiter {
    Parentheses,
    Brackets
};

struct NumberLiteralSyntax final {
    std::string text;
};

struct StringLiteralSyntax final {
    std::string text;
};

struct IdentifierSyntax final {
    std::string name;
};

struct HistoryReferenceSyntax final {
    std::size_t depth = 1;
};

struct ArrayLiteralSyntax final {
    std::vector<SyntaxNodePtr> elements;
};

struct CallSyntax final {
    std::string name;
    std::vector<SyntaxNodePtr> arguments;
    CallDelimiter delimiter = CallDelimiter::Parentheses;
};

struct GroupSyntax final {
    SyntaxNodePtr expression;
};

struct UnarySyntax final {
    UnaryOperator operation = UnaryOperator::Plus;
    SyntaxNodePtr operand;
};

struct BinarySyntax final {
    BinaryOperator operation = BinaryOperator::Add;
    SyntaxNodePtr left;
    SyntaxNodePtr right;
};

struct PostfixSyntax final {
    PostfixOperator operation = PostfixOperator::Factorial;
    SyntaxNodePtr operand;
};

struct ComparisonSyntax final {
    std::vector<SyntaxNodePtr> operands;
    std::vector<ComparisonOperator> operations;
};

struct AssignmentSyntax final {
    SyntaxNodePtr target;
    SyntaxNodePtr value;
};

struct FunctionSignatureSyntax final {
    std::string name;
    std::vector<std::string> parameters;
    CallDelimiter delimiter = CallDelimiter::Parentheses;
};

struct UnitAppliedSyntax final {
    SyntaxNodePtr value;
    std::string unit;
};

using SyntaxNodeData = std::variant<
    NumberLiteralSyntax,
    StringLiteralSyntax,
    IdentifierSyntax,
    HistoryReferenceSyntax,
    ArrayLiteralSyntax,
    CallSyntax,
    GroupSyntax,
    UnarySyntax,
    BinarySyntax,
    PostfixSyntax,
    ComparisonSyntax,
    AssignmentSyntax,
    FunctionSignatureSyntax,
    UnitAppliedSyntax>;

struct SyntaxNode final {
    source::SourceSpan span;
    SyntaxNodeData data;
};

template <class NodeData>
[[nodiscard]] SyntaxNodePtr makeSyntaxNode(source::SourceSpan span, NodeData data) {
    return std::make_shared<SyntaxNode>(SyntaxNode{span, SyntaxNodeData{std::move(data)}});
}

[[nodiscard]] std::string_view syntaxKindName(const SyntaxNode& node) noexcept;

} // namespace mmcal::syntax
