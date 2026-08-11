// 構文解析perser
#include "parser.hpp"

#include "error/error_message.hpp"
#include "symbols/symbol_registry.hpp"

#include <algorithm>
#include <string>
#include <unordered_set>
#include <utility>

namespace mmcal::syntax {
namespace {

[[nodiscard]] source::SourceSpan combinedSpan(
    const SyntaxNodePtr& left,
    const SyntaxNodePtr& right) noexcept {
    return source::combine(left->span, right->span);
}

} // namespace

ParserOptions ParserOptions::defaults() {
    const auto constants = symbols::defaultSymbolRegistry().sourcePredefinedNames();
    return ParserOptions{
        constants,
        constants,
        {"Deg", "Rad", "Grad", "deg", "rad", "grad", "mm", "cm", "m", "inch"}
    };
}

bool ParserOptions::isConstant(std::string_view name) const {
    return constants.contains(std::string{name});
}

bool ParserOptions::isProtected(std::string_view name) const {
    return protectedNames.contains(std::string{name});
}

bool ParserOptions::isUnit(std::string_view name) const {
    return units.contains(std::string{name});
}

Parser::Parser(
    std::shared_ptr<const std::string> sourceText,
    std::vector<Token> tokens,
    ParserOptions options)
    : sourceText_(std::move(sourceText)),
      tokens_(std::move(tokens)),
      options_(std::move(options)) {
    if (tokens_.empty() || tokens_.back().kind != TokenKind::End)
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "Parser requires an end-of-input token");
}

SyntaxTree Parser::parse() {
    SyntaxNodePtr root = parseAssignment();
    consume(TokenKind::End, "Unexpected token after the expression");
    return SyntaxTree{sourceText_, std::move(tokens_), std::move(root)};
}

const Token& Parser::current() const noexcept {
    return tokens_[index_];
}

const Token& Parser::previous() const noexcept {
    return tokens_[index_ - 1];
}

bool Parser::check(TokenKind kind) const noexcept {
    return current().kind == kind;
}

bool Parser::match(TokenKind kind) noexcept {
    if (!check(kind))
        return false;

    ++index_;
    return true;
}

const Token& Parser::consume(TokenKind kind, std::string_view message) {
    if (check(kind)) {
        ++index_;
        return previous();
    }

    error::throwCalcError(
        error::CalcErrorType::Syntax,
        std::string{message} + "; found " + std::string{tokenKindName(current().kind)},
        current().span);
}

std::string_view Parser::tokenText(const Token& token) const noexcept {
    return token.text(*sourceText_);
}

SyntaxNodePtr Parser::parseAssignment() {
    SyntaxNodePtr left = parseComparison();
    if (!match(TokenKind::Assign))
        return left;

    SyntaxNodePtr target = makeAssignmentTarget(left);
    SyntaxNodePtr value = parseAssignment();
    const source::SourceSpan span{target->span.begin, value->span.end};
    return makeSyntaxNode(
        span,
        AssignmentSyntax{std::move(target), std::move(value)});
}

SyntaxNodePtr Parser::parseComparison() {
    std::vector<SyntaxNodePtr> operands;
    std::vector<ComparisonOperator> operations;
    operands.push_back(parseExpression());

    while (check(TokenKind::Less)
        || check(TokenKind::LessEqual)
        || check(TokenKind::Greater)
        || check(TokenKind::GreaterEqual)
        || check(TokenKind::EqualEqual)
        || check(TokenKind::BangEqual)) {
        const TokenKind operation = current().kind;
        ++index_;
        operations.push_back(comparisonOperator(operation));
        operands.push_back(parseExpression());
    }

    if (operations.empty())
        return operands.front();

    const source::SourceSpan span{operands.front()->span.begin, operands.back()->span.end};
    return makeSyntaxNode(
        span,
        ComparisonSyntax{std::move(operands), std::move(operations)});
}

SyntaxNodePtr Parser::parseExpression() {
    SyntaxNodePtr left = parseTerm();

    while (check(TokenKind::Plus) || check(TokenKind::Minus)) {
        const TokenKind operation = current().kind;
        ++index_;
        SyntaxNodePtr right = parseTerm();
        const source::SourceSpan span = combinedSpan(left, right);
        left = makeSyntaxNode(
            span,
            BinarySyntax{
                operation == TokenKind::Plus
                    ? BinaryOperator::Add
                    : BinaryOperator::Subtract,
                std::move(left),
                std::move(right)});
    }

    return left;
}

SyntaxNodePtr Parser::parseTerm() {
    SyntaxNodePtr left = parseUnary();

    while (true) {
        if (match(TokenKind::Star)) {
            SyntaxNodePtr right = parseUnary();
            const source::SourceSpan span = combinedSpan(left, right);
            left = makeSyntaxNode(
                span,
                BinarySyntax{BinaryOperator::Multiply, std::move(left), std::move(right)});
            continue;
        }

        if (match(TokenKind::Slash)) {
            SyntaxNodePtr right = parseUnary();
            const source::SourceSpan span = combinedSpan(left, right);
            left = makeSyntaxNode(
                span,
                BinarySyntax{BinaryOperator::Divide, std::move(left), std::move(right)});
            continue;
        }

        if (!canStartImplicitFactor())
            break;

        if (check(TokenKind::Identifier) && options_.isUnit(tokenText(current()))) {
            const Token unit = current();
            ++index_;
            const source::SourceSpan span{left->span.begin, unit.span.end};
            left = makeSyntaxNode(
                span,
                UnitAppliedSyntax{std::move(left), std::string{tokenText(unit)}});
            continue;
        }

        if (std::holds_alternative<IdentifierSyntax>(left->data)
            && check(TokenKind::Number)
            && !isConstantIdentifier(left))
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "An identifier followed by a number is not implicit multiplication",
                current().span);

        SyntaxNodePtr right = parseUnary();
        const source::SourceSpan span = combinedSpan(left, right);
        left = makeSyntaxNode(
            span,
            BinarySyntax{
                BinaryOperator::ImplicitMultiply,
                std::move(left),
                std::move(right)});
    }

    return left;
}

SyntaxNodePtr Parser::parseUnary() {
    if (match(TokenKind::Plus)) {
        const Token operation = previous();
        SyntaxNodePtr operand = parseUnary();
        const source::SourceSpan span{operation.span.begin, operand->span.end};
        return makeSyntaxNode(
            span,
            UnarySyntax{UnaryOperator::Plus, std::move(operand)});
    }

    if (match(TokenKind::Minus)) {
        const Token operation = previous();
        SyntaxNodePtr operand = parseUnary();
        const source::SourceSpan span{operation.span.begin, operand->span.end};
        return makeSyntaxNode(
            span,
            UnarySyntax{UnaryOperator::Minus, std::move(operand)});
    }

    return parsePower();
}

SyntaxNodePtr Parser::parsePower() {
    SyntaxNodePtr left = parsePostfix();
    if (!match(TokenKind::Caret))
        return left;

    SyntaxNodePtr right = parseUnary();
    const source::SourceSpan span = combinedSpan(left, right);
    return makeSyntaxNode(
        span,
        BinarySyntax{BinaryOperator::Power, std::move(left), std::move(right)});
}

SyntaxNodePtr Parser::parsePostfix() {
    SyntaxNodePtr value = parsePrimary();

    while (match(TokenKind::Bang)) {
        const Token operation = previous();
        const source::SourceSpan span{value->span.begin, operation.span.end};
        value = makeSyntaxNode(
            span,
            PostfixSyntax{PostfixOperator::Factorial, std::move(value)});
    }

    return value;
}

SyntaxNodePtr Parser::parsePrimary() {
    if (match(TokenKind::Number)) {
        const Token token = previous();
        return makeSyntaxNode(
            token.span,
            NumberLiteralSyntax{std::string{tokenText(token)}});
    }

    if (match(TokenKind::String)) {
        const Token token = previous();
        return makeSyntaxNode(
            token.span,
            StringLiteralSyntax{std::string{tokenText(token)}});
    }

    if (check(TokenKind::Percent))
        return parseHistoryReference();
    if (check(TokenKind::LBrace))
        return parseArray();
    if (check(TokenKind::Identifier))
        return parseIdentifierOrCall();
    if (check(TokenKind::LParen))
        return parseGroup();
    if (check(TokenKind::LBracket))
        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "Square brackets are only valid for function calls",
            current().span);

    error::throwCalcError(
        error::CalcErrorType::Syntax,
        "Expected a value",
        current().span);
}

SyntaxNodePtr Parser::parseArray() {
    const Token opening = consume(TokenKind::LBrace, "Expected '{'");
    std::vector<SyntaxNodePtr> elements;

    if (!check(TokenKind::RBrace)) {
        do {
            elements.push_back(parseAssignment());
        } while (match(TokenKind::Comma));
    }

    const Token closing = consume(TokenKind::RBrace, "Expected '}' after array elements");
    return makeSyntaxNode(
        source::SourceSpan{opening.span.begin, closing.span.end},
        ArrayLiteralSyntax{std::move(elements)});
}

SyntaxNodePtr Parser::parseHistoryReference() {
    const Token first = consume(TokenKind::Percent, "Expected '%'");
    std::size_t depth = 1;
    while (match(TokenKind::Percent))
        ++depth;

    return makeSyntaxNode(
        source::SourceSpan{first.span.begin, previous().span.end},
        HistoryReferenceSyntax{depth});
}

SyntaxNodePtr Parser::parseIdentifierOrCall() {
    const Token identifier = consume(TokenKind::Identifier, "Expected identifier");
    const std::string name{tokenText(identifier)};

    if (!check(TokenKind::LParen) && !check(TokenKind::LBracket))
        return makeSyntaxNode(identifier.span, IdentifierSyntax{name});

    const bool parentheses = match(TokenKind::LParen);
    if (!parentheses)
        consume(TokenKind::LBracket, "Expected '['");

    const TokenKind closingKind = parentheses ? TokenKind::RParen : TokenKind::RBracket;
    const std::string_view closingText = parentheses ? ")" : "]";
    std::vector<SyntaxNodePtr> arguments;

    if (!check(closingKind)) {
        do {
            arguments.push_back(parseAssignment());
        } while (match(TokenKind::Comma));
    }

    const Token closing = consume(
        closingKind,
        std::string{"Expected '"} + std::string{closingText} + "' after function arguments");

    return makeSyntaxNode(
        source::SourceSpan{identifier.span.begin, closing.span.end},
        CallSyntax{
            name,
            std::move(arguments),
            parentheses ? CallDelimiter::Parentheses : CallDelimiter::Brackets});
}

SyntaxNodePtr Parser::parseGroup() {
    const Token opening = consume(TokenKind::LParen, "Expected '('");
    SyntaxNodePtr expression = parseAssignment();
    const Token closing = consume(TokenKind::RParen, "Expected ')' after grouped expression");

    return makeSyntaxNode(
        source::SourceSpan{opening.span.begin, closing.span.end},
        GroupSyntax{std::move(expression)});
}

bool Parser::canStartImplicitFactor() const noexcept {
    switch (current().kind) {
    case TokenKind::Number:
    case TokenKind::Identifier:
    case TokenKind::LParen:
    case TokenKind::Percent:
        return true;
    default:
        return false;
    }
}

bool Parser::isConstantIdentifier(const SyntaxNodePtr& node) const {
    const auto* identifier = std::get_if<IdentifierSyntax>(&node->data);
    return identifier && options_.isConstant(identifier->name);
}

SyntaxNodePtr Parser::makeAssignmentTarget(const SyntaxNodePtr& node) const {
    if (const auto* identifier = std::get_if<IdentifierSyntax>(&node->data)) {
        if (options_.isProtected(identifier->name))
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "Protected names cannot be assignment targets",
                node->span);

        return node;
    }

    const auto* call = std::get_if<CallSyntax>(&node->data);
    if (!call)
        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "Assignment target must be an identifier or function signature",
            node->span);

    if (options_.isProtected(call->name))
        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "Protected names cannot be function names",
            node->span);

    std::vector<std::string> parameters;
    parameters.reserve(call->arguments.size());
    std::unordered_set<std::string> uniqueParameters;

    for (const SyntaxNodePtr& argument : call->arguments) {
        const auto* parameter = std::get_if<IdentifierSyntax>(&argument->data);
        if (!parameter)
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "Function parameters must be identifiers",
                argument->span);

        if (options_.isProtected(parameter->name))
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "Protected names cannot be function parameters",
                argument->span);

        if (!uniqueParameters.insert(parameter->name).second)
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "Function parameters must be unique",
                argument->span);

        parameters.push_back(parameter->name);
    }

    return makeSyntaxNode(
        node->span,
        FunctionSignatureSyntax{call->name, std::move(parameters), call->delimiter});
}

ComparisonOperator Parser::comparisonOperator(TokenKind kind) {
    switch (kind) {
    case TokenKind::Less:
        return ComparisonOperator::Less;
    case TokenKind::LessEqual:
        return ComparisonOperator::LessEqual;
    case TokenKind::Greater:
        return ComparisonOperator::Greater;
    case TokenKind::GreaterEqual:
        return ComparisonOperator::GreaterEqual;
    case TokenKind::EqualEqual:
        return ComparisonOperator::Equal;
    case TokenKind::BangEqual:
        return ComparisonOperator::NotEqual;
    default:
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "Token is not a comparison operator");
    }
}

} // namespace mmcal::syntax
