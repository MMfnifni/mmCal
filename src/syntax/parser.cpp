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
        {},
        {"Deg", "Rad", "Grad", "deg", "rad", "grad", "mm", "cm", "m", "inch"}
    };
}

bool ParserOptions::isConstant(std::string_view name) const {
    return constants.contains(std::string{name});
}

bool ParserOptions::isProtected(std::string_view name) const {
    return protectedNames.contains(std::string{name});
}

bool ParserOptions::isFunction(std::string_view name) const {
    return functions.contains(std::string{name});
}

bool ParserOptions::isUnit(std::string_view name) const {
    return units.contains(std::string{name});
}

Parser::Parser(
    std::shared_ptr<const std::string> sourceText,
    std::vector<Token> tokens,
    ParserOptions options,
    ParseBudget* budget)
    : sourceText_(std::move(sourceText)),
      tokens_(std::move(tokens)),
      options_(std::move(options)),
      budget_(budget ? budget : &ownedBudget_) {
    if (tokens_.empty() || tokens_.back().kind != TokenKind::End)
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "Parser requires an end-of-input token");

    budget_->checkTokenCount(tokens_.size(), tokens_.back().span);
}

SyntaxTree Parser::parse() {
    const auto depthGuard = budget_->enter(current().span);
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
    const auto depthGuard = budget_->enter(current().span);
    std::vector<SyntaxNodePtr> targets;
    SyntaxNodePtr value = parseComparison();
    std::size_t operationCount = 0;

    while (match(TokenKind::Assign)) {
        budget_->checkOperatorChain(++operationCount, previous().span);
        targets.push_back(makeAssignmentTarget(value));
        value = parseComparison();
    }

    // v1.5.3では右辺をparseAssignment()で再帰していた。ここでは同じ右結合を
    // 末尾から組み立て，長い a:=b:=... がC++スタックを消費しないようにする。
    for (auto iterator = targets.rbegin(); iterator != targets.rend(); ++iterator) {
        SyntaxNodePtr target = *iterator;
        const source::SourceSpan span{target->span.begin, value->span.end};
        value = makeNode(
            span,
            AssignmentSyntax{std::move(target), std::move(value)});
    }

    return value;
}

SyntaxNodePtr Parser::parseComparison() {
    const auto depthGuard = budget_->enter(current().span);
    std::vector<SyntaxNodePtr> operands;
    std::vector<ComparisonOperator> operations;
    std::size_t operationCount = 0;
    operands.push_back(parseExpression());

    while (check(TokenKind::Less)
        || check(TokenKind::LessEqual)
        || check(TokenKind::Greater)
        || check(TokenKind::GreaterEqual)
        || check(TokenKind::EqualEqual)
        || check(TokenKind::BangEqual)) {
        const TokenKind operation = current().kind;
        const source::SourceSpan operationSpan = current().span;
        ++index_;
        budget_->checkOperatorChain(++operationCount, operationSpan);
        operations.push_back(comparisonOperator(operation));
        operands.push_back(parseExpression());
    }

    if (operations.empty())
        return operands.front();

    const source::SourceSpan span{operands.front()->span.begin, operands.back()->span.end};
    return makeNode(
        span,
        ComparisonSyntax{std::move(operands), std::move(operations)});
}

SyntaxNodePtr Parser::parseExpression() {
    const auto depthGuard = budget_->enter(current().span);
    SyntaxNodePtr left = parseTerm();
    std::size_t operationCount = 0;

    while (check(TokenKind::Plus) || check(TokenKind::Minus)) {
        const TokenKind operation = current().kind;
        const source::SourceSpan operationSpan = current().span;
        ++index_;
        budget_->checkOperatorChain(++operationCount, operationSpan);
        SyntaxNodePtr right = parseTerm();
        const source::SourceSpan span = combinedSpan(left, right);
        left = makeNode(
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
    const auto depthGuard = budget_->enter(current().span);
    SyntaxNodePtr left = parseUnary();
    std::size_t operationCount = 0;

    while (true) {
        if (match(TokenKind::Star)) {
            budget_->checkOperatorChain(++operationCount, previous().span);
            SyntaxNodePtr right = parseUnary();
            const source::SourceSpan span = combinedSpan(left, right);
            left = makeNode(
                span,
                BinarySyntax{BinaryOperator::Multiply, std::move(left), std::move(right)});
            continue;
        }

        if (match(TokenKind::Slash)) {
            budget_->checkOperatorChain(++operationCount, previous().span);
            SyntaxNodePtr right = parseUnary();
            const source::SourceSpan span = combinedSpan(left, right);
            left = makeNode(
                span,
                BinarySyntax{BinaryOperator::Divide, std::move(left), std::move(right)});
            continue;
        }

        if (!canStartImplicitFactor())
            break;

        if (check(TokenKind::Identifier) && options_.isUnit(tokenText(current()))) {
            const Token unit = current();
            ++index_;
            budget_->checkOperatorChain(++operationCount, unit.span);
            const source::SourceSpan span{left->span.begin, unit.span.end};
            left = makeNode(
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

        const source::SourceSpan operationSpan = current().span;
        budget_->checkOperatorChain(++operationCount, operationSpan);
        SyntaxNodePtr right = parseUnary();
        const source::SourceSpan span = combinedSpan(left, right);
        left = makeNode(
            span,
            BinarySyntax{
                BinaryOperator::ImplicitMultiply,
                std::move(left),
                std::move(right)});
    }

    return left;
}

SyntaxNodePtr Parser::parseUnary() {
    const auto depthGuard = budget_->enter(current().span);

    struct PowerOperand final {
        std::vector<Token> prefixes;
        SyntaxNodePtr base;
    };

    std::vector<PowerOperand> operands;
    std::size_t operationCount = 0;

    while (true) {
        PowerOperand operand;
        while (check(TokenKind::Plus) || check(TokenKind::Minus)) {
            operand.prefixes.push_back(current());
            budget_->checkOperatorChain(++operationCount, current().span);
            ++index_;
        }

        operand.base = parsePostfix();
        operands.push_back(std::move(operand));

        if (!match(TokenKind::Caret))
            break;
        budget_->checkOperatorChain(++operationCount, previous().span);
    }

    const auto applyPrefixes = [this](PowerOperand& operand, SyntaxNodePtr value) {
        for (auto iterator = operand.prefixes.rbegin();
             iterator != operand.prefixes.rend();
             ++iterator) {
            const source::SourceSpan span{iterator->span.begin, value->span.end};
            value = makeNode(
                span,
                UnarySyntax{
                    iterator->kind == TokenKind::Plus
                        ? UnaryOperator::Plus
                        : UnaryOperator::Minus,
                    std::move(value)});
        }
        return value;
    };

    SyntaxNodePtr value = applyPrefixes(operands.back(), operands.back().base);
    for (std::size_t index = operands.size() - 1; index != 0; --index) {
        PowerOperand& operand = operands[index - 1];
        const source::SourceSpan span = combinedSpan(operand.base, value);
        value = makeNode(
            span,
            BinarySyntax{BinaryOperator::Power, operand.base, std::move(value)});
        value = applyPrefixes(operand, std::move(value));
    }

    // v1.5.3のparseUnary()/parsePower()相互再帰と同じ優先順位・右結合を保つ。
    // 配列へ平坦化したのは，++++1 や 2^2^... でC++スタックを使わないためである。
    return value;
}

SyntaxNodePtr Parser::parsePostfix() {
    const auto depthGuard = budget_->enter(current().span);
    SyntaxNodePtr value = parsePrimary();
    std::size_t operationCount = 0;

    while (match(TokenKind::Bang)) {
        const Token operation = previous();
        budget_->checkOperatorChain(++operationCount, operation.span);
        const source::SourceSpan span{value->span.begin, operation.span.end};
        value = makeNode(
            span,
            PostfixSyntax{PostfixOperator::Factorial, std::move(value)});
    }

    return value;
}

SyntaxNodePtr Parser::parsePrimary() {
    const auto depthGuard = budget_->enter(current().span);
    if (match(TokenKind::Number)) {
        const Token token = previous();
        return makeNode(
            token.span,
            NumberLiteralSyntax{std::string{tokenText(token)}});
    }

    if (match(TokenKind::String)) {
        const Token token = previous();
        return makeNode(
            token.span,
            StringLiteralSyntax{std::string{tokenText(token)}});
    }

    if (check(TokenKind::Percent))
        return parseHistoryReference();
    if (check(TokenKind::At))
        return parseInputHistoryReference();
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
    const auto depthGuard = budget_->enter(current().span);
    const Token opening = consume(TokenKind::LBrace, "Expected '{'");
    std::vector<SyntaxNodePtr> elements;

    if (!check(TokenKind::RBrace)) {
        do {
            elements.push_back(parseAssignment());
            budget_->checkArrayElements(elements.size(), elements.back()->span);
        } while (match(TokenKind::Comma));
    }

    const Token closing = consume(TokenKind::RBrace, "Expected '}' after array elements");
    return makeNode(
        source::SourceSpan{opening.span.begin, closing.span.end},
        ArrayLiteralSyntax{std::move(elements)});
}

SyntaxNodePtr Parser::parseHistoryReference() {
    const auto depthGuard = budget_->enter(current().span);
    const Token first = consume(TokenKind::Percent, "Expected '%'");
    std::size_t depth = 1;
    while (match(TokenKind::Percent))
        ++depth;

    return makeNode(
        source::SourceSpan{first.span.begin, previous().span.end},
        HistoryReferenceSyntax{HistoryReferenceKind::Output, depth});
}

SyntaxNodePtr Parser::parseInputHistoryReference() {
    const auto depthGuard = budget_->enter(current().span);
    const Token first = consume(TokenKind::At, "Expected '@'");
    std::size_t depth = 1;
    while (match(TokenKind::At))
        ++depth;

    return makeNode(
        source::SourceSpan{first.span.begin, previous().span.end},
        HistoryReferenceSyntax{HistoryReferenceKind::Input, depth});
}

SyntaxNodePtr Parser::parseIdentifierOrCall() {
    const auto depthGuard = budget_->enter(current().span);
    const Token identifier = consume(TokenKind::Identifier, "Expected identifier");
    const std::string name{tokenText(identifier)};

    // v1.5.2以降、函数呼出しは name[...] に限定する。
    // name(expr) は通常のidentifierなら暗黙乗算としてparseTermへ返すが、
    // 既知の函数名については旧構文の打ち間違いを黙って乗算へ変えない。
    if (check(TokenKind::LParen) && options_.isFunction(name))
        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "Function calls require square brackets; use " + name + "[...]",
            current().span);

    if (!check(TokenKind::LBracket))
        return makeNode(identifier.span, IdentifierSyntax{name});

    if (name == "cases")
        return parseCasesCall(identifier);

    consume(TokenKind::LBracket, "Expected '['");
    std::vector<SyntaxNodePtr> arguments;

    if (!check(TokenKind::RBracket)) {
        do {
            arguments.push_back(parseAssignment());
            budget_->checkCallArguments(arguments.size(), arguments.back()->span);
        } while (match(TokenKind::Comma));
    }

    const Token closing = consume(
        TokenKind::RBracket,
        "Expected ']' after function arguments");

    return makeNode(
        source::SourceSpan{identifier.span.begin, closing.span.end},
        CallSyntax{name, std::move(arguments)});
}

SyntaxNodePtr Parser::parseCasesCall(const Token& identifier) {
    const auto depthGuard = budget_->enter(current().span);
    consume(TokenKind::LBracket, "Expected '[' after cases");
    std::vector<CasesBranchSyntax> branches;
    bool defaultSeen = false;

    if (check(TokenKind::RBracket))
        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "cases requires at least one conditional branch",
            current().span);

    while (true) {
        const bool previousDelimiter = casesIfDelimiter_;
        casesIfDelimiter_ = true;
        SyntaxNodePtr value = parseAssignment();
        casesIfDelimiter_ = previousDelimiter;

        SyntaxNodePtr condition;
        if (check(TokenKind::Identifier) && tokenText(current()) == "if") {
            if (defaultSeen)
                error::throwCalcError(
                    error::CalcErrorType::Syntax,
                    "cases default branch must be last",
                    current().span);
            ++index_;
            condition = parseAssignment();
        }
        else {
            defaultSeen = true;
        }

        branches.push_back(CasesBranchSyntax{std::move(value), std::move(condition)});
        budget_->checkCallArguments(branches.size(), branches.back().value->span);

        if (!match(TokenKind::Semicolon))
            break;
        if (defaultSeen)
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "cases default branch must be last",
                previous().span);
        if (check(TokenKind::RBracket))
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "Expected cases branch after ';'",
                current().span);
    }

    const Token closing = consume(TokenKind::RBracket, "Expected ']' after cases branches");
    bool hasConditional = false;
    for (const CasesBranchSyntax& branch : branches)
        hasConditional = hasConditional || static_cast<bool>(branch.condition);
    if (!hasConditional)
        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "cases requires at least one 'value if condition' branch",
            source::SourceSpan{identifier.span.begin, closing.span.end});

    return makeNode(
        source::SourceSpan{identifier.span.begin, closing.span.end},
        CasesSyntax{std::move(branches)});
}

SyntaxNodePtr Parser::parseGroup() {
    const auto depthGuard = budget_->enter(current().span);
    const Token opening = consume(TokenKind::LParen, "Expected '('");
    SyntaxNodePtr expression = parseAssignment();
    const Token closing = consume(TokenKind::RParen, "Expected ')' after grouped expression");

    return makeNode(
        source::SourceSpan{opening.span.begin, closing.span.end},
        GroupSyntax{std::move(expression)});
}

bool Parser::canStartImplicitFactor() const noexcept {
    if (casesIfDelimiter_ && check(TokenKind::Identifier) && tokenText(current()) == "if")
        return false;

    switch (current().kind) {
    case TokenKind::Number:
    case TokenKind::Identifier:
    case TokenKind::LParen:
    case TokenKind::Percent:
    case TokenKind::At:
        return true;
    default:
        return false;
    }
}

bool Parser::isConstantIdentifier(const SyntaxNodePtr& node) const {
    const auto* identifier = std::get_if<IdentifierSyntax>(&node->data);
    return identifier && options_.isConstant(identifier->name);
}

SyntaxNodePtr Parser::makeAssignmentTarget(const SyntaxNodePtr& node) {
    if (const auto* identifier = std::get_if<IdentifierSyntax>(&node->data)) {
        if (options_.isProtected(identifier->name))
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "Protected names cannot be assignment targets",
                node->span);

        return node;
    }

    const auto* call = std::get_if<CallSyntax>(&node->data);
    if (!call) {
        // name(expr) := ... はv1.5.2以前の函数定義構文に見えるため、
        // 一般的なassignment target errorより移行方法を直接示す。
        if (const auto* multiply = std::get_if<BinarySyntax>(&node->data);
            multiply && multiply->operation == BinaryOperator::ImplicitMultiply) {
            const auto* name = std::get_if<IdentifierSyntax>(&multiply->left->data);
            if (name && std::holds_alternative<GroupSyntax>(multiply->right->data))
                error::throwCalcError(
                    error::CalcErrorType::Syntax,
                    "Function definitions require square brackets; use "
                        + name->name + "[...] := ...",
                    node->span);
        }

        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "Assignment target must be an identifier or function signature",
            node->span);
    }

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

    return makeNode(
        node->span,
        FunctionSignatureSyntax{call->name, std::move(parameters)});
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
