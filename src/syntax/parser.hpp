#pragma once

#include "parse_budget.hpp"
#include "syntax_tree.hpp"

#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mmcal::syntax {

struct ParserOptions final {
    std::unordered_set<std::string> constants;
    std::unordered_set<std::string> protectedNames;
    std::unordered_set<std::string> functions;
    std::unordered_set<std::string> units;

    [[nodiscard]] static ParserOptions defaults();
    [[nodiscard]] bool isConstant(std::string_view name) const;
    [[nodiscard]] bool isProtected(std::string_view name) const;
    [[nodiscard]] bool isFunction(std::string_view name) const;
    [[nodiscard]] bool isUnit(std::string_view name) const;
};

// トークン列から構文ASTを構築する。数値化や組み込み函数の評価は行わない。
class Parser final {
public:
    Parser(
        std::shared_ptr<const std::string> sourceText,
        std::vector<Token> tokens,
        ParserOptions options = ParserOptions::defaults(),
        ParseBudget* budget = nullptr);

    [[nodiscard]] SyntaxTree parse();

private:
    std::shared_ptr<const std::string> sourceText_;
    std::vector<Token> tokens_;
    ParserOptions options_;
    ParseBudget ownedBudget_;
    ParseBudget* budget_ = nullptr;
    std::size_t index_ = 0;
    bool casesIfDelimiter_ = false;

    [[nodiscard]] const Token& current() const noexcept;
    [[nodiscard]] const Token& previous() const noexcept;
    [[nodiscard]] bool check(TokenKind kind) const noexcept;
    bool match(TokenKind kind) noexcept;
    const Token& consume(TokenKind kind, std::string_view message);
    [[nodiscard]] std::string_view tokenText(const Token& token) const noexcept;

    [[nodiscard]] SyntaxNodePtr parseAssignment();
    [[nodiscard]] SyntaxNodePtr parseComparison();
    [[nodiscard]] SyntaxNodePtr parseExpression();
    [[nodiscard]] SyntaxNodePtr parseTerm();
    [[nodiscard]] SyntaxNodePtr parseUnary();
    [[nodiscard]] SyntaxNodePtr parsePostfix();
    [[nodiscard]] SyntaxNodePtr parsePrimary();
    [[nodiscard]] SyntaxNodePtr parseArray();
    [[nodiscard]] SyntaxNodePtr parseHistoryReference();
    [[nodiscard]] SyntaxNodePtr parseInputHistoryReference();
    [[nodiscard]] SyntaxNodePtr parseIdentifierOrCall();
    [[nodiscard]] SyntaxNodePtr parseCasesCall(const Token& identifier);
    [[nodiscard]] SyntaxNodePtr parseGroup();

    [[nodiscard]] bool canStartImplicitFactor() const noexcept;
    [[nodiscard]] bool isConstantIdentifier(const SyntaxNodePtr& node) const;
    [[nodiscard]] SyntaxNodePtr makeAssignmentTarget(const SyntaxNodePtr& node);
    [[nodiscard]] static ComparisonOperator comparisonOperator(TokenKind kind);

    template <class NodeData>
    [[nodiscard]] SyntaxNodePtr makeNode(source::SourceSpan span, NodeData data) {
        budget_->consumeNode(span);
        return makeSyntaxNode(span, std::move(data));
    }
};

} // namespace mmcal::syntax
