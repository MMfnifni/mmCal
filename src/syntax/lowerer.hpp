#pragma once

#include "expression/expr.hpp"
#include "expression/origin_map.hpp"
#include "syntax_tree.hpp"
#include "source/source_document.hpp"
#include "symbols/symbol_registry.hpp"
#include "symbols/symbol_table.hpp"

#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>

namespace mmcal::syntax {

struct LoweringOptions final {
    std::unordered_set<std::string> constants;
    std::unordered_set<std::string> functions;
    std::unordered_set<std::string> variables;

    [[nodiscard]] static LoweringOptions defaults();
    [[nodiscard]] bool isConstant(std::string_view name) const;
    [[nodiscard]] bool isFunction(std::string_view name) const;
    [[nodiscard]] bool isVariable(std::string_view name) const;
};

struct LoweringResult final {
    expression::Expr expression;
    expression::OriginMap origins;
};

// 構文ASTを評価器用Exprへ変換し、構文と意味表現の責務を分離する。
class Lowerer final {
public:
    explicit Lowerer(
        LoweringOptions options = LoweringOptions::defaults(),
        const symbols::SymbolRegistry& symbolRegistry = symbols::defaultSymbolRegistry(),
        symbols::SymbolTable& symbolTable = symbols::defaultSymbolTable());

    [[nodiscard]] expression::Expr lower(const SyntaxTree& tree) const;
    [[nodiscard]] LoweringResult lowerTracked(const SyntaxTree& tree) const;
    [[nodiscard]] LoweringResult lowerTracked(
        const SyntaxTree& tree,
        std::shared_ptr<const source::SourceDocument> document) const;
    [[nodiscard]] const LoweringOptions& options() const noexcept;
    [[nodiscard]] LoweringOptions& options() noexcept;

private:
    LoweringOptions options_;
    const symbols::SymbolRegistry& symbolRegistry_;
    symbols::SymbolTable& symbolTable_;

    [[nodiscard]] expression::Expr lowerNode(
        const SyntaxNode& node,
        expression::OriginMap* origins) const;
    [[nodiscard]] expression::Expr lowerNumber(
        const NumberLiteralSyntax& number,
        source::SourceSpan span) const;
    [[nodiscard]] expression::Expr lowerString(
        const StringLiteralSyntax& string,
        source::SourceSpan span) const;
    [[nodiscard]] expression::Expr lowerIdentifier(const IdentifierSyntax& identifier) const;
    [[nodiscard]] expression::Expr lowerArray(
        const ArrayLiteralSyntax& array,
        source::SourceSpan span,
        expression::OriginMap* origins) const;
    [[nodiscard]] expression::Expr lowerCall(
        const CallSyntax& call,
        expression::OriginMap* origins) const;
    [[nodiscard]] expression::Expr lowerComparison(
        const ComparisonSyntax& comparison,
        expression::OriginMap* origins) const;
    [[nodiscard]] expression::Expr lowerAssignment(
        const AssignmentSyntax& assignment,
        expression::OriginMap* origins) const;
};

} // namespace mmcal::syntax
