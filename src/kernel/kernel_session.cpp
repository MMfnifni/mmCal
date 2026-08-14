// 定義・履歴・角度・診断を持つセッション
#include "kernel_session.hpp"

#include "syntax/lexer.hpp"
#include "syntax/parser.hpp"
#include "error/error_message.hpp"
#include "source/source_document.hpp"

#include <memory>
#include <string>
#include <utility>
#include <variant>

namespace mmcal::kernel {
namespace {

[[nodiscard]] syntax::LoweringOptions makeLoweringOptions(
    const symbols::SymbolRegistry& symbols,
    const evaluation::BuiltinRegistry& registry) {
    syntax::LoweringOptions options;
    options.constants = symbols.sourcePredefinedNames();
    options.functions = registry.sourceFunctionNames();
    return options;
}

[[nodiscard]] syntax::ParserOptions makeParserOptions(
    const symbols::SymbolRegistry& symbols,
    const evaluation::BuiltinRegistry& registry,
    const evaluation::UserFunctionRegistry& userFunctions) {
    syntax::ParserOptions options = syntax::ParserOptions::defaults();
    options.constants = symbols.sourcePredefinedNames();
    options.protectedNames = options.constants;
    for (const std::string& name : registry.sourceFunctionNames()) {
        options.protectedNames.insert(name);
        options.functions.insert(name);
    }
    for (const std::string& name : userFunctions.names())
        options.functions.insert(name);
    return options;
}

} // namespace

KernelSession::KernelSession()
    : symbolRegistry_(symbols::SymbolRegistry::defaults(symbolTable_)),
      registry_(evaluation::BuiltinRegistry::defaults(symbolTable_)),
      mathRegistry_(mathematics::MathRegistry::defaults(symbolTable_, registry_)),
      angleSemantics_(mathematics::AngleUnit::Radian),
      lowerer_(makeLoweringOptions(symbolRegistry_, registry_), symbolRegistry_, symbolTable_),
      evaluator_(
          environment_,
          registry_,
          &userFunctions_,
          symbolRegistry_,
          mathRegistry_,
          angleSemantics_) {}

expression::Expr KernelSession::evaluate(std::string_view sourceText) {
    exitRequested_ = false;
    clearRequested_ = false;
    definitionsChanged_ = false;
    ++inputCount_;
    // clearHistory後も画面の絶対入力番号とvector indexを一致させるため、
    // emplaceではなく現在の入力番号までslotを確保する。
    inputHistory_.resize(inputCount_);
    outputHistory_.resize(inputCount_);
    inputHistory_.back().reset();
    outputHistory_.back().reset();
    diagnostics_.clear();

    auto source = std::make_shared<const std::string>(sourceText);
    auto document = std::make_shared<const source::SourceDocument>(inputCount_, source);

    try {
        syntax::Lexer lexer{*source};
        syntax::Parser parser{
            source,
            lexer.tokenize(),
            makeParserOptions(symbolRegistry_, registry_, userFunctions_)};
        syntax::SyntaxTree tree = parser.parse();
        syntax::LoweringResult lowered = lowerer_.lowerTracked(tree, document);
        inputHistory_.back() = lowered.expression;

        const evaluation::EvaluationContext context{
            history_, inputHistory_, outputHistory_, &diagnostics_,
            &exitRequested_, &clearRequested_, &definitionsChanged_, &angleSemantics_};
        expression::Expr result = evaluator_.evaluate(
            lowered.expression,
            lowered.origins,
            context);

        // Clear[]は現在の入力を履歴へ残さず、履歴番号も1へ戻す。
        // 定義と履歴だけを消し、角度設定やRNG stateは保持する。
        if (clearRequested_) {
            clearDefinitions();
            history_.clear();
            inputHistory_.clear();
            outputHistory_.clear();
            diagnostics_.clear();
            inputCount_ = 0;
            return result;
        }

        if (definitionsChanged_)
            resetKnownNames();
        else
            rememberSuccessfulDefinition(tree);

        history_.push_back(result);
        outputHistory_.back() = result;
        return result;
    }
    catch (error::CalcError& exception) {
        // Lexer/Parser段階のエラーにも入力文書を結び付け、表示側を一貫させる。
        if (exception.span())
            exception.attachDocumentIfMissing(document);
        throw;
    }
}

const expression::Expr* KernelSession::history(std::size_t depth) const noexcept {
    if (depth == 0 || depth > history_.size())
        return nullptr;

    return &history_[history_.size() - depth];
}

const expression::Expr* KernelSession::inputHistory(std::size_t index) const noexcept {
    if (index == 0 || index > inputHistory_.size() || !inputHistory_[index - 1])
        return nullptr;
    return &*inputHistory_[index - 1];
}

const expression::Expr* KernelSession::outputHistory(std::size_t index) const noexcept {
    if (index == 0 || index > outputHistory_.size() || !outputHistory_[index - 1])
        return nullptr;
    return &*outputHistory_[index - 1];
}

std::span<const evaluation::EvaluationDiagnostic> KernelSession::diagnostics() const noexcept {
    return diagnostics_;
}

std::size_t KernelSession::historySize() const noexcept {
    return history_.size();
}

std::size_t KernelSession::inputCount() const noexcept {
    return inputCount_;
}

std::size_t KernelSession::nextInputNumber() const noexcept {
    return inputCount_ + 1;
}

bool KernelSession::exitRequested() const noexcept {
    return exitRequested_;
}

bool KernelSession::clearRequested() const noexcept {
    return clearRequested_;
}

void KernelSession::clearHistory() noexcept {
    history_.clear();
    inputHistory_.clear();
    outputHistory_.clear();
    diagnostics_.clear();
}

void KernelSession::clearDefinitions() {
    environment_.clear();
    userFunctions_.clear();
    resetKnownNames();
}

void KernelSession::resetForIndependentEvaluation() {
    clearDefinitions();
    clearHistory();
    inputCount_ = 0;
    exitRequested_ = false;
    clearRequested_ = false;
    definitionsChanged_ = false;
}

void KernelSession::reset() {
    resetForIndependentEvaluation();
    evaluator_.reseedRandomFromEntropy();
}

const evaluation::Environment& KernelSession::environment() const noexcept {
    return environment_;
}

const evaluation::BuiltinRegistry& KernelSession::builtinRegistry() const noexcept {
    return registry_;
}

const symbols::SymbolTable& KernelSession::symbolTable() const noexcept {
    return symbolTable_;
}

const symbols::SymbolRegistry& KernelSession::symbolRegistry() const noexcept {
    return symbolRegistry_;
}

const mathematics::MathRegistry& KernelSession::mathRegistry() const noexcept {
    return mathRegistry_;
}

mathematics::AngleUnit KernelSession::defaultAngleUnit() const noexcept {
    return angleSemantics_.defaultUnit();
}

void KernelSession::setDefaultAngleUnit(mathematics::AngleUnit unit) noexcept {
    angleSemantics_.setDefaultUnit(unit);
}

const evaluation::UserFunctionRegistry& KernelSession::userFunctions() const noexcept {
    return userFunctions_;
}

void KernelSession::setEvaluationDepthLimit(std::size_t limit) {
    evaluator_.setDepthLimit(limit);
}

std::size_t KernelSession::evaluationDepthLimit() const noexcept {
    return evaluator_.depthLimit();
}

void KernelSession::rememberSuccessfulDefinition(const syntax::SyntaxTree& tree) {
    const auto* assignment = std::get_if<syntax::AssignmentSyntax>(&tree.root().data);
    if (!assignment)
        return;

    if (const auto* identifier = std::get_if<syntax::IdentifierSyntax>(
        &assignment->target->data)) {
        lowerer_.options().variables.insert(identifier->name);
        return;
    }

    if (const auto* signature = std::get_if<syntax::FunctionSignatureSyntax>(
        &assignment->target->data))
        lowerer_.options().functions.insert(signature->name);
}

void KernelSession::resetKnownNames() {
    lowerer_.options().variables.clear();
    for (const auto& [symbol, value] : environment_.definitions()) {
        static_cast<void>(value);
        lowerer_.options().variables.insert(symbol.name());
    }

    lowerer_.options().functions = registry_.sourceFunctionNames();
    if (userFunctions_.size() != 0)
        for (const std::string& name : userFunctions_.names())
            lowerer_.options().functions.insert(name);
}

} // namespace mmcal::kernel
