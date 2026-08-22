#pragma once

#include "evaluation/builtin_registry.hpp"
#include "evaluation/environment.hpp"
#include "evaluation/evaluator.hpp"
#include "evaluation/diagnostic.hpp"
#include "evaluation/user_function_registry.hpp"
#include "expression/expr.hpp"
#include "syntax/lowerer.hpp"
#include "symbols/symbol_registry.hpp"
#include "symbols/symbol_table.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace mmcal::kernel {

// 一連の入力で共有する定義、履歴、Frontend、評価器を所有する実行単位。
class KernelSession final {
public:
    KernelSession();

    [[nodiscard]] expression::Expr evaluate(std::string_view sourceText);
    [[nodiscard]] expression::Expr evaluate(
        std::string_view sourceText,
        const evaluation::EvaluationCancellationToken& cancellation);

    [[nodiscard]] const expression::Expr* history(std::size_t depth = 1) const noexcept;
    [[nodiscard]] const expression::Expr* inputHistory(std::size_t index) const noexcept;
    [[nodiscard]] const expression::Expr* outputHistory(std::size_t index) const noexcept;
    [[nodiscard]] std::span<const evaluation::EvaluationDiagnostic> diagnostics() const noexcept;
    [[nodiscard]] std::size_t historySize() const noexcept;
    [[nodiscard]] std::size_t inputCount() const noexcept;
    [[nodiscard]] std::size_t nextInputNumber() const noexcept;
    [[nodiscard]] bool exitRequested() const noexcept;
    [[nodiscard]] bool clearRequested() const noexcept;
    void clearHistory() noexcept;
    void clearDefinitions();
    // 独立評価用に定義・履歴・入力番号等だけを戻し，角度設定とRNG streamは保持する。
    void resetForIndependentEvaluation();
    void reset();

    [[nodiscard]] const evaluation::Environment& environment() const noexcept;
    [[nodiscard]] const evaluation::BuiltinRegistry& builtinRegistry() const noexcept;
    [[nodiscard]] const symbols::SymbolTable& symbolTable() const noexcept;
    [[nodiscard]] const symbols::SymbolRegistry& symbolRegistry() const noexcept;
    [[nodiscard]] const mathematics::MathRegistry& mathRegistry() const noexcept;
    [[nodiscard]] mathematics::AngleUnit defaultAngleUnit() const noexcept;
    void setDefaultAngleUnit(mathematics::AngleUnit unit) noexcept;
    [[nodiscard]] const evaluation::UserFunctionRegistry& userFunctions() const noexcept;

    void setEvaluationDepthLimit(std::size_t limit);
    [[nodiscard]] std::size_t evaluationDepthLimit() const noexcept;
    void setEvaluationLimits(evaluation::EvaluationLimits limits);
    [[nodiscard]] const evaluation::EvaluationLimits& evaluationLimits() const noexcept;
    [[nodiscard]] const evaluation::EvaluationUsage& lastEvaluationUsage() const noexcept;

private:
    symbols::SymbolTable symbolTable_;
    symbols::SymbolRegistry symbolRegistry_;
    evaluation::BuiltinRegistry registry_;
    mathematics::MathRegistry mathRegistry_;
    mathematics::AngleSemantics angleSemantics_;
    evaluation::Environment environment_;
    evaluation::UserFunctionRegistry userFunctions_;
    syntax::Lowerer lowerer_;
    evaluation::Evaluator evaluator_;
    std::vector<expression::Expr> history_;
    std::vector<std::optional<expression::Expr>> inputHistory_;
    std::vector<std::optional<expression::Expr>> outputHistory_;
    std::vector<evaluation::EvaluationDiagnostic> diagnostics_;
    std::size_t inputCount_ = 0;
    bool exitRequested_ = false;
    bool clearRequested_ = false;
    bool definitionsChanged_ = false;
    evaluation::EvaluationUsage lastEvaluationUsage_{};

    void rememberSuccessfulDefinition(const syntax::SyntaxTree& tree);
    void resetKnownNames();
    [[nodiscard]] expression::Expr evaluateImpl(
        std::string_view sourceText,
        const evaluation::EvaluationCancellationToken* cancellation);
};

} // namespace mmcal::kernel
