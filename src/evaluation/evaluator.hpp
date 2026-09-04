#pragma once

#include "builtin_registry.hpp"
#include "approximation/approximation_context.hpp"
#include "builtins/signal_processing.hpp"
#include "environment.hpp"
#include "evaluation_context.hpp"
#include "evaluation_budget.hpp"
#include "expression/origin_map.hpp"
#include "user_function_registry.hpp"
#include "symbols/symbol_registry.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "random/random_engine.hpp"

#include <cstddef>
#include <optional>
#include <string_view>
#include <span>
#include <vector>

namespace mmcal::evaluation {

// 不変な式を入力とし、明示的な評価スタック上で簡約済みの新しい式を返す評価器。
class Evaluator final {
public:
    explicit Evaluator(
        Environment& environment,
        const BuiltinRegistry& registry = defaultBuiltinRegistry(),
        UserFunctionRegistry* userFunctions = nullptr,
        const symbols::SymbolRegistry& symbolRegistry = symbols::defaultSymbolRegistry(),
        const mathematics::MathRegistry& mathematics = mathematics::defaultMathRegistry(),
        const mathematics::AngleSemantics& angleSemantics = mathematics::defaultAngleSemantics());

    [[nodiscard]] expression::Expr evaluate(const expression::Expr& expression);
    [[nodiscard]] expression::Expr evaluate(
        const expression::Expr& expression,
        const expression::OriginMap& origins);
    [[nodiscard]] expression::Expr evaluate(
        const expression::Expr& expression,
        const expression::OriginMap& origins,
        EvaluationContext context);

    void setDepthLimit(std::size_t limit);
    [[nodiscard]] std::size_t depthLimit() const noexcept;
    void setEvaluationLimits(EvaluationLimits limits);
    [[nodiscard]] const EvaluationLimits& evaluationLimits() const noexcept;
    void reseedRandomFromEntropy();

private:
    struct ActiveFunctionCall final {
        expression::Symbol name;
        std::vector<expression::Expr> arguments;

        [[nodiscard]] bool operator==(const ActiveFunctionCall&) const = default;
    };

    struct ActiveUserFunctionFrame final {
        ActiveFunctionCall call;
        const UserFunctionDefinition* definition = nullptr;
        std::optional<source::SourceReference> callOrigin;
    };

    Environment& environment_;
    const BuiltinRegistry& registry_;
    UserFunctionRegistry* userFunctions_ = nullptr;
    const symbols::SymbolRegistry& symbolRegistry_;
    const mathematics::MathRegistry& mathematics_;
    const mathematics::AngleSemantics& angleSemantics_;
    random::RandomEngine randomEngine_;
    builtins::FourierTransformCache fourierTransformCache_;
    const expression::OriginMap* origins_ = nullptr;
    const EvaluationContext* context_ = nullptr;
    EvaluationLimits limits_{};
    std::vector<expression::Symbol> resolvingSymbols_;
    std::vector<ActiveUserFunctionFrame> activeUserFunctions_;
    // 外側のsymbolic binderを内側frontendの先行materializeから保護する。
    // Environmentそのものを書き換えず，この評価器の再入時だけ自由記号として扱う。
    std::vector<expression::Symbol> materializationProtectedSymbols_;
    // N[...] の評価中だけ有効な要求精度スタック。
    // 外側のexact評価規則は変えず、FFT等のprecision-aware builtinだけがこの情報を参照する。
    std::vector<approximation::ApproximationContext> approximationContexts_;

    [[nodiscard]] expression::Expr evaluateMachine(
        const expression::Expr& expression,
        const expression::OriginMap* origins,
        const EvaluationContext* context);

    [[nodiscard]] expression::Expr dispatchBuiltin(
        const BuiltinDefinition& definition,
        const expression::CallExpr& call,
        std::span<const expression::Expr> arguments);
    [[nodiscard]] expression::Expr evaluateSet(
        std::span<const expression::Expr> arguments);
    [[nodiscard]] expression::Expr evaluateSetDelayed(
        const expression::CallExpr& call,
        std::span<const expression::Expr> arguments);
    [[nodiscard]] expression::Expr evaluateHistory(
        std::span<const expression::Expr> arguments);
    [[nodiscard]] expression::Expr evaluateIndexedHistory(
        std::span<const expression::Expr> arguments,
        bool input);
    [[nodiscard]] expression::Expr resolveHeldHistoryReferences(
        const expression::Expr& expression);
    [[nodiscard]] expression::Expr materializeSafeHeldFrontends(
        const expression::Expr& expression);
    [[nodiscard]] expression::Expr materializeSafeHeldFrontends(
        const expression::Expr& expression,
        std::span<const expression::Symbol> protectedSymbols);
    [[nodiscard]] expression::Expr evaluateDefinitions() const;
    [[nodiscard]] expression::Expr evaluateUndefine(std::span<const expression::Expr> arguments);
    void emitWarning(std::string_view code, std::string message);
    void emitInfo(std::string_view code, std::string message,
        std::optional<expression::Expr> previousExpression = std::nullopt);
    [[nodiscard]] expression::Expr finalizeNumericalApproximation(
        const expression::Expr& value,
        std::size_t precisionDigits,
        bool warnOnFailure = true);
    [[nodiscard]] const approximation::ApproximationContext* currentApproximationContext() const noexcept;
    [[nodiscard]] std::optional<source::SourceReference> originOf(
        const expression::Expr& expression) const;
};

} // namespace mmcal::evaluation
