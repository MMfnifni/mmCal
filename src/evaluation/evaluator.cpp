// 非再帰タスクスタック式評価器
#include "evaluator.hpp"

#include "expression/array_utils.hpp"

#include "approximation/approximation_context.hpp"
#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/expression_interval.hpp"
#include "builtins/names.hpp"
#include "builtins/iteration.hpp"
#include "error/error_message.hpp"
#include "evaluation/iterator_spec.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/integer_algorithms.hpp"
#include "simplification/simplifier.hpp"
#include "simplification/expression_cost.hpp"
#include "solver/solution_set.hpp"
#include "symbolic/series.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

namespace mmcal::evaluation {
namespace {

struct EvaluateTask final {
    expression::Expr expression;
    const expression::OriginMap* origins = nullptr;
    std::size_t depth = 1;
};

struct PushResultTask final {
    expression::Expr expression;
};

struct FinishSymbolTask final {
    expression::Symbol symbol;
};

struct BuildArrayTask final {
    expression::Expr sourceExpression;
    std::vector<std::size_t> expressionIndices;
    const expression::OriginMap* origins = nullptr;
};

struct BuildListTask final {
    expression::Expr sourceExpression;
    std::size_t elementCount = 0;
    const expression::OriginMap* origins = nullptr;
};

struct DispatchBuiltinTask final {
    expression::Expr expression;
    const BuiltinDefinition* definition = nullptr;
    const expression::OriginMap* origins = nullptr;
    std::size_t depth = 1;
};

struct IfConditionTask final {
    expression::Expr expression;
    const expression::OriginMap* origins = nullptr;
    std::size_t childDepth = 1;
};

struct CasesConditionTask final {
    expression::Expr expression;
    std::size_t branchIndex = 0;
    const expression::OriginMap* origins = nullptr;
    std::size_t childDepth = 1;
};

struct EnterUserFunctionTask final {
    expression::Expr expression;
    const UserFunctionDefinition* definition = nullptr;
    const expression::OriginMap* origins = nullptr;
    std::size_t childDepth = 1;
};

struct FinishUserFunctionTask final {};

struct BeginLocalScopeTask final {
    expression::Symbol symbol;
    expression::Expr value;
};

struct EndLocalScopeTask final {};

struct BuildTableTask final {
    std::size_t elementCount = 0;
};

// N[...] は第1引数を先にexact評価しない。precisionだけを確定してから子式を評価し、
// precision-aware builtinへ要求精度を伝播した後、最後に従来どおりcertified decimalへ落とす。
struct BeginNumericalApproximationTask final {
    expression::Expr expression;
    const expression::OriginMap* origins = nullptr;
    std::size_t childDepth = 1;
};

struct FinishNumericalApproximationTask final {
    std::size_t precisionDigits = approximation::ApproximationContext::defaultDecimalDigits;
    std::size_t warningCountBefore = 0;
};

using EvaluationTask = std::variant<
    EvaluateTask,
    PushResultTask,
    FinishSymbolTask,
    BuildArrayTask,
    BuildListTask,
    DispatchBuiltinTask,
    IfConditionTask,
    CasesConditionTask,
    EnterUserFunctionTask,
    FinishUserFunctionTask,
    BeginLocalScopeTask,
    EndLocalScopeTask,
    BuildTableTask,
    BeginNumericalApproximationTask,
    FinishNumericalApproximationTask>;

struct TaskSource final {
    const expression::Expr* expression = nullptr;
    const expression::OriginMap* origins = nullptr;
};

template <class... Ts>
struct Overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
Overloaded(Ts...) -> Overloaded<Ts...>;

[[nodiscard]] TaskSource taskSource(const EvaluationTask& task) noexcept {
    return std::visit(Overloaded{
        [](const EvaluateTask& value) noexcept {
            return TaskSource{&value.expression, value.origins};
        },
        [](const BuildArrayTask& value) noexcept {
            return TaskSource{&value.sourceExpression, value.origins};
        },
        [](const BuildListTask& value) noexcept {
            return TaskSource{&value.sourceExpression, value.origins};
        },
        [](const DispatchBuiltinTask& value) noexcept {
            return TaskSource{&value.expression, value.origins};
        },
        [](const IfConditionTask& value) noexcept {
            return TaskSource{&value.expression, value.origins};
        },
        [](const CasesConditionTask& value) noexcept {
            return TaskSource{&value.expression, value.origins};
        },
        [](const EnterUserFunctionTask& value) noexcept {
            return TaskSource{&value.expression, value.origins};
        },
        [](const BeginNumericalApproximationTask& value) noexcept {
            return TaskSource{&value.expression, value.origins};
        },
        [](const auto&) noexcept {
            return TaskSource{};
        }
    }, task);
}

[[nodiscard]] std::optional<source::SourceReference> findOrigin(
    const TaskSource& source) {
    if (!source.expression || !source.origins)
        return std::nullopt;
    return source.origins->find(*source.expression);
}

// {variable, lower, upper} 形式のiterator specではvariableだけを保持し、lower/upperは通常の評価規則へ通す。
// integrate/sum/product等のbinderで共用する。
void scheduleIteratorSpecArgument(
    std::vector<EvaluationTask>& tasks,
    const expression::Expr& argument,
    const expression::OriginMap* origins,
    std::size_t depth) {
    const expression::ArrayExpr* spec = rangeIteratorArray(argument);
    if (!spec) {
        tasks.emplace_back(PushResultTask{argument});
        return;
    }

    tasks.emplace_back(BuildArrayTask{argument, {1, 2}, origins});
    tasks.emplace_back(EvaluateTask{spec->element(2), origins, depth});
    tasks.emplace_back(EvaluateTask{spec->element(1), origins, depth});
}

// table iteratorではvariableだけを保持し、残りのrange引数だけを通常評価する。
void scheduleTableIteratorSpecArgument(
    std::vector<EvaluationTask>& tasks,
    const expression::Expr& argument,
    const expression::OriginMap* origins,
    std::size_t depth) {
    const expression::ArrayExpr* spec = tableIteratorArray(argument);
    if (!spec) {
        tasks.emplace_back(PushResultTask{argument});
        return;
    }

    std::vector<std::size_t> indices;
    indices.reserve(spec->size() - 1);
    for (std::size_t i = 1; i < spec->size(); ++i)
        indices.push_back(i);
    tasks.emplace_back(BuildArrayTask{argument, indices, origins});
    for (std::size_t i = spec->size(); i-- > 1;)
        tasks.emplace_back(EvaluateTask{spec->element(i), origins, depth});
}

[[nodiscard]] std::optional<std::vector<expression::Expr>> tableIteratorSequence(
    const expression::Expr& expression) {
    std::vector<expression::Expr> specs;

    if (expression.isArray()) {
        const expression::ArrayExpr& array = expression.asArray();
        if (array.rank() != 2 || array.shape[0] == 0 || array.shape[1] < 2 || array.shape[1] > 4)
            return std::nullopt;

        specs.reserve(array.shape[0]);
        for (std::size_t row = 0; row < array.shape[0]; ++row) {
            std::vector<expression::Expr> elements;
            elements.reserve(array.shape[1]);
            for (std::size_t column = 0; column < array.shape[1]; ++column)
                elements.push_back(array.element(row * array.shape[1] + column));
            expression::Expr spec = expression::Expr::array({array.shape[1]}, std::move(elements));
            if (!parseTableIteratorSpec(spec))
                return std::nullopt;
            specs.push_back(std::move(spec));
        }
        return specs;
    }

    if (!expression.isList() || expression.asList().elements.empty())
        return std::nullopt;

    specs.reserve(expression.asList().elements.size());
    for (const expression::Expr& element : expression.asList().elements) {
        if (!parseTableIteratorSpec(element))
            return std::nullopt;
        specs.push_back(element);
    }
    return specs;
}

[[nodiscard]] expression::Expr mapLeafCalls(
    const expression::Expr& container,
    const expression::Symbol& function) {
    if (container.isArray()) {
        const auto& array = container.asArray();
        std::vector<expression::Expr> mapped;
        mapped.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i)
            mapped.push_back(expression::Expr::call(function, {array.element(i)}));
        return expression::Expr::array(array.shape, std::move(mapped));
    }

    if (container.isList()) {
        std::vector<expression::Expr> mapped;
        mapped.reserve(container.asList().elements.size());
        for (const auto& element : container.asList().elements) {
            if (element.isArray() || element.isList())
                mapped.push_back(mapLeafCalls(element, function));
            else
                mapped.push_back(expression::Expr::call(function, {element}));
        }
        return expression::Expr::list(std::move(mapped));
    }

    error::throwCalcError(
        error::CalcErrorType::Type,
        "map expects an Array or brace value as its second argument");
}

[[nodiscard]] std::vector<expression::Expr> takeResults(
    std::vector<expression::Expr>& results,
    std::size_t count) {
    if (results.size() < count)
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "Evaluation result stack is inconsistent");

    const auto begin = results.end() - static_cast<std::ptrdiff_t>(count);
    std::vector<expression::Expr> values{begin, results.end()};
    results.erase(begin, results.end());
    return values;
}

[[nodiscard]] std::string arityMessage(const BuiltinDefinition& definition) {
    if (definition.minimumArguments == definition.maximumArguments)
        return std::string{definition.name()} + " expects "
            + std::to_string(definition.minimumArguments) + " argument(s)";
    if (definition.maximumArguments == BuiltinDefinition::unlimited)
        return std::string{definition.name()} + " expects at least "
            + std::to_string(definition.minimumArguments) + " argument(s)";

    return std::string{definition.name()} + " expects between "
        + std::to_string(definition.minimumArguments) + " and "
        + std::to_string(definition.maximumArguments) + " arguments";
}

[[nodiscard]] std::string userFunctionArityMessage(
    const UserFunctionRegistry& registry,
    const expression::Symbol& name,
    std::size_t actualArity) {
    const std::vector<std::size_t> available = registry.arities(name);
    std::string message = name.name() + " is not defined for "
        + std::to_string(actualArity) + " argument(s)";

    if (available.empty())
        return message;

    message += "; available arities: ";
    for (std::size_t i = 0; i < available.size(); ++i) {
        if (i != 0)
            message += ", ";
        message += std::to_string(available[i]);
    }

    return message;
}

struct HistoryIndex final {
    bool relative = false;
    std::size_t magnitude = 0;
};

[[nodiscard]] std::optional<HistoryIndex> historyIndexValue(
    const expression::Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;

    const numeric::BigInt& value = expression.asNumber().asReal().asInteger();
    if (value.isZero())
        return std::nullopt;

    const auto magnitude = numeric::tryToUint64(value.abs());
    if (!magnitude || *magnitude > std::numeric_limits<std::size_t>::max())
        return std::nullopt;

    return HistoryIndex{value.isNegative(), static_cast<std::size_t>(*magnitude)};
}

[[nodiscard]] std::optional<std::size_t> positiveSizeValue(
    const expression::Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;

    const numeric::BigInt& value = expression.asNumber().asReal().asInteger();
    if (value.isNegative() || value.isZero())
        return std::nullopt;

    const std::string text = value.toString();
    std::size_t result = 0;
    const auto conversion = std::from_chars(text.data(), text.data() + text.size(), result);
    if (conversion.ec != std::errc{} || conversion.ptr != text.data() + text.size())
        return std::nullopt;

    return result;
}

constexpr std::size_t maximumNumericalApproximationRefinements = 16;

[[nodiscard]] std::size_t nextApproximationGuardDigits(std::size_t current) {
    const std::size_t growth = std::max<std::size_t>(8, current / 2);
    if (growth > std::numeric_limits<std::size_t>::max() - current)
        throw std::overflow_error("Approximation precision is too large");
    return current + growth;
}

[[nodiscard]] numeric::DecimalApproximation reduceApproximationPrecision(
    const numeric::DecimalApproximation& value,
    std::size_t precisionDigits) {
    const std::size_t available = value.requestedSignificantDigits();
    if (available != 0 && precisionDigits >= available)
        return value;

    // VerifiedApproximationは数値算法の候補点でありrigorous enclosureではない。
    // 精度を下げてもcertified constructorへ通して証明済みに昇格させず，
    // 元の数値点を低い表示精度へ再量子化したVerified値として保持する。
    if (!value.hasRigorousEnclosure())
        return numeric::DecimalApproximation::fromVerifiedValueSignificant(
            numeric::RealNumber{value.certifiedLower()}, precisionDigits);

    // 桁数を下げる場合も元のInformationEnclosureを保持し，
    // 新しい表示丸め量子だけを追加で情報量上限へ反映する。
    if (const auto rounded = numeric::DecimalApproximation::fromCertifiedIntervalWithInformationSignificant(
        value.certifiedLower(),
        value.certifiedUpper(),
        value.informationLower(),
        value.informationUpper(),
        precisionDigits))
        return *rounded;

    // 元の保証区間が粗く、より低い桁への丸め境界を跨ぐ特殊caseでは、
    // 情報を捨てて誤った桁へ丸めず既存の近似値を保持する。
    return value;
}

[[nodiscard]] expression::Expr reduceApproximationPrecision(
    const numeric::ComplexDecimalApproximation& value,
    std::size_t precisionDigits) {
    const auto real = reduceApproximationPrecision(value.real(), precisionDigits);
    const auto imaginary = reduceApproximationPrecision(value.imaginary(), precisionDigits);
    return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
        real, imaginary, value.realExactlyZero(), value.imaginaryExactlyZero())};
}

[[nodiscard]] bool isBooleanExpression(const expression::Expr& expression) {
    if (expression.isBoolean())
        return true;
    if (!expression.isCall())
        return false;

    const std::string_view head = expression.asCall().head.view();
    return head == builtins::names::less
        || head == builtins::names::lessEqual
        || head == builtins::names::greater
        || head == builtins::names::greaterEqual
        || head == builtins::names::equal
        || head == builtins::names::notEqual
        || head == builtins::names::logicalAnd;
}

void checkExpressionBigIntegerBits(const expression::Expr& expression) {
    if (!currentEvaluationBudget())
        return;

    const auto realBits = [](const numeric::RealNumber& value) {
        if (value.isInteger())
            return value.asInteger().bitLength();
        const numeric::Rational& rational = value.asRational();
        return std::max(
            rational.numerator().bitLength(),
            rational.denominator().bitLength());
    };
    const auto checkNumber = [&](const numeric::Number& value) {
        const std::size_t bits = value.isReal()
            ? realBits(value.asReal())
            : std::max(
                realBits(value.asComplex().real),
                realBits(value.asComplex().imaginary));
        checkEvaluationBigIntegerBits(bits);
    };

    std::vector<expression::Expr> pending{expression};
    while (!pending.empty()) {
        expression::Expr current = std::move(pending.back());
        pending.pop_back();
        if (current.isNumber()) {
            checkNumber(current.asNumber());
            continue;
        }
        if (current.isCall()) {
            for (const expression::Expr& argument : current.asCall().arguments)
                pending.push_back(argument);
            continue;
        }
        if (current.isArray()) {
            const expression::ArrayExpr& array = current.asArray();
            for (std::size_t i = 0; i < array.size(); ++i) {
                switch (array.storedKindAt(i)) {
                case expression::ArrayStorageKind::Integer:
                case expression::ArrayStorageKind::Rational:
                case expression::ArrayStorageKind::Number:
                    checkNumber(array.exactNumber(i));
                    break;
                case expression::ArrayStorageKind::Generic:
                    pending.push_back(array.element(i));
                    break;
                case expression::ArrayStorageKind::DecimalApproximation:
                case expression::ArrayStorageKind::ComplexDecimalApproximation:
                    break;
                }
            }
            continue;
        }
        if (current.isList()) {
            for (const expression::Expr& element : current.asList().elements)
                pending.push_back(element);
            continue;
        }
        if (current.isSolutionSet()) {
            const solver::SolutionSet& solutions = current.asSolutionSet();
            for (const solver::SolutionBranch& branch : solutions.branches())
                for (const solver::SolutionBinding& binding : branch.bindings)
                    pending.push_back(binding.value);
            for (const solver::SolutionCase& solutionCase : solutions.cases())
                for (const solver::SolutionBranch& branch : solutionCase.branches)
                    for (const solver::SolutionBinding& binding : branch.bindings)
                        pending.push_back(binding.value);
        }
    }
}

void consumeGeneratedExpression(const expression::Expr& expression) {
    if (!currentEvaluationBudget())
        return;

    checkExpressionBigIntegerBits(expression);
    consumeEvaluationBudget(
        EvaluationResource::GeneratedNode,
        simplification::measureExpressionCost(expression).nodes);
}

} // namespace

Evaluator::Evaluator(
    Environment& environment,
    const BuiltinRegistry& registry,
    UserFunctionRegistry* userFunctions,
    const symbols::SymbolRegistry& symbolRegistry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angleSemantics)
    : environment_(environment),
      registry_(registry),
      userFunctions_(userFunctions),
      symbolRegistry_(symbolRegistry),
      mathematics_(mathematics),
      angleSemantics_(angleSemantics) {}

expression::Expr Evaluator::evaluate(const expression::Expr& expression) {
    EvaluationBudget budget{limits_};
    EvaluationContext context{};
    context.budget = &budget;
    return evaluateMachine(expression, nullptr, &context);
}

expression::Expr Evaluator::evaluate(
    const expression::Expr& expression,
    const expression::OriginMap& origins) {
    EvaluationBudget budget{limits_};
    EvaluationContext context{};
    context.budget = &budget;
    return evaluateMachine(expression, &origins, &context);
}

expression::Expr Evaluator::evaluate(
    const expression::Expr& expression,
    const expression::OriginMap& origins,
    EvaluationContext context) {
    if (!context.budget) {
        EvaluationBudget budget{limits_};
        context.budget = &budget;
        return evaluateMachine(expression, &origins, &context);
    }
    return evaluateMachine(expression, &origins, &context);
}

void Evaluator::setDepthLimit(std::size_t limit) {
    if (limit == 0)
        throw std::invalid_argument("Evaluation depth limit must be greater than zero");

    limits_.maxDepth = limit;
}

std::size_t Evaluator::depthLimit() const noexcept {
    return limits_.maxDepth;
}

void Evaluator::setEvaluationLimits(EvaluationLimits limits) {
    if (limits.maxDepth == 0)
        throw std::invalid_argument("Evaluation depth limit must be greater than zero");
    limits_ = std::move(limits);
}

const EvaluationLimits& Evaluator::evaluationLimits() const noexcept {
    return limits_;
}

void Evaluator::reseedRandomFromEntropy() {
    static_cast<void>(randomEngine_.reseedFromEntropy());
}

expression::Expr Evaluator::evaluateMachine(
    const expression::Expr& expression,
    const expression::OriginMap* origins,
    const EvaluationContext* context) {
    const std::size_t initialLocalDepth = environment_.localDepth();
    EvaluationBudgetScope budgetScope{context ? context->budget : nullptr};
    context_ = context;
    origins_ = origins;
    resolvingSymbols_.clear();
    activeUserFunctions_.clear();
    approximationContexts_.clear();

    std::vector<EvaluationTask> tasks;
    std::vector<expression::Expr> results;
    tasks.reserve(64);
    results.reserve(64);
    tasks.emplace_back(EvaluateTask{expression, origins, 1});

    auto cleanup = [&] {
        while (environment_.localDepth() > initialLocalDepth)
            environment_.popScope();
        resolvingSymbols_.clear();
        activeUserFunctions_.clear();
        approximationContexts_.clear();
        origins_ = nullptr;
        context_ = nullptr;
    };

    try {
        while (!tasks.empty()) {
            consumeEvaluationBudget(EvaluationResource::EvaluationStep);
            EvaluationTask task = std::move(tasks.back());
            tasks.pop_back();
            const TaskSource source = taskSource(task);
            origins_ = source.origins;

            try {
                std::visit(Overloaded{
                    [&](const EvaluateTask& current) {
                        if (EvaluationBudget* budget = currentEvaluationBudget())
                            budget->checkDepth(current.depth);

                        switch (current.expression.kind()) {
                        case expression::ExprKind::Number:
                            checkExpressionBigIntegerBits(current.expression);
                            results.push_back(current.expression);
                            return;

                        case expression::ExprKind::DecimalApproximation:
                        case expression::ExprKind::ComplexDecimalApproximation:
                        case expression::ExprKind::Boolean:
                        case expression::ExprKind::String:
                            results.push_back(current.expression);
                            return;

                        case expression::ExprKind::SolutionSet:
                            checkExpressionBigIntegerBits(current.expression);
                            results.push_back(current.expression);
                            return;

                        case expression::ExprKind::Symbol: {
                            const expression::Symbol symbol = current.expression.asSymbol();
                            if (std::find(materializationProtectedSymbols_.begin(),
                                    materializationProtectedSymbols_.end(), symbol)
                                != materializationProtectedSymbols_.end()) {
                                results.push_back(current.expression);
                                return;
                            }
                            const expression::Expr* binding = environment_.find(symbol);
                            if (!binding) {
                                // 未束縛SymbolはCASの自由記号として保持する。
                                // simplify / D / integrate / solve等のsymbolic evaluationはこの意味論を共有する。
                                results.push_back(current.expression);
                                return;
                            }

                            if (std::find(resolvingSymbols_.begin(), resolvingSymbols_.end(), symbol)
                                != resolvingSymbols_.end())
                                error::throwCalcError(
                                    error::CalcErrorType::Evaluation,
                                    "Cyclic symbol definition: " + symbol.name());

                            resolvingSymbols_.push_back(symbol);
                            tasks.emplace_back(FinishSymbolTask{symbol});
                            tasks.emplace_back(EvaluateTask{*binding, current.origins, current.depth + 1});
                            return;
                        }

                        case expression::ExprKind::Array: {
                            const expression::ArrayExpr& array = current.expression.asArray();
                            consumeEvaluationBudget(
                                EvaluationResource::DenseArrayElement, array.size());
                            checkExpressionBigIntegerBits(current.expression);
                            const auto entries = array.expressionEntries();
                            if (entries.empty()) {
                                results.push_back(current.expression);
                                return;
                            }

                            std::vector<std::size_t> indices;
                            indices.reserve(entries.size());
                            for (const auto& entry : entries)
                                indices.push_back(entry.index);
                            tasks.emplace_back(BuildArrayTask{
                                current.expression, std::move(indices), current.origins});
                            for (auto iterator = entries.rbegin(); iterator != entries.rend(); ++iterator)
                                tasks.emplace_back(EvaluateTask{
                                    iterator->expression, current.origins, current.depth + 1});
                            return;
                        }

                        case expression::ExprKind::List: {
                            const expression::ListExpr& list = current.expression.asList();
                            tasks.emplace_back(BuildListTask{
                                current.expression,
                                list.elements.size(),
                                current.origins
                            });
                            for (auto iterator = list.elements.rbegin(); iterator != list.elements.rend(); ++iterator)
                                tasks.emplace_back(EvaluateTask{*iterator, current.origins, current.depth + 1});
                            return;
                        }

                        case expression::ExprKind::Call: {
                            const expression::CallExpr& call = current.expression.asCall();
                            if (const BuiltinDefinition* definition = registry_.find(call.head)) {
                                if (!definition->acceptsArity(call.arguments.size()))
                                    error::throwCalcError(
                                        error::CalcErrorType::Type,
                                        arityMessage(*definition));

                                // Ifだけは条件結果を見てから片側の枝だけを積む。
                                if (definition->id == BuiltinId::If) {
                                    tasks.emplace_back(IfConditionTask{
                                        current.expression,
                                        current.origins,
                                        current.depth + 1
                                    });
                                    tasks.emplace_back(EvaluateTask{
                                        call.arguments[0],
                                        current.origins,
                                        current.depth + 1
                                    });
                                    return;
                                }

                                if (definition->id == BuiltinId::Cases) {
                                    const auto scheduleBranch = [&](std::size_t index) {
                                        if (index >= call.arguments.size()) {
                                            const auto* indeterminate = symbolRegistry_.find(
                                                symbols::PredefinedSymbolId::Indeterminate);
                                            if (!indeterminate)
                                                error::throwCalcError(
                                                    error::CalcErrorType::Internal,
                                                    "Indeterminate symbol is not registered");
                                            results.emplace_back(indeterminate->symbol);
                                            return;
                                        }
                                        const expression::Expr& branchExpression = call.arguments[index];
                                        if (!branchExpression.isCall()
                                            || !branchExpression.asCall().head.sameIdentity(
                                                registry_.symbol(BuiltinId::CaseBranch))
                                            || branchExpression.asCall().arguments.empty()
                                            || branchExpression.asCall().arguments.size() > 2)
                                            error::throwCalcError(
                                                error::CalcErrorType::Type,
                                                "cases contains an invalid branch");
                                        const auto& branch = branchExpression.asCall().arguments;
                                        if (branch.size() == 1) {
                                            tasks.emplace_back(EvaluateTask{
                                                branch[0], current.origins, current.depth + 1});
                                            return;
                                        }
                                        tasks.emplace_back(CasesConditionTask{
                                            current.expression, index, current.origins, current.depth + 1});
                                        tasks.emplace_back(EvaluateTask{
                                            branch[1], current.origins, current.depth + 1});
                                    };
                                    scheduleBranch(0);
                                    return;
                                }

                                // Nだけは値を先にexact評価すると、FFT等が巨大なexact中間式を構築した後でしか
                                // 近似要求を知れない。第1引数を保持し、precisionを先に評価してから子式へ伝播する。
                                if (definition->id == BuiltinId::NumericalApproximation) {
                                    tasks.emplace_back(BeginNumericalApproximationTask{
                                        current.expression,
                                        current.origins,
                                        current.depth + 1
                                    });
                                    if (call.arguments.size() == 2)
                                        tasks.emplace_back(EvaluateTask{
                                            call.arguments[1],
                                            current.origins,
                                            current.depth + 1
                                        });
                                    return;
                                }

                                tasks.emplace_back(DispatchBuiltinTask{
                                    current.expression,
                                    definition,
                                    current.origins,
                                    current.depth
                                });

                                for (std::size_t index = call.arguments.size(); index-- > 0;) {
                                    if (definition->argumentEvaluation == ArgumentEvaluation::HoldFirstAndIteratorSpec
                                        && index == 1) {
                                        scheduleIteratorSpecArgument(
                                            tasks,
                                            call.arguments[index],
                                            current.origins,
                                            current.depth + 1);
                                        continue;
                                    }
                                    if (definition->argumentEvaluation == ArgumentEvaluation::HoldFirstAndTableIteratorSpec
                                        && index == 1) {
                                        scheduleTableIteratorSpecArgument(
                                            tasks,
                                            call.arguments[index],
                                            current.origins,
                                            current.depth + 1);
                                        continue;
                                    }

                                    const bool held = definition->argumentEvaluation == ArgumentEvaluation::HoldAll
                                        || (definition->argumentEvaluation == ArgumentEvaluation::HoldFirst && index == 0)
                                        || (definition->argumentEvaluation == ArgumentEvaluation::HoldFirstTwo && index < 2)
                                        || (definition->argumentEvaluation == ArgumentEvaluation::HoldFirstAndIteratorSpec && index == 0)
                                        || (definition->argumentEvaluation == ArgumentEvaluation::HoldFirstAndTableIteratorSpec && index == 0);
                                    if (held)
                                        tasks.emplace_back(PushResultTask{call.arguments[index]});
                                    else
                                        tasks.emplace_back(EvaluateTask{
                                            call.arguments[index],
                                            current.origins,
                                            current.depth + 1
                                        });
                                }
                                return;
                            }

                            if (userFunctions_) {
                                if (const UserFunctionDefinition* definition = userFunctions_->find(
                                    call.head,
                                    call.arguments.size())) {
                                    tasks.emplace_back(EnterUserFunctionTask{
                                        current.expression,
                                        definition,
                                        current.origins,
                                        current.depth + 1
                                    });
                                    for (auto iterator = call.arguments.rbegin(); iterator != call.arguments.rend(); ++iterator)
                                        tasks.emplace_back(EvaluateTask{
                                            *iterator,
                                            current.origins,
                                            current.depth + 1
                                        });
                                    return;
                                }

                                if (userFunctions_->contains(call.head))
                                    error::throwCalcError(
                                        error::CalcErrorType::Type,
                                        userFunctionArityMessage(*userFunctions_, call.head, call.arguments.size()));
                            }

                            error::throwCalcError(
                                error::CalcErrorType::Name,
                                "Unknown function: " + call.head.name());
                        }
                        }

                        error::throwCalcError(
                            error::CalcErrorType::Internal,
                            "Unknown expression kind");
                    },
                    [&](const PushResultTask& current) {
                        results.push_back(current.expression);
                    },
                    [&](const FinishSymbolTask& current) {
                        if (resolvingSymbols_.empty() || resolvingSymbols_.back() != current.symbol)
                            error::throwCalcError(
                                error::CalcErrorType::Internal,
                                "Symbol resolution stack is inconsistent");
                        resolvingSymbols_.pop_back();
                    },
                    [&](const BuildArrayTask& current) {
                        std::vector<expression::Expr> elements = takeResults(
                            results, current.expressionIndices.size());
                        expression::Expr rebuilt = expression::rebuildEvaluatedArray(
                            current.sourceExpression.asArray(),
                            current.expressionIndices,
                            std::move(elements));
                        consumeGeneratedExpression(rebuilt);
                        results.push_back(std::move(rebuilt));
                    },
                    [&](const BuildListTask& current) {
                        std::vector<expression::Expr> elements = takeResults(results, current.elementCount);
                        expression::Expr rebuilt = expression::braceValue(std::move(elements));
                        consumeGeneratedExpression(rebuilt);
                        results.push_back(std::move(rebuilt));
                    },
                    [&](const DispatchBuiltinTask& current) {
                        const expression::CallExpr& call = current.expression.asCall();
                        std::vector<expression::Expr> arguments = takeResults(results, call.arguments.size());

                        if (current.definition->id == BuiltinId::Map) {
                            if (!arguments[0].isSymbol())
                                error::throwCalcError(
                                    error::CalcErrorType::Type,
                                    "map first argument must be a function symbol");
                            expression::Expr mapped = mapLeafCalls(
                                arguments[1], arguments[0].asSymbol());
                            tasks.emplace_back(EvaluateTask{
                                std::move(mapped), current.origins, current.depth + 1});
                            return;
                        }

                        if (current.definition->id == BuiltinId::Table) {
                            // {{i,...},{j,...},...} は左から外側iteratorとして解釈する。
                            // 既存の単一iterator評価器へ nested table callとして落とすことで，
                            // scope・range・budget semanticsを二重実装しない。
                            if (const auto sequence = tableIteratorSequence(arguments[1])) {
                                expression::Expr nested = arguments[0];
                                for (auto iterator = sequence->rbegin(); iterator != sequence->rend(); ++iterator)
                                    nested = expression::Expr::call(call.head, {std::move(nested), *iterator});
                                tasks.emplace_back(EvaluateTask{
                                    std::move(nested), current.origins, current.depth + 1});
                                return;
                            }

                            const auto spec = parseTableIteratorSpec(arguments[1]);
                            if (!spec)
                                error::throwCalcError(
                                    error::CalcErrorType::Type,
                                    "table iterator must be {symbol, end}, {symbol, lower, upper}, {symbol, lower, upper, step}, or a brace value of such iterators");
                            std::vector<expression::Expr> values = builtins::exactRangeValues(
                                spec->rangeArguments, "table");
                            tasks.emplace_back(BuildTableTask{values.size()});
                            for (std::size_t i = values.size(); i-- > 0;) {
                                tasks.emplace_back(EndLocalScopeTask{});
                                tasks.emplace_back(EvaluateTask{
                                    arguments[0], current.origins, current.depth + 1});
                                tasks.emplace_back(BeginLocalScopeTask{
                                    spec->variable, std::move(values[i])});
                            }
                            return;
                        }

                        expression::Expr dispatched = [&]() {
                            const bool hasDedicatedApproxArithmetic =
                                current.definition->id == BuiltinId::Add
                                || current.definition->id == BuiltinId::Subtract
                                || current.definition->id == BuiltinId::Multiply
                                || current.definition->id == BuiltinId::Divide
                                || current.definition->id == BuiltinId::Negate;
                            if (!hasDedicatedApproxArithmetic) {
                                expression::Expr approximateCall = expression::Expr::call(
                                    call.head,
                                    std::vector<expression::Expr>{arguments.begin(), arguments.end()});
                                if (const auto approximate = approximation::evaluateApproximateExpression(
                                        approximateCall, registry_, mathematics_, angleSemantics_))
                                    return *approximate;
                            }

                            expression::Expr result = dispatchBuiltin(*current.definition, call, arguments);
                            // log2/log10/fract等は通常Evaluatorでprimitiveへrewriteされる。
                            // rewrite後の式にも近似値が残る場合、同じcertified経路へ一度だけ再投入し、
                            // wrapper builtinごとのDecimalApproximation分岐を増やさない。
                            if (!hasDedicatedApproxArithmetic && result.isCall()) {
                                if (const auto approximate = approximation::evaluateApproximateExpression(
                                        result, registry_, mathematics_, angleSemantics_))
                                    return *approximate;
                            }
                            return result;
                        }();

                        // 組み込み函数が局所的に評価した後、共通Simplifierへ一度通す。
                        // これにより Exp[Log[Pi]] のような複数函数にまたがる安全な書換えを各builtinへ重複実装せず、同じ数学知識へ集約する。
                        const bool pureForPostSimplification =
                            current.definition->id != BuiltinId::Set
                            && current.definition->id != BuiltinId::SetDelayed
                            && current.definition->id != BuiltinId::If
                            && current.definition->id != BuiltinId::Cases
                            && current.definition->id != BuiltinId::CaseBranch
                            && current.definition->id != BuiltinId::History
                            && current.definition->id != BuiltinId::InputHistory
                            && current.definition->id != BuiltinId::OutputHistory
                            && current.definition->id != BuiltinId::Explain
                            && current.definition->id != BuiltinId::Exit
                            && current.definition->id != BuiltinId::Clear
                            && current.definition->id != BuiltinId::Definitions
                            && current.definition->id != BuiltinId::Undefine
                            // 明示的な正規形変換は各engineが返した形そのものに意味がある。
                            // automatic Simplifierで直後に並べ替えない。
                            && current.definition->id != BuiltinId::Expand
                            && current.definition->id != BuiltinId::Factor
                            && current.definition->id != BuiltinId::Collect
                            && current.definition->id != BuiltinId::Solve;
                        if (pureForPostSimplification) {
                            simplification::SimplificationContext simplificationContext{
                                registry_, mathematics_, angleSemantics_};
                            simplificationContext.predefinedSymbols = &symbolRegistry_;
                            dispatched = simplification::Simplifier{}.simplify(
                                dispatched, simplificationContext);
                        }

                        // In[n]は保存した入力Exprを「貼り戻す」意味とし、取得したExprを現在の環境で通常評価する。Out[n]は保存済み結果なので再評価しない。
                        if (current.definition->id == BuiltinId::InputHistory) {
                            consumeGeneratedExpression(dispatched);
                            tasks.emplace_back(EvaluateTask{
                                std::move(dispatched), nullptr, current.depth + 1});
                            return;
                        }
                        consumeGeneratedExpression(dispatched);
                        results.push_back(std::move(dispatched));
                    },
                    [&](const IfConditionTask& current) {
                        std::vector<expression::Expr> conditionResult = takeResults(results, 1);
                        expression::Expr condition = conditionResult.front();
                        const expression::CallExpr& call = current.expression.asCall();

                        if (condition.isBoolean()) {
                            tasks.emplace_back(EvaluateTask{
                                call.arguments[condition.asBoolean() ? 1 : 2],
                                current.origins,
                                current.childDepth
                            });
                            return;
                        }

                        // 真偽未確定の記号命題だけはIf式のまま保持する。
                        if (isBooleanExpression(condition)) {
                            results.push_back(expression::Expr::call(
                                call.head,
                                {std::move(condition), call.arguments[1], call.arguments[2]}));
                            return;
                        }

                        if (current.origins) {
                            if (const auto origin = current.origins->find(call.arguments[0]))
                                error::throwCalcError(
                                    error::CalcErrorType::Type,
                                    "If condition must evaluate to True or False",
                                    *origin);
                        }

                        error::throwCalcError(
                            error::CalcErrorType::Type,
                            "If condition must evaluate to True or False");
                    },
                    [&](const CasesConditionTask& current) {
                        std::vector<expression::Expr> conditionResult = takeResults(results, 1);
                        expression::Expr condition = conditionResult.front();
                        const expression::CallExpr& call = current.expression.asCall();
                        if (current.branchIndex >= call.arguments.size())
                            error::throwCalcError(
                                error::CalcErrorType::Internal,
                                "cases branch index is out of range");
                        const expression::Expr& branchExpression = call.arguments[current.branchIndex];
                        if (!branchExpression.isCall()
                            || !branchExpression.asCall().head.sameIdentity(
                                registry_.symbol(BuiltinId::CaseBranch))
                            || branchExpression.asCall().arguments.size() != 2)
                            error::throwCalcError(
                                error::CalcErrorType::Internal,
                                "cases condition task references an invalid branch");
                        const auto& branch = branchExpression.asCall().arguments;

                        if (condition.isBoolean()) {
                            if (condition.asBoolean()) {
                                tasks.emplace_back(EvaluateTask{
                                    branch[0], current.origins, current.childDepth});
                                return;
                            }

                            const std::size_t next = current.branchIndex + 1;
                            if (next >= call.arguments.size()) {
                                const auto* indeterminate = symbolRegistry_.find(
                                    symbols::PredefinedSymbolId::Indeterminate);
                                if (!indeterminate)
                                    error::throwCalcError(
                                        error::CalcErrorType::Internal,
                                        "Indeterminate symbol is not registered");
                                results.emplace_back(indeterminate->symbol);
                                return;
                            }
                            const expression::Expr& nextExpression = call.arguments[next];
                            if (!nextExpression.isCall()
                                || !nextExpression.asCall().head.sameIdentity(
                                    registry_.symbol(BuiltinId::CaseBranch))
                                || nextExpression.asCall().arguments.empty()
                                || nextExpression.asCall().arguments.size() > 2)
                                error::throwCalcError(
                                    error::CalcErrorType::Type,
                                    "cases contains an invalid branch");
                            const auto& nextBranch = nextExpression.asCall().arguments;
                            if (nextBranch.size() == 1) {
                                tasks.emplace_back(EvaluateTask{
                                    nextBranch[0], current.origins, current.childDepth});
                                return;
                            }
                            tasks.emplace_back(CasesConditionTask{
                                current.expression, next, current.origins, current.childDepth});
                            tasks.emplace_back(EvaluateTask{
                                nextBranch[1], current.origins, current.childDepth});
                            return;
                        }

                        if (isBooleanExpression(condition)) {
                            std::vector<expression::Expr> remaining;
                            remaining.reserve(call.arguments.size() - current.branchIndex);
                            remaining.push_back(expression::Expr::call(
                                registry_.symbol(BuiltinId::CaseBranch),
                                {branch[0], std::move(condition)}));
                            for (std::size_t i = current.branchIndex + 1;
                                 i < call.arguments.size(); ++i)
                                remaining.push_back(call.arguments[i]);
                            results.push_back(expression::Expr::call(
                                call.head, std::move(remaining)));
                            return;
                        }

                        error::throwCalcError(
                            error::CalcErrorType::Type,
                            "cases condition must evaluate to True or False");
                    },
                    [&](const EnterUserFunctionTask& current) {
                        const expression::CallExpr& call = current.expression.asCall();
                        std::vector<expression::Expr> arguments = takeResults(results, call.arguments.size());
                        const ActiveFunctionCall activeCall{current.definition->name, arguments};

                        const auto duplicate = std::find_if(
                            activeUserFunctions_.begin(),
                            activeUserFunctions_.end(),
                            [&](const ActiveUserFunctionFrame& frame) {
                                return frame.call == activeCall;
                            });
                        if (duplicate != activeUserFunctions_.end())
                            error::throwCalcError(
                                error::CalcErrorType::Evaluation,
                                "Cyclic function call: " + current.definition->name.name());

                        environment_.pushScope();
                        for (std::size_t i = 0; i < current.definition->parameters.size(); ++i)
                            environment_.setLocal(current.definition->parameters[i], arguments[i]);

                        std::optional<source::SourceReference> callOrigin;
                        if (current.origins)
                            callOrigin = current.origins->find(current.expression);

                        activeUserFunctions_.push_back(ActiveUserFunctionFrame{
                            activeCall,
                            current.definition,
                            std::move(callOrigin)
                        });

                        tasks.emplace_back(FinishUserFunctionTask{});
                        tasks.emplace_back(EvaluateTask{
                            current.definition->body,
                            &current.definition->origins,
                            current.childDepth
                        });
                    },
                    [&](const FinishUserFunctionTask&) {
                        if (activeUserFunctions_.empty() || environment_.localDepth() <= initialLocalDepth)
                            error::throwCalcError(
                                error::CalcErrorType::Internal,
                                "User-function evaluation stack is inconsistent");

                        environment_.popScope();
                        activeUserFunctions_.pop_back();
                    },
                    [&](const BeginLocalScopeTask& current) {
                        environment_.pushScope();
                        environment_.setLocal(current.symbol, current.value);
                    },
                    [&](const EndLocalScopeTask&) {
                        if (environment_.localDepth() <= initialLocalDepth)
                            error::throwCalcError(
                                error::CalcErrorType::Internal,
                                "Local evaluation scope stack is inconsistent");
                        environment_.popScope();
                    },
                    [&](const BuildTableTask& current) {
                        expression::Expr rebuilt = expression::braceValue(
                            takeResults(results, current.elementCount));
                        consumeGeneratedExpression(rebuilt);
                        results.push_back(std::move(rebuilt));
                    },
                    [&](const BeginNumericalApproximationTask& current) {
                        constexpr std::size_t defaultPrecisionDigits = 16;
                        const expression::CallExpr& call = current.expression.asCall();
                        std::size_t precisionDigits = defaultPrecisionDigits;

                        if (call.arguments.size() == 2) {
                            std::vector<expression::Expr> precisionResult = takeResults(results, 1);
                            const auto requestedDigits = positiveSizeValue(precisionResult.front());
                            if (!requestedDigits) {
                                if (const auto origin = current.origins
                                    ? current.origins->find(call.arguments[1])
                                    : std::nullopt)
                                    error::throwCalcError(
                                        error::CalcErrorType::Type,
                                        "N precision must be a positive integer",
                                        *origin);
                                error::throwCalcError(
                                    error::CalcErrorType::Type,
                                    "N precision must be a positive integer");
                            }
                            precisionDigits = *requestedDigits;
                        }
                        if (EvaluationBudget* budget = currentEvaluationBudget())
                            budget->checkRequestedPrecisionDigits(precisionDigits);

                        std::size_t warningCountBefore = 0;
                        if (context_ && context_->diagnostics)
                            warningCountBefore = static_cast<std::size_t>(std::count_if(
                                context_->diagnostics->begin(), context_->diagnostics->end(),
                                [](const EvaluationDiagnostic& diagnostic) {
                                    return diagnostic.severity == DiagnosticSeverity::Warning;
                                }));

                        approximationContexts_.emplace_back(precisionDigits);
                        tasks.emplace_back(FinishNumericalApproximationTask{
                            precisionDigits, warningCountBefore});
                        tasks.emplace_back(EvaluateTask{
                            call.arguments.front(),
                            current.origins,
                            current.childDepth
                        });
                    },
                    [&](const FinishNumericalApproximationTask& current) {
                        if (approximationContexts_.empty())
                            error::throwCalcError(
                                error::CalcErrorType::Internal,
                                "Numerical approximation context stack is inconsistent");

                        std::size_t warningCountAfter = current.warningCountBefore;
                        if (context_ && context_->diagnostics)
                            warningCountAfter = static_cast<std::size_t>(std::count_if(
                                context_->diagnostics->begin(), context_->diagnostics->end(),
                                [](const EvaluationDiagnostic& diagnostic) {
                                    return diagnostic.severity == DiagnosticSeverity::Warning;
                                }));

                        std::vector<expression::Expr> valueResult = takeResults(results, 1);
                        expression::Expr approximated = finalizeNumericalApproximation(
                            valueResult.front(), current.precisionDigits,
                            warningCountAfter == current.warningCountBefore);
                        approximationContexts_.pop_back();
                        consumeGeneratedExpression(approximated);
                        results.push_back(std::move(approximated));
                    }
                }, task);
            }
            catch (error::CalcError& exception) {
                if (const auto origin = findOrigin(source))
                    exception.attachSourceIfMissing(*origin);
                throw;
            }
            catch (const std::domain_error& exception) {
                if (const auto origin = findOrigin(source))
                    error::throwCalcError(error::CalcErrorType::Domain, exception.what(), *origin);
                error::throwCalcError(error::CalcErrorType::Domain, exception.what());
            }
            catch (const std::overflow_error& exception) {
                if (const auto origin = findOrigin(source))
                    error::throwCalcError(error::CalcErrorType::Overflow, exception.what(), *origin);
                error::throwCalcError(error::CalcErrorType::Overflow, exception.what());
            }
            catch (const std::length_error& exception) {
                if (const auto origin = findOrigin(source))
                    error::throwCalcError(error::CalcErrorType::Overflow, exception.what(), *origin);
                error::throwCalcError(error::CalcErrorType::Overflow, exception.what());
            }
            catch (const std::invalid_argument& exception) {
                if (const auto origin = findOrigin(source))
                    error::throwCalcError(error::CalcErrorType::Type, exception.what(), *origin);
                error::throwCalcError(error::CalcErrorType::Type, exception.what());
            }
        }

        if (results.size() != 1)
            error::throwCalcError(
                error::CalcErrorType::Internal,
                "Evaluation result stack did not finish with one value");
        if (environment_.localDepth() != initialLocalDepth)
            error::throwCalcError(
                error::CalcErrorType::Internal,
                "Evaluation left local scopes active");

        expression::Expr result = results.back();
        cleanup();
        return result;
    }
    catch (error::CalcError& exception) {
        // 函数本体由来のエラーへ、内側から外側の呼出元を付加する。
        if (!activeUserFunctions_.empty()) {
            if (exception.document())
                exception.setSourceLabel("Defined at");

            for (auto iterator = activeUserFunctions_.rbegin();
                iterator != activeUserFunctions_.rend(); ++iterator) {
                if (iterator->callOrigin)
                    exception.addTrace("Called from", *iterator->callOrigin);
            }
        }

        cleanup();
        throw;
    }
    catch (...) {
        cleanup();
        throw;
    }
}

expression::Expr Evaluator::evaluateSet(
    std::span<const expression::Expr> arguments) {
    if (!arguments.front().isSymbol())
        error::throwCalcError(
            error::CalcErrorType::Type,
            "Set requires a symbol as its first argument");

    const expression::Symbol symbol = arguments.front().asSymbol();
    if (registry_.contains(symbol))
        error::throwCalcError(
            error::CalcErrorType::Name,
            "Cannot assign to builtin function name: " + symbol.name());
    if (symbolRegistry_.isProtected(symbol))
        error::throwCalcError(
            error::CalcErrorType::Name,
            "Cannot assign to protected symbol: " + symbol.name());

    const expression::Expr value = arguments.back();
    std::optional<expression::Expr> previous;
    if (!environment_.containsLocal(symbol))
        if (const expression::Expr* existing = environment_.find(symbol))
            previous = *existing;

    environment_.assign(symbol, value);
    if (context_ && context_->definitionsChanged)
        *context_->definitionsChanged = true;
    if (previous && *previous != value)
        emitInfo("definition::redefined", symbol.name() + " redefined", *previous);
    return value;
}

expression::Expr Evaluator::evaluateSetDelayed(
    const expression::CallExpr& call,
    std::span<const expression::Expr> arguments) {
    if (!userFunctions_)
        error::throwCalcError(
            error::CalcErrorType::Evaluation,
            "User-function definitions are not available in this evaluator");

    const expression::Expr& signatureExpression = arguments.front();
    if (!signatureExpression.isCall()
        || signatureExpression.asCall().head.view() != builtins::names::functionSignature)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "SetDelayed requires a function signature as its first argument");

    const expression::CallExpr& signature = signatureExpression.asCall();
    if (signature.arguments.empty() || !signature.arguments.front().isSymbol())
        error::throwCalcError(
            error::CalcErrorType::Type,
            "Function signature is invalid");

    const expression::Symbol name = signature.arguments.front().asSymbol();
    if (registry_.contains(name))
        error::throwCalcError(
            error::CalcErrorType::Name,
            "Cannot redefine builtin function: " + name.name());
    if (symbolRegistry_.isProtected(name))
        error::throwCalcError(
            error::CalcErrorType::Name,
            "Cannot define a function with a protected symbol name: " + name.name());

    std::vector<expression::Symbol> parameters;
    parameters.reserve(signature.arguments.size() - 1);
    for (std::size_t i = 1; i < signature.arguments.size(); ++i) {
        if (!signature.arguments[i].isSymbol())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Function parameters must be symbols");
        parameters.push_back(signature.arguments[i].asSymbol());
    }

    expression::OriginMap definitionOrigins;
    if (origins_)
        definitionOrigins = *origins_;

    std::optional<expression::Expr> previousDefinition;
    if (const UserFunctionDefinition* previous = userFunctions_->find(name, parameters.size())) {
        std::vector<expression::Expr> signatureArguments;
        signatureArguments.reserve(previous->parameters.size() + 1);
        signatureArguments.emplace_back(previous->name);
        for (const expression::Symbol& parameter : previous->parameters)
            signatureArguments.emplace_back(parameter);
        previousDefinition = expression::Expr::call(
            expression::Symbol{builtins::names::setDelayed},
            {expression::Expr::call(
                expression::Symbol{builtins::names::functionSignature},
                std::move(signatureArguments)), previous->body});
    }

    userFunctions_->define(UserFunctionDefinition{
        name,
        std::move(parameters),
        arguments.back(),
        std::move(definitionOrigins)
    });

    if (context_ && context_->definitionsChanged)
        *context_->definitionsChanged = true;
    if (previousDefinition)
        emitInfo("definition::redefined", name.name() + " redefined", *previousDefinition);

    // 定義内容を確認できるよう、現段階ではSetDelayed式自身を返す。
    return expression::Expr::call(call.head, call.arguments);
}

expression::Expr Evaluator::evaluateHistory(
    std::span<const expression::Expr> arguments) {
    const auto depth = positiveSizeValue(arguments.front());
    if (!depth)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "History depth must be a positive integer");
    if (!context_ || *depth > context_->history.size())
        error::throwCalcError(
            error::CalcErrorType::Evaluation,
            "History entry is not available at depth " + std::to_string(*depth));

    return context_->history[context_->history.size() - *depth];
}

expression::Expr Evaluator::evaluateIndexedHistory(
    std::span<const expression::Expr> arguments,
    bool input) {
    const auto index = historyIndexValue(arguments.front());
    if (!index)
        error::throwCalcError(
            error::CalcErrorType::Type,
            input ? "In index must be a non-zero integer" : "Out index must be a non-zero integer");

    if (!context_)
        error::throwCalcError(error::CalcErrorType::Evaluation, "Session history is not available");

    if (index->relative) {
        // In[-n]は入力履歴そのものを相対参照する。現在評価中の入力slotは除外するため、
        // In[-1] / @ は必ず直前の入力を指す。評価失敗した入力でもlower済みExprがあれば再評価できる。
        if (input) {
            const std::size_t previousInputCount = context_->inputs.empty()
                ? 0
                : context_->inputs.size() - 1;
            if (index->magnitude > previousInputCount)
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "Input entry is not available at relative index -"
                        + std::to_string(index->magnitude));

            const std::size_t entryIndex = previousInputCount - index->magnitude;
            if (!context_->inputs[entryIndex])
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "Input entry is not available at relative index -"
                        + std::to_string(index->magnitude));
            return *context_->inputs[entryIndex];
        }

        // Out[-n]は%/%%と同じく「成功した出力」の相対履歴を参照する。
        // これにより評価失敗した入力slotを挟んでも % == Out[-1] が常に成立する。
        if (index->magnitude > context_->history.size())
            error::throwCalcError(
                error::CalcErrorType::Evaluation,
                "Output entry is not available at relative index -"
                    + std::to_string(index->magnitude));
        return context_->history[context_->history.size() - index->magnitude];
    }

    const std::size_t number = index->magnitude;
    const auto& entries = input ? context_->inputs : context_->outputs;
    if (input && number == context_->inputs.size()
        && number <= context_->outputs.size() && !context_->outputs[number - 1])
        error::throwCalcError(
            error::CalcErrorType::Evaluation,
            "In cannot reference the input currently being evaluated");
    if (number > entries.size() || !entries[number - 1])
        error::throwCalcError(
            error::CalcErrorType::Evaluation,
            std::string{input ? "Input" : "Output"} + " entry is not available at index "
                + std::to_string(number));
    return *entries[number - 1];
}

void Evaluator::emitWarning(std::string_view code, std::string message) {
    if (!context_ || !context_->diagnostics)
        return;

    const auto duplicate = std::find_if(
        context_->diagnostics->begin(), context_->diagnostics->end(),
        [&](const EvaluationDiagnostic& diagnostic) {
            return diagnostic.severity == DiagnosticSeverity::Warning
                && diagnostic.code == code && diagnostic.message == message;
        });
    if (duplicate != context_->diagnostics->end())
        return;

    context_->diagnostics->push_back(EvaluationDiagnostic{
        DiagnosticSeverity::Warning, std::string{code}, std::move(message), std::nullopt});
}

void Evaluator::emitInfo(
    std::string_view code,
    std::string message,
    std::optional<expression::Expr> previousExpression) {
    if (!context_ || !context_->diagnostics)
        return;
    context_->diagnostics->push_back(EvaluationDiagnostic{
        DiagnosticSeverity::Info, std::string{code}, std::move(message),
        std::move(previousExpression)});
}

expression::Expr Evaluator::evaluateDefinitions() const {
    std::vector<expression::Expr> definitions;
    for (const auto& [symbol, value] : environment_.definitions())
        definitions.push_back(expression::Expr::call(
            registry_.symbol(BuiltinId::Set),
            {expression::Expr{symbol}, value}));

    if (userFunctions_) {
        for (const UserFunctionDefinition& definition : userFunctions_->definitions()) {
            std::vector<expression::Expr> signatureArguments;
            signatureArguments.reserve(definition.parameters.size() + 1);
            signatureArguments.emplace_back(definition.name);
            for (const expression::Symbol& parameter : definition.parameters)
                signatureArguments.emplace_back(parameter);
            definitions.push_back(expression::Expr::call(
                registry_.symbol(BuiltinId::SetDelayed),
                {expression::Expr::call(
                    expression::Symbol{builtins::names::functionSignature},
                    std::move(signatureArguments)), definition.body}));
        }
    }

    const std::size_t count = definitions.size();
    return expression::Expr::array({count}, std::move(definitions));
}

expression::Expr Evaluator::evaluateUndefine(std::span<const expression::Expr> arguments) {
    std::size_t changed = 0;
    for (const expression::Expr& argument : arguments) {
        if (!argument.isSymbol())
            error::throwCalcError(error::CalcErrorType::Type, "UnDef expects symbol arguments");
        const expression::Symbol symbol = argument.asSymbol();
        if (registry_.contains(symbol) || symbolRegistry_.isProtected(symbol))
            error::throwCalcError(
                error::CalcErrorType::Name,
                "Cannot undefine protected symbol: " + symbol.name());

        const bool removedVariable = environment_.erase(symbol);
        const bool removedFunction = userFunctions_ && userFunctions_->erase(symbol);
        if (removedVariable || removedFunction)
            ++changed;
    }

    if (changed != 0 && context_ && context_->definitionsChanged)
        *context_->definitionsChanged = true;
    return expression::Expr{numeric::Number{numeric::BigInt::parse(std::to_string(changed))}};
}

expression::Expr Evaluator::finalizeNumericalApproximation(
    const expression::Expr& value,
    std::size_t precisionDigits,
    bool warnOnFailure) {
    // free symbolを含む式は「数値化失敗」ではなく部分数値化の対象である。
    // protected constant/domainはfree symbolに数えず，ユーザー未知量だけを検出する。
    std::function<bool(const expression::Expr&)> containsUnknownSymbol;
    containsUnknownSymbol = [&](const expression::Expr& current) -> bool {
        if (current.isSymbol()) {
            const auto& symbol = current.asSymbol();
            return !symbolRegistry_.contains(symbol)
                && !registry_.contains(symbol)
                && mathematics_.findConstant(symbol) == nullptr;
        }
        if (current.isCall()) {
            for (const auto& argument : current.asCall().arguments)
                if (containsUnknownSymbol(argument))
                    return true;
            return false;
        }
        if (current.isArray()) {
            for (std::size_t i = 0; i < current.asArray().size(); ++i)
                if (containsUnknownSymbol(current.asArray().element(i)))
                    return true;
            return false;
        }
        if (current.isList()) {
            for (const auto& element : current.asList().elements)
                if (containsUnknownSymbol(element))
                    return true;
        }
        return false;
    };

    // InformationEnclosureの幅が有限precision入力そのものに由来するかを判別する。
    // exact式でも有限working precisionでは一時的にbranch cutを跨ぎ得るため，
    // その場合はInputInformation例外を即停止条件にせずguard増加で再試行する。
    std::function<bool(const expression::Expr&)> containsFinitePrecisionInput;
    containsFinitePrecisionInput = [&](const expression::Expr& current) -> bool {
        if (current.isDecimalApproximation() || current.isComplexDecimalApproximation())
            return true;
        if (current.isCall()) {
            const auto& call = current.asCall();
            if (const auto* definition = registry_.find(call.head);
                definition && definition->id == BuiltinId::NumericalApproximation)
                return true;
            for (const auto& argument : call.arguments)
                if (containsFinitePrecisionInput(argument))
                    return true;
            return false;
        }
        if (current.isArray()) {
            for (std::size_t i = 0; i < current.asArray().size(); ++i)
                if (containsFinitePrecisionInput(current.asArray().element(i)))
                    return true;
            return false;
        }
        if (current.isList()) {
            for (const auto& element : current.asList().elements)
                if (containsFinitePrecisionInput(element))
                    return true;
        }
        return false;
    };

    // Nはscalarだけでなく配列・SolutionSet・一般symbolic expressionへ部分的に作用する。
    // whole-expressionのcertificationを最優先し，失敗したときだけnumeric subpart traversalへ落とす。
    std::function<expression::Expr(const expression::Expr&, bool, bool)> approximate;
    approximate = [&](
        const expression::Expr& current,
        bool allowWarning,
        bool preserveExactInteger) -> expression::Expr {
        // precision-aware builtinが既に近似値を返した場合、外側Nがより低い桁を要求するなら
        // certified enclosureから安全に丸め直す。より高い桁は元情報以上に増やせないため保持する。
        if (current.isDecimalApproximation())
            return expression::Expr{reduceApproximationPrecision(
                current.asDecimalApproximation(), precisionDigits)};
        if (current.isComplexDecimalApproximation())
            return reduceApproximationPrecision(
                current.asComplexDecimalApproximation(), precisionDigits);

        // Boolean/StringやInfinity・domain symbol・自由変数は数値ではないため，Nの失敗ではない。
        // Pi/E/PhiのようなMathRegistry定数だけは後段のcertified evaluatorへ送る。
        if (current.isBoolean() || current.isString())
            return current;
        if (current.isSymbol() && mathematics_.findConstant(current.asSymbol()) == nullptr)
            return current;

        if (current.isArray()) {
            const auto& array = current.asArray();
            if (array.storageKind() == expression::ArrayStorageKind::DecimalApproximation) {
                std::vector<numeric::DecimalApproximation> values;
                values.reserve(array.size());
                for (std::size_t i = 0; i < array.size(); ++i)
                    values.push_back(reduceApproximationPrecision(array.decimalAt(i), precisionDigits));
                return expression::Expr::decimalArray(array.shape, std::move(values));
            }
            if (array.storageKind() == expression::ArrayStorageKind::ComplexDecimalApproximation) {
                std::vector<numeric::ComplexDecimalApproximation> values;
                values.reserve(array.size());
                for (std::size_t i = 0; i < array.size(); ++i) {
                    const expression::Expr reduced = reduceApproximationPrecision(
                        array.complexDecimalAt(i), precisionDigits);
                    values.push_back(reduced.asComplexDecimalApproximation());
                }
                return expression::Expr::complexDecimalArray(array.shape, std::move(values));
            }
            std::vector<expression::Expr> elements;
            elements.reserve(array.size());
            for (std::size_t i = 0; i < array.size(); ++i)
                elements.push_back(approximate(array.element(i), allowWarning, false));
            return expression::Expr::array(array.shape, std::move(elements));
        }
        if (current.isList()) {
            const auto& list = current.asList();
            std::vector<expression::Expr> elements;
            elements.reserve(list.elements.size());
            for (const expression::Expr& element : list.elements)
                elements.push_back(approximate(element, allowWarning, false));
            return expression::braceValue(std::move(elements));
        }

        // SolutionSetは未知変数を数値化せず，各bindingの右辺だけへNを作用させる。
        // 数値化できないparameterized solutionはexactのまま保持するため，binding内部の失敗は警告しない。
        if (current.isSolutionSet()) {
            const solver::SolutionSet& solutions = current.asSolutionSet();
            std::vector<solver::SolverVariable> variables{
                solutions.variables().begin(), solutions.variables().end()};

            const auto approximateBranch = [&](const solver::SolutionBranch& source) {
                solver::SolutionBranch result = source;
                for (solver::SolutionBinding& binding : result.bindings)
                    binding.value = approximate(binding.value, false, false);
                return result;
            };

            solver::SolutionSet transformed = [&]() {
                switch (solutions.kind()) {
                case solver::SolutionSetKind::Empty:
                    return solver::SolutionSet::empty(std::move(variables));
                case solver::SolutionSetKind::Finite: {
                    std::vector<solver::SolutionBranch> branches;
                    branches.reserve(solutions.branches().size());
                    for (const solver::SolutionBranch& branch : solutions.branches())
                        branches.push_back(approximateBranch(branch));
                    return solver::SolutionSet::finite(std::move(variables), std::move(branches));
                }
                case solver::SolutionSetKind::Universal:
                    return solver::SolutionSet::universal(
                        std::move(variables), solutions.conditions());
                case solver::SolutionSetKind::Conditional: {
                    std::vector<solver::SolutionCase> cases{
                        solutions.cases().begin(), solutions.cases().end()};
                    for (solver::SolutionCase& solutionCase : cases) {
                        if (solutionCase.outcome != solver::SolutionSetKind::Finite)
                            continue;
                        for (solver::SolutionBranch& branch : solutionCase.branches)
                            branch = approximateBranch(branch);
                    }
                    return solver::SolutionSet::conditional(std::move(variables), std::move(cases));
                }
                case solver::SolutionSetKind::Unresolved:
                    return solver::SolutionSet::unresolved(
                        std::move(variables), solutions.conditions());
                }
                throw std::logic_error("Unknown SolutionSet kind");
            }();

            if (solutions.kind() == solver::SolutionSetKind::Finite
                || solutions.kind() == solver::SolutionSetKind::Conditional)
                transformed = transformed.withAdditionalConditions(solutions.conditions());
            return expression::Expr::solutionSet(std::move(transformed));
        }

        // SeriesDataの指数格子metadataは構造情報なのでexact整数のまま保持する。
        // centerと係数だけへNを作用させ，変換後もparseSeriesData可能な形を維持する。
        if (const auto series = symbolic::parseSeriesData(current, registry_)) {
            symbolic::SeriesData transformed = *series;
            transformed.center = approximate(transformed.center, false, false);
            for (expression::Expr& coefficient : transformed.coefficients)
                coefficient = approximate(coefficient, false, false);
            for (auto& layer : transformed.logarithmicCoefficients)
                for (expression::Expr& coefficient : layer)
                    coefficient = approximate(coefficient, false, false);
            return symbolic::makeSeriesData(std::move(transformed), registry_);
        }

        // casesは条件をexactのまま保持し，各branchの値だけへNを作用させる。
        if (current.isCall()
            && current.asCall().head.sameIdentity(registry_.symbol(BuiltinId::Cases))) {
            std::vector<expression::Expr> branches;
            branches.reserve(current.asCall().arguments.size());
            for (const expression::Expr& branchExpression : current.asCall().arguments) {
                if (!branchExpression.isCall()
                    || !branchExpression.asCall().head.sameIdentity(
                        registry_.symbol(BuiltinId::CaseBranch))
                    || branchExpression.asCall().arguments.empty()
                    || branchExpression.asCall().arguments.size() > 2)
                    return current;
                const auto& branch = branchExpression.asCall().arguments;
                std::vector<expression::Expr> branchArguments{
                    approximate(branch[0], false, false)};
                if (branch.size() == 2)
                    branchArguments.push_back(branch[1]);
                branches.push_back(expression::Expr::call(
                    registry_.symbol(BuiltinId::CaseBranch), std::move(branchArguments)));
            }
            return expression::Expr::call(
                registry_.symbol(BuiltinId::Cases), std::move(branches));
        }

        // UnitAppliedは単位文字列そのものを数値化せず、値の部分だけへNを作用させる。
        // arg等が返す明示Radも、この経路で近似値と単位を両立できる。
        if (current.isCall()) {
            const auto& currentCall = current.asCall();
            const auto* definition = registry_.find(currentCall.head);
            if (definition && definition->id == BuiltinId::UnitApplied
                && currentCall.arguments.size() == 2 && currentCall.arguments[1].isString()) {
                return expression::Expr::call(currentCall.head, {
                    approximate(currentCall.arguments[0], allowWarning, false),
                    currentCall.arguments[1]
                });
            }
        }

        // symbolic callを部分数値化するとき，exact integerは係数だけでなく指数・branch番号・
        // 個数などの構造parameterにも使われる。意味を2.0へ変えないようatomのまま保持する。
        // N[2,p]のように値そのものを要求された場合は従来どおりapproximationへ変換する。
        if (current.isNumber() && preserveExactInteger) {
            const numeric::Number& number = current.asNumber();
            if (number.isReal() && number.asReal().isInteger())
                return current;
        }

        // exactなNumberだけは区間算法へ送る必要がない。有限小数なら必要最小桁で表示し、
        // 循環小数だけ要求桁へ丸めるという従来のNの表示規則を保つ。
        if (current.isNumber()) {
            const numeric::Number& number = current.asNumber();
            if (number.isReal())
                return expression::Expr{numeric::DecimalApproximation::fromRealSignificant(
                    number.asReal(), precisionDigits)};

            const auto& complex = number.asComplex();
            return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
                numeric::DecimalApproximation::fromRealSignificant(complex.real, precisionDigits),
                numeric::DecimalApproximation::fromRealSignificant(complex.imaginary, precisionDigits),
                complex.real.isZero(),
                complex.imaginary.isZero())};
        }

        approximation::CertifiedEvaluator certified{registry_, mathematics_, angleSemantics_};
        approximation::ApproximationContext context{precisionDigits};
        std::string lastPrecisionFailure;
        for (std::size_t refinement = 0;
             refinement < maximumNumericalApproximationRefinements;
             ++refinement) {
            consumeEvaluationBudget(EvaluationResource::CertifiedRefinement);
            try {
                const std::size_t bits = context.workingBinaryBits();
                // finite-precision入力を含む式では，値そのもののCertifiedEnclosureより先に
                // InformationEnclosureでdomain/branchを検査する。exact入力だけなら両者は
                // 同一なので，重い特殊函数を二度評価せずCertifiedEnclosureをそのまま再利用する。
                const bool finitePrecisionInput = containsFinitePrecisionInput(current);
                std::optional<approximation::CertifiedValue> information;
                if (finitePrecisionInput)
                    information = certified.enclose(
                        current, bits, approximation::CertifiedEvaluator::EnclosureKind::Information);
                const auto enclosed = certified.enclose(
                    current, bits, approximation::CertifiedEvaluator::EnclosureKind::Certified);
                if (!finitePrecisionInput && enclosed)
                    information = *enclosed;
                if (!enclosed || !information) {
                    // whole expressionをcertifyできない場合だけ，通常評価型のcallのnumeric subpartへNを作用させる。
                    // HoldAll/HoldFirst系は変数指定やiterator等の構文的引数を持つため勝手に書き換えない。
                    if (current.isCall()) {
                        const auto& call = current.asCall();
                        const auto* definition = registry_.find(call.head);
                        const bool structural = !definition
                            || definition->argumentEvaluation == ArgumentEvaluation::All;
                        if (structural) {
                            std::vector<expression::Expr> arguments;
                            arguments.reserve(call.arguments.size());
                            bool changed = false;
                            for (const auto& argument : call.arguments) {
                                expression::Expr transformed = approximate(argument, false, true);
                                changed = changed || !(transformed == argument);
                                arguments.push_back(std::move(transformed));
                            }
                            if (changed) {
                                expression::Expr rebuilt =
                                    expression::Expr::rebuildCall(call, std::move(arguments));
                                // numeric childだけが変わっても，親が閉じた未対応式のままなら
                                // 「部分的に何か変わった」ことを成功扱いしない。再度whole-expressionを試し，
                                // free symbolがなければ最終的に適切なN warningへ落とす。
                                return approximate(rebuilt, allowWarning, preserveExactInteger);
                            }
                        }
                    }

                    // 自由記号を含む式はpartial Nとして正常に保持する。完全に数値閉包なのに
                    // backendが値を作れない場合だけgeneric warningを出す。
                    if (allowWarning && !containsUnknownSymbol(current))
                        emitWarning("N::unevaluated",
                            "N could not certify a numerical value for part of the expression; it remains unevaluated");
                    return current;
                }
                if (const auto decimal = approximation::finalizeCertifiedApproximation(
                        *enclosed, *information, precisionDigits))
                    return *decimal;
            }
            catch (const approximation::CertifiedBackendUnsupported& exception) {
                if (allowWarning)
                    emitWarning("N::unsupported",
                        std::string{exception.what()} + "; the expression remains unevaluated");
                return current;
            }
            catch (const approximation::PrecisionInsufficient& exception) {
                // 数学的domain errorではなく，現在の区間幅では分岐・非零性等を証明できない。
                // finite-precision入力自身のInformationEnclosureが原因ならguard増加では
                // 改善しないため，局所retryを即座に終了する。
                lastPrecisionFailure = exception.what();
                if (!exception.refinable() && containsFinitePrecisionInput(current))
                    break;
            }

            if (refinement + 1 < maximumNumericalApproximationRefinements)
                context.setGuardDigits(nextApproximationGuardDigits(context.guardDigits()));
        }

        if (allowWarning && !containsUnknownSymbol(current)) {
            std::string message =
                "N could not certify the requested decimal precision after bounded refinement";
            if (!lastPrecisionFailure.empty())
                message += ": " + lastPrecisionFailure;
            message += "; the expression remains unevaluated";
            emitWarning("N::precision", std::move(message));
        }
        return current;
    };

    return approximate(value, warnOnFailure, false);
}

const approximation::ApproximationContext* Evaluator::currentApproximationContext() const noexcept {
    return approximationContexts_.empty() ? nullptr : &approximationContexts_.back();
}

std::optional<source::SourceReference> Evaluator::originOf(
    const expression::Expr& expression) const {
    if (!origins_)
        return std::nullopt;

    return origins_->find(expression);
}

} // namespace mmcal::evaluation
