// 非再帰タスクスタック式評価器
#include "evaluator.hpp"

#include "expression/array_utils.hpp"

#include "approximation/approximation_context.hpp"
#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "evaluation/iterator_spec.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/integer_algorithms.hpp"
#include "simplification/simplifier.hpp"

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
    std::vector<std::size_t> shape;
    std::size_t elementCount = 0;
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

struct EnterUserFunctionTask final {
    expression::Expr expression;
    const UserFunctionDefinition* definition = nullptr;
    const expression::OriginMap* origins = nullptr;
    std::size_t childDepth = 1;
};

struct FinishUserFunctionTask final {};

// N[...] は第1引数を先にexact評価しない。precisionだけを確定してから子式を評価し、
// precision-aware builtinへ要求精度を伝播した後、最後に従来どおりcertified decimalへ落とす。
struct BeginNumericalApproximationTask final {
    expression::Expr expression;
    const expression::OriginMap* origins = nullptr;
    std::size_t childDepth = 1;
};

struct FinishNumericalApproximationTask final {
    std::size_t fractionalDigits = approximation::ApproximationContext::defaultDecimalDigits;
};

using EvaluationTask = std::variant<
    EvaluateTask,
    PushResultTask,
    FinishSymbolTask,
    BuildArrayTask,
    BuildListTask,
    DispatchBuiltinTask,
    IfConditionTask,
    EnterUserFunctionTask,
    FinishUserFunctionTask,
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

    tasks.emplace_back(BuildArrayTask{argument, spec->shape, spec->elements.size(), origins});
    tasks.emplace_back(EvaluateTask{spec->elements[2], origins, depth});
    tasks.emplace_back(EvaluateTask{spec->elements[1], origins, depth});
    tasks.emplace_back(PushResultTask{spec->elements[0]});
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

[[nodiscard]] std::size_t nextApproximationGuardDigits(std::size_t current) {
    const std::size_t growth = std::max<std::size_t>(8, current / 2);
    if (growth > std::numeric_limits<std::size_t>::max() - current)
        throw std::overflow_error("Approximation precision is too large");
    return current + growth;
}

[[nodiscard]] bool intervalIsExactZero(
    const approximation::RealInterval& interval) noexcept {
    return interval.isPoint() && interval.lower().isZero();
}

[[nodiscard]] std::optional<expression::Expr> certifiedDecimalExpression(
    const approximation::CertifiedValue& value,
    std::size_t fractionalDigits) {
    if (value.isReal()) {
        const auto decimal = numeric::DecimalApproximation::fromCertifiedInterval(
            value.asReal().lower().toRational(),
            value.asReal().upper().toRational(),
            fractionalDigits);
        return decimal ? std::optional<expression::Expr>{expression::Expr{*decimal}}
                       : std::nullopt;
    }

    const approximation::ComplexInterval& complex = value.asComplex();
    const auto real = numeric::DecimalApproximation::fromCertifiedInterval(
        complex.real().lower().toRational(),
        complex.real().upper().toRational(),
        fractionalDigits);
    const auto imaginary = numeric::DecimalApproximation::fromCertifiedInterval(
        complex.imaginary().lower().toRational(),
        complex.imaginary().upper().toRational(),
        fractionalDigits);
    if (!real || !imaginary)
        return std::nullopt;

    return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
        *real,
        *imaginary,
        intervalIsExactZero(complex.real()),
        intervalIsExactZero(complex.imaginary()))};
}

[[nodiscard]] numeric::DecimalApproximation reduceApproximationDigits(
    const numeric::DecimalApproximation& value,
    std::size_t fractionalDigits) {
    if (fractionalDigits >= value.requestedFractionalDigits())
        return value;

    // exact point由来なら従来の有限小数最小表記を維持する。
    if (value.enclosureIsPoint())
        return numeric::DecimalApproximation::fromReal(
            numeric::RealNumber{value.certifiedLower()}, fractionalDigits);

    if (const auto rounded = numeric::DecimalApproximation::fromCertifiedInterval(
        value.certifiedLower(), value.certifiedUpper(), fractionalDigits))
        return *rounded;

    // 元の保証区間が粗く、より低い桁への丸め境界を跨ぐ特殊caseでは、
    // 情報を捨てて誤った桁へ丸めず既存の近似値を保持する。
    return value;
}

[[nodiscard]] expression::Expr reduceApproximationDigits(
    const numeric::ComplexDecimalApproximation& value,
    std::size_t fractionalDigits) {
    const auto real = reduceApproximationDigits(value.real(), fractionalDigits);
    const auto imaginary = reduceApproximationDigits(value.imaginary(), fractionalDigits);
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
    return evaluateMachine(expression, nullptr, nullptr);
}

expression::Expr Evaluator::evaluate(
    const expression::Expr& expression,
    const expression::OriginMap& origins) {
    return evaluateMachine(expression, &origins, nullptr);
}

expression::Expr Evaluator::evaluate(
    const expression::Expr& expression,
    const expression::OriginMap& origins,
    EvaluationContext context) {
    return evaluateMachine(expression, &origins, &context);
}

void Evaluator::setDepthLimit(std::size_t limit) {
    if (limit == 0)
        throw std::invalid_argument("Evaluation depth limit must be greater than zero");

    depthLimit_ = limit;
}

std::size_t Evaluator::depthLimit() const noexcept {
    return depthLimit_;
}

void Evaluator::reseedRandomFromEntropy() {
    static_cast<void>(randomEngine_.reseedFromEntropy());
}

expression::Expr Evaluator::evaluateMachine(
    const expression::Expr& expression,
    const expression::OriginMap* origins,
    const EvaluationContext* context) {
    const std::size_t initialLocalDepth = environment_.localDepth();
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
            EvaluationTask task = std::move(tasks.back());
            tasks.pop_back();
            const TaskSource source = taskSource(task);
            origins_ = source.origins;

            try {
                std::visit(Overloaded{
                    [&](const EvaluateTask& current) {
                        if (current.depth > depthLimit_)
                            error::throwCalcError(
                                error::CalcErrorType::Evaluation,
                                "Evaluation depth limit exceeded");

                        switch (current.expression.kind()) {
                        case expression::ExprKind::Number:
                        case expression::ExprKind::DecimalApproximation:
                        case expression::ExprKind::ComplexDecimalApproximation:
                        case expression::ExprKind::Boolean:
                        case expression::ExprKind::String:
                        case expression::ExprKind::SolutionSet:
                            results.push_back(current.expression);
                            return;

                        case expression::ExprKind::Symbol: {
                            const expression::Symbol symbol = current.expression.asSymbol();
                            const expression::Expr* binding = environment_.find(symbol);
                            if (!binding) {
                                // 暫定方針: 未束縛Symbolは自由記号としてそのまま通す。
                                // Simplify/Expand/Factor/Collect/SolveをCLIから直接検証できるようにするため。
                                // 将来SymbolicEvaluationContextを導入したら、通常評価でNameErrorに戻すモードと明示的なsymbolic modeをここで分離できる。
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
                            tasks.emplace_back(BuildArrayTask{
                                current.expression,
                                array.shape,
                                array.elements.size(),
                                current.origins
                            });
                            for (auto iterator = array.elements.rbegin(); iterator != array.elements.rend(); ++iterator)
                                tasks.emplace_back(EvaluateTask{*iterator, current.origins, current.depth + 1});
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

                                    const bool held = definition->argumentEvaluation == ArgumentEvaluation::HoldAll
                                        || (definition->argumentEvaluation == ArgumentEvaluation::HoldFirst && index == 0)
                                        || (definition->argumentEvaluation == ArgumentEvaluation::HoldFirstTwo && index < 2)
                                        || (definition->argumentEvaluation == ArgumentEvaluation::HoldFirstAndIteratorSpec && index == 0);
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
                        std::vector<expression::Expr> elements = takeResults(results, current.elementCount);
                        results.push_back(expression::rebuildEvaluatedArray(current.shape, std::move(elements)));
                    },
                    [&](const BuildListTask& current) {
                        std::vector<expression::Expr> elements = takeResults(results, current.elementCount);
                        results.push_back(expression::braceValue(std::move(elements)));
                    },
                    [&](const DispatchBuiltinTask& current) {
                        const expression::CallExpr& call = current.expression.asCall();
                        std::vector<expression::Expr> arguments = takeResults(results, call.arguments.size());

                        // 評価後に0になった除数でも、元の除数範囲をエラー位置として使う。
                        if (current.definition->id == BuiltinId::Divide
                            && arguments[1].isNumber()
                            && arguments[1].asNumber().isZero()) {
                            if (current.origins) {
                                if (const auto origin = current.origins->find(call.arguments[1]))
                                    error::throwCalcError(
                                        error::CalcErrorType::Domain,
                                        "Division by zero",
                                        *origin);
                            }
                            error::throwCalcError(
                                error::CalcErrorType::Domain,
                                "Division by zero");
                        }

                        expression::Expr dispatched = dispatchBuiltin(*current.definition, call, arguments);

                        // 組み込み函数が局所的に評価した後、共通Simplifierへ一度通す。
                        // これにより Exp[Log[Pi]] のような複数函数にまたがる安全な書換えを各builtinへ重複実装せず、同じ数学知識へ集約する。
                        const bool pureForPostSimplification =
                            current.definition->id != BuiltinId::Set
                            && current.definition->id != BuiltinId::SetDelayed
                            && current.definition->id != BuiltinId::If
                            && current.definition->id != BuiltinId::History
                            && current.definition->id != BuiltinId::InputHistory
                            && current.definition->id != BuiltinId::OutputHistory
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
                            const simplification::SimplificationContext simplificationContext{
                                registry_, mathematics_, angleSemantics_};
                            dispatched = simplification::Simplifier{}.simplify(
                                dispatched, simplificationContext);
                        }

                        // In[n]は保存した入力Exprを「貼り戻す」意味とし、取得したExprを現在の環境で通常評価する。Out[n]は保存済み結果なので再評価しない。
                        if (current.definition->id == BuiltinId::InputHistory) {
                            tasks.emplace_back(EvaluateTask{
                                std::move(dispatched), nullptr, current.depth + 1});
                            return;
                        }
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
                    [&](const BeginNumericalApproximationTask& current) {
                        constexpr std::size_t defaultFractionalDigits = 16;
                        const expression::CallExpr& call = current.expression.asCall();
                        std::size_t fractionalDigits = defaultFractionalDigits;

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
                            fractionalDigits = *requestedDigits;
                        }

                        approximationContexts_.emplace_back(fractionalDigits);
                        tasks.emplace_back(FinishNumericalApproximationTask{fractionalDigits});
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

                        std::vector<expression::Expr> valueResult = takeResults(results, 1);
                        expression::Expr approximated = finalizeNumericalApproximation(
                            valueResult.front(), current.fractionalDigits);
                        approximationContexts_.pop_back();
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
    std::size_t fractionalDigits) {
    // Nはscalarだけでなく配列へも要素単位に作用する。
    // FFT/行列等のexact配列を表示用近似へ落とす際に、各builtinが独自のdigits引数を持つ必要をなくす。
    std::function<expression::Expr(const expression::Expr&)> approximate;
    approximate = [&](const expression::Expr& current) -> expression::Expr {
        // precision-aware builtinが既に近似値を返した場合、外側Nがより低い桁を要求するなら
        // certified enclosureから安全に丸め直す。より高い桁は元情報以上に増やせないため保持する。
        if (current.isDecimalApproximation())
            return expression::Expr{reduceApproximationDigits(
                current.asDecimalApproximation(), fractionalDigits)};
        if (current.isComplexDecimalApproximation())
            return reduceApproximationDigits(
                current.asComplexDecimalApproximation(), fractionalDigits);

        if (current.isArray()) {
            const auto& array = current.asArray();
            std::vector<expression::Expr> elements;
            elements.reserve(array.elements.size());
            for (const expression::Expr& element : array.elements)
                elements.push_back(approximate(element));
            return expression::Expr::array(array.shape, std::move(elements));
        }
        if (current.isList()) {
            const auto& list = current.asList();
            std::vector<expression::Expr> elements;
            elements.reserve(list.elements.size());
            for (const expression::Expr& element : list.elements)
                elements.push_back(approximate(element));
            return expression::braceValue(std::move(elements));
        }

        // UnitAppliedは単位文字列そのものを数値化せず、値の部分だけへNを作用させる。
        // arg等が返す明示Radも、この経路で近似値と単位を両立できる。
        if (current.isCall()) {
            const auto& currentCall = current.asCall();
            const auto* definition = registry_.find(currentCall.head);
            if (definition && definition->id == BuiltinId::UnitApplied
                && currentCall.arguments.size() == 2 && currentCall.arguments[1].isString()) {
                return expression::Expr::call(currentCall.head, {
                    approximate(currentCall.arguments[0]),
                    currentCall.arguments[1]
                });
            }
        }

        // exactなNumberだけは区間算法へ送る必要がない。有限小数なら必要最小桁で表示し、
        // 循環小数だけ要求桁へ丸めるという従来のNの表示規則を保つ。
        if (current.isNumber()) {
            const numeric::Number& number = current.asNumber();
            if (number.isReal())
                return expression::Expr{numeric::DecimalApproximation::fromReal(
                    number.asReal(), fractionalDigits)};

            const auto& complex = number.asComplex();
            return expression::Expr{numeric::ComplexDecimalApproximation::fromComponents(
                numeric::DecimalApproximation::fromReal(complex.real, fractionalDigits),
                numeric::DecimalApproximation::fromReal(complex.imaginary, fractionalDigits),
                complex.real.isZero(),
                complex.imaginary.isZero())};
        }

        approximation::CertifiedEvaluator certified{registry_, mathematics_, angleSemantics_};
        approximation::ApproximationContext context{fractionalDigits};
        for (;;) {
            try {
                const auto enclosed = certified.enclose(current, context.workingBinaryBits());
                if (!enclosed) {
                    emitWarning("N::unevaluated",
                        "N could not certify a numerical value for part of the expression; it remains unevaluated");
                    return current;
                }
                if (const auto decimal = certifiedDecimalExpression(*enclosed, fractionalDigits))
                    return *decimal;
            }
            catch (const approximation::PrecisionInsufficient&) {
                // 数学的domain errorではなく、現在の区間幅では分岐を証明できない。
            }
            context.setGuardDigits(nextApproximationGuardDigits(context.guardDigits()));
        }
    };

    return approximate(value);
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
