// 非再帰タスクスタック式評価器
#include "evaluator.hpp"

#include "mathematics/angle.hpp"
#include "mathematics/assumption_parser.hpp"

#include "builtins/arithmetic.hpp"
#include "builtins/approximation_utilities.hpp"
#include "builtins/aggregate.hpp"
#include "builtins/statistics.hpp"
#include "builtins/signal_processing.hpp"
#include "builtins/array_vector.hpp"
#include "builtins/comparison.hpp"
#include "builtins/combinatorics.hpp"
#include "builtins/discrete_math.hpp"
#include "builtins/elementary_utilities.hpp"
#include "builtins/complex_functions.hpp"
#include "builtins/hyperbolic.hpp"
#include "builtins/linear_algebra.hpp"
#include "builtins/numerical_calculus.hpp"
#include "builtins/trigonometric.hpp"
#include "builtins/stable_elementary.hpp"
#include "builtins/special_functions.hpp"
#include "builtins/random_functions.hpp"
#include "simplification/simplifier.hpp"
#include "builtins/transcendental.hpp"
#include "approximation/approximation_context.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/certification_error.hpp"
#include "builtins/names.hpp"
#include "solver/polynomial_solver.hpp"
#include "solver/solve_constraints.hpp"
#include "solver/transcendental_solver.hpp"
#include "symbolic/algebra_transforms.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/integration.hpp"
#include "symbolic/limit.hpp"
#include "simplification/full_simplifier.hpp"
#include "error/error_message.hpp"
#include "evaluation/iterator_spec.hpp"
#include "numeric/integer_algorithms.hpp"

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

using EvaluationTask = std::variant<
    EvaluateTask,
    PushResultTask,
    FinishSymbolTask,
    BuildArrayTask,
    DispatchBuiltinTask,
    IfConditionTask,
    EnterUserFunctionTask,
    FinishUserFunctionTask>;

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
        [](const DispatchBuiltinTask& value) noexcept {
            return TaskSource{&value.expression, value.origins};
        },
        [](const IfConditionTask& value) noexcept {
            return TaskSource{&value.expression, value.origins};
        },
        [](const EnterUserFunctionTask& value) noexcept {
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

[[nodiscard]] bool containsBuiltinCall(
    const expression::Expr& root,
    const expression::Symbol& head) {
    std::vector<const expression::Expr*> pending{&root};
    while (!pending.empty()) {
        const expression::Expr* current = pending.back();
        pending.pop_back();
        if (current->isCall()) {
            const auto& call = current->asCall();
            if (call.head.sameIdentity(head))
                return true;
            for (const expression::Expr& argument : call.arguments)
                pending.push_back(&argument);
        }
        else if (current->isArray()) {
            for (const expression::Expr& element : current->asArray().elements)
                pending.push_back(&element);
        }
    }
    return false;
}

[[nodiscard]] bool containsUnresolvedSolution(const solver::SolutionSet& solutions) {
    if (solutions.kind() == solver::SolutionSetKind::Unresolved)
        return true;
    if (solutions.kind() != solver::SolutionSetKind::Conditional)
        return false;
    return std::any_of(solutions.cases().begin(), solutions.cases().end(),
        [](const solver::SolutionCase& item) {
            return item.outcome == solver::SolutionSetKind::Unresolved;
        });
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
                        results.push_back(expression::Expr::array(current.shape, std::move(elements)));
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

expression::Expr Evaluator::dispatchBuiltin(
    const BuiltinDefinition& definition,
    const expression::CallExpr& call,
    std::span<const expression::Expr> arguments) {
    switch (definition.id) {
    case BuiltinId::Add:
        return builtins::evaluateAdd(arguments, registry_);
    case BuiltinId::Subtract:
        return builtins::evaluateSubtract(arguments, registry_);
    case BuiltinId::Multiply:
        return builtins::evaluateMultiply(arguments, registry_);
    case BuiltinId::Divide:
        return builtins::evaluateDivide(arguments, registry_);
    case BuiltinId::Power:
        return builtins::evaluatePower(arguments, registry_, mathematics_);
    case BuiltinId::Negate:
        return builtins::evaluateNegate(arguments, registry_);
    case BuiltinId::Factorial:
        return builtins::evaluateFactorial(arguments);
    case BuiltinId::Derivative: {
        if (arguments.size() < 2)
            error::throwCalcError(error::CalcErrorType::Type,
                "D expects an expression followed by one or more derivative specifications");

        expression::Expr result = arguments[0];
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            expression::Symbol variable;
            std::uint64_t order = 1;
            if (arguments[i].isSymbol()) {
                variable = arguments[i].asSymbol();
            }
            else if (arguments[i].isArray()) {
                const auto& spec = arguments[i].asArray();
                if (spec.shape.size() != 1 || spec.shape[0] != 2 || spec.elements.size() != 2
                    || !spec.elements[0].isSymbol() || !spec.elements[1].isNumber()
                    || !spec.elements[1].asNumber().isReal()
                    || !spec.elements[1].asNumber().asReal().isInteger()) {
                    error::throwCalcError(error::CalcErrorType::Type,
                        "D derivative specification must be a symbol or {symbol, nonnegative integer}");
                }
                const auto parsed = numeric::tryToUint64(
                    spec.elements[1].asNumber().asReal().asInteger());
                if (!parsed)
                    error::throwCalcError(error::CalcErrorType::Domain,
                        "D derivative order must be a nonnegative integer that fits in uint64");
                if (*parsed > 4096)
                    error::throwCalcError(error::CalcErrorType::Overflow,
                        "D derivative order is too large");
                variable = spec.elements[0].asSymbol();
                order = *parsed;
            }
            else {
                error::throwCalcError(error::CalcErrorType::Type,
                    "D derivative specification must be a symbol or {symbol, nonnegative integer}");
            }

            for (std::uint64_t derivative = 0; derivative < order; ++derivative)
                result = symbolic::differentiateExpression(
                    result, variable, registry_, mathematics_, angleSemantics_);
        }

        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Derivative)))
            emitWarning("D::unevaluated",
                "D could not fully evaluate the derivative; unevaluated D[...] remains");
        return result;
    }
    case BuiltinId::SymbolicIntegral: {
        if (arguments.size() < 2 || arguments.size() > 3)
            error::throwCalcError(error::CalcErrorType::Type,
                "integrate expects integrate[expression, variable, optional assumptions] or integrate[expression, {variable, lower, upper}, optional assumptions]");

        mathematics::AssumptionSet assumptions;
        if (arguments.size() == 3)
            assumptions = mathematics::parseAssumptions(arguments[2], registry_, mathematics_);

        std::optional<expression::Expr> result;
        if (arguments[1].isSymbol()) {
            result = symbolic::integrateExpression(
                arguments[0], arguments[1].asSymbol(), registry_, mathematics_,
                angleSemantics_, assumptions);
        }
        else if (const auto iterator = parseRangeIteratorSpec(arguments[1])) {
            const auto* infinity = symbolRegistry_.find("Infinity");
            if (!infinity)
                error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
            result = symbolic::integrateExpression(
                arguments[0], iterator->variable, iterator->lower, iterator->upper,
                registry_, mathematics_, angleSemantics_, infinity->symbol, assumptions);
        }
        else {
            error::throwCalcError(error::CalcErrorType::Type,
                "integrate expects a symbol or {variable, lower, upper} as the second argument");
        }

        if (containsBuiltinCall(*result, registry_.symbol(BuiltinId::SymbolicIntegral)))
            emitWarning("integrate::unevaluated",
                "integrate could not fully prove the symbolic antiderivative or definite integral; unevaluated integrate[...] remains");
        return *result;
    }
    case BuiltinId::Limit: {
        if (arguments.size() < 3 || arguments.size() > 4 || !arguments[1].isSymbol())
            error::throwCalcError(error::CalcErrorType::Type,
                "limit expects limit[expression, variable, point] or limit[expression, variable, point, direction]");

        symbolic::LimitDirection direction = symbolic::LimitDirection::TwoSided;
        if (arguments.size() == 4) {
            if (!arguments[3].isNumber() || !arguments[3].asNumber().isReal()
                || !arguments[3].asNumber().asReal().isInteger())
                error::throwCalcError(error::CalcErrorType::Type,
                    "limit direction must be -1 for left or 1 for right");
            const auto& value = arguments[3].asNumber().asReal().asInteger();
            if (value == numeric::BigInt{-1})
                direction = symbolic::LimitDirection::Left;
            else if (value == numeric::BigInt{1})
                direction = symbolic::LimitDirection::Right;
            else
                error::throwCalcError(error::CalcErrorType::Domain,
                    "limit direction must be -1 for left or 1 for right");
        }

        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
        expression::Expr result = symbolic::limitExpression(
            arguments[0], arguments[1].asSymbol(), arguments[2], direction,
            registry_, mathematics_, angleSemantics_, infinity->symbol);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Limit)))
            emitWarning("limit::unevaluated",
                "limit could not prove the requested limit; unevaluated limit[...] remains");
        return result;
    }
    case BuiltinId::Floor:
        return builtins::evaluateFloor(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Ceil:
        return builtins::evaluateCeil(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Trunc:
        return builtins::evaluateTrunc(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Round:
        return builtins::evaluateRound(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Frac:
        return builtins::evaluateFrac(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Gcd:
        return builtins::evaluateGcd(arguments, registry_);
    case BuiltinId::Lcm:
        return builtins::evaluateLcm(arguments, registry_);
    case BuiltinId::Mod:
        return builtins::evaluateMod(arguments, registry_);
    case BuiltinId::Rem:
        return builtins::evaluateRem(arguments, registry_);
    case BuiltinId::Quotient:
        return builtins::evaluateQuotient(arguments, registry_);
    case BuiltinId::Permutation:
        return builtins::evaluatePermutation(arguments, registry_);
    case BuiltinId::Combination:
        return builtins::evaluateCombination(arguments, registry_);
    case BuiltinId::Fibonacci:
        return builtins::evaluateFibonacci(arguments, registry_);
    case BuiltinId::DiscreteFourierTransform:
        return builtins::evaluateDft(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::FastFourierTransform:
        return builtins::evaluateFft(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::InverseFourierTransform:
        return builtins::evaluateIfft(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Convolution:
        return builtins::evaluateConvolution(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Transpose:
        return builtins::evaluateTranspose(arguments, registry_);
    case BuiltinId::MatrixAdd:
        return builtins::evaluateMatrixAdd(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::MatrixMultiply:
        return builtins::evaluateMatrixMultiply(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Determinant:
        return builtins::evaluateDeterminant(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Inverse:
        return builtins::evaluateMatrixInverse(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Rref: {
        expression::Expr result = builtins::evaluateRref(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Rref)))
            emitWarning("rref::unevaluated",
                "rref could not determine the symbolic pivots; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::Rank: {
        expression::Expr result = builtins::evaluateMatrixRank(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Rank)))
            emitWarning("rank::unevaluated",
                "rank could not determine the symbolic pivots; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::NumericDerivative:
        return builtins::evaluateNumericDerivative(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::NumericIntegral:
        return builtins::evaluateNumericIntegral(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Cbrt:
        return builtins::evaluateCbrt(arguments, registry_, mathematics_);
    case BuiltinId::Hypot:
        return builtins::evaluateHypot(arguments, registry_);
    case BuiltinId::Cis:
        return builtins::evaluateCis(arguments, registry_);
    case BuiltinId::Polar:
        return builtins::evaluatePolar(arguments, registry_);
    case BuiltinId::NextPow2:
        return builtins::evaluateNextPow2(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::DegreeToRadian:
        return builtins::evaluateDegreeToRadian(arguments, registry_, mathematics_);
    case BuiltinId::DegreeToGradian:
        return builtins::evaluateDegreeToGradian(arguments, registry_);
    case BuiltinId::RadianToDegree:
        return builtins::evaluateRadianToDegree(arguments, registry_, mathematics_);
    case BuiltinId::RadianToGradian:
        return builtins::evaluateRadianToGradian(arguments, registry_, mathematics_);
    case BuiltinId::GradianToDegree:
        return builtins::evaluateGradianToDegree(arguments, registry_);
    case BuiltinId::GradianToRadian:
        return builtins::evaluateGradianToRadian(arguments, registry_, mathematics_);
    case BuiltinId::Sum:
        return builtins::evaluateSum(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Product:
        return builtins::evaluateProduct(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Min:
        return builtins::evaluateMin(arguments, registry_);
    case BuiltinId::Max:
        return builtins::evaluateMax(arguments, registry_);
    case BuiltinId::Mean:
        return builtins::evaluateMean(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Median:
    case BuiltinId::Mode:
    case BuiltinId::Quantile:
    case BuiltinId::Percentile:
    case BuiltinId::VariancePopulation:
    case BuiltinId::VarianceSample:
    case BuiltinId::StddevPopulation:
    case BuiltinId::StddevSample:
    case BuiltinId::GeometricMean:
    case BuiltinId::HarmonicMean:
    case BuiltinId::Rms:
    case BuiltinId::MedianAbsoluteDeviation:
    case BuiltinId::MeanAbsoluteDeviation:
    case BuiltinId::Skewness:
    case BuiltinId::KurtosisPopulation:
    case BuiltinId::KurtosisSample:
    case BuiltinId::CoefficientVariation:
    case BuiltinId::StandardError:
    case BuiltinId::ZScore:
    case BuiltinId::Iqr:
    case BuiltinId::TrimMean:
    case BuiltinId::WinsorMean:
    case BuiltinId::Winsorized:
    case BuiltinId::Covariance:
    case BuiltinId::Correlation:
    case BuiltinId::SpearmanCorrelation:
    case BuiltinId::PercentRank:
        return builtins::evaluateStatistic(definition.id, arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Identity:
        return builtins::evaluateIdentity(arguments);
    case BuiltinId::Zeros:
        return builtins::evaluateZeros(arguments);
    case BuiltinId::MatrixGet:
        return builtins::evaluateMatrixGet(arguments);
    case BuiltinId::Trace:
        return builtins::evaluateTrace(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Rows:
        return builtins::evaluateRows(arguments);
    case BuiltinId::Cols:
        return builtins::evaluateCols(arguments);
    case BuiltinId::Diag:
        return builtins::evaluateDiag(arguments);
    case BuiltinId::VectorAdd:
        return builtins::evaluateVectorAdd(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorSubtract:
        return builtins::evaluateVectorSubtract(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorScale:
        return builtins::evaluateVectorScale(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorDot:
        return builtins::evaluateVectorDot(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorCross:
        return builtins::evaluateVectorCross(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorNorm:
        return builtins::evaluateVectorNorm(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorManhattan:
        return builtins::evaluateVectorManhattan(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorEuclidean:
        return builtins::evaluateVectorEuclidean(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorNormalize:
        return builtins::evaluateVectorNormalize(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorProject:
        return builtins::evaluateVectorProject(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorAngle:
        return builtins::evaluateVectorAngle(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorReflect:
        return builtins::evaluateVectorReflect(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorReflectAxis:
        return builtins::evaluateVectorReflectAxis(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorSum:
        return builtins::evaluateVectorSum(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Expm1:
    case BuiltinId::Log1p:
    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc:
    case BuiltinId::Sinhc:
    case BuiltinId::Tanhc:
    case BuiltinId::Expc:
        return builtins::evaluateStableElementary(
            definition.id, arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Log2:
    case BuiltinId::Log10:
    case BuiltinId::Gamma:
    case BuiltinId::LogGamma:
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
    case BuiltinId::Beta:
    case BuiltinId::BetaLog:
    case BuiltinId::GeneralizedBinomial:
    case BuiltinId::FallingFactorial:
    case BuiltinId::RisingFactorial:
        return builtins::evaluateSpecialFunction(
            definition.id, arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::RandSeed:
        return builtins::evaluateRandSeed(arguments, randomEngine_);
    case BuiltinId::Rand:
        return builtins::evaluateRand(arguments, randomEngine_);
    case BuiltinId::RandInt:
        return builtins::evaluateRandInt(arguments, randomEngine_);
    case BuiltinId::Choice:
        return builtins::evaluateChoice(arguments, randomEngine_);
    case BuiltinId::RandN:
        return builtins::evaluateRandN(
            arguments, randomEngine_, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Sqrt:
        return builtins::evaluateSqrt(arguments, registry_, mathematics_);
    case BuiltinId::Abs:
        return builtins::evaluateAbs(arguments, registry_, mathematics_);
    case BuiltinId::Sign:
        return builtins::evaluateSign(arguments, registry_, mathematics_);
    case BuiltinId::Re:
        return builtins::evaluateRe(arguments, registry_, mathematics_);
    case BuiltinId::Im:
        return builtins::evaluateIm(arguments, registry_, mathematics_);
    case BuiltinId::Conj:
        return builtins::evaluateConj(arguments, registry_, mathematics_);
    case BuiltinId::Sin:
        return builtins::evaluateSin(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Cos:
        return builtins::evaluateCos(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Tan:
        return builtins::evaluateTan(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Cot:
        return builtins::evaluateCot(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Sec:
        return builtins::evaluateSec(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Csc:
        return builtins::evaluateCsc(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Asin:
        return builtins::evaluateAsin(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Acos:
        return builtins::evaluateAcos(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Atan:
        return builtins::evaluateAtan(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Atan2:
        return builtins::evaluateAtan2(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Sinh:
        return builtins::evaluateSinh(arguments, registry_, mathematics_);
    case BuiltinId::Cosh:
        return builtins::evaluateCosh(arguments, registry_, mathematics_);
    case BuiltinId::Tanh:
        return builtins::evaluateTanh(arguments, registry_, mathematics_);
    case BuiltinId::Asinh:
        return builtins::evaluateAsinh(arguments, registry_, mathematics_);
    case BuiltinId::Acosh:
        return builtins::evaluateAcosh(arguments, registry_, mathematics_);
    case BuiltinId::Atanh:
        return builtins::evaluateAtanh(arguments, registry_, mathematics_);
    case BuiltinId::Csch:
        return builtins::evaluateCsch(arguments, registry_, mathematics_);
    case BuiltinId::Sech:
        return builtins::evaluateSech(arguments, registry_, mathematics_);
    case BuiltinId::Coth:
        return builtins::evaluateCoth(arguments, registry_, mathematics_);
    case BuiltinId::Arg:
        return builtins::evaluateArg(arguments, registry_, mathematics_);
    case BuiltinId::Log:
        return builtins::evaluateLog(arguments, registry_, mathematics_);
    case BuiltinId::Exp:
        return builtins::evaluateExp(arguments, registry_, mathematics_);
    case BuiltinId::NumericalApproximation:
        return evaluateNumericalApproximation(call, arguments);
    case BuiltinId::Precision: {
        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
        if (const auto result = builtins::evaluatePrecision(arguments, infinity->symbol))
            return *result;
        emitWarning("precision::unevaluated",
            "precision could not determine the guaranteed precision; the expression remains unevaluated");
        return expression::Expr::call(call.head, {arguments.front()});
    }
    case BuiltinId::Accuracy: {
        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
        if (const auto result = builtins::evaluateAccuracy(arguments, infinity->symbol))
            return *result;
        emitWarning("accuracy::unevaluated",
            "accuracy could not determine the guaranteed accuracy; the expression remains unevaluated");
        return expression::Expr::call(call.head, {arguments.front()});
    }
    case BuiltinId::Rationalize:
        if (const auto result = builtins::evaluateRationalize(arguments))
            return *result;
        emitWarning("rationalize::unevaluated",
            "rationalize could not convert part of the expression; it remains unevaluated");
        return expression::Expr::call(call.head, std::vector<expression::Expr>{arguments.begin(), arguments.end()});
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify: {
        if (arguments.empty() || arguments.size() > 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "simplify expects an expression and optional assumptions");
        mathematics::AssumptionSet assumptions;
        if (arguments.size() == 2)
            assumptions = mathematics::parseAssumptions(arguments[1], registry_, mathematics_);
        const simplification::SimplificationContext context{
            registry_, mathematics_, angleSemantics_, std::move(assumptions)};
        if (definition.id == BuiltinId::FullSimplify)
            return simplification::fullSimplify(arguments.front(), context);
        return simplification::Simplifier{}.simplify(arguments.front(), context);
    }
    case BuiltinId::Expand:
        return symbolic::expandExpression(
            arguments.front(), registry_, mathematics_, angleSemantics_);
    case BuiltinId::Factor:
        return symbolic::factorExpression(
            arguments.front(), registry_, mathematics_, angleSemantics_);
    case BuiltinId::Collect: {
        if (arguments.size() != 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "collect expects an expression and a symbol or symbol array");
        std::vector<expression::Symbol> variables;
        if (arguments[1].isSymbol())
            variables.push_back(arguments[1].asSymbol());
        else if (arguments[1].isArray() && arguments[1].asArray().rank() == 1) {
            for (const expression::Expr& item : arguments[1].asArray().elements) {
                if (!item.isSymbol())
                    error::throwCalcError(
                        error::CalcErrorType::Type,
                        "collect variable array must contain only symbols");
                variables.push_back(item.asSymbol());
            }
        }
        else
            error::throwCalcError(
                error::CalcErrorType::Type,
                "collect expects an expression and a symbol or symbol array");
        return symbolic::collectExpression(
            arguments[0], variables, registry_, mathematics_, angleSemantics_);
    }
    case BuiltinId::Solve: {
        if (arguments.size() < 2 || arguments.size() > 3)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "solve expects equation(s), variable(s), and optional constraints");

        std::vector<expression::Symbol> variables;
        if (arguments[1].isSymbol())
            variables.push_back(arguments[1].asSymbol());
        else if (arguments[1].isArray() && arguments[1].asArray().rank() == 1) {
            for (const expression::Expr& item : arguments[1].asArray().elements) {
                if (!item.isSymbol())
                    error::throwCalcError(
                        error::CalcErrorType::Type,
                        "solve variable array must contain only symbols");
                if (std::find(variables.begin(), variables.end(), item.asSymbol()) != variables.end())
                    error::throwCalcError(
                        error::CalcErrorType::Type,
                        "solve variable array contains a duplicate symbol");
                variables.push_back(item.asSymbol());
            }
        }
        else
            error::throwCalcError(
                error::CalcErrorType::Type,
                "solve expects a symbol or symbol array as its second argument");

        solver::SolveConstraints constraints;
        if (arguments.size() == 3)
            constraints = solver::parseSolveConstraints(
                arguments[2], variables, registry_, mathematics_, angleSemantics_);

        solver::SolutionSet solutions = [&]() {
            if (variables.size() == 1 && !arguments[0].isArray()) {
                const bool realDomain = constraints.domain
                    && mathematics::isSubdomainOf(
                        *constraints.domain, mathematics::NumericDomain::Real);
                if (realDomain) {
                    if (auto transcendental = solver::solveRealInjectiveFunctionRelation(
                            arguments[0], variables.front(), registry_, mathematics_,
                            angleSemantics_, constraints.assumptions))
                        return *transcendental;
                }
                return solver::solveUnivariatePolynomialRelation(
                    arguments[0], variables.front(), registry_, mathematics_, angleSemantics_);
            }

            std::vector<expression::Expr> equations;
            if (arguments[0].isArray() && arguments[0].asArray().rank() == 1)
                equations.assign(
                    arguments[0].asArray().elements.begin(),
                    arguments[0].asArray().elements.end());
            else
                equations.push_back(arguments[0]);

            // 一変数のrelation配列は論理積として扱う。等式があれば先に解いて有限候補を作り、残りをexact constraintとして絞る。
            // 等式がなければ最初の不等式からReal領域branchを作り、残りの不等式を条件として交差させる。
            if (variables.size() == 1 && !equations.empty()) {
                auto first = equations.begin();
                const auto equality = std::find_if(
                    equations.begin(), equations.end(), [&](const expression::Expr& item) {
                        return item.isCall()
                            && item.asCall().head.sameIdentity(registry_.symbol(BuiltinId::Equal));
                    });
                if (equality != equations.end())
                    first = equality;

                solver::SolutionSet result = solver::solveUnivariatePolynomialRelation(
                    *first, variables.front(), registry_, mathematics_, angleSemantics_);
                for (auto iterator = equations.begin(); iterator != equations.end(); ++iterator) {
                    if (iterator == first)
                        continue;
                    const std::array<expression::Symbol, 1> oneVariable{variables.front()};
                    const solver::SolveConstraints relationConstraint =
                        solver::parseSolveConstraints(
                            *iterator, oneVariable, registry_, mathematics_, angleSemantics_);
                    result = solver::applySolveConstraints(
                        std::move(result), relationConstraint,
                        registry_, mathematics_, angleSemantics_);
                }
                return result;
            }

            return solver::solveLinearPolynomialSystem(
                equations, variables, registry_, mathematics_, angleSemantics_);
        }();

        if (constraints.domain
            && *constraints.domain == mathematics::NumericDomain::Complex
            && !solutions.variables().empty()
            && solutions.variables().front().domain == mathematics::NumericDomain::Real) {
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "Ordered inequalities are defined only over Real or a subdomain");
        }

        solutions = solver::applySolveConstraints(
            std::move(solutions), constraints, registry_, mathematics_, angleSemantics_);
        if (containsUnresolvedSolution(solutions))
            emitWarning("solve::unresolved",
                "solve could not determine a complete solution set; unresolved cases remain");
        return expression::Expr::solutionSet(std::move(solutions));
    }
    case BuiltinId::Set:
        return evaluateSet(arguments);
    case BuiltinId::SetDelayed:
        return evaluateSetDelayed(call, arguments);
    case BuiltinId::Less:
    case BuiltinId::LessEqual:
    case BuiltinId::Greater:
    case BuiltinId::GreaterEqual:
    case BuiltinId::Equal:
    case BuiltinId::NotEqual:
        return builtins::evaluateComparison(registry_.symbol(definition.id), arguments);
    case BuiltinId::LogicalAnd:
        return builtins::evaluateLogicalAnd(arguments, registry_);
    case BuiltinId::Element:
        return builtins::evaluateElement(arguments, registry_, mathematics_);
    case BuiltinId::If:
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "If must be handled by the evaluation machine");
    case BuiltinId::History:
        return evaluateHistory(arguments);
    case BuiltinId::InputHistory:
        return evaluateAbsoluteHistory(arguments, true);
    case BuiltinId::OutputHistory:
        return evaluateAbsoluteHistory(arguments, false);
    case BuiltinId::Exit:
        if (context_ && context_->exitRequested)
            *context_->exitRequested = true;
        return expression::Expr{true};
    case BuiltinId::Clear:
        if (context_ && context_->clearRequested)
            *context_->clearRequested = true;
        return expression::Expr{true};
    case BuiltinId::Definitions:
        return evaluateDefinitions();
    case BuiltinId::Undefine:
        return evaluateUndefine(arguments);
    case BuiltinId::AngleMode: {
        const auto angleSymbol = [&](mathematics::AngleUnit unit) -> expression::Expr {
            std::string_view name;
            switch (unit) {
            case mathematics::AngleUnit::Degree: name = "Deg"; break;
            case mathematics::AngleUnit::Radian: name = "Rad"; break;
            case mathematics::AngleUnit::Gradian: name = "Grad"; break;
            }

            const auto* predefined = symbolRegistry_.find(name);
            if (!predefined)
                error::throwCalcError(
                    error::CalcErrorType::Internal,
                    "Angle-mode symbol is not registered");
            return expression::Expr{predefined->symbol};
        };

        if (arguments.empty())
            return angleSymbol(angleSemantics_.defaultUnit());

        if (arguments.size() != 1 || !arguments.front().isSymbol())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "angleMode expects no argument or one of Rad, Deg, Grad");

        const auto* predefined = symbolRegistry_.find(arguments.front().asSymbol());
        if (!predefined)
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "angleMode expects one of Rad, Deg, Grad");

        mathematics::AngleUnit unit;
        switch (predefined->id) {
        case symbols::PredefinedSymbolId::DegreeUnit:
            unit = mathematics::AngleUnit::Degree;
            break;
        case symbols::PredefinedSymbolId::RadianUnit:
            unit = mathematics::AngleUnit::Radian;
            break;
        case symbols::PredefinedSymbolId::GradianUnit:
            unit = mathematics::AngleUnit::Gradian;
            break;
        default:
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "angleMode expects one of Rad, Deg, Grad");
        }

        if (!context_ || !context_->angleSemantics)
            error::throwCalcError(
                error::CalcErrorType::Internal,
                "angleMode requires a mutable kernel session");
        context_->angleSemantics->setDefaultUnit(unit);
        return angleSymbol(unit);
    }
    case BuiltinId::UnitApplied: {
        if (arguments.size() != 2 || !arguments[1].isString())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "UnitApplied requires a value and unit name");

        // 角度単位だけは数学層で意味を持つため、綴りを正規化する。
        // 長さ等の単位はまだ演算しないが、将来のunit systemへ渡せるようUnitApplied式として保持し、評価エラーにはしない。
        std::string unit = arguments[1].asString();
        if (const auto angleUnit = mathematics::AngleSemantics::parseUnit(unit))
            unit = std::string{mathematics::AngleSemantics::canonicalName(*angleUnit)};
        return expression::Expr::call(
            registry_.symbol(BuiltinId::UnitApplied),
            {arguments[0], expression::Expr{std::move(unit)}});
    }
    }

    error::throwCalcError(
        error::CalcErrorType::Internal,
        "Builtin dispatch is incomplete");
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

expression::Expr Evaluator::evaluateAbsoluteHistory(
    std::span<const expression::Expr> arguments,
    bool input) {
    const auto number = positiveSizeValue(arguments.front());
    if (!number)
        error::throwCalcError(
            error::CalcErrorType::Type,
            input ? "In index must be a positive integer" : "Out index must be a positive integer");

    if (!context_)
        error::throwCalcError(error::CalcErrorType::Evaluation, "Session history is not available");

    const auto& entries = input ? context_->inputs : context_->outputs;
    if (input && *number == context_->inputs.size()
        && *number <= context_->outputs.size() && !context_->outputs[*number - 1])
        error::throwCalcError(
            error::CalcErrorType::Evaluation,
            "In cannot reference the input currently being evaluated");
    if (*number > entries.size() || !entries[*number - 1])
        error::throwCalcError(
            error::CalcErrorType::Evaluation,
            std::string{input ? "Input" : "Output"} + " entry is not available at index "
                + std::to_string(*number));
    return *entries[*number - 1];
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

expression::Expr Evaluator::evaluateNumericalApproximation(
    const expression::CallExpr& call,
    std::span<const expression::Expr> arguments) {
    constexpr std::size_t defaultFractionalDigits = 16;
    std::size_t fractionalDigits = defaultFractionalDigits;

    if (arguments.size() == 2) {
        const auto requestedDigits = positiveSizeValue(arguments[1]);
        if (!requestedDigits) {
            if (const auto origin = originOf(call.arguments[1]))
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

    const expression::Expr& value = arguments.front();

    // Nはscalarだけでなく配列へも要素単位に作用する。
    // FFT/行列等のexact配列を表示用近似へ落とす際に、各builtinが独自のdigits引数を持つ必要をなくす。
    std::function<expression::Expr(const expression::Expr&)> approximate;
    approximate = [&](const expression::Expr& current) -> expression::Expr {
        if (current.isArray()) {
            const auto& array = current.asArray();
            std::vector<expression::Expr> elements;
            elements.reserve(array.elements.size());
            for (const expression::Expr& element : array.elements)
                elements.push_back(approximate(element));
            return expression::Expr::array(array.shape, std::move(elements));
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

std::optional<source::SourceReference> Evaluator::originOf(
    const expression::Expr& expression) const {
    if (!origins_)
        return std::nullopt;

    return origins_->find(expression);
}

} // namespace mmcal::evaluation
