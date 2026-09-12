// 極限limit
#include "symbolic/cases.hpp"
#include "limit.hpp"

#include "builtins/array_vector.hpp"
#include "builtins/polynomial_ideal.hpp"
#include "error/error_message.hpp"
#include "mathematics/assumption_parser.hpp"
#include "mathematics/definedness.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "evaluation/iterator_spec.hpp"
#include "expression/array_utils.hpp"
#include "expression/exact_value.hpp"
#include "mathematics/knowledge_context.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "simplification/full_simplifier.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "symbolic/algebra_transforms.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/integration.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/series.hpp"
#include "symbolic/substitution.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::symbolic {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using expression::exact::integer;
using expression::exact::isZero;
using expression::exact::rational;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

constexpr std::size_t maximumLimitDepth = 24;
constexpr std::size_t maximumLHopitalSteps = 12;

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return builtins.isCallTo(expression, id);
}

[[nodiscard]] Expr call(
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id,
    std::vector<Expr> arguments) {
    return Expr::call(builtins.symbol(id), std::move(arguments));
}

[[nodiscard]] bool sameVariable(
    const Expr& expression,
    const expression::Symbol& variable) {
    return expression.isSymbol() && expression.asSymbol().sameIdentity(variable);
}

[[nodiscard]] bool iteratorBinds(
    const Expr& expression,
    const expression::Symbol& variable) {
    const auto iterator = evaluation::parseRangeIteratorSpec(expression);
    return iterator && iterator->variable.sameIdentity(variable);
}

[[nodiscard]] bool tableIteratorBinds(
    const Expr& expression,
    const expression::Symbol& variable) {
    const auto iterator = evaluation::parseTableIteratorSpec(expression);
    return iterator && iterator->variable.sameIdentity(variable);
}

[[nodiscard]] bool derivativeSpecBinds(
    const Expr& expression,
    const expression::Symbol& variable) {
    if (sameVariable(expression, variable))
        return true;
    if (!expression.isArray() || expression.asArray().rank() != 1
        || expression.asArray().size() < 1)
        return false;
    return sameVariable(expression.asArray().element(0), variable);
}

[[nodiscard]] bool solveVariableSpecBinds(
    const Expr& expression,
    const expression::Symbol& variable) {
    if (sameVariable(expression, variable))
        return true;
    if (!expression.isArray() || expression.asArray().rank() != 1)
        return false;
    for (std::size_t i = 0; i < expression.asArray().size(); ++i)
        if (sameVariable(expression.asArray().element(i), variable))
            return true;
    return false;
}

[[nodiscard]] std::optional<Expr> materializeHeldDerivative(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Derivative))
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() < 2)
        return std::nullopt;

    Expr result = arguments[0];
    for (std::size_t i = 1; i < arguments.size(); ++i) {
        expression::Symbol variable;
        std::uint64_t order = 1;
        if (arguments[i].isSymbol()) {
            variable = arguments[i].asSymbol();
        }
        else if (arguments[i].isArray()) {
            const auto& spec = arguments[i].asArray();
            if (spec.rank() != 1 || spec.size() != 2 || !spec.element(0).isSymbol()
                || !spec.element(1).isNumber() || !spec.element(1).asNumber().isReal()
                || !spec.element(1).asNumber().asReal().isInteger())
                return std::nullopt;
            const auto parsed = numeric::tryToUint64(
                spec.element(1).asNumber().asReal().asInteger());
            if (!parsed)
                return std::nullopt;
            variable = spec.element(0).asSymbol();
            order = *parsed;
        }
        else {
            return std::nullopt;
        }

        if (order > 1) {
            if (auto known = differentiateKnownRepeatedExpression(
                    result, variable, order, builtins, mathematics, angles)) {
                result = std::move(*known);
                continue;
            }
        }
        for (std::uint64_t derivative = 0; derivative < order; ++derivative) {
            // D本体と同様，固定orderではなくrequest-scoped budgetで実作業量を制御する。
            // exact 0へ到達した後は残りの高階微分を反復する必要がない。
            evaluation::consumeEvaluationBudget(evaluation::EvaluationResource::EvaluationStep);
            result = differentiateExpression(result, variable, builtins, mathematics, angles);
            if (order > 1)
                result = canonicalizeDerivativeOutput(result, builtins, mathematics, angles);
            if (result.isNumber() && result.asNumber().isZero())
                break;
        }
    }
    return result;
}

[[nodiscard]] std::optional<Expr> materializeHeldVectorCalculus(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;

    switch (definition->id) {
    case BuiltinId::Gradient:
        return arguments.size() == 2
            ? std::optional<Expr>{mmcal::builtins::evaluateGradient(arguments, builtins, mathematics, angles)}
            : std::nullopt;
    case BuiltinId::Divergence:
        return arguments.size() == 2
            ? std::optional<Expr>{mmcal::builtins::evaluateDivergence(arguments, builtins, mathematics, angles)}
            : std::nullopt;
    case BuiltinId::Curl:
        return arguments.size() == 2
            ? std::optional<Expr>{mmcal::builtins::evaluateCurl(arguments, builtins, mathematics, angles)}
            : std::nullopt;
    case BuiltinId::Laplacian:
        return arguments.size() == 2
            ? std::optional<Expr>{mmcal::builtins::evaluateLaplacian(arguments, builtins, mathematics, angles)}
            : std::nullopt;
    case BuiltinId::Jacobian:
        return arguments.size() == 2
            ? std::optional<Expr>{mmcal::builtins::evaluateJacobian(arguments, builtins, mathematics, angles)}
            : std::nullopt;
    case BuiltinId::Hessian:
        return arguments.size() == 2
            ? std::optional<Expr>{mmcal::builtins::evaluateHessian(arguments, builtins, mathematics, angles)}
            : std::nullopt;
    case BuiltinId::DirectionalDerivative:
        return arguments.size() == 3
            ? std::optional<Expr>{mmcal::builtins::evaluateDirectionalDerivative(arguments, builtins, mathematics, angles)}
            : std::nullopt;
    default:
        return std::nullopt;
    }
}

// Limitの直接点代入専用のfree-symbol判定。
// evaluatorのHold/binderを無視した通常containsSymbolを使うと，内側LimitやTableの
// control variableまで外側のxとして数えてしまい，後段でbinderを破壊する。
[[nodiscard]] bool containsFreeLimitSymbol(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (sameVariable(expression, variable))
        return true;
    if (expression.isCall()) {
        const auto* definition = builtins.find(expression.asCall().head);
        const auto& arguments = expression.asCall().arguments;
        if (definition) {
            if (definition->id == BuiltinId::Limit && arguments.size() >= 2
                && sameVariable(arguments[1], variable)) {
                for (std::size_t i = 2; i < arguments.size(); ++i)
                    if (containsFreeLimitSymbol(arguments[i], variable, builtins))
                        return true;
                return false;
            }
            if ((definition->id == BuiltinId::SymbolicIntegral
                    || definition->id == BuiltinId::NumericIntegral)
                && arguments.size() >= 2 && iteratorBinds(arguments[1], variable)) {
                const auto iterator = evaluation::parseRangeIteratorSpec(arguments[1]);
                return containsFreeLimitSymbol(iterator->lower, variable, builtins)
                    || containsFreeLimitSymbol(iterator->upper, variable, builtins);
            }
            if (definition->id == BuiltinId::Table && arguments.size() == 2
                && tableIteratorBinds(arguments[1], variable)) {
                const auto iterator = evaluation::parseTableIteratorSpec(arguments[1]);
                for (const Expr& bound : iterator->rangeArguments)
                    if (containsFreeLimitSymbol(bound, variable, builtins))
                        return true;
                return false;
            }
            if (definition->id == BuiltinId::Solve && arguments.size() >= 2
                && solveVariableSpecBinds(arguments[1], variable))
                return false;
            if (definition->id == BuiltinId::NumericDerivative && arguments.size() >= 2
                && sameVariable(arguments[1], variable)) {
                for (std::size_t i = 2; i < arguments.size(); ++i)
                    if (containsFreeLimitSymbol(arguments[i], variable, builtins))
                        return true;
                return false;
            }
        }
        for (const Expr& argument : arguments)
            if (containsFreeLimitSymbol(argument, variable, builtins))
                return true;
        return false;
    }
    if (expression.isArray()) {
        for (std::size_t i = 0; i < expression.asArray().size(); ++i)
            if (containsFreeLimitSymbol(expression.asArray().element(i), variable, builtins))
                return true;
        return false;
    }
    if (expression.isList()) {
        for (const Expr& element : expression.asList().elements)
            if (containsFreeLimitSymbol(element, variable, builtins))
                return true;
    }
    return false;
}

[[nodiscard]] std::optional<Expr> substituteLimitFreeSymbol(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins) {
    if (sameVariable(expression, variable))
        return value;

    if (expression.isCall()) {
        const auto* definition = builtins.find(expression.asCall().head);
        const auto& source = expression.asCall().arguments;
        if (definition) {
            // D・不定積分・Seriesはcontrol variableと結果の自由変数が同一である。
            // 未評価のまま点代入すると D[f[0],0] 等になるため，ここでは変形せず
            // Limit本体を未解決として残す。
            if (definition->id == BuiltinId::Derivative && source.size() >= 2
                && derivativeSpecBinds(source[1], variable))
                return std::nullopt;
            if (definition->id == BuiltinId::SymbolicIntegral && source.size() >= 2
                && sameVariable(source[1], variable))
                return std::nullopt;
            if (definition->id == BuiltinId::Series && source.size() >= 2
                && iteratorBinds(source[1], variable))
                return std::nullopt;

            std::vector<Expr> arguments = source;
            if (definition->id == BuiltinId::Limit && source.size() >= 2
                && sameVariable(source[1], variable)) {
                for (std::size_t i = 2; i < source.size(); ++i) {
                    auto replaced = substituteLimitFreeSymbol(source[i], variable, value, builtins);
                    if (!replaced)
                        return std::nullopt;
                    arguments[i] = std::move(*replaced);
                }
                return Expr::rebuildCall(expression.asCall(), std::move(arguments));
            }
            if ((definition->id == BuiltinId::SymbolicIntegral
                    || definition->id == BuiltinId::NumericIntegral)
                && source.size() >= 2 && iteratorBinds(source[1], variable)) {
                const auto iterator = evaluation::parseRangeIteratorSpec(source[1]);
                auto lower = substituteLimitFreeSymbol(iterator->lower, variable, value, builtins);
                auto upper = substituteLimitFreeSymbol(iterator->upper, variable, value, builtins);
                if (!lower || !upper)
                    return std::nullopt;
                arguments[1] = Expr::array({3}, {
                    Expr{iterator->variable}, std::move(*lower), std::move(*upper)});
                return Expr::rebuildCall(expression.asCall(), std::move(arguments));
            }
            if (definition->id == BuiltinId::Table && source.size() == 2
                && tableIteratorBinds(source[1], variable)) {
                const auto iterator = evaluation::parseTableIteratorSpec(source[1]);
                std::vector<Expr> spec;
                spec.reserve(iterator->rangeArguments.size() + 1);
                spec.push_back(Expr{iterator->variable});
                for (const Expr& bound : iterator->rangeArguments) {
                    auto replaced = substituteLimitFreeSymbol(bound, variable, value, builtins);
                    if (!replaced)
                        return std::nullopt;
                    spec.push_back(std::move(*replaced));
                }
                const std::size_t specSize = spec.size();
                arguments[1] = Expr::array({specSize}, std::move(spec));
                return Expr::rebuildCall(expression.asCall(), std::move(arguments));
            }
            if (definition->id == BuiltinId::Solve && source.size() >= 2
                && solveVariableSpecBinds(source[1], variable)) {
                // solveの未知変数はrelationとconstraintsの双方を束縛する。外側Limitの
                // 同名変数を代入すると solve[0==1,0] のようにbinderを破壊する。
                return expression;
            }
            if (definition->id == BuiltinId::NumericDerivative && source.size() >= 2
                && sameVariable(source[1], variable)) {
                // ndの変数は被微分式だけを束縛する。評価点・digits側に同名の自由変数が
                // あればそこだけ外側Limitの点代入を許し，control slotは保護する。
                for (std::size_t i = 2; i < source.size(); ++i) {
                    auto replaced = substituteLimitFreeSymbol(source[i], variable, value, builtins);
                    if (!replaced)
                        return std::nullopt;
                    arguments[i] = std::move(*replaced);
                }
                return Expr::rebuildCall(expression.asCall(), std::move(arguments));
            }
        }

        std::vector<Expr> arguments;
        arguments.reserve(source.size());
        for (const Expr& argument : source) {
            auto replaced = substituteLimitFreeSymbol(argument, variable, value, builtins);
            if (!replaced)
                return std::nullopt;
            arguments.push_back(std::move(*replaced));
        }
        return Expr::rebuildCall(expression.asCall(), std::move(arguments));
    }

    if (expression.isArray()) {
        const auto& array = expression.asArray();
        std::vector<Expr> elements;
        elements.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i) {
            auto replaced = substituteLimitFreeSymbol(array.element(i), variable, value, builtins);
            if (!replaced)
                return std::nullopt;
            elements.push_back(std::move(*replaced));
        }
        return Expr::array(array.shape, std::move(elements));
    }
    if (expression.isList()) {
        std::vector<Expr> elements;
        elements.reserve(expression.asList().elements.size());
        for (const Expr& element : expression.asList().elements) {
            auto replaced = substituteLimitFreeSymbol(element, variable, value, builtins);
            if (!replaced)
                return std::nullopt;
            elements.push_back(std::move(*replaced));
        }
        return expression::braceValue(std::move(elements));
    }
    return expression;
}


[[nodiscard]] Expr simplify(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
}

[[nodiscard]] Expr fullSimplify(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplification::fullSimplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions},
        simplification::FullSimplificationOptions{48});
}

[[nodiscard]] std::optional<std::vector<expression::Symbol>> transformVariables(
    const Expr& expression) {
    if (expression.isSymbol())
        return std::vector<expression::Symbol>{expression.asSymbol()};
    if (!expression.isArray() || expression.asArray().rank() != 1
        || expression.asArray().size() == 0)
        return std::nullopt;
    std::vector<expression::Symbol> variables;
    variables.reserve(expression.asArray().size());
    for (std::size_t i = 0; i < expression.asArray().size(); ++i) {
        const Expr item = expression.asArray().element(i);
        if (!item.isSymbol())
            return std::nullopt;
        variables.push_back(item.asSymbol());
    }
    return variables;
}

[[nodiscard]] std::optional<Expr> materializeHeldAlgebraTransform(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;

    switch (definition->id) {
    case BuiltinId::Expand:
        return arguments.size() == 1
            ? std::optional<Expr>{expandExpression(arguments[0], builtins, mathematics, angles)}
            : std::nullopt;
    case BuiltinId::Factor:
        return arguments.size() == 1
            ? std::optional<Expr>{factorExpression(arguments[0], builtins, mathematics, angles)}
            : std::nullopt;
    case BuiltinId::Collect: {
        if (arguments.size() != 2)
            return std::nullopt;
        const auto variables = transformVariables(arguments[1]);
        if (!variables)
            return std::nullopt;
        return collectExpression(arguments[0], *variables, builtins, mathematics, angles);
    }
    case BuiltinId::GroebnerBasis:
        return arguments.size() >= 2 && arguments.size() <= 3
            ? std::optional<Expr>{mmcal::builtins::evaluateGroebnerBasis(arguments, builtins)}
            : std::nullopt;
    case BuiltinId::PolynomialReduce:
        return arguments.size() >= 3 && arguments.size() <= 4
            ? std::optional<Expr>{mmcal::builtins::evaluatePolynomialReduce(arguments, builtins)}
            : std::nullopt;
    default:
        return std::nullopt;
    }
}

[[nodiscard]] std::optional<Expr> materializeHeldSeriesCall(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isHead(expression, builtins, BuiltinId::Series))
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() < 2 || arguments.size() > 3 || !arguments[1].isArray())
        return std::nullopt;
    const auto& spec = arguments[1].asArray();
    if (spec.rank() != 1 || spec.size() != 3 || !spec.element(0).isSymbol())
        return std::nullopt;
    const Expr orderExpression = spec.element(2);
    if (!orderExpression.isNumber() || !orderExpression.asNumber().isReal()
        || !orderExpression.asNumber().asReal().isInteger())
        return std::nullopt;
    const auto order = numeric::tryToUint64(orderExpression.asNumber().asReal().asInteger());
    if (!order || *order > 1024)
        return std::nullopt;

    mathematics::AssumptionSet seriesAssumptions;
    if (arguments.size() == 3)
        seriesAssumptions = mathematics::parseAssumptions(arguments[2], builtins, mathematics);
    return seriesExpression(
        arguments[0], spec.element(0).asSymbol(), spec.element(1),
        static_cast<std::size_t>(*order), builtins, mathematics, angles, seriesAssumptions);
}

[[nodiscard]] std::optional<Expr> materializeHeldNormal(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || (definition->id != BuiltinId::Normal
        && definition->id != BuiltinId::ToNormal)
        || expression.asCall().arguments.size() != 1)
        return std::nullopt;

    Expr argument = expression.asCall().arguments.front();
    if (auto series = materializeHeldSeriesCall(argument, builtins, mathematics, angles))
        argument = std::move(*series);

    if (definition->id == BuiltinId::ToNormal)
        return toNormalExpression(argument, builtins, mathematics, angles);
    if (const auto series = parseSeriesData(argument, builtins))
        return normalSeriesExpression(*series, builtins, mathematics, angles);
    return argument;
}

[[nodiscard]] Expr negate(
    Expr value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplify(call(builtins, BuiltinId::Negate, {std::move(value)}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr divide(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplify(call(builtins, BuiltinId::Divide, {std::move(lhs), std::move(rhs)}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr unresolved(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& point,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Expr> arguments{expression, Expr{variable}, point};
    if (direction == LimitDirection::Left)
        arguments.push_back(integer(-1));
    else if (direction == LimitDirection::Right)
        arguments.push_back(integer(1));
    return call(builtins, BuiltinId::Limit, std::move(arguments));
}

[[nodiscard]] bool isInfinity(const Expr& expression, const expression::Symbol& infinity) {
    return expression.isSymbol() && expression.asSymbol().sameIdentity(infinity);
}

[[nodiscard]] bool isUnaryFunctionOfVariable(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return isHead(expression, builtins, id)
        && expression.asCall().arguments.size() == 1
        && expression.asCall().arguments[0].isSymbol()
        && expression.asCall().arguments[0].asSymbol().sameIdentity(variable);
}

[[nodiscard]] bool isNegativeInfinity(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    return isHead(expression, builtins, BuiltinId::Negate)
        && expression.asCall().arguments.size() == 1
        && isInfinity(expression.asCall().arguments[0], infinity);
}

[[nodiscard]] Expr positiveInfinity(const expression::Symbol& infinity) {
    return Expr{infinity};
}

[[nodiscard]] Expr signedInfinity(
    int sign,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    if (sign >= 0)
        return positiveInfinity(infinity);
    return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<Expr> finiteLimitFromLocalSeries(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& point,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    constexpr std::size_t kSeriesOrder = 6;
    auto expanded = seriesExpression(
        expression, variable, point, kSeriesOrder,
        builtins, mathematics, angles, assumptions);
    if (!expanded) return std::nullopt;
    auto data = parseSeriesData(*expanded, builtins);
    if (!data || data->exponentDenominator == 0)
        return std::nullopt;

    for (std::size_t i = 0; i < data->coefficients.size(); ++i) {
        const std::int64_t exponentNumerator =
            data->minimumExponent + static_cast<std::int64_t>(i);
        const bool ordinaryNonZero = !isZero(data->coefficients[i]);
        bool logarithmicNonZero = false;
        for (const auto& layer : data->logarithmicCoefficients) {
            if (!isZero(layer[i])) {
                logarithmicNonZero = true;
                break;
            }
        }
        if (!ordinaryNonZero && !logarithmicNonZero)
            continue;
        if (exponentNumerator < 0)
            return std::nullopt;
        if (exponentNumerator == 0 && logarithmicNonZero)
            return std::nullopt;
        if (exponentNumerator == 0 && ordinaryNonZero)
            return fullSimplify(
                data->coefficients[i], builtins, mathematics, angles, assumptions);
        if (exponentNumerator > 0)
            return integer(0);
    }
    return integer(0);
}

[[nodiscard]] bool domainConditionsHold(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const auto conditions = mathematics::expressionDomainConditions(expression, builtins, mathematics);
    if (!conditions)
        return true;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    for (const auto& predicate : conditions->predicates()) {
        if (knowledge.prove(predicate) == mathematics::TruthValue::False)
            return false;
    }
    return true;
}

[[nodiscard]] bool recoverableSpeculativeError(error::CalcErrorType type) {
    // Limit内部の点代入・Seriesは候補生成であり，特異点に当たっただけなら
    // 公開評価の失敗にせず別の証明経路へ回す。資源制限と内部不変条件は伝播する。
    switch (type) {
    case error::CalcErrorType::Domain:
    case error::CalcErrorType::Type:
    case error::CalcErrorType::Overflow:
    case error::CalcErrorType::Evaluation:
        return true;
    case error::CalcErrorType::Syntax:
    case error::CalcErrorType::Name:
    case error::CalcErrorType::Internal:
    case error::CalcErrorType::ResourceLimit:
        return false;
    }
    return false;
}

[[nodiscard]] std::optional<Expr> substitutedLimitCandidate(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& point,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    try {
        auto substituted = substituteLimitFreeSymbol(
            expression, variable, point, builtins);
        if (!substituted
            || !domainConditionsHold(*substituted, builtins, mathematics, assumptions))
            return std::nullopt;
        return simplify(
            std::move(*substituted), builtins, mathematics, angles, assumptions);
    }
    catch (const error::CalcError& exception) {
        if (!recoverableSpeculativeError(exception.type()))
            throw;
        return std::nullopt;
    }
}

struct LocalPolynomialBehavior final {
    std::size_t zeroOrder = 0;
    Rational leading{};
};

[[nodiscard]] std::optional<LocalPolynomialBehavior> localBehavior(
    RationalPolynomial polynomial,
    const Rational& point) {
    if (polynomial.isZero())
        return std::nullopt;

    std::size_t order = 0;
    while (evaluatePolynomial(polynomial, point).isZero()) {
        const auto quotient = divideByLinearFactor(polynomial, point);
        if (!quotient)
            return std::nullopt;
        polynomial = *quotient;
        ++order;
    }
    return LocalPolynomialBehavior{order, evaluatePolynomial(polynomial, point)};
}

[[nodiscard]] int rationalSign(const Rational& value) {
    if (value.isZero()) return 0;
    return value.numerator().isNegative() ? -1 : 1;
}

[[nodiscard]] std::optional<Expr> rationalFunctionFiniteLimit(
    const Expr& expression,
    const expression::Symbol& variable,
    const Rational& point,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    Expr numerator = expression;
    Expr denominator = integer(1);
    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        numerator = expression.asCall().arguments[0];
        denominator = expression.asCall().arguments[1];
    }

    const auto np = toRationalPolynomial(numerator, variable, builtins, {256, 2048});
    const auto dp = toRationalPolynomial(denominator, variable, builtins, {256, 2048});
    if (!np || !dp || dp->isZero())
        return std::nullopt;

    const auto nb = localBehavior(*np, point);
    const auto db = localBehavior(*dp, point);
    if (!nb || !db)
        return std::nullopt;

    if (nb->zeroOrder > db->zeroOrder)
        return integer(0);

    if (nb->zeroOrder == db->zeroOrder)
        return rational(nb->leading / db->leading);

    const std::size_t poleOrder = db->zeroOrder - nb->zeroOrder;
    int sign = rationalSign(nb->leading / db->leading);
    if (direction == LimitDirection::TwoSided && (poleOrder % 2) != 0)
        return std::nullopt;
    if (direction == LimitDirection::Left && (poleOrder % 2) != 0)
        sign = -sign;
    return signedInfinity(sign, builtins, mathematics, angles, infinity, assumptions);
}

[[nodiscard]] bool isRationalExpressionCandidate(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (sameVariable(expression, variable))
        return true;
    if (expression::exact::realRational(expression))
        return true;
    if (!expression.isCall())
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return false;
    const auto& arguments = expression.asCall().arguments;
    switch (definition->id) {
    case BuiltinId::Negate:
        return arguments.size() == 1
            && isRationalExpressionCandidate(arguments[0], variable, builtins);
    case BuiltinId::Add:
    case BuiltinId::Subtract:
    case BuiltinId::Multiply:
    case BuiltinId::Divide:
        if (arguments.empty()) return false;
        for (const Expr& argument : arguments)
            if (!isRationalExpressionCandidate(argument, variable, builtins))
                return false;
        return true;
    case BuiltinId::Power:
        if (arguments.size() != 2
            || !isRationalExpressionCandidate(arguments[0], variable, builtins))
            return false;
        if (const auto exponent = expression::exact::realRational(arguments[1]))
            return exponent->isInteger();
        return false;
    default:
        return false;
    }
}

[[nodiscard]] std::optional<Expr> rationalFunctionInfiniteLimit(
    const Expr& expression,
    const expression::Symbol& variable,
    bool negativeInfinityPoint,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    // 1 + 1/x のようにAST上で単一Divideになっていない有理式も，
    // 既存のexact rational normalizerで通分してから最高次項を比較する。
    // limit専用に別の通分算法を持たず，積分器と同じ正規形を共有する。
    Expr normalized = expression;
    if (isRationalExpressionCandidate(expression, variable, builtins)) {
        if (auto rationalized = normalizeRationalExpression(
                expression, variable, builtins, mathematics, angles))
            normalized = std::move(*rationalized);
    }

    Expr numerator = normalized;
    Expr denominator = integer(1);
    if (isHead(normalized, builtins, BuiltinId::Divide)
        && normalized.asCall().arguments.size() == 2) {
        numerator = normalized.asCall().arguments[0];
        denominator = normalized.asCall().arguments[1];
    }

    const auto np = toRationalPolynomial(numerator, variable, builtins, {256, 2048});
    const auto dp = toRationalPolynomial(denominator, variable, builtins, {256, 2048});
    if (!np || !dp || np->isZero() || dp->isZero())
        return std::nullopt;

    if (np->degree() < dp->degree())
        return integer(0);

    const Rational ratio = np->coefficient(np->degree()) / dp->coefficient(dp->degree());
    if (np->degree() == dp->degree())
        return rational(ratio);

    const std::size_t degreeDifference = np->degree() - dp->degree();
    int sign = rationalSign(ratio);
    if (negativeInfinityPoint && (degreeDifference % 2) != 0)
        sign = -sign;
    return signedInfinity(sign, builtins, mathematics, angles, infinity, assumptions);
}

[[nodiscard]] std::optional<Expr> quadraticRadicalInfiniteLimit(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    if (!isHead(expression, builtins, BuiltinId::Sqrt)
        || expression.asCall().arguments.size() != 1)
        return std::nullopt;
    const auto polynomial = toRationalPolynomial(
        expression.asCall().arguments[0], variable, builtins, {2, 16});
    if (!polynomial || polynomial->degree() != 2
        || rationalSign(polynomial->coefficient(2)) <= 0)
        return std::nullopt;
    // 正の二次先頭係数ならP(x)>0が十分遠方で保証され，principal sqrtは
    // ±Infinityの双方で非負の大きさInfinityへ発散する。
    return Expr{infinity};
}

[[nodiscard]] std::optional<Expr> quadraticRadicalRatioInfiniteLimit(
    const Expr& expression,
    const expression::Symbol& variable,
    bool negativeInfinityPoint,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!isHead(expression, builtins, BuiltinId::Divide)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;

    const auto quadraticRadicand = [&](const Expr& candidate)
        -> std::optional<RationalPolynomial> {
        if (!isHead(candidate, builtins, BuiltinId::Sqrt)
            || candidate.asCall().arguments.size() != 1)
            return std::nullopt;
        auto polynomial = toRationalPolynomial(
            candidate.asCall().arguments[0], variable, builtins, {2, 16});
        if (!polynomial || polynomial->degree() != 2
            || rationalSign(polynomial->coefficient(2)) <= 0)
            return std::nullopt;
        return polynomial;
    };

    bool radicalInNumerator = false;
    std::optional<RationalPolynomial> quadratic = quadraticRadicand(arguments[0]);
    std::optional<RationalPolynomial> affine;
    if (quadratic) {
        radicalInNumerator = true;
        affine = toRationalPolynomial(arguments[1], variable, builtins, {1, 8});
    }
    else {
        quadratic = quadraticRadicand(arguments[1]);
        if (quadratic)
            affine = toRationalPolynomial(arguments[0], variable, builtins, {1, 8});
    }
    if (!quadratic || !affine || affine->degree() != 1
        || affine->coefficient(1).isZero())
        return std::nullopt;

    Expr leadingRoot = simplify(
        call(builtins, BuiltinId::Sqrt,
            {rational(quadratic->coefficient(2))}),
        builtins, mathematics, angles, assumptions);
    Expr result = radicalInNumerator
        ? divide(
            std::move(leadingRoot), rational(affine->coefficient(1)),
            builtins, mathematics, angles, assumptions)
        : divide(
            rational(affine->coefficient(1)), std::move(leadingRoot),
            builtins, mathematics, angles, assumptions);
    if (negativeInfinityPoint)
        result = negate(
            std::move(result), builtins, mathematics, angles, assumptions);
    return fullSimplify(
        std::move(result), builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr imaginaryPi(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
    Expr imaginaryUnit{Number::complex(
        numeric::RealNumber{BigInt{0}}, numeric::RealNumber{BigInt{1}})};
    return simplify(call(builtins, BuiltinId::Multiply, {
        std::move(imaginaryUnit), Expr{pi->symbol}}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr inverseHalfTurn(
    bool negative,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr value = [&]() -> Expr {
        switch (angles.defaultUnit()) {
        case mathematics::AngleUnit::Degree:
            return integer(90);
        case mathematics::AngleUnit::Gradian:
            return integer(100);
        case mathematics::AngleUnit::Radian: {
            const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
            return divide(Expr{pi->symbol}, integer(2),
                builtins, mathematics, angles, assumptions);
        }
        }
        return integer(0);
    }();
    return negative
        ? negate(std::move(value), builtins, mathematics, angles, assumptions)
        : value;
}

[[nodiscard]] std::optional<Rational> affineSlopeAtInfinity(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    const auto polynomial = toRationalPolynomial(expression, variable, builtins, {1, 4});
    if (!polynomial || polynomial->degree() != 1)
        return std::nullopt;
    return polynomial->coefficient(1);
}

[[nodiscard]] std::optional<int> polynomialSignAtInfinity(
    const Expr& expression,
    const expression::Symbol& variable,
    bool negativeInfinityPoint,
    const evaluation::BuiltinRegistry& builtins) {
    const auto polynomial = toRationalPolynomial(
        expression, variable, builtins, {256, 2048});
    if (!polynomial || polynomial->isZero() || polynomial->degree() == 0)
        return std::nullopt;
    int sign = rationalSign(polynomial->coefficient(polynomial->degree()));
    if (negativeInfinityPoint && (polynomial->degree() % 2) != 0)
        sign = -sign;
    return sign;
}


[[nodiscard]] bool isExactRealRationalFunction(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins);

[[nodiscard]] bool isSignedInfinity(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    return isInfinity(expression, infinity)
        || isNegativeInfinity(expression, builtins, infinity);
}

[[nodiscard]] bool isNamedLimitSentinel(
    const Expr& expression,
    const expression::Symbol* symbol) {
    return symbol && expression.isSymbol()
        && expression.asSymbol().sameIdentity(*symbol);
}

[[nodiscard]] bool isExceptionalLimitValue(
    const Expr& expression,
    const expression::Symbol* complexInfinity,
    const expression::Symbol* indeterminate) {
    return isNamedLimitSentinel(expression, complexInfinity)
        || isNamedLimitSentinel(expression, indeterminate);
}

[[nodiscard]] bool containsLimitSentinel(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity,
    const expression::Symbol* complexInfinity,
    const expression::Symbol* indeterminate) {
    if (isSignedInfinity(expression, builtins, infinity)
        || isExceptionalLimitValue(expression, complexInfinity, indeterminate))
        return true;
    if (expression.isCall()) {
        for (const Expr& argument : expression.asCall().arguments) {
            if (containsLimitSentinel(
                    argument, builtins, infinity, complexInfinity, indeterminate))
                return true;
        }
    }
    if (expression.isArray()) {
        const auto& array = expression.asArray();
        for (std::size_t i = 0; i < array.size(); ++i) {
            if (containsLimitSentinel(
                    array.element(i), builtins, infinity, complexInfinity, indeterminate))
                return true;
        }
    }
    if (expression.isList()) {
        for (const Expr& element : expression.asList().elements) {
            if (containsLimitSentinel(
                    element, builtins, infinity, complexInfinity, indeterminate))
                return true;
        }
    }
    return false;
}

[[nodiscard]] bool isExactOne(const Expr& expression) {
    const auto value = expression::exact::realRational(expression);
    return value && *value == Rational{BigInt{1}};
}

[[nodiscard]] std::optional<Expr> infinitySeriesLimit(
    const Expr& expression,
    const expression::Symbol& variable,
    bool negativeInfinityPoint,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions) {
    Expr prepared = expression;
    if (negativeInfinityPoint) {
        Expr reflectedVariable = simplify(
            call(builtins, BuiltinId::Negate, {Expr{variable}}),
            builtins, mathematics, angles, assumptions);
        auto reflected = substituteLimitFreeSymbol(
            expression, variable, reflectedVariable, builtins);
        if (!reflected)
            return std::nullopt;
        prepared = simplify(
            std::move(*reflected), builtins, mathematics, angles, assumptions);
    }

    constexpr std::size_t kSeriesOrder = 6;
    auto expanded = seriesExpression(
        prepared, variable, Expr{infinity}, kSeriesOrder,
        builtins, mathematics, angles, assumptions);
    if (!expanded) return std::nullopt;
    auto data = parseSeriesData(*expanded, builtins);
    if (!data || data->exponentDenominator == 0)
        return std::nullopt;

    for (std::size_t i = 0; i < data->coefficients.size(); ++i) {
        const std::int64_t exponentNumerator =
            data->minimumExponent + static_cast<std::int64_t>(i);
        const bool ordinaryNonZero = !isZero(data->coefficients[i]);
        bool logarithmicNonZero = false;
        for (const auto& layer : data->logarithmicCoefficients) {
            if (i < layer.size() && !isZero(layer[i])) {
                logarithmicNonZero = true;
                break;
            }
        }
        if (!ordinaryNonZero && !logarithmicNonZero)
            continue;

        // t=1/x -> 0+なので，負次数は無限大へ発散する。
        // log層が同じ先頭次数にある場合は符号解析が必要になるため保守的に未解決とする。
        if (exponentNumerator < 0) {
            if (logarithmicNonZero)
                return std::nullopt;
            const auto coefficient = expression::exact::realRational(data->coefficients[i]);
            if (!coefficient || coefficient->isZero())
                return std::nullopt;
            return signedInfinity(
                rationalSign(*coefficient),
                builtins, mathematics, angles, infinity, assumptions);
        }
        if (exponentNumerator == 0 && logarithmicNonZero)
            return std::nullopt;
        if (exponentNumerator == 0 && ordinaryNonZero)
            return fullSimplify(
                data->coefficients[i], builtins, mathematics, angles, assumptions);
        if (exponentNumerator > 0)
            return integer(0);
    }
    return integer(0);
}

[[nodiscard]] bool isOscillatoryPeriodicBuiltin(BuiltinId id) {
    return id == BuiltinId::Sin || id == BuiltinId::Cos || id == BuiltinId::Tan;
}

[[nodiscard]] bool isExactRealRationalFunction(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        const auto numerator = toRationalPolynomial(
            expression.asCall().arguments[0], variable, builtins, {256, 2048});
        const auto denominator = toRationalPolynomial(
            expression.asCall().arguments[1], variable, builtins, {256, 2048});
        return numerator && denominator && !denominator->isZero();
    }
    return toRationalPolynomial(expression, variable, builtins, {256, 2048}).has_value();
}

[[nodiscard]] std::optional<Expr> reduceExponentiallyScaledRationalProductAtInfinity(
    const Expr& expression,
    const expression::Symbol& variable,
    bool negativeInfinityPoint,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Multiply))
        return std::nullopt;

    std::vector<Expr> exponents;
    std::vector<Expr> rationalFactors;
    for (const Expr& factor : expression.asCall().arguments) {
        if (isHead(factor, builtins, BuiltinId::Exp)
            && factor.asCall().arguments.size() == 1) {
            exponents.push_back(factor.asCall().arguments[0]);
            continue;
        }
        // exp(q(x))の減衰は任意のexact実有理函数より速い。分母も
        // 非零多項式と証明できるものだけを許し，函数係数や近似値へは広げない。
        if (!isExactRealRationalFunction(factor, variable, builtins))
            return std::nullopt;
        rationalFactors.push_back(factor);
    }
    if (exponents.empty())
        return std::nullopt;

    const std::size_t exponentCount = exponents.size();
    Expr combinedExponent = exponentCount == 1
        ? exponents.front()
        : call(builtins, BuiltinId::Add, std::move(exponents));
    const auto sign = polynomialSignAtInfinity(
        combinedExponent, variable, negativeInfinityPoint, builtins);
    if (sign && *sign < 0)
        return integer(0);

    // exp(p)exp(q)=exp(p+q)はbranch条件のないentire函数恒等式である。
    // 指数の変数依存性がexactに相殺した場合は，0*Infinityという人工的な
    // 不定形を避け，残る有理函数（および定数exp）だけを再度limitへ渡す。
    if (exponentCount < 2)
        return std::nullopt;
    const auto exponentPolynomial = toRationalPolynomial(
        combinedExponent, variable, builtins, {256, 2048});
    if (!exponentPolynomial || exponentPolynomial->degree() != 0)
        return std::nullopt;
    if (!exponentPolynomial->coefficient(0).isZero())
        rationalFactors.push_back(call(builtins, BuiltinId::Exp, {
            Expr{Number{exponentPolynomial->coefficient(0)}}}));
    if (rationalFactors.empty())
        return integer(1);
    if (rationalFactors.size() == 1)
        return rationalFactors.front();
    return call(builtins, BuiltinId::Multiply, std::move(rationalFactors));
}

[[nodiscard]] bool isBoundedRealTrigFactor(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall() || expression.asCall().arguments.size() != 1)
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || (definition->id != BuiltinId::Sin && definition->id != BuiltinId::Cos))
        return false;
    return isExactRealRationalFunction(expression.asCall().arguments[0], variable, builtins);
}

[[nodiscard]] std::optional<Rational> affinePlusBoundedTrigSlope(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (const auto polynomial = toRationalPolynomial(
            expression, variable, builtins, {1, 8});
        polynomial && polynomial->degree() <= 1)
        return polynomial->coefficient(1);
    if (isBoundedRealTrigFactor(expression, variable, builtins))
        return Rational{BigInt{0}};
    if (!expression.isCall())
        return std::nullopt;

    const auto* definition = builtins.find(expression.asCall().head);
    const auto& arguments = expression.asCall().arguments;
    if (!definition)
        return std::nullopt;
    if (definition->id == BuiltinId::Negate && arguments.size() == 1) {
        if (auto slope = affinePlusBoundedTrigSlope(arguments[0], variable, builtins))
            return -*slope;
        return std::nullopt;
    }
    if ((definition->id == BuiltinId::Add || definition->id == BuiltinId::Subtract)
        && !arguments.empty()) {
        Rational slope{BigInt{0}};
        for (std::size_t i = 0; i < arguments.size(); ++i) {
            const auto termSlope = affinePlusBoundedTrigSlope(
                arguments[i], variable, builtins);
            if (!termSlope)
                return std::nullopt;
            slope = slope + (definition->id == BuiltinId::Subtract && i == 1
                ? -*termSlope : *termSlope);
        }
        return slope;
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> rationalizeQuadraticRadicalAtPositiveInfinity(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    bool scaleForLimit) {
    const auto signedRadicand = [&](const Expr& term)
        -> std::optional<std::pair<int, Expr>> {
        if (isHead(term, builtins, BuiltinId::Sqrt)
            && term.asCall().arguments.size() == 1)
            return std::pair<int, Expr>{1, term.asCall().arguments[0]};
        if (isHead(term, builtins, BuiltinId::Negate)
            && term.asCall().arguments.size() == 1) {
            const Expr& inner = term.asCall().arguments[0];
            if (isHead(inner, builtins, BuiltinId::Sqrt)
                && inner.asCall().arguments.size() == 1)
                return std::pair<int, Expr>{-1, inner.asCall().arguments[0]};
        }
        return std::nullopt;
    };

    int radicalSign = 0;
    std::optional<Expr> radicand;
    std::optional<Expr> affine;
    if (isHead(expression, builtins, BuiltinId::Subtract)
        && expression.asCall().arguments.size() == 2) {
        const Expr& lhs = expression.asCall().arguments[0];
        const Expr& rhs = expression.asCall().arguments[1];
        if (auto lhsRadical = signedRadicand(lhs)) {
            radicalSign = lhsRadical->first;
            radicand = lhsRadical->second;
            affine = simplify(
                call(builtins, BuiltinId::Negate, {rhs}),
                builtins, mathematics, angles, assumptions);
        }
        else if (auto rhsRadical = signedRadicand(rhs)) {
            radicalSign = -rhsRadical->first;
            radicand = rhsRadical->second;
            affine = lhs;
        }
        else {
            return std::nullopt;
        }
    }
    else if (isHead(expression, builtins, BuiltinId::Add)
        && expression.asCall().arguments.size() == 2) {
        const Expr& first = expression.asCall().arguments[0];
        const Expr& second = expression.asCall().arguments[1];
        if (auto firstRadical = signedRadicand(first)) {
            radicalSign = firstRadical->first;
            radicand = firstRadical->second;
            affine = second;
        }
        else if (auto secondRadical = signedRadicand(second)) {
            radicalSign = secondRadical->first;
            radicand = secondRadical->second;
            affine = first;
        }
        else {
            return std::nullopt;
        }
    }
    else {
        return std::nullopt;
    }

    if (!radicand || !affine)
        return std::nullopt;
    const Expr& radicandExpr = *radicand;
    const Expr& affineExpr = *affine;
    const auto p = toRationalPolynomial(radicandExpr, variable, builtins, {2, 16});
    const auto a = toRationalPolynomial(affineExpr, variable, builtins, {1, 8});
    if (!p || !a || p->degree() != 2 || a->degree() != 1)
        return std::nullopt;

    const Rational affineLeading = a->coefficient(1);
    if (p->coefficient(2) != affineLeading * affineLeading
        || rationalSign(affineLeading) != -radicalSign)
        return std::nullopt;

    // (-x)^2のような未簡約Powerを残すとInfinity-Infinityが再発するため，
    // 既に証明済みの有理多項式係数上でP-affine^2をexactに作る。
    std::vector<Rational> residualCoefficients(3, Rational{BigInt{0}});
    for (std::size_t exponent = 0; exponent < residualCoefficients.size(); ++exponent) {
        Rational affineSquare{BigInt{0}};
        for (std::size_t left = 0; left <= exponent; ++left)
            affineSquare = affineSquare
                + a->coefficient(left) * a->coefficient(exponent - left);
        residualCoefficients[exponent] = p->coefficient(exponent) - affineSquare;
    }
    RationalPolynomial residualPolynomial{std::move(residualCoefficients)};
    // P=affine^2で符号条件も満たす場合，principal sqrtは十分大きい正のxで
    // -radicalSign*affineそのものとなる。0/xという人工的な定義域穴を作らず，
    // このeventual identityから極限0を直接返す。
    if (residualPolynomial.isZero())
        return integer(0);
    Expr residual = polynomialToExpandedExpr(
        residualPolynomial, variable, builtins);

    Expr conjugateRadical = call(builtins, BuiltinId::Sqrt, {radicandExpr});
    if (radicalSign < 0)
        conjugateRadical = call(
            builtins, BuiltinId::Negate, {std::move(conjugateRadical)});
    Expr conjugate = simplify(
        call(builtins, BuiltinId::Subtract,
            {std::move(conjugateRadical), affineExpr}),
        builtins, mathematics, angles, assumptions);
    if (!scaleForLimit)
        return call(
            builtins, BuiltinId::Divide,
            {std::move(residual), std::move(conjugate)});

    const Expr variableExpr{variable};
    Expr variableSquared = call(
        builtins, BuiltinId::Power, {variableExpr, integer(2)});
    Expr numerator = simplify(
        divide(std::move(residual), variableExpr,
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);

    Expr scaledRadicand = divide(
        radicandExpr, std::move(variableSquared),
        builtins, mathematics, angles, assumptions);
    Expr scaledRadical = call(
        builtins, BuiltinId::Sqrt, {std::move(scaledRadicand)});
    if (radicalSign < 0)
        scaledRadical = call(
            builtins, BuiltinId::Negate, {std::move(scaledRadical)});
    Expr scaledAffine = divide(
        affineExpr, variableExpr,
        builtins, mathematics, angles, assumptions);
    Expr denominator = simplify(
        call(builtins, BuiltinId::Subtract,
            {std::move(scaledRadical), std::move(scaledAffine)}),
        builtins, mathematics, angles, assumptions);
    return call(
        builtins, BuiltinId::Divide,
        {std::move(numerator), std::move(denominator)});
}

[[nodiscard]] std::optional<Expr> rationalizeQuadraticRadicalSubexpressionAtPositiveInfinity(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    bool scaleDirect) {
    if (auto direct = rationalizeQuadraticRadicalAtPositiveInfinity(
            expression, variable, builtins, mathematics, angles, assumptions, scaleDirect))
        return direct;
    if (!expression.isCall())
        return std::nullopt;

    const auto& source = expression.asCall();
    for (std::size_t i = 0; i < source.arguments.size(); ++i) {
        auto transformed = rationalizeQuadraticRadicalSubexpressionAtPositiveInfinity(
            source.arguments[i], variable,
            builtins, mathematics, angles, assumptions, false);
        if (!transformed)
            continue;
        std::vector<Expr> arguments = source.arguments;
        arguments[i] = std::move(*transformed);
        return simplify(
            Expr::rebuildCall(source, std::move(arguments)),
            builtins, mathematics, angles, assumptions);
    }
    return std::nullopt;
}

[[nodiscard]] Expr limitCore(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& point,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinity,
    const mathematics::AssumptionSet& assumptions,
    const expression::Symbol* complexInfinity,
    const expression::Symbol* indeterminate,
    std::size_t depth) {
    if (depth > maximumLimitDepth)
        return unresolved(expression, variable, point, direction, builtins);
    if (!containsFreeLimitSymbol(expression, variable, builtins))
        return expression;

    // brace値の極限は成分ごとに独立である。ragged Listをgeneric substitutionへ
    // 落とすとshape判定やdefinedness境界で元式へ戻るため，Array/Listともここで明示的に扱う。
    if (expression.isArray()) {
        const auto& array = expression.asArray();
        std::vector<Expr> elements;
        elements.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i)
            elements.push_back(limitCore(
                array.element(i), variable, point, direction, builtins, mathematics, angles,
                infinity, assumptions, complexInfinity, indeterminate, depth + 1));
        return Expr::array(array.shape, std::move(elements));
    }
    if (expression.isList()) {
        std::vector<Expr> elements;
        elements.reserve(expression.asList().elements.size());
        for (const Expr& element : expression.asList().elements)
            elements.push_back(limitCore(
                element, variable, point, direction, builtins, mathematics, angles,
                infinity, assumptions, complexInfinity, indeterminate, depth + 1));
        return expression::braceValue(std::move(elements));
    }

    // 片側極限のassumptionでabs等が簡約されると，held Dの被微分式だけが
    // 書き換わってD自体は再dispatchされない。D[x,x]のように安全に確定できる
    // 場合だけここでmaterializeし，未対応の微分は従来どおりheldのまま残す。
    if (isHead(expression, builtins, BuiltinId::Derivative)) {
        if (auto derivative = materializeHeldDerivative(
                expression, builtins, mathematics, angles);
            derivative && *derivative != expression)
            return limitCore(*derivative, variable, point, direction, builtins, mathematics, angles,
                infinity, assumptions, complexInfinity, indeterminate, depth + 1);
    }

    // vector calculusもHoldAllなので，直接点代入より先に既存kernelでmaterializeする。
    // 先にx->0すると grad[x^2,{x}] が grad[0,{0}] へ壊れる。
    if (auto calculus = materializeHeldVectorCalculus(
            expression, builtins, mathematics, angles);
        calculus && *calculus != expression)
        return limitCore(*calculus, variable, point, direction, builtins, mathematics, angles,
            infinity, assumptions, complexInfinity, indeterminate, depth + 1);

    // algebra変形のcontrol variableやSeriesのiteratorへ点代入する前に，
    // 公開frontendを既存kernelで一度materializeする。collect[0,0] 等の構文破壊と，
    // normal[series[...]] がLimit内だけ未接続になる非対称性を防ぐ。
    if (auto transformed = materializeHeldAlgebraTransform(
            expression, builtins, mathematics, angles);
        transformed && *transformed != expression)
        return limitCore(*transformed, variable, point, direction, builtins, mathematics, angles,
            infinity, assumptions, complexInfinity, indeterminate, depth + 1);
    if (auto normalized = materializeHeldNormal(
            expression, builtins, mathematics, angles);
        normalized && *normalized != expression)
        return limitCore(*normalized, variable, point, direction, builtins, mathematics, angles,
            infinity, assumptions, complexInfinity, indeterminate, depth + 1);

    // casesの条件が極限変数に依存しない場合だけ，各branchの極限へ分配できる。
    // 変数依存条件ではapproach directionとbranch境界の解析が必要であり，
    // 点代入でbranch値を先に評価すると0/0等を誤ってDomainErrorへ落とすため未解決で保持する。
    if (isHead(expression, builtins, BuiltinId::Cases)) {
        const auto& sourceBranches = expression.asCall().arguments;
        std::vector<Expr> branches;
        branches.reserve(sourceBranches.size());
        for (const Expr& branchExpression : sourceBranches) {
            if (!isHead(branchExpression, builtins, BuiltinId::CaseBranch)
                || branchExpression.asCall().arguments.empty()
                || branchExpression.asCall().arguments.size() > 2)
                return unresolved(expression, variable, point, direction, builtins);
            const auto& branch = branchExpression.asCall().arguments;
            if (branch.size() == 2 && containsSymbol(branch[1], variable))
                return unresolved(expression, variable, point, direction, builtins);
            Expr value = limitCore(
                branch[0], variable, point, direction, builtins, mathematics, angles,
                infinity, assumptions, complexInfinity, indeterminate, depth + 1);
            if (branch.size() == 2)
                branches.push_back(detail::makeCaseBranch(builtins, std::move(value), branch[1]));
            else
                branches.push_back(detail::makeCaseBranch(builtins, std::move(value)));
        }
        return detail::makeCases(builtins, std::move(branches));
    }

    // 実引数が無限大へ走る周期三角函数は単一の極限値を持たない。
    // これは「未実装」ではなく不存在を証明できる場合なので，Evaluatorから
    // Indeterminate atomが渡されていればそれを返す。有限点の二側極限では，
    // どちらか一方で無限振動を証明できれば二側極限の不存在も確定する。
    if (expression.isCall() && expression.asCall().arguments.size() == 1) {
        const auto* definition = builtins.find(expression.asCall().head);
        if (definition && isOscillatoryPeriodicBuiltin(definition->id) && indeterminate) {
            const Expr& argument = expression.asCall().arguments[0];
            if (direction == LimitDirection::TwoSided
                && !isInfinity(point, infinity)
                && !isNegativeInfinity(point, builtins, infinity)) {
                Expr left = limitCore(argument, variable, point, LimitDirection::Left,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity,
                    indeterminate, depth + 1);
                if (isSignedInfinity(left, builtins, infinity))
                    return Expr{*indeterminate};
                Expr right = limitCore(argument, variable, point, LimitDirection::Right,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity,
                    indeterminate, depth + 1);
                if (isSignedInfinity(right, builtins, infinity))
                    return Expr{*indeterminate};
            }
            else {
                Expr inner = limitCore(argument, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity,
                    indeterminate, depth + 1);
                if (isSignedInfinity(inner, builtins, infinity))
                    return Expr{*indeterminate};
            }
        }
    }

    const bool atPositiveInfinity = isInfinity(point, infinity);
    const bool atNegativeInfinity = isNegativeInfinity(point, builtins, infinity);

    // base->1なら，十分近くではprincipal Logのbranch cutと0を避ける。
    // 有限点でもbase^exponent=exp(exponent Log(base))へ写し，exponent単独の
    // 極限が存在しない1^Infinity型を積の極限として証明する。
    if (isHead(expression, builtins, BuiltinId::Power)
        && expression.asCall().arguments.size() == 2) {
        const Expr& base = expression.asCall().arguments[0];
        const Expr& exponent = expression.asCall().arguments[1];
        Expr baseLimit = limitCore(
            base, variable, point, direction,
            builtins, mathematics, angles, infinity, assumptions,
            complexInfinity, indeterminate, depth + 1);
        if (isExactOne(baseLimit)) {
            Expr logarithmicExponent = simplify(
                call(builtins, BuiltinId::Multiply,
                    {exponent, call(builtins, BuiltinId::Log, {base})}),
                builtins, mathematics, angles, assumptions);
            Expr exponentResult = limitCore(
                logarithmicExponent, variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions,
                complexInfinity, indeterminate, depth + 1);
            if (isInfinity(exponentResult, infinity))
                return Expr{infinity};
            if (isNegativeInfinity(exponentResult, builtins, infinity))
                return integer(0);
            if (!isHead(exponentResult, builtins, BuiltinId::Limit)
                && !isExceptionalLimitValue(
                    exponentResult, complexInfinity, indeterminate)) {
                return fullSimplify(
                    call(builtins, BuiltinId::Exp, {std::move(exponentResult)}),
                    builtins, mathematics, angles, assumptions);
            }
        }
    }

    // 線形な外側構造は先に分解する。improper integralの原始函数で
    // -exp[-x] や atan[x]-atan[0] を無用に未評価へ落とさない。
    if (expression.isCall()) {
        const auto* outer = builtins.find(expression.asCall().head);
        const auto& arguments = expression.asCall().arguments;
        if (outer && outer->id == BuiltinId::Negate && arguments.size() == 1) {
            Expr inner = limitCore(arguments[0], variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
            if (isExceptionalLimitValue(inner, complexInfinity, indeterminate))
                return inner;
            if (!isHead(inner, builtins, BuiltinId::Limit))
                return negate(std::move(inner), builtins, mathematics, angles, assumptions);
        }
        if (outer && (outer->id == BuiltinId::Add || outer->id == BuiltinId::Subtract)) {
            std::vector<Expr> values;
            values.reserve(arguments.size());
            bool decomposable = true;
            bool positiveInfinitySeen = false;
            bool negativeInfinitySeen = false;
            bool boundedUnresolvedSeen = false;
            for (std::size_t i = 0; i < arguments.size(); ++i) {
                const Expr& argument = arguments[i];
                Expr value = limitCore(argument, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
                if (isHead(value, builtins, BuiltinId::Limit)
                    || isExceptionalLimitValue(value, complexInfinity, indeterminate)) {
                    if ((atPositiveInfinity || atNegativeInfinity)
                        && isBoundedRealTrigFactor(argument, variable, builtins)) {
                        boundedUnresolvedSeen = true;
                        continue;
                    }
                    decomposable = false;
                    break;
                }
                if (outer->id == BuiltinId::Subtract && i == 1)
                    value = negate(std::move(value), builtins, mathematics, angles, assumptions);
                positiveInfinitySeen |= isInfinity(value, infinity);
                negativeInfinitySeen |= isNegativeInfinity(value, builtins, infinity);
                values.push_back(std::move(value));
            }
            if (boundedUnresolvedSeen && !positiveInfinitySeen && !negativeInfinitySeen)
                decomposable = false;
            if (decomposable && !(positiveInfinitySeen && negativeInfinitySeen)) {
                if (positiveInfinitySeen)
                    return Expr{infinity};
                if (negativeInfinitySeen)
                    return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);
                const BuiltinId combinedId = outer->id == BuiltinId::Subtract
                    ? BuiltinId::Add : outer->id;
                return simplify(call(builtins, combinedId, std::move(values)),
                    builtins, mathematics, angles, assumptions);
            }
            // 未解決項やInfinity-Infinityは，後段の局所Seriesで相殺を調べる余地を残す。
        }
        if (outer && outer->id == BuiltinId::Divide && arguments.size() == 2) {
            Expr numerator = limitCore(arguments[0], variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
            Expr denominator = limitCore(arguments[1], variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
            const bool numeratorFinite = !isHead(numerator, builtins, BuiltinId::Limit)
                && !isInfinity(numerator, infinity)
                && !isNegativeInfinity(numerator, builtins, infinity)
                && !isExceptionalLimitValue(numerator, complexInfinity, indeterminate);
            const bool denominatorFinite = !isHead(denominator, builtins, BuiltinId::Limit)
                && !isInfinity(denominator, infinity)
                && !isNegativeInfinity(denominator, builtins, infinity)
                && !isExceptionalLimitValue(denominator, complexInfinity, indeterminate);
            if (numeratorFinite && denominatorFinite && !isZero(denominator))
                return divide(std::move(numerator), std::move(denominator),
                    builtins, mathematics, angles, assumptions);

            const bool denominatorInfinite = isSignedInfinity(
                denominator, builtins, infinity);
            if (denominatorInfinite
                && (numeratorFinite
                    || isBoundedRealTrigFactor(arguments[0], variable, builtins)))
                return integer(0);

            if (isSignedInfinity(numerator, builtins, infinity)
                && denominatorFinite) {
                if (const auto finiteDenominator = expression::exact::realRational(denominator);
                    finiteDenominator && !finiteDenominator->isZero()) {
                    int sign = isNegativeInfinity(numerator, builtins, infinity) ? -1 : 1;
                    sign *= rationalSign(*finiteDenominator);
                    return signedInfinity(
                        sign, builtins, mathematics, angles, infinity, assumptions);
                }
            }

            // 共通分母が±Infinityへ走る和は項別の商へ分配し，各項の極限が
            // すべて証明できた場合だけ再結合する。bounded/xのsqueezeもここで効く。
            if (denominatorInfinite
                && (isHead(arguments[0], builtins, BuiltinId::Add)
                    || isHead(arguments[0], builtins, BuiltinId::Subtract))) {
                const bool subtraction = isHead(
                    arguments[0], builtins, BuiltinId::Subtract);
                std::vector<Expr> values;
                bool decomposable = true;
                bool positiveInfinitySeen = false;
                bool negativeInfinitySeen = false;
                for (const Expr& term : arguments[0].asCall().arguments) {
                    Expr quotient = call(
                        builtins, BuiltinId::Divide, {term, arguments[1]});
                    Expr value = limitCore(
                        quotient, variable, point, direction,
                        builtins, mathematics, angles, infinity, assumptions,
                        complexInfinity, indeterminate, depth + 1);
                    if (isHead(value, builtins, BuiltinId::Limit)
                        || isExceptionalLimitValue(
                            value, complexInfinity, indeterminate)) {
                        decomposable = false;
                        break;
                    }
                    if (subtraction && values.size() == 1)
                        value = negate(
                            std::move(value), builtins, mathematics, angles, assumptions);
                    positiveInfinitySeen |= isInfinity(value, infinity);
                    negativeInfinitySeen |= isNegativeInfinity(value, builtins, infinity);
                    values.push_back(std::move(value));
                }
                if (decomposable && !(positiveInfinitySeen && negativeInfinitySeen)) {
                    if (positiveInfinitySeen)
                        return Expr{infinity};
                    if (negativeInfinitySeen)
                        return negate(
                            Expr{infinity}, builtins, mathematics, angles, assumptions);
                    return simplify(
                        call(builtins, BuiltinId::Add, std::move(values)),
                        builtins, mathematics, angles, assumptions);
                }
            }
        }
        if (outer && outer->id == BuiltinId::Multiply) {
            // 実有理函数を引数に取るsin/cosは実軸上で絶対値1以下なので，
            // もう一方のfactorが0へ収束する2因子積はsqueeze theoremで0となる。
            if (arguments.size() == 2) {
                for (std::size_t boundedIndex = 0; boundedIndex < 2; ++boundedIndex) {
                    if (!isBoundedRealTrigFactor(arguments[boundedIndex], variable, builtins))
                        continue;
                    Expr vanishing = limitCore(arguments[1 - boundedIndex], variable, point, direction,
                        builtins, mathematics, angles, infinity, assumptions, complexInfinity,
                        indeterminate, depth + 1);
                    if (isZero(vanishing))
                        return integer(0);
                }
            }

            // すべてのfactorが有限極限へ収束する場合は通常の積へ戻す。
            // また，±Infinityとexact real rationalな非零有限因子だけから成る場合は，
            // 符号をexactに合成して無限大を確定する。0*Infinityや符号不明な有限因子は
            // 不定形として後段のSeries等へ残す。
            std::vector<Expr> values;
            values.reserve(arguments.size());
            bool unresolvedFactor = false;
            std::size_t infinityFactors = 0;
            int productSign = 1;
            bool finiteFactorsHaveKnownRealSign = true;
            bool zeroFactor = false;
            const mathematics::KnowledgeContext knowledge{
                builtins, mathematics, assumptions};
            for (const Expr& argument : arguments) {
                Expr value = limitCore(argument, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
                if (isHead(value, builtins, BuiltinId::Limit)
                    || isExceptionalLimitValue(value, complexInfinity, indeterminate)) {
                    unresolvedFactor = true;
                    break;
                }
                if (isInfinity(value, infinity)) {
                    ++infinityFactors;
                    continue;
                }
                if (isNegativeInfinity(value, builtins, infinity)) {
                    ++infinityFactors;
                    productSign = -productSign;
                    continue;
                }
                if (const auto rationalValue = expression::exact::realRational(value)) {
                    if (rationalValue->isZero())
                        zeroFactor = true;
                    else
                        productSign *= rationalSign(*rationalValue);
                }
                else {
                    const auto positive = knowledge.prove(mathematics::relation(
                        mathematics::RelationKind::Greater, value, integer(0)));
                    const auto negative = knowledge.prove(mathematics::relation(
                        mathematics::RelationKind::Less, value, integer(0)));
                    if (negative == mathematics::TruthValue::True)
                        productSign = -productSign;
                    else if (positive != mathematics::TruthValue::True)
                        finiteFactorsHaveKnownRealSign = false;
                }
                values.push_back(std::move(value));
            }
            if (!unresolvedFactor && infinityFactors == 0)
                return simplify(call(builtins, BuiltinId::Multiply, std::move(values)),
                    builtins, mathematics, angles, assumptions);
            if (!unresolvedFactor && infinityFactors > 0 && !zeroFactor
                && finiteFactorsHaveKnownRealSign)
                return signedInfinity(productSign,
                    builtins, mathematics, angles, infinity, assumptions);
        }
    }

    if (atPositiveInfinity || atNegativeInfinity) {
        if (auto rationalLimit = rationalFunctionInfiniteLimit(
                expression, variable, atNegativeInfinity,
                builtins, mathematics, angles, infinity, assumptions))
            return *rationalLimit;
        if (auto exponentialProduct = reduceExponentiallyScaledRationalProductAtInfinity(
                expression, variable, atNegativeInfinity, builtins)) {
            Expr result = limitCore(
                *exponentialProduct, variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions,
                complexInfinity, indeterminate, depth + 1);
            if (!isHead(result, builtins, BuiltinId::Limit))
                return result;
        }
        if (isHead(expression, builtins, BuiltinId::Divide)
            && expression.asCall().arguments.size() == 2) {
            const auto numeratorSlope = affinePlusBoundedTrigSlope(
                expression.asCall().arguments[0], variable, builtins);
            const auto denominatorSlope = affinePlusBoundedTrigSlope(
                expression.asCall().arguments[1], variable, builtins);
            if (numeratorSlope && denominatorSlope && !denominatorSlope->isZero())
                return rational(*numeratorSlope / *denominatorSlope);
        }
        if (auto radicalLimit = quadraticRadicalInfiniteLimit(
                expression, variable, builtins, infinity))
            return *radicalLimit;
        if (auto radicalRatio = quadraticRadicalRatioInfiniteLimit(
                expression, variable, atNegativeInfinity,
                builtins, mathematics, angles, assumptions))
            return *radicalRatio;

        // sqrt[quadratic] +/- affine のInfinity-Infinity型は共役で有理化し，
        // xでscaleしてから通常のlimit kernelへ戻す。+Infinityではx>0が最終的に
        // 保証されるため sqrt(P)/x = sqrt(P/x^2) をbranch安全に使える。
        Expr radicalCandidate = expression;
        if (atNegativeInfinity) {
            Expr reflectedVariable = simplify(
                call(builtins, BuiltinId::Negate, {Expr{variable}}),
                builtins, mathematics, angles, assumptions);
            if (auto reflected = substituteLimitFreeSymbol(
                    expression, variable, reflectedVariable, builtins))
                radicalCandidate = simplify(
                    std::move(*reflected), builtins, mathematics, angles, assumptions);
        }
        if (auto rationalized = rationalizeQuadraticRadicalSubexpressionAtPositiveInfinity(
                radicalCandidate, variable,
                builtins, mathematics, angles, assumptions, true)) {
            Expr result = limitCore(
                *rationalized, variable, Expr{infinity}, direction,
                builtins, mathematics, angles, infinity, assumptions,
                complexInfinity, indeterminate, depth + 1);
            if (!isHead(result, builtins, BuiltinId::Limit))
                return result;
        }

        // +Infinity seriesはt=1/x, t->0+へ写す既存Series kernelを共有する。
        // 先頭Laurent/Puiseux項の相殺が有限値へ落ちる場合だけexactに確定する。
        // -Infinityはx->-xで同じkernelへ送る。
        try {
            if (auto seriesLimit = infinitySeriesLimit(
                    expression, variable, atNegativeInfinity,
                    builtins, mathematics, angles, infinity, assumptions))
                return *seriesLimit;
        }
        catch (const error::CalcError& exception) {
            if (!recoverableSpeculativeError(exception.type()))
                throw;
        }

        // Classical integral functions have branch-sensitive but exact real-axis asymptotics.
        // Ei(x) ~ exp(x)/x for x->+Infinity and tends to zero for x->-Infinity.
        if (isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::ExponentialIntegralEi))
            return atNegativeInfinity ? integer(0) : Expr{infinity};

        // Principal Ci(x) tends to zero on the positive real axis. On the negative real
        // branch cut Ci(-r)=Ci(r)+I Pi (r>0), hence x->-Infinity tends exactly to I Pi.
        if (isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::CosineIntegralCi))
            return atNegativeInfinity
                ? imaginaryPi(builtins, mathematics, angles, assumptions)
                : integer(0);
        // DLMF 6.2.14: Si(x) -> +/- Pi/2 on the real axis.
        if (isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::SineIntegralSi)) {
            Expr halfPi = divide(Expr{mathematics.findConstant(mathematics::ConstantId::Pi)->symbol},
                integer(2), builtins, mathematics, angles, assumptions);
            return atNegativeInfinity
                ? negate(std::move(halfPi), builtins, mathematics, angles, assumptions)
                : halfPi;
        }

        // DLMF 7.2.9: both Fresnel C and S tend to 1/2 at +Infinity; both are odd.
        if (isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::FresnelC)
            || isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::FresnelS))
            return rational(atNegativeInfinity
                ? Rational{BigInt{-1}, BigInt{2}}
                : Rational{BigInt{1}, BigInt{2}});


        // li(x)=Ei(Log(x)). On the positive real axis li(x)->+Infinity.
        // Along x->-Infinity on the principal branch, Log(x)=log|x|+I Pi and
        // |Ei(Log(x))| ~ |x|/|Log(x)| -> Infinity while the value is complex.
        // mmCal has no DirectedInfinity yet, so ComplexInfinity is the exact conservative
        // extended-complex result: it records unbounded magnitude without inventing a real value.
        if (isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::LogarithmicIntegralLi)) {
            if (!atNegativeInfinity)
                return Expr{infinity};
            if (complexInfinity)
                return Expr{*complexInfinity};
        }

        if (expression.isCall() && expression.asCall().arguments.size() == 1) {
            const auto* definition = builtins.find(expression.asCall().head);
            const Expr& argument = expression.asCall().arguments[0];
            if (definition) {
                if (definition->id == BuiltinId::Exp) {
                    if (const auto sign = polynomialSignAtInfinity(
                            argument, variable, atNegativeInfinity, builtins))
                        return *sign > 0 ? Expr{infinity} : integer(0);
                }
                const auto slope = affineSlopeAtInfinity(argument, variable, builtins);
                const int infinitySign = atNegativeInfinity ? -1 : 1;
                if (slope && !slope->isZero()) {
                    const int argumentSign = rationalSign(*slope) * infinitySign;
                    switch (definition->id) {
                    case BuiltinId::Exp:
                        return argumentSign > 0 ? Expr{infinity} : integer(0);
                    case BuiltinId::Log:
                        if (argumentSign > 0)
                            return Expr{infinity};
                        break;
                    case BuiltinId::Atan:
                        return inverseHalfTurn(argumentSign < 0,
                            builtins, mathematics, angles, assumptions);
                    case BuiltinId::Tanh:
                    case BuiltinId::Erf:
                        return integer(argumentSign > 0 ? 1 : -1);
                    case BuiltinId::Erfc:
                        return integer(argumentSign > 0 ? 0 : 2);
                    case BuiltinId::Asinh:
                        return signedInfinity(argumentSign,
                            builtins, mathematics, angles, infinity, assumptions);
                    case BuiltinId::Sqrt:
                        // principal sqrtは正実軸上で非負であり，
                        // 引数が+Infinityへ走る場合だけ実の+Infinityへ発散する。
                        // 負実軸側はprincipal branch上で虚方向へ発散するため，
                        // DirectedInfinityを持たない現在は未解決に残す。
                        if (argumentSign > 0)
                            return Expr{infinity};
                        break;
                    default:
                        break;
                    }
                }
            }
        }

        // ここで確定しない商は，下段の安全な点代入とL'Hopital候補へ回す。
    }

    const auto rationalPoint = expression::exact::realRational(point);
    if (rationalPoint) {
        if (rationalPoint->isZero()) {
            // Ei(x)=gamma+log|x|+O(x) on the real axis, so both real one-sided
            // limits at its logarithmic singularity are -Infinity.
            if (isUnaryFunctionOfVariable(
                    expression, variable, builtins, BuiltinId::ExponentialIntegralEi))
                return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);

            // Ci(z)=gamma+Log(z)+O(z^2). From the right this is real -Infinity;
            // from the left the principal branch adds the bounded term I Pi. In either
            // case the directed limit is -Infinity, so the two-sided real limit agrees.
            if (isUnaryFunctionOfVariable(
                    expression, variable, builtins, BuiltinId::CosineIntegralCi))
                return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);

            // li(x)=Ei(Log(x)) tends to zero at the origin from either real side,
            // despite x=0 being a branch point of the complex principal function.
            if (isUnaryFunctionOfVariable(
                    expression, variable, builtins, BuiltinId::LogarithmicIntegralLi))
                return integer(0);
            // Si, Fresnel C and Fresnel S are entire odd functions and vanish at zero.
            if (isUnaryFunctionOfVariable(
                    expression, variable, builtins, BuiltinId::SineIntegralSi)
                || isUnaryFunctionOfVariable(
                    expression, variable, builtins, BuiltinId::FresnelC)
                || isUnaryFunctionOfVariable(
                    expression, variable, builtins, BuiltinId::FresnelS))
                return integer(0);
        }

        // li(x)=Ei(Log(x)) and Ei(t)->-Infinity as t->0 from either real side.
        // Hence both one-sided limits, and therefore the two-sided real limit, at x=1
        // are -Infinity even though li(1) itself is undefined.
        if (*rationalPoint == Rational{BigInt{1}}
            && isUnaryFunctionOfVariable(
                expression, variable, builtins, BuiltinId::LogarithmicIntegralLi))
            return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);

        // DLMF 4.37.24: principal atanh(x)=1/2 Log((1+x)/(1-x))。
        // 実定義域(-1,1)の内側からbranch pointへ近づく一側極限だけを
        // real Infinityとして確定する。反対側はprincipal branch cut上で
        // 虚部を伴うため，現在のInfinity sentinelでは表現しない。
        if (isUnaryFunctionOfVariable(expression, variable, builtins, BuiltinId::Atanh)) {
            if (*rationalPoint == Rational{BigInt{1}}
                && direction == LimitDirection::Left)
                return Expr{infinity};
            if (*rationalPoint == Rational{BigInt{-1}}
                && direction == LimitDirection::Right)
                return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);
        }

        if (auto rationalLimit = rationalFunctionFiniteLimit(
                expression, variable, *rationalPoint, direction,
                builtins, mathematics, angles, infinity, assumptions))
            return *rationalLimit;

        // principal Logの実軸右極限。左側は -Infinity + I Pi となるため、
        // 現在のInfinity sentinelで複素無限大を捏造せず未評価に残す。
        if (rationalPoint->isZero() && direction == LimitDirection::Right
            && isHead(expression, builtins, BuiltinId::Log)
            && expression.asCall().arguments.size() == 1
            && expression.asCall().arguments[0].isSymbol()
            && expression.asCall().arguments[0].asSymbol().sameIdentity(variable))
            return negate(Expr{infinity}, builtins, mathematics, angles, assumptions);

        // x^a log[x] -> 0 (x->0+, a>0)。improper integralのendpointで頻出し、
        // 実軸上の標準極限としてbranch安全に扱える。
        if (rationalPoint->isZero() && direction == LimitDirection::Right
            && isHead(expression, builtins, BuiltinId::Multiply)) {
            const auto& factors = expression.asCall().arguments;
            if (factors.size() == 2) {
                for (std::size_t logIndex = 0; logIndex < 2; ++logIndex) {
                    const Expr& logarithm = factors[logIndex];
                    const Expr& vanishing = factors[1 - logIndex];
                    if (!isHead(logarithm, builtins, BuiltinId::Log)
                        || logarithm.asCall().arguments.size() != 1
                        || !logarithm.asCall().arguments[0].isSymbol()
                        || !logarithm.asCall().arguments[0].asSymbol().sameIdentity(variable))
                        continue;
                    Rational exponent{BigInt{0}};
                    bool matched = false;
                    if (vanishing.isSymbol() && vanishing.asSymbol().sameIdentity(variable)) {
                        exponent = Rational{BigInt{1}};
                        matched = true;
                    }
                    else if (isHead(vanishing, builtins, BuiltinId::Power)
                        && vanishing.asCall().arguments.size() == 2
                        && vanishing.asCall().arguments[0].isSymbol()
                        && vanishing.asCall().arguments[0].asSymbol().sameIdentity(variable)) {
                        if (const auto value = expression::exact::realRational(vanishing.asCall().arguments[1])) {
                            exponent = *value;
                            matched = !exponent.isZero() && !exponent.numerator().isNegative();
                        }
                    }
                    if (matched)
                        return integer(0);
                }
            }
        }
    }

    // 既存の個別規則で決まらない有限点極限だけを，局所Seriesの先頭項で補完する。
    // 発散項やlog^kの定数次数が残る場合は推測せず，従来kernelへ処理を戻す。
    if (!isInfinity(point, infinity) && !isNegativeInfinity(point, builtins, infinity)) {
        try {
            if (auto seriesLimit = finiteLimitFromLocalSeries(
                    expression, variable, point,
                    builtins, mathematics, angles, assumptions))
                return *seriesLimit;
        }
        catch (const error::CalcError& exception) {
            if (!recoverableSpeculativeError(exception.type()))
                throw;
        }
    }

    if (auto substituted = substitutedLimitCandidate(
            expression, variable, point,
            builtins, mathematics, angles, assumptions)) {
        // Infinityを部分式に残した形式代入は極限の証明ではない。
        // exp[Infinity]/Infinity等を値として漏らさず，不定形解析へ回す。
        if (!containsFreeLimitSymbol(*substituted, variable, builtins)
            && !containsLimitSentinel(
                *substituted, builtins, infinity, complexInfinity, indeterminate))
            return fullSimplify(
                std::move(*substituted), builtins, mathematics, angles, assumptions);
    }

    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        Expr numerator = expression.asCall().arguments[0];
        Expr denominator = expression.asCall().arguments[1];
        for (std::size_t step = 0; step < maximumLHopitalSteps; ++step) {
            try {
                // 形式代入ではexp[Infinity]等がsentinelへ閉じないため，分子・分母の
                // 極限を既存kernelで証明し，0/0またはsigned Infinity/Infinityだけに適用する。
                Expr numeratorAt = limitCore(
                    numerator, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions,
                    complexInfinity, indeterminate, depth + 1);
                Expr denominatorAt = limitCore(
                    denominator, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions,
                    complexInfinity, indeterminate, depth + 1);
                const bool zeroOverZero = isZero(numeratorAt) && isZero(denominatorAt);
                const bool infinityOverInfinity =
                    isSignedInfinity(numeratorAt, builtins, infinity)
                    && isSignedInfinity(denominatorAt, builtins, infinity);
                if (!zeroOverZero && !infinityOverInfinity)
                    break;

                numerator = differentiateExpression(
                    numerator, variable, builtins, mathematics, angles);
                denominator = differentiateExpression(
                    denominator, variable, builtins, mathematics, angles);
                if (isHead(numerator, builtins, BuiltinId::Derivative)
                    || isHead(denominator, builtins, BuiltinId::Derivative))
                    break;
                Expr quotient = divide(numerator, denominator,
                    builtins, mathematics, angles, assumptions);
                Expr result = limitCore(
                    quotient, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions,
                    complexInfinity, indeterminate, depth + 1);
                if (!isHead(result, builtins, BuiltinId::Limit)
                    && !isExceptionalLimitValue(
                        result, complexInfinity, indeterminate))
                    return result;
            }
            catch (const error::CalcError& exception) {
                if (!recoverableSpeculativeError(exception.type()))
                    throw;
                break;
            }
        }
    }

    if (expression.isCall()) {
        const auto* definition = builtins.find(expression.asCall().head);
        const auto& arguments = expression.asCall().arguments;
        if (definition && definition->id == BuiltinId::Multiply && arguments.size() == 2) {
            // x*log(x) -> 0 (x->0+) など、0*Infinity型の代表的な形は商へ直してL'Hopitalへ送る。
            for (std::size_t logIndex = 0; logIndex < 2; ++logIndex) {
                const Expr& logarithm = arguments[logIndex];
                const Expr& factor = arguments[1 - logIndex];
                if (!isHead(logarithm, builtins, BuiltinId::Log)
                    || logarithm.asCall().arguments.size() != 1)
                    continue;
                Expr rewritten = divide(logarithm,
                    divide(integer(1), factor, builtins, mathematics, angles, assumptions),
                    builtins, mathematics, angles, assumptions);
                Expr result = limitCore(rewritten, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
                if (!isHead(result, builtins, BuiltinId::Limit))
                    return result;
            }
        }
    }

    return unresolved(expression, variable, point, direction, builtins);
}

} // namespace

Expr limitExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& point,
    LimitDirection direction,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions,
    const expression::Symbol* complexInfinitySymbol,
    const expression::Symbol* indeterminateSymbol) {
    if (direction == LimitDirection::TwoSided
        && isHead(expression, builtins, BuiltinId::Cases)
        && !isInfinity(point, infinitySymbol)
        && !isNegativeInfinity(point, builtins, infinitySymbol)) {
        // x依存Casesは点そのもののbranchを選ぶのではなく，左右から別々に
        // branchを確定して極限を比較する。境界値のclosed/open差を二側極限へ
        // 混ぜず，両側が同じと証明できた場合だけ値を返す。
        const auto directedLimit = [&](LimitDirection side) {
            mathematics::AssumptionSet directed = assumptions;
            directed.add(mathematics::elementOf(Expr{variable}, mathematics::NumericDomain::Real));
            directed.add(mathematics::relation(
                side == LimitDirection::Right
                    ? mathematics::RelationKind::Greater
                    : mathematics::RelationKind::Less,
                Expr{variable}, point));
            const Expr preparedSide = simplification::Simplifier{}.simplify(
                expression,
                simplification::SimplificationContext{builtins, mathematics, angles, directed});
            return limitCore(preparedSide, variable, point, side,
                builtins, mathematics, angles, infinitySymbol, directed,
                complexInfinitySymbol, indeterminateSymbol, 0);
        };

        Expr left = directedLimit(LimitDirection::Left);
        Expr right = directedLimit(LimitDirection::Right);
        if (!isHead(left, builtins, BuiltinId::Limit)
            && !isHead(right, builtins, BuiltinId::Limit)) {
            if (left == right)
                return left;
            if (indeterminateSymbol)
                return Expr{*indeterminateSymbol};
        }
    }

    mathematics::AssumptionSet local = assumptions;
    if (!isInfinity(point, infinitySymbol)
        && !isNegativeInfinity(point, builtins, infinitySymbol)
        && direction != LimitDirection::TwoSided) {
        local.add(mathematics::elementOf(Expr{variable}, mathematics::NumericDomain::Real));
        local.add(mathematics::relation(
            direction == LimitDirection::Right
                ? mathematics::RelationKind::Greater
                : mathematics::RelationKind::Less,
            Expr{variable}, point));
    }
    const Expr prepared = simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles, local});
    return limitCore(prepared, variable, point, direction,
        builtins, mathematics, angles, infinitySymbol, local, complexInfinitySymbol,
        indeterminateSymbol, 0);
}

} // namespace mmcal::symbolic
