// 極限limit
#include "symbolic/cases.hpp"
#include "limit.hpp"

#include "builtins/array_vector.hpp"
#include "builtins/polynomial_ideal.hpp"
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
                arguments[1] = Expr::array({spec.size()}, std::move(spec));
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

[[nodiscard]] std::optional<Expr> rationalFunctionInfiniteLimit(
    const Expr& expression,
    const expression::Symbol& variable,
    bool negativeInfinityPoint,
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


[[nodiscard]] bool isSignedInfinity(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const expression::Symbol& infinity) {
    return isInfinity(expression, infinity)
        || isNegativeInfinity(expression, builtins, infinity);
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

    // 線形な外側構造は先に分解する。improper integralの原始函数で
    // -exp[-x] や atan[x]-atan[0] を無用に未評価へ落とさない。
    if (expression.isCall()) {
        const auto* outer = builtins.find(expression.asCall().head);
        const auto& arguments = expression.asCall().arguments;
        if (outer && outer->id == BuiltinId::Negate && arguments.size() == 1) {
            Expr inner = limitCore(arguments[0], variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
            if (!isHead(inner, builtins, BuiltinId::Limit))
                return negate(std::move(inner), builtins, mathematics, angles, assumptions);
        }
        if (outer && (outer->id == BuiltinId::Add || outer->id == BuiltinId::Subtract)) {
            std::vector<Expr> values;
            values.reserve(arguments.size());
            bool decomposable = true;
            bool positiveInfinitySeen = false;
            bool negativeInfinitySeen = false;
            for (const Expr& argument : arguments) {
                Expr value = limitCore(argument, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
                if (isHead(value, builtins, BuiltinId::Limit)) {
                    decomposable = false;
                    break;
                }
                if (outer->id == BuiltinId::Subtract && values.size() == 1)
                    value = negate(std::move(value), builtins, mathematics, angles, assumptions);
                positiveInfinitySeen |= isInfinity(value, infinity);
                negativeInfinitySeen |= isNegativeInfinity(value, builtins, infinity);
                values.push_back(std::move(value));
            }
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
                && !isNegativeInfinity(numerator, builtins, infinity);
            const bool denominatorFinite = !isHead(denominator, builtins, BuiltinId::Limit)
                && !isInfinity(denominator, infinity)
                && !isNegativeInfinity(denominator, builtins, infinity);
            if (numeratorFinite && denominatorFinite && !isZero(denominator))
                return divide(std::move(numerator), std::move(denominator),
                    builtins, mathematics, angles, assumptions);
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

            // すべてのfactorが有限極限へ収束する場合だけ積を合成する。
            // 0*Infinity等の不定形はここで決めず、既存の専用ruleへ残す。
            std::vector<Expr> values;
            values.reserve(arguments.size());
            bool finite = true;
            for (const Expr& argument : arguments) {
                Expr value = limitCore(argument, variable, point, direction,
                    builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
                if (isHead(value, builtins, BuiltinId::Limit)
                    || isInfinity(value, infinity)
                    || isNegativeInfinity(value, builtins, infinity)) {
                    finite = false;
                    break;
                }
                values.push_back(std::move(value));
            }
            if (finite)
                return simplify(call(builtins, BuiltinId::Multiply, std::move(values)),
                    builtins, mathematics, angles, assumptions);
        }
    }

    if (atPositiveInfinity || atNegativeInfinity) {
        if (auto rationalLimit = rationalFunctionInfiniteLimit(
                expression, variable, atNegativeInfinity,
                builtins, mathematics, angles, infinity, assumptions))
            return *rationalLimit;

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

        return unresolved(expression, variable, point, direction, builtins);
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
        if (auto seriesLimit = finiteLimitFromLocalSeries(
                expression, variable, point,
                builtins, mathematics, angles, assumptions))
            return *seriesLimit;
    }

    if (auto rawSubstituted = substituteLimitFreeSymbol(
            expression, variable, point, builtins);
        rawSubstituted && domainConditionsHold(*rawSubstituted, builtins, mathematics, assumptions)) {
        Expr substituted = simplify(
            std::move(*rawSubstituted), builtins, mathematics, angles, assumptions);
        if (!containsFreeLimitSymbol(substituted, variable, builtins))
            return fullSimplify(std::move(substituted), builtins, mathematics, angles, assumptions);
    }

    if (isHead(expression, builtins, BuiltinId::Divide)
        && expression.asCall().arguments.size() == 2) {
        Expr numerator = expression.asCall().arguments[0];
        Expr denominator = expression.asCall().arguments[1];
        for (std::size_t step = 0; step < maximumLHopitalSteps; ++step) {
            const auto numeratorSubstituted = substituteLimitFreeSymbol(
                numerator, variable, point, builtins);
            const auto denominatorSubstituted = substituteLimitFreeSymbol(
                denominator, variable, point, builtins);
            if (!numeratorSubstituted || !denominatorSubstituted)
                break;
            Expr numeratorAt = simplify(*numeratorSubstituted,
                builtins, mathematics, angles, assumptions);
            Expr denominatorAt = simplify(*denominatorSubstituted,
                builtins, mathematics, angles, assumptions);
            if (!(isZero(numeratorAt) && isZero(denominatorAt)))
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
            Expr result = limitCore(quotient, variable, point, direction,
                builtins, mathematics, angles, infinity, assumptions, complexInfinity, indeterminate, depth + 1);
            if (!isHead(result, builtins, BuiltinId::Limit))
                return result;
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
