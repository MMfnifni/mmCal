// builtin dispatchを非再帰評価機械から分離する。
#include "evaluator.hpp"
#include "evaluation/builtin_categories.hpp"
#include "expression/array_utils.hpp"

#include "mathematics/assumption_parser.hpp"
#include "mathematics/value_facts.hpp"

#include "builtins/arithmetic.hpp"
#include "builtins/approximation_utilities.hpp"
#include "builtins/aggregate.hpp"
#include "builtins/statistics.hpp"
#include "builtins/signal_processing.hpp"
#include "builtins/array.hpp"
#include "builtins/array_vector.hpp"
#include "builtins/comparison.hpp"
#include "builtins/combinatorics.hpp"
#include "builtins/discrete_math.hpp"
#include "builtins/elementary_utilities.hpp"
#include "builtins/explain.hpp"
#include "builtins/complex_functions.hpp"
#include "builtins/hyperbolic.hpp"
#include "builtins/iteration.hpp"
#include "builtins/linear_algebra.hpp"
#include "builtins/numerical_calculus.hpp"
#include "builtins/polynomial_ideal.hpp"
#include "builtins/trigonometric.hpp"
#include "builtins/stable_elementary.hpp"
#include "builtins/special_functions.hpp"
#include "builtins/random_functions.hpp"
#include "builtins/transcendental.hpp"
#include "error/error_message.hpp"
#include "evaluation/iterator_spec.hpp"
#include "simplification/full_simplifier.hpp"
#include "simplification/simplifier.hpp"
#include "solver/absolute_value_solver.hpp"
#include "solver/polynomial_solver.hpp"
#include "solver/radical_solver.hpp"
#include "solver/real_equation_prover.hpp"
#include "solver/solve_constraints.hpp"
#include "solver/solve_normalization.hpp"
#include "solver/solver_support.hpp"
#include "solver/transcendental_solver.hpp"
#include "symbolic/algebraic_expression.hpp"
#include "symbolic/algebra_transforms.hpp"
#include "symbolic/algebraic_number.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/integration.hpp"
#include "symbolic/limit.hpp"
#include "symbolic/series.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/real_number.hpp"
#include "graphics/graphics_backend.hpp"
#include "formatting/expr_formatter.hpp"
#include "plot/plot_pipeline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace mmcal::evaluation {
namespace {

enum class InfiniteValue {
    None,
    Positive,
    Negative,
    Undirected
};

[[nodiscard]] bool isHeldSeriesPipeline(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins) {
    if (builtins.isCallTo(expression, BuiltinId::Series)
        || builtins.isCallTo(expression, BuiltinId::Derivative))
        return true;
    if (builtins.isCallTo(expression, BuiltinId::NumericalApproximation)) {
        const auto& arguments = expression.asCall().arguments;
        return !arguments.empty() && builtins.isCallTo(arguments.front(), BuiltinId::Series);
    }
    if ((!builtins.isCallTo(expression, BuiltinId::Normal)
            && !builtins.isCallTo(expression, BuiltinId::ToNormal))
        || expression.asCall().arguments.size() != 1)
        return false;
    return builtins.isCallTo(expression.asCall().arguments.front(), BuiltinId::Series);
}

[[nodiscard]] bool isPredefinedAtom(
    const expression::Expr& expression,
    symbols::PredefinedSymbolId id,
    const symbols::SymbolRegistry& symbolRegistry) noexcept {
    if (!expression.isSymbol())
        return false;
    if (const auto* definition = symbolRegistry.find(expression.asSymbol()); definition)
        return definition->id == id;
    const auto* canonical = symbolRegistry.find(id);
    return canonical && expression.asSymbol() == canonical->symbol;
}

[[nodiscard]] expression::Expr predefinedValue(
    symbols::PredefinedSymbolId id,
    const symbols::SymbolRegistry& symbolRegistry) {
    const auto* definition = symbolRegistry.find(id);
    if (!definition)
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "Required predefined exceptional value is not registered");
    return expression::Expr{definition->symbol};
}

[[nodiscard]] bool isIndeterminate(
    const expression::Expr& expression,
    const symbols::SymbolRegistry& symbolRegistry) noexcept {
    return isPredefinedAtom(
        expression, symbols::PredefinedSymbolId::Indeterminate, symbolRegistry);
}

[[nodiscard]] InfiniteValue infiniteValue(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins,
    const symbols::SymbolRegistry& symbolRegistry) noexcept {
    if (isPredefinedAtom(
            expression, symbols::PredefinedSymbolId::Infinity, symbolRegistry))
        return InfiniteValue::Positive;
    if (isPredefinedAtom(
            expression, symbols::PredefinedSymbolId::ComplexInfinity, symbolRegistry))
        return InfiniteValue::Undirected;

    if (!builtins.isCallTo(expression, BuiltinId::Negate)
        || expression.asCall().arguments.size() != 1)
        return InfiniteValue::None;
    const expression::Expr& operand = expression.asCall().arguments.front();
    if (isPredefinedAtom(
            operand, symbols::PredefinedSymbolId::Infinity, symbolRegistry))
        return InfiniteValue::Negative;
    if (isPredefinedAtom(
            operand, symbols::PredefinedSymbolId::ComplexInfinity, symbolRegistry))
        return InfiniteValue::Undirected;
    return InfiniteValue::None;
}

[[nodiscard]] bool isExactlyZero(const expression::Expr& value) noexcept {
    if (value.isNumber())
        return value.asNumber().isZero();
    if (value.isDecimalApproximation())
        return value.asDecimalApproximation().informationExactlyZero();
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        return complex.realInformationExactlyZero()
            && complex.imaginaryInformationExactlyZero();
    }
    return false;
}

[[nodiscard]] bool isExactlyNonZero(const expression::Expr& value) noexcept {
    if (value.isNumber())
        return !value.asNumber().isZero();
    if (value.isDecimalApproximation()) {
        const auto& decimal = value.asDecimalApproximation();
        return decimal.informationUpper() < numeric::Rational{}
            || decimal.informationLower() > numeric::Rational{};
    }
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        const bool realExcludesZero = complex.realInformationUpper() < numeric::Rational{}
            || complex.realInformationLower() > numeric::Rational{};
        const bool imaginaryExcludesZero = complex.imaginaryInformationUpper() < numeric::Rational{}
            || complex.imaginaryInformationLower() > numeric::Rational{};
        return realExcludesZero || imaginaryExcludesZero;
    }
    return false;
}

[[nodiscard]] bool hasExactUnitMagnitude(const expression::Expr& value) {
    if (!value.isNumber())
        return false;
    const numeric::RealNumber real = value.asNumber().realPart();
    const numeric::RealNumber imaginary = value.asNumber().imaginaryPart();
    return real * real + imaginary * imaginary
        == numeric::RealNumber{numeric::BigInt{1}};
}

[[nodiscard]] bool propagatesIndeterminate(
    const BuiltinDefinition& definition,
    const mathematics::MathRegistry& mathRegistry) noexcept {
    if (mathRegistry.findFunction(definition.symbol))
        return true;
    switch (definition.id) {
    case BuiltinId::Add:
    case BuiltinId::Subtract:
    case BuiltinId::Multiply:
    case BuiltinId::Divide:
    case BuiltinId::Power:
    case BuiltinId::Negate:
    case BuiltinId::Factorial:
    case BuiltinId::Floor:
    case BuiltinId::Ceil:
    case BuiltinId::Trunc:
    case BuiltinId::Round:
    case BuiltinId::Frac:
    case BuiltinId::Fma:
    case BuiltinId::Clamp:
    case BuiltinId::Proj:
        return true;
    default:
        return false;
    }
}

[[nodiscard]] std::optional<expression::Expr> exceptionalArithmeticResult(
    const BuiltinDefinition& definition,
    std::span<const expression::Expr> arguments,
    const BuiltinRegistry& builtins,
    const symbols::SymbolRegistry& symbolRegistry,
    const mathematics::MathRegistry& mathRegistry) {
    const bool hasIndeterminate = std::any_of(
        arguments.begin(), arguments.end(),
        [&](const expression::Expr& argument) {
            return isIndeterminate(argument, symbolRegistry);
        });
    if (hasIndeterminate) {
        if (definition.id == BuiltinId::Equal)
            return expression::Expr{false};
        if (definition.id == BuiltinId::NotEqual)
            return expression::Expr{true};
        if (definition.id == BuiltinId::Element)
            return expression::Expr{false};
        if (propagatesIndeterminate(definition, mathRegistry))
            return predefinedValue(
                symbols::PredefinedSymbolId::Indeterminate, symbolRegistry);
    }

    const auto indeterminate = [&] {
        return predefinedValue(
            symbols::PredefinedSymbolId::Indeterminate, symbolRegistry);
    };
    const auto complexInfinity = [&] {
        return predefinedValue(
            symbols::PredefinedSymbolId::ComplexInfinity, symbolRegistry);
    };

    switch (definition.id) {
    case BuiltinId::Add: {
        bool positiveInfinity = false;
        bool negativeInfinity = false;
        for (const expression::Expr& argument : arguments) {
            const InfiniteValue infinity = infiniteValue(
                argument, builtins, symbolRegistry);
            positiveInfinity = positiveInfinity || infinity == InfiniteValue::Positive;
            negativeInfinity = negativeInfinity || infinity == InfiniteValue::Negative;
        }
        if (positiveInfinity && negativeInfinity)
            return indeterminate();
        break;
    }
    case BuiltinId::Subtract:
        if (arguments.size() == 2) {
            const InfiniteValue lhs = infiniteValue(
                arguments[0], builtins, symbolRegistry);
            const InfiniteValue rhs = infiniteValue(
                arguments[1], builtins, symbolRegistry);
            if (lhs != InfiniteValue::None && rhs != InfiniteValue::None
                && (lhs == rhs || lhs == InfiniteValue::Undirected
                    || rhs == InfiniteValue::Undirected))
                return indeterminate();
        }
        break;
    case BuiltinId::Multiply: {
        const bool hasZero = std::any_of(
            arguments.begin(), arguments.end(), isExactlyZero);
        const bool hasInfinity = std::any_of(
            arguments.begin(), arguments.end(),
            [&](const expression::Expr& argument) {
                return infiniteValue(argument, builtins, symbolRegistry)
                    != InfiniteValue::None;
            });
        if (hasZero && hasInfinity)
            return indeterminate();
        break;
    }
    case BuiltinId::Divide:
        if (arguments.size() == 2) {
            const InfiniteValue numeratorInfinity =
                infiniteValue(arguments[0], builtins, symbolRegistry);
            const InfiniteValue denominatorInfinity =
                infiniteValue(arguments[1], builtins, symbolRegistry);
            if (numeratorInfinity != InfiniteValue::None
                && denominatorInfinity != InfiniteValue::None)
                return indeterminate();
            if (isExactlyZero(arguments[1])) {
                if (arguments[0].isArray()) {
                    const auto& array = arguments[0].asArray();
                    std::vector<expression::Expr> elements;
                    elements.reserve(array.size());
                    for (std::size_t i = 0; i < array.size(); ++i) {
                        const std::array<expression::Expr, 2> pair{
                            array.element(i), arguments[1]};
                        const auto element = exceptionalArithmeticResult(
                            definition, pair, builtins, symbolRegistry, mathRegistry);
                        elements.push_back(element ? *element : indeterminate());
                    }
                    return expression::Expr::array(array.shape, std::move(elements));
                }
                if (isExactlyZero(arguments[0]))
                    return indeterminate();
                if (isExactlyNonZero(arguments[0])
                    || numeratorInfinity != InfiniteValue::None) {
                    return complexInfinity();
                }

                const mathematics::ValueFacts facts = mathematics::inferValueFacts(
                    arguments[0], builtins, mathRegistry);
                if (facts.sign == mathematics::RealSign::Positive
                    || facts.sign == mathematics::RealSign::Negative
                    || facts.sign == mathematics::RealSign::NonZero)
                    return complexInfinity();
                return indeterminate();
            }
        }
        break;
    case BuiltinId::Power:
        if (arguments.size() == 2) {
            const expression::Expr& base = arguments[0];
            const expression::Expr& exponent = arguments[1];
            const InfiniteValue baseInfinity = infiniteValue(
                base, builtins, symbolRegistry);
            const InfiniteValue exponentInfinity = infiniteValue(
                exponent, builtins, symbolRegistry);

            if (exponentInfinity == InfiniteValue::Undirected)
                return indeterminate();
            if (baseInfinity != InfiniteValue::None && isExactlyZero(exponent))
                return indeterminate();
            if (isExactlyZero(base) && isExactlyZero(exponent))
                return indeterminate();
            if (isExactlyZero(base) && exponent.isNumber()) {
                const numeric::Number& number = exponent.asNumber();
                const numeric::RealNumber realPart = number.realPart();
                if (number.isComplex() && realPart.isZero())
                    return indeterminate();
                if (realPart.isNegative())
                    return complexInfinity();
            }
            if (exponentInfinity != InfiniteValue::None
                && hasExactUnitMagnitude(base))
                return indeterminate();
        }
        break;
    case BuiltinId::Negate:
        if (arguments.size() == 1
            && infiniteValue(arguments.front(), builtins, symbolRegistry)
                == InfiniteValue::Undirected)
            return complexInfinity();
        break;
    default:
        break;
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<mathematics::NumericDomain> solveDomainSymbol(
    const expression::Expr& expression,
    const symbols::SymbolRegistry& symbols) {
    if (!expression.isSymbol())
        return std::nullopt;
    const auto* definition = symbols.find(expression.asSymbol());
    if (!definition || definition->kind != symbols::PredefinedSymbolKind::MathematicalDomain)
        return std::nullopt;
    switch (definition->id) {
    case symbols::PredefinedSymbolId::IntegerDomain: return mathematics::NumericDomain::Integer;
    case symbols::PredefinedSymbolId::RationalDomain: return mathematics::NumericDomain::Rational;
    case symbols::PredefinedSymbolId::RealDomain: return mathematics::NumericDomain::Real;
    case symbols::PredefinedSymbolId::ComplexDomain: return mathematics::NumericDomain::Complex;
    default: return std::nullopt;
    }
}

void collectSolveUnknowns(
    const expression::Expr& expression,
    const symbols::SymbolRegistry& symbols,
    const BuiltinRegistry& builtins,
    std::vector<expression::Symbol>& result) {
    if (expression.isSymbol()) {
        const expression::Symbol& symbol = expression.asSymbol();
        if (symbols.contains(symbol) || builtins.contains(symbol))
            return;
        if (std::find(result.begin(), result.end(), symbol) == result.end())
            result.push_back(symbol);
        return;
    }
    if (expression.isCall()) {
        for (const expression::Expr& argument : expression.asCall().arguments)
            collectSolveUnknowns(argument, symbols, builtins, result);
        return;
    }
    if (expression.isArray()) {
        for (const expression::Expr& element : expression.asArray().storedExpressions())
            collectSolveUnknowns(element, symbols, builtins, result);
        return;
    }
    if (expression.isList())
        for (const expression::Expr& element : expression.asList().elements)
            collectSolveUnknowns(element, symbols, builtins, result);
}

void validateSolveVariable(
    const expression::Symbol& variable,
    const symbols::SymbolRegistry& symbols,
    const BuiltinRegistry& builtins) {
    if (symbols.contains(variable) || builtins.contains(variable))
        error::throwCalcError(
            error::CalcErrorType::Type,
            "solve variable must be an unprotected user symbol");
}


[[nodiscard]] bool containsBuiltinCall(
    const expression::Expr& root,
    const expression::Symbol& head) {
    std::vector<expression::Expr> pending{root};
    while (!pending.empty()) {
        expression::Expr current = std::move(pending.back());
        pending.pop_back();
        if (current.isCall()) {
            const auto& call = current.asCall();
            if (call.head.sameIdentity(head))
                return true;
            for (const expression::Expr& argument : call.arguments)
                pending.push_back(argument);
        }
        else if (current.isArray()) {
            for (const expression::Expr& element : current.asArray().storedExpressions())
                pending.push_back(element);
        }
        else if (current.isList()) {
            for (const expression::Expr& element : current.asList().elements)
                pending.push_back(element);
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

[[nodiscard]] std::vector<expression::Symbol> collectVariables(
    const expression::Expr& specification) {
    std::vector<expression::Symbol> variables;
    if (specification.isSymbol()) {
        variables.push_back(specification.asSymbol());
        return variables;
    }
    if (specification.isArray() && specification.asArray().rank() == 1) {
        const auto& array = specification.asArray();
        variables.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i) {
            const expression::Expr item = array.element(i);
            if (!item.isSymbol())
                error::throwCalcError(
                    error::CalcErrorType::Type,
                    "collect variable array must contain only symbols");
            variables.push_back(item.asSymbol());
        }
        return variables;
    }
    error::throwCalcError(
        error::CalcErrorType::Type,
        "collect expects an expression and a symbol or symbol array");
}

void consumeSolverSolutionBranches(const solver::SolutionSet& solutions) {
    std::size_t branches = solutions.branches().size();
    for (const solver::SolutionCase& item : solutions.cases())
        branches += std::max<std::size_t>(1, item.branches.size());
    consumeEvaluationBudget(
        EvaluationResource::SolverBranch,
        std::max<std::size_t>(1, branches));
}

[[nodiscard]] std::optional<std::pair<expression::Expr, expression::Expr>>
plotRangePair(const expression::Expr& expression) {
    if (expression.isArray()) {
        const auto& array = expression.asArray();
        if (array.rank() != 1 || array.size() != 2)
            return std::nullopt;
        return std::pair{array.element(0), array.element(1)};
    }
    if (expression.isList() && expression.asList().elements.size() == 2)
        return std::pair{
            expression.asList().elements[0], expression.asList().elements[1]};
    return std::nullopt;
}

[[nodiscard]] std::optional<bool> plotBooleanOptionValue(
    const expression::Expr& expression) noexcept {
    if (!expression.isBoolean())
        return std::nullopt;
    return expression.asBoolean();
}

[[nodiscard]] plot::PlotRequestOptions parsePlotOptions(
    std::span<const expression::Expr> arguments,
    const BuiltinRegistry& builtins) {
    plot::PlotRequestOptions options;
    for (std::size_t i = 2; i < arguments.size(); ++i) {
        const auto& option = arguments[i];
        if (!builtins.isCallTo(option, BuiltinId::Rule)
            || option.asCall().arguments.size() != 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Plot options must be rules such as PlotRange -> {min,max}");

        const auto& ruleArguments = option.asCall().arguments;
        if (!ruleArguments[0].isSymbol())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Plot option names must be symbols");

        const std::string_view name = ruleArguments[0].asSymbol().view();
        if (name == "PlotRange") {
            if (options.plotRange)
                error::throwCalcError(
                    error::CalcErrorType::Domain,
                    "Duplicate option 'PlotRange' for Plot");
            const auto pair = plotRangePair(ruleArguments[1]);
            if (!pair)
                error::throwCalcError(
                    error::CalcErrorType::Type,
                    "PlotRange expects {minimum,maximum}");
            options.plotRange = plot::PlotRangeOption{pair->first, pair->second};
        }
        else if (name == "AspectRatio") {
            if (options.aspectRatio)
                error::throwCalcError(
                    error::CalcErrorType::Domain,
                    "Duplicate option 'AspectRatio' for Plot");
            options.aspectRatio = ruleArguments[1];
        }
        else if (name == "Ticks") {
            if (options.ticks)
                error::throwCalcError(
                    error::CalcErrorType::Domain,
                    "Duplicate option 'Ticks' for Plot");
            const auto value = plotBooleanOptionValue(ruleArguments[1]);
            if (!value)
                error::throwCalcError(
                    error::CalcErrorType::Type,
                    "Ticks expects True or False");
            options.ticks = *value;
        }
        else if (name == "PlotPoints") {
            if (options.plotPoints)
                error::throwCalcError(
                    error::CalcErrorType::Domain,
                    "Duplicate option 'PlotPoints' for Plot");
            const auto& value = ruleArguments[1];
            if (!value.isNumber() || !value.asNumber().asReal().isInteger())
                error::throwCalcError(
                    error::CalcErrorType::Type,
                    "PlotPoints expects an integer from 100 through 1024");
            const auto parsed = numeric::tryToUint64(value.asNumber().asReal().asInteger());
            if (!parsed || *parsed < 100 || *parsed > 1024)
                error::throwCalcError(
                    error::CalcErrorType::Domain,
                    "PlotPoints expects an integer from 100 through 1024");
            options.plotPoints = static_cast<std::size_t>(*parsed);
        }
        else
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "Unknown option '" + std::string{name} + "' for Plot");
    }
    return options;
}

[[nodiscard]] std::optional<plot::PlotRequestSet> parsePlotRequest(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins) {
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::Plot)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() < 2)
        return std::nullopt;

    const auto iterator = parseRangeIteratorSpec(arguments[1]);
    if (!iterator)
        return std::nullopt;

    std::vector<expression::Expr> curves;
    if (arguments[0].isArray()) {
        const auto& array = arguments[0].asArray();
        if (array.rank() != 1 || array.size() == 0)
            return std::nullopt;
        curves = array.materialize();
    }
    else if (arguments[0].isList()) {
        if (arguments[0].asList().elements.empty())
            return std::nullopt;
        curves = arguments[0].asList().elements;
    }
    else
        curves.push_back(arguments[0]);

    return plot::PlotRequestSet{
        std::move(curves), iterator->variable, iterator->lower, iterator->upper,
        parsePlotOptions(arguments, builtins)};
}

[[nodiscard]] bool appendPlotRequests(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins,
    std::vector<plot::PlotRequest>& requests) {
    if (const auto requestSet = parsePlotRequest(expression, builtins)) {
        auto split = plot::splitPlotRequests(*requestSet);
        requests.insert(requests.end(), split.begin(), split.end());
        return true;
    }

    if (!expression.isCall())
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::Show)
        return false;

    for (const auto& argument : expression.asCall().arguments)
        if (!appendPlotRequests(argument, builtins, requests))
            return false;
    return !requests.empty();
}

[[nodiscard]] std::optional<std::vector<plot::PlotRequest>> parsePlotRequests(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins) {
    std::vector<plot::PlotRequest> requests;
    if (!appendPlotRequests(expression, builtins, requests) || requests.empty())
        return std::nullopt;
    return requests;
}


[[nodiscard]] std::optional<std::pair<expression::Expr, expression::Expr>>
parametricCoordinatePair(const expression::Expr& expression) {
    if (expression.isArray()) {
        const auto& array = expression.asArray();
        if (array.rank() != 1 || array.size() != 2)
            return std::nullopt;
        return std::pair{array.element(0), array.element(1)};
    }
    if (expression.isList() && expression.asList().elements.size() == 2)
        return std::pair{expression.asList().elements[0], expression.asList().elements[1]};
    return std::nullopt;
}

[[nodiscard]] std::optional<plot::ParametricPlotRequestSet> parseParametricPlotRequest(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins) {
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::ParametricPlot)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() < 2)
        return std::nullopt;
    const auto iterator = parseRangeIteratorSpec(arguments[1]);
    if (!iterator)
        return std::nullopt;

    std::vector<std::pair<expression::Expr, expression::Expr>> curves;
    if (arguments[0].isArray() && arguments[0].asArray().rank() == 2) {
        const auto& array = arguments[0].asArray();
        if (array.extent(0) == 0 || array.extent(1) != 2)
            return std::nullopt;
        curves.reserve(array.extent(0));
        for (std::size_t row = 0; row < array.extent(0); ++row)
            curves.emplace_back(array.element(row * 2), array.element(row * 2 + 1));
    }
    else {
        std::vector<expression::Expr> entries;
        if (arguments[0].isList())
            entries = arguments[0].asList().elements;
        bool allEntriesArePairs = !entries.empty();
        if (allEntriesArePairs)
            for (const auto& entry : entries)
                if (!parametricCoordinatePair(entry)) {
                    allEntriesArePairs = false;
                    break;
                }
        if (allEntriesArePairs) {
            for (const auto& entry : entries)
                curves.push_back(*parametricCoordinatePair(entry));
        }
        else if (const auto single = parametricCoordinatePair(arguments[0]))
            curves.push_back(*single);
        else
            return std::nullopt;
    }

    auto options = parsePlotOptions(arguments, builtins);
    if (options.plotRange)
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "ParametricPlot does not yet support PlotRange; use Automatic range");
    return plot::ParametricPlotRequestSet{
        std::move(curves), iterator->variable, iterator->lower, iterator->upper,
        std::move(options)};
}


struct ParsedListPlot final {
    std::vector<std::pair<expression::Expr, expression::Expr>> points;
    bool joined = false;
    bool joinedSpecified = false;
    std::optional<expression::Expr> aspectRatio;
    std::optional<bool> ticks;
};

[[nodiscard]] expression::Expr exactZeroExpr() {
    return expression::Expr{numeric::Number{numeric::BigInt{0}}};
}

[[nodiscard]] std::optional<ParsedListPlot> parseListPlotRequest(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins) {
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::ListPlot)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.empty())
        return std::nullopt;

    ParsedListPlot request;
    const auto zero = exactZeroExpr();
    const auto appendRows = [&](const std::vector<std::vector<expression::Expr>>& rows)
        -> bool {
        if (rows.empty())
            return false;
        const std::size_t width = rows.front().size();
        if (width != 1 && width != 2)
            return false;
        for (const auto& row : rows) {
            if (row.size() != width)
                return false;
            request.points.emplace_back(row[0], width == 1 ? zero : row[1]);
        }
        return true;
    };

    const auto& data = arguments[0];
    if (data.isArray()) {
        const auto& array = data.asArray();
        if (array.rank() == 1) {
            if (array.size() == 0)
                return std::nullopt;
            request.points.reserve(array.size());
            for (std::size_t i = 0; i < array.size(); ++i)
                request.points.emplace_back(array.element(i), zero);
        }
        else if (array.rank() == 2) {
            const std::size_t rows = array.extent(0);
            const std::size_t columns = array.extent(1);
            if (rows == 0 || (columns != 1 && columns != 2))
                return std::nullopt;
            request.points.reserve(rows);
            for (std::size_t row = 0; row < rows; ++row)
                request.points.emplace_back(
                    array.element(row * columns),
                    columns == 1 ? zero : array.element(row * columns + 1));
        }
        else
            return std::nullopt;
    }
    else if (data.isList()) {
        const auto& elements = data.asList().elements;
        if (elements.empty())
            return std::nullopt;
        bool nested = false;
        for (const auto& element : elements)
            nested = nested || element.isArray() || element.isList();
        if (!nested) {
            request.points.reserve(elements.size());
            for (const auto& element : elements)
                request.points.emplace_back(element, zero);
        }
        else {
            std::vector<std::vector<expression::Expr>> rows;
            rows.reserve(elements.size());
            for (const auto& element : elements) {
                if (element.isArray()) {
                    const auto& row = element.asArray();
                    if (row.rank() != 1)
                        return std::nullopt;
                    rows.push_back(row.materialize());
                }
                else if (element.isList())
                    rows.push_back(element.asList().elements);
                else
                    return std::nullopt;
            }
            if (!appendRows(rows))
                return std::nullopt;
        }
    }
    else
        return std::nullopt;

    for (std::size_t i = 1; i < arguments.size(); ++i) {
        const auto& option = arguments[i];
        if (!builtins.isCallTo(option, BuiltinId::Rule)
            || option.asCall().arguments.size() != 2
            || !option.asCall().arguments[0].isSymbol())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "ListPlot options must be rules such as Joined -> True");
        const auto& rule = option.asCall().arguments;
        const std::string_view name = rule[0].asSymbol().view();
        if (name == "Joined") {
            if (request.joinedSpecified)
                error::throwCalcError(error::CalcErrorType::Domain,
                    "Duplicate option 'Joined' for ListPlot");
            const auto value = plotBooleanOptionValue(rule[1]);
            if (!value)
                error::throwCalcError(error::CalcErrorType::Type,
                    "Joined expects True or False");
            request.joined = *value;
            request.joinedSpecified = true;
        }
        else if (name == "AspectRatio") {
            if (request.aspectRatio)
                error::throwCalcError(error::CalcErrorType::Domain,
                    "Duplicate option 'AspectRatio' for ListPlot");
            request.aspectRatio = rule[1];
        }
        else if (name == "Ticks") {
            if (request.ticks)
                error::throwCalcError(error::CalcErrorType::Domain,
                    "Duplicate option 'Ticks' for ListPlot");
            const auto value = plotBooleanOptionValue(rule[1]);
            if (!value)
                error::throwCalcError(error::CalcErrorType::Type,
                    "Ticks expects True or False");
            request.ticks = *value;
        }
        else if (name == "PlotPoints")
            error::throwCalcError(error::CalcErrorType::Domain,
                "ListPlot does not support PlotPoints because input points are explicit");
        else if (name == "PlotRange")
            error::throwCalcError(error::CalcErrorType::Domain,
                "ListPlot does not yet support PlotRange");
        else
            error::throwCalcError(error::CalcErrorType::Domain,
                "Unknown option '" + std::string{name} + "' for ListPlot");
    }
    return request;
}

[[nodiscard]] expression::Expr listPlotNormalCoordinates(const ParsedListPlot& request) {
    std::vector<expression::Expr> values;
    values.reserve(request.points.size() * 2);
    for (const auto& [x, y] : request.points) {
        values.push_back(x);
        values.push_back(y);
    }
    return expression::Expr::array({request.points.size(), 2}, std::move(values));
}

[[nodiscard]] std::optional<numeric::BigFloat> listPlotCoordinateValue(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::size_t precisionBits = 64) {
    const expression::Symbol dummy{"__mmcal_listplot_parameter"};
    const auto compiled = plot::compilePlotProgram(expression, dummy, builtins, mathematics);
    if (!compiled || !compiled.program)
        return std::nullopt;
    plot::BigFloatPlotExecutor executor(*compiled.program, precisionBits, angles);
    const auto zero = numeric::BigFloat::fromBigInt(
        numeric::BigInt{0}, precisionBits, numeric::RoundingMode::NearestEven);
    const auto value = executor.evaluate(zero);
    if (!value.finite())
        return std::nullopt;
    return value.value;
}

[[nodiscard]] std::optional<double> listPlotPositiveDouble(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto value = listPlotCoordinateValue(expression, builtins, mathematics, angles);
    if (!value || !value->isPositive())
        return std::nullopt;
    const auto rational = value->toRational();
    try {
        const long double numerator = std::stold(rational.numerator().toString());
        const long double denominator = std::stold(rational.denominator().toString());
        const double result = static_cast<double>(numerator / denominator);
        if (!std::isfinite(result) || !(result > 0.0))
            return std::nullopt;
        return result;
    }
    catch (...) {
        return std::nullopt;
    }
}

struct PreparedListPlotGraphics final {
    plot::PlotScene plotScene;
    plot::CurveViewRange2D range;
    plot::PlotViewportMm viewport;
    std::size_t precisionBits = 64;
    bool ticks = true;
};

[[nodiscard]] std::optional<PreparedListPlotGraphics> prepareListPlotGraphics(
    const ParsedListPlot& request,
    const BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    constexpr std::size_t precisionBits = 64;
    std::vector<std::pair<numeric::BigFloat, numeric::BigFloat>> points;
    points.reserve(request.points.size());
    for (const auto& [xExpression, yExpression] : request.points) {
        auto x = listPlotCoordinateValue(xExpression, builtins, mathematics, angles, precisionBits);
        auto y = listPlotCoordinateValue(yExpression, builtins, mathematics, angles, precisionBits);
        if (!x || !y)
            return std::nullopt;
        points.emplace_back(std::move(*x), std::move(*y));
    }

    const auto zeroExpression = exactZeroExpr();
    plot::SampledCurve2D sampled;
    sampled.precisionBits = precisionBits;
    plot::SampledCurveSegment2D sampledSegment{
        plot::PlotInterval{zeroExpression, plot::PlotEndpointInclusion::Closed,
            zeroExpression, plot::PlotEndpointInclusion::Closed},
        {}, plot::PlotSegmentGeometryKind::Polyline,
        std::nullopt, std::nullopt, std::nullopt, std::nullopt,
        false, false, std::nullopt, std::nullopt,
        std::nullopt, std::nullopt, false, false};
    sampledSegment.samples.reserve(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto parameter = numeric::BigFloat::fromBigInt(
            numeric::BigInt{static_cast<std::int64_t>(i)}, precisionBits,
            numeric::RoundingMode::NearestEven);
        sampledSegment.samples.emplace_back(
            parameter, points[i].first, points[i].second,
            plot::PlotNumericStatus::Finite, plot::PlotPointTag::Coarse);
    }
    sampled.segments.push_back(std::move(sampledSegment));

    const auto xRange = plot::estimateCurveCoordinateRange(
        sampled, plot::CurveCoordinate2D::X);
    const auto yRange = plot::estimateCurveCoordinateRange(
        sampled, plot::CurveCoordinate2D::Y);
    if (!xRange || !yRange)
        return std::nullopt;

    plot::PlotViewportMm viewport;
    if (request.aspectRatio) {
        const auto ratio = listPlotPositiveDouble(
            *request.aspectRatio, builtins, mathematics, angles);
        if (!ratio)
            return std::nullopt;
        viewport.heightMm = viewport.widthMm * *ratio;
    }
    plot::PlotScene scene;
    plot::PlotSceneCurve curve;
    curve.id = 1;
    if (request.joined && points.size() >= 2) {
        plot::PlotSceneSegment segment{
            plot::PlotInterval{zeroExpression, plot::PlotEndpointInclusion::Closed,
                zeroExpression, plot::PlotEndpointInclusion::Closed},
            plot::PlotSegmentGeometryKind::Polyline, {},
            std::nullopt, std::nullopt, std::nullopt, std::nullopt,
            false, false, std::nullopt, std::nullopt};
        segment.vertices.reserve(points.size());
        for (const auto& [x, y] : points)
            segment.vertices.push_back(plot::PlotSceneVertex{x, y, std::nullopt});
        curve.segments.push_back(std::move(segment));
    }
    for (const auto& [x, y] : points) {
        plot::PlotSceneSegment marker{
            plot::PlotInterval{zeroExpression, plot::PlotEndpointInclusion::Closed,
                zeroExpression, plot::PlotEndpointInclusion::Closed},
            plot::PlotSegmentGeometryKind::Polyline, {},
            std::nullopt, std::nullopt, std::nullopt, std::nullopt,
            false, false, std::nullopt, std::nullopt};
        marker.vertices.push_back(plot::PlotSceneVertex{x, y, std::nullopt});
        marker.markLowerEndpoint = true;
        curve.segments.push_back(std::move(marker));
    }
    scene.curves.push_back(std::move(curve));

    const bool ticks = request.ticks.value_or(true);
    return PreparedListPlotGraphics{
        std::move(scene), plot::CurveViewRange2D{*xRange, *yRange},
        viewport, precisionBits, ticks};
}

[[nodiscard]] std::optional<graphics::GraphicsScene> buildListPlotGraphics(
    const ParsedListPlot& request,
    const BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto prepared = prepareListPlotGraphics(request, builtins, mathematics, angles);
    if (!prepared)
        return std::nullopt;
    const auto transform = plot::makePlotViewTransform(
        prepared->range, prepared->precisionBits, prepared->viewport);
    if (!transform)
        return std::nullopt;

    plot::PlotAxisLayoutOptions axisOptions;
    axisOptions.generateMajorTicks = prepared->ticks;
    const auto axes = plot::layoutPlotAxes(*transform, axisOptions);
    if (!axes || !axes.layout)
        return std::nullopt;
    plot::PlotGraphicsLoweringOptions loweringOptions;
    loweringOptions.showMajorTicks = prepared->ticks;
    loweringOptions.showMajorTickLabels = prepared->ticks;
    const auto lowered = plot::lowerPlotToGraphics(
        prepared->plotScene, *axes.layout, *transform, loweringOptions);
    if (!lowered || !lowered.scene)
        return std::nullopt;
    return *lowered.scene;
}

enum class CompositeFirstDrawable {
    Plot,
    ListPlot
};

struct ParsedCompositeDrawable final {
    std::vector<plot::PlotRequest> plots;
    std::vector<ParsedListPlot> listPlots;
    std::optional<CompositeFirstDrawable> first;
};

[[nodiscard]] bool appendCompositeDrawable(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins,
    ParsedCompositeDrawable& result) {
    if (const auto request = parsePlotRequest(expression, builtins)) {
        if (!result.first)
            result.first = CompositeFirstDrawable::Plot;
        auto split = plot::splitPlotRequests(*request);
        result.plots.insert(result.plots.end(), split.begin(), split.end());
        return true;
    }
    if (const auto request = parseListPlotRequest(expression, builtins)) {
        if (!result.first)
            result.first = CompositeFirstDrawable::ListPlot;
        result.listPlots.push_back(*request);
        return true;
    }
    if (!expression.isCall())
        return false;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::Show
        || expression.asCall().arguments.empty())
        return false;
    for (const auto& argument : expression.asCall().arguments)
        if (!appendCompositeDrawable(argument, builtins, result))
            return false;
    return true;
}

[[nodiscard]] std::optional<ParsedCompositeDrawable> parseCompositeDrawable(
    const expression::Expr& expression,
    const BuiltinRegistry& builtins) {
    ParsedCompositeDrawable result;
    if (!appendCompositeDrawable(expression, builtins, result)
        || (!result.first)
        || (result.plots.empty() && result.listPlots.empty()))
        return std::nullopt;
    return result;
}

[[nodiscard]] plot::CurveRangeEstimate1D transformAxisRange(
    const numeric::BigFloat& minimum,
    const numeric::BigFloat& maximum) {
    return plot::CurveRangeEstimate1D{
        minimum, maximum, minimum, maximum, minimum, maximum,
        0, 0, false, false};
}

[[nodiscard]] plot::CurveViewRange2D unionCompositeRanges(
    const std::vector<plot::CurveViewRange2D>& ranges) {
    plot::CurveViewRange2D result = ranges.front();
    const auto merge = [](plot::CurveRangeEstimate1D& destination,
                           const plot::CurveRangeEstimate1D& source) {
        if (source.observedMinimum < destination.observedMinimum)
            destination.observedMinimum = source.observedMinimum;
        if (source.observedMaximum > destination.observedMaximum)
            destination.observedMaximum = source.observedMaximum;
        if (source.dataMinimum < destination.dataMinimum)
            destination.dataMinimum = source.dataMinimum;
        if (source.dataMaximum > destination.dataMaximum)
            destination.dataMaximum = source.dataMaximum;
        if (source.viewMinimum < destination.viewMinimum)
            destination.viewMinimum = source.viewMinimum;
        if (source.viewMaximum > destination.viewMaximum)
            destination.viewMaximum = source.viewMaximum;
        destination.includedFiniteSamples += source.includedFiniteSamples;
        destination.excludedFiniteSamples += source.excludedFiniteSamples;
        destination.boundaryTrimmed =
            destination.boundaryTrimmed || source.boundaryTrimmed;
        destination.flatExpanded = destination.flatExpanded || source.flatExpanded;
    };
    for (std::size_t i = 1; i < ranges.size(); ++i) {
        merge(result.x, ranges[i].x);
        merge(result.y, ranges[i].y);
    }
    return result;
}

struct CompositeGraphicsResult final {
    plot::PlotPipelineStatus status = plot::PlotPipelineStatus::EmptyRequest;
    std::optional<graphics::GraphicsScene> scene;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == plot::PlotPipelineStatus::Success && scene.has_value();
    }
};

[[nodiscard]] CompositeGraphicsResult buildCompositeGraphics(
    const ParsedCompositeDrawable& request,
    const BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol) {
    std::optional<plot::PlotPipelineOutput> plotOutput;
    if (!request.plots.empty()) {
        auto built = plot::buildPlotPipeline(
            request.plots, builtins, mathematics, angles, infinitySymbol);
        if (!built || !built.output)
            return CompositeGraphicsResult{built.status, std::nullopt};
        plotOutput = std::move(*built.output);
    }

    std::vector<PreparedListPlotGraphics> listOutputs;
    listOutputs.reserve(request.listPlots.size());
    for (const ParsedListPlot& listPlot : request.listPlots) {
        auto prepared = prepareListPlotGraphics(listPlot, builtins, mathematics, angles);
        if (!prepared)
            return CompositeGraphicsResult{
                plot::PlotPipelineStatus::SamplingFailed, std::nullopt};
        listOutputs.push_back(std::move(*prepared));
    }

    std::vector<plot::CurveViewRange2D> ranges;
    std::size_t precisionBits = 64;
    plot::PlotViewportMm viewport;
    bool ticks = true;
    plot::PlotScene combinedScene;
    plot::PlotSceneCurveId maximumCurveId = 0;
    if (plotOutput) {
        precisionBits = plotOutput->transform.precisionBits();
        ranges.push_back(plot::CurveViewRange2D{
            transformAxisRange(
                plotOutput->transform.xMinimum(), plotOutput->transform.xMaximum()),
            transformAxisRange(
                plotOutput->transform.yMinimum(), plotOutput->transform.yMaximum())});
        viewport = plotOutput->transform.viewport();
        ticks = !plotOutput->axes.xAxis.majorTicks.empty()
            || !plotOutput->axes.yAxis.majorTicks.empty();
        combinedScene = plotOutput->plotScene;
        for (const auto& curve : combinedScene.curves)
            maximumCurveId = std::max(maximumCurveId, curve.id);
    }
    for (PreparedListPlotGraphics& listOutput : listOutputs) {
        ranges.push_back(listOutput.range);
        if (!plotOutput && ranges.size() == 1) {
            precisionBits = listOutput.precisionBits;
            viewport = listOutput.viewport;
            ticks = listOutput.ticks;
        }
        for (auto& curve : listOutput.plotScene.curves) {
            curve.id = ++maximumCurveId;
            combinedScene.curves.push_back(std::move(curve));
        }
    }
    if (ranges.empty() || combinedScene.curves.empty())
        return CompositeGraphicsResult{
            plot::PlotPipelineStatus::SceneBuildFailed, std::nullopt};

    // Showの先頭がListPlotなら，先頭graphicsの物理縦横比とtick指定を継承する。
    if (request.first == CompositeFirstDrawable::ListPlot && !listOutputs.empty()) {
        viewport = listOutputs.front().viewport;
        ticks = listOutputs.front().ticks;
    }
    const plot::CurveViewRange2D range = unionCompositeRanges(ranges);
    const auto transform = plot::makePlotViewTransform(range, precisionBits, viewport);
    if (!transform)
        return CompositeGraphicsResult{
            plot::PlotPipelineStatus::ViewTransformFailed, std::nullopt};
    plot::PlotAxisLayoutOptions axisOptions;
    axisOptions.generateMajorTicks = ticks;
    const auto axes = plot::layoutPlotAxes(*transform, axisOptions);
    if (!axes || !axes.layout)
        return CompositeGraphicsResult{
            axes.status == plot::PlotAxisLayoutStatus::InvalidOptions
                ? plot::PlotPipelineStatus::InvalidOptions
                : plot::PlotPipelineStatus::AxisLayoutFailed,
            std::nullopt};
    plot::PlotGraphicsLoweringOptions loweringOptions;
    loweringOptions.showMajorTicks = ticks;
    loweringOptions.showMajorTickLabels = ticks;
    const auto lowered = plot::lowerPlotToGraphics(
        combinedScene, *axes.layout, *transform, loweringOptions);
    if (!lowered || !lowered.scene)
        return CompositeGraphicsResult{
            lowered.status == plot::PlotGraphicsLoweringStatus::InvalidOptions
                ? plot::PlotPipelineStatus::InvalidOptions
                : plot::PlotPipelineStatus::GraphicsLoweringFailed,
            std::nullopt};
    return CompositeGraphicsResult{
        plot::PlotPipelineStatus::Success, std::move(*lowered.scene)};
}

[[noreturn]] void throwPlotPipelineFailure(
    plot::PlotPipelineStatus status,
    std::string evaluationMessage) {
    if (status == plot::PlotPipelineStatus::InvalidOptions)
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "Plot options are invalid for the requested plot");
    error::throwCalcError(
        error::CalcErrorType::Evaluation, std::move(evaluationMessage));
}

struct ExportGraphicsTarget final {
    std::filesystem::path path;
    graphics::GraphicsFormat format = graphics::GraphicsFormat::Svg;
};

[[nodiscard]] std::optional<ExportGraphicsTarget> exportGraphicsTarget(
    const std::filesystem::path& path,
    const std::optional<std::string>& explicitFormat) {
    if (!explicitFormat) {
        const auto format = graphics::graphicsFormatFromExtension(path.extension().string());
        if (!format)
            return std::nullopt;
        return ExportGraphicsTarget{path, *format};
    }

    const auto format = graphics::parseGraphicsFormat(*explicitFormat);
    if (!format)
        return std::nullopt;

    std::filesystem::path outputPath = path;
    const auto extensionFormat = graphics::graphicsFormatFromExtension(path.extension().string());
    if (!extensionFormat || *extensionFormat != *format)
        outputPath += graphics::graphicsFormatExtension(*format);
    return ExportGraphicsTarget{std::move(outputPath), *format};
}

[[nodiscard]] std::size_t plotCoordinateDigits(std::size_t precisionBits) noexcept {
    // Plot内部の64-bit BigFloatは作業値であり，その全bitを利用者向け精度として主張しない。
    // 約8 guard bitを落として10進化するため，既定64 bitでは16桁になる。
    if (precisionBits <= 8)
        return 1;
    constexpr double decimalDigitsPerBit = 0.3010299956639812;
    return std::max<std::size_t>(
        1, static_cast<std::size_t>(
            static_cast<double>(precisionBits - 8) * decimalDigitsPerBit));
}

[[nodiscard]] numeric::DecimalApproximation plotCoordinateValue(
    const numeric::BigFloat& value,
    std::size_t significantDigits) {
    return numeric::DecimalApproximation::fromVerifiedValueSignificant(
        numeric::RealNumber{value.toRational()}, significantDigits);
}

[[nodiscard]] expression::Expr sampledSegmentCoordinates(
    const plot::SampledCurveSegment2D& segment,
    std::size_t precisionBits) {
    const std::size_t significantDigits = plotCoordinateDigits(precisionBits);
    std::vector<numeric::DecimalApproximation> values;
    values.reserve(segment.samples.size() * 2);
    std::size_t pointCount = 0;
    for (const auto& sample : segment.samples) {
        if (!sample.finite())
            continue;
        values.push_back(plotCoordinateValue(sample.x, significantDigits));
        values.push_back(plotCoordinateValue(sample.y, significantDigits));
        ++pointCount;
    }
    return expression::Expr::decimalArray({pointCount, 2}, std::move(values));
}

[[nodiscard]] expression::Expr sampledCurveCoordinates(
    const plot::SampledCurve2D& curve) {
    if (curve.segments.size() == 1)
        return sampledSegmentCoordinates(curve.segments.front(), curve.precisionBits);

    std::vector<expression::Expr> segments;
    segments.reserve(curve.segments.size());
    for (const auto& segment : curve.segments)
        segments.push_back(sampledSegmentCoordinates(segment, curve.precisionBits));
    return expression::Expr::list(std::move(segments));
}

[[nodiscard]] expression::Expr sampledCurvesCoordinates(
    const std::vector<plot::SampledCurve2D>& curves) {
    if (curves.size() == 1)
        return sampledCurveCoordinates(curves.front());

    std::vector<expression::Expr> result;
    result.reserve(curves.size());
    for (const auto& curve : curves)
        result.push_back(sampledCurveCoordinates(curve));
    return expression::Expr::list(std::move(result));
}

} // namespace


expression::Expr Evaluator::resolveHeldHistoryReferences(
    const expression::Expr& expression) {
    if (!expression.isCall())
        return expression;

    const expression::CallExpr& call = expression.asCall();
    if (const BuiltinDefinition* definition = registry_.find(call.head)) {
        if (definition->id == BuiltinId::History)
            return resolveHeldHistoryReferences(evaluateHistory(call.arguments));
        if (definition->id == BuiltinId::OutputHistory)
            return resolveHeldHistoryReferences(evaluateIndexedHistory(call.arguments, false));
        if (definition->id == BuiltinId::InputHistory) {
            const expression::Expr input = evaluateIndexedHistory(call.arguments, true);

            // In[n]は保存入力を現在の環境で再評価する。外側のHoldを解除して
            // 入力全体を評価すると他の項まで副作用を受けるため，履歴参照だけを
            // 独立Evaluatorで再評価し，PRNGとsession contextは現在のものを引き継ぐ。
            Evaluator nested{
                environment_, registry_, userFunctions_, symbolRegistry_, mathematics_, angleSemantics_};
            nested.setEvaluationLimits(limits_);
            nested.randomEngine_ = randomEngine_;
            expression::Expr evaluated = nested.evaluateMachine(input, nullptr, context_);
            randomEngine_ = nested.randomEngine_;
            return evaluated;
        }
    }

    std::vector<expression::Expr> arguments;
    arguments.reserve(call.arguments.size());
    bool changed = false;
    for (const expression::Expr& argument : call.arguments) {
        expression::Expr resolved = resolveHeldHistoryReferences(argument);
        changed = changed || !(resolved == argument);
        arguments.push_back(std::move(resolved));
    }
    return changed
        ? expression::Expr::rebuildCall(call, std::move(arguments))
        : expression;
}

expression::Expr Evaluator::materializeSafeHeldFrontends(
    const expression::Expr& expression,
    std::span<const expression::Symbol> protectedSymbols) {
    const std::size_t previousSize = materializationProtectedSymbols_.size();
    for (const expression::Symbol& symbol : protectedSymbols)
        if (std::find(materializationProtectedSymbols_.begin(),
                materializationProtectedSymbols_.end(), symbol)
            == materializationProtectedSymbols_.end())
            materializationProtectedSymbols_.push_back(symbol);

    try {
        expression::Expr result = materializeSafeHeldFrontends(expression);
        materializationProtectedSymbols_.resize(previousSize);
        return result;
    }
    catch (...) {
        materializationProtectedSymbols_.resize(previousSize);
        throw;
    }
}

expression::Expr Evaluator::materializeSafeHeldFrontends(
    const expression::Expr& expression) {
    // D / Series / 明示的なNormal[Series]等は，外側frontendがHoldしていても
    // 利用者が評価を要求した独立したsymbolic objectとしてmaterializeできる。
    if (isHeldSeriesPipeline(expression, registry_))
        return evaluate(expression);
    if (expression.isCall()) {
        if (const BuiltinDefinition* definition = registry_.find(expression.asCall().head); definition) {
            if (evaluation::isVectorCalculusBuiltin(definition->id))
                return evaluate(expression);

            // 内側のlimit/integrateは自身のbinder規則で変数を保護できる。
            // 外側D等のHoldだけを理由に未dispatchへ残すと，単体では計算できる式が
            // 合成した途端に未評価になるため，binder自身を一度だけ閉じてから外側へ渡す。
            if (definition->id == BuiltinId::SymbolicIntegral) {
                const auto& arguments = expression.asCall().arguments;
                // D[integrate[f,x],x]には微分器側の基本定理ruleがある。ここで原始函数を
                // 先に展開すると，同値ではあってもSimplifierへ不要な負荷を押し付けるため，
                // 外側binderと同じ変数の不定積分だけはheldのまま渡す。
                if (arguments.size() >= 2 && arguments[1].isSymbol()
                    && std::find(materializationProtectedSymbols_.begin(),
                            materializationProtectedSymbols_.end(), arguments[1].asSymbol())
                        != materializationProtectedSymbols_.end())
                    return expression;
                return evaluate(expression);
            }
            if (definition->id == BuiltinId::Limit)
                return evaluate(expression);
        }
    }

    if (expression.isArray()) {
        const auto& array = expression.asArray();
        std::vector<expression::Expr> elements;
        elements.reserve(array.size());
        bool changed = false;
        for (std::size_t i = 0; i < array.size(); ++i) {
            expression::Expr materialized = materializeSafeHeldFrontends(array.element(i));
            changed = changed || materialized != array.element(i);
            elements.push_back(std::move(materialized));
        }
        return changed ? expression::Expr::array(array.shape, std::move(elements)) : expression;
    }

    if (expression.isList()) {
        std::vector<expression::Expr> elements;
        elements.reserve(expression.asList().elements.size());
        bool changed = false;
        for (const expression::Expr& element : expression.asList().elements) {
            expression::Expr materialized = materializeSafeHeldFrontends(element);
            changed = changed || materialized != element;
            elements.push_back(std::move(materialized));
        }
        return changed ? expression::braceValue(std::move(elements)) : expression;
    }

    if (!expression.isCall())
        return expression;

    const expression::CallExpr& call = expression.asCall();
    const BuiltinDefinition* definition = registry_.find(call.head);
    // HoldAll/HoldFirst系の内部へ先回りすると，LimitやCases等が与える局所assumptionや
    // binder scopeを無視してしまう。通常評価で全引数を評価するbuiltinだけを降りる。
    if (!definition || definition->argumentEvaluation != ArgumentEvaluation::All)
        return expression;

    std::vector<expression::Expr> arguments;
    arguments.reserve(call.arguments.size());
    bool changed = false;
    for (const expression::Expr& argument : call.arguments) {
        expression::Expr materialized = materializeSafeHeldFrontends(argument);
        changed = changed || materialized != argument;
        arguments.push_back(std::move(materialized));
    }

    // expand/factor/simplify/collectは純粋な式変形frontendであり，外側Holdのためだけに
    // 未dispatchへ残す理由はない。一方generic evaluate()を呼ぶと制御変数のsession定義まで
    // 解決し得るため，held引数を既存kernelへ直接渡してbinder意味論を保つ。
    switch (definition->id) {
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify: {
        if (arguments.empty() || arguments.size() > 2)
            break;
        mathematics::AssumptionSet assumptions;
        if (arguments.size() == 2)
            assumptions = mathematics::parseAssumptions(arguments[1], registry_, mathematics_);
        simplification::SimplificationContext context{
            registry_, mathematics_, angleSemantics_, std::move(assumptions)};
        context.predefinedSymbols = &symbolRegistry_;
        expression::Expr result = definition->id == BuiltinId::FullSimplify
            ? simplification::fullSimplify(arguments.front(), context)
            : arguments.front();
        return simplification::simplifyExplicitLinearCombination(result, context);
    }
    case BuiltinId::Expand:
        if (arguments.size() == 1)
            return symbolic::expandExpression(
                arguments.front(), registry_, mathematics_, angleSemantics_);
        break;
    case BuiltinId::Factor:
        if (arguments.size() == 1)
            return symbolic::factorExpression(
                arguments.front(), registry_, mathematics_, angleSemantics_);
        break;
    case BuiltinId::Collect:
        if (arguments.size() == 2)
            return symbolic::collectExpression(
                arguments[0], collectVariables(arguments[1]),
                registry_, mathematics_, angleSemantics_);
        break;
    default:
        break;
    }

    return changed
        ? expression::Expr::rebuildCall(call, std::move(arguments))
        : expression;
}

expression::Expr Evaluator::dispatchBuiltin(
    const BuiltinDefinition& definition,
    const expression::CallExpr& call,
    std::span<const expression::Expr> arguments) {
    if (const auto exceptional = exceptionalArithmeticResult(
            definition, arguments, registry_, symbolRegistry_, mathematics_))
        return *exceptional;

    if (evaluation::requiresRectangularArray(definition.id)) {
        const bool nonRectangular = std::any_of(arguments.begin(), arguments.end(),
            [](const expression::Expr& value) { return value.isList(); });
        if (nonRectangular) {
            emitWarning("Array::nonRectangular",
                std::string{definition.name()}
                    + " requires a rectangular dense array; the brace value remains unevaluated");
            return expression::Expr::rebuildCall(
                call, std::vector<expression::Expr>{arguments.begin(), arguments.end()});
        }
    }

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

        std::vector<expression::Symbol> derivativeVariables;
        derivativeVariables.reserve(arguments.size() - 1);
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            if (arguments[i].isSymbol())
                derivativeVariables.push_back(arguments[i].asSymbol());
            else if (arguments[i].isArray()) {
                const auto& spec = arguments[i].asArray();
                if (spec.rank() == 1 && spec.size() == 2 && spec.element(0).isSymbol())
                    derivativeVariables.push_back(spec.element(0).asSymbol());
            }
        }

        // DはHoldAllだが，明示的なIn/Out/%参照はhistory snapshotとして先に解決する。
        expression::Expr result = resolveHeldHistoryReferences(arguments[0]);
        // 外側Dの変数は内側limit/integrateの点・境界でも自由記号として扱う。
        // これを保護しないと y:=5 の下で D[limit[...,x,y],y] が5を取り込む。
        result = materializeSafeHeldFrontends(result, derivativeVariables);
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            expression::Symbol variable;
            std::uint64_t order = 1;
            if (arguments[i].isSymbol()) {
                variable = arguments[i].asSymbol();
            }
            else if (arguments[i].isArray()) {
                const auto& spec = arguments[i].asArray();
                const expression::Expr specVariable = spec.size() > 0 ? spec.element(0) : arguments[i];
                const expression::Expr specOrder = spec.size() > 1 ? spec.element(1) : arguments[i];
                if (spec.shape.size() != 1 || spec.shape[0] != 2 || spec.size() != 2
                    || !specVariable.isSymbol() || !specOrder.isNumber()
                    || !specOrder.asNumber().isReal()
                    || !specOrder.asNumber().asReal().isInteger()) {
                    error::throwCalcError(error::CalcErrorType::Type,
                        "D derivative specification must be a symbol or {symbol, nonnegative integer}");
                }
                const auto parsed = numeric::tryToUint64(
                    specOrder.asNumber().asReal().asInteger());
                if (!parsed)
                    error::throwCalcError(error::CalcErrorType::Domain,
                        "D derivative order must be a nonnegative integer that fits in uint64");
                variable = specVariable.asSymbol();
                order = *parsed;
            }
            else {
                error::throwCalcError(error::CalcErrorType::Type,
                    "D derivative specification must be a symbol or {symbol, nonnegative integer}");
            }

            if (order > 1) {
                if (auto knownRepeated = symbolic::differentiateKnownRepeatedExpression(
                        result, variable, order, registry_, mathematics_, angleSemantics_)) {
                    result = std::move(*knownRepeated);
                    continue;
                }
            }
            for (std::uint64_t derivative = 0; derivative < order; ++derivative) {
                // 高階Dは固定order上限ではなくrequest budgetで制御する。
                // exact 0へ到達した後は残りorderを反復する必要がない。
                consumeEvaluationBudget(EvaluationResource::EvaluationStep);
                result = symbolic::differentiateExpression(
                    result, variable, registry_, mathematics_, angleSemantics_);
                if (const auto normalized = symbolic::normalizeRationalExpression(
                        result, variable, registry_, mathematics_, angleSemantics_))
                    result = *normalized;

                if (order > 1)
                    result = symbolic::canonicalizeDerivativeOutput(
                        result, registry_, mathematics_, angleSemantics_);
                if (result.isNumber() && result.asNumber().isZero())
                    break;
            }
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

        std::vector<expression::Symbol> integrationVariables;
        if (arguments[1].isSymbol())
            integrationVariables.push_back(arguments[1].asSymbol());
        else if (const auto iterator = parseRangeIteratorSpec(arguments[1]))
            integrationVariables.push_back(iterator->variable);

        expression::Expr integrand = materializeSafeHeldFrontends(
            resolveHeldHistoryReferences(arguments[0]), integrationVariables);
        // 外側の積分変数を保護したまま，内側の安全なsymbolic frontendだけを先に閉じる。
        std::optional<expression::Expr> result;
        std::optional<symbolic::IntegrationDisposition> integrationDisposition;
        if (arguments[1].isSymbol()) {
            symbolic::IntegrationResult detailed = symbolic::integrateExpressionDetailed(
                integrand, arguments[1].asSymbol(), registry_, mathematics_,
                angleSemantics_, assumptions);
            integrationDisposition = detailed.disposition;
            result = std::move(detailed.expression);
        }
        else if (const auto iterator = parseRangeIteratorSpec(arguments[1])) {
            const auto* infinity = symbolRegistry_.find("Infinity");
            if (!infinity)
                error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
            result = symbolic::integrateExpression(
                integrand, iterator->variable, iterator->lower, iterator->upper,
                registry_, mathematics_, angleSemantics_, infinity->symbol, assumptions);
        }
        else {
            error::throwCalcError(error::CalcErrorType::Type,
                "integrate expects a symbol or {variable, lower, upper} as the second argument");
        }

        if (containsBuiltinCall(*result, registry_.symbol(BuiltinId::SymbolicIntegral))) {
            if (!integrationDisposition) {
                emitWarning("integrate::conditionsRequired",
                    "integrate kept the definite integral unevaluated because a safe symbolic result could not be established on the requested interval");
            }
            else {
                switch (*integrationDisposition) {
                case symbolic::IntegrationDisposition::Partial:
                    emitWarning("integrate::partial",
                        "integrate partially evaluated the expression; remaining subintegral(s) are outside the current symbolic rule set");
                    break;
                case symbolic::IntegrationDisposition::KnownNoFiniteClosedForm:
                    emitWarning("integrate::noKnownClosedForm",
                        "integrate recognized a family with no known finite closed form in mmCal's supported standard-function vocabulary; the integral remains unevaluated");
                    break;
                case symbolic::IntegrationDisposition::ConditionsRequired:
                    emitWarning("integrate::conditionsRequired",
                        "integrate needs additional domain or branch assumptions before it can choose a safe symbolic antiderivative");
                    break;
                case symbolic::IntegrationDisposition::UnsupportedByEngine:
                    emitWarning("integrate::unsupported",
                        "mmCal has no implemented symbolic integration rule for this expression; this does not imply that no closed form exists");
                    break;
                case symbolic::IntegrationDisposition::Solved:
                    emitWarning("integrate::unsupported",
                        "integrate left an unexpected unevaluated subintegral; this does not imply that no closed form exists");
                    break;
                }
            }
        }
        return *result;
    }
    case BuiltinId::Limit: {
        if (arguments.size() < 2 || arguments.size() > 4)
            error::throwCalcError(error::CalcErrorType::Type,
                "limit expects limit[expression, variable, point], limit[expression, variable, point, direction], or limit[expression, {variable, point, direction}]");

        expression::Symbol variable;
        std::optional<expression::Expr> point;
        std::optional<expression::Expr> directionExpression;
        if (arguments.size() == 2) {
            const auto spec = parseRangeIteratorSpec(arguments[1]);
            if (!spec)
                error::throwCalcError(error::CalcErrorType::Type,
                    "limit iterator must be {variable, point, direction}");
            variable = spec->variable;
            point = spec->lower;
            directionExpression = spec->upper;
        }
        else {
            if (!arguments[1].isSymbol())
                error::throwCalcError(error::CalcErrorType::Type,
                    "limit variable must be a symbol");
            variable = arguments[1].asSymbol();
            point = arguments[2];
            if (arguments.size() == 4)
                directionExpression = arguments[3];
        }

        symbolic::LimitDirection direction = symbolic::LimitDirection::TwoSided;
        if (directionExpression) {
            const auto& directionValue = *directionExpression;
            if (!directionValue.isNumber() || !directionValue.asNumber().isReal()
                || !directionValue.asNumber().asReal().isInteger())
                error::throwCalcError(error::CalcErrorType::Type,
                    "limit direction must be -1 for left or 1 for right");
            const auto& value = directionValue.asNumber().asReal().asInteger();
            if (value == numeric::BigInt{-1})
                direction = symbolic::LimitDirection::Left;
            else if (value == numeric::BigInt{1})
                direction = symbolic::LimitDirection::Right;
            else
                error::throwCalcError(error::CalcErrorType::Domain,
                    "limit direction must be -1 for left or 1 for right");
        }

        const auto* infinity = symbolRegistry_.find("Infinity");
        const auto* complexInfinity = symbolRegistry_.find("ComplexInfinity");
        const auto* indeterminate = symbolRegistry_.find("Indeterminate");
        if (!infinity || !complexInfinity || !indeterminate)
            error::throwCalcError(error::CalcErrorType::Internal,
                "Limit exceptional symbols are not registered");
        expression::Expr result = symbolic::limitExpression(
            arguments[0], variable, *point, direction,
            registry_, mathematics_, angleSemantics_, infinity->symbol, {},
            &complexInfinity->symbol, &indeterminate->symbol);

        // 点代入後に integrate[0,x] や別変数の内側limitが残る場合，
        // そのbinder自身は既に外側極限の局所assumptionを通過しているのでここで閉じる。
        // 同じ変数の未解決limitを再評価すると自己再帰になるため，その形だけは保持する。
        bool mayMaterializeResidual = true;
        if (registry_.isCallTo(result, BuiltinId::Limit)) {
            const auto& residualArguments = result.asCall().arguments;
            if (residualArguments.size() >= 2) {
                expression::Symbol residualVariable;
                bool parsed = false;
                if (residualArguments[1].isSymbol()) {
                    residualVariable = residualArguments[1].asSymbol();
                    parsed = true;
                }
                else if (const auto residualSpec = parseRangeIteratorSpec(residualArguments[1])) {
                    residualVariable = residualSpec->variable;
                    parsed = true;
                }
                if (parsed && residualVariable.sameIdentity(variable))
                    mayMaterializeResidual = false;
            }
        }
        if (mayMaterializeResidual) {
            const std::array<expression::Symbol, 1> limitVariables{variable};
            result = materializeSafeHeldFrontends(result, limitVariables);
        }

        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Limit)))
            emitWarning("limit::unevaluated",
                "limit could not prove the requested limit; unevaluated limit[...] remains");
        return result;
    }
    case BuiltinId::Series: {
        if (arguments.size() < 2 || arguments.size() > 3)
            error::throwCalcError(error::CalcErrorType::Type,
                "series expects series[expression, {variable, center, order}] with optional assumptions");

        if (!arguments[1].isArray())
            error::throwCalcError(error::CalcErrorType::Type,
                "series specification must be {variable, center, nonnegative integer order}");
        const auto& spec = arguments[1].asArray();
        if (spec.rank() != 1 || spec.size() != 3 || !spec.element(0).isSymbol())
            error::throwCalcError(error::CalcErrorType::Type,
                "series specification must be {variable, center, nonnegative integer order}");
        const expression::Expr orderExpression = spec.element(2);
        if (!orderExpression.isNumber() || !orderExpression.asNumber().isReal()
            || !orderExpression.asNumber().asReal().isInteger())
            error::throwCalcError(error::CalcErrorType::Type,
                "series order must be a nonnegative integer");
        const auto order = numeric::tryToUint64(
            orderExpression.asNumber().asReal().asInteger());
        if (!order)
            error::throwCalcError(error::CalcErrorType::Domain,
                "series order must be a nonnegative integer that fits in uint64");
        if (*order > 1024)
            error::throwCalcError(error::CalcErrorType::ResourceLimit,
                "series order exceeds the current limit of 1024");

        mathematics::AssumptionSet assumptions;
        if (arguments.size() == 3)
            assumptions = mathematics::parseAssumptions(arguments[2], registry_, mathematics_);

        const std::array<expression::Symbol, 1> seriesVariables{spec.element(0).asSymbol()};
        expression::Expr held = materializeSafeHeldFrontends(
            resolveHeldHistoryReferences(arguments[0]), seriesVariables);
        // 外側Seriesの変数を保護し，内側binderの点・境界へsession値が漏れないようにする。
        if (auto result = symbolic::seriesExpression(
                held, spec.element(0).asSymbol(), spec.element(1),
                static_cast<std::size_t>(*order), registry_, mathematics_,
                angleSemantics_, assumptions))
            return *result;

        emitWarning("series::unsupported",
            "series could not construct a supported local expansion; the request remains unevaluated");
        return expression::Expr::call(registry_.symbol(BuiltinId::Series),
            std::vector<expression::Expr>{arguments.begin(), arguments.end()});
    }
    case BuiltinId::SeriesData: {
        expression::Expr value = expression::Expr::call(
            registry_.symbol(BuiltinId::SeriesData),
            std::vector<expression::Expr>{arguments.begin(), arguments.end()});
        if (!symbolic::parseSeriesData(value, registry_))
            error::throwCalcError(error::CalcErrorType::Type,
                "seriesData expects six arguments, optionally followed by logarithmic coefficient layers");
        return value;
    }
    case BuiltinId::Normal: {
        if (arguments.size() != 1)
            error::throwCalcError(error::CalcErrorType::Type,
                "normal expects one expression");
        if (const auto series = symbolic::parseSeriesData(arguments.front(), registry_))
            return symbolic::normalSeriesExpression(
                *series, registry_, mathematics_, angleSemantics_);
        return arguments.front();
    }
    case BuiltinId::ToNormal: {
        if (arguments.size() != 1)
            error::throwCalcError(error::CalcErrorType::Type,
                "toNormal expects one expression");

        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(
                error::CalcErrorType::Internal,
                "Infinity symbol is not registered");

        plot::PlotPipelineOptions plotOptions;
        plotOptions.forceSampledGeometry = true;
        if (const auto requests = parsePlotRequests(arguments.front(), registry_)) {
            const auto built = plot::buildPlotPipeline(
                *requests, registry_, mathematics_, angleSemantics_, infinity->symbol,
                {}, plotOptions);
            if (!built || !built.output)
                throwPlotPipelineFailure(
                    built.status,
                    "toNormal could not materialize Plot coordinates");
            return sampledCurvesCoordinates(built.output->curves);
        }

        if (registry_.isCallTo(arguments.front(), BuiltinId::Show)) {
            if (const auto composite = parseCompositeDrawable(arguments.front(), registry_)) {
                std::vector<expression::Expr> coordinates;
                if (!composite->plots.empty()) {
                    const auto built = plot::buildPlotPipeline(
                        composite->plots, registry_, mathematics_, angleSemantics_,
                        infinity->symbol, {}, plotOptions);
                    if (!built || !built.output)
                        throwPlotPipelineFailure(
                            built.status,
                            "toNormal could not materialize Show Plot coordinates");
                    for (const auto& curve : built.output->curves)
                        coordinates.push_back(sampledCurveCoordinates(curve));
                }
                for (const auto& listPlot : composite->listPlots)
                    coordinates.push_back(listPlotNormalCoordinates(listPlot));
                if (coordinates.size() == 1)
                    return coordinates.front();
                return expression::Expr::list(std::move(coordinates));
            }
        }

        if (const auto request = parseParametricPlotRequest(arguments.front(), registry_)) {
            const auto built = plot::buildParametricPlotPipeline(
                *request, registry_, mathematics_, angleSemantics_, infinity->symbol,
                {}, plotOptions);
            if (!built || !built.output)
                throwPlotPipelineFailure(
                    built.status,
                    "toNormal could not materialize ParametricPlot coordinates");
            return sampledCurvesCoordinates(built.output->curves);
        }

        if (const auto request = parseListPlotRequest(arguments.front(), registry_))
            return listPlotNormalCoordinates(*request);

        return symbolic::toNormalExpression(
            arguments.front(), registry_, mathematics_, angleSemantics_);
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
    case BuiltinId::BitAnd:
        return builtins::evaluateBitAnd(arguments, registry_);
    case BuiltinId::BitOr:
        return builtins::evaluateBitOr(arguments, registry_);
    case BuiltinId::BitXor:
        return builtins::evaluateBitXor(arguments, registry_);
    case BuiltinId::BitNot:
        return builtins::evaluateBitNot(arguments, registry_);
    case BuiltinId::BitShiftLeft:
        return builtins::evaluateBitShiftLeft(arguments, registry_);
    case BuiltinId::BitShiftRight:
        return builtins::evaluateBitShiftRight(arguments, registry_);
    case BuiltinId::BitLength:
        return builtins::evaluateBitLength(arguments, registry_);
    case BuiltinId::BitCount:
        return builtins::evaluateBitCount(arguments, registry_);
    case BuiltinId::BitGet:
        return builtins::evaluateBitGet(arguments, registry_);
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
    case BuiltinId::IsPrime:
        return builtins::evaluateIsPrime(arguments, registry_);
    case BuiltinId::NextPrime:
        return builtins::evaluateNextPrime(arguments, registry_);
    case BuiltinId::PreviousPrime:
        return builtins::evaluatePreviousPrime(arguments, registry_);
    case BuiltinId::FactorInteger:
        return builtins::evaluateFactorInteger(arguments, registry_);
    case BuiltinId::Totient:
        return builtins::evaluateTotient(arguments, registry_);
    case BuiltinId::Permutation:
        return builtins::evaluatePermutation(arguments, registry_);
    case BuiltinId::Combination:
        return builtins::evaluateCombination(arguments, registry_);
    case BuiltinId::Fibonacci:
        return builtins::evaluateFibonacci(arguments, registry_);
    case BuiltinId::DiscreteFourierTransform:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateDft(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateDft(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::FastFourierTransform:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateFft(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateFft(arguments, registry_, mathematics_, angleSemantics_, fourierTransformCache_);
    case BuiltinId::InverseFourierTransform:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateIfft(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateIfft(arguments, registry_, mathematics_, angleSemantics_, fourierTransformCache_);
    case BuiltinId::Convolution:
        return builtins::evaluateConvolution(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Transpose:
        return builtins::evaluateTranspose(arguments, registry_);
    case BuiltinId::ConjugateTranspose:
        return builtins::evaluateConjugateTranspose(
            arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::MatrixAdd:
        return builtins::evaluateMatrixAdd(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::MatrixMultiply:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateDot(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateDot(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Determinant:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateDeterminant(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateDeterminant(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Inverse:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateInverse(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateMatrixInverse(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Rref: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateRref(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateRref(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Rref)))
            emitWarning("rref::unevaluated",
                "rref could not determine the symbolic pivots; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::Rank: {
        // rankは不連続量なので、exact入力ではまずexact eliminationを優先する。
        // exactに決まらない場合だけcertified intervalへ降ろし、epsilon判定は導入しない。
        expression::Expr result = builtins::evaluateMatrixRank(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Rank))) {
            if (const auto* approximation = currentApproximationContext())
                if (const auto approximateResult = builtins::evaluateApproximateMatrixRank(
                    arguments, registry_, mathematics_, angleSemantics_, *approximation))
                    return *approximateResult;
            emitWarning("matrixRank::unevaluated",
                "matrixRank could not determine the symbolic pivots; the expression remains unevaluated");
        }
        return result;
    }
    case BuiltinId::SolveLinear: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateSolveLinear(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateSolveLinear(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::SolveLinear)))
            emitWarning("solveLinear::unevaluated",
                "solveLinear could not determine the symbolic pivots; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::NullSpace: {
        // nullSpaceもrankと同じくrank deficiencyに依存する不連続量なので，
        // exact入力ではまずexact eliminationを優先する。未解決時だけcertified intervalへ降ろす。
        expression::Expr result = builtins::evaluateNullSpace(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::NullSpace))) {
            if (const auto* approximation = currentApproximationContext())
                if (const auto approximateResult = builtins::evaluateApproximateNullSpace(
                    arguments, registry_, mathematics_, angleSemantics_, *approximation))
                    return *approximateResult;
            emitWarning("nullSpace::unevaluated",
                "nullSpace could not certify the pivot structure; the expression remains unevaluated");
        }
        return result;
    }
    case BuiltinId::LuDecomposition: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateLuDecomposition(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateLuDecomposition(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::LuDecomposition)))
            emitWarning("luDecomposition::unevaluated",
                "luDecomposition could not prove a required symbolic pivot nonzero; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::QrDecomposition: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateQrDecomposition(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateQrDecomposition(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::QrDecomposition)))
            emitWarning("qrDecomposition::unevaluated",
                "qrDecomposition requires an exact real matrix for the exact fraction-free backend; N[...] uses the certified Householder backend for exact matrices, while existing finite-precision matrix leaves remain conservative");
        return result;
    }
    case BuiltinId::SingularValueDecomposition: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateSingularValueDecomposition(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateSingularValueDecomposition(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::SingularValueDecomposition)))
            emitWarning("svd::unevaluated",
                "svd exact form is only emitted for natural exact cases; N[...] uses the numerical SVD backend for exact matrices, while existing finite-precision matrix leaves remain conservative");
        return result;
    }
    case BuiltinId::ConditionNumber: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateConditionNumber(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateConditionNumber(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::ConditionNumber)))
            emitWarning("conditionNumber::unevaluated",
                "conditionNumber exact form is unavailable; N[...] uses the numerical SVD backend for exact matrices, while existing finite-precision matrix leaves remain conservative");
        return result;
    }
    case BuiltinId::PseudoInverse: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximatePseudoInverse(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluatePseudoInverse(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::PseudoInverse)))
            emitWarning("pseudoInverse::unevaluated",
                "pseudoInverse could not determine a safe exact or numerical rank; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::LeastSquares: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateLeastSquares(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateLeastSquares(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::LeastSquares)))
            emitWarning("leastSquares::unevaluated",
                "leastSquares could not determine a safe exact or numerical rank; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::Eigenvalues: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateEigenvalues(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateEigenvalues(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Eigenvalues)))
            emitWarning("eigenvalues::unevaluated",
                "eigenvalues exact form is unavailable; N[...] uses the numerical Schur backend for exact matrices, while existing finite-precision matrix leaves remain conservative");
        return result;
    }
    case BuiltinId::Eigenvectors: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateEigenvectors(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateEigenvectors(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Eigenvectors)))
            emitWarning("eigenvectors::unevaluated",
                "eigenvectors exact form is only emitted for natural exact cases; N[...] uses the numerical Schur backend for exact matrices, while existing finite-precision matrix leaves remain conservative");
        return result;
    }
    case BuiltinId::Eigensystem: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateEigensystem(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateEigensystem(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Eigensystem)))
            emitWarning("eigensystem::unevaluated",
                "eigensystem exact form is only emitted for natural exact cases; N[...] uses the numerical Schur backend for exact matrices, while existing finite-precision matrix leaves remain conservative");
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
    case BuiltinId::Fma:
        return builtins::evaluateFma(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Clamp:
        return builtins::evaluateClamp(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Proj:
        return builtins::evaluateProj(arguments, registry_);
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
    case BuiltinId::Range:
        return builtins::evaluateRange(arguments);
    case BuiltinId::Map:
    case BuiltinId::Table:
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "map/table must be handled by the evaluation machine");
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
        return builtins::evaluateStatistic(
            definition.id, arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Dimensions:
        return builtins::evaluateDimensions(arguments);
    case BuiltinId::ArrayRank:
        return builtins::evaluateArrayRank(arguments);
    case BuiltinId::Length:
        return builtins::evaluateLength(arguments);
    case BuiltinId::ArrayGet:
        return builtins::evaluateArrayGet(arguments);
    case BuiltinId::Reshape:
        return builtins::evaluateReshape(arguments);
    case BuiltinId::Identity:
        return builtins::evaluateIdentity(arguments);
    case BuiltinId::Zeros:
        return builtins::evaluateZeros(arguments);
    case BuiltinId::Trace:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateTrace(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
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
    case BuiltinId::VectorCross:
        return builtins::evaluateVectorCross(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorNorm:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateNorm(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateNorm(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorManhattan:
        return builtins::evaluateVectorManhattan(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorEuclidean:
        return builtins::evaluateVectorEuclidean(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorNormalize:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateNormalize(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateNormalize(arguments, registry_, mathematics_, angleSemantics_);
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
    case BuiltinId::VectorInner:
        return builtins::evaluateInner(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorOuter:
        return builtins::evaluateOuter(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorRejection:
        return builtins::evaluateRejection(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::OrthogonalQ:
        return builtins::evaluateOrthogonalQ(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::OrthonormalQ:
        return builtins::evaluateOrthonormalQ(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::LinearIndependentQ:
        return builtins::evaluateLinearIndependentQ(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::GramSchmidt:
        return builtins::evaluateGramSchmidt(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Gradient:
    case BuiltinId::Divergence:
    case BuiltinId::Curl:
    case BuiltinId::Laplacian:
    case BuiltinId::Jacobian:
    case BuiltinId::Hessian:
    case BuiltinId::DirectionalDerivative: {
        const std::vector<expression::Symbol> calculusVariables = collectVariables(
            definition.id == BuiltinId::DirectionalDerivative ? arguments[2] : arguments[1]);
        expression::Expr field = materializeSafeHeldFrontends(
            resolveHeldHistoryReferences(arguments[0]), calculusVariables);
        if (definition.id == BuiltinId::DirectionalDerivative) {
            std::array<expression::Expr, 3> resolved{
                std::move(field),
                evaluate(resolveHeldHistoryReferences(arguments[1])),
                resolveHeldHistoryReferences(arguments[2])
            };
            return builtins::evaluateDirectionalDerivative(
                resolved, registry_, mathematics_, angleSemantics_);
        }

        std::array<expression::Expr, 2> resolved{
            std::move(field),
            resolveHeldHistoryReferences(arguments[1])
        };
        switch (definition.id) {
        case BuiltinId::Gradient:
            return builtins::evaluateGradient(resolved, registry_, mathematics_, angleSemantics_);
        case BuiltinId::Divergence:
            return builtins::evaluateDivergence(resolved, registry_, mathematics_, angleSemantics_);
        case BuiltinId::Curl:
            return builtins::evaluateCurl(resolved, registry_, mathematics_, angleSemantics_);
        case BuiltinId::Laplacian:
            return builtins::evaluateLaplacian(resolved, registry_, mathematics_, angleSemantics_);
        case BuiltinId::Jacobian:
            return builtins::evaluateJacobian(resolved, registry_, mathematics_, angleSemantics_);
        case BuiltinId::Hessian:
            return builtins::evaluateHessian(resolved, registry_, mathematics_, angleSemantics_);
        default:
            break;
        }
        throw std::logic_error("Unknown vector-calculus builtin");
    }
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
    case BuiltinId::LambertW:
    case BuiltinId::Zeta:
    case BuiltinId::Digamma:
    case BuiltinId::Trigamma:
    case BuiltinId::IncompleteBeta:
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS:
    case BuiltinId::Hypergeometric1F1:
    case BuiltinId::Hypergeometric2F1:
    case BuiltinId::EllipticF:
    case BuiltinId::EllipticE:
    case BuiltinId::EllipticPi:
    case BuiltinId::ExponentialIntegralEi:
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
    case BuiltinId::LogarithmicIntegralLi:
    case BuiltinId::Polylog:
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
        error::throwCalcError(error::CalcErrorType::Internal,
            "N must be evaluated through the precision-aware evaluation path");
    case BuiltinId::Precision: {
        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
        if (const auto result = builtins::evaluatePrecision(arguments, infinity->symbol))
            return *result;
        emitWarning("precision::unevaluated",
            "precision could not determine the guaranteed precision; the expression remains unevaluated");
        return expression::Expr::rebuildCall(call, {arguments.front()});
    }
    case BuiltinId::Explain:
        return builtins::evaluateExplain(arguments, registry_, symbolRegistry_, mathematics_);
    case BuiltinId::Accuracy: {
        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
        if (const auto result = builtins::evaluateAccuracy(arguments, infinity->symbol))
            return *result;
        emitWarning("accuracy::unevaluated",
            "accuracy could not determine the guaranteed accuracy; the expression remains unevaluated");
        return expression::Expr::rebuildCall(call, {arguments.front()});
    }
    case BuiltinId::Rationalize:
        if (const auto result = builtins::evaluateRationalize(arguments))
            return *result;
        emitWarning("rationalize::unevaluated",
            "rationalize could not convert part of the expression; it remains unevaluated");
        return expression::Expr::rebuildCall(
            call, std::vector<expression::Expr>{arguments.begin(), arguments.end()});
    case BuiltinId::Root: {
        if (arguments.size() < 2 || arguments.size() > 3)
            error::throwCalcError(error::CalcErrorType::Type,
                "root expects root[coefficients,index] or root[coefficients,index,Complex]");
        const auto coefficients = symbolic::rootPolynomialCoefficients(arguments[0]);
        const auto index = symbolic::positiveRootIndex(arguments[1]);
        if (!coefficients || coefficients->size() < 2 || !index)
            error::throwCalcError(error::CalcErrorType::Type,
                "root expects exact real Rational coefficients and a positive integer index");
        const bool complexDomain = arguments.size() == 3;
        if (complexDomain && (!arguments[2].isSymbol() || arguments[2].asSymbol().view() != "Complex"))
            error::throwCalcError(error::CalcErrorType::Type,
                "root third argument must be Complex");
        if (coefficients->size() == 2) {
            if (*index != 1)
                error::throwCalcError(error::CalcErrorType::Domain,
                    "root index exceeds the number of roots");
            return expression::Expr{numeric::Number{-(*coefficients)[0] / (*coefficients)[1]}};
        }
        if (complexDomain) {
            const auto algebraic = symbolic::ComplexAlgebraicNumber::create(*coefficients, *index);
            if (!algebraic)
                error::throwCalcError(error::CalcErrorType::Domain,
                    "complex root index is invalid or root isolation exceeded the current budget");
            return symbolic::makeCanonicalRootExpression(*algebraic, registry_);
        }
        const auto algebraic = symbolic::RealAlgebraicNumber::create(*coefficients, *index);
        if (!algebraic)
            error::throwCalcError(error::CalcErrorType::Domain,
                "root index is invalid or the polynomial exceeds the current algebraic degree limit");
        return symbolic::makeCanonicalRootExpression(*algebraic, registry_);
    }
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify: {
        if (arguments.empty() || arguments.size() > 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "simplify expects an expression and optional assumptions");
        mathematics::AssumptionSet assumptions;
        if (arguments.size() == 2)
            assumptions = mathematics::parseAssumptions(arguments[1], registry_, mathematics_);
        simplification::SimplificationContext context{
            registry_, mathematics_, angleSemantics_, std::move(assumptions)};
        context.predefinedSymbols = &symbolRegistry_;
        if (definition.id == BuiltinId::FullSimplify) {
            expression::Expr result = simplification::fullSimplify(arguments.front(), context);
            return simplification::simplifyExplicitLinearCombination(result, context);
        }
        return simplification::simplifyExplicitLinearCombination(arguments.front(), context);
    }
    case BuiltinId::Expand:
        return symbolic::expandExpression(
            arguments.front(), registry_, mathematics_, angleSemantics_);
    case BuiltinId::Factor:
        return symbolic::factorExpression(
            arguments.front(), registry_, mathematics_, angleSemantics_);
    case BuiltinId::Collect:
        if (arguments.size() != 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "collect expects an expression and a symbol or symbol array");
        return symbolic::collectExpression(
            arguments[0], collectVariables(arguments[1]),
            registry_, mathematics_, angleSemantics_);
    case BuiltinId::Solve: {
        consumeEvaluationBudget(EvaluationResource::SolverBranch);
        if (arguments.size() < 2 || arguments.size() > 3)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "solve expects equation(s), variable(s), and optional constraints");

        std::vector<expression::Symbol> variables;
        solver::SolveConstraints constraints;

        // solve[equation, Real] 等は、方程式中の未知symbolが一意のときだけ
        // domain指定の短縮形として受理する。複数候補から変数を推測しない。
        const auto shorthandDomain = arguments.size() == 2
            ? solveDomainSymbol(arguments[1], symbolRegistry_)
            : std::nullopt;
        if (shorthandDomain) {
            collectSolveUnknowns(arguments[0], symbolRegistry_, registry_, variables);
            if (variables.size() != 1) {
                const std::string message = variables.empty()
                    ? "solve[equation, domain] requires exactly one unknown symbol; none was found"
                    : "solve[equation, domain] requires exactly one unknown symbol; multiple candidates were found";
                error::throwCalcError(error::CalcErrorType::Type, message);
            }
            constraints.domain = *shorthandDomain;
        }
        else {
            if (arguments[1].isSymbol()) {
                validateSolveVariable(arguments[1].asSymbol(), symbolRegistry_, registry_);
                variables.push_back(arguments[1].asSymbol());
            }
            else if (arguments[1].isArray() && arguments[1].asArray().rank() == 1) {
                const auto& array = arguments[1].asArray();
                for (std::size_t i = 0; i < array.size(); ++i) {
                    const expression::Expr item = array.element(i);
                    if (!item.isSymbol())
                        error::throwCalcError(
                            error::CalcErrorType::Type,
                            "solve variable array must contain only symbols");
                    validateSolveVariable(item.asSymbol(), symbolRegistry_, registry_);
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
                    "solve expects a symbol, symbol array, or domain as its second argument");

            if (arguments.size() == 3)
                constraints = solver::parseSolveConstraints(
                    arguments[2], variables, registry_, mathematics_, angleSemantics_);
        }

        mathematics::AssumptionSet normalizationAssumptions = constraints.assumptions;
        if (constraints.domain
            && mathematics::isSubdomainOf(*constraints.domain, mathematics::NumericDomain::Real)) {
            // Solveのambient domainは正規化にも有効な仮定である。ここを落とすと
            // sqrt[x^2]がReal solveでもabs[x]へ還元されず，後段の既存abs solverへ
            // 接続できない。Complex solveではprincipal branchを変えない。
            for (const expression::Symbol& variable : variables)
                normalizationAssumptions.add(mathematics::elementOf(
                    expression::Expr{variable}, *constraints.domain));
        }
        const expression::Expr materializedSolveInput = materializeSafeHeldFrontends(
            resolveHeldHistoryReferences(arguments[0]), variables);
        const expression::Expr solveInput = solver::normalizeForSolve(
            materializedSolveInput, registry_, mathematics_, angleSemantics_, normalizationAssumptions);

        const mathematics::NumericDomain solveAmbientDomain = constraints.domain.value_or(
            mathematics::NumericDomain::Complex);
        std::vector<solver::SolverVariable> solverVariables;
        solverVariables.reserve(variables.size());
        for (const expression::Symbol& variable : variables)
            solverVariables.push_back(solver::SolverVariable{variable, solveAmbientDomain});

        const bool realSolveDomain = constraints.domain
            && mathematics::isSubdomainOf(
                *constraints.domain, mathematics::NumericDomain::Real);

        // 単一relationのdispatchをscalar入力とrelation配列で共有する。
        // ここが分裂すると，solve[sin[x]==0,x,Real]だけ解けて
        // solve[{sin[x]==0,x>0},x,Real]がpolynomial solverへ落ちる非対称性を生む。
        const auto solveSingleRelation = [&](const expression::Expr& relation) {
            const std::vector<solver::SolverVariable> singleVariables{solverVariables.front()};
            if (auto identical = solver::solveIdenticalEquality(
                    relation, singleVariables, registry_, mathematics_, angleSemantics_))
                return *identical;
            if (auto direct = solver::solveDirectAlgebraicBindingRelation(
                    relation, variables.front(), registry_, mathematics_))
                return *direct;
            if (auto absolute = solver::solveRealAbsoluteValueRelation(
                    relation, variables.front(), realSolveDomain,
                    registry_, mathematics_, angleSemantics_, constraints.assumptions))
                return *absolute;
            if (auto radical = solver::solveRadicalRelation(
                    relation, variables.front(), registry_, mathematics_,
                    angleSemantics_, constraints.assumptions))
                return *radical;
            if (auto principalLambert = solver::solvePrincipalLambertRelation(
                    relation, variables.front(), registry_, mathematics_,
                    angleSemantics_, constraints.assumptions))
                return *principalLambert;
            if (!realSolveDomain) {
                if (auto complexPeriodic = solver::solveComplexExponentialLogRelation(
                        relation, variables.front(), registry_, mathematics_,
                        angleSemantics_, constraints.assumptions))
                    return *complexPeriodic;
            }
            if (realSolveDomain) {
                consumeEvaluationBudget(EvaluationResource::SolverBranch, 3);
                if (auto exponential = solver::solveRealExponentialRelation(
                        relation, variables.front(), registry_, mathematics_,
                        angleSemantics_, constraints.assumptions))
                    return *exponential;
                if (auto periodic = solver::solveRealPeriodicFunctionRelation(
                        relation, variables.front(), registry_, mathematics_,
                        angleSemantics_, constraints.assumptions))
                    return *periodic;
                if (auto transcendental = solver::solveRealInjectiveFunctionRelation(
                        relation, variables.front(), registry_, mathematics_,
                        angleSemantics_, constraints.assumptions))
                    return *transcendental;
                const auto* infinity = symbolRegistry_.find("Infinity");
                if (!infinity)
                    error::throwCalcError(
                        error::CalcErrorType::Internal,
                        "Infinity symbol is not registered");
                // 高次exact polynomial等式は一般函数proofより先に共通algebraic Root経路へ渡す。
                // Complex/Realで同じcertified complex isolationを使い，Realでは証明済み実根だけを返す。
                if (auto algebraic = solver::solveRealAlgebraicPolynomialEquation(
                        relation, variables.front(), registry_, mathematics_, angleSemantics_))
                    return *algebraic;
                if (auto proved = solver::solveRealEquationByProof(
                        relation, variables.front(), registry_, mathematics_,
                        angleSemantics_, infinity->symbol, constraints.assumptions))
                    return *proved;
            }
            solver::SolutionSet polynomial = solver::solveUnivariatePolynomialRelation(
                relation, variables.front(), registry_, mathematics_, angleSemantics_);
            consumeSolverSolutionBranches(polynomial);
            return polynomial;
        };

        solver::SolutionSet solutions = [&]() {
            if (variables.size() == 1 && !solveInput.isArray())
                return solveSingleRelation(solveInput);

            std::vector<expression::Expr> equations;
            if (solveInput.isArray() && solveInput.asArray().rank() == 1)
                equations = solveInput.asArray().materialize();
            else
                equations.push_back(solveInput);

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

                solver::SolutionSet result = solveSingleRelation(*first);
                consumeSolverSolutionBranches(result);
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
                    consumeSolverSolutionBranches(result);
                }
                return result;
            }

            const mathematics::NumericDomain polynomialSystemDomain = constraints.domain
                && mathematics::isSubdomainOf(
                    *constraints.domain, mathematics::NumericDomain::Real)
                ? mathematics::NumericDomain::Real
                : mathematics::NumericDomain::Complex;
            if (auto polynomialSystem = solver::solvePolynomialSystem(
                    equations, variables, polynomialSystemDomain,
                    registry_, mathematics_, angleSemantics_))
                return *polynomialSystem;

            return solver::solveLinearPolynomialSystem(
                equations, variables, registry_, mathematics_, angleSemantics_);
        }();

        consumeSolverSolutionBranches(solutions);

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
        consumeSolverSolutionBranches(solutions);
        if (containsUnresolvedSolution(solutions))
            emitWarning("solve::unresolved",
                "solve could not determine a complete solution set; unresolved cases remain");
        return expression::Expr::solutionSet(std::move(solutions));
    }
    case BuiltinId::GroebnerBasis:
        return builtins::evaluateGroebnerBasis(arguments, registry_);
    case BuiltinId::PolynomialReduce:
        return builtins::evaluatePolynomialReduce(arguments, registry_);
    case BuiltinId::Cases:
    case BuiltinId::CaseBranch:
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "cases must be handled by the evaluation machine");
    case BuiltinId::Set:
        return evaluateSet(arguments);
    case BuiltinId::SetDelayed:
        return evaluateSetDelayed(call, arguments);
    case BuiltinId::Rule:
        return expression::Expr::rebuildCall(
            call, {arguments[0], arguments[1]});
    case BuiltinId::Less:
    case BuiltinId::LessEqual:
    case BuiltinId::Greater:
    case BuiltinId::GreaterEqual:
    case BuiltinId::Equal:
    case BuiltinId::NotEqual:
        return builtins::evaluateComparison(
            registry_.symbol(definition.id), arguments, registry_, mathematics_);
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
        return evaluateIndexedHistory(arguments, true);
    case BuiltinId::OutputHistory:
        return evaluateIndexedHistory(arguments, false);
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
    case BuiltinId::ListPlot: {
        if (arguments.empty())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "ListPlot expects ListPlot[data, options...]");
        const expression::Expr held = expression::Expr::rebuildCall(
            call, std::vector<expression::Expr>{arguments.begin(), arguments.end()});
        if (!parseListPlotRequest(held, registry_))
            error::throwCalcError(
                error::CalcErrorType::Type,
                "ListPlot expects a non-empty 1D array, Nx1 array, or Nx2 array");
        return held;
    }
    case BuiltinId::ParametricPlot: {
        if (arguments.size() < 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "ParametricPlot expects ParametricPlot[{x,y},{parameter,lower,upper},options...]");
        const auto iterator = parseRangeIteratorSpec(arguments[1]);
        if (!iterator)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "ParametricPlot expects ParametricPlot[{x,y},{parameter,lower,upper},options...]");

        const std::size_t previousProtectedSize = materializationProtectedSymbols_.size();
        if (std::find(materializationProtectedSymbols_.begin(),
                materializationProtectedSymbols_.end(), iterator->variable)
            == materializationProtectedSymbols_.end())
            materializationProtectedSymbols_.push_back(iterator->variable);

        std::optional<expression::Expr> materializedExpression;
        std::optional<expression::Expr> materializedIterator;
        try {
            materializedExpression = evaluate(resolveHeldHistoryReferences(arguments[0]));
            materializedIterator = evaluate(resolveHeldHistoryReferences(arguments[1]));
            materializationProtectedSymbols_.resize(previousProtectedSize);
        }
        catch (...) {
            materializationProtectedSymbols_.resize(previousProtectedSize);
            throw;
        }

        std::vector<expression::Expr> heldArguments;
        heldArguments.reserve(arguments.size());
        heldArguments.push_back(std::move(*materializedExpression));
        heldArguments.push_back(std::move(*materializedIterator));
        heldArguments.insert(heldArguments.end(), arguments.begin() + 2, arguments.end());
        const expression::Expr held = expression::Expr::rebuildCall(call, std::move(heldArguments));
        if (!parseParametricPlotRequest(held, registry_))
            error::throwCalcError(
                error::CalcErrorType::Type,
                "ParametricPlot expects {x,y} or {{x1,y1},{x2,y2},...} with {parameter,lower,upper}");
        return held;
    }
    case BuiltinId::Plot: {
        if (arguments.size() < 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Plot expects Plot[expression,{variable,lower,upper},options...]");
        const auto iterator = parseRangeIteratorSpec(arguments[1]);
        if (!iterator)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Plot expects Plot[expression,{variable,lower,upper},options...]");

        // Plotはbinder変数だけをsymbolicに保持し，parameter・user function・endpointは
        // 現在のsession環境で一度materializeする。x:=5の下でもplot[a*x,{x,...}]のxは
        // 保護される一方，a:=2は2へ解決される。
        const std::size_t previousProtectedSize = materializationProtectedSymbols_.size();
        if (std::find(materializationProtectedSymbols_.begin(),
                materializationProtectedSymbols_.end(), iterator->variable)
            == materializationProtectedSymbols_.end())
            materializationProtectedSymbols_.push_back(iterator->variable);

        std::optional<expression::Expr> materializedExpression;
        std::optional<expression::Expr> materializedIterator;
        try {
            materializedExpression = evaluate(resolveHeldHistoryReferences(arguments[0]));
            materializedIterator = evaluate(resolveHeldHistoryReferences(arguments[1]));
            materializationProtectedSymbols_.resize(previousProtectedSize);
        }
        catch (...) {
            materializationProtectedSymbols_.resize(previousProtectedSize);
            throw;
        }

        std::vector<expression::Expr> heldArguments;
        heldArguments.reserve(arguments.size());
        heldArguments.push_back(std::move(*materializedExpression));
        heldArguments.push_back(std::move(*materializedIterator));
        heldArguments.insert(heldArguments.end(), arguments.begin() + 2, arguments.end());
        const expression::Expr heldPlot = expression::Expr::rebuildCall(
            call, std::move(heldArguments));
        if (!parsePlotRequest(heldPlot, registry_))
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Plot expects Plot[expression,{variable,lower,upper},options...]");
        // PlotSceneをExprへ押し込まず，公開Plotは描画要求として保持する。
        // Export等のconsumerが必要なbackendへloweringする。
        return heldPlot;
    }
    case BuiltinId::Show: {
        if (arguments.empty())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Show expects one or more Plot or ListPlot expressions");

        // Showは描画済みmm geometryではなくPlot requestを合成する。各Plotの指定区間は
        // 保持し，viewportだけを全区間の包絡へ統一するため，ここではflattenだけ行う。
        std::vector<expression::Expr> resolvedArguments;
        resolvedArguments.reserve(arguments.size());
        for (const auto& argument : arguments)
            resolvedArguments.push_back(resolveHeldHistoryReferences(argument));

        std::vector<expression::Expr> flattened;
        for (const auto& resolved : resolvedArguments) {
            expression::Expr evaluated = evaluate(resolved);
            if (evaluated.isCall()) {
                const auto* nested = registry_.find(evaluated.asCall().head);
                if (nested && nested->id == BuiltinId::Show) {
                    flattened.insert(
                        flattened.end(),
                        evaluated.asCall().arguments.begin(),
                        evaluated.asCall().arguments.end());
                    continue;
                }
            }
            flattened.push_back(std::move(evaluated));
        }

        const expression::Expr heldShow = expression::Expr::call(
            registry_.symbol(BuiltinId::Show), flattened);
        if (!parseCompositeDrawable(heldShow, registry_))
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Show expects only Plot or ListPlot expressions");
        return heldShow;
    }
    case BuiltinId::Export: {
        const EvaluationContext* exportEvaluationContext = context_;
        if (arguments.size() < 2 || arguments.size() > 3)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Export expects Export[object,file] or Export[object,file,format]");

        // HoldAllの引数を一つずつnested evaluateすると現在のEvaluationContextが
        // 一時的に置き換わる。%/Out参照はcontextが生きているうちに全て解決する。
        std::vector<expression::Expr> resolvedArguments;
        resolvedArguments.reserve(arguments.size());
        for (const auto& argument : arguments)
            resolvedArguments.push_back(resolveHeldHistoryReferences(argument));

        expression::Expr drawableExpression = evaluate(resolvedArguments[0]);
        // nested evaluateは独立contextを使うため，Export自身のdiagnostic sinkへ戻す。
        context_ = exportEvaluationContext;
        const auto requests = parsePlotRequests(drawableExpression, registry_);
        const auto parametricRequest = parseParametricPlotRequest(drawableExpression, registry_);
        const auto listPlotRequest = parseListPlotRequest(drawableExpression, registry_);
        const auto compositeRequest = registry_.isCallTo(drawableExpression, BuiltinId::Show)
            ? parseCompositeDrawable(drawableExpression, registry_)
            : std::nullopt;
        if (!requests && !parametricRequest && !listPlotRequest && !compositeRequest)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Export currently expects a Plot, ParametricPlot, ListPlot, or Show expression as its first argument");

        expression::Expr destination = evaluate(resolvedArguments[1]);
        context_ = exportEvaluationContext;
        if (!destination.isString())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Export file name must evaluate to a string");

        std::optional<std::string> format;
        if (resolvedArguments.size() == 3) {
            expression::Expr evaluatedFormat = evaluate(resolvedArguments[2]);
            context_ = exportEvaluationContext;
            if (!evaluatedFormat.isString())
                error::throwCalcError(
                    error::CalcErrorType::Type,
                    "Export format must evaluate to a string");
            format = evaluatedFormat.asString();
        }

        const std::filesystem::path requestedPath{destination.asString()};
        const auto target = exportGraphicsTarget(requestedPath, format);
        if (!target)
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "Export currently supports SVG, EPS, and PDF plot output");

        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(
                error::CalcErrorType::Internal,
                "Infinity symbol is not registered");

        std::optional<graphics::GraphicsScene> graphicsScene;
        plot::PlotPipelineStatus buildStatus = plot::PlotPipelineStatus::EmptyRequest;
        if (compositeRequest) {
            auto built = buildCompositeGraphics(
                *compositeRequest, registry_, mathematics_, angleSemantics_, infinity->symbol);
            buildStatus = built.status;
            if (built)
                graphicsScene = std::move(built.scene);
        }
        else if (requests) {
            const auto built = plot::buildPlotPipeline(
                *requests, registry_, mathematics_, angleSemantics_, infinity->symbol);
            buildStatus = built.status;
            if (built && built.output)
                graphicsScene = built.output->graphicsScene;
        }
        else if (parametricRequest) {
            const auto built = plot::buildParametricPlotPipeline(
                *parametricRequest, registry_, mathematics_, angleSemantics_, infinity->symbol);
            buildStatus = built.status;
            if (built && built.output) {
                graphicsScene = built.output->graphicsScene;
                for (const auto& reduction : built.output->periodReductions) {
                    if (!reduction)
                        continue;
                    emitInfo(
                        "ParametricPlot::period",
                        "ParametricPlot reduced a repeated periodic parameter interval from ["
                            + formatting::formatExpr(reduction->originalLower) + ", "
                            + formatting::formatExpr(reduction->originalUpper) + "] to ["
                            + formatting::formatExpr(reduction->originalLower) + ", "
                            + formatting::formatExpr(reduction->effectiveUpper)
                            + "] (proven period " + formatting::formatExpr(reduction->period)
                            + ", repetitions " + std::to_string(reduction->repetitions) + ").");
                }
            }
        }
        else if (listPlotRequest) {
            graphicsScene = buildListPlotGraphics(
                *listPlotRequest, registry_, mathematics_, angleSemantics_);
            buildStatus = graphicsScene
                ? plot::PlotPipelineStatus::Success
                : plot::PlotPipelineStatus::SamplingFailed;
        }
        if (!graphicsScene)
            throwPlotPipelineFailure(
                buildStatus, "Export could not build the plot graphics scene");

        const auto rendered = graphics::renderGraphics(
            *graphicsScene, target->format);
        if (!rendered || !rendered.data)
            error::throwCalcError(
                error::CalcErrorType::Evaluation,
                rendered.status == graphics::GraphicsRenderStatus::UnsupportedFeature
                    ? "Export format does not support a graphics feature used by this plot"
                    : "Export could not render the plot in the requested format");

        std::ofstream output(target->path, std::ios::binary | std::ios::trunc);
        if (!output)
            error::throwCalcError(
                error::CalcErrorType::Evaluation,
                "Export could not open the output file");
        output.write(rendered.data->data(), static_cast<std::streamsize>(rendered.data->size()));
        if (!output)
            error::throwCalcError(
                error::CalcErrorType::Evaluation,
                "Export could not write the output file");
        return expression::Expr{target->path.string()};
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

} // namespace mmcal::evaluation
