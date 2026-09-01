#pragma once

#include "builtins/names.hpp"
#include "cli/output_layout.hpp"
#include "expression/expr.hpp"
#include "expression/array_utils.hpp"
#include "formatting/expr_formatter.hpp"
#include "solver/solution_set.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstddef>
#include <numeric>
#include <string>
#include <string_view>
#include <span>
#include <vector>

namespace mmcal::cli {

[[nodiscard]] inline bool isCasesCall(const expression::Expr& expression) noexcept {
    return expression.isCall()
        && expression.asCall().head.view() == builtins::names::cases;
}

[[nodiscard]] inline bool isSeriesDataCall(const expression::Expr& expression) noexcept {
    return expression.isCall()
        && expression.asCall().head.view() == builtins::names::seriesData
        && (expression.asCall().arguments.size() == 6
            || expression.asCall().arguments.size() == 7);
}

[[nodiscard]] inline bool prefersMultiline(const expression::Expr& expression) noexcept {
    if (expression.isArray())
        return expression.asArray().rank() >= 2 && expression.asArray().size() != 0;
    if (expression.isList()) {
        const auto& elements = expression.asList().elements;
        return std::any_of(elements.begin(), elements.end(), [](const expression::Expr& item) {
            return item.isArray() || item.isList() || item.isSolutionSet() || isCasesCall(item);
        });
    }
    if (expression.isSolutionSet()) {
        const auto& solutions = expression.asSolutionSet();
        return solutions.kind() == solver::SolutionSetKind::Conditional
            || (solutions.kind() == solver::SolutionSetKind::Finite
                && solutions.branches().size() > 1);
    }
    return isCasesCall(expression) && expression.asCall().arguments.size() > 1;
}

namespace detail {

inline void appendIndent(std::string& output, std::size_t indent) {
    output.append(indent, ' ');
}

[[nodiscard]] inline std::string canonical(const expression::Expr& expression) {
    return formatting::formatExpr(expression);
}

inline void appendPretty(
    std::string& output,
    const expression::Expr& expression,
    std::size_t indent,
    std::size_t width,
    bool forceMulti);

[[nodiscard]] inline std::vector<expression::Expr> seriesCoefficients(
    const expression::Expr& expression) {
    if (expression.isArray())
        return expression.asArray().materialize();
    if (expression.isList())
        return expression.asList().elements;
    return {};
}

inline void appendSeriesData(
    std::string& output,
    const expression::CallExpr& call) {
    const auto& arguments = call.arguments;
    if ((arguments.size() != 6 && arguments.size() != 7) || !arguments[0].isSymbol()) {
        output += canonical(expression::Expr::call(call.head, call.arguments));
        return;
    }

    const std::vector<expression::Expr> coefficients = seriesCoefficients(arguments[2]);
    if (coefficients.empty() || !arguments[3].isNumber()
        || !arguments[4].isNumber() || !arguments[5].isNumber()) {
        output += canonical(expression::Expr::call(call.head, call.arguments));
        return;
    }

    const auto parseInteger = [](const expression::Expr& value) -> std::optional<std::int64_t> {
        const std::string text = canonical(value);
        std::int64_t result = 0;
        const auto parsed = std::from_chars(text.data(), text.data() + text.size(), result);
        if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
            return std::nullopt;
        return result;
    };
    const auto minimumExponent = parseInteger(arguments[3]);
    const auto order = parseInteger(arguments[4]);
    const auto denominator = parseInteger(arguments[5]);
    if (!minimumExponent || !order || !denominator || *denominator <= 0) {
        output += canonical(expression::Expr::call(call.head, call.arguments));
        return;
    }

    std::vector<std::vector<expression::Expr>> logarithmicLayers;
    if (arguments.size() == 7) {
        if (arguments[6].isArray()) {
            const auto& array = arguments[6].asArray();
            if (array.shape.size() != 2 || array.shape[1] != coefficients.size()) {
                output += canonical(expression::Expr::call(call.head, call.arguments));
                return;
            }
            const auto flat = array.materialize();
            logarithmicLayers.reserve(array.shape[0]);
            for (std::size_t row = 0; row < array.shape[0]; ++row) {
                const auto begin = flat.begin()
                    + static_cast<std::ptrdiff_t>(row * array.shape[1]);
                logarithmicLayers.emplace_back(
                    begin, begin + static_cast<std::ptrdiff_t>(array.shape[1]));
            }
        }
        else if (arguments[6].isList()) {
            logarithmicLayers.reserve(arguments[6].asList().elements.size());
            for (const auto& layerExpression : arguments[6].asList().elements) {
                auto layer = seriesCoefficients(layerExpression);
                if (layer.size() != coefficients.size()) {
                    output += canonical(expression::Expr::call(call.head, call.arguments));
                    return;
                }
                logarithmicLayers.push_back(std::move(layer));
            }
        }
        else {
            output += canonical(expression::Expr::call(call.head, call.arguments));
            return;
        }
    }

    const std::string variable = arguments[0].asSymbol().name();
    const std::string center = canonical(arguments[1]);
    const bool atInfinity = center == "Infinity";
    std::string base;
    if (atInfinity)
        base = variable;
    else if (center == "0")
        base = variable;
    else if (!center.empty() && center.front() == '-')
        base = "(" + variable + "+" + center.substr(1) + ")";
    else
        base = "(" + variable + "-" + center + ")";

    const auto exponentText = [&](std::int64_t numerator) {
        if (*denominator == 1)
            return std::to_string(numerator);
        const std::int64_t divisor = std::gcd(
            numerator < 0 ? -numerator : numerator, *denominator);
        const std::int64_t reducedNumerator = numerator / divisor;
        const std::int64_t reducedDenominator = *denominator / divisor;
        if (reducedDenominator == 1)
            return std::to_string(reducedNumerator);
        return std::to_string(reducedNumerator) + "/" + std::to_string(reducedDenominator);
    };
    const auto powerFactor = [&](std::int64_t exponentNumerator) {
        if (exponentNumerator == 0) return std::string{};
        const std::int64_t displayNumerator =
            atInfinity ? -exponentNumerator : exponentNumerator;
        if (displayNumerator == *denominator) return base;
        const std::string exponent = exponentText(displayNumerator);
        if (*denominator == 1)
            return base + "^" + exponent;
        return base + "^(" + exponent + ")";
    };
    const auto logFactor = [&](std::size_t degree) {
        if (degree == 0) return std::string{};
        std::string factor = atInfinity
            ? "log[1/" + variable + "]"
            : "log[" + base + "]";
        if (degree != 1)
            factor += "^" + std::to_string(degree);
        return factor;
    };

    bool wroteTerm = false;
    const auto appendTerm = [&](const expression::Expr& coefficientExpression,
                                std::int64_t exponentNumerator,
                                std::size_t logDegree) {
        const std::string coefficient = canonical(coefficientExpression);
        if (coefficient == "0") return;

        const std::string power = powerFactor(exponentNumerator);
        const std::string logarithm = logFactor(logDegree);
        std::string factor;
        if (!power.empty()) factor = power;
        if (!logarithm.empty()) {
            if (!factor.empty()) factor += "*";
            factor += logarithm;
        }

        std::string term;
        if (factor.empty())
            term = coefficient;
        else if (coefficient == "1")
            term = factor;
        else if (coefficient == "-1")
            term = "-" + factor;
        else
            term = coefficient + "*" + factor;

        if (!wroteTerm) {
            output += term;
            wroteTerm = true;
        }
        else if (!term.empty() && term.front() == '-') {
            output += " - ";
            output += term.substr(1);
        }
        else {
            output += " + ";
            output += term;
        }
    };

    for (std::size_t index = 0; index < coefficients.size(); ++index) {
        const std::int64_t exponentNumerator =
            *minimumExponent + static_cast<std::int64_t>(index);
        appendTerm(coefficients[index], exponentNumerator, 0);
        for (std::size_t layer = 0; layer < logarithmicLayers.size(); ++layer)
            appendTerm(logarithmicLayers[layer][index], exponentNumerator, layer + 1);
    }

    if (!wroteTerm)
        output += "0";

    output += " + O[";
    output += base;
    const std::int64_t displayOrder = atInfinity ? -*order : *order;
    if (displayOrder != *denominator) {
        output.push_back('^');
        const std::string exponent = exponentText(displayOrder);
        if (*denominator == 1)
            output += exponent;
        else
            output += "(" + exponent + ")";
    }
    output.push_back(']');
}

inline void appendSequence(
    std::string& output,
    std::span<const expression::Expr> elements,
    std::size_t indent,
    std::size_t width,
    bool forceMulti) {
    const std::string flat = canonical(expression::Expr::list(
        std::vector<expression::Expr>{elements.begin(), elements.end()}));
    if (!forceMulti && indent + flat.size() <= width) {
        output += flat;
        return;
    }

    output.push_back('{');
    if (elements.empty()) {
        output.push_back('}');
        return;
    }
    output.push_back('\n');
    for (std::size_t i = 0; i < elements.size(); ++i) {
        appendIndent(output, indent + 2);
        appendPretty(output, elements[i], indent + 2, width, false);
        if (i + 1 != elements.size())
            output.push_back(',');
        output.push_back('\n');
    }
    appendIndent(output, indent);
    output.push_back('}');
}

inline void appendArrayDimension(
    std::string& output,
    const expression::ArrayExpr& array,
    std::size_t dimension,
    std::size_t offset,
    std::size_t indent,
    std::size_t width,
    bool forceMulti) {
    if (dimension >= array.shape.size()) {
        output += canonical(array.element(offset));
        return;
    }

    const std::size_t extent = array.shape[dimension];
    std::size_t block = 1;
    for (std::size_t i = dimension + 1; i < array.shape.size(); ++i)
        block *= array.shape[i];

    if (dimension + 1 == array.shape.size()) {
        std::string flat{"{"};
        for (std::size_t i = 0; i < extent; ++i) {
            if (i != 0)
                flat += ", ";
            flat += canonical(array.element(offset + i));
        }
        flat.push_back('}');
        if (!forceMulti && indent + flat.size() <= width) {
            output += flat;
            return;
        }
    }

    output.push_back('{');
    if (extent == 0) {
        output.push_back('}');
        return;
    }
    output.push_back('\n');
    for (std::size_t i = 0; i < extent; ++i) {
        appendIndent(output, indent + 2);
        appendArrayDimension(
            output, array, dimension + 1, offset + i * block,
            indent + 2, width, false);
        if (i + 1 != extent)
            output.push_back(',');
        output.push_back('\n');
    }
    appendIndent(output, indent);
    output.push_back('}');
}

inline void appendCases(
    std::string& output,
    const expression::CallExpr& call,
    std::size_t indent,
    std::size_t width,
    bool forceMulti) {
    const std::string flat = canonical(expression::Expr::call(call.head, call.arguments));
    if (!forceMulti && indent + flat.size() <= width) {
        output += flat;
        return;
    }

    output += "cases[";
    if (call.arguments.empty()) {
        output.push_back(']');
        return;
    }
    output.push_back('\n');
    for (std::size_t i = 0; i < call.arguments.size(); ++i) {
        appendIndent(output, indent + 2);
        const auto& branchExpression = call.arguments[i];
        if (branchExpression.isCall()
            && branchExpression.asCall().head.view() == builtins::names::caseBranch
            && !branchExpression.asCall().arguments.empty()
            && branchExpression.asCall().arguments.size() <= 2) {
            const auto& branch = branchExpression.asCall();
            appendPretty(output, branch.arguments[0], indent + 2, width, false);
            if (branch.arguments.size() == 2) {
                output += " if ";
                output += canonical(branch.arguments[1]);
            }
        }
        else {
            appendPretty(output, branchExpression, indent + 2, width, false);
        }
        if (i + 1 != call.arguments.size())
            output.push_back(';');
        output.push_back('\n');
    }
    appendIndent(output, indent);
    output.push_back(']');
}

[[nodiscard]] inline std::string_view domainName(
    mathematics::NumericDomain domain) noexcept {
    switch (domain) {
    case mathematics::NumericDomain::Integer: return "Integer";
    case mathematics::NumericDomain::Rational: return "Rational";
    case mathematics::NumericDomain::Real: return "Real";
    case mathematics::NumericDomain::Complex: return "Complex";
    case mathematics::NumericDomain::Unknown: return "Unknown";
    }
    return "Unknown";
}

inline void appendSolutionBranch(
    std::string& output,
    const solver::SolutionBranch& branch,
    std::size_t indent,
    std::size_t width) {
    const bool multipleBindings = branch.bindings.size() > 1;
    if (multipleBindings)
        output.push_back('{');

    for (std::size_t i = 0; i < branch.bindings.size(); ++i) {
        if (i != 0)
            output += ", ";
        output += branch.bindings[i].variable.name();
        output += " == ";
        appendPretty(output, branch.bindings[i].value, indent, width, false);
    }

    if (multipleBindings)
        output.push_back('}');

    if (!branch.freeVariables.empty()) {
        const bool regionBranch = branch.bindings.empty();
        if (!regionBranch)
            output += " where ";
        for (std::size_t i = 0; i < branch.freeVariables.size(); ++i) {
            if (i != 0)
                output += ", ";
            output += branch.freeVariables[i].symbol.name();
            output += " in ";
            output += domainName(branch.freeVariables[i].domain);
        }
    }

    if (!branch.conditions.empty()) {
        output += " if ";
        output += formatting::formatAssumptions(branch.conditions);
    }
    if (branch.multiplicity && *branch.multiplicity > 1) {
        output += " (multiplicity ";
        output += std::to_string(*branch.multiplicity);
        output.push_back(')');
    }
}

inline void appendSolutionBranches(
    std::string& output,
    std::span<const solver::SolutionBranch> branches,
    std::size_t indent,
    std::size_t width,
    bool forceMulti) {
    output.push_back('{');
    if (branches.empty()) {
        output.push_back('}');
        return;
    }

    if (!forceMulti) {
        for (std::size_t i = 0; i < branches.size(); ++i) {
            if (i != 0)
                output += ", ";
            appendSolutionBranch(output, branches[i], indent, width);
        }
        output.push_back('}');
        return;
    }

    output.push_back('\n');
    for (std::size_t i = 0; i < branches.size(); ++i) {
        appendIndent(output, indent + 2);
        appendSolutionBranch(output, branches[i], indent + 2, width);
        if (i + 1 != branches.size())
            output.push_back(',');
        output.push_back('\n');
    }
    appendIndent(output, indent);
    output.push_back('}');
}

inline void appendSolutionOutcome(
    std::string& output,
    solver::SolutionSetKind outcome,
    std::span<const solver::SolutionBranch> branches,
    std::size_t indent,
    std::size_t width) {
    using solver::SolutionSetKind;
    switch (outcome) {
    case SolutionSetKind::Empty:
        output += "{}";
        return;
    case SolutionSetKind::Finite:
        appendSolutionBranches(output, branches, indent, width, branches.size() > 1);
        return;
    case SolutionSetKind::Universal:
        output += "All";
        return;
    case SolutionSetKind::Unresolved:
    case SolutionSetKind::Conditional:
        output += "Unresolved";
        return;
    }
}

inline void appendGlobalSolutionConditions(
    std::string& output,
    const mathematics::AssumptionSet& conditions) {
    if (conditions.empty())
        return;
    output += " if ";
    output += formatting::formatAssumptions(conditions);
}

inline void appendSolutionSet(
    std::string& output,
    const solver::SolutionSet& solutions,
    std::size_t indent,
    std::size_t width,
    bool forceMulti) {
    using solver::SolutionSetKind;
    switch (solutions.kind()) {
    case SolutionSetKind::Empty:
        output += "{}";
        appendGlobalSolutionConditions(output, solutions.conditions());
        return;
    case SolutionSetKind::Universal:
        output += "All";
        appendGlobalSolutionConditions(output, solutions.conditions());
        return;
    case SolutionSetKind::Unresolved:
        output += "UnresolvedSolutionSet[";
        for (std::size_t i = 0; i < solutions.variables().size(); ++i) {
            if (i != 0)
                output += ", ";
            output += solutions.variables()[i].symbol.name();
        }
        output.push_back(']');
        appendGlobalSolutionConditions(output, solutions.conditions());
        return;
    case SolutionSetKind::Finite:
        appendSolutionBranches(
            output, solutions.branches(), indent, width,
            forceMulti || solutions.branches().size() > 1);
        appendGlobalSolutionConditions(output, solutions.conditions());
        return;
    case SolutionSetKind::Conditional:
        output += "cases[";
        if (solutions.cases().empty()) {
            output.push_back(']');
            appendGlobalSolutionConditions(output, solutions.conditions());
            return;
        }
        output.push_back('\n');
        for (std::size_t i = 0; i < solutions.cases().size(); ++i) {
            const auto& solutionCase = solutions.cases()[i];
            appendIndent(output, indent + 2);
            appendSolutionOutcome(
                output, solutionCase.outcome, solutionCase.branches,
                indent + 2, width);
            if (!solutionCase.conditions.empty()) {
                output += " if ";
                output += formatting::formatAssumptions(solutionCase.conditions);
            }
            if (i + 1 != solutions.cases().size())
                output.push_back(';');
            output.push_back('\n');
        }
        appendIndent(output, indent);
        output.push_back(']');
        appendGlobalSolutionConditions(output, solutions.conditions());
        return;
    }
}

inline void appendPretty(
    std::string& output,
    const expression::Expr& expression,
    std::size_t indent,
    std::size_t width,
    bool forceMulti) {
    if (expression.isArray()) {
        const auto& array = expression.asArray();
        if (!expression::braceLiteralPreservesShape(array.shape)) {
            output += canonical(expression);
            return;
        }
        appendArrayDimension(output, array, 0, 0, indent, width, forceMulti);
        return;
    }
    if (expression.isList()) {
        appendSequence(output, expression.asList().elements, indent, width, forceMulti);
        return;
    }
    if (isSeriesDataCall(expression)) {
        appendSeriesData(output, expression.asCall());
        return;
    }
    if (expression.isSolutionSet()) {
        const std::string flat = canonical(expression);
        if (!forceMulti && indent + flat.size() <= width) {
            output += flat;
            return;
        }
        appendSolutionSet(
            output, expression.asSolutionSet(), indent, width, forceMulti);
        return;
    }
    if (isCasesCall(expression)) {
        appendCases(output, expression.asCall(), indent, width, forceMulti);
        return;
    }
    output += canonical(expression);
}

} // namespace detail

[[nodiscard]] inline std::string formatReplExpression(
    const expression::Expr& expression,
    OutputLayout layout,
    std::size_t width) {
    const std::string flat = formatting::formatExpr(expression);
    if (layout == OutputLayout::Single)
        return flat;

    width = std::max<std::size_t>(width, 20);
    if (isSeriesDataCall(expression)) {
        std::string output;
        detail::appendPretty(output, expression, 0, width, false);
        return output;
    }
    const bool structured = prefersMultiline(expression);
    const bool forceMulti = layout == OutputLayout::Multi
        || (layout == OutputLayout::Auto && structured);
    if (layout == OutputLayout::Auto && !structured && flat.size() <= width)
        return flat;

    std::string output;
    detail::appendPretty(output, expression, 0, width, forceMulti);
    return output;
}

} // namespace mmcal::cli
