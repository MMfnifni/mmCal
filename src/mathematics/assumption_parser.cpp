// 仮定式の解析
#include "assumption_parser.hpp"

#include "error/error_message.hpp"
#include "knowledge_context.hpp"

#include <optional>
#include <string_view>

namespace mmcal::mathematics {
namespace {

using evaluation::BuiltinId;
using expression::Expr;

[[nodiscard]] std::optional<NumericDomain> domainFromExpr(const Expr& expression) {
    if (!expression.isSymbol())
        return std::nullopt;
    const std::string_view name = expression.asSymbol().view();
    if (name == "Integer") return NumericDomain::Integer;
    if (name == "Rational") return NumericDomain::Rational;
    if (name == "Real") return NumericDomain::Real;
    if (name == "Complex") return NumericDomain::Complex;
    return std::nullopt;
}

[[nodiscard]] std::optional<RelationKind> relationKind(BuiltinId id) noexcept {
    switch (id) {
    case BuiltinId::Equal: return RelationKind::Equal;
    case BuiltinId::NotEqual: return RelationKind::NotEqual;
    case BuiltinId::Less: return RelationKind::Less;
    case BuiltinId::LessEqual: return RelationKind::LessEqual;
    case BuiltinId::Greater: return RelationKind::Greater;
    case BuiltinId::GreaterEqual: return RelationKind::GreaterEqual;
    default: return std::nullopt;
    }
}

void addChecked(
    Predicate predicate,
    AssumptionSet& result,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    const KnowledgeContext knowledge{builtins, mathematics, result};
    const TruthValue truth = knowledge.prove(predicate);
    if (truth == TruthValue::False)
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "Assumptions are inconsistent");
    if (truth != TruthValue::True)
        result.add(std::move(predicate));
}

void appendAssumptions(
    const Expr& expression,
    AssumptionSet& result,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    if (expression.isBoolean()) {
        if (expression.asBoolean())
            return;
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "Assumptions are inconsistent");
    }

    if (expression.isArray()) {
        if (expression.asArray().rank() != 1)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Assumptions array must be one-dimensional");
        const auto& array = expression.asArray();
        for (std::size_t i = 0; i < array.size(); ++i)
            appendAssumptions(array.element(i), result, builtins, mathematics);
        return;
    }

    if (!expression.isCall())
        error::throwCalcError(
            error::CalcErrorType::Type,
            "Assumptions must be comparisons, element[expr, domain], or an array");

    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "Unknown assumption predicate");

    const auto& arguments = expression.asCall().arguments;
    if (definition->id == BuiltinId::LogicalAnd) {
        for (const Expr& item : arguments)
            appendAssumptions(item, result, builtins, mathematics);
        return;
    }

    if (definition->id == BuiltinId::Element) {
        if (arguments.size() != 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "element expects an expression and a numeric domain");
        const auto domain = domainFromExpr(arguments[1]);
        if (!domain)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "element domain must be Integer, Rational, Real, or Complex");
        addChecked(elementOf(arguments[0], *domain), result, builtins, mathematics);
        return;
    }

    if (const auto relation = relationKind(definition->id)) {
        if (arguments.size() != 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "Assumption comparison must have two arguments");
        addChecked(mathematics::relation(*relation, arguments[0], arguments[1]),
            result, builtins, mathematics);
        return;
    }

    error::throwCalcError(
        error::CalcErrorType::Type,
        "Unsupported assumption predicate");
}

} // namespace

AssumptionSet parseAssumptions(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    AssumptionSet result;
    appendAssumptions(expression, result, builtins, mathematics);
    return result;
}

} // namespace mmcal::mathematics
