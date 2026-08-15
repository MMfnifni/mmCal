// 比較演算
#include "comparison.hpp"

#include "error/error_message.hpp"
#include "names.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/rational.hpp"

#include <optional>

#include <compare>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using expression::Expr;
using expression::Symbol;

void requireBinary(std::span<const Expr> arguments, std::string_view name) {
    if (arguments.size() != 2)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects 2 arguments");
}

[[nodiscard]] Expr unresolved(
    const Symbol& head,
    const Expr& lhs,
    const Expr& rhs) {
    return Expr::call(head, {lhs, rhs});
}

struct RealBounds final {
    numeric::Rational lower;
    numeric::Rational upper;
};

struct ComplexBounds final {
    RealBounds real;
    RealBounds imaginary;
};

[[nodiscard]] std::optional<RealBounds> informationBounds(const Expr& value) {
    if (value.isNumber() && value.asNumber().isReal()) {
        const numeric::Rational exact = value.asNumber().asReal().toRational();
        return RealBounds{exact, exact};
    }
    if (value.isDecimalApproximation()) {
        const auto& approximate = value.asDecimalApproximation();
        return RealBounds{approximate.informationLower(), approximate.informationUpper()};
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<ComplexBounds> complexInformationBounds(const Expr& value) {
    if (value.isNumber()) {
        const auto& number = value.asNumber();
        if (number.isReal()) {
            const numeric::Rational exact = number.asReal().toRational();
            return ComplexBounds{RealBounds{exact, exact}, RealBounds{numeric::Rational{}, numeric::Rational{}}};
        }
        const auto& complex = number.asComplex();
        return ComplexBounds{
            RealBounds{complex.real.toRational(), complex.real.toRational()},
            RealBounds{complex.imaginary.toRational(), complex.imaginary.toRational()}};
    }
    if (value.isDecimalApproximation()) {
        const auto& approximate = value.asDecimalApproximation();
        return ComplexBounds{
            RealBounds{approximate.informationLower(), approximate.informationUpper()},
            RealBounds{numeric::Rational{}, numeric::Rational{}}};
    }
    if (value.isComplexDecimalApproximation()) {
        const auto& approximate = value.asComplexDecimalApproximation();
        return ComplexBounds{
            RealBounds{approximate.real().informationLower(), approximate.real().informationUpper()},
            RealBounds{approximate.imaginary().informationLower(), approximate.imaginary().informationUpper()}};
    }
    return std::nullopt;
}

[[nodiscard]] bool disjoint(const RealBounds& lhs, const RealBounds& rhs) {
    return lhs.upper < rhs.lower || rhs.upper < lhs.lower;
}

[[nodiscard]] bool samePoint(const RealBounds& lhs, const RealBounds& rhs) {
    return lhs.lower == lhs.upper
        && rhs.lower == rhs.upper
        && lhs.lower == rhs.lower;
}

} // namespace

Expr evaluateComparison(
    const Symbol& head,
    std::span<const Expr> arguments) {
    requireBinary(arguments, head.view());

    const Expr& lhs = arguments[0];
    const Expr& rhs = arguments[1];

    if (head.view() == names::equal || head.view() == names::notEqual) {
        bool determined = false;
        bool equal = false;

        if (lhs == rhs) {
            determined = true;
            equal = true;
        }
        else if (lhs.isNumber() && rhs.isNumber()) {
            determined = true;
            equal = lhs.asNumber() == rhs.asNumber();
        }
        else if (lhs.isBoolean() && rhs.isBoolean()) {
            determined = true;
            equal = lhs.asBoolean() == rhs.asBoolean();
        }
        else if (lhs.isString() && rhs.isString()) {
            determined = true;
            equal = lhs.asString() == rhs.asString();
        }
        else {
            const auto leftBounds = complexInformationBounds(lhs);
            const auto rightBounds = complexInformationBounds(rhs);
            if (leftBounds && rightBounds) {
                if (disjoint(leftBounds->real, rightBounds->real)
                    || disjoint(leftBounds->imaginary, rightBounds->imaginary)) {
                    determined = true;
                    equal = false;
                }
                else if (samePoint(leftBounds->real, rightBounds->real)
                    && samePoint(leftBounds->imaginary, rightBounds->imaginary)) {
                    determined = true;
                    equal = true;
                }
            }
        }

        // InformationEnclosureだけで証明できない近似値や記号式の不一致を、誤ってfalseとは断定しない。
        if (!determined)
            return unresolved(head, lhs, rhs);

        return Expr{head.view() == names::equal ? equal : !equal};
    }

    const auto left = informationBounds(lhs);
    const auto right = informationBounds(rhs);
    if (!left || !right)
        return unresolved(head, lhs, rhs);

    if (head.view() == names::less) {
        if (left->upper < right->lower) return Expr{true};
        if (left->lower >= right->upper) return Expr{false};
        return unresolved(head, lhs, rhs);
    }
    if (head.view() == names::lessEqual) {
        if (left->upper <= right->lower) return Expr{true};
        if (left->lower > right->upper) return Expr{false};
        return unresolved(head, lhs, rhs);
    }
    if (head.view() == names::greater) {
        if (left->lower > right->upper) return Expr{true};
        if (left->upper <= right->lower) return Expr{false};
        return unresolved(head, lhs, rhs);
    }
    if (head.view() == names::greaterEqual) {
        if (left->lower >= right->upper) return Expr{true};
        if (left->upper < right->lower) return Expr{false};
        return unresolved(head, lhs, rhs);
    }

    error::throwCalcError(
        error::CalcErrorType::Internal,
        "Unknown comparison operator");
}

Expr evaluateLogicalAnd(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    std::vector<Expr> unresolvedArguments;

    for (const Expr& argument : arguments) {
        if (!argument.isBoolean()) {
            unresolvedArguments.push_back(argument);
            continue;
        }

        if (!argument.asBoolean())
            return Expr{false};
    }

    if (unresolvedArguments.empty())
        return Expr{true};
    if (unresolvedArguments.size() == 1)
        return unresolvedArguments.front();

    return Expr::call(registry.symbol(evaluation::BuiltinId::LogicalAnd), std::move(unresolvedArguments));
}

Expr evaluateElement(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireBinary(arguments, names::element);
    if (!arguments[1].isSymbol())
        error::throwCalcError(
            error::CalcErrorType::Type,
            "element domain must be Integer, Rational, Real, or Complex");

    mathematics::NumericDomain domain = mathematics::NumericDomain::Unknown;
    const std::string_view name = arguments[1].asSymbol().view();
    if (name == "Integer") domain = mathematics::NumericDomain::Integer;
    else if (name == "Rational") domain = mathematics::NumericDomain::Rational;
    else if (name == "Real") domain = mathematics::NumericDomain::Real;
    else if (name == "Complex") domain = mathematics::NumericDomain::Complex;
    else
        error::throwCalcError(
            error::CalcErrorType::Type,
            "element domain must be Integer, Rational, Real, or Complex");

    const mathematics::AssumptionSet assumptions;
    const mathematics::KnowledgeContext knowledge{registry, mathematics, assumptions};
    const mathematics::TruthValue truth = knowledge.prove(
        mathematics::elementOf(arguments[0], domain));
    if (truth == mathematics::TruthValue::True)
        return Expr{true};
    if (truth == mathematics::TruthValue::False)
        return Expr{false};
    return Expr::call(registry.symbol(evaluation::BuiltinId::Element),
        {arguments[0], arguments[1]});
}

} // namespace mmcal::builtins
