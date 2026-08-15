//sum・prod・meanなどの集約
#include "aggregate.hpp"

#include "builtins/exact_operations.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/rational.hpp"

#include <optional>

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace mmcal::builtins {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

struct RealBounds final {
    numeric::Rational lower;
    numeric::Rational upper;
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

[[nodiscard]] std::vector<Expr> aggregateItems(
    std::span<const Expr> arguments,
    std::string_view name,
    bool allowEmpty) {
    if (arguments.size() == 1 && arguments.front().isArray())
        return arguments.front().asArray().materialize();

    if (!allowEmpty && arguments.empty())
        error::throwCalcError(error::CalcErrorType::Type,
            std::string{name} + " requires at least one value");

    for (const Expr& argument : arguments)
        if (argument.isArray())
            error::throwCalcError(error::CalcErrorType::Type,
                std::string{name} + " accepts either one array or scalar arguments");

    return {arguments.begin(), arguments.end()};
}

[[nodiscard]] Expr extrema(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    BuiltinId id,
    bool minimum) {
    std::vector<Expr> items = aggregateItems(arguments,
        minimum ? "min" : "max", false);
    if (items.empty())
        error::throwCalcError(error::CalcErrorType::Domain,
            minimum ? "min of an empty array is undefined" : "max of an empty array is undefined");

    std::vector<RealBounds> bounds;
    bounds.reserve(items.size());
    for (const Expr& item : items) {
        const auto valueBounds = informationBounds(item);
        if (!valueBounds)
            return Expr::call(registry.symbol(id), std::move(items));
        bounds.push_back(*valueBounds);
    }

    std::size_t best = 0;
    for (std::size_t i = 1; i < items.size(); ++i) {
        if (items[i] == items[best])
            continue;
        if (minimum) {
            if (bounds[i].upper < bounds[best].lower) {
                best = i;
                continue;
            }
            if (bounds[best].upper <= bounds[i].lower)
                continue;
        }
        else {
            if (bounds[i].lower > bounds[best].upper) {
                best = i;
                continue;
            }
            if (bounds[best].lower >= bounds[i].upper)
                continue;
        }
        return Expr::call(registry.symbol(id), std::move(items));
    }
    return items[best];
}

} // namespace

Expr evaluateSum(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> items = aggregateItems(arguments, "sum", true);
    if (items.empty())
        return integer(0);
    return exact::add(std::move(items), registry, mathematics, angles);
}

Expr evaluateProduct(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> items = aggregateItems(arguments, "prod", true);
    if (items.empty())
        return integer(1);
    return exact::multiply(std::move(items), registry, mathematics, angles);
}

Expr evaluateMin(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    return extrema(arguments, registry, BuiltinId::Min, true);
}

Expr evaluateMax(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    return extrema(arguments, registry, BuiltinId::Max, false);
}

Expr evaluateMean(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> items = aggregateItems(arguments, "mean", false);
    if (items.empty())
        error::throwCalcError(error::CalcErrorType::Domain, "mean of an empty array is undefined");
    const std::size_t count = items.size();
    Expr total = exact::add(std::move(items), registry, mathematics, angles);
    return exact::divide(std::move(total), Expr{Number{BigInt::parse(std::to_string(count))}},
        registry, mathematics, angles);
}

} // namespace mmcal::builtins
