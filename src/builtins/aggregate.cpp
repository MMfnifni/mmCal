//sum・prod・meanなどの集約
#include "aggregate.hpp"

#include "builtins/exact_operations.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

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

    bool allExactReal = true;
    for (const Expr& item : items)
        allExactReal = allExactReal && item.isNumber() && item.asNumber().isReal();

    if (!allExactReal)
        return Expr::call(registry.symbol(id), std::move(items));

    Expr result = items.front();
    for (std::size_t i = 1; i < items.size(); ++i) {
        const auto ordering = items[i].asNumber().asReal() <=> result.asNumber().asReal();
        if ((minimum && ordering < 0) || (!minimum && ordering > 0))
            result = items[i];
    }
    return result;
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
