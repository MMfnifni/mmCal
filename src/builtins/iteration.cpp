// exactな有限列生成
#include "iteration.hpp"

#include "error/error_message.hpp"
#include "expression/array_utils.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>

namespace mmcal::builtins {
namespace {

using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] std::optional<Rational> exactRealRational(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    return expression.asNumber().asReal().toRational();
}

[[nodiscard]] BigInt floorNonNegative(const Rational& value) {
    return value.numerator() / value.denominator();
}

[[nodiscard]] std::size_t checkedCount(const BigInt& count, std::string_view caller) {
    const auto converted = numeric::tryToUint64(count);
    if (!converted || *converted > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()))
        error::throwCalcError(
            error::CalcErrorType::Overflow,
            std::string{caller} + " result is too large");
    return static_cast<std::size_t>(*converted);
}

} // namespace

std::vector<Expr> exactRangeValues(
    std::span<const Expr> arguments,
    std::string_view caller) {
    if (arguments.empty() || arguments.size() > 3)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{caller} + " expects 1 to 3 range arguments");

    Expr defaultOne = integer(1);
    const Expr& startExpr = arguments.size() == 1 ? defaultOne : arguments[0];
    const Expr& endExpr = arguments.size() == 1 ? arguments[0] : arguments[1];
    const Expr& stepExpr = arguments.size() == 3 ? arguments[2] : defaultOne;

    const auto start = exactRealRational(startExpr);
    const auto end = exactRealRational(endExpr);
    const auto step = exactRealRational(stepExpr);
    if (!start || !end || !step)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{caller} + " requires exact real integer or rational bounds");
    if (step->isZero())
        error::throwCalcError(
            error::CalcErrorType::Domain,
            std::string{caller} + " step must be nonzero");

    const bool increasing = step->numerator().isPositive();
    if ((increasing && *start > *end) || (!increasing && *start < *end))
        return {};

    const Rational quotient = (*end - *start) / *step;
    if (quotient.numerator().isNegative())
        return {};

    const BigInt countBig = floorNonNegative(quotient) + BigInt{1};
    const std::size_t count = checkedCount(countBig, caller);

    std::vector<Expr> result;
    result.reserve(count);
    Rational current = *start;
    for (std::size_t i = 0; i < count; ++i) {
        result.emplace_back(Number{current});
        current += *step;
    }
    return result;
}

Expr evaluateRange(std::span<const Expr> arguments) {
    return expression::braceValue(exactRangeValues(arguments, "range"));
}

} // namespace mmcal::builtins
