// 組合せ・Fibonacci
#include "combinatorics.hpp"
#include "builtin_helpers.hpp"

#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"

#include <algorithm>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>

namespace mmcal::builtins {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(BigInt value) {
    return Expr{Number{std::move(value)}};
}

[[nodiscard]] const BigInt* exactInteger(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return nullptr;
    return &expression.asNumber().asReal().asInteger();
}

[[nodiscard]] Expr hold(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    return Expr::call(registry.symbol(id), {arguments.begin(), arguments.end()});
}

[[nodiscard]] std::uint64_t requireNonnegativeCount(
    const Expr& expression,
    std::string_view function,
    std::string_view argumentName) {
    const BigInt* integerValue = exactInteger(expression);
    if (!integerValue)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{function} + " requires integer arguments");
    if (integerValue->isNegative())
        error::throwCalcError(
            error::CalcErrorType::Domain,
            std::string{function} + " requires non-negative " + std::string{argumentName});

    const auto converted = numeric::tryToUint64(*integerValue);
    if (!converted)
        error::throwCalcError(
            error::CalcErrorType::Overflow,
            std::string{function} + " argument is too large for exact evaluation");
    return *converted;
}

[[nodiscard]] BigInt permutation(std::uint64_t n, std::uint64_t r) {
    if (r > n)
        return BigInt{};

    BigInt result{1};
    for (std::uint64_t i = 0; i < r; ++i)
        result *= BigInt::parse(std::to_string(n - i));
    return result;
}

[[nodiscard]] BigInt combination(std::uint64_t n, std::uint64_t r) {
    if (r > n)
        return BigInt{};

    r = std::min(r, n - r);
    BigInt result{1};
    for (std::uint64_t i = 1; i <= r; ++i) {
        result *= BigInt::parse(std::to_string(n - r + i));
        result /= BigInt::parse(std::to_string(i));
    }
    return result;
}

[[nodiscard]] std::pair<BigInt, BigInt> fibonacciPair(std::uint64_t n) {
    if (n == 0)
        return {BigInt{}, BigInt{1}};

    auto [a, b] = fibonacciPair(n >> 1U);
    const BigInt c = a * (BigInt{2} * b - a);
    const BigInt d = a * a + b * b;
    if ((n & 1U) == 0)
        return {c, d};
    return {d, c + d};
}

} // namespace

Expr evaluatePermutation(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::permutation);
    if (!exactInteger(arguments[0]) || !exactInteger(arguments[1])) {
        if (!arguments[0].isNumber() || !arguments[1].isNumber())
            return hold(BuiltinId::Permutation, arguments, registry);
        error::throwCalcError(error::CalcErrorType::Type, "perm requires integer arguments");
    }

    const std::uint64_t n = requireNonnegativeCount(arguments[0], "perm", "n");
    const std::uint64_t r = requireNonnegativeCount(arguments[1], "perm", "r");
    return integer(permutation(n, r));
}

Expr evaluateCombination(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::combination);
    if (!exactInteger(arguments[0]) || !exactInteger(arguments[1])) {
        if (!arguments[0].isNumber() || !arguments[1].isNumber())
            return hold(BuiltinId::Combination, arguments, registry);
        error::throwCalcError(error::CalcErrorType::Type, "comb requires integer arguments");
    }

    const std::uint64_t n = requireNonnegativeCount(arguments[0], "comb", "n");
    const std::uint64_t r = requireNonnegativeCount(arguments[1], "comb", "r");
    return integer(combination(n, r));
}

Expr evaluateFibonacci(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::fibonacci);
    if (!exactInteger(arguments.front())) {
        if (!arguments.front().isNumber())
            return hold(BuiltinId::Fibonacci, arguments, registry);
        error::throwCalcError(error::CalcErrorType::Type, "fib requires an integer argument");
    }

    const std::uint64_t n = requireNonnegativeCount(arguments.front(), "fib", "n");
    return integer(fibonacciPair(n).first);
}

} // namespace mmcal::builtins
