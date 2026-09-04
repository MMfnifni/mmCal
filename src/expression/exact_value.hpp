#pragma once

#include "expression/expr.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"

#include <charconv>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <system_error>
#include <utility>

namespace mmcal::expression::exact {

// exact scalarのAST化と判定を一箇所へ寄せる。
// 各symbolic subsystemがNumberの内部表現へ個別依存すると，Integer/Rationalのpromotion規則がずれやすい。
[[nodiscard]] inline Expr integer(std::int64_t value) {
    return Expr{numeric::Number{numeric::BigInt{value}}};
}

[[nodiscard]] inline Expr integer(numeric::BigInt value) {
    return Expr{numeric::Number{std::move(value)}};
}

[[nodiscard]] inline Expr rational(numeric::Rational value) {
    return Expr{numeric::Number{std::move(value)}};
}

[[nodiscard]] inline Expr rational(std::int64_t numerator, std::int64_t denominator) {
    return rational(numeric::Rational{numeric::BigInt{numerator}, numeric::BigInt{denominator}});
}

[[nodiscard]] inline bool isZero(const Expr& expression) noexcept {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] inline bool isOne(const Expr& expression) {
    return expression.isNumber()
        && expression.asNumber().isReal()
        && expression.asNumber().asReal().toRational()
            == numeric::Rational{numeric::BigInt{1}};
}

[[nodiscard]] inline std::optional<numeric::Rational> realRational(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    return expression.asNumber().asReal().toRational();
}

[[nodiscard]] inline std::optional<std::size_t> positiveSize(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    const numeric::BigInt& value = expression.asNumber().asReal().asInteger();
    if (value.isNegative() || value.isZero())
        return std::nullopt;

    // BigIntからsize_tへの暗黙の切詰めは許さず，桁数指定などの上限超過も失敗として扱う。
    const std::string text = value.toString();
    std::size_t result = 0;
    const auto converted = std::from_chars(text.data(), text.data() + text.size(), result);
    if (converted.ec != std::errc{} || converted.ptr != text.data() + text.size())
        return std::nullopt;
    return result;
}

} // namespace mmcal::expression::exact
