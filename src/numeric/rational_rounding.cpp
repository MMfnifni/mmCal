#include "rational_rounding.hpp"

#include "integer_algorithms.hpp"

#include <compare>

namespace mmcal::numeric {

BigInt truncateToInteger(const Rational& value) {
    return value.numerator() / value.denominator();
}

BigInt floorToInteger(const Rational& value) {
    auto result = divmod(value.numerator(), value.denominator());
    if (!result.remainder.isZero() && value.numerator().isNegative())
        result.quotient -= BigInt{1};
    return result.quotient;
}

BigInt ceilToInteger(const Rational& value) {
    auto result = divmod(value.numerator(), value.denominator());
    if (!result.remainder.isZero() && value.numerator().isPositive())
        result.quotient += BigInt{1};
    return result.quotient;
}

BigInt roundToNearestEvenInteger(const Rational& value) {
    auto result = divmod(value.numerator(), value.denominator());
    if (result.remainder.isZero())
        return result.quotient;

    const BigInt twiceRemainder = result.remainder.abs() * BigInt{2};
    const auto comparison = twiceRemainder <=> value.denominator();
    bool awayFromZero = comparison == std::strong_ordering::greater;
    if (comparison == std::strong_ordering::equal)
        awayFromZero = !(result.quotient.abs() % BigInt{2}).isZero();

    if (awayFromZero)
        result.quotient += value.numerator().isNegative() ? BigInt{-1} : BigInt{1};
    return result.quotient;
}

} // namespace mmcal::numeric
