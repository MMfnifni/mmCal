// 実区間
#include "real_interval.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rounding_mode.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <utility>

namespace mmcal::approximation {
namespace {

using numeric::BigFloat;
using numeric::RoundingMode;

[[nodiscard]] BigFloat minimum(std::array<BigFloat, 4>& values) {
    return *std::min_element(values.begin(), values.end());
}

[[nodiscard]] BigFloat maximum(std::array<BigFloat, 4>& values) {
    return *std::max_element(values.begin(), values.end());
}

} // namespace

RealInterval::RealInterval(BigFloat lower, BigFloat upper)
    : lower_(std::move(lower)), upper_(std::move(upper)) {
    if (lower_ > upper_)
        throw std::invalid_argument("RealInterval lower bound exceeds upper bound");
}

RealInterval RealInterval::point(BigFloat value) {
    // 同一objectを同一初期化式でcopy/moveしない。MSVCを含め評価順の差に依存せず、
    // point intervalは必ず同一の2端点から構築する。
    BigFloat lower = value;
    BigFloat upper = value;
    return RealInterval{std::move(lower), std::move(upper)};
}

RealInterval RealInterval::fromRational(
    const numeric::Rational& value,
    std::size_t precisionBits) {
    return RealInterval{
        BigFloat::fromRational(value, precisionBits, RoundingMode::TowardNegative),
        BigFloat::fromRational(value, precisionBits, RoundingMode::TowardPositive)
    };
}

RealInterval RealInterval::fromRationalBounds(
    const numeric::Rational& lower,
    const numeric::Rational& upper,
    std::size_t precisionBits) {
    if (lower > upper)
        throw std::invalid_argument("RealInterval rational bounds are reversed");

    return RealInterval{
        BigFloat::fromRational(lower, precisionBits, RoundingMode::TowardNegative),
        BigFloat::fromRational(upper, precisionBits, RoundingMode::TowardPositive)
    };
}

const BigFloat& RealInterval::lower() const noexcept {
    return lower_;
}

const BigFloat& RealInterval::upper() const noexcept {
    return upper_;
}

bool RealInterval::isPoint() const noexcept {
    return lower_ == upper_;
}

bool RealInterval::containsZero() const noexcept {
    const BigFloat zero;
    return lower_ <= zero && zero <= upper_;
}

bool RealInterval::contains(const numeric::Rational& value) const {
    return lower_.toRational() <= value && value <= upper_.toRational();
}

RealInterval RealInterval::roundedOutward(std::size_t precisionBits) const {
    return RealInterval{
        lower_.rounded(precisionBits, RoundingMode::TowardNegative),
        upper_.rounded(precisionBits, RoundingMode::TowardPositive)
    };
}

RealInterval hull(const RealInterval& lhs, const RealInterval& rhs) {
    const BigFloat lower = lhs.lower() < rhs.lower() ? lhs.lower() : rhs.lower();
    const BigFloat upper = lhs.upper() > rhs.upper() ? lhs.upper() : rhs.upper();
    return RealInterval{lower, upper};
}

RealInterval negate(const RealInterval& value) {
    return RealInterval{-value.upper(), -value.lower()};
}

RealInterval add(
    const RealInterval& lhs,
    const RealInterval& rhs,
    std::size_t precisionBits) {
    return RealInterval{
        numeric::add(
            lhs.lower(), rhs.lower(), precisionBits, RoundingMode::TowardNegative),
        numeric::add(
            lhs.upper(), rhs.upper(), precisionBits, RoundingMode::TowardPositive)
    };
}

RealInterval subtract(
    const RealInterval& lhs,
    const RealInterval& rhs,
    std::size_t precisionBits) {
    return RealInterval{
        numeric::subtract(
            lhs.lower(), rhs.upper(), precisionBits, RoundingMode::TowardNegative),
        numeric::subtract(
            lhs.upper(), rhs.lower(), precisionBits, RoundingMode::TowardPositive)
    };
}

RealInterval multiply(
    const RealInterval& lhs,
    const RealInterval& rhs,
    std::size_t precisionBits) {
    // 積 xy は長方形 [a,b]x[c,d] 上の双線形函数なので、極値は4隅のいずれかにある。
    // 各隅を下向き/上向きに別々に丸め、その最小/最大を採ることで包含を保証する。
    std::array<BigFloat, 4> lowerProducts{
        numeric::multiply(lhs.lower(), rhs.lower(), precisionBits, RoundingMode::TowardNegative),
        numeric::multiply(lhs.lower(), rhs.upper(), precisionBits, RoundingMode::TowardNegative),
        numeric::multiply(lhs.upper(), rhs.lower(), precisionBits, RoundingMode::TowardNegative),
        numeric::multiply(lhs.upper(), rhs.upper(), precisionBits, RoundingMode::TowardNegative)
    };
    std::array<BigFloat, 4> upperProducts{
        numeric::multiply(lhs.lower(), rhs.lower(), precisionBits, RoundingMode::TowardPositive),
        numeric::multiply(lhs.lower(), rhs.upper(), precisionBits, RoundingMode::TowardPositive),
        numeric::multiply(lhs.upper(), rhs.lower(), precisionBits, RoundingMode::TowardPositive),
        numeric::multiply(lhs.upper(), rhs.upper(), precisionBits, RoundingMode::TowardPositive)
    };

    return RealInterval{minimum(lowerProducts), maximum(upperProducts)};
}

RealInterval divide(
    const RealInterval& lhs,
    const RealInterval& rhs,
    std::size_t precisionBits) {
    if (rhs.containsZero())
        throw std::domain_error("RealInterval division by an interval containing zero");

    // 分母区間が0を跨がなければ x/y は長方形上で連続であり、極値は4隅で得られる。
    std::array<BigFloat, 4> lowerQuotients{
        numeric::divide(lhs.lower(), rhs.lower(), precisionBits, RoundingMode::TowardNegative),
        numeric::divide(lhs.lower(), rhs.upper(), precisionBits, RoundingMode::TowardNegative),
        numeric::divide(lhs.upper(), rhs.lower(), precisionBits, RoundingMode::TowardNegative),
        numeric::divide(lhs.upper(), rhs.upper(), precisionBits, RoundingMode::TowardNegative)
    };
    std::array<BigFloat, 4> upperQuotients{
        numeric::divide(lhs.lower(), rhs.lower(), precisionBits, RoundingMode::TowardPositive),
        numeric::divide(lhs.lower(), rhs.upper(), precisionBits, RoundingMode::TowardPositive),
        numeric::divide(lhs.upper(), rhs.lower(), precisionBits, RoundingMode::TowardPositive),
        numeric::divide(lhs.upper(), rhs.upper(), precisionBits, RoundingMode::TowardPositive)
    };

    return RealInterval{minimum(lowerQuotients), maximum(upperQuotients)};
}

} // namespace mmcal::approximation
