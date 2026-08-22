// 実平方根の保証付き評価
#include "certified_sqrt.hpp"
#include "certified_precision.hpp"

#include "approximation_context.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "numeric/detail/binary_scale.hpp"
#include "numeric/integer_algorithms.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace mmcal::approximation {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::DecimalApproximation;
using numeric::Rational;
using numeric::RealNumber;
using numeric::RoundingMode;


[[nodiscard]] std::size_t nextGuardDigits(std::size_t guardDigits) {
    const std::size_t growth = std::max<std::size_t>(8, guardDigits / 2);
    return checkedPrecisionAdd(guardDigits, growth, "Certified sqrt precision is too large");
}

[[nodiscard]] std::optional<DecimalApproximation> tryCertifiedDecimal(
    const RealInterval& interval,
    std::size_t fractionalDigits) {
    return DecimalApproximation::fromCertifiedInterval(
        interval.lower().toRational(),
        interval.upper().toRational(),
        fractionalDigits);
}

[[nodiscard]] std::int64_t floorHalf(std::int64_t value) noexcept {
    // C++の整数除算は0方向へ丸める。負の奇数だけ1小さい整数へ補正すればfloor(value/2)。
    std::int64_t result = value / 2;
    if (value < 0 && value % 2 != 0)
        --result;
    return result;
}

[[nodiscard]] std::size_t positiveShiftMagnitude(std::int64_t value) {
    // value<0 が前提。INT64_MINでも符号付きnegateを行わない。
    const std::uint64_t magnitude = std::uint64_t{0} - static_cast<std::uint64_t>(value);
    if (magnitude > std::numeric_limits<std::size_t>::max())
        throw std::length_error("Certified sqrt shift distance is too large");
    return static_cast<std::size_t>(magnitude);
}

[[nodiscard]] BigFloat::exponent_type subtractScaleBits(
    std::int64_t exponent,
    std::size_t scaleBits) {
    if (scaleBits > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::overflow_error("Certified sqrt exponent underflow");

    const auto scale = static_cast<std::int64_t>(scaleBits);
    if (exponent < std::numeric_limits<std::int64_t>::min() + scale)
        throw std::overflow_error("Certified sqrt exponent underflow");
    return exponent - scale;
}

} // namespace

CertifiedSqrtEnclosure encloseSqrt(
    const Rational& value,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified sqrt precision must be at least one bit");
    if (value.numerator().isNegative())
        throw std::domain_error("Certified real sqrt requires a non-negative value");

    if (value.isZero()) {
        const BigFloat zero = BigFloat::fromBigInt(BigInt{}, precisionBits, RoundingMode::TowardZero);
        return CertifiedSqrtEnclosure{RealInterval::point(zero), precisionBits};
    }

    // -------------------------------------------------------------------------
    // 1. まず x を 4^k * z, 1 <= z < 4 の形へexactに正規化する。
    // -------------------------------------------------------------------------
    // e = floor(log2(x)) とすると k = floor(e/2) で
    //     z = x / 2^(2k)
    // は必ず 1 <= z < 4 に入る。したがって sqrt(z) は [1,2) にあり、「何bitが有効桁なのか」を扱いやすい。
    // log2自体もBigIntのbit長と比較だけで厳密に求め、machine floatは使わない。
    const std::int64_t binaryExponent = numeric::detail::floorLog2PositiveRatio(
        value.numerator(), value.denominator());
    const std::int64_t scaleExponent = floorHalf(binaryExponent);

    BigInt normalizedNumerator = value.numerator();
    BigInt normalizedDenominator = value.denominator();

    if (scaleExponent >= 0) {
        const auto shift64 = static_cast<std::uint64_t>(scaleExponent) * 2U;
        if (shift64 > std::numeric_limits<std::size_t>::max())
            throw std::length_error("Certified sqrt normalization shift is too large");
        normalizedDenominator <<= static_cast<std::size_t>(shift64);
    }
    else {
        // 2*scaleExponent の負の大きさを符号付きoverflowなしで求める。
        const std::size_t oneShift = positiveShiftMagnitude(scaleExponent);
        if (oneShift > std::numeric_limits<std::size_t>::max() / 2)
            throw std::length_error("Certified sqrt normalization shift is too large");
        normalizedNumerator <<= oneShift * 2;
    }

    // -------------------------------------------------------------------------
    // 2. sqrt(z) を固定小数点の整数平方根へ変換する。
    // -------------------------------------------------------------------------
    // S bitだけ小数部を持つ整数 m を
    //     m = floor(sqrt(z) * 2^S)
    // としたい。両辺を二乗すると
    //     m^2 <= z * 2^(2S)
    // なので、z = p/q に対して
    //     Q = floor(p * 2^(2S) / q)
    //     m = floor(sqrt(Q))
    // と完全な整数演算だけで求められる。
    // floor(sqrt(floor(A))) == floor(sqrt(A)) なので、途中のfloorで下界を失わない。
    //
    // precisionBitsより数bit余分にmを作り、最後にBigFloatのdirected roundingで下端/上端へ丸める。
    // 正しさはこの余分bit数そのものには依存せず、最終的なdecimal丸め一致までapproximateSqrt側が作業精度を増やす。
    const std::size_t fixedFractionBits = checkedPrecisionAdd(
        precisionBits, 4, "Certified sqrt precision is too large");
    if (fixedFractionBits > std::numeric_limits<std::size_t>::max() / 2)
        throw std::overflow_error("Certified sqrt precision is too large");

    const std::size_t doubledBits = fixedFractionBits * 2;
    const BigInt scaledNumerator = normalizedNumerator << doubledBits;
    const auto quotientAndRemainder = numeric::divmod(
        scaledNumerator, normalizedDenominator);
    const auto rootResult = numeric::integerSqrt(quotientAndRemainder.quotient);
    const BigInt& lowerInteger = rootResult.root;

    // m/2^S が真値そのものかをexactに判定する。
    // quotientの余りだけを見るのでは不十分なので、元のp/qへ戻して
    //     m^2 q == p 2^(2S)
    // をBigIntで直接比較する。
    const BigInt exactLeft = lowerInteger * lowerInteger * normalizedDenominator;
    const bool exactDyadicRoot = exactLeft == scaledNumerator;
    const BigInt upperInteger = exactDyadicRoot
        ? lowerInteger
        : lowerInteger + BigInt{1};

    // sqrt(x) = sqrt(z) * 2^k なので、固定小数点m/2^Sの最終指数は k-S。
    const auto outputExponent = subtractScaleBits(scaleExponent, fixedFractionBits);

    const BigFloat lower = BigFloat::fromDyadic(
        lowerInteger,
        outputExponent,
        precisionBits,
        RoundingMode::TowardNegative);
    const BigFloat upper = BigFloat::fromDyadic(
        upperInteger,
        outputExponent,
        precisionBits,
        RoundingMode::TowardPositive);

    return CertifiedSqrtEnclosure{RealInterval{lower, upper}, precisionBits};
}

CertifiedSqrtEnclosure encloseSqrt(
    const RealInterval& value,
    std::size_t precisionBits) {
    if (value.lower().isNegative())
        throw std::domain_error("Certified real sqrt interval cannot include negative values");

    // sqrtは[0,+infinity)で単調増加なので、入力区間[a,b]の像は厳密に[sqrt(a), sqrt(b)]。
    // 各端点はBigFloat自身がexactなdyadic Rationalなので、Rational版のcertified sqrtへそのまま渡せる。
    const CertifiedSqrtEnclosure lower = encloseSqrt(
        value.lower().toRational(), precisionBits);
    const CertifiedSqrtEnclosure upper = encloseSqrt(
        value.upper().toRational(), precisionBits);

    return CertifiedSqrtEnclosure{
        RealInterval{lower.interval.lower(), upper.interval.upper()},
        precisionBits
    };
}

DecimalApproximation approximateSqrt(
    const RealNumber& value,
    std::size_t fractionalDigits) {
    if (fractionalDigits == 0)
        throw std::invalid_argument("Certified sqrt precision must be greater than zero");
    if (value.isNegative())
        throw std::domain_error("Certified real sqrt requires a non-negative value");

    ApproximationContext context{fractionalDigits};
    const Rational exact = value.toRational();
    for (;;) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        const CertifiedSqrtEnclosure enclosure = encloseSqrt(
            exact, context.workingBinaryBits());
        if (const auto decimal = tryCertifiedDecimal(enclosure.interval, fractionalDigits))
            return *decimal;

        // 「何回Newtonしたか」ではなく、上下端の要求桁丸めが一致したかだけを最終certification条件にする。
        // 丸め境界近傍なら必要なだけprecisionを増やす。
        context.setGuardDigits(nextGuardDigits(context.guardDigits()));
    }
}

} // namespace mmcal::approximation
