#pragma once

#include "real_interval.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/rational.hpp"
#include "numeric/real_number.hpp"

#include <cstddef>

namespace mmcal::approximation {

// sin/cos/tan の証明付き包含区間。
struct CertifiedTrigEnclosure final {
    RealInterval interval;
    std::size_t termsUsed = 0;
    std::size_t precisionBits = 0;
};

// 最終的に10進丸めまで確定した結果。
struct CertifiedTrigResult final {
    numeric::DecimalApproximation value;
    std::size_t termsUsed = 0;
};

// exactなRationalラジアン値をTaylorの剰余上界で包含する。
// reference実装として、argument reductionを必要としないため任意の有限Rationalに使える。
[[nodiscard]] CertifiedTrigEnclosure encloseSinRadian(
    const numeric::Rational& argument,
    std::size_t precisionBits);

[[nodiscard]] CertifiedTrigEnclosure encloseCosRadian(
    const numeric::Rational& argument,
    std::size_t precisionBits);

[[nodiscard]] CertifiedTrigEnclosure encloseTanRadian(
    const numeric::Rational& argument,
    std::size_t precisionBits);

// 入力ラジアン自身がPi等に由来するcertified区間である場合のAPI。
// 区間中点をTaylor評価し、|sin'|,|cos'|<=1 によるLipschitz誤差を加える。
[[nodiscard]] CertifiedTrigEnclosure encloseSinRadianInterval(
    const RealInterval& argument,
    std::size_t precisionBits);

[[nodiscard]] CertifiedTrigEnclosure encloseCosRadianInterval(
    const RealInterval& argument,
    std::size_t precisionBits);

[[nodiscard]] CertifiedTrigEnclosure encloseTanRadianInterval(
    const RealInterval& argument,
    std::size_t precisionBits);

// exactなturn値を周期・象限対称性で縮約し、certified Piを通してラジアンへ変換する。
// Degree/Gradおよび q*Pi Rad の数値評価はこの経路を使う。
[[nodiscard]] CertifiedTrigEnclosure encloseSinTurns(
    const numeric::Rational& turns,
    std::size_t precisionBits);

[[nodiscard]] CertifiedTrigEnclosure encloseCosTurns(
    const numeric::Rational& turns,
    std::size_t precisionBits);

[[nodiscard]] CertifiedTrigEnclosure encloseTanTurns(
    const numeric::Rational& turns,
    std::size_t precisionBits);

// 明示ラジアン用の既存API。
[[nodiscard]] CertifiedTrigResult approximateSin(
    const numeric::RealNumber& argument,
    std::size_t fractionalDigits);

[[nodiscard]] CertifiedTrigResult approximateCos(
    const numeric::RealNumber& argument,
    std::size_t fractionalDigits);

[[nodiscard]] CertifiedTrigResult approximateTan(
    const numeric::RealNumber& argument,
    std::size_t fractionalDigits);

// Degree/Grad/default-angleや q*Pi Rad をturnとしてexactに表せる場合のAPI。
[[nodiscard]] CertifiedTrigResult approximateSinTurns(
    const numeric::Rational& turns,
    std::size_t fractionalDigits);

[[nodiscard]] CertifiedTrigResult approximateCosTurns(
    const numeric::Rational& turns,
    std::size_t fractionalDigits);

[[nodiscard]] CertifiedTrigResult approximateTanTurns(
    const numeric::Rational& turns,
    std::size_t fractionalDigits);

} // namespace mmcal::approximation
