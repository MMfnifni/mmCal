#pragma once

#include "real_interval.hpp"
#include "numeric/rational.hpp"

#include <cstddef>

namespace mmcal::approximation {

// 実Gamma函数の包含区間。非正整数poleを含む入力区間は評価しない。
[[nodiscard]] RealInterval encloseGammaReal(
    const RealInterval& input,
    std::size_t precisionBits);

// C/POSIXのlgammaと同じく、実軸上の log(|Gamma(x)|)。
[[nodiscard]] RealInterval encloseLogGammaReal(
    const RealInterval& input,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseErfReal(
    const RealInterval& input,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseErfcReal(
    const RealInterval& input,
    std::size_t precisionBits);

// Fresnel C(x)=integral_0^x cos(Pi t^2/2)dt, S(x)=integral_0^x sin(Pi t^2/2)dt。
// 実軸ではTaylor/部分積分漸近展開を切り替え、剰余を明示的に包含する。
[[nodiscard]] RealInterval encloseFresnelCReal(
    const RealInterval& input,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseFresnelSReal(
    const RealInterval& input,
    std::size_t precisionBits);

// Kummerの合流型超幾何函数 M(a,b,z)=1F1(a;b;z)。
// 現backendはexact Rationalのpoint引数を級数＋厳密tail boundで保証評価する。
[[nodiscard]] RealInterval encloseHypergeometric1F1Real(
    const numeric::Rational& a,
    const numeric::Rational& b,
    const numeric::Rational& z,
    std::size_t precisionBits);

// a,b>0 に対するBetaとlog Beta。Gammaの比ではなくlog-domainで評価する。
[[nodiscard]] RealInterval encloseBetaPositive(
    const RealInterval& a,
    const RealInterval& b,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseBetaLogPositive(
    const RealInterval& a,
    const RealInterval& b,
    std::size_t precisionBits);

} // namespace mmcal::approximation
