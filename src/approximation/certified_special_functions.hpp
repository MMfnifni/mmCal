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

// Gauss 2F1。現backendはprincipal branchのうち |z|<1 のexact Rational pointを
// Gauss級数＋厳密な幾何tail boundで保証評価する。
[[nodiscard]] RealInterval encloseHypergeometric2F1Real(
    const numeric::Rational& a,
    const numeric::Rational& b,
    const numeric::Rational& c,
    const numeric::Rational& z,
    std::size_t precisionBits);

// Legendre不完全楕円積分。amplitudeは常にRadian。
// 現backendは |m|<1（Piは加えて|n|<1）のexact Rational pointを
// m/n級数とsin偶数冪積分の漸化式で保証評価する。
[[nodiscard]] RealInterval encloseEllipticFReal(
    const numeric::Rational& phi,
    const numeric::Rational& m,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseEllipticEReal(
    const numeric::Rational& phi,
    const numeric::Rational& m,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseEllipticPiReal(
    const numeric::Rational& n,
    const numeric::Rational& phi,
    const numeric::Rational& m,
    std::size_t precisionBits);

// 古典的積分函数。実数backendではprincipal real valueを保証区間で返す。
[[nodiscard]] RealInterval encloseExponentialIntegralEiReal(
    const RealInterval& input,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseSineIntegralSiReal(
    const RealInterval& input,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseCosineIntegralCiPositive(
    const RealInterval& input,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseLogarithmicIntegralLiPositive(
    const RealInterval& input,
    std::size_t precisionBits);

// 現backendは正整数sと|z|<1のexact Rational pointを級数で保証評価する。
[[nodiscard]] RealInterval enclosePolylogReal(
    std::uint64_t order,
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
