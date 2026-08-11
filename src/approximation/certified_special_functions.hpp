#pragma once

#include "real_interval.hpp"

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
