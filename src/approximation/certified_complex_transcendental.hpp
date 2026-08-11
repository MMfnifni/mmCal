#pragma once

#include "complex_interval.hpp"

#include <cstddef>

namespace mmcal::approximation {

struct CertifiedArgumentResult final {
    RealInterval interval;
    std::size_t termsUsed = 0;
};

struct CertifiedComplexTranscendentalResult final {
    ComplexInterval interval;
    std::size_t termsUsed = 0;
};

// principal Arg(z) をラジアンで囲う。値域は (-Pi, Pi]。
// 負実軸は +Pi 側を採る。入力長方形がbranch cutを跨ぐ場合は、一意な連続区間で表現するため必要に応じて [-Pi,Pi] まで広げる。
[[nodiscard]] CertifiedArgumentResult enclosePrincipalArgument(
    const ComplexInterval& value,
    std::size_t precisionBits);

// principal Log(z) = ln|z| + I Arg(z)。z=0は定義域外。
[[nodiscard]] CertifiedComplexTranscendentalResult enclosePrincipalComplexLog(
    const ComplexInterval& value,
    std::size_t precisionBits);

// Exp(x+yI) = exp(x)(cos(y)+I sin(y))。ここでyは数学上のラジアン値であり、ユーザー向けAngleSemanticsは介在しない。
[[nodiscard]] CertifiedComplexTranscendentalResult encloseComplexExp(
    const ComplexInterval& value,
    std::size_t precisionBits);

// principal Power(z,w) := Exp(w * principal Log(z))。z=0の特殊規則は上位exact evaluatorで処理し、ここは非零baseを前提とする。
[[nodiscard]] CertifiedComplexTranscendentalResult enclosePrincipalPower(
    const ComplexInterval& base,
    const ComplexInterval& exponent,
    std::size_t precisionBits);

} // namespace mmcal::approximation
