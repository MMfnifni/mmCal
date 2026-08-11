#pragma once

#include "complex_interval.hpp"
#include "real_interval.hpp"

#include <cstddef>

namespace mmcal::approximation {

// ラジアン引数の複素三角函数。ユーザー向け角度単位変換はCertifiedEvaluator側で行う。
[[nodiscard]] ComplexInterval encloseComplexSinRadian(
    const ComplexInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseComplexCosRadian(
    const ComplexInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseComplexTanRadian(
    const ComplexInterval& value, std::size_t precisionBits);

// principal inverse trig。返り値は数学上のラジアン。
[[nodiscard]] RealInterval encloseAsinRealRadian(
    const RealInterval& value, std::size_t precisionBits);
[[nodiscard]] RealInterval encloseAcosRealRadian(
    const RealInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval enclosePrincipalComplexAsinRadian(
    const ComplexInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval enclosePrincipalComplexAcosRadian(
    const ComplexInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval enclosePrincipalComplexAtanRadian(
    const ComplexInterval& value, std::size_t precisionBits);

// 実立方根。real cube rootなので負実数も実数へ写す。
[[nodiscard]] RealInterval encloseRealCubeRoot(
    const RealInterval& value, std::size_t precisionBits);

// 双曲線函数。実数専用経路は虚部0を厳密に保つ。
[[nodiscard]] RealInterval encloseSinhReal(
    const RealInterval& value, std::size_t precisionBits);
[[nodiscard]] RealInterval encloseCoshReal(
    const RealInterval& value, std::size_t precisionBits);
[[nodiscard]] RealInterval encloseTanhReal(
    const RealInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseComplexSinh(
    const ComplexInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseComplexCosh(
    const ComplexInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseComplexTanh(
    const ComplexInterval& value, std::size_t precisionBits);

[[nodiscard]] RealInterval encloseAsinhReal(
    const RealInterval& value, std::size_t precisionBits);
[[nodiscard]] RealInterval encloseAcoshReal(
    const RealInterval& value, std::size_t precisionBits);
[[nodiscard]] RealInterval encloseAtanhReal(
    const RealInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval enclosePrincipalComplexAsinh(
    const ComplexInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval enclosePrincipalComplexAcosh(
    const ComplexInterval& value, std::size_t precisionBits);
[[nodiscard]] ComplexInterval enclosePrincipalComplexAtanh(
    const ComplexInterval& value, std::size_t precisionBits);

} // namespace mmcal::approximation
