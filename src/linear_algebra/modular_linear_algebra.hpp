#pragma once

#include "linear_algebra/fraction_free_elimination.hpp"
#include "numeric/rational.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace mmcal::linear_algebra {

// word-sized prime fieldを使うexact整数線形代数の実行統計。
// algorithm選択とbenchmark監査用であり，数学的結果の一部ではない。
struct ModularLinearAlgebraStats final {
    std::size_t primesTried = 0;
    std::size_t primesAccepted = 0;
    std::size_t reconstructedModulusBits = 0;
};

// Bareissとmodular backendのautomatic dispatch policy。
// 小行列では変換・CRT overheadを避け，大きい次数または高い係数heightでmodularへ送る。
[[nodiscard]] bool preferModularDeterminant(const IntegerMatrixBuffer& matrix) noexcept;
[[nodiscard]] bool preferModularSolve(
    const IntegerMatrixBuffer& augmented,
    std::size_t variables) noexcept;
[[nodiscard]] bool preferModularInverse(const IntegerMatrixBuffer& matrix) noexcept;

// det(A) mod pを複数の31-bit primeで計算し，Hadamard上界を越えるまでCRT再構成する。
[[nodiscard]] numeric::BigInt modularDeterminant(
    const IntegerMatrixBuffer& matrix,
    ModularLinearAlgebraStats* stats = nullptr);

// [A|B]を有限体上で解き，CRT + rational reconstruction後に元の整数系でexact verificationする。
// unique solutionを再構成できた場合だけvariables x rhsColumnsのrow-major Rationalを返す。
// singular/inconsistent/bad-primeのみのときはnulloptとし，呼出側がBareissへfallbackする。
[[nodiscard]] std::optional<std::vector<numeric::Rational>> modularSolve(
    const IntegerMatrixBuffer& augmented,
    std::size_t variables,
    ModularLinearAlgebraStats* stats = nullptr);

// A^-1 = adj(A)/det(A)を利用するinverse専用backend。detとadjugateを整数CRTで
// 再構成するため，n^2個の独立したrational reconstructionを避ける。
[[nodiscard]] std::optional<std::vector<numeric::Rational>> modularInverse(
    const IntegerMatrixBuffer& matrix,
    ModularLinearAlgebraStats* stats = nullptr);

// benchmarkとdispatch監査用のrigorous reconstruction bound。
[[nodiscard]] std::size_t determinantReconstructionBits(
    const IntegerMatrixBuffer& matrix) noexcept;
[[nodiscard]] std::size_t solutionReconstructionBits(
    const IntegerMatrixBuffer& augmented,
    std::size_t variables) noexcept;

} // namespace mmcal::linear_algebra
