#pragma once

#include "approximation/approximation_context.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "linear_algebra/exact_matrix.hpp"
#include "linear_algebra/matrix.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <optional>

namespace mmcal::linear_algebra {

// 固有値。exactでは三角行列と2x2を扱い，一般行列はprecision-aware Schur backendへ委ねる。
[[nodiscard]] std::optional<expression::Expr> eigenvalues(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

// 固有vectorを列に並べたn×n行列。exactは対角行列のみ。
[[nodiscard]] std::optional<expression::Expr> eigenvectors(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

// {values, vectors}。vectorsの各列が対応する固有vector。
[[nodiscard]] std::optional<expression::Expr> eigensystem(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

// 一般実/複素正方行列のprecision-aware eigen backend。
// Hessenberg化 + implicit shifted complex QRでSchur形を作り，必要なら三角back substitutionで
// eigenvectorをmaterializeする。非正規行列では各rootの包含区間ではなく，Schur/eigenpair関係を監査する。
[[nodiscard]] std::optional<expression::Expr> approximateEigenvalues(
    const expression::ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateEigenvectors(
    const expression::ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateEigensystem(
    const expression::ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

} // namespace mmcal::linear_algebra
