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

// reduced SVD。k=min(m,n)としてU:m×k, S:k×k, V:n×kを一般braceで返す。
// 実数はA=U S Transpose[V]，複素数はA=U S conjugateTranspose[V]。
// exactは自然に閉じる実対角行列だけを扱い，一般行列はapproximate backendへ委ねる。
[[nodiscard]] std::optional<expression::Expr> singularValueDecomposition(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

// 一般実/複素行列のprecision-aware numerical SVD。A^H Aは形成しない。
// Householder bidiagonalizationの後，一側Jacobiで直交化し，reconstruction/orthogonalityを
// intervalで監査してから要求桁へdecimal化する。
[[nodiscard]] std::optional<expression::Expr> approximateSingularValueDecomposition(
    const expression::ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

} // namespace mmcal::linear_algebra
