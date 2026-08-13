#pragma once

#include "approximation/approximation_context.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "linear_algebra/exact_matrix.hpp"
#include "linear_algebra/matrix.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <cstddef>
#include <optional>

namespace mmcal::linear_algebra {

// PA=LU。返値はshape {3,n,n}で、順にP,L,U。
[[nodiscard]] std::optional<expression::Expr> luDecomposition(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

// reduced A=QR。k=min(m,n)としてQ:m×k, R:k×nを一般brace {Q,R}で返す。
// factor shapeが同じ場合だけbraceValueがdense Arrayへ自動最適化する。
[[nodiscard]] std::optional<expression::Expr> qrDecomposition(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

[[nodiscard]] std::optional<expression::Expr> approximateLuDecomposition(
    const expression::ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] std::optional<expression::Expr> approximateQrDecomposition(
    const expression::ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

// benchmark用。blockColumns=1がunblocked相当。
[[nodiscard]] std::optional<expression::Expr> approximateQrDecompositionWithBlockSize(
    const expression::ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context,
    std::size_t blockColumns);

} // namespace mmcal::linear_algebra
