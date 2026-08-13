#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "linear_algebra/matrix.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <cstddef>
#include <optional>

namespace mmcal::linear_algebra {

struct ExactMatrixContext final {
    const evaluation::BuiltinRegistry& builtins;
    const mathematics::MathRegistry& mathematics;
    const mathematics::AngleSemantics& angles;
};

[[nodiscard]] bool allExactNumbers(const MatrixView& matrix) noexcept;

[[nodiscard]] std::optional<expression::Expr> determinant(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

[[nodiscard]] std::optional<expression::Expr> inverse(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

[[nodiscard]] std::optional<MatrixBuffer> rref(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

[[nodiscard]] std::optional<std::size_t> matrixRank(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

[[nodiscard]] std::optional<expression::Expr> solveLinear(
    const MatrixView& matrix,
    const expression::ArrayExpr& rhs,
    const ExactMatrixContext& context);

[[nodiscard]] std::optional<expression::Expr> nullSpace(
    const MatrixView& matrix,
    const ExactMatrixContext& context);

} // namespace mmcal::linear_algebra
