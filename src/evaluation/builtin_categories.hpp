#pragma once

#include "evaluation/builtin_registry.hpp"

namespace mmcal::evaluation {

// BuiltinIdのfamily判定を各subsystemのswitchへ複製しないための共通predicate。
// 新しい同系函数を追加した場合は，registryとこの分類だけを更新すればよい。
[[nodiscard]] constexpr bool isStatisticBuiltin(BuiltinId id) noexcept {
    switch (id) {
    case BuiltinId::Median:
    case BuiltinId::Mode:
    case BuiltinId::Quantile:
    case BuiltinId::Percentile:
    case BuiltinId::VariancePopulation:
    case BuiltinId::VarianceSample:
    case BuiltinId::StddevPopulation:
    case BuiltinId::StddevSample:
    case BuiltinId::GeometricMean:
    case BuiltinId::HarmonicMean:
    case BuiltinId::Rms:
    case BuiltinId::MedianAbsoluteDeviation:
    case BuiltinId::MeanAbsoluteDeviation:
    case BuiltinId::Skewness:
    case BuiltinId::KurtosisPopulation:
    case BuiltinId::KurtosisSample:
    case BuiltinId::CoefficientVariation:
    case BuiltinId::StandardError:
    case BuiltinId::ZScore:
    case BuiltinId::Iqr:
    case BuiltinId::TrimMean:
    case BuiltinId::WinsorMean:
    case BuiltinId::Winsorized:
    case BuiltinId::Covariance:
    case BuiltinId::Correlation:
    case BuiltinId::SpearmanCorrelation:
    case BuiltinId::PercentRank:
        return true;
    default:
        return false;
    }
}

[[nodiscard]] constexpr bool isVectorBuiltin(BuiltinId id) noexcept {
    switch (id) {
    case BuiltinId::VectorAdd:
    case BuiltinId::VectorSubtract:
    case BuiltinId::VectorScale:
    case BuiltinId::VectorCross:
    case BuiltinId::VectorNorm:
    case BuiltinId::VectorManhattan:
    case BuiltinId::VectorEuclidean:
    case BuiltinId::VectorNormalize:
    case BuiltinId::VectorProject:
    case BuiltinId::VectorAngle:
    case BuiltinId::VectorReflect:
    case BuiltinId::VectorReflectAxis:
    case BuiltinId::VectorSum:
    case BuiltinId::VectorInner:
    case BuiltinId::VectorOuter:
    case BuiltinId::VectorRejection:
    case BuiltinId::OrthogonalQ:
    case BuiltinId::OrthonormalQ:
    case BuiltinId::LinearIndependentQ:
    case BuiltinId::GramSchmidt:
        return true;
    default:
        return false;
    }
}

[[nodiscard]] constexpr bool isVectorCalculusBuiltin(BuiltinId id) noexcept {
    switch (id) {
    case BuiltinId::Gradient:
    case BuiltinId::Divergence:
    case BuiltinId::Curl:
    case BuiltinId::Laplacian:
    case BuiltinId::Jacobian:
    case BuiltinId::Hessian:
    case BuiltinId::DirectionalDerivative:
        return true;
    default:
        return false;
    }
}

[[nodiscard]] constexpr bool isMatrixBuiltin(BuiltinId id) noexcept {
    switch (id) {
    case BuiltinId::Transpose:
    case BuiltinId::ConjugateTranspose:
    case BuiltinId::MatrixAdd:
    case BuiltinId::MatrixMultiply:
    case BuiltinId::Determinant:
    case BuiltinId::Inverse:
    case BuiltinId::Rref:
    case BuiltinId::Rank:
    case BuiltinId::SolveLinear:
    case BuiltinId::NullSpace:
    case BuiltinId::LuDecomposition:
    case BuiltinId::QrDecomposition:
    case BuiltinId::SingularValueDecomposition:
    case BuiltinId::ConditionNumber:
    case BuiltinId::LeastSquares:
    case BuiltinId::PseudoInverse:
    case BuiltinId::Eigenvalues:
    case BuiltinId::Eigenvectors:
    case BuiltinId::Eigensystem:
    case BuiltinId::Trace:
    case BuiltinId::Rows:
    case BuiltinId::Cols:
    case BuiltinId::Diag:
        return true;
    default:
        return false;
    }
}

[[nodiscard]] constexpr bool requiresRectangularArray(BuiltinId id) noexcept {
    return isMatrixBuiltin(id) || isVectorBuiltin(id) || isVectorCalculusBuiltin(id);
}

} // namespace mmcal::evaluation
