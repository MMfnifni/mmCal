#pragma once

#include "approximation_context.hpp"
#include "real_interval.hpp"
#include "mathematics/math_ids.hpp"
#include "numeric/decimal_approximation.hpp"

#include <cstddef>
#include <optional>

namespace mmcal::approximation {

struct CertifiedConstantResult final {
    RealInterval interval;
    std::size_t termsUsed = 0;
    std::size_t workingBits = 0;
};

// 指定した2進作業精度でPiを必ず含む区間を作る。
// 初期実装は検証しやすさを優先し、Machin公式 + atan交代級数を使う。
[[nodiscard]] CertifiedConstantResult enclosePi(std::size_t precisionBits);

// 現在certified providerを持つ数学定数だけを返す。未実装定数はnullopt。
[[nodiscard]] std::optional<CertifiedConstantResult> encloseConstant(
    mathematics::ConstantId id,
    std::size_t precisionBits);

// 要求小数桁へ最近接・偶数丸めした結果が区間全体で一意になるまで、作業精度を自動的に増やしてcertifyする。
[[nodiscard]] numeric::DecimalApproximation approximatePi(
    std::size_t fractionalDigits);

[[nodiscard]] std::optional<numeric::DecimalApproximation> approximateConstant(
    mathematics::ConstantId id,
    std::size_t fractionalDigits);

} // namespace mmcal::approximation
