#pragma once

#include "complex_interval.hpp"
#include "real_interval.hpp"

namespace mmcal::approximation {

// 特異集合に対する区間の関係。ExactSingularityは区間全体がその特異点に
// 一致する場合だけを表し，MayContainSingularityは有限precision等で候補を
// 含み得るが一意には決められない場合を表す。
enum class SingularityRelation {
    Clear,
    ExactSingularity,
    MayContainSingularity
};

[[nodiscard]] bool intervalIsExactZero(const RealInterval& value) noexcept;
[[nodiscard]] bool intervalIsExactZero(const ComplexInterval& value) noexcept;

// {0,-1,-2,...} のpole集合との関係を分類する。
[[nodiscard]] SingularityRelation classifyNonPositiveIntegerPole(
    const RealInterval& value);
[[nodiscard]] SingularityRelation classifyNonPositiveIntegerPole(
    const ComplexInterval& value);

// principal branchのInformationEnclosure判定。いずれも「cut上の規約値と
// 不連続な片側候補を同時に含むか」を返す。
[[nodiscard]] bool informationCrossesPrincipalNegativeRealCut(
    const ComplexInterval& value);
[[nodiscard]] bool informationIsAmbiguousAtOuterRealCuts(
    const ComplexInterval& value);
[[nodiscard]] bool informationIsAmbiguousAtAcoshCut(
    const ComplexInterval& value);
[[nodiscard]] bool informationIsAmbiguousAtOuterImaginaryCuts(
    const ComplexInterval& value);
[[nodiscard]] bool informationIsAmbiguousAtPositiveRealCut(
    const ComplexInterval& value);

[[nodiscard]] bool informationMayContainComplexZero(
    const ComplexInterval& value) noexcept;

} // namespace mmcal::approximation
