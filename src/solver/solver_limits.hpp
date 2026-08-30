#pragma once

#include <cstddef>

namespace mmcal::solver::limits {

// solver内部の計算量制限。数学的な可解性の境界ではない。
inline constexpr std::size_t realProofNodes = 192;
inline constexpr std::size_t realAnalysisNodes = 256;
inline constexpr std::size_t realDomainPieces = 24;
inline constexpr std::size_t realCriticalPoints = 16;

} // namespace mmcal::solver::limits
