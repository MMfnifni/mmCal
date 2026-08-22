#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>

namespace mmcal::benchmarks {

struct CertificationBoundaryFuzzerOptions final {
    std::uint64_t seed = 0;
    std::uint64_t cases = 10000;
    std::optional<std::uint64_t> singleCase;
    std::uint64_t reportEvery = 1000;
    std::size_t threads = 1;
    std::uint64_t caseTimeoutMilliseconds = 2000;
    bool loop = false;
    bool noStopLoop = false;
};

// branch cut・pole・backend boundary近傍のclosed numeric inputを生成し，
// certified evaluationの分類契約とbounded-workを監査する。
[[nodiscard]] bool runCertificationBoundaryFuzzer(
    const CertificationBoundaryFuzzerOptions& options);

} // namespace mmcal::benchmarks
