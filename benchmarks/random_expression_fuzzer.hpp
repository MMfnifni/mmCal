#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>

namespace mmcal::benchmarks {

struct RandomExpressionFuzzerOptions final {
    std::uint64_t seed = 0;
    std::uint64_t cases = 10000;
    std::optional<std::uint64_t> singleCase;
    std::size_t maxDepth = 16;
    std::uint64_t reportEvery = 10000;
    std::size_t threads = 1;
    bool loop = false;
    bool noStopLoop = false;
};

// 文法と型をある程度理解した式を生成し，数学的不変量を検査する。
// falseはFAILを1件検出して詳細を表示済みであることを表す。
[[nodiscard]] bool runRandomExpressionFuzzer(const RandomExpressionFuzzerOptions& options);

} // namespace mmcal::benchmarks
