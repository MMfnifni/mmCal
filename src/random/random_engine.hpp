#pragma once

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <array>
#include <cstdint>

namespace mmcal::random {

// KernelSession内で共有する再現可能なPRNG。
// xoshiro256**本体とseed展開を自前で固定し、標準libraryのdistribution差に依存しない。
class RandomEngine final {
public:
    RandomEngine();

    void seed(const numeric::BigInt& seedValue);
    [[nodiscard]] numeric::BigInt reseedFromEntropy();

    [[nodiscard]] std::uint64_t nextUint64() noexcept;
    [[nodiscard]] numeric::BigInt uniformBelow(const numeric::BigInt& exclusiveUpper);
    [[nodiscard]] numeric::Rational unit53();
    [[nodiscard]] numeric::Rational positiveUnit53();

private:
    std::array<std::uint64_t, 4> state_{};
};

} // namespace mmcal::random
