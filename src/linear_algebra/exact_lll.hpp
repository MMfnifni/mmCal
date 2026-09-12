#pragma once

#include "numeric/big_int.hpp"

#include <cstddef>
#include <vector>

namespace mmcal::linear_algebra {

using IntegerLatticeRow = std::vector<numeric::BigInt>;
using IntegerLatticeBasis = std::vector<IntegerLatticeRow>;

struct ExactLllOptions final {
    // このkernelは小規模なexact lattice専用である。dimensionだけでなく，
    // Rational Gram-Schmidt中の係数膨張と反復回数にも独立した防壁を置く。
    std::size_t maximumRank = 24;
    std::size_t maximumColumns = 192;
    std::size_t maximumIntermediateBits = 8192;
    std::size_t maximumSwaps = 4096;
    std::size_t maximumSizeReductions = 32768;
};

enum class ExactLllStatus {
    Reduced,
    DimensionMismatch,
    RankLimitExceeded,
    ColumnLimitExceeded,
    RankDeficient,
    IntermediateBitLimitExceeded,
    SwapLimitExceeded,
    SizeReductionLimitExceeded
};

struct ExactLllResult final {
    IntegerLatticeBasis basis;
    // basis == transformation * input。budget停止時にもこの不変条件を保つ。
    IntegerLatticeBasis transformation;
    ExactLllStatus status = ExactLllStatus::DimensionMismatch;
    std::size_t swaps = 0;
    std::size_t sizeReductions = 0;
    std::size_t maximumObservedBits = 0;

    [[nodiscard]] bool reduced() const noexcept {
        return status == ExactLllStatus::Reduced;
    }
};

// delta=3/4，nearest-even size reductionのexact row-LLL。
// 入力rowは一次独立でなければならない。近似浮動小数点は用いない。
[[nodiscard]] ExactLllResult exactLllReduceRows(
    IntegerLatticeBasis basis,
    ExactLllOptions options = {});

} // namespace mmcal::linear_algebra
