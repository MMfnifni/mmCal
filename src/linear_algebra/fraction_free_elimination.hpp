#pragma once

#include "numeric/big_int.hpp"

#include <cstddef>
#include <vector>

namespace mmcal::linear_algebra {

// Bareiss法で使う整数専用row-major作業領域。
// exact Rational行列は呼出側で行ごとに分母を払ってここへliftする。
class IntegerMatrixBuffer final {
public:
    IntegerMatrixBuffer(std::size_t rows, std::size_t columns);
    IntegerMatrixBuffer(
        std::size_t rows,
        std::size_t columns,
        std::vector<numeric::BigInt> elements);

    [[nodiscard]] std::size_t rows() const noexcept;
    [[nodiscard]] std::size_t columns() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;

    [[nodiscard]] numeric::BigInt& operator()(
        std::size_t row, std::size_t column) noexcept;
    [[nodiscard]] const numeric::BigInt& operator()(
        std::size_t row, std::size_t column) const noexcept;

    void swapRows(std::size_t lhs, std::size_t rhs) noexcept;

    [[nodiscard]] std::vector<numeric::BigInt>& elements() noexcept;
    [[nodiscard]] const std::vector<numeric::BigInt>& elements() const noexcept;

private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<numeric::BigInt> elements_;
};

struct BareissEchelonResult final {
    IntegerMatrixBuffer matrix;
    std::vector<std::size_t> pivotColumns;
    std::size_t rowSwaps = 0;
};

// pivotColumnLimitより左だけをpivot候補とし、右側列も同じfraction-free行操作で更新する。
// augmented matrixのinverse/linear solveでも同じkernelを再利用できる。
[[nodiscard]] BareissEchelonResult bareissEchelon(
    IntegerMatrixBuffer matrix,
    std::size_t pivotColumnLimit);

[[nodiscard]] numeric::BigInt bareissDeterminant(IntegerMatrixBuffer matrix);

} // namespace mmcal::linear_algebra
