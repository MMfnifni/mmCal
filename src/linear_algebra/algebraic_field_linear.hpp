#pragma once

#include "symbolic/number_field.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace mmcal::linear_algebra {

// 単一のNumberFieldContext上で行うexact Gaussian elimination。
// 行列要素はすべて同じembedded fieldに属することを前提とし，pivot判定は
// AlgebraicElement::isZero()でexactに行う。
struct AlgebraicFieldRref final {
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::vector<symbolic::AlgebraicElement> elements;
    std::vector<std::size_t> pivotColumns;

    [[nodiscard]] const symbolic::AlgebraicElement& operator()(
        std::size_t row, std::size_t column) const noexcept {
        return elements[row * columns + column];
    }
};

[[nodiscard]] std::optional<AlgebraicFieldRref> algebraicFieldRref(
    std::size_t rows,
    std::size_t columns,
    std::span<const symbolic::AlgebraicElement> elements);

// kernel basisをfree-column昇順で返す。各vectorはlength=columns。
[[nodiscard]] std::optional<std::vector<std::vector<symbolic::AlgebraicElement>>>
algebraicFieldNullSpace(
    std::size_t rows,
    std::size_t columns,
    std::span<const symbolic::AlgebraicElement> elements);

} // namespace mmcal::linear_algebra
