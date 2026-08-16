#pragma once

#include "numeric/rational.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace mmcal::symbolic {

struct RationalRootInterval final {
    numeric::Rational lower;
    numeric::Rational upper;

    [[nodiscard]] bool isPoint() const noexcept { return lower == upper; }
};

// 有理係数多項式の1つの実根をexactに表すspecialized IR。
// polynomial[i]はx^iの係数。rootIndexは異なる実根を小さい順に数えた1-based index。
// isolating intervalはその根だけを含むexact Rational区間である。
class RealAlgebraicNumber final {
public:
    [[nodiscard]] static std::optional<RealAlgebraicNumber> create(
        std::span<const numeric::Rational> polynomial,
        std::size_t rootIndex);

    [[nodiscard]] static std::optional<std::vector<RealAlgebraicNumber>> isolateAll(
        std::span<const numeric::Rational> polynomial);

    [[nodiscard]] std::span<const numeric::Rational> polynomial() const noexcept;
    [[nodiscard]] std::size_t degree() const noexcept;
    [[nodiscard]] std::size_t rootIndex() const noexcept;
    [[nodiscard]] const RationalRootInterval& isolatingInterval() const noexcept;

    // unique-root intervalを要求bit数に応じて相対的に細分化する。
    [[nodiscard]] RationalRootInterval refined(std::size_t precisionBits) const;

private:
    std::vector<numeric::Rational> polynomial_;
    std::size_t rootIndex_ = 0;
    RationalRootInterval interval_;

    RealAlgebraicNumber(
        std::vector<numeric::Rational> polynomial,
        std::size_t rootIndex,
        RationalRootInterval interval);
};

} // namespace mmcal::symbolic
