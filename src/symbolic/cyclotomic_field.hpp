#pragma once

#include "numeric/rational.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace mmcal::symbolic {

// exact Fourier backend向けのcyclotomic quotient field Q[t]/Phi_n(t)。
// embedding/root isolationは持たず，tを選択したprimitive n-th root of unityとして
// FFT側が解釈する。計算本体はpower-basis係数だけで完結する。
class CyclotomicFieldContext final {
public:
    [[nodiscard]] static std::shared_ptr<const CyclotomicFieldContext> create(
        std::size_t conductor);

    [[nodiscard]] std::size_t conductor() const noexcept;
    [[nodiscard]] std::size_t degree() const noexcept;
    [[nodiscard]] std::span<const numeric::Rational> minimalPolynomial() const noexcept;
    [[nodiscard]] std::span<const numeric::Rational> power(std::size_t exponent) const;

    [[nodiscard]] std::vector<numeric::Rational> multiply(
        std::span<const numeric::Rational> lhs,
        std::span<const numeric::Rational> rhs) const;

    // a+b I をQ(zeta_n)へexactに埋め込む。imaginary!=0なら4|nを要求する。
    [[nodiscard]] std::optional<std::vector<numeric::Rational>> embedGaussianRational(
        const numeric::Rational& real,
        const numeric::Rational& imaginary) const;
    [[nodiscard]] std::optional<std::pair<numeric::Rational, numeric::Rational>>
        exactGaussianRational(std::span<const numeric::Rational> value) const;

private:
    std::size_t conductor_ = 0;
    std::vector<numeric::Rational> polynomial_;
    std::vector<numeric::Rational> reduction_;
    std::vector<std::vector<numeric::Rational>> powers_;

    CyclotomicFieldContext(
        std::size_t conductor,
        std::vector<numeric::Rational> polynomial,
        std::vector<numeric::Rational> reduction,
        std::vector<std::vector<numeric::Rational>> powers);
};

[[nodiscard]] std::size_t cyclotomicDegree(std::size_t conductor) noexcept;

} // namespace mmcal::symbolic
