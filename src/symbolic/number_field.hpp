#pragma once

#include "algebraic_number.hpp"

#include <cstddef>
#include <memory>
#include <mutex>
#include <optional>
#include <span>
#include <vector>

namespace mmcal::symbolic {

enum class AlgebraicSign {
    Negative,
    Zero,
    Positive
};

// 選択されたembeddingまで固定したQ上のsimple extension。
// generator()が表すthetaを用い，power basis 1,theta,...,theta^(d-1) を暗黙に採用する。
class NumberFieldContext final {
public:
    // generatorはQ上既約と証明済みのminimal polynomialを持つことを前提とする。
    // 同じembedded generator identityはbounded weak internerでContextを共有する。
    [[nodiscard]] static std::shared_ptr<const NumberFieldContext> create(
        AlgebraicNumber generator);

    [[nodiscard]] const AlgebraicNumber& generator() const noexcept;
    [[nodiscard]] AlgebraicRootDomain domain() const noexcept;
    [[nodiscard]] std::span<const numeric::Rational> minimalPolynomial() const noexcept;
    [[nodiscard]] std::span<const numeric::Rational> reduction() const noexcept;
    [[nodiscard]] std::size_t degree() const noexcept;

    [[nodiscard]] std::vector<numeric::Rational> multiply(
        std::span<const numeric::Rational> lhs,
        std::span<const numeric::Rational> rhs) const;
    [[nodiscard]] std::optional<std::vector<numeric::Rational>> reciprocal(
        std::span<const numeric::Rational> value) const;
    [[nodiscard]] std::optional<std::vector<numeric::Rational>> minimalPolynomialOf(
        std::span<const numeric::Rational> value) const;

private:
    AlgebraicNumber generator_;
    std::vector<numeric::Rational> minimalPolynomial_;
    std::vector<numeric::Rational> reduction_;

    struct ReciprocalCacheEntry final {
        std::vector<numeric::Rational> value;
        std::vector<numeric::Rational> reciprocal;
    };
    mutable std::mutex reciprocalCacheMutex_;
    mutable std::vector<ReciprocalCacheEntry> reciprocalCache_;

    struct MinimalPolynomialCacheEntry final {
        std::vector<numeric::Rational> value;
        std::vector<numeric::Rational> polynomial;
    };
    mutable std::mutex minimalPolynomialCacheMutex_;
    mutable std::vector<MinimalPolynomialCacheEntry> minimalPolynomialCache_;

    NumberFieldContext(
        AlgebraicNumber generator,
        std::vector<numeric::Rational> minimalPolynomial,
        std::vector<numeric::Rational> reduction);

    [[nodiscard]] std::optional<std::vector<numeric::Rational>> findCachedReciprocal(
        std::span<const numeric::Rational> value) const;
    [[nodiscard]] std::optional<std::vector<numeric::Rational>> reciprocalUncached(
        std::span<const numeric::Rational> value) const;
    void publishReciprocal(
        std::vector<numeric::Rational> value,
        std::vector<numeric::Rational> reciprocal) const;

    [[nodiscard]] std::optional<std::vector<numeric::Rational>> findCachedMinimalPolynomial(
        std::span<const numeric::Rational> value) const;
    [[nodiscard]] std::optional<std::vector<numeric::Rational>> minimalPolynomialUncached(
        std::span<const numeric::Rational> value) const;
    void publishMinimalPolynomial(
        std::vector<numeric::Rational> value,
        std::vector<numeric::Rational> polynomial) const;
};

// NumberFieldContextのpower basis上のexact座標。
class AlgebraicElement final {
public:
    [[nodiscard]] static std::optional<AlgebraicElement> create(
        std::shared_ptr<const NumberFieldContext> field,
        std::vector<numeric::Rational> coefficients);
    [[nodiscard]] static std::optional<AlgebraicElement> generator(
        std::shared_ptr<const NumberFieldContext> field);

    [[nodiscard]] const std::shared_ptr<const NumberFieldContext>& field() const noexcept;
    [[nodiscard]] std::span<const numeric::Rational> coefficients() const noexcept;
    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] std::optional<bool> exactEquals(const AlgebraicElement& rhs) const noexcept;
    // chosen real embedding上の値をexact Rational intervalで囲む。Complex fieldではnullopt。
    [[nodiscard]] std::optional<RationalRootInterval> refinedRealInterval(
        std::size_t precisionBits) const;
    // 実embedding上でのexact符号。現在のrefinement budget内で分離できない場合はnullopt。
    [[nodiscard]] std::optional<AlgebraicSign> exactSign() const;

    [[nodiscard]] std::optional<AlgebraicElement> add(const AlgebraicElement& rhs) const;
    [[nodiscard]] std::optional<AlgebraicElement> subtract(const AlgebraicElement& rhs) const;
    [[nodiscard]] std::optional<AlgebraicElement> multiply(const AlgebraicElement& rhs) const;
    [[nodiscard]] std::optional<AlgebraicElement> divide(const AlgebraicElement& rhs) const;

    [[nodiscard]] std::optional<std::vector<numeric::Rational>> minimalPolynomial() const;

private:
    std::shared_ptr<const NumberFieldContext> field_;
    std::vector<numeric::Rational> coefficients_;

    AlgebraicElement(
        std::shared_ptr<const NumberFieldContext> field,
        std::vector<numeric::Rational> coefficients);

    [[nodiscard]] bool hasSameField(const AlgebraicElement& rhs) const noexcept;
};

} // namespace mmcal::symbolic
