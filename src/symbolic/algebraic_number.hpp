#pragma once

#include "numeric/rational.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace mmcal::symbolic {

struct RationalRootInterval final {
    numeric::Rational lower;
    numeric::Rational upper;

    [[nodiscard]] bool isPoint() const noexcept { return lower == upper; }
};

struct RationalComplexDisk final {
    numeric::Rational real;
    numeric::Rational imaginary;
    numeric::Rational radius;
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

// square-free有理係数多項式の1つの複素根をexactに識別するIR。
// rootIndexは Re(z)+Pi Im(z) の昇順で全複素根を並べた1-based index。
// Piが超越数で各rootの実部・虚部が代数数なので，異なる代数根のkeyは一致しない。
// diskはRational中心・Rational半径を持ち，Rouche判定でその内部に根がちょうど1個あることを証明する。
class ComplexAlgebraicNumber final {
public:
    [[nodiscard]] static std::optional<ComplexAlgebraicNumber> create(
        std::span<const numeric::Rational> polynomial,
        std::size_t rootIndex);

    [[nodiscard]] static std::optional<std::vector<ComplexAlgebraicNumber>> isolateAll(
        std::span<const numeric::Rational> polynomial);

    [[nodiscard]] std::span<const numeric::Rational> polynomial() const noexcept;
    [[nodiscard]] std::size_t degree() const noexcept;
    [[nodiscard]] std::size_t rootIndex() const noexcept;
    [[nodiscard]] const RationalComplexDisk& isolatingDisk() const noexcept;

    [[nodiscard]] RationalComplexDisk refined(std::size_t precisionBits) const;

private:
    std::vector<numeric::Rational> polynomial_;
    std::size_t rootIndex_ = 0;
    RationalComplexDisk disk_;

    ComplexAlgebraicNumber(
        std::vector<numeric::Rational> polynomial,
        std::size_t rootIndex,
        RationalComplexDisk disk);
};

enum class AlgebraicRootDomain {
    Real,
    Complex
};

enum class AlgebraicBinaryOperation {
    Add,
    Subtract,
    Multiply,
    Divide
};

// Real/Complex Rootの共通exact代数数view。定義多項式とisolating regionでroot identityを保持し，
// bounded resultant arithmeticでは演算後の候補多項式から正しいrootを再分離する。
// primitive-element最小多項式化までは要求せず，結果多項式はsquare-free canonical formで保持する。
class AlgebraicNumber final {
public:
    [[nodiscard]] static std::optional<AlgebraicNumber> create(
        std::span<const numeric::Rational> polynomial,
        std::size_t rootIndex,
        AlgebraicRootDomain domain);

    [[nodiscard]] static std::optional<AlgebraicNumber> fromRational(
        const numeric::Rational& value);
    [[nodiscard]] static std::optional<AlgebraicNumber> fromComplexRational(
        const numeric::Rational& real,
        const numeric::Rational& imaginary);
    [[nodiscard]] static std::optional<AlgebraicNumber> combine(
        const AlgebraicNumber& lhs,
        const AlgebraicNumber& rhs,
        AlgebraicBinaryOperation operation);

    [[nodiscard]] AlgebraicRootDomain domain() const noexcept;
    [[nodiscard]] std::span<const numeric::Rational> polynomial() const noexcept;
    [[nodiscard]] std::size_t rootIndex() const noexcept;
    [[nodiscard]] const RealAlgebraicNumber* asReal() const noexcept;
    [[nodiscard]] const ComplexAlgebraicNumber* asComplex() const noexcept;
    // root自体がRationalまたはRational+i Rationalへ退化する場合だけexact成分を返す。
    [[nodiscard]] std::optional<std::pair<numeric::Rational, numeric::Rational>>
        exactRationalParts() const;

private:
    std::variant<RealAlgebraicNumber, ComplexAlgebraicNumber> value_;

    explicit AlgebraicNumber(RealAlgebraicNumber value);
    explicit AlgebraicNumber(ComplexAlgebraicNumber value);
};

} // namespace mmcal::symbolic
