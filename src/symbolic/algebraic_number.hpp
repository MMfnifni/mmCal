#pragma once

#include "numeric/rational.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace mmcal::symbolic {

class AlgebraicElement;

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
    // 既約minimal polynomialと，その根をただ1つ含む候補区間から直接構築する内部fast path。
    // 区間の一意根性とroot indexはSturm列で再証明する。
    [[nodiscard]] static std::optional<RealAlgebraicNumber> createFromMinimalPolynomialInterval(
        std::span<const numeric::Rational> polynomial,
        RationalRootInterval interval);

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

enum class AlgebraicOrder {
    Less,
    Equal,
    Greater
};

// Real/Complex Rootの共通exact代数数view。canonical root identityとisolating regionを保持し，
// 証明できる場合はpersistent AlgebraicElementを内部算術表現として併置する。user-visibleな
// root[minpoly,k] identityとfield座標は分離し，budget超過時は従来resultant経路へfallbackする。
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
    [[nodiscard]] static AlgebraicNumber fromRealRoot(RealAlgebraicNumber value);
    [[nodiscard]] static AlgebraicNumber fromComplexRoot(ComplexAlgebraicNumber value);
    [[nodiscard]] static std::optional<AlgebraicNumber> combine(
        const AlgebraicNumber& lhs,
        const AlgebraicNumber& rhs,
        AlgebraicBinaryOperation operation);

    [[nodiscard]] AlgebraicRootDomain domain() const noexcept;
    [[nodiscard]] std::span<const numeric::Rational> polynomial() const noexcept;
    [[nodiscard]] std::size_t rootIndex() const noexcept;
    [[nodiscard]] const RealAlgebraicNumber* asReal() const noexcept;
    [[nodiscard]] const ComplexAlgebraicNumber* asComplex() const noexcept;
    [[nodiscard]] bool hasSameRootIdentity(const AlgebraicNumber& rhs) const noexcept;
    // exact algebraic equality。budget超過や証明不能時はnullopt。
    [[nodiscard]] std::optional<bool> exactEquals(const AlgebraicNumber& rhs) const;
    // 実代数数だけに定義するexact order。Complex root列挙順とは無関係。
    [[nodiscard]] std::optional<AlgebraicOrder> exactRealCompare(const AlgebraicNumber& rhs) const;
    [[nodiscard]] const AlgebraicElement* arithmeticElement() const noexcept;
    [[nodiscard]] AlgebraicNumber withArithmeticElement(
        std::shared_ptr<const AlgebraicElement> element) const;
    [[nodiscard]] AlgebraicNumber withGeneratorField() const;
    // root自体がRationalまたはRational+i Rationalへ退化する場合だけexact成分を返す。
    [[nodiscard]] std::optional<std::pair<numeric::Rational, numeric::Rational>>
        exactRationalParts() const;

private:
    std::variant<RealAlgebraicNumber, ComplexAlgebraicNumber> value_;
    std::shared_ptr<const AlgebraicElement> arithmeticElement_;

    explicit AlgebraicNumber(RealAlgebraicNumber value);
    explicit AlgebraicNumber(ComplexAlgebraicNumber value);
};

} // namespace mmcal::symbolic
