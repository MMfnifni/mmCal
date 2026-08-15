#pragma once

#include "real_number.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>

namespace mmcal::numeric {

enum class ApproximationOrigin {
    ExactValue,
    CertifiedInterval
};

// 厳密実数やcertified区間を、ユーザーへ提示する10進表現として保持する値型。
// 反復算法の作業値ではなく「確定した表示結果 + certified enclosure」であり、BigFloat/RealIntervalとは責務を分ける。保存enclosureは後続のcertified四則演算へ再利用できる。
class DecimalApproximation final {
public:
    [[nodiscard]] static DecimalApproximation fromReal(
        const RealNumber& value,
        std::size_t repeatingFractionalDigits = 16);

    // N[expr,p]用。pは小数部桁数ではなく有効10進桁数を表す。
    [[nodiscard]] static DecimalApproximation fromRealSignificant(
        const RealNumber& value,
        std::size_t significantDigits);

    // 指定した小数部桁数へ必ず固定桁で最近接・偶数丸めする。
    // 証明済み区間の両端比較など、桁数を揃える必要がある内部処理向け。
    [[nodiscard]] static DecimalApproximation fromRealFixed(
        const RealNumber& value,
        std::size_t fractionalDigits);

    // 真値が [lower, upper] に含まれることが別途証明されている場合に使う。
    // 両端を同じ桁数へ丸めた結果が一致したときだけ、その10進表現を確定値として返す。
    // 要求桁数はmetadataに保持し，表示上は連続する末尾0を1個まで圧縮する。
    [[nodiscard]] static std::optional<DecimalApproximation> fromCertifiedInterval(
        const Rational& lower,
        const Rational& upper,
        std::size_t fractionalDigits);

    [[nodiscard]] static std::optional<DecimalApproximation> fromCertifiedIntervalSignificant(
        const Rational& lower,
        const Rational& upper,
        std::size_t significantDigits);

    // 真値保証区間とは別に，この近似値から後続計算で利用してよい情報量の区間を指定する。
    // information enclosureはcertified enclosureを必ず包含し，表示丸めの量子幅も内部で包含させる。
    [[nodiscard]] static std::optional<DecimalApproximation> fromCertifiedIntervalWithInformation(
        const Rational& certifiedLower,
        const Rational& certifiedUpper,
        const Rational& informationLower,
        const Rational& informationUpper,
        std::size_t fractionalDigits);

    [[nodiscard]] static std::optional<DecimalApproximation> fromCertifiedIntervalWithInformationSignificant(
        const Rational& certifiedLower,
        const Rational& certifiedUpper,
        const Rational& informationLower,
        const Rational& informationUpper,
        std::size_t significantDigits);

    [[nodiscard]] std::string_view text() const noexcept;
    [[nodiscard]] std::size_t fractionalDigits() const noexcept;
    [[nodiscard]] std::size_t requestedFractionalDigits() const noexcept;
    [[nodiscard]] std::size_t requestedSignificantDigits() const noexcept;
    [[nodiscard]] bool isRounded() const noexcept;
    [[nodiscard]] ApproximationOrigin origin() const noexcept;
    // 表示文字列が表す10進値そのものをexact Rationalで保持する。
    // accuracy/precisionは文字列を再parseせずInformationEnclosureを直接使う。
    [[nodiscard]] const Rational& displayedValue() const noexcept;
    [[nodiscard]] const Rational& certifiedLower() const noexcept;
    [[nodiscard]] const Rational& certifiedUpper() const noexcept;
    [[nodiscard]] const Rational& informationLower() const noexcept;
    [[nodiscard]] const Rational& informationUpper() const noexcept;
    [[nodiscard]] bool certifiedEnclosureIsPoint() const noexcept;
    [[nodiscard]] bool informationEnclosureIsPoint() const noexcept;
    [[nodiscard]] bool operator==(const DecimalApproximation&) const = default;

private:
    std::string text_;
    std::size_t fractionalDigits_ = 0;
    std::size_t requestedFractionalDigits_ = 0;
    std::size_t requestedSignificantDigits_ = 0;
    bool rounded_ = false;
    ApproximationOrigin origin_ = ApproximationOrigin::ExactValue;
    Rational displayedValue_;
    Rational certifiedLower_;
    Rational certifiedUpper_;
    Rational informationLower_;
    Rational informationUpper_;

    DecimalApproximation(
        std::string text,
        std::size_t fractionalDigits,
        std::size_t requestedFractionalDigits,
        std::size_t requestedSignificantDigits,
        bool rounded,
        ApproximationOrigin origin,
        Rational displayedValue,
        Rational certifiedLower,
        Rational certifiedUpper,
        Rational informationLower,
        Rational informationUpper);
};

} // namespace mmcal::numeric
