#pragma once

#include "decimal_approximation.hpp"

#include <string>
#include <string_view>

namespace mmcal::numeric {

// certifiedな複素近似をユーザーへ提示するための出力表現。
// DecimalApproximationと同様、反復算法のworking valueではなく「確定済み表示結果 + certified enclosure」を保持し、後続のcertified四則演算へ再利用できる。
class ComplexDecimalApproximation final {
public:
    [[nodiscard]] static ComplexDecimalApproximation fromComponents(
        DecimalApproximation real,
        DecimalApproximation imaginary,
        bool realExactlyZero,
        bool imaginaryExactlyZero);

    [[nodiscard]] std::string_view text() const noexcept;
    [[nodiscard]] const DecimalApproximation& real() const noexcept;
    [[nodiscard]] const DecimalApproximation& imaginary() const noexcept;
    // legacy accessor。意味は「後続計算でexact zeroとして再利用してよいcomponent」。
    [[nodiscard]] bool realExactlyZero() const noexcept;
    [[nodiscard]] bool imaginaryExactlyZero() const noexcept;
    [[nodiscard]] bool realCertifiedExactlyZero() const noexcept;
    [[nodiscard]] bool imaginaryCertifiedExactlyZero() const noexcept;
    [[nodiscard]] bool realInformationExactlyZero() const noexcept;
    [[nodiscard]] bool imaginaryInformationExactlyZero() const noexcept;
    // Complex値として再利用可能なinformation bounds。information-exact-zero証明済み成分は
    // 表示用DecimalApproximationのzero-centered量子幅ではなくpoint zeroを返す。
    [[nodiscard]] const Rational& realInformationLower() const noexcept;
    [[nodiscard]] const Rational& realInformationUpper() const noexcept;
    [[nodiscard]] const Rational& imaginaryInformationLower() const noexcept;
    [[nodiscard]] const Rational& imaginaryInformationUpper() const noexcept;

    // exact-zero component flagを含むprovenanceを保ったまま符号反転する。
    [[nodiscard]] ComplexDecimalApproximation negated() const;

    [[nodiscard]] bool operator==(const ComplexDecimalApproximation&) const = default;

private:
    std::string text_;
    DecimalApproximation real_;
    DecimalApproximation imaginary_;
    bool realExactlyZero_ = false;
    bool imaginaryExactlyZero_ = false;

    ComplexDecimalApproximation(
        std::string text,
        DecimalApproximation real,
        DecimalApproximation imaginary,
        bool realExactlyZero,
        bool imaginaryExactlyZero);
};

} // namespace mmcal::numeric
