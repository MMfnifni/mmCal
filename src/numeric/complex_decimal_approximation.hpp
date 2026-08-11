#pragma once

#include "decimal_approximation.hpp"

#include <string>
#include <string_view>

namespace mmcal::numeric {

// certifiedな複素近似をユーザーへ提示するための出力表現。
// DecimalApproximationと同様、内部演算値ではなく「確定済みの表示結果」を保持する。
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
    [[nodiscard]] bool realExactlyZero() const noexcept;
    [[nodiscard]] bool imaginaryExactlyZero() const noexcept;
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
