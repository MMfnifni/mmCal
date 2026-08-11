#pragma once

#include <optional>
#include <string_view>

namespace mmcal::mathematics {

// 三角函数の入力角度をどの単位として解釈するかを表す。
// mmCalは数学函数との整合性を優先し、未指定時はRadianを既定とする。
enum class AngleUnit {
    Degree,
    Radian,
    Gradian
};

class AngleSemantics final {
public:
    constexpr AngleSemantics() noexcept = default;
    explicit constexpr AngleSemantics(AngleUnit defaultUnit) noexcept
        : defaultUnit_(defaultUnit) {}

    [[nodiscard]] constexpr AngleUnit defaultUnit() const noexcept { return defaultUnit_; }
    void setDefaultUnit(AngleUnit unit) noexcept { defaultUnit_ = unit; }

    // 入力互換のため大文字・小文字の短縮名を受け付ける。
    // 出力時はDeg / Rad / Gradへ正規化する。
    [[nodiscard]] static std::optional<AngleUnit> parseUnit(std::string_view name) noexcept;
    [[nodiscard]] static std::string_view canonicalName(AngleUnit unit) noexcept;

private:
    AngleUnit defaultUnit_ = AngleUnit::Radian;
};

// 単体Evaluator等が利用する不変の既定設定。KernelSessionは専用設定を所有する。
[[nodiscard]] const AngleSemantics& defaultAngleSemantics();

} // namespace mmcal::mathematics
