#pragma once

namespace mmcal::mathematics {

// 記号計算では「真でも偽でもあるとまだ証明できない」が普通に起こる。
// boolへ潰すと未知をfalseと誤認するため、知識層では三値を使う。
enum class TruthValue {
    False,
    Unknown,
    True
};

[[nodiscard]] constexpr TruthValue logicalNot(TruthValue value) noexcept {
    if (value == TruthValue::True)
        return TruthValue::False;
    if (value == TruthValue::False)
        return TruthValue::True;
    return TruthValue::Unknown;
}

} // namespace mmcal::mathematics
