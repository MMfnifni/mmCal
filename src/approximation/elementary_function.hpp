// 初等函数識別の公開
#pragma once

namespace mmcal::approximation {

// 近似数学バックエンドへ渡す函数ID。Evaluatorはcmathへ直接依存しない。
enum class ElementaryFunction {
    Sin,
    Cos,
    Tan,
    Asin,
    Acos,
    Atan,
    Exp,
    Log,
    Log10,
    Sqrt
};

} // namespace mmcal::approximation
