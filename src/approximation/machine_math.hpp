#pragma once

#include "elementary_function.hpp"

namespace mmcal::approximation {

// 将来のApproximateRealから利用する、標準double数学函数の薄い境界層。
// 数学的な定義域判定や実数/複素数への昇格は上位層の責務とする。
class MachineMath final {
public:
    [[nodiscard]] static double evaluate(ElementaryFunction function, double value);
    [[nodiscard]] static double atan2(double y, double x);
    [[nodiscard]] static double pow(double base, double exponent);
};

} // namespace mmcal::approximation
