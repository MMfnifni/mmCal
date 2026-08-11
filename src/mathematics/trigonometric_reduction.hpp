#pragma once

#include "math_ids.hpp"
#include "numeric/rational.hpp"

namespace mmcal::mathematics {

// sin/cos/tan の周期性と象限対称性だけを使って、任意のturn値を第1象限の基準角 [0, 1/4 turn] へ厳密に写す。
//
// ここでは数値近似を一切行わない。
// 入力も出力も Rational なので、exact simplifier と certified numerical evaluator の双方が同じ縮約規則を共有できる。
struct ReducedTrigAngle final {
    numeric::Rational referenceTurns;
    bool negative = false;
};

[[nodiscard]] numeric::Rational normalizeTurns(numeric::Rational turns);

[[nodiscard]] ReducedTrigAngle reduceTrigTurns(
    FunctionId function,
    numeric::Rational turns);

} // namespace mmcal::mathematics
