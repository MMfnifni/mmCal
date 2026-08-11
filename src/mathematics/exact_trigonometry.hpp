#pragma once

#include "angle.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "math_registry.hpp"
#include "numeric/rational.hpp"

#include <optional>

namespace mmcal::mathematics {

struct ExactAngle final {
    // 1 turn = 360 Deg = 2 Pi Rad = 400 Grad。
    numeric::Rational turns;
};

// 式が「Piの有理数倍」である場合だけ係数qを返す。
// 例: Pi -> 1, 2 Pi -> 2, Pi/6 -> 1/6。
[[nodiscard]] std::optional<numeric::Rational> extractRationalPiMultiple(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics);

// 三角函数へ渡された角度式を、厳密に可能な場合だけturnへ変換する。
// 未指定角度はAngleSemantics::defaultUnit()に従う。
[[nodiscard]] std::optional<ExactAngle> extractExactAngle(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AngleSemantics& angleSemantics);

// sin/cos/tan/cot/sec/cscの数学知識によるexact簡約。
// 簡約不能ならstd::nulloptを返し、呼出側は元の記号式を保持する。
[[nodiscard]] std::optional<expression::Expr> simplifyExactTrig(
    FunctionId function,
    const expression::Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AngleSemantics& angleSemantics);


// principal inverse trigをexactに決定できる特殊値だけ簡約する。
// 戻り値はAngleSemanticsの現在単位で表す。Radian既定なら asin[1/2] -> Pi/6。
[[nodiscard]] std::optional<expression::Expr> simplifyExactInverseTrig(
    FunctionId function,
    const expression::Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AngleSemantics& angleSemantics);

// atan2(y,x) のprincipal値 (-1/2,1/2] turn をexactに決定できる場合だけ返す。
[[nodiscard]] std::optional<expression::Expr> simplifyExactAtan2(
    const expression::Expr& y,
    const expression::Expr& x,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AngleSemantics& angleSemantics);

// 既存のRational Taylor backendは「ラジアン値」を受け取る。
// Nから利用するため、明示Radまたは既定Radianのexact realだけを取り出す。
[[nodiscard]] std::optional<numeric::RealNumber> extractExactRadianValue(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const AngleSemantics& angleSemantics);

} // namespace mmcal::mathematics
