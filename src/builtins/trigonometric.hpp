#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <span>

namespace mmcal::builtins {

// 三角函数は裸の角度をAngleSemanticsの既定単位として解釈する。
// 逆三角函数も同じ設定の単位で値を返すが、内部branch定義はラジアンで固定する。
[[nodiscard]] expression::Expr evaluateSin(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);
[[nodiscard]] expression::Expr evaluateCos(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);
[[nodiscard]] expression::Expr evaluateTan(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);
[[nodiscard]] expression::Expr evaluateCot(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);
[[nodiscard]] expression::Expr evaluateSec(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);
[[nodiscard]] expression::Expr evaluateCsc(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);
[[nodiscard]] expression::Expr evaluateAsin(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);
[[nodiscard]] expression::Expr evaluateAcos(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);
[[nodiscard]] expression::Expr evaluateAtan(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);
[[nodiscard]] expression::Expr evaluateAtan2(
    std::span<const expression::Expr>, const evaluation::BuiltinRegistry&,
    const mathematics::MathRegistry&, const mathematics::AngleSemantics&);

} // namespace mmcal::builtins
