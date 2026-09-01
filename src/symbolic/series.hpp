#pragma once

#include "expression/expr.hpp"
#include "evaluation/builtin_registry.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"
#include "mathematics/angle.hpp"

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace mmcal::symbolic {

// SeriesDataは公開ExprKindを増やさず，内部builtin head seriesData[...]で保持する。
// exponentDenominatorはPuiseux級数の指数格子1/qを表す。
// minimumExponent/orderNumeratorはその格子上の整数numeratorとして保持する。
// logarithmicCoefficients[k]はlog(x-center)^(k+1)に掛かる同一指数格子の係数列である。
// center==Infinityでは局所変数を1/xとし，指数格子・log層も(1/x)^(n/q)，log[1/x]^kとして解釈する。
struct SeriesData final {
    expression::Symbol variable;
    expression::Expr center;
    std::vector<expression::Expr> coefficients;
    std::int64_t minimumExponent = 0;
    std::int64_t orderNumerator = 1;
    std::uint32_t exponentDenominator = 1;
    std::vector<std::vector<expression::Expr>> logarithmicCoefficients;
};

[[nodiscard]] std::optional<SeriesData> parseSeriesData(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins);

[[nodiscard]] expression::Expr makeSeriesData(
    SeriesData data,
    const evaluation::BuiltinRegistry& builtins);

// v1.5.5 WIP: truncated Taylor/Laurent/Puiseux series arithmetic。
// 四則・整数冪に加え，exp/log/sin/cos/sinh/cosh，principal sqrt / exact有理冪，
// erf/Si/Ei/Ciの局所展開をTPSA係数漸化式・primitive compositionで扱う。
// logはprincipal branch上で展開中心の正則性を証明できる場合だけ展開し，
// 非整数冪もprincipal branch上の正則性を証明できる場合だけ展開する。
// direct trigはsession angle modeと明示Rad/Deg/Gradを同じ意味論で扱う。
[[nodiscard]] std::optional<expression::Expr> seriesExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const expression::Expr& center,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

[[nodiscard]] expression::Expr normalSeriesExpression(
    const SeriesData& series,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] expression::Expr toNormalExpression(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] std::optional<expression::Expr> differentiateSeriesExpression(
    const SeriesData& series,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] std::optional<expression::Expr> integrateSeriesExpression(
    const SeriesData& series,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::symbolic
