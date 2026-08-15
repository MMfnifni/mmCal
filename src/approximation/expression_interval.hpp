#pragma once

#include "approximation_context.hpp"
#include "certified_evaluator.hpp"
#include "complex_interval.hpp"
#include "expression/expr.hpp"
#include "real_interval.hpp"

#include <cstddef>
#include <optional>
#include <span>

namespace mmcal::approximation {

// DecimalApproximationの既存保証区間も失わず、任意Exprを複素区間へ持ち上げる。
[[nodiscard]] std::optional<ComplexInterval> encloseComplexExpression(
    const expression::Expr& expression,
    std::size_t precisionBits,
    const CertifiedEvaluator& certified);

// certified区間を要求小数桁へ一意に丸められる場合だけExpr化する。
// exact pointはN[exact,p]と同じ最小表記を維持する。
[[nodiscard]] std::optional<expression::Expr> decimalExpression(
    const RealInterval& value,
    std::size_t fractionalDigits);
[[nodiscard]] std::optional<expression::Expr> decimalExpression(
    const ComplexInterval& value,
    std::size_t fractionalDigits);

// 既存Approximationを含む入力から、最も低い要求桁をbackend精度として推定する。
[[nodiscard]] std::optional<ApproximationContext> inferredApproximationContext(
    std::span<const expression::Expr> expressions);

// DecimalApproximation / ComplexDecimalApproximationを含むscalar四則演算。
// exact Numberはpoint intervalとして混在できる。入力approximationが一つもない場合や，
// scalar数値以外を含む場合はnulloptを返し，通常のexact/symbolic経路へ委ねる。
[[nodiscard]] std::optional<expression::Expr> addApproximateScalars(
    std::span<const expression::Expr> expressions);
[[nodiscard]] std::optional<expression::Expr> subtractApproximateScalars(
    const expression::Expr& lhs,
    const expression::Expr& rhs);
[[nodiscard]] std::optional<expression::Expr> multiplyApproximateScalars(
    std::span<const expression::Expr> expressions);
[[nodiscard]] std::optional<expression::Expr> divideApproximateScalars(
    const expression::Expr& lhs,
    const expression::Expr& rhs);
[[nodiscard]] std::optional<expression::Expr> negateApproximateScalar(
    const expression::Expr& value);

[[nodiscard]] std::size_t nextGuardDigits(std::size_t current);

} // namespace mmcal::approximation
