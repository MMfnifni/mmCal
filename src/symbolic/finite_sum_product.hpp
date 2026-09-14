#pragma once

#include "evaluation/builtin_registry.hpp"
#include "evaluation/iterator_spec.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <optional>

namespace mmcal::symbolic {

// 記号有限和・積の閉形式kernel。
// 戻り値なしは「証明可能な閉形式をこのkernelでは構成できない」を意味し，
// 呼出側が有限展開または未評価保持へfallbackする。
[[nodiscard]] std::optional<expression::Expr> finiteSymbolicSum(
    const expression::Expr& body,
    const evaluation::TableIteratorSpec& iterator,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] std::optional<expression::Expr> finiteSymbolicProduct(
    const expression::Expr& body,
    const evaluation::TableIteratorSpec& iterator,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// hypergeometric termの基礎認識。Factorial / Combination / rising/falling factorial
// 等のshift quotientを構造的に正規化し，f(k+1)/f(k)がkの有理函数ならratioを返す。
struct HypergeometricTermRecognition final {
    expression::Expr ratio;
};

[[nodiscard]] std::optional<HypergeometricTermRecognition> recognizeHypergeometricTerm(
    const expression::Expr& body,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::symbolic
