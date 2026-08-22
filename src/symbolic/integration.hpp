#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"

#include <optional>

namespace mmcal::symbolic {


// 未評価積分の理由を、数学的な「閉形式がない」と現在の実装能力不足で混同しないための分類。
// KnownNoFiniteClosedForm は一般的不可能性の証明ではなく、mmCalが明示的に認識している
// familyについて「現在採用している有限個の標準函数では閉形式を持たない」と判断した場合だけ使う。
enum class IntegrationDisposition {
    Solved,
    Partial,
    UnsupportedByEngine,
    KnownNoFiniteClosedForm,
    ConditionsRequired
};

struct IntegrationResult final {
    expression::Expr expression;
    IntegrationDisposition disposition = IntegrationDisposition::Solved;
};

[[nodiscard]] IntegrationResult integrateExpressionDetailed(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// 不定積分の原始函数代表元を返す。積分定数は表示しない。
// assumptionsは積分前の安全な簡約・definedness証明に利用する。
// 未対応の場合は integrate[expression, variable] をそのまま返す。
// 1変数の有理函数として認識できる式を通分・多項式GCDでexact正規化する。
// Dの後処理等でも共有し，数値近似による同値判定は行わない。
[[nodiscard]] std::optional<expression::Expr> normalizeRationalExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] expression::Expr integrateExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// 実区間の記号定積分。Infinity endpointと証明可能なimproper integralも扱う。
// 安全性を証明できない場合は integrate[...] を保持する。
[[nodiscard]] expression::Expr integrateExpression(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const expression::Expr& lower,
    const expression::Expr& upper,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::symbolic
