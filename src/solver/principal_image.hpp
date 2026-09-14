#pragma once

#include "expression/expr.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"

#include <optional>
#include <vector>

namespace mmcal::evaluation {
class BuiltinRegistry;
}
namespace mmcal::mathematics {
class AngleSemantics;
}

namespace mmcal::solver {

// principal function の入力domainを制限したときの像を指定する。
// Complexは通常のprincipal image全体，Realは実入力だけから到達できる部分像。
enum class PrincipalImageInputDomain {
    Complex,
    Real
};

// principal branch の像を DNF（OR of AND conditions）で保持する。
// alternatives.empty() は像外，{{}} は無条件に像内を表す。
struct PrincipalImageAnalysis final {
    std::vector<mathematics::AssumptionSet> alternatives;

    [[nodiscard]] bool rejected() const noexcept { return alternatives.empty(); }
    [[nodiscard]] bool unconditional() const noexcept {
        return alternatives.size() == 1 && alternatives.front().empty();
    }
};

// FunctionBranchRuleで定義されたprincipal branchの値域を共通判定する。
// 現在はsqrt / Log系 / asin / acos / atan / asinh / acosh / atanhを扱う。
// Real指定ではprincipal函数を実入力へ制限した像を返し，Real solveでもbranch cut越しの
// 複素主値（例: log[-1]=I Pi, asin[x>1]=Pi/2-I a）を安全に反転できる。
// 記号値では像を有限個の条件領域へ分解し，既知のassumptionとcertified比較で
// 不可能な領域を除去する。未対応branch ruleはnulloptを返す。
[[nodiscard]] std::optional<PrincipalImageAnalysis> analyzePrincipalImage(
    mathematics::FunctionBranchRule branchRule,
    const expression::Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {},
    PrincipalImageInputDomain inputDomain = PrincipalImageInputDomain::Complex);

} // namespace mmcal::solver
