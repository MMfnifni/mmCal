#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "solution_set.hpp"

#include <optional>
#include <span>

namespace mmcal::solver {

// 一変数exact polynomial relation solver。
// == に加え、Real上の < <= > >= を扱う。不等式は一変数Rational係数Polynomialから
// exactな符号表を構成し、解区間をfree-variable branchとして保持する。
[[nodiscard]] SolutionSet solveUnivariatePolynomialRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// x == exact algebraic value の直接bindingだけを拾う低コスト経路。
// Real/Rational制約付きsolveで一般函数proofへ入る前に使い，最終domain判定はapplySolveConstraintsへ任せる。
[[nodiscard]] std::optional<SolutionSet> solveDirectAlgebraicBindingRelation(
    const expression::Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics);

// 互換用の等式solver entry point。ordered relationも渡せるが、新規コードは上を使う。
// 既存radical/binomial solverで閉じない有理係数等式を、実代数数Rootでexactに補完する。
// 非多項式・現degree budget外はnulloptとして既存Unresolved semanticsを維持する。
[[nodiscard]] std::optional<SolutionSet> solveRealAlgebraicPolynomialEquation(
    const expression::Expr& equation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// Q上で複数因子へ分かれるReal polynomialだけを因子単位で解く高速経路。
// 一般の既約多項式はnulloptとして既存のproof / radical経路を優先する。
[[nodiscard]] std::optional<SolutionSet> solveFactoredRealAlgebraicPolynomialEquation(
    const expression::Expr& equation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// 高次入力がrepeated factorを持つ場合だけReal Algebraic Root fallbackを優先する。
// minimal factorが二次へ落ちても，元が高次fallbackだったという表現契約を保持する。
[[nodiscard]] std::optional<SolutionSet> solveRepeatedRealAlgebraicPolynomialEquation(
    const expression::Expr& equation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] SolutionSet solvePolynomialEquation(
    const expression::Expr& equation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// 多変数の一次方程式系。係数はexact Rationalに限定し、
// 一意解/矛盾を完全に判定できない場合はUnresolvedを返す。
// exact Rational係数の多変数多項式等式系をgeneral Gröbner basisで処理する。
// 現solve bridgeはzero-dimensional shape-position basisに加え，eliminant rootごとの
// bounded specializationで一部の非shape-position系も完全列挙する。積=0をexactに
// 分岐できるpositive-dimensional系は自由変数branchとして返す。1方程式の次数1/2
// projectionはComplex/Real free parameterへliftし，Realでは半代数的parameter条件を保持する。
// 一般多様体の推測parameterizationは行わずUnresolvedを保つ。
[[nodiscard]] std::optional<SolutionSet> solvePolynomialSystem(
    std::span<const expression::Expr> equations,
    std::span<const expression::Symbol> variables,
    mathematics::NumericDomain ambientDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] SolutionSet solveLinearPolynomialSystem(
    std::span<const expression::Expr> equations,
    std::span<const expression::Symbol> variables,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::solver
