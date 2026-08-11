#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "solution_set.hpp"

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

// 互換用の等式solver entry point。ordered relationも渡せるが、新規コードは上を使う。
[[nodiscard]] SolutionSet solvePolynomialEquation(
    const expression::Expr& equation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// 多変数の一次方程式系。係数はexact Rationalに限定し、
// 一意解/矛盾を完全に判定できない場合はUnresolvedを返す。
[[nodiscard]] SolutionSet solveLinearPolynomialSystem(
    std::span<const expression::Expr> equations,
    std::span<const expression::Symbol> variables,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::solver
