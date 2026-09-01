// 解集合表現の回帰テスト
#include "solution_set_tests.hpp"

#include "builtins/array.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "solver/solution_set.hpp"
#include "solver/polynomial_solver.hpp"
#include "formatting/expr_formatter.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "symbols/symbol_table.hpp"
#include "symbolic/series.hpp"
#include "test_framework.hpp"

#include <array>
#include <stdexcept>

namespace mmcal::tests {
namespace {

using expression::Expr;
using mathematics::NumericDomain;
using mathematics::RelationKind;
using numeric::BigInt;
using numeric::Number;
using solver::SolutionBinding;
using solver::SolutionBranch;
using solver::SolverVariable;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

} // namespace

void runSolutionSetTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const auto builtins = evaluation::BuiltinRegistry::defaults(symbols);
    const auto x = symbols.intern("x");
    const auto a = symbols.intern("a");

    // x^3 == -8 のComplex上の全解を、Power[-8,1/3]という「函数値」とは別に
    // 3本のbranchとして表現できることを確認する。CubeRoot/Solveそのものはまだ
    // 実装しないが、Solver結果型は最初から多価性を失わない。
    const Expr sqrt3 = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sqrt),
        {integer(3)});
    const Expr imaginary{Number::complex(BigInt{0}, BigInt{1})};
    const Expr iSqrt3 = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Multiply),
        {imaginary, sqrt3});
    const Expr upperRoot = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Add),
        {integer(1), iSqrt3});
    const Expr lowerRoot = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Subtract),
        {integer(1), iSqrt3});

    const auto roots = solver::SolutionSet::finite(
        {SolverVariable{x, NumericDomain::Complex}},
        {
            SolutionBranch{{SolutionBinding{x, integer(-2)}}, {}, std::size_t{1}, {}, std::nullopt},
            SolutionBranch{{SolutionBinding{x, upperRoot}}, {}, std::size_t{1}, {}, std::nullopt},
            SolutionBranch{{SolutionBinding{x, lowerRoot}}, {}, std::size_t{1}, {}, std::nullopt}
        });
    tests.expect(
        roots.kind() == solver::SolutionSetKind::Finite
            && roots.branches().size() == 3
            && roots.variables().front().domain == NumericDomain::Complex,
        "SolutionSet: finite solver result keeps all branches and ambient domain");

    mathematics::AssumptionSet nonZeroA;
    nonZeroA.add(mathematics::relation(RelationKind::NotEqual, Expr{a}, integer(0)));
    const auto conditionalBranch = solver::SolutionSet::finite(
        {SolverVariable{x, NumericDomain::Complex}},
        {SolutionBranch{{SolutionBinding{x, Expr{a}}}, nonZeroA, std::nullopt, {}, std::nullopt}});
    tests.expect(
        !conditionalBranch.branches().front().unconditional()
            && conditionalBranch.branches().front().conditions.size() == 1,
        "SolutionSet: a branch preserves mathematical side conditions");

    const auto k = symbols.intern("k");
    const solver::SolutionBranch metadataBranch{
        {SolutionBinding{x, Expr{a}}},
        nonZeroA,
        std::size_t{2},
        {SolverVariable{k, NumericDomain::Integer}},
        NumericDomain::Real};
    const auto metadataSet = solver::SolutionSet::finite(
        {SolverVariable{x, NumericDomain::Real}}, {metadataBranch});
    const std::array<Expr, 2> branchSelection{
        Expr::solutionSet(metadataSet), integer(0)};
    const Expr selectedBranch = builtins::evaluateArrayGet(branchSelection);
    tests.expect(
        selectedBranch.isSolutionSet()
            && selectedBranch.asSolutionSet().kind() == solver::SolutionSetKind::Finite
            && selectedBranch.asSolutionSet().variables().size() == 1
            && selectedBranch.asSolutionSet().variables().front().domain == NumericDomain::Real
            && selectedBranch.asSolutionSet().branches().size() == 1
            && selectedBranch.asSolutionSet().branches().front() == metadataBranch,
        "at: selecting a SolutionSet branch preserves branch metadata and variable domains");

    const std::array<Expr, 3> bindingSelection{
        Expr::solutionSet(metadataSet), integer(0), Expr{x}};
    tests.expect(
        builtins::evaluateArrayGet(bindingSelection) == Expr{a},
        "at: selecting a SolutionSet binding returns its right-hand-side expression");

    const auto allReal = solver::SolutionSet::universal(
        {SolverVariable{x, NumericDomain::Real}});
    tests.expect(
        allReal.kind() == solver::SolutionSetKind::Universal
            && allReal.branches().empty(),
        "SolutionSet: universal and finite solution sets are distinct states");


    mathematics::AssumptionSet zeroA;
    zeroA.add(mathematics::relation(RelationKind::Equal, Expr{a}, integer(0)));
    const auto conditionalSet = solver::SolutionSet::conditional(
        {SolverVariable{x, NumericDomain::Complex}},
        {
            solver::SolutionCase{
                nonZeroA,
                solver::SolutionSetKind::Finite,
                {SolutionBranch{{SolutionBinding{x, integer(1)}}, {}, std::nullopt, {}, std::nullopt}}},
            solver::SolutionCase{zeroA, solver::SolutionSetKind::Universal, {}}
        });
    tests.expect(
        conditionalSet.kind() == solver::SolutionSetKind::Conditional
            && conditionalSet.cases().size() == 2,
        "SolutionSet: coefficient-dependent solver results preserve set-level cases");
    tests.expectEqual(
        formatting::formatExpr(Expr::solutionSet(conditionalSet)),
        std::string{"cases[{x == 1} if a != 0; All if a == 0]"},
        "SolutionSet: conditional cases have a readable exact representation");

    const auto math = mathematics::MathRegistry::defaults(symbols, builtins);
    const auto& angles = mathematics::defaultAngleSemantics();

    symbolic::SeriesData bindingSeries{
        a, integer(0), {integer(1), integer(1)}, 0, 2, 1, {}};
    const solver::SolutionBranch normalizableBranch{
        {SolutionBinding{x, symbolic::makeSeriesData(bindingSeries, builtins)}},
        nonZeroA, std::size_t{2}, {SolverVariable{k, NumericDomain::Integer}}, NumericDomain::Real};
    const auto normalizableSet = solver::SolutionSet::finite(
        {SolverVariable{x, NumericDomain::Real}}, {normalizableBranch});
    const Expr normalizedSolutions = symbolic::toNormalExpression(
        Expr::solutionSet(normalizableSet), builtins, math, angles);
    tests.expect(
        normalizedSolutions.isSolutionSet()
            && normalizedSolutions.asSolutionSet().kind() == solver::SolutionSetKind::Finite
            && normalizedSolutions.asSolutionSet().variables().size() == 1
            && normalizedSolutions.asSolutionSet().variables().front()
                == normalizableSet.variables().front()
            && normalizedSolutions.asSolutionSet().branches().size() == 1
            && normalizedSolutions.asSolutionSet().branches().front().conditions == nonZeroA
            && normalizedSolutions.asSolutionSet().branches().front().multiplicity == std::optional<std::size_t>{2}
            && normalizedSolutions.asSolutionSet().branches().front().freeVariables
                == std::vector<SolverVariable>{{k, NumericDomain::Integer}}
            && normalizedSolutions.asSolutionSet().branches().front().bindingsCertifiedDomain
                == std::optional<NumericDomain>{NumericDomain::Real},
        "toNormal: SolutionSet conversion preserves branch metadata and solver domains");
    tests.expectEqual(
        formatting::formatExpr(normalizedSolutions),
        std::string{"{x == a+1 where k in Integer if a != 0 (multiplicity 2)}"},
        "toNormal: SolutionSet conversion recursively normalizes binding values");

    const Expr quadraticEquation = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Equal),
        {
            Expr::call(builtins.symbol(evaluation::BuiltinId::Power), {Expr{x}, integer(2)}),
            integer(1)
        });
    const auto quadratic = solver::solvePolynomialEquation(
        quadraticEquation, x, builtins, math, angles);
    tests.expect(
        quadratic.kind() == solver::SolutionSetKind::Finite
            && quadratic.branches().size() == 2,
        "PolynomialSolver: quadratic equation returns both complex-domain roots");
    tests.expectEqual(
        formatting::formatExpr(Expr::solutionSet(quadratic)),
        std::string{"{x == 1, x == -1}"},
        "PolynomialSolver: solution sets have a readable runtime representation");

    const Expr cubicEquation = Expr::call(
        builtins.symbol(evaluation::BuiltinId::Equal),
        {
            Expr::call(builtins.symbol(evaluation::BuiltinId::Power), {Expr{x}, integer(3)}),
            integer(-8)
        });
    const auto cubic = solver::solvePolynomialEquation(
        cubicEquation, x, builtins, math, angles);
    tests.expect(
        cubic.kind() == solver::SolutionSetKind::Finite
            && cubic.branches().size() == 3,
        "PolynomialSolver: perfect-cube binomial keeps all three roots");

    tests.expectThrows<std::invalid_argument>(
        [&] {
            static_cast<void>(solver::SolutionSet::finite(
                {SolverVariable{x, NumericDomain::Complex}},
                {SolutionBranch{{SolutionBinding{a, integer(1)}}, {}, std::nullopt, {}, std::nullopt}}));
        },
        "SolutionSet: a branch cannot bind a symbol outside the solver variable list");
}

} // namespace mmcal::tests
