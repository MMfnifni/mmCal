// 安全な標準式簡約の回帰テスト
#include "simplifier_tests.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "formatting/expr_formatter.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"
#include "mathematics/numeric_domain.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "simplification/expression_cost.hpp"
#include "simplification/expression_ordering.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "solver/solution_set.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <vector>

namespace mmcal::tests {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr call(
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id,
    std::initializer_list<Expr> arguments) {
    return Expr::call(builtins.symbol(id), std::vector<Expr>{arguments});
}

[[nodiscard]] simplification::SimplificationContext context(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& math,
    mathematics::AssumptionSet assumptions = {}) {
    return simplification::SimplificationContext{
        builtins,
        math,
        mathematics::defaultAngleSemantics(),
        std::move(assumptions),
        32
    };
}

} // namespace

void runSimplifierTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const auto builtins = evaluation::BuiltinRegistry::defaults(symbols);
    const auto math = mathematics::MathRegistry::defaults(symbols, builtins);
    const simplification::Simplifier simplifier;

    const Expr x{symbols.intern("x")};
    const Expr y{symbols.intern("y")};
    const Expr z{symbols.intern("z")};

    // Addのcanonical orderは入力順やhash/pointer値に依存しないstrict total orderを使う。
    const simplification::ExpressionLess less;
    std::vector<Expr> orderedSamples{
        integer(-2), integer(0), integer(3), Expr{false}, Expr{true},
        Expr{std::string{"a"}}, x, y,
        Expr::array({2}, {x, integer(1)}),
        call(builtins, BuiltinId::Sin, {x}),
        call(builtins, BuiltinId::Add, {x, y}),
        Expr::solutionSet(solver::SolutionSet::unresolved({
            solver::SolverVariable{symbols.intern("q"), mathematics::NumericDomain::Real}}))
    };
    bool totalOrder = true;
    for (std::size_t i = 0; i < orderedSamples.size(); ++i) {
        if (less(orderedSamples[i], orderedSamples[i]))
            totalOrder = false;
        for (std::size_t j = i + 1; j < orderedSamples.size(); ++j) {
            if (orderedSamples[i] == orderedSamples[j]
                || less(orderedSamples[i], orderedSamples[j]) == less(orderedSamples[j], orderedSamples[i]))
                totalOrder = false;
        }
    }
    for (const Expr& a : orderedSamples)
        for (const Expr& b : orderedSamples)
            for (const Expr& c : orderedSamples)
                if (less(a, b) && less(b, c) && !less(a, c))
                    totalOrder = false;
    tests.expect(totalOrder,
        "Simplifier: expression ordering is a strict total order across Expr kinds");

    const auto addPermutation = [&](std::vector<Expr> terms) {
        return formatting::formatExpr(simplifier.simplify(
            Expr::call(builtins.symbol(BuiltinId::Add), std::move(terms)),
            context(builtins, math)));
    };
    const std::string canonicalAddOrder = addPermutation({z, x, integer(2), y});
    tests.expectEqual(addPermutation({y, integer(2), z, x}), canonicalAddOrder,
        "Simplifier: Add canonical order is independent of input permutation");
    tests.expectEqual(addPermutation({x, z, y, integer(2)}), canonicalAddOrder,
        "Simplifier: Add canonical order is deterministic across permutations");

    const Expr productA = simplifier.simplify(
        call(builtins, BuiltinId::Multiply, {
            call(builtins, BuiltinId::Divide, {x, y}), z}),
        context(builtins, math));
    const Expr productB = simplifier.simplify(
        call(builtins, BuiltinId::Divide, {
            call(builtins, BuiltinId::Multiply, {z, x}), y}),
        context(builtins, math));
    tests.expect(productA == productB,
        "Simplifier: product/division normal form removes tree-shape differences");

    const Expr productC = simplifier.simplify(
        call(builtins, BuiltinId::Multiply, {
            call(builtins, BuiltinId::Divide, {x, y}),
            call(builtins, BuiltinId::Divide, {integer(2), z})}),
        context(builtins, math));
    const Expr productD = simplifier.simplify(
        call(builtins, BuiltinId::Divide, {
            call(builtins, BuiltinId::Multiply, {integer(2), x}),
            call(builtins, BuiltinId::Multiply, {y, z})}),
        context(builtins, math));
    tests.expect(productC == productD,
        "Simplifier: nested symbolic products and divisions share one normal form");

    const Expr collected = simplifier.simplify(
        call(builtins, BuiltinId::Add, {
            call(builtins, BuiltinId::Multiply, {integer(2), x}),
            call(builtins, BuiltinId::Multiply, {integer(3), x})
        }),
        context(builtins, math));
    tests.expectEqual(
        formatting::formatExpr(collected),
        std::string{"5x"},
        "Simplifier: rational coefficients of identical atoms are collected");

    const Expr cancellation = simplifier.simplify(
        call(builtins, BuiltinId::Add, {
            x,
            call(builtins, BuiltinId::Negate, {x})
        }),
        context(builtins, math));
    tests.expectEqual(
        formatting::formatExpr(cancellation),
        std::string{"0"},
        "Simplifier: x + (-x) cancels structurally");

    const Expr quotient = call(builtins, BuiltinId::Divide, {x, x});
    const Expr unknownQuotient = simplifier.simplify(quotient, context(builtins, math));
    tests.expectEqual(
        formatting::formatExpr(unknownQuotient),
        std::string{"x/x"},
        "Simplifier: x/x is retained when x != 0 is unknown");

    mathematics::AssumptionSet nonzero;
    nonzero.add(mathematics::relation(
        mathematics::RelationKind::NotEqual, x, integer(0)));
    const Expr knownQuotient = simplifier.simplify(
        quotient, context(builtins, math, nonzero));
    tests.expectEqual(
        formatting::formatExpr(knownQuotient),
        std::string{"1"},
        "Simplifier: x/x -> 1 only when x != 0 is proven");

    const Expr expX = call(builtins, BuiltinId::Exp, {x});
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            call(builtins, BuiltinId::Divide, {expX, expX}),
            context(builtins, math))),
        std::string{"1"},
        "Simplifier: zero-free MathKnowledge permits exp[x]/exp[x] cancellation");

    const Expr gammaX = call(builtins, BuiltinId::Gamma, {x});
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            call(builtins, BuiltinId::Divide, {gammaX, gammaX}),
            context(builtins, math))),
        std::string{"gamma[x]/gamma[x]"},
        "Simplifier: zero-free knowledge does not erase unresolved gamma poles");

    const Expr powerZero = call(builtins, BuiltinId::Power, {x, integer(0)});
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(powerZero, context(builtins, math))),
        std::string{"x^0"},
        "Simplifier: x^0 is retained when x may be zero");
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            powerZero, context(builtins, math, nonzero))),
        std::string{"1"},
        "Simplifier: x^0 -> 1 when x != 0 is proven");

    mathematics::AssumptionSet positive;
    positive.add(mathematics::elementOf(x, mathematics::NumericDomain::Real));
    positive.add(mathematics::relation(
        mathematics::RelationKind::GreaterEqual, x, integer(0)));
    const Expr sqrtSquare = call(builtins, BuiltinId::Sqrt, {
        call(builtins, BuiltinId::Power, {x, integer(2)})
    });
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            sqrtSquare, context(builtins, math, positive))),
        std::string{"x"},
        "Simplifier: sqrt[x^2] -> x under x >= 0");
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            sqrtSquare, context(builtins, math))),
        std::string{"sqrt[x^2]"},
        "Simplifier: sqrt[x^2] is retained without a sign proof");

    const Expr sinNegative = call(builtins, BuiltinId::Sin, {
        call(builtins, BuiltinId::Negate, {x})
    });
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            sinNegative, context(builtins, math))),
        std::string{"-sin[x]"},
        "Simplifier: sin is odd without numerically evaluating x");

    const Expr cosNegative = call(builtins, BuiltinId::Cos, {
        call(builtins, BuiltinId::Negate, {x})
    });
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            cosNegative, context(builtins, math))),
        std::string{"cos[x]"},
        "Simplifier: cos is even without numerically evaluating x");

    const Expr expLog = call(builtins, BuiltinId::Exp, {
        call(builtins, BuiltinId::Log, {x})
    });
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(expLog, context(builtins, math))),
        std::string{"exp[log[x]]"},
        "Simplifier: Exp[Log[x]] is retained when x != 0 is unknown");
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            expLog, context(builtins, math, nonzero))),
        std::string{"x"},
        "Simplifier: Exp[Log[x]] -> x when x != 0 is proven");

    const Expr logExp = call(builtins, BuiltinId::Log, {
        call(builtins, BuiltinId::Exp, {x})
    });
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(logExp, context(builtins, math))),
        std::string{"log[exp[x]]"},
        "Simplifier: Log[Exp[x]] is retained for unconstrained complex x");

    mathematics::AssumptionSet real;
    real.add(mathematics::elementOf(x, mathematics::NumericDomain::Real));
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            logExp, context(builtins, math, real))),
        std::string{"x"},
        "Simplifier: Log[Exp[x]] -> x when x is provably real");


    const Expr angle = integer(1);
    const Expr sinSquare = call(builtins, BuiltinId::Power, {
        call(builtins, BuiltinId::Sin, {angle}), integer(2)
    });
    const Expr cosSquare = call(builtins, BuiltinId::Power, {
        call(builtins, BuiltinId::Cos, {angle}), integer(2)
    });
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            call(builtins, BuiltinId::Add, {sinSquare, cosSquare}),
            context(builtins, math))),
        std::string{"1"},
        "Simplifier: sin[u]^2 + cos[u]^2 uses the exact Pythagorean identity");

    const Expr symbolicSinSquare = call(builtins, BuiltinId::Power, {
        call(builtins, BuiltinId::Sin, {x}), integer(2)
    });
    const Expr symbolicCosSquare = call(builtins, BuiltinId::Power, {
        call(builtins, BuiltinId::Cos, {x}), integer(2)
    });
    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            call(builtins, BuiltinId::Add, {symbolicSinSquare, symbolicCosSquare}),
            context(builtins, math))),
        std::string{"1"},
        "Simplifier: Pythagorean identity does not require a numeric angle");

    tests.expectEqual(
        formatting::formatExpr(simplifier.simplify(
            call(builtins, BuiltinId::Add, {
                call(builtins, BuiltinId::Multiply, {integer(3), symbolicSinSquare}),
                call(builtins, BuiltinId::Multiply, {integer(3), symbolicCosSquare})
            }),
            context(builtins, math))),
        std::string{"3"},
        "Simplifier: equal exact coefficients factor through the Pythagorean identity");

    const Expr nested = call(builtins, BuiltinId::Add, {
        x,
        call(builtins, BuiltinId::Multiply, {integer(2), y})
    });
    const auto cost = simplification::measureExpressionCost(nested);
    tests.expect(
        cost.nodes == 5 && cost.leaves == 3 && cost.depth == 3,
        "Simplifier: expression cost is independent infrastructure for FullSimplify");
}

} // namespace mmcal::tests
