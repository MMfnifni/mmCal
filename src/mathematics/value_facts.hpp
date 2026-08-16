#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "math_registry.hpp"
#include "numeric_domain.hpp"

namespace mmcal::mathematics {

class AssumptionSet;

struct ValueFacts final {
    NumericDomain domain = NumericDomain::Unknown;
    RealSign sign = RealSign::Unknown;
    bool exact = false;
    bool provablyNonReal = false;
    // domainは「属すること」の上界だけでは Integer⊂Rational を否定できない。
    // exact non-integer / irrational knowledgeを別に保持し，Element/Solver条件のFalse証明へ使う。
    bool provablyNonInteger = false;
    bool provablyNonRational = false;

    [[nodiscard]] bool isNumeric() const noexcept;
    [[nodiscard]] bool isProvablyReal() const noexcept;
    [[nodiscard]] bool isProvablyComplex() const noexcept;
    [[nodiscard]] bool isProvablyNonNegativeReal() const noexcept;
    [[nodiscard]] bool isProvablyNegativeReal() const noexcept;
};

// Exprを評価値へ落とさず、式構造とMathRegistryの知識だけから数学的domainを推論する。
// native C++ recursionは使わず、深い記号式でもOS stackを消費しない。
[[nodiscard]] ValueFacts inferValueFacts(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics);

// Assumption付きの推論。Solver/Refineは未束縛Symbolを勝手にRealとせず、x ∈ Real や x > 0 のような明示前提がある場合だけ知識を追加する。
[[nodiscard]] ValueFacts inferValueFacts(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics,
    const AssumptionSet& assumptions);

} // namespace mmcal::mathematics
