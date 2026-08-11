#pragma once

#include "assumption_set.hpp"
#include "truth_value.hpp"
#include "value_facts.hpp"

namespace mmcal::mathematics {

// MathRegistryの恒久的な数学知識と、その場限りのAssumptionをまとめて問い合わせる窓口。
// Solver/Simplifierが独自に「xは実数だったはず」と状態を持たないための共通層。
class KnowledgeContext final {
public:
    KnowledgeContext(
        const evaluation::BuiltinRegistry& builtins,
        const MathRegistry& mathematics,
        const AssumptionSet& assumptions);

    [[nodiscard]] ValueFacts facts(const expression::Expr& expression) const;
    [[nodiscard]] TruthValue prove(const Predicate& predicate) const;

private:
    const evaluation::BuiltinRegistry& builtins_;
    const MathRegistry& mathematics_;
    const AssumptionSet& assumptions_;
};

} // namespace mmcal::mathematics
