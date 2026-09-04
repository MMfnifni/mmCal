#pragma once

#include "expression/expr.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/predicate.hpp"
#include "solution_set.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace mmcal::evaluation {
class BuiltinRegistry;
enum class BuiltinId;
}
namespace mmcal::mathematics {
class AngleSemantics;
class MathRegistry;
}

namespace mmcal::solver {

[[nodiscard]] expression::Expr integerExpr(std::int64_t value);

[[nodiscard]] expression::Expr builtinCall(
    const evaluation::BuiltinRegistry& builtins,
    evaluation::BuiltinId id,
    std::vector<expression::Expr> arguments);

[[nodiscard]] expression::Expr simplifyForSolve(
    expression::Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

[[nodiscard]] std::optional<mathematics::RelationKind> relationKindOf(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins);

[[nodiscard]] mathematics::RelationKind reversedRelation(
    mathematics::RelationKind relation) noexcept;

[[nodiscard]] evaluation::BuiltinId builtinForRelation(
    mathematics::RelationKind relation) noexcept;

[[nodiscard]] expression::Expr relationExpr(
    mathematics::RelationKind relation,
    expression::Expr lhs,
    expression::Expr rhs,
    const evaluation::BuiltinRegistry& builtins);


[[nodiscard]] std::optional<SolutionSet> solveIdenticalEquality(
    const expression::Expr& relation,
    std::vector<SolverVariable> variables,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] mathematics::AssumptionSet withRealVariable(
    const mathematics::AssumptionSet& assumptions,
    const expression::Symbol& variable);

[[nodiscard]] std::size_t expressionNodeCount(
    const expression::Expr& expression,
    std::size_t limit);

[[nodiscard]] bool containsBuiltinCall(
    const expression::Expr& expression,
    evaluation::BuiltinId id,
    const evaluation::BuiltinRegistry& builtins);

} // namespace mmcal::solver
