#pragma once

#include "expression/expr.hpp"
#include "numeric_domain.hpp"

#include <variant>

namespace mmcal::mathematics {

enum class RelationKind {
    Equal,
    NotEqual,
    Less,
    LessEqual,
    Greater,
    GreaterEqual
};

struct RelationPredicate final {
    RelationKind relation = RelationKind::Equal;
    expression::Expr lhs;
    expression::Expr rhs;

    [[nodiscard]] bool operator==(const RelationPredicate&) const = default;
};

struct DomainPredicate final {
    expression::Expr expression;
    NumericDomain domain = NumericDomain::Unknown;

    [[nodiscard]] bool operator==(const DomainPredicate&) const = default;
};

// Solverの条件、Assumption、知識問い合わせを同じ構造で表す。
// 現時点では関係式と数体系所属だけに絞り、論理式をExprへ無理に埋め込まない。
using Predicate = std::variant<RelationPredicate, DomainPredicate>;

[[nodiscard]] inline Predicate relation(
    RelationKind kind,
    expression::Expr lhs,
    expression::Expr rhs) {
    return RelationPredicate{kind, std::move(lhs), std::move(rhs)};
}

[[nodiscard]] inline Predicate elementOf(
    expression::Expr expression,
    NumericDomain domain) {
    return DomainPredicate{std::move(expression), domain};
}

} // namespace mmcal::mathematics
