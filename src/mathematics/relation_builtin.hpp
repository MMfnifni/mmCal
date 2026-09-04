#pragma once

#include "evaluation/builtin_registry.hpp"
#include "mathematics/predicate.hpp"

#include <optional>

namespace mmcal::mathematics {

// 関係演算子のBuiltinId対応はAssumptionとSolverで共通である。
// ここを一箇所に固定し，新しい関係演算子を追加した際の解釈ずれを防ぐ。
[[nodiscard]] inline std::optional<RelationKind> relationKindForBuiltin(
    evaluation::BuiltinId id) noexcept {
    using evaluation::BuiltinId;
    switch (id) {
    case BuiltinId::Equal: return RelationKind::Equal;
    case BuiltinId::NotEqual: return RelationKind::NotEqual;
    case BuiltinId::Less: return RelationKind::Less;
    case BuiltinId::LessEqual: return RelationKind::LessEqual;
    case BuiltinId::Greater: return RelationKind::Greater;
    case BuiltinId::GreaterEqual: return RelationKind::GreaterEqual;
    default: return std::nullopt;
    }
}

[[nodiscard]] inline evaluation::BuiltinId builtinForRelation(
    RelationKind relation) noexcept {
    using evaluation::BuiltinId;
    switch (relation) {
    case RelationKind::Equal: return BuiltinId::Equal;
    case RelationKind::NotEqual: return BuiltinId::NotEqual;
    case RelationKind::Less: return BuiltinId::Less;
    case RelationKind::LessEqual: return BuiltinId::LessEqual;
    case RelationKind::Greater: return BuiltinId::Greater;
    case RelationKind::GreaterEqual: return BuiltinId::GreaterEqual;
    }
    return BuiltinId::Equal;
}

} // namespace mmcal::mathematics
