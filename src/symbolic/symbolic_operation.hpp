// 記号演算の共通結果型
#pragma once

namespace mmcal::symbolic {

// 公開予定の記号処理を一つの巨大Simplifierへ混ぜないための分類。
// Simplify/FullSimplifyは同値簡約、Expand/Factor/Collectは指示された正規形変換、
// SolveはSolutionSetを返す探索処理であり、それぞれ別アルゴリズムとして実装する。
enum class SymbolicOperation {
    Simplify,
    FullSimplify,
    Expand,
    Factor,
    Collect,
    Solve
};

enum class SymbolicResultKind {
    Expression,
    SolutionSet
};

[[nodiscard]] constexpr SymbolicResultKind resultKind(SymbolicOperation operation) noexcept {
    return operation == SymbolicOperation::Solve
        ? SymbolicResultKind::SolutionSet
        : SymbolicResultKind::Expression;
}

} // namespace mmcal::symbolic
