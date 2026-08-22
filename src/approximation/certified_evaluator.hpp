#pragma once

#include "certified_value.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"

#include <cstddef>
#include <optional>
#include <span>

namespace mmcal::approximation {

// Expr全体をcertifiedな実/複素区間へ写す数値評価器。
//
// 通常Evaluatorが exact simplification を担当し、このクラスは N[...] が明示されたときだけ呼ばれる。
// したがって Pi や sqrt[2] のexact identityを壊さず、必要な瞬間に初めてBigFloat/RealIntervalへ変換する。
struct CertifiedBinding final {
    expression::Symbol symbol;
    CertifiedValue value;
};

class CertifiedEvaluator final {
public:
    enum class EnclosureKind {
        Certified,
        Information
    };

    CertifiedEvaluator(
        const evaluation::BuiltinRegistry& builtins,
        const mathematics::MathRegistry& mathematics,
        const mathematics::AngleSemantics& angleSemantics);

    // 現在実装済みの数学函数・定数だけを評価する。
    // 未対応のSymbol/函数/branchに遭遇した場合はnulloptを返し、呼出側は元のExprを保持する。
    // 再帰評価は内部depth budgetでも制限し、病的な式でOSのstack overflowへ到達しない。
    [[nodiscard]] std::optional<CertifiedValue> enclose(
        const expression::Expr& expression,
        std::size_t precisionBits,
        EnclosureKind enclosureKind = EnclosureKind::Certified) const;

    // 数値解析用。自由変数をcertified区間へ束縛したまま式全体を評価する。
    // 点束縛だけでなく区間束縛を許すため、積分の導函数上界なども同じ評価器で扱える。
    [[nodiscard]] std::optional<CertifiedValue> enclose(
        const expression::Expr& expression,
        std::size_t precisionBits,
        std::span<const CertifiedBinding> bindings,
        EnclosureKind enclosureKind = EnclosureKind::Certified) const;

private:
    const evaluation::BuiltinRegistry& builtins_;
    const mathematics::MathRegistry& mathematics_;
    const mathematics::AngleSemantics& angleSemantics_;

    [[nodiscard]] std::optional<CertifiedValue> encloseBound(
        const expression::Expr& expression,
        std::size_t precisionBits,
        std::span<const CertifiedBinding> bindings,
        EnclosureKind enclosureKind,
        std::size_t recursionDepth) const;
    [[nodiscard]] std::optional<CertifiedValue> encloseCall(
        const expression::CallExpr& call,
        std::size_t precisionBits,
        std::span<const CertifiedBinding> bindings,
        EnclosureKind enclosureKind,
        std::size_t recursionDepth) const;
};

} // namespace mmcal::approximation
