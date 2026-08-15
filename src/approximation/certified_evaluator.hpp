#pragma once

#include "complex_interval.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <variant>

namespace mmcal::approximation {

// certified evaluator内部の数値domain。
// Realは必要時だけComplexへ昇格し、Complexの虚部が「小さい」という理由ではRealへ戻さない。虚部区間が厳密に[0,0]と証明できた場合だけRealとして扱える。
class CertifiedValue final {
public:
    CertifiedValue(RealInterval real);
    CertifiedValue(ComplexInterval complex);

    [[nodiscard]] bool isReal() const noexcept;
    [[nodiscard]] bool isComplex() const noexcept;
    [[nodiscard]] const RealInterval& asReal() const;
    [[nodiscard]] const ComplexInterval& asComplex() const;
    [[nodiscard]] ComplexInterval toComplex() const;

private:
    std::variant<RealInterval, ComplexInterval> value_;
};

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
        EnclosureKind enclosureKind) const;
    [[nodiscard]] std::optional<CertifiedValue> encloseCall(
        const expression::CallExpr& call,
        std::size_t precisionBits,
        std::span<const CertifiedBinding> bindings,
        EnclosureKind enclosureKind) const;
};

} // namespace mmcal::approximation
