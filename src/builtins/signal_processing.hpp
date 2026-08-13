#pragma once

#include "approximation/approximation_context.hpp"
#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

namespace mmcal::builtins {

// Evaluator単位で保持するexact FFT plan cache。
// Symbol identityを含むExprのtwiddleをprocess-globalへ漏らさず、同一session内のtransform間だけ再利用する。
class FourierTransformCache final {
public:
    struct Stage final {
        std::size_t length = 0;
        std::vector<expression::Expr> forwardRoots;
        std::vector<expression::Expr> inverseRoots;
    };

    struct Plan final {
        std::vector<std::size_t> bitReversed;
        std::vector<Stage> stages;
    };

    void clear() noexcept;
    [[nodiscard]] std::size_t planCount() const noexcept;
    [[nodiscard]] Plan& plan(std::size_t size);

private:
    std::unordered_map<std::size_t, Plan> plans_;
};

[[nodiscard]] expression::Expr evaluateDft(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);
[[nodiscard]] expression::Expr evaluateFft(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    FourierTransformCache& cache);
[[nodiscard]] expression::Expr evaluateIfft(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    FourierTransformCache& cache);

// N[fft[...],p] 等から呼ばれるprecision-aware backend。
// exact Exprを巨大な記号式へ展開せず、BigFloat端点のComplexInterval上で直接変換する。
// 入力をcertified数値へ写せない場合はnulloptを返し、Evaluatorは従来のexact経路へfallbackする。
[[nodiscard]] std::optional<expression::Expr> evaluateApproximateDft(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);
[[nodiscard]] std::optional<expression::Expr> evaluateApproximateFft(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);
[[nodiscard]] std::optional<expression::Expr> evaluateApproximateIfft(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context);

[[nodiscard]] expression::Expr evaluateConvolution(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

} // namespace mmcal::builtins
