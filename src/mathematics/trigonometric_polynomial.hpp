#pragma once

#include "expression/expr.hpp"
#include "numeric/rational.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace mmcal::evaluation {
class BuiltinRegistry;
}

namespace mmcal::mathematics {

struct TrigFourierTerm final {
    numeric::Rational coefficient;
    std::size_t frequency = 0;
    bool sine = false;
};

struct TrigFourierExpansion final {
    expression::Expr argument;
    std::vector<TrigFourierTerm> terms;
    std::size_t totalDegree = 0;
};

// sin[u]^m cos[u]^n の有限Fourier係数を厳密に生成する。
// frequency=0 は定数項を表す。Expr化せず積分器等から直接利用できる。
[[nodiscard]] expression::Expr scaledTrigArgumentForFrequency(
    const expression::Expr& source,
    std::size_t multiplier,
    const evaluation::BuiltinRegistry& builtins);

[[nodiscard]] std::optional<TrigFourierExpansion> expandTrigMonomial(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    std::size_t maximumTotalDegree = 64);

// sin[u]^m cos[u]^n を有限Fourier和へ厳密変換する。
// 恒等式として大域的に成立する変換だけを扱い、積分器や将来のtrigReduceから共有する。
[[nodiscard]] std::optional<expression::Expr> reduceTrigMonomial(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    std::size_t maximumTotalDegree = 64);

// 1次のsin/cos 2因子を積和公式で厳密変換する。
// 明示角度単位がある場合は両因子の単位が同一のときだけ変換する。
[[nodiscard]] std::optional<expression::Expr> reduceTrigProduct(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins);

} // namespace mmcal::mathematics
