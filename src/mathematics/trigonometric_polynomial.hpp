#pragma once

#include "expression/expr.hpp"

#include <cstddef>
#include <optional>

namespace mmcal::evaluation {
class BuiltinRegistry;
}

namespace mmcal::mathematics {

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
