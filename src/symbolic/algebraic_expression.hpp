#pragma once

#include "symbolic/algebraic_number.hpp"

#include <optional>

namespace mmcal::evaluation {
class BuiltinRegistry;
}

namespace mmcal::expression {
class Expr;
}

namespace mmcal::mathematics {
class MathRegistry;
}

namespace mmcal::symbolic {

// Exprとして保持されているexact代数数を，可能な範囲でAlgebraicNumberへ持ち上げる。
// root[...]の内部cacheだけでなく，Rational/complex Rational，sqrt/cbrt，Phi，および
// それらからなるboundedな四則演算・小整数冪を同じexact algebraic arithmeticへ接続する。
// 数値近似から代数性を推測せず，証明できない式はnulloptを返す。
[[nodiscard]] std::optional<AlgebraicNumber> exactAlgebraicValue(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics);

} // namespace mmcal::symbolic
