#pragma once

#include "symbolic/algebraic_number.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

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


// root[...]の構文表現とAlgebraicNumber backendの共有境界。Evaluator / Solver /
// algebraic expression bridgeで同じ係数・index変換とcanonical Root生成を使う。
[[nodiscard]] std::optional<std::vector<numeric::Rational>> rootPolynomialCoefficients(
    const expression::Expr& expression);
[[nodiscard]] std::optional<std::size_t> positiveRootIndex(
    const expression::Expr& expression);
[[nodiscard]] expression::Expr makeCanonicalRootExpression(
    const RealAlgebraicNumber& algebraic,
    const evaluation::BuiltinRegistry& builtins);
[[nodiscard]] expression::Expr makeCanonicalRootExpression(
    const ComplexAlgebraicNumber& algebraic,
    const evaluation::BuiltinRegistry& builtins);
[[nodiscard]] expression::Expr makeCanonicalAlgebraicExpression(
    const AlgebraicNumber& algebraic,
    const evaluation::BuiltinRegistry& builtins);

// Exprとして保持されているexact代数数を，可能な範囲でAlgebraicNumberへ持ち上げる。
// root[...]の内部cacheだけでなく，Rational/complex Rational，sqrt/cbrt，Phi，および
// それらからなるboundedな四則演算・小整数冪を同じexact algebraic arithmeticへ接続する。
// 数値近似から代数性を推測せず，証明できない式はnulloptを返す。
[[nodiscard]] std::optional<AlgebraicNumber> exactAlgebraicValue(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics);

} // namespace mmcal::symbolic
