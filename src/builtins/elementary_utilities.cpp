// abs・sign
#include "elementary_utilities.hpp"

#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "mathematics/exact_roots.hpp"
#include "mathematics/exact_trigonometry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace mmcal::builtins {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;

void requireArity(
    std::span<const Expr> arguments,
    std::size_t arity,
    std::string_view name) {
    if (arguments.size() != arity)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects " + std::to_string(arity) + " argument(s)");
}

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr pi(const mathematics::MathRegistry& mathematics) {
    const auto* definition = mathematics.findConstant(mathematics::ConstantId::Pi);
    if (!definition)
        error::throwCalcError(error::CalcErrorType::Internal, "Pi is not registered");
    return Expr{definition->symbol};
}

[[nodiscard]] Expr multiply(
    const evaluation::BuiltinRegistry& registry,
    std::vector<Expr> arguments) {
    return Expr::call(registry.symbol(BuiltinId::Multiply), std::move(arguments));
}

[[nodiscard]] Expr divide(
    const evaluation::BuiltinRegistry& registry,
    Expr numerator,
    Expr denominator) {
    return Expr::call(
        registry.symbol(BuiltinId::Divide),
        {std::move(numerator), std::move(denominator)});
}

[[nodiscard]] Expr holdUnary(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    BuiltinId id,
    std::string_view name) {
    requireArity(arguments, 1, name);
    return Expr::call(registry.symbol(id), {arguments.front()});
}

} // namespace

Expr evaluateCbrt(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireArity(arguments, 1, names::cbrt);
    static_cast<void>(mathematics);

    if (arguments.front().isNumber()) {
        const Number& value = arguments.front().asNumber();
        if (!value.isReal())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "cbrt requires a real argument");
        if (const auto root = mathematics::exactRealCubeRoot(value.asReal()))
            return Expr{Number{*root}};
    }

    return Expr::call(registry.symbol(BuiltinId::Cbrt), {arguments.front()});
}

Expr evaluateHypot(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::hypot);
    for (const Expr& argument : arguments) {
        if (argument.isNumber() && !argument.asNumber().isReal())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "hypot requires real arguments");
    }

    if (arguments[0].isNumber() && arguments[1].isNumber()) {
        const auto& x = arguments[0].asNumber().asReal();
        const auto& y = arguments[1].asNumber().asReal();
        const auto squared = x * x + y * y;
        if (const auto root = mathematics::exactSquareRoot(squared))
            return Expr{Number{*root}};
        // 無理数になる場合も、既にexactに求まった平方和へ落としてからsqrtを保持する。
        // hypot[1,1] -> sqrt[2] のように余分な x^2+y^2 構造を残さない。
        return Expr::call(registry.symbol(BuiltinId::Sqrt), {Expr{Number{squared}}});
    }

    return Expr::call(registry.symbol(BuiltinId::Hypot), {arguments[0], arguments[1]});
}

Expr evaluateCis(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    return holdUnary(arguments, registry, BuiltinId::Cis, names::cis);
}

Expr evaluatePolar(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::polar);
    // polar[r,theta] は r cis[theta] の薄いexact wrapperとして定義する。
    return multiply(
        registry,
        {arguments[0], Expr::call(registry.symbol(BuiltinId::Cis), {arguments[1]})});
}

Expr evaluateDegreeToRadian(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireArity(arguments, 1, names::degreeToRadian);
    return divide(registry,
        multiply(registry, {arguments[0], pi(mathematics)}), integer(180));
}

Expr evaluateDegreeToGradian(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::degreeToGradian);
    return divide(registry, multiply(registry, {integer(10), arguments[0]}), integer(9));
}

Expr evaluateRadianToDegree(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireArity(arguments, 1, names::radianToDegree);
    if (const auto coefficient = mathematics::extractRationalPiMultiple(
            arguments[0], registry, mathematics))
        return Expr{Number{*coefficient * numeric::Rational{BigInt{180}}}};
    return divide(registry,
        multiply(registry, {integer(180), arguments[0]}), pi(mathematics));
}

Expr evaluateRadianToGradian(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireArity(arguments, 1, names::radianToGradian);
    if (const auto coefficient = mathematics::extractRationalPiMultiple(
            arguments[0], registry, mathematics))
        return Expr{Number{*coefficient * numeric::Rational{BigInt{200}}}};
    return divide(registry,
        multiply(registry, {integer(200), arguments[0]}), pi(mathematics));
}

Expr evaluateGradianToDegree(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::gradianToDegree);
    return divide(registry, multiply(registry, {integer(9), arguments[0]}), integer(10));
}

Expr evaluateGradianToRadian(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireArity(arguments, 1, names::gradianToRadian);
    return divide(registry,
        multiply(registry, {arguments[0], pi(mathematics)}), integer(200));
}

} // namespace mmcal::builtins
