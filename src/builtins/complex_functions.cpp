// 複素数函数
#include "complex_functions.hpp"

#include "error/error_message.hpp"

#include <string>
#include <string_view>

namespace mmcal::builtins {
namespace {

[[nodiscard]] expression::Expr holdUnary(
    std::span<const expression::Expr> arguments,
    std::string_view name,
    evaluation::BuiltinId id,
    const evaluation::BuiltinRegistry& registry) {
    if (arguments.size() != 1)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects 1 argument(s)");
    return expression::Expr::call(registry.symbol(id), {arguments.front()});
}

} // namespace

expression::Expr evaluateAbs(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnary(arguments, "abs", evaluation::BuiltinId::Abs, registry);
}

expression::Expr evaluateSign(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnary(arguments, "sign", evaluation::BuiltinId::Sign, registry);
}

expression::Expr evaluateRe(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnary(arguments, "re", evaluation::BuiltinId::Re, registry);
}

expression::Expr evaluateIm(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnary(arguments, "im", evaluation::BuiltinId::Im, registry);
}

expression::Expr evaluateConj(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnary(arguments, "conj", evaluation::BuiltinId::Conj, registry);
}

} // namespace mmcal::builtins
