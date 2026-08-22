// exp・log系函数
#include "transcendental.hpp"
#include "builtin_helpers.hpp"

#include "error/error_message.hpp"

#include <vector>

namespace mmcal::builtins {

expression::Expr evaluateArg(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnaryBuiltin(arguments, registry, evaluation::BuiltinId::Arg, "arg");
}

expression::Expr evaluateLog(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    if (arguments.size() < 1 || arguments.size() > 2)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "log expects 1 or 2 argument(s)");
    return expression::Expr::call(
        registry.symbol(evaluation::BuiltinId::Log),
        std::vector<expression::Expr>{arguments.begin(), arguments.end()});
}

expression::Expr evaluateExp(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnaryBuiltin(arguments, registry, evaluation::BuiltinId::Exp, "exp");
}

} // namespace mmcal::builtins
