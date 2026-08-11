// exp・log系函数
#include "transcendental.hpp"

#include "error/error_message.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace mmcal::builtins {
namespace {

void requireUnary(
    std::span<const expression::Expr> arguments,
    std::string_view name) {
    if (arguments.size() != 1)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects 1 argument(s)");
}

[[nodiscard]] expression::Expr holdUnary(
    std::span<const expression::Expr> arguments,
    std::string_view name,
    evaluation::BuiltinId id,
    const evaluation::BuiltinRegistry& registry) {
    requireUnary(arguments, name);
    return expression::Expr::call(registry.symbol(id), {arguments.front()});
}

} // namespace

expression::Expr evaluateArg(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnary(arguments, "arg", evaluation::BuiltinId::Arg, registry);
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
    return holdUnary(arguments, "exp", evaluation::BuiltinId::Exp, registry);
}

} // namespace mmcal::builtins
