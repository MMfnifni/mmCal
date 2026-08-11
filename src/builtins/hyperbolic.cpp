// 双曲線函数
#include "hyperbolic.hpp"

#include "error/error_message.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace mmcal::builtins {
namespace {

[[nodiscard]] expression::Expr holdUnary(
    evaluation::BuiltinId id,
    std::string_view name,
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    if (arguments.size() != 1)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects 1 argument(s)");
    return expression::Expr::call(registry.symbol(id), {arguments.front()});
}

} // namespace

#define MMCAL_HOLD_HYPERBOLIC(NAME, ID, TEXT) \
expression::Expr NAME( \
    std::span<const expression::Expr> arguments, \
    const evaluation::BuiltinRegistry& registry, \
    const mathematics::MathRegistry& mathematics) { \
    static_cast<void>(mathematics); \
    return holdUnary(evaluation::BuiltinId::ID, TEXT, arguments, registry); \
}

MMCAL_HOLD_HYPERBOLIC(evaluateSinh, Sinh, "sinh")
MMCAL_HOLD_HYPERBOLIC(evaluateCosh, Cosh, "cosh")
MMCAL_HOLD_HYPERBOLIC(evaluateTanh, Tanh, "tanh")
MMCAL_HOLD_HYPERBOLIC(evaluateAsinh, Asinh, "asinh")
MMCAL_HOLD_HYPERBOLIC(evaluateAcosh, Acosh, "acosh")
MMCAL_HOLD_HYPERBOLIC(evaluateAtanh, Atanh, "atanh")
MMCAL_HOLD_HYPERBOLIC(evaluateCsch, Csch, "csch")
MMCAL_HOLD_HYPERBOLIC(evaluateSech, Sech, "sech")
MMCAL_HOLD_HYPERBOLIC(evaluateCoth, Coth, "coth")

#undef MMCAL_HOLD_HYPERBOLIC

} // namespace mmcal::builtins
