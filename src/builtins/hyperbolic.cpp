// 双曲線函数
#include "hyperbolic.hpp"
#include "builtin_helpers.hpp"



namespace mmcal::builtins {

#define MMCAL_HOLD_HYPERBOLIC(NAME, ID, TEXT) \
expression::Expr NAME( \
    std::span<const expression::Expr> arguments, \
    const evaluation::BuiltinRegistry& registry, \
    const mathematics::MathRegistry& mathematics) { \
    static_cast<void>(mathematics); \
    return holdUnaryBuiltin(arguments, registry, evaluation::BuiltinId::ID, TEXT); \
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
