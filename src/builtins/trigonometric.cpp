// 三角函数、逆三角函数
#include "trigonometric.hpp"
#include "builtin_helpers.hpp"



namespace mmcal::builtins {

#define MMCAL_HOLD_TRIG_UNARY(NAME, ID, TEXT) \
expression::Expr NAME( \
    std::span<const expression::Expr> arguments, \
    const evaluation::BuiltinRegistry& registry, \
    const mathematics::MathRegistry& mathematics, \
    const mathematics::AngleSemantics& angleSemantics) { \
    static_cast<void>(mathematics); \
    static_cast<void>(angleSemantics); \
    return holdUnaryBuiltin(arguments, registry, evaluation::BuiltinId::ID, TEXT); \
}

MMCAL_HOLD_TRIG_UNARY(evaluateSin, Sin, "sin")
MMCAL_HOLD_TRIG_UNARY(evaluateCos, Cos, "cos")
MMCAL_HOLD_TRIG_UNARY(evaluateTan, Tan, "tan")
MMCAL_HOLD_TRIG_UNARY(evaluateCot, Cot, "cot")
MMCAL_HOLD_TRIG_UNARY(evaluateSec, Sec, "sec")
MMCAL_HOLD_TRIG_UNARY(evaluateCsc, Csc, "csc")
MMCAL_HOLD_TRIG_UNARY(evaluateAsin, Asin, "asin")
MMCAL_HOLD_TRIG_UNARY(evaluateAcos, Acos, "acos")
MMCAL_HOLD_TRIG_UNARY(evaluateAtan, Atan, "atan")

#undef MMCAL_HOLD_TRIG_UNARY

expression::Expr evaluateAtan2(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angleSemantics) {
    static_cast<void>(mathematics);
    static_cast<void>(angleSemantics);
    return holdBuiltin(arguments, registry, evaluation::BuiltinId::Atan2, "atan2", 2);
}

} // namespace mmcal::builtins
