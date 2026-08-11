// 三角函数、逆三角函数
#include "trigonometric.hpp"

#include "error/error_message.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace mmcal::builtins {
namespace {

void requireArity(
    std::span<const expression::Expr> arguments,
    std::size_t arity,
    std::string_view name) {
    if (arguments.size() != arity)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects " + std::to_string(arity) + " argument(s)");
}

[[nodiscard]] expression::Expr holdFunction(
    evaluation::BuiltinId id,
    std::string_view name,
    std::span<const expression::Expr> arguments,
    std::size_t arity,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, arity, name);
    return expression::Expr::call(
        registry.symbol(id),
        std::vector<expression::Expr>{arguments.begin(), arguments.end()});
}

} // namespace

#define MMCAL_HOLD_TRIG_UNARY(NAME, ID, TEXT) \
expression::Expr NAME( \
    std::span<const expression::Expr> arguments, \
    const evaluation::BuiltinRegistry& registry, \
    const mathematics::MathRegistry& mathematics, \
    const mathematics::AngleSemantics& angleSemantics) { \
    static_cast<void>(mathematics); \
    static_cast<void>(angleSemantics); \
    return holdFunction(evaluation::BuiltinId::ID, TEXT, arguments, 1, registry); \
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
    return holdFunction(evaluation::BuiltinId::Atan2, "atan2", arguments, 2, registry);
}

} // namespace mmcal::builtins
