// 複素数函数
#include "complex_functions.hpp"
#include "builtin_helpers.hpp"



namespace mmcal::builtins {

expression::Expr evaluateAbs(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnaryBuiltin(arguments, registry, evaluation::BuiltinId::Abs, "abs");
}

expression::Expr evaluateSign(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnaryBuiltin(arguments, registry, evaluation::BuiltinId::Sign, "sign");
}

expression::Expr evaluateRe(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnaryBuiltin(arguments, registry, evaluation::BuiltinId::Re, "re");
}

expression::Expr evaluateIm(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnaryBuiltin(arguments, registry, evaluation::BuiltinId::Im, "im");
}

expression::Expr evaluateConj(
    std::span<const expression::Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    static_cast<void>(mathematics);
    return holdUnaryBuiltin(arguments, registry, evaluation::BuiltinId::Conj, "conj");
}

} // namespace mmcal::builtins
