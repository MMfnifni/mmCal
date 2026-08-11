// expm1・log1p・cardinal函数
#include "stable_elementary.hpp"

#include "builtins/exact_operations.hpp"
#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "mathematics/exact_trigonometry.hpp"
#include "numeric/number.hpp"

#include <optional>
#include <string>
#include <utility>

namespace mmcal::builtins {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(std::int64_t value) { return Expr{Number{BigInt{value}}}; }

[[nodiscard]] bool exactZero(const Expr& expression) {
    if (expression.isNumber())
        return expression.asNumber().isZero();
    if (!expression.isCall())
        return false;
    const auto& call = expression.asCall();
    return call.arguments.size() == 2 && call.arguments[0].isNumber()
        && call.arguments[0].asNumber().isZero() && call.arguments[1].isString();
}


[[nodiscard]] std::optional<Expr> exactCardinalTrig(
    BuiltinId id,
    const Expr& x,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto angle = mathematics::extractExactAngle(x, registry, mathematics, angles);
    if (!angle)
        return std::nullopt;
    if (angle->turns.isZero())
        return id == BuiltinId::Cosc ? integer(0) : integer(1);

    const mathematics::FunctionId trig = id == BuiltinId::Sinc
        ? mathematics::FunctionId::Sin
        : id == BuiltinId::Cosc ? mathematics::FunctionId::Cos : mathematics::FunctionId::Tan;
    auto trigValue = mathematics::simplifyExactTrig(trig, x, registry, mathematics, angles);
    if (!trigValue)
        return std::nullopt;

    Expr numerator = id == BuiltinId::Cosc
        ? exact::subtract(integer(1), std::move(*trigValue), registry, mathematics, angles)
        : std::move(*trigValue);

    const numeric::Rational piCoefficient = angle->turns * numeric::Rational{BigInt{2}};
    const numeric::Rational reciprocalCoefficient = numeric::Rational{BigInt{1}} / piCoefficient;
    Expr scaledNumerator = exact::multiply(
        {std::move(numerator), Expr{Number{reciprocalCoefficient}}}, registry, mathematics, angles);
    const mathematics::ConstantDefinition* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
    if (!pi)
        error::throwCalcError(error::CalcErrorType::Internal, "Pi is not registered");
    return exact::divide(std::move(scaledNumerator), Expr{pi->symbol}, registry, mathematics, angles);
}
[[nodiscard]] bool exactMinusOne(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isReal()
        && expression.asNumber().asReal().isInteger()
        && expression.asNumber().asReal().asInteger() == BigInt{-1};
}

} // namespace

Expr evaluateStableElementary(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() != 1)
        error::throwCalcError(error::CalcErrorType::Type, "stable elementary function expects one argument");

    const Expr& x = arguments.front();
    if (id == BuiltinId::Log1p && exactMinusOne(x))
        error::throwCalcError(error::CalcErrorType::Domain, "log1p is undefined at -1");

    if (id == BuiltinId::Sinc || id == BuiltinId::Cosc || id == BuiltinId::Tanc)
        if (auto exact = exactCardinalTrig(id, x, registry, mathematics, angles))
            return std::move(*exact);

    if (exactZero(x)) {
        switch (id) {
        case BuiltinId::Expm1:
        case BuiltinId::Log1p:
        case BuiltinId::Cosc:
            return integer(0);
        case BuiltinId::Sinc:
        case BuiltinId::Tanc:
        case BuiltinId::Sinhc:
        case BuiltinId::Tanhc:
        case BuiltinId::Expc:
            return integer(1);
        default:
            break;
        }
    }

    return Expr::call(registry.symbol(id), {x});
}

} // namespace mmcal::builtins
