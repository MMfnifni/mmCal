// 式全体の保証付き区間評価
#include "certified_evaluator.hpp"

#include "certified_constants.hpp"
#include "certification_error.hpp"
#include "certified_atan.hpp"
#include "certified_complex_transcendental.hpp"
#include "certified_exponential.hpp"
#include "certified_elementary_functions.hpp"
#include "certified_logarithm.hpp"
#include "certified_complex_sqrt.hpp"
#include "certified_sqrt.hpp"
#include "certified_special_functions.hpp"
#include "certified_trigonometry.hpp"
#include "interval_math.hpp"
#include "mathematics/exact_trigonometry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <charconv>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace mmcal::approximation {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::AngleUnit;
using mathematics::ConstantId;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

constexpr std::size_t maximumCertifiedExpressionDepth = 96;

[[nodiscard]] bool exceedsCertifiedExpressionDepth(const Expr& root) {
    struct Pending final {
        const Expr* expression = nullptr;
        std::size_t depth = 0;
    };

    std::vector<Pending> pending;
    pending.push_back(Pending{&root, 1});
    while (!pending.empty()) {
        const Pending current = pending.back();
        pending.pop_back();
        if (current.depth > maximumCertifiedExpressionDepth)
            return true;

        if (current.expression->isCall()) {
            const auto& arguments = current.expression->asCall().arguments;
            for (const Expr& argument : arguments)
                pending.push_back(Pending{&argument, current.depth + 1});
        }
        else if (current.expression->isArray()) {
            const auto& elements = current.expression->asArray().elements;
            for (const Expr& element : elements)
                pending.push_back(Pending{&element, current.depth + 1});
        }
    }
    return false;
}

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

[[nodiscard]] RealInterval exactRealInterval(
    const RealNumber& value,
    std::size_t precisionBits) {
    return RealInterval::fromRational(value.toRational(), precisionBits);
}

[[nodiscard]] CertifiedValue normalizeComplex(ComplexInterval value) {
    if (value.isProvablyReal())
        return CertifiedValue{value.real()};
    return CertifiedValue{std::move(value)};
}

[[nodiscard]] CertifiedValue addValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{approximation::add(lhs.asReal(), rhs.asReal(), precisionBits)};

    return normalizeComplex(approximation::add(
        lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue subtractValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{approximation::subtract(lhs.asReal(), rhs.asReal(), precisionBits)};

    return normalizeComplex(approximation::subtract(
        lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue multiplyValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{approximation::multiply(lhs.asReal(), rhs.asReal(), precisionBits)};

    return normalizeComplex(approximation::multiply(
        lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue divideValues(
    const CertifiedValue& lhs,
    const CertifiedValue& rhs,
    std::size_t precisionBits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{approximation::divide(lhs.asReal(), rhs.asReal(), precisionBits)};

    return normalizeComplex(approximation::divide(
        lhs.toComplex(), rhs.toComplex(), precisionBits));
}

[[nodiscard]] CertifiedValue negateValue(const CertifiedValue& value) {
    if (value.isReal())
        return CertifiedValue{approximation::negate(value.asReal())};
    return normalizeComplex(approximation::negate(value.asComplex()));
}

[[nodiscard]] std::optional<std::uint64_t> toUint64(const BigInt& value) {
    if (value.isNegative())
        return std::nullopt;

    const std::string text = value.toString();
    std::uint64_t result = 0;
    const auto conversion = std::from_chars(
        text.data(), text.data() + text.size(), result);
    if (conversion.ec != std::errc{} || conversion.ptr != text.data() + text.size())
        return std::nullopt;
    return result;
}

[[nodiscard]] std::optional<Rational> exactRealRational(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    return expression.asNumber().asReal().toRational();
}

[[nodiscard]] bool isOneHalf(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return false;

    return expression.asNumber().asReal().toRational() == rational(1, 2);
}

[[nodiscard]] std::optional<BigInt> exactIntegerExponent(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    return expression.asNumber().asReal().asInteger();
}

[[nodiscard]] CertifiedValue integerPower(
    CertifiedValue base,
    std::uint64_t exponent,
    std::size_t precisionBits) {
    const RealInterval one = RealInterval::fromRational(rational(1), precisionBits);
    CertifiedValue result{one};

    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = multiplyValues(result, base, precisionBits);

        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyValues(base, base, precisionBits);
    }

    return result;
}

[[nodiscard]] CertifiedValue principalSqrtReal(
    const RealInterval& value,
    std::size_t precisionBits) {
    const numeric::BigFloat zero;

    // 全区間が非負なら通常の実平方根。
    if (value.lower() >= zero)
        return CertifiedValue{encloseSqrt(value, precisionBits).interval};

    // 全区間が負なら principal sqrt は純虚数。z in [a,b] < 0 なら sqrt(z) = I*sqrt(-z) であり、-z は [-b,-a] > 0。
    if (value.upper() < zero) {
        const RealInterval magnitude = approximation::negate(value);
        const RealInterval imaginary = encloseSqrt(magnitude, precisionBits).interval;
        return CertifiedValue{ComplexInterval{
            RealInterval::fromRational(rational(0), precisionBits),
            imaginary
        }};
    }

    // 区間が0を跨ぐ場合、真値が負か正かを現precisionでは一意に決められない。
    // principal sqrt の像は負側では正の虚軸、正側では正の実軸に乗る。その両方を含む第1象限の長方形を返せば包含は保証できる。
    const RealInterval zeroInterval = RealInterval::fromRational(rational(0), precisionBits);
    const RealInterval positiveInput{zeroInterval.lower(), value.upper()};
    const RealInterval negativeInput{value.lower(), zeroInterval.upper()};
    const RealInterval realPart = encloseSqrt(positiveInput, precisionBits).interval;
    const RealInterval imaginaryPart = encloseSqrt(
        approximation::negate(negativeInput), precisionBits).interval;

    return CertifiedValue{ComplexInterval{
        RealInterval{zeroInterval.lower(), realPart.upper()},
        RealInterval{zeroInterval.lower(), imaginaryPart.upper()}
    }};
}

struct AngleOperand final {
    Expr value;
    AngleUnit unit = AngleUnit::Radian;
};

[[nodiscard]] std::optional<AngleOperand> splitAngleOperand(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::AngleSemantics& angleSemantics) {
    const auto* definition = expression.isCall()
        ? builtins.find(expression.asCall().head)
        : nullptr;

    if (!definition || definition->id != BuiltinId::UnitApplied)
        return AngleOperand{expression, angleSemantics.defaultUnit()};

    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() != 2 || !arguments[1].isString())
        return std::nullopt;

    const auto unit = mathematics::AngleSemantics::parseUnit(arguments[1].asString());
    if (!unit)
        return std::nullopt;

    return AngleOperand{arguments[0], *unit};
}

[[nodiscard]] CertifiedValue scaleValueByReal(
    const CertifiedValue& value,
    const RealInterval& factor,
    std::size_t precisionBits) {
    if (value.isReal())
        return CertifiedValue{approximation::multiply(value.asReal(), factor, precisionBits)};
    return normalizeComplex(approximation::multiply(
        value.asComplex(), ComplexInterval::fromReal(factor), precisionBits));
}

[[nodiscard]] CertifiedValue angleValueToRadians(
    const CertifiedValue& value,
    AngleUnit unit,
    const mathematics::MathRegistry& mathematics,
    std::size_t precisionBits) {
    if (unit == AngleUnit::Radian)
        return value;

    const auto* piDefinition = mathematics.findConstant(ConstantId::Pi);
    if (!piDefinition)
        throw std::logic_error("Pi is not registered in MathRegistry");
    const auto pi = encloseConstant(piDefinition->id, precisionBits);
    if (!pi)
        throw std::logic_error("Pi does not have a certified provider");

    const Rational divisor = unit == AngleUnit::Degree ? rational(180) : rational(200);
    const RealInterval factor = approximation::divide(
        pi->interval,
        RealInterval::fromRational(divisor, precisionBits),
        precisionBits);
    return scaleValueByReal(value, factor, precisionBits);
}

[[nodiscard]] CertifiedValue radiansToAngleValue(
    const CertifiedValue& value,
    AngleUnit unit,
    const mathematics::MathRegistry& mathematics,
    std::size_t precisionBits) {
    if (unit == AngleUnit::Radian)
        return value;

    const auto* piDefinition = mathematics.findConstant(ConstantId::Pi);
    if (!piDefinition)
        throw std::logic_error("Pi is not registered in MathRegistry");
    const auto pi = encloseConstant(piDefinition->id, precisionBits);
    if (!pi)
        throw std::logic_error("Pi does not have a certified provider");

    const Rational multiplier = unit == AngleUnit::Degree ? rational(180) : rational(200);
    const RealInterval factor = approximation::divide(
        RealInterval::fromRational(multiplier, precisionBits),
        pi->interval,
        precisionBits);
    return scaleValueByReal(value, factor, precisionBits);
}

[[nodiscard]] CertifiedValue reciprocalValue(
    const CertifiedValue& value,
    std::size_t precisionBits,
    const char* message) {
    if (value.isReal()) {
        if (value.asReal().containsZero())
            throw PrecisionInsufficient{message};
        return CertifiedValue{approximation::divide(
            RealInterval::fromRational(rational(1), precisionBits),
            value.asReal(), precisionBits)};
    }

    if (value.asComplex().containsZero())
        throw PrecisionInsufficient{message};
    return normalizeComplex(approximation::divide(
        ComplexInterval::fromReal(RealInterval::fromRational(rational(1), precisionBits)),
        value.asComplex(), precisionBits));
}

} // namespace

CertifiedValue::CertifiedValue(RealInterval real)
    : value_(std::move(real)) {}

CertifiedValue::CertifiedValue(ComplexInterval complex)
    : value_(std::move(complex)) {}

bool CertifiedValue::isReal() const noexcept {
    return std::holds_alternative<RealInterval>(value_);
}

bool CertifiedValue::isComplex() const noexcept {
    return std::holds_alternative<ComplexInterval>(value_);
}

const RealInterval& CertifiedValue::asReal() const {
    if (!isReal())
        throw std::logic_error("CertifiedValue does not contain a real interval");
    return std::get<RealInterval>(value_);
}

const ComplexInterval& CertifiedValue::asComplex() const {
    if (!isComplex())
        throw std::logic_error("CertifiedValue does not contain a complex interval");
    return std::get<ComplexInterval>(value_);
}

ComplexInterval CertifiedValue::toComplex() const {
    if (isComplex())
        return asComplex();
    return ComplexInterval::fromReal(asReal());
}

CertifiedEvaluator::CertifiedEvaluator(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angleSemantics)
    : builtins_(builtins),
      mathematics_(mathematics),
      angleSemantics_(angleSemantics) {}

std::optional<CertifiedValue> CertifiedEvaluator::enclose(
    const Expr& expression,
    std::size_t precisionBits) const {
    if (exceedsCertifiedExpressionDepth(expression))
        return std::nullopt;
    return encloseBound(expression, precisionBits, {});
}

std::optional<CertifiedValue> CertifiedEvaluator::enclose(
    const Expr& expression,
    std::size_t precisionBits,
    std::span<const CertifiedBinding> bindings) const {
    if (exceedsCertifiedExpressionDepth(expression))
        return std::nullopt;
    return encloseBound(expression, precisionBits, bindings);
}

std::optional<CertifiedValue> CertifiedEvaluator::encloseBound(
    const Expr& expression,
    std::size_t precisionBits,
    std::span<const CertifiedBinding> bindings) const {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified evaluation precision must be at least one bit");

    if (expression.isNumber()) {
        const Number& number = expression.asNumber();
        if (number.isReal())
            return CertifiedValue{exactRealInterval(number.asReal(), precisionBits)};

        const auto& complex = number.asComplex();
        return CertifiedValue{ComplexInterval{
            exactRealInterval(complex.real, precisionBits),
            exactRealInterval(complex.imaginary, precisionBits)
        }};
    }

    if (expression.isSymbol()) {
        for (const CertifiedBinding& binding : bindings) {
            if (binding.symbol.sameIdentity(expression.asSymbol()))
                return binding.value;
        }

        const auto* constant = mathematics_.findConstant(expression.asSymbol());
        if (!constant)
            return std::nullopt;

        const auto enclosed = encloseConstant(constant->id, precisionBits);
        if (!enclosed)
            return std::nullopt;
        return CertifiedValue{enclosed->interval};
    }

    if (expression.isCall())
        return encloseCall(expression.asCall(), precisionBits, bindings);

    return std::nullopt;
}

std::optional<CertifiedValue> CertifiedEvaluator::encloseCall(
    const expression::CallExpr& call,
    std::size_t precisionBits,
    std::span<const CertifiedBinding> bindings) const {
    const auto* definition = builtins_.find(call.head);
    if (!definition)
        return std::nullopt;

    const auto encloseArgument = [&](std::size_t index) -> std::optional<CertifiedValue> {
        if (index >= call.arguments.size())
            return std::nullopt;
        return encloseBound(call.arguments[index], precisionBits, bindings);
    };

    switch (definition->id) {
    case BuiltinId::Add: {
        CertifiedValue result{RealInterval::fromRational(rational(0), precisionBits)};
        for (const Expr& argument : call.arguments) {
            const auto enclosed = encloseBound(argument, precisionBits, bindings);
            if (!enclosed)
                return std::nullopt;
            result = addValues(result, *enclosed, precisionBits);
        }
        return result;
    }

    case BuiltinId::Subtract: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        const auto lhs = encloseArgument(0);
        const auto rhs = encloseArgument(1);
        if (!lhs || !rhs)
            return std::nullopt;
        return subtractValues(*lhs, *rhs, precisionBits);
    }

    case BuiltinId::Multiply: {
        CertifiedValue result{RealInterval::fromRational(rational(1), precisionBits)};
        for (const Expr& argument : call.arguments) {
            const auto enclosed = encloseBound(argument, precisionBits, bindings);
            if (!enclosed)
                return std::nullopt;
            result = multiplyValues(result, *enclosed, precisionBits);
        }
        return result;
    }

    case BuiltinId::Divide: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        const auto lhs = encloseArgument(0);
        const auto rhs = encloseArgument(1);
        if (!lhs || !rhs)
            return std::nullopt;
        return divideValues(*lhs, *rhs, precisionBits);
    }

    case BuiltinId::Negate: {
        const auto value = encloseArgument(0);
        return value ? std::optional<CertifiedValue>{negateValue(*value)} : std::nullopt;
    }

    case BuiltinId::Cbrt: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        if (!value->isReal())
            throw std::domain_error("cbrt requires a real argument");
        return CertifiedValue{encloseRealCubeRoot(value->asReal(), precisionBits)};
    }

    case BuiltinId::Hypot: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        const auto x = encloseArgument(0);
        const auto y = encloseArgument(1);
        if (!x || !y)
            return std::nullopt;
        if (!x->isReal() || !y->isReal())
            throw std::domain_error("hypot requires real arguments");
        const RealInterval sum = approximation::add(
            squareInterval(x->asReal(), precisionBits),
            squareInterval(y->asReal(), precisionBits),
            precisionBits);
        return CertifiedValue{encloseSqrt(sum, precisionBits).interval};
    }

    case BuiltinId::Cis: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto angleOperand = splitAngleOperand(
            call.arguments[0], builtins_, angleSemantics_);
        if (!angleOperand)
            return std::nullopt;
        const auto scalar = encloseBound(angleOperand->value, precisionBits, bindings);
        if (!scalar)
            return std::nullopt;
        const CertifiedValue radians = angleValueToRadians(
            *scalar, angleOperand->unit, mathematics_, precisionBits);
        const ComplexInterval value = radians.toComplex();
        const ComplexInterval sine = encloseComplexSinRadian(value, precisionBits);
        const ComplexInterval cosine = encloseComplexCosRadian(value, precisionBits);
        const ComplexInterval iSine{
            approximation::negate(sine.imaginary()),
            sine.real()};
        return normalizeComplex(approximation::add(cosine, iSine, precisionBits));
    }

    case BuiltinId::DegreeToRadian:
    case BuiltinId::DegreeToGradian:
    case BuiltinId::RadianToDegree:
    case BuiltinId::RadianToGradian:
    case BuiltinId::GradianToDegree:
    case BuiltinId::GradianToRadian: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        switch (definition->id) {
        case BuiltinId::DegreeToRadian:
            return angleValueToRadians(*value, AngleUnit::Degree, mathematics_, precisionBits);
        case BuiltinId::DegreeToGradian:
            return radiansToAngleValue(
                angleValueToRadians(*value, AngleUnit::Degree, mathematics_, precisionBits),
                AngleUnit::Gradian, mathematics_, precisionBits);
        case BuiltinId::RadianToDegree:
            return radiansToAngleValue(*value, AngleUnit::Degree, mathematics_, precisionBits);
        case BuiltinId::RadianToGradian:
            return radiansToAngleValue(*value, AngleUnit::Gradian, mathematics_, precisionBits);
        case BuiltinId::GradianToDegree:
            return radiansToAngleValue(
                angleValueToRadians(*value, AngleUnit::Gradian, mathematics_, precisionBits),
                AngleUnit::Degree, mathematics_, precisionBits);
        case BuiltinId::GradianToRadian:
            return angleValueToRadians(*value, AngleUnit::Gradian, mathematics_, precisionBits);
        default:
            return std::nullopt;
        }
    }

    case BuiltinId::Expm1: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        const CertifiedValue one{RealInterval::fromRational(rational(1), precisionBits)};
        if (value->isReal())
            return subtractValues(
                CertifiedValue{encloseExp(value->asReal(), precisionBits).interval},
                one, precisionBits);
        return subtractValues(
            normalizeComplex(encloseComplexExp(value->asComplex(), precisionBits).interval),
            one, precisionBits);
    }

    case BuiltinId::Log1p: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        const CertifiedValue shifted = addValues(
            *value, CertifiedValue{RealInterval::fromRational(rational(1), precisionBits)},
            precisionBits);
        if (shifted.isReal()) {
            const numeric::BigFloat zero;
            if (shifted.asReal().lower() > zero)
                return CertifiedValue{encloseLogPositive(shifted.asReal(), precisionBits).interval};
        }
        return normalizeComplex(enclosePrincipalComplexLog(
            shifted.toComplex(), precisionBits).interval);
    }

    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto angleOperand = splitAngleOperand(call.arguments[0], builtins_, angleSemantics_);
        if (!angleOperand)
            return std::nullopt;
        const auto scalar = encloseBound(angleOperand->value, precisionBits, bindings);
        if (!scalar)
            return std::nullopt;
        const CertifiedValue radians = angleValueToRadians(
            *scalar, angleOperand->unit, mathematics_, precisionBits);
        const RealInterval oneReal = RealInterval::fromRational(rational(1), precisionBits);

        if (radians.isReal()) {
            const RealInterval& x = radians.asReal();
            if (x.isPoint() && x.lower().isZero()) {
                return CertifiedValue{RealInterval::fromRational(
                    definition->id == BuiltinId::Cosc ? rational(0) : rational(1), precisionBits)};
            }
            if (x.containsZero())
                throw PrecisionInsufficient{"Cardinal trigonometric denominator cannot yet be proven nonzero"};
            const auto sine = encloseSinRadianInterval(x, precisionBits).interval;
            const auto cosine = encloseCosRadianInterval(x, precisionBits).interval;
            if (definition->id == BuiltinId::Sinc)
                return CertifiedValue{approximation::divide(sine, x, precisionBits)};
            if (definition->id == BuiltinId::Cosc)
                return CertifiedValue{approximation::divide(
                    approximation::subtract(oneReal, cosine, precisionBits), x, precisionBits)};
            if (cosine.containsZero())
                throw PrecisionInsufficient{"tanc pole cannot yet be excluded"};
            return CertifiedValue{approximation::divide(
                sine, approximation::multiply(x, cosine, precisionBits), precisionBits)};
        }

        const ComplexInterval x = radians.asComplex();
        const bool pointZero = x.real().isPoint() && x.real().lower().isZero()
            && x.imaginary().isPoint() && x.imaginary().lower().isZero();
        if (pointZero)
            return CertifiedValue{RealInterval::fromRational(
                definition->id == BuiltinId::Cosc ? rational(0) : rational(1), precisionBits)};
        if (x.containsZero())
            throw PrecisionInsufficient{"Complex cardinal trigonometric denominator cannot yet be proven nonzero"};
        const ComplexInterval sine = encloseComplexSinRadian(x, precisionBits);
        const ComplexInterval cosine = encloseComplexCosRadian(x, precisionBits);
        if (definition->id == BuiltinId::Sinc)
            return normalizeComplex(approximation::divide(sine, x, precisionBits));
        if (definition->id == BuiltinId::Cosc) {
            const ComplexInterval one = ComplexInterval::fromReal(oneReal);
            return normalizeComplex(approximation::divide(
                approximation::subtract(one, cosine, precisionBits), x, precisionBits));
        }
        if (cosine.containsZero())
            throw PrecisionInsufficient{"Complex tanc pole cannot yet be excluded"};
        return normalizeComplex(approximation::divide(
            sine, approximation::multiply(x, cosine, precisionBits), precisionBits));
    }

    case BuiltinId::Sinhc:
    case BuiltinId::Tanhc:
    case BuiltinId::Expc: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        const bool pointZero = value->isReal()
            ? value->asReal().isPoint() && value->asReal().lower().isZero()
            : value->asComplex().real().isPoint() && value->asComplex().real().lower().isZero()
                && value->asComplex().imaginary().isPoint() && value->asComplex().imaginary().lower().isZero();
        if (pointZero)
            return CertifiedValue{RealInterval::fromRational(rational(1), precisionBits)};
        if ((value->isReal() && value->asReal().containsZero())
            || (value->isComplex() && value->asComplex().containsZero()))
            throw PrecisionInsufficient{"Cardinal function denominator cannot yet be proven nonzero"};

        if (definition->id == BuiltinId::Expc) {
            const CertifiedValue one{RealInterval::fromRational(rational(1), precisionBits)};
            const CertifiedValue exponential = value->isReal()
                ? CertifiedValue{encloseExp(value->asReal(), precisionBits).interval}
                : normalizeComplex(encloseComplexExp(value->asComplex(), precisionBits).interval);
            return divideValues(subtractValues(exponential, one, precisionBits), *value, precisionBits);
        }

        if (value->isReal()) {
            const RealInterval sinhValue = encloseSinhReal(value->asReal(), precisionBits);
            if (definition->id == BuiltinId::Sinhc)
                return CertifiedValue{approximation::divide(sinhValue, value->asReal(), precisionBits)};
            const RealInterval coshValue = encloseCoshReal(value->asReal(), precisionBits);
            return CertifiedValue{approximation::divide(
                sinhValue, approximation::multiply(value->asReal(), coshValue, precisionBits), precisionBits)};
        }

        const ComplexInterval z = value->asComplex();
        const ComplexInterval sinhValue = encloseComplexSinh(z, precisionBits);
        if (definition->id == BuiltinId::Sinhc)
            return normalizeComplex(approximation::divide(sinhValue, z, precisionBits));
        const ComplexInterval coshValue = encloseComplexCosh(z, precisionBits);
        if (coshValue.containsZero())
            throw PrecisionInsufficient{"Complex tanhc pole cannot yet be excluded"};
        return normalizeComplex(approximation::divide(
            sinhValue, approximation::multiply(z, coshValue, precisionBits), precisionBits));
    }

    case BuiltinId::Gamma:
    case BuiltinId::LogGamma:
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        if (!value->isReal())
            return std::nullopt; // complex special-function backend is intentionally not implemented yet.
        switch (definition->id) {
        case BuiltinId::Gamma:
            return CertifiedValue{encloseGammaReal(value->asReal(), precisionBits)};
        case BuiltinId::LogGamma:
            return CertifiedValue{encloseLogGammaReal(value->asReal(), precisionBits)};
        case BuiltinId::Erf:
            return CertifiedValue{encloseErfReal(value->asReal(), precisionBits)};
        case BuiltinId::Erfc:
            return CertifiedValue{encloseErfcReal(value->asReal(), precisionBits)};
        case BuiltinId::FresnelC:
            return CertifiedValue{encloseFresnelCReal(value->asReal(), precisionBits)};
        case BuiltinId::FresnelS:
            return CertifiedValue{encloseFresnelSReal(value->asReal(), precisionBits)};
        default:
            return std::nullopt;
        }
    }

    case BuiltinId::Hypergeometric1F1: {
        if (call.arguments.size() != 3)
            return std::nullopt;
        const auto a = exactRealRational(call.arguments[0]);
        const auto b = exactRealRational(call.arguments[1]);
        const auto z = exactRealRational(call.arguments[2]);
        if (!a || !b || !z)
            return std::nullopt; // 現backendはexact Rationalのparameter/pointだけを保証評価する。
        return CertifiedValue{encloseHypergeometric1F1Real(
            *a, *b, *z, precisionBits)};
    }

    case BuiltinId::Beta:
    case BuiltinId::BetaLog: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        const auto a = encloseArgument(0);
        const auto b = encloseArgument(1);
        if (!a || !b || !a->isReal() || !b->isReal())
            return std::nullopt;
        if (definition->id == BuiltinId::Beta)
            return CertifiedValue{encloseBetaPositive(a->asReal(), b->asReal(), precisionBits)};
        return CertifiedValue{encloseBetaLogPositive(a->asReal(), b->asReal(), precisionBits)};
    }

    case BuiltinId::Log2:
    case BuiltinId::Log10:
    case BuiltinId::GeneralizedBinomial:
    case BuiltinId::FallingFactorial:
    case BuiltinId::RisingFactorial:
    case BuiltinId::RandSeed:
    case BuiltinId::Rand:
    case BuiltinId::RandInt:
    case BuiltinId::Choice:
    case BuiltinId::RandN:
        // 通常Evaluatorで既存primitiveへexact rewriteされる。ここへ残るものは現段階ではcertified backendを持たない。
        return std::nullopt;

    case BuiltinId::Sqrt: {
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;

        if (value->isReal())
            return principalSqrtReal(value->asReal(), precisionBits);

        // 一般複素数でもprincipal branchをcertifiedに評価する。
        // branch cut（負実軸）を入力区間が跨ぐ場合はcertified_complex_sqrt側が虚部を両符号へ広げ、真値を落とさない。
        return normalizeComplex(enclosePrincipalComplexSqrt(
            value->asComplex(), precisionBits));
    }

    case BuiltinId::Abs: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        if (value->isReal())
            return CertifiedValue{absoluteInterval(value->asReal(), precisionBits)};

        const ComplexInterval complex = value->asComplex();
        const RealInterval magnitudeSquared = approximation::add(
            squareInterval(complex.real(), precisionBits),
            squareInterval(complex.imaginary(), precisionBits),
            precisionBits);
        return CertifiedValue{encloseSqrt(magnitudeSquared, precisionBits).interval};
    }

    case BuiltinId::Sign: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        const RealInterval zero = RealInterval::fromRational(rational(0), precisionBits);
        const RealInterval one = RealInterval::fromRational(rational(1), precisionBits);
        const RealInterval minusOne = RealInterval::fromRational(rational(-1), precisionBits);
        if (value->isReal()) {
            const RealInterval& real = value->asReal();
            if (real.lower() > zero.upper())
                return CertifiedValue{one};
            if (real.upper() < zero.lower())
                return CertifiedValue{minusOne};
            if (real.isPoint() && real.lower().isZero())
                return CertifiedValue{zero};
            throw PrecisionInsufficient{"Sign cannot yet be resolved around zero"};
        }

        const ComplexInterval complex = value->asComplex();
        if (complex.real().isPoint() && complex.real().lower().isZero()
            && complex.imaginary().isPoint() && complex.imaginary().lower().isZero())
            return CertifiedValue{zero};
        const RealInterval magnitudeSquared = approximation::add(
            squareInterval(complex.real(), precisionBits),
            squareInterval(complex.imaginary(), precisionBits),
            precisionBits);
        const RealInterval magnitude = encloseSqrt(magnitudeSquared, precisionBits).interval;
        if (magnitude.containsZero())
            throw PrecisionInsufficient{"Complex sign magnitude cannot yet be proven nonzero"};
        return normalizeComplex(approximation::divide(
            complex, ComplexInterval::fromReal(magnitude), precisionBits));
    }

    case BuiltinId::Re: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        return value->isReal()
            ? *value
            : CertifiedValue{value->asComplex().real()};
    }

    case BuiltinId::Im: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        if (value->isReal())
            return CertifiedValue{RealInterval::fromRational(rational(0), precisionBits)};
        return CertifiedValue{value->asComplex().imaginary()};
    }

    case BuiltinId::Conj: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        if (value->isReal())
            return *value;
        return normalizeComplex(ComplexInterval{
            value->asComplex().real(),
            approximation::negate(value->asComplex().imaginary())});
    }

    case BuiltinId::Power: {
        if (call.arguments.size() != 2)
            return std::nullopt;

        const auto base = encloseArgument(0);
        if (!base)
            return std::nullopt;

        // x^(1/2) は principal sqrt と同義。0.5もexact Rational 1/2なのでここへ来る。
        if (isOneHalf(call.arguments[1])) {
            if (base->isReal())
                return principalSqrtReal(base->asReal(), precisionBits);
            return normalizeComplex(enclosePrincipalComplexSqrt(
                base->asComplex(), precisionBits));
        }

        const auto exponent = exactIntegerExponent(call.arguments[1]);
        if (exponent) {
            if (exponent->isZero())
                return CertifiedValue{RealInterval::fromRational(rational(1), precisionBits)};

            const bool negative = exponent->isNegative();
            const auto magnitude = toUint64(exponent->abs());
            if (!magnitude)
                return std::nullopt;

            CertifiedValue result = integerPower(*base, *magnitude, precisionBits);
            if (!negative)
                return result;

            const CertifiedValue one{RealInterval::fromRational(rational(1), precisionBits)};
            return divideValues(one, result, precisionBits);
        }

        // 一般の非整数・複素指数は principal Power の定義へ一本化する。
        // Power(z,w) := Exp(w * principal Log(z))。
        // exact evaluatorでは式を不用意に展開せずPower表記を保つが、Nの内部意味論はこの定義を使用するため、Solverの「全ての根」とは明確に別物になる。
        const auto exponentValue = encloseArgument(1);
        if (!exponentValue)
            return std::nullopt;
        return normalizeComplex(enclosePrincipalPower(
            base->toComplex(), exponentValue->toComplex(), precisionBits).interval);
    }

    case BuiltinId::Exp: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;

        if (value->isReal())
            return CertifiedValue{encloseExp(value->asReal(), precisionBits).interval};

        return normalizeComplex(encloseComplexExp(
            value->asComplex(), precisionBits).interval);
    }

    case BuiltinId::Log: {
        if (call.arguments.size() < 1 || call.arguments.size() > 2)
            return std::nullopt;

        const auto principalLog = [&](const CertifiedValue& value) -> CertifiedValue {
            // 正実数は専用の単調real Logが最も鋭い。負実数や一般複素数は principal Log = ln|z| + I Arg(z) の共通backendへ送る。
            if (value.isReal()) {
                const numeric::BigFloat zero;
                if (value.asReal().lower() > zero)
                    return CertifiedValue{
                        encloseLogPositive(value.asReal(), precisionBits).interval};
            }
            return normalizeComplex(enclosePrincipalComplexLog(
                value.toComplex(), precisionBits).interval);
        };

        if (call.arguments.size() == 1) {
            const auto value = encloseArgument(0);
            if (!value)
                return std::nullopt;
            return principalLog(*value);
        }

        const auto base = encloseArgument(0);
        const auto value = encloseArgument(1);
        if (!base || !value)
            return std::nullopt;

        const CertifiedValue baseLog = principalLog(*base);
        const CertifiedValue valueLog = principalLog(*value);
        // base=1ならLog[base]=0で数学的に未定義。exactな1はSimplifierがDomainErrorにし、区間が0を含むだけならworking precisionを上げて再試行する。
        const CertifiedValue reciprocalBaseLog = reciprocalValue(
            baseLog, precisionBits, "Logarithm base could not be certified away from one");
        return multiplyValues(valueLog, reciprocalBaseLog, precisionBits);
    }

    case BuiltinId::Arg: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        return CertifiedValue{enclosePrincipalArgument(
            value->toComplex(), precisionBits).interval};
    }

    case BuiltinId::Sin:
    case BuiltinId::Cos:
    case BuiltinId::Tan:
    case BuiltinId::Cot:
    case BuiltinId::Sec:
    case BuiltinId::Csc: {
        if (call.arguments.size() != 1)
            return std::nullopt;

        const Expr& argument = call.arguments.front();
        const auto evaluateRealTrig = [&](const RealInterval& radians) -> CertifiedValue {
            /*
            旧実装ではsinだけを求める場合でもsin/cosを両方certified評価していた。
            巨大radian reductionではPi評価まで二重になるため、必要な函数だけ遅延生成する。
            tan/cotだけは分子・分母の双方が必要なので従来どおり両方を使う。
            */
            switch (definition->id) {
            case BuiltinId::Sin:
                return CertifiedValue{
                    encloseSinRadianInterval(radians, precisionBits).interval};
            case BuiltinId::Cos:
                return CertifiedValue{
                    encloseCosRadianInterval(radians, precisionBits).interval};
            case BuiltinId::Tan: {
                const CertifiedTrigEnclosure sine = encloseSinRadianInterval(radians, precisionBits);
                const CertifiedTrigEnclosure cosine = encloseCosRadianInterval(radians, precisionBits);
                if (cosine.interval.containsZero())
                    throw PrecisionInsufficient{"Tangent denominator cannot yet be proven nonzero"};
                return CertifiedValue{approximation::divide(
                    sine.interval, cosine.interval, precisionBits)};
            }
            case BuiltinId::Cot: {
                const CertifiedTrigEnclosure sine = encloseSinRadianInterval(radians, precisionBits);
                const CertifiedTrigEnclosure cosine = encloseCosRadianInterval(radians, precisionBits);
                if (sine.interval.containsZero())
                    throw PrecisionInsufficient{"Cotangent denominator cannot yet be proven nonzero"};
                return CertifiedValue{approximation::divide(
                    cosine.interval, sine.interval, precisionBits)};
            }
            case BuiltinId::Sec: {
                const CertifiedTrigEnclosure cosine = encloseCosRadianInterval(radians, precisionBits);
                return reciprocalValue(
                    CertifiedValue{cosine.interval}, precisionBits,
                    "Secant denominator cannot yet be proven nonzero");
            }
            case BuiltinId::Csc: {
                const CertifiedTrigEnclosure sine = encloseSinRadianInterval(radians, precisionBits);
                return reciprocalValue(
                    CertifiedValue{sine.interval}, precisionBits,
                    "Cosecant denominator cannot yet be proven nonzero");
            }
            default:
                throw std::logic_error("Unexpected trigonometric builtin");
            }
        };

        // exact turn化できる実角は、巨大角でもPi近似前に周期縮約する。
        if (const auto exactAngle = mathematics::extractExactAngle(
                argument, builtins_, mathematics_, angleSemantics_)) {
            const CertifiedTrigEnclosure sine = encloseSinTurns(exactAngle->turns, precisionBits);
            const CertifiedTrigEnclosure cosine = encloseCosTurns(exactAngle->turns, precisionBits);
            switch (definition->id) {
            case BuiltinId::Sin:
                return CertifiedValue{sine.interval};
            case BuiltinId::Cos:
                return CertifiedValue{cosine.interval};
            case BuiltinId::Tan:
                return CertifiedValue{encloseTanTurns(exactAngle->turns, precisionBits).interval};
            case BuiltinId::Cot:
                if (sine.interval.containsZero())
                    throw PrecisionInsufficient{"Cotangent denominator cannot yet be proven nonzero"};
                return CertifiedValue{approximation::divide(
                    cosine.interval, sine.interval, precisionBits)};
            case BuiltinId::Sec:
                return reciprocalValue(
                    CertifiedValue{cosine.interval}, precisionBits,
                    "Secant denominator cannot yet be proven nonzero");
            case BuiltinId::Csc:
                return reciprocalValue(
                    CertifiedValue{sine.interval}, precisionBits,
                    "Cosecant denominator cannot yet be proven nonzero");
            default:
                break;
            }
        }

        const auto angleOperand = splitAngleOperand(argument, builtins_, angleSemantics_);
        if (!angleOperand)
            return std::nullopt;
        const auto scalar = encloseBound(angleOperand->value, precisionBits, bindings);
        if (!scalar)
            return std::nullopt;

        const CertifiedValue radians = angleValueToRadians(
            *scalar, angleOperand->unit, mathematics_, precisionBits);
        if (radians.isReal())
            return evaluateRealTrig(radians.asReal());

        const ComplexInterval sine = encloseComplexSinRadian(
            radians.asComplex(), precisionBits);
        const ComplexInterval cosine = encloseComplexCosRadian(
            radians.asComplex(), precisionBits);
        switch (definition->id) {
        case BuiltinId::Sin:
            return normalizeComplex(sine);
        case BuiltinId::Cos:
            return normalizeComplex(cosine);
        case BuiltinId::Tan:
            if (cosine.containsZero())
                throw PrecisionInsufficient{"Complex tangent denominator cannot yet be proven nonzero"};
            return normalizeComplex(approximation::divide(sine, cosine, precisionBits));
        case BuiltinId::Cot:
            if (sine.containsZero())
                throw PrecisionInsufficient{"Complex cotangent denominator cannot yet be proven nonzero"};
            return normalizeComplex(approximation::divide(cosine, sine, precisionBits));
        case BuiltinId::Sec:
            return reciprocalValue(normalizeComplex(cosine), precisionBits,
                "Complex secant denominator cannot yet be proven nonzero");
        case BuiltinId::Csc:
            return reciprocalValue(normalizeComplex(sine), precisionBits,
                "Complex cosecant denominator cannot yet be proven nonzero");
        default:
            return std::nullopt;
        }
    }

    case BuiltinId::Asin:
    case BuiltinId::Acos:
    case BuiltinId::Atan: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto input = encloseArgument(0);
        if (!input)
            return std::nullopt;

        CertifiedValue radians = [&]() -> CertifiedValue {
            if (input->isReal()) {
                const RealInterval& real = input->asReal();
                const Rational lower = real.lower().toRational();
                const Rational upper = real.upper().toRational();
                if (definition->id == BuiltinId::Atan)
                    return CertifiedValue{encloseAtan(real, precisionBits).interval};

                const bool insideUnitInterval = lower >= rational(-1) && upper <= rational(1);
                const bool outsideUnitInterval = upper < rational(-1) || lower > rational(1);
                if (insideUnitInterval) {
                    return definition->id == BuiltinId::Asin
                        ? CertifiedValue{encloseAsinRealRadian(real, precisionBits)}
                        : CertifiedValue{encloseAcosRealRadian(real, precisionBits)};
                }
                if (!outsideUnitInterval)
                    throw PrecisionInsufficient{
                        "Inverse trigonometric branch cannot yet be resolved near +/-1"};
            }

            const ComplexInterval complex = input->toComplex();
            if (definition->id == BuiltinId::Asin)
                return normalizeComplex(enclosePrincipalComplexAsinRadian(complex, precisionBits));
            if (definition->id == BuiltinId::Acos)
                return normalizeComplex(enclosePrincipalComplexAcosRadian(complex, precisionBits));
            return normalizeComplex(enclosePrincipalComplexAtanRadian(complex, precisionBits));
        }();

        return radiansToAngleValue(
            radians, angleSemantics_.defaultUnit(), mathematics_, precisionBits);
    }

    case BuiltinId::Atan2: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        const auto y = encloseArgument(0);
        const auto x = encloseArgument(1);
        if (!y || !x)
            return std::nullopt;
        if (!y->isReal() || !x->isReal())
            throw std::domain_error("atan2 expects real arguments");

        const RealInterval radians = enclosePrincipalArgument(
            ComplexInterval{x->asReal(), y->asReal()}, precisionBits).interval;
        return radiansToAngleValue(
            CertifiedValue{radians}, angleSemantics_.defaultUnit(), mathematics_, precisionBits);
    }

    case BuiltinId::Sinh:
    case BuiltinId::Cosh:
    case BuiltinId::Tanh:
    case BuiltinId::Csch:
    case BuiltinId::Sech:
    case BuiltinId::Coth: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto input = encloseArgument(0);
        if (!input)
            return std::nullopt;

        if (input->isReal()) {
            const RealInterval sinhValue = encloseSinhReal(input->asReal(), precisionBits);
            const RealInterval coshValue = encloseCoshReal(input->asReal(), precisionBits);
            switch (definition->id) {
            case BuiltinId::Sinh:
                return CertifiedValue{sinhValue};
            case BuiltinId::Cosh:
                return CertifiedValue{coshValue};
            case BuiltinId::Tanh:
                return CertifiedValue{approximation::divide(
                    sinhValue, coshValue, precisionBits)};
            case BuiltinId::Csch:
                return reciprocalValue(CertifiedValue{sinhValue}, precisionBits,
                    "Hyperbolic cosecant denominator cannot yet be proven nonzero");
            case BuiltinId::Sech:
                return reciprocalValue(CertifiedValue{coshValue}, precisionBits,
                    "Hyperbolic secant denominator cannot yet be proven nonzero");
            case BuiltinId::Coth:
                if (sinhValue.containsZero())
                    throw PrecisionInsufficient{
                        "Hyperbolic cotangent denominator cannot yet be proven nonzero"};
                return CertifiedValue{approximation::divide(
                    coshValue, sinhValue, precisionBits)};
            default:
                break;
            }
        }

        const ComplexInterval complex = input->toComplex();
        const ComplexInterval sinhValue = encloseComplexSinh(complex, precisionBits);
        const ComplexInterval coshValue = encloseComplexCosh(complex, precisionBits);
        switch (definition->id) {
        case BuiltinId::Sinh:
            return normalizeComplex(sinhValue);
        case BuiltinId::Cosh:
            return normalizeComplex(coshValue);
        case BuiltinId::Tanh:
            if (coshValue.containsZero())
                throw PrecisionInsufficient{"Complex tanh denominator cannot yet be proven nonzero"};
            return normalizeComplex(approximation::divide(sinhValue, coshValue, precisionBits));
        case BuiltinId::Csch:
            return reciprocalValue(normalizeComplex(sinhValue), precisionBits,
                "Complex csch denominator cannot yet be proven nonzero");
        case BuiltinId::Sech:
            return reciprocalValue(normalizeComplex(coshValue), precisionBits,
                "Complex sech denominator cannot yet be proven nonzero");
        case BuiltinId::Coth:
            if (sinhValue.containsZero())
                throw PrecisionInsufficient{"Complex coth denominator cannot yet be proven nonzero"};
            return normalizeComplex(approximation::divide(coshValue, sinhValue, precisionBits));
        default:
            return std::nullopt;
        }
    }

    case BuiltinId::Asinh:
    case BuiltinId::Acosh:
    case BuiltinId::Atanh: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto input = encloseArgument(0);
        if (!input)
            return std::nullopt;

        if (input->isReal()) {
            const RealInterval& real = input->asReal();
            const Rational lower = real.lower().toRational();
            const Rational upper = real.upper().toRational();
            if (definition->id == BuiltinId::Asinh)
                return CertifiedValue{encloseAsinhReal(real, precisionBits)};
            if (definition->id == BuiltinId::Acosh) {
                if (lower >= rational(1))
                    return CertifiedValue{encloseAcoshReal(real, precisionBits)};
                if (upper >= rational(1))
                    throw PrecisionInsufficient{
                        "Acosh branch cannot yet be resolved near one"};
            }
            if (definition->id == BuiltinId::Atanh) {
                if (lower > rational(-1) && upper < rational(1))
                    return CertifiedValue{encloseAtanhReal(real, precisionBits)};
                const bool outside = upper < rational(-1) || lower > rational(1);
                if (!outside)
                    throw PrecisionInsufficient{
                        "Atanh branch cannot yet be resolved near +/-1"};
            }
        }

        const ComplexInterval complex = input->toComplex();
        if (definition->id == BuiltinId::Asinh)
            return normalizeComplex(enclosePrincipalComplexAsinh(complex, precisionBits));
        if (definition->id == BuiltinId::Acosh)
            return normalizeComplex(enclosePrincipalComplexAcosh(complex, precisionBits));
        return normalizeComplex(enclosePrincipalComplexAtanh(complex, precisionBits));
    }

    case BuiltinId::Polar:
    case BuiltinId::NextPow2:
    case BuiltinId::Sum:
    case BuiltinId::Product:
    case BuiltinId::Min:
    case BuiltinId::Max:
    case BuiltinId::Mean:
    case BuiltinId::Median:
    case BuiltinId::Mode:
    case BuiltinId::Quantile:
    case BuiltinId::Percentile:
    case BuiltinId::VariancePopulation:
    case BuiltinId::VarianceSample:
    case BuiltinId::StddevPopulation:
    case BuiltinId::StddevSample:
    case BuiltinId::GeometricMean:
    case BuiltinId::HarmonicMean:
    case BuiltinId::Rms:
    case BuiltinId::MedianAbsoluteDeviation:
    case BuiltinId::MeanAbsoluteDeviation:
    case BuiltinId::Skewness:
    case BuiltinId::KurtosisPopulation:
    case BuiltinId::KurtosisSample:
    case BuiltinId::CoefficientVariation:
    case BuiltinId::StandardError:
    case BuiltinId::ZScore:
    case BuiltinId::Iqr:
    case BuiltinId::TrimMean:
    case BuiltinId::WinsorMean:
    case BuiltinId::Winsorized:
    case BuiltinId::Covariance:
    case BuiltinId::Correlation:
    case BuiltinId::SpearmanCorrelation:
    case BuiltinId::PercentRank:
    case BuiltinId::Identity:
    case BuiltinId::Zeros:
    case BuiltinId::MatrixGet:
    case BuiltinId::Trace:
    case BuiltinId::Rows:
    case BuiltinId::Cols:
    case BuiltinId::Diag:
    case BuiltinId::VectorAdd:
    case BuiltinId::VectorSubtract:
    case BuiltinId::VectorScale:
    case BuiltinId::VectorDot:
    case BuiltinId::VectorCross:
    case BuiltinId::VectorNorm:
    case BuiltinId::VectorManhattan:
    case BuiltinId::VectorEuclidean:
    case BuiltinId::VectorNormalize:
    case BuiltinId::VectorProject:
    case BuiltinId::VectorAngle:
    case BuiltinId::VectorReflect:
    case BuiltinId::VectorReflectAxis:
    case BuiltinId::VectorSum:
    case BuiltinId::UnitApplied:
    case BuiltinId::Factorial:
    case BuiltinId::Derivative:
    case BuiltinId::SymbolicIntegral:
    case BuiltinId::Limit:
    case BuiltinId::Floor:
    case BuiltinId::Ceil:
    case BuiltinId::Trunc:
    case BuiltinId::Round:
    case BuiltinId::Frac:
    case BuiltinId::Gcd:
    case BuiltinId::Lcm:
    case BuiltinId::Mod:
    case BuiltinId::Rem:
    case BuiltinId::Quotient:
    case BuiltinId::Permutation:
    case BuiltinId::Combination:
    case BuiltinId::Fibonacci:
    case BuiltinId::DiscreteFourierTransform:
    case BuiltinId::FastFourierTransform:
    case BuiltinId::InverseFourierTransform:
    case BuiltinId::Convolution:
    case BuiltinId::Transpose:
    case BuiltinId::MatrixAdd:
    case BuiltinId::MatrixMultiply:
    case BuiltinId::Determinant:
    case BuiltinId::Inverse:
    case BuiltinId::Rref:
    case BuiltinId::Rank:
    case BuiltinId::NumericDerivative:
    case BuiltinId::NumericIntegral:
    case BuiltinId::NumericalApproximation:
    case BuiltinId::Precision:
    case BuiltinId::Accuracy:
    case BuiltinId::Rationalize:
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify:
    case BuiltinId::Expand:
    case BuiltinId::Factor:
    case BuiltinId::Collect:
    case BuiltinId::Solve:
    case BuiltinId::Set:
    case BuiltinId::SetDelayed:
    case BuiltinId::Less:
    case BuiltinId::LessEqual:
    case BuiltinId::Greater:
    case BuiltinId::GreaterEqual:
    case BuiltinId::Equal:
    case BuiltinId::NotEqual:
    case BuiltinId::LogicalAnd:
    case BuiltinId::Element:
    case BuiltinId::If:
    case BuiltinId::History:
    case BuiltinId::InputHistory:
    case BuiltinId::OutputHistory:
    case BuiltinId::Exit:
    case BuiltinId::Clear:
    case BuiltinId::Definitions:
    case BuiltinId::Undefine:
    case BuiltinId::AngleMode:
        return std::nullopt;
    }

    return std::nullopt;
}

} // namespace mmcal::approximation
