// 式全体の保証付き区間評価
#include "certified_evaluator.hpp"
#include "expression_interval.hpp"

#include "certified_constants.hpp"
#include "certified_domain.hpp"
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
#include "certified_value_math.hpp"
#include "interval_math.hpp"
#include "mathematics/exact_trigonometry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational_rounding.hpp"
#include "numeric/detail/binary_scale.hpp"
#include "symbolic/algebraic_number.hpp"

#include <charconv>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
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

// CertifiedEvaluatorは式木を再帰評価する。Windowsの既定stackでもOS例外へ到達しないよう、
// 構文上の深さとは別に実評価の再帰budgetを持つ。
constexpr std::size_t maximumCertifiedExpressionDepth = 96;
constexpr std::size_t maximumCertifiedRecursionDepth = 48;

[[nodiscard]] Rational absRational(const Rational& value) {
    return value.numerator().isNegative() ? -value : value;
}

[[nodiscard]] Rational intervalMagnitudeUpper(const RealInterval& value) {
    const Rational lower = absRational(value.lower().toRational());
    const Rational upper = absRational(value.upper().toRational());
    return lower < upper ? upper : lower;
}

[[nodiscard]] Rational valueComponentMagnitudeUpper(const CertifiedValue& value) {
    if (value.isReal())
        return intervalMagnitudeUpper(value.asReal());
    const Rational real = intervalMagnitudeUpper(value.asComplex().real());
    const Rational imaginary = intervalMagnitudeUpper(value.asComplex().imaginary());
    return real < imaginary ? imaginary : real;
}

[[nodiscard]] std::size_t cancellationGuardBits(const CertifiedValue& value) {
    const Rational magnitude = valueComponentMagnitudeUpper(value);
    if (magnitude.isZero() || magnitude >= Rational{BigInt{1}})
        return 0;

    const std::int64_t exponent = numeric::detail::floorLog2PositiveRatio(
        magnitude.numerator().abs(), magnitude.denominator());
    if (exponent >= 0)
        return 0;

    const std::uint64_t loss = static_cast<std::uint64_t>(-(exponent + 1)) + 1;
    constexpr std::size_t margin = 24;
    if (loss > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max() - margin))
        throw std::overflow_error("Stable elementary cancellation guard is too large");
    return static_cast<std::size_t>(loss) + margin;
}

[[nodiscard]] std::size_t checkedPrecisionWithGuard(
    std::size_t precisionBits,
    std::size_t guardBits) {
    if (guardBits > std::numeric_limits<std::size_t>::max() - precisionBits)
        throw std::overflow_error("Stable elementary working precision is too large");
    return precisionBits + guardBits;
}

[[nodiscard]] std::optional<std::size_t> positivePrecisionDigits(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    const BigInt& value = expression.asNumber().asReal().asInteger();
    if (value.isNegative() || value.isZero())
        return std::nullopt;
    const std::string text = value.toString();
    std::size_t result = 0;
    const auto converted = std::from_chars(text.data(), text.data() + text.size(), result);
    if (converted.ec != std::errc{} || converted.ptr != text.data() + text.size())
        return std::nullopt;
    return result;
}

[[nodiscard]] bool exceedsCertifiedExpressionDepth(const Expr& root) {
    struct Pending final {
        Expr expression;
        std::size_t depth = 0;
    };

    std::vector<Pending> pending;
    pending.push_back(Pending{root, 1});
    while (!pending.empty()) {
        Pending current = std::move(pending.back());
        pending.pop_back();
        if (current.depth > maximumCertifiedExpressionDepth)
            return true;

        if (current.expression.isCall()) {
            const auto& arguments = current.expression.asCall().arguments;
            for (const Expr& argument : arguments)
                pending.push_back(Pending{argument, current.depth + 1});
        }
        else if (current.expression.isArray()) {
            for (const Expr& element : current.expression.asArray().storedExpressions())
                pending.push_back(Pending{element, current.depth + 1});
        }
        else if (current.expression.isList()) {
            const auto& elements = current.expression.asList().elements;
            for (const Expr& element : elements)
                pending.push_back(Pending{element, current.depth + 1});
        }
    }
    return false;
}

[[nodiscard]] std::optional<symbolic::RealAlgebraicNumber> algebraicRoot(
    const expression::CallExpr& call) {
    if (call.arguments.size() != 2 || !call.arguments[0].isArray()
        || call.arguments[0].asArray().rank() != 1)
        return std::nullopt;

    const auto& coefficientsArray = call.arguments[0].asArray();
    std::vector<Rational> coefficients;
    coefficients.reserve(coefficientsArray.size());
    for (std::size_t i = 0; i < coefficientsArray.size(); ++i) {
        const Expr value = coefficientsArray.element(i);
        if (!value.isNumber() || !value.asNumber().isReal())
            return std::nullopt;
        coefficients.push_back(value.asNumber().asReal().toRational());
    }

    const Expr& indexExpression = call.arguments[1];
    if (!indexExpression.isNumber() || !indexExpression.asNumber().isReal()
        || !indexExpression.asNumber().asReal().isInteger())
        return std::nullopt;
    const BigInt& indexInteger = indexExpression.asNumber().asReal().asInteger();
    if (indexInteger.isNegative() || indexInteger.isZero())
        return std::nullopt;
    const auto index = numeric::tryToUint64(indexInteger);
    if (!index || *index > std::numeric_limits<std::size_t>::max())
        return std::nullopt;
    return symbolic::RealAlgebraicNumber::create(
        coefficients, static_cast<std::size_t>(*index));
}


[[nodiscard]] std::optional<symbolic::ComplexAlgebraicNumber> complexAlgebraicRoot(
    const expression::CallExpr& call) {
    if (call.arguments.size() != 3 || !call.arguments[2].isSymbol()
        || call.arguments[2].asSymbol().view() != "Complex"
        || !call.arguments[0].isArray() || call.arguments[0].asArray().rank() != 1)
        return std::nullopt;

    const auto& coefficientsArray = call.arguments[0].asArray();
    std::vector<Rational> coefficients;
    coefficients.reserve(coefficientsArray.size());
    for (std::size_t i = 0; i < coefficientsArray.size(); ++i) {
        const Expr value = coefficientsArray.element(i);
        if (!value.isNumber() || !value.asNumber().isReal())
            return std::nullopt;
        coefficients.push_back(value.asNumber().asReal().toRational());
    }

    const Expr& indexExpression = call.arguments[1];
    if (!indexExpression.isNumber() || !indexExpression.asNumber().isReal()
        || !indexExpression.asNumber().asReal().isInteger())
        return std::nullopt;
    const BigInt& indexInteger = indexExpression.asNumber().asReal().asInteger();
    if (indexInteger.isNegative() || indexInteger.isZero())
        return std::nullopt;
    const auto index = numeric::tryToUint64(indexInteger);
    if (!index || *index > std::numeric_limits<std::size_t>::max())
        return std::nullopt;
    return symbolic::ComplexAlgebraicNumber::create(
        coefficients, static_cast<std::size_t>(*index));
}

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

[[nodiscard]] Rational rationalPower(Rational base, std::size_t exponent) {
    Rational result{BigInt{1}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result *= base;
        exponent >>= 1U;
        if (exponent != 0)
            base *= base;
    }
    return result;
}

[[nodiscard]] RealInterval scaleInterval(
    const RealInterval& value,
    const Rational& scale,
    std::size_t precisionBits) {
    return approximation::multiply(
        value, RealInterval::fromRational(scale, precisionBits), precisionBits);
}

[[nodiscard]] RealInterval addSymmetricRemainder(
    const RealInterval& value,
    const Rational& radius,
    std::size_t precisionBits) {
    return approximation::add(value,
        RealInterval::fromRationalBounds(-radius, radius, precisionBits), precisionBits);
}

// 可除特異点を含む小区間では商を直接作らず，Taylor多項式と剰余上界で包む。
// |x|<=1だけを対象とし，それ以外は従来backendへ委ねる。
[[nodiscard]] std::optional<RealInterval> stableCardinalRealNearZero(
    BuiltinId id,
    const RealInterval& x,
    std::size_t precisionBits) {
    const Rational radius = intervalMagnitudeUpper(x);
    if (radius > rational(1))
        return std::nullopt;

    const RealInterval one = RealInterval::fromRational(rational(1), precisionBits);
    const RealInterval x2 = approximation::squareInterval(x, precisionBits);
    const RealInterval x3 = approximation::multiply(x2, x, precisionBits);
    const RealInterval x4 = approximation::squareInterval(x2, precisionBits);
    const RealInterval x5 = approximation::multiply(x4, x, precisionBits);

    switch (id) {
    case BuiltinId::Sinc: {
        RealInterval result = approximation::subtract(
            one, scaleInterval(x2, rational(1, 6), precisionBits), precisionBits);
        result = approximation::add(
            result, scaleInterval(x4, rational(1, 120), precisionBits), precisionBits);
        return addSymmetricRemainder(
            result, rationalPower(radius, 6) * rational(1, 5040), precisionBits);
    }
    case BuiltinId::Cosc: {
        RealInterval result = scaleInterval(x, rational(1, 2), precisionBits);
        result = approximation::subtract(
            result, scaleInterval(x3, rational(1, 24), precisionBits), precisionBits);
        result = approximation::add(
            result, scaleInterval(x5, rational(1, 720), precisionBits), precisionBits);
        return addSymmetricRemainder(
            result, rationalPower(radius, 7) * rational(1, 40320), precisionBits);
    }
    case BuiltinId::Sinhc: {
        RealInterval result = approximation::add(
            one, scaleInterval(x2, rational(1, 6), precisionBits), precisionBits);
        result = approximation::add(
            result, scaleInterval(x4, rational(1, 120), precisionBits), precisionBits);
        // k>=3の正項tailは|x|<=1で初項r^6/7!と比率<=1/72の幾何級数で抑える。
        return addSymmetricRemainder(
            result, rationalPower(radius, 6) * rational(1, 4970), precisionBits);
    }
    case BuiltinId::Expc: {
        RealInterval result = approximation::add(
            one, scaleInterval(x, rational(1, 2), precisionBits), precisionBits);
        result = approximation::add(
            result, scaleInterval(x2, rational(1, 6), precisionBits), precisionBits);
        result = approximation::add(
            result, scaleInterval(x3, rational(1, 24), precisionBits), precisionBits);
        result = approximation::add(
            result, scaleInterval(x4, rational(1, 120), precisionBits), precisionBits);
        // expのLagrange剰余でe^|x|<3を使う。
        return addSymmetricRemainder(
            result, rationalPower(radius, 5) * rational(1, 240), precisionBits);
    }
    default:
        return std::nullopt;
    }
}

[[nodiscard]] RealInterval exactRealInterval(
    const RealNumber& value,
    std::size_t precisionBits) {
    return RealInterval::fromRational(value.toRational(), precisionBits);
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
            result = multiplyCertifiedValues(result, base, precisionBits);

        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyCertifiedValues(base, base, precisionBits);
    }

    return result;
}

[[nodiscard]] std::optional<RealInterval> provablyReal(const CertifiedValue& value) {
    if (value.isReal())
        return value.asReal();
    if (approximation::intervalIsExactZero(value.asComplex().imaginary()))
        return value.asComplex().real();
    return std::nullopt;
}

[[nodiscard]] bool informationStraddlesRealOnlyBackend(const CertifiedValue& value) {
    return value.isComplex()
        && value.asComplex().imaginary().containsZero()
        && !approximation::intervalIsExactZero(value.asComplex().imaginary());
}

void requireNoNonPositiveIntegerPole(
    const CertifiedValue& value,
    std::string_view functionName,
    bool informationEnclosure = false) {
    const SingularityRelation relation = value.isReal()
        ? classifyNonPositiveIntegerPole(value.asReal())
        : classifyNonPositiveIntegerPole(value.asComplex());
    switch (relation) {
    case SingularityRelation::Clear:
        return;
    case SingularityRelation::ExactSingularity:
        throw std::domain_error(
            std::string{functionName} + " is undefined at a non-positive integer");
    case SingularityRelation::MayContainSingularity:
        throw PrecisionInsufficient(
            std::string{functionName}
                + " InformationEnclosure may contain a non-positive-integer pole",
            informationEnclosure
                ? PrecisionInsufficientKind::InputInformation
                : PrecisionInsufficientKind::Refinable);
    }
}

void requireNoPointSingularity(
    const CertifiedValue& value,
    const Rational& point,
    std::string_view domainMessage,
    std::string_view precisionMessage,
    bool informationEnclosure = false) {
    const ComplexInterval complex = value.toComplex();
    const Rational realLower = complex.real().lower().toRational();
    const Rational realUpper = complex.real().upper().toRational();
    if (point < realLower || point > realUpper || !complex.imaginary().containsZero())
        return;

    const bool exactPoint = complex.real().isPoint()
        && realLower == point
        && intervalIsExactZero(complex.imaginary());
    if (exactPoint)
        throw std::domain_error(std::string{domainMessage});

    throw PrecisionInsufficient(
        std::string{precisionMessage},
        informationEnclosure
            ? PrecisionInsufficientKind::InputInformation
            : PrecisionInsufficientKind::Refinable);
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
    return normalizeCertifiedComplex(approximation::multiply(
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
    return normalizeCertifiedComplex(approximation::divide(
        ComplexInterval::fromReal(RealInterval::fromRational(rational(1), precisionBits)),
        value.asComplex(), precisionBits));
}

} // namespace

CertifiedEvaluator::CertifiedEvaluator(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angleSemantics)
    : builtins_(builtins),
      mathematics_(mathematics),
      angleSemantics_(angleSemantics) {}

std::optional<CertifiedValue> CertifiedEvaluator::enclose(
    const Expr& expression,
    std::size_t precisionBits,
    EnclosureKind enclosureKind) const {
    if (exceedsCertifiedExpressionDepth(expression))
        return std::nullopt;
    return encloseBound(expression, precisionBits, {}, enclosureKind, 0);
}

std::optional<CertifiedValue> CertifiedEvaluator::enclose(
    const Expr& expression,
    std::size_t precisionBits,
    std::span<const CertifiedBinding> bindings,
    EnclosureKind enclosureKind) const {
    if (exceedsCertifiedExpressionDepth(expression))
        return std::nullopt;
    return encloseBound(expression, precisionBits, bindings, enclosureKind, 0);
}

std::optional<CertifiedValue> CertifiedEvaluator::encloseBound(
    const Expr& expression,
    std::size_t precisionBits,
    std::span<const CertifiedBinding> bindings,
    EnclosureKind enclosureKind,
    std::size_t recursionDepth) const {
    if (recursionDepth > maximumCertifiedRecursionDepth)
        return std::nullopt;
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

    if (expression.isDecimalApproximation()) {
        const auto& value = expression.asDecimalApproximation();
        const Rational& lower = enclosureKind == EnclosureKind::Information
            ? value.informationLower() : value.certifiedLower();
        const Rational& upper = enclosureKind == EnclosureKind::Information
            ? value.informationUpper() : value.certifiedUpper();
        return CertifiedValue{RealInterval::fromRationalBounds(lower, upper, precisionBits)};
    }

    if (expression.isComplexDecimalApproximation()) {
        const auto& value = expression.asComplexDecimalApproximation();
        const auto component = [&](const numeric::DecimalApproximation& part,
                                   bool exactlyZero,
                                   const Rational& informationLower,
                                   const Rational& informationUpper) {
            if (enclosureKind == EnclosureKind::Information && exactlyZero)
                return RealInterval::fromRational(Rational{}, precisionBits);
            const Rational& lower = enclosureKind == EnclosureKind::Information
                ? informationLower : part.certifiedLower();
            const Rational& upper = enclosureKind == EnclosureKind::Information
                ? informationUpper : part.certifiedUpper();
            return RealInterval::fromRationalBounds(lower, upper, precisionBits);
        };
        return CertifiedValue{ComplexInterval{
            component(value.real(), value.realExactlyZero(),
                value.realInformationLower(), value.realInformationUpper()),
            component(value.imaginary(), value.imaginaryExactlyZero(),
                value.imaginaryInformationLower(), value.imaginaryInformationUpper())}};
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
        return encloseCall(expression.asCall(), precisionBits, bindings, enclosureKind, recursionDepth);

    return std::nullopt;
}

std::optional<CertifiedValue> CertifiedEvaluator::encloseCall(
    const expression::CallExpr& call,
    std::size_t precisionBits,
    std::span<const CertifiedBinding> bindings,
    EnclosureKind enclosureKind,
    std::size_t recursionDepth) const {
    if (recursionDepth > maximumCertifiedRecursionDepth)
        return std::nullopt;
    const auto* definition = builtins_.find(call.head);
    if (!definition)
        return std::nullopt;

    const auto encloseArgument = [&](std::size_t index) -> std::optional<CertifiedValue> {
        if (index >= call.arguments.size())
            return std::nullopt;
        return encloseBound(call.arguments[index], precisionBits, bindings, enclosureKind, recursionDepth + 1);
    };

    switch (definition->id) {
    case BuiltinId::NumericalApproximation: {
        if (call.arguments.empty() || call.arguments.size() > 2)
            return std::nullopt;
        constexpr std::size_t defaultPrecisionDigits = 16;
        const std::size_t digits = call.arguments.size() == 2
            ? positivePrecisionDigits(call.arguments[1]).value_or(0)
            : defaultPrecisionDigits;
        if (digits == 0)
            return std::nullopt;

        ApproximationContext context{digits};
        for (;;) {
            const std::size_t nestedBits = context.workingBinaryBits();
            const auto information = encloseBound(
                call.arguments[0], nestedBits, bindings, EnclosureKind::Information, recursionDepth + 1);
            const auto certified = encloseBound(
                call.arguments[0], nestedBits, bindings, EnclosureKind::Certified, recursionDepth + 1);
            if (!information || !certified)
                return std::nullopt;
            if (const auto materialized = finalizeCertifiedApproximation(
                    *certified, *information, digits))
                return encloseBound(
                    *materialized, precisionBits, bindings, enclosureKind, recursionDepth + 1);

            const std::size_t growth = std::max<std::size_t>(8, context.guardDigits() / 2);
            if (growth > std::numeric_limits<std::size_t>::max() - context.guardDigits())
                return std::nullopt;
            context.setGuardDigits(context.guardDigits() + growth);
            if (context.guardDigits() > digits + 256)
                return std::nullopt;
        }
    }

    case BuiltinId::Root: {
        if (call.arguments.size() == 3) {
            const auto algebraic = complexAlgebraicRoot(call);
            if (!algebraic)
                return std::nullopt;
            if (const auto general = symbolic::AlgebraicNumber::create(
                    algebraic->polynomial(), algebraic->rootIndex(),
                    symbolic::AlgebraicRootDomain::Complex)) {
                if (const auto exact = general->exactRationalParts())
                    return CertifiedValue{ComplexInterval{
                        RealInterval::fromRational(exact->first, precisionBits),
                        RealInterval::fromRational(exact->second, precisionBits)}};
            }
            const symbolic::RationalComplexDisk disk = algebraic->refined(precisionBits);
            return CertifiedValue{ComplexInterval{
                RealInterval::fromRationalBounds(
                    disk.real - disk.radius, disk.real + disk.radius, precisionBits),
                RealInterval::fromRationalBounds(
                    disk.imaginary - disk.radius, disk.imaginary + disk.radius, precisionBits)}};
        }
        const auto algebraic = algebraicRoot(call);
        if (!algebraic)
            return std::nullopt;
        const symbolic::RationalRootInterval interval = algebraic->refined(precisionBits);
        return CertifiedValue{RealInterval::fromRationalBounds(
            interval.lower, interval.upper, precisionBits)};
    }

    case BuiltinId::Add: {
        CertifiedValue result{RealInterval::fromRational(rational(0), precisionBits)};
        for (const Expr& argument : call.arguments) {
            const auto enclosed = encloseBound(argument, precisionBits, bindings, enclosureKind, recursionDepth + 1);
            if (!enclosed)
                return std::nullopt;
            result = addCertifiedValues(result, *enclosed, precisionBits);
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
        return subtractCertifiedValues(*lhs, *rhs, precisionBits);
    }

    case BuiltinId::Multiply: {
        CertifiedValue result{RealInterval::fromRational(rational(1), precisionBits)};
        for (const Expr& argument : call.arguments) {
            const auto enclosed = encloseBound(argument, precisionBits, bindings, enclosureKind, recursionDepth + 1);
            if (!enclosed)
                return std::nullopt;
            result = multiplyCertifiedValues(result, *enclosed, precisionBits);
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
        return divideCertifiedValues(*lhs, *rhs, precisionBits);
    }

    case BuiltinId::Negate: {
        const auto value = encloseArgument(0);
        return value ? std::optional<CertifiedValue>{negateCertifiedValue(*value)} : std::nullopt;
    }

    case BuiltinId::Frac: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value || !value->isReal())
            return std::nullopt;
        const Rational lower = value->asReal().lower().toRational();
        const Rational upper = value->asReal().upper().toRational();
        const BigInt lowerFloor = numeric::floorToInteger(lower);
        const BigInt upperFloor = numeric::floorToInteger(upper);
        if (lowerFloor != upperFloor)
            throw PrecisionInsufficient(
                "fract InformationEnclosure crosses an integer boundary",
                enclosureKind == EnclosureKind::Information
                    ? PrecisionInsufficientKind::InputInformation
                    : PrecisionInsufficientKind::Refinable);
        return CertifiedValue{RealInterval::fromRationalBounds(
            lower - Rational{lowerFloor}, upper - Rational{lowerFloor}, precisionBits)};
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
        const auto scalar = encloseBound(angleOperand->value, precisionBits, bindings, enclosureKind, recursionDepth + 1);
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
        return normalizeCertifiedComplex(approximation::add(cosine, iSine, precisionBits));
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
        const auto initial = encloseArgument(0);
        if (!initial)
            return std::nullopt;

        // exp(x)-1 は x≈0 で桁落ちする。入力の2進桁位置だけ余分に作業精度を
        // 与え，同じ式を再包含してから減算する。有限precision入力では
        // InformationEnclosure自体は狭まらないため，隠れた精度を発明しない。
        const std::size_t workBits = checkedPrecisionWithGuard(
            precisionBits, cancellationGuardBits(*initial));
        const auto value = workBits == precisionBits
            ? initial
            : encloseBound(call.arguments[0], workBits, bindings, enclosureKind, recursionDepth + 1);
        if (!value)
            return std::nullopt;
        const CertifiedValue one{RealInterval::fromRational(rational(1), workBits)};
        CertifiedValue result = value->isReal()
            ? subtractCertifiedValues(
                CertifiedValue{encloseExp(value->asReal(), workBits).interval}, one, workBits)
            : subtractCertifiedValues(
                normalizeCertifiedComplex(encloseComplexExp(value->asComplex(), workBits).interval),
                one, workBits);
        return result.isReal()
            ? CertifiedValue{result.asReal().roundedOutward(precisionBits)}
            : CertifiedValue{result.asComplex().roundedOutward(precisionBits)};
    }

    case BuiltinId::Log1p: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto initial = encloseArgument(0);
        if (!initial)
            return std::nullopt;

        // log(1+x) は x≈0 の加算とlog評価，x≈-1 の 1+x の双方で
        // cancellationを起こし得る。初回区間から必要guardを見積もり，
        // shifted側がさらに小さい場合はその桁落ち分も加える。
        const CertifiedValue initialOne{
            RealInterval::fromRational(rational(1), precisionBits)};
        const CertifiedValue initialShifted = addCertifiedValues(
            *initial, initialOne, precisionBits);
        const std::size_t guardBits = std::max(
            cancellationGuardBits(*initial), cancellationGuardBits(initialShifted));
        const std::size_t workBits = checkedPrecisionWithGuard(precisionBits, guardBits);
        const auto value = workBits == precisionBits
            ? initial
            : encloseBound(call.arguments[0], workBits, bindings, enclosureKind, recursionDepth + 1);
        if (!value)
            return std::nullopt;
        const CertifiedValue shifted = addCertifiedValues(
            *value, CertifiedValue{RealInterval::fromRational(rational(1), workBits)}, workBits);
        if (shifted.isReal()) {
            const numeric::BigFloat zero;
            if (shifted.asReal().lower() > zero)
                return CertifiedValue{
                    encloseLogPositive(shifted.asReal(), workBits).interval.roundedOutward(precisionBits)};
        }
        return normalizeCertifiedComplex(enclosePrincipalComplexLog(
            shifted.toComplex(), workBits).interval.roundedOutward(precisionBits));
    }

    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto angleOperand = splitAngleOperand(call.arguments[0], builtins_, angleSemantics_);
        if (!angleOperand)
            return std::nullopt;
        const auto scalar = encloseBound(angleOperand->value, precisionBits, bindings, enclosureKind, recursionDepth + 1);
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
            if (x.containsZero()) {
                if (definition->id == BuiltinId::Sinc || definition->id == BuiltinId::Cosc) {
                    if (const auto stable = stableCardinalRealNearZero(
                            definition->id, x, precisionBits))
                        return CertifiedValue{*stable};
                }
                else if (definition->id == BuiltinId::Tanc) {
                    if (const auto stableSinc = stableCardinalRealNearZero(
                            BuiltinId::Sinc, x, precisionBits)) {
                        const RealInterval cosine = encloseCosRadianInterval(
                            x, precisionBits).interval;
                        if (!cosine.containsZero())
                            return CertifiedValue{approximation::divide(
                                *stableSinc, cosine, precisionBits)};
                    }
                }
                throw PrecisionInsufficient{
                    "Cardinal trigonometric input information crosses the removable singularity",
                    enclosureKind == EnclosureKind::Information
                        ? PrecisionInsufficientKind::InputInformation
                        : PrecisionInsufficientKind::Refinable};
            }
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
            return normalizeCertifiedComplex(approximation::divide(sine, x, precisionBits));
        if (definition->id == BuiltinId::Cosc) {
            const ComplexInterval one = ComplexInterval::fromReal(oneReal);
            return normalizeCertifiedComplex(approximation::divide(
                approximation::subtract(one, cosine, precisionBits), x, precisionBits));
        }
        if (cosine.containsZero())
            throw PrecisionInsufficient{"Complex tanc pole cannot yet be excluded"};
        return normalizeCertifiedComplex(approximation::divide(
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
        if (value->isReal() && value->asReal().containsZero()) {
            if (definition->id == BuiltinId::Sinhc || definition->id == BuiltinId::Expc) {
                if (const auto stable = stableCardinalRealNearZero(
                        definition->id, value->asReal(), precisionBits))
                    return CertifiedValue{*stable};
            }
            else if (definition->id == BuiltinId::Tanhc) {
                if (const auto stableSinhc = stableCardinalRealNearZero(
                        BuiltinId::Sinhc, value->asReal(), precisionBits)) {
                    const RealInterval coshValue = encloseCoshReal(
                        value->asReal(), precisionBits);
                    if (!coshValue.containsZero())
                        return CertifiedValue{approximation::divide(
                            *stableSinhc, coshValue, precisionBits)};
                }
            }
            throw PrecisionInsufficient{
                "Cardinal function input information crosses the removable singularity",
                enclosureKind == EnclosureKind::Information
                    ? PrecisionInsufficientKind::InputInformation
                    : PrecisionInsufficientKind::Refinable};
        }
        if (value->isComplex() && value->asComplex().containsZero())
            throw PrecisionInsufficient{
                "Complex cardinal function denominator cannot yet be proven nonzero",
                enclosureKind == EnclosureKind::Information
                    ? PrecisionInsufficientKind::InputInformation
                    : PrecisionInsufficientKind::Refinable};

        if (definition->id == BuiltinId::Expc) {
            const CertifiedValue one{RealInterval::fromRational(rational(1), precisionBits)};
            const CertifiedValue exponential = value->isReal()
                ? CertifiedValue{encloseExp(value->asReal(), precisionBits).interval}
                : normalizeCertifiedComplex(encloseComplexExp(value->asComplex(), precisionBits).interval);
            return divideCertifiedValues(subtractCertifiedValues(exponential, one, precisionBits), *value, precisionBits);
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
            return normalizeCertifiedComplex(approximation::divide(sinhValue, z, precisionBits));
        const ComplexInterval coshValue = encloseComplexCosh(z, precisionBits);
        if (coshValue.containsZero())
            throw PrecisionInsufficient{"Complex tanhc pole cannot yet be excluded"};
        return normalizeCertifiedComplex(approximation::divide(
            sinhValue, approximation::multiply(z, coshValue, precisionBits), precisionBits));
    }

    case BuiltinId::Gamma:
    case BuiltinId::LogGamma:
    case BuiltinId::Zeta:
    case BuiltinId::Digamma:
    case BuiltinId::Trigamma:
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto exactRational = exactRealRational(call.arguments[0]);
        if (exactRational && definition->id == BuiltinId::Gamma)
            return CertifiedValue{encloseGammaRational(*exactRational, precisionBits)};
        if (exactRational && definition->id == BuiltinId::LogGamma)
            return CertifiedValue{encloseLogGammaRational(*exactRational, precisionBits)};

        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;

        if (definition->id == BuiltinId::Gamma
            || definition->id == BuiltinId::LogGamma
            || definition->id == BuiltinId::Digamma
            || definition->id == BuiltinId::Trigamma)
            requireNoNonPositiveIntegerPole(
                *value, definition->name(), enclosureKind == EnclosureKind::Information);
        if (definition->id == BuiltinId::Zeta)
            requireNoPointSingularity(
                *value, rational(1),
                "zeta has a pole at s = 1",
                "zeta InformationEnclosure may contain the pole at s = 1",
                enclosureKind == EnclosureKind::Information);

        if (!value->isReal()) {
            switch (definition->id) {
            case BuiltinId::Gamma:
                return normalizeCertifiedComplex(encloseGammaComplex(value->asComplex(), precisionBits));
            case BuiltinId::Zeta:
                return normalizeCertifiedComplex(encloseZetaComplex(value->asComplex(), precisionBits));
            case BuiltinId::Digamma:
                return normalizeCertifiedComplex(encloseDigammaComplex(value->asComplex(), precisionBits));
            case BuiltinId::Trigamma:
                return normalizeCertifiedComplex(encloseTrigammaComplex(value->asComplex(), precisionBits));
            case BuiltinId::Erf:
                return normalizeCertifiedComplex(encloseErfComplex(value->asComplex(), precisionBits));
            case BuiltinId::Erfc:
                return normalizeCertifiedComplex(encloseErfcComplex(value->asComplex(), precisionBits));
            case BuiltinId::FresnelC:
                return normalizeCertifiedComplex(encloseFresnelCComplex(value->asComplex(), precisionBits));
            case BuiltinId::FresnelS:
                return normalizeCertifiedComplex(encloseFresnelSComplex(value->asComplex(), precisionBits));
            default:
                return std::nullopt;
            }
        }
        switch (definition->id) {
        case BuiltinId::Gamma:
            return CertifiedValue{encloseGammaReal(value->asReal(), precisionBits)};
        case BuiltinId::LogGamma:
            return CertifiedValue{encloseLogGammaReal(value->asReal(), precisionBits)};
        case BuiltinId::Zeta:
            return CertifiedValue{encloseZetaReal(value->asReal(), precisionBits)};
        case BuiltinId::Digamma:
            if (value->asReal().lower().toRational() > Rational{})
                return CertifiedValue{encloseDigammaPositive(value->asReal(), precisionBits)};
            // digammaは非正整数のpoleを除けば負の実軸上でも実数値を持つ。
            // positive-only real backendへ誤送せず，複素recurrenceでpoleを分類して実部へ射影する。
            return CertifiedValue{encloseDigammaComplex(
                ComplexInterval::fromReal(value->asReal()), precisionBits).real()};
        case BuiltinId::Trigamma:
            if (value->asReal().lower().toRational() > Rational{})
                return CertifiedValue{encloseTrigammaPositive(value->asReal(), precisionBits)};
            return CertifiedValue{encloseTrigammaComplex(
                ComplexInterval::fromReal(value->asReal()), precisionBits).real()};
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

    case BuiltinId::LambertW: {
        if (call.arguments.empty() || call.arguments.size() > 2)
            return std::nullopt;

        int branch = 0;
        std::size_t valueIndex = 0;
        if (call.arguments.size() == 2) {
            const auto branchInteger = exactIntegerExponent(call.arguments[0]);
            if (!branchInteger)
                return std::nullopt;
            if (*branchInteger == BigInt{0})
                branch = 0;
            else if (*branchInteger == BigInt{-1})
                branch = -1;
            else
                throw CertifiedBackendUnsupported{
                    "Certified complex Lambert W branches other than 0 and -1 are not implemented"};
            valueIndex = 1;
        }

        const auto value = encloseArgument(valueIndex);
        if (!value)
            return std::nullopt;
        if (!value->isReal())
            throw CertifiedBackendUnsupported{
                "Certified complex Lambert W evaluation is not implemented"};
        return CertifiedValue{encloseLambertWReal(value->asReal(), branch, precisionBits)};
    }

    case BuiltinId::Hypergeometric1F1: {
        if (call.arguments.size() != 3)
            return std::nullopt;
        const auto a = exactRealRational(call.arguments[0]);
        const auto b = exactRealRational(call.arguments[1]);
        const auto z = exactRealRational(call.arguments[2]);
        if (a && b && z)
            return CertifiedValue{encloseHypergeometric1F1Real(
                *a, *b, *z, precisionBits)};

        const auto complexA = encloseArgument(0);
        const auto complexB = encloseArgument(1);
        const auto complexZ = encloseArgument(2);
        if (!complexA || !complexB || !complexZ)
            return std::nullopt;
        requireNoNonPositiveIntegerPole(
            *complexB, "hypergeometric1F1 denominator parameter",
            enclosureKind == EnclosureKind::Information);
        const ComplexInterval result = encloseHypergeometric1F1Complex(
            complexA->toComplex(), complexB->toComplex(),
            complexZ->toComplex(), precisionBits);
        // 実parameter/実argumentでは1F1は実数値。complex interval backendの
        // 虚部に丸め由来のzero-centered幅が残っても，実数値であること自体は
        // 数学的に保証されるためreal projectionしてbackend経路差を露出させない。
        if (complexA->isReal() && complexB->isReal() && complexZ->isReal())
            return CertifiedValue{result.real()};
        return normalizeCertifiedComplex(result);
    }

    case BuiltinId::Hypergeometric2F1: {
        if (call.arguments.size() != 4)
            return std::nullopt;
        const auto a = exactRealRational(call.arguments[0]);
        const auto b = exactRealRational(call.arguments[1]);
        const auto c = exactRealRational(call.arguments[2]);
        const auto z = exactRealRational(call.arguments[3]);
        if (a && b && c && z) {
            const Rational absZ = z->numerator().isNegative() ? -*z : *z;
            if (absZ <= rational(9, 10))
                return CertifiedValue{encloseHypergeometric2F1Real(
                    *a, *b, *c, *z, precisionBits)};
        }

        const auto complexA = encloseArgument(0);
        const auto complexB = encloseArgument(1);
        const auto complexC = encloseArgument(2);
        const auto complexZ = encloseArgument(3);
        if (!complexA || !complexB || !complexC || !complexZ)
            return std::nullopt;
        requireNoNonPositiveIntegerPole(
            *complexC, "hypergeometric2F1 denominator parameter",
            enclosureKind == EnclosureKind::Information);
        if (enclosureKind == EnclosureKind::Information
            && informationIsAmbiguousAtPositiveRealCut(complexZ->toComplex()))
            throw PrecisionInsufficient{
                "hypergeometric2F1 InformationEnclosure cannot determine the principal branch-cut side",
                PrecisionInsufficientKind::InputInformation};
        const ComplexInterval result = encloseHypergeometric2F1Complex(
            complexA->toComplex(), complexB->toComplex(),
            complexC->toComplex(), complexZ->toComplex(), precisionBits);
        // principal branch cutはz>=1。入力区間全体が実軸上かつ1未満なら
        // Gauss continuationも実数値なので，complex backendの虚部roundoffを捨てられる。
        if (complexA->isReal() && complexB->isReal() && complexC->isReal()
            && complexZ->isReal()
            && complexZ->asReal().upper().toRational() < rational(1))
            return CertifiedValue{result.real()};
        return normalizeCertifiedComplex(result);
    }

    case BuiltinId::EllipticF:
    case BuiltinId::EllipticE: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        const auto exactPhi = exactRealRational(call.arguments[0]);
        const auto exactM = exactRealRational(call.arguments[1]);
        if (exactPhi && exactM)
            return CertifiedValue{definition->id == BuiltinId::EllipticF
                ? encloseEllipticFReal(*exactPhi, *exactM, precisionBits)
                : encloseEllipticEReal(*exactPhi, *exactM, precisionBits)};

        const auto phi = encloseArgument(0);
        const auto m = encloseArgument(1);
        if (!phi || !m)
            return std::nullopt;
        const auto realPhi = provablyReal(*phi);
        const auto realM = provablyReal(*m);
        if (!realPhi || !realM) {
            if (enclosureKind == EnclosureKind::Information
                && (informationStraddlesRealOnlyBackend(*phi)
                    || informationStraddlesRealOnlyBackend(*m)))
                throw PrecisionInsufficient{
                    "Elliptic integral InformationEnclosure straddles the real-only certified backend boundary",
                    PrecisionInsufficientKind::InputInformation};
            throw CertifiedBackendUnsupported{
                "Certified complex elliptic F/E backend is not implemented"};
        }
        return CertifiedValue{definition->id == BuiltinId::EllipticF
            ? encloseEllipticFReal(*realPhi, *realM, precisionBits)
            : encloseEllipticEReal(*realPhi, *realM, precisionBits)};
    }

    case BuiltinId::EllipticPi: {
        if (call.arguments.size() != 3)
            return std::nullopt;
        const auto exactN = exactRealRational(call.arguments[0]);
        const auto exactPhi = exactRealRational(call.arguments[1]);
        const auto exactM = exactRealRational(call.arguments[2]);
        if (exactN && exactPhi && exactM)
            return CertifiedValue{encloseEllipticPiReal(
                *exactN, *exactPhi, *exactM, precisionBits)};

        const auto n = encloseArgument(0);
        const auto phi = encloseArgument(1);
        const auto m = encloseArgument(2);
        if (!n || !phi || !m)
            return std::nullopt;
        const auto realN = provablyReal(*n);
        const auto realPhi = provablyReal(*phi);
        const auto realM = provablyReal(*m);
        if (!realN || !realPhi || !realM) {
            if (enclosureKind == EnclosureKind::Information
                && (informationStraddlesRealOnlyBackend(*n)
                    || informationStraddlesRealOnlyBackend(*phi)
                    || informationStraddlesRealOnlyBackend(*m)))
                throw PrecisionInsufficient{
                    "Elliptic Pi InformationEnclosure straddles the real-only certified backend boundary",
                    PrecisionInsufficientKind::InputInformation};
            throw CertifiedBackendUnsupported{
                "Certified complex elliptic Pi backend is not implemented"};
        }
        return CertifiedValue{encloseEllipticPiReal(
            *realN, *realPhi, *realM, precisionBits)};
    }

    case BuiltinId::ExponentialIntegralEi:
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
    case BuiltinId::LogarithmicIntegralLi: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;

        const bool informationInput = enclosureKind == EnclosureKind::Information;
        if (definition->id == BuiltinId::ExponentialIntegralEi)
            requireNoPointSingularity(
                *value, rational(0),
                "Ei is singular at zero",
                "Ei InformationEnclosure may contain the singular point zero",
                informationInput);
        else if (definition->id == BuiltinId::CosineIntegralCi)
            requireNoPointSingularity(
                *value, rational(0),
                "Ci is singular at zero",
                "Ci InformationEnclosure may contain the singular point zero",
                informationInput);
        else if (definition->id == BuiltinId::LogarithmicIntegralLi) {
            requireNoPointSingularity(
                *value, rational(1),
                "li is singular at x=1",
                "li InformationEnclosure may contain the singular point x=1",
                informationInput);
            if (informationInput
                && informationCrossesPrincipalNegativeRealCut(value->toComplex()))
                throw PrecisionInsufficient{
                    "li InformationEnclosure cannot determine the principal branch-cut side",
                    PrecisionInsufficientKind::InputInformation};
        }

        if (!value->isReal()) {
            switch (definition->id) {
            case BuiltinId::ExponentialIntegralEi:
                return normalizeCertifiedComplex(encloseExponentialIntegralEiComplex(
                    value->asComplex(), precisionBits));
            case BuiltinId::SineIntegralSi:
                return normalizeCertifiedComplex(encloseSineIntegralSiComplex(
                    value->asComplex(), precisionBits));
            case BuiltinId::CosineIntegralCi:
                return normalizeCertifiedComplex(encloseCosineIntegralCiComplex(
                    value->asComplex(), precisionBits));
            case BuiltinId::LogarithmicIntegralLi:
                return normalizeCertifiedComplex(encloseLogarithmicIntegralLiComplex(
                    value->asComplex(), precisionBits));
            default:
                return std::nullopt;
            }
        }
        switch (definition->id) {
        case BuiltinId::ExponentialIntegralEi:
            return CertifiedValue{encloseExponentialIntegralEiReal(value->asReal(), precisionBits)};
        case BuiltinId::SineIntegralSi:
            return CertifiedValue{encloseSineIntegralSiReal(value->asReal(), precisionBits)};
        case BuiltinId::CosineIntegralCi: {
            const RealInterval& real = value->asReal();
            const Rational lower = real.lower().toRational();
            const Rational upper = real.upper().toRational();
            if (informationInput && lower > Rational{}
                && lower <= rational(96) && upper > rational(96))
                throw PrecisionInsufficient{
                    "Ci InformationEnclosure crosses the certified series boundary x=96",
                    PrecisionInsufficientKind::InputInformation};
            if (lower > Rational{})
                return CertifiedValue{encloseCosineIntegralCiPositive(real, precisionBits)};
            if (upper < Rational{}) {
                // principal branchでは x<0 に対して Ci(x)=Ci(-x)+i Pi。
                // 負の実軸をgeneric complex-logへ送るとbranch cut crossingと区別できないため，
                // 実軸上であることが保証された経路だけこの恒等式でprincipal sideを固定する。
                const RealInterval realPart = encloseCosineIntegralCiPositive(
                    approximation::negate(real), precisionBits);
                return CertifiedValue{ComplexInterval{
                    realPart, enclosePi(precisionBits).interval}};
            }
            if (real.isPoint())
                throw std::domain_error("Ci is singular at zero");
            throw PrecisionInsufficient{
                "Ci real input interval crosses the singular point zero"};
        }
        case BuiltinId::LogarithmicIntegralLi:
            if (value->asReal().upper().toRational() < Rational{})
                return normalizeCertifiedComplex(encloseLogarithmicIntegralLiComplex(
                    ComplexInterval::fromReal(value->asReal()), precisionBits));
            return CertifiedValue{encloseLogarithmicIntegralLiPositive(
                value->asReal(), precisionBits)};
        default:
            return std::nullopt;
        }
    }

    case BuiltinId::Polylog: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        const auto order = exactRealRational(call.arguments[0]);
        if (!order || !order->isInteger() || !order->numerator().isPositive()) {
            // PolyLogは一般の数値orderにも数学的には定義されるが，現certified backendは
            // exact positive integer order専用。数値orderをgeneric unevaluatedへ落とさない。
            if (encloseArgument(0))
                throw CertifiedBackendUnsupported{
                    "Certified polylog currently requires an exact positive integer order"};
            return std::nullopt;
        }
        const auto count = numeric::tryToUint64(order->numerator());
        if (!count)
            throw CertifiedBackendUnsupported{
                "Certified polylog order exceeds the bounded-work backend range"};
        if (const auto z = exactRealRational(call.arguments[1]))
            return CertifiedValue{enclosePolylogReal(*count, *z, precisionBits)};
        const auto z = encloseArgument(1);
        if (!z)
            return std::nullopt;
        if (enclosureKind == EnclosureKind::Information && z->isComplex()
            && informationIsAmbiguousAtPositiveRealCut(z->asComplex()))
            throw PrecisionInsufficient{
                "Polylog InformationEnclosure cannot determine the principal branch-cut side",
                PrecisionInsufficientKind::InputInformation};
        if (z->isReal()) {
            const ComplexInterval result = enclosePolylogComplex(
                *count, z->toComplex(), precisionBits);
            // 現bounded-work領域は|z|<=49/50でprincipal cutの手前にある。
            // 実入力ならpolylogは実数値なのでreal componentだけを返す。
            return CertifiedValue{result.real()};
        }
        return normalizeCertifiedComplex(enclosePolylogComplex(
            *count, z->asComplex(), precisionBits));
    }

    case BuiltinId::IncompleteBeta: {
        if (call.arguments.size() != 3)
            return std::nullopt;

        const auto exactA = exactRealRational(call.arguments[0]);
        const auto exactB = exactRealRational(call.arguments[1]);
        const auto a = encloseArgument(0);
        const auto b = encloseArgument(1);
        const auto x = encloseArgument(2);
        if (!a || !b || !x)
            return std::nullopt;
        if (!a->isReal() || !b->isReal() || !x->isReal())
            throw CertifiedBackendUnsupported{
                "Certified complex ibeta parameters are not implemented"};

        const numeric::BigFloat zero;
        const auto classifyPositiveParameter = [&](
            const RealInterval& value, std::string_view name) {
            if (value.upper() <= zero)
                throw std::domain_error(
                    std::string{"ibeta requires positive "} + std::string{name});
            if (value.lower() <= zero)
                throw PrecisionInsufficient{
                    std::string{"ibeta "} + std::string{name}
                        + " InformationEnclosure crosses the positive-real domain boundary",
                    PrecisionInsufficientKind::InputInformation};
        };
        classifyPositiveParameter(a->asReal(), "a");
        classifyPositiveParameter(b->asReal(), "b");

        const RealInterval& xInterval = x->asReal();
        const Rational xLower = xInterval.lower().toRational();
        const Rational xUpper = xInterval.upper().toRational();
        if (xUpper < rational(0) || xLower > rational(1))
            throw std::domain_error("ibeta requires x in [0,1]");
        if (xLower < rational(0) || xUpper > rational(1))
            throw PrecisionInsufficient{
                "ibeta x InformationEnclosure crosses the x in [0,1] domain boundary",
                PrecisionInsufficientKind::InputInformation};

        if (!exactA || !exactB)
            throw CertifiedBackendUnsupported{
                "Certified ibeta currently requires exact a and b parameters"};
        return CertifiedValue{encloseIncompleteBetaRegularized(
            *exactA, *exactB, xInterval, precisionBits)};
    }

    case BuiltinId::Beta:
    case BuiltinId::BetaLog: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        const auto exactA = exactRealRational(call.arguments[0]);
        const auto exactB = exactRealRational(call.arguments[1]);
        if (exactA && exactB) {
            if (definition->id == BuiltinId::Beta)
                return CertifiedValue{encloseBetaRational(*exactA, *exactB, precisionBits)};
            return CertifiedValue{encloseBetaLogRational(*exactA, *exactB, precisionBits)};
        }

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
    case BuiltinId::IsPrime:
    case BuiltinId::NextPrime:
    case BuiltinId::PreviousPrime:
    case BuiltinId::FactorInteger:
    case BuiltinId::Totient:
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
        // finite-precision入力のInformationEnclosureが負実軸のbranch sideを
        // 決められない場合，両側を含む粗い像を数値として返してはいけない。
        if (enclosureKind == EnclosureKind::Information
            && informationCrossesPrincipalNegativeRealCut(value->asComplex()))
            throw PrecisionInsufficient{
                "sqrt input InformationEnclosure crosses the principal branch cut",
                PrecisionInsufficientKind::InputInformation};
        return normalizeCertifiedComplex(enclosePrincipalComplexSqrt(
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
        return normalizeCertifiedComplex(approximation::divide(
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
        return normalizeCertifiedComplex(ComplexInterval{
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
            if (enclosureKind == EnclosureKind::Information
                && informationCrossesPrincipalNegativeRealCut(base->asComplex()))
                throw PrecisionInsufficient{
                    "Power input InformationEnclosure crosses the principal branch cut",
                    PrecisionInsufficientKind::InputInformation};
            return normalizeCertifiedComplex(enclosePrincipalComplexSqrt(
                base->asComplex(), precisionBits));
        }

        const auto exponent = exactIntegerExponent(call.arguments[1]);
        if (exponent) {
            if (exponent->isZero()) {
                const bool containsZero = base->isReal()
                    ? base->asReal().containsZero()
                    : base->asComplex().containsZero();
                if (enclosureKind == EnclosureKind::Information && containsZero)
                    throw PrecisionInsufficient{
                        "Power base InformationEnclosure may contain zero for exponent zero",
                        PrecisionInsufficientKind::InputInformation};
                return CertifiedValue{RealInterval::fromRational(rational(1), precisionBits)};
            }

            const bool negative = exponent->isNegative();
            const auto magnitude = toUint64(exponent->abs());
            if (!magnitude)
                return std::nullopt;

            CertifiedValue result = integerPower(*base, *magnitude, precisionBits);
            if (!negative)
                return result;

            const CertifiedValue one{RealInterval::fromRational(rational(1), precisionBits)};
            return divideCertifiedValues(one, result, precisionBits);
        }

        // 一般の非整数・複素指数は principal Power の定義へ一本化する。
        // Power(z,w) := Exp(w * principal Log(z))。
        // exact evaluatorでは式を不用意に展開せずPower表記を保つが、Nの内部意味論はこの定義を使用するため、Solverの「全ての根」とは明確に別物になる。
        const auto exponentValue = encloseArgument(1);
        if (!exponentValue)
            return std::nullopt;
        if (enclosureKind == EnclosureKind::Information
            && informationCrossesPrincipalNegativeRealCut(base->toComplex()))
            throw PrecisionInsufficient{
                "Power base InformationEnclosure crosses the principal Log branch cut",
                PrecisionInsufficientKind::InputInformation};
        return normalizeCertifiedComplex(enclosePrincipalPower(
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

        return normalizeCertifiedComplex(encloseComplexExp(
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
            if (enclosureKind == EnclosureKind::Information
                && informationCrossesPrincipalNegativeRealCut(value.toComplex()))
                throw PrecisionInsufficient{
                    "Log input InformationEnclosure crosses the principal branch cut",
                    PrecisionInsufficientKind::InputInformation};
            return normalizeCertifiedComplex(enclosePrincipalComplexLog(
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
        // base=1ならLog[base]=0で数学的に未定義。finite-precision入力の
        // InformationEnclosureが1を含み得る場合はguardを増やしても入力情報自体は
        // 狭まらないため，persistentなprecision不足として即座に返す。
        const bool baseLogMayContainZero = baseLog.isReal()
            ? baseLog.asReal().containsZero()
            : baseLog.asComplex().containsZero();
        if (enclosureKind == EnclosureKind::Information && baseLogMayContainZero)
            throw PrecisionInsufficient{
                "Logarithm base InformationEnclosure may contain one",
                PrecisionInsufficientKind::InputInformation};
        const CertifiedValue reciprocalBaseLog = reciprocalValue(
            baseLog, precisionBits, "Logarithm base could not be certified away from one");
        return multiplyCertifiedValues(valueLog, reciprocalBaseLog, precisionBits);
    }

    case BuiltinId::Arg: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = encloseArgument(0);
        if (!value)
            return std::nullopt;
        if (enclosureKind == EnclosureKind::Information) {
            const ComplexInterval complex = value->toComplex();
            if (informationMayContainComplexZero(complex))
                throw PrecisionInsufficient{
                    "Arg input InformationEnclosure may contain the undefined origin",
                    PrecisionInsufficientKind::InputInformation};
            if (informationCrossesPrincipalNegativeRealCut(complex))
                throw PrecisionInsufficient{
                    "Arg input InformationEnclosure crosses the principal branch cut",
                    PrecisionInsufficientKind::InputInformation};
        }
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
        const auto scalar = encloseBound(angleOperand->value, precisionBits, bindings, enclosureKind, recursionDepth + 1);
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
            return normalizeCertifiedComplex(sine);
        case BuiltinId::Cos:
            return normalizeCertifiedComplex(cosine);
        case BuiltinId::Tan:
            if (cosine.containsZero())
                throw PrecisionInsufficient{"Complex tangent denominator cannot yet be proven nonzero"};
            return normalizeCertifiedComplex(approximation::divide(sine, cosine, precisionBits));
        case BuiltinId::Cot:
            if (sine.containsZero())
                throw PrecisionInsufficient{"Complex cotangent denominator cannot yet be proven nonzero"};
            return normalizeCertifiedComplex(approximation::divide(cosine, sine, precisionBits));
        case BuiltinId::Sec:
            return reciprocalValue(normalizeCertifiedComplex(cosine), precisionBits,
                "Complex secant denominator cannot yet be proven nonzero");
        case BuiltinId::Csc:
            return reciprocalValue(normalizeCertifiedComplex(sine), precisionBits,
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

        if (enclosureKind == EnclosureKind::Information && input->isComplex()) {
            const ComplexInterval complex = input->asComplex();
            const bool crossesCut = definition->id == BuiltinId::Atan
                ? informationIsAmbiguousAtOuterImaginaryCuts(complex)
                : informationIsAmbiguousAtOuterRealCuts(complex);
            if (crossesCut)
                throw PrecisionInsufficient{
                    "Inverse trigonometric InformationEnclosure cannot determine the principal branch-cut side",
                    PrecisionInsufficientKind::InputInformation};
        }

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
                return normalizeCertifiedComplex(enclosePrincipalComplexAsinRadian(complex, precisionBits));
            if (definition->id == BuiltinId::Acos)
                return normalizeCertifiedComplex(enclosePrincipalComplexAcosRadian(complex, precisionBits));
            return normalizeCertifiedComplex(enclosePrincipalComplexAtanRadian(complex, precisionBits));
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

        const ComplexInterval cartesian{x->asReal(), y->asReal()};
        if (enclosureKind == EnclosureKind::Information) {
            if (informationMayContainComplexZero(cartesian))
                throw PrecisionInsufficient{
                    "atan2 InformationEnclosure may contain the undefined origin",
                    PrecisionInsufficientKind::InputInformation};
            if (informationCrossesPrincipalNegativeRealCut(cartesian))
                throw PrecisionInsufficient{
                    "atan2 InformationEnclosure cannot determine the principal branch-cut side",
                    PrecisionInsufficientKind::InputInformation};
        }

        const RealInterval radians = enclosePrincipalArgument(
            cartesian, precisionBits).interval;
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
            return normalizeCertifiedComplex(sinhValue);
        case BuiltinId::Cosh:
            return normalizeCertifiedComplex(coshValue);
        case BuiltinId::Tanh:
            if (coshValue.containsZero())
                throw PrecisionInsufficient{"Complex tanh denominator cannot yet be proven nonzero"};
            return normalizeCertifiedComplex(approximation::divide(sinhValue, coshValue, precisionBits));
        case BuiltinId::Csch:
            return reciprocalValue(normalizeCertifiedComplex(sinhValue), precisionBits,
                "Complex csch denominator cannot yet be proven nonzero");
        case BuiltinId::Sech:
            return reciprocalValue(normalizeCertifiedComplex(coshValue), precisionBits,
                "Complex sech denominator cannot yet be proven nonzero");
        case BuiltinId::Coth:
            if (sinhValue.containsZero())
                throw PrecisionInsufficient{"Complex coth denominator cannot yet be proven nonzero"};
            return normalizeCertifiedComplex(approximation::divide(coshValue, sinhValue, precisionBits));
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

        if (enclosureKind == EnclosureKind::Information && input->isComplex()) {
            const ComplexInterval complex = input->asComplex();
            bool crossesCut = false;
            if (definition->id == BuiltinId::Asinh)
                crossesCut = informationIsAmbiguousAtOuterImaginaryCuts(complex);
            else if (definition->id == BuiltinId::Acosh)
                crossesCut = informationIsAmbiguousAtAcoshCut(complex);
            else
                crossesCut = informationIsAmbiguousAtOuterRealCuts(complex);
            if (crossesCut)
                throw PrecisionInsufficient{
                    "Inverse hyperbolic InformationEnclosure cannot determine the principal branch-cut side",
                    PrecisionInsufficientKind::InputInformation};
        }

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
            return normalizeCertifiedComplex(enclosePrincipalComplexAsinh(complex, precisionBits));
        if (definition->id == BuiltinId::Acosh)
            return normalizeCertifiedComplex(enclosePrincipalComplexAcosh(complex, precisionBits));
        return normalizeCertifiedComplex(enclosePrincipalComplexAtanh(complex, precisionBits));
    }

    case BuiltinId::Fma: {
        if (call.arguments.size() != 3)
            return std::nullopt;
        const auto a = encloseArgument(0);
        const auto b = encloseArgument(1);
        const auto c = encloseArgument(2);
        if (!a || !b || !c)
            return std::nullopt;
        return addCertifiedValues(multiplyCertifiedValues(*a, *b, precisionBits), *c, precisionBits);
    }

    case BuiltinId::Clamp: {
        if (call.arguments.size() != 3)
            return std::nullopt;
        const auto x = encloseArgument(0);
        const auto lo = encloseArgument(1);
        const auto hi = encloseArgument(2);
        if (!x || !lo || !hi || !x->isReal() || !lo->isReal() || !hi->isReal())
            return std::nullopt;
        if (hi->asReal().lower() < lo->asReal().upper())
            throw PrecisionInsufficient{"clamp bounds cannot yet be ordered"};
        const auto clampValue = [](const numeric::BigFloat& value,
                                   const numeric::BigFloat& lower,
                                   const numeric::BigFloat& upper) {
            if (value < lower)
                return lower;
            if (upper < value)
                return upper;
            return value;
        };
        numeric::BigFloat lower = clampValue(
            x->asReal().lower(), lo->asReal().lower(), hi->asReal().lower());
        numeric::BigFloat upper = clampValue(
            x->asReal().upper(), lo->asReal().upper(), hi->asReal().upper());
        return CertifiedValue{RealInterval{std::move(lower), std::move(upper)}};
    }

    case BuiltinId::Proj:
        if (call.arguments.size() != 1)
            return std::nullopt;
        return encloseArgument(0);

    case BuiltinId::Polar:
    case BuiltinId::NextPow2:
    case BuiltinId::BitAnd:
    case BuiltinId::BitOr:
    case BuiltinId::BitXor:
    case BuiltinId::BitNot:
    case BuiltinId::BitShiftLeft:
    case BuiltinId::BitShiftRight:
    case BuiltinId::BitLength:
    case BuiltinId::BitCount:
    case BuiltinId::BitGet:
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
    case BuiltinId::Dimensions:
    case BuiltinId::ArrayRank:
    case BuiltinId::ArrayGet:
    case BuiltinId::Reshape:
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
    case BuiltinId::SolveLinear:
    case BuiltinId::NullSpace:
    case BuiltinId::LuDecomposition:
    case BuiltinId::QrDecomposition:
    case BuiltinId::SingularValueDecomposition:
    case BuiltinId::ConditionNumber:
    case BuiltinId::LeastSquares:
    case BuiltinId::PseudoInverse:
    case BuiltinId::Eigenvalues:
    case BuiltinId::Eigenvectors:
    case BuiltinId::Eigensystem:
    case BuiltinId::ConjugateTranspose:
    case BuiltinId::Length:
    case BuiltinId::NumericDerivative:
    case BuiltinId::NumericIntegral:
    case BuiltinId::Precision:
    case BuiltinId::Accuracy:
    case BuiltinId::Explain:
    case BuiltinId::Map:
    case BuiltinId::Range:
    case BuiltinId::Table:
    case BuiltinId::Rationalize:
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify:
    case BuiltinId::Expand:
    case BuiltinId::Factor:
    case BuiltinId::Collect:
    case BuiltinId::Solve:
    case BuiltinId::GroebnerBasis:
    case BuiltinId::PolynomialReduce:
    case BuiltinId::Cases:
    case BuiltinId::CaseBranch:
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
