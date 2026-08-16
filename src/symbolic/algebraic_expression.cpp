// Expr表現とexact algebraic-number backendの橋渡し
#include "algebraic_expression.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/math_ids.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace mmcal::symbolic {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Rational;

constexpr std::size_t maximumBridgeNodes = 64;
constexpr std::uint64_t maximumIntegerPowerMagnitude = 32;

[[nodiscard]] std::optional<AlgebraicNumber> rationalValue(const Rational& value) {
    return AlgebraicNumber::fromRational(value);
}

[[nodiscard]] std::optional<AlgebraicNumber> combine(
    const AlgebraicNumber& lhs,
    const AlgebraicNumber& rhs,
    AlgebraicBinaryOperation operation) {
    return AlgebraicNumber::combine(lhs, rhs, operation);
}


[[nodiscard]] std::optional<AlgebraicNumber> positiveRationalSquareRoot(
    const Rational& value) {
    if (value <= Rational{BigInt{0}})
        return std::nullopt;

    const std::vector<Rational> polynomial{
        -value,
        Rational{},
        Rational{BigInt{1}}
    };
    return AlgebraicNumber::create(polynomial, 2, AlgebraicRootDomain::Real);
}

[[nodiscard]] std::optional<AlgebraicNumber> positiveRationalCubeRoot(
    const Rational& value) {
    if (value <= Rational{BigInt{0}})
        return std::nullopt;

    const std::vector<Rational> polynomial{
        -value,
        Rational{},
        Rational{},
        Rational{BigInt{1}}
    };
    return AlgebraicNumber::create(polynomial, 1, AlgebraicRootDomain::Real);
}

[[nodiscard]] std::optional<AlgebraicNumber> phiValue() {
    // Phi=(1+sqrt(5))/2 は x^2-x-1 の正の根。
    const std::vector<Rational> polynomial{
        Rational{BigInt{-1}},
        Rational{BigInt{-1}},
        Rational{BigInt{1}}
    };
    return AlgebraicNumber::create(polynomial, 2, AlgebraicRootDomain::Real);
}

[[nodiscard]] std::optional<AlgebraicNumber> integerPower(
    const AlgebraicNumber& base,
    const BigInt& exponent) {
    const auto magnitude = numeric::tryToUint64(exponent.abs());
    if (!magnitude || *magnitude > maximumIntegerPowerMagnitude)
        return std::nullopt;

    auto one = rationalValue(Rational{BigInt{1}});
    if (!one)
        return std::nullopt;
    if (*magnitude == 0) {
        // held式をbridgeするときも0^0を1へ決め打ちしない。
        // baseがexact nonzeroと証明できる場合だけ通常のa^0=1を適用する。
        const auto zero = rationalValue(Rational{});
        if (!zero)
            return std::nullopt;
        const auto equalsZero = base.exactEquals(*zero);
        if (!equalsZero || *equalsZero)
            return std::nullopt;
        return one;
    }

    AlgebraicNumber factor = base;
    AlgebraicNumber result = *one;
    std::uint64_t power = *magnitude;
    while (power != 0) {
        if ((power & 1U) != 0U) {
            auto multiplied = combine(result, factor, AlgebraicBinaryOperation::Multiply);
            if (!multiplied)
                return std::nullopt;
            result = std::move(*multiplied);
        }
        power >>= 1U;
        if (power == 0)
            break;
        auto squared = combine(factor, factor, AlgebraicBinaryOperation::Multiply);
        if (!squared)
            return std::nullopt;
        factor = std::move(*squared);
    }

    if (!exponent.isNegative())
        return result;
    return combine(*one, result, AlgebraicBinaryOperation::Divide);
}

[[nodiscard]] std::optional<AlgebraicNumber> exactAlgebraicValueImpl(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    std::size_t& remainingNodes) {
    if (remainingNodes == 0)
        return std::nullopt;
    --remainingNodes;

    if (expression.isCall() && expression.asCall().algebraicValue)
        return *expression.asCall().algebraicValue;

    if (expression.isNumber()) {
        const auto& number = expression.asNumber();
        if (number.isReal())
            return AlgebraicNumber::fromRational(number.asReal().toRational());
        return AlgebraicNumber::fromComplexRational(
            number.asComplex().real.toRational(),
            number.asComplex().imaginary.toRational());
    }

    if (expression.isSymbol()) {
        const auto* constant = mathematics.findConstant(expression.asSymbol());
        if (constant && constant->id == mathematics::ConstantId::Phi)
            return phiValue();
        return std::nullopt;
    }

    if (!expression.isCall())
        return std::nullopt;

    const auto& call = expression.asCall();
    const auto* definition = builtins.find(call.head);
    if (!definition)
        return std::nullopt;

    const auto child = [&](std::size_t index) {
        return exactAlgebraicValueImpl(
            call.arguments[index], builtins, mathematics, remainingNodes);
    };

    switch (definition->id) {
    case BuiltinId::Root: {
        if (call.arguments.size() < 2 || call.arguments.size() > 3)
            return std::nullopt;
        const auto coefficients = rootPolynomialCoefficients(call.arguments[0]);
        const auto index = positiveRootIndex(call.arguments[1]);
        if (!coefficients || coefficients->size() < 2 || !index)
            return std::nullopt;
        if (coefficients->size() == 2) {
            if (*index != 1 || (*coefficients)[1].isZero())
                return std::nullopt;
            return rationalValue(-(*coefficients)[0] / (*coefficients)[1]);
        }

        AlgebraicRootDomain domain = AlgebraicRootDomain::Real;
        if (call.arguments.size() == 3) {
            if (!call.arguments[2].isSymbol()
                || call.arguments[2].asSymbol().view() != "Complex")
                return std::nullopt;
            domain = AlgebraicRootDomain::Complex;
        }
        return AlgebraicNumber::create(*coefficients, *index, domain);
    }

    case BuiltinId::Negate: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        const auto value = child(0);
        const auto minusOne = rationalValue(Rational{BigInt{-1}});
        if (!value || !minusOne)
            return std::nullopt;
        return combine(*minusOne, *value, AlgebraicBinaryOperation::Multiply);
    }

    case BuiltinId::Add:
    case BuiltinId::Multiply: {
        if (call.arguments.empty())
            return std::nullopt;
        auto result = child(0);
        if (!result)
            return std::nullopt;
        const AlgebraicBinaryOperation operation = definition->id == BuiltinId::Add
            ? AlgebraicBinaryOperation::Add
            : AlgebraicBinaryOperation::Multiply;
        for (std::size_t i = 1; i < call.arguments.size(); ++i) {
            auto rhs = child(i);
            if (!rhs)
                return std::nullopt;
            auto next = combine(*result, *rhs, operation);
            if (!next)
                return std::nullopt;
            result = std::move(next);
        }
        return result;
    }

    case BuiltinId::Subtract:
    case BuiltinId::Divide: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        const auto lhs = child(0);
        const auto rhs = child(1);
        if (!lhs || !rhs)
            return std::nullopt;
        return combine(
            *lhs,
            *rhs,
            definition->id == BuiltinId::Subtract
                ? AlgebraicBinaryOperation::Subtract
                : AlgebraicBinaryOperation::Divide);
    }

    case BuiltinId::Sqrt: {
        if (call.arguments.size() != 1
            || !call.arguments[0].isNumber()
            || !call.arguments[0].asNumber().isReal())
            return std::nullopt;
        const Rational value = call.arguments[0].asNumber().asReal().toRational();
        return positiveRationalSquareRoot(value);
    }

    case BuiltinId::Cbrt: {
        if (call.arguments.size() != 1
            || !call.arguments[0].isNumber()
            || !call.arguments[0].asNumber().isReal())
            return std::nullopt;
        const Rational value = call.arguments[0].asNumber().asReal().toRational();
        return positiveRationalCubeRoot(value);
    }

    case BuiltinId::Power: {
        if (call.arguments.size() != 2
            || !call.arguments[1].isNumber()
            || !call.arguments[1].asNumber().isReal()
            || !call.arguments[1].asNumber().asReal().isInteger())
            return std::nullopt;
        const auto base = child(0);
        if (!base)
            return std::nullopt;
        return integerPower(
            *base,
            call.arguments[1].asNumber().asReal().asInteger());
    }

    default:
        return std::nullopt;
    }
}

} // namespace

std::optional<std::vector<Rational>> rootPolynomialCoefficients(
    const Expr& expression) {
    if (!expression.isArray() || expression.asArray().rank() != 1)
        return std::nullopt;

    const auto& array = expression.asArray();
    std::vector<Rational> coefficients;
    coefficients.reserve(array.size());
    for (std::size_t i = 0; i < array.size(); ++i) {
        const Expr value = array.element(i);
        if (!value.isNumber() || !value.asNumber().isReal())
            return std::nullopt;
        coefficients.push_back(value.asNumber().asReal().toRational());
    }
    while (!coefficients.empty() && coefficients.back().isZero())
        coefficients.pop_back();
    return coefficients;
}

std::optional<std::size_t> positiveRootIndex(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    const BigInt& integer = expression.asNumber().asReal().asInteger();
    if (integer.isNegative() || integer.isZero())
        return std::nullopt;
    const auto value = numeric::tryToUint64(integer);
    if (!value || *value > std::numeric_limits<std::size_t>::max())
        return std::nullopt;
    return static_cast<std::size_t>(*value);
}

Expr makeCanonicalRootExpression(
    const RealAlgebraicNumber& algebraic,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Rational> coefficients(
        algebraic.polynomial().begin(), algebraic.polynomial().end());
    const std::size_t coefficientCount = coefficients.size();
    const AlgebraicNumber cached =
        AlgebraicNumber::fromRealRoot(algebraic).withGeneratorField();
    return Expr::call(builtins.symbol(BuiltinId::Root), {
        Expr::rationalArray({coefficientCount}, std::move(coefficients)),
        Expr{numeric::Number{BigInt::fromUnsigned(algebraic.rootIndex())}}
    }, std::make_shared<const AlgebraicNumber>(cached));
}

Expr makeCanonicalRootExpression(
    const ComplexAlgebraicNumber& algebraic,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Rational> coefficients(
        algebraic.polynomial().begin(), algebraic.polynomial().end());
    const std::size_t coefficientCount = coefficients.size();
    const AlgebraicNumber cached =
        AlgebraicNumber::fromComplexRoot(algebraic).withGeneratorField();
    return Expr::call(builtins.symbol(BuiltinId::Root), {
        Expr::rationalArray({coefficientCount}, std::move(coefficients)),
        Expr{numeric::Number{BigInt::fromUnsigned(algebraic.rootIndex())}},
        Expr{expression::Symbol{"Complex"}}
    }, std::make_shared<const AlgebraicNumber>(cached));
}

std::optional<AlgebraicNumber> exactAlgebraicValue(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    std::size_t remainingNodes = maximumBridgeNodes;
    return exactAlgebraicValueImpl(
        expression, builtins, mathematics, remainingNodes);
}

} // namespace mmcal::symbolic
