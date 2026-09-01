// 四則演算、冪
#include "arithmetic.hpp"
#include "builtin_helpers.hpp"

#include "approximation/expression_interval.hpp"
#include "error/error_message.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "names.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/integer_algorithms.hpp"
#include "symbolic/algebraic_number.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] std::size_t realNumberBits(const RealNumber& value) {
    if (value.isInteger())
        return value.asInteger().bitLength();
    const Rational& rational = value.asRational();
    return std::max(
        rational.numerator().bitLength(),
        rational.denominator().bitLength());
}

[[nodiscard]] std::size_t numberBits(const Number& value) {
    if (value.isReal())
        return realNumberBits(value.asReal());
    return std::max(
        realNumberBits(value.asComplex().real),
        realNumberBits(value.asComplex().imaginary));
}

void checkNumberBudget(const Number& number) {
    evaluation::checkEvaluationBigIntegerBits(numberBits(number));
}

[[nodiscard]] Expr numberExpr(Number number) {
    checkNumberBudget(number);
    return Expr{std::move(number)};
}

[[nodiscard]] Expr numberArray(
    std::vector<std::size_t> shape,
    std::vector<Number> values) {
    for (const Number& value : values)
        checkNumberBudget(value);
    return Expr::numberArray(std::move(shape), std::move(values));
}

[[nodiscard]] Expr integerExpr(std::int64_t value) {
    return numberExpr(Number{BigInt{value}});
}

[[nodiscard]] bool isExactZero(const Expr& value) noexcept {
    return value.isNumber() && value.asNumber().isZero();
}

[[nodiscard]] bool isExactOne(const Expr& value) noexcept {
    return value.isNumber() && value.asNumber() == Number{BigInt{1}};
}

[[nodiscard]] std::optional<std::uint64_t> toUint64(const BigInt& value) {
    if (value.isNegative())
        return std::nullopt;

    const std::string text = value.toString();
    std::uint64_t result = 0;
    const auto conversion = std::from_chars(text.data(), text.data() + text.size(), result);
    if (conversion.ec != std::errc{} || conversion.ptr != text.data() + text.size())
        return std::nullopt;

    return result;
}



[[nodiscard]] std::optional<symbolic::AlgebraicNumber> algebraicValue(
    const Expr& expression,
    const evaluation::BuiltinRegistry& registry) {
    if (expression.isCall() && expression.asCall().algebraicValue)
        return expression.asCall().algebraicValue->withGeneratorField();
    if (expression.isNumber()) {
        const Number& number = expression.asNumber();
        if (number.isReal())
            return symbolic::AlgebraicNumber::fromRational(number.asReal().toRational());
        return symbolic::AlgebraicNumber::fromComplexRational(
            number.asComplex().real.toRational(), number.asComplex().imaginary.toRational());
    }
    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = registry.find(expression.asCall().head);
    if (!definition || definition->id != evaluation::BuiltinId::Root)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if ((arguments.size() != 2 && arguments.size() != 3)
        || !arguments[0].isArray() || arguments[0].asArray().rank() != 1)
        return std::nullopt;

    std::vector<Rational> coefficients;
    coefficients.reserve(arguments[0].asArray().size());
    for (std::size_t i = 0; i < arguments[0].asArray().size(); ++i) {
        const Expr value = arguments[0].asArray().element(i);
        if (!value.isNumber() || !value.asNumber().isReal())
            return std::nullopt;
        coefficients.push_back(value.asNumber().asReal().toRational());
    }
    if (!arguments[1].isNumber() || !arguments[1].asNumber().isReal()
        || !arguments[1].asNumber().asReal().isInteger())
        return std::nullopt;
    const BigInt& indexInteger = arguments[1].asNumber().asReal().asInteger();
    const auto index = numeric::tryToUint64(indexInteger);
    if (indexInteger.isNegative() || indexInteger.isZero() || !index
        || *index > std::numeric_limits<std::size_t>::max())
        return std::nullopt;
    symbolic::AlgebraicRootDomain domain = symbolic::AlgebraicRootDomain::Real;
    if (arguments.size() == 3) {
        if (!arguments[2].isSymbol() || arguments[2].asSymbol().view() != "Complex")
            return std::nullopt;
        domain = symbolic::AlgebraicRootDomain::Complex;
    }
    auto algebraic = symbolic::AlgebraicNumber::create(
        coefficients, static_cast<std::size_t>(*index), domain);
    if (!algebraic)
        return std::nullopt;
    return algebraic->withGeneratorField();
}

[[nodiscard]] Expr algebraicExpr(
    const symbolic::AlgebraicNumber& value,
    const evaluation::BuiltinRegistry& registry) {
    if (const auto exact = value.exactRationalParts()) {
        if (exact->second.isZero())
            return Expr{Number{exact->first}};
        return Expr{Number::complex(RealNumber{exact->first}, RealNumber{exact->second})};
    }
    const auto polynomial = value.polynomial();
    if (polynomial.size() == 2)
        return Expr{Number{-polynomial[0] / polynomial[1]}};

    std::vector<Rational> coefficients(polynomial.begin(), polynomial.end());
    const std::size_t coefficientCount = coefficients.size();
    std::vector<Expr> arguments;
    arguments.reserve(value.domain() == symbolic::AlgebraicRootDomain::Complex ? 3 : 2);
    arguments.push_back(Expr::rationalArray({coefficientCount}, std::move(coefficients)));
    arguments.push_back(Expr{Number{BigInt::fromUnsigned(value.rootIndex())}});
    if (value.domain() == symbolic::AlgebraicRootDomain::Complex)
        arguments.push_back(Expr{expression::Symbol{"Complex"}});
    const symbolic::AlgebraicNumber cached = value.withGeneratorField();
    return Expr::call(
        registry.symbol(evaluation::BuiltinId::Root),
        std::move(arguments),
        std::make_shared<const symbolic::AlgebraicNumber>(cached));
}

[[nodiscard]] std::optional<Expr> algebraicBinary(
    const Expr& lhs,
    const Expr& rhs,
    symbolic::AlgebraicBinaryOperation operation,
    const evaluation::BuiltinRegistry& registry) {
    const auto left = algebraicValue(lhs, registry);
    if (!left)
        return std::nullopt;
    const auto right = algebraicValue(rhs, registry);
    if (!right)
        return std::nullopt;
    const auto result = symbolic::AlgebraicNumber::combine(*left, *right, operation);
    if (!result)
        return std::nullopt;
    return algebraicExpr(*result, registry);
}

[[nodiscard]] std::optional<Expr> algebraicFold(
    const std::vector<Expr>& arguments,
    symbolic::AlgebraicBinaryOperation operation,
    const evaluation::BuiltinRegistry& registry) {
    if (arguments.empty())
        return std::nullopt;
    bool containsRoot = false;
    for (const Expr& argument : arguments) {
        if (argument.isCall()) {
            const auto* definition = registry.find(argument.asCall().head);
            containsRoot = containsRoot || (definition && definition->id == evaluation::BuiltinId::Root);
        }
        if (!algebraicValue(argument, registry))
            return std::nullopt;
    }
    if (!containsRoot)
        return std::nullopt;

    auto accumulated = algebraicValue(arguments.front(), registry);
    for (std::size_t i = 1; i < arguments.size(); ++i) {
        const auto next = algebraicValue(arguments[i], registry);
        accumulated = symbolic::AlgebraicNumber::combine(*accumulated, *next, operation);
        if (!accumulated)
            return std::nullopt;
    }
    return algebraicExpr(*accumulated, registry);
}

[[nodiscard]] Expr scalarAdd(
    const std::vector<Expr>& arguments,
    const evaluation::BuiltinRegistry& registry) {
    Number sum{BigInt{0}};
    bool numeric = true;
    for (const Expr& argument : arguments) {
        if (!argument.isNumber()) {
            numeric = false;
            break;
        }
        sum += argument.asNumber();
        checkNumberBudget(sum);
    }
    if (numeric)
        return numberExpr(std::move(sum));

    const Expr* soleNonZero = nullptr;
    bool additiveIdentityOnly = true;
    for (const Expr& argument : arguments) {
        if (isExactZero(argument))
            continue;
        if (soleNonZero) {
            additiveIdentityOnly = false;
            break;
        }
        soleNonZero = &argument;
    }
    if (additiveIdentityOnly && soleNonZero)
        return *soleNonZero;

    if (const auto algebraic = algebraicFold(arguments, symbolic::AlgebraicBinaryOperation::Add, registry))
        return *algebraic;
    if (const auto approximate = approximation::addApproximateScalars(arguments))
        return *approximate;
    return Expr::call(registry.symbol(evaluation::BuiltinId::Add), arguments);
}

[[nodiscard]] Expr scalarMultiply(
    const std::vector<Expr>& arguments,
    const evaluation::BuiltinRegistry& registry) {
    Number product{BigInt{1}};
    bool numeric = true;
    for (const Expr& argument : arguments) {
        if (!argument.isNumber()) {
            numeric = false;
            break;
        }
        product *= argument.asNumber();
        checkNumberBudget(product);
    }
    if (numeric)
        return numberExpr(std::move(product));

    const Expr* soleNonOne = nullptr;
    bool multiplicativeIdentityOnly = true;
    for (const Expr& argument : arguments) {
        if (isExactOne(argument))
            continue;
        if (soleNonOne) {
            multiplicativeIdentityOnly = false;
            break;
        }
        soleNonOne = &argument;
    }
    if (multiplicativeIdentityOnly && soleNonOne)
        return *soleNonOne;

    if (const auto algebraic = algebraicFold(arguments, symbolic::AlgebraicBinaryOperation::Multiply, registry))
        return *algebraic;
    if (const auto approximate = approximation::multiplyApproximateScalars(arguments))
        return *approximate;
    return Expr::call(registry.symbol(evaluation::BuiltinId::Multiply), arguments);
}

[[nodiscard]] Expr scalarSubtract(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& registry) {
    if (lhs.isNumber() && rhs.isNumber())
        return numberExpr(lhs.asNumber() - rhs.asNumber());
    if (isExactZero(rhs))
        return lhs;
    if (const auto algebraic = algebraicBinary(
        lhs, rhs, symbolic::AlgebraicBinaryOperation::Subtract, registry))
        return *algebraic;
    if (const auto approximate = approximation::subtractApproximateScalars(lhs, rhs))
        return *approximate;
    return Expr::call(registry.symbol(evaluation::BuiltinId::Subtract), {lhs, rhs});
}

[[nodiscard]] Expr scalarNegate(
    const Expr& value,
    const evaluation::BuiltinRegistry& registry) {
    if (value.isNumber())
        return numberExpr(-value.asNumber());
    if (const auto zero = symbolic::AlgebraicNumber::fromRational(Rational{})) {
        if (const auto algebraic = algebraicValue(value, registry)) {
            if (const auto result = symbolic::AlgebraicNumber::combine(
                *zero, *algebraic, symbolic::AlgebraicBinaryOperation::Subtract))
                return algebraicExpr(*result, registry);
        }
    }
    if (const auto approximate = approximation::negateApproximateScalar(value))
        return *approximate;
    return Expr::call(registry.symbol(evaluation::BuiltinId::Negate), {value});
}

[[nodiscard]] bool isInformationExactZero(const Expr& value) noexcept {
    if (value.isNumber())
        return value.asNumber().isZero();
    if (value.isDecimalApproximation())
        return value.asDecimalApproximation().informationExactlyZero();
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        return complex.realInformationExactlyZero()
            && complex.imaginaryInformationExactlyZero();
    }
    return false;
}

[[noreturn]] void arrayArithmeticError(std::string message) {
    error::throwCalcError(error::CalcErrorType::Type, std::move(message));
}

} // namespace

Expr evaluateAdd(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    const bool hasArray = std::any_of(arguments.begin(), arguments.end(),
        [](const Expr& argument) { return argument.isArray(); });
    if (!hasArray)
        return scalarAdd(std::vector<Expr>{arguments.begin(), arguments.end()}, registry);

    for (const Expr& argument : arguments)
        if (!argument.isArray())
            arrayArithmeticError("Array addition requires arrays with identical shapes");

    const auto& first = arguments.front().asArray();
    for (const Expr& argument : arguments)
        if (argument.asArray().shape != first.shape)
            error::throwCalcError(error::CalcErrorType::Domain,
                "Array addition requires identical shapes");

    const bool allExact = std::all_of(arguments.begin(), arguments.end(),
        [](const Expr& argument) { return argument.asArray().hasExactNumberStorage(); });
    if (allExact) {
        std::vector<Number> values(first.size(), Number{BigInt{0}});
        for (const Expr& argument : arguments) {
            const auto& array = argument.asArray();
            for (std::size_t i = 0; i < values.size(); ++i)
                values[i] += array.exactNumber(i);
        }
        return numberArray(first.shape, std::move(values));
    }

    std::vector<Expr> elements;
    elements.reserve(first.size());
    for (std::size_t i = 0; i < first.size(); ++i) {
        std::vector<Expr> terms;
        terms.reserve(arguments.size());
        for (const Expr& argument : arguments)
            terms.push_back(argument.asArray().element(i));
        elements.push_back(scalarAdd(terms, registry));
    }
    return Expr::array(first.shape, std::move(elements));
}

Expr evaluateSubtract(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::subtract);
    const Expr& lhs = arguments[0];
    const Expr& rhs = arguments[1];

    if (lhs.isArray() || rhs.isArray()) {
        if (!lhs.isArray() || !rhs.isArray())
            arrayArithmeticError("Array subtraction requires two arrays with identical shapes");
        if (lhs.asArray().shape != rhs.asArray().shape)
            error::throwCalcError(error::CalcErrorType::Domain,
                "Array subtraction requires identical shapes");

        if (lhs.asArray().hasExactNumberStorage() && rhs.asArray().hasExactNumberStorage()) {
            std::vector<Number> values;
            values.reserve(lhs.asArray().size());
            for (std::size_t i = 0; i < lhs.asArray().size(); ++i)
                values.push_back(lhs.asArray().exactNumber(i) - rhs.asArray().exactNumber(i));
            return numberArray(lhs.asArray().shape, std::move(values));
        }

        std::vector<Expr> elements;
        elements.reserve(lhs.asArray().size());
        for (std::size_t i = 0; i < lhs.asArray().size(); ++i) {
            const Expr left = lhs.asArray().element(i);
            const Expr right = rhs.asArray().element(i);
            elements.push_back(scalarSubtract(left, right, registry));
        }
        return Expr::array(lhs.asArray().shape, std::move(elements));
    }

    return scalarSubtract(lhs, rhs, registry);
}

Expr evaluateMultiply(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    std::size_t arrayCount = 0;
    std::size_t arrayIndex = 0;
    for (std::size_t i = 0; i < arguments.size(); ++i)
        if (arguments[i].isArray()) {
            ++arrayCount;
            arrayIndex = i;
        }

    if (arrayCount == 0)
        return scalarMultiply(std::vector<Expr>{arguments.begin(), arguments.end()}, registry);
    if (arrayCount > 1)
        arrayArithmeticError(
            "Array multiplication is scalar-only; use dot[...] for vector or matrix contraction");

    const auto& array = arguments[arrayIndex].asArray();
    std::vector<Expr> scalarFactors;
    scalarFactors.reserve(arguments.size() - 1);
    for (std::size_t i = 0; i < arguments.size(); ++i)
        if (i != arrayIndex)
            scalarFactors.push_back(arguments[i]);

    if (array.hasExactNumberStorage()
        && std::all_of(scalarFactors.begin(), scalarFactors.end(),
            [](const Expr& value) { return value.isNumber(); })) {
        Number scalar{BigInt{1}};
        for (const Expr& factor : scalarFactors)
            scalar *= factor.asNumber();
        std::vector<Number> values;
        values.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i)
            values.push_back(array.exactNumber(i) * scalar);
        return numberArray(array.shape, std::move(values));
    }

    std::vector<Expr> elements;
    elements.reserve(array.size());
    for (std::size_t i = 0; i < array.size(); ++i) {
        const Expr element = array.element(i);
        std::vector<Expr> factors;
        factors.reserve(scalarFactors.size() + 1);
        factors.push_back(element);
        factors.insert(factors.end(), scalarFactors.begin(), scalarFactors.end());
        elements.push_back(scalarMultiply(factors, registry));
    }
    return Expr::array(array.shape, std::move(elements));
}

Expr evaluateDivide(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 2, names::divide);

    const Expr& numerator = arguments[0];
    const Expr& denominator = arguments[1];
    if (isInformationExactZero(denominator))
        error::throwCalcError(error::CalcErrorType::Domain, "Division by zero");
    if (isExactOne(denominator))
        return numerator;

    if (numerator.isArray() || denominator.isArray()) {
        if (!numerator.isArray() || denominator.isArray())
            arrayArithmeticError(
                "Array division is scalar-only; only array/scalar is supported");

        const auto& array = numerator.asArray();
        if (array.hasExactNumberStorage() && denominator.isNumber()) {
            std::vector<Number> values;
            values.reserve(array.size());
            for (std::size_t i = 0; i < array.size(); ++i)
                values.push_back(array.exactNumber(i) / denominator.asNumber());
            return numberArray(array.shape, std::move(values));
        }

        std::vector<Expr> elements;
        elements.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i) {
            const std::array<Expr, 2> pair{array.element(i), denominator};
            elements.push_back(evaluateDivide(pair, registry));
        }
        return Expr::array(array.shape, std::move(elements));
    }
    if (numerator.isNumber() && denominator.isNumber())
        return numberExpr(numerator.asNumber() / denominator.asNumber());
    if (const auto algebraic = algebraicBinary(
        numerator, denominator, symbolic::AlgebraicBinaryOperation::Divide, registry))
        return *algebraic;
    if (const auto approximate = approximation::divideApproximateScalars(
        numerator, denominator))
        return *approximate;

    return Expr::call(
        registry.symbol(evaluation::BuiltinId::Divide),
        {numerator, denominator});
}

Expr evaluatePower(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireArity(arguments, 2, names::power);

    const Expr& base = arguments[0];
    const Expr& exponent = arguments[1];

    if (base.isNumber() && base.asNumber().isZero() && exponent.isNumber()
        && exponent.asNumber().isReal()) {
        const RealNumber& realExponent = exponent.asNumber().asReal();
        if (realExponent.isZero())
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "Zero to the zero power is indeterminate");
        if (realExponent.isNegative())
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "Zero cannot be raised to a negative power");

        return integerExpr(0);
    }

    if (base.isNumber() && base.asNumber() == Number{BigInt{1}})
        return integerExpr(1);

    // mmCalでは有限小数もexact Rationalとして読むため、0.5 は厳密に 1/2。
    // Power[x, 1/2] はprincipal square rootの定義と一致させる。これにより
    //     (-2)^0.5  -> I * sqrt[2]
    //     4^0.5     -> 2
    //     x^0.5     -> sqrt[x]
    // となり、「小数で書いたからmachine-real powerへ落ちる」という別意味を作らない。
    if (exponent.isNumber() && exponent.asNumber().isReal()
        && exponent.asNumber().asReal().toRational()
            == Rational{BigInt{1}, BigInt{2}}) {
        return evaluateSqrt(std::span<const Expr>{&base, 1}, registry, mathematics);
    }

    if (!exponent.isNumber() || !exponent.asNumber().isReal()
        || !exponent.asNumber().asReal().isInteger())
        return Expr::call(registry.symbol(evaluation::BuiltinId::Power), {base, exponent});

    const BigInt& integerExponent = exponent.asNumber().asReal().asInteger();
    if (integerExponent.isZero()) {
        // 0^0をIndeterminateとしている以上、未知のsymbolic baseに対してx^0 -> 1 と無条件簡約するのは安全ではない。
        // baseが非零と証明できる場合だけ1へ畳み込み、未知ならPower式を保持する。
        if (base.isNumber())
            return integerExpr(1);

        const mathematics::ValueFacts baseFacts = mathematics::inferValueFacts(
            base, registry, mathematics);
        if (baseFacts.sign == mathematics::RealSign::Positive
            || baseFacts.sign == mathematics::RealSign::Negative
            || baseFacts.sign == mathematics::RealSign::NonZero)
            return integerExpr(1);

        return Expr::call(registry.symbol(evaluation::BuiltinId::Power), {base, exponent});
    }
    if (integerExponent == BigInt{1})
        return base;

    if (!base.isNumber()) {
        const auto algebraicBase = algebraicValue(base, registry);
        const auto magnitude = numeric::tryToUint64(integerExponent.abs());
        if (algebraicBase && magnitude && *magnitude <= 32) {
            if (auto fieldPower = symbolic::AlgebraicNumber::integerPowerInField(
                    *algebraicBase, *magnitude, integerExponent.isNegative()))
                return algebraicExpr(*fieldPower, registry);

            auto result = symbolic::AlgebraicNumber::fromRational(Rational{BigInt{1}});
            auto factor = algebraicBase;
            std::uint64_t power = *magnitude;
            while (result && factor && power != 0) {
                if ((power & 1U) != 0)
                    result = symbolic::AlgebraicNumber::combine(
                        *result, *factor, symbolic::AlgebraicBinaryOperation::Multiply);
                power >>= 1U;
                if (power != 0)
                    factor = symbolic::AlgebraicNumber::combine(
                        *factor, *factor, symbolic::AlgebraicBinaryOperation::Multiply);
            }
            // 中間代数演算がbudget/証明不能で失敗した場合，部分結果を確定値として返さない。
            if (power != 0)
                result = std::nullopt;
            if (result) {
                if (integerExponent.isNegative()) {
                    const auto one = symbolic::AlgebraicNumber::fromRational(Rational{BigInt{1}});
                    result = one ? symbolic::AlgebraicNumber::combine(
                        *one, *result, symbolic::AlgebraicBinaryOperation::Divide) : std::nullopt;
                }
                if (result)
                    return algebraicExpr(*result, registry);
            }
        }
        return Expr::call(registry.symbol(evaluation::BuiltinId::Power), {base, exponent});
    }

    const bool negativeExponent = integerExponent.isNegative();
    const auto magnitude = toUint64(integerExponent.abs());
    if (!magnitude)
        error::throwCalcError(
            error::CalcErrorType::Overflow,
            "Exponent is too large for exact evaluation");

    if (negativeExponent && base.asNumber().isZero())
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "Zero cannot be raised to a negative power");

    // 巨大値を構築してから検査するのでは遅い。既約なexact realの整数冪では，
    // numerator/denominatorの少なくとも一方が概ね exponent*(bits-1) bitになる。
    if (base.asNumber().isReal()) {
        if (evaluation::EvaluationBudget* budget = evaluation::currentEvaluationBudget()) {
            const std::size_t bits = realNumberBits(base.asNumber().asReal());
            const std::size_t limit = budget->limits().maxBigIntegerBits;
            if (bits > 1 && *magnitude > (limit == 0 ? 0 : (limit - 1) / (bits - 1)))
                budget->checkBigIntegerBits(limit + (limit != std::numeric_limits<std::size_t>::max()));
        }
    }

    Number result = numeric::integerPower(base.asNumber(), *magnitude);
    if (negativeExponent)
        result = Number{BigInt{1}} / result;

    return numberExpr(std::move(result));
}

Expr evaluateNegate(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    requireArity(arguments, 1, names::negate);
    if (arguments.front().isArray()) {
        const auto& array = arguments.front().asArray();
        if (array.hasExactNumberStorage()) {
            std::vector<Number> values;
            values.reserve(array.size());
            for (std::size_t i = 0; i < array.size(); ++i)
                values.push_back(-array.exactNumber(i));
            return numberArray(array.shape, std::move(values));
        }
        std::vector<Expr> elements;
        elements.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i) {
            elements.push_back(scalarNegate(array.element(i), registry));
        }
        return Expr::array(array.shape, std::move(elements));
    }
    return scalarNegate(arguments.front(), registry);
}

Expr evaluateFactorial(std::span<const Expr> arguments) {
    requireArity(arguments, 1, names::factorial);

    const Expr& argument = arguments.front();
    if (!argument.isNumber() || !argument.asNumber().isReal()
        || !argument.asNumber().asReal().isInteger())
        error::throwCalcError(
            error::CalcErrorType::Type,
            "Factorial requires a non-negative integer");

    const BigInt& value = argument.asNumber().asReal().asInteger();
    if (value.isNegative())
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "Factorial is undefined for negative integers");

    const auto count = toUint64(value);
    if (!count)
        error::throwCalcError(
            error::CalcErrorType::Overflow,
            "Factorial argument is too large for exact evaluation");

    if (evaluation::EvaluationBudget* budget = evaluation::currentEvaluationBudget()) {
        const std::size_t limit = budget->limits().maxBigIntegerBits;
        if (*count >= 4) {
            const std::uint64_t half = *count / 2;
            const std::size_t lowerBitsPerFactor = BigInt::fromUnsigned(half).bitLength() - 1;
            if (lowerBitsPerFactor != 0
                && half > limit / lowerBitsPerFactor)
                budget->checkBigIntegerBits(limit + (limit != std::numeric_limits<std::size_t>::max()));
        }
    }

    return numberExpr(Number{numeric::factorial(*count)});
}

Expr evaluateSqrt(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    requireArity(arguments, 1, names::sqrt);
    static_cast<void>(mathematics);

    // principal sqrtのexact/conditional簡約は共通Simplifierに集約する。
    return Expr::call(
        registry.symbol(evaluation::BuiltinId::Sqrt),
        {arguments.front()});
}

} // namespace mmcal::builtins
