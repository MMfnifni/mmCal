// 三角函数の整数冪・積を有限Fourier和へ落とす厳密恒等式
#include "trigonometric_polynomial.hpp"

#include "evaluation/builtin_registry.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::mathematics {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return expression.isCall()
        && expression.asCall().head.sameIdentity(builtins.symbol(id));
}

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr rational(Rational value) {
    return Expr{Number{std::move(value)}};
}

[[nodiscard]] Expr call(
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id,
    std::vector<Expr> arguments) {
    return Expr::call(builtins.symbol(id), std::move(arguments));
}

[[nodiscard]] std::optional<std::size_t> nonnegativeSmallInteger(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    const Rational value = expression.asNumber().asReal().toRational();
    if (!value.isInteger() || value.numerator().isNegative())
        return std::nullopt;
    const auto converted = numeric::tryToUint64(value.numerator());
    if (!converted || *converted > static_cast<std::uint64_t>(SIZE_MAX))
        return std::nullopt;
    return static_cast<std::size_t>(*converted);
}

struct TrigFactor final {
    BuiltinId id = BuiltinId::Sin;
    Expr argument;
    std::size_t exponent = 0;
};

[[nodiscard]] std::optional<TrigFactor> trigFactor(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    const Expr* base = &expression;
    std::size_t exponent = 1;
    if (isHead(expression, builtins, BuiltinId::Power)) {
        const auto& a = expression.asCall().arguments;
        if (a.size() != 2)
            return std::nullopt;
        const auto parsed = nonnegativeSmallInteger(a[1]);
        if (!parsed)
            return std::nullopt;
        exponent = *parsed;
        base = &a[0];
    }

    if (!base->isCall() || base->asCall().arguments.size() != 1)
        return std::nullopt;
    if (isHead(*base, builtins, BuiltinId::Sin))
        return TrigFactor{BuiltinId::Sin, base->asCall().arguments.front(), exponent};
    if (isHead(*base, builtins, BuiltinId::Cos))
        return TrigFactor{BuiltinId::Cos, base->asCall().arguments.front(), exponent};
    return std::nullopt;
}

[[nodiscard]] Expr scaledAngleArgument(
    const Expr& source,
    std::size_t multiplier,
    const evaluation::BuiltinRegistry& builtins) {
    if (multiplier == 1)
        return source;

    const Expr factor{Number{BigInt::fromUnsigned(static_cast<std::uint64_t>(multiplier))}};
    if (isHead(source, builtins, BuiltinId::UnitApplied)) {
        const auto& a = source.asCall().arguments;
        if (a.size() == 2)
            return call(builtins, BuiltinId::UnitApplied, {
                call(builtins, BuiltinId::Multiply, {factor, a[0]}), a[1]});
    }
    return call(builtins, BuiltinId::Multiply, {factor, source});
}

[[nodiscard]] Expr scaledTerm(
    Rational coefficient,
    Expr atom,
    const evaluation::BuiltinRegistry& builtins) {
    if (coefficient == Rational{BigInt{1}})
        return atom;
    if (coefficient == Rational{BigInt{-1}})
        return call(builtins, BuiltinId::Negate, {std::move(atom)});
    return call(builtins, BuiltinId::Multiply, {
        rational(std::move(coefficient)), std::move(atom)});
}

[[nodiscard]] std::optional<std::pair<Expr, Expr>> compatibleAnglePair(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins) {
    const bool lhsUnit = isHead(lhs, builtins, BuiltinId::UnitApplied);
    const bool rhsUnit = isHead(rhs, builtins, BuiltinId::UnitApplied);
    if (lhsUnit != rhsUnit)
        return std::nullopt;
    if (!lhsUnit)
        return std::pair<Expr, Expr>{lhs, rhs};

    const auto& la = lhs.asCall().arguments;
    const auto& ra = rhs.asCall().arguments;
    if (la.size() != 2 || ra.size() != 2 || !(la[1] == ra[1]))
        return std::nullopt;
    return std::pair<Expr, Expr>{la[0], ra[0]};
}

[[nodiscard]] Expr combineAngleArguments(
    const Expr& lhs,
    const Expr& rhs,
    bool subtract,
    const Expr& originalLhs,
    const evaluation::BuiltinRegistry& builtins) {
    Expr combined = call(builtins, subtract ? BuiltinId::Subtract : BuiltinId::Add, {lhs, rhs});
    if (!isHead(originalLhs, builtins, BuiltinId::UnitApplied))
        return combined;
    return call(builtins, BuiltinId::UnitApplied, {
        std::move(combined), originalLhs.asCall().arguments[1]});
}

} // namespace

Expr scaledTrigArgumentForFrequency(
    const Expr& source,
    std::size_t multiplier,
    const evaluation::BuiltinRegistry& builtins) {
    return scaledAngleArgument(source, multiplier, builtins);
}

std::optional<TrigFourierExpansion> expandTrigMonomial(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    std::size_t maximumTotalDegree) {
    std::vector<TrigFactor> factors;
    if (isHead(expression, builtins, BuiltinId::Multiply)) {
        for (const Expr& factor : expression.asCall().arguments) {
            auto parsed = trigFactor(factor, builtins);
            if (!parsed)
                return std::nullopt;
            factors.push_back(std::move(*parsed));
        }
    }
    else {
        auto parsed = trigFactor(expression, builtins);
        if (!parsed)
            return std::nullopt;
        factors.push_back(std::move(*parsed));
    }

    if (factors.empty())
        return std::nullopt;

    Expr argument = factors.front().argument;
    std::size_t sinePower = 0;
    std::size_t cosinePower = 0;
    for (const TrigFactor& factor : factors) {
        if (!(factor.argument == argument))
            return std::nullopt;
        if (factor.exponent > maximumTotalDegree)
            return std::nullopt;
        if (factor.id == BuiltinId::Sin)
            sinePower += factor.exponent;
        else
            cosinePower += factor.exponent;
        if (sinePower > maximumTotalDegree || cosinePower > maximumTotalDegree
            || sinePower + cosinePower > maximumTotalDegree)
            return std::nullopt;
    }

    const std::size_t total = sinePower + cosinePower;
    if (total <= 1)
        return std::nullopt;

    // (y-y^-1)^m (y+y^-1)^n の整数係数を直接畳み込む。
    // 旧来の個別sin^2/cos^2規則を高次数へコピーせず，全整数冪を同じ恒等式から生成する。
    std::vector<BigInt> coefficients{BigInt{1}};
    for (std::size_t factorIndex = 0; factorIndex < total; ++factorIndex) {
        const bool sineFactor = factorIndex < sinePower;
        std::vector<BigInt> next(coefficients.size() + 1, BigInt{});
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            next[i] += coefficients[i];
            next[i + 1] += sineFactor ? -coefficients[i] : coefficients[i];
        }
        coefficients = std::move(next);
    }

    BigInt denominator{1};
    denominator <<= total;
    const BigInt pairDenominator = denominator >> 1U;
    std::vector<TrigFourierTerm> terms;

    if ((sinePower & 1U) == 0) {
        const bool negativePhase = ((sinePower / 2U) & 1U) != 0;
        for (std::size_t j = 0; j * 2U < total; ++j) {
            BigInt numerator = coefficients[j];
            if (negativePhase)
                numerator = -numerator;
            if (numerator.isZero())
                continue;
            terms.push_back(TrigFourierTerm{
                Rational{std::move(numerator), pairDenominator}, total - 2U * j, false});
        }
        if ((total & 1U) == 0) {
            BigInt numerator = coefficients[total / 2U];
            if (negativePhase)
                numerator = -numerator;
            if (!numerator.isZero())
                terms.push_back(TrigFourierTerm{
                    Rational{std::move(numerator), denominator}, 0, false});
        }
    }
    else {
        const bool negativePhase = (((sinePower - 1U) / 2U) & 1U) != 0;
        for (std::size_t j = 0; j * 2U < total; ++j) {
            BigInt numerator = coefficients[j];
            if (negativePhase)
                numerator = -numerator;
            if (numerator.isZero())
                continue;
            terms.push_back(TrigFourierTerm{
                Rational{std::move(numerator), pairDenominator}, total - 2U * j, true});
        }
    }

    return TrigFourierExpansion{std::move(argument), std::move(terms), total};
}

std::optional<Expr> reduceTrigMonomial(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    std::size_t maximumTotalDegree) {
    auto expansion = expandTrigMonomial(expression, builtins, maximumTotalDegree);
    if (!expansion)
        return std::nullopt;

    std::vector<Expr> terms;
    terms.reserve(expansion->terms.size());
    for (TrigFourierTerm& term : expansion->terms) {
        if (term.frequency == 0) {
            terms.push_back(rational(std::move(term.coefficient)));
            continue;
        }
        Expr atom = call(builtins, term.sine ? BuiltinId::Sin : BuiltinId::Cos, {
            scaledTrigArgumentForFrequency(expansion->argument, term.frequency, builtins)});
        terms.push_back(scaledTerm(std::move(term.coefficient), std::move(atom), builtins));
    }

    if (terms.empty())
        return integer(0);
    if (terms.size() == 1)
        return terms.front();
    return call(builtins, BuiltinId::Add, std::move(terms));
}

std::optional<Expr> reduceTrigProduct(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::Multiply)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;

    const Expr& lhs = expression.asCall().arguments[0];
    const Expr& rhs = expression.asCall().arguments[1];
    auto left = trigFactor(lhs, builtins);
    auto right = trigFactor(rhs, builtins);
    if (!left || !right || left->exponent != 1 || right->exponent != 1)
        return std::nullopt;

    const auto pair = compatibleAnglePair(left->argument, right->argument, builtins);
    if (!pair)
        return std::nullopt;

    const Expr sum = combineAngleArguments(
        pair->first, pair->second, false, left->argument, builtins);
    const Expr difference = combineAngleArguments(
        pair->first, pair->second, true, left->argument, builtins);
    const Expr two = integer(2);

    if (left->id == BuiltinId::Cos && right->id == BuiltinId::Cos) {
        return call(builtins, BuiltinId::Divide, {
            call(builtins, BuiltinId::Add, {
                call(builtins, BuiltinId::Cos, {difference}),
                call(builtins, BuiltinId::Cos, {sum})}), two});
    }
    if (left->id == BuiltinId::Sin && right->id == BuiltinId::Sin) {
        return call(builtins, BuiltinId::Divide, {
            call(builtins, BuiltinId::Subtract, {
                call(builtins, BuiltinId::Cos, {difference}),
                call(builtins, BuiltinId::Cos, {sum})}), two});
    }

    // sin[a]cos[b] = (sin[a+b] + sin[a-b])/2。
    // cos[a]sin[b] は引数を交換して同じ恒等式へ落とす。
    if (left->id == BuiltinId::Cos && right->id == BuiltinId::Sin)
        return reduceTrigProduct(
            call(builtins, BuiltinId::Multiply, {rhs, lhs}), builtins);

    return call(builtins, BuiltinId::Divide, {
        call(builtins, BuiltinId::Add, {
            call(builtins, BuiltinId::Sin, {sum}),
            call(builtins, BuiltinId::Sin, {difference})}), two});
}

} // namespace mmcal::mathematics
