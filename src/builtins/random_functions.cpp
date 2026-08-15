// 乱数組込み函数
#include "random_functions.hpp"

#include "builtins/exact_operations.hpp"
#include "error/error_message.hpp"
#include "numeric/number.hpp"
#include "numeric/real_number.hpp"
#include "numeric/integer_algorithms.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] Expr number(RealNumber value) {
    return Expr{Number{std::move(value)}};
}

[[nodiscard]] Expr integer(BigInt value) {
    return number(RealNumber{std::move(value)});
}

[[nodiscard]] const RealNumber* exactReal(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return nullptr;
    return &expression.asNumber().asReal();
}

[[nodiscard]] const BigInt* exactInteger(const Expr& expression) {
    const RealNumber* real = exactReal(expression);
    if (!real || !real->isInteger())
        return nullptr;
    return &real->asInteger();
}

[[nodiscard]] Expr unitSample(random::RandomEngine& engine) {
    return number(RealNumber{engine.unit53()});
}

[[nodiscard]] BigInt inclusiveIntegerSample(
    const BigInt& lower,
    const BigInt& upper,
    random::RandomEngine& engine) {
    if (lower > upper)
        error::throwCalcError(error::CalcErrorType::Domain, "randint requires lower <= upper");

    const BigInt width = upper - lower + BigInt{1};
    return lower + engine.uniformBelow(width);
}

} // namespace

Expr evaluateRandSeed(
    std::span<const Expr> arguments,
    random::RandomEngine& engine) {
    if (arguments.empty())
        return integer(engine.reseedFromEntropy());

    const BigInt* seedValue = exactInteger(arguments.front());
    if (!seedValue)
        error::throwCalcError(error::CalcErrorType::Type, "randSeed requires an integer seed");

    engine.seed(*seedValue);
    return integer(*seedValue);
}

Expr evaluateRand(
    std::span<const Expr> arguments,
    random::RandomEngine& engine) {
    if (arguments.empty())
        return unitSample(engine);

    if (arguments.size() == 1) {
        const RealNumber* upper = exactReal(arguments.front());
        if (!upper)
            error::throwCalcError(error::CalcErrorType::Type, "rand requires exact real bounds");
        if (upper->isNegative())
            error::throwCalcError(error::CalcErrorType::Domain, "rand upper bound must be non-negative");
        if (upper->isZero())
            return number(*upper);

        RealNumber value{engine.unit53()};
        value *= *upper;
        return number(std::move(value));
    }

    const RealNumber* lower = exactReal(arguments[0]);
    const RealNumber* upper = exactReal(arguments[1]);
    if (!lower || !upper)
        error::throwCalcError(error::CalcErrorType::Type, "rand requires exact real bounds");
    if (*lower > *upper)
        error::throwCalcError(error::CalcErrorType::Domain, "rand requires lower <= upper");
    if (*lower == *upper)
        return number(*lower);

    RealNumber width = *upper - *lower;
    width *= RealNumber{engine.unit53()};
    width += *lower;
    return number(std::move(width));
}

Expr evaluateRandInt(
    std::span<const Expr> arguments,
    random::RandomEngine& engine) {
    if (arguments.empty())
        return integer(engine.uniformBelow(BigInt{2}));

    if (arguments.size() == 1) {
        const BigInt* bound = exactInteger(arguments.front());
        if (!bound)
            error::throwCalcError(error::CalcErrorType::Type, "randint requires integer bounds");
        if (bound->isNegative())
            return integer(inclusiveIntegerSample(*bound, BigInt{}, engine));
        return integer(inclusiveIntegerSample(BigInt{}, *bound, engine));
    }

    const BigInt* lower = exactInteger(arguments[0]);
    const BigInt* upper = exactInteger(arguments[1]);
    if (!lower || !upper)
        error::throwCalcError(error::CalcErrorType::Type, "randint requires integer bounds");
    return integer(inclusiveIntegerSample(*lower, *upper, engine));
}

Expr evaluateChoice(
    std::span<const Expr> arguments,
    random::RandomEngine& engine) {
    if (arguments.size() == 1 && arguments.front().isArray()) {
        const auto& array = arguments.front().asArray();
        if (array.rank() != 1)
            error::throwCalcError(error::CalcErrorType::Type, "choice array must have rank 1");
        if (array.empty())
            error::throwCalcError(error::CalcErrorType::Domain, "choice requires a non-empty array");
        const BigInt index = engine.uniformBelow(BigInt::parse(std::to_string(array.size())));
        const auto converted = numeric::tryToUint64(index);
        return array.element(static_cast<std::size_t>(*converted));
    }

    if (arguments.empty())
        error::throwCalcError(error::CalcErrorType::Domain, "choice requires at least one value");

    const BigInt index = engine.uniformBelow(BigInt::parse(std::to_string(arguments.size())));
    const auto converted = numeric::tryToUint64(index);
    return arguments[static_cast<std::size_t>(*converted)];
}

Expr evaluateRandN(
    std::span<const Expr> arguments,
    random::RandomEngine& engine,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    RealNumber mean{BigInt{}};
    RealNumber sigma{BigInt{1}};

    if (!arguments.empty()) {
        const RealNumber* value = exactReal(arguments[0]);
        if (!value)
            error::throwCalcError(error::CalcErrorType::Type, "randn requires exact real parameters");
        mean = *value;
    }
    if (arguments.size() == 2) {
        const RealNumber* value = exactReal(arguments[1]);
        if (!value)
            error::throwCalcError(error::CalcErrorType::Type, "randn requires exact real parameters");
        if (value->isNegative())
            error::throwCalcError(error::CalcErrorType::Domain, "randn sigma must be non-negative");
        sigma = *value;
    }
    if (sigma.isZero())
        return number(std::move(mean));

    // Box-Mullerをexactな53bit dyadic uniform sampleへ適用する。u1は(0,1]としてlog[0]を避け、角度はsession既定に依存しない明示Radで構成する。
    const Expr u1 = number(RealNumber{engine.positiveUnit53()});
    const Expr u2 = number(RealNumber{engine.unit53()});
    const Expr two{Number{BigInt{2}}};
    const Expr minusTwo{Number{BigInt{-2}}};
    const Expr logU1 = Expr::call(registry.symbol(evaluation::BuiltinId::Log), {u1});
    const Expr radius = exact::sqrt(
        exact::multiply({minusTwo, logU1}, registry, mathematics, angles),
        registry, mathematics, angles);

    const auto* piDefinition = mathematics.findConstant(mathematics::ConstantId::Pi);
    if (!piDefinition)
        error::throwCalcError(error::CalcErrorType::Internal, "Pi is not registered");
    const Expr piSymbol{piDefinition->symbol};
    const Expr theta = exact::multiply({two, piSymbol, u2}, registry, mathematics, angles);
    const Expr explicitRad = Expr::call(
        registry.symbol(evaluation::BuiltinId::UnitApplied),
        {theta, Expr{std::string{"Rad"}}});
    const Expr cosine = Expr::call(registry.symbol(evaluation::BuiltinId::Cos), {explicitRad});
    Expr sample = exact::multiply({radius, cosine}, registry, mathematics, angles);

    sample = exact::multiply({number(sigma), std::move(sample)}, registry, mathematics, angles);
    return exact::add({number(std::move(mean)), std::move(sample)}, registry, mathematics, angles);
}

} // namespace mmcal::builtins
