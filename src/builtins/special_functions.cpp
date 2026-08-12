// Gamma・erf・Betaなどの特殊函数
#include "special_functions.hpp"

#include "builtins/exact_operations.hpp"
#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

void requireArity(std::span<const Expr> arguments, std::size_t expected, std::string_view name) {
    if (arguments.size() != expected)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects " + std::to_string(expected) + " argument(s)");
}

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr integer(BigInt value) {
    return Expr{Number{std::move(value)}};
}

[[nodiscard]] Expr rationalExpr(Rational value) {
    return Expr{Number{RealNumber{std::move(value)}}};
}

[[nodiscard]] Expr hold(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry) {
    return Expr::call(registry.symbol(id), {arguments.begin(), arguments.end()});
}

[[nodiscard]] const Rational* exactRealRational(const Expr& expression, Rational& storage) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return nullptr;
    storage = expression.asNumber().asReal().toRational();
    return &storage;
}


[[nodiscard]] std::optional<std::uint64_t> nonnegativeIntegerCount(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    const BigInt& value = expression.asNumber().asReal().asInteger();
    if (value.isNegative())
        return std::nullopt;
    return numeric::tryToUint64(value);
}

[[nodiscard]] Expr piExpr(const mathematics::MathRegistry& mathematics) {
    const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
    if (!pi)
        error::throwCalcError(error::CalcErrorType::Internal, "Pi is not registered");
    return Expr{pi->symbol};
}

[[nodiscard]] Expr gammaHalfInteger(
    const Rational& value,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    // denominator=2、numeratorは奇数を想定。
    const BigInt numerator = value.numerator();
    const bool positive = numerator.isPositive();

    BigInt nBig;
    Rational coefficient;
    if (positive) {
        // x=n+1/2: Gamma(x)=(2n)!/(4^n n!) sqrt(Pi)
        nBig = (numerator - BigInt{1}) / BigInt{2};
        const auto n = numeric::tryToUint64(nBig);
        if (!n || *n > std::numeric_limits<std::uint64_t>::max() / 2U
            || *n > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) / 2U)
            error::throwCalcError(error::CalcErrorType::Overflow,
                "gamma half-integer argument is too large for exact evaluation");
        const std::uint64_t twiceN = 2U * *n;
        BigInt fourPower{1};
        fourPower <<= static_cast<std::size_t>(twiceN);
        coefficient = Rational{numeric::factorial(twiceN), fourPower * numeric::factorial(*n)};
    }
    else {
        // x=1/2-n: Gamma(x)=(-4)^n n!/(2n)! sqrt(Pi)
        nBig = (BigInt{1} - numerator) / BigInt{2};
        const auto n = numeric::tryToUint64(nBig);
        if (!n || *n > std::numeric_limits<std::uint64_t>::max() / 2U
            || *n > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) / 2U)
            error::throwCalcError(error::CalcErrorType::Overflow,
                "gamma half-integer argument is too large for exact evaluation");
        const std::uint64_t twiceN = 2U * *n;
        BigInt fourPower{1};
        fourPower <<= static_cast<std::size_t>(twiceN);
        if ((*n & 1U) != 0)
            fourPower = -fourPower;
        coefficient = Rational{fourPower * numeric::factorial(*n), numeric::factorial(twiceN)};
    }

    Expr rootPi = exact::sqrt(piExpr(mathematics), registry, mathematics, angles);
    if (coefficient == Rational{BigInt{1}})
        return rootPi;
    return exact::multiply(
        {rationalExpr(std::move(coefficient)), std::move(rootPi)},
        registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateGamma(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::gamma);
    Rational value;
    if (!exactRealRational(arguments.front(), value))
        return hold(BuiltinId::Gamma, arguments, registry);

    if (value.isInteger()) {
        const BigInt& n = value.numerator();
        if (!n.isPositive())
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "gamma is undefined at non-positive integers");
        const BigInt nMinusOne = n - BigInt{1};
        const auto count = numeric::tryToUint64(nMinusOne);
        if (!count)
            error::throwCalcError(
                error::CalcErrorType::Overflow,
                "gamma integer argument is too large for exact evaluation");
        return integer(numeric::factorial(*count));
    }

    if (value.denominator() == BigInt{2})
        return gammaHalfInteger(value, registry, mathematics, angles);

    return hold(BuiltinId::Gamma, arguments, registry);
}

[[nodiscard]] Expr evaluateLogGamma(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::logGamma);
    Rational value;
    if (exactRealRational(arguments.front(), value) && value.isInteger()
        && !value.numerator().isPositive())
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "lgamma is undefined at non-positive integers");

    Expr gamma = evaluateGamma(arguments, registry, mathematics, angles);
    if (gamma.isCall() && gamma.asCall().head.sameIdentity(registry.symbol(BuiltinId::Gamma)))
        return hold(BuiltinId::LogGamma, arguments, registry);

    Expr magnitude = exact::call(BuiltinId::Abs, {std::move(gamma)}, registry, mathematics, angles);
    return exact::call(BuiltinId::Log, {std::move(magnitude)}, registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateErfLike(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, id == BuiltinId::Erf ? names::erf : names::erfc);
    Rational value;
    if (!exactRealRational(arguments.front(), value))
        return hold(id, arguments, registry);

    if (value.isZero())
        return integer(id == BuiltinId::Erf ? 0 : 1);

    if (value.numerator().isNegative()) {
        Expr positive = rationalExpr(-value);
        Expr reflected = Expr::call(registry.symbol(id), {std::move(positive)});
        if (id == BuiltinId::Erf)
            return exact::negate(std::move(reflected), registry, mathematics, angles);
        return exact::subtract(integer(2), std::move(reflected), registry, mathematics, angles);
    }
    return hold(id, arguments, registry);
}


[[nodiscard]] Expr evaluateFresnel(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::string_view name = id == BuiltinId::FresnelC ? names::fresnelC : names::fresnelS;
    requireArity(arguments, 1, name);

    Rational value;
    if (!exactRealRational(arguments.front(), value))
        return hold(id, arguments, registry);
    if (value.isZero())
        return integer(0);

    // C(-x)=-C(x), S(-x)=-S(x)。entireな奇函数なのでbranch条件なしで安全に使える。
    if (value.numerator().isNegative()) {
        Expr positive = rationalExpr(-value);
        return exact::negate(
            Expr::call(registry.symbol(id), {std::move(positive)}),
            registry, mathematics, angles);
    }
    return hold(id, arguments, registry);
}

[[nodiscard]] Expr evaluateBeta(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 2, names::beta);

    Rational a;
    Rational b;
    const Rational* exactA = exactRealRational(arguments[0], a);
    const Rational* exactB = exactRealRational(arguments[1], b);
    if (!exactA || !exactB)
        return hold(BuiltinId::Beta, arguments, registry);
    if (a <= Rational{BigInt{0}} || b <= Rational{BigInt{0}})
        error::throwCalcError(error::CalcErrorType::Domain,
            "beta currently requires positive real arguments");

    // B(x,1)=1/x, B(1,y)=1/y。正実数域ではpoleを持ち込まず安全にexact化できる。
    if (a == Rational{BigInt{1}})
        return rationalExpr(Rational{BigInt{1}} / b);
    if (b == Rational{BigInt{1}})
        return rationalExpr(Rational{BigInt{1}} / a);

    // 正整数はfactorial比を直接使う。
    const auto m = nonnegativeIntegerCount(arguments[0]);
    const auto n = nonnegativeIntegerCount(arguments[1]);
    if (m && n) {
        if (*m > std::numeric_limits<std::uint64_t>::max() - *n)
            error::throwCalcError(error::CalcErrorType::Overflow, "beta arguments are too large");
        const BigInt numerator = numeric::factorial(*m - 1) * numeric::factorial(*n - 1);
        const BigInt denominator = numeric::factorial(*m + *n - 1);
        return rationalExpr(Rational{numerator, denominator});
    }

    // 半整数等、Gamma(a), Gamma(b), Gamma(a+b)が全てexactに畳める場合だけGamma恒等式を使う。
    // 一般symbolic Betaを無条件にGamma比へ展開しない。
    const Expr sum = rationalExpr(a + b);
    const std::array<Expr, 1> aArg{arguments[0]};
    const std::array<Expr, 1> bArg{arguments[1]};
    const std::array<Expr, 1> sumArg{sum};
    Expr ga = evaluateGamma(aArg, registry, mathematics, angles);
    Expr gb = evaluateGamma(bArg, registry, mathematics, angles);
    Expr gs = evaluateGamma(sumArg, registry, mathematics, angles);
    const auto isHeldGamma = [&](const Expr& value) {
        return value.isCall()
            && value.asCall().head.sameIdentity(registry.symbol(BuiltinId::Gamma));
    };
    if (isHeldGamma(ga) || isHeldGamma(gb) || isHeldGamma(gs))
        return hold(BuiltinId::Beta, arguments, registry);

    return exact::divide(
        exact::multiply({std::move(ga), std::move(gb)}, registry, mathematics, angles),
        std::move(gs), registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateBetaLog(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 2, names::betaLog);
    Expr beta = evaluateBeta(arguments, registry, mathematics, angles);
    if (beta.isCall()
        && beta.asCall().head.sameIdentity(registry.symbol(BuiltinId::Beta)))
        return hold(BuiltinId::BetaLog, arguments, registry);
    return exact::call(BuiltinId::Log, {std::move(beta)}, registry, mathematics, angles);
}

[[nodiscard]] Expr finiteFactorialProduct(
    const Expr& x,
    std::uint64_t n,
    bool rising,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (n == 0)
        return integer(1);

    std::vector<Expr> factors;
    factors.reserve(static_cast<std::size_t>(n));
    for (std::uint64_t k = 0; k < n; ++k) {
        const Expr offset = integer(BigInt::parse(std::to_string(k)));
        factors.push_back(rising
            ? exact::add({x, offset}, registry, mathematics, angles)
            : exact::subtract(x, offset, registry, mathematics, angles));
    }
    return exact::multiply(std::move(factors), registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateFiniteFactorial(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::string_view name = id == BuiltinId::FallingFactorial
        ? names::fallingFactorial : names::risingFactorial;
    requireArity(arguments, 2, name);
    if (!arguments[1].isNumber())
        return hold(id, arguments, registry);
    const auto n = nonnegativeIntegerCount(arguments[1]);
    if (!n) {
        if (arguments[1].isNumber() && arguments[1].asNumber().isReal()
            && arguments[1].asNumber().asReal().isInteger())
            error::throwCalcError(error::CalcErrorType::Domain,
                std::string{name} + " requires a non-negative integer order");
        error::throwCalcError(error::CalcErrorType::Type,
            std::string{name} + " requires an integer order");
    }
    return finiteFactorialProduct(
        arguments[0], *n, id == BuiltinId::RisingFactorial,
        registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateGeneralizedBinomial(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 2, names::generalizedBinomial);
    if (!arguments[1].isNumber())
        return hold(BuiltinId::GeneralizedBinomial, arguments, registry);
    const auto n = nonnegativeIntegerCount(arguments[1]);
    if (!n)
        return hold(BuiltinId::GeneralizedBinomial, arguments, registry);

    Expr numerator = finiteFactorialProduct(
        arguments[0], *n, false, registry, mathematics, angles);
    return exact::divide(
        std::move(numerator), integer(numeric::factorial(*n)),
        registry, mathematics, angles);
}

} // namespace

Expr evaluateSpecialFunction(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (id) {
    case BuiltinId::Log2:
    case BuiltinId::Log10:
        requireArity(arguments, 1, id == BuiltinId::Log2 ? names::log2 : names::log10);
        return exact::call(
            BuiltinId::Log,
            {integer(id == BuiltinId::Log2 ? 2 : 10), arguments.front()},
            registry, mathematics, angles);
    case BuiltinId::Gamma:
        return evaluateGamma(arguments, registry, mathematics, angles);
    case BuiltinId::LogGamma:
        return evaluateLogGamma(arguments, registry, mathematics, angles);
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
        return evaluateErfLike(id, arguments, registry, mathematics, angles);
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS:
        return evaluateFresnel(id, arguments, registry, mathematics, angles);
    case BuiltinId::Beta:
        return evaluateBeta(arguments, registry, mathematics, angles);
    case BuiltinId::BetaLog:
        return evaluateBetaLog(arguments, registry, mathematics, angles);
    case BuiltinId::GeneralizedBinomial:
        return evaluateGeneralizedBinomial(arguments, registry, mathematics, angles);
    case BuiltinId::FallingFactorial:
    case BuiltinId::RisingFactorial:
        return evaluateFiniteFactorial(id, arguments, registry, mathematics, angles);
    default:
        error::throwCalcError(error::CalcErrorType::Internal, "Unexpected special-function builtin");
    }
}

} // namespace mmcal::builtins
