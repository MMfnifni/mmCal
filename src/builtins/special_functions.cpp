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


[[nodiscard]] bool isNonPositiveInteger(const Rational& value) {
    return value.isInteger() && !value.numerator().isPositive();
}

[[nodiscard]] Expr evaluateHypergeometric1F1(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 3, names::hypergeometric1F1);

    Rational a;
    Rational b;
    Rational z;
    if (!exactRealRational(arguments[0], a)
        || !exactRealRational(arguments[1], b)
        || !exactRealRational(arguments[2], z))
        return hold(BuiltinId::Hypergeometric1F1, arguments, registry);

    // b=0,-1,-2,... は通常parameter pole。ただし a=-m で級数がpole到達前に
    // terminateする場合だけ有限多項式として安全に評価できる。
    std::optional<std::uint64_t> terminatingOrder;
    if (a.isInteger() && !a.numerator().isPositive()) {
        const auto order = numeric::tryToUint64(-a.numerator());
        if (order && *order <= 4096)
            terminatingOrder = *order;
    }

    if (isNonPositiveInteger(b)) {
        if (!terminatingOrder)
            return hold(BuiltinId::Hypergeometric1F1, arguments, registry);
        const auto poleIndex = numeric::tryToUint64(-b.numerator());
        if (!poleIndex || *terminatingOrder > *poleIndex)
            return hold(BuiltinId::Hypergeometric1F1, arguments, registry);
    }

    if (z.isZero())
        return integer(1);

    // M(a,a,z)=exp(z)。parameter poleを跨がない場合のみdefinednessを保って簡約する。
    if (a == b && !isNonPositiveInteger(b))
        return exact::call(BuiltinId::Exp, {arguments[2]}, registry, mathematics, angles);

    // a=-m は有限級数なので、exact Rational入力なら完全にexact評価する。
    if (terminatingOrder) {
        Rational term{BigInt{1}};
        Rational sum{BigInt{1}};
        for (std::uint64_t k = 0; k < *terminatingOrder; ++k) {
            const Rational numeratorFactor = a + Rational{BigInt::fromUnsigned(k)};
            const Rational denominatorFactor = b + Rational{BigInt::fromUnsigned(k)};
            if (denominatorFactor.isZero())
                return hold(BuiltinId::Hypergeometric1F1, arguments, registry);
            term *= numeratorFactor * z;
            term /= denominatorFactor * Rational{BigInt::fromUnsigned(k + 1)};
            sum += term;
        }
        return rationalExpr(std::move(sum));
    }

    return hold(BuiltinId::Hypergeometric1F1, arguments, registry);
}


[[nodiscard]] std::optional<std::uint64_t> terminatingHypergeometricOrder(const Rational& value) {
    if (!value.isInteger() || value.numerator().isPositive())
        return std::nullopt;
    const auto order = numeric::tryToUint64(-value.numerator());
    if (!order || *order > 4096)
        return std::nullopt;
    return *order;
}

[[nodiscard]] Expr evaluateHypergeometric2F1(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 4, names::hypergeometric2F1);
    static_cast<void>(mathematics);
    static_cast<void>(angles);

    Rational a;
    Rational b;
    Rational c;
    Rational z;
    if (!exactRealRational(arguments[0], a)
        || !exactRealRational(arguments[1], b)
        || !exactRealRational(arguments[2], c)
        || !exactRealRational(arguments[3], z))
        return hold(BuiltinId::Hypergeometric2F1, arguments, registry);

    const auto orderA = terminatingHypergeometricOrder(a);
    const auto orderB = terminatingHypergeometricOrder(b);
    std::optional<std::uint64_t> terminatingOrder;
    if (orderA && orderB)
        terminatingOrder = std::min(*orderA, *orderB);
    else if (orderA)
        terminatingOrder = orderA;
    else if (orderB)
        terminatingOrder = orderB;

    // c=0,-1,-2,... は通常parameter pole。有限級数がpole到達前に停止する場合だけ
    // exact polynomialとして安全に受理する。1F1と同じdefinedness方針。
    if (isNonPositiveInteger(c)) {
        if (!terminatingOrder)
            return hold(BuiltinId::Hypergeometric2F1, arguments, registry);
        const auto poleIndex = numeric::tryToUint64(-c.numerator());
        if (!poleIndex || *terminatingOrder > *poleIndex)
            return hold(BuiltinId::Hypergeometric2F1, arguments, registry);
    }

    if (z.isZero())
        return integer(1);

    if (terminatingOrder) {
        Rational term{BigInt{1}};
        Rational sum{BigInt{1}};
        for (std::uint64_t k = 0; k < *terminatingOrder; ++k) {
            const Rational ka{BigInt::fromUnsigned(k)};
            const Rational denominatorFactor = c + ka;
            if (denominatorFactor.isZero())
                return hold(BuiltinId::Hypergeometric2F1, arguments, registry);
            term *= (a + ka) * (b + ka) * z;
            term /= denominatorFactor * Rational{BigInt::fromUnsigned(k + 1)};
            sum += term;
        }
        return rationalExpr(std::move(sum));
    }

    return hold(BuiltinId::Hypergeometric2F1, arguments, registry);
}

[[nodiscard]] bool exactZero(const Expr& expression) {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] Expr evaluateEllipticF(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 2, names::ellipticF);
    if (exactZero(arguments[0]))
        return integer(0);
    // F(phi|0)=phi。phiがsymbolicでも成立するためSolverにも安全に還元できる。
    if (exactZero(arguments[1]))
        return arguments[0];

    Rational phi;
    if (exactRealRational(arguments[0], phi) && phi.numerator().isNegative())
        return exact::negate(
            Expr::call(registry.symbol(BuiltinId::EllipticF), {
                rationalExpr(-phi), arguments[1]}), registry, mathematics, angles);
    return hold(BuiltinId::EllipticF, arguments, registry);
}

[[nodiscard]] Expr evaluateEllipticE(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 2, names::ellipticE);
    if (exactZero(arguments[0]))
        return integer(0);
    if (exactZero(arguments[1]))
        return arguments[0];

    Rational phi;
    if (exactRealRational(arguments[0], phi) && phi.numerator().isNegative())
        return exact::negate(
            Expr::call(registry.symbol(BuiltinId::EllipticE), {
                rationalExpr(-phi), arguments[1]}), registry, mathematics, angles);
    return hold(BuiltinId::EllipticE, arguments, registry);
}

[[nodiscard]] Expr evaluateEllipticPi(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 3, names::ellipticPi);
    if (exactZero(arguments[1]))
        return integer(0);
    // Pi(0;phi|m)=F(phi|m)。特殊値を既存函数へ落とすことでD/N/Solveの知識も共有する。
    if (exactZero(arguments[0])) {
        const std::array<Expr, 2> fArguments{arguments[1], arguments[2]};
        return evaluateEllipticF(fArguments, registry, mathematics, angles);
    }

    Rational phi;
    if (exactRealRational(arguments[1], phi) && phi.numerator().isNegative())
        return exact::negate(
            Expr::call(registry.symbol(BuiltinId::EllipticPi), {
                arguments[0], rationalExpr(-phi), arguments[2]}),
            registry, mathematics, angles);
    return hold(BuiltinId::EllipticPi, arguments, registry);
}

[[nodiscard]] Expr evaluateClassicalIntegralFunction(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::string_view name;
    switch (id) {
    case BuiltinId::ExponentialIntegralEi: name = names::exponentialIntegralEi; break;
    case BuiltinId::SineIntegralSi: name = names::sineIntegralSi; break;
    case BuiltinId::CosineIntegralCi: name = names::cosineIntegralCi; break;
    case BuiltinId::LogarithmicIntegralLi: name = names::logarithmicIntegralLi; break;
    default:
        error::throwCalcError(error::CalcErrorType::Internal, "Unexpected classical integral function");
    }
    requireArity(arguments, 1, name);

    Rational value;
    if (!exactRealRational(arguments.front(), value))
        return hold(id, arguments, registry);

    if ((id == BuiltinId::ExponentialIntegralEi || id == BuiltinId::CosineIntegralCi)
        && value.isZero())
        error::throwCalcError(error::CalcErrorType::Domain,
            std::string{name} + " is undefined at zero");
    if (id == BuiltinId::LogarithmicIntegralLi && value == Rational{BigInt{1}})
        error::throwCalcError(error::CalcErrorType::Domain, "li is undefined at one");

    // Siはentireな奇函数。Ei/Ci/liはprincipal branchを持つため、同じ反射を一般化しない。
    if (id == BuiltinId::SineIntegralSi) {
        if (value.isZero())
            return integer(0);
        if (value.numerator().isNegative())
            return exact::negate(
                Expr::call(registry.symbol(id), {rationalExpr(-value)}),
                registry, mathematics, angles);
    }
    return hold(id, arguments, registry);
}

[[nodiscard]] Expr evaluatePolylog(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 2, names::polylog);

    if (exactZero(arguments[1]))
        return integer(0);

    Rational order;
    if (!exactRealRational(arguments[0], order) || !order.isInteger())
        return hold(BuiltinId::Polylog, arguments, registry);

    if (order.numerator() == BigInt{0}) {
        // Li_0(z)=z/(1-z)。z=1ではpoleなのでexact divisionのDomainErrorをそのまま使う。
        return exact::divide(
            arguments[1],
            exact::subtract(integer(1), arguments[1], registry, mathematics, angles),
            registry, mathematics, angles);
    }
    if (order.numerator() == BigInt{1}) {
        // Li_1(z)=-Log(1-z)。zがsymbolicでも成立するprincipal-branch identity。
        Expr oneMinusZ = exact::subtract(integer(1), arguments[1], registry, mathematics, angles);
        return exact::negate(
            exact::call(BuiltinId::Log, {std::move(oneMinusZ)}, registry, mathematics, angles),
            registry, mathematics, angles);
    }

    Rational z;
    if (!exactRealRational(arguments[1], z))
        return hold(BuiltinId::Polylog, arguments, registry);

    if (order.numerator() == BigInt{2} && (z == Rational{BigInt{1}} || z == Rational{BigInt{-1}})) {
        Expr piSquared = exact::call(
            BuiltinId::Power, {piExpr(mathematics), integer(2)}, registry, mathematics, angles);
        return exact::divide(
            z.numerator().isNegative()
                ? exact::negate(std::move(piSquared), registry, mathematics, angles)
                : std::move(piSquared),
            integer(z.numerator().isNegative() ? 12 : 6),
            registry, mathematics, angles);
    }

    return hold(BuiltinId::Polylog, arguments, registry);
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
    case BuiltinId::Hypergeometric1F1:
        return evaluateHypergeometric1F1(arguments, registry, mathematics, angles);
    case BuiltinId::Hypergeometric2F1:
        return evaluateHypergeometric2F1(arguments, registry, mathematics, angles);
    case BuiltinId::EllipticF:
        return evaluateEllipticF(arguments, registry, mathematics, angles);
    case BuiltinId::EllipticE:
        return evaluateEllipticE(arguments, registry, mathematics, angles);
    case BuiltinId::EllipticPi:
        return evaluateEllipticPi(arguments, registry, mathematics, angles);
    case BuiltinId::ExponentialIntegralEi:
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
    case BuiltinId::LogarithmicIntegralLi:
        return evaluateClassicalIntegralFunction(id, arguments, registry, mathematics, angles);
    case BuiltinId::Polylog:
        return evaluatePolylog(arguments, registry, mathematics, angles);
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
