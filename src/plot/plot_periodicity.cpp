#include "plot_periodicity.hpp"

#include "numeric/integer_algorithms.hpp"
#include "simplification/simplifier.hpp"
#include "symbolic/polynomial.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::plot {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Rational;

struct PeriodKnowledge final {
    bool dependsOnParameter = false;
    std::optional<Rational> periodTurns;
    PeriodProofKind proofKind = PeriodProofKind::Composition;
};

[[nodiscard]] Rational positive(Rational value) {
    if (value < Rational{BigInt{0}})
        return -value;
    return value;
}

[[nodiscard]] Rational rationalLcm(const Rational& lhs, const Rational& rhs) {
    const Rational a = positive(lhs);
    const Rational b = positive(rhs);
    return Rational{
        numeric::lcm(a.numerator(), b.numerator()),
        numeric::gcd(a.denominator(), b.denominator())};
}

// 正の有理周波数f_iに対し，f_i/gが全て整数となる最大の有理g。
[[nodiscard]] Rational rationalGcd(const Rational& lhs, const Rational& rhs) {
    const Rational a = positive(lhs);
    const Rational b = positive(rhs);
    return Rational{
        numeric::gcd(a.numerator(), b.numerator()),
        numeric::lcm(a.denominator(), b.denominator())};
}

[[nodiscard]] std::optional<Rational> commonPeriod(
    const std::optional<Rational>& lhs,
    const std::optional<Rational>& rhs) {
    if (!lhs)
        return rhs;
    if (!rhs)
        return lhs;
    return rationalLcm(*lhs, *rhs);
}

[[nodiscard]] std::optional<Rational> exactRealRational(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    return expression.asNumber().asReal().toRational();
}

[[nodiscard]] std::optional<Rational> affineCoefficient(
    const Expr& expression,
    const expression::Symbol& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto polynomial = symbolic::toExpressionPolynomial(
        expression, parameter, builtins, mathematics, angles);
    if (!polynomial || polynomial->degree() != 1)
        return std::nullopt;
    auto coefficient = exactRealRational(polynomial->coefficient(1));
    if (!coefficient || coefficient->isZero())
        return std::nullopt;
    return coefficient;
}

[[nodiscard]] std::optional<std::size_t> nonnegativeSmallInteger(const Expr& expression) {
    const auto rational = exactRealRational(expression);
    if (!rational || !rational->isInteger() || rational->numerator().isNegative())
        return std::nullopt;
    const auto value = numeric::tryToUint64(rational->numerator());
    if (!value || *value > static_cast<std::uint64_t>(SIZE_MAX))
        return std::nullopt;
    return static_cast<std::size_t>(*value);
}

void normalizeFrequencies(std::vector<Rational>& frequencies) {
    for (auto& frequency : frequencies)
        frequency = positive(std::move(frequency));
    std::sort(frequencies.begin(), frequencies.end());
    frequencies.erase(std::unique(frequencies.begin(), frequencies.end()), frequencies.end());
}

[[nodiscard]] std::optional<HarmonicSpectrum> mergeSpectra(
    HarmonicSpectrum lhs,
    const HarmonicSpectrum& rhs,
    std::size_t maximumFrequencies) {
    lhs.frequencies.insert(lhs.frequencies.end(), rhs.frequencies.begin(), rhs.frequencies.end());
    normalizeFrequencies(lhs.frequencies);
    if (lhs.frequencies.size() > maximumFrequencies)
        return std::nullopt;
    return lhs;
}

[[nodiscard]] std::optional<HarmonicSpectrum> multiplySpectra(
    const HarmonicSpectrum& lhs,
    const HarmonicSpectrum& rhs,
    std::size_t maximumFrequencies) {
    std::vector<Rational> result;
    if (lhs.frequencies.empty() || rhs.frequencies.empty())
        return HarmonicSpectrum{};
    if (lhs.frequencies.size() > maximumFrequencies
        || rhs.frequencies.size() > maximumFrequencies)
        return std::nullopt;

    // sin/cos積の周波数supportは |a-b| と a+b の部分集合。
    // 位相・係数を捨てるためsupportは過大評価し得るが，得られる共通周期は常に安全。
    for (const auto& a : lhs.frequencies) {
        for (const auto& b : rhs.frequencies) {
            result.push_back(positive(a - b));
            result.push_back(a + b);
            if (result.size() > maximumFrequencies * 4U)
                return std::nullopt;
        }
    }
    normalizeFrequencies(result);
    if (result.size() > maximumFrequencies)
        return std::nullopt;
    return HarmonicSpectrum{std::move(result)};
}

[[nodiscard]] std::optional<HarmonicSpectrum> analyzeSpectrumImpl(
    const Expr& expression,
    const expression::Symbol& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const HarmonicSpectrumOptions& options) {
    if (!symbolic::containsSymbol(expression, parameter))
        return HarmonicSpectrum{{Rational{BigInt{0}}}};
    if (!expression.isCall())
        return std::nullopt;

    const auto& call = expression.asCall();
    const auto* builtin = builtins.find(call.head);
    if (!builtin)
        return std::nullopt;

    if ((builtin->id == BuiltinId::Sin || builtin->id == BuiltinId::Cos)
        && call.arguments.size() == 1) {
        const auto coefficient = affineCoefficient(
            call.arguments[0], parameter, builtins, mathematics, angles);
        if (!coefficient)
            return std::nullopt;
        return HarmonicSpectrum{{positive(*coefficient)}};
    }

    if (builtin->id == BuiltinId::Negate && call.arguments.size() == 1)
        return analyzeSpectrumImpl(
            call.arguments[0], parameter, builtins, mathematics, angles, options);

    if ((builtin->id == BuiltinId::Add || builtin->id == BuiltinId::Subtract)
        && !call.arguments.empty()) {
        auto result = analyzeSpectrumImpl(
            call.arguments.front(), parameter, builtins, mathematics, angles, options);
        if (!result)
            return std::nullopt;
        for (std::size_t i = 1; i < call.arguments.size(); ++i) {
            const auto child = analyzeSpectrumImpl(
                call.arguments[i], parameter, builtins, mathematics, angles, options);
            if (!child)
                return std::nullopt;
            result = mergeSpectra(std::move(*result), *child, options.maximumFrequencies);
            if (!result)
                return std::nullopt;
        }
        return result;
    }

    if (builtin->id == BuiltinId::Multiply && !call.arguments.empty()) {
        HarmonicSpectrum result{{Rational{BigInt{0}}}};
        for (const auto& argument : call.arguments) {
            const auto child = analyzeSpectrumImpl(
                argument, parameter, builtins, mathematics, angles, options);
            if (!child)
                return std::nullopt;
            auto product = multiplySpectra(result, *child, options.maximumFrequencies);
            if (!product)
                return std::nullopt;
            result = std::move(*product);
        }
        return result;
    }

    if (builtin->id == BuiltinId::Divide && call.arguments.size() == 2
        && !symbolic::containsSymbol(call.arguments[1], parameter)) {
        return analyzeSpectrumImpl(
            call.arguments[0], parameter, builtins, mathematics, angles, options);
    }

    if (builtin->id == BuiltinId::Power && call.arguments.size() == 2) {
        const auto exponent = nonnegativeSmallInteger(call.arguments[1]);
        if (!exponent || *exponent > options.maximumPower)
            return std::nullopt;
        if (*exponent == 0)
            return HarmonicSpectrum{{Rational{BigInt{0}}}};
        const auto base = analyzeSpectrumImpl(
            call.arguments[0], parameter, builtins, mathematics, angles, options);
        if (!base)
            return std::nullopt;
        HarmonicSpectrum result{{Rational{BigInt{0}}}};
        for (std::size_t i = 0; i < *exponent; ++i) {
            auto product = multiplySpectra(result, *base, options.maximumFrequencies);
            if (!product)
                return std::nullopt;
            result = std::move(*product);
        }
        return result;
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<Rational> periodFromSpectrum(const HarmonicSpectrum& spectrum) {
    std::optional<Rational> fundamentalFrequency;
    for (const auto& frequency : spectrum.frequencies) {
        if (frequency.isZero())
            continue;
        fundamentalFrequency = fundamentalFrequency
            ? std::optional<Rational>{rationalGcd(*fundamentalFrequency, frequency)}
            : std::optional<Rational>{positive(frequency)};
    }
    if (!fundamentalFrequency || fundamentalFrequency->isZero())
        return std::nullopt;
    return Rational{BigInt{1}} / *fundamentalFrequency;
}

[[nodiscard]] PeriodKnowledge analyzeExpressionPeriod(
    const Expr& expression,
    const expression::Symbol& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!symbolic::containsSymbol(expression, parameter))
        return PeriodKnowledge{};
    if (expression.isSymbol())
        return PeriodKnowledge{true, std::nullopt, PeriodProofKind::Composition};
    if (!expression.isCall())
        return PeriodKnowledge{true, std::nullopt, PeriodProofKind::Composition};

    const auto& call = expression.asCall();
    const auto* definition = mathematics.findFunction(call.head);
    if (definition && definition->periodTurns && call.arguments.size() == 1) {
        if (const auto coefficient = affineCoefficient(
                call.arguments[0], parameter, builtins, mathematics, angles)) {
            return PeriodKnowledge{
                true, positive(*definition->periodTurns / positive(*coefficient)),
                PeriodProofKind::Registry};
        }
    }

    bool anyDependent = false;
    std::optional<Rational> period;
    bool registryOnly = true;
    for (const auto& argument : call.arguments) {
        const auto child = analyzeExpressionPeriod(
            argument, parameter, builtins, mathematics, angles);
        if (!child.dependsOnParameter)
            continue;
        anyDependent = true;
        if (!child.periodTurns)
            return PeriodKnowledge{true, std::nullopt, PeriodProofKind::Composition};
        period = commonPeriod(period, child.periodTurns);
        registryOnly = registryOnly && child.proofKind == PeriodProofKind::Registry;
    }
    return PeriodKnowledge{
        anyDependent, period,
        registryOnly ? PeriodProofKind::Registry : PeriodProofKind::Composition};
}

[[nodiscard]] Expr makePeriodExpression(
    const Rational& turns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Degree:
        return Expr{numeric::Number{turns * Rational{BigInt{360}}}};
    case mathematics::AngleUnit::Gradian:
        return Expr{numeric::Number{turns * Rational{BigInt{400}}}};
    case mathematics::AngleUnit::Radian: {
        const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
        const Rational coefficient = turns * Rational{BigInt{2}};
        if (!pi)
            return Expr{numeric::Number{BigInt{0}}};
        if (coefficient == Rational{BigInt{1}})
            return Expr{pi->symbol};
        return Expr::call(
            builtins.symbol(BuiltinId::Multiply),
            {Expr{numeric::Number{coefficient}}, Expr{pi->symbol}});
    }
    }
    return Expr{numeric::Number{BigInt{0}}};
}

[[nodiscard]] std::optional<std::size_t> exactPositiveInteger(const Expr& expression) {
    const auto rational = exactRealRational(expression);
    if (!rational || !rational->isInteger() || rational->numerator().isNegative()
        || rational->numerator().isZero())
        return std::nullopt;
    try {
        const auto text = rational->numerator().toString();
        const auto value = std::stoull(text);
        if (value > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()))
            return std::nullopt;
        return static_cast<std::size_t>(value);
    }
    catch (...) {
        return std::nullopt;
    }
}

} // namespace

std::optional<HarmonicSpectrum> analyzeHarmonicSpectrum(
    const expression::Expr& expression,
    const expression::Symbol& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const HarmonicSpectrumOptions& options) {
    if (options.maximumFrequencies == 0)
        return std::nullopt;
    auto result = analyzeSpectrumImpl(
        expression, parameter, builtins, mathematics, angles, options);
    if (!result)
        return std::nullopt;
    normalizeFrequencies(result->frequencies);
    return result;
}

expression::Expr parametricPeriodExpression(
    const numeric::Rational& turns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return makePeriodExpression(turns, builtins, mathematics, angles);
}

std::optional<std::size_t> exactParametricPeriodRepetitions(
    const ParametricCurveRequest& request,
    const numeric::Rational& periodTurns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    simplification::SimplificationContext context{builtins, mathematics, angles, assumptions};
    const Expr period = simplification::Simplifier{}.simplify(
        parametricPeriodExpression(periodTurns, builtins, mathematics, angles), context);
    const Expr span = Expr::call(
        builtins.symbol(BuiltinId::Subtract), {request.upper, request.lower});
    const Expr ratio = simplification::Simplifier{}.simplify(
        Expr::call(builtins.symbol(BuiltinId::Divide), {span, period}), context);
    return exactPositiveInteger(ratio);
}

std::optional<PeriodCertificate> analyzeParametricCurvePeriod(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const HarmonicSpectrumOptions& harmonicOptions) {
    const auto xSpectrum = analyzeHarmonicSpectrum(
        request.xExpression, request.parameter, builtins, mathematics, angles, harmonicOptions);
    const auto ySpectrum = analyzeHarmonicSpectrum(
        request.yExpression, request.parameter, builtins, mathematics, angles, harmonicOptions);
    if (xSpectrum && ySpectrum) {
        HarmonicSpectrum combined = *xSpectrum;
        const auto merged = mergeSpectra(
            std::move(combined), *ySpectrum, harmonicOptions.maximumFrequencies);
        if (merged) {
            if (const auto period = periodFromSpectrum(*merged))
                return PeriodCertificate{
                    *period, false, PeriodProofKind::HarmonicSpectrum, *merged};
        }
    }

    const auto x = analyzeExpressionPeriod(
        request.xExpression, request.parameter, builtins, mathematics, angles);
    const auto y = analyzeExpressionPeriod(
        request.yExpression, request.parameter, builtins, mathematics, angles);

    if (x.dependsOnParameter && !x.periodTurns)
        return std::nullopt;
    if (y.dependsOnParameter && !y.periodTurns)
        return std::nullopt;
    const auto period = commonPeriod(x.periodTurns, y.periodTurns);
    if (!period || period->isZero())
        return std::nullopt;
    const PeriodProofKind proof = x.proofKind == PeriodProofKind::Registry
            && y.proofKind == PeriodProofKind::Registry
        ? PeriodProofKind::Registry
        : PeriodProofKind::Composition;
    return PeriodCertificate{*period, false, proof, std::nullopt};
}

std::optional<ParametricPeriodReduction> reduceParametricCurvePeriod(
    ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    const HarmonicSpectrumOptions& harmonicOptions) {
    const auto analysis = analyzeParametricCurvePeriod(
        request, builtins, mathematics, angles, harmonicOptions);
    if (!analysis)
        return std::nullopt;

    simplification::SimplificationContext context{builtins, mathematics, angles, assumptions};
    const Expr period = simplification::Simplifier{}.simplify(
        parametricPeriodExpression(analysis->turns, builtins, mathematics, angles), context);
    const auto repetitions = exactParametricPeriodRepetitions(
        request, analysis->turns, builtins, mathematics, angles, assumptions);
    if (!repetitions || *repetitions < 2)
        return std::nullopt;

    const Expr effectiveUpper = simplification::Simplifier{}.simplify(
        Expr::call(builtins.symbol(BuiltinId::Add), {request.lower, period}), context);
    ParametricPeriodReduction reduction{
        request.lower, request.upper, effectiveUpper, period, *repetitions, *analysis};
    request.upper = effectiveUpper;
    return reduction;
}

} // namespace mmcal::plot
