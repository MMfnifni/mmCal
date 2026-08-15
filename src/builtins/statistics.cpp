// 統計
#include "statistics.hpp"

#include "builtins/exact_operations.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
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

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr natural(std::size_t value) {
    return Expr{Number{BigInt::parse(std::to_string(value))}};
}

[[nodiscard]] Rational rationalOf(const Expr& expression, std::string_view name) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        error::throwCalcError(error::CalcErrorType::Type,
            std::string{name} + " requires exact real values");
    return expression.asNumber().asReal().toRational();
}

void requireExactRealData(std::span<const Expr> values, std::string_view name) {
    for (const Expr& value : values)
        static_cast<void>(rationalOf(value, name));
}

[[nodiscard]] std::vector<Expr> scalarOrVector(
    std::span<const Expr> arguments,
    std::string_view name,
    std::size_t minimumCount = 1) {
    if (arguments.size() == 1 && arguments.front().isArray()) {
        const auto& array = arguments.front().asArray();
        if (array.rank() != 1)
            error::throwCalcError(error::CalcErrorType::Type,
                std::string{name} + " requires a rank-1 array");
        if (array.size() < minimumCount)
            error::throwCalcError(error::CalcErrorType::Domain,
                std::string{name} + " requires more observations");
        return array.materialize();
    }

    if (arguments.size() < minimumCount)
        error::throwCalcError(error::CalcErrorType::Domain,
            std::string{name} + " requires more observations");
    for (const Expr& argument : arguments)
        if (argument.isArray())
            error::throwCalcError(error::CalcErrorType::Type,
                std::string{name} + " accepts either one vector or scalar observations");
    return {arguments.begin(), arguments.end()};
}

[[nodiscard]] std::vector<Expr> tailData(
    std::span<const Expr> arguments,
    std::string_view name,
    std::size_t minimumCount = 1) {
    if (arguments.size() == 2 && arguments[1].isArray()) {
        const auto& array = arguments[1].asArray();
        if (array.rank() != 1)
            error::throwCalcError(error::CalcErrorType::Type,
                std::string{name} + " requires a rank-1 data array");
        if (array.size() < minimumCount)
            error::throwCalcError(error::CalcErrorType::Domain,
                std::string{name} + " requires more observations");
        return array.materialize();
    }

    if (arguments.size() - 1 < minimumCount)
        error::throwCalcError(error::CalcErrorType::Domain,
            std::string{name} + " requires more observations");
    std::vector<Expr> result;
    result.reserve(arguments.size() - 1);
    for (std::size_t i = 1; i < arguments.size(); ++i) {
        if (arguments[i].isArray())
            error::throwCalcError(error::CalcErrorType::Type,
                std::string{name} + " accepts p/x followed by one vector or scalar observations");
        result.push_back(arguments[i]);
    }
    return result;
}

struct PairedData final {
    std::vector<Expr> x;
    std::vector<Expr> y;
};

[[nodiscard]] PairedData pairedData(std::span<const Expr> arguments, std::string_view name) {
    if (arguments.size() == 2 && arguments[0].isArray() && arguments[1].isArray()) {
        const auto& x = arguments[0].asArray();
        const auto& y = arguments[1].asArray();
        if (x.rank() != 1 || y.rank() != 1)
            error::throwCalcError(error::CalcErrorType::Type,
                std::string{name} + " requires rank-1 arrays");
        if (x.size() != y.size() || x.empty())
            error::throwCalcError(error::CalcErrorType::Domain,
                std::string{name} + " requires equal non-empty sample lengths");
        return {x.materialize(), y.materialize()};
    }

    if (arguments.size() < 2 || arguments.size() % 2 != 0)
        error::throwCalcError(error::CalcErrorType::Type,
            std::string{name} + " requires two vectors or an even scalar list split in half");
    const std::size_t half = arguments.size() / 2;
    for (const Expr& argument : arguments)
        if (argument.isArray())
            error::throwCalcError(error::CalcErrorType::Type,
                std::string{name} + " cannot mix arrays and scalar observations");
    return {{arguments.begin(), arguments.begin() + static_cast<std::ptrdiff_t>(half)},
            {arguments.begin() + static_cast<std::ptrdiff_t>(half), arguments.end()}};
}

[[nodiscard]] bool exactRealLess(const Expr& lhs, const Expr& rhs) {
    return lhs.asNumber().asReal() < rhs.asNumber().asReal();
}

[[nodiscard]] std::vector<Expr> sortedExactReal(std::vector<Expr> values, std::string_view name) {
    requireExactRealData(values, name);
    std::sort(values.begin(), values.end(), exactRealLess);
    return values;
}

[[nodiscard]] Expr addAll(
    std::vector<Expr> values,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return exact::add(std::move(values), registry, mathematics, angles);
}

[[nodiscard]] Expr meanOf(
    const std::vector<Expr>& values,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> copy = values;
    return exact::divide(addAll(std::move(copy), registry, mathematics, angles), natural(values.size()),
        registry, mathematics, angles);
}

[[nodiscard]] Expr powerInteger(
    Expr value,
    std::int64_t exponent,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (exponent == 0)
        return integer(1);
    if (exponent > 0 && exponent <= 8) {
        std::vector<Expr> factors(static_cast<std::size_t>(exponent), value);
        return exact::multiply(std::move(factors), registry, mathematics, angles);
    }
    return exact::call(BuiltinId::Power, {std::move(value), integer(exponent)},
        registry, mathematics, angles);
}

[[nodiscard]] Expr populationCentralMoment(
    const std::vector<Expr>& values,
    std::int64_t exponent,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Expr mean = meanOf(values, registry, mathematics, angles);
    std::vector<Expr> terms;
    terms.reserve(values.size());
    for (const Expr& value : values) {
        Expr deviation = exact::subtract(value, mean, registry, mathematics, angles);
        terms.push_back(powerInteger(std::move(deviation), exponent, registry, mathematics, angles));
    }
    return exact::divide(addAll(std::move(terms), registry, mathematics, angles), natural(values.size()),
        registry, mathematics, angles);
}

[[nodiscard]] Expr varianceOf(
    const std::vector<Expr>& values,
    bool sample,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (sample && values.size() < 2)
        error::throwCalcError(error::CalcErrorType::Domain, "sample variance requires at least two observations");
    const Expr mean = meanOf(values, registry, mathematics, angles);
    std::vector<Expr> squares;
    squares.reserve(values.size());
    for (const Expr& value : values) {
        Expr deviation = exact::subtract(value, mean, registry, mathematics, angles);
        squares.push_back(powerInteger(std::move(deviation), 2, registry, mathematics, angles));
    }
    const std::size_t denominator = sample ? values.size() - 1 : values.size();
    return exact::divide(addAll(std::move(squares), registry, mathematics, angles), natural(denominator),
        registry, mathematics, angles);
}

[[nodiscard]] bool isExactZero(const Expr& value) {
    return value.isNumber() && value.asNumber().isZero();
}

[[nodiscard]] Expr medianSorted(
    const std::vector<Expr>& sorted,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::size_t n = sorted.size();
    if (n % 2 != 0)
        return sorted[n / 2];
    return exact::divide(exact::add({sorted[n / 2 - 1], sorted[n / 2]}, registry, mathematics, angles),
        integer(2), registry, mathematics, angles);
}

[[nodiscard]] std::size_t floorNonnegativeRational(const Rational& value) {
    if (value.numerator().isNegative())
        error::throwCalcError(error::CalcErrorType::Domain, "internal nonnegative rational expected");
    const BigInt q = value.numerator() / value.denominator();
    return static_cast<std::size_t>(std::stoull(q.toString()));
}

[[nodiscard]] Expr quantileSorted(
    const std::vector<Expr>& sorted,
    const Rational& p,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Rational zero{BigInt{0}};
    const Rational one{BigInt{1}};
    if (p < zero || p > one)
        error::throwCalcError(error::CalcErrorType::Domain, "quantile probability must be in [0, 1]");
    if (sorted.size() == 1)
        return sorted.front();

    const Rational h = p * Rational{BigInt::parse(std::to_string(sorted.size() - 1))};
    const std::size_t lower = floorNonnegativeRational(h);
    if (lower >= sorted.size() - 1)
        return sorted.back();
    const Rational fraction = h - Rational{BigInt::parse(std::to_string(lower))};
    if (fraction.isZero())
        return sorted[lower];

    Expr delta = exact::subtract(sorted[lower + 1], sorted[lower], registry, mathematics, angles);
    Expr weight{Number{fraction}};
    return exact::add({sorted[lower], exact::multiply({std::move(weight), std::move(delta)}, registry, mathematics, angles)},
        registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateMedian(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = sortedExactReal(scalarOrVector(arguments, "median"), "median");
    return medianSorted(values, registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateMode(std::span<const Expr> arguments) {
    const auto values = scalarOrVector(arguments, "mode");
    std::vector<std::size_t> counts(values.size(), 0);
    std::size_t maximum = 0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        for (std::size_t j = 0; j < values.size(); ++j)
            if (values[i] == values[j])
                ++counts[i];
        maximum = std::max(maximum, counts[i]);
    }
    if (maximum <= 1)
        return Expr::array({0}, {});

    std::vector<Expr> modes;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (counts[i] != maximum)
            continue;
        bool duplicate = false;
        for (const Expr& existing : modes)
            duplicate = duplicate || existing == values[i];
        if (!duplicate)
            modes.push_back(values[i]);
    }
    if (modes.size() == 1)
        return modes.front();
    const std::size_t modeCount = modes.size();
    return Expr::array({modeCount}, std::move(modes));
}

[[nodiscard]] Expr evaluateQuantileLike(
    std::span<const Expr> arguments,
    bool percentile,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Rational p = rationalOf(arguments.front(), percentile ? "percentile" : "quantile");
    if (percentile)
        p /= Rational{BigInt{100}};
    auto values = sortedExactReal(tailData(arguments, percentile ? "percentile" : "quantile"),
        percentile ? "percentile" : "quantile");
    return quantileSorted(values, p, registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateVarianceLike(
    std::span<const Expr> arguments,
    bool sample,
    bool standardDeviation,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, sample ? (standardDeviation ? "stddevs" : "vars")
                                                   : (standardDeviation ? "stddev" : "var"), sample ? 2 : 1);
    requireExactRealData(values, sample ? "sample statistic" : "population statistic");
    Expr variance = varianceOf(values, sample, registry, mathematics, angles);
    if (!standardDeviation)
        return variance;
    return exact::sqrt(std::move(variance), registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateGeometricMean(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, "geomean");
    requireExactRealData(values, "geomean");
    for (const Expr& value : values)
        if (value.asNumber().asReal().isNegative())
            error::throwCalcError(error::CalcErrorType::Domain, "geomean requires nonnegative real observations");
    Expr product = exact::multiply(values, registry, mathematics, angles);
    if (values.size() == 1)
        return product;
    if (values.size() == 2)
        return exact::sqrt(std::move(product), registry, mathematics, angles);
    if (values.size() == 3)
        return exact::call(BuiltinId::Cbrt, {std::move(product)}, registry, mathematics, angles);
    const Rational exponent{BigInt{1}, BigInt::parse(std::to_string(values.size()))};
    return exact::call(BuiltinId::Power, {std::move(product), Expr{Number{exponent}}},
        registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateHarmonicMean(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, "harmmean");
    requireExactRealData(values, "harmmean");
    std::vector<Expr> reciprocals;
    reciprocals.reserve(values.size());
    for (const Expr& value : values) {
        if (value.asNumber().isZero())
            error::throwCalcError(error::CalcErrorType::Domain, "harmmean is undefined for zero observations");
        reciprocals.push_back(exact::divide(integer(1), value, registry, mathematics, angles));
    }
    Expr denominator = addAll(std::move(reciprocals), registry, mathematics, angles);
    if (isExactZero(denominator))
        error::throwCalcError(error::CalcErrorType::Domain, "harmmean reciprocal sum is zero");
    return exact::divide(natural(values.size()), std::move(denominator), registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateRms(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, "rms");
    requireExactRealData(values, "rms");
    std::vector<Expr> squares;
    squares.reserve(values.size());
    for (const Expr& value : values)
        squares.push_back(powerInteger(value, 2, registry, mathematics, angles));
    Expr meanSquare = exact::divide(addAll(std::move(squares), registry, mathematics, angles), natural(values.size()),
        registry, mathematics, angles);
    return exact::sqrt(std::move(meanSquare), registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateAbsoluteDeviation(
    std::span<const Expr> arguments,
    bool fromMedian,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, fromMedian ? "mad" : "madR");
    requireExactRealData(values, fromMedian ? "mad" : "madR");
    Expr center = fromMedian
        ? medianSorted(sortedExactReal(values, "mad"), registry, mathematics, angles)
        : meanOf(values, registry, mathematics, angles);
    std::vector<Expr> deviations;
    deviations.reserve(values.size());
    for (const Expr& value : values) {
        Expr difference = exact::subtract(value, center, registry, mathematics, angles);
        deviations.push_back(exact::call(BuiltinId::Abs, {std::move(difference)}, registry, mathematics, angles));
    }
    if (fromMedian) {
        auto sorted = sortedExactReal(std::move(deviations), "mad");
        return medianSorted(sorted, registry, mathematics, angles);
    }
    return meanOf(deviations, registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateSkewness(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, "skew", 2);
    requireExactRealData(values, "skew");
    Expr m2 = populationCentralMoment(values, 2, registry, mathematics, angles);
    if (isExactZero(m2))
        error::throwCalcError(error::CalcErrorType::Domain, "skew is undefined for zero variance");
    Expr m3 = populationCentralMoment(values, 3, registry, mathematics, angles);
    Expr denominator = exact::multiply({m2, exact::sqrt(m2, registry, mathematics, angles)},
        registry, mathematics, angles);
    return exact::divide(std::move(m3), std::move(denominator), registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateKurtosisPopulation(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, "kurtp", 2);
    requireExactRealData(values, "kurtp");
    Expr m2 = populationCentralMoment(values, 2, registry, mathematics, angles);
    if (isExactZero(m2))
        error::throwCalcError(error::CalcErrorType::Domain, "kurtp is undefined for zero variance");
    Expr m4 = populationCentralMoment(values, 4, registry, mathematics, angles);
    Expr ratio = exact::divide(std::move(m4), powerInteger(m2, 2, registry, mathematics, angles),
        registry, mathematics, angles);
    return exact::subtract(std::move(ratio), integer(3), registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateKurtosisSample(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, "kurts", 4);
    requireExactRealData(values, "kurts");
    const std::size_t n = values.size();
    Expr mean = meanOf(values, registry, mathematics, angles);
    std::vector<Expr> fourthPowers;
    fourthPowers.reserve(n);
    for (const Expr& value : values) {
        Expr d = exact::subtract(value, mean, registry, mathematics, angles);
        fourthPowers.push_back(powerInteger(std::move(d), 4, registry, mathematics, angles));
    }
    Expr s2 = varianceOf(values, true, registry, mathematics, angles);
    if (isExactZero(s2))
        error::throwCalcError(error::CalcErrorType::Domain, "kurts is undefined for zero variance");
    Expr normalizedFourth = exact::divide(addAll(std::move(fourthPowers), registry, mathematics, angles),
        powerInteger(s2, 2, registry, mathematics, angles), registry, mathematics, angles);

    const BigInt nBig = BigInt::parse(std::to_string(n));
    const BigInt nm1 = BigInt::parse(std::to_string(n - 1));
    const BigInt nm2 = BigInt::parse(std::to_string(n - 2));
    const BigInt nm3 = BigInt::parse(std::to_string(n - 3));
    const Rational a{nBig * (nBig + BigInt{1}), nm1 * nm2 * nm3};
    const Rational b{BigInt{3} * nm1 * nm1, nm2 * nm3};
    Expr first = exact::multiply({Expr{Number{a}}, std::move(normalizedFourth)}, registry, mathematics, angles);
    return exact::subtract(std::move(first), Expr{Number{b}}, registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateCv(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, "cv");
    requireExactRealData(values, "cv");
    Expr mean = meanOf(values, registry, mathematics, angles);
    if (isExactZero(mean))
        error::throwCalcError(error::CalcErrorType::Domain, "cv is undefined for zero mean");
    Expr sd = exact::sqrt(varianceOf(values, false, registry, mathematics, angles), registry, mathematics, angles);
    return exact::divide(std::move(sd), std::move(mean), registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateStandardError(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto values = scalarOrVector(arguments, "stderr", 2);
    requireExactRealData(values, "stderr");
    Expr sd = exact::sqrt(varianceOf(values, true, registry, mathematics, angles), registry, mathematics, angles);
    Expr rootN = exact::sqrt(natural(values.size()), registry, mathematics, angles);
    return exact::divide(std::move(sd), std::move(rootN), registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateZScore(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    for (const Expr& argument : arguments)
        static_cast<void>(rationalOf(argument, "zscore"));
    if ((arguments[2].asNumber().asReal() <=> numeric::RealNumber{BigInt{0}}) != std::strong_ordering::greater)
        error::throwCalcError(error::CalcErrorType::Domain, "zscore requires sigma > 0");
    return exact::divide(exact::subtract(arguments[0], arguments[1], registry, mathematics, angles),
        arguments[2], registry, mathematics, angles);
}

[[nodiscard]] Expr evaluateIqr(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto sorted = sortedExactReal(scalarOrVector(arguments, "iqr", 2), "iqr");
    const Expr q1 = quantileSorted(sorted, Rational{BigInt{1}, BigInt{4}}, registry, mathematics, angles);
    const Expr q3 = quantileSorted(sorted, Rational{BigInt{3}, BigInt{4}}, registry, mathematics, angles);
    return exact::subtract(q3, q1, registry, mathematics, angles);
}

[[nodiscard]] Rational trimmingFraction(const Expr& value, std::string_view name) {
    Rational p = rationalOf(value, name);
    const Rational zero{BigInt{0}};
    const Rational half{BigInt{1}, BigInt{2}};
    if (p < zero || !(p < half))
        error::throwCalcError(error::CalcErrorType::Domain,
            std::string{name} + " requires 0 <= p < 1/2");
    return p;
}

[[nodiscard]] std::size_t trimCount(const Rational& p, std::size_t n) {
    return floorNonnegativeRational(p * Rational{BigInt::parse(std::to_string(n))});
}

[[nodiscard]] Expr evaluateTrimMean(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Rational p = trimmingFraction(arguments.front(), "trimmean");
    auto sorted = sortedExactReal(tailData(arguments, "trimmean"), "trimmean");
    const std::size_t k = trimCount(p, sorted.size());
    if (2 * k >= sorted.size())
        error::throwCalcError(error::CalcErrorType::Domain, "trimmean removes all observations");
    std::vector<Expr> kept(sorted.begin() + static_cast<std::ptrdiff_t>(k),
        sorted.end() - static_cast<std::ptrdiff_t>(k));
    return meanOf(kept, registry, mathematics, angles);
}

[[nodiscard]] std::vector<Expr> winsorizedValues(
    std::span<const Expr> arguments,
    std::string_view name) {
    const Rational p = trimmingFraction(arguments.front(), name);
    auto original = tailData(arguments, name);
    requireExactRealData(original, name);
    auto sorted = sortedExactReal(original, name);
    const std::size_t k = trimCount(p, sorted.size());
    if (k == 0)
        return original;
    const Expr low = sorted[k];
    const Expr high = sorted[sorted.size() - k - 1];
    for (Expr& value : original) {
        if (exactRealLess(value, low))
            value = low;
        else if (exactRealLess(high, value))
            value = high;
    }
    return original;
}

[[nodiscard]] Expr covarianceOf(
    const std::vector<Expr>& x,
    const std::vector<Expr>& y,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr mx = meanOf(x, registry, mathematics, angles);
    Expr my = meanOf(y, registry, mathematics, angles);
    std::vector<Expr> products;
    products.reserve(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        Expr dx = exact::subtract(x[i], mx, registry, mathematics, angles);
        Expr dy = exact::subtract(y[i], my, registry, mathematics, angles);
        products.push_back(exact::multiply({std::move(dx), std::move(dy)}, registry, mathematics, angles));
    }
    return exact::divide(addAll(std::move(products), registry, mathematics, angles), natural(x.size()),
        registry, mathematics, angles);
}

[[nodiscard]] Expr correlationOf(
    const std::vector<Expr>& x,
    const std::vector<Expr>& y,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr vx = varianceOf(x, false, registry, mathematics, angles);
    Expr vy = varianceOf(y, false, registry, mathematics, angles);
    if (isExactZero(vx) || isExactZero(vy))
        error::throwCalcError(error::CalcErrorType::Domain, "correlation is undefined for zero variance");
    Expr denominator = exact::sqrt(exact::multiply({vx, vy}, registry, mathematics, angles),
        registry, mathematics, angles);
    return exact::divide(covarianceOf(x, y, registry, mathematics, angles), std::move(denominator),
        registry, mathematics, angles);
}

[[nodiscard]] std::vector<Expr> averageRanks(
    const std::vector<Expr>& values) {
    std::vector<std::size_t> order(values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
        order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return exactRealLess(values[a], values[b]);
    });

    std::vector<Expr> ranks(values.size(), integer(0));
    std::size_t begin = 0;
    while (begin < order.size()) {
        std::size_t end = begin + 1;
        while (end < order.size() && values[order[begin]] == values[order[end]])
            ++end;
        const BigInt numerator = BigInt::parse(std::to_string(begin + 1 + end));
        const Expr rank{Number{Rational{numerator, BigInt{2}}}};
        for (std::size_t i = begin; i < end; ++i)
            ranks[order[i]] = rank;
        begin = end;
    }
    return ranks;
}

[[nodiscard]] Expr evaluatePercentRank(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Expr x = arguments.front();
    static_cast<void>(rationalOf(x, "percentrank"));
    auto sorted = sortedExactReal(tailData(arguments, "percentrank"), "percentrank");
    if (sorted.size() == 1) {
        if (x == sorted.front())
            return integer(0);
        error::throwCalcError(error::CalcErrorType::Domain, "percentrank x is outside the data range");
    }
    if (exactRealLess(x, sorted.front()) || exactRealLess(sorted.back(), x))
        error::throwCalcError(error::CalcErrorType::Domain, "percentrank x is outside the data range");

    for (std::size_t i = 0; i < sorted.size(); ++i) {
        if (!(x == sorted[i]))
            continue;
        std::size_t end = i + 1;
        while (end < sorted.size() && x == sorted[end])
            ++end;
        const BigInt numerator = BigInt::parse(std::to_string(i + end - 1));
        const BigInt denominator = BigInt::parse(std::to_string(2 * (sorted.size() - 1)));
        return Expr{Number{Rational{numerator, denominator}}};
    }

    for (std::size_t i = 0; i + 1 < sorted.size(); ++i) {
        if (exactRealLess(sorted[i], x) && exactRealLess(x, sorted[i + 1])) {
            Expr fraction = exact::divide(exact::subtract(x, sorted[i], registry, mathematics, angles),
                exact::subtract(sorted[i + 1], sorted[i], registry, mathematics, angles),
                registry, mathematics, angles);
            Expr position = exact::add({natural(i), std::move(fraction)}, registry, mathematics, angles);
            return exact::divide(std::move(position), natural(sorted.size() - 1), registry, mathematics, angles);
        }
    }
    return integer(1);
}

} // namespace

Expr evaluateStatistic(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (id) {
    case BuiltinId::Median:
        return evaluateMedian(arguments, registry, mathematics, angles);
    case BuiltinId::Mode:
        return evaluateMode(arguments);
    case BuiltinId::Quantile:
        return evaluateQuantileLike(arguments, false, registry, mathematics, angles);
    case BuiltinId::Percentile:
        return evaluateQuantileLike(arguments, true, registry, mathematics, angles);
    case BuiltinId::VariancePopulation:
        return evaluateVarianceLike(arguments, false, false, registry, mathematics, angles);
    case BuiltinId::VarianceSample:
        return evaluateVarianceLike(arguments, true, false, registry, mathematics, angles);
    case BuiltinId::StddevPopulation:
        return evaluateVarianceLike(arguments, false, true, registry, mathematics, angles);
    case BuiltinId::StddevSample:
        return evaluateVarianceLike(arguments, true, true, registry, mathematics, angles);
    case BuiltinId::GeometricMean:
        return evaluateGeometricMean(arguments, registry, mathematics, angles);
    case BuiltinId::HarmonicMean:
        return evaluateHarmonicMean(arguments, registry, mathematics, angles);
    case BuiltinId::Rms:
        return evaluateRms(arguments, registry, mathematics, angles);
    case BuiltinId::MedianAbsoluteDeviation:
        return evaluateAbsoluteDeviation(arguments, true, registry, mathematics, angles);
    case BuiltinId::MeanAbsoluteDeviation:
        return evaluateAbsoluteDeviation(arguments, false, registry, mathematics, angles);
    case BuiltinId::Skewness:
        return evaluateSkewness(arguments, registry, mathematics, angles);
    case BuiltinId::KurtosisPopulation:
        return evaluateKurtosisPopulation(arguments, registry, mathematics, angles);
    case BuiltinId::KurtosisSample:
        return evaluateKurtosisSample(arguments, registry, mathematics, angles);
    case BuiltinId::CoefficientVariation:
        return evaluateCv(arguments, registry, mathematics, angles);
    case BuiltinId::StandardError:
        return evaluateStandardError(arguments, registry, mathematics, angles);
    case BuiltinId::ZScore:
        return evaluateZScore(arguments, registry, mathematics, angles);
    case BuiltinId::Iqr:
        return evaluateIqr(arguments, registry, mathematics, angles);
    case BuiltinId::TrimMean:
        return evaluateTrimMean(arguments, registry, mathematics, angles);
    case BuiltinId::Winsorized: {
        auto values = winsorizedValues(arguments, "winsorR");
        const std::size_t count = values.size();
        return Expr::array({count}, std::move(values));
    }
    case BuiltinId::WinsorMean: {
        auto values = winsorizedValues(arguments, "winsor");
        return meanOf(values, registry, mathematics, angles);
    }
    case BuiltinId::Covariance: {
        PairedData data = pairedData(arguments, "cov");
        requireExactRealData(data.x, "cov");
        requireExactRealData(data.y, "cov");
        return covarianceOf(data.x, data.y, registry, mathematics, angles);
    }
    case BuiltinId::Correlation: {
        PairedData data = pairedData(arguments, "corr");
        requireExactRealData(data.x, "corr");
        requireExactRealData(data.y, "corr");
        return correlationOf(data.x, data.y, registry, mathematics, angles);
    }
    case BuiltinId::SpearmanCorrelation: {
        PairedData data = pairedData(arguments, "corrspearman");
        requireExactRealData(data.x, "corrspearman");
        requireExactRealData(data.y, "corrspearman");
        return correlationOf(averageRanks(data.x), averageRanks(data.y), registry, mathematics, angles);
    }
    case BuiltinId::PercentRank:
        return evaluatePercentRank(arguments, registry, mathematics, angles);
    default:
        break;
    }
    error::throwCalcError(error::CalcErrorType::Type, "Unsupported statistics builtin");
}

} // namespace mmcal::builtins
