// DFT・FFT・畳み込み
#include "signal_processing.hpp"

#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/certified_trigonometry.hpp"
#include "approximation/complex_interval.hpp"
#include "approximation/expression_interval.hpp"
#include "builtins/exact_operations.hpp"
#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "symbolic/cyclotomic_field.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
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

void requireArity(std::span<const Expr> arguments, std::size_t expected, std::string_view name) {
    if (arguments.size() != expected)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects " + std::to_string(expected) + " argument(s)");
}

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] BigInt sizeInteger(std::size_t value) {
    return BigInt::parse(std::to_string(value));
}

[[nodiscard]] std::vector<Expr> vectorArgument(const Expr& expression, std::string_view name) {
    if (!expression.isArray() || expression.asArray().rank() != 1)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " requires a rank-1 array");
    return expression.asArray().materialize();
}

[[nodiscard]] Expr zero() {
    return integer(0);
}

[[nodiscard]] Expr multiply(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() * rhs.asNumber()};
    return exact::multiply({std::move(lhs), std::move(rhs)}, registry, mathematics, angles);
}

[[nodiscard]] Expr add(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() + rhs.asNumber()};
    return exact::add({std::move(lhs), std::move(rhs)}, registry, mathematics, angles);
}

[[nodiscard]] Expr subtract(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() - rhs.asNumber()};
    return exact::subtract(std::move(lhs), std::move(rhs), registry, mathematics, angles);
}

[[nodiscard]] Expr divideBySize(
    Expr value,
    std::size_t denominator,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (denominator == 0)
        return value;
    return exact::divide(
        std::move(value), Expr{Number{sizeInteger(denominator)}}, registry, mathematics, angles);
}

[[nodiscard]] Expr pi(const mathematics::MathRegistry& mathematics) {
    const auto* definition = mathematics.findConstant(mathematics::ConstantId::Pi);
    if (!definition)
        error::throwCalcError(error::CalcErrorType::Internal, "Pi is not registered");
    return Expr{definition->symbol};
}

// Fourier位相はセッションの既定角度に依存しない。常にradを明示する。
[[nodiscard]] Expr twiddle(
    std::size_t numerator,
    std::size_t denominator,
    bool inverse,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (denominator == 0)
        error::throwCalcError(error::CalcErrorType::Internal, "Fourier transform size is zero");

    BigInt signedNumerator = sizeInteger(numerator) * BigInt{2};
    if (!inverse)
        signedNumerator = -signedNumerator;
    const Rational coefficient{std::move(signedNumerator), sizeInteger(denominator)};

    Expr phase = exact::multiply(
        {Expr{Number{coefficient}}, pi(mathematics)}, registry, mathematics, angles);
    Expr radians = Expr::call(
        registry.symbol(BuiltinId::UnitApplied), {std::move(phase), Expr{std::string{"Rad"}}});
    return exact::call(BuiltinId::Cis, {std::move(radians)}, registry, mathematics, angles);
}

[[nodiscard]] std::vector<Expr> directTransform(
    const std::vector<Expr>& input,
    bool inverse,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::size_t n = input.size();
    if (n == 0)
        return {};

    std::vector<Expr> roots;
    roots.reserve(n);
    for (std::size_t m = 0; m < n; ++m)
        roots.push_back(twiddle(m, n, inverse, registry, mathematics, angles));

    std::vector<Expr> output;
    output.reserve(n);
    for (std::size_t k = 0; k < n; ++k) {
        Expr sum = zero();
        std::size_t rootIndex = 0;
        for (std::size_t j = 0; j < n; ++j) {
            sum = add(std::move(sum), multiply(input[j], roots[rootIndex],
                    registry, mathematics, angles), registry, mathematics, angles);

            if (j + 1 < n && k != 0) {
                if (rootIndex >= n - k)
                    rootIndex -= n - k;
                else
                    rootIndex += k;
            }
        }
        if (inverse)
            sum = divideBySize(std::move(sum), n, registry, mathematics, angles);
        output.push_back(std::move(sum));
    }
    return output;
}

[[nodiscard]] bool isPowerOfTwo(std::size_t value) noexcept {
    return value != 0 && (value & (value - 1)) == 0;
}

[[nodiscard]] FourierTransformCache::Plan& transformPlan(
    std::size_t n,
    FourierTransformCache& cache,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    auto& plan = cache.plan(n);
    if (!plan.bitReversed.empty())
        return plan;

    plan.bitReversed.resize(n);
    for (std::size_t i = 0, j = 0; i < n; ++i) {
        plan.bitReversed[i] = j;
        if (i + 1 == n)
            break;
        std::size_t bit = n >> 1;
        while ((j & bit) != 0) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
    }

    for (std::size_t length = 2; length <= n; length <<= 1) {
        const std::size_t half = length >> 1;
        FourierTransformCache::Stage stage;
        stage.length = length;
        stage.forwardRoots.reserve(half);
        stage.inverseRoots.reserve(half);
        for (std::size_t j = 0; j < half; ++j) {
            stage.forwardRoots.push_back(twiddle(
                j, length, false, registry, mathematics, angles));
            stage.inverseRoots.push_back(twiddle(
                j, length, true, registry, mathematics, angles));
        }
        plan.stages.push_back(std::move(stage));
        if (length == n)
            break;
    }
    return plan;
}

[[nodiscard]] std::vector<Expr> radix2Transform(
    const std::vector<Expr>& input,
    bool inverse,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    FourierTransformCache& cache) {
    const std::size_t n = input.size();
    if (!isPowerOfTwo(n))
        return directTransform(input, inverse, registry, mathematics, angles);

    auto& plan = transformPlan(n, cache, registry, mathematics, angles);
    std::vector<Expr> data(n, zero());
    for (std::size_t i = 0; i < n; ++i)
        data[plan.bitReversed[i]] = input[i];

    for (const auto& stage : plan.stages) {
        const std::size_t length = stage.length;
        const std::size_t half = length >> 1;
        const auto& roots = inverse ? stage.inverseRoots : stage.forwardRoots;
        for (std::size_t block = 0; block < n; block += length) {
            for (std::size_t j = 0; j < half; ++j) {
                Expr even = data[block + j];
                Expr odd = multiply(data[block + j + half], roots[j],
                    registry, mathematics, angles);
                data[block + j] = add(even, odd, registry, mathematics, angles);
                data[block + j + half] = subtract(std::move(even), std::move(odd),
                    registry, mathematics, angles);
            }
        }
    }

    if (inverse)
        for (Expr& value : data)
            value = divideBySize(std::move(value), n, registry, mathematics, angles);
    return data;
}

[[nodiscard]] Expr cyclotomicGeneratorExpr(
    const symbolic::CyclotomicFieldContext& field,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    // 表示は従来のexact Fourierと同じcisを使う。計算内部だけQ[t]/Phi_n(t)へ写す。
    return twiddle(1, field.conductor(), false, registry, mathematics, angles);
}

[[nodiscard]] std::vector<Rational> addCoordinates(
    std::span<const Rational> lhs,
    std::span<const Rational> rhs) {
    if (lhs.size() != rhs.size())
        return {};
    std::vector<Rational> result(lhs.begin(), lhs.end());
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] += rhs[i];
    return result;
}

[[nodiscard]] std::vector<Rational> negateCoordinates(std::span<const Rational> value) {
    std::vector<Rational> result(value.begin(), value.end());
    for (Rational& coefficient : result)
        coefficient = -coefficient;
    return result;
}

[[nodiscard]] std::optional<std::vector<Rational>> powerCoordinates(
    const symbolic::CyclotomicFieldContext& field,
    std::vector<Rational> base,
    std::uint64_t exponent) {
    std::vector<Rational> result(field.degree());
    result[0] = Rational{BigInt{1}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0U)
            result = field.multiply(result, base);
        exponent >>= 1U;
        if (exponent != 0)
            base = field.multiply(base, base);
    }
    return result;
}

[[nodiscard]] std::optional<std::vector<Rational>> cyclotomicCoordinatesImpl(
    const Expr& expression,
    const symbolic::CyclotomicFieldContext& field,
    const Expr& generator,
    const evaluation::BuiltinRegistry& registry,
    std::size_t& remainingNodes) {
    if (remainingNodes == 0)
        return std::nullopt;
    --remainingNodes;

    if (expression == generator)
        return std::vector<Rational>(field.power(1).begin(), field.power(1).end());

    if (expression.isNumber()) {
        const Number& number = expression.asNumber();
        if (number.isReal())
            return field.embedGaussianRational(number.asReal().toRational(), Rational{});
        return field.embedGaussianRational(
            number.asComplex().real.toRational(),
            number.asComplex().imaginary.toRational());
    }
    if (!expression.isCall())
        return std::nullopt;

    const auto& call = expression.asCall();
    const auto* definition = registry.find(call.head);
    if (!definition)
        return std::nullopt;
    const auto child = [&](std::size_t index) {
        return cyclotomicCoordinatesImpl(
            call.arguments[index], field, generator, registry, remainingNodes);
    };

    switch (definition->id) {
    case BuiltinId::Negate: {
        if (call.arguments.size() != 1)
            return std::nullopt;
        auto value = child(0);
        return value ? std::optional<std::vector<Rational>>{negateCoordinates(*value)}
                     : std::nullopt;
    }
    case BuiltinId::Add: {
        std::vector<Rational> result(field.degree());
        for (std::size_t i = 0; i < call.arguments.size(); ++i) {
            auto value = child(i);
            if (!value)
                return std::nullopt;
            result = addCoordinates(result, *value);
        }
        return result;
    }
    case BuiltinId::Subtract: {
        if (call.arguments.size() != 2)
            return std::nullopt;
        auto lhs = child(0);
        auto rhs = child(1);
        if (!lhs || !rhs)
            return std::nullopt;
        return addCoordinates(*lhs, negateCoordinates(*rhs));
    }
    case BuiltinId::Multiply: {
        std::vector<Rational> result(field.degree());
        result[0] = Rational{BigInt{1}};
        for (std::size_t i = 0; i < call.arguments.size(); ++i) {
            auto value = child(i);
            if (!value)
                return std::nullopt;
            result = field.multiply(result, *value);
        }
        return result;
    }
    case BuiltinId::Power: {
        if (call.arguments.size() != 2
            || !call.arguments[1].isNumber()
            || !call.arguments[1].asNumber().isReal()
            || !call.arguments[1].asNumber().asReal().isInteger())
            return std::nullopt;
        const BigInt& exponent = call.arguments[1].asNumber().asReal().asInteger();
        if (exponent.isNegative())
            return std::nullopt;
        const auto magnitude = numeric::tryToUint64(exponent);
        if (!magnitude || *magnitude > field.conductor())
            return std::nullopt;
        auto base = child(0);
        if (!base)
            return std::nullopt;
        return powerCoordinates(field, std::move(*base), *magnitude);
    }
    default:
        return std::nullopt;
    }
}

[[nodiscard]] std::optional<std::vector<Rational>> cyclotomicCoordinates(
    const Expr& expression,
    const symbolic::CyclotomicFieldContext& field,
    const Expr& generator,
    const evaluation::BuiltinRegistry& registry) {
    std::size_t remainingNodes = 256;
    return cyclotomicCoordinatesImpl(
        expression, field, generator, registry, remainingNodes);
}

[[nodiscard]] Expr cyclotomicPolynomialExpr(
    std::span<const Rational> coordinates,
    const symbolic::CyclotomicFieldContext& field,
    const Expr& generator,
    const evaluation::BuiltinRegistry& registry) {
    if (const auto gaussian = field.exactGaussianRational(coordinates)) {
        if (gaussian->second.isZero())
            return Expr{Number{gaussian->first}};
        return Expr{Number::complex(
            numeric::RealNumber{gaussian->first}, numeric::RealNumber{gaussian->second})};
    }

    std::vector<Expr> terms;
    for (std::size_t exponent = 0; exponent < coordinates.size(); ++exponent) {
        const Rational& coefficient = coordinates[exponent];
        if (coefficient.isZero())
            continue;
        if (exponent == 0) {
            terms.push_back(Expr{Number{coefficient}});
            continue;
        }

        Expr power = exponent == 1
            ? generator
            : Expr::call(registry.symbol(BuiltinId::Power),
                {generator, Expr{Number{BigInt::fromUnsigned(exponent)}}});
        if (coefficient == Rational{BigInt{1}})
            terms.push_back(std::move(power));
        else if (coefficient == Rational{BigInt{-1}})
            terms.push_back(Expr::call(
                registry.symbol(BuiltinId::Negate), {std::move(power)}));
        else
            terms.push_back(Expr::call(
                registry.symbol(BuiltinId::Multiply),
                {Expr{Number{coefficient}}, std::move(power)}));
    }

    if (terms.empty())
        return zero();
    if (terms.size() == 1)
        return std::move(terms.front());
    return Expr::call(registry.symbol(BuiltinId::Add), std::move(terms));
}

[[nodiscard]] std::size_t multiplyModulo(
    std::size_t lhs,
    std::size_t rhs,
    std::size_t modulus) noexcept {
    if (modulus == 0)
        return 0;
    lhs %= modulus;
    std::size_t result = 0;
    while (rhs != 0) {
        if ((rhs & 1U) != 0U)
            result = result >= modulus - lhs ? result - (modulus - lhs) : result + lhs;
        rhs >>= 1U;
        if (rhs == 0)
            break;
        lhs = lhs >= modulus - lhs ? lhs - (modulus - lhs) : lhs + lhs;
    }
    return result;
}

[[nodiscard]] std::optional<std::vector<std::vector<Rational>>> coordinatesForField(
    const std::vector<Expr>& input,
    const symbolic::CyclotomicFieldContext& field,
    const Expr& generator,
    const evaluation::BuiltinRegistry& registry) {
    std::vector<std::vector<Rational>> values;
    values.reserve(input.size());
    for (const Expr& expression : input) {
        auto coordinates = cyclotomicCoordinates(expression, field, generator, registry);
        if (!coordinates)
            return std::nullopt;
        values.push_back(std::move(*coordinates));
    }
    return values;
}

[[nodiscard]] std::shared_ptr<const symbolic::CyclotomicFieldContext> cachedOrCreateCyclotomicField(
    std::size_t conductor,
    FourierTransformCache& cache) {
    if (symbolic::cyclotomicDegree(conductor) > 64)
        return {};
    if (auto existing = cache.cyclotomicField(conductor))
        return existing;
    auto created = symbolic::CyclotomicFieldContext::create(conductor);
    if (created)
        cache.rememberCyclotomicField(conductor, created);
    return created;
}

[[nodiscard]] std::optional<std::vector<Expr>> exactCyclotomicTransform(
    const std::vector<Expr>& input,
    bool inverse,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    FourierTransformCache& cache) {
    const std::size_t n = input.size();
    if (n < 5 || isPowerOfTwo(n))
        return std::nullopt;
    if (n > std::numeric_limits<std::size_t>::max() / 4)
        return std::nullopt;

    bool allNumbers = true;
    bool hasImaginaryRational = false;
    for (const Expr& expression : input) {
        if (!expression.isNumber()) {
            allNumbers = false;
            continue;
        }
        if (!expression.asNumber().isReal()
            && !expression.asNumber().asComplex().imaginary.isZero())
            hasImaginaryRational = true;
    }

    std::vector<std::size_t> conductors;
    const std::size_t gaussianConductor = std::lcm(n, std::size_t{4});
    if (allNumbers)
        conductors.push_back(hasImaginaryRational ? gaussianConductor : n);
    else {
        for (const std::size_t conductor : {n, gaussianConductor})
            if (std::find(conductors.begin(), conductors.end(), conductor) == conductors.end())
                conductors.push_back(conductor);
    }

    std::shared_ptr<const symbolic::CyclotomicFieldContext> cyclotomic;
    std::optional<std::vector<std::vector<Rational>>> values;
    std::optional<Expr> generator;
    for (const std::size_t conductor : conductors) {
        auto candidate = cachedOrCreateCyclotomicField(conductor, cache);
        if (!candidate)
            continue;
        Expr candidateGenerator = cyclotomicGeneratorExpr(
            *candidate, registry, mathematics, angles);
        auto candidateValues = coordinatesForField(
            input, *candidate, candidateGenerator, registry);
        if (!candidateValues)
            continue;
        cyclotomic = std::move(candidate);
        values = std::move(candidateValues);
        generator = std::move(candidateGenerator);
        break;
    }
    if (!cyclotomic || !values || !generator)
        return std::nullopt;

    const std::size_t conductor = cyclotomic->conductor();
    if (conductor % n != 0)
        return std::nullopt;
    const std::size_t degree = cyclotomic->degree();
    const std::size_t rootStep = conductor / n;
    std::vector<std::vector<Rational>> transformed(
        n, std::vector<Rational>(degree));
    for (std::size_t k = 0; k < n; ++k) {
        auto& sum = transformed[k];
        for (std::size_t j = 0; j < n; ++j) {
            const std::size_t jk = multiplyModulo(j, k, n);
            std::size_t exponent = rootStep * jk;
            if (inverse && exponent != 0)
                exponent = conductor - exponent;
            const auto term = cyclotomic->multiply(
                (*values)[j], cyclotomic->power(exponent));
            for (std::size_t i = 0; i < degree; ++i)
                sum[i] += term[i];
        }
        if (inverse) {
            const Rational scale{BigInt{1}, sizeInteger(n)};
            for (Rational& coefficient : sum)
                coefficient *= scale;
        }
    }

    std::vector<Expr> output;
    output.reserve(n);
    for (const auto& coordinates : transformed)
        output.push_back(cyclotomicPolynomialExpr(coordinates, *cyclotomic, *generator, registry));
    return output;
}

/*
旧 exact FFT dispatch（v1.5.3のcyclotomic backend導入前）。

    return vectorExpr(radix2Transform(
        input, inverse, registry, mathematics, angles, cache));

非2冪長ではradix2Transform()がdirectTransform()へfallbackし，twiddleを
cis[2 Pi k/n]のgeneric Exprとして構築していた。この方式はsymbolic入力のfallbackとして
現在も残すが，exact Rational / Gaussian Rational と同一cyclotomic field由来の入力では，
Q(zeta_n)のpower-basis座標演算へ先に落とす。変更理由は，5/7/10/12点等で同じ
root-of-unity恒等式をSimplifierへ何度も再証明させ，ifft[fft[...]]が巨大式になるためである。
*/

[[nodiscard]] Expr vectorExpr(std::vector<Expr> elements);

[[nodiscard]] approximation::ComplexInterval exactComplexRational(
    const Rational& real,
    const Rational& imaginary,
    std::size_t precisionBits) {
    return approximation::ComplexInterval{
        approximation::RealInterval::fromRational(real, precisionBits),
        approximation::RealInterval::fromRational(imaginary, precisionBits)};
}

[[nodiscard]] std::optional<approximation::ComplexInterval> exactQuarterTurnRoot(
    std::size_t exponent,
    std::size_t length,
    bool inverse,
    std::size_t precisionBits) {
    if (length == 0 || exponent > std::numeric_limits<std::size_t>::max() / 4)
        return std::nullopt;
    const std::size_t scaled = exponent * 4;
    if (scaled % length != 0)
        return std::nullopt;

    std::size_t quarter = (scaled / length) & 3U;
    if (!inverse && quarter != 0)
        quarter = (4 - quarter) & 3U;
    switch (quarter) {
    case 0: return exactComplexRational(Rational{BigInt{1}}, Rational{}, precisionBits);
    case 1: return exactComplexRational(Rational{}, Rational{BigInt{1}}, precisionBits);
    case 2: return exactComplexRational(Rational{BigInt{-1}}, Rational{}, precisionBits);
    case 3: return exactComplexRational(Rational{}, Rational{BigInt{-1}}, precisionBits);
    default: return std::nullopt;
    }
}

[[nodiscard]] approximation::ComplexInterval approximatePrimitiveRoot(
    std::size_t length,
    bool inverse,
    std::size_t precisionBits) {
    // radix-2の最初の二段は根が有理複素数なので、Pi/trig区間へ落とさずexact pointを使う。
    // これによりN[fft[{1,2,3,4}],p]でも整数・Gaussian integer成分を従来どおり最小表記で保てる。
    if (length == 2)
        return exactComplexRational(Rational{BigInt{-1}}, Rational{}, precisionBits);
    if (length == 4)
        return exactComplexRational(
            Rational{}, Rational{BigInt{inverse ? 1 : -1}}, precisionBits);

    Rational turns{BigInt{1}, sizeInteger(length)};
    if (!inverse)
        turns = -turns;
    return approximation::ComplexInterval{
        approximation::encloseCosTurns(turns, precisionBits).interval,
        approximation::encloseSinTurns(turns, precisionBits).interval};
}

[[nodiscard]] std::vector<approximation::ComplexInterval> approximateDirectTransform(
    const std::vector<approximation::ComplexInterval>& input,
    bool inverse,
    std::size_t precisionBits) {
    const std::size_t n = input.size();
    if (n == 0)
        return {};

    const auto one = exactComplexRational(Rational{BigInt{1}}, Rational{}, precisionBits);
    const auto primitive = approximatePrimitiveRoot(n, inverse, precisionBits);
    std::vector<approximation::ComplexInterval> roots;
    roots.reserve(n);
    roots.push_back(one);
    for (std::size_t i = 1; i < n; ++i) {
        if (const auto exact = exactQuarterTurnRoot(i, n, inverse, precisionBits))
            roots.push_back(*exact);
        else
            roots.push_back(approximation::multiply(roots.back(), primitive, precisionBits));
    }

    std::vector<approximation::ComplexInterval> output;
    output.reserve(n);
    for (std::size_t k = 0; k < n; ++k) {
        approximation::ComplexInterval sum = exactComplexRational(Rational{}, Rational{}, precisionBits);
        std::size_t rootIndex = 0;
        for (std::size_t j = 0; j < n; ++j) {
            sum = approximation::add(sum,
                approximation::multiply(input[j], roots[rootIndex], precisionBits),
                precisionBits);
            if (j + 1 < n && k != 0) {
                if (rootIndex >= n - k)
                    rootIndex -= n - k;
                else
                    rootIndex += k;
            }
        }
        if (inverse) {
            const auto scale = approximation::ComplexInterval::fromReal(
                approximation::RealInterval::fromRational(
                    Rational{BigInt{1}, sizeInteger(n)}, precisionBits));
            sum = approximation::multiply(sum, scale, precisionBits);
        }
        output.push_back(std::move(sum));
    }
    return output;
}

[[nodiscard]] std::vector<std::size_t> bitReversedIndices(std::size_t n) {
    std::vector<std::size_t> result(n);
    for (std::size_t i = 0, j = 0; i < n; ++i) {
        result[i] = j;
        if (i + 1 == n)
            break;
        std::size_t bit = n >> 1;
        while ((j & bit) != 0) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
    }
    return result;
}

[[nodiscard]] std::vector<approximation::ComplexInterval> approximateRadix2PowerOfTwo(
    const std::vector<approximation::ComplexInterval>& input,
    bool inverse,
    std::size_t precisionBits) {
    const std::size_t n = input.size();
    if (n == 0)
        return {};
    if (!isPowerOfTwo(n))
        throw std::invalid_argument("Radix-2 Fourier backend requires a power-of-two length");

    const std::vector<std::size_t> bitReversed = bitReversedIndices(n);
    const auto zero = exactComplexRational(Rational{}, Rational{}, precisionBits);
    std::vector<approximation::ComplexInterval> data(n, zero);
    for (std::size_t i = 0; i < n; ++i)
        data[bitReversed[i]] = input[i];

    for (std::size_t length = 2; length <= n; length <<= 1) {
        const std::size_t half = length >> 1;
        const auto primitive = approximatePrimitiveRoot(length, inverse, precisionBits);
        std::vector<approximation::ComplexInterval> roots;
        roots.reserve(half);
        roots.push_back(exactComplexRational(Rational{BigInt{1}}, Rational{}, precisionBits));
        for (std::size_t j = 1; j < half; ++j) {
            if (const auto exact = exactQuarterTurnRoot(j, length, inverse, precisionBits))
                roots.push_back(*exact);
            else
                roots.push_back(approximation::multiply(roots.back(), primitive, precisionBits));
        }

        for (std::size_t block = 0; block < n; block += length) {
            for (std::size_t j = 0; j < half; ++j) {
                const auto even = data[block + j];
                const auto odd = approximation::multiply(
                    data[block + j + half], roots[j], precisionBits);
                data[block + j] = approximation::add(even, odd, precisionBits);
                data[block + j + half] = approximation::subtract(even, odd, precisionBits);
            }
        }
        if (length == n)
            break;
    }

    if (inverse) {
        const auto scale = approximation::ComplexInterval::fromReal(
            approximation::RealInterval::fromRational(
                Rational{BigInt{1}, sizeInteger(n)}, precisionBits));
        for (auto& value : data)
            value = approximation::multiply(value, scale, precisionBits);
    }
    return data;
}

[[nodiscard]] approximation::ComplexInterval conjugate(
    const approximation::ComplexInterval& value) {
    return approximation::ComplexInterval{
        value.real(), approximation::negate(value.imaginary())};
}

[[nodiscard]] approximation::ComplexInterval approximateTurnRoot(
    const Rational& turns,
    std::size_t precisionBits) {
    return approximation::ComplexInterval{
        approximation::encloseCosTurns(turns, precisionBits).interval,
        approximation::encloseSinTurns(turns, precisionBits).interval};
}

[[nodiscard]] std::vector<approximation::ComplexInterval> approximateChirp(
    std::size_t n,
    bool inverse,
    std::size_t precisionBits) {
    if (n == 0)
        return {};

    // c_k = exp(sign*pi*i*k^2/n)。各kでtrigを再評価せず、
    // c_{k+1}/c_k = exp(sign*pi*i*(2k+1)/n) をさらに一定比で更新する。
    // これによりBluestein用chirp生成のtranscendental評価は2回だけで済む。
    BigInt sign{inverse ? 1 : -1};
    const BigInt denominator = sizeInteger(n);
    const Rational initialTurns{sign, denominator * BigInt{2}};
    const Rational ratioStepTurns{sign, denominator};
    approximation::ComplexInterval ratio = approximateTurnRoot(initialTurns, precisionBits);
    const approximation::ComplexInterval ratioStep = approximateTurnRoot(
        ratioStepTurns, precisionBits);

    std::vector<approximation::ComplexInterval> chirp;
    chirp.reserve(n);
    chirp.push_back(exactComplexRational(Rational{BigInt{1}}, Rational{}, precisionBits));
    for (std::size_t k = 1; k < n; ++k) {
        chirp.push_back(approximation::multiply(chirp.back(), ratio, precisionBits));
        ratio = approximation::multiply(ratio, ratioStep, precisionBits);
    }
    return chirp;
}

[[nodiscard]] std::size_t convolutionLength(std::size_t n) {
    if (n == 0)
        return 0;
    const std::size_t maximumInput = std::numeric_limits<std::size_t>::max() / 2 + 1;
    if (n > maximumInput)
        throw std::overflow_error("Fourier transform size is too large");
    const std::size_t required = n * 2 - 1;
    std::size_t length = 1;
    while (length < required) {
        if (length > std::numeric_limits<std::size_t>::max() / 2)
            throw std::overflow_error("Fourier convolution size is too large");
        length <<= 1;
    }
    return length;
}

[[nodiscard]] std::vector<approximation::ComplexInterval> approximateBluesteinTransform(
    const std::vector<approximation::ComplexInterval>& input,
    bool inverse,
    std::size_t precisionBits) {
    const std::size_t n = input.size();
    if (n == 0)
        return {};

    const std::size_t m = convolutionLength(n);
    const auto zero = exactComplexRational(Rational{}, Rational{}, precisionBits);
    const std::vector<approximation::ComplexInterval> chirp = approximateChirp(
        n, inverse, precisionBits);

    std::vector<approximation::ComplexInterval> a(m, zero);
    std::vector<approximation::ComplexInterval> b(m, zero);
    for (std::size_t k = 0; k < n; ++k) {
        a[k] = approximation::multiply(input[k], chirp[k], precisionBits);
        const auto opposite = conjugate(chirp[k]);
        b[k] = opposite;
        if (k != 0)
            b[m - k] = opposite;
    }

    auto spectrumA = approximateRadix2PowerOfTwo(a, false, precisionBits);
    auto spectrumB = approximateRadix2PowerOfTwo(b, false, precisionBits);
    for (std::size_t k = 0; k < m; ++k)
        spectrumA[k] = approximation::multiply(
            spectrumA[k], spectrumB[k], precisionBits);
    auto convolution = approximateRadix2PowerOfTwo(spectrumA, true, precisionBits);

    std::vector<approximation::ComplexInterval> output;
    output.reserve(n);
    for (std::size_t k = 0; k < n; ++k)
        output.push_back(approximation::multiply(
            convolution[k], chirp[k], precisionBits));

    if (inverse) {
        const auto scale = approximation::ComplexInterval::fromReal(
            approximation::RealInterval::fromRational(
                Rational{BigInt{1}, sizeInteger(n)}, precisionBits));
        for (auto& value : output)
            value = approximation::multiply(value, scale, precisionBits);
    }
    return output;
}

[[nodiscard]] std::vector<approximation::ComplexInterval> approximateFastTransform(
    const std::vector<approximation::ComplexInterval>& input,
    bool inverse,
    std::size_t precisionBits) {
    if (isPowerOfTwo(input.size()))
        return approximateRadix2PowerOfTwo(input, inverse, precisionBits);

    // ごく小さい非2冪は直接DFTの方が定数項が小さい。それ以上はBluesteinで
    // O(N^2) fallbackを避け、次の2冪長のconvolutionへ還元する。
    // GCC Releaseの強制比較では319点でdirect、335点でBluesteinが僅差で逆転した。
    // MSVC実測では257点でdirect、509点でBluesteinが優位だったため、primary環境へ
    // 保守的に寄せた384点をpolicy境界とする。--fft-thresholdで再測定可能。
    if (input.size() < approximateFftBluesteinThreshold)
        return approximateDirectTransform(input, inverse, precisionBits);
    return approximateBluesteinTransform(input, inverse, precisionBits);
}

template <class Transform>
[[nodiscard]] std::optional<Expr> approximateTransform(
    std::span<const Expr> arguments,
    std::string_view name,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context,
    Transform&& transform) {
    requireArity(arguments, 1, name);
    const std::vector<Expr> inputExpressions = vectorArgument(arguments.front(), name);
    approximation::CertifiedEvaluator certified{registry, mathematics, angles};

    constexpr std::size_t maximumRefinements = 12;
    for (std::size_t refinement = 0; refinement < maximumRefinements; ++refinement) {
        evaluation::consumeEvaluationBudget(evaluation::EvaluationResource::CertifiedRefinement);
        const std::size_t bits = context.workingBinaryBits();
        try {
            std::vector<approximation::ComplexInterval> input;
            std::vector<approximation::ComplexInterval> informationInput;
            input.reserve(inputExpressions.size());
            informationInput.reserve(inputExpressions.size());
            for (const Expr& expression : inputExpressions) {
                const auto information = approximation::encloseComplexExpression(
                    expression, bits, certified,
                    approximation::CertifiedEvaluator::EnclosureKind::Information);
                const auto value = approximation::encloseComplexExpression(
                    expression, bits, certified,
                    approximation::CertifiedEvaluator::EnclosureKind::Certified);
                if (!value || !information)
                    return std::nullopt;
                input.push_back(*value);
                informationInput.push_back(*information);
            }

            const auto transformed = transform(input, bits);
            const auto informationTransformed = transform(informationInput, bits);
            if (transformed.size() != informationTransformed.size())
                return std::nullopt;
            std::vector<Expr> output;
            output.reserve(transformed.size());
            bool rounded = true;
            for (std::size_t i = 0; i < transformed.size(); ++i) {
                const auto decimal = approximation::finalizeCertifiedApproximation(
                    approximation::CertifiedValue{transformed[i]},
                    approximation::CertifiedValue{informationTransformed[i]},
                    context.decimalDigits());
                if (!decimal) {
                    rounded = false;
                    break;
                }
                output.push_back(*decimal);
            }
            if (rounded)
                return vectorExpr(std::move(output));
        }
        catch (const approximation::PrecisionInsufficient&) {
            // 現作業精度では象限や丸めを証明できない。guardを増やして同じ式を再評価する。
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            return std::nullopt;
        }

        context.setGuardDigits(approximation::nextGuardDigits(context.guardDigits()));
    }
    return std::nullopt;
}

[[nodiscard]] Expr vectorExpr(std::vector<Expr> elements) {
    const std::size_t count = elements.size();
    return Expr::array({count}, std::move(elements));
}

} // namespace

void FourierTransformCache::clear() noexcept {
    plans_.clear();
    cyclotomicFields_.clear();
}

std::size_t FourierTransformCache::planCount() const noexcept {
    return plans_.size();
}

std::size_t FourierTransformCache::cyclotomicFieldCount() const noexcept {
    return cyclotomicFields_.size();
}

FourierTransformCache::Plan& FourierTransformCache::plan(std::size_t size) {
    return plans_[size];
}

std::shared_ptr<const symbolic::CyclotomicFieldContext> FourierTransformCache::cyclotomicField(
    std::size_t conductor) const {
    const auto iterator = cyclotomicFields_.find(conductor);
    return iterator == cyclotomicFields_.end() ? nullptr : iterator->second;
}

void FourierTransformCache::rememberCyclotomicField(
    std::size_t conductor,
    std::shared_ptr<const symbolic::CyclotomicFieldContext> field) {
    cyclotomicFields_[conductor] = std::move(field);
}

Expr evaluateDft(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::dft);
    const std::vector<Expr> input = vectorArgument(arguments.front(), names::dft);
    if (const auto context = approximation::inferredApproximationContext(input))
        if (const auto result = evaluateApproximateDft(
            arguments, registry, mathematics, angles, *context))
            return *result;
    return vectorExpr(directTransform(input, false, registry, mathematics, angles));
}

Expr evaluateFft(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    FourierTransformCache& cache) {
    requireArity(arguments, 1, names::fft);
    const std::vector<Expr> input = vectorArgument(arguments.front(), names::fft);
    if (const auto context = approximation::inferredApproximationContext(input))
        if (const auto result = evaluateApproximateFft(
            arguments, registry, mathematics, angles, *context))
            return *result;
    if (const auto cyclotomic = exactCyclotomicTransform(
            input, false, registry, mathematics, angles, cache))
        return vectorExpr(*cyclotomic);
    return vectorExpr(radix2Transform(
        input, false, registry, mathematics, angles, cache));
}

Expr evaluateIfft(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    FourierTransformCache& cache) {
    requireArity(arguments, 1, names::ifft);
    const std::vector<Expr> input = vectorArgument(arguments.front(), names::ifft);
    if (const auto context = approximation::inferredApproximationContext(input))
        if (const auto result = evaluateApproximateIfft(
            arguments, registry, mathematics, angles, *context))
            return *result;
    if (const auto cyclotomic = exactCyclotomicTransform(
            input, true, registry, mathematics, angles, cache))
        return vectorExpr(*cyclotomic);
    return vectorExpr(radix2Transform(
        input, true, registry, mathematics, angles, cache));
}

std::optional<Expr> evaluateApproximateDft(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return approximateTransform(arguments, names::dft, registry, mathematics, angles,
        std::move(context), [](const auto& input, std::size_t bits) {
            return approximateDirectTransform(input, false, bits);
        });
}

std::optional<Expr> evaluateApproximateFft(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return approximateTransform(arguments, names::fft, registry, mathematics, angles,
        std::move(context), [](const auto& input, std::size_t bits) {
            return approximateFastTransform(input, false, bits);
        });
}

std::optional<Expr> evaluateApproximateBluesteinFftForBenchmark(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return approximateTransform(arguments, names::fft, registry, mathematics, angles,
        std::move(context), [](const auto& input, std::size_t bits) {
            return approximateBluesteinTransform(input, false, bits);
        });
}

std::optional<Expr> evaluateApproximateIfft(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return approximateTransform(arguments, names::ifft, registry, mathematics, angles,
        std::move(context), [](const auto& input, std::size_t bits) {
            return approximateFastTransform(input, true, bits);
        });
}

Expr evaluateConvolution(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 2, names::convolve);
    const std::vector<Expr> lhs = vectorArgument(arguments[0], names::convolve);
    const std::vector<Expr> rhs = vectorArgument(arguments[1], names::convolve);
    if (lhs.empty() || rhs.empty())
        return vectorExpr({});

    if (lhs.size() > std::numeric_limits<std::size_t>::max() - rhs.size() + 1)
        error::throwCalcError(error::CalcErrorType::Overflow, "convolution result is too large");

    std::vector<Expr> output(lhs.size() + rhs.size() - 1, zero());
    for (std::size_t i = 0; i < lhs.size(); ++i)
        for (std::size_t j = 0; j < rhs.size(); ++j)
            output[i + j] = add(std::move(output[i + j]),
                multiply(lhs[i], rhs[j], registry, mathematics, angles),
                registry, mathematics, angles);
    return vectorExpr(std::move(output));
}

} // namespace mmcal::builtins
