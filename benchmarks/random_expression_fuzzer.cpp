#include "random_expression_fuzzer.hpp"

#include "approximation/approximation_context.hpp"
#include "approximation/certified_evaluator.hpp"
#include "formatting/expr_formatter.hpp"
#include "error/error_message.hpp"
#include "kernel/kernel_session.hpp"
#include "mathematics/definedness.hpp"
#include "mathematics/knowledge_context.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "solver/solution_set.hpp"

#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace mmcal::benchmarks {
namespace {

using Clock = std::chrono::steady_clock;
using mmcal::expression::Expr;

struct GeneratedExpr final {
    std::string source;
    std::vector<GeneratedExpr> children;
};

enum class Invariant {
    ExactRoundTrip,
    FullSimplifyPreservesValue,
    ExpandPreservesPolynomial,
    FactorPreservesPolynomial,
    DoubleTranspose,
    DeterminantTranspose,
    DerivativeOfIntegral,
    SolveKnownRoots,
    InverseReconstruction,
    FourierRoundTrip,
    DomainAndBranchBoundary,
    PolynomialLimitMatchesSubstitution,
    CasesShortCircuit,
    NestedNDoesNotInventPrecision,
    AlgebraicRootIdentity,
    GroebnerGeneratorsReduceToZero,
    ArrayReshapeRoundTrip,
    SimplifyPreservesDomainHole,
    RemovableLimitPreservesDomainHole
};

struct GeneratedCase final {
    Invariant invariant = Invariant::ExactRoundTrip;
    GeneratedExpr expression;
    std::size_t targetDepth = 1;
    std::int64_t substitution = 0;
    std::int64_t firstRoot = 0;
    std::int64_t secondRoot = 1;
    std::size_t boundaryCase = 0;
    std::size_t precisionLow = 20;
    std::size_t precisionHigh = 40;
    std::size_t rows = 1;
    std::size_t columns = 1;
    std::string termOrder{"GrevLex"};
};

struct Failure final {
    std::string reason;
    std::string expected;
    std::string actual;
    std::string budget{};
    bool inconclusive = false;
};

[[nodiscard]] std::uint64_t splitMix64(std::uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

[[nodiscard]] std::uint64_t caseSeed(std::uint64_t seed, std::uint64_t caseIndex) noexcept {
    return splitMix64(seed ^ splitMix64(caseIndex));
}

[[nodiscard]] std::size_t randomBetween(
    std::mt19937_64& rng,
    std::size_t lower,
    std::size_t upper) {
    if (upper <= lower)
        return lower;
    std::uniform_int_distribution<std::size_t> distribution(lower, upper);
    return distribution(rng);
}

[[nodiscard]] std::size_t chooseTargetDepth(std::mt19937_64& rng, std::size_t maxDepth) {
    if (maxDepth <= 1)
        return 1;

    std::uniform_int_distribution<unsigned> roll(0, 99);
    const unsigned bucket = roll(rng);

    std::size_t low = 1;
    std::size_t high = std::min<std::size_t>(4, maxDepth);
    if (bucket >= 60 && bucket < 85) {
        low = std::min<std::size_t>(5, maxDepth);
        high = std::min<std::size_t>(7, maxDepth);
    }
    else if (bucket >= 85 && bucket < 95) {
        low = std::min<std::size_t>(8, maxDepth);
        high = std::min<std::size_t>(10, maxDepth);
    }
    else if (bucket >= 95 && bucket < 99) {
        low = std::min<std::size_t>(11, maxDepth);
        high = std::min<std::size_t>(13, maxDepth);
    }
    else if (bucket >= 99) {
        low = std::min<std::size_t>(14, maxDepth);
        high = maxDepth;
    }

    if (low > high)
        low = high;
    return randomBetween(rng, low, high);
}

[[nodiscard]] GeneratedExpr leaf(std::mt19937_64& rng) {
    static constexpr std::array<std::string_view, 12> leaves{
        "0", "1", "-1", "2", "-2", "3", "5", "7", "11",
        "1/2", "-1/3", "2/5"};
    return GeneratedExpr{std::string{leaves[randomBetween(rng, 0, leaves.size() - 1)]}, {}};
}

[[nodiscard]] GeneratedExpr generateExactScalar(
    std::mt19937_64& rng,
    std::size_t depth,
    std::size_t targetDepth) {
    if (depth >= targetDepth)
        return leaf(rng);

    // 深くなるほどleaf率を上げ，target depthの全枝展開によるAST爆発を避ける。
    const unsigned leafPercent = static_cast<unsigned>(15 + 55 * depth / std::max<std::size_t>(1, targetDepth));
    std::uniform_int_distribution<unsigned> percent(0, 99);
    if (depth != 0 && percent(rng) < leafPercent)
        return leaf(rng);

    const unsigned operation = static_cast<unsigned>(randomBetween(rng, 0, 10));
    if (operation <= 4) {
        auto lhs = generateExactScalar(rng, depth + 1, targetDepth);
        auto rhs = generateExactScalar(rng, depth + 1, targetDepth);
        static constexpr std::array<std::string_view, 5> operators{"+", "-", "*", "+", "*"};
        const std::string_view op = operators[operation];
        return GeneratedExpr{
            "(" + lhs.source + std::string{op} + rhs.source + ")",
            {std::move(lhs), std::move(rhs)}};
    }

    if (operation == 5) {
        auto child = generateExactScalar(rng, depth + 1, targetDepth);
        const unsigned exponent = static_cast<unsigned>(randomBetween(rng, 1, 5));
        return GeneratedExpr{
            "(" + child.source + ")^" + std::to_string(exponent),
            {std::move(child)}};
    }

    if (operation == 6) {
        // radicandを平方にしてdomain errorを避ける。
        auto child = generateExactScalar(rng, depth + 1, targetDepth);
        return GeneratedExpr{
            "sqrt[(" + child.source + ")^2]",
            {std::move(child)}};
    }

    if (operation == 7) {
        auto child = generateExactScalar(rng, depth + 1, targetDepth);
        return GeneratedExpr{"abs[" + child.source + "]", {std::move(child)}};
    }

    if (operation == 8) {
        static constexpr std::array<std::string_view, 8> trig{
            "sin[0]", "cos[0]", "sin[Pi/6]", "cos[Pi/3]",
            "sin[Pi/2]", "cos[Pi]", "tan[0]", "exp[0]"};
        return GeneratedExpr{std::string{trig[randomBetween(rng, 0, trig.size() - 1)]}, {}};
    }

    if (operation == 9) {
        auto child = generateExactScalar(rng, depth + 1, targetDepth);
        return GeneratedExpr{"fullSimplify[" + child.source + "]", {std::move(child)}};
    }

    auto child = generateExactScalar(rng, depth + 1, targetDepth);
    return GeneratedExpr{"simplify[" + child.source + "]", {std::move(child)}};
}

[[nodiscard]] GeneratedExpr polynomialLeaf(std::mt19937_64& rng) {
    static constexpr std::array<std::string_view, 9> leaves{
        "x", "0", "1", "-1", "2", "-2", "3", "5", "7"};
    return GeneratedExpr{std::string{leaves[randomBetween(rng, 0, leaves.size() - 1)]}, {}};
}

[[nodiscard]] GeneratedExpr generatePolynomial(
    std::mt19937_64& rng,
    std::size_t depth,
    std::size_t targetDepth) {
    if (depth >= targetDepth)
        return polynomialLeaf(rng);

    std::uniform_int_distribution<unsigned> percent(0, 99);
    const unsigned leafPercent = static_cast<unsigned>(20 + 60 * depth / std::max<std::size_t>(1, targetDepth));
    if (depth != 0 && percent(rng) < leafPercent)
        return polynomialLeaf(rng);

    const unsigned operation = static_cast<unsigned>(randomBetween(rng, 0, 4));
    if (operation <= 2) {
        auto lhs = generatePolynomial(rng, depth + 1, targetDepth);
        auto rhs = generatePolynomial(rng, depth + 1, targetDepth);
        static constexpr std::array<std::string_view, 3> operators{"+", "-", "*"};
        return GeneratedExpr{
            "(" + lhs.source + std::string{operators[operation]} + rhs.source + ")",
            {std::move(lhs), std::move(rhs)}};
    }

    auto child = generatePolynomial(rng, depth + 1, targetDepth);
    const unsigned exponent = static_cast<unsigned>(randomBetween(rng, 1, 4));
    return GeneratedExpr{
        "(" + child.source + ")^" + std::to_string(exponent),
        {std::move(child)}};
}

[[nodiscard]] GeneratedExpr generateMatrix(std::mt19937_64& rng) {
    const std::size_t rows = randomBetween(rng, 1, 4);
    const std::size_t cols = randomBetween(rng, 1, 4);
    std::string source = "{";
    std::vector<GeneratedExpr> elements;
    elements.reserve(rows * cols);
    for (std::size_t row = 0; row < rows; ++row) {
        if (row != 0)
            source += ',';
        source += '{';
        for (std::size_t col = 0; col < cols; ++col) {
            if (col != 0)
                source += ',';
            const std::int64_t value = static_cast<std::int64_t>(randomBetween(rng, 0, 10)) - 5;
            source += std::to_string(value);
            elements.push_back(GeneratedExpr{std::to_string(value), {}});
        }
        source += '}';
    }
    source += '}';
    return GeneratedExpr{std::move(source), std::move(elements)};
}

[[nodiscard]] GeneratedExpr generateInvertibleMatrix(
    std::mt19937_64& rng,
    std::size_t n) {
    // 非零対角の上三角整数行列に限定し，inverse側の失敗をgeneratorのsingular
    // blind spotと混同しない。非対角項は0以外も十分含める。
    std::string source = "{";
    std::vector<GeneratedExpr> elements;
    elements.reserve(n * n);
    for (std::size_t row = 0; row < n; ++row) {
        if (row != 0)
            source += ',';
        source += '{';
        for (std::size_t column = 0; column < n; ++column) {
            if (column != 0)
                source += ',';
            std::int64_t value = 0;
            if (column == row) {
                value = static_cast<std::int64_t>(randomBetween(rng, 1, 5));
                if (randomBetween(rng, 0, 1) == 0)
                    value = -value;
            }
            else if (column > row)
                value = static_cast<std::int64_t>(randomBetween(rng, 0, 8)) - 4;
            source += std::to_string(value);
            elements.push_back(GeneratedExpr{std::to_string(value), {}});
        }
        source += '}';
    }
    source += '}';
    return GeneratedExpr{std::move(source), std::move(elements)};
}

[[nodiscard]] std::string identityMatrixSource(std::size_t n) {
    std::string source = "{";
    for (std::size_t row = 0; row < n; ++row) {
        if (row != 0)
            source += ',';
        source += '{';
        for (std::size_t column = 0; column < n; ++column) {
            if (column != 0)
                source += ',';
            source += row == column ? '1' : '0';
        }
        source += '}';
    }
    source += '}';
    return source;
}

[[nodiscard]] GeneratedExpr generateFourierVector(std::mt19937_64& rng) {
    // semantic fuzzerは1 caseをboundedに保つ。11/12点級のexact non-power-of-two FFTは
    // cyclotomic式のfullSimplifyが秒単位になることがあり，専用FFT回帰/benchmarkで扱う。
    const std::size_t length = randomBetween(rng, 1, 10);
    std::string source = "{";
    std::vector<GeneratedExpr> elements;
    elements.reserve(length);
    for (std::size_t i = 0; i < length; ++i) {
        if (i != 0)
            source += ',';
        const std::int64_t real = static_cast<std::int64_t>(randomBetween(rng, 0, 10)) - 5;
        const std::int64_t imaginary = static_cast<std::int64_t>(randomBetween(rng, 0, 6)) - 3;
        std::string element = std::to_string(real);
        if (imaginary != 0)
            element += imaginary > 0
                ? "+" + std::to_string(imaginary) + "I"
                : std::to_string(imaginary) + "I";
        source += element;
        elements.push_back(GeneratedExpr{std::move(element), {}});
    }
    source += '}';
    return GeneratedExpr{std::move(source), std::move(elements)};
}

[[nodiscard]] GeneratedExpr generateReshapeVector(
    std::mt19937_64& rng,
    std::size_t rows,
    std::size_t columns) {
    const std::size_t length = rows * columns;
    std::string source = "{";
    std::vector<GeneratedExpr> elements;
    elements.reserve(length);
    for (std::size_t i = 0; i < length; ++i) {
        if (i != 0)
            source += ',';
        const std::int64_t value = static_cast<std::int64_t>(randomBetween(rng, 0, 20)) - 10;
        source += std::to_string(value);
        elements.push_back(GeneratedExpr{std::to_string(value), {}});
    }
    source += '}';
    return GeneratedExpr{std::move(source), std::move(elements)};
}

[[nodiscard]] std::string replaceSymbolX(std::string_view source, std::int64_t value) {
    // このfuzzerの多項式generatorが導入するsymbolはxだけである。
    // Formatterは125xのような暗黙乗算を出すため，数字をidentifier境界とは扱わない。
    const std::string replacement = "(" + std::to_string(value) + ")";
    std::string result;
    result.reserve(source.size() + 16);
    for (const char character : source) {
        if (character == 'x')
            result += replacement;
        else
            result.push_back(character);
    }
    return result;
}

[[nodiscard]] std::string invariantName(Invariant invariant) {
    switch (invariant) {
    case Invariant::ExactRoundTrip:
        return "exact-round-trip";
    case Invariant::FullSimplifyPreservesValue:
        return "fullSimplify-preserves-value";
    case Invariant::ExpandPreservesPolynomial:
        return "expand-preserves-polynomial";
    case Invariant::FactorPreservesPolynomial:
        return "factor-preserves-polynomial";
    case Invariant::DoubleTranspose:
        return "double-transpose";
    case Invariant::DeterminantTranspose:
        return "det-transpose";
    case Invariant::DerivativeOfIntegral:
        return "D-integrate-derivative-back";
    case Invariant::SolveKnownRoots:
        return "solve-known-roots";
    case Invariant::InverseReconstruction:
        return "matrix-inverse-reconstruction";
    case Invariant::FourierRoundTrip:
        return "fft-ifft-round-trip";
    case Invariant::DomainAndBranchBoundary:
        return "domain-branch-boundary";
    case Invariant::PolynomialLimitMatchesSubstitution:
        return "polynomial-limit-matches-substitution";
    case Invariant::CasesShortCircuit:
        return "cases-short-circuit";
    case Invariant::NestedNDoesNotInventPrecision:
        return "nested-N-does-not-invent-precision";
    case Invariant::AlgebraicRootIdentity:
        return "algebraic-root-identity";
    case Invariant::GroebnerGeneratorsReduceToZero:
        return "groebner-generators-reduce-to-zero";
    case Invariant::ArrayReshapeRoundTrip:
        return "array-reshape-round-trip";
    case Invariant::SimplifyPreservesDomainHole:
        return "simplify-preserves-domain-hole";
    case Invariant::RemovableLimitPreservesDomainHole:
        return "removable-limit-preserves-domain-hole";
    }
    return "unknown";
}

[[nodiscard]] GeneratedCase generateCase(
    std::uint64_t masterSeed,
    std::uint64_t index,
    std::size_t maxDepth) {
    std::mt19937_64 rng{caseSeed(masterSeed, index)};
    GeneratedCase result;
    result.targetDepth = chooseTargetDepth(rng, std::max<std::size_t>(1, maxDepth));
    result.substitution = static_cast<std::int64_t>(randomBetween(rng, 0, 10)) - 5;

    const unsigned choice = static_cast<unsigned>(randomBetween(rng, 0, 103));
    if (choice < 16) {
        result.invariant = Invariant::ExactRoundTrip;
        result.expression = generateExactScalar(rng, 0, result.targetDepth);
    }
    else if (choice < 29) {
        result.invariant = Invariant::FullSimplifyPreservesValue;
        result.expression = generateExactScalar(rng, 0, result.targetDepth);
    }
    else if (choice < 39) {
        result.invariant = Invariant::ExpandPreservesPolynomial;
        result.expression = generatePolynomial(rng, 0, std::min<std::size_t>(result.targetDepth, 6));
    }
    else if (choice < 47) {
        result.invariant = Invariant::FactorPreservesPolynomial;
        result.expression = generatePolynomial(rng, 0, std::min<std::size_t>(result.targetDepth, 5));
    }
    else if (choice < 52) {
        result.invariant = Invariant::DoubleTranspose;
        result.expression = generateMatrix(rng);
    }
    else if (choice < 56) {
        result.invariant = Invariant::DeterminantTranspose;
        const std::size_t n = randomBetween(rng, 1, 4);
        std::string source = "{";
        std::vector<GeneratedExpr> elements;
        for (std::size_t row = 0; row < n; ++row) {
            if (row != 0)
                source += ',';
            source += '{';
            for (std::size_t col = 0; col < n; ++col) {
                if (col != 0)
                    source += ',';
                const std::int64_t value = static_cast<std::int64_t>(randomBetween(rng, 0, 10)) - 5;
                source += std::to_string(value);
                elements.push_back(GeneratedExpr{std::to_string(value), {}});
            }
            source += '}';
        }
        source += '}';
        result.expression = GeneratedExpr{std::move(source), std::move(elements)};
    }
    else if (choice < 62) {
        result.invariant = Invariant::DerivativeOfIntegral;
        result.expression = generatePolynomial(
            rng, 0, std::min<std::size_t>(result.targetDepth, 4));
    }
    else if (choice < 67) {
        result.invariant = Invariant::SolveKnownRoots;
        result.firstRoot = static_cast<std::int64_t>(randomBetween(rng, 0, 10)) - 5;
        do {
            result.secondRoot = static_cast<std::int64_t>(randomBetween(rng, 0, 10)) - 5;
        } while (result.secondRoot == result.firstRoot);
        result.expression = GeneratedExpr{
            "(x-(" + std::to_string(result.firstRoot) + "))*(x-("
                + std::to_string(result.secondRoot) + "))==0",
            {}};
    }
    else if (choice < 71) {
        result.invariant = Invariant::InverseReconstruction;
        result.expression = generateInvertibleMatrix(rng, randomBetween(rng, 1, 4));
    }
    else if (choice < 74) {
        result.invariant = Invariant::FourierRoundTrip;
        result.expression = generateFourierVector(rng);
    }
    else if (choice < 78) {
        static constexpr std::array<std::string_view, 7> boundaries{
            "1/0", "0^0", "cot[0]", "atan2[0,0]", "atanh[1]",
            "sqrt[-1]", "sqrt[(-3)^2]"};
        result.invariant = Invariant::DomainAndBranchBoundary;
        result.boundaryCase = randomBetween(rng, 0, boundaries.size() - 1);
        result.expression = GeneratedExpr{
            std::string{boundaries[result.boundaryCase]}, {}};
    }
    else if (choice < 83) {
        result.invariant = Invariant::PolynomialLimitMatchesSubstitution;
        result.substitution = static_cast<std::int64_t>(randomBetween(rng, 0, 8)) - 4;
        result.expression = generatePolynomial(
            rng, 0, std::min<std::size_t>(result.targetDepth, 4));
    }
    else if (choice < 87) {
        result.invariant = Invariant::CasesShortCircuit;
        result.expression = generateExactScalar(
            rng, 0, std::min<std::size_t>(result.targetDepth, 5));
    }
    else if (choice < 91) {
        static constexpr std::array<std::string_view, 5> exactValues{
            "Pi", "E", "sqrt[2]", "root[{-2,0,1},2]", "root[{-3,0,1},2]"};
        result.invariant = Invariant::NestedNDoesNotInventPrecision;
        result.precisionLow = randomBetween(rng, 12, 32);
        result.precisionHigh = result.precisionLow + randomBetween(rng, 8, 32);
        result.expression = GeneratedExpr{
            std::string{exactValues[randomBetween(rng, 0, exactValues.size() - 1)]}, {}};
    }
    else if (choice < 94) {
        static constexpr std::array<std::int64_t, 7> radicands{2, 3, 5, 6, 7, 10, 11};
        const std::int64_t radicand = radicands[randomBetween(rng, 0, radicands.size() - 1)];
        result.invariant = Invariant::AlgebraicRootIdentity;
        result.substitution = radicand;
        result.expression = GeneratedExpr{
            "root[{-" + std::to_string(radicand) + ",0,1},2]", {}};
    }
    else if (choice < 97) {
        const std::int64_t a = static_cast<std::int64_t>(randomBetween(rng, 0, 6)) - 3;
        const std::int64_t b = static_cast<std::int64_t>(randomBetween(rng, 1, 6));
        static constexpr std::array<std::string_view, 3> orders{"Lex", "GrLex", "GrevLex"};
        result.invariant = Invariant::GroebnerGeneratorsReduceToZero;
        result.termOrder = std::string{orders[randomBetween(rng, 0, orders.size() - 1)]};
        GeneratedExpr first{"x-y-(" + std::to_string(a) + ")", {}};
        GeneratedExpr second{"y^2-(" + std::to_string(b) + ")", {}};
        result.expression = GeneratedExpr{
            "{" + first.source + "," + second.source + "}",
            {std::move(first), std::move(second)}};
    }
    else if (choice < 99) {
        result.invariant = Invariant::ArrayReshapeRoundTrip;
        result.rows = randomBetween(rng, 1, 4);
        result.columns = randomBetween(rng, 1, 4);
        result.expression = generateReshapeVector(rng, result.rows, result.columns);
    }
    else if (choice < 101) {
        static constexpr std::array<std::string_view, 5> holes{
            "1/x-1/x",
            "0*(1/x)",
            "exp[1/x]/exp[1/x]",
            "exp[1/x]^0",
            "sin[1/x]^2+cos[1/x]^2"};
        result.invariant = Invariant::SimplifyPreservesDomainHole;
        result.boundaryCase = randomBetween(rng, 0, holes.size() - 1);
        result.expression = GeneratedExpr{std::string{holes[result.boundaryCase]}, {}};
    }
    else {
        result.invariant = Invariant::RemovableLimitPreservesDomainHole;
        result.firstRoot = static_cast<std::int64_t>(randomBetween(rng, 0, 6)) - 3;
        result.secondRoot = static_cast<std::int64_t>(randomBetween(rng, 1, 6));
        const std::string a = std::to_string(result.firstRoot);
        const std::string b = std::to_string(result.secondRoot);
        result.expression = GeneratedExpr{
            "((x-(" + a + "))*(x+(" + b + ")))/(x-(" + a + "))", {}};
    }
    return result;
}

[[nodiscard]] Expr evaluate(mmcal::kernel::KernelSession& session, std::string_view source) {
    return session.evaluate(source);
}

[[nodiscard]] const numeric::Rational& fourierOracleTolerance() {
    static const numeric::Rational tolerance = [] {
        numeric::BigInt denominator{1};
        for (std::size_t i = 0; i < 40; ++i)
            denominator *= numeric::BigInt{10};
        return numeric::Rational{numeric::BigInt{1}, std::move(denominator)};
    }();
    return tolerance;
}

[[nodiscard]] bool intervalCertifiesExact(
    const approximation::RealInterval& enclosure,
    const numeric::Rational& expected) {
    if (!enclosure.contains(expected))
        return false;

    const numeric::Rational& tolerance = fourierOracleTolerance();
    const numeric::Rational lower = enclosure.lower().toRational();
    const numeric::Rational upper = enclosure.upper().toRational();
    return lower >= expected - tolerance && upper <= expected + tolerance;
}

[[nodiscard]] bool certifiedValueContainsExact(
    const approximation::CertifiedValue& enclosure,
    const Expr& exact) {
    if (!exact.isNumber())
        return false;

    const numeric::Number& expected = exact.asNumber();
    const numeric::Rational expectedReal = expected.realPart().toRational();
    const numeric::Rational expectedImaginary = expected.imaginaryPart().toRational();

    if (enclosure.isReal())
        return expectedImaginary.isZero() && intervalCertifiesExact(enclosure.asReal(), expectedReal);

    const auto& complex = enclosure.asComplex();
    return intervalCertifiesExact(complex.real(), expectedReal)
        && intervalCertifiesExact(complex.imaginary(), expectedImaginary);
}

[[nodiscard]] Expr substituteSymbolX(const Expr& expression, const Expr& replacement) {
    if (expression.isSymbol() && expression.asSymbol().view() == "x")
        return replacement;
    if (!expression.isCall())
        return expression;

    std::vector<Expr> arguments;
    arguments.reserve(expression.asCall().arguments.size());
    for (const Expr& argument : expression.asCall().arguments)
        arguments.push_back(substituteSymbolX(argument, replacement));
    return Expr::rebuildCall(expression.asCall(), std::move(arguments));
}

[[nodiscard]] mathematics::Predicate substitutePredicateX(
    const mathematics::Predicate& predicate,
    const Expr& replacement) {
    if (const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate)) {
        return mathematics::relation(
            relation->relation,
            substituteSymbolX(relation->lhs, replacement),
            substituteSymbolX(relation->rhs, replacement));
    }
    const auto& domain = std::get<mathematics::DomainPredicate>(predicate);
    return mathematics::elementOf(
        substituteSymbolX(domain.expression, replacement), domain.domain);
}

[[nodiscard]] std::optional<bool> domainIsDefinedAt(
    const mathematics::AssumptionSet& domain,
    std::int64_t point,
    const mmcal::kernel::KernelSession& session) {
    const Expr replacement{numeric::Number{numeric::BigInt{point}}};
    const mathematics::AssumptionSet none;
    const mathematics::KnowledgeContext knowledge{
        session.builtinRegistry(), session.mathRegistry(), none};

    for (const mathematics::Predicate& predicate : domain.predicates()) {
        const mathematics::TruthValue truth = knowledge.prove(
            substitutePredicateX(predicate, replacement));
        if (truth == mathematics::TruthValue::False)
            return false;
        if (truth != mathematics::TruthValue::True)
            return std::nullopt;
    }
    return true;
}

[[nodiscard]] std::optional<std::string> compareDomainOnIntegerWindow(
    const mathematics::AssumptionSet& expected,
    const mathematics::AssumptionSet& actual,
    std::int64_t center,
    const mmcal::kernel::KernelSession& session) {
    for (std::int64_t point = center - 6; point <= center + 6; ++point) {
        const auto expectedDefined = domainIsDefinedAt(expected, point, session);
        const auto actualDefined = domainIsDefinedAt(actual, point, session);
        if (!expectedDefined || !actualDefined)
            return "domain predicate remained undecidable at x=" + std::to_string(point);
        if (*expectedDefined != *actualDefined)
            return "definedness differs at x=" + std::to_string(point)
                + " expected=" + (*expectedDefined ? "defined" : "undefined")
                + " actual=" + (*actualDefined ? "defined" : "undefined");
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Failure> checkCase(
    mmcal::kernel::KernelSession& session,
    const GeneratedCase& testCase,
    std::string_view overrideSource = {}) {
    const std::string source = overrideSource.empty()
        ? testCase.expression.source
        : std::string{overrideSource};

    try {
        if (testCase.invariant == Invariant::ExactRoundTrip) {
            session.resetForIndependentEvaluation();
            const Expr original = evaluate(session, source);
            const std::string formatted = mmcal::formatting::formatExpr(original);
            session.resetForIndependentEvaluation();
            const Expr reparsed = evaluate(session, formatted);
            const std::string reformatted = mmcal::formatting::formatExpr(reparsed);
            if (formatted != reformatted)
                return Failure{"format(parse(format(expr))) is not stable", formatted, reformatted};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::FullSimplifyPreservesValue) {
            session.resetForIndependentEvaluation();
            const Expr original = evaluate(session, source);
            session.resetForIndependentEvaluation();
            const Expr simplified = evaluate(session, "fullSimplify[" + source + "]");
            if (!(original == simplified))
                return Failure{"fullSimplify changed an exact numeric value",
                    mmcal::formatting::formatExpr(original),
                    mmcal::formatting::formatExpr(simplified)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::ExpandPreservesPolynomial
            || testCase.invariant == Invariant::FactorPreservesPolynomial) {
            session.resetForIndependentEvaluation();
            const std::string transformedSource = testCase.invariant == Invariant::ExpandPreservesPolynomial
                ? "expand[" + source + "]"
                : "factor[expand[" + source + "]]";
            const Expr transformed = evaluate(session, transformedSource);
            const std::string transformedText = mmcal::formatting::formatExpr(transformed);

            // v1.5.3までは一点代入だけだったため，異なる多項式が偶然同じ値になる
            // blind spotがあった。まず形式差を展開し，さらに独立な三点で照合する。
            session.resetForIndependentEvaluation();
            const Expr residual = evaluate(
                session,
                "expand[(" + source + ")-(" + transformedText + ")]");
            session.resetForIndependentEvaluation();
            const Expr zero = evaluate(session, "0");
            if (!(residual == zero))
                return Failure{
                    testCase.invariant == Invariant::ExpandPreservesPolynomial
                        ? "expand left a nonzero symbolic polynomial residual"
                        : "factor[expand[...]] left a nonzero symbolic polynomial residual",
                    "0",
                    mmcal::formatting::formatExpr(residual)};

            const std::array<std::int64_t, 3> substitutions{
                testCase.substitution,
                testCase.substitution + 11,
                testCase.substitution - 13};
            for (const std::int64_t substitution : substitutions) {
                const std::string originalAtPoint = replaceSymbolX(source, substitution);
                const std::string transformedAtPoint = replaceSymbolX(transformedText, substitution);
                session.resetForIndependentEvaluation();
                const Expr lhs = evaluate(session, originalAtPoint);
                session.resetForIndependentEvaluation();
                const Expr rhs = evaluate(session, transformedAtPoint);
                if (!(lhs == rhs))
                    return Failure{
                        (testCase.invariant == Invariant::ExpandPreservesPolynomial
                            ? "expand changed polynomial value at x="
                            : "factor[expand[...]] changed polynomial value at x=")
                            + std::to_string(substitution),
                        mmcal::formatting::formatExpr(lhs),
                        mmcal::formatting::formatExpr(rhs)};
            }
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::DoubleTranspose) {
            session.resetForIndependentEvaluation();
            const Expr original = evaluate(session, source);
            session.resetForIndependentEvaluation();
            const Expr transformed = evaluate(session, "transpose[transpose[" + source + "]]" );
            if (!(original == transformed))
                return Failure{"transpose[transpose[A]] != A",
                    mmcal::formatting::formatExpr(original),
                    mmcal::formatting::formatExpr(transformed)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::DeterminantTranspose) {
            session.resetForIndependentEvaluation();
            const Expr lhs = evaluate(session, "det[" + source + "]");
            session.resetForIndependentEvaluation();
            const Expr rhs = evaluate(session, "det[transpose[" + source + "]]" );
            if (!(lhs == rhs))
                return Failure{"det[A] != det[transpose[A]]",
                    mmcal::formatting::formatExpr(lhs),
                    mmcal::formatting::formatExpr(rhs)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::DerivativeOfIntegral) {
            session.resetForIndependentEvaluation();
            const Expr residual = evaluate(
                session,
                "fullSimplify[D[integrate[(" + source + "),x],x]-(" + source + ")]" );
            session.resetForIndependentEvaluation();
            const Expr zero = evaluate(session, "0");
            if (!(residual == zero))
                return Failure{"D[integrate[p,x],x] left a nonzero residual",
                    "0", mmcal::formatting::formatExpr(residual)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::SolveKnownRoots) {
            session.resetForIndependentEvaluation();
            const Expr solved = evaluate(session, "solve[" + source + ",x]");
            if (!solved.isSolutionSet())
                return Failure{"solve did not return a SolutionSet",
                    "two finite roots", mmcal::formatting::formatExpr(solved)};

            const solver::SolutionSet& solutions = solved.asSolutionSet();
            if (solutions.kind() != solver::SolutionSetKind::Finite)
                return Failure{"solve did not return a finite solution set",
                    "Finite", std::to_string(static_cast<int>(solutions.kind()))};
            std::vector<std::string> actualRoots;
            for (const solver::SolutionBranch& branch : solutions.branches()) {
                if (branch.bindings.size() != 1)
                    return Failure{"solve returned a branch with unexpected arity",
                        "one binding", std::to_string(branch.bindings.size())};
                actualRoots.push_back(
                    mmcal::formatting::formatExpr(branch.bindings.front().value));
            }
            std::vector<std::string> expectedRoots{
                std::to_string(testCase.firstRoot),
                std::to_string(testCase.secondRoot)};
            std::sort(actualRoots.begin(), actualRoots.end());
            std::sort(expectedRoots.begin(), expectedRoots.end());
            if (actualRoots != expectedRoots) {
                const auto join = [](const std::vector<std::string>& values) {
                    std::string text;
                    for (const std::string& value : values) {
                        if (!text.empty())
                            text += ", ";
                        text += value;
                    }
                    return text;
                };
                return Failure{"solve roots differ from the constructed roots",
                    join(expectedRoots), join(actualRoots)};
            }
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::InverseReconstruction) {
            session.resetForIndependentEvaluation();
            const Expr matrix = evaluate(session, source);
            if (!matrix.isArray() || !matrix.asArray().isMatrix()
                || matrix.asArray().shape[0] != matrix.asArray().shape[1])
                return Failure{"invertible-matrix generator produced a non-square matrix",
                    "square matrix", mmcal::formatting::formatExpr(matrix)};
            const std::size_t size = matrix.asArray().shape[0];
            session.resetForIndependentEvaluation();
            const Expr reconstructed = evaluate(
                session, "dot[" + source + ",inverse[" + source + "]]" );
            session.resetForIndependentEvaluation();
            const Expr identity = evaluate(session, identityMatrixSource(size));
            if (!(reconstructed == identity))
                return Failure{"dot[A,inverse[A]] != I",
                    mmcal::formatting::formatExpr(identity),
                    mmcal::formatting::formatExpr(reconstructed)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::FourierRoundTrip) {
            session.resetForIndependentEvaluation();
            const Expr original = evaluate(session, source);
            session.resetForIndependentEvaluation();
            // Bluestein/Cooley--Tukeyの中間式には完全に相殺するexact radicalが
            // 残り得るため，数学的不変量の比較前に共通のcanonicalizerを通す。
            const Expr roundTrip = evaluate(
                session, "fullSimplify[ifft[fft[" + source + "]]]" );
            if (roundTrip == original)
                return std::nullopt;

            // 非2冪長ではexact root-of-unityの相殺形がcanonicalizerに残ることがある。
            // actual-expectedを先に構成すると，数学的には0でもradical cancellationが
            // bounded refinementを使い切るため，exact round-trip Exprそのものを
            // CertifiedEvaluatorへ渡す。各enclosureが元のexact Gaussian整数点を含むことを
            // 確認し，Formatter/reparseやsubtractive cancellationをoracleへ持ち込まない。
            if (roundTrip.isArray() && original.isArray()
                && roundTrip.asArray().shape == original.asArray().shape) {
                const mathematics::AngleSemantics angles{session.defaultAngleUnit()};
                const approximation::CertifiedEvaluator certified{
                    session.builtinRegistry(), session.mathRegistry(), angles};
                const approximation::ApproximationContext context{50};

                bool certifiedEqual = true;
                bool certifiedAvailable = true;
                for (std::size_t index = 0; index < original.asArray().size(); ++index) {
                    const auto enclosed = certified.enclose(
                        roundTrip.asArray().element(index), context.workingBinaryBits());
                    if (!enclosed) {
                        certifiedAvailable = false;
                        break;
                    }
                    if (!certifiedValueContainsExact(
                            *enclosed, original.asArray().element(index))) {
                        certifiedEqual = false;
                        break;
                    }
                }
                if (certifiedAvailable && certifiedEqual)
                    return std::nullopt;
                if (!certifiedAvailable)
                    return Failure{
                        "FFT exact oracle exceeded the certified evaluation work/depth budget",
                        mmcal::formatting::formatExpr(original),
                        mmcal::formatting::formatExpr(roundTrip), {}, true};
            }

            return Failure{"ifft[fft[v]] != v",
                mmcal::formatting::formatExpr(original),
                mmcal::formatting::formatExpr(roundTrip)};
        }

        if (testCase.invariant == Invariant::DomainAndBranchBoundary) {
            session.resetForIndependentEvaluation();
            if (testCase.boundaryCase == 0 || testCase.boundaryCase == 1) {
                const Expr actual = evaluate(session, source);
                session.resetForIndependentEvaluation();
                const Expr expected = evaluate(
                    session, testCase.boundaryCase == 0 ? "ComplexInfinity" : "Indeterminate");
                if (!(actual == expected))
                    return Failure{"special undefined-value boundary changed",
                        mmcal::formatting::formatExpr(expected),
                        mmcal::formatting::formatExpr(actual)};
                return std::nullopt;
            }

            if (testCase.boundaryCase < 5) {
                try {
                    const Expr unexpected = evaluate(session, source);
                    return Failure{"domain boundary did not raise DomainError",
                        "DomainError", mmcal::formatting::formatExpr(unexpected)};
                }
                catch (const error::CalcError& exception) {
                    if (exception.type() != error::CalcErrorType::Domain)
                        return Failure{"domain boundary raised the wrong diagnostic class",
                            "DomainError", std::string{error::calcErrorTypeName(exception.type())}};
                    return std::nullopt;
                }
            }

            const Expr actual = evaluate(session, source);
            session.resetForIndependentEvaluation();
            const Expr expected = evaluate(
                session, testCase.boundaryCase == 5 ? "I" : "3");
            if (!(actual == expected))
                return Failure{"principal sqrt branch boundary changed",
                    mmcal::formatting::formatExpr(expected),
                    mmcal::formatting::formatExpr(actual)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::PolynomialLimitMatchesSubstitution) {
            session.resetForIndependentEvaluation();
            const Expr limited = evaluate(
                session,
                "limit[(" + source + "),x," + std::to_string(testCase.substitution) + "]");
            session.resetForIndependentEvaluation();
            const Expr substituted = evaluate(
                session, replaceSymbolX(source, testCase.substitution));
            if (!(limited == substituted))
                return Failure{"polynomial limit differs from direct substitution",
                    mmcal::formatting::formatExpr(substituted),
                    mmcal::formatting::formatExpr(limited)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::CasesShortCircuit) {
            session.resetForIndependentEvaluation();
            const Expr expected = evaluate(session, source);
            session.resetForIndependentEvaluation();
            const Expr actual = evaluate(
                session, "cases[1/0 if False;" + source + " if True;1/0]");
            if (!(actual == expected))
                return Failure{"cases evaluated an inactive branch or changed the active value",
                    mmcal::formatting::formatExpr(expected),
                    mmcal::formatting::formatExpr(actual)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::NestedNDoesNotInventPrecision) {
            const std::string low = std::to_string(testCase.precisionLow);
            const std::string high = std::to_string(testCase.precisionHigh);
            session.resetForIndependentEvaluation();
            const Expr first = evaluate(session, "N[" + source + "," + low + "]");
            session.resetForIndependentEvaluation();
            const Expr nested = evaluate(
                session, "N[N[" + source + "," + low + "]," + high + "]");
            if (!(nested == first))
                return Failure{"nested N invented information beyond the first approximation",
                    mmcal::formatting::formatExpr(first),
                    mmcal::formatting::formatExpr(nested)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::AlgebraicRootIdentity) {
            session.resetForIndependentEvaluation();
            const Expr residual = evaluate(
                session,
                "fullSimplify[(" + source + ")^2-(" + std::to_string(testCase.substitution) + ")]");
            session.resetForIndependentEvaluation();
            const Expr zero = evaluate(session, "0");
            if (!(residual == zero))
                return Failure{"algebraic square root does not satisfy its defining polynomial",
                    "0", mmcal::formatting::formatExpr(residual)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::GroebnerGeneratorsReduceToZero) {
            for (const GeneratedExpr& generator : testCase.expression.children) {
                session.resetForIndependentEvaluation();
                const Expr remainder = evaluate(
                    session,
                    "at[polynomialReduce[" + generator.source + ",groebnerBasis["
                        + source + ",{x,y}," + testCase.termOrder + "],{x,y},"
                        + testCase.termOrder + "],1]");
                if (!remainder.isNumber() || !remainder.asNumber().isZero())
                    return Failure{"Groebner basis did not reduce an input generator to zero",
                        "0", mmcal::formatting::formatExpr(remainder)};
            }
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::ArrayReshapeRoundTrip) {
            const std::size_t length = testCase.rows * testCase.columns;
            session.resetForIndependentEvaluation();
            const Expr original = evaluate(session, source);
            session.resetForIndependentEvaluation();
            const Expr roundTrip = evaluate(
                session,
                "reshape[reshape[" + source + ",{" + std::to_string(testCase.rows)
                    + "," + std::to_string(testCase.columns) + "}],{"
                    + std::to_string(length) + "}]");
            if (!(roundTrip == original))
                return Failure{"reshape matrix/vector round-trip changed array values or order",
                    mmcal::formatting::formatExpr(original),
                    mmcal::formatting::formatExpr(roundTrip)};
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::SimplifyPreservesDomainHole) {
            session.resetForIndependentEvaluation();
            const Expr original = evaluate(session, source);
            const auto originalDomain = mathematics::expressionDomainConditions(
                original, session.builtinRegistry(), session.mathRegistry());
            if (!originalDomain)
                return Failure{"definedness collector could not describe the generated source",
                    "known domain conditions", mmcal::formatting::formatExpr(original)};

            for (const std::string_view functionName : {"simplify", "fullSimplify"}) {
                session.resetForIndependentEvaluation();
                const Expr simplified = evaluate(
                    session, std::string{functionName} + "[" + source + "]");
                const auto simplifiedDomain = mathematics::expressionDomainConditions(
                    simplified, session.builtinRegistry(), session.mathRegistry());
                if (!simplifiedDomain)
                    return Failure{std::string{functionName}
                            + " produced an expression with unknown definedness",
                        "known domain conditions", mmcal::formatting::formatExpr(simplified)};
                if (const auto mismatch = compareDomainOnIntegerWindow(
                        *originalDomain, *simplifiedDomain, 0, session))
                    return Failure{std::string{functionName} + " changed symbolic domain conditions",
                        "same exact definedness on generated rational domain", *mismatch};

                session.resetForIndependentEvaluation();
                const Expr underAssumption = evaluate(
                    session, std::string{functionName} + "[" + source + ",x!=0]");
                session.resetForIndependentEvaluation();
                const Expr expected = evaluate(
                    session, testCase.boundaryCase <= 1 ? "0" : "1");
                if (!(underAssumption == expected))
                    return Failure{std::string{functionName}
                            + " did not exploit an explicit nonzero assumption",
                        mmcal::formatting::formatExpr(expected),
                        mmcal::formatting::formatExpr(underAssumption)};
            }
            return std::nullopt;
        }

        if (testCase.invariant == Invariant::RemovableLimitPreservesDomainHole) {
            const std::int64_t point = testCase.firstRoot;
            const std::int64_t expectedLimit = point + testCase.secondRoot;

            session.resetForIndependentEvaluation();
            const Expr limited = evaluate(
                session, "limit[" + source + ",x," + std::to_string(point) + "]");
            session.resetForIndependentEvaluation();
            const Expr expected = evaluate(session, std::to_string(expectedLimit));
            if (!(limited == expected))
                return Failure{"removable singularity limit is incorrect",
                    mmcal::formatting::formatExpr(expected),
                    mmcal::formatting::formatExpr(limited)};

            session.resetForIndependentEvaluation();
            const Expr original = evaluate(session, source);
            const auto originalDomain = mathematics::expressionDomainConditions(
                original, session.builtinRegistry(), session.mathRegistry());
            if (!originalDomain)
                return Failure{"definedness collector could not describe removable singularity",
                    "known domain conditions", mmcal::formatting::formatExpr(original)};

            session.resetForIndependentEvaluation();
            const Expr simplified = evaluate(session, "fullSimplify[" + source + "]");
            const auto simplifiedDomain = mathematics::expressionDomainConditions(
                simplified, session.builtinRegistry(), session.mathRegistry());
            if (!simplifiedDomain)
                return Failure{"fullSimplify produced unknown definedness at removable singularity",
                    "known domain conditions", mmcal::formatting::formatExpr(simplified)};
            if (const auto mismatch = compareDomainOnIntegerWindow(
                    *originalDomain, *simplifiedDomain, point, session))
                return Failure{"fullSimplify erased or changed a removable singularity",
                    "same exact definedness on generated rational domain", *mismatch};
            return std::nullopt;
        }

        return Failure{"unknown invariant", "implemented invariant", invariantName(testCase.invariant)};
    }
    catch (const std::exception& exception) {
        return Failure{"unexpected exception", "successful invariant evaluation", exception.what()};
    }
    catch (...) {
        return Failure{"unexpected non-standard exception", "successful invariant evaluation", "unknown exception"};
    }
}

[[nodiscard]] std::string shrinkFailure(
    mmcal::kernel::KernelSession& session,
    const GeneratedCase& testCase) {
    // 行列ではshapeを壊す縮約を避ける。scalar/polynomialだけ部分式へ縮める。
    if (testCase.invariant == Invariant::DoubleTranspose
        || testCase.invariant == Invariant::DeterminantTranspose
        || testCase.invariant == Invariant::SolveKnownRoots
        || testCase.invariant == Invariant::InverseReconstruction
        || testCase.invariant == Invariant::FourierRoundTrip
        || testCase.invariant == Invariant::DomainAndBranchBoundary
        || testCase.invariant == Invariant::NestedNDoesNotInventPrecision
        || testCase.invariant == Invariant::AlgebraicRootIdentity
        || testCase.invariant == Invariant::GroebnerGeneratorsReduceToZero
        || testCase.invariant == Invariant::ArrayReshapeRoundTrip
        || testCase.invariant == Invariant::SimplifyPreservesDomainHole
        || testCase.invariant == Invariant::RemovableLimitPreservesDomainHole)
        return testCase.expression.source;

    std::string best = testCase.expression.source;
    std::vector<const GeneratedExpr*> frontier;
    frontier.reserve(testCase.expression.children.size());
    for (const auto& child : testCase.expression.children)
        frontier.push_back(&child);

    static constexpr std::array<std::string_view, 3> simple{"0", "1", "-1"};
    for (const std::string_view candidate : simple) {
        if (checkCase(session, testCase, candidate)) {
            best = std::string{candidate};
            return best;
        }
    }

    while (!frontier.empty()) {
        const GeneratedExpr* candidate = frontier.back();
        frontier.pop_back();
        if (checkCase(session, testCase, candidate->source)) {
            best = candidate->source;
            frontier.clear();
            for (const auto& child : candidate->children)
                frontier.push_back(&child);
        }
    }
    return best;
}

[[nodiscard]] std::string budgetSummary(const evaluation::EvaluationUsage& usage) {
    std::ostringstream output;
    output
        << "input-bytes=" << usage.inputBytes
        << " work=" << usage.evaluationSteps
        << " depth=" << usage.maximumDepth
        << " nodes=" << usage.generatedNodes
        << " simplify=" << usage.simplificationCandidates
        << " solve=" << usage.solverBranches
        << " integrate=" << usage.integrationCandidates
        << " refine=" << usage.certifiedRefinements
        << " arrays=" << usage.denseArrayElements
        << " matrix-temp=" << usage.temporaryMatrixElements
        << " bigint-bits=" << usage.maximumBigIntegerBits
        << " precision-digits=" << usage.maximumRequestedPrecisionDigits
        << " algebraic-degree=" << usage.maximumAlgebraicDegree
        << " algebraic-refine=" << usage.algebraicRefinements;
    return output.str();
}

void printInconclusive(
    const RandomExpressionFuzzerOptions& options,
    std::uint64_t index,
    const GeneratedCase& testCase,
    const Failure& result) {
    std::cerr
        << "\n[INCONCLUSIVE] Random expression invariant\n"
        << "Seed      : " << options.seed << '\n'
        << "Case      : " << index << '\n'
        << "Case seed : " << caseSeed(options.seed, index) << '\n'
        << "Depth     : " << testCase.targetDepth << '\n'
        << "Invariant : " << invariantName(testCase.invariant) << '\n'
        << "Expression: " << testCase.expression.source << '\n'
        << "Reason    : " << result.reason << '\n'
        << "Reproduce : mmCal.Benchmarks --random-expressions --seed "
        << options.seed << " --case " << index << '\n';
}

void printFailure(
    const RandomExpressionFuzzerOptions& options,
    std::uint64_t index,
    const GeneratedCase& testCase,
    const Failure& failure,
    std::string_view reduced) {
    std::cerr
        << "\n[FAIL] Random expression invariant\n"
        << "Seed      : " << options.seed << '\n'
        << "Case      : " << index << '\n'
        << "Case seed : " << caseSeed(options.seed, index) << '\n'
        << "Depth     : " << testCase.targetDepth << '\n'
        << "Invariant : " << invariantName(testCase.invariant) << '\n'
        << "Expression: " << testCase.expression.source << '\n';
    if (!reduced.empty() && reduced != testCase.expression.source)
        std::cerr << "Reduced   : " << reduced << '\n';
    std::cerr
        << "Reason    : " << failure.reason << '\n'
        << "Expected  : " << failure.expected << '\n'
        << "Actual    : " << failure.actual << '\n';
    if (!failure.budget.empty())
        std::cerr << "Budget    : " << failure.budget << '\n';
    std::cerr
        << "Reproduce : mmCal.Benchmarks --random-expressions --seed "
        << options.seed << " --case " << index << '\n';
}

} // namespace

bool runRandomExpressionFuzzer(const RandomExpressionFuzzerOptions& options) {
    const auto started = Clock::now();
    const std::uint64_t firstCase = options.singleCase.value_or(1);
    const std::uint64_t lastCase = options.singleCase
        ? firstCase
        : options.loop ? std::numeric_limits<std::uint64_t>::max() : options.cases;
    const std::size_t workerCount = options.singleCase ? 1 : std::max<std::size_t>(1, options.threads);

    std::cout
        << "Random expression fuzzing\n"
        << "Seed      : " << options.seed << '\n'
        << "Mode      : " << (options.singleCase
            ? "single case"
            : options.noStopLoop
                ? "infinite (FAIL is reported and fuzzing continues)"
                : options.loop
                    ? "infinite (FAIL stops immediately)"
                    : "finite") << '\n'
        << "Threads   : " << workerCount << '\n'
        << "Max depth : " << options.maxDepth << '\n';
    if (!options.loop && !options.singleCase)
        std::cout << "Cases     : " << options.cases << '\n';
    std::cout << std::flush;

    // single caseはthreadを作らず，failure reproductionを最短経路にする。
    if (options.singleCase) {
        mmcal::kernel::KernelSession session;
        const GeneratedCase testCase = generateCase(options.seed, firstCase, options.maxDepth);
        if (auto failure = checkCase(session, testCase)) {
            failure->budget = budgetSummary(session.lastEvaluationUsage());
            if (failure->inconclusive) {
                printInconclusive(options, firstCase, testCase, *failure);
                return true;
            }
            const std::string reduced = shrinkFailure(session, testCase);
            printFailure(options, firstCase, testCase, *failure, reduced);
            return false;
        }

        const double seconds = std::chrono::duration<double>(Clock::now() - started).count();
        std::cout << "Random expression case: PASS  " << seconds << " s\n";
        return true;
    }

    std::atomic<std::uint64_t> nextCase{firstCase};
    std::atomic<std::uint64_t> completed{0};
    std::atomic<std::uint64_t> failures{0};
    std::atomic<std::uint64_t> inconclusive{0};
    std::atomic<std::size_t> observedMaxDepth{0};
    std::atomic<bool> stop{false};
    std::mutex outputMutex;

    const auto updateMaxDepth = [&](std::size_t depth) {
        std::size_t observed = observedMaxDepth.load(std::memory_order_relaxed);
        while (observed < depth
            && !observedMaxDepth.compare_exchange_weak(
                observed, depth, std::memory_order_relaxed)) {}
    };

    const auto worker = [&]() {
        // KernelSessionはmutable stateを持つためworker間で共有しない。
        mmcal::kernel::KernelSession session;

        while (!stop.load(std::memory_order_relaxed)) {
            const std::uint64_t index = nextCase.fetch_add(1, std::memory_order_relaxed);
            if (!options.loop && index > lastCase)
                break;
            if (index == 0) // uint64_t wraparound。実運用では到達しないが無限loopを安全に終える。
                break;

            const GeneratedCase testCase = generateCase(options.seed, index, options.maxDepth);
            updateMaxDepth(testCase.targetDepth);

            if (auto failure = checkCase(session, testCase)) {
                failure->budget = budgetSummary(session.lastEvaluationUsage());
                if (failure->inconclusive) {
                    inconclusive.fetch_add(1, std::memory_order_relaxed);
                }
                else {
                    failures.fetch_add(1, std::memory_order_relaxed);

                    if (!options.noStopLoop) {
                        bool expected = false;
                        if (stop.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
                            const std::string reduced = shrinkFailure(session, testCase);
                            std::scoped_lock lock{outputMutex};
                            printFailure(options, index, testCase, *failure, reduced);
                        }
                        break;
                    }

                    const std::string reduced = shrinkFailure(session, testCase);
                    std::scoped_lock lock{outputMutex};
                    printFailure(options, index, testCase, *failure, reduced);
                }
            }

            const std::uint64_t checked = completed.fetch_add(1, std::memory_order_relaxed) + 1;
            if (options.reportEvery != 0 && checked % options.reportEvery == 0) {
                const double seconds = std::chrono::duration<double>(Clock::now() - started).count();
                std::scoped_lock lock{outputMutex};
                std::cout << '[' << checked << "] "
                          << (failures.load(std::memory_order_relaxed) == 0 ? "PASS" : "RUN")
                          << "  " << seconds << " s  max-depth="
                          << observedMaxDepth.load(std::memory_order_relaxed);
                if (options.noStopLoop)
                    std::cout << "  failures=" << failures.load(std::memory_order_relaxed);
                const auto inconclusiveCount = inconclusive.load(std::memory_order_relaxed);
                if (inconclusiveCount != 0)
                    std::cout << "  inconclusive=" << inconclusiveCount;
                std::cout << '\n';
            }
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(workerCount);
    for (std::size_t i = 0; i < workerCount; ++i)
        workers.emplace_back(worker);
    for (auto& thread : workers)
        thread.join();

    const double seconds = std::chrono::duration<double>(Clock::now() - started).count();
    const std::uint64_t failureCount = failures.load(std::memory_order_relaxed);
    if (!options.loop) {
        std::cout << "Random expression fuzzing: "
                  << (failureCount == 0 ? "PASS  " : "FAIL  ")
                  << completed.load(std::memory_order_relaxed) << " cases  "
                  << seconds << " s";
        if (failureCount != 0)
            std::cout << "  failures=" << failureCount;
        const auto inconclusiveCount = inconclusive.load(std::memory_order_relaxed);
        if (inconclusiveCount != 0)
            std::cout << "  inconclusive=" << inconclusiveCount;
        std::cout << '\n';
    }
    return failureCount == 0;
}

} // namespace mmcal::benchmarks
