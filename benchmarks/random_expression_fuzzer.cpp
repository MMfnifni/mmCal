#include "random_expression_fuzzer.hpp"

#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"

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
    DeterminantTranspose
};

struct GeneratedCase final {
    Invariant invariant = Invariant::ExactRoundTrip;
    GeneratedExpr expression;
    std::size_t targetDepth = 1;
    std::int64_t substitution = 0;
};

struct Failure final {
    std::string reason;
    std::string expected;
    std::string actual;
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

    const unsigned choice = static_cast<unsigned>(randomBetween(rng, 0, 99));
    if (choice < 35) {
        result.invariant = Invariant::ExactRoundTrip;
        result.expression = generateExactScalar(rng, 0, result.targetDepth);
    }
    else if (choice < 60) {
        result.invariant = Invariant::FullSimplifyPreservesValue;
        result.expression = generateExactScalar(rng, 0, result.targetDepth);
    }
    else if (choice < 75) {
        result.invariant = Invariant::ExpandPreservesPolynomial;
        result.expression = generatePolynomial(rng, 0, std::min<std::size_t>(result.targetDepth, 6));
    }
    else if (choice < 87) {
        result.invariant = Invariant::FactorPreservesPolynomial;
        result.expression = generatePolynomial(rng, 0, std::min<std::size_t>(result.targetDepth, 5));
    }
    else if (choice < 95) {
        result.invariant = Invariant::DoubleTranspose;
        result.expression = generateMatrix(rng);
    }
    else {
        result.invariant = Invariant::DeterminantTranspose;
        // detには正方行列だけを生成する。
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
    return result;
}

[[nodiscard]] Expr evaluate(mmcal::kernel::KernelSession& session, std::string_view source) {
    return session.evaluate(source);
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

            const std::string originalAtPoint = replaceSymbolX(source, testCase.substitution);
            const std::string transformedAtPoint = replaceSymbolX(transformedText, testCase.substitution);
            session.resetForIndependentEvaluation();
            const Expr lhs = evaluate(session, originalAtPoint);
            session.resetForIndependentEvaluation();
            const Expr rhs = evaluate(session, transformedAtPoint);
            if (!(lhs == rhs))
                return Failure{
                    testCase.invariant == Invariant::ExpandPreservesPolynomial
                        ? "expand changed polynomial value"
                        : "factor[expand[...]] changed polynomial value",
                    mmcal::formatting::formatExpr(lhs),
                    mmcal::formatting::formatExpr(rhs)};
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
        || testCase.invariant == Invariant::DeterminantTranspose)
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
        << "Actual    : " << failure.actual << '\n'
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
        if (const auto failure = checkCase(session, testCase)) {
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

            if (const auto failure = checkCase(session, testCase)) {
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
        std::cout << '\n';
    }
    return failureCount == 0;
}

} // namespace mmcal::benchmarks
