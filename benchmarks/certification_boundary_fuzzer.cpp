#include "certification_boundary_fuzzer.hpp"

#include "error/error_message.hpp"
#include "evaluation/diagnostic.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <random>
#include <sstream>
#include <span>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace mmcal::benchmarks {
namespace {

using Clock = std::chrono::steady_clock;

constexpr std::size_t kPrecisionRefinementCeiling = 4096;

enum class BoundaryOutcome {
    Value,
    DomainError,
    PrecisionInsufficient,
    BackendUnsupported,
    Unevaluated,
    ResourceLimit,
    Timeout,
    OtherError,
    ConflictingDiagnostics
};

struct ProbeTemplate final {
    std::string_view family;
    std::string_view name;
    std::string_view source;
    BoundaryOutcome expected = BoundaryOutcome::Value;
};

struct GeneratedCase final {
    std::string family;
    std::string name;
    std::string source;
    BoundaryOutcome expected = BoundaryOutcome::Value;
    std::size_t requestedDigits = 20;
    std::size_t inputDigits = 5;
    std::size_t epsilonPower = 40;
};

struct CaseObservation final {
    BoundaryOutcome outcome = BoundaryOutcome::OtherError;
    std::string formattedResult;
    std::string detail;
    evaluation::EvaluationUsage usage{};
    double elapsedMilliseconds = 0.0;
};

struct Failure final {
    std::uint64_t index = 0;
    GeneratedCase testCase;
    CaseObservation observation;
    std::string reason;
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

[[nodiscard]] std::size_t choose(
    std::mt19937_64& rng,
    std::span<const std::size_t> values) {
    std::uniform_int_distribution<std::size_t> distribution(0, values.size() - 1);
    return values[distribution(rng)];
}

[[nodiscard]] std::string replaceAll(
    std::string text,
    std::string_view token,
    std::string_view replacement) {
    std::size_t offset = 0;
    while ((offset = text.find(token, offset)) != std::string::npos) {
        text.replace(offset, token.size(), replacement);
        offset += replacement.size();
    }
    return text;
}

[[nodiscard]] std::string instantiate(
    std::string_view source,
    std::size_t requestedDigits,
    std::size_t inputDigits,
    std::size_t epsilonPower) {
    std::string result{source};
    result = replaceAll(std::move(result), "{p}", std::to_string(requestedDigits));
    result = replaceAll(std::move(result), "{q}", std::to_string(inputDigits));
    result = replaceAll(std::move(result), "{k}", std::to_string(epsilonPower));
    return result;
}

[[nodiscard]] constexpr std::array<ProbeTemplate, 46> probeTemplates() {
    return {{
        // Principal branch cuts: exact side/cut points remain evaluable, finite side ambiguity does not.
        {"elementary-branch", "log-exact-upper-side", "N[log[-1+I/10^{k}],{p}]", BoundaryOutcome::Value},
        {"elementary-branch", "log-cut-ambiguity", "N[log[-1+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"elementary-branch", "sqrt-cut-ambiguity", "N[sqrt[-1+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"elementary-branch", "arg-origin-ambiguity", "N[arg[N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"elementary-branch", "atan2-cut-ambiguity", "N[atan2[N[0,{q}],-1],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"elementary-branch", "atan2-origin-ambiguity", "N[atan2[0,N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"elementary-branch", "power-cut-ambiguity", "N[(-1+I*N[0,{q}])^(1/3),{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"elementary-branch", "arbitrary-log-value-cut", "N[log[2,-1+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"elementary-branch", "arbitrary-log-base-cut", "N[log[-2+I*N[0,{q}],8],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"elementary-branch", "arbitrary-log-base-one", "N[log[N[1,{q}],8],{p}]", BoundaryOutcome::PrecisionInsufficient},

        {"inverse-branch", "asin-real-cut", "N[asin[2+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"inverse-branch", "acos-real-cut", "N[acos[2+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"inverse-branch", "atan-imaginary-cut", "N[atan[2I+N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"inverse-branch", "asinh-imaginary-cut", "N[asinh[2I+N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"inverse-branch", "acosh-real-cut", "N[acosh[1/2+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"inverse-branch", "atanh-real-cut", "N[atanh[2+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"inverse-branch", "asin-safe-side", "N[asin[2+I/10^{k}],{p}]", BoundaryOutcome::Value},
        {"inverse-branch", "acosh-safe-side", "N[acosh[1/2+I/10^{k}],{p}]", BoundaryOutcome::Value},

        // Exact singularities are domain errors; finite information merely containing them is precision-insufficient.
        {"special-pole", "gamma-exact-pole", "gamma[0]", BoundaryOutcome::DomainError},
        {"special-pole", "gamma-finite-pole", "N[gamma[N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"special-pole", "digamma-exact-pole", "digamma[0]", BoundaryOutcome::DomainError},
        {"special-pole", "digamma-finite-pole", "N[digamma[N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"special-pole", "trigamma-exact-pole", "trigamma[-1]", BoundaryOutcome::DomainError},
        {"special-pole", "trigamma-finite-pole", "N[trigamma[N[-1,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"special-pole", "zeta-exact-pole", "zeta[1]", BoundaryOutcome::DomainError},
        {"special-pole", "zeta-finite-pole", "N[zeta[N[1,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"special-pole", "Ei-exact-singularity", "Ei[0]", BoundaryOutcome::DomainError},
        {"special-pole", "Ei-finite-singularity", "N[Ei[N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"special-pole", "Ci-exact-singularity", "Ci[0]", BoundaryOutcome::DomainError},
        {"special-pole", "Ci-finite-singularity", "N[Ci[N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"special-pole", "li-exact-singularity", "li[1]", BoundaryOutcome::DomainError},
        {"special-pole", "li-finite-singularity", "N[li[N[1,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},

        // Parameter poles/domain boundaries.
        {"parameter-boundary", "1F1-denominator-pole", "N[hypergeometric1F1[1,N[0,{q}],2],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"parameter-boundary", "2F1-denominator-pole", "N[hypergeometric2F1[1,2,N[0,{q}],2],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"parameter-boundary", "ibeta-positive-domain-boundary", "N[ibeta[N[0,{q}],1,1/2],{p}]", BoundaryOutcome::PrecisionInsufficient},

        // Branch/backend classification must happen before bounded-work fallback.
        {"special-branch", "2F1-cut-ambiguity", "N[hypergeometric2F1[1/2,1/3,5/4,2+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"special-branch", "li-cut-ambiguity", "N[li[-2+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"special-branch", "polylog-cut-ambiguity", "N[polylog[2,2+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},

        // Existing mathematical values for which the current certified backend is intentionally unavailable.
        {"backend-boundary", "lambert-nonreal-branch", "N[lambertw[2,1],{p}]", BoundaryOutcome::BackendUnsupported},
        {"backend-boundary", "polylog-finite-order", "N[polylog[N[2,{q}],1/2],{p}]", BoundaryOutcome::BackendUnsupported},
        {"backend-boundary", "ibeta-finite-parameter", "N[ibeta[N[1,{q}],1,1/2],{p}]", BoundaryOutcome::BackendUnsupported},
        {"backend-boundary", "elliptic-complex-input", "N[ellipticF[1/2+I/10,1/3],{p}]", BoundaryOutcome::BackendUnsupported},
        {"backend-boundary", "elliptic-real-complex-boundary", "N[ellipticF[1/2+I*N[0,{q}],1/3],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"backend-boundary", "elliptic-series-threshold", "N[ellipticF[1/2,19/20],{p}]", BoundaryOutcome::BackendUnsupported},
        {"backend-boundary", "Ci-series-threshold", "N[Ci[97],{p}]", BoundaryOutcome::BackendUnsupported},
        {"backend-boundary", "Ci-threshold-crossing", "N[Ci[N[96,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient}
    }};
}

[[nodiscard]] std::string_view outcomeName(BoundaryOutcome outcome) noexcept {
    switch (outcome) {
    case BoundaryOutcome::Value: return "Value";
    case BoundaryOutcome::DomainError: return "DomainError";
    case BoundaryOutcome::PrecisionInsufficient: return "PrecisionInsufficient";
    case BoundaryOutcome::BackendUnsupported: return "CertifiedBackendUnsupported";
    case BoundaryOutcome::Unevaluated: return "Unevaluated";
    case BoundaryOutcome::ResourceLimit: return "ResourceLimit";
    case BoundaryOutcome::Timeout: return "Timeout";
    case BoundaryOutcome::OtherError: return "OtherError";
    case BoundaryOutcome::ConflictingDiagnostics: return "ConflictingDiagnostics";
    }
    return "Unknown";
}

[[nodiscard]] bool hasWarning(
    const kernel::KernelSession& session,
    std::string_view code) {
    return std::ranges::any_of(session.diagnostics(), [&](const auto& diagnostic) {
        return diagnostic.severity == evaluation::DiagnosticSeverity::Warning
            && diagnostic.code == code;
    });
}

[[nodiscard]] std::string diagnosticsText(const kernel::KernelSession& session) {
    std::ostringstream stream;
    bool first = true;
    for (const auto& diagnostic : session.diagnostics()) {
        if (!first)
            stream << "; ";
        first = false;
        stream << diagnostic.code << ": " << diagnostic.message;
    }
    return stream.str();
}

[[nodiscard]] BoundaryOutcome classifySuccessfulEvaluation(
    const kernel::KernelSession& session,
    const expression::Expr& result) {
    const bool precision = hasWarning(session, "N::precision");
    const bool unsupported = hasWarning(session, "N::unsupported");
    const bool unevaluated = hasWarning(session, "N::unevaluated");
    const unsigned classificationCount = static_cast<unsigned>(precision)
        + static_cast<unsigned>(unsupported) + static_cast<unsigned>(unevaluated);
    if (classificationCount > 1)
        return BoundaryOutcome::ConflictingDiagnostics;
    if (precision)
        return BoundaryOutcome::PrecisionInsufficient;
    if (unsupported)
        return BoundaryOutcome::BackendUnsupported;
    if (unevaluated)
        return BoundaryOutcome::Unevaluated;
    if (result.isNumber() || result.isDecimalApproximation()
        || result.isComplexDecimalApproximation())
        return BoundaryOutcome::Value;
    return BoundaryOutcome::Unevaluated;
}

[[nodiscard]] CaseObservation observe(const GeneratedCase& testCase, std::uint64_t timeoutMilliseconds) {
    kernel::KernelSession session;
    auto limits = session.evaluationLimits();
    // A branch/pole ambiguity that cannot improve with more guard digits should stop locally.
    // Keep enough room for legitimate special-function work while making refinement cliffs visible.
    limits.maxCertifiedRefinements = 20'000;
    session.setEvaluationLimits(limits);

    evaluation::EvaluationCancellationToken cancellation;
    std::mutex watchdogMutex;
    std::condition_variable watchdogCondition;
    bool finished = false;
    std::atomic<bool> timedOut{false};
    std::thread watchdog{[&] {
        std::unique_lock lock{watchdogMutex};
        if (!watchdogCondition.wait_for(
                lock,
                std::chrono::milliseconds{timeoutMilliseconds},
                [&] { return finished; })) {
            timedOut.store(true, std::memory_order_relaxed);
            cancellation.requestCancellation();
        }
    }};
    const auto stopWatchdog = [&] {
        {
            std::lock_guard lock{watchdogMutex};
            finished = true;
        }
        watchdogCondition.notify_one();
        watchdog.join();
    };

    const auto start = Clock::now();
    try {
        const expression::Expr result = session.evaluate(testCase.source, cancellation);
        stopWatchdog();
        const auto end = Clock::now();
        return CaseObservation{
            classifySuccessfulEvaluation(session, result),
            formatting::formatExpr(result),
            diagnosticsText(session),
            session.lastEvaluationUsage(),
            std::chrono::duration<double, std::milli>(end - start).count()};
    }
    catch (const error::CalcError& exception) {
        stopWatchdog();
        const auto end = Clock::now();
        BoundaryOutcome outcome = BoundaryOutcome::OtherError;
        if (timedOut.load(std::memory_order_relaxed))
            outcome = BoundaryOutcome::Timeout;
        else if (exception.type() == error::CalcErrorType::Domain)
            outcome = BoundaryOutcome::DomainError;
        else if (exception.type() == error::CalcErrorType::ResourceLimit)
            outcome = BoundaryOutcome::ResourceLimit;
        return CaseObservation{
            outcome,
            {},
            std::string{error::calcErrorTypeName(exception.type())} + ": " + exception.what(),
            session.lastEvaluationUsage(),
            std::chrono::duration<double, std::milli>(end - start).count()};
    }
    catch (const std::exception& exception) {
        stopWatchdog();
        const auto end = Clock::now();
        return CaseObservation{
            timedOut.load(std::memory_order_relaxed) ? BoundaryOutcome::Timeout : BoundaryOutcome::OtherError,
            {},
            exception.what(),
            session.lastEvaluationUsage(),
            std::chrono::duration<double, std::milli>(end - start).count()};
    }
}

[[nodiscard]] GeneratedCase generateCase(std::uint64_t masterSeed, std::uint64_t index) {
    std::mt19937_64 rng{caseSeed(masterSeed, index)};
    constexpr std::array<std::size_t, 5> requestedDigits{2, 5, 20, 50, 100};
    constexpr std::array<std::size_t, 4> inputDigits{2, 5, 10, 20};
    constexpr std::array<std::size_t, 5> epsilonPowers{12, 24, 40, 80, 160};
    constexpr auto probes = probeTemplates();

    std::uniform_int_distribution<std::size_t> probeDistribution(0, probes.size() - 1);
    const auto& probe = probes[probeDistribution(rng)];
    const std::size_t p = choose(rng, requestedDigits);
    const std::size_t q = choose(rng, inputDigits);
    const std::size_t k = std::max(choose(rng, epsilonPowers), p + 4);
    return GeneratedCase{
        std::string{probe.family},
        std::string{probe.name},
        instantiate(probe.source, p, q, k),
        probe.expected,
        p,
        q,
        k};
}

[[nodiscard]] std::string budgetText(const evaluation::EvaluationUsage& usage) {
    std::ostringstream stream;
    stream << "steps=" << usage.evaluationSteps
           << " depth=" << usage.maximumDepth
           << " nodes=" << usage.generatedNodes
           << " simplify=" << usage.simplificationCandidates
           << " refine=" << usage.certifiedRefinements
           << " precision-digits=" << usage.maximumRequestedPrecisionDigits;
    return stream.str();
}

[[nodiscard]] std::optional<Failure> checkCase(
    std::uint64_t masterSeed,
    std::uint64_t index,
    std::uint64_t timeoutMilliseconds) {
    GeneratedCase testCase = generateCase(masterSeed, index);
    CaseObservation observation = observe(testCase, timeoutMilliseconds);

    if (observation.outcome != testCase.expected) {
        return Failure{
            index, std::move(testCase), std::move(observation),
            "classification mismatch"};
    }

    if (testCase.expected == BoundaryOutcome::PrecisionInsufficient
        && observation.usage.certifiedRefinements > kPrecisionRefinementCeiling) {
        return Failure{
            index, std::move(testCase), std::move(observation),
            "persistent information ambiguity consumed excessive certified refinements"};
    }

    return std::nullopt;
}

void printFailure(std::uint64_t masterSeed, const Failure& failure) {
    std::cerr
        << "\n[FAIL] Certification boundary invariant\n"
        << "Seed      : " << masterSeed << '\n'
        << "Case      : " << failure.index << '\n'
        << "Case seed : " << caseSeed(masterSeed, failure.index) << '\n'
        << "Family    : " << failure.testCase.family << '\n'
        << "Probe     : " << failure.testCase.name << '\n'
        << "Expression: " << failure.testCase.source << '\n'
        << "Expected  : " << outcomeName(failure.testCase.expected) << '\n'
        << "Actual    : " << outcomeName(failure.observation.outcome) << '\n'
        << "Reason    : " << failure.reason << '\n';
    if (!failure.observation.formattedResult.empty())
        std::cerr << "Result    : " << failure.observation.formattedResult << '\n';
    if (!failure.observation.detail.empty())
        std::cerr << "Detail    : " << failure.observation.detail << '\n';
    std::cerr
        << "Elapsed   : " << failure.observation.elapsedMilliseconds << " ms\n"
        << "Budget    : " << budgetText(failure.observation.usage) << '\n'
        << "Reproduce : mmCal.Benchmarks --certification-boundaries --seed "
        << masterSeed << " --case " << failure.index << '\n';
}

[[nodiscard]] bool printSingleCase(
    std::uint64_t masterSeed,
    std::uint64_t index,
    std::uint64_t timeoutMilliseconds) {
    GeneratedCase testCase = generateCase(masterSeed, index);
    std::cout
        << "Certification boundary case\n"
        << "Seed      : " << masterSeed << '\n'
        << "Case      : " << index << '\n'
        << "Case seed : " << caseSeed(masterSeed, index) << '\n'
        << "Family    : " << testCase.family << '\n'
        << "Probe     : " << testCase.name << '\n'
        << "Expression: " << testCase.source << '\n'
        << "Expected  : " << outcomeName(testCase.expected) << '\n'
        << std::flush;
    CaseObservation observation = observe(testCase, timeoutMilliseconds);
    std::cout << "Actual    : " << outcomeName(observation.outcome) << '\n';
    if (!observation.formattedResult.empty())
        std::cout << "Result    : " << observation.formattedResult << '\n';
    if (!observation.detail.empty())
        std::cout << "Detail    : " << observation.detail << '\n';
    std::cout
        << "Elapsed   : " << observation.elapsedMilliseconds << " ms\n"
        << "Budget    : " << budgetText(observation.usage) << '\n';

    const bool excessiveRefinement = testCase.expected == BoundaryOutcome::PrecisionInsufficient
        && observation.usage.certifiedRefinements > kPrecisionRefinementCeiling;
    const bool passed = observation.outcome == testCase.expected && !excessiveRefinement;
    if (passed)
        std::cout << "Status    : PASS\n";
    else if (excessiveRefinement)
        std::cout << "Status    : FAIL (excessive certified refinement)\n";
    else
        std::cout << "Status    : FAIL\n";
    if (!passed)
        std::cout << "Reproduce : mmCal.Benchmarks --certification-boundaries --seed "
                  << masterSeed << " --case " << index << '\n';
    return passed;
}

} // namespace

bool runCertificationBoundaryFuzzer(const CertificationBoundaryFuzzerOptions& options) {
    if (options.singleCase)
        return printSingleCase(
            options.seed, *options.singleCase, options.caseTimeoutMilliseconds);

    std::cout
        << "Certification boundary fuzzing\n"
        << "Seed      : " << options.seed << '\n'
        << "Mode      : " << (options.loop
            ? (options.noStopLoop ? "infinite (FAIL is reported and fuzzing continues)" : "infinite (stop on first FAIL)")
            : "finite") << '\n'
        << "Threads   : " << options.threads << '\n'
        << "Probes    : " << probeTemplates().size() << '\n'
        << "Timeout   : " << options.caseTimeoutMilliseconds << " ms/case\n";
    if (!options.loop)
        std::cout << "Cases     : " << options.cases << '\n';

    const auto start = Clock::now();
    std::atomic<std::uint64_t> nextIndex{1};
    std::atomic<std::uint64_t> completed{0};
    std::atomic<std::uint64_t> failures{0};
    std::atomic<bool> stop{false};
    std::mutex outputMutex;

    const auto worker = [&] {
        while (!stop.load(std::memory_order_relaxed)) {
            const std::uint64_t index = nextIndex.fetch_add(1, std::memory_order_relaxed);
            if (!options.loop && index > options.cases)
                break;

            if (auto failure = checkCase(options.seed, index, options.caseTimeoutMilliseconds)) {
                failures.fetch_add(1, std::memory_order_relaxed);
                {
                    std::lock_guard lock{outputMutex};
                    printFailure(options.seed, *failure);
                }
                if (!options.noStopLoop) {
                    stop.store(true, std::memory_order_relaxed);
                    break;
                }
            }

            const std::uint64_t done = completed.fetch_add(1, std::memory_order_relaxed) + 1;
            if (options.reportEvery != 0 && done % options.reportEvery == 0) {
                const auto now = Clock::now();
                std::lock_guard lock{outputMutex};
                std::cout << '[' << done << "] PASS  "
                          << std::chrono::duration<double>(now - start).count() << " s"
                          << "  failures=" << failures.load(std::memory_order_relaxed) << '\n';
            }
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(options.threads);
    for (std::size_t i = 0; i < options.threads; ++i)
        workers.emplace_back(worker);
    for (auto& thread : workers)
        thread.join();

    const auto end = Clock::now();
    const std::uint64_t failureCount = failures.load(std::memory_order_relaxed);
    std::cout
        << "Certification boundary fuzzing: "
        << (failureCount == 0 ? "PASS" : "FAIL")
        << "  cases=" << completed.load(std::memory_order_relaxed)
        << " failures=" << failureCount
        << " elapsed=" << std::chrono::duration<double>(end - start).count() << " s\n";
    return failureCount == 0;
}

} // namespace mmcal::benchmarks
