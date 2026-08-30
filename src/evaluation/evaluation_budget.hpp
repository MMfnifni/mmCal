#pragma once

#include "error/error_message.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <csignal>
#include <string>
#include <string_view>
#include <utility>

namespace mmcal::evaluation {

// 1回の評価要求で共有する上限。局所engineが持つ、より小さいalgorithm上限は
// この上限の子budgetとして併用し、既存の停止性を失わない。
struct EvaluationLimits final {
    std::size_t maxInputBytes = 16 * 1024 * 1024;
    std::size_t maxEvaluationSteps = 10'000'000;
    std::size_t maxDepth = 1024;
    std::size_t maxGeneratedNodes = 5'000'000;
    std::size_t maxSimplificationCandidates = 2'000'000;
    std::size_t maxSolverBranches = 50'000;
    std::size_t maxIntegrationCandidates = 100'000;
    std::size_t maxCertifiedRefinements = 1'000'000;
    std::size_t maxDenseArrayElements = 10'000'000;
    std::size_t maxTemporaryMatrixElements = 25'000'000;
    std::size_t maxBigIntegerBits = 8'000'000;
    std::size_t maxRequestedPrecisionDigits = 100'000;
    std::size_t maxAlgebraicDegree = 96;
    std::size_t maxAlgebraicRefinements = 100'000;

    [[nodiscard]] bool operator==(const EvaluationLimits&) const = default;
};

// 成功・失敗のどちらでも取得できる、直前評価の決定論的telemetry。
// deadlineやwall clockはCoreへ入れず、frontend側の取消しと分離する。
struct EvaluationUsage final {
    std::size_t inputBytes = 0;
    std::size_t evaluationSteps = 0;
    std::size_t maximumDepth = 0;
    std::size_t generatedNodes = 0;
    std::size_t simplificationCandidates = 0;
    std::size_t solverBranches = 0;
    std::size_t integrationCandidates = 0;
    std::size_t certifiedRefinements = 0;
    std::size_t denseArrayElements = 0;
    std::size_t temporaryMatrixElements = 0;
    std::size_t maximumBigIntegerBits = 0;
    std::size_t maximumRequestedPrecisionDigits = 0;
    std::size_t maximumAlgebraicDegree = 0;
    std::size_t algebraicRefinements = 0;
    std::size_t modularPrimes = 0;

    [[nodiscard]] bool operator==(const EvaluationUsage&) const = default;
};

enum class EvaluationResource {
    EvaluationStep,
    GeneratedNode,
    SimplificationCandidate,
    SolverBranch,
    IntegrationCandidate,
    CertifiedRefinement,
    DenseArrayElement,
    TemporaryMatrixElement,
    AlgebraicRefinement
};

// wall clockをCoreへ持ち込まず，deadline・UI cancel・process policyをfrontend側で
// 実装するための共有token。token自体はbudgetを所有せず，要求終了後も再利用できる。
class EvaluationCancellationToken final {
public:
    void requestCancellation() noexcept {
        requested_.store(true, std::memory_order_relaxed);
    }

    // POSIX signal handlerのようなasync-signal contextではstd::atomic操作を要求せず，
    // sig_atomic_tへの単純storeだけでcancel要求を伝える。通常threadからは上のAPIを使う。
    void requestSignalCancellation() noexcept {
        signalRequested_ = 1;
    }

    void reset() noexcept {
        requested_.store(false, std::memory_order_relaxed);
        signalRequested_ = 0;
    }

    [[nodiscard]] bool cancellationRequested() const noexcept {
        return requested_.load(std::memory_order_relaxed) || signalRequested_ != 0;
    }

private:
    std::atomic<bool> requested_{false};
    volatile std::sig_atomic_t signalRequested_ = 0;
};

class EvaluationBudget final {
public:
    explicit EvaluationBudget(
        EvaluationLimits limits = {},
        const EvaluationCancellationToken* cancellation = nullptr) noexcept
        : limits_(std::move(limits)), cancellation_(cancellation) {}

    [[nodiscard]] const EvaluationLimits& limits() const noexcept {
        return limits_;
    }

    [[nodiscard]] const EvaluationUsage& usage() const noexcept {
        return usage_;
    }

    void consume(EvaluationResource resource, std::size_t amount = 1) {
        checkCancellation();
        switch (resource) {
        case EvaluationResource::EvaluationStep:
            consumeCounter(usage_.evaluationSteps, amount,
                limits_.maxEvaluationSteps, "Evaluation step budget");
            return;
        case EvaluationResource::GeneratedNode:
            consumeCounter(usage_.generatedNodes, amount,
                limits_.maxGeneratedNodes, "Generated expression node budget");
            return;
        case EvaluationResource::SimplificationCandidate:
            consumeCounter(usage_.simplificationCandidates, amount,
                limits_.maxSimplificationCandidates, "Simplification candidate budget");
            return;
        case EvaluationResource::SolverBranch:
            consumeCounter(usage_.solverBranches, amount,
                limits_.maxSolverBranches, "Solver branch budget");
            return;
        case EvaluationResource::IntegrationCandidate:
            consumeCounter(usage_.integrationCandidates, amount,
                limits_.maxIntegrationCandidates, "Integration candidate budget");
            return;
        case EvaluationResource::CertifiedRefinement:
            consumeCounter(usage_.certifiedRefinements, amount,
                limits_.maxCertifiedRefinements, "Certified refinement budget");
            return;
        case EvaluationResource::DenseArrayElement:
            consumeCounter(usage_.denseArrayElements, amount,
                limits_.maxDenseArrayElements, "Dense array element budget");
            return;
        case EvaluationResource::TemporaryMatrixElement:
            consumeCounter(usage_.temporaryMatrixElements, amount,
                limits_.maxTemporaryMatrixElements, "Temporary matrix element budget");
            return;
        case EvaluationResource::AlgebraicRefinement:
            consumeCounter(usage_.algebraicRefinements, amount,
                limits_.maxAlgebraicRefinements, "Algebraic refinement budget");
            return;
        }
    }

    void checkDepth(std::size_t depth) {
        checkCancellation();
        usage_.maximumDepth = std::max(usage_.maximumDepth, depth);
        checkMaximum(depth, limits_.maxDepth, "Evaluation depth budget");
    }

    void checkBigIntegerBits(std::size_t bits) {
        checkCancellation();
        usage_.maximumBigIntegerBits = std::max(usage_.maximumBigIntegerBits, bits);
        checkMaximum(bits, limits_.maxBigIntegerBits, "BigInt bit-length budget");
    }

    void checkRequestedPrecisionDigits(std::size_t digits) {
        checkCancellation();
        usage_.maximumRequestedPrecisionDigits =
            std::max(usage_.maximumRequestedPrecisionDigits, digits);
        checkMaximum(digits, limits_.maxRequestedPrecisionDigits,
            "Requested precision budget");
    }

    void checkAlgebraicDegree(std::size_t degree) {
        checkCancellation();
        usage_.maximumAlgebraicDegree = std::max(usage_.maximumAlgebraicDegree, degree);
        checkMaximum(degree, limits_.maxAlgebraicDegree,
            "Algebraic construction degree budget");
    }

    void checkInputBytes(std::size_t bytes) {
        checkCancellation();
        usage_.inputBytes = bytes;
        checkMaximum(bytes, limits_.maxInputBytes, "Frontend input byte budget");
    }

    void pollCancellation() const {
        checkCancellation();
    }

    void recordModularPrime() {
        checkCancellation();
        ++usage_.modularPrimes;
    }

private:
    EvaluationLimits limits_;
    EvaluationUsage usage_;
    const EvaluationCancellationToken* cancellation_ = nullptr;

    void checkCancellation() const {
        if (cancellation_ && cancellation_->cancellationRequested())
            error::throwCalcError(
                error::CalcErrorType::ResourceLimit,
                "Evaluation cancelled by frontend");
    }

    [[noreturn]] static void exceeded(std::string_view name, std::size_t limit) {
        error::throwCalcError(
            error::CalcErrorType::ResourceLimit,
            std::string{name} + " exceeded (limit " + std::to_string(limit) + ')');
    }

    static void consumeCounter(
        std::size_t& used,
        std::size_t amount,
        std::size_t limit,
        std::string_view name) {
        if (amount > limit || used > limit - amount)
            exceeded(name, limit);
        used += amount;
    }

    static void checkMaximum(
        std::size_t value,
        std::size_t limit,
        std::string_view name) {
        if (value > limit)
            exceeded(name, limit);
    }
};

// 深いsubsystemへ引数を延々追加せず同じ要求budgetを渡すための、非所有scope。
// thread_localなので並列KernelSession間では共有されず、nested scope終了時に必ず復元する。
inline thread_local EvaluationBudget* activeEvaluationBudget = nullptr;

class EvaluationBudgetScope final {
public:
    explicit EvaluationBudgetScope(EvaluationBudget* budget) noexcept
        : previous_(std::exchange(activeEvaluationBudget, budget)) {}

    EvaluationBudgetScope(const EvaluationBudgetScope&) = delete;
    EvaluationBudgetScope& operator=(const EvaluationBudgetScope&) = delete;

    ~EvaluationBudgetScope() {
        activeEvaluationBudget = previous_;
    }

private:
    EvaluationBudget* previous_ = nullptr;
};

[[nodiscard]] inline EvaluationBudget* currentEvaluationBudget() noexcept {
    return activeEvaluationBudget;
}

inline void consumeEvaluationBudget(
    EvaluationResource resource,
    std::size_t amount = 1) {
    if (activeEvaluationBudget)
        activeEvaluationBudget->consume(resource, amount);
}

inline void checkEvaluationCancellation() {
    if (activeEvaluationBudget)
        activeEvaluationBudget->pollCancellation();
}

inline void recordEvaluationModularPrime() {
    if (activeEvaluationBudget)
        activeEvaluationBudget->recordModularPrime();
}

inline void checkEvaluationBigIntegerBits(std::size_t bits) {
    if (activeEvaluationBudget)
        activeEvaluationBudget->checkBigIntegerBits(bits);
}

inline void checkEvaluationRequestedPrecisionDigits(std::size_t digits) {
    if (activeEvaluationBudget)
        activeEvaluationBudget->checkRequestedPrecisionDigits(digits);
}

inline void checkEvaluationAlgebraicDegree(std::size_t degree) {
    if (activeEvaluationBudget)
        activeEvaluationBudget->checkAlgebraicDegree(degree);
}

} // namespace mmcal::evaluation
