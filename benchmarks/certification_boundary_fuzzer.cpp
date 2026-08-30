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
#include <optional>
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

enum class ProbeKind {
    Classification,
    Metamorphic
};

enum class MetamorphicRelation {
    SameValue,
    InformationDoesNotNarrow,
    InformationContainsExactRoute,
    StablePrecisionInsufficient,
    SameOutcome
};

struct ProbeTemplate final {
    std::string_view family;
    std::string_view name;
    std::string_view source;
    BoundaryOutcome expected = BoundaryOutcome::Value;
};

struct MetamorphicTemplate final {
    std::string_view family;
    std::string_view name;
    std::string_view baselineSource;
    std::string_view variantSource;
    BoundaryOutcome expected = BoundaryOutcome::Value;
    MetamorphicRelation relation = MetamorphicRelation::SameValue;
};

struct GeneratedCase final {
    ProbeKind kind = ProbeKind::Classification;
    std::string family;
    std::string name;
    std::string source;
    std::string variantSource;
    BoundaryOutcome expected = BoundaryOutcome::Value;
    MetamorphicRelation relation = MetamorphicRelation::SameValue;
    std::size_t requestedDigits = 20;
    std::size_t inputDigits = 5;
    std::size_t epsilonPower = 40;
};

struct CaseObservation final {
    BoundaryOutcome outcome = BoundaryOutcome::OtherError;
    std::optional<expression::Expr> result;
    std::string formattedResult;
    std::string detail;
    evaluation::EvaluationUsage usage{};
    double elapsedMilliseconds = 0.0;
};

struct Failure final {
    std::uint64_t index = 0;
    GeneratedCase testCase;
    CaseObservation observation;
    std::optional<CaseObservation> variantObservation;
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
    const std::size_t lowDigits = std::min(requestedDigits, inputDigits);
    const std::size_t highDigits = std::max(requestedDigits, inputDigits + 1);
    std::string result{source};
    result = replaceAll(std::move(result), "{p}", std::to_string(requestedDigits));
    result = replaceAll(std::move(result), "{q}", std::to_string(inputDigits));
    result = replaceAll(std::move(result), "{lo}", std::to_string(lowDigits));
    result = replaceAll(std::move(result), "{hi}", std::to_string(highDigits));
    result = replaceAll(std::move(result), "{k}", std::to_string(epsilonPower));
    return result;
}

[[nodiscard]] constexpr std::array<ProbeTemplate, 59> probeTemplates() {
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
        {"special-branch", "Ei-large-cut-ambiguity", "N[Ei[-1000+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"special-branch", "Ci-large-cut-ambiguity", "N[Ci[-1000+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"backend-extension", "Ei-positive-axis-finite-imaginary", "N[Ei[1000+I*N[0,{q}]],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "Ci-pure-imaginary-finite-input", "N[Ci[I*N[513,{q}]],{p}]", BoundaryOutcome::Value},

        // Existing mathematical values for which the current certified backend is intentionally unavailable.
        {"backend-extension", "lambert-nonreal-branch", "N[lambertw[2,1],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "lambert-principal-branchpoint-upper", "N[lambertw[-1/E+I/10^{k}],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "lambert-principal-branchpoint-lower", "N[lambertw[-1/E-I/10^{k}],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "lambert-minus-one-branchpoint-upper", "N[lambertw[-1,-1/E+I/10^{k}],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "lambert-plus-one-branchpoint-lower", "N[lambertw[1,-1/E-I/10^{k}],{p}]", BoundaryOutcome::Value},
        {"special-branch", "lambert-branchpoint-cut-ambiguity", "N[lambertw[-1,-1/E+I*N[0,{q}]],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"backend-boundary", "polylog-finite-order", "N[polylog[N[2,{q}],1/2],{p}]", BoundaryOutcome::BackendUnsupported},
        {"backend-extension", "ibeta-finite-parameter", "N[ibeta[N[1,{q}],1,1/2],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "zeta-critical-strip-complex", "N[zeta[1+I],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "zeta-critical-strip-finite-realpart", "N[zeta[N[1,{q}]+I],{p}]", BoundaryOutcome::Value},
        {"backend-boundary", "elliptic-complex-input", "N[ellipticF[1/2+I/10,1/3],{p}]", BoundaryOutcome::BackendUnsupported},
        {"backend-boundary", "elliptic-real-complex-boundary", "N[ellipticF[1/2+I*N[0,{q}],1/3],{p}]", BoundaryOutcome::PrecisionInsufficient},
        {"backend-extension", "elliptic-former-series-threshold", "N[ellipticF[1/2,19/20],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "Ci-former-series-threshold", "N[Ci[97],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "Ci-former-threshold-finite-input", "N[Ci[N[96,{q}]],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "Ei-former-complex-threshold", "N[Ei[513I],{p}]", BoundaryOutcome::Value},
        {"backend-extension", "Ci-former-complex-threshold", "N[Ci[140+I],{p}]", BoundaryOutcome::Value}
    }};
}

[[nodiscard]] constexpr std::array<MetamorphicTemplate, 36> metamorphicTemplates() {
    using Relation = MetamorphicRelation;
    return {{
        // Increasing the outer requested precision must never recover information hidden by the inner N.
        {"N-provenance", "nested-N-pi-no-gain",
            "N[Pi,{lo}]", "N[N[Pi,{lo}],{hi}]", BoundaryOutcome::Value, Relation::SameValue},
        {"N-provenance", "nested-N-rational-no-gain",
            "N[1/3,{lo}]", "N[N[1/3,{lo}],{hi}]", BoundaryOutcome::Value, Relation::SameValue},
        {"N-provenance", "nested-N-log-no-gain",
            "N[log[2],{lo}]", "N[N[log[2],{lo}],{hi}]", BoundaryOutcome::Value, Relation::SameValue},
        {"N-provenance", "nested-N-complex-no-gain",
            "N[1+I/3,{lo}]", "N[N[1+I/3,{lo}],{hi}]", BoundaryOutcome::Value, Relation::SameValue},

        // Reducing display precision may discard information, but may not make the InformationEnclosure narrower.
        {"N-provenance", "nested-N-pi-degrades-only",
            "N[Pi,{hi}]", "N[N[Pi,{hi}],{lo}]", BoundaryOutcome::Value, Relation::InformationDoesNotNarrow},
        {"N-provenance", "nested-N-log-degrades-only",
            "N[log[2],{hi}]", "N[N[log[2],{hi}],{lo}]", BoundaryOutcome::Value, Relation::InformationDoesNotNarrow},
        {"N-provenance", "nested-N-complex-degrades-only",
            "N[1+I/3,{hi}]", "N[N[1+I/3,{hi}],{lo}]", BoundaryOutcome::Value, Relation::InformationDoesNotNarrow},

        // Exact identities must not requantize or otherwise rewrite finite-precision provenance.
        {"N-provenance", "exact-zero-addition-preserves-value",
            "N[Pi,{q}]", "N[Pi,{q}]+0", BoundaryOutcome::Value, Relation::SameValue},
        {"N-provenance", "exact-one-multiplication-preserves-value",
            "N[I/10^{k},{q}]", "N[I/10^{k},{q}]*1", BoundaryOutcome::Value, Relation::SameValue},
        {"N-provenance", "double-negation-preserves-value",
            "N[1+I/3,{q}]", "-(-N[1+I/3,{q}])", BoundaryOutcome::Value, Relation::SameValue},

        // Persistent finite-information ambiguity must survive harmless rewrites and outer N refinement.
        {"ambiguity-stability", "log-cut-survives-outer-N",
            "N[log[-1+I*N[0,{lo}]],{hi}]",
            "N[log[-1+I*N[N[0,{lo}],{hi}]],{hi}]",
            BoundaryOutcome::PrecisionInsufficient, Relation::StablePrecisionInsufficient},
        {"ambiguity-stability", "sqrt-cut-survives-exact-zero-addition",
            "N[sqrt[-1+I*N[0,{lo}]],{hi}]",
            "N[sqrt[-1+I*(N[0,{lo}]+0)],{hi}]",
            BoundaryOutcome::PrecisionInsufficient, Relation::StablePrecisionInsufficient},
        {"ambiguity-stability", "power-cut-survives-exact-one-multiplication",
            "N[(-1+I*N[0,{lo}])^(1/3),{hi}]",
            "N[(-1+I*(N[0,{lo}]*1))^(1/3),{hi}]",
            BoundaryOutcome::PrecisionInsufficient, Relation::StablePrecisionInsufficient},
        {"ambiguity-stability", "gamma-pole-survives-outer-N",
            "N[gamma[N[0,{lo}]],{hi}]",
            "N[gamma[N[N[0,{lo}],{hi}]],{hi}]",
            BoundaryOutcome::PrecisionInsufficient, Relation::StablePrecisionInsufficient},
        {"ambiguity-stability", "zeta-pole-survives-outer-N",
            "N[zeta[N[1,{lo}]],{hi}]",
            "N[zeta[N[N[1,{lo}],{hi}]],{hi}]",
            BoundaryOutcome::PrecisionInsufficient, Relation::StablePrecisionInsufficient},
        {"ambiguity-stability", "li-pole-survives-outer-N",
            "N[li[N[1,{lo}]],{hi}]",
            "N[li[N[N[1,{lo}],{hi}]],{hi}]",
            BoundaryOutcome::PrecisionInsufficient, Relation::StablePrecisionInsufficient},
        {"ambiguity-stability", "atan2-origin-survives-outer-N",
            "N[atan2[0,N[0,{lo}]],{hi}]",
            "N[atan2[0,N[N[0,{lo}],{hi}]],{hi}]",
            BoundaryOutcome::PrecisionInsufficient, Relation::StablePrecisionInsufficient},
        {"ambiguity-stability", "2F1-cut-survives-outer-N",
            "N[hypergeometric2F1[1/2,1/3,5/4,2+I*N[0,{lo}]],{hi}]",
            "N[hypergeometric2F1[1/2,1/3,5/4,2+I*N[N[0,{lo}],{hi}]],{hi}]",
            BoundaryOutcome::PrecisionInsufficient, Relation::StablePrecisionInsufficient},
        {"ambiguity-stability", "elliptic-boundary-survives-outer-N",
            "N[ellipticF[1/2+I*N[0,{lo}],1/3],{hi}]",
            "N[ellipticF[1/2+I*N[N[0,{lo}],{hi}],1/3],{hi}]",
            BoundaryOutcome::PrecisionInsufficient, Relation::StablePrecisionInsufficient},

        // Linear algebra must not recover a hidden exact zero from a nested finite-precision leaf.
        {"linear-algebra-provenance", "rank-finite-zero-stays-undetermined",
            "matrixRank[{{N[1,{lo}],0},{0,N[0,{lo}]}}]",
            "matrixRank[{{N[1,{lo}],0},{0,N[N[0,{lo}],{hi}]}}]",
            BoundaryOutcome::Unevaluated, Relation::SameOutcome},
        {"linear-algebra-provenance", "nullspace-finite-zero-stays-undetermined",
            "nullSpace[{{N[1,{lo}],0},{0,N[0,{lo}]}}]",
            "nullSpace[{{N[1,{lo}],0},{0,N[N[0,{lo}],{hi}]}}]",
            BoundaryOutcome::Unevaluated, Relation::SameOutcome},
        {"linear-algebra-provenance", "rref-finite-zero-stays-undetermined",
            "rref[{{N[1,{lo}],0},{0,N[0,{lo}]}}]",
            "rref[{{N[1,{lo}],0},{0,N[N[0,{lo}],{hi}]}}]",
            BoundaryOutcome::Unevaluated, Relation::SameOutcome},
        {"linear-algebra-provenance", "solve-finite-zero-stays-undetermined",
            "solveLinear[{{N[1,{lo}],0},{0,N[0,{lo}]}},{1,0}]",
            "solveLinear[{{N[1,{lo}],0},{0,N[N[0,{lo}],{hi}]}},{1,0}]",
            BoundaryOutcome::Unevaluated, Relation::SameOutcome},
        {"linear-algebra-provenance", "inverse-finite-zero-stays-undetermined",
            "inverse[{{N[1,{lo}],0},{0,N[0,{lo}]}}]",
            "inverse[{{N[1,{lo}],0},{0,N[N[0,{lo}],{hi}]}}]",
            BoundaryOutcome::Unevaluated, Relation::SameOutcome},

        // f[N[x,p]] may lose information relative to N[f[x],p], but its InformationEnclosure
        // must still contain the exact-route certified truth.
        {"function-provenance", "sqrt-finite-route-contains-exact-route",
            "sqrt[N[2,{q}]]", "N[sqrt[2],{q}]",
            BoundaryOutcome::Value, Relation::InformationContainsExactRoute},
        {"function-provenance", "log-finite-route-contains-exact-route",
            "log[N[2,{q}]]", "N[log[2],{q}]",
            BoundaryOutcome::Value, Relation::InformationContainsExactRoute},
        {"function-provenance", "exp-finite-route-contains-exact-route",
            "exp[N[1/3,{q}]]", "N[exp[1/3],{q}]",
            BoundaryOutcome::Value, Relation::InformationContainsExactRoute},
        {"function-provenance", "sin-finite-route-contains-exact-route",
            "sin[N[1,{q}]]", "N[sin[1],{q}]",
            BoundaryOutcome::Value, Relation::InformationContainsExactRoute},
        {"function-provenance", "gamma-finite-route-contains-exact-route",
            "gamma[N[1/3,{q}]]", "N[gamma[1/3],{q}]",
            BoundaryOutcome::Value, Relation::InformationContainsExactRoute},
        {"function-provenance", "zeta-finite-route-contains-exact-route",
            "zeta[N[2,{q}]]", "N[zeta[2],{q}]",
            BoundaryOutcome::Value, Relation::InformationContainsExactRoute},

        // Nested N on array leaves and numerical-calculus points must be semantically inert.
        {"array-provenance", "fft-nested-N-preserves-output",
            "fft[{N[1,{lo}],N[-1+1/10^{k},{lo}]}]",
            "fft[{N[N[1,{lo}],{hi}],N[N[-1+1/10^{k},{lo}],{hi}]}]",
            BoundaryOutcome::Value, Relation::SameValue},
        {"array-provenance", "dot-nested-N-preserves-output",
            "dot[{N[1/3,{lo}],N[2/3,{lo}]},{2,3}]",
            "dot[{N[N[1/3,{lo}],{hi}],N[N[2/3,{lo}],{hi}]},{2,3}]",
            BoundaryOutcome::Value, Relation::SameValue},
        {"numerical-calculus-provenance", "diff-nested-point-preserves-output",
            "diff[x^2,x,N[1,{lo}],{hi}]",
            "diff[x^2,x,N[N[1,{lo}],{hi}],{hi}]",
            BoundaryOutcome::Value, Relation::SameValue},
        {"numerical-calculus-provenance", "nintegrate-nested-bound-preserves-output",
            "nintegrate[1,{x,N[1,{lo}],2},{hi}]",
            "nintegrate[1,{x,N[N[1,{lo}],{hi}],2},{hi}]",
            BoundaryOutcome::Value, Relation::SameValue},
        {"linear-algebra-provenance", "det-nested-N-preserves-finite-zero-information",
            "det[{{N[1,{lo}],0},{0,N[0,{lo}]}}]",
            "det[{{N[N[1,{lo}],{hi}],0},{0,N[N[0,{lo}],{hi}]}}]",
            BoundaryOutcome::Value, Relation::SameValue},

        // An existing finite parameter remains a finite parameter after outer N; it must not become exact.
        {"backend-provenance", "polylog-finite-order-stays-unsupported",
            "N[polylog[N[2,{lo}],1/2],{hi}]",
            "N[polylog[N[N[2,{lo}],{hi}],1/2],{hi}]",
            BoundaryOutcome::BackendUnsupported, Relation::SameOutcome}
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

[[nodiscard]] std::string_view relationName(MetamorphicRelation relation) noexcept {
    switch (relation) {
    case MetamorphicRelation::SameValue: return "SameValue";
    case MetamorphicRelation::InformationDoesNotNarrow: return "InformationDoesNotNarrow";
    case MetamorphicRelation::InformationContainsExactRoute: return "InformationContainsExactRoute";
    case MetamorphicRelation::StablePrecisionInsufficient: return "StablePrecisionInsufficient";
    case MetamorphicRelation::SameOutcome: return "SameOutcome";
    }
    return "Unknown";
}

struct InformationBox final {
    numeric::Rational realLower;
    numeric::Rational realUpper;
    numeric::Rational imaginaryLower;
    numeric::Rational imaginaryUpper;
};

[[nodiscard]] std::optional<InformationBox> informationBox(const expression::Expr& value) {
    if (value.isDecimalApproximation()) {
        const auto& approximation = value.asDecimalApproximation();
        return InformationBox{
            approximation.informationLower(),
            approximation.informationUpper(),
            numeric::Rational{},
            numeric::Rational{}};
    }

    if (value.isComplexDecimalApproximation()) {
        const auto& approximation = value.asComplexDecimalApproximation();
        return InformationBox{
            approximation.realInformationLower(),
            approximation.realInformationUpper(),
            approximation.imaginaryInformationLower(),
            approximation.imaginaryInformationUpper()};
    }

    return std::nullopt;
}

[[nodiscard]] bool contains(const InformationBox& outer, const InformationBox& inner) {
    return outer.realLower <= inner.realLower
        && outer.realUpper >= inner.realUpper
        && outer.imaginaryLower <= inner.imaginaryLower
        && outer.imaginaryUpper >= inner.imaginaryUpper;
}

[[nodiscard]] std::optional<InformationBox> certifiedBox(const expression::Expr& value) {
    if (value.isDecimalApproximation()) {
        const auto& approximation = value.asDecimalApproximation();
        return InformationBox{
            approximation.certifiedLower(),
            approximation.certifiedUpper(),
            numeric::Rational{},
            numeric::Rational{}};
    }

    if (value.isComplexDecimalApproximation()) {
        const auto& approximation = value.asComplexDecimalApproximation();
        return InformationBox{
            approximation.real().certifiedLower(),
            approximation.real().certifiedUpper(),
            approximation.imaginary().certifiedLower(),
            approximation.imaginary().certifiedUpper()};
    }

    return std::nullopt;
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
        || result.isComplexDecimalApproximation() || result.isArray())
        return BoundaryOutcome::Value;
    return BoundaryOutcome::Unevaluated;
}

[[nodiscard]] CaseObservation observeSource(
    std::string_view source,
    std::uint64_t timeoutMilliseconds) {
    kernel::KernelSession session;
    // 正常値のcertified backendは内部級数・Euler-Maclaurin等で多数のwork unitを
    // 消費し得るため，fuzzer固有の低いglobal上限は課さない。persistent ambiguityは
    // 評価後の4096 refinement上限とfrontend timeoutで別に検出する。

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
        const expression::Expr result = session.evaluate(source, cancellation);
        stopWatchdog();
        const auto end = Clock::now();
        return CaseObservation{
            classifySuccessfulEvaluation(session, result),
            result,
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
            std::nullopt,
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
            std::nullopt,
            {},
            exception.what(),
            session.lastEvaluationUsage(),
            std::chrono::duration<double, std::milli>(end - start).count()};
    }
}

[[nodiscard]] CaseObservation observe(
    const GeneratedCase& testCase,
    std::uint64_t timeoutMilliseconds) {
    return observeSource(testCase.source, timeoutMilliseconds);
}

[[nodiscard]] GeneratedCase generateCase(std::uint64_t masterSeed, std::uint64_t index) {
    std::mt19937_64 rng{caseSeed(masterSeed, index)};
    constexpr std::array<std::size_t, 5> requestedDigits{2, 5, 20, 50, 100};
    constexpr std::array<std::size_t, 4> inputDigits{2, 5, 10, 20};
    constexpr std::array<std::size_t, 5> epsilonPowers{12, 24, 40, 80, 160};
    constexpr auto classificationProbes = probeTemplates();
    constexpr auto metamorphicProbes = metamorphicTemplates();

    const std::size_t p = choose(rng, requestedDigits);
    const std::size_t q = choose(rng, inputDigits);
    const std::size_t k = std::max(choose(rng, epsilonPowers), p + 4);
    std::uniform_int_distribution<std::size_t> probeDistribution(
        0, classificationProbes.size() + metamorphicProbes.size() - 1);
    const std::size_t selected = probeDistribution(rng);

    if (selected < classificationProbes.size()) {
        const auto& probe = classificationProbes[selected];
        GeneratedCase result;
        result.kind = ProbeKind::Classification;
        result.family = probe.family;
        result.name = probe.name;
        result.source = instantiate(probe.source, p, q, k);
        result.expected = probe.expected;
        result.requestedDigits = p;
        result.inputDigits = q;
        result.epsilonPower = k;
        return result;
    }

    const auto& probe = metamorphicProbes[selected - classificationProbes.size()];
    GeneratedCase result;
    result.kind = ProbeKind::Metamorphic;
    result.family = probe.family;
    result.name = probe.name;
    result.source = instantiate(probe.baselineSource, p, q, k);
    result.variantSource = instantiate(probe.variantSource, p, q, k);
    result.expected = probe.expected;
    result.relation = probe.relation;
    result.requestedDigits = p;
    result.inputDigits = q;
    result.epsilonPower = k;
    return result;
}

[[nodiscard]] std::string budgetText(const evaluation::EvaluationUsage& usage) {
    std::ostringstream stream;
    stream << "work=" << usage.evaluationSteps
           << " depth=" << usage.maximumDepth
           << " nodes=" << usage.generatedNodes
           << " simplify=" << usage.simplificationCandidates
           << " refine=" << usage.certifiedRefinements
           << " precision-digits=" << usage.maximumRequestedPrecisionDigits;
    return stream.str();
}

[[nodiscard]] std::optional<std::string> invariantFailureReason(
    const GeneratedCase& testCase,
    const CaseObservation& observation,
    const std::optional<CaseObservation>& variantObservation) {
    if (observation.outcome != testCase.expected)
        return "baseline classification mismatch";

    if (testCase.kind == ProbeKind::Classification) {
        if (testCase.expected == BoundaryOutcome::PrecisionInsufficient
            && observation.usage.certifiedRefinements > kPrecisionRefinementCeiling)
            return "persistent information ambiguity consumed excessive certified refinements";
        return std::nullopt;
    }

    if (!variantObservation)
        return "metamorphic variant was not evaluated";
    const auto& variant = *variantObservation;

    if (variant.outcome != testCase.expected)
        return "metamorphic variant classification mismatch";

    switch (testCase.relation) {
    case MetamorphicRelation::SameValue:
        if (!observation.result || !variant.result)
            return "value relation did not produce two values";
        if (!(*observation.result == *variant.result))
            return "semantics-preserving transform changed the certified value or provenance";
        break;

    case MetamorphicRelation::InformationDoesNotNarrow: {
        if (!observation.result || !variant.result)
            return "information relation did not produce two values";
        const auto baselineBox = informationBox(*observation.result);
        const auto variantBox = informationBox(*variant.result);
        if (!baselineBox || !variantBox)
            return "information relation produced a non-approximate value";
        if (!contains(*variantBox, *baselineBox))
            return "outer N narrowed InformationEnclosure and invented information";
        break;
    }

    case MetamorphicRelation::InformationContainsExactRoute: {
        if (!observation.result || !variant.result)
            return "function provenance relation did not produce two values";
        const auto finiteInformation = informationBox(*observation.result);
        const auto exactCertified = certifiedBox(*variant.result);
        if (!finiteInformation || !exactCertified)
            return "function provenance relation produced a non-approximate value";
        if (!contains(*finiteInformation, *exactCertified))
            return "finite-input InformationEnclosure excluded the exact-route certified truth";
        break;
    }

    case MetamorphicRelation::StablePrecisionInsufficient:
        if (observation.usage.certifiedRefinements > kPrecisionRefinementCeiling
            || variant.usage.certifiedRefinements > kPrecisionRefinementCeiling)
            return "persistent ambiguity transform consumed excessive certified refinements";
        break;

    case MetamorphicRelation::SameOutcome:
        break;
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<Failure> checkCase(
    std::uint64_t masterSeed,
    std::uint64_t index,
    std::uint64_t timeoutMilliseconds) {
    GeneratedCase testCase = generateCase(masterSeed, index);
    CaseObservation observation = observe(testCase, timeoutMilliseconds);
    std::optional<CaseObservation> variantObservation;
    if (testCase.kind == ProbeKind::Metamorphic)
        variantObservation = observeSource(testCase.variantSource, timeoutMilliseconds);

    if (auto reason = invariantFailureReason(testCase, observation, variantObservation)) {
        return Failure{
            index,
            std::move(testCase),
            std::move(observation),
            std::move(variantObservation),
            std::move(*reason)};
    }
    return std::nullopt;
}

void printObservation(std::string_view prefix, const CaseObservation& observation) {
    std::cerr << prefix << " outcome : " << outcomeName(observation.outcome) << '\n';
    if (!observation.formattedResult.empty())
        std::cerr << prefix << " result  : " << observation.formattedResult << '\n';
    if (!observation.detail.empty())
        std::cerr << prefix << " detail  : " << observation.detail << '\n';
    std::cerr << prefix << " elapsed : " << observation.elapsedMilliseconds << " ms\n"
              << prefix << " budget  : " << budgetText(observation.usage) << '\n';
}

void printFailure(std::uint64_t masterSeed, const Failure& failure) {
    std::cerr
        << "\n[FAIL] Certification boundary invariant\n"
        << "Seed      : " << masterSeed << '\n'
        << "Case      : " << failure.index << '\n'
        << "Case seed : " << caseSeed(masterSeed, failure.index) << '\n'
        << "Family    : " << failure.testCase.family << '\n'
        << "Probe     : " << failure.testCase.name << '\n'
        << "Kind      : " << (failure.testCase.kind == ProbeKind::Metamorphic
            ? "Metamorphic" : "Classification") << '\n'
        << "Expression: " << failure.testCase.source << '\n';
    if (failure.testCase.kind == ProbeKind::Metamorphic) {
        std::cerr
            << "Variant   : " << failure.testCase.variantSource << '\n'
            << "Relation  : " << relationName(failure.testCase.relation) << '\n';
    }
    std::cerr
        << "Expected  : " << outcomeName(failure.testCase.expected) << '\n'
        << "Reason    : " << failure.reason << '\n';
    printObservation("Base", failure.observation);
    if (failure.variantObservation)
        printObservation("Variant", *failure.variantObservation);
    std::cerr
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
        << "Kind      : " << (testCase.kind == ProbeKind::Metamorphic
            ? "Metamorphic" : "Classification") << '\n'
        << "Expression: " << testCase.source << '\n';
    if (testCase.kind == ProbeKind::Metamorphic) {
        std::cout
            << "Variant   : " << testCase.variantSource << '\n'
            << "Relation  : " << relationName(testCase.relation) << '\n';
    }
    std::cout << "Expected  : " << outcomeName(testCase.expected) << '\n' << std::flush;

    CaseObservation observation = observe(testCase, timeoutMilliseconds);
    std::optional<CaseObservation> variantObservation;
    if (testCase.kind == ProbeKind::Metamorphic)
        variantObservation = observeSource(testCase.variantSource, timeoutMilliseconds);

    const auto printToStdout = [](std::string_view prefix, const CaseObservation& value) {
        std::cout << prefix << " outcome : " << outcomeName(value.outcome) << '\n';
        if (!value.formattedResult.empty())
            std::cout << prefix << " result  : " << value.formattedResult << '\n';
        if (!value.detail.empty())
            std::cout << prefix << " detail  : " << value.detail << '\n';
        std::cout << prefix << " elapsed : " << value.elapsedMilliseconds << " ms\n"
                  << prefix << " budget  : " << budgetText(value.usage) << '\n';
    };
    printToStdout("Base", observation);
    if (variantObservation)
        printToStdout("Variant", *variantObservation);

    const auto reason = invariantFailureReason(testCase, observation, variantObservation);
    if (!reason) {
        std::cout << "Status    : PASS\n";
        return true;
    }

    std::cout
        << "Status    : FAIL (" << *reason << ")\n"
        << "Reproduce : mmCal.Benchmarks --certification-boundaries --seed "
        << masterSeed << " --case " << index << '\n';
    return false;
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
        << "Classify  : " << probeTemplates().size() << " probes\n"
        << "Metamorph : " << metamorphicTemplates().size() << " probes\n"
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
