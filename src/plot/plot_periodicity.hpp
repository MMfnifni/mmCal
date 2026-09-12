#pragma once

#include "plot_request.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/rational.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace mmcal::plot {

// 有限trigonometric polynomialのparameter周波数support。
// frequency=0は定数成分を表す。係数を保持しないためcancellation後の最小supportは
// 主張しないが，ここに現れる全周波数の共通周期は元式の厳密な周期として使える。
struct HarmonicSpectrum final {
    std::vector<numeric::Rational> frequencies;
};

struct HarmonicSpectrumOptions final {
    std::size_t maximumFrequencies = 256;
    std::size_t maximumPower = 64;
};

enum class PeriodProofKind {
    Registry,
    Composition,
    HarmonicSpectrum
};

// 「periodであること」と「最小periodであること」を分離するcertificate。
// 現段階のHarmonicSpectrumはsupportを保守的に過大評価し得るため，通常は
// provenMinimal=falseのままでもsampling短縮には安全に利用できる。
struct PeriodCertificate final {
    numeric::Rational turns;
    bool provenMinimal = false;
    PeriodProofKind proofKind = PeriodProofKind::Composition;
    std::optional<HarmonicSpectrum> spectrum;
};

struct ParametricPeriodReduction final {
    expression::Expr originalLower;
    expression::Expr originalUpper;
    expression::Expr effectiveUpper;
    expression::Expr period;
    std::size_t repetitions = 1;
    PeriodCertificate certificate;
};

// sin/cosのaffine phaseと四則・非負整数冪から有限周波数supportを厳密に構成する。
// 解析不能はnulloptであり，非周期の断定ではない。
[[nodiscard]] std::optional<HarmonicSpectrum> analyzeHarmonicSpectrum(
    const expression::Expr& expression,
    const expression::Symbol& parameter,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const HarmonicSpectrumOptions& options = {});

// x/y双方から証明できる共通周期をcertificateとして返す。
[[nodiscard]] std::optional<PeriodCertificate> analyzeParametricCurvePeriod(
    const ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const HarmonicSpectrumOptions& harmonicOptions = {});

[[nodiscard]] expression::Expr parametricPeriodExpression(
    const numeric::Rational& turns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] std::optional<std::size_t> exactParametricPeriodRepetitions(
    const ParametricCurveRequest& request,
    const numeric::Rational& periodTurns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {});

// 元spanがperiodのexact整数倍(>=2)のときだけupperをlower+periodへ短縮する。
[[nodiscard]] std::optional<ParametricPeriodReduction> reduceParametricCurvePeriod(
    ParametricCurveRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions = {},
    const HarmonicSpectrumOptions& harmonicOptions = {});

} // namespace mmcal::plot
