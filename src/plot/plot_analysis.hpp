#pragma once

#include "plot_request.hpp"

#include "evaluation/builtin_registry.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"

#include <optional>
#include <vector>

namespace mmcal::plot {

// 現在分かっているdomain情報が，要求区間内の実定義域をどこまで記述しているか。
enum class PlotDomainCoverage {
    Unknown,
    Partial,
    Complete
};

enum class PlotEndpointInclusion {
    Open,
    Closed
};

// Plotの数学座標上にある一つの連続区間。
// endpointはPi/2やRoot等のexact Exprをそのまま保持し，数値化はsampling直前まで遅延する。
struct PlotInterval final {
    expression::Expr lower;
    PlotEndpointInclusion lowerInclusion = PlotEndpointInclusion::Closed;
    expression::Expr upper;
    PlotEndpointInclusion upperInclusion = PlotEndpointInclusion::Closed;

    [[nodiscard]] bool operator==(const PlotInterval&) const = default;
};

// intervalsは「現時点でsampling対象として残っている区間」を表す。
// coverageがCompleteのときだけ，そのunionが要求区間内の実定義域全体であることを保証する。
struct PlotDomain final {
    PlotDomainCoverage coverage = PlotDomainCoverage::Unknown;
    std::vector<PlotInterval> intervals;

    [[nodiscard]] bool operator==(const PlotDomain&) const = default;
};

enum class PlotLandmarkKind {
    Pole,
    RemovableSingularity,
    JumpDiscontinuity,
    BranchPoint,
    DomainBoundary,
    UndefinedPoint,
    NonSmoothPoint,
    VerticalAsymptote
};

enum class PlotLandmarkConfidence {
    Proven,
    Candidate
};

// 解析上重要な一点に関する一つの事実。
// 同じpositionがPoleとVerticalAsymptoteの両方である場合は，別Landmarkとして併存させる。
struct PlotLandmark final {
    expression::Expr position;
    PlotLandmarkKind kind;
    PlotLandmarkConfidence confidence = PlotLandmarkConfidence::Proven;
    // removable singularity等で有限極限値までexactに証明できた場合だけ保持する。
    // 未定義点そのものの函数値ではなく，片側から共通に近づく描画上のsemantic値である。
    std::optional<expression::Expr> finiteLimit;

    [[nodiscard]] bool operator==(const PlotLandmark&) const = default;
};

enum class PlotSamplingSafety {
    Safe,
    UnsupportedDiscontinuity
};

struct PlotAnalysis final {
    PlotDomain domain;
    std::vector<PlotLandmark> landmarks;
    PlotSamplingSafety samplingSafety = PlotSamplingSafety::Safe;

    [[nodiscard]] bool operator==(const PlotAnalysis&) const = default;
};

// symbolic prepass前の解析結果。要求区間全体をsampling候補として保持するが，
// まだdomain completenessは主張しない。
[[nodiscard]] PlotAnalysis makeInitialPlotAnalysis(const PlotRequest& request);

// exact knowledgeだけを使うbounded symbolic prepass。
// 証明できない部分はUnknown/Candidateへ残し，sampling側が安全にfallbackできるようにする。
[[nodiscard]] PlotAnalysis analyzePlotRequest(
    const PlotRequest& request,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const expression::Symbol& infinitySymbol,
    const mathematics::AssumptionSet& assumptions = {});

} // namespace mmcal::plot
