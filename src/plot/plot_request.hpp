#pragma once

#include "expression/expr.hpp"
#include "expression/symbol.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace mmcal::plot {

struct PlotRangeOption final {
    expression::Expr minimum;
    expression::Expr maximum;

    [[nodiscard]] bool operator==(const PlotRangeOption&) const = default;
};

// PlotRange/AspectRatio/Ticks/PlotPointsは予約語ではなく通常のSymbolであり，PlotがRuleの左辺名を
// consumer側で解釈する。ここには評価済みoption値だけを保持する。
struct PlotRequestOptions final {
    std::optional<PlotRangeOption> plotRange;
    std::optional<expression::Expr> aspectRatio;
    std::optional<bool> ticks;
    // 100を現行既定密度とする相対sampling倍率。公開範囲は100..1024。
    std::optional<std::size_t> plotPoints;

    [[nodiscard]] bool operator==(const PlotRequestOptions&) const = default;
};

// Plot呼び出しから切り出した不変の意味情報。
// sampling/style詳細は後続層の責務とし，ここには数学的入力と公開view optionだけを保持する。
struct PlotRequest final {
    expression::Expr expression;
    expression::Symbol variable;
    expression::Expr lower;
    expression::Expr upper;
    PlotRequestOptions options{};

    [[nodiscard]] bool operator==(const PlotRequest&) const = default;
};

// 複数curve Plotの共通数学domain。各expressionは既存PlotRequestへ展開して
// analysis/compile/samplingできるため，単一curve pipelineを重複実装しない。
struct PlotRequestSet final {
    std::vector<expression::Expr> expressions;
    expression::Symbol variable;
    expression::Expr lower;
    expression::Expr upper;
    PlotRequestOptions options{};
};

[[nodiscard]] inline std::vector<PlotRequest> splitPlotRequests(const PlotRequestSet& request) {
    std::vector<PlotRequest> result;
    result.reserve(request.expressions.size());
    for (const auto& expression : request.expressions)
        result.push_back(PlotRequest{
            expression, request.variable, request.lower, request.upper, request.options});
    return result;
}

// ParametricPlotの1本の曲線。parameterはgeometry座標から独立したbinderとして保持する。
struct ParametricCurveRequest final {
    expression::Expr xExpression;
    expression::Expr yExpression;
    expression::Symbol parameter;
    expression::Expr lower;
    expression::Expr upper;
    PlotRequestOptions options{};

    [[nodiscard]] bool operator==(const ParametricCurveRequest&) const = default;
};

struct ParametricPlotRequestSet final {
    std::vector<std::pair<expression::Expr, expression::Expr>> curves;
    expression::Symbol parameter;
    expression::Expr lower;
    expression::Expr upper;
    PlotRequestOptions options{};
};

[[nodiscard]] inline std::vector<ParametricCurveRequest> splitParametricPlotRequests(
    const ParametricPlotRequestSet& request) {
    std::vector<ParametricCurveRequest> result;
    result.reserve(request.curves.size());
    for (const auto& [x, y] : request.curves)
        result.push_back(ParametricCurveRequest{
            x, y, request.parameter, request.lower, request.upper, request.options});
    return result;
}

} // namespace mmcal::plot
