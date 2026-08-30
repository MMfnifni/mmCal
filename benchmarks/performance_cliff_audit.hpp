#pragma once

#include <cstddef>

namespace mmcal::benchmarks {

struct PerformanceCliffAuditOptions final {
    std::size_t iterations = 1;
    double cliffRatio = 4.0;
};

// 公開KernelSession経由で特殊函数のargument/precision sweepを行い，wall timeと
// EvaluationBudget telemetryの隣接比からbounded-work境界・performance cliffを抽出する。
void runSpecialFunctionPerformanceCliffAudit(
    const PerformanceCliffAuditOptions& options = {});

// Complex Algebraic Rootの初期分離と単根refinementを分けて測定する。
// 対称疎多項式・generic dense・近接根を同じrunnerで追跡する。
void runAlgebraicRootPerformanceCliffAudit(
    const PerformanceCliffAuditOptions& options = {});

} // namespace mmcal::benchmarks
