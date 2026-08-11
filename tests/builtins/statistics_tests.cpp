// 記述統計の回帰テスト
#include "statistics_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

std::string eval(kernel::KernelSession& session, std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

error::CalcError evalError(kernel::KernelSession& session, std::string_view source) {
    try { static_cast<void>(session.evaluate(source)); }
    catch (const error::CalcError& e) { return e; }
    throw std::logic_error("Expected CalcError");
}

} // namespace

void runStatisticsTests(TestRunner& tests) {
    kernel::KernelSession session;

    tests.expectEqual(eval(session, "median[1,3,5]"), std::string{"3"},
        "odd median is exact");
    tests.expectEqual(eval(session, "median[{1,2,3,4}]"), std::string{"5/2"},
        "even median averages the middle pair exactly");
    tests.expectEqual(eval(session, "mode[1,2,2,3]"), std::string{"2"},
        "unique mode returns the modal value");
    tests.expectEqual(eval(session, "mode[1,1,2,2]"), std::string{"{1, 2}"},
        "tied modes are represented explicitly");
    tests.expectEqual(eval(session, "mode[1,2,3]"), std::string{"{}"},
        "a sample with no repeated value has no mode");

    tests.expectEqual(eval(session, "quantile[1/4,1,2,3,4,5,6,7]"), std::string{"5/2"},
        "quantile uses exact Type-7 linear interpolation");
    tests.expectEqual(eval(session, "percentile[50,{1,3,5}]"), std::string{"3"},
        "percentile maps 0..100 onto the quantile definition");
    tests.expectEqual(eval(session, "iqr[1,2,3,4]"), std::string{"3/2"},
        "IQR consistently uses Type-7 Q3-Q1");

    tests.expectEqual(eval(session, "var[1,2,3]"), std::string{"2/3"},
        "population variance remains rational");
    tests.expectEqual(eval(session, "vars[{1,2,3}]"), std::string{"1"},
        "sample variance accepts one vector");
    tests.expectEqual(eval(session, "stddev[1,2,3]"), std::string{"sqrt[6]/3"},
        "population standard deviation preserves a radical");
    tests.expectEqual(eval(session, "stddevs[{1,2,3}]"), std::string{"1"},
        "sample standard deviation is exact when possible");
    tests.expectEqual(eval(session, "geomean[1,4,1/32]"), std::string{"1/2"},
        "geometric mean uses exact roots");
    tests.expectEqual(eval(session, "harmmean[1,2,6]"), std::string{"9/5"},
        "harmonic mean remains rational");
    tests.expectEqual(eval(session, "rms[1,-1,1,-1]"), std::string{"1"},
        "RMS remains exact");

    tests.expectEqual(eval(session, "mad[1,1,2,2,4]"), std::string{"1"},
        "mad is median absolute deviation");
    tests.expectEqual(eval(session, "madR[1,2,3]"), std::string{"2/3"},
        "madR is mean absolute deviation about the arithmetic mean");
    tests.expectEqual(eval(session, "skew[1,2,3,4,5]"), std::string{"0"},
        "symmetric data has zero population moment skewness");
    tests.expectEqual(eval(session, "kurtp[-2,-1,0,1,2]"), std::string{"-13/10"},
        "kurtp is population excess kurtosis");
    tests.expectEqual(eval(session, "kurts[-2,-1,0,1,2]"), std::string{"-6/5"},
        "kurts uses unbiased Fisher sample excess kurtosis");

    tests.expectEqual(eval(session, "cv[10,10,10]"), std::string{"0"},
        "coefficient of variation is zero for a constant nonzero sample");
    tests.expectEqual(eval(session, "stderr[{1,2,3}]"), std::string{"sqrt[3]/3"},
        "standard error uses sample standard deviation divided by sqrt(n)");
    tests.expectEqual(eval(session, "zscore[5,3,1]"), std::string{"2"},
        "zscore is exact");

    tests.expectEqual(eval(session, "trimmean[1/5,1,2,100,3,4]"), std::string{"3"},
        "trimmean removes floor(p*n) observations from each tail");
    tests.expectEqual(eval(session, "winsor[1/5,1,2,100,3,4]"), std::string{"3"},
        "winsor returns the mean after symmetric winsorization");
    tests.expectEqual(eval(session, "winsorR[1/5,1,2,100,3,4]"), std::string{"{2, 2, 4, 3, 4}"},
        "winsorR returns clipped data in original order");

    tests.expectEqual(eval(session, "cov[{1,2,3},{2,4,6}]"), std::string{"4/3"},
        "population covariance is exact");
    tests.expectEqual(eval(session, "cov[1,2,3,2,4,6]"), std::string{"4/3"},
        "legacy even scalar covariance syntax splits observations in half");
    tests.expectEqual(eval(session, "corr[{1,2,3},{2,4,6}]"), std::string{"1"},
        "Pearson correlation simplifies exactly");
    tests.expectEqual(eval(session, "corrspearman[{1,2,3},{10,20,30}]"), std::string{"1"},
        "Spearman correlation uses average ranks");
    tests.expectEqual(eval(session, "corrspearman[{1,1,2},{10,10,20}]"), std::string{"1"},
        "Spearman ranking handles ties exactly");
    tests.expectEqual(eval(session, "percentrank[3,1,2,3,4,5]"), std::string{"1/2"},
        "percentrank is exact at an observed value");
    tests.expectEqual(eval(session, "percentrank[5/2,1,2,3,4]"), std::string{"1/2"},
        "percentrank linearly interpolates between observations");
    tests.expectEqual(eval(session, "percentrank[2,1,2,2,3]"), std::string{"1/2"},
        "percentrank averages tied observed ranks");

    tests.expect(evalError(session, "quantile[-1/10,1,2,3]").type() == error::CalcErrorType::Domain,
        "quantile rejects p outside [0,1]");
    tests.expect(evalError(session, "geomean[-1,2]").type() == error::CalcErrorType::Domain,
        "real geometric mean rejects negative observations");
    tests.expect(evalError(session, "corr[{1,1},{2,3}]").type() == error::CalcErrorType::Domain,
        "correlation rejects zero variance");
    tests.expect(evalError(session, "zscore[1,0,0]").type() == error::CalcErrorType::Domain,
        "zscore requires a positive sigma");
}

} // namespace mmcal::tests
