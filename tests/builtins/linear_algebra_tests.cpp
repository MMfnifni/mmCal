// 行列・線形代数の回帰テスト
#include "linear_algebra_tests.hpp"

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

void runLinearAlgebraTests(TestRunner& tests) {
    kernel::KernelSession session;
    tests.expectEqual(eval(session, "transpose[{{1,2,3},{4,5,6}}]"),
        std::string{"{{1, 4}, {2, 5}, {3, 6}}"},
        "Linear algebra: transpose preserves exact elements");
    tests.expectEqual(eval(session, "madd[{{1,2},{3,4}},{{5,6},{7,8}}]"),
        std::string{"{{6, 8}, {10, 12}}"},
        "Linear algebra: matrix addition is exact");
    tests.expectEqual(eval(session, "matmul[{{1,2},{3,4}},{{5,6},{7,8}}]"),
        std::string{"{{19, 22}, {43, 50}}"},
        "Linear algebra: matrix multiplication is exact");
    tests.expectEqual(eval(session, "matmul[{1,2,3},{4,5,6}]"), std::string{"32"},
        "Linear algebra: vector dot product uses matmul");
    tests.expectEqual(eval(session, "det[{{1,2},{3,4}}]"), std::string{"-2"},
        "Linear algebra: determinant is exact");
    tests.expectEqual(eval(session, "det[{{a,b},{c,d}}]"), std::string{"a d-b c"},
        "Linear algebra: symbolic determinant remains exact");
    tests.expectEqual(eval(session, "inverse[{{1,2},{3,4}}]"),
        std::string{"{{-2, 1}, {3/2, -1/2}}"},
        "Linear algebra: exact rational matrix inverse");
    tests.expectEqual(eval(session, "inverse[{{a,b},{c,d}}]"),
        std::string{"{{d/(a d-b c), -b/(a d-b c)}, {-c/(a d-b c), a/(a d-b c)}}"},
        "Linear algebra: symbolic inverse keeps determinant denominators");
    tests.expectEqual(eval(session, "rref[{{1,2},{3,4}}]"),
        std::string{"{{1, 0}, {0, 1}}"},
        "Linear algebra: rref uses exact pivots");
    tests.expectEqual(eval(session, "rank[{{1,2},{2,4}}]"), std::string{"1"},
        "Linear algebra: rank is exact for rational matrices");
    tests.expectEqual(eval(session, "rank[{{x,0},{0,1}}]"),
        std::string{"rank[{{x, 0}, {0, 1}}]"},
        "Linear algebra: symbolic rank does not guess an undecidable pivot");
    const auto singular = evalError(session, "inverse[{{1,2},{2,4}}]");
    tests.expect(singular.type() == error::CalcErrorType::Domain,
        "Linear algebra: singular inverse is a domain error");
    const auto mismatch = evalError(session, "matmul[{{1,2}},{{1,2}}]");
    tests.expect(mismatch.type() == error::CalcErrorType::Domain,
        "Linear algebra: incompatible dimensions are rejected");
}

} // namespace mmcal::tests
