// ユーザー変数環境の回帰テスト
#include "environment_tests.hpp"

#include "evaluation/environment.hpp"
#include "formatting/expr_formatter.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "test_framework.hpp"

namespace mmcal::tests {

void runEnvironmentTests(TestRunner& tests) {
    using evaluation::Environment;
    using expression::Expr;
    using expression::Symbol;
    using formatting::formatExpr;
    using numeric::BigInt;
    using numeric::Number;

    Environment environment;
    const Symbol x{"x"};
    const Symbol y{"y"};

    tests.expectEqual(environment.size(), std::size_t{0}, "Environment: starts empty");
    tests.expect(!environment.contains(x), "Environment: missing symbol is not contained");
    tests.expect(environment.find(x) == nullptr, "Environment: missing symbol has no binding");

    environment.set(x, Expr{Number{BigInt{3}}});
    tests.expect(environment.contains(x), "Environment: contains assigned symbol");
    tests.expectEqual(environment.size(), std::size_t{1}, "Environment: counts assigned symbol");
    tests.expectEqual(formatExpr(*environment.find(x)), "3", "Environment: retrieves assigned value");

    environment.set(x, Expr{Number{BigInt{5}}});
    tests.expectEqual(environment.size(), std::size_t{1}, "Environment: assignment replaces existing value");
    tests.expectEqual(formatExpr(*environment.find(x)), "5", "Environment: retrieves replacement value");

    environment.set(y, Expr{Symbol{"x"}});
    tests.expectEqual(environment.size(), std::size_t{2}, "Environment: stores multiple symbols");
    tests.expect(environment.erase(x), "Environment: erases existing symbol");
    tests.expect(!environment.erase(x), "Environment: reports missing erase target");


    environment.set(x, Expr{Number{BigInt{10}}});
    environment.pushScope();
    environment.setLocal(x, Expr{Number{BigInt{20}}});
    tests.expectEqual(formatExpr(*environment.find(x)), "20",
        "Environment: local binding shadows global value");
    tests.expect(environment.containsLocal(x),
        "Environment: reports active local binding");
    environment.assign(x, Expr{Number{BigInt{30}}});
    tests.expectEqual(formatExpr(*environment.find(x)), "30",
        "Environment: assign updates nearest local binding");
    environment.popScope();
    tests.expectEqual(formatExpr(*environment.find(x)), "10",
        "Environment: leaving scope restores global value");
    tests.expectEqual(environment.localDepth(), std::size_t{0},
        "Environment: local scope depth returns to zero");

    environment.clear();
    tests.expectEqual(environment.size(), std::size_t{0}, "Environment: clear removes all bindings");
}

} // namespace mmcal::tests
