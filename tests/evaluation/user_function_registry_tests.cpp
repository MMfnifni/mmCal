// ユーザー定義函数の登録の回帰テスト
#include "user_function_registry_tests.hpp"

#include "evaluation/user_function_registry.hpp"
#include "formatting/expr_formatter.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "test_framework.hpp"

#include <string>
#include <vector>

namespace mmcal::tests {
namespace {

[[nodiscard]] expression::Expr integer(std::int64_t value) {
    return expression::Expr{numeric::Number{numeric::BigInt{value}}};
}

} // namespace

void runUserFunctionRegistryTests(TestRunner& tests) {
    using evaluation::UserFunctionDefinition;
    using evaluation::UserFunctionRegistry;
    using expression::Symbol;

    UserFunctionRegistry registry;
    const Symbol f{"f"};
    tests.expectEqual(registry.size(), std::size_t{0},
        "UserFunctionRegistry: starts empty");

    registry.define(UserFunctionDefinition{
        Symbol{"f"},
        {Symbol{"x"}},
        integer(1)
    });

    tests.expect(registry.contains(f),
        "UserFunctionRegistry: contains defined name");
    tests.expect(registry.find(f, 1) != nullptr,
        "UserFunctionRegistry: finds matching arity");
    tests.expect(registry.find(f, 2) == nullptr,
        "UserFunctionRegistry: separates arities");

    registry.define(UserFunctionDefinition{
        Symbol{"f"},
        {Symbol{"x"}, Symbol{"y"}},
        integer(2)
    });
    tests.expectEqual(registry.size(), std::size_t{2},
        "UserFunctionRegistry: stores overloads by arity");
    tests.expect(registry.arities(f) == std::vector<std::size_t>({1, 2}),
        "UserFunctionRegistry: reports available arities");

    registry.define(UserFunctionDefinition{
        Symbol{"f"},
        {Symbol{"value"}},
        integer(9)
    });
    tests.expectEqual(registry.size(), std::size_t{2},
        "UserFunctionRegistry: replaces same-arity definition");
    tests.expectEqual(
        formatting::formatExpr(registry.find(f, 1)->body),
        std::string{"9"},
        "UserFunctionRegistry: replacement becomes active");

    tests.expect(registry.erase(f),
        "UserFunctionRegistry: erases all overloads by name");
    tests.expect(!registry.contains(f),
        "UserFunctionRegistry: erased name is absent");

    registry.define(UserFunctionDefinition{Symbol{"g"}, {}, integer(3)});
    registry.clear();
    tests.expectEqual(registry.size(), std::size_t{0},
        "UserFunctionRegistry: clear removes all definitions");
}

} // namespace mmcal::tests
