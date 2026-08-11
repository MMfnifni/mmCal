// 非再帰タスクスタック式評価器の回帰テスト
#include "evaluator_tests.hpp"

#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "evaluation/environment.hpp"
#include "evaluation/evaluator.hpp"
#include "formatting/expr_formatter.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "syntax/lexer.hpp"
#include "syntax/lowerer.hpp"
#include "syntax/parser.hpp"
#include "test_framework.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mmcal::tests {
namespace {

using expression::Expr;
using expression::Symbol;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr rational(std::int64_t numerator, std::int64_t denominator) {
    return Expr{Number{Rational{BigInt{numerator}, BigInt{denominator}}}};
}

[[nodiscard]] Expr call(std::string_view head, std::vector<Expr> arguments) {
    return Expr::call(Symbol{head}, std::move(arguments));
}

[[nodiscard]] std::string evaluateText(
    std::string text,
    evaluation::Environment& environment) {
    auto sourceText = std::make_shared<const std::string>(std::move(text));
    syntax::Lexer lexer{*sourceText};
    syntax::Parser parser{sourceText, lexer.tokenize()};
    const syntax::SyntaxTree tree = parser.parse();
    syntax::Lowerer lowerer;
    syntax::LoweringResult lowered = lowerer.lowerTracked(tree);
    evaluation::Evaluator evaluator{environment};
    return formatting::formatExpr(evaluator.evaluate(lowered.expression, lowered.origins));
}

[[nodiscard]] error::CalcError evaluateTextError(
    std::string text,
    evaluation::Environment& environment) {
    auto sourceText = std::make_shared<const std::string>(std::move(text));

    try {
        syntax::Lexer lexer{*sourceText};
        syntax::Parser parser{sourceText, lexer.tokenize()};
        const syntax::SyntaxTree tree = parser.parse();
        syntax::Lowerer lowerer;
        syntax::LoweringResult lowered = lowerer.lowerTracked(tree);
        evaluation::Evaluator evaluator{environment};
        static_cast<void>(evaluator.evaluate(lowered.expression, lowered.origins));
    }
    catch (const error::CalcError& exception) {
        return exception;
    }

    throw std::logic_error("Expected CalcError was not thrown");
}

} // namespace

void runEvaluatorTests(TestRunner& tests) {
    using builtins::names::add;
    using builtins::names::divide;
    using builtins::names::factorial;
    using builtins::names::multiply;
    using builtins::names::negate;
    using builtins::names::power;
    using builtins::names::set;
    using builtins::names::sqrt;
    using builtins::names::subtract;
    using evaluation::Environment;
    using evaluation::Evaluator;
    using formatting::formatExpr;

    Environment environment;
    Evaluator evaluator{environment};

    tests.expectEqual(formatExpr(evaluator.evaluate(integer(42))), "42", "Evaluator: leaves numbers unchanged");
    tests.expectEqual(formatExpr(evaluator.evaluate(Expr{true})), "True", "Evaluator: leaves booleans unchanged");
    tests.expectEqual(formatExpr(evaluator.evaluate(Expr{Symbol{"Pi"}})), "Pi", "Evaluator: preserves known symbolic constants");
    tests.expectEqual(formatExpr(evaluator.evaluate(Expr{Symbol{"unknown"}})), "unknown",
        "Evaluator: temporarily preserves unbound symbols for symbolic evaluation");

    environment.set(Symbol{"x"}, integer(3));
    environment.set(Symbol{"y"}, Expr{Symbol{"x"}});
    environment.set(Symbol{"q"}, Expr{Symbol{"Pi"}});
    environment.set(Symbol{"a"}, Expr{Symbol{"Pi"}});
    environment.set(Symbol{"b"}, Expr{Symbol{"E"}});
    environment.set(Symbol{"c"}, integer(5));
    environment.set(Symbol{"d"}, integer(2));
    tests.expectEqual(formatExpr(evaluator.evaluate(Expr{Symbol{"y"}})), "3", "Evaluator: resolves chained symbols");

    const Expr assigned = evaluator.evaluate(call(set, {Expr{Symbol{"z"}}, rational(1, 3)}));
    tests.expectEqual(formatExpr(assigned), "1/3", "Evaluator: Set returns assigned value");
    tests.expectEqual(formatExpr(evaluator.evaluate(Expr{Symbol{"z"}})), "1/3", "Evaluator: Set stores assigned value");

    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(call(set, {integer(1), integer(2)})));
    }, "Evaluator: Set requires symbol first argument");

    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(call(set, {Expr{Symbol{"Pi"}}, integer(3)})));
    }, "Evaluator: Set rejects protected predefined symbols");

    tests.expectEqual(
        formatExpr(evaluator.evaluate(call(add, {integer(2), integer(3), rational(1, 2)}))),
        "11/2",
        "Evaluator: adds exact numeric arguments");
    tests.expectEqual(
        formatExpr(evaluator.evaluate(call(add, {Expr{Symbol{"q"}}, integer(0)}))),
        "Pi",
        "Evaluator: removes additive zero");
    tests.expectEqual(
        formatExpr(evaluator.evaluate(call(add, {
            call(add, {Expr{Symbol{"a"}}, integer(2)}),
            integer(3)
        }))),
        "5+Pi",
        "Evaluator: flattens nested addition and combines numbers");

    Environment simplificationEnvironment;
    tests.expectEqual(
        evaluateText("exp[log[Pi]]", simplificationEnvironment),
        std::string{"Pi"},
        "Evaluator: shared Simplifier combines Exp and principal Log safely");
    tests.expectEqual(
        evaluateText("log[exp[2]]", simplificationEnvironment),
        std::string{"2"},
        "Evaluator: shared Simplifier uses real-domain knowledge for Log[Exp[x]]");

    tests.expectEqual(
        formatExpr(evaluator.evaluate(call(subtract, {integer(7), rational(1, 2)}))),
        "13/2",
        "Evaluator: subtracts exact numbers");
    tests.expectEqual(
        formatExpr(evaluator.evaluate(call(subtract, {Expr{Symbol{"q"}}, Expr{Symbol{"q"}}}))),
        "0",
        "Evaluator: simplifies equal-expression subtraction");

    tests.expectEqual(
        formatExpr(evaluator.evaluate(call(multiply, {integer(2), rational(3, 4)}))),
        "3/2",
        "Evaluator: multiplies exact numbers");
    tests.expectEqual(
        formatExpr(evaluator.evaluate(call(multiply, {Expr{Symbol{"q"}}, integer(1)}))),
        "Pi",
        "Evaluator: removes multiplicative one");
    tests.expectEqual(
        formatExpr(evaluator.evaluate(call(multiply, {Expr{Symbol{"q"}}, integer(0)}))),
        "0",
        "Evaluator: short-circuits multiplication by zero");

    tests.expectEqual(
        formatExpr(evaluator.evaluate(call(divide, {integer(7), integer(2)}))),
        "7/2",
        "Evaluator: keeps division exact");
    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(call(divide, {integer(1), integer(0)})));
    }, "Evaluator: rejects division by zero");

    tests.expectEqual(formatExpr(evaluator.evaluate(call(power, {integer(2), integer(100)}))),
        "1267650600228229401496703205376", "Evaluator: exact positive integer power");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(power, {integer(2), integer(-3)}))),
        "1/8", "Evaluator: exact negative integer power");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(power, {rational(2, 3), integer(4)}))),
        "16/81", "Evaluator: exact rational power");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(power, {
        Expr{Number::complex(numeric::RealNumber{}, numeric::RealNumber{BigInt{1}})},
        integer(4)
    }))), "1", "Evaluator: exact imaginary-unit power");
    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(call(power, {integer(0), integer(0)})));
    }, "Evaluator: treats zero power zero as indeterminate");
    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(call(power, {integer(0), integer(-1)})));
    }, "Evaluator: rejects negative powers of zero");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(power, {integer(2), rational(1, 2)}))),
        "sqrt[2]", "Evaluator: half power uses principal square root semantics");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(power, {integer(-2), rational(1, 2)}))),
        "I sqrt[2]", "Evaluator: negative half power promotes to principal complex root");

    tests.expectEqual(formatExpr(evaluator.evaluate(call(factorial, {integer(0)}))), "1", "Evaluator: zero factorial");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(factorial, {integer(20)}))),
        "2432902008176640000", "Evaluator: exact factorial");
    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(call(factorial, {integer(-1)})));
    }, "Evaluator: rejects negative factorial");
    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(call(factorial, {rational(1, 2)})));
    }, "Evaluator: factorial requires integer");

    tests.expectEqual(formatExpr(evaluator.evaluate(call(negate, {integer(5)}))), "-5", "Evaluator: negates number");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(sqrt, {integer(144)}))), "12", "Evaluator: exact integer sqrt");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(sqrt, {rational(9, 16)}))), "3/4", "Evaluator: exact rational sqrt");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(sqrt, {integer(-9)}))), "3I", "Evaluator: exact negative sqrt");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(sqrt, {integer(2)}))), "sqrt[2]", "Evaluator: preserves irrational sqrt");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(sqrt, {
        Expr{Number::complex(numeric::RealNumber{BigInt{3}}, numeric::RealNumber{BigInt{4}})}
    }))), "2+I", "Evaluator: exact principal sqrt closes over rational complex components");
    tests.expectEqual(formatExpr(evaluator.evaluate(call(sqrt, {
        Expr{Number::complex(numeric::RealNumber{BigInt{-3}}, numeric::RealNumber{BigInt{-4}})}
    }))), "1-2I", "Evaluator: exact principal sqrt chooses lower-half-plane imaginary sign");

    tests.expectEqual(evaluateText("1 < 2 <= 2", environment), "True", "Evaluator: exact comparison chain");
    tests.expectEqual(evaluateText("2 == 2", environment), "True", "Evaluator: exact equality");
    tests.expectEqual(evaluateText("2 != 2", environment), "False", "Evaluator: exact inequality");
    tests.expectEqual(evaluateText("Pi ^ 2 == Pi * Pi", environment),
        "True", "Evaluator: structural multiplication normalization can prove exact equality");

    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(call("unknownFunction", {integer(1)})));
    }, "Evaluator: rejects unknown functions");

    environment.set(Symbol{"cycleA"}, Expr{Symbol{"cycleB"}});
    environment.set(Symbol{"cycleB"}, Expr{Symbol{"cycleA"}});
    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(Expr{Symbol{"cycleA"}}));
    }, "Evaluator: detects cyclic symbol definitions");

    const error::CalcError divisionError = evaluateTextError("1 + 8 / (3 - 3)", environment);
    tests.expect(divisionError.type() == error::CalcErrorType::Domain,
        "Evaluator: division error keeps domain type");
    tests.expect(divisionError.span().has_value(),
        "Evaluator: division error has source span");
    if (divisionError.span()) {
        tests.expectEqual(divisionError.span()->begin.column, std::size_t{10},
            "Evaluator: division error points to denominator");
        tests.expectEqual(divisionError.span()->end.column, std::size_t{15},
            "Evaluator: division error covers denominator expression");
    }

    tests.expectEqual(evaluateText("2+missing", environment), "2+missing",
        "Evaluator: unbound symbols are temporarily preserved as free symbols");

    tests.expectEqual(formatExpr(call("f", {integer(2)})), "f[2]", "Formatter: uses brackets for function calls");
    tests.expectEqual(formatExpr(call(sqrt, {integer(2)})), "sqrt[2]", "Formatter: uses brackets for builtin functions");

    evaluator.setDepthLimit(2);
    tests.expectThrows<error::CalcError>([&] {
        static_cast<void>(evaluator.evaluate(call(add, {call(add, {integer(1), integer(2)}), integer(3)})));
    }, "Evaluator: enforces depth limit");
    evaluator.setDepthLimit(1024);

    tests.expectThrows<std::invalid_argument>([&] {
        evaluator.setDepthLimit(0);
    }, "Evaluator: rejects zero depth limit");
}

} // namespace mmcal::tests
