// 構文木からExprへのloweringの回帰テスト
#include "lowerer_tests.hpp"

#include "error/error_message.hpp"
#include "evaluation/environment.hpp"
#include "evaluation/evaluator.hpp"
#include "formatting/expr_formatter.hpp"
#include "syntax/lexer.hpp"
#include "syntax/lowerer.hpp"
#include "syntax/parser.hpp"
#include "test_framework.hpp"

#include <memory>
#include <string>

namespace mmcal::tests {
namespace {

[[nodiscard]] expression::Expr lower(
    std::string text,
    syntax::LoweringOptions options = syntax::LoweringOptions::defaults()) {
    auto sourceText = std::make_shared<const std::string>(std::move(text));
    syntax::Lexer lexer{*sourceText};
    syntax::Parser parser{sourceText, lexer.tokenize()};
    const syntax::SyntaxTree tree = parser.parse();
    return syntax::Lowerer{std::move(options)}.lower(tree);
}

[[nodiscard]] std::string lowerAndFormat(std::string text) {
    return formatting::formatExpr(lower(std::move(text)));
}

[[nodiscard]] std::string evaluateAndFormat(std::string text) {
    evaluation::Environment environment;
    evaluation::Evaluator evaluator{environment};
    return formatting::formatExpr(evaluator.evaluate(lower(std::move(text))));
}

} // namespace

void runLowererTests(TestRunner& tests) {
    tests.expectEqual(evaluateAndFormat("1 + 2 * 3"), std::string{"7"},
        "frontend evaluates arithmetic precedence");
    tests.expectEqual(evaluateAndFormat("0.1 + 0.2"), std::string{"3/10"},
        "frontend keeps decimal arithmetic exact");
    tests.expectEqual(evaluateAndFormat("2#101.01"), std::string{"21/4"},
        "frontend lowers general radix fraction");
    tests.expectEqual(evaluateAndFormat("0xFF + 1"), std::string{"256"},
        "frontend lowers prefixed integer");
    tests.expectEqual(evaluateAndFormat("sqrt[9/16]"), std::string{"3/4"},
        "frontend supports bracket function call");
    tests.expectEqual(evaluateAndFormat("I * I"), std::string{"-1"},
        "uppercase I is the imaginary unit");
    tests.expectEqual(lowerAndFormat("i"), std::string{"i"},
        "lowercase i remains a symbol");
    tests.expectEqual(lowerAndFormat("2Pi"), std::string{"2Pi"},
        "frontend lowers implicit multiplication");
    tests.expectEqual(lowerAndFormat("Pi(2)"), std::string{"Pi*2"},
        "constant followed by parentheses is multiplication");
    tests.expectEqual(lowerAndFormat("f[2]"), std::string{"f[2]"},
        "square bracket syntax lowers to a function call");
    tests.expectEqual(lowerAndFormat("f[x]:=x+1"), std::string{"f[x]:=x+1"},
        "function definition syntax is preserved for later evaluation");
    tests.expectEqual(lowerAndFormat("1 < x <= 3"), std::string{"And[1<x, x<=3]"},
        "comparison chain lowers without losing operands");
    tests.expectEqual(lowerAndFormat("%%"), std::string{"%%"},
        "history reference depth is preserved");
    tests.expectEqual(lowerAndFormat("30deg"), std::string{"30 deg"},
        "unit suffix is represented explicitly");
    tests.expectEqual(lowerAndFormat("\"a\\nb\""), std::string{"\"a\\nb\""},
        "string literal round trip");
    tests.expectEqual(lowerAndFormat("{{1, 2}, {3, 4}}"), std::string{"{{1, 2}, {3, 4}}"},
        "rectangular array shape is inferred");

    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(lower("{1, {2, 3}}")); },
        "lowerer rejects ragged array");
    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(lower("2#102")); },
        "lowerer rejects a digit outside the radix");
}

} // namespace mmcal::tests
