// 構文解析の回帰テスト
#include "parser_tests.hpp"

#include "error/error_message.hpp"
#include "syntax/lexer.hpp"
#include "syntax/parser.hpp"
#include "test_framework.hpp"

#include <memory>
#include <string>

namespace mmcal::tests {
namespace {

[[nodiscard]] syntax::SyntaxTree parse(std::string text) {
    auto sourceText = std::make_shared<const std::string>(std::move(text));
    syntax::Lexer lexer{*sourceText};
    syntax::Parser parser{sourceText, lexer.tokenize()};
    return parser.parse();
}

} // namespace

void runParserTests(TestRunner& tests) {
    {
        const auto tree = parse("1 + 2 * 3");
        const auto* add = std::get_if<syntax::BinarySyntax>(&tree.root().data);
        tests.expect(add && add->operation == syntax::BinaryOperator::Add,
            "parser preserves additive precedence");
        const auto* multiply = add
            ? std::get_if<syntax::BinarySyntax>(&add->right->data)
            : nullptr;
        tests.expect(multiply && multiply->operation == syntax::BinaryOperator::Multiply,
            "parser nests multiplication under addition");
    }

    {
        const auto tree = parse("2^3^4");
        const auto* outer = std::get_if<syntax::BinarySyntax>(&tree.root().data);
        const auto* inner = outer
            ? std::get_if<syntax::BinarySyntax>(&outer->right->data)
            : nullptr;
        tests.expect(
            outer && outer->operation == syntax::BinaryOperator::Power
                && inner && inner->operation == syntax::BinaryOperator::Power,
            "parser makes power right associative");
    }

    {
        const auto tree = parse("-2^2");
        const auto* unary = std::get_if<syntax::UnarySyntax>(&tree.root().data);
        const auto* power = unary
            ? std::get_if<syntax::BinarySyntax>(&unary->operand->data)
            : nullptr;
        tests.expect(
            unary && unary->operation == syntax::UnaryOperator::Minus
                && power && power->operation == syntax::BinaryOperator::Power,
            "parser applies power before unary minus");
    }

    {
        const auto tree = parse("sqrt[2]");
        const auto* call = std::get_if<syntax::CallSyntax>(&tree.root().data);
        tests.expect(
            call && call->name == "sqrt" && call->arguments.size() == 1,
            "parser recognizes square-bracket function calls");
    }

    {
        const auto tree = parse("f[x, y] := x + y");
        const auto* assignment = std::get_if<syntax::AssignmentSyntax>(&tree.root().data);
        const auto* signature = assignment
            ? std::get_if<syntax::FunctionSignatureSyntax>(&assignment->target->data)
            : nullptr;
        tests.expect(
            signature && signature->name == "f" && signature->parameters.size() == 2,
            "parser recognizes square-bracket function signature assignment");
    }

    {
        const auto tree = parse("x(x + 1)");
        const auto* multiply = std::get_if<syntax::BinarySyntax>(&tree.root().data);
        tests.expect(
            multiply && multiply->operation == syntax::BinaryOperator::ImplicitMultiply,
            "parser treats identifier followed by a group as implicit multiplication");
    }

    {
        const auto tree = parse("@");
        const auto* history = std::get_if<syntax::HistoryReferenceSyntax>(&tree.root().data);
        tests.expect(
            history && history->kind == syntax::HistoryReferenceKind::Input
                && history->depth == 1,
            "parser recognizes previous-input shorthand");
    }

    {
        const auto tree = parse("{1, {2, 3}}");
        tests.expect(
            std::holds_alternative<syntax::ArrayLiteralSyntax>(tree.root().data),
            "parser stores nested array syntax");
        tests.expect(tree.root().span.begin.column == 1 && tree.root().span.end.column == 12,
            "parser stores source span on root");
    }

    {
        const auto tree = parse("@@@");
        const auto* history = std::get_if<syntax::HistoryReferenceSyntax>(&tree.root().data);
        tests.expect(
            history && history->kind == syntax::HistoryReferenceKind::Input
                && history->depth == 3,
            "parser recognizes multi-depth input-history shorthand");
    }
    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(parse("[1 + 2]")); },
        "parser rejects square bracket grouping");
    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(parse("Pi := 3")); },
        "parser rejects constant assignment");
    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(parse("f[x, x] := x")); },
        "parser rejects duplicate function parameters");
    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(parse("f(x) := x")); },
        "parser rejects parenthesized function definition syntax");
    tests.expect(
        std::holds_alternative<syntax::UnitAppliedSyntax>(parse("(30)rad").root().data),
        "parser accepts an angle-unit suffix on a grouped expression");
    tests.expect(
        std::holds_alternative<syntax::UnitAppliedSyntax>(parse("Pi/6 Rad").root().data),
        "parser applies a unit suffix to the complete preceding term");
}

} // namespace mmcal::tests
