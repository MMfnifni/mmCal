// 構文解析の回帰テスト
#include "parser_tests.hpp"

#include "error/error_message.hpp"
#include "syntax/lexer.hpp"
#include "syntax/parser.hpp"
#include "test_framework.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

[[nodiscard]] syntax::SyntaxTree parse(std::string text) {
    auto sourceText = std::make_shared<const std::string>(std::move(text));
    syntax::Lexer lexer{*sourceText};
    syntax::Parser parser{sourceText, lexer.tokenize()};
    return parser.parse();
}

[[nodiscard]] syntax::SyntaxTree parseWithLimits(
    std::string text,
    syntax::ParseLimits limits) {
    auto sourceText = std::make_shared<const std::string>(std::move(text));
    syntax::ParseBudget budget{limits};
    syntax::Lexer lexer{*sourceText, &budget};
    syntax::Parser parser{
        sourceText,
        lexer.tokenize(),
        syntax::ParserOptions::defaults(),
        &budget};
    return parser.parse();
}

[[nodiscard]] error::CalcError parseError(
    std::string text,
    syntax::ParseLimits limits = {}) {
    try {
        static_cast<void>(parseWithLimits(std::move(text), limits));
    }
    catch (const error::CalcError& exception) {
        return exception;
    }
    throw std::logic_error("Expected CalcError was not thrown");
}

[[nodiscard]] std::string nestedInput(
    std::string_view opening,
    std::string_view closing,
    std::size_t depth) {
    std::string result;
    result.reserve((opening.size() + closing.size()) * depth + 1);
    for (std::size_t index = 0; index < depth; ++index)
        result += opening;
    result += '1';
    for (std::size_t index = 0; index < depth; ++index)
        result += closing;
    return result;
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

    {
        syntax::ParseLimits limits;
        limits.maxOperatorChain = 2;
        static_cast<void>(parseWithLimits("1+2+3", limits));
        const error::CalcError boundaryError = parseError("1+2+3+4", limits);
        tests.expect(
            boundaryError.type() == error::CalcErrorType::ResourceLimit
                && boundaryError.span()
                && boundaryError.span()->begin.column == 6,
            "parser gives an explicit resource error at the operator-chain boundary");
    }

    {
        syntax::ParseLimits limits;
        limits.maxTokens = 3;
        tests.expect(
            parseError("1+2", limits).type() == error::CalcErrorType::ResourceLimit,
            "lexer token budget includes the end-of-input token");

        limits = {};
        limits.maxLiteralDigits = 3;
        tests.expect(
            parseError("1234", limits).type() == error::CalcErrorType::ResourceLimit,
            "lexer rejects an oversized numeric literal before numeric conversion");

        limits = {};
        limits.maxNodes = 2;
        tests.expect(
            parseError("1+2", limits).type() == error::CalcErrorType::ResourceLimit,
            "parser enforces the AST node budget");
    }

    {
        syntax::ParseLimits limits;
        limits.maxCallArguments = 2;
        tests.expect(
            parseError("f[1,2,3]", limits).type() == error::CalcErrorType::ResourceLimit,
            "parser enforces the function argument budget");

        limits = {};
        limits.maxArrayElements = 2;
        tests.expect(
            parseError("{1,2,3}", limits).type() == error::CalcErrorType::ResourceLimit,
            "parser enforces the array element budget");

        limits = {};
        limits.maxRecursionDepth = 32;
        tests.expect(
            parseError(nestedInput("(", ")", 100), limits).type()
                == error::CalcErrorType::ResourceLimit,
            "parser enforces its nesting budget with a source-spanned error");
    }

    // 旧実装がC++再帰で落ち得た形を，十分大きい動的corpusとして固定する。
    const auto expectHostileResourceLimit = [&](std::string input, std::string_view name) {
        tests.expect(
            parseError(std::move(input)).type() == error::CalcErrorType::ResourceLimit,
            name);
    };

    expectHostileResourceLimit(
        nestedInput("(", ")", 4000),
        "parser survives deeply nested groups without stack overflow");
    expectHostileResourceLimit(
        nestedInput("f[", "]", 4000),
        "parser survives deeply nested calls without stack overflow");
    expectHostileResourceLimit(
        nestedInput("{", "}", 4000),
        "parser survives deeply nested arrays without stack overflow");

    std::string unaryChain(50'000, '-');
    unaryChain += '1';
    expectHostileResourceLimit(
        std::move(unaryChain),
        "parser rejects a huge unary chain without recursive descent");

    std::string powerChain;
    powerChain.reserve(100'001);
    for (std::size_t index = 0; index < 50'000; ++index)
        powerChain += "1^";
    powerChain += '1';
    expectHostileResourceLimit(
        std::move(powerChain),
        "parser rejects a huge power chain without recursive descent");

    std::string additiveChain;
    additiveChain.reserve(100'001);
    for (std::size_t index = 0; index < 50'000; ++index)
        additiveChain += "1+";
    additiveChain += '1';
    expectHostileResourceLimit(
        std::move(additiveChain),
        "parser bounds a huge left-associative chain before lowering");

    std::string assignmentChain;
    assignmentChain.reserve(60'001);
    for (std::size_t index = 0; index < 20'000; ++index)
        assignmentChain += "x:=";
    assignmentChain += '1';
    expectHostileResourceLimit(
        std::move(assignmentChain),
        "parser rejects a huge assignment chain without recursive descent");

    std::string postfixChain{"1"};
    postfixChain.append(50'000, '!');
    expectHostileResourceLimit(
        std::move(postfixChain),
        "parser bounds huge postfix chains before lowering");
}

} // namespace mmcal::tests
