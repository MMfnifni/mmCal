// 構文木からExprへのloweringの回帰テスト
#include "lowerer_tests.hpp"

#include "error/error_message.hpp"
#include "evaluation/builtin_registry.hpp"
#include "evaluation/environment.hpp"
#include "evaluation/evaluator.hpp"
#include "formatting/expr_formatter.hpp"
#include "syntax/lexer.hpp"
#include "syntax/lowerer.hpp"
#include "syntax/parser.hpp"
#include "test_framework.hpp"

#include <cstdint>
#include <memory>
#include <vector>
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

class DeterministicGenerator final {
public:
    [[nodiscard]] std::uint64_t next() noexcept {
        state_ ^= state_ << 13;
        state_ ^= state_ >> 7;
        state_ ^= state_ << 17;
        return state_;
    }

    [[nodiscard]] std::size_t choose(std::size_t count) noexcept {
        return static_cast<std::size_t>(next() % count);
    }

private:
    std::uint64_t state_ = 0x6d6d43616c150001ULL;
};

[[nodiscard]] expression::Expr generatedExpression(
    DeterministicGenerator& random,
    int depth,
    const evaluation::BuiltinRegistry& builtins,
    symbols::SymbolTable& symbols,
    bool allowArray = true) {
    using evaluation::BuiltinId;
    using expression::Expr;
    using numeric::BigInt;
    using numeric::Number;

    const auto atom = [&]() -> Expr {
        switch (random.choose(7)) {
        case 0: return Expr{Number{BigInt{static_cast<std::int64_t>(random.choose(19)) - 9}}};
        case 1: return Expr{numeric::Rational{
            BigInt{static_cast<std::int64_t>(random.choose(17)) - 8},
            BigInt{static_cast<std::int64_t>(random.choose(8)) + 1}}};
        case 2: return Expr{symbols.intern("x")};
        case 3: return Expr{symbols.intern("y")};
        case 4: return Expr{symbols.intern("e")};
        case 5: return Expr{symbols.intern("Pi")};
        default: return Expr{symbols.intern("E")};
        }
    };

    if (depth <= 0 || random.choose(5) == 0)
        return atom();

    const auto recurse = [&]() {
        return generatedExpression(random, depth - 1, builtins, symbols, allowArray);
    };
    const auto call = [&](BuiltinId id, std::vector<Expr> arguments) {
        return Expr::call(builtins.symbol(id), std::move(arguments));
    };

    switch (random.choose(allowArray ? 12 : 11)) {
    case 0:
        return call(BuiltinId::Add, {recurse(), recurse(), recurse()});
    case 1:
        return call(BuiltinId::Subtract, {recurse(), recurse()});
    case 2:
        return call(BuiltinId::Multiply, {recurse(), recurse(), recurse()});
    case 3:
        return call(BuiltinId::Divide, {recurse(), recurse()});
    case 4:
        return call(BuiltinId::Power, {recurse(), Expr{Number{BigInt{
            static_cast<std::int64_t>(random.choose(5)) - 2}}}});
    case 5:
        return call(BuiltinId::Negate, {recurse()});
    case 6:
        return call(BuiltinId::Sin, {recurse()});
    case 7:
        return call(BuiltinId::Cos, {recurse()});
    case 8:
        return call(BuiltinId::Exp, {recurse()});
    case 9:
        return call(BuiltinId::Sqrt, {recurse()});
    case 10:
        return call(BuiltinId::Log, {recurse()});
    default:
        return Expr::array({3}, {
            generatedExpression(random, depth - 1, builtins, symbols, false),
            generatedExpression(random, depth - 1, builtins, symbols, false),
            generatedExpression(random, depth - 1, builtins, symbols, false)});
    }
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
    tests.expectEqual(lowerAndFormat("x(x+1)"), std::string{"x*(x+1)"},
        "formatter makes identifier/group multiplication explicit");
    tests.expectEqual(lowerAndFormat("(x^2)^3"), std::string{"(x^2)^3"},
        "formatter preserves left-nested Power associativity");
    tests.expectEqual(lowerAndFormat("x^(2^3)"), std::string{"x^2^3"},
        "formatter uses the parser's right-associative Power notation");
    tests.expectEqual(lowerAndFormat("f[2]"), std::string{"f[2]"},
        "square bracket syntax lowers to a function call");
    tests.expectEqual(lowerAndFormat("f[x]:=x+1"), std::string{"f[x]:=x+1"},
        "function definition syntax is preserved for later evaluation");
    tests.expectEqual(lowerAndFormat("1 < x <= 3"), std::string{"And[1<x, x<=3]"},
        "comparison chain lowers without losing operands");
    tests.expectEqual(lowerAndFormat("%%"), std::string{"%%"},
        "output history reference depth is preserved");
    tests.expectEqual(lowerAndFormat("@"), std::string{"In[-1]"},
        "input history shorthand lowers to the formal relative reference");
    tests.expectEqual(lowerAndFormat("30deg"), std::string{"30 deg"},
        "unit suffix is represented explicitly");
    tests.expectEqual(lowerAndFormat("\"a\\nb\""), std::string{"\"a\\nb\""},
        "string literal round trip");
    tests.expectEqual(lowerAndFormat("{{1, 2}, {3, 4}}"), std::string{"{{1, 2}, {3, 4}}"},
        "rectangular array shape is inferred");

    tests.expectEqual(lowerAndFormat("{1, {2, 3}}"), std::string{"{1, {2, 3}}"},
        "lowerer preserves a non-rectangular brace value");
    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(lower("2#102")); },
        "lowerer rejects a digit outside the radix");

    // formatterが出した式をparser/lowererへ戻し、再formatして固定点になることを生成式で検証する。
    auto& symbolTable = symbols::defaultSymbolTable();
    const auto builtins = evaluation::BuiltinRegistry::defaults(symbolTable);
    DeterministicGenerator random;
    bool generatedRoundTrip = true;
    std::string failedInput;
    std::string failedOutput;
    for (std::size_t i = 0; i < 3000; ++i) {
        const expression::Expr generated = generatedExpression(random, 4, builtins, symbolTable);
        const std::string formatted = formatting::formatExpr(generated);
        try {
            const std::string reformatted = formatting::formatExpr(lower(formatted));
            if (reformatted != formatted) {
                generatedRoundTrip = false;
                failedInput = formatted;
                failedOutput = reformatted;
                break;
            }
        }
        catch (const error::CalcError& error) {
            generatedRoundTrip = false;
            failedInput = formatted;
            failedOutput = error.what();
            break;
        }
    }
    tests.expect(generatedRoundTrip,
        "formatter/parser: 3000 generated expressions reach a format-parse-format fixed point");
    if (!generatedRoundTrip) {
        tests.expectEqual(failedOutput, failedInput,
            "formatter/parser: generated round-trip failure detail");
    }

    const std::vector<std::string> lexicalBoundaries{
        "2exp[x]", "2E", "2E^x", "2e3", "2e-3", "-2^2", "(-2)^2",
        "2sqrt[2]", "(x+1)(x-1)", "x!^2", "1/(x/y)", "{{1,2},{3,4}}"};
    bool lexicalRoundTrip = true;
    for (const std::string& source : lexicalBoundaries) {
        const std::string once = formatting::formatExpr(lower(source));
        if (formatting::formatExpr(lower(once)) != once) {
            lexicalRoundTrip = false;
            break;
        }
    }
    tests.expect(lexicalRoundTrip,
        "formatter/parser: lexical ambiguity boundary cases round-trip canonically");
}

} // namespace mmcal::tests
