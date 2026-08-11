// トークン化の回帰テスト
#include "lexer_tests.hpp"

#include "error/error_message.hpp"
#include "syntax/lexer.hpp"
#include "test_framework.hpp"

#include <string>
#include <vector>

namespace mmcal::tests {
namespace {

[[nodiscard]] std::vector<syntax::TokenKind> tokenKinds(std::string_view sourceText) {
    syntax::Lexer lexer{sourceText};
    const auto tokens = lexer.tokenize();
    std::vector<syntax::TokenKind> kinds;
    kinds.reserve(tokens.size());
    for (const syntax::Token& token : tokens)
        kinds.push_back(token.kind);
    return kinds;
}

} // namespace

void runLexerTests(TestRunner& tests) {
    using syntax::TokenKind;

    tests.expect(
        tokenKinds("x := 1/3") == std::vector<TokenKind>{
            TokenKind::Identifier,
            TokenKind::Assign,
            TokenKind::Number,
            TokenKind::Slash,
            TokenKind::Number,
            TokenKind::End},
        "lexer tokenizes assignment");

    tests.expect(
        tokenKinds("sqrt[16#A.F]") == std::vector<TokenKind>{
            TokenKind::Identifier,
            TokenKind::LBracket,
            TokenKind::Number,
            TokenKind::RBracket,
            TokenKind::End},
        "lexer tokenizes bracket call and radix number");

    tests.expect(
        tokenKinds(".5e-2 + 0xFF + 0b101") == std::vector<TokenKind>{
            TokenKind::Number,
            TokenKind::Plus,
            TokenKind::Number,
            TokenKind::Plus,
            TokenKind::Number,
            TokenKind::End},
        "lexer tokenizes exact numeric forms");

    tests.expect(
        tokenKinds("2exp[x]+2E^x") == std::vector<TokenKind>{
            TokenKind::Number,
            TokenKind::Identifier,
            TokenKind::LBracket,
            TokenKind::Identifier,
            TokenKind::RBracket,
            TokenKind::Plus,
            TokenKind::Number,
            TokenKind::Identifier,
            TokenKind::Caret,
            TokenKind::Identifier,
            TokenKind::End},
        "lexer separates e/E identifiers from numbers when no decimal exponent follows");

    tests.expect(
        tokenKinds("a<=b!=c==d>=e") == std::vector<TokenKind>{
            TokenKind::Identifier,
            TokenKind::LessEqual,
            TokenKind::Identifier,
            TokenKind::BangEqual,
            TokenKind::Identifier,
            TokenKind::EqualEqual,
            TokenKind::Identifier,
            TokenKind::GreaterEqual,
            TokenKind::Identifier,
            TokenKind::End},
        "lexer tokenizes comparison operators");

    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(syntax::Lexer{"x = 1"}.tokenize()); },
        "lexer rejects single equals");
    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(syntax::Lexer{"\"unterminated"}.tokenize()); },
        "lexer rejects unterminated string");
    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(syntax::Lexer{"1e+"}.tokenize()); },
        "lexer rejects missing exponent digits");
    tests.expectThrows<error::CalcError>(
        [] { static_cast<void>(syntax::Lexer{"1.2.3"}.tokenize()); },
        "lexer rejects multiple decimal points");
}

} // namespace mmcal::tests
