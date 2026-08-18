// トークン化lexer
#include "lexer.hpp"

#include "error/error_message.hpp"

#include <cctype>
#include <string>
#include <utility>

namespace mmcal::syntax {
namespace {

[[nodiscard]] bool isAsciiLetter(char value) noexcept {
    return (value >= 'A' && value <= 'Z') || (value >= 'a' && value <= 'z');
}

[[nodiscard]] bool isIdentifierStart(char value) noexcept {
    return isAsciiLetter(value) || value == '_';
}

[[nodiscard]] bool isIdentifierContinue(char value) noexcept {
    return isIdentifierStart(value) || (value >= '0' && value <= '9');
}

[[nodiscard]] bool isRadixDigit(char value) noexcept {
    return (value >= '0' && value <= '9')
        || (value >= 'A' && value <= 'Z')
        || (value >= 'a' && value <= 'z');
}

[[nodiscard]] std::size_t numericLiteralDigitCount(std::string_view text) noexcept {
    const bool hashRadix = text.find('#') != std::string_view::npos;
    const bool prefixedRadix = text.size() > 2 && text.front() == '0'
        && (text[1] == 'b' || text[1] == 'B'
            || text[1] == 'o' || text[1] == 'O'
            || text[1] == 'x' || text[1] == 'X');

    std::size_t count = 0;
    for (std::size_t index = 0; index < text.size(); ++index) {
        const char value = text[index];
        if (value >= '0' && value <= '9') {
            ++count;
            continue;
        }
        if ((hashRadix || (prefixedRadix && index >= 2)) && isAsciiLetter(value))
            ++count;
    }
    return count;
}

} // namespace

Lexer::Lexer(std::string_view sourceText, ParseBudget* budget)
    : sourceText_(sourceText),
      budget_(budget ? budget : &ownedBudget_) {}

std::vector<Token> Lexer::tokenize() {
    std::vector<Token> tokens;

    const auto append = [&](Token token) {
        budget_->checkTokenCount(tokens.size() + 1, token.span);
        if (token.kind == TokenKind::Number)
            budget_->checkLiteralDigits(
                numericLiteralDigitCount(token.text(sourceText_)), token.span);
        tokens.push_back(std::move(token));
    };

    while (true) {
        skipWhitespace();
        const source::SourcePosition begin = position();

        if (atEnd()) {
            append(makeToken(TokenKind::End, begin));
            return tokens;
        }

        const char value = current();
        if ((value >= '0' && value <= '9')
            || (value == '.' && peek() >= '0' && peek() <= '9')) {
            append(lexNumber());
            continue;
        }

        if (isIdentifierStart(value)) {
            append(lexIdentifier());
            continue;
        }

        if (value == '"') {
            append(lexString());
            continue;
        }

        advance();
        switch (value) {
        case '+':
            append(makeToken(TokenKind::Plus, begin));
            break;
        case '-':
            append(makeToken(TokenKind::Minus, begin));
            break;
        case '*':
            append(makeToken(TokenKind::Star, begin));
            break;
        case '/':
            append(makeToken(TokenKind::Slash, begin));
            break;
        case '^':
            append(makeToken(TokenKind::Caret, begin));
            break;
        case '(':
            append(makeToken(TokenKind::LParen, begin));
            break;
        case ')':
            append(makeToken(TokenKind::RParen, begin));
            break;
        case '[':
            append(makeToken(TokenKind::LBracket, begin));
            break;
        case ']':
            append(makeToken(TokenKind::RBracket, begin));
            break;
        case '{':
            append(makeToken(TokenKind::LBrace, begin));
            break;
        case '}':
            append(makeToken(TokenKind::RBrace, begin));
            break;
        case ',':
            append(makeToken(TokenKind::Comma, begin));
            break;
        case '%':
            append(makeToken(TokenKind::Percent, begin));
            break;
        case '@':
            append(makeToken(TokenKind::At, begin));
            break;
        case ':':
            if (current() != '=')
                error::throwCalcError(
                    error::CalcErrorType::Syntax,
                    "Expected '=' after ':'",
                    makeToken(TokenKind::End, begin).span);
            advance();
            append(makeToken(TokenKind::Assign, begin));
            break;
        case '<':
            if (current() == '=') {
                advance();
                append(makeToken(TokenKind::LessEqual, begin));
            }
            else
                append(makeToken(TokenKind::Less, begin));
            break;
        case '>':
            if (current() == '=') {
                advance();
                append(makeToken(TokenKind::GreaterEqual, begin));
            }
            else
                append(makeToken(TokenKind::Greater, begin));
            break;
        case '=':
            if (current() != '=')
                error::throwCalcError(
                    error::CalcErrorType::Syntax,
                    "Assignment uses ':='; equality uses '=='",
                    makeToken(TokenKind::End, begin).span);
            advance();
            append(makeToken(TokenKind::EqualEqual, begin));
            break;
        case '!':
            if (current() == '=') {
                advance();
                append(makeToken(TokenKind::BangEqual, begin));
            }
            else
                append(makeToken(TokenKind::Bang, begin));
            break;
        default:
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                std::string{"Unexpected character: "} + value,
                makeToken(TokenKind::End, begin).span);
        }
    }
}

bool Lexer::atEnd() const noexcept {
    return index_ >= sourceText_.size();
}

char Lexer::current() const noexcept {
    return atEnd() ? '\0' : sourceText_[index_];
}

char Lexer::peek(std::size_t distance) const noexcept {
    const std::size_t target = index_ + distance;
    return target >= sourceText_.size() ? '\0' : sourceText_[target];
}

source::SourcePosition Lexer::position() const noexcept {
    return source::SourcePosition{index_, line_, column_};
}

char Lexer::advance() noexcept {
    const char value = current();
    if (atEnd())
        return value;

    ++index_;
    if (value == '\n') {
        ++line_;
        column_ = 1;
    }
    else
        ++column_;

    return value;
}

void Lexer::skipWhitespace() noexcept {
    while (!atEnd() && std::isspace(static_cast<unsigned char>(current())) != 0)
        advance();
}

Token Lexer::lexNumber() {
    const source::SourcePosition begin = position();

    if (current() == '.') {
        advance();
        while (current() >= '0' && current() <= '9')
            advance();
    }
    else {
        while (current() >= '0' && current() <= '9')
            advance();

        if (current() == '#') {
            advance();
            bool hasDigit = false;
            bool hasPoint = false;

            while (isRadixDigit(current()) || current() == '.') {
                if (current() == '.') {
                    if (hasPoint)
                        error::throwCalcError(
                            error::CalcErrorType::Syntax,
                            "Number contains multiple radix points",
                            source::SourceSpan{begin, position()});
                    hasPoint = true;
                }
                else
                    hasDigit = true;

                advance();
            }

            if (!hasDigit)
                error::throwCalcError(
                    error::CalcErrorType::Syntax,
                    "Radix literal requires at least one digit",
                    source::SourceSpan{begin, position()});

            return makeToken(TokenKind::Number, begin);
        }

        if (index_ - begin.offset == 1
            && sourceText_[begin.offset] == '0'
            && (current() == 'b' || current() == 'B'
                || current() == 'o' || current() == 'O'
                || current() == 'x' || current() == 'X')) {
            advance();
            const std::size_t digitsBegin = index_;
            while (isRadixDigit(current()))
                advance();

            if (index_ == digitsBegin)
                error::throwCalcError(
                    error::CalcErrorType::Syntax,
                    "Prefixed integer requires at least one digit",
                    source::SourceSpan{begin, position()});

            return makeToken(TokenKind::Number, begin);
        }

        if (current() == '.') {
            advance();
            while (current() >= '0' && current() <= '9')
                advance();
        }
    }

    if (current() == 'e' || current() == 'E') {
        const char next = peek();
        const bool hasSign = next == '+' || next == '-';
        const bool startsExponent = (next >= '0' && next <= '9')
            || (hasSign && peek(2) >= '0' && peek(2) <= '9');

        // e/E は指数部が実際に続く場合だけ数値へ取り込む。
        // 2exp[x] や 2E^x は暗黙乗算として後続identifierへ分離する。
        if (startsExponent) {
            advance();
            if (current() == '+' || current() == '-')
                advance();
            while (current() >= '0' && current() <= '9')
                advance();
        }
        else if (hasSign) {
            advance();
            advance();
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "Exponent requires at least one digit",
                source::SourceSpan{begin, position()});
        }
    }

    if (current() == '.')
        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "Malformed number literal",
            source::SourceSpan{begin, position()});

    return makeToken(TokenKind::Number, begin);
}

Token Lexer::lexIdentifier() {
    const source::SourcePosition begin = position();
    while (isIdentifierContinue(current()))
        advance();
    return makeToken(TokenKind::Identifier, begin);
}

Token Lexer::lexString() {
    const source::SourcePosition begin = position();
    advance();

    bool escaped = false;
    while (!atEnd()) {
        const char value = advance();

        if (escaped) {
            escaped = false;
            continue;
        }

        if (value == '\\') {
            escaped = true;
            continue;
        }

        if (value == '"')
            return makeToken(TokenKind::String, begin);
    }

    error::throwCalcError(
        error::CalcErrorType::Syntax,
        "Unterminated string literal",
        source::SourceSpan{begin, position()});
}

Token Lexer::makeToken(TokenKind kind, source::SourcePosition begin) const noexcept {
    return Token{kind, source::SourceSpan{begin, position()}};
}

} // namespace mmcal::syntax
