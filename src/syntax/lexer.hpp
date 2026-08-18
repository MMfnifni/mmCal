#pragma once

#include "parse_budget.hpp"
#include "token.hpp"

#include <string_view>
#include <vector>

namespace mmcal::syntax {

// 文字列を意味解釈せず、位置情報付きトークン列へ分解する。
class Lexer final {
public:
    explicit Lexer(std::string_view sourceText, ParseBudget* budget = nullptr);

    [[nodiscard]] std::vector<Token> tokenize();

private:
    std::string_view sourceText_;
    ParseBudget ownedBudget_;
    ParseBudget* budget_ = nullptr;
    std::size_t index_ = 0;
    std::size_t line_ = 1;
    std::size_t column_ = 1;

    [[nodiscard]] bool atEnd() const noexcept;
    [[nodiscard]] char current() const noexcept;
    [[nodiscard]] char peek(std::size_t distance = 1) const noexcept;
    [[nodiscard]] source::SourcePosition position() const noexcept;

    char advance() noexcept;
    void skipWhitespace() noexcept;
    [[nodiscard]] Token lexNumber();
    [[nodiscard]] Token lexIdentifier();
    [[nodiscard]] Token lexString();
    [[nodiscard]] Token makeToken(TokenKind kind, source::SourcePosition begin) const noexcept;
};

} // namespace mmcal::syntax
