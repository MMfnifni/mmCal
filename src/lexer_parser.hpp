#pragma once

#include "ast.hpp"
#include "core.hpp"
#include "math_util.hpp"
#include "value_arithmetic.hpp"
#include "value_compare.hpp"
#include <cctype>
#include <charconv>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mm::cal {
 /* ============================
Token / Lexer
============================ */

 struct Token {
   TokenType type = TokenType::End;
   std::string text;
   size_t pos = 0;
   double number = 0.0;
   Token() = default;

   Token(TokenType t, std::string txt, size_t p) : type(t), text(std::move(txt)), number(0.0), pos(p) {}

   Token(TokenType t, std::string txt, double num, size_t p) : type(t), text(std::move(txt)), number(num), pos(p) {}
 };
 class Lexer {
  public:
   std::string_view src;
   size_t p = 0;

   Token lookahead{};
   bool hasPeek = false;
   explicit Lexer(std::string_view s) : src(s) {}
   const Token &peek();
   Token get();
   bool eof();

  private:
   void skipSpaces();
   Token nextToken();
 };

 /* ============================
   Parser（暗黙乗算）
   ============================ */

 class Parser {
  public:
   Lexer lex;
   Token cur;
   Token prev{TokenType::End, "", 0};
   SystemConfig &cfg;

   Parser(SystemConfig &cfg, const std::string &s);
   void advance();
   // const Token &peek();

   bool accept(TokenType t);
   void expect(TokenType t);

   bool startsPrimary() const;

   bool isFunctionName(const std::string &name) const;
   bool isConstantName(const std::string &name) const;
   bool isUnitName(const std::string &s) const;
   bool isImplicitMul(const Token &prev, const Token &cur) const;

   bool isValueEnd(TokenType t);
   bool isValueStart(TokenType t);

   //  parse
   // └ parseCompare
   //    └ parseExpression
   //      └ parseTerm
   //        └ parseUnary
   //          └ parsePower
   //            └ parsePostfix
   //              └ parsePrimar

   std::unique_ptr<ASTNode> parse(); // パースの入口
   std::unique_ptr<ASTNode> parseCompare();
   std::unique_ptr<ASTNode> parseExpression();
   std::unique_ptr<ASTNode> parseTerm();
   std::unique_ptr<ASTNode> parsePower();
   std::unique_ptr<ASTNode> parseUnary();
   std::unique_ptr<ASTNode> parsePrimary();
   std::unique_ptr<ASTNode> parsePostfix();
   std::unique_ptr<ASTNode> parseAssignment();

 }; // class Parser
 Value evalCompare(const Value &lhs, const Value &rhs, CmpOp op, FunctionContext &ctx);

} // namespace mm::cal