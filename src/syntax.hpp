#pragma once

#include "core.hpp"
#include "functions.hpp"
#include "math_util.hpp"
#include <cctype>
#include <charconv>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mm::cal {
 /* ============================
Token / Lexer
============================ */

 enum class TokenType {
  End,
  Number,
  Identifier,
  Plus,  // +
  Minus, // -
  Mul,   // *
  Div,   // /
  Pow,   // ^
  // Star,
  // Slash,
  // Caret,
  LParen,    // (
  RParen,    // )
  LBracket,  // [
  RBracket,  // ]
  Comma,     // ,
  Bang,      // !
  Percent,   // %
  Less,      // <
  LessEq,    // <=
  Greater,   // >
  GreaterEq, // >=
  Equal,     // ==
  NotEqual,  // !=
 };

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
   struct ASTNode;
   Lexer lex;
   Token cur;
   Token prev{TokenType::End, "", 0};
   SystemConfig &cfg;

   Parser(SystemConfig &cfg, const std::string &s);
   void advance();
   const Token &peek();

   bool accept(TokenType t);
   void expect(TokenType t);

   bool startsPrimary() const;

   bool isFunctionName(const std::string &name) const;
   bool isConstantName(const std::string &name) const;
   bool isImplicitMul(const Token &prev, const Token &cur) const;

   std::unique_ptr<ASTNode> parseCompare();
   std::unique_ptr<ASTNode> parseExpression();
   std::unique_ptr<ASTNode> parseTerm();
   std::unique_ptr<ASTNode> parsePower();
   std::unique_ptr<ASTNode> parseUnary();
   std::unique_ptr<ASTNode> parsePrimary();
   std::unique_ptr<ASTNode> parsePostfix();

   /* ============================
      AST
      ============================ */

   enum class UnaryOp { Plus, Minus };
   enum class BinOp { Add, Sub, Mul, Div, Pow };
   enum class CmpOp { Less, LessEq, Greater, GreaterEq, Equal, NotEqual };

   static CmpOp tokenToCmpOp(TokenType t);
   static BinOp tokenToBinOp(TokenType t);

   struct ASTNode {
     size_t pos = 0;
     virtual ~ASTNode() = default;

     // 外から呼ばれる唯一の eval
     virtual Value eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const final {
      Value v = evalImpl(cfg, hist, base);
      checkFiniteValue(v, pos); // nan / inf は必ずここで捕まる
      return v;
     }

    protected:
     // 派生クラスはこれだけ実装する
     virtual Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const = 0;

     // 派生クラスからも使える
     static void checkFinite(const Value &v, size_t pos) {
      std::visit(
          [&](auto &&x) {
           using T = std::decay_t<decltype(x)>;

           if constexpr (std::is_same_v<T, double>) {
            if (std::isnan(x)) throw CalcError(CalcErrorType::NaNDetected, errorMessage(CalcErrorType::NaNDetected), pos);
            if (!std::isfinite(x)) throw CalcError(CalcErrorType::InfinityDetected, errorMessage(CalcErrorType::InfinityDetected), pos);
           } else if constexpr (std::is_same_v<T, std::complex<double>>) {
            if (std::isnan(x.real()) || std::isnan(x.imag())) throw CalcError(CalcErrorType::NaNDetected, errorMessage(CalcErrorType::NaNDetected), pos);
            if (!std::isfinite(x.real()) || !std::isfinite(x.imag())) throw CalcError(CalcErrorType::InfinityDetected, errorMessage(CalcErrorType::InfinityDetected), pos);
           }
          },
          v);
     }
   };

   /* ---------- Number ---------- */
   struct NumberNode : ASTNode {
     Value value;

     NumberNode(Value v, size_t p);
     Value evalImpl(SystemConfig &, const std::vector<InputEntry> &, int) const override;
   };

   /* ---------- Binary ---------- */

   struct BinaryNode : ASTNode {
     BinOp op;
     std::unique_ptr<ASTNode> lhs, rhs;

     BinaryNode(BinOp o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p);

     Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override;
   };

   /* ---------- Unary ---------- */
   struct UnaryNode : ASTNode {
     UnaryOp op;
     std::unique_ptr<ASTNode> rhs;
     UnaryNode(UnaryOp o, std::unique_ptr<ASTNode> r, size_t p);

     Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override;
   };

   /* ---------- Postfix ---------- */
   struct PostfixUnaryNode : ASTNode {
     char op;
     std::unique_ptr<ASTNode> expr;

     PostfixUnaryNode(char op, std::unique_ptr<ASTNode> e, size_t p);

     Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override;
   };

   /* ---------- Compare ---------- */
   struct CompareNode : ASTNode {
     CmpOp op;
     std::unique_ptr<ASTNode> lhs, rhs;

     CompareNode(CmpOp o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p);

     Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override;
   };

   /* ---------- Function ---------- */

   struct FunctionCallNode : ASTNode {
     std::string name;
     std::vector<std::unique_ptr<ASTNode>> args;

     Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override;
   };

   /* ---------- History % ---------- */
   //
   // struct OutRelativeNode : ASTNode {
   //  int offset; // % の個数
   //
   //  OutRelativeNode(int n, size_t p) : offset(n) ;
   //
   //  Value eval(SystemConfig &, const std::vector<InputEntry> &hist, int base) const override {
   //   if (offset <= 0 || offset > (int)hist.size()) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos);
   //
   //   // 相対参照
   //   return {hist[hist.size() - offset].value, AngleUnit::None};
   //  }
   //};
   //

   /* ---------- History Out[n] ---------- */
   struct OutNode : ASTNode {
     int index; // 正: Out[n], 負: % / %%
     OutNode(int idx, size_t p);

     Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override;
   };

   /* ============================
      ユーティリティ先輩
      ============================ */

   inline bool isValueEnd(TokenType t);
   inline bool isValueStart(TokenType t);

 }; // class Parser
 Value evalCompare(const Value &lhs, const Value &rhs, Parser::CmpOp op, FunctionContext &ctx);
} // namespace mm::cal