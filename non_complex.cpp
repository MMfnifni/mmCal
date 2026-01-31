#include <cassert>
#include <cctype>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

/* ============================
   エラー
   ============================ */

enum class CalcErrorType {
 Syntax,
 UnknownIdentifier,
 InvalidCharacter,
 MismatchedParen,
 OperandMissing,
 FunctionMissing,
 DivisionByZero,
 AngleWithoutdouble,
 DomainError,
 NaNResult,
 InfiniteResult,
 OutOfRange,
 NeedInteger,
 Overflow,
 InvalidNumber,
 InvalidOperation,
};

struct CalcError : std::runtime_error {
  CalcErrorType type;
  size_t pos;
  CalcError(CalcErrorType t, const std::string &m, size_t p) : std::runtime_error(m), type(t), pos(p) {}
};

std::string errorMessage(CalcErrorType t) {
 switch (t) {
  case CalcErrorType::Syntax: return "syntax error";
  case CalcErrorType::UnknownIdentifier: return "unknown identifier";
  case CalcErrorType::InvalidCharacter: return "invalid character";
  case CalcErrorType::MismatchedParen: return "mismatched parentheses";
  case CalcErrorType::OperandMissing: return "operand missing";
  case CalcErrorType::FunctionMissing: return "function argument missing";
  case CalcErrorType::DivisionByZero: return "division by zero";
  case CalcErrorType::AngleWithoutdouble: return "angle without value";
  case CalcErrorType::DomainError: return "domain error";
  case CalcErrorType::NaNResult: return "result is NaN";
  case CalcErrorType::InfiniteResult: return "result is infinite";
  case CalcErrorType::OutOfRange: return "history out of range";
  case CalcErrorType::NeedInteger: return "this func need integer";
  case CalcErrorType::Overflow: return "overflow detected";
  case CalcErrorType::InvalidNumber: return "invalid number";
  case CalcErrorType::InvalidOperation: return "onvalid operation";
 }
 return "error";
}

/* ============================
   基本定義
   ============================ */
struct SystemConfig;
struct FunctionContext;
struct FunctionDef;
struct InputEntry;

using FuncImpl = std::function<double(const std::vector<double> &, FunctionContext &)>;

struct FunctionDef {
  int minArgc = 0;
  int maxArgc = -1; // -1 = unlimited
  FuncImpl f;

  bool validArgc(int argc) const {
   if (argc < minArgc) return false;
   if (maxArgc >= 0 && argc > maxArgc) return false;
   return true;
  }
};

inline double round12(double v) {
 constexpr double scale = 1e12;
 return std::round(v * scale) / scale;
}

struct SystemConfig {
  int fixedDigits = -1;
  std::unordered_map<std::string, FunctionDef> functions;
  bool shouldExit = false;
  bool shouldClear = false;
};

struct FunctionContext {
  SystemConfig &cfg;
  const std::vector<InputEntry> &hist;
  int base;
  size_t pos;
};

double factorial(double x, size_t pos) {
 if (!std::isfinite(x)) throw CalcError(CalcErrorType::DomainError, "factorial: invalid value", pos);

 // 整数チェック
 double ix;
 if (std::modf(x, &ix) != 0.0) throw CalcError(CalcErrorType::DomainError, "factorial: not an integer", pos);

 if (ix < 0) throw CalcError(CalcErrorType::DomainError, "factorial: negative value", pos);

 if (ix > 170) // 171! は double で overflow
  throw CalcError(CalcErrorType::Overflow, "factorial: overflow", pos);

 double r = 1.0;
 for (int i = 2; i <= (int)ix; ++i)
  r *= i;
 return r;
}

class ExitRequest {};

/* ============================
   Token / Lexer
   ============================ */

enum class TokenType {
 End,
 Number,
 Identifier,
 Plus,
 Minus,
 Mul,
 Div,
 Pow,
 Star,
 Slash,
 Caret,
 LParen,
 RParen,
 LBracket,
 RBracket,
 Comma,
 Bang,
 Percent,
};

struct Token {
  TokenType type;
  std::string text;
  double number = 0.0;
  size_t pos = 0;
};

class Lexer {
 public:
  std::string_view src;
  size_t p = 0;

  Token lookahead{};
  bool hasPeek = false;

  explicit Lexer(std::string_view s) : src(s) {}

  const Token &peek() {
   if (!hasPeek) {
    lookahead = nextToken();
    hasPeek = true;
   }
   return lookahead;
  }

  Token get() {
   if (hasPeek) {
    hasPeek = false;
    return lookahead;
   }
   return nextToken();
  }

  bool eof() { return peek().type == TokenType::End; }

 private:
  void skipSpaces() {
   while (p < src.size() && std::isspace((unsigned char)src[p]))
    ++p;
  }

  Token nextToken() {
   skipSpaces();
   size_t start = p;

   if (p >= src.size()) return {TokenType::End, "", 0.0, start};

   char c = src[p++];

   // number
   if (std::isdigit(c) || c == '.') {
    --p;
    char *end;
    double v = std::strtod(src.data() + p, &end);
    if (end == src.data() + p) { throw CalcError(CalcErrorType::InvalidNumber, "invalid number", p); }
    p = end - src.data();
    return {TokenType::Number, "", v, start};
   }

   // identifier
   if (std::isalpha(c) || c == '_') {
    std::string s(1, c);
    while (p < src.size() && (std::isalnum(src[p]) || src[p] == '_'))
     s += src[p++];
    return {TokenType::Identifier, s, 0.0, start};
   }

   switch (c) {
    case '+': return {TokenType::Plus, "+", 0, start};
    case '-': return {TokenType::Minus, "-", 0, start};
    case '*': return {TokenType::Mul, "*", 0, start};
    case '/': return {TokenType::Div, "/", 0, start};
    case '^': return {TokenType::Pow, "^", 0, start};
    case '(': return {TokenType::LParen, "(", 0, start};
    case ')': return {TokenType::RParen, ")", 0, start};
    case '[': return {TokenType::LBracket, "[", 0, start};
    case ']': return {TokenType::RBracket, "]", 0, start};
    case ',': return {TokenType::Comma, ",", 0, start};
    case '!': return {TokenType::Bang, "!", 0, start};
    case '%': return {TokenType::Percent, "%", 0, start};
   }

   throw CalcError(CalcErrorType::InvalidCharacter, std::string("unexpected character: ") + c, start);
  }
};

/* ============================
   定数・補助
   ============================ */

const std::unordered_map<std::string, double> constants = {
    { "Pi", 3.14159265358979323846},
    {  "E", 2.71828182845904523536},
    {"Teu", 6.28318530717958647692},
    {"Phi", 1.61803398874989484820},
    { "NA",          6.02214076e23},
    //{"teu", 6.28318530717958647692},
    //{"teu", 6.28318530717958647692},
    //{"teu", 6.28318530717958647692},
};

struct InputEntry {
  std::string expr; // In[n]
  double value;     // Out[n]
};

inline void checkFinite(double v, size_t pos) {
 if (std::isnan(v)) throw CalcError(CalcErrorType::NaNResult, errorMessage(CalcErrorType::NaNResult), pos);
 if (std::isinf(v)) throw CalcError(CalcErrorType::InfiniteResult, errorMessage(CalcErrorType::InfiniteResult), pos);
}

double radToDeg(double r) { return r * 180.0 / constants.at("Pi"); }

inline long long requireInt(double x, size_t pos) {
 if (!std::isfinite(x) || std::floor(x) != x) throw CalcError(CalcErrorType::NeedInteger, errorMessage(CalcErrorType::NeedInteger), pos);
 return (long long)x;
}

long long permLL(long long n, long long r) {
 if (r < 0 || r > n) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), 0);
 long long res = 1;
 for (long long i = 0; i < r; ++i)
  res *= (n - i);
 return res;
}

long long combLL(long long n, long long r) {
 if (r < 0 || r > n) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), 0);
 r = std::min(r, n - r);
 long long res = 1;
 for (long long i = 1; i <= r; ++i)
  res = res * (n - r + i) / i;
 return res;
}

long long gcdLL(long long a, long long b) {
 if (a < 0) a = -a;
 if (b < 0) b = -b;
 while (b != 0) {
  long long r = a % b;
  a = b;
  b = r;
 }
 return a;
}

double evaluate(const std::string &src, SystemConfig &cfg, const std::vector<InputEntry> &hist, int base);

/* ============================
   AST
   ============================ */

struct ASTNode {
  size_t pos;
  virtual ~ASTNode() = default;
  virtual double eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const = 0;
};

/* ---------- Number ---------- */

struct NumberNode : ASTNode {
  double value;

  NumberNode(double v, size_t p) : value(v) { pos = p; }
  double eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override { return value; }
};

/* ---------- Unary ---------- */

struct UnaryNode : ASTNode {
  char op;
  std::unique_ptr<ASTNode> rhs;
  UnaryNode(char o, std::unique_ptr<ASTNode> r, size_t p) : op(o), rhs(std::move(r)) { pos = p; }
  double eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
   auto v = rhs->eval(cfg, hist, base);
   if (op == '-') return -v;
   return v;
  }
};

struct PostfixUnaryNode : ASTNode {
  char op;
  std::unique_ptr<ASTNode> expr;

  PostfixUnaryNode(char op, std::unique_ptr<ASTNode> e, size_t p) : op(op), expr(std::move(e)) { pos = p; }

  double eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
   double v = expr->eval(cfg, hist, base);

   double r;
   switch (op) {
    case '!': r = factorial(v, pos); break;
    case '%': r = v / 100.0; break;
    default: throw CalcError(CalcErrorType::Syntax, "unknown postfix operator", pos);
   }

   checkFinite(r, pos);
   return r;
  }
};

/* ---------- Binary ---------- */

struct BinaryNode : ASTNode {
  char op;
  std::unique_ptr<ASTNode> lhs, rhs;
  BinaryNode(char o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p) : op(o), lhs(std::move(l)), rhs(std::move(r)) { pos = p; }
  double eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
   auto a = lhs->eval(cfg, hist, base);
   auto b = rhs->eval(cfg, hist, base);
   double r = 0.0;
   switch (op) {
    case '+': r = a + b; break;
    case '-': r = a - b; break;
    case '*': r = a * b; break;
    case '/':
     if (b == 0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
     r = a / b;
     break;
    case '^': r = std::pow(a, b); break;
   }
   checkFinite(r, pos);
   return r;
  }
};

/* ---------- Function ---------- */

struct FunctionCallNode : ASTNode {
  std::string name;
  std::vector<std::unique_ptr<ASTNode>> args;

  double eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
   auto it = cfg.functions.find(name);
   if (it == cfg.functions.end()) throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), pos);

   auto &f = it->second;
   int argc = (int)args.size();

   if (!f.validArgc(argc)) throw CalcError(CalcErrorType::FunctionMissing, errorMessage(CalcErrorType::FunctionMissing), pos);

   std::vector<double> v;
   for (auto &a : args)
    v.push_back(a->eval(cfg, hist, base));

   FunctionContext ctx{cfg, hist, base, pos};

   double r = f.f(v, ctx);
   checkFinite(r, pos);
   return r;
  }
};

/* ---------- History % ---------- */
//
// struct OutRelativeNode : ASTNode {
//  int offset; // % の個数
//
//  OutRelativeNode(int n, size_t p) : offset(n) { pos = p; }
//
//  double eval(SystemConfig &, const std::vector<InputEntry> &hist, int base) const override {
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
  OutNode(int idx, size_t p) : index(idx) { pos = p; }

  double eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
   if (index == 0) { throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos); }

   int i;
   if (index > 0) {
    i = index;
   } else {
    i = (int)hist.size() + index + 1; // -1 → last
   }

   if (i <= 0 || i > (int)hist.size()) { throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos); }

   return hist[i - 1].value;
  }
};

/* ============================
   Parser（暗黙乗算）
   ============================ */

class Parser {
 public:
  Lexer lex;
  Token cur;

  Parser(const std::string &s) : lex(s) { cur = lex.get(); }

  void advance() { cur = lex.get(); }

  const Token &peek() { return lex.peek(); }

  bool accept(TokenType t) {
   if (cur.type == t) {
    advance();
    return true;
   }
   return false;
  }

  void expect(TokenType t) {
   if (!accept(t)) throw CalcError(CalcErrorType::Syntax, errorMessage(CalcErrorType::Syntax), cur.pos);
  }

  bool startsPrimary() const { return cur.type == TokenType::Number || cur.type == TokenType::Identifier || cur.type == TokenType::LParen || cur.type == TokenType::Percent; }

  std::unique_ptr<ASTNode> parseExpression();

 private:
  std::unique_ptr<ASTNode> parseTerm();
  std::unique_ptr<ASTNode> parsePower();
  std::unique_ptr<ASTNode> parseUnary();
  std::unique_ptr<ASTNode> parsePrimary();
  std::unique_ptr<ASTNode> parsePostfix();
  bool acceptParen() {
   if (cur.type == TokenType::LParen || cur.type == TokenType::LBracket) {
    advance();
    return true;
   }
   return false;
  }

  void expectParenClose(TokenType open) {
   if (open == TokenType::LParen) expect(TokenType::RParen);
   else expect(TokenType::RBracket);
  }
};

std::unique_ptr<ASTNode> Parser::parsePostfix() {
 auto node = parsePrimary();

 while (cur.type == TokenType::Bang) {
  size_t p = cur.pos;
  advance();
  node = std::make_unique<PostfixUnaryNode>('!', std::move(node), p);
 }

 return node;
}

std::unique_ptr<ASTNode> Parser::parseExpression() {
 auto n = parseTerm();
 while (cur.type == TokenType::Plus || cur.type == TokenType::Minus) {
  char op = cur.text[0];
  size_t opPos = cur.pos;
  advance();
  auto r = parseTerm();
  n = std::make_unique<BinaryNode>(op, std::move(n), std::move(r), cur.pos);
 }
 return n;
}

std::unique_ptr<ASTNode> Parser::parseTerm() {
 auto n = parsePower();
 while (true) {
  if (cur.type == TokenType::Mul || cur.type == TokenType::Div) {
   char op = cur.text[0];
   size_t opPos = cur.pos;
   advance();
   auto r = parsePower();
   n = std::make_unique<BinaryNode>(op, std::move(n), std::move(r), cur.pos);
  }
  // 暗黙乗算
  else if (startsPrimary()) {
   size_t opPos = cur.pos;
   auto r = parsePower();
   n = std::make_unique<BinaryNode>('*', std::move(n), std::move(r), cur.pos);
  } else break;
 }
 return n;
}

std::unique_ptr<ASTNode> Parser::parsePower() {
 auto n = parseUnary();
 if (cur.type == TokenType::Pow) {
  size_t opPos = cur.pos;
  advance();
  auto r = parsePower();
  n = std::make_unique<BinaryNode>('^', std::move(n), std::move(r), cur.pos);
 }
 return n;
}

std::unique_ptr<ASTNode> Parser::parseUnary() {
 if (cur.type == TokenType::Plus || cur.type == TokenType::Minus) {
  char op = cur.text[0];
  size_t p = cur.pos;
  advance();
  return std::make_unique<UnaryNode>(op, parseUnary(), p);
 }
 return parsePostfix();
}

std::unique_ptr<ASTNode> Parser::parsePrimary() {
 size_t p = cur.pos;

 // ---- number ----
 if (cur.type == TokenType::Number) {
  double v = cur.number;
  advance();
  return std::make_unique<NumberNode>(v, p);
 }

 // ---- history % ----
 if (cur.type == TokenType::Percent) {
  int c = 0;
  while (cur.type == TokenType::Percent) {
   advance();
   ++c;
  }
  return std::make_unique<OutNode>(-c, p);
 }

 // ---- identifier ----
 if (cur.type == TokenType::Identifier) {
  std::string name = cur.text;
  advance();

  // function call: () or []
  TokenType open;
  if (cur.type == TokenType::LParen || cur.type == TokenType::LBracket) {
   open = cur.type;
   advance();
  } else {
   // constant
   if (constants.count(name)) return std::make_unique<NumberNode>(constants.at(name), p);

   throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), p);
  }

  auto f = std::make_unique<FunctionCallNode>();
  f->name = name;
  f->pos = p;

  // no args
  if ((open == TokenType::LParen && accept(TokenType::RParen)) || (open == TokenType::LBracket && accept(TokenType::RBracket))) { return f; }

  // args
  while (true) {
   f->args.push_back(parseExpression());
   if (accept(TokenType::Comma)) continue;

   if (open == TokenType::LParen) expect(TokenType::RParen);
   else expect(TokenType::RBracket);
   break;
  }

  return f;
 }

 // ---- grouped expression ----
 if (cur.type == TokenType::LParen || cur.type == TokenType::LBracket) {
  TokenType open = cur.type;
  advance();

  auto n = parseExpression();

  if (open == TokenType::LParen) expect(TokenType::RParen);
  else expect(TokenType::RBracket);

  return n;
 }

 throw CalcError(CalcErrorType::Syntax, errorMessage(CalcErrorType::Syntax), cur.pos);
}

/* ============================
   評価
   ============================ */

double evaluate(const std::string &src, SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) {
 Parser p(src);
 auto ast = p.parseExpression();
 if (p.cur.type != TokenType::End) throw CalcError(CalcErrorType::Syntax, errorMessage(CalcErrorType::Syntax), p.cur.pos);
 return ast->eval(cfg, hist, base);
}

/* ============================
   表示
   ============================ */

std::string formatResult(double v, const SystemConfig &cfg) {
 std::ostringstream oss;
 if (cfg.fixedDigits >= 0) {
  oss << std::fixed << std::setprecision(cfg.fixedDigits) << v;
 } else {
  oss << std::setprecision(15) << v;
 }
 return oss.str();
}

/* ============================
   関数定義
   ============================ */

void initFunctions(SystemConfig &cfg) {

 const double deg2rad = constants.at("Pi") / 180.0;
 const double rad2deg = 180.0 / constants.at("Pi");

 // ---- basic math ----
 cfg.functions["abs"] = {1, 1, [](auto &v, auto &) { return std::fabs(v[0]); }};
 cfg.functions["sqrt"] = {1, 1, [](auto &v, auto &) {
                           if (v[0] < 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), 0);
                           return std::sqrt(v[0]);
                          }};
 cfg.functions["floor"] = {1, 1, [](auto &v, auto &) { return std::floor(v[0]); }};
 cfg.functions["ceil"] = {1, 1, [](auto &v, auto &) { return std::ceil(v[0]); }};
 cfg.functions["round"] = {1, 1, [](auto &v, auto &) { return std::round(v[0]); }};
 cfg.functions["pow"] = {2, 2, [](auto &v, auto &) { return std::pow(v[0], v[1]); }};
 cfg.functions["mod"] = {2, 2, [](auto &v, auto &ctx) {
                          if (v[1] == 0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), ctx.pos);
                          return std::fmod(v[0], v[1]);
                         }};

 // ---- exp / log ----
 cfg.functions["exp"] = {1, 1, [](auto &v, auto &) { return std::exp(v[0]); }};
 cfg.functions["log"] = {1, 1, [](auto &v, auto &) {
                          if (v[0] <= 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), 0);
                          return std::log(v[0]);
                         }};

 // ---- trig (degree based) ----
 cfg.functions["sin"] = {1, 1, [=](auto &v, auto &) { return std::sin(v[0] * deg2rad); }};
 cfg.functions["cos"] = {1, 1, [=](auto &v, auto &) { return std::cos(v[0] * deg2rad); }};
 cfg.functions["tan"] = {1, 1, [=](auto &v, auto &) { return std::tan(v[0] * deg2rad); }};
 cfg.functions["asin"] = {1, 1, [=](auto &v, auto &) { return std::asin(v[0]) * rad2deg; }};
 cfg.functions["acos"] = {1, 1, [=](auto &v, auto &) { return std::acos(v[0]) * rad2deg; }};
 cfg.functions["atan"] = {1, 1, [=](auto &v, auto &) { return std::atan(v[0]) * rad2deg; }}; // ---- hyperbolic ----
 cfg.functions["sinh"] = {1, 1, [](auto &v, auto &) { return std::sinh(v[0]); }};
 cfg.functions["cosh"] = {1, 1, [](auto &v, auto &) { return std::cosh(v[0]); }};
 cfg.functions["tanh"] = {1, 1, [](auto &v, auto &) { return std::tanh(v[0]); }};

 // ---- inverse hyperbolic ----
 cfg.functions["asinh"] = {1, 1, [](auto &v, auto &) { return std::asinh(v[0]); }};
 cfg.functions["acosh"] = {1, 1, [](auto &v, auto &) {
                            if (v[0] < 1.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), 0);
                            return std::acosh(v[0]);
                           }};
 cfg.functions["atanh"] = {1, 1, [](auto &v, auto &) {
                            if (std::fabs(v[0]) >= 1.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), 0);
                            return std::atanh(v[0]);
                           }};

 // ---- misc ----
 cfg.functions["hypot"] = {2, 2, [](auto &v, auto &) { return std::hypot(v[0], v[1]); }};
 cfg.functions["cbrt"] = {1, 1, [](auto &v, auto &) { return std::cbrt(v[0]); }};
 cfg.functions["sign"] = {1, 1, [](auto &v, auto &) {
                           if (v[0] > 0) return 1.0;
                           if (v[0] < 0) return -1.0;
                           return 0.0;
                          }};

 // ---- angle conversion ----
 cfg.functions["torad"] = {1, 1, [=](auto &v, auto &) { return v[0] * deg2rad; }};
 cfg.functions["todeg"] = {1, 1, [=](auto &v, auto &) { return v[0] * rad2deg; }};

 // ---- number theory / statistics ----
 cfg.functions["gcd"] = {2, 2, [](auto &v, auto &ctx) {
                          int a = (int)requireInt(v[0], ctx.pos);
                          int b = (int)requireInt(v[1], ctx.pos);
                          return (double)gcdLL(a, b);
                         }};
 cfg.functions["comb"] = {2, 2, [](auto &v, auto &ctx) {
                           long long n = requireInt(v[0], ctx.pos);
                           long long r = requireInt(v[1], ctx.pos);
                           return (double)combLL(n, r);
                          }};
 cfg.functions["perm"] = {2, 2, [](auto &v, auto &ctx) {
                           long long n = requireInt(v[0], ctx.pos);
                           long long r = requireInt(v[1], ctx.pos);
                           return (double)permLL(n, r);
                          }};
 cfg.functions["lcm"] = {2, 2, [](auto &v, auto &ctx) {
                          long long a = requireInt(v[0], ctx.pos);
                          long long b = requireInt(v[1], ctx.pos);
                          return (double)(a / gcdLL(a, b) * b);
                         }};
 cfg.functions["sum"] = {1, -1, [](auto &v, auto &) { return std::accumulate(v.begin(), v.end(), 0.0); }};
 cfg.functions["mean"] = {1, -1, [](auto &v, auto &) { return std::accumulate(v.begin(), v.end(), 0.0) / v.size(); }};
 cfg.functions["min"] = {1, -1, [](auto &v, auto &) { return *std::min_element(v.begin(), v.end()); }};
 cfg.functions["max"] = {1, -1, [](auto &v, auto &) { return *std::max_element(v.begin(), v.end()); }};
 cfg.functions["gamma"] = {1, 1, [](auto &v, auto &ctx) { return std::tgamma(v[0]); }};
 cfg.functions["erf"] = {1, 1, [](auto &v, auto &ctx) { return std::erf(v[0]); }};

 // ---- programmierfeld ----

 cfg.functions["if"] = {3, 3, [](auto &v, auto &) { return v[0] != 0.0 ? v[1] : v[2]; }};
 cfg.functions["rand"] = {0, 0, [](auto &, auto &) { return (double)std::rand() / RAND_MAX; }};

 // ---- history / state ----
 cfg.functions["In"] = {1, 1, [](auto &v, auto &ctx) {
                         int i = (int)requireInt(v[0], ctx.pos);
                         if (i <= 0 || i > (int)ctx.hist.size()) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), ctx.pos);
                         return evaluate(ctx.hist[i - 1].expr, ctx.cfg, ctx.hist, ctx.base);
                        }};
 cfg.functions["Out"] = {1, 1, [](auto &v, auto &ctx) {
                          int i = (int)requireInt(v[0], ctx.pos);
                          if (i <= 0 || i > (int)ctx.hist.size()) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), ctx.pos);
                          return ctx.hist[i - 1].value;
                         }};
 cfg.functions["Prev"] = {1, 1, [](auto &v, auto &ctx) {
                           int k = (int)requireInt(v[0], ctx.pos);
                           int idx = ctx.base - k + 1;
                           if (idx <= 0 || idx > (int)ctx.hist.size()) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), ctx.pos);
                           return ctx.hist[idx - 1].value;
                          }};
 cfg.functions["Clear"] = {0, 0, [](auto &, auto &ctx) {
                            ctx.cfg.shouldClear = true;
                            return 0.0;
                           }};
 cfg.functions["Exit"] = {0, 0, [](auto &, auto &ctx) {
                           ctx.cfg.shouldExit = true;
                           return 0.0;
                          }};
}

/* ============================
  main
  ============================ */

int main() {

 SystemConfig cfg;
 std::vector<InputEntry> history;
 initFunctions(cfg);

 std::string line;
 std::cout << "================================\n  mm Calculator\n        mmKreutzef 2026\n================================\n";

 while (true) {
  std::cout << "\nIn [" << history.size() + 1 << "] := ";
  if (!std::getline(std::cin, line)) break;

  // 空行は戻りなさい
  if (line.find_first_not_of(" \t\r\n") == std::string::npos) { continue; }

  if (line.starts_with("SetFix")) {
   cfg.fixedDigits = std::stoi(line.substr(6));
   continue;
  }

  try {
   double v = evaluate(line, cfg, history, history.size());
   double v_out = round12(v);
   if (cfg.shouldClear) {
    history.clear();
    cfg.shouldClear = false;
    continue;
   }
   history.push_back({line, v_out});
   if (cfg.shouldExit) break;
   std::cout << "Out[" << history.size() << "] := " << formatResult(v_out, cfg) << "\n\n";
  } catch (const CalcError &e) {
   std::cout << "Error: " << e.what() << "\n" << line << "\n";
   std::cout << std::string(e.pos, ' ') << "^\n\n";
   history.push_back({line, 0.0});
  } catch (const ExitRequest &) { break; }
 }
 std::cout << "bye...nara";
}
