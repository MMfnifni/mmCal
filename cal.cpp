#include <cassert>
#include <cctype>
#include <cmath>
#include <complex>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

/*
Expr        ::= Add
Add         ::= Mul (("+" | "-") Mul)*
Mul         ::= Unary (("*" | "/" | ImplicitMul) Unary)*
Unary       ::= ("+" | "-") Unary
              | Power
Power       ::= Postfix ("^" Power)?
Postfix     ::= Primary ("!")*
Primary     ::= Number
              | "(" Expr ")"
              | Identifier
*/

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
 TypeError,
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
  case CalcErrorType::InvalidOperation: return "invalid operation";
  case CalcErrorType::TypeError: return "double or complex type error";
 }
 return "error";
}

struct CalcError : std::runtime_error {
  CalcErrorType type;
  size_t pos;
  CalcError(CalcErrorType t, const std::string &m, size_t p) : std::runtime_error(m), type(t), pos(p) {}
};

/* ============================
   基本定義
   ============================ */

inline constexpr double PI = 3.14159265358979323846264338327950288;

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
enum class CmpOp { Less, LessEq, Greater, GreaterEq, Equal, NotEqual };
CmpOp tokenToCmpOp(TokenType t) {
 switch (t) {
  case TokenType::Less: return CmpOp::Less;
  case TokenType::LessEq: return CmpOp::LessEq;
  case TokenType::Greater: return CmpOp::Greater;
  case TokenType::GreaterEq: return CmpOp::GreaterEq;
  case TokenType::Equal: return CmpOp::Equal;
  case TokenType::NotEqual: return CmpOp::NotEqual;
  default: throw std::logic_error("not a compare token");
 }
}

struct ExitRequest {};
struct ClearRequest {};
struct InvalidValue {};
struct SystemConfig;
struct FunctionContext;
struct FunctionDef;
struct InputEntry;

using Value = std::variant<double, std::complex<double>, InvalidValue>;
using FuncImpl = std::function<Value(const std::vector<Value> &, FunctionContext &)>;
using Complex = std::complex<double>;

enum class EvalKind { None, Value, Clear, Exit };

struct EvalResult {
  EvalKind kind = EvalKind::None;
  Value value{};
};
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
struct RuntimeState {
  bool shouldExit = false;
  bool shouldClear = false;
};
const int cnst_precision = 12;
struct SystemConfig {
  int precision = cnst_precision; // 表示精度（有効桁）
  bool trimTrailingZeros = true;  // 末尾の 0 を削るか
  std::unordered_map<std::string, FunctionDef> functions;
};
struct FunctionContext {
  SystemConfig &cfg;
  const std::vector<InputEntry> &hist;
  int base;
  size_t pos;
};

template <class> inline constexpr bool always_false = false;

bool eq(double a, double b) { return std::abs(a - b) < std::pow(10.0, -cnst_precision); }
inline bool isInvalid(const Value &v) { return std::holds_alternative<InvalidValue>(v); }
inline Complex toComplex(const Value &v) {
 if (std::holds_alternative<double>(v)) return Complex(std::get<double>(v), 0.0);
 return std::get<Complex>(v);
}
inline double inf(int sign = +1) { return sign >= 0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity(); }
inline double asReal(const Value &v, size_t pos) {
 if (!std::holds_alternative<double>(v)) throw CalcError(CalcErrorType::TypeError, errorMessage(CalcErrorType::TypeError), pos);
 return std::get<double>(v);
}

inline Complex asComplex(const Value &v) {
 if (std::holds_alternative<Complex>(v)) return std::get<Complex>(v);
 return Complex(std::get<double>(v), 0.0);
}

inline double round12(double v) {
 constexpr double scale = 1e12;
 return std::round(v * scale) / scale;
}

inline bool isReal(const Value &v) { return std::holds_alternative<double>(v); }

inline bool isComplex(const Value &v) { return std::holds_alternative<std::complex<double>>(v); }
inline bool isDouble(const Value &v) { return std::holds_alternative<double>(v); }

inline double asDouble(const Value &v, size_t pos = 0) {
 if (!isDouble(v)) throw CalcError(CalcErrorType::TypeError, errorMessage(CalcErrorType::TypeError), pos);
 return std::get<double>(v);
}

Value div(const Value &a, const Value &b, size_t pos) {
 if (isDouble(a) && isDouble(b)) {
  double r = asDouble(b, pos);
  if (r == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
  return asDouble(a) / r;
 }
 return asComplex(a) / asComplex(b);
}

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

Value roundValue(const Value &v, const SystemConfig &cfg) {
 if (isInvalid(v)) return v;

 auto roundDouble = [&](double d) {
  double scale = std::pow(10.0, cfg.precision);
  return std::round(d * scale) / scale;
 };

 if (std::holds_alternative<double>(v)) { return roundDouble(std::get<double>(v)); }

 if (std::holds_alternative<std::complex<double>>(v)) {
  auto c = std::get<std::complex<double>>(v);
  return std::complex<double>(roundDouble(c.real()), roundDouble(c.imag()));
 }

 return InvalidValue{};
}

/* ============================
   Token / Lexer
   ============================ */

struct Token {
  TokenType type = TokenType::End;
  std::string text;
  size_t pos = 0;
  Token() = default;

  Token(TokenType t, std::string s, size_t p) : type(t), text(std::move(s)), pos(p) {}
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

   if (p >= src.size()) return {TokenType::End, "", start};

   char c = src[p++];

   // number
   if (std::isdigit(c) || c == '.') {
    --p;
    char *end;
    double v = std::strtod(src.data() + p, &end);
    if (end == src.data() + p) { throw CalcError(CalcErrorType::InvalidNumber, "invalid number", p); }
    p = end - src.data();
    std::string s(src.data() + start, p - start);
    return {TokenType::Number, s, start};
   }

   // identifier
   if (std::isalpha(c) || c == '_') {
    std::string s(1, c);
    while (p < src.size() && (std::isalnum(src[p]) || src[p] == '_'))
     s += src[p++];
    return {TokenType::Identifier, s, start};
   }

   if (c == '<') {
    if (p < src.size() && src[p] == '=') {
     ++p;
     return {TokenType::LessEq, "<=", start};
    }
    return {TokenType::Less, "<", start};
   }

   if (c == '>') {
    if (p < src.size() && src[p] == '=') {
     ++p;
     return {TokenType::GreaterEq, ">=", start};
    }
    return {TokenType::Greater, ">", start};
   }

   if (c == '=') {
    if (p < src.size() && src[p] == '=') {
     ++p;
     return {TokenType::Equal, "==", start};
    }
    throw CalcError(CalcErrorType::InvalidCharacter, "single '=' not allowed", start);
   }

   switch (c) {
    case '+': return {TokenType::Plus, "+", start};
    case '-': return {TokenType::Minus, "-", start};
    case '*': return {TokenType::Mul, "*", start};
    case '/': return {TokenType::Div, "/", start};
    case '^': return {TokenType::Pow, "^", start};
    case '(': return {TokenType::LParen, "(", start};
    case ')': return {TokenType::RParen, ")", start};
    case '[': return {TokenType::LBracket, "[", start};
    case ']': return {TokenType::RBracket, "]", start};
    case ',': return {TokenType::Comma, ",", start};
    case '!': return {TokenType::Bang, "!", start};
    case '%': return {TokenType::Percent, "%", start};
   }

   throw CalcError(CalcErrorType::InvalidCharacter, std::string("unexpected character: ") + c, start);
  }
};

bool isValueEnd(TokenType t) {
 switch (t) {
  case TokenType::Number:
  case TokenType::Identifier:
  case TokenType::RParen:
  case TokenType::RBracket: return true;
  default: return false;
 }
}
bool isValueStart(TokenType t) {
 switch (t) {
  case TokenType::Number:
  case TokenType::Identifier:
  case TokenType::LParen:
  case TokenType::LBracket: return true;
  default: return false;
 }
}

bool isImplicitMul(const Token &a, const Token &b) {
 if (!isValueEnd(a.type)) return false;
 if (!isValueStart(b.type)) return false;

 // 例外: 関数呼び出し foo(...)
 if (a.type == TokenType::Identifier && b.type == TokenType::LParen) { return false; }

 return true;
}

/* ============================
   定数・補助
   ============================ */

const std::unordered_map<std::string, Value> constants = {
    {"Pi", 3.14159265358979323846},
    {"E", 2.71828182845904523536},
    {"Teu", 6.28318530717958647692},
    {"Phi", 1.61803398874989484820},
    {"NA", 6.02214076e23},
    {"I", Complex(0.0, 1.0)},
    {"ESP", std::pow(10.0, -cnst_precision)},
};

struct InputEntry {
  std::string expr; // In[n]
  Value value;      // Out[n]
};

inline void checkFinite(double v, size_t pos) {
 if (std::isnan(v)) throw CalcError(CalcErrorType::NaNResult, errorMessage(CalcErrorType::NaNResult), pos);
 if (std::isinf(v)) throw CalcError(CalcErrorType::InfiniteResult, errorMessage(CalcErrorType::InfiniteResult), pos);
}
inline void checkFinite(const Value &v, size_t pos) {
 if (std::holds_alternative<double>(v)) {
  double x = std::get<double>(v);
  if (std::isnan(x)) throw CalcError(CalcErrorType::NaNResult, errorMessage(CalcErrorType::NaNResult), pos);
  if (std::isinf(x)) throw CalcError(CalcErrorType::InfiniteResult, errorMessage(CalcErrorType::InfiniteResult), pos);
  return;
 }

 // complex
 const auto &z = std::get<std::complex<double>>(v);
 if (std::isnan(z.real()) || std::isnan(z.imag())) throw CalcError(CalcErrorType::NaNResult, errorMessage(CalcErrorType::NaNResult), pos);

 if (std::isinf(z.real()) || std::isinf(z.imag())) throw CalcError(CalcErrorType::InfiniteResult, errorMessage(CalcErrorType::InfiniteResult), pos);
}

inline Value makeValue(double v) { return Value{v}; }
inline Value makeValue(std::complex<double> v) { return Value{v}; }
double radToDeg(double r) { return r * 180.0 / PI; }
inline bool tolerantEqual(double a, double b, int precision) {
 double eps = std::pow(10.0, -precision);
 return std::abs(a - b) <= eps * std::max({1.0, std::abs(a), std::abs(b)});
}

inline long long requireInt(const Value &v, size_t pos) {
 double d = asDouble(v, pos);
 if (std::floor(d) != d) throw CalcError(CalcErrorType::NeedInteger, errorMessage(CalcErrorType::NeedInteger), pos);
 return static_cast<int64_t>(d);
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

Value evalCompare(const Value &lhs, const Value &rhs, CmpOp op, FunctionContext &ctx) {
 // 複素数は大小比較不可
 if (isComplex(lhs) || isComplex(rhs)) {
  if (op != CmpOp::Equal && op != CmpOp::NotEqual) throw CalcError(CalcErrorType::TypeError, "complex comparison", ctx.pos);
 }

 double a = asDouble(lhs, ctx.pos);
 double b = asDouble(rhs, ctx.pos);

 bool r = false;
 int prec = ctx.cfg.precision;

 switch (op) {
  case CmpOp::Equal: r = tolerantEqual(a, b, prec); break;

  case CmpOp::NotEqual: r = !tolerantEqual(a, b, prec); break;

  case CmpOp::Less: r = a < b && !tolerantEqual(a, b, prec); break;

  case CmpOp::LessEq: r = a < b || tolerantEqual(a, b, prec); break;

  case CmpOp::Greater: r = a > b && !tolerantEqual(a, b, prec); break;

  case CmpOp::GreaterEq: r = a > b || tolerantEqual(a, b, prec); break;
 }

 return r ? 1.0 : 0.0;
}

Value evaluate(const std::string &src, SystemConfig &cfg, const std::vector<InputEntry> &hist, int base);

/* ============================
   AST
   ============================ */

struct ASTNode {
  size_t pos = 0;
  virtual ~ASTNode() = default;
  virtual Value eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const = 0;
};

/* ---------- Number ---------- */

struct NumberNode : ASTNode {
  Value value;

  NumberNode(Value v, size_t p) : value(std::move(v)) { pos = p; }
  Value eval(SystemConfig &, const std::vector<InputEntry> &, int) const override { return value; }
};

/* ---------- Unary ---------- */

enum class UnaryOp { Plus, Minus };

struct UnaryNode : ASTNode {
  UnaryOp op;
  std::unique_ptr<ASTNode> rhs;

  UnaryNode(UnaryOp o, std::unique_ptr<ASTNode> r, size_t p) : op(o), rhs(std::move(r)) { pos = p; }

  Value eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
   auto v = rhs->eval(cfg, hist, base);
   if (op == UnaryOp::Minus) {
    if (isComplex(v)) return -toComplex(v);
    return -asDouble(v, pos);
   }
   return v;
  }
};

struct PostfixUnaryNode : ASTNode {
  char op;
  std::unique_ptr<ASTNode> expr;

  PostfixUnaryNode(char op, std::unique_ptr<ASTNode> e, size_t p) : op(op), expr(std::move(e)) { pos = p; }

  Value eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {

   Value v = expr->eval(cfg, hist, base);

   // 後置演算子は実数のみ対応
   if (!std::holds_alternative<double>(v)) { throw CalcError(CalcErrorType::TypeError, errorMessage(CalcErrorType::TypeError), pos); }

   double d = std::get<double>(v);
   double r;

   switch (op) {
    case '!': r = factorial(d, pos); break;
    default: throw CalcError(CalcErrorType::Syntax, "unknown postfix operator", pos);
   }

   Value out = r;
   checkFinite(out, pos);
   return out;
  }
};

/* ---------- Compare ---------- */
struct CompareNode : ASTNode {
  CmpOp op;
  std::unique_ptr<ASTNode> lhs, rhs;

  CompareNode(CmpOp o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p) : op(o), lhs(std::move(l)), rhs(std::move(r)) { pos = p; }

  Value eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
   Value a = lhs->eval(cfg, hist, base);
   Value b = rhs->eval(cfg, hist, base);

   FunctionContext ctx{cfg, hist, base, pos};
   return evalCompare(a, b, op, ctx);
  }
};

/* ---------- Binary ---------- */
enum class BinOp { Add, Sub, Mul, Div, Pow };
inline BinOp tokenToBinOp(TokenType t) {
 switch (t) {
  case TokenType::Plus: return BinOp::Add;
  case TokenType::Minus: return BinOp::Sub;
  case TokenType::Mul: return BinOp::Mul;
  case TokenType::Div: return BinOp::Div;
  case TokenType::Pow: return BinOp::Pow;
  default: throw CalcError(CalcErrorType::InvalidOperation, errorMessage(CalcErrorType::InvalidOperation), 0);
 }
}

struct BinaryNode : ASTNode {
  BinOp op;
  std::unique_ptr<ASTNode> lhs, rhs;

  BinaryNode(BinOp o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p) : op(o), lhs(std::move(l)), rhs(std::move(r)) { pos = p; }

  Value eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {

   auto a = lhs->eval(cfg, hist, base);
   auto b = rhs->eval(cfg, hist, base);

   if (isComplex(a) || isComplex(b)) {
    Complex x = toComplex(a);
    Complex y = toComplex(b);
    switch (op) {
     case BinOp::Add: return x + y;
     case BinOp::Sub: return x - y;
     case BinOp::Mul: return x * y;
     case BinOp::Div: return x / y;
     case BinOp::Pow: return std::pow(x, y);
    }
   }

   double x = asDouble(a, pos);
   double y = asDouble(b, pos);

   switch (op) {
    case BinOp::Add: return x + y;
    case BinOp::Sub: return x - y;
    case BinOp::Mul: return x * y;
    case BinOp::Div:
     if (y == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
     return x / y;
    case BinOp::Pow: return std::pow(x, y);
   }

   throw CalcError(CalcErrorType::InvalidOperation, "invalid op", pos);
  }
};

/* ---------- Function ---------- */

struct FunctionCallNode : ASTNode {
  std::string name;
  std::vector<std::unique_ptr<ASTNode>> args;

  Value eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {

   auto it = cfg.functions.find(name);
   if (it == cfg.functions.end()) throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), pos);

   auto &f = it->second;
   int argc = static_cast<int>(args.size());

   if (!f.validArgc(argc)) throw CalcError(CalcErrorType::FunctionMissing, errorMessage(CalcErrorType::FunctionMissing), pos);

   std::vector<Value> v;
   v.reserve(args.size());
   for (auto &a : args)
    v.push_back(a->eval(cfg, hist, base));

   FunctionContext ctx{cfg, hist, base, pos};

   Value r = f.f(v, ctx);

   // double のときだけ有限チェック
   if (std::holds_alternative<double>(r)) checkFinite(std::get<double>(r), pos);

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
  OutNode(int idx, size_t p) : index(idx) { pos = p; }

  Value eval(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
   if (index == 0) { throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos); }

   int i;
   if (index > 0) {
    i = index;
   } else {
    i = static_cast<int>(hist.size()) + index + 1; // -1 → last
   }

   if (i <= 0 || i > static_cast<int>(hist.size())) { throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos); }

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
  Token prev;
  SystemConfig &cfg;

  Parser(SystemConfig &cfg, const std::string &s) : cfg(cfg), lex(s) { cur = lex.get(); }

  void advance() {
   prev = cur;
   cur = lex.get();
  }

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

  bool isFunctionName(const std::string &name) const { return cfg.functions.find(name) != cfg.functions.end(); }
  bool isConstantName(const std::string &name) const { return constants.find(name) != constants.end(); }

  bool allowImplicitMul(const Token &prev, const Token &cur, const SystemConfig &cfg) {
   // 右が値を開始しないならダメ
   if (!isValueStart(cur.type)) return false;

   // 左が値を完結していないならダメ
   if (!isValueEnd(prev.type)) return false;

   // Identifier + Number
   if (prev.type == TokenType::Identifier && cur.type == TokenType::Number) {
    return isConstantName(prev.text); // ← sin30 を永久封印
   }

   // Identifier + '('
   if (prev.type == TokenType::Identifier && cur.type == TokenType::LParen) {
    return !isFunctionName(prev.text); // Pi(Pi) OK, sin(30) は関数呼び出し
   }

   return true;
  }

  bool isImplicitMul(const Token &prev, const Token &cur) const {
   // 右が「値の開始」でなければダメ
   auto isValueStart = [&](TokenType t) { return t == TokenType::Number || t == TokenType::Identifier || t == TokenType::LParen; };

   // 左が「値の終了」でなければダメ
   auto isValueEnd = [&](TokenType t) { return t == TokenType::Number || t == TokenType::Identifier || t == TokenType::RParen || t == TokenType::Bang; };

   if (!isValueEnd(prev.type) || !isValueStart(cur.type)) return false;

   // Identifier + Number → 定数のみ許可（sin30 永久封印）
   if (prev.type == TokenType::Identifier && cur.type == TokenType::Number) { return isConstantName(prev.text); }

   // Identifier + '('
   if (prev.type == TokenType::Identifier && cur.type == TokenType::LParen) {
    // 関数なら implicit mul しない（= 関数呼び出しに任せる）
    return !isFunctionName(prev.text);
   }

   return true;
  }

  std::unique_ptr<ASTNode> parseExpression();
  std::unique_ptr<ASTNode> parseCompare();

 private:
  std::unique_ptr<ASTNode> parseTerm();
  std::unique_ptr<ASTNode> parsePower();
  std::unique_ptr<ASTNode> parseUnary();
  std::unique_ptr<ASTNode> parsePrimary();
  std::unique_ptr<ASTNode> parsePostfix();
  // bool acceptParen() {
  //  if (cur.type == TokenType::LParen || cur.type == TokenType::LBracket) {
  //   advance();
  //   return true;
  //  }
  //  return false;
  // }

  // void expectParenClose(TokenType open) {
  //  if (open == TokenType::LParen) expect(TokenType::RParen);
  //  else expect(TokenType::RBracket);
  // }
};
std::unique_ptr<ASTNode> Parser::parseCompare() {
 auto n = parseExpression();

 while (cur.type == TokenType::Less || cur.type == TokenType::LessEq || cur.type == TokenType::Greater || cur.type == TokenType::GreaterEq || cur.type == TokenType::Equal || cur.type == TokenType::NotEqual) {

  CmpOp op = tokenToCmpOp(cur.type);
  size_t p = cur.pos;
  advance();
  auto r = parseExpression();
  n = std::make_unique<CompareNode>(op, std::move(n), std::move(r), p);
 }

 return n;
}

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
  BinOp op = tokenToBinOp(cur.type);
  size_t p = cur.pos;
  advance();
  auto r = parseTerm();
  n = std::make_unique<BinaryNode>(op, std::move(n), std::move(r), p);
 }
 return n;
}

std::unique_ptr<ASTNode> Parser::parseTerm() {
 auto n = parseUnary();

 while (true) {
  if (cur.type == TokenType::Mul || cur.type == TokenType::Div) {
   BinOp op = tokenToBinOp(cur.type);
   size_t p = cur.pos;
   advance();
   auto r = parseUnary();
   n = std::make_unique<BinaryNode>(op, std::move(n), std::move(r), p);
  }
  // 暗黙乗算
  else if (isImplicitMul(prev, cur)) {
   size_t p = cur.pos;
   auto r = parseUnary();
   n = std::make_unique<BinaryNode>(BinOp::Mul, std::move(n), std::move(r), p);
  } else {
   break;
  }
 }
 return n;
}

std::unique_ptr<ASTNode> Parser::parsePower() {
 auto n = parsePostfix();
 if (cur.type == TokenType::Pow) {
  size_t p = cur.pos;
  advance();
  auto r = parseUnary();
  n = std::make_unique<BinaryNode>(BinOp::Pow, std::move(n), std::move(r), p);
 }
 return n;
}

std::unique_ptr<ASTNode> Parser::parseUnary() {
 if (cur.type == TokenType::Plus || cur.type == TokenType::Minus) {
  UnaryOp op = (cur.type == TokenType::Minus) ? UnaryOp::Minus : UnaryOp::Plus;
  size_t p = cur.pos;
  advance();
  return std::make_unique<UnaryNode>(op, parseUnary(), p);
 }
 return parsePower();
}

std::unique_ptr<ASTNode> Parser::parsePrimary() {
 size_t p = cur.pos;

 // ---- number ----
 if (cur.type == TokenType::Number) {
  double v;
  try {
   v = std::stod(cur.text);
  } catch (...) { throw CalcError(CalcErrorType::InvalidNumber, "invalid number", p); }
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

  // ----- followed by ( or [ -----
  if (cur.type == TokenType::LParen || cur.type == TokenType::LBracket) {
   TokenType open = cur.type;
   size_t callPos = cur.pos;
   advance();

   // 関数なら普通に call
   if (cfg.functions.count(name)) {
    auto f = std::make_unique<FunctionCallNode>();
    f->name = name;
    f->pos = p;

    // no args
    if ((open == TokenType::LParen && accept(TokenType::RParen)) || (open == TokenType::LBracket && accept(TokenType::RBracket))) { return f; }

    while (true) {
     f->args.push_back(parseExpression());
     if (accept(TokenType::Comma)) continue;
     if (open == TokenType::LParen) expect(TokenType::RParen);
     else expect(TokenType::RBracket);
     break;
    }
    return f;
   }

   // ★ 定数なら「定数 × (expr)」
   if (constants.count(name)) {
    auto lhs = std::make_unique<NumberNode>(constants.at(name), p);
    auto rhs = parseExpression();
    if (open == TokenType::LParen) expect(TokenType::RParen);
    else expect(TokenType::RBracket);
    return std::make_unique<BinaryNode>(BinOp::Mul, std::move(lhs), std::move(rhs), callPos);
   }

   throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), p);
  }

  // ----- identifier only -----
  if (constants.count(name)) return std::make_unique<NumberNode>(constants.at(name), p);

  throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), p);
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

Value evaluate(const std::string &src, SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) {
 Parser p(cfg, src);
 auto ast = p.parseCompare();
 if (p.cur.type != TokenType::End) throw CalcError(CalcErrorType::Syntax, errorMessage(CalcErrorType::Syntax), p.cur.pos);
 return ast->eval(cfg, hist, base);
}

/* ============================
   表示
   ============================ */
inline void setupStream(std::ostringstream &oss, const SystemConfig &cfg) {
 if (cfg.precision >= 0) {
  oss << std::fixed << std::setprecision(cfg.precision);
 } else {
  oss << std::setprecision(15);
 }
}

// std::string formatResult(double v, const SystemConfig &cfg) {
//  std::ostringstream oss;
//  if (cfg.fixedDigits >= 0) {
//   oss << std::fixed << std::setprecision(cfg.fixedDigits) << v;
//  } else {
//   oss << std::setprecision(15) << v;
//  }
//  return oss.str();
// }
// std::string formatResult(const std::complex<double> &z, const SystemConfig &cfg) {
//  std::ostringstream oss;
//  setupStream(oss, cfg);
//
//  oss << z.real();
//
//  if (z.imag() >= 0) oss << "+";
//  oss << z.imag() << "I";
//
//  return oss.str();
// }
std::string formatReal(double x, const SystemConfig &cfg) {
 std::ostringstream oss;

 if (cfg.trimTrailingZeros) {
  oss << std::fixed << std::setprecision(cfg.precision);
 } else {
  oss << std::setprecision(cfg.precision);
 }

 oss << x;
 std::string s = oss.str();

 // fixed + trimTrailingZeros のとき末尾処理
 if (cfg.trimTrailingZeros && s.find('.') != std::string::npos) {
  while (!s.empty() && s.back() == '0')
   s.pop_back();
  if (!s.empty() && s.back() == '.') s.pop_back();
 }

 return s;
}

std::string formatComplex(const std::complex<double> &c, const SystemConfig &cfg) {
 double re = c.real();
 double im = c.imag();

 // 純実数
 if (std::abs(im) < std::pow(10.0, -cnst_precision)) return formatReal(re, cfg);

 std::ostringstream oss;

 if (re != 0.0) oss << formatReal(re, cfg);

 if (im > 0 && re != 0.0) oss << "+";

 if (im == 1.0) oss << "I";
 else if (im == -1.0) oss << "-I";
 else oss << formatReal(im, cfg) << "I";

 return oss.str();
}

// std::string formatResult(const Value &v, const SystemConfig &cfg) {
//  if (isReal(v)) return formatReal(std::get<double>(v), cfg);
//  return formatComplex(std::get<std::complex<double>>(v), cfg);
// }

std::string formatResult(const Value &v, const SystemConfig &cfg) {
 return std::visit(
     [&](auto &&x) -> std::string {
      using T = std::decay_t<decltype(x)>;

      if constexpr (std::is_same_v<T, double>) {
       if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
       if (std::isnan(x)) return "nan";
       return formatReal(x, cfg);

      } else if constexpr (std::is_same_v<T, std::complex<double>>) {
       double re = x.real(), im = x.imag();

       if (std::isinf(re) || std::isinf(im)) {
        std::string s;
        if (std::isinf(re)) s += (re > 0 ? "inf" : "-inf");
        if (std::isinf(im)) {
         if (!s.empty() && im > 0) s += "+";
         s += (im > 0 ? "inf" : "-inf");
         s += "i";
        }
        return s;
       }

       return formatComplex(x, cfg);

      } else if constexpr (std::is_same_v<T, InvalidValue>) {
       return "Invalid";

      } else {
       static_assert(always_false<T>, "Unhandled Value type");
      }
     },
     v);
}

inline Value mulReal(const Value &v, double k, size_t pos) {
 if (isComplex(v)) return asComplex(v) * k;
 return asDouble(v, pos) * k;
}

inline Value divReal(const Value &v, double k, size_t pos) {
 if (k == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
 if (isComplex(v)) return asComplex(v) / k;
 return asDouble(v, pos) / k;
}

/* ============================
   関数定義
   ============================ */
Value deg2rad_v(const Value &v, size_t pos) {
 constexpr double k = PI / 180.0;
 return mulReal(v, k, pos);
}

Value rad2deg_v(const Value &v, size_t pos) {
 constexpr double k = 180.0 / PI;
 return mulReal(v, k, pos);
}

Value deg2grad_v(const Value &v, size_t pos) { return mulReal(v, 400.0 / 360.0, pos); }

Value grad2deg_v(const Value &v, size_t pos) { return mulReal(v, 360.0 / 400.0, pos); }

Value rad2grad_v(const Value &v, size_t pos) {
 constexpr double k = 200.0 / PI;
 return mulReal(v, k, pos);
}

Value grad2rad_v(const Value &v, size_t pos) {
 constexpr double k = PI / 200.0;
 return mulReal(v, k, pos);
}
void initFunctions(SystemConfig &cfg) {

 const double deg2rad = PI / 180.0;
 const double rad2deg = 180.0 / PI;

 // ---- basic math ----
 cfg.functions["abs"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                          if (isComplex(v[0])) return std::abs(asComplex(v[0]));
                          return std::fabs(asReal(v[0], ctx.pos));
                         }};

 cfg.functions["sqrt"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           if (isComplex(v[0])) return std::sqrt(asComplex(v[0]));
                           double x = asReal(v[0], ctx.pos);
                           if (x < 0) return Complex(0, std::sqrt(-x));
                           return std::sqrt(x);
                          }};

 cfg.functions["floor"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::floor(asReal(v[0], ctx.pos)); }};
 cfg.functions["ceil"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::ceil(asReal(v[0], ctx.pos)); }};
 cfg.functions["round"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::round(asReal(v[0], ctx.pos)); }};

 cfg.functions["pow"] = {2, 2, [](auto &v, auto &ctx) -> Value { return std::pow(asComplex(v[0]), asComplex(v[1])); }};

 // ---- exp / log ----
 cfg.functions["exp"] = {1, 1, [](auto &v, auto &) -> Value { return std::exp(asComplex(v[0])); }};

 cfg.functions["log"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                          if (isComplex(v[0])) return std::log(asComplex(v[0]));
                          double x = asReal(v[0], ctx.pos);
                          if (x <= 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                          return std::log(x);
                         }};
 cfg.functions["log"] = {1, 2, [](auto &v, auto &ctx) -> Value {
                          if (isComplex(v[0]) || (v.size() == 2 && isComplex(v[1]))) {
                           if (v.size() == 1) {
                            return std::log(asComplex(v[0]));
                           } else {
                            auto base = asComplex(v[0]);
                            auto x = asComplex(v[1]);
                            return std::log(x) / std::log(base);
                           }
                          }
                          if (v.size() == 1) {
                           double x = asReal(v[0], ctx.pos);
                           if (x <= 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                           return std::log(x);
                          }
                          // v.size() == 2
                          double base = asReal(v[0], ctx.pos);
                          double x = asReal(v[1], ctx.pos);

                          if (base <= 0 || base == 1 || x <= 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                          return std::log(x) / std::log(base);
                         }};

 // ---- trig (degree based) ----
 cfg.functions["sin"] = {1, 1, [=](auto &v, auto &) -> Value { return std::sin(asComplex(v[0]) * deg2rad); }};
 cfg.functions["cos"] = {1, 1, [=](auto &v, auto &) -> Value { return std::cos(asComplex(v[0]) * deg2rad); }};
 // cfg.functions["tan"] = {1, 1, [=](auto &v, auto &) -> Value { return std::tan(asComplex(v[0]) * deg2rad); }};
 cfg.functions["tan"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                          double x = asReal(v[0], ctx.pos);
                          double r = x * deg2rad;
                          double c = std::cos(r);
                          if (std::abs(c) < asDouble(constants.at("ESP"))) {
                           double s = std::sin(r);
                           return inf(s >= 0 ? +1 : -1);
                          }
                          return std::tan(r);
                         }};
 cfg.functions["cot"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                          double x = asReal(v[0], ctx.pos);
                          double r = x * deg2rad;

                          double s = std::sin(r);
                          if (std::abs(s) < asDouble(constants.at("ESP"))) {
                           double c = std::cos(r);
                           return inf(c >= 0 ? +1 : -1);
                          }
                          return std::cos(r) / s;
                         }};
 cfg.functions["sec"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                          double x = asReal(v[0], ctx.pos);
                          double r = x * deg2rad;

                          double c = std::cos(r);
                          if (std::abs(c) < asDouble(constants.at("ESP"))) {
                           double s = std::sin(r);
                           return inf(s >= 0 ? +1 : -1);
                          }
                          return 1.0 / c;
                         }};
 cfg.functions["csc"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                          double x = asReal(v[0], ctx.pos);
                          double r = x * deg2rad;

                          double s = std::sin(r);
                          if (std::abs(s) < asDouble(constants.at("ESP"))) {
                           double c = std::cos(r);
                           return inf(c >= 0 ? +1 : -1);
                          }
                          return 1.0 / s;
                         }};
 cfg.functions["asin"] = {1, 1, [=](auto &v, auto &) -> Value { return std::asin(asComplex(v[0])) * rad2deg; }};
 cfg.functions["acos"] = {1, 1, [=](auto &v, auto &) -> Value { return std::acos(asComplex(v[0])) * rad2deg; }};
 cfg.functions["atan"] = {1, 1, [=](auto &v, auto &) -> Value { return std::atan(asComplex(v[0])) * rad2deg; }};

 // ---- hyperbolic ----
 cfg.functions["sinh"] = {1, 1, [](auto &v, auto &) -> Value { return std::sinh(asComplex(v[0])); }};
 cfg.functions["cosh"] = {1, 1, [](auto &v, auto &) -> Value { return std::cosh(asComplex(v[0])); }};
 cfg.functions["tanh"] = {1, 1, [](auto &v, auto &) -> Value { return std::tanh(asComplex(v[0])); }};

 // ---- coonvert angle ----
 cfg.functions["DtoR"] = {1, 1, [](auto &v, auto &ctx) -> Value { return deg2rad_v(v[0], ctx.pos); }};
 cfg.functions["DtoG"] = {1, 1, [](auto &v, auto &ctx) -> Value { return deg2grad_v(v[0], ctx.pos); }};
 cfg.functions["RtoD"] = {1, 1, [](auto &v, auto &ctx) -> Value { return rad2deg_v(v[0], ctx.pos); }};
 cfg.functions["RtoG"] = {1, 1, [](auto &v, auto &ctx) -> Value { return rad2grad_v(v[0], ctx.pos); }};
 cfg.functions["GtoD"] = {1, 1, [](auto &v, auto &ctx) -> Value { return grad2deg_v(v[0], ctx.pos); }};
 cfg.functions["GtoR"] = {1, 1, [](auto &v, auto &ctx) -> Value { return grad2rad_v(v[0], ctx.pos); }};

 // ---- inverse hyperbolic ----
 cfg.functions["asinh"] = {1, 1, [](auto &v, auto &) -> Value { return std::asinh(asComplex(v[0])); }};
 cfg.functions["acosh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = asReal(v[0], ctx.pos);
                            if (x < 1.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                            return std::acosh(x);
                           }};
 cfg.functions["atanh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = asReal(v[0], ctx.pos);
                            if (std::fabs(x) >= 1.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                            return std::atanh(x);
                           }};

 // ---- misc ----
 cfg.functions["hypot"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            double x = asReal(v[0], ctx.pos);
                            double y = asReal(v[1], ctx.pos);
                            return std::hypot(x, y);
                           }};
 cfg.functions["cbrt"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::cbrt(asReal(v[0], ctx.pos)); }};
 cfg.functions["sign"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           double x = asReal(v[0], ctx.pos);
                           return (x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0;
                          }};

 // ---- angle conversion ----
 cfg.functions["torad"] = {1, 1, [=](auto &v, auto &) -> Value { return asComplex(v[0]) * deg2rad; }};
 cfg.functions["todeg"] = {1, 1, [=](auto &v, auto &) -> Value { return asComplex(v[0]) * rad2deg; }};

 // ---- number theory / statistics ----
 cfg.functions["gcd"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                          long long a = requireInt(v[0], ctx.pos);
                          long long b = requireInt(v[1], ctx.pos);
                          return (double)gcdLL(a, b);
                         }};
 cfg.functions["comb"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           long long n = requireInt(v[0], ctx.pos);
                           long long r = requireInt(v[1], ctx.pos);
                           return (double)combLL(n, r);
                          }};
 cfg.functions["perm"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           long long n = requireInt(v[0], ctx.pos);
                           long long r = requireInt(v[1], ctx.pos);
                           return (double)permLL(n, r);
                          }};
 cfg.functions["lcm"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                          long long a = requireInt(v[0], ctx.pos);
                          long long b = requireInt(v[1], ctx.pos);
                          return (double)(a / gcdLL(a, b) * b);
                         }};

 cfg.functions["sum"] = {1, -1, [](auto &v, auto &) -> Value {
                          Complex acc{0, 0};
                          for (auto &x : v)
                           acc += asComplex(x);
                          return acc.imag() == 0 ? Value(acc.real()) : Value(acc);
                         }};
 cfg.functions["mean"] = {1, -1, [](auto &v, auto &) -> Value {
                           Complex acc{0, 0};
                           for (auto &x : v)
                            acc += asComplex(x);
                           acc /= (double)v.size();
                           return acc.imag() == 0 ? Value(acc.real()) : Value(acc);
                          }};
 cfg.functions["mod"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                          Complex a = asComplex(v[0]), b = asComplex(v[1]);
                          if (a.imag() != 0 || b.imag() != 0) throw CalcError(CalcErrorType::DomainError, "mod: complex argument", ctx.pos);
                          double x = a.real(), y = b.real();
                          if (y == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), ctx.pos);
                          double r = x - y * std::floor(x / y);
                          return r;
                         }};
 cfg.functions["min"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                          double m = asReal(v[0], ctx.pos);
                          for (auto &x : v)
                           m = std::min(m, asReal(x, ctx.pos));
                          return m;
                         }};
 cfg.functions["max"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                          double m = asReal(v[0], ctx.pos);
                          for (auto &x : v)
                           m = std::max(m, asReal(x, ctx.pos));
                          return m;
                         }};
 cfg.functions["gamma"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = asReal(v[0], ctx.pos);
                            return std::tgamma(x);
                           }};
 cfg.functions["erf"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                          double x = asReal(v[0], ctx.pos);
                          return std::erf(x);
                         }};

 // ---- control ----
 cfg.functions["if"] = {3, 3, [](auto &v, auto &ctx) -> Value { return asReal(v[0], ctx.pos) != 0.0 ? v[1] : v[2]; }};
 cfg.functions["rand"] = {0, 0, [](auto &, auto &) -> Value { return (double)std::rand() / RAND_MAX; }};

 // ---- history / state ----
 cfg.functions["In"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                         int i = (int)requireInt(v[0], ctx.pos);
                         if (i <= 0 || i > (int)ctx.hist.size()) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), ctx.pos);
                         return evaluate(ctx.hist[i - 1].expr, ctx.cfg, ctx.hist, ctx.base);
                        }};
 cfg.functions["Out"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                          int i = (int)requireInt(v[0], ctx.pos);
                          if (i <= 0 || i > (int)ctx.hist.size()) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), ctx.pos);
                          return ctx.hist[i - 1].value;
                         }};
 cfg.functions["Prev"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           int k = (int)requireInt(v[0], ctx.pos);
                           int idx = ctx.base - k + 1;
                           if (idx <= 0 || idx > (int)ctx.hist.size()) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), ctx.pos);
                           return ctx.hist[idx - 1].value;
                          }};
 cfg.functions["Clear"] = {0, 0, [](auto &, auto &) -> Value { throw ClearRequest{}; }};
 cfg.functions["Exit"] = {0, 0, [](auto &, auto &) -> Value { throw ExitRequest{}; }};
}

/* ============================
  評価マン
  ============================ */

EvalResult evalLine(const std::string &line, SystemConfig &cfg, RuntimeState &rt, std::vector<InputEntry> &history) {
 EvalResult res{};

 try {
  Value v = evaluate(line, cfg, history, static_cast<int>(history.size()));
  res.kind = EvalKind::Value;
  res.value = roundValue(v, cfg);
  return res;
 } catch (const ClearRequest &) {
  rt.shouldClear = true;
  res.kind = EvalKind::Clear;
  return res;
 } catch (const ExitRequest &) {
  rt.shouldExit = true;
  res.kind = EvalKind::Exit;
  return res;
 }
}

/* ============================
  main
  ============================ */
int main(int argc, char *argv[]) {
 SystemConfig syscfg;
 RuntimeState rtmstt;
 std::vector<InputEntry> history;

 initFunctions(syscfg);

 // ================================
 // CLI モード
 // ================================
 if (argc >= 2) {
  std::string line = argv[1];

  try {
   EvalResult res = evalLine(line, syscfg, rtmstt, history);

   if (res.kind == EvalKind::Value) { std::cout << formatResult(res.value, syscfg) << "\n"; }
   // Clear / None / Exit は何も出さず終了
  } catch (const CalcError &e) {
   std::cout << "Error: " << e.what() << "\n" << line << "\n" << std::string(e.pos, ' ') << "^\n";
   return 1;
  }

  return 0;
 }

 // ================================
 // REPL モード
 // ================================
 std::cout << "================================\n"
              "  mm Calculator\n"
              "        mmKreutzef 2021-2026\n"
              "================================\n";

 std::string line;

 while (true) {
  std::cout << "In [" << history.size() + 1 << "] := ";
  if (!std::getline(std::cin, line)) break;

  if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

  // ---- メタコマンド ----
  if (line.starts_with("SetFix")) {
   syscfg.precision = std::stoi(line.substr(6));
   continue;
  }

  try {
   EvalResult res = evalLine(line, syscfg, rtmstt, history);

   switch (res.kind) {
    case EvalKind::Clear: history.clear(); break;

    case EvalKind::Exit: return 0;

    case EvalKind::Value:
     history.push_back({line, res.value});
     std::cout << "Out[" << history.size() << "] := " << formatResult(res.value, syscfg) << "\n\n";
     break;

    case EvalKind::None: break;
   }
  } catch (const CalcError &e) { std::cout << "Error: " << e.what() << "\n" << line << "\n" << std::string(e.pos, ' ') << "^\n\n"; }
 }
 std::cout << "bye...nara\n";
}