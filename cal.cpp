#include <algorithm>
#include <cassert>
#include <cctype>
#include <charconv>
#include <cmath>
#include <complex>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>
#ifdef _WIN32
# define DLL_EXPORT __declspec(dllexport)
#else
# define DLL_EXPORT
#endif

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

namespace mm::cal {

 /* ============================
    エラー
    ============================ */

 enum class CalcErrorType {
  NotImplemented,
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
  InvalidArgument,
  NaNDetected,
  InfinityDetected,
  InternalError,
 };

 constexpr const char *errorMessage(CalcErrorType t) {
  switch (t) {
   case CalcErrorType::NotImplemented: return "this feature is not yet implemented";
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
   case CalcErrorType::InvalidArgument: return "invalid argument count";
   case CalcErrorType::NaNDetected: return "NaN detected";
   case CalcErrorType::InfinityDetected: return "Infinity detected";
   case CalcErrorType::InternalError: return "InternalError";
  }
  return "unknown calculation error";
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
 // 以下のほうが高速になる。可読性とのトレードオフ
 // struct Value {
 //  enum { Real, Complex } type;
 //  double r;
 //  std::complex<double> c;
 //};
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

 inline void calcWarn(const FunctionContext &ctx, const std::string &msg) { std::cerr << "[WARN] pos=" << ctx.pos << " " << msg << "\n"; }
 inline void calcWarn(SystemConfig &cfg, size_t pos, const std::string &msg) { std::cerr << "[WARN] pos=" << pos << " " << msg << "\n"; }

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

 inline int64_t asInt64(const Value &v, size_t pos) {
  if (isComplex(v)) throw CalcError(CalcErrorType::TypeError, "int required", pos);
  double x = std::get<double>(v);
  if (!std::isfinite(x)) throw CalcError(CalcErrorType::DomainError, "nan/inf", pos);

  if (std::floor(x) != x) throw CalcError(CalcErrorType::TypeError, "integer required", pos);

  if (x < (double)LLONG_MIN || x > (double)LLONG_MAX) throw CalcError(CalcErrorType::OutOfRange, "int64 overflow", pos);

  return (int64_t)x;
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
 inline long long checkedMul(long long a, long long b, size_t pos) {
  if (a == 0 || b == 0) return 0;

  if (a > 0) {
   if (b > 0) {
    if (a > LLONG_MAX / b) throw CalcError(CalcErrorType::Overflow, errorMessage(CalcErrorType::Overflow), pos);
   } else {
    if (b < LLONG_MIN / a) throw CalcError(CalcErrorType::Overflow, errorMessage(CalcErrorType::Overflow), pos);
   }
  } else {
   if (b > 0) {
    if (a < LLONG_MIN / b) throw CalcError(CalcErrorType::Overflow, errorMessage(CalcErrorType::Overflow), pos);
   } else {
    if (a != 0 && b < LLONG_MAX / a) throw CalcError(CalcErrorType::Overflow, errorMessage(CalcErrorType::Overflow), pos);
   }
  }

  return a * b;
 }
 static long double factLD(long long n, int pos) {
  if (n < 0) throw CalcError(CalcErrorType::DomainError, "fact: n < 0", pos);
  long double acc = 1.0L;
  for (long long i = 2; i <= n; ++i) {
   acc *= (long double)i;
   if (!std::isfinite((double)acc)) throw CalcError(CalcErrorType::Overflow, "fact: overflow", pos);
  }
  return acc;
 }
 static unsigned long long fibULL(long long n, int pos) {
  if (n < 0) throw CalcError(CalcErrorType::DomainError, "fib: n < 0", pos);
  if (n == 0) return 0;
  if (n == 1) return 1;

  unsigned long long a = 0, b = 1;
  for (long long i = 2; i <= n; ++i) {
   unsigned long long c = a + b;
   if (c < b) throw CalcError(CalcErrorType::Overflow, "fib: overflow", pos);
   a = b;
   b = c;
  }
  return b;
 }
 static bool isPrimeLL(long long n) {
  if (n <= 1) return false;
  if (n <= 3) return true;
  if (n % 2 == 0) return false;
  if (n % 3 == 0) return false;

  for (long long i = 5; i <= n / i; i += 6) {
   if (n % i == 0) return false;
   if (n % (i + 2) == 0) return false;
  }
  return true;
 }
 static inline std::vector<double> collectReals(const std::vector<Value> &v, FunctionContext &ctx) {
  std::vector<double> a;
  a.reserve(v.size());
  for (auto &x : v)
   a.push_back(asReal(x, ctx.pos));
  return a;
 }
 static inline double medianOfSorted(const std::vector<double> &a) {
  size_t n = a.size();
  if (n == 0) return 0.0;
  if (n & 1) return a[n / 2];
  return 0.5 * (a[n / 2 - 1] + a[n / 2]);
 }

 auto makeSincLike = [&](double scale) {
  return [=](auto &v, auto &ctx) -> Value {
   // sinc(x) = sin(scale*x)/x
   // small t 対策: t = scale*x, |t| < 1e-8 で sin(t)/t をテイラー展開

   if (isComplex(v[0])) {
    Complex x = asComplex(v[0]);
    if (x == Complex(0, 0)) return 1.0;

    Complex t = x * scale;
    double at = std::abs(t);

    if (at < 1e-8) {
     Complex t2 = t * t;
     Complex t4 = t2 * t2;
     Complex t6 = t4 * t2;
     Complex t8 = t4 * t4;

     // sin(t)/t
     Complex s_over_t = Complex(1, 0) - t2 / 6.0 + t4 / 120.0 - t6 / 5040.0 + t8 / 362880.0;

     // sinc(x) = (sin(t)/t) * (t/x)
     return s_over_t * (t / x);
    }

    return std::sin(t) / x;
   }

   double x = asReal(v[0], ctx.pos);
   if (x == 0.0) return 1.0;

   double t = x * scale;
   double at = std::abs(t);

   if (at < 1e-8) {
    double t2 = t * t;
    double t4 = t2 * t2;
    double t6 = t4 * t2;
    double t8 = t4 * t4;

    double s_over_t = 1.0 - t2 / 6.0 + t4 / 120.0 - t6 / 5040.0 + t8 / 362880.0;
    return s_over_t * (t / x);
   }

   return std::sin(t) / x;
  };
 };
 auto makeCoscLike = [&](double scale) {
  return [=](auto &v, auto &ctx) -> Value {
   // cosc(x) = (1 - cos(scale*x)) / x
   // small t 対策: 1 - cos(t) = 2*sin(t/2)^2 を使う
   // t = scale*x

   if (isComplex(v[0])) {
    Complex x = asComplex(v[0]);
    if (x == Complex(0, 0)) return 0.0;

    Complex t = x * scale;             // radians
    Complex h = t * Complex(0.5, 0.0); // t/2

    Complex s = std::sin(h);
    return (Complex(2, 0) * s * s) / x;
   }

   double x = asReal(v[0], ctx.pos);
   if (x == 0.0) return 0.0;

   double t = x * scale;
   double h = 0.5 * t;

   double s = std::sin(h);
   return (2.0 * s * s) / x;
  };
 };
 auto makeTancLike = [&](double scale) {
  return [=](auto &v, auto &ctx) -> Value {
   // tanc(x) = tan(scale*x) / x
   // small t 対策: tan(t)/t をテイラーで計算して (t/x) を掛ける
   // t = scale*x

   if (isComplex(v[0])) {
    Complex x = asComplex(v[0]);
    if (x == Complex(0, 0)) return 1.0;

    Complex t = x * scale;
    double at = std::abs(t);

    if (at < 1e-8) {
     Complex t2 = t * t;
     Complex t4 = t2 * t2;
     Complex t6 = t4 * t2;
     Complex t8 = t4 * t4;

     // tan(t)/t
     Complex tan_over_t = Complex(1, 0) + t2 / 3.0 + (Complex(2, 0) * t4) / 15.0 + (Complex(17, 0) * t6) / 315.0 + (Complex(62, 0) * t8) / 2835.0;

     return tan_over_t * (t / x);
    }

    // 元の仕様: cos が小さいなら inf
    Complex c = std::cos(t);
    if (std::abs(c) < std::pow(10.0, -cnst_precision)) {
     Complex s = std::sin(t);
     return inf(std::real(s) >= 0 ? +1 : -1);
    }

    return std::tan(t) / x;
   }

   double x = asReal(v[0], ctx.pos);
   if (x == 0.0) return 1.0;

   double t = x * scale;
   double at = std::abs(t);

   if (at < 1e-8) {
    double t2 = t * t;
    double t4 = t2 * t2;
    double t6 = t4 * t2;
    double t8 = t4 * t4;

    double tan_over_t = 1.0 + t2 / 3.0 + (2.0 * t4) / 15.0 + (17.0 * t6) / 315.0 + (62.0 * t8) / 2835.0;

    return tan_over_t * (t / x);
   }

   double c = std::cos(t);
   if (std::abs(c) < std::pow(10.0, -cnst_precision)) {
    double s = std::sin(t);
    return inf(s >= 0 ? +1 : -1);
   }

   return std::tan(t) / x;
  };
 };
 static inline double requireFinite(double x, size_t pos) {
  if (!std::isfinite(x)) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), pos);
  return x;
 }

 static inline double toDeg(double rad) { return rad * (180.0 / PI); }
 static inline double toRad(double deg) { return deg * (PI / 180.0); }

 static std::vector<double> gatherReals(const std::vector<Value> &v, size_t pos) {
  std::vector<double> a;
  a.reserve(v.size());
  for (auto &x : v) {
   if (isComplex(x)) {
    Complex z = asComplex(x);
    if (std::imag(z) != 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), pos);
    a.push_back(std::real(z));
   } else {
    a.push_back(asReal(x, pos));
   }
  }
  return a;
 }

 static double meanOf(const std::vector<double> &a) {
  long double s = 0;
  for (double x : a)
   s += (long double)x;
  return (double)(s / (long double)a.size());
 }

 static double variancePopulation(const std::vector<double> &a, double mu) {
  long double ss = 0;
  for (double x : a) {
   long double d = (long double)x - (long double)mu;
   ss += d * d;
  }
  return (double)(ss / (long double)a.size());
 }

 static double stddevPopulation(const std::vector<double> &a, double mu) { return std::sqrt(variancePopulation(a, mu)); }

 static double quantileLinear(std::vector<double> a, double p, size_t pos) {
  if (!(p >= 0.0 && p <= 1.0)) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), pos);
  if (a.empty()) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), pos);

  std::sort(a.begin(), a.end());
  if (a.size() == 1) return a[0];

  double idx = p * (double)(a.size() - 1);
  size_t i = (size_t)std::floor(idx);
  size_t j = std::min(i + 1, a.size() - 1);
  double t = idx - (double)i;
  return a[i] * (1.0 - t) + a[j] * t;
 }

 static std::vector<double> rankAverageTies(const std::vector<double> &x) {
  // average rank for ties, ranks start at 1
  size_t n = x.size();
  std::vector<size_t> idx(n);
  for (size_t i = 0; i < n; ++i)
   idx[i] = i;

  std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });

  std::vector<double> r(n);
  size_t i = 0;
  while (i < n) {
   size_t j = i;
   while (j + 1 < n && x[idx[j + 1]] == x[idx[i]])
    ++j;

   // ranks i..j correspond to 1-based rank (i+1 .. j+1)
   double avg = 0.5 * ((double)(i + 1) + (double)(j + 1));
   for (size_t k = i; k <= j; ++k)
    r[idx[k]] = avg;

   i = j + 1;
  }
  return r;
 }

 static double pearsonCorr(const std::vector<double> &x, const std::vector<double> &y, size_t pos) {
  if (x.size() != y.size() || x.empty()) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), pos);

  double mx = meanOf(x);
  double my = meanOf(y);

  long double sxx = 0, syy = 0, sxy = 0;
  for (size_t i = 0; i < x.size(); ++i) {
   long double dx = (long double)x[i] - (long double)mx;
   long double dy = (long double)y[i] - (long double)my;
   sxx += dx * dx;
   syy += dy * dy;
   sxy += dx * dy;
  }

  if (sxx == 0 || syy == 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), pos);

  return (double)(sxy / std::sqrt(sxx * syy));
 }
 Value roundValue(const Value &v, const SystemConfig &cfg) {
  if (isInvalid(v)) return v;

  auto roundDouble = [&](double d) {
   double scale = std::pow(10.0, (cfg.precision));
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

     const char *begin = src.data() + p;
     const char *end = src.data() + src.size();

     double value;
     auto [ptr, ec] = std::from_chars(begin, end, value);

     if (ec == std::errc::invalid_argument) throw CalcError(CalcErrorType::InvalidNumber, "invalid number", p);

     if (ec == std::errc::result_out_of_range) throw CalcError(CalcErrorType::Overflow, "number out of range", p);

     // 1..1, 1.2.3潰し
     if (ptr < end && *ptr == '.') { throw CalcError(CalcErrorType::InvalidNumber, "multiple '.' in number", static_cast<size_t>(ptr - src.data())); }

     size_t len = static_cast<size_t>(ptr - begin);
     std::string s(src.substr(p, len));

     p += len;
     return {TokenType::Number, s, value, start};
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
     {"EPS", std::pow(10.0, -cnst_precision)},
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

 inline void checkFinite(const std::complex<double> &z, size_t pos) {
  checkFinite(z.real(), pos);
  checkFinite(z.imag(), pos);
 }

 inline void checkFiniteValue(const Value &v, size_t pos) {
  std::visit(
      [&](auto &&x) {
       using T = std::decay_t<decltype(x)>;
       if constexpr (std::is_same_v<T, double>) {
        checkFinite(x, pos);
       } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        checkFinite(x, pos);
       } else if constexpr (std::is_same_v<T, InvalidValue>) {
        // Invalidは別ルートで処理される前提
       } else {
        static_assert(always_false<T>);
       }
      },
      v);
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

 long long permLL(long long n, long long r, size_t pos) {
  if (r < 0 || r > n) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), pos);

  long long res = 1;
  for (long long i = 0; i < r; ++i) {
   res = checkedMul(res, n - i, pos);
  }
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

 long long combLL(long long n, long long r, size_t pos) {
  if (r < 0 || r > n) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), pos);

  r = std::min(r, n - r);
  if (r == 0) return 1;

  long long res = 1;
  for (long long i = 1; i <= r; ++i) {
   long long num = n - r + i; // numerator
   long long den = i;         // denominator

   // まず num と den を約分
   long long g1 = gcdLL(num, den);
   num /= g1;
   den /= g1;

   // 次に res と den を約分
   long long g2 = gcdLL(res, den);
   res /= g2;
   den /= g2;

   // ここで den は 1 になってるはず（整数解なので）
   if (den != 1) {
    // 念のため（理論上は起きない）
    throw CalcError(CalcErrorType::InternalError, "comb: internal reduction failed", pos);
   }

   res = checkedMul(res, num, pos);
  }
  return res;
 }
 Value evalCompare(const Value &lhs, const Value &rhs, CmpOp op, FunctionContext &ctx) {
  // 複素数は大小比較不可
  if (isComplex(lhs) || isComplex(rhs)) throw CalcError(CalcErrorType::NotImplemented, "complex comparison is not implemented", ctx.pos);

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

   NumberNode(Value v, size_t p) : value(std::move(v)) { pos = p; }
   Value evalImpl(SystemConfig &, const std::vector<InputEntry> &, int) const override { return value; }
 };

 /* ---------- Unary ---------- */

 enum class UnaryOp { Plus, Minus };

 struct UnaryNode : ASTNode {
   UnaryOp op;
   std::unique_ptr<ASTNode> rhs;

   UnaryNode(UnaryOp o, std::unique_ptr<ASTNode> r, size_t p) : op(o), rhs(std::move(r)) { pos = p; }

   Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
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

   Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {

    Value v = expr->eval(cfg, hist, base);

    // 後置演算子は実数のみ対応
    if (!std::holds_alternative<double>(v)) { throw CalcError(CalcErrorType::TypeError, errorMessage(CalcErrorType::TypeError), pos); }

    double d = std::get<double>(v);
    double r;

    switch (op) {
     case '!': r = factorial(d, pos); break;
     default: throw CalcError(CalcErrorType::Syntax, "unknown postfix operator", pos);
    }

    return r;
   }
 };

 /* ---------- Compare ---------- */
 struct CompareNode : ASTNode {
   CmpOp op;
   std::unique_ptr<ASTNode> lhs, rhs;

   CompareNode(CmpOp o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p) : op(o), lhs(std::move(l)), rhs(std::move(r)) { pos = p; }

   Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
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

   Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
    auto a = lhs->eval(cfg, hist, base);
    auto b = rhs->eval(cfg, hist, base);

    // ===== complex が含まれる場合 =====
    if (isComplex(a) || isComplex(b)) {
     Complex x = toComplex(a);
     Complex y = toComplex(b);

     // nan注意
     if (!std::isfinite(x.real()) || !std::isfinite(x.imag()) || !std::isfinite(y.real()) || !std::isfinite(y.imag())) { throw CalcError(CalcErrorType::DomainError, "nan or inf", pos); }

     switch (op) {
      case BinOp::Add: return x + y;
      case BinOp::Sub: return x - y;
      case BinOp::Mul: return x * y;
      case BinOp::Div:
       if (std::abs(y) == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
       return x / y;
      case BinOp::Pow: return std::pow(x, y);
     }
    }

    // ===== real 同士 =====
    double x = asDouble(a, pos);
    double y = asDouble(b, pos);

    // nan注意（real）
    if (!std::isfinite(x) || !std::isfinite(y)) throw CalcError(CalcErrorType::DomainError, "nan or inf", pos);

    switch (op) {
     case BinOp::Add: return x + y;
     case BinOp::Sub: return x - y;
     case BinOp::Mul: return x * y;
     case BinOp::Div:
      if (y == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
      return x / y;
     case BinOp::Pow:
      // real→complex 昇格判定はここだけ
      if (x < 0.0 && std::floor(y) != y) {
       // 警告（stderr）
       calcWarn(cfg, pos, "(-a)^(p/q) : principal value only; other branches may exist.if need real, use pow.");
       return std::pow(Complex(x, 0.0), Complex(y, 0.0));
      }
      return std::pow(x, y);
    }

    throw CalcError(CalcErrorType::InvalidOperation, "invalid op", pos);
   }
 };

 /* ---------- Function ---------- */

 struct FunctionCallNode : ASTNode {
   std::string name;
   std::vector<std::unique_ptr<ASTNode>> args;

   Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {

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

    // double のときだけ有限チェック -> 一括捕捉に変更
    /* if (std::holds_alternative<double>(r)) checkFinite(std::get<double>(r), pos);*/

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

   Value evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const override {
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
   Token prev{TokenType::End, "", 0};
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
   // auto r = parsePower(); // 厳格な右辺結合ならpsersePowerがいいが[2^-2]の-が読めなくなる
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

  double eps = std::pow(10.0, -(cnst_precision));
  bool re0 = std::abs(re) < eps;
  bool im0 = std::abs(im) < eps;

  // 純実数
  if (im0) return formatReal(re, cfg);

  std::ostringstream oss;

  // 実部
  if (!re0) oss << formatReal(re, cfg);

  // 虚部の符号
  if (!re0 && im > 0) oss << "+";

  // 虚部の大きさ
  if (std::abs(im - 1.0) < eps) oss << "I";
  else if (std::abs(im + 1.0) < eps) oss << "-I";
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

 Value evaluateFunction(const std::string &name, const std::vector<Value> &args, FunctionContext &ctx) {
  auto it = ctx.cfg.functions.find(name);
  if (it == ctx.cfg.functions.end()) throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), ctx.pos);

  const FunctionDef &f = it->second;

  if (!f.validArgc((int)args.size())) throw CalcError(CalcErrorType::InvalidArgument, errorMessage(CalcErrorType::InvalidArgument), ctx.pos);

  Value r = f.f(args, ctx);

  // double のみ有限性チェック
  if (std::holds_alternative<double>(r)) checkFinite(std::get<double>(r), ctx.pos);

  return r;
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
  cfg.functions["trunc"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::trunc(asReal(v[0], ctx.pos)); }};
  // cfg.functions["round"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::round(asReal(v[0], ctx.pos)); }};

  cfg.functions["pow"] = {2, 2, [](auto &v, auto &ctx) -> Value { return std::pow(asComplex(v[0]), asComplex(v[1])); }};

  // ---- exp / log ----
  cfg.functions["exp"] = {1, 1, [](auto &v, auto &) -> Value { return std::exp(asComplex(v[0])); }};

  // cfg.functions["log"] = {1, 1, [](auto &v, auto &ctx) -> Value {
  //                          if (isComplex(v[0])) return std::log(asComplex(v[0]));
  //                          double x = asReal(v[0], ctx.pos);
  //                          if (x <= 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
  //                          return std::log(x);
  //                         }};
  cfg.functions["log10"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             double x = asReal(v[0], ctx.pos);
                             if (x <= 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                             return std::log10(x);
                            }};

  cfg.functions["log2"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = asReal(v[0], ctx.pos);
                            if (x <= 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                            return std::log2(x);
                           }};

  cfg.functions["log"] = {1, 2, [](auto &v, auto &ctx) -> Value {
                           // ---- complex case ----
                           if (isComplex(v[0]) || (v.size() == 2 && isComplex(v[1]))) {
                            if (v.size() == 1) return std::log(asComplex(v[0]));
                            auto base = asComplex(v[0]);
                            auto x = asComplex(v[1]);
                            Complex lb = std::log(base);
                            // base=1 などで log(base)=0 になると破綻する
                            double eps = asDouble(constants.at("EPS"));
                            if (std::abs(lb) < eps) throw CalcError(CalcErrorType::DomainError, "log: invalid base", ctx.pos);
                            return std::log(x) / lb;
                           }

                           // ---- real case ----
                           if (v.size() == 1) {
                            double x = asReal(v[0], ctx.pos);
                            if (x <= 0) throw CalcError(CalcErrorType::DomainError, "log: x <= 0", ctx.pos);
                            return std::log(x);
                           }

                           // v.size() == 2
                           double base = asReal(v[0], ctx.pos);
                           double x = asReal(v[1], ctx.pos);

                           if (base <= 0 || base == 1 || x <= 0) throw CalcError(CalcErrorType::DomainError, "log: invalid base or x", ctx.pos);

                           return std::log(x) / std::log(base);
                          }};
  cfg.functions["ln"] = cfg.functions["log"];
  cfg.functions["log1p"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (isComplex(v[0])) return std::log(Complex(1, 0) + asComplex(v[0]));

                             double x = asReal(v[0], ctx.pos);
                             if (x <= -1.0) throw CalcError(CalcErrorType::DomainError, "log1p: x <= -1", ctx.pos);
                             return std::log1p(x);
                            }};

  cfg.functions["expm1"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (isComplex(v[0])) return std::exp(asComplex(v[0])) - Complex(1, 0);
                             return std::expm1(asReal(v[0], ctx.pos));
                            }};

  // ---- trig (degree based) ----
  cfg.functions["sin"] = {1, 1, [=](auto &v, auto &) -> Value { return std::sin(asComplex(v[0]) * deg2rad); }};
  cfg.functions["cos"] = {1, 1, [=](auto &v, auto &) -> Value { return std::cos(asComplex(v[0]) * deg2rad); }};
  // cfg.functions["tan"] = {1, 1, [=](auto &v, auto &) -> Value { return std::tan(asComplex(v[0]) * deg2rad); }};
  cfg.functions["tan"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           // complex case
                           if (isComplex(v[0])) { return std::tan(asComplex(v[0]) * deg2rad); }

                           // real case
                           double x = asReal(v[0], ctx.pos);
                           double r = x * deg2rad;
                           double c = std::cos(r);

                           if (std::abs(c) < asDouble(constants.at("EPS"))) {
                            double s = std::sin(r);
                            return inf(s >= 0 ? +1 : -1);
                           }
                           return std::tan(r);
                          }};
  cfg.functions["cot"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (isComplex(v[0])) {
                            auto z = asComplex(v[0]) * deg2rad;
                            return std::cos(z) / std::sin(z);
                           }

                           double x = asReal(v[0], ctx.pos);
                           double r = x * deg2rad;

                           double s = std::sin(r);
                           if (std::abs(s) < asDouble(constants.at("EPS"))) {
                            double c = std::cos(r);
                            return inf(c >= 0 ? +1 : -1);
                           }
                           return std::cos(r) / s;
                          }};

  cfg.functions["sec"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (isComplex(v[0])) {
                            auto z = asComplex(v[0]) * deg2rad;
                            return 1.0 / std::cos(z);
                           }

                           double x = asReal(v[0], ctx.pos);
                           double r = x * deg2rad;

                           double c = std::cos(r);
                           if (std::abs(c) < asDouble(constants.at("EPS"))) {
                            double s = std::sin(r);
                            return inf(s >= 0 ? +1 : -1);
                           }
                           return 1.0 / c;
                          }};

  cfg.functions["csc"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (isComplex(v[0])) {
                            auto z = asComplex(v[0]) * deg2rad;
                            return 1.0 / std::sin(z);
                           }

                           double x = asReal(v[0], ctx.pos);
                           double r = x * deg2rad;

                           double s = std::sin(r);
                           if (std::abs(s) < asDouble(constants.at("EPS"))) {
                            double c = std::cos(r);
                            return inf(c >= 0 ? +1 : -1);
                           }
                           return 1.0 / s;
                          }};

  cfg.functions["asin"] = {1, 1, [=](auto &v, auto &) -> Value { return std::asin(asComplex(v[0])) * rad2deg; }};
  cfg.functions["acos"] = {1, 1, [=](auto &v, auto &) -> Value { return std::acos(asComplex(v[0])) * rad2deg; }};
  cfg.functions["atan"] = {1, 1, [=](auto &v, auto &) -> Value { return std::atan(asComplex(v[0])) * rad2deg; }};
  cfg.functions["atan2"] = {2, 2, [=](auto &v, auto &ctx) -> Value {
                             double y = asReal(v[0], ctx.pos);
                             double x = asReal(v[1], ctx.pos);
                             if (x == 0.0 && y == 0.0) throw CalcError(CalcErrorType::DomainError, "atan2(0,0)", ctx.pos);
                             return std::atan2(y, x) * rad2deg;
                            }};

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
                             if (isComplex(v[0])) return std::acosh(asComplex(v[0]));
                             double x = asReal(v[0], ctx.pos);
                             // 実数範囲なら実数
                             if (x >= 1.0) return std::acosh(x);
                             // それ以外は complex に落とす
                             return std::acosh(Complex(x, 0.0));
                            }};
  cfg.functions["atanh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (isComplex(v[0])) return std::atanh(asComplex(v[0]));
                             double x = asReal(v[0], ctx.pos);
                             // 実数範囲なら実数
                             if (std::abs(x) < 1.0) return std::atanh(x);
                             // |x|>=1 は complex
                             return std::atanh(Complex(x, 0.0));
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

  // ---- sinc / cosc / tanc /sinhc / tanhc / expc  (degree based) ----
  // sinc(x) = sin(x)/x
  // cosc(x) = (1 - cos(x))/x
  // tanc(x) = tan(x)/x
  //
  // この電卓は trig が degree 基準なので、ここも degree 基準。
  //
  // sinhc(x) = sinh(x)/x
  // tanhc(x) = tanh(x)/x
  // expc(x)  = (exp(x)-1)/x = expm1(x)/x
  //
  // x=0 は極限値で定義する
  // sinhc(0)=1, tanhc(0)=1, expc(0)=1

  // cfg.functions["sinc"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
  //                           if (isComplex(v[0])) {
  //                            Complex x = asComplex(v[0]);
  //                            if (x == Complex(0, 0)) return 1.0;
  //                            return std::sin(x * deg2rad) / x;
  //                           }

  //                          double x = asReal(v[0], ctx.pos);
  //                          if (x == 0.0) return 1.0;
  //                          return std::sin(x * deg2rad) / x;
  //                         }};

  cfg.functions["sinc"] = {1, 1, makeSincLike(deg2rad)};
  cfg.functions["sincR"] = {1, 1, makeSincLike(1.0)};

  cfg.functions["cosc"] = {1, 1, makeCoscLike(deg2rad)};
  cfg.functions["coscR"] = {1, 1, makeCoscLike(1.0)};

  cfg.functions["tanc"] = {1, 1, makeTancLike(deg2rad)};
  cfg.functions["tancR"] = {1, 1, makeTancLike(1.0)};

  cfg.functions["sinhc"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (isComplex(v[0])) {
                              Complex x = asComplex(v[0]);
                              if (x == Complex(0, 0)) return 1.0;
                              return std::sinh(x) / x;
                             }

                             double x = asReal(v[0], ctx.pos);
                             if (x == 0.0) return 1.0;
                             return std::sinh(x) / x;
                            }};

  cfg.functions["tanhc"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (isComplex(v[0])) {
                              Complex x = asComplex(v[0]);
                              if (x == Complex(0, 0)) return 1.0;
                              return std::tanh(x) / x;
                             }

                             double x = asReal(v[0], ctx.pos);
                             if (x == 0.0) return 1.0;
                             return std::tanh(x) / x;
                            }};

  cfg.functions["expc"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (isComplex(v[0])) {
                             Complex x = asComplex(v[0]);
                             if (x == Complex(0, 0)) return 1.0;
                             return (std::exp(x) - Complex(1, 0)) / x;
                            }

                            double x = asReal(v[0], ctx.pos);
                            if (x == 0.0) return 1.0;

                            // exp(x)-1 の桁落ち回避
                            return std::expm1(x) / x;
                           }};

  // ---- number theory / statistics ----
  cfg.functions["fact"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            long long n = requireInt(v[0], ctx.pos);
                            return (double)factLD(n, ctx.pos);
                           }};
  cfg.functions["gcd"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           long long a = requireInt(v[0], ctx.pos);
                           long long b = requireInt(v[1], ctx.pos);
                           return (double)gcdLL(a, b);
                          }};
  cfg.functions["comb"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            long long n = requireInt(v[0], ctx.pos);
                            long long r = requireInt(v[1], ctx.pos);
                            return (double)combLL(n, r, ctx.pos);
                           }};
  cfg.functions["perm"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            long long n = requireInt(v[0], ctx.pos);
                            long long r = requireInt(v[1], ctx.pos);
                            return (double)permLL(n, r, ctx.pos);
                           }};
  cfg.functions["lcm"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           long long a = requireInt(v[0], ctx.pos);
                           long long b = requireInt(v[1], ctx.pos);
                           if (a == 0 || b == 0) return 0.0;
                           long long g = gcdLL(a, b);
                           // lcm = |a/g| * |b|
                           long long x = a / g;
                           long long r = checkedMul(std::llabs(x), std::llabs(b), ctx.pos);
                           return (double)r;
                          }};

  cfg.functions["sum"] = {0, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.empty()) throw CalcError(CalcErrorType::DomainError, "sum: no elements", ctx.pos);
                           Complex acc{0, 0};
                           for (auto &x : v)
                            acc += asComplex(x);
                           return (std::abs(acc.imag()) < std::pow(10.0, -cnst_precision)) ? Value(acc.real()) : Value(acc);
                          }};

  cfg.functions["prod"] = {0, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.empty()) throw CalcError(CalcErrorType::DomainError, "prod: no elements", ctx.pos);
                            Complex acc{1, 0};
                            for (auto &x : v)
                             acc *= asComplex(x);

                            // 実数ならdouble返し
                            return (std::abs(acc.imag()) < std::pow(10.0, -cnst_precision)) ? Value(acc.real()) : Value(acc);
                           }};

  cfg.functions["mean"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.empty()) throw CalcError(CalcErrorType::DomainError, "mean: no elements", ctx.pos);
                            Complex acc{0, 0};
                            for (auto &x : v)
                             acc += asComplex(x);
                            acc /= (double)v.size();
                            return (std::abs(acc.imag()) < std::pow(10.0, -cnst_precision)) ? Value(acc.real()) : Value(acc);
                           }};
  cfg.functions["mod"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           Complex a = asComplex(v[0]), b = asComplex(v[1]);
                           if (a.imag() != 0 || b.imag() != 0) throw CalcError(CalcErrorType::DomainError, "mod: complex argument", ctx.pos);
                           double x = a.real(), y = b.real();
                           if (y == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), ctx.pos);
                           double r = x - y * std::floor(x / y);
                           return r;
                          }};
  cfg.functions["geomean"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                               int n = (int)v.size();
                               long double sumLog = 0.0L;

                               for (auto &x : v) {
                                if (isComplex(x)) throw CalcError(CalcErrorType::DomainError, "geomean: complex", ctx.pos);

                                double a = asReal(x, ctx.pos);
                                if (a < 0) throw CalcError(CalcErrorType::DomainError, "geomean: negative", ctx.pos);
                                if (a == 0) return 0.0;

                                sumLog += std::log((long double)a);
                               }

                               long double r = std::exp(sumLog / (long double)n);
                               if (!std::isfinite((double)r)) throw CalcError(CalcErrorType::Overflow, "geomean: overflow", ctx.pos);
                               return (double)r;
                              }};
  cfg.functions["harmmean"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                                int n = (int)v.size();
                                long double acc = 0.0L;

                                for (auto &x : v) {
                                 if (isComplex(x)) throw CalcError(CalcErrorType::DomainError, "harmmean: complex", ctx.pos);

                                 double a = asReal(x, ctx.pos);
                                 if (a == 0.0) throw CalcError(CalcErrorType::DivisionByZero, "harmmean: zero element", ctx.pos);
                                 acc += 1.0L / (long double)a;
                                }

                                long double r = (long double)n / acc;
                                if (!std::isfinite((double)r)) throw CalcError(CalcErrorType::Overflow, "harmmean: overflow", ctx.pos);
                                return (double)r;
                               }};
  cfg.functions["quantile"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                double p = asReal(v[0], ctx.pos);
                                if (!(p >= 0.0 && p <= 1.0)) throw CalcError(CalcErrorType::DomainError, "quantile: p out of range", ctx.pos);

                                if (v.size() < 2) throw CalcError(CalcErrorType::DomainError, "quantile: no samples", ctx.pos);

                                std::vector<double> a;
                                a.reserve(v.size() - 1);
                                for (size_t i = 1; i < v.size(); ++i)
                                 a.push_back(asReal(v[i], ctx.pos));

                                std::sort(a.begin(), a.end());

                                double idx = p * (a.size() - 1);
                                size_t i0 = (size_t)std::floor(idx);
                                size_t i1 = (size_t)std::ceil(idx);

                                if (i0 == i1) return a[i0];
                                return a[i0] * (i1 - idx) + a[i1] * (idx - i0);
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
  cfg.functions["clamp"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                             double x = asReal(v[0], ctx.pos);
                             double lo = asReal(v[1], ctx.pos);
                             double hi = asReal(v[2], ctx.pos);
                             if (lo > hi) throw CalcError(CalcErrorType::DomainError, "clamp: lo > hi", ctx.pos);
                             return std::min(std::max(x, lo), hi);
                            }};
  cfg.functions["fract"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             double x = asReal(v[0], ctx.pos);
                             return x - std::floor(x);
                            }};
  cfg.functions["gamma"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             double x = asReal(v[0], ctx.pos);
                             double r = std::tgamma(x);
                             if (std::isnan(r)) throw CalcError(CalcErrorType::DomainError, "gamma: domain error", ctx.pos);
                             if (!std::isfinite(r)) throw CalcError(CalcErrorType::Overflow, "gamma: overflow", ctx.pos);
                             return r;
                            }};
  cfg.functions["lgamma"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              double x = asReal(v[0], ctx.pos);
                              double r = std::lgamma(x);
                              if (std::isnan(r)) throw CalcError(CalcErrorType::DomainError, "lgamma: domain error", ctx.pos);
                              if (!std::isfinite(r)) throw CalcError(CalcErrorType::Overflow, "lgamma: overflow", ctx.pos);

                              return r;
                             }};
  cfg.functions["mode"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                            auto a = collectReals(v, ctx);
                            std::sort(a.begin(), a.end());

                            // double比較用トレランス
                            const double tol = std::pow(10.0, -cnst_precision);

                            auto eq = [&](double x, double y) -> bool { return std::abs(x - y) <= tol; };

                            double best = a[0];
                            int bestCnt = 1;

                            double cur = a[0];
                            int curCnt = 1;

                            for (size_t i = 1; i < a.size(); ++i) {
                             if (eq(a[i], cur)) {
                              curCnt++;
                             } else {
                              if (curCnt > bestCnt) {
                               bestCnt = curCnt;
                               best = cur;
                              }
                              cur = a[i];
                              curCnt = 1;
                             }
                            }

                            // last
                            if (curCnt > bestCnt) {
                             bestCnt = curCnt;
                             best = cur;
                            }

                            // 全部1回ずつなら「最小値」が返る（仕様としていいのかな?）
                            return best;
                           }};

  cfg.functions["erf"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           double x = asReal(v[0], ctx.pos);
                           return std::erf(x);
                          }};
  cfg.functions["erfc"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = asReal(v[0], ctx.pos);
                            return std::erfc(x);
                           }};
  cfg.functions["round"] = {1, 2, [](auto &v, auto &ctx) -> Value {
                             double x = asReal(v[0], ctx.pos);

                             if (v.size() == 1) return std::round(x);

                             int n = (int)requireInt(v[1], ctx.pos);

                             // double の安全域（10^±15くらいが現実的）
                             if (n < -15 || n > 15) throw CalcError(CalcErrorType::OutOfRange, "round: n out of range (-15..15)", ctx.pos);

                             double s = std::pow(10.0, (double)n);
                             return std::round(x * s) / s;
                            }};

  cfg.functions["var"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                           int n = (int)v.size();
                           if (n == 0) throw CalcError(CalcErrorType::DomainError, "var: no elements", ctx.pos);
                           if (n == 1) return 0.0;

                           double mean = 0;
                           for (auto &x : v)
                            mean += asReal(x, ctx.pos);
                           mean /= n;

                           double acc = 0;
                           for (auto &x : v) {
                            double d = asReal(x, ctx.pos) - mean;
                            acc += d * d;
                           }
                           return acc / n;
                          }};
  cfg.functions["vars"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            int n = (int)v.size();
                            if (n < 2) throw CalcError(CalcErrorType::DomainError, "vars: too few elements", ctx.pos);

                            double mean = 0;
                            for (auto &x : v)
                             mean += asReal(x, ctx.pos);
                            mean /= n;

                            double acc = 0;
                            for (auto &x : v) {
                             double d = asReal(x, ctx.pos) - mean;
                             acc += d * d;
                            }
                            return acc / (n - 1);
                           }};
  cfg.functions["stddev"] = {1, -1, [](auto &v, auto &ctx) -> Value { return std::sqrt(asReal(evaluateFunction("var", v, ctx), ctx.pos)); }};

  cfg.functions["stddevs"] = {1, -1, [](auto &v, auto &ctx) -> Value { return std::sqrt(asReal(evaluateFunction("vars", v, ctx), ctx.pos)); }};

  cfg.functions["median"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                              if (v.empty()) throw CalcError(CalcErrorType::DomainError, "median: no elements", ctx.pos);

                              std::vector<double> a;
                              a.reserve(v.size());
                              for (auto &x : v)
                               a.push_back(asReal(x, ctx.pos));

                              std::sort(a.begin(), a.end());
                              int n = (int)a.size();
                              if (n % 2 == 1) return a[n / 2];
                              return 0.5 * (a[n / 2 - 1] + a[n / 2]);
                             }};

  cfg.functions["mad"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                           auto a = collectReals(v, ctx);
                           std::sort(a.begin(), a.end());

                           double med = medianOfSorted(a);

                           std::vector<double> dev;
                           dev.reserve(a.size());
                           for (double x : a)
                            dev.push_back(std::abs(x - med));

                           std::sort(dev.begin(), dev.end());
                           return medianOfSorted(dev);
                          }};
  cfg.functions["skew"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                            auto a = collectReals(v, ctx);
                            size_t n = a.size();

                            double mean = 0.0;
                            for (double x : a)
                             mean += x;
                            mean /= (double)n;

                            double m2 = 0.0, m3 = 0.0;
                            for (double x : a) {
                             double d = x - mean;
                             double d2 = d * d;
                             m2 += d2;
                             m3 += d2 * d;
                            }
                            m2 /= (double)n;
                            m3 /= (double)n;

                            if (m2 == 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                            return m3 / std::pow(m2, 1.5);
                           }};

  cfg.functions["kurt"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                            auto a = collectReals(v, ctx);
                            size_t n = a.size();

                            double mean = 0.0;
                            for (double x : a)
                             mean += x;
                            mean /= (double)n;

                            double m2 = 0.0, m4 = 0.0;
                            for (double x : a) {
                             double d = x - mean;
                             double d2 = d * d;
                             m2 += d2;
                             m4 += d2 * d2;
                            }
                            m2 /= (double)n;
                            m4 /= (double)n;

                            if (m2 == 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                            return (m4 / (m2 * m2)) - 3.0;
                           }};
  cfg.functions["percentile"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                  double p = asReal(v[0], ctx.pos);
                                  if (p < 0 || p > 100) throw CalcError(CalcErrorType::DomainError, "percentile: p out of range", ctx.pos);

                                  std::vector<double> a;
                                  for (size_t i = 1; i < v.size(); ++i)
                                   a.push_back(asReal(v[i], ctx.pos));

                                  std::sort(a.begin(), a.end());
                                  double idx = (p / 100.0) * (a.size() - 1);
                                  size_t i0 = (size_t)std::floor(idx);
                                  size_t i1 = (size_t)std::ceil(idx);
                                  if (i0 == i1) return a[i0];
                                  return a[i0] * (i1 - idx) + a[i1] * (idx - i0);
                                 }};
  cfg.functions["cov"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.size() % 2 != 0) throw CalcError(CalcErrorType::DomainError, "cov: invalid argument count", ctx.pos);

                           int n = (int)v.size() / 2;
                           if (n < 2) throw CalcError(CalcErrorType::DomainError, "cov: too few samples", ctx.pos);

                           double mx = 0, my = 0;
                           for (int i = 0; i < n; ++i) {
                            mx += asReal(v[i], ctx.pos);
                            my += asReal(v[i + n], ctx.pos);
                           }
                           mx /= n;
                           my /= n;

                           double acc = 0;
                           for (int i = 0; i < n; ++i)
                            acc += (asReal(v[i], ctx.pos) - mx) * (asReal(v[i + n], ctx.pos) - my);

                           return acc / n;
                          }};
  cfg.functions["corr"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                            int total = (int)v.size();
                            if (total % 2 != 0) throw CalcError(CalcErrorType::InvalidArgument, "corr: argument count must be even", ctx.pos);

                            int n = total / 2;
                            std::vector<Value> vx(v.begin(), v.begin() + n);
                            std::vector<Value> vy(v.begin() + n, v.end());
                            double c = asReal(evaluateFunction("cov", v, ctx), ctx.pos);
                            double sx = asReal(evaluateFunction("var", vx, ctx), ctx.pos);
                            double sy = asReal(evaluateFunction("var", vy, ctx), ctx.pos);
                            if (sx == 0.0 || sy == 0.0) throw CalcError(CalcErrorType::DomainError, "corr: zero variance", ctx.pos);
                            return c / std::sqrt(sx * sy);
                           }};

  // ---- geometric vector ----
  cfg.functions["norm"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.empty()) throw CalcError(CalcErrorType::DomainError, "norm: no elements", ctx.pos);

                            double acc = 0.0;
                            for (auto &x : v) {
                             double r = asReal(x, ctx.pos);
                             acc += r * r;
                            }
                            return std::sqrt(acc);
                           }};
  cfg.functions["dot"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.size() % 2 != 0) throw CalcError(CalcErrorType::DomainError, "dot: dimension mismatch", ctx.pos);

                           int n = (int)v.size() / 2;
                           if (n == 0) throw CalcError(CalcErrorType::DomainError, "dot: no elements", ctx.pos);

                           double acc = 0.0;
                           for (int i = 0; i < n; ++i) {
                            acc += asReal(v[i], ctx.pos) * asReal(v[i + n], ctx.pos);
                           }
                           return acc;
                          }};

  // cfg.functions["cross"] = {6, 6, [](auto &v, auto &ctx) -> Value {
  //                            double x1 = asReal(v[0], ctx.pos);
  //                            double y1 = asReal(v[1], ctx.pos);
  //                            double z1 = asReal(v[2], ctx.pos);
  //                            double x2 = asReal(v[3], ctx.pos);
  //                            double y2 = asReal(v[4], ctx.pos);
  //                            double z2 = asReal(v[5], ctx.pos);
  //
  //                            // (x1,y1,z1) × (x2,y2,z2)
  //                            // → ベクトルとして返したいが、現状 Value が scalar 前提なので
  //                            // テストに合わせて "(x,y,z)" 形式の文字列 Value を返す想定
  //                            return VectorValue{y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2};
  //                           }};

  cfg.functions["lerp"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                            double a = asReal(v[0], ctx.pos);
                            double b = asReal(v[1], ctx.pos);
                            double t = asReal(v[2], ctx.pos);
                            return a * (1 - t) + b * t;
                           }};
  cfg.functions["distance"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                if (v.size() % 2 != 0) throw CalcError(CalcErrorType::DomainError, "distance: dimension mismatch", ctx.pos);

                                int n = (int)v.size() / 2;
                                if (n == 0) throw CalcError(CalcErrorType::DomainError, "distance: no elements", ctx.pos);

                                double acc = 0.0;
                                for (int i = 0; i < n; ++i) {
                                 double d = asReal(v[i + n], ctx.pos) - asReal(v[i], ctx.pos);
                                 acc += d * d;
                                }
                                return std::sqrt(acc);
                               }};
  cfg.functions["fib"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           long long n = requireInt(v[0], ctx.pos);
                           return (double)fibULL(n, ctx.pos);
                          }};
  cfg.functions["isprime"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                               long long n = requireInt(v[0], ctx.pos);
                               return isPrimeLL(n) ? 1.0 : 0.0;
                              }};

  // ---- control ----
  cfg.functions["if"] = {3, 3, [](auto &v, auto &ctx) -> Value { return asReal(v[0], ctx.pos) != 0.0 ? v[1] : v[2]; }};
  cfg.functions["rand"] = {0, 2, [](auto &v, auto &ctx) -> Value {
                            static thread_local std::mt19937_64 rng{std::random_device{}()};

                            auto make = [&](double lo, double hi) -> Value {
                             // NaN/Inf guard
                             checkFinite(lo, ctx.pos);
                             checkFinite(hi, ctx.pos);

                             if (lo > hi) std::swap(lo, hi);

                             // degenerate case
                             if (lo == hi) return lo;

                             // [lo, hi)
                             std::uniform_real_distribution<double> dist(lo, hi);
                             return dist(rng);
                            };
                            if (v.empty()) { return make(0.0, 1.0); }
                            if (v.size() == 1) {
                             double hi = asReal(v[0], ctx.pos);
                             return make(0.0, hi);
                            }
                            double lo = asReal(v[0], ctx.pos);
                            double hi = asReal(v[1], ctx.pos);
                            return make(lo, hi);
                           }};
  cfg.functions["randint"] = {0, 2, [](auto &v, auto &ctx) -> Value {
                               static thread_local std::mt19937_64 rng{std::random_device{}()};
                               auto gen = [&](long long lo, long long hi) -> Value {
                                if (lo > hi) std::swap(lo, hi);

                                // 退化ケース
                                if (lo == hi) return (double)lo;
                                std::uniform_int_distribution<long long> dist(lo, hi);
                                return (double)dist(rng);
                               };
                               // 0 or 1
                               if (v.empty()) { return gen(0, 1); }
                               if (v.size() == 1) {
                                long long n = requireInt(v[0], ctx.pos);
                                return gen(0, n); // n<0 でも gen が swap する
                               }
                               long long a = requireInt(v[0], ctx.pos);
                               long long b = requireInt(v[1], ctx.pos);
                               return gen(a, b);
                              }};
  cfg.functions["randn"] = {0, 2, [](auto &v, auto &ctx) -> Value {
                             static thread_local std::mt19937_64 rng{std::random_device{}()};

                             double mu = 0.0;
                             double sigma = 1.0;

                             if (v.size() >= 1) mu = asReal(v[0], ctx.pos);
                             if (v.size() >= 2) sigma = asReal(v[1], ctx.pos);

                             if (sigma < 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                             std::normal_distribution<double> dist(mu, sigma);
                             return dist(rng);
                            }};
  cfg.functions["choice"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                              static thread_local std::mt19937_64 rng{std::random_device{}()};

                              if (v.empty()) throw CalcError(CalcErrorType::InvalidArgument, "choice: no arguments", ctx.pos);

                              std::uniform_int_distribution<size_t> dist(0, v.size() - 1);
                              return v[dist(rng)];
                             }};
  // fma(a,b,c)
  cfg.functions["fma"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                           double a = asReal(v[0], ctx.pos);
                           double b = asReal(v[1], ctx.pos);
                           double c = asReal(v[2], ctx.pos);
                           return std::fma(a, b, c);
                          }};

  // ave(...) = mean(...)
  cfg.functions["ave"] = cfg.functions["mean"];

  // rms(...)
  cfg.functions["rms"] = {1, 9999, [](auto &v, auto &ctx) -> Value {
                           auto a = gatherReals(v, ctx.pos);
                           if (a.empty()) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                           long double ss = 0;
                           for (double x : a)
                            ss += (long double)x * (long double)x;
                           return std::sqrt((double)(ss / (long double)a.size()));
                          }};

  // cv(...)
  cfg.functions["cv"] = {1, 9999, [](auto &v, auto &ctx) -> Value {
                          auto a = gatherReals(v, ctx.pos);
                          if (a.size() < 2) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                          double mu = meanOf(a);
                          if (mu == 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                          double sd = stddevPopulation(a, mu);
                          return sd / mu;
                         }};

  // stderr(...)
  cfg.functions["stderr"] = {1, 9999, [](auto &v, auto &ctx) -> Value {
                              auto a = gatherReals(v, ctx.pos);
                              if (a.size() < 2) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                              double mu = meanOf(a);
                              double sd = stddevPopulation(a, mu);
                              return sd / std::sqrt((double)a.size());
                             }};

  // zscore(x, mu, sigma)
  cfg.functions["zscore"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                              double x = asReal(v[0], ctx.pos);
                              double mu = asReal(v[1], ctx.pos);
                              double sigma = asReal(v[2], ctx.pos);
                              if (sigma == 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                              return (x - mu) / sigma;
                             }};

  // iqr(...)
  cfg.functions["iqr"] = {1, 9999, [](auto &v, auto &ctx) -> Value {
                           auto a = gatherReals(v, ctx.pos);
                           if (a.size() < 2) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                           double q1 = quantileLinear(a, 0.25, ctx.pos);
                           double q3 = quantileLinear(a, 0.75, ctx.pos);
                           return q3 - q1;
                          }};

  // trimmean(p, ...)
  cfg.functions["trimmean"] = {2, 9999, [](auto &v, auto &ctx) -> Value {
                                double p = asReal(v[0], ctx.pos);
                                if (!(p >= 0.0 && p < 0.5)) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                                std::vector<Value> rest(v.begin() + 1, v.end());
                                auto a = gatherReals(rest, ctx.pos);
                                if (a.empty()) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                                std::sort(a.begin(), a.end());
                                size_t n = a.size();
                                size_t k = (size_t)std::floor(p * (double)n);

                                if (2 * k >= n) k = (n - 1) / 2;

                                long double s = 0;
                                size_t cnt = 0;
                                for (size_t i = k; i < n - k; ++i) {
                                 s += (long double)a[i];
                                 ++cnt;
                                }
                                if (cnt == 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                                return (double)(s / (long double)cnt);
                               }};

  // winsor(p, ...)  -> returns mean of winsorized data
  cfg.functions["winsor"] = {2, 9999, [](auto &v, auto &ctx) -> Value {
                              double p = asReal(v[0], ctx.pos);
                              if (!(p >= 0.0 && p < 0.5)) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                              std::vector<Value> rest(v.begin() + 1, v.end());
                              auto a = gatherReals(rest, ctx.pos);
                              if (a.empty()) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                              std::sort(a.begin(), a.end());
                              size_t n = a.size();
                              size_t k = (size_t)std::floor(p * (double)n);
                              if (k >= n) k = n - 1;

                              double lo = a[k];
                              double hi = a[n - 1 - k];

                              for (double &x : a) {
                               if (x < lo) x = lo;
                               if (x > hi) x = hi;
                              }
                              return meanOf(a);
                             }};

  // corrspearman(x..., y...)
  cfg.functions["corrspearman"] = {2, 9999, [](auto &v, auto &ctx) -> Value {
                                    // split by half
                                    if (v.size() % 2 != 0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                                    size_t n = v.size() / 2;
                                    if (n < 2) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);

                                    std::vector<Value> vx(v.begin(), v.begin() + n);
                                    std::vector<Value> vy(v.begin() + n, v.end());

                                    auto x = gatherReals(vx, ctx.pos);
                                    auto y = gatherReals(vy, ctx.pos);

                                    auto rx = rankAverageTies(x);
                                    auto ry = rankAverageTies(y);

                                    return pearsonCorr(rx, ry, ctx.pos);
                                   }};

  // re(z)
  cfg.functions["re"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                          if (isComplex(v[0])) return std::real(asComplex(v[0]));
                          return asReal(v[0], ctx.pos);
                         }};

  // im(z)
  cfg.functions["im"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                          if (isComplex(v[0])) return std::imag(asComplex(v[0]));
                          return 0.0;
                         }};

  // arg(z) -> degree
  cfg.functions["arg"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           Complex z = isComplex(v[0]) ? asComplex(v[0]) : Complex(asReal(v[0], ctx.pos), 0.0);
                           if (z == Complex(0, 0)) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                           return toDeg(std::atan2(std::imag(z), std::real(z)));
                          }};

  // conj(z)
  cfg.functions["conj"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            Complex z = isComplex(v[0]) ? asComplex(v[0]) : Complex(asReal(v[0], ctx.pos), 0.0);
                            return std::conj(z);
                           }};

  // polar(r, theta)  theta: degree
  cfg.functions["polar"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             double r = asReal(v[0], ctx.pos);
                             double th = asReal(v[1], ctx.pos);
                             double t = toRad(th);
                             return Complex(r * std::cos(t), r * std::sin(t));
                            }};

  // cis(x) x: degree
  cfg.functions["cis"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           double th = asReal(v[0], ctx.pos);
                           double t = toRad(th);
                           return Complex(std::cos(t), std::sin(t));
                          }};

  // stress(F, A) = F/A
  cfg.functions["stress"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              double F = asReal(v[0], ctx.pos);
                              double A = asReal(v[1], ctx.pos);
                              if (A == 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                              return F / A;
                             }};

  // strain(dL, L) = dL/L
  cfg.functions["strain"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              double dL = asReal(v[0], ctx.pos);
                              double L = asReal(v[1], ctx.pos);
                              if (L == 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                              return dL / L;
                             }};

  // young(sigma, eps) = sigma/eps
  cfg.functions["young"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             double sigma = asReal(v[0], ctx.pos);
                             double eps = asReal(v[1], ctx.pos);
                             if (eps == 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                             return sigma / eps;
                            }};

  // moment_rect(b,h) = b*h^3/12
  cfg.functions["moment_rect"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                   double b = asReal(v[0], ctx.pos);
                                   double h = asReal(v[1], ctx.pos);
                                   return b * h * h * h / 12.0;
                                  }};

  // moment_circle(d) = pi*d^4/64
  cfg.functions["moment_circle"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                     double d = asReal(v[0], ctx.pos);
                                     return PI * d * d * d * d / 64.0;
                                    }};

  // sectionmod_rect(b,h) = b*h^2/6
  cfg.functions["sectionmod_rect"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                       double b = asReal(v[0], ctx.pos);
                                       double h = asReal(v[1], ctx.pos);
                                       return b * h * h / 6.0;
                                      }};

  // sectionmod_circle(d) = pi*d^3/32
  cfg.functions["sectionmod_circle"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                         double d = asReal(v[0], ctx.pos);
                                         return PI * d * d * d / 32.0;
                                        }};

  // torsion_J_circle(d) = pi*d^4/32
  cfg.functions["torsion_J_circle"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                        double d = asReal(v[0], ctx.pos);
                                        return PI * d * d * d * d / 32.0;
                                       }};

  // polarZ_circle(d) = pi*d^3/16
  cfg.functions["polarZ_circle"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                     double d = asReal(v[0], ctx.pos);
                                     return PI * d * d * d / 16.0;
                                    }};

  // bolt_stress(F, d) = F / (pi d^2 / 4)
  cfg.functions["bolt_stress"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                   double F = asReal(v[0], ctx.pos);
                                   double d = asReal(v[1], ctx.pos);
                                   if (d == 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                                   double A = PI * d * d / 4.0;
                                   return F / A;
                                  }};

  // torque_from_preload(F, d, K) = K*F*d
  cfg.functions["torque_from_preload"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                           double F = asReal(v[0], ctx.pos);
                                           double d = asReal(v[1], ctx.pos);
                                           double K = asReal(v[2], ctx.pos);
                                           return K * F * d;
                                          }};

  // preload_from_torque(T, d, K) = T/(K*d)
  cfg.functions["preload_from_torque"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                           double T = asReal(v[0], ctx.pos);
                                           double d = asReal(v[1], ctx.pos);
                                           double K = asReal(v[2], ctx.pos);
                                           if (K == 0.0 || d == 0.0) throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), ctx.pos);
                                           return T / (K * d);
                                          }};

  // friction(mu, N) = mu*N
  cfg.functions["friction"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                double mu = asReal(v[0], ctx.pos);
                                double N = asReal(v[1], ctx.pos);
                                return mu * N;
                               }};

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
   res.value = v;
   return res;
  } catch (const ClearRequest &) {
   rt.shouldClear = true;
   res.kind = EvalKind::Clear;
   return res;
  } catch (const ExitRequest &) {
   rt.shouldExit = true;
   std::cout << "bye...nara\n";
   res.kind = EvalKind::Exit;
   return res;
  }
 }
} // namespace mm::cal

using namespace mm::cal;
/* ============================
  DLL
  ============================ */
struct CalcResult {
  bool ok;
  std::string output;  // 成功時: 計算結果 / 失敗時: エラーメッセージ
  size_t errorPos = 0; // エラー時のみ
};

CalcResult calcEval(const std::string &expr, SystemConfig &cfg, std::vector<InputEntry> &history, int base = 10);

CalcResult calcEval(const std::string &expr, SystemConfig &cfg, std::vector<InputEntry> &history, int base) {
 try {
  Value v = evaluate(expr, cfg, history, base);

  // 表示用に丸め
  Value vr = roundValue(v, cfg);
  std::string out = formatResult(vr, cfg);

  history.push_back({expr, vr});

  return {true, out, 0};

 } catch (const CalcError &e) { return {false, e.what(), e.pos}; } catch (const std::exception &e) {
  return {false, e.what(), 0};

 } catch (...) { return {false, "unknown error", 0}; }
}

extern "C" __declspec(dllexport) const char *CalcEvalDLL(const char *expr, SystemConfig *cfg, std::vector<InputEntry> *history, int base) {
 static thread_local std::string result; // DLL境界安全

 try {
  Value v = evaluate(expr, *cfg, *history, base);

  Value vr = roundValue(v, *cfg);
  result = formatResult(vr, *cfg);

  history->push_back({expr, vr});
  return result.c_str();

 } catch (const ExitRequest &) {
  result = "ERROR: Exit requested";
  return result.c_str();

 } catch (const CalcError &e) {
  result = "ERROR: ";
  result += e.what();
  return result.c_str();

 } catch (const std::exception &e) {
  result = "ERROR: ";
  result += e.what();
  return result.c_str();

 } catch (...) {
  result = "ERROR: unknown error";
  return result.c_str();
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
              "================================\n\n";

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
     std::cout << "\nOut[" << history.size() << "] := " << formatResult(res.value, syscfg) << "\n\n";
     break;

    case EvalKind::None: break;
   }
  } catch (const CalcError &e) { std::cout << "\nError: " << e.what() << "\n" << line << "\n" << std::string(e.pos, ' ') << "^\n\n"; }
 }
 std::cout << "bye...nara\n";
}