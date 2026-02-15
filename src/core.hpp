#pragma once

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

namespace mm::cal {
 inline constexpr double PI = 3.14159265358979323846264338327950288419716939937510;
 inline constexpr int cnst_precision = 12;
 inline constexpr double cnst_precision_inv = 1e-12;

 struct InvalidValue {};
 struct ExitRequest {};
 struct ClearRequest {};
 struct SystemConfig;
 struct FunctionContext;
 struct FunctionDef;
 struct InputEntry;
 using Value = std::variant<double, std::complex<double>, std::string, InvalidValue>;
 // 以下のほうが高速になる。可読性とのトレードオフ
 // struct Value {
 //  enum { Real, Complex } type;
 //  double r;
 //  std::complex<double> c;
 //};
 using Complex = std::complex<double>;

 struct FunctionContext;
 struct RuntimeState {
   bool shouldExit = false;
   bool shouldClear = false;
 };
 struct InputEntry {
   std::string expr; // In[n]
   Value value;      // Out[n]
 };
 struct CalcResult {
   bool ok;
   std::string output;  // 成功時: 計算結果 / 失敗時: エラーメッセージ
   size_t errorPos = 0; // エラー時のみ
 };
 enum class EvalKind { None, Value, Clear, Exit };
 struct EvalResult {
   EvalKind kind = EvalKind::None;
   Value value{};
 };

 using FuncImpl = std::function<Value(const std::vector<Value> &, FunctionContext &)>;
 // using Fn = std::function<Value(std::span<const Value>, FunctionContext &)>;  // vectorは重いのでそのうち固定長バッファに切り替える
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
   const double deg2rad = PI / 180.0;
   const double rad2deg = 180.0 / PI;
 };

 inline const std::unordered_map<std::string, Value> constants = {
     {"Pi", PI},
     {"E", 2.71828182845904523536028747135266249775724709369995},
     {"Tau", 2 * PI},
     {"Phi", 1.61803398874989484820458683436563811772030917980576},
     {"NA", 6.02214076e23},
     {"I", Complex(0.0, 1.0)},
     {"EPS", cnst_precision_inv},
 };

 template <class> inline constexpr bool always_false = false;

 EvalResult evalLine(const std::string &line, SystemConfig &cfg, RuntimeState &rt, std::vector<InputEntry> &history);
 CalcResult calcEval(const std::string &expr, SystemConfig &cfg, std::vector<InputEntry> &history, int base = 10);

 /* ============================
    エラー
    ============================ */

 enum class CalcErrorType {
  NotImplemented,
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
  NonConvergence,
  Overflow,
  InvalidNumber,
  InvalidOperation,
  TypeError,
  InvalidArgument,
  NaNDetected,
  InfinityDetected,
  InternalError,
  SyntaxError,
  RuntimeError,
 };

 constexpr const char *errorMessage(CalcErrorType t) {
  switch (t) {
   case CalcErrorType::NotImplemented: return "this feature is not yet implemented";
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
   case CalcErrorType::NonConvergence: return "did not converge";
   case CalcErrorType::Overflow: return "overflow detected";
   case CalcErrorType::InvalidNumber: return "invalid number";
   case CalcErrorType::InvalidOperation: return "invalid operation";
   case CalcErrorType::TypeError: return "double or complex type error";
   case CalcErrorType::InvalidArgument: return "invalid argument count";
   case CalcErrorType::NaNDetected: return "NaN detected";
   case CalcErrorType::InfinityDetected: return "Infinity detected";
   case CalcErrorType::InternalError: return "InternalError";
   case CalcErrorType::SyntaxError: return "SyntaxError";
   case CalcErrorType::RuntimeError: return "RuntimeError";
  }
  return "unknown calculation error";
 }

 struct CalcError : std::runtime_error {
   CalcErrorType type;
   size_t pos;
   CalcError(CalcErrorType t, const std::string &m, size_t p) : std::runtime_error(m), type(t), pos(p) {}
 };

 [[noreturn]] inline void throwOverflow(size_t pos) { throw CalcError(CalcErrorType::Overflow, errorMessage(CalcErrorType::Overflow), pos); }
 [[noreturn]] inline void throwDomain(size_t pos) { throw CalcError(CalcErrorType::DomainError, errorMessage(CalcErrorType::DomainError), pos); }
 [[noreturn]] inline void throwDomain(size_t pos, const char *msg) { throw CalcError(CalcErrorType::DomainError, msg, pos); }
 [[noreturn]] inline void throwInvalid(size_t pos, const char *msg) { throw CalcError(CalcErrorType::InvalidArgument, msg, pos); }

 /* ============================
    補助関数(ヘルパーマン)
    ============================ */

 inline bool isReal(const Value &v) { return std::holds_alternative<double>(v); }
 inline bool isComplex(const Value &v) { return std::holds_alternative<std::complex<double>>(v); }
 inline bool isDouble(const Value &v) { return std::holds_alternative<double>(v); }
 inline bool isInvalid(const Value &v) { return std::holds_alternative<InvalidValue>(v); }
 inline bool isString(const Value &v) { return std::holds_alternative<std::string>(v); }
 inline const std::string &asString(const Value &v, size_t pos) {
  if (!std::holds_alternative<std::string>(v)) { throw CalcError(CalcErrorType::TypeError, "expected string", pos); }
  return std::get<std::string>(v);
 }

 inline std::complex<double> toComplex(const Value &v) {
  if (std::holds_alternative<double>(v)) return std::complex<double>(std::get<double>(v), 0.0);
  return std::get<std::complex<double>>(v);
 }

 inline double asDouble(const Value &v, size_t pos = 0) {
  if (!isDouble(v)) throw CalcError(CalcErrorType::TypeError, errorMessage(CalcErrorType::TypeError), pos);
  return std::get<double>(v);
 }
 inline double asReal(const Value &v, size_t pos) {
  if (!std::holds_alternative<double>(v)) throw CalcError(CalcErrorType::TypeError, errorMessage(CalcErrorType::TypeError), pos);
  return std::get<double>(v);
 }
 inline Complex asComplex(const Value &v) {
  if (std::holds_alternative<Complex>(v)) return std::get<Complex>(v);
  return Complex(std::get<double>(v), 0.0);
 }
 inline int64_t asInt64(const Value &v, size_t pos) {
  if (isComplex(v)) throw CalcError(CalcErrorType::TypeError, "int required", pos);
  double x = std::get<double>(v);
  if (!std::isfinite(x)) throw CalcError(CalcErrorType::DomainError, "nan/inf", pos);
  if (std::floor(x) != x) throw CalcError(CalcErrorType::TypeError, "integer required", pos);
  if (x < (double)LLONG_MIN || x > (double)LLONG_MAX) throw CalcError(CalcErrorType::OutOfRange, "int64 overflow", pos);

  return (int64_t)x;
 }

 inline Value makeValue(double v) { return Value{v}; }
 inline Value makeValue(std::complex<double> v) { return Value{v}; }

 inline Value mulReal(const Value &v, double k, size_t pos) {
  if (isComplex(v)) return asComplex(v) * k;
  return asDouble(v, pos) * k;
 }

 inline Value divReal(const Value &v, double k, size_t pos) {
  if (k == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
  if (isComplex(v)) return asComplex(v) / k;
  return asDouble(v, pos) / k;
 }

 // 有限性チェックたち
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
       } else if constexpr (std::is_same_v<T, std::string>) {
        // string は有限チェック不要（タグ/オプション用途）
       } else if constexpr (std::is_same_v<T, InvalidValue>) {
        // Invalidは別ルートで処理される前提
       } else {
        static_assert(always_false<T>);
       }
      },
      v);
 }

 inline double requireFinite(double x, size_t pos) {
  if (!std::isfinite(x)) throwDomain(pos);
  return x;
 }

 // 警告ウマーマン
 inline void calcWarn(const FunctionContext &ctx, const std::string &msg) { std::cerr << "[WARN] pos=" << ctx.pos << " " << msg << "\n"; }
 inline void calcWarn(SystemConfig &cfg, size_t pos, const std::string &msg) { std::cerr << "[WARN] pos=" << pos << " " << msg << "\n"; }

 Value roundValue(const Value &v, const SystemConfig &cfg);

 inline double round12(double v) { return std::round(v / cnst_precision_inv) * cnst_precision_inv; }
} // namespace mm::cal
