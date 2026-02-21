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
 struct MultiValue;

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
  DimensionMismatch,
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
   case CalcErrorType::FunctionMissing: return "function argument missing(many or few args)";
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
   case CalcErrorType::DimensionMismatch: return "Dimension mismatch";
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
 //[[noreturn]] inline void throwRuntime(const char *msg) { throw CalcError(CalcErrorType::RuntimeError, msg, 0); }

 /* ============================
   値だぞっ
   ============================ */
 class Value {
  public:
   using MultiPtr = std::shared_ptr<MultiValue>;

  private:
   using Storage = std::variant<std::monostate, double, std::complex<double>, std::string, MultiPtr>;

   Storage data_;

  public:
   // constructors
   Value() noexcept : data_(std::monostate{}) {}
   Value(double v) : data_(v) {}
   Value(std::complex<double> v) : data_(v) {}
   Value(const std::string &s) : data_(s) {}
   Value(std::string &&s) : data_(std::move(s)) {}
   Value(MultiPtr mv) : data_(std::move(mv)) {}

   // type checks
   bool isInvalid() const noexcept { return std::holds_alternative<std::monostate>(data_); }
   bool isScalar() const noexcept { return std::holds_alternative<double>(data_); }
   bool isComplex() const noexcept { return std::holds_alternative<std::complex<double>>(data_); }
   bool isString() const noexcept { return std::holds_alternative<std::string>(data_); }
   bool isMulti() const noexcept { return std::holds_alternative<MultiPtr>(data_); }
   bool isNumeric() const noexcept { return isScalar() || isComplex(); }

   // getters
   double asScalar(size_t pos) const {
    if (!isScalar()) throw CalcError(CalcErrorType::TypeError, "asScalar: expected scalar", pos);
    return std::get<double>(data_);
   }

   std::complex<double> asComplex(size_t pos) const {
    if (auto p = std::get_if<std::complex<double>>(&data_)) return *p;
    if (auto p = std::get_if<double>(&data_)) return std::complex<double>(*p, 0.0);
    throw CalcError(CalcErrorType::TypeError, "asComplex: expected complex", pos);
   }

   const std::string &asString(size_t pos) const {
    if (!isString()) throw CalcError(CalcErrorType::TypeError, "asString: expected string", pos);
    return std::get<std::string>(data_);
   }

   MultiPtr asMulti(size_t pos) const {
    if (!isMulti()) throw CalcError(CalcErrorType::TypeError, "asMulti: expected multivalue", pos);
    return std::get<MultiPtr>(data_);
   }

   // conversions
   std::complex<double> toComplex(size_t pos) const {
    if (isScalar()) return {std::get<double>(data_), 0.0};
    if (isComplex()) return std::get<std::complex<double>>(data_);

    throw CalcError(CalcErrorType::TypeError, "toComplex: numeric required", pos);
   }

   int64_t toInt64(size_t pos) const {
    double x = asScalar(pos);

    if (!std::isfinite(x)) throw CalcError(CalcErrorType::DomainError, "nan/inf", pos);
    if (std::floor(x) != x) throw CalcError(CalcErrorType::TypeError, "integer required", pos);
    if (x < (double)LLONG_MIN || x > (double)LLONG_MAX) throw CalcError(CalcErrorType::OutOfRange, "int64 overflow", pos);

    return static_cast<int64_t>(x);
   }

  private:
   static inline void checkFiniteImpl(const std::monostate &, size_t) {
    // ignore
   }
   static inline void checkFiniteImpl(double x, size_t pos) {
    if (!std::isfinite(x)) throwDomain(pos);
   }
   static void checkFiniteImpl(const std::complex<double> &z, size_t pos) {
    if (!std::isfinite(z.real()) || !std::isfinite(z.imag())) throwDomain(pos);
   }
   static void checkFiniteImpl(const std::string &, size_t) {
    // string has no numeric domain将来
   }
   static void checkFiniteImpl(const MultiPtr &mv, size_t pos);

  public:
   // Finite guarantee layer

   static void checkFinite(const Value &v, size_t pos) {
    std::visit([&](auto &&x) { checkFiniteImpl(x, pos); }, v.data_);
   }
   static Value requireFinite(Value v, size_t pos) {
    checkFinite(v, pos);
    return v;
   }
   // finite
   void validateFinite(size_t pos) const;

  public:
   // visitor access
   const Storage &storage() const noexcept { return data_; }
   Storage &storage() noexcept { return data_; }

   template <class Visitor> decltype(auto) visit(Visitor &&v) { return std::visit(std::forward<Visitor>(v), data_); }
   template <class Visitor> decltype(auto) visit(Visitor &&v) const { return std::visit(std::forward<Visitor>(v), data_); }
 };

 struct MultiValue {
   std::vector<Value> elems_;
   MultiValue() = default;

   explicit MultiValue(std::vector<Value> elems) : elems_(std::move(elems)) {}

   const std::vector<Value> &elems() const noexcept { return elems_; }

   std::size_t size() const noexcept { return elems_.size(); }

   const Value &operator[](std::size_t i) const { return elems_[i]; }

   auto begin() const noexcept { return elems_.begin(); }
   auto end() const noexcept { return elems_.end(); }
 };
 // 以下のほうが高速になる。可読性とのトレードオフ
 // struct Value {
 //  enum { Real, Complex } type;
 //  double r;
 //  std::complex<double> c;
 //};
 //

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
    補助関数(ヘルパーマン)
    ============================ */

 template <class... Ts> struct Overloaded : Ts... {
   using Ts::operator()...;
 };

 template <class... Ts> Overloaded(Ts...) -> Overloaded<Ts...>;

 template <typename RealFn, typename ComplexFn> inline Value applyUnaryNumeric(const Value &v, RealFn rf, ComplexFn cf, size_t pos) {
  return std::visit(Overloaded{

                        [&](double d) -> Value { return rf(d); },

                        [&](std::complex<double> c) -> Value { return cf(c); },

                        [&](std::shared_ptr<MultiValue> m) -> Value {
                         auto result = std::make_shared<MultiValue>();
                         result->elems_.reserve(m->elems_.size());

                         for (const auto &e : m->elems_)
                          result->elems_.push_back(applyUnaryNumeric(e, rf, cf, pos));

                         return result;
                        },

                        [&](auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "numeric argument required", pos); }

                    },
                    v.storage());
 }

 inline bool isScaler(const Value &v) { return v.isScalar(); }
 inline bool isComplex(const Value &v) { return v.isComplex(); }
 inline bool isDouble(const Value &v) { return v.isScalar(); }
 inline bool isInvalid(const Value &v) { return v.isInvalid(); }
 inline bool isString(const Value &v) { return v.isString(); }
 inline bool isMulti(const Value &v) { return v.isMulti(); }

 inline Value makeValue(double v) { return Value{v}; }
 inline Value makeValue(std::complex<double> v) { return Value{v}; }

 inline Value mulReal(const Value &v, double k, size_t pos) {
  return std::visit(Overloaded{

                        [&](double d) -> Value { return d * k; },

                        [&](std::complex<double> c) -> Value { return c * k; },

                        [&](std::shared_ptr<MultiValue> m) -> Value {
                         auto result = std::make_shared<MultiValue>();
                         result->elems_.reserve(m->elems_.size());

                         for (const auto &e : m->elems_)
                          result->elems_.push_back(mulReal(e, k, pos));

                         return result;
                        },

                        [&](auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "invalid operand for real multiplication", pos); }

                    },
                    v.storage());
 }

 inline Value divReal(const Value &v, double k, size_t pos) {
  if (k == 0.0) throw CalcError(CalcErrorType::DivisionByZero, "division by zero", pos);

  return std::visit(Overloaded{

                        [&](double d) -> Value { return d / k; },

                        [&](std::complex<double> c) -> Value { return c / k; },

                        [&](std::shared_ptr<MultiValue> m) -> Value {
                         auto result = std::make_shared<MultiValue>();
                         result->elems_.reserve(m->elems_.size());

                         for (const auto &e : m->elems_)
                          result->elems_.push_back(divReal(e, k, pos));

                         return result;
                        },

                        [&](auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "invalid operand for real division", pos); }

                    },
                    v.storage());
 }

 // Finite
 inline void Value::validateFinite(size_t pos) const {
  std::visit(
      [&](auto &&x) {
       using T = std::decay_t<decltype(x)>;

       if constexpr (std::is_same_v<T, std::monostate>) {
        return;
       } else if constexpr (std::is_same_v<T, double>) {
        checkFinite(x, pos);
       } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        checkFinite(x.real(), pos);
        checkFinite(x.imag(), pos);
       } else if constexpr (std::is_same_v<T, std::string>) {
        return;
       } else if constexpr (std::is_same_v<T, MultiPtr>) {
        if (!x) return;
        for (const auto &v : x->elems())
         v.validateFinite(pos);
       }
      },
      data_);

  // MultiValue 拡張時もここに追加すること(絶対未来の私は忘れてる)
 }

 // 有限性チェックたち
 inline void checkFinite(const std::complex<double> &z, size_t pos) { Value::checkFinite(Value(z), pos); }

 inline double requireFinite(double x, size_t pos) {
  if (!std::isfinite(x)) throwDomain(pos);
  return x;
 }

 inline void checkFinite(double x, size_t pos) {
  if (!std::isfinite(x)) throwDomain(pos);
 }

 inline void checkFinite(const Value &v, size_t pos) { Value::checkFinite(v, pos); }

 // 警告ウマーマン
 inline void calcWarn(const FunctionContext &ctx, const std::string &msg) { std::cerr << "[WARN] pos=" << ctx.pos << " " << msg << "\n"; }
 inline void calcWarn(SystemConfig &cfg, size_t pos, const std::string &msg) { std::cerr << "[WARN] pos=" << pos << " " << msg << "\n"; }
 inline void calcWarn(size_t pos, const std::string &msg) { std::cerr << "[WARN] pos=" << pos << " " << msg << "\n"; }

 Value roundValue(const Value &v, const SystemConfig &cfg);

} // namespace mm::cal
