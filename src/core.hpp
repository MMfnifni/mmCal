#pragma once

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
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

 using Complex = std::complex<double>;

 struct ControlRequest : std::exception {
   enum class Kind { Exit, Clear };
   Kind kind;

   explicit ControlRequest(Kind k) : kind(k) {}
 };
 struct SystemConfig;
 struct FunctionContext;
 struct FunctionDef;
 struct InputEntry;
 struct MultiValue;
 struct SideEffect;
 struct EvalResult;
 struct EvaluationContext;
 struct ASTNode; // 実装はsyntax.hppへ
 struct UserFunction {
   std::vector<std::string> params;
   std::unique_ptr<ASTNode> body;
   UserFunction();
   ~UserFunction();
   UserFunction(UserFunction &&) noexcept;
   UserFunction &operator=(UserFunction &&) noexcept;
   UserFunction(const UserFunction &) = delete;
   UserFunction &operator=(const UserFunction &) = delete;
 };
 enum class AngleMode {
  Deg,
  Rad,
  Grad,
 };
 /* ============================
    エラー
    ============================ */

 enum class CalcErrorType {
  UnknownIdentifier,
  SyntaxError,
  TypeError,
  DomainError,
  DivisionByZero,
  OutOfRange,
  Overflow,
  NonConvergence,
  IOError,
  DefinitionError,
  InternalError,
 };

 constexpr const char *errorMessage(CalcErrorType t) {
  switch (t) {
   case CalcErrorType::UnknownIdentifier: return "unknown identifier";
   case CalcErrorType::SyntaxError: return "syntax error";
   case CalcErrorType::TypeError: return "type error";
   case CalcErrorType::DomainError: return "domain error";
   case CalcErrorType::DivisionByZero: return "division by zero";
   case CalcErrorType::OutOfRange: return "OutOfRange: history out of range";
   case CalcErrorType::Overflow: return "overflow detected";
   case CalcErrorType::NonConvergence: return "did not converge";
   case CalcErrorType::IOError: return "I/O error";
   case CalcErrorType::DefinitionError: return "definition error";
   case CalcErrorType::InternalError: return "internal error";
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
 [[noreturn]] inline void throwType(size_t pos) { throw CalcError(CalcErrorType::TypeError, errorMessage(CalcErrorType::TypeError), pos); }
 [[noreturn]] inline void throwType(size_t pos, const char *msg) { throw CalcError(CalcErrorType::TypeError, msg, pos); }
 [[noreturn]] inline void throwSyntax(size_t pos) { throw CalcError(CalcErrorType::SyntaxError, errorMessage(CalcErrorType::SyntaxError), pos); }
 [[noreturn]] inline void throwSyntax(size_t pos, const char *msg) { throw CalcError(CalcErrorType::SyntaxError, msg, pos); }
 [[noreturn]] inline void throwUnknown(size_t pos, const std::string &msg) { throw CalcError(CalcErrorType::UnknownIdentifier, msg, pos); }
 [[noreturn]] inline void throwInternal(size_t pos, const char *msg) { throw CalcError(CalcErrorType::InternalError, msg, pos); }

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
    if (!isScalar()) throw CalcError(CalcErrorType::TypeError, "TypeError: asScalar: expected scalar", pos);
    return std::get<double>(data_);
   }

   std::complex<double> asComplex(size_t pos) const {
    if (auto p = std::get_if<std::complex<double>>(&data_)) return *p;
    if (auto p = std::get_if<double>(&data_)) return std::complex<double>(*p, 0.0);
    throw CalcError(CalcErrorType::TypeError, "TypeError: asComplex: expected complex", pos);
   }

   const std::string &asString(size_t pos) const {
    if (!isString()) throw CalcError(CalcErrorType::TypeError, "TypeError: asString: expected string", pos);
    return std::get<std::string>(data_);
   }

   MultiPtr asMulti(size_t pos) const {
    if (!isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: asMulti: expected multivalue", pos);
    return std::get<MultiPtr>(data_);
   }
   const MultiValue &asMultiRef(size_t pos) const;

   // conversions
   std::complex<double> toComplex(size_t pos) const {
    if (isScalar()) return {std::get<double>(data_), 0.0};
    if (isComplex()) return std::get<std::complex<double>>(data_);

    throw CalcError(CalcErrorType::TypeError, "TypeError: toComplex: numeric required", pos);
   }

   int64_t toInt64(size_t pos) const {
    double x = asScalar(pos);

    if (!std::isfinite(x)) throw CalcError(CalcErrorType::DomainError, "DomainError: nan/inf", pos);
    if (std::floor(x) != x) throw CalcError(CalcErrorType::TypeError, "TypeError: integer required", pos);
    if (x < (double)LLONG_MIN || x > (double)LLONG_MAX) throw CalcError(CalcErrorType::OutOfRange, "OutOfRange: int64 overflow", pos);

    return static_cast<int64_t>(x);
   }

   // multi util
   std::size_t multiSize(size_t pos) const;
   bool multiEmpty(size_t pos) const;
   const Value &multiAt(std::size_t i, size_t pos) const;
   bool isEmptyMulti() const noexcept;

   template <class F> void forEachMulti(size_t pos, F &&f) const;

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

   bool hasNestedMulti(size_t pos) const;

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

   // visitor access
   const Storage &storage() const noexcept { return data_; }
   Storage &storage() noexcept { return data_; }
   explicit Value(const EvalResult &r);

   template <class Visitor> decltype(auto) visit(Visitor &&v) { return std::visit(std::forward<Visitor>(v), data_); }
   template <class Visitor> decltype(auto) visit(Visitor &&v) const { return std::visit(std::forward<Visitor>(v), data_); }
   std::size_t size() const;

   static Value None() { return Value(); } // monostate
   bool isNone() const { return std::holds_alternative<std::monostate>(data_); };
 };
 struct MultiValue {
   std::vector<Value> elems_;
   MultiValue() = default;

   explicit MultiValue(std::vector<Value> elems) : elems_(std::move(elems)) {}

   const std::vector<Value> &elems() const noexcept { return elems_; }

   std::size_t size() const noexcept { return elems_.size(); }

   const Value &operator[](std::size_t i) const { return elems_[i]; }

   bool empty() const noexcept { return elems_.empty(); }

   template <class F> void forEach(F &&f) const {
    for (const auto &v : elems_)
     f(v);
   }

   template <class F> Value map(F &&f) const;

   template <class F> Value mapMulti(size_t pos, F &&f) const {
    if (!Value::isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: multivalue required", pos);
    return Value::asMultiRef(pos).map(std::forward<F>(f));
   }

   auto begin() const noexcept { return elems_.begin(); }
   auto end() const noexcept { return elems_.end(); }

   bool hasNested() const noexcept;

   std::vector<std::vector<double>> toMatrix() const;
 };

 struct InputEntry {
   std::string expr; // In[n]
   Value value;      // Out[n]
 };

 struct EvalResult {
   Value value{};
   std::string displayOverride{};
   std::string explain{};
   bool suppressDisplay = false;

   bool hasDisplayOverride() const { return !displayOverride.empty(); }
 };

 struct SideEffect {
   enum class Kind { FileWrite, ClipboardCopy, Message, Explain, SuppressDisplay };

   Kind kind;
   std::string a; // path or text
   std::string b; // file content (for write)
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
   std::size_t maxDisplayElements = 10000;
   AngleMode angleMode = AngleMode::Deg;
   bool allowVariableRedefinition = true;
   bool allowFunctionRedefinition = false;
 };
 struct SessionState {
   SystemConfig &cfg;
   std::vector<InputEntry> &history;
   int &base;
   std::unordered_map<std::string, Value> &globals;
   std::unordered_map<std::string, UserFunction> &userFunctions;

   SessionState(SystemConfig &c, std::vector<InputEntry> &h, int &b, std::unordered_map<std::string, Value> &g, std::unordered_map<std::string, UserFunction> &uf) : cfg(c), history(h), base(b), globals(g), userFunctions(uf) {}
 };

 struct FunctionContext {
   SessionState &session;
   size_t pos;
   std::vector<SideEffect> &sideEffects;
 };

 struct EvaluationContext {
   SessionState &session;
   std::unordered_map<std::string, Value> locals;
   std::vector<SideEffect> sideEffects;
   std::vector<std::string> callStack;

   EvaluationContext(SessionState &s) : session(s) {}
 };

 inline bool hasVariable(const EvaluationContext &ctx, const std::string &name) { return ctx.locals.contains(name) || ctx.session.globals.contains(name); }

 inline const Value *findVariable(const EvaluationContext &ctx, const std::string &name) {
  if (auto it = ctx.locals.find(name); it != ctx.locals.end()) return &it->second;
  if (auto it = ctx.session.globals.find(name); it != ctx.session.globals.end()) return &it->second;
  return nullptr;
 }

 inline Value *findVariable(EvaluationContext &ctx, const std::string &name) {
  if (auto it = ctx.locals.find(name); it != ctx.locals.end()) return &it->second;
  if (auto it = ctx.session.globals.find(name); it != ctx.session.globals.end()) return &it->second;
  return nullptr;
 }

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

                        [&](auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "TypeError: numeric argument required", pos); }

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

                        [&](auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "TypeError: invalid operand for real multiplication", pos); }

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

                        [&](auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "TypeError: invalid operand for real division", pos); }

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

 inline std::size_t Value::size() const {
  if (isMulti()) { return std::get<MultiPtr>(data_)->size(); }
  if (isString()) { return std::get<std::string>(data_).size(); }
  if (isInvalid()) { return 0; }
  // scalar / complex → 要素1とみなす
  if (isNumeric()) { return 1; }
  throw CalcError(CalcErrorType::InternalError, "RuntimeError: Value->size()???", 0);
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
 inline void calcWarn(const std::string &msg) { std::cerr << "[WARN] " << msg << "\n"; }

 Value roundValue(const Value &v, const SystemConfig &cfg);

 std::vector<std::vector<double>> toMatrix(const Value &val, const FunctionContext &ctx);

 // その他ユーティリティ
 std::string whatByteUnit(size_t filesize);
 void serializeValueImpl(const Value &v, std::ostream &os, const SystemConfig &cfg, int depth);
 inline void serializeValue(const Value &v, std::ostream &os, const SystemConfig &cfg) { serializeValueImpl(v, os, cfg, 0); }
 void applySideEffects(EvaluationContext &ectx, EvalResult &result);
 // 内部的にイコールとする
 inline constexpr double numeric_epsilon() noexcept { return cnst_precision_inv; }
 inline bool nearly_equal(double a, double b, double eps = cnst_precision_inv) noexcept { return std::abs(a - b) <= eps * std::max({1.0, std::abs(a), std::abs(b)}); }
 inline bool nearly_zero(double x, double eps = cnst_precision_inv) noexcept { return std::abs(x) <= eps; }
 inline bool nearly_equal(std::complex<double> a, std::complex<double> b, double eps = cnst_precision_inv) noexcept { return nearly_equal(a.real(), b.real(), eps) && nearly_equal(a.imag(), b.imag(), eps); }
 inline bool nearly_zero(std::complex<double> z, double eps = cnst_precision_inv) noexcept { return nearly_zero(z.real(), eps) && nearly_zero(z.imag(), eps); }
 inline bool tolerantEqual(double a, double b, int /*precision*/) noexcept { return nearly_equal(a, b); }          // 互換用
 inline bool tolerantEqual(std::complex<double> a, std::complex<double> b) noexcept { return nearly_equal(a, b); } // 互換用

 template <class F> inline void Value::forEachMulti(size_t pos, F &&f) const {
  if (!isMulti()) return;
  const auto &mv = asMultiRef(pos);
  for (const auto &v : mv.elems())
   f(v);
 }

} // namespace mm::cal
