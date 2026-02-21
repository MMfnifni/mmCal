#pragma once

#include "core.hpp"
namespace mm::cal {
 template <class Op> inline Value numericBinaryOp(const Value &a, const Value &b, size_t pos, Op op) {
  return std::visit(Overloaded{

                        // --------------------------
                        // scalar / complex
                        // --------------------------

                        [&](double x, double y) -> Value { return op(x, y); },

                        [&](double x, std::complex<double> y) -> Value { return op(std::complex<double>(x, 0.0), y); },

                        [&](std::complex<double> x, double y) -> Value { return op(x, std::complex<double>(y, 0.0)); },

                        [&](std::complex<double> x, std::complex<double> y) -> Value { return op(x, y); },

                        // --------------------------
                        // multi + scalar/complex
                        // --------------------------

                        [&](std::shared_ptr<MultiValue> m, double y) -> Value {
                         auto result = std::make_shared<MultiValue>();

                         for (const auto &e : m->elems_)
                          result->elems_.push_back(numericBinaryOp(e, Value(y), pos, op));

                         return result;
                        },

                        [&](std::shared_ptr<MultiValue> m, std::complex<double> y) -> Value {
                         auto result = std::make_shared<MultiValue>();

                         for (const auto &e : m->elems_)
                          result->elems_.push_back(numericBinaryOp(e, Value(y), pos, op));

                         return result;
                        },

                        // --------------------------
                        // scalar/complex + multi
                        // --------------------------

                        [&](double x, std::shared_ptr<MultiValue> m) -> Value {
                         auto result = std::make_shared<MultiValue>();

                         for (const auto &e : m->elems_)
                          result->elems_.push_back(numericBinaryOp(Value(x), e, pos, op));

                         return result;
                        },

                        [&](std::complex<double> x, std::shared_ptr<MultiValue> m) -> Value {
                         auto result = std::make_shared<MultiValue>();

                         for (const auto &e : m->elems_)
                          result->elems_.push_back(numericBinaryOp(Value(x), e, pos, op));

                         return result;
                        },

                        // --------------------------
                        // multi + multi
                        // --------------------------

                        [&](std::shared_ptr<MultiValue> m1, std::shared_ptr<MultiValue> m2) -> Value {
                         if (m1->elems_.size() != m2->elems_.size()) throw CalcError(CalcErrorType::DimensionMismatch, "multi-value size mismatch", pos);

                         auto result = std::make_shared<MultiValue>();
                         result->elems_.reserve(m1->elems_.size());

                         for (size_t i = 0; i < m1->elems_.size(); ++i)
                          result->elems_.push_back(numericBinaryOp(m1->elems_[i], m2->elems_[i], pos, op));

                         return result;
                        },

                        // --------------------------
                        // invalid
                        // --------------------------

                        [&](auto &&, auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "invalid operand types", pos); }

                    },
                    a.storage(), b.storage());
 }

 inline Value add(const Value &a, const Value &b, size_t pos) {
  return numericBinaryOp(a, b, pos, [](auto x, auto y) { return x + y; });
 }

 inline Value sub(const Value &a, const Value &b, size_t pos) {
  return numericBinaryOp(a, b, pos, [](auto x, auto y) { return x - y; });
 }

 inline Value mul(const Value &a, const Value &b, size_t pos) {
  return numericBinaryOp(a, b, pos, [](auto x, auto y) { return x * y; });
 }

 inline Value div(const Value &a, const Value &b, size_t pos) {
  return numericBinaryOp(a, b, pos, [&](auto x, auto y) {
   using Y = std::decay_t<decltype(y)>;

   if constexpr (std::is_same_v<Y, double>) {
    if (y == 0.0) throw CalcError(CalcErrorType::DivisionByZero, "division by zero", pos);
   } else if constexpr (std::is_same_v<Y, std::complex<double>>) {
    if (y == std::complex<double>(0.0, 0.0)) throw CalcError(CalcErrorType::DivisionByZero, "division by zero", pos);
   }

   return x / y;
  });
 }

 inline bool isInteger(double x) {
  double r = std::round(x);
  return std::abs(x - r) < 1e-12;
 }

 inline Value power(const Value &a, const Value &b, size_t pos) {
  if (!a.isNumeric() || !b.isNumeric()) throw CalcError(CalcErrorType::TypeError, "numeric required", pos);

  // ===== real Å~ real fast path =====
  if (a.isScalar() && b.isScalar()) {
   double base = a.asScalar(pos);
   double exp = b.asScalar(pos);

   // 0^0
   if (base == 0.0 && exp == 0.0) throwDomain(pos);

   // 0^negative
   if (base == 0.0) {
    if (exp > 0.0) return 0.0;
    throwDomain(pos);
   }

   // ïâêîÇÃîÒêÆêîèÊ Å® complex principal branch
   if (base < 0.0 && !isInteger(exp)) {
    calcWarn(pos, "(-a)^(p/q) : principal value only; other branches may exist");

    std::complex<double> z(base, 0.0);
    std::complex<double> w(exp, 0.0);
    auto r = std::pow(z, w);

    if (std::abs(r.imag()) < cnst_precision_inv) return r.real();

    return r;
   }

   // í èÌ real pow
   return std::pow(base, exp);
  }

  // ===== complex path =====
  std::complex<double> z = a.toComplex(pos);
  std::complex<double> w = b.toComplex(pos);

  if (std::abs(z) == 0.0 && std::abs(w) == 0.0) throw CalcError(CalcErrorType::DomainError, "0^0 is indeterminate", pos);

  std::complex<double> result = std::pow(z, w);

  if (std::abs(result.imag()) < cnst_precision_inv) return result.real();

  return result;
 }

 inline Value negate(const Value &v, size_t pos) {
  return std::visit(Overloaded{

                        [](double x) -> Value { return -x; },

                        [](std::complex<double> x) -> Value { return -x; },

                        [&](std::shared_ptr<MultiValue> m) -> Value {
                         auto result = std::make_shared<MultiValue>();
                         result->elems_.reserve(m->elems_.size());

                         for (const auto &e : m->elems_)
                          result->elems_.push_back(negate(e, pos));

                         return result;
                        },

                        [&](auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "invalid operand for unary -", pos); }

                    },
                    v.storage());
 }

 inline double factorialScalar(double x, size_t pos) {
  if (x < 0.0) throw CalcError(CalcErrorType::DomainError, "factorial undefined for negative values", pos);

  double intpart;
  if (std::modf(x, &intpart) != 0.0) throw CalcError(CalcErrorType::DomainError, "factorial requires integer value", pos);

  if (x > 170.0) // double overflow safety
   throw CalcError(CalcErrorType::Overflow, "factorial overflow", pos);

  double result = 1.0;
  for (int i = 1; i <= static_cast<int>(x); ++i)
   result *= i;

  return result;
 }

 inline Value factorial(const Value &v, size_t pos) {
  return std::visit(Overloaded{// double
                               [&](double x) -> Value { return factorialScalar(x, pos); },

                               // multi (vectorized)
                               [&](std::shared_ptr<MultiValue> m) -> Value {
                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(m->elems_.size());

                                for (const auto &e : m->elems_)
                                 result->elems_.push_back(factorial(e, pos));

                                return result;
                               },

                               // invalid
                               [&](auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "invalid operand for factorial", pos); }

                    },
                    v.storage());
 }
} // namespace mm::cal