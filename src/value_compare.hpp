#pragma once

#include "core.hpp"
#include "value_arithmetic.hpp"

namespace mm::cal {

 enum class CompareOp { Eq, Ne, Lt, Le, Gt, Ge };

 // ---------------------------------------------
 // 内部的にイコールとする
 // ---------------------------------------------

 inline bool tolerantEqual(double a, double b) { return std::abs(a - b) <= cnst_precision_inv; }

 inline bool tolerantEqual(std::complex<double> a, std::complex<double> b) { return tolerantEqual(a.real(), b.real()) && tolerantEqual(a.imag(), b.imag()); }

 // ---------------------------------------------
 // 数値比較のdispatcher
 // ---------------------------------------------

 inline bool compareNumeric(double x, double y, CompareOp op) {
  switch (op) {
   case CompareOp::Eq: return tolerantEqual(x, y);
   case CompareOp::Ne: return !tolerantEqual(x, y);
   case CompareOp::Lt: return x < y && !tolerantEqual(x, y);
   case CompareOp::Le: return x < y || tolerantEqual(x, y);
   case CompareOp::Gt: return x > y && !tolerantEqual(x, y);
   case CompareOp::Ge: return x > y || tolerantEqual(x, y);
  }
  return false;
 }

 inline bool compareNumeric(std::complex<double> x, std::complex<double> y, CompareOp op, size_t pos) {
  switch (op) {

   case CompareOp::Eq: return tolerantEqual(x, y);

   case CompareOp::Ne: return !tolerantEqual(x, y);

   case CompareOp::Lt:
   case CompareOp::Le:
   case CompareOp::Gt:
   case CompareOp::Ge: throw CalcError(CalcErrorType::TypeError, "ordering comparison not defined for complex numbers", pos);
  }

  return false;
 }

 // ---------------------------------------------
 // main
 // ---------------------------------------------
 inline Value compare(const Value &a, const Value &b, CompareOp op, size_t pos) {

  auto result = std::visit(Overloaded{

                               [&](double x, double y) -> bool { return compareNumeric(x, y, op); },

                               [&](std::complex<double> x, std::complex<double> y) -> bool { return compareNumeric(x, y, op, pos); },

                               [&](double x, std::complex<double> y) -> bool { return compareNumeric(std::complex<double>(x, 0), y, op, pos); },

                               [&](std::complex<double> x, double y) -> bool { return compareNumeric(x, std::complex<double>(y, 0), op, pos); },

                               [&](const std::string &x, const std::string &y) -> bool {
                                switch (op) {
                                 case CompareOp::Eq: return x == y;
                                 case CompareOp::Ne: return x != y;
                                 case CompareOp::Lt: return x < y;
                                 case CompareOp::Le: return x <= y;
                                 case CompareOp::Gt: return x > y;
                                 case CompareOp::Ge: return x >= y;
                                }
                                return false;
                               },

                               [&](auto &&, auto &&) -> bool { throw CalcError(CalcErrorType::TypeError, "invalid operand types for comparison", pos); }

                           },
                           a.storage(), b.storage());

  return result ? Value(1.0) : Value(0.0);
 }

} // namespace mm::cal
