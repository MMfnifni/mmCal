#pragma once

#include "core.hpp"
#include <algorithm>
#include <cmath>
#include <intrin.h>
#include <limits>
#include <numeric>
#include <vector>

namespace mm::cal {

 /* ============================
   数学補助系
   ============================ */
 inline constexpr double toDeg(double rad) { return rad * (180.0 / PI); }
 inline constexpr double toRad(double deg) { return deg * (PI / 180.0); }
 inline constexpr double radToDeg(double r) { return r * 180.0 / PI; }

 inline bool tolerantEqual(double a, double b, int precision) {
  double eps = std::pow(10.0, -precision);
  return std::abs(a - b) <= eps * std::max({1.0, std::abs(a), std::abs(b)});
 }

 bool eq(double a, double b);

 inline double inf(int sign = +1) { return sign >= 0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity(); }

 /* ============================
   整数系
   ============================ */

 inline long long requireInt(const Value &v, size_t pos) {
  double d = asDouble(v, pos);
  if (std::floor(d) != d) throw CalcError(CalcErrorType::NeedInteger, errorMessage(CalcErrorType::NeedInteger), pos);
  return static_cast<int64_t>(d);
 }
 inline long long checkedAdd(long long a, long long b, size_t pos) {
  if ((b > 0 && a > LLONG_MAX - b) || (b < 0 && a < LLONG_MIN - b)) throw CalcError(CalcErrorType::Overflow, errorMessage(CalcErrorType::Overflow), pos);
  return a + b;
 }
 inline long long checkedSub(long long a, long long b, size_t pos) {
  if ((b < 0 && a > LLONG_MAX + b) || (b > 0 && a < LLONG_MIN + b)) throw CalcError(CalcErrorType::Overflow, errorMessage(CalcErrorType::Overflow), pos);
  return a - b;
 }
 // 実験的にx64/x86ごとに実装してみる
 inline long long checkedMul(long long a, long long b, size_t pos) {
#if defined(_M_X64)
  long long hi;
  const long long lo = _mul128(a, b, &hi);
  if (hi != (lo >> 63)) throwOverflow(pos);
  return lo;
#else
  // x86 / ARM / etc: portable fallback
  if (a == 0 || b == 0) return 0;
  const bool overflow = (a > 0) ? ((b > 0) ? (a > LLONG_MAX / b) : (b < LLONG_MIN / a)) : ((b > 0) ? (a < LLONG_MIN / b) : (a < LLONG_MAX / b));
  if (overflow) throwOverflow(pos);
  return a * b;
#endif
 }
 inline long long checkedPow(long long base, long long exp, size_t pos) {
  if (exp < 0) throw CalcError(CalcErrorType::DomainError, "pow: negative exponent", pos);

  long long result = 1;
  long long x = base;

  while (exp) {
   if (exp & 1) result = checkedMul(result, x, pos);
   exp >>= 1;
   if (exp) x = checkedMul(x, x, pos);
  }
  return result;
 }
 inline constexpr long long gcdLL(long long a, long long b) noexcept { return std::gcd(a, b); }
 long long combLL(long long n, long long r, size_t pos);
 long long permLL(long long n, long long r, size_t pos);
 double factorial(double x, size_t pos);
 long double factLD(long long n, size_t pos);
 unsigned long long fibULL(unsigned long long n, size_t pos);
 bool isPrimeLL(long long n);

 /* ============================
   統計系
   ============================ */
 std::vector<double> gatherReals(const std::vector<Value> &v, size_t pos);
 std::vector<double> collectReals(const std::vector<Value> &v, FunctionContext &ctx);
 double meanOf(const std::vector<double> &a);
 double stddevPopulation(const std::vector<double> &a, double mu);
 double quantileLinear(std::vector<double> a, double p, size_t pos);
 std::vector<double> rankAverageTies(const std::vector<double> &x, size_t pos);
 double pearsonCorr(const std::vector<double> &x, const std::vector<double> &y, size_t pos);
 double medianOfSorted(const std::vector<double> &a);
 double variancePopulation(const std::vector<double> &a, double mu);

} // namespace mm::cal
