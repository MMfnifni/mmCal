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

 static constexpr double deg2rad = PI / 180.0;
 static constexpr double rad2deg = 180.0 / PI;
 inline constexpr double eps = cnst_precision_inv;

 inline constexpr double toDeg(double rad) { return rad * (180.0 / PI); }
 inline constexpr double toRad(double deg) { return deg * (PI / 180.0); }
 inline constexpr double radToDeg(double r) { return r * 180.0 / PI; }

 inline bool tolerantEqual(double a, double b, int precision) {
  double eps = std::pow(10.0, -precision);
  return std::abs(a - b) <= eps * std::max({1.0, std::abs(a), std::abs(b)});
 }

 inline double inf(int sign = +1) { return sign >= 0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity(); }

 // Complex(x, 0.0)
 static inline Complex cx(double x) { return {x, 0.0}; }

 // 「1/denom」を安全に返す関数
 Value safeInv(double denom, double sign);
 static inline Value signedInfBy(double sign) { return inf(sign >= 0 ? +1 : -1); }
 // denom が 0 に近いときに ±inf を返す（符号は signSource から取る）
 static inline Value invOrSignedInf(double denom, double signSource) { return (std::abs(denom) < eps) ? signedInfBy(signSource) : (1.0 / denom); }

 bool eq(double a, double b);

 // -------------------------------
 // 基本昇格
 // -------------------------------

 inline Complex toComplex(const Value &v, size_t pos) {
  return std::visit(Overloaded{[](double d) { return Complex(d, 0.0); }, [](Complex c) { return c; }, [&](auto &&) -> Complex { throw CalcError(CalcErrorType::TypeError, "numeric value required", pos); }}, v.storage());
 }

 inline bool isZero(Complex z) { return z.real() == 0.0 && z.imag() == 0.0; }

 // -------------------------------
 // 実数限定取得
 // -------------------------------

 inline double requireReal(const Value &v, size_t pos) {
  return std::visit(Overloaded{[](double d) { return d; },
                               [&](Complex c) -> double {
                                if (std::abs(c.imag()) > cnst_precision_inv) throw CalcError(CalcErrorType::DomainError, "real value required", pos);
                                return c.real();
                               },
                               [&](auto &&) -> double { throw CalcError(CalcErrorType::TypeError, "real value required", pos); }},
                    v.storage());
 }

 // -------------------------------
 // realに戻す
 // -------------------------------

 inline Value realIfPossible(Complex z) {
  if (std::abs(z.imag()) < cnst_precision_inv) return Value(z.real());
  return Value(z);
 }

 // -------------------------------
 // 単項適用
 // -------------------------------

 template <class F> inline Value applyUnaryNumeric(const Value &v, size_t pos, F f) {
  return std::visit(
      Overloaded{[&](double d) -> Value { return realIfPossible(f(Complex(d, 0.0))); }, [&](Complex c) -> Value { return realIfPossible(f(c)); }, [&](auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "numeric operand required", pos); }},
      v.storage());
 }

 // -------------------------------
 // 二項適用
 // -------------------------------

 template <class F> inline Value applyBinaryNumeric(const Value &a, const Value &b, size_t pos, F f) {
  return std::visit(Overloaded{[&](double x, double y) -> Value { return realIfPossible(f(Complex(x, 0), Complex(y, 0))); }, [&](double x, Complex y) -> Value { return realIfPossible(f(Complex(x, 0), y)); },
                               [&](Complex x, double y) -> Value { return realIfPossible(f(x, Complex(y, 0))); }, [&](Complex x, Complex y) -> Value { return realIfPossible(f(x, y)); },
                               [&](auto &&, auto &&) -> Value { throw CalcError(CalcErrorType::TypeError, "numeric operands required", pos); }},
                    a.storage(), b.storage());
 }

 /* ============================
   整数系
   ============================ */

 inline long long requireInt(const Value &v, size_t pos) {
  double d = v.asScalar(pos);
  if (std::floor(d) != d) throw CalcError(CalcErrorType::NeedInteger, errorMessage(CalcErrorType::NeedInteger), pos);
  return static_cast<int64_t>(d);
 }
 inline std::complex<double> requireComplex(const Value &v, size_t pos) { return v.asComplex(pos); }
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
 inline double varianceRealWelford(const std::vector<Value> &v, FunctionContext &ctx, int ddof) {
  const int n = (int)v.size();
  if (n <= ddof) throw CalcError(CalcErrorType::DomainError, "var: too few elements", ctx.pos);
  if (n == 1) return 0.0;

  long double mean = 0.0L;
  long double m2 = 0.0L;
  int k = 0;

  for (auto &x : v) {
   const long double a = (long double)x.asScalar(ctx.pos);
   ++k;
   const long double d = a - mean;
   mean += d / k;
   m2 += d * (a - mean);
  }

  return (double)(m2 / (long double)(n - ddof));
 }
 // -----------------------------
 // Welfordで mean + M2 を返す変態設計
 // -----------------------------
 struct MeanM2 {
   long double mean = 0.0L;
   long double m2 = 0.0L;
   int n = 0;
 };

 inline MeanM2 welfordMeanM2Real(const std::vector<Value> &v, FunctionContext &ctx) {
  MeanM2 s;
  for (const auto &x : v) {
   const long double a = static_cast<long double>(requireReal(x, ctx.pos));

   ++s.n;
   const long double d = a - s.mean;
   s.mean += d / static_cast<long double>(s.n);
   s.m2 += d * (a - s.mean);
  }
  return s;
 }

 inline double medianInplace(std::vector<double> &a) {
  const size_t n = a.size();
  if (n == 0) return 0.0;
  const size_t mid = n / 2;
  std::nth_element(a.begin(), a.begin() + mid, a.end());
  const double hi = a[mid];
  if (n & 1) return hi;
  const double lo = *std::max_element(a.begin(), a.begin() + mid);
  return 0.5 * (lo + hi);
 }

 inline double skewnessReal(const std::vector<Value> &v, FunctionContext &ctx) {
  const size_t n = v.size();
  if (n == 0) throw CalcError(CalcErrorType::DomainError, "skew: no elements", ctx.pos);

  long double mean = 0.0L;
  long double m2 = 0.0L;
  long double m3 = 0.0L;
  size_t k = 0;

  for (const auto &xv : v) {

   const long double x = static_cast<long double>(requireReal(xv, ctx.pos));

   const long double k1 = static_cast<long double>(k);

   ++k;

   const long double delta = x - mean;
   const long double delta_n = delta / k;
   const long double delta_n2 = delta_n * delta_n;
   const long double term1 = delta * delta_n * k1;

   mean += delta_n;
   m3 += term1 * delta_n * (k1 - 1.0L) - 3.0L * delta_n * m2;
   m2 += term1;
  }

  if (m2 == 0.0L) throwDomain(ctx.pos);

  return static_cast<double>(m3 / std::pow(m2, 1.5L));
 }

 inline double kurtosisExcessReal(const std::vector<Value> &v, FunctionContext &ctx) {
  const size_t n = v.size();
  if (n == 0) throw CalcError(CalcErrorType::DomainError, "kurt: no elements", ctx.pos);

  long double mean = 0.0L;

  for (const auto &xv : v)
   mean += static_cast<long double>(requireReal(xv, ctx.pos));

  mean /= static_cast<long double>(n);

  long double m2 = 0.0L;
  long double m4 = 0.0L;

  for (const auto &xv : v) {

   const long double x = static_cast<long double>(requireReal(xv, ctx.pos));

   const long double d = x - mean;
   const long double d2 = d * d;

   m2 += d2;
   m4 += d2 * d2;
  }

  if (m2 == 0.0L) throwDomain(ctx.pos);

  const long double nld = static_cast<long double>(n);

  return static_cast<double>(nld * m4 / (m2 * m2) - 3.0L);
 }

 inline double kurtosisExcessSample(const std::vector<Value> &v, FunctionContext &ctx) {
  const size_t n = v.size();
  if (n < 4) throw CalcError(CalcErrorType::DomainError, "kurts: need at least 4 elements", ctx.pos);

  long double mean = 0.0L;
  long double m2 = 0.0L, m3 = 0.0L, m4 = 0.0L;
  size_t k = 0;

  for (const auto &xv : v) {

   const long double x = static_cast<long double>(requireReal(xv, ctx.pos));

   const long double k1 = static_cast<long double>(k);

   ++k;

   const long double delta = x - mean;
   const long double delta_n = delta / k;
   const long double delta_n2 = delta_n * delta_n;
   const long double term1 = delta * delta_n * k1;

   mean += delta_n;

   m4 += term1 * delta_n2 * (k1 * k1 - 3.0L * k1 + 3.0L) + 6.0L * delta_n2 * m2 - 4.0L * delta_n * m3;

   m3 += term1 * delta_n * (k1 - 1.0L) - 3.0L * delta_n * m2;

   m2 += term1;
  }

  if (m2 == 0.0L) throwDomain(ctx.pos);

  const long double g2 = m4 / (m2 * m2) - 3.0L;

  const long double nn = static_cast<long double>(n);

  return static_cast<double>(((nn - 1.0L) / ((nn - 2.0L) * (nn - 3.0L))) * ((nn + 1.0L) * g2 + 6.0L));
 }

 inline double quantileLinearSorted(const std::vector<double> &a, double p01, size_t pos) {
  if (!(p01 >= 0.0 && p01 <= 1.0)) throw CalcError(CalcErrorType::DomainError, "quantile: p out of range", pos);
  if (a.empty()) throw CalcError(CalcErrorType::DomainError, "quantile: no samples", pos);

  // ソート済み前提
  if (a.size() == 1) return a[0];

  const double idx = p01 * (double)(a.size() - 1);
  const size_t i0 = (size_t)idx;
  const size_t i1 = std::min(i0 + 1, a.size() - 1);
  const double t = idx - (double)i0;
  return a[i0] * (1.0 - t) + a[i1] * t;
 }

 inline double covariancePopulationReal(const std::vector<Value> &v, FunctionContext &ctx) {
  const int total = static_cast<int>(v.size());

  if (total % 2 != 0) throw CalcError(CalcErrorType::DomainError, "cov: invalid argument count", ctx.pos);

  const int n = total / 2;
  if (n < 2) throw CalcError(CalcErrorType::DomainError, "cov: too few samples", ctx.pos);

  long double mx = 0.0L, my = 0.0L;
  long double c = 0.0L;
  int k = 0;

  for (int i = 0; i < n; ++i) {

   const long double x = static_cast<long double>(requireReal(v[i], ctx.pos));

   const long double y = static_cast<long double>(requireReal(v[i + n], ctx.pos));

   ++k;

   const long double dx = x - mx;
   const long double dy = y - my;

   mx += dx / k;
   my += dy / k;

   c += dx * (y - my);
  }

  return static_cast<double>(c / static_cast<long double>(n));
 }

 static double corrPopulationReal(const std::vector<Value> &v, FunctionContext &ctx) {
  const int total = (int)v.size();
  if (total % 2 != 0) throw CalcError(CalcErrorType::InvalidArgument, "corr: argument count must be even", ctx.pos);

  const int n = total / 2;
  if (n < 2) throw CalcError(CalcErrorType::DomainError, "corr: too few samples", ctx.pos);

  long double mx = 0.0L, my = 0.0L;
  long double sxx = 0.0L, syy = 0.0L, sxy = 0.0L;
  int k = 0;

  for (int i = 0; i < n; ++i) {
   const long double x = static_cast<long double>(requireReal(v[i], ctx.pos));

   const long double y = static_cast<long double>(requireReal(v[i + n], ctx.pos));

   ++k;
   const long double dx = x - mx;
   const long double dy = y - my;

   mx += dx / k;
   my += dy / k;

   sxx += dx * (x - mx);
   syy += dy * (y - my);
   sxy += dx * (y - my);
  }

  if (sxx == 0.0L || syy == 0.0L) throw CalcError(CalcErrorType::DomainError, "corr: zero variance", ctx.pos);

  return static_cast<double>(sxy / std::sqrt(sxx * syy));
 }
 Value areaPolygon(const std::vector<Value> &v, FunctionContext &ctx);

 inline constexpr long long gcdLL(long long a, long long b) noexcept { return std::gcd(a, b); }
 long long combLL(long long n, long long r, size_t pos);
 long long permLL(long long n, long long r, size_t pos);

 double factorial(double x, size_t pos);
 long double factLD(long long n, size_t pos);
 unsigned long long fibULL(unsigned long long n, size_t pos);
 inline bool isPrimeLL(long long n) {
  if (n <= 1) return false;                   // 1以下は素数ではない
  if (n <= 3) return true;                    // 2と3は素数
  if (n % 2 == 0 || n % 3 == 0) return false; // 2または3で割り切れる数は素数ではない

  // 5から始めて、iを6ずつ増やしていく
  for (long long i = 5; i * i <= n; i += 6) {
   if (n % i == 0 || n % (i + 2) == 0) return false; // iかi+2で割り切れる場合は素数ではない
  }

  return true; // それ以外は素数
 }
 inline Value fn_log(const std::vector<Value> &v, FunctionContext &ctx) {
  if (v.size() == 1) {
   return applyUnaryNumeric(v[0], ctx.pos, [&](Complex x) {
    if (x == Complex(0, 0)) throwDomain(ctx.pos);
    return std::log(x);
   });
  }

  return applyBinaryNumeric(v[1], v[0], ctx.pos, [&](Complex x, Complex base) {
   if (base == Complex(1, 0)) throwDomain(ctx.pos);
   return std::log(x) / std::log(base);
  });
 }

 inline Complex sin_over_t_series(Complex t) {
  // sin(t)/t = 1 - t^2/6 + t^4/120 - t^6/5040 + t^8/362880
  const Complex t2 = t * t;
  const Complex t4 = t2 * t2;
  const Complex t6 = t4 * t2;
  const Complex t8 = t4 * t4;
  return Complex(1, 0) - t2 / 6.0 + t4 / 120.0 - t6 / 5040.0 + t8 / 362880.0;
 }

 inline Complex tan_over_t_series(Complex t) {
  // tan(t)/t = 1 + t^2/3 + 2 t^4/15 + 17 t^6/315 + 62 t^8/2835
  const Complex t2 = t * t;
  const Complex t4 = t2 * t2;
  const Complex t6 = t4 * t2;
  const Complex t8 = t4 * t4;
  return Complex(1, 0) + t2 / 3.0 + (Complex(2, 0) * t4) / 15.0 + (Complex(17, 0) * t6) / 315.0 + (Complex(62, 0) * t8) / 2835.0;
 }

 inline double absC(Complex z) { return std::abs(z); } // 読みやすさ用

 inline FuncImpl makeSincLike(double scale) {
  return [=](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   const Complex x = toComplex(v[0], ctx.pos);

   if (x == Complex(0, 0)) return 1.0;

   const Complex t = x * scale;

   const Complex r = (std::abs(t) < 1e-8) ? (sin_over_t_series(t) * scale) : (std::sin(t) / x);

   return realIfPossible(r);
  };
 }

 inline FuncImpl makeCoscLike(double scale) {
  return [=](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   const Complex x = toComplex(v[0], ctx.pos);

   if (x == Complex(0, 0)) return 0.0;

   const Complex t = x * scale;
   const Complex h = t * 0.5;

   const Complex s = std::sin(h);
   const Complex r = (Complex(2, 0) * s * s) / x;

   return realIfPossible(r);
  };
 }

 inline FuncImpl makeTancLike(double scale) {
  return [=](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   const Complex x = toComplex(v[0], ctx.pos);

   if (x == Complex(0, 0)) return 1.0;

   const Complex t = x * scale;

   if (std::abs(t) < 1e-8) return realIfPossible(tan_over_t_series(t) * scale);

   const Complex c = std::cos(t);

   if (std::abs(c) < cnst_precision_inv) {
    const Complex s = std::sin(t);
    return inf(std::real(s) >= 0 ? +1 : -1);
   }

   const Complex r = std::tan(t) / x;

   return realIfPossible(r);
  };
 }

 inline FuncImpl makeDivXCmplxReal(std::function<Complex(Complex)> fc, std::function<double(double)> fr, double x0Value) {
  return [=](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   const Value &arg = v[0];

   if (arg.isComplex()) {
    Complex x = toComplex(arg, ctx.pos);

    if (x == Complex(0, 0)) return x0Value;

    return realIfPossible(fc(x) / x);
   } else {
    double x = requireReal(arg, ctx.pos);

    if (x == 0.0) return x0Value;

    return fr(x) / x;
   }
  };
 }

 inline FuncImpl makeExpc() {
  return [](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   const Value &arg = v[0];

   if (arg.isComplex()) {
    Complex x = toComplex(arg, ctx.pos);

    if (x == Complex(0, 0)) return 1.0;

    double ax = std::abs(x);

    if (ax < 1e-8) {
     Complex x2 = x * x;
     Complex x3 = x2 * x;
     Complex x4 = x2 * x2;

     return realIfPossible(Complex(1, 0) + x / 2.0 + x2 / 6.0 + x3 / 24.0 + x4 / 120.0);
    }

    return realIfPossible((std::exp(x) - Complex(1, 0)) / x);
   } else {
    double x = requireReal(arg, ctx.pos);

    if (x == 0.0) return 1.0;

    return std::expm1(x) / x;
   }
  };
 }

 inline long long absLL(long long x, size_t pos) {
  if (x == LLONG_MIN) throwOverflow(pos);
  return std::llabs(x);
 }

 //// IRR計算用ヘルパー
 // double fin_irr = [](const std::vector<double> &cf, auto &ctx) -> double {
 //  if (cf.size() < 2) throwDomain(ctx.pos);
 //  double x0 = 0.1; // 初期利率
 //  double eps = 1e-10;
 //  int maxIter = 200;

 // for (int i = 0; i < maxIter; ++i) {
 //  double f = 0, df = 0;
 //  for (size_t t = 0; t < cf.size(); ++t) {
 //   f += cf[t] / std::pow(1 + x0, t);
 //   if (t > 0) df -= t * cf[t] / std::pow(1 + x0, t + 1);
 //  }
 //  if (std::abs(f) < eps) return x0;
 //  if (df == 0) break; // ゼロ除算防止
 //  x0 -= f / df;
 // }
 // throwDomain(ctx.pos); // 収束しない場合
 //};

 //  IRR / RATE (Newton-Raphson + 二分法)
 // 精度が悪いので廃止(早いんだけどね)
 // auto newton_raphson = [&](auto f, auto df, double guess, double tol = 1e-10, int maxIter = 200) -> double {
 // double x = guess;
 // for (int i = 0; i < maxIter; ++i) {
 //  double fx = f(x);
 //  double dfx = df(x);
 //  if (std::abs(dfx) < 1e-14) break;
 //  double x1 = x - fx / dfx;
 //  if (std::abs(x1 - x) < tol) return x1;
 //  x = x1;
 // }
 // throwDomain(-1); // 収束失敗
 //};

 //  NPV
 inline double npv(const std::vector<double> &cf, double r) {
  double sum = 0.0;
  for (size_t t = 0; t < cf.size(); ++t)
   sum += cf[t] / std::pow(1.0 + r, t);
  return sum;
 }

 //  根のブラケット探索
 template <typename Func> std::pair<double, double> bracketRoot(Func f, auto &ctx, double start = -0.999, double end = 10.0, double step = 0.001, int maxSteps = 10000) {
  double prev = f(start);
  for (int i = 1; i <= maxSteps && start + i * step <= end; ++i) {
   double x = start + i * step;
   double curr = f(x);
   if (prev * curr <= 0.0) return {start + (i - 1) * step, x};
   prev = curr;
  }
  throw CalcError(CalcErrorType::NonConvergence, "IRR root not bracketed in search range", ctx.pos);
 }

 // --- Brent法 ---
 double brent(std::function<double(double)> f, double a, double b, FunctionContext &ctx, double tol = cnst_precision_inv, int maxIter = 100);

 // 全区間スキャン(短区間探索)
 std::pair<double, double> fullScanBracket(std::function<double(double)> f, double start, double end, double step, FunctionContext &ctx);
 double zetaEulerMaclaurin(double s);
 inline double betacf(double a, double b, double x) {
  const int MAXIT = 200;
  const double EPS = 3.0e-14;
  const double FPMIN = 1e-30;

  double qab = a + b;
  double qap = a + 1.0;
  double qam = a - 1.0;

  double c = 1.0;
  double d = 1.0 - qab * x / qap;
  if (std::fabs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  double h = d;

  for (int m = 1; m <= MAXIT; ++m) {
   int m2 = 2 * m;

   double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
   d = 1.0 + aa * d;
   if (std::fabs(d) < FPMIN) d = FPMIN;
   c = 1.0 + aa / c;
   if (std::fabs(c) < FPMIN) c = FPMIN;
   d = 1.0 / d;
   h *= d * c;

   aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
   d = 1.0 + aa * d;
   if (std::fabs(d) < FPMIN) d = FPMIN;
   c = 1.0 + aa / c;
   if (std::fabs(c) < FPMIN) c = FPMIN;
   d = 1.0 / d;
   double del = d * c;
   h *= del;

   if (std::fabs(del - 1.0) < EPS) break;
  }

  return h;
 }

 /* ============================
   統計系
   ============================ */
 std::vector<double> gatherReals(const std::vector<Value> &v, size_t pos);
 std::vector<double> collectReals(const std::vector<Value> &v, FunctionContext &ctx);
 std::vector<Complex> collectComplex(const std::vector<Value> &v, FunctionContext &ctx);
 double meanOf(const std::vector<double> &a);
 double stddevPopulation(const std::vector<double> &a, double mu);
 double quantileLinear(std::vector<double> a, double p, size_t pos);
 std::vector<double> rankAverageTies(const std::vector<double> &x, size_t pos);
 double pearsonCorr(const std::vector<double> &x, const std::vector<double> &y, size_t pos);
 double medianOfSorted(const std::vector<double> &a);
 double variancePopulation(const std::vector<double> &a, double mu);
 Value area_polygon_impl(const std::vector<Value> &v, FunctionContext &ctx);
 inline std::vector<double> collectNumericVector(const std::vector<Value> &v, const FunctionContext &ctx);
} // namespace mm::cal
