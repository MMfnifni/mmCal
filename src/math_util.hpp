#pragma once

#include "core.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace mm::cal {

 /* ============================
   数学補助系
   ============================ */

 static constexpr double deg2rad = PI / 180.0;
 static constexpr double rad2deg = 180.0 / PI;

 inline constexpr double toDeg(double rad) { return rad * (180.0 / PI); }
 inline constexpr double toRad(double deg) { return deg * (PI / 180.0); }
 inline constexpr double radToDeg(double r) { return r * 180.0 / PI; }

 inline double inf(int sign = +1) { return sign >= 0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity(); }

 // Complex(x, 0.0)
 static inline Complex cx(double x) { return {x, 0.0}; }

 // 「1/denom」を安全に返す関数
 Value safeInv(double denom, double sign);

 static inline Value signedInfBy(double sign) { return inf(sign >= 0 ? +1 : -1); }
 // denom が 0 に近いときに ±inf を返す（符号は signSource から取る）
 static inline Value invOrSignedInf(double denom, double signSource) { return (std::abs(denom) < cnst_precision_inv) ? signedInfBy(signSource) : (1.0 / denom); }

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
  if (std::floor(d) != d) throw CalcError(CalcErrorType::TypeError, errorMessage(CalcErrorType::TypeError), pos);
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
  if (n < 4) throw CalcError(CalcErrorType::DomainError, "DomainError: kurts: need at least 4 elements", ctx.pos);

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

 inline double corrPopulationReal(const std::vector<Value> &v, FunctionContext &ctx) {
  const int total = (int)v.size();
  if (total % 2 != 0) throw CalcError(CalcErrorType::SyntaxError, "corr: argument count must be even", ctx.pos);

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

  if (sxx == 0.0L || syy == 0.0L) throw CalcError(CalcErrorType::DomainError, "DomainError: corr: zero variance", ctx.pos);

  return static_cast<double>(sxy / std::sqrt(sxx * syy));
 }

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

 inline long double abs_cnst_precision_inv_ld() { return 1e-18L; }

 inline bool nearly_zero_ld(long double x, long double cnst_precision_inv = 1e-18L) { return std::fabsl(x) <= cnst_precision_inv; }

 inline long double pow1p_int_ld(long double r, int n) {
  // (1+r)^n を log1p / exp で安定に計算
  return std::expl((long double)n * std::log1pl(r));
 }

 inline long double kahan_sum_ld(const std::vector<long double> &xs) {
  long double sum = 0.0L;
  long double c = 0.0L;
  for (long double x : xs) {
   long double y = x - c;
   long double t = sum + y;
   c = (t - sum) - y;
   sum = t;
  }
  return sum;
 }

 inline long double fin_pmt_ld(long double rate, int nper, long double pv, int type) {
  if (nper <= 0) throw std::runtime_error("fin_pmt_ld: nper <= 0");

  if (nearly_zero_ld(rate)) {
   return pv / (long double)nper; // 正の支払額
  }

  const long double A = pow1p_int_ld(rate, nper);
  long double pmt = pv * rate * A / (A - 1.0L); // 期末払い
  if (type == 1) {
   pmt /= (1.0L + rate); // 期首払い補正
  }
  return pmt; // 正の支払額
 }

 //  NPV
 inline double npv(const std::vector<double> &cf, double r) {
  double sum = 0.0;
  for (size_t t = 0; t < cf.size(); ++t)
   sum += cf[t] / std::pow(1.0 + r, t);
  return sum;
 }

 //  根のブラケット探索
 template <typename Func> std::pair<double, double> bracketRoot(Func f, auto &ctx, double start = -0.999, double end = 10.0, double step = 0.001, int maxStcnst_precision_inv = 10000) {
  double prev = f(start);
  for (int i = 1; i <= maxStcnst_precision_inv && start + i * step <= end; ++i) {
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
  const double cnst_precision_inv = 3.0e-14;
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

   if (std::fabs(del - 1.0) < cnst_precision_inv) break;
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

 /* ============================
   ベクトル系
   ============================ */
 Value fromVector(const std::vector<double> &xs);
 Value fromIndexVector(const std::vector<size_t> &xs);
 Value fromMatrix(const std::vector<std::vector<double>> &A);
 std::string requireFunctionName(const Value &v, size_t pos);
 double defaultDiffStep(double x);
 double defaultDiffStep2(double x);
 std::shared_ptr<MultiValue> makeFactorList(long long n, size_t pos);
 double binomGeneral(double x, double y, size_t pos);
 double fallfactGeneral(double x, double n, size_t pos);
 double risefactGeneral(double x, double n, size_t pos);
 size_t nrows(const std::vector<std::vector<double>> &A);
 size_t ncols(const std::vector<std::vector<double>> &A);
 std::vector<std::vector<double>> identity(size_t n);
 bool isSymmetricMatrix(const std::vector<std::vector<double>> &A, double tol);
 std::vector<double> mat_vec(const std::vector<std::vector<double>> &A, const std::vector<double> &x);
 std::vector<std::vector<double>> mat_mul(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);
 std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &A);
 double norm(const std::vector<double> &x);
 double dot(const std::vector<double> &a, const std::vector<double> &b);
 double power_method(const std::vector<std::vector<double>> &A, std::vector<double> &eigenvec, int max_iter = 1000);
 void hessenberg(std::vector<std::vector<double>> &A);
 void qr_decompose(const std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &Q, std::vector<std::vector<double>> &R);
 std::vector<double> eigenvalues(std::vector<std::vector<double>> A, int max_iter = 1000);
 void lu_decompose(std::vector<std::vector<double>> A, std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U, std::vector<size_t> &P);
 std::vector<double> inverse_iteration(const std::vector<std::vector<double>> &A, double lambda);
 std::vector<double> singular_values(const std::vector<std::vector<double>> &A);
 double condition_number(const std::vector<std::vector<double>> &A);
 void bidiagonalize(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &U, std::vector<std::vector<double>> &V);
 void bidiagonal_qr(std::vector<std::vector<double>> &B, std::vector<std::vector<double>> &U, std::vector<std::vector<double>> &V);
 void svd_jacobi(const std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &U, std::vector<double> &S, std::vector<std::vector<double>> &V);
 bool isPowerOfTwo(size_t n);
 std::vector<std::complex<double>> toComplexVector1D(const Value &v, FunctionContext &ctx);
 Value fromComplexVector(const std::vector<std::complex<double>> &vec);
 std::vector<std::complex<double>> dft_impl(const std::vector<std::complex<double>> &x, bool inverse);
 void fft_rec(std::vector<std::complex<double>> &a);
 std::vector<std::complex<double>> fft_dispatch(const std::vector<std::complex<double>> &x, bool inverse, FunctionContext &ctx);
 Value fft2D(const Value &v, bool inverse, FunctionContext &ctx);

 inline void validate_matrix_shape(const std::vector<std::vector<double>> &A, const char *name) {
  if (A.empty()) return;

  const size_t cols = A[0].size();
  for (size_t i = 1; i < A.size(); ++i) {
   if (A[i].size() != cols) { throw std::runtime_error(std::string(name) + ": irregular matrix"); }
  }
 }

 inline void ensureFlatVectorValue(const Value &v, FunctionContext &ctx, const char *name) {
  if (!v.isMulti()) { throw CalcError(CalcErrorType::TypeError, std::string(name) + " requires vector argument", ctx.pos); }

  const auto &mv = v.asMultiRef(ctx.pos);
  for (const auto &e : mv.elems()) {
   if (e.isMulti()) { throw CalcError(CalcErrorType::DomainError, std::string("DomainError: ") + name + " requires 1-dimensional vector", ctx.pos); }
  }
 }

 inline std::vector<double> toVectorChecked(const Value &v, FunctionContext &ctx, const char *name) {
  ensureFlatVectorValue(v, ctx, name);
  const auto &mv = v.asMultiRef(ctx.pos);

  std::vector<double> out;
  out.reserve(mv.size());
  for (const auto &e : mv.elems()) {
   out.push_back(e.asScalar(ctx.pos));
  }
  return out;
 }

 inline Value fromVectorValue(const std::vector<double> &xs) {
  auto mv = std::make_shared<MultiValue>();
  mv->elems_.reserve(xs.size());
  for (double x : xs) {
   mv->elems_.emplace_back(x);
  }
  return Value(mv);
 }

 inline void ensureMatrixValue(const Value &v, FunctionContext &ctx, const char *name) {
  if (!v.isMulti()) { throw CalcError(CalcErrorType::TypeError, std::string(name) + " requires matrix argument", ctx.pos); }

  const auto &outer = v.asMultiRef(ctx.pos);
  for (const auto &row : outer.elems()) {
   if (!row.isMulti()) { throw CalcError(CalcErrorType::DomainError, std::string("DomainError: ") + name + " requires 2-dimensional matrix", ctx.pos); }
  }
 }

 inline std::vector<std::vector<double>> toMatrixChecked(const Value &v, FunctionContext &ctx, const char *name) {
  ensureMatrixValue(v, ctx, name);
  auto A = toMatrix(v, ctx);
  validate_matrix_shape(A, name);
  return A;
 }

 inline bool isGammaPole(double z) { return z <= 0.0 && std::floor(z) == z; }

 inline double gammaSign(double z) {
  if (z > 0.0) return 1.0;
  const double k = std::floor(-z);
  return (static_cast<long long>(k) % 2 == 0) ? -1.0 : 1.0;
 }

} // namespace mm::cal
