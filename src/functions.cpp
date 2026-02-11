#include "functions.hpp"
#include "constants.hpp"
#include "core.hpp"
#include "evaluate.hpp"
#include "math_util.hpp"

namespace mm::cal {

 /* ============================
   地獄の無限関数定義編
   ============================ */

 static const double deg2rad = PI / 180.0;
 static const double rad2deg = 180.0 / PI;
 static const double eps = asDouble(constants.at("EPS"));

 ///!!!!
 static void registerBasicMath(SystemConfig &cfg);   // 基本関数
 static void registerExpLog(SystemConfig &cfg);      // log系
 static void registerTrig(SystemConfig &cfg);        // 三角関数
 static void registerHyperbolic(SystemConfig &cfg);  // hyper三角関数
 static void registerGeoVec(SystemConfig &cfg);      // 図形, ベクトル
 static void registerSpecialTrig(SystemConfig &cfg); // hyper三角関数
 static void registerStatistics(SystemConfig &cfg);  // 統計系
 static void registerComplex(SystemConfig &cfg);     // 複素
 static void registerRandom(SystemConfig &cfg);      // 乱数
 static void registerEngineering(SystemConfig &cfg); // 機械設計
 static void registerOthers(SystemConfig &cfg);      // アミューズ

 // Complex(x, 0.0)
 static inline Complex cx(double x) { return {x, 0.0}; }

 // 「1/denom」を安全に返す関数
 auto safeInv = [=](double denom, double sign) -> Value { return (std::abs(denom) < eps) ? inf(sign >= 0 ? +1 : -1) : (1.0 / denom); };
 static inline Value signedInfBy(double sign) { return inf(sign >= 0 ? +1 : -1); }
 // denom が 0 に近いときに ±inf を返す（符号は signSource から取る）
 static inline Value invOrSignedInf(double denom, double signSource) { return (std::abs(denom) < eps) ? signedInfBy(signSource) : (1.0 / denom); }

 static Value fn_log(const std::vector<Value> &v, FunctionContext &ctx) {
  const double eps = asDouble(constants.at("EPS"));
  auto warnPrincipal = [&] { calcWarn(ctx.cfg, ctx.pos, "log: complex principal value only; other branches may exist"); };
  auto safeLogBase = [&](Complex base) -> Complex {
   const auto lb = std::log(base);
   if (std::abs(lb) < eps) throwDomain(ctx.pos); // base == 1 (or too close)
   return lb;
  };

  auto logBase = [&](Complex base, Complex x) -> Complex { return std::log(x) / safeLogBase(base); };
  const bool complexMode = isComplex(v[0]) || (v.size() == 2 && isComplex(v[1]));
  // ---- complex mode (explicit) ----
  if (complexMode) {
   if (v.size() == 1) return std::log(asComplex(v[0]));
   return logBase(asComplex(v[0]), asComplex(v[1]));
  }

  // ---- real mode ----
  if (v.size() == 1) {
   const double x = asReal(v[0], ctx.pos);
   if (x == 0.0) throwDomain(ctx.pos); // log(0) は -inf だが電卓としてはエラーが自然
   if (x > 0.0) return std::log(x);
   // x < 0 -> complex principal value
   calcWarn(ctx.cfg, ctx.pos, "log(x<0): complex principal value only; branch cut on negative real axis");
   return std::log(cx(x));
  }
  // ---- base + x (both real) ----
  const double base = asReal(v[0], ctx.pos);
  const double x = asReal(v[1], ctx.pos);
  if (base == 0.0) throw CalcError(CalcErrorType::DomainError, "log: base must not be zero", ctx.pos);
  if (x == 0.0) throw CalcError(CalcErrorType::DomainError, "log: log(0) is undefined", ctx.pos);
  // real-safe region
  if (base > 0.0 && base != 1.0 && x > 0.0) return std::log(x) / std::log(base);
  // otherwise promote to complex (principal value)
  warnPrincipal();
  return logBase(Complex(base, 0.0), cx(x));
 }
 // sum, prod用
 static inline Value complexToValue(Complex z) { return (std::abs(z.imag()) < cnst_precision_inv) ? Value(z.real()) : Value(z); }
 // 実数チェックの御老体
 static inline double requireRealNoComplex(const Value &x, size_t pos, const char *msg) {
  if (isComplex(x)) {
   const auto z = asComplex(x);
   if (z.imag() != 0.0) throw CalcError(CalcErrorType::DomainError, msg, pos);
   return z.real();
  }
  return asReal(x, pos);
 }
 static inline bool isNegativeInteger(double x, double eps) {
  if (!(x < 0.0)) return false;
  double r = std::round(x);
  return std::abs(x - r) < eps;
 }
 static double varianceRealWelford(const std::vector<Value> &v, FunctionContext &ctx, int ddof) {
  const int n = (int)v.size();
  if (n <= ddof) throw CalcError(CalcErrorType::DomainError, "var: too few elements", ctx.pos);
  if (n == 1) return 0.0;

  long double mean = 0.0L;
  long double m2 = 0.0L;
  int k = 0;

  for (auto &x : v) {
   const long double a = (long double)asReal(x, ctx.pos);
   ++k;
   const long double d = a - mean;
   mean += d / k;
   m2 += d * (a - mean);
  }

  return (double)(m2 / (long double)(n - ddof));
 }
 // Welfordで mean + M2 を返す変態設計
 struct MeanM2 {
   long double mean = 0.0L;
   long double m2 = 0.0L;
   int n = 0;
 };

 static MeanM2 welfordMeanM2Real(const std::vector<Value> &v, FunctionContext &ctx) {
  MeanM2 s;
  for (const auto &x : v) {
   const long double a = (long double)asReal(x, ctx.pos);
   ++s.n;
   const long double d = a - s.mean;
   s.mean += d / (long double)s.n;
   s.m2 += d * (a - s.mean);
  }
  return s;
 }

 static double medianInplace(std::vector<double> &a) {
  const size_t n = a.size();
  if (n == 0) return 0.0;
  const size_t mid = n / 2;
  std::nth_element(a.begin(), a.begin() + mid, a.end());
  const double hi = a[mid];
  if (n & 1) return hi;
  const double lo = *std::max_element(a.begin(), a.begin() + mid);
  return 0.5 * (lo + hi);
 }

 static double skewnessReal(const std::vector<Value> &v, FunctionContext &ctx) {
  const size_t n = v.size();
  if (n == 0) throw CalcError(CalcErrorType::DomainError, "skew: no elements", ctx.pos);

  long double mean = 0.0L;
  long double m2 = 0.0L;
  long double m3 = 0.0L;
  size_t k = 0;

  for (auto &xv : v) {
   const long double x = (long double)asReal(xv, ctx.pos);
   const long double k1 = (long double)k;
   ++k;

   const long double delta = x - mean;
   const long double delta_n = delta / (long double)k;
   const long double delta_n2 = delta_n * delta_n;
   const long double term1 = delta * delta_n * k1;

   mean += delta_n;
   m3 += term1 * delta_n * (k1 - 1.0L) - 3.0L * delta_n * m2;
   m2 += term1;
  }

  if (m2 == 0.0L) throwDomain(ctx.pos);

  // population skewness: (m3/n) / ( (m2/n)^(3/2) )  == m3 / m2^(3/2)
  const long double denom = std::pow(m2, 1.5L);
  return (double)(m3 / denom);
 }

 static double kurtosisExcessReal(const std::vector<Value> &v, FunctionContext &ctx) {
  const size_t n = v.size();
  if (n == 0) throw CalcError(CalcErrorType::DomainError, "kurt: no elements", ctx.pos);
  // 1st pass: mean
  long double mean = 0.0L;
  for (const auto &xv : v)
   mean += (long double)asReal(xv, ctx.pos);
  mean /= (long double)n;
  // 2nd pass: central moments
  long double m2 = 0.0L;
  long double m4 = 0.0L;
  for (const auto &xv : v) {
   const long double x = (long double)asReal(xv, ctx.pos);
   const long double d = x - mean;
   const long double d2 = d * d;
   m2 += d2;
   m4 += d2 * d2;
  }
  if (m2 == 0.0L) throwDomain(ctx.pos);
  // population excess kurtosis:
  // (m4/n) / ((m2/n)^2) - 3  ==  n*m4/(m2*m2) - 3
  const long double nld = (long double)n;
  return (double)(nld * m4 / (m2 * m2) - 3.0L);
 }

 static double kurtosisExcessSample(const std::vector<Value> &v, FunctionContext &ctx) {
  const size_t n = v.size();
  if (n < 4) throw CalcError(CalcErrorType::DomainError, "kurts: need at least 4 elements", ctx.pos);

  long double mean = 0.0L;
  long double m2 = 0.0L, m3 = 0.0L, m4 = 0.0L;
  size_t k = 0;

  for (auto &xv : v) {
   const long double x = (long double)asReal(xv, ctx.pos);
   const long double k1 = (long double)k;
   ++k;

   const long double delta = x - mean;
   const long double delta_n = delta / (long double)k;
   const long double delta_n2 = delta_n * delta_n;
   const long double term1 = delta * delta_n * k1;

   mean += delta_n;
   m4 += term1 * delta_n2 * (k1 * k1 - 3.0L * k1 + 3.0L) + 6.0L * delta_n2 * m2 - 4.0L * delta_n * m3;
   m3 += term1 * delta_n * (k1 - 1.0L) - 3.0L * delta_n * m2;
   m2 += term1;
  }

  if (m2 == 0.0L) throwDomain(ctx.pos);

  // g2: population excess kurtosis
  const long double g2 = m4 / (m2 * m2) - 3.0L;

  // G2: unbiased / bias-corrected sample excess kurtosis
  const long double nn = (long double)n;
  const long double G2 = ((nn - 1.0L) / ((nn - 2.0L) * (nn - 3.0L))) * ((nn + 1.0L) * g2 + 6.0L);

  return (double)G2;
 }

 static double quantileLinearSorted(const std::vector<double> &a, double p01, size_t pos) {
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
 static double covariancePopulationReal(const std::vector<Value> &v, FunctionContext &ctx) {
  const int total = (int)v.size();
  if (total % 2 != 0) throw CalcError(CalcErrorType::DomainError, "cov: invalid argument count", ctx.pos);

  const int n = total / 2;
  if (n < 2) throw CalcError(CalcErrorType::DomainError, "cov: too few samples", ctx.pos);

  long double mx = 0.0L, my = 0.0L;
  long double c = 0.0L;
  int k = 0;

  for (int i = 0; i < n; ++i) {
   const long double x = (long double)asReal(v[i], ctx.pos);
   const long double y = (long double)asReal(v[i + n], ctx.pos);

   ++k;
   const long double dx = x - mx;
   const long double dy = y - my;

   mx += dx / k;
   my += dy / k;

   c += dx * (y - my); // Welford型
  }

  return (double)(c / (long double)n); // population
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
   const long double x = (long double)asReal(v[i], ctx.pos);
   const long double y = (long double)asReal(v[i + n], ctx.pos);

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

  return (double)(sxy / std::sqrt(sxx * syy));
 }

 // 実数に落とせるなら落とす
 static inline Value realIfPossible(Complex z) { return (std::abs(z.imag()) < cnst_precision_inv) ? Value(z.real()) : Value(z); }

 static inline bool isZero(Complex z) { return z.real() == 0.0 && z.imag() == 0.0; }

 static inline Complex sin_over_t_series(Complex t) {
  // sin(t)/t = 1 - t^2/6 + t^4/120 - t^6/5040 + t^8/362880
  const Complex t2 = t * t;
  const Complex t4 = t2 * t2;
  const Complex t6 = t4 * t2;
  const Complex t8 = t4 * t4;
  return Complex(1, 0) - t2 / 6.0 + t4 / 120.0 - t6 / 5040.0 + t8 / 362880.0;
 }

 static inline Complex tan_over_t_series(Complex t) {
  // tan(t)/t = 1 + t^2/3 + 2 t^4/15 + 17 t^6/315 + 62 t^8/2835
  const Complex t2 = t * t;
  const Complex t4 = t2 * t2;
  const Complex t6 = t4 * t2;
  const Complex t8 = t4 * t4;
  return Complex(1, 0) + t2 / 3.0 + (Complex(2, 0) * t4) / 15.0 + (Complex(17, 0) * t6) / 315.0 + (Complex(62, 0) * t8) / 2835.0;
 }

 static inline double absC(Complex z) { return std::abs(z); } // 読みやすさ用

 static FuncImpl makeSincLike(double scale) {
  return [=](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   const Complex x = asComplex(v[0]);
   if (isZero(x)) return 1.0;

   const Complex t = x * scale;

   // sinc(x) = sin(scale*x)/x = sin(t)/x
   // small t: (sin(t)/t) * (t/x) = (sin(t)/t) * scale
   const Complex r = (absC(t) < 1e-8) ? (sin_over_t_series(t) * scale) : (std::sin(t) / x);

   return realIfPossible(r);
  };
 }
 static FuncImpl makeCoscLike(double scale) {
  return [=](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   const Complex x = asComplex(v[0]);
   if (isZero(x)) return 0.0;

   const Complex t = x * scale;
   const Complex h = t * 0.5;

   // cosc(x) = (1 - cos(t))/x = 2*sin(t/2)^2 / x
   const Complex s = std::sin(h);
   const Complex r = (Complex(2, 0) * s * s) / x;

   return realIfPossible(r);
  };
 }
 static FuncImpl makeTancLike(double scale) {
  return [=](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   const Complex x = asComplex(v[0]);
   if (isZero(x)) return 1.0;

   const Complex t = x * scale;

   // small t: tan(t)/t * scale
   if (absC(t) < 1e-8) { return realIfPossible(tan_over_t_series(t) * scale); }

   // 元仕様: cos が小さいなら inf（複素でも適用）
   const Complex c = std::cos(t);
   if (std::abs(c) < cnst_precision_inv) {
    const Complex s = std::sin(t);
    return inf(std::real(s) >= 0 ? +1 : -1);
   }

   const Complex r = std::tan(t) / x;
   return realIfPossible(r);
  };
 }

 static FuncImpl makeDivXCmplxReal(std::function<Complex(Complex)> fc, std::function<double(double)> fr, double x0Value) {
  return [=](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   if (isComplex(v[0])) {
    Complex x = asComplex(v[0]);
    if (x == Complex(0, 0)) return x0Value;
    return fc(x) / x;
   }
   double x = asReal(v[0], ctx.pos);
   if (x == 0.0) return x0Value;
   return fr(x) / x;
  };
 }

 static FuncImpl makeExpc() {
  return [](const std::vector<Value> &v, FunctionContext &ctx) -> Value {
   if (isComplex(v[0])) {
    Complex x = asComplex(v[0]);
    if (x == Complex(0, 0)) return 1.0;

    // exp(x)-1 をなるべく桁落ちしないようにする
    // 複素には std::expm1 が無いので小さい時はテイラー
    double ax = std::abs(x);
    if (ax < 1e-8) {
     // (exp(x)-1)/x = 1 + x/2 + x^2/6 + x^3/24 + ...
     Complex x2 = x * x;
     Complex x3 = x2 * x;
     Complex x4 = x2 * x2;
     return Complex(1, 0) + x / 2.0 + x2 / 6.0 + x3 / 24.0 + x4 / 120.0;
    }

    return (std::exp(x) - Complex(1, 0)) / x;
   }

   double x = asReal(v[0], ctx.pos);
   if (x == 0.0) return 1.0;
   return std::expm1(x) / x;
  };
 }

 static inline long long absLL(long long x, size_t pos) {
  if (x == LLONG_MIN) throwOverflow(pos);
  return std::llabs(x);
 }

 void registerBasicMath(SystemConfig &cfg) { // 基本関数
  cfg.functions["abs"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           if (isComplex(v[0])) return std::abs(asComplex(v[0]));
                           return std::fabs(asReal(v[0], ctx.pos));
                          }};
  cfg.functions["sign"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = asReal(v[0], ctx.pos);
                            return (x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0;
                           }};
  cfg.functions["sqrt"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (isComplex(v[0])) return std::sqrt(asComplex(v[0]));
                            double x = asReal(v[0], ctx.pos);
                            if (x < 0) return Complex(0, std::sqrt(-x));
                            return std::sqrt(x);
                           }};
  cfg.functions["cbrt"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::cbrt(asReal(v[0], ctx.pos)); }};
  cfg.functions["floor"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::floor(asReal(v[0], ctx.pos)); }};
  cfg.functions["ceil"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::ceil(asReal(v[0], ctx.pos)); }};
  cfg.functions["trunc"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::trunc(asReal(v[0], ctx.pos)); }};
  cfg.functions["pow"] = {2, 2, [](auto &v, auto &ctx) -> Value { return std::pow(asComplex(v[0]), asComplex(v[1])); }};

  // ---- coonvert angle ----
  cfg.functions["DtoR"] = {1, 1, [](auto &v, auto &ctx) -> Value { return deg2rad_v(v[0], ctx.pos); }};
  cfg.functions["DtoG"] = {1, 1, [](auto &v, auto &ctx) -> Value { return deg2grad_v(v[0], ctx.pos); }};
  cfg.functions["RtoD"] = {1, 1, [](auto &v, auto &ctx) -> Value { return rad2deg_v(v[0], ctx.pos); }};
  cfg.functions["RtoG"] = {1, 1, [](auto &v, auto &ctx) -> Value { return rad2grad_v(v[0], ctx.pos); }};
  cfg.functions["GtoD"] = {1, 1, [](auto &v, auto &ctx) -> Value { return grad2deg_v(v[0], ctx.pos); }};
  cfg.functions["GtoR"] = {1, 1, [](auto &v, auto &ctx) -> Value { return grad2rad_v(v[0], ctx.pos); }};
 }

 void registerExpLog(SystemConfig &cfg) { // log, exp系
  cfg.functions["exp"] = {1, 1, [](auto &v, auto &) -> Value { return std::exp(asComplex(v[0])); }};
  cfg.functions["log10"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             double x = asReal(v[0], ctx.pos);
                             if (x <= 0) throwDomain(ctx.pos);
                             return std::log10(x);
                            }};
  cfg.functions["log2"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = asReal(v[0], ctx.pos);
                            if (x <= 0) throwDomain(ctx.pos);
                            return std::log2(x);
                           }};
  cfg.functions["log"] = {1, 2, fn_log};
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
 }

 void registerTrig(SystemConfig &cfg) {
  cfg.functions["sin"] = {1, 1, [=](auto &v, auto &) -> Value { return std::sin(asComplex(v[0]) * deg2rad); }};
  cfg.functions["cos"] = {1, 1, [=](auto &v, auto &) -> Value { return std::cos(asComplex(v[0]) * deg2rad); }};
  // cfg.functions["tan"] = {1, 1, [=](auto &v, auto &) -> Value { return std::tan(asComplex(v[0]) * deg2rad); }};
  cfg.functions["tan"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           if (isComplex(v[0])) return std::tan(asComplex(v[0]) * deg2rad);
                           const double r = asReal(v[0], ctx.pos) * deg2rad;
                           const double c = std::cos(r);
                           return (std::abs(c) < eps) ? signedInfBy(std::sin(r)) : std::tan(r);
                          }};
  cfg.functions["cot"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           if (isComplex(v[0])) {
                            const auto z = asComplex(v[0]) * deg2rad;
                            return std::cos(z) / std::sin(z);
                           }
                           const double r = asReal(v[0], ctx.pos) * deg2rad;
                           const double s = std::sin(r);
                           return (std::abs(s) < eps) ? signedInfBy(std::cos(r)) : (std::cos(r) / s);
                          }};
  cfg.functions["sec"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           if (isComplex(v[0])) return 1.0 / std::cos(asComplex(v[0]) * deg2rad);
                           const double r = asReal(v[0], ctx.pos) * deg2rad;
                           const double c = std::cos(r);
                           return (std::abs(c) < eps) ? signedInfBy(std::sin(r)) : (1.0 / c);
                          }};
  cfg.functions["csc"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           if (isComplex(v[0])) return 1.0 / std::sin(asComplex(v[0]) * deg2rad);

                           const double r = asReal(v[0], ctx.pos) * deg2rad;
                           const double s = std::sin(r);
                           return (std::abs(s) < eps) ? signedInfBy(std::cos(r)) : (1.0 / s);
                          }};
  cfg.functions["asin"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (isComplex(v[0])) return std::asin(asComplex(v[0])) * rad2deg;
                            const double x = asReal(v[0], ctx.pos);
                            if (std::abs(x) <= 1.0) return std::asin(x) * rad2deg;
                            calcWarn(ctx.cfg, ctx.pos, "asin(|x|>1): complex principal value only");
                            return std::asin(cx(x)) * rad2deg;
                           }};
  cfg.functions["acos"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (isComplex(v[0])) return std::acos(asComplex(v[0])) * rad2deg;
                            const double x = asReal(v[0], ctx.pos);
                            if (std::abs(x) <= 1.0) return std::acos(x) * rad2deg;
                            calcWarn(ctx.cfg, ctx.pos, "acos(|x|>1): complex principal value only");
                            return std::acos(cx(x)) * rad2deg;
                           }};
  cfg.functions["atan"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (isComplex(v[0])) return std::atan(asComplex(v[0])) * rad2deg;
                            return std::atan(asReal(v[0], ctx.pos)) * rad2deg;
                           }};
  cfg.functions["atan2"] = {2, 2, [=](auto &v, auto &ctx) -> Value {
                             double y = asReal(v[0], ctx.pos);
                             double x = asReal(v[1], ctx.pos);
                             if (x == 0.0 && y == 0.0) throwDomain(ctx.pos);
                             return std::atan2(y, x) * rad2deg;
                            }};
 }

 void registerHyperbolic(SystemConfig &cfg) { // hyper三角関数
  cfg.functions["sinh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (isComplex(v[0])) return std::sinh(asComplex(v[0]));
                            return std::sinh(asReal(v[0], ctx.pos));
                           }};
  cfg.functions["cosh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (isComplex(v[0])) return std::cosh(asComplex(v[0]));
                            return std::cosh(asReal(v[0], ctx.pos));
                           }};
  cfg.functions["tanh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (isComplex(v[0])) return std::tanh(asComplex(v[0]));
                            return std::tanh(asReal(v[0], ctx.pos));
                           }};

  // ---- inverse hyperbolic ----
  cfg.functions["asinh"] = {1, 1, [](auto &v, auto &ctx) -> Value { return isComplex(v[0]) ? Value(std::asinh(asComplex(v[0]))) : Value(std::asinh(asReal(v[0], ctx.pos))); }};

  cfg.functions["acosh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (isComplex(v[0])) return std::acosh(asComplex(v[0]));

                             const double x = asReal(v[0], ctx.pos);
                             if (x >= 1.0) return std::acosh(x);

                             calcWarn(ctx.cfg, ctx.pos, "acosh(x<1): complex principal value only");
                             return std::acosh(cx(x));
                            }};

  cfg.functions["atanh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (isComplex(v[0])) return std::atanh(asComplex(v[0]));

                             const double x = asReal(v[0], ctx.pos);
                             if (std::abs(x) < 1.0) return std::atanh(x);

                             calcWarn(ctx.cfg, ctx.pos, "atanh(|x|>=1): complex principal value only");
                             return std::atanh(cx(x));
                            }};
 }

 void registerStatistics(SystemConfig &cfg) {
  cfg.functions["hypot"] = {2, 2, [](auto &v, auto &ctx) -> Value { return std::hypot(asReal(v[0], ctx.pos), asReal(v[1], ctx.pos)); }};
  cfg.functions["fact"] = {1, 1, [](auto &v, auto &ctx) -> Value { return (double)factLD(requireInt(v[0], ctx.pos), ctx.pos); }};
  cfg.functions["gcd"] = {2, 2, [](auto &v, auto &ctx) -> Value { return (double)gcdLL(requireInt(v[0], ctx.pos), requireInt(v[1], ctx.pos)); }};
  cfg.functions["comb"] = {2, 2, [](auto &v, auto &ctx) -> Value { return (double)combLL(requireInt(v[0], ctx.pos), requireInt(v[1], ctx.pos), ctx.pos); }};
  cfg.functions["perm"] = {2, 2, [](auto &v, auto &ctx) -> Value { return (double)permLL(requireInt(v[0], ctx.pos), requireInt(v[1], ctx.pos), ctx.pos); }};
  cfg.functions["lcm"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           const long long a = requireInt(v[0], ctx.pos);
                           const long long b = requireInt(v[1], ctx.pos);
                           if (a == 0 || b == 0) return 0.0;
                           const long long g = gcdLL(a, b);
                           const long long x = a / g;
                           const auto ax = absLL(x, ctx.pos);
                           const auto ab = absLL(b, ctx.pos);
                           return (double)checkedMul(ax, ab, ctx.pos);
                          }};
  cfg.functions["sum"] = {0, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.empty()) throwDomain(ctx.pos);
                           Complex acc{0, 0};
                           for (const auto &x : v)
                            acc += asComplex(x);
                           return complexToValue(acc);
                          }};
  cfg.functions["prod"] = {0, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.empty()) throwDomain(ctx.pos);
                            Complex acc{1, 0};
                            for (const auto &x : v)
                             acc *= asComplex(x);
                            return complexToValue(acc);
                           }};
  cfg.functions["mean"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.empty()) throwDomain(ctx.pos);
                            Complex acc{0, 0};
                            for (const auto &x : v)
                             acc += asComplex(x);
                            acc /= (double)v.size();
                            return complexToValue(acc);
                           }};
  cfg.functions["mod"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           const double x = requireRealNoComplex(v[0], ctx.pos, "mod: complex argument");
                           const double y = requireRealNoComplex(v[1], ctx.pos, "mod: complex argument");
                           if (y == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), ctx.pos);
                           return x - y * std::floor(x / y);
                          }};
  cfg.functions["geomean"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                               if (v.empty()) throwDomain(ctx.pos);
                               long double sumLog = 0.0L;
                               for (const auto &x : v) {
                                if (isComplex(x)) throw CalcError(CalcErrorType::DomainError, "geomean: complex", ctx.pos);
                                const double a = asReal(x, ctx.pos);
                                if (a < 0.0) throw CalcError(CalcErrorType::DomainError, "geomean: negative", ctx.pos);
                                if (a == 0.0) return 0.0;
                                sumLog += std::log((long double)a);
                               }
                               const long double r = std::exp(sumLog / (long double)v.size());
                               if (!std::isfinite((double)r)) throwOverflow(ctx.pos);
                               return (double)r;
                              }};
  cfg.functions["harmmean"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                                if (v.empty()) throwDomain(ctx.pos);
                                long double acc = 0.0L;
                                for (const auto &x : v) {
                                 if (isComplex(x)) throw CalcError(CalcErrorType::DomainError, "harmmean: complex", ctx.pos);
                                 const double a = asReal(x, ctx.pos);
                                 if (a == 0.0) throw CalcError(CalcErrorType::DivisionByZero, "harmmean: zero element", ctx.pos);
                                 acc += 1.0L / (long double)a;
                                }
                                const long double r = (long double)v.size() / acc;
                                if (!std::isfinite((double)r)) throwOverflow(ctx.pos);
                                return (double)r;
                               }};

  cfg.functions["quantile"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                const double p = asReal(v[0], ctx.pos);
                                if (!(p >= 0.0 && p <= 1.0)) throw CalcError(CalcErrorType::DomainError, "quantile: p out of range", ctx.pos);
                                if (v.size() < 2) throw CalcError(CalcErrorType::DomainError, "quantile: no samples", ctx.pos);
                                std::vector<double> a;
                                a.reserve(v.size() - 1);
                                for (size_t i = 1; i < v.size(); ++i)
                                 a.push_back(asReal(v[i], ctx.pos));
                                const size_t n = a.size();
                                const double idx = p * (double)(n - 1);
                                size_t i0 = (size_t)idx;
                                size_t i1 = i0 + 1;
                                if (i0 >= n - 1) return *std::max_element(a.begin(), a.end()); // p==1 対策
                                std::nth_element(a.begin(), a.begin() + i0, a.end());          // i0 番目
                                const double x0 = a[i0];
                                std::nth_element(a.begin() + i0 + 1, a.begin() + i1, a.end()); // i1 番目（後ろ側だけでOK）
                                const double x1 = a[i1];
                                const double t = idx - (double)i0;
                                return std::lerp(x0, x1, t);
                               }};

  cfg.functions["min"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.empty()) throwDomain(ctx.pos);
                           double m = asReal(v[0], ctx.pos);
                           for (size_t i = 1; i < v.size(); ++i)
                            m = std::min(m, asReal(v[i], ctx.pos));
                           return m;
                          }};
  cfg.functions["max"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.empty()) throwDomain(ctx.pos);
                           double m = asReal(v[0], ctx.pos);
                           for (size_t i = 1; i < v.size(); ++i)
                            m = std::max(m, asReal(v[i], ctx.pos));
                           return m;
                          }};
  cfg.functions["clamp"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                             const double x = asReal(v[0], ctx.pos);
                             const double lo = asReal(v[1], ctx.pos);
                             const double hi = asReal(v[2], ctx.pos);
                             if (lo > hi) throw CalcError(CalcErrorType::DomainError, "clamp: lo > hi", ctx.pos);
                             return std::clamp(x, lo, hi);
                            }};
  cfg.functions["fract"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             double x = asReal(v[0], ctx.pos);
                             return x - std::floor(x);
                            }};
  cfg.functions["gamma"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             const double x = asReal(v[0], ctx.pos);
                             // 負の整数は極
                             if (x < 0 && std::floor(x) == x) throw CalcError(CalcErrorType::DomainError, "gamma: pole at negative integer", ctx.pos);
                             const double r = std::tgamma(x);
                             if (std::isnan(r)) throwDomain(ctx.pos);
                             if (!std::isfinite(r)) throwOverflow(ctx.pos);
                             return r;
                            }};
  cfg.functions["lgamma"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              const double x = asReal(v[0], ctx.pos);
                              const double r = std::lgamma(x);
                              if (std::isnan(r)) throwDomain(ctx.pos);
                              if (!std::isfinite(r)) throwOverflow(ctx.pos);
                              return r;
                             }};

  cfg.functions["mode"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                            auto a = collectReals(v, ctx);
                            if (a.empty()) throwDomain(ctx.pos);
                            std::sort(a.begin(), a.end());
                            auto eq = [&](double x, double y) { return std::abs(x - y) <= cnst_precision_inv; };
                            double best = a[0], cur = a[0];
                            int bestCnt = 1, curCnt = 1;
                            for (size_t i = 1; i < a.size(); ++i) {
                             if (eq(a[i], cur)) {
                              ++curCnt;
                              continue;
                             }
                             if (curCnt > bestCnt) {
                              bestCnt = curCnt;
                              best = cur;
                             }
                             cur = a[i];
                             curCnt = 1;
                            }
                            if (curCnt > bestCnt) best = cur; // 最後の run
                            return best;                      // 同率なら最小値が残る（仕様として自然）
                           }};

  cfg.functions["erf"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::erf(asReal(v[0], ctx.pos)); }};
  cfg.functions["erfc"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::erfc(asReal(v[0], ctx.pos)); }};

  cfg.functions["round"] = {1, 2, [](auto &v, auto &ctx) -> Value {
                             const double x = asReal(v[0], ctx.pos);
                             if (v.size() == 1) return std::round(x);
                             const int n = (int)requireInt(v[1], ctx.pos);
                             if (n < -15 || n > 15) throw CalcError(CalcErrorType::OutOfRange, "round: n out of range (-15..15)", ctx.pos);
                             if (n == 0) return std::round(x);
                             // 冪乗テーブル
                             static constexpr double p10[] = {1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15};
                             const double s = p10[n + 15];
                             return std::round(x * s) / s;
                            }};
  cfg.functions["var"] = {1, -1, [](auto &v, auto &ctx) -> Value { return varianceRealWelford(v, ctx, 0); }};
  cfg.functions["vars"] = {1, -1, [](auto &v, auto &ctx) -> Value { return varianceRealWelford(v, ctx, 1); }};
  cfg.functions["stddev"] = {1, -1, [](auto &v, auto &ctx) -> Value { return std::sqrt(varianceRealWelford(v, ctx, 0)); }};
  cfg.functions["stddevs"] = {1, -1, [](auto &v, auto &ctx) -> Value { return std::sqrt(varianceRealWelford(v, ctx, 1)); }};
  cfg.functions["median"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                              if (v.empty()) throw CalcError(CalcErrorType::DomainError, "median: no elements", ctx.pos);
                              auto a = collectReals(v, ctx);
                              return medianInplace(a);
                             }};
  cfg.functions["mad"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                           if (v.empty()) throw CalcError(CalcErrorType::DomainError, "mad: no elements", ctx.pos);
                           auto a = collectReals(v, ctx);
                           const double med = medianInplace(a);
                           for (double &x : a)
                            x = std::abs(x - med);
                           return medianInplace(a);
                          }};
  cfg.functions["skew"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value { return skewnessReal(v, ctx); }};
  cfg.functions["kurtp"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value { return kurtosisExcessReal(v, ctx); }};
  cfg.functions["kurts"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value { return kurtosisExcessSample(v, ctx); }};
  cfg.functions["percentile"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                  const double p = asReal(v[0], ctx.pos);
                                  if (!(p >= 0.0 && p <= 100.0)) throw CalcError(CalcErrorType::DomainError, "percentile: p out of range", ctx.pos);
                                  if (v.size() < 2) throw CalcError(CalcErrorType::DomainError, "percentile: no samples", ctx.pos);
                                  std::vector<double> a;
                                  a.reserve(v.size() - 1);
                                  for (size_t i = 1; i < v.size(); ++i)
                                   a.push_back(asReal(v[i], ctx.pos));
                                  return quantileLinearSorted(std::move(a), p / 100.0, ctx.pos);
                                 }};
  cfg.functions["cov"] = {2, -1, [](auto &v, auto &ctx) -> Value { return covariancePopulationReal(v, ctx); }};
  cfg.functions["corr"] = {2, -1, [](auto &v, auto &ctx) -> Value { return corrPopulationReal(v, ctx); }};
  cfg.functions["ave"] = cfg.functions["mean"];
  cfg.functions["rms"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                           long double ss = 0.0L;
                           size_t n = 0;
                           for (auto &xv : v) {
                            double x = asReal(xv, ctx.pos);
                            ss += (long double)x * (long double)x;
                            ++n;
                           }
                           if (n == 0) throwDomain(ctx.pos);
                           return std::sqrt((double)(ss / (long double)n));
                          }};
  cfg.functions["cv"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                          if (v.size() < 2) throwDomain(ctx.pos);
                          const auto s = welfordMeanM2Real(v, ctx);
                          const double mu = (double)s.mean;
                          if (!std::isfinite(mu) || mu == 0.0) throwDomain(ctx.pos);

                          // population variance
                          const long double var = s.m2 / (long double)s.n;
                          if (!(var >= 0.0L)) throw CalcError(CalcErrorType::DomainError, "cv: invalid variance", ctx.pos);

                          const double sd = std::sqrt((double)var);
                          if (!std::isfinite(sd)) throwOverflow(ctx.pos);

                          return sd / mu;
                         }};

  cfg.functions["stderr"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                              if (v.size() < 2) throwDomain(ctx.pos);
                              const auto s = welfordMeanM2Real(v, ctx);
                              // population variance
                              const long double var = s.m2 / (long double)s.n;
                              if (!(var >= 0.0L)) throw CalcError(CalcErrorType::DomainError, "stderr: invalid variance", ctx.pos);
                              const double sd = std::sqrt((double)var);
                              const double se = sd / std::sqrt((double)s.n);

                              if (!std::isfinite(se)) throwOverflow(ctx.pos);
                              return se;
                             }};

  cfg.functions["zscore"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                              double x = asReal(v[0], ctx.pos);
                              double mu = asReal(v[1], ctx.pos);
                              double sigma = asReal(v[2], ctx.pos);
                              if (sigma == 0.0) throwDomain(ctx.pos);
                              return (x - mu) / sigma;
                             }};
  cfg.functions["iqr"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                           auto a = gatherReals(v, ctx.pos);
                           if (a.size() < 2) throwDomain(ctx.pos);
                           std::sort(a.begin(), a.end());
                           return quantileLinearSorted(a, 0.75, ctx.pos) - quantileLinearSorted(a, 0.25, ctx.pos); // quantileLinearでソート回数を減らすやつ
                          }};
  cfg.functions["trimmean"] = {2, INT_MAX, [](auto &v, auto &ctx) -> Value {
                                const double p = asReal(v[0], ctx.pos);
                                if (!(p >= 0.0 && p < 0.5)) throwDomain(ctx.pos);
                                // v[1..] を直接 gather できるならそれが理想
                                std::vector<double> a;
                                a.reserve(v.size() - 1);
                                for (size_t i = 1; i < v.size(); ++i)
                                 a.push_back(asReal(v[i], ctx.pos));

                                if (a.empty()) throwDomain(ctx.pos);

                                std::sort(a.begin(), a.end());
                                const size_t n = a.size();
                                size_t k = (size_t)std::floor(p * (double)n);
                                if (2 * k >= n) k = (n - 1) / 2; // pが0.4999...みたいな時にfloor(p*n)がずれるのを防止

                                const size_t cnt = n - 2 * k;
                                if (cnt == 0) throwDomain(ctx.pos);

                                long double s = 0.0L;
                                for (size_t i = k; i < n - k; ++i)
                                 s += (long double)a[i];

                                return (double)(s / (long double)cnt);
                               }};
  cfg.functions["winsor"] = {2, INT_MAX, [](auto &v, auto &ctx) -> Value {
                              const double p = asReal(v[0], ctx.pos);
                              if (!(p >= 0.0 && p < 0.5)) throwDomain(ctx.pos);

                              std::vector<double> a;
                              a.reserve(v.size() - 1);
                              for (size_t i = 1; i < v.size(); ++i)
                               a.push_back(asReal(v[i], ctx.pos));

                              if (a.empty()) throwDomain(ctx.pos);

                              std::sort(a.begin(), a.end());
                              const size_t n = a.size();
                              size_t k = (size_t)std::floor(p * (double)n);
                              if (2 * k >= n) k = (n - 1) / 2;
                              const double lo = a[k];
                              const double hi = a[n - 1 - k];
                              long double s = 0.0L;
                              // 先頭k個 → lo
                              s += (long double)k * (long double)lo;
                              // 中央
                              for (size_t i = k; i < n - k; ++i)
                               s += (long double)a[i];
                              // 末尾k個 → hi
                              s += (long double)k * (long double)hi;
                              return (double)(s / (long double)n);
                             }};
  cfg.functions["corrspearman"] = {2, INT_MAX, [](auto &v, auto &ctx) -> Value {
                                    if (v.size() % 2 != 0) throwDomain(ctx.pos);
                                    const size_t n = v.size() / 2;
                                    if (n < 2) throwDomain(ctx.pos);
                                    std::vector<double> x;
                                    std::vector<double> y;
                                    x.reserve(n);
                                    y.reserve(n);
                                    for (size_t i = 0; i < n; ++i)
                                     x.push_back(asReal(v[i], ctx.pos));
                                    for (size_t i = 0; i < n; ++i)
                                     y.push_back(asReal(v[i + n], ctx.pos));
                                    auto rx = rankAverageTies(x, ctx.pos);
                                    auto ry = rankAverageTies(y, ctx.pos);
                                    return pearsonCorr(rx, ry, ctx.pos);
                                   }};
 }

 void registerSpecialTrig(SystemConfig &cfg) {

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
  cfg.functions["sinhc"] = {1, 1, makeDivXCmplxReal([](Complex x) { return std::sinh(x); }, [](double x) { return std::sinh(x); }, 1.0)};
  cfg.functions["tanhc"] = {1, 1, makeDivXCmplxReal([](Complex x) { return std::tanh(x); }, [](double x) { return std::tanh(x); }, 1.0)};
  cfg.functions["expc"] = {1, 1, makeExpc()};
 }

 void registerGeoVec(SystemConfig &cfg) {
  cfg.functions["norm"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.empty()) throw CalcError(CalcErrorType::DomainError, "norm: no elements", ctx.pos);
                            double acc = 0.0;
                            // for (const auto &x : v)
                            //  acc += std::pow(asReal(x, ctx.pos), 2);
                            for (const auto &x : v) {
                             const double r = asReal(x, ctx.pos);
                             acc += r * r; // pow(r,2) より圧倒的に速い
                            }
                            return std::sqrt(acc);
                           }};
  cfg.functions["dot"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.size() % 2 != 0) throw CalcError(CalcErrorType::DomainError, "dot: dimension mismatch", ctx.pos);
                           size_t n = v.size() / 2;
                           if (n == 0) throw CalcError(CalcErrorType::DomainError, "dot: no elements", ctx.pos);

                           double acc = 0.0;
                           // for (size_t i = 0; i < n; ++i)
                           //  acc += asReal(v[i], ctx.pos) * asReal(v[i + n], ctx.pos);
                           for (size_t i = 0; i < n; ++i)
                            acc += asReal(v[i], ctx.pos) * asReal(v[i + n], ctx.pos);

                           return acc;
                          }};
  cfg.functions["lerp"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                            double a = asReal(v[0], ctx.pos);
                            double b = asReal(v[1], ctx.pos);
                            double t = asReal(v[2], ctx.pos);
                            return a * (1 - t) + b * t;
                            // return std::lerp(a, b, t); // leapは丸めが微妙なことがあった
                           }};
  cfg.functions["distance"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                if (v.size() % 2 != 0) throw CalcError(CalcErrorType::DomainError, "distance: dimension mismatch", ctx.pos);
                                size_t n = v.size() / 2;
                                if (n == 0) throw CalcError(CalcErrorType::DomainError, "distance: no elements", ctx.pos);
                                double acc = 0.0;
                                // for (size_t i = 0; i < n; ++i) {
                                //  acc += std::pow(asReal(v[i + n], ctx.pos) - asReal(v[i], ctx.pos), 2);
                                // }
                                for (size_t i = 0; i < n; ++i) {
                                 const double d = asReal(v[i + n], ctx.pos) - asReal(v[i], ctx.pos);
                                 acc += d * d; // pow(d,2) より圧倒的に速い
                                }
                                return std::sqrt(acc);
                               }};
 }
 void registerComplex(SystemConfig &cfg) {
  cfg.functions["re"] = {1, 1, [](auto &v, auto &ctx) -> Value { return isComplex(v[0]) ? std::real(asComplex(v[0])) : asReal(v[0], ctx.pos); }};
  cfg.functions["real"] = cfg.functions["re"];
  cfg.functions["im"] = {1, 1, [](auto &v, auto &ctx) -> Value { return isComplex(v[0]) ? std::imag(asComplex(v[0])) : 0.0; }};
  cfg.functions["imag"] = cfg.functions["im"];
  cfg.functions["arg"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           Complex z = isComplex(v[0]) ? asComplex(v[0]) : Complex(asReal(v[0], ctx.pos), 0.0);
                           if (z == Complex(0, 0)) throwDomain(ctx.pos);
                           return toDeg(std::atan2(std::imag(z), std::real(z)));
                          }};
  cfg.functions["conj"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            Complex z = isComplex(v[0]) ? asComplex(v[0]) : Complex(asReal(v[0], ctx.pos), 0.0);
                            return std::conj(z);
                           }};
  cfg.functions["polar"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             double r = asReal(v[0], ctx.pos);
                             double th = asReal(v[1], ctx.pos);
                             return Complex(r * std::cos(toRad(th)), r * std::sin(toRad(th)));
                            }};
  cfg.functions["rect"] = cfg.functions["polar"];
  cfg.functions["cis"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           double th = asReal(v[0], ctx.pos);
                           return Complex(std::cos(toRad(th)), std::sin(toRad(th)));
                          }};
  cfg.functions["proj"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            Complex z = isComplex(v[0]) ? asComplex(v[0]) : Complex(asReal(v[0], ctx.pos), 0.0);
                            return std::proj(z);
                           }};
  cfg.functions["unit"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            Complex z = isComplex(v[0]) ? asComplex(v[0]) : Complex(asReal(v[0], ctx.pos), 0.0);
                            double a = std::abs(z);
                            if (a == 0.0) throw CalcError(CalcErrorType::DomainError, "unit(0) is undefined", ctx.pos);
                            return z / a;
                           }};

  cfg.functions["csgn"] = cfg.functions["unit"];
  cfg.functions["mag"] = cfg.functions["abs"];
 }
 void registerRandom(SystemConfig &cfg) {
  cfg.functions["rand"] = {0, 2, [](auto &v, auto &ctx) -> Value {
                            static thread_local std::mt19937_64 rng{std::random_device{}()};
                            auto make = [&](double lo, double hi) -> Value {
                             checkFinite(lo, ctx.pos);
                             checkFinite(hi, ctx.pos);
                             if (lo > hi) std::swap(lo, hi);
                             if (lo == hi) return lo; // 退化ケース
                             std::uniform_real_distribution<double> dist(lo, hi);
                             return dist(rng);
                            };
                            if (v.empty()) return make(0.0, 1.0);                       // 引数なし → [0, 1)
                            if (v.size() == 1) return make(0.0, asReal(v[0], ctx.pos)); // 引数1つ → [0, hi)
                            return make(asReal(v[0], ctx.pos), asReal(v[1], ctx.pos));  // 引数2つ → [lo, hi)
                           }};

  cfg.functions["randint"] = {0, 2, [](auto &v, auto &ctx) -> Value {
                               static thread_local std::mt19937_64 rng{std::random_device{}()};
                               auto gen = [&](long long lo, long long hi) -> Value {
                                if (lo > hi) std::swap(lo, hi);
                                if (lo == hi) return static_cast<double>(lo);
                                std::uniform_int_distribution<long long> dist(lo, hi);
                                return static_cast<double>(dist(rng));
                               };

                               if (v.empty()) return gen(0, 1);
                               if (v.size() == 1) return gen(0, requireInt(v[0], ctx.pos));
                               return gen(requireInt(v[0], ctx.pos), requireInt(v[1], ctx.pos));
                              }};

  cfg.functions["randn"] = {0, 2, [](auto &v, auto &ctx) -> Value {
                             static thread_local std::mt19937_64 rng{std::random_device{}()};
                             // 正規分布に従う乱数
                             double mu = (v.size() >= 1) ? asReal(v[0], ctx.pos) : 0.0;
                             double sigma = (v.size() >= 2) ? asReal(v[1], ctx.pos) : 1.0;
                             if (sigma < 0.0) throwDomain(ctx.pos); // sigmaは負にできない
                             std::normal_distribution<double> dist(mu, sigma);
                             return dist(rng);
                            }};

  cfg.functions["choice"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                              static thread_local std::mt19937_64 rng{std::random_device{}()};
                              if (v.empty()) throw CalcError(CalcErrorType::InvalidArgument, "choice: no arguments", ctx.pos);
                              std::uniform_int_distribution<size_t> dist(0, v.size() - 1);
                              return v[dist(rng)];
                             }};

  // fma(a, b, c)
  cfg.functions["fma"] = {3, 3, [](auto &v, auto &ctx) -> Value { return std::fma(asReal(v[0], ctx.pos), asReal(v[1], ctx.pos), asReal(v[2], ctx.pos)); }};
 }

 void registerEngineering(SystemConfig &cfg) { // stress(F, A) = F/A
  cfg.functions["stress"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              double F = asReal(v[0], ctx.pos);
                              double A = asReal(v[1], ctx.pos);
                              if (A == 0.0) throwDomain(ctx.pos);
                              return F / A;
                             }};

  // strain(dL, L) = dL/L
  cfg.functions["strain"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              double dL = asReal(v[0], ctx.pos);
                              double L = asReal(v[1], ctx.pos);
                              if (L == 0.0) throwDomain(ctx.pos);
                              return dL / L;
                             }};

  // young(sigma, eps) = sigma/eps
  cfg.functions["young"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             double sigma = asReal(v[0], ctx.pos);
                             double eps = asReal(v[1], ctx.pos);
                             if (eps == 0.0) throwDomain(ctx.pos);
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
                                   if (d == 0.0) throwDomain(ctx.pos);
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
                                           if (K == 0.0 || d == 0.0) throwDomain(ctx.pos);
                                           return T / (K * d);
                                          }};

  // friction(mu, N) = mu*N
  cfg.functions["friction"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                double mu = asReal(v[0], ctx.pos);
                                double N = asReal(v[1], ctx.pos);
                                return mu * N;
                               }};
 }
 void registerOthers(SystemConfig &cfg) {
  cfg.functions["cnst"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            // v[0] が string であることを要求
                            const std::string &name = asString(v[0], ctx.pos);
                            auto it = constants_dic.find(name);
                            if (it == constants_dic.end()) { throw CalcError(CalcErrorType::UnknownIdentifier, "unknown constant: " + name, ctx.pos); }
                            return it->second;
                           }};
  cfg.functions["fib"] = {1, 1, [](auto &v, auto &ctx) -> Value { return (double)fibULL(requireInt(v[0], ctx.pos), ctx.pos); }};
  cfg.functions["isprime"] = {1, 1, [](auto &v, auto &ctx) -> Value { return isPrimeLL(requireInt(v[0], ctx.pos)) ? 1.0 : 0.0; }};
  cfg.functions["if"] = {3, 3, [](auto &v, auto &ctx) -> Value { return asReal(v[0], ctx.pos) != 0.0 ? v[1] : v[2]; }};
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
 ////!!!

 void initFunctions(SystemConfig &cfg) {
  registerBasicMath(cfg);   // 基本関数
  registerExpLog(cfg);      // log系
  registerTrig(cfg);        // 三角関数
  registerHyperbolic(cfg);  // hyper三角関数
  registerGeoVec(cfg);      // 図形, ベクトル
  registerSpecialTrig(cfg); // hyper三角関数
  registerStatistics(cfg);  // 統計系
  registerComplex(cfg);     // 複素
  registerRandom(cfg);      // 乱数
  registerEngineering(cfg); // 機械設計
  registerOthers(cfg);      // アミューズ
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

} // namespace mm::cal