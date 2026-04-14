#include "functions/functions.hpp"

namespace mm::cal {

 void registerStatistics(SystemConfig &cfg) {
  cfg.functions["fact"] = {1, 1, [](auto &v, auto &ctx) -> Value { return (double)factLD(requireInt(v[0], ctx.pos), ctx.pos); }};
  cfg.functions["gcd"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.size() < 2) throwDomain(ctx.pos, "need at least 2 arguments");
                           long long g = requireInt(v[0], ctx.pos);
                           for (size_t i = 1; i < v.size(); ++i) {
                            g = gcdLL(g, requireInt(v[i], ctx.pos));
                           }
                           return (double)g;
                          }};
  cfg.functions["comb"] = {2, 2, [](auto &v, auto &ctx) -> Value { return (double)combLL(requireInt(v[0], ctx.pos), requireInt(v[1], ctx.pos), ctx.pos); }};
  cfg.functions["perm"] = {2, 2, [](auto &v, auto &ctx) -> Value { return (double)permLL(requireInt(v[0], ctx.pos), requireInt(v[1], ctx.pos), ctx.pos); }};
  cfg.functions["lcm"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.size() < 2) throwDomain(ctx.pos, "need at least 2 arguments");
                           long long l = requireInt(v[0], ctx.pos);
                           for (size_t i = 1; i < v.size(); ++i) {
                            long long b = requireInt(v[i], ctx.pos);
                            if (l == 0 || b == 0) return 0.0;
                            l = checkedMul(l / gcdLL(l, b), b, ctx.pos);
                           }
                           return (double)l;
                          }};
  cfg.functions["sum"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                           double sum_real = 0.0;
                           double sum_imag = 0.0;

                           for (const auto &arg : v) { // MultiValue 展開
                            if (arg.isMulti()) {
                             const auto &mv = arg.asMultiRef(ctx.pos);
                             for (const auto &elem : mv.elems()) {
                              auto s = elem.asComplex(ctx.pos);
                              sum_real += s.real();
                              sum_imag += s.imag();
                             }
                            } else {
                             auto s = arg.asComplex(ctx.pos);
                             sum_real += s.real();
                             sum_imag += s.imag();
                            }
                           }
                           if (std::abs(sum_imag) < 1e-12) return Value(sum_real);
                           return Value(std::complex<double>(sum_real, sum_imag));
                          }};
  cfg.functions["prod"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.size() == 1 && v[0].isMulti()) {
                             const auto &mv = v[0].asMultiRef(ctx.pos);
                             Complex acc = 1;
                             for (const auto &elem : mv.elems()) {
                              acc *= elem.toComplex(ctx.pos);
                             }
                             return acc;
                            } else {
                             auto x = collectComplex(v, ctx);
                             Complex acc = 1;
                             for (auto d : x)
                              acc *= d;
                             return acc;
                            }
                           }};
  cfg.functions["mean"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.size() == 1 && v[0].isMulti()) {
                             const auto &mv = v[0].asMultiRef(ctx.pos);
                             Complex acc = 0;
                             for (const auto &elem : mv.elems()) {
                              acc += elem.toComplex(ctx.pos);
                             }
                             return (acc / static_cast<double>(mv.size()));
                            } else {
                             auto x = collectComplex(v, ctx);
                             if (x.empty()) throwDomain(ctx.pos);
                             Complex acc = 0.0;
                             for (auto d : x)
                              acc += d;
                             return (acc / static_cast<double>(x.size()));
                            }
                           }};
  cfg.functions["mod"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           const double x = v[0].asScalar(ctx.pos);
                           const double y = v[1].asScalar(ctx.pos);
                           if (y == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), ctx.pos);
                           return x - y * std::floor(x / y);
                          }};

  cfg.functions["geomean"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                               if (v.empty()) throwDomain(ctx.pos);
                               long double sumLog = 0.0L;
                               for (const auto &x : v) {
                                if (isComplex(x)) throw CalcError(CalcErrorType::DomainError, "DomainError: geomean: complex", ctx.pos);
                                const double a = x.asScalar(ctx.pos);
                                if (a < 0.0) throw CalcError(CalcErrorType::DomainError, "DomainError: geomean: negative", ctx.pos);
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
                                 if (isComplex(x)) throw CalcError(CalcErrorType::DomainError, "DomainError: harmmean: complex", ctx.pos);
                                 const double a = x.asScalar(ctx.pos);
                                 if (a == 0.0) throw CalcError(CalcErrorType::DivisionByZero, "DivisionByZero: harmmean: zero element", ctx.pos);
                                 acc += 1.0L / (long double)a;
                                }
                                const long double r = (long double)v.size() / acc;
                                if (!std::isfinite((double)r)) throwOverflow(ctx.pos);
                                return (double)r;
                               }};

  cfg.functions["quantile"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                const double p = v[0].asScalar(ctx.pos);
                                if (!(p >= 0.0 && p <= 1.0)) throw CalcError(CalcErrorType::DomainError, "DomainError: quantile: p out of range", ctx.pos);
                                if (v.size() < 2) throw CalcError(CalcErrorType::DomainError, "DomainError: quantile: no samples", ctx.pos);
                                std::vector<double> a;
                                a.reserve(v.size() - 1);
                                for (size_t i = 1; i < v.size(); ++i)
                                 a.push_back(v[i].asScalar(ctx.pos));
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
                           double m = v[0].asScalar(ctx.pos);
                           for (size_t i = 1; i < v.size(); ++i)
                            m = std::min(m, v[i].asScalar(ctx.pos));
                           return m;
                          }};
  cfg.functions["max"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                           if (v.empty()) throwDomain(ctx.pos);
                           double m = v[0].asScalar(ctx.pos);
                           for (size_t i = 1; i < v.size(); ++i)
                            m = std::max(m, v[i].asScalar(ctx.pos));
                           return m;
                          }};
  cfg.functions["clamp"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                             const double x = v[0].asScalar(ctx.pos);
                             const double lo = v[1].asScalar(ctx.pos);
                             const double hi = v[2].asScalar(ctx.pos);
                             if (lo > hi) throw CalcError(CalcErrorType::DomainError, "DomainError: clamp: lo > hi", ctx.pos);
                             return std::clamp(x, lo, hi);
                            }};
  cfg.functions["fract"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             double x = v[0].asScalar(ctx.pos);
                             return x - std::floor(x);
                            }};
  cfg.functions["gamma"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             const double x = v[0].asScalar(ctx.pos);
                             // 負の整数は極
                             if (x < 0 && std::floor(x) == x) throw CalcError(CalcErrorType::DomainError, "DomainError: gamma: pole at negative integer", ctx.pos);
                             const double r = std::tgamma(x);
                             if (std::isnan(r)) throwDomain(ctx.pos);
                             if (!std::isfinite(r)) throwOverflow(ctx.pos);
                             return r;
                            }};
  cfg.functions["lgamma"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              const double x = v[0].asScalar(ctx.pos);
                              const double r = std::lgamma(x);
                              if (std::isnan(r)) throwDomain(ctx.pos);
                              if (!std::isfinite(r)) throwOverflow(ctx.pos);
                              return r;
                             }};
  // 反射対応
  cfg.functions["digamma"] = {1, 1, [&](auto &v, auto &ctx) -> Value {
                               double x = v[0].asScalar(ctx.pos);
                               if (x <= 0.0 && std::floor(x) == x) throwDomain(ctx.pos, "digamma: pole");

                               // reflection formula
                               if (x < 0.0) {
                                double rec = evaluateFunction("digamma", {Value(1.0 - x)}, ctx).asScalar(ctx.pos);
                                return rec - PI / std::tan(PI * x);
                               }
                               double result = 0.0;
                               while (x < 6.0) {
                                result -= 1.0 / x;
                                x += 1.0;
                               }
                               double f = 1.0 / (x * x);
                               result += std::log(x) - 0.5 / x - f * (1.0 / 12 - f * (1.0 / 120 - f / 252));
                               return result;
                              }};

  cfg.functions["zeta"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double s = v[0].asScalar(ctx.pos);
                            if (s <= 1.0) throwDomain(ctx.pos, "zeta: s must be > 1");
                            double r = zetaEulerMaclaurin(s);
                            if (!std::isfinite(r)) throwOverflow(ctx.pos);
                            return r;
                           }};

  cfg.functions["trigamma"] = {1, 1, [&](auto &v, auto &ctx) -> Value {
                                double x = v[0].asScalar(ctx.pos);
                                if (x <= 0.0 && std::floor(x) == x) throwDomain(ctx.pos, "trigamma: pole");
                                // reflection formula
                                if (x < 0.0) {
                                 Value rec = evaluateFunction("trigamma", {Value(1.0 - x)}, ctx);
                                 double s = std::sin(PI * x);
                                 double term = (PI * PI) / (s * s);
                                 return Value(rec.asScalar(ctx.pos) + term);
                                }
                                double result = 0.0;
                                // asymptotic reduction
                                while (x < 6.0) {
                                 result += 1.0 / (x * x);
                                 x += 1.0;
                                }
                                result += 1.0 / (2.0 * x * x) + (1.0 + (1.0 / (6.0 * x * x)) * (1.0 - 1.0 / (30.0 * x * x))) / x;
                                return result;
                               }};

  // std::tgamma(x) * std::tgamma(y) / std::tgamma(x + y)だとすぐ漏れちゃうから対策
  cfg.functions["beta"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            double x = v[0].asScalar(ctx.pos);
                            double y = v[1].asScalar(ctx.pos);
                            if (x <= 0.0 || y <= 0.0) throwDomain(ctx.pos);
                            double logb = std::lgamma(x) + std::lgamma(y) - std::lgamma(x + y);
                            double r = std::exp(logb);
                            if (!std::isfinite(r)) throwOverflow(ctx.pos);
                            return r;
                           }};

  cfg.functions["ibeta"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                             double a = v[0].asScalar(ctx.pos);
                             double b = v[1].asScalar(ctx.pos);
                             double x = v[2].asScalar(ctx.pos);

                             if (x < 0.0 || x > 1.0 || a <= 0.0 || b <= 0.0) throwDomain(ctx.pos);
                             if (x == 0.0) return 0.0;
                             if (x == 1.0) return 1.0;
                             double bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + a * std::log(x) + b * std::log(1.0 - x));
                             double result;
                             if (x < (a + 1.0) / (a + b + 2.0)) result = bt * betacf(a, b, x) / a;
                             else result = 1.0 - bt * betacf(b, a, 1.0 - x) / b;
                             if (!std::isfinite(result)) throwOverflow(ctx.pos);
                             return result;
                            }};

  cfg.functions["numdiff"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                               auto x = collectReals(v, ctx);
                               if (x.size() < 2) throwDomain(ctx.pos, "need at least 2 samples");

                               double sumsq = 0;
                               for (size_t i = 0; i < x.size() - 1; ++i) {
                                double d = x[i + 1] - x[i];
                                sumsq += d * d;
                               }

                               return std::sqrt(sumsq / (x.size() - 1));
                              }};

  cfg.functions["mode"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                            auto a = collectReals(v, ctx);
                            if (a.empty()) throwDomain(ctx.pos);
                            std::sort(a.begin(), a.end());
                            double best = a[0], cur = a[0];
                            int bestCnt = 1, curCnt = 1;
                            for (size_t i = 1; i < a.size(); ++i) {
                             if (nearly_equal(a[i], cur)) {
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
                            return best;                      // 同率なら最小値が残る
                           }};

  cfg.functions["erf"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::erf(v[0].asScalar(ctx.pos)); }};
  cfg.functions["erfc"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::erfc(v[0].asScalar(ctx.pos)); }};

  cfg.functions["round"] = {1, 2, [](auto &v, auto &ctx) -> Value {
                             const double x = v[0].asScalar(ctx.pos);
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
                              if (v.empty()) throw CalcError(CalcErrorType::DomainError, "DomainError: median: no elements", ctx.pos);
                              auto a = collectReals(v, ctx);
                              return medianInplace(a);
                             }};
  // MAD (中央値絶対偏差)
  cfg.functions["mad"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                           if (v.empty()) throw CalcError(CalcErrorType::DomainError, "DomainError: mad: no elements", ctx.pos);
                           auto a = collectReals(v, ctx);
                           const double med = medianInplace(a);
                           for (double &x : a)
                            x = std::abs(x - med);
                           return medianInplace(a);
                          }};
  // MAD (平均絶対偏差)
  cfg.functions["madR"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            auto x = collectReals(v, ctx);
                            if (x.empty()) throwDomain(ctx.pos);
                            double mean = 0;
                            for (auto d : x)
                             mean += d;
                            mean /= x.size();
                            double mad = 0;
                            for (auto d : x)
                             mad += std::abs(d - mean);
                            return mad / x.size();
                           }};
  cfg.functions["skew"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value { return skewnessReal(v, ctx); }};
  cfg.functions["kurtp"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value { return kurtosisExcessReal(v, ctx); }};
  cfg.functions["kurts"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value { return kurtosisExcessSample(v, ctx); }};
  cfg.functions["percentile"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                  const double p = v[0].asScalar(ctx.pos);
                                  if (!(p >= 0.0 && p <= 100.0)) throw CalcError(CalcErrorType::DomainError, "DomainError: percentile: p out of range", ctx.pos);
                                  if (v.size() < 2) throw CalcError(CalcErrorType::DomainError, "DomainError: percentile: no samples", ctx.pos);
                                  std::vector<double> a;
                                  a.reserve(v.size() - 1);
                                  for (size_t i = 1; i < v.size(); ++i)
                                   a.push_back(v[i].asScalar(ctx.pos));
                                  return quantileLinearSorted(std::move(a), p / 100.0, ctx.pos);
                                 }};
  cfg.functions["cov"] = {2, -1, [](auto &v, auto &ctx) -> Value { return covariancePopulationReal(v, ctx); }};
  cfg.functions["corr"] = {2, -1, [](auto &v, auto &ctx) -> Value { return corrPopulationReal(v, ctx); }};
  cfg.functions["ave"] = cfg.functions["mean"];
  cfg.functions["rms"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                           long double ss = 0.0L;
                           size_t n = 0;
                           for (auto &xv : v) {
                            double x = xv.asScalar(ctx.pos);
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
                          if (!std::isfinite(mu) || mu == 0.0) throwDomain(ctx.pos, "DomainError: cv: mean must be nonzero");

                          // population variance
                          const long double var = s.m2 / (long double)s.n;
                          if (!(var >= 0.0L)) throw CalcError(CalcErrorType::DomainError, "DomainError: cv: invalid variance", ctx.pos);

                          const double sd = std::sqrt((double)var);
                          if (!std::isfinite(sd)) throwOverflow(ctx.pos);

                          return sd / mu;
                         }};

  cfg.functions["stderr"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                              if (v.size() < 2) throwDomain(ctx.pos, "DomainError: stderr: need at least 2 samples");
                              const auto s = welfordMeanM2Real(v, ctx);
                              // population variance
                              const long double var = s.m2 / (long double)s.n;
                              if (!(var >= 0.0L)) throw CalcError(CalcErrorType::DomainError, "DomainError: stderr: invalid variance", ctx.pos);
                              const double sd = std::sqrt((double)var);
                              const double se = sd / std::sqrt((double)s.n);

                              if (!std::isfinite(se)) throwOverflow(ctx.pos);
                              return se;
                             }};

  cfg.functions["zscore"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                              double x = v[0].asScalar(ctx.pos);
                              double mu = v[1].asScalar(ctx.pos);
                              double sigma = v[2].asScalar(ctx.pos);
                              if (sigma == 0.0) throwDomain(ctx.pos, "sigma must be nonzero");
                              return (x - mu) / sigma;
                             }};
  cfg.functions["iqr"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                           auto a = gatherReals(v, ctx.pos);
                           if (a.size() < 2) throwDomain(ctx.pos, "need at least 2 samples");
                           std::sort(a.begin(), a.end());
                           return quantileLinearSorted(a, 0.75, ctx.pos) - quantileLinearSorted(a, 0.25, ctx.pos); // quantileLinearでソート回数を減らすやつ
                          }};
  cfg.functions["trimmean"] = {2, INT_MAX, [](auto &v, auto &ctx) -> Value {
                                const double p = v[0].asScalar(ctx.pos);
                                if (!(p >= 0.0 && p < 0.5)) throwDomain(ctx.pos);
                                std::vector<double> a;
                                a.reserve(v.size() - 1);
                                for (size_t i = 1; i < v.size(); ++i)
                                 a.push_back(v[i].asScalar(ctx.pos));

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
                              const double p = v[0].asScalar(ctx.pos);
                              if (!(p >= 0.0 && p < 0.5)) throwDomain(ctx.pos);

                              std::vector<double> a;
                              a.reserve(v.size() - 1);
                              for (size_t i = 1; i < v.size(); ++i)
                               a.push_back(v[i].asScalar(ctx.pos));

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
                                     x.push_back(v[i].asScalar(ctx.pos));
                                    for (size_t i = 0; i < n; ++i)
                                     y.push_back(v[i + n].asScalar(ctx.pos));
                                    auto rx = rankAverageTies(x, ctx.pos);
                                    auto ry = rankAverageTies(y, ctx.pos);
                                    return pearsonCorr(rx, ry, ctx.pos);
                                   }};
  cfg.functions["polylog"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                               double s = v[0].asScalar(ctx.pos);
                               double z = v[1].asScalar(ctx.pos);
                               if (std::abs(z) >= 1.0) throwDomain(ctx.pos, "polylog: |z| must be < 1");
                               const int MAXIT = 100000;
                               double sum = 0.0;
                               double term = z;
                               for (int n = 1; n < MAXIT; ++n) {
                                double add = term / std::pow((double)n, s);
                                sum += add;
                                if (std::abs(add) < cnst_precision_inv) break;
                                term *= z;
                               }
                               if (!std::isfinite(sum)) throwOverflow(ctx.pos);
                               return sum;
                              }};
 }
} // namespace mm::cal
