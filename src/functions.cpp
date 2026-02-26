#include "functions.hpp"
#include "constants.hpp"
#include "core.hpp"
#include "evaluate.hpp"
#include "math_util.hpp"

namespace mm::cal {

 /* ============================
   地獄の無限関数定義編
   ============================ */

 class FunctionBuilder {
   SystemConfig &cfg;

  public:
   FunctionBuilder(SystemConfig &c) : cfg(c) {}

   template <typename F> void def(const std::string &name, int min, int max, F f) { cfg.functions[name] = {min, max, f}; }
 };

 ///!!!!
 static void registerBasicMath(SystemConfig &cfg);     // 基本関数
 static void registerExpLog(SystemConfig &cfg);        // log系
 static void registerTrig(SystemConfig &cfg);          // 三角関数
 static void registerHyperbolic(SystemConfig &cfg);    // hyper三角関数
 static void registerGeoVec(SystemConfig &cfg);        // 図形, ベクトル
 static void registerVector(SystemConfig &cfg);        // 行列・線形代数
 static void registerSpecialTrig(SystemConfig &cfg);   // hyper三角関数
 static void registerStatistics(SystemConfig &cfg);    // 統計系
 static void registerComplex(SystemConfig &cfg);       // 複素
 static void registerRandom(SystemConfig &cfg);        // 乱数
 static void registerAreaVol(SystemConfig &cfg);       // 面積，体積
 static void registerEngineering(SystemConfig &cfg);   // 機械設計
 static void registerMoldInjection(SystemConfig &cfg); // 金型系
 static void registerFinance(SystemConfig &cfg);       // 財務系
 static void registerOthers(SystemConfig &cfg);        // アミューズ

 void registerBasicMath(SystemConfig &cfg) { // 基本関数
  cfg.functions["abs"] = {1, 1, [](auto &v, auto &ctx) -> Value { return applyUnaryNumeric(v[0], [](double d) { return std::fabs(d); }, [](Complex c) { return std::abs(c); }, ctx.pos); }};

  cfg.functions["sign"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = requireReal(v[0], ctx.pos);
                            return (x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0;
                           }};

  cfg.functions["sqrt"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            return applyUnaryNumeric(
                                v[0],
                                [](double d) -> Value {
                                 if (d < 0) return Complex(0, std::sqrt(-d));
                                 return std::sqrt(d);
                                },
                                [](Complex c) { return std::sqrt(c); }, ctx.pos);
                           }};

  cfg.functions["cbrt"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::cbrt(requireReal(v[0], ctx.pos)); }};
  cfg.functions["floor"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::floor(requireReal(v[0], ctx.pos)); }};
  cfg.functions["ceil"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::ceil(requireReal(v[0], ctx.pos)); }};
  cfg.functions["trunc"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::trunc(requireReal(v[0], ctx.pos)); }};
  cfg.functions["pow"] = {2, 2, [](auto &v, auto &ctx) -> Value { return realIfPossible(std::pow(toComplex(v[0], ctx.pos), toComplex(v[1], ctx.pos))); }};

  cfg.functions["nextpow2"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                double x = requireReal(v[0], ctx.pos);
                                if (x <= 0.0) throwDomain(ctx.pos);

                                int n = 0;
                                double val = 1.0;
                                while (val < x) {
                                 val *= 2.0;
                                 ++n;
                                }
                                return n;
                               }};

  // ---- convert angle ----
  cfg.functions["DtoR"] = {1, 1, [](auto &v, auto &ctx) -> Value { return deg2rad_v(v[0], ctx.pos); }};
  cfg.functions["DtoG"] = {1, 1, [](auto &v, auto &ctx) -> Value { return deg2grad_v(v[0], ctx.pos); }};
  cfg.functions["RtoD"] = {1, 1, [](auto &v, auto &ctx) -> Value { return rad2deg_v(v[0], ctx.pos); }};
  cfg.functions["RtoG"] = {1, 1, [](auto &v, auto &ctx) -> Value { return rad2grad_v(v[0], ctx.pos); }};
  cfg.functions["GtoD"] = {1, 1, [](auto &v, auto &ctx) -> Value { return grad2deg_v(v[0], ctx.pos); }};
  cfg.functions["GtoR"] = {1, 1, [](auto &v, auto &ctx) -> Value { return grad2rad_v(v[0], ctx.pos); }};
 }

 void registerExpLog(SystemConfig &cfg) { // log, exp系
  cfg.functions["exp"] = {1, 1, [](auto &v, auto &ctx) -> Value { return realIfPossible(std::exp(toComplex(v[0], ctx.pos))); }};
  cfg.functions["log10"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             double x = requireReal(v[0], ctx.pos);
                             if (x <= 0.0) throwDomain(ctx.pos);
                             return std::log10(x);
                            }};
  cfg.functions["log2"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = requireReal(v[0], ctx.pos);
                            if (x <= 0.0) throwDomain(ctx.pos);
                            return std::log2(x);
                           }};
  cfg.functions["log"] = {1, 2, [](const std::vector<Value> &args, const FunctionContext &ctx) -> Value {
                           auto ln_real = [](double x) { return std::log(x); };
                           auto ln_complex = [](std::complex<double> z) { return std::log(z); };

                           // --------------------------
                           // 1引数: ln(x)
                           // --------------------------
                           if (args.size() == 1) {
                            const Value &v = args[0];

                            // realの場合
                            if (v.isScalar()) {
                             double x = v.asScalar(ctx.pos);

                             if (x > 0.0) return ln_real(x);

                             if (x == 0.0) throwDomain(ctx.pos); // ln(0)はエラー

                             // x < 0 → complex必要
                             std::complex<double> z(x, 0.0);
                             return ln_complex(z);
                            }

                            // complex
                            std::complex<double> z = v.toComplex(ctx.pos);
                            return ln_complex(z);
                           }

                           // --------------------------
                           // 2引数: log(base, x)
                           // --------------------------

                           std::complex<double> base = args[0].toComplex(ctx.pos);
                           std::complex<double> x = args[1].toComplex(ctx.pos);

                           // log(?,0) は未定義
                           if (x == std::complex<double>(0.0, 0.0)) throwDomain(ctx.pos);
                           if (base == std::complex<double>(0.0, 0.0) || base == std::complex<double>(1.0, 0.0)) throwDomain(ctx.pos);

                           return std::log(x) / std::log(base);
                          }};
  cfg.functions["ln"] = cfg.functions["log"];
  cfg.functions["log1p"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (v[0].isComplex()) {
                              Complex x = toComplex(v[0], ctx.pos);
                              return realIfPossible(std::log(Complex(1, 0) + x));
                             }
                             double x = requireReal(v[0], ctx.pos);
                             if (x <= -1.0) throw CalcError(CalcErrorType::DomainError, "log1p: x <= -1", ctx.pos);

                             return std::log1p(x);
                            }};
  cfg.functions["expm1"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (v[0].isComplex()) {
                              Complex x = toComplex(v[0], ctx.pos);
                              return realIfPossible(std::exp(x) - Complex(1, 0));
                             }
                             double x = requireReal(v[0], ctx.pos);
                             return std::expm1(x);
                            }};
 }

 void registerTrig(SystemConfig &cfg) {
  cfg.functions["sin"] = {1, 1, [=](auto &v, auto &ctx) -> Value { return realIfPossible(std::sin(toComplex(v[0], ctx.pos) * deg2rad)); }};
  cfg.functions["cos"] = {1, 1, [=](auto &v, auto &ctx) -> Value { return realIfPossible(std::cos(toComplex(v[0], ctx.pos) * deg2rad)); }};
  // cfg.functions["tan"] = {1, 1, [=](auto &v, auto &) -> Value { return std::tan(asComplex(v[0]) * deg2rad); }};
  cfg.functions["tan"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (v[0].isComplex()) return realIfPossible(std::tan(toComplex(v[0], ctx.pos) * deg2rad));
                           double r = requireReal(v[0], ctx.pos) * deg2rad;
                           double c = std::cos(r);
                           return (std::abs(c) < eps) ? signedInfBy(std::sin(r)) : std::tan(r);
                          }};
  cfg.functions["cot"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (v[0].isComplex()) {
                            auto z = toComplex(v[0], ctx.pos) * deg2rad;
                            return realIfPossible(std::cos(z) / std::sin(z));
                           }
                           double r = requireReal(v[0], ctx.pos) * deg2rad;
                           double s = std::sin(r);
                           return (std::abs(s) < eps) ? signedInfBy(std::cos(r)) : std::cos(r) / s;
                          }};
  cfg.functions["sec"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (v[0].isComplex()) return realIfPossible(Complex(1, 0) / std::cos(toComplex(v[0], ctx.pos) * deg2rad));
                           double r = requireReal(v[0], ctx.pos) * deg2rad;
                           double c = std::cos(r);
                           return (std::abs(c) < eps) ? signedInfBy(std::sin(r)) : 1.0 / c;
                          }};
  cfg.functions["csc"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (v[0].isComplex()) return realIfPossible(Complex(1, 0) / std::sin(toComplex(v[0], ctx.pos) * deg2rad));
                           double r = requireReal(v[0], ctx.pos) * deg2rad;
                           double s = std::sin(r);
                           return (std::abs(s) < eps) ? signedInfBy(std::cos(r)) : 1.0 / s;
                          }};
  cfg.functions["asin"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (v[0].isComplex()) return realIfPossible(std::asin(toComplex(v[0], ctx.pos)) * rad2deg);
                            double x = requireReal(v[0], ctx.pos);
                            if (std::abs(x) <= 1.0) return std::asin(x) * rad2deg;
                            calcWarn(ctx.cfg, ctx.pos, "asin(|x|>1): complex principal value only");
                            return realIfPossible(std::asin(Complex(x, 0)) * rad2deg);
                           }};
  cfg.functions["acos"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (v[0].isComplex()) return realIfPossible(std::acos(toComplex(v[0], ctx.pos)) * rad2deg);
                            double x = requireReal(v[0], ctx.pos);
                            if (std::abs(x) <= 1.0) return std::acos(x) * rad2deg;
                            calcWarn(ctx.cfg, ctx.pos, "acos(|x|>1): complex principal value only");
                            return realIfPossible(std::acos(Complex(x, 0)) * rad2deg);
                           }};
  cfg.functions["atan"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (v[0].isComplex()) return realIfPossible(std::atan(toComplex(v[0], ctx.pos)) * rad2deg);
                            return std::atan(requireReal(v[0], ctx.pos)) * rad2deg;
                           }};

  cfg.functions["atan2"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             double y = requireReal(v[0], ctx.pos);
                             double x = requireReal(v[1], ctx.pos);
                             if (x == 0.0 && y == 0.0) throwDomain(ctx.pos);
                             return std::atan2(y, x) * rad2deg;
                            }};
 }

 void registerHyperbolic(SystemConfig &cfg) { // hyper三角関数
  cfg.functions["sinh"] = {1, 1, [](auto &v, auto &ctx) -> Value { return realIfPossible(std::sinh(toComplex(v[0], ctx.pos))); }};
  cfg.functions["cosh"] = {1, 1, [](auto &v, auto &ctx) -> Value { return realIfPossible(std::cosh(toComplex(v[0], ctx.pos))); }};
  cfg.functions["tanh"] = {1, 1, [](auto &v, auto &ctx) -> Value { return realIfPossible(std::tanh(toComplex(v[0], ctx.pos))); }};

  // ---- inverse hyperbolic ----
  cfg.functions["asinh"] = {1, 1, [](auto &v, auto &ctx) -> Value { return realIfPossible(std::asinh(toComplex(v[0], ctx.pos))); }};
  cfg.functions["acosh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (v[0].isComplex()) return realIfPossible(std::acosh(toComplex(v[0], ctx.pos)));
                             double x = requireReal(v[0], ctx.pos);
                             if (x >= 1.0) return std::acosh(x);
                             calcWarn(ctx.cfg, ctx.pos, "acosh(x<1): complex principal value only");
                             return realIfPossible(std::acosh(Complex(x, 0)));
                            }};

  cfg.functions["atanh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (v[0].isComplex()) return realIfPossible(std::atanh(toComplex(v[0], ctx.pos)));
                             double x = requireReal(v[0], ctx.pos);
                             if (std::abs(x) < 1.0) return std::atanh(x);
                             calcWarn(ctx.cfg, ctx.pos, "atanh(|x|>=1): complex principal value only");
                             return realIfPossible(std::atanh(Complex(x, 0)));
                            }};
  cfg.functions["csch"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = requireReal(v[0], ctx.pos);
                            if (x == 0.0) throwDomain(ctx.pos, "must be non-zero");
                            return 1.0 / std::sinh(x);
                           }};

  cfg.functions["sech"] = {1, 1, [](auto &v, auto &ctx) -> Value { return 1.0 / std::cosh(requireReal(v[0], ctx.pos)); }};
  cfg.functions["coth"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = requireReal(v[0], ctx.pos);
                            if (x == 0.0) throwDomain(ctx.pos, "must be non-zero");
                            return std::cosh(x) / std::sinh(x);
                           }};
 }

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
                                if (isComplex(x)) throw CalcError(CalcErrorType::DomainError, "geomean: complex", ctx.pos);
                                const double a = x.asScalar(ctx.pos);
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
                                 const double a = x.asScalar(ctx.pos);
                                 if (a == 0.0) throw CalcError(CalcErrorType::DivisionByZero, "harmmean: zero element", ctx.pos);
                                 acc += 1.0L / (long double)a;
                                }
                                const long double r = (long double)v.size() / acc;
                                if (!std::isfinite((double)r)) throwOverflow(ctx.pos);
                                return (double)r;
                               }};

  cfg.functions["quantile"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                const double p = v[0].asScalar(ctx.pos);
                                if (!(p >= 0.0 && p <= 1.0)) throw CalcError(CalcErrorType::DomainError, "quantile: p out of range", ctx.pos);
                                if (v.size() < 2) throw CalcError(CalcErrorType::DomainError, "quantile: no samples", ctx.pos);
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
                             if (lo > hi) throw CalcError(CalcErrorType::DomainError, "clamp: lo > hi", ctx.pos);
                             return std::clamp(x, lo, hi);
                            }};
  cfg.functions["fract"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             double x = v[0].asScalar(ctx.pos);
                             return x - std::floor(x);
                            }};
  cfg.functions["gamma"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             const double x = v[0].asScalar(ctx.pos);
                             // 負の整数は極
                             if (x < 0 && std::floor(x) == x) throw CalcError(CalcErrorType::DomainError, "gamma: pole at negative integer", ctx.pos);
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
                                double s = std::sin(PI * x);
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
                                double f = 1.0 / (x * x);
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
                              if (v.empty()) throw CalcError(CalcErrorType::DomainError, "median: no elements", ctx.pos);
                              auto a = collectReals(v, ctx);
                              return medianInplace(a);
                             }};
  // MAD (中央値絶対偏差)
  cfg.functions["mad"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                           if (v.empty()) throw CalcError(CalcErrorType::DomainError, "mad: no elements", ctx.pos);
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
                                  if (!(p >= 0.0 && p <= 100.0)) throw CalcError(CalcErrorType::DomainError, "percentile: p out of range", ctx.pos);
                                  if (v.size() < 2) throw CalcError(CalcErrorType::DomainError, "percentile: no samples", ctx.pos);
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
                          if (!std::isfinite(mu) || mu == 0.0) throwDomain(ctx.pos, "mean must be nonzero");

                          // population variance
                          const long double var = s.m2 / (long double)s.n;
                          if (!(var >= 0.0L)) throw CalcError(CalcErrorType::DomainError, "cv: invalid variance", ctx.pos);

                          const double sd = std::sqrt((double)var);
                          if (!std::isfinite(sd)) throwOverflow(ctx.pos);

                          return sd / mu;
                         }};

  cfg.functions["stderr"] = {1, INT_MAX, [](auto &v, auto &ctx) -> Value {
                              if (v.size() < 2) throwDomain(ctx.pos, "need at least 2 samples");
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
                                // v[1..] を直接 gather できるならそれが理想
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

  //                          double x = v[0].asScalar(ctx.pos);
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
  cfg.functions["hypot"] = {2, 2, [](auto &v, auto &ctx) -> Value { return std::hypot(requireReal(v[0], ctx.pos), requireReal(v[1], ctx.pos)); }};
  cfg.functions["norm"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.empty()) throw CalcError(CalcErrorType::DomainError, "no elements", ctx.pos);
                            double acc = 0.0;
                            // for (const auto &x : v)
                            //  acc += std::pow(x.asScalar(ctx.pos), 2);
                            for (const auto &x : v) {
                             const double r = x.asScalar(ctx.pos);
                             acc += r * r; // pow(r,2) より圧倒的に速い
                            }
                            return std::sqrt(acc);
                           }};

  // ベクトル和
  cfg.functions["vadd"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vadd requires two multivalue arguments", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            const auto &b = v[1].asMultiRef(ctx.pos);

                            if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);

                            auto result = std::make_shared<MultiValue>();
                            result->elems_.reserve(a.size());
                            for (size_t i = 0; i < a.size(); ++i) {
                             result->elems_.push_back(Value(a[i].asScalar(ctx.pos) + b[i].asScalar(ctx.pos)));
                            }
                            return Value(result);
                           }};

  // ベクトル差
  cfg.functions["vsub"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vsub requires two multivalue arguments", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            const auto &b = v[1].asMultiRef(ctx.pos);

                            if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);

                            auto result = std::make_shared<MultiValue>();
                            result->elems_.reserve(a.size());
                            for (size_t i = 0; i < a.size(); ++i) {
                             result->elems_.push_back(Value(a[i].asScalar(ctx.pos) - b[i].asScalar(ctx.pos)));
                            }
                            return Value(result);
                           }};

  // スカラ倍
  cfg.functions["vscalar"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                               if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vscalar requires multivalue argument", ctx.pos);
                               double scalar = v[1].asScalar(ctx.pos);
                               const auto &vec = v[0].asMultiRef(ctx.pos);

                               auto result = std::make_shared<MultiValue>();
                               result->elems_.reserve(vec.size());
                               for (const auto &elem : vec.elems()) {
                                result->elems_.push_back(Value(elem.asScalar(ctx.pos) * scalar));
                               }
                               return Value(result);
                              }};
  // 内積
  cfg.functions["vdot"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vdot requires two multivalue arguments", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            const auto &b = v[1].asMultiRef(ctx.pos);

                            if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);

                            double acc = 0.0;
                            for (size_t i = 0; i < a.size(); ++i) {
                             acc += a[i].asScalar(ctx.pos) * b[i].asScalar(ctx.pos);
                            }
                            return Value(acc);
                           }};

  // 外積 (3次元限定)
  cfg.functions["vcross"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vcross requires two multivalue arguments", ctx.pos);
                              const auto &a = v[0].asMultiRef(ctx.pos);
                              const auto &b = v[1].asMultiRef(ctx.pos);

                              if (a.size() != 3 || b.size() != 3) throw CalcError(CalcErrorType::DomainError, "vcross requires 3-dimensional vectors", ctx.pos);

                              double x = a[1].asScalar(ctx.pos) * b[2].asScalar(ctx.pos) - a[2].asScalar(ctx.pos) * b[1].asScalar(ctx.pos);
                              double y = a[2].asScalar(ctx.pos) * b[0].asScalar(ctx.pos) - a[0].asScalar(ctx.pos) * b[2].asScalar(ctx.pos);
                              double z = a[0].asScalar(ctx.pos) * b[1].asScalar(ctx.pos) - a[1].asScalar(ctx.pos) * b[0].asScalar(ctx.pos);

                              auto result = std::make_shared<MultiValue>(std::vector<Value>{Value(x), Value(y), Value(z)});
                              return Value(result);
                             }};
  cfg.functions["vnorm"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vnorm requires multivalue argument", ctx.pos);
                             const auto &vec = v[0].asMultiRef(ctx.pos);

                             double sum = 0.0;
                             for (const auto &elem : vec.elems()) {
                              double val = elem.asScalar(ctx.pos);
                              sum += val * val;
                             }
                             return Value(std::sqrt(sum));
                            }};

  // マンハッタン距離
  cfg.functions["vmanhattan"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vmanhattan requires two multivalue arguments", ctx.pos);
                                  const auto &a = v[0].asMultiRef(ctx.pos);
                                  const auto &b = v[1].asMultiRef(ctx.pos);

                                  if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);

                                  double sum = 0.0;
                                  for (size_t i = 0; i < a.size(); ++i) {
                                   sum += std::abs(a[i].asScalar(ctx.pos) - b[i].asScalar(ctx.pos));
                                  }
                                  return Value(sum);
                                 }};

  // ユークリッド距離
  cfg.functions["veuclidean"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "veuclidean requires two multivalue arguments", ctx.pos);
                                  const auto &a = v[0].asMultiRef(ctx.pos);
                                  const auto &b = v[1].asMultiRef(ctx.pos);

                                  if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);

                                  double sum = 0.0;
                                  for (size_t i = 0; i < a.size(); ++i) {
                                   double diff = a[i].asScalar(ctx.pos) - b[i].asScalar(ctx.pos);
                                   sum += diff * diff;
                                  }
                                  return Value(std::sqrt(sum));
                                 }};

  cfg.functions["dot"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "dot requires two multivalue arguments", ctx.pos);

                           const auto &a = v[0].asMultiRef(ctx.pos);
                           const auto &b = v[1].asMultiRef(ctx.pos);
                           if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);
                           if (a.empty()) throw CalcError(CalcErrorType::DomainError, "dot of empty multivalue", ctx.pos);
                           double acc = 0.0;
                           for (size_t i = 0; i < a.size(); ++i) { // nested Multi は禁止するよん
                            if (a[i].isMulti() || b[i].isMulti()) throw CalcError(CalcErrorType::TypeError, "nested multivalue not allowed", ctx.pos);
                            acc += a[i].asScalar(ctx.pos) * b[i].asScalar(ctx.pos);
                           }
                           return Value(acc);
                          }};
  cfg.functions["vsum"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vsum requires multivalue argument", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            double sum = 0.0;
                            for (const auto &elem : a.elems()) {
                             sum += elem.asScalar(ctx.pos);
                            }
                            return Value(sum);
                           }};
  // ベクトルの正規化
  cfg.functions["vnormalize"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vnormalize requires multivalue argument", ctx.pos);
                                  const auto &vec = v[0].asMultiRef(ctx.pos);

                                  double norm = 0.0;
                                  for (const auto &elem : vec.elems()) {
                                   double val = elem.asScalar(ctx.pos);
                                   norm += val * val;
                                  }
                                  norm = std::sqrt(norm);

                                  if (norm == 0.0) throw CalcError(CalcErrorType::DomainError, "cannot normalize zero vector", ctx.pos);

                                  auto result = std::make_shared<MultiValue>();
                                  result->elems_.reserve(vec.size());
                                  for (const auto &elem : vec.elems()) {
                                   result->elems_.push_back(Value(elem.asScalar(ctx.pos) / norm));
                                  }
                                  return Value(result);
                                 }};

  // ベクトルの射影
  cfg.functions["vproject"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vproject requires two multivalue arguments", ctx.pos);
                                const auto &a = v[0].asMultiRef(ctx.pos);
                                const auto &b = v[1].asMultiRef(ctx.pos);

                                if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);

                                // 内積計算
                                double dot_product = 0.0;
                                for (size_t i = 0; i < a.size(); ++i) {
                                 dot_product += a[i].asScalar(ctx.pos) * b[i].asScalar(ctx.pos);
                                }

                                // bのノルムの2乗
                                double b_norm_sq = 0.0;
                                for (size_t i = 0; i < b.size(); ++i) {
                                 double val = b[i].asScalar(ctx.pos);
                                 b_norm_sq += val * val;
                                }

                                if (std::abs(b_norm_sq) < 1e-15) throw CalcError(CalcErrorType::DomainError, "cannot project onto zero vector", ctx.pos);

                                double scalar = dot_product / b_norm_sq;

                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(b.size());
                                for (size_t i = 0; i < b.size(); ++i) {
                                 result->elems_.push_back(Value(b[i].asScalar(ctx.pos) * scalar));
                                }
                                return Value(result);
                               }};

  // 2ベクトル間の角度 (度数法)
  cfg.functions["vangle"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vangle requires two multivalue arguments", ctx.pos);
                              const auto &a = v[0].asMultiRef(ctx.pos);
                              const auto &b = v[1].asMultiRef(ctx.pos);

                              if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);

                              // 内積計算
                              double dot_product = 0.0;
                              for (size_t i = 0; i < a.size(); ++i) {
                               dot_product += a[i].asScalar(ctx.pos) * b[i].asScalar(ctx.pos);
                              }

                              // norm計算
                              double norm_a = 0.0, norm_b = 0.0;
                              for (size_t i = 0; i < a.size(); ++i) {
                               double val_a = a[i].asScalar(ctx.pos);
                               double val_b = b[i].asScalar(ctx.pos);
                               norm_a += val_a * val_a;
                               norm_b += val_b * val_b;
                              }
                              norm_a = std::sqrt(norm_a);
                              norm_b = std::sqrt(norm_b);

                              if (norm_a == 0.0 || norm_b == 0.0) throw CalcError(CalcErrorType::DomainError, "cannot calculate angle with zero vector", ctx.pos);

                              // cos計算
                              double cos_angle = dot_product / (norm_a * norm_b);
                              // cosが-1以上1以下になるようにクランプ
                              cos_angle = std::max(-1.0, std::min(1.0, cos_angle));

                              double angle_rad = std::acos(cos_angle);
                              return Value(angle_rad * 180.0 / PI);
                             }};

  // ベクトルの反射
  cfg.functions["vreflect"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vreflect requires two multivalue arguments", ctx.pos);

                                const auto &a = v[0].asMultiRef(ctx.pos);
                                const auto &n = v[1].asMultiRef(ctx.pos);

                                if (a.size() != n.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);

                                double dot = 0.0;
                                double nn = 0.0;

                                for (size_t i = 0; i < a.size(); ++i) {
                                 double ai = a[i].asScalar(ctx.pos);
                                 double ni = n[i].asScalar(ctx.pos);
                                 dot += ai * ni;
                                 nn += ni * ni;
                                }

                                if (std::abs(nn) < 1e-15) throw CalcError(CalcErrorType::DomainError, "cannot reflect with zero normal vector", ctx.pos);

                                double factor = 2.0 * dot / nn;

                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(a.size());

                                for (size_t i = 0; i < a.size(); ++i) {
                                 double ai = a[i].asScalar(ctx.pos);
                                 double ni = n[i].asScalar(ctx.pos);
                                 result->elems_.push_back(Value(ai - factor * ni));
                                }

                                return Value(result);
                               }};
  // 軸対称反射
  cfg.functions["vreflect_axis"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                     if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vreflect_axis requires two multivalue arguments", ctx.pos);

                                     const auto &a = v[0].asMultiRef(ctx.pos);
                                     const auto &n = v[1].asMultiRef(ctx.pos);

                                     if (a.size() != n.size()) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);

                                     double dot = 0.0;
                                     double nn = 0.0;

                                     for (size_t i = 0; i < a.size(); ++i) {
                                      double ai = a[i].asScalar(ctx.pos);
                                      double ni = n[i].asScalar(ctx.pos);
                                      dot += ai * ni;
                                      nn += ni * ni;
                                     }

                                     if (std::abs(nn) < 1e-15) throw CalcError(CalcErrorType::DomainError, "cannot reflect with zero axis vector", ctx.pos);

                                     double factor = 2.0 * dot / nn;

                                     auto result = std::make_shared<MultiValue>();
                                     result->elems_.reserve(a.size());

                                     for (size_t i = 0; i < a.size(); ++i) {
                                      double ai = a[i].asScalar(ctx.pos);
                                      double ni = n[i].asScalar(ctx.pos);
                                      result->elems_.push_back(Value(factor * ni - ai));
                                     }

                                     return Value(result);
                                    }};

  cfg.functions["vlength"] = cfg.functions["vnorm"];
  cfg.functions["vdistance"] = cfg.functions["veuclidean"];

  // 単位ベクトル
  cfg.functions["vunit"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vunit requires multivalue argument", ctx.pos);
                             const auto &vec = v[0].asMultiRef(ctx.pos);

                             double norm = 0.0;
                             for (const auto &elem : vec.elems()) {
                              double val = elem.asScalar(ctx.pos);
                              norm += val * val;
                             }
                             norm = std::sqrt(norm);

                             if (norm == 0.0) throw CalcError(CalcErrorType::DomainError, "cannot create unit vector from zero vector", ctx.pos);

                             auto result = std::make_shared<MultiValue>();
                             result->elems_.reserve(vec.size());
                             for (const auto &elem : vec.elems()) {
                              result->elems_.push_back(Value(elem.asScalar(ctx.pos) / norm));
                             }
                             return Value(result);
                            }};

  cfg.functions["scalar"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "scalar requires multivalue argument", ctx.pos);
                              double scalar = v[1].asScalar(ctx.pos);
                              const auto &a = v[0].asMultiRef(ctx.pos);

                              auto result = std::make_shared<MultiValue>();
                              result->elems_.reserve(a.size());
                              for (const auto &elem : a.elems()) {
                               result->elems_.push_back(Value(elem.asScalar(ctx.pos) * scalar));
                              }
                              return Value(result);
                             }};

  cfg.functions["lerp"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                            double a = v[0].asScalar(ctx.pos);
                            double b = v[1].asScalar(ctx.pos);
                            double t = v[2].asScalar(ctx.pos);
                            return a * (1 - t) + b * t;
                            // return std::lerp(a, b, t); // leapは丸めが微妙なことがあった
                           }};
  cfg.functions["distance"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                if (v.size() % 2 != 0) throw CalcError(CalcErrorType::DomainError, "dimension mismatch", ctx.pos);
                                size_t n = v.size() / 2;
                                if (n == 0) throw CalcError(CalcErrorType::DomainError, "no elements", ctx.pos);
                                double acc = 0.0;
                                // for (size_t i = 0; i < n; ++i) {
                                //  acc += std::pow(asReal(v[i + n], ctx.pos) - v[i].asScalar(ctx.pos), 2);
                                // }
                                for (size_t i = 0; i < n; ++i) {
                                 const double d = v[i + n].asScalar(ctx.pos) - v[i].asScalar(ctx.pos);
                                 acc += d * d; // pow(d,2) より圧倒的に速い
                                }
                                return std::sqrt(acc);
                               }};
  cfg.functions["totient"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                               int n = (int)v[0].asScalar(ctx.pos);
                               if (n < 1) throwDomain(ctx.pos, "must be >= 1");
                               int phi = n;
                               for (int p = 2; p * p <= n; ++p) {
                                if (n % p == 0) {
                                 while (n % p == 0)
                                  n /= p;
                                 phi -= phi / p;
                                }
                               }
                               if (n > 1) phi -= phi / n;
                               return (double)phi;
                              }};
  cfg.functions["covmatrix"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                 auto x = collectReals(v, ctx);
                                 if (x.size() < 2) throwDomain(ctx.pos);
                                 double mean = 0.0;
                                 for (auto d : x)
                                  mean += d;
                                 mean /= x.size();
                                 double cov = 0.0;
                                 for (auto d : x)
                                  cov += (d - mean) * (d - mean);
                                 return cov / (x.size() - 1);
                                }};
  cfg.functions["corrmatrix"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                  if (v.size() % 2 != 0) throwDomain(ctx.pos, "dimension mismatch");

                                  size_t n = v.size() / 2;
                                  if (n < 2) throwDomain(ctx.pos, "need at least 2 samples");

                                  std::vector<double> x, y;
                                  x.reserve(n);
                                  y.reserve(n);

                                  for (size_t i = 0; i < n; ++i)
                                   x.push_back(requireReal(v[i], ctx.pos));

                                  for (size_t i = 0; i < n; ++i)
                                   y.push_back(requireReal(v[i + n], ctx.pos));

                                  double mx = 0, my = 0;
                                  for (double d : x)
                                   mx += d;
                                  for (double d : y)
                                   my += d;

                                  mx /= n;
                                  my /= n;

                                  double num = 0, dx = 0, dy = 0;

                                  for (size_t i = 0; i < n; ++i) {
                                   double dxi = x[i] - mx;
                                   double dyi = y[i] - my;
                                   num += dxi * dyi;
                                   dx += dxi * dxi;
                                   dy += dyi * dyi;
                                  }

                                  if (dx == 0.0 || dy == 0.0) throwDomain(ctx.pos, "zero variance");

                                  return num / std::sqrt(dx * dy);
                                 }};
  cfg.functions["percentrank"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                   if (v.size() < 2) throwDomain(ctx.pos);

                                   const double val = requireReal(v[0], ctx.pos);

                                   std::vector<double> x;
                                   for (size_t i = 1; i < v.size(); ++i)
                                    x.push_back(requireReal(v[i], ctx.pos));

                                   if (x.empty()) throwDomain(ctx.pos);

                                   size_t count = 0;
                                   for (double d : x)
                                    if (d <= val) ++count;

                                   return (double)count / x.size() * 100.0;
                                  }};

  cfg.functions["winsorR"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                               if (v.size() < 3) throwDomain(ctx.pos);

                               const double p = requireReal(v[0], ctx.pos);
                               if (!(p >= 0.0 && p <= 0.5)) throwDomain(ctx.pos);

                               std::vector<double> x;
                               for (size_t i = 1; i < v.size(); ++i)
                                x.push_back(requireReal(v[i], ctx.pos));

                               if (x.empty()) throwDomain(ctx.pos);

                               std::vector<double> y = x;
                               std::sort(y.begin(), y.end());

                               size_t k = (size_t)(p * y.size());
                               double low = y[k];
                               double high = y[y.size() - 1 - k];

                               long double sum = 0;
                               for (double d : x)
                                sum += std::min(std::max(d, low), high);

                               return (double)(sum / x.size());
                              }};

  // --- Fourier transforms (簡易版) ---
  cfg.functions["fft"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                           auto x = collectReals(v, ctx);
                           size_t N = x.size();
                           if (N == 0) throwDomain(ctx.pos);

                           // 簡易FFT: DFT の振幅平均
                           double amp_sum = 0;
                           for (size_t k = 0; k < N; ++k) {
                            std::complex<double> sum = 0;
                            for (size_t n = 0; n < N; ++n)
                             sum += x[n] * std::exp(std::complex<double>(0, -2.0 * PI * k * n / N));
                            amp_sum += std::abs(sum);
                           }

                           return amp_sum / N;
                          }};
  cfg.functions["ifft"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            auto x = collectReals(v, ctx);
                            size_t N = x.size();
                            if (N == 0) throwDomain(ctx.pos);
                            double amp_sum = 0;
                            for (size_t n = 0; n < N; ++n) {
                             std::complex<double> sum = 0;
                             for (size_t k = 0; k < N; ++k)
                              sum += x[k] * std::exp(std::complex<double>(0, 2.0 * PI * k * n / N));
                             sum /= static_cast<double>(N);
                             amp_sum += std::abs(sum);
                            }
                            return amp_sum / N; // スカラー
                           }};

  cfg.functions["dftR"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            auto x = collectReals(v, ctx);
                            size_t N = x.size();
                            if (N == 0) throwDomain(ctx.pos);

                            double amp_sum = 0;
                            for (size_t k = 0; k < N; ++k) {
                             std::complex<double> sum = 0;
                             for (size_t n = 0; n < N; ++n)
                              sum += x[n] * std::exp(std::complex<double>(0, -2.0 * PI * k * n / N));
                             amp_sum += std::abs(sum); // 振幅の合計
                            }

                            return amp_sum / N;
                           }};

  // --- Hilbert transform & convolution（簡易スカラ版） ---
  cfg.functions["hilbertR"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                                auto x = collectReals(v, ctx);
                                size_t N = x.size();
                                if (N == 0) throwDomain(ctx.pos);
                                // DFT
                                std::vector<std::complex<double>> X(N);
                                for (size_t k = 0; k < N; ++k) {
                                 std::complex<double> sum = 0;
                                 for (size_t n = 0; n < N; ++n)
                                  sum += x[n] * std::exp(std::complex<double>(0, -2.0 * PI * k * n / N));
                                 X[k] = sum;
                                }
                                // ハーフスペクトルをゼロ化
                                for (size_t k = 1; k < N / 2; ++k)
                                 X[k] *= 2;
                                for (size_t k = N / 2 + 1; k < N; ++k)
                                 X[k] = 0;
                                // IFFT
                                double amp_sum = 0;
                                for (size_t n = 0; n < N; ++n) {
                                 std::complex<double> sum = 0;
                                 for (size_t k = 0; k < N; ++k)
                                  sum += X[k] * std::exp(std::complex<double>(0, 2.0 * PI * k * n / N));
                                 sum /= static_cast<double>(N);
                                 amp_sum += std::abs(sum);
                                }
                                return amp_sum / N;
                               }};
  cfg.functions["convolve"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                if (v.size() < 2) throwDomain(ctx.pos);

                                size_t n = v.size() / 2;
                                if (n == 0 || v.size() % 2 != 0) throwDomain(ctx.pos, "dimension mismatch");

                                std::vector<double> x, y;
                                x.reserve(n);
                                y.reserve(n);

                                for (size_t i = 0; i < n; ++i)
                                 x.push_back(requireReal(v[i], ctx.pos));

                                for (size_t i = 0; i < n; ++i)
                                 y.push_back(requireReal(v[i + n], ctx.pos));

                                long double sum = 0;

                                for (size_t i = 0; i < n; ++i)
                                 for (size_t j = 0; j <= i; ++j)
                                  sum += (long double)x[j] * y[i - j];

                                return (double)sum;
                               }};
 }

 void registerVector(SystemConfig &cfg) {
  // 任意サイズ行列の作成
  cfg.functions["matrix"] = {0, -1, [](auto &v, auto &ctx) -> Value {
                              if (v.empty()) throw CalcError(CalcErrorType::TypeError, "matrix requires at least one row", ctx.pos);

                              // 1引数で既に2Dならそのまま返す
                              if (v.size() == 1 && v[0].isMulti()) {

                               const auto &first = v[0].asMultiRef(ctx.pos);

                               if (!first.empty() && first[0].isMulti()) {
                                return v[0]; // すでに2D
                               }
                              }

                              size_t rows = v.size();
                              size_t cols = v[0].asMultiRef(ctx.pos).size();

                              auto result = std::make_shared<MultiValue>();
                              result->elems_.reserve(rows);

                              for (size_t i = 0; i < rows; ++i) {

                               const auto &row = v[i].asMultiRef(ctx.pos);

                               if (row.size() != cols) throw CalcError(CalcErrorType::DomainError, "all rows must have the same number of columns", ctx.pos);

                               result->elems_.push_back(v[i]);
                              }

                              return Value(result);
                             }};

  // 任意サイズ行列の和
  cfg.functions["madd"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "madd requires two matrix arguments", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            const auto &b = v[1].asMultiRef(ctx.pos);

                            if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "matrix dimension mismatch", ctx.pos);

                            auto result = std::make_shared<MultiValue>();
                            result->elems_.reserve(a.size());

                            for (size_t i = 0; i < a.size(); ++i) {
                             const auto &row_a = a[i].asMultiRef(ctx.pos);
                             const auto &row_b = b[i].asMultiRef(ctx.pos);

                             if (row_a.size() != row_b.size()) { throw CalcError(CalcErrorType::DomainError, "row dimension mismatch", ctx.pos); }

                             auto row_result = std::make_shared<MultiValue>();
                             row_result->elems_.reserve(row_a.size());

                             for (size_t j = 0; j < row_a.size(); ++j) {
                              row_result->elems_.push_back(Value(row_a[j].asScalar(ctx.pos) + row_b[j].asScalar(ctx.pos)));
                             }
                             result->elems_.push_back(Value(row_result));
                            }
                            return Value(result);
                           }};

  // 任意サイズ行列の積
  cfg.functions["mmul"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "mmul requires two matrix arguments", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            const auto &b = v[1].asMultiRef(ctx.pos);

                            if (a[0].asMultiRef(ctx.pos).size() != b.size()) { throw CalcError(CalcErrorType::DomainError, "matrix dimension mismatch for multiplication", ctx.pos); }

                            auto result = std::make_shared<MultiValue>();
                            result->elems_.reserve(a.size());

                            for (size_t i = 0; i < a.size(); ++i) {
                             const auto &row_a = a[i].asMultiRef(ctx.pos);
                             auto row_result = std::make_shared<MultiValue>();
                             row_result->elems_.reserve(b[0].asMultiRef(ctx.pos).size());

                             for (size_t j = 0; j < b[0].asMultiRef(ctx.pos).size(); ++j) {
                              double sum = 0.0;
                              for (size_t k = 0; k < row_a.size(); ++k) {
                               sum += row_a[k].asScalar(ctx.pos) * b[k].asMultiRef(ctx.pos)[j].asScalar(ctx.pos);
                              }
                              row_result->elems_.push_back(Value(sum));
                             }
                             result->elems_.push_back(Value(row_result));
                            }
                            return Value(result);
                           }};

  // 行列の転置
  cfg.functions["mtranspose"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mtranspose requires matrix argument", ctx.pos);
                                  const auto &matrix = v[0].asMultiRef(ctx.pos);

                                  if (matrix.empty()) { throw CalcError(CalcErrorType::DomainError, "cannot transpose empty matrix", ctx.pos); }

                                  auto result = std::make_shared<MultiValue>();
                                  result->elems_.reserve(matrix[0].asMultiRef(ctx.pos).size());

                                  for (size_t j = 0; j < matrix[0].asMultiRef(ctx.pos).size(); ++j) {
                                   auto row_result = std::make_shared<MultiValue>();
                                   row_result->elems_.reserve(matrix.size());
                                   for (size_t i = 0; i < matrix.size(); ++i) {
                                    row_result->elems_.push_back(matrix[i].asMultiRef(ctx.pos)[j]);
                                   }
                                   result->elems_.push_back(Value(row_result));
                                  }
                                  return Value(result);
                                 }};

  // 行列のトレース（対角成分の和）
  cfg.functions["mtrace"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mtrace requires matrix argument", ctx.pos);
                              const auto &matrix = v[0].asMultiRef(ctx.pos);

                              if (matrix.empty() || matrix[0].asMultiRef(ctx.pos).empty()) { throw CalcError(CalcErrorType::DomainError, "cannot calculate trace of empty matrix", ctx.pos); }

                              size_t min_dim = std::min(matrix.size(), matrix[0].asMultiRef(ctx.pos).size());

                              double trace = 0.0;
                              for (size_t i = 0; i < min_dim; ++i) {
                               trace += matrix[i].asMultiRef(ctx.pos)[i].asScalar(ctx.pos);
                              }
                              return Value(trace);
                             }};

  // 行列のランク
  cfg.functions["mrank"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mrank requires matrix argument", ctx.pos);
                             const auto &matrix = v[0].asMultiRef(ctx.pos);

                             if (matrix.empty() || matrix[0].asMultiRef(ctx.pos).empty()) { return Value(0.0); }

                             // シンプルなランク計算（ガウス消去法の簡易版）
                             size_t rows = matrix.size();
                             size_t cols = matrix[0].asMultiRef(ctx.pos).size();

                             // 行列をコピー
                             std::vector<std::vector<double>> temp(rows, std::vector<double>(cols));
                             for (size_t i = 0; i < rows; ++i) {
                              for (size_t j = 0; j < cols; ++j) {
                               temp[i][j] = matrix[i].asMultiRef(ctx.pos)[j].asScalar(ctx.pos);
                              }
                             }

                             size_t rank = 0;
                             for (size_t col = 0; col < cols && rank < rows; ++col) {
                              size_t pivot = rank;
                              for (size_t row = rank + 1; row < rows; ++row) {
                               if (std::abs(temp[row][col]) > std::abs(temp[pivot][col])) { pivot = row; }
                              }

                              if (std::abs(temp[pivot][col]) < 1e-10) continue;

                              // 行の交換
                              if (pivot != rank) { std::swap(temp[rank], temp[pivot]); }

                              // 前進消去
                              for (size_t row = rank + 1; row < rows; ++row) {
                               double factor = temp[row][col] / temp[rank][col];
                               for (size_t k = col; k < cols; ++k) {
                                temp[row][k] -= factor * temp[rank][k];
                               }
                              }
                              rank++;
                             }
                             return Value(static_cast<double>(rank));
                            }};

  // 行列の決定因子（行列式）
  cfg.functions["mdet"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mdet requires matrix argument", ctx.pos);
                            const auto &matrix = v[0].asMultiRef(ctx.pos);

                            if (matrix.empty() || matrix[0].asMultiRef(ctx.pos).empty()) { return Value(0.0); }

                            size_t rows = matrix.size();
                            size_t cols = matrix[0].asMultiRef(ctx.pos).size();

                            if (rows != cols) { throw CalcError(CalcErrorType::DomainError, "determinant requires square matrix", ctx.pos); }

                            // ガウス消去法で行列式を計算
                            std::vector<std::vector<double>> temp(rows, std::vector<double>(cols));
                            for (size_t i = 0; i < rows; ++i) {
                             for (size_t j = 0; j < cols; ++j) {
                              temp[i][j] = matrix[i].asMultiRef(ctx.pos)[j].asScalar(ctx.pos);
                             }
                            }

                            double det = 1.0;
                            for (size_t i = 0; i < rows; ++i) {
                             size_t pivot = i;
                             for (size_t j = i + 1; j < rows; ++j) {
                              if (std::abs(temp[j][i]) > std::abs(temp[pivot][i])) { pivot = j; }
                             }

                             if (std::abs(temp[pivot][i]) < 1e-10) { return Value(0.0); }

                             if (pivot != i) {
                              std::swap(temp[i], temp[pivot]);
                              det *= -1.0;
                             }

                             det *= temp[i][i];

                             for (size_t j = i + 1; j < rows; ++j) {
                              double factor = temp[j][i] / temp[i][i];
                              for (size_t k = i; k < cols; ++k) {
                               temp[j][k] -= factor * temp[i][k];
                              }
                             }
                            }
                            return Value(det);
                           }};

  // 単位行列の作成
  cfg.functions["identity"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                int size = static_cast<int>(v[0].asScalar(ctx.pos));
                                if (size <= 0) throw CalcError(CalcErrorType::DomainError, "matrix size must be positive", ctx.pos);

                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(size);

                                for (int i = 0; i < size; ++i) {
                                 auto row = std::make_shared<MultiValue>();
                                 row->elems_.reserve(size);
                                 for (int j = 0; j < size; ++j) {
                                  row->elems_.push_back(Value(i == j ? 1.0 : 0.0));
                                 }
                                 result->elems_.push_back(Value(row));
                                }
                                return Value(result);
                               }};

  // ゼロ行列の作成
  cfg.functions["zeros"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             int rows = static_cast<int>(v[0].asScalar(ctx.pos));
                             int cols = static_cast<int>(v[1].asScalar(ctx.pos));
                             if (rows <= 0 || cols <= 0) throw CalcError(CalcErrorType::DomainError, "matrix dimensions must be positive", ctx.pos);

                             auto result = std::make_shared<MultiValue>();
                             result->elems_.reserve(rows);

                             for (int i = 0; i < rows; ++i) {
                              auto row = std::make_shared<MultiValue>();
                              row->elems_.reserve(cols);
                              for (int j = 0; j < cols; ++j) {
                               row->elems_.push_back(Value(0.0));
                              }
                              result->elems_.push_back(Value(row));
                             }
                             return Value(result);
                            }};

  // 行列の要素アクセス
  cfg.functions["mget"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mget requires matrix argument", ctx.pos);
                            const auto &matrix = v[0].asMultiRef(ctx.pos);
                            int row = static_cast<int>(v[1].asScalar(ctx.pos));
                            int col = static_cast<int>(v[2].asScalar(ctx.pos));

                            if (row < 0 || row >= static_cast<int>(matrix.size())) { throw CalcError(CalcErrorType::DomainError, "row index out of bounds", ctx.pos); }
                            if (col < 0 || col >= static_cast<int>(matrix[0].asMultiRef(ctx.pos).size())) { throw CalcError(CalcErrorType::DomainError, "column index out of bounds", ctx.pos); }

                            return Value(matrix[row].asMultiRef(ctx.pos)[col].asScalar(ctx.pos));
                           }};

  // 行列の行数取得
  cfg.functions["mrows"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mrows requires matrix argument", ctx.pos);
                             const auto &matrix = v[0].asMultiRef(ctx.pos);
                             return Value(static_cast<double>(matrix.size()));
                            }};

  // 行列の列数取得
  cfg.functions["mcols"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mcols requires matrix argument", ctx.pos);
                             const auto &matrix = v[0].asMultiRef(ctx.pos);
                             if (matrix.empty()) return Value(0.0);
                             return Value(static_cast<double>(matrix[0].asMultiRef(ctx.pos).size()));
                            }};

  // 行列の対角成分取得
  cfg.functions["mdiag"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mdiag requires matrix argument", ctx.pos);
                             const auto &matrix = v[0].asMultiRef(ctx.pos);

                             if (matrix.empty() || matrix[0].asMultiRef(ctx.pos).empty()) {
                              auto result = std::make_shared<MultiValue>();
                              return Value(result);
                             }

                             size_t min_dim = std::min(matrix.size(), matrix[0].asMultiRef(ctx.pos).size());

                             auto result = std::make_shared<MultiValue>();
                             result->elems_.reserve(min_dim);

                             for (size_t i = 0; i < min_dim; ++i) {
                              result->elems_.push_back(matrix[i].asMultiRef(ctx.pos)[i]);
                             }
                             return Value(result);
                            }};

  // 行列のLU分解
  cfg.functions["mlu"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mlu requires matrix argument", ctx.pos);

                           const auto &M = v[0].asMultiRef(ctx.pos);

                           size_t n = M.size();
                           if (n == 0) throw CalcError(CalcErrorType::DomainError, "empty matrix", ctx.pos);

                           size_t m = M[0].asMultiRef(ctx.pos).size();

                           if (n != m) throw CalcError(CalcErrorType::DomainError, "LU requires square matrix", ctx.pos);

                           std::vector<std::vector<double>> A(n, std::vector<double>(n));
                           for (size_t i = 0; i < n; ++i) {
                            const auto &row = M[i].asMultiRef(ctx.pos);
                            if (row.size() != n) throw CalcError(CalcErrorType::DomainError, "irregular matrix", ctx.pos);

                            for (size_t j = 0; j < n; ++j)
                             A[i][j] = row[j].asScalar(ctx.pos);
                           }

                           // L と U
                           std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
                           std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));

                           for (size_t i = 0; i < n; ++i) {
                            // --- U 計算 ---
                            for (size_t j = i; j < n; ++j) {
                             double sum = 0.0;
                             for (size_t k = 0; k < i; ++k)
                              sum += L[i][k] * U[k][j];

                             U[i][j] = A[i][j] - sum;
                            }

                            if (std::fabs(U[i][i]) < 1e-14) throw CalcError(CalcErrorType::DomainError, "singular matrix in LU", ctx.pos);

                            // --- L 計算 ---
                            L[i][i] = 1.0;

                            for (size_t j = i + 1; j < n; ++j) {
                             double sum = 0.0;
                             for (size_t k = 0; k < i; ++k)
                              sum += L[j][k] * U[k][i];

                             L[j][i] = (A[j][i] - sum) / U[i][i];
                            }
                           }

                           // MultiValue化
                           auto l_mv = std::make_shared<MultiValue>();
                           auto u_mv = std::make_shared<MultiValue>();

                           for (size_t i = 0; i < n; ++i) {
                            auto l_row = std::make_shared<MultiValue>();
                            auto u_row = std::make_shared<MultiValue>();

                            for (size_t j = 0; j < n; ++j) {
                             l_row->elems_.emplace_back(L[i][j]);
                             u_row->elems_.emplace_back(U[i][j]);
                            }

                            l_mv->elems_.emplace_back(Value(l_row));
                            u_mv->elems_.emplace_back(Value(u_row));
                           }

                           auto result = std::make_shared<MultiValue>();
                           result->elems_.emplace_back(Value(l_mv));
                           result->elems_.emplace_back(Value(u_mv));

                           return Value(result);
                          }};

  // 行列の固有値
  cfg.functions["meigenvals"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  std::vector<std::vector<double>> A = toMatrix(v[0], ctx);
                                  if (A.size() == 1) {
                                   auto result = std::make_shared<MultiValue>();
                                   result->elems_.emplace_back(A[0][0]);
                                   return Value(result);
                                  }

                                  auto eig = eigenvalues(A);

                                  auto result = std::make_shared<MultiValue>();
                                  result->elems_.reserve(eig.size());

                                  for (double lamba : eig)
                                   result->elems_.emplace_back(lamba);
                                  std::sort(eig.begin(), eig.end(), std::greater<double>());
                                  return Value(result);
                                 }};

  // 行列の固有ベクトル
  cfg.functions["meigenvecs"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "meigenvecs requires matrix argument", ctx.pos);

                                  const auto &matrix = v[0].asMultiRef(ctx.pos);

                                  // 行列の次元確認
                                  if (matrix.empty()) throw CalcError(CalcErrorType::InvalidOperation, "std::vector<std::vector<double>> is empty", ctx.pos);

                                  size_t n = matrix.size(); // 行列の次元

                                  // 正方行列でない場合はエラー
                                  for (size_t i = 0; i < n; ++i) {
                                   if (matrix[i].size() != n) throw CalcError(CalcErrorType::InvalidOperation, "std::vector<std::vector<double>> must be square", ctx.pos);
                                  }

                                  // 行列のコピーを作成（元の行列を変更しない）
                                  std::vector<std::vector<double>> A(n, std::vector<double>(n));
                                  for (size_t i = 0; i < n; ++i) {
                                   const auto &row = matrix[i].asMultiRef(ctx.pos);
                                   for (size_t j = 0; j < n; ++j) {
                                    A[i][j] = row[j].asScalar(ctx.pos);
                                   }
                                  }

                                  // 固有値を求める（QR反復法）
                                  std::vector<double> eigenvals(n);
                                  std::vector<std::vector<double>> V(n, std::vector<double>(n, 0.0));

                                  // 初期化
                                  for (size_t i = 0; i < n; ++i) {
                                   V[i][i] = 1.0;
                                  }

                                  // QR反復法による固有値計算
                                  const int max_iterations = 1000;
                                  const double epsilon = 1e-12;

                                  for (int iter = 0; iter < max_iterations; ++iter) {
                                   bool converged = true;

                                   // QR分解を計算
                                   std::vector<std::vector<double>> Q(n, std::vector<double>(n, 0.0));
                                   std::vector<std::vector<double>> R = A;

                                   // 単位行列で初期化
                                   for (size_t i = 0; i < n; ++i) {
                                    Q[i][i] = 1.0;
                                   }

                                   // ガウス・ホルツ変換でQR分解
                                   for (size_t k = 0; k < n - 1; ++k) {
                                    // ベクトルv = R[k:n][k] を取得
                                    std::vector<double> v(n - k);
                                    for (size_t i = 0; i < n - k; ++i) {
                                     v[i] = R[k + i][k];
                                    }

                                    // ノルムを計算
                                    double norm = 0.0;
                                    for (double val : v) {
                                     norm += val * val;
                                    }
                                    norm = std::sqrt(norm);

                                    if (norm < epsilon) continue;

                                    // 通常の符号を持つ値を選択
                                    double sign = (v[0] >= 0) ? 1.0 : -1.0;
                                    double alpha = -sign * norm;

                                    // w = v + alpha * e_1
                                    std::vector<double> w(v);
                                    w[0] += alpha;

                                    // ノルムを正規化
                                    double w_norm = 0.0;
                                    for (double val : w) {
                                     w_norm += val * val;
                                    }
                                    w_norm = std::sqrt(w_norm);

                                    if (w_norm < epsilon) continue;

                                    // ハウスホルダー変換の計算
                                    std::vector<double> u(w);
                                    for (double &val : u) {
                                     val /= w_norm;
                                    }

                                    // H = I - 2 * u * u^T を計算
                                    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0.0));
                                    for (size_t i = 0; i < n; ++i) {
                                     H[i][i] = 1.0;
                                    }
                                    for (size_t i = 0; i < n - k; ++i) {
                                     for (size_t j = 0; j < n - k; ++j) {
                                      H[k + i][k + j] -= 2.0 * u[i] * u[j];
                                     }
                                    }

                                    // 行列の積を計算
                                    std::vector<std::vector<double>> temp_R(n, std::vector<double>(n));
                                    for (size_t i = 0; i < n; ++i) {
                                     for (size_t j = 0; j < n; ++j) {
                                      double sum = 0.0;
                                      for (size_t k_idx = 0; k_idx < n; ++k_idx) {
                                       sum += H[i][k_idx] * R[k_idx][j];
                                      }
                                      temp_R[i][j] = sum;
                                     }
                                    }
                                    R = temp_R;

                                    // Q = Q * H^T
                                    std::vector<std::vector<double>> temp_Q(n, std::vector<double>(n));
                                    for (size_t i = 0; i < n; ++i) {
                                     for (size_t j = 0; j < n; ++j) {
                                      double sum = 0.0;
                                      for (size_t k_idx = 0; k_idx < n; ++k_idx) {
                                       sum += Q[i][k_idx] * H[j][k_idx];
                                      }
                                      temp_Q[i][j] = sum;
                                     }
                                    }
                                    Q = temp_Q;
                                   }

                                   // 収束判定
                                   for (size_t i = 0; i < n - 1; ++i) {
                                    if (std::abs(R[i][i + 1]) > epsilon) {
                                     converged = false;
                                     break;
                                    }
                                   }

                                   if (converged) {
                                    // 固有値を対角成分から取得
                                    for (size_t i = 0; i < n; ++i) {
                                     eigenvals[i] = R[i][i];
                                    }
                                    break;
                                   }
                                  }

                                  // 固有ベクトルの計算（固有値ごとに求める）
                                  auto result = std::make_shared<MultiValue>();

                                  for (size_t k = 0; k < n; ++k) {
                                   double eigenval = eigenvals[k];
                                   std::vector<std::vector<double>> I_minus_lambda_A(n, std::vector<double>(n, 0.0));

                                   // (A - λI) を計算
                                   for (size_t i = 0; i < n; ++i) {
                                    for (size_t j = 0; j < n; ++j) {
                                     if (i == j) {
                                      I_minus_lambda_A[i][j] = A[i][j] - eigenval;
                                     } else {
                                      I_minus_lambda_A[i][j] = A[i][j];
                                     }
                                    }
                                   }

                                   // 連立方程式 (A - λI) x = 0 を解く（簡単な方法）
                                   std::vector<double> eigenvector(n, 0.0);
                                   bool found = false;

                                   // 単純なガウス・ザイデル法で解を求める
                                   const int gauss_iter = 100;
                                   for (int iter = 0; iter < gauss_iter; ++iter) {
                                    std::vector<double> new_eigenvector(n, 0.0);

                                    for (size_t i = 0; i < n; ++i) {
                                     double sum = 0.0;
                                     for (size_t j = 0; j < n; ++j) {
                                      if (i != j) { sum += I_minus_lambda_A[i][j] * eigenvector[j]; }
                                     }
                                     if (std::abs(I_minus_lambda_A[i][i]) > epsilon) { new_eigenvector[i] = -sum / I_minus_lambda_A[i][i]; }
                                    }

                                    // 正規化
                                    double norm = 0.0;
                                    for (double val : new_eigenvector) {
                                     norm += val * val;
                                    }
                                    norm = std::sqrt(norm);

                                    if (norm > epsilon) {
                                     for (double &val : new_eigenvector) {
                                      val /= norm;
                                     }
                                     eigenvector = new_eigenvector;
                                    } else {
                                     // ベクトルがゼロになる場合は単位ベクトルに設定
                                     eigenvector = std::vector<double>(n, 0.0);
                                     eigenvector[0] = 1.0;
                                     break;
                                    }
                                   }

                                   // 固有ベクトルを結果に追加
                                   auto vec = std::make_shared<MultiValue>();
                                   vec->elems_.reserve(n);
                                   for (double val : eigenvector) {
                                    vec->elems_.push_back(Value(val));
                                   }
                                   result->elems_.push_back(Value(vec));
                                  }

                                  return Value(result);
                                 }};

  // 行列の逆行列（ガウス・ジョルダン法）
  cfg.functions["minverse"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "minverse requires matrix argument", ctx.pos);
                                const auto &matrix = v[0].asMultiRef(ctx.pos);

                                size_t rows = matrix.size();
                                size_t cols = matrix[0].asMultiRef(ctx.pos).size();

                                if (rows != cols) { throw CalcError(CalcErrorType::DomainError, "inverse requires square matrix", ctx.pos); }

                                // 単位行列を作成
                                std::vector<std::vector<double>> augmented(rows, std::vector<double>(2 * cols));

                                // 元の行列をaugmented行列にコピー
                                for (size_t i = 0; i < rows; ++i) {
                                 for (size_t j = 0; j < cols; ++j) {
                                  augmented[i][j] = matrix[i].asMultiRef(ctx.pos)[j].asScalar(ctx.pos);
                                 }
                                 augmented[i][i + cols] = 1.0; // 単位行列の部分
                                }

                                // ガウス・ジョルダン法
                                for (size_t i = 0; i < rows; ++i) {
                                 // ピボット選択
                                 size_t pivot_row = i;
                                 for (size_t j = i + 1; j < rows; ++j) {
                                  if (std::abs(augmented[j][i]) > std::abs(augmented[pivot_row][i])) { pivot_row = j; }
                                 }

                                 // 行の交換
                                 if (pivot_row != i) { std::swap(augmented[i], augmented[pivot_row]); }

                                 // 原点化
                                 double pivot = augmented[i][i];
                                 if (std::abs(pivot) < 1e-10) { throw CalcError(CalcErrorType::DomainError, "matrix is singular, cannot invert", ctx.pos); }

                                 for (size_t j = i; j < 2 * cols; ++j) {
                                  augmented[i][j] /= pivot;
                                 }

                                 // 前進消去
                                 for (size_t j = 0; j < rows; ++j) {
                                  if (j != i) {
                                   double factor = augmented[j][i];
                                   for (size_t k = i; k < 2 * cols; ++k) {
                                    augmented[j][k] -= factor * augmented[i][k];
                                   }
                                  }
                                 }
                                }

                                // 結果を返す
                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(rows);

                                for (size_t i = 0; i < rows; ++i) {
                                 auto row = std::make_shared<MultiValue>();
                                 row->elems_.reserve(cols);
                                 for (size_t j = 0; j < cols; ++j) {
                                  row->elems_.push_back(Value(augmented[i][j + cols]));
                                 }
                                 result->elems_.push_back(Value(row));
                                }
                                return Value(result);
                               }};

  // 行列のQR分解
  cfg.functions["mqr"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mqr requires matrix argument", ctx.pos);

                           const auto &matrix = v[0].asMultiRef(ctx.pos);

                           // 行列の次元確認
                           if (matrix.empty()) throw CalcError(CalcErrorType::InvalidOperation, "std::vector<std::vector<double>> is empty", ctx.pos);

                           size_t m = matrix.size();    // 行数
                           size_t n = matrix[0].size(); // 列数

                           // 行列のコピーを作成（元の行列を変更しない）
                           std::vector<std::vector<double>> A(m, std::vector<double>(n));
                           for (size_t i = 0; i < m; ++i) {
                            if (matrix[i].size() != n) throw CalcError(CalcErrorType::InvalidOperation, "std::vector<std::vector<double>> has inconsistent row lengths", ctx.pos);
                            const auto &row = matrix[i].asMultiRef(ctx.pos);
                            for (size_t j = 0; j < n; ++j) {
                             A[i][j] = row[j].asScalar(ctx.pos);
                            }
                           }

                           // QR分解の実装（ハウスホルダー変換）
                           std::vector<std::vector<double>> Q(m, std::vector<double>(m, 0.0));
                           std::vector<std::vector<double>> R = A;

                           // 単位行列で初期化
                           for (size_t i = 0; i < m; ++i) {
                            Q[i][i] = 1.0;
                           }

                           // ハウスホルダー変換を適用
                           for (size_t k = 0; k < std::min(m, n); ++k) {
                            // ベクトルv = R[k:m][k] を取得
                            std::vector<double> v(m - k);
                            for (size_t i = 0; i < m - k; ++i) {
                             v[i] = R[k + i][k];
                            }

                            // ノルムを計算
                            double norm = 0.0;
                            for (double val : v) {
                             norm += val * val;
                            }
                            norm = std::sqrt(norm);

                            // 通常の符号を持つ値を選択
                            double sign = (v[0] >= 0) ? 1.0 : -1.0;
                            /*double alpha = -sign * norm;*/
                            double alpha = (v[0] >= 0) ? -norm : norm;

                            // w = v + alpha * e_1
                            std::vector<double> w(v);
                            w[0] += alpha;

                            // ノルムを正規化
                            double w_norm = 0.0;
                            for (double val : w) {
                             w_norm += val * val;
                            }
                            w_norm = std::sqrt(w_norm);

                            if (w_norm < 1e-12) continue; // 0ベクトルならスキップ

                            // ハウスホルダー変換の計算
                            std::vector<double> u(w);
                            for (double &val : u) {
                             val /= w_norm;
                            }

                            // H = I - 2 * u * u^T を計算
                            std::vector<std::vector<double>> H(m, std::vector<double>(m, 0.0));
                            for (size_t i = 0; i < m; ++i) {
                             H[i][i] = 1.0;
                            }
                            for (size_t i = 0; i < m - k; ++i) {
                             for (size_t j = 0; j < m - k; ++j) {
                              H[k + i][k + j] -= 2.0 * u[i] * u[j];
                             }
                            }

                            // 行列の積を計算
                            // R = H * R
                            std::vector<std::vector<double>> temp_R(m, std::vector<double>(n));
                            for (size_t i = 0; i < m; ++i) {
                             for (size_t j = 0; j < n; ++j) {
                              double sum = 0.0;
                              for (size_t k_idx = 0; k_idx < m; ++k_idx) {
                               sum += H[i][k_idx] * R[k_idx][j];
                              }
                              temp_R[i][j] = sum;
                             }
                            }
                            R = temp_R;

                            // Q = Q * H^T
                            std::vector<std::vector<double>> temp_Q(m, std::vector<double>(m));
                            for (size_t i = 0; i < m; ++i) {
                             for (size_t j = 0; j < m; ++j) {
                              double sum = 0.0;
                              for (size_t k_idx = 0; k_idx < m; ++k_idx) {
                               sum += Q[i][k_idx] * H[j][k_idx];
                              }
                              temp_Q[i][j] = sum;
                             }
                            }
                            Q = temp_Q;
                           }

                           // 結果の構築
                           auto result = std::make_shared<MultiValue>();

                           // Q行列の構築
                           auto q_matrix = std::make_shared<MultiValue>();
                           q_matrix->elems_.reserve(m);
                           for (size_t i = 0; i < m; ++i) {
                            auto row = std::make_shared<MultiValue>();
                            row->elems_.reserve(m);
                            for (size_t j = 0; j < m; ++j) {
                             row->elems_.push_back(Value(Q[i][j]));
                            }
                            q_matrix->elems_.push_back(Value(row));
                           }
                           result->elems_.push_back(Value(q_matrix));

                           // R行列の構築
                           auto r_matrix = std::make_shared<MultiValue>();
                           r_matrix->elems_.reserve(m);
                           for (size_t i = 0; i < m; ++i) {
                            auto row = std::make_shared<MultiValue>();
                            row->elems_.reserve(n);
                            for (size_t j = 0; j < n; ++j) {
                             row->elems_.push_back(Value(R[i][j]));
                            }
                            r_matrix->elems_.push_back(Value(row));
                           }
                           result->elems_.push_back(Value(r_matrix));

                           return Value(result);
                          }};

  // 行列の特異値分解（SVD）
  cfg.functions["msvd"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "msvd requires matrix argument", ctx.pos);

                            const auto &matrix = v[0].asMultiRef(ctx.pos);
                            if (matrix.empty()) throw CalcError(CalcErrorType::InvalidOperation, "matrix is empty", ctx.pos);

                            size_t m = matrix.size();
                            size_t n = matrix[0].asMultiRef(ctx.pos).size();

                            if (m < n) throw CalcError(CalcErrorType::InvalidOperation, "msvd currently requires m >= n", ctx.pos);

                            // --- A を double に変換 ---
                            std::vector<std::vector<double>> A(m, std::vector<double>(n));
                            for (size_t i = 0; i < m; ++i) {
                             const auto &row = matrix[i].asMultiRef(ctx.pos);
                             if (row.size() != n) throw CalcError(CalcErrorType::InvalidOperation, "inconsistent row size", ctx.pos);
                             for (size_t j = 0; j < n; ++j)
                              A[i][j] = row[j].asScalar(ctx.pos);
                            }

                            // --- AtA 計算 ---
                            std::vector<std::vector<double>> AtA(n, std::vector<double>(n, 0.0));
                            for (size_t i = 0; i < n; ++i)
                             for (size_t j = 0; j < n; ++j)
                              for (size_t k = 0; k < m; ++k)
                               AtA[i][j] += A[k][i] * A[k][j];

                            // --- Jacobi 固有値分解（対称行列）---
                            std::vector<std::vector<double>> V(n, std::vector<double>(n, 0.0));
                            for (size_t i = 0; i < n; ++i)
                             V[i][i] = 1.0;

                            const int max_iter = 100;
                            const double eps = 1e-12;

                            for (int iter = 0; iter < max_iter; ++iter) {
                             bool converged = true;
                             for (size_t p = 0; p < n - 1; ++p) {
                              for (size_t q = p + 1; q < n; ++q) {

                               if (std::abs(AtA[p][q]) < eps) continue;

                               converged = false;

                               double theta = 0.5 * std::atan2(2.0 * AtA[p][q], AtA[q][q] - AtA[p][p]);

                               double c = std::cos(theta);
                               double s = std::sin(theta);

                               for (size_t i = 0; i < n; ++i) {
                                double ip = AtA[i][p];
                                double iq = AtA[i][q];
                                AtA[i][p] = c * ip - s * iq;
                                AtA[i][q] = s * ip + c * iq;
                               }

                               for (size_t i = 0; i < n; ++i) {
                                double pi = AtA[p][i];
                                double qi = AtA[q][i];
                                AtA[p][i] = c * pi - s * qi;
                                AtA[q][i] = s * pi + c * qi;
                               }

                               for (size_t i = 0; i < n; ++i) {
                                double vip = V[i][p];
                                double viq = V[i][q];
                                V[i][p] = c * vip - s * viq;
                                V[i][q] = s * vip + c * viq;
                               }
                              }
                             }
                             if (converged) break;
                            }

                            // --- 特異値 ---
                            std::vector<double> sigma(n);
                            for (size_t i = 0; i < n; ++i)
                             sigma[i] = std::sqrt(std::max(0.0, AtA[i][i]));

                            // --- 降順ソート ---
                            for (size_t i = 0; i < n - 1; ++i) {
                             for (size_t j = i + 1; j < n; ++j) {
                              if (sigma[i] < sigma[j]) {
                               std::swap(sigma[i], sigma[j]);
                               for (size_t k = 0; k < n; ++k)
                                std::swap(V[k][i], V[k][j]);
                              }
                             }
                            }

                            // --- U = A V Σ^-1 ---
                            std::vector<std::vector<double>> U(m, std::vector<double>(n, 0.0));

                            for (size_t i = 0; i < n; ++i) {
                             if (sigma[i] < eps) continue;

                             for (size_t r = 0; r < m; ++r)
                              for (size_t c = 0; c < n; ++c)
                               U[r][i] += A[r][c] * V[c][i];

                             for (size_t r = 0; r < m; ++r)
                              U[r][i] /= sigma[i];
                            }

                            // --- MultiValue構築 ---
                            auto result = std::make_shared<MultiValue>();

                            // U
                            auto u_mat = std::make_shared<MultiValue>();
                            for (size_t i = 0; i < m; ++i) {
                             auto row = std::make_shared<MultiValue>();
                             for (size_t j = 0; j < n; ++j)
                              row->elems_.push_back(Value(U[i][j]));
                             u_mat->elems_.push_back(Value(row));
                            }
                            result->elems_.push_back(Value(u_mat));

                            // Σ
                            auto s_vec = std::make_shared<MultiValue>();
                            for (double s : sigma)
                             s_vec->elems_.push_back(Value(s));
                            result->elems_.push_back(Value(s_vec));

                            // V^T
                            auto vt_mat = std::make_shared<MultiValue>();
                            for (size_t i = 0; i < n; ++i) {
                             auto row = std::make_shared<MultiValue>();
                             for (size_t j = 0; j < n; ++j)
                              row->elems_.push_back(Value(V[j][i]));
                             vt_mat->elems_.push_back(Value(row));
                            }
                            result->elems_.push_back(Value(vt_mat));

                            return Value(result);
                           }};

  // 行列の軌道（行列式の絶対値）
  cfg.functions["mcond"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             auto A = toMatrix(v[0], ctx);

                             auto S = singular_values(A);

                             if (S.empty()) throw CalcError(CalcErrorType::DomainError, "matrix is empty", ctx.pos);

                             double smax = *std::max_element(S.begin(), S.end());
                             double smin = *std::min_element(S.begin(), S.end());

                             if (smin == 0.0) return Value(std::numeric_limits<double>::infinity());

                             return Value(smax / smin);
                            }};

  // 行列の条件数
  cfg.functions["mcondition"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mcondition requires matrix argument", ctx.pos);

                                  const auto &matrix = v[0].asMultiRef(ctx.pos);

                                  // 行列の次元確認
                                  if (matrix.empty()) throw CalcError(CalcErrorType::InvalidOperation, "std::vector<std::vector<double>> is empty", ctx.pos);

                                  size_t m = matrix.size();    // 行数
                                  size_t n = matrix[0].size(); // 列数

                                  // 行列のコピーを作成（元の行列を変更しない）
                                  std::vector<std::vector<double>> A(m, std::vector<double>(n));
                                  for (size_t i = 0; i < m; ++i) {
                                   if (matrix[i].size() != n) throw CalcError(CalcErrorType::InvalidOperation, "std::vector<std::vector<double>> has inconsistent row lengths", ctx.pos);
                                   const auto &row = matrix[i].asMultiRef(ctx.pos);
                                   for (size_t j = 0; j < n; ++j) {
                                    A[i][j] = row[j].asScalar(ctx.pos);
                                   }
                                  }

                                  // 特異値を求める（SVDの実装）
                                  std::vector<double> singular_values;

                                  if (m >= n) {
                                   // A^T * A を計算
                                   std::vector<std::vector<double>> ATA(n, std::vector<double>(n, 0.0));
                                   for (size_t i = 0; i < n; ++i) {
                                    for (size_t j = 0; j < n; ++j) {
                                     double sum = 0.0;
                                     for (size_t k = 0; k < m; ++k) {
                                      sum += A[k][i] * A[k][j];
                                     }
                                     ATA[i][j] = sum;
                                    }
                                   }
                                   // 固有値を求める（QR反復法）
                                   singular_values = compute_eigenvalues(ATA, ctx);
                                  } else {
                                   // A * A^T を計算
                                   std::vector<std::vector<double>> AAT(m, std::vector<double>(m, 0.0));
                                   for (size_t i = 0; i < m; ++i) {
                                    for (size_t j = 0; j < m; ++j) {
                                     double sum = 0.0;
                                     for (size_t k = 0; k < n; ++k) {
                                      sum += A[i][k] * A[j][k];
                                     }
                                     AAT[i][j] = sum;
                                    }
                                   }
                                   // 固有値を求める（QR反復法）
                                   singular_values = compute_eigenvalues(AAT, ctx);
                                  }

                                  // 特異値は固有値の平方根
                                  std::vector<double> svd_values;
                                  for (double val : singular_values) {
                                   if (val < 0) val = 0; // 数値誤差で負になる場合がある
                                   svd_values.push_back(std::sqrt(val));
                                  }

                                  // 条件数を計算（最大特異値 / 最小特異値）
                                  if (svd_values.empty()) { return Value(std::numeric_limits<double>::infinity()); }

                                  double max_sv = *std::max_element(svd_values.begin(), svd_values.end());
                                  double min_sv = *std::min_element(svd_values.begin(), svd_values.end());

                                  // 最小特異値がゼロに近い場合（数値誤差）
                                  if (min_sv < 1e-12) { return Value(std::numeric_limits<double>::infinity()); }

                                  return Value(max_sv / min_sv);
                                 }};

  // 行列の最小二乗法解(LU分解だよ，SVDやQRは無理ゲ)
  // cfg.functions["mlsqr"] = {2, 2, [](auto &v, auto &ctx) -> Value {
  //                           if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "mlsqr requires matrix and vector arguments", ctx.pos);

  //                           const auto &matrix_A = v[0].asMultiRef(ctx.pos);
  //                           const auto &vector_b = v[1].asMultiRef(ctx.pos);

  //                           if (matrix_A.empty()) throw CalcError(CalcErrorType::InvalidOperation, "matrix A is empty", ctx.pos);

  //                           size_t m = matrix_A.size();
  //                           size_t n = matrix_A[0].asMultiRef(ctx.pos).size();

  //                           if (vector_b.size() != m) throw CalcError(CalcErrorType::InvalidOperation, "A and b dimension mismatch", ctx.pos);

  //                           // --- A, b を double に変換 ---
  //                           std::vector<std::vector<double>> A(m, std::vector<double>(n));
  //                           for (size_t i = 0; i < m; ++i) {
  //                            const auto &row = matrix_A[i].asMultiRef(ctx.pos);
  //                            for (size_t j = 0; j < n; ++j)
  //                             A[i][j] = row[j].asScalar(ctx.pos);
  //                           }

  //                           std::vector<double> b(m);
  //                           for (size_t i = 0; i < m; ++i)
  //                            b[i] = vector_b[i].asScalar(ctx.pos);

  //                           // --- A^T A ---
  //                           std::vector<std::vector<double>> ATA(n, std::vector<double>(n, 0.0));
  //                           for (size_t i = 0; i < n; ++i)
  //                            for (size_t j = 0; j < n; ++j)
  //                             for (size_t k = 0; k < m; ++k)
  //                              ATA[i][j] += A[k][i] * A[k][j];

  //                           // --- A^T b ---
  //                           std::vector<double> ATb(n, 0.0);
  //                           for (size_t i = 0; i < n; ++i)
  //                            for (size_t k = 0; k < m; ++k)
  //                             ATb[i] += A[k][i] * b[k];

  //                           // --- LU分解 ---
  //                           std::vector<std::vector<double>> L, U;
  //                           std::vector<size_t> P;
  //                           lu_decompose(ATA, L, U, P);

  //                           // --- 前進代入 Ly = ATb ---
  //                           std::vector<double> y(n);
  //                           for (size_t i = 0; i < n; ++i) {
  //                            y[i] = ATb[P[i]];
  //                            for (size_t j = 0; j < i; ++j)
  //                             y[i] -= L[i][j] * y[j];
  //                           }

  //                           // --- 後退代入 Ux = y ---
  //                           std::vector<double> x(n);
  //                           for (int i = (int)n - 1; i >= 0; --i) {
  //                            x[i] = y[i];
  //                            for (size_t j = i + 1; j < n; ++j)
  //                             x[i] -= U[i][j] * x[j];

  //                            if (std::abs(U[i][i]) < 1e-12) throw CalcError(CalcErrorType::InvalidOperation, "matrix is singular", ctx.pos);

  //                            x[i] /= U[i][i];
  //                           }

  //                           auto result = std::make_shared<MultiValue>();
  //                           result->elems_.reserve(n);
  //                           for (double val : x)
  //                            result->elems_.push_back(Value(val));

  //                           return Value(result);
  //                          }};

  cfg.functions["mlsqr"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "mlsqr requires matrix and vector arguments", ctx.pos);

                             const auto &matrix_A = v[0].asMultiRef(ctx.pos);
                             const auto &vector_b = v[1].asMultiRef(ctx.pos);

                             if (matrix_A.empty()) throw CalcError(CalcErrorType::InvalidOperation, "matrix A is empty", ctx.pos);

                             size_t m = matrix_A.size();
                             size_t n = matrix_A[0].asMultiRef(ctx.pos).size();

                             if (vector_b.size() != m) throw CalcError(CalcErrorType::InvalidOperation, "A and b dimension mismatch", ctx.pos);

                             // --- A, b を double に変換 ---
                             std::vector<std::vector<double>> A(m, std::vector<double>(n));
                             for (size_t i = 0; i < m; ++i) {
                              const auto &row = matrix_A[i].asMultiRef(ctx.pos);
                              for (size_t j = 0; j < n; ++j)
                               A[i][j] = row[j].asScalar(ctx.pos);
                             }

                             std::vector<double> b(m);
                             for (size_t i = 0; i < m; ++i)
                              b[i] = vector_b[i].asScalar(ctx.pos);

                             // --- QR分解 ---
                             if (m < n) throw CalcError(CalcErrorType::InvalidOperation, "mlsqr requires m >= n (overdetermined system)", ctx.pos);
                             std::vector<std::vector<double>> Q, R;
                             qr_decompose(A, Q, R);

                             // --- y = Q^T b ---
                             std::vector<double> y(n, 0.0);
                             for (size_t i = 0; i < n; ++i)
                              for (size_t k = 0; k < m; ++k)
                               y[i] += Q[k][i] * b[k];

                             // --- 後退代入 R x = y ---
                             std::vector<double> x(n);
                             for (int i = (int)n - 1; i >= 0; --i) {

                              if (std::abs(R[i][i]) < 1e-12) throw CalcError(CalcErrorType::InvalidOperation, "matrix is rank deficient", ctx.pos);

                              x[i] = y[i];
                              for (size_t j = i + 1; j < n; ++j)
                               x[i] -= R[i][j] * x[j];

                              x[i] /= R[i][i];
                             }

                             // --- 結果をMultiValueへ ---
                             auto result = std::make_shared<MultiValue>();
                             result->elems_.reserve(n);

                             for (double val : x)
                              result->elems_.push_back(Value(val));

                             return Value(result);
                            }};

  // 行列のフロベニウスノルム
  cfg.functions["mnormf"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mnormf requires matrix argument", ctx.pos);
                              const auto &matrix = v[0].asMultiRef(ctx.pos);

                              double sum = 0.0;
                              for (size_t i = 0; i < matrix.size(); ++i) {
                               for (size_t j = 0; j < matrix[i].asMultiRef(ctx.pos).size(); ++j) {
                                double val = matrix[i].asMultiRef(ctx.pos)[j].asScalar(ctx.pos);
                                sum += val * val;
                               }
                              }
                              return Value(std::sqrt(sum));
                             }};

  // 行列のフロベニウスノルム（別名）
  cfg.functions["mnorm"] = cfg.functions["mnormf"];

  // 行列の最小固有値
  cfg.functions["mmaxeigen_power"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                       std::vector<std::vector<double>> A = toMatrix(v[0], ctx);

                                       size_t m = A.size();
                                       if (m == 0 || m != A[0].size()) throw CalcError(CalcErrorType::DomainError, "matrix must be square", ctx.pos);

                                       // 1×1 特殊ケース
                                       if (m == 1) {
                                        if (!std::isfinite(A[0][0])) throw CalcError(CalcErrorType::DomainError, "invalid number", ctx.pos);

                                        return Value(A[0][0]);
                                       }

                                       std::vector<double> eigenvec;
                                       double lambda = power_method(A, eigenvec);
                                       return Value(lambda);
                                      }};
  cfg.functions["mmaxeigen"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                 std::vector<std::vector<double>> A = toMatrix(v[0], ctx);

                                 size_t m = A.size();
                                 if (m == 0 || m != A[0].size()) throw CalcError(CalcErrorType::DomainError, "matrix must be square", ctx.pos);

                                 // 1×1 特殊ケース
                                 if (m == 1) {
                                  if (!std::isfinite(A[0][0])) throw CalcError(CalcErrorType::DomainError, "invalid number", ctx.pos);

                                  return Value(A[0][0]);
                                 }

                                 auto eig = eigenvalues(A);
                                 if (eig.empty()) throw CalcError(CalcErrorType::DomainError, "eigenvalues empty", ctx.pos);
                                 return *std::max_element(eig.begin(), eig.end());
                                }};

  // 行列の特異値の和
  cfg.functions["msumsv"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              auto A = toMatrix(v[0], ctx);

                              if (A.empty()) throw CalcError(CalcErrorType::DomainError, "matrix is empty", ctx.pos);

                              auto S = singular_values(A);

                              double sum = 0.0;
                              for (double s : S)
                               sum += s;

                              return Value(sum);
                             }};

  // 行列の対称性判定
  cfg.functions["misse_symmetric"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                       if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "misse symmetric requires matrix argument", ctx.pos);
                                       const auto &matrix = v[0].asMultiRef(ctx.pos);

                                       size_t rows = matrix.size();
                                       size_t cols = matrix[0].asMultiRef(ctx.pos).size();

                                       if (rows != cols) return Value(false);

                                       for (size_t i = 0; i < rows; ++i) {
                                        for (size_t j = i + 1; j < cols; ++j) {
                                         double a = matrix[i].asMultiRef(ctx.pos)[j].asScalar(ctx.pos);
                                         double b = matrix[j].asMultiRef(ctx.pos)[i].asScalar(ctx.pos);
                                         if (std::abs(a - b) > 1e-10) { return Value(false); }
                                        }
                                       }
                                       return Value(true);
                                      }};

  // 行列の正定性判定
  cfg.functions["mispd"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mispd requires matrix argument", ctx.pos);

                             const auto &matrix = v[0].asMultiRef(ctx.pos);

                             // 行列の次元確認
                             if (matrix.empty()) throw CalcError(CalcErrorType::InvalidOperation, "std::vector<std::vector<double>> is empty", ctx.pos);

                             size_t n = matrix.size(); // 行列の次元

                             // 正方行列でない場合はfalseを返す
                             for (size_t i = 0; i < n; ++i) {
                              if (matrix[i].size() != n) throw CalcError(CalcErrorType::InvalidOperation, "std::vector<std::vector<double>> must be square", ctx.pos);
                             }

                             // 行列のコピーを作成（元の行列を変更しない）
                             std::vector<std::vector<double>> A(n, std::vector<double>(n));
                             for (size_t i = 0; i < n; ++i) {
                              const auto &row = matrix[i].asMultiRef(ctx.pos);
                              for (size_t j = 0; j < n; ++j) {
                               A[i][j] = row[j].asScalar(ctx.pos);
                              }
                             }

                             // Cholesky分解を試行する
                             try {
                              std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

                              for (size_t i = 0; i < n; ++i) {
                               double sum = 0.0;
                               for (size_t k = 0; k < i; ++k) {
                                sum += L[i][k] * L[i][k];
                               }
                               double diag = A[i][i] - sum;

                               if (diag < -1e-12) { return Value(false); }
                               diag = std::max(0.0, diag);

                               L[i][i] = std::sqrt(diag);

                               if (L[i][i] < 1e-12) { return Value(false); }

                               for (size_t j = i + 1; j < n; ++j) {
                                sum = 0.0;
                                for (size_t k = 0; k < i; ++k) {
                                 sum += L[j][k] * L[i][k];
                                }
                                L[j][i] = (A[j][i] - sum) / L[i][i];
                               }
                              }

                              return Value(true);

                             } catch (const std::exception &) { return Value(false); }
                            }};
 }

 void registerComplex(SystemConfig &cfg) {
  cfg.functions["re"] = {1, 1, [](auto &v, auto &ctx) -> Value { return v[0].isComplex() ? v[0].asComplex(ctx.pos).real() : v[0].asScalar(ctx.pos); }};
  cfg.functions["real"] = cfg.functions["re"];
  cfg.functions["im"] = {1, 1, [](auto &v, auto &ctx) -> Value { return v[0].isComplex() ? v[0].asComplex(ctx.pos).imag() : 0.0; }};
  cfg.functions["imag"] = cfg.functions["im"];
  cfg.functions["arg"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           auto z = v[0].toComplex(ctx.pos);
                           if (std::abs(z) == 0) throwDomain(ctx.pos);
                           return std::arg(z * 180.0 / PI) * rad2deg;
                          }};
  cfg.functions["conj"] = {1, 1, [](auto &v, auto &ctx) -> Value { return std::conj(requireComplex(v[0], ctx.pos)); }};
  cfg.functions["polar"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             double r = v[0].asScalar(ctx.pos);
                             double th = v[1].asScalar(ctx.pos);
                             return Complex(r * std::cos(toRad(th)), r * std::sin(toRad(th)));
                            }};
  cfg.functions["rect"] = cfg.functions["polar"];
  cfg.functions["cis"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           double th = v[0].asScalar(ctx.pos);
                           return Complex(std::cos(toRad(th)), std::sin(toRad(th)));
                          }};
  cfg.functions["proj"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            Complex z = v[0].isComplex() ? v[0].asComplex(ctx.pos) : Complex(v[0].asScalar(ctx.pos), 0.0);
                            return std::proj(z);
                           }};
  cfg.functions["unit"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            Complex z = v[0].isComplex() ? v[0].asComplex(ctx.pos) : Complex(v[0].asScalar(ctx.pos), 0.0);
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
                            if (v.empty()) return make(0.0, 1.0);                        // 引数なし → [0, 1)
                            if (v.size() == 1) return make(0.0, v[0].asScalar(ctx.pos)); // 引数1つ → [0, hi)
                            return make(v[0].asScalar(ctx.pos), v[1].asScalar(ctx.pos)); // 引数2つ → [lo, hi)
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
                             double mu = (v.size() >= 1) ? v[0].asScalar(ctx.pos) : 0.0;
                             double sigma = (v.size() >= 2) ? v[1].asScalar(ctx.pos) : 1.0;
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
  cfg.functions["fma"] = {3, 3, [](auto &v, auto &ctx) -> Value { return std::fma(v[0].asScalar(ctx.pos), v[1].asScalar(ctx.pos), v[2].asScalar(ctx.pos)); }};
 }

 void registerAreaVol(SystemConfig &cfg) {
  cfg.functions["area_circle"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                   double d = v[0].asScalar(ctx.pos);
                                   if (d < 0) throwDomain(ctx.pos, "diameter must be >= 0");
                                   return PI * d * d / 4.0;
                                  }};

  cfg.functions["area_triangle"] = {2, 3, [](auto &v, auto &ctx) -> Value {
                                     if (v.size() == 2) return 0.5 * v[0].asScalar(ctx.pos) * v[1].asScalar(ctx.pos);
                                     double a = v[0].asScalar(ctx.pos);
                                     double b = v[1].asScalar(ctx.pos);
                                     double c = v[2].asScalar(ctx.pos);

                                     if (a <= 0 || b <= 0 || c <= 0 || a + b <= c || b + c <= a || c + a <= b) { throw CalcError(CalcErrorType::DomainError, "area_triangle: invalid triangle sides", ctx.pos); }

                                     double s = 0.5 * (a + b + c);
                                     double area = std::sqrt(s * (s - a) * (s - b) * (s - c));
                                     return area;
                                    }};

  cfg.functions["area_trapezoid"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                      double a = v[0].asScalar(ctx.pos);
                                      double b = v[1].asScalar(ctx.pos);
                                      double h = v[2].asScalar(ctx.pos);
                                      return 0.5 * (a + b) * h;
                                     }};

  // 多角形の座標から面積（Shoelace formula）
  cfg.functions["area_polygon"] = {6, -1, area_polygon_impl};
  cfg.functions["vol_cylinder"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                    double d = v[0].asScalar(ctx.pos);
                                    double h = v[1].asScalar(ctx.pos);
                                    if (d < 0 || h < 0) throwDomain(ctx.pos, "must be >= 0");
                                    return PI * (d * d / 4.0) * h;
                                   }};

  cfg.functions["vol_cone"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                double d = v[0].asScalar(ctx.pos);
                                double h = v[1].asScalar(ctx.pos);
                                if (d < 0 || h < 0) throwDomain(ctx.pos);
                                return PI * (d * d / 4.0) * h / 3.0;
                               }};

  cfg.functions["vol_sphere"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  double d = v[0].asScalar(ctx.pos);
                                  double r = d / 2.0;
                                  return 4.0 / 3.0 * PI * r * r * r;
                                 }};

  cfg.functions["vol_prism"] = {6, -1, [](auto &v, auto &ctx) -> Value {
                                 if (v.size() < 7 || (v.size() - 1) % 2 != 0) throw CalcError(CalcErrorType::DomainError, "vol_prism: need at least 3 base points + height", ctx.pos);

                                 double h = v.back().asScalar(ctx.pos);
                                 std::vector<Value> coords(v.begin(), v.end() - 1);

                                 return area_polygon_impl(coords, ctx).asScalar(ctx.pos) * h;
                                }};
 }

 void registerEngineering(SystemConfig &cfg) { // stress(F, A) = F/A
  cfg.functions["stress"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              double F = v[0].asScalar(ctx.pos);
                              double A = v[1].asScalar(ctx.pos);
                              if (A == 0.0) throwDomain(ctx.pos);
                              return F / A;
                             }};

  // strain(dL, L) = dL/L
  cfg.functions["strain"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              double dL = v[0].asScalar(ctx.pos);
                              double L = v[1].asScalar(ctx.pos);
                              if (L == 0.0) throwDomain(ctx.pos);
                              return dL / L;
                             }};

  // young(sigma, eps) = sigma/eps
  cfg.functions["young"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             double sigma = v[0].asScalar(ctx.pos);
                             double eps = v[1].asScalar(ctx.pos);
                             if (eps == 0.0) throwDomain(ctx.pos);
                             return sigma / eps;
                            }};

  // moment_rect(b,h) = b*h^3/12
  cfg.functions["moment_rect"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                   double b = v[0].asScalar(ctx.pos);
                                   double h = v[1].asScalar(ctx.pos);
                                   return b * h * h * h / 12.0;
                                  }};

  // moment_circle(d) = pi*d^4/64
  cfg.functions["moment_circle"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                     double d = v[0].asScalar(ctx.pos);
                                     return PI * d * d * d * d / 64.0;
                                    }};

  // sectionmod_rect(b,h) = b*h^2/6
  cfg.functions["sectionmod_rect"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                       double b = v[0].asScalar(ctx.pos);
                                       double h = v[1].asScalar(ctx.pos);
                                       return b * h * h / 6.0;
                                      }};

  // sectionmod_circle(d) = pi*d^3/32
  cfg.functions["sectionmod_circle"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                         double d = v[0].asScalar(ctx.pos);
                                         return PI * d * d * d / 32.0;
                                        }};

  // torsion_J_circle(d) = pi*d^4/32
  cfg.functions["torsion_J_circle"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                        double d = v[0].asScalar(ctx.pos);
                                        return PI * d * d * d * d / 32.0;
                                       }};

  // polarZ_circle(d) = pi*d^3/16
  cfg.functions["polarZ_circle"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                     double d = v[0].asScalar(ctx.pos);
                                     return PI * d * d * d / 16.0;
                                    }};

  // bolt_stress(F, d) = F / (pi d^2 / 4)
  cfg.functions["bolt_stress"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                   double F = v[0].asScalar(ctx.pos);
                                   double d = v[1].asScalar(ctx.pos);
                                   if (d <= 0.0) throwDomain(ctx.pos);
                                   double A = PI * d * d / 4.0;
                                   return F / A;
                                  }};

  // torque_from_preload(F, d, K) = K*F*d
  cfg.functions["torque_from_preload"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                           double F = v[0].asScalar(ctx.pos);
                                           double d = v[1].asScalar(ctx.pos);
                                           double K = v[2].asScalar(ctx.pos);
                                           return K * F * d;
                                          }};

  // preload_from_torque(T, d, K) = T/(K*d)
  cfg.functions["preload_from_torque"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                           double T = v[0].asScalar(ctx.pos);
                                           double d = v[1].asScalar(ctx.pos);
                                           double K = v[2].asScalar(ctx.pos);
                                           if (K == 0.0 || d == 0.0) throwDomain(ctx.pos);
                                           return T / (K * d);
                                          }};

  // friction(mu, N) = mu*N
  cfg.functions["friction"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                double mu = v[0].asScalar(ctx.pos);
                                double N = v[1].asScalar(ctx.pos);
                                return mu * N;
                               }};
 }

 void registerMoldInjection(SystemConfig &cfg) {
  cfg.functions["mold_clamp"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                  double P_avg = v[0].asScalar(ctx.pos);  // 平均型内圧
                                  double A_proj = v[1].asScalar(ctx.pos); // 投影面積
                                  return P_avg * A_proj;                  // F_clamp = P_avg * A_proj
                                 }};

  cfg.functions["mold_clamp_safe"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                       double SF = v[0].asScalar(ctx.pos);
                                       double P_avg = v[1].asScalar(ctx.pos);
                                       double A_proj = v[2].asScalar(ctx.pos);
                                       return SF * P_avg * A_proj; // F_machine >= SF * P_avg * A_proj
                                      }};
  cfg.functions["mold_Pinj"] = {4, 4, [](auto &v, auto &ctx) -> Value {
                                 double P_cav = v[0].asScalar(ctx.pos);
                                 double dP_runner = v[1].asScalar(ctx.pos);
                                 double dP_gate = v[2].asScalar(ctx.pos);
                                 double dP_nozzle = v[3].asScalar(ctx.pos);
                                 return P_cav + dP_runner + dP_gate + dP_nozzle;
                                }};
  cfg.functions["mold_flowrate"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                     double V = v[0].asScalar(ctx.pos);      // 製品体積
                                     double t_fill = v[1].asScalar(ctx.pos); // 充填時間
                                     if (t_fill <= 0.0) throwDomain(ctx.pos);
                                     return V / t_fill; // Q = V / t_fill
                                    }};

  cfg.functions["mold_gate_velocity"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                          double Q = v[0].asScalar(ctx.pos);      // 体積流量
                                          double A_gate = v[1].asScalar(ctx.pos); // ゲート断面積
                                          return Q / A_gate;                      // v_gate = Q / A_gate
                                         }};
  cfg.functions["mold_shear_gate"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                       double Q = v[0].asScalar(ctx.pos);
                                       double b = v[1].asScalar(ctx.pos);
                                       double h = v[2].asScalar(ctx.pos);
                                       return 6 * Q / (b * h * h);
                                      }};

  cfg.functions["mold_shear_runner"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                         double Q = v[0].asScalar(ctx.pos);
                                         double D = v[1].asScalar(ctx.pos);
                                         return 32 * Q / (PI * D * D * D);
                                        }};
  cfg.functions["mold_pressure_loss_runner"] = {4, 4, [](auto &v, auto &ctx) -> Value {
                                                 double mu = v[0].asScalar(ctx.pos); // 見かけ粘度
                                                 double L = v[1].asScalar(ctx.pos);  // 流路長さ
                                                 double Q = v[2].asScalar(ctx.pos);  // 流量
                                                 double D = v[3].asScalar(ctx.pos);  // ランナー径
                                                 return 128 * mu * L * Q / (PI * std::pow(D, 4));
                                                }};
  cfg.functions["mold_eject_friction"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                           double mu = v[0].asScalar(ctx.pos);
                                           double N = v[1].asScalar(ctx.pos);
                                           return mu * N;
                                          }};

  cfg.functions["mold_eject_contact"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                          double p_contact = v[0].asScalar(ctx.pos);
                                          double A_contact = v[1].asScalar(ctx.pos);
                                          return p_contact * A_contact;
                                         }};

  cfg.functions["mold_eject_total"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                        double mu = v[0].asScalar(ctx.pos);
                                        double p_contact = v[1].asScalar(ctx.pos);
                                        double A_contact = v[2].asScalar(ctx.pos);
                                        return mu * p_contact * A_contact;
                                       }};
  cfg.functions["mold_mu_tex"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                   double k = v[0].asScalar(ctx.pos);
                                   double h = v[1].asScalar(ctx.pos); // シボ深さ [mm]
                                   return k * h;
                                  }};
  cfg.functions["mold_pin_stress"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                       double F_eject = v[0].asScalar(ctx.pos);
                                       double n = v[1].asScalar(ctx.pos);
                                       double A_pin = v[2].asScalar(ctx.pos);
                                       return (F_eject / n) / A_pin;
                                      }};
  cfg.functions["mold_plate_deflection"] = {5, 5, [](auto &v, auto &ctx) -> Value {
                                             double K = v[0].asScalar(ctx.pos); // 境界条件＋形状込み係数
                                             double P = v[1].asScalar(ctx.pos); // 型内圧
                                             double a = v[2].asScalar(ctx.pos); // 支点間スパン
                                             double E = v[3].asScalar(ctx.pos); // ヤング率
                                             double t = v[4].asScalar(ctx.pos); // 板厚
                                             return K * P * std::pow(a, 4) / (E * std::pow(t, 3));
                                            }};
 }

 void registerFinance(SystemConfig &cfg) {
  // ------------------ 基本キャッシュフロー ------------------

  // 将来価値 FV
  cfg.functions["fin_fv"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                              double rate = v[0].asScalar(ctx.pos);
                              double nper = v[1].asScalar(ctx.pos);
                              double pmt = v[2].asScalar(ctx.pos);
                              if (rate == 0) return pmt * nper;
                              return pmt * (std::pow(1.0 + rate, nper) - 1.0) / rate;
                             }};

  // 現在価値 PV
  cfg.functions["fin_pv"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                              double rate = v[0].asScalar(ctx.pos);
                              double nper = v[1].asScalar(ctx.pos);
                              double pmt = v[2].asScalar(ctx.pos);
                              if (rate == 0) return pmt * nper;
                              return pmt * (1.0 - std::pow(1.0 + rate, -nper)) / rate;
                             }};

  // 支払額 PMT
  cfg.functions["fin_pmt"] = {3, 3, [](const auto &v, auto &ctx) -> Value {
                               double rate = v[0].asScalar(ctx.pos);
                               double nper = v[1].asScalar(ctx.pos);
                               double pv = v[2].asScalar(ctx.pos);
                               if (nper <= 0) throwDomain(ctx.pos, "nper must be > 0");
                               if (rate == 0.0) return pv / nper;
                               double factor = 1.0;
                               // 高精度累乗（pow誤差対策）
                               for (int i = 0; i < int(nper); ++i)
                                factor *= (1.0 + rate);

                               return pv * rate * factor / (factor - 1.0);
                              }};

  // 支払総額
  cfg.functions["fin_total_payment"] = {3, 3, [](const auto &v, auto &ctx) -> Value {
                                         double rate = v[0].asScalar(ctx.pos);
                                         double nper = v[1].asScalar(ctx.pos);
                                         double pv = v[2].asScalar(ctx.pos);
                                         if (nper <= 0) throwDomain(ctx.pos);

                                         double pmt;
                                         if (rate == 0.0) pmt = pv / nper;
                                         else {
                                          double factor = 1.0;
                                          for (int i = 0; i < int(nper); ++i)
                                           factor *= (1.0 + rate);
                                          pmt = pv * rate * factor / (factor - 1.0);
                                         }

                                         return pmt * nper;
                                        }};

  // 元金返済額 PPMT（期指定）
  cfg.functions["fin_ppmt"] = {4, 4, [](auto &v, auto &ctx) -> Value {
                                double rate = v[0].asScalar(ctx.pos);
                                double nper = v[1].asScalar(ctx.pos);
                                double per = v[2].asScalar(ctx.pos);
                                double pv = v[3].asScalar(ctx.pos);
                                if (per < 1 || per > nper) throwDomain(ctx.pos);
                                double pmt = (rate == 0) ? pv / nper : pv * rate / (1.0 - std::pow(1.0 + rate, -nper));
                                double interest = (pv - pmt * (per - 1.0)) * rate;
                                return pmt - interest;
                               }};

  // 利息返済額 IPMT（期指定）
  cfg.functions["fin_ipmt"] = {4, 4, [](auto &v, auto &ctx) -> Value {
                                double rate = v[0].asScalar(ctx.pos);
                                double nper = v[1].asScalar(ctx.pos);
                                double per = v[2].asScalar(ctx.pos);
                                double pv = v[3].asScalar(ctx.pos);
                                if (per < 1 || per > nper) throwDomain(ctx.pos);
                                double pmt = (rate == 0) ? pv / nper : pv * rate / (1.0 - std::pow(1.0 + rate, -nper));
                                double principal_paid = pmt * (per - 1.0);
                                double remaining = pv - principal_paid;
                                return remaining * rate;
                               }};

  // ------------------ NPV / IRR ------------------

  // NPV
  cfg.functions["fin_npv"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                               if (v.size() <= 2) throwDomain(ctx.pos);
                               double rate = v[0].asScalar(ctx.pos);
                               std::vector<double> cf = collectReals(v, ctx);
                               cf.erase(cf.begin()); // rateを除去
                               double sum = 0;
                               for (size_t t = 0; t < cf.size(); ++t)
                                sum += cf[t] / std::pow(1 + rate, t + 1);
                               return sum;
                              }};

  // IRR
  cfg.functions["fin_irr"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                               if (v.size() < 2) throw CalcError(CalcErrorType::DomainError, "fin_irr: at least 2 cash flows required", ctx.pos);

                               std::vector<double> cf;
                               for (auto &x : v)
                                cf.push_back(x.asScalar(ctx.pos));

                               auto f = [&](double r) {
                                double sum = 0.0;
                                for (size_t t = 0; t < cf.size(); ++t)
                                 sum += cf[t] / std::pow(1.0 + r, t);
                                return sum;
                               };

                               // 全区間スキャンで符号変化を見つける
                               double step = 0.0001;
                               double start = -0.9999, end = 10.0;
                               std::vector<std::pair<double, double>> brackets;

                               double prev = f(start);
                               for (double x = start + step; x <= end; x += step) {
                                double curr = f(x);
                                if (prev * curr <= 0) brackets.push_back({x - step, x});
                                prev = curr;
                               }

                               if (brackets.empty()) throw CalcError(CalcErrorType::NonConvergence, "fin_irr: no root bracketed", ctx.pos);

                               // 最初の正の IRR を返す
                               auto brent = [&](double a, double b) {
                                double fa = f(a), fb = f(b);
                                if (fa * fb > 0) throw CalcError(CalcErrorType::NonConvergence, "Brent method not bracketed", ctx.pos);

                                double c = a, fc = fa, d = 0;
                                bool mflag = true;
                                double tol = 1e-12;

                                for (int iter = 0; iter < 512; ++iter) {
                                 double s;
                                 if (fa != fc && fb != fc) s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb));
                                 else s = b - fb * (b - a) / (fb - fa);

                                 double midpoint = (3 * a + b) / 4;
                                 bool needBisect =
                                     !((s > midpoint && s < b) || (s < midpoint && s > b)) || (mflag ? std::abs(s - b) >= std::abs(b - c) / 2 : std::abs(s - b) >= std::abs(c - d) / 2) || (mflag ? std::abs(b - c) < tol : std::abs(c - d) < tol);

                                 if (needBisect) s = (a + b) / 2, mflag = true;
                                 else mflag = false;

                                 double fs = f(s);
                                 d = c;
                                 c = b;
                                 fc = fb;

                                 if (fa * fs < 0) {
                                  b = s;
                                  fb = fs;
                                 } else {
                                  a = s;
                                  fa = fs;
                                 }

                                 if (std::abs(fa) < std::abs(fb)) {
                                  std::swap(a, b);
                                  std::swap(fa, fb);
                                 }

                                 if (std::abs(b - a) < tol) return b;
                                }
                                throw CalcError(CalcErrorType::NonConvergence, "fin_irr: Brent did not converge", ctx.pos);
                               };

                               for (auto &[low, high] : brackets) {
                                try {
                                 double r = brent(low, high);
                                 if (r >= -1.0) return r; // 正の IRR 優先
                                } catch (...) { continue; }
                               }
                               throw CalcError(CalcErrorType::NonConvergence, "fin_irr: root not found", ctx.pos);
                              }};

  // MIRR
  cfg.functions["fin_mirr"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                if (v.size() < 3) throw CalcError(CalcErrorType::DomainError, "fin_mirr: at least reinvestRate + 2 cash flows required", ctx.pos);

                                double reinvestRate = v[0].asScalar(ctx.pos);
                                std::vector<double> cf;
                                for (size_t i = 1; i < v.size(); ++i)
                                 cf.push_back(v[i].asScalar(ctx.pos));
                                size_t n = cf.size();

                                double positiveFV = 0.0, negativePV = 0.0;
                                for (size_t t = 0; t < n; ++t) {
                                 if (cf[t] >= 0) positiveFV += cf[t] * std::pow(1.0 + reinvestRate, n - t - 1);
                                 else negativePV += cf[t] / std::pow(1.0 + reinvestRate, t);
                                }

                                if (positiveFV <= 0.0 || negativePV >= 0.0) return reinvestRate;

                                return std::pow(-positiveFV / negativePV, 1.0 / n) - 1.0;
                               }};

  // NPER (期間)
  cfg.functions["fin_nper"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                double rate = v[0].asScalar(ctx.pos);
                                double pmt = v[1].asScalar(ctx.pos);
                                double pv = v[2].asScalar(ctx.pos);
                                if (pmt == 0) throwDomain(ctx.pos);
                                if (rate == 0) return pv / pmt;
                                double denom = 1 - pv * rate / pmt;
                                if (denom <= 0) throwDomain(ctx.pos);
                                return std::log(denom) / std::log(1 + rate);
                               }};

  // RATE (利率)
  cfg.functions["fin_rate"] = {3, 3, [&](const auto &v, auto &ctx) -> Value {
                                double nper = v[0].asScalar(ctx.pos);
                                double pmt = v[1].asScalar(ctx.pos);
                                double pv = v[2].asScalar(ctx.pos);

                                if (nper <= 0) throwDomain(ctx.pos);
                                if (pmt == 0.0 && pv == 0.0) return 0.0;
                                if (pmt == 0.0 && pv != 0.0) throwDomain(ctx.pos);

                                auto f = [&](double r) {
                                 if (r == 0.0) return pv + pmt * nper;
                                 return pv * std::pow(1.0 + r, nper) + pmt * (std::pow(1.0 + r, nper) - 1.0) / r;
                                };

                                return brent(f, -0.999, 1e6, ctx);
                               }};

  // CUMIPMT (累積利息)
  // 1. その期の利息を計算
  // 2. 支払額から元本返済を計算
  // 3. 元本返済で残高を減らす
  // 4. type = 1 の場合のみ期首補正
  cfg.functions["fin_cumipmt"] = {6, 6, [](auto &v, auto &ctx) -> Value {
                                   double rate = v[0].asScalar(ctx.pos);
                                   int nper = (int)v[1].asScalar(ctx.pos);
                                   double pv = v[2].asScalar(ctx.pos);
                                   int start = (int)v[3].asScalar(ctx.pos);
                                   int end = (int)v[4].asScalar(ctx.pos);
                                   int type = (int)v[5].asScalar(ctx.pos);

                                   if (rate < 0 || nper <= 0) throwDomain(ctx.pos);
                                   if (start < 1 || end < start || end > nper) throwDomain(ctx.pos);

                                   const double EPS = 1e-15;

                                   double factor = std::pow(1.0 + rate, nper);
                                   double pmt = (std::fabs(rate) < EPS) ? -pv / nper : -pv * rate * factor / (factor - 1.0);

                                   double balance = pv;

                                   double sum = 0.0;
                                   double c = 0.0;

                                   auto kahan_add = [&](double x) {
                                    double y = x - c;
                                    double t = sum + y;
                                    c = (t - sum) - y;
                                    sum = t;
                                   };

                                   for (int k = 1; k <= end; k++) {

                                    double interest = 0.0;

                                    if (!(type == 1 && k == 1)) { interest = balance * rate; }

                                    if (k >= start) kahan_add(interest);

                                    balance = balance * (1.0 + rate) + pmt;
                                   }

                                   return sum;
                                  }};

  // CUMPRINC (累積元金)
  cfg.functions["fin_cumprinc"] = {6, 6, [](auto &v, auto &ctx) -> Value {
                                    double rate = v[0].asScalar(ctx.pos);
                                    int nper = (int)v[1].asScalar(ctx.pos);
                                    double pv = v[2].asScalar(ctx.pos);
                                    int start = (int)v[3].asScalar(ctx.pos);
                                    int end = (int)v[4].asScalar(ctx.pos);
                                    int type = (int)v[5].asScalar(ctx.pos);

                                    if (rate < 0 || nper <= 0) throwDomain(ctx.pos);
                                    if (start < 1 || end < start || end > nper) throwDomain(ctx.pos);

                                    double pmt = (std::fabs(rate) < 1e-15) ? -pv / nper : -pv * rate * std::pow(1 + rate, nper) / (std::pow(1 + rate, nper) - 1);

                                    double balance = pv;
                                    double sum = 0;

                                    for (int k = 1; k <= end; k++) {

                                     double prev_balance = balance;

                                     if (type == 1 && k == 1) balance = balance + pmt;
                                     else balance = balance * (1 + rate) + pmt;

                                     double principal = prev_balance - balance;

                                     if (k >= start) sum += principal;
                                    }

                                    return sum;
                                   }};

  // EFFECTIVE RATE / NOMINAL RATE / CAGR
  cfg.functions["fin_effective_rate"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                          double nominal = v[0].asScalar(ctx.pos);
                                          double npery = v[1].asScalar(ctx.pos);
                                          return std::pow(1 + nominal / npery, npery) - 1.0;
                                         }};
  cfg.functions["fin_nominal_rate"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                        double effective = v[0].asScalar(ctx.pos);
                                        double npery = v[1].asScalar(ctx.pos);
                                        return npery * (std::pow(1 + effective, 1.0 / npery) - 1.0);
                                       }};

  cfg.functions["fin_cagr"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                double start = v[0].asScalar(ctx.pos);
                                double end = v[1].asScalar(ctx.pos);
                                double n = v[2].asScalar(ctx.pos);
                                if (start <= 0 || end <= 0 || n <= 0) throwDomain(ctx.pos); // 値がゼロ以下や期間が0以下は定義不能
                                // CAGR = (終値/初値)^(1/n) - 1
                                return std::pow(end / start, 1.0 / n) - 1.0;
                               }};
 }

 void registerOthers(SystemConfig &cfg) {
  cfg.functions["cnst"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            // v[0] が string であることを要求
                            const std::string &name = v[0].asString(ctx.pos);
                            auto it = constants_dic.find(name);
                            if (it == constants_dic.end()) { throw CalcError(CalcErrorType::UnknownIdentifier, "unknown constant: " + name, ctx.pos); }
                            return it->second;
                           }};
  cfg.functions["fib"] = {1, 1, [](auto &v, auto &ctx) -> Value { return (double)fibULL(requireInt(v[0], ctx.pos), ctx.pos); }};
  cfg.functions["isprime"] = {1, 1, [](auto &v, auto &ctx) -> Value { return isPrimeLL(requireInt(v[0], ctx.pos)) ? 1.0 : 0.0; }};
  cfg.functions["nextprime"] = {1, 1, [](auto &v, auto &ctx) {
                                 long long n = requireInt(v[0], ctx.pos) + 1;
                                 while (!isPrimeLL(n))
                                  ++n;
                                 return (double)n;
                                }};
  cfg.functions["prevprime"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                 int n = (int)std::floor(v[0].asScalar(ctx.pos));
                                 while (--n >= 2)
                                  if (isPrimeLL(n)) return (double)n;
                                 throwDomain(ctx.pos);
                                }};
  cfg.functions["if"] = {3, 3, [](auto &v, auto &ctx) -> Value { return (v[0].asScalar(ctx.pos) != 0.0 ? v[1] : v[2]); }};

  cfg.functions["convert"] = {2, 3, [](auto &v, auto &ctx) -> Value {
                               double value = v[0].asScalar(ctx.pos);
                               if (v.size() == 2) {
                                // 引数2つ: SI単位系に変換
                                std::string unit = v[1].asString(ctx.pos);
                                // SI単位系に変換
                                if (unit == "m") return Value(value);
                                else if (unit == "kg") return Value(value);
                                else if (unit == "s") return Value(value);
                                else if (unit == "K") return Value(value);
                                else if (unit == "A") return Value(value);
                                else if (unit == "mol") return Value(value);
                                else if (unit == "cd") return Value(value);
                                else {
                                 // 他の単位もSI単位系に変換
                                 for (const auto &[quantity, units] : unitTable) {
                                  auto it = units.find(unit);
                                  if (it != units.end()) {
                                   // SI単位系に変換
                                   double si_value = value * it->second;
                                   // SI単位系の値を返す
                                   return Value(si_value);
                                  }
                                 }
                                 throw CalcError(CalcErrorType::UnknownIdentifier, "unknown unit: " + unit, ctx.pos);
                                }
                               } else {
                                // 引数3つ: 通常の変換
                                std::string from_unit = v[1].asString(ctx.pos);
                                std::string to_unit = v[2].asString(ctx.pos);

                                // 温度特殊処理(羊のせいでえらい目にあった)
                                if ((from_unit == "C" || from_unit == "F" || from_unit == "K") && (to_unit == "C" || to_unit == "F" || to_unit == "K")) {

                                 // ケルビンに変換
                                 double kelvin;
                                 if (from_unit == "C") {
                                  kelvin = value + 273.15;
                                 } else if (from_unit == "F") {
                                  kelvin = (value - 32) * 5 / 9 + 273.15;
                                 } else {
                                  kelvin = value;
                                 }

                                 // ターゲットに変換
                                 if (to_unit == "C") {
                                  return Value(kelvin - 273.15);
                                 } else if (to_unit == "F") {
                                  return Value((kelvin - 273.15) * 9 / 5 + 32);
                                 } else {
                                  return Value(kelvin);
                                 }
                                }

                                // 他の単位の変換
                                for (const auto &[quantity, units] : unitTable) {
                                 auto from_it = units.find(from_unit);
                                 auto to_it = units.find(to_unit);

                                 if (from_it != units.end() && to_it != units.end()) {
                                  // SI単位系に変換してから変換
                                  double si_value = value * from_it->second;
                                  double result = si_value / to_it->second;
                                  return Value(result);
                                 }
                                }

                                throw CalcError(CalcErrorType::UnknownIdentifier, "unknown unit combination: " + from_unit + " to " + to_unit, ctx.pos);
                               }
                              }};

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
  registerBasicMath(cfg);     // 基本関数
  registerExpLog(cfg);        // log系
  registerTrig(cfg);          // 三角関数
  registerHyperbolic(cfg);    // hyper三角関数
  registerGeoVec(cfg);        // 図形, ベクトル
  registerVector(cfg);        // 行列・線形代数
  registerSpecialTrig(cfg);   // hyper三角関数
  registerStatistics(cfg);    // 統計系
  registerComplex(cfg);       // 複素
  registerRandom(cfg);        // 乱数
  registerAreaVol(cfg);       // 面積, 体積
  registerEngineering(cfg);   // 機械設計
  registerMoldInjection(cfg); // 金型
  registerFinance(cfg);       // 財務系
  registerOthers(cfg);        // アミューズ
 }

 Value evaluateFunction(const std::string &name, const std::vector<Value> &args, FunctionContext &ctx) {
  auto it = ctx.cfg.functions.find(name);
  if (it == ctx.cfg.functions.end()) throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), ctx.pos);

  const FunctionDef &f = it->second;

  if (!f.validArgc(static_cast<int>(args.size()))) throw CalcError(CalcErrorType::InvalidArgument, errorMessage(CalcErrorType::InvalidArgument), ctx.pos);

  Value result = f.f(args, ctx);

  // 型を知らずに一括チェック
  result.validateFinite(ctx.pos);

  return result;
 }

} // namespace mm::cal