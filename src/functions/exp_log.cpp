#include "functions/functions.hpp"

namespace mm::cal {

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

                           // 1引数: ln(x)
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

                           // 2引数: log(base, x)
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
                             if (x <= -1.0) throw CalcError(CalcErrorType::DomainError, "DomainError: log1p: x <= -1", ctx.pos);

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

} // namespace mm::cal
