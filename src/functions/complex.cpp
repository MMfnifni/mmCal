#include "functions/functions.hpp"

namespace mm::cal {

 void registerComplex(SystemConfig &cfg) {
  cfg.functions["re"] = {1, 1, [](auto &v, auto &ctx) -> Value { return v[0].isComplex() ? v[0].asComplex(ctx.pos).real() : v[0].asScalar(ctx.pos); }};
  cfg.functions["real"] = cfg.functions["re"];
  cfg.functions["im"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                          if (v[0].isComplex()) return v[0].asComplex(ctx.pos).imag();
                          (void)v[0].asScalar(ctx.pos); // scalar以外はここでTypeError
                          return 0.0;
                         }};
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
                            if (a == 0.0) throw CalcError(CalcErrorType::DomainError, "DomainError: unit(0) is undefined", ctx.pos);
                            return z / a;
                           }};

  cfg.functions["csgn"] = cfg.functions["unit"];
  cfg.functions["mag"] = cfg.functions["abs"];
 }

} // namespace mm::cal
