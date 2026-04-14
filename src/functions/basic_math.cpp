#include "functions/functions.hpp"
#include "unit_conv.hpp"

namespace mm::cal {

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

} // namespace mm::cal
