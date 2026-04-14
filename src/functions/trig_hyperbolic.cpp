#include "functions/functions.hpp"

namespace mm::cal {

 void registerTrig(SystemConfig &cfg) {
  cfg.functions["sin"] = {1, 1, [=](auto &v, auto &ctx) -> Value { return realIfPossible(std::sin(toComplex(v[0], ctx.pos) * deg2rad)); }};
  cfg.functions["cos"] = {1, 1, [=](auto &v, auto &ctx) -> Value { return realIfPossible(std::cos(toComplex(v[0], ctx.pos) * deg2rad)); }};
  // cfg.functions["tan"] = {1, 1, [=](auto &v, auto &) -> Value { return std::tan(asComplex(v[0]) * deg2rad); }};
  cfg.functions["tan"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (v[0].isComplex()) return realIfPossible(std::tan(toComplex(v[0], ctx.pos) * deg2rad));
                           double r = requireReal(v[0], ctx.pos) * deg2rad;
                           double c = std::cos(r);
                           return (std::abs(c) < cnst_precision_inv) ? signedInfBy(std::sin(r)) : std::tan(r);
                          }};
  cfg.functions["cot"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (v[0].isComplex()) {
                            auto z = toComplex(v[0], ctx.pos) * deg2rad;
                            return realIfPossible(std::cos(z) / std::sin(z));
                           }
                           double r = requireReal(v[0], ctx.pos) * deg2rad;
                           double s = std::sin(r);
                           return (std::abs(s) < cnst_precision_inv) ? signedInfBy(std::cos(r)) : std::cos(r) / s;
                          }};
  cfg.functions["sec"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (v[0].isComplex()) return realIfPossible(Complex(1, 0) / std::cos(toComplex(v[0], ctx.pos) * deg2rad));
                           double r = requireReal(v[0], ctx.pos) * deg2rad;
                           double c = std::cos(r);
                           return (std::abs(c) < cnst_precision_inv) ? signedInfBy(std::sin(r)) : 1.0 / c;
                          }};
  cfg.functions["csc"] = {1, 1, [=](auto &v, auto &ctx) -> Value {
                           if (v[0].isComplex()) return realIfPossible(Complex(1, 0) / std::sin(toComplex(v[0], ctx.pos) * deg2rad));
                           double r = requireReal(v[0], ctx.pos) * deg2rad;
                           double s = std::sin(r);
                           return (std::abs(s) < cnst_precision_inv) ? signedInfBy(std::cos(r)) : 1.0 / s;
                          }};
  cfg.functions["asin"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (v[0].isComplex()) return realIfPossible(std::asin(toComplex(v[0], ctx.pos)) * rad2deg);
                            double x = requireReal(v[0], ctx.pos);
                            if (std::abs(x) <= 1.0) return std::asin(x) * rad2deg;
                            calcWarn(ctx.session.cfg, ctx.pos, "asin(|x|>1): complex principal value only");
                            return realIfPossible(std::asin(Complex(x, 0)) * rad2deg);
                           }};
  cfg.functions["acos"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (v[0].isComplex()) return realIfPossible(std::acos(toComplex(v[0], ctx.pos)) * rad2deg);
                            double x = requireReal(v[0], ctx.pos);
                            if (std::abs(x) <= 1.0) return std::acos(x) * rad2deg;
                            calcWarn(ctx.session.cfg, ctx.pos, "acos(|x|>1): complex principal value only");
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
                             calcWarn(ctx.session.cfg, ctx.pos, "acosh(x<1): complex principal value only");
                             return realIfPossible(std::acosh(Complex(x, 0)));
                            }};

  cfg.functions["atanh"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (v[0].isComplex()) return realIfPossible(std::atanh(toComplex(v[0], ctx.pos)));
                             double x = requireReal(v[0], ctx.pos);
                             if (std::abs(x) < 1.0) return std::atanh(x);
                             calcWarn(ctx.session.cfg, ctx.pos, "atanh(|x|>=1): complex principal value only");
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

} // namespace mm::cal
