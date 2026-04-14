#include "functions/functions.hpp"

namespace mm::cal {

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

  // young(sigma, cnst_precision_inv) = sigma/cnst_precision_inv
  cfg.functions["young"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             double sigma = v[0].asScalar(ctx.pos);
                             double cnst_precision_inv = v[1].asScalar(ctx.pos);
                             if (cnst_precision_inv == 0.0) throwDomain(ctx.pos);
                             return sigma / cnst_precision_inv;
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

} // namespace mm::cal
