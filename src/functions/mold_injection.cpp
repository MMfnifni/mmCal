#include "functions/functions.hpp"

namespace mm::cal {

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

} // namespace mm::cal
