#include "functions/functions.hpp"

namespace mm::cal {

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

                                     if (a <= 0 || b <= 0 || c <= 0 || a + b <= c || b + c <= a || c + a <= b) { throw CalcError(CalcErrorType::DomainError, "DomainError: area_triangle: invalid triangle sides", ctx.pos); }

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
                                 if (v.size() < 7 || (v.size() - 1) % 2 != 0) throw CalcError(CalcErrorType::DomainError, "DomainError: vol_prism: need at least 3 base points + height", ctx.pos);

                                 double h = v.back().asScalar(ctx.pos);
                                 std::vector<Value> coords(v.begin(), v.end() - 1);

                                 return area_polygon_impl(coords, ctx).asScalar(ctx.pos) * h;
                                }};
 }

} // namespace mm::cal
