#include "functions/functions.hpp"

namespace mm::cal {

 void registerFinance(SystemConfig &cfg) {
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

  // NPV
  cfg.functions["fin_npv"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                               if (v.size() <= 2) throwDomain(ctx.pos);

                               const long double rate = (long double)v[0].asScalar(ctx.pos);

                               std::vector<long double> terms;
                               terms.reserve(v.size() - 1);

                               for (size_t i = 1; i < v.size(); ++i) {
                                const long double cf = (long double)v[i].asScalar(ctx.pos);
                                const long double disc = std::expl((long double)i * std::log1pl(rate));
                                terms.push_back(cf / disc);
                               }

                               return (double)kahan_sum_ld(terms);
                              }};

  // IRR
  cfg.functions["fin_irr"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                               if (v.size() < 2) { throw CalcError(CalcErrorType::DomainError, "DomainError: fin_irr: at least 2 cash flows required", ctx.pos); }

                               std::vector<long double> cf;
                               for (auto &x : v) {
                                cf.push_back((long double)x.asScalar(ctx.pos));
                               }

                               auto f = [&](double r_in) -> double {
                                const long double r = (long double)r_in;
                                long double sum = 0.0L;
                                long double c = 0.0L;

                                for (size_t t = 0; t < cf.size(); ++t) {
                                 const long double denom = std::expl((long double)t * std::log1pl(r));
                                 const long double term = cf[t] / denom;

                                 long double y = term - c;
                                 long double tmp = sum + y;
                                 c = (tmp - sum) - y;
                                 sum = tmp;
                                }

                                return (double)sum;
                               };

                               const double step = 1e-4;
                               const double start = -0.9999;
                               const double end = 10.0;

                               std::vector<std::pair<double, double>> brackets;
                               double prev = f(start);

                               for (double x = start + step; x <= end; x += step) {
                                double curr = f(x);
                                if (prev * curr <= 0.0) { brackets.push_back({x - step, x}); }
                                prev = curr;
                               }

                               if (brackets.empty()) { throw CalcError(CalcErrorType::NonConvergence, "fin_irr: no root bracketed", ctx.pos); }

                               for (auto &[low, high] : brackets) {
                                try {
                                 double r = brent(f, low, high, ctx);
                                 if (r >= -1.0) return r;
                                } catch (...) { continue; }
                               }

                               throw CalcError(CalcErrorType::NonConvergence, "fin_irr: root not found", ctx.pos);
                              }};

  // MIRR
  cfg.functions["fin_mirr"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                if (v.size() < 3) { throw CalcError(CalcErrorType::DomainError, "DomainError: fin_mirr: at least reinvestRate + 2 cash flows required", ctx.pos); }

                                const long double reinvestRate = (long double)v[0].asScalar(ctx.pos);

                                std::vector<long double> cf;
                                for (size_t i = 1; i < v.size(); ++i) {
                                 cf.push_back((long double)v[i].asScalar(ctx.pos));
                                }

                                const int n = (int)cf.size();

                                long double positiveFV = 0.0L;
                                long double negativePV = 0.0L;

                                for (int t = 0; t < n; ++t) {
                                 if (cf[t] >= 0.0L) {
                                  positiveFV += cf[t] * pow1p_int_ld(reinvestRate, n - t - 1);
                                 } else {
                                  negativePV += cf[t] / pow1p_int_ld(reinvestRate, t);
                                 }
                                }

                                if (positiveFV <= 0.0L || negativePV >= 0.0L) { return (double)reinvestRate; }

                                return (double)(std::powl(-positiveFV / negativePV, 1.0L / (long double)n) - 1.0L);
                               }};

  // NPER (期間)
  cfg.functions["fin_nper"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                                const long double rate = (long double)v[0].asScalar(ctx.pos);
                                const long double pmt = (long double)v[1].asScalar(ctx.pos);
                                const long double pv = (long double)v[2].asScalar(ctx.pos);

                                if (nearly_zero_ld(pmt)) throwDomain(ctx.pos);

                                if (nearly_zero_ld(rate)) { return (double)(pv / pmt); }

                                const long double denom = 1.0L - pv * rate / pmt;
                                if (denom <= 0.0L) throwDomain(ctx.pos);

                                return (double)(std::logl(denom) / std::log1pl(rate));
                               }};

  // RATE (利率)
  cfg.functions["fin_rate"] = {3, 3, [&](const auto &v, auto &ctx) -> Value {
                                const long double nper = (long double)v[0].asScalar(ctx.pos);
                                const long double pmt = (long double)v[1].asScalar(ctx.pos);
                                const long double pv = (long double)v[2].asScalar(ctx.pos);

                                if (nper <= 0.0L) throwDomain(ctx.pos);
                                if (nearly_zero_ld(pmt) && nearly_zero_ld(pv)) return 0.0;
                                if (nearly_zero_ld(pmt) && !nearly_zero_ld(pv)) throwDomain(ctx.pos);

                                auto f = [&](double r_in) -> double {
                                 const long double r = (long double)r_in;
                                 if (nearly_zero_ld(r)) { return (double)(pv + pmt * nper); }

                                 const long double A = std::expl(nper * std::log1pl(r));
                                 return (double)(pv * A + pmt * (A - 1.0L) / r);
                                };

                                return brent(f, -0.999, 1e6, ctx);
                               }};

  // CUMIPMT (累積利息)
  // 1. その期の利息を計算
  // 2. 支払額から元本返済を計算
  // 3. 元本返済で残高を減らす
  // 4. type = 1 の場合のみ期首補正
  cfg.functions["fin_cumipmt"] = {6, 6, [](auto &v, auto &ctx) -> Value {
                                   const long double rate = (long double)v[0].asScalar(ctx.pos);
                                   const int nper = (int)v[1].asScalar(ctx.pos);
                                   const long double pv = (long double)v[2].asScalar(ctx.pos);
                                   const int start = (int)v[3].asScalar(ctx.pos);
                                   const int end = (int)v[4].asScalar(ctx.pos);
                                   const int type = (int)v[5].asScalar(ctx.pos);

                                   if (rate < 0.0L || nper <= 0) throwDomain(ctx.pos);
                                   if (start < 1 || end < start || end > nper) throwDomain(ctx.pos);
                                   if (!(type == 0 || type == 1)) throwDomain(ctx.pos);

                                   if (nearly_zero_ld(rate)) { return 0.0; }

                                   const long double pmt = fin_pmt_ld(rate, nper, pv, type); // 正の支払額
                                   long double balance = pv;

                                   long double sum = 0.0L;
                                   long double c = 0.0L;

                                   auto kahan_add = [&](long double x) {
                                    long double y = x - c;
                                    long double t = sum + y;
                                    c = (t - sum) - y;
                                    sum = t;
                                   };

                                   for (int k = 1; k <= end; ++k) {
                                    long double interest = 0.0L;
                                    long double principal = 0.0L;

                                    if (type == 1 && k == 1) {
                                     interest = 0.0L;
                                     principal = pmt;
                                     balance -= principal;
                                    } else {
                                     interest = balance * rate;
                                     principal = pmt - interest;
                                     balance -= principal;
                                    }

                                    if (k >= start) { kahan_add(interest); }
                                   }

                                   return (double)sum;
                                  }};

  // CUMPRINC (累積元金)
  cfg.functions["fin_cumprinc"] = {6, 6, [](auto &v, auto &ctx) -> Value {
                                    const long double rate = (long double)v[0].asScalar(ctx.pos);
                                    const int nper = (int)v[1].asScalar(ctx.pos);
                                    const long double pv = (long double)v[2].asScalar(ctx.pos);
                                    const int start = (int)v[3].asScalar(ctx.pos);
                                    const int end = (int)v[4].asScalar(ctx.pos);
                                    const int type = (int)v[5].asScalar(ctx.pos);

                                    if (rate < 0.0L || nper <= 0) throwDomain(ctx.pos);
                                    if (start < 1 || end < start || end > nper) throwDomain(ctx.pos);
                                    if (!(type == 0 || type == 1)) throwDomain(ctx.pos);

                                    if (nearly_zero_ld(rate)) { return (double)(pv * (long double)(end - start + 1) / (long double)nper); }

                                    const long double pmt = fin_pmt_ld(rate, nper, pv, type); // 正の支払額
                                    long double balance = pv;

                                    long double sum = 0.0L;
                                    long double c = 0.0L;

                                    auto kahan_add = [&](long double x) {
                                     long double y = x - c;
                                     long double t = sum + y;
                                     c = (t - sum) - y;
                                     sum = t;
                                    };

                                    for (int k = 1; k <= end; ++k) {
                                     long double interest = 0.0L;
                                     long double principal = 0.0L;

                                     if (type == 1 && k == 1) {
                                      interest = 0.0L;
                                      principal = pmt;
                                      balance -= principal;
                                     } else {
                                      interest = balance * rate;
                                      principal = pmt - interest;
                                      balance -= principal;
                                     }

                                     if (k >= start) { kahan_add(principal); }
                                    }

                                    return (double)sum;
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

} // namespace mm::cal
