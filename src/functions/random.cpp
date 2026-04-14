#include "functions/functions.hpp"

namespace mm::cal {

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
                              if (v.empty()) throw CalcError(CalcErrorType::DomainError, "DomainError: choice: no arguments", ctx.pos);
                              std::uniform_int_distribution<size_t> dist(0, v.size() - 1);
                              return v[dist(rng)];
                             }};

  // fma(a, b, c)
  cfg.functions["fma"] = {3, 3, [](auto &v, auto &ctx) -> Value { return std::fma(v[0].asScalar(ctx.pos), v[1].asScalar(ctx.pos), v[2].asScalar(ctx.pos)); }};
 }

} // namespace mm::cal
