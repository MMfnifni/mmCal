#include "functions/functions.hpp"

namespace mm::cal {

 static double evalRealFunction(const std::string &name, double x, FunctionContext &ctx) {
  Value y = evaluateFunction(name, {Value(x)}, ctx);
  if (!y.isScalar()) { throw CalcError(CalcErrorType::TypeError, "TypeError: calculus requires a real-valued scalar function", ctx.pos); }

  double r = y.asScalar(ctx.pos);
  if (!std::isfinite(r)) throw CalcError(CalcErrorType::DomainError, "DomainError: function evaluation produced non-finite value", ctx.pos);
  return r;
 }

 void registerCalculus(SystemConfig &cfg) {
  cfg.functions["diff"] = {2, 3, [](auto &v, auto &ctx) -> Value {
                            const std::string name = requireFunctionName(v[0], ctx.pos);
                            const double x = requireReal(v[1], ctx.pos);
                            double h = (v.size() >= 3) ? requireReal(v[2], ctx.pos) : defaultDiffStep(x);

                            if (!std::isfinite(x)) throw CalcError(CalcErrorType::DomainError, "DomainError: diff: x must be finite", ctx.pos);
                            if (!std::isfinite(h) || h == 0.0) throw CalcError(CalcErrorType::DomainError, "DomainError: diff: h must be finite and nonzero", ctx.pos);
                            h = std::abs(h);

                            const double fp = evalRealFunction(name, x + h, ctx);
                            const double fm = evalRealFunction(name, x - h, ctx);
                            const double d = (fp - fm) / (2.0 * h);

                            if (!std::isfinite(d)) throw CalcError(CalcErrorType::DomainError, "DomainError: diff: non-finite derivative estimate", ctx.pos);
                            return d;
                           }};

  cfg.functions["diff2"] = {2, 3, [](auto &v, auto &ctx) -> Value {
                             const std::string name = requireFunctionName(v[0], ctx.pos);
                             const double x = requireReal(v[1], ctx.pos);
                             double h = (v.size() >= 3) ? requireReal(v[2], ctx.pos) : defaultDiffStep2(x);

                             if (!std::isfinite(x)) throw CalcError(CalcErrorType::DomainError, "DomainError: diff2: x must be finite", ctx.pos);
                             if (!std::isfinite(h) || h == 0.0) throw CalcError(CalcErrorType::DomainError, "DomainError: diff2: h must be finite and nonzero", ctx.pos);
                             h = std::abs(h);

                             const double fp = evalRealFunction(name, x + h, ctx);
                             const double f0 = evalRealFunction(name, x, ctx);
                             const double fm = evalRealFunction(name, x - h, ctx);
                             const double d2 = (fp - 2.0 * f0 + fm) / (h * h);

                             if (!std::isfinite(d2)) throw CalcError(CalcErrorType::DomainError, "DomainError: diff2: non-finite derivative estimate", ctx.pos);
                             return d2;
                            }};

  cfg.functions["simpson"] = {3, 4, [](auto &v, auto &ctx) -> Value {
                               const std::string name = requireFunctionName(v[0], ctx.pos);
                               double a = requireReal(v[1], ctx.pos);
                               double b = requireReal(v[2], ctx.pos);
                               long long n = (v.size() >= 4) ? requireInt(v[3], ctx.pos) : 100;

                               if (!std::isfinite(a) || !std::isfinite(b)) throw CalcError(CalcErrorType::DomainError, "DomainError: simpson: bounds must be finite", ctx.pos);
                               if (n <= 0) throw CalcError(CalcErrorType::DomainError, "DomainError: simpson: n must be positive", ctx.pos);
                               if ((n % 2) != 0) throw CalcError(CalcErrorType::DomainError, "DomainError: simpson: n must be even", ctx.pos);
                               if (a == b) return 0.0;

                               double sign = 1.0;
                               if (a > b) {
                                std::swap(a, b);
                                sign = -1.0;
                               }

                               const double h = (b - a) / static_cast<double>(n);
                               long double sum = evalRealFunction(name, a, ctx) + evalRealFunction(name, b, ctx);

                               for (long long i = 1; i < n; ++i) {
                                const double x = a + h * static_cast<double>(i);
                                const double fx = evalRealFunction(name, x, ctx);
                                sum += ((i % 2) == 0 ? 2.0L : 4.0L) * static_cast<long double>(fx);
                               }

                               const long double result = static_cast<long double>(sign) * static_cast<long double>(h) * sum / 3.0L;
                               if (!std::isfinite(static_cast<double>(result))) throw CalcError(CalcErrorType::Overflow, "Overflow: simpson overflow", ctx.pos);
                               return static_cast<double>(result);
                              }};

  cfg.functions["trapz"] = {3, 4, [](auto &v, auto &ctx) -> Value {
                             const std::string name = requireFunctionName(v[0], ctx.pos);
                             double a = requireReal(v[1], ctx.pos);
                             double b = requireReal(v[2], ctx.pos);
                             long long n = (v.size() >= 4) ? requireInt(v[3], ctx.pos) : 100;

                             if (!std::isfinite(a) || !std::isfinite(b)) throw CalcError(CalcErrorType::DomainError, "DomainError: trapz: bounds must be finite", ctx.pos);
                             if (n <= 0) throw CalcError(CalcErrorType::DomainError, "DomainError: trapz: n must be positive", ctx.pos);
                             if (a == b) return 0.0;

                             double sign = 1.0;
                             if (a > b) {
                              std::swap(a, b);
                              sign = -1.0;
                             }

                             const double h = (b - a) / static_cast<double>(n);
                             long double sum = 0.5L * static_cast<long double>(evalRealFunction(name, a, ctx) + evalRealFunction(name, b, ctx));

                             for (long long i = 1; i < n; ++i) {
                              const double x = a + h * static_cast<double>(i);
                              sum += static_cast<long double>(evalRealFunction(name, x, ctx));
                             }

                             const long double result = static_cast<long double>(sign) * static_cast<long double>(h) * sum;
                             if (!std::isfinite(static_cast<double>(result))) throw CalcError(CalcErrorType::Overflow, "Overflow: trapz overflow", ctx.pos);
                             return static_cast<double>(result);
                            }};

  cfg.functions["integrate"] = {3, 4, [](auto &v, auto &ctx) -> Value { return evaluateFunction("simpson", v, ctx); }};
 }

} // namespace mm::cal
