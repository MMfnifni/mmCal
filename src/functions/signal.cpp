#include "functions/functions.hpp"

namespace mm::cal {
 void registerSignal(SystemConfig &cfg) {

  cfg.functions["sinc"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = requireReal(v[0], ctx.pos);
                            double r = angleIn(x, ctx.session.cfg);
                            if (std::abs(r) < cnst_precision_inv) return 1.0;
                            return std::sin(r) / r;
                           }};
  cfg.functions["cosc"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = requireReal(v[0], ctx.pos);
                            double r = angleIn(x, ctx.session.cfg);
                            if (std::abs(r) < cnst_precision_inv) return 0.0;
                            return (1.0 - std::cos(r)) / r;
                           }};
  cfg.functions["tanc"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            double x = requireReal(v[0], ctx.pos);
                            double r = angleIn(x, ctx.session.cfg);
                            if (std::abs(r) < cnst_precision_inv) return 1.0;
                            double c = std::cos(r);
                            return (std::abs(c) < cnst_precision_inv) ? signedInfBy(std::sin(r)) : std::tan(r) / r;
                           }};
  cfg.functions["sinhc"] = {1, 1, makeDivXCmplxReal([](Complex x) { return std::sinh(x); }, [](double x) { return std::sinh(x); }, 1.0)};
  cfg.functions["tanhc"] = {1, 1, makeDivXCmplxReal([](Complex x) { return std::tanh(x); }, [](double x) { return std::tanh(x); }, 1.0)};
  cfg.functions["expc"] = {1, 1, makeExpc()};
  cfg.functions["fft"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: fft requires vector/matrix", ctx.pos);
                           const auto &mv = v[0].asMultiRef(ctx.pos);
                           if (mv.empty()) return v[0];
                           bool is2D = true;
                           for (const auto &row : mv)
                            if (!row.isMulti()) {
                             is2D = false;
                             break;
                            }

                           if (is2D) return fft2D(v[0], false, ctx);

                           auto x = toComplexVector1D(v[0], ctx);
                           return fromComplexVector(fft_dispatch(x, false, ctx));
                          }};
  cfg.functions["ifft"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: ifft requires vector/matrix", ctx.pos);

                            const auto &mv = v[0].asMultiRef(ctx.pos);

                            if (!mv.empty() && mv[0].isMulti()) return fft2D(v[0], true, ctx);

                            auto x = toComplexVector1D(v[0], ctx);
                            return fromComplexVector(fft_dispatch(x, true, ctx));
                           }};
  cfg.functions["dft"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           auto x = toComplexVector1D(v[0], ctx);
                           return fromComplexVector(dft_impl(x, false));
                          }};
  cfg.functions["hilbert"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                               auto x = toComplexVector1D(v[0], ctx);
                               size_t N = x.size();

                               auto X = fft_dispatch(x, false, ctx);

                               std::vector<double> h(N, 0.0);

                               if (N % 2 == 0) {
                                h[0] = 1;
                                h[N / 2] = 1;
                                for (size_t i = 1; i < N / 2; i++)
                                 h[i] = 2;
                               } else {
                                h[0] = 1;
                                for (size_t i = 1; i <= (N - 1) / 2; i++)
                                 h[i] = 2;
                               }

                               for (size_t i = 0; i < N; i++)
                                X[i] *= h[i];

                               X = fft_dispatch(X, true, ctx);

                               return fromComplexVector(X);
                              }};
 }

} // namespace mm::cal
