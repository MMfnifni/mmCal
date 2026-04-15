#include "functions/functions.hpp"

namespace mm::cal {

 void registerGeoVec(SystemConfig &cfg) {
  cfg.functions["hypot"] = {2, 2, [](auto &v, auto &ctx) -> Value { return std::hypot(requireReal(v[0], ctx.pos), requireReal(v[1], ctx.pos)); }};
  cfg.functions["norm"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                            if (v.empty()) throw CalcError(CalcErrorType::DomainError, "DomainError: no elements", ctx.pos);
                            double acc = 0.0;
                            for (const auto &x : v) {
                             const double r = x.asScalar(ctx.pos);
                             acc += r * r; // pow(r,2) より圧倒的に速い
                            }
                            return std::sqrt(acc);
                           }};

  // ベクトル和
  cfg.functions["vadd"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vadd requires two multivalue arguments", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            const auto &b = v[1].asMultiRef(ctx.pos);

                            if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);

                            auto result = std::make_shared<MultiValue>();
                            result->elems_.reserve(a.size());
                            for (size_t i = 0; i < a.size(); ++i) {
                             result->elems_.push_back(Value(a[i].asScalar(ctx.pos) + b[i].asScalar(ctx.pos)));
                            }
                            return Value(result);
                           }};

  // ベクトル差
  cfg.functions["vsub"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vsub requires two multivalue arguments", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            const auto &b = v[1].asMultiRef(ctx.pos);

                            if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);

                            auto result = std::make_shared<MultiValue>();
                            result->elems_.reserve(a.size());
                            for (size_t i = 0; i < a.size(); ++i) {
                             result->elems_.push_back(Value(a[i].asScalar(ctx.pos) - b[i].asScalar(ctx.pos)));
                            }
                            return Value(result);
                           }};

  // スカラ倍
  cfg.functions["vscalar"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                               if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vscalar requires multivalue argument", ctx.pos);
                               double scalar = v[1].asScalar(ctx.pos);
                               const auto &vec = v[0].asMultiRef(ctx.pos);

                               auto result = std::make_shared<MultiValue>();
                               result->elems_.reserve(vec.size());
                               for (const auto &elem : vec.elems()) {
                                result->elems_.push_back(Value(elem.asScalar(ctx.pos) * scalar));
                               }
                               return Value(result);
                              }};
  // 内積
  cfg.functions["vdot"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vdot requires two multivalue arguments", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            const auto &b = v[1].asMultiRef(ctx.pos);

                            if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);

                            double acc = 0.0;
                            for (size_t i = 0; i < a.size(); ++i) {
                             acc += a[i].asScalar(ctx.pos) * b[i].asScalar(ctx.pos);
                            }
                            return Value(acc);
                           }};

  // 外積 (3次元限定)
  cfg.functions["vcross"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vcross requires two multivalue arguments", ctx.pos);
                              const auto &a = v[0].asMultiRef(ctx.pos);
                              const auto &b = v[1].asMultiRef(ctx.pos);

                              if (a.size() != 3 || b.size() != 3) throw CalcError(CalcErrorType::DomainError, "DomainError: vcross requires 3-dimensional vectors", ctx.pos);

                              double x = a[1].asScalar(ctx.pos) * b[2].asScalar(ctx.pos) - a[2].asScalar(ctx.pos) * b[1].asScalar(ctx.pos);
                              double y = a[2].asScalar(ctx.pos) * b[0].asScalar(ctx.pos) - a[0].asScalar(ctx.pos) * b[2].asScalar(ctx.pos);
                              double z = a[0].asScalar(ctx.pos) * b[1].asScalar(ctx.pos) - a[1].asScalar(ctx.pos) * b[0].asScalar(ctx.pos);

                              auto result = std::make_shared<MultiValue>(std::vector<Value>{Value(x), Value(y), Value(z)});
                              return Value(result);
                             }};
  cfg.functions["vnorm"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vnorm requires multivalue argument", ctx.pos);
                             const auto &vec = v[0].asMultiRef(ctx.pos);

                             double sum = 0.0;
                             for (const auto &elem : vec.elems()) {
                              double val = elem.asScalar(ctx.pos);
                              sum += val * val;
                             }
                             return Value(std::sqrt(sum));
                            }};

  // マンハッタン距離
  cfg.functions["vmanhattan"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vmanhattan requires two multivalue arguments", ctx.pos);
                                  const auto &a = v[0].asMultiRef(ctx.pos);
                                  const auto &b = v[1].asMultiRef(ctx.pos);

                                  if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);

                                  double sum = 0.0;
                                  for (size_t i = 0; i < a.size(); ++i) {
                                   sum += std::abs(a[i].asScalar(ctx.pos) - b[i].asScalar(ctx.pos));
                                  }
                                  return Value(sum);
                                 }};

  // ユークリッド距離
  cfg.functions["veuclidean"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "veuclidean requires two multivalue arguments", ctx.pos);
                                  const auto &a = v[0].asMultiRef(ctx.pos);
                                  const auto &b = v[1].asMultiRef(ctx.pos);

                                  if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);

                                  double sum = 0.0;
                                  for (size_t i = 0; i < a.size(); ++i) {
                                   double diff = a[i].asScalar(ctx.pos) - b[i].asScalar(ctx.pos);
                                   sum += diff * diff;
                                  }
                                  return Value(std::sqrt(sum));
                                 }};

  cfg.functions["dot"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                           if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "dot requires two multivalue arguments", ctx.pos);

                           const auto &a = v[0].asMultiRef(ctx.pos);
                           const auto &b = v[1].asMultiRef(ctx.pos);
                           if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);
                           if (a.empty()) throw CalcError(CalcErrorType::DomainError, "DomainError: dot of empty multivalue", ctx.pos);
                           double acc = 0.0;
                           for (size_t i = 0; i < a.size(); ++i) { // nested Multi は禁止するよん
                            if (a[i].isMulti() || b[i].isMulti()) throw CalcError(CalcErrorType::TypeError, "nested multivalue not allowed", ctx.pos);
                            acc += a[i].asScalar(ctx.pos) * b[i].asScalar(ctx.pos);
                           }
                           return Value(acc);
                          }};
  cfg.functions["vsum"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vsum requires multivalue argument", ctx.pos);
                            const auto &a = v[0].asMultiRef(ctx.pos);
                            double sum = 0.0;
                            for (const auto &elem : a.elems()) {
                             sum += elem.asScalar(ctx.pos);
                            }
                            return Value(sum);
                           }};
  // ベクトルの正規化
  cfg.functions["vnormalize"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vnormalize requires multivalue argument", ctx.pos);
                                  const auto &vec = v[0].asMultiRef(ctx.pos);

                                  double norm = 0.0;
                                  for (const auto &elem : vec.elems()) {
                                   double val = elem.asScalar(ctx.pos);
                                   norm += val * val;
                                  }
                                  norm = std::sqrt(norm);

                                  if (norm == 0.0) throw CalcError(CalcErrorType::DomainError, "DomainError: cannot normalize zero vector", ctx.pos);

                                  auto result = std::make_shared<MultiValue>();
                                  result->elems_.reserve(vec.size());
                                  for (const auto &elem : vec.elems()) {
                                   result->elems_.push_back(Value(elem.asScalar(ctx.pos) / norm));
                                  }
                                  return Value(result);
                                 }};

  // ベクトルの射影
  cfg.functions["vproject"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vproject requires two multivalue arguments", ctx.pos);
                                const auto &a = v[0].asMultiRef(ctx.pos);
                                const auto &b = v[1].asMultiRef(ctx.pos);

                                if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);

                                // 内積計算
                                double dot_product = 0.0;
                                for (size_t i = 0; i < a.size(); ++i) {
                                 dot_product += a[i].asScalar(ctx.pos) * b[i].asScalar(ctx.pos);
                                }

                                // bのノルムの2乗
                                double b_norm_sq = 0.0;
                                for (size_t i = 0; i < b.size(); ++i) {
                                 double val = b[i].asScalar(ctx.pos);
                                 b_norm_sq += val * val;
                                }

                                if (std::abs(b_norm_sq) < 1e-15) throw CalcError(CalcErrorType::DomainError, "DomainError: cannot project onto zero vector", ctx.pos);

                                double scalar = dot_product / b_norm_sq;

                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(b.size());
                                for (size_t i = 0; i < b.size(); ++i) {
                                 result->elems_.push_back(Value(b[i].asScalar(ctx.pos) * scalar));
                                }
                                return Value(result);
                               }};

  // 2ベクトル間の角度 (度数法)
  cfg.functions["vangle"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "vangle requires two multivalue arguments", ctx.pos);
                              const auto &a = v[0].asMultiRef(ctx.pos);
                              const auto &b = v[1].asMultiRef(ctx.pos);

                              if (a.size() != b.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);

                              // 内積計算
                              double dot_product = 0.0;
                              for (size_t i = 0; i < a.size(); ++i) {
                               dot_product += a[i].asScalar(ctx.pos) * b[i].asScalar(ctx.pos);
                              }

                              // norm計算
                              double norm_a = 0.0, norm_b = 0.0;
                              for (size_t i = 0; i < a.size(); ++i) {
                               double val_a = a[i].asScalar(ctx.pos);
                               double val_b = b[i].asScalar(ctx.pos);
                               norm_a += val_a * val_a;
                               norm_b += val_b * val_b;
                              }
                              norm_a = std::sqrt(norm_a);
                              norm_b = std::sqrt(norm_b);

                              if (norm_a == 0.0 || norm_b == 0.0) throw CalcError(CalcErrorType::DomainError, "DomainError: cannot calculate angle with zero vector", ctx.pos);

                              // cos計算
                              double cos_angle = dot_product / (norm_a * norm_b);
                              // cosが-1以上1以下になるようにクランプ
                              cos_angle = std::max(-1.0, std::min(1.0, cos_angle));

                              double angle_rad = std::acos(cos_angle);
                              return angleOut(angle_rad, ctx.session.cfg);
                             }};

  // ベクトルの反射
  cfg.functions["vreflect"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: vreflect requires two multivalue arguments", ctx.pos);

                                const auto &a = v[0].asMultiRef(ctx.pos);
                                const auto &n = v[1].asMultiRef(ctx.pos);

                                if (a.size() != n.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);

                                double dot = 0.0;
                                double nn = 0.0;

                                for (size_t i = 0; i < a.size(); ++i) {
                                 double ai = a[i].asScalar(ctx.pos);
                                 double ni = n[i].asScalar(ctx.pos);
                                 dot += ai * ni;
                                 nn += ni * ni;
                                }

                                if (std::abs(nn) < 1e-15) throw CalcError(CalcErrorType::DomainError, "DomainError: cannot reflect with zero normal vector", ctx.pos);

                                double factor = 2.0 * dot / nn;

                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(a.size());

                                for (size_t i = 0; i < a.size(); ++i) {
                                 double ai = a[i].asScalar(ctx.pos);
                                 double ni = n[i].asScalar(ctx.pos);
                                 result->elems_.push_back(Value(ai - factor * ni));
                                }

                                return Value(result);
                               }};
  // 軸対称反射
  cfg.functions["vreflect_axis"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                     if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: vreflect_axis requires two multivalue arguments", ctx.pos);

                                     const auto &a = v[0].asMultiRef(ctx.pos);
                                     const auto &n = v[1].asMultiRef(ctx.pos);

                                     if (a.size() != n.size()) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);

                                     double dot = 0.0;
                                     double nn = 0.0;

                                     for (size_t i = 0; i < a.size(); ++i) {
                                      double ai = a[i].asScalar(ctx.pos);
                                      double ni = n[i].asScalar(ctx.pos);
                                      dot += ai * ni;
                                      nn += ni * ni;
                                     }

                                     if (std::abs(nn) < 1e-15) throw CalcError(CalcErrorType::DomainError, "DomainError: cannot reflect with zero axis vector", ctx.pos);

                                     double factor = 2.0 * dot / nn;

                                     auto result = std::make_shared<MultiValue>();
                                     result->elems_.reserve(a.size());

                                     for (size_t i = 0; i < a.size(); ++i) {
                                      double ai = a[i].asScalar(ctx.pos);
                                      double ni = n[i].asScalar(ctx.pos);
                                      result->elems_.push_back(Value(factor * ni - ai));
                                     }

                                     return Value(result);
                                    }};

  cfg.functions["vlength"] = cfg.functions["vnorm"];
  cfg.functions["vdistance"] = cfg.functions["veuclidean"];

  // 単位ベクトル
  cfg.functions["vunit"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "vunit requires multivalue argument", ctx.pos);
                             const auto &vec = v[0].asMultiRef(ctx.pos);

                             double norm = 0.0;
                             for (const auto &elem : vec.elems()) {
                              double val = elem.asScalar(ctx.pos);
                              norm += val * val;
                             }
                             norm = std::sqrt(norm);

                             if (norm == 0.0) throw CalcError(CalcErrorType::DomainError, "DomainError: cannot create unit vector from zero vector", ctx.pos);

                             auto result = std::make_shared<MultiValue>();
                             result->elems_.reserve(vec.size());
                             for (const auto &elem : vec.elems()) {
                              result->elems_.push_back(Value(elem.asScalar(ctx.pos) / norm));
                             }
                             return Value(result);
                            }};

  cfg.functions["scalar"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "scalar requires multivalue argument", ctx.pos);
                              double scalar = v[1].asScalar(ctx.pos);
                              const auto &a = v[0].asMultiRef(ctx.pos);

                              auto result = std::make_shared<MultiValue>();
                              result->elems_.reserve(a.size());
                              for (const auto &elem : a.elems()) {
                               result->elems_.push_back(Value(elem.asScalar(ctx.pos) * scalar));
                              }
                              return Value(result);
                             }};

  cfg.functions["lerp"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                            double a = v[0].asScalar(ctx.pos);
                            double b = v[1].asScalar(ctx.pos);
                            double t = v[2].asScalar(ctx.pos);
                            return a * (1 - t) + b * t;
                            // return std::lerp(a, b, t); // leapは丸めが微妙なことがあった
                           }};
  cfg.functions["distance"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                if (v.size() % 2 != 0) throw CalcError(CalcErrorType::DomainError, "DomainError: dimension mismatch", ctx.pos);
                                size_t n = v.size() / 2;
                                if (n == 0) throw CalcError(CalcErrorType::DomainError, "DomainError: no elements", ctx.pos);
                                double acc = 0.0;
                                // for (size_t i = 0; i < n; ++i) {
                                //  acc += std::pow(asReal(v[i + n], ctx.pos) - v[i].asScalar(ctx.pos), 2);
                                // }
                                for (size_t i = 0; i < n; ++i) {
                                 const double d = v[i + n].asScalar(ctx.pos) - v[i].asScalar(ctx.pos);
                                 acc += d * d; // pow(d,2) より圧倒的に速い
                                }
                                return std::sqrt(acc);
                               }};
  cfg.functions["totient"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                               int n = (int)v[0].asScalar(ctx.pos);
                               if (n < 1) throwDomain(ctx.pos, "must be >= 1");
                               int phi = n;
                               for (int p = 2; p * p <= n; ++p) {
                                if (n % p == 0) {
                                 while (n % p == 0)
                                  n /= p;
                                 phi -= phi / p;
                                }
                               }
                               if (n > 1) phi -= phi / n;
                               return (double)phi;
                              }};
  cfg.functions["covmatrix"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                                 if (v.empty()) throwDomain(ctx.pos, "need at least 1 vector");

                                 std::vector<std::vector<double>> vars;
                                 vars.reserve(v.size());

                                 size_t n = 0;
                                 for (const auto &arg : v) {
                                  vars.push_back(requireRealVectorArg(arg, ctx.pos));
                                  if (vars.back().size() < 2) throwDomain(ctx.pos, "each vector must have at least 2 samples");
                                  if (n == 0) n = vars.back().size();
                                  else if (vars.back().size() != n) throwDomain(ctx.pos, "sample length mismatch");
                                 }

                                 const size_t p = vars.size();
                                 std::vector<double> mean(p, 0.0);
                                 for (size_t i = 0; i < p; ++i) {
                                  long double acc = 0.0;
                                  for (double x : vars[i])
                                   acc += x;
                                  mean[i] = static_cast<double>(acc / n);
                                 }

                                 std::vector<std::vector<double>> C(p, std::vector<double>(p, 0.0));
                                 for (size_t i = 0; i < p; ++i) {
                                  for (size_t j = i; j < p; ++j) {
                                   long double acc = 0.0;
                                   for (size_t k = 0; k < n; ++k) {
                                    acc += (long double)(vars[i][k] - mean[i]) * (vars[j][k] - mean[j]);
                                   }
                                   double cij = static_cast<double>(acc / (n - 1));
                                   C[i][j] = cij;
                                   C[j][i] = cij;
                                  }
                                 }

                                 return makeRowMajorMatrixValue(C);
                                }};
  cfg.functions["corrmatrix"] = {1, -1, [](auto &v, auto &ctx) -> Value {
                                  if (v.empty()) throwDomain(ctx.pos, "need at least 1 vector");

                                  std::vector<std::vector<double>> vars;
                                  vars.reserve(v.size());

                                  size_t n = 0;
                                  for (const auto &arg : v) {
                                   vars.push_back(requireRealVectorArg(arg, ctx.pos));
                                   if (vars.back().size() < 2) throwDomain(ctx.pos, "each vector must have at least 2 samples");
                                   if (n == 0) n = vars.back().size();
                                   else if (vars.back().size() != n) throwDomain(ctx.pos, "sample length mismatch");
                                  }

                                  const size_t p = vars.size();
                                  std::vector<double> mean(p, 0.0), sd(p, 0.0);
                                  for (size_t i = 0; i < p; ++i) {
                                   long double acc = 0.0;
                                   for (double x : vars[i])
                                    acc += x;
                                   mean[i] = static_cast<double>(acc / n);

                                   long double ss = 0.0;
                                   for (double x : vars[i]) {
                                    long double d = x - mean[i];
                                    ss += d * d;
                                   }
                                   sd[i] = std::sqrt(static_cast<double>(ss / (n - 1)));
                                   if (sd[i] == 0.0) throwDomain(ctx.pos, "zero variance");
                                  }

                                  std::vector<std::vector<double>> R(p, std::vector<double>(p, 0.0));
                                  for (size_t i = 0; i < p; ++i) {
                                   R[i][i] = 1.0;
                                   for (size_t j = i + 1; j < p; ++j) {
                                    long double acc = 0.0;
                                    for (size_t k = 0; k < n; ++k) {
                                     acc += (long double)(vars[i][k] - mean[i]) * (vars[j][k] - mean[j]);
                                    }
                                    double rij = static_cast<double>(acc / ((n - 1) * sd[i] * sd[j]));
                                    R[i][j] = rij;
                                    R[j][i] = rij;
                                   }
                                  }

                                  return makeRowMajorMatrixValue(R);
                                 }};
  cfg.functions["percentrank"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                                   if (v.size() < 2) throwDomain(ctx.pos);

                                   const double val = requireReal(v[0], ctx.pos);

                                   std::vector<double> x;
                                   for (size_t i = 1; i < v.size(); ++i)
                                    x.push_back(requireReal(v[i], ctx.pos));

                                   if (x.empty()) throwDomain(ctx.pos);

                                   size_t count = 0;
                                   for (double d : x)
                                    if (d <= val) ++count;

                                   return (double)count / x.size() * 100.0;
                                  }};

  cfg.functions["winsorR"] = {2, -1, [](auto &v, auto &ctx) -> Value {
                               if (v.size() < 3) throwDomain(ctx.pos);

                               const double p = requireReal(v[0], ctx.pos);
                               if (!(p >= 0.0 && p <= 0.5)) throwDomain(ctx.pos);

                               std::vector<double> x;
                               for (size_t i = 1; i < v.size(); ++i)
                                x.push_back(requireReal(v[i], ctx.pos));

                               if (x.empty()) throwDomain(ctx.pos);

                               std::vector<double> y = x;
                               std::sort(y.begin(), y.end());

                               size_t k = (size_t)(p * y.size());
                               double low = y[k];
                               double high = y[y.size() - 1 - k];

                               long double sum = 0;
                               for (double d : x)
                                sum += std::min(std::max(d, low), high);

                               return (double)(sum / x.size());
                              }};

  cfg.functions["convolve"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                                auto x = requireRealVectorArg(v[0], ctx.pos);
                                auto y = requireRealVectorArg(v[1], ctx.pos);
                                if (x.empty() || y.empty()) throwDomain(ctx.pos, "empty vector");

                                std::vector<long double> acc(x.size() + y.size() - 1, 0.0L);
                                for (size_t i = 0; i < x.size(); ++i) {
                                 for (size_t j = 0; j < y.size(); ++j) {
                                  acc[i + j] += (long double)x[i] * y[j];
                                 }
                                }

                                auto out = std::make_shared<MultiValue>();
                                out->elems_.reserve(acc.size());
                                for (auto z : acc)
                                 out->elems_.push_back(Value((double)z));
                                return Value(out);
                               }};
 }

} // namespace mm::cal
