#include "math_util.hpp"

namespace mm::cal {

 /* ============================
   数学補助系
   ============================ */

 /* ============================
   整数系
   ============================ */

 long long combLL(long long n, long long r, size_t pos) {
  if (r < 0) throwDomain(pos);
  if (r > n) return 0;

  r = std::min(r, n - r);
  long long res = 1;

  for (long long i = 1; i <= r; ++i) {
   long long num = n - r + i;
   long long den = i;

   const auto g1 = gcdLL(num, den);
   num /= g1;
   den /= g1;

   const auto g2 = gcdLL(res, den);
   res /= g2;
   den /= g2;

   // den は必ず 1 になる（整数の二項係数なので）
   res = checkedMul(res, num, pos);
  }
  return res;
 }

 long long permLL(long long n, long long r, size_t pos) {
  if (n < 0 || r < 0) throwDomain(pos);
  if (r > n) return 0;

  long long res = 1;
  for (long long k = 0; k < r; ++k)
   res = checkedMul(res, n - k, pos);

  return res;
 }

 double factorial(double x, size_t pos) {
  if (!std::isfinite(x)) throwDomain(pos);
  if (x < 0) throwDomain(pos);
  if (std::floor(x) != x) throwDomain(pos);
  if (x > 170) // 171! は double overflow
   throwOverflow(pos);

  double r = 1.0;
  for (int i = 2; i <= (int)x; ++i)
   r *= i;
  return r;
 }

 long double factLD(long long n, size_t pos) {
  if (n < 0) throw CalcError(CalcErrorType::DomainError, "fact: n < 0", pos);

  long double acc = 1.0L;
  for (long long i = 2; i <= n; ++i) {
   acc *= (long double)i;
   if (!std::isfinite(acc)) throwOverflow(pos);
  }
  return acc;
 }

 unsigned long long fibULL(unsigned long long n, size_t pos) {

  auto add = [&](unsigned long long x, unsigned long long y) -> unsigned long long {
   if (ULLONG_MAX - x < y) throwOverflow(pos);
   return x + y;
  };

  auto sub = [&](unsigned long long x, unsigned long long y) -> unsigned long long {
   // 本来 fib の内部では起きないが、式を安全にするため入れておく
   if (x < y) throw CalcError(CalcErrorType::Overflow, "fib: underflow", pos);
   return x - y;
  };

  auto mul = [&](unsigned long long x, unsigned long long y) -> unsigned long long {
   if (x && ULLONG_MAX / x < y) throwOverflow(pos);
   return x * y;
  };

  // 戻り値: {F(k), F(k+1)}
  auto rec = [&](auto &&self, unsigned long long k) -> std::pair<unsigned long long, unsigned long long> {
   if (k == 0) return {0, 1};

   auto [a, b] = self(self, k >> 1); // a=F(m), b=F(m+1)

   // c = F(2m)   = a * (2b - a)
   // d = F(2m+1) = a^2 + b^2
   const auto two_b_minus_a = add(b, sub(b, a)); // = 2b - a を安全に
   const auto c = mul(a, two_b_minus_a);
   const auto d = add(mul(a, a), mul(b, b));

   if ((k & 1) == 0) return {c, d}; // {F(2m), F(2m+1)}
   return {d, add(c, d)};           // {F(2m+1), F(2m+2)}
  };

  return rec(rec, (unsigned long long)n).first;
 }

 // static unsigned long long fibULL(long long n, int pos) {
 //  if (n < 0) throw CalcError(CalcErrorType::DomainError, "fib: n < 0", pos);
 // unsigned long long a = 0, b = 1;
 // while (n--) {
 //  if (ULLONG_MAX - a < b) throwOverflow;
 //  const auto c = a + b;
 //  a = b;
 //  b = c;
 // }
 // return a;
 //}

 Value area_polygon_impl(const std::vector<Value> &v, FunctionContext &ctx) {
  size_t n = v.size() / 2;

  if (n < 3) throw CalcError(CalcErrorType::DomainError, "area_polygon: need at least 3 points", ctx.pos);

  double area = 0.0;

  for (size_t i = 0; i < n; ++i) {
   double x0 = v[2 * i].asScalar(ctx.pos);
   double y0 = v[2 * i + 1].asScalar(ctx.pos);

   double x1 = v[2 * ((i + 1) % n)].asScalar(ctx.pos);
   double y1 = v[2 * ((i + 1) % n) + 1].asScalar(ctx.pos);

   area += x0 * y1 - x1 * y0;
  }

  return std::abs(area) / 2.0;
 }

 /* ============================
   統計系
   ============================ */
 std::vector<double> gatherReals(const std::vector<Value> &v, size_t pos) {
  std::vector<double> a;
  a.reserve(v.size());
  for (const auto &x : v) {
   if (isComplex(x)) {
    const auto z = x.asComplex(pos);
    if (std::abs(std::imag(z)) > cnst_precision_inv) throwDomain(pos);
    a.push_back(std::real(z));
   } else {
    a.push_back(x.asScalar(pos));
   }
  }
  return a;
 }

 std::vector<double> collectReals(const std::vector<Value> &v, FunctionContext &ctx) {
  std::vector<double> a;
  a.reserve(v.size());
  for (const auto &x : v)
   a.push_back(x.asScalar(ctx.pos));
  return a;
 }
 std::vector<Complex> collectComplex(const std::vector<Value> &v, FunctionContext &ctx) {
  std::vector<Complex> a;
  a.reserve(v.size());
  for (const auto &x : v)
   a.push_back(x.asComplex(ctx.pos));
  return a;
 }

 double meanOf(const std::vector<double> &a) {
  const long double s = std::accumulate(a.begin(), a.end(), 0.0L);
  return (double)(s / (long double)a.size());
 }

 double stddevPopulation(const std::vector<double> &a, double mu) { return std::sqrt(variancePopulation(a, mu)); }

 double quantileLinear(std::vector<double> a, double p, size_t pos) {
  if (!(0.0 <= p && p <= 1.0) || a.empty()) throwDomain(pos);
  std::sort(a.begin(), a.end());
  if (a.size() == 1) return a[0];
  const double x = p * (double)(a.size() - 1);
  const auto i = (size_t)x;
  const auto j = std::min(i + 1, a.size() - 1);
  const double t = x - (double)i;
  return std::lerp(a[i], a[j], t);
 }

 std::vector<double> rankAverageTies(const std::vector<double> &x, size_t pos) {
  const size_t n = x.size();
  std::vector<double> r(n);
  if (n == 0) return r;

  // NaN は統計として扱えない（sortの前提も壊れる）
  for (double v : x) {
   if (std::isnan(v)) throwDomain(pos);
  }

  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);

  // strict weak ordering を守る比較
  std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });

  // ties 判定（相対+絶対）
  auto eq = [](double a, double b) {
   const double d = std::abs(a - b);
   const double s = std::max({1.0, std::abs(a), std::abs(b)});
   return d <= cnst_precision_inv * s;
  };

  for (size_t i = 0; i < n;) {
   size_t j = i;
   while (j + 1 < n && eq(x[idx[j + 1]], x[idx[i]]))
    ++j;

   // rank は 1-based
   const double avg = 0.5 * ((double)(i + 1) + (double)(j + 1));
   for (size_t k = i; k <= j; ++k)
    r[idx[k]] = avg;

   i = j + 1;
  }

  return r;
 }

 double pearsonCorr(const std::vector<double> &x, const std::vector<double> &y, size_t pos) {
  if (x.size() != y.size() || x.empty()) throwDomain(pos);
  long double mx = 0.0L, my = 0.0L;
  long double sxx = 0.0L, syy = 0.0L, sxy = 0.0L;
  for (size_t i = 0; i < x.size(); ++i) {
   const long double xi = x[i];
   const long double yi = y[i];
   const long double dx = xi - mx;
   const long double dy = yi - my;
   mx += dx / (long double)(i + 1);
   my += dy / (long double)(i + 1);
   sxx += dx * (xi - mx);
   syy += dy * (yi - my);
   sxy += dx * (yi - my);
  }
  if (sxx == 0.0L || syy == 0.0L) throwDomain(pos);
  return (double)(sxy / std::sqrt(sxx * syy));
 }

 double medianOfSorted(const std::vector<double> &a) {
  const size_t n = a.size();
  return n == 0 ? 0.0 : (n & 1) ? a[n / 2] : 0.5 * (a[n / 2 - 1] + a[n / 2]);
 }

 double variancePopulation(const std::vector<double> &a, double mu) {
  long double ss = 0;
  for (double x : a) {
   long double d = (long double)x - (long double)mu;
   ss += d * d;
  }
  return (double)(ss / (long double)a.size());
 }

 double zetaEulerMaclaurin(double s) {
  const int N = 12; // 打ち切り位置
  const int M = 4;  // Bernoulli 項数

  // Bernoulli numbers B2, B4, B6, B8
  static const double B[] = {1.0 / 6.0, -1.0 / 30.0, 1.0 / 42.0, -1.0 / 30.0};

  double sum = 0.0;
  for (int n = 1; n < N; ++n)
   sum += 1.0 / std::pow(n, s);

  double Ns = std::pow(N, -s);
  sum += std::pow(N, 1.0 - s) / (s - 1.0);
  sum += 0.5 * Ns;

  for (int k = 1; k <= M; ++k) {
   double term = B[k - 1] * std::tgamma(s + 2 * k - 1) / std::tgamma(s) * std::pow(N, -s - 2 * k + 1);
   sum += term;
  }

  return sum;
 }

 double brent(std::function<double(double)> f, double a, double b, FunctionContext &ctx, double tol, int maxIter) {
  double fa = f(a), fb = f(b);
  if (fa * fb > 0) throw CalcError(CalcErrorType::NonConvergence, "Brent method did not converge (root not bracketed)", ctx.pos);
  if (std::abs(fa) < std::abs(fb)) {
   std::swap(a, b);
   std::swap(fa, fb);
  }

  double c = a, fc = fa, d = 0;
  bool mflag = true;
  for (int iter = 0; iter < maxIter; ++iter) {
   double s;
   if (fa != fc && fb != fc) { // インターポレーション（逆二次）
    s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb));
   } else { // 線形補間
    s = b - fb * (b - a) / (fb - fa);
   }

   double midpoint = (3 * a + b) / 4;
   bool bisect = !((s > midpoint && s < b) || (s < midpoint && s > b)) || (mflag ? std::abs(s - b) >= std::abs(b - c) / 2 : std::abs(s - b) >= std::abs(c - d) / 2) || (mflag ? std::abs(b - c) < tol : std::abs(c - d) < tol);

   if (bisect) s = (a + b) / 2, mflag = true;
   else mflag = false;

   double fs = f(s);
   d = c;
   c = b;
   fc = fb;

   if (fa * fs < 0) b = s, fb = fs;
   else a = s, fa = fs;

   if (std::abs(fa) < std::abs(fb)) {
    std::swap(a, b);
    std::swap(fa, fb);
   }

   if (std::abs(b - a) < tol) return b;
  }

  throw CalcError(CalcErrorType::NonConvergence, "Brent method did not converge", ctx.pos);
 }

 std::pair<double, double> fullScanBracket(std::function<double(double)> f, double start, double end, double step, FunctionContext &ctx) {
  double prev = f(start);
  for (double x = start + step; x <= end; x += step) {
   double curr = f(x);
   if (prev * curr <= 0) return {x - step, x};
   prev = curr;
  }
  throw CalcError(CalcErrorType::NonConvergence, "IRR root not bracketed in full scan", ctx.pos);
 }

 Value safeInv(double denom, double sign) {
  if (std::abs(denom) < cnst_precision_inv) return std::numeric_limits<double>::infinity() * (sign >= 0 ? 1.0 : -1.0);
  return 1.0 / denom;
 }

 bool eq(double a, double b) { return std::abs(a - b) < cnst_precision_inv; }

 inline Value areaPolygon(const std::vector<Value> &v, FunctionContext &ctx) {
  size_t n = v.size() / 2;

  if (n < 3) throw CalcError(CalcErrorType::DomainError, "area_polygon: need at least 3 points", ctx.pos);

  double area = 0.0;

  for (size_t i = 0; i < n; ++i) {

   double x0 = requireReal(v[2 * i], ctx.pos);
   double y0 = requireReal(v[2 * i + 1], ctx.pos);

   double x1 = requireReal(v[2 * ((i + 1) % n)], ctx.pos);

   double y1 = requireReal(v[2 * ((i + 1) % n) + 1], ctx.pos);

   area += x0 * y1 - x1 * y0;
  }

  return std::abs(area) / 2.0;
 }
 inline std::vector<double> collectNumericVector(const std::vector<Value> &v, FunctionContext &ctx) {
  std::vector<Value> tmp = v;
  return collectReals(tmp, ctx);
 }

 // 固有値を求める関数（QR反復法）
 std::vector<double> compute_eigenvalues(std::vector<std::vector<double>> &matrix, const FunctionContext &ctx) {
  size_t n = matrix.size();
  std::vector<double> eigenvals(n);

  // 対角成分を初期化
  for (size_t i = 0; i < n; ++i) {
   eigenvals[i] = matrix[i][i];
  }

  // QR反復法の実装
  const int max_iterations = 1000;
  const double epsilon = 1e-12;

  for (int iter = 0; iter < max_iterations; ++iter) {
   bool converged = true;

   // QR分解を計算
   std::vector<std::vector<double>> Q(n, std::vector<double>(n, 0.0));
   std::vector<std::vector<double>> R = matrix;

   // 単位行列で初期化
   for (size_t i = 0; i < n; ++i) {
    Q[i][i] = 1.0;
   }

   // ハウスホルダー変換でQR分解
   for (size_t k = 0; k < n - 1; ++k) {
    // ベクトルv = R[k:n][k] を取得
    std::vector<double> v(n - k);
    for (size_t i = 0; i < n - k; ++i) {
     v[i] = R[k + i][k];
    }

    // ノルムを計算
    double norm = 0.0;
    for (double val : v) {
     norm += val * val;
    }
    norm = std::sqrt(norm);

    if (norm < epsilon) continue;

    // 通常の符号を持つ値を選択
    double sign = (v[0] >= 0) ? 1.0 : -1.0;
    double alpha = -sign * norm;

    // w = v + alpha * e_1
    std::vector<double> w(v);
    w[0] += alpha;

    // ノルムを正規化
    double w_norm = 0.0;
    for (double val : w) {
     w_norm += val * val;
    }
    w_norm = std::sqrt(w_norm);

    if (w_norm < epsilon) continue;

    // ハウスホルダー変換の計算
    std::vector<double> u(w);
    for (double &val : u) {
     val /= w_norm;
    }

    // H = I - 2 * u * u^T を計算
    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
     H[i][i] = 1.0;
    }
    for (size_t i = 0; i < n - k; ++i) {
     for (size_t j = 0; j < n - k; ++j) {
      H[k + i][k + j] -= 2.0 * u[i] * u[j];
     }
    }

    // 行列の積を計算
    std::vector<std::vector<double>> temp_R(n, std::vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
     for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k_idx = 0; k_idx < n; ++k_idx) {
       sum += H[i][k_idx] * R[k_idx][j];
      }
      temp_R[i][j] = sum;
     }
    }
    R = temp_R;

    // Q = Q * H^T
    std::vector<std::vector<double>> temp_Q(n, std::vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
     for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k_idx = 0; k_idx < n; ++k_idx) {
       sum += Q[i][k_idx] * H[j][k_idx];
      }
      temp_Q[i][j] = sum;
     }
    }
    Q = temp_Q;
   }

   // 収束判定
   for (size_t i = 0; i < n - 1; ++i) {
    if (std::abs(R[i][i + 1]) > epsilon) {
     converged = false;
     break;
    }
   }

   if (converged) {
    // 固有値を対角成分から取得
    for (size_t i = 0; i < n; ++i) {
     eigenvals[i] = R[i][i];
    }
    break;
   }
  }

  return eigenvals;
 };

 // 固有値計算関数
 std::vector<double> meigenvals(const std::vector<std::vector<double>> &matrix) {
  size_t n = matrix.size();

  // 行列の次元確認
  if (n == 0) { return {}; }

  // 正方行列でない場合はエラー
  for (size_t i = 0; i < n; ++i) {
   if (matrix[i].size() != n) { throw std::runtime_error("Matrix must be square"); }
  }

  // 行列のコピーを作成（元の行列を変更しない）
  std::vector<std::vector<double>> A = matrix;

  // 固有値を求める（QR反復法）
  std::vector<double> eigenvals(n);

  // QR反復法の実装
  const int max_iterations = 1000;
  const double epsilon = 1e-12;

  for (int iter = 0; iter < max_iterations; ++iter) {
   bool converged = true;

   // QR分解を計算
   std::vector<std::vector<double>> Q(n, std::vector<double>(n, 0.0));
   std::vector<std::vector<double>> R = A;

   // 単位行列で初期化
   for (size_t i = 0; i < n; ++i) {
    Q[i][i] = 1.0;
   }

   // ガウス・ホルツ変換でQR分解
   for (size_t k = 0; k < n - 1; ++k) {
    // ベクトルv = R[k:n][k] を取得
    std::vector<double> v(n - k);
    for (size_t i = 0; i < n - k; ++i) {
     v[i] = R[k + i][k];
    }

    // ノルムを計算
    double norm = 0.0;
    for (double val : v) {
     norm += val * val;
    }
    norm = std::sqrt(norm);

    if (norm < epsilon) continue;

    // 通常の符号を持つ値を選択
    double sign = (v[0] >= 0) ? 1.0 : -1.0;
    double alpha = -sign * norm;

    // w = v + alpha * e_1
    std::vector<double> w(v);
    w[0] += alpha;

    // ノルムを正規化
    double w_norm = 0.0;
    for (double val : w) {
     w_norm += val * val;
    }
    w_norm = std::sqrt(w_norm);

    if (w_norm < epsilon) continue;

    // ハウスホルダー変換の計算
    std::vector<double> u(w);
    for (double &val : u) {
     val /= w_norm;
    }

    // H = I - 2 * u * u^T を計算
    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
     H[i][i] = 1.0;
    }
    for (size_t i = 0; i < n - k; ++i) {
     for (size_t j = 0; j < n - k; ++j) {
      H[k + i][k + j] -= 2.0 * u[i] * u[j];
     }
    }

    // 行列の積を計算
    std::vector<std::vector<double>> temp_R(n, std::vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
     for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k_idx = 0; k_idx < n; ++k_idx) {
       sum += H[i][k_idx] * R[k_idx][j];
      }
      temp_R[i][j] = sum;
     }
    }
    R = temp_R;

    // Q = Q * H^T
    std::vector<std::vector<double>> temp_Q(n, std::vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
     for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k_idx = 0; k_idx < n; ++k_idx) {
       sum += Q[i][k_idx] * H[j][k_idx];
      }
      temp_Q[i][j] = sum;
     }
    }
    Q = temp_Q;
   }

   // 収束判定
   for (size_t i = 0; i < n - 1; ++i) {
    if (std::abs(R[i][i + 1]) > epsilon) {
     converged = false;
     break;
    }
   }

   if (converged) {
    // 固有値を対角成分から取得
    for (size_t i = 0; i < n; ++i) {
     eigenvals[i] = R[i][i];
    }
    break;
   }
  }

  return eigenvals;
 }

 size_t nrows(const std::vector<std::vector<double>> &A) { return A.size(); }
 size_t ncols(const std::vector<std::vector<double>> &A) { return A.empty() ? 0 : A[0].size(); }

 std::vector<std::vector<double>> identity(size_t n) {
  std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
  for (size_t i = 0; i < n; ++i)
   I[i][i] = 1.0;
  return I;
 }

 std::vector<double> mat_vec(const std::vector<std::vector<double>> &A, const std::vector<double> &x) {
  size_t n = nrows(A);
  std::vector<double> r(n, 0.0);
  for (size_t i = 0; i < n; ++i)
   for (size_t j = 0; j < x.size(); ++j)
    r[i] += A[i][j] * x[j];
  return r;
 }

 std::vector<std::vector<double>> mat_mul(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B) {
  size_t n = nrows(A), m = ncols(B), k = ncols(A);
  std::vector<std::vector<double>> C(n, std::vector<double>(m, 0.0));
  for (size_t i = 0; i < n; ++i)
   for (size_t j = 0; j < m; ++j)
    for (size_t t = 0; t < k; ++t)
     C[i][j] += A[i][t] * B[t][j];
  return C;
 }

 std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &A) {
  std::vector<std::vector<double>> T(ncols(A), std::vector<double>(nrows(A)));
  for (size_t i = 0; i < nrows(A); ++i)
   for (size_t j = 0; j < ncols(A); ++j)
    T[j][i] = A[i][j];
  return T;
 }

 double norm(const std::vector<double> &x) {
  double s = 0;
  for (double v : x)
   s += v * v;
  return std::sqrt(s);
 }

 double dot(const std::vector<double> &a, const std::vector<double> &b) {
  double s = 0;
  for (size_t i = 0; i < a.size(); ++i)
   s += a[i] * b[i];
  return s;
 }

 double power_method(const std::vector<std::vector<double>> &A, std::vector<double> &eigenvec, int max_iter) {
  size_t n = nrows(A);
  eigenvec.assign(n, 1.0);

  double lambda = 0;

  for (int k = 0; k < max_iter; ++k) {
   std::vector<double> y = mat_vec(A, eigenvec);
   double new_lambda = dot(y, eigenvec) / dot(eigenvec, eigenvec);

   double nrm = norm(y);
   for (size_t i = 0; i < n; ++i)
    eigenvec[i] = y[i] / nrm;

   if (std::fabs(new_lambda - lambda) < cnst_precision_inv) break;

   lambda = new_lambda;
  }

  return lambda;
 }

 void hessenberg(std::vector<std::vector<double>> &A) {
  size_t n = nrows(A);

  for (size_t k = 0; k < n - 2; ++k) {
   std::vector<double> x(n - k - 1);
   for (size_t i = k + 1; i < n; ++i)
    x[i - k - 1] = A[i][k];

   double alpha = norm(x);
   if (alpha < cnst_precision_inv) continue;

   if (x[0] >= 0) alpha = -alpha;
   x[0] -= alpha;

   double beta = dot(x, x);

   for (size_t j = k; j < n; ++j) {
    double s = 0;
    for (size_t i = 0; i < x.size(); ++i)
     s += x[i] * A[k + 1 + i][j];
    s /= beta;

    for (size_t i = 0; i < x.size(); ++i)
     A[k + 1 + i][j] -= 2 * s * x[i];
   }

   for (size_t i = 0; i < n; ++i) {
    double s = 0;
    for (size_t j = 0; j < x.size(); ++j)
     s += x[j] * A[i][k + 1 + j];
    s /= beta;

    for (size_t j = 0; j < x.size(); ++j)
     A[i][k + 1 + j] -= 2 * s * x[j];
   }
  }
 }

 void qr_decompose(const std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &Q, std::vector<std::vector<double>> &R) {
  size_t n = nrows(A);
  Q = identity(n);
  R = A;

  for (size_t j = 0; j < n - 1; ++j) {
   for (size_t i = n - 1; i > j; --i) {
    double a = R[i - 1][j];
    double b = R[i][j];
    if (std::fabs(b) < cnst_precision_inv) continue;

    double r = std::hypot(a, b);
    double c = a / r;
    double s = -b / r;

    for (size_t k = j; k < n; ++k) {
     double t1 = R[i - 1][k];
     double t2 = R[i][k];
     R[i - 1][k] = c * t1 - s * t2;
     R[i][k] = s * t1 + c * t2;
    }

    for (size_t k = 0; k < n; ++k) {
     double t1 = Q[k][i - 1];
     double t2 = Q[k][i];
     Q[k][i - 1] = c * t1 - s * t2;
     Q[k][i] = s * t1 + c * t2;
    }
   }
  }
 }
 std::vector<double> eigenvalues(std::vector<std::vector<double>> A, int max_iter) {
  size_t n = nrows(A);
  hessenberg(A);

  for (int iter = 0; iter < max_iter; ++iter) {
   bool done = true;
   for (size_t i = 1; i < n; ++i)
    if (std::fabs(A[i][i - 1]) > cnst_precision_inv) done = false;
   if (done) break;

   double d = (A[n - 2][n - 2] - A[n - 1][n - 1]) / 2.0;
   double mu = A[n - 1][n - 1] - std::copysign(1.0, d) * A[n - 1][n - 2] * A[n - 1][n - 2] / (std::fabs(d) + std::sqrt(d * d + A[n - 1][n - 2] * A[n - 1][n - 2]));

   for (size_t i = 0; i < n; ++i)
    A[i][i] -= mu;

   std::vector<std::vector<double>> Q, R;
   qr_decompose(A, Q, R);
   A = mat_mul(R, Q);

   for (size_t i = 0; i < n; ++i)
    A[i][i] += mu;
  }

  std::vector<double> eig(n);
  for (size_t i = 0; i < n; ++i)
   eig[i] = A[i][i];

  return eig;
 }
 void lu_decompose(std::vector<std::vector<double>> A, std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U, std::vector<size_t> &P) {
  size_t n = nrows(A);
  L = identity(n);
  U = A;
  P.resize(n);
  for (size_t i = 0; i < n; ++i)
   P[i] = i;

  for (size_t k = 0; k < n; ++k) {
   size_t pivot = k;
   double maxv = std::fabs(U[k][k]);
   for (size_t i = k + 1; i < n; ++i)
    if (std::fabs(U[i][k]) > maxv) {
     maxv = std::fabs(U[i][k]);
     pivot = i;
    }

   if (maxv < cnst_precision_inv) throw std::runtime_error("Singular std::vector<std::vector<double>>");

   if (pivot != k) {
    std::swap(U[k], U[pivot]);
    std::swap(P[k], P[pivot]);
    for (size_t j = 0; j < k; ++j)
     std::swap(L[k][j], L[pivot][j]);
   }

   for (size_t i = k + 1; i < n; ++i) {
    L[i][k] = U[i][k] / U[k][k];
    for (size_t j = k; j < n; ++j)
     U[i][j] -= L[i][k] * U[k][j];
   }
  }
 }
 std::vector<double> inverse_iteration(const std::vector<std::vector<double>> &A, double lambda) {
  size_t n = nrows(A);
  std::vector<std::vector<double>> B = A;
  for (size_t i = 0; i < n; ++i)
   B[i][i] -= lambda;

  std::vector<std::vector<double>> L, U;
  std::vector<size_t> P;
  lu_decompose(B, L, U, P);

  std::vector<double> x(n, 1.0);

  for (int iter = 0; iter < 100; ++iter) {
   std::vector<double> y(n);
   for (size_t i = 0; i < n; ++i)
    y[i] = x[P[i]];

   for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < i; ++j)
     y[i] -= L[i][j] * y[j];

   for (size_t i = n - 1; i >= 0; --i) {
    for (size_t j = i + 1; j < n; ++j)
     y[i] -= U[i][j] * x[j];
    y[i] /= U[i][i];
   }

   double nrm = norm(y);
   for (size_t i = 0; i < n; ++i)
    x[i] = y[i] / nrm;
  }

  return x;
 }

 std::vector<double> singular_values(const std::vector<std::vector<double>> &A) {
  size_t m = A.size();
  if (m == 0) return {};

  size_t n = A[0].size();
  if (n == 0) return {};

  // 🔹 1x1 special case
  if (m == 1 && n == 1) return {std::fabs(A[0][0])};

  // 🔹 一般ケース
  std::vector<std::vector<double>> U, V;
  std::vector<double> S;

  svd_jacobi(A, U, S, V);

  return S;
 }
 double condition_number(const std::vector<std::vector<double>> &A) {
  auto sv = singular_values(A);
  double maxv = *std::max_element(sv.begin(), sv.end());
  double minv = *std::min_element(sv.begin(), sv.end());
  return maxv / minv;
 }

 void bidiagonalize(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &U, std::vector<std::vector<double>> &V) {
  size_t m = nrows(A);
  size_t n = ncols(A);

  U = identity(m);
  V = identity(n);

  for (size_t k = 0; k < std::min(m, n); ++k) {
   // --- 左Householder ---
   std::vector<double> x(m - k);
   for (size_t i = k; i < m; ++i)
    x[i - k] = A[i][k];

   double alpha = norm(x);
   if (alpha > cnst_precision_inv) {
    if (x[0] >= 0) alpha = -alpha;
    x[0] -= alpha;

    double beta = dot(x, x);

    for (size_t j = k; j < n; ++j) {
     double s = 0;
     for (size_t i = 0; i < x.size(); ++i)
      s += x[i] * A[k + i][j];
     s /= beta;
     for (size_t i = 0; i < x.size(); ++i)
      A[k + i][j] -= 2 * s * x[i];
    }

    for (size_t j = 0; j < m; ++j) {
     double s = 0;
     for (size_t i = 0; i < x.size(); ++i)
      s += x[i] * U[j][k + i];
     s /= beta;
     for (size_t i = 0; i < x.size(); ++i)
      U[j][k + i] -= 2 * s * x[i];
    }
   }

   if (k + 1 >= n) continue;

   // --- 右Householder ---
   std::vector<double> y(n - k - 1);
   for (size_t j = k + 1; j < n; ++j)
    y[j - k - 1] = A[k][j];

   alpha = norm(y);
   if (alpha > cnst_precision_inv) {
    if (y[0] >= 0) alpha = -alpha;
    y[0] -= alpha;

    double beta = dot(y, y);

    for (size_t i = k; i < m; ++i) {
     double s = 0;
     for (size_t j = 0; j < y.size(); ++j)
      s += y[j] * A[i][k + 1 + j];
     s /= beta;
     for (size_t j = 0; j < y.size(); ++j)
      A[i][k + 1 + j] -= 2 * s * y[j];
    }

    for (size_t j = 0; j < n; ++j) {
     double s = 0;
     for (size_t i = 0; i < y.size(); ++i)
      s += y[i] * V[j][k + 1 + i];
     s /= beta;
     for (size_t i = 0; i < y.size(); ++i)
      V[j][k + 1 + i] -= 2 * s * y[i];
    }
   }
  }
 }
 void bidiagonal_qr(std::vector<std::vector<double>> &B, std::vector<std::vector<double>> &U, std::vector<std::vector<double>> &V) {
  size_t n = ncols(B);

  for (int iter = 0; iter < 1000; ++iter) {
   bool converged = true;
   for (size_t i = 1; i < n; ++i)
    if (std::fabs(B[i][i - 1]) > cnst_precision_inv) converged = false;
   if (converged) break;

   for (size_t k = 0; k < n - 1; ++k) {
    double a = B[k][k];
    double b = B[k][k + 1];
    double r = std::hypot(a, b);
    if (r < cnst_precision_inv) continue;
    double c = a / r, s = -b / r;

    for (size_t j = k; j < n; ++j) {
     double t1 = B[k][j], t2 = B[k + 1][j];
     B[k][j] = c * t1 - s * t2;
     B[k + 1][j] = s * t1 + c * t2;
    }

    for (size_t j = 0; j < n; ++j) {
     double t1 = V[j][k], t2 = V[j][k + 1];
     V[j][k] = c * t1 - s * t2;
     V[j][k + 1] = s * t1 + c * t2;
    }
   }
  }
 }

 void svd_jacobi(const std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &U, std::vector<double> &S, std::vector<std::vector<double>> &V) {
  size_t m = A.size();
  if (m == 0) throw std::runtime_error("empty matrix");

  size_t n = A[0].size();
  if (n == 0) throw std::runtime_error("empty matrix");

  // ---- 1x1 special case ----
  if (m == 1 && n == 1) {
   U = {{1.0}};
   V = {{1.0}};
   S = {std::fabs(A[0][0])};
   return;
  }

  // ---- m < n の場合は転置 ----
  bool transposed = false;
  std::vector<std::vector<double>> B = A;

  if (m < n) {
   transposed = true;
   B.assign(n, std::vector<double>(m));
   for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
     B[j][i] = A[i][j];
   std::swap(m, n);
  }

  // ---- AtA 計算 ----
  std::vector<std::vector<double>> AtA(n, std::vector<double>(n, 0.0));
  for (size_t i = 0; i < n; ++i)
   for (size_t j = 0; j < n; ++j)
    for (size_t k = 0; k < m; ++k)
     AtA[i][j] += B[k][i] * B[k][j];

  // ---- Jacobi 固有値分解 ----
  V.assign(n, std::vector<double>(n, 0.0));
  for (size_t i = 0; i < n; ++i)
   V[i][i] = 1.0;

  const int max_iter = 100;
  const double eps = 1e-12;

  for (int iter = 0; iter < max_iter; ++iter) {
   bool converged = true;

   for (size_t p = 0; p < n - 1; ++p) {
    for (size_t q = p + 1; q < n; ++q) {

     if (std::abs(AtA[p][q]) < eps) continue;
     converged = false;

     double theta = 0.5 * std::atan2(2.0 * AtA[p][q], AtA[q][q] - AtA[p][p]);

     double c = std::cos(theta);
     double s = std::sin(theta);

     for (size_t i = 0; i < n; ++i) {
      double ip = AtA[i][p];
      double iq = AtA[i][q];
      AtA[i][p] = c * ip - s * iq;
      AtA[i][q] = s * ip + c * iq;
     }

     for (size_t i = 0; i < n; ++i) {
      double pi = AtA[p][i];
      double qi = AtA[q][i];
      AtA[p][i] = c * pi - s * qi;
      AtA[q][i] = s * pi + c * qi;
     }

     for (size_t i = 0; i < n; ++i) {
      double vip = V[i][p];
      double viq = V[i][q];
      V[i][p] = c * vip - s * viq;
      V[i][q] = s * vip + c * viq;
     }
    }
   }

   if (converged) break;
  }

  // ---- 特異値 ----
  S.resize(n);
  for (size_t i = 0; i < n; ++i)
   S[i] = std::sqrt(std::max(0.0, AtA[i][i]));

  // ---- 降順ソート ----
  for (size_t i = 0; i < n - 1; ++i) {
   for (size_t j = i + 1; j < n; ++j) {
    if (S[i] < S[j]) {
     std::swap(S[i], S[j]);
     for (size_t k = 0; k < n; ++k)
      std::swap(V[k][i], V[k][j]);
    }
   }
  }

  // ---- U = B V Σ^-1 ----
  U.assign(m, std::vector<double>(n, 0.0));

  for (size_t i = 0; i < n; ++i) {
   if (S[i] < eps) continue;

   for (size_t r = 0; r < m; ++r)
    for (size_t c = 0; c < n; ++c)
     U[r][i] += B[r][c] * V[c][i];

   for (size_t r = 0; r < m; ++r)
    U[r][i] /= S[i];
  }

  // ---- 転置して戻す ----
  if (transposed) {
   // A = U Σ V^T だったが，元は A^T を分解しているので，U と V を入れ替える
   std::swap(U, V);
  }
 }

 //// svdだけど現状ゴミ
 // void svd(const std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &U, std::vector<double> &S, std::vector<std::vector<double>> &V) {
 //  std::vector<std::vector<double>> B = A;
 //  bidiagonalize(B, U, V);
 //  bidiagonal_qr(B, U, V);

 // size_t n = std::min(nrows(B), ncols(B));
 // S.resize(n);

 // for (size_t i = 0; i < n; ++i)
 //  S[i] = std::fabs(B[i][i]);
 //}

 // 信号系

 bool isPowerOfTwo(size_t n) { return n && !(n & (n - 1)); }

 std::vector<std::complex<double>> toComplexVector1D(const Value &v, FunctionContext &ctx) {
  if (!v.isMulti()) throw CalcError(CalcErrorType::TypeError, "vector required", ctx.pos);

  const auto &mv = v.asMultiRef(ctx.pos);
  std::vector<std::complex<double>> out;
  out.reserve(mv.size());

  for (const auto &e : mv.elems())
   out.push_back(e.toComplex(ctx.pos));

  return out;
 }

 Value fromComplexVector(const std::vector<std::complex<double>> &vec) {
  std::vector<Value> elems;
  elems.reserve(vec.size());

  for (auto &c : vec) {
   if (std::abs(c.imag()) < cnst_precision_inv) elems.emplace_back(c.real());
   else elems.emplace_back(c);
  }
  return Value(std::make_shared<MultiValue>(std::move(elems)));
 }

 std::vector<std::complex<double>> dft_impl(const std::vector<std::complex<double>> &x, bool inverse) {

  size_t N = x.size();
  std::vector<std::complex<double>> X(N);

  double sign = inverse ? 1.0 : -1.0;

  for (size_t k = 0; k < N; ++k) {
   std::complex<double> sum = 0;
   for (size_t n = 0; n < N; ++n) {
    double angle = 2.0 * PI * k * n / N;
    sum += x[n] * std::exp(std::complex<double>(0, sign * angle));
   }
   if (inverse) sum /= static_cast<double>(N);
   X[k] = sum;
  }

  return X;
 }

 void fft_rec(std::vector<std::complex<double>> &a) {

  size_t N = a.size();
  if (N <= 1) return;

  size_t half = N / 2;

  std::vector<std::complex<double>> even(half), odd(half);

  for (size_t i = 0; i < half; ++i) {
   even[i] = a[2 * i];
   odd[i] = a[2 * i + 1];
  }

  fft_rec(even);
  fft_rec(odd);

  std::complex<double> w = 1.0;
  std::complex<double> wn = std::exp(std::complex<double>(0, -2.0 * PI / N));

  for (size_t k = 0; k < half; ++k) {
   std::complex<double> t = w * odd[k];

   a[k] = even[k] + t;
   a[k + half] = even[k] - t;

   w *= wn;
  }
 }

 std::vector<std::complex<double>> fft_dispatch(const std::vector<std::complex<double>> &x, bool inverse, FunctionContext &ctx) {

  if (!isPowerOfTwo(x.size())) {
   calcWarn(ctx, "FFT fallback to O(N^2) DFT (non power-of-two length)");
   return dft_impl(x, inverse);
  }

  std::vector<std::complex<double>> out = x;

  if (inverse) {
   for (auto &c : out)
    c = std::conj(c);
  }

  fft_rec(out);

  if (inverse) {
   for (auto &c : out)
    c = std::conj(c);

   double N = static_cast<double>(out.size());
   for (auto &c : out)
    c /= N;
  }

  return out;
 }

 Value fft2D(const Value &v, bool inverse, FunctionContext &ctx) {

  const auto &mv = v.asMultiRef(ctx.pos);

  size_t rows = mv.size();
  if (rows == 0) return v;

  size_t cols = mv[0].asMultiRef(ctx.pos).size();

  std::vector<std::vector<std::complex<double>>> mat(rows, std::vector<std::complex<double>>(cols));

  for (size_t i = 0; i < rows; i++) {
   const auto &row = mv[i].asMultiRef(ctx.pos);
   for (size_t j = 0; j < cols; j++)
    mat[i][j] = row[j].toComplex(ctx.pos);
  }

  // 行方向
  for (size_t i = 0; i < rows; i++)
   mat[i] = fft_dispatch(mat[i], inverse, ctx);

  // 列方向
  for (size_t j = 0; j < cols; j++) {
   std::vector<std::complex<double>> col(rows);
   for (size_t i = 0; i < rows; i++)
    col[i] = mat[i][j];

   col = fft_dispatch(col, inverse, ctx);

   for (size_t i = 0; i < rows; i++)
    mat[i][j] = col[i];
  }

  // 戻す
  std::vector<Value> outRows;
  outRows.reserve(rows);

  for (size_t i = 0; i < rows; i++) {
   std::vector<Value> rowVals;
   rowVals.reserve(cols);

   for (size_t j = 0; j < cols; j++) {
    auto &c = mat[i][j];
    if (std::abs(c.imag()) < cnst_precision_inv) rowVals.emplace_back(c.real());
    else rowVals.emplace_back(c);
   }

   outRows.emplace_back(std::make_shared<MultiValue>(std::move(rowVals)));
  }

  return Value(std::make_shared<MultiValue>(std::move(outRows)));
 }

} // namespace mm::cal