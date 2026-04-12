#include "math_util.hpp"

#include <intrin.h>

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

  for (size_t i = 0; i < n;) {
   size_t j = i;
   while (j + 1 < n && nearly_equal(x[idx[j + 1]], x[idx[i]]))
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
  const size_t m = nrows(A);
  const size_t n = ncols(A);

  if (m == 0 || n == 0) {
   Q = {};
   R = A;
   return;
  }

  Q = identity(m); // Q は m×m
  R = A;           // R は m×n

  for (size_t j = 0; j < std::min(m, n); ++j) {
   for (size_t i = m - 1; i > j; --i) {
    const double a = R[i - 1][j];
    const double b = R[i][j];

    if (std::fabs(b) < cnst_precision_inv) continue;

    const double r = std::hypot(a, b);
    const double c = a / r;
    const double s = -b / r;

    // R = G * R
    for (size_t k = j; k < n; ++k) {
     const double t1 = R[i - 1][k];
     const double t2 = R[i][k];
     R[i - 1][k] = c * t1 - s * t2;
     R[i][k] = s * t1 + c * t2;
    }

    // Q = Q * G^T
    for (size_t k = 0; k < m; ++k) {
     const double t1 = Q[k][i - 1];
     const double t2 = Q[k][i];
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

   if (maxv < cnst_precision_inv) throw std::runtime_error("Singular std::vector<std::vector<double>]");

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

   for (size_t i = n; i-- > 0;) {
    for (size_t j = i + 1; j < n; ++j)
     y[i] -= U[i][j] * y[j];
    y[i] /= U[i][i];
   }

   double nrm = norm(y);
   for (size_t i = 0; i < n; ++i)
    if (nrm > 1e-14) x[i] = y[i] / nrm;
  }

  return x;
 }

 std::vector<double> singular_values(const std::vector<std::vector<double>> &A) {
  const size_t m = A.size();
  if (m == 0) return {};

  const size_t n = A[0].size();
  if (n == 0) return {};

  for (size_t i = 1; i < m; ++i) {
   if (A[i].size() != n) { throw std::runtime_error("irregular matrix"); }
  }

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

 static std::vector<std::vector<long double>> complete_orthonormal_basis(const std::vector<std::vector<long double>> &Qinit, size_t rows) {
  using ld = long double;

  const ld eps = std::numeric_limits<ld>::epsilon() * 256.0L;

  std::vector<std::vector<ld>> Q(rows, std::vector<ld>(rows, 0.0L));

  auto dot_col = [&](const std::vector<ld> &a, const std::vector<ld> &b) -> ld {
   ld s = 0.0L;
   ld c = 0.0L;
   for (size_t i = 0; i < rows; ++i) {
    const ld x = a[i] * b[i];
    const ld y = x - c;
    const ld t = s + y;
    c = (t - s) - y;
    s = t;
   }
   return s;
  };

  auto norm2 = [&](const std::vector<ld> &v) -> ld { return dot_col(v, v); };

  auto normalize = [&](std::vector<ld> &v) -> bool {
   ld n2 = norm2(v);
   if (n2 <= eps) return false;
   ld inv = 1.0L / std::sqrt(n2);
   for (size_t i = 0; i < rows; ++i)
    v[i] *= inv;
   return true;
  };

  size_t col = 0;

  if (!Qinit.empty()) {
   const size_t init_cols = Qinit[0].size();

   for (size_t j = 0; j < init_cols && col < rows; ++j) {
    std::vector<ld> v(rows, 0.0L);
    for (size_t i = 0; i < rows; ++i)
     v[i] = Qinit[i][j];

    for (size_t k = 0; k < col; ++k) {
     std::vector<ld> qk(rows, 0.0L);
     for (size_t i = 0; i < rows; ++i)
      qk[i] = Q[i][k];
     ld proj = dot_col(v, qk);
     for (size_t i = 0; i < rows; ++i)
      v[i] -= proj * qk[i];
    }

    if (!normalize(v)) continue;

    for (size_t i = 0; i < rows; ++i)
     Q[i][col] = v[i];
    ++col;
   }
  }

  for (size_t base = 0; base < rows && col < rows; ++base) {
   std::vector<ld> v(rows, 0.0L);
   v[base] = 1.0L;

   for (size_t k = 0; k < col; ++k) {
    std::vector<ld> qk(rows, 0.0L);
    for (size_t i = 0; i < rows; ++i)
     qk[i] = Q[i][k];
    ld proj = dot_col(v, qk);
    for (size_t i = 0; i < rows; ++i)
     v[i] -= proj * qk[i];
   }

   if (!normalize(v)) continue;

   for (size_t i = 0; i < rows; ++i)
    Q[i][col] = v[i];
   ++col;
  }

  for (size_t tries = 0; col < rows && tries < rows * 4; ++tries) {
   std::vector<ld> v(rows, 0.0L);
   for (size_t i = 0; i < rows; ++i) {
    v[i] = static_cast<ld>((i + 1) * (tries + 2));
   }

   for (size_t k = 0; k < col; ++k) {
    std::vector<ld> qk(rows, 0.0L);
    for (size_t i = 0; i < rows; ++i)
     qk[i] = Q[i][k];
    ld proj = dot_col(v, qk);
    for (size_t i = 0; i < rows; ++i)
     v[i] -= proj * qk[i];
   }

   if (!normalize(v)) continue;

   for (size_t i = 0; i < rows; ++i)
    Q[i][col] = v[i];
   ++col;
  }

  if (col != rows) { throw std::runtime_error("failed to complete orthonormal basis"); }

  return Q;
 }

 void svd_jacobi(const std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &U, std::vector<double> &S, std::vector<std::vector<double>> &V) {
  using ld = long double;

  const size_t m0 = A.size();
  if (m0 == 0) throw std::runtime_error("empty matrix");

  const size_t n0 = A[0].size();
  if (n0 == 0) throw std::runtime_error("empty matrix");

  for (size_t i = 1; i < m0; ++i) {
   if (A[i].size() != n0) throw std::runtime_error("irregular matrix");
  }

  const size_t r0 = std::min(m0, n0);

  if (m0 == 1 && n0 == 1) {
   U = {{1.0}};
   V = {{1.0}};
   S = {std::fabs(A[0][0])};
   return;
  }

  bool transposed = false;

  size_t m = m0;
  size_t n = n0;
  std::vector<std::vector<ld>> B;

  if (m0 >= n0) {
   B.assign(m0, std::vector<ld>(n0, 0.0L));
   for (size_t i = 0; i < m0; ++i) {
    for (size_t j = 0; j < n0; ++j) {
     B[i][j] = static_cast<ld>(A[i][j]);
    }
   }
  } else {
   transposed = true;
   m = n0;
   n = m0;
   B.assign(m, std::vector<ld>(n, 0.0L));
   for (size_t i = 0; i < m0; ++i) {
    for (size_t j = 0; j < n0; ++j) {
     B[j][i] = static_cast<ld>(A[i][j]);
    }
   }
  }

  const size_t r = n; // 内部では常に m >= n

  std::vector<std::vector<ld>> Vthin(n, std::vector<ld>(n, 0.0L));
  for (size_t i = 0; i < n; ++i) {
   Vthin[i][i] = 1.0L;
  }

  auto col_dot = [&](size_t p, size_t q) -> ld {
   ld s = 0.0L;
   ld c = 0.0L;
   for (size_t i = 0; i < m; ++i) {
    const ld x = B[i][p] * B[i][q];
    const ld y = x - c;
    const ld t = s + y;
    c = (t - s) - y;
    s = t;
   }
   return s;
  };

  auto col_norm2 = [&](size_t p) -> ld {
   ld s = 0.0L;
   ld c = 0.0L;
   for (size_t i = 0; i < m; ++i) {
    const ld x = B[i][p] * B[i][p];
    const ld y = x - c;
    const ld t = s + y;
    c = (t - s) - y;
    s = t;
   }
   return s;
  };

  const int max_sweeps = 100;
  const ld eps = std::numeric_limits<ld>::epsilon();

  for (int sweep = 0; sweep < max_sweeps; ++sweep) {
   bool changed = false;

   for (size_t p = 0; p + 1 < n; ++p) {
    for (size_t q = p + 1; q < n; ++q) {
     const ld app = col_norm2(p);
     const ld aqq = col_norm2(q);
     const ld apq = col_dot(p, q);

     const ld scale = std::sqrt(std::max<ld>(0.0L, app * aqq));
     if (scale == 0.0L) continue;

     if (std::abs(apq) <= eps * scale * 64.0L) continue;

     changed = true;

     const ld tau = (aqq - app) / (2.0L * apq);
     const ld t = (tau >= 0.0L) ? 1.0L / (tau + std::sqrt(1.0L + tau * tau)) : -1.0L / (-tau + std::sqrt(1.0L + tau * tau));
     const ld c = 1.0L / std::sqrt(1.0L + t * t);
     const ld s = c * t;

     for (size_t i = 0; i < m; ++i) {
      const ld bip = B[i][p];
      const ld biq = B[i][q];
      B[i][p] = c * bip - s * biq;
      B[i][q] = s * bip + c * biq;
     }

     for (size_t i = 0; i < n; ++i) {
      const ld vip = Vthin[i][p];
      const ld viq = Vthin[i][q];
      Vthin[i][p] = c * vip - s * viq;
      Vthin[i][q] = s * vip + c * viq;
     }
    }
   }

   if (!changed) break;
  }

  std::vector<ld> Sld(r, 0.0L);
  for (size_t j = 0; j < r; ++j) {
   Sld[j] = std::sqrt(std::max<ld>(0.0L, col_norm2(j)));
  }

  for (size_t i = 0; i + 1 < r; ++i) {
   size_t best = i;
   for (size_t j = i + 1; j < r; ++j) {
    if (Sld[j] > Sld[best]) best = j;
   }
   if (best != i) {
    std::swap(Sld[i], Sld[best]);

    for (size_t row = 0; row < m; ++row) {
     std::swap(B[row][i], B[row][best]);
    }
    for (size_t row = 0; row < n; ++row) {
     std::swap(Vthin[row][i], Vthin[row][best]);
    }
   }
  }

  ld smax = 0.0L;
  for (ld s : Sld)
   smax = std::max(smax, s);

  const ld sval_tol = eps * std::max<ld>(1.0L, smax) * static_cast<ld>(std::max(m, n)) * 128.0L;

  std::vector<std::vector<ld>> Useed(m, std::vector<ld>(r, 0.0L));

  for (size_t j = 0; j < r; ++j) {
   if (Sld[j] <= sval_tol) continue;
   for (size_t i = 0; i < m; ++i) {
    Useed[i][j] = B[i][j] / Sld[j];
   }
  }

  std::vector<std::vector<ld>> Ufull = complete_orthonormal_basis(Useed, m);
  std::vector<std::vector<ld>> Vfull = complete_orthonormal_basis(Vthin, n);

  if (!transposed) {
   U.assign(m0, std::vector<double>(m0, 0.0));
   V.assign(n0, std::vector<double>(n0, 0.0));
   S.assign(r0, 0.0);

   for (size_t i = 0; i < m0; ++i) {
    for (size_t j = 0; j < m0; ++j) {
     U[i][j] = static_cast<double>(Ufull[i][j]);
    }
   }

   for (size_t i = 0; i < n0; ++i) {
    for (size_t j = 0; j < n0; ++j) {
     V[i][j] = static_cast<double>(Vfull[i][j]);
    }
   }

   for (size_t i = 0; i < r0; ++i) {
    S[i] = static_cast<double>(Sld[i]);
   }
  } else {
   // B = A^T = U_b * S * V_b^T
   // A = V_b * S * U_b^T
   // 元のAに対して:
   // U = V_b (m0 x m0)
   // V = U_b (n0 x n0)

   U.assign(m0, std::vector<double>(m0, 0.0));
   V.assign(n0, std::vector<double>(n0, 0.0));
   S.assign(r0, 0.0);

   for (size_t i = 0; i < m0; ++i) {
    for (size_t j = 0; j < m0; ++j) {
     U[i][j] = static_cast<double>(Vfull[i][j]);
    }
   }

   for (size_t i = 0; i < n0; ++i) {
    for (size_t j = 0; j < n0; ++j) {
     V[i][j] = static_cast<double>(Ufull[i][j]);
    }
   }

   for (size_t i = 0; i < r0; ++i) {
    S[i] = static_cast<double>(Sld[i]);
   }
  }
 }

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

  if (!isPowerOfTwo(x.size())) { return dft_impl(x, inverse); }

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
  if (!isPowerOfTwo(rows) || !isPowerOfTwo(cols)) { calcWarn(ctx, "2D FFT fallback to O(N^2) DFT (non power-of-two dimension)"); }

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
 Value fromVector(const std::vector<double> &xs) {
  auto mv = std::make_shared<MultiValue>();
  mv->elems_.reserve(xs.size());

  for (double x : xs) {
   mv->elems_.emplace_back(x);
  }

  return Value(mv);
 }

 Value fromIndexVector(const std::vector<size_t> &xs) {
  auto mv = std::make_shared<MultiValue>();
  mv->elems_.reserve(xs.size());

  for (size_t x : xs) {
   mv->elems_.emplace_back(static_cast<double>(x));
  }

  return Value(mv);
 }

 Value fromMatrix(const std::vector<std::vector<double>> &A) {
  auto outer = std::make_shared<MultiValue>();
  outer->elems_.reserve(A.size());

  for (const auto &rowv : A) {
   auto row = std::make_shared<MultiValue>();
   row->elems_.reserve(rowv.size());

   for (double x : rowv) {
    row->elems_.emplace_back(x);
   }

   outer->elems_.emplace_back(Value(row));
  }

  return Value(outer);
 }
} // namespace mm::cal