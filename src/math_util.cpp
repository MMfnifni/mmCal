#include "math_util.hpp"

namespace mm::cal {

 /* ============================
   数学補助系
   ============================ */

 bool eq(double a, double b) { return std::abs(a - b) < cnst_precision_inv; }

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
  if (n < 0) throw CalcError(CalcErrorType::DomainError, "fib: n < 0", pos);

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

 bool isPrimeLL(long long n) {
  if (n <= 1) return false;                   // 1以下は素数ではない
  if (n <= 3) return true;                    // 2と3は素数
  if (n % 2 == 0 || n % 3 == 0) return false; // 2または3で割り切れる数は素数ではない

  // 5から始めて、iを6ずつ増やしていく
  for (long long i = 5; i * i <= n; i += 6) {
   if (n % i == 0 || n % (i + 2) == 0) return false; // iかi+2で割り切れる場合は素数ではない
  }

  return true; // それ以外は素数
 }

 /* ============================
   統計系
   ============================ */
 std::vector<double> gatherReals(const std::vector<Value> &v, size_t pos) {
  std::vector<double> a;
  a.reserve(v.size());
  for (const auto &x : v) {
   if (isComplex(x)) {
    const auto z = asComplex(x);
    if (std::abs(std::imag(z)) > cnst_precision_inv) throwDomain(pos);
    a.push_back(std::real(z));
   } else {
    a.push_back(asReal(x, pos));
   }
  }
  return a;
 }

 std::vector<double> collectReals(const std::vector<Value> &v, FunctionContext &ctx) {
  std::vector<double> a;
  a.reserve(v.size());
  for (const auto &x : v)
   a.push_back(asReal(x, ctx.pos));
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

} // namespace mm::cal