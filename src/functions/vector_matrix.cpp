#include "functions/functions.hpp"

namespace mm::cal {

 void registerVector(SystemConfig &cfg) {
  // 任意サイズ行列の作成
  cfg.functions["matrix"] = {0, -1, [](auto &v, auto &ctx) -> Value {
                              if (v.empty()) { throw CalcError(CalcErrorType::TypeError, "TypeError: matrix requires at least one row", ctx.pos); }

                              // 1引数で既に2D matrix なら検証してそのまま返す
                              if (v.size() == 1 && v[0].isMulti()) {
                               const auto &outer = v[0].asMultiRef(ctx.pos);

                               if (outer.empty()) { return v[0]; }

                               bool all_rows_are_multi = true;
                               for (const auto &row : outer.elems()) {
                                if (!row.isMulti()) {
                                 all_rows_are_multi = false;
                                 break;
                                }
                               }

                               if (all_rows_are_multi) {
                                auto A = toMatrixChecked(v[0], ctx, "matrix");
                                (void)A;
                                return v[0];
                               }
                              }

                              // 複数引数: 各引数が row vector であることを要求
                              size_t cols = 0;
                              bool first = true;

                              auto result = std::make_shared<MultiValue>();
                              result->elems_.reserve(v.size());

                              for (const auto &arg : v) {
                               ensureFlatVectorValue(arg, ctx, "matrix");
                               const auto &row = arg.asMultiRef(ctx.pos);

                               if (first) {
                                cols = row.size();
                                first = false;
                               } else if (row.size() != cols) {
                                throw CalcError(CalcErrorType::DomainError, "DomainError: all rows must have the same number of columns", ctx.pos);
                               }

                               result->elems_.push_back(arg);
                              }

                              return Value(result);
                             }};

  // 任意サイズ行列の和
  cfg.functions["madd"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) { throw CalcError(CalcErrorType::TypeError, "madd requires two matrix arguments", ctx.pos); }

                            auto A = toMatrixChecked(v[0], ctx, "madd");
                            auto B = toMatrixChecked(v[1], ctx, "madd");

                            if (A.size() != B.size()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix dimension mismatch", ctx.pos); }

                            if (A.empty()) { return fromMatrix({}); }

                            if (A[0].size() != B[0].size()) { throw CalcError(CalcErrorType::DomainError, "DomainError: row dimension mismatch", ctx.pos); }

                            const size_t rows = A.size();
                            const size_t cols = A[0].size();

                            std::vector<std::vector<double>> C(rows, std::vector<double>(cols, 0.0));
                            for (size_t i = 0; i < rows; ++i) {
                             for (size_t j = 0; j < cols; ++j) {
                              C[i][j] = A[i][j] + B[i][j];
                             }
                            }

                            return fromMatrix(C);
                           }};

  // 任意サイズ行列の積
  cfg.functions["mmul"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti() || !v[1].isMulti()) { throw CalcError(CalcErrorType::TypeError, "mmul requires two matrix arguments", ctx.pos); }

                            auto A = toMatrixChecked(v[0], ctx, "mmul");
                            auto B = toMatrixChecked(v[1], ctx, "mmul");

                            if (A.empty() || B.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos); }

                            if (A[0].size() != B.size()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix dimension mismatch for multiplication", ctx.pos); }

                            auto C = mat_mul(A, B);
                            return fromMatrix(C);
                           }};

  // 行列の転置
  cfg.functions["mtranspose"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  if (!v[0].isMulti()) { throw CalcError(CalcErrorType::TypeError, "mtranspose requires matrix argument", ctx.pos); }

                                  auto A = toMatrixChecked(v[0], ctx, "mtranspose");

                                  if (A.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: cannot transpose empty matrix", ctx.pos); }

                                  const size_t rows = A.size();
                                  const size_t cols = A[0].size();

                                  auto result = std::make_shared<MultiValue>();
                                  result->elems_.reserve(cols);

                                  for (size_t j = 0; j < cols; ++j) {
                                   auto row_result = std::make_shared<MultiValue>();
                                   row_result->elems_.reserve(rows);

                                   for (size_t i = 0; i < rows; ++i) {
                                    row_result->elems_.push_back(Value(A[i][j]));
                                   }

                                   result->elems_.push_back(Value(row_result));
                                  }

                                  return Value(result);
                                 }};

  // 行列のトレース（対角成分の和）
  cfg.functions["mtrace"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "mtrace requires matrix argument", ctx.pos);
                              const auto &matrix = v[0].asMultiRef(ctx.pos);

                              if (matrix.empty() || matrix[0].asMultiRef(ctx.pos).empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: cannot calculate trace of empty matrix", ctx.pos); }

                              size_t min_dim = std::min(matrix.size(), matrix[0].asMultiRef(ctx.pos).size());

                              double trace = 0.0;
                              for (size_t i = 0; i < min_dim; ++i) {
                               trace += matrix[i].asMultiRef(ctx.pos)[i].asScalar(ctx.pos);
                              }
                              return Value(trace);
                             }};

  // 行列のランク
  cfg.functions["mrank"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) { throw CalcError(CalcErrorType::TypeError, "mrank requires matrix argument", ctx.pos); }

                             auto A = toMatrixChecked(v[0], ctx, "mrank");

                             if (A.empty()) { return Value(0.0); }

                             const size_t rows = A.size();
                             const size_t cols = A[0].size();

                             if (cols == 0) { return Value(0.0); }

                             std::vector<std::vector<double>> M = A;
                             const double eps = 1e-10;

                             size_t rank = 0;
                             size_t pivot_col = 0;

                             while (rank < rows && pivot_col < cols) {
                              size_t pivot = rank;
                              double best = std::abs(M[pivot][pivot_col]);

                              for (size_t i = rank + 1; i < rows; ++i) {
                               const double cand = std::abs(M[i][pivot_col]);
                               if (cand > best) {
                                best = cand;
                                pivot = i;
                               }
                              }

                              if (best < eps) {
                               ++pivot_col;
                               continue;
                              }

                              if (pivot != rank) { std::swap(M[rank], M[pivot]); }

                              const double diag = M[rank][pivot_col];

                              for (size_t i = rank + 1; i < rows; ++i) {
                               const double factor = M[i][pivot_col] / diag;
                               if (std::abs(factor) < eps) continue;

                               M[i][pivot_col] = 0.0;
                               for (size_t j = pivot_col + 1; j < cols; ++j) {
                                M[i][j] -= factor * M[rank][j];
                               }
                              }

                              ++rank;
                              ++pivot_col;
                             }

                             return Value(static_cast<double>(rank));
                            }};

  // 行列の決定因子（行列式）
  cfg.functions["mdet"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            if (!v[0].isMulti()) { throw CalcError(CalcErrorType::TypeError, "mdet requires matrix argument", ctx.pos); }

                            auto A = toMatrixChecked(v[0], ctx, "mdet");

                            if (A.empty()) { return Value(0.0); }

                            const size_t rows = A.size();
                            const size_t cols = A[0].size();

                            if (rows != cols) { throw CalcError(CalcErrorType::DomainError, "DomainError: determinant requires square matrix", ctx.pos); }

                            std::vector<std::vector<double>> M = A;
                            const double eps = 1e-10;
                            long double det = 1.0L;
                            int sign = 1;

                            for (size_t k = 0; k < rows; ++k) {
                             size_t pivot = k;
                             double best = std::abs(M[k][k]);

                             for (size_t i = k + 1; i < rows; ++i) {
                              const double cand = std::abs(M[i][k]);
                              if (cand > best) {
                               best = cand;
                               pivot = i;
                              }
                             }

                             if (best < eps) { return Value(0.0); }

                             if (pivot != k) {
                              std::swap(M[k], M[pivot]);
                              sign = -sign;
                             }

                             const double diag = M[k][k];
                             det *= static_cast<long double>(diag);

                             for (size_t i = k + 1; i < rows; ++i) {
                              const double factor = M[i][k] / diag;
                              if (std::abs(factor) < eps) continue;

                              M[i][k] = 0.0;
                              for (size_t j = k + 1; j < cols; ++j) {
                               M[i][j] -= factor * M[k][j];
                              }
                             }
                            }

                            det *= static_cast<long double>(sign);

                            double out = static_cast<double>(det);
                            if (std::abs(out) < eps) out = 0.0;
                            return Value(out);
                           }};

  // 単位行列の作成
  cfg.functions["identity"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                int size = static_cast<int>(v[0].asScalar(ctx.pos));
                                if (size <= 0) throw CalcError(CalcErrorType::DomainError, "DomainError: matrix size must be positive", ctx.pos);

                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(size);

                                for (int i = 0; i < size; ++i) {
                                 auto row = std::make_shared<MultiValue>();
                                 row->elems_.reserve(size);
                                 for (int j = 0; j < size; ++j) {
                                  row->elems_.push_back(Value(i == j ? 1.0 : 0.0));
                                 }
                                 result->elems_.push_back(Value(row));
                                }
                                return Value(result);
                               }};

  // ゼロ行列の作成
  cfg.functions["zeros"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             int rows = static_cast<int>(v[0].asScalar(ctx.pos));
                             int cols = static_cast<int>(v[1].asScalar(ctx.pos));
                             if (rows <= 0 || cols <= 0) throw CalcError(CalcErrorType::DomainError, "DomainError: matrix dimensions must be positive", ctx.pos);

                             auto result = std::make_shared<MultiValue>();
                             result->elems_.reserve(rows);

                             for (int i = 0; i < rows; ++i) {
                              auto row = std::make_shared<MultiValue>();
                              row->elems_.reserve(cols);
                              for (int j = 0; j < cols; ++j) {
                               row->elems_.push_back(Value(0.0));
                              }
                              result->elems_.push_back(Value(row));
                             }
                             return Value(result);
                            }};

  // 行列の要素アクセス
  cfg.functions["mget"] = {3, 3, [](auto &v, auto &ctx) -> Value {
                            auto A = toMatrixChecked(v[0], ctx, "mget");

                            const int row = static_cast<int>(requireInt(v[1], ctx.pos));
                            const int col = static_cast<int>(requireInt(v[2], ctx.pos));

                            if (row < 0 || row >= static_cast<int>(A.size())) { throw CalcError(CalcErrorType::DomainError, "DomainError: row index out of bounds", ctx.pos); }

                            const int cols = A.empty() ? 0 : static_cast<int>(A[0].size());
                            if (col < 0 || col >= cols) { throw CalcError(CalcErrorType::DomainError, "DomainError: column index out of bounds", ctx.pos); }

                            return Value(A[row][col]);
                           }};

  // 行列の行数取得
  cfg.functions["mrows"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             auto A = toMatrixChecked(v[0], ctx, "mrows");
                             return Value(static_cast<double>(A.size()));
                            }};

  // 行列の列数取得
  cfg.functions["mcols"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             auto A = toMatrixChecked(v[0], ctx, "mcols");
                             if (A.empty()) { return Value(0.0); }
                             return Value(static_cast<double>(A[0].size()));
                            }};

  // 行列の対角成分取得
  cfg.functions["mdiag"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             auto A = toMatrixChecked(v[0], ctx, "mdiag");

                             if (A.empty()) { return fromVectorValue({}); }

                             const size_t d = std::min(A.size(), A[0].size());
                             std::vector<double> diag;
                             diag.reserve(d);

                             for (size_t i = 0; i < d; ++i) {
                              diag.push_back(A[i][i]);
                             }

                             return fromVectorValue(diag);
                            }};

  // 行列のLU分解
  cfg.functions["mlu"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           auto A = toMatrix(v[0], ctx);

                           if (A.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: empty matrix", ctx.pos); }
                           if (A.size() != A[0].size()) { throw CalcError(CalcErrorType::DomainError, "DomainError: LU requires square matrix", ctx.pos); }

                           std::vector<std::vector<double>> L, U;
                           std::vector<size_t> P;
                           lu_decompose(A, L, U, P);

                           auto result = std::make_shared<MultiValue>();
                           result->elems_.emplace_back(fromIndexVector(P));
                           result->elems_.emplace_back(fromMatrix(L));
                           result->elems_.emplace_back(fromMatrix(U));
                           return Value(result);
                          }};

  // 行列の固有値
  cfg.functions["meigenvals"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  auto A = toMatrix(v[0], ctx);

                                  if (A.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos); }
                                  if (A.size() != A[0].size()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix must be square", ctx.pos); }
                                  if (!isSymmetricMatrix(A, 1e-10)) { throw CalcError(CalcErrorType::DomainError, "DomainError: meigenvals requires a symmetric real matrix", ctx.pos); }

                                  if (A.size() == 1 && A[0].size() == 1) { return fromVector({A[0][0]}); }

                                  auto eig = eigenvalues(A);
                                  std::sort(eig.begin(), eig.end(), std::greater<double>());
                                  return fromVector(eig);
                                 }};

  // 行列の固有ベクトル
  cfg.functions["meigenvecs"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                  auto A = toMatrix(v[0], ctx);

                                  if (A.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos); }
                                  if (A.size() != A[0].size()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix must be square", ctx.pos); }

                                  const size_t n = A.size();

                                  if (n == 1) { return fromMatrix({{1.0}}); }

                                  if (!isSymmetricMatrix(A, 1e-10)) { throw CalcError(CalcErrorType::DomainError, "DomainError: meigenvecs currently supports only matrices with real eigenvectors reliably; use a symmetric matrix", ctx.pos); }

                                  auto eig = eigenvalues(A);
                                  if (eig.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: failed to compute eigenvector", ctx.pos); }

                                  std::sort(eig.begin(), eig.end(), std::greater<double>());

                                  std::vector<std::vector<double>> vecs;
                                  vecs.reserve(eig.size());

                                  for (double lambda : eig) {
                                   std::vector<std::vector<double>> M = A;
                                   for (size_t i = 0; i < n; ++i) {
                                    M[i][i] -= lambda;
                                   }

                                   std::vector<std::vector<double>> U, V;
                                   std::vector<double> S;
                                   svd_jacobi(M, U, S, V);

                                   if (S.empty() || V.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: failed to compute eigenvector", ctx.pos); }

                                   size_t min_idx = 0;
                                   for (size_t i = 1; i < S.size(); ++i) {
                                    if (S[i] < S[min_idx]) { min_idx = i; }
                                   }

                                   std::vector<double> x(n, 0.0);
                                   for (size_t i = 0; i < n; ++i) {
                                    x[i] = V[i][min_idx];
                                   }

                                   double nrm = 0.0;
                                   for (double xi : x) {
                                    nrm += xi * xi;
                                   }
                                   nrm = std::sqrt(nrm);

                                   if (nearly_zero(nrm)) { throw CalcError(CalcErrorType::DomainError, "DomainError: failed to compute eigenvector", ctx.pos); }

                                   for (double &xi : x) {
                                    xi /= nrm;
                                   }

                                   for (double xi : x) {
                                    if (!nearly_zero(xi)) {
                                     if (xi < 0.0) {
                                      for (double &xj : x) {
                                       xj = -xj;
                                      }
                                     }
                                     break;
                                    }
                                   }

                                   vecs.push_back(std::move(x));
                                  }

                                  return fromMatrix(vecs);
                                 }};

  // 行列の逆行列
  cfg.functions["minverse"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                auto A = toMatrix(v[0], ctx);

                                if (A.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos); }

                                const size_t rows = A.size();
                                const size_t cols = A[0].size();

                                if (rows != cols) { throw CalcError(CalcErrorType::DomainError, "DomainError: inverse requires square matrix", ctx.pos); }

                                for (size_t i = 1; i < rows; ++i) {
                                 if (A[i].size() != cols) { throw CalcError(CalcErrorType::DomainError, "DomainError: irregular matrix", ctx.pos); }
                                }

                                std::vector<std::vector<double>> U, V;
                                std::vector<double> S;
                                svd_jacobi(A, U, S, V);

                                if (S.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos); }

                                const double smax = *std::max_element(S.begin(), S.end());
                                const double tol = std::numeric_limits<double>::epsilon() * std::max(1.0, smax) * static_cast<double>(rows) * 128.0;

                                for (double s : S) {
                                 if (s <= tol) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is singular, cannot invert", ctx.pos); }
                                }

                                // A^{-1} = V * diag(1/S) * U^T
                                std::vector<std::vector<double>> Ainv(rows, std::vector<double>(cols, 0.0));

                                for (size_t i = 0; i < rows; ++i) {
                                 for (size_t j = 0; j < cols; ++j) {
                                  long double sum = 0.0L;
                                  for (size_t k = 0; k < S.size(); ++k) {
                                   sum += static_cast<long double>(V[i][k]) * (1.0L / static_cast<long double>(S[k])) * static_cast<long double>(U[j][k]);
                                  }
                                  Ainv[i][j] = static_cast<double>(sum);
                                 }
                                }

                                return fromMatrix(Ainv);
                               }};

  // 行列のQR分解
  cfg.functions["mqr"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           auto A = toMatrix(v[0], ctx);

                           if (A.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos); }

                           const size_t m = A.size();
                           const size_t n = A[0].size();
                           for (size_t i = 1; i < m; ++i) {
                            if (A[i].size() != n) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix has inconsistent row lengths", ctx.pos); }
                           }

                           std::vector<std::vector<double>> Q, R;
                           qr_decompose(A, Q, R);

                           auto result = std::make_shared<MultiValue>();
                           result->elems_.emplace_back(fromMatrix(Q));
                           result->elems_.emplace_back(fromMatrix(R));
                           return Value(result);
                          }};

  // 行列の特異値分解（SVD）
  cfg.functions["msvd"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            auto A = toMatrix(v[0], ctx);

                            if (A.empty() || A[0].empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos); }

                            std::vector<std::vector<double>> U, V;
                            std::vector<double> S;
                            svd_jacobi(A, U, S, V);

                            // API仕様として msvd は {U, S, Vt} を返す
                            std::vector<std::vector<double>> Vt(V.empty() ? 0 : V[0].size(), std::vector<double>(V.size(), 0.0));

                            for (size_t i = 0; i < V.size(); ++i) {
                             for (size_t j = 0; j < V[i].size(); ++j) {
                              Vt[j][i] = V[i][j];
                             }
                            }

                            auto result = std::make_shared<MultiValue>();
                            result->elems_.emplace_back(fromMatrix(U));
                            result->elems_.emplace_back(fromVector(S));
                            result->elems_.emplace_back(fromMatrix(Vt));
                            return Value(result);
                           }};

  // 行列の軌道（行列式の絶対値）
  cfg.functions["mcond"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) { throw CalcError(CalcErrorType::TypeError, "TypeError: mcond requires matrix argument", ctx.pos); }

                             auto A = toMatrix(v[0], ctx);
                             if (A.empty()) { throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos); }

                             auto S = singular_values(A);
                             if (S.empty()) { return Value(std::numeric_limits<double>::infinity()); }

                             if (A.size() == 1 || A[0].size() == 1) {
                              bool all_zero = true;
                              double nrm2 = 0.0;
                              for (const auto &row : A) {
                               for (double x : row) {
                                nrm2 += x * x;
                                if (!nearly_zero(x)) all_zero = false;
                               }
                              }
                              if (all_zero) { return Value(std::numeric_limits<double>::infinity()); }
                              return Value(1.0);
                             }

                             const double smax = *std::max_element(S.begin(), S.end());
                             const double smin = *std::min_element(S.begin(), S.end());

                             const double tol = std::numeric_limits<double>::epsilon() * std::max(1.0, smax) * static_cast<double>(std::max(A.size(), A[0].size())) * 128.0;

                             if (smin <= tol) { return Value(std::numeric_limits<double>::infinity()); }

                             return Value(smax / smin);
                            }};

  cfg.functions["mcondition"] = cfg.functions["mcond"];

  cfg.functions["mlsqr"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti() || !v[1].isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: mlsqr requires matrix and vector arguments", ctx.pos);

                             const auto &matrix_A = v[0].asMultiRef(ctx.pos);
                             const auto &vector_b = v[1].asMultiRef(ctx.pos);

                             if (matrix_A.empty()) throw CalcError(CalcErrorType::DomainError, "DomainError: mlsqr: matrix is empty", ctx.pos);

                             size_t m = matrix_A.size();
                             size_t n = matrix_A[0].asMultiRef(ctx.pos).size();

                             if (vector_b.size() != m) throw CalcError(CalcErrorType::DomainError, "DomainError: mlsqr: dimension mismatch", ctx.pos);

                             // --- A, b を double に変換 ---
                             std::vector<std::vector<double>> A(m, std::vector<double>(n));
                             for (size_t i = 0; i < m; ++i) {
                              const auto &row = matrix_A[i].asMultiRef(ctx.pos);
                              for (size_t j = 0; j < n; ++j)
                               A[i][j] = row[j].asScalar(ctx.pos);
                             }

                             std::vector<double> b(m);
                             for (size_t i = 0; i < m; ++i)
                              b[i] = vector_b[i].asScalar(ctx.pos);

                             // --- QR分解 ---
                             if (m < n) throw CalcError(CalcErrorType::DomainError, "DomainError: mlsqr requires m >= n (overdetermined system)", ctx.pos);
                             std::vector<std::vector<double>> Q, R;
                             qr_decompose(A, Q, R);

                             // --- y = Q^T b ---
                             std::vector<double> y(n, 0.0);
                             for (size_t i = 0; i < n; ++i)
                              for (size_t k = 0; k < m; ++k)
                               y[i] += Q[k][i] * b[k];

                             // --- 後退代入 R x = y ---
                             std::vector<double> x(n);
                             for (int i = (int)n - 1; i >= 0; --i) {

                              if (nearly_zero(R[i][i])) throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is rank deficient", ctx.pos);

                              x[i] = y[i];
                              for (size_t j = i + 1; j < n; ++j)
                               x[i] -= R[i][j] * x[j];

                              x[i] /= R[i][i];
                             }

                             // --- 結果をMultiValueへ ---
                             auto result = std::make_shared<MultiValue>();
                             result->elems_.reserve(n);

                             for (double val : x)
                              result->elems_.push_back(Value(val));

                             return Value(result);
                            }};

  // 行列のフロベニウスノルム
  cfg.functions["mnormf"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: mnormf requires matrix argument", ctx.pos);
                              const auto &matrix = v[0].asMultiRef(ctx.pos);

                              double sum = 0.0;
                              for (size_t i = 0; i < matrix.size(); ++i) {
                               for (size_t j = 0; j < matrix[i].asMultiRef(ctx.pos).size(); ++j) {
                                double val = matrix[i].asMultiRef(ctx.pos)[j].asScalar(ctx.pos);
                                sum += val * val;
                               }
                              }
                              return Value(std::sqrt(sum));
                             }};

  // 行列のフロベニウスノルム（別名）
  cfg.functions["mnorm"] = cfg.functions["mnormf"];

  // 行列の最大固有値（対称実行列専用, 冪乗法）
  cfg.functions["mmaxeigen_power"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                       std::vector<std::vector<double>> A = toMatrix(v[0], ctx);

                                       size_t m = A.size();
                                       if (m == 0 || m != A[0].size()) throw CalcError(CalcErrorType::DomainError, "DomainError: matrix must be square", ctx.pos);
                                       if (!isSymmetricMatrix(A, 1e-10)) { throw CalcError(CalcErrorType::DomainError, "DomainError: mmaxeigen_power requires a symmetric real matrix", ctx.pos); }

                                       // 1×1 特殊ケース
                                       if (m == 1) {
                                        if (!std::isfinite(A[0][0])) throw CalcError(CalcErrorType::DomainError, "DomainError: invalid number", ctx.pos);

                                        return Value(A[0][0]);
                                       }

                                       std::vector<double> eigenvec;
                                       double lambda = power_method(A, eigenvec);
                                       return Value(lambda);
                                      }};
  cfg.functions["mmaxeigen"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                 std::vector<std::vector<double>> A = toMatrix(v[0], ctx);

                                 size_t m = A.size();
                                 if (m == 0 || m != A[0].size()) throw CalcError(CalcErrorType::DomainError, "DomainError: matrix must be square", ctx.pos);
                                 if (!isSymmetricMatrix(A, 1e-10)) { throw CalcError(CalcErrorType::DomainError, "DomainError: mmaxeigen requires a symmetric real matrix", ctx.pos); }

                                 // 1×1 特殊ケース
                                 if (m == 1) {
                                  if (!std::isfinite(A[0][0])) throw CalcError(CalcErrorType::DomainError, "DomainError: invalid number", ctx.pos);

                                  return Value(A[0][0]);
                                 }

                                 auto eig = eigenvalues(A);
                                 if (eig.empty()) throw CalcError(CalcErrorType::DomainError, "DomainError: eigenvalues empty", ctx.pos);
                                 return *std::max_element(eig.begin(), eig.end());
                                }};

  // 行列の特異値の和
  cfg.functions["msumsv"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              auto A = toMatrix(v[0], ctx);

                              if (A.empty()) throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos);

                              auto S = singular_values(A);

                              double sum = 0.0;
                              for (double s : S)
                               sum += s;

                              return Value(sum);
                             }};

  // 行列の対称性判定
  cfg.functions["misse_symmetric"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                       if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: misse symmetric requires matrix argument", ctx.pos);
                                       const auto &matrix = v[0].asMultiRef(ctx.pos);

                                       size_t rows = matrix.size();
                                       size_t cols = matrix[0].asMultiRef(ctx.pos).size();

                                       if (rows != cols) return Value(false);

                                       for (size_t i = 0; i < rows; ++i) {
                                        for (size_t j = i + 1; j < cols; ++j) {
                                         double a = matrix[i].asMultiRef(ctx.pos)[j].asScalar(ctx.pos);
                                         double b = matrix[j].asMultiRef(ctx.pos)[i].asScalar(ctx.pos);
                                         if (std::abs(a - b) > 1e-10) { return Value(false); }
                                        }
                                       }
                                       return Value(true);
                                      }};

  // 行列の正定性判定
  cfg.functions["mispd"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                             if (!v[0].isMulti()) throw CalcError(CalcErrorType::TypeError, "TypeError: mispd requires matrix argument", ctx.pos);

                             const auto &matrix = v[0].asMultiRef(ctx.pos);

                             // 行列の次元確認
                             if (matrix.empty()) throw CalcError(CalcErrorType::DomainError, "DomainError: matrix is empty", ctx.pos);

                             size_t n = matrix.size(); // 行列の次元

                             // 正方行列でない場合はfalseを返す
                             for (size_t i = 0; i < n; ++i) {
                              if (matrix[i].size() != n) throw CalcError(CalcErrorType::DomainError, "DomainError: matrix must be square", ctx.pos);
                             }

                             // 行列のコピーを作成（元の行列を変更しない）
                             std::vector<std::vector<double>> A(n, std::vector<double>(n));
                             for (size_t i = 0; i < n; ++i) {
                              const auto &row = matrix[i].asMultiRef(ctx.pos);
                              for (size_t j = 0; j < n; ++j) {
                               A[i][j] = row[j].asScalar(ctx.pos);
                              }
                             }

                             // Cholesky分解を試行する
                             try {
                              std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

                              for (size_t i = 0; i < n; ++i) {
                               double sum = 0.0;
                               for (size_t k = 0; k < i; ++k) {
                                sum += L[i][k] * L[i][k];
                               }
                               double diag = A[i][i] - sum;

                               if (diag < -1e-12) { return Value(false); }
                               diag = std::max(0.0, diag);

                               L[i][i] = std::sqrt(diag);

                               if (L[i][i] < 1e-12) { return Value(false); }

                               for (size_t j = i + 1; j < n; ++j) {
                                sum = 0.0;
                                for (size_t k = 0; k < i; ++k) {
                                 sum += L[j][k] * L[i][k];
                                }
                                L[j][i] = (A[j][i] - sum) / L[i][i];
                               }
                              }

                              return Value(true);

                             } catch (const std::exception &) { return Value(false); }
                            }};
 }
} // namespace mm::cal
