# Changelog

## v1.5.2 — in development

- Array / Matrix + Bareiss / `solveLinear` / `nullSpace` / LU / Householder QR / SVD / Eigen regression: internal `2027 / 2027 PASS`, black-box `1504 / 1504 PASS`. Fixed-seed Matrix / FFT checks in `mmCal.Benchmarks --random-only` also pass.

- Reworked `RealInterval::point` so the same `BigFloat` is no longer copied and moved within one initialization expression. This prevents negative point intervals from becoming reversed on MSVC, fixing `RealInterval lower bound exceeds upper bound` failures in SVD/Eigen interval audits, and adds a dedicated negative-point regression test.

### Syntax (breaking change)

- Function-call syntax is now exclusively `name[...]`; parenthesized `name(...)` calls were removed
- Parentheses `()` are grouping-only, while `x(x+1)` for an ordinary identifier is parsed as implicit multiplication
- Legacy `sin(x)` syntax on a known function name raises SyntaxError instead of being silently treated as multiplication
- The Formatter always emits `name[...]` for calls and explicit `x*(...)` for identifier/group multiplication, preserving unambiguous round trips
- Removed function-call delimiter branching from Parser/AST/Lowerer
- Unified formal history access around `In[n]` / `Out[n]` and added negative relative indices; zero is invalid
- Added repeated `@` shorthand for `In[-1]`, `In[-2]`, ... and documented the matching repeated `%` shorthand for `Out[-n]`
- Negative `Out[-n]` counts successful outputs, while negative `In[-n]` counts input slots

### Array / basic linear algebra

- Consolidated Array invariants around shape + row-major flat storage; evaluated nested Arrays are rectangularity-checked and flattened at the `Expr::array` normalization boundary
- Preserved zero-length dimensions explicitly and format otherwise ambiguous empty shapes through `reshape[{}, {...}]` for round-trip safety
- Added `dimensions`, `arrayRank`, `at`, and `reshape`, with shared zero-based Array indexing
- Added `MatrixView` / `MatrixBuffer` so matrix algorithms operate directly on flat Arrays instead of repeatedly materializing `vector<vector<Expr>>`
- Unified canonical names around `dot`, `matrixRank`, `norm`, and `normalize`; legacy matrix/vector names remain compatibility aliases
- Integrated same-shape Array `+` / `-` and scalar×Array into ordinary arithmetic, while Array×Array `*` is rejected in favor of explicit `dot[...]` contraction
- Added a flat exact-Number Gaussian / Gauss-Jordan backend that keeps Expr/Simplifier out of pivot loops for `det`, `inverse`, `rref`, and `matrixRank`
- Added Stage 3 Bareiss fraction-free elimination for integer/Rational matrices. Rational inputs are lifted by per-row denominator LCMs and `det`, `rref`, `matrixRank`, and `inverse` share the `IntegerMatrixBuffer` kernel
- Bareiss divisions verify exactness through `divmod`; pivot selection prefers the smallest nonzero bit length to limit BigInt growth. Exact complex matrices retain the Gaussian fallback
- Added `solveLinear[A,b]`: square systems and consistent overdetermined systems with full column rank return a unique solution, while inconsistent or free-variable systems are Domain errors
- Added `nullSpace[A]`: free columns in ascending order define a canonical RREF basis; full column rank preserves the empty-basis vector dimension as `reshape[{}, {0,n}]`. Integer/Rational inputs share Bareiss forward elimination, exact complex uses Gaussian fallback, and symbolic input proceeds only with provably nonzero pivots
- Exact integer/Rational `solveLinear` reuses per-row denominator clearing and augmented Bareiss `[A|b]`, limiting pivots to coefficient columns for consistency/uniqueness checks; exact complex uses Gaussian fallback and symbolic systems proceed only with provably nonzero pivots
- On the same Release benchmark, Stage 3 improves `det` by about 5.4–11.4x and `rref` by about 7.7–15.5x for orders 8–16
- Added a shared symbolic-expansion budget plus triangular fast path for `det` / `inverse` to prevent factorial expression growth on dense high-order symbolic matrices
- Extracted certified expression→interval conversion, decimalization, and guard-digit refinement into a shared approximation layer used by both FFT and Matrix backends
- Added direct precision-aware certified Matrix dispatch for `N[dot[...],p]`, `N[det[...],p]`, `N[inverse[...],p]`, `N[rref[...],p]`, `N[matrixRank[...],p]`, `N[norm[...],p]`, and related operations
- `N[solveLinear[...],p]` likewise dispatches directly to certified augmented interval elimination without constructing the exact solution first, and never substitutes epsilon guesses when pivots or consistency cannot be certified
- `matrixRank` prefers exact elimination for exact inputs and uses no epsilon threshold in the approximate backend; only certified nonzero pivots advance the rank, and rank deficiency is never guessed from a tolerance
- `nullSpace` is likewise discontinuous at rank deficiency: exact inputs prefer exact pivot structure, while approximate inputs return a basis only when interval elimination certifies the pivot structure
- Extended `at` with zero-based prefix indexing so rank-3 decomposition results can expose factor subarrays through calls such as `at[result,0]`
- Added `luDecomposition[A]` for square matrices with exact-first row-pivoted `P A = L U`; certified approximate LU uses interval-based partial pivoting and `N[...]` dispatches directly to certified interval LU
- Added Householder `qrDecomposition[A]` with `A = Q R`; `N[...]` dispatches directly to certified real/complex interval Householder QR. A measured 4x4 exact expression explosion led to a policy limit of 3x3 for general exact QR, while upper-triangular matrices use an any-size fast path
- Prototyped a column-block Householder-application kernel and benchmarked block sizes 1/8/16/32. Orders 8/16/24 showed no consistent win, so automatic blocking is not enabled; the kernel/benchmark remain for future backend work
- Compacted `N` output without weakening certification: exact terminating decimals remain unpadded, while certified fixed-digit results retain only one zero from a redundant trailing-zero run (`1.000... -> 1.0`, `1.500... -> 1.50`). Requested digits and the certified enclosure remain in metadata
- Added fixed-seed Matrix random checks (exact/Rational `solveLinear` round trip, `nullSpace` basis verification through `A.v==0`, exact inverse round trip, determinant transpose invariance, RREF, certified determinant) and Matrix timings to `mmCal.Benchmarks`
- Generalized `{...}` to a finite brace container: equal-shape children auto-optimize to dense Arrays, while heterogeneous-shape results such as `{Q,R}` / `{U,S,V}` remain general braces and are audited at Matrix boundaries
- Extended `qrDecomposition` to rectangular reduced QR and added real/complex reduced `svd` / `singularValueDecomposition`; the numerical SVD avoids `A^H A` and uses Householder bidiagonalization + one-sided Jacobi
- Added `eigenvalues`, `eigenvectors`, and `eigensystem`. Exact triangular/diagonal/distinct-root Number 2x2 cases remain exact; general `N[...]` dispatches directly to a Complex BigFloat Hessenberg + implicit shifted QR Schur backend with interval audits of Schur/eigenpair relations
- Added an explicit `<stdexcept>` include to `expression_interval.cpp`, where `std::overflow_error` is used directly, removing reliance on a transitive declaration
- Added `mmCal.Benchmarks --matrix-large`. Ten-decimal random matrices model the supplied generator workload; order-1024 parsing alone measured about 10.9 s / 1.99 GB RSS, while direct Expr construction was already about 0.69 GB, exposing a storage bottleneck before dense cubic algorithms


## v1.5.1 — 2026-08-12

v1.5.1 preserves the exact-first CAS foundation of v1.5.0 while concentrating on canonicalization, verification, large-integer arithmetic, high-precision numerical evaluation, and reproducible benchmarking.

### Symbolic / CAS

- Added deterministic strict total ordering over the AST for canonical `Add` ordering
- Added definedness-preserving product/division normal forms
- Strengthened MathKnowledge nonzero facts and reversed-relation inference
- Added generated formatter/parser round-trip tests and fixed radix-prefix, Array-adjacency, and negative-Rational formatting failures
- Added a derivative-back harness across integration rule families
- Added automatic Reference ↔ builtin-registry consistency checks
- Expanded extreme-value approximation tests
- Added FFT plan/twiddle caching across transforms

### Multiprecision / high precision

- Adaptive schoolbook / Karatsuba / Toom-3 BigUInt multiplication
- Dedicated squaring path
- Retained balanced-product-tree factorial with faster leaf construction and one-limb paths
- Added Burnikel–Ziegler division and power-of-two division fast paths
- Added `10^9` chunked and divide-and-conquer decimal conversion
- Added early oversized-value rejection in `tryToUint64`
- Added a directed-rounding-safe fast path for extreme BigFloat exponent gaps
- Switched `Pi` to binary-splitting Chudnovsky
- Switched `exp/E` and `log` to binary-splitting based certified evaluation with range reduction
- Added certified argument reduction for huge-radian `sin/cos/tan`

### Benchmark / test

- Added `mmCal.Benchmarks` as a separate Visual Studio project
- Added fixed-seed randomized invariants, threshold sweeps, factorial/decimal-I/O benchmarks, and high-precision `Pi/exp/log` benchmarks
- v1.5.1 release state: internal 1691 / 1691, black-box 1337 / 1337

### Benchmarked but not selected

- Prime-Swing factorial
- binary GCD
- Karatsuba vector-pool workspace
- Karatsuba recursion-depth scratch workspace
- dedicated Toom-3 squaring
- low Toom-3 thresholds
- machine-`fmod` huge-trig reduction

See `docs/performance_optimization.md` for the measured rationale behind each decision.

---

## v1.5.0

A near-complete reconstruction of the old calculator: numerical model, Lexer/Parser/AST/Evaluator, Simplifier, Solver, CertifiedEvaluator, CLI, formatter, tests, and documentation were reorganized around an exact-first CLI calculator / compact CAS architecture.
