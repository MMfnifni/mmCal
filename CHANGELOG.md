# Changelog

## Unreleased

### Development and validation

- Added a grammar-aware semantic expression fuzzer to `mmCal.Benchmarks --random-expressions`.
- `--loop` runs without a case limit and stops immediately on the first FAIL after shrinking the expression and printing seed/case reproduction information.
- Every case is derived independently from the master seed and a 1-based case index, so `--seed N --case M` reproduces the exact case directly.
- Expression depth uses a weighted distribution: shallow expressions dominate, while deeper cases are sampled occasionally up to the configurable `--max-depth` (16 by default).
- Initial invariants cover exact formatter round-trip stability, value preservation by `fullSimplify`, polynomial value preservation by `expand`/`factor`, double transpose, and `det[A] == det[transpose[A]]`.
- The generator focuses on semantically valid exact arithmetic, polynomials, and small matrices instead of flooding the run with random TypeErrors.
- Added `--threads N` for case-level parallel fuzzing with one independent `KernelSession` per worker; master-seed + case replay remains independent of thread count.
- Added `--nostop-loop`; normal `--loop` stops on the first FAIL, while `--nostop-loop` shrinks/reports failures and continues.
- Split fuzzer case isolation into `KernelSession::resetForIndependentEvaluation()`, which clears definitions, histories, input numbering, and diagnostics while preserving the RNG stream and angle settings; long-running loops no longer accumulate session history and no longer entropy-reseed on every reset.
- Replaced the maximum-payload-sized `Expr::Node` `std::variant` with an `ExprKind` header plus kind-specific typed payload nodes as a standalone representation refactor; the public `Expr` API, node identity, and structural equality are preserved.
- In the same GCC Release/LTO-off `--matrix-large transpose 1024 16` run, maximum RSS fell from `693312 KiB` (~677.1 MiB) with the legacy variant node to `299668 KiB` (~292.6 MiB) with typed nodes, a reduction of about 56.8% / 384.4 MiB.
- Reworked persistent dense `ArrayExpr` storage around immutable shared pages of 1024 elements. Integer, Rational, Number, DecimalApproximation, ComplexDecimalApproximation, and Generic values are packed per page while Arrays remain a single public `ExprKind::Array`.
- `ArrayBuilder` limits promotion to the current page and never rebuilds completed pages. If the last element of a million-element numeric Array becomes symbolic, only the final page of at most 1024 elements becomes Generic; preceding packed pages remain shared.
- Added immutable shape/offset/stride views to `ArrayExpr`. `transpose` now swaps shape/strides without copying the backing storage; contiguous reshape and eligible prefix slices also preserve the backing.
- Rectangular brace literals now stream leaves through one `ArrayBuilder`, and numeric literals enter packed pages without first allocating one `Expr` node per element. Ragged braces still lower to general `ListExpr` values.
- Rejected the naive single-`vector<Rational>` packed design because transpose deep-copied one million Rational/BigInt values and took roughly 650–675 ms. The adopted shared-page/view design measured about 0.059 ms for `--matrix-large transpose 1024 16` with about 133.4 MiB maximum RSS. A 13.63 MB 1024x1024 decimal-literal CLI `dimensions[...]` run measured about 4.55 s / 431 MiB.
- BigUInt / BigInt SBO remains deferred so its maintainability and benefit can be evaluated independently.
- Canonicalize nested positive exact integer powers with the branch-safe rule `(a^m)^n -> a^(mn)` in the normal Simplifier, fixing random-fuzzer `expand`/`factor` invariant failures.
- Formatter now preserves left-nested Power associativity as `(a^b)^c` while retaining right-associative `a^b^c` notation.
- Fold consecutive exact-rational subtraction constants so nested powers of `x-1-3` normalize to `(x-4)^144`.
- For positive exact-integer powers `(c*a)^n`, safely raise only the exact numeric coefficient `c` and extract it, allowing `(861(1-x)^64)^3 -> 638277381(1-x)^192` without applying product-power distribution to general complex exponents.

## v1.5.2 — 2026-08-13

v1.5.2 preserves the exact-first numerical foundation of v1.5.1 while substantially expanding symbolic calculus, special functions, precision-aware evaluation, Arrays, and linear algebra. Release state: internal `2027 / 2027 PASS`, black-box `1504 / 1504 PASS`; fixed-seed BigInt, special-function, Matrix, and FFT invariants in `mmCal.Benchmarks --random-only` also pass.

### Syntax and REPL (including breaking changes)

- Function calls are now exclusively `name[...]`; `name(...)` is removed and `()` is grouping only.
- Ordinary identifier adjacency such as `x(x+1)` remains implicit multiplication and formats canonically as `x*(x+1)`.
- Parser/AST/Lowerer call-delimiter branching was removed; legacy `sin(x)` on a known function is a SyntaxError.
- History access is standardized on `In[n]` / `Out[n]`: positive indices are absolute, negative indices are relative, and zero is invalid.
- `@` / `@@` / ... map to `In[-1]` / `In[-2]` / ...; `%` / `%%` / ... map to `Out[-1]` / `Out[-2]` / .... Relative `In` counts input slots, while relative `Out` counts successful outputs.

### Symbolic calculus and integration Knowledge

- Added shared finite-Fourier reduction for nonnegative integer powers of `sin[u]^m cos[u]^n`, avoiding per-power rules.
- Added reciprocal-trigonometric recurrence handling for negative integer sine/cosine powers, including results such as `integrate[sin[2x]^(-2),x] -> -cot[2x]/2`.
- Added integer-power reduction for `tan/cot/sec/csc` and cross-frequency trigonometric product-to-sum rules.
- Added bounded Weierstrass substitution `t=tan[x/2]` for rational expressions in a common `sin/cos` argument, feeding the transformed expression into the exact rational integrator.
- Strengthened structural inverse-chain/substitution matching and concrete quadratic-radical families.
- Retained derivative-back auditing while explicitly allowing `ResolutionOnly` cases so proof-engine limitations do not remove useful primitives.
- Integration diagnostics now distinguish `unsupported`, `partial`, `conditionsRequired`, and `noKnownClosedForm`, avoiding the false implication that an unimplemented method proves mathematical impossibility.
- The broad `docs/memorandum/integralCatalog.md` catalog was used to audit trigonometric, rational, special-function, and branch-sensitive families.

### Special functions

- Added `fresnelc` / `fresnels` with exact special values, odd symmetry, differentiation, certified real `N`, and quadratic-phase integration support.
- Added `hypergeometric1F1[a,b,z]` with terminating/safe exact reductions, differentiation, certified real `N`, and branch-safer `integrate[exp[x^n],x]` representations.
- Added `hypergeometric2F1[a,b,c,z]` with terminating exact series, differentiation, certified real `N`, and binomial-power integration families.
- Added incomplete `ellipticF` / `ellipticE` / `ellipticPi` with principal-branch semantics, amplitude derivatives, certified real `N`, and standard-kernel integration.
- Added `Ei` / `Si` / `Ci` / `li` / `polylog`, including representative exact reductions, derivatives, certified real `N`, and integration Knowledge.
- General inverse special functions are not invented for `solve`; only safe exact degeneracies fall through to existing solvers.

### Precision-aware `N` and FFT

- Extended `N[expr,p]` into a precision-aware evaluation entry point for opted-in builtins rather than always materializing a complete exact result first.
- `N[fft[data],p]` dispatches directly to certified BigFloat/`ComplexInterval` radix-2 transforms instead of constructing huge exact Fourier expressions.
- Certified non-power-of-two FFT uses direct DFT for smaller cases and Bluestein for larger ones; the current measured policy crossover is around 96 points and remains environment-dependent.
- FFT and Matrix paths share expression-to-interval conversion, decimalization, and guard-digit refinement helpers.
- Certified fixed-digit display compresses redundant trailing-zero runs while preserving precision/enclosure metadata (`1.000... -> 1.0`, `1.500... -> 1.50`); exact finite decimals remain compact (`N[1/2,10] -> 0.5`).

### Arrays and linear algebra

- `{...}` is now a general finite brace container. Rectangular children are automatically optimized into dense row-major `ArrayExpr`; differently shaped factors such as `{Q,R}` / `{U,S,V}` remain general brace values.
- Zero-length dimensions are preserved; empty shapes that cannot round-trip through braces alone format via `reshape[{}, {...}]`.
- Added/regularized `dimensions`, `arrayRank`, `length`, prefix-aware `at`, and `reshape`, using zero-based indexing.
- Added `MatrixView` / `MatrixBuffer` so row-major Arrays are not unnecessarily copied into nested vectors.
- Canonicalized core APIs around `dot`, `matrixRank`, `norm`, and `normalize`, retaining legacy aliases.
- Added Bareiss fraction-free elimination for Integer/Rational matrices. Per-row denominator clearing lifts to integer buffers shared by `det`, `rref`, `matrixRank`, `inverse`, `solveLinear`, and `nullSpace`.
- Added `solveLinear[A,b]` for unique solutions, including consistent overdetermined full-column-rank systems; inconsistent or underdetermined systems produce Domain errors rather than invented parameterizations.
- Added `nullSpace[A]` with a deterministic free-column basis while preserving `{0,n}` shape for full-column-rank empty bases.
- Added `luDecomposition[A]` with exact row-pivoted `P A = L U` and direct certified interval partial pivoting under `N[...]`.
- Added rectangular reduced Householder `qrDecomposition[A]`. General exact QR is policy-limited to 3x3 after measured radical-expression explosion, with triangular/trapezoidal fast paths retained.
- A column-block Householder kernel was benchmarked at block sizes 1/8/16/32; no stable winner emerged, so automatic blocking is not enabled.
- Added real/complex reduced `svd` / `singularValueDecomposition`, deliberately avoiding `A^H A`; the numerical backend uses Householder bidiagonalization plus one-sided Jacobi and interval-audits reconstruction/orthogonality.
- Added `conjugateTranspose` for shared Hermitian relations.
- Added `eigenvalues` / `eigenvectors` / `eigensystem`: exact triangular/diagonal and distinct-root exact Number 2x2 cases remain exact, while general `N[...]` uses Complex BigFloat Hessenberg reduction, implicit shifted QR, and Schur relations.
- Symbolic `det` / `inverse` use triangular fast paths and a shared expansion budget to avoid factorial expression growth.
- Discontinuous rank/null-space decisions use exact pivots or interval-certified structure; no arbitrary epsilon threshold is introduced.

### Performance, stability, and development infrastructure

- Bareiss measured about 5.4–11.4x faster for order-8–16 determinants and 7.7–15.5x faster for RREF than the Stage-2 Gaussian/Gauss-Jordan path.
- Added `mmCal.Benchmarks --matrix-large` for order-32/64 Matrix timings and order-1024 storage/parse stress measurements.
- Directly constructing a 1024x1024 ten-decimal Rational Expr matrix reached about 0.69 GB maximum RSS; parsing generator-style text and evaluating only `dimensions` took about 10.9 s and 1.99 GB, showing that representation cost precedes algorithmic cost at this scale.
- Recorded next-version performance ToDos: typed `Expr::Node` storage instead of the large variant fixed cost, BigUInt/BigInt SBO, packed numeric Array/approximate Matrix storage, and fewer temporary allocations while parsing huge braces.
- Fixed evaluation-order-sensitive `RealInterval::point` construction that could reverse negative point intervals under MSVC and fail SVD/Eigen audits; added a direct negative-point regression.
- Added an explicit `<stdexcept>` include where `std::overflow_error` is used instead of relying on transitive declarations.
- Expanded fixed-seed Matrix invariants through Bareiss, inverse, solve, nullSpace, LU, QR, real/complex SVD, and Eigen, alongside FFT round-trip checks.

### Documentation and licensing

- Synchronized README, Reference, Architecture, Roadmap, and performance documentation with the v1.5.2 release state.
- Aligned the BSD 3-Clause copyright notice with project metadata and explicitly separated trademark/brand policy into `TRADEMARKS.md` / `TRADEMARKS.ja.md` without changing the copyright permissions granted by BSD-3-Clause.

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
