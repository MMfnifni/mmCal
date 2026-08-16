# Changelog

## Unreleased

- Add post-v1.5.3 changes here.

## v1.5.3 — 2026-08-16

v1.5.3 builds on the symbolic-calculus, certified-`N`, Array, and linear-algebra foundation of v1.5.2 and connects **exact algebraic number fields, Solver semantic normalization, high-precision special functions, exact Cyclotomic FFT, and representation-level optimization** into one release. The exact-first, proof-only policy remains unchanged: unsupported or unproven cases are not promoted to guessed results. Release validation: internal `2335 / 2335 PASS`, black-box `1715 / 1715 PASS`, Random Expression Fuzzer `100000 / 100000 PASS` with 8 threads.

### Algebraic numbers and number fields

- Unified Real/Complex `root[...]` under `AlgebraicNumber`, with bounded minimal-polynomial reduction, primitive-element reduction, resultant fallback, and certified root re-identification.
- Added immutable `NumberFieldContext` / `AlgebraicElement` storage that keeps the chosen embedding and exact Rational power-basis coordinates across arithmetic. Same-field `+ - * /`, small integer powers, and exact inversion stay inside the field.
- Added bounded weak interning for identical embedded generators, compositum/embedding reuse, reciprocal LRU reuse, minimal-polynomial LRU reuse, and incremental exact Krylov elimination for minimal-polynomial derivation.
- Added exact algebraic equality and Real ordering. Deterministic Complex Root enumeration remains separate from mathematical ordering.
- Added `symbolic::exactAlgebraicValue` to bridge canonical `root[...]`, exact Rational/complex Rational, `sqrt` / `cbrt`, `Phi`, and bounded exact arithmetic without forcing a new display form. Comparison, domain reasoning, and direct Solve bindings share this view.
- Consolidated Root coefficient/index parsing and canonical Root Expr construction in the `algebraic_expression` boundary, removing duplicate Evaluator/Solver helpers.

### Solver and semantic coherence

- Added `solve[equation,Integer|Rational|Real|Complex]`, inferring the unknown only when exactly one user symbol is eligible and rejecting protected symbols as solve variables.
- Added solve-safe normalization so representation differences such as `E^x` / `exp[x]`, `ln` / `log`, and `log2` / `log10` do not create Solver capability differences. Proven positive constant-base exponentials are inverted through logarithms.
- Added exact symbolic `lambertw[z]` / `lambertw[k,z]` and certified Real `N` for branches `k=0/-1`. The supported `a^x==x^2` Real family is closed through Lambert W only when branch conditions are proven.
- Added negative domain facts `provablyNonInteger` / `provablyNonRational` and propagated them into membership and Solve constraints.
- `N[SolutionSet]` preserves solution structure while certifiedly approximating only numeric-closed binding right-hand sides.
- Formatter output now places one space around comparison operators `== != < <= > >=` while retaining compact arithmetic formatting.

### Certified numerical evaluation and special functions

- Promoted `DecimalApproximation` / `ComplexDecimalApproximation` to first-class certified numeric leaves with separate `CertifiedEnclosure` and `InformationEnclosure` propagation.
- Standardized `N[expr,p]` as significant decimal digits. If whole-expression certification does not close, numeric subtrees may be approximated structurally without violating Hold attributes; duplicate outer `N` warnings are suppressed when an inner diagnostic already explains the failure.
- Added `zeta`, `digamma`, `trigamma`, and regularized `ibeta`, connected to exact reductions, derivatives, and certified `N` on supported real domains.
- Optimized Gamma/Beta evaluation through exact-Rational argument retention, balanced rising products, a static exact Bernoulli table, high-precision BigInt remainder planning, and point/shared-normalization `ibeta` paths.
- Added `CertifiedBackendUnsupported` to distinguish missing certified backends from true mathematical DomainErrors.

### Exact FFT and Array representation

- Stabilized the exact Cyclotomic FFT backend. Supported non-power-of-two exact inputs are transformed in Rational power-basis coordinates over `Q[t]/Phi_n(t)`, allowing exact 5/7/10/12-point round trips without large `cis[...]` expression growth. Budget or membership failures fall back to the previous generic exact DFT.
- Replaced the large inline `Expr::Node` variant with kind-specific typed nodes while preserving the public `Expr` API and structural semantics.
- Reworked dense `ArrayExpr` storage to immutable paged packed backing plus shape/offset/strides views. Transpose is a zero-copy stride view, and rectangular numeric braces lower directly through `ArrayBuilder`.

### Functions, utilities, and diagnostics

- Added deterministic exact `isprime`, `nextprime`, `prevprime`, `factorint`, and `totient` over the `uint64` range. Larger BigInts are not promoted from probable-prime evidence to exact truth.
- Added BigInt bit utilities, `round[x,n]`, `fma`, `clamp`, and `proj`.
- Added `range`, `table`, `map`, `explain`, and parameterized Real periodic solution families.
- Syntax failures before evaluation no longer consume `In[n]`. Indeterminate Infinity forms such as `Infinity-Infinity` and `0*Infinity` remain conservatively unevaluated instead of simplifying incorrectly.
- Removed the unused Machine/double shim from Core, keeping Exact, BigFloat, and certified backends explicitly separated.

### Performance and validation

- Added the persistent `mmCal.Benchmarks --random-expressions` semantic fuzzer with seed+case reproducibility independent of thread count and session reset support for long-running burn-in.
- Added `test_set/tester.py --timings` and dedicated algebraic-field, special-function, and exact-cyclotomic FFT benchmarks. Adoption and rejection rationale is recorded in `docs/performance_optimization.*`.
- Persistent multiplication-matrix caching, eager generation of larger Bernoulli tables, and unconditional high-precision Stirling-K expansion were rejected after measurement showed insufficient benefit or regressions.
- Release preparation consolidated Root-expression helpers and removed constructor shadowing; Visual Studio project XML validation and GCC shadow-warning checks are part of the release checklist.

### Documentation and compatibility

- Synchronized README, Reference, Architecture, Roadmap, and performance documentation with the v1.5.3 implementation state. Intentional omissions and search budgets are tracked separately in `docs/roadmap.*` and the v1.5.3 intentional-limits memorandum.
- User-visible algebraic canonical output remains `root[minpoly,k]`; internal primitive elements and cache representations are not exposed by the Formatter.

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

- Added the Stage 7-7 exact Cyclotomic FFT backend. Exact Rational inputs of non-power-of-two length 5 or greater are transformed in Rational power-basis coordinates of the quotient `Q[t]/Phi_n(t)` instead of repeatedly simplifying generic `cis[...]` expressions; Gaussian Rational inputs are embedded exactly in `Q(zeta_lcm(n,4))` when needed. Exact 5/7/10/12-point `ifft[fft[v]]` now closes back to `v` without leaving large root-of-unity expressions. The internal generator does not require root isolation and is materialized using the existing single `cis[-2 Pi/n Rad]` vocabulary. Inputs exceeding the degree-64 budget or symbolic inputs that cannot be certified into the quotient field retain the pre-Stage-7-7 generic exact-DFT fallback. The previous dispatch is preserved next to the replacement as a commented reference with the reason for the change. Added `--exact-cyclotomic-fft`; GCC Release/LTO-off warm round trips measured about 0.49/1.74/1.33/1.27 ms for lengths 5/7/10/12.
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
