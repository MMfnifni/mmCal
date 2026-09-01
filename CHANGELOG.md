# Changelog

## v1.5.5 — Unreleased

### Series and asymptotics

- Added `series[...]` on a `SeriesData` / TPSA core. Exact Taylor, Laurent, Puiseux, and logarithmic coefficient layers share one representation and connect directly to `D`, `integrate`, and `normal`. In addition to elementary functions, supported local providers cover `erf` / `erfc` / `Si` / `Ei` / `Ci` / `li`, Fresnel functions, principal inverse trigonometric functions, Lambert W, `gamma` / `lgamma` / `digamma` / `trigamma`, and supported `polylog` neighborhoods using exact coefficient recurrences and DLMF formulas. `tan/cot/sec/csc`, hyperbolic variants, `expm1/log1p`, cardinal functions, and `log2/log10` lower into the same TPSA algebra.
- Added `series[expr,{x,Infinity,n}]`, where `Infinity` denotes real `+Infinity` and is mapped to the local variable `t=1/x`. Rational functions, polynomial growth, reciprocal composition, Puiseux terms, and logarithmic asymptotics reuse the existing Series algebra. When a positive leading direction is proved, `log[A(x)]` is factored as `log[c]+r log[1/x]+log[1+h]`. Oscillatory forms, essential growth, negative logarithmic powers, and other transseries that need a richer representation remain unevaluated rather than guessed.
- Added a limited finite-point `limit` fallback that uses local Series only when the existing limit kernel is unresolved, including exact cancellation of singular terms. `toNormal[expr]` recursively normalizes `SeriesData` inside lists, arrays, and calls; for finite and conditional `SolutionSet` values it transforms only binding right-hand sides while preserving conditions, free variables, multiplicity, and domains.
- Fixed `N[SeriesData]` incorrectly approximating structural exponent-grid metadata and thereby breaking later `normal` / `D`. Only coefficients and the center are approximated; `minimumExponent`, `orderNumerator`, and `exponentDenominator` remain exact integers.
- Fixed false `Division by zero` failures in positive-infinity reciprocal composition such as `tan[1/x]`, `sec[1/x]`, `cot[1/x]`, and `csc[1/x]`. Reciprocal normalization is restricted to the local-variable substitution path and does not change global Simplifier semantics.

### Arrays, vectors, and vector calculus

- Added `Array / scalar` as elementwise scalar scaling. `scalar / Array` and `Array / Array` remain rejected rather than silently acquiring elementwise division semantics.
- Added `inner[a,b]`, `outer[a,b]`, `distance[a,b]`, and `projection[a,b]`. `dot` remains a bilinear contraction, while `inner` is Hermitian and conjugates its first argument. `norm` remains consistent with the Hermitian convention; `vdistance` / `veuclidean` / `vproject` route to the canonical kernels. `vmanhattan` now naturally supports complex components as a sum of component magnitudes.
- Added Cartesian `grad`, `divergence`, `curl`, `laplacian`, `jacobian`, and `hessian`. Coordinates are explicit and the operators reuse the existing exact symbolic `D` kernel. `curl` is restricted to three dimensions; cylindrical or spherical scale factors are not inferred implicitly. Nested vector-calculus operators can be evaluated compositionally.
- Extended `at` to finite `SolutionSet` values. `at[solutions,i]` uses zero-based indexing and returns a one-branch `SolutionSet` while preserving conditions, free variables, multiplicity, and solver-variable domains. `at[solutions,i,x]` extracts the requested binding right-hand side. Non-finite sets, out-of-range indices, and unbound variables are rejected explicitly.

### Symbolic evaluation, CLI, and FFT

- Fixed `D[cases[...]]` so explicit default branches do not invent derivative values at boundaries, compacted additive output from higher derivatives, and extended exact-rational linear-combination and rational-affine assumption handling in `simplify` / `fullSimplify`.
- Added REPL `:quit` / `:exit` commands and interactive `:layout` composition.
- For power-of-two exact `ifft[fft[v]]` within the current degree budget, FFT-generated root-of-unity expressions are re-embedded into cyclotomic coordinates to reduce inverse-expression growth without changing the forward representation.

## v1.5.4 - 2026-08-31

### Solver and symbolic computation

- Expanded multivariate polynomial `solve`. Gröbner elimination and recursive specialization now completely enumerate more supported zero-dimensional systems, while positive-dimensional systems such as `x y==0` or `x^2+y^2==1` can return free parameters with exact conditions when this can be proved. mmCal still does not invent incomplete parametrizations for general varieties.
- Extended Real `solve` for `u exp[u]=a`, affine exponential equations, positive constant-base exponentials, selected `x^x==r` cases, `lambertw[u]==r`, `cosh[u]==r`, and related forms. New range, monotonicity, and convexity proofs return an empty set or a unique exact solution only when the result is proved.
- Added Real solving for `abs` and principal radical equations. Extraneous roots introduced by squaring are removed exactly, while Complex loci such as circles are not collapsed into finite point sets.
- Added first-class mathematical case expressions, `cases[value if condition; ...]`, integrated with `simplify`, `N`, `D`, `integrate`, and `limit`. Unknown conditions are retained and false branches are not evaluated.
- Added `groebnerBasis[...]` / `polynomialReduce[...]` with exact Rational Gröbner bases in Lex / GrLex / GrevLex orderings, also used by nonlinear polynomial `solve`.
- Strengthened `simplify` / `fullSimplify` to preserve definedness. Rewrites such as `F/F -> 1`, `F^0 -> 1`, and special-function degeneracies are applied only when the required conditions are proved. Exact complex results with zero imaginary part now normalize back to real values.

### Differentiation, integration, and limits

- Generalized exact rational-function integration with Hermite reduction and algebraic `Root` fallbacks, extending support to higher-degree and repeated denominators while preserving simpler `Log/atan/atanh` and elementary forms when available.
- Expanded assumption-aware definite and improper integration, including Gamma/Beta/Mellin, Frullani-type, logarithmic moments, `1+x^q` families, and Dirichlet/Fresnel integrals. Results remain unevaluated when convergence or internal singularities cannot be proved safe.
- Strengthened `D` across `cases[...]` and for higher derivatives. Variable-dependent boundaries no longer receive spurious derivative values, and Lambert W, `polylog`, and quadratic-exponential families gain compact exact higher-order formulas.
- Added `limit[expr,{x,a,direction}]` and expanded known limits for `Ei` / `Ci` / `li`, real-function endpoints, oscillatory cases, and squeeze arguments. Integration derivative-back checks now verify identities on their common domain.

### Algebraic roots and high-degree polynomials

- Removed major performance cliffs in high-degree Complex `Root` and general polynomial `solve`. Root isolation/refinement and irreducibility work are reused more aggressively, including difficult symmetric and multi-scale polynomials; the audited range now extends through degree 96.
- Improved AlgebraicNumber / number-field arithmetic by avoiding unnecessary intermediate `Root` materialization while preserving exactness.

### Certified `N` and special functions

- Greatly expanded certified complex `N`. Principal-branch, enclosure-aware evaluation now covers broader regions for `erf/erfc`, `Ei/Si/Ci`, Fresnel functions, `1F1`, `2F1`, `polylog`, `zeta`, `gamma`, `digamma/trigamma`, Lambert W, and elliptic integrals.
- Removed many artificial fixed-magnitude limits. Real `Ei/Si/Ci`, complex `Ei/Ci`, Fresnel functions, `1F1`, `2F1`, `polylog`, and `ellipticF/E/Pi` now select algorithms from convergence and resource limits instead. Certified `zeta` now covers the finite complex plane except `s=1`, and `ibeta` accepts positive real interval parameters including finite-precision inputs.
- Added certified complex Lambert W near the `-1/e` branch point and for arbitrary integer branches. If finite-precision input does not determine the branch side, mmCal returns `N::precision` instead of guessing.
- Clarified `CertifiedEnclosure` versus `InformationEnclosure` so hidden guard digits from finite-precision inputs cannot be promoted into unsupported output precision. Pole, branch-cut, and algorithm-boundary ambiguity is separated from mathematical DomainError and reported as `N::precision` or `N::unsupported` as appropriate.
- Simplified approximate CLI rendering while retaining approximation provenance. Nested `N`, `diff[...,digits]`, and `nintegrate[...,digits]` now share the same precision and information limits.
- Fixed cancellation and excessive refinement near the `acosh` branch cut, around zero in `expm1/log1p`, and in several `2F1` / complex `li` paths.

### Exact linear algebra

- Added `conditionNumber`, `pseudoInverse`, and `leastSquares`. Exact matrices use exact rank decisions, while outer `N[...]` paths use certified SVD where supported.
- Added a modular + CRT path for exact Integer/Rational `det` and `solveLinear`, selected against Bareiss with exact verification. Also fixed avoidable precision loss in `leastSquares` intermediates.

### Evaluation resources and CLI

- Added shared per-evaluation `EvaluationBudget` / `EvaluationLimits` covering large integers, Arrays/Matrices, Solver and integration work, algebraic construction, and requested precision. Exhaustion is reported as `ResourceLimitError`.
- Added cooperative cancellation with Ctrl-C / Ctrl-Break / SIGINT.
- Added `--batch` for line-oriented automation. REPL `:help` now provides a complete callable catalog, constants, input rules, examples, and nearby-name suggestions.

### Validation and compatibility

- Expanded public-CLI black-box audits and the Random Expression Fuzzer across `D`, `integrate`, `limit`, `N`, `cases`, Gröbner operations, Solve, exact FFT, linear algebra, and precision provenance.
- Fixed regressions in formatting, history references, calculus across `cases[...]`, principal branches, removable singularities, and DomainError classification while preserving the public syntax and exact-first semantics.

## v1.5.3 — 2026-08-16

v1.5.3 connects **exact algebraic number fields, Solver semantic unification, high-precision special functions, exact Cyclotomic FFT, and representation optimizations** on top of the symbolic calculus, certified `N`, and Array/linear-algebra foundation from v1.5.2. The exact-first, proof-only policy remains: unproved or unsupported cases are not promoted to guessed results.

### Algebraic numbers and number fields

- Unified Real/Complex `root[...]` under `AlgebraicNumber` with bounded minimal-polynomial reduction, primitive elements, resultants, and certified root re-identification.
- Added immutable `NumberFieldContext` / `AlgebraicElement` with selected embeddings and Rational power-basis coordinates. Same-field arithmetic, small integer powers, and exact inverse stay inside the field.
- Added reuse caches for fields, reciprocals, and minimal polynomials; minimal polynomials use exact Krylov elimination. Exact algebraic equality and Real ordering are supported.
- Connected `root[...]`, exact Rational/complex Rational, `sqrt` / `cbrt`, `Phi`, and bounded exact arithmetic through a shared algebraic-value layer used by comparison, domain reasoning, and Solve. User-visible representation remains based on `root[minpoly,k]`.

### Solver and semantic coherence

- Added `solve[equation,Integer|Rational|Real|Complex]`, with safe variable inference only when a single unknown user symbol exists; protected symbols are rejected as solve variables.
- Normalized equivalent forms such as `E^x` / `exp[x]` and `ln` / `log`, and added proof-safe logarithmic inversion for positive constant-base exponentials. Added exact symbolic `lambertw` plus certified real branches `k=0/-1`.
- Propagated provable noninteger/non-Rational facts through `ValueFacts`. `N[SolutionSet]` preserves solution-set structure while approximating only numerically closed right-hand sides; comparison formatting was also cleaned up.

### Certified numerical evaluation and special functions

- Promoted `DecimalApproximation` / `ComplexDecimalApproximation` to first-class certified numeric leaves carrying both `CertifiedEnclosure` and `InformationEnclosure`.
- Standardized `N[expr,p]` on significant decimal digits. If whole-expression certification cannot close, numeric subtrees can still be approximated without violating Hold semantics, and duplicate diagnostics are suppressed.
- Added `zeta`, `digamma`, `trigamma`, and regularized `ibeta`, plus high-precision Gamma/Beta optimizations. Existing-but-unsupported certified values use `CertifiedBackendUnsupported` rather than DomainError.

### Exact FFT and Array representation

- Formalized the exact Cyclotomic FFT for supported non-power-of-two exact inputs using Rational coordinates in `Q[t]/Phi_n(t)`. Round-trips such as lengths 5/7/10/12 close exactly without huge `cis[...]` expressions; unsupported/budget-exceeded cases fall back to the generic exact DFT.
- Replaced the large inline `Expr::Node` variant with kind-specific typed nodes while preserving the public `Expr` API.
- Reworked dense `ArrayExpr` into immutable paged packed storage with shape/offset/strides views. Transpose is a zero-copy view, and rectangular numeric braces build directly through `ArrayBuilder`.

### Functions, utilities, and diagnostics

- Added deterministic exact `isprime` / `nextprime` / `prevprime` / `factorint` / `totient` across `uint64`; unsupported BigInt primality is not promoted from probable-prime evidence to `True`.
- Added bit operations, `round[x,n]`, `fma`, `clamp`, `proj`, `range` / `table` / `map`, `explain`, and parameterized periodic Real solution families.
- Parse/lower failures no longer consume `In[n]`; indeterminate forms such as `Infinity-Infinity` and `0*Infinity` are not incorrectly simplified. Unused Machine/double shims were removed.

### Performance and validation

- Added a seed+case reproducible parallel semantic fuzzer to `mmCal.Benchmarks --random-expressions`, with independent-evaluation reset support for long runs.
- Added algebraic-field, special-function, and exact-Cyclotomic-FFT benchmarks and recorded measured optimization decisions in `docs/performance_optimization.*`. Caches/expansions that did not show a reliable gain remain disabled.

### Documentation and compatibility

- Synchronized README, Reference, Architecture, Roadmap, and performance documentation with v1.5.3; intentional limitations and exploration budgets are kept in the roadmap and intentional-limits memorandum.
- Internal primitive elements and cache representation remain hidden from Formatter output; the canonical algebraic display stays `root[minpoly,k]`.

## v1.5.2 — 2026-08-13

### Syntax and REPL (including breaking changes)

- Standardized function calls on `name[...]`; `name(...)` was removed and `()` is grouping-only. Using the old syntax for a known function is a SyntaxError.
- Ordinary `x(x+1)` is accepted as implicit multiplication and normalizes to `x*(x+1)`.
- Standardized history on `In[n]` / `Out[n]`: positive indices are absolute, negative indices relative, zero invalid; `@` forms abbreviate `In[-n]` and `%` forms abbreviate `Out[-n]`.

### Symbolic calculus and integration knowledge

- Added shared finite-Fourier reduction for nonnegative integer powers of `sin[u]^m cos[u]^n`, reciprocal-trigonometric recurrence for negative powers, integer-power reductions for `tan/cot/sec/csc`, and cross-frequency product-to-sum rules.
- Added bounded Weierstrass substitution `t=tan[x/2]` for rational expressions in a common `sin/cos` argument, feeding the transformed result to the exact rational integrator.
- Strengthened inverse-chain/substitution matching and concrete quadratic-radical families, covering patterns such as `2x(1+x^2)^5` and `x/(1+x^4)`.
- Kept derivative-back auditing while introducing `ResolutionOnly` so proof-engine limitations do not discard correct primitives. Integration failures are classified as `unsupported`, `partial`, `conditionsRequired`, or `noKnownClosedForm`.
- Used `docs/memorandum/integralCatalog.md` to audit trigonometric, rational, special-function, and branch-sensitive integration families.

### Special functions

- Added `fresnelc` / `fresnels`, `hypergeometric1F1` / `hypergeometric2F1`, and incomplete `ellipticF` / `ellipticE` / `ellipticPi`, with exact reductions, differentiation, certified real `N`, and corresponding integration support.
- Added `Ei` / `Si` / `Ci` / `li` / `polylog` with representative exact values, derivatives, certified real `N`, and integration rules.
- `integrate[exp[x^n],x]` and related cases prefer branch-safer 1F1 representations where appropriate.
- Solver does not invent nonexistent general inverse special functions; only safe exact degeneracies fall through to existing solving paths.

### Precision-aware `N` and FFT

- Added exact Cyclotomic FFT for supported non-power-of-two exact Rational inputs of length at least 5, operating in `Q[t]/Phi_n(t)`. Exact round-trips avoid huge root-of-unity expressions; cases outside the degree/field proof budget fall back to the generic exact DFT.
- Extended `N[expr,p]` into a precision-aware evaluation entry point that propagates `ApproximationContext`. `N[fft[data],p]` dispatches directly to certified BigFloat/`ComplexInterval` FFT instead of first constructing a huge exact transform.
- Certified non-power-of-two FFT uses direct DFT for smaller sizes and Bluestein for larger ones. FFT and Matrix paths share interval conversion, decimalization, and guard-digit refinement infrastructure.
- Certified approximate display compresses redundant trailing zeros while preserving enclosure metadata; terminating values produced by `N` retain at least one provenance zero to remain visually distinct from exact numbers.

### Arrays and linear algebra

- Generalized `{...}` into a finite brace container; rectangular children are optimized into dense row-major `ArrayExpr`. Zero-length dimensions and empty-array shapes are preserved.
- Regularized `dimensions`, `arrayRank`, `length`, `at`, and `reshape`; added `MatrixView` / `MatrixBuffer` to reduce copying. Canonical APIs are `dot`, `matrixRank`, `norm`, and `normalize`, with legacy aliases retained.
- Added Bareiss fraction-free elimination for Integer/Rational matrices, shared by `det`, `rref`, `matrixRank`, `inverse`, `solveLinear`, and `nullSpace`. `solveLinear` returns unique solutions only; `nullSpace` has deterministic bases and preserves empty shape.
- Added `luDecomposition[A]`: exact row-pivoted `P A = L U`, with certified interval partial pivoting under `N[...]`.
- Added rectangular reduced Householder `qrDecomposition[A]`. General exact QR is limited to 3x3 after measured expression growth; triangular/trapezoidal fast paths remain. Automatic blocking was benchmarked but not enabled because no block size won consistently.
- Replaced the general exact `qrDecomposition[A]` path from Expr-level Householder expansion with fraction-free orthogonalization. Orthogonalization creates neither square roots nor divisions, keeps primitive integer vectors with GCD content reduction, and uses a symmetric Bareiss Gram factorization (fraction-free LDLᵀ-equivalent) as the full-rank fast path. The former 3x3 hard cap is removed; rank-deficient inputs fall back to direct fraction-free orthogonalization.
- Re-audited exact Matrix performance cliffs and replaced the post-Bareiss generic Rational RREF used by `inverse` / unique `solveLinear` with common-denominator BigInt back substitution. Full-column-rank `rref` also skips the Rational backward phase. Modular inverse was remeasured through 32x32 / 256-bit and remains slower than Bareiss, so automatic dispatch stays disabled.
- Added real/complex reduced `svd` without forming `A^H A`, using Householder bidiagonalization plus one-sided Jacobi; added `conjugateTranspose`.
- Added `eigenvalues` / `eigenvectors` / `eigensystem`: simple exact matrices stay exact, while general `N[...]` uses Complex BigFloat Hessenberg reduction, shifted QR, and Schur processing.
- Symbolic `det` / `inverse` use triangular fast paths and expansion budgets. Discontinuous rank/null-space decisions require exact pivots or interval certification rather than arbitrary epsilon thresholds.

### Performance, stability, and development infrastructure

- Bareiss measured about 5.4–11.4x faster for order-8–16 determinants and 7.7–15.5x for RREF. Added `--matrix-large` for order-32/64 timings and order-1024 storage/parse stress.
- Large 1024x1024 Rational matrices showed Expr/Rational representation as the main memory/parse bottleneck, motivating typed `Expr::Node`, BigUInt/BigInt SBO, and packed Array work for later versions.
- Fixed MSVC-sensitive negative `RealInterval::point` construction and a transitive-include dependency on `<stdexcept>`, and expanded fixed-seed Matrix/FFT invariants.

### Documentation and licensing

- Synchronized README, Reference, Architecture, Roadmap, and performance documentation with v1.5.2.
- Aligned BSD 3-Clause copyright metadata and documented trademark/brand use separately in `TRADEMARKS.md` / `TRADEMARKS.ja.md`.

## v1.5.1 — 2026-08-12

v1.5.1 retains the exact-first CAS foundation from v1.5.0 while focusing on canonicalization, validation, large integers, high-precision numerical evaluation, and benchmark infrastructure.

### Symbolic / CAS

- Added deterministic total ordering for `Add`, definedness-preserving product/division normal forms, and stronger `MathKnowledge` nonzero/reversed-relation inference.
- Added/expanded formatter-parser round-trip tests, integration derivative-back checks, Reference-to-builtin consistency checks, and extreme-value approximation tests.
- Fixed radix-prefix, Array-adjacency, and negative-Rational formatting issues, and added FFT plan/twiddle caching across transforms.

### Multiprecision / high precision

- Added adaptive schoolbook / Karatsuba / Toom-3 BigUInt multiplication, dedicated squaring, Burnikel-Ziegler division, power-of-two division fast paths, and divide-and-conquer decimal conversion.
- Improved balanced-product-tree factorial, oversized `tryToUint64` rejection, and BigFloat addition/subtraction for extreme exponent gaps.
- Switched `Pi` to binary-splitting Chudnovsky, `exp/E` and `log` to binary-splitting certified evaluation with range reduction, and added certified argument reduction for huge-radian `sin/cos/tan`.

### Benchmark / test

- Added `mmCal.Benchmarks` to the Visual Studio solution with fixed-seed random invariants, algorithm-threshold sweeps, factorial/decimal-I/O benchmarks, and high-precision `Pi/exp/log` benchmarks.
- v1.5.1 release state: internal 1691 / 1691, black-box 1337 / 1337.

### Benchmarked but not selected

- Prime-Swing factorial, binary GCD, Karatsuba workspace variants, dedicated/low-threshold Toom-3 squaring, and machine-`fmod` huge-trig reduction were not adopted because measurements did not show sufficient benefit.

See `docs/performance_optimization.*` for the measurement details behind these decisions.

---

## v1.5.0

Rebuilt most of the numerical model, Lexer/Parser/AST/Evaluator, Simplifier, Solver, CertifiedEvaluator, CLI, Formatter, tests, and documentation, redefining mmCal as an exact-first CLI calculator / compact CAS.
