# Changelog

## Unreleased

### Semantic hardening and build

- Normalized exact complex arithmetic back to real `Number` storage whenever the imaginary component becomes zero. Also supplied the missing certified finite-precision elliptic-integral and complex `li` definitions, restoring successful CMake/GCC linking for all targets.
- Added a certification-boundary fuzzer to `mmCal.Benchmarks`. It generates 46 families covering branch cuts, exact poles, finite-precision boundary crossings, real/complex backend boundaries, and work limits; classifies `Value`, `DomainError`, `N::precision`, `N::unsupported`, unevaluated, and resource failures separately; and is reproducible through `--seed --case`. Timeout/cancellation failures are kept distinct from mathematical failures, and the Visual Studio/CMake source sets are synchronized.
- Made `simplify` / `fullSimplify` definedness-aware. Reductions such as `F-F -> 0`, `F/F -> 1`, `F^0 -> 1`, and special-function degeneracies are applied only when the required domain conditions are provable. This includes the `0^0` convention, negative Rational powers, `zeta[1]`, and Gamma-family poles.
- Changed integration derivative-back regression checks to prove identities only on the common domain. The Random Expression Fuzzer now also covers `limit`, `cases`, nested `N`, AlgebraicNumber, Gröbner reduction, Array reshape, domain-hole preservation, and removable singularities.
- Standardized the MSVC stack reserve at 16 MiB for the CLI, tests, and benchmarks. Shared builtin arity/unevaluated-call handling and Rational rounding helpers were consolidated. Exact complex `digamma` / `trigamma` also gain an exact recurrence path for positive-integer real parts to avoid unnecessary interval cancellation.

### Certified `N` and precision provenance

- Avoided duplicate certified-function work for exact `N` inputs, where CertifiedEnclosure and InformationEnclosure are identical. Internal guard precision is also no longer compounded unnecessarily across complex `li` (`Log -> Ei`) and the principal `2F1` `1/z` connection, preserving the same principal values and interval guarantees while substantially reducing representative 20-digit runtimes.
- Fixed `leastSquares` losing an extra output digit because its intermediate pseudoinverse was rounded to the requested display precision before the final product. The intermediate now retains working precision and only the final result is rounded. Finite-precision `nullSpace` regressions now also require both pivot existence and free-column status to be certified from InformationEnclosure rather than recovering hidden exact-zero truth.
- Fixed an `acosh` performance cliff immediately above and below the interior branch cut `-1<x<1`, where tiny imaginary components could trigger excessive refinement and frontend timeouts. A stable path is used only after the half-plane is certified, so cases as small as `N[acosh[1/2+I/10^160],100]` retain the principal branch without pathological refinement.
- Fixed cancellation in `expm1` / `log1p` near zero by deriving extra working precision from the binary scale of the input instead of subtracting at the requested precision. Exact real `2F1` inputs with `z>1` now use the defined principal-cut continuation value, while finite-precision inputs spanning both sides still return `N::precision`. Complex Fresnel series evaluation now has a bounded-work limit of `|z|<8`.
- Preserved the historical `diff[...,digits]` / `nintegrate[...,digits]` contract in which `digits` counts fractional decimal places. Finite-input `InformationEnclosure` limits still cap claimed precision, but the display mode is no longer accidentally changed to the significant-digit rule used by `N`.
- Separated singularity, branch-cut, and algorithm-boundary decisions between `CertifiedEnclosure` and `InformationEnclosure`. Exact poles yield DomainError; ambiguity caused by finite-precision input information yields `N::precision`; mathematically valid values outside the implemented certified algorithms yield `N::unsupported`. The same rule now covers the Gamma family, `Ei` / `Ci` / `li`, `zeta`, `1F1` / `2F1`, `ibeta`, `log` / `sqrt` / noninteger powers, inverse trigonometric and hyperbolic functions, `Arg` / `atan2`, `polylog`, and elliptic integrals.
- Refactored the certified numerical layer into value representation, zero/pole/branch-cut classification, Real/Complex arithmetic, and precision helpers. This is an internal separation of responsibilities with no intended mathematical semantic change.
- Clarified the contract between `CertifiedEnclosure` and `InformationEnclosure` in `DecimalApproximation` / `ComplexDecimalApproximation`. Proven exact zero is no longer automatically treated as reusable exact input information. `precision`, `accuracy`, `explain`, comparisons, and certified arithmetic now share the same information-quality model; `explain` adds `PrecisionDigits` / `AccuracyDigits`.
- Audited hidden-guard reuse across finite-precision values. Zero/nonzero decisions such as `N[0,p]^0` and `1/N[0,p]`, display capping after complex-to-real projection, removable cardinal singularities, FFT cancellation, and matrix pivot/rank decisions now respect the InformationEnclosure. Continuous quantities propagate Certified/Information enclosures in parallel. `lu/qr/svd/conditionNumber/pseudoInverse/leastSquares/eigen*` remain conservative for matrices that already contain finite-precision leaves until their input-perturbation certificates are complete, rather than recovering results from hidden certified points.
- Bounded top-level certified refinement in `N` to 16 local attempts. Finite-precision inputs that continue to straddle a branch cut, pole, or algorithm boundary return `N::precision` instead of consuming unbounded guard digits, while exact inputs may still refine when additional precision can prove the correct side.
- Extended InformationEnclosure branch gating to inverse trigonometric/hyperbolic functions, `Arg` / `atan2`, `polylog`, and elliptic `F/E/Pi`. Complex inverse functions use a derivative-bound enclosure away from branch points to reduce dependency blow-up without inventing precision.
- Simplified CLI rendering of certified approximations while keeping precision metadata internal. Near-zero non-point values render compactly as `0.0`; trailing zeros are shortened; and `N`-produced terminating values retain a provenance zero, e.g. `N[1/2,20] -> 0.50`, `N[2,20] -> 2.0`, `N[I,20] -> 1.0I`, and `N[log[-1+I/10^1000],20] -> 0.0+3.1415926535897932385I`. Neutral exact operations preserve metadata, and zero-centered formatting uses absolute accuracy rather than relative precision.
- Aligned real/complex certified paths, including safe real projection for all-real `1F1` / `2F1` / `polylog`, `RealInterval` elliptic evaluation, and bounded interval evaluation for `2F1` near denominator-parameter poles. Numerical `2F1` with exact Rational parameters now accumulates the Gauss series in outward-rounded intervals instead of growing huge exact Rational numerators and denominators. `diff` / `nintegrate` propagate `CertifiedEnclosure` and `InformationEnclosure` in parallel, including held `N[...]` subexpressions, so they cannot manufacture output precision beyond the input information. Persistent singularity or branch-side ambiguity from finite input information fails locally instead of triggering pointless refinement. Dedicated N black-box/property regressions were added.

### Black-box validation

- Added public-CLI-only audits for `D`, `integrate`, `limit`, `N`, `cases`, Gröbner expressions, and exact algebraic Solve, plus property suites for derivative-back identities, Gröbner invariants, solution counts, precision provenance, recurrences, and independently generated high-precision references.
- Fixed calculus across `cases[...]`: `D`, `integrate`, and `limit` distribute branch-wise only when conditions are independent of the calculus variable. Variable-dependent or singular branches remain unevaluated or are handled by the correct limit rather than producing spurious DomainErrors.
- Allowed safe composition of `groebnerBasis[...]` inside held polynomial builtins, and updated pre-`cases` black-box expectations/reference examples to the current scalar-case representation.

### Scalar cases, polynomial ideals, and complex polygamma

- Added first-class scalar `cases[value if condition; ...]`, distinct from control-flow `if[...]` and `SolutionSet`. Unknown conditions are retained, false branches are not evaluated, and the construct is integrated with `simplify`, `N`, `D`, `integrate`, and `limit`.
- Added a general multivariate polynomial core over `Q[x1,...,xn]` with Lex / GrLex / GrevLex orderings, division/normal forms, S-polynomials, and Buchberger reduction. Exact Rational Gröbner bases are exposed through `groebnerBasis[...]` / `polynomialReduce[...]` under existing evaluation budgets.
- Connected nonlinear polynomial `solve[{...},{...}]` to Lex Gröbner elimination. Contradictory ideals return the empty set; supported zero-dimensional systems use exact univariate roots and back-substitution with exact verification; positive-dimensional cases remain `UnresolvedSolutionSet` rather than inventing parametrizations.
- Added certified complex `N` for `digamma` / `trigamma` using recurrence and Bernoulli/Stirling asymptotics. General `polygamma[n,x]` remains future work. The Visual Studio project was synchronized with the new polynomial sources.

### Symbolic and formatting fixes

- Improved conditional-solution `&&` spacing, nested-iterator `table`, and exact Complex Rational-imaginary formatting (`2I/29`).
- Added exact perfect-power recognition to univariate polynomial `factor`, and fixed HoldAll `D` so `%` / `Out[n]` resolves the stored output before differentiation without changing `In[n]` re-evaluation semantics.
- Added the general symbolic power rule for `integrate[x^n,x]`, `integrate[log[log[x]],x]`, and a branch-safe primitive for `integrate[li[x],x]`. The convergent improper integral `integrate[log[log[x]],{x,1,E}]` is now handled exactly.
- Fixed Formatter radix-prefix collision handling at lexical boundaries, avoiding unnecessary `*` in forms such as `380x^9`. Univariate polynomials are displayed in descending degree without changing internal Expr ordering.
- Added `limit[expr,{x,a,direction}]`, expanded known principal-branch limits for `Ei` / `Ci` / `li`, and simplified `li[0] -> 0`. Periodic oscillation of real `sin` / `cos` / `tan` now returns `Indeterminate` when no single limit exists, while a bounded squeeze rule closes cases such as `limit[x sin[1/x],x,0] -> 0` exactly.

### `D`, integration, and certified complex `N` audit

- Audited derivative formulas, primitives, special-function identities, principal branches, and removable singularities. Rules were added or strengthened for Beta-family functions, finite combinatorial functions, `Ei/Si/Ci/li/digamma/trigamma/LambertW`, `polylog`, `1F1/2F1`, `ibeta`, and `sinc/cosc/expc` families.
- Derivatives of `sinc/cosc/tanc/sinhc/tanhc/expc` now retain continuous-extension values at zero through `cases[...]`; `Si'`, principal `LambertW`, and `polylog` likewise preserve finite zero-point derivatives where appropriate.
- Extended certified complex `N` on `ComplexInterval` to `erf/erfc`, `Ei/Si/Ci`, Fresnel functions, `1F1`, `2F1`, `polylog`, `zeta`, and `gamma` using convergent series, Euler-Maclaurin, Stirling, and explicit remainder bounds. `2F1` uses the principal `1/z` connection formula only when its safety conditions are provable.
- Negative noninteger Rational `digamma/trigamma` and negative Rational `zeta` are routed through exact recurrence/functional equations. `lgamma` remains real-axis `log|Gamma|`. Derivative-back auditing avoids promoting local principal-branch identities into unsafe global identities.

### Certified backend and algorithm thresholds

- Clarified `PrecisionInsufficient` versus `CertifiedBackendUnsupported`: fixed series ranges, term limits, and planning limits that cannot be resolved by extra precision now return `N::unsupported`; interval-width ambiguity at branch/work boundaries may still refine.
- Retuned special-function series boundaries and cancellation planning. Real `Ei/Si/Ci` remain bounded at 96, complex `Ei` at `|z|<=512`, and complex `Ci` at `|z|<=128`; `Si/Ci` and complex `Ei/Ci` now add argument-dependent working precision for cancellation. The existing `1F1: |z|<=160`, `2F1: |z|<=9/10`, elliptic `F/E: |m|<=9/10` (plus `|n|<=9/10` for `Pi`), and positive-order `polylog: |z|<=49/50` limits remain.
- Corrected `zeta` classification so only `s=1` is a DomainError; regions unsupported by the present Euler-Maclaurin backend are no longer mislabeled as mathematically undefined.
- Bounded certified FFT precision retries to 12 and kept the direct-DFT/Bluestein policy boundary at 384 points. Regression tests now explicitly cross exact-linear-algebra and BigUInt Karatsuba / Toom-3 / Burnikel-Ziegler / decimal-conversion thresholds.

### Exact linear algebra and cancellation

- Added `conditionNumber`, `pseudoInverse`, and `leastSquares`. Exact matrices use exact rank decisions; the pseudoinverse uses rank factorization to construct the Moore-Penrose inverse exactly, including rank-deficient Rational and complex matrices. The outer-`N` path for exact input uses certified SVD. Matrices that already contain finite-precision leaves are currently kept conservative rather than recovering rank or singular subspaces from hidden guard digits. Zero-sized matrix shapes are preserved.
- Added a 31-bit prime-field + CRT modular path for exact Integer/Rational matrices. `det` reconstructs to a Hadamard-bound certificate; `solveLinear` uses rational reconstruction followed by exact `A X = B` verification; `inverse` verifies the reconstructed adjugate and falls back to Bareiss when needed.
- Added measured Bareiss/modular dispatch for `det` / `solveLinear`. `inverse`, `rref`, `matrixRank`, and `nullSpace` remain on Bareiss for now. Modular-prime telemetry and `--exact-linear-algebra` / `--budget-telemetry` benchmarks were added.
- Connected `EvaluationCancellationToken` to top-level CLI evaluation. Ctrl-C / Ctrl-Break on Windows and SIGINT on POSIX request cooperative cancellation, reported as `ResourceLimitError`; noninteractive exit code `3` is preserved.

### Evaluation resource policy

- Added per-top-level-evaluation `EvaluationBudget` / `EvaluationLimits` / `EvaluationUsage`, covering evaluation steps/depth, generated Expr nodes, Simplifier/Solve/integration candidates, certified refinements, Array/Matrix elements, BigInt size, requested precision, and algebraic construction work.
- Added `KernelSession::setEvaluationLimits`, `evaluationLimits`, and `lastEvaluationUsage`; legacy `setEvaluationDepthLimit` maps into the common limits. Resource exhaustion is reported as `ResourceLimitError` with the resource name and limit, separate from DomainError and ordinary unevaluated results.
- Large integer powers/factorials, Arrays, Matrix workspaces, `N[expr,p]`, and source text are checked as early as practical before expensive allocation or computation. Core remains free of wall-clock deadlines; deterministic operation budgets are separate from frontend cancellation/timeout policy.

### Fuzzing and validation

- Extended the Random Expression Fuzzer with derivative-back, Solve, `A inverse[A] == I`, `ifft[fft[v]] == v`, DomainError classification, and principal-`sqrt` boundary checks. Polynomial identities use structural residuals plus independent substitutions, and Solve compares `SolutionSet` bindings as sets.
- Exact FFT round-trips first use structural equality; remaining root-of-unity cancellation forms are passed directly to `CertifiedEvaluator` as a second oracle. Deep cases stop under a recursion-depth budget, and indeterminate cases are counted separately as `inconclusive` rather than mathematical failures.
- Failure output now includes seed/case/reduced expression and budget telemetry; regression coverage was added for limit crossings, reset behavior, and success/failure telemetry.

### CLI and compatibility

- Added `--batch` for line-oriented automation.
- Kept startup `--help` concise while expanding REPL `:help` into a complete callable catalog with descriptions, input rules, and examples. Added `:help Pi` / `:help constants` and deterministic nearby-name suggestions; help remains frontend-only and does not consume evaluation/history slots.
- Added `mmCal.Benchmarks --fft-threshold [iterations]` to compare direct DFT and forced Bluestein from 65–509 points. GCC/MSVC measurements retuned the certified non-power-of-two FFT policy boundary to 384 points, retained by forced-comparison tests.

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
