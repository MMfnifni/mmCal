# Changelog

## v1.5.5 — Unreleased

### Series and asymptotics

- Added `series[...]` on a `SeriesData` / TPSA core. Exact Taylor, Laurent, Puiseux, and logarithmic coefficient layers share one representation and connect directly to `D`, `integrate`, and `normal`. In addition to elementary functions, supported local providers cover `erf` / `erfc` / `Si` / `Ei` / `Ci` / `li`, Fresnel functions, principal inverse trigonometric functions, Lambert W, `gamma` / `lgamma` / `digamma` / `trigamma`, and supported `polylog` neighborhoods using exact coefficient recurrences and DLMF formulas. `tan/cot/sec/csc`, hyperbolic variants, `expm1/log1p`, cardinal functions, and `log2/log10` lower into the same TPSA algebra.
- Added `series[expr,{x,Infinity,n}]`, where `Infinity` denotes real `+Infinity` and is mapped to the local variable `t=1/x`. Rational functions, polynomial growth, reciprocal composition, Puiseux terms, and logarithmic asymptotics reuse the existing Series algebra. When a positive leading direction is proved, `log[A(x)]` is factored as `log[c]+r log[1/x]+log[1+h]`. Oscillatory forms, essential growth, negative logarithmic powers, and other transseries that need a richer representation remain unevaluated rather than guessed.
- Added a limited finite-point `limit` fallback that uses local Series only when the existing limit kernel is unresolved, including exact cancellation of singular terms. `toNormal[expr]` recursively normalizes `SeriesData` inside lists, arrays, and calls; for finite and conditional `SolutionSet` values it transforms only binding right-hand sides while preserving conditions, free variables, multiplicity, and domains.
- Fixed `N[SeriesData]` incorrectly approximating structural exponent-grid metadata and thereby breaking later `normal` / `D`. Only coefficients and the center are approximated; `minimumExponent`, `orderNumerator`, and `exponentDenominator` remain exact integers.
- Fixed false `Division by zero` failures in positive-infinity reciprocal composition such as `tan[1/x]`, `sec[1/x]`, `cot[1/x]`, and `csc[1/x]`. Reciprocal normalization is restricted to the local-variable substitution path and does not change global Simplifier semantics.
- Reduced a performance cliff in regular `gamma` Series expansion. The `lgamma`-to-exponential coefficient recurrence now builds each raw coefficient sum first and simplifies once per coefficient instead of repeatedly simplifying every intermediate term, preserving low-order canonical output while substantially reducing expression growth.
- Added coefficient recurrences for `Ei` / `li` Series with exact Rational affine arguments, avoiding generic Series products, reciprocals, and repeated simplification. Puiseux compositions such as `li[2+sqrt[x]]` reuse the affine form in the local variable instead of rebuilding large intermediate expressions.

### Arrays, vectors, and vector calculus

- Added `Array / scalar` as elementwise scalar scaling. Division by exact zero preserves scalar exceptional-value semantics componentwise, so `{0,1}/0 -> {Indeterminate, ComplexInfinity}`. `scalar / Array` and `Array / Array` remain rejected rather than inventing an implicit linear-algebra meaning.
- Added `inner[a,b]`, `outer[a,b]`, `distance[a,b]`, `manhattanDistance[a,b]`, `projection[a,b]`, and `rejection[a,b]`. `dot` remains bilinear, while `inner` conjugates its first argument. Norm, distance, projection, rejection, and reflection now share the same Hermitian semantics.
- Added `orthogonalQ`, `orthonormalQ`, `linearIndependentQ`, and `gramSchmidt` for row-vector sets stored in rank-2 Arrays. Orthogonalization is exact-first and keeps an unnormalized orthogonal basis internally to avoid radical-heavy intermediate growth; provably dependent rows are dropped without guessing unresolved symbolic zero/nonzero decisions.
- Removed the old vector convenience spellings `vadd` / `vsub` / `vscalar` / `vsum`; Array `+` / `-` / scalar multiplication and `sum` are now the direct API. Also removed `vcross` / `vmanhattan` / `veuclidean` / `vproject` / `vangle` / `vreflect` / `vreflect_axis` in favor of the primary names `cross` / `manhattanDistance` / `distance` / `projection` / `vectorAngle` / `reflectNormal` / `reflectAxis`.
- Added `vectorAngle[a,b,assumptions]` so symbolic realness can be stated explicitly. Geometric angle remains restricted to real vectors rather than being guessed for complex vectors.
- Added Cartesian `grad`, `divergence`, `curl`, `laplacian`, `jacobian`, `hessian`, and `directionalDerivative`. Two-dimensional `curl` returns the scalar curl, three-dimensional `curl` returns the vector curl, and `laplacian` now applies componentwise to vector fields. Cylindrical/spherical metric factors are never inferred implicitly.
- Extended `at` to finite `SolutionSet` values. `at[solutions,i]` returns the zero-based selected branch while preserving conditions, free variables, multiplicity, and solver-variable domains; `at[solutions,i,x]` returns the requested binding value. Non-finite sets, out-of-range indices, and non-solver selectors are rejected explicitly.

### Symbolic evaluation, CLI, and FFT

- Principal `sqrt` equations now keep a complete solution set when the right-hand side is a generic complex parameter. `sqrt[x]==y` retains `x==y^2` together with the principal-square-root range condition `re[y] > 0` or `re[y] == 0 && im[y] >= 0` on its solution branches. Exact complex right-hand sides are reduced with bounded exact arithmetic before range testing so values such as `-2+3I` do not survive as false conditional branches.
- Real `solve` now inverts principal `asin` / `acos` / `atan` / `acosh` equations with their exact output-range conditions. Inverse-trigonometric ranges follow the active angle mode and distinguish open from closed endpoints. Symbolic `log2[x]==y` / `log10[x]==y` inversion now retains `y in Real` as a branch condition.
- Strengthened exact endpoint handling for principal inverse functions. Closed endpoints such as `solve[asin[x]==-Pi/2,x,Real]` no longer retain tautological range conditions, while open endpoints are rejected before evaluating the inverse transform itself.
- Certified complex `Ci` now rejects asymptotic routes that cannot reach the requested precision before doing expensive `E1` work, avoiding repeated guard retries near backend thresholds such as `N[Ci[140+I],100]`.
- Exact-input `N[...]` no longer evaluates the same expensive special-function backend twice merely to construct `CertifiedEnclosure` and `InformationEnclosure`; without finite-precision provenance those input enclosures are identical. Expressions containing finite-precision values, nested `N`, or external bindings still propagate the two enclosures separately. Exact-real `2F1` sufficiently close to `z=1` now uses a guarded DLMF 15.8.4 `1-z` connection, and exact-Rational connection coefficients are sent directly through the Rational Gamma / LogGamma backend, avoiding both near-unit Gauss-series cliffs and generic interval-Gamma overhead.
- Tightened certified real dispatch for several special functions. Exact-Rational `zeta` / `digamma` / `trigamma` preserve the original rational argument instead of evaluating both dyadic enclosure endpoints, finite-precision negative-real `digamma` shifts to the positive-real backend through a real recurrence, and exact half-integer-`Pi` amplitudes in `ellipticF` / `ellipticE` / `ellipticPi` reuse one complete Carlson integral rather than recomputing it through period reduction.
- Reduced major algebraic-log rational-integration cliffs by keeping residues `P(r)/Q'(r)` as exact expressions over the same `Root` instead of repeatedly canonicalizing intermediate AlgebraicNumbers. Added a low-cost factorization for `x^(4m)+x^(2m)+1` families to expose simpler factors earlier.
- Higher derivatives `D[exp[q(x)],{x,n}]` with quadratic-or-lower `q` now use a polynomial-coefficient recurrence, and nested `D` materializes the inner derivative before continuing. Polynomial products with linear `exp/sin/cos/sinh/cosh` are integrated by finite coefficient recurrences rather than deep repeated integration by parts. Simple mixed-frequency `sin/cos` products are sent through product-to-sum before generic candidate search; when this provably yields a shallower expression tree, flat primitives such as `-cos[2x]/4` are preferred over equivalent power forms. Higher-degree pure-binomial radicals / reciprocals are dispatched early to their existing `2F1` primitives, while `c/log[a*x+b]` connects directly to `li` whenever the affine chain is proved, avoiding generic candidate-search overhead.
- Real polynomial `solve` now shares real-root isolation across proper factors of reducible polynomials. Rationally factorable `x^(2m)+b x^m+c` forms also use a low-cost quadratic-in-`x^m` factorization before entering the existing factored fast path. When all real roots are Rational they are returned directly, preserving canonical forms such as `solve[x^4-1==0,x,Real] -> {x == -1, x == 1}` while retaining minimal-polynomial `root[...]` output for genuinely algebraic factors. Direct bindings are domain-filtered before general algebraic proof work, avoiding needless AlgebraicNumber construction in cases such as `solve[x==Phi,x,Real|Rational]`.
- Arithmetic between distinct pure-quadratic `Root` values now constructs the known degree-2/4 annihilating polynomial and exact conjugate index directly, bypassing general resultant / primitive-element construction while preserving canonical `root[...]` output for forms such as `sqrt[2]±sqrt[3]`.
- Fixed `D[cases[...]]` so explicit default branches do not invent derivative values at boundaries, compacted additive output from higher derivatives, and extended exact-rational linear-combination and rational-affine assumption handling in `simplify` / `fullSimplify`.
- Fixed cross-frontend composition around held symbolic operators. `limit` now protects inner binders/control variables in `D`, vector calculus, `solve`, `collect` / Gröbner operations, and `Series` / `Normal`, materializing supported inner calls before outer point substitution. Two-sided limits of variable-dependent `cases` evaluate the two directions separately and return a value only when they agree; disagreement is `Indeterminate`. Brace/List limits are handled componentwise.
- Fixed Solver definedness and dispatch asymmetries. Identical-expression equations can reduce to `All` while retaining finite representable definedness predicates, relation arrays reuse the same transcendental/radical dispatch as scalar relations, and the Real ambient domain is available during normalization. This connects forms such as `sqrt[x^2]==±x`, `abs[x]==±x`, `exp[log[x]]==x`, and `log[exp[x]]==x` through existing branch/domain knowledge.
- Regularized safe principal-inverse composition. Finite `sin[asin[x]]`, `cos[acos[x]]`, `sinh[asinh[x]]`, and `cosh[acosh[x]]` reduce in the function-after-inverse direction, while `tan[atan[x]]` and `tanh[atanh[x]]` retain the definedness needed to exclude exceptional points. Reverse compositions such as `asin[sin[x]]` remain unsimplified in general. Assumptions such as `x>1` now propagate sign/nonzero facts for `log/log2/log10`, and zero-direction errors from `rejection`, `reflectNormal`, and `reflectAxis` report the public function name rather than an internal `projection` delegate.
- Connected held symbolic frontends across calculus. `integrate[D[...]]`, `Series[D[...]]`, `solve[D[...]==...]`, `integrate[Normal[Series[...]]]`, and nested `Normal[Series]` pipelines now use narrowly scoped safe materialization instead of remaining unevaluated. The same safe materialization works when those frontends are nested inside normally evaluated expressions, Vector Calculus operators, and explicit `simplify` / `fullSimplify` / `expand` / `factor` / `collect` transforms. Transform kernels consume held arguments directly instead of resolving session definitions, while arbitrary held calls are still not eagerly evaluated, preserving binder semantics.
- Strengthened composition of nested binder frontends. When outer `D` / `integrate` / `series` / `solve` / vector-calculus operators materialize inner `limit` / `integrate`, outer binder variables are temporarily protected from current session definitions. `D[integrate[f,x],x]` still prefers the existing fundamental-theorem rule instead of expanding an antiderivative first. Residual `integrate[0,x]` and inner limits over different variables after limit substitution are materialized only when safe.
- Added REPL `:quit` / `:exit` commands and interactive `:layout` composition.
- For power-of-two exact `ifft[fft[v]]` within the current degree budget, FFT-generated root-of-unity expressions are re-embedded into cyclotomic coordinates to reduce inverse-expression growth without changing the forward representation.
- Extended cyclotomic re-embedding for non-power-of-two exact FFTs by recognizing `sqrt[3] = zeta_12 + zeta_12^(-1)` inside `Q(zeta_12)`. Six-point Gaussian-integer `ifft[fft[v]]` no longer falls back to generic Expr simplification and now closes its exact round trip in milliseconds instead of roughly one second.
- Random Expression Fuzzer single-case replay now prints depth, invariant, and generated expression even on PASS, making slow or timeout-prone cases directly diagnosable without first turning them into failures.

### Internal refactoring

- Consolidated exact-scalar AST construction and predicates, integer powers of `Rational`, relation-operator `BuiltinId` mapping, and `Cases/CaseBranch` construction into shared helpers. Series, integration, limits, Solver, and Assumption handling now share the same exact-value and relation semantics, reducing synchronization risks when these facilities are extended.
- Consolidated exact and approximate linear-algebra internals. Exact `NumberMatrix` storage and matrix-operation wrappers are shared, while SVD and eigen code use one NearestEven BigFloat point-arithmetic layer for shape handling, rounding mode, and elementary operations.
- Centralized Builtin-family predicates for Vector, Matrix, statistics, and related groups. Algorithm-specific certified retry loops and small local aliases remain explicit where abstraction would obscure control boundaries rather than remove semantic duplication.
- Added Japanese implementation comments around non-obvious certified retry/fallback decisions, including why `PrecisionInsufficient` may be retried with more guard digits while `CertifiedBackendUnsupported` must fall back immediately. Removed dead integration helpers left behind by earlier algebraic-integration paths. Public APIs, mathematical semantics, and Formatter contracts are unchanged.

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

- The internal test runner reports per-suite elapsed time only when invoked with `--timings`, preserving normal test output by default.
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
