# mmCal Performance and Algorithm Selection Notes

This document records **which performance techniques were adopted, which were benchmarked and rejected, and why**.

Performance is never allowed to override mmCal's core contracts: exact arithmetic, branch/domain correctness, directed rounding, certified enclosures, and reparsable canonical output.

Thresholds are environment-dependent. The measurements below are representative results from the development GCC environment, not universal guarantees. Re-run `mmCal.Benchmarks` on MSVC or other hardware before retuning thresholds.

---

## 1. Selection policy

An optimization is selected only when it:

1. preserves the old mathematical semantics;
2. passes fixed-seed randomized/invariant checks;
3. remains correct at boundaries and extreme values;
4. provides a meaningful measured benefit;
5. can isolate small-size regressions behind a threshold; and
6. justifies its implementation complexity.

An algorithm is not selected merely because it has a better asymptotic bound or is common in other libraries.

---

# 2. BigInt multiplication

## Schoolbook — retained

The simple double loop remains the lowest-overhead base case.

## Karatsuba — adopted

Threshold sweeps showed regressions when Karatsuba started at 8–16 limbs. Gains became stable around 32–48 limbs, so v1.5.1 uses a crossover near **48 32-bit limbs**. Strongly unbalanced operands fall back to schoolbook multiplication.

Representative results:

| operand | schoolbook | adaptive | speedup |
|---:|---:|---:|---:|
| 128 limbs | ~9.5 µs | ~6.9 µs | ~1.4x |
| 256 limbs | ~38 µs | ~21 µs | ~1.8x |
| 512 limbs | ~151 µs | ~63 µs | ~2.4x |

The same backend accelerated the existing balanced factorial product tree.

## Toom-3 — adopted

For still larger balanced operands, Toom-3 becomes beneficial. Representative v1.5.1 thresholds are approximately:

```text
top-level Toom-3  ~1280 limbs
recursive Toom-3  ~448 limbs
```

Starting Toom-3 too early loses to Karatsuba because of evaluation/interpolation overhead. Representative gains over Karatsuba-only were about 1.17x at 4096 limbs and about 1.3x at 6144 limbs.

## Toom-4 / FFT / NTT — deferred

These remain candidates for much larger operands. They are not added until a useful crossover is demonstrated by the benchmark project.

---

# 3. Dedicated squaring

`x*x` uses symmetry rather than the generic multiplication path:

- symmetric schoolbook square for small values;
- Karatsuba square for larger values.

Representative results:

```text
512 limbs   ~63 µs generic multiply → ~38 µs square
1024 limbs  ~189 µs                 → ~111 µs
```

This directly benefits exponentiation by squaring.

### Dedicated Toom-3 square — rejected

It was implemented and benchmarked but remained slower than the Karatsuba square in the relevant size range. Evaluation/interpolation overhead moved its crossover too far upward, so it is not a default v1.5.1 path.

---

# 4. Karatsuba workspace reuse — rejected

Two allocation-reuse designs were benchmarked:

- vector-pool workspace;
- recursion-depth scratch buffers.

Both regressed by up to roughly 5–10% around 512–1024 limbs in the GCC environment. Pool/scratch management, resize work, and cache effects outweighed allocator savings.

This may be worth retesting under MSVC or a different allocator.

---

# 5. Factorial

## Balanced product tree — retained

The existing product tree combines similarly sized integers and works well with adaptive multiplication. Leaves were improved by constructing `BigInt` directly from unsigned machine integers and using one-limb multiplication fast paths.

## Prime-Swing — rejected for the current backend

Prime-Swing was implemented, including fixes for unnecessary repeated prime-exponent work, but still lost to the current product tree.

Representative result:

```text
320000!
product tree  ~0.59 s
Prime-Swing   ~1.5 s
```

This is not a general statement that Prime-Swing is inferior; it is the result for mmCal's current backend and may be revisited after major backend changes.

---

# 6. Division

## Knuth normalized long division — retained

It remains a fast low-overhead base case for small/medium operands and for Burnikel–Ziegler recursion.

## Burnikel–Ziegler — adopted

Large balanced divisions use Burnikel–Ziegler. The GCC benchmark showed stable gains from roughly 32 limbs upward when the quotient is also large.

Representative results:

```text
1024-limb divisor/quotient  ~1.18 ms → ~0.12 ms
2048-limb divisor/quotient  ~5.0  ms → ~0.35 ms
```

This also accelerates `%`, Euclidean GCD, Rational normalization, integer roots, decimal conversion, and BigFloat division.

## Power-of-two division — adopted

Division by `2^k` uses shifts and low-bit extraction instead of general long division. At 4096 limbs this reduced millisecond-scale work to a few microseconds.

---

# 7. GCD

## Euclidean GCD — retained

It benefits automatically from faster division.

## Binary GCD — rejected

A Stein/binary GCD implementation was benchmarked and was roughly 4–30x slower in some current-backend cases. The next plausible candidate is Lehmer GCD rather than replacing Euclid with binary GCD.

---

# 8. Decimal conversion

Once arithmetic became faster, decimal output became the dominant cost for huge factorials.

## `10^9` chunks — adopted

Replacing digit-at-a-time `/10` conversion with 9-digit chunks reduced `40000!` conversion from roughly 6 seconds to roughly 0.5 seconds. Parsing uses the same chunking idea.

## Divide-and-conquer decimal conversion — adopted

Large powers `10^(9*2^k)` split the value roughly in half recursively, with `10^9` chunk conversion at the leaves.

Representative `40000!` (~166,714 decimal digits):

```text
digit-wise       ~6.0 s
10^9 chunks      ~0.53 s
divide-conquer   ~0.19 s
```

Formatter overhead was small compared with BigInt-to-decimal conversion itself.

---

# 9. `tryToUint64`

Clearly oversized integers are rejected by bit length before any decimal conversion. This removed a pathological path where a 100,000-bit value was rendered to decimal only to fail `uint64_t` conversion.

---

# 10. Extreme BigFloat exponent gaps

The old addition path always aligned exact exponents, so expressions such as `1 + 2^-5000000` created huge shifts.

v1.5.1 adds a fast path only when the requested precision, signs, and distance from a rounding boundary are sufficient to prove the result for all four rounding modes:

```text
NearestEven
TowardPositive
TowardNegative
TowardZero
```

Ambiguous cases fall back to the original exact-alignment path.

Representative 53-bit result: about 1.4 ms to sub-microsecond scale for `1 + 2^-5,000,000`.

---

# 11. Pi

The old certified Machin-series implementation remained mathematically clear but scaled poorly, reaching roughly 12 seconds for 10,000 digits.

Binary-splitting Chudnovsky is now the default:

```text
N[Pi,10000]  ~12.3 s → ~0.4 s
```

The old method is retained as reference code; the certified-enclosure contract is unchanged.

---

# 12. `exp` / `E`

Sequential RealInterval Taylor evaluation became dominated by interval-object and large-intermediate costs at high precision.

v1.5.1 uses binary splitting with certified range reduction. Small exact rationals use exact binary splitting; large numerators/denominators switch to fixed-precision interval binary splitting to avoid pathological exact intermediates.

Representative `N[E,5000]`: multi-second scale to roughly 0.1-second scale.

---

# 13. `log`

Sequential atanh-series evaluation was similarly expensive (`N[log[2],3000]` reached roughly 6–7 seconds).

v1.5.1 uses binary splitting. Large-bit mantissas are repeatedly reduced toward 1 with certified square roots before interval binary splitting.

Representative results:

```text
N[log[2],3000]                     ~6.7 s → ~0.1 s scale
N[log[123456789/987654321],5000] ~22 s   → ~3 s scale
```

Bit-burst or AGM-based logarithms remain future candidates for substantially higher precision.

---

# 14. Huge-radian trigonometric functions

The old path could send a huge raw radian rational almost directly to Taylor evaluation; `N[sin[10^6],20]` could exceed a five-second timeout.

v1.5.1 performs certified argument reduction using a certified Pi enclosure. It proves the quadrant integer for `x/(Pi/2)` before reducing to a small interval. It does not use machine `fmod` or a `double` Pi value.

The same method extends to interval-valued arguments such as expressions involving `sqrt[2]`.

---

# 15. FFT plan cache

Exact radix-2 transforms cache size-dependent bit-reversal and twiddle information across transforms within the evaluator/session. Input/output data are not cached, and the cache is not process-global.

---

# 15.4. Exact Cyclotomic FFT — quotient-field backend

## Previous generic `cis` / Expr direct DFT — retained as fallback

Before the quotient-field backend, exact non-power-of-two FFTs fell from `radix2Transform()` into `directTransform()`, constructing every twiddle as a generic `cis[-2 Pi k/n Rad]` expression. Even at lengths 5/7/10/12, forward transforms accumulated root-of-unity expressions and inverse transforms asked the generic Simplifier to rediscover the same cyclotomic identities, making exact round trips a major part of `test5_matrix.txt`.

The old dispatch and generic DFT implementation are not deleted. They remain the fallback for symbolic inputs and cases outside the new backend budget, and the old dispatch is preserved next to the replacement as a commented reference explaining the change.

## Canonical `NumberFieldContext` with ζ_n — rejected

The first prototype reused the general algebraic backend and represented the primitive root `zeta_n` as a canonical Complex Root inside `NumberFieldContext`. That representation is mathematically natural, but FFT arithmetic only needs repeated linear operations in one root-of-unity field; root isolation and canonical Root materialization dominated the work. The prototype measured roughly 2.5 s for a 5-point forward transform and 19.8 s for 7 points, with first construction of the ζ_7 context alone around 0.4 s. It was rejected.

The adopted design therefore performs the transform in a quotient representation without requiring an embedded root identity during the kernel. This is consistent with exact-DFT designs that operate over cyclotomic fields/quotients and evaluate the root of unity only at representation boundaries.

## `Q[t]/Phi_n(t)` quotient backend — selected

`CyclotomicFieldContext` stores only the conductor, exact cyclotomic polynomial, power-basis reduction, and coordinates of `t^k`. It has no root isolation or embedding object. Exact Rational inputs of non-power-of-two length at least 5 are transformed directly in Rational coordinates. Gaussian Rational inputs enlarge the conductor to `lcm(n,4)` when required, using `I=t^(3n/4)` exactly. Only the output boundary materializes the existing single generator vocabulary `cis[-2 Pi/n Rad]`.

The current budget is `phi(n)<=64`. Inputs above the budget, or symbolic expressions that cannot be proven to lie in the same quotient field, fall back to the legacy generic exact DFT. The approximate/certified FFT path is unchanged.

Representative GCC Release / LTO-off timings from `--exact-cyclotomic-fft 5`:

| length | first round trip | warm round trip |
|---:|---:|---:|
| 5 | ~0.72 ms | ~0.49 ms |
| 7 | ~1.79 ms | ~1.74 ms |
| 10 | ~1.50 ms | ~1.33 ms |
| 12 | ~1.44 ms | ~1.27 ms |
| 15 | ~9.11 ms | ~9.39 ms |
| 21 | ~41.1 ms | ~37.3 ms |

`tester.py --timings` reduced `test5_matrix.txt` from about 1.46 s in the earlier algebraic-field audit to about 0.38 s after the quotient-field backend.

## Persistent mixed-radix plan — deferred

The current non-power-of-two cyclotomic kernel is still a direct O(n^2) DFT over quotient coordinates. At 21 points the measured round trip is about 34 ms and the primary generic-expression explosion has already been removed. Mixed-radix Cooley–Tukey and exact Rader/Bluestein variants remain plausible, but should be added only after larger supported lengths demonstrate a crossover that justifies the extra implementation complexity.

# 15.5. Precision-aware `N` and certified FFT

## Old path — construct exact FFT first, approximate afterward

Previously `N` received already-evaluated arguments, so

```text
N[fft[data],16]
```

first built the complete exact Fourier expression and only then approximated its components. Each butterfly therefore paid the full generic `Expr` multiply/add/simplification cost even when the caller only wanted decimal output.

## Precision-aware evaluation — selected

`N` now holds its first argument, resolves the requested precision first, and keeps a precision context active while evaluating the child expression. FFT consumes that context and performs the transform directly on certified `ComplexInterval`/BigFloat endpoints. Ordinary exact `fft[...]` is unchanged, and no machine `double` backend is introduced.

Representative benchmark in the same GCC Release environment, measured under the older `N` semantics at 16 fractional digits:

```text
32 points   exact ~11.6 ms   certified ~2.7 ms
64 points   exact ~63.5 ms   certified ~6.2 ms
128 points  exact ~327.6 ms  certified ~13.1 ms
```

For non-power-of-two certified transforms, `--fft-threshold` compares direct DFT with forced Bluestein on the same 16-digit input. A GCC Release rerun measured 319 points at 1533 ms direct versus 1609 ms Bluestein, and 335 points at 1714 ms direct versus 1638 ms Bluestein, placing that environment's crossover between them. The supplied MSVC `--full` run instead measured 257 points at 1823/5002 ms and 509 points at 14007/7896 ms. The current policy conservatively favors the primary MSVC environment: direct below 384 points and Bluestein from 384 upward. This is an environment-dependent implementation threshold, not a mathematical constant.

The dedicated sweep covers 65 / 95 / 127 / 191 / 255 / 257 / 319 / 335 / 351 / 367 / 383 / 384 / 385 / 447 / 509 points. It forces each algorithm instead of comparing two calls routed through the existing policy. Results must have exactly equal finalized decimal real/imaginary component values; this ignores only representation differences such as `-12` versus `-12+0...I`.

A 2026-08-22 GCC Release re-audit measured 383 points at about 1990/1437 ms direct/Bluestein, 384 at 1001/1434 ms, and 385 at 2010/1414 ms. The winner is therefore not monotone around the boundary. Certified direct evaluation can change refinement and argument-reduction cost with the arithmetic structure of the transform length, so retuning the global threshold to one local crossover such as 335 would overfit this compiler run. The policy remains 384, while 383/384/385 are now permanent sweep points for future environment-specific audits.

The same audit separates mathematical convergence from practical bounded work for certified special functions. On this GCC Release build, `polylog[2,0.999]`, `2F1[...,0.98]`, and `ellipticF[...,0.98]` are mathematically inside their series convergence regions but become multi-second or worse with the current exact-majorant implementations. At the time, conservative work boundaries were `|z|<=49/50` for positive-order polylog, `<=9/10` for the 2F1/elliptic series, `|z|<=160` for 1F1, and 96 for the real Ei/Si/Ci series. Failures caused by a fixed backend range, term cap, or planner cap are classified as `CertifiedBackendUnsupported`, preventing meaningless guard-precision retries.

Why selected:

- preserves the exact-first public API while accelerating explicit approximation;
- keeps arbitrary precision and outward rounding without introducing machine `double`;
- makes precision propagation reusable by other expensive builtins;
- removes the approximate-path O(N^2) cliff for larger non-power-of-two sizes.

Exact symbolic FFT expression growth is a separate problem and is intentionally not changed here.

---

# 15.6. Array / Matrix storage and exact elimination

## Flat Array + exact Number backend — selected

v1.5.2 keeps matrices on the shared Array `shape + row-major flat storage` representation instead of introducing a separate nested matrix Value. `MatrixView` reads an Array without copying; algorithms that mutate their workspace use a flat `MatrixBuffer`. This describes the v1.5.2 release representation; v1.5.3 replaces the persistent physical storage with immutable paged backing plus strided views.

For exact Number matrices, pivot loops operate directly on `Number` values without constructing Expr nodes or invoking the Simplifier. Numeric `dot` likewise performs each cell accumulation directly in `Number`; symbolic products are built only when necessary.

General symbolic `det` / `inverse` use a shared expansion-work budget. Triangular matrices take an order-independent diagonal-product fast path and sufficiently sparse cases can still evaluate, while dense higher-order cases remain unevaluated before factorial expression growth occurs.

## Precision-aware certified Matrix — selected

Operations such as `N[det[A],p]` do not first finish an exact determinant. Like FFT, they receive the active `ApproximationContext` and dispatch directly to a `ComplexInterval` backend. Expression-to-interval conversion, decimalization, and guard-digit refinement are shared with the FFT approximation path.

`matrixRank` prefers exact elimination for exact inputs. The approximate backend uses no machine epsilon: an interval that still contains zero triggers guard-precision refinement, and rank deficiency is never inferred merely from a tolerance.

Example Release / LTO-off measurements from 2026-08-13:

| size | exact `dot` | exact `det` | exact `rref` | `N[det,16]` |
|---:|---:|---:|---:|---:|
| 8 | — | 0.283 ms | 0.434 ms | 1.410 ms |
| 12 | — | 1.290 ms | 2.233 ms | 10.144 ms |
| 16 | 0.328 ms | 3.735 ms | 5.799 ms | 22.799 ms |
| 32 | 1.702 ms | — | — | — |
| 64 | 16.553 ms | — | — | — |

The initial implementation intentionally retained ordinary Gaussian/Gauss-Jordan exact elimination; integer/Rational Bareiss elimination was introduced later as a separate fraction-free optimization.

# 15.7. Bareiss / fraction-free exact Matrix — selected

The fraction-free backend lifts exact real matrices to integers by clearing denominators independently per row, then runs a shared Bareiss kernel on `IntegerMatrixBuffer`. Integer inputs require no lift; Rational inputs are scaled only for the elimination workspace.

- `det` uses Bareiss forward elimination and restores the product of row denominator scales once at the end.
- `rref` remains fraction-free through the forward phase. Full-column-rank results are materialized directly as identity columns, while only rank-deficient cases enter canonical Rational backward normalization.
- `matrixRank` uses the number of Bareiss pivots without materializing a full RREF.
- `inverse` writes `B=D A`, performs fraction-free forward elimination on `[B|D]`, then uses the final pivot as a common denominator for BigInt-only back substitution instead of generic Rational Gauss-Jordan.
- `solveLinear` shares the same BigInt back-substitution path for unique full-column-rank systems and constructs Rational values only at the end.
- Exact complex matrices retain the previous `Number` Gaussian/Gauss-Jordan fallback until a dedicated exact complex integer-domain representation is justified.

Pivot selection prefers the nonzero candidate with the smallest bit length to limit intermediate BigInt growth. Every Bareiss division is checked with `BigInt::divmod`; a nonzero remainder is treated as an invariant failure rather than silently truncating.

Release / LTO-off measurements on 2026-08-13 using the same benchmark matrices:

| size | Gaussian `det` | Bareiss `det` | speedup | Gauss-Jordan `rref` | Bareiss `rref` | speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 0.256 ms | 0.047 ms | 5.4x | 0.447 ms | 0.058 ms | 7.7x |
| 12 | 1.241 ms | 0.109 ms | 11.4x | 2.058 ms | 0.146 ms | 14.1x |
| 16 | 3.570 ms | 0.385 ms | 9.3x | 5.798 ms | 0.375 ms | 15.5x |

`N[det[...],p]`, `N[inverse[...],p]`, and related operations do not build these exact Bareiss results first. They continue to dispatch directly to the precision-aware certified Matrix backend shared with the FFT approximation infrastructure.

A 2026-08-25 performance-cliff audit showed that exact inverse cost was concentrated after Bareiss forward elimination in the generic Rational RREF phase. Full-rank square inverse and unique solve now use the final pivot `D` as a common denominator and perform exact BigInt back substitution

```text
n_i = (b_i D - sum_{j>i} U_ij n_j) / U_ii.
```

Full-column-rank `rref` similarly stops after the forward phase and materializes the known identity-column result directly. Rank-deficient `rref/nullSpace` retains the canonical Rational backward phase.

On the same public-path benchmark, 16-bit 32x32 measured about 103 ms for `inverse`, 5.2 ms for `rref`, 5.2 ms for `matrixRank`, and 5.8 ms for a nullity-one 32x33 `nullSpace`. At 96 bits, the 32x32 figures were about 0.70 s / 34.7 ms / 34.1 ms / 37.6 ms; at 256 bits, about 3.46 s / 147 ms / 145 ms / 144 ms. Rank-deficient 48x49 / 256-bit `rref/rank/nullSpace` measured about 1.25 / 1.13 / 1.10 s. No abrupt algorithmic cliff was observed in these structural paths; growth tracks coefficient height. High-bit inverse is increasingly dominated by canonicalizing the final 1024 huge Rational elements, so shared-denominator persistent storage is a future representation problem rather than a reason to select the slower modular inverse backend.

# 15.7.1. Modular / CRT exact Matrix — selected

Post-v1.5.3 adds a modular backend for larger dense integer/Rational matrices so that Bareiss intermediate BigInts do not dominate bit complexity. The finite-field kernel uses 31-bit primes, keeping products safely inside `uint64_t` and avoiding an MSVC-specific `__int128` dependency. Rational inputs continue to clear denominators per row before entering the integer workspace.

- `det`: performs Gaussian elimination over several prime fields, incrementally combines images with CRT, and stops once an integer-only Hadamard bound guarantees a unique centered reconstruction.
- `solveLinear`: combines finite-field solutions by CRT, applies rational reconstruction, and accepts a candidate only after exact verification of `A X = B` in the original integer system. Bad primes are skipped; unsuccessful reconstruction falls back to Bareiss.
- `inverse`: uses exact `det(A)` as the common denominator, reconstructs integer `adj(A)` from modular inverse images, and verifies `A adj(A)=det(A)I` exactly. The backend exists, but automatic dispatch does not use it because Bareiss remains faster throughout the measured range.
- `rref`, `matrixRank`, and `nullSpace` remain on Bareiss pending a separate modular-certificate design.

The automatic dispatcher requires coefficient density of at least 25% and uses the following conservative GCC Release crossover policy. `height` is the maximum coefficient bit length in the integer workspace.

| operation | modular condition |
| --- | --- |
| `det` | order>=48, or order>=32 & height>=64, order>=24 & height>=192 |
| `solveLinear` | variables>=24, or variables>=12 & height>=96, variables>=8 & height>=256, variables>=6 & height>=512 |
| `inverse` | never selected automatically |

Representative GCC Release / LTO-off timings from 2026-08-20, averaged over three iterations, in milliseconds:

| workload | Bareiss | modular |
| --- | ---: | ---: |
| `det` 32x32, 96-bit | 37.691 | 26.270 |
| `det` 24x24, 256-bit | 41.221 | 29.889 |
| `det` 32x32, 256-bit | 174.095 | 58.399 |
| `det` 20x20, 512-bit | 55.836 | 56.766 |
| `det` 24x24, 512-bit | 136.994 | 70.267 |
| `det` 32x32, 512-bit | 497.769 | 143.100 |
| `solveLinear` 12x12, 96-bit | 0.453 | 0.234 |
| `solveLinear` 8x8, 256-bit | 0.297 | 0.179 |
| `solveLinear` 24x24, 256-bit | 47.954 | 0.686 |
| `solveLinear` 6x6, 512-bit | 0.285 | 0.113 |
| `solveLinear` 24x24, 512-bit | 156.020 | 0.703 |
| `solveLinear` 32x32, 512-bit | 528.836 | 1.586 |

At 20x20 with 512-bit coefficients Bareiss remained marginally faster in this sweep, so automatic dispatch deliberately stays on Bareiss there. The 24x24 512-bit case clearly favors modular and is already covered by the 24+/192-bit rule.

A 2026-08-25 remeasurement still found no inverse crossover: at 32x32 the 16-bit case measured `15.2 / 84.5`, 96-bit `74.2 / 657.6`, and 256-bit `357 / 3339` ms (Bareiss forward core / modular inverse). Automatic modular inverse therefore remains disabled. The public exact `inverse` path has an additional output-materialization cost; 32x32 / 256-bit is about 3.46 s because the final 1024 canonical Rational values are themselves large.

These thresholds are performance policy, not mathematical semantics. Compiler, BigInt, prime-kernel, or CPU changes should be remeasured with `mmCal.Benchmarks --exact-linear-algebra 1`. Certified `N[det[...],p]` / `N[solveLinear[...],p]` paths still dispatch directly to precision-aware interval kernels without first performing exact reconstruction.

# 15.8. LU / fraction-free exact QR / certified Householder QR — selected

Decomposition code is centralized in `linear_algebra/decomposition.*`. `luDecomposition[A]` keeps row-pivoted `P A = L U`. QR deliberately uses different exact and approximate algorithms. `N[qrDecomposition[A],p]` continues to dispatch directly to certified Householder without first constructing exact factors, while exact real matrices were moved from Expr-level Householder expansion to fraction-free orthogonalization on 2026-08-25.

## Exact QR: delay normalization until materialization — selected

The old exact Householder path repeatedly performed `norm -> sqrt -> reflector -> Expr arithmetic -> simplify`. Once the first radical entered later column norms, expression growth became extreme: roughly 1.2 ms at 2x2, 59 ms at 3x3, but about 18 s and roughly 677 KB of formatted output at 4x4. The former `maximumExactQrOrder = 3` was a policy guard against this blow-up.

The replacement first lifts Rational columns to primitive integer vectors and creates neither square roots nor Rational divisions during orthogonalization. With integer orthogonal vector `p_i` and `d_i=p_i^T p_i`, projection removal uses

```text
p <- d_i v - (p_i^T v) p_i
```

with GCD content reduction. Only final materialization builds

```text
Q[:,i] = p_i / sqrt(d_i)
R[i,j] = (p_i^T a_j) / sqrt(d_i).
```

Radical decomposition is performed once per column and the same radical node is shared by all Q/R elements. This removes the old `maximumExactQrOrder` hard cap. The upper-triangular/trapezoidal `{I,A}` fast path remains.

## Gram + symmetric Bareiss (fraction-free LDL^T-equivalent) — selected

For full-rank input, primitive columns `C` form `G=C^T C`; symmetric Bareiss elimination yields the principal-determinant sequence and lower coefficients from which an orthogonal integer basis is recovered without square roots. This measured roughly 3–5x faster than direct fraction-free Gram-Schmidt over representative orders 8–64. If a leading principal minor vanishes or the matrix is rank deficient, direct fraction-free orthogonalization is used instead and standard-basis directions are completed by the same kernel.

A conventional Rational LDL^T implementation was measured and rejected: Rational normalization made it roughly 6–35x slower than the direct fraction-free path. Forming `A^T A` does not introduce numerical-conditioning problems in exact arithmetic; the relevant costs here are BigInt growth and measured runtime.

A second experiment skipped square-factor extraction for large norms and emitted the exact but less canonical identity `1/sqrt(n)=sqrt(n)/n`. It gave no material improvement at workloads such as 32x32 / 256-bit and sometimes regressed slightly, while weakening canonical radical form, so it is rejected. The current implementation keeps one canonical radical decomposition per column.

Representative GCC Release / LTO-off public-path timings from `--exact-linear-algebra 1` on 2026-08-25, including final Rational/radical Expr materialization:

| coefficient | size | exact QR |
|---:|---:|---:|
| 16 bit | 8x8 | 2.22 ms |
| 16 bit | 16x16 | 14.3 ms |
| 16 bit | 24x24 | 54.0 ms |
| 16 bit | 32x32 | 147 ms |
| 96 bit | 16x16 | 96.4 ms |
| 96 bit | 24x24 | 418 ms |
| 96 bit | 32x32 | 1.30 s |
| 256 bit | 16x16 | 415 ms |
| 256 bit | 24x24 | 2.00 s |
| 256 bit | 32x32 | 6.85 s |

The old 4x4 expression cliff is gone at low and medium coefficient heights. High-order/high-bit workloads are now dominated by the final huge Rational/radical output rather than by the orthogonalization kernel. The implementation therefore does not restore an order hard cap; EvaluationBudget and natural output cost remain the boundary. A future improvement would require a shared-denominator or delayed-radical internal representation rather than silently switching exact work to machine approximation.

## Certified Householder — retained

The approximate Householder kernel also has a column-block experiment that processes multiple columns in one row-major scan. Release measurements of block sizes 1/8/16/32 over orders 8–24 stayed within a few percent and no block size won consistently because BigFloat/interval arithmetic dominates. The default therefore remains effectively block=1, while the block kernel and benchmark are retained.

Representative 2026-08-13 Release / LTO-off timings:

| size | exact `LU` | `N[LU,16]` | `N[QR,16]` |
|---:|---:|---:|---:|
| 8 | 0.290 ms | 1.811 ms | 7.062 ms |
| 12 | 1.315 ms | 5.355 ms | 21.218 ms |
| 16 | 3.548 ms | 10.562 ms | 46.718 ms |

# 15.9. Reduced SVD — selected

The numerical SVD backend deliberately avoids forming `A^H A`, which would square the condition number. It performs Householder bidiagonalization followed by one-sided Jacobi column orthogonalization. Real and complex inputs share the same precision-aware policy; candidate factors are accepted only after interval verification of reconstruction and U/V orthogonality. Exact SVD is limited to natural closed cases.

Release / LTO-off repeated measurements on 2026-08-13 for `N[svd[A],16]`:

| size | time |
|---:|---:|
| 4x4 | about 3.8 ms |
| 8x8 | about 17.1 ms |
| 12x12 | about 46.6 ms |
| 16x16 | about 82.3 ms |

The observed range is smooth and consistent with the expected cubic regime at these sizes. The dominant cost is the Jacobi iteration and BigFloat arithmetic. No extra SVD blocking policy is enabled yet: adding a cache block around the smaller bidiagonalization fraction would not justify itself without a measured win. The existing QR block kernel remains benchmarked independently.

# 15.10. Eigen / complex Schur — selected

The general eigenproblem is reduced directly on Complex BigFloat values: Hessenberg reduction → implicit shifted QR → complex Schur form. Eigenvectors, when requested, are recovered by back substitution on the triangular Schur matrix and mapped through the accumulated Schur vectors. The exact path prefers triangular/diagonal matrices and distinct-root exact Number 2x2 cases.

A fixed-seed 3x3 case exposed that using essentially the same precision for QR stopping and the final relation certificate left too little accumulation margin. The internal QR stopping target is therefore ten decimal digits stricter than the requested display precision. The original `ComplexInterval` input is used to audit `A Q-Q T`, `A v-lambda v`, and Schur-vector unitarity. For non-normal matrices, this is a relation certificate rather than an unconditional claim of unique componentwise eigenvalue/eigenvector enclosures.

Release / LTO-off measurements on 2026-08-13 using random decimal matrices in [-1,1] with ten fractional digits:

| size | `N[eigenvalues,16]` | `N[eigensystem,16]` |
|---:|---:|---:|
| 4 | 13.4 ms | 14.0 ms |
| 8 | 68.1 ms | 77.0 ms |
| 16 | 364.8 ms | 466.3 ms |
| 32 | 2826 ms | 3658 ms |
| 64 | 19436 ms | >35 s measurement cap |

# 15.11. Large dense Matrix audit

`mmCal.Benchmarks` now provides `--matrix-large <op> <size> [digits]`. Algorithm timings exclude parser cost by directly constructing exact Rational elements with the same [-1,1] range and ten-decimal granularity as the supplied random-matrix generator; the RNG stream is not intended to be Python-identical. CLI parse/lowering cost is measured separately with the Python text format.

Representative order-32/order-64 results:

| op | 32x32 | 64x64 |
|---|---:|---:|
| `N[dot,16]` | 167 ms | 1.21 s |
| `N[det,16]` | 294 ms | 2.96 s |
| `N[inverse,16]` | 1.22 s | 8.68 s |
| `N[matrixRank,16]` | 400 ms | 4.16 s |
| `N[solveLinear,16]` | 549 ms | 5.42 s |
| `N[nullSpace,16]` | 401 ms | 4.01 s |
| `N[LU,16]` | 147 ms | 1.41 s |
| `N[QR,16]` | 1.47 s | 12.88 s |
| `N[SVD,16]` | 1.40 s | 11.60 s |
| `N[eigenvalues,16]` | 2.83 s | 19.44 s |

At 1024x1024, representation cost becomes a primary limit before the cubic algorithms themselves. Direct construction of 1,048,576 ten-decimal Rational Expr elements reached about 0.69 GB maximum RSS; the timed transpose itself took about 107 ms and trace about 60 ms. Parsing roughly 14.16 MB of generator-style text and evaluating only `dimensions[...]` took about 10.9 s wall time and about 1.99 GB maximum RSS. A 1024-order `N[dot,16]` run did not complete within a 10 s cap and reached about 0.96 GB RSS; `N[LU,16]` likewise exceeded 10 s and reached about 1.59 GB. Further 1024 QR/SVD/Eigen runs were stopped to avoid unnecessary memory pressure.

After the v1.5.3 `Expr::Node` typed-node refactor, an apples-to-apples x86-64 GCC Release/LTO-off rerun of `--matrix-large transpose 1024 16` measured `693312 KiB` (~677.1 MiB) maximum RSS for the legacy variant source and `299668 KiB` (~292.6 MiB) for the typed-node source: about 384.4 MiB / **56.8% less RSS**. The one-shot transpose timing changed from 181.7 ms to 157.4 ms, but timing noise is not the adoption criterion; the memory reduction plus unchanged regression semantics are.

The persistent representation packs `ArrayExpr` values into fixed immutable pages and separates shape/offset/strides from the backing. A first naive single-`vector<Rational>` design reduced storage overhead but made transpose deep-copy one million Rational/BigInt values; direct-packed order-1024 transpose took roughly 650–675 ms, so that design was rejected. The adopted design shares immutable 1024-element pages and makes transpose a stride-only view operation.

`ArrayBuilder` promotes only the current page. Completed pages are immutable, so a symbolic value near the end of a huge numeric Array does not force the preceding data through a Generic-Expr rebuild. A dedicated 1,048,576-element test with only the final value changed to `x` was essentially identical to the all-integer case at about 0.30 s / 69.8 MiB, with only the final page becoming Generic. Rectangular brace lowering also streams numeric leaves directly into one builder instead of first allocating an Expr node for every scalar.

With the benchmark fixture likewise constructing exact Rationals directly through `ArrayBuilder`, `--matrix-large transpose 1024 16` now measured `136576 KiB` (~133.4 MiB) maximum RSS and about 0.059 ms for the transpose view itself. `trace 1024` measured about 133.5 MiB and 10.95 ms for trace. Because the fixture construction path is intentionally part of the new representation work, the 292.6→133.4 MiB change is an end-to-end storage+builder improvement rather than the same kind of single-change A/B used for typed nodes.

A 13.63 MB 1024x1024 ten-decimal literal fed through the CLI and evaluated only as `dimensions[...]` measured about 4.55 s wall time and `441564 KiB` (~431 MiB) maximum RSS. The older 14.16 MB / 10.9 s / 1.99 GB measurement did not use byte-identical input, so this is not presented as a strict A/B comparison, but it confirms that avoiding per-scalar Expr allocation in lowering materially reduces the post-parse representation cost.

Even a pure cubic extrapolation from order 64 suggests roughly 1.6 h for `N[LU]`, 1.4 h for `N[dot]`, 3.4 h for `N[det]`, 6.2 h for `N[solveLinear]`, 9.9 h for `N[inverse]`, 13 h for `N[SVD]`, 15 h for `N[QR]`, and 22 h for `N[eigenvalues]`. Extrapolating the observed 32→64 exponent instead gives a broad roughly 1–21 h range depending on the operation. These are projections, not 1024 completion measurements, and cache/allocation/guard-precision/iteration effects can make them worse.

A dense order of 1024 is not intrinsically huge in a machine-double + BLAS setting, but it remains a stress regime for mmCal's certified arbitrary-precision dense algorithms. Persistent exact-Array representation cost is now substantially lower after typed nodes, paged packed backing, and direct builder lowering. Approximate SVD/Eigen paths already use dedicated contiguous working buffers, so blocking and selective threading should now be remeasured against the new storage balance instead of introducing another persistent-matrix type first.


# 15.12. Persistent algebraic fields and bounded caches

## Compositum / embedding reuse — adopted

The compositum cache retains a proven primitive-element compositum and both operand power-basis embeddings for a Root pair in a bounded cache. Real `AlgebraicElement` materialization also evaluates the coordinate polynomial over the chosen generator interval and uses an exact Sturm root count to identify the canonical `root[minpoly,k]` directly. This avoids duplicate all-root isolation and impossible Rational / `Q+iQ` degeneration probes.

For the representative expression

```text
(root[{-2,0,1},2]+root[{-3,0,0,1},1])
*(root[{-2,0,1},2]-root[{-3,0,0,1},1])
```

The earlier implementation took about 1.23 s per evaluation. Persistent compositum caching reduced this to roughly 0.14–0.15 s on the first evaluation and about 0.02 s warm in the same GCC Release / LTO-off environment.

## Reciprocal reuse — adopted

The reciprocal cache memoizes exact reciprocals computed by extended Euclid in `Q[t]/(m)` in a thread-safe per-`NumberFieldContext` LRU capped at 16 reciprocal pairs. One entry is bidirectional because `inverse(inverse(x)) = x`. Rational constant coordinates invert directly through the canonical embedding of `Q`.

Representative measurements:

| degree | operation | before cache | cached reciprocal warm |
|---:|---|---:|---:|
| 6 | reciprocal | ~178 us | ~0.24 us |
| 6 | divide | ~243 us | ~84 us |
| 12 | reciprocal | ~525 us | ~0.44 us |
| 12 | divide | ~773 us | ~211–233 us |

### Persistent multiplication-matrix cache — rejected

For a degree-12 field, a representative ordinary multiplication took about 112 us, matrix-vector multiplication about 101 us, while constructing the left-multiplication matrix cost about 228 us. The roughly 10% per-multiply saving is not enough to amortize construction unless the same multiplier is reused many times. The memory and complexity cost of attaching such a cache to general `AlgebraicElement` values is therefore not justified at present.

## Incremental Krylov minimal polynomial — adopted

The old `AlgebraicElement::minimalPolynomial()` rebuilt a Rational matrix and reran Gauss-Jordan from scratch for every candidate degree in `1,a,...,a^k`, repeatedly discarding the same independence information.

Incremental minimal-polynomial derivation appends the Krylov sequence `1,a,a^2,...` one column at a time and reduces only the new column against a persistent exact row-echelon state. The first dependence

```text
c0 + c1 a + ... + a^k = 0
```

is the minimal polynomial because all earlier powers are linearly independent. The primitive-element tensor path uses the same incremental basis both for the minimal polynomial of `theta` and for conversion of `alpha` / `beta` into the established power basis.

Representative GCC Release / LTO-off results:

| field degree | repeated Gauss-Jordan | incremental Krylov first |
|---:|---:|---:|
| 6 | ~651 us | ~475 us |
| 12 | ~18.9 ms | ~7.7–7.9 ms |

Derived minimal polynomials are also retained by exact power-basis coordinate key in a thread-safe per-field LRU capped at 16 entries. A warm degree-12 hit is about 0.59 us. Cache misses and eviction only cause recomputation and cannot change the mathematical result.

### Multiplication-matrix minpoly / modular reconstruction — deferred

Libraries such as FLINT provide exact Rational-matrix minimal-polynomial backends, and representing a finite algebra element by its multiplication matrix is standard. mmCal already has bounded-degree power-basis coordinates, however, and incremental Krylov gives a substantial improvement with much less machinery in the current degree range. A permanent multiplication-matrix-to-matrix-minpoly path is therefore not added. Modular images with rational reconstruction and fraction-free matrix minpoly remain candidates if higher degrees or coefficient heights make first derivation dominant again.

# 15.13. Black-box workload audit

`test_set/tester.py` now accepts `--timings [N]`, recording mmCal process wall time for each test file and printing the slowest files. Test-file I/O and parsing remain outside the overall timed region as before.

Representative incremental-Krylov / GCC Release / LTO-off results for all 1652 black-box cases on 2026-08-16:

| test file | tests | wall time |
|---|---:|---:|
| `test16_exact_calculus_solver.txt` | 85 | ~2927 ms |
| `test5_matrix.txt` | 105 | ~1457 ms |
| `test9_special_func.txt` | 75 | ~297 ms |
| `test8_calculus.txt` | 46 | ~173 ms |
| `test22_number_field_interning.txt` | 1 | ~152 ms |

`test16` is mainly a stress set for integration, high-degree Solve, and algebraic Root construction. Despite its name, much of the time in `test5_matrix` is in exact FFT/DFT non-power-of-two round trips rather than small matrix operations; exact `ifft[fft[...]]` around 7, 12, and related sizes remains a future performance target. In `test9`, `N[ibeta[1/3,2/3,1/4],20]` and `N[gamma[1/3],20]` stand out comparatively.

These timings never affect PASS/FAIL semantics; they are profiling signals used only to select optimization targets.

# 15.14. Certified `gamma` / `ibeta` — interval and cache optimizations

The `tester.py --timings` audit identified `N[ibeta[1/3,2/3,1/4],20]` and `N[gamma[1/3],20]` as clear hotspots inside `test9_special_func.txt`, so both certified backends were measured directly.

## `ibeta` point/shared normalization — adopted

The old interval wrapper evaluated both endpoints even when `lower == upper`, and each endpoint rebuilt `Beta(a,b)`. Since the 2F1 series itself was only about 1.9 ms at the representative 20-digit workload while Beta/Gamma normalization dominated, The optimized path now:

- evaluates an exact point once;
- computes `Beta(a,b)` once and shares it across interval endpoints;
- reuses the same normalization through the complement identity because `B(a,b)=B(b,a)`;
- uses only +40 guard bits for non-complement point evaluation and +80 when complement/interval evaluation requires it;
- returns exact `x=0,1` before constructing the normalization.

Representative warm/direct timings moved from roughly 116→37–43 ms at 80 bits, 255→98–117 ms at 160 bits, 835→387–392 ms at 320 bits, and about 10.3→5.3 s for the first 640-bit call. Repeated 640-bit calls fall to about 1.0 s once the Gamma plan cache is warm.

## Gamma lazy Bernoulli / Horner / plan reuse — adopted

The old Gamma backend eagerly generated `B0...B128` on first use. The revised Gamma backend preserves the Akiyama–Tanigawa state and extends it only to the requested even Bernoulli order under a mutex. It also limits low/mid-precision plan search to `min(64,max(16,ceil(bits/5)))`, evaluates the Stirling polynomial in Horner form, builds exact-point recurrence products as balanced Rational products, adds exact `Gamma(1)=Gamma(2)=1` / `logGamma(1)=logGamma(2)=0` fast paths, and keeps a thread-local bounded cache of exact Stirling plans.

Separate-process cold timings for `gamma[1/3]` improved from about 79→9.5 ms at 80 bits and 89→43 ms at 160 bits. The 320/640-bit first-call cost is approximately unchanged; the remaining high-precision cost is in the actual Stirling/recurrence work. Repeated 640-bit `lgamma[1/3]` falls from roughly 1.36 s to about 0.34 s when the plan cache is warm.

### Eager Bernoulli generation through `B256` — measured and rejected

Simply extending the old Akiyama–Tanigawa initialization through `B256` was also evaluated as a prerequisite for longer Stirling sums. Exact Rational generation alone rose from roughly 68 ms through `B128` to roughly 560 ms through `B256`, imposing that fixed first-use cost even on low-precision calls. The stateful lazy cache is therefore retained instead.

### `maximumK > 64` to reduce shift — measured and rejected

Allowing `K=96` does reduce the large recurrence shift at roughly 640 bits, but the larger exact Bernoulli coefficients, Rational/interval conversion and longer Stirling sum outweighed that gain: the representative workload regressed from about 1.37 s to about 2.7 s. Keep `K<=64` until a rectangular/binary-splitting Stirling kernel and a cheaper high-order Bernoulli backend justify reevaluation.

### fixed-k exact binary-search planning — measured and rejected

A fixed `k=64` binary search over the exact remainder bound was also tested. Repeated construction of huge Rational powers `x^(2k-1)` made the planner itself slower than the existing scan, so it is not retained.

The replaced implementations remain commented next to the new code in `certified_special_functions.cpp`, together with the reason for replacement, specifically for this algorithm-comparison cycle.

Note: the stateful lazy Bernoulli generator described in this section was subsequently replaced by a static exact `B_2...B_128` table by a static exact `B_2...B_128` table; the generator remains only as commented comparison code.

# 15.15. Exact-Rational `Gamma` / high-precision Stirling planning

Reprofiling the 640-bit and higher `gamma[1/3]` / `ibeta[1/3,2/3,1/4]` paths after the interval/cache optimization showed that major costs remained at representation boundaries and in parameter planning, before any need for a fundamentally different Gamma formula. The design was cross-checked against Fredrik Johansson, *Arbitrary-precision computation of the gamma function* (arXiv:2109.08392), especially its treatment of rational rising factorials and Stirling parameter selection.

## Preserve exact Rational identity through the certified backend — adopted

The evaluator knew `1/3` exactly, but the old special-function path first converted it with `RealInterval::fromRational`. A non-dyadic Rational is not a point interval, so the exact rising-factorial path added in the preceding optimization was effectively bypassed for the representative `gamma[1/3]` workload. The high-precision path now:

- detects exact Rational arguments before interval conversion for `Gamma`, `LogGamma`, `Beta`, and `BetaLog`;
- preserves the original `p/q` through positive-Rational Gamma argument shifting;
- forms `(p/q)_n` as a balanced binary product `prod(p+qk)/q^n` and canonicalizes one final Rational;
- passes exact `a`, `b`, and `a+b` into the three LogGamma evaluations used by Beta;
- preserves exact Rational identity under negative-argument reflection and evaluates `sin(Pi x)` as exact turns `sinTurns(x/2)`.

## Static exact table for `B_2...B_128` — adopted

The stateful lazy Akiyama–Tanigawa generator from the preceding implementation avoided low-precision eager initialization, but the first call that reached the highest Bernoulli orders still paid a large exact-Rational state-update cost. The current Stirling kernel only needs the fixed constants `B_2...B_128`, so they are stored as exact decimal numerator/denominator literals and parsed into `BigInt/Rational` only on first reference. The replaced lazy generator remains commented beside the new implementation with its replacement rationale. This does **not** revive the rejected idea of runtime eager generation through `B256`; higher dynamic Bernoulli generation remains deferred.

## Certified BigInt high-precision planner — adopted

The old planner advanced the shift in increments of eight and tested `k=1...64` with exact Rational arithmetic at every candidate. At 1000-bit precision this planning work itself became a large cold-start cost. For positive Rational `x=(p+qs)/q`, coefficient `c=A/B`, and `d=2k-1`, the test

```text
|c| / x^d <= 2^-P
```

is exactly equivalent to

```text
|A| q^d 2^P <= B (p+q s)^d.
```

Above 768 bits, the current planner tests this integer inequality for `k=64` and finds a sufficient shift with doubling plus binary search. The final remainder bound is rebuilt as an exact Rational, so no floating heuristic weakens the certified contract.

This is distinct from the fixed-k Rational-power binary search rejected in the earlier planner: that version repeatedly constructed normalized Rational powers `x^(2k-1)`. The adopted planner removes those GCD/normalization costs and compares BigInt cross products directly.

## Representative measurements

Representative GCC Release / LTO-off timings in the same environment:

| workload | earlier analysis baseline | current planner |
|---|---:|---:|
| `gamma[1/3]`, 640 bit | about 0.69 s | about 65 ms first / 38 ms warm average |
| `gamma[1/3]`, 1280 bit first | about 1.5 s | about 0.11 s |
| `ibeta[1/3,2/3,1/4]`, 640 bit | about 5.3 s | about 0.21 s |
| `ibeta[1/3,2/3,1/4]`, 1280 bit | about 4.5 s | about 0.49 s |
| `gamma[-1/3]`, 1280 bit | generic interval reflection | about 0.15 s |

A normal `--special-functions 3` run gives roughly 4.5/7.2/13.7/37.6 ms for `gamma[1/3]` at 80/160/320/640 bits and about 15.8/32.2/44.3/210 ms for the corresponding `ibeta` workload. The permanent benchmark now also includes 1280 bits.

## Improved Stirling main sum / Algorithm 6 — next candidate

Johansson's Theorem 3.5 / Algorithm 6 splits the Stirling main sum into low-index Bernoulli terms and a re-expanded high-index hypergeometric tail, reducing the number of Bernoulli values needed at high precision. The FLINT/Arb Gamma backend likewise documents an improved Stirling sum using rectangular splitting for low-index terms and high-index re-expansion. The current planner removes the dominant representation/planner overhead through the 1280-bit range; future work above this range should therefore implement this improved main-sum direction instead of simply increasing `K` or returning to runtime `B256` generation.

# 16. `mmCal.Benchmarks`

v1.5.1 adds a separate console project to the Visual Studio solution:

```text
mmCal
mmCal.Core
mmCal.Tests
mmCal.Benchmarks
```

It keeps benchmark workloads out of normal regression-test timing while reusing `mmCal.Core`.

It includes:

- fixed-seed randomized BigInt division invariants;
- decimal round trips;
- certified `exp(x)exp(-x)` invariants;
- certified `log(x)+log(1/x)` invariants;
- multiplication/square/division threshold benchmarks;
- factorial and decimal I/O benchmarks;
- high-precision `Pi/exp/log` benchmarks;
- certified `gamma/ibeta` precision-scaling benchmark (`--special-functions`);
- fixed-seed certified Matrix invariants, including Bareiss / LU / QR / solve / nullSpace, real/complex SVD, and Eigen;
- exact/certified FFT benchmarks and direct/Bluestein crossover measurements;
- fixed-seed certified FFT round-trip invariants.

Modes:

```text
mmCal.Benchmarks
mmCal.Benchmarks --full
mmCal.Benchmarks --random-only
mmCal.Benchmarks --benchmark-only
mmCal.Benchmarks --matrix-large nsvd 64 16
mmCal.Benchmarks --special-functions 1
```

Correctness checks should run before accepting any new threshold solely because it benchmarks faster.

---

# 17. Explicitly not selected for v1.5.1

| Candidate | Decision | Reason |
|---|---|---|
| Prime-Swing factorial | rejected | slower than the current product tree |
| binary GCD | rejected | 4–30x slower in some current cases |
| Karatsuba vector pool | rejected | ~5–10% regression |
| Karatsuba depth scratch | rejected | management cost outweighed savings |
| dedicated Toom-3 square | rejected | slower than Karatsuba square |
| low Toom-3 threshold | rejected | overhead wins around 512–1024 limbs |
| machine `fmod` for huge trig reduction | rejected | loses certified semantics |
| silently converting ordinary evaluation to Machine/double | policy rejection | changes exact-first semantics |
| persistent algebraic multiplication-matrix cache | rejected | degree-12 build cost ~228 us while multiply improved only ~112→101 us; poor amortization |
| multiplication-matrix minpoly / modular reconstruction | deferred | incremental Krylov is sufficient in the current bounded degree range and needs much less machinery |
| runtime eager Bernoulli generation through `B256` | rejected | exact Rational generation alone costs about 560 ms. the current implementation uses a fixed exact `B_2...B_128` table, but does not dynamically generate a larger range |
| Gamma Stirling `maximumK>64` | rejected | reduced shift but high-order Bernoulli/Rational and longer Stirling work regressed the 640-bit case from ~1.37 s to ~2.7 s |
| Gamma fixed-k **Rational-power** binary search | rejected | normalized Rational `x^(2k-1)` probes are expensive; the current implementation adopts a different BigInt cross-product formulation of the same certified test |

---

# 18. Next candidates

The `Expr::Node` typed-node refactor is adopted in v1.5.3. It was kept as a standalone public-API-preserving change, passed the regression/fuzzer checks, and reduced order-1024 Matrix RSS by about 56.8% in the same-environment comparison above.

1. Remeasure blocked LU / QR against the new paged packed Array cost balance.
2. Evaluate threading thresholds only for pure numeric working kernels.
3. Audit parser-AST temporary allocation for huge braces if it remains material.
4. Keep BigUInt / BigInt SBO as an isolated experiment and adopt it only if the complexity can be confined to a storage abstraction.
5. Measure Toom-4 / higher-Toom crossovers and consider FFT/NTT multiplication for still larger integers.
6. Lehmer GCD.
7. bit-burst / AGM logarithm backends.
8. A Cyclotomic exact FFT backend. The `tester.py --timings` audit shows non-power-of-two exact FFT/iFFT round trips as a clear interactive black-box hotspot, so this priority should be revisited.

The naive flat packed-Array design is rejected; immutable paged backing plus stride views is adopted. Approximate Matrix algorithms keep their existing dedicated contiguous working buffers rather than forcing persistent Array storage and algorithm temporaries into one type. BigUInt SBO is explicitly deferred for now.

Future changes should continue to record both adoption and rejection rationale in this document.

# 19. Cross-function special-function performance-cliff audit

A new `mmCal.Benchmarks --performance-cliffs [iterations]` mode sweeps representative argument and precision ranges through the public `KernelSession`, recording wall time and `EvaluationUsage` together. Numeric-to-held/error transitions are reported as `BOUNDARY`, while only sharp growth between successful neighboring points is reported as a `CLIFF` candidate. This prevents intentional bounded-work exits from being mistaken for slow evaluation.

A 2026-08-25 GCC Release / LTO-off audit (two runs, warm values) found the main current hotspots to be: `1F1[1/2,5/4,160]` at about 1.08 s, `polylog[2,49/50]` at about 139 ms, elliptic F/E/Pi near the 9/10 parameter boundary at about 162/206/307 ms, complex Fresnel `7+I` at about 1.97 s, representative complex 2F1 continuation at about 0.98 s, and complex `li[-2+I]` at about 1.23 s. General real zeta uses roughly 42k certified refinements at representative non-integer points, while complex zeta remains around 0.45–0.55 s even with only about 8k–17k refinements, indicating substantial complex interval/Euler–Maclaurin kernel cost.

Gamma, ibeta, and the ordinary 2F1 precision sweeps showed comparatively smooth growth from 10 to 80 decimal digits and are lower-priority cliff targets in this pass. The next likely optimization targets are complex Fresnel near `|z|=8`, large positive 1F1, elliptic near-boundary evaluation, and a finer decomposition of zeta and complex continuation costs.
# 20. Removing exact-majorant cliffs in complex Fresnel and large positive 1F1

The performance-cliff runner showed that the first-order problem in both complex Fresnel and large positive-real `1F1` was not simply the number of series terms. The dominant cost came from **exact `Rational` proof bookkeeping growing on every term**.

The complex Fresnel Maclaurin terms were already stored as fixed-precision `ComplexInterval` values, but the Taylor-tail majorant was multiplied as an exact `Rational` at every iteration. Near `|z|=8`, the numerator and denominator inherited the dyadic Pi upper bound and grew repeatedly, pushing the representative 20-digit `7+I` case to about 1.97 s. The tail majorant now propagates as an outward-rounded positive `RealInterval`; the fixed `z^4*pi^2/4` complex factor is materialized once and reused; and the expensive tail certificate is tested every eight terms rather than on every term. When the tail is accepted, the interval upper endpoint is converted back to a Rational and used for the same explicit symmetric inflation as before, so the certification contract is unchanged.

For large positive-real `1F1`, the old backend retained the term and partial sum as exact `Rational` values. At `z=160`, the value is on the order of `e^z`, and Rational normalization, gcd work, and allocation dominated the roughly one thousand series iterations. For moderate positive parameters on `16<=z<=160`, all terms are positive, so the backend now reserves an argument-dependent working-precision guard and carries term/sum values as outward-rounded `RealInterval`s.

The old tail proof also waited until `n>=6|z|` to establish a fixed one-half future-ratio bound. For positive parameters, with `N=n+1`, every future ratio satisfies

```text
r_j = z (a+j) / ((b+j)(j+1))
    <= z (1+a/N) / (N+1) = Q.
```

Once `Q<1`, the omitted tail including the not-yet-added `next` term is bounded by `|next|/(1-Q)`. Testing this certificate every sixteen terms reduces the representative `z=160` CertifiedRefinement count from roughly 964 to roughly 513 while avoiding a new Rational-heavy proof loop.

A 2026-08-25 GCC Release / LTO-off remeasurement gave the following representative results. Absolute times remain machine/compiler dependent; the ratios and work counts are the important adoption signal.

| workload | previous warm value | optimized representative | note |
| --- | ---: | ---: | --- |
| `N[fresnelc[7+I],20]` | about 1.97 s | about 14 ms | existing `\|z\|<8` limit retained |
| `N[1F1[1/2,5/4,64],20]` | about 80 ms | about 2.5 ms | positive-real interval fast path |
| `N[1F1[1/2,5/4,128],20]` | about 619 ms | about 4.5 ms | same |
| `N[1F1[1/2,5/4,160],20]` | about 1.08 s | about 6 ms | existing `\|z\|<=160` limit retained |

The 1F1 `z=160` precision sweep remains roughly 5/6/7/10 ms at 10/20/40/80 digits, while complex Fresnel `7+I` remains in the tens of milliseconds across the same sweep. The intrinsic term-count increase toward the Fresnel boundary still exists, but the pathological representation-driven cliff is gone.

The positive 1F1 path was experimentally able to certify `z=512` in only tens of milliseconds, but widening only the exact-point backend would create a new asymmetry against finite-precision and complex inputs. This pass therefore keeps the public `|z|<=160` bounded-work policy unchanged; any future boundary expansion should first be audited with InformationEnclosure and the metamorphic certification suite.

# 21. Removing the elliptic F/E/Pi near-boundary series cliff

The performance-cliff audit measured roughly 162 ms for `ellipticF[1/2,9/10]`, 206 ms for `ellipticE[1/2,9/10]`, and 307 ms for `ellipticPi[9/10,1/2,1/3]` at 20 digits. AGM / Carlson symmetric forms were initially considered, but profiling showed that the first-order problem was not the Legendre series itself. It was the repeated growth of **exact `Rational` certification state** across hundreds of terms.

The old implementation retained `c_k`, `m^k`, the Pi combined coefficient `q_k`, and the tail power `r^k` as exact Rationals. As `r=max(|m|,|n|)` approached 9/10, numerator/denominator sizes, GCD normalization, and allocation costs grew with the number of terms. The series value and trigonometric quantities were already certified with `RealInterval`, so exact fractions were not required for these intermediate states.

The coefficient, parameter-power, Pi combined-coefficient, and tail-power recurrences now use outward `RealInterval` arithmetic at fixed working precision. The tail certificate only needs a rigorous upper bound; when it closes, the interval upper endpoint is converted back to a Rational and added as the same explicit symmetric error inflation used before. Tail certificates are checked every eight terms rather than every term.

When `|phi|<=Pi/2` can be proved, the tail bound also uses

```text
I_k(phi) = Integral[sin(t)^(2k), {t,0,phi}]
|I_k(phi)| <= |phi| sin(|phi|)^(2k)
```

instead of the previous unconditional `|I_k|<=|phi|`. The effective geometric ratio therefore becomes `rho = r sin(phi)^2`. At `phi=1/2`, `sin(phi)^2` is about 0.23, so even `r=0.9` produces `rho` near 0.21. If the amplitude is not provably inside `[-Pi/2,Pi/2]`, the implementation falls back to the original `rho=r` proof.

On 2026-08-26 with GCC Release / LTO off at 20 digits:

| workload | audit baseline | optimized | CertifiedRefinement |
| --- | ---: | ---: | ---: |
| `ellipticF[1/2,9/10]` | ~162 ms | ~2.3 ms | 799 -> 102 |
| `ellipticE[1/2,9/10]` | ~206 ms | ~2.6 ms | 799 -> 102 |
| `ellipticPi[9/10,1/2,1/3]` | ~307 ms | ~2.8 ms | 862 -> 110 |

Near `phi=3/2`, where `sin(phi)^2` is close to one, F/E are about 26 ms and Pi about 28 ms at 20 digits. This is a genuine convergence cost near the complete case rather than the former Rational-representation cliff. The performance-cliff runner now includes amplitude sweeps at `m/n=9/10` and 10/20/40/80-digit precision sweeps.

Carlson symmetric forms remain a strong future backend for Legendre elliptic integrals, but they were not needed to remove the present `phi=1/2, m/n->9/10` cliff: repairing the certification representation reduced it to a few milliseconds. Because a new Carlson backend would also require new proof machinery and amplitude/branch reduction policy, it is deferred until complete-case behavior or an expansion beyond the current `|m|<=9/10` / `|n|<=9/10` work boundary justifies that complexity.

# 22. Removing zeta Euler–Maclaurin certification / representation cliffs

The performance-cliff audit found that general non-integer real zeta spent roughly 42k `CertifiedRefinement` units at 20 digits, while representative complex zeta cases were around half a second. A separate high-precision cliff appeared near 80 digits, where the complex path reached roughly 2.5 seconds. Profiling showed that the Euler–Maclaurin formula itself was not the primary problem. The dominant cost came from **recomputing the same mathematical quantities through transcendental evaluation and growing exact-Rational proof state in the planner and finite Dirichlet sum**.

The real planner previously rebuilt the rising factorial, `(2k)!`, and `N^(-s-2k+1)` from scratch for each `k`. These now advance by recurrence, with a single certified `N^-2` factor reused between adjacent corrections. In the finite Dirichlet sum, positive integer bases use complete multiplicativity

```text
(ab)^(-s) = a^(-s) b^(-s)    (a,b>0)
```

so only prime bases require a fresh `Log -> Exp` certification; composites are assembled from already certified interval factors. Integer exponents retain the exact-power path, while half-integer exponents use exact integer powers plus a certified square root, avoiding a transcendental round trip for cases such as `n^(-3/2)=1/(n sqrt(n))`.

The complex high-precision cliff was even more representation-driven. The Bernoulli-remainder planner converted a high-precision dyadic Pi bound to an exact Rational and repeatedly formed `(2Pi)^(-2k)` with `rationalPower`. In the representative 80-digit case, nearly all wall time was spent in this planner. Since the planner needs a rigorous upper bound rather than an exact symbolic value, the recurrences for

```text
rising bound
N^(1-sigma-2k)
(2Pi)^(-2k)
```

now propagate as outward-rounded `RealInterval`s at fixed working precision. Only the accepted upper endpoint is converted back to a Rational for the existing explicit remainder inflation. The complex Euler–Maclaurin correction itself likewise advances the rising factorial, factorial, and `N^-2` power recurrence incrementally.

With the previous `K<=40` correction cap, the representative 80-digit case narrowly failed at `N=64` and jumped to `N=128,K=32`, doubling much of the expensive finite-sum work. Once the planner recurrence became cheap, a modest increase from **40 to 48 corrections** allowed the same case to close at `N=64,K=44`. This is a work-balance adjustment inside the existing certified backend, not an expansion of its mathematical domain.

Representative GCC Release / LTO-off measurements on 2026-08-26 were:

| workload | audit value | optimized |
| --- | ---: | ---: |
| `N[zeta[3/2],20]` | about 155 ms | about 10–12 ms |
| `N[zeta[12/5],20]` | about 275–321 ms | about 40–46 ms |
| `N[zeta[3/2+I],20]` | about 0.7–1.2 s | about 14 ms |
| `N[zeta[3/2+I],80]` | about 2.5 s even after the first pass | about 73 ms |

The 10/20/40/80-digit complex `zeta[3/2+I]` sweep is roughly 11/14/29/73 ms and no longer triggers the runner's 4x cliff threshold. General real points also avoid the previous roughly 42k-refinement plateau.

When this optimization was introduced, the certified region remained `Re(s)>1`, with the principal-branch convention, the `s=1` pole, and the `PrecisionInsufficient` / `CertifiedBackendUnsupported` classification unchanged. A later capability audit extended certified continuation through the critical strip and left half-plane; the complete-multiplicativity reuse introduced here remains limited to positive integer bases and still does not add any new complex-log branch transformation.

Across the Fresnel, 1F1, elliptic, and zeta passes, the recurring lesson is that **certified algorithms should distinguish quantities that need mathematical exactness from quantities that only need a rigorous enclosure or majorant**. Keeping the latter as normalized exact Rationals can create a representation cliff long before the underlying numerical algorithm becomes difficult. Future special-function work should profile planner / majorant / tail-certificate representation before adding a more complicated backend.

# 23. Removing complex 2F1 continuation / li Arg and Gamma representation cliffs

The remaining performance-cliff hotspots included roughly 0.96 s for the representative 20-digit complex `2F1` principal `1/z` continuation and roughly 1.03 s for `li[-2+I]`. Decomposition showed that the common dominant cost was not simply the number of Gamma evaluations or the `li=Ei(Log(z))` composition. It was the **exact-Rational certified `atan` series used by principal complex Log / Arg**.

After quadrant reduction, `enclosePrincipalArgument` reduces to expressions such as `atan(x/y)`. The certified atan backend then reduces to `|x|<=1/2` and uses the alternating series

```text
atan(x) = x - x^3/3 + x^5/5 - ...
```

with adjacent partial sums as the proof enclosure. The proof is sound, but interval division produces high-precision dyadic Rational endpoints. Repeatedly forming exact Rational powers of those endpoints grows denominator bit lengths term by term. Non-dyadic inputs therefore produced hundred-millisecond or second-scale Arg work even though simple cases such as `Log(1+I)` remained fast.

The input endpoint remains exact, but the atan series power, partial sum, and next term now propagate as outward-rounded `RealInterval`s at `precision+32` bits. The exact alternating-series value lies between two adjacent exact partial sums; since the two interval states enclose those partial sums, their hull is still a rigorous atan enclosure. Termination requires the upper endpoint of the next-term interval to be below the existing threshold, so the branch/Arg proof contract is unchanged while normalized-Rational denominator growth disappears.

Complex Ei tail-majorant bookkeeping was aligned with the same rule. The exact magnitude is not required; only a rigorous upper bound is. Its majorant now advances as an outward `RealInterval`, with the expensive tail certificate checked every eight terms. The atan change is the dominant li improvement, but this removes the same representation hazard from Ei itself.

Two additional costs were removed from the 2F1 continuation coefficient path:

1. Of the seven Gamma arguments in the connection coefficients, provably real arguments such as `a`, `b`, `b-a`, and `a-b` now use the real Gamma backend and are lifted back to `ComplexInterval` only afterward.
2. The complex Gamma Stirling remainder previously raised the post-shift real lower endpoint, often a high-precision dyadic Rational inherited from a non-dyadic input, to `lower^(2n-1)` exactly. The remainder proof only needs a positive lower bound, so the implementation now uses `floor(lower)` as a conservative integer bound. This can only enlarge the remainder bound while eliminating the huge dyadic denominator powers that dominated high-precision continuation cases.

On 2026-08-26 with GCC Release / LTO off, three-run warm measurements were:

| workload | audit baseline | optimized | improvement |
| --- | ---: | ---: | ---: |
| `N[2F1[3.4,5.6,4+I,4.6+2I],20]` | ~0.96 s | ~100 ms | ~9.6x |
| `N[li[-2+I],20]` | ~1.03 s | ~13.5 ms | ~76x |
| `N[li[-2],20]` | ~575 ms | ~9 ms | ~64x |

The corresponding precision sweeps are approximately:

```text
2F1 continuation : 10/20/40/80 digits ~= 77 / 100 / 177 / 485 ms
complex li        : 10/20/40/80 digits ~= 10 / 13.5 / 19 / 36.5 ms
```

The 2F1 continuation remains intrinsically more expensive than the ordinary unit-disk Gauss series because the connection formula contains seven Gamma factors, two inner 2F1 evaluations, and two principal powers. That fixed cost is now distinct from the former Rational/Arg cliff. Growth from 40 to 80 digits remains below the runner's 4x threshold, so a more complex dedicated Gamma-ratio backend is not yet justified by the measured crossover.

The public 2F1 bounded-work region, principal `1/z` connection conditions, the principal `Log -> Ei` definition of li, and all branch-cut / pole / `PrecisionInsufficient` classifications are unchanged. The performance-cliff runner now includes 10/20/40/80-digit sweeps for both paths.


### Removing the former 2F1 9/10 work boundary (2026-08-26)

The first performance-cliff audit restricted the certified `2F1` Gauss series to `|z|<=9/10`. This was not a mathematical convergence boundary: it was a work policy introduced because the former exact-`Rational` certification state grew pathologically as `z` approached one. After moving the Gauss-series term and tail proof state to fixed-working-precision outward intervals, that threshold was re-audited.

On GCC Release / LTO-off at 20 digits, representative `2F1[1/2,1/3,5/4,z]` cases certify in roughly 0.07 s at `z=19/20`, 0.18 s at `49/50`, and 0.41 s at `99/100`. At `999/1000` the cost still rises to roughly 5.3 s, so approaching the unit circle retains a genuine convergence cost; however, that growth is continuous and there is no algorithmic discontinuity at 9/10.

The fixed 9/10 threshold is therefore removed. The Gauss series now attempts the full convergence region `|z|<1`, while the term cap and shared `EvaluationBudget` provide bounded work. A finite-precision magnitude interval straddling `|z|=1` is `PrecisionInsufficient`; a provable `z=1` point uses Gauss summation when `Re(c-a-b)>0`; other `|z|=1` points remain `CertifiedBackendUnsupported` until their parameter-dependent boundary formulas are implemented. The principal `1/z` continuation for `|z|>1` is unchanged.

This establishes a broader rule: **hard thresholds introduced as performance workarounds should be re-measured after the underlying representation cliff is removed, and should be moved back to the mathematical boundary whenever bounded-work machinery is sufficient.**

## High-precision contraction for complex Lambert W branches

After arbitrary integer branches were added, `N[lambertw[2,1],100]` exposed a cliff: rectangle inclusion for `Log_k(z)-Log(w)` repeatedly recomputed a certified complex logarithm until the box reached 100-digit width, taking about 2.7 s. Candidate generation now uses Newton iteration on `w exp(w)-z=0` instead of the linearly convergent logarithmic fixed point. The proof remains the branch-explicit Banach disk inclusion for the logarithmic map; the Newton candidate is only a search aid, and the fast path returns only when the certified disk radius satisfies the requested-width gate.

On GCC Release / LTO-off, the 100-digit `lambertw[2,1]` case falls from about 2.7 s to about 0.07 s. The requested-width gate is essential: returning an enclosing disk merely because it proves existence would otherwise expose too few digits.


# 24. Removing Complex Root symmetric-seed and re-isolation cliffs

A re-audit of `N[root[...,k,Complex],p]` found a substantially larger cliff in Complex Root than in the real Sturm path. The worst cases depended more on polynomial symmetry than on degree itself.

The former Durand-Kerner candidate generator placed all roots on one radius with a perfectly regular angular spacing. That is normally adequate, but for rotationally symmetric polynomials such as `x^n-a` it can keep the iteration trapped in a symmetric trajectory. `root[{-2,0,0,0,1},1,Complex]` could exceed ten seconds, with `x^8-2` and `x^10-2` also reaching multi-second or timeout behavior.

The candidate seeds now receive small deterministic angular and radial perturbations. Seeds are only search aids: every accepted root disk is still proved to contain exactly one root by the existing exact-Rational Rouche test, so this does not weaken the mathematical contract.

Exact Root calls also already carried a canonical `AlgebraicNumber` and isolating disk in their internal `algebraicValue` cache. The certified evaluator previously ignored that cache and rebuilt the ComplexAlgebraicNumber from polynomial coefficients. `N[Out[-1],p]` and equivalent paths now reuse the exact-evaluation cache directly.

`ComplexAlgebraicNumber::refined` formerly called `isolateComplexDisks` again even when only one selected root needed additional bits, repeating candidate generation, certification, and ordering for every root. Refinement now starts from the existing certified disk, Newton-refines only its center, and accepts a local disk only when

```text
new disk is certified unique by Rouche
AND
new disk is contained in the old certified disk
```

so the selected root identity is preserved. Full re-isolation remains the fallback when the local proof cannot close.

A second representation cliff appeared inside the local Rouche proof. Converting the full-working-precision Newton center directly to an exact Rational can preserve a mathematically zero component as a tiny dyadic with a huge denominator. For the purely imaginary selected root of `x^8-2`, that meaningless component made exact Taylor expansion expensive enough to dominate refinement. The proof center is now quantized to half the working precision, and components smaller than that absolute proof resolution are normalized to zero. The resulting disk is still re-certified by Rouche, so this changes only the choice of proof center, not the accepted enclosure.

On 2026-08-27, GCC Release / LTO off, two-run averages from `--algebraic-root-cliffs 2` were:

| case | initial isolation | refine 80 bit | refine 320 bit |
| --- | ---: | ---: | ---: |
| `x^2-2` | 4.7 ms | 0.8 ms | 1.0 ms |
| `x^4-2` | 13.9 ms | 3.2 ms | 9.7 ms |
| `x^5-2` | 46.8 ms | 5.9 ms | 27.1 ms |
| `x^8-2` | 68.7 ms | 10.1 ms | 47.3 ms |
| `x^10-2` | 226.8 ms | 52.7 ms | 169.4 ms |
| `x^12-2` | 332.3 ms | 95.2 ms | 286.8 ms |
| `x^16-2` | 696.9 ms | 163.0 ms | 771.9 ms |
| generic dense degree 5 | 18.8 ms | 7.2 ms | 26.0 ms |
| close-root degree 4 | 28.4 ms | 0.3 ms | 2.1 ms |

The former >10 s `x^4-2` pathology is gone, and behavior is continuous through degree 16. At this point in the audit, higher-degree isolation still measured roughly 2.2 s for `x^20-2`, 4.4 s for degree 24, and 15.7 s for degree 32, so `maximumSupportedAlgebraicDegree=64` was temporarily retained as a backend/resource safety cap. A later high-degree Rouché and candidate-generation audit removed the underlying cliff and raised the measured budget to 96.

`mmCal.Benchmarks --algebraic-root-cliffs [iterations]` now measures initial all-root isolation separately from selected-root 80/320-bit local refinement and keeps symmetric sparse, generic dense, and close-root cases as permanent regression probes.

# 25. Removing the fixed `|z|<=160` limit from `1F1`

After the Complex Root audit, the remaining certified-special-function hard limits were reviewed again. The `|z|<=160` rule in `hypergeometric1F1` was not a mathematical boundary; it was a policy threshold left over from the former exact-`Rational` series representation. The earlier large-positive-real optimization had already moved terms and partial sums to outward `RealInterval`s, but the negative-real and complex paths still needed a common proof before the fixed threshold could be removed.

The large-`|z|` real path now uses a sign-independent future-term majorant. Once `N>=2|a|,2|b|`, all later term ratios satisfy

```text
|t_(j+1)/t_j| <= 3 |z| / (N+1)
```

so when this bound `Q` drops below one, the tail including the next not-yet-added term is enclosed by `|next|/(1-Q)`. Terms and partial sums remain guarded dyadic intervals, so large negative-real cancellation does not reintroduce exact-Rational denominator growth.

The complex backend uses the same majorant. Its temporary `|z|<=8000` cancellation-guard cutoff has also been removed. Real and complex `1F1` therefore have no fixed magnitude threshold. Instead, the backend checks whether the conservative `2|a|`, `2|b|`, and `3|z|` tail-certification start can fit within the 250000-term algorithm budget, while the shared `EvaluationBudget::CertifiedRefinement` remains the request-wide safety net. Bounded work is thus determined by estimated work, not by a magic argument magnitude.

Finite-precision arguments use the same path without collapsing their `InformationEnclosure`. For example, `N[hypergeometric1F1[1/2,5/4,N[161,5]],20]` now evaluates beyond the former boundary, but does not manufacture twenty digits from a five-digit input; its resulting precision remains about two digits.

On 2026-08-27, GCC Release / LTO off, three-run warm values from `--performance-cliffs 3` were:

| workload | warm |
| --- | ---: |
| `N[1F1[1/2,5/4,64],20]` | about 2.1 ms |
| `N[1F1[1/2,5/4,160],20]` | about 6.0 ms |
| `N[1F1[1/2,5/4,161],20]` | about 5.0 ms |
| `N[1F1[1/2,5/4,256],20]` | about 9.7 ms |
| `N[1F1[1/2,5/4,512],20]` | about 26.6 ms |
| `N[1F1[1/2,5/4,1000],20]` | about 98 ms |
| `N[1F1[1/2,5/4,-512],20]` | about 27 ms |
| `N[1F1[1/2,5/4,256+I],20]` | about 53 ms |
| `N[1F1[1/2,5/4,-256+I],20]` | about 44 ms |

The `z=512` precision sweep was about 25/30/28/32 ms at 10/20/40/80 digits, while complex `256+I` was about 41/43/65/117 ms. No four-times cliff remained at the former 160 boundary or in the tested high-precision range.

The fixed `|z|<=160` rule is therefore removed rather than replaced by a larger magic number. If genuinely large arguments expose a new asymptotic cliff later, Kummer transformation or a dedicated asymptotic backend should be considered; the current measurements do not justify restoring a magnitude cutoff.

# 27. Removing the fixed real Ei / Si / Ci boundary with certified asymptotics

After removing the fixed 1F1 boundary, the remaining real `Ei/Si/Ci` boundary at 96 was re-audited. Once the Taylor state had been moved to guarded intervals, arguments just above 96 already evaluated in milliseconds, so the threshold had become a legacy workaround for exact-Rational cancellation and representation cost. Merely extending the Taylor range was still the wrong design: `Ci[1000]` drove the `gamma+log(x)` cancellation path to roughly 2000-bit EulerGamma work and exposed a different backend limit. The fix is therefore a distinct large-argument backend rather than another magnitude threshold.

For positive real `Si/Ci`, the implementation uses the DLMF 6.12 auxiliary expansions

```text
f(x) ~ 1/x (1 - 2!/x^2 + 4!/x^4 - ...)
g(x) ~ 1/x^2 (1 - 3!/x^2 + 5!/x^4 - ...)

Si(x) = Pi/2 - f(x) cos(x) - g(x) sin(x)
Ci(x) = f(x) sin(x) - g(x) cos(x).
```

On the positive real axis the remainders of `f` and `g` are bounded by the first neglected terms and have their signs. The backend therefore propagates the asymptotic terms with outward `RealInterval` arithmetic and adds the first omitted terms as one-sided remainder intervals. If the terms start increasing before the requested width is certified, it falls back to the convergent Taylor path. This removes EulerGamma completely from large-argument Ci value construction.

Positive real `Ei` uses `Ei(x) ~ exp(x)/x (1 + 1!/x + 2!/x^2 + ...)` with the DLMF 6.12.2 remainder bound `(1+chi(n+1))*nextTerm`. For integer indices the required `chi` values are updated from `chi(2)=2`, `chi(3)=3 Pi/4`, and `chi(t+2)=chi(t)(t+2)/(t+1)`, avoiding a separate Gamma-ratio evaluation. Certifying the bracket before multiplying by `exp(x)/x` preserves relative/significant-digit accuracy even for enormous values. Negative real Ei keeps the existing certified `Ei(-x)=-E1(x)` asymptotic path.

On GCC Release / LTO-off, `--performance-cliffs 2` gives representative 20-digit warm values of about 6.0 ms for `Ei[96]`, 1.8 ms for `Ei[256]`, 1.5 ms for `Ei[512]`, 1.4 ms for `Si[1000]`, and 2.0 ms for `Ci[1000]`. `Si[10000]`, `Ci[10000]`, and `Ei[10000]/exp[10000]` also complete quickly with certified results. The former fixed 96 threshold is not replaced by another magnitude constant: backend choice is controlled by whether the asymptotic remainder reaches the requested precision, the Taylor term cap, and the shared `EvaluationBudget`.

The audit also exposed the existing `N[Ei[-64],20] -> 0.0` loss of significant digits. Scaling the negative-Ei Taylor tail target with the rigorous lower bound `E1(x)>=exp(-x)/(x+1)` preserves values such as `-2.46796855945...e-30` instead of rounding them away.

# 28. Removing the fixed `|z|<8` boundary from complex Fresnel

After eliminating the exact-majorant representation cliff in complex Fresnel, the historical `|z|<8` rule was re-audited as a capability limit. Simply removing the boundary lets the Maclaurin path certify `8+I` in about 0.02 s, `10+I` in about 0.08 s, and `20+I` in about 1.2 s, but `32+I` still exceeds 8 s because the oscillatory series develops a genuine cancellation cliff. Thus 8 itself was stale, but large near-axis complex arguments still needed a different certified representation.

Rather than moving the cutoff, mmCal now uses the DLMF 7.12 complex Fresnel auxiliary-function asymptotics,

```text
C(z) = 1/2 + f(z) sin(Pi z^2/2) - g(z) cos(Pi z^2/2)
S(z) = 1/2 - f(z) cos(Pi z^2/2) - g(z) sin(Pi z^2/2).
```

For `|arg z|<Pi/8`, the `f/g` remainders are bounded in magnitude by the first neglected terms. The implementation never guesses this sector from an approximate angle: it uses the stronger directly provable interval condition `Re(z)>0` and `4|Im(z)|<=Re(z)`, which lies strictly inside that DLMF sector. Quarter-turn identities `C(i z)=i C(z)` and `S(i z)=-i S(z)` map neighborhoods of all four coordinate axes into the same certified wedge. Inputs that do not prove wedge membership fall back to the guarded entire Maclaurin series.

The `f/g` terms are propagated as `ComplexInterval`s. First-neglected-term magnitudes are carried through the certified complex sine/cosine factors before accepting the requested output width; if the asymptotic terms begin increasing before the proof closes, evaluation falls back to the series rather than extending a divergent asymptotic expansion.

On 2026-08-27, GCC Release / LTO off, three-run warm values are roughly 10 ms for `N[fresnelc[32+I],20]`, 12 ms for `1+32I`, and 13 ms for diagonal `8+8I`. The `32+I` 10/20/40/80-digit sweep is about 8/10/16/35 ms. Finite-precision probes also retain their input information limits, so the new dispatch does not recover hidden exact centers.

The fixed complex-Fresnel magnitude cutoff is therefore removed. Bounded work is now controlled by the asymptotic remainder certificate, the Maclaurin term cap, and the shared `EvaluationBudget`.



# 29. Removing the elliptic 9/10 capability boundary

`ellipticF/E/Pi` retained fixed `|m|<=9/10` limits, with `Pi` also requiring `|n|<=9/10`. After moving the series certification state away from exact-Rational growth, inputs immediately below that boundary were already in the millisecond range, so the cutoff had become a stale performance workaround.

A certified Carlson symmetric-form backend now evaluates the Legendre forms through `RF/RD/RJ`, using the DLMF 19.25 transformations and the DLMF 19.26 duplication identities. Positive-real monotonicity encloses the Carlson residuals. The `RC` correction required by `RJ` switches to a regular near-equal series where the elementary representation would suffer interval cancellation. Exact `q Pi` amplitudes are period-reduced from the rational coefficient rather than by dividing independent interval approximations to `Pi`.

Dispatch no longer uses a fixed 0.9 parameter cutoff. The Legendre series remains the fast path when the certified effective tail ratio is `<=9/10`; otherwise evaluation falls back to Carlson. This removes the artificial performance discontinuity at the former capability boundary.

Representative 20-digit warm timings on GCC Release / LTO off, 2026-08-28:

| workload | warm |
| --- | ---: |
| `ellipticF[1/2,0.90]` | about 2.2 ms |
| `ellipticF[1/2,0.95]` | about 2.6 ms |
| `ellipticF[1/2,0.99]` | about 2.4 ms |
| `ellipticF[1/2,0.999]` | about 2.5 ms |
| `ellipticPi[0.99,1/2,0.99]` | about 3.0 ms |
| `ellipticF[3/2,0.99]` | about 26.7 ms |
| `ellipticE[3/2,0.99]` | about 53.4 ms |
| `ellipticPi[0.99,3/2,0.99]` | about 75.2 ms |

The 10/20/40/80-digit Carlson sweeps are roughly 16/27/54/143 ms for F, 32/53/115/317 ms for E, and 46/75/155/414 ms for Pi, with continuous growth as precision increases.

The capability boundary is now mathematical rather than a fixed parameter magnitude. Local real values with `m>1` or `n>1` are certified when the reduced real integration path is proved to stay before the branch point or pole. Paths crossing a singularity and general complex continuation remain unsupported. Exact `m=1` for `ellipticE` uses its finite real degeneration.

# 30. Removing the positive-order polylog 49/50 boundary with DLMF continuation and a near-one backend

Positive-order `polylog` still carried a fixed `|z|<=49/50` bounded-work cutoff inherited from the old exact-Rational series. Removing the cutoff for measurement showed that, even within the mathematical convergence region `|z|<1`, keeping `z^k` as a growing exact Rational was a major part of the near-unit-point cost.

The direct series now propagates

```text
t_(k+1) = t_k z (k/(k+1))^s
```

with outward `RealInterval` / `ComplexInterval` recurrence. The tail remains certified by the uniform future-ratio bound `|z|`, giving `|next|/(1-|z|)`. This removes exact numerator/denominator growth, although the natural slow convergence of cases such as `Li_2(0.999)` remains if the direct series is used alone.

For `Li_2`, the principal DLMF 25.12.3, 25.12.4, and 25.12.6 connection formulas are now available. A transformation is accepted only when the complete input enclosure proves that the relevant branch cut is avoided and the transformed argument gives a useful magnitude contraction.

```text
Li_2(z) + Li_2(z/(z-1)) = -1/2 Log(1-z)^2
Li_2(z) + Li_2(1/z) = -Pi^2/6 - 1/2 Log(-z)^2
Li_2(x) + Li_2(1-x) = Pi^2/6 - log(x)log(1-x), 0<x<1
```

This covers the positive-real unit-point neighborhood, the negative real axis, `z=I`, and selected `|z|>1` regions away from the cut by mapping them to smaller arguments. Exact positive-real `z>1` lies on the principal branch cut itself, so mmCal deliberately keeps such a point unevaluated rather than inventing an upper- or lower-side boundary value.

For higher positive integer orders, `z=1` is finite but the direct series still became second-scale near `z=0.999`. Taking the positive-integer `s -> n` limit of the noninteger DLMF 25.12.12 formula gives, for positive real `0<z<1` and `mu=log(z)<0`,

```text
Li_n(exp(mu))
 = mu^(n-1)/(n-1)! (H_(n-1) - log(-mu))
   + sum_{k>=0, k!=n-1} zeta(n-k) mu^k/k!
```

and this is now used as a certified near-one backend. Negative-integer zeta values reduce exactly through

```text
zeta(1-2r) = -B_(2r)/(2r).
```

The remainder is not accepted merely because terms appear small. Using

```text
|B_(2r)| = 2 (2r)! zeta(2r)/(2Pi)^(2r),  zeta(2r) < 2,
```

a majorant is constructed whose successive tail ratio is at most `(|mu|/(2Pi))^2`; the omitted tail is then enclosed geometrically. This supplies an explicit remainder certificate for the integer-limit expansion.

The near-one path is an internal fast-path choice, not a capability boundary. GCC Release / LTO-off crossover measurements currently select it for orders 3 through 12 with `|mu|<=1/20`; all other points fall back to the existing unit-disk series. At high order the `1/k^s` factor already makes the direct series cheaper. Finite-precision real inputs use the same expansion on their full `RealInterval`, never by recovering a hidden midpoint.

Representative 20-digit measurements on 2026-08-28, GCC Release / LTO off, are:

| workload | previous audit | new backend |
| --- | ---: | ---: |
| `Li_2(0.999)` | ~1.2 s | ~2.9 ms |
| `Li_3(0.999)` | ~0.96 s | ~11 ms |
| `Li_4(0.999)` | ~0.80 s | ~13 ms |
| `Li_8(0.999)` | ~0.24 s | ~24 ms |
| `Li_2(I)` | outside series region | ~35 ms |
| `Li_2(-2)` | outside series region | ~5 ms |
| `Li_2(2+I)` | outside series region | ~36 ms |

The 10/20/40/80-digit `Li_2(0.999)` sweep is about 2.3/2.9/4.2/7.3 ms, with continuous scaling. `Li_3(0.999)` at 80 digits is about 0.09 s including CLI startup. `N[polylog[3,N[999/1000,8]],20]` fell from roughly 2.9 s on the old interval series to about 0.01 s on the near-one interval path while still displaying only `1.2004154`, so the optimization does not manufacture precision beyond the eight-digit input.

The fixed `|z|<=49/50` cutoff is therefore removed. Bounded work is now defined by the **unit-disk series term cap, shared `EvaluationBudget`, branch proofs for the DLMF connection formulas, and the Bernoulli remainder proof for the near-one expansion**.


# 31. Removing the fixed 512/128 boundaries from complex Ei/Ci

Certified complex `Ei` and `Ci` retained historical bounded-work limits `|z|<=512` and `|z|<=128`. Simply deleting those checks exposed the underlying Taylor cancellation cliff: at 20 digits, `Ei[513+I]` took about 0.7 s, `Ci[129+I]` about 1.5 s, and `1000+I` cases could exhaust the certified-refinement budget.

A large-argument certified kernel now evaluates the DLMF 6.12.1 expansion

```text
E1(z) ~ exp(-z)/z (1 - 1!/z + 2!/z^2 - ...).
```

Terms are propagated as `ComplexInterval`s. In `|arg z|<=Pi/2` the first omitted term bounds the remainder; in the left half-plane the DLMF sector bound adds the `csc(|arg z|)` factor. If the requested whole-complex width cannot be certified before the asymptotic terms begin increasing, the result is rejected and evaluation falls back to the convergent series.

Complex `Ei` uses the principal identity `Ei(z)=-E1(-z)+Log(z)-Log(-z)`, concentrating branch correction in the principal Log implementation. Because the positive real axis is not an Ei cut even though `-z` lies on the E1 cut, a finite-precision enclosure that straddles zero imaginary part there is handled from a certified real-Ei anchor plus a vertical-path bound using `Ei'(z)=exp(z)/z`. This preserves the input width without recovering a hidden exact zero.

Complex `Ci` uses `Ci(z)=-1/2(E1(i z)+E1(-i z))` in the right half-plane, reflects the left half-plane with `Ci(z)=Ci(-z)+Log(z)-Log(-z)`, and uses the Chi degeneration on the pure imaginary axis. When optimal asymptotic truncation is insufficient at high precision, the guarded series remains a fallback; its tail-majorant bookkeeping has also been changed from growing exact Rationals to outward `RealInterval` recurrence.

On 2026-08-28, GCC Release / LTO off, representative 20-digit timings are about 13.2/8.2/9.4/9.2 ms for `Ei[128+I]`, `Ei[512+I]`, `Ei[513+I]`, and `Ei[1000+I]`; about 49.7/44.2/45.7/41.9/18.4 ms for `Ci[120+I]`, `Ci[128+I]`, `Ci[129+I]`, `Ci[140+I]`, and `Ci[1000+I]`. The `Ei[1000+I]` 10/20/40/80-digit sweep is about 7.3/8.5/12.3/28.6 ms, while `Ci[1000+I]` is about 11.7/19.6/40.9/114.1 ms. A fallback-heavy `Ci[140+I]` remains around 0.04/0.28/0.36 s at 20/50/100 digits.

The fixed 512/128 limits are therefore removed. Bounded work is now governed by **branch-cut proofs, the E1 asymptotic remainder certificate, guarded-series term caps, and the shared `EvaluationBudget`**.


## Local contraction at the Lambert W `-1/e` branch point

The general complex Lambert W backend uses Newton candidates for `w exp(w)-z=0` and a Banach proof for `Log_k(z)-Log(w)`, but the derivative degenerates at `w=-1`, making boxes difficult to close near `-1/e`. Using the square-root local structure of DLMF 4.13.9_1, set `u=W+1`, `q=e z+1` and write `q=u^2 A(u)/2` with `A(0)=1`. The local proof fixes `p=sqrt(2q)` once and iterates `u=±p/sqrt(A(u))`; `A` and `A'` are enclosed by outward complex interval series, which provides both mapping inclusion and a uniform contraction bound.

An earlier attempt to reevaluate `sqrt(2q/A(u))` at every iteration was rejected because rectangular complex-sqrt dependency around its cut prevented the enclosure from shrinking below roughly `1e-8`. Fixing `p` and evaluating `1/sqrt(A)` with a certified binomial series in `A-1` removes that dependency. `W_0` follows principal `p`, `W_-1` uses `-p` on the upper-side local branch, and `W_1` uses the symmetric lower-side local branch. The first three Puiseux terms are used only for candidate generation; final acceptance is by the contraction proof.

For extreme offsets, evaluating all of `z` at finite precision before forming `e z+1` loses the distance to the branch point through cancellation. The evaluator now simplifies `z+1/E` while the input is still exact and passes that offset directly to the local kernel. Termination is likewise based on component-wise relative width of the nonzero components of `u=W+1`, rather than on absolute width of the full value `W≈-1`. The Puiseux seed is refined by Newton iteration on `q=u^2 A(u)/2`, then a dynamically small proof box is constructed from the candidate accuracy. The Newton value is never accepted as the proof result; final acceptance still requires `T(B) subset B` together with a uniform contraction bound.

The local `A(u)`, `A'(u)`, and `A(u)^(-1/2)` tail majorants also no longer accumulate growing exact Rationals. They certify a geometric tail from the outward interval norm of the current term and a uniform term-ratio bound. On the 2026-08-28 GCC Release / LTO-off runner, representative 100-digit timings are about 1.51 s for `-1/e+i*10^-2`, 0.31 s for `-1/e+i*10^-12`, and 0.031 s for `-1/e-i*10^-160`. The real point `-1/e+10^-12` is about 0.25 s including CLI startup (`W_-1` about 0.29 s), so the shallow and extreme branch-point regimes no longer create multi-second timeouts.

The audit also exposed a common performance bug outside Lambert W. `algebraicBinary` attempted to convert a huge complex Rational right operand to `AlgebraicNumber` even when the left operand had already been proved nonalgebraic. Short-circuiting after the failed left conversion avoids unnecessary algebraic-field construction for expressions such as `x-I/10^160` and `1/E-I/10^160`.

# 32. Newton-polygon multi-radius candidates for Complex Root isolation

Complex Root certification already separated approximate candidate generation from exact proof: Durand-Kerner supplies centers, Newton refines each center locally, and an exact Rational Rouché test proves that the final disk contains exactly one root. The remaining weakness was the initial geometry. Every root candidate started near one perturbed Cauchy radius, which is reasonable for roots of comparable magnitude but pathological when one polynomial contains several widely separated root-radius groups.

Candidate initialization now computes the upper Newton polygon of the nonzero coefficient points

```text
(k, log2(|a_k|)).
```

For consecutive polygon vertices `k0 < k1`, `m=k1-k0` candidates are placed on a circle with approximate radius

```text
r = (|a_k0 / a_k1|)^(1/m).
```

The logarithms and radii are deliberately approximate: they affect only deterministic starting points. The existing exact Rouché certification remains the sole acceptance criterion, so an inaccurate polygon caused by finite candidate arithmetic can at worst cost performance and cannot certify a wrong root. Polynomials with a zero constant term or unusable approximate geometry fall back to the previous Cauchy-radius initialization.

This follows the same initialization principle documented for FLINT/Arb `acb_poly_find_roots`, which uses Newton-polygon circles before Durand-Kerner and then rigorously validates the roots. mmCal retains its own exact Rational disk proof and deterministic `Re(z)+Pi Im(z)` ordering.

A permanent benchmark case uses

```text
(x^8 - 2^80)(x^8 - 2^-80)
= x^16 - (2^80 + 2^-80)x^8 + 1,
```

whose roots split into two groups of eight with radii `2^-10` and `2^10`. On the 2026-08-28 GCC Release / LTO-off build, the old single-Cauchy-radius initialization took about 10.79 s for initial all-root isolation, while the Newton-polygon initialization takes about 0.36 s. The existing `x^16-2` case remains around 0.91 s, so the optimization removes a scale-separation pathology without regressing the symmetric single-radius family.

The Newton-polygon seed change alone did not remove the degree-64 Complex Root backend limit, because certification, pairwise disk separation, and deterministic ordering still grew materially with degree. A later audit replaced the dominant Rational Rouché proof path and improved Durand-Kerner termination, after which the measured defining-polynomial budget was raised to 96. The `--algebraic-root-cliffs` runner keeps the symmetric `x^n-2` family, the two-radius stress case, and high-degree Solve cases.

# 33. Direct finite-Fourier integration for Advanced Integration

Case-level timing of `runAdvancedIntegrationTests` showed that the small-degree `sin[u]^m cos[u]^n` grid repeatedly spent roughly 0.5--1.1 s per integration case. The derivative-back checks were only a few milliseconds; the cost was in `integrate[...]` itself. The old path first materialized a finite Fourier Expr, ran the general simplifier, and then re-entered `integrateCore`, after already paying for generic reverse-chain, Weierstrass, and integration-by-parts candidate searches.

`trigonometric_polynomial` can now return a structured `TrigFourierExpansion` containing the common argument, exact Rational coefficient, frequency, and sine/cosine kind. For same-argument integer sine/cosine powers of total degree at least three, integration attempts this representation before generic candidate exploration. If the common derivative is provably nonzero and independent of the integration variable, it is computed once and the primitive is built directly from

```text
c cos(k u) -> c sin(k u)/(k u')
c sin(k u) -> -c cos(k u)/(k u')
constant   -> constant*x.
```

The existing `TrigArgument` scale is retained, so Rad/Deg/Grad session semantics do not change. Nonlinear arguments are not forced through the fast path. For example, `sin[2x^2]^4` falls back to the former Fourier-Expr/general-integrator path and therefore still reaches the Fresnel transformation. Existing compact square and `sin[x]cos[x]` forms also retain their earlier routes.

On 2026-08-28, GCC Release / LTO off, `runAdvancedIntegrationTests` fell from roughly 14.9 s to about 2.2 s while retaining every test, including the degree-256 Fourier case and nonlinear Fresnel composition.

# 34. Algebraic field persistence and Lambert W disk certificates in CalculusKnowledge

The two dominant `runCalculusKnowledgeTests` costs were a chained algebraic power equivalent to `(sqrt(2)+sqrt(3))^3` at roughly 1.9 s, and `N[lambertw[-1,-1/E-I/10^8],20]` at roughly 2.0--2.3 s.

For the algebraic case, the current irreducibility prover could not prove that

```text
x^4 - 10x^2 + 1
```

is irreducible over Q, so the persistent number-field representation was not attached and later powers fell back to resultant/minimal-polynomial reconstruction. The existing modular criterion is only sufficient: a V4-type quartic need not remain irreducible modulo any of the small primes being tried.

For a monic even quartic `x^4+b x^2+d`, write any monic quadratic factorization as

```text
(x^2+p x+q)(x^2-p x+s).
```

The linear coefficient gives `p(s-q)=0`. Reducibility is therefore completely decided by either `b^2-4d` being a rational square (`p=0`), or `d` being a rational square and `2q-b` being a rational square for one of `q=+-sqrt(d)` (`q=s`). This exact proof is now used after the modular sufficient tests.

Integer powers of a real `AlgebraicNumber` that already has field coordinates are also computed by binary exponentiation directly in the `AlgebraicElement` Q(theta) coordinates. Intermediate powers are no longer materialized as visible `Root[minpoly,k]` values; minimal-polynomial construction and root re-identification occur only once at the end. The `(sqrt(2)+sqrt(3))^2/^3` equivalents consequently fell from about 0.6/1.9 s to roughly 0.14--0.17 s.

For nonlocal Lambert W branches, the logarithmic backend already obtained a valid Banach disk around the Newton candidate. It nevertheless discarded that certified disk when its width missed an internal guarded target by only a few digits, then paid for as many as 256 interval-Log rectangle contractions. The backend now returns any disk that has proved mapping inclusion and contraction. The outer `CertifiedEvaluator` already decides whether the returned enclosure carries enough precision and requests a higher-precision reevaluation when necessary, so repeating that target requirement inside the backend was redundant. The mathematical certificate is unchanged.

`N[lambertw[-1,-1/E-I/10^8],20]` drops from about 2.3 s to about 0.05 s. After these changes the largest per-case intervals in CalculusKnowledge are around 0.31 s, mainly complex Root isolation, rather than multi-second outliers. The whole group fell from about 7.8 s to about 3.2 s.

No coverage was removed from `runAdvancedIntegrationTests` or `runCalculusKnowledgeTests`. Final 2026-08-28 GCC Release / LTO-off validation recorded 2820/2820 internal C++ regressions, 710/710 relevant black-box cases, and 1000/1000 Certification Boundary Fuzzer cases (59 classification + 36 metamorphic probes).

# 35. Modular irreducibility proof and batch root canonicalization for general polynomial `solve`

Auditing `solve[x^16+x+1==0,x]` found two structural cliffs around Complex Root isolation rather than in Durand-Kerner or the Rouché certificates themselves.

First, the former `provenIrreducibleOverQ` only used the sufficient test that one small-prime reduction remain irreducible at the full degree. The polynomial `x^16+x+1` is irreducible over Q, but every reduction in the former small-prime set split, so the prover gave up and degree-16 Kronecker factor search was entered.

The prover now uses several square-free good-prime reductions jointly. For each prime, exact Frobenius/gcd computations

```text
gcd(f, x^(p^k)-x)
```

recover the number of irreducible mod-p factors of each degree. A Q-factor of degree `d` must preserve degree at every good prime and therefore reduce to a subset of those finite-field factors whose degrees sum to `d`. mmCal builds the possible proper-factor degree set `1..floor(n/2)` for each prime and intersects the sets. An empty intersection proves irreducibility over Q.

For `x^16+x+1`, representative reductions give

```text
mod 2 : factor degrees 8 + 8   -> possible proper degree {8}
mod 3 : factor degrees 1 + 15  -> possible proper degree {1}
```

so no Q-factor degree can exist. This is not a probabilistic irreducibility guess: only exact good-prime information is used, and primes for which square-freeness or degree preservation cannot be established are skipped.

Second, Complex `solve` previously called `ComplexAlgebraicNumber::isolateAll` once and then called `ComplexAlgebraicNumber::create` separately for every returned root in order to obtain a canonical minimal-polynomial Root. `create` re-isolates every root of the same defining polynomial, making a degree-`n` solve effectively pay for all-root isolation `1+n` times.

`ComplexAlgebraicNumber::canonicalizeAll` now canonicalizes an already certified all-root set in one batch. If the defining polynomial is proved irreducible, the existing disks and global ordering are reused directly. For reducible polynomials up to the current degree-16 minimal-polynomial reduction limit, factorization is performed once, each proven irreducible factor is isolated once, and intersection with the original ordered disks identifies the minimal factor and the factor-local root index. If that correspondence cannot be proved uniquely, the solver falls back to the former per-root path. Above the minimal-polynomial reduction limit, where the old `create` path would retain the original polynomial anyway, the already certified disks are reused instead of repeating all-root isolation.

Representative 2026-08-28 GCC Release / LTO-off timings are:

| workload | wall time |
| --- | ---: |
| `solve[x^6-3x^5-x^4+2x^3+2x^2-2x-1==0,x]` | about 0.046 s |
| `solve[(x^4-2)(x^4+1)==0,x]` | about 0.17--0.20 s |
| `solve[x^16+x+1==0,x]` | about 1.38--1.45 s |
| `solve[x^18+x+1==0,x]` | about 2.37 s |
| `solve[x^20+x+1==0,x]` | about 3.70 s |

Instrumentation of the degree-16 case measured roughly 70 ms for candidate generation and about 1.5 s for the initial all-root isolation including all sixteen Rouché certificates. The pathological cost was therefore the irreducibility/factorization fallback and repeated materialization, not an unusually difficult root geometry. `--algebraic-root-cliffs` now permanently includes a general degree-16 Solve case, a degree-6 multi-prime factor-degree proof case, and a reducible degree-8 batch-canonicalization case so that Solve materialization is monitored together with raw Root isolation.


# 36. High-degree Complex Root Rouché proof, candidate iteration, and degree-budget re-audit

After batch canonicalization removed repeated all-root isolation from general polynomial `solve`, `solve[x^32-x+1==0,x]` still took about 25.5 s. Phase measurements showed about 0.40 s in Durand-Kerner candidate generation and about 0.05 s in disk separation/ordering, while Rouché certification of all 32 roots consumed about 11.05 s. The old proof built exact Rational Taylor coefficients around every disk center, repeatedly normalizing large Rationals and creating a representation cliff at higher degrees.

The current fast proof first evaluates `p(c)` and `p'(c)` with directed BigFloat intervals and proves the sufficient condition

```text
|p(c)| + (r^2/2) max_{|z-c|<=r}|p''(z)| < |p'(c)| r
```

for a candidate disk of center `c` and radius `r`. The `p''` bound is obtained by outward Horner evaluation of the absolute-coefficient polynomial at an upper bound for `|c|+r`. Every comparison uses directed upper/lower bounds; approximate candidates are never accepted as roots by themselves. If this inexpensive bound does not close, the implementation falls back to the exact Rational Taylor proof, whose translation is now formed by an O(n^2) Horner shift.

This reduced internal Rouché time from about 11.05 s to 0.14 s at degree 32 and from about 20.40 s to 0.47 s at degree 64. Once certification became cheap, Durand-Kerner candidate generation became the dominant cost. The old loop continued for at least roughly `degree` iterations even after corrections were already tiny. Candidate generation now stops after the fourth iteration once the correction criterion is met; if certification fails, the normal higher-precision retry still regenerates the candidates. Since certification, not Durand-Kerner convergence, is the acceptance condition, this does not weaken exactness.

Root materialization also no longer attempts to construct generator fields above the algebraic-field arithmetic candidate-degree limit of 16. Such fields could not be used by subsequent arithmetic and caused repeated high-degree irreducibility work, particularly around degree 65.

Representative GCC Release / LTO-off measurements from `--algebraic-root-cliffs 1` on 2026-08-28 are:

| workload | wall time |
| --- | ---: |
| `x^16-2` all-root create | about 0.158 s |
| two-radius degree-16 create | about 0.223 s |
| `solve[x^16+x+1==0,x]` | about 0.087 s |
| `solve[x^32-x+1==0,x]` | about 0.361 s |
| `solve[x^64+x+1==0,x]` | about 1.562 s |
| `solve[x^65+x+1==0,x]` | about 1.521 s |

Additional probes showed no algorithmic discontinuity at the former 64-degree boundary, with roughly 2.9 s at degree 80 and 4.1 s at degree 96. `maximumSupportedAlgebraicDegree` and the default `EvaluationLimits::maxAlgebraicDegree` are therefore raised from 64 to 96. This is a measured resource-safety boundary, not a mathematical domain restriction.

Aberth-Ehrlich iteration remains a future candidate-generation experiment. Any comparison should reuse the same Newton-polygon multi-radius initialization and measure ordinary, clustered, and scale-separated roots. If adopted, Aberth would replace only approximate candidate generation; Rouché certification, disk separation, and deterministic ordering remain the exact acceptance contract.


# 35. Direct construction of repeated symbolic derivatives

Repeatedly applying first-derivative rules can cause substantial expression growth for `LambertW`, `polylog`, and `exp[q(x)]`, because piecewise forms, quotient/product rules, and common exponential factors are rebuilt at every order. Since 2026-08-29, direct-variable requests with `n<=64` first use exact family-specific recurrences: the DLMF 4.13.4_1--4.13.4_2 polynomial recurrence for Lambert W, the `theta=xD` / signed-Stirling identity for polylogarithms, and `P_(n+1)=P'_n+q'P_n` for `exp[q]` with quadratic `q`.

These are exact Expr-construction paths, not approximate shortcuts. Orders above 64 or unmatched forms fall back to the generic `D` implementation, leaving the public order budget of 4096 and derivative semantics unchanged. No long benchmark was run for this batch; only compact representative outputs plus targeted compilation/smoke checks were used.

# 36. Hermite reduction and algebraic-log fallback for rational integration

The exact rational integrator was already fast and compact when a denominator reduced to linear and irreducible quadratic factors, but it stopped at higher-degree factors such as `x^3+x+1` and at repeated powers of those factors. Running a full Complex Root decomposition too early is also undesirable: it would replace elementary results such as `x/(1+x^4) -> atan[x^2]/2` with a much larger Root/Log sum. The higher-degree algebraic path therefore runs only after the elementary search has failed.

Multiplicity extraction no longer depends on the factor engine. A Yun square-free decomposition over Q[x] derives `Q=product f_i^i` exactly from `gcd(Q,Q')`. Powers of each square-free `f_i` are lowered one step at a time by the Hermite identity

```text
A/f^k = D[B/f^(k-1)] + C/f^(k-1)
B = -A (f')^(-1)/(k-1) mod f
```

where `(f')^(-1) mod f` is computed exactly by the extended Euclidean algorithm in Q[x]. No numerical roots participate in this reduction. For the final square-free part `P/Q`, all complex roots `r` are identified by the existing certified `ComplexAlgebraicNumber` backend and

```text
P(x)/Q(x) = sum_r P(r)/(Q'(r)(x-r))
```

produces `sum_r P(r)/Q'(r) Log[x-r]`. Residues are materialized through persistent exact algebraic-number arithmetic; approximate roots are not used to justify the identity.

This fallback remains inside the existing degree-12 specialized rational work budget. The bound is a resource-safety policy, not a mathematical domain boundary: it limits all-root isolation, algebraic residue materialization, and output size. Representative GCC Release checks were about 0.25 s for `integrate[1/(x^3+x+1),x]`, 0.56 s for `1/(x^5+x+1)`, and 4.0 s for `1/(x^8+x+1)`. No long-running benchmark suite was executed for this change.

Rothstein-Trager / Lazard-Rioboo-Trager residue grouping remains a future optimization. The current kernel is already exact but explicitly enumerates roots; grouping equal residues and conjugate roots could reduce expression size and recover more compact real `Log/atan` forms without replacing the Hermite/square-free capability added here.
