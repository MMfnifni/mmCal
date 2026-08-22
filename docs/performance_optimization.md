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

# 15.4. Exact Cyclotomic FFT — Stage 7-7

## Previous generic `cis` / Expr direct DFT — retained as fallback

Before Stage 7-7, exact non-power-of-two FFTs fell from `radix2Transform()` into `directTransform()`, constructing every twiddle as a generic `cis[-2 Pi k/n Rad]` expression. Even at lengths 5/7/10/12, forward transforms accumulated root-of-unity expressions and inverse transforms asked the generic Simplifier to rediscover the same cyclotomic identities, making exact round trips a major part of `test5_matrix.txt`.

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

`tester.py --timings` reduced `test5_matrix.txt` from about 1.46 s at the Stage-5-3 audit to about 0.38 s after Stage 7-7.

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

The same audit separates mathematical convergence from practical bounded work for certified special functions. On this GCC Release build, `polylog[2,0.999]`, `2F1[...,0.98]`, and `ellipticF[...,0.98]` are mathematically inside their series convergence regions but become multi-second or worse with the current exact-majorant implementations. Conservative work boundaries are therefore `|z|<=49/50` for positive-order polylog, `<=9/10` for the 2F1/elliptic series, `|z|<=160` for 1F1, and 96 for the real Ei/Si/Ci series. Failures caused by a fixed backend range, term cap, or planner cap are classified as `CertifiedBackendUnsupported`, preventing meaningless guard-precision retries.

Why selected:

- preserves the exact-first public API while accelerating explicit approximation;
- keeps arbitrary precision and outward rounding without introducing machine `double`;
- makes precision propagation reusable by other expensive builtins;
- removes the approximate-path O(N^2) cliff for larger non-power-of-two sizes.

Exact symbolic FFT expression growth is a separate problem and is intentionally not changed here.

---

# 15.6. Array / Matrix Stage 1–2

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

This stage intentionally retained ordinary Gaussian/Gauss-Jordan exact elimination; integer/Rational Bareiss elimination was kept as a separate Stage 3 optimization.

# 15.7. Bareiss / fraction-free exact Matrix — selected

Stage 3 lifts exact real matrices to integers by clearing denominators independently per row, then runs a shared Bareiss kernel on `IntegerMatrixBuffer`. Integer inputs require no lift; Rational inputs are scaled only for the elimination workspace.

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

The second representation stage packs persistent `ArrayExpr` values into fixed immutable pages and separates shape/offset/strides from the backing. A first naive single-`vector<Rational>` design reduced storage overhead but made transpose deep-copy one million Rational/BigInt values; direct-packed order-1024 transpose took roughly 650–675 ms, so that design was rejected. The adopted design shares immutable 1024-element pages and makes transpose a stride-only view operation.

`ArrayBuilder` promotes only the current page. Completed pages are immutable, so a symbolic value near the end of a huge numeric Array does not force the preceding data through a Generic-Expr rebuild. A dedicated 1,048,576-element test with only the final value changed to `x` was essentially identical to the all-integer case at about 0.30 s / 69.8 MiB, with only the final page becoming Generic. Rectangular brace lowering also streams numeric leaves directly into one builder instead of first allocating an Expr node for every scalar.

With the benchmark fixture likewise constructing exact Rationals directly through `ArrayBuilder`, `--matrix-large transpose 1024 16` now measured `136576 KiB` (~133.4 MiB) maximum RSS and about 0.059 ms for the transpose view itself. `trace 1024` measured about 133.5 MiB and 10.95 ms for trace. Because the fixture construction path is intentionally part of the new representation work, the 292.6→133.4 MiB change is an end-to-end storage+builder improvement rather than the same kind of single-change A/B used for typed nodes.

A 13.63 MB 1024x1024 ten-decimal literal fed through the CLI and evaluated only as `dimensions[...]` measured about 4.55 s wall time and `441564 KiB` (~431 MiB) maximum RSS. The older 14.16 MB / 10.9 s / 1.99 GB measurement did not use byte-identical input, so this is not presented as a strict A/B comparison, but it confirms that avoiding per-scalar Expr allocation in lowering materially reduces the post-parse representation cost.

Even a pure cubic extrapolation from order 64 suggests roughly 1.6 h for `N[LU]`, 1.4 h for `N[dot]`, 3.4 h for `N[det]`, 6.2 h for `N[solveLinear]`, 9.9 h for `N[inverse]`, 13 h for `N[SVD]`, 15 h for `N[QR]`, and 22 h for `N[eigenvalues]`. Extrapolating the observed 32→64 exponent instead gives a broad roughly 1–21 h range depending on the operation. These are projections, not 1024 completion measurements, and cache/allocation/guard-precision/iteration effects can make them worse.

A dense order of 1024 is not intrinsically huge in a machine-double + BLAS setting, but it remains a stress regime for mmCal's certified arbitrary-precision dense algorithms. Persistent exact-Array representation cost is now substantially lower after typed nodes, paged packed backing, and direct builder lowering. Approximate SVD/Eigen paths already use dedicated contiguous working buffers, so blocking and selective threading should now be remeasured against the new storage balance instead of introducing another persistent-matrix type first.


# 15.12. Persistent algebraic fields, Stages 5-1 through 5-3

## Compositum / embedding reuse — adopted

Stage 5-1 retains a proven primitive-element compositum and both operand power-basis embeddings for a Root pair in a bounded cache. Real `AlgebraicElement` materialization also evaluates the coordinate polynomial over the chosen generator interval and uses an exact Sturm root count to identify the canonical `root[minpoly,k]` directly. This avoids duplicate all-root isolation and impossible Rational / `Q+iQ` degeneration probes.

For the representative expression

```text
(root[{-2,0,1},2]+root[{-3,0,0,1},1])
*(root[{-2,0,1},2]-root[{-3,0,0,1},1])
```

Stage 4 took about 1.23 s per evaluation. Stage 5-1 reduced this to roughly 0.14–0.15 s on the first evaluation and about 0.02 s warm in the same GCC Release / LTO-off environment.

## Reciprocal reuse — adopted

Stage 5-2 memoizes exact reciprocals computed by extended Euclid in `Q[t]/(m)` in a thread-safe per-`NumberFieldContext` LRU capped at 16 reciprocal pairs. One entry is bidirectional because `inverse(inverse(x)) = x`. Rational constant coordinates invert directly through the canonical embedding of `Q`.

Representative measurements:

| degree | operation | before cache | Stage 5-2 warm |
|---:|---|---:|---:|
| 6 | reciprocal | ~178 us | ~0.24 us |
| 6 | divide | ~243 us | ~84 us |
| 12 | reciprocal | ~525 us | ~0.44 us |
| 12 | divide | ~773 us | ~211–233 us |

### Persistent multiplication-matrix cache — rejected

For a degree-12 field, a representative ordinary multiplication took about 112 us, matrix-vector multiplication about 101 us, while constructing the left-multiplication matrix cost about 228 us. The roughly 10% per-multiply saving is not enough to amortize construction unless the same multiplier is reused many times. The memory and complexity cost of attaching such a cache to general `AlgebraicElement` values is therefore not justified at present.

## Incremental Krylov minimal polynomial — adopted

The old `AlgebraicElement::minimalPolynomial()` rebuilt a Rational matrix and reran Gauss-Jordan from scratch for every candidate degree in `1,a,...,a^k`, repeatedly discarding the same independence information.

Stage 5-3 appends the Krylov sequence `1,a,a^2,...` one column at a time and reduces only the new column against a persistent exact row-echelon state. The first dependence

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

Representative Stage 5-3 / GCC Release / LTO-off results for all 1652 black-box cases on 2026-08-16:

| test file | tests | wall time |
|---|---:|---:|
| `test16_exact_calculus_solver.txt` | 85 | ~2927 ms |
| `test5_matrix.txt` | 105 | ~1457 ms |
| `test9_special_func.txt` | 75 | ~297 ms |
| `test8_calculus.txt` | 46 | ~173 ms |
| `test22_number_field_interning.txt` | 1 | ~152 ms |

`test16` is mainly a stress set for integration, high-degree Solve, and algebraic Root construction. Despite its name, much of the time in `test5_matrix` is in exact FFT/DFT non-power-of-two round trips rather than small matrix operations; exact `ifft[fft[...]]` around 7, 12, and related sizes remains a future performance target. In `test9`, `N[ibeta[1/3,2/3,1/4],20]` and `N[gamma[1/3],20]` stand out comparatively.

These timings never affect PASS/FAIL semantics; they are profiling signals used only to select optimization targets.

# 15.14. Certified `gamma` / `ibeta` — Step 6-1 / 6-2

The `tester.py --timings` audit identified `N[ibeta[1/3,2/3,1/4],20]` and `N[gamma[1/3],20]` as clear hotspots inside `test9_special_func.txt`, so both certified backends were measured directly.

## `ibeta` point/shared normalization — adopted

The old interval wrapper evaluated both endpoints even when `lower == upper`, and each endpoint rebuilt `Beta(a,b)`. Since the 2F1 series itself was only about 1.9 ms at the representative 20-digit workload while Beta/Gamma normalization dominated, Step 6 now:

- evaluates an exact point once;
- computes `Beta(a,b)` once and shares it across interval endpoints;
- reuses the same normalization through the complement identity because `B(a,b)=B(b,a)`;
- uses only +40 guard bits for non-complement point evaluation and +80 when complement/interval evaluation requires it;
- returns exact `x=0,1` before constructing the normalization.

Representative warm/direct timings moved from roughly 116→37–43 ms at 80 bits, 255→98–117 ms at 160 bits, 835→387–392 ms at 320 bits, and about 10.3→5.3 s for the first 640-bit call. Repeated 640-bit calls fall to about 1.0 s once the Gamma plan cache is warm.

## Gamma lazy Bernoulli / Horner / plan reuse — adopted

The old Gamma backend eagerly generated `B0...B128` on first use. Step 6 preserves the Akiyama–Tanigawa state and extends it only to the requested even Bernoulli order under a mutex. It also limits low/mid-precision plan search to `min(64,max(16,ceil(bits/5)))`, evaluates the Stirling polynomial in Horner form, builds exact-point recurrence products as balanced Rational products, adds exact `Gamma(1)=Gamma(2)=1` / `logGamma(1)=logGamma(2)=0` fast paths, and keeps a thread-local bounded cache of exact Stirling plans.

Separate-process cold timings for `gamma[1/3]` improved from about 79→9.5 ms at 80 bits and 89→43 ms at 160 bits. The 320/640-bit first-call cost is approximately unchanged; the remaining high-precision cost is in the actual Stirling/recurrence work. Repeated 640-bit `lgamma[1/3]` falls from roughly 1.36 s to about 0.34 s when the plan cache is warm.

### Eager Bernoulli generation through `B256` — measured and rejected

Simply extending the old Akiyama–Tanigawa initialization through `B256` was also evaluated as a prerequisite for longer Stirling sums. Exact Rational generation alone rose from roughly 68 ms through `B128` to roughly 560 ms through `B256`, imposing that fixed first-use cost even on low-precision calls. The stateful lazy cache is therefore retained instead.

### `maximumK > 64` to reduce shift — measured and rejected

Allowing `K=96` does reduce the large recurrence shift at roughly 640 bits, but the larger exact Bernoulli coefficients, Rational/interval conversion and longer Stirling sum outweighed that gain: the representative workload regressed from about 1.37 s to about 2.7 s. Keep `K<=64` until a rectangular/binary-splitting Stirling kernel and a cheaper high-order Bernoulli backend justify reevaluation.

### fixed-k exact binary-search planning — measured and rejected

A fixed `k=64` binary search over the exact remainder bound was also tested. Repeated construction of huge Rational powers `x^(2k-1)` made the planner itself slower than the existing scan, so it is not retained.

The replaced implementations remain commented next to the new code in `certified_special_functions.cpp`, together with the reason for replacement, specifically for this algorithm-comparison cycle.

Note: the stateful lazy Bernoulli generator described in this Step 6-2 section was subsequently replaced in Step 6-3 by a static exact `B_2...B_128` table; the generator remains only as commented comparison code.

# 15.15. Exact-Rational `Gamma` / high-precision Stirling planning — Step 6-3

Reprofiling the 640-bit and higher `gamma[1/3]` / `ibeta[1/3,2/3,1/4]` paths after Step 6-2 showed that major costs remained at representation boundaries and in parameter planning, before any need for a fundamentally different Gamma formula. The design was cross-checked against Fredrik Johansson, *Arbitrary-precision computation of the gamma function* (arXiv:2109.08392), especially its treatment of rational rising factorials and Stirling parameter selection.

## Preserve exact Rational identity through the certified backend — adopted

The evaluator knew `1/3` exactly, but the old special-function path first converted it with `RealInterval::fromRational`. A non-dyadic Rational is not a point interval, so the exact rising-factorial path added in Step 6-2 was effectively bypassed for the representative `gamma[1/3]` workload. Step 6-3 now:

- detects exact Rational arguments before interval conversion for `Gamma`, `LogGamma`, `Beta`, and `BetaLog`;
- preserves the original `p/q` through positive-Rational Gamma argument shifting;
- forms `(p/q)_n` as a balanced binary product `prod(p+qk)/q^n` and canonicalizes one final Rational;
- passes exact `a`, `b`, and `a+b` into the three LogGamma evaluations used by Beta;
- preserves exact Rational identity under negative-argument reflection and evaluates `sin(Pi x)` as exact turns `sinTurns(x/2)`.

## Static exact table for `B_2...B_128` — adopted

The stateful lazy Akiyama–Tanigawa generator from Step 6-2 avoided low-precision eager initialization, but the first call that reached the highest Bernoulli orders still paid a large exact-Rational state-update cost. The current Stirling kernel only needs the fixed constants `B_2...B_128`, so they are stored as exact decimal numerator/denominator literals and parsed into `BigInt/Rational` only on first reference. The replaced lazy generator remains commented beside the new implementation with its replacement rationale. This does **not** revive the rejected idea of runtime eager generation through `B256`; higher dynamic Bernoulli generation remains deferred.

## Certified BigInt high-precision planner — adopted

The old planner advanced the shift in steps of eight and tested `k=1...64` with exact Rational arithmetic at every candidate. At 1000-bit precision this planning work itself became a large cold-start cost. For positive Rational `x=(p+qs)/q`, coefficient `c=A/B`, and `d=2k-1`, the test

```text
|c| / x^d <= 2^-P
```

is exactly equivalent to

```text
|A| q^d 2^P <= B (p+q s)^d.
```

Above 768 bits, Step 6-3 tests this integer inequality for `k=64` and finds a sufficient shift with doubling plus binary search. The final remainder bound is rebuilt as an exact Rational, so no floating heuristic weakens the certified contract.

This is distinct from the fixed-k Rational-power binary search rejected in Step 6-2: that version repeatedly constructed normalized Rational powers `x^(2k-1)`. The adopted planner removes those GCD/normalization costs and compares BigInt cross products directly.

## Representative measurements

Representative GCC Release / LTO-off timings in the same environment:

| workload | Step 6-2 / analysis baseline | Step 6-3 |
|---|---:|---:|
| `gamma[1/3]`, 640 bit | about 0.69 s | about 65 ms first / 38 ms warm average |
| `gamma[1/3]`, 1280 bit first | about 1.5 s | about 0.11 s |
| `ibeta[1/3,2/3,1/4]`, 640 bit | about 5.3 s | about 0.21 s |
| `ibeta[1/3,2/3,1/4]`, 1280 bit | about 4.5 s | about 0.49 s |
| `gamma[-1/3]`, 1280 bit | generic interval reflection | about 0.15 s |

A normal `--special-functions 3` run gives roughly 4.5/7.2/13.7/37.6 ms for `gamma[1/3]` at 80/160/320/640 bits and about 15.8/32.2/44.3/210 ms for the corresponding `ibeta` workload. The permanent benchmark now also includes 1280 bits.

## Improved Stirling main sum / Algorithm 6 — next candidate

Johansson's Theorem 3.5 / Algorithm 6 splits the Stirling main sum into low-index Bernoulli terms and a re-expanded high-index hypergeometric tail, reducing the number of Bernoulli values needed at high precision. The FLINT/Arb Gamma backend likewise documents an improved Stirling sum using rectangular splitting for low-index terms and high-index re-expansion. Step 6-3 removes the dominant representation/planner overhead through the 1280-bit range; future work above this range should therefore implement this improved main-sum direction instead of simply increasing `K` or returning to runtime `B256` generation.

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
| runtime eager Bernoulli generation through `B256` | rejected | exact Rational generation alone costs about 560 ms. Step 6-3 uses a fixed exact `B_2...B_128` table, but does not dynamically generate a larger range |
| Gamma Stirling `maximumK>64` | rejected | reduced shift but high-order Bernoulli/Rational and longer Stirling work regressed the 640-bit case from ~1.37 s to ~2.7 s |
| Gamma fixed-k **Rational-power** binary search | rejected | normalized Rational `x^(2k-1)` probes are expensive; Step 6-3 adopts a different BigInt cross-product formulation of the same certified test |

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
