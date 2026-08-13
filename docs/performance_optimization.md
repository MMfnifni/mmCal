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

# 15.5. Precision-aware `N` and certified FFT

## Old path — construct exact FFT first, approximate afterward

Previously `N` received already-evaluated arguments, so

```text
N[fft[data],16]
```

first built the complete exact Fourier expression and only then approximated its components. Each butterfly therefore paid the full generic `Expr` multiply/add/simplification cost even when the caller only wanted decimal output.

## Precision-aware evaluation — selected

`N` now holds its first argument, resolves the requested precision first, and keeps a precision context active while evaluating the child expression. FFT consumes that context and performs the transform directly on certified `ComplexInterval`/BigFloat endpoints. Ordinary exact `fft[...]` is unchanged, and no machine `double` backend is introduced.

Representative benchmark in the same GCC Release environment at 16 fractional digits:

```text
32 points   exact ~11.6 ms   certified ~2.7 ms
64 points   exact ~63.5 ms   certified ~6.2 ms
128 points  exact ~327.6 ms  certified ~13.1 ms
```

For non-power-of-two certified transforms, direct DFT and Bluestein were measured against each other. Around 65 points direct evaluation remains about 60 ms, while at 127 points Bluestein is about 199 ms versus about 217 ms direct. The current policy therefore keeps direct evaluation below 96 points and uses Bluestein above it. This is a measured implementation threshold, not a mathematical constant, and should be remeasured on MSVC.

Why selected:

- preserves the exact-first public API while accelerating explicit approximation;
- keeps arbitrary precision and outward rounding without introducing machine `double`;
- makes precision propagation reusable by other expensive builtins;
- removes the approximate-path O(N^2) cliff for larger non-power-of-two sizes.

Exact symbolic FFT expression growth is a separate problem and is intentionally not changed here.

---

# 15.6. Array / Matrix Stage 1–2

## Flat Array + exact Number backend — selected

v1.5.2 keeps matrices on the shared Array `shape + row-major flat storage` representation instead of introducing a separate nested matrix Value. `MatrixView` reads an Array without copying; algorithms that mutate their workspace use a flat `MatrixBuffer`.

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
- `rref` remains fraction-free through the forward phase and introduces Rational values only during backward normalization.
- `matrixRank` uses the number of Bareiss pivots without materializing a full RREF.
- `inverse` writes `B=D A` and fraction-free eliminates the augmented matrix `[B|D]`, whose reduced right half is `A^-1`.
- Exact complex matrices retain the previous `Number` Gaussian/Gauss-Jordan fallback until a dedicated exact complex integer-domain representation is justified.

Pivot selection prefers the nonzero candidate with the smallest bit length to limit intermediate BigInt growth. Every Bareiss division is checked with `BigInt::divmod`; a nonzero remainder is treated as an invariant failure rather than silently truncating.

Release / LTO-off measurements on 2026-08-13 using the same benchmark matrices:

| size | Gaussian `det` | Bareiss `det` | speedup | Gauss-Jordan `rref` | Bareiss `rref` | speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 0.256 ms | 0.047 ms | 5.4x | 0.447 ms | 0.058 ms | 7.7x |
| 12 | 1.241 ms | 0.109 ms | 11.4x | 2.058 ms | 0.146 ms | 14.1x |
| 16 | 3.570 ms | 0.385 ms | 9.3x | 5.798 ms | 0.375 ms | 15.5x |

`N[det[...],p]`, `N[inverse[...],p]`, and related operations do not build these exact Bareiss results first. They continue to dispatch directly to the precision-aware certified Matrix backend shared with the FFT approximation infrastructure.

# 15.8. LU / Householder QR — selected

Decomposition code is grouped in `linear_algebra/decomposition.*`. `luDecomposition[A]` provides row-pivoted `P A = L U`, while `qrDecomposition[A]` uses Householder reflectors for `A = Q R`. Under `N[...]`, factors are not built exactly first; the active `ApproximationContext` dispatches directly to a certified `ComplexInterval` backend shared with the precision-aware Matrix/FFT infrastructure.

A column-block Householder application kernel was also tested, processing several columns in one row-major sweep. Repeated Release measurements with block sizes 1/8/16/32 on orders 8 through 24 produced only a few percent difference with no stable winning block; BigFloat/interval arithmetic dominated cache effects. The default therefore remains block=1-equivalent, while the block kernel and benchmark are retained for future backend changes.

Representative Release / LTO-off measurements from 2026-08-13:

| size | exact `LU` | `N[LU,16]` | `N[QR,16]` |
|---:|---:|---:|---:|
| 8 | 0.290 ms | 1.811 ms | 7.062 ms |
| 12 | 1.315 ms | 5.355 ms | 21.218 ms |
| 16 | 3.548 ms | 10.562 ms | 46.718 ms |

General exact Householder QR grows radical expressions rapidly: the same benchmark family measured about 1.2 ms at 2x2 and 59 ms at 3x3, while a 4x4 case took about 18 seconds and formatted to roughly 677 KB. General exact QR is therefore policy-limited to 3x3; only the upper-triangular `{I,A}` fast path remains exact at arbitrary order. General order 4+ use should go through `N[qrDecomposition[A],p]`.

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

Even a pure cubic extrapolation from order 64 suggests roughly 1.6 h for `N[LU]`, 1.4 h for `N[dot]`, 3.4 h for `N[det]`, 6.2 h for `N[solveLinear]`, 9.9 h for `N[inverse]`, 13 h for `N[SVD]`, 15 h for `N[QR]`, and 22 h for `N[eigenvalues]`. Extrapolating the observed 32→64 exponent instead gives a broad roughly 1–21 h range depending on the operation. These are projections, not 1024 completion measurements, and cache/allocation/guard-precision/iteration effects can make them worse.

A dense order of 1024 is not intrinsically huge in a machine-double + BLAS setting, but it is currently a stress regime for mmCal's exact Decimal→Rational→Expr representation and certified arbitrary-precision dense algorithms. If order-1024 dense work becomes an explicit target, compact numeric Array storage, fewer parser/lowering Expr allocations, and packed approximate-Matrix storage should be considered before more sophisticated blocking.

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

---

# 18. Next candidates

1. Replace the large fixed-size `Expr::Node` `std::variant` payload with kind-specific typed nodes and remeasure order-1024 dense Matrix RSS.
2. Benchmark BigUInt / BigInt small-object optimization independently to determine whether removing heap allocation for small integers pays off.
3. Investigate compact packed storage for numeric Arrays / approximate Matrices and reduce temporary allocation while parsing/lowering huge braces.
4. Remeasure blocked LU / QR only after the storage work changes the cost balance.
5. Measure Toom-4 / higher-Toom crossovers and consider FFT/NTT multiplication for still larger integers.
6. Lehmer GCD.
7. bit-burst / AGM logarithm backends.
8. 5,000–10,000-digit benchmark coverage for Gamma, erf, and other special functions.
9. A Cyclotomic exact FFT backend.

Items 1–3 specifically target the representation bottleneck exposed by the v1.5.2 order-1024 audit. The `Expr::Node` change should be a standalone refactor with the public `Expr` API and full regression corpus frozen, rather than being mixed with unrelated performance work.

Future changes should continue to record both adoption and rejection rationale in this document.
