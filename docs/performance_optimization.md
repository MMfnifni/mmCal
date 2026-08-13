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
- exact/certified FFT benchmarks and direct/Bluestein crossover measurements;
- fixed-seed certified FFT round-trip invariants.

Modes:

```text
mmCal.Benchmarks
mmCal.Benchmarks --full
mmCal.Benchmarks --random-only
mmCal.Benchmarks --benchmark-only
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

1. Toom-4 / higher-Toom crossover measurement;
2. FFT/NTT integer multiplication for still larger operands;
3. Lehmer GCD;
4. bit-burst / AGM logarithm backends;
5. 5,000–10,000-digit benchmark coverage for Gamma, erf, and other special functions;
6. allocator/SBO work only if benchmarks show a real gain.

Future changes should continue to record both adoption and rejection rationale in this document.
