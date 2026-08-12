# Changelog

## v1.5.2 — in development

### Syntax (breaking change)

- Function-call syntax is now exclusively `name[...]`; parenthesized `name(...)` calls were removed
- Parentheses `()` are grouping-only, while `x(x+1)` for an ordinary identifier is parsed as implicit multiplication
- Legacy `sin(x)` syntax on a known function name raises SyntaxError instead of being silently treated as multiplication
- The Formatter always emits `name[...]` for calls and explicit `x*(...)` for identifier/group multiplication, preserving unambiguous round trips
- Removed function-call delimiter branching from Parser/AST/Lowerer


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
