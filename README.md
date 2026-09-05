# mmCalculator – Mathematical Machinery Calculator

An exact-first CLI calculator and compact CAS for engineering, research, and manufacturing.

© 2021–2026 mmKreutzef (aka Daiki.NIIMI)  
Licensed under the BSD 3-Clause License

**Latest release: v1.5.5 — Superior BugFix**

[English](README.md) | [日本語](README.ja.md)

## Overview

mmCalculator (mmCal) is an exact-first CLI mathematical calculator / compact CAS designed for technical work such as research, engineering, and manufacturing.

While it can be used as a general-purpose calculator, it also provides:

- **Computation with integers, fractions, and symbolic expressions kept in exact form whenever possible**
- Continuous calculations using previous results
- Variables and user-defined functions
- Complex numbers, vectors, and matrices
- Expansion, factorization, simplification, Taylor/Laurent/Puiseux/logarithmic series expansion, differentiation, integration, limits, and equation solving
- Arbitrary-precision numerical approximation, special functions, statistics, and signal processing

Unlike heavyweight systems, this tool aims to be:

> **Lightweight and ready to use — a practical tool for real work**

It also does not immediately convert input to `double` as many ordinary calculators do. Numerical approximation is performed only when needed.

## Quick examples — details later

```text
In [1]> 999999999999999999999999999999^2
Out[1]> 999999999999999999999999999998000000000000000000000000000001

In [2]> 0.1+0.2==0.3
Out[2]> True

In [3]> sqrt[72]+sin[Pi/6]
Out[3]> 1/2+6sqrt[2]

In [4]> factor[expand[(x+1)^3]-1]
Out[4]> x*(x^2+3x+3)

In [5]> fullSimplify[(x^2-1)/(x-1),x!=1]
Out[5]> x+1

In [6]> simplify[sqrt[x^2],element[x,Real]]
Out[6]> abs[x]

In [7]> D[exp[x^2],x]
Out[7]> 2x exp[x^2]

In [8]> integrate[x^2+sin[x],{x,0,Pi}]
Out[8]> 2+Pi^3/3

In [9]> limit[(1-cos[x])/x^2,x,0]
Out[9]> 1/2

In [10]> solve[1.1^x==x^2,x,Real]
Out[10]> {x == -2lambertw[log[11/10]/2]/log[11/10], x == -2lambertw[-log[11/10]/2]/log[11/10], x == -2lambertw[-1, -log[11/10]/2]/log[11/10]}

In [11]> N[%,20]
Out[11]> {x == -0.95548727594562198165, x == 1.0513800237472769374, x == 95.716830168405222740}

In [12]> inverse[{{1,2},{3,4}}]
Out[12]> {{-2, 1}, {3/2, -1/2}}

In [13]> ifft[fft[{1,2,3,4,5,6,7}]]
Out[13]> {1, 2, 3, 4, 5, 6, 7}
```

The following sections provide an overview only.
For function specifications and implementation details, see the [reference](docs/reference.md) or the Markdown documents in the `docs` directory.

Finite-precision results retain provenance and information metadata rather than only display text. `ExactValue` and `CertifiedInterval` approximations carry rigorous enclosures and can participate in certified arithmetic, while residual-verified SVD / eigen candidates are explicitly tagged `VerifiedApproximation` instead of being promoted to mathematical point enclosures.

## v1.5.5

v1.5.5 adds Series and Vector Calculus while focusing primarily on **cross-cutting bug fixes, composition between existing symbolic frontends, finite-precision semantics, performance cliffs, and removal of obsolete fixed limits**.

Major changes include:

- **Series and asymptotics**: added `SeriesData` / TPSA-based Taylor, Laurent, Puiseux, and logarithmic series, `+Infinity` expansions, and `toNormal`, connected directly to `D`, `integrate`, and `limit`
- **Array and Vector**: regularized the Vector API around Hermitian inner products and added Cartesian Vector Calculus including `grad`, `divergence`, `curl`, `laplacian`, `jacobian`, and `hessian`
- **Symbolic composition**: fixed Solver definedness, principal inverses, `cases` boundaries, Limit binder capture, nested `D` / `integrate` / `series` / `solve`, and transform-front-end combinations that previously produced unevaluated or semantically damaged expressions
- **Certified `N` and formatting**: fixed finite-precision signs, parentheses, and polynomial ordering; preserved structural exact integers; and split approximation provenance into `ExactValue`, `CertifiedInterval`, and `VerifiedApproximation` so rigorous enclosures are not conflated with residual-verified numerical candidates
- **Integration, differentiation, and algebra**: extended higher-degree Rational-function integration through Yun decomposition, polynomial CRT, Hermite reduction, and exact residues, while strengthening dedicated recurrences for higher derivatives and polynomial-times-elementary integrals
- **Performance and limits**: removed perfect-power false-positive cliffs such as `factor[x^257-1]`, and replaced obsolete fixed 64 / 128 / 256 / 4096 boundaries with closed forms, structural fast paths, or shared `EvaluationBudget` limits where safe
- **CLI and validation**: added `:quit` / `:exit` / `:layout` and expanded certification-boundary, performance-cliff, and cross-feature regression coverage

See **v1.5.5** in [`CHANGELOG.md`](CHANGELOG.md) for the detailed release history. README and Reference describe the v1.5.5 behavior; release-by-release history remains in the changelog.

## 1. Getting started

On Windows, simply launch `mmCal.exe`.

The display precision and default angle unit can also be specified at startup.

```text
mmCal --fix 16 --angle deg
mmCal --angle rad
mmCal --angle grad --fix 8
mmCal --layout multi
mmCal --eval "factor[x^2-1]"
mmCal --batch < expressions.txt
```

- `--fix 16`: Display results with **up to 16(0..1000) digits after the decimal point**. Unnecessary trailing zeros are omitted
- `--angle deg`: Treat trigonometric inputs without an explicit angle unit as degrees
- `--angle rad`: Radians. This is the default
- `--angle grad`: Gradians
- `--layout auto|single|multi`: Interactive REPL layout. The default `auto` structurally expands Arrays/Lists/`cases`/solution sets on a TTY and falls back to one-line output through pipes or redirects
- `--eval expr`: Evaluate one expression and write only its value to standard output
- `--batch`: Evaluate standard input one line at a time in one session.
- `--help`, `-h`: Show concise startup usage. Use commands such as `:help sin` after startup for function details

`--eval` and `--batch` are automation modes: they do not emit the banner, prompts, `Out[...]` labels, or farewell. Values go to standard output; warnings and errors go to standard error. Stable exit codes are `0` for success, `2` for arguments, `3` for syntax/resource limits, `4` for evaluation errors, and `5` for internal errors. Batch mode continues after an error and returns the greatest exit code observed.

On Linux and similar systems, mmCal can be built from source using CMake 3.20 or later with GCC or Clang.

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

For CMake-generated MSVC builds, `MMCAL_PARALLEL_COMPILE=ON` is the default and adds `/MP`. Disable it at configure time with `-DMMCAL_PARALLEL_COMPILE=OFF` when necessary. GCC and Clang do not receive a compiler-specific parallel-build flag; parallelism is delegated to Ninja, Make, or the selected build tool through `cmake --build ... --parallel`. Unity builds are intentionally not enabled by default because the existing large translation units would increase memory pressure and make incremental rebuilds coarser.

## 2. Exact values and decimal display are different things

mmCal is exact-first. Integers, rationals, finite decimal literals, radicals, complex values, and symbolic expressions remain exact whenever practical.

```text
1/3      -> 1/3
0.125    -> 1/8
sqrt[2]  -> sqrt[2]
```

Use `N[expr,p]` when a decimal approximation is wanted. `p` means significant decimal digits, not a fixed number of digits after the decimal point.

```text
N[1/3,20]
-> 0.33333333333333333333

N[Pi,30]
-> 3.14159265358979323846264338328
```

Exact integer components that remain mathematically exact after numerical evaluation are displayed without a redundant `.0`. `:fix` changes presentation only and does not replace the stored value.

Numerical values distinguish `ExactValue`, `CertifiedInterval`, and `VerifiedApproximation`. Use `precision`, `accuracy`, `rationalize`, and `explain` to inspect them. Detailed enclosure, provenance, stopping, and branch semantics are documented in the Reference.

## 3. Basic syntax

Function calls use square brackets `[]`; parentheses `()` are reserved for grouping.

```text
sin[Pi/6]
log[10,1000]
(x+1)^2
2Pi
```

Arrays and general brace containers use `{...}`.

```text
{1,2,3}
{{1,2},{3,4}}
```

See the Cheatsheet / Reference for assignment, comparisons, conditions, implicit multiplication, literals, and operator precedence.

## 4. Angles

Radians are the default. `angleMode` changes the session default, while `Deg` / `Rad` / `Grad` can be written explicitly in an expression.

```text
sin[Pi/6]
angleMode[Deg]
sin[30]
sin[30 Deg]
```

Angle-conversion functions: `DtoR`, `DtoG`, `RtoD`, `RtoG`, `GtoD`, `GtoR`

## 5. Variables and user-defined functions

```text
x:=2
f[t]:=t^2+1
f[4]
-> 17
```

Definition/session operations: `Defs`, `UnDef`, `Clear`, `Exit`

See the Reference for evaluation and scoping rules.

## 6. Input/output history

`%` refers to the most recent successful output and `@` to the most recent input. Repeated forms such as `%%` / `@@` move further back.

The formal history functions are `In[n]` / `Out[n]`; positive indices are absolute and negative indices are relative.

```text
Out[-1]
In[-1]
```

`Out` returns a stored output snapshot, while `In` re-evaluates the stored input in the current definition environment. See the Reference for details.

## 7. Main mathematical features

The names below are the main canonical function names. Use `:help functions` for the complete current list including callable aliases, the Cheatsheet for examples, and the Reference for exact forms and domains.

### Numeric and complex

`N`, `precision`, `accuracy`, `rationalize`, `explain`, `sqrt`, `cbrt`, `abs`, `sign`, `re`, `im`, `conj`, `arg`, `cis`, `polar`, `proj`, `hypot`, `fma`, `clamp`

### Elementary functions

`exp`, `expm1`, `expc`, `log`, `log1p`, `log2`, `log10`

`sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `asin`, `acos`, `atan`, `atan2`

`sinh`, `cosh`, `tanh`, `csch`, `sech`, `coth`, `asinh`, `acosh`, `atanh`

`sinc`, `cosc`, `tanc`, `sinhc`, `tanhc`

### Integer and discrete mathematics

`floor`, `ceil`, `trunc`, `round`, `frac`, `gcd`, `lcm`, `mod`, `rem`, `quotient`

`bitand`, `bitor`, `bitxor`, `bitnot`, `bitshiftl`, `bitshiftr`, `bitlength`, `bitcount`, `bitget`

`isprime`, `nextprime`, `prevprime`, `factorint`, `totient`, `perm`, `comb`, `fib`, `nextpow2`

### Symbolic algebra and calculus

`simplify`, `fullSimplify`, `expand`, `factor`, `collect`, `cases`

`series`, `normal`, `toNormal`

`D`, `diff`, `integrate`, `nintegrate`, `limit`

`solve`, `root`, `groebnerBasis`, `polynomialReduce`, `element`, `if`

### Special functions

`gamma`, `lgamma`, `beta`, `betaln`, `ibeta`, `binom`, `fallingfact`, `risingfact`

`erf`, `erfc`, `zeta`, `digamma`, `trigamma`, `lambertw`

`fresnelc`, `fresnels`, `hypergeometric1F1`, `hypergeometric2F1`

`ellipticF`, `ellipticE`, `ellipticPi`

`Ei`, `Si`, `Ci`, `li`, `polylog`

When a symbolic result cannot be established safely, mmCal may return a Warning and leave the expression unevaluated rather than inventing a result.

## 8. Arrays, matrices, vectors, and statistics

### Arrays and sequence operations

`dimensions`, `arrayRank`, `length`, `at`, `reshape`, `zeros`, `identity`

`range`, `table`, `map`, `sum`, `prod`, `min`, `max`

Indices are zero-based.

### Matrices and linear algebra

`transpose`, `conjugateTranspose`, `madd`, `dot`, `det`, `inverse`, `rref`, `matrixRank`

`solveLinear`, `nullSpace`, `luDecomposition`, `qrDecomposition`, `svd`

`conditionNumber`, `leastSquares`, `pseudoInverse`

`eigenvalues`, `eigenvectors`, `eigensystem`

`trace`, `rows`, `cols`, `diag`

### Vectors

`cross`, `norm`, `manhattanDistance`, `distance`, `normalize`

`projection`, `rejection`, `vectorAngle`, `reflectNormal`, `reflectAxis`

`inner`, `outer`, `orthogonalQ`, `orthonormalQ`, `linearIndependentQ`, `gramSchmidt`

### Vector calculus

`grad`, `divergence`, `curl`, `laplacian`, `jacobian`, `hessian`, `directionalDerivative`

### Statistics

`mean`, `median`, `mode`, `quantile`, `percentile`

`var`, `vars`, `stddev`, `stddevs`, `geomean`, `harmmean`, `rms`

`mad`, `madR`, `skew`, `kurtp`, `kurts`, `cv`, `stderr`, `zscore`, `iqr`

`trimmean`, `winsor`, `winsorR`, `cov`, `corr`, `corrspearman`, `percentrank`

See the Reference for exact / finite-precision behavior and shape requirements.

## 9. FFT and random numbers

Signal processing: `dft`, `fft`, `ifft`, `convolve`

Random numbers: `randSeed`, `rand`, `randint`, `choice`, `randn`

Random generators keep session state and reproduce the same sequence when reseeded with the same value. They are not cryptographic generators.

## 10. About Warnings

`D`, `integrate`, `limit`, `solve`, `N`, and other operations may return a Warning and leave an expression unevaluated when the current implementation cannot establish a safe result.

A Warning does not mean that the unevaluated expression is the mathematical answer. It reports conditions such as unresolved evaluation, insufficient assumptions or precision, or an unsupported backend. See the Reference for the detailed classification.

## 11. CLI help and display settings

```text
:help
:help sin
:help functions
:help constants
:fix 16
:fix off
:layout auto
:layout single
:layout multi
:status
:quit
:exit
```

`:help functions` lists the current canonical function names and callable aliases. `:help <function>` shows accepted forms and examples.

`:fix` changes decimal presentation only, `:layout` changes interactive REPL composition only, and `:status` reports current session state.

## 12. Naming

Ordinary mathematical functions use lowercase canonical names.

```text
sin cos log sqrt integrate solve
```

Short symbolic and session operations use proper names.

```text
D N In Out Exit Clear Defs UnDef
```

Use `:help functions` for the current complete name list including aliases.

## 13. Detailed documentation

- `docs/reference.md` — Detailed function, syntax, and current-behavior reference
- `docs/mathematics.md` — Mathematical policy for domains, principal values, and numerical evaluation
- `docs/architecture.md` — Internal architecture for developers
- `docs/grammar.ebnf` — Machine-readable overview of the grammar
- `docs/performance_optimization.md` — Performance work adopted or rejected, with benchmark rationale
- `docs/evaluation_budget.ja.md` — Request-scoped evaluation limits, cancellation, telemetry, and diagnostic contract (Japanese)
- `docs/multiprecision_implementation.ja.md` — Detailed Japanese notes on the multiprecision / certified numerical backend
- `CHANGELOG.md` — Major changes by release

## 14. License and trademarks

The source code is distributed under the **BSD 3-Clause License**. See `LICENSE` for the copyright permissions governing commercial use, modification, redistribution, and incorporation into other software.

**Use of the mmCal name, official logos, and related branding is addressed separately by `TRADEMARKS.md` .** The trademark policy does not prohibit independent GUIs, forks, commercial products, or other uses permitted by the BSD license; it is intended to avoid confusion about whether a third-party product is mmCal itself or is officially endorsed.

- [BSD 3-Clause License](LICENSE)
- [Trademark and Brand Policy](TRADEMARKS.md)

If mmCal is used in an academic publication or product, attribution beyond the BSD requirements is not mandatory, but a factual acknowledgement is appreciated.

## 15. Tests and development environment

v1.5.5 has been verified with **3428 / 3428** internal regression tests and **2465 / 2465** black-box tests passing.
They focus especially on exact arithmetic, boundary values, domains, error classification, formatter round-trip parsing, and certified numerical enclosures. A separate `mmCal.Benchmarks` project provides fixed-seed randomized correctness checks, algorithm-threshold sweeps, and performance comparisons without mixing benchmark workloads into the ordinary test suite.

Primary Windows development environment:

- Windows 10 1909
- Intel Core i7-9800X
- Visual Studio Community 2026 / MSVC

Builds and tests are also performed on Linux with GCC and Clang.

LLMs are used as an auxiliary tool for documentation organization, implementation-policy review, test design, and organizing MathKnowledge such as integration rules.

## 16. Notes

This project aims to combine rigorous operator semantics with practical expression evaluation.
It also aims to remain easy to run as a lightweight CLI in research, design, and manufacturing environments.

And, ultimately, I am building the tool I want to use.

## 17. Disclaimer

This software is provided **AS IS**, as stated in the BSD 3-Clause License.
No express or implied warranties are made, including warranties of merchantability, fitness for a particular purpose, or non-infringement.
The author or copyright holders shall not be liable for any claim, damages, or other liability arising from or related to the software, whether in contract, tort, or otherwise.

By using this software, you acknowledge that all risks associated with its use are your own responsibility. The author shall not be liable for data loss, system malfunction, or other damage resulting from use of the software.

See `LICENSE` for the authoritative terms and disclaimer.

## Acknowledgements

I would like to thank my past self for starting this tool, the university and professors who provided a place to learn, and my current workplace for providing a real-world environment in which to use it.

In addition, the numerical calculations were implemented with reference to [DLMF](https://dlmf.nist.gov/). We would like to express our gratitude to the authors.

## Requests and feedback

Please submit requests, bug reports, and implementation proposals on GitHub. They are very welcome 🍀

## Future ideas

- It would be nice to have something like `for`
- A `plot` function for graphing
