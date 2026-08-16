# mmCalculator – Mathematical Machinery Calculator

An exact-first CLI calculator and compact CAS for engineering, research, and manufacturing.

© 2021–2026 mmKreutzef (aka Daiki.NIIMI)  
Licensed under the BSD 3-Clause License

**Current release: v1.5.2**

[English](README.md) | [日本語](README.ja.md)

## Overview

mmCalculator (mmCal) is an exact-first CLI mathematical calculator / compact CAS designed for technical work such as research, engineering, and manufacturing.

While it can be used as a general-purpose calculator, it also provides:

- **Computation with integers, fractions, and symbolic expressions kept in exact form whenever possible**
- Continuous calculations using previous results
- Variables and user-defined functions
- Complex numbers, vectors, and matrices
- Expansion, factorization, simplification, differentiation, integration, limits, and equation solving
- Arbitrary-precision numerical approximation, special functions, statistics, and signal processing

Unlike heavyweight systems, this tool aims to be:

> **Lightweight and ready to use — a practical tool for real work**

It also does not immediately convert input to `double` as many ordinary calculators do. Numerical approximation is performed only when needed.

## Quick examples — details later

```text
In [1]> 999999999999999999999999999999^2
Out[1]> 999999999999999999999999999998000000000000000000000000000001

In [2]> 0.1+0.2
Out[2]> 3/10

In [3]> 0.1+0.2==0.3
Out[3]> True

In [4]> 1/3+1/6
Out[4]> 1/2

In [5]> sqrt[72]
Out[5]> 6sqrt[2]

In [6]> sin[Pi/6]
Out[6]> 1/2

In [7]> expand[(x+1)^3]
Out[7]> x^3+3x^2+3x+1

In [8]> factor[x^2-1]
Out[8]> (x-1)(x+1)

In [9]> fullSimplify[(x^2-1)/(x-1),x!=1]
Out[9]> 1+x

In [10]> simplify[sqrt[x^2],element[x,Real]]
Out[10]> abs[x]

In [11]> D[exp[x^2],x]
Out[11]> 2x exp[x^2]

In [12]> integrate[x^2+sin[x],x]
Out[12]> x^3/3-cos[x]

In [13]> integrate[sin[x],{x,0,Pi}]
Out[13]> 2

In [14]> limit[(1-cos[x])/x^2,x,0]
Out[14]> 1/2

In [15]> solve[x^2+1==0,x,Complex]
Out[15]> {x==I,x==-I}

In [16]> (1+I)/(1-I)
Out[16]> I

In [17]> sqrt[-8]
Out[17]> 2I sqrt[2]

In [18]> a:=12
Out[18]> 12

In [19]> a^2+1
Out[19]> 145

In [20]> @+a
Out[20]> 157

In [21]> A:={{1,2},{3,4}}
Out[21]> {{1,2},{3,4}}

In [22]> det[A]
Out[22]> -2

In [23]> inverse[A]
Out[23]> {{-2,1},{3/2,-1/2}}

In [24]> solveLinear[A,{5,11}]
Out[24]> {1,2}

In [25]> dot[A,inverse[A]]
Out[25]> {{1,0},{0,1}}

In [26]> Pi
Out[26]> Pi

In [27]> N[%,30]
Out[27]> 3.14159265358979323846264338328
```

The following sections provide an overview only.
For function specifications and implementation details, see the [reference](docs/reference.md) or the Markdown documents in the `docs` directory.

`N` results retain certified enclosures rather than only display text. Certified decimal approximations can therefore participate in ordinary `+ - * /` together with exact Numbers; propagated uncertainty may reduce the reported accuracy, and an outer `N` never reconstructs digits that were not guaranteed by the input approximation.

## v1.5.2

v1.5.2 keeps the exact-first numerical foundation of v1.5.1 while substantially expanding **symbolic calculus, special functions, precision-aware numerical evaluation, Arrays, and linear algebra**.

Highlights:

- function calls are now canonicalized to `name[...]`; `()` is grouping only
- `In [n]` / `Out[n]` and the `@` / `%` shorthands provide absolute and relative history access
- integration Knowledge now covers broader trigonometric powers, reciprocal powers, product-to-sum, bounded Weierstrass substitution, inverse-chain matching, and related families
- added `fresnelc/fresnels`, `hypergeometric1F1/2F1`, incomplete elliptic integrals, and `Ei/Si/Ci/li/polylog`, connected to `D`, `integrate`, and certified `N` where supported
- `N[expr,p]` is precision-aware, allowing FFT and Matrix operations to enter BigFloat/interval backends without first materializing huge exact intermediate expressions
- `{...}` is a general brace container; rectangular values are automatically optimized into dense row-major Arrays
- added Bareiss fraction-free elimination, `solveLinear`, `nullSpace`, LU, rectangular Householder QR, real/complex SVD, and Eigen/Schur support
- large dense Matrix benchmarks show that at order 1024 the current Expr/Rational representation and parse/storage cost become major bottlenecks before the algorithms themselves
- release state: internal **2027 / 2027** PASS, black-box **1504 / 1504** PASS, with fixed-seed Matrix/FFT invariants also passing

See [`CHANGELOG.md`](CHANGELOG.md) for release details and [`docs/performance_optimization.md`](docs/performance_optimization.md) for measured adoption/rejection rationale.

## 1. Getting started

On Windows, simply launch `mmCal.exe`.

The display precision and default angle unit can also be specified at startup.

```text
mmCal --fix 16 --angle deg
mmCal --angle rad
mmCal --angle grad --fix 8
```

- `--fix 16`: Display results with **up to 16 digits after the decimal point**. Unnecessary trailing zeros are omitted
- `--angle deg`: Treat trigonometric inputs without an explicit angle unit as degrees
- `--angle rad`: Radians. This is the default
- `--angle grad`: Gradians
- `--help`, `-h`: Show startup options

On Linux and similar systems, mmCal can be built from source using CMake 3.20 or later with GCC or Clang.

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## 2. Exact values and decimal display are different things

`0.1` is not immediately converted to a binary floating-point value; it is treated as the exact value `1/10`.
Large integers, rational numbers, algebraic expressions containing radicals, complex numbers, and symbolic expressions are kept exact whenever possible. Expansion, factorization, simplification, differentiation, integration, limits, and equation solving all operate within the same expression system.

```text
In [1]> 1/3
Out[1]> 1/3

In [2]> sqrt[2]
Out[2]> sqrt[2]
```

Values such as `Pi`, `E`, `Phi`, and `sqrt[2]` are not treated as pre-stored machine-precision floating-point constants. mmCal distinguishes exact symbolic expressions from numerical approximations.

When a numerical approximation is required, use `N[expr,p]`. `p` is the number of **significant decimal digits**, not a fixed number of digits after the decimal point. Fixed-decimal presentation belongs to `:fix` / `--fix`.

```text
In [3]> N[1/3,20]
Out[3]> 0.33333333333333333333

In [4]> N[Pi,30]
Out[4]> 3.14159265358979323846264338328
```

`N[expression,p]` evaluates the expression to a certified approximation at `p` significant digits, so requested relative precision follows the scale of the value.
By contrast, `:fix` changes **only the number of displayed fractional digits** and does not change the exact value stored internally or the Precision semantics of `N`.

```text
:fix 6
Display: Fixed(6)

In [5]> 1/3
Out[5]> 0.333333

:fix off
Display: Exact

In [6]> Out[5]
Out[6]> 1/3
```

In other words, even after using `:fix`, the value stored in `Out[5]` remains `1/3`.
The value being calculated and its presentation are kept separate.

### precision / accuracy / rationalize

Approximate values produced by `N` retain a truth-certifying CertifiedEnclosure and a separate InformationEnclosure that limits how much information later computations may reuse. Approximate values can be fed back into ordinary arithmetic and certified scalar functions such as `sin`, `exp`, `log`, and `sqrt`; cancellation and error propagation naturally reduce the resulting Accuracy / Precision, while an outer `N` cannot recover undeclared guard digits.

- `accuracy[x]`: An integer lower bound on the **guaranteed number of absolute decimal digits** relative to the true value
- `precision[x]`: An integer lower bound on the **guaranteed number of relative decimal digits** relative to the true value
- `rationalize[x]`: Recover an exact rational number from the certified interval of an approximate value

```text
In [7]> accuracy[N[1/3,20]]
Out[7]> 20

In [8]> precision[N[1/3,20]]
Out[8]> 19

In [9]> rationalize[N[1/3,20]]
Out[9]> 1/3

In [10]> accuracy[1/3]
Out[10]> Infinity
```

`precision` and `accuracy` do not simply return the `n` supplied to `N[...,n]`.
They are conservatively derived from guaranteed error bounds, so the result may be smaller than the requested number of digits.
Exact values and exact symbolic expressions return `Infinity` in this sense.

Numerical approximation does not merely repeat a calculation until the result "looks close enough."
The basic approach is to compute an interval containing the true value and confirm that the required decimal representation is uniquely determined.

## 3. Basic syntax

Function calls use **square brackets `[]` only**. Parentheses `()` are reserved for expression grouping and are not function-call delimiters.

```text
sin[Pi/6]
sqrt[2]
log[10,1000]
```

Therefore `sin(Pi/6)` is not a function call. Using the legacy `()` syntax on a known function name is a SyntaxError; write `sin[Pi/6]`. Adjacency such as `x(x+1)` for an ordinary identifier is accepted as implicit multiplication, while the Formatter canonicalizes it to the explicit `x*(x+1)`.

Ordinary notation is used for arithmetic and powers.

```text
2+3*4
(x+1)^2
1/(x+1)
2Pi
2sqrt[2]
```

Finite decimal literals are parsed as exact fractions.

```text
0.125
-> 1/8
```

Arrays use braces.

```text
{1,2,3}
{{1,2},{3,4}}
```

## 4. Angles

The default angle unit is radians.

```text
sin[Pi/6]
-> 1/2
```

The session default can be inspected or changed with `angleMode`.

```text
angleMode[]
-> Rad

angleMode[Deg]
-> Deg

sin[30]
-> 1/2
```

An angle unit can also be specified explicitly for part of an expression.

```text
sin[30 Deg]
sin[Pi/6 Rad]
```

Explicit `Deg` / `Rad` / `Grad` units take precedence over the session default.

## 5. Variables and user-defined functions

```text
x:=2
-> 2

f[t]:=t^2+1
f[4]
-> 17
```

Redefining an existing definition produces a notification.

```text
x:=4
INFO: x redefined (was 2)
-> 4
```

List current definitions:

```text
Defs[]
```

Remove definitions:

```text
UnDef[x]
UnDef[x,y,f]
```

Remove definitions and history together:

```text
Clear[]
```

Exit:

```text
Exit[]
```

## 6. Input/output history

The most recent successful output is referenced with `%`, while the immediately previous input expression is referenced with `@`. Repeating the shorthand moves further back: `%%` / `%%%` for earlier successful outputs and `@@` / `@@@` for earlier inputs. There is no fixed shorthand depth limit.

```text
In [1]> 2+3
Out[1]> 5

In [2]> %*2
Out[2]> 10

In [3]> Pi
Out[3]> Pi

In [4]> N[@,30]
Out[4]> 3.14159265358979323846264338328
```

The formal history interface is `In [n]` / `Out[n]`: positive indices are absolute and negative indices are relative. Zero is invalid.

```text
In [1]
Out[1]
In [-1]    // previous input; equivalent to @
Out[-1]   // previous successful output; equivalent to %
```

`Out[n]` returns a stored output snapshot without reevaluation. Negative `Out[-n]` counts successful outputs only, so `% == Out[-1]` remains true across failed evaluations.
`In [n]` retrieves a previous lowered input expression and **evaluates it again in the current definition environment**. Negative `In [-n]` counts input slots, so `In [-1]` / `@` reevaluates the immediately previous input.

## 7. Main mathematical features

### Basic functions and complex numbers

```text
abs[-3]              -> 3
abs[3+4I]            -> 5
re[3+4I]             -> 3
im[3+4I]             -> 4
conj[3+4I]           -> 3-4I
arg[-1]              -> Pi Rad
```

`sqrt`, `log`, general powers, inverse trigonometric functions, inverse hyperbolic functions, and similar operations use principal values, including on the complex domain.

```text
sqrt[-1]             -> I
log[-1]              -> I Pi
```

### Trigonometric and hyperbolic functions

Implemented functions include `sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `asin`, `acos`, `atan`, `atan2`, as well as `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, and others.

### Symbolic differentiation

```text
D[x^3+sin[x],x]
-> 3x^2+cos[x]

D[x^5,{x,3}]
-> 60x^2
```

### Symbolic integration

Indefinite integrals do not display `+C`; one representative antiderivative is returned.

```text
integrate[x^2+sin[x],x]
-> x^3/3-cos[x]

integrate[x cos[x],x]
-> x sin[x]+cos[x]

integrate[1/(1+x^2),{x,0,1}]
-> Pi/4

integrate[exp[-x],{x,0,Infinity}]
-> 1
```

v1.5.2 expands integration as shared Knowledge rather than as isolated rules: finite-Fourier reduction for integer trigonometric powers, sec/csc recurrences for negative powers, cross-frequency product-to-sum, bounded Weierstrass substitution, inverse-chain candidates, quadratic radicals, and related families. Failure diagnostics distinguish unsupported rules, partial evaluation, missing conditions, and recognized families with no known finite closed form in the current standard-function vocabulary. Proof limitations alone do not remove a useful primitive.

### Special functions added in v1.5.2

The following special-function foundations are connected to symbolic differentiation, integration, and certified numerical evaluation where supported:

```text
fresnelc[x]  fresnels[x]
hypergeometric1F1[a,b,z]
hypergeometric2F1[a,b,c,z]
ellipticF[phi,m]  ellipticE[phi,m]  ellipticPi[n,phi,m]
Ei[x]  Si[x]  Ci[x]  li[x]  polylog[s,z]
```

Unreleased also adds `zeta`, `digamma`, `trigamma`, and the regularized incomplete beta `ibeta`, connected to representative exact reductions, derivative relations, and certified `N` on their supported real domains. Lightweight exact number theory now includes `isprime`, `nextprime`, `prevprime`, `factorint`, and `totient`, deterministic over the `uint64` range; larger BigInts are not promoted from probable-prime evidence to certified truth.

Representative reductions include:

```text
integrate[exp[-x^2],x]
-> erf[x]sqrt[Pi]/2

integrate[exp[-x^2],{x,-Infinity,Infinity}]
-> sqrt[Pi]

integrate[exp[x^6],x]
-> x hypergeometric1F1[1/6,7/6,x^6]

integrate[1/(1+x^5),x]
-> x hypergeometric2F1[1,1/5,6/5,-x^5]

integrate[sin[x]/x,x] -> Si[x]
integrate[cos[x]/x,x] -> Ci[x]
integrate[exp[x]/x,x] -> Ei[x]
integrate[1/log[x],x] -> li[x]
integrate[log[1-x]/x,x] -> -polylog[2,x]
```

mmCal does not invent general inverse special functions for `solve`; when branch structure or injectivity cannot be established, the equation remains unresolved.

### Limits

```text
limit[sin[x]/x,x,0]       -> 1
limit[1/x,x,0,1]          -> Infinity
limit[1/x,x,0,-1]         -> -Infinity
```

### Equations and inequalities

```text
solve[x^2-2==0,x]
-> {x==sqrt[2], x==-sqrt[2]}

solve[x^2<4,x,Real]
-> {x in Real if x>-2&&x<2}

solve[exp[x]==2,x,Real]
-> {x==log[2]}

solve[sin[x]==0,x,Real]
-> {x==Pi k where k in Integer}
```

Real periodic `sin/cos/tan` equations can now return integer-parameter families when the argument is affine in the solve variable with an exact nonzero linear coefficient. Nonlinear arguments and complete Complex-domain families are not fabricated from a principal inverse.

For higher-degree Rational-coefficient polynomials over Real, the solver can fall back to exact `root` values when the existing radical-oriented solver does not close naturally. `root[{a0,...,an},k]` denotes the `k`th distinct increasing real root of the ascending-power coefficient polynomial, and `N` certifiedly refines its Sturm isolating interval.

```text
solve[x^5-x+1==0,x,Real]
-> {x==root[{1,-1,0,0,0,1},1]}

N[root[{-2,0,1},2],30]
-> 1.41421356237309504880168872421
```

Complex roots are represented by `root[{a0,...,an},k,Complex]`, using certified isolating disks and exact Root fallback for high-degree Rational-polynomial Complex Solve. Individual Root construction now reduces the defining polynomial to the selected proven Rational irreducible factor when bounded factorization succeeds (degree at most 16). For Root-to-Root `+ - * /`, a proof-based primitive-element path is used when the operand minimal polynomials are proven irreducible and `theta=alpha+c beta` is certified to generate an irreducible extension of product degree; otherwise arithmetic falls back to exact resultants plus certified root re-identification. General algebraic equality/ordering, `rootApproximant`, and complete Q-factorization beyond the current budget remain deliberate future work.

When mmCal cannot guarantee a complete solution set, it does not return an arbitrary convenient solution as though it were complete.
Instead, it reports the unresolved state using a Warning and the result representation.

## 8. Arrays, matrices, vectors, and statistics

`{...}` is a general finite brace container rather than matrix-only syntax. Values whose children share one shape are automatically optimized to the dense **shape + row-major flat storage** `ArrayExpr`; heterogeneous-shape values such as `{Q,R}` remain general braces. Matrix functions audit rectangularity at their boundary and leave non-rectangular values unevaluated with a Warning. Only empty shapes that cannot be preserved by braces alone are formatted through `reshape`.

```text
dimensions[{{1,2,3},{4,5,6}}]
-> {2, 3}

arrayRank[{{1,2},{3,4}}]
-> 2

at[{{1,2},{3,4}},1]
-> {3, 4}

at[{{1,2},{3,4}},1,0]
-> 3

reshape[{1,2,3,4},{2,2}]
-> {{1, 2}, {3, 4}}

zeros[0,3]
-> reshape[{}, {0, 3}]
```

Indices are zero-based. `arrayRank[A]` means the number of Array dimensions, while `matrixRank[A]` is the linear-algebra rank.

Use `range` / `table` / `map` for finite sequence generation and explicit element-wise application. `table` iterators are locally scoped and do not overwrite an outer definition. `map` visits scalar leaves explicitly, so ordinary calls such as `exp[A]` are not silently given element-wise semantics.

```text
range[0,1,1/3] -> {0,1/3,2/3,1}
table[i^2,{i,5}] -> {1,4,9,16,25}
map[sin,{0,Pi/2,Pi}] -> {0,1,0}
```

The canonical basic linear-algebra API is:

```text
dot[{{1,2},{3,4}},{{5,6},{7,8}}]
-> {{19, 22}, {43, 50}}

dot[{1,2},{3,4}]
-> 11

det[{{1,2},{3,4}}]
-> -2

inverse[{{1,2},{3,4}}]
-> {{-2, 1}, {3/2, -1/2}}

rref[{{1,2},{3,4}}]
-> {{1, 0}, {0, 1}}

matrixRank[{{1,2},{2,4}}]
-> 1

nullSpace[{{1,2},{2,4}}]
-> {{-2, 1}}

solveLinear[{{2,1},{1,-1}},{5,1}]
-> {2, 1}

luDecomposition[{{0,2},{3,4}}]
-> {{{0,1},{1,0}},{{1,0},{0,1}},{{3,4},{0,2}}}

qrDecomposition[{{3,0},{4,0}}]
-> {{{-3/5,-4/5},{-4/5,3/5}},{{-5,0},{0,0}}}

svd[{{3,0},{0,4}}]
-> {{{0,1},{1,0}},{{4,0},{0,3}},{{0,1},{1,0}}}

eigenvalues[{{0,-1},{1,0}}]
-> {I, -I}

norm[{3+4I}]
-> 5

normalize[{3,4}]
-> {3/5, 4/5}
```

`dot` covers vector-vector, matrix-vector, vector-matrix, and matrix-matrix contraction.
`A*B` is deliberately not matrix multiplication: ordinary arithmetic supports same-shape Array `+` / `-` and scalar×Array, while matrix multiplication and inner products remain explicit as `dot[A,B]`.

Integer and Rational matrices remain exact rather than being converted to BigFloat. `det`, `rref`, `matrixRank`, `nullSpace`, `inverse`, and `solveLinear` share per-row denominator clearing and Bareiss fraction-free elimination to suppress intermediate Rational growth; exact complex matrices fall back to the `Number` Gaussian backend. `solveLinear[A,b]` returns only unique solutions, including consistent overdetermined systems with full column rank. Inconsistent systems and systems with free variables are Domain errors rather than invented parametric answers. General symbolic determinant/inverse expansion has a work budget: triangular and sufficiently sparse cases are still evaluated, while potentially explosive dense cases remain unevaluated instead of constructing factorial-size expressions. `luDecomposition[A]` returns `{P,L,U}` for square matrices. `qrDecomposition[A]` is a rectangular reduced Householder QR: for m×n input with `k=min[m,n]`, it returns `{Q,R}` with `Q:m×k` and `R:k×n`. `svd[A]` likewise returns rectangular reduced `{U,S,V}`; its general numerical backend avoids forming `A^H A` and instead uses Householder bidiagonalization plus one-sided Jacobi. Factors are extracted with prefix indexing such as `at[result,0]`. General exact QR is limited to 3x3 to prevent expression explosion, while upper-triangular/trapezoidal cases retain an any-size fast path. `eigenvalues`, `eigenvectors`, and `eigensystem` target square matrices: exact triangular/diagonal and distinct-root exact Number 2x2 cases remain exact, while general `N[...]` uses Hessenberg reduction plus implicit shifted complex QR and audits Schur/eigenpair relations against the certified input intervals. Defective or near-multiple cases are not supplied with guessed independent eigenvectors.

Under `N`, matrix operations dispatch directly to a precision-aware backend just like FFT:

```text
N[det[{{Pi,0},{0,2}}],12]
-> 6.28318530718

N[inverse[{{Pi,0},{0,2}}],12]
-> {{0.318309886184, 0}, {0, 0.5}}

N[solveLinear[{{Pi,0},{0,2}},{Pi,4}],12]
-> {1, 2}

N[qrDecomposition[{{1,2},{3,4}}],8]
-> {{{-0.31622777,-0.94868330},{-0.94868330,0.31622777}},{{-3.16227766,-4.42718872},{0,-0.63245553}}}

N[eigenvalues[{{1,2},{3,4}}],8]
-> {-0.37228132, 5.3722813}
```

This avoids first constructing a huge exact result. Requested precision is handled directly with certified BigFloat/interval operations. `solveLinear` likewise tries direct augmented interval elimination and never substitutes an epsilon guess when pivots or consistency cannot be certified; dependent overdetermined rows can be intrinsically difficult to certify because interval evaluation loses correlation. `matrixRank` and `nullSpace` are discontinuous with respect to rank deficiency, so exact inputs use exact elimination first. Approximate inputs still use no arbitrary epsilon: the interval backend returns a result only when the pivot structure is certified rather than guessing rank deficiency.

For `N` output, exact terminating decimals remain compact (`N[1/2,10] -> 0.5`). Certified-interval results compact only redundant runs of trailing zeros while retaining one visible zero, so `1.000000000000 -> 1.0` and `1.500000000000 -> 1.50`. Requested digits and the certified enclosure remain in metadata.

Legacy `matmul` / `mmul` / `vdot`, `rank` / `mrank`, `vnorm`, `vnormalize`, and `mget` remain compatibility aliases.

Statistical operations also retain exact numbers and symbolic expressions whenever possible.

```text
mean[{1,2,4}]
-> 7/3

var[{1,2,3}]
-> 2/3

stddev[{1,2,3}]
-> sqrt[6]/3

corr[{1,2,3},{2,4,6}]
-> 1
```

## 9. FFT and random numbers

```text
dft[{1,2,3,4}]
fft[{1,2,3,4}]

ifft[fft[{1,2,3,4}]]
-> {1, 2, 3, 4}

convolve[{1,2},{3,4}]
-> {3, 10, 8}
```

Random-number state is maintained per session.

```text
randSeed[42]
rand[]
randint[1,6]
choice[{a,b,c}]
randn[]
```

Setting the same seed reproduces the same sequence.
The random-number generator is **not intended for cryptographic use**.

## 10. About Warnings

When `D`, `integrate`, `limit`, `solve`, or similar operations cannot prove a complete evaluation, mmCal may return the expression unevaluated together with a Warning rather than inventing an incorrect value.

```text
integrate[gamma[x],x]
WARN: integrate could not fully prove the symbolic antiderivative or definite integral; unevaluated integrate[...] remains
-> integrate[gamma[x], x]
```

This does not mean that the expression itself is the mathematical result. It means that **the current implementation could not safely complete the calculation**.
Unsupported regions will continue to be expanded. I'm working on it.

## 11. CLI display settings

```text
:fix 16
:fix off
:status
```

`:fix n` only rounds the **display** to at most `n` digits after the decimal point; it does not change the stored value or the semantics of `precision` / `accuracy`.
`:status` shows the current angle mode, display mode, number of definitions, history count, and similar state.
The console title also shows the angle and display mode as auxiliary information.

## 12. Naming

Ordinary mathematical functions use lowercase canonical names.

```text
sin cos log sqrt integrate solve
```

Short symbolic operations and session operations use the following proper names.

```text
D N In Out Exit Clear Defs UnDef
```

## 13. Detailed documentation

- `docs/reference.md` — Detailed function, syntax, and current-behavior reference
- `docs/mathematics.md` — Mathematical policy for domains, principal values, and numerical evaluation
- `docs/architecture.md` — Internal architecture for developers
- `docs/grammar.ebnf` — Machine-readable overview of the grammar
- `docs/performance_optimization.md` — Performance work adopted or rejected for v1.5.1–v1.5.2, with benchmark rationale
- `docs/multiprecision_implementation.ja.md` — Detailed Japanese notes on the multiprecision / certified numerical backend
- `CHANGELOG.md` — Major changes by release

## 14. License and trademarks

The source code is distributed under the **BSD 3-Clause License**. See `LICENSE` for the copyright permissions governing commercial use, modification, redistribution, and incorporation into other software.

**Use of the mmCal name, official logos, and related branding is addressed separately by `TRADEMARKS.md` .** The trademark policy does not prohibit independent GUIs, forks, commercial products, or other uses permitted by the BSD license; it is intended to avoid confusion about whether a third-party product is mmCal itself or is officially endorsed.

- [BSD 3-Clause License](LICENSE)
- [Trademark and Brand Policy](TRADEMARKS.md)

If mmCal is used in an academic publication or product, attribution beyond the BSD requirements is not mandatory, but a factual acknowledgement is appreciated.

## 15. Tests and development environment

At the v1.5.2 release point, this project contains 2,027 internal regression tests and 1,504 black-box tests.
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

## Requests and feedback

Please submit requests, bug reports, and implementation proposals on GitHub. They are very welcome 🍀

## Future ideas

- It would be nice to have something like `for`
- A `plot` function for graphing
