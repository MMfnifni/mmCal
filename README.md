# mmCalculator – Mathematical Machinery Calculator

An exact-first CLI calculator and compact CAS for engineering, research, and manufacturing.

© 2021–2026 mmKreutzef (aka Daiki.NIIMI)  
Licensed under the BSD 3-Clause License

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

In [3]>  0.1+0.2==0.3
Out[3]> True

In [4]> 1/3+1/6
Out[4]> 1/2

In [5]> sqrt[72]
Out[5]> 6sqrt[2]

In [6]> sIn [Pi/6]
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

In [12]> integrate[x^2+sIn [x],x]
Out[12]> x^3/3-cos[x]

In [13]> integrate[sIn [x],{x,0,Pi}]
Out[13]> 2

In [14]> limit[(1-cos[x])/x^2,x,0]
Out[14]> 1/2

In [15]> solve[x^2+1==0,x,Complex]
Out[15]> {x==I, x==-I}

In [16]> (1+I)/(1-I)
Out[16]> I

In [17]> sqrt[-8]
Out[17]> 2I sqrt[2]

In [18]> Pi
Out[18]> Pi

In [19]> N[%,30]
Out[19]> 3.141592653589793238462643383280

In [20]> N[%%%,20]
Out[20]> 2.82842712474619009760I
```

The following sections provide an overview only.
For function specifications and implementation details, see the [reference](docs/reference.md) or the Markdown documents in the `docs` directory.

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

When a decimal value is required, use `N[expr,n]` to request an arbitrary-precision numerical approximation.

```text
In [3]> N[1/3,20]
Out[3]> 0.33333333333333333333

In [4]> N[Pi,30]
Out[4]> 3.141592653589793238462643383280
```

`N[expression,digits]` returns a numerical approximation of the expression itself.
By contrast, `:fix` changes **only how values are displayed** and does not change the exact value stored internally.

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

Approximate values produced by `N` retain not only a display string but also precision metadata such as a certified enclosure containing the true value.

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

Square brackets are the standard notation for function calls. Parentheses can also be used for calls.

```text
sIn [Pi/6]
sqrt[2]
log[10,1000]
```

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
sIn [Pi/6]
-> 1/2
```

The session default can be inspected or changed with `angleMode`.

```text
angleMode[]
-> Rad

angleMode[Deg]
-> Deg

sIn [30]
-> 1/2
```

An angle unit can also be specified explicitly for part of an expression.

```text
sIn [30 Deg]
sIn [Pi/6 Rad]
```

Explicit `Deg` / `Rad` / `Grad` units take precedence over the session default.

## 5. Variables and user-defined functions

```text
x:=2
-> 2

f(t):=t^2+1
f(4)
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

The most recent successful result can be referenced with `%`, the one before that with `%%`, and earlier results by adding more `%` characters such as `%%%`.

```text
In [1]> 2+3
Out[1]> 5

In [2]> %*2
Out[2]> 10
```

Absolute indices are also available.

```text
In [1]
Out[1]
```

`Out[n]` returns the output stored at that point.
`In [n]` retrieves the previous input expression and **evaluates it again in the current definition environment**.

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
D[x^3+sIn [x],x]
-> 3x^2+cos[x]

D[x^5,{x,3}]
-> 60x^2
```

### Symbolic integration

Indefinite integrals do not display `+C`; one representative antiderivative is returned.

```text
integrate[x^2+sIn [x],x]
-> x^3/3-cos[x]

integrate[x cos[x],x]
-> x sIn [x]+cos[x]

integrate[1/(1+x^2),{x,0,1}]
-> Pi/4

integrate[exp[-x],{x,0,Infinity}]
-> 1
```

### Limits

```text
limit[sIn [x]/x,x,0]       -> 1
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
```

When mmCal cannot guarantee a complete solution set, it does not return an arbitrary convenient solution as though it were complete.
Instead, it reports the unresolved state using a Warning and the result representation.

## 8. Matrices, vectors, and statistics

```text
matmul[{{1,2},{3,4}},{{5,6},{7,8}}]
-> {{19, 22}, {43, 50}}

det[{{1,2},{3,4}}]
-> -2

inverse[{{1,2},{3,4}}]
-> {{-2, 1}, {3/2, -1/2}}

rref[{{1,2},{3,4}}]
-> {{1, 0}, {0, 1}}

vdot[{1,2},{3,4}]
-> 11

vcross[{1,0,0},{0,1,0}]
-> {0, 0, 1}

vnorm[{3,4}]
-> 5

mean[{1,2,4}]
-> 7/3

var[{1,2,3}]
-> 2/3

stddev[{1,2,3}]
-> sqrt[6]/3

corr[{1,2,3},{2,4,6}]
-> 1
```

Matrix and statistical operations also retain exact numbers and symbolic expressions whenever possible.
`N[{...}]` recursively approximates array elements numerically.

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
- `docs/roadmap.md` — Major currently unsupported features and future candidates
- `docs/grammar.ebnf` — Machine-readable overview of the grammar

## 14. License

BSD 3-Clause

Copyright (c) 2021–2026 mmKreutzef

See `LICENSE` for the full terms.

If mmCal is used in an academic paper, product, or similar work, I would appreciate a note in the documentation or publication acknowledging its use, although **this is not an additional license requirement**.
I'd be even happier if you let me know.

Reselling the source with little more than a different platform wrapper or UI is permitted by BSD 3-Clause, but it does make the author a little sad.

Embedding mmCal as part of another software system is very welcome.

## 15. Tests and development environment

This project contains more than 1,500 internal regression tests and more than 1,000 black-box tests.
They focus especially on exact arithmetic, boundary values, domains, error classification, formatter round-trip parsing, and certified numerical enclosures.

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
