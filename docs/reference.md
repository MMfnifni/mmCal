# mmCal Specification and Function Reference

This document is the detailed specification of **mmCal as implemented**.
For a user-oriented introduction, see the root-level `README.md`.

> Target: **v1.5.4 development tree**

This document changes with each version; retrieve older versions from the Git history when needed.

## 0. Unchanging Principles

### Design Principles

Be rigorous.
Do everything from scratch.
Explicitly state any approximations.
Admit when you don’t know something.

### Distribution Principles

Deliver the system as a single executable.
Open source under the BSD 3-Clause License.

## 1. Current design philosophy

mmCal is not a calculator that immediately lowers every input to `double`. It is designed as an **exact-first compact CAS / numerical-computation kernel**.

Priorities are:

1. Mathematical correctness
2. Preserving principal branches, definedness, and domains
3. Never conflating exact and approximate values
4. Reusing shared mathematical knowledge across the Solver, Simplifier, and CertifiedEvaluator
5. Improving performance only where the above semantics remain intact

Representative examples:

```text
0.1 + 0.2
-> 3/10

sin[Pi/6]
-> 1/2

sqrt[-8]
-> 2I sqrt[2]

log[-1]
-> I Pi

N[Pi,30]
-> 3.14159265358979323846264338328
```

`Pi` and `sqrt[2]` are not decimals stored internally. They remain exact expressions, and certified numerical evaluation is performed only when `N[...]` is requested.

---

# 2. Internal model for numbers and expressions

## 2.1 Integer

Arbitrary-length integer `BigInt`.

```text
123456789012345678901234567890
-> 123456789012345678901234567890

100!
-> exact BigInt
```

Values are not rounded to fixed-width 64-bit integers.

## 2.2 Rational

Finite decimal literals are parsed directly as exact Rational values.

```text
0.1
-> 1/10

1.25
-> 5/4

0.1 + 0.2 == 0.3
-> True
```

Therefore the usual IEEE 754 `0.30000000000000004` issue does not occur in ordinary exact evaluation.

## 2.3 Exact complex

Both the real and imaginary parts are stored as exact `Number` values.

```text
I^2
-> -1

(3 + 4I) / 5
-> 3/5 + 4/5I
```

## 2.4 Symbolic expressions

Expressions that do not close to exact numeric values remain as AST expressions.

```text
sqrt[2]
-> sqrt[2]

sin[1]
-> sin[1]

gamma[1/3]
-> gamma[1/3]
```

Here, the `1` in `sin[1]` means **1 radian**. The default angle unit is radians.

## 2.5 Certified approximation

`BigFloat` is used for arbitrary-precision working values, while `RealInterval` / `ComplexInterval` provide certified enclosures.

For `N[expr,p]`, mmCal confirms that both endpoints of an interval containing the true value round to the same `p`-significant-digit decimal result before returning a `DecimalApproximation`. Near zero, where relative Precision may be unavailable, a zero-centered approximation can still be returned when its InformationEnclosure proves useful absolute Accuracy.

A `DecimalApproximation` is not merely a display string. It retains the requested digit count, provenance (exact input / certified interval), the displayed decimal value as an exact Rational, and two exact Rational intervals: a **CertifiedEnclosure** and an **InformationEnclosure**. The CertifiedEnclosure proves containment of the true value. The InformationEnclosure limits how much information later computations may legitimately reuse, with the invariant `CertifiedEnclosure ⊆ InformationEnclosure`. `ComplexDecimalApproximation` keeps the same metadata independently for its real and imaginary components. `precision`, `accuracy`, and `rationalize` use the InformationEnclosure directly rather than reparsing the display string and guessing its quality.

```text
N[sqrt[2],30]
-> 1.41421356237309504880168872421
```

The evaluator does not rely solely on heuristic stopping conditions such as "the difference became sufficiently small."

---

# 3. Predefined symbols

Current protected predefined symbols:

| Name            | Meaning                                                                                                                                          |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `Pi`            | Circle constant. Exact transcendental constant                                                                                                   |
| `E`             | Base of the natural logarithm. Exact transcendental constant                                                                                     |
| `Phi`           | Golden ratio. Exact algebraic constant                                                                                                           |
| `I`             | Imaginary unit                                                                                                                                   |
| `True`, `False` | Boolean values                                                                                                                                   |
| `Integer`       | Integer domain                                                                                                                                   |
| `Rational`      | Rational domain                                                                                                                                  |
| `Real`          | Real domain                                                                                                                                      |
| `Complex`       | Complex domain                                                                                                                                   |
| `Infinity`      | Positive extended-real infinity; also returned by `precision/accuracy` for exact values                                                          |
| `ComplexInfinity` | Infinite magnitude with undetermined direction, produced by a provably nonzero value divided by exact zero                                     |
| `Indeterminate` | Protected nonnumeric result for an expression whose numerical value is not unambiguously defined                                                  |

The legacy constants `Tau`, `NA`, and `ESP` are no longer predefined.

---

# 4. Input syntax

## 4.1 Function calls

Function calls use **square brackets `[]` only**. Parentheses `()` are reserved for grouping and are not function-call delimiters.

```text
sin[Pi/6]
sqrt[2]
f[x]
```

Using the legacy parenthesized form such as `sin(x)` on a known function name is a SyntaxError. For an ordinary identifier, `x(x+1)` is accepted as implicit multiplication, but the Formatter emits the explicit canonical form `x*(x+1)`. Parentheses are used for grouping, a standalone `[x+1]` is not a grouping expression, and user-defined function signatures use `f[x] := ...`.

## 4.2 Arrays

```text
{1,2,3}
{{1,2},{3,4}}
```

Internally, dense values use `ArrayExpr`. Numeric values may live in immutable packed pages while shape/offset/strides form a separate layout, allowing transpose and some reshape/slice operations to share the backing storage. This is an implementation detail and does not create new user-visible Array types.

## 4.3 Variables and user-defined functions

```text
x := 3
-> 3

f[t] := t^2 + 1
f[4]
-> 17
```

For variables, `:=` behaves like `Set`: the right-hand side is evaluated and the result is stored.
Function definitions behave like `SetDelayed`: the function body itself is retained.

Overwriting an existing definition with different content emits an Info diagnostic separately from the evaluation result.

```text
x := 2
-> 2

x := 4
INFO: x redefined (was 2)
-> 4
```

The same applies when redefining a user function with the same arity. Reassigning exactly the same value does not emit an Info diagnostic.

Built-in names and predefined symbols cannot be redefined.

Use the following to inspect or remove current definitions:

```text
Defs[]
UnDef[x]
UnDef[x,y,f]
```

`Defs[]` returns current global user variables and user-function definitions as an array of expressions. `UnDef[...]` removes both variable and function definitions for the specified names and returns the number of names actually changed.

## 4.4 History

The shorthand forms are:

```text
@       // In [-1]
@@      // In [-2]
@@@     // In [-3]
%       // Out[-1]
%%      // Out[-2]
%%%     // Out[-3]
```

The formal interface is `In [n]` / `Out[n]`. `n > 0` is an absolute prompt index, `n < 0` is relative, and `n = 0` raises TypeError.

```text
In [1]
Out[1]
In [-1]
Out[-1]
```

`In [n]` retrieves the lowered Expr for the target input and then **evaluates it normally in the current session environment**. Positive `In [n]` uses an absolute input number. Negative `In [-n]` counts previous input slots while excluding the input currently being evaluated. Therefore `@` / `@@` / `@@@` / ... mean `In [-1]` / `In [-2]` / `In [-3]` / ... respectively, with no fixed shorthand depth limit. An input that reached parse/lower but failed during evaluation can therefore be retried through `In [-1]`. By contrast, an input rejected by the Lexer, Parser, or Lowerer is never committed to history and does not consume an `In[n]` number. The error still refers to the pending input number, and the corrected next input reuses that same prompt number.

`Out[n]` returns a stored result snapshot without reevaluation. Positive `Out[n]` is indexed by the absolute input number, while negative `Out[-n]` counts **successful outputs only** from the most recent one. Therefore `% == Out[-1]`, `%% == Out[-2]`, `%%% == Out[-3]`, ... remain true even when failed evaluations occur between successful outputs. Repeated `%` also has no fixed shorthand depth limit.
Explicit history references are also resolved inside held symbolic operators such as `D`, `integrate`, and `limit`. `Out[n]` / `%` splice the stored result snapshot into the held expression, while `In[n]` / `@` re-evaluate the stored input in the current session environment. Thus `integrate[Out[1],x]` and nested forms such as `integrate[2*In[1],x]` are not treated as unknown functions named `Out` or `In`.

```text
In [1]> fft[{1,2,3}]
Out[1]> {6, ...}
In [2]> N[@,30]
Out[2]> {6, -1.50+0.866025403784...I, ...}
```

Here `N[@,30]` does not merely approximate the stored `Out[1]`; it reevaluates the `fft[...]` from `In [1]` inside a 30-digit approximation context.

Absolute-reference example:

```text
In [1]> 1+1
Out[1]> 2
In [2]> 2+2
Out[2]> 4
In [3]> In [1]+In [2]
Out[3]> 6
In [4]> In [3]
Out[4]> 6
```

Thus, `In [n]` means "paste the previous input back into the current environment and execute it again." If the previous input contains variable references, assignment, randomness, or other stateful behavior, the current definitions and RNG state are used. This is deliberately separate from displaying a raw historical AST. A positive absolute reference to the input currently being evaluated is forbidden to prevent direct self-recursion.

If evaluation fails but parsing/lowering succeeded, the absolute input slot remains, but no corresponding positive `Out[n]` snapshot exists.

## 4.5 Comparisons

```text
<  <=  >  >=  ==  !=
```

When the result can be proven, evaluation returns `True` or `False`. If the result is symbolically undecidable, the Predicate expression is retained.

## 4.6 Operator precedence

Important implementation-level precedence rules:

1. postfix `!`
2. power `^` — right-associative
3. unary `+ -`
4. `* /` and implicit multiplication — processed left-to-right at the same term level
5. `+ -`
6. comparison
7. assignment `:=`

```text
2^3^2
-> 512

-2^2
-> -4
```

## 4.7 Implicit multiplication

```text
2Pi
2(x+1)
(x+1)(x-1)
2 x
2exp[x]
2E^x
```

A function name directly following a number is not parsed as part of a function-call token; it is interpreted as implicit multiplication. An `e/E` immediately following a number is consumed as scientific notation only when actual exponent digits follow.

## 4.8 Radix-prefixed numeric literals

The current Lowerer supports:

```text
0b1010
0o17
0xFF
2#1010
16#FF
```

For `base#digits`, Rational literals are also parsed where supported.
Lexer-level infix bitwise operators such as `&`, `|`, `<<`, and `>>` are not implemented. Use the function APIs described later, including `bitAnd`, `bitOr`, `bitXor`, `bitNot`, `bitShiftLeft`, and `bitShiftRight`.

---

# 5. Conditional evaluation and mathematical cases

```text
if[condition,trueExpr,falseExpr]
cases[value1 if condition1; value2 if condition2; ...]
```

`if[...]` is evaluation control: the condition is evaluated first and only the selected branch is evaluated. Therefore DomainError or random-number consumption in the unselected branch does not occur.

`cases[...]` is a first-class scalar piecewise mathematical expression. Proven-false branches are removed without evaluating their values; a proven-true branch short-circuits; undecidable predicates remain symbolic. `simplify` may select branches from assumptions and `N` approximates branch values without numericalizing predicates. `integrate` / `limit` distribute only when conditions are independent of the calculus variable. `D` additionally differentiates safe open interiors of variable-dependent inequality branches, while retaining an unevaluated boundary `D[...]` unless differentiability at the closed boundary can be proved.

---

# 6. Angle semantics

**The default is radians.**

```text
sin[Pi/6]
-> 1/2

asin[1/2]
-> Pi/6
```

Explicit units:

```text
sin[30 Deg]
-> 1/2

sin[Pi/6 Rad]
-> 1/2

sin[100 Grad]
-> 1
```

Angle-unit suffixes can be applied not only to numeric literals but also to general expressions.

```text
x Deg
Pi/6 Rad
(2x + 1) Grad
```

The session default angle unit can be inspected or changed with `angleMode`. The legacy `:angle` command will not be restored.

```text
angleMode[]
-> Rad

angleMode[Deg]
-> Deg

angleMode[Grad]
-> Grad
```

`angleMode` changes mathematical evaluation state, not presentation state. Explicit `Deg/Rad/Grad` units always override the session default. The Kernel API function `setDefaultAngleUnit()` manipulates the same state.

Angle conversion:

```text
DtoR[180] -> Pi
DtoG[90]  -> 100
RtoD[Pi]  -> 180
RtoG[Pi]  -> 200
GtoD[200] -> 180
GtoR[200] -> Pi
```

---

# 7. `N` — numerical approximation

```text
N[expr]
N[expr,p]
```

`p` is the number of **significant decimal digits**; the default is 16. It is not a fixed number of digits after the decimal point. Fixed-decimal presentation is controlled separately by `:fix` / `--fix`.
`N` applies recursively to Arrays. For explicit angle-unit values such as those returned by `arg`, only the numeric component is approximated and the unit is retained.

`N` is also the entry point for precision-aware evaluation. It resolves the requested precision before evaluating its first argument and keeps that precision context active while the child expression is evaluated. Ordinary builtins still follow exact-first evaluation; only explicitly supported builtins such as FFT consume the context and evaluate directly in a certified approximate domain.

If whole-expression certification is unavailable, ordinary evaluated Calls / Arrays / Lists may still be traversed structurally: only numerically closed subexpressions are approximated, while free symbols and unresolved symbolic function parts remain exact. Calls with `HoldAll` / `HoldFirst`-style semantics are not blindly rebuilt, preserving their evaluation contract.

```text
N[x+Pi,20]
-> 3.1415926535897932385+x

N[sin[x]+Pi,20]
-> 3.1415926535897932385+sin[x]

N[True,20]          -> True
N[Infinity,20]      -> Infinity
N[Indeterminate,20] -> Indeterminate
```

Leaving a free symbol, Boolean, `Infinity`, `ComplexInfinity`, or `Indeterminate` exact is not itself a Warning. If an inner operation such as `D`, `limit`, `solve`, or `rref` has already emitted a more specific Warning, the enclosing `N` suppresses a redundant generic `N::unevaluated`. A mathematically existing value whose certified backend is unavailable produces `N::unsupported` and remains symbolic rather than being reported as DomainError. If increasing guard precision cannot overcome the width of an already finite-precision `InformationEnclosure`, `N::precision` returns the expression after bounded refinement instead of consuming the global refinement budget. Therefore exact `gamma[0]` or a `2F1` denominator parameter known exactly to be a non-positive integer is a DomainError, while `gamma[N[0,5]]` or `hypergeometric2F1[1,2,N[0,5],2]` yields `N::precision` because the finite input information still permits nearby nonsingular values. Likewise, `log[-1+I*N[0,5]]` is not collapsed to one side of its branch cut. Complex Lambert W cases for which the current contraction proof cannot construct a certified box remain held with `N::unsupported`.

```text
N[Pi,20]
-> 3.1415926535897932385
N[Phi,20]
-> 1.6180339887498948482
N[fft[{1,2,3,4}],20]
N[arg[-1],20]
-> 3.1415926535897932385 Rad

N[Pi*10^20,20]
-> 314159265358979323850
precision[N[Pi*10^20,20]]
-> 19
accuracy[N[Pi/10^20,20]]
-> 39
```

When an exact Rational has a terminating decimal representation, unnecessary trailing zeros are not displayed; for example `N[1/2,10] -> 0.5`. For a certified-interval result produced at a requested significant precision, only a run of redundant trailing zeros is compacted, with one trailing zero retained: a certified `1.000000000000` is displayed as `1.0`, while `1.500000000000` is displayed as `1.50`. The requested digit count, CertifiedEnclosure, and InformationEnclosure remain intact in metadata, so the number of visible zeros is not itself the precision guarantee.

## 7.1 CertifiedEnclosure / InformationEnclosure

`DecimalApproximation` and `ComplexDecimalApproximation` retain two different intervals for every approximate component.

- **CertifiedEnclosure** — an interval that the backend has proved contains the true value. Internal guard digits may make it much narrower than the precision declared to the user. Truth-containment checks and unique decimal-rounding checks use this interval.
- **InformationEnclosure** — an interval describing how much information later computation is allowed to reuse from the approximation. For nonzero `N[x,p]`, if `e=floor(log10(|d|))` for displayed value `d`, it contains both the CertifiedEnclosure and at least `d ± 0.5*10^(e-p+1)`. The information contract therefore follows the value scale and prevents hidden guard digits from later reappearing as user-visible Accuracy. For zero-centered approximations, the InformationEnclosure directly expresses absolute Accuracy instead of relative Precision.

The invariant is always

```text
CertifiedEnclosure ⊆ InformationEnclosure
```

`Infinity` is positive extended-real infinity. `ComplexInfinity` records infinite magnitude without a determined real or complex direction, while `Indeterminate` records that no unambiguous numerical value exists. These are protected atoms rather than finite algebraic symbols. The evaluator currently canonicalizes the following proof-safe exceptional forms:

| Form | Result |
| --- | --- |
| `0/0`, `Infinity/Infinity`, `Infinity-Infinity`, `0*Infinity` | `Indeterminate` |
| `0^0`, `1^Infinity`, `Infinity^0`, `(-1)^Infinity`, `0^I` | `Indeterminate` |
| `z^Infinity` when `abs[z]==1` is proved | `Indeterminate` |
| `a/0` when `a` is proved nonzero | `ComplexInfinity` |
| `x^(1/0)` for any `x` | `Indeterminate` |

For an exact numeric base, unit magnitude is checked with exact Rational real and imaginary parts. For a symbolic base it is used only under an explicit proof such as `simplify[z^Infinity,abs[z]==1]`; a numerical approximation near the unit circle is never guessed to be exactly on it. `Indeterminate` propagates through arithmetic and registered scalar mathematical functions, `N[Indeterminate,p]` preserves it, and `Indeterminate==Indeterminate` is `False`. This is a deliberately bounded exceptional-value contract, not yet a complete algebra of all directed infinities.

The InformationEnclosure is not a probability distribution or a statistical confidence interval, nor does it claim that the backend considers every point in the wider interval mathematically possible. Truth certification belongs exclusively to the CertifiedEnclosure. The InformationEnclosure is an **information contract**: a narrower internal certificate alone does not grant later code permission to recover undeclared digits.

Whenever a later operation must make a discrete semantic decision—proving an approximate value exactly zero/nonzero, selecting a branch side, excluding a pole, or certifying a matrix pivot/rank—it uses the **InformationEnclosure alone**. A CertifiedEnclosure that happens to collapse to an internal point does not authorize such a decision when the InformationEnclosure still crosses the boundary. Consequently `N[0,5]^0` and `1/N[0,5]` are not silently reinterpreted as exact `0^0` or `1/0`.

Ordinary `+ - * /` and unary `-` propagate both intervals independently. Exact `Number` operands enter both paths as identical point intervals.

For complex approximations, mmCal does not compute whole-value quality by taking the minimum of two component-wise relative precisions. Precision/Accuracy are derived from the whole complex InformationEnclosure and value magnitude, and a component proved exactly zero contributes the point interval `{0,0}` to that effective enclosure. An irrelevant exact zero component therefore cannot collapse the precision of a pure-imaginary value.

```text
precision[N[I,20]]
-> 19

precision[N[I/10^100,20]]
-> 19
accuracy[N[I/10^100,20]]
-> 119
```

Exact identities `+0`, `-0`, `*1`, `/1`, and unary `-` do not consume information, so they preserve or mirror the existing Certified/Information metadata instead of re-quantizing the approximation. Exact powers-of-ten scaling likewise does not lose an extra relative digit. Zero-centered results select their display quantum from absolute Accuracy carried by the InformationEnclosure rather than from relative Precision.

```text
precision[N[Pi,20]*10]
-> 19
accuracy[N[Pi,20]*10]
-> 18

precision[sin[N[Pi,20]]]
-> 0
accuracy[sin[N[Pi,20]]]
-> 19
```

If an operation proves `q` significant digits from its InformationEnclosure, formatting may retain up to `q+1` significant display digits. This does not invent information; it prevents the decimal output quantum from being made one decade coarser than the already-guaranteed half-quantum.

`DecimalApproximation` / `ComplexDecimalApproximation` are first-class leaves for the certified scalar evaluator. Functions with interval backends, including `sin`, `exp`, `log`, `sqrt`, hyperbolic/inverse functions, `gamma`, `erf`, `Ei`, `Si`, and `Ci`, plus supported regions of Fresnel C/S, `1F1`, `2F1`, `zeta`, and `polylog`, propagate both enclosures independently on `ComplexInterval`. Functions such as `log2`, `log10`, and `fract` that rewrite to supported primitives re-enter the same path after rewriting. Ordered comparisons and discrete selectors such as `min` / `max` are resolved only when the **InformationEnclosure alone** proves the result, preventing hidden guard digits from leaking through Boolean decisions. Backends that currently require exact Rational parameters, including parts of `1F1` / `2F1`, elliptic functions, and `polylog`, remain conservatively unevaluated for unsupported approximate parameters.

```text
N[Pi,20] + 1/3

Certified:    C(Pi) + {1/3}
Information:  I(Pi) + {1/3}
```

The output decimal is justified from the CertifiedEnclosure, while the number of digits the result may declare is capped by the InformationEnclosure. Scaling and cancellation therefore reduce `accuracy` / `precision` naturally. The propagated InformationEnclosure is stored in the result itself rather than being compressed back to a digit count after every operation.

An outer `N` cannot narrow an existing InformationEnclosure merely by increasing working precision. Thus

```text
N[N[Pi,20],100]
-> 3.1415926535897932385
```

retains the original 20-digit guarantee. Asking for fewer digits is allowed to discard information by adding the coarser output-rounding interval. Both enclosures are propagated through `RealInterval` / `ComplexInterval` with outward rounding; no `double` or machine-real fallback is used.

---

# 8. precision / accuracy / rationalize

## 8.1 `accuracy[x]`

For a `DecimalApproximation`, returns an **integer lower bound on the guaranteed number of absolute decimal digits** relative to the true value.

```text
accuracy[N[1/3,20]]
-> 20
```

Given displayed value `d` and InformationEnclosure `[iL,iU]`, the implementation uses

```text
max(|d-iL|, |d-iU|)
```

as its absolute error bound. The InformationEnclosure created by `N[...,p]` already contains the scale-dependent half-quantum implied by significant-digit rounding, so even a terminating decimal whose CertifiedEnclosure is an exact point cannot recover hidden guard information as additional Accuracy.

Exact numbers and exact symbolic expressions return `Infinity`.

```text
accuracy[1/3] -> Infinity
accuracy[Pi]  -> Infinity
```

## 8.2 `precision[x]`

Uses the same absolute error bound divided by a positive lower bound on the value magnitude derived from the InformationEnclosure, returning an **integer lower bound on the guaranteed number of relative decimal digits**.

```text
precision[N[1/3,20]]
-> 19
```

This function does not mechanically return the requested 20 significant digits. Around `1/3`, 20-significant-digit rounding has quantum `10^-20`; the resulting relative uncertainty permits an integer lower bound of 19 guaranteed relative decimal digits. If the InformationEnclosure contains zero, no positive lower bound on the value magnitude is available, so the result is 0. Exact expressions return `Infinity`. Near-cancellation can therefore preserve substantial absolute Accuracy while losing many relative Precision digits.

## 8.3 `rationalize[x]`

Finds the **exact Rational with the smallest denominator** contained in the InformationEnclosure of an approximate value. The search uses continued-fraction-style interval recursion over exact Rational values and never converts the interval to `double`. Using the InformationEnclosure prevents `rationalize` from recovering hidden guard digits or a hidden exact point that was not declared by the approximation.

```text
rationalize[N[1/3,20]]
-> 1/3
```

`rationalize[x,tol]` chooses the smallest-denominator Rational within `[x-tol,x+tol]` centered on the displayed value. `tol` must be a non-negative exact real.

```text
rationalize[N[Pi,20],1/1000]
-> 201/64
```

`201/64` lies within 0.001 of `Pi` and has a smaller denominator than `355/113`, so it is the correct result under this specification.

With `tol=0`, the displayed finite decimal itself is converted back to an exact Rational.

```text
rationalize[N[1/3,20],0]
-> 33333333333333333333/100000000000000000000
```

Because the source literal `0.1` is already parsed as exact `1/10` in mmCal, `rationalize[0.1]` simply remains `1/10`. `DecimalApproximation` values nested inside Arrays or expressions are also rationalized recursively.

## 8.4 `explain[value]`

A lightweight introspection function that returns only information **already carried by the evaluated value**. It does not run `det`, mathematical `matrixRank`, LU, Eigen, or other derived computations, and it does not scan a Generic Array merely to infer additional properties.

```text
explain[{{1,2},{3,4}}]
-> {{"Kind","Array"},
    {"Domain","Integer"},
    {"Exactness","Exact"},
    {"ArrayRank",2},
    {"Dimensions",{2,2}},
    {"ElementCount",4},
    {"Rectangular",True},
    {"Empty",False},
    {"Matrix",True},
    {"Square",True},
    {"Order",2}}
```

Arguments are evaluated normally before introspection, so `explain[1+2]` describes `3`. History outputs can be inspected directly with `explain[Out[n]]` or `explain[%]`.

The result is displayed as a brace sequence of property/value pairs. Because values such as `Dimensions` and the enclosure properties may themselves be Arrays or brace values, the internal result is a general `ListExpr` rather than a dense `ArrayExpr`; the structured values are not flattened into strings.

`Exactness` is a classification rather than a boolean. Current principal values are `"Exact"`, `"CertifiedApproximation"`, and `"Unknown"`.

Built-in mathematical constants and predefined symbols are not collapsed into ordinary unknown symbols. `Pi/E/Phi` use the mathematical metadata already registered in MathRegistry, while values such as `Infinity` use their predefined SymbolRegistry semantics; both are O(1) lookups.

```text
explain[Pi]
-> {{"Kind","Constant"},
    {"Domain","Real"},
    {"Exactness","Exact"},
    {"Name","Pi"},
    {"Real",True},
    {"Positive",True},
    {"Irrational",True},
    {"ArithmeticClass","Transcendental"}}

explain[Infinity]
-> {{"Kind","Constant"},
    {"Domain","ExtendedReal"},
    {"Exactness","Exact"},
    {"Name","Infinity"},
    {"Infinite",True},
    {"Finite",False},
    {"Sign","Positive"}}
```

`I` is lowered by normal evaluation to an exact complex `Number`, so `explain[I]` describes the evaluated complex value rather than the input token.

Builtin function symbols are also described directly from existing BuiltinRegistry / MathRegistry metadata in O(1) time. This can expose arity, held-argument rules, and mathematical metadata such as domain, parity, period, principal inverse, and real range without executing the function.

```text
explain[sin]
-> {{"Kind","BuiltinFunction"},
    {"Domain","Function"},
    {"Exactness","Exact"},
    {"Name","sin"},
    {"Arity",1},
    {"ArgumentEvaluation","All"},
    {"FunctionDomain","ComplexToComplexRealPreserving"},
    {"Parity","Odd"},
    {"PeriodTurns",1},
    {"PrincipalInverse","asin"}, ...}

explain[table]
-> ... {"ArgumentEvaluation","HoldFirstAndTableIteratorSpec"} ...
```

This is registry introspection, not an attempt to execute the function to discover additional properties.

Certified decimal approximations expose both the exact Rational CertifiedEnclosure and the InformationEnclosure in addition to their requested significant-digit precision. The former is the truth certificate; the latter is the information limit propagated by later approximate arithmetic. Arrays expose O(1) storage-domain/exactness metadata plus essentially free shape properties such as `ArrayRank`, `Dimensions`, `ElementCount`, `Vector`, `Matrix`, `Square`, `Order`, and `Empty`. Determinant, mathematical rank, invertibility, eigenvalues, and similar derived properties are intentionally omitted.

Integers expose sign, zero, and `BitLength`. Decimal digit count is not computed automatically because huge integers would require decimal conversion; Rationals instead expose numerator/denominator bit lengths.

```text
explain[value,"internal"]
```

adds development/performance diagnostics such as `Representation`, Array `Storage` / `Contiguous` / `StoredExpressions`, and approximation `ApproximationOrigin`. **`"internal"` property names and values are not a compatibility-stable API.** Unknown modes are errors; no hidden expensive `"full"` mode is executed.

---

# 9. Basic arithmetic and algebra

Standard operators:

```text
+  -  *  /  ^  !
```

`pow[x,y]` is a source-level alias for `Power[x,y]`, and `fact[x]` is an alias for factorial.

Representative exact simplifications:

```text
sqrt[8]
-> 2 sqrt[2]

sqrt[2/3]
-> sqrt[6] / 3

sqrt[-8]
-> 2I sqrt[2]

sqrt[z]^2
-> z
```

Transformations that would violate principal-branch semantics are not applied.

```text
sqrt[x^2]
-> sqrt[x ^ 2]       // x is unconstrained

simplify[sqrt[x^2], element[x,Real]]
-> abs[x]

simplify[sqrt[x^2], x >= 0]
-> x
```

---

# 10. Basic mathematical and complex functions

| Function      | Description                        | Example                          |
| ------------- | ---------------------------------- | -------------------------------- |
| `sqrt[x]`     | Principal square root              | `sqrt[-4] -> 2I`                 |
| `cbrt[x]`     | Real cube root, real domain        | `cbrt[-8] -> -2`                 |
| `abs[z]`      | Absolute value / complex magnitude | `abs[3+4I] -> 5`                 |
| `sign[z]`     | Real sign / complex `z/abs[z]`     | `sign[3+4I] -> 3/5+4/5I`         |
| `re[z]`       | Real part                          | `re[3+4I] -> 3`                  |
| `im[z]`       | Imaginary part                     | `im[3+4I] -> 4`                  |
| `conj[z]`     | Complex conjugate                  | `conj[3+4I] -> 3-4I`             |
| `arg[z]`      | Principal argument                 | `arg[-1] -> Pi Rad`              |
| `hypot[x,y]`  | Exact `sqrt[x^2+y^2]`              | `hypot[3,4] -> 5`                |
| `cis[x]`      | `cos[x]+I sin[x]`                  | `cis[Pi/3] -> 1/2 + I sqrt[3]/2` |
| `polar[r,t]`  | `r cis[t]`                         | `polar[2,Pi/3]`                  |
| `nextpow2[x]` | Smallest `n` such that `2^n >= x`  | `nextpow2[9] -> 4`               |

Compatibility aliases:

```text
real -> re
imag -> im
mag  -> abs
unit,csgn -> sign
rect -> polar
```

---

# 11. Exponential and logarithmic functions

| Function   | Semantics                                 |
| ---------- | ----------------------------------------- |
| `exp[x]`   | Entire complex exponential                |
| `log[x]`   | Principal natural logarithm               |
| `log[b,x]` | Principal `Log[x]/Log[b]`                 |
| `log2[x]`  | `log[2,x]` frontend                       |
| `log10[x]` | `log[10,x]` frontend                      |
| `expm1[x]` | Stable evaluation of `exp[x]-1` near zero |
| `log1p[x]` | Stable evaluation of `log[1+x]` near zero |

`expm1` and `log1p` do not form `exp[x]-1` or `1+x` only at the requested precision when `x` is tiny. Extra working precision is derived from the binary scale of the input, and the result is rounded only after the cancellation-sensitive operation.

Examples:

```text
exp[1]
-> E

log[E]
-> 1

log[-1]
-> I Pi

log[10,1000]
-> 3

log[1,10]
-> DomainError       // base 1
```

`log[0]` produces a DomainError rather than being replaced by Infinity.

---

# 12. Trigonometric functions

Implemented:

```text
sin cos tan cot sec csc
asin acos atan atan2
```

With radians as the default:

```text
sin[Pi/6] -> 1/2
cos[Pi/3] -> 1/2
tan[Pi/4] -> 1
asin[1/2] -> Pi/6
atan2[1,-1] -> 3Pi/4
```

Poles of `tan`, `sec`, `cot`, and `csc` are treated as definedness conditions. Finite values are not fabricated at poles.

---

# 13. Hyperbolic functions

Implemented:

```text
sinh cosh tanh
asinh acosh atanh
csch sech coth
```

Inverse hyperbolic functions with principal complex branches carry branch metadata in `MathRegistry`.

---

# 14. Cardinal / stable elementary functions

```text
sinc[x]
cosc[x]
tanc[x]
sinhc[x]
tanhc[x]
expc[x]
```

Removable singularities are filled exactly.

```text
sinc[0]  -> 1
cosc[0]  -> 0
tanc[0]  -> 1
sinhc[0] -> 1
tanhc[0] -> 1
expc[0]  -> 1
```

Trigonometric cardinal functions normalize angle expressions to radian quantities before forming the ratio. If a finite-precision InformationEnclosure crosses zero, `sinc/cosc/sinhc/expc` use a local Taylor polynomial with an explicit remainder bound and `tanc/tanhc` use safe continuous-extension forms, so a removable singularity does not by itself force the result to remain unevaluated.

```text
sinc[Pi/2]
-> 2 / Pi

sinc[90 Deg]
-> 2 / Pi
```

---

# 15. Rounding and integer utilities

```text
floor ceil trunc round frac
gcd lcm mod rem quotient
bitand bitor bitxor bitnot
bitshiftl bitshiftr bitlength bitcount bitget
fma clamp proj
```

Examples:

```text
floor[-3/2] -> -2
ceil[-3/2]  -> -1
trunc[-3/2] -> -1
round[5/2]  -> 2       // nearest-even
round[125,-1] -> 120
round[135,-1] -> 140
frac[-3/2]  -> 1/2

gcd[84,126,210] -> 42
lcm[6,8,9] -> 72

quotient[-5,3] -> -1
rem[-5,3]      -> -2
mod[-5,3]      -> 1

bitand[-1,5] -> 5
bitor[-8,3]  -> -5
bitxor[-1,5] -> -6
bitnot[5]    -> -6
bitshiftr[-3,1] -> -2
bitget[-2,100]  -> 1
```

`round[x,n]` performs nearest-even rounding to a quantum of `10^-n`; negative `n` is allowed. Approximate inputs resolve only when the whole InformationEnclosure rounds to the same value.

The bitwise functions use **infinite two's-complement semantics** over arbitrary BigInt values. Negative operands therefore sign-extend with ones. `bitshiftr` is the user-facing arithmetic right shift, while `bitcount` accepts only nonnegative integers because negative infinite two's-complement values have infinitely many one bits. This version deliberately keeps the semantics in function form rather than adding lexer-level `& | << >>` syntax.

`fma[a,b,c]` remains exact for exact operands and, for certified approximations, evaluates `a*b+c` without an intermediate DecimalApproximation rounding. `clamp[x,lo,hi]` is real-valued and uses InformationEnclosure to prove approximate ordering. `proj[z]` is currently the identity on finite exact/certified complex values; full complex-infinity/Riemann-sphere semantics remain tied to the deferred extended-real model.

`mod` corresponds to a floor quotient, while `rem` corresponds to a truncate-toward-zero quotient.

---

# 16. Combinatorics and lightweight number theory

```text
perm[n,r]
comb[n,r]
fib[n]
```

Examples:

```text
perm[10,3] -> 720
comb[10,3] -> 120
fib[100]   -> 354224848179261915075
```

`fib` uses fast doubling in O(log n).
`comb` uses the symmetry `r=min[r,n-r]`.

The lightweight exact number-theory set also includes:

```text
isprime[n]
nextprime[n]
prevprime[n]
factorint[n]
totient[n]
```

```text
isprime[97]    -> True
nextprime[14]  -> 17
prevprime[14]  -> 13
factorint[360] -> {2, 2, 2, 3, 3, 5}
totient[9]     -> 6
```

`isprime` is deterministic on `0 <= n <= 2^64-1` using strong Miller-Rabin bases whose proven range covers all `uint64` values. `factorint` / `totient` use deterministic Pollard-Rho together with exact prime verification in the same `uint64` domain. Values beyond the current proof backend remain unevaluated rather than returning a probable-prime result as `True`. `factorint[-n]` prefixes `-1` to the flat prime-factor list; `factorint[0]` is a DomainError.

---

# 17. Special functions

## 17.1 Gamma / LogGamma

```text
gamma[5]
-> 24

gamma[1/2]
-> sqrt[Pi]

gamma[-1/2]
-> -2 sqrt[Pi]

N[gamma[1/3],20]
-> 2.6789385347077476337
```

General real arguments use Stirling–Bernoulli with a rigorous remainder bound; negative real arguments use reflection.
Non-positive integer poles produce DomainError. For exact complex input, `N` applies recurrence shifts into a suitable right half-plane and evaluates a `ComplexInterval` Stirling expansion with an explicit remainder bound.

```text
N[gamma[1+I],20]
-> 0.49801566811835604271-0.15494982830181068512I
```

`lgamma[x]` currently means **`log[abs[gamma[x]]]` on the real axis**. It is kept separate from complex `LogGamma`.

## 17.2 Erf

```text
erf[x]
erfc[x]
```

```text
erf[0]  -> 0
erfc[0] -> 1
N[erf[1],20] -> 0.84270079294971486934
N[erf[1+I],20]
-> 1.3161512816979476449+0.19045346923783468628I
```

Complex input is evaluated by the entire power series on `ComplexInterval` with an explicit tail bound.

## 17.3 Beta

Currently restricted to positive real arguments.

```text
beta[2,3] -> 1/12
beta[1/2,1/2] -> Pi
betaln[1/2,1/2] -> log[Pi]
```

The implementation does not unconditionally expand to a general Gamma ratio when doing so could break pole cancellation.

## 17.4 Zeta / Digamma / Trigamma / regularized incomplete Beta

```text
zeta[s]
digamma[x]
trigamma[x]
ibeta[a,b,x]
```

`zeta` denotes the Riemann zeta function. Representative exact values and trivial zeros are simplified exactly. The certified `N` backend covers the finite complex plane except `s=1`: for `Re[s]>=0` it uses Euler-Maclaurin summation with the remainder condition checked at each correction order, while the left half-plane is mapped through the functional equation. Only `s=1` is a pole; a finite-precision enclosure that may contain the pole yields `N::precision`.

```text
zeta[0]  -> -1/2
zeta[-2] -> 0
zeta[2]  -> Pi^2/6
N[zeta[3],20] -> 1.2020569031595942854
N[zeta[2+I],20]
-> 1.1503557032549026717-0.43753086591960788112I
N[zeta[1/2],20] -> -1.4603545088095868129
N[zeta[1/2+I],20]
-> 0.14393642707718906032-0.72209974353167308913I
```

`digamma[x]` is the derivative of `lgamma[x]`; `trigamma[x]` is the derivative of `digamma[x]`. Non-positive integer poles are DomainErrors. The certified backend covers positive real inputs and general complex inputs by recurrence into the right half-plane followed by Bernoulli asymptotics on `ComplexInterval` with an explicit remainder bound. Intervals containing poles and inputs beyond the bounded planner are distinguished as domain/backend failures. Positive-integer trigamma values reduce exactly to `Pi^2/6` minus a finite second-order harmonic sum.

```text
N[digamma[1],20]  -> -0.57721566490153286061
N[trigamma[1],20] -> 1.6449340668482264365
N[digamma[1+I],30]  -> 0.0946503206224769772718784827219+1.07667404746858117413405079475I
N[trigamma[1+I],30] -> 0.463000096622763786298326518184-0.794233542759318865583013617157I
trigamma[2]        -> Pi^2/6-1
D[gamma[x],x]      -> digamma[x]gamma[x]
D[lgamma[x],x]     -> digamma[x]
D[digamma[x],x]    -> trigamma[x]
```

`ibeta[a,b,x]` is the **regularized incomplete beta function** `I_x(a,b)`. Its real contract is `a>0`, `b>0`, `0<=x<=1`. Positive-integer `a,b` with exact Rational `x` reduce to a finite binomial sum. Certified `N` also propagates finite-precision positive-real `a,b` and real `x` as intervals; monotonicity of `I_x(a,b)` in `a`, `b`, and `x` encloses the whole parameter box from certified endpoint evaluations.

```text
ibeta[1,1,1/4] -> 1/4
ibeta[2,3,1/2] -> 11/16
N[ibeta[1/3,2/3,1/4],20] -> 0.53302858123542523627
N[ibeta[N[1/3,8],N[2/3,8],1/4],8] -> 0.53302858
```

Higher polygamma and complex-parameter `ibeta` remain intentionally deferred.

## 17.5 Generalized factorial family

```text
binom[x,n]
fallingfact[x,n]
risingfact[x,n]
```

At present, these functions primarily construct exact finite products when the order is a non-negative integer.

```text
binom[1/2,2] -> -1/8
fallingfact[5,3] -> 60
risingfact[5,3] -> 210
```

## 17.6 Fresnel C / S

mmCal uses the standard Fresnel integrals corresponding to

```text
fresnelc[x] = integral_0^x cos[Pi t^2/2] dt
fresnels[x] = integral_0^x sin[Pi t^2/2] dt
```

as entire odd functions. Arguments that do not close exactly remain symbolic, while `N` certifies both real and complex inputs on `ComplexInterval`. The complex path has no fixed `|z|` boundary: moderate and diagonal-sector arguments use the entire Maclaurin series, while large near-axis arguments use the DLMF 7.12 `f/g` asymptotic expansions with first-neglected-term remainder bounds. Quarter-turn identities `C[i z]=i C[z]` and `S[i z]=-i S[z]` map all coordinate-axis neighborhoods into the same certified wedge; if the asymptotic proof cannot close, evaluation falls back to the Maclaurin path. Work is bounded by term caps and the shared `EvaluationBudget`. The real path likewise has a certified large-argument asymptotic backend.

```text
fresnelc[0] -> 0
fresnels[0] -> 0
N[fresnelc[1],20] -> 0.77989340037682282947
N[fresnels[1],20] -> 0.43825914739035476608
N[fresnelc[1+I],20] -> 2.5557937781024390246+2.5557937781024390246I
N[fresnels[1+I],20] -> -2.0618882191948404681+2.0618882191948404681I

D[fresnelc[x],x] -> cos[Pi x^2/2 Rad]
D[fresnels[x],x] -> sin[Pi x^2/2 Rad]
```

`Rad` is explicit in the derivatives because the Fresnel definitions themselves must not depend on the session's default angle unit.

## 17.7 Confluent hypergeometric 1F1

Kummer's confluent hypergeometric function is written as

```text
hypergeometric1F1[a,b,z]
```

It is entire in `z`; in general `b = 0,-1,-2,...` is a parameter pole, so those cases are not unconditionally simplified to finite values. Exact evaluation currently handles safely terminating series, `z=0`, `a=b`, and related closed cases. The certified `N` backend accepts exact Rational parameters `a,b` and real or complex `z`, evaluating the entire series on `ComplexInterval` with an explicit tail bound. Approximate parameters themselves remain unsupported.

```text
hypergeometric1F1[0,3,2] -> 1
hypergeometric1F1[-2,3,2] -> 0
hypergeometric1F1[2,2,1] -> E
N[hypergeometric1F1[1/6,7/6,1],20]
-> 1.1920688079818883008
N[hypergeometric1F1[1/2,5/4,1+I],20]
-> 1.2988692086674201067+0.72862412303674683434I
```

When the parameters do not depend on the differentiation variable,

```text
D[hypergeometric1F1[a,b,z],z]
= a hypergeometric1F1[a+1,b+1,z]/b
```

is used. For integration, mmCal prefers the 1F1 form when an upper-incomplete-Gamma representation would introduce principal-branch structure or a removable hole at the origin. For example,

```text
integrate[exp[-x^2],x]
-> erf[x]sqrt[Pi]/2

integrate[exp[-x^2],{x,0,Infinity}]
-> sqrt[Pi]/2

integrate[exp[-x^2],{x,-Infinity,Infinity}]
-> sqrt[Pi]

integrate[exp[x^6],x]
-> x hypergeometric1F1[1/6, 7/6, x^6]
```

The same family handles `exp[c x^n]` for positive integer `n`.

## 17.8 Gauss hypergeometric 2F1

The Gauss hypergeometric function is written as

```text
hypergeometric2F1[a,b,c,z]
```

In general `c = 0,-1,-2,...` is a parameter pole, and mmCal uses the principal branch in `z`. Exact evaluation currently handles terminating series generated by non-positive-integer numerator parameters and other safe degenerations such as `a=0` or `b=0`. The certified `N` backend lifts supported exact numeric parameters to `ComplexInterval`. For `|z|<1` it uses the Gauss series with a rigorous tail bound. At `z=1`, Gauss summation `Gamma[c] Gamma[c-a-b]/(Gamma[c-a] Gamma[c-b])` is used whenever `Re(c-a-b)>0` is certified. For provable `|z|>1`, it may use the principal `1/z` connection formula only when parameter degeneracies and branch-cut hazards are provably excluded and `|1/z|<1` is certified. An exact real `z>1` is evaluated using the defined principal-cut continuation value. By contrast, a finite-precision input such as `2+I*N[0,p]` that still permits both sides of the cut returns `N::precision` instead of selecting one side. Other unproved boundary or degenerate cases remain unevaluated.

```text
hypergeometric2F1[-2,1,3,1/2] -> 17/24
hypergeometric2F1[0,2,3,x] -> 1
N[hypergeometric2F1[1/2,1/2,3/2,1/4],20]
-> 1.0471975511965977462
N[hypergeometric2F1[1/2,1/3,5/4,1/2+I/4],20]
-> 1.0768624682230010329+0.057258816434281232216I
N[hypergeometric2F1[3.4,5.6,4+I,4.6+2I],20]
-> 0.0046136876612922014955+0.0019659119401108294965I
```

When the parameters do not depend on the differentiation variable,

```text
D[hypergeometric2F1[a,b,c,z],z]
= a b hypergeometric2F1[a+1,b+1,c+1,z]/c
```

is used. Integration uses this family for binomial powers, for example

```text
integrate[sqrt[1+2x^3],x]
-> x hypergeometric2F1[-1/2, 1/3, 4/3, -2x^3]

integrate[1/(1+x^5),x]
-> x hypergeometric2F1[1, 1/5, 6/5, -x^5]
```

There is intentionally no general `Solve` inversion rule for 2F1 because global injectivity is not available in general. Only exact degenerations that reduce to existing algebraic expressions are passed on to the ordinary solver.

## 17.9 Incomplete elliptic integrals F / E / Pi

mmCal writes the Legendre incomplete elliptic integrals as

```text
ellipticF[phi,m]
ellipticE[phi,m]
ellipticPi[n,phi,m]
```

The second argument `m` is the parameter, and the amplitude `phi` is **always interpreted in radians**, independent of the session's `Deg/Rad/Grad` angle mode. Principal branches are used; the branch cuts and poles of general complex parameters are not collapsed into an unsafe "defined everywhere" rule.

```text
ellipticF[x,0] -> x
ellipticE[x,0] -> x
ellipticPi[0,x,0] -> x

N[ellipticF[1/2,1/3],20]
-> 0.5068477562654311092
N[ellipticE[1/2,1/3],20]
-> 0.49331536201475850521
N[ellipticPi[1/5,1/2,1/3],20]
-> 0.51520338216141386085
N[ellipticF[1/2,99/100],20]
-> 0.52198775871658283077
N[ellipticE[Pi/2,1],20]
-> 1.0
N[ellipticF[1/2,2],20]
-> 0.55135887907967981413
```

When the certified effective tail ratio is sufficiently small, the real backend keeps the guarded Legendre series as a fast path. Outside that region it reduces real amplitudes modulo `Pi` and evaluates the Legendre forms through certified Carlson symmetric integrals `RF`, `RD`, and `RJ`. There is no longer a fixed `|m|<=9/10` or `|n|<=9/10` capability boundary. Local real values with `m>1` or `n>1` are accepted when the whole reduced integration path can be proved to stay before the corresponding branch point or pole; period-crossing values require those singularities to be excluded. Exact `m=1` for `ellipticE` uses its finite real degeneration. General complex elliptic continuation is still unsupported. Amplitude derivatives are

```text
D[ellipticF[phi,m],phi]
= 1/sqrt[1-m sin[phi Rad]^2]

D[ellipticE[phi,m],phi]
= sqrt[1-m sin[phi Rad]^2]

D[ellipticPi[n,phi,m],phi]
= 1/((1-n sin[phi Rad]^2)sqrt[1-m sin[phi Rad]^2])
```

so the standard kernels integrate directly.

```text
integrate[1/sqrt[1-(1/3)sin[x]^2],x]
-> ellipticF[x, 1/3]

integrate[sqrt[1-(1/3)sin[x]^2],x]
-> ellipticE[x, 1/3]

integrate[1/((1-(1/5)sin[x]^2)sqrt[1-(1/3)sin[x]^2]),x]
-> ellipticPi[1/5, x, 1/3]

integrate[1/sqrt[1-x^4],x]
-> ellipticF[asin[x], -1]
```

The final quartic reduction is a correct local primitive, but the current `fullSimplify` cannot always prove the corresponding `sin[asin[x]]` and principal-square-root product identity globally. The derivative-back harness therefore monitors it in ResolutionOnly mode instead of reducing integration capability because of a proof-engine limitation. General elliptic equations also remain unresolved by `Solve` until a principled inverse-elliptic function family exists; exact degenerations such as `m=0` are solved by the existing solver.

## 17.10 Ei / Si / Ci / li / Polylogarithm

The principal special functions commonly required by symbolic integration are exposed as

```text
Ei[x]
Si[x]
Ci[x]
li[x]
polylog[s,z]
```

`Ei`, `Ci`, `li`, and `polylog` generally have branch structure and are registered as principal-branch functions. `Si` is entire and odd. `Ei/Si/Ci` certify supported complex series regions directly on `ComplexInterval`. `polylog[n,z]` supports positive integer order in the provable `|z|<1` region through a series plus tail bound; for `n=2`, the principal DLMF 25.12.3/25.12.4/25.12.6 connection formulas are also used where the branch cut can be excluded. Near `z=1` on the positive real axis, orders 3 through 12 use a certified positive-integer-limit `mu=log(z)` expansion as a fast path, including finite-precision real intervals without collapsing them to hidden point values. Branch cuts and convergence boundaries are not filled with heuristic values.

```text
Si[0] -> 0
Si[-1] -> -Si[1]
li[0] -> 0
polylog[0,z] -> z/(1-z)
polylog[1,z] -> -log[1-z]
polylog[2,1] -> Pi^2/6
polylog[2,-1] -> -Pi^2/12
polylog[3,1] -> zeta[3]
polylog[3,-1] -> -3zeta[3]/4

N[Ei[1],20] -> 1.8951178163559367555
N[Si[1],20] -> 0.94608307036718301494
N[Ci[1],20] -> 0.33740392290096813466
N[li[2],20] -> 1.0451637801174927848
N[polylog[2,1/2],20] -> 0.58224052646501250590
N[polylog[2,999/1000],20] -> 1.6370226052761177427
N[polylog[3,999/1000],20] -> 1.2004153539954643452
N[polylog[2,-2],20] -> -1.4367463668836809464
N[Ei[1+I],20] -> 1.7646259855638540684+2.3877698515105224193I
N[Si[1+I],20] -> 1.1042226582355817396+0.88245380500791774338I
N[Ci[1+I],20] -> 0.88217218055593632505+0.28724913351995593953I
N[polylog[2,1/2+I/4],20]
-> 0.54586750496407962676+0.33913769923976904082I
```

Certified special-function backends also have bounded-work implementation limits that are separate from their mathematical domains. `1F1` no longer has a fixed `|z|` threshold: both real and complex series are attempted while a rigorous future-term ratio can enter the convergent regime within the 250000-term series budget. The `2F1` Gauss series uses its mathematical convergence region `|z|<1` directly, with no fixed interior threshold and with term/global budgets providing bounded work; `z=1` is additionally supported by Gauss summation when `Re(c-a-b)>0`. Real elliptic `F/E/Pi` likewise no longer use the former `9/10` parameter work boundary: the series remains a fast path when the effective certified tail ratio is `<=9/10`, while a Carlson `RF/RD/RJ` backend covers additional real regions subject to explicit branch/pole proofs, duplication refinement limits, and the shared `EvaluationBudget`. Complex `Ei/Ci` likewise have no fixed magnitude boundary. Moderate arguments use guarded interval series, while large arguments use the DLMF 6.12 `E1` asymptotic expansion plus principal connection formulas; the negative-real-axis cut is decided from the InformationEnclosure. If the asymptotic remainder proof does not close, evaluation falls back to the series, with the shared `EvaluationBudget` bounding work. The former positive-order `polylog` cutoff at `|z|<=49/50` has been removed: the `|z|<1` series is bounded by its 1000000-term cap and the shared `EvaluationBudget`. `Li_2` extends to the negative real axis, near the unit circle, and selected `|z|>1` regions away from the principal cut through certified connection formulas; exact positive-real points on the cut remain held rather than choosing an upper or lower boundary value. Higher positive integer orders use the `mu=log(z)` integer-limit expansion as a fast path near positive-real `z=1`. Real `Ei/Si/Ci` no longer have a fixed magnitude boundary: guarded Taylor series handle moderate arguments, certified asymptotic expansions handle large arguments, and remainder proofs, term caps, plus the shared `EvaluationBudget` provide bounded work. Complex `2F1` additionally uses the principal `1/z` connection only when `|z|>1` is certified and degenerate parameters and branch-cut hazards can be excluded.

A value existing mathematically outside one of these thresholds does not imply that the current backend can certify it. Fixed series ranges, fixed term caps, and fixed planner budgets return `N::unsupported` with the original expression instead of repeatedly increasing precision without changing the applicable algorithm. For zeta, only `s=1` is a pole; the critical strip and left half-plane are now covered by the certified continuation backend. Guard precision is increased when interval width may be responsible for straddling a backend boundary or branch cut, but top-level `N` is bounded to 16 local refinement attempts. If an existing finite-precision input enclosure keeps the decision ambiguous and the requested digits cannot be certified, `N::precision` preserves the unevaluated expression instead of exhausting the global resource budget.

When the argument and order parameters are independent of the differentiation variable, the derivative knowledge includes

```text
D[Ei[x],x] -> exp[x]/x
D[Si[x],x] -> sinc[x Rad]
D[Ci[x],x] -> cos[x]/x
D[li[x],x] -> 1/log[x]
D[polylog[s,x],x] -> cases[polylog[s-1, x]/x if x != 0; 1 if x == 0]
```

The exact degeneration `polylog[1,x] -> -log[1-x]` gives

```text
D[polylog[2,x],x] -> cases[-log[1-x]/x if x != 0; 1 if x == 0]
```

For direct-variable repeated derivatives `D[polylog[s,x],{x,n}]` with `n<=64`, mmCal uses the Euler operator `theta=x D` and signed Stirling numbers instead of building nested `D[cases[...]]`; at `x=0`, the exact series coefficient gives `n!/n^s`.

The same shared knowledge closes

```text
integrate[exp[x]/x,x] -> Ei[x]
integrate[sin[x]/x,x] -> Si[x]
integrate[cos[x]/x,x] -> Ci[x]
integrate[1/log[x],x] -> li[x]
integrate[li[x],x] -> x li[x]-Ei[2log[x]]
integrate[log[1-x]/x,x] -> -polylog[2, x]
```

No general `Solve` rule invents a single principal inverse for `Ei/Si/Ci/li/polylog`: global injectivity and branch structure are not generally available. Only exact degenerations such as `polylog[0,z]` and `polylog[1,z]` are passed to the existing algebraic/logarithmic Solver.

## 17.11 Lambert W

```text
lambertw[z]
lambertw[k,z]
```

`lambertw` denotes the Lambert W function satisfying `w exp[w] == z`. The one-argument form is the principal branch `k=0`; the two-argument form specifies an integer branch `k`. It remains an exact symbolic function for representative exact values, differentiation, and Real-domain exponential-equation solving, while `N` also provides certified evaluation for the real branches and arbitrary integer complex branches.

```text
lambertw[0] -> 0
lambertw[E] -> 1
lambertw[-1/E] -> -1
lambertw[-1,-1/E] -> -1
D[lambertw[x],x] -> exp[-lambertw[x]]/(1+lambertw[x])
D[lambertw[x],{x,2}] -> (-2-lambertw[x])exp[-2lambertw[x]]/(1+lambertw[x])^3
```

On the real axis the Solver distinguishes the real `k=0` and `k=-1` branches where required and prefers the monotone real inverse backend whenever the branch value is real. The complex certified backend combines a principal-branch Maclaurin/contraction path with the branch-explicit fixed-point form `Log[z]+2 Pi I k-Log[w]`. Any integer branch index is accepted and the proof keeps the requested `k` rather than inferring a nearby branch. Negative-real complex values are also handled when the branch-cut side is exact; only residual regions where the present contraction proof cannot close remain held with `N::unsupported`.

```text
N[lambertw[1],20] -> 0.5671432904097838730
N[lambertw[-1,-1/10],20] -> -3.5771520639572972184
N[lambertw[1+I],20] -> 0.65696606923043640587+0.32545033941341502999I
N[lambertw[2,1],20] -> -2.4015851048680028842+10.776299516115070898I
N[lambertw[-1/E+I/10^8],20] -> -0.99983512787429915685+0.00016485400656056139308I
N[lambertw[-1,-1/E+I/10^8],20] -> -1.0001648721257003724-0.00016489025031827418034I
```

Near `-1/E`, complex evaluation switches to a square-root local coordinate using `u=W+1` and `q=E z+1`. `W_0` follows the principal-square-root side, `W_-1` is the local branch approached from the upper side of the negative real axis, and `W_1` is the symmetric lower-side branch. If finite input information cannot determine the cut side for a nonprincipal branch, mmCal returns `N::precision` instead of choosing a side.

---

# 18. Arrays

## 18.1 Array foundation

At the language level, `{...}` is a general finite brace container rather than a matrix-only literal. When every child has the same shape, the value is automatically promoted to a dense `ArrayExpr`; heterogeneous-shape values such as `{Q,R}` and ragged braces remain general brace values. Numeric dense Arrays may internally use shared packed Integer / Rational / Number pages and strided views, but those storage choices are not user-visible types. Matrix operations still accept only dense rectangular Arrays and audit this at their boundary.

```text
dimensions[A]
arrayRank[A]
length[A]
at[A,i,...]
reshape[A,{d1,d2,...}]
identity[n]
zeros[rows,cols]
rows[A]
cols[A]
diag[A]
trace[A]
```

Indices are zero-based. `at` also accepts a prefix shorter than the Array rank and returns the remaining subarray; only a full-rank index returns a scalar. Finite `SolutionSet` values use the same zero-based convention. `at[solutions,i]` returns a one-branch `SolutionSet`, preserving branch conditions, free variables, multiplicity, and solver-variable domains. `at[solutions,i,x]` returns only the right-hand side bound to `x` in that branch. Because that three-argument form does not carry condition metadata, use the two-argument branch form when conditional information must be preserved. `Conditional`, `Universal`, and `Unresolved` sets are not indexable because they do not define one unambiguous explicit branch sequence.

```text
dimensions[{{1,2,3},{4,5,6}}] -> {2, 3}
arrayRank[{{1,2},{3,4}}] -> 2
at[{{1,2},{3,4}},1] -> {3, 4}
at[{{1,2},{3,4}},1,0] -> 3
at[solve[x^2==1,x],0] -> {x == 1}
at[solve[x^2==1,x],1,x] -> -1
at[solve[{x+y==3,x*y==2},{x,y}],1,y] -> 1
reshape[{1,2,3,4},{2,2}] -> {{1, 2}, {3, 4}}
```

For a non-rectangular brace, `dimensions` / `arrayRank` report only the rectangular prefix common to every child, while `length` always reports the outer element count. A general brace is valid by itself; a Matrix function given such a non-rectangular value emits a Warning and remains unevaluated.

```text
dimensions[{{1,2},{3}}] -> {2}
arrayRank[{{1,2},{3}}] -> 1
length[{{1,2},{3}}] -> 2
at[{{1,2},{3}},0] -> {1, 2}
transpose[{{1,2},{3}}] -> Warning + unevaluated
```

`mget[A,row,col]` is a compatibility alias of `at` and uses the same zero-based indexing.

A brace literal cannot preserve trailing shape information after a leading zero-length dimension. The formatter therefore uses `reshape` only when necessary for round-trip safety.

```text
zeros[0,3]
-> reshape[{}, {0, 3}]

dimensions[zeros[0,3]]
-> {0,3}
```

If evaluation turns Array elements into Arrays, equal child shapes are flattened into the common representation. Mixed scalar/Array leaves or inconsistent child shapes are TypeErrors.

# 19. Aggregate functions

```text
sum
prod
min
max
mean
```

They accept either variadic scalar arguments or a single Array.

```text
sum[1,2,3] -> 6
sum[{1,2,3}] -> 6
prod[] -> 1
sum[] -> 0
mean[1,2,4] -> 7/3
```

`min/max` do not arbitrarily order symbolic values whose ordering cannot be proven.

```text
min[x,3]
-> min[x, 3]
```

Symbolic finite sums of the form `sum[f,{k,a,b}]` are not yet implemented. Explicit finite sequences can be generated with `table` and then aggregated with `sum`.

## 19.1 `range` / `table` / `map`

Use the following functions for exact finite sequence generation and explicit element-wise application.

```text
range[n]
range[a,b]
range[a,b,increment]

table[expr,{i,n}]
table[expr,{i,a,b}]
table[expr,{i,a,b,increment}]

map[f,arrayOrBrace]
```

`range` accepts exact real Integer / Rational bounds and increments. The endpoint is included when it lies in the generated sequence. Floating increments are not silently rounded into a sequence.

```text
range[5] -> {1, 2, 3, 4, 5}
range[0,1,1/3] -> {0, 1/3, 2/3, 1}
range[5,1,-2] -> {5, 3, 1}
```

`table` holds its body and binds only the iterator variable in a local scope for each iteration. An outer definition of the same name is restored after iteration, and nested tables have independent scopes.

```text
table[i^2,{i,5}] -> {1, 4, 9, 16, 25}
table[i/2,{i,0,2,1/2}] -> {0, 1/4, 1/2, 3/4, 1}
```

`map[f,value]` explicitly applies `f[...]` to the **scalar leaves** of an Array or general brace while preserving dense shape or ragged-brace structure. Ordinary calls such as `exp[A]` are deliberately not made element-wise automatically, leaving room for future matrix-function semantics.

```text
map[sin,{0,Pi/2,Pi}] -> {0, 1, 0}
```

---

# 20. Descriptive statistics

Statistical functions generally accept **exact real data** and preserve Rational results when the quantity closes rationally.
Many functions accept either a single rank-1 Array or a scalar argument list.

## 20.1 Order statistics

```text
median
mode
quantile
percentile
iqr
percentrank
```

`quantile` uses Hyndman–Fan Type 7.

```text
median[1,2,3,4] -> 5/2
quantile[1/4,1,2,3,4,5,6,7] -> 5/2
iqr[1,2,3,4] -> 3/2
```

If `mode` has multiple modes, it returns an Array. If every value occurs exactly once, it returns an empty Array.

## 20.2 Variance and standard deviation

```text
var      // population variance
vars     // sample unbiased variance
stddev   // sqrt[var]
stddevs  // sqrt[vars]
```

```text
var[1,2,3] -> 2/3
vars[1,2,3] -> 1
stddev[1,2,3] -> sqrt[6]/3
stddevs[1,2,3] -> 1
```

## 20.3 Other statistics

```text
geomean harmmean rms
mad madR
skew kurtp kurts
cv stderr zscore
trimmean winsor winsorR
cov corr corrspearman
```

Recommended form for `cov/corr/corrspearman`:

```text
cov[{1,2,3},{2,4,6}] -> 4/3
corr[{1,2,3},{2,4,6}] -> 1
```

For compatibility, an even number of scalar arguments may also be split into first-half / second-half datasets.

---

# 21. Assumptions and domains

```text
element[x,Real]
element[x,Integer]
```

The second argument of `simplify/fullSimplify` may be a Predicate, an Array, or an `And`-like condition set.

```text
simplify[sqrt[x^2], element[x,Real]]
-> abs[x]

simplify[abs[x], x >= 0]
-> x
```

Contradictory assumptions produce DomainError.

`element` can prove negative membership as well as positive membership. A non-integral exact Rational is known not to be an Integer; known irrational/transcendental constants are known not to be Rational/Integer; and an algebraic Root whose minimal degree is proven greater than one is likewise non-Rational/non-Integer.

```text
element[1/2,Integer] -> False
element[Pi,Rational] -> False
element[Phi,Rational] -> False
element[root[{-2,0,1},2],Rational] -> False
```

Failure to prove membership is never converted into `False`.

---

# 22. Expression transformation

```text
simplify[expr]
simplify[expr,assumptions]
fullSimplify[expr]
fullSimplify[expr,assumptions]
expand[expr]
factor[expr]
collect[expr,x]
```

Examples:

```text
simplify[sin[x]^2 + cos[x]^2]
-> 1

fullSimplify[x^2 + 2x + 1]
-> (1 + x)^2

expand[(x+1)^3]
-> x^3 + 3 x^2 + 3 x + 1

factor[x^2-1]
-> (x-1)(x+1)
```

`fullSimplify` performs bounded candidate search. A shorter expression is not preferred if obtaining it changes the domain.

```text
fullSimplify[(x^2-1)/(x-1)]
-> original hole preserved

fullSimplify[(x^2-1)/(x-1), x != 1]
-> 1 + x

simplify[1/x-1/x]
-> 1/x-1/x

simplify[1/x-1/x, x != 0]
-> 0

simplify[x^0]
-> x^0

simplify[x^0, x != 0]
-> 1
```

Rules that collapse an expression to a constant while potentially erasing a hole—such as `F-F -> 0`, `0*F -> 0`, `F/F -> 1`, and `F^0 -> 1`—are applied only when `F` is proven defined under the current assumptions. Because mmCal defines `0^0 -> Indeterminate`, symbolic `x^0` is likewise preserved unless `x != 0` is proven.

The same rule applies when a special/combinatorial function degenerates to a constant and a parameter disappears from the value. For example, `hypergeometric2F1[0,2,3,z]`, `polylog[s,0]`, and `binom[x,0]` have constant degenerations, but an undefined disappearing argument is not silently erased.

```text
simplify[hypergeometric2F1[0,2,3,1/x]]
-> hypergeometric2F1[0, 2, 3, 1/x]

simplify[hypergeometric2F1[0,2,3,1/x], x != 0]
-> 1
```

`Power` definedness distinguishes exact Rational exponents: a positive non-integer Rational exponent permits a zero base where the evaluator defines it, while a negative Rational exponent retains a nonzero-base requirement. `zeta[s]` is not treated as globally unknown for definedness; its sole pole is represented by the finite condition `s != 1`.

## 22.1 `series` / `normal` / `toNormal` (v1.5.5 WIP)

```text
series[expr,{x,a,n}]
series[expr,{x,a,n},assumptions]
normal[seriesExpr]
toNormal[expr]
```

`series` constructs a local expansion about `x=a` and retains it as internal `seriesData[...]`. `normal` converts only a top-level `SeriesData` object, discarding the remainder order and returning the retained truncated expression. `toNormal` recursively walks the expression tree and converts supported structured objects nested inside it to ordinary expressions. It handles `SeriesData` inside lists, arrays, and calls, and now also recurses into binding right-hand sides of finite and conditional `SolutionSet` objects while preserving set structure, branch conditions, free variables, multiplicity, and domain metadata. Ordinary expressions and unsupported structures are preserved, so the recursive behavior remains distinct from the compatibility-oriented top-level `normal`. The TPSA kernel composes exact constants, the expansion variable, sums, differences, products, division, integer powers, `exp` / `log` / `sin` / `cos` / `sinh` / `cosh`, `tan/cot/sec/csc`, `tanh/coth/sech/csch`, `expm1/log1p`, `sinc/cosc/tanc`, `sinhc/tanhc/expc`, `log2/log10`, and principal `sqrt` / exact rational powers, supporting Taylor series, Laurent series with a finite principal part, and Puiseux series on an exact rational exponent grid. `log2/log10` lower through the general `log[base,x] = log[x]/log[base]` form for both finite logarithmic Series and the `+Infinity` logarithmic layers. Different Puiseux denominators are re-embedded into an exact LCM grid. Analytic functions compose through coefficient recurrences rather than repeated higher differentiation. A symbolic leading coefficient is inverted only when nonzero status is proved, while principal `log` is expanded only at a proved positive-real or nonreal regular center. Direct trigonometric series honor the current angle mode and explicit `Rad` / `Deg` / `Grad`. At a branch point, non-integer rational powers are restricted to a positive leading coefficient with a simple zero/pole, or to an already branched Puiseux expression. Higher-multiplicity cases such as `sqrt[x^2]` and uncertified negative leading directions remain unevaluated rather than selecting a branch by guesswork.

`Infinity` is also accepted as an expansion center. Here `Infinity` means real `+Infinity`, not a general point on the Riemann sphere; internally the expansion is mapped to `t->0+` with `t=1/x`. Therefore exponent `r` in `seriesData[x,Infinity,...]` denotes `(1/x)^r`, and logarithmic layers denote `log[1/x]^k`. For example, `series[1/(x+1),{x,Infinity,4}]` represents `x^-1-x^-2+x^-3-x^-4+O[x^-5]`. `normal` / `toNormal` emit ordinary powers of `x` rather than leaving intermediate forms such as `(1/x)^(-m)`. `D` incorporates `dt/dx=-t^2`, while `integrate` uses `dx=-t^-2 dt`; integrating an explicit `1/x` term therefore closes into the logarithmic layer as `-log[1/x]`. If the truncation remainder is exactly `O(1/x)`, integration is left unevaluated because the unknown remainder can generate a logarithm that the current `O(t^r)` metadata cannot describe. The initial scope includes rational functions, polynomial growth, `exp[1/x]`, Puiseux powers, and `log[1/x]`. Oscillatory `sin[x]`, essential growth `exp[x]`, and direct `log[x]` require dedicated asymptotic providers and remain unevaluated.

A logarithmic asymptotic provider is also available specifically for real `+Infinity`. When the leading form can be proved to be `A(x)~c(1/x)^r` with positive `c`, `log[A(x)]` is factored as `log[c]+r log[1/x]+log[1+h]` and composed through the existing TPSA/logarithmic layers. This supports `series[log[x],{x,Infinity,n}]`, `log[2x]`, `log[x+1]`, `log[x^2+1]`, `log[sqrt[x]+1]`, powers of `log[x]`, and products such as `log[x]/x^m` in the same `SeriesData` representation. The internal logarithmic basis remains `log[1/x]`, but `normal` / `toNormal` use the known positive-infinity direction to emit ordinary `log[x]`. Transseries requiring negative logarithmic powers such as `1/log[x]`, and negative or complex leading directions that need an additional principal-branch decision, remain unevaluated. Oscillatory `sin[x]` and essential growth `exp[x]` remain outside this provider.

```text
series[log[x+1],{x,Infinity,4}]
-> seriesData[x, Infinity, {0, 1, -1/2, 1/3, -1/4}, 0, 5, 1, {{-1, 0, 0, 0, 0}}]

normal[series[log[x+1],{x,Infinity,3}]]
-> x^(-1)-x^(-2)/2+x^(-3)/3+log[x]

series[1/log[x],{x,Infinity,3}]
-> series[1/log[x], {x, Infinity, 3}]
```

Local special-function Series use a primitive-composition provider. Rather than storing higher-derivative tables per function, mmCal expands the known first-derivative kernel with the existing TPSA machinery, integrates the coefficients, and restores the function value at the center as the constant term. Supported local providers include `erf` / `Si` / `Ei` / `Ci`, `erfc` / `fresnelc` / `fresnels`, `li`, regular-center principal `asin` / `acos` / `atan`, and regular-center `lambertw`. `li` uses `li'(z)=1/log[z]`; inverse trigonometric Series use `asin'(z)=(1-z^2)^(-1/2)`, `acos'(z)=-(1-z^2)^(-1/2)`, and `atan'(z)=1/(1+z^2)` through the same TPSA machinery. `erfc` shares the Gaussian kernel with the complementary sign implied by `erfc[z]=1-erf[z]`; Fresnel C/S pass the DLMF 7.2.7–7.2.8 kernels `cos[Pi z^2/2]` / `sin[Pi z^2/2]` through TPSA. Zero-centered coefficients agree with DLMF 7.6.1, 7.6.4, 7.6.6, and 6.6.5. `SeriesData` also carries coefficient layers for powers of `log(x-center)`, so the logarithmic singularities of `Ei` / `Ci` at zero are now retained exactly using the local DLMF 6.6.1 / 6.6.6 series. The same representation closes over `log[x]`, `x log[x]`, powers of the local logarithm, products with Puiseux factors, and integration of `x^-1 log[x]^k`. For origin-centered composition, if an argument can be certified as `A(t)=c t^r(1+h)` with `0<r<=1` and `c>0`, mmCal factors `log A=log c+r log t+log(1+h)` and composes `log` / `Ei` / `Ci` through Taylor/Puiseux zeros. Cases that can change the principal winding, such as `r>1`, or negative/complex leading coefficients remain unevaluated rather than choosing a branch. Centers on the principal cut described in DLMF 6.2 also remain unevaluated. For principal `li(z)=Ei(Log(z))`, the provider deliberately accepts only proved real centers `x>1` as in DLMF 6.2.8, or provably nonreal centers. Real centers at `0`, `1`, `0<x<1`, and on the negative real axis remain unevaluated rather than inventing a two-sided principal neighborhood. `asin` / `acos` accept real centers only when `-1<a<1` is provable, or centers provably off the real-axis cuts. `atan` accepts real centers, complex centers with provably nonzero real part, and imaginary-axis centers provably between `-I` and `I`; its branch points and the outward principal cuts remain unevaluated. Inverse-trigonometric return values follow the session angle mode, so their Series coefficients include the Radian/Degree/Gradian output scale. By contrast, the defining kernels of `Si` / `Ci` / Fresnel C/S are intrinsically radian-based and do not depend on session angle mode. Lambert W regular-center expansion generates the integer-coefficient derivative polynomials `p_n(W)` from DLMF 4.13.4 and composes regular-center Taylor coefficients through the existing TPSA/Puiseux kernel. At regular centers principal `W_0` avoids the cut `(-Infinity,-1/E]`, while explicit integer branches `k!=0` avoid `(-Infinity,0]`. At the branch point `z=-1/E`, the DLMF 4.13.9_1–4.13.9_2 expansion in `s=sqrt[E z+1]` represents exact-center `z=-1/E` expansions of `W_0` and `W_-1` on a square-root Puiseux grid. The `d_n` recurrence is evaluated as rational even terms and rational-times-`sqrt[2]` odd terms, avoiding general algebraic simplification inside the coefficient-generation loop. The principal sheet is selected only when the leading direction certifies the principal square root; `W_-1` uses the opposite sign of the same local variable. Other branches, directions entering the cut, and higher-multiplicity contacts that would require choosing a value for `sqrt[x^2]` remain unevaluated. `gamma` / `lgamma` use the same local framework. It composes the DLMF 5.7 coefficients for `log Gamma(1+z)` and the half-integer base through the existing TPSA machinery, then transports them to exact integer and half-integer centers using the logarithmic derivative of `Gamma(z+1)=z Gamma(z)`. It does not obtain coefficients by repeated higher derivatives. In particular, `lgamma` builds the local logarithmic increment directly rather than first constructing a Gamma Series and taking its logarithm, avoiding large cancellation expressions at higher orders. `gamma` can compose regular supported centers through complex local directions and Puiseux arguments. Since mmCal's current `lgamma[x]` means real-axis `log[abs[gamma[x]]]`, not complex `LogGamma`, its Series provider accepts only real-coefficient local directions. Gamma poles at non-positive integers, unsupported exact centers such as `1/3`, and complex-direction `lgamma` remain unevaluated. `digamma` / `trigamma` use the same local basis. Digamma coefficients matching DLMF 5.7.4 are obtained by differentiating the same log-Gamma coefficients once; trigamma differentiates them a second time, so no duplicate coefficient table is maintained. Transport from the base centers 1 and 1/2 to exact integer and half-integer centers uses `psi(z+1)=psi(z)+1/z` and `psi1(z+1)=psi1(z)-1/z^2` directly in Series arithmetic. Supported regular centers compose through complex directions and Puiseux arguments, while non-positive integer poles and exact centers without the current coefficient basis remain unevaluated. Origin-centered `polylog[s,z]` is constructed directly from the defining DLMF 25.12.10 power series. When the order `s` does not depend on the expansion variable it may remain symbolic; coefficients `n^(-s)` are preserved as exact expressions while TPSA/Puiseux arguments are composed. The simple zero `Li_s(z)=z+O[z^2]` at the origin is exposed to valuation, so products, quotients, and Laurent inversion compose naturally. Positive-integer orders also support nonzero regular centers. It forms Taylor coefficients directly from `D^n Li_s(z)=z^(-n) sum_k s(n,k) Li_(s-k)(z)` using signed Stirling numbers. Terms reaching nonpositive integer order are reduced exactly with `Li_0(z)=z/(1-z)` and the Eulerian-polynomial identity `Li_{-m}(z)=z A_m(z)/(1-z)^(m+1)`, avoiding both higher-derivative tables and negative-order `polylog` heads in the result. Real centers are accepted only when `a<1` is proved, while provably nonreal centers are also allowed, keeping the principal cut `[1,Infinity)` excluded. For positive-integer order with proved `0<a<1`, positivity of the defining series also supplies nonzero evidence to reciprocal and negative-integer-power inversion. Noninteger orders at nonzero centers, cut centers, and symbolic centers whose regularity cannot be proved remain unevaluated. Known valuations for special functions, inverse trigonometric functions, Lambert W, the supported Gamma/psi family, and polylog are exposed to product, quotient, and Laurent inversion composition.

```text
series[log[x],{x,0,4}]
-> seriesData[x, 0, {0, 0, 0, 0, 0}, 0, 5, 1, {{1, 0, 0, 0, 0}}]

series[Ei[x],{x,0,4}]
-> seriesData[x, 0, {-digamma[1], 1, 1/4, 1/18, 1/96}, 0, 5, 1, {{1, 0, 0, 0, 0}}]

series[Ci[x],{x,0,6}]
-> seriesData[x, 0, {-digamma[1], 0, -1/4, 0, 1/96, 0, -1/4320}, 0, 7, 1, {{1, 0, 0, 0, 0, 0, 0}}]

series[log[2*x],{x,0,4}]
-> seriesData[x, 0, {log[2], 0, 0, 0, 0}, 0, 5, 1, {{1, 0, 0, 0, 0}}]

series[log[sqrt[x]],{x,0,4}]
-> seriesData[x, 0, {0, 0, 0, 0, 0, 0, 0, 0, 0}, 0, 9, 2, {{1/2, 0, 0, 0, 0, 0, 0, 0, 0}}]

series[log[x^2],{x,0,4}]
-> series[log[x^2], {x, 0, 4}]
```

```text
series[lgamma[1+x],{x,0,5}]
-> seriesData[x, 0, {digamma[1], Pi^2/12, -zeta[3]/3, Pi^4/360, -zeta[5]/5}, 1, 6, 1]

series[gamma[1/2+x],{x,0,2}]
-> seriesData[x, 0, {sqrt[Pi], (digamma[1]-2log[2])sqrt[Pi], (Pi^2/2+(digamma[1]-2log[2])^2)sqrt[Pi]/2}, 0, 3, 1]

series[gamma[x],{x,0,3}]
-> series[gamma[x], {x, 0, 3}]

series[lgamma[1+I*x],{x,0,3}]
-> series[lgamma[I x+1], {x, 0, 3}]

series[digamma[1+x],{x,0,5}]
-> seriesData[x, 0, {digamma[1], Pi^2/6, -zeta[3], Pi^4/90, -zeta[5], zeta[6]}, 0, 6, 1]

series[trigamma[1/2+x],{x,0,3}]
-> seriesData[x, 0, {Pi^2/2, -14zeta[3], Pi^4/2, -124zeta[5]}, 0, 4, 1]

series[polylog[2,x],{x,0,6}]
-> seriesData[x, 0, {1, 1/4, 1/9, 1/16, 1/25, 1/36}, 1, 7, 1]

series[polylog[a,sqrt[x]],{x,0,2}]
-> seriesData[x, 0, {1, 2^(-a), 3^(-a), 4^(-a)}, 1, 5, 2]

series[polylog[2,1/2+x],{x,0,3}]
-> seriesData[x, 0, {polylog[2, 1/2], -2log[1/2], 2(1+log[1/2]), 4(-1-2log[1/2])/3}, 0, 4, 1]

series[1/polylog[2,1/2+x],{x,0,2}]
-> seriesData[x, 0, {1/polylog[2, 1/2], 2log[1/2]/polylog[2, 1/2]^2, -(2(1+log[1/2])/polylog[2, 1/2]-4log[1/2]^2/polylog[2, 1/2]^2)/polylog[2, 1/2]}, 0, 3, 1]
```

```text
series[(1+x)^3,{x,0,5}]
-> seriesData[x, 0, {1, 3, 3, 1, 0, 0}, 0, 6, 1]

normal[%]
-> x^3+3x^2+3x+1
```

```text
series[1/(1-x),{x,0,4}]
-> seriesData[x, 0, {1, 1, 1, 1, 1}, 0, 5, 1]

series[1/x,{x,0,3}]
-> seriesData[x, 0, {1, 0, 0, 0, 0}, -1, 4, 1]

normal[%]
-> x^(-1)
```

```text
toNormal[{series[(1+x)^2,{x,0,3}],series[log[x],{x,0,2}]}]
-> {x^2+2x+1, log[x]}

normal[{series[(1+x)^2,{x,0,3}]}]
-> {seriesData[x, 0, {1, 2, 1, 0}, 0, 4, 1]}
```

Inside `SolutionSet`, only bindings are recursively normalized while solution metadata is retained. Low-cost elementary functions are lowered to the existing TPSA algebra rather than acquiring independent coefficient tables.

```text
series[tan[x],{x,0,5}]
-> seriesData[x, 0, {1, 0, 1/3, 0, 2/15}, 1, 6, 1]

series[log2[x],{x,Infinity,3}]
-> seriesData[x, Infinity, {0, 0, 0, 0}, 0, 4, 1, {{-1/log[2], 0, 0, 0}}]
```

`toNormal` is recursive and idempotent for converted objects. For `SolutionSet`, only binding right-hand sides are normalized; conditions, free variables, multiplicity, and domains are preserved. Additional structured representations can be added to the same frontend in the future.


```text
series[exp[x],{x,0,5}]
-> seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120}, 0, 6, 1]

series[1/sin[x],{x,0,5}]
-> seriesData[x, 0, {1, 0, 1/6, 0, 7/360, 0, 31/15120}, -1, 6, 1]

series[log[x],{x,I,3}]
-> seriesData[x, I, {I Pi/2, -I, 1/2, I/3}, 0, 4, 1]

series[sqrt[1+x],{x,0,6}]
-> seriesData[x, 0, {1, 1/2, -1/8, 1/16, -5/128, 7/256, -21/1024}, 0, 7, 1]

series[(1+x)^(3/2),{x,0,6}]
-> seriesData[x, 0, {1, 3/2, 3/8, -1/16, 3/128, -3/256, 7/1024}, 0, 7, 1]
```

```text
series[sqrt[x],{x,0,5}]
-> seriesData[x, 0, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 1, 11, 2]

series[sqrt[x]*(1+x),{x,0,4}]
-> seriesData[x, 0, {1, 0, 1, 0, 0, 0, 0, 0}, 1, 9, 2]

series[exp[sqrt[x]],{x,0,3}]
-> seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120, 1/720}, 0, 7, 2]
```

```text
series[erf[x],{x,0,7}]
-> seriesData[x, 0, {2/sqrt[Pi], 0, -2/sqrt[Pi]/3, 0, 1/(5sqrt[Pi]), 0, -1/(3sqrt[Pi])/7}, 1, 8, 1]

series[Si[x],{x,0,7}]
-> seriesData[x, 0, {1, 0, -1/18, 0, 1/600, 0, -1/35280}, 1, 8, 1]

series[Ei[1+x],{x,0,4}]
-> seriesData[x, 0, {Ei[1], E, 0, E/6, -E/12}, 0, 5, 1]

series[Ci[1+x],{x,0,4}]
-> seriesData[x, 0, {Ci[1], cos[1 Rad], (-cos[1 Rad]-sin[1 Rad])/2, (cos[1 Rad]/2+sin[1 Rad])/3, (-cos[1 Rad]/2-5sin[1 Rad]/6)/4}, 0, 5, 1]

series[erfc[x],{x,0,5}]
-> seriesData[x, 0, {1, -2/sqrt[Pi], 0, 2/(3sqrt[Pi]), 0, -1/sqrt[Pi]/5}, 0, 6, 1]

series[fresnelc[x],{x,0,5}]
-> seriesData[x, 0, {1, 0, 0, 0, -Pi^2/40}, 1, 6, 1]

series[fresnels[x],{x,0,7}]
-> seriesData[x, 0, {Pi/6, 0, 0, 0, -Pi Pi^2/336}, 3, 8, 1]

series[li[2+x],{x,0,2}]
-> seriesData[x, 0, {li[2], 1/log[2], -1/(2log[2]^2)/2}, 0, 3, 1]

series[li[a+x],{x,0,2},a>1]
-> seriesData[x, 0, {li[a], 1/log[a], -1/(a log[a]^2)/2}, 0, 3, 1]

series[asin[x],{x,0,7}]
-> seriesData[x, 0, {1, 0, 1/6, 0, 3/40, 0, 5/112}, 1, 8, 1]

series[atan[1+I+x],{x,0,3}]
-> seriesData[x, 0, {atan[1+I], 1/5-2I/5, -1/25+7I/25, -1/375-68I/375}, 0, 4, 1]

series[asin[a+x],{x,0,2},-1<a<1]
-> seriesData[x, 0, {asin[a], (1-a^2)^(-1/2), a*(1-a^2)^(-1/2)/(2(1-a^2))}, 0, 3, 1]

series[lambertw[x],{x,0,7}]
-> seriesData[x, 0, {1, -1, 3/2, -8/3, 125/24, -54/5, 16807/720}, 1, 8, 1]

series[lambertw[E+x],{x,0,4}]
-> seriesData[x, 0, {1, exp[-1]/2, -3*exp[-2]/16, 19*exp[-3]/192, -185*exp[-4]/3072}, 0, 5, 1]

series[1/lambertw[x],{x,0,5}]
-> seriesData[x, 0, {1, 1, -1/2, 2/3, -9/8, 32/15, -625/144}, -1, 6, 1]

series[lambertw[-1/E+x],{x,0,3}]
-> seriesData[x, 0, {-1, sqrt[2]sqrt[E], -2*E/3, 11*E sqrt[2]sqrt[E]/36, -43*exp[2]/135, 769*exp[2]sqrt[2]sqrt[E]/4320, -1768*E exp[2]/8505}, 0, 7, 2]

series[lambertw[-1,-1/E+x],{x,0,3}]
-> seriesData[x, 0, {-1, -sqrt[2]sqrt[E], -2*E/3, -11*E sqrt[2]sqrt[E]/36, -43*exp[2]/135, -769*exp[2]sqrt[2]sqrt[E]/4320, -1768*E exp[2]/8505}, 0, 7, 2]
```

In interactive auto/multi layout, `series[sqrt[x]*(1+x),{x,0,4}]` is displayed as `x^(1/2) + x^(3/2) + O[x^(9/2)]`. `single` and machine-facing output retain canonical `seriesData[...]`.

Requests centered on the principal-log branch cut, such as `series[log[x],{x,-1,n}]`, and requests such as `series[sqrt[x^2],{x,0,n}]` that cannot be certified as one principal local branch remain unevaluated.

`SeriesData` supports direct coefficient arithmetic under `D` / `integrate` with respect to its expansion variable. Ordinary coefficients and every `log(x-center)^k` layer are transformed by the product rule and the exact integration recurrence, so Taylor/Laurent/Puiseux exponent grids remain intact while expressions such as `D[log[x]^k]` and `integrate[x^-1 log[x]^k,x]` stay inside SeriesData.

```text
D[series[exp[x],{x,0,5}],x]
-> seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24}, 0, 5, 1]

integrate[series[sqrt[x],{x,0,4}],x]
-> seriesData[x, 0, {2/3, 0, 0, 0, 0, 0, 0, 0}, 3, 11, 2]
```

Integrating an `x^(-1)` term moves it into a logarithmic layer; for example, `integrate[series[1/x,{x,0,4}],x]` returns SeriesData representing `log[x]`. Likewise `x^-1 log[x]^k` maps exactly to `log[x]^(k+1)/(k+1)`.

`seriesData[variable,center,coefficients,minExponent,orderNumerator,exponentDenominator]` remains the ordinary internal representation. When logarithmic terms are present, a seventh argument `logarithmicCoefficientLayers` is appended; layer `k` is the coefficient sequence multiplying `log(x-center)^(k+1)` on the same exponent grid. If all log layers are zero, the six-argument canonical form is retained. Normal use should go through `series` / `normal`. Coefficient `i` has exponent `(minExponent+i)/exponentDenominator`, and the remainder order is `orderNumerator/exponentDenominator`.

---

# 23. Symbolic and numerical differentiation

## 23.1 `D`

```text
D[expr,x]
D[expr,{x,n}]
D[expr,x,y,...]
```

`{x,n}` means the non-negative integer `n`-th derivative. Multiple specifications are applied from left to right.

```text
D[sin[x],{x,4}]
-> sin[x]

D[x^2 y^3,x,y]
-> 6 x y^2
```

`D` has HoldAll semantics, so existing variable values do not replace the symbols being differentiated.

```text
D[x^3 + 2x,x]
-> 2 + 3 x^2

D[exp[x^2],x]
-> 2 x exp[x^2]
```

With radians as the default:

```text
D[sin[x],x]
-> cos[x]
```

Functions such as `abs/sign/re/im/conj/arg`, which do not in general have an ordinary holomorphic derivative with respect to a complex variable, are not given fabricated derivatives; unevaluated `D[...]` is retained instead.

Stable functions whose removable singularity is filled at zero retain that point after differentiation. For example:

```text
D[sinc[x],x]
-> cases[(x cos[x]-sin[x])/x^2 if x != 0; 0 if x == 0]

D[cosc[x],x]
-> cases[(-1+x sin[x]+cos[x])/x^2 if x != 0; 1/2 if x == 0]
```

These are first-class scalar `cases[...]` values, not evaluation-control `if[...]` and not `SolutionSet` objects. `if[...]` selects one branch during evaluation; `cases[...]` preserves mathematical branch conditions symbolically; `SolutionSet` remains reserved for equation solution sets.

Integration reuses the same derivative knowledge.

```text
D[integrate[f[x],x],x]
-> f[x]

D[integrate[t^2,{t,0,x}],x]
-> x^2
```

For general forms with variable endpoints or where the differentiation variable occurs inside the integral, the Leibniz rule is constructed formally. A quotient whose denominator is independent of the differentiation variable uses `f'/c` directly rather than expanding into the general quotient rule.

## 23.2 `diff`

```text
diff[expr,x,at]
diff[expr,x,at,digits]
```

This is not a separate finite-difference formula. mmCal first constructs an exact derivative Expr using `D`, then passes that expression to the `CertifiedEvaluator` at the specified point. For finite-precision inputs, `CertifiedEnclosure` and `InformationEnclosure` are evaluated separately, so the result cannot claim more precision than the input carries. Held `N[...]` subexpressions obey the same rule. If the available input information cannot determine a singularity or branch-cut side, `diff` returns an EvaluationError rather than guessing a value.
The fourth argument `digits` counts fractional decimal places; it is not the significant-digit convention used by `N[...,p]`. If a finite-precision input carries less information, the `InformationEnclosure` remains the upper bound on the precision that can be claimed.

```text
diff[x^2,x,3]
-> 6.0
```

---

# 24. Symbolic integration and certified numerical integration

## 24.1 `integrate` — exact / symbolic integration

```text
integrate[expr,x]
integrate[expr,{x,a,b}]
integrate[expr,x,assumptions]
integrate[expr,{x,a,b},assumptions]
```

The first form is an indefinite integral, and the second is an exact/symbolic definite integral. The integration variable is held as a binder and is not replaced by a global definition with the same name.

An indefinite integral returns one representative of the antiderivative modulo an additive constant. Therefore `+ C` is not displayed.

```text
integrate[x^2,x]
-> x^3 / 3

integrate[(2x+3)^5,x]
-> (3 + 2 x)^6 / 12

integrate[1/(2x+3),x]
-> log[3 + 2 x] / 2
```

Major exact rules currently implemented:

- Constants, `x`, and arbitrary finite polynomials
- Rational powers of affine bases; exponent `-1` is mapped to Log
- Rational functions with Rational coefficients. Linear/quadratic factors use exact partial fractions, including repeated irreducible quadratics `(a x^2+b x+c)^k` through an exact completing-the-square recurrence. Denominators containing factors of degree three or higher use a Q[x] square-free decomposition plus Hermite reduction; the remaining square-free part is represented as a finite algebraic-log sum over certified complex `root[...,k,Complex]` values with exact residues `P(r)/Q'(r)`. Elementary reverse-chain rules run first so compact `atan/asin` forms are preferred when available. The specialized rational path currently uses a degree-12 work budget. Forms such as `x^m/(1+x^n)` may still fall through to a `hypergeometric2F1` primitive when appropriate
- Quadratic inverse-square-root forms with provably positive Rational scale, plus `sqrt[q(x)]` primitives for exact Rational quadratics `q(x)`
- Finite Fourier reduction for `sin^m/cos^n`; the integrator may explicitly expand positive integer total degree up to 256
- `sin[u]^(-n)` / `cos[u]^(-n)` (`1<=n<=256`) through the standard `csc/sec` reduction recurrences
- Positive integer powers (`2<=n<=256`) of `tan/cot/sec/csc` through the standard reduction formulas
- Linearity over sums, differences, negation, and factors independent of the integration variable
- Safe standard primitives for `exp/sin/cos/tan/cot/sec/csc`
- Safe standard primitives for `sinh/cosh/tanh/coth/sech/csch`
- `log/log1p/expm1/sqrt/cbrt`; logarithmic-derivative knowledge recognizes forms such as `log[x]/x` and `1/(x log[x])`
- `asin/acos/atan/asinh/acosh/atanh`
- `erf/erfc`
- `fresnelc/fresnels`; exact Rational-coefficient `sin/cos[a x^2+b x+c]` phases are completed to a square and reduced to standard Fresnel integrals; the defining `Pi*x^2/2` kernels are recognized directly
- Gaussian family: for exact positive Rational `a`, `exp[-a x^2]` prefers the canonical `erf` primitive over generic 1F1, and improper endpoints close through the the exact `erf` endpoint limits at ±Infinity
- `hypergeometric1F1`; outside the preferred Gaussian family, `exp[c x^n]` with positive integer `n>=2` reduces to an entire 1F1 primitive at the origin
- Exact inverse chain rule; `f'(x) f(x)^p` is also recognized structurally instead of depending on the accidental post-`D` expression shape
- Finite integration by parts for polynomial × `exp/sin/cos/sinh/cosh`
- Exact integration of `exp[a x+b] sin/cos[c x+d]` forms by solving a linear system
- Local principal-branch `t=sqrt[x]` substitution for rational forms involving `sqrt[x]`, including the quadratic nested-radical class
- Bounded Weierstrass substitution `t=tan(theta/2)` for rational expressions in a common `sin(theta),cos(theta)` argument; the transformed expression is handed to the exact rational integrator
- Dilogarithm knowledge reducing `log[1+beta*x^n]/x` to `polylog[2,-beta*x^n]`
- Bounded distribution of products over short sums, while exact whole-expression chain rules are tried first to avoid capability regressions
- `x^n log[x]` for non-negative integer `n`

Examples:

```text
integrate[exp[2x+1],x]
-> exp[1 + 2 x] / 2

integrate[sin[2x],x]
-> -cos[2 x] / 2

integrate[tan[x],x]
-> -log[cos[x]]

integrate[sech[x],x]
-> atan[sinh[x]]

integrate[1/(1-x^2),x]
-> atanh[x]

integrate[1/(x^2-1),x]
-> -atanh[x]

integrate[(x+1)/(x+2),x]
-> x - log[2 + x]

integrate[(x+1)/(x^2+4),x]
-> log[4 + x^2] / 2 + atan[x/2] / 2

integrate[1/sqrt[4-x^2],x]
-> asin[x/2]

integrate[1/sqrt[x^2+4],x]
-> asinh[x/2]

integrate[sin[x]^2,x]
-> (x - sin[2 x] / 2) / 2

integrate[log[x]/x,x]
-> log[x]^2/2

integrate[sec[x]^3,x]
-> sec[x]tan[x]/2+log[sec[x]+tan[x]]/2

integrate[exp[x^6],x]
-> x hypergeometric1F1[1/6, 7/6, x^6]

integrate[sin[2x]^(-2),x]
-> -cot[2x]/2

integrate[cos[4x^2],x]
-> fresnelc[x sqrt[8/Pi]]/sqrt[8/Pi]

integrate[sin[2x^2]^4,x]
-> 3x/8+fresnelc[x sqrt[16/Pi]]/(8sqrt[16/Pi])-fresnelc[x sqrt[8/Pi]]/sqrt[8/Pi]/2

integrate[asin[x],x]
-> x asin[x] + sqrt[1 - x^2]

integrate[erf[x],x]
-> x erf[x] + exp[-x^2] / sqrt[Pi]

integrate[x exp[x],x]
-> x exp[x] - exp[x]

integrate[x^2 log[x],x]
-> x^3 log[x] / 3 - x^3 / 9

integrate[x^2+sin[x],x]
-> x^3 / 3 - cos[x]

integrate[E^x cos[x],x]
-> (cos[x] + sin[x]) / 2 * exp[x]

integrate[1/(x^3+1),{x,0,1}]
-> log[2] / 3 + Pi sqrt[3] / 9
```

Representative substitution and connection cases added in this batch:

```text
integrate[sqrt[4-x^2],x]
-> x sqrt[4-x^2]/2+2asin[x/2]

integrate[2*x*(1+x^2)^5,x]
-> (1+x^2)^6/6

integrate[x/(1+x^4),x]
-> atan[x^2]/2

integrate[cos[Pi*x^2/2],x]
-> fresnelc[x]

integrate[log[1+x^2]/x,x]
-> -polylog[2,-x^2]/2

integrate[1/(1+sin[x]),x]
-> -2/(1+tan[x/2])
```

Angle semantics are shared with the existing `AngleSemantics` / `D` infrastructure.

```text
integrate[sin[x Deg],x]
-> -(180 / Pi cos[x Deg])
```

When a candidate antiderivative is discovered structurally, for example through an inverse-chain rule, the existing `D` implementation is used as a proof engine. The candidate is accepted only after proving an exact proportional relationship between its derivative and the original integrand. Agreement at a finite set of numerical sample points is not sufficient.

**Global simplifications** that would violate branch or definedness semantics are not performed. An antiderivative, however, need not obey the same standard as a global algebraic identity. A locally valid primitive on a common analytic region may be accepted as an integration-specific rule.

For example, the principal-square-root integral

```text
integrate[1/sqrt[x^2-1],x]
-> log[sqrt[x^2-1] + x]
```

may be returned, but the Simplifier is not given the globally unsafe rule

```text
sqrt[x^2-1] == sqrt[x-1] sqrt[x+1]
```

Integration formulas and global algebraic identities are therefore treated under separate validity criteria.

Linearity over sums also supports partial evaluation.

```text
integrate[x^2 + gamma[x],x]
-> x^3 / 3 + integrate[gamma[x],x]
WARN: integrate partially evaluated the expression; remaining subintegral(s) are outside the current symbolic rule set
```

Only the unresolved part is retained; successfully integrated terms are not rolled back.

### Unevaluated-integration diagnostics

Unevaluated results are no longer collapsed into one warning. The current diagnostic codes are:

- `integrate::unsupported` — no current symbolic rule matched. **This does not mean no closed form exists.**
- `integrate::partial` — part of the expression was integrated, but remaining subintegrals are outside the current rules.
- `integrate::conditionsRequired` — additional domain or branch assumptions are required before choosing a safe primitive.
- `integrate::noKnownClosedForm` — an explicitly recognized family has no known finite closed form in mmCal's currently supported standard-function vocabulary.

The last category is deliberately narrower than a claim of absolute mathematical impossibility. A series, a newly introduced special function, or a broader function class may still represent the antiderivative.

### Nested square-root substitution

The following class is currently supported:

```text
integrate[sqrt[x + sqrt[x]],x]
-> 2 (x + sqrt[x]) sqrt[x + sqrt[x]] / 3
   - ((1 + 2 sqrt[x]) sqrt[x + sqrt[x]] / 4
      - log[1 + 2 sqrt[x + sqrt[x]] + 2 sqrt[x]] / 8)
```

This is not a one-off formula. It handles the class that becomes `2 t sqrt[q(t)]` under `t=sqrt[x]`, where `q` is an exact Rational-coefficient quadratic with a positive leading coefficient. It is a local substitution rule on the principal branch, not a general radical-substitution search.

### Exact / symbolic definite integrals

When a safe antiderivative is available, exact substitution is performed at the endpoints and the difference is taken. For functions with domain constraints, the CertifiedEvaluator is used where possible to preflight the **entire interval**, so intermediate poles or branch violations are not missed.

```text
integrate[x^2,{x,0,1}]
-> 1/3

integrate[sin[x],{x,0,Pi}]
-> 2

integrate[1/x,{x,1,2}]
-> log[2]

integrate[log[x],{x,1,Pi}]
-> 1 + Pi log[Pi] - Pi
```

If the interval crosses a singularity, endpoint substitution alone is not used to manufacture a result.

```text
integrate[1/x,{x,-1,1}]
-> WARN + unevaluated integrate[...]

integrate[tan[x],{x,0,2}]
-> WARN + unevaluated integrate[...]
```

Cauchy principal value is not assumed automatically.

### Assumptions and improper integrals

An assumption may be supplied as the third argument. It is integrated into the existing `KnowledgeContext` and can be used for branch-sensitive simplification such as `abs` and `sqrt[x^2]`.

```text
integrate[abs[x],x,x>=0]
-> x ^ 2 / 2

integrate[abs[x],x,x<=0]
-> -x ^ 2 / 2

integrate[sqrt[x^2],x,x>=0]
-> x ^ 2 / 2
```

If an endpoint is `Infinity` or ordinary substitution is undefined there, the corresponding one-sided or infinite limit of the antiderivative is used to evaluate the improper integral. Only classes for which the absence of internal singularities can be proven are accepted; Cauchy principal values are not guessed.

```text
integrate[exp[-x],{x,0,Infinity}]
-> 1

integrate[1/x^2,{x,1,Infinity}]
-> 1

integrate[1/(1+x^2),{x,-Infinity,Infinity}]
-> Pi

integrate[1/sqrt[x],{x,0,1}]
-> 2

integrate[log[x],{x,0,1}]
-> -1

integrate[exp[-a*x],{x,0,Infinity},a>0]
-> 1/a

integrate[x^(s-1)*exp[-x],{x,0,Infinity},s>0]
-> gamma[s]

integrate[x^(a-1)*(1-x)^(b-1),{x,0,1},{a>0,b>0}]
-> beta[a,b]

integrate[1/(1+x^4),{x,0,Infinity}]
-> Pi/(2sqrt[2])

integrate[log[x]^2,{x,0,1}]
-> 2

integrate[sin[x]/x,{x,0,Infinity}]
-> Pi/2

integrate[cos[Pi*x^2/2],{x,0,Infinity}]
-> 1/2
```

These parameterized families reduce to Gamma/Beta/Mellin forms only when the assumptions prove the required convergence and real-branch conditions. A symbolic finite-interval pole is likewise accepted only when assumptions prove it lies outside the interval.

If an internal pole cannot be excluded, the integral remains unevaluated.

```text
integrate[1/(x-2),{x,1,Infinity}]
-> WARN + unevaluated integrate[...]
```

## 24.2 `limit` — exact / symbolic limits

```text
limit[expr,x,a]
limit[expr,x,a,-1]   // left-hand limit
limit[expr,x,a,1]    // right-hand limit
limit[expr,{x,a,direction}]
```

The fourth argument, or the third element of the brace form, specifies direction: `-1` means left and `1` means right. If omitted, a two-sided limit is requested. `limit[expr,{x,a,direction}]` is equivalent to the four-argument form. Direction is not merely a display option; it is supplied to `KnowledgeContext` as a temporary assumption `x<a` / `x>a`.

```text
limit[sin[x]/x,x,0]
-> 1

limit[(1-cos[x])/x^2,x,0]
-> 1/2

limit[1/x,x,0,1]
-> Infinity

limit[1/x,x,0,-1]
-> -Infinity

limit[abs[x]/x,x,0,1]
-> 1

limit[abs[x]/x,x,0,-1]
-> -1

limit[atan[x],x,Infinity]
-> Pi / 2

limit[exp[-x],x,Infinity]
-> 0

limit[sin[1/x],x,0]
-> Indeterminate

limit[sin[1/x],x,0,1]
-> Indeterminate

limit[x sin[1/x],x,0]
-> 0

limit[x*Ei[x],x,0,1]
-> 0

limit[sqrt[x]*log[x],x,0,1]
-> 0

limit[Ei[x]-log[x],x,0,1]
-> -digamma[1]

limit[Ci[x]-log[x],x,0,1]
-> -digamma[1]

limit[log[2*x]-log[x],x,0,1]
-> log[2]

limit[Ei[x],x,0]
-> -Infinity

limit[Ei[x],x,Infinity]
-> Infinity

limit[Ei[x],x,-Infinity]
-> 0

limit[Ci[x],x,0]
-> -Infinity

limit[Ci[x],x,Infinity]
-> 0

limit[Ci[x],x,-Infinity]
-> I Pi

limit[li[x],x,0]
-> 0

limit[li[x],x,1]
-> -Infinity

limit[li[x],{x,1,1}]
-> -Infinity

limit[li[x],x,Infinity]
-> Infinity

limit[li[x],x,-Infinity]
-> ComplexInfinity
```

`Ei`, `Ci`, and `li` are interpreted on their principal branches. Along the real axis, `Ei[x]` tends to `-Infinity` at zero and to 0 as `x -> -Infinity`. `Ci[x]` tends to `-Infinity` at zero and to 0 at positive infinity, while the negative real axis lies on its branch cut and `Ci[x] -> I Pi` as `x -> -Infinity`. `li[x]` tends to 0 at the origin. For `li[x]` along the negative real axis as `x -> -Infinity`, the current direction representation does not collapse the unbounded complex value to real `-Infinity`; mmCal returns the direction-unspecified complex infinity `ComplexInfinity`.

For indeterminate `0/0` forms, repeated l'Hôpital evaluation using the existing `D` implementation is available with safety limits. Local zero/pole orders of Rational functions and degree comparisons at infinity are handled exactly. For finite points that remain unresolved by these specialized rules, a supported local `SeriesData` may be used as a supplemental backend to prove a finite constant or zero from the leading nonzero term. This also resolves cancellations such as `Ei[x]-log[x]` that look like `Infinity-Infinity` termwise. Negative powers or constant-order `log^k` terms are kept unevaluated rather than guessing a divergence direction, so the Series path does not replace the existing limit kernel. For `sin` / `cos` / `tan`, when the real argument is proved to run to `+/-Infinity` on a one-sided approach or at infinity, periodic oscillation proves that no single limiting value exists and the result is `Indeterminate`; for a finite two-sided limit, oscillation on either side is already sufficient. The engine also uses the real bound `abs(sin[u]), abs(cos[u]) <= 1` for exact Rational-function arguments to close two-factor squeeze cases such as `x sin[1/x] -> 0`. Unresolved two-sided limits are not collapsed into principal values or other guessed results.

```text
limit[1/x,x,0]
-> WARN + limit[1 / x, x, 0]
```

## 24.3 `nintegrate` — certified numerical integration

```text
nintegrate[expr,{x,a,b}]
nintegrate[expr,{x,a,b},digits]
```

```text
nintegrate[x^2,{x,0,1},12]
-> 0.333333333333
```

Internally, the interval is normalized with `x=a+(b-a)t`. Newton–Cotes quadrature on exact Rational grids is combined with certified bounds on higher derivatives to enclose the quadrature error.

This is not a `double`-based Simpson calculation that returns a value because it merely appears close. A result is returned only when the final rounding is uniquely determined. For finite-precision bounds or integrands, `CertifiedEnclosure` and `InformationEnclosure` are integrated in parallel, preventing the result from claiming more information than its inputs. Held `N[...]` subexpressions retain the same information limit.
The third argument `digits` counts fractional decimal places. Thus, when the input information is sufficient, a large integer part does not reduce the requested number of decimal places. If finite-precision input information is the tighter limit, the output is capped by that information instead.

Before constructing higher derivatives, the original integrand is preflighted over the full interval so obvious singularities can be rejected early. No accidental cancellation across a singularity is accepted.

---

# 25. Solver

```text
solve[equation,x]
solve[equation,x,domainOrConstraint]
solve[equation,domain]
solve[{equations...},{variables...}]
```

The default ambient domain for equation systems is Complex.
Ordered inequalities are handled over Real or a real subdomain.

`solve[equation,domain]` is a shorthand for `Integer` / `Rational` / `Real` / `Complex` domains. It infers the solve variable only when the equation contains **exactly one** unknown user symbol. Zero or multiple candidates produce TypeError rather than a guess. Explicit-variable forms likewise reject protected constants and builtin/domain symbols such as `Pi` or `Real` as solve variables.

```text
solve[x^2 == 1,x]
-> {x == 1, x == -1}

solve[x^2 + 1 == 0,x,Real]
-> {}

solve[x^2 + 1 == 0,x,Complex]
-> {x == I, x == -I}

solve[x^2 < 4,x]
-> {x in Real if x > -2 && x < 2}

solve[{2x+3y==5,x-2y==9},{x,y}]
-> {{x==37/7, y==-13/7}}
```

`SolutionSet` distinguishes Empty / Finite / Universal / Conditional / Unresolved.
Unsupported expressions are not misreported as having no solutions.

### Multivariate polynomials and Gröbner bases

Exact rational multivariate polynomials use a dedicated `Q[x1,...,xn]` ring representation.

```text
groebnerBasis[polys,{x,y,...}]
groebnerBasis[polys,{x,y,...},Lex]
groebnerBasis[polys,{x,y,...},GrLex]
groebnerBasis[polys,{x,y,...},GrevLex]
polynomialReduce[f,G,{x,y,...}]
polynomialReduce[f,G,{x,y,...},order]
```

The default term order is `GrevLex`. The exact backend implements multivariate division, normal forms, S-polynomials, Buchberger product/chain criteria, sugar pair selection, interreduction, and reduced Gröbner bases. Variable and degree counts are not hard-coded, but pair, reduction, basis-size, and term-count work is bounded.

```text
groebnerBasis[{x y-1,y^2-x},{x,y},Lex]
-> {x-y^2, y^3-1}

polynomialReduce[x^2+y^2,{x-y,y^2-1},{x,y},Lex]
-> {{x+y, 2}, 2}

polynomialReduce[x y-1,groebnerBasis[{x y-1,y^2-x},{x,y},Lex],{x,y},Lex]
-> {{y, 1}, 0}
```

A direct `groebnerBasis[...]` result may be composed inside `polynomialReduce[...]` (or another `groebnerBasis[...]`) without evaluating the held polynomial variables against current session bindings.

For nonlinear polynomial equality systems, `solve[{...},{...}]` attempts Lex Gröbner elimination. Inconsistent ideals return `{}`; zero-dimensional shape-position bases are enumerated through exact univariate roots plus back substitution and exact remainder verification. If a basis is not in shape position but an exact univariate eliminant exists, mmCal may instead specialize the original system at each exact eliminant root and recursively solve every lower-dimensional branch; the result is returned only when all branches close completely within the bounded-work limits.

Positive-dimensional systems are not assigned guessed algebraic-variety parameterizations. Exact `f*g==0` structure may be split into the complete union `f==0` / `g==0`, yielding free-variable branches such as `solve[{x*y==0},{x,y}] -> {x == 0 where y in Complex, y == 0 where x in Complex}`. Multi-equation systems may additionally eliminate a solver variable when it occurs with a provably nonzero constant linear coefficient; the exact binding is substituted into the remaining system, which is then solved recursively. For example, `solve[{z-x-y==0,x^2+y^2==1},{x,y,z},Real]` reduces to the circle projection while preserving `x in [-1,1]`. When the reduced problem is a single polynomial equation, it may be projected through one degree-1/2 symbolic-coefficient variable and lifted with the remaining solver variables free. In the default Complex domain, for example, `solve[{x*y==1},{x,y}] -> {y == 1/x where x in Complex if x != 0}`. In the Real domain, free parameters are kept Real; linear projections retain denominator nonzero conditions, while quadratic projections retain the nonnegative discriminant/radicand required by the principal `sqrt`. With one free parameter, that condition is normalized through the exact univariate polynomial-inequality solver, so `solve[{x^2+y^2==1},{x,y},Real]` returns the complete branches `y==±sqrt[1-x^2]` with `-1<=x<=1`. With several free parameters, the semialgebraic predicate is preserved explicitly. Nonconstant linear coefficients, higher-degree function-field algebraic equations, and general rational parameterizations remain unresolved.

Because `solve` is `HoldAll`, its input is not sent through the general Evaluator before classification. A dedicated **solve-safe normalization** layer canonicalizes builtin aliases and applies proof-safe Simplifier rewrites only. This avoids evaluation side effects while making equivalent spellings such as `E^x` / `exp[x]`, `ln` / `log`, and `log2` / `log10` share solver capabilities.

Denominator zeros, Log definedness, rational-function holes/poles, and related constraints are retained as global conditions where possible.

When current Predicate representation cannot completely express a condition, such as the pole set of Gamma, the solver leaves the result unresolved rather than fabricating an incomplete condition.

### Global inverse solving on the real axis

`MathRegistry` stores not only principal inverses but also metadata about global injectivity, monotonicity, and real-valued ranges. On **Real or a real subdomain**, `solve` safely applies inverses only to functions that can be proven globally one-to-one.

```text
solve[exp[x]==2,x,Real]
-> {x == log[2]}

solve[log[x]==2,x,Real]
-> {x == exp[2]}

solve[sinh[3x]==2,x,Real]
-> {x == asinh[2] / 3}

solve[tanh[x]==1/2,x,Real]
-> {x == atanh[1/2]}

solve[tanh[x]==2,x,Real]
-> {}

solve[exp[x]==a,x,Real]
-> {x == log[a] if a in Real && a > 0}

solve[E^x==8,x,Real]
-> {x == log[8]}

solve[ln[x]==2,x,Real]
-> {x == exp[2]}

solve[log2[x]==3,x,Real]
-> {x == 8}
```

### Real-domain nonexistence / uniqueness proof

Real equalities pass through an exact proof layer after the more specific inverse-function solvers. This layer never treats failure to find a numerical root as evidence of nonexistence. It currently attempts, in increasing cost order:

- global strict sign/nonzero facts for the residual `f(x)=lhs-rhs` from `ValueFacts`;
- exclusion by a registered real function range, or exact-anchor inversion for a globally injective real function even when no inverse builtin exists;
- strict monotonicity on a connected real-domain piece where `f` and `f'` are real and defined, with an exactly verified root anchor; a non-strict derivative sign `f'>=0` / `f'<=0` is also promoted to strict monotonicity when the complete zero set of `f'` is proved finite or at most countable through Integer-parameter families, hence contains no real interval;
- strict convexity/concavity when `f''` has a strict global sign and the unique critical point can be constructed exactly, allowing the global extremum value to prove either nonexistence or a tangent unique root.

The univariate real analyzer decomposes algebraically provable domain conditions into connected intervals and splits those intervals again at exact critical points. On each piece it can certify strict monotonicity, one-sided endpoint limits, and the resulting range. For example, `log[x]` is increasing on `(0,Infinity)` with range `(-Infinity,Infinity)`, `sqrt[x]` is increasing on `[0,Infinity)` with the same nonnegative range, and `atanh[x]` is increasing on `(-1,1)` with range `(-Infinity,Infinity)`. A disconnected domain such as `1/(x^2-1)` is never joined across its poles.

This decomposition is currently bounded to finitely many algebraic boundaries and exactly solvable critical points. It returns Unknown rather than inventing completeness for domains with infinite periodic pole sets or expressions where multiple complex-valued subexpressions could cancel back to real values. Conversely, proving uniqueness is not enough to manufacture a symbolic root: `erf[x]==1/2` remains `UnresolvedSolutionSet[x]` because mmCal currently has neither inverse-erf nor a general transcendental Root representation.

The real monotonicity metadata for `erf` uses DLMF 7.10.1, `erf'(x)=2 exp(-x^2)/sqrt(Pi)>0`. The layer is conceptually informed by the exact set/range/optimization roles of Wolfram Language `Reduce`, `FunctionRange`, and global optimization, but mmCal accepts only certificates supported by its own exact knowledge and symbolic differentiation infrastructure.

```text
solve[exp[x]==x,x,Real] -> {}
solve[exp[x]+x^2+1==0,x,Real] -> {}
solve[x+exp[x]-1==0,x,Real] -> {x == 0}
solve[exp[x]==x+1,x,Real] -> {x == 0}
solve[exp[x]-x+5==0,x,Real] -> {}
solve[log[x]-x-1==0,x,Real] -> {}
solve[log[x]-x+1==0,x,Real] -> {x == 1}
solve[1/x+x==0,x,Real] -> {}
solve[1/x-x==0,x,Real] -> {x == -1, x == 1}
solve[sin[x]==x,x,Real] -> {x == 0}
solve[erf[x]==0,x,Real] -> {x == 0}
solve[erf[x^2-1]==0,x,Real] -> {x == 1, x == -1}
solve[erf[x]==1,x,Real] -> {}
solve[erfc[x]==1,x,Real] -> {x == 0}
solve[erfc[x]==0,x,Real] -> {}
solve[erfc[x]==2,x,Real] -> {}
solve[erf[x]==1/2,x,Real] -> UnresolvedSolutionSet[x]
solve[exp[x]==x+2,x,Real] -> {x == -lambertw[-exp[-2]]-2, x == -lambertw[-1, -exp[-2]]-2}
solve[2^x==x,x,Real] -> {}
solve[(4/3)^x==x,x,Real] -> {x == -lambertw[-log[4/3]]/log[4/3], x == -lambertw[-1, -log[4/3]]/log[4/3]}
solve[x^x==1,x,Real] -> {x == 1}
solve[x^x==2,x,Real] -> {x == exp[lambertw[log[2]]]}
```

Equations of the form `exp[p x+q]==c x+d`, and `a^(m x+n)==c x+d` for a positive constant base, are reduced to Lambert W when coefficient reality and nonzero conditions are proved exactly. The transformed argument `z` is compared with `-1/E` by certified ordering so the complete real branch count is classified as zero, one, or two; the coincident `W_0/W_-1` branches at the branch point are merged.

The principal equation `x^x==r` is more delicate on the negative real axis because the value is generally complex and targets inside the unit interval can also receive discrete negative-integer solutions. The current complete classifier therefore handles `r>1` together with the exact cases `r=1,0,-1` and `r<-1`, where negative real roots can be excluded completely. For example, `solve[x^x==1/4,x,Real]` remains `UnresolvedSolutionSet[x]` rather than returning only the positive-real branch.

These proofs are not reused in the Complex domain. For example, a Real nonexistence proof alone does not turn `solve[exp[x]-x+5==0,x,Complex]` into `{}`.

### Absolute-value equations and inequalities

For Real relations, `abs[u]` is handled with its exact range `[0,Infinity)` intact. Ordered inequalities use Real semantics and, once the sign of the right-hand side is proved, may reduce safely to the corresponding polynomial inequality in `u^2`. An equality `abs[u]==a` is split into `u==a` or `u==-a` only under an **explicit Real domain**. In the default Complex domain, `abs[z]==a` generally describes a locus such as a circle, so mmCal keeps it unresolved rather than collapsing it to two points.

```text
solve[abs[x]<2,x]
-> {x in Real if x > -2 && x < 2}

solve[abs[x-1]>=3,x]
-> {x in Real if x <= -2, x in Real if x >= 4}

solve[abs[2x-1]<=3,x]
-> {x in Real if x >= -1 && x <= 2}

solve[abs[x]==2,x,Real]
-> {x == 2, x == -2}

solve[abs[x]==-2,x,Real]
-> {}

solve[abs[x]!=2,x,Real]
-> {x in Real if x != 2 && x != -2}

solve[abs[x]==2,x]
-> UnresolvedSolutionSet[x]
```

### Real exponentials and Lambert W

When the Real-domain solver can prove `a>0` and the exponent real, the principal power `a^u = exp[u log[a]]` is strictly positive. Zero equations therefore close to the empty set without numerical search.

```text
solve[1.1^x == 0,x,Real] -> {}
solve[1.1^x == 0,Real]   -> {}
solve[2^x == 8,x,Real]   -> {x == 3}
solve[2^(2x+1) == 8,x,Real] -> {x == 1}
solve[2^x == -1,x,Real]  -> {}
```

When the right-hand side is independent of the solve variable and the solver can prove `a>0`, `a!=1`, and `r>0`, a constant-base equation `a^u==r` is safely inverted to `u==log[a,r]` and passed to the existing polynomial solver. The rewrite is not used when its branch/domain requirements cannot be proved.

The Lambert W normalization layer directly recognizes `u exp[u]==a`. On the Real domain it enumerates `W_0(a)` under `a>=-1/E` and also `W_-1(a)` under `-1/E<a<0`, without duplicating the branch point. Equations such as `exp[-u]==u` are normalized to the same form. Principal `lambertw[u]==r` is inverted to `u==r exp[r]` only when the real principal range `r>=-1` is proved; this does not guess a general Complex branch family.

```text
solve[x*exp[x]==1,x,Real] -> {x == lambertw[1]}
solve[exp[-x]==x,x,Real] -> {x == lambertw[1]}
solve[x+log[x]==0,x,Real] -> {x == lambertw[1]}
solve[exp[x]+x==0,x,Real] -> {x == -lambertw[1]}
solve[lambertw[x]==1,x] -> {x == E}
solve[cosh[x]==2,x,Real] -> {x == acosh[2], x == -acosh[2]}
```

The existing `a^x==x^2` classifier remains alongside this normalization. Writing `L=log[a]`, one real root always comes from the principal branch; the two negative-argument branches `W_0` / `W_-1` are added only when `|L|<=2/E` is certified. At the branch point they coincide and are not duplicated.

```text
solve[1.1^x == x^2,x,Real]
-> {x == -2lambertw[log[11/10]/2]/log[11/10],
    x == -2lambertw[-log[11/10]/2]/log[11/10],
    x == -2lambertw[-1,-log[11/10]/2]/log[11/10]}

N[solve[1.1^x == x^2,x,Real],20]
-> {x == -0.95548727594562198165,
    x == 1.0513800237472769374,
    x == 95.71683016840522274}
```

General `a^(b x+c)==P(x)`, complete Complex branch families, and inequalities involving Lambert W remain deferred. If branch conditions cannot be proved, the solver keeps conditional branches or an `UnresolvedSolutionSet` rather than guessing.

### Principal radical equation solve

Equations involving `sqrt` / `cbrt` are not solved by merely raising both sides to a power. Candidate generation from the transformed polynomial is separated from an exact range check for the original radical.

For the principal square root,

```text
sqrt[A] == B
    <=> A == B^2 and B lies in the image of the principal sqrt
```

is used. Roots introduced by squaring are therefore rejected by an exact range proof. For exact numbers the principal image `Re(B)>0` or `Re(B)==0 && Im(B)>=0` is checked directly; when `B` is proved real this reduces to `B>=0`. If a symbolic complex parameter would require a disjunctive half-plane condition that the current predicate representation cannot express completely, the solver keeps `UnresolvedSolutionSet` rather than weakening the condition.

```text
solve[sqrt[x]==2,x]       -> {x == 4}
solve[sqrt[x]==-2,x]      -> {}
solve[sqrt[x]==I,x]       -> {x == -1}
solve[sqrt[x+1]==x-1,x]   -> {x == 3}
solve[sqrt[x^2]==2,x]     -> {x == 2, x == -2}
```

`cbrt` is the real cube root, so `cbrt[A]==B` is equivalent to `A==B^3` together with `B in Real`. When the transformed equation is a higher-degree Rational polynomial and an affine real right-hand side proves that every original solution has a real solve variable, the solver uses real algebraic isolation.

```text
solve[cbrt[x]==2,x]         -> {x == 8}
solve[cbrt[x]==-2,x]        -> {x == -8}
solve[cbrt[x+1]==x-1,x]     -> {x == root[{-2, 2, -3, 1}, 1]}
solve[cbrt[x]==a,x]         -> {x == a^3 if a in Real}
```

The Simplifier also applies `cbrt[z]^3 -> z` only when `z` is proved real, preserving the real-domain requirement of `cbrt`.

### Parameterized real solution families for periodic functions

`sin/cos/tan` are not globally injective, so mmCal does not collapse them to one principal inverse. On the real axis it can instead preserve periodicity with an integer formal parameter. The parameter is locally bound by `where k in Integer`; if `k` already occurs in the relation, a fresh name such as `k1` is selected.

```text
solve[sin[x]==0,x,Real]
-> {x==Pi k where k in Integer}

solve[cos[x]==0,x,Real]
-> {x==Pi/2+Pi k where k in Integer}

solve[tan[x]==1,x,Real]
-> {x==Pi/4+Pi k where k in Integer}
```

The first implementation is deliberately limited to one-argument `sin/cos/tan` equations whose argument is affine in the solve variable with an **exact nonzero linear coefficient**. Existing Knowledge checks the real target range, so equations such as `sin[x]==2` reduce to the empty set. Nonlinear arguments and complete Complex-domain periodic families remain unresolved rather than being guessed.

The period follows the current session angle mode. The first periodic Solver does not yet normalize an explicit `Rad` / `Deg` / `Grad` suffix around the entire argument into the affine polynomial matcher, so that form remains unresolved for now.

```text
angleMode[Deg]
solve[sin[x]==0,x,Real]
-> {x==180k where k in Integer}
```

The Complex domain likewise does not fabricate complete solution sets from principal inverses alone.

### Exact algebraic roots: `root` / `AlgebraicNumber`

`root` represents exact real and complex algebraic roots without forcing radical expansions.

```text
root[{a0,a1,...,an},k]
root[{a0,a1,...,an},k,Complex]
```

The coefficient list stores the exact Rational polynomial `a0 + a1 x + ... + an x^n` in ascending power order. The two-argument form denotes the 1-based `k`th distinct real root in increasing order. The three-argument `Complex` form isolates all complex roots and assigns deterministic 1-based indices by increasing `Re(z)+Pi Im(z)`. Because the real and imaginary parts of algebraic roots are algebraic while Pi is transcendental, two distinct algebraic roots cannot share that exact ordering key. Root certification itself does not depend on this approximate ordering procedure: a unique Rational-center/Rational-radius disk is first proven by an exact Rouche test, then ordering is certified.

Defining polynomials normalize exactly to monic square-free form.

```text
root[{4,0,-4,0,1},2]
-> root[{-2,0,1},2]
```

`RealAlgebraicNumber` uses exact Rational Sturm sequences and isolating intervals. `ComplexAlgebraicNumber` stores exact Rational-center isolating disks; numerical root candidates are used only to propose disks that must subsequently pass exact certification. Candidate initialization uses Newton-polygon radius groups when available, while exact Rouché certification remains the sole acceptance criterion. `N[root[...,k],p]` reuses the existing certified interval/disk and locally refines only the selected root.

```text
N[root[{-2,0,1},2],30]
-> 1.41421356237309504880168872421

N[root[{1,0,1},2,Complex],30]
-> I

solve[x^5-x+1==0,x,Real]
-> {x==root[{1,-1,0,0,0,1},1]}

solve[x^5-x+1==0,x]
-> {x==root[{1,-1,0,0,0,1},1,Complex], ...}
```

Existing linear, quadratic, binomial, and Rational-root-deflation solvers remain preferred when they close naturally. Remaining Rational-polynomial equations use Sturm real Root fallback in a Real domain and certified complex Root isolation in the default/Complex domain. The defining-polynomial degree budget is currently 96.

`AlgebraicNumber` unifies Real/Complex Root values and performs bounded exact `+ - * /` arithmetic with other Root values and exact Rational/complex-Rational operands. It forms a resultant candidate polynomial and then uses the operands' isolating intervals/disks to certify exactly one result root before replacing the expression.

```text
root[{-2,0,1},2]*root[{-2,0,1},2]
-> 2

root[{1,0,1},2,Complex]+I
-> 2I
```

Individual Root construction no longer stops at square-free monic normalization when a stronger statement can be proven: **the defining polynomial is reduced to the Rational irreducible factor containing the selected root only when that factor is certified exactly**. The bounded factorization backend currently covers degree at most 16, combining finite-field irreducibility proofs with exact Kronecker factor search. Real roots identify the selected factor by exact Sturm root counts; Complex roots require a unique certified match between the selected isolating disk and the factor's root disks. If the proof fails, the original square-free defining polynomial is retained rather than guessing a minimal polynomial.

```text
root[{6,0,-5,0,1},1]
-> root[{-3,0,1},1]

root[{2,0,3,0,1},1,Complex]
-> root[{2,0,1},1,Complex]
```

`isolateAll` intentionally preserves the original polynomial and global root index so enumeration remains stable for one polynomial; minimal-polynomial canonicalization is applied when an individual `root[...]` is constructed or a Solve result is materialized.

For Root-to-Root field arithmetic, if both operand minimal polynomials are proven irreducible over Q and a candidate `theta=alpha+c beta` has a first exact power dependence of degree `deg(alpha)deg(beta)` whose polynomial is also proven irreducible, theta is accepted as a **primitive element**. Only then is the tensor-product basis converted exactly into the theta power basis, and the result's minimal polynomial is derived from exact linear dependence inside the simple extension. Overlapping extensions, failed certificates, or budget overflow fall back to the previous resultant plus isolating-region re-identification path.

```text
root[{-2,0,1},2]+root[{-3,0,0,1},1]
-> root[{1,-36,12,-6,-6,0,1},2]
```

Once primitive-element reduction certifies a simple extension, the result does not discard that work after formatting it as `root[minpoly,k]`. Internally, `NumberFieldContext` retains the generator minimal polynomial, selected Real/Complex embedding, and reduction of `theta^d` modulo `m(theta)`, while `AlgebraicElement` retains exact Rational coordinates in the power basis. Later same-Context `+ - * /` therefore use coefficient arithmetic modulo `m(theta)` rather than rebuilding a resultant or primitive element. Individual Roots also receive a generator-field representation when Q-irreducibility is proven, so repeated arithmetic on the same Root can remain inside its original simple extension.

```text
root[{1,-1,0,0,0,1},1,Complex]^2
-> root[{-1,1,0,-2,0,1},3,Complex]
```

The public canonical form remains `root[minpoly,k]`; field coordinates are not exposed to formatting or structural equality. `Expr::rebuildCall` preserves the internal Algebraic cache only when a Call's arguments remain structurally unchanged and invalidates it when they change, allowing Simplifier, substitution, and constraint-processing paths to retain the same field lineage safely.

A bounded weak interner is used for `NumberFieldContext`. Independently constructed contexts share one immutable object only when they have the **same embedded generator identity**: the same minimal polynomial, Root domain, and root index. The interner owns no Contexts, stores only `weak_ptr`s, prunes expired entries, and uses a 256-entry LRU bound; misses and eviction affect performance only, never mathematical identity. Roots with different selected embeddings are not merged even when their minimal polynomials agree, and isomorphic fields expressed through different primitive generators or subfield relations are not guessed equivalent.

This lets independently derived elements of the same simple extension enter the pointer-level same-field fast path. For example, the two sides below are built through separate primitive-element reductions but now multiply inside their shared field instead of re-entering a higher-degree construction:

```text
(root[{-2,0,1},2]+root[{-3,0,0,1},1])
*(root[{-2,0,1},2]-root[{-3,0,0,1},1])
-> root[{1,12,-6,1},1]
```

The common-field construction itself is also reused. After primitive-element reduction succeeds for a Root pair, a bounded cache of at most 64 entries retains the compositum `NumberFieldContext` and the power-basis embeddings of both operands. Reversed operand order is recognized, so sequences such as `alpha+beta` followed by `alpha-beta` do not repeat the same tensor-product / primitive-element search. The output field is referenced weakly and entries whose field has expired are discarded.

For a Real field, an `AlgebraicElement` can evaluate its coordinate polynomial over the chosen generator's certified isolating interval using exact Rational interval arithmetic. A Sturm certificate proves that this interval contains exactly one root of the result minimal polynomial, and the number of roots below the interval determines the canonical root index directly. This avoids isolating every real root of the result polynomial and then repeating all-root isolation while constructing `root[minpoly,k]`; the previous isolating-region path remains a fallback when the direct certificate does not resolve the root.

A value carrying a persistent field representation also already has a certified minimal polynomial. Therefore a Real value of degree greater than 1 cannot be Rational, and a Complex value of degree greater than 2 cannot lie in `Q+iQ`; `exactRationalParts` skips its former 192-bit refinement in those cases. These are proof-reuse optimizations only and do not alter canonical output or exact semantics.

Reciprocals inside one `NumberFieldContext` are also reused. On the first miss, the inverse of a power-basis coordinate vector `u` is computed exactly by extended Euclid in `Q[t]/(m)` and stored in a thread-safe per-field LRU capped at 16 entries. Since `inverse(inverse(u)) = u`, both directions of the pair are published together; a hit only copies the cached Rational coefficient vector. Eviction can only cause recomputation and cannot affect the mathematical result. Constant coordinates `{q,0,...}` use the canonical embedding of `Q` and return `{1/q,0,...}` directly without polynomial Euclid. A persistent multiplication-matrix cache was also measured, but on a degree-12 field it reduced a representative multiplication only from about 112 us to 101 us while costing about 228 us to build, so it is deliberately not enabled in the current implementation.

A dedicated development benchmark is available:

```text
mmCal.Benchmarks --algebraic-field [iterations]
```

It measures first-run and warm average timings, in one session, for the compositum-reuse expression corresponding to `(sqrt[2]+cuberoot[3])*(sqrt[2]-cuberoot[3])`, and also microbenchmarks first/warm reciprocal lookup, warm division, and first/warm minimal-polynomial derivation in a degree-12 simple extension. Timings depend on compiler, build configuration, and CPU and are intended for same-environment regression monitoring rather than absolute performance guarantees.

This representation is also connected to exact comparisons. `==` / `!=` compare power-basis coordinates directly inside one `NumberFieldContext`, accept an identical canonical Root identity immediately, and certify distinct root indices of one polynomial or distinct proven irreducible minimal polynomials as unequal. Remaining bounded cases may construct the exact difference through the existing primitive-element/resultant path and test zero from field coordinates or certified root isolation. Proof failure or budget overflow never becomes `False`.

Mathematical `< <= > >=` is defined only for real algebraic values. In one field, the sign of `a-b` is certified by exact Rational interval evaluation of its power-basis polynomial at the chosen real embedding. Across different real fields, the certified isolating intervals are refined until they separate; exact difference construction is available as a fallback. The deterministic `Re(z)+Pi Im(z)` ordering used to enumerate Complex Roots is not a mathematical order, so `< <= > >=` on Complex algebraic values remains unevaluated.

```text
root[{-2,0,1},2] > 1
-> True

root[{-2,0,1},2] != root[{-3,0,1},2]
-> True

root[{1,0,1},1,Complex] < root[{1,0,1},2,Complex]
-> root[{1,0,1},1,Complex] < root[{1,0,1},2,Complex]
```

Expressions that are not syntactically Root values are also bridged when their exact algebraic meaning can be certified. The current bridge covers canonical `root[...]`, exact Rational / exact complex-Rational values, `sqrt[q]` / `cbrt[q]` for safe exact-Rational real cases, `Phi`, and bounded `+ - * /` or small integer powers built from them. Formatting is not forced into Root form; comparisons, domain proofs, and Solve use the shared `AlgebraicNumber` view internally while preserving the original user-visible radical or constant expression.

```text
root[{-2,0,1},2] == sqrt[2]
-> True

root[{-2,0,0,1},1] == cbrt[2]
-> True

Phi == root[{-1,-1,1},2]
-> True

element[sqrt[2]+sqrt[3],Rational]
-> False

solve[x == sqrt[2],x,Rational]
-> {}

solve[x == sqrt[2],x,Real]
-> {x == sqrt[2]}
```

The bridge has node and small-integer-power budgets and falls back to the previous symbolic path when conversion cannot be certified. It never guesses a minimal polynomial from an approximation. Undefined forms such as `0^0` become `Indeterminate` before entering the algebraic bridge; `a^0=1` is used only when exact nonzeroness of the base can be proven.

To bound resultant and primitive-element growth, the current **algebraic-field candidate-degree budget is 16**. Complete arbitrary-degree Q-factorization, general number-field merging/canonicalization across different primitive generators or subfield relations, general reduction of overlapping extensions when the current simple-extension certificate fails, complete cross-context equality for unproven/nonminimal representations, ordering beyond the current refinement/algebraic-construction budgets, and `rootApproximant` remain deferred. The current layer is therefore a persistent-lineage bounded exact algebraic-field backend with exact comparison support, not yet a complete number-field canonicalizer.

---

# 26. Linear algebra

## 26.1 Exact-first linear algebra

Canonical API:

```text
transpose[A]
conjugateTranspose[A]
dot[A,B]
det[A]
inverse[A]
rref[A]
matrixRank[A]
nullSpace[A]
solveLinear[A,b]
luDecomposition[A]
qrDecomposition[A]
svd[A]
conditionNumber[A]
pseudoInverse[A]
leastSquares[A,b]
eigenvalues[A]
eigenvectors[A]
eigensystem[A]
norm[v]
normalize[v]
trace[A]
```

`dot` supports rank-1 and rank-2 Arrays.

```text
dot[{1,2,3},{4,5,6}] -> 32
dot[{{1,2},{3,4}},{5,6}] -> {17, 39}
dot[{5,6},{{1,2},{3,4}}] -> {23, 34}
dot[{{1,2},{3,4}},{{5,6},{7,8}}] -> {{19, 22}, {43, 50}}
```

Array-by-Array `*` is not matrix multiplication. `*` accepts scalar×Array multiplication; matrix multiplication and vector contraction use explicit `dot`. `+/-` are element-wise for equal shapes.

Exact real/Rational matrices clear row denominators and use Bareiss fraction-free elimination on an integer work matrix, avoiding Rational construction at every pivot. Exact complex matrices fall back to the flat `Number` Gaussian backend. Symbolic elimination never guesses a pivot whose nonzero status cannot be proven.

```text
det[{{1,2},{3,4}}] -> -2
inverse[{{1,2},{3,4}}] -> {{-2, 1}, {3/2, -1/2}}
rref[{{1,2},{3,4}}] -> {{1, 0}, {0, 1}}
matrixRank[{{1,2},{2,4}}] -> 1
nullSpace[{{1,2},{2,4}}] -> {{-2, 1}}
solveLinear[{{2,1},{1,-1}},{5,1}] -> {2, 1}
```

`nullSpace[A]` returns a canonical RREF basis by taking free columns in ascending order and setting each corresponding free variable to one. Its result shape is `{nullity, columns}`; full column rank therefore formats as `reshape[{}, {0,n}]` so the vector dimension of the empty basis is not lost. Exact integer/Rational inputs share the Bareiss forward elimination path, exact complex matrices use the Gaussian fallback, and symbolic matrices produce a basis only when pivot nonzero status is provable.
For a matrix that already contains finite-precision elements, both pivot existence and the absence of a pivot in a free column must be certified from the InformationEnclosure before nullity is fixed. Thus `nullSpace[N[{{Pi,1}},12]]` may evaluate when the pivot structure is certified by the declared input information, while a finite-precision zero column is not treated as exact zero merely because its hidden CertifiedEnclosure is a point.

`solveLinear[A,b]` treats `A` as an m×n matrix and `b` as a length-m vector. It returns a length-n vector only when the solution is unique. The matrix need not be square: a consistent overdetermined system is accepted when it has full column rank. Inconsistent systems and systems with free variables are Domain errors; this function does not invent a parametric solution.

Exact integer/Rational matrices clear denominators per row into an integer workspace and automatically choose Bareiss or the 31-bit modular backend according to order, coefficient height, and density. Modular `det` reconstructs a unique integer by CRT through an integer Hadamard bound. Modular `solveLinear` applies rational reconstruction after CRT and returns a candidate only after exact verification against the original integer system; bad primes or unsuccessful reconstruction fall back to Bareiss. A modular inverse backend also exists, but automatic `inverse` remains on Bareiss because Bareiss still wins throughout the currently measured GCC range.

`luDecomposition[A]` currently targets square matrices and returns shape `{3,n,n}` containing `{P,L,U}`, with the convention `P A = L U`. Certified approximate LU uses partial pivoting among provably nonzero candidates, maximizing the certified lower bound of `|pivot|^2`; no epsilon threshold is used. Row pivoting is used, exact Number input remains exact for Rational and complex values, triangular symbolic matrices avoid unnecessary division, and a general symbolic decomposition remains unevaluated when a required pivot cannot be proved nonzero. Prefix indexing extracts each factor.

```text
lu = luDecomposition[A]
at[lu,0] -> P
at[lu,1] -> L
at[lu,2] -> U
```

`qrDecomposition[A]` uses reduced QR for rectangular m×n matrices. With `k=min(m,n)` it returns the general brace `{Q,R}` with `Q:m×k`, `R:k×n`, and `A = Q R`. Equal-shape square factors may be optimized internally to a dense Array, but the user representation is the same. Exact real matrices first reduce each column to a primitive integer direction and use division-free fraction-free projection; square roots are not generated during orthogonalization and are introduced only when Q/R are materialized as final expressions. For full-rank leading columns, a symmetric Bareiss factorization of the Gram matrix (a fraction-free LDLᵀ-equivalent path) reconstructs the orthogonal integer basis; rank-deficient cases fall back to direct fraction-free orthogonalization. The former 3x3 hard cap is removed. Safe upper-triangular/trapezoidal cases retain their fast path. `N[qrDecomposition[A],p]` bypasses exact expansion and dispatches directly to a certified interval Householder backend for both real and complex matrices. QR factor column signs are not mathematically unique, so componentwise sign equality between the exact fraction-free backend and the certified Householder backend is not part of the contract. The contract is `A = Q R` with orthonormal Q columns; the exact backend uses a deterministic orientation induced by its primitive integer directions.

```text
qr = qrDecomposition[A]
at[qr,0] -> Q
at[qr,1] -> R
```

Householder application also has a column-block kernel that processes multiple columns during one row-major scan. On the current no-BLAS BigFloat/interval backend, measurements at orders 8/16/24 did not show a consistent speedup, so automatic blocking is not enabled; the unblocked-equivalent path remains the default and the block kernel/benchmark are retained for later backend optimization.

`svd[A]` returns reduced `{U,S,V}`. For m×n input, `k=min(m,n)`, `U:m×k`, `S:k×k`, and `V:n×k`. Real input satisfies `A = U S Transpose[V]`; complex input satisfies `A = U S conjugateTranspose[V]`. The general numerical backend deliberately does not form `A^H A`: it uses Householder bidiagonalization followed by one-sided Jacobi column orthogonalization. Candidate factors are returned only after interval checks validate reconstruction residual and `U^H U` / `V^H V` orthogonality more strictly than the requested output digits; otherwise guard digits are increased and the calculation is retried. Exact SVD is restricted to natural closed cases such as exact real diagonal matrices. Singular vectors are not unique inside repeated-singular-value subspaces, so the certificate concerns reconstruction and orthogonality rather than a unique componentwise vector.

`conditionNumber[A]` returns the spectral 2-norm condition number `sigma_max/sigma_min`. Proven exact rank deficiency returns `Infinity`; a nonzero rectangular matrix with only one singular value returns `1`, and exact real diagonal matrices return an exact ratio. A general exact nondiagonal matrix is not forced into a large singular-value expression; use `N[...]` to dispatch it to certified SVD. The condition number of an empty matrix is a DomainError.

`pseudoInverse[A]` returns the Moore-Penrose pseudoinverse. Exact numeric matrices use rank factorization `A=FG` and evaluate `A^+=G^H(GG^H)^-1(F^H F)^-1F^H` with exact arithmetic, so rank-deficient Rational and complex matrices remain exact. Zero-by-n and n-by-zero inputs return an empty matrix with transposed shape. For `N[pseudoInverse[A],p]` with exact `A`, the requested-precision path uses certified SVD. A matrix that already contains finite-precision elements is currently kept unevaluated because the SVD backend does not yet certify the full input perturbation through singular subspaces; hidden CertifiedEnclosure points are not used to reconstruct rank or singular values.

`leastSquares[A,b]` returns the minimum-norm least-squares solution `A^+ b`. The length of `b` must equal the row count of `A`. Exact numeric inputs remain exact, including rank-deficient cases. The outer-`N` path for an exact matrix may use certified SVD, while a matrix that already contains finite-precision elements follows the same conservative rule as `pseudoInverse` and remains unevaluated.

```text
conditionNumber[{{3,0},{0,4}}] -> 4/3
conditionNumber[{{1,2},{2,4}}] -> Infinity
pseudoInverse[{{1,2},{2,4}}] -> {{1/25, 2/25}, {2/25, 4/25}}
pseudoInverse[{{I,0},{0,2I}}] -> {{-I, 0}, {0, -I/2}}
leastSquares[{{1,0},{0,1},{1,1}},{1,2,4}] -> {4/3, 7/3}
dimensions[pseudoInverse[zeros[0,3]]] -> {3, 0}
```

`eigenvalues[A]` / `eigenvectors[A]` / `eigensystem[A]` handle eigenvalues, eigenvectors, and the paired result for square matrices. Eigenvectors are returned as **columns**, and `eigensystem[A]` returns `{values,vectors}`. The exact path handles diagonal entries of upper-triangular matrices, the standard basis of diagonal matrices, and exact Number 2x2 matrices with distinct eigenvalues. A nondiagonal repeated-root 2x2 matrix is not given duplicated vectors merely to fill a basis. General `N[...]` uses Complex BigFloat Hessenberg reduction followed by implicit shifted QR to obtain a Schur relation `A Q ≈ Q T`; eigenvectors are then recovered from triangular back substitution. The original certified input intervals are used to audit the Schur relation, `A v ≈ λ v` residuals, and Schur-vector unitarity more strictly than the requested display digits, retrying with more guard digits when necessary. For a general non-normal matrix, individual eigenvalue/eigenvector components can be perturbation-sensitive, so mmCal does not claim that every displayed component is a unique componentwise enclosure of a mathematically distinguished exact value; the certificate concerns the computed Schur/eigenpair relations. Near-multiple or defective cases that do not yield a stable independent eigenvector basis remain unevaluated rather than being guessed.

`conjugateTranspose[A]` computes the Hermitian transpose used by complex SVD and complex orthogonality checks. Rank-1 input is conjugated componentwise; rank-2 input is transposed and conjugated.

`norm` is a Hermitian norm for complex vectors.

```text
norm[{3,4}] -> 5
norm[{3+4I}] -> 5
normalize[{3,4}] -> {3/5, 4/5}
```

### Precision-aware `N`

Matrix backends share the same `ApproximationContext` and certified Expr/interval conversion layer as FFT. Therefore calls such as

```text
N[dot[A,B],100]
N[det[A],100]
N[inverse[A],100]
N[rref[A],100]
N[solveLinear[A,b],100]
N[luDecomposition[A],100]
N[qrDecomposition[A],100]
N[svd[A],100]
N[conditionNumber[A],100]
N[pseudoInverse[A],100]
N[leastSquares[A,b],100]
N[eigenvalues[A],100]
N[eigensystem[A],100]
N[norm[v],100]
```

can pass the requested precision directly to the BigFloat/interval backend instead of first constructing a potentially huge exact intermediate expression. When an exact matrix is sent into an outer `N`, its input enclosure can be recomputed at the requested working precision. When matrix elements already contain `DecimalApproximation` / `ComplexDecimalApproximation` leaves, value computation uses the CertifiedEnclosure while zero/nonzero, pivot, and rank decisions use the InformationEnclosure, so undeclared guard information cannot reappear. `solveLinear`, `inverse`, `rref`, `matrixRank`, and `nullSpace` return results only when InformationEnclosures certify the required pivot structure; uncertain cases are not completed with epsilon thresholds or hidden point values. Continuous quantities such as `det`, `dot`, and `norm` propagate Certified/Information enclosures in parallel so cancellation naturally reduces output Precision. The current `luDecomposition`, `qrDecomposition`, `svd`, `conditionNumber`, `pseudoInverse`, `leastSquares`, and `eigen*` perturbation certificates are not yet defined for matrices that already contain finite-precision leaves, so those inputs remain conservatively unevaluated. Exact matrices evaluated through `N[...,p]` continue to use the certified numerical backends.

## 26.2 Vector helpers and compatibility aliases

The vector-oriented API contains both independent canonical functions and compatibility aliases.

Canonical functions:

```text
madd
vadd vsub vscalar
vcross
inner outer
vproject vangle
vmanhattan veuclidean
vreflect vreflect_axis
vsum
grad divergence curl laplacian jacobian hessian
```

Compatibility aliases:

```text
matmul mmul vdot -> dot
rank mrank       -> matrixRank
mget             -> at
vnorm vlength    -> norm
vnormalize vunit -> normalize
vdistance distance -> veuclidean
cross              -> vcross
projection         -> vproject
gradient           -> grad
singularValueDecomposition -> svd
```

Aliases have no separate algorithm; they resolve to the same `BuiltinId`. Canonical vector helpers such as `vadd` are public APIs in their own right, not compatibility names. `dot` remains a bilinear contraction, while `inner[a,b]` is the Hermitian inner product that conjugates its first argument and therefore matches `norm[v]`. `projection[a,b]` / `vproject[a,b]` uses `b inner[b,a]/inner[b,b]`.

Vector-calculus operators use explicit Cartesian coordinates.

```text
grad[f,{x,y,z}]
divergence[{P,Q,R},{x,y,z}]
curl[{P,Q,R},{x,y,z}]
laplacian[f,{x,y,z}]
jacobian[{f1,f2,...},{x1,x2,...}]
hessian[f,{x1,x2,...}]
```

The coordinate specification must be a rank-1 Array of distinct symbols. `curl` is restricted to three-dimensional Cartesian fields. No curvilinear scale factors or metric are assumed implicitly.

---

# 27. Signal processing

```text
dft[v]
fft[v]
ifft[v]
convolve[a,b]
```

Fourier phase is explicitly evaluated in radians and does not depend on the session's default angle unit.

```text
dft[{1,2,3,4}]
-> {10,-2+2I,-2,-2-2I}

ifft[fft[{1+I,2-I,3+2I,4-3I}]]
-> {1+I,2-I,3+2I,4-3I}

convolve[{1,2},{3,4}]
-> {3,10,8}
```

For exact inputs, power-of-two FFTs use radix-2 Cooley–Tukey and retain the existing public forward representation. For power-of-two `ifft` calls within the current degree budget (lengths 16–128), when the input contains FFT-generated root-of-unity structure, mmCal re-embeds that structure into Rational power-basis coordinates of `Q[t]/Phi_n(t)` and performs the radix-2 inverse there. For power-of-two cyclotomic fields, `Phi_(2^m)(t)=t^(2^(m-1))+1` reduces twiddle multiplication to coefficient shifts with sign changes, avoiding the huge symbolic expansion that previously appeared in `ifft[fft[v]]`. Purely numeric spectra or expressions whose membership cannot be proved fall back to the previous path. For non-power-of-two lengths of at least 5, inputs that can be certified as exact Rational/Gaussian Rational values or expressions in the same cyclotomic quotient are transformed in the same coordinate representation. Gaussian Rational inputs extend the conductor to `lcm(n,4)` when needed so that `I` lies in the same cyclotomic field. If the current cyclotomic-degree budget of 64 is exceeded, or symbolic inputs cannot be proven to lie in the quotient field, the legacy generic exact DFT remains the fallback. Ordinary `fft[...]` remains exact-first and never silently converts to machine `double`.

`N[fft[v],p]` does not first expand the full exact Fourier expression. `N` propagates the requested precision into the FFT call, which performs butterflies directly on certified `ComplexInterval`/BigFloat endpoints and returns decimal components only after their requested rounding is proven unique. For approximate operands, CertifiedEnclosure and InformationEnclosure are transformed independently through the same FFT, so cancellation cannot resurrect hidden guard digits. A difference smaller than the information carried by finite-precision inputs therefore becomes a low-Precision zero-centered result rather than a high-precision tiny value.

The approximate path uses radix-2 for power-of-two sizes and Bluestein convolution for sufficiently large non-power-of-two sizes. Small non-power-of-two transforms keep direct DFT because its constant factor wins there. The current conservative policy uses direct evaluation below 384 points, based on a forced-algorithm GCC sweep together with the MSVC `--full` results. This is not a mathematical boundary; rerun `mmCal.Benchmarks --fft-threshold 1` for the current compiler and machine.

---

# 28. Random numbers

Random-number functions are stateful built-ins.
The RNG state is independent for each `KernelSession`.

## 28.1 Seed

```text
randSeed[42]
-> 42
```

Using the same seed restores the same sequence.

```text
randSeed[42]
a := rand[]
randSeed[42]
rand[] == a
-> True
```

`randSeed[]` reseeds from entropy and returns the integer seed that can be used to reproduce the sequence.

## 28.2 Uniform real

```text
rand[]
rand[hi]
rand[lo,hi]
```

`rand[]` returns an exact Rational on a 53-bit dyadic lattice.
For example, the first sample for seed 42 is:

```text
randSeed[42]
rand[]
-> 227930101193189/1125899906842624
```

Ranges:

```text
rand[]       : [0,1)
rand[hi]     : [0,hi), hi >= 0
rand[lo,hi]  : [lo,hi), lo <= hi
```

## 28.3 Integer

```text
randint[]
randint[a]
randint[a,b]
```

```text
randint[]    : {0,1}
randint[5]   : [0,5] inclusive
randint[-5]  : [-5,0] inclusive
randint[1,6] : [1,6] inclusive
```

Arbitrary `BigInt` ranges are supported. Rejection sampling is used to avoid modulo bias.

## 28.4 Choice

```text
choice[2,3,5,7]
choice[{2,3,5,7}]
```

Selects one element from either a rank-1 Array or a variadic list.

## 28.5 Normal

```text
randn[]
randn[mu]
randn[mu,sigma]
```

Box–Muller is applied to exact dyadic uniform samples and returns a symbolic expression.

```text
randn[5,0]
-> 5

N[randn[],8]
-> example: -0.80379286
```

`randn` uses explicit radians and does not depend on the session's default angle unit.

**Not suitable for cryptographic use.**

---

# 29. Major aliases

| Alias                    | Canonical    |
| ------------------------ | ------------ |
| `pow`                    | `Power`      |
| `fact`                   | `Factorial`  |
| `fract`                  | `frac`       |
| `ln`                     | `log`        |
| `real`                   | `re`         |
| `imag`                   | `im`         |
| `mag`                    | `abs`        |
| `unit`, `csgn`           | `sign`       |
| `rect`                   | `polar`      |
| `ave`                    | `mean`       |
| `matmul`, `mmul`, `vdot` | `dot`        |
| `mtranspose`             | `transpose`  |
| `mget`                   | `at`         |
| `singularValueDecomposition` | `svd` |
| `mdet`                   | `det`        |
| `minverse`               | `inverse`    |
| `rank`, `mrank`          | `matrixRank` |
| `mtrace`                 | `trace`      |
| `mrows`                  | `rows`       |
| `mcols`                  | `cols`       |
| `mdiag`                  | `diag`       |
| `vnorm`, `vlength`       | `norm`       |
| `vdistance`              | `veuclidean` |
| `vnormalize`, `vunit`    | `normalize`  |

Aliases are not separate implementations; they resolve to the same `BuiltinId`. Mathematical metadata and Solver rules are therefore not duplicated.

In mmCal 1.5.0, capitalized aliases added only for Mathematica compatibility (`Sin`, `ArcTan`, `Integrate`, `Solve`, etc.) were removed. Lowercase canonical names are the default for mathematical functions. `D`, `N`, `In`, `Out`, `Exit`, `Clear`, `Defs`, and `UnDef` remain as intentional proper names for symbolic and Kernel operations. If compatibility syntax becomes necessary, it should be implemented as a separate import/compatibility layer rather than by adding aliases to the default namespace.

---

# 30. Current source-callable function list

The current development tree contains **274 registered builtin/alias names / 254 source-callable names**. Internal heads are not included in the source-callable count.

```text
Clear, D, Defs, DtoG, DtoR, Exit, GtoD, GtoR, In, N,
Out, RtoD, RtoG, UnDef, abs, accuracy, acos, acosh, angleMode, arg,
arrayRank, asin, asinh, at, atan, atan2, atanh, ave, beta, betaln, binom, cbrt, cases,
ceil, choice, cis, collect, cols, comb, conditionNumber, conj, conjugateTranspose, convolve, corr, corrspearman,
cos, cosc, cosh, cot, coth, cov, cross, csc, csch, csgn, curl, cv,
det, dft, diag, digamma, diff, dimensions, distance, divergence, dot, eigenvalues, eigenvectors, eigensystem, element, erf, erfc, exp, explain, expand, expc,
Ei, Si, Ci, li, polylog, fresnelc, fresnels, hypergeometric1F1, hypergeometric2F1, ellipticF, ellipticE, ellipticPi,
expm1, fact, factor, factorint, fallingfact, fft, fib, floor, frac, fract, fullSimplify,
gamma, gcd, geomean, grad, gradient, groebnerBasis, harmmean, hessian, hypot, ibeta, identity, if, ifft, im, imag, inner,
integrate, inverse, iqr, isprime, jacobian, kurtp, kurts, laplacian, lcm, leastSquares, length, lgamma, lambertw, limit, ln, log,
log10, log1p, log2, mad, madR, madd, mag, map, matmul, max, mcols,
mdet, mdiag, matrixRank, mean, median, mget, min, minverse, mmul, mod, mode,
luDecomposition, mrank, mrows, mtrace, mtranspose, nextpow2, nextprime, nintegrate, norm, normal, toNormal, normalize, nullSpace, percentile, percentrank, perm, polar, prevprime,
outer, polynomialReduce, pow, precision, prod, projection, pseudoInverse, quantile, quotient, rand, randSeed, randint, randn, range, rank,
qrDecomposition, rationalize, re, real, rect, rem, reshape, risingfact, rms, root, round, rows, rref,
sec, sech, series, sign, simplify, sin, sinc, sinh, sinhc, skew, solve, solveLinear,
singularValueDecomposition, sqrt, stddev, stddevs, stderr, sum, svd, table, tan, tanc, tanh, tanhc, trace,
totient, transpose, trigamma, trimmean, trunc, unit, vadd, vangle, var, vars, vcross, vdistance,
vdot, veuclidean, vlength, vmanhattan, vnorm, vnormalize, vproject, vreflect, vreflect_axis, vscalar,
vsub, vsum, vunit, winsor, winsorR, zeros, zscore, zeta,
bitand, bitor, bitxor, bitnot, bitshiftl, bitshiftr, bitlength, bitcount, bitget, fma, clamp, proj
```

---

# 31. Error / Warning policy

Main `CalcError` categories:

- Syntax
- Domain
- Type
- Overflow
- Name
- Evaluation
- Internal
- ResourceLimit

When evaluation itself succeeds but an algorithmic built-in cannot complete its work, a Warning is returned separately from the result Expr. This applies to `D`, `solve`, `solveLinear`, `N`, `rref`, `matrixRank`, `nullSpace`, `luDecomposition`, `qrDecomposition`, `svd`, `conditionNumber`, `pseudoInverse`, `leastSquares`, `eigenvalues`, `eigenvectors`, `eigensystem` (including the compatibility alias `rank`), and `integrate`, as well as cases where `precision/accuracy/rationalize` retain unsupported input unevaluated.

Info diagnostics are used for normal state-change notifications. Currently this includes variable and function redefinition notices.

```text
D[abs[x],x]
WARN: D could not fully evaluate the derivative; unevaluated D[...] remains
-> D[abs[x], x]
```

An expression such as `sin[x]` that is correctly retained symbolically does not produce a Warning.

Mathematical singular values are preserved as explicit exceptional values such as `ComplexInfinity` or `Indeterminate` when their meaning is defined. Domain violations, type errors, and other failures that cannot be represented as values remain Errors.

Examples:

```text
1/0
-> ComplexInfinity

0/0
-> Indeterminate

log[0]
-> DomainError

gamma[-2]
-> DomainError

randint[5,1]
-> DomainError
```

The Parser/Evaluator retains source spans and documents, allowing Errors that propagate through user-function definitions to include call traces.

---

# 32. Performance policy

Exact and certified computation is substantially more expensive than native `double`.
Historical microbenchmarks have shown costs ranging from roughly 100× to more than 100,000× that of `double`, depending on the operation.

Even so, operations taking tens of microseconds to a few milliseconds remain practical in an interactive CLI, so ordinary semantics are not lowered to `double` merely for speed.

Implemented optimizations include:

- Number real-real fast path
- Removal of redundant GCD operations in Rational multiplication/division
- Minimization of reduction ranges in Rational addition
- Elimination of repeated Simplifier structural-key computation
- Indexing of like terms in Add
- Newton method for BigInt cube root
- Balanced product tree for factorial
- Certified Log range reduction / shared log(2) enclosure
- Radix-2 FFT

Future features such as `for/Plot`, which may evaluate expressions thousands or millions of times, are expected to use an explicit Machine evaluator separate from Exact/Certified semantics.

---

# 33. Major currently unimplemented / deferred features

Representative items:

- General parametric linear systems
- `hilbert` (legacy naming/specification still to be confirmed)
- Engineering functions, financial functions, and unit conversion
- Legacy colon commands such as `:defs`, `:unset`, and `:undef` (`Defs[]/UnDef[]` function forms are implemented). `:angle` has been replaced by `angleMode[]`; `:help` / `:fix` / `:layout` / `:status` are frontend commands
- `for`, `plot`
- General Machine/double evaluation mode

Further candidates:

complete arbitrary-degree Q-factorization, general number-field merging/canonicalization across different primitive generators, complete cross-context comparison for unproven representations, `rootApproximant`, more general parameterized solution families, and a Machine evaluator / `for` / `plot`.

---

# 34. Current CLI

The CLI separates mathematical Kernel state from frontend presentation state. The following startup options are available:

```text
mmCal --fix 16 --angle deg
mmCal --angle rad
mmCal --angle grad --fix 8
mmCal --layout multi
mmCal --eval "expand[(x+1)^3]"
mmCal --batch < expressions.txt
```

- `--fix n`: Set the startup limit on decimal display digits. Internal values are unchanged and unnecessary trailing zeros are omitted
- `--angle deg|rad|grad`: Set the default angle unit at startup
- `--layout auto|single|multi`: Set interactive REPL composition only; it cannot be combined with non-interactive modes
- `--eval expr`: Evaluate one expression non-interactively
- `--batch`: Evaluate standard input one expression per line in one session.
- `--help`, `-h`: Show usage

`--layout` is presentation-only and cannot be combined with `--eval` or `--batch`. The default `auto` uses terminal width and expression structure on a TTY, but falls back to canonical one-line output through pipes or redirects. `single` always uses one line; `multi` structurally expands supported compound results. Kernel Expr values, history, and canonical `formatExpr()` output are unchanged.

Automation modes emit no banner, prompt, `Out[...]` label, or farewell. Successful values go to stdout and warnings/errors to stderr. Exit codes are `0` for success, `2` for argument errors, `3` for `SyntaxError` / `ResourceLimitError`, `4` for evaluation errors, and `5` for `InternalError`. Batch processing continues after line errors and returns the greatest code observed. `--eval` and `--batch` are mutually exclusive.

During an active top-level evaluation, the CLI connects an `EvaluationCancellationToken` to console interruption. Windows handles Ctrl-C / Ctrl-Break and POSIX handles SIGINT. A long-running kernel stops when it next polls cancellation and reports `ResourceLimitError: Evaluation cancelled by frontend`; non-interactive `--eval` therefore returns exit code `3`. Cancellation is cooperative rather than an asynchronous forced stop of every subsystem.

The Lexer/Parser has independent budgets for tokens, AST nodes, nesting, operator chains, numeric-literal digits, function arguments, and Array elements. Exceeding a budget stops before AST lowering with a source-spanned `ResourceLimitError`.

```text
In [1]> 1/3
Out[1]> 1/3
```

- Parse/evaluate one line at a time
- Prompts are fixed as `In [n]>` / `Out[n]>` with no extra spaces
- `Exit[]` is the mathematical function for ending a session. CLI compatibility commands `:quit` / `:exit` are also accepted; bare `exit` / `quit` receive no special treatment
- `Clear[]`: Remove user definitions and all history, resetting the next input number to 1
- `Defs[]`, `UnDef[...]`: Inspect and remove user definitions
- History references `@`, `%`, `%%`, ... together with signed-index reevaluating `In [n]` and snapshot `Out[n]`

## 34.1 `:help`

```text
:help
:help sin
:help functions
:help constants
:help Pi
```

With no argument, `:help` prints a concise REPL-command summary. `:help function` uses `BuiltinRegistry` as the source of truth for canonical names, callable aliases, and arity. A separate user-facing catalog supplies a curated description, explicit input rules, notes where necessary, and one or more examples for **every** registered source-callable builtin. Multi-form operations such as `D`, `integrate`, `root`, `qrDecomposition`, `svd`, and `solve` display multiple usage lines and examples rather than an arity-only placeholder. The internal test suite enumerates `BuiltinRegistry::sourceFunctionNames()` and rejects a missing detailed entry.

`:help Pi` and the `:help constants` index cover protected mathematical constants, Boolean values, numeric domains, `Infinity`, and the `Rad` / `Deg` / `Grad` unit symbols. `:help functions` prints a sorted index of callable canonical names and aliases. An unknown topic remains a successful frontend query and may offer one deterministic nearby topic; adjacent transpositions such as `sdv` suggest `svd`, while common discoverability shorthand such as `qr` suggests the callable `qrDecomposition` without registering `qr` as a new function alias.

Help lookup is handled before parsing or evaluation. Successful, constant, and unknown lookups therefore neither advance `In[n]` nor enter history; an unknown name also points to the function and constant indexes.

## 34.2 `:fix` — presentation-only decimal display

```text
:fix 16
Display: Fixed(16)

In [1]> 1/3
Out[1]> 0.3333333333333333

:fix off
Display: Exact

In [2]> Out[1]
Out[2]> 1/3
```

`:fix n` changes **display only**, rounding to at most `n` digits after the decimal point. Unnecessary trailing zeros are omitted, so for example `31/10` is displayed as `3.1` under `:fix 5`. Stored Expr values, `Out[n]`, and the semantics of `precision/accuracy` are unchanged. This is not a conversion from exact values to Machine/double. The current range for `n` is 0..1000. Entering `:fix` with no argument displays the current mode.

An entire expression that can be certified numerically is approximated only for display. Symbolic expressions containing free variables retain exact notation.

## 34.3 `:layout` — interactive REPL composition

```text
:layout
Layout: Auto

:layout single
Layout: Single

:layout multi
Layout: Multi
```

`:layout` changes **presentation composition only** in the ordinary REPL. `single` preserves the canonical one-line representation. `multi` structurally expands Arrays/Lists/`cases`/finite or conditional solution sets. `auto` chooses from terminal width and expression structure on a TTY, while pipe or redirect output falls back to one line. It does not alter the stored Expr, `Out[n]`, the reparsable canonical formatter, or automation output.

For example, `multi` can display:

```text
Out[1]> {
          {1, 2},
          {3, 4}
        }

Out[2]> cases[
          x^2 if x >= 0;
          -x if x < 0;
          0
        ]
```

With no argument, `:layout` reports the current mode.

## 34.4 `:status`

```text
:status
Angle: Rad
Display: Exact
Layout: Auto
Evaluation: Exact-first
Definitions: 0
History: 0
```

`:status` is also a CLI command and is not stored in history. The design boundary is maintained: mathematical state changes use Kernel functions such as `angleMode[...]`, while frontend queries and presentation changes use CLI commands such as `:help`, `:fix`, and `:layout`.

## 34.5 `:quit` / `:exit`

```text
:quit
:exit
```

Both commands terminate the current CLI session successfully. They are compatibility commands with the same purpose as the expression-level `Exit[]`; they do not parse/evaluate an expression or consume history. Bare `quit` / `exit` are not special. `--batch` also recognizes these commands and stops successfully at that line.

## 34.6 Console title

As auxiliary information, the title is updated to forms such as:

```text
mmCal <version> - Rad - Exact - Layout(Auto)
mmCal <version> - Deg - Fixed(16) - Layout(Multi)
```

- Windows: `SetConsoleTitleA`
- Linux/macOS: ANSI OSC title sequence only when attached to a TTY
- Other platforms: no-op

Failure to change the title is never treated as a calculation Error. `:status` is authoritative for state inspection; terminal software overriding the title has no effect on semantics.

## 34.7 Canonical formatter

Ordinary `Out[n]` uses compact, reparsable mathematical notation rather than an AST dump.

```text
x^2+sin[x]
A-B+C
2(x+sqrt[x])sqrt[x+sqrt[x]]/3
```

- Do not print `+ -` / `+-`; negative terms are rendered with `-`
- Addition/subtraction forms such as `A-(B-C)` may be flattened for display to `A-B+C`
- Do not add unnecessary spaces around `+`, `-`, `*`, `/`, or `^`
- Render relation operators `==`, `!=`, `<`, `<=`, `>`, `>=` with one space on each side for readability
- Implicit multiplication is concatenated only when lexically safe (`2x`, `2sqrt[x]`). Forms that would collide with exponent notation, such as `2exp[x]` or `2E`, are rendered explicitly as `2*exp[x]`, `2*E`
- Preserve necessary spacing when adjacent identifiers would merge into another token (`I Pi`, `x y`)
- When juxtaposition would be ambiguous, such as adjacent numeric tokens, use explicit `*` rather than whitespace
- Preserve precedence and associativity, and regression-test that format → parse → format does not change meaning

The canonical formatter remains the one-line serialization contract. Interactive `:layout` composition is a separate presentation layer and does not alter this representation. Debug/full-form display of internal structure remains a separate future feature.

---

# 35. Major implementation layers

```text
Lexer / Parser
    ↓
Lowerer
    ↓
Expr AST
    ↓
Evaluator explicit task stack
    ↓
Builtin / user function / symbolic operation
    ↓
Simplifier + MathRegistry + KnowledgeContext
    ↓
Exact result
      or
CertifiedEvaluator -> interval -> DecimalApproximation
```

Major responsibilities:

- `SymbolTable`: intern / identity
- `SymbolRegistry`: protected symbol / constant / domain name
- `BuiltinRegistry`: name / alias / arity / Hold attributes
- `MathRegistry`: domain / parity / branch / definedness
- `ValueFacts`: conservative numeric-domain/sign inference
- `KnowledgeContext`: permanent facts + assumptions
- `Simplifier`: safe local rewriting
- `FullSimplifier`: bounded candidate search
- `CertifiedEvaluator`: interval evaluation of the complete expression
- `SolutionSet`: Solver solution-set representation
- `RandomEngine`: session-local stateful PRNG

This separation is maintained so that adding a function does not require duplicating the same mathematical knowledge independently across the Solver, Simplifier, and numerical backend.

## CertifiedEvaluator safety limit

Because CertifiedEvaluator recursively evaluates expressions as intervals, pathologically deep ASTs are treated as unsupported before reaching an OS stack overflow. The current depth limit is 96 levels. This is a safety guard on AST nesting depth, not on the number of terms in ordinary n-ary `Add` / `Multiply` expressions.

`nintegrate` preflights the integrand over the entire interval before constructing higher derivatives, allowing obvious singularities to be rejected as DomainError first.
