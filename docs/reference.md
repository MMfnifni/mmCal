# mmCal Specification and Function Reference

This document is the detailed specification of **mmCal as implemented**.
For a user-oriented introduction, see the root-level `README.md`.
This document changes with each version; retrieve older versions from the Git history when needed.

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
-> 3.141592653589793238462643383280
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

For `N[expr,n]`, mmCal confirms that both endpoints of an interval containing the true value round to the same `n`-digit result before returning a `DecimalApproximation`.

A `DecimalApproximation` is not merely a display string. It retains the requested digit count, provenance (exact input / certified interval), the displayed decimal value as an exact Rational, and an exact Rational enclosure containing the true value. `ComplexDecimalApproximation` similarly retains metadata for the real and imaginary components. `precision`, `accuracy`, and `rationalize` use this metadata directly rather than reparsing the display string and guessing its quality.

```text
N[sqrt[2],30]
-> 1.414213562373095048801688724210
```

The evaluator does not rely solely on heuristic stopping conditions such as "the difference became sufficiently small."

---

# 3. Predefined symbols

Current protected predefined symbols:

| Name | Meaning |
|---|---|
| `Pi` | Circle constant. Exact transcendental constant |
| `E` | Base of the natural logarithm. Exact transcendental constant |
| `Phi` | Golden ratio. Exact algebraic constant |
| `I` | Imaginary unit |
| `True`, `False` | Boolean values |
| `Integer` | Integer domain |
| `Rational` | Rational domain |
| `Real` | Real domain |
| `Complex` | Complex domain |
| `Infinity` | Infinite-precision sentinel returned by `precision/accuracy` for exact values. Automatic extended-real arithmetic simplification remains limited |

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

Internally, arrays are represented as `ArrayExpr` values containing shape information and flattened elements.

## 4.3 Variables and user-defined functions

```text
x := 3
-> 3

f(t) := t^2 + 1
f(4)
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
@       // In[-1]
@@      // In[-2]
@@@     // In[-3]
%       // Out[-1]
%%      // Out[-2]
%%%     // Out[-3]
```

The formal interface is `In[n]` / `Out[n]`. `n > 0` is an absolute prompt index, `n < 0` is relative, and `n = 0` raises TypeError.

```text
In[1]
Out[1]
In[-1]
Out[-1]
```

`In[n]` retrieves the lowered Expr for the target input and then **evaluates it normally in the current session environment**. Positive `In[n]` uses an absolute input number. Negative `In[-n]` counts previous input slots while excluding the input currently being evaluated. Therefore `@` / `@@` / `@@@` / ... mean `In[-1]` / `In[-2]` / `In[-3]` / ... respectively, with no fixed shorthand depth limit. An input that reached parse/lower but failed during evaluation can therefore be retried through `In[-1]`; a slot that never produced a lowered Expr is unavailable.

`Out[n]` returns a stored result snapshot without reevaluation. Positive `Out[n]` is indexed by the absolute input number, while negative `Out[-n]` counts **successful outputs only** from the most recent one. Therefore `% == Out[-1]`, `%% == Out[-2]`, `%%% == Out[-3]`, ... remain true even when failed evaluations occur between successful outputs. Repeated `%` also has no fixed shorthand depth limit.

```text
In[1]> fft[{1,2,3}]
Out[1]> {6, ...}
In[2]> N[@,30]
Out[2]> {6, -1.500000000000000000000000000000+0.866025403784...I, ...}
```

Here `N[@,30]` does not merely approximate the stored `Out[1]`; it reevaluates the `fft[...]` from `In[1]` inside a 30-digit approximation context.

Absolute-reference example:

```text
In[1]> 1+1
Out[1]> 2
In[2]> 2+2
Out[2]> 4
In[3]> In[1]+In[2]
Out[3]> 6
In[4]> In[3]
Out[4]> 6
```

Thus, `In[n]` means "paste the previous input back into the current environment and execute it again." If the previous input contains variable references, assignment, randomness, or other stateful behavior, the current definitions and RNG state are used. This is deliberately separate from displaying a raw historical AST. A positive absolute reference to the input currently being evaluated is forbidden to prevent direct self-recursion.

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
Bitwise operators themselves are not yet implemented.

---

# 5. Angle semantics

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

# 6. Basic arithmetic and algebra

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

# 7. Basic mathematical and complex functions

| Function | Description | Example |
|---|---|---|
| `sqrt[x]` | Principal square root | `sqrt[-4] -> 2I` |
| `cbrt[x]` | Real cube root, real domain | `cbrt[-8] -> -2` |
| `abs[z]` | Absolute value / complex magnitude | `abs[3+4I] -> 5` |
| `sign[z]` | Real sign / complex `z/abs[z]` | `sign[3+4I] -> 3/5+4/5I` |
| `re[z]` | Real part | `re[3+4I] -> 3` |
| `im[z]` | Imaginary part | `im[3+4I] -> 4` |
| `conj[z]` | Complex conjugate | `conj[3+4I] -> 3-4I` |
| `arg[z]` | Principal argument | `arg[-1] -> Pi Rad` |
| `hypot[x,y]` | Exact `sqrt[x^2+y^2]` | `hypot[3,4] -> 5` |
| `cis[x]` | `cos[x]+I sin[x]` | `cis[Pi/3] -> 1/2 + I sqrt[3]/2` |
| `polar[r,t]` | `r cis[t]` | `polar[2,Pi/3]` |
| `nextpow2[x]` | Smallest `n` such that `2^n >= x` | `nextpow2[9] -> 4` |

Compatibility aliases:

```text
real -> re
imag -> im
mag  -> abs
unit,csgn -> sign
rect -> polar
```

---

# 8. Exponential and logarithmic functions

| Function | Semantics |
|---|---|
| `exp[x]` | Entire complex exponential |
| `log[x]` | Principal natural logarithm |
| `log[b,x]` | Principal `Log[x]/Log[b]` |
| `log2[x]` | `log[2,x]` frontend |
| `log10[x]` | `log[10,x]` frontend |
| `expm1[x]` | Stable evaluation of `exp[x]-1` near zero |
| `log1p[x]` | Stable evaluation of `log[1+x]` near zero |

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

# 9. Trigonometric functions

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
atan2[1,-1] -> 3 Pi / 4
```

Poles of `tan`, `sec`, `cot`, and `csc` are treated as definedness conditions. Finite values are not fabricated at poles.

---

# 10. Hyperbolic functions

Implemented:

```text
sinh cosh tanh
asinh acosh atanh
csch sech coth
```

Inverse hyperbolic functions with principal complex branches carry branch metadata in `MathRegistry`.

---

# 11. Cardinal / stable elementary functions

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

Trigonometric cardinal functions normalize angle expressions to radian quantities before forming the ratio.

```text
sinc[Pi/2]
-> 2 / Pi

sinc[90 Deg]
-> 2 / Pi
```

---

# 12. Rounding and integer utilities

```text
floor ceil trunc round frac
gcd lcm mod rem quotient
```

Examples:

```text
floor[-3/2] -> -2
ceil[-3/2]  -> -1
trunc[-3/2] -> -1
round[5/2]  -> 2       // nearest-even
frac[-3/2]  -> 1/2

gcd[84,126,210] -> 42
lcm[6,8,9] -> 72

quotient[-5,3] -> -1
rem[-5,3]      -> -2
mod[-5,3]      -> 1
```

`mod` corresponds to a floor quotient, while `rem` corresponds to a truncate-toward-zero quotient.

The current `round` accepts one argument only. The legacy `round[x,n]` form has not yet been restored.

---

# 13. Combinatorics and lightweight number theory

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

`isprime/nextprime/prevprime/factorint/totient` are tracked as not yet implemented.

---

# 14. Special functions

## 14.1 Gamma / LogGamma

```text
gamma[5]
-> 24

gamma[1/2]
-> sqrt[Pi]

gamma[-1/2]
-> -2 sqrt[Pi]

N[gamma[1/3],20]
-> 2.67893853470774763366
```

General real arguments use Stirling–Bernoulli with a rigorous remainder bound; negative real arguments use reflection.
Non-positive integer poles produce DomainError.

`lgamma[x]` currently means **`log[abs[gamma[x]]]` on the real axis**. It is kept separate from complex `LogGamma`.

## 14.2 Erf

```text
erf[x]
erfc[x]
```

```text
erf[0]  -> 0
erfc[0] -> 1
N[erf[1],20] -> 0.84270079294971486934
```

## 14.3 Beta

Currently restricted to positive real arguments.

```text
beta[2,3] -> 1/12
beta[1/2,1/2] -> Pi
betaln[1/2,1/2] -> log[Pi]
```

The implementation does not unconditionally expand to a general Gamma ratio when doing so could break pole cancellation.

## 14.4 Generalized factorial family

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

## 14.5 Fresnel C / S

mmCal uses the standard Fresnel integrals corresponding to

```text
fresnelc[x] = integral_0^x cos[Pi t^2/2] dt
fresnels[x] = integral_0^x sin[Pi t^2/2] dt
```

as entire odd functions. Arguments that do not close exactly remain symbolic, while `N` performs certified real evaluation.

```text
fresnelc[0] -> 0
fresnels[0] -> 0
N[fresnelc[1],20] -> 0.77989340037682282947
N[fresnels[1],20] -> 0.43825914739035476608

D[fresnelc[x],x] -> cos[Pi x^2/2 Rad]
D[fresnels[x],x] -> sin[Pi x^2/2 Rad]
```

`Rad` is explicit in the derivatives because the Fresnel definitions themselves must not depend on the session's default angle unit.

## 14.6 Confluent hypergeometric 1F1

Kummer's confluent hypergeometric function is written as

```text
hypergeometric1F1[a,b,z]
```

It is entire in `z`; in general `b = 0,-1,-2,...` is a parameter pole, so those cases are not unconditionally simplified to finite values. Exact evaluation currently handles safely terminating series, `z=0`, `a=b`, and related closed cases. The certified real `N` backend currently accepts exact Rational `a,b,z`.

```text
hypergeometric1F1[0,3,2] -> 1
hypergeometric1F1[-2,3,2] -> 0
hypergeometric1F1[2,2,1] -> E
N[hypergeometric1F1[1/6,7/6,1],20]
-> 1.19206880798188830082
```

When the parameters do not depend on the differentiation variable,

```text
D[hypergeometric1F1[a,b,z],z]
= a hypergeometric1F1[a+1,b+1,z]/b
```

is used. For integration, mmCal prefers the 1F1 form when an upper-incomplete-Gamma representation would introduce principal-branch structure or a removable hole at the origin. For example,

```text
integrate[exp[x^6],x]
-> x hypergeometric1F1[1/6, 7/6, x^6]
```

The same family handles `exp[c x^n]` for positive integer `n`.

## 14.7 Gauss hypergeometric 2F1

The Gauss hypergeometric function is written as

```text
hypergeometric2F1[a,b,c,z]
```

In general `c = 0,-1,-2,...` is a parameter pole, and mmCal uses the principal branch in `z`. Exact evaluation currently handles terminating series generated by non-positive-integer numerator parameters and other safe degenerations such as `a=0` or `b=0`. The certified real `N` backend currently accepts exact Rational parameters with real `|z|<1`.

```text
hypergeometric2F1[-2,1,3,1/2] -> 17/24
hypergeometric2F1[0,2,3,x] -> 1
N[hypergeometric2F1[1/2,1/2,3/2,1/4],20]
-> 1.04719755119659774615
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

## 14.8 Incomplete elliptic integrals F / E / Pi

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
-> 0.50684775626543110920
N[ellipticE[1/2,1/3],20]
-> 0.49331536201475850521
N[ellipticPi[1/5,1/2,1/3],20]
-> 0.51520338216141386085
```

The current certified real backend accepts exact Rational amplitudes with `|m|<1` for `F/E`, and additionally `|n|<1` for `Pi`. Amplitude derivatives are

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


## 14.9 Ei / Si / Ci / li / Polylogarithm

The principal special functions commonly required by symbolic integration are exposed as

```text
Ei[x]
Si[x]
Ci[x]
li[x]
polylog[s,z]
```

`Ei`, `Ci`, `li`, and `polylog` generally have branch structure and are registered as principal-branch functions. `Si` is entire and odd. The current certified real backends deliberately cover only regions where a rigorous tail bound is available; unsupported regions are not filled with heuristic numeric values.

```text
Si[0] -> 0
Si[-1] -> -Si[1]
polylog[0,z] -> z/(1-z)
polylog[1,z] -> -log[1-z]
polylog[2,1] -> Pi^2/6
polylog[2,-1] -> -Pi^2/12

N[Ei[1],20] -> 1.89511781635593675547
N[Si[1],20] -> 0.94608307036718301494
N[Ci[1],20] -> 0.33740392290096813466
N[li[2],20] -> 1.04516378011749278484
N[polylog[2,1/2],20] -> 0.58224052646501250590
```

When the argument and order parameters are independent of the differentiation variable, the derivative knowledge includes

```text
D[Ei[x],x] -> exp[x]/x
D[Si[x],x] -> sin[x]/x
D[Ci[x],x] -> cos[x]/x
D[li[x],x] -> 1/log[x]
D[polylog[s,x],x] -> polylog[s-1,x]/x
```

The exact degeneration `polylog[1,x] -> -log[1-x]` gives

```text
D[polylog[2,x],x] -> -log[1-x]/x
```

and the same shared knowledge closes

```text
integrate[exp[x]/x,x] -> Ei[x]
integrate[sin[x]/x,x] -> Si[x]
integrate[cos[x]/x,x] -> Ci[x]
integrate[1/log[x],x] -> li[x]
integrate[log[1-x]/x,x] -> -polylog[2,x]
```

No general `Solve` rule invents a single principal inverse for `Ei/Si/Ci/li/polylog`: global injectivity and branch structure are not generally available. Only exact degenerations such as `polylog[0,z]` and `polylog[1,z]` are passed to the existing algebraic/logarithmic Solver.

---

# 15. Aggregate functions

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

Symbolic finite sums of the form `sum[f,{k,a,b}]` are not yet implemented.

---

# 16. Descriptive statistics

Statistical functions generally accept **exact real data** and preserve Rational results when the quantity closes rationally.
Many functions accept either a single rank-1 Array or a scalar argument list.

## 16.1 Order statistics

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

## 16.2 Variance and standard deviation

```text
var      // population variance
vars     // sample unbiased variance
stddev   // sqrt[var]
stddevs  // sqrt[vars]
```

```text
var[1,2,3] -> 2/3
vars[1,2,3] -> 1
stddev[1,2,3] -> sqrt[6] / 3
stddevs[1,2,3] -> 1
```

## 16.3 Other statistics

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

# 17. Array / Vector / Matrix

## 17.1 Array utilities

```text
identity[n]
zeros[rows,cols]
mget[A,row,col]
rows[A]
cols[A]
diag[A]
trace[A]
```

```text
identity[2]
-> {{1,0},{0,1}}

trace[{{1,2},{3,4}}]
-> 5
```

## 17.2 Exact-first linear algebra

```text
transpose[A]
madd[A,B,...]
matmul[A,B]
det[A]
inverse[A]
rref[A]
rank[A]
```

```text
matmul[{{1,2},{3,4}},{{5,6},{7,8}}]
-> {{19,22},{43,50}}

det[{{1,2},{3,4}}]
-> -2

inverse[{{1,2},{3,4}}]
-> {{-2,1},{3/2,-1/2}}
```

If nonzero symbolic pivots cannot be proven, symbolic rank/pivot routines do not choose pivots arbitrarily.

## 17.3 Vector

```text
vadd vsub vscalar
vdot vcross
vnorm vnormalize
vproject vangle
vmanhattan veuclidean
vreflect vreflect_axis
vsum
```

```text
vdot[{1,2},{3,4}] -> 11
vnormalize[{3,4}] -> {3/5,4/5}
vangle[{1,0},{0,1}] -> Pi/2
```

---

# 18. Signal processing

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

For exact inputs, power-of-two FFTs use radix-2 Cooley–Tukey and non-power-of-two lengths fall back to exact DFT. Ordinary `fft[...]` remains exact-first and never silently converts to machine `double`.

`N[fft[v],p]` does not first expand the full exact Fourier expression. `N` propagates the requested precision into the FFT call, which performs butterflies directly on certified `ComplexInterval`/BigFloat endpoints and returns decimal components only after their requested rounding is proven unique. `fft[v]` with approximate operands dispatches to the same backend.

The approximate path uses radix-2 for power-of-two sizes and Bluestein convolution for sufficiently large non-power-of-two sizes. Small non-power-of-two transforms keep direct DFT because its constant factor wins there; the current measured crossover policy uses direct evaluation below 96 points and should be remeasured with `mmCal.Benchmarks` on MSVC.

---
# 19. Symbolic and numerical differentiation

## 19.1 `D`

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

Integration reuses the same derivative knowledge.

```text
D[integrate[f[x],x],x]
-> f[x]

D[integrate[t^2,{t,0,x}],x]
-> x^2
```

For general forms with variable endpoints or where the differentiation variable occurs inside the integral, the Leibniz rule is constructed formally. A quotient whose denominator is independent of the differentiation variable uses `f'/c` directly rather than expanding into the general quotient rule.

## 19.2 `diff`

```text
diff[expr,x,at]
diff[expr,x,at,digits]
```

This is not a separate finite-difference formula. mmCal first constructs an exact derivative Expr using `D`, then passes that expression to the CertifiedEvaluator at the specified point.

```text
diff[x^2,x,3]
-> 6.0000000000000000
```

---

# 20. Symbolic integration and certified numerical integration

## 20.1 `integrate` — exact / symbolic integration

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
- Rational functions with Rational coefficients. In addition to linear and quadratic denominators, exact partial fractions are used when Rational roots reduce the denominator to linear factors with at most an irreducible quadratic remainder. Repeated linear factors are supported
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
- `hypergeometric1F1`; `exp[c x^n]` with positive integer `n>=2` reduces to an entire 1F1 primitive at the origin
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
```

If an internal pole cannot be excluded, the integral remains unevaluated.

```text
integrate[1/(x-2),{x,1,Infinity}]
-> WARN + unevaluated integrate[...]
```

## 20.2 `limit` — exact / symbolic limits

```text
limit[expr,x,a]
limit[expr,x,a,-1]   // left-hand limit
limit[expr,x,a,1]    // right-hand limit
```

The fourth argument specifies direction: `-1` means left and `1` means right. If omitted, a two-sided limit is requested. Direction is not merely a display option; it is supplied to `KnowledgeContext` as a temporary assumption `x<a` / `x>a`.

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
```

For indeterminate `0/0` forms, repeated l'Hôpital evaluation using the existing `D` implementation is available with safety limits. Local zero/pole orders of Rational functions and degree comparisons at infinity are handled exactly. Unresolved two-sided limits are not collapsed into principal values or other guessed results.

```text
limit[1/x,x,0]
-> WARN + limit[1 / x, x, 0]
```

## 20.3 `nintegrate` — certified numerical integration

```text
nintegrate[expr,{x,a,b}]
nintegrate[expr,{x,a,b},digits]
```

```text
nintegrate[x^2,{x,0,1},12]
-> 0.333333333333
```

Internally, the interval is normalized with `x=a+(b-a)t`. Newton–Cotes quadrature on exact Rational grids is combined with certified bounds on higher derivatives to enclose the quadrature error.

This is not a `double`-based Simpson calculation that returns a value because it merely appears close. A result is returned only when the final rounding is uniquely determined.

Before constructing higher derivatives, the original integrand is preflighted over the full interval so obvious singularities can be rejected early. No accidental cancellation across a singularity is accepted.

---

# 21. Expression transformation

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
```

---

# 22. Assumptions and domains

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

---

# 23. Solver

```text
solve[equation,x]
solve[equation,x,domainOrConstraint]
solve[{equations...},{variables...}]
```

The default ambient domain for equation systems is Complex.
Ordered inequalities are handled over Real or a real subdomain.

```text
solve[x^2 == 1,x]
-> {x == 1, x == -1}

solve[x^2 + 1 == 0,x,Real]
-> {}

solve[x^2 + 1 == 0,x,Complex]
-> {x == I, x == -I}

solve[x^2 < 4,x]
-> {x in Real if x > -2 && x < 2}
```

`SolutionSet` distinguishes Empty / Finite / Universal / Conditional / Unresolved.
Unsupported expressions are not misreported as having no solutions.

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
```

Inverse, range, and period metadata are also registered for `sin/cos/tan`, but these functions are not globally injective on the real axis. Because integer-parameter families of solutions are not yet implemented, mmCal does not collapse the full solution set to a single principal inverse.

```text
solve[sin[x]==0,x,Real]
-> WARN + UnresolvedSolutionSet[x]
```

The Complex domain likewise does not fabricate complete solution sets from principal inverses alone.

---

# 24. `N` — numerical approximation

```text
N[expr]
N[expr,digits]
```

The default is 16 fractional digits.
`N` applies recursively to Arrays. For explicit angle-unit values such as those returned by `arg`, only the numeric component is approximated and the unit is retained.

Since v1.5.2, `N` is also the entry point for precision-aware evaluation. It resolves the requested precision before evaluating its first argument and keeps that precision context active while the child expression is evaluated. Ordinary builtins still follow exact-first evaluation; only explicitly supported builtins such as FFT consume the context and evaluate directly in a certified approximate domain.

```text
N[Pi,20]
-> 3.14159265358979323846
N[Phi,20]
-> 1.61803398874989484820
N[fft[{1,2,3,4}],20]
N[arg[-1],20]
-> 3.14159265358979323846 Rad
```

When an exact Rational has a terminating decimal representation, unnecessary trailing zeros are not displayed. The current implementation nevertheless retains `requestedFractionalDigits` separately as metadata.

Approximate values are not yet a general Machine/ApproximateReal arithmetic domain. The displayed value, requested digit count, and certified enclosure are stored separately.

---

# 25. precision / accuracy / rationalize

## 25.1 `accuracy[x]`

For a `DecimalApproximation`, returns an **integer lower bound on the guaranteed number of absolute decimal digits** relative to the true value.

```text
accuracy[N[1/3,20]]
-> 20
```

Given displayed value `d` and certified source enclosure `[l,u]`, the implementation computes

```text
max(|d-l|, |d-u|)
```

and also applies `0.5*10^-n` as a semantic error floor so that an `N[...,n]` approximation retains the meaning of an approximate value. Therefore, even if a terminating decimal happens to coincide exactly with the true value, it does not acquire infinite accuracy beyond the requested digits.

Exact numbers and exact symbolic expressions return `Infinity`.

```text
accuracy[1/3] -> Infinity
accuracy[Pi]  -> Infinity
```

## 25.2 `precision[x]`

Uses the same absolute error bound divided by a lower bound on the absolute value of the true value derived from the certified enclosure, returning an **integer lower bound on the guaranteed number of relative decimal digits**.

```text
precision[N[1/3,20]]
-> 19
```

This function does not mechanically return the requested 20 digits. Around `1/3`, an absolute error floor of `0.5*10^-20` corresponds to a relative error of approximately `1.5*10^-20`, so only 19 integer decimal digits can be guaranteed. If the enclosure contains zero, no positive lower bound on relative accuracy is available, so the result is 0. Exact expressions return `Infinity`.

## 25.3 `rationalize[x]`

Finds the **exact Rational with the smallest denominator** contained in the certified enclosure of an approximate value. The search uses continued-fraction-style interval recursion over exact Rational values and never converts the interval to `double`.

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

---

# 26. Random numbers

Random-number functions are stateful built-ins.
The RNG state is independent for each `KernelSession`.

## 26.1 Seed

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

## 26.2 Uniform real

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

## 26.3 Integer

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

## 26.4 Choice

```text
choice[2,3,5,7]
choice[{2,3,5,7}]
```

Selects one element from either a rank-1 Array or a variadic list.

## 26.5 Normal

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

# 27. Conditional evaluation

```text
if[condition,trueExpr,falseExpr]
```

As a special form, `if` evaluates the condition first and evaluates only the selected branch.
Therefore a DomainError or random-number consumption in the unselected branch does not occur.

---

# 28. Major aliases

| Alias | Canonical |
|---|---|
| `pow` | `Power` |
| `fact` | `Factorial` |
| `fract` | `frac` |
| `ln` | `log` |
| `real` | `re` |
| `imag` | `im` |
| `mag` | `abs` |
| `unit`, `csgn` | `sign` |
| `rect` | `polar` |
| `ave` | `mean` |
| `mmul` | `matmul` |
| `mtranspose` | `transpose` |
| `mdet` | `det` |
| `minverse` | `inverse` |
| `mrank` | `rank` |
| `mtrace` | `trace` |
| `mrows` | `rows` |
| `mcols` | `cols` |
| `mdiag` | `diag` |
| `vlength` | `vnorm` |
| `vdistance` | `veuclidean` |
| `vunit` | `vnormalize` |

Aliases are not separate implementations; they resolve to the same `BuiltinId`. Mathematical metadata and Solver rules are therefore not duplicated.

In mmCal 1.5.0, capitalized aliases added only for Mathematica compatibility (`Sin`, `ArcTan`, `Integrate`, `Solve`, etc.) were removed. Lowercase canonical names are the default for mathematical functions. `D`, `N`, `In`, `Out`, `Exit`, `Clear`, `Defs`, and `UnDef` remain as intentional proper names for symbolic and Kernel operations. If compatibility syntax becomes necessary, it should be implemented as a separate import/compatibility layer rather than by adding aliases to the default namespace.

---

# 29. Current source-callable function list

The current development tree contains **217 built-in definitions / 199 source-callable names**. Internal heads are not included in the source-callable count.

```text
Clear, D, Defs, DtoG, DtoR, Exit, GtoD, GtoR, In, N,
Out, RtoD, RtoG, UnDef, abs, accuracy, acos, acosh, angleMode, arg,
asin, asinh, atan, atan2, atanh, ave, beta, betaln, binom, cbrt,
ceil, choice, cis, collect, cols, comb, conj, convolve, corr, corrspearman,
cos, cosc, cosh, cot, coth, cov, csc, csch, csgn, cv,
det, dft, diag, diff, element, erf, erfc, exp, expand, expc,
Ei, Si, Ci, li, polylog, fresnelc, fresnels, hypergeometric1F1, hypergeometric2F1, ellipticF, ellipticE, ellipticPi,
expm1, fact, factor, fallingfact, fft, fib, floor, frac, fract, fullSimplify,
gamma, gcd, geomean, harmmean, hypot, identity, if, ifft, im, imag,
integrate, inverse, iqr, kurtp, kurts, lcm, lgamma, limit, ln, log,
log10, log1p, log2, mad, madR, madd, mag, matmul, max, mcols,
mdet, mdiag, mean, median, mget, min, minverse, mmul, mod, mode,
mrank, mrows, mtrace, mtranspose, nextpow2, nintegrate, percentile, percentrank, perm, polar,
pow, precision, prod, quantile, quotient, rand, randSeed, randint, randn, rank,
rationalize, re, real, rect, rem, risingfact, rms, round, rows, rref,
sec, sech, sign, simplify, sin, sinc, sinh, sinhc, skew, solve,
sqrt, stddev, stddevs, stderr, sum, tan, tanc, tanh, tanhc, trace,
transpose, trimmean, trunc, unit, vadd, vangle, var, vars, vcross, vdistance,
vdot, veuclidean, vlength, vmanhattan, vnorm, vnormalize, vproject, vreflect, vreflect_axis, vscalar,
vsub, vsum, vunit, winsor, winsorR, zeros, zscore
```

---

# 30. Error / Warning policy

Main `CalcError` categories:

- Syntax
- Domain
- Type
- Overflow
- Name
- Evaluation
- Internal

When evaluation itself succeeds but an algorithmic built-in cannot complete its work, a Warning is returned separately from the result Expr. This applies to `D`, `solve`, `N`, `rref`, `rank`, and `integrate`, as well as cases where `precision/accuracy/rationalize` retain unsupported input unevaluated.

Info diagnostics are used for normal state-change notifications. Currently this includes variable and function redefinition notices.

```text
D[abs[x],x]
WARN: D could not fully evaluate the derivative; unevaluated D[...] remains
-> D[abs[x], x]
```

An expression such as `sin[x]` that is correctly retained symbolically does not produce a Warning.

Values that are mathematically undefined generally produce an Error immediately rather than flowing through evaluation as NaN/Inf.

Examples:

```text
1/0
-> DomainError

log[0]
-> DomainError

gamma[-2]
-> DomainError

randint[5,1]
-> DomainError
```

The Parser/Evaluator retains source spans and documents, allowing Errors that propagate through user-function definitions to include call traces.

---

# 31. Performance policy

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

# 32. Major currently unimplemented / deferred features

See `docs/roadmap.md` for future candidates and the reasons they are deferred.

Representative items:

- `digamma`, `trigamma`, `zeta`, `ibeta`
- `isprime`, `nextprime`, `prevprime`, `factorint`, `totient`
- Advanced LU/QR/SVD/eigen/condition number/least squares
- `hilbert` (legacy naming/specification still to be confirmed)
- `fma`, `clamp`, `proj`
- Engineering functions, financial functions, and unit conversion
- Legacy colon commands such as `:defs`, `:help`, `:unset`, and `:undef` (`Defs[]/UnDef[]` function forms are implemented). `:angle` has been replaced by `angleMode[]`; `:fix` / `:status` remain as presentation commands
- `for`, `plot`
- General Machine/double evaluation mode

Further candidates:

AlgebraicNumber / Root infrastructure, `rootApproximant`, integer-parameter solution families for periodic functions, and a Machine evaluator / `for` / `plot`.

---

# 33. Current CLI

The CLI separates mathematical Kernel state from frontend presentation state. The following startup options are available:

```text
mmCal --fix 16 --angle deg
mmCal --angle rad
mmCal --angle grad --fix 8
```

- `--fix n`: Set the startup limit on decimal display digits. Internal values are unchanged and unnecessary trailing zeros are omitted
- `--angle deg|rad|grad`: Set the default angle unit at startup
- `--help`, `-h`: Show usage

```text
In[1]> 1/3
Out[1]> 1/3
```

- Parse/evaluate one line at a time
- Prompts are fixed as `In[n]>` / `Out[n]>` with no extra spaces
- Exit is unified under `Exit[]`; bare `exit` / `quit` receive no special treatment
- `Clear[]`: Remove user definitions and all history, resetting the next input number to 1
- `Defs[]`, `UnDef[...]`: Inspect and remove user definitions
- History references `@`, `%`, `%%`, ... together with signed-index reevaluating `In[n]` and snapshot `Out[n]`

## 33.1 `:fix` — presentation-only decimal display

```text
:fix 16
Display: Fixed(16)

In[1]> 1/3
Out[1]> 0.3333333333333333

:fix off
Display: Exact

In[2]> Out[1]
Out[2]> 1/3
```

`:fix n` changes **display only**, rounding to at most `n` digits after the decimal point. Unnecessary trailing zeros are omitted, so for example `31/10` is displayed as `3.1` under `:fix 5`. Stored Expr values, `Out[n]`, and the semantics of `precision/accuracy` are unchanged. This is not a conversion from exact values to Machine/double. The current range for `n` is 0..1000. Entering `:fix` with no argument displays the current mode.

An entire expression that can be certified numerically is approximated only for display. Symbolic expressions containing free variables retain exact notation.

## 33.2 `:status`

```text
:status
Angle: Rad
Display: Exact
Evaluation: Exact-first
Definitions: 0
History: 0
```

`:status` is also a CLI command and is not stored in history. The design boundary is maintained: mathematical state changes use Kernel functions such as `angleMode[...]`, while presentation state changes use CLI commands such as `:fix`.

## 33.3 Console title

As auxiliary information, the title is updated to forms such as:

```text
mmCal 1.5.0 - Rad - Exact
mmCal 1.5.0 - Deg - Fixed(16)
```

- Windows: `SetConsoleTitleA`
- Linux/macOS: ANSI OSC title sequence only when attached to a TTY
- Other platforms: no-op

Failure to change the title is never treated as a calculation Error. `:status` is authoritative for state inspection; terminal software overriding the title has no effect on semantics.

## 33.4 Canonical formatter

Ordinary `Out[n]` uses compact, reparsable mathematical notation rather than an AST dump.

```text
x^2+sin[x]
A-B+C
2(x+sqrt[x])sqrt[x+sqrt[x]]/3
```

- Do not print `+ -` / `+-`; negative terms are rendered with `-`
- Addition/subtraction forms such as `A-(B-C)` may be flattened for display to `A-B+C`
- Do not add unnecessary spaces around `+`, `-`, `*`, `/`, `^`, or comparison operators
- Implicit multiplication is concatenated only when lexically safe (`2x`, `2sqrt[x]`). Forms that would collide with exponent notation, such as `2exp[x]` or `2E`, are rendered explicitly as `2*exp[x]`, `2*E`
- Preserve necessary spacing when adjacent identifiers would merge into another token (`I Pi`, `x y`)
- When juxtaposition would be ambiguous, such as adjacent numeric tokens, use explicit `*` rather than whitespace
- Preserve precedence and associativity, and regression-test that format → parse → format does not change meaning

Debug/full-form display of internal structure is intended to remain separate from the normal formatter as a future feature.

---

# 34. Major implementation layers

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
