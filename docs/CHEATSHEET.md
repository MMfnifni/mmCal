# mmCal Cheatsheet

Based on the current development tree. For details, see `docs/reference.md`. For quick in-session lookup, use `:help <function-name>`.

## 1. Input syntax

| Purpose | Syntax | Example |
|---|---|---|
| Function call | `name[...]` | `sin[Pi/6]` |
| Grouping | `(...)` | `(x+1)^2` |
| Array | `{...}` | `{1,2,3}` |
| Matrix | `{{...},...}` | `{{1,2},{3,4}}` |
| Assignment | `:=` | `x:=3` |
| User-defined function | `f[x]:=...` | `f[x]:=x^2+1` |
| Numerical approximation | `N[expr,p]` | `N[Pi,30]` |

**Function calls use `[]` only. Parentheses `()` are reserved for grouping.**

```text
sin[Pi/6]
(x+1)^2
```

`sin(Pi/6)` is invalid legacy-style syntax for a known function and raises `SyntaxError`. With an ordinary identifier, `x(x+1)` is interpreted as implicit multiplication, `x*(x+1)`.

Implicit multiplication is supported.

```text
2Pi
2sqrt[2]
2(x+1)
(x+1)(x-1)
```

Operators:

```text
+  -  *  /  ^  !
```

```text
2^3^2   -> 512
-2^2    -> -4
```

Comparisons:

```text
==  !=  <  <=  >  >=
```

## 2. Numbers and constants

```text
123
0.125          -> 1/8
1/3
3+4I
Pi
E
Phi
Infinity
```

Finite decimal literals are exact `Rational` values by default. The imaginary unit is uppercase `I`.

Base-prefixed literals:

```text
0b1010
0o17
0xFF
2#1010
16#FF
```

Domain symbols:

```text
Integer  Rational  Real  Complex
```

## 3. Angles

**The default angle unit is Radian (`Rad`).**

```text
sin[Pi/6]       -> 1/2
sin[30 Deg]     -> 1/2
sin[Pi/6 Rad]   -> 1/2
sin[100 Grad]   -> 1
```

Session setting:

```text
angleMode[]        -> Rad
angleMode[Deg]     -> Deg
angleMode[Rad]     -> Rad
angleMode[Grad]    -> Grad
```

An explicit angle unit overrides the session setting. Units can also be attached to general expressions.

```text
sin[x Deg]
sin[(2x+1) Grad]
```

Conversions:

```text
DtoR[180] -> Pi    DtoG[90]  -> 100
RtoD[Pi]  -> 180   RtoG[Pi]  -> 200
GtoD[200] -> 180   GtoR[200] -> Pi
```

## 4. Exact values and numerical approximation

```text
1/3              -> 1/3
sqrt[2]          -> sqrt[2]
0.1+0.2==0.3     -> True
```

```text
N[1/3,20]
N[Pi,30]
N[sqrt[2],50]
```

In `N[...,p]`, `p` is the number of **significant decimal digits**.

```text
accuracy[N[1/3,20]]
precision[N[1/3,20]]
rationalize[N[1/3,20]]   -> 1/3
explain[N[Pi,20]]
```

`:fix n` changes only the **maximum number of digits after the decimal point used for display**; it does not change the computed value.

## 5. Variables, user-defined functions, and history

```text
x:=3
f[t]:=t^2+1
f[4]              -> 17
```

```text
Defs[]
UnDef[x]
UnDef[x,y,f]
Clear[]
Exit[]
```

History:

| Input | Meaning |
|---|---|
| `%`, `%%`, ... | Previous, second previous, ... successful output |
| `@`, `@@`, ... | Re-evaluate the previous, second previous, ... input in the current environment |
| `Out[n]` | Stored output snapshot |
| `In [n]` | Re-evaluate a previous input in the current environment |

```text
N[@,30]
Out[-1]
In [-1]
```

## 6. Conditionals and assumptions

Evaluation control:

```text
if[x>0,1,-1]
```

Mathematical piecewise expression:

```text
cases[1/x if x!=0; 0 if x==0]
```

Assumptions:

```text
element[x,Real]
element[n,Integer]
simplify[sqrt[x^2],element[x,Real]]   -> abs[x]
simplify[sqrt[x^2],x>=0]             -> x
```

## 7. Basic mathematics

| Category | Main forms |
|---|---|
| Roots / absolute value | `sqrt[x]` `cbrt[x]` `abs[x]` `sign[x]` |
| Complex numbers | `re[z]` `im[z]` `conj[z]` `arg[z]` |
| Exponential / logarithm | `exp[x]` `log[x]` `log[b,x]` `log2[x]` `log10[x]` |
| Stable forms | `expm1[x]` `log1p[x]` `sinc[x]` `cosc[x]` `tanc[x]` |
| Trigonometric | `sin[x]` `cos[x]` `tan[x]` `cot[x]` `sec[x]` `csc[x]` |
| Inverse trigonometric | `asin[x]` `acos[x]` `atan[x]` `atan2[y,x]` |
| Hyperbolic | `sinh[x]` `cosh[x]` `tanh[x]` `csch[x]` `sech[x]` `coth[x]` |
| Inverse hyperbolic | `asinh[x]` `acosh[x]` `atanh[x]` |

```text
sqrt[8]       -> 2sqrt[2]
abs[3+4I]     -> 5
re[3+4I]      -> 3
im[3+4I]      -> 4
conj[3+4I]    -> 3-4I
arg[-1]       -> Pi Rad
log[10,1000]  -> 3
```

## 8. Rounding, integers, and number theory

```text
floor[-3/2]    -> -2
ceil[-3/2]     -> -1
trunc[-3/2]    -> -1
round[5/2]     -> 2
frac[-3/2]     -> 1/2
```

```text
gcd[84,126,210] -> 42
lcm[6,8,9]      -> 72
mod[-5,3]       -> 1
rem[-5,3]       -> -2
quotient[-5,3]  -> -1
```

```text
10!
perm[10,3]      -> 720
comb[10,3]      -> 120
fib[100]
isprime[97]     -> True
factorint[360]
totient[9]      -> 6
```

## 9. Expression transformation

```text
simplify[expr]
simplify[expr,assumptions]
fullSimplify[expr]
fullSimplify[expr,assumptions]
expand[expr]
factor[expr]
collect[expr,x]
```

```text
expand[(x+1)^3]
factor[x^2-1]
fullSimplify[(x^2-1)/(x-1),x!=1] -> x+1
```

## 10. Series

```text
series[expr,{x,a,n}]
series[expr,{x,a,n},assumptions]
```

```text
series[exp[x],{x,0,5}]
series[1/(1+x),{x,0,5}]
series[1/(x+1),{x,Infinity,4}]
```

Convert back to an ordinary expression:

```text
normal[series[(1+x)^3,{x,0,5}]]
-> x^3+3x^2+3x+1
```

`normal` acts only on the top level. To recursively convert supported `Series` objects inside an expression:

```text
toNormal[expr]
```

## 11. Differentiation, integration, and limits

```text
D[x^3,x]                    -> 3x^2
D[x^5,{x,3}]
D[x^2*y^3,x,y]
D[exp[x^2],x]               -> 2x exp[x^2]
```

Numerical differentiation at a point:

```text
diff[x^2,x,3]               -> 6.0
diff[sin[x],x,1,20]
```

Integration:

```text
integrate[x^2,x]                    -> x^3/3
integrate[sin[x],x]                 -> -cos[x]
integrate[x^2,{x,0,1}]              -> 1/3
integrate[exp[-x],{x,0,Infinity}]   -> 1
integrate[expr,x,assumptions]
```

Certified numerical integration:

```text
nintegrate[x^2,{x,0,1},12]
```

Limits:

```text
limit[sin[x]/x,x,0]   -> 1
limit[1/x,x,0,1]      -> Infinity
limit[1/x,x,0,-1]     -> -Infinity
```

Direction `1` means from the right; `-1` means from the left.

## 12. Solve

```text
solve[x^2==1,x]
-> {x == 1, x == -1}

solve[x^2+1==0,x,Real]
-> {}

solve[x^2+1==0,x,Complex]
-> {x == I, x == -I}

solve[x^2<4,x]
solve[{2x+3y==5,x-2y==9},{x,y}]
```

For equation systems, the default ambient domain is `Complex`. Ordered inequalities use the `Real` domain family.

Finite-solution branch selection is **0-based**.

```text
s:=solve[x^2==1,x]
at[s,0]     -> {x == 1}
at[s,1]     -> {x == -1}
at[s,1,x]   -> -1
```

Exact roots:

```text
root[{-2,0,1},2]
N[root[{-2,0,1},2],30]
```

For `root[{a0,a1,...,an},k]`, polynomial coefficients are listed in ascending powers.

## 13. Arrays and sequence generation

```text
{1,2,3}
{{1,2},{3,4}}
```

```text
dimensions[{{1,2,3},{4,5,6}}] -> {2,3}
arrayRank[{{1,2},{3,4}}]       -> 2
length[{10,20,30}]             -> 3
at[{{1,2},{3,4}},1]            -> {3,4}
at[{{1,2},{3,4}},1,0]          -> 3
```

**`at` uses 0-based indexing.**

```text
reshape[{1,2,3,4},{2,2}]
zeros[2,3]
identity[3]
transpose[{{1,2},{3,4}}]
```

Sequence generation and mapping:

```text
range[5]             -> {1,2,3,4,5}
range[0,1,1/3]       -> {0,1/3,2/3,1}
table[i^2,{i,5}]     -> {1,4,9,16,25}
map[sin,{0,Pi/2,Pi}] -> {0,1,0}
```

`sin[A]`, `exp[A]`, and similar forms are not automatically element-wise. Use `map` explicitly.

Aggregations:

```text
sum[{1,2,3}]    -> 6
prod[{1,2,3}]   -> 6
min[{1,2,3}]
max[{1,2,3}]
mean[{1,2,4}]   -> 7/3
```

## 14. Matrices

**`A*B` is not matrix multiplication. Use `dot[A,B]` for matrix multiplication and contraction.**

```text
dot[{{1,2},{3,4}},{{5,6},{7,8}}]
-> {{19,22},{43,50}}
```

Main functions:

```text
transpose[A]             conjugateTranspose[A]
det[A]                   inverse[A]
rref[A]                  matrixRank[A]
nullSpace[A]             solveLinear[A,b]
luDecomposition[A]       qrDecomposition[A]
svd[A]                   conditionNumber[A]
pseudoInverse[A]         leastSquares[A,b]
eigenvalues[A]           eigenvectors[A]
eigensystem[A]           trace[A]
```

```text
det[{{1,2},{3,4}}]                -> -2
solveLinear[{{2,1},{1,-1}},{5,1}] -> {2,1}
```

Use `at` to extract decomposition results as well.

```text
lu:=luDecomposition[A]
at[lu,0]
at[lu,1]
at[lu,2]
```

## 15. Vectors

Rank-1 arrays are used as vectors.

```text
norm[{3,4}]       -> 5
normalize[{3,4}]  -> {3/5,4/5}
cross[{1,0,0},{0,1,0}]
distance[a,b]
outer[a,b]
projection[a,b]
rejection[a,b]
```

For complex vectors, `dot` and `inner` are distinct.

```text
inner[{1,I},{1,I}] -> 2
dot[{1,I},{1,I}]   -> 0
```

- `dot[a,b]`: bilinear; no conjugation.
- `inner[a,b]`: Hermitian inner product; conjugates the first argument.

Orthogonalization:

```text
orthogonalQ[vectors]
orthonormalQ[vectors]
linearIndependentQ[vectors]
gramSchmidt[vectors]
```

Cartesian vector calculus:

```text
grad[f,{x,y,z}]
divergence[{P,Q,R},{x,y,z}]
curl[{P,Q,R},{x,y,z}]
laplacian[f,{x,y,z}]
jacobian[{f1,f2},{x,y}]
hessian[f,{x,y}]
```

## 16. Statistics, signals, and random numbers

Statistics:

```text
mean[{1,2,4}]                  -> 7/3
median[1,2,3,4]                -> 5/2
quantile[1/4,1,2,3,4,5,6,7]   -> 5/2
var[1,2,3]                     -> 2/3
stddev[1,2,3]                  -> sqrt[6]/3
corr[{1,2,3},{2,4,6}]          -> 1
```

Signals:

```text
dft[v]
fft[v]
ifft[v]
convolve[a,b]
ifft[fft[{1,2,3,4}]] -> {1,2,3,4}
```

Fourier phases are always in Radians, independent of the session angle mode.

Random numbers:

```text
randSeed[42]
rand[]
randint[1,6]
choice[{a,b,c}]
randn[]
```

## 17. Special functions

Representative input forms. For branches, domains, and detailed behavior, use `:help <function-name>`.

```text
gamma[x]                  lgamma[x]
digamma[x]                trigamma[x]
erf[x]                    erfc[x]
beta[a,b]                 ibeta[a,b,x]       betaln[a,b]
zeta[s]
lambertw[x]               lambertw[k,x]
Ei[x]                     Si[x]               Ci[x]
li[x]                     polylog[s,z]
fresnelc[x]               fresnels[x]
hypergeometric1F1[a,b,z]
hypergeometric2F1[a,b,c,z]
ellipticF[phi,m]          ellipticE[phi,m]    ellipticPi[n,phi,m]
```

The amplitude `phi` of `ellipticF/E/Pi` is always interpreted in Radians, regardless of the session angle mode.

## 18. CLI quick reference

```text
mmCal --angle rad
mmCal --angle deg
mmCal --angle grad
mmCal --fix 16
mmCal --layout auto
mmCal --layout single
mmCal --layout multi
mmCal --eval "factor[x^2-1]"
mmCal --batch < expressions.txt
```

REPL:

```text
:help
:help sin
:help functions
:help constants
:status
:fix 16
:fix off
:layout auto
:layout single
:layout multi
:quit
:exit
```

## 19. Common mistakes

| Incorrect | Correct |
|---|---|
| `sin(x)` | `sin[x]` |
| `[x+1]` for grouping | `(x+1)` |
| Treating `sin[30]` as 30 degrees | `sin[30 Deg]` |
| `:angle deg` | `angleMode[Deg]` |
| Treating `A*B` as matrix multiplication | `dot[A,B]` |
| Using `at[A,1]` for the first element | `at[A,0]` |
| Assuming `sin[A]` is element-wise | `map[sin,A]` |

When in doubt:

```text
:help <name>
:help functions
:status
```
