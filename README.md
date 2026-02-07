# mmCalculator – Mechanical & Manufacturing Calculator

© 2021–2026 mmKreuzef (aka Daiki.NIIMI) Project mme
Licensed under the BSD 3-Clause License

[English](README.md) | [日本語](README.ja.md)

### Overview

This project is a **mathematical expression evaluation engine** featuring **numerical computation, complex number arithmetic, and high-precision mathematical functions**.

While it can be used as a calculator, it is primarily designed for **CAD, Mechanical processing, and practical use at manufacturing and design sites**, where both numerical accuracy and usability are required.

The design is inspired by _Mathematica_ in spirit, but **variable assignment and symbolic computation are intentionally not implemented**.

---

### Numeric Types

- **Real numbers**: IEEE 754 `double` (displayed with up to 12 decimal digits)
- **Complex numbers**: `a + bI`

---

### Constants

All constants start with an capital letter.

| Name  | Description                    | Expansion value                         |
| ----- | ------------------------------ | --------------------------------------- |
| `Pi`  | π (pi)                         | `3.14159265358979323846264338327950288` |
| `Tue` | 2π (τ, tau)                    | `6.283185307179586476925286766559006`   |
| `E`   | Base of natural logarithm      | `2.7182818284590452353602874713526625`  |
| `Phi` | Golden ratio                   | `1.618033988749894848204586834365638`   |
| `NA`  | Avogadro constant              | ` 6.02214076e23`                        |
| `I`   | Imaginary unit                 | `0+I`                                   |
| `ESP` | Internal convergence tolerance | `1e-12`                                 |

---

### Operators

#### Arithmetic, Power, and Factorial

| Operator | Description                            | Associativity |
| -------- | -------------------------------------- | ------------- |
| `+ -`    | Addition / subtraction                 | Left          |
| `* /`    | Multiplication / division              | Left          |
| `^`      | Power                                  | **Right**     |
| `!`      | Factorial (non-negative integers only) | Postfix       |

```text
2^3^2 = 2^(3^2)
-2^2 = -(2^2)
```

---

### Unary Operators

- Unary `+` and `-` may be chained arbitrarily.

```text
---5 = -5
```

---

### Implicit Multiplication

The following patterns are **automatically interpreted as multiplication**:

```text
2Pi
2(Pi)
Pi(2)
(2+3)(4+5)
2 4
```

But, please be mindful of the following types of input.

| Imput   | Calculator's Thinking | Output                          | Expected result (or input error) |
| ------- | --------------------- | ------------------------------- | -------------------------------- |
| `1..1`  | `1. .<-???`           | `Error: multiple '.' in number` | `1. * .1 = 0.1`                  |
| `1.2.1` | `1.2 .<-???`          | `Error: multiple '.' in number` | `1.2 * .1 = 0.12`                |

**Function names and numeric values do not implicitly combine.**

```text
[NG]
sin30        (Error)
sin2(30)    (Error)

[OK]
sin(30)
cos(15+30)
```

Functions must always be enclosed in `()` or `[]`.

---

### Functions

Below is a complete list of all currently implemented functions. For each function, the function name, a brief description, and representative input/output examples are provided.

---

### Basic Math (Algebra, Exponentials, Logarithms, Rounding)

| Function Name    | Description                               | Example Inputs/Outputs           |
| ---------------- | ----------------------------------------- | -------------------------------- |
| `abs(x)`         | Absolute value (supports complex numbers) | `abs(-3) = 3`<br>`abs(3+4I) = 5` |
| `sign(x)`        | Sign function                             | `sign(-5) = -1`<br>`sign(0) = 0` |
| `sqrt(x)`        | Square root (negative → complex)          | `sqrt(4) = 2`<br>`sqrt(-4) = 2I` |
| `cbrt(x)`        | Cube root                                 | `cbrt(8) = 2`                    |
| `exp(x)`         | Exponential function                      | `exp(1) = 2.71828...`            |
| `expm1(x)`       | `exp(x) - 1`                              | `expm1(1)=1.718281828459`        |
| `log(x)`,`ln(x)` | Natural logarithm                         | `log(E) = 1`                     |
| `log(b, x)`      | Logarithm with arbitrary base             | `log(10, 1000) = 3`              |
| `log10(x)`       | Base-10 logarithm                         | `log10(100) = 2`                 |
| `log2(x)`        | Base-2 logarithm                          | `log2(8) = 3`                    |
| `log1p(x)`       | `log(1 + x)`                              | `log1p(E-1)=1`                   |
| `floor(x)`       | Floor                                     | `floor(3.7) = 3`                 |
| `ceil(x)`        | Ceiling                                   | `ceil(3.1) = 4`                  |
| `trunc(x)`       | Truncation toward zero                    | `trunc(-1.9)=-1`                 |
| `round(x)`       | Rounding (nearest integer)                | `round(2.6) = 3`                 |
| `round(x,n)`     | Rounding to _n_ digits                    | `round(1.2345,2)=1.23`           |
| `fract(x)`       | Fractional part                           | `fract(3.14)=0.14`               |
| `fact(x)`        | Factorial                                 | `fact(10)=3628800`               |

---

### Special Functions

| Function Name | Description            | Example Inputs/Outputs        |
| ------------- | ---------------------- | ----------------------------- |
| `gamma(x)`    | Gamma function         | `gamma(5)=24`                 |
| `lgamma(x)`   | Log-gamma function     | `lgamma(5)=3.178053830348...` |
| `erf(x)`      | Error function         | `erf(1)=0.842700792949...`    |
| `erfc(x)`     | Complementary error fn | `erfc(1)=0.157299207050...`   |

---

#### Angle Conversion

| Function Name | Description     | Example Inputs/Outputs |
| ------------- | --------------- | ---------------------- |
| `DtoR(x)`     | Degree → Radian | `DtoR(180)=Pi`         |
| `DtoG(x)`     | Degree → Grad   | `DtoG(90)=100`         |
| `RtoD(x)`     | Radian → Degree | `RtoD(Pi)=180`         |
| `RtoG(x)`     | Radian → Grad   | `RtoG(Pi)=200`         |
| `GtoD(x)`     | Grad → Degree   | `GtoD(200)=180`        |
| `GtoR(x)`     | Grad → Radian   | `GtoR(200)=Pi`         |

---

### Trigonometric Functions (Degree Mode / degree)

| Function Name | Description               | Example Inputs/Outputs |
| ------------- | ------------------------- | ---------------------- |
| `sin(x)`      | Sine                      | `sin(30)=0.5`          |
| `cos(x)`      | Cosine                    | `cos(60)=0.5`          |
| `tan(x)`      | Tangent                   | `tan(45)=1`            |
| `cot(x)`      | Cotangent                 | `cot(45)=1`            |
| `sec(x)`      | Secant                    | `sec(60)=2`            |
| `csc(x)`      | Cosecant                  | `csc(30)=2`            |
| `asin(x)`     | Inverse sine              | `asin(0.5)=30`         |
| `acos(x)`     | Inverse cosine            | `acos(0.5)=60`         |
| `atan(x)`     | Inverse tangent           | `atan(1)=45`           |
| `atan2(y,x)`  | Quadrant-aware arctangent | `atan2(1,1)=45`        |

---

### Hyperbolic Functions

| Function Name | Description             | Example Inputs/Outputs |
| ------------- | ----------------------- | ---------------------- |
| `sinh(x)`     | Hyperbolic sine         | `sinh(1)=1.175...`     |
| `cosh(x)`     | Hyperbolic cosine       | `cosh(0)=1`            |
| `tanh(x)`     | Hyperbolic tangent      | `tanh(0)=0`            |
| `asinh(x)`    | Inverse hyperbolic sine | `asinh(1)=0.881...`    |
| `acosh(x)`    | Inverse hyperbolic cos  | `acosh(1)=0`           |
| `atanh(x)`    | Inverse hyperbolic tan  | `atanh(0.5)=0.549...`  |

---

### Special Ratio Functions (c-series / degree)

| Function Name | Description    | Example Inputs/Outputs   |
| ------------- | -------------- | ------------------------ |
| `sinc(x)`     | `sin(x)/x`     | `sinc(90)=0.111...`      |
| `cosc(x)`     | `(1-cos(x))/x` | `cosc(60)=0.008333...`   |
| `tanc(x)`     | `tan(x)/x`     | `tanc(0)=1`              |
| `sinhc(x)`    | `sinh(x)/x`    | `sinhc(3)=3.339...`      |
| `tanhc(x)`    | `tanh(x)/x`    | `tanhc(30)=0.03333...`   |
| `expc(x)`     | `(exp(x)-1)/x` | `expc(1)=1.718281828...` |

- Since sin/cos/tan assume degree input, these also take degree input.
- When `x` is extremely small, loss-of-significance countermeasures are applied.

| Function Name | Description    | Example Inputs/Outputs |
| ------------- | -------------- | ---------------------- |
| `sincR(x)`    | `sin(x)/x`     | `sincR(1e-8)≈1`        |
| `coscR(x)`    | `(1-cos(x))/x` | `coscR(1e-8)≈5e-9`     |
| `tancR(x)`    | `tan(x)/x`     | `tancR(1e-8)≈1`        |

- These interpret `x` as radians.
- For the signal-processing-style “sinc”, this is generally the one you want.

---

### Number Theory / Combinatorics

| Function Name | Description                                            | Example Inputs/Outputs      |
| ------------- | ------------------------------------------------------ | --------------------------- |
| `gcd(a,b)`    | Greatest common divisor                                | `gcd(12,18)=6`              |
| `lcm(a,b)`    | Least common multiple                                  | `lcm(6,8)=24`               |
| `perm(n,r)`   | Permutations [(returns 0 if r<0 or r>n)]               | `perm(5,2)=20`              |
| `comb(n,r)`   | Combinations [(returns 0 if r<0 or r>n)]               | `comb(5,2)=10`              |
| `isprime(x)`  | Prime number determination<br>(0: non-prime, 1: prime) | `isprime(67280421310721)=1` |
| `fib(x)`      | Fibonacci sequence                                     | `fib(25)=75025`             |

---

### Statistics (Descriptive Statistics / Aggregation)

| Function Name          | Description                                           | Example Inputs/Outputs                                                                                                                                                             |
| ---------------------- | ----------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sum(...)`             | Sum                                                   | `sum(1,2,3)=6`                                                                                                                                                                     |
| `mean(...)`, `ave()`   | Mean                                                  | `mean(1,2,3)=2`                                                                                                                                                                    |
| `min(...)`             | Minimum                                               | `min(3,1,2)=1`                                                                                                                                                                     |
| `max(...)`             | Maximum                                               | `max(3,1,2)=3`                                                                                                                                                                     |
| `prod(...)`            | Product                                               | `prod(2,3,4)=24`                                                                                                                                                                   |
| `median(...)`          | Median                                                | `median(1,3,5)=3`                                                                                                                                                                  |
| `mode(...)`            | Mode                                                  | `mode(1,2,2,3)=2`                                                                                                                                                                  |
| `percentile(p,...)`    | Percentile                                            | `percentile(50,1,3,5)=3`                                                                                                                                                           |
| `quantile(p,...)`      | Quantile (`p` in \[0,1\])                             | `quantile(1/4,1,2,3,4,5,6,7)=2.5`                                                                                                                                                  |
| `var(x...)`            | Population variance (divide by _n_)                   | `var(1,2,3)=0.666666666667`                                                                                                                                                        |
| `vars(x...)`           | Sample variance (unbiased, divide by _n-1_)           | `vars(1,2,3)=1`                                                                                                                                                                    |
| `stddev(x...)`         | Standard deviation                                    | `stddev(1,2,3)=0.816496580928`                                                                                                                                                     |
| `stddevs(x...)`        | Sample standard deviation (unbiased)                  | `stddevs(1,2,3)=1`                                                                                                                                                                 |
| `geomean(...)`         | Geometric mean                                        | `geomean(1,4,1/32)=0.5`                                                                                                                                                            |
| `harmmean(...)`        | Harmonic mean                                         | `harmmean(1,2,6)=1.8`                                                                                                                                                              |
| `mad(...)`             | Median Absolute Deviation (MAD)                       | `mad(1,1,2,2,4)=1`                                                                                                                                                                 |
| `skew(...)`            | Skewness (symmetry measure)                           | `skew(1,2,3,4,5)=0` (symmetric)<br>`skew(0,0,0,0,10) > 0` (right-tailed)<br>`skew(-10,0,0,0,0) < 0` (left-tailed)<br>`skew(-2,-1,0,1,2)=0`                                         |
| `kurtp(...)`           | Population Excess kurtosis (normal distribution -> 0) | `kurtp(-1,-0.5,0,0.5,1) < 0 (-1.3)` (light tails)<br>`kurtp(0,0,0,0,0,10) >> 0` (heavy tails)<br>`kurtp(-2,-1,0,1,2) < 0` (uniform-like)<br>`kurtp(0,0,0,0,0,-2,2) > 0` (outliers) |
| `kurts(...)`           | Sample Excess kurtosis, `n\geq 5`                     | `kurts(-1,-0.5,0,0.5,1)= -7.09`                                                                                                                                                    |
| `cv(...)`              | Coefficient of variation = stddev/mean                | `cv(10,10,10)=0`                                                                                                                                                                   |
| `stderr(...)`          | Standard error = stddev/sqrt(n), `n\geq 2`            | `stderr(1,2,3)=0.471...`                                                                                                                                                           |
| `zscore(x, mu, sigma)` | (x - mu) / sigma                                      | `zscore(5,3,1)=2`                                                                                                                                                                  |
| `iqr(...)`             | Interquartile range (Q3 - Q1)                         | `iqr(1,2,3,4)=2`                                                                                                                                                                   |
| `trimmean(p,...)`      | Trimmed mean (discard fraction _p_ from both ends)    | `trimmean(0.2,1,2,100,3,4)=2.5`                                                                                                                                                    |
| `winsor(p,...)`        | Winsorization (clip fraction _p_ from both ends)      | `winsor(0.2,1,2,100,3,4)=?`                                                                                                                                                        |

---

### Correlation / Covariance

| Function Name             | Description                             | Example Inputs/Outputs            |
| ------------------------- | --------------------------------------- | --------------------------------- |
| `cov(x...,y...)`          | Covariance                              | `cov(1,2,3,2,4,6)=0.666...`       |
| `corr(x...,y...)`         | Correlation coefficient                 | `corr(1,2,3,2,4,6)=1`             |
| `corrspearman(x...,y...)` | Spearman correlation (rank correlation) | `corrspearman(1,2,3, 10,20,30)=1` |

---

### Numerical Computation (Floating-Point Utilities)

| Function Name | Description                                  | Example Inputs/Outputs |
| ------------- | -------------------------------------------- | ---------------------- |
| `hypot(x,y)`  | √(x² + y²) (robust against cancellation)     | `hypot(3,4)=5`         |
| `fma(a,b,c)`  | Compute a\*b + c with a single rounding step | `fma(10,10,2)=102`     |
| `rms(...)`    | Root mean square                             | `rms(1,-1,1,-1)=1`     |

---

### Complex Numbers

| Function Name                       | Description                    | Example                |
| ----------------------------------- | ------------------------------ | ---------------------- |
| `re(z)`, `real(z)`                  | Real part                      | `re(3+4I)=3`           |
| `im(z)`, `imag(z)`                  | Imaginary part                 | `im(3+4I)=4`           |
| `arg(z)`                            | Argument (degree)              | `arg(1+I)=45`          |
| `conj(z)`                           | Complex conjugate              | `conj(3+4I)=3-4I`      |
| `polar(r,theta)`<br>`rect(r,theta)` | Polar form r∠θ (degree)        | `polar(2,60)=1+1.732I` |
| `cis(x)`                            | cos(x) + i sin(x) (degree)     | `cis(60)=0.5+0.866I`   |
| `mag(z)`                            | Alias of `abs(z)`              | `mag(3+4I)=5`          |
| `proj(z)`                           | Riemann sphere projection      | `proj(3+4I)=3+4I`      |
| `unit(z)`, `csgn(z)`                | Unit complex number `z/abs(z)` | `unit(3+4I)=0.6+0.8I`  |

---

#### Geometry / Vectors (2D)

| Function Name           | Description                 | Example Inputs/Outputs |
| ----------------------- | --------------------------- | ---------------------- |
| `norm(x,y)`             | Vector magnitude            | `norm(3,4)=5`          |
| `dot(x1,y1,x2,y2)`      | Dot product                 | `dot(1,0,0,1)=0`       |
| `cross(x1,y1,x2,y2)`    | Cross product (Z component) | `cross(1,0,0,1)=1`     |
| `lerp(a,b,t)`           | Linear interpolation        | `lerp(0,10,0.5)=5`     |
| `distance(x1,y1,x2,y2)` | Distance                    | `distance(0,0,3,4)=5`  |

---

### Random Numbers

| Function Name      | Description                                           | Example Inputs/Outputs                |
| ------------------ | ----------------------------------------------------- | ------------------------------------- |
| `rand()`           | Uniform random in \[0,1)                              | `rand()=0.37...`                      |
| `rand(hi)`         | Uniform random in \[0,hi)                             | `rand(5)=3.37...` (returns 0 if hi=0) |
| `rand(lo,hi)`      | Uniform random in \[lo,hi)                            | `rand(-2,2)=-0.71...`                 |
| `randint()`        | Random integer 0 or 1                                 | `randint()=0`                         |
| `randint(a)`       | Random integer in \[0,a] if `a>0`, or \[a,0] if `a<0` | `randint(5)=2`                        |
| `randint(a,b)`     | Random integer in \[a,b]                              | `randint(1,6)=3` (dice)               |
| `choice(a, b,...)` | Randomly select one element                           | `choice(2,3,5,7,9,11,13,17,19)=17`    |
| `randn()`          | Normal random N(0,1)                                  | `randn()=-0.23...`                    |
| `randn(mu)`        | Normal random N(mu,1)                                 | `randn(10)=9.61...`                   |
| `randn(mu,sigma)`  | Normal random N(mu,sigma)                             | `randn(0,2)=1.74...`                  |

---

### Strength of Materials (Stress, Strain, Elasticity)

| Function Name      | Description            | Example                   |
| ------------------ | ---------------------- | ------------------------- |
| `stress(F,A)`      | Stress: σ = F/A        | `stress(1000,10)=100`     |
| `strain(dL,L)`     | Strain: ε = dL/L       | `strain(0.1,100)=0.001`   |
| `young(sigma,eps)` | Young’s modulus: E=σ/ε | `young(200,0.001)=200000` |

---

### Section Properties (I, Z, J)

| Function Name          | Description                            | Example                        |
| ---------------------- | -------------------------------------- | ------------------------------ |
| `moment_rect(b,h)`     | Rectangular area moment: I=b\*h³/12    | `moment_rect(10,20)=6666.6`    |
| `moment_circle(d)`     | Circular area moment: I=π\*d⁴/64       | `moment_circle(10)=490.87`     |
| `sectionmod_rect(b,h)` | Rectangular section modulus: Z=b\*h²/6 | `sectionmod_rect(10,20)=666.6` |
| `sectionmod_circle(d)` | Circular section modulus: Z=π\*d³/32   | `sectionmod_circle(10)=98.17`  |
| `torsion_J_circle(d)`  | Polar moment of inertia: J=π\*d⁴/32    | `torsion_J_circle(10)=981.74`  |
| `polarZ_circle(d)`     | Polar section modulus: Zp=π\*d³/16     | `polarZ_circle(10)=196.35`     |

---

### Fasteners / Friction (Bolts, Torque)

| Function Name                | Description                     | Example                                  |
| ---------------------------- | ------------------------------- | ---------------------------------------- |
| `bolt_stress(F,d)`           | Bolt stress: σ = F/(πd²/4)      | `bolt_stress(10000,10)=127.3`            |
| `torque_from_preload(F,d,K)` | Tightening torque: T = K F d    | `torque_from_preload(10000,0.01,0.2)=20` |
| `preload_from_torque(T,d,K)` | Inverse conversion: F = T/(K d) | `preload_from_torque(20,0.01,0.2)=10000` |
| `friction(mu,N)`             | Friction force: F = μ N         | `friction(0.2,100)=20`                   |

---

### Utilities

| Function Name      | Description | Example Inputs/Outputs |
| ------------------ | ----------- | ---------------------- |
| `clamp(x, lo, hi)` | Clamp       | `clamp(5,0,10)=5`      |

---

#### History

```text
In[n], Out[n], %, %%, %%%...
```

---

#### System Functions

| Function Name | Description                                   |
| ------------- | --------------------------------------------- |
| `Clear[]`     | clears history and resets`In[n] / Out[n]`to 1 |
| `Exit`        | terminates the program                        |

---

## Operator Precedence Table

The following table shows **operator precedence (high → low)** in this engine:

| Priority | Operator / Syntax       | Type                      | Associativity | Example             |
| -------- | ----------------------- | ------------------------- | ------------- | ------------------- |
| 1        | `!`                     | Factorial                 | Postfix       | `3! = 6`            |
| 2        | Unary `+ -`             | Unary                     | Right         | `--5 = 5`           |
| 3        | `^`                     | Power                     | **Right**     | `2^3^2 = 512`       |
| 4        | Implicit multiplication | Multiplication            | Left          | `2Pi`, `(2+3)(4+5)` |
| 5        | `* /`                   | Multiplication / division | Left          | `10/2*3 = 15`       |
| 6        | `+ -`                   | Addition / subtraction    | Left          | `1-2-3 = -4`        |
| 7        | `< <= > >= ==`          | Comparison                | Left          | `2 < 3 = 1`         |

### Notes on Precedence

- `!` is the **highest-priority postfix operator**
- Unary operators bind tighter than power, but power is right-associative
- Implicit multiplication has higher precedence than explicit `* /`
- Comparison operators are always evaluated last

---

### Comparison Operators

| Operator    | Meaning               |
| ----------- | --------------------- |
| `< <= > >=` | Relational comparison |
| `==`        | Equality comparison   |

- Results are **1 (true) / 0 (false)**

---

### Complex Number Behavior

- `sqrt(-1) = I`
- `exp(I)` follows Euler’s formula
- Real numbers are automatically promoted to complex numbers when required
- `+0I` and `-0I` both evaluate to zero. Therefore, `arg(-1+0I) = 180` and `arg(-1-0I) = 180` are obtained respectively.

---

### Error Handling

The following conditions are explicitly treated as **Errors**:

- Mismatched parentheses
- Invalid number of arguments
- Domain errors (`log(-1)`, `tan(90)`, etc.)
- Invalid syntax (`2**3`, `1+*2`, etc.)

#### Error Categories

| Category           | Description                   | Example        | Error Message               |
| ------------------ | ----------------------------- | -------------- | --------------------------- |
| Syntax error       | Invalid token order           | `2**3`, `1+*2` | `syntax error`              |
| Parenthesis error  | Unbalanced parentheses        | `(1+2`         | `syntax error`              |
| Argument error     | Invalid argument count        | `sin()`        | `function argument missing` |
| Type error         | Complex value in real-only op | `(-3)!`        | `factorial: negative value` |
| Domain error       | Mathematically undefined      | `log(-1)`      | `domain error`              |
| Division by zero   | Denominator is zero           | `1/0`          | `division by zero`          |
| Undefined behavior | Unsupported syntax            | `sin30`        | `unknown identifier`        |

#### Error Policy

- Evaluation stops immediately upon error
- Errors do not propagate as values
- IEEE754 results (`inf`, `-inf`) are distinguished from errors
- Distinction between returning `inf` and `result is infinite` is currently ambiguous (needs improvement)

---

#### Singularities and Exceptions

- When an error occurs, subsequent evaluations are not performed
- Errors are not propagated as values but are displayed immediately

---

## Usage Examples

```text
In [1] := sin(30)^2 + cos(30)^2
Out[1] := 1

In [2] := sqrt(-4)
Out[2] := 2I

In [3] :=

In [3] := (2+3)(4+5)
Out[3] := 45

In [4] := 10/2+3
Out[4] := 8

In [5] := --3!-(-(3!))
Out[5] := 12

In [6] := gcd(8,4)
Out[6] := 4

In [7] := pow(4,0.5)
Out[7] := 2

In [8] := sin(90)
Out[8] := 1

In [9] := Out[3]+2
Out[9] := 47

In [10] := In[2+3]*5
Out[10] := 60


In [11] := %%
Out[11] := 47

In [12] := Clear[]

In [1] :=

```

## Design Notes

### Angle System

- All trigonometric functions use **degrees** by default
- Degrees were chosen for usability in practical environments
- Use `RtoD` if needed
- Native radian support is a future extension

---

### Numerical Precision

- Internally uses `double`
- Since the internal representation is double precision, it can only handle values within the range of 1.7E±308.
  Values exceeding ±1.7E+308 will result in an overflow or infinity error.

- Subject to IEEE754 rounding errors

```text
atan(0.57735026919) = 30.000000000016
```

Displayed values and history are rounded to 12 decimal places.

---

### Intentional Limitations

- No variable assignment
- No user-defined functions
- No symbolic computation

---

### Known Behavior

- `0^-1` → `±inf` (IEEE754 compliant)
- `-0` may be displayed
- Expressions inside `In[] / Out[]` are undefined

---

### Specification behavior

- Spaces are generally ignored, but between numbers they are treated as implicit multiplication. `3 2` -> `3*2`
- In principle, complex solutions are returned as complex numbers, but `cbrt` returns real solutions. (`(-8)^(1/3)=1+1.732050807569I`, `cbrt(-8)=-2`)

---

### Future Plans

- [N]Base Support (0x, 0b, add(), shift())
- AngleMode Support (sin(30deg), setAngle=RADIAN)

---

## Grammar Definition (EBNF)

```ebnf
expression ::= comparison ;

comparison ::= additive
| additive ( ("<" | "<=" | ">" | ">=" | "==") additive ) ;

additive ::= multiplicative { ("+" | "-") multiplicative } ;

multiplicative ::= power { ("*" | "/") power | implicit_mul power } ;

implicit_mul ::= /* adjacency rule */ ;

power ::= unary [ "^" power ] ; // right associative

unary ::= { "+" | "-" } postfix ;

postfix ::= primary { "!" } ;

primary ::= number
| constant
| function_call
| "(" expression ")" ;

function_call ::= identifier "(" [ arguments ] ")" ;

arguments ::= expression { "," expression } ;

number ::= integer | float ;

constant ::= "Pi" | "E" | "Phi" | "I" ;

identifier ::= letter { letter | digit } ;

integer ::= digit { digit } ;

float ::= digit { digit } "." digit { digit } ;

letter ::= "A"…"Z" | "a"…"z" ;

digit ::= "0"…"9" ;
```

---

## License

BSD 3-Clause License

Copyright (c) 2021–2026 mmKreutzef

If you use this software in academic papers or commercial products, please clearly state that fact in your documentation or publications.

A short notice would be appreciated.(and feel delight)

---

## Tests

This project includes **600 automated tests**, covering:

- Operator precedence and associativity
- Error handling
- Complex arithmetic
- Boundary and edge cases

---

## Notes

This project aims to balance **strict operator semantics** with **practical expression evaluation**, ensuring it works reliably in design and manufacturing environments of any scale.

---

### **Disclaimer**

This software is provided “as is” without any warranty, express or implied, including but not limited to warranties of merchantability, fitness for a particular purpose, or non-infringement.  
The author or copyright holder shall not be liable for any claims, damages, or other liabilities arising out of or in connection with the software, whether based on contract, tort, or otherwise.

By using this software, you acknowledge and automatically agree that all risks associated with its use are your responsibility.
The author assumes no responsibility for any damage resulting from data loss, system malfunction, or other consequences arising from the use of this software.

---

## Requests / Contributions

Please submit them via Github.I love to hear any suggestions for implementation, especially those needed in manufacturing or design environments.
