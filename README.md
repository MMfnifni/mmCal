# mmCalculator – Mechanical & Manufacturing Calculator

© 2021–2026 mmKreutzef (aka Daiki.NIIMI) Licensed under the BSD 3-Clause License

[English](README.md) | [日本語](README.ja.md)

# Overview

mmCalculator is a **high-performance mathematical calculator** designed for engineering, manufacturing, and technical work.

It goes beyond a simple calculator by supporting:

- Continuous calculations using previous results
- Variables and user-defined functions
- Complex numbers, vectors, and matrices
- Advanced numerical and engineering functions

Unlike heavyweight systems, it remains:

> **Fast, lightweight, and immediately usable — a practical tool for real work**

## Enough theory — let’s dive straight into practical examples!

```txt
In [1] := (2+3)(4+5) + 2Pi
Out[1] := 51.28318530718

In [2] := a := hypot(3,4)
Out[2] := <defined: a := 5>

In [3] := f(x) := x^2 + 2x + 1
Out[3] := <defined: f(x)>

In [4] := f(a)
Out[4] := 36

In [5] := round(ln(E^5) + log10(1000) + log2(8), 10)
Out[5] := 11

In [6] := fma(10^16, 3, -10^16*3) + expm1(1) - (E-1)
Out[6] := 0

In [7] := sin(30)^2 + cos(30)^2
Out[7] := 1

In [8] := tan(45) + cot(Pi/4 rad) + sec(60deg) + csc(Pi/6 rad)
Out[8] := 6

In [9] := atan2(1,1) + atan2(1,-1) + atan2(-1,-1) + atan2(-1,1)
Out[9] := 360

In [10] := Out[5]*a - %%
Out[10] := 19

In [11] := %% / %
Out[11] := 18.947368421053

In [12] := 3!^2 / round(%)
Out[12] := 2

In [13] := (-8)^(1/3)
[WARN] pos=4 (-a)^(p/q) : principal value only; other branches may exist
Out[13] := 1+1.732050807569I

In [14] := exp(I Pi) + 1
Out[14] := 0

In [15] := polar(2,60) + cis(60)
Out[15] := 1.5+2.598076211353I

In [16] := abs(%) + arg(%) + re(%) + im(%)
Out[16] := 67.098076211353

In [17] := conj(%%) + unit(%%) - proj(%%)
Out[17] := 0.5-4.330127018922I

In [18] := mean(3,5,7,11,13,17) + median(1,3,5) + mode(1,2,2,3)
Out[18] := 14.3333333333333

In [19] := stddevs(3,5,7,11,13,17) + vars(1,2,3) + iqr(1,2,3,4)
Out[19] := 7.778888771954

In [20] := gcd(84,120) + lcm(6,8) + perm(10,3) + comb(10,3)
Out[20] := 876

In [21] := ifft(fft({1+I,2-I,3+2I,4-3I}))
Out[21] := {1+I,2-I,3+2I,4-3I}

In [22] := mlsqr({{1,2},{3,4}},{1,2})
Out[22] := {0, 0.5}

In [23] := vdot({1,2},{3,4}) + vmanhattan({0,0},{3,4}) + veuclidean({0,0},{3,4})
Out[23] := 23
```

## Key Points Express

- Implicit multiplication works (natural input)
- Support for user-defined variables and functions
- Use `%` or `Out[n]` to reuse past history in calculations
- Supports complex numbers
- Angles default to degrees but can be specified with `rad`, etc.
- Extensive set of functions and constants!(Matrices, Vectors also acceptable)

# Detailed implementation

The following is closer to a **technical note** and is not organized. Please bear with me for a moment.

## Numeric Types

- **Real numbers**: IEEE 754 `double` (displayed with up to 12 decimal digits)
- **Complex numbers**: `a + bI`
- **Array** `{}` Arrays can contain real numbers, complex numbers, or other arrays.

### Constants

All constants start with a capital letter.

| Name  | Description                    | Expansion value                         |
| ----- | ------------------------------ | --------------------------------------- |
| `Pi`  | π (pi)                         | `3.14159265358979323846264338327950288` |
| `Tau` | 2π (τ)                         | `6.283185307179586476925286766559006`   |
| `E`   | Base of natural logarithm      | `2.7182818284590452353602874713526625`  |
| `Phi` | Golden ratio                   | `1.618033988749894848204586834365638`   |
| `NA`  | Avogadro constant              | `6.02214076e23`                         |
| `I`   | Imaginary unit                 | `0+I`                                   |
| `ESP` | Internal convergence tolerance | `1e-12`                                 |

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

### Unary Operators

- Unary `+` and `-` may be chained arbitrarily.

```text
---5 = -5
```

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

### User Variables and Function Definitions

You can define variables or functions using `:=`.
Variable and function names must comply with the following constraints:

- They must not conflict with built-in functions or variables (e.g., `Pi`, `ave(x)`)
- You cannot use a name defined as a variable for a function, or vice versa
- Names must consist only of letters, digits, and underscores, and must begin with a letter

Examples

```
In [1] := x:=2
Out[1] := <defined: x := 2>

In [2] := f(y):=y^y
Out[2] := <defined: f(y)>

In [3] := f (3)
Out[3] := 27

In [4] := f(x)
Out[4] := 4
```

Regarding the function definition `f(x)`, `x` accepts any string that follows the naming conventions.
However, if `x` has been defined previously, it is recognized as a function variable rather than a user-defined variable within the function definition.
Additionally, by default, variables are overwritable, while functions are not.

```
In [6] := h(y):=x^y
Out[6] := <defined: h(y)>

In [7] := h(2)
Out[7] := 4

In [8] := x:=3
Out[8] := <redefined: x := 3 (was 2)>

In [9] := h(2)
Out[9] := 9
```

Lists and removal of user-defined variables and functions are provided by commands such as `:unset`.
For details, see “System Functions”.

### Functions

Below is a complete list of all currently implemented functions. For each function, the function name, a brief description, and representative input/output examples are provided.

### Basic Math (Algebra, Exponentials, Logarithms, Rounding)

| Function Name    | Description                               | Example Inputs/Outputs           |
| ---------------- | ----------------------------------------- | -------------------------------- |
| `abs(x)`         | Absolute value (supports complex numbers) | `abs(-3) = 3`<br>`abs(3+4I) = 5` |
| `pow(x,y)`       | Power function                            | `pow(2,10)=1024`                 |
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
| `fact(x)`        | Factorial (non-negative integers only)    | `fact(5)=120`                    |
| `mod(x,y)`       | Surplus                                   | `mod(11,3)=2`                    |
| `nextpow2(x)`    | The smallest n such that 2^n >= x         | `nextpow2(9)=4`                  |

### Special Functions

| Function Name   | Description                          | Example Inputs/Outputs        |
| --------------- | ------------------------------------ | ----------------------------- |
| `gamma(x)`      | Gamma function                       | `gamma(5)=24`                 |
| `lgamma(x)`     | Log-gamma function                   | `lgamma(5)=3.178053830348...` |
| `erf(x)`        | Error function                       | `erf(1)=0.842700792949...`    |
| `erfc(x)`       | Complementary error fn               | `erfc(1)=0.157299207050...`   |
| `digamma(x)`    | Digamma function                     | `digamma(1)=-0.57721...`      |
| `trigamma(x)`   | Trigamma function                    | `trigamma(1)=1.64493...`      |
| `zeta(s)`       | Riemann zeta function (`s>1`)        | `zeta(2)=1.64493...`          |
| `beta(x,y)`     | Beta function                        | `beta(2,3)=0.0833333...`      |
| `ibeta(a,b,x)`  | Regularized incomplete beta function | `ibeta(2,3,0.5)=0.6875`       |
| `betaln(x,y)`   | Beta-log function                    | `betaln(2,3)=-2.48490...`     |
| `binom(x,y)`    | Generalized binomial coefficient     | `binom(1/2,2)=-0.125`         |
| `fallfact(x,n)` | Falling factorial                    | `fallfact(5,3)=60`            |
| `risefact(x,n)` | Rising factorial                     | `risefact(5,3)=210`           |

#### Angle Conversion

| Function Name | Description     | Example Inputs/Outputs |
| ------------- | --------------- | ---------------------- |
| `DtoR(x)`     | Degree → Radian | `DtoR(180)=Pi`         |
| `DtoG(x)`     | Degree → Grad   | `DtoG(90)=100`         |
| `RtoD(x)`     | Radian → Degree | `RtoD(Pi)=180`         |
| `RtoG(x)`     | Radian → Grad   | `RtoG(Pi)=200`         |
| `GtoD(x)`     | Grad → Degree   | `GtoD(200)=180`        |
| `GtoR(x)`     | Grad → Radian   | `GtoR(200)=Pi`         |

### Trigonometric Functions (Degree Mode / degree)

| Function Name | Description                                     | Example Inputs/Outputs |
| ------------- | ----------------------------------------------- | ---------------------- |
| `sin(x)`      | Sine                                            | `sin(30)=0.5`          |
| `cos(x)`      | Cosine                                          | `cos(60)=0.5`          |
| `tan(x)`      | Tangent                                         | `tan(45)=1`            |
| `cot(x)`      | Cotangent                                       | `cot(45)=1`            |
| `sec(x)`      | Secant                                          | `sec(60)=2`            |
| `csc(x)`      | Cosecant                                        | `csc(30)=2`            |
| `asin(x)`     | Inverse sine                                    | `asin(0.5)=30`         |
| `acos(x)`     | Inverse cosine                                  | `acos(0.5)=60`         |
| `atan(x)`     | Inverse tangent                                 | `atan(1)=45`           |
| `atan2(y,x)`  | Quadrant-aware arctangent                       | `atan2(1,1)=45`        |
| `csch(x)`     | Hyperbolic secant function `1/sinh(x)`          | `csch(1)=0.8509...`    |
| `sech(x)`     | Hyperbolic tangent function `1/cosh(x)`         | `sech(0)=1`            |
| `coth(x)`     | Hyperbolic cotangent function `cosh(x)/sinh(x)` | `coth(1)=1.313...`     |

### Hyperbolic Functions

| Function Name | Description             | Example Inputs/Outputs |
| ------------- | ----------------------- | ---------------------- |
| `sinh(x)`     | Hyperbolic sine         | `sinh(1)=1.175...`     |
| `cosh(x)`     | Hyperbolic cosine       | `cosh(0)=1`            |
| `tanh(x)`     | Hyperbolic tangent      | `tanh(0)=0`            |
| `asinh(x)`    | Inverse hyperbolic sine | `asinh(1)=0.881...`    |
| `acosh(x)`    | Inverse hyperbolic cos  | `acosh(1)=0`           |
| `atanh(x)`    | Inverse hyperbolic tan  | `atanh(0.5)=0.549...`  |

### Signal

| Function Name  | Description                    | Example Inputs/Outputs                     |
| -------------- | ------------------------------ | ------------------------------------------ |
| `fft(...)`     | Fast Fourier Transform         | `fft({1,1,1,1,1,1,1,1})={8,0,0,0,0,0,0,0}` |
| `ifft(...)`    | Inverse Fast Fourier Transform | `ifft({4,0,0,0})={1, 1, 1, 1}`             |
| `dft(...)`     | Discrete Fourier Transform     | `dft({1,2,3,4})={10,-2+2I,-2,-2-2I}`       |
| `hilbert(...)` | Hilbert transform              | `hilbert({1,0,-1,0})={1,I,-1,-I}`          |
| `sinc(x)`      | `sin(x)/x`                     | `sinc(90)=0.111...`                        |
| `cosc(x)`      | `(1-cos(x))/x`                 | `cosc(60)=0.008333...`                     |
| `tanc(x)`      | `tan(x)/x`                     | `tanc(0)=1`                                |
| `sinhc(x)`     | `sinh(x)/x`                    | `sinhc(3)=3.339...`                        |
| `tanhc(x)`     | `tanh(x)/x`                    | `tanhc(30)=0.03333...`                     |
| `expc(x)`      | `(exp(x)-1)/x`                 | `expc(1)=1.7182818...`                     |

- Since sin/cos/tan assume degree input, these also take degree input.
- When `x` is extremely small, loss-of-significance countermeasures are applied.
- When `fft` and `ifft` are not `N^2`, they are processed as `dft`.

| Function Name | Description    | Example Inputs/Outputs |
| ------------- | -------------- | ---------------------- |
| `sincR(x)`    | `sin(x)/x`     | `sincR(1e-8)≈1`        |
| `coscR(x)`    | `(1-cos(x))/x` | `coscR(1e-8)≈5e-9`     |
| `tancR(x)`    | `tan(x)/x`     | `tancR(1e-8)≈1`        |

- These interpret `x` as radians.
- For the signal-processing-style “sinc”, this is generally the one you want.

### Number Theory / Combinatorics

| Function Name                     | Description                                            | Example Inputs/Outputs         |
| --------------------------------- | ------------------------------------------------------ | ------------------------------ |
| `gcd(a,b,...)`                    | Greatest common divisor                                | `gcd(12,18,24)=6`              |
| `lcm(a,b,...)`                    | Least common multiple                                  | `lcm(6,8)=24`                  |
| `perm(n,r)`                       | Permutations [(returns 0 if r<0 or r>n)]               | `perm(5,2)=20`                 |
| `comb(n,r)`                       | Combinations [(returns 0 if r<0 or r>n)]               | `comb(5,2)=10`                 |
| `isprime(x)`                      | Prime number determination<br>(0: non-prime, 1: prime) | `isprime(9973)=1`              |
| `nextprime(n)`                    | The smallest prime number exceeding                    | `nextprime(11) = 13`           |
| `prevprime(n)`                    | The largest prime number less than _n_                 | `prevprime(11) = 7`            |
| `fib(x)`                          | Fibonacci sequence                                     | `fib(25)=75025`                |
| `factorint(n)`, `primefactors(n)` | Prime factorization (returns a list of prime factors)  | `factorint(84) = {2, 2, 3, 7}` |

### Statistics (Descriptive Statistics / Aggregation)

| Function Name          | Description                                           | Example Inputs/Outputs                                                                                                                                                                                |
| ---------------------- | ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sum(...)`             | Sum                                                   | `sum(1,2,3)=6`                                                                                                                                                                                        |
| `mean(...)`, `ave()`   | Mean                                                  | `mean(1,2,3)=2`                                                                                                                                                                                       |
| `min(...)`             | Minimum                                               | `min(3,1,2)=1`                                                                                                                                                                                        |
| `max(...)`             | Maximum                                               | `max(3,1,2)=3`                                                                                                                                                                                        |
| `prod(...)`            | Product                                               | `prod(2,3,4)=24`                                                                                                                                                                                      |
| `median(...)`          | Median                                                | `median(1,3,5)=3`                                                                                                                                                                                     |
| `mode(...)`            | Mode                                                  | `mode(1,2,2,3)=2`                                                                                                                                                                                     |
| `percentile(p,...)`    | Percentile                                            | `percentile(50,1,3,5)=3`                                                                                                                                                                              |
| `quantile(p,...)`      | Quantile (`p` in \[0,1\])                             | `quantile(1/4,1,2,3,4,5,6,7)=2.5`                                                                                                                                                                     |
| `var(x...)`            | Population variance (divide by _n_)                   | `var(1,2,3)=0.666666666667`                                                                                                                                                                           |
| `vars(x...)`           | Sample variance (unbiased, divide by _n-1_)           | `vars(1,2,3)=1`                                                                                                                                                                                       |
| `stddev(x...)`         | Standard deviation                                    | `stddev(1,2,3)=0.816496580928`                                                                                                                                                                        |
| `stddevs(x...)`        | Sample standard deviation (unbiased)                  | `stddevs(1,2,3)=1`                                                                                                                                                                                    |
| `geomean(...)`         | Geometric mean                                        | `geomean(1,4,1/32)=0.5`                                                                                                                                                                               |
| `harmmean(...)`        | Harmonic mean                                         | `harmmean(1,2,6)=1.8`                                                                                                                                                                                 |
| `mad(...)`             | Median Absolute Deviation                             | `mad(1,1,2,2,4)=1`                                                                                                                                                                                    |
| `madR(v)`              | Mean Absolute Deviation                               | `madR(1,2,3) = 0.666...`                                                                                                                                                                              |
| `skew(...)`            | Skewness (symmetry measure)                           | `skew(1,2,3,4,5)=0` (symmetric)<br>`skew(0,0,0,0,10) > 0` (right-tailed)<br>`skew(-10,0,0,0,0) < 0` (left-tailed)<br>`skew(-2,-1,0,1,2)=0`                                                            |
| `kurtp(...)`           | Population Excess kurtosis (normal distribution -> 0) | `kurtp(-1,-0.5,0,0.5,1) < 0 (-1.3)` (light tails)<br>`kurtp(0,0,0,0,0,10) >> 0 (1.2)` (heavy tails)<br>`kurtp(-2,-1,0,1,2) < 0 (-1.3)` (uniform-like)<br>`kurtp(0,0,0,0,0,-2,2) > 0 (0.5)` (outliers) |
| `kurts(...)`           | Sample Excess kurtosis, `n\geq 5`                     | `kurts(-1,-0.5,0,0.5,1)= -7.09`                                                                                                                                                                       |
| `cv(...)`              | Coefficient of variation = stddev/mean                | `cv(10,10,10)=0`                                                                                                                                                                                      |
| `stderr(...)`          | Standard error = stddev/sqrt(n), `n\geq 2`            | `stderr(1,2,3)=0.471...`                                                                                                                                                                              |
| `zscore(x, mu, sigma)` | (x - mu) / sigma                                      | `zscore(5,3,1)=2`                                                                                                                                                                                     |
| `iqr(...)`             | Interquartile range (Q3 - Q1)                         | `iqr(1,2,3,4)=2`                                                                                                                                                                                      |
| `trimmean(p,...)`      | Trimmed mean (discard fraction _p_ from both ends)    | `trimmean(0.2,1,2,100,3,4)=2.5`                                                                                                                                                                       |
| `winsor(p,...)`        | Winsorization (clip fraction _p_ from both ends)      | `winsor(0.2,1,2,100,3,4)=3`                                                                                                                                                                           |

### Numerical Calculus (Preliminary Implementation)

| Function Name                            | Description                                                | Input and Output Examples         |
| ---------------------------------------- | ---------------------------------------------------------- | --------------------------------- |
| `diff(f,x)`                              | Numerical first-order derivative (central difference)      | `diff(sin,30)=0.015114...`        |
| `diff(f,x,h)`                            | Numerical first-order derivative with specified step size  | `diff(exp,1,1e-6)=2.71828...`     |
| `diff2(f,x)`                             | Numerical second-order derivative                          | `diff2(sin,0)=-0.000304... `      |
| `diff2(f,x,h)`                           | Numerical second-order derivative with specified step size | `diff2(exp,1,1e-4)=2.71828...`    |
| `simpson(f,a,b)`, `integrate(f,a,b)`     | Numerical integration using Simpson's rule                 | `simpson(sin,0,180)=114.59...`    |
| `simpson(f,a,b,n)`, `integrate(f,a,b,n)` | Simpson's rule with specified number of segments           | `simpson(exp,0,1,200)=1.718...`   |
| `trapz(f,a,b)`                           | Numerical integration using the trapezoidal rule           | `trapz(exp,0,1)=1.718...`         |
| `trapz(f,a,b,n)`                         | Trapezoidal rule with specified number of segments         | `trapz(sin,0,180,1000)=114.59...` |

### Financial Functions

| Function Name                           | Description                                           | Input and Output Example                      |
| --------------------------------------- | ----------------------------------------------------- | --------------------------------------------- |
| `fin_fv(rate,n,pmt)`                    | Future Value (FV)                                     | `fin_fv(0.05,10,100)=1257.789...`             |
| `fin_pv(rate,nper,pmt)`                 | Present Value (PV)                                    | `fin_pv(0.05,10,100)=772.173...`              |
| `fin_pmt(rate,nper,pv)`                 | Payment Amount (PMT)                                  | `fin_pmt(0.05,10,1000)=129.504...`            |
| `fin_total_payment(rate,nper,pv)`       | Total Payment                                         | `fin_total_payment(0.05,10,1000)=1295.045...` |
| `fin_ppmt(rate,nper,per,pv)`            | Principal Payment Amount (PPMT, Period-specified)     | `fin_ppmt(0.05,10,3,1000)=92.455...`          |
| `fin_ipmt(rate,nper,per,pv)`            | Interest Payment Amount (IPMT, Period Specified)      | `fin_ipmt(0.05,10,3,1000)=37.049...`          |
| `fin_npv(rate,cf...)`                   | Net Present Value (NPV)                               | `fin_npv(0.1,100,-50,-50)=12.021... `         |
| `fin_irr(cf...)`                        | Internal Rate of Return (IRR, Brent)                  | `fin_irr(-100,50,60)=0.063...`                |
| `fin_mirr(rate,cf...)`                  | Modified Internal Rate of Return (MIRR)               | `fin_mirr(0.1,-100,60,60)=...`                |
| `fin_rate(nper,pmt,pv)`                 | Interest Rate Calculation (RATE, Brent)               | `fin_rate(10,100,1000)=0.05`                  |
| `fin_nper(rate,pmt,pv)`                 | Period Calculation (NPER)                             | `fin_nper(0.05,100,1000)=10`                  |
| `fin_cumipmt(rate,nper,pmt,start,end)`  | Cumulative interest over specified period (CUMIPMT)   | `fin_cumipmt(0.05,10,100,1,5)=210.562...`     |
| `fin_cumprinc(rate,nper,pmt,start,end)` | Cumulative principal over specified period (CUMPRINC) | `fin_cumprinc(0.05,10,100,1,5)=289.437...`    |
| `fin_effective_rate(nominal,npery)`     | Effective interest rate (EFFECTIVE RATE)              | `fin_effective_rate(0.06,12)=0.061...`        |
| `fin_nominal_rate(effective,npery)`     | Nominal Rate (NOMINAL RATE)                           | `fin_nominal_rate(0.0617,12)=0.060...`        |
| `fin_cagr(start,end,n)`                 | Compound Annual Growth Rate (CAGR)                    | `fin_cagr(1000,2000,10)=0.071...`             |

### Correlation / Covariance / Uncategorized(ToDo: To be organized)

| Function Name             | Description                                             | Example Inputs/Outputs            |
| ------------------------- | ------------------------------------------------------- | --------------------------------- |
| `cov(x...,y...)`          | Covariance                                              | `cov(1,2,3,2,4,6)=0.666...`       |
| `corr(x...,y...)`         | Correlation coefficient                                 | `corr(1,2,3,2,4,6)=1`             |
| `corrspearman(x...,y...)` | Spearman correlation (rank correlation)                 | `corrspearman(1,2,3, 10,20,30)=1` |
| `polylog(s, z)`           | Polylog calculation using a multiple series (`\|z\|<1`) | `polylog(2, 0.5)=0.582240526465`  |
| `totient(n)`              | Euler's totient function φ(n)                           | `totient(9)=6`                    |
| `winsorR(p,...)`          | Returns the winsorized data series                      | `winsorR(0.2,1,2,100,3,4)=...`    |
| `convolve(a,b)`           | Discrete convolution                                    | `convolve({1,2},{3,4})={3,10,8}`  |

### Numerical Computation (Floating-Point Utilities)

| Function Name | Description                                  | Example Inputs/Outputs |
| ------------- | -------------------------------------------- | ---------------------- |
| `hypot(x,y)`  | √(x² + y²) (robust against cancellation)     | `hypot(3,4)=5`         |
| `fma(a,b,c)`  | Compute a\*b + c with a single rounding step | `fma(10,10,2)=102`     |
| `rms(...)`    | Root mean square                             | `rms(1,-1,1,-1)=1`     |

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

#### Geometry / Vectors

| Function Name                       | Description                           | Example                           |
| ----------------------------------- | ------------------------------------- | --------------------------------- |
| `norm(x,y)`                         | Vector length                         | `norm(3,4)=5`                     |
| `lerp(a,b,t) `                      | Linear interpolation                  | `lerp(0,10,0.5)=5`                |
| `distance(x1,y1,x2,y2)`             | Distance                              | `distance(0,0,3,4)=5`             |
| `vadd(a,b)`                         | Vector addition                       | `vadd({1,2},{3,4})={4,6}`         |
| `vsub(a,b)`                         | Vector Subtraction                    | `vsub({5,6},{1,2})={4,4}`         |
| `vscalar(a,s)`                      | Scalar Multiplication                 | `vscalar({1,2},3)={3,6}`          |
| `vdot(a,b)`                         | Dot product                           | `vdot({1,2},{3,4})=11`            |
| `vcross(a,b)`                       | Cross product (3D only)               | `vcross({1,0,0},{0,1,0})={0,0,1}` |
| `vnorm(a)`, `vlength(a)`            | Vector length                         | `vnorm({3,4})=5`                  |
| `vmanhattan(a,b)`                   | Manhattan distance                    | `vmanhattan({0,0},{3,4})=7`       |
| `veuclidean(a,b)`, `vdistance(a,b)` | Euclidean distance                    | `veuclidean({0,0},{3,4})=5`       |
| `vnormalize(a)`                     | Vector normalization (unit vector)    | `vnormalize({3,4})={0.6,0.8}`     |
| `vproject(a,b)`                     | Project vector a onto vector b        | `vproject({1,2},{0,1})={0,2}`     |
| `vangle(a,b)`                       | Angle between two vectors (degrees)   | `vangle({1,0},{0,1})=90`          |
| `vreflect(a,n)`                     | Reflect vector a with normal vector n | `vreflect({1,1},{0,1})={1,-1}`    |
| `vunit(a)`                          | Unit vector                           | `vunit({3,4})={0.6,0.8}`          |

### Matrix

- Prioritizing rigor over execution speed, these operations can become extremely slow with large matrices.

| Function Name               | Description                                                                         | Input and Output Examples                                                       |
| --------------------------- | ----------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| `matrix(a,b,c,...)`         | Create a matrix of arbitrary size                                                   | `matrix({1,2},{3,4})`                                                           |
| `identity(n)`               | Create an n×n identity matrix                                                       | `identity(3)`                                                                   |
| `zeros(rows,cols)`          | Create a zero matrix of size rows x cols                                            | `zeros(2,3)`                                                                    |
| `mget(A,row,col)`           | Access matrix elements                                                              | `mget({{1,2},{3,4}},0,1)=2`                                                     |
| `madd(A,B)`                 | Matrix addition                                                                     | `madd({{1,2},{3,4}}, {{5,6},{7,8}})`                                            |
| `mmul(A,B)`                 | Matrix multiplication                                                               | `mmul({{1,2},{3,4}}, {{5,6},{7,8}})`                                            |
| `mtranspose(A)`             | Matrix transpose                                                                    | `mtranspose({{1,2},{3,4}})`                                                     |
| `mtrace(A)`                 | Matrix trace (sum of diagonal entries)                                              | `mtrace({{1,2},{3,4}})=5`                                                       |
| `mrank(A)`                  | Matrix rank                                                                         | `mrank({{1,2},{3,4}})=2`                                                        |
| `mdet(A)`                   | Determinant of matrix                                                               | `mdet({{1,2},{3,4}})=-2`                                                        |
| `mrows(A)`                  | Get number of rows in matrix                                                        | `mrows({{1,2},{3,4}})=2`                                                        |
| `mcols(A)`                  | Get number of columns in matrix                                                     | `mcols({{1,2},{3,4}})=2`                                                        |
| `mdiag(A)`                  | Get the diagonal elements of the matrix                                             | `mdiag({{1,2},{3,4}})={1,4}`                                                    |
| `mlu(A)`                    | LU decomposition of a matrix (with partial pivoting). The return value is `{P,L,U}` | `mlu({{2,1},{4,3}})= {{1,0},{...}}`                                             |
| `meigenvals(A)`             | Matrix eigenvalues                                                                  | `meigenvals({{1,2},{3,4}})={5.372.., -0.372..}`                                 |
| `meigenvecs(A)`             | Matrix eigenvectors                                                                 | `meigenvecs({{1,2},{3,4}})={{0.415.., 0.909..}, {0.824.., -0.565..}}`           |
| `minverse(A)`               | Inverse matrix                                                                      | `minverse({{1,2},{3,4}})={{-2,1},{1.5,-0.5}}={{-2,1},{1.5,-0.5}}`               |
| `mqr(A)`                    | QR decomposition (Householder transformation)                                       | `mqr({{1,2},{3,4}})={{{0.316227766017, ..},  {.., ..}},  {{.., ..}, {.., ..}}}` |
| `msvd(A)`                   | Singular value decomposition (Jacobi method)                                        | `msvd({{1,2},{3,4}})`                                                           |
| `mcond(A)`, `mcondition(A)` | Condition number of matrix (LU decomposition)                                       | `mcond({{1,2},{3,4}})`                                                          |
| `mlsqr(A,b)`                | Least squares solution (QR method)                                                  | `mlsqr({{1,2},{3,4}},{1,2})`                                                    |
| `mnormf(A)`, `mnorm(A)`     | Matrix Frobenius norm                                                               | `mnormf({{1,2},{3,4}})`                                                         |
| `mmaxeigen(A)`              | Matrix maximum eigenvalue                                                           | `mmaxeigen({{1,2},{3,4}})`                                                      |
| `mmaxeigen_power(A)`        | Maximum eigenvalue of a matrix (power method)                                       | `mmaxeigen_power({{1,2},{3,4}})=5.37228132327`                                  |
| `msumsv(A)`                 | Sum of Singular Values of Matrix (QR Iterative SVD)                                 | `msumsv({{1,2},{3,4}})`                                                         |
| `misse_symmetric(A)`        | Determines matrix symmetry                                                          | `misse_symmetric({{1,2},{2,1}})`                                                |
| `mispd(A)`                  | Determines matrix positive definiteness (Cholesky decomposition)                    | `mispd({{2,1},{1,2}})`                                                          |
| `covmatrix(M)`              | Covariance matrix                                                                   | `covmatrix({{1,2},{3,4}})=...`                                                  |
| `corrmatrix(...)`           | Correlation matrix                                                                  | `corrmatrix({1,2,3},{2,4,6})=...`                                               |
| `percentrank(x,...)`        | Percentile rank                                                                     | `percentrank(3,1,2,3,4,5)=0.5`                                                  |

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

### Area and Volume

| Function Name                         | Description                                  | Input and Output Example          |
| ------------------------------------- | -------------------------------------------- | --------------------------------- |
| `area_circle(d)`                      | Area of a circle (diameter specified)        | `area_circle(2)=3.14159...`       |
| `area_triangle(b,h)`                  | Triangle area (base and height)              | `area_triangle(3,4)=6`            |
| `area_triangle(a,b,c)`                | Triangle area (three side lengths)           | `area_triangle(3,4,5)=6`          |
| `area_trapezoid(a,b,h)`               | Trapezoid area                               | `area_trapezoid(3,5,4)=16`        |
| `area_polygon(x0,y0,x1,y1,...,xn,yn)` | Polygon area (coordinate specification)      | `area_polygon(0,0,1,0,1,1,0,1)=1` |
| `vol_cylinder(d,h)`                   | Cylinder volume (diameter specified)         | `vol_cylinder(2,3)=9.42477...`    |
| `vol_cone(d,h)`                       | Volume of a cone (diameter specified)        | `vol_cone(2,3)=3.14159...`        |
| `vol_sphere(d)`                       | Volume of a sphere (diameter specified)      | `vol_sphere(2)=4.18879...`        |
| `vol_prism(x0,y0,...,xn,yn,h)`        | Volume of prism with base polygon and height | `vol_prism(0,0,1,0,1,1,0,1,2)=2`  |

### Strength of Materials (Stress, Strain, Elasticity)

| Function Name      | Description            | Example                   |
| ------------------ | ---------------------- | ------------------------- |
| `stress(F,A)`      | Stress: σ = F/A        | `stress(1000,10)=100`     |
| `strain(dL,L)`     | Strain: ε = dL/L       | `strain(0.1,100)=0.001`   |
| `young(sigma,eps)` | Young’s modulus: E=σ/ε | `young(200,0.001)=200000` |

### Section Properties (I, Z, J)

| Function Name          | Description                            | Example                        |
| ---------------------- | -------------------------------------- | ------------------------------ |
| `moment_rect(b,h)`     | Rectangular area moment: I=b\*h³/12    | `moment_rect(10,20)=6666.6`    |
| `moment_circle(d)`     | Circular area moment: I=π\*d⁴/64       | `moment_circle(10)=490.87`     |
| `sectionmod_rect(b,h)` | Rectangular section modulus: Z=b\*h²/6 | `sectionmod_rect(10,20)=666.6` |
| `sectionmod_circle(d)` | Circular section modulus: Z=π\*d³/32   | `sectionmod_circle(10)=98.17`  |
| `torsion_J_circle(d)`  | Polar moment of inertia: J=π\*d⁴/32    | `torsion_J_circle(10)=981.74`  |
| `polarZ_circle(d)`     | Polar section modulus: Zp=π\*d³/16     | `polarZ_circle(10)=196.35`     |

### Mold & Injection Molding

| Function Name                                     | Description                                          | Example                                                              | Output Unit (MKS System) |
| ------------------------------------------------- | ---------------------------------------------------- | -------------------------------------------------------------------- | ------------------------ |
| `mold_clamp(P_avg, A_proj)`                       | Mold clamping force (clamp force)                    | `mold_clamp(5e6, 0.01)=50000`                                        | N (kg·m/s²)              |
| `mold_clamp_safe(SF, P_avg, A_proj)`              | Clamping force including safety factor               | `mold_clamp_safe(1.2, 5e6, 0.01)=60000`                              | N                        |
| `mold_Pinj(P_cav, ΔP_runner, ΔP_gate, ΔP_nozzle)` | Injection pressure                                   | `mold_Pinj(50e6, 1e6, 0.2e6, 0.1e6)=51300000`                        | Pa (kg/m·s²)             |
| `mold_flowrate(V, t_fill)`                        | Fill Flow Rate (Volume Flow Rate)                    | `mold_flowrate(0.0001, 0.5)=0.0002`                                  | m³/s                     |
| `mold_gate_velocity(Q, A_gate)`                   | Gate flow velocity                                   | `mold_gate_velocity(2e-4, 2e-5)=10`                                  | m/s                      |
| `mold_shear_gate(Q, b, h)`                        | Rectangular gate shear velocity                      | `mold_shear_gate(2e-4, 0.002, 0.001)=60000`                          | 1/s                      |
| `mold_shear_runner(Q, D)`                         | Shear velocity of runner circular pipe               | `mold_shear_runner(2e-4, 0.004)=31830.988...`                        | 1/                       |
| `mold_pressure_loss_runner(mu, L, Q, D)`          | Runner pressure loss (Newtonian fluid approximation) | `mold_pressure_loss_runner(500, 0.2, 2e-4, 0.004)=3183098861.837...` | Pa                       |
| `mold_eject_friction(mu, N)`                      | Mold ejection friction force                         | `mold_eject_friction(0.2, 1e4)=200`                                  | N                        |
| `mold_eject_contact(p_contact, A_contact)`        | Contact force                                        | `mold_eject_contact(5e6, 0.002)=10000`                               | N                        |
| `mold_eject_total(mu, p_contact, A_contact)`      | Total Ejection Force                                 | `mold_eject_total(0.2, 5e6, 0.002)=2000`                             | N                        |
| `mold_mu_tex(k, h)`                               | Texture friction correction factor                   | `mold_mu_tex(0.5, 0.001)=0.0005`                                     | Dimensionless            |
| `mold_pin_stress(F_eject, n, A_pin)`              | Ejector pin surface pressure                         | `mold_pin_stress(1000, 4, 1e-6)=250000000`                           | Pa                       |
| `mold_plate_deflection(K, P, a, E, t)`            | Maximum mold plate deflection                        | `mold_plate_deflection(0.005, 5e6, 0.2, 2e11, 0.01)=0.0002`          | m                        |

- This is an experimental implementation. In the future, self-evident and unnecessary functions may be removed, and the unit system may be subject to change.

### Fasteners / Friction (Bolts, Torque)

| Function Name                | Description                     | Example                                  |
| ---------------------------- | ------------------------------- | ---------------------------------------- |
| `bolt_stress(F,d)`           | Bolt stress: σ = F/(πd²/4)      | `bolt_stress(10000,10)=127.3`            |
| `torque_from_preload(F,d,K)` | Tightening torque: T = K F d    | `torque_from_preload(10000,0.01,0.2)=20` |
| `preload_from_torque(T,d,K)` | Inverse conversion: F = T/(K d) | `preload_from_torque(20,0.01,0.2)=10000` |
| `friction(mu,N)`             | Friction force: F = μ N         | `friction(0.2,100)=20`                   |

### Utilities

| Function Name              | Description  | Example Inputs/Outputs                                         |
| -------------------------- | ------------ | -------------------------------------------------------------- |
| `clamp(x, lo, hi)`         | Clamp        | `clamp(5,0,10)=5`                                              |
| `cnst("name")`             | Constant     | `cnst("celeritas")=299792458`                                  |
| `convert(x, "from", "to")` | Convert unit | `convert(10,"mm","cm")=1`<br>`convert(10,"mm")=0.01` (SI unit) |

Refer to the following for the list of constants-> [【List of Constants】](constants_list.md)

#### History

```text
In[n], Out[n], %, %%, %%%...
```

#### System Functions

System functions are not intended for use with other functions or formulas.

| Function Name               | Description                                                                               |
| --------------------------- | ----------------------------------------------------------------------------------------- |
| `Clear[]`                   | Clear erases history and sets `In[n] / Out[n]` to `n=1`                                   |
| `Exit[]`                    | Terminates the program                                                                    |
| `filein("filename")`        | Reads input from an external file and evaluates it. The result is stored in `Out[n]`      |
| `fileout(data, "filename")` | Outputs to a file. `[n]` does not increment.                                              |
| `clip(data)`                | Copies to the clipboard. `[n]` does not increment.                                        |
| `explain(data)`             | Displays information such as data type and internal values. `[n]` does not increment.     |
| `if(cond,a,b)`              | Conditional branching (experimental)                                                      |
| `silent(x)`                 | Utility for suppressing output (useful for suppressing the output of long matrices, etc.) |

Additionally, functions related to the calculator system are provided by a set of functions beginning with `:`. Some of these are compatible with the above.

| Function Name | Description                                            | Input         | Output Example                                                                                         |
| ------------- | ------------------------------------------------------ | ------------- | ------------------------------------------------------------------------------------------------------ |
| `:clear`      | Same as `Clear[]`                                      | `:clear`      | None                                                                                                   |
| `:exit`       | Same as `Exit[]`                                       | `:exit`       | None                                                                                                   |
| `:defs`       | Display a list of user-defined variables and functions | `:defs`       | `definitions => Variables:  x := 2   Functions:  f(x)`                                                 |
| `:unset`      | Delete user-defined variables                          | `:unset x`    | `<unset: x>`                                                                                           |
| `:undef`      | Delete user-defined functions                          | `:undef f`    | `<undefined: f>`                                                                                       |
| `:help`       | Function help                                          | `:help ibeta` | `ibeta[a, b, x]: regularized incomplete beta function I_x[a, b], where a > 0, b > 0, and 0 <= x <= 1.` |

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

### Comparison Operators

| Operator    | Meaning               |
| ----------- | --------------------- |
| `< <= > >=` | Relational comparison |
| `==`        | Equality comparison   |

- Results are **1 (true) / 0 (false)**

### Complex Number Behavior

- `sqrt(-1) = I`
- `exp(I)` follows Euler’s formula
- Real numbers are automatically promoted to complex numbers when required
- `+0I` and `-0I` both evaluate to zero. Therefore, `arg(-1+0I) = 180` and `arg(-1-0I) = 180` are obtained respectively.

### Error Handling

The following conditions are explicitly treated as **Errors**:

- Mismatched parentheses
- Invalid number of arguments
- Domain errors (`log(-1)`, `tan(90)`, etc.)
- Invalid syntax (`2**3`, `1+*2`, etc.)

#### Error Categories

| Error Type                                 | Example Cases                                                                                                                            | Comments                                                                                             |
| ------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| **Division by zero**                       | `1/0`, `0/0`, `I/0`, `mod(1,0)`                                                                                                          | Occurs when dividing by zero                                                                         |
| **Syntax error**                           | `1+*2`, `1/`, `()`, `(()`, `(`, `[)`                                                                                                     | Invalid expression syntax                                                                            |
| **Domain error**                           | `log(0)`, `(-1)!`, `1.5!`, `preload_from_torque(1,0,0.2)`,<br> `bolt_stress(1,0)`, `trimmean(0.5,1,2,3,4)`, `stddevs(1)`, `kurts(1,2,3)` | Calculation outside the function's domain (negative factorial, square root of negative number, etc.) |
| **Unknown identifier**                     | `foo`, `foo(1)`, `unknown(1)`, `cnst("")`, `cnst("GFN")`                                                                                 | Undefined function or constant                                                                       |
| **Result is infinite**                     | `tan(90)`, `cot(0)`, `sec(90)`, `csc(0)`, `atanh(1)`, `atanh(-1)`, `log(0+0I)`                                                           | Result becomes infinite                                                                              |
| **Result is NaN**                          | `0^I`, `(-0)^I`                                                                                                                          | Undefined numeric calculation                                                                        |
| **History out of range**                   | `Out[0]`, `Out[-100]`, `Out[100]`,                                                                                                       | History access is out of range                                                                       |
| **Multiple '.' in number**                 | `1..1`                                                                                                                                   | Syntax error in numeric notation                                                                     |
| **Unknown constant / unterminated string** | `cnst(")`, `cnst("GFN)`                                                                                                                  | Constant undefined or string not properly terminated                                                 |
| **Double or complex type error**           | `var(1,I)`, `median(1,I)`                                                                                                                | Calculation with unsupported type                                                                    |
| **Error in log with base**                 | `log(0,10)`, `log(1,0)`                                                                                                                  | Invalid logarithm base                                                                               |
| **Overflow detected**                      | `fib(-1)`                                                                                                                                | Value exceeds calculation range                                                                      |
| **Other domain-related**                   | `acos(2)`, `asin(2)`, `acosh(0.5)`, `atan2(0,0)`, `log1p(-1)`, `log1p(-2)`                                                               | Out of domain for trigonometric or logarithmic functions                                             |

#### Error Policy

- Evaluation stops immediately upon error
- Errors do not propagate as values
- IEEE754 results (`inf`, `-inf`) are distinguished from errors
- Distinction between returning `inf` and `result is infinite` is currently ambiguous (needs improvement)

#### Singularities and Exceptions

- When an error occurs, subsequent evaluations are not performed
- Errors are not propagated as values but are displayed immediately

### Key Points Express

- Implicit multiplication works (natural input)
- Use `%` or `Out[n]` to reuse past history in calculations
- Supports complex numbers
- Angles default to degrees but can be specified with `rad`, etc.
- Loads of functions and constants!
- Works anywhere!

# Design Notes

### Angle System

- All trigonometric functions use **degrees** by default

You can specify radians when explicitly stated as follows.

```text
sin(Pi/2 rad) = 1
```

### Numerical Precision

- Internally uses `double`
- Since the internal representation is double precision, it can only handle values within the range of 1.7E±308.
  Values exceeding ±1.7E+308 will result in an overflow or infinity error.

- Subject to IEEE754 rounding errors

```text
atan(0.57735026919) = 30.000000000016
```

Displayed values and history are rounded to 12 decimal places.

### Intentional Limitations

- No symbolic computation

### Known Behavior

- `0^-1` → `±inf` (IEEE754 compliant)
- `-0` may be displayed
- Expressions inside `In[] / Out[]` are undefined

### Specification behavior

- Spaces are generally ignored, but between numbers they are treated as implicit multiplication. `3 2` -> `3*2`
- In principle, complex solutions are returned as complex numbers, but `cbrt` returns real solutions. (`(-8)^(1/3)=1+1.732050807569I`, `cbrt(-8)=-2`)

## Grammar Definition (EBNF)

```ebnf
(* Lexical elements *)

digit         ::= "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9" ;
letter        ::= "A" | ... | "Z" | "a" | ... | "z" ;
underscore    ::= "_" ;
identifier    ::= ( letter | underscore ) , { letter | digit | underscore } ;

integer       ::= digit , { digit } ;
float         ::= digit , { digit } , "." , { digit }
                | "." , digit , { digit } ;
number        ::= integer | float ;

string_char   ::= ? any character except '"' ? ;
string        ::= '"' , { string_char } , '"' ;

constant      ::= "Pi" | "Tau" | "E" | "Phi" | "NA" | "I" | "ESP" ;
unit_name     ::= "deg" | "rad" | "grad" | "mm" | "cm" | "m" | "inch" ;

plus          ::= "+" ;
minus         ::= "-" ;
mul           ::= "*" ;
div           ::= "/" ;
pow           ::= "^" ;
bang          ::= "!" ;

assign        ::= ":=" ;

lt            ::= "<" ;
le            ::= "<=" ;
gt            ::= ">" ;
ge            ::= ">=" ;
eq            ::= "==" ;
neq           ::= "!=" ;

lparen        ::= "(" ;
rparen        ::= ")" ;
lbracket      ::= "[" ;
rbracket      ::= "]" ;
lbrace        ::= "{" ;
rbrace        ::= "}" ;
comma         ::= "," ;
percent       ::= "%" ;

(* Entry point *)

Input         ::= Assignment ;

(*  Assignment / definition right-associative *)

Assignment    ::= Compare
                | Assignable , assign , Assignment ;

Assignable    ::= identifier
                | FunctionSignature ;

FunctionSignature
              ::= identifier , CallOpen , [ ParamList ] , CallClose ;

ParamList      ::= identifier , { comma , identifier } ;

(* Comparison left-associative chain *)

Compare       ::= Expression ,
                  { CompareOp , Expression } ;

CompareOp     ::= lt | le | gt | ge | eq | neq ;

(* Arithmetic *)

Expression    ::= Term ,
                  { ( plus | minus ) , Term } ;

Term          ::= Unary ,
                  { ( mul | div | ImplicitMul ) , Unary } ;

Unary         ::= { plus | minus } , Power ;

Power         ::= Postfix ,
                  [ pow , Unary ] ;
                  (* parsePostfix "^" parseUnary *)

Postfix       ::= Primary , { bang } ;

(* Primary *)

Primary       ::= number
                | string
                | HistoryRef
                | MultiLiteral
                | IdentifierExpr
                | Group ;

HistoryRef    ::= percent , { percent } ;
                (* %, %%, %%% ... *)

MultiLiteral  ::= lbrace ,
                  [ Assignment , { comma , Assignment } ] ,
                  rbrace ;

IdentifierExpr
              ::= constant
                | identifier , [ CallSuffix ] ;

CallSuffix    ::= CallOpen ,
                  [ Assignment , { comma , Assignment } ] ,
                  CallClose ;

CallOpen      ::= lparen | lbracket ;
CallClose     ::= rparen | rbracket ;

Group         ::= lparen , Expression , rparen
                | lbracket , Expression , rbracket ;

(* Implicit multiplication *)

ImplicitMul   ::= /* no token: adjacency-based multiplication */ ;

UnitApplied   ::= ( number | constant ) , unit_name ;
```

## Command-Line Arguments

If a formula is passed as a command-line argument, the solution is output and the session ends.
The history function cannot be used.

Adding `--batch` to the command-line arguments enables a mode where `In[n]` and `Out[n]` are not displayed.
This is intended for external program processing.

## UI Version

This is a UI wrapper around the CLI version. Usage is the same as the CLI version.  
However, it requires the following runtime:

[.NET8.0 Runtime](https://dotnet.microsoft.com/ja-jp/download/dotnet/8.0)

## License

BSD 3-Clause License

Copyright (c) 2021–2026 mmKreutzef

If you use this software in academic papers or commercial products, please clearly state that fact in your documentation or publications.

A short notice would be appreciated.(and feel delight)

However, **simply repackaging** this DLL or source code with a different platform or UI for commercial use is technically permissible under the license, but it's a real letdown.

We'd be thrilled if you integrate it as part of some other software.

## Tests

This project includes **600 automated tests**, covering:

- Operator precedence and associativity
- Error handling
- Complex, Vector, Matrix arithmetic
- Boundary and edge cases
- Perform FFT->IFFT on a random array and measure maximum absolute error, mean absolute error, maximum relative error, and energy conservation.
  - 2048×2048 complex random: no observable round-trip error
  - 4096×4096 complex random: no observable round-trip error
  - 997×997 complex random: max absolute error ≈ 1.41e-12, energy relative error ≈ 1.76e-16

Refer to the test_set folder for test details.

## Notes

This project aims to reconcile rigorous operator semantics with practical formula evaluation.(But, due to incremental additions to the implementation, there is some overlap in responsibilities; I hereby pledge to streamline this.)

It also seeks to operate with ease in any environment, whether for designers or on the production floor.

And I'm just making what I want.

### **Disclaimer**

This software is provided “as is” without any warranty, express or implied, including but not limited to warranties of merchantability, fitness for a particular purpose, or non-infringement.
The author or copyright holder shall not be liable for any claims, damages, or other liabilities arising out of or in connection with the software, whether based on contract, tort, or otherwise.

By using this software, you acknowledge and automatically agree that all risks associated with its use are your responsibility.
The author assumes no responsibility for any damage resulting from data loss, system malfunction, or other consequences arising from the use of this software.

### Acknowledgments

I express my gratitude to my past self who inspired the development of this tool, to my university and professors who served as my place of learning, and to my current workplace as the environment where it is actually used.

## Requests / Contributions

Please submit them via Github.I love to hear any suggestions for implementation, especially those needed in manufacturing or design environments.

### Future Plans

- Support for different bases (_N_ notation, 0x, 0b, add(), shift(), >>)
- AngleMode support (setAngle = RADIAN)
- Arbitrary precision (setFix=5)
- `eval` now supports `expect`
- Basic calculus stuff -> Placeholder implementation
- Variables aren’t really supported on purpose, but loop-like stuff (like `for`) would be nice
- `plot` function for graphing
- BigInt? (still depends on double, so very large numbers can get a bit messy at the lower digits)
