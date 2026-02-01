# mm Calculator – Manufacturing-Oriented Mathematical Calculator

© 2021–2026 mmKreuzef (aka Daiki.NIIMI) Project mme
Licensed under the BSD 3-Clause License

[English](README.md) | [日本語](README.ja.md)

## Language Specification

### Overview

This project is a **mathematical expression evaluation engine** featuring **numerical computation, complex number arithmetic, and high-precision mathematical functions**.

While it can be used as a calculator, it is primarily designed for **CAD, mathematical processing, and practical use at manufacturing and design sites**, where both numerical accuracy and usability are required.

The design is inspired by _Mathematica_ in spirit, but **variable assignment and symbolic computation are intentionally not implemented**.

---

### Numeric Types

- **Real numbers**: IEEE 754 `double` (displayed with up to 12 decimal digits)
- **Complex numbers**: `a + bI`

---

### Constants

All constants start with an uppercase letter.

| Name  | Description                              |
| ----- | ---------------------------------------- |
| `Pi`  | π (pi)                                   |
| `Tue` | 2π (τ, tau)                              |
| `E`   | Base of natural logarithm                |
| `Phi` | Golden ratio                             |
| `NA`  | Avogadro constant                        |
| `I`   | Imaginary unit                           |
| `ESP` | Internal convergence tolerance (`1e-12`) |

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
```

**Implicit binding between function names and numbers is prohibited.**

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

#### Basic Mathematical Functions

| Function       | Description                                           | Input / Output Examples                      |
| -------------- | ----------------------------------------------------- | -------------------------------------------- |
| `abs(x)`       | Absolute value (supports complex numbers)             | `abs(-5) = 5`<br>`abs(3+4I) = 5`             |
| `sign(x)`      | Sign of a real number                                 | `sign(10)=1`<br>`sign(-3)=-1`<br>`sign(0)=0` |
| `sqrt(x)`      | Square root (returns complex for negative real input) | `sqrt(4)=2`<br>`sqrt(-9)=3I`                 |
| `cbrt(x)`      | Cube root                                             | `cbrt(-8)=-2`                                |
| `exp(x)`       | Exponential function                                  | `exp(1)=E`<br>`exp(I)=cos(1)+sin(1)I`        |
| `log(x)`       | Natural logarithm                                     | `log(E)=1`                                   |
| `log(base, x)` | Logarithm with arbitrary base                         | `log(2,8)=3`                                 |
| `log10(x)`     | Base-10 logarithm                                     | `log10(100)=2`                               |
| `log2(x)`      | Base-2 logarithm                                      | `log2(8)=3`                                  |
| `floor(x)`     | Floor function                                        | `floor(3.7)=3`                               |
| `ceil(x)`      | Ceiling function                                      | `ceil(3.2)=4`                                |
| `round(x)`     | Round to nearest integer                              | `round(1.5)=2`                               |
| `round(x, n)`  | Round with decimal precision                          | `round(1.234,2)=1.23`                        |
| `fract(x)`     | Fractional part                                       | `fract(1.25)=0.25`                           |
| `pow(x, y)`    | Power function                                        | `pow(2,3)=8`                                 |

#### Angle Conversion

| Function  | Description     | Example         |
| --------- | --------------- | --------------- |
| `DtoR(x)` | Degree → Radian | `DtoR(180)=Pi`  |
| `DtoG(x)` | Degree → Grad   | `DtoG(90)=100`  |
| `RtoD(x)` | Radian → Degree | `RtoD(Pi)=180`  |
| `RtoG(x)` | Radian → Grad   | `RtoG(Pi)=200`  |
| `GtoD(x)` | Grad → Degree   | `GtoD(200)=180` |
| `GtoR(x)` | Grad → Radian   | `GtoR(200)=Pi`  |

- D: Degree, R: Radian, G: Grad
- Converts the input value into the target angle unit
- `RtoD(Pi) = 180`
- `RtoG(4+10I) = 254.647908947033+636.619772367581I`

---

#### Trigonometric Functions (Degree-based)

| Function      | Description                      | Example         |
| ------------- | -------------------------------- | --------------- |
| `sin(x)`      | Sine (degree)                    | `sin(30)=0.5`   |
| `cos(x)`      | Cosine (degree)                  | `cos(60)=0.5`   |
| `tan(x)`      | Tangent (degree)                 | `tan(45)=1`     |
| `cot(x)`      | Cotangent                        | `cot(45)=1`     |
| `sec(x)`      | Secant                           | `sec(60)=2`     |
| `csc(x)`      | Cosecant                         | `csc(30)=2`     |
| `asin(x)`     | Inverse sine (returns degree)    | `asin(0.5)=30`  |
| `acos(x)`     | Inverse cosine (returns degree)  | `acos(0.5)=60`  |
| `atan(x)`     | Inverse tangent (returns degree) | `atan(1)=45`    |
| `atan2(y, x)` | Quadrant-aware inverse tangent   | `atan2(1,1)=45` |

- Both input and output are in **degrees**
- Other angle units can be used explicitly via `RtoD`, etc.
- Out-of-domain inputs result in an Error  
  (e.g. `tan(90)` → `Error: result is infinite`)

---

#### Hyperbolic Functions

| Function   | Description                | Example             |
| ---------- | -------------------------- | ------------------- |
| `sinh(x)`  | Hyperbolic sine            | `sinh(1)=1.1752`    |
| `cosh(x)`  | Hyperbolic cosine          | `cosh(1)=1.5431`    |
| `tanh(x)`  | Hyperbolic tangent         | `tanh(1)=0.7616`    |
| `asinh(x)` | Inverse hyperbolic sine    | `asinh(1)=0.8814`   |
| `acosh(x)` | Inverse hyperbolic cosine  | `acosh(2)=1.3170`   |
| `atanh(x)` | Inverse hyperbolic tangent | `atanh(0.5)=0.5493` |

---

#### Number Theory & Combinatorics

| Function    | Description             | Example        |
| ----------- | ----------------------- | -------------- |
| `gcd(a,b)`  | Greatest common divisor | `gcd(8,4)=4`   |
| `lcm(a,b)`  | Least common multiple   | `lcm(6,8)=24`  |
| `perm(n,r)` | Permutations            | `perm(5,2)=20` |
| `comb(n,r)` | Combinations            | `comb(5,2)=10` |

---

#### Aggregation & Control

| Function             | Description                         | Example                        |
| -------------------- | ----------------------------------- | ------------------------------ |
| `sum(...)`           | Sum of values                       | `sum(1,2,3)=6`                 |
| `mean(...)`          | Arithmetic mean                     | `mean(2,4)=3`                  |
| `min(...)`           | Minimum                             | `min(3,1,5)=1`                 |
| `max(...)`           | Maximum                             | `max(3,1,5)=5`                 |
| `var(x...)`          | Population deviation (divided by n) | `var(1,2,3)=0.666666666667`    |
| `vars(x...)`         | Unbiased variance (divided by n-1)  | `vars(1,2,3)=1`                |
| `stddev(x...)`       | standard deviation                  | `stddev(1,2,3)=0.816496580928` |
| `stddevs(x...)`      | Unbiased standard deviation         | `stddevs(1,2,3)=1`             |
| `median(...)`        | Median                              | `median(1,2,3,4)=2.5`          |
| `percentile(p, ...)` | Percentile (0–100)                  | `percentile(25,1,2,3,4)=1.75`  |
| `cov(x..., y...)`    | Covariance                          | `cov(1,2,3, 3,2,1)=-2/3`       |
| `corr(x..., y...)`   | Correlation coefficient             | `corr(1,2,3, 1,2,3)=1`         |

---

#### Geometry / Vector Operations

| Function                | Description               | Example               |
| ----------------------- | ------------------------- | --------------------- |
| `norm(x,y)`             | 2D vector length          | `norm(3,4)=5`         |
| `dot(x1,y1,x2,y2)`      | Dot product               | `dot(1,2,3,4)=11`     |
| `cross(x1,y1,x2,y2)`    | 2D cross product (scalar) | `cross(1,0,0,1)=1`    |
| `lerp(a,b,t)`           | Linear interpolation      | `lerp(0,10,0.5)=5`    |
| `distance(x1,y1,x2,y2)` | Distance between points   | `distance(0,0,3,4)=5` |

---

#### Control & Utility

| Function           | Description          | Example           |
| ------------------ | -------------------- | ----------------- |
| `clamp(x, lo, hi)` | Clamp value to range | `clamp(5,0,10)=5` |
| `rand()`           | Random number [0,1)  | `rand()`          |

---

#### History

```text
In[n], Out[n], %, %%, %%%...
```

---

#### System Functions

```text
Clear[], Exit[]
```

- `Clear[]` clears history and resets `In[n] / Out[n]` to 1
- `Exit[]` terminates the program

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

- エラー発生時，それ以降の評価は行われません
- エラーは値として伝播せず，直ちに表示されます
- IEEE754 に基づく結果（`inf`, `-inf`）はエラーと区別されます
- `inf`, `-inf`を返すか`result is infinite`は曖昧です(要改修です)

---

## Usage Examples

```text
In [1] := sin(30)^2 + cos(30)^2
Out[1] := 1

In [2] := sqrt(-4)
Out[2] := 2I
```

---

## Design Notes

### Angle System

- All trigonometric functions use **degrees** by default
- Degrees were chosen for usability in practical environments
- Use `RtoD` if needed
- Native radian support is a future extension

---

### Numerical Precision

- Internally uses `double`
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
A short notice would be appreciated.

---

## Tests

This project includes **300 automated tests**, covering:

- Operator precedence and associativity
- Error handling
- Complex arithmetic
- Boundary and edge cases

---

## Notes

This project aims to balance **strict operator semantics** with **practical expression evaluation**, ensuring it works reliably in design and manufacturing environments of any scale.

## Requests / Contributions

Please submit them via Git.
