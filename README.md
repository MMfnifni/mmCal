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

```text
abs, sign, sqrt, cbrt, exp, log
```

- `log(x)` : Natural logarithm
- `log(base, x)` : Logarithm with arbitrary base

---

#### Angle Conversion

```text
DtoR, DtoG, RtoD, RtoG
GtoD, GtoR
```

- D: Degree, R: Radian, G: Grad
- Converts the input value into the target angle unit
- `RtoD(Pi) = 180`
- `RtoG(4+10I) = 254.647908947033+636.619772367581I`

---

#### Trigonometric Functions (Degree-based)

```text
sin, cos, tan, cot, sec, csc
asin, acos, atan
```

- Both input and output are in **degrees**
- Other angle units can be used explicitly via `RtoD`, etc.
- Out-of-domain inputs result in an Error  
  (e.g. `tan(90)` → `Error: result is infinite`)

---

#### Hyperbolic Functions

```text
sinh, cosh, tanh
asinh, acosh, atanh
```

---

#### Number Theory & Combinatorics

```text
gcd, lcm
perm, comb
```

---

#### Aggregation & Control

```text
sum, mean, min, max
```

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
