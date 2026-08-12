# 不定積分・代表例カタログ

## 0. 記号

- **E** — 初等函数だけで表現可能
- **S** — 標準的な特殊函数が必要
- **H** — 楕円積分・超幾何函数などが必要
- **P** — 定義域・仮定・区分函数への注意が必要
- **U** — 一般的な有限個の標準函数では扱いにくく、未評価でも妥当

以下は原則として

```text
integrate[f[x], x]
```

形式を想定する。

実数範囲では `log[abs[x]]` を用いる。複素解析では局所的に `log[x]` として扱うなど、枝の選択が別途必要。

---

# 1. 定数・冪・多項式

| 種別 | 入力 | 結果 |
|---|---|---|
| E | `integrate[0, x]` | `0` |
| E | `integrate[1, x]` | `x` |
| E | `integrate[3, x]` | `3x` |
| E | `integrate[x, x]` | `x^2/2` |
| E | `integrate[x^2, x]` | `x^3/3` |
| E | `integrate[x^5, x]` | `x^6/6` |
| E | `integrate[x^(-2), x]` | `-1/x` |
| E | `integrate[x^(-3), x]` | `-1/(2x^2)` |
| E/P | `integrate[1/x, x]` | `log[abs[x]]` |
| E | `integrate[sqrt[x], x]` | `2x^(3/2)/3` |
| E | `integrate[1/sqrt[x], x]` | `2sqrt[x]` |

一般則：

```text
integrate[x^n, x]
    = x^(n+1)/(n+1)       n != -1
```

特殊点

```text
integrate[x^(-1), x]
    = log[abs[x]]
```

は別規則となる。

### 多項式

```text
integrate[x^3 + 2x^2 - 5x + 7, x]
= x^4/4 + 2x^3/3 - 5x^2/2 + 7x
```

```text
integrate[3x^5 - 4x^2 + 1, x]
= x^6/2 - 4x^3/3 + x
```

---

# 2. 線形函数の冪

```text
integrate[(2x + 1)^3, x]
= (2x + 1)^4/8
```

```text
integrate[(3x - 7)^(-2), x]
= -1/(3(3x - 7))
```

一般則：

```text
integrate[(a x + b)^n, x]
= (a x + b)^(n+1) / (a(n+1))
```

ただし

```text
a != 0
n != -1
```

の場合。

`n = -1` は

```text
integrate[1/(a x + b), x]
= log[abs[a x + b]]/a
```

---

# 3. 基本的な置換積分

これは積分器にとって極めて重要な規則群。

## 3.1 f'(x)/f(x)

```text
integrate[2x/(x^2 + 1), x]
= log[x^2 + 1]
```

```text
integrate[(3x^2 + 2)/(x^3 + 2x + 5), x]
= log[abs[x^3 + 2x + 5]]
```

一般形：

```text
integrate[D[f[x],x] / f[x], x]
= log[abs[f[x]]]
```

---

## 3.2 f'(x) f(x)^n

```text
integrate[2x (x^2 + 1)^5, x]
= (x^2 + 1)^6/6
```

```text
integrate[3x^2 sqrt[x^3 + 1], x]
= 2(x^3 + 1)^(3/2)/3
```

一般形：

```text
integrate[f'[x] f[x]^n, x]
= f[x]^(n+1)/(n+1)
```

---

## 3.3 exp[f(x)]

```text
integrate[2x exp[x^2], x]
= exp[x^2]
```

```text
integrate[(3x^2 + 1) exp[x^3 + x], x]
= exp[x^3 + x]
```

---

## 3.4 三角函数との合成

```text
integrate[2x cos[x^2], x]
= sin[x^2]
```

```text
integrate[3x^2 sin[x^3], x]
= -cos[x^3]
```

```text
integrate[2x sec[x^2]^2, x]
= tan[x^2]
```

---

# 4. 有理函数

有理函数

\[
\frac{P(x)}{Q(x)}
\]

は初等函数で必ず積分可能であり、CASの基本中の基本。

## 4.1 単純なもの

```text
integrate[1/(x + 1), x]
= log[abs[x + 1]]
```

```text
integrate[x/(x^2 + 1), x]
= log[x^2 + 1]/2
```

```text
integrate[(2x + 3)/(x^2 + 3x + 5), x]
= log[x^2 + 3x + 5]
```

---

## 4.2 部分分数

```text
integrate[1/(x(x + 1)), x]
= log[abs[x]] - log[abs[x + 1]]
```

```text
integrate[1/((x - 1)(x + 1)), x]
= log[abs[(x - 1)/(x + 1)]]/2
```

```text
integrate[1/(x^2(x + 1)), x]
= -log[abs[x]] - 1/x + log[abs[x + 1]]
```

---

## 4.3 二次式

```text
integrate[1/(x^2 + 1), x]
= atan[x]
```

```text
integrate[1/(x^2 + a^2), x]
= atan[x/a]/a
```

`a > 0` と仮定。

```text
integrate[1/(x^2 + 2x + 2), x]
= atan[x + 1]
```

```text
integrate[1/(a^2 - x^2), x]
= log[abs[(a + x)/(a - x)]]/(2a)
```

```text
integrate[1/(x^2 - a^2), x]
= log[abs[(x - a)/(x + a)]]/(2a)
```

---

# 5. 根号

## 5.1 基本

```text
integrate[sqrt[x], x]
= 2x^(3/2)/3
```

```text
integrate[x sqrt[x^2 + 1], x]
= (x^2 + 1)^(3/2)/3
```

```text
integrate[x/sqrt[x^2 + 1], x]
= sqrt[x^2 + 1]
```

---

## 5.2 円型

```text
integrate[1/sqrt[1 - x^2], x]
= asin[x]
```

一般に

```text
integrate[1/sqrt[a^2 - x^2], x]
= asin[x/a]
```

`a > 0`。

```text
integrate[sqrt[a^2 - x^2], x]
= x sqrt[a^2 - x^2]/2
  + a^2 asin[x/a]/2
```

---

## 5.3 双曲型

```text
integrate[1/sqrt[x^2 + a^2], x]
= asinh[x/a]
```

あるいは

```text
= log[x + sqrt[x^2 + a^2]]
```

定数差を除いて同値。

```text
integrate[sqrt[x^2 + a^2], x]
= x sqrt[x^2 + a^2]/2
  + a^2 asinh[x/a]/2
```

```text
integrate[1/sqrt[x^2 - a^2], x]
= log[abs[x + sqrt[x^2 - a^2]]]
```

```text
integrate[sqrt[x^2 - a^2], x]
= x sqrt[x^2 - a^2]/2
  - a^2 log[abs[x + sqrt[x^2 - a^2]]]/2
```

---

## 5.4 根号置換

```text
integrate[1/(sqrt[x](1 + x)), x]
= 2atan[sqrt[x]]
```

```text
integrate[sqrt[x]/(1 + x), x]
= 2sqrt[x] - 2atan[sqrt[x]]
```

---

# 6. 指数函数

```text
integrate[exp[x], x]
= exp[x]
```

```text
integrate[exp[2x], x]
= exp[2x]/2
```

```text
integrate[exp[a x], x]
= exp[a x]/a
```

```text
integrate[2^x, x]
= 2^x/log[2]
```

一般に

```text
integrate[a^x, x]
= a^x/log[a]
```

---

## 6.1 多項式 × exp

```text
integrate[x exp[x], x]
= exp[x](x - 1)
```

```text
integrate[x^2 exp[x], x]
= exp[x](x^2 - 2x + 2)
```

```text
integrate[x^3 exp[x], x]
= exp[x](x^3 - 3x^2 + 6x - 6)
```

---

## 6.2 exp × 三角函数

```text
integrate[exp[a x] sin[b x], x]
= exp[a x](a sin[b x] - b cos[b x])/(a^2 + b^2)
```

```text
integrate[exp[a x] cos[b x], x]
= exp[a x](a cos[b x] + b sin[b x])/(a^2 + b^2)
```

特に

```text
integrate[exp[x] sin[x], x]
= exp[x](sin[x] - cos[x])/2
```

```text
integrate[exp[x] cos[x], x]
= exp[x](sin[x] + cos[x])/2
```

---

# 7. 対数

```text
integrate[log[x], x]
= x log[x] - x
```

```text
integrate[x log[x], x]
= x^2 log[x]/2 - x^2/4
```

```text
integrate[log[x]^2, x]
= x(log[x]^2 - 2log[x] + 2)
```

```text
integrate[log[x]/x, x]
= log[x]^2/2
```

```text
integrate[1/(x log[x]), x]
= log[abs[log[x]]]
```

より一般に

```text
integrate[log[x]^n/x, x]
= log[x]^(n+1)/(n+1)
```

---

# 8. 基本三角函数

## 8.1 sin / cos

```text
integrate[sin[x], x]
= -cos[x]
```

```text
integrate[cos[x], x]
= sin[x]
```

```text
integrate[sin[a x], x]
= -cos[a x]/a
```

```text
integrate[cos[a x], x]
= sin[a x]/a
```

---

# 9. tan / cot / sec / csc

```text
integrate[tan[x], x]
= -log[abs[cos[x]]]
```

```text
integrate[cot[x], x]
= log[abs[sin[x]]]
```

```text
integrate[sec[x], x]
= log[abs[sec[x] + tan[x]]]
```

```text
integrate[csc[x], x]
= log[abs[csc[x] - cot[x]]]
```

```text
integrate[sec[x]^2, x]
= tan[x]
```

```text
integrate[csc[x]^2, x]
= -cot[x]
```

```text
integrate[sec[x] tan[x], x]
= sec[x]
```

```text
integrate[csc[x] cot[x], x]
= -csc[x]
```

---

# 10. 目的の例：負冪三角函数

```text
integrate[sin[2x]^(-2), x]
```

は

```text
integrate[csc[2x]^2, x]
```

なので

```text
= -cot[2x]/2
```

つまりこれは特殊函数など不要で、かなり基本的な **E** クラス。

一般に

```text
integrate[sin[a x]^(-2), x]
= -cot[a x]/a
```

```text
integrate[cos[a x]^(-2), x]
= tan[a x]/a
```

---

# 11. 三角函数の整数冪

## 11.1 二乗

```text
integrate[sin[x]^2, x]
= x/2 - sin[2x]/4
```

```text
integrate[cos[x]^2, x]
= x/2 + sin[2x]/4
```

---

## 11.2 三乗

```text
integrate[sin[x]^3, x]
= -cos[x] + cos[x]^3/3
```

```text
integrate[cos[x]^3, x]
= sin[x] - sin[x]^3/3
```

---

## 11.3 四乗

```text
integrate[sin[x]^4, x]
= 3x/8 - sin[2x]/4 + sin[4x]/32
```

```text
integrate[cos[x]^4, x]
= 3x/8 + sin[2x]/4 + sin[4x]/32
```

---

## 11.4 混合

```text
integrate[sin[x]^2 cos[x]^2, x]
= x/8 - sin[4x]/32
```

```text
integrate[sin[x] cos[x], x]
= sin[x]^2/2
```

```text
integrate[sin[x]^3 cos[x], x]
= sin[x]^4/4
```

```text
integrate[sin[x] cos[x]^3, x]
= -cos[x]^4/4
```

---

# 12. tan / cot の冪

```text
integrate[tan[x]^2, x]
= tan[x] - x
```

```text
integrate[cot[x]^2, x]
= -cot[x] - x
```

```text
integrate[sec[x]^3, x]
= (sec[x] tan[x] + log[abs[sec[x] + tan[x]]])/2
```

```text
integrate[csc[x]^3, x]
= (-csc[x] cot[x] + log[abs[csc[x] - cot[x]]])/2
```

これらから reduction formula によって高冪へ拡張できる。

---

# 13. 高冪三角函数

例えば

```text
integrate[sin[2x]^64, x]
```

も原理的には完全に初等的。

偶数冪なので

\[
\sin^{64}(2x)
\]

を有限Fourier級数

\[
c_0+\sum_{k=1}^{32}c_k\cos(4kx)
\]

へ展開し、項別積分すればよい。

したがって

```text
integrate[sin[2x]^256, x]
```

も同様に有限個の項で厳密に求まる。

計算量・出力サイズは大きくなるが、**数学的困難はない**。

---

# 14. 積和公式を必要とするもの

```text
integrate[sin[a x] cos[b x], x]
```

`a != ±b` なら

```text
= -cos[(a+b)x]/(2(a+b))
  -cos[(a-b)x]/(2(a-b))
```

```text
integrate[sin[a x] sin[b x], x]
= sin[(a-b)x]/(2(a-b))
  - sin[(a+b)x]/(2(a+b))
```

```text
integrate[cos[a x] cos[b x], x]
= sin[(a-b)x]/(2(a-b))
  + sin[(a+b)x]/(2(a+b))
```

`a = b` などは別ケースに退化するので条件分岐が必要。

---

# 15. 三角有理式

## 15.1 簡単な恒等変形

```text
integrate[1/(1 + sin[x]), x]
= tan[x] - sec[x]
```

```text
integrate[1/(1 - sin[x]), x]
= tan[x] + sec[x]
```

```text
integrate[1/(1 + cos[x]), x]
= tan[x/2]
```

```text
integrate[1/(1 - cos[x]), x]
= -cot[x/2]
```

---

## 15.2 Weierstrass置換

一般に

```text
t = tan[x/2]
```

とすれば

```text
sin[x] = 2t/(1+t^2)
cos[x] = (1-t^2)/(1+t^2)
dx = 2dt/(1+t^2)
```

なので

```text
integrate[R[sin[x], cos[x]], x]
```

は有理函数積分へ変換できる。

例えば `a > abs[b]` なら

```text
integrate[1/(a + b cos[x]), x]
=
2/sqrt[a^2-b^2]
atan[
    sqrt[(a-b)/(a+b)] tan[x/2]
]
```

---

# 16. 部分積分

## 16.1 x × trig

```text
integrate[x sin[x], x]
= sin[x] - x cos[x]
```

```text
integrate[x cos[x], x]
= x sin[x] + cos[x]
```

```text
integrate[x^2 sin[x], x]
= -x^2 cos[x] + 2x sin[x] + 2cos[x]
```

```text
integrate[x^2 cos[x], x]
= x^2 sin[x] + 2x cos[x] - 2sin[x]
```

---

## 16.2 x × inverse trig

```text
integrate[x atan[x], x]
=
(x^2 + 1)atan[x]/2 - x/2
```

---

# 17. 逆三角函数

```text
integrate[asin[x], x]
= x asin[x] + sqrt[1 - x^2]
```

```text
integrate[acos[x], x]
= x acos[x] - sqrt[1 - x^2]
```

```text
integrate[atan[x], x]
= x atan[x] - log[1 + x^2]/2
```

---

# 18. 双曲線函数

```text
integrate[sinh[x], x]
= cosh[x]
```

```text
integrate[cosh[x], x]
= sinh[x]
```

```text
integrate[tanh[x], x]
= log[cosh[x]]
```

```text
integrate[coth[x], x]
= log[abs[sinh[x]]]
```

```text
integrate[sech[x]^2, x]
= tanh[x]
```

```text
integrate[csch[x]^2, x]
= -coth[x]
```

```text
integrate[sech[x] tanh[x], x]
= -sech[x]
```

```text
integrate[csch[x] coth[x], x]
= -csch[x]
```

---

## 18.1 sech / csch

```text
integrate[sech[x], x]
= atan[sinh[x]]
```

同値な形として

```text
= 2atan[tanh[x/2]]
```

```text
integrate[csch[x], x]
= log[abs[tanh[x/2]]]
```

---

# 19. 逆双曲線函数

```text
integrate[asinh[x], x]
= x asinh[x] - sqrt[x^2 + 1]
```

```text
integrate[acosh[x], x]
= x acosh[x] - sqrt[x^2 - 1]
```

`x > 1` を想定。

```text
integrate[atanh[x], x]
= x atanh[x] + log[1 - x^2]/2
```

適切な実数領域内での式。

---

# 20. Gaussian ― erf / erfi

ここから「初等函数では積分不能」になる。

## 20.1 erf

**S**

```text
integrate[exp[-x^2], x]
= sqrt[Pi]/2 erf[x]
```

一般に `a > 0` なら

```text
integrate[exp[-a x^2], x]
= sqrt[Pi]/(2sqrt[a]) erf[sqrt[a] x]
```

---

## 20.2 erfi

```text
integrate[exp[x^2], x]
= sqrt[Pi]/2 erfi[x]
```

したがって `erf` だけでなく **erfi** も積分システムにはかなり有用。

---

# 21. 指数積分 Ei

**S**

```text
integrate[exp[x]/x, x]
= Ei[x]
```

```text
integrate[exp[2x]/x, x]
= Ei[2x]
```

```text
integrate[exp[-x]/x, x]
= Ei[-x]
```

さらに面白い例：

```text
integrate[exp[exp[x]], x]
= Ei[exp[x]]
```

---

# 22. 対数積分 li

**S**

```text
integrate[1/log[x], x]
= li[x]
```

かつ

```text
li[x] = Ei[log[x]]
```

なので、独立函数を持たず `Ei` に還元する設計も可能。

---

# 23. Si / Ci

## 正弦積分

```text
integrate[sin[x]/x, x]
= Si[x]
```

## 余弦積分

```text
integrate[cos[x]/x, x]
= Ci[x]
```

これらは非常に代表的な「初等積分不能」の例。

---

# 24. Shi / Chi

双曲線版。

```text
integrate[sinh[x]/x, x]
= Shi[x]
```

```text
integrate[cosh[x]/x, x]
= Chi[x]
```

したがって特殊函数体系としては

```text
Si
Ci
Shi
Chi
Ei
```

はかなり自然な一群になる。

---

# 25. Fresnel積分

標準定義を

\[
C(z)=\int_0^z \cos\left(\frac{\pi t^2}{2}\right)dt
\]

\[
S(z)=\int_0^z \sin\left(\frac{\pi t^2}{2}\right)dt
\]

とする。

すると

```text
integrate[cos[x^2], x]
=
sqrt[Pi/2] fresnelc[sqrt[2/Pi] x]
```

```text
integrate[sin[x^2], x]
=
sqrt[Pi/2] fresnels[sqrt[2/Pi] x]
```

これは `fresnelc/fresnels` を入れるなら必須級のテスト。

さらに

```text
integrate[cos[Pi x^2/2], x]
= fresnelc[x]
```

```text
integrate[sin[Pi x^2/2], x]
= fresnels[x]
```

は定義そのものなので最初に通したい。

---

# 26. Dilogarithm / Polylogarithm

**S**

```text
integrate[log[1 - x]/x, x]
= -polylog[2, x]
```

```text
integrate[log[1 + x]/x, x]
= -polylog[2, -x]
```

```text
integrate[log[1 + x^2]/x, x]
= -polylog[2, -x^2]/2
```

一般則：

```text
integrate[polylog[s, x]/x, x]
= polylog[s + 1, x]
```

この規則は非常に美しい。

---

# 27. 不完全Gamma函数

**S/H**

```text
integrate[x^(s - 1) exp[-x], x]
= lowerGamma[s, x]
```

あるいは上側不完全Gammaを用いるなら

```text
= -upperGamma[s, x]
```

一般化して

```text
integrate[exp[-x^n], x]
=
lowerGamma[1/n, x^n]/n
```

適切な枝・領域を仮定。

例えば

```text
integrate[exp[-x^4], x]
= lowerGamma[1/4, x^4]/4
```

Gaussian

```text
exp[-x^2]
```

も本質的にはこの系列に属する。

---

# 28. 楕円積分 EllipticF

ここから通常の初等函数・erf等より一段上。

## 第一種

```text
integrate[
    1/sqrt[1 - m sin[x]^2],
    x
]
=
ellipticF[x, m]
```

---

## 第二種

```text
integrate[
    sqrt[1 - m sin[x]^2],
    x
]
=
ellipticE[x, m]
```

---

## 第三種

```text
integrate[
    1 / (
        (1 - n sin[x]^2)
        sqrt[1 - m sin[x]^2]
    ),
    x
]
=
ellipticPi[n, x, m]
```

---

# 29. 代数函数から楕円積分へ

非常に重要な境界例。

```text
integrate[
    1/sqrt[1 - x^4],
    x
]
=
ellipticF[asin[x], -1]
```

なぜなら

\[
1-x^4=(1-x^2)(1+x^2)
\]

だからである。

より一般に

```text
integrate[
  1/sqrt[(1-x^2)(1-k^2 x^2)],
  x
]
=
ellipticF[asin[x], k^2]
```

また

```text
integrate[
  sqrt[(1-k^2 x^2)/(1-x^2)],
  x
]
=
ellipticE[asin[x], k^2]
```

---

# 30. Hypergeometric 2F1

さらに一般化すると、多数の積分を `hypergeometric2F1` 一つで表せる。

代表公式：

\[
\int x^m(1+\beta x^n)^p dx
\]

は

```text
x^(m+1)/(m+1)
*
hypergeometric2F1[
    -p,
    (m+1)/n,
    1 + (m+1)/n,
    -beta x^n
]
```

となる。

したがって例えば

```text
integrate[1/(1 + x^n), x]
=
x hypergeometric2F1[
    1,
    1/n,
    1 + 1/n,
    -x^n
]
```

一般の記号 `n` に対して有用。

整数 `n` が固定なら、因数分解して初等函数へ落とせる場合も多い。

---

# 31. Bessel函数を含む積分

特殊函数を入力にも許す場合。

Bessel函数には

\[
\frac{d}{dx}\left(x^\nu J_\nu(x)\right)
=
x^\nu J_{\nu-1}(x)
\]

があるので

```text
integrate[
    x^v besselJ[v - 1, x],
    x
]
=
x^v besselJ[v, x]
```

特例：

```text
integrate[x besselJ[0, x], x]
= x besselJ[1, x]
```

また

```text
integrate[besselJ[1, x], x]
= -besselJ[0, x]
```

このあたりは将来的な特殊函数Knowledge向け。

---

# 32. 絶対値

実数領域限定。

```text
integrate[abs[x], x]
=
x abs[x]/2
```

これは `x = 0` でも微分可能なので、綺麗な大域的原始函数になる。

一方

```text
integrate[sign[x], x]
```

に対して単純に

```text
abs[x]
```

と返すのは厳密には問題がある。

`abs[x]` は `x=0` で微分不能だから

```text
D[abs[x],x] = sign[x]
```

は全実数上で成立しない。

**不連続函数は微分のDarboux性により、大域的な原始函数そのものが存在しない場合がある。**

これは厳密CASでは重要なテスト。

---

# 33. Piecewiseが必要になる積分

例えば

```text
integrate[1/abs[x], x]
```

は `x != 0` の各区間で

```text
sign[x] log[abs[x]]
```

と書けるが、0を跨いだ一個の滑らかな原始函数とは扱えない。

同様にパラメータ付き

```text
integrate[1/(x^2 + a^2), x]
```

も

- `a > 0`
- `a < 0`
- `a = 0`
- `a` が複素数

で適切な表現が変わる。

厳密CASではこの条件管理が非常に重要。

---

# 34. 「見た目は複雑だが簡単」なテスト

Knowledgeやsimplifyとの連携確認に有効。

```text
integrate[
    (1 + tan[x]^2),
    x
]
= tan[x]
```

```text
integrate[
    sin[x]^2 + cos[x]^2,
    x
]
= x
```

```text
integrate[
    (cos[x]^2 - sin[x]^2),
    x
]
= sin[2x]/2
```

```text
integrate[
    2sin[x]cos[x],
    x
]
= sin[x]^2
```

```text
integrate[
    exp[x](sin[x] + cos[x]),
    x
]
= exp[x] sin[x]
```

```text
integrate[
    (1 + x^2)^(-1),
    x
]
= atan[x]
```

---

# 35. 「置換を見抜けば一発」のテスト

```text
integrate[
    x/(1 + x^4),
    x
]
= atan[x^2]/2
```

```text
integrate[
    x^3/(1 + x^4),
    x
]
= log[1 + x^4]/4
```

```text
integrate[
    x/sqrt[1 - x^4],
    x
]
= asin[x^2]/2
```

```text
integrate[
    cos[x]/(1 + sin[x]),
    x
]
= log[abs[1 + sin[x]]]
```

```text
integrate[
    sin[x]/(1 + cos[x]),
    x
]
= -log[abs[1 + cos[x]]]
```

```text
integrate[
    exp[x]/(1 + exp[x]),
    x
]
= log[1 + exp[x]]
```

---

# 36. 初等積分できそうでできない代表例

以下はCAS積分器にとって重要な境界。

## Gaussian

```text
integrate[exp[-x^2], x]
```

初等函数では不能 → `erf`

## Fresnel

```text
integrate[sin[x^2], x]
```

初等函数では不能 → `fresnels`

```text
integrate[cos[x^2], x]
```

初等函数では不能 → `fresnelc`

## 指数積分

```text
integrate[exp[x]/x, x]
```

初等函数では不能 → `Ei`

## 正弦積分

```text
integrate[sin[x]/x, x]
```

初等函数では不能 → `Si`

## 対数積分

```text
integrate[1/log[x], x]
```

初等函数では不能 → `li`

## 楕円積分

```text
integrate[1/sqrt[1 - x^4], x]
```

初等函数では不能 → `ellipticF`

この6系列は「初等積分器から特殊函数積分器への境界」を試す非常に良いセット。

---

# 37. 有限個の普通の特殊函数へ落としにくい例

一般的なCASでも、どこまで特殊函数を許すかで結果が変わる。

```text
integrate[x^x, x]
```

**U**

通常は簡単な閉形式なし。

```text
integrate[sin[sin[x]], x]
```

**U**

Fourier/Bessel級数展開は可能だが、有限個の通常の初等函数では表現しにくい。

```text
integrate[exp[sin[x]], x]
```

**U**

Bessel級数などによる表現は可能だが、単純な有限閉形式ではない。

このようなものは

```text
integrate[...]
```

を未評価のまま返す判断も数学的に妥当。

---

# 38. 積分Knowledgeとして特に重要な一般規則

個々の積分表を何千件登録するより、この規則を持つ方が強い。

## 線形性

```text
integrate[a f[x] + b g[x], x]
=
a integrate[f[x],x]
+
b integrate[g[x],x]
```

---

## 定数因子

```text
integrate[a f[x], x]
=
a integrate[f[x], x]
```

`a` が `x` に依存しない場合。

---

## 合成函数

既知の

```text
integrate[f[u], u] = F[u]
```

に対して

```text
integrate[
    f[g[x]] g'[x],
    x
]
=
F[g[x]]
```

---

## 対数微分

```text
integrate[f'[x]/f[x], x]
=
log[abs[f[x]]]
```

---

## 指数微分

```text
integrate[f'[x] exp[f[x]], x]
=
exp[f[x]]
```

---

## sin / cos

```text
integrate[f'[x] cos[f[x]], x]
=
sin[f[x]]
```

```text
integrate[f'[x] sin[f[x]], x]
=
-cos[f[x]]
```

---

## tan

```text
integrate[f'[x] sec[f[x]]^2, x]
=
tan[f[x]]
```

---

## inverse trig

```text
integrate[
    f'[x]/sqrt[1 - f[x]^2],
    x
]
=
asin[f[x]]
```

```text
integrate[
    f'[x]/(1 + f[x]^2),
    x
]
=
atan[f[x]]
```

この「微分された内部函数を探す」機構は積分器の中核になる。

---

# 39. 回帰試験として特に良い最小セット

積分器を段階的に育てる場合は、まず次を全部通せるとかなり強い。

```text
integrate[x^5, x]
integrate[x^(-2), x]
integrate[1/x, x]

integrate[1/(2x + 1), x]
integrate[2x/(x^2 + 1), x]
integrate[1/(x^2 + 1), x]
integrate[1/(1 - x^2), x]

integrate[sqrt[x], x]
integrate[1/sqrt[1 - x^2], x]

integrate[exp[2x], x]
integrate[x exp[x], x]
integrate[exp[x] sin[x], x]

integrate[log[x], x]
integrate[log[x]/x, x]

integrate[sin[2x], x]
integrate[cos[2x], x]
integrate[tan[x], x]

integrate[sin[2x]^(-2), x]
integrate[cos[2x]^(-2), x]

integrate[sin[x]^2, x]
integrate[sin[x]^3, x]
integrate[sin[x]^4, x]
integrate[sin[x]^64, x]

integrate[sin[x] cos[x], x]
integrate[x sin[x], x]

integrate[sinh[x], x]
integrate[tanh[x], x]

integrate[asin[x], x]
integrate[atan[x], x]

integrate[exp[-x^2], x]
integrate[exp[x^2], x]

integrate[sin[x]/x, x]
integrate[cos[x]/x, x]

integrate[sin[x^2], x]
integrate[cos[x^2], x]

integrate[exp[x]/x, x]
integrate[1/log[x], x]

integrate[log[1-x]/x, x]

integrate[1/sqrt[1-x^4], x]

integrate[x^x, x]
```

最後の

```text
integrate[x^x, x]
```

まで「何ができるか」ではなく、

```text
これは既知の安全な閉形式へ変換できない
→ 未評価
```

と正しく判断できれば、積分器としてむしろ健全である。

---

# 40. 積分器全体の分類

実装体系として見るなら、大雑把には次の順序になる。

1. **線形性・定数抽出**
2. **多項式・冪**
3. **基本函数の直接積分表**
4. **線形引数 `f[a x+b]`**
5. **微分された内部函数の認識**
6. **有理函数・部分分数**
7. **三角恒等式**
8. **三角整数冪 reduction**
9. **Weierstrass置換**
10. **部分積分**
11. **根号・二次式**
12. **指数函数 × 多項式・三角函数**
13. **erf / erfi**
14. **Ei / Si / Ci / Shi / Chi**
15. **fresnelc / fresnels**
16. **gamma / incomplete gamma**
17. **polylog**
18. **ellipticF/E/Pi**
19. **hypergeometric**
20. **条件・枝・Piecewise管理**
21. **解けない積分を安全に未評価で返す**

この順序なら、単なる巨大な「積分表」ではなく、かなりCASらしい積分器へ発展できる。