# mmCal 仕様・函数リファレンス

この文書は **mmCalの実装そのもの** を基準にした詳細仕様書である。
ユーザー向けの導入はルートの`README.ja.md`を参照する。
この文書はバージョンごとに常に変動するため，過去バージョンはgitより引っ張り出してください。

## 1. 現在の設計思想

mmCalは、入力を最初から`double`へ落とす電卓ではなく **exact-firstの小型CAS / 数値計算kernel** とする。

優先順位は次の通り。

1. 数学的厳密性
2. principal branch / definedness / domainを失わないこと
3. exactとapproximateを混同しないこと
4. 共通知識をSolver・Simplifier・CertifiedEvaluatorで再利用すること
5. 上記を壊さない範囲で高速化すること

代表例:

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

`Pi`や`sqrt[2]`は「内部に保存した小数」ではない。exactな数式として保持し、`N[...]`が指定されたときだけcertified numerical evaluationへ進む。

---

# 2. 数値・式の内部モデル

## 2.1 Integer

任意長整数`BigInt`。

```text
123456789012345678901234567890
-> 123456789012345678901234567890

100!
-> exact BigInt
```

固定64bit整数へ丸めない。

## 2.2 Rational

有限小数も最初からexact Rationalとしてparseする。

```text
0.1
-> 1/10

1.25
-> 5/4

0.1 + 0.2 == 0.3
-> True
```

したがってIEEE 754由来の`0.30000000000000004`問題は通常評価には存在しない。

## 2.3 exact complex

実部・虚部ともexactな`Number`として保持する。

```text
I^2
-> -1

(3 + 4I) / 5
-> 3/5 + 4/5I
```

## 2.4 Symbolic expression

exactに閉じない式はASTのまま保持する。

```text
sqrt[2]
-> sqrt[2]

sin[1]
-> sin[1]

gamma[1/3]
-> gamma[1/3]
```

ここで`sin[1]`の1は **1 radian**。既定角度はRadianである。

## 2.5 certified approximation

任意精度の作業値は`BigFloat`、証明付き区間は`RealInterval` / `ComplexInterval`。

`N[expr,n]`では、真値を含む区間の両端が同じn桁丸めへ入ることを確認してから`DecimalApproximation`を返す。

現在の`DecimalApproximation`は表示文字列だけではなく、要求桁数、由来（exact入力 / certified interval）、表示10進値そのもののexact Rational、真値を含むexact Rational enclosureを保持する。`ComplexDecimalApproximation`も実部・虚部のmetadataを保持する。`precision/accuracy/rationalize`はこのmetadataを直接使い、表示文字列を再parseして精度を推測しない。

```text
N[sqrt[2],30]
-> 1.414213562373095048801688724210
```

「差が小さくなったので終了」という経験的停止条件だけには依存しない。

---

# 3. predefined symbols

現在の保護されたpredefined symbol:

| 名前 | 意味 |
|---|---|
| `Pi` | 円周率。exact transcendental constant |
| `E` | 自然対数の底。exact transcendental constant |
| `Phi` | 黄金比。exact algebraic constant |
| `I` | 虚数単位 |
| `True`, `False` | Boolean |
| `Integer` | 整数domain |
| `Rational` | 有理数domain |
| `Real` | 実数domain |
| `Complex` | 複素数domain |
| `Infinity` | exact値の`precision/accuracy`が返す無限精度sentinel。拡張実数算術の自動簡約はまだ限定的 |

旧版の`Tau`, `NA`, `ESP`は現在predefined constantではない。

---

# 4. 入力構文

## 4.1 function call

丸括弧と角括弧の両方を函数呼出に使用できる。

```text
sin(Pi/6)
sin[Pi/6]
```

グルーピングは丸括弧を使う。単独の`[x+1]`はgroupではない。

## 4.2 配列

```text
{1,2,3}
{{1,2},{3,4}}
```

内部ではshapeとflatten済みelementsを持つ`ArrayExpr`。

## 4.3 変数・ユーザー函数

```text
x := 3
-> 3

f(t) := t^2 + 1
f(4)
-> 17
```

変数への`:=`は右辺を評価して保存する`Set`相当。
函数定義は本体を保持する`SetDelayed`相当。

既存定義を別内容で上書きすると、評価結果とは別にInfo diagnosticを返す。

```text
x := 2
-> 2

x := 4
INFO: x redefined (was 2)
-> 4
```

函数の同一arity再定義でも同様に以前の定義を表示する。同じ値を再代入しただけならInfoは出さない。

組込み名・predefined symbolは再定義できない。

現在の定義確認・削除は次を使う。

```text
Defs[]
UnDef[x]
UnDef[x,y,f]
```

`Defs[]`は現在のglobal user variableとuser function definitionを式の配列として返す。`UnDef[...]`は指定した名前について変数定義と函数定義を削除し、実際に変更した名前数を返す。

## 4.4 履歴

```text
%
%%
%%%
```

成功した出力を直前から参照する。

絶対番号では次を使う。

```text
In[1]
Out[1]
```

`In[n]`は入力nのlowered Exprを取得した後、**現在のsession環境で通常評価する**。`Out[n]`は入力nが成功した現在の保存済み出力snapshotを返し、再評価しない。

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

したがって`In[n]`は「過去入力を現在環境へ貼り戻して再実行する」意味である。過去入力が変数参照・代入・乱数等を含めば現在の定義やRNG stateを使う。生の入力ASTを表示する用途とは分離する。現在評価中の入力自身を`In[n]`で参照することは禁止し、自己再帰によるstack overflowを防ぐ。

評価エラーになった入力でもparse/lowerまで成功していれば履歴slot自体は残るが、対応する`Out[n]`は存在しない。

## 4.5 比較

```text
<  <=  >  >=  ==  !=
```

確定できれば`True/False`、symbolicに未確定ならPredicate式を保持する。

## 4.6 演算子優先順位

実装上の重要点:

1. postfix `!`
2. power `^` — 右結合
3. unary `+ -`
4. `* /` と暗黙乗算 — 同じtermレベルで左から処理
5. `+ -`
6. comparison
7. assignment `:=`

```text
2^3^2
-> 512

-2^2
-> -4
```

## 4.7 暗黙乗算

```text
2Pi
2(x+1)
(x+1)(x-1)
2 x
2exp[x]
2E^x
```

函数名と数値の直接結合は函数呼出とは解釈せず、暗黙乗算として扱う。数値直後の`e/E`は、その後ろに指数の数字が実際に続く場合だけ科学表記へ取り込む。

## 4.8 基数付き数値literal

現在のLowererは次を扱う。

```text
0b1010
0o17
0xFF
2#1010
16#FF
```

`base#digits`形式ではRational literalもparse可能な範囲で扱う。
bit演算子自体はまだ未実装。

---

# 5. 角度仕様

**既定はRadian。**

```text
sin[Pi/6]
-> 1/2

asin[1/2]
-> Pi/6
```

明示単位:

```text
sin[30 Deg]
-> 1/2

sin[Pi/6 Rad]
-> 1/2

sin[100 Grad]
-> 1
```

単位接尾辞は数値literalだけでなく一般式へ適用できる。

```text
x Deg
Pi/6 Rad
(2x + 1) Grad
```

session既定角度は`angleMode`で確認・変更する。旧`:angle` commandは復活させない。

```text
angleMode[]
-> Rad

angleMode[Deg]
-> Deg

angleMode[Grad]
-> Grad
```

`angleMode`は表示設定ではなく数学的評価状態である。明示単位`Deg/Rad/Grad`は常にsession既定より優先する。Kernel APIの`setDefaultAngleUnit()`も同じ状態を操作する。

角度変換:

```text
DtoR[180] -> Pi
DtoG[90]  -> 100
RtoD[Pi]  -> 180
RtoG[Pi]  -> 200
GtoD[200] -> 180
GtoR[200] -> Pi
```

---

# 6. 基本演算・代数

通常の演算子:

```text
+  -  *  /  ^  !
```

`pow[x,y]`は`Power[x,y]`のsource alias、`fact[x]`はfactorial alias。

主なexact simplification:

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

ただしprincipal branchを壊す変形はしない。

```text
sqrt[x^2]
-> sqrt[x ^ 2]       // xが何者か不明

simplify[sqrt[x^2], element[x,Real]]
-> abs[x]

simplify[sqrt[x^2], x >= 0]
-> x
```

---

# 7. 基本数学・複素函数

| 函数 | 概要 | 例 |
|---|---|---|
| `sqrt[x]` | principal square root | `sqrt[-4] -> 2I` |
| `cbrt[x]` | 実立方根。実数domain | `cbrt[-8] -> -2` |
| `abs[z]` | 絶対値・複素magnitude | `abs[3+4I] -> 5` |
| `sign[z]` | 実符号 / 複素`z/abs[z]` | `sign[3+4I] -> 3/5+4/5I` |
| `re[z]` | 実部 | `re[3+4I] -> 3` |
| `im[z]` | 虚部 | `im[3+4I] -> 4` |
| `conj[z]` | 複素共役 | `conj[3+4I] -> 3-4I` |
| `arg[z]` | principal argument | `arg[-1] -> Pi Rad` |
| `hypot[x,y]` | exact `sqrt[x^2+y^2]` | `hypot[3,4] -> 5` |
| `cis[x]` | `cos[x]+I sin[x]` | `cis[Pi/3] -> 1/2 + I sqrt[3]/2` |
| `polar[r,t]` | `r cis[t]` | `polar[2,Pi/3]` |
| `nextpow2[x]` | 最小nで`2^n >= x` | `nextpow2[9] -> 4` |

互換alias:

```text
real -> re
imag -> im
mag  -> abs
unit,csgn -> sign
rect -> polar
```

---

# 8. 指数・対数

| 函数 | 仕様 |
|---|---|
| `exp[x]` | complex entire exponential |
| `log[x]` | principal natural logarithm |
| `log[b,x]` | principal `Log[x]/Log[b]` |
| `log2[x]` | `log[2,x]` frontend |
| `log10[x]` | `log[10,x]` frontend |
| `expm1[x]` | `exp[x]-1`を0近傍で安定評価 |
| `log1p[x]` | `log[1+x]`を0近傍で安定評価 |

例:

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
-> DomainError       // 底1
```

`log[0]`はInfinityへ置換せずDomainError。

---

# 9. 三角函数

実装済み:

```text
sin cos tan cot sec csc
asin acos atan atan2
```

既定Radなので:

```text
sin[Pi/6] -> 1/2
cos[Pi/3] -> 1/2
tan[Pi/4] -> 1
asin[1/2] -> Pi/6
atan2[1,-1] -> 3 Pi / 4
```

`tan`, `sec`, `cot`, `csc`のpoleはdefinednessとして扱い、有限値を捏造しない。

---

# 10. 双曲線函数

実装済み:

```text
sinh cosh tanh
asinh acosh atanh
csch sech coth
```

principal complex branchesを持つ逆双曲線函数は`MathRegistry`へbranch metadataを持つ。

---

# 11. cardinal / 安定初等函数

```text
sinc[x]
cosc[x]
tanc[x]
sinhc[x]
tanhc[x]
expc[x]
```

removable singularityはexactに埋める。

```text
sinc[0]  -> 1
cosc[0]  -> 0
tanc[0]  -> 1
sinhc[0] -> 1
tanhc[0] -> 1
expc[0]  -> 1
```

三角cardinal函数は角度表現をRadian量へ正規化してから比を取る。

```text
sinc[Pi/2]
-> 2 / Pi

sinc[90 Deg]
-> 2 / Pi
```

---

# 12. 丸め・整数utility

```text
floor ceil trunc round frac
gcd lcm mod rem quotient
```

例:

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

`mod`はfloor quotient、`rem`はtruncate-toward-zero quotientに対応する。

現在の`round`は1引数版のみ。旧`round[x,n]`はまだ戻していない。

---

# 13. 組合せ・軽量数論

```text
perm[n,r]
comb[n,r]
fib[n]
```

例:

```text
perm[10,3] -> 720
comb[10,3] -> 120
fib[100]   -> 354224848179261915075
```

`fib`はfast doubling O(log n)。
`comb`は対称性`r=min[r,n-r]`を利用。

`isprime/nextprime/prevprime/factorint/totient`は未実装一覧で追跡中。

---

# 14. 特殊函数

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

一般実数はStirling–Bernoulli + rigorous remainder、負実数はreflectionを使用。
非正整数poleはDomainError。

`lgamma[x]`は現在 **実軸上の`log[abs[gamma[x]]]`**。複素`LogGamma`とは分離している。

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

現在は正実数域に限定。

```text
beta[2,3] -> 1/12
beta[1/2,1/2] -> Pi
betaln[1/2,1/2] -> log[Pi]
```

一般Gamma比へ無条件展開してpole cancellationを壊さない。

## 14.4 generalized factorial系

```text
binom[x,n]
fallingfact[x,n]
risingfact[x,n]
```

現段階では非負整数次数へexact finite productを構成できる範囲が中心。

```text
binom[1/2,2] -> -1/8
fallingfact[5,3] -> 60
risingfact[5,3] -> 210
```

---

# 15. 集約函数

```text
sum
prod
min
max
mean
```

scalar variadicまたは1個のArrayを受ける。

```text
sum[1,2,3] -> 6
sum[{1,2,3}] -> 6
prod[] -> 1
sum[] -> 0
mean[1,2,4] -> 7/3
```

`min/max`は順序を証明できないsymbolic値を勝手に並べない。

```text
min[x,3]
-> min[x, 3]
```

`sum[f,{k,a,b}]`型の記号有限和はまだ未実装。

---

# 16. 記述統計

統計函数は原則 **exact real data** を受け、Rationalで閉じる量はRationalのまま返す。
1個のrank-1 Arrayまたはscalar列を受けるものが多い。

## 16.1 順序統計

```text
median
mode
quantile
percentile
iqr
percentrank
```

`quantile`はHyndman–Fan Type 7。

```text
median[1,2,3,4] -> 5/2
quantile[1/4,1,2,3,4,5,6,7] -> 5/2
iqr[1,2,3,4] -> 3/2
```

`mode`が複数ならArrayを返し、全値が1回ずつなら空Array。

## 16.2 分散・標準偏差

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

## 16.3 その他

```text
geomean harmmean rms
mad madR
skew kurtp kurts
cv stderr zscore
trimmean winsor winsorR
cov corr corrspearman
```

`cov/corr/corrspearman`の推奨形:

```text
cov[{1,2,3},{2,4,6}] -> 4/3
corr[{1,2,3},{2,4,6}] -> 1
```

互換用に偶数個scalarを前半/後半へ分割する形式も受ける。

---

# 17. Array / Vector / Matrix

## 17.1 Array utility

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

## 17.2 exact-first線形代数

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

symbolic rank/pivotの非零性を証明できない場合、勝手にpivotを選ばない。

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

Fourier位相はsessionの角度既定に依存せず、内部で明示Radian。

```text
dft[{1,2,3,4}]
-> {10,-2+2I,-2,-2-2I}

ifft[fft[{1+I,2-I,3+2I,4-3I}]]
-> {1+I,2-I,3+2I,4-3I}

convolve[{1,2},{3,4}]
-> {3,10,8}
```

2冪長FFTはradix-2 Cooley–Tukey。非2冪長はexact DFTへfallback。
現在はmachine FFTではない。

---

# 19. 記号微分と数値微分

## 19.1 `D`

```text
D[expr,x]
D[expr,{x,n}]
D[expr,x,y,...]
```

`{x,n}`は非負整数`n`階微分、複数specは左から順に適用する。

```text
D[sin[x],{x,4}]
-> sin[x]

D[x^2 y^3,x,y]
-> 6 x y^2
```

HoldAllなので既存の変数値に置換せず式を微分する。

```text
D[x^3 + 2x,x]
-> 2 + 3 x^2

D[exp[x^2],x]
-> 2 x exp[x^2]
```

既定Radでは:

```text
D[sin[x],x]
-> cos[x]
```

`abs/sign/re/im/conj/arg`など一般複素変数で通常のholomorphic derivativeを持たないものは、偽の微分を返さず未評価`D[...]`を保持する。

積分でも同じ微分知識を再利用する。

```text
D[integrate[f[x],x],x]
-> f[x]

D[integrate[t^2,{t,0,x}],x]
-> x^2
```

可変上下端や積分内部に微分変数が現れる一般形ではLeibniz ruleを形式的に構築する。分母が微分変数に依存しない商は、一般quotient ruleへ膨張させず`f'/c`を直接使う。

## 19.2 `diff`

```text
diff[expr,x,at]
diff[expr,x,at,digits]
```

有限差分の別公式ではなく、まず`D`でexact derivative Exprを作り、それを指定点でCertifiedEvaluatorへ渡す。

```text
diff[x^2,x,3]
-> 6.0000000000000000
```

---

# 20. 記号積分とcertified数値積分

## 20.1 `integrate` — exact/symbolic積分

```text
integrate[expr,x]
integrate[expr,{x,a,b}]
integrate[expr,x,assumptions]
integrate[expr,{x,a,b},assumptions]
```

第1形式は不定積分、第2形式はexact/symbolic定積分。積分変数はbinderとしてholdされ、同名のglobal定義に置換されない。

不定積分は「加法定数を法として選んだ原始函数の代表元」を返す。したがって`+ C`は表示しない。

```text
integrate[x^2,x]
-> x^3 / 3

integrate[(2x+3)^5,x]
-> (3 + 2 x)^6 / 12

integrate[1/(2x+3),x]
-> log[3 + 2 x] / 2
```

現在の主なexact規則:

- 定数、`x`、任意の有限多項式
- affine baseの有理冪。指数`-1`はLogへ送る
- Rational係数の有理函数。1次・2次分母に加え、Rational rootで1次因子へ分解でき、残余が高々既約2次となる場合はexact partial fractionへ分解する。重複1次因子にも対応
- 正のRational scaleを証明できる二次平方根型の`asin/asinh` primitive
- `sin^2/cos^2/tan^2/...`の安全な倍角・恒等式reduction
- 和・差・符号反転、積分変数に依存しない係数の線形性
- `exp/sin/cos/tan/cot/sec/csc`の安全な標準原始函数
- `sinh/cosh/tanh/coth/sech/csch`の安全な標準原始函数
- `log/log1p/expm1/sqrt/cbrt`
- `asin/acos/atan/asinh/acosh/atanh`
- `erf/erfc`
- exactな逆chain rule
- 多項式×`exp/sin/cos/sinh/cosh`に対する有限回のintegration by parts
- `exp[a x+b] sin/cos[c x+d]`型を連立一次式としてexact積分
- `sqrt[q(sqrt[x])]`で`q`がexact Rational係数2次式かつ正leadingの場合の` t=sqrt[x] `局所置換
- `x^n log[x]` (`n`が非負整数)

例:

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

角度単位も既存のAngleSemantics/Dと共有する。

```text
integrate[sin[x Deg],x]
-> -(180 / Pi cos[x Deg])
```

逆chain rule等で構造から候補原始函数を発見した場合は、既存の`D`をproof engineとして使い、候補の微分と元 integrand のexactな比例関係を証明してから採用する。単に数値点で一致した候補は採用しない。

branch/definednessを壊す**global simplification**は行わない。一方、原始函数は大域恒等式と同じ基準である必要はない。共通の解析領域上で正しい局所原始函数は、integrate専用の規則として採用できる。

例えばprincipal square rootについて

```text
integrate[1/sqrt[x^2-1],x]
-> log[sqrt[x^2-1] + x]
```

を返すが、Simplifierへ

```text
sqrt[x^2-1] == sqrt[x-1] sqrt[x+1]
```

という危険な大域規則は追加しない。積分公式と代数的恒等式の正当性を分離している。

また和の線形性は部分評価する。

```text
integrate[x^2 + gamma[x],x]
-> x^3 / 3 + integrate[gamma[x],x]
WARN: ... unevaluated integrate[...] remains
```

解けない部分だけを保持し、既に求まった項まで巻き戻さない。

### nested square-root substitution

現在は次のclassも扱う。

```text
integrate[sqrt[x + sqrt[x]],x]
-> 2 (x + sqrt[x]) sqrt[x + sqrt[x]] / 3
   - ((1 + 2 sqrt[x]) sqrt[x + sqrt[x]] / 4
      - log[1 + 2 sqrt[x + sqrt[x]] + 2 sqrt[x]] / 8)
```

これは単発公式ではなく、`t=sqrt[x]`により`2 t sqrt[q(t)]`へ落ちる、`q`がexact Rational係数2次式で正leadingのclassを処理する。principal branch上の局所置換則であり、一般のradical substitution探索ではない。

### exact/symbolic定積分

原始函数を安全に得られた場合、上下端へexact substitutionして差を取る。定義域条件がある函数は、可能ならCertifiedEvaluatorで**区間全体**をpreflightし、途中のpoleやbranch上の不成立を見逃さない。

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

特異点を跨ぐ場合はendpoint代入だけで値を作らない。

```text
integrate[1/x,{x,-1,1}]
-> WARN + unevaluated integrate[...]

integrate[tan[x],{x,0,2}]
-> WARN + unevaluated integrate[...]
```

Cauchy principal valueを自動的に意味することもない。

### assumptions と improper integral

第3引数にassumptionを渡せる。既存の`KnowledgeContext`へ統合され、`abs`や`sqrt[x^2]`のbranch-sensitive簡約に利用する。

```text
integrate[abs[x],x,x>=0]
-> x ^ 2 / 2

integrate[abs[x],x,x<=0]
-> -x ^ 2 / 2

integrate[sqrt[x^2],x,x>=0]
-> x ^ 2 / 2
```

endpointがInfinityまたは通常代入で定義されない場合は、原始函数の対応する片側/無限遠極限を使ってimproper integralを評価する。内部特異点がないことを証明できるclassだけを受理し、Cauchy principal valueは推測しない。

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

途中にpoleを含む可能性を排除できない場合は未評価保持する。

```text
integrate[1/(x-2),{x,1,Infinity}]
-> WARN + unevaluated integrate[...]
```

## 20.2 `limit` — exact/symbolic limit

```text
limit[expr,x,a]
limit[expr,x,a,-1]   // 左極限
limit[expr,x,a,1]    // 右極限
```

第4引数は方向を表し、`-1`が左、`1`が右。省略時は二側極限。方向は単なる表示指定ではなく一時的なassumption `x<a` / `x>a` としてKnowledgeContextへ渡される。

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

0/0型では既存`D`を使った反復l'Hopitalを安全弁付きで利用する。Rational functionの局所zero/pole次数や無限遠次数比較はexactに処理する。未解決二側極限をprincipal value等へ潰さない。

```text
limit[1/x,x,0]
-> WARN + limit[1 / x, x, 0]
```

## 20.3 `nintegrate` — certified数値積分

```text
nintegrate[expr,{x,a,b}]
nintegrate[expr,{x,a,b},digits]
```

```text
nintegrate[x^2,{x,0,1},12]
-> 0.333333333333
```

内部では区間を`x=a+(b-a)t`へ正規化し、exact Rational格子上のNewton–Cotesと高階導函数のcertified boundから求積誤差を包含する。

単なるdouble Simpsonの「近そうな値」ではなく、最終丸めが一意になった場合だけ返す。

高階導函数を構築する前に元の被積分函数を区間全体でpreflightするため、明白な特異点を早期に拒否する。特異点を跨いで偶然相殺する処理は行わない。

---

# 21. 式変形

```text
simplify[expr]
simplify[expr,assumptions]
fullSimplify[expr]
fullSimplify[expr,assumptions]
expand[expr]
factor[expr]
collect[expr,x]
```

例:

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

`fullSimplify`はbounded candidate search。短い式を選ぶために定義域を変えてよいわけではない。

```text
fullSimplify[(x^2-1)/(x-1)]
-> 元のholeを保持

fullSimplify[(x^2-1)/(x-1), x != 1]
-> 1 + x
```

---

# 22. Assumption / domain

```text
element[x,Real]
element[x,Integer]
```

`simplify/fullSimplify`の第2引数にはPredicate、Array、`And`相当の条件を渡せる。

```text
simplify[sqrt[x^2], element[x,Real]]
-> abs[x]

simplify[abs[x], x >= 0]
-> x
```

矛盾したassumptionはDomainError。

---

# 23. Solver

```text
solve[equation,x]
solve[equation,x,domainOrConstraint]
solve[{equations...},{variables...}]
```

既定ambient domainは等式系でComplex。
ordered inequalityはRealまたはそのsubdomainで扱う。

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

`SolutionSet`はEmpty / Finite / Universal / Conditional / Unresolvedを区別する。
対応外の式を「解なし」と誤認しない。

分母zero、Logのdefinedness、rational-functionのhole/pole等を可能な範囲でglobal conditionとして保持する。

Gammaのpole集合のように現Predicateで完全表現できない条件は、不完全な条件を捏造せずunresolvedのまま扱う。

### 実軸global inverse solve

MathRegistryはprincipal inverseだけでなく、実軸上のglobal injectivity・単調性・実値域をmetadataとして持つ。`solve`は**Realまたはそのsubdomain**で、globalに一対一であることを証明できる函数だけを安全に反転する。

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

`sin/cos/tan`にもinverse・range・period metadataは登録しているが、実軸全体ではglobal injectiveではない。整数parameterを含む全解family表現が未実装なので、principal inverse一個へ潰さない。

```text
solve[sin[x]==0,x,Real]
-> WARN + UnresolvedSolutionSet[x]
```

Complex領域でもprincipal inverseだけから全解を捏造しない。

---

# 24. `N` — numerical approximation

```text
N[expr]
N[expr,digits]
```

既定は16 fractional digits。
Arrayへ再帰的に適用できるほか、`arg`などが返す明示角度単位では値の部分だけを近似し、単位は保持する。

```text
N[Pi,20]
-> 3.14159265358979323846
N[Phi,20]
-> 1.61803398874989484820
N[fft[{1,2,3,4}],20]
N[arg[-1],20]
-> 3.14159265358979323846 Rad
```

exact Rationalが有限10進になる場合、表示は必要以上に0埋めしない。ただし現在は`requestedFractionalDigits`を別metadataとして保持する。

近似値は現在まだ一般の四則演算用Machine/ApproximateReal domainではない。表示値・要求桁・certified enclosureを別々に保持する。

---

# 25. precision / accuracy / rationalize

## 25.1 `accuracy[x]`

`DecimalApproximation`について、真値に対する**保証可能な絶対10進桁数の整数下限**を返す。

```text
accuracy[N[1/3,20]]
-> 20
```

表示値`d`とcertified source enclosure `[l,u]`から

```text
max(|d-l|, |d-u|)
```

を求め、さらに`N[...,n]`が近似値型であるという意味を失わないよう`0.5*10^-n`をsemantic error floorとして加味する。したがって、有限小数が真値と偶然完全一致しても要求桁を越えて無限accuracyとはしない。

exactな数・exact symbolic expressionは`Infinity`を返す。

```text
accuracy[1/3] -> Infinity
accuracy[Pi]  -> Infinity
```

## 25.2 `precision[x]`

同じabsolute error boundを、certified enclosureから得られる真値絶対値の下限で割り、**保証可能な相対10進桁数の整数下限**を返す。

```text
precision[N[1/3,20]]
-> 19
```

これは要求した20桁を機械的に返す函数ではない。`1/3`近傍では`0.5*10^-20`の絶対誤差floorが相対的には約`1.5*10^-20`となるため、厳密に保証できる整数桁数は19になる。0を含むenclosureでは相対誤差を正に下から評価できないため0を返す。exact expressionは`Infinity`。

## 25.3 `rationalize[x]`

近似値が持つcertified enclosure内から、**分母が最小になるexact Rational**を求める。探索はexact Rational上のcontinued-fraction型interval recursionで行い、doubleへ変換しない。

```text
rationalize[N[1/3,20]]
-> 1/3
```

`rationalize[x,tol]`は表示値を中心とする`[x-tol,x+tol]`内から最小分母Rationalを選ぶ。`tol`は非負exact real。

```text
rationalize[N[Pi,20],1/1000]
-> 201/64
```

`201/64`は`Pi`から0.001以内で、`355/113`より小さい分母を持つため、この仕様ではこちらが正しい。

`tol=0`は表示された有限10進値そのものをexact Rationalへ戻す。

```text
rationalize[N[1/3,20],0]
-> 33333333333333333333/100000000000000000000
```

mmCalではソースの`0.1`自体が最初から`1/10`なので、`rationalize[0.1]`は単に`1/10`のままである。Arrayや式内部のDecimalApproximationも再帰的にRational化する。

---

# 26. 乱数

乱数だけはstateful builtin。
RNG stateはKernelSessionごとに独立する。

## 25.1 seed

```text
randSeed[42]
-> 42
```

同じseedなら同じ列へ戻る。

```text
randSeed[42]
a := rand[]
randSeed[42]
rand[] == a
-> True
```

`randSeed[]`はentropyから再seedし、その再現用整数seedを返す。

## 25.2 uniform real

```text
rand[]
rand[hi]
rand[lo,hi]
```

`rand[]`は53bit dyadic lattice上のexact Rational。
例えばseed 42の先頭sampleは:

```text
randSeed[42]
rand[]
-> 227930101193189/1125899906842624
```

範囲:

```text
rand[]       : [0,1)
rand[hi]     : [0,hi), hi >= 0
rand[lo,hi]  : [lo,hi), lo <= hi
```

## 25.3 integer

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

BigInt範囲に対応。modulo biasを避けるrejection sampling。

## 25.4 choice

```text
choice[2,3,5,7]
choice[{2,3,5,7}]
```

rank-1 arrayまたはvariadic listから1要素。

## 25.5 normal

```text
randn[]
randn[mu]
randn[mu,sigma]
```

Box-Mullerをexact dyadic uniform sampleへ適用したsymbolic expressionを返す。

```text
randn[5,0]
-> 5

N[randn[],8]
-> 例: -0.80379286
```

`randn`の角度は明示Radで、sessionの角度既定に依存しない。

**暗号用途ではない。**

---

# 27. 条件分岐

```text
if[condition,trueExpr,falseExpr]
```

特殊形式として条件を先に評価し、選択されたbranchだけを評価する。
したがって非選択branch内のDomainErrorや乱数消費は発生しない。

---

# 28. 主要alias

| alias | canonical |
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

aliasは別実装ではなく同一`BuiltinId`へ束ねる。数学metadataやSolver規則を二重管理しない。

mmCal 1.5.0では、Mathematica互換だけを目的とした大文字始まりalias（`Sin`, `ArcTan`, `Integrate`, `Solve`等）を削除した。数学函数はlowercase canonicalを原則とする。`D`, `N`, `In`, `Out`, `Exit`, `Clear`, `Defs`, `UnDef`は記号演算・Kernel操作の固有名として例外的に維持する。互換構文が必要になった場合は、default namespaceへaliasを増やすのではなく独立したimport/compatibility層として検討する。

---

# 29. 現在のsource-callable函数一覧

mmCal 1.5.0では **205 builtin definitions / 187 source-callable names**。内部headはsource-callable数に含めない。

```text
Clear, D, Defs, DtoG, DtoR, Exit, GtoD, GtoR, In, N,
Out, RtoD, RtoG, UnDef, abs, accuracy, acos, acosh, angleMode, arg,
asin, asinh, atan, atan2, atanh, ave, beta, betaln, binom, cbrt,
ceil, choice, cis, collect, cols, comb, conj, convolve, corr, corrspearman,
cos, cosc, cosh, cot, coth, cov, csc, csch, csgn, cv,
det, dft, diag, diff, element, erf, erfc, exp, expand, expc,
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

# 30. エラー / Warning方針

主なCalcError種別:

- Syntax
- Domain
- Type
- Overflow
- Name
- Evaluation
- Internal

評価自体は成功したがalgorithmic builtinが処理を完了できない場合は、結果Exprと別にWarningを返す。`D`, `solve`, `N`, `rref`, `rank`, `integrate`に加え、`precision/accuracy/rationalize`が対象外入力を未評価保持する場合もWarningになる。

正常な状態変更の補足にはInfo diagnosticを使う。現在は変数・函数の再定義通知が対象。

```text
D[abs[x],x]
WARN: D could not fully evaluate the derivative; unevaluated D[...] remains
-> D[abs[x], x]
```

`sin[x]`のように記号函数として保持すること自体が正しい場合はWarningにしない。

数学的に未定義な値をNaN/Infへ流して継続するのではなく、原則その時点でError。

例:

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

Parser/Evaluatorはsource spanとdocumentを保持し、函数定義経由のErrorにはcall traceを付与できる。

---

# 31. 性能方針

exact/certifiedはCPUのnative doubleより大幅に重い。
過去のmicrobenchmarkでは、対象によりdouble比で約100倍〜10万倍超の差がある。

それでも対話型CLIで数十µs〜数msの処理は実用上問題になりにくいため、通常意味論をdoubleへ落として速度を稼がない。

実施済み高速化例:

- Number real-real fast path
- Rational乗除算の重複GCD除去
- Rational加算の縮約範囲最小化
- Simplifier structural key再計算削減
- Add同類項索引
- BigInt cube root Newton法
- factorial balanced product tree
- certified Log range reduction / log(2) enclosure共有
- FFT radix-2

将来`for/Plot`のように数千〜数百万回の評価を行う処理では、Exact/Certifiedとは別に明示的Machine evaluatorを追加する予定。

---

# 32. 現在未実装・保留の主な項目

今後の候補と保留理由は`docs/roadmap.md`を参照する。

代表:

- `digamma`, `trigamma`, `zeta`, `ibeta`, `polylog`
- `isprime`, `nextprime`, `prevprime`, `factorint`, `totient`
- advanced LU/QR/SVD/eigen/condition number/least squares
- `hilbert`（旧仕様の名称再確認）
- `fma`, `clamp`, `proj`
- 工学函数、財務函数、単位変換
- 旧colon command `:defs`, `:help`, `:unset`, `:undef` 等（函数版`Defs[]/UnDef[]`は実装済み）。`:angle`は`angleMode[]`へ置換し、表示設定として`:fix`/`:status`を実装済み
- `for`, `plot`
- general Machine/double evaluation mode

今後追加検討:

今後はAlgebraicNumber / Root基盤、`rootApproximant`、周期函数の整数parameter付き解集合、Machine evaluator / `for` / `plot`を候補とする。

---

# 33. 現在のCLI

CLIはKernelの数学状態とfrontendの表示状態を分離する。起動時には次のオプションを指定できる。

```text
mmCal --fix 16 --angle deg
mmCal --angle rad
mmCal --angle grad --fix 8
```

- `--fix n`: 起動時の小数表示桁数上限。内部値は変更せず、末尾の不要な0は省略する
- `--angle deg|rad|grad`: 起動時の既定角度
- `--help`, `-h`: 使用法を表示

```text
In[1]> 1/3
Out[1]> 1/3
```

- 1行ごとにparse/evaluate
- promptは`In[n]>` / `Out[n]>`で固定し、余分な空白を入れない
- 終了は`Exit[]`に一本化。裸の`exit` / `quit`特別扱いはない
- `Clear[]`: user definitionsと全履歴を消し、次の入力番号を1へ戻す
- `Defs[]`, `UnDef[...]`: user definitionsの確認・削除
- 計算履歴 `%`, `%%`, ... および再評価型`In[n]`, snapshot型`Out[n]`

## 33.1 `:fix` — presentation-only小数表示

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

`:fix n`は小数点以下最大`n`桁へ丸める**表示だけ**を変更する。末尾の不要な0は省略するため、例えば`:fix 5`で`31/10`は`3.1`と表示する。保存されるExpr、`Out[n]`、`precision/accuracy`の意味論は変更しない。exact値をMachine/doubleへ変換する機能ではない。`n`は現在0..1000。`:fix`のみなら現在の表示modeを表示する。

数値としてcertifyできる式全体は表示時だけ近似する。自由変数を含むsymbolic expressionはexact表記を維持する。

## 33.2 `:status`

```text
:status
Angle: Rad
Display: Exact
Evaluation: Exact-first
Definitions: 0
History: 0
```

`:status`も履歴へ入らないCLI commandである。数学状態の変更は`angleMode[...]`等のKernel函数、presentation状態の変更は`:fix`等のCLI command、という境界を維持する。

## 33.3 console title

タイトルは補助情報として、例えば次の形式に更新する。

```text
mmCal 1.5.0 - Rad - Exact
mmCal 1.5.0 - Deg - Fixed(16)
```

- Windows: `SetConsoleTitleA`
- Linux/macOS: TTY時のみANSI OSC title sequence
- その他: no-op

タイトル変更失敗は計算Errorにしない。状態確認の正本は`:status`であり、terminalがtitleを上書きしても意味論には影響しない。

## 33.4 canonical formatter

通常`Out[n]`はAST dumpではなく、再parse可能なcompact数学表記とする。

```text
x^2+sin[x]
A-B+C
2(x+sqrt[x])sqrt[x+sqrt[x]]/3
```

- `+ -` / `+-`は出さず、負項を`-`として表示
- `A-(B-C)`のような加減算は表示時だけ`A-B+C`へflattenできる
- `+`, `-`, `*`, `/`, `^`, 比較演算子の前後に不要な空白を置かない
- implicit multiplicationは字句上安全な場合だけ連結する（`2x`, `2sqrt[x]`）。`2exp[x]`や`2E`のように指数表記と衝突する連結は`2*exp[x]`, `2*E`と明示する
- identifier同士など連結で別tokenになる場合は必要な空白を残す（`I Pi`, `x y`）
- 数字同士など曖昧になる場合は空白ではなく明示`*`を使う
- precedence/associativityを守り、format → parse → formatで意味が変わらないことを回帰テストする

内部構造を見せるdebug/full-form表示は、通常formatterとは将来別機能に分離する。

---

# 34. 実装上の主要層

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

主要責務:

- `SymbolTable`: intern / identity
- `SymbolRegistry`: protected symbol / constant / domain name
- `BuiltinRegistry`: name / alias / arity / Hold属性
- `MathRegistry`: domain / parity / branch / definedness
- `ValueFacts`: conservative numeric-domain/sign inference
- `KnowledgeContext`: permanent facts + assumptions
- `Simplifier`: safe local rewrite
- `FullSimplifier`: bounded candidate search
- `CertifiedEvaluator`: expression全体のinterval evaluation
- `SolutionSet`: Solverの解集合表現
- `RandomEngine`: session-local stateful PRNG

この分離を維持し、函数追加ごとにSolver・Simplifier・数値backendへ同じ知識を重複記述しないことを基本方針とする。


## CertifiedEvaluator の安全限界

CertifiedEvaluatorは式を再帰的に区間評価するため、病的に深いASTについてはOSのstack overflowへ到達する前に評価対象外として扱う。現在の深さ上限は96段。通常のn-ary `Add` / `Multiply` の項数ではなく、ASTの入れ子深さに対する安全弁である。

`nintegrate` は高階導函数を作る前に被積分函数を区間全体でpreflightし、明白な特異点を先にDomainErrorへ落とす。
