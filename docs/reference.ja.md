# mmCal 仕様・函数リファレンス

この文書は **mmCalの実装そのもの** を基準にした詳細仕様書である。
ユーザー向けの導入はルートの`README.ja.md`を参照する。
この文書はバージョンごとに常に変動するため，過去バージョンはgitより引っ張り出してください。

## 0. 不変の理念

### 設計理念
厳密に。全て自前で。近似は明示的に。解らないものは解らないと言う。

### 配布理念
唯一本の実行ファイルに。
Open source under the BSD 3-Clause License.

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
-> 3.14159265358979323846264338328
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

任意精度の作業値は`BigFloat`，証明付き区間は`RealInterval` / `ComplexInterval`。

`N[expr,p]`では，真値を含む区間の両端が同じ`p`有効桁の10進丸めへ入ることを確認してから`DecimalApproximation`を返す。0近傍で相対Precisionを定義できない場合でも，InformationEnclosureから絶対Accuracyを保証できるならzero-centered approximationを返せる。

現在の`DecimalApproximation`は表示文字列だけではなく，要求桁数，由来（exact入力 / certified interval），表示10進値そのもののexact Rationalに加えて，**CertifiedEnclosure**と**InformationEnclosure**の2種類のexact Rational区間を保持する。CertifiedEnclosureは真値包含を証明する区間，InformationEnclosureはその近似値から後続計算で利用してよい情報量を表す区間であり，常に`CertifiedEnclosure ⊆ InformationEnclosure`を満たす。`ComplexDecimalApproximation`も実部・虚部ごとに同じmetadataを保持する。`precision/accuracy/rationalize`はInformationEnclosureを直接使い，表示文字列を再parseして精度を推測しない。

```text
N[sqrt[2],30]
-> 1.4142135623730954881688724210
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

函数呼出には**角括弧 `[]` だけ**を使用する。丸括弧 `()` はgrouping専用であり、函数呼出delimiterにはしない。

```text
sin[Pi/6]
sqrt[2]
f[x]
```

既知の函数名に `sin(x)` のような旧丸括弧構文を使うとSyntaxError。通常identifierの `x(x+1)` は暗黙乗算として受理するが、Formatterは `x*(x+1)` と明示する。グルーピングは丸括弧を使い、単独の`[x+1]`はgroupではない。ユーザー函数の定義も `f[x] := ...` の形式に限定する。

## 4.2 配列

```text
{1,2,3}
{{1,2},{3,4}}
```

内部ではdense `ArrayExpr`として扱う。Unreleased実装ではnumeric値をimmutableなpacked pageへ保持し，shape / offset / stridesを別に持つため，transposeや一部reshape/sliceはbackingを共有できる。これは内部最適化であり，ユーザーからは通常のArrayとして見える。

## 4.3 変数・ユーザー函数

```text
x := 3
-> 3

f[t] := t^2 + 1
f[4]
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

短縮記法は次のとおり。

```text
@       // In [-1]
@@      // In [-2]
@@@     // In [-3]
%       // Out[-1]
%%      // Out[-2]
%%%     // Out[-3]
```

正式な参照は`In [n]` / `Out[n]`を使う。`n > 0`は画面上の絶対入力番号，`n < 0`は相対参照，`n = 0`はTypeErrorとする。

```text
In [1]
Out[1]
In [-1]
Out[-1]
```

`In [n]`は対象入力のlowered Exprを取得した後，**現在のsession環境で通常評価する**。正の`In [n]`は絶対入力番号，負の`In [-n]`は現在評価中の入力slotを除いて過去の入力slotを数える。したがって`@` / `@@` / `@@@` / ... はそれぞれ`In [-1]` / `In [-2]` / `In [-3]` / ... を意味する。連続個数に固定上限はない。評価エラーになった入力でもparse/lowerまで成功していれば`In [-1]`で再評価できる。一方，Lexer / Parser / Lowererで式として成立しなかった入力は履歴へcommitせず，`In[n]`番号も消費しない。エラー表示にはその時点のpending入力番号を使うため，修正後の次入力は同じ`In[n]`から再開する。

`Out[n]`は保存済み結果snapshotを返し，再評価しない。正の`Out[n]`は絶対入力番号に対応するsnapshot，負の`Out[-n]`は**成功した出力だけ**を直前から数える。このため評価失敗を挟んでも常に`% == Out[-1]`，`%% == Out[-2]`，`%%% == Out[-3]`，... となる。`%`も連続個数に固定上限はない。

```text
In [1]> fft[{1,2,3}]
Out[1]> {6, ...}
In [2]> N[@,30]
Out[2]> {6, -1.50+0.866025403784...I, ...}
```

上の`N[@,30]`では，`Out[1]`を後から30桁化するのではなく，`In [1]`の`fft[...]`を30桁のapproximation contextで再評価できる。

絶対参照の例：

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

したがって`In [n]`は「過去入力を現在環境へ貼り戻して再実行する」意味である。過去入力が変数参照・代入・乱数等を含めば現在の定義やRNG stateを使う。生の入力ASTを表示する用途とは分離する。正の添字で現在評価中の入力自身を参照することは禁止し，自己再帰によるstack overflowを防ぐ。

評価エラーになった入力でもparse/lowerまで成功していれば絶対入力slot自体は残るが，対応する正の`Out[n]`は存在しない。

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
-> 2.6789385347077476337
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

## 14.5 Fresnel C / S

mmCalでは標準Fresnel積分を

```text
fresnelc[x] = integral_0^x cos[Pi t^2/2] dt
fresnels[x] = integral_0^x sin[Pi t^2/2] dt
```

に対応するentireな奇函数として扱う。exactに閉じない引数は記号式を保持し、`N`では誤差保証付き実数評価を行う。

```text
fresnelc[0] -> 0
fresnels[0] -> 0
N[fresnelc[1],20] -> 0.77989340037682282947
N[fresnels[1],20] -> 0.43825914739035476608

D[fresnelc[x],x] -> cos[Pi x^2/2 Rad]
D[fresnels[x],x] -> sin[Pi x^2/2 Rad]
```

微分の位相には`Rad`を明示する。Fresnel函数の定義自体はsessionの既定角度単位に依存しないためである。

## 14.6 合流型超幾何函数 1F1

Kummerの合流型超幾何函数を

```text
hypergeometric1F1[a,b,z]
```

で表す。`z`についてentireであり、`b = 0,-1,-2,...`には一般にparameter poleがあるため、その場合を無条件に有限値へ簡約しない。現段階のexact評価は停止する級数、`z=0`、`a=b`等の安全に閉じる場合を扱う。`N`の誤差保証付き実数backendはexact Rationalの`a,b,z`を対象とする。

```text
hypergeometric1F1[0,3,2] -> 1
hypergeometric1F1[-2,3,2] -> 0
hypergeometric1F1[2,2,1] -> E
N[hypergeometric1F1[1/6,7/6,1],20]
-> 1.1920688079818883008
```

parameterが微分変数に依存しないとき、

```text
D[hypergeometric1F1[a,b,z],z]
= a hypergeometric1F1[a+1,b+1,z]/b
```

を使う。積分器では、上側不完全Gammaによる局所式がprincipal branchや原点のremovable holeを持つ場合に、原点を含めてentireな1F1表現を優先する。例えば、

```text
integrate[exp[x^6],x]
-> x hypergeometric1F1[1/6, 7/6, x^6]
```

より一般に正整数`n`について`exp[c x^n]`を同じ系列へ還元できる。

## 14.7 Gauss超幾何函数 2F1

Gaussの超幾何函数を

```text
hypergeometric2F1[a,b,c,z]
```

で表す。`c = 0,-1,-2,...`には一般にparameter poleがあり、`z`についてはprincipal branchを採用する。現段階のexact評価は、上側parameterが非正整数で停止する有限級数や、`a=0` / `b=0`など安全に閉じる場合を扱う。誤差保証付き`N` backendはexact Rational parameterと`|z|<1`の実数引数を対象とする。

```text
hypergeometric2F1[-2,1,3,1/2] -> 17/24
hypergeometric2F1[0,2,3,x] -> 1
N[hypergeometric2F1[1/2,1/2,3/2,1/4],20]
-> 1.0471975511965977462
```

parameterが微分変数に依存しない場合、

```text
D[hypergeometric2F1[a,b,c,z],z]
= a b hypergeometric2F1[a+1,b+1,c+1,z]/c
```

を使う。積分器では、例えば

```text
integrate[sqrt[1+2x^3],x]
-> x hypergeometric2F1[-1/2, 1/3, 4/3, -2x^3]

integrate[1/(1+x^5),x]
-> x hypergeometric2F1[1, 1/5, 6/5, -x^5]
```

のようなbinomial-power familyへ利用する。一般の2F1を`Solve`で逆函数化する規則は持たない。大域単射性を証明できないためであり、停止級数やexact退化で既存代数式へ落ちた場合だけ通常のSolverへ渡す。

## 14.8 不完全楕円積分 F / E / Pi

mmCalではLegendre形の不完全楕円積分を

```text
ellipticF[phi,m]
ellipticE[phi,m]
ellipticPi[n,phi,m]
```

で表す。第2引数`m`はparameterであり、振幅`phi`は**常にRadian**として解釈する。sessionの`Deg/Rad/Grad`設定には依存しない。principal branchを採用し、一般complex parameterのbranch cutやpoleを単純な「everywhere defined」には扱わない。

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
```

現段階のcertified real backendは、exact Rational振幅について`F/E`では`|m|<1`、`Pi`ではさらに`|n|<1`の安全な領域を扱う。振幅微分は

```text
D[ellipticF[phi,m],phi]
= 1/sqrt[1-m sin[phi Rad]^2]

D[ellipticE[phi,m],phi]
= sqrt[1-m sin[phi Rad]^2]

D[ellipticPi[n,phi,m],phi]
= 1/((1-n sin[phi Rad]^2)sqrt[1-m sin[phi Rad]^2])
```

である。したがって標準kernelは直接積分できる。

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

最後のquartic reductionは正しい局所primitiveだが、現在の`fullSimplify`は`sin[asin[x]]`とprincipal square rootの積を一般に安全な恒等式へ潰し切れない。そのためderivative-back harnessではResolutionOnlyとして監視し、証明器不足を理由に積分能力を削らない。一般の楕円函数方程式も、逆楕円函数族をまだ持たないため`Solve`は未解決を保持する。`m=0`等でexactに通常式へ退化した場合だけ既存Solverが解く。


## 14.9 Ei / Si / Ci / li / Polylogarithm

積分で頻出するprincipal special functionsを次の名前で表す。

```text
Ei[x]
Si[x]
Ci[x]
li[x]
polylog[s,z]
```

`Ei`, `Ci`, `li`, `polylog`は一般にbranchを持つため、MathRegistryではprincipal branchとして扱う。`Si`はentireな奇函数である。現在のcertified real backendは、安全にtail boundを証明できる領域に限定し、対応外を推測値で埋めない。

```text
Si[0] -> 0
Si[-1] -> -Si[1]
polylog[0,z] -> z/(1-z)
polylog[1,z] -> -log[1-z]
polylog[2,1] -> Pi^2/6
polylog[2,-1] -> -Pi^2/12

N[Ei[1],20] -> 1.8951178163559367555
N[Si[1],20] -> 0.94608307036718301494
N[Ci[1],20] -> 0.33740392290096813466
N[li[2],20] -> 1.0451637801174927848
N[polylog[2,1/2],20] -> 0.5822405264650125059
```

現在の微分Knowledgeは、引数・order parameterが微分変数に依存しない範囲で

```text
D[Ei[x],x] -> exp[x]/x
D[Si[x],x] -> sin[x]/x
D[Ci[x],x] -> cos[x]/x
D[li[x],x] -> 1/log[x]
D[polylog[s,x],x] -> polylog[s-1,x]/x
```

を使う。`polylog[2,x]`は`polylog[1,x]`を`-log[1-x]`へexact退化させるため、

```text
D[polylog[2,x],x] -> -log[1-x]/x
```

まで閉じる。積分器ではこの共有Knowledgeにより、

```text
integrate[exp[x]/x,x] -> Ei[x]
integrate[sin[x]/x,x] -> Si[x]
integrate[cos[x]/x,x] -> Ci[x]
integrate[1/log[x],x] -> li[x]
integrate[log[1-x]/x,x] -> -polylog[2,x]
```

を返す。一般の`Ei/Si/Ci/li/polylog`方程式にprincipal inverseを一個だけ返す`Solve`規則は持たない。大域単射性・branchを証明できないためであり、`polylog[0,z]`や`polylog[1,z]`のように既存の代数函数・`log`へexact退化した場合だけ既存Solverへ渡す。

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

## 17.1 Array基盤

Arrayはrankごとに別Value型を増やさず，共通のdense `ArrayExpr`を使う。Unreleasedではphysical storageをimmutable packed page，logical layoutをshape / offset / stridesへ分離している。
公開リテラルは従来どおり `{...}` / `{{...},...}` とし，Matrix演算へ渡せるArrayは常にdense rectangularである。

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

indexは0始まり。`at`はrank未満のprefix indexも受け取り，残り次元を保持したsubarrayを返す。全rank分を指定した場合だけscalarになる。

```text
dimensions[{{1,2,3},{4,5,6}}] -> {2,3}
arrayRank[{{1,2},{3,4}}] -> 2
at[{{1,2},{3,4}},1] -> {3,4}
at[{{1,2},{3,4}},1,0] -> 3
reshape[{1,2,3,4},{2,2}] -> {{1,2},{3,4}}
```

`{...}`はユーザー意味論として一般の有限brace containerである。child shapeが全て一致する矩形値は内部でdense `ArrayExpr`へ自動昇格し，`{{1,2},{3}}`や分解結果の`{Q,R}`のようにshapeが揃わない値は一般braceのまま保持する。numeric dense Arrayは内部でInteger / Rational / Number等のpacked pageを共有する場合があるが，storage種別はユーザー意味論へ露出しない。一般brace自体は正常な値であり，Matrix函数へ渡した時点で矩形性監査が入り，非矩形ならWarningを出して未評価保持する。`dimensions` / `arrayRank`は非矩形値では全childに共通するrectangular prefixだけを返す。

```text
dimensions[{{1,2},{3}}] -> {2}
arrayRank[{{1,2},{3}}] -> 1
length[{{1,2},{3}}] -> 2
at[{{1,2},{3}},0] -> {1,2}
transpose[{{1,2},{3}}] -> Warning + unevaluated
```

`mget[A,row,col]` は互換aliasとして `at` と同じ0始まりindexを使う。

先頭側に0長次元を持つArrayはbrace literalだけではshapeを復元できないため、Formatterは必要な場合だけ `reshape` を使う。

```text
zeros[0,3]
-> reshape[{}, {0, 3}]

dimensions[zeros[0,3]]
-> {0,3}
```

評価後にArray要素がArrayへ変わる場合も、同一shapeなら自動的に一段flattenして共通Arrayへ正規化する。scalar/Array混在またはchild shape不一致はTypeError。

## 17.2 exact-first線形代数

canonical API:

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
eigenvalues[A]
eigenvectors[A]
eigensystem[A]
norm[v]
normalize[v]
trace[A]
```

`dot` はStage 2ではrank-1/rank-2を扱う。

```text
dot[{1,2,3},{4,5,6}] -> 32
dot[{{1,2},{3,4}},{5,6}] -> {17,39}
dot[{5,6},{{1,2},{3,4}}] -> {23,34}
dot[{{1,2},{3,4}},{{5,6},{7,8}}] -> {{19,22},{43,50}}
```

Array同士の `*` は行列積にしない。`*` はscalar×Arrayだけを許し、行列積・vector contractionは明示的に `dot` を使う。同shapeの `+/-` はelement-wise。

exact実数/Rational行列は各行の分母を払って整数行列へliftし、Bareiss fraction-free eliminationを使う。これによりpivotごとのRational生成を避ける。exact complexはflat `Number` Gaussian backendへfallbackする。symbolic行列は非零性を証明できないpivotを勝手に選ばない。

```text
det[{{1,2},{3,4}}] -> -2
inverse[{{1,2},{3,4}}] -> {{-2,1},{3/2,-1/2}}
rref[{{1,2},{3,4}}] -> {{1,0},{0,1}}
matrixRank[{{1,2},{2,4}}] -> 1
nullSpace[{{1,2},{2,4}}] -> {{-2,1}}
solveLinear[{{2,1},{1,-1}},{5,1}] -> {2,1}
```

`nullSpace[A]`はRREFのfree columnを昇順に取り，各free variableを1としたcanonical basisを返す。返り値shapeは `{nullity, columns}` であり，full column rankでは `reshape[{}, {0,n}]` として空basisのvector次元を保持する。exact整数/RationalではBareiss forward eliminationを共有し，exact complexはGaussian fallback，symbolicではpivotの非零性を証明できる場合だけbasisを構成する。

`solveLinear[A,b]` は `A` を m×n 行列、`b` を長さmのvectorとして扱う。一意解が存在すれば長さnのvectorを返す。正方行列に限定せず、整合した過剰決定系もfull column rankなら解ける。不整合系、または自由変数が残る系はDomain error。一般parametric solutionはこの函数では捏造しない。

`luDecomposition[A]` は現在正方行列を対象とし，shape `{3,n,n}` の `{P,L,U}` を返す。規約は `P A = L U`。certified approximate LUでは，非零を証明できた候補のうち`|pivot|^2`の区間下限が最大の行を選ぶpartial pivotingを使い，epsilon判定は行わない。row pivotingを行い，exact NumberではRational/complexをexactに保持する。三角symbolic行列は不要な除算を行わずそのまま分解でき，非零性を証明できないpivotが必要な一般symbolic行列は未評価に留める。factorはprefix indexingで取り出せる。

```text
lu = luDecomposition[A]
at[lu,0] -> P
at[lu,1] -> L
at[lu,2] -> U
```

`qrDecomposition[A]` はHouseholder reflectorを使うreduced QRで，矩形m×nにも対応する。`k=min(m,n)`として `Q:m×k`，`R:k×n` を一般brace `{Q,R}` で返し，規約は `A = Q R`。factor shapeが同じ正方caseでは内部的にdense Arrayへ自動最適化されるが，ユーザー構文は同じ`{Q,R}`である。exact実数行列の一般Householder展開はexpression growthを避けるため3×3以下に制限し，安全な上三角/上台形caseはfast pathを使う。`N[qrDecomposition[A],p]`はexact QRを先に展開せず，実/複素ともcertified interval Householder backendへ直接入る。

```text
qr = qrDecomposition[A]
at[qr,0] -> Q
at[qr,1] -> R
```

Householder適用には複数列を一度のrow-major走査で処理できるcolumn-block kernelも実装している。ただし外部BLASを使わない現backendでは8/16/24次の実測で一貫した高速化が得られなかったため，自動block化は採用せずunblocked相当を既定とする。block kernelとbenchmarkは今後のBigFloat/Matrix backend最適化用に残す。

`svd[A]`はreduced SVDを `{U,S,V}` で返す。m×n入力に対して`k=min(m,n)`，`U:m×k`，`S:k×k`，`V:n×k`。実数なら `A = U S Transpose[V]`，複素数なら `A = U S conjugateTranspose[V]`。一般数値backendは条件数を二乗する`A^H A`を形成せず，Householder bidiagonalizationの後にone-sided Jacobiで列を直交化する。候補factorはreconstruction residualと`U^H U` / `V^H V`の直交性を区間演算で要求表示桁より厳しく監査し，証明できなければguard digitsを増やして再試行する。exact SVDは自然に閉じる実対角等へ限定する。重複特異値の部分空間ではsingular vector basisは一意ではないため，componentごとの「唯一の真値」を主張せず，再構成・直交性を保証する。

`eigenvalues[A]` / `eigenvectors[A]` / `eigensystem[A]` は正方行列の固有値・固有vector・組を扱う。`eigenvectors`の各**列**が対応する固有vectorであり，`eigensystem[A]`は `{values,vectors}` を返す。exact pathは上三角行列の対角固有値，対角行列の標準基底，およびdistinct eigenvalueを持つexact Number 2×2を明示処理する。重根を持つ非対角2×2では不足する固有vectorを複製せず未評価に留める。一般行列の `N[...]` はComplex BigFloat Hessenberg reduction + implicit shifted QRからSchur形 `A Q ≈ Q T` を求め，Schur三角行列からback substitutionで固有vectorを構成する。元入力のcertified intervalに対するSchur relationと `A v ≈ λ v` residual，およびSchur vectorのunitarityを要求表示桁より厳しく区間監査し，証明できなければguard digitsを増やして再試行する。一般非正規行列では固有値・固有vectorは摂動に敏感であり，返した各componentが唯一の真値を個別区間包含するとは主張しない。保証対象は計算されたSchur/eigenpair relationである。近接重根・defective caseで独立固有vectorを安定に構成できない場合，`eigenvectors` / `eigensystem`は推測せず未評価に留める。

`conjugateTranspose[A]`はHermitian transposeであり，complex SVDの`V^H`や複素直交性の検証に使う。rank-1では成分の共役だけを行い，rank-2では転置と共役を同時に行う。

`norm` は複素vectorに対してHermitian normを使う。

```text
norm[{3,4}] -> 5
norm[{3+4I}] -> 5
normalize[{3,4}] -> {3/5,4/5}
```

### precision-aware `N`

FFTと同じ `ApproximationContext` / certified interval変換を共有する。したがって例えば

```text
N[dot[A,B],100]
N[det[A],100]
N[inverse[A],100]
N[rref[A],100]
N[solveLinear[A,b],100]
N[luDecomposition[A],100]
N[qrDecomposition[A],100]
N[svd[A],100]
N[eigenvalues[A],100]
N[eigensystem[A],100]
N[norm[v],100]
```

は、巨大なexact中間式を完成させてから近似するのではなく、対応するBigFloat/interval backendへ要求精度を渡して直接評価できる。`solveLinear`はaugmented interval eliminationでpivotと整合性を証明し、証明不能なcaseをepsilonで補わない。依存した過剰決定系ではinterval相関の消失により直接証明できない場合がある。`matrixRank`と`nullSpace`はrank deficiencyに依存する不連続演算なので，exact入力ではexact eliminationを優先する。近似入力では浮動小数の任意thresholdを使わず、区間からpivot構造を証明できる場合だけ結果を返し，rank deficiencyを推測しない。

## 17.3 互換Vector / Matrix函数

従来名は互換のため維持する。

```text
madd
matmul mmul
rank mrank
mget
vadd vsub vscalar
vdot vcross
vnorm vnormalize
vproject vangle
vmanhattan veuclidean
vreflect vreflect_axis
vsum
```

`matmul/mmul/vdot` は `dot`、`rank/mrank` は `matrixRank`、`vnorm` は `norm`、`vnormalize` は `normalize` へ束ねる。

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

exact入力では2冪長FFTはradix-2 Cooley–Tukey、非2冪長はexact DFTへfallbackする。通常の`fft[...]`は引き続きexact-firstであり、machine `double`へ暗黙変換しない。

`N[fft[v],p]`では`N`が第1引数を先にexact展開せず、要求精度`p`をFFTへ伝播する。FFT側は`ComplexInterval`/BigFloat端点で直接butterflyを行い、各出力成分が要求桁へ一意に丸められることを証明してから`DecimalApproximation`を返す。近似入力を含む`fft[v]`も同じbackendへdispatchする。

近似FFTでは2冪長をradix-2、十分大きい非2冪長をBluestein convolutionへ還元する。小さい非2冪はdirect DFTの定数項が小さいため、現在のbenchmarkでは96点未満をdirectとしている。この閾値はMSVC環境で`mmCal.Benchmarks`から再測定する前提の実装値である。

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
- 正のRational scaleを証明できる二次逆平方根型の`asin/asinh` primitive、およびexact Rational係数二次式`q(x)`の`sqrt[q(x)]` primitive
- `sin^m/cos^n`の有限Fourier reduction。積分器は正整数総次数256までを明示的に展開可能
- `sin[u]^(-n)` / `cos[u]^(-n)` (`1<=n<=256`) を `csc/sec` の標準漸化式で積分
- `tan/cot/sec/csc`の正整数冪 (`2<=n<=256`) を標準reduction formulaで積分
- 和・差・符号反転、積分変数に依存しない係数の線形性
- `exp/sin/cos/tan/cot/sec/csc`の安全な標準原始函数
- `sinh/cosh/tanh/coth/sech/csch`の安全な標準原始函数
- `log/log1p/expm1/sqrt/cbrt`。`log[x]/x`や`1/(x log[x])`は対数微分Knowledgeから認識
- `asin/acos/atan/asinh/acosh/atanh`
- `erf/erfc`
- `fresnelc/fresnels`。exact Rational係数の`sin/cos[a x^2+b x+c]`を平方完成して標準Fresnel積分へ還元。`Pi*x^2/2`の定義核も直接認識
- `hypergeometric1F1`。正整数`n>=2`の`exp[c x^n]`を原点でentireな1F1 primitiveへ還元
- exactな逆chain rule。`f'(x) f(x)^p`はD後の偶然の式形に依存せず構造的にも認識
- 多項式×`exp/sin/cos/sinh/cosh`に対する有限回のintegration by parts
- `exp[a x+b] sin/cos[c x+d]`型を連立一次式としてexact積分
- principal `sqrt[x]`を含む有理的な形への` t=sqrt[x] `局所置換、および`sqrt[q(sqrt[x])]`の二次根号class
- 共通引数を持つ`R(sin(theta),cos(theta))`に対するbounded Weierstrass置換`t=tan(theta/2)`。変換後は既存exact有理積分器へ渡す
- `log[1+beta*x^n]/x`を`polylog[2,-beta*x^n]`へ還元するdilogarithm Knowledge
- boundedな積×短い和の分配。ただし全体がexact chain-ruleで一発に閉じる場合はchain-ruleを優先し、能力退行を防ぐ
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

今回追加した置換・接続の代表例:

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
WARN: integrate partially evaluated the expression; remaining subintegral(s) are outside the current symbolic rule set
```

解けない部分だけを保持し、既に求まった項まで巻き戻さない。

### 未評価理由の診断

未評価を単一のWARNへ潰さず、現在は次のdiagnostic codeを区別する。

- `integrate::unsupported` — 現在のsymbolic rule setに解法がない。**閉形式が存在しないことを意味しない**。
- `integrate::partial` — 一部は積分済みだが、残るsubintegralが現在のrule外。
- `integrate::conditionsRequired` — domain / branch仮定不足で安全なprimitiveを選択できない。
- `integrate::noKnownClosedForm` — mmCalの現在の標準函数語彙で有限閉形式がない代表familyとして明示的に認識したもの。

例えば、

```text
integrate[gamma[x],x]
-> integrate[gamma[x],x]
WARN: mmCal has no implemented symbolic integration rule for this expression; this does not imply that no closed form exists

integrate[abs[x],x]
-> integrate[abs[x],x]
WARN: integrate needs additional domain or branch assumptions before it can choose a safe symbolic antiderivative

integrate[x^x,x]
-> integrate[x^x,x]
WARN: integrate recognized a family with no known finite closed form in mmCal's supported standard-function vocabulary; the integral remains unevaluated
```

`noKnownClosedForm`も「あらゆる数学的表現で不可能」という判定ではない。級数、新しい特殊函数、より広い函数classを許せば表現できる場合がある。mmCalが主張する範囲を現在サポートする有限標準函数語彙に限定する。

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

solve[{2x+3y==5,x-2y==9},{x,y}]
-> {{x==37/7, y==-13/7}}
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
N[expr,p]
```

`p`は**有効10進桁数(significant decimal digits)**であり，既定は16桁。小数点以下の表示桁数ではない。固定小数表示は`:fix` / `--fix`が担当する。
Arrayへ再帰的に適用できるほか、`arg`などが返す明示角度単位では値の部分だけを近似し、単位は保持する。

v1.5.2では`N`をprecision-aware evaluationの入口として扱う。第2引数の要求精度を先に確定し、第1引数の評価中はそのprecision contextを保持する。通常builtinは従来どおりexact評価され、FFTなど明示的に対応したbuiltinだけが要求精度を受け取って直接certified backendへ降りる。したがってexact-firstの意味論を全体へ暗黙に変更しない。

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

exact Rationalが有限10進になる場合，表示は必要以上に0埋めしない。例えば `N[1/2,10] -> 0.5` である。certified interval由来の有効桁結果では，要求桁に対応する末尾0の連続だけを圧縮し，最後に1個の0を残す。したがって内部の12桁保証が `1.000000000000` を確定していても表示は `1.0`，`1.500000000000` なら `1.50` とする。要求桁数，CertifiedEnclosure，InformationEnclosureはmetadataに全て保持し，表示上の0の個数を精度保証そのものとして扱わない。

## 24.1 CertifiedEnclosure / InformationEnclosure

`DecimalApproximation` / `ComplexDecimalApproximation`は，近似値ごとに2種類の区間を保持する。

- **CertifiedEnclosure** — 真値が必ず含まれることをbackendが証明した区間。内部guard桁により，ユーザーへ宣言した桁数より大幅に狭い場合がある。数学的な正しさ，zero判定，要求桁への一意丸め判定にはこちらを使う。
- **InformationEnclosure** — その近似値から後続計算で利用してよい情報量を表す区間。非zeroの`N[x,p]`で表示値`d`の10進指数を`e=floor(log10(|d|))`とすると，少なくとも`d ± 0.5*10^(e-p+1)`とCertifiedEnclosureの双方を包含する。したがって情報量は値のscaleに追従し，内部guard桁をユーザー可視のAccuracyとして後から回収しない。zero-centered approximationでは相対Precisionではなく，InformationEnclosureが直接absolute Accuracyを表す。

常に次を不変条件とする。

```text
CertifiedEnclosure ⊆ InformationEnclosure
```

`Infinity`は拡張実数sentinelであり，一般symbolの有限代数則を適用しない。特に`Infinity-Infinity`，`0*Infinity`，`Infinity/Infinity`は現段階では未評価に留め，`0`等を捏造しない。完全な拡張実数算術は別仕様として扱う。

InformationEnclosureは確率分布や統計的confidence intervalではない。また「真値がこの広い区間のどこにでもあり得る」とbackendが主張するものでもない。真値保証そのものはCertifiedEnclosureが担当し，InformationEnclosureは**現在の値から利用してよい情報量の契約**を表す。したがってbackendがより狭いCertifiedEnclosureを内部に持っていても，それだけを理由に既存近似値の情報量は増えない。

通常の`+ - * /`と単項`-`では2区間を独立に伝播する。exact `Number`は両方について同じpoint intervalとして混在できる。

Unreleasedでは`DecimalApproximation` / `ComplexDecimalApproximation`をcertified numerical evaluatorのfirst-class leafとして扱う。`sin` / `exp` / `log` / `sqrt` / 双曲線・逆函数・`gamma` / `erf` / `Ei` / `Si` / `Ci`等，interval backendを持つscalar函数では両enclosureを独立に伝播する。`log2` / `log10` / `fract`のようにprimitiveへrewriteされる函数もrewrite後に同じ経路へ入る。ordered comparison，`min` / `max`等の離散的判定は**InformationEnclosureだけで結論を証明できる場合**に限って確定し，内部guard桁をBoolean結果から漏らさない。exact Rational parameterだけを受ける現行`1F1` / `2F1` / elliptic / `polylog`の一部backend等は，approximate parameterへ無理に拡張せず未評価に留める。

```text
N[Pi,20] + 1/3

Certified:    C(Pi) + {1/3}
Information:  I(Pi) + {1/3}
```

出力10進値の正当性はCertifiedEnclosureから決定し，出力として宣言できる桁数はInformationEnclosureを越えない範囲へ制限する。scale拡大や近接減算ではInformationEnclosureも演算されるため，`accuracy` / `precision`は自然に低下し得る。演算結果自身にも伝播後のInformationEnclosureを保存するので，複数回の演算を跨いでも単なる「要求桁数」へ情報を圧縮し直さない。

外側の`N`はInformationEnclosureを狭めて情報を発明しない。したがって

```text
N[N[Pi,20],100]
-> 3.1415926535897932385
```

は元の20桁保証を保持する。一方，より低い桁を要求した場合は表示丸めに対応するInformationEnclosureを追加して安全に情報を捨てられる。これらは`double`等のmachine arithmeticへ変換せず，両enclosureを`RealInterval` / `ComplexInterval`へ持ち上げて外向き丸めで計算する。

---

# 25. precision / accuracy / rationalize

## 25.1 `accuracy[x]`

`DecimalApproximation`について，真値に対する**保証可能な絶対10進桁数の整数下限**を返す。

```text
accuracy[N[1/3,20]]
-> 20
```

表示値`d`とInformationEnclosure `[iL,iU]`から

```text
max(|d-iL|, |d-iU|)
```

をabsolute error boundとして使う。`N[...,p]`生成時のInformationEnclosureには有効桁丸めに対応するscale依存の半量子が既に含まれるため，有限小数がCertifiedEnclosure上で真値と偶然完全一致していても，要求桁を越えたhidden guard情報をAccuracyとして回収しない。

exactな数・exact symbolic expressionは`Infinity`を返す。

```text
accuracy[1/3] -> Infinity
accuracy[Pi]  -> Infinity
```

## 25.2 `precision[x]`

同じabsolute error boundを，InformationEnclosureから得られる値絶対値の正の下限で割り，**保証可能な相対10進桁数の整数下限**を返す。

```text
precision[N[1/3,20]]
-> 19
```

これは要求した20有効桁を機械的に返す函数ではない。`1/3`近傍では有効20桁の丸め量子が`10^-20`なので，InformationEnclosureの相対的不確かさから保証できる整数桁数は19になる。InformationEnclosureが0を含む場合は値絶対値の正の下限を得られないため0を返す。exact expressionは`Infinity`。近接減算ではabsolute Accuracyを多く残したまま結果scaleだけが小さくなるため，Precisionだけが大きく落ちることがある。

## 25.3 `rationalize[x]`

近似値が持つInformationEnclosure内から，**分母が最小になるexact Rational**を求める。探索はexact Rational上のcontinued-fraction型interval recursionで行い，doubleへ変換しない。これにより，CertifiedEnclosureだけが保持しているhidden guard桁やhidden exact pointから，ユーザーへ宣言していない情報を`rationalize`で掘り返さない。

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

## 25.4 `explain[value]`

評価済みの値が**既に保持している情報だけ**を構造化して返す軽量introspection函数。対象を改めて`det` / `matrixRank` / LU / Eigen等へ掛けたり，Generic Arrayを全走査して性質を推論したりはしない。したがって`explain`自体のために高価な数学計算は開始しない。

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

通常の引数評価を先に行うため，`explain[1+2]`は`3`を説明する。履歴結果もそのまま対象にできる。

```text
explain[Out[25]]
explain[%]
```

返値はbrace構文で表示されるproperty/value pair列。`Dimensions`や各種Enclosureの値自体がArray/braceになり得るため，内部表現はdense `ArrayExpr`ではなく一般`ListExpr`を用いる。表示上は`{{"Kind",...},{...}}`であり，情報を文字列へ潰さない。

`Exactness`はbooleanではなく分類値を返す。現在の主な値は`"Exact"` / `"CertifiedApproximation"` / `"Unknown"`。近似値を単に`Exact=False`とは表現しない。

組込み数学定数・予約symbolは一般の未知symbolへ落とさない。`Pi/E/Phi`はMathRegistryに既に登録された数学metadataをO(1)で参照し，`Infinity`等はSymbolRegistryの予約意味を用いる。

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

`I`は通常評価でexact complex `Number`へloweringされるため，`explain[I]`は入力tokenではなく評価後の複素数値を説明する。

certified decimal approximationでは，要求有効桁数に加えてCertifiedEnclosureとInformationEnclosureを別々に確認できる。前者は真値保証，後者は後続計算で利用してよい情報量である。

```text
explain[N[Pi,20]]
-> {{"Kind","DecimalApproximation"},
    {"Domain","Real"},
    {"Exactness","CertifiedApproximation"},
    {"RequestedPrecisionDigits",20},
    ...
    {"CertifiedEnclosure",{certifiedLower,certifiedUpper}},
    {"InformationEnclosure",{informationLower,informationUpper}}}
```

Arrayではstorage metadataからO(1)で分かる`Domain` / `Exactness`と，shapeからほぼ無料で分かる`ArrayRank` / `Dimensions` / `ElementCount` / `Vector` / `Matrix` / `Square` / `Order` / `Empty`を返す。`det`，数学的`matrixRank`，invertibility，eigenvalue等は返さない。必要なら既存函数を明示的に呼ぶ。

整数では符号・zero・`BitLength`を返す。10進桁数は巨大整数で10進変換を必要とするため，軽量性を優先して自動計算しない。Rationalでは分子・分母のbit長を返す。

```text
explain[value,"internal"]
```

は開発・性能調査用で，`Representation`，Arrayの`Storage` / `Contiguous` / `StoredExpressions`，近似値の`ApproximationOrigin`等を追加する。**`"internal"`のproperty名・値は互換性保証対象ではない。** 未知modeはerrorとし，将来高コストな`"full"`相当を暗黙に実行しない。

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
| `matmul`, `mmul`, `vdot` | `dot` |
| `mtranspose` | `transpose` |
| `mget` | `at` |
| `mdet` | `det` |
| `minverse` | `inverse` |
| `rank`, `mrank` | `matrixRank` |
| `mtrace` | `trace` |
| `mrows` | `rows` |
| `mcols` | `cols` |
| `mdiag` | `diag` |
| `vnorm`, `vlength` | `norm` |
| `vdistance` | `veuclidean` |
| `vnormalize`, `vunit` | `normalize` |

aliasは別実装ではなく同一`BuiltinId`へ束ねる。数学metadataやSolver規則を二重管理しない。

mmCal 1.5.0では、Mathematica互換だけを目的とした大文字始まりalias（`Sin`, `ArcTan`, `Integrate`, `Solve`等）を削除した。数学函数はlowercase canonicalを原則とする。`D`, `N`, `In`, `Out`, `Exit`, `Clear`, `Defs`, `UnDef`は記号演算・Kernel操作の固有名として例外的に維持する。互換構文が必要になった場合は、default namespaceへaliasを増やすのではなく独立したimport/compatibility層として検討する。

---

# 29. 現在のsource-callable函数一覧

現在の開発treeでは **builtin/alias登録名237個 / sourceから呼出可能な名前219個**。内部headはsource-callable数に含めない。

```text
Clear, D, Defs, DtoG, DtoR, Exit, GtoD, GtoR, In, N,
Out, RtoD, RtoG, UnDef, abs, accuracy, acos, acosh, angleMode, arg,
arrayRank, asin, asinh, at, atan, atan2, atanh, ave, beta, betaln, binom, cbrt,
ceil, choice, cis, collect, cols, comb, conj, conjugateTranspose, convolve, corr, corrspearman,
cos, cosc, cosh, cot, coth, cov, csc, csch, csgn, cv,
det, dft, diag, diff, dimensions, dot, eigenvalues, eigenvectors, eigensystem, element, erf, erfc, exp, explain, expand, expc,
Ei, Si, Ci, li, polylog, fresnelc, fresnels, hypergeometric1F1, hypergeometric2F1, ellipticF, ellipticE, ellipticPi,
expm1, fact, factor, fallingfact, fft, fib, floor, frac, fract, fullSimplify,
gamma, gcd, geomean, harmmean, hypot, identity, if, ifft, im, imag,
integrate, inverse, iqr, kurtp, kurts, lcm, length, lgamma, limit, ln, log,
log10, log1p, log2, mad, madR, madd, mag, matmul, max, mcols,
mdet, mdiag, matrixRank, mean, median, mget, min, minverse, mmul, mod, mode,
luDecomposition, mrank, mrows, mtrace, mtranspose, nextpow2, nintegrate, norm, normalize, nullSpace, percentile, percentrank, perm, polar,
pow, precision, prod, quantile, quotient, rand, randSeed, randint, randn, rank,
qrDecomposition, rationalize, re, real, rect, rem, reshape, risingfact, rms, round, rows, rref,
sec, sech, sign, simplify, sin, sinc, sinh, sinhc, skew, solve, solveLinear,
singularValueDecomposition, sqrt, stddev, stddevs, stderr, sum, svd, tan, tanc, tanh, tanhc, trace,
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

評価自体は成功したがalgorithmic builtinが処理を完了できない場合は、結果Exprと別にWarningを返す。`D`, `solve`, `solveLinear`, `N`, `rref`, `matrixRank`, `nullSpace`, `luDecomposition`, `qrDecomposition`, `svd`, `eigenvalues`, `eigenvectors`, `eigensystem`（互換alias `rank`を含む）, `integrate`に加え、`precision/accuracy/rationalize`が対象外入力を未評価保持する場合もWarningになる。

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

- `digamma`, `trigamma`, `zeta`, `ibeta`
- `isprime`, `nextprime`, `prevprime`, `factorint`, `totient`
- condition number / least squares / 一般parametric linear system
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
In [1]> 1/3
Out[1]> 1/3
```

- 1行ごとにparse/evaluate
- promptは`In [n]>` / `Out[n]>`で固定し、余分な空白を入れない
- 終了は`Exit[]`に一本化。裸の`exit` / `quit`特別扱いはない
- `Clear[]`: user definitionsと全履歴を消し、次の入力番号を1へ戻す
- `Defs[]`, `UnDef[...]`: user definitionsの確認・削除
- 計算履歴 `@`, `%`, `%%`, ... および正負添字を持つ再評価型`In [n]`, snapshot型`Out[n]`

## 33.1 `:fix` — presentation-only小数表示

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
