# mmCal Cheatsheet

現行開発tree基準。詳細は `docs/reference.ja.md`。その場での確認は `:help <函数名>`。

## 1. 入力規則

| 用途 | 書式 | 例 |
|---|---|---|
| 函数呼び出し | `name[...]` | `sin[Pi/6]` |
| grouping | `(...)` | `(x+1)^2` |
| Array | `{...}` | `{1,2,3}` |
| Matrix | `{{...},...}` | `{{1,2},{3,4}}` |
| 代入 | `:=` | `x:=3` |
| ユーザー函数 | `f[x]:=...` | `f[x]:=x^2+1` |
| 数値近似 | `N[expr,p]` | `N[Pi,30]` |

**函数呼び出しは `[]` のみ。`()` はgrouping専用。**

```text
sin[Pi/6]
(x+1)^2
```

`sin(Pi/6)` は既知函数名に対する旧構文なので `SyntaxError`。通常identifierの `x(x+1)` は暗黙乗算 `x*(x+1)` と解釈される。

暗黙乗算も使える。

```text
2Pi
2sqrt[2]
2(x+1)
(x+1)(x-1)
```

演算子:

```text
+  -  *  /  ^  !
```

```text
2^3^2   -> 512
-2^2    -> -4
```

比較:

```text
==  !=  <  <=  >  >=
```

## 2. 数値・定数

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

有限小数は原則exact Rational。虚数単位は大文字 `I`。

基数付きliteral:

```text
0b1010
0o17
0xFF
2#1010
16#FF
```

定義域symbol:

```text
Integer  Rational  Real  Complex
```

## 3. 角度

**既定は Radian (`Rad`)。**

```text
sin[Pi/6]       -> 1/2
sin[30 Deg]     -> 1/2
sin[Pi/6 Rad]   -> 1/2
sin[100 Grad]   -> 1
```

session設定:

```text
angleMode[]        -> Rad
angleMode[Deg]     -> Deg
angleMode[Rad]     -> Rad
angleMode[Grad]    -> Grad
```

明示単位はsession設定より優先する。一般式にも付けられる。

```text
sin[x Deg]
sin[(2x+1) Grad]
```

変換:

```text
DtoR[180] -> Pi    DtoG[90]  -> 100
RtoD[Pi]  -> 180   RtoG[Pi]  -> 200
GtoD[200] -> 180   GtoR[200] -> Pi
```

## 4. exact / numerical approximation

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

`N[...,p]` の `p` は**有効10進桁数**。

```text
accuracy[N[1/3,20]]
precision[N[1/3,20]]
rationalize[N[1/3,20]]   -> 1/3
explain[N[Pi,20]]
```

`:fix n` は計算値ではなく**表示上の小数点以下最大桁数**だけを変える。

## 5. 変数・ユーザー函数・履歴

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

履歴:

| 入力 | 意味 |
|---|---|
| `%`, `%%`, ... | 直前，2つ前，…の成功出力 |
| `@`, `@@`, ... | 直前，2つ前，…の入力を現在環境で再評価 |
| `Out[n]` | 保存済み出力snapshot |
| `In [n]` | 過去入力を現在環境で再評価 |

```text
N[@,30]
Out[-1]
In [-1]
```

## 6. 条件分岐・assumption

評価制御:

```text
if[x>0,1,-1]
```

数学的piecewise:

```text
cases[1/x if x!=0; 0 if x==0]
```

assumption:

```text
element[x,Real]
element[n,Integer]
simplify[sqrt[x^2],element[x,Real]]   -> abs[x]
simplify[sqrt[x^2],x>=0]             -> x
```

## 7. 基本数学

| 分類 | 主な入力 |
|---|---|
| 根・絶対値 | `sqrt[x]` `cbrt[x]` `abs[x]` `sign[x]` |
| 複素数 | `re[z]` `im[z]` `conj[z]` `arg[z]` |
| 指数・対数 | `exp[x]` `log[x]` `log[b,x]` `log2[x]` `log10[x]` |
| 安定形 | `expm1[x]` `log1p[x]` `sinc[x]` `cosc[x]` `tanc[x]` |
| 三角 | `sin[x]` `cos[x]` `tan[x]` `cot[x]` `sec[x]` `csc[x]` |
| 逆三角 | `asin[x]` `acos[x]` `atan[x]` `atan2[y,x]` |
| 双曲線 | `sinh[x]` `cosh[x]` `tanh[x]` `csch[x]` `sech[x]` `coth[x]` |
| 逆双曲線 | `asinh[x]` `acosh[x]` `atanh[x]` |

```text
sqrt[8]       -> 2sqrt[2]
abs[3+4I]     -> 5
re[3+4I]      -> 3
im[3+4I]      -> 4
conj[3+4I]    -> 3-4I
arg[-1]       -> Pi Rad
log[10,1000]  -> 3
```

## 8. 丸め・整数・数論

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

## 9. 式変形

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

通常式へ戻す:

```text
normal[series[(1+x)^3,{x,0,5}]]
-> x^3+3x^2+3x+1
```

`normal` はtop-levelだけ。式中の対応済みSeriesも再帰変換するなら:

```text
toNormal[expr]
```

## 11. 微分・積分・極限

```text
D[x^3,x]                    -> 3x^2
D[x^5,{x,3}]
D[x^2*y^3,x,y]
D[exp[x^2],x]               -> 2x exp[x^2]
```

点で数値微分:

```text
diff[x^2,x,3]               -> 6.0
diff[sin[x],x,1,20]
```

積分:

```text
integrate[x^2,x]                    -> x^3/3
integrate[sin[x],x]                 -> -cos[x]
integrate[x^2,{x,0,1}]              -> 1/3
integrate[exp[-x],{x,0,Infinity}]   -> 1
integrate[expr,x,assumptions]
```

保証付き数値積分:

```text
nintegrate[x^2,{x,0,1},12]
```

極限:

```text
limit[sin[x]/x,x,0]   -> 1
limit[1/x,x,0,1]      -> Infinity
limit[1/x,x,0,-1]     -> -Infinity
```

方向は `1` が右，`-1` が左。

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

等式系の既定ambient domainは `Complex`。ordered inequalityは `Real` 系。

有限解のbranch選択は**0-based**。

```text
s:=solve[x^2==1,x]
at[s,0]     -> {x == 1}
at[s,1]     -> {x == -1}
at[s,1,x]   -> -1
```

exact root:

```text
root[{-2,0,1},2]
N[root[{-2,0,1},2],30]
```

`root[{a0,a1,...,an},k]` の係数は昇冪順。

## 13. Array・列生成

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

**`at` は0-based。**

```text
reshape[{1,2,3,4},{2,2}]
zeros[2,3]
identity[3]
transpose[{{1,2},{3,4}}]
```

列生成・map:

```text
range[5]             -> {1,2,3,4,5}
range[0,1,1/3]       -> {0,1/3,2/3,1}
table[i^2,{i,5}]     -> {1,4,9,16,25}
map[sin,{0,Pi/2,Pi}] -> {0,1,0}
```

`sin[A]` や `exp[A]` は自動element-wiseではない。明示的に `map` を使う。

集約:

```text
sum[{1,2,3}]    -> 6
prod[{1,2,3}]   -> 6
min[{1,2,3}]
max[{1,2,3}]
mean[{1,2,4}]   -> 7/3
```

## 14. 行列

**`A*B` は行列積ではない。行列積・contractionは `dot[A,B]`。**

```text
dot[{{1,2},{3,4}},{{5,6},{7,8}}]
-> {{19,22},{43,50}}
```

主な函数:

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

分解結果も `at` で取得する。

```text
lu:=luDecomposition[A]
at[lu,0]
at[lu,1]
at[lu,2]
```

## 15. Vector

rank-1 Arrayをvectorとして使う。

```text
norm[{3,4}]       -> 5
normalize[{3,4}]  -> {3/5,4/5}
cross[{1,0,0},{0,1,0}]
distance[a,b]
outer[a,b]
projection[a,b]
rejection[a,b]
```

複素vectorでは `dot` と `inner` を区別する。

```text
inner[{1,I},{1,I}] -> 2
dot[{1,I},{1,I}]   -> 0
```

- `dot[a,b]`: bilinear，共役なし。
- `inner[a,b]`: Hermitian内積，第1引数を共役。

直交化:

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

## 16. 統計・Signal・乱数

統計:

```text
mean[{1,2,4}]                  -> 7/3
median[1,2,3,4]                -> 5/2
quantile[1/4,1,2,3,4,5,6,7]   -> 5/2
var[1,2,3]                     -> 2/3
stddev[1,2,3]                  -> sqrt[6]/3
corr[{1,2,3},{2,4,6}]          -> 1
```

Signal:

```text
dft[v]
fft[v]
ifft[v]
convolve[a,b]
ifft[fft[{1,2,3,4}]] -> {1,2,3,4}
```

Fourier位相はsessionの角度modeに依存せずRadian。

乱数:

```text
randSeed[42]
rand[]
randint[1,6]
choice[{a,b,c}]
randn[]
```

## 17. 特殊函数

代表的な入力形。枝・定義域等は `:help <函数名>`。

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

`ellipticF/E/Pi` の振幅 `phi` はsession設定にかかわらずRadianとして扱う。

## 18. CLI早見

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

## 19. 間違えやすい点

| 誤り | 正しい入力 |
|---|---|
| `sin(x)` | `sin[x]` |
| `[x+1]` でgrouping | `(x+1)` |
| `sin[30]` を30度と思う | `sin[30 Deg]` |
| `:angle deg` | `angleMode[Deg]` |
| `A*B` を行列積と思う | `dot[A,B]` |
| 最初の要素を `at[A,1]` | `at[A,0]` |
| `sin[A]` をelement-wiseと思う | `map[sin,A]` |

迷ったら:

```text
:help <name>
:help functions
:status
```
