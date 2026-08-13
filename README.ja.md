# mmCalculator – Mathematical Machinery Calculator

© 2021–2026 mmKreutzef (aka Daiki.NIIMI)  
Licensed under the BSD 3-Clause License

**Current release: v1.5.1**

[English](README.md) | [日本語](README.ja.md)

## 概要

mmCalculator（以下mmCal）は，研究・設計・製造などの技術用途を想定した，厳密計算を優先する(exact-first)CLI数式計算機／小型CASである。

一般的な電卓として使える一方で，次の機能を備える。

- **整数・分数・記号式を，可能な限り厳密(exact)な形のまま計算**
- 過去の計算結果を用いた連続計算
- 変数およびユーザー定義函数
- 複素数，ベクトル，行列
- 展開・因数分解・簡約・微分・積分・極限・方程式求解
- 任意精度の数値近似，特殊函数，統計，信号処理

鈍重なシステムとは違い，本ツールは：

> **軽量、そしてすぐに使える — 実務に役立つ実用的なツール**

また一般的な電卓のように入力直後から`double`へ変換せず，必要なときだけ数値近似を求める。


## 細かいことは置いといて実例超特急

```text
In[1]> 999999999999999999999999999999^2
Out[1]> 999999999999999999999999999998000000000000000000000000000001

In[2]> 0.1+0.2
Out[2]> 3/10

In[3]>  0.1+0.2==0.3
Out[3]> True

In[4]> 1/3+1/6
Out[4]> 1/2

In[5]> sqrt[72]
Out[5]> 6sqrt[2]

In[6]> sin[Pi/6]
Out[6]> 1/2

In[7]> expand[(x+1)^3]
Out[7]> x^3+3x^2+3x+1

In[8]> factor[x^2-1]
Out[8]> (x-1)(x+1)

In[9]> fullSimplify[(x^2-1)/(x-1),x!=1]
Out[9]> 1+x

In[10]> simplify[sqrt[x^2],element[x,Real]]
Out[10]> abs[x]

In[11]> D[exp[x^2],x]
Out[11]> 2x exp[x^2]

In[12]> integrate[x^2+sin[x],x]
Out[12]> x^3/3-cos[x]

In[13]> integrate[sin[x],{x,0,Pi}]
Out[13]> 2

In[14]> limit[(1-cos[x])/x^2,x,0]
Out[14]> 1/2

In[15]> solve[x^2+1==0,x,Complex]
Out[15]> {x==I, x==-I}

In[16]> (1+I)/(1-I)
Out[16]> I

In[17]> sqrt[-8]
Out[17]> 2I sqrt[2]

In[18]> Pi
Out[18]> Pi

In[19]> N[%,30]
Out[19]> 3.141592653589793238462643383280

In[20]> N[%%%,20]
Out[20]> 2.82842712474619009760I
```

以下では概要のみを示す。
各函数の仕様や内部の詳細は[リファレンス](docs/reference.ja.md)または`docs`フォルダ内の文書を参照。

## v1.5.1

v1.5.1は，v1.5.0のexact-first CAS基盤を保ったまま，canonicalization，検証基盤，多倍長整数，高精度数値評価を重点的に固めたreleaseである。

主な内部改善:

- BigInt乗算をschoolbook / Karatsuba / Toom-3の適応dispatchへ変更し，専用squareも追加
- Burnikel–Ziegler除算，2冪除算fast path，divide-and-conquer 10進変換
- `Pi`をbinary-splitting Chudnovsky，`exp` / `log`をbinary splitting中心のcertified算法へ更新
- 巨大Radianの三角函数を保証付きargument reductionで高速化
- BigFloatの極端なexponent gap加減算をdirected roundingを保ったままfast-path化
- formatter/parser生成型round-trip，積分derivative-back，Reference↔registry照合，極端値近似testを追加
- `mmCal.Benchmarks`を独立projectとして追加

高速化の採用理由，棄却した算法，代表benchmarkは[`docs/performance_optimization.ja.md`](docs/performance_optimization.ja.md)を参照。

## 1. まず使う

Windowsでは`mmCal.exe`を起動するだけ。

起動時に表示桁数や既定角度を指定できる。

```text
mmCal --fix 16 --angle deg
mmCal --angle rad
mmCal --angle grad --fix 8
```

- `--fix 16`: 結果を小数点以下**最大16桁**で表示する。末尾の不要な0は省略
- `--angle deg`: 角度指定のない三角函数を度として扱う
- `--angle rad`: ラジアン。既定値
- `--angle grad`: グラード
- `--help`, `-h`: 起動オプションを表示する

Linux等ではCMake 3.20以上とGCCまたはClangを用いてソースからビルドできる。

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## 2. 「正確な値」と「小数表示」は別物

`0.1`も最初から二進浮動小数へ変換せず，exactな`1/10`として扱う。
巨大整数・有理数・根号を含む代数的表現・複素数・記号式も可能な限りexactな形を保ち，展開・因数分解・簡約・微分・積分・極限・方程式求解などを同じ式体系上で処理する。

```text
In[1]> 1/3
Out[1]> 1/3

In[2]> sqrt[2]
Out[2]> sqrt[2]
```

`Pi`，`E`，`Phi`，`sqrt[2]`なども，保存済みの機械精度浮動小数点数として扱うのではなく，exactな記号式と数値近似を区別する。

小数値が必要な場合だけ`N[expr,n]`で任意精度の数値近似を要求する。

```text
In[3]> N[1/3,20]
Out[3]> 0.33333333333333333333

In[4]> N[Pi,30]
Out[4]> 3.141592653589793238462643383280
```

`N[式,桁数]`は式そのものを数値近似した値を返す。
一方，`:fix`は**画面上の見え方だけ**を変更し，保存されるexact値は変更しない。

```text
:fix 6
Display: Fixed(6)

In[5]> 1/3
Out[5]> 0.333333

:fix off
Display: Exact

In[6]> Out[5]
Out[6]> 1/3
```

つまり`:fix`を使っても，保存されている`Out[5]`は`1/3`のままである。
計算値と表示形式は分離している。

### precision / accuracy / rationalize

`N`で得た近似値は，単なる表示文字列だけでなく，真値を必ず含む保証区間(certified enclosure)などの精度情報を保持する。

- `accuracy[x]`: 真値に対して保証できる**絶対10進桁数**の整数下限
- `precision[x]`: 真値に対して保証できる**相対10進桁数**の整数下限
- `rationalize[x]`: 近似値の保証区間からexactな有理数を復元する

```text
In[7]> accuracy[N[1/3,20]]
Out[7]> 20

In[8]> precision[N[1/3,20]]
Out[8]> 19

In[9]> rationalize[N[1/3,20]]
Out[9]> 1/3

In[10]> accuracy[1/3]
Out[10]> Infinity
```

`precision`や`accuracy`は，`N[...,n]`の`n`をそのまま返す函数ではない。
保証可能な誤差境界から保守的に求めるため，値によっては要求桁数より小さくなる。
exactな値やexactな記号式は，この意味では`Infinity`を返す。

数値近似では，単に計算を何回か繰り返して「たぶん合っている」とはしない。
真値を含む範囲を計算し，必要な10進表示が一意に確定したことを確認する方式を基本としている。

## 3. 基本構文

函数呼び出しは**角括弧 `[]` のみ**を使用する。丸括弧 `()` は数式のグルーピング専用であり、函数呼び出しには使用しない。

```text
sin[Pi/6]
sqrt[2]
log[10,1000]
```

したがって `sin(Pi/6)` は函数呼び出しではない。既知の函数名に旧 `()` 構文を使った場合はSyntaxErrorとし、`sin[Pi/6]` のように書く。`x(x+1)` のような通常identifierと丸括弧の隣接は暗黙乗算として受理し、Formatterは明確さのため `x*(x+1)` と正規化する。

四則演算や冪は通常の記法を使う。

```text
2+3*4
(x+1)^2
1/(x+1)
2Pi
2sqrt[2]
```

有限小数literalはexactな分数として読み取る。

```text
0.125
-> 1/8
```

配列は波括弧。

```text
{1,2,3}
{{1,2},{3,4}}
```

## 4. 角度

既定はラジアンである。

```text
sin[Pi/6]
-> 1/2
```

セッション中の既定角度は`angleMode`で変更できる。

```text
angleMode[]
-> Rad

angleMode[Deg]
-> Deg

sin[30]
-> 1/2
```

式の一部だけ明示指定することもできる。

```text
sin[30 Deg]
sin[Pi/6 Rad]
```

明示した`Deg` / `Rad` / `Grad`はセッション設定より優先される。

## 5. 変数・ユーザー函数

```text
x:=2
-> 2

f(t):=t^2+1
f(4)
-> 17
```

既存定義を変更すると通知される。

```text
x:=4
INFO: x redefined (was 2)
-> 4
```

現在の定義一覧:

```text
Defs[]
```

定義の削除:

```text
UnDef[x]
UnDef[x,y,f]
```

定義と履歴をまとめて削除:

```text
Clear[]
```

終了:

```text
Exit[]
```

## 6. 入出力履歴

直前の成功結果は`%`，直前の入力式は`@`で参照できる。`%%`，`%%%`はさらに前の成功出力，`@@`，`@@@`はさらに前の入力式を参照する。連続個数に固定上限はない。

```text
In[1]> 2+3
Out[1]> 5

In[2]> %*2
Out[2]> 10

In[3]> Pi
Out[3]> Pi

In[4]> N[@,30]
Out[4]> 3.141592653589793238462643383280
```

正式な履歴参照は`In[n]` / `Out[n]`で，正の添字は絶対番号，負の添字は相対参照とする。`0`は無効。

```text
In[1]
Out[1]
In[-1]    // 直前の入力。@ と同義
Out[-1]   // 直前の成功出力。% と同義
```

`Out[n]`は保存済み出力snapshotを返し，再評価しない。負の`Out[-n]`は成功出力だけを数えるため，評価失敗を挟んでも`% == Out[-1]`が成立する。
`In[n]`は過去のlowered入力式を取得し，**現在の定義環境でもう一度評価する**。負の`In[-n]`は入力slotを数え，`In[-1]` / `@`は直前の入力を再評価する。

## 7. 主な数学機能

### 基本・複素数

```text
abs[-3]              -> 3
abs[3+4I]            -> 5
re[3+4I]             -> 3
im[3+4I]             -> 4
conj[3+4I]           -> 3-4I
arg[-1]              -> Pi Rad
```

`sqrt`，`log`，一般の冪，逆三角函数，逆双曲線函数などは，複素数まで含めた主値を基準にする。

```text
sqrt[-1]             -> I
log[-1]              -> I Pi
```

### 三角・双曲線

`sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `asin`, `acos`, `atan`, `atan2`，および`sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`などを実装している。

### 記号微分

```text
D[x^3+sin[x],x]
-> 3x^2+cos[x]

D[x^5,{x,3}]
-> 60x^2
```

### 記号積分

不定積分の`+C`は表示しない。原始函数の代表を1つ返す。

```text
integrate[x^2+sin[x],x]
-> x^3/3-cos[x]

integrate[x cos[x],x]
-> x sin[x]+cos[x]

integrate[1/(1+x^2),{x,0,1}]
-> Pi/4

integrate[exp[-x],{x,0,Infinity}]
-> 1
```

### 極限

```text
limit[sin[x]/x,x,0]       -> 1
limit[1/x,x,0,1]          -> Infinity
limit[1/x,x,0,-1]         -> -Infinity
```

### 方程式・不等式

```text
solve[x^2-2==0,x]
-> {x==sqrt[2], x==-sqrt[2]}

solve[x^2<4,x,Real]
-> {x in Real if x>-2&&x<2}

solve[exp[x]==2,x,Real]
-> {x==log[2]}
```

完全な解集合を保証できない場合，都合のよい1解だけを返さない。
未解決であることをWarningと結果で示す。

## 8. Array・行列・ベクトル・統計

`{...}`は一般の有限brace containerであり，vector / matrix / tensor専用の構文ではない。child shapeが全て一致する矩形値は**shape + row-major flat storage**のdense `ArrayExpr`へ自動最適化し，`{Q,R}`のようにshapeが異なる値は一般braceとして保持する。Matrix函数へ入る境界では矩形性を監査し，非矩形値はWarning + 未評価とする。shapeをbraceだけでは保存できない空Arrayのみ`reshape`を用いて表示する。

```text
dimensions[{{1,2,3},{4,5,6}}]
-> {2, 3}

arrayRank[{{1,2},{3,4}}]
-> 2

at[{{1,2},{3,4}},1]
-> {3, 4}

at[{{1,2},{3,4}},1,0]
-> 3

reshape[{1,2,3,4},{2,2}]
-> {{1, 2}, {3, 4}}

zeros[0,3]
-> reshape[{}, {0, 3}]
```

添字は0-based。`arrayRank[A]`はArrayの次元数，`matrixRank[A]`は線形代数上の階数であり，意味を分離している。

基本線形代数のcanonical APIは次のとおり。

```text
dot[{{1,2},{3,4}},{{5,6},{7,8}}]
-> {{19, 22}, {43, 50}}

dot[{1,2},{3,4}]
-> 11

det[{{1,2},{3,4}}]
-> -2

inverse[{{1,2},{3,4}}]
-> {{-2, 1}, {3/2, -1/2}}

rref[{{1,2},{3,4}}]
-> {{1, 0}, {0, 1}}

matrixRank[{{1,2},{2,4}}]
-> 1

nullSpace[{{1,2},{2,4}}]
-> {{-2, 1}}

solveLinear[{{2,1},{1,-1}},{5,1}]
-> {2, 1}

luDecomposition[{{0,2},{3,4}}]
-> {{{0,1},{1,0}},{{1,0},{0,1}},{{3,4},{0,2}}}

qrDecomposition[{{3,0},{4,0}}]
-> {{{-3/5,-4/5},{-4/5,3/5}},{{-5,0},{0,0}}}

svd[{{3,0},{0,4}}]
-> {{{0,1},{1,0}},{{4,0},{0,3}},{{0,1},{1,0}}}

eigenvalues[{{0,-1},{1,0}}]
-> {I, -I}

norm[{3+4I}]
-> 5

normalize[{3,4}]
-> {3/5, 4/5}
```

`dot`はvector-vector / matrix-vector / vector-matrix / matrix-matrixを扱う。
`A*B`を行列積にはせず，同shape Arrayの`+` / `-`とscalar×Arrayのみを通常算術へ統合する。行列積・内積は`dot[A,B]`で明示する。

整数・Rational行列は不用意にBigFloatへ変換せずexactに処理する。`det` / `rref` / `matrixRank` / `nullSpace` / `inverse` / `solveLinear`は行ごとの分母除去とBareiss fraction-free eliminationを共有し，中間Rationalの増殖を抑える。exact complex行列は`Number` Gaussian backendへfallbackする。`solveLinear[A,b]`は一意解だけを返し，整合した過剰決定系もfull column rankなら扱う。不整合系や自由変数が残る系はDomain errorとし，parametric solutionを捏造しない。一般symbolic determinant / inverseには式爆発を防ぐ展開budgetを設け，三角・疎行列など安全に処理できる場合を除き，巨大な式を作る前に未評価で保持する。`luDecomposition[A]`は正方行列に`{P,L,U}`を返す。`qrDecomposition[A]`は矩形にも対応するreduced Householder QRで，m×nに対し`k=min[m,n]`，`Q:m×k`，`R:k×n`の`{Q,R}`を返す。`svd[A]`も矩形reduced `{U,S,V}`を返し，一般数値backendは条件数を二乗する`A^H A`を形成せずHouseholder bidiagonalization + one-sided Jacobiを使う。`at[result,0]`等でfactorを取り出せる。一般exact QRは式爆発を避けるため3×3以下に制限し，上三角/上台形caseだけ任意次数のexact fast pathを許す。`eigenvalues` / `eigenvectors` / `eigensystem`は正方行列を対象とし，exact pathは三角・対角・distinct-root exact Number 2×2を処理，一般`N[...]`はHessenberg + implicit shifted complex QRでSchur形を作り，元入力intervalに対するSchur/eigenpair relationを監査する。defective/近接重根で独立vectorを安全に構成できない場合は推測しない。

`N`の下ではFFTと同様にprecision-aware backendへ直接dispatchする。

```text
N[det[{{Pi,0},{0,2}}],12]
-> 6.283185307180

N[inverse[{{Pi,0},{0,2}}],12]
-> {{0.318309886184, 0}, {0, 0.5}}

N[solveLinear[{{Pi,0},{0,2}},{Pi,4}],12]
-> {1.0, 2}

N[qrDecomposition[{{1,2},{3,4}}],8]
-> {{{-0.31622777,-0.94868330},{-0.94868330,0.31622777}},{{-3.16227766,-4.42718872},{0,-0.63245553}}}

N[eigenvalues[{{1,2},{3,4}}],8]
-> {-0.37228132, 5.37228132}
```

この経路はexactな巨大中間式を完成させてから丸めるのではなく，BigFloat/interval系のcertified演算で要求精度を直接処理する。`solveLinear`もaugmented interval eliminationを直接試し，pivotと整合性を証明できない場合にepsilonで推測しない。特に依存した過剰決定系ではinterval相関のため整合性証明が難しい場合がある。`matrixRank`と`nullSpace`はrank deficiencyに依存する不連続量なのでexact入力ではexact eliminationを先に使い，近似入力・未解決caseでもepsilon閾値は導入しない。interval backendはpivot構造を証明できる場合だけ結果を返し，rank deficiencyを推測しない。

`N`の表示ではexact有限小数は `N[1/2,10] -> 0.5` のように不要な0埋めをしない。certified interval由来で固定桁の末尾0が連続する場合は1個だけ残し，`1.000000000000 -> 1.0`，`1.500000000000 -> 1.50` と表示する。要求桁数と保証区間はmetadataに保持される。

旧`matmul` / `mmul` / `vdot`，`rank` / `mrank`，`vnorm`，`vnormalize`，`mget`は互換aliasとして残している。

統計処理でも，可能な範囲では要素をexactな数や記号式のまま扱う。

```text
mean[{1,2,4}]
-> 7/3

var[{1,2,3}]
-> 2/3

stddev[{1,2,3}]
-> sqrt[6]/3

corr[{1,2,3},{2,4,6}]
-> 1
```

## 9. FFT・乱数

```text
dft[{1,2,3,4}]
fft[{1,2,3,4}]

ifft[fft[{1,2,3,4}]]
-> {1, 2, 3, 4}

convolve[{1,2},{3,4}]
-> {3, 10, 8}
```

乱数はセッションごとの状態を保つ。

```text
randSeed[42]
rand[]
randint[1,6]
choice[{a,b,c}]
randn[]
```

同じseedを設定すれば同じ列を再現する。
暗号用途の乱数ではない。

## 10. Warningについて

`D`，`integrate`，`limit`，`solve`などで完全な評価を証明できない場合，誤った値を作らず，Warningとともに式を未評価のまま返すことがある。

```text
integrate[gamma[x],x]
WARN: integrate could not fully prove the symbolic antiderivative or definite integral; unevaluated integrate[...] remains
-> integrate[gamma[x], x]
```

これは「計算結果がその式」という意味ではなく，**現在の実装では安全に計算し切れなかった**という通知である。
未解決領域は今後も拡充予定。がんばります。

## 11. CLI表示設定

```text
:fix 16
:fix off
:status
```

`:fix n`は小数点以下最大`n`桁へ丸めて**表示するだけ**で，保存されている値や`precision` / `accuracy`の意味は変更しない。
`:status`では現在の角度，表示形式，定義数，履歴数などを確認できる。
コンソールタイトルにも角度と表示形式を補助表示する。

## 12. 名前について

通常の数学函数は小文字を標準とする。

```text
sin cos log sqrt integrate solve
```

一方，短い記号演算やセッション操作には次の固有名を使う。

```text
D N In Out Exit Clear Defs UnDef
```

## 13. 詳細資料

- `docs/reference.ja.md` — 函数・構文・現在仕様の詳細
- `docs/mathematics.md` — 定義域，主値，数値計算の数学方針
- `docs/architecture.md` — 開発者向け内部構造
- `docs/roadmap.md` — 現在未実装の主な機能と今後の候補
- `docs/grammar.ebnf` — 文法の機械可読な概要
- `docs/multiprecision_implementation.ja.md` — 多倍長整数・任意精度・保証付き評価の実装詳細
- `docs/performance_optimization.ja.md` — v1.5.1で採用・棄却した高速化と実測根拠
- `CHANGELOG.ja.md` — releaseごとの主要変更

## 14. ライセンス

BSD 3-Clause

Copyright (c) 2021–2026 mmKreutzef

詳細な条件は`LICENSE`を参照。

学術論文や製品等でmmCalを利用した場合，**ライセンス上の追加義務ではないが**，ドキュメントや出版物等で使用した旨を記載していただけると嬉しい。
お知らせいただければ，さらに喜びます。

また，BSD 3-Clause上は許可されるが，ソースへ別プラットフォームやUIを被せただけの再販は作者が少し萎えます。

何かのソフトウェアの一部として組み込む利用は大歓迎です。

## 15. テスト・制作環境

v1.5.1時点で，本プロジェクトには1691件の内部回帰テストと1337件のブラックボックステストが含まれる。
exact算術，境界値，定義域，エラー分類，formatterの再入力性，数値近似の保証区間などを重点的に検証している。さらに`mmCal.Benchmarks`を独立projectとして用意し，固定seedのランダム正当性試験，算法threshold sweep，巨大数・高精度函数の性能比較を通常testから分離して実行できる。

主なWindows開発環境:

- Windows 10 1909
- Intel Core i7-9800X
- Visual Studio Community 2026 / MSVC

LinuxではGCCおよびClangによるビルド・テストも行っている。

ドキュメント整理，実装方針の検討，テスト設計，および積分規則などのMathKnowledge整理にはLLMを補助的に利用している。

## 16. 備考

本プロジェクトは，厳密な演算子意味論と実用的な数式評価の両立を目指している。
研究・設計・製造現場などで，軽量なCLIとして簡便に利用できることも目的としている。

そして，私が欲しいものを作っているだけである。

## 17. 免責事項

本ソフトウェアはBSD 3-Clause Licenseに定めるとおり，**現状のまま（AS IS）**提供される。
商業的利用の適合性，特定目的への適合性，非侵害など，明示または黙示の保証は一切ありません。
著者または著作権者は，契約，不法行為，その他の理由で発生する，または発生したソフトウェアに関連するすべての請求，損害，またはその他の責任に対して，一切責任を負いません。

本ソフトウェアを使用することにより，ソフトウェアの使用に関して発生するすべてのリスクは自己責任であることを認め，自動的に同意したことになります。著者は，データの損失，システムの不具合，その他ソフトウェアの使用によって生じた損害について一切責任を負いません。

正式な条件および免責事項は`LICENSE`を参照。

## 謝辞

このツール開発のきっかけとなった過去の自分，学びの園であった大学と教授殿，そして実使用環境となっている現職場に感謝の意を表します。

## 要望等

要望・不具合報告・実装提案等はGitHubへ投げてください。大歓迎です🍀

## 将来のお話

- `for`的なものは欲しいよね
- `plot`函数（グラフ描画）
