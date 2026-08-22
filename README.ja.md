# mmCalculator – Mathematical Machinery Calculator

© 2021–2026 mmKreutzef (aka Daiki.NIIMI)  
Licensed under the BSD 3-Clause License

**Current release: v1.5.3**

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

> **軽量，そしてすぐに使える — 実務に役立つ実用的なツール**

また一般的な電卓のように入力直後から`double`へ変換せず，必要なときだけ数値近似を求める。

## 細かいことは置いといて実例超特急

```text
In [1]> 999999999999999999999999999999^2
Out[1]> 999999999999999999999999999998000000000000000000000000000001

In [2]> 0.1+0.2==0.3
Out[2]> True

In [3]> 1/3+1/6
Out[3]> 1/2

In [4]> sqrt[72]+sin[Pi/6]
Out[4]> 1/2+6sqrt[2]

In [5]> factor[expand[(x+1)^3]-1]
Out[5]> x*(x^2+3x+3)

In [6]> fullSimplify[(x^2-1)/(x-1),x!=1]
Out[6]> 1+x

In [7]> simplify[sqrt[x^2],element[x,Real]]
Out[7]> abs[x]

In [8]> D[exp[x^2],x]
Out[8]> 2x exp[x^2]

In [9]> integrate[x^2+sin[x],{x,0,Pi}]
Out[9]> 2+Pi^3/3

In [10]> limit[(1-cos[x])/x^2,x,0]
Out[10]> 1/2

In [11]> root[{-2,0,1},2]==sqrt[2]
Out[11]> True

In [12]> element[sqrt[2]+sqrt[3],Rational]
Out[12]> False

In [13]> solve[1.1^x==x^2,x,Real]
Out[13]> {x == -2lambertw[log[11/10]/2]/log[11/10], x == -2lambertw[-log[11/10]/2]/log[11/10], x == -2lambertw[-1, -log[11/10]/2]/log[11/10]}

In [14]> N[%,20]
Out[14]> {x == -0.95548727594562198165, x == 1.0513800237472769374, x == 95.71683016840522274}

In [15]> solve[x^2+1==0,x,Complex]
Out[15]> {x == I, x == -I}

In [16]> A:={{1,2},{3,4}}
Out[16]> {{1,2},{3,4}}

In [17]> det[dot[A,inverse[A]]]
Out[17]> 1

In [18]> ifft[fft[{1,2,3,4,5,6,7}]]
Out[18]> {1,2,3,4,5,6,7}

In [19]> Pi
Out[19]> Pi

In [20]> N[@+x,30]
Out[20]> 3.14159265358979323846264338328+x

In [21]> explain[N[@@,5]]
Out[21]> {{"Kind", "DecimalApproximation"}, {"Domain", "Real"}, {"Exactness", "CertifiedApproximation"}, {"RequestedPrecisionDigits", 5}, {"DisplayedFractionalDigits", 4}, {"Rounded", True}, {"CertifiedEnclosure", {6908435304715/2199023255552, 13816870609431/4398046511104}}, {"InformationEnclosure", {62831/20000, 62833/20000}}}
```

以下では概要のみを示す。
各函数の仕様や内部の詳細は[リファレンス](docs/reference.ja.md)または`docs`フォルダ内の文書を参照。

## v1.5.3

v1.5.3は，v1.5.2までのexact-first CAS基盤を保ったまま，**代数数をpersistent number fieldへ進め，Solver・domain knowledge・certified `N`・exact FFTまで表現横断で接続したrelease**である。

主な変更:

- Real/Complex `root[...]`，minimal polynomial，primitive element，persistent `NumberFieldContext` / `AlgebraicElement`，exact equality / Real orderingを統合
- `root[...]`と`sqrt` / `cbrt` / `Phi`等を内部exact algebraic viewで接続し，表示形を変えず比較・domain・Solveの知識を共有
- `solve[equation,Real]`等のdomain短縮形，solve-safe normalization，Lambert Wによる証明付き指数方程式解法，`N[solve[...]]`を追加
- certified近似値をfirst-class数値として扱い，CertifiedEnclosure / InformationEnclosureを分離。`N[x+Pi,p]`のようなstructural partial numericalizationにも対応
- `zeta` / `digamma` / `trigamma` / `ibeta`を追加し，Gamma/Beta高精度backendを大幅に高速化
- first-class `cases[...]`，general `Q[x1,...,xn]` Gröbner basis（Lex / GrLex / GrevLex），`polynomialReduce`，zero-dimensional polynomial Solve統合を追加
- `digamma` / `trigamma`のcertified `N`をcomplex入力まで拡張
- exact非2冪FFTへCyclotomic quotient backendを追加し，5/7/10/12点等のround-tripを巨大なroot-of-unity式なしでexactに閉じる
- `Expr::Node` typed-node化とimmutable paged packed Array + stride viewにより，大規模Array/Matrixのmemory固定費とtranspose costを削減
- `isprime` / `factorint`等の軽量数論，bit utility，`range` / `table` / `map` / `explain`を追加
- semantic fuzzer，black-box timing，algebraic / special-function / cyclotomic FFT benchmarkを常設し，採用しなかった最適化も実測理由付きで記録

正式化時点の検証値はinternal **2335 / 2335**，black-box **1715 / 1715**，Random Expression Fuzzer **100000 / 100000 PASS（8 threads）**。

詳細な変更履歴は[`CHANGELOG.ja.md`](CHANGELOG.ja.md)を参照のこと。

## 1. まず使う

Windowsでは`mmCal.exe`を起動するだけ。

起動時に表示桁数や既定角度を指定できる。

```text
mmCal --fix 16 --angle deg
mmCal --angle rad
mmCal --angle grad --fix 8
mmCal --eval "factor[x^2-1]"
mmCal --batch < expressions.txt
```

- `--fix 16`: 結果を小数点以下**最大16桁(0..1000)**で表示する。末尾の不要な0は省略
- `--angle deg`: 角度指定のない三角函数を度として扱う
- `--angle rad`: ラジアン。既定値
- `--angle grad`: グラード
- `--eval expr`: 1式だけ評価し，値だけを標準出力へ出す
- `--batch`（互換alias: `--bach`）: 標準入力を1行1式として同一sessionで順に評価する。
- `--help`, `-h`: 短い起動usageを表示する。函数の詳細は起動後に`:help sin`等で確認する

`--eval` / `--batch`（`--bach`）は自動処理用であり，banner・prompt・`Out[...]`・終了挨拶を出さない。値は標準出力，Warning / Errorは標準エラーへ分離する。終了codeは成功`0`，引数`2`，Syntax / ResourceLimit`3`，評価`4`，内部error`5`である。batch modeはerror後も次行を処理し，発生した最大の終了codeを返す。

Linux等ではCMake 3.20以上とGCCまたはClangを用いてソースからビルドできる。

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

CMake生成のMSVC buildでは`MMCAL_PARALLEL_COMPILE=ON`が既定であり，compilerへ`/MP`を付与する。無効化する場合はconfigure時に`-DMMCAL_PARALLEL_COMPILE=OFF`を指定する。GCC / Clangではcompiler固有の並列optionを埋め込まず，`cmake --build ... --parallel`でNinja / Make等のbuild toolへ並列性を委ねる。Unity buildは巨大translation unitのmemory消費とincremental rebuild粒度を悪化させるため既定では使用しない。

## 2. 「正確な値」と「小数表示」は別物

`0.1`も最初から二進浮動小数へ変換せず，exactな`1/10`として扱う。
巨大整数・有理数・根号を含む代数的表現・複素数・記号式も可能な限りexactな形を保ち，展開・因数分解・簡約・微分・積分・極限・方程式求解などを同じ式体系上で処理する。

```text
In [1]> 1/3
Out[1]> 1/3

In [2]> sqrt[2]
Out[2]> sqrt[2]
```

`Pi`，`E`，`Phi`，`sqrt[2]`なども，保存済みの機械精度浮動小数点数として扱うのではなく，exactな記号式と数値近似を区別する。

数値近似が必要な場合だけ`N[expr,p]`を使う。`p`は**有効10進桁数(significant decimal digits)**であり，小数点以下の固定表示桁数ではない。固定小数表示は`:fix` / `--fix`が担当する。

```text
In [3]> N[1/3,20]
Out[3]> 0.33333333333333333333

In [4]> N[Pi,30]
Out[4]> 3.14159265358979323846264338328
```

`N[式,p]`は式そのものを`p`有効桁のcertified近似値へ評価する。値のscaleが変わっても相対Precisionを基準にする。branch cutやbackend境界を区間幅だけが跨ぐ場合はguard precisionを増やして再判定するが，top-level `N`は局所16回で必ず停止する。既存の有限precision入力enclosureが原因で要求桁を確定できない場合は`N::precision`，数学値は存在するが現certified backend未対応なら`N::unsupported`として式を保持し，真のDomainErrorとは区別する。
一方，`:fix`は**画面上の小数点以下表示桁数だけ**を変更し，保存されるexact値や`N`のPrecision意味論は変更しない。

```text
:fix 6
Display: Fixed(6)

In [5]> 1/3
Out[5]> 0.333333

:fix off
Display: Exact

In [6]> Out[5]
Out[6]> 1/3
```

つまり`:fix`を使っても，保存されている`Out[5]`は`1/3`のままである。
計算値と表示形式は分離している。

### precision / accuracy / rationalize

`N`で得た近似値は，単なる表示文字列ではなく，真値保証用のCertifiedEnclosureと後続計算で利用してよい情報量を表すInformationEnclosureを保持する。
近似値同士，または近似値とexact Numberの四則演算に加え，interval backendを持つ`sin` / `exp` / `log` / `sqrt`等のscalar函数へも近似値をそのまま再投入できる。
誤差伝播や相殺で保証可能なAccuracy / Precisionは自然に低下し，`N[N[Pi,20],100]`のように外側の`N`を増桁しても，元の近似値に存在しない情報やbackend内部のguard桁を新しいAccuracyとして復元しない。

- `accuracy[x]`: 真値に対して保証できる**絶対10進桁数**の整数下限
- `precision[x]`: 真値に対して保証できる**相対10進桁数**の整数下限
- `rationalize[x]`: 近似値の保証区間からexactな有理数を復元する

```text
In [7]> accuracy[N[1/3,20]]
Out[7]> 20

In [8]> precision[N[1/3,20]]
Out[8]> 19

In [9]> rationalize[N[1/3,20]]
Out[9]> 1/3

In [10]> accuracy[1/3]
Out[10]> Infinity
```

`precision`や`accuracy`は，`N[...,n]`の`n`をそのまま返す函数ではない。 純虚数ではexact-zero成分を誤差源として数えず，`+0` / `*1`等のexact identityや単項符号反転もInformationEnclosureを再量子化しない。zero-centered値はrelative Precisionではなくabsolute Accuracyで情報量を追跡する。
保証可能な誤差境界から保守的に求めるため，値によっては要求桁数より小さくなる。
exactな値やexactな記号式は，この意味では`Infinity`を返す。

数値近似では，単に計算を何回か繰り返して「たぶん合っている」とはしない。
真値を含む範囲を計算し，必要な10進表示が一意に確定したことを確認する方式を基本としている。

## 3. 基本構文

函数呼び出しは**角括弧 `[]` のみ**を使用する。丸括弧 `()` は数式のグルーピング専用であり，函数呼び出しには使用しない。

```text
sin[Pi/6]
sqrt[2]
log[10,1000]
```

したがって `sin(Pi/6)` は函数呼び出しではない。既知の函数名に旧 `()` 構文を使った場合はSyntaxErrorとし，`sin[Pi/6]` のように書く。`x(x+1)` のような通常identifierと丸括弧の隣接は暗黙乗算として受理し，Formatterは明確さのため `x*(x+1)` と正規化する。

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

f[t]:=t^2+1
f[4]
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
In [1]> 2+3
Out[1]> 5

In [2]> %*2
Out[2]> 10

In [3]> Pi
Out[3]> Pi

In [4]> N[@,30]
Out[4]> 3.14159265358979323846264338328
```

正式な履歴参照は`In [n]` / `Out[n]`で，正の添字は絶対番号，負の添字は相対参照とする。`0`は無効。

```text
In [1]
Out[1]
In [-1]    // 直前の入力。@ と同義
Out[-1]   // 直前の成功出力。% と同義
```

`Out[n]`は保存済み出力snapshotを返し，再評価しない。負の`Out[-n]`は成功出力だけを数えるため，評価失敗を挟んでも`% == Out[-1]`が成立する。
`In [n]`は過去のlowered入力式を取得し，**現在の定義環境でもう一度評価する**。負の`In [-n]`は入力slotを数え，`In [-1]` / `@`は直前の入力を再評価する。

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

v1.5.2では，個別規則の羅列だけでなく三角整数冪の有限Fourier還元，負整数冪のsec/csc reduction，異周波数product-to-sum，bounded Weierstrass置換，inverse-chain候補，二次根号等を共通知識として強化した。
積分できない場合も「未実装」「一部のみ解決」「条件不足」「現在の標準函数語彙では既知の有限閉形式なし」を区別してWarningを出し，証明器不足だけを理由に既存のprimitiveを削らない。

### 特殊函数

積分・微分・数値評価を共有するため，次の特殊函数基盤を追加した。

```text
fresnelc[x]  fresnels[x]
hypergeometric1F1[a,b,z]
hypergeometric2F1[a,b,c,z]
ellipticF[phi,m]  ellipticE[phi,m]  ellipticPi[n,phi,m]
Ei[x]  Si[x]  Ci[x]  li[x]  polylog[s,z]
```

v1.5.3ではさらに，`zeta` / `digamma` / `trigamma` / 正則化不完全Beta `ibeta`を追加し，代表exact値・微分関係・当初対応していた実領域のcertified `N`へ接続した。現Unreleased treeでは`digamma` / `trigamma`のcertified評価をcomplex入力まで拡張している。
また軽量数論として`isprime` / `nextprime` / `prevprime` / `factorint` / `totient`を`uint64`範囲で決定的に評価する。
証明backendを超えるBigIntをprobable-primeとして確定しない。

代表例:

```text
integrate[exp[-x^2],x]
-> erf[x]sqrt[Pi]/2

integrate[exp[-x^2],{x,-Infinity,Infinity}]
-> sqrt[Pi]

integrate[exp[x^6],x]
-> x hypergeometric1F1[1/6,7/6,x^6]

integrate[1/(1+x^5),x]
-> x hypergeometric2F1[1,1/5,6/5,-x^5]

integrate[sin[x]/x,x] -> Si[x]
integrate[cos[x]/x,x] -> Ci[x]
integrate[exp[x]/x,x] -> Ei[x]
integrate[1/log[x],x] -> li[x]
integrate[log[1-x]/x,x] -> -polylog[2,x]
```

一般の特殊函数方程式について逆函数を捏造せず，branchや単射性を証明できない場合は`solve`を未解決のまま保持する。

### 極限

```text
limit[sin[x]/x,x,0]       -> 1
limit[1/x,x,0,1]          -> Infinity
limit[1/x,x,0,-1]         -> -Infinity
```

### 方程式・不等式

```text
solve[x^2-2==0,x]
-> {x == sqrt[2], x == -sqrt[2]}

solve[x^2<4,x,Real]
-> {x in Real if x>-2&&x<2}

solve[exp[x]==2,x,Real]
-> {x == log[2]}

solve[sin[x]==0,x,Real]
-> {x == Pi k where k in Integer}
```

`sin/cos/tan`の実軸周期解は，exact非零一次係数を持つaffine argumentから段階的にinteger-parameter familyへ対応している。
非線形argumentやComplex全解をprincipal inverseだけから捏造しない。

一般高次Rational係数多項式のReal解は，既存radical solverで自然に閉じない場合にexact `root`へfallbackできる。
`root[{a0,...,an},k]`は昇冪係数多項式の異なる実根を小さい順に数えた第`k`根で，`N`ではSturm分離区間をcertifiedに細分化する。

```text
solve[x^5-x+1==0,x,Real]
-> {x == root[{1,-1,0,0,0,1},1]}

N[root[{-2,0,1},2],30]
-> 1.41421356237309504880168872421
```

Complex側は`root[{a0,...,an},k,Complex]`で全複素根をcertified isolating diskへ分離し，高次Rational係数多項式のComplex Solveもexact Rootへfallbackできる。
個々のRoot生成時には，次数16以下で証明可能な場合に有理既約因子をexactに選択して定義多項式をminimal polynomialへ縮約する。
Root同士の`+ - * /`では，operandのminimal polynomialが既約と証明でき，`theta=alpha+c beta`が積次数の既約拡大を生成すると証明できる場合にprimitive-element reductionを使い，それ以外はresultantとisolating regionによる再同定へfallbackする。
v1.5.3ではbounded exact algebraic equalityとReal orderingまで実装済みである。
異なるprimitive generator / subfieldを含む一般field merge，`rootApproximant`，次数budgetを超える完全なQ因子分解は意図的に後続へ残す。

完全な解集合を保証できない場合，都合のよい1解だけを返さない。
未解決であることをWarningと結果で示す。

## 8. Array・行列・ベクトル・統計

`{...}`は一般の有限brace containerであり，vector / matrix / tensor専用の構文ではない。
child shapeが全て一致する矩形値はdense `ArrayExpr`へ自動最適化し，v1.5.3ではphysical storageを固定1024要素のimmutable packed page，logical layoutをshape / offset / stridesとして保持する。
transposeや一部のview操作はbackingを共有し，`{Q,R}`のようにshapeが異なる値は一般braceとして保持する。
Matrix函数へ入る境界では矩形性を監査し，非矩形値はWarning + 未評価とする。
shapeをbraceだけでは保存できない空Arrayのみ`reshape`を用いて表示する。

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

有限列の生成と明示的なelement-wise適用には`range` / `table` / `map`を使う。`table`のiteratorはlocal scopeであり，外側の同名定義を汚さない。`map`はscalar leafへ明示適用するため，通常の`exp[A]`等を自動element-wise化しない。

```text
range[0,1,1/3] -> {0,1/3,2/3,1}
table[i^2,{i,5}] -> {1,4,9,16,25}
map[sin,{0,Pi/2,Pi}] -> {0,1,0}
```

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

整数・Rational行列は不用意にBigFloatへ変換せずexactに処理する。
`det` / `rref` / `matrixRank` / `nullSpace` / `inverse` / `solveLinear`は行ごとの分母除去で整数workspaceへliftする。`rref` / `matrixRank` / `nullSpace`と小規模・疎な`det` / `solveLinear`はBareiss fraction-free eliminationを使い，大きいdense exact整数/Rationalの`det`は31-bit prime + CRT，`solveLinear`はCRT + rational reconstructionへ自動dispatchする。modular solveの候補は元の整数系でexact verificationした場合だけ採用し，bad primeや復元失敗時はBareissへfallbackする。modular inverse backendも実装しているが，現GCC benchmarkではautomatic pathはBareissのままとする。
exact complex行列は`Number` Gaussian backendへfallbackする。
`solveLinear[A,b]`は一意解だけを返し，整合した過剰決定系もfull column rankなら扱う。
不整合系や自由変数が残る系はDomain errorとし，parametric solutionを捏造しない。
一般symbolic determinant / inverseには式爆発を防ぐ展開budgetを設け，三角・疎行列など安全に処理できる場合を除き，巨大な式を作る前に未評価で保持する。
`luDecomposition[A]`は正方行列に`{P,L,U}`を返す。`qrDecomposition[A]`は矩形にも対応するreduced Householder QRで，m×nに対し`k=min[m,n]`，`Q:m×k`，`R:k×n`の`{Q,R}`を返す。
`svd[A]`も矩形reduced `{U,S,V}`を返し，一般数値backendは条件数を二乗する`A^H A`を形成せずHouseholder bidiagonalization + one-sided Jacobiを使う。
`at[result,0]`等でfactorを取り出せる。一般exact QRは式爆発を避けるため3×3以下に制限し，上三角/上台形caseだけ任意次数のexact fast pathを許す。
`eigenvalues` / `eigenvectors` / `eigensystem`は正方行列を対象とし，exact pathは三角・対角・distinct-root exact Number 2×2を処理，一般`N[...]`はHessenberg + implicit shifted complex QRでSchur形を作り，元入力intervalに対するSchur/eigenpair relationを監査する。
defective/近接重根で独立vectorを安全に構成できない場合は推測しない。

`N`の下ではFFTと同様にprecision-aware backendへ直接dispatchする。

```text
N[det[{{Pi,0},{0,2}}],12]
-> 6.28318530718

N[inverse[{{Pi,0},{0,2}}],12]
-> {{0.318309886184, 0}, {0, 0.5}}

N[solveLinear[{{Pi,0},{0,2}},{Pi,4}],12]
-> {1, 2}

N[qrDecomposition[{{1,2},{3,4}}],8]
-> {{{-0.31622777,-0.94868330},{-0.94868330,0.31622777}},{{-3.16227766,-4.42718872},{0,-0.63245553}}}

N[eigenvalues[{{1,2},{3,4}}],8]
-> {-0.37228132, 5.3722813}
```

この経路はexactな巨大中間式を完成させてから丸めるのではなく，BigFloat/interval系のcertified演算で要求精度を直接処理する。
`solveLinear`もaugmented interval eliminationを直接試し，pivotと整合性を証明できない場合にepsilonで推測しない。
特に依存した過剰決定系ではinterval相関のため整合性証明が難しい場合がある。
`matrixRank`と`nullSpace`はrank deficiencyに依存する不連続量なのでexact入力ではexact eliminationを先に使い，近似入力・未解決caseでもepsilon閾値は導入しない。
interval backendはpivot構造を証明できる場合だけ結果を返し，rank deficiencyを推測しない。

`N`の表示ではexact有限小数は `N[1/2,10] -> 0.5` のように不要な0埋めをしない。certified interval由来の有効桁結果で末尾0が連続する場合も冗長な0列は圧縮し，要求Precisionと両enclosureはmetadataに保持する。固定小数点以下桁数が必要なら`:fix`を使う。

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

## 11. CLI help / 表示設定

```text
:help
:help sin
:help functions
:help constants
:fix 16
:fix off
:status
```

`:help 函数名`はBuiltinRegistryの函数名・alias・引数個数を正本とし，全callable builtinに個別の説明，明示的な入力規則，一つ以上の例を表示する。`integrate` / `root` / `qrDecomposition` / `svd` / `solve`等の複数形式を持つ函数では，各形式と追加note・複数例も示す。`:help Pi`と`:help constants`は保護された定数，domain，角度単位symbolを扱う。未知名が明確な近傍topicを持つ場合は，`sdv` -> `svd`，`qr` -> `qrDecomposition`のように候補を提示する。
`:help functions`は利用可能なcanonical名とcallable aliasを一覧する。
`:help` / `:fix` / `:status`は評価式ではないため`In[n]`を進めず，履歴にも入らない。
`:fix n`は小数点以下最大`n`桁へ丸めて**表示するだけ**で，保存されている値や`precision` / `accuracy`の意味は変更しない。
`:status`では現在の角度，表示形式，定義数，履歴数などを確認できる。
コンソールタイトルにも角度と表示形式を補助表示する。

Parserはtoken数，AST node数，入れ子深さ，演算子鎖，数値literal桁数，函数引数数，Array要素数を独立に制限する。上限超過は`ResourceLimitError`であり，OSのstack overflowや曖昧な`InternalError`にはしない。

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
- `docs/grammar.ebnf` — 文法の機械可読な概要
- `docs/multiprecision_implementation.ja.md` — 多倍長整数・任意精度・保証付き評価の実装詳細
- `docs/performance_optimization.ja.md` — v1.5.1–v1.5.3で採用・棄却した高速化と実測根拠
- `docs/evaluation_budget.ja.md` — 評価資源上限，cancellation，telemetry，diagnostic契約
- `docs/roadmap.md` — 意図的未実装，探索上限，次の優先候補
- `CHANGELOG.ja.md` — releaseごとの主要変更

## 14. ライセンスと商標

ソースコードは**BSD 3-Clause License**で提供する。商用利用，改変，再配布，組込み利用を含む著作権上の許諾条件は`LICENSE`を参照。

**mmCalの名称・公式ロゴ等のブランド利用は，ソースコードのライセンスとは別に`TRADEMARKS.ja.md`で扱う。**
独立したGUI，fork，商用製品等を作ること自体を制限するものではなく，第三者製品を公式mmCalそのもの・公式認定品であるかのように表示しないための方針である。

- [BSD 3-Clause License](LICENSE)
- [商標・ブランドポリシー](TRADEMARKS.ja.md)

学術論文や製品等でmmCalを利用した場合，ライセンス上の追加義務ではないが，使用した旨を記載していただけると嬉しい。

## 15. テスト・制作環境

v1.5.3正式化時点で，本プロジェクトには2335件の内部回帰テストと1715件のブラックボックステストが含まれる。
exact算術，境界値，定義域，エラー分類，formatterの再入力性，数値近似の保証区間などを重点的に検証している。
さらに`mmCal.Benchmarks`を独立projectとして用意し，固定seedのランダム正当性試験，算法threshold sweep，巨大数・高精度函数の性能比較を通常testから分離して実行できる。

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

本ソフトウェアを使用することにより，ソフトウェアの使用に関して発生するすべてのリスクは自己責任であることを認め，自動的に同意したことになります。
著者は，データの損失，システムの不具合，その他ソフトウェアの使用によって生じた損害について一切責任を負いません。

正式な条件および免責事項は`LICENSE`を参照。

## 謝辞

このツール開発のきっかけとなった過去の自分，学びの園であった大学と教授殿，そして実使用環境となっている現職場に感謝の意を表します。

## 要望等

要望・不具合報告・実装提案等はGitHubへ投げてください。大歓迎です🍀

## 将来のお話

- `for`的なものは欲しいよね
- `plot`函数（グラフ描画）
