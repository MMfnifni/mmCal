# mmCalculator – Mathematical Machinery Calculator

© 2021–2026 mmKreutzef (aka Daiki.NIIMI)  
Licensed under the BSD 3-Clause License

**最新リリース: v1.5.5**

[English](README.md) | [日本語](README.ja.md)

## 概要

mmCalculator（以下mmCal）は，研究・設計・製造などの技術用途を想定した，厳密計算を優先する(exact-first)CLI数式計算機／小型CASである。

一般的な電卓として使える一方で，次の機能を備える。

- **整数・分数・記号式を，可能な限り厳密(exact)な形のまま計算**
- 過去の計算結果を用いた連続計算
- 変数およびユーザー定義函数
- 複素数，ベクトル，行列
- 展開・因数分解・簡約・Taylor/Laurent/Puiseux/対数級数展開・微分・積分・極限・方程式求解
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

In [3]> sqrt[72]+sin[Pi/6]
Out[3]> 1/2+6sqrt[2]

In [4]> factor[expand[(x+1)^3]-1]
Out[4]> x*(x^2+3x+3)

In [5]> fullSimplify[(x^2-1)/(x-1),x!=1]
Out[5]> x+1

In [6]> simplify[sqrt[x^2],element[x,Real]]
Out[6]> abs[x]

In [7]> D[exp[x^2],x]
Out[7]> 2x exp[x^2]

In [8]> integrate[x^2+sin[x],{x,0,Pi}]
Out[8]> 2+Pi^3/3

In [9]> limit[(1-cos[x])/x^2,x,0]
Out[9]> 1/2

In [10]> solve[1.1^x==x^2,x,Real]
Out[10]> {x == -2lambertw[log[11/10]/2]/log[11/10], x == -2lambertw[-log[11/10]/2]/log[11/10], x == -2lambertw[-1, -log[11/10]/2]/log[11/10]}

In [11]> N[%,20]
Out[11]> {x == -0.95548727594562198165, x == 1.0513800237472769374, x == 95.716830168405222740}

In [12]> inverse[{{1,2},{3,4}}]
Out[12]> {{-2, 1}, {3/2, -1/2}}

In [13]> ifft[fft[{1,2,3,4,5,6,7}]]
Out[13]> {1, 2, 3, 4, 5, 6, 7}
```

以下では概要のみを示す。
各函数の仕様や内部の詳細は[リファレンス](docs/reference.ja.md)または`docs`フォルダ内の文書を参照。

## v1.5.5

v1.5.5は，SeriesとVector Calculusを追加しつつ，主眼を**横断的なbug fix，既存frontend間の接続，有限precision意味論，性能の崖，古い固定制限の整理**へ置いた版である。

主な変更点:

- **Series・漸近展開**: `SeriesData` / TPSAを核にTaylor / Laurent / Puiseux / logarithmic Series，`+Infinity`展開，`toNormal`を追加し，`D` / `integrate` / `limit`へ接続
- **Array・Vector**: Hermitian内積を基準にVector APIを整理し，`grad` / `divergence` / `curl` / `laplacian` / `jacobian` / `hessian`等のCartesian Vector Calculusを追加
- **記号計算の接続**: Solve definedness，principal inverse，`cases`境界，Limit binder capture，nested `D` / `integrate` / `series` / `solve`，式変形frontendの組合せで生じる未評価・意味論破壊を横断修正
- **保証付き `N` と表示**: finite-precision formatterの符号・括弧・多項式順を修正し，exact integerの構造parameterを保護。近似値由来を`ExactValue` / `CertifiedInterval` / `VerifiedApproximation`へ分離し，証明用enclosureと残差検証済み数値を区別
- **積分・微分・代数**: 高次有理函数積分をYun分解・polynomial CRT・Hermite reduction・exact residueへ拡張し，高階`D`や多項式×初等函数積分の専用漸化式を強化
- **性能と制限**: `factor[x^257-1]`等のperfect-power空振りを除去し，古い64 / 128 / 256 / 4096等の固定境界を，閉形式・構造fast path・共通`EvaluationBudget`へ置換できる範囲で撤廃
- **CLI・検証**: `:quit` / `:exit` / `:layout`を追加し，certification境界・performance cliff・cross-feature regressionを拡充

詳細な変更履歴は[`CHANGELOG.ja.md`](CHANGELOG.ja.md)の**v1.5.5**を参照。READMEとReferenceはv1.5.5の現在仕様を記述し，旧版固有の変更説明はCHANGELOGへ集約する。

## 1. まず使う

Windowsでは`mmCal.exe`を起動するだけ。

起動時に表示桁数や既定角度を指定できる。

```text
mmCal --fix 16 --angle deg
mmCal --angle rad
mmCal --angle grad --fix 8
mmCal --layout multi
mmCal --eval "factor[x^2-1]"
mmCal --batch < expressions.txt
```

- `--fix 16`: 結果を小数点以下**最大16桁(0..1000)**で表示する。末尾の不要な0は省略
- `--angle deg`: 角度指定のない三角函数を度として扱う
- `--angle rad`: ラジアン。既定値
- `--angle grad`: グラード
- `--layout auto|single|multi`: 通常REPLの表示組版。既定`auto`はTTY上でArray/List/`cases`/解集合を構造的に改行し，pipe/redirect時は1行へ退避する
- `--eval expr`: 1式だけ評価し，値だけを標準出力へ出す
- `--batch`: 標準入力を1行1式として同一sessionで順に評価する。
- `--help`, `-h`: 短い起動usageを表示する。函数の詳細は起動後に`:help sin`等で確認する

`--eval` / `--batch`は自動処理用であり，banner・prompt・`Out[...]`・終了挨拶を出さない。値は標準出力，Warning / Errorは標準エラーへ分離する。終了codeは成功`0`，引数`2`，Syntax / ResourceLimit`3`，評価`4`，内部error`5`である。batch modeはerror後も次行を処理し，発生した最大の終了codeを返す。

Linux等ではCMake 3.20以上とGCCまたはClangを用いてソースからビルドできる。

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

CMake生成のMSVC buildでは`MMCAL_PARALLEL_COMPILE=ON`が既定であり，compilerへ`/MP`を付与する。無効化する場合はconfigure時に`-DMMCAL_PARALLEL_COMPILE=OFF`を指定する。GCC / Clangではcompiler固有の並列optionを埋め込まず，`cmake --build ... --parallel`でNinja / Make等のbuild toolへ並列性を委ねる。Unity buildは巨大translation unitのmemory消費とincremental rebuild粒度を悪化させるため既定では使用しない。

## 2. 「正確な値」と「小数表示」は別物

mmCalはexact-firstを基本とする。整数，有理数，有限小数，根号，複素数，記号式は，可能な限りexactな形を保持する。

```text
1/3      -> 1/3
0.125    -> 1/8
sqrt[2]  -> sqrt[2]
```

小数近似が必要な場合は`N[expr,p]`を使う。`p`は有効10進桁数であり，小数点以下の固定桁数ではない。

```text
N[1/3,20]
-> 0.33333333333333333333

N[Pi,30]
-> 3.14159265358979323846264338328
```

exact整数など，近似後も値が厳密に確定している成分は不要な`.0`を付けずに表示する。`:fix`は保存値を変えず，画面上の小数表示だけを変更する。

数値近似は`ExactValue`，`CertifiedInterval`，`VerifiedApproximation`を区別し，`precision`，`accuracy`，`rationalize`，`explain`で状態を確認できる。詳細なenclosure，provenance，停止条件，branch処理はReferenceへ記載する。

## 3. 基本構文

函数呼び出しは角括弧`[]`を使い，丸括弧`()`はグルーピング専用である。

```text
sin[Pi/6]
log[10,1000]
(x+1)^2
2Pi
```

Arrayや一般のbrace containerは`{...}`を使う。

```text
{1,2,3}
{{1,2},{3,4}}
```

代入・比較・条件式，暗黙乗算，literal，演算子優先順位などの詳細はCheatsheet / Referenceを参照。

## 4. 角度

既定はラジアンである。`angleMode`でセッション既定を変更でき，`Deg` / `Rad` / `Grad`を式中で明示することもできる。

```text
sin[Pi/6]
angleMode[Deg]
sin[30]
sin[30 Deg]
```

角度変換函数: `DtoR`, `DtoG`, `RtoD`, `RtoG`, `GtoD`, `GtoR`

## 5. 変数・ユーザー函数

```text
x:=2
f[t]:=t^2+1
f[4]
-> 17
```

定義・セッション操作: `Defs`, `UnDef`, `Clear`, `Exit`

詳細な評価規則とscopeはReferenceを参照。

## 6. 入出力履歴

直前の成功出力は`%`，直前の入力は`@`で参照できる。`%%` / `@@`のように遡ることもできる。

正式な履歴函数は`In[n]` / `Out[n]`であり，正の添字は絶対番号，負の添字は相対参照である。

```text
Out[-1]
In[-1]
```

`Out`は保存済み出力を返し，`In`は保存済み入力を現在の定義環境で再評価する。詳細はReferenceを参照。

## 7. 主な数学機能

以下は主なcanonical函数名である。aliasを含む完全な一覧は`:help functions`，使用例はCheatsheet，正確な引数形式と定義域はReferenceを参照。

### 数値・複素数

`N`, `precision`, `accuracy`, `rationalize`, `explain`, `sqrt`, `cbrt`, `abs`, `sign`, `re`, `im`, `conj`, `arg`, `cis`, `polar`, `proj`, `hypot`, `fma`, `clamp`

### 初等函数

`exp`, `expm1`, `expc`, `log`, `log1p`, `log2`, `log10`

`sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `asin`, `acos`, `atan`, `atan2`

`sinh`, `cosh`, `tanh`, `csch`, `sech`, `coth`, `asinh`, `acosh`, `atanh`

`sinc`, `cosc`, `tanc`, `sinhc`, `tanhc`

### 整数・離散数学

`floor`, `ceil`, `trunc`, `round`, `frac`, `gcd`, `lcm`, `mod`, `rem`, `quotient`

`bitand`, `bitor`, `bitxor`, `bitnot`, `bitshiftl`, `bitshiftr`, `bitlength`, `bitcount`, `bitget`

`isprime`, `nextprime`, `prevprime`, `factorint`, `totient`, `perm`, `comb`, `fib`, `nextpow2`

### 記号計算・微積分

`simplify`, `fullSimplify`, `expand`, `factor`, `collect`, `cases`

`series`, `normal`, `toNormal`

`D`, `diff`, `integrate`, `nintegrate`, `limit`

`solve`, `root`, `groebnerBasis`, `polynomialReduce`, `element`, `if`

### 特殊函数

`gamma`, `lgamma`, `beta`, `betaln`, `ibeta`, `binom`, `fallingfact`, `risingfact`

`erf`, `erfc`, `zeta`, `digamma`, `trigamma`, `lambertw`

`fresnelc`, `fresnels`, `hypergeometric1F1`, `hypergeometric2F1`

`ellipticF`, `ellipticE`, `ellipticPi`

`Ei`, `Si`, `Ci`, `li`, `polylog`

完全に解けない記号計算では，安全でない推測をせずWarningと未評価式を返す場合がある。

## 8. Array・行列・ベクトル・統計

### Array・列操作

`dimensions`, `arrayRank`, `length`, `at`, `reshape`, `zeros`, `identity`

`range`, `table`, `map`, `sum`, `prod`, `min`, `max`

添字は0-basedである。

### 行列・線形代数

`transpose`, `conjugateTranspose`, `madd`, `dot`, `det`, `inverse`, `rref`, `matrixRank`

`solveLinear`, `nullSpace`, `luDecomposition`, `qrDecomposition`, `svd`

`conditionNumber`, `leastSquares`, `pseudoInverse`

`eigenvalues`, `eigenvectors`, `eigensystem`

`trace`, `rows`, `cols`, `diag`

### ベクトル

`cross`, `norm`, `manhattanDistance`, `distance`, `normalize`

`projection`, `rejection`, `vectorAngle`, `reflectNormal`, `reflectAxis`

`inner`, `outer`, `orthogonalQ`, `orthonormalQ`, `linearIndependentQ`, `gramSchmidt`

### Vector Calculus

`grad`, `divergence`, `curl`, `laplacian`, `jacobian`, `hessian`, `directionalDerivative`

### 統計

`mean`, `median`, `mode`, `quantile`, `percentile`

`var`, `vars`, `stddev`, `stddevs`, `geomean`, `harmmean`, `rms`

`mad`, `madR`, `skew`, `kurtp`, `kurts`, `cv`, `stderr`, `zscore`, `iqr`

`trimmean`, `winsor`, `winsorR`, `cov`, `corr`, `corrspearman`, `percentrank`

Array・線形代数・統計のexact / finite-precision挙動やshape条件はReferenceを参照。

## 9. FFT・乱数

信号処理: `dft`, `fft`, `ifft`, `convolve`

乱数: `randSeed`, `rand`, `randint`, `choice`, `randn`

乱数はセッション状態を持ち，同じseedで同じ列を再現できる。暗号用途ではない。

## 10. Warningについて

`D`, `integrate`, `limit`, `solve`, `N`などは，現在の実装で安全に結果を確定できない場合，誤った値を作らずWarningと未評価式を返すことがある。

Warningは「その式が答え」という意味ではなく，未解決・条件不足・precision不足・backend未対応などを利用者へ通知するためのものである。分類の詳細はReferenceを参照。

## 11. CLI help / 表示設定

```text
:help
:help sin
:help functions
:help constants
:fix 16
:fix off
:layout auto
:layout single
:layout multi
:status
:quit
:exit
```

`:help functions`で現在利用可能なcanonical函数名とcallable aliasを一覧できる。`:help <function>`は入力形式と例を表示する。

`:fix`は表示上の小数桁，`:layout`はREPLの組版だけを変更する。`:status`は現在のsession状態を表示する。

## 12. 名前について

通常の数学函数は小文字をcanonical名とする。

```text
sin cos log sqrt integrate solve
```

短い記号演算・セッション操作には固有名を使う。

```text
D N In Out Exit Clear Defs UnDef
```

aliasを含む現在の名前一覧は`:help functions`を参照。

## 13. 詳細資料

- `docs/reference.ja.md` — 函数・構文・現在仕様の詳細
- `docs/mathematics.md` — 定義域，主値，数値計算の数学方針
- `docs/architecture.md` — 開発者向け内部構造
- `docs/grammar.ebnf` — 文法の機械可読な概要
- `docs/multiprecision_implementation.ja.md` — 多倍長整数・任意精度・保証付き評価の実装詳細
- `docs/performance_optimization.ja.md` — 採用・棄却した高速化と実測根拠
- `docs/evaluation_budget.ja.md` — 評価資源上限，cancellation，telemetry，diagnostic契約
- `CHANGELOG.ja.md` — releaseごとの主要変更

## 14. ライセンスと商標

ソースコードは**BSD 3-Clause License**で提供する。商用利用，改変，再配布，組込み利用を含む著作権上の許諾条件は`LICENSE`を参照。

**mmCalの名称・公式ロゴ等のブランド利用は，ソースコードのライセンスとは別に`TRADEMARKS.ja.md`で扱う。**
独立したGUI，fork，商用製品等を作ること自体を制限するものではなく，第三者製品を公式mmCalそのもの・公式認定品であるかのように表示しないための方針である。

- [BSD 3-Clause License](LICENSE)
- [商標・ブランドポリシー](TRADEMARKS.ja.md)

学術論文や製品等でmmCalを利用した場合，ライセンス上の追加義務ではないが，使用した旨を記載していただけると嬉しい。

## 15. テスト・制作環境

v1.5.5では，内部回帰テスト **3428 / 3428**，ブラックボックステスト **2465 / 2465** の通過を確認している。
exact算術，境界値，定義域，エラー分類，formatterの再入力性，数値近似の保証区間などを重点的に検証している。
さらに`mmCal.Benchmarks`を独立projectとして用意し，固定seedのランダム正当性試験，算法閾値 sweep，巨大数・高精度函数の性能比較を通常testから分離して実行できる。

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

本ソフトウェアはBSD 3-Clause Licenseに定めるとおり，**現状のまま（AS IS）** 提供される。
商業的利用の適合性，特定目的への適合性，非侵害など，明示または黙示の保証は一切ありません。
著者または著作権者は，契約，不法行為，その他の理由で発生する，または発生したソフトウェアに関連するすべての請求，損害，またはその他の責任に対して，一切責任を負いません。

本ソフトウェアを使用することにより，ソフトウェアの使用に関して発生するすべてのリスクは自己責任であることを認め，自動的に同意したことになります。
著者は，データの損失，システムの不具合，その他ソフトウェアの使用によって生じた損害について一切責任を負いません。

正式な条件および免責事項は`LICENSE`を参照。

## 謝辞

このツール開発のきっかけとなった過去の自分，学びの園であった大学と教授殿，そして実使用環境となっている現職場に感謝の意を表します。

また、数値計算は一部[DLMF](https://dlmf.nist.gov/)を参考に実装してる。著者の皆様に深く感謝申し上げる。

## 要望等

要望・不具合報告・実装提案等はGitHubへ投げてください。大歓迎です🍀

## 将来のお話

- `for`的なものは欲しいよね
- `plot`函数（グラフ描画）
