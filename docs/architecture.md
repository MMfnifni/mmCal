# mmCal 内部構造

## 方針

mmCalは、入力をすぐ機械精度浮動小数へ変換せず、整数・有理数・記号式を可能な限り正確に保つ。近似計算、表示上の固定小数、数学的なセッション設定は別の責務として扱う。

依存方向は原則として次の順序を守る。

```text
numeric / symbols / expression
          ↓
mathematics
          ↓
simplification / symbolic / solver / approximation
          ↓
evaluation
          ↓
kernel
          ↓
frontend
```

CoreからCLIへ依存しない。数学層からKernelSessionへ依存しない。表示都合をExprやNumberへ逆流させない。

## 主要module

### `numeric`

`BigInt`, `Rational`, `RealNumber`, `Number`, `BigFloat`と近似値metadataを所有する。exactな整数・有理数の演算はこの層で閉じる。v1.5.1ではBigUIntの乗算をschoolbook/Karatsuba/Toom-3で適応dispatchし，巨大除算にはBurnikel–Ziegler，巨大10進変換にはdivide-and-conquerを使用する。

### `expression`

評価対象となるExpr、Call、Array、Symbolなどの構造を所有する。表示方法やユーザー入力位置は本体の数学値から分離する。

### `mathematics`

函数の定義域、逆函数、単調性、値域、周期、definednessなど、評価文脈に依存しない数学知識を`MathRegistry`へ集約する。`KnowledgeContext`はユーザー仮定や片側極限など局所的な事実を保持し，関係式の左右反転や安全な非零知識も共通推論として利用する。

### `simplification`

principal branchや定義域を壊さない範囲で式を標準化する。`Add`はAST全域のstrict total orderingで決定的に並べ，積・除算はdefinednessを保つ範囲で係数・分子因子・分母因子へ正規化する。`FullSimplify`は複数候補を探索し、式コストで選択する。

### `symbolic`

`D`, `integrate`, `limit`, 代数変形、多項式、置換を実装する。積分候補の一部は`D`を検証器として利用する。v1.5.1のderivative-back harnessはrule familyを横断して検証するが，証明器の能力不足だけで既存積分を拒否しないようStrict/ResolutionOnlyを分ける。

### `solver`

多項式，Rational function，制約付き解集合，実軸で安全な超越函数反転を扱う。完全解を証明できない場合は未解決状態を保持する。Unreleasedでは`SolutionBranch::freeVariables`をsolver変数の自由変数だけでなく，`where k in Integer`のような**branch-local formal parameter**にも一般化した。実軸`sin/cos/tan`の周期解はMathRegistryのperiod / principal inverse / real rangeを利用し，solve変数に対するexact affine argumentだけをinteger-parameter familyへ展開する。formal parameterはEnvironmentのユーザー変数ではなくSolutionSet内で局所束縛され，非線形argumentやComplex全解は安全な表現がない限り未解決のまま保持する。

### `approximation`

任意精度作業値と区間演算を使い、必要桁が保証できる数値近似を生成する。`Pi`はbinary-splitting Chudnovsky，`exp/log`はbinary splittingと保証付きrange reduction，巨大Radianの三角函数はPi保証区間によるargument reductionを使う。深すぎるASTはOSのstack overflowへ到達する前に拒否する。

v1.5.2では`N`の要求精度を子builtinへ伝播できるprecision-aware経路を追加した。これは全評価を近似化するモードではなく，明示対応したbuiltinだけが利用する。FFTではexact Exprを展開せず，`ComplexInterval`上のradix-2/Bluestein backendへ降りる。Matrixでも同じ`ApproximationContext`を使い，expression→certified interval変換，decimalization，guard-digit refinementを共通helperへ集約する。表示ではexact有限小数を不要に0埋めせず，certified fixed-digit結果の末尾0列は1桁だけ残して圧縮する。

Unreleasedの`DecimalApproximation` / `ComplexDecimalApproximation`は，真値保証用の**CertifiedEnclosure**と，後続計算で利用してよい情報量を表す**InformationEnclosure**を別々のexact Rational boundsとして保持し，常に`CertifiedEnclosure ⊆ InformationEnclosure`を保つ。`N[x,p]`の`p`は有効10進桁数であり，非zero表示値`d`の10進指数を`e=floor(log10(|d|))`とすると，丸め半量子`0.5*10^(e-p+1)`をInformationEnclosureへ含める。backend内部のguard桁は後からAccuracyとして回収しない。通常四則演算とcertified対応scalar函数では両enclosureを独立にinterval伝播し，exact `Number`は双方へ同じpoint intervalとして混在させる。ordered comparisonや`min/max`の離散判定はInformationEnclosureだけで証明できる場合に限る。出力値の正当性はCertifiedEnclosure，`accuracy` / `precision` / default `rationalize`と外側`N`の情報量制限はInformationEnclosureを基準にする。結果自身へ伝播済みInformationEnclosureを保存するため，複数演算を跨いでも単なる要求桁数へ情報を圧縮し直さない。`:fix`はこれとは独立した固定小数表示である。

### `linear_algebra`

rank-2 Array上の線形代数algorithmをbuiltin dispatchから分離する。`MatrixView`は`ArrayExpr`の論理要素列を参照し，書換えが必要なalgorithmだけ`MatrixBuffer`へ複製する。exact Number行列はpivot loop内でExpr/Simplifierを使わない専用backend，`N`配下ではBigFloat/`ComplexInterval`系のcertified backendへ分岐する。approximate Matrix側には既に`PointMatrix` / `IntervalMatrix` / `ComplexMatrix`等の連続working bufferがあるため，永続Array storageとalgorithm temporaryは分離する。

ユーザー構文の`{...}`は一般の有限brace containerとする。child shapeが一致する矩形値はdense `ArrayExpr`へ自動最適化し，shapeが異なる`{Q,R}`やragged値は`ListExpr`として保持する。Unreleasedの`ArrayExpr`は固定1024要素のimmutable pageをshared backingとして持ち，各pageをInteger / Rational / Number / DecimalApproximation / ComplexDecimalApproximation / Genericの最狭表現でpacked保持する。Array全体の意味論上の型は常に`ExprKind::Array`であり，page種別はstorage detailである。shape / offset / stridesをbackingから分離したため，transposeはstride交換だけのzero-copy view，contiguous reshapeや一部sliceもbackingを共有する。Matrix algorithmは`ArrayExpr`だけを受け取り，Evaluator dispatchの矩形性監査で`ListExpr`をWarning + 未評価へ戻す。zero-length dimensionはdense Arrayのshapeとして保持する。`dimensions`は一般braceに対して全childに共通するrectangular prefixを返し，`length` / `at`はArray/List双方を扱う。

`ArrayBuilder`は完成済みpageをimmutable化し，promotionを現在の最大1024要素page内に限定する。したがって大規模numeric Arrayの末尾でsymbolic値が現れても，全要素をGeneric `Expr`へ再構築しない。矩形brace literalはLowererから単一builderへleafを直接投入し，numeric literalを一度`Expr`化してからpackする二重表現を避ける。

一般symbolic `det` / `inverse`はexact-firstを維持する一方，無制限Laplace/adjugate展開は行わない。三角行列fast pathと共有展開budgetにより，式爆発が見込まれる場合は未評価式へ戻す。

分解系のLU / QRは`linear_algebra/decomposition.*`へ集約する。`luDecomposition`は正方行列でrow-pivoted `P A = L U`（certified approximateでは非零候補の|pivot|保証下限を比較するpartial pivoting），`qrDecomposition`は矩形を含むreduced Householder QRとして`k=min(m,n)`，`Q:m×k`，`R:k×n`を返す。公開結果は一般brace `{P,L,U}` / `{Q,R}`であり，同shapeなら内部dense Array，異shapeならList表現になる。`at[result,i]`でfactorを切り出す。`N[...]`ではexact分解を構築せずcertified `ComplexInterval` backendへ直接dispatchする。Householder適用にはcolumn-block kernelも持つが，no-BLASのBigFloat/interval backendではblock=8/16/32の実測優位が一貫しなかったため既定はunblocked相当とし，benchmarkだけ常設する。

SVDは`linear_algebra/svd.*`へ分離し，reduced `{U,S,V}`を返す。一般数値backendは`A^H A`を形成せず，Householder bidiagonalization + one-sided JacobiをBigFloat中心値で行い，元入力のcertified intervalと候補factorからreconstruction residualおよび`U^H U` / `V^H V` orthogonalityを区間監査する。実数・複素数双方を扱い，複素点演算は`linear_algebra/complex_point.*`へ共通化する。監査が要求桁を満たさなければguard digitsを増やして再試行し，exact SVDは自然に閉じるcaseだけを返す。重複特異値ではvector basisが一意でないため，componentwiseな唯一性ではなく再構成・直交性を保証対象とする。

Eigenは`linear_algebra/eigen.*`へ分離する。exact pathは上三角行列の対角固有値，対角行列の標準基底，distinct-root exact Number 2×2を扱い，一般数値pathはComplex BigFloat上でHessenberg reduction → implicit shifted QR → complex Schur形へ進む。Schur vectorを蓄積し，固有vectorは上三角Schur行列からback substitutionして列として返す。停止精度は出力桁より十分厳しく設定し，元入力の`ComplexInterval`に対して`A Q-Q T`および`A v-λv`をinterval演算で監査する。一般非正規行列では固有量のcomponentwise enclosureを安易に主張せず，Schur/eigenpair relationをcertificate境界とする。重根・defective/near-defective caseで安定な独立固有vectorを作れない場合は未評価へ戻す。

v1.5.2の1024×1024 dense Matrix監査では，算法より先に`Expr::Node` / Rational / parse-loweringの固定費がmemory bottleneckになることを確認した。Unreleasedの第一段階では公開`Expr` APIを維持したtyped-node化により，同一GCC Release/LTO-offの1024×1024 Rational Matrix `transpose`で最大RSSを`693312 KiB`から`299668 KiB`へ約56.8%削減した。第二段階ではpersistent Arrayをimmutable paged packed backing + stride viewへ移行し，benchmark入力も`ArrayBuilder`から直接packed構築する経路へ変更した結果，同負荷の最大RSSは約`136576 KiB`（133.4 MiB），transpose本体は約0.059 msとなった。単一flat packed vector案はtransposeのdeep copyで退行したため採用していない。BigUInt/BigInt SBOは保守性を優先して今回は見送り，必要なら別experimentとして測定する。

Stage 3ではexact実数（整数/Rational）行列を行ごとの分母LCMで整数行列へliftし，`IntegerMatrixBuffer`上のBareiss fraction-free eliminationへdispatchする。分母除去は共通`liftRealRows` helperへ集約し，`det`はBareissの最終pivotから復元，`rref` / `matrixRank` / `nullSpace`はfraction-free forward eliminationを共有する。`nullSpace`はfree columnを昇順に選ぶRREF basisを構成し，full column rankでもshape `{0,n}` を保持する。`inverse`は `B=D A` に対するaugmented matrix `[B|D]`，`solveLinear[A,b]`は `[A|b]` を同じkernelへ渡し，後者ではpivot候補を係数列だけに制限して整合性と一意性を判定する。これにより中間Rational生成をpivot loopからほぼ排除する。exact複素行列は現在も`Number` Gaussian backendへfallbackする。`N[solveLinear[...],p]`はexact解を先に構築せず，certified interval augmented eliminationを直接試す。一方`matrixRank` / `nullSpace`はrank deficiencyに不連続なので，exact入力ではexact pivot structureを優先し，近似入力ではintervalでpivot構造を証明できる場合だけ結果を確定する。

### `evaluation`

Builtin属性、Hold規則、iterator、代入、ユーザー函数、履歴参照、診断を統合する。評価器は深い通常式でC++再帰stackを消費しにくい明示task-stack方式を維持する。`N`は特殊taskとして第1引数を保持し、precisionを先に確定してから子式を評価する。precision contextはstack管理されるためnested `N`でも外側の要求精度を破壊しない。

### `kernel`

1セッションのユーザー定義、履歴、角度設定、乱数状態、Warning/Infoを所有する。通常の`reset()`はRNGをentropy reseedする一方，benchmark/fuzzer等で独立評価を連続実行する用途には`resetForIndependentEvaluation()`を使い，定義・履歴・入力番号・diagnosticだけを初期化してRNG streamと角度設定を保持する。

### `cli`

標準入出力、`:fix`, `:status`, 起動時引数、console titleを担当する。`:fix`は表示だけを変え、KernelのExprや履歴を書き換えない。

## 数学知識の共有

`MathRegistry`は函数固有の安定したmetadata、`KnowledgeContext`は評価時の仮定、`ValueFacts`はExprから導出した符号・実数性などの問い合わせを担当する。

`D`, `integrate`, `limit`, `solve`は別algorithmのまま維持するが、domain・inverse・range・periodicity等は可能な限り共通metadataを参照する。

## 表示

通常formatterは再parse可能なcompact表記を生成する。ASTをそのままdumpすることは目的にしない。

- `a+(-b)`は`a-b`
- 不要な演算子空白は出さない
- `xy`が1つのSymbolへ読まれる場合など、字句境界が必要な積では空白を残す
- `2exp[x]`，`Pi^0x`，Array隣接などlexer/parserと衝突する境界では明示`*`を使う
- 優先順位・結合規則を守るため必要な括弧は残す

小数表示指定はCLI presentationであり、内部の正確値を近似値へ置換しない。指定桁数へ丸めた後、表示上不要な小数部末尾の0だけを除去する。

## stack安全性

特に監査対象とする処理:

- CertifiedEvaluatorのAST深さ
- Simplifierの再入
- 高階微分
- 積分候補検証
- l'Hopital反復
- solver recursion
- formatter traversal

病的入力は通常の評価失敗・未対応として処理し、OS-level stack overflowへ到達させないことを目標とする。

## Benchmark project

`mmCal.Benchmarks`は通常の回帰testとは分離したConsole projectである。`mmCal.Core`へだけ依存し，次を担当する。

- 固定seedの巨大BigInt商余り・10進round-trip・exact/certified Matrix・FFT等のランダム正当性試験
- grammar-aware semantic Random Expression Fuzzer。caseごとに再現可能なseedを持ち，各workerの`KernelSession`を独立resetして長時間探索する
- Karatsuba / Toom-3 / Burnikel–Ziegler等のthreshold sweep
- factorial，decimal conversion，高精度`Pi/exp/log`，exact/certified Matrix，FFT等の速度比較
- `--full`による大規模case，`--random-only` / `--benchmark-only`による用途分離

性能測定をUnit testのPASS/FAIL時間へ混ぜず，算法選定の根拠を再現可能に残すことが目的である。採用・棄却履歴は`performance_optimization.ja.md`を参照する。
