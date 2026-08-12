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

多項式、Rational function、制約付き解集合、実軸で安全な超越函数反転を扱う。完全解を証明できない場合は未解決状態を保持する。

### `approximation`

任意精度作業値と区間演算を使い、必要桁が保証できる数値近似を生成する。`Pi`はbinary-splitting Chudnovsky，`exp/log`はbinary splittingと保証付きrange reduction，巨大Radianの三角函数はPi保証区間によるargument reductionを使う。深すぎるASTはOSのstack overflowへ到達する前に拒否する。

### `evaluation`

Builtin属性、Hold規則、iterator、代入、ユーザー函数、履歴参照、診断を統合する。評価器は深い通常式でC++再帰stackを消費しにくい明示task-stack方式を維持する。

### `kernel`

1セッションのユーザー定義、履歴、角度設定、乱数状態、Warning/Infoを所有する。

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

- 固定seedの巨大BigInt商余り・10進round-trip等のランダム正当性試験
- Karatsuba / Toom-3 / Burnikel–Ziegler等のthreshold sweep
- factorial，decimal conversion，高精度`Pi/exp/log`等の速度比較
- `--full`による大規模case，`--random-only` / `--benchmark-only`による用途分離

性能測定をUnit testのPASS/FAIL時間へ混ぜず，算法選定の根拠を再現可能に残すことが目的である。採用・棄却履歴は`performance_optimization.ja.md`を参照する。
