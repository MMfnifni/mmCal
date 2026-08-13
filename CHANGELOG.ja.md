# Changelog

## v1.5.2 — 開発中

- Array / Matrix + Bareiss / `solveLinear` / `nullSpace` / LU / Householder QR / SVD / Eigen 回帰: internal `2027 / 2027 PASS`、black-box `1504 / 1504 PASS`。`mmCal.Benchmarks --random-only` の Matrix / FFT を含む固定seed検証も PASS。

- `RealInterval::point`で同一`BigFloat`を一つの初期化式内でcopy/moveしていた構築を分離し，負のpoint intervalがMSVCで反転し得る問題を修正。SVD/Eigenのinterval監査で発生していた`RealInterval lower bound exceeds upper bound`を防止し，負値pointの専用回帰を追加。

### 構文（breaking change）

- 函数呼び出しを `name[...]` に一本化し、`name(...)` 呼び出しを廃止
- 丸括弧 `()` はgrouping専用とし、通常identifierの `x(x+1)` は暗黙乗算として扱う
- 既知の函数名に旧 `sin(x)` 構文を使った場合は、誤って乗算へ解釈せずSyntaxErrorを返す
- Formatterは函数を常に `name[...]`、identifierとgroupの積を `x*(...)` と明示してround-tripを一意化
- Parser/AST/Lowererから函数呼出delimiterの分岐を削除
- 履歴参照を`In[n]` / `Out[n]`へ整理し，負添字による相対参照を追加（`0`は無効）
- `@` / `@@` / ... を`In[-1]` / `In[-2]` / ... の短縮記法として追加し，`%` / `%%` / ... と`Out[-n]`の対応も明文化
- 負の`Out[-n]`は成功出力，負の`In[-n]`は入力slotを基準に数える

### Array / 基本線形代数

- Arrayをshape + row-major flat storageの共通表現として整理し，評価後に生じるnested Arrayも矩形性を検証して自動flattenする不変条件を`Expr::array`へ集約
- zero-length dimensionを正式に保持し，brace表記だけではshapeを保存できない空Arrayを`reshape[{}, {...}]`でround-trip可能にした
- `dimensions` / `arrayRank` / `at` / `reshape`を追加し，0-based indexingを共通Array APIへ統合
- `MatrixView` / `MatrixBuffer`を追加し，row-major Arrayを不要に`vector<vector<Expr>>`へ複製しない線形代数基盤へ整理
- canonical APIを`dot` / `matrixRank` / `norm` / `normalize`へ統一。旧`matmul` / `mmul` / `vdot` / `rank` / `mrank` / `vnorm` / `vnormalize` / `mget`は互換aliasとして維持
- 同shape Arrayの`+` / `-`とscalar×Arrayを通常算術へ統合し，Array×Arrayの`*`は拒否して`dot[...]`を明示的な contraction とした
- exact Number行列用にExpr/Simplifierをpivot loopへ持ち込まないflat Gaussian / Gauss-Jordan backendを追加し，`det` / `inverse` / `rref` / `matrixRank`を共通化
- Stage 3として整数/Rational行列へBareiss fraction-free eliminationを追加。Rationalは行ごとの分母LCMで整数liftし，`det` / `rref` / `matrixRank` / `inverse`が共通`IntegerMatrixBuffer` kernelを利用
- Bareiss除算は`divmod`でexactnessを検証し，pivot候補は中間BigInt growthを抑えるためbit length最小を優先。exact complexは従来Gaussian backendへfallback
- `solveLinear[A,b]`を追加。正方系だけでなくfull column rankの整合した過剰決定系も一意解として扱い，不整合系・自由変数を含む系はDomain errorとする
- `nullSpace[A]`を追加。free column昇順のcanonical RREF basisを返し，full column rankでは`reshape[{}, {0,n}]`として空basisのvector次元を保持する。整数/RationalはBareiss forward eliminationを共有し，exact complexはGaussian，symbolicはpivot非零性を証明できる場合だけ処理
- exact整数/Rationalの`solveLinear`は共通の行別分母除去とaugmented Bareiss `[A|b]`を利用し，pivot候補を係数列へ制限して整合性・一意性を判定。exact complexはGaussian，symbolicはpivot非零性を証明できる範囲だけ処理
- 同一Release benchmarkでStage 2比、`det`は8–16次で約5.4–11.4倍、`rref`は約7.7–15.5倍高速化
- symbolic `det` / `inverse`へ共有展開budgetと三角行列fast pathを追加し，dense高次行列の階乗級expression explosionを抑止
- FFTとMatrixで共有するcertified expression→interval / decimalization / guard-digit管理を`approximation`層へ抽出
- `N[dot[...],p]` / `N[det[...],p]` / `N[inverse[...],p]` / `N[rref[...],p]` / `N[matrixRank[...],p]` / `N[norm[...],p]`等をprecision-aware certified Matrix backendへ直接dispatch
- `N[solveLinear[...],p]`もexact解を先に構築せずcertified augmented interval eliminationへ直接dispatchし，pivot/整合性を証明できない場合はepsilon判定を行わない
- `matrixRank`は不連続量のためexact入力ではexact eliminationを優先し，近似backendでもepsilon閾値を使わない。非零pivotをintervalで証明できる場合だけrankを進め，近似値からrank deficiencyを推測しない
- `nullSpace`もrank deficiencyに不連続なためexact入力ではexact pivot structureを優先し，近似済み入力ではintervalでpivot構造を証明できる場合だけbasisを返す
- `at`をprefix indexingへ拡張し，rank-3 decomposition結果から`at[result,0]`のようにfactor subarrayをzero-basedで取得可能にした
- `luDecomposition[A]`を追加。正方行列に対してrow-pivoted `P A = L U`をexact-firstで構成し，certified approximate側は`|pivot|^2`保証下限最大のpartial pivotingを使い，`N[...]`ではcertified interval LUへ直接dispatch
- `qrDecomposition[A]`を追加。Householder reflectorによる`A = Q R`をexact-firstで構成し，`N[...]`では実/複素のcertified interval Householder QRへ直接dispatch。一般exact QRは4×4で式爆発が実測されたため3×3以下へpolicy制限し，上三角行列は任意サイズでfast path
- Householder適用へcolumn-block kernelを試作し，block=1/8/16/32をRelease計測。8/16/24次で一貫した優位がなかったため自動block化は不採用とし，kernel/benchmarkのみ保持
- `N`表示を整理。exact有限小数は従来どおり不要な0埋めをせず，certified fixed-digit結果の末尾0列は1個だけ残して圧縮する（`1.000... -> 1.0`，`1.500... -> 1.50`）。要求桁数とcertified enclosureはmetadataに保持する
- `mmCal.Benchmarks`へexact/Rational `solveLinear` round-trip / `nullSpace` basisの`A.v==0`検証 / inverse round-trip / determinant transpose invariance / RREF / certified determinantを含むfixed-seed Matrix random checkとMatrix timingを追加
- `{...}`を一般有限brace containerへ拡張し，矩形childはdense Arrayへ自動最適化，shapeの異なる`{Q,R}` / `{U,S,V}`等は一般braceとして保持。Matrix境界で矩形性を監査
- `qrDecomposition`を矩形reduced QRへ拡張し，`svd` / `singularValueDecomposition`を実・複素reduced SVDとして追加。`A^H A`を形成せずHouseholder bidiagonalization + one-sided Jacobiを採用
- `eigenvalues` / `eigenvectors` / `eigensystem`を追加。exactは三角/対角/distinct-root Number 2×2を優先し，一般`N[...]`はComplex BigFloat Hessenberg + implicit shifted QR + Schur backendへ直接dispatch。元入力intervalに対するSchur/eigenpair relationを監査
- `expression_interval.cpp`が直接使用する`std::overflow_error`の宣言元として`<stdexcept>`を明示includeし，transitive include依存を除去
- `mmCal.Benchmarks --matrix-large`を追加。添付generator相当の10桁random decimal行列で32/64次の主要Matrix函数を測定し，1024次ではCLI parseだけで約10.9秒/1.99 GB RSS，direct Expr入力でも約0.69 GBとなるstorage bottleneckを確認


## v1.5.1 — 2026-08-12

v1.5.0のexact-first CAS基盤を維持しつつ、canonicalization、検証、巨大整数、高精度数値評価、benchmark基盤を重点的に改善した。

### 数式・CAS

- `Add`へAST全域の決定的strict total orderingを導入
- 積・除算をdefinednessを保つcanonical normal formへ整理
- `MathKnowledge`の非零知識と関係式左右反転推論を強化
- formatter/parser生成型round-trip testを追加し、radix prefix衝突、Array隣接、負Rationalの表記不安定を修正
- 積分rule/familyを横断するderivative-back harnessを追加
- Referenceとbuiltin registryの自動照合を追加
- approximation extreme-value testを拡充
- FFT plan/twiddleのtransform間cacheを追加

### 多倍長・高精度

- BigUInt乗算をschoolbook / Karatsuba / Toom-3の適応dispatchへ変更
- 専用squareを追加
- factorialはbalanced product treeを維持し、leaf構築・1-limb経路を高速化
- Burnikel–Ziegler divisionと`2^k`除算fast pathを追加
- decimal parse/toStringを`10^9` chunk + divide-and-conquer化
- `tryToUint64`の巨大値早期棄却を追加
- BigFloat extreme exponent-gap加減算をdirected roundingを保ったままfast-path化
- `Pi`をbinary-splitting Chudnovskyへ変更
- `exp/E`と`log`をbinary splitting + certified range reductionへ変更
- 巨大Radianの`sin/cos/tan`へcertified argument reductionを追加

### Benchmark / test

- Visual Studio solutionへ`mmCal.Benchmarks`を追加
- fixed-seed random invariant、threshold sweep、factorial、decimal I/O、高精度`Pi/exp/log`のbenchmarkを常設化
- v1.5.1確定時点: internal 1691 / 1691、black-box 1337 / 1337

### 比較したが採用しなかったもの

- Prime-Swing factorial
- binary GCD
- Karatsuba vector-pool workspace
- Karatsuba recursion-depth scratch workspace
- Toom-3専用square
- 低threshold Toom-3
- machine `fmod`による巨大trig縮約

各判断の実測根拠は`docs/performance_optimization.ja.md`を参照。

---

## v1.5.0

旧版から数値モデル、Lexer/Parser/AST/Evaluator、Simplifier、Solver、CertifiedEvaluator、CLI、formatter、tests、documentationをほぼ全面再構築し、mmCalをexact-first CLI calculator / compact CASとして再定義したrelease。
