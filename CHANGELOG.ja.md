# Changelog

## Unreleased

### 意味論・ビルド

- exact複素数の四則演算後に虚部が0になった場合，`Number`を実数表現へ正規化するよう修正した。あわせて宣言だけ存在していた有限精度の楕円積分と複素`li`の保証付き実装を補い，CMake/GCCでも全対象をリンクできる状態へ戻した。
- `mmCal.Benchmarks`へ保証付き数値計算の境界fuzzerを追加した。分岐切断，exactな極，有限精度の境界交差，実数／複素数の計算経路，計算量閾値を46種の試験から生成し，`Value` / `DomainError` / `N::precision` / `N::unsupported` / 未評価 / 資源超過を区別する。`--seed --case`で再現でき，停止監視も数学的失敗とは分離して扱う。Visual StudioプロジェクトとCMakeのソース構成も同期した。
- `simplify` / `fullSimplify`を定義性を考慮するよう強化した。`F-F -> 0`，`F/F -> 1`，`F^0 -> 1`や特殊函数の0/1への退化は，必要な定義条件を証明できる場合だけ適用する。`0^0`，負のRational指数，`zeta[1]`，Gamma系の極なども同じ原則で扱い，未証明の定義域の穴を消さない。
- 積分のderivative-back検証を，共通定義域上での恒等式確認へ変更した。Random Expression Fuzzerにも`limit`，`cases`，入れ子の`N`，AlgebraicNumber，Gröbner還元，Array reshape，除外点・可除特異点の意味論的不変条件を追加した。
- MSVCのstack reserveを`mmCal` / `mmCal.Tests` / `mmCal.Benchmarks`で16 MiBへ統一した。builtinの引数個数検査・未評価保持，Rational丸め等の共通処理も集約した。exact複素数の`digamma` / `trigamma`には，整数実部をexactな漸化式で戻す経路を追加し，不要な区間相殺を避ける。

### 保証付き `N`・精度情報

- exact入力に対する`N`ではCertifiedEnclosureとInformationEnclosureが同一なので，重い保証付き函数を二重評価しないようにした。複素`li`の`Log -> Ei`合成と`2F1`の主値`1/z`接続公式でも内部ガード精度の二重加算を整理し，主値・区間保証を維持したまま代表20桁入力を大幅に高速化した。
- `leastSquares`が内部の擬似逆行列を要求表示桁へ一度丸めてから積を取っていたため，有効桁を余計に1桁失う問題を修正した。中間擬似逆行列は作業桁まで保持し，最終結果だけを要求桁へ丸める。有限precisionの`nullSpace`については，pivotの存在だけでなくfree columnであることもInformationEnclosureから証明できる場合だけbasisを返す契約を回帰テストへ明記した。
- `acosh`の分岐切断 `-1<x<1` の直上・直下で，極小虚部を含む入力が過剰精密化して時間切れになる問題を修正した。分岐切断のどちら側にあるかを区間で証明できた場合だけ安定な式へ切り替え，`N[acosh[1/2+I/10^160],100]`級でも主値枝を保ったまま短時間で確定する。
- `expm1` / `log1p`の0近傍で，内部の`exp[x]-1` / `log[1+x]`が要求桁より先に相殺する問題を修正した。入力の2進桁位置から必要な作業精度を増やし，極小exact入力でも有効桁を維持する。`2F1`はexactな実数`z>1`について主値の分岐切断上で定めた規約値を許可し，有限精度入力が分岐切断の上下を含む場合だけ`N::precision`とする。複素Fresnel級数には`|z|<8`の有界計算境界を設けた。
- `diff[...,digits]` / `nintegrate[...,digits]`の`digits`を従来どおり小数部桁数として保持した。有限精度入力の`InformationEnclosure`による精度上限は維持しつつ，`N`の有効桁表示規則へ誤って変更されないよう分離した。
- 特異点・分岐切断・算法境界の判定を`CertifiedEnclosure`と`InformationEnclosure`で分離した。exactな極はDomainError，有限精度入力の情報幅のため極や分岐側を確定できない場合は`N::precision`，値は存在するが現行算法で扱えない場合は`N::unsupported`とする。Gamma系，`Ei` / `Ci` / `li`，`zeta`，`1F1` / `2F1`，`ibeta`，`log` / `sqrt` / 非整数冪，逆三角・逆双曲線，`Arg` / `atan2`，`polylog`，楕円積分を同じ規則へ揃えた。
- 保証付き数値層を整理し，値表現，零・極・分岐切断判定，実数／複素数の四則演算，精度計算を分離した。数学的仕様は変えず，重複処理と旧comment-out実装を削除した。
- `DecimalApproximation` / `ComplexDecimalApproximation`で`CertifiedEnclosure`と`InformationEnclosure`の役割を明確化した。真値がexact zeroであることと，入力情報としてexact zeroを再利用できることを区別し，`precision` / `accuracy` / `explain` / 比較 / 数値演算は同じ情報量評価を共有する。`explain`には`PrecisionDigits` / `AccuracyDigits`を追加した。
- 有限precision値の隠れた保護桁の情報を再利用する経路を横断監査した。`N[0,p]^0`や`1/N[0,p]`の零／非零判定，複素値が実数へ射影された場合の表示桁，cardinal函数の0近傍，FFTの相殺，行列pivot・階数判定をInformationEnclosure基準へ統一した。連続量はCertified/Information両区間を並行伝播し，`lu/qr/svd/conditionNumber/pseudoInverse/leastSquares/eigen*`は既に有限precisionの行列に対する摂動保証が未完成なため，内部の点値から結果を復元せず保守的に未評価へ戻す。
- top-level `N`のrefinementを局所16回までに制限した。有限精度入力が分岐切断，極，算法境界を跨ぎ続ける場合は無制限に保護桁を増やさず`N::precision`へ戻す一方，exact入力から追加精度で側を証明できる場合は従来どおり継続する。
- 逆三角・逆双曲線，`Arg` / `atan2`，`polylog`，楕円函数`F/E/Pi`の主値枝と実数／複素数境界をInformationEnclosureで判定するよう統一した。複素逆函数では分岐点から十分離れた場合に導関数上界を使って区間の過剰拡大を抑え，条件数に由来する正当な精度低下は残す。
- 近似値のCLI表示を内部精度情報から分離して簡潔化した。0近傍は`0.0`，末尾0列は必要最小限とし，`N[1/2,20] -> 0.50`，`N[2,20] -> 2.0`，`N[I,20] -> 1.0I`，`N[log[-1+I/10^1000],20] -> 0.0+3.1415926535897932385I`のように表示する。`+0` / `-0` / `*1` / `/1`等では精度情報を再量子化せず，0近傍の表示量子は相対精度ではなく絶対精度から決める。
- 実数／複素数の評価経路を揃え，全引数が実数の`1F1` / `2F1` / `polylog`の安全な実数射影，楕円積分の`RealInterval`対応，極近傍`2F1`の有界な区間級数を追加した。exact Rationalの`2F1`数値級数も区間累積へ切り替え，極に近い分母パラメータで巨大なRational分子・分母を生成しない。`diff` / `nintegrate`も`CertifiedEnclosure`と`InformationEnclosure`を並行して伝播し，holdされた`N[...]`を含めて入力以上の精度を作らない。入力情報だけでは特異点や分岐側を決められない場合は，無意味な再精密化をせず評価不能として返す。`N`専用のブラックボックス・性質回帰も追加した。

### ブラックボックス検証

- public CLIだけを介する`D` / `integrate` / `limit` / `N` / `cases` / Gröbner / exact代数方程式Solveの監査と，derivative-back恒等式，Gröbner不変条件，Solve解数，精度情報，漸化式，高精度参照値を照合する性質検証を追加した。
- `cases[...]`を跨ぐ`D` / `integrate` / `limit`は，条件が微積分変数へ依存しない場合だけ枝ごとに分配する。変数依存条件や可除・非可除特異点では，無効な枝から偽の`DomainError`を出さず，未評価または正しい極限を保つ。
- Holdされた多項式函数内で`groebnerBasis[...]`を安全に合成できるようにし，旧`if[...]`表現を使っていた期待値・Reference例を`cases[...]`へ同期した。

### 場合分け・多項式イデアル・複素ポリガンマ

- 数学的場合分けを評価制御`if[...]`や`SolutionSet`から分離し，scalar `cases[value if condition; ...]`を追加した。未確定条件を保持し，`simplify` / `N` / `D` / `integrate` / `limit`へ接続した。
- `Q[x1,...,xn]`上の一般多変数多項式基盤と，Lex / GrLex / GrevLex，除算・標準形・S多項式・Buchberger法を用いるexact Rational係数Gröbner基底を追加した。公開函数は`groebnerBasis[...]` / `polynomialReduce[...]`で，既存の評価資源上限内で処理する。
- 非線形多項式`solve[{...},{...}]`をLex Gröbner消去へ接続した。矛盾系は空集合，0次元の適合系は一変数exact rootと後退代入で列挙し，元方程式をexactに再検証する。正次元多様体は偽のパラメータ表示を作らず`UnresolvedSolutionSet`へ戻す。
- 複素数`digamma` / `trigamma`の保証付き`N`を追加し，漸化式とBernoulli/Stirling型漸近展開で評価する。一般`polygamma[n,x]`は将来課題のままとした。Visual Studioプロジェクトにも関連ソースを登録した。

### 記号計算・表示修正

- 条件付き解集合の`&&`表示，入れ子iterator形式の`table`，exact複素数のRational虚部表示（`2I/29`）を整えた。
- `factor`へexact一変数多項式の完全冪認識を追加し，`D`内の`%` / `Out[n]`を副作用なく履歴値へ解決してから微分するよう修正した。`In[n]`の再評価意味論は変更しない。
- 積分へ一般冪則`integrate[x^n,x]`，`integrate[log[log[x]],x]`，主値枝を保つ`integrate[li[x],x]`を追加した。`integrate[log[log[x]],{x,1,E}]`も収束する広義積分としてexactに扱う。
- Formatterのradix-prefix衝突判定を字句境界基準へ修正し，`380x^9`等の不要な`*`を除いた。単変数多項式は表示時だけ次数降順へ並べ，内部Exprの順序は変更しない。
- `limit[expr,{x,a,direction}]`を追加し，`Ei` / `Ci` / `li`の主値枝上の既知極限と`li[0] -> 0`を補強した。`sin` / `cos` / `tan`の周期振動は極限不存在を`Indeterminate`で表し，限定的なはさみうち則により`limit[x sin[1/x],x,0] -> 0`をexactに処理する。

### 微分・積分・複素数 `N` 監査

- 微分公式，原始函数，特殊函数恒等式，主値枝，可除特異点を横断監査し，Beta系，有限組合せ函数，`Ei/Si/Ci/li/digamma/trigamma/LambertW`，`polylog`，`1F1/2F1`，`ibeta`，`sinc/cosc/expc`等の規則を補強した。
- `sinc/cosc/tanc/sinhc/tanhc/expc`等は0で偽の`0/0`を作らないよう`cases[...]`で連続延長値を保持する。`Si'`，主値`LambertW`，`polylog`も0での有限導函数値を維持する。
- 複素数の保証付き`N`を`ComplexInterval`上へ拡張し，`erf/erfc`，`Ei/Si/Ci`，Fresnel函数，`1F1`，`2F1`，`polylog`，`zeta`，`gamma`を級数・Euler-Maclaurin・Stirling等で評価する。`2F1`は条件を証明できる場合のみ安全な`1/z`接続公式を用いる。
- 負の非整数Rational `digamma/trigamma`と負Rational `zeta`はexactな漸化式・函数等式を経由して既存領域へ送る。`lgamma`は引き続き実軸の`log|Gamma|`である。derivative-back監査も局所恒等式を全域恒等式へ誤拡張しない。

### 保証付き数値計算・算法閾値監査

- `PrecisionInsufficient`と`CertifiedBackendUnsupported`の使い分けを整理し，追加精度で解決しない固定級数範囲・項数・計画上限は`N::unsupported`へ戻す。区間幅だけが分岐切断や計算量境界を跨ぐ場合は精度の精密化を許す。
- 特殊函数の級数境界と相殺対策を再調整した。実数`Ei/Si/Ci`は96，複素`Ei`は`|z|<=512`，複素`Ci`は`|z|<=128`を計算量上限とし，`Si/Ci`と複素`Ei/Ci`は引数の大きさから相殺に必要な作業精度を追加する。`1F1: |z|<=160`，`2F1: |z|<=9/10`，楕円積分`F/E: |m|<=9/10`（`Pi`は`|n|<=9/10`も要求），正整数位数`polylog: |z|<=49/50`も維持する。
- `zeta`は`s=1`だけをDomainErrorとし，現Euler-Maclaurin実装が扱わない領域を数学的な定義域外と誤分類しない。
- 保証付きFFTの精度再試行を12回に制限し，非2冪FFTの直接DFT / Bluestein境界は384点を維持した。exact線形代数およびBigUIntのKaratsuba / Toom-3 / Burnikel-Ziegler / 10進変換についても，算法切替境界を回帰テストで直接検証する。

### Exact線形代数・計算取消し

- `conditionNumber`，`pseudoInverse`，`leastSquares`を追加した。exact行列では階数を厳密に判定し，擬似逆行列は階数分解からMoore–Penrose逆行列をexactに構成する。exact入力を外側`N`から評価する一般数値経路は保証付きSVDを使う。既に有限precisionの行列では階数・特異部分空間を隠れた保護桁から推測せず，現状は保守的に未評価へ戻す。0次元行列の形状も保持する。
- exact整数/Rational行列へ31-bit素数体 + CRTによるmodular経路を追加した。`det`はHadamard上界まで復元し，`solveLinear`は有理数再構成後に`A X = B`をexact検証する。`inverse`もadjugateをexact検証し，失敗時はBareissへ戻す。
- GCC実測に基づき`det` / `solveLinear`のBareiss・modular自動選択を追加した。`inverse` / `rref` / `matrixRank` / `nullSpace`は現状Bareissを維持する。使用素数数のtelemetryと`--exact-linear-algebra` / `--budget-telemetry` benchmarkも追加した。
- 最上位評価へ`EvaluationCancellationToken`を接続し，WindowsのCtrl-C / Ctrl-BreakとPOSIXのSIGINTで協調的に停止できるようにした。取消しは`ResourceLimitError`として扱い，非対話時の終了codeは`3`を維持する。

### 評価資源制限

- 1回の最上位評価で共有する`EvaluationBudget` / `EvaluationLimits` / `EvaluationUsage`を追加した。評価step・depth，生成Expr，Simplifier/Solve/積分候補，certified refinement，Array/Matrix要素，BigInt bit長，要求精度，代数的構成などを共通に計測・制限する。
- `KernelSession::setEvaluationLimits` / `evaluationLimits` / `lastEvaluationUsage`を追加し，従来の`setEvaluationDepthLimit`も互換APIとして接続した。資源超過は`DomainError`や未評価とは分離し，`ResourceLimitError`で資源名と上限を報告する。
- 巨大整数冪・factorial・Array・Matrix・`N[expr,p]`・入力text等は，大規模確保や計算を始める前に可能な範囲で上限を検査する。実時間上限はCoreへ入れず，決定論的な演算量budgetとフロントエンド側の取消し／タイムアウトを分離した。

### ファジング・検証

- Random Expression Fuzzerをderivative-back，Solve，`A inverse[A] == I`，`ifft[fft[v]] == v`，DomainError分類，主値`sqrt`境界まで拡張した。多項式は構造残差と独立代入，Solveは`SolutionSet`の集合比較で検証する。
- exact FFTは構造一致を優先し，相殺形が残る場合だけ`CertifiedEvaluator`へ直接渡す二段目判定を使う。深い式は再帰depth上限で停止し，判定不能を数学的FAILと混同せず`inconclusive`として数える。
- FAIL時にseed / case / 簡約式と各budget使用量を表示し，資源境界・reset・telemetryの回帰テストを追加した。

### CLI・互換性

- 行単位自動化向けの`--batch`を追加した。
- `--help`は短い起動案内のまま，REPLの`:help`を全callable函数の説明・入力規則・例を持つ一覧へ拡張した。`:help Pi` / `:help constants`，未知項目の近傍候補提示も追加し，help参照では評価・履歴を消費しない。
- `mmCal.Benchmarks --fft-threshold [iterations]`を追加し，直接DFTとBluesteinを65～509点で比較した。GCC / MSVCの結果から保証付き非2冪FFTの方針境界を384点へ再調整し，強制比較テストでも維持した。

## v1.5.3 — 2026-08-16

v1.5.3は，v1.5.2の記号微積分，保証付き`N`，Array・線形代数を基礎に，**exact代数数体，Solverの意味論統一，高精度特殊函数，exact Cyclotomic FFT，内部表現の最適化**を接続した版である。exact-first・proof-onlyの方針を維持し，証明不能や未対応を推測値として確定しない。

### 代数数・数体

- 実数／複素数`root[...]`を共通`AlgebraicNumber`へ統合し，最小多項式簡約，原始元，終結式，保証付きの根再同定を有界計算で実装した。
- 不変な`NumberFieldContext` / `AlgebraicElement`を導入し，選択した埋め込みとRational冪基底座標を保持する。同一数体の四則演算，小整数冪，exact逆元を数体内で処理する。
- 数体・逆元・最小多項式等の再利用キャッシュを追加し，最小多項式はexact Krylov消去で導出する。exact代数数の等値判定と実数上の大小比較も追加した。
- `root[...]`，exact Rational/複素Rational，`sqrt` / `cbrt`，`Phi`等を共通の代数数評価へ接続し，比較・定義域判定・Solveで同一のexact値として扱う。表示は従来どおり`root[minpoly,k]`を基本とする。

### Solver・意味論

- `solve[equation,Integer|Rational|Real|Complex]`を追加し，未知user symbolが一意な場合だけ変数を推定する。protected symbolはsolve変数にしない。
- `E^x` / `exp[x]`，`ln` / `log`等の表現差を正規化し，証明可能な正の底の指数方程式を対数で反転する。exact記号`lambertw`と実数枝`k=0/-1`の保証付き`N`も追加した。
- 非整数・非Rationalであることを`ValueFacts`へ伝播し，`N[SolutionSet]`は解集合構造を保ったまま数値化できる右辺だけを近似する。比較演算子の表示も読みやすく整えた。

### 保証付き数値評価・特殊函数

- `DecimalApproximation` / `ComplexDecimalApproximation`を第一級の保証付き数値として扱い，`CertifiedEnclosure`と`InformationEnclosure`を後続演算へ伝播する。
- `N[expr,p]`を有効10進桁数として統一し，全体を保証できない場合もHold属性を壊さず数値部分だけを近似する。重複診断も抑制した。
- `zeta` / `digamma` / `trigamma` / regularized `ibeta`を追加し，Gamma/Beta系の高精度計算を高速化した。値は存在するが保証付き計算未対応の場合は`CertifiedBackendUnsupported`としてDomainErrorから分離する。

### Exact FFT・Array表現

- exact Cyclotomic FFTを正式化し，対応する非2冪exact入力を`Q[t]/Phi_n(t)`上で計算する。5/7/10/12点等の往復を巨大な`cis[...]`式なしでexactに閉じ，適用不能時は従来のexact DFTへ戻す。
- `Expr::Node`をkind別typed nodeへ変更し，公開`Expr` APIを保ったまま固定payloadを削減した。
- dense `ArrayExpr`をimmutableなpaged packed storage + shape/offset/strides viewへ変更した。transposeはzero-copy viewとなり，矩形数値braceは`ArrayBuilder`へ直接構築する。

### 函数・補助機能・診断

- `isprime` / `nextprime` / `prevprime` / `factorint` / `totient`を`uint64`全域でdeterministic exactに追加した。証明できないBigIntをprobable-primeだけで`True`にはしない。
- bit演算群，`round[x,n]`，`fma`，`clamp`，`proj`，`range` / `table` / `map`，`explain`，周期的な実数解族を追加した。
- 構文解析に失敗した入力は`In[n]`を消費せず，`Infinity-Infinity`や`0*Infinity`等も誤簡約しない。未使用のMachine/double shimも削除した。

### 性能・検証

- `mmCal.Benchmarks --random-expressions`へseed+caseで再現可能な並列の意味論fuzzerを追加し，長時間実行時の履歴蓄積も防止した。
- 代数数体，特殊函数，exact Cyclotomic FFT等のベンチマークを追加し，採用・棄却した最適化の実測根拠を`docs/performance_optimization.*`へ記録した。利益のないキャッシュや無条件な高精度展開は採用しない。

### 文書・互換性

- README / Reference / Architecture / Roadmap / 性能文書をv1.5.3へ同期し，意図的な未実装・探索上限はroadmapとintentional-limits memorandumへ分離した。
- 内部原始元やcache表現は利用者向けFormatterへ露出せず，代数数の標準表示は`root[minpoly,k]`を維持する。

## v1.5.2 — 2026-08-13

### 構文・REPL（破壊的変更を含む）

- 函数呼び出しを`name[...]`へ一本化し，`name(...)`を廃止した。`()`はgrouping専用で，既知函数へ旧形式を使うとSyntaxErrorとなる。
- `x(x+1)`は暗黙乗算として受理し，Formatterは`x*(x+1)`へ正規化する。
- 履歴参照を`In[n]` / `Out[n]`へ統一した。正添字は絶対番号，負添字は相対参照，`0`は無効で，`@`群は`In[-n]`，`%`群は`Out[-n]`の短縮形である。

### 記号微積分・積分規則

- `sin[u]^m cos[u]^n`の非負整数冪を有限Fourier多項式へ還元し，負整数冪は`csc/sec`漸化式，`tan/cot/sec/csc`整数冪は各還元公式で処理する。異周波数積のproduct-to-sumも追加した。
- 共通引数を持つ`R(sin[x],cos[x])`へ有界なWeierstrass置換`t=tan[x/2]`を追加し，exact有理函数積分へ接続した。
- inverse-chain / substitution認識と二次根号族を拡張し，`2x(1+x^2)^5`，`x/(1+x^4)`等を構造的に処理する。
- derivative-back監査を維持しつつ，証明器不足だけで正しい原始函数を捨てない`ResolutionOnly`を明文化した。積分失敗も`unsupported` / `partial` / `conditionsRequired` / `noKnownClosedForm`へ分類した。
- `docs/memorandum/integralCatalog.md`を用い，三角・有理・特殊函数・分岐に敏感な積分族を横断監査した。

### 特殊函数

- `fresnelc` / `fresnels`，`hypergeometric1F1` / `hypergeometric2F1`，`ellipticF` / `ellipticE` / `ellipticPi`を追加し，exact退化，微分，保証付き実数`N`，対応する積分規則へ接続した。
- `Ei` / `Si` / `Ci` / `li` / `polylog`を追加し，代表的なexact値・微分・保証付き実数`N`・積分を実装した。
- `integrate[exp[x^n],x]`等では主値枝を保ちやすい1F1表現を優先する。
- 存在しない一般的な特殊函数の逆函数をSolverで捏造せず，安全なexact退化で既存Solverへ落ちる場合だけ解く。

### 精度対応 `N`・FFT

- exact Cyclotomic FFTを追加し，5点以上の対応する非2冪exact Rational入力を`Q[t]/Phi_n(t)`上で処理する。`ifft[fft[v]]`を巨大なroot-of-unity式なしでexactに閉じ，次数上限や数体所属を証明できない場合は従来のexact DFTへ戻す。
- `N[expr,p]`を精度対応の評価入口へ拡張し，対応函数へ`ApproximationContext`を伝播する。`N[fft[data],p]`は巨大exact式を作らずBigFloat/`ComplexInterval`の保証付きFFTへ直接進む。
- 保証付き非2冪FFTは小サイズで直接DFT，大サイズでBluesteinを使う。FFT/Matrix間で区間変換，10進化，保護桁の精密化等を共有した。
- 保証付き近似値は保証情報を保ったまま末尾0列を短く表示し，`N`由来の有限小数・整数にはexact値との区別用に最低1個の0を残す。

### Array・線形代数

- `{...}`を一般有限brace containerへ拡張し，矩形要素はrow-majorのdense `ArrayExpr`へ自動最適化した。0長次元と空Arrayのshapeも保持する。
- `dimensions` / `arrayRank` / `length` / `at` / `reshape`を整備し，`MatrixView` / `MatrixBuffer`で不要な複製を減らした。主要APIは`dot` / `matrixRank` / `norm` / `normalize`へ統一し，旧名は互換aliasとして残す。
- integer/Rational行列へBareiss fraction-free eliminationを導入し，`det` / `rref` / `matrixRank` / `inverse` / `solveLinear` / `nullSpace`で共有した。`solveLinear`は一意解だけを返し，`nullSpace`は決定的な基底と空shapeを保持する。
- `luDecomposition[A]`を追加し，exactでは行ピボット付き`P A = L U`，`N[...]`では保証付きinterval partial pivotingを使う。
- `qrDecomposition[A]`を矩形reduced Householder QRとして追加した。一般exact QRは式爆発のため3×3以下に制限し，上三角／上台形には高速経路を使う。自動block化は実測で優位性が安定せず不採用とした。
- exact `qrDecomposition[A]`の一般経路をExpr上のHouseholder展開からfraction-free直交化へ置換した。直交化中は平方根・除算を生成せず，primitive整数vectorとGCD content reductionで進め，full-rank caseはGram行列の対称Bareiss（fraction-free LDLᵀ相当）をfast pathに使う。従来の3×3 hard capを撤去し，rank落ちは直接fraction-free経路へfallbackする。
- exact Matrixのperformance-cliff再監査で，`inverse` / unique `solveLinear`のBareiss後段をgeneric Rational RREFから共通分母BigInt back-substitutionへ置換した。full-column-rank `rref`もRational backward phaseを省略する。modular inverseは32×32・256-bitまで再計測して引き続きBareissより遅く，自動dispatchは有効化しない。
- 実・複素reduced `svd`を追加し，`A^H A`を形成せずHouseholder bidiagonalization + one-sided Jacobiを使う。`conjugateTranspose`も追加した。
- `eigenvalues` / `eigenvectors` / `eigensystem`を追加した。簡単なexact行列はexactのまま，一般`N[...]`はComplex BigFloat Hessenberg + shifted QR + Schur経路で処理する。
- 記号`det` / `inverse`には三角行列の高速経路と展開上限を設け，rank/null-space等の不連続量は任意epsilonではなくexact pivotまたは区間で証明できる場合だけ確定する。

### 性能・安定性・開発基盤

- Bareissにより8～16次で`det`約5.4～11.4倍，`rref`約7.7～15.5倍の改善を実測した。`--matrix-large`も追加し，32/64次行列と1024次の保存・parse負荷を測定した。
- 1024×1024のRational行列ではExpr/Rational表現が先にメモリ・構文解析のボトルネックになることを確認し，型別`Expr::Node`，BigUInt/BigInt SBO，packed Array等を次版候補として記録した。
- MSVCで負の`RealInterval::point`が反転し得る評価順依存と`<stdexcept>`の推移的include依存を修正し，Matrix / FFTのfixed-seed invariantを拡充した。

### 文書・ライセンス

- README / Reference / Architecture / Roadmap / 性能文書をv1.5.2へ同期した。
- BSD 3-Clauseの著作権表示をプロジェクトmetadataと揃え，商標・ブランド利用は`TRADEMARKS.md` / `TRADEMARKS.ja.md`の別方針として明記した。

## v1.5.1 — 2026-08-12

v1.5.0のexact-first CAS基盤を維持しつつ，正規化，検証，多倍長整数，高精度数値評価，ベンチマーク基盤を重点的に改善した。

### 数式・CAS

- `Add`へAST全域の決定的な全順序を導入し，積・除算も定義性を保つ正規形へ整理した。`MathKnowledge`の非零知識と関係式の左右反転推論も強化した。
- Formatter/Parser 往復，積分derivative-back，Referenceとbuiltin registryの照合，極端値近似のテストを追加・拡充した。
- radix prefix衝突，Array隣接，負Rational表示を修正し，FFT plan/twiddleのtransform間cacheを追加した。

### 多倍長・高精度

- BigUInt乗算をschoolbook / Karatsuba / Toom-3の適応選択へ変更し，専用square，Burnikel–Ziegler除算，`2^k`除算，decimal divide-and-conquer変換を追加した。
- factorialのbalanced product treeを高速化し，`tryToUint64`の巨大値早期棄却，BigFloatの巨大指数差加減算の高速経路を追加した。
- `Pi`をbinary-splitting Chudnovskyへ，`exp/E`・`log`をbinary splitting + 保証付き範囲縮約へ変更し，巨大Radianの`sin/cos/tan`にも保証付き引数縮約を追加した。

### ベンチマーク・テスト

- Visual Studio solutionへ`mmCal.Benchmarks`を追加し，固定seedのrandom invariant，算法閾値，factorial，decimal I/O，高精度`Pi/exp/log`のベンチマークを常設化した。
- v1.5.1確定時点でinternal 1691 / 1691，ブラックボックス 1337 / 1337を確認した。

### 比較したが採用しなかったもの

- Prime-Swing factorial，binary GCD，Karatsuba workspace各種，Toom-3専用square／低閾値化，machine `fmod`による巨大trig縮約は実測上の利点が不足したため採用しなかった。

各判断の実測根拠は`docs/performance_optimization.ja.md`を参照。

---

## v1.5.0

旧版から数値モデル，Lexer/Parser/AST/Evaluator，Simplifier，Solver，CertifiedEvaluator，CLI，Formatter，テスト，文書をほぼ全面再構築し，mmCalをexact-firstなCLI電卓 / 小規模CASとして再定義した版である。
