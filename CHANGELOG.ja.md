# Changelog

## v1.5.5 — 2026-09-05

### Series・解析

- `SeriesData` / TPSAを核とする`series[...]`を追加し，Taylor / Laurent / Puiseux / logarithmic展開を`D` / `integrate` / `limit` / `toNormal`へ接続した。`series[...,{x,Infinity,n}]`にも対応する。
- 初等函数に加え，`Ei` / `Si` / `Ci` / `li`，Fresnel函数，Lambert W，Gamma系，`polylog`等の局所展開を拡張した。
- Seriesのmetadata，Infinity展開，nested frontend合成，Cases境界等で発生していた未評価・誤評価を修正した。

### Array・Vector

- `inner` / `outer` / `distance` / `projection` / `rejection`，`gramSchmidt`，直交・一次独立判定等を追加し，複素vectorではHermitian内積へ意味論を統一した。
- `grad` / `divergence` / `curl` / `laplacian` / `jacobian` / `hessian` / `directionalDerivative`を追加した。
- Vectorの旧convenience名を整理し，Array演算と正式APIへ統合した。`at`もfinite `SolutionSet`のbranch / binding取得へ拡張した。

### 記号計算・Solve・積分

- `solve`のdefinedness，principal inverse，Real domain，relation array，radical / transcendental dispatch等を横断的に修正した。
- `limit` / `D` / `integrate` / `series` / `solve` / Vector Calculus / 式変形frontend間のheld binder処理を整理し，組合せた場合だけ未評価・構文破壊になる問題を多数修正した。
- 高次有理函数積分，代数的log積分，高階微分，多項式×初等函数積分，三角冪積分等を高速化・一般化した。
- `factor[x^257-1]`を代表とするperfect-power探索の性能崖を有限体precheckで解消した。

### 数値近似・表示

- 数値近似のprovenanceを`ExactValue` / `CertifiedInterval` / `VerifiedApproximation`へ分離し，SVD / eigen等の検証済み数値をrigorous point enclosureとして誤用しないようにした。
- `N`の部分数値化で指数，branch番号，離散parameter等のexact integer構造を保持するようにした。
- finite-precision値の符号・括弧・多項式順序をformatterで正規化し，`+-` / `--`，複素係数のprecedence崩れ，不要な`+0.0I`等を修正した。exact-source integer pointはtop-level値表示で不要な`.0`を省略する。
- certified特殊函数評価の重複計算，境界dispatch，guard再試行等を整理し，Ci，2F1，Gamma系，elliptic函数等の性能と安定性を改善した。

### 制限撤廃・内部整理

- 古い固定次数・指数・元数制限のうち，共通`EvaluationBudget`や構造fast pathで安全に置換できるものを撤廃した。高階`D`，binomial solve，Series整数冪，三角冪・指数函数積分，有理函数積分等が対象である。
- exact scalar，relation，Cases，linear algebra，Builtin family判定等の重複実装を共通化し，certified retry / fallback境界を整理した。
- REPLへ`:quit` / `:exit`と`:layout`を追加し，exact FFTのcyclotomic再埋込みも拡張した。

## v1.5.4 — 2026-08-31

### Solver・記号計算

- 多変数polynomial `solve`を強化した。Gröbner消去と再帰的specializationを拡張し，対応可能なzero-dimensional系をexactに完全列挙するほか，`x y==0`や`x^2+y^2==1`等のpositive-dimensional系も，証明可能な範囲で自由パラメータとexact条件を持つ解集合として返す。一般多様体の不完全なparameterizationは引き続き推測しない。
- Real `solve`を拡張し，`u exp[u]=a`，`exp[p x+q]==c x+d`，正定数底の指数方程式，`x^x==r`の一部，`lambertw[u]==r`，`cosh[u]==r`等をLambert Wや既存inverseへ接続した。値域・単調性・凸性等による不存在／一意性証明も追加し，証明できた場合だけ空集合またはexact解を返す。
- `abs`および主値`radical`方程式をReal solverへ追加した。平方等で導入される偽根をexactに除外し，Complex上の円周等を有限個の解へ誤って潰さない。
- 数学的場合分け`cases[value if condition; ...]`を追加し，`simplify`，`N`，`D`，`integrate`，`limit`へ接続した。未知条件は保持し，false branchは評価しない。
- `groebnerBasis[...]` / `polynomialReduce[...]`を追加した。Lex / GrLex / GrevLexのexact Rational Gröbner基底を扱い，非線形多項式`solve`にも利用する。
- `simplify` / `fullSimplify`を定義性に配慮するよう強化し，`F/F -> 1`，`F^0 -> 1`，特殊函数の退化等は必要な定義条件を証明できる場合だけ適用する。exact複素数は演算後に虚部が0なら実数へ正規化する。

### 微分・積分・極限

- exact有理函数積分を一般化し，Hermite reductionとalgebraic `Root`を用いて高次・重複分母を含む有理函数まで扱える範囲を拡張した。既存の簡潔な`Log/atan/atanh`表現や初等的置換を優先する。
- 仮定付き定積分・広義積分を拡張し，Gamma/Beta/Mellin，Frullani型，対数moment，`1+x^q`型，Dirichlet/Fresnel等のexact条件付き評価を追加した。収束条件や内部特異点を証明できない場合は未評価を保持する。
- `D`の`cases[...]`処理と高階微分を強化した。変数依存の境界では非微分可能点へ偽の値を割り当てず，Lambert W，`polylog`，二次式指数函数等の高階微分をよりcompactなexact式で構成する。
- `limit[expr,{x,a,direction}]`を追加し，`Ei` / `Ci` / `li`や実函数の既知極限，周期振動，squeeze可能な極限を補強した。積分のderivative-back検証も共通定義域上の恒等性を確認するよう修正した。

### Algebraic Root・高次多項式

- 高次Complex `Root`と一般高次polynomial `solve`の性能段差を大幅に改善した。root isolation / refinementと既約性判定の重複計算を減らし，対称多項式やroot半径が大きく異なる場合の停滞も解消した。96次までの監査範囲を拡張している。
- AlgebraicNumber / number field内の演算を高速化し，中間`Root`生成を抑えてexact性を保ったまま高次冪等の負荷を低減した。

### 保証付き `N`・特殊函数

- 複素数を含む保証付き`N`を大幅に拡張した。`erf/erfc`，`Ei/Si/Ci`，Fresnel函数，`1F1`，`2F1`，`polylog`，`zeta`，`gamma`，`digamma/trigamma`，Lambert W，楕円積分等で主値枝と誤差保証を維持した評価範囲を広げた。
- 従来の固定的な大きさ境界を多数撤去し，実数`Ei/Si/Ci`，複素`Ei/Ci`，Fresnel，`1F1`，`2F1`，`polylog`，`ellipticF/E/Pi`等を収束性と資源上限に基づいて選択的に評価するようにした。`zeta`は`s=1`を除く有限複素平面へ，`ibeta`はfinite-precisionを含む正実パラメータへ対応範囲を拡張した。
- Lambert Wは`-1/e`近傍と任意整数branchの保証付き複素評価を追加した。finite-precision入力で分岐側を証明できない場合は推測せず`N::precision`を返す。
- `CertifiedEnclosure`と`InformationEnclosure`の役割を整理し，有限precision入力の隠れた保護桁から過剰な精度を捏造しないよう統一した。特異点・分岐切断・算法境界の曖昧性は，数学的なDomainErrorと区別して`N::precision`または`N::unsupported`として扱う。
- `N`の表示を簡潔化し，近零値や末尾0を整理しつつapproximation provenanceを保持する。nested `N`，`diff[...,digits]`，`nintegrate[...,digits]`も共通のprecision制約と情報量制約へ統一した。
- `acosh`の分岐切断近傍，`expm1/log1p`の0近傍，`2F1`の接続公式，複素`li`等の相殺・過剰refinementによる性能問題を修正した。

### Exact線形代数

- `conditionNumber`，`pseudoInverse`，`leastSquares`を追加した。exact行列ではrankをexactに判定し，外側`N[...]`では保証付きSVDを利用する。
- exact Integer/Rational行列の`det` / `solveLinear`へmodular + CRT経路を追加し，検証付きでBareissと自動選択する。`leastSquares`の中間丸めによる精度損失も修正した。

### 評価資源・CLI

- 評価ごとの共通`EvaluationBudget` / `EvaluationLimits`を追加した。巨大整数，Array/Matrix，Solver，積分，AlgebraicNumber，要求precision等を計算前または計算中に制限し，超過時は`ResourceLimitError`として報告する。
- Ctrl-C / Ctrl-Break / SIGINTによる協調的な評価取消しを追加した。
- 行単位自動化向け`--batch`を追加した。REPLの`:help`も全callable函数，定数，入力規則，例，近似候補を含む一覧へ拡張した。

### 検証・互換性

- public CLIだけを使うblack-box監査とRandom Expression Fuzzerを拡張し，`D` / `integrate` / `limit` / `N` / `cases` / Gröbner / Solve / exact FFT / 線形代数 / precision provenance等の性質検証を追加した。
- Formatter，履歴参照，`cases[...]`を跨ぐ微積分，主値分岐，可除特異点，DomainError分類等の回帰を修正した。公開構文とexact-first方針は維持する。

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
- `D` / `integrate` / `limit`等のheld symbolic operator内でも`In[n]` / `Out[n]` / `%`を履歴参照として再帰的に解決し，`integrate[Out[1],x]`等が未知函数扱いになる不整合を修正した。
- 多変数polynomialの正次元代替経路に定数係数の一次変数消去を追加し，複数方程式を低次元systemへexactに落として既存のRealパラメータ-定義域 射影へ接続した。

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
- exact `qrDecomposition[A]`の一般経路をExpr上のHouseholder展開からfraction-free直交化へ置換した。直交化中は平方根・除算を生成せず，primitive整数vectorとGCD content reductionで進め，full-rank caseはGram行列の対称Bareiss（fraction-free LDLᵀ相当）をfast pathに使う。従来の3×3 hard capを撤去し，rank落ちは直接fraction-free経路へ切り替える。
- exact Matrixのperformance-cliff再監査で，`inverse` / unique `solveLinear`のBareiss後段をgeneric Rational RREFから共通分母BigInt back-substitutionへ置換した。full-column-rank `rref`もRational backward phaseを省略する。modular inverseは32×32・256-bitまで再計測して引き続きBareissより遅く，自動dispatchは有効化しない。
- 実・複素reduced `svd`を追加し，`A^H A`を形成せずHouseholder bidiagonalization + one-sided Jacobiを使う。`conjugateTranspose`も追加した。
- `eigenvalues` / `eigenvectors` / `eigensystem`を追加した。簡単なexact行列はexactのまま，一般`N[...]`はComplex BigFloat Hessenberg + shifted QR + Schur経路で処理する。
- 記号`det` / `inverse`には三角行列の高速経路と展開上限を設け，rank/null-space等の不連続量は任意epsilonではなくexact pivotまたは区間で証明できる場合だけ確定する。

### 性能・安定性・開発基盤

- internal test runnerは`--timings`指定時だけsuite別実時間を表示し，通常実行の出力契約は維持する。
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
