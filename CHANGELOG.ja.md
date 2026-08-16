# Changelog

## Unreleased

- ﾏﾀﾞﾅｲﾖ

## v1.5.3 — 2026-08-16

v1.5.3は，v1.5.2で追加したsymbolic calculus / certified `N` / Array・線形代数の上に，**exact algebraic number field，Solverの意味論統一，高精度特殊函数，exact Cyclotomic FFT，representation最適化**をまとめて接続したreleaseである。exact-first・proof-onlyの方針は維持し，証明不能・backend未対応を推測結果として確定しない。正式化時点の検証値はinternal `2335 / 2335 PASS`，black-box `1715 / 1715 PASS`，Random Expression Fuzzer `100000 / 100000 PASS`（8 threads）。

### Algebraic number / number field

- Real/Complex `root[...]`を共通`AlgebraicNumber`へ統合し，bounded minimal-polynomial reduction，primitive-element reduction，resultant fallback，certified root re-identificationを実装。
- immutable `NumberFieldContext` / `AlgebraicElement`を導入し，選択embeddingとpower-basis Rational座標を演算列でpersistentに保持。同一fieldの`+ - * /`，小整数冪，exact inverseをfield内演算へ閉じる。
- 同一embedded generatorのbounded weak interning，primitive compositum/embedding reuse，reciprocal LRU，minimal-polynomial LRUを追加。minimal polynomialはincremental exact Krylov eliminationで導出する。
- exact algebraic equalityとReal orderingを追加。Complex Rootの決定的列挙順は数学的大小順へ流用しない。
- `root[...]`，exact Rational/complex Rational，`sqrt` / `cbrt`，`Phi`およびbounded exact arithmeticを`symbolic::exactAlgebraicValue`で接続し，表示形を強制変更せず比較・domain判定・Solveで同じ代数数として扱う。
- Root係数/indexの解釈とcanonical Root Expr生成を`algebraic_expression`境界へ集約し，Evaluator / Solver間の重複を解消。

### Solver / semantic coherence

- `solve[equation,Integer|Rational|Real|Complex]`を追加し，未知user symbolが一意の場合だけ安全に変数を推定する。protected symbolをsolve変数として誤受理しない。
- solve-safe normalizationを追加し，`E^x` / `exp[x]`，`ln` / `log`，`log2` / `log10`等の表現差をSolver能力差へしない。正baseのconstant-base exponentialを証明付きで対数反転する。
- exact symbolic `lambertw[z]` / `lambertw[k,z]`とReal `k=0/-1`のcertified `N`を追加。`a^x==x^2`の実解をbranch条件を証明できる範囲でLambert Wへ閉じる。
- `ValueFacts`へ`provablyNonInteger` / `provablyNonRational`を追加し，domain membershipとSolve constraintへ否定知識を伝播。
- `N[SolutionSet]`は解集合構造を保持したまま数値閉包なbinding右辺だけをcertifiedに近似する。
- Formatterは比較演算子`== != < <= > >=`だけ前後1空白で表示し，算術演算子のcompact形式は維持する。

### Certified numerical evaluation / special functions

- `DecimalApproximation` / `ComplexDecimalApproximation`をfirst-class certified numeric leafへ昇格し，真値保証用`CertifiedEnclosure`と情報量制限用`InformationEnclosure`を分離して後続演算へ伝播する。
- `N[expr,p]`を有効10進桁数として統一し，whole-expression certificationが閉じない場合はHold属性を壊さずnumeric subtreeだけをstructuralに近似する。inner diagnosticが既に原因を示した場合は外側`N`の重複WARNを抑止する。
- `zeta` / `digamma` / `trigamma` / regularized `ibeta`を追加。対応実領域のexact reduction・微分・certified `N`へ接続する。
- Gamma/Beta系を最適化。exact Rational引数を保持したrising product，static exact Bernoulli table，高精度BigInt remainder planner，`ibeta` point/shared-normalization経路により代表的な高精度負荷を大幅に短縮。
- Complex値自体は存在するがcertified backend未実装な場合を`CertifiedBackendUnsupported`としてDomainErrorから分離する。

### Exact FFT / Array representation

- exact Cyclotomic FFT backendを正式化。対応する非2冪exact入力を`Q[t]/Phi_n(t)`のRational power-basis座標で変換し，5/7/10/12点等の`ifft[fft[v]]`を巨大な`cis[...]`式なしでexactに閉じる。budget外・field membership未証明時は従来generic exact DFTへfallbackする。
- `Expr::Node`を巨大inline variantからkind別typed nodeへ変更し，公開`Expr` APIとstructural semanticsを維持したまま固定payloadを削減。
- dense `ArrayExpr`をimmutable paged packed backing + shape/offset/strides viewへ変更。transposeはzero-copy stride viewとなり，Lowererは矩形numeric braceを`ArrayBuilder`へ直接流す。

### Functions / utilities / diagnostics

- `isprime` / `nextprime` / `prevprime` / `factorint` / `totient`を`uint64`全域でdeterministic exactに追加。proof backend外のBigIntをprobable-primeから`True`へ昇格しない。
- `bitand` / `bitor` / `bitxor` / `bitnot` / bit shift/length/count/get，`round[x,n]`，`fma`，`clamp`，`proj`を追加。
- `range` / `table` / `map`，`explain`，parameterized periodic Real Solve familyを追加。
- Lexer/Parser/Lowererで成立しない入力は`In[n]`を消費しない。`Infinity-Infinity`，`0*Infinity`等は誤簡約せず保守的に未評価へ残す。
- Coreの未使用Machine/double shimを削除し，Exact / BigFloat / certified backendの境界を明確化。

### Performance / validation

- `mmCal.Benchmarks --random-expressions`を常設し，seed+case番号でthread数に依存せず再現できるcase-level parallel semantic fuzzerを追加。`resetForIndependentEvaluation()`で長時間運転時の履歴蓄積を防ぐ。
- `test_set/tester.py --timings`，`--algebraic-field`，`--special-functions`，`--exact-cyclotomic-fft`等のbenchmarkを追加し，採用・棄却した最適化を`docs/performance_optimization.*`へ実測根拠付きで記録。
- persistent multiplication-matrix cache，単純な高次数Bernoulli eager生成，高精度Stirling `K`の無条件拡大等は実測で利益が不足または退行したため採用しない。
- v1.5.3 release preparationでRoot表現helperを共通化し，constructor shadowingを除去。Visual Studio project XMLとGCC shadow warningをrelease checklistへ含める。

### Documentation / compatibility

- README / Reference / Architecture / Roadmap / performance文書をv1.5.3の実装状態へ同期。意図的未実装と探索budgetは`docs/roadmap.*`およびv1.5.3 intentional-limits memorandumへ分離して記録する。
- user-visibleなAlgebraic canonical outputは引き続き`root[minpoly,k]`を基本とし，内部primitive elementやcache representationをFormatterへ露出しない。

## v1.5.2 — 2026-08-13

v1.5.2は，v1.5.1までのexact-first数値基盤を維持しつつ，記号微積分・特殊函数・precision-aware評価・Array/線形代数をまとめて拡張したrelease。正式化時点でinternal `2027 / 2027 PASS`，black-box `1504 / 1504 PASS`。`mmCal.Benchmarks --random-only`のBigInt / special-function / Matrix / FFT fixed-seed invariantもPASS。

### 構文・REPL（breaking changeを含む）

- 函数呼び出しを`name[...]`へ一本化し，`name(...)`を廃止。`()`はgrouping専用
- 通常identifierの`x(x+1)`は暗黙乗算として受理し，Formatterは`x*(x+1)`へ正規化
- Parser/AST/Lowererから函数call delimiter分岐を削除し，既知函数へ旧`sin(x)`を使った場合はSyntaxError
- 履歴参照を`In[n]` / `Out[n]`へ整理。正添字は絶対番号，負添字は相対参照，`0`は無効
- `@` / `@@` / ... = `In[-1]` / `In[-2]` / ...，`%` / `%%` / ... = `Out[-1]` / `Out[-2]` / ...。`In[-n]`は入力slot，`Out[-n]`は成功出力を数える

### 記号微積分・積分Knowledge

- `sin[u]^m cos[u]^n`の非負整数冪を有限Fourier多項式へ還元する共通Knowledgeを追加し，高次数も個別ruleなしで処理
- `sin/cos`の負整数冪を`csc/sec` recurrenceへ還元し，`integrate[sin[2x]^(-2),x] -> -cot[2x]/2`等へ対応
- `tan/cot/sec/csc`整数冪のreduction formula，異周波数`sin/cos`積のproduct-to-sumを追加
- 共通引数を持つ有理式`R(sin[x],cos[x])`へbounded Weierstrass substitution `t=tan[x/2]`を追加し，変換後をexact rational integratorへ接続
- inverse-chain/substitution matcherを強化し，`2x(1+x^2)^5`，`3x^2 sqrt[1+x^3]`，`x/(1+x^4)`等を構造的に認識
- `sqrt[a-x^2]` / `sqrt[x^2+a]` / `sqrt[x^2-a]`等の具体係数二次根号familyを拡張
- 積分結果のderivative-back監査を維持し，証明器不足だけで正しいprimitiveを削らない`ResolutionOnly`方針を明文化
- 積分失敗diagnosticを`unsupported` / `partial` / `conditionsRequired` / `noKnownClosedForm`へ分類し，「未実装」と「閉形式がない」を混同しない
- `docs/memorandum/integralCatalog.md`の広範なcatalogを使って三角・有理・特殊函数・branch-sensitive familyを横断監査

### 特殊函数

- `fresnelc` / `fresnels`を追加。exact特殊値，奇対称，`D`，certified real `N`，二次位相積分へ接続
- `hypergeometric1F1[a,b,z]`を追加し，停止級数・安全なexact退化，`D`，certified real `N`を実装。`integrate[exp[x^n],x]`は原点を含めbranch-safeな1F1表現を優先
- `hypergeometric2F1[a,b,c,z]`を追加し，停止級数・`D`・certified real `N`とbinomial-power積分へ接続
- `ellipticF` / `ellipticE` / `ellipticPi`を追加。principal branch，amplitude derivative，certified real `N`，標準kernel積分を実装
- `Ei` / `Si` / `Ci` / `li` / `polylog`を追加し，代表的なexact退化，`D`，certified real `N`，積分Knowledgeを実装
- 一般特殊函数方程式に存在しないprincipal inverseを捏造せず，exact退化で既存Solverへ落ちる場合だけ解く

### precision-aware `N` / FFT

- Stage 7-7としてexact Cyclotomic FFT backendを追加。5点以上の非2冪exact Rational入力を，generic `cis[...]` Expr演算ではなくcyclotomic quotient `Q[t]/Phi_n(t)`のRational power-basis座標へ写して変換し，Gaussian Rational入力は必要に応じ`Q(zeta_lcm(n,4))`へexactに埋め込む。これにより5/7/10/12点等の`ifft[fft[v]]`が巨大なroot-of-unity式を残さずexactに`v`へ閉じる。内部generatorはroot isolationを行わず，既存表示と連続する単一`cis[-2 Pi/n Rad]`としてmaterializeする。degree budget 64を超える場合やsymbolic入力をfield座標へ証明付きで写せない場合は，Stage 7-7以前のgeneric exact DFTへfallbackする。旧dispatchは変更理由付きコメントとして隣接保存。専用`--exact-cyclotomic-fft` benchmarkを追加し，GCC Release/LTO-offで5/7/10/12点round-trip warmを約0.49/1.74/1.33/1.27 msで確認した。
- `N[expr,p]`を「完成したexact式を後からp桁化」するだけでなく，対応builtinへ`ApproximationContext`を伝播するprecision-aware評価入口へ拡張
- `N[fft[data],p]`は巨大exact Fourier式を先に構築せず，BigFloat/`ComplexInterval`上のcertified radix-2 backendへ直接dispatch
- certified非2冪FFTは小サイズdirect DFT，大サイズBluesteinへ分岐。現benchmark環境のpolicy thresholdは約96点
- FFT/Matrixでexpression→interval，decimalization，guard-digit refinement等の共通approximation helperを共有
- certified fixed-digit結果の表示は保証metadataを保持したまま冗長な末尾0列を圧縮（`1.000... -> 1.0`，`1.500... -> 1.50`）。exact有限小数は`N[1/2,10] -> 0.5`のまま

### Array・線形代数

- `{...}`を一般有限brace containerへ拡張。矩形childはshape + row-major flat storageのdense `ArrayExpr`へ自動最適化し，shapeが異なる`{Q,R}` / `{U,S,V}`等は一般braceとして保持
- zero-length dimensionを正式保持し，braceだけでshapeを復元できない空Arrayを`reshape[{}, {...}]`でround-trip
- `dimensions` / `arrayRank` / `length` / prefix対応`at` / `reshape`を整備。indexは0-based
- `MatrixView` / `MatrixBuffer`を追加し，row-major Arrayを不要にnested vectorへ複製しないalgorithm基盤へ整理
- canonical APIを`dot` / `matrixRank` / `norm` / `normalize`へ統一し，旧`matmul/mmul/vdot/rank/mrank/vnorm/vnormalize/mget`は互換aliasとして維持
- integer/Rational行列へBareiss fraction-free eliminationを導入。行別分母LCMで整数liftし，`det` / `rref` / `matrixRank` / `inverse` / `solveLinear` / `nullSpace`でkernelを共有
- `solveLinear[A,b]`は一意解だけを返し，full column rankの整合した過剰決定系にも対応。不整合・自由変数系はDomain error
- `nullSpace[A]`はfree column昇順のcanonical basisを返し，full column rankでもshape `{0,n}`を保持
- `luDecomposition[A]`を追加。exactはrow-pivoted `P A = L U`，`N[...]`はcertified interval partial pivotingへ直接dispatch
- `qrDecomposition[A]`を矩形reduced Householder QRとして追加。一般exact QRは式爆発実測により3×3以下へ制限し，上三角/上台形caseはfast path
- Householder column-block kernelを試作・計測したがblock=1/8/16/32で一貫した勝者がなく，自動block化は不採用。benchmarkは保持
- `svd` / `singularValueDecomposition`を実・複素reduced SVDとして追加。`A^H A`を形成せずHouseholder bidiagonalization + one-sided Jacobiを採用し，再構成・直交性をinterval監査
- `conjugateTranspose`を追加し，複素SVD等のHermitian relationを共通表現
- `eigenvalues` / `eigenvectors` / `eigensystem`を追加。exactは三角/対角/distinct-root exact Number 2×2，一般`N[...]`はComplex BigFloat Hessenberg + implicit shifted QR + Schur backendへ直接dispatch
- general symbolic `det` / `inverse`へ三角fast pathと展開budgetを追加し，dense高次の階乗級expression explosionを抑止
- rank/null-spaceのような不連続量はepsilon判定を導入せず，exact pivot structureまたはintervalで証明できる場合だけ結果を確定

### 性能・安定性・開発基盤

- Bareiss導入によりStage 2 Gaussian比で8–16次`det`約5.4–11.4倍，`rref`約7.7–15.5倍の改善を実測
- `mmCal.Benchmarks --matrix-large`を追加し，32/64次の主要Matrix函数と1024次storage/parse負荷を測定
- 1024×1024の10桁Rational Exprを直接構築すると最大RSS約0.69 GB，generator形式textのparse + `dimensions`だけで約10.9 s / 1.99 GBとなり，large dense MatrixではExpr/Rational storageが先にbottleneckになることを確認
- 次版ToDoとして`Expr::Node`巨大variantのtyped-node化，BigUInt/BigInt SBO，numeric Array/approximate Matrix packed storage，巨大brace parse allocation削減を記録
- `RealInterval::point`で同一`BigFloat`を一初期化式内でcopy/moveしていた評価順依存を除去し，MSVCで負point intervalが反転し得るSVD/Eigen failureを修正。負値point専用回帰を追加
- `expression_interval.cpp`が直接使う`std::overflow_error`の宣言元として`<stdexcept>`を明示includeし，transitive include依存を除去
- fixed-seed Matrix invariantをBareiss / inverse / solve / nullSpace / LU / QR / real+complex SVD / Eigenまで拡張し，FFT random round-tripと併せてbenchmark projectへ常設

### 文書・ライセンス

- README / Reference / Architecture / Roadmap / performance docsをv1.5.2正式版へ同期
- BSD 3-Clauseの著作権表示をproject metadataと整合させ，商標・ブランド利用は`TRADEMARKS.md` / `TRADEMARKS.ja.md`の別ポリシーで扱うことを明記

## v1.5.1 — 2026-08-12

v1.5.0のexact-first CAS基盤を維持しつつ，canonicalization，検証，巨大整数，高精度数値評価，benchmark基盤を重点的に改善した。

### 数式・CAS

- `Add`へAST全域の決定的strict total orderingを導入
- 積・除算をdefinednessを保つcanonical normal formへ整理
- `MathKnowledge`の非零知識と関係式左右反転推論を強化
- formatter/parser生成型round-trip testを追加し，radix prefix衝突，Array隣接，負Rationalの表記不安定を修正
- 積分rule/familyを横断するderivative-back harnessを追加
- Referenceとbuiltin registryの自動照合を追加
- approximation extreme-value testを拡充
- FFT plan/twiddleのtransform間cacheを追加

### 多倍長・高精度

- BigUInt乗算をschoolbook / Karatsuba / Toom-3の適応dispatchへ変更
- 専用squareを追加
- factorialはbalanced product treeを維持し，leaf構築・1-limb経路を高速化
- Burnikel–Ziegler divisionと`2^k`除算fast pathを追加
- decimal parse/toStringを`10^9` chunk + divide-and-conquer化
- `tryToUint64`の巨大値早期棄却を追加
- BigFloat extreme exponent-gap加減算をdirected roundingを保ったままfast-path化
- `Pi`をbinary-splitting Chudnovskyへ変更
- `exp/E`と`log`をbinary splitting + certified range reductionへ変更
- 巨大Radianの`sin/cos/tan`へcertified argument reductionを追加

### Benchmark / test

- Visual Studio solutionへ`mmCal.Benchmarks`を追加
- fixed-seed random invariant，threshold sweep，factorial，decimal I/O，高精度`Pi/exp/log`のbenchmarkを常設化
- v1.5.1確定時点: internal 1691 / 1691，black-box 1337 / 1337

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

旧版から数値モデル，Lexer/Parser/AST/Evaluator，Simplifier，Solver，CertifiedEvaluator，CLI，formatter，tests，documentationをほぼ全面再構築し，mmCalをexact-first CLI calculator / compact CASとして再定義したrelease。
