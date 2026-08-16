# Changelog

## Unreleased

### 開発・検証基盤

- Visual Studioの`mmCal.Core.vcxproj` / `.vcxproj.filters`で欠落していた`</ItemGroup>`を修復し，全`.vcxproj` / `.filters`をXML parserで検証。プロジェクト読込不能を解消
- AlgebraicNumberのbounded canonicalizationを拡張。次数16以下で証明できる場合，選択Rootをexact有理既約因子＝minimal polynomialへ縮約し，実RootはSturm根数，Complex Rootはcertified isolating diskで因子所属を証明する。さらにoperand minimal polynomialの既約性と積次数simple extensionを証明できる場合は`theta=alpha+c beta`のprimitive-element reductionを使い，power basis上のexact線形従属から結果minimal polynomialを導出。証明失敗時はresultantへfallbackする。また検証中に既存resultant arithmeticの共役root誤選択を発見し，root enumerationとindividual canonicalizationを分離して修正
- Core内で未使用だった`machine_math.cpp/.hpp`と専用`elementary_function.hpp`を削除。現在の`src/`実装はMachine/double数値計算へ依存せず，Exact / BigFloat / certified backendだけで閉じる。将来`plot`等でMachine evaluatorが必要になった場合は，既存Exact/Certified経路へ混入させず明示的な別backendとして再導入する
- BigInt向けbit utilityとして`bitand` / `bitor` / `bitxor` / `bitnot` / `bitshiftl` / `bitshiftr` / `bitlength` / `bitcount` / `bitget`を追加。負数を含むbit演算は無限長2の補数semanticsへ統一し，公開右shiftはarithmetic shiftとする。`round[x,n]`はnearest-evenの10進桁丸め（負の`n`も可），`fma`は近似値で中間Decimal丸めを挟まない一括certified評価，`clamp`はInformationEnclosureで境界を証明できる場合だけ確定，`proj`は現在の有限複素値に対する射影恒等写像として追加
- `root[poly,k,Complex]`と`ComplexAlgebraicNumber`を追加し，square-free monic Rational多項式の全複素根をcertified isolating diskへ分離。高次Rational係数多項式のComplex/default `solve`もexact Complex Rootへfallbackする。Root / exact Rational / exact complex Rational間の`+ - * /`と小整数冪はresultant候補多項式とoperand/result isolating regionの再同定でbounded exact algebraic-field演算へ接続し，候補次数budget 16を超える場合や結果rootを一意に証明できない場合はsymbolic式を保持。後続のbounded canonicalizationでminimal-polynomial factor reductionとprimitive-element reductionまで実装し，一般algebraic equality・ordering，persistent number-field，budget超過の完全Q因子分解は引き続き意図的に保留
- 特殊函数へ`zeta` / `digamma` / `trigamma` / 正則化不完全Beta `ibeta`を追加。`zeta`は代表exact値・自明零点と実軸`s>1`のEuler-Maclaurin certified backend，`digamma/trigamma`は正実数のBernoulli漸近＋recurrence，`ibeta`は`a,b>0, 0<=x<=1`のreal contractとexact integer-parameter有限和・exact Rational parameterのcertified 2F1/Beta経路を実装。`D[gamma]` / `D[lgamma]` / `D[digamma]`も新しい函数語彙へ接続
- 軽量数論として`isprime` / `nextprime` / `prevprime` / `factorint` / `totient`を追加。`uint64`全域では証明範囲が既知のdeterministic strong Miller-Rabinを使い，因数分解はprime verification付きPollard-Rhoでexactに返す。現在の証明backendを超えるBigIntはprobable-primeを`True`として返さず未評価を保持
- `docs/roadmap.md`を正式化し，意図的未実装，backend未対応，探索budgetを区別して整理。Complex Root/algebraic-field，特殊函数の複素解析接続，一般parametric Solve，arbitrary-BigInt数論，Machine evaluator，BigUInt SBO，Core全体multi-thread化等の保留理由を明文化
- Lexer / Parser / Lowererで式として成立しなかった入力は履歴へcommitせず，`In[n]`番号を消費しないよう変更。SyntaxError表示はpending番号を示し，修正後の次入力は同じ`In[n]`を再利用する。parse/lower成功後の評価エラーは従来どおり入力履歴へ残す
- `mmCal.Benchmarks --random-expressions`へgrammar-awareなsemantic expression fuzzerを追加
- `--loop`ではcase数を制限せず連続実行し，最初のFAILを検出した時点でshrinking・seed/case再現情報を表示して即停止する
- 各caseはmaster seedと1-based case番号から独立生成され，`--seed N --case M`だけで該当caseを直接再現できる
- 生成depthは浅い式を主体にしつつ稀に深い式を混ぜる重み付き分布とし，既定`--max-depth 16`の範囲でcaseごとに変化する
- 初期invariantとしてexact formatter round-trip，`fullSimplify`値保存，`expand`/`factor`の多項式値保存，double transpose，`det[A]==det[transpose[A]]`を検証
- generatorは合法なexact算術・多項式・小行列を主体とし，ランダムなTypeErrorをFAILとして量産しないsemantic fuzzingを優先
- `--threads N`を追加し，workerごとに独立`KernelSession`を持つcase-level parallel fuzzingへ対応。master seed + case番号による再現性はthread数に依存しない
- `--nostop-loop`を追加。通常`--loop`は最初のFAILで停止する一方，`--nostop-loop`はFAILをshrinking・表示した後も継続する
- fuzzerのcase間初期化を`KernelSession::resetForIndependentEvaluation()`へ分離し，定義・履歴・入力番号・diagnosticを破棄しつつRNG streamと角度設定を保持。長時間loopでsession履歴を蓄積せず，`reset()`のentropy reseedも回避
- `Expr::Node`の最大payload依存`std::variant`を廃止し，`ExprKind` header + kind別typed payloadへ単独refactor。公開`Expr` API，node identity，structural equalityは維持
- 同一GCC Release/LTO-offの`--matrix-large transpose 1024 16`で最大RSSを旧variant版`693312 KiB`（約677.1 MiB）からtyped-node版`299668 KiB`（約292.6 MiB）へ削減し，約56.8% / 384.4 MiB削減を確認
- dense `ArrayExpr`の永続表現を固定1024要素pageのimmutable shared backingへ変更し，Integer / Rational / Number / DecimalApproximation / ComplexDecimalApproximation / Genericをpage単位でpacked保持。Array全体は従来どおり単一の`ExprKind::Array`として扱う
- `ArrayBuilder`はpromotionを現在page内だけへ限定し，完成済みpageを再構築しない。100万要素の末尾にsymbolic値が現れてもGeneric化するのは最後の最大1024要素だけで，既存packed pageは共有したまま保持
- `ArrayExpr`へshape / offset / stridesのimmutable viewを導入し，`transpose`はbackingを複製せずshape/stride交換だけで実行。contiguous `reshape`とprefix sliceも可能な範囲でbackingを共有
- 矩形brace literalはLowererから単一`ArrayBuilder`へleafを直接流し，numeric literalでは要素ごとの`Expr` nodeを先に構築しない。ragged braceは従来どおり一般`ListExpr`へlowering
- naiveな単一`vector<Rational>` packed案は，transposeで100万個のRational/BigIntをdeep copyして約650–675 msとなるため棄却。shared paged backing + view方式では`--matrix-large transpose 1024 16`が約0.059 ms，最大RSS約133.4 MiB。13.63 MBの1024×1024 decimal literalをCLIで`dimensions[...]`した測定は約4.55 s / 431 MiB
- BigUInt / BigInt SBOは今回混ぜず，保守性と単独benchmark可能性を優先して見送り
- `explain[value]`を追加。評価済み値が既に保持する安価なmetadataだけを`{{"Property", value}, ...}`形式のbrace値で返し，追加の数学計算やArray全走査は行わない。Arrayはshape / element count / square判定等，certified近似はexact Rational enclosureを公開し，`explain[value,"internal"]`では表現・storage等の非互換保証debug metadataも返す。`Pi/E/Phi`はMathRegistryの定数metadataを，`Infinity`等の予約symbolはSymbolRegistryの既知意味を直接参照し，通常の未定義symbolと区別して説明する
- `DecimalApproximation` / `ComplexDecimalApproximation`を通常の`+ - * /`と単項`-`へ再投入可能にし，exact `Number`をpoint intervalとして混在できるようにした。近似値metadataを真値保証用`CertifiedEnclosure`と情報量制限用`InformationEnclosure`へ分離し，`CertifiedEnclosure ⊆ InformationEnclosure`を不変条件として両方を独立伝播する。`accuracy` / `precision` / tolerance省略時`rationalize`はInformationEnclosureを基準にし，`N[N[Pi,20],100]`や後続演算で内部guard桁・hidden exact pointを新しい情報として回収しない
- `N[expr,p]`の第2引数を固定小数点以下桁数ではなく**有効10進桁数(significant decimal digits)**として定義。値のscaleに追従するInformationEnclosureを生成し，`N[Pi*10^20,20]`や微小値でも相対Precisionを保つ。固定小数表示は`:fix` / `--fix`へ明確に分離し，0近傍では相対Precisionが0でもabsolute Accuracyを保持するzero-centered certified approximationを返せるようにした
- `DecimalApproximation` / `ComplexDecimalApproximation`をcertified対応scalar函数のfirst-class入力へ拡張。`sin` / `exp` / `log` / `sqrt` / 双曲線・逆函数 / `gamma` / `erf` / `Ei` / `Si` / `Ci`等でCertifiedEnclosureとInformationEnclosureを独立伝播し，primitiveへrewriteされる`log2` / `log10` / `fract`も同経路へ再投入する。ordered comparisonと`min` / `max`はInformationEnclosureだけで結論を証明できる場合に限って確定し，hidden guard桁をBoolean結果から漏らさない。exact Rational parameterのみを受ける一部特殊函数backendは無理に近似parameter対応せず未評価を維持
- `Infinity`を通常の有限symbolと同じ相殺・0積簡約へ流さないよう修正。`Infinity-Infinity`，`0*Infinity`，`Infinity/Infinity`等の未定義・未仕様形は現段階では保守的に未評価へ残し，`0`等を捏造しない
- InformationEnclosureを近似値の順序比較だけでなく`==` / `!=`にも接続。実部または虚部のInformationEnclosureが互いに素なら不等を確定し，双方が同一pointと証明できる場合だけ等値を確定する。重なる区間はhidden guard桁を使わず未評価へ残す。`floor` / `ceil` / `trunc` / `round`および`min` / `max`も既存どおりInformationEnclosureのみで離散結果を証明する
- exact finite sequence / iteration APIとして`range` / `table` / `map`を追加。`range`はexact Integer / Rational stepを用い，`table`はheld bodyとlocal iterator scopeでnested反復・外側定義の復元を保証する。`map`はArray / ragged braceのscalar leafへ函数を明示適用し，通常函数を暗黙element-wise化しない。`explain[builtin]`はBuiltinRegistry / MathRegistryからarity・held属性・domain・parity・period・principal inverse・real range等をO(1)で照会可能にした
- Real上の`sin/cos/tan`周期方程式へinteger formal parameter付きSolutionSetを追加。引数がsolve変数に対するexact非零一次係数のaffine式である範囲を安全に扱い，`solve[sin[x]==0,x,Real] -> {x==Pi k where k in Integer}`等を返す。target range，session角度mode，endpoint branch重複除去，parameter name衝突を監査し，非線形argumentやComplex全解は未解決のまま保持する
- Gaussian積分のpreferred canonical formを改善。exact positive Rational `a`に対する`integrate[exp[-a x^2],x]`を一般1F1より`erf`へ優先し，有限積/商の安全なInfinity極限を追加して`integrate[exp[-x^2],{x,0,Infinity}] -> sqrt[Pi]/2`および全実軸Gaussian `-> sqrt[Pi]`までexactに閉じる。`erf[0]` / `erfc[0]`もSimplifierでcanonical化
- exact実代数根`root[{a0,...,an},k]`と`RealAlgebraicNumber` specialized IRを追加。係数は昇冪exact Rational，`k`は異なる実根を昇順に数える1-based index。定義多項式をmonic square-freeへ正規化し，Rational Sturm列で根数をexactに数え，Rational isolating intervalを保持して`N[root[...,k],p]`をcertified refinementする。既存Solverで閉じないRational係数高次多項式はまずReal領域のRoot fallbackへ接続した。後続のUnreleased変更でComplex Root isolationとbounded algebraic-field演算まで拡張し，後続のbounded実装でminimal-polynomial factor reduction / primitive-element reductionまで接続
- nested positive exact integer Powerを通常Simplifierで`(a^m)^n -> a^(mn)`へ安全に正規化し，random fuzzerが検出した`expand`/`factor`不変量違反を修正
- Formatterは`^`の右結合性を明示し，左nested Powerを`(a^b)^c`と括弧付きで出力してASTの意味を保持
- exact Rational定数を連続減算する`(a-b)-c`を`a-(b+c)`へ畳み，`((((x-1)-3)^3)^4...) -> (x-4)^144`のcanonicalizationを改善
- 正のexact integer冪`(c*a)^n`ではexact numeric係数`c`だけを安全に冪乗して外へ出し，`(861(1-x)^64)^3 -> 638277381(1-x)^192`までcanonical化。一般複素指数への積の冪分配は行わない

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
