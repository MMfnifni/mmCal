# mmCal v1.5.3 意図的な制限・探索budget・実装閾値一覧

> 対象: v1.5.3正式化ソース
>
> 本文では「数学上の制限」「安全弁」「探索budget」「実装閾値」を区別する。閾値は入力を拒否する上限ではない。

## 1. Random Expression Fuzzer

| 項目                       |         現在値 | 意味                                              |
| -------------------------- | -------------: | ------------------------------------------------- |
| `maxDepth`既定値           |             16 | hard limitではない。`--max-depth N`で変更可能     |
| `--max-depth` CLI hard max |           なし | 1以上だけ検査。上げすぎるとAST/式サイズが急増する |
| target depth 1–4           |            60% | 通常case                                          |
| target depth 5–7           |            25% | 中程度                                            |
| target depth 8–10          |            10% | 深め                                              |
| target depth 11–13         |             4% | stress                                            |
| target depth 14–max        |             1% | rare stress                                       |
| `expand`用generator depth  |          最大6 | 多項式爆発を避けるfuzzer側制限                    |
| `factor`用generator depth  |          最大5 | 同上                                              |
| scalar exponent            |           1..5 | generator側                                       |
| polynomial exponent        |           1..4 | generator側                                       |
| Matrix shape               |    1..4 × 1..4 | generator側                                       |
| Matrix integer entry       |          -5..5 | generator側                                       |
| report interval            |          10000 | `--report-every`で変更可能                        |
| finite cases               |          10000 | `--cases`で変更可能                               |
| `--loop`                   | case数上限なし | 最初のFAILで停止                                  |

### loopのsession状態

v1.5.3では，各workerが独立`KernelSession`を所有し，case間・比較評価間に`resetForIndependentEvaluation()`を呼ぶ。これは`history_` / `inputHistory_` / `outputHistory_` / diagnostic / 定義 / `inputCount_`を初期化する一方，角度設定とRNG streamを保持する。通常の`reset()`は同じtransient resetの後にentropy reseedするため，fuzzerだけが不必要なreseedingを避けられる。

したがって，旧実装で問題だった「`clearHistory()`しても`inputCount_`が増え続け，次回`evaluate()`が巨大vectorへresizeする」経路は解消済みである。参考測定では同一Release/LTO-off build，`--threads 8 --seed 1234`で100,000 caseの最大RSSが`48136 KiB`，200,000 caseが`51364 KiB`であり，case数に比例する履歴蓄積は見られなかった。overnight burn-inではallocator high-water等を含め，process RSSの定常性を引き続き監視する。

現generatorはrandom builtinを生成しない。将来これを追加する場合は，case番号からrandom builtin用seedも明示的に導出し，thread schedulingから再現性を分離する。

## 2. 評価・簡約・式木

| 項目                            |          現在値 | 動作                                                  |
| ------------------------------- | --------------: | ----------------------------------------------------- |
| Evaluator depth                 |            1024 | 通常評価の病的nesting/recurison安全弁                 |
| CertifiedEvaluator AST depth    |              96 | OS stack overflow前にcertified evaluation対象外へする |
| Simplifier pass                 |              32 | canonical simplificationの反復上限                    |
| FullSimplify candidates         |              96 | exact同値候補探索の上限                               |
| `expand` expanded terms         |            4096 | 超過時は巨大展開を作らず元部分式を保持                |
| Polynomial conversion degree    |            4096 | 多項式化による式爆発防止                              |
| Polynomial conversion terms     |            4096 | 同上                                                  |
| 一般trig monomial reduction既定 | total degree 64 | 共通helperの既定budget                                |

## 3. 微分・積分・極限

| 項目                                  |        現在値 | 動作                                                   |
| ------------------------------------- | ------------: | ------------------------------------------------------ |
| `D[...,{x,n}]` order                  |      最大4096 | それ以上はOverflow error                               |
| symbolic integrate recursion depth    |            24 | それ以上は探索打切り                                   |
| integration substitution candidates   |            32 | 候補爆発防止                                           |
| trig integer-power integration        |   主に最大256 | `sin/cos/tan/cot/sec/csc`等の有限展開/reduction policy |
| specialized rational denominator path |  degree最大12 | 高次数の無制限partial-fraction探索を避ける             |
| direct Lambert/polylog/quadratic-exp repeated D | `n<=64` | compact exact fast path。超過時は一般`D`反復へ代替経路し，全体の`D` order上限4096は不変 |
| `limit` recursion depth               |            24 | 未解決式へ戻す                                         |
| l'Hopital iterations                  |            12 | それ以上は無制限反復しない                             |
| `nintegrate` maximum subintervals     | 2^18 = 262144 | adaptive refinement安全弁                              |
| Newton–Cotes degree                   |             8 | certified numerical integrationの固定rule              |
| numerical calculus decimal digits     |    既定100000 | 共通`maxRequestedPrecisionDigits`。独立した局所hard capは持たない |

有理函数積分は，一次・二次因子のcompact partial fractionを優先した後，Q[x]上のYun square-free decompositionとHermite reductionで重複高次因子を処理し，次数3以上のsquare-free部分をcertified Complex Rootの留数分解`P(r)/Q'(r) Log[x-r]`へ落とす。specialized rational denominator pathの現work budgetはdegree 12であり，これは数学的定義域境界ではなく式サイズ・全根isolation・algebraic residue materializationを無制限化しないためのresource-safety境界である。Rothstein–Trager / Lazard–Rioboo–Trager型の留数groupingは未実装で，現出力よりcompactな実`log/atan`表現へまとめる余地がある。

## 4. Solve / polynomial

| 項目                                  |             現在値 | 動作                                                                                                                                              |
| ------------------------------------- | -----------------: | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| symbolic-coefficient linear system    |          最大4変数 | Cramer式爆発防止。Rational coefficient系の一般Gauss-Jordanとは別                                                                                  |
| binomial polynomial solve             |      degree最大256 | `a*x^n+b`系のexact root展開budget                                                                                                                 |
| rational-function polynomial exponent |            -64..64 | solver内部変換の式爆発防止                                                                                                                        |
| polynomial conversion                 | degree/terms各4096 | 共通budget                                                                                                                                        |
| algebraic `root` defining polynomial  |       degree最大96 | Real Sturm分離・Complex certified isolationを無制限化しない。2026-08-28にoutward interval Rouché証明，Durand–Kerner早期終了，高次Root materializationの不要なfield構築除去を行った後，`solve[x^65+x+1==0,x]`約1.5秒，80次約2.9秒，96次約4.1秒まで連続的に測定できたため，旧64上限を96へ再設定した。96超は未測定領域として共通budgetで停止する |
| algebraic-field candidate             |       degree最大16 | resultant / minimal-polynomial factor reduction / primitive-element reductionの次数爆発を抑える。証明不能・超過時は代替経路またはsymbolic式を保持 |

## 5. Linear algebra

| 項目                                   |             現在値 | 動作                                                                                                              |
| -------------------------------------- | -----------------: | ----------------------------------------------------------------------------------------------------------------- |
| general exact real QR                  |      固定次数上限なし | fraction-free直交化 + full-rank Gram/symmetric Bareiss。高bitでは最終radical materializationが実用上の壁になり得る |
| upper-triangular exact QR              |   明示次数上限なし | `{I,A}` fast path                                                                                                 |
| LU                                     |  square matrixのみ | 現仕様                                                                                                            |
| determinant / inverse / eigen          |  square matrixのみ | 数学的shape要件                                                                                                   |
| symbolic determinant expansion budget  |                512 | 超える場合は無理に式を作らない                                                                                    |
| symbolic inverse expansion budget      |                256 | adjugate/minor式爆発防止                                                                                          |
| certified Matrix/LU/QR precision retry |           最大12回 | guard precisionを増やしても証明不能なら未解決                                                                     |
| SVD precision retry                    |           最大10回 | 同上                                                                                                              |
| Eigen precision retry                  |           最大10回 | 同上                                                                                                              |
| real SVD Jacobi sweeps                 |  `max(32, 8*n+16)` | 収束反復budget                                                                                                    |
| complex SVD Jacobi sweeps              | `max(32, 10*n+20)` | 同上                                                                                                              |
| Eigen QR sweeps without deflation      |   `max(256, 96*n)` | 収束しないcaseで無限反復しない                                                                                    |
| dense Matrix dimension hard max        |               なし | `size_t`/allocation overflow検査のみ。1024級は現状practical stress regime                                         |
| `ArrayExpr` packed page容量            |           1024要素 | v1.5.3内部実装policy。Arrayのdimension/element hard maxではなく，builder promotion costを最大1 pageへ制限する単位 |

## 6. `N` / numerical display / certified 計算基盤

| 項目                              |         現在値 | 動作                                              |
| --------------------------------- | -------------: | ------------------------------------------------- |
| user-requested certified decimal digits | 既定100000 | 共通`EvaluationLimits::maxRequestedPrecisionDigits`。設定可能な資源policyであり数学的固定上限ではない |
| default N digits                  |             16 | 未指定時                                          |
| default guard digits              |              8 | certified working precision                       |
| top-level `N` local refinements   |             16 | branch/極/work-boundary ambiguityでglobal budgetまで無制限retryしない |
| `:fix` / `--fix`                  |        0..1000 | 表示のみ。内部precisionではない                   |
| CertifiedEvaluator AST depth      |             96 | certified numerical pathのnesting safety          |

## 7. Certified special-function 計算基盤の代表的制限

これらは特殊函数そのものの数学的定義域ではなく，現計算基盤が保証付き評価を有限時間・有限資源で行うための実装範囲である。境界外で数学的に値が存在する場合は`DomainError`にせず，原則として`CertifiedBackendUnsupported`から`N::unsupported`へ戻して未評価式を保持する。`CertifiedBackendUnsupported`は`std::domain_error`とは独立した内部例外である。区間幅だけが境界を跨いでいる場合は，再精密化で解決し得るため`PrecisionInsufficient`を使うが，top-level `N`は局所16回で停止し，既存入力enclosureの幅が原因で解決しない場合は`N::precision`へ戻す。

| 計算基盤 | 現在の代表制限 |
| --- | --- |
| `1F1` real / complex series | 固定`\|z\|`境界なし。将来項比majorantが収束域へ入るまでを含めseries最大250000 terms，共通`EvaluationBudget`で停止 |
| `2F1` real Gauss series | `\|z\| < 1`，series最大250000 terms。固定の内部閾値は設けず，真の収束境界とresource budgetで停止 |
| `2F1` complex | `\|z\|<1`はGauss series，`\|z\|>1`は退化パラメータ・分岐切断を安全に除外できる場合だけ主値 `1/z` 解析接続。`z=1`かつ`Re(c-a-b)>0`はGauss summationで対応。他の`\|z\|=1`点はパラメータ依存の境界公式が必要なため現計算基盤では未対応 |
| elliptic `F/E` real | 保証付き実効tail ratio `<=9/10`なら最大200000-term Legendre seriesをfast pathとして使用。それ以外はCarlson `RF/RD` duplication＋分岐証明へ接続。旧`\|m\|<=9/10` capability boundaryは撤去済み。一般complex 解析接続は未対応 |
| elliptic `Pi` real | 保証付き実効tail ratio `<=9/10`なら最大200000-term Legendre seriesをfast pathとして使用。それ以外はCarlson `RF/RJ` duplication＋極/分岐証明へ接続。旧`\|m\|,\|n\|<=9/10` capability boundaryは撤去済み。極を跨ぐ主値 計算基盤と一般complex 解析接続は未対応 |
| `Ei` real | 固定magnitude境界なし。moderate argumentは最大200000-term Taylor，正/負大引数は保証付き漸近展開，共通`EvaluationBudget`で停止 |
| `Si/Ci` real | 固定magnitude境界なし。moderate argumentは最大200000-term Taylor，大きな正実数は`f/g`保証付き漸近展開，共通`EvaluationBudget`で停止 |
| `Ei` complex | 固定`\|z\|`境界なし。moderate argumentはguard付きcomplex series，大引数はDLMF 6.12のcertified `E1`漸近展開と主値 Log connectionを使用。負実軸cutを跨ぐfinite-precision enclosureは`N::precision`，剰余証明が閉じなければseriesへ代替経路 |
| `Ci` complex | 固定`\|z\|`境界なし。moderate argumentはguard付きcomplex series，大引数は`Ci(z)=-1/2(E1(iz)+E1(-iz))`と左半平面reflection，純虚軸は`Chi`退化を使用。cut/極ではなく分岐側の証明と共通`EvaluationBudget`で停止 |
| Fresnel C/S complex | 固定`\|z\|`境界なし。moderate/対角方向は最大200000-term Maclaurin，軸近傍大引数は90度回転で`\|tan(arg z)\|<=1/4`の保証wedgeへ写して`f/g`漸近展開を使用。最初の未使用項による剰余証明が閉じなければseriesへ代替経路し，共通`EvaluationBudget`で停止 |
| positive-order `polylog` | 固定`49/50`境界なし。`\|z\|<1`は最大1000000-term interval series，`Li_2`は分岐を証明できるDLMF connection formula，正実`z≈1`のorder 3～12は`mu=log(z)`整数極限展開をfast pathとして使用。exact正実分岐切断上の`Li_2`は未対応 |
| EulerGamma 計算基盤 | internal `n < 2^20` safety bound |
| `zeta` certified real | `s>1`は単調real fast path。それ以外も`s=1`を除きcomplex Euler-Maclaurin / functional equationへ接続。planner/budgetで証明できない場合だけ計算量制限付き扱い |
| `zeta` certified complex | `Re(s)>=0`はEuler-Maclaurin `N<=128`, correction `k<=48`，左半平面はfunctional equation。`s=1`のみ極 |
| `digamma/trigamma` certified real | 現在`x>0`，recurrence shift `<=1000000`，Bernoulli asymptotic `k<=64` |
| exact integer `trigamma[n]` | `n<=100000`で有限二次調和和へexact還元 |
| `ibeta[a,b,x]` exact finite sum | positive integer `a,b<=4096`，`0<=x<=1` |
| `ibeta` certified real | finite-precisionを含むreal `a,b>0`, `x in [0,1]`。パラメータに関する単調性で区間端点をcertifyし，2F1/Beta 計算基盤のbudgetを継承 |
| complex `LambertW` | 任意整数branch番号を受理。主値級数/縮小写像，分岐付きLog縮小写像に加え，`-1/e`近傍は平方根局所座標で`k=0/-1`（および対称な`k=1`下側）をcertifyする。分岐点から離れた負実軸上など，現保証boxを閉じられない残存領域は`N::unsupported` |

旧`49/50`や複素`Ei/Ci`の`512/128`等は数学的な定義域・収束半径そのものではなく，性能・保証計算量のpolicyであった。複素`Ei/Ci`はcertified `E1`漸近計算基盤と主値 connectionを導入し，旧`512/128`固定閾値を撤去した。`polylog`もtermを巨大exact Rationalで保持しないinterval recurrenceへ変更し，`Li_2`のconnection formulaと高位整数orderのnear-one展開を追加した後に旧`49/50`固定閾値を撤去した。`2F1`はcertification stateの`RealInterval`化後に再測定し，旧`9/10`固定閾値を撤去した。楕円積分もseries状態のinterval化後にCarlson symmetric forms `RF/RD/RJ`を第二計算基盤として導入し，旧`9/10`固定閾値を撤去した。現在はパラメータ値ではなく実効tail ratioでseries fast pathを選び，それ以外は実積分路上のbranch/極判定とCarlson duplicationの保証収束で閉じる。`2F1`はGauss級数の真の収束域`|z|<1`をそのまま試行し，unit circleへの接近による自然な収束悪化はseries term capと共通resource budgetで制御する。`z=1, Re(c-a-b)>0`はGauss summationで閉じる。他の`|z|=1`点はパラメータ依存の境界公式を別途要する。

## 8. Real equation proof layerの現在境界

Real equalityのnonexistence / uniqueness proofは，数値samplingではなくexact certificateだけを使う。`real_function_analysis`は有限個のalgebraic 定義域境界をconnected intervalsへ分解し，exact 臨界点で再分割してpieceごとのstrict 単調性・one-sided endpoint limit・値域を証明するため，半直線や有限個の極で分断された定義域も一般proofへ供給できる。

現在の境界は，定義域 piece最大24，臨界点最大16，解析式最大256 nodeの計算量制限付きである。`tan`等の無限periodic 極集合，一般特殊函数のbranch/極集合，および複数のcomplex-valued subexpressionが相殺して実数へ戻り得る式ではreal-valued 定義域の完全性を捏造せずUnknownへ戻す。臨界点をexactに列挙できない場合も同様である。さらに，一意性だけを証明できてもexact root表現が無ければ`UnresolvedSolutionSet`を維持する。

また「一意に存在する」と「mmCalのexact式としてそのrootを構成できる」は別である。`erf[x]==1/2`や`cos[x]==x`のように一意性を証明できても，対応するinverse builtinまたは一般transcendental Root表現がない場合は`UnresolvedSolutionSet`を維持する。`f'>=0` / `f'<=0`からのstrictnessは，`f'=0`の完全解集合が有限またはIntegerパラメータ族として高々可算で，実intervalを含まないことをexactに証明できる場合まで実装した。今後の候補は，periodic 極を含む定義域表現，より一般のendpoint/asymptotic sign，zero setが離散であることの追加certificate，および必要ならalgebraic `Root`とは区別したtranscendental root certificate表現である。

Real exponential classifierは`exp[p x+q]==c x+d`とpositive constant baseの`a^(m x+n)==c x+d`をLambert Wへexactに変換し，変換argumentと`-1/E`のcertified orderingで実分岐数を完全分類できる場合まで扱う。主値 `x^x==r`は負実軸で一般に複素値となり，`0<r<1`では負の偶整数解が混在し得るため，現段階では`r>1`とexactな`r=1,0,-1`，`r<-1`だけをcompleteとして扱う。`x^x==1/4`等は正実分岐だけを返さずUnresolvedを維持する。

## 9. 軽量数論計算基盤

| 項目 | 現在値 | 動作 |
| --- | ---: | --- |
| `isprime`証明範囲                | `0..2^64-1` | first 12 prime basesのdeterministic strong Miller-Rabin                  |
| `nextprime/prevprime`            |  `uint64`内 | 候補をdeterministic primalityで検証                                      |
| `factorint/totient`              | `\|n\|<=2^64-1` | Pollard-Rho分割 + deterministic primality verification |
| Pollard-Rho polynomial パラメータ |  `c=1..127` | 1 パラメータにつき最大2,000,000 iteration。分割できなければ誤答せず未評価 |

`uint64`を超えるBigInt自体はmmCalで表現できるが，現在の数論proof 計算基盤の対象外である。probable-prime判定をexact `True`へ昇格しない。

## 10. Algorithm selection 閾値 — 上限ではない

以下は性能測定で選んだdispatch境界であり，数学的・意味論的な制限ではない。

| 項目                           | 現在値 | 意味 |
| ------------------------------ | -----: | ---- |
| BigUInt Karatsuba              | 約48 limbs | multiplication dispatch |
| BigUInt top-level Toom-3       | 約1280 limbs | top-level multiplication dispatch |
| Toom-3 recursive               | 約448 limbs | recursive multiplication dispatch |
| square Karatsuba               | 約48 limbs | squaring dispatch |
| Burnikel–Ziegler division      | 約32 limbs + offset policy | division dispatch |
| decimal divide-and-conquer     | 約128 limbs | decimal conversion dispatch |
| certified non-power-of-two FFT | 384点未満direct DFT，それ以上Bluestein | `--fft-threshold`で環境別再測定可能 |
| exact Cyclotomic FFT           | `phi(n) <= 64` | 対応exact入力を`Q[t]/Phi_n(t)`へ写せる場合の現在計算基盤 budget。超過時はgeneric exact DFTへ代替経路 |
| QR column block default        | 1（unblocked相当） | 現在の実測dispatch |

## 11. 「上限がない」もの

- BigUInt / BigIntの桁数に32/64/128-bitの固定数学上限はない。
- Rationalの分子・分母も任意精度である。
- BigFloatのprecisionにIEEE754型のような固定53-bit上限はない。
- brace / Arrayのshapeに1024等の固定次数上限はない。
- `@`, `@@`, ... / `%`, `%%`, ... の連続個数に固定上限はない。
- `N`の一般precisionに1000桁という制限はない。1000は`:fix`の表示上限である。要求precisionには既定100000桁の共通`EvaluationBudget`上限があるが，これは設定可能な資源policyであり数値表現の固定上限ではない。

実際の上限はmemory，時間，`size_t`，共通`EvaluationBudget`，および各symbolic/certified algorithmが設ける個別budgetで決まる。

## 12. Benchmark multi-thread方針

Benchmark executableだけをmulti-thread化することは可能である。ただし用途を分離する。

### 推奨

- `--random-expressions --threads N`: 推奨。各workerが独立`KernelSession`を所有する。
- fixed-seed random invariant: caseをworkerへ分配可能。
- throughput stress benchmark: 独立operationをN workerで実行してops/sを測る。

### 既存timingでは非推奨

現在の`benchmarkMultiply`, `benchmarkExactFft`, `N[SVD]`等はsingle-operation latencyを比較するための値である。これを複数threadで同時実行するとallocator contention，CPU boost，cache，memory bandwidthを含むthroughput値になり，過去のsingle-thread benchmarkと比較できなくなる。

### 1個のMatrix演算を速くしたい場合

Benchmark側だけでN threadsを作っても，一つの`svd[A]`や`eigenvalues[A]`の内部計算は速くならない。Householder/Jacobi/QR/Matrix multiply等のCore kernel自体をparallelizeする必要があり，これはBenchmark限定変更ではなくCoreの算法変更である。

### Thread safety

同じ`KernelSession` / `Evaluator` / `FourierTransformCache`を複数threadから共有してはならない。workerごとに完全に分離する。現在の設計はsession-local stateが多いためこの方式とは相性が良いが，正式にthread-safeを保証する前にはThreadSanitizer等でrace監査する。

多変数polynomial `solve`のnon-shape 代替経路は，solver変数最大4，eliminant root最大64の特殊化 budgetで低次元systemを再帰処理する。positive-dimensional系はexact factor component，定数係数の一次変数消去，または1方程式・次数1/2 射影として完全性を構造的に証明できる場合だけパラメータ化する。一次消去は最大6変数までとし，非定数係数で追加case splitが必要な場合は現段階では使わない。Realの次数1/2 射影は自由パラメータのReal 定義域と，分母非零・二次radicand非負の半代数条件をexactに保持し，1自由パラメータなら一変数polynomial inequalityへ戻してinterval化する。複数パラメータでは証明済みPredicateをそのまま保持する。function field上の高次algebraic equation，一般多様体のrational parameterization，完全な半代数分解を必要とする一般多方程式正次元系は現段階ではUnresolvedを維持する。
