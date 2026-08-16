# mmCal v1.5.2 意図的な制限・探索budget・実装threshold一覧

> 対象: v1.5.2 + Random Expression Fuzzer追加後のソース
>
> 本文では「数学上の制限」「安全弁」「探索budget」「実装threshold」を区別する。thresholdは入力を拒否する上限ではない。

## 1. Random Expression Fuzzer

| 項目 | 現在値 | 意味 |
|---|---:|---|
| `maxDepth`既定値 | 16 | hard limitではない。`--max-depth N`で変更可能 |
| `--max-depth` CLI hard max | なし | 1以上だけ検査。上げすぎるとAST/式サイズが急増する |
| target depth 1–4 | 60% | 通常case |
| target depth 5–7 | 25% | 中程度 |
| target depth 8–10 | 10% | 深め |
| target depth 11–13 | 4% | stress |
| target depth 14–max | 1% | rare stress |
| `expand`用generator depth | 最大6 | 多項式爆発を避けるfuzzer側制限 |
| `factor`用generator depth | 最大5 | 同上 |
| scalar exponent | 1..5 | generator側 |
| polynomial exponent | 1..4 | generator側 |
| Matrix shape | 1..4 × 1..4 | generator側 |
| Matrix integer entry | -5..5 | generator側 |
| report interval | 10000 | `--report-every`で変更可能 |
| finite cases | 10000 | `--cases`で変更可能 |
| `--loop` | case数上限なし | 最初のFAILで停止 |

### loopのsession状態

Unreleasedでは，各workerが独立`KernelSession`を所有し，case間・比較評価間に`resetForIndependentEvaluation()`を呼ぶ。これは`history_` / `inputHistory_` / `outputHistory_` / diagnostic / 定義 / `inputCount_`を初期化する一方，角度設定とRNG streamを保持する。通常の`reset()`は同じtransient resetの後にentropy reseedするため，fuzzerだけが不必要なreseedingを避けられる。

したがって，旧実装で問題だった「`clearHistory()`しても`inputCount_`が増え続け，次回`evaluate()`が巨大vectorへresizeする」経路は解消済みである。参考測定では同一Release/LTO-off build，`--threads 8 --seed 1234`で100,000 caseの最大RSSが`48136 KiB`，200,000 caseが`51364 KiB`であり，case数に比例する履歴蓄積は見られなかった。overnight burn-inではallocator high-water等を含め，process RSSの定常性を引き続き監視する。

現generatorはrandom builtinを生成しない。将来これを追加する場合は，case番号からrandom builtin用seedも明示的に導出し，thread schedulingから再現性を分離する。

## 2. 評価・簡約・式木

| 項目 | 現在値 | 動作 |
|---|---:|---|
| Evaluator depth | 1024 | 通常評価の病的nesting/recurison安全弁 |
| CertifiedEvaluator AST depth | 96 | OS stack overflow前にcertified evaluation対象外へする |
| Simplifier pass | 32 | canonical simplificationの反復上限 |
| FullSimplify candidates | 96 | exact同値候補探索の上限 |
| `expand` expanded terms | 4096 | 超過時は巨大展開を作らず元部分式を保持 |
| Polynomial conversion degree | 4096 | 多項式化による式爆発防止 |
| Polynomial conversion terms | 4096 | 同上 |
| 一般trig monomial reduction既定 | total degree 64 | 共通helperの既定budget |

## 3. 微分・積分・極限

| 項目 | 現在値 | 動作 |
|---|---:|---|
| `D[...,{x,n}]` order | 最大4096 | それ以上はOverflow error |
| symbolic integrate recursion depth | 24 | それ以上は探索打切り |
| integration substitution candidates | 32 | 候補爆発防止 |
| trig integer-power integration | 主に最大256 | `sin/cos/tan/cot/sec/csc`等の有限展開/reduction policy |
| specialized rational denominator path | degree最大12 | 高次数の無制限partial-fraction探索を避ける |
| `limit` recursion depth | 24 | 未解決式へ戻す |
| l'Hopital steps | 12 | それ以上は無制限反復しない |
| `nintegrate` maximum subintervals | 2^18 = 262144 | adaptive refinement安全弁 |
| Newton–Cotes degree | 8 | certified numerical integrationの固定rule |
| numerical calculus decimal digits | 最大100000 | tolerance構築の明示上限 |

## 4. Solve / polynomial

| 項目 | 現在値 | 動作 |
|---|---:|---|
| symbolic-coefficient linear system | 最大4変数 | Cramer式爆発防止。Rational coefficient系の一般Gauss-Jordanとは別 |
| binomial polynomial solve | degree最大256 | `a*x^n+b`系のexact root展開budget |
| rational-function polynomial exponent | -64..64 | solver内部変換の式爆発防止 |
| polynomial conversion | degree/terms各4096 | 共通budget |
| algebraic `root` defining polynomial | degree最大64 | Real Sturm分離・Complex certified isolation/refinementを無制限化しない |
| algebraic-field candidate | degree最大16 | resultant / minimal-polynomial factor reduction / primitive-element reductionの次数爆発を抑える。証明不能・超過時はfallbackまたはsymbolic式を保持 |

## 5. Linear algebra

| 項目 | 現在値 | 動作 |
|---|---:|---|
| general exact Householder QR | `min(m,n) <= 3` | 4次以上のradical式爆発を防ぐ |
| upper-triangular exact QR | 明示次数上限なし | `{I,A}` fast path |
| LU | square matrixのみ | 現仕様 |
| determinant / inverse / eigen | square matrixのみ | 数学的shape要件 |
| symbolic determinant expansion budget | 512 | 超える場合は無理に式を作らない |
| symbolic inverse expansion budget | 256 | adjugate/minor式爆発防止 |
| certified Matrix/LU/QR precision retry | 最大12回 | guard precisionを増やしても証明不能なら未解決 |
| SVD precision retry | 最大10回 | 同上 |
| Eigen precision retry | 最大10回 | 同上 |
| real SVD Jacobi sweeps | `max(32, 8*n+16)` | 収束反復budget |
| complex SVD Jacobi sweeps | `max(32, 10*n+20)` | 同上 |
| Eigen QR sweeps without deflation | `max(256, 96*n)` | 収束しないcaseで無限反復しない |
| dense Matrix dimension hard max | なし | `size_t`/allocation overflow検査のみ。1024級は現状practical stress regime |
| `ArrayExpr` packed page容量 | 1024要素 | Unreleased内部実装policy。Arrayのdimension/element hard maxではなく，builder promotion costを最大1 pageへ制限する単位 |

## 6. `N` / numerical display / certified backend

| 項目 | 現在値 | 動作 |
|---|---:|---|
| `N[expr,p]`の一般decimal hard max | 明示固定値なし | memory, `size_t`, 個別backendのbudgetで制約される |
| default N digits | 16 | 未指定時 |
| default guard digits | 8 | certified working precision |
| `:fix` / `--fix` | 0..1000 | 表示のみ。内部precisionではない |
| CertifiedEvaluator AST depth | 96 | certified numerical pathのnesting safety |

## 7. Certified special-function backendの代表的制限

これらは特殊函数そのものの数学的domainではなく，現backendが保証付き級数で安全に処理するための実装範囲である。

| backend | 現在の代表制限 |
|---|---|
| `1F1` | series最大200000 terms |
| `2F1` | series最大250000 terms |
| elliptic series | 最大200000 terms |
| `Ei` | 現series pathは `|x| <= 8`，最大200000 terms |
| `Si` | 現series pathは `|x| <= 8`，最大200000 terms |
| `Ci` | real certified pathは `0 < x <= 8`，最大200000 terms |
| positive-order `polylog` | 現series pathは `|z| < 1`，最大1000000 terms |
| EulerGamma backend | internal `n < 2^20` safety bound |
| `zeta` certified real | 現在`x>1`，Euler-Maclaurin tail start `N<=4096`，Bernoulli correction `k<=64` |
| `digamma/trigamma` certified real | 現在`x>0`，recurrence shift `<=1000000`，Bernoulli asymptotic `k<=64` |
| exact integer `trigamma[n]` | `n<=100000`で有限二次調和和へexact還元 |
| `ibeta[a,b,x]` exact finite sum | positive integer `a,b<=4096`，`0<=x<=1` |
| `ibeta` certified real | exact Rational `a,b>0` + certified real `x in [0,1]`。2F1/Beta backendのbudgetを継承 |

これらを超えた入力に対して，mmCalは誤った近似を返すのではなく`PrecisionInsufficient`等で「現在のbackendでは保証できない」とする。

## 8. 軽量数論backend

| 項目 | 現在値 | 動作 |
|---|---:|---|
| `isprime`証明範囲 | `0..2^64-1` | first 12 prime basesのdeterministic strong Miller-Rabin |
| `nextprime/prevprime` | `uint64`内 | 候補をdeterministic primalityで検証 |
| `factorint/totient` | `|n|<=2^64-1` | Pollard-Rho分割 + deterministic primality verification |
| Pollard-Rho polynomial parameter | `c=1..127` | 1 parameterにつき最大2,000,000 iteration。分割できなければ誤答せず未評価 |

`uint64`を超えるBigInt自体はmmCalで表現できるが，現在の数論proof backendの対象外である。probable-prime判定をexact `True`へ昇格しない。

## 9. Algorithm selection threshold — 上限ではない

以下は性能測定で選んだdispatch境界であり，数学的・意味論的な制限ではない。

| 項目 | 現在値 |
|---|---:|
| BigUInt Karatsuba | 約48 limbs |
| BigUInt top-level Toom-3 | 約1280 limbs |
| Toom-3 recursive | 約448 limbs |
| square Karatsuba | 約48 limbs |
| Burnikel–Ziegler division | 約32 limbs + offset policy |
| decimal divide-and-conquer | 約128 limbs |
| certified non-power-of-two FFT | 96点未満direct DFT，それ以上Bluestein |
| QR column block default | 1（unblocked相当） |

## 10. 「上限がない」もの

- BigUInt / BigIntの桁数に32/64/128-bitの固定数学上限はない。
- Rationalの分子・分母も任意精度である。
- BigFloatのprecisionにIEEE754型のような固定53-bit上限はない。
- brace / Arrayのshapeに1024等の固定次数上限はない。
- `@`, `@@`, ... / `%`, `%%`, ... の連続個数に固定上限はない。
- `N`の一般precisionに1000桁という制限はない。1000は`:fix`の表示上限である。

実際の上限はmemory，時間，`size_t`，および各symbolic/certified algorithmが設ける個別budgetで決まる。

## 11. Benchmark multi-thread方針

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
