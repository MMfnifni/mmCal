# mmCal.Benchmarks

`mmCal.Benchmarks`は，mmCal本体の通常回帰testとは分離した，**性能測定・固定seedランダム不変量試験・大規模Matrix負荷試験・grammar-aware数式fuzzing**のための開発用Console programである。

通常の`mmCal.Tests`や`test_set`が「既知の仕様を壊していないか」を確認するのに対し，`mmCal.Benchmarks`は主に次を担当する。

- BigInt / BigFloat / certified算法の性能変化を測る
- Karatsuba，Toom-3，Burnikel–Ziegler，FFT等のthresholdを実測する
- fixed-seed random inputで数学的不変量を検査する
- Matrix，FFT，特殊函数等の数値backendがcertified relationを維持しているか確認する
- 大きなMatrixで時間・メモリ特性を調べる
- 合法な数式を大量生成し，函数同士の未知の組合せ不具合を探索する

benchmark値はCPU，compiler，最適化設定，allocator，OS，cache状態等に依存する。**異なる環境の数値を直接性能保証として比較しないこと。** thresholdを変更する場合は，実際に採用するcompilerとbuild configurationで再測定する。

---

## 1. Build

### Visual Studio / MSVC

`mmCal.sln`には`mmCal.Benchmarks` projectが含まれている。

性能測定では原則として，

```text
Release | x64
```

を使用する。

Debug buildはassertion，最適化不足，iterator/debug runtime等の影響を大きく受けるため，性能比較には使用しない。

### CMake

CMakeでは`mmCal.Benchmarks`は`EXCLUDE_FROM_ALL` targetであるため，明示的にbuildする。

```text
cmake --build build --config Release --target mmCal.Benchmarks
```

生成物の場所はgenerator / build環境に依存する。

---

## 2. まず使うcommand

```text
mmCal.Benchmarks
```

fixed-seed random invariant checkを実行した後，通常規模のbenchmark一式を実行する。

より重い測定：

```text
mmCal.Benchmarks --full
```

random invariantだけ：

```text
mmCal.Benchmarks --random-only
```

性能測定だけ：

```text
mmCal.Benchmarks --benchmark-only
```

大きなMatrixを一種類だけ測定：

```text
mmCal.Benchmarks --matrix-large nsvd 64 16
```

random expression fuzzerを有限回実行：

```text
mmCal.Benchmarks --random-expressions --cases 50000 --seed 1234
```

FAILが出るまで無期限に回す：

```text
mmCal.Benchmarks --random-expressions --loop --threads 8
```

FAILを報告しつつ無期限に回し続ける：

```text
mmCal.Benchmarks --random-expressions --nostop-loop --threads 8
```

特定caseを直接再現：

```text
mmCal.Benchmarks --random-expressions --seed 1234 --case 48172
```

help：

```text
mmCal.Benchmarks --help
```

---

# 3. 通常mode

## 3.1 Default

```text
mmCal.Benchmarks
```

次のrandom invariantを先に実行する。

| 分類 | 既定case数 |
|---|---:|
| Random BigInt | 300 |
| certified exp / log | 40 |
| 2F1 / elliptic | 40 |
| certified Matrix | 40 |
| certified FFT | 40 |

すべてPASSした後，通常benchmarkを実行する。

## 3.2 `--full`

```text
mmCal.Benchmarks --full
```

random試験数を増やし，通常benchmarkより大きなoperand，高精度，高次Matrix / FFTまで測る。

| 分類 | `--full` case数 |
|---|---:|
| Random BigInt | 2000 |
| certified exp / log | 200 |
| 2F1 / elliptic | 200 |
| certified Matrix | 200 |
| certified FFT | 200 |

日常的なcommit前確認より，算法threshold変更，高速化patch，release前の性能監査向けである。

## 3.3 `--random-only`

```text
mmCal.Benchmarks --random-only
```

fixed-seed invariantだけを実行し，timing benchmarkは行わない。

算法実装を変更した直後に，まず「速さ以前に数学的に壊れていないか」を確認する用途を想定している。

## 3.4 `--benchmark-only`

```text
mmCal.Benchmarks --benchmark-only
```

random invariantを省略し，timingだけを測る。

性能を何度も測るときに使用する。正当性確認を置き換えるoptionではない。

`--random-only`と`--benchmark-only`は同時指定できない。

---

# 4. 通常benchmarkの内容

## 4.1 BigInt multiply / square / division

32-bit limb数を変えながら，

- multiplication
- dedicated square
- division

を測る。

通常mode：

```text
64, 128, 256, 512, 1024 limbs
```

`--full`：

```text
64, 128, 256, 512, 1024, 2048, 4096 limbs
```

出力単位は1回あたり`us`である。

Karatsuba / Toom-3 / Burnikel–Ziegler等のcrossoverは数学定数ではなく実装環境依存である。threshold変更時はこの領域を特に比較する。

## 4.2 Factorial

product tree等を含む巨大factorial本体の時間を測る。

通常：

```text
10000!, 20000!, 40000!
```

`--full`では`80000!`まで測る。

## 4.3 Decimal conversion

factorialで作った巨大BigIntについて，

```text
BigInt -> decimal text -> BigInt
```

の`toString`と`parse`を別々に測る。

巨大整数では計算本体より10進変換が支配的になることがあるため，factorial timingとは分離している。

## 4.4 Exact / certified Matrix

主に以下を測る。

- exact `dot`
- exact `det`
- exact `rref`
- exact `solveLinear`
- exact `nullSpace`
- exact LU
- certified `N[det,16]`
- certified `N[LU,16]`
- certified `N[QR,16]`

exact整数 / Rational pathとprecision-aware certified pathを同じ表の中で比較できる。

## 4.5 Exact Householder QR

一般exact QRはradical式が急激に膨張するため，通常benchmarkでは2×2，3×3だけを測る。

これは性能制約を意図的に監視するtestでもある。exact-firstは「巨大exact式を無条件に展開する」という意味ではない。

## 4.6 Householder QR block sweep

certified QRについて，column block幅

```text
1, 8, 16, 32
```

を比較する。

現時点ではBigFloat / interval演算コストが支配的で，一つのblock幅が全サイズで安定して勝つとは確認されていない。そのため自動block化の採否はbenchmarkで再評価する前提である。

## 4.7 SVD

certified reduced SVDを，主に4，8，12，16次で測る。

SVDはreconstruction / orthogonality relationを監査する数値backendであり，単純な`double` LAPACK benchmarkとは性質が異なる。

## 4.8 Eigen / Eigensystem

certified approximate eigenvalue / eigensystem backendを4，8，12，16次で測る。

一般非正規行列では固有値問題そのものが摂動に敏感であるため，時間だけでなくrandom invariant側のrelation checkと併せて評価する。

## 4.9 FFT

2冪長についてexact FFTとcertified approximate FFTを比較する。

通常：

```text
32, 64, 128 points
```

`--full`では512点まで拡大する。

また非2冪長について，direct DFTとFFT/Bluesteinの実測を比較する。現在のdirect/Bluestein thresholdも，この測定から再評価する前提である。

## 4.10 Certified exp / log

`exp[1]`，`log[2]`を100～数千桁で測る。

`--full`では10000桁まで拡大し，さらに一般Rational

```text
123456789/987654321
```

に対するcertified exp / logも測る。

---

# 5. Fixed-seed random invariant

通常のrandom checkは**再現性を優先して固定seed**を使用する。

対象は概ね次の通りである。

### BigInt

巨大random整数について，算術，division，変換等のinvariantを検証する。

### certified exp / log

exact Rational値や函数恒等式がcertified intervalに包含されることを検査する。

例として，正の`x`について

```text
log[x] + log[1/x]
```

のintervalが0を包含すること等を確認する。

### 2F1 / elliptic

既知の退化式・恒等式をrandomなsafe-domain pointへ適用し，certified enclosureがexact値を包含することを確認する。

### Matrix

Matrix random checkでは，単に期待文字列と比較するのではなく，函数ごとの数学的不変量を確認する。

例：

- `A . inverse[A] == I`
- RREF
- `solveLinear`
- `A . v == 0` for `v` in `nullSpace[A]`
- `P A == L U`
- QR reconstruction / orthogonality
- SVD reconstruction / orthogonality
- Eigen relation

整数・Rational・certified approximate pathを跨いで監視する。

### FFT

FFT round-tripやcertified approximate transformが期待関係を満たすことを固定seedで検査する。

---

# 6. Large Matrix benchmark

```text
mmCal.Benchmarks --matrix-large <op> <size> [digits]
```

一種類のMatrix operationだけを指定サイズで測定する。

例：

```text
mmCal.Benchmarks --matrix-large nlu 64 16
mmCal.Benchmarks --matrix-large nqr 64 16
mmCal.Benchmarks --matrix-large nsvd 64 16
mmCal.Benchmarks --matrix-large neigen 64 16
```

`digits`省略時は16である。

## 6.1 operation一覧

| op | 内容 | `digits` |
|---|---|---|
| `transpose` | exact transpose | 不使用 |
| `trace` | exact trace | 不使用 |
| `ndot` | certified approximate matrix product | 使用 |
| `det` | exact determinant | 不使用 |
| `ndet` | certified approximate determinant | 使用 |
| `ninv` | certified approximate inverse | 使用 |
| `rref` | exact RREF | 不使用 |
| `rank` | exact matrix rank | 不使用 |
| `nrank` | certified approximate matrix rank | 使用 |
| `nsolve` | certified approximate linear solve | 使用 |
| `nnull` | certified approximate null space | 使用 |
| `lu` | exact LU decomposition | 不使用 |
| `nlu` | certified approximate LU | 使用 |
| `nqr` | certified approximate Householder QR | 使用 |
| `nsvd` | certified approximate SVD | 使用 |
| `neigen` | certified approximate eigenvalues | 使用 |
| `neigensystem` | certified approximate eigensystem | 使用 |

## 6.2 入力Matrix

large Matrix benchmarkはparser costを測るものではない。

添付random matrix generatorと同系統の，

```text
[-1,1]
小数10桁相当
```

の値をC++側でexact Rationalとして直接構築する。

したがって，

```text
--matrix-large nsvd 64 16
```

のtimingには，巨大な`{{...},{...}}`文字列をLexer / Parser / Lowererへ通す時間は含まれない。

Python generatorとRNG列そのものを一致させる目的でもない。**分布と桁幅を近づけ，算法時間を独立測定するためのfixture**である。

巨大text inputのparse / lowering性能を調べる場合は，CLIへ実際のbrace形式を入力して別途測定する。

## 6.3 1024×1024について

1024×1024 dense matrixは1,048,576要素を持つ。

`double + BLAS`では特別巨大とは限らないが，mmCalではexact Rational / Expr，あるいはBigFloat / intervalを各要素に持つため，メモリと演算コストは大きく異なる。

特にO(n^3)算法では，64次から1024次へ上げるだけで理想的な三次則でも約4096倍の仕事量になる。

したがって1024次を試す場合は，

1. `transpose` / `trace`等のO(n)～O(n²)処理
2. LU / dot等
3. QR / SVD / Eigen

の順に段階的に上げ，process memoryも同時に監視することを推奨する。

巨大Matrixでtime outした結果を「算法が壊れた」と即断しない。現在は多倍長算術，parser / lowering，allocation，working precisionも重要な性能要因である。

Unreleasedの第一段階typed-node化では，同一GCC Release/LTO-offの`--matrix-large transpose 1024 16`で旧variant版`693312 KiB`→typed-node版`299668 KiB`となり，最大RSSを約56.8%削減した。第二段階では`ArrayExpr`を固定1024要素のimmutable packed page + shape/offset/stridesへ移行し，benchmark fixtureも`ArrayBuilder`からexact Rationalを直接構築する。これにより同負荷は最大RSS約`136576 KiB`（133.4 MiB），transpose本体約0.059 msとなった。transposeはデータcopyではなくlayout view生成なので，旧実装のtranspose timingと算法量そのものが異なる点に注意する。

`ArrayBuilder`はpromotionを現在page内へ限定する。1,048,576要素の最後だけsymbolic値にした専用測定でもall-integer版とほぼ同じ約0.30 s / 約69.8 MiBで，先行pageはInteger packedのまま，最終pageだけGenericとなった。これは「最後の1要素で全ArrayをGenericへ作り直す」実装を避けるための設計である。

巨大text inputは別測定する。13.63 MBの1024×1024・10桁decimal literalをCLIへ渡して`dimensions[...]`だけを評価した現Unreleased測定は約4.55 s / 最大RSS `441564 KiB`（431 MiB）。v1.5.2文書の約14.16 MB / 10.9 s / 1.99 GBとは入力textが完全同一ではないため厳密なA/B値には使わない。

---

# 7. Random Expression Fuzzer

```text
mmCal.Benchmarks --random-expressions
```

文法と型をある程度理解した**grammar-aware semantic fuzzer**である。

完全なrandom token列を生成するのではなく，原則として評価可能なexact算術，多項式，小型Matrix等を生成し，数学的不変量を検査する。

目的は，

> 固定testでは思いつかなかった函数・構文・簡約の組合せを自動探索すること

である。

## 7.1 有限実行

既定は10000 caseである。

```text
mmCal.Benchmarks --random-expressions
```

件数指定：

```text
mmCal.Benchmarks --random-expressions --cases 50000
```

再現可能なmaster seedを固定：

```text
mmCal.Benchmarks --random-expressions --cases 50000 --seed 1234
```

## 7.2 無限loop

```text
mmCal.Benchmarks --random-expressions --loop
```

case上限を設けず継続する。

通常の`--loop`は**最初のFAILを検出した時点で必ず停止する。** 同じ根本原因から派生したfailureを大量保存するより，最初の再現可能な反例を小さくすることを優先する。

FAIL後もburn-inを継続したい場合は，

```text
mmCal.Benchmarks --random-expressions --nostop-loop
```

を使用する。こちらは各FAILをshrinking・表示した後，次caseへ進む。通常終了は`Ctrl+C`で行う。

## 7.3 複数thread

```text
mmCal.Benchmarks --random-expressions --cases 100000 --threads 8
mmCal.Benchmarks --random-expressions --loop --threads 8
```

`--threads N`は**case単位のthroughput並列化**である。1個のSVDや1個の式評価を内部でN thread化するoptionではない。各workerは独立した`KernelSession`を所有し，mutableな履歴・定義・diagnostic stateを共有しない。

各caseはmaster seed + 1-based case番号から独立生成されるため，thread schedulingや`--threads`値が変わっても，同じseed/case番号は同じ式を生成する。FAIL再現時は`--case`を使えばよく，通常1 threadで直接再現される。

worker内ではcase間・比較評価間に`KernelSession::resetForIndependentEvaluation()`を使用する。これは定義・履歴・input/output history・diagnostic・入力番号を初期化するが，通常`reset()`と異なりentropy reseedを行わない。そのため長時間loopで履歴vectorがcase数に比例して成長せず，既存のRNG streamも不用意に変更しない。現generatorはrandom builtinを生成しないため，random builtinを将来追加する場合はcase seedからの明示seed policyを別途定める。

参考として，同一Release/LTO-off build，`--threads 8 --seed 1234`で有限実行した最大RSSは100,000 caseで`48136 KiB`，200,000 caseで`51364 KiB`だった。case数を2倍にしても履歴保持に相当する線形増加は見られない。これは長時間`--loop`の完全な定常性証明ではないため，overnight burn-inでは引き続きprocess RSSを監視する。

## 7.4 FAIL時

FAILすると，概ね次の情報を表示する。

```text
[FAIL] Random expression invariant
Seed      : 92847163
Case      : 48172
Case seed : ...
Depth     : ...
Invariant : ...
Expression: ...
Reduced   : ...
Reason    : ...
Expected  : ...
Actual    : ...
Reproduce : mmCal.Benchmarks --random-expressions --seed 92847163 --case 48172
```

FAIL後は簡易shrinkingを行い，可能ならより小さい反例を`Reduced`として表示する。

通常の`--loop`またはfinite modeでFAILした場合，processはexit code 1で終了する。`--nostop-loop`ではFAILは継続中の反例として扱われ，その場ではprocessを終了しない。

## 7.5 Caseの直接再現

```text
mmCal.Benchmarks --random-expressions --seed 92847163 --case 48172
```

case番号は**1-based**である。

各caseはmaster seedとcase番号から独立したcase seedを導出して生成する。このため，48172件目を再現するために1～48171件目を再実行する必要はない。

長時間`--loop`を回す上で重要な仕様である。

## 7.6 Depth

既定最大depth：

```text
16
```

変更：

```text
mmCal.Benchmarks --random-expressions --max-depth 24
```

depthは全caseで固定ではなく，caseごとに重み付きrandomで選択される。

浅い式を主体としつつ，稀に`--max-depth`付近の深い式を生成する。すべての枝を最大depthまで伸ばして巨大ASTを作ることが目的ではない。

## 7.7 Progress表示

既定では10000 caseごとに進捗を表示する。

```text
[10000] PASS  2.34 s  max-depth=16
[20000] PASS  4.69 s  max-depth=16
```

変更：

```text
mmCal.Benchmarks --random-expressions --loop --threads 8 --report-every 100000
```

`--report-every 0`を指定した場合は定期progress表示を行わない。

## 7.8 現在のInvariant

初期実装では主に次を検証する。

### Formatter / Parser round-trip

```text
format(parse(format(expr))) == format(expr)
```

内部ASTが完全同型であることまでは要求しない。Formatterがcanonical textを安定して再生成できることを見る。

### `fullSimplify`

exact numeric expressionについて，`fullSimplify`が値を変えていないことを確認する。

### `expand`

random polynomialを数値評価し，展開前後で値が一致することを確認する。

### `factor`

random polynomialについてfactor / expand後の数学的値を確認する。

### Matrix transpose

```text
transpose[transpose[A]] == A
```

### Determinant

```text
det[A] == det[transpose[A]]
```

今後，`D` / `integrate` / `N` / `inverse` / `solveLinear` / `nullSpace`等のmetamorphic invariantを追加できる。

---

# 8. Fuzzerと通常random checkの違い

両者は目的が異なる。

### Fixed-seed random invariant

```text
mmCal.Benchmarks --random-only
```

- test内容そのものは開発者が設計する
- 入力値だけをrandom化する
- release regression向き
- 同じseedで毎回同じ集合を確認する

### Random Expression Fuzzer

```text
mmCal.Benchmarks --random-expressions --loop
```

- 数式構造そのものをrandom生成する
- 函数同士の未知の組合せを探索する
- 長時間burn-in向き
- FAILしたcaseはseed + case番号で固定testへ昇格できる

どちらか一方で代用するものではない。

---

# 9. Benchmark結果の読み方

## 9.1 1回の結果だけでthresholdを変えない

CPU boost，background process，cache，allocator等で数%程度は容易に変動する。

threshold変更では，

1. 同じRelease binary
2. 同じmachine
3. 同じ電源設定
4. 同じoperand範囲
5. 複数回測定

を基本とする。

僅差なら複雑な算法を採用しない。

mmCalでは実際に，理論上有望でもbenchmarkで優位性が得られず棄却・保留した最適化が複数存在する。

## 9.2 correctnessとspeedを混同しない

速い結果は正しい証明ではない。

高速化変更では原則として，

```text
通常unit / black-box test
↓
fixed-seed random invariant
↓
benchmark
```

の順に確認する。

## 9.3 exactとapproximateを分けて読む

例えば，

```text
exact det
N[det,16]
```

は同じ数学函数でもbackendが異なる。

exact pathはBigInt / Rational / symbolic expressionを維持し，approximate pathは`ApproximationContext`からBigFloat / intervalへ直接dispatchする場合がある。

したがって単純に「Nを付けると何倍速い」と一般化せず，どのbackendを測っているかを見る。

---

# 10. 推奨workflow

## 通常の変更

```text
mmCal.Tests
black-box tests
mmCal.Benchmarks --random-only
```

## 数値backendを変更した場合

```text
mmCal.Benchmarks --random-only
mmCal.Benchmarks --benchmark-only
```

## BigInt算法thresholdを変更した場合

```text
mmCal.Benchmarks --full --benchmark-only
```

変更前後のBigInt multiply / square / divisionを保存して比較する。

## Matrix算法を変更した場合

まず通常random invariantを通す。

```text
mmCal.Benchmarks --random-only
```

次に対象だけを段階的に測る。

```text
mmCal.Benchmarks --matrix-large nlu 32 16
mmCal.Benchmarks --matrix-large nlu 64 16
mmCal.Benchmarks --matrix-large nlu 128 16
```

いきなり1024次から始めない。

## 放置burn-in

```text
mmCal.Benchmarks --random-expressions --loop --threads 8
```

FAILが出たら，表示されたcommandをそのまま使用する。

```text
mmCal.Benchmarks --random-expressions --seed <seed> --case <case>
```

再現後，原因を修正し，可能ならその反例を通常unit / black-box regressionへ昇格する。

---

# 11. Exit code

| code | 意味 |
|---:|---|
| `0` | 指定した試験・benchmarkが正常終了 |
| `1` | random invariant / random-expression fuzzer等でFAIL |
| `2` | option不足，未知option，矛盾したoption等のCLI error |

benchmark中の内部invariant違反等ではprogramがabortする場合がある。これはtiming値を誤って正常結果として採用しないためである。

---

# 12. 現在の設計上の注意

`mmCal.Benchmarks`は一般的なbenchmark frameworkではなく，**mmCal固有の算法選定と正当性監視のための開発tool**である。

特に次を意図している。

- 理論上速そう，という理由だけで算法を採用しない
- 棄却した算法も再測定可能な形をできるだけ残す
- thresholdをCPU/compiler非依存の普遍値だと思わない
- arbitrary precision / intervalのコストを`double` benchmarkと混同しない
- 大規模Matrixでは算法時間だけでなくrepresentation / allocation / parser costも監視する
- fuzzerで見つけた未知のfailureを，再現可能な固定testへ変換する

mmCalの性能方針は，

> **正当性を確認し，random testで揺さぶり，benchmarkで測り，速いものだけを採用する。**

である。
