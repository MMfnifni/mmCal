# mmCal.Benchmarks

`mmCal.Benchmarks`は，mmCal本体の通常回帰testとは分離した，**性能測定・固定seedランダム不変量試験・大規模Matrix負荷試験・grammar-aware数式fuzzing**のための開発用Console programである。

通常の`mmCal.Tests`や`test_set`が「既知の仕様を壊していないか」を確認するのに対し，`mmCal.Benchmarks`は主に次を担当する。

- BigInt / BigFloat / certified算法の性能変化を測る
- Karatsuba，Toom-3，Burnikel–Ziegler，FFT等の閾値を実測する
- fixed-seed random inputで数学的不変量を検査する
- Matrix，FFT，特殊函数等の数値計算基盤がcertified relationを維持しているか確認する
- 大きなMatrixで時間・メモリ特性を調べる
- 合法な数式を大量生成し，函数同士の未知の組合せ不具合を探索する

benchmark値はCPU，compiler，最適化設定，allocator，OS，cache状態等に依存する。**異なる環境の数値を直接性能保証として比較しないこと。** 閾値を変更する場合は，実際に採用するcompilerとbuild configurationで再測定する。

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

exact MatrixのBareiss / modular crossoverを単独測定：

```text
mmCal.Benchmarks --exact-linear-algebra 1
```

EvaluationBudgetの代表負荷telemetryを表示：

```text
mmCal.Benchmarks --budget-telemetry
```

大きなMatrixを一種類だけ測定：

```text
mmCal.Benchmarks --matrix-large nsvd 64 16
```

代数体compositum / embedding再利用とsame-field reciprocal再利用を単独測定：

```text
mmCal.Benchmarks --algebraic-field 8
```

前半は同一sessionで，`(root[{-2,0,1},2]+root[{-3,0,0,1},1])*(root[{-2,0,1},2]-root[{-3,0,0,1},1])`の初回時間と2回目以降の平均を測定する。後半は12次simple extensionでextended-Euclid reciprocalの初回，cache hit時のwarm reciprocal，warm divisionに加え，minimal polynomialのfirst derivationとwarm cache hitを測定する。persistent field / primitive-element embedding cache / canonical Root materialization / reciprocal reuse / incremental Krylov minimal-polynomial derivationの性能退行を監視する用途であり，絶対性能値そのものを保証するものではない。

certified Gamma / regularized incomplete Betaのprecision scalingを単独測定：

```text
mmCal.Benchmarks --special-functions 1
```

80 / 160 / 320 / 640 / 1280 bitで`gamma[1/3]`と`ibeta[1/3,2/3,1/4]`のcertified 計算基盤を直接測る。Gammaのexact Rational dispatch，static exact Bernoulli table，高精度Stirling planner，ibetaのpoint fast path / shared Beta normalizationの退行監視用である。Stirling-planや定数cacheはprocess内でwarmになるため，完全なcold-start比較では各precisionを別processでも測定する。

random expression fuzzerを有限回実行：

```text
mmCal.Benchmarks --random-expressions --cases 50000 --seed 1234
```

Certification境界だけを集中監査する場合：

```text
mmCal.Benchmarks --certification-boundaries --cases 10000 --seed 1234
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

日常的なcommit前確認より，算法閾値変更，高速化patch，release前の性能監査向けである。

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

## 3.5 `--exact-linear-algebra` / `--budget-telemetry`

```text
mmCal.Benchmarks --exact-linear-algebra 1
```

exact整数行列の係数height 16 / 96 / 256 / 512 bitと複数次数について，`det`と`solveLinear`のBareiss / modular 計算基盤を強制比較する。`inverse`は16 / 96 / 256 bitを別sweepし，automatic採用に十分なcrossoverがあるか確認する。dispatcher 閾値は数学仕様ではなくcompiler・CPU依存の性能policyなので，変更時はこのrunnerを正本とする。

```text
mmCal.Benchmarks --budget-telemetry
```

代数式，modular determinant，modular solve，certified評価，積分の代表入力を1回ずつ評価し，work count，depth，generated node，探索candidate，dense Array / Matrix temporary，BigInt，precision，algebraic refinementに加えて`modular-primes`を表示する。これは既定budget値の校正用telemetryであり，PASS/FAIL判定や数学的意味論を変更しない。

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

Karatsuba / Toom-3 / Burnikel–Ziegler等のcrossoverは数学定数ではなく実装環境依存である。閾値変更時はこの領域を特に比較する。

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

## 4.5 Exact fraction-free QR

exact実数QRはprimitive整数vector上のfraction-free直交化を使い，平方根をQ/Rの最終materializationまで遅延する。full-rankではGram行列のsymmetric Bareiss経路，rank-deficientではdirect fraction-free 代替経路を使う。旧3×3 hard capは撤去済みで，通常benchmarkも2/4/8/16次を測る。

`--exact-linear-algebra`では係数bit長を変えながら，final Rational/radical materialization込みの`inverse/rref/rank/nullSpace/QR`と，rank-deficient 48/64次の構造系pathも測る。

## 4.6 Householder QR block sweep

certified QRについて，column block幅

```text
1, 8, 16, 32
```

を比較する。

現時点ではBigFloat / interval演算コストが支配的で，一つのblock幅が全サイズで安定して勝つとは確認されていない。そのため自動block化の採否はbenchmarkで再評価する前提である。

## 4.7 SVD

certified reduced SVDを，主に4，8，12，16次で測る。

SVDはreconstruction / orthogonality relationを監査する数値計算基盤であり，単純な`double` LAPACK benchmarkとは性質が異なる。

## 4.8 Eigen / Eigensystem

certified approximate eigenvalue / eigensystem 計算基盤を4，8，12，16次で測る。

timing入力は，非対角成分の絶対値を2以下，隣接する対角値の間隔を`4n+1`とした密実対称行列である。Gershgorin円板が互いに素になるためsimple spectrumが構成的に保証され，完全な固有vector基底を持たない入力をtiming failureと混同しない。

計算基盤が結果を返さない場合は，途中までのtiming行を残して`abort`するのではなく，operation，size，digitsをstderrへ出してexit code 1で終了する。

同じ4，8，12，16次の固定入力は`--random-only`でも固有関係と非零vectorを検査する。これにより，長いtiming列を開始する前に入力fixtureと計算基盤契約の不整合を検出する。

一般非正規行列では固有値問題そのものが摂動に敏感であるため，時間だけでなくrandom invariant側のrelation checkと併せて評価する。

## 4.9 FFT

2冪長についてexact FFTとcertified approximate FFTを比較する。

通常：

```text
32, 64, 128 points
```

`--full`では512点まで拡大する。

また非2冪長について，direct DFTとFFT/Bluesteinの実測を比較する。現在のdirect/Bluestein 閾値も，この測定から再評価する前提である。

境界だけを旧policyから独立して再測定する場合：

```text
mmCal.Benchmarks --fft-threshold 1
```

65 / 95 / 127 / 191 / 255 / 257 / 319 / 335 / 351 / 367 / 383 / 384 / 385 / 447 / 509点について，direct DFTと強制Bluesteinを同一入力・16桁で比較し，速い側も表示する。反復数は省略時1。両算法の確定済み十進表示値が成分ごとに一致しなければbenchmark failureとする。現在の保守的policy境界は384点であり，compiler / CPUが変わればこのmodeで再測定する。

GCC Releaseの2026-08-22再測定では383点でBluestein，384点でdirect，385点で再びBluesteinが優位となり，勝敗は境界近傍で単調ではなかった。このため単一crossoverへpolicyを過適合させず，384点を維持する。383/384/385の3点は，今後のcompiler / CPU変更時にこの非単調性も含めて確認するための境界監視点である。

## 4.9.1 Exact Cyclotomic FFT

非2冪exact Cyclotomic 計算基盤だけを短時間で測る場合：

```text
mmCal.Benchmarks --exact-cyclotomic-fft 5
```

5 / 7 / 10 / 12 / 15 / 21点について，最初の`ifft[fft[v]]` round-tripと同一`FourierTransformCache`を使ったwarm平均を表示する。結果は必ず入力Arrayとstructural equalityで一致することをbenchmark内部で確認する。

このbenchmarkは`Q[t]/Phi_n(t)` quotient kernelとfield cacheの退行監視用であり，certified approximate FFTのdirect/Bluestein crossoverとは別に扱う。

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

既知の退化式・恒等式をrandomなsafe-定義域 pointへ適用し，certified enclosureがexact値を包含することを確認する。

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

整数・Rational・certified approximate pathを跨いで監視する。可逆性だけでは完全な固有vector基底の存在を保証しないため，Eigen relationにはdistinct eigenvalueを持つ上三角行列を別生成する。defective / near-defective行列で`eigensystem`が未評価に留まることは仕様どおりであり，failureとして扱わない。

失敗時は実装が出力する`phase` field，1-originの`case`，0-originの`case-index`，size，反例行列を出力する。固定seedなので，同一version・同一case数なら反例を再現できる。

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

v1.5.3の第一段階typed-node化では，同一GCC Release/LTO-offの`--matrix-large transpose 1024 16`で旧variant版`693312 KiB`→typed-node版`299668 KiB`となり，最大RSSを約56.8%削減した。第二段階では`ArrayExpr`を固定1024要素のimmutable packed page + shape/offset/stridesへ移行し，benchmark fixtureも`ArrayBuilder`からexact Rationalを直接構築する。これにより同負荷は最大RSS約`136576 KiB`（133.4 MiB），transpose本体約0.059 msとなった。transposeはデータcopyではなくlayout view生成なので，旧実装のtranspose timingと算法量そのものが異なる点に注意する。

`ArrayBuilder`はpromotionを現在page内へ限定する。1,048,576要素の最後だけsymbolic値にした専用測定でもall-integer版とほぼ同じ約0.30 s / 約69.8 MiBで，先行pageはInteger packedのまま，最終pageだけGenericとなった。これは「最後の1要素で全ArrayをGenericへ作り直す」実装を避けるための設計である。

巨大text inputは別測定する。13.63 MBの1024×1024・10桁decimal literalをCLIへ渡して`dimensions[...]`だけを評価したv1.5.3測定は約4.55 s / 最大RSS `441564 KiB`（431 MiB）。v1.5.2文書の約14.16 MB / 10.9 s / 1.99 GBとは入力textが完全同一ではないため厳密なA/B値には使わない。

---


# 7. Certification Boundary Fuzzer

```text
mmCal.Benchmarks --certification-boundaries
```

`N`のcertified evaluationについて，**分岐切断・極・定義域境界・certified 計算基盤境界**の近傍と，有限precision値に対する**metamorphic変形**を重点的に生成する専用fuzzerである。通常のRandom Expression Fuzzerが式全体のsemantic invariantを広く探索するのに対し，こちらはclosed numeric expressionを中心に，固定分類と変形前後の関係を直接監査する。

```text
Value
DomainError
PrecisionInsufficient       # N::precision
CertifiedBackendUnsupported # N::unsupported
Unevaluated
ResourceLimit
Timeout
OtherError
```

主要な不変条件は次である。

- exact 極 / exact singularityは`DomainError`であり，有限precisionのInformationEnclosureが極を**含み得るだけ**なら`N::precision`である。
- 主値分岐 sideをInformationEnclosureから一意に決定できない場合，片側の値を捏造せず`N::precision`へ戻る。
- 数学的な値が存在するが現certified 計算基盤の外側なら`N::unsupported`であり，`DomainError`へ誤分類しない。
- finite information由来の境界曖昧性はguard digitを増やしても改善しないため，過剰なcertified refinementを消費しない。
- generic `N::unevaluated`，resource exhaustion，case timeoutは通常の数学的分類とは別のfailureとして扱う。
- `N[N[x,p],q]`で外側の要求桁を増やしても，内側の近似値より狭いInformationEnclosureやhidden guard digitsを復元しない。
- 外側の`N`で要求桁を下げる場合，InformationEnclosureは同じか広がるだけであり，狭くならない。
- `+0`，`*1`，二重符号反転等のexact identityは有限precision provenanceを再量子化しない。
- `f[N[x,p]]`は`N[f[x],p]`と同じ表示値になる必要はないが，前者のInformationEnclosureはexact-routeのCertifiedEnclosureを包含しなければならない。
- branch/極曖昧性，finite zeroを含むmatrix pivot判定，FFTや数値微積分の有限precision入力は，nested `N`等の意味保存変形によってhidden exact truthを回収しない。

固定classification probeは59種類で，`log/sqrt/Arg/atan2/Power`，逆三角・逆双曲線，Gamma/digamma/trigamma/zeta，`Ei/Ci/li/polylog`，`1F1/2F1`，`ibeta`，Lambert W，elliptic F等を横断する。さらに36種類のmetamorphic probeで，nested `N`，InformationEnclosure単調性，exact identity，persistent ambiguity，finite-precision線形代数，FFT，数値微積分，`f[N[x,p]]`と`N[f[x],p]`の関係を検査する。各caseではprobeを選んだうえで，要求precisionを`2/5/20/50/100`桁，有限precision入力を`2/5/10/20`桁，exact branch-side epsilonの桁をseedから変化させる。

## 7.1 有限実行

既定は10000 caseである。

```text
mmCal.Benchmarks --certification-boundaries
mmCal.Benchmarks --certification-boundaries --cases 50000 --seed 1234
```

case生成はmaster seedと1-based case番号だけで決まり，thread schedulingには依存しない。

## 7.2 Caseの直接再現

```text
mmCal.Benchmarks --certification-boundaries --seed 1234 --case 817
```

FAIL時にはfamily，probe名，生成式，期待分類，実分類，diagnostic，EvaluationBudget使用量を表示し，同じcommandを`Reproduce`として出す。`--case`では評価開始前に式をflushして表示するため，重いcaseやtimeout候補でも何を実行しているか確認できる。

## 7.3 Timeoutと計算量制限付き

境界fuzzerは各caseをfrontend cancellation token付きで評価し，既定では**2000 ms/case**を超えるとcancelする。

```text
mmCal.Benchmarks --certification-boundaries --timeout-ms 5000
```

これはCoreへwall-clock deadlineを持ち込むものではない。benchmark frontendのwatchdogが既存`EvaluationCancellationToken`を要求するだけであり，timeoutは`DomainError`や`N::precision`とは別の`Timeout` failureとして報告する。

また，`N::precision`が期待されるpersistent ambiguityで4096回を超えるcertified refinementを消費した場合もfailureとする。入力InformationEnclosureが境界を跨いだままなら，global budget近くまでguardを増やす挙動は性能退行とみなす。一方，正常値の特殊函数計算基盤はEuler-Maclaurinや級数評価で多数のwork unitを正当に消費し得るため，fuzzer側から低いglobal `CertifiedRefinement`上限は上書きしない。正常値のcliffはcase timeoutで分離して検出する。

## 7.4 無限loop・複数thread

```text
mmCal.Benchmarks --certification-boundaries --loop --threads 8
mmCal.Benchmarks --certification-boundaries --nostop-loop --threads 8
```

`--loop`は最初のFAILで停止する。`--nostop-loop`はFAILを表示して次caseへ進み，burn-inを継続する。各workerは独立した`KernelSession`を使用し，caseごとに独立したwatchdogを持つ。

`--seed` / `--case` / `--cases` / `--threads` / `--report-every`はRandom Expression Fuzzerと同じ名前を使うが，`--certification-boundaries`が指定されている場合は境界fuzzerへrouteされる。両fuzzerを同一processで同時指定することはできない。

# 8. Random Expression Fuzzer

```text
mmCal.Benchmarks --random-expressions
```

文法と型をある程度理解した**grammar-aware semantic fuzzer**である。

完全なrandom token列を生成するのではなく，原則として評価可能なexact算術，多項式，小型Matrix等を生成し，数学的不変量を検査する。

目的は，

> 固定testでは思いつかなかった函数・構文・簡約の組合せを自動探索すること

である。

## 8.1 有限実行

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

## 8.2 無限loop

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

## 8.3 複数thread

```text
mmCal.Benchmarks --random-expressions --cases 100000 --threads 8
mmCal.Benchmarks --random-expressions --loop --threads 8
```

`--threads N`は**case単位のthroughput並列化**である。1個のSVDや1個の式評価を内部でN thread化するoptionではない。各workerは独立した`KernelSession`を所有し，mutableな履歴・定義・diagnostic stateを共有しない。

各caseはmaster seed + 1-based case番号から独立生成されるため，thread schedulingや`--threads`値が変わっても，同じseed/case番号は同じ式を生成する。FAIL再現時は`--case`を使えばよく，通常1 threadで直接再現される。

worker内ではcase間・比較評価間に`KernelSession::resetForIndependentEvaluation()`を使用する。これは定義・履歴・input/output history・diagnostic・入力番号を初期化するが，通常`reset()`と異なりentropy reseedを行わない。そのため長時間loopで履歴vectorがcase数に比例して成長せず，既存のRNG streamも不用意に変更しない。現generatorはrandom builtinを生成しないため，random builtinを将来追加する場合はcase seedからの明示seed policyを別途定める。

参考として，同一Release/LTO-off build，`--threads 8 --seed 1234`で有限実行した最大RSSは100,000 caseで`48136 KiB`，200,000 caseで`51364 KiB`だった。case数を2倍にしても履歴保持に相当する線形増加は見られない。これは長時間`--loop`の完全な定常性証明ではないため，overnight burn-inでは引き続きprocess RSSを監視する。

## 8.4 FAIL時

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
Budget    : work=... depth=... nodes=... simplify=... solve=... integrate=...
Reproduce : mmCal.Benchmarks --random-expressions --seed 92847163 --case 48172
```

FAIL後は簡易shrinkingを行い，可能ならより小さい反例を`Reduced`として表示する。

通常の`--loop`またはfinite modeでFAILした場合，processはexit code 1で終了する。`--nostop-loop`ではFAILは継続中の反例として扱われ，その場ではprocessを終了しない。

## 8.5 Caseの直接再現

```text
mmCal.Benchmarks --random-expressions --seed 92847163 --case 48172
```

case番号は**1-based**である。

各caseはmaster seedとcase番号から独立したcase seedを導出して生成する。このため，48172件目を再現するために1～48171件目を再実行する必要はない。

単一case再現ではFAILだけでなくPASS時も`Depth`，`Invariant`，`Expression`を表示する。数学的不変条件は通るが異常に遅いcaseを，生成式を失わず直接診断できる。

長時間`--loop`を回す上で重要な仕様である。

## 8.6 Depth

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

## 8.7 Progress表示

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

## 8.8 現在のInvariant

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

### 微積分 derivative-back

random polynomial `p`について，次をsymbolic residualとして検査する。

```text
fullSimplify[D[integrate[p,x],x]-p] == 0
```

### Solve

異なる二つの整数根から方程式を構成し，`solve`が返した`SolutionSet` bindingを順序非依存で比較する。Formatterのbranch表示順をcorrectness条件にはしない。

### Matrix inverse

非零対角を持つ小型上三角整数行列を構成し，generator段階で非特異性を保証した上で検査する。

```text
dot[A,inverse[A]] == I
```

### FFT

1～12点のexact Gaussian整数vectorについて，2冪・非2冪計算基盤を跨いで検査する。

```text
ifft[fft[v]] == v
```

まず`fullSimplify`後のcanonical構造一致を検査する。generic exact DFTでroot-of-unityの相殺形が残る場合は，`actual-expected`を先に作らない。exact round-trip `Expr`の各成分を`CertifiedEvaluator`へ直接渡し，50桁のcertified enclosureが元のexact Gaussian整数点を含み，さらに実部・虚部ともその点から絶対`10^-40`以内に収まることを確認する。これにより`(-I/2+sqrt[3]/2)^3`や相殺する`sqrt[3]`線形項を含む式で，subtractive residual側だけがbounded refinementを使い切る偽FAILを避ける。Formatterによる文字列化・再parseもoracleへ挟まず，近似FFTへの置換ではなくexact FFT結果そのものを検証する二段目oracleである。

二段目oracle自体もboundedである。`CertifiedEvaluator`は実再帰depth budgetを持ち，巨大なgeneric exact DFT式をWindowsのnative stack overflowまで再帰させない。oracleがこのdepth/work上限で判定できなかったcaseはFFTの数学的FAILとは区別して`inconclusive`として計数し，`--nostop-loop`の進捗表示にも件数を出す。単一`--case`再現では`[INCONCLUSIVE]`を表示する。

### 定義域 / branch境界

`1/0`，`0^0`，`cot[0]`，`atan2[0,0]`，`atanh[1]`が必ず`DomainError`へ分類されること，`sqrt[-1] == I`，`sqrt[(-3)^2] == 3`となる主値 square-root 分岐を検査する。`ResourceLimitError`や一般exceptionへの分類退行もFAILである。

FAIL時の`Budget`行は，直前評価の`EvaluationUsage`である。反例が資源上限そのものなのか，低い使用量で生じた意味論failureなのかをseed / caseと同時に判断できる。

---

# 9. Fuzzerと通常random checkの違い

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

# 10. Benchmark結果の読み方

## 10.1 1回の結果だけで閾値を変えない

CPU boost，background process，cache，allocator等で数%程度は容易に変動する。

閾値変更では，

1. 同じRelease binary
2. 同じmachine
3. 同じ電源設定
4. 同じoperand範囲
5. 複数回測定

を基本とする。

僅差なら複雑な算法を採用しない。

mmCalでは実際に，理論上有望でもbenchmarkで優位性が得られず棄却・保留した最適化が複数存在する。

## 10.2 correctnessとspeedを混同しない

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

## 10.3 exactとapproximateを分けて読む

例えば，

```text
exact det
N[det,16]
```

は同じ数学函数でも計算基盤が異なる。

exact pathはBigInt / Rational / symbolic expressionを維持し，approximate pathは`ApproximationContext`からBigFloat / intervalへ直接dispatchする場合がある。

したがって単純に「Nを付けると何倍速い」と一般化せず，どの計算基盤を測っているかを見る。

---

# 11. 推奨workflow

## 通常の変更

```text
mmCal.Tests
black-box tests
mmCal.Benchmarks --random-only
```

## 数値計算基盤を変更した場合

```text
mmCal.Benchmarks --random-only
mmCal.Benchmarks --benchmark-only
```

## BigInt算法閾値を変更した場合

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

# 12. Exit code

| code | 意味 |
|---:|---|
| `0` | 指定した試験・benchmarkが正常終了 |
| `1` | random invariant / random-expression fuzzer等でFAIL |
| `2` | option不足，未知option，矛盾したoption等のCLI error |

benchmark中の内部invariant違反等ではprogramがabortする場合がある。これはtiming値を誤って正常結果として採用しないためである。ただし，certified Eigen / Eigensystem 計算基盤が結果を返せない場合は診断を表示し，code 1で終了する。

---

# 13. 現在の設計上の注意

`mmCal.Benchmarks`は一般的なbenchmark frameworkではなく，**mmCal固有の算法選定と正当性監視のための開発tool**である。

特に次を意図している。

- 理論上速そう，という理由だけで算法を採用しない
- 棄却した算法も再測定可能な形をできるだけ残す
- 閾値をCPU/compiler非依存の普遍値だと思わない
- arbitrary precision / intervalのコストを`double` benchmarkと混同しない
- 大規模Matrixでは算法時間だけでなくrepresentation / allocation / parser costも監視する
- fuzzerで見つけた未知のfailureを，再現可能な固定testへ変換する

mmCalの性能方針は，

> **正当性を確認し，random testで揺さぶり，benchmarkで測り，速いものだけを採用する。**

である。

---

# 14. 特殊函数 performance-cliff 監査

特殊函数の計算量制限付き境界と精度増加時の計算量を同じ尺度で比較するには，次を使う。

```text
mmCal.Benchmarks --performance-cliffs [iterations]
```

各pointについて公開`KernelSession`経由の評価時間と`EvaluationUsage`を表示する。`iterations >= 2`では最初の評価を`first-ms`，以後の平均を`warm-ms`として分離する。static table/cacheのfirst-use taxと定常costを混同しないためである。

主な列は次の意味を持つ。

- `first-ms` / `warm-ms`: 公開評価経路全体のwall time
- `cert-ref`: `CertifiedRefinement`
- `work`: evaluator work count
- `nodes`: generated Expr node
- `CLIFF(time ...)`: 同一sweep内の隣接成功点に対してwall timeが既定4倍以上へ増えた候補
- `CLIFF(work ...)`: wall timeより先にcertified workが4倍以上へ増えた候補
- `BOUNDARY`: numericからheld/error等へ結果分類が切り替わった点

`BOUNDARY`はperformance cliffではない。例えば分岐切断や極を跨ぐため現計算基盤が未対応の入力を即時に未評価へ戻すことは，病的な反復へ入らないための正常な計算量制限付き policyである。固定閾値は算法改善後に再監査し，数学的理由がなくなったものは撤去する。

現runnerは，Gamma，LogGamma，erf/erfc，Lambert W，Fresnel，`1F1` / `2F1`，楕円積分，`Ei/Si/Ci/li`，polylog，zeta，digamma/trigamma，Beta/ibetaを一巡する。Lambert Wは`-1/e`に対する`10^-2`～`10^-160`のcomplex/real offsetと20/40/80/100桁sweepを持ち，平方根branch-point 計算基盤の浅い近傍・極端な近傍の双方を監視する。楕円積分は`m/n`境界だけでなく`m/n=9/10`でのamplitude sweepと10/20/40/80桁precision sweepも持ち，small-amplitude fast pathの改善後にcomplete case近傍へcliffが移っていないかを継続監視する。real/complexの代表計算基盤もsnapshotとして測る。

2026-08-25のGCC Release / LTO offでの初回横断監査では，次が主要hotspotだった。値はalgorithm選定用のsignalであり，他CPU/compilerへ固定閾値として移植しない。

| workload | warm代表値 | 主なsignal |
|---|---:|---|
| `N[1F1[1/2,5/4,160],20]` | 約1.08 s | `z=64 -> 128`で約7.8倍，境界160までは増加 |
| `N[polylog[2,49/50],20]` | 約139 ms | unit-circle境界へ近づき連続的に増加 |
| `N[ellipticF[1/2,9/10],20]` | 約162 ms | `m`増加に伴い級数cost増大 |
| `N[ellipticE[1/2,9/10],20]` | 約206 ms | 同上 |
| `N[ellipticPi[9/10,1/2,1/3],20]` | 約307 ms | characteristic境界側が特に重い |
| `N[fresnelc[7+I],20]` | 約1.97 s | `\|z\|<8`境界直前で急増 |
| complex `2F1` 解析接続代表 | 約0.98 s | `1/z` 解析接続経路が重い |
| `N[li[-2+I],20]` | 約1.23 s | complex `Log -> Ei`経路が重い |
| `N[zeta[3/2],20]` | 約155 ms | 約4.2万 certified refinement |
| complex `zeta`代表 | 約0.45–0.55 s | refinement数だけでは説明できない高いinterval算術cost |

一方，`gamma[1/3]`，`ibeta[1/3,2/3,1/4]`，20桁の通常`2F1`はprecisionを10→80桁へ上げても今回の4倍判定で急なcliffは出ず，ほぼ連続的に増加した。

### polylog 49/50境界撤去後の再監査（2026-08-28）

旧`|z|<=49/50`はexact-Rational級数のrepresentation cliffを避けるためのpolicyだった。direct termをoutward interval recurrenceへ移し，`Li_2`へDLMF connection formula，order 3～12の正実unit-point近傍へ`mu=log(z)`展開を追加した後，固定閾値を撤去した。

GCC Release / LTO off，20桁，`--performance-cliffs 1`では`Li_2(z)`の`z=1/2,4/5,9/10,19/20,49/50,99/100,199/200,999/1000` sweepが約1.3～3.7 msで推移し，旧49/50位置に不連続はない。`z=999/1000`は約3.0 msである。runnerにはこのargument sweepに加え，10/20/40/80桁precision，`I/-2/2+I`等の解析接続，`Li_3/Li_4/Li_8` near-one sweepを常設する。
### complex Fresnel / large positive 1F1 最適化後

上表の初回値を受けて，complex Fresnelのtail majorantをexact Rational反復からoutward intervalへ変更し，large positive `1F1`のterm/partial sumもguard付きintervalへ移した。この段階では公開計算量制限付き境界を`Fresnel |z|<8`，`1F1 |z|<=160`のまま維持していた。

同じGCC Release / LTO off環境の再測定では，`N[fresnelc[7+I],20]`は約14 ms，`N[1F1[1/2,5/4,160],20]`は約6 msまで低下した。後続のcapability監査では1F1のinterval級数を符号非依存majorantへ一般化し，real/complexとも固定`160`境界を撤去した。さらにcomplex Fresnelも軸近傍大引数へ保証付き`f/g`漸近計算基盤を追加し，旧`|z|<8`境界を撤去した。runnerは`1F1 magnitude scaling`，`z=512` / complex `256+I` precision sweepに加え，complex Fresnel `32+I` precisionとargument sweepを常設する。


### complex 2F1 解析接続 / li の再測定（2026-08-26）

後続監査で，両者の共通ボトルネックが主値 complex Log/Arg内部のcertified `atan` exact-Rational級数であることを特定した。atanの隣接部分和をoutward `RealInterval`で保持し，2F1側ではprovably-real Gammaをreal 計算基盤へ分流，complex GammaのStirling remainderも高precision dyadic下端のexact冪乗を避けた。

同じGCC Release / LTO off系統の3回warm測定では，代表complex `2F1` 解析接続は約0.96 sから約100 ms，`li[-2+I]`は約1.03 sから約13.5 msへ低下した。`li[-2]`も約575 msから約9 msまで低下した。runnerには`2F1 continuation precision`と`complex li precision`を追加し，10/20/40/80桁の伸びを継続監視する。

---

# 15. Complex Root performance-cliff 監査

Complex Algebraic Rootの全根分離と，既存isolating diskからの単根refinementを分けて測るには次を使う。

```text
mmCal.Benchmarks --algebraic-root-cliffs [iterations]
```

出力列は`create-ms`が初期の全根isolation＋決定的ordering，`refine80-ms` / `refine320-ms`が選択済みrootだけのcertified local refinementである。対称seed pathologyを拾う`x^n-2`系列，generic dense 5次，近接複素根を持つ4次caseに加え，根半径が`2^-10`と`2^10`へ分離するtwo-radius 16次caseを常設する。

2026-08-28のGCC Release / LTO offでは，Newton-polygon multi-radius seed導入時点で`x^8-2`約68 ms，`x^16-2`約0.91 sを維持し，two-radius 16次caseを旧単一Cauchy半径seed約10.79 sから約0.36 sへ低下させた。その後の高次数監査でRouché証明をdirected BigFloat intervalの十分条件へ先行移行し，閉じない場合だけexact Rational Taylor boundへ切り替える構成へ変更した。Durand-Kernerも候補補正が十分小さくなればdegree回数を待たず終了し，最終Rouchéが失敗した場合だけprecisionを上げて再試行する。最終runnerでは`x^16-2`約0.16 s，一般Solveは16次約0.09 s，32次約0.36 s，64次約1.56 s，65次約1.52 sである。96次までの追加監査でも性能が連続したためRoot defining-polynomial budgetは64から96へ再設定した。今後の候補生成改善としてAberth–Ehrlich法を比較候補に残すが，採用するとしてもapproximate seed/refinementに限定し，Rouché certificationを置換しない。

### real Ei / Si / Ci の大引数監査

旧`|x|<=96`境界撤去後は，`--performance-cliffs`のreal `Ei/Si/Ci` sweepで96/97だけでなく512/1000/10000を継続監視する。moderate argumentのTaylorと大引数asymptoticの切替は固定magnitude contractではなく，漸近剰余が要求幅へ届くか，Taylor term cap，共通`EvaluationBudget`で決まる。`Ci[1000]`等がEulerGamma cancellation pathへ戻ること，またpositive `Ei`の巨大値でabsolute-error基準へ退行することをregressionとして検出する。

### complex Fresnel の固定 `|z|<8` 境界撤去

Maclaurinのexact-majorant cliffを除去した後に旧`8`境界を再監査すると，`8+I`から`20+I`まではそのままcertifyできる一方，`32+I`付近では振動的な巨大中間項の相殺によりseriesが再び秒級へ伸びた。固定境界を別の値へ移す代わりに，DLMF 7.12の複素Fresnel補助函数`f/g`漸近展開を追加した。90度回転対称性により入力を正実軸近傍へ写し，`|tan(arg z)|<=1/4`を区間端点から証明できる場合だけ，最初の未使用項で剰余を押さえられる安全なwedgeとして漸近計算基盤を使う。それ以外はentire Maclaurinへ切り替える。

2026-08-27 GCC Release / LTO off，3回warmでは`32+I`が約10 msで，10/20/40/80桁は約8/10/16/35 msだった。`1+32I`もquarter-turn経路で約12 ms，`8+8I`のような対角方向はseriesのまま約13 msである。旧`|z|=8`は現在の計算基盤 capabilityを表さないため削除した。


### complex Ei/Ci の固定512/128境界撤去

`--performance-cliffs`は旧complex `Ei |z|<=512` / `Ci |z|<=128`位置の前後と10/20/40/80桁precision sweepを常設する。大引数はcertified `E1`漸近＋主値 connection，漸近証明が閉じない高precision入力はguard付きseriesへ切り替える。分岐切断を跨ぐfinite-precision enclosureは性能境界ではなく`N::precision`分類として別に監視する。

2026-08-28 GCC Release / LTO offの20桁代表値は`Ei[512+I]`約8.2 ms，`Ei[513+I]`約9.4 ms，`Ci[128+I]`約44.2 ms，`Ci[129+I]`約45.7 msで，旧境界位置に不連続はない。`Ei[1000+I]`は約9.2 ms，`Ci[1000+I]`は約18.4 ms。`Ci[140+I]`は高precisionでseries 代替経路を使うため，20/50/100桁の時間も併せて監視する。
