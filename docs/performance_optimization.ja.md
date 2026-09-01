# mmCal 高速化・算法選定記録

この文書は，各バージョンで行った高速化について，**何を採用したか，何を比較したが棄却したか，なぜその判断をしたか**を残すための記録である。

速度だけでなく，mmCalの基本契約である

- exact arithmetic
- 主値分岐 / definedness / 定義域
- directed rounding
- certified enclosure
- formatter round-trip safety

を壊さないことを採用条件とする。

benchmark値はCPU/compiler/allocator/cacheに依存する。以下は開発時のGCC系環境で得た代表値であり，普遍的な性能保証ではない。Visual Studio / MSVCを含む別環境では`mmCal.Benchmarks`で再測定する。

---

## 1. 採用判断の原則

高速化は次の順で判断する。

1. 旧実装と数学的意味論が一致すること
2. fixed-seed random / invariant testで差分がないこと
3. 境界値・extreme valueで正しさを維持すること
4. 実測で十分な利益があること
5. 小サイズを犠牲にする場合は閾値 dispatchで隔離できること
6. コード複雑性に見合わない微小改善は採用しないこと

したがって「理論上漸近的に速い」「有名libraryでも使われる」だけでは採用理由にしない。

---

# 2. BigInt multiplication

## 2.1 Schoolbook — 維持

小さいoperandでは二重loopのschoolbook multiplicationが最も低overheadだったためbase caseとして維持した。

## 2.2 Karatsuba — 採用

初期閾値 sweepでは8～16 limbsからKaratsubaへ入れると明確に退行した。一方，32～48 limbs以降で利益が安定した。

v1.5.1既定値:

```text
Karatsuba crossover ≈ 48 limbs
1 limb = 32 bit
```

極端にunbalancedなoperandはschoolbookへ戻す。

代表値:

|   operand | schoolbook | adaptive | speedup |
| --------: | ---------: | -------: | ------: |
| 128 limbs |   約9.5 µs | 約6.9 µs |  約1.4x |
| 256 limbs |    約38 µs |  約21 µs |  約1.8x |
| 512 limbs |   約151 µs |  約63 µs |  約2.4x |

factorialにもproduct tree経由で波及し，代表測定では

```text
10000!  約4.4 ms → 約2.0 ms
20000!  約19.5 ms → 約6.8 ms
40000!  約88 ms → 約23 ms
```

程度まで改善した。

## 2.3 Toom-3 — 採用

Karatsubaよりさらに巨大なbalanced multiplicationではToom-3が有利になった。

v1.5.1の代表閾値:

```text
top-level Toom-3   ≈ 1280 limbs
recursive Toom-3   ≈ 448 limbs
```

512～1024 limbs付近ではToom-3のevaluation/interpolation overheadが勝つcaseがあるため，早過ぎるdispatchは避けた。

代表値ではKaratsuba-only比で

```text
4096 limbs  約1.17x
6144 limbs  約1.3x
```

程度の改善が得られた。

## 2.4 Toom-4 / FFT / NTT — 保留

v1.5.1では未導入。

Toom-3よりさらに巨大な領域では候補になるが，interactive用途でcrossoverが十分に現れるまで複雑化しない。将来は`mmCal.Benchmarks`でToom-4とFFT/NTTの閾値を測って判断する。

---

# 3. Dedicated square

`x*x`は一般乗算と異なりcross termが対称なので，専用squareを採用した。

- small: symmetric schoolbook square
- large: Karatsuba square

代表値:

```text
512 limbs   一般乗算 約63 µs → square 約38 µs
1024 limbs  一般乗算 約189 µs → square 約111 µs
```

`pow`のrepeated squaringへ直接効く。

## Toom-3 square — 棄却

専用Toom-3 squareも実装して比較したが，現在の閾値域ではKaratsuba squareより遅かった。

理由:

- evaluation/interpolation cost
- square専用Karatsubaが既に対称性を十分利用している
- crossoverが一般Toom-3 multiplicationより後ろへ寄る

そのためv1.5.1では採用していない。

---

# 4. Karatsuba workspace

temporary vector allocationを減らす目的で二方式を比較した。

1. vector pool型
2. recursion depthごとのscratch型

結果は512～1024 limbs付近で最大約5～10%の退行。

allocator削減より

- pool管理
- resize/clear
- scratch管理
- cache locality

のcostが勝ったため棄却した。

MSVC allocatorでは結果が変わる可能性があるため，将来再測定は可能。

---

# 5. Factorial

## Balanced product tree — 維持

既存の`productRange()`はbalanced treeで同程度の大きさのBigInt同士を掛けるため，adaptive multiplicationとの相性が良かった。

追加した軽量最適化:

- machine integer leafをdecimal string経由でparseせず`BigInt::fromUnsigned()`へ直接構築
- 1-limb multiplicationを`multiplySmall()`へ落とす

## Prime-Swing — 棄却

Prime-Swingも実装し，初版で見つかったprime exponent重複計算等も修正した上で比較した。

それでも代表値は

```text
320000!
product tree  約0.59 s
Prime-Swing   約1.5 s
```

で，現BigInt 計算基盤では既存product treeが速かった。

したがって「一般に高級な算法だから」という理由では採用しなかった。

将来，prime処理や乗算計算基盤が変われば再評価対象。

---

# 6. Division

## 6.1 Knuth normalized long division — 維持

小～中サイズでは低overheadで，Burnikel–Zieglerのbase caseとして優秀なため削除しない。

## 6.2 Burnikel–Ziegler — 採用

巨大でbalancedなdivisionに導入。

GCC環境では約32 limbs以降で利益が安定した。商が小さいcaseはKnuthへ残す。

代表値:

```text
1024-limb divisor/quotient  約1.18 ms → 約0.12 ms
2048-limb divisor/quotient  約5.0  ms → 約0.35 ms
```

波及対象:

- `/`, `%`
- Euclidean GCD
- Rational normalization
- integerSqrt / integerCubeRoot
- decimal divide-and-conquer
- BigFloat division

## 6.3 Power-of-two division — 採用

`x / 2^k`を一般divisionへ流さず，

```text
quotient  = x >> k
remainder = low k bits
```

とする。

4096 limbs級でms級から数µs級まで短縮した。

---

# 7. GCD

## Euclidean GCD — 維持

Burnikel–Zieglerの高速化をそのまま利用できるため，既存Euclidean `%`を維持した。

## Binary GCD — 棄却

Stein binary GCDを実装・比較したが，現BigIntでは約4～30倍遅いcaseがあった。

shift/subtraction回数と大きなtemporary処理が多く，現在の高速divisionを使うEuclidに勝てなかった。

次候補はLehmer GCD。

---

# 8. Decimal conversion

巨大factorialの計算本体を高速化すると，次にbinary-limb→decimal conversionが支配的になった。

## 8.1 `10^9` chunk — 採用

旧実装の`/10` 1桁ずつを，`/10^9` 9桁ずつへ変更。

`40000!`で約6秒から約0.5秒級まで改善。

parse側も9桁chunk化し，巨大decimal parseを大幅に短縮した。

## 8.2 Divide-and-conquer decimal conversion — 採用

さらに巨大な`10^(9*2^k)`でほぼ半分へ分割して再帰変換する。

`40000!`（約166,714桁）の代表値:

```text
旧 digit-wise        約6.0 s
10^9 chunk           約0.53 s
D&C                  約0.19 s
```

元実装比で約30倍級。

formatter固有overheadはこの規模でも小さく，主costはBigInt→decimal自体だった。

---

# 9. `tryToUint64`

巨大BigIntを一度10進文字列へ変換してから`from_chars`で失敗する経路を修正した。

64bit超がbit lengthで明白なら即時`nullopt`。

100000-bit級で十数msから測定限界近くまで短縮した。

---

# 10. BigFloat extreme exponent gap

旧加算は常にcommon exponentへexact alignmentし，`1 + 2^-5000000`でも巨大left shiftを作った。

v1.5.1では，要求precisionと符号から結果を一意に証明できるcaseだけfast pathへ入れる。

全rounding mode:

```text
NearestEven
TowardPositive
TowardNegative
TowardZero
```

を扱い，曖昧な場合は旧exact alignmentへ切り替える。

代表値:

```text
53-bit, 1 + 2^-5,000,000
約1.4 ms → sub-µs級
```

速度のためにsmall operandを単純破棄しているわけではない。

---

# 11. Pi

## Machin formula — 旧referenceへ降格

v1.5.0のMachin公式は保証構造が単純だったが，10000桁級で約12秒まで伸びた。

## Binary-splitting Chudnovsky — 採用

代表値:

```text
N[Pi,10000]
約12.3 s → 約0.4 s
```

高桁既知値と照合し，certified enclosureの契約を維持した。

---

# 12. `exp` / `E`

## 逐次RealInterval Taylor — 旧referenceへ降格

高桁でinterval object生成と巨大Rational中間値が支配的になった。

## Binary splitting + certified 値域 reduction — 採用

小さいRationalはexact binary splitting，大きい分子・分母はfixed-precision interval binary splittingへ切り替える。

代表値:

```text
N[E,5000]  数秒級 → 約0.1 s級
```

一般Rationalでも深めの値域 reductionを行い，中間項の成長を抑える。

---

# 13. `log`

## 逐次atanh級数 — 旧referenceへ降格

`N[log[2],3000]`が約6～7秒まで伸びた。

## Binary splitting + sqrt 値域 reduction — 採用

- `log[2]`等の小係数: exact binary splitting
- 高bit mantissa: certified sqrtを複数回行って1近傍へ縮約し，interval binary splitting

代表値:

```text
N[log[2],3000]                 約6.7 s → 約0.1 s級
N[log[123456789/987654321],5000]  約22 s → 約3 s級
```

さらに高桁ではbit-burst / AGM系が将来候補。

---

# 14. 巨大Radianの三角函数

旧経路では巨大な生Radian Rationalを小区間へ落とさずTaylor評価し，

```text
N[sin[10^6],20]
```

が5秒timeout級になった。

## Certified argument reduction — 採用

`double fmod`は使わず，Piの保証区間から`x/(Pi/2)`の象限integerを一意に証明し，小区間へ縮約する。

point Rationalだけでなく，`10^6 sqrt[2]`のような保証区間入力にも適用する。

変更後は`10^6`～`10^12`級の代表入力が対話的時間へ戻った。

---

# 15. FFT plan cache

exact radix-2 FFTでは同じtransform sizeでbit-reversalとtwiddleを何度も構築するため，Evaluator/session内でplanをcacheする。

cache対象は入力結果ではなく，size依存のplan情報。

process-globalにはせずSymbol/session lifetimeを安全に保つ。

---

# 15.4. exact Cyclotomic FFT —

## 旧 generic `cis` / Expr direct DFT — 代替経路へ降格

以前の非2冪exact FFTは，`radix2Transform()`から`directTransform()`へ代替経路し，各twiddleを`cis[-2 Pi k/n Rad]`のgeneric `Expr`として構築していた。5/7/10/12点程度でもforwardでroot-of-unity式が増え，inverseでは同じcyclotomic恒等式をSimplifierへ再証明させるため，`ifft[fft[v]]`が巨大式になり`tester.py`の`test5_matrix.txt`を支配していた。

旧dispatchとgeneric DFT本体は削除せず，symbolic入力や新計算基盤のbudget外で使う代替経路として維持する。dispatch旧形は変更理由付きコメントで`signal_processing.cpp`に隣接保存する。

## canonical `NumberFieldContext`でζ_nを持つ案 — 棄却

最初の試作では一般algebraic 計算基盤をそのまま利用し，primitive root `zeta_n`をcanonical Complex Rootとして`NumberFieldContext`へ載せた。しかしFFTで必要なのは同じroot-of-unity field内の線形演算であり，generatorのroot isolation / canonical Root materializationが支配的になった。試作では5点forwardが約2.5 s，7点forwardが約19.8 sまで悪化し，ζ_7のContext初回構築だけでも約0.4 s級だったため棄却した。

この結果から，FFT内部ではembedded algebraic numberとしてのroot identityを毎回要求せず，quotient ring/field arithmeticだけを使う方針へ変更した。これはSageのexact DFTがCyclotomicField上の列を扱う設計や，FLINTがcyclotomic polynomial quotient上で先に計算してからroot of unityへ評価する実装方針とも整合する。

## `Q[t]/Phi_n(t)` quotient 計算基盤 — 採用

`CyclotomicFieldContext`はembedding/root isolationを持たず，次だけを保持する。

- conductor `n`
- exact cyclotomic polynomial `Phi_n(t)`
- `t^degree mod Phi_n(t)` reduction
- `t^k`のpower-basis座標

5点以上の非2冪exact Rational入力はこの座標へ直接写し，DFTの加減乗算をRational vectorだけで行う。Gaussian Rationalを含む場合は必要に応じconductorを`lcm(n,4)`へ拡張し，`I=t^(3n/4)`としてexactに埋め込む。出力境界だけで従来表示と連続する単一generator `cis[-2 Pi/n Rad]`へ戻す。

現budgetは`phi(n)<=64`。budget超過，またはsymbolic expressionを同じquotient fieldへ証明付きで写せない場合は旧generic exact DFTへ切り替える。近似/certified FFT経路は変更しない。

GCC Release / LTO offの専用`--exact-cyclotomic-fft 5`代表値：

| length | first round-trip | warm round-trip |
| -----: | ---------------: | --------------: |
|      5 |        約0.72 ms |       約0.49 ms |
|      7 |        約1.79 ms |       約1.74 ms |
|     10 |        約1.50 ms |       約1.33 ms |
|     12 |        約1.44 ms |       約1.27 ms |
|     15 |        約9.11 ms |       約9.39 ms |
|     21 |        約41.1 ms |       約37.3 ms |

`tester.py --timings`では`test5_matrix.txt`が約1.46 sから約0.38 sへ低下した。

## 2冪inverseの円分体再埋込み — 採用

2冪長のforward FFTは既存radix-2表現の方が読みやすいため変更しない。一方，16点以上の`ifft[fft[v]]`ではradix-2が生成した`cis[rational Pi]`，Rational scale，`sqrt[2]`等を同じ`Q(zeta_n)`のpower-basis座標へ再認識し，inverseだけを座標上で実行する。純粋な数値spectrumやmembershipを証明できない式には適用せず，従来経路へfallbackする。

2冪では`Phi_(2^m)(t)=t^(2^(m-1))+1`なので，`multiplyByPower()`は一般多項式積を使わず係数shiftと符号反転で`t^k`倍を処理する。forward出力は16点で旧実装とbyte-for-byte一致し，GCC RelWithDebInfoの統合確認ではexact round-tripが16点約0.01 s，32点約0.08 s，64点約0.37 s，128点約6.1 sで元vectorへ完全復元した。128点はなお重いが，旧generic Exprの式爆発は解消している。

## mixed-radix常設化 — 保留

現在の非2冪cyclotomic transform本体は座標上のdirect O(n^2) DFTである。21点でもround-trip約34 msで，今回の主問題だったgeneric Expr explosionは既に解消した。mixed-radix Cooley–Tukeyやprime長Rader/Bluesteinをexact quotient座標へ追加することは可能だが，実装複雑度を増やす前により大きい対応長でcrossoverを測る。

# 15.5. precision-aware `N` と certified FFT

## 旧経路 — exact FFT完成後に`N`を適用

従来の`N`は引数を通常評価してから呼ばれていたため，

```text
N[fft[data],16]
```

でもまず巨大なexact Fourier式を構築し，その後に各成分を近似していた。FFT本体が`Expr`のexact multiply/add/Simplifierをbutterflyごとに通るため，近似値しか要らない場合にもsymbolic costを全額支払っていた。

## precision-aware evaluation — 採用

`N`の第1引数を保持し，precisionを先に確定する。FFT dispatch時にprecision contextが存在すれば，`double`ではなく`ComplexInterval`/BigFloat端点で直接transformする。exactな`fft[...]`の経路は変更しない。

代表benchmark（同一GCC Release環境。当時の旧`N`仕様で16 fractional digits）:

```text
32 points   exact ~11.6 ms   certified ~2.7 ms
64 points   exact ~63.5 ms   certified ~6.2 ms
128 points  exact ~327.6 ms  certified ~13.1 ms
```

非2冪のcertified FFTは`--fft-threshold`でdirect DFTと強制Bluesteinを同じ入力・16桁・指定反復数で比較する。GCC Release再測定では319点がdirect 1533 ms / Bluestein 1609 ms，335点がdirect 1714 ms / Bluestein 1638 msで，crossoverはこの間だった。一方，MSVC `--full`では257点がdirect 1823 ms / Bluestein 5002 ms，509点がdirect 14007 ms / Bluestein 7896 msだった。primaryのMSVC環境へ保守的に寄せ，現在は384点未満をdirect，それ以上をBluesteinへ送る。これは数学定数ではなく環境依存のpolicy値である。

専用sweepは65 / 95 / 127 / 191 / 255 / 257 / 319 / 335 / 351 / 367 / 383 / 384 / 385 / 447 / 509点を測る。旧policyを通る`evaluateApproximateFft`同士の比較ではなく，directとBluesteinを明示的に強制する。両結果は各実部・虚部の確定済み十進表示値をexact比較し，`-12`と`-12+0…I`のようなzero-component表現差だけを無視する。

2026-08-22のGCC Release再監査では，383点がdirect 1990 ms / Bluestein 1437 ms，384点が1001 / 1434 ms，385点が2010 / 1414 msとなり，境界近傍の勝敗が単調ではなかった。certified direct側のargument reduction・refinement回数や長さの算術構造が定数項へ効くため，単一の局所crossoverだけに合わせて閾値を335等へ動かすことは過学習になる。primary MSVC実測も考慮し，policyは384を維持する。benchmarkには383/384/385を恒久的に含め，環境変更時に再監査する。

同じ監査でcertified特殊函数も「数学的収束」と「実用的計算量制限付き」を分離した。GCC Releaseでは`polylog[2,0.999]`，`2F1[...,0.98]`，`ellipticF[...,0.98]`が収束域内でも現exact-majorant seriesでは数秒から10秒超へ急増した。一方`polylog`は`|z|<=49/50`，`2F1`/ellipticは`<=9/10`，`1F1`は当時`|z|<=160`，`Ei/Si/Ci`は96までを保守的な計算量制限付き policyとした。固定計算基盤範囲やplanner上限を超えた失敗は`CertifiedBackendUnsupported`とし，guard precisionだけを増やす無意味なretryを禁止する。

採用理由:

- exact-first APIを維持したまま近似要求だけを高速化できる
- machine `double`を導入せず任意精度・外向き丸めを維持できる
- `N`のprecision伝播は将来ほかの高cost builtinにも再利用できる
- 非2冪のO(N^2)崖をapproximate pathではBluesteinで回避できる

Cyclotomic quotient 計算基盤導入以前はexact FFT自体のsymbolic expression explosionを意図的に残していた。現在は対応非2冪exact入力を上記Cyclotomic quotient 計算基盤へ送ることでこの問題を解消し，budget外・symbolic membership未証明だけgeneric Expr 代替経路を維持する。

---

# 15.6. Array / Matrix

## flat Array + exact Number 計算基盤 — 採用

v1.5.2では行列専用のnested containerを増やさず，既存Arrayの`shape + row-major flat storage`を基盤とした。`MatrixView`はArrayをzero-copy参照し，Gaussian / Gauss-Jordan等で書換えが必要な場合だけflat `MatrixBuffer`へ複製する。これはv1.5.2 release時点の設計であり，v1.5.3ではpersistent Arrayのphysical storageをimmutable paged backing + stride viewへ置換している。

exact Number行列ではpivot loopからExpr生成とSimplifier呼出しを外し，`Number`を直接累積・消去する。`dot`も全要素がNumberならcellごとの積和を`Number`だけで処理し，symbolicの場合のみExprを構築する。

一般symbolic `det` / `inverse`はLaplace/adjugate展開の仕事量を共有budgetで制限する。三角行列は次数に依存せず対角積へ落とし，疎行列はbudget内なら処理するが，dense高次行列は階乗級の式を作る前に未評価へ戻す。

## precision-aware certified Matrix — 採用

`N[det[A],p]`等はexact結果を完成してから近似せず，FFTと同じ`ApproximationContext`を受けて`ComplexInterval` 計算基盤へ直接dispatchする。expression→interval変換，decimalization，guard-digit増加はFFTと共通化した。

`matrixRank`はexact入力ではexact eliminationを優先する。近似計算基盤ではmachine epsilonを用いず，intervalが0を含むがexact zeroでもないpivotは`PrecisionInsufficient`としてguard precisionを増やす。full rank等を非零pivotから証明できる場合は確定するが，近似値だけからrank deficiencyを推測しない。

2026-08-13のRelease / LTO off計測例:

| size | exact `dot` | exact `det` | exact `rref` | `N[det,16]` |
| ---: | ----------: | ----------: | -----------: | ----------: |
|    8 |           — |    0.283 ms |     0.434 ms |    1.410 ms |
|   12 |           — |    1.290 ms |     2.233 ms |   10.144 ms |
|   16 |    0.328 ms |    3.735 ms |     5.799 ms |   22.799 ms |
|   32 |    1.702 ms |           — |            — |           — |
|   64 |   16.553 ms |           — |            — |           — |

# 15.7. Bareiss / fraction-free exact Matrix — 採用

ではexact実数行列を行ごとの分母LCMで整数行列へliftし，`IntegerMatrixBuffer`上のBareiss eliminationを共通kernelとして追加した。整数行列はそのまま，Rational行列は各行を非零整数倍してから処理する。

- `det`: Bareissのfraction-free forward eliminationで計算し，Rational入力では行scale積を最後に一度だけ戻す。
- `rref`: Bareissでinteger echelon formまで進める。full-column-rankならRREFが単位列で確定するためRational backward phaseを省略し，rank-deficient caseだけcanonical Rational RREFを構築する。
- `matrixRank`: echelonのpivot数だけで決定し，RREF全体を構築しない。
- `inverse`: `B=D A` として `[B|D]` をfraction-free forward eliminationし，最終pivotを共通分母とするBigInt back-substitutionで右側を直接解く。generic Rational Gauss-Jordanへ戻さない。
- `solveLinear`: unique full-column-rank caseはinverseと同じBigInt back-substitutionを使い，最後にだけRationalを生成する。
- exact complex: `Q(i)`等へ整数liftする専用環をまだ持たないため，従来`Number` Gaussian/Gauss-Jordanを代替経路として保持する。

pivotは数値安定性のためではなく中間BigInt growthを抑えるため，候補中でbit lengthが小さい非零値を優先する。Bareissの各除算は`BigInt::divmod`で余り0を検証し，fraction-free invariantが壊れた場合は黙ってtruncationしない。

2026-08-13 Release / LTO off，同一benchmark入力でGaussianと比較:

| size | `det` Gaussian | `det` Bareiss | speedup | `rref` Gauss-Jordan | `rref` Bareiss | speedup |
| ---: | -------------: | ------------: | ------: | ------------------: | -------------: | ------: |
|    8 |       0.256 ms |      0.047 ms |    5.4x |            0.447 ms |       0.058 ms |    7.7x |
|   12 |       1.241 ms |      0.109 ms |   11.4x |            2.058 ms |       0.146 ms |   14.1x |
|   16 |       3.570 ms |      0.385 ms |    9.3x |            5.798 ms |       0.375 ms |   15.5x |

`N[det[...],p]` / `N[inverse[...],p]`等はこのexact Bareiss結果を先に作らず，FFT共通のprecision-aware certified Matrix 計算基盤へ直接dispatchする。したがってBareiss採用はexact pathの改善であり，`N`の近似経路を後退させない。

2026-08-25のperformance-cliff再監査では，inverseの重さはBareiss forward kernelより後段のRational RREFに集中していた。そこでfull-rank square/unique solveでは最終pivot `D` を共通分母とし，後退代入を

```text
n_i = (b_i D - sum_{j>i} U_ij n_j) / U_ii
```

というexact BigInt divisionだけで行う経路へ変更した。`rref`もfull-column-rankならforward elimination終了時点で結果を直接構築する。rank-deficient `rref/nullSpace`は従来のcanonical Rational backward phaseを維持する。

同じpublic-path benchmarkでは16-bit 32×32で`inverse`約103 ms，`rref`約5.2 ms，`matrixRank`約5.2 ms，nullity 1の32×33 `nullSpace`約5.8 ms。96-bit 32×32ではそれぞれ約0.70 s / 34.7 ms / 34.1 ms / 37.6 ms，256-bit 32×32では約3.46 s / 147 ms / 145 ms / 144 msだった。48×49・256-bitのrank-deficient `rref/rank/nullSpace`も約1.25 / 1.13 / 1.10 sで，この範囲では突然のalgorithmic cliffではなく係数bit growthに沿った増加である。inverseの高bit側は最終1024個の巨大canonical Rational生成自体が支配し始めるため，共通分母をpersistent Array storageで共有する等は将来のrepresentation課題として分離する。

# 15.7.1. Modular / CRT exact Matrix — 採用

post-v1.5.3では，大きいdense整数/Rational行列でBareiss中間BigIntのbit growthを避けるため，31-bit prime field上のmodular 計算基盤を追加した。31-bit primeを使うことで積は`uint64_t`へ安全に収め，MSVC固有の`__int128`へ依存しない。Rational入力は従来どおり行ごとの分母除去で整数workspaceへliftしてから処理する。

- `det`: 各prime上でGaussian eliminationし，整数演算だけで得るHadamard上界を満たすまでincremental CRTを進め，centered representativeを一意復元する。
- `solveLinear`: 有限体解をCRTし，rational reconstructionした候補を元の整数系`A X = B`でexact verificationする。rank drop等を起こすbad primeはskipし，保証域までに復元できなければBareissへ切り替える。
- `inverse`: exact `det(A)`を共通分母として有限体`A^-1`から`adj(A)`像を作り，CRTした整数adjugateを`A adj(A)=det(A)I`でexact verificationする。計算基盤は実装済みだが，現測定範囲ではBareissが速いためautomatic pathへは送らない。
- `rref` / `matrixRank` / `nullSpace`: rank certificate設計を別問題として扱い，現段階ではBareissのままとする。

automatic dispatcherは係数部density 25%以上を前提とし，GCC Release測定から次の保守的policyを採用した。heightはworkspace中の最大係数bit長である。

| operation     | modularを選ぶ条件                                                                                       |
| ------------- | ------------------------------------------------------------------------------------------------------- |
| `det`         | order>=48，またはorder>=32 & height>=64，order>=24 & height>=192                                        |
| `solveLinear` | variables>=24，またはvariables>=12 & height>=96，variables>=8 & height>=256，variables>=6 & height>=512 |
| `inverse`     | automaticでは選ばない                                                                                   |

2026-08-20 GCC Release / LTO-off，3 iteration平均の代表値。単位はmsであり，絶対性能保証ではなくcrossover policyの資料である。

| workload                     | Bareiss | modular |
| ---------------------------- | ------: | ------: |
| `det` 32×32, 96-bit          |  37.691 |  26.270 |
| `det` 24×24, 256-bit         |  41.221 |  29.889 |
| `det` 32×32, 256-bit         | 174.095 |  58.399 |
| `det` 20×20, 512-bit         |  55.836 |  56.766 |
| `det` 24×24, 512-bit         | 136.994 |  70.267 |
| `det` 32×32, 512-bit         | 497.769 | 143.100 |
| `solveLinear` 12×12, 96-bit  |   0.453 |   0.234 |
| `solveLinear` 8×8, 256-bit   |   0.297 |   0.179 |
| `solveLinear` 24×24, 256-bit |  47.954 |   0.686 |
| `solveLinear` 6×6, 512-bit   |   0.285 |   0.113 |
| `solveLinear` 24×24, 512-bit | 156.020 |   0.703 |
| `solveLinear` 32×32, 512-bit | 528.836 |   1.586 |

20×20・512-bit determinantはこの測定でBareissが僅かに優位だったため，自動dispatchはここでmodularへ切り替えない。24×24・512-bitではmodularが明確に優位であり，既存の24次/192-bit条件に包含される。

`inverse`は2026-08-25に再測定してもcrossoverがなく，32×32で16-bit `15.2 / 84.5`，96-bit `74.2 / 657.6`，256-bit `357 / 3339` ms（Bareiss forward core / modular inverse）だった。したがって計算基盤の存在とautomatic採用を分離し続ける。なおpublic `inverse`全体では最終Rational materializationが別途支配し，256-bit 32×32で約3.46 sとなる。

閾値は数学的意味論ではなく性能policyである。compiler，BigInt実装，prime kernel，CPUが変われば`mmCal.Benchmarks --exact-linear-algebra 1`で再測定する。`N[det[...],p]` / `N[solveLinear[...],p]`等のcertified pathはexact reconstructionを経由せず，従来どおりprecision-aware interval 計算基盤へ直接dispatchする。

# 15.8. LU / fraction-free exact QR / certified Householder QR — 採用

分解処理は`linear_algebra/decomposition.*`へまとめ，`luDecomposition[A]`はrow-pivoted `P A = L U`を維持する。QRはexactとapproximateで算法を分離した。`N[qrDecomposition[A],p]`は従来どおりcertified Householderを直接使い，exact factorを先に展開しない。一方exact実数行列は2026-08-25にExpr-level Householderからfraction-free直交化へ置換した。

## exact QR: 正規化を最後まで遅延 — 採用

旧exact Householderでは，各反復で`norm -> sqrt -> reflector -> Expr arithmetic -> simplify`を繰り返した。このため最初のradicalが次列のnormへ入り，2×2約1.2 ms，3×3約59 msに対し4×4では約18秒，formatted output約677 KBまで膨張した。旧`maximumExactQrOrder = 3`はこの式爆発を防ぐpolicyだった。

新経路では，Rational列をまずprimitive整数vectorへliftし，直交化中は平方根もRational除算も作らない。直交basis `p_i` とその整数norm `d_i=p_i^T p_i`を保持し，射影除去は

```text
p <- d_i v - (p_i^T v) p_i
```

をGCDで約分して進める。最終materializationでのみ

```text
Q[:,i] = p_i / sqrt(d_i)
R[i,j] = (p_i^T a_j) / sqrt(d_i)
```

をExpr化する。平方根分解も列ごとに一度だけ行い，同じradical nodeをQ/R全要素で共有する。これにより旧`maximumExactQrOrder`は撤去した。上三角／上台形の`{I,A}` fast pathはそのまま残す。

## Gram + symmetric Bareiss（fraction-free LDL^T相当） — 採用

full-rank caseではprimitive列`C`からGram行列`G=C^T C`を作り，対称Bareiss消去で主値 determinant列とlower係数を得る。これから平方根を作らず直交整数basisを復元する。この経路はdirect fraction-free Gram-Schmidtより代表8～64次で概ね3～5倍速かった。先頭主値 minorが0になるcaseやrank-deficient caseではdirect fraction-free直交化へ代替経路し，標準basisを同じkernelで直交補完する。

通常のRational LDL^Tも比較したが，Rational正規化が支配しdirect fraction-freeより約6～35倍遅かったため不採用とした。`A^T A`を作ることによる数値条件数悪化はexact arithmeticでは問題にならず，ここでの判断基準はBigInt growthと実測時間である。

一般巨大normについて平方因子試行を省略し，`1/sqrt(n)=sqrt(n)/n`の未簡約radicalをそのまま返す案も実測した。しかし32×32・256-bit等で有意な改善がなく，caseによっては僅かに退行し，canonical radicalも弱めるため棄却した。現在は列単位にradical decompositionを一度だけ行う。

2026-08-25 GCC Release / LTO-off，`--exact-linear-algebra 1`のpublic exact path代表値。係数heightを含めて測っており，最後のRational/radical Expr materializationも時間に含む。

| coefficient |  size | exact QR |
| ----------: | ----: | -------: |
|      16 bit |   8×8 |  2.22 ms |
|      16 bit | 16×16 |  14.3 ms |
|      16 bit | 24×24 |  54.0 ms |
|      16 bit | 32×32 |   147 ms |
|      96 bit | 16×16 |  96.4 ms |
|      96 bit | 24×24 |   418 ms |
|      96 bit | 32×32 |   1.30 s |
|     256 bit | 16×16 |   415 ms |
|     256 bit | 24×24 |   2.00 s |
|     256 bit | 32×32 |   6.85 s |

低～中bitでは旧4×4 expression cliffは消え，次数に対して滑らかに伸びる。高bit・高次では直交化kernelより最終的な巨大Rational/radical出力そのものが支配する。この領域はhard capへ戻さず，EvaluationBudgetと自然な出力costに任せる。将来さらに詰めるならshared-denominator / delayed-radicalの内部表現を別途設計すべきであり，exact semanticsを弱めるmachine近似へは逃がさない。

## certified Householder — 維持

Householderのapproximate kernelでは複数列を一度のrow-major走査で処理するcolumn-block版も試作した。block=1/8/16/32を複数回Release計測したが，8～24次で差は概ね数%以内かつ最速blockが安定せず，BigFloat/interval演算costが支配的だった。このため既定はblock=1相当とし，block kernelとbenchmarkのみ残した。

2026-08-13 Release / LTO offの代表値:

| size | exact `LU` | `N[LU,16]` | `N[QR,16]` |
| ---: | ---------: | ---------: | ---------: |
|    8 |   0.290 ms |   1.811 ms |   7.062 ms |
|   12 |   1.315 ms |   5.355 ms |  21.218 ms |
|   16 |   3.548 ms |  10.562 ms |  46.718 ms |

# 15.9. reduced SVD — 採用

数値SVDは条件数を二乗する`A^H A`を形成せず，Householder bidiagonalizationの後にone-sided Jacobiで列を直交化する。実数・複素数で同じprecision-aware policyを使い，候補factorはreconstructionとU/V orthogonalityを区間監査してから返す。exact SVDは自然に閉じるcaseへ限定する。

2026-08-13 Release / LTO offで`N[svd[A],16]`を複数回計測した代表値:

|  size |      time |
| ----: | --------: |
|   4×4 |  約3.8 ms |
|   8×8 | 約17.1 ms |
| 12×12 | 約46.6 ms |
| 16×16 | 約82.3 ms |

この範囲では急激な悪化はなく，概ね三次成長に沿う。支配costはJacobi反復とBigFloat演算であり，bidiagonalization側だけをcache block化しても寄与率が小さいため，SVD専用block policyは現時点で追加しない。QRのblock kernelは独立benchmarkとして残し，将来計算基盤が変わった時に再測定する。

# 15.10. Eigen / complex Schur — 採用

一般固有値問題は`A^H A`等へ変形せず，Complex BigFloat上でHessenberg reduction → implicit shifted QR → complex Schur形へ進む。固有vectorが必要な場合はSchur三角行列からback substitutionし，Schur vectorを掛け戻す。exact pathは上三角/対角とdistinct-root exact Number 2×2を優先する。

反復停止精度と最終certificateを同じ桁へ置くと行列積で誤差余裕を使い切ることがfixed-seed 3×3で判明したため，内部QR停止精度は表示要求より10桁相当厳しくする。元入力の`ComplexInterval`に対する`A Q-Q T`および`A v-λv`を区間監査し，Schur vectorのunitarityも同時に監視する。非正規行列では固有値・固有vectorのcomponentwise enclosureを無条件には主張せず，Schur/eigenpair relationをcertification境界とする。

2026-08-13 Release / LTO off，random decimal Matrix（[-1,1]，小数10桁相当）:

| size | `N[eigenvalues,16]` | `N[eigensystem,16]` |
| ---: | ------------------: | ------------------: |
|    4 |             13.4 ms |             14.0 ms |
|    8 |             68.1 ms |             77.0 ms |
|   16 |            364.8 ms |            466.3 ms |
|   32 |             2826 ms |             3658 ms |
|   64 |            19436 ms |   >35 s（計測上限） |

# 15.11. large dense Matrix監査

添付のrandom matrix generatorと同じ「[-1,1]，小数10桁」という入力特性を再現する`--matrix-large <op> <size> [digits]`を`mmCal.Benchmarks`へ追加した。算法benchmarkではparser costを除くため，同じ範囲・10桁量子化をexact RationalとしてC++から直接構築する。PythonとRNG列そのものは一致させず，数値分布と桁幅を合わせる。一方CLI負荷はPython generatorと同じ出力形式を使って別測定する。

32/64次のrepresentative timing:

| op                  |  32×32 |   64×64 |
| ------------------- | -----: | ------: |
| `N[dot,16]`         | 167 ms |  1.21 s |
| `N[det,16]`         | 294 ms |  2.96 s |
| `N[inverse,16]`     | 1.22 s |  8.68 s |
| `N[matrixRank,16]`  | 400 ms |  4.16 s |
| `N[solveLinear,16]` | 549 ms |  5.42 s |
| `N[nullSpace,16]`   | 401 ms |  4.01 s |
| `N[LU,16]`          | 147 ms |  1.41 s |
| `N[QR,16]`          | 1.47 s | 12.88 s |
| `N[SVD,16]`         | 1.40 s | 11.60 s |
| `N[eigenvalues,16]` | 2.83 s | 19.44 s |

1024×1024では算法より先にrepresentation costが目立つ。C++から直接1,048,576個の10桁Rational Exprを構築したbenchmark processは入力だけで最大RSS約0.69 GB。`transpose`本体は約107 ms，`trace`本体は約60 msだった。Python形式の約14.16 MBテキストをCLIへ渡し`dimensions[...]`だけを評価した測定ではwall約10.9 s，最大RSS約1.99 GBだった。さらに1024次`N[dot,16]`は10秒上限で未完了（最大RSS約0.96 GB），`N[LU,16]`も10秒上限で未完了（最大RSS約1.59 GB）だったため，QR/SVD/Eigenの1024実走はメモリ圧迫を避けて中止した。

v1.5.3の`Expr::Node` typed-node refactor後，同一x86-64 GCC / Release / LTO offで旧variant sourceと新sourceを同じ`--matrix-large transpose 1024 16`へ掛けて再比較した。旧variant版は最大RSS `693312 KiB`（約677.1 MiB），typed-node版は`299668 KiB`（約292.6 MiB）で，約384.4 MiB / **56.8%削減**。単発`transpose` timingは181.7 ms→157.4 msだったが，timingはnoiseを含むため採用根拠はRSS削減と全regression維持を主とする。

第二段階ではpersistent `ArrayExpr`をfixed-size immutable pageへpackedし，shape / offset / stridesをbackingから分離した。最初に試した単一`vector<Rational>`方式は保存時のRSSは下がるものの，transposeでRational/BigIntを100万要素deep copyし，1024×1024 direct-packed transposeが約650～675 msへ退行したため棄却した。採用版は1024要素pageを`shared_ptr<const page>`で共有し，transposeをstride交換だけのview生成にした。

`ArrayBuilder`はpromotionを現在page内だけへ限定する。完成済みpageは不変なので，大規模Arrayの末尾にsymbolic値が出ても全量Generic化しない。1,048,576要素の最後だけ`x`にした専用測定はall-integer版とほぼ同じ約0.30 s / 約69.8 MiBで，最後のpageだけGenericだった。また矩形brace Lowererはleafを単一builderへ直接流し，numeric literalを一度`Expr` nodeへ包んでからpackする二重表現を避けた。

現`--matrix-large transpose 1024 16`はbenchmark fixtureも`ArrayBuilder`からexact Rationalを直接構築し，最大RSS `136576 KiB`（約133.4 MiB），transpose本体約0.059 ms。`trace 1024`は最大RSS約133.5 MiB，trace本体約10.95 msだった。P1 typed-nodeの292.6 MiBからさらに減っているが，fixture construction pathも同時に現実装へ合わせているため，この差はpersistent storage + builderの総合改善であって単一変更A/Bではない。

CLI側では13.63 MBの1024×1024・10桁decimal literalを`dimensions[...]`へ入力した測定がwall約4.55 s，最大RSS `441564 KiB`（約431 MiB）だった。v1.5.2監査の約14.16 MB / 10.9 s / 1.99 GBとは入力textが完全同一ではないため厳密比較ではないが，parse後段のExpr allocation削減が大きく効いていることを示す。

64次値から純粋なO(n^3)を仮定した1024次の粗い外挿でも，`N[LU]`約1.6時間，`N[dot]`約1.4時間，`N[det]`約3.4時間，`N[solveLinear]`約6.2時間，`N[inverse]`約9.9時間，`N[SVD]`約13時間，`N[QR]`約15時間，`N[eigenvalues]`約22時間となる。32→64の実測指数をそのまま延長すると約1～21時間程度へ揺れるため，これらは予測値であって1024実測ではない。cache・allocator・guard precision・反復回数により悪化し得る。

結論として，1024 dense自体はmachine double + BLASの世界では特別巨大な次数ではないが，mmCalのcertified arbitrary-precision dense算法にとっては依然stress領域である。一方，persistent exact Arrayのrepresentation固定費はtyped-node + paged packed backing + direct builderで大きく下がった。approximate SVD/Eigen等には既に連続working bufferがあるため，次はstorage改善後のcost balanceでblock化・threadingを再評価する。

# 15.12. persistent algebraic field

## compositum / embedding再利用 — 採用

primitive-element compositum cacheでは，同じRoot pairから一度証明したprimitive-element compositumと両operandのpower-basis embeddingをbounded cacheへ保持する。Real `AlgebraicElement`からcanonical `root[minpoly,k]`へ戻す際も，chosen generator interval上でexact Rational interval評価し，Sturm root countで該当root indexを直接証明する。全根isolationの重複と，不可能なRational / `Q+iQ`退化probeを避ける。

代表式：

```text
(root[{-2,0,1},2]+root[{-3,0,0,1},1])
*(root[{-2,0,1},2]-root[{-3,0,0,1},1])
```

旧実装では約1.23 s/回だったが，persistent compositum cache導入後は同一GCC Release / LTO off環境で初回約0.14～0.15 s，warm約0.02 sまで短縮した。

## reciprocal reuse — 採用

reciprocal cacheでは`Q[t]/(m)`上のextended Euclidで求めたexact reciprocalを，`NumberFieldContext`ごと最大16組のthread-safe LRUへ保持する。`inverse(inverse(x))=x`なので1 entryを双方向pairとして扱う。有理定数座標は`Q`の埋め込みを使い直接逆数化する。

代表測定：

| degree | 処理       |  cache前 | reciprocal cache warm |
| -----: | ---------- | -------: | -------------: |
|      6 | reciprocal | 約178 us |      約0.24 us |
|      6 | divide     | 約243 us |        約84 us |
|     12 | reciprocal | 約525 us |      約0.44 us |
|     12 | divide     | 約773 us |  約211～233 us |

### persistent multiplication matrix cache — 棄却

12次体で左乗算matrixを試したところ，通常乗算約112 usに対してmatrix-vector約101 usで，matrix構築自体が約228 usだった。1回あたりの短縮が約10%に留まり，同じelementによる多数回乗算がなければ償却できない。全`AlgebraicElement`へmatrix cacheを持たせる複雑性とmemory retentionに見合わないため現段階では採用しない。

## incremental Krylov minimal polynomial — 採用

旧`AlgebraicElement::minimalPolynomial()`は，`1,a,...,a^k`について各`k`ごとに新しいRational matrixを構築し，Gauss-Jordanを最初からやり直していた。次数`d`まで進むと，同じ独立性情報を繰り返し計算する。

incremental Krylov minimal-polynomial導出では`1,a,a^2,...`を1列ずつ追加し，既存のexact row-echelon stateで新列だけをreduceする。最初に線形従属した

```text
c0 + c1 a + ... + a^k = 0
```

は，それ以前の`1,a,...,a^(k-1)`が独立なのでそのままminimal polynomialである。primitive-element tensor algebraでも同じincremental basisを使い，`theta`のminimal polynomialだけでなく，確立済みpower basisへの`alpha` / `beta`座標変換にも消去状態を再利用する。

同一GCC Release / LTO offの代表値：

| field degree | 旧repeated Gauss-Jordan | incremental Krylov first |
| -----------: | ----------------------: | -----------------------: |
|            6 |                約651 us |                 約475 us |
|           12 |               約18.9 ms |            約7.7～7.9 ms |

さらに導出済みminimal polynomialをexact power-basis座標keyでfieldごと最大16 entryのthread-safe LRUへ保持する。12次のwarm hitは約0.59 usである。cache miss/evictionは再計算を増やすだけで，数学的結果は変えない。

### multiplication-matrix minpoly / modular reconstruction — 保留

FLINT等にはexact Rational matrixのminimal-polynomial 計算基盤があり，有限次元代数をmultiplication matrixとして扱う方法自体は標準的である。しかし現在のmmCalはalgebraic-field候補次数をboundedに保ち，power-basis座標を既に持つ。現状の次数域ではincremental Krylovが小さな実装で十分な改善を出したため，常時multiplication matrixを作ってmatrix-minpolyへ渡す経路は追加しない。modular image + rational reconstructionやfraction-free matrix minpolyも，将来次数・係数heightが増えてfirst derivationが再び支配的になった時の候補とする。

# 15.13. black-box test workload監査

`test_set/tester.py`へ`--timings [N]`を追加した。各test fileについてmmCal processのwall timeを計測し，遅い順に表示する。test記述のfile I/O / parseは従来どおり総elapsed前に完了する。

2026-08-16のincremental Krylov実装 / GCC Release / LTO offで全1652 black-boxを測定した代表値：

| test file                           | tests | wall time |
| ----------------------------------- | ----: | --------: |
| `test16_exact_calculus_solver.txt`  |    85 | 約2927 ms |
| `test5_matrix.txt`                  |   105 | 約1457 ms |
| `test9_special_func.txt`            |    75 |  約297 ms |
| `test8_calculus.txt`                |    46 |  約173 ms |
| `test22_number_field_interning.txt` |     1 |  約152 ms |

`test16`はintegration / high-degree Solve / algebraic Root constructionが主なstress集合であり，`test5_matrix`は名称に反して小Matrix演算よりexact FFT/DFTのround-tripが大きな比率を占めていた。非2冪はcyclotomic quotient，16点以上の2冪inverseは円分体再埋込みで主要な式爆発を解消した。今後は128点超の2冪と，より大きい非2冪長でcrossoverを再測定する。`test9`では`N[ibeta[1/3,2/3,1/4],20]`と`N[gamma[1/3],20]`が相対的に重い。

このtimingはtest correctnessのPASS/FAIL判定には使わず，optimization対象を選ぶためのprofiling signalとしてのみ利用する。

# 15.14. certified `gamma` / `ibeta`

`tester.py --timings`で`test9_special_func.txt`を分解すると，`N[ibeta[1/3,2/3,1/4],20]`と`N[gamma[1/3],20]`が明確なhotspotだったため，両計算基盤を直接計測して最適化した。

## `ibeta` point/shared normalization — 採用

旧`encloseIncompleteBetaRegularized`はexact pointでも`lower==upper`を別々に評価し，各endpointで`Beta(a,b)`まで再構築していた。2F1自体は20桁級で約1.9 msなのに対しBeta/Gamma normalizationが支配的だったため，次へ変更した。

- `lower==upper`ならpointを1回だけ評価
- interval endpoint間で`Beta(a,b)`を1回だけ構築して共有
- complement `I_x(a,b)=1-I_{1-x}(b,a)`でも`B(a,b)=B(b,a)`を利用して同じnormalizationを共有
- pointが`x<=1/2`なら必要guardを`+40 bit`に抑え，complementを使うpoint / intervalだけ`+80 bit`
- `x=0,1`はnormalizationを構築せずexactに返す

incremental Krylov実装のdirect warm probeとの代表比較：

|           precision | 旧 `ibeta[1/3,2/3,1/4]` |       修正後 |
| ------------------: | ----------------------: | -----------: |
|              80 bit |                約116 ms |   約37–43 ms |
|             160 bit |                約255 ms |  約98–117 ms |
|             320 bit |                約835 ms | 約387–392 ms |
|       640 bit first |                約10.3 s |      約5.3 s |
| 640 bit warm repeat |                約10 s級 |      約1.0 s |

warm repeatの追加短縮は後述のStirling-plan cacheも受ける。算法・branch contractは変更していない。

## Gamma lazy Bernoulli / Horner / plan reuse — 採用

旧Gamma 計算基盤ではAkiyama–Tanigawa法で`B0...B128`をfirst use時に一括生成していた。また低精度でも最初の小さいshift候補で`k=1...64`を総当たりするため，最終planが小さい`k`でも高次Bernoulliまで生成しやすかった。

- Akiyama–Tanigawaの内部状態を保持し，要求された`B_2k`までだけ逐次延長するthread-safe lazy cacheへ変更
- 低～中精度のStirling探索上限を`min(64,max(16,ceil(bits/5)))`として，不要な高次Bernoulli生成を避ける
- Stirling和を`x^-1(c1+x^-2(c2+...))`のHorner形へ変更
- exact pointのrecurrence productをbalanced exact Rational productとして構築してからlogを1回だけ取る
- `logGamma(1)=logGamma(2)=0`，`Gamma(1)=Gamma(2)=1`をcertified 計算基盤でも即時処理
- `(inputLower,precisionBits)`から得たexact Stirling planをthread-local最大16 entryで再利用

別process cold-startの代表値：

|           precision | 旧実装 |     修正後 |
| ------------------: | --------: | ---------: |
| 80 bit `gamma[1/3]` |   約79 ms |   約9.5 ms |
|             160 bit |   約89 ms |    約43 ms |
|             320 bit |  約256 ms | 約260 ms級 |
|             640 bit |  約1.42 s | 約1.42 s級 |

したがってlazy化は主に低～中精度のfirst-use taxを除去する。高精度first callの主要costは依然としてStirling/recurrence本体に残る。一方plan cacheが効く同一threadの反復では640 bit `lgamma[1/3]`が約1.36 s級から約0.34 sまで低下した。

### Bernoulliを`B256`までeager生成 — 実測棄却

Stirling項数を増やす前提として，旧Akiyama–Tanigawa実装を単純に`B256`まで一括拡張する案も検討した。しかしexact Rational生成だけで代表測定は`B128`約68 msに対し`B256`約560 msまで増え，低精度を含む全first useへ大きな固定費を課す。必要次数まで状態を逐次延長するlazy cacheの方が明確に有利なので，eager拡張は採用しない。

### `maximumK > 64`によるshift削減 — 実測棄却

分析段階では640 bit級で旧planが`shift=272, omittedK=64`へ達するため，Bernoulliをlazyに拡張して`K=80..128`を許せばshiftを減らせると予想した。実際に`K=96`まで拡張するとshiftは縮んだが，代表測定は`gamma[1/3]`約1.37 sから約2.7 sへ悪化した。

原因は，高次Bernoulli係数のexact Rational生成・interval conversion・より長いStirling和のcostが，recurrence shift削減を上回ったためである。よって現計算基盤では`K<=64`を維持する。将来rectangular splitting，binary splitting，より効率的なBernoulli 計算基盤を導入した時だけ再評価する。

### fixed-k exact binary-search plan — 実測棄却

高精度で`k=64`を固定し，remainder boundをshiftに対して二分探索する案も試した。しかし巨大Rationalの`x^(2k-1)`を各probeで構築するcostが大きく，linear scanより悪化したため採用しない。現状はlow/mid precisionの`K` budget削減とplan cacheの方が効果が大きい。

旧実装は`certified_special_functions.cpp`内にコメントアウトで残し，置換理由を隣接記述している。通常はdead codeを残さない方針だが，今回は算法比較と将来の再評価用に意図的に保存した。

注：この節のstateful lazy Bernoulli generatorは旧の実装であり，現在では`B_2...B_128`のstatic exact tableへ置換された。旧generator自体は比較用コメントとして残している。

# 15.15. exact Rational `Gamma` / high-precision Stirling planner

640 bit以上の`gamma[1/3]` / `ibeta[1/3,2/3,1/4]`を再計測したところ，算法以前に入力表現とplannerに大きな無駄が残っていた。Johansson, _Arbitrary-precision computation of the gamma function_ (arXiv:2109.08392)のrational rising-factorial / Stirling パラメータ-selectionの整理も参照している。

## exact Rational identityをcertified 計算基盤まで保持 — 採用

Evaluatorでは`1/3`をexact Rationalとして保持しているが，旧経路はspecial-function 計算基盤へ入る前に`RealInterval::fromRational`へ変換していた。非dyadic Rationalはpoint intervalにならないため，Sexact rising-factorial経路が`gamma[1/3]`では実質使われていなかった。

- `Gamma` / `LogGamma` / `Beta` / `BetaLog`でexact Rationalをinterval化より先に検出する
- positive Rational `Gamma`ではoriginal `p/q`をStirling shiftまで保持する
- `(p/q)_n`を`prod(p+qk)/q^n`としてbalanced binary productし，最後に一度だけRational化する
- `Beta(a,b)`では`a,b,a+b`をexact Rationalのまま3つの`LogGamma`へ渡す
- negative exact Rationalはreflectionで`1-x`をexact Rationalのまま保持し，`sin(Pi x)=sinTurns(x/2)`としてexact turn reductionを使う

これにより「non-dyadic Rationalをintervalへ落としたためexact用高速路が死ぬ」という表現境界の損失を除去した。

## Bernoulli `B_2...B_128` static exact table — 採用

stateful lazy Akiyama–Tanigawaは低精度のeager taxを避けたが，1000 bit級で初めて高いBernoulli次数へ到達した際，内部Rational stateを逐次更新するcold-startが再び大きくなった。Stirlingが現状使用する`B_2...B_128`は固定された厳密有理定数なので，numerator/denominatorのstatic decimal tableへ置換し，参照された値だけ`BigInt/Rational`へlazy parseする。旧stateful generatorは変更理由付きコメントとして隣接保存している。

これは`B256`までruntime eager生成する案とは別である。高次Bernoulliを動的に大量生成する案は引き続き棄却し，現行`K<=64`に必要な既知定数だけをtable化した。

## high-precision planをBigInt inequalityで直接証明 — 採用

旧plannerはshiftを8ずつ増やし，各shiftで`k=1...64`のremainder boundをexact Rationalで評価していた。1000 bit級ではplanner自体が大きなcold costになる。

positive Rational `x=(p+qs)/q`，Stirling coefficient `c=A/B`，`d=2k-1`に対し，

```text
|c| / x^d <= 2^-P
```

は正の整数だけを使って

```text
|A| q^d 2^P <= B (p+q s)^d
```

と同値である。768 bit超では`k=64`についてこの不等式をBigInt cross multiplicationで評価し，doubling + binary searchで十分なshiftを直接求める。最後のremainder boundもexact Rationalで再構築するため，plannerの高速化にfloating-point heuristicを使わずcertified contractを維持する。

過去に棄却した「fixed-k exact Rational二分探索」は，各probeでnormalized Rational `x^(2k-1)`を構築した実装である。今回採用したのはGCD/normalizationを挟まずBigInt cross productだけを使う別実装であり，旧棄却理由と矛盾しない。

## 代表測定

GCC Release / LTO off，同一環境の代表値：

| workload                       |                      旧 |                               新 |
| ------------------------------ | ----------------------: | -------------------------------: |
| `gamma[1/3]`, 640 bit          |              約0.69 s級 | first 約65 ms / warm平均 約38 ms |
| `gamma[1/3]`, 1280 bit first   |                 約1.5 s |                         約0.11 s |
| `ibeta[1/3,2/3,1/4]`, 640 bit  |                 約5.3 s |                         約0.21 s |
| `ibeta[1/3,2/3,1/4]`, 1280 bit |                 約4.5 s |                         約0.49 s |
| `gamma[-1/3]`, 1280 bit        | interval reflection経路 |                         約0.15 s |

通常の`--special-functions 3`では80/160/320/640 bitの`gamma[1/3]`が約4.5/7.2/13.7/37.6 ms，対応する`ibeta`が約15.8/32.2/44.3/210 msだった。1280 bitも恒久benchmarkへ追加する。

## 改良Stirling主和 / Algorithm 6 — 次候補

JohanssonのTheorem 3.5 / Algorithm 6は，低indexのBernoulli項とhigh-index hypergeometric tailへStirling主和を分割し，高精度で必要なBernoulli数を減らしながら高速化する。FLINT/ArbのGamma 計算基盤もimproved Stirling sumでrectangular splittingとhigh-index re-expansionを使う。1280 bit級の主要な表現/planner overheadは取れたため，次に1000 bit超～さらに高精度を伸ばす場合は単純な`K>64`や`B256` runtime生成へ戻らず，この方向を実装候補とする。

# 16. `mmCal.Benchmarks`

v1.5.1でVisual Studio solutionへ独立Console projectとして追加した。

```text
mmCal
mmCal.Core
mmCal.Tests
mmCal.Benchmarks
```

通常testへbenchmark時間を混ぜず，CoreへProjectReferenceして次を行う。

- fixed-seed random BigInt division invariant
- decimal round-trip
- certified `exp(x)exp(-x)` invariant
- certified `log(x)+log(1/x)` invariant
- multiply / square / divide 閾値 benchmark
- factorial benchmark
- decimal parse/toString benchmark
- high-precision `Pi/exp/log` benchmark
- certified `gamma/ibeta` precision-scaling benchmark (`--special-functions`)
- exact/certified FFT benchmark + direct/Bluestein crossover
- fixed-seed certified Matrix invariant（Bareiss / LU / QR / solve / nullSpace / 実・複素SVD / Eigen）
- fixed-seed certified FFT round-trip invariant

実行例:

```text
mmCal.Benchmarks
mmCal.Benchmarks --full
mmCal.Benchmarks --random-only
mmCal.Benchmarks --benchmark-only
mmCal.Benchmarks --matrix-large nsvd 64 16
mmCal.Benchmarks --special-functions 1
```

閾値変更時は速度だけでなくrandom invariantを先に通す。

---

# 17. v1.5.1で意図的に採用しなかった一覧

| 候補                                                   | 判断             | 理由                                                                                                                          |
| ------------------------------------------------------ | ---------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| Prime-Swing factorial                                  | 棄却             | 現product treeより巨大factorialで遅い                                                                                         |
| binary GCD                                             | 棄却             | Euclidean+B/Zより4～30倍遅いcase                                                                                              |
| Karatsuba vector pool                                  | 棄却             | 5～10%程度退行                                                                                                                |
| Karatsuba depth scratch                                | 棄却             | 同様に管理costが勝る                                                                                                          |
| Toom-3 square                                          | 棄却             | Karatsuba squareより遅い                                                                                                      |
| 低閾値 Toom-3                                     | 棄却             | 512～1024 limbsでoverheadが勝つ                                                                                               |
| machine `fmod`による巨大trig縮約                       | 棄却             | certified semanticsを失う                                                                                                     |
| 全体をMachine/double化                                 | 方針として不採用 | exact-firstの意味論を変える                                                                                                   |
| persistent algebraic multiplication-matrix cache       | 棄却             | 12次体で構築約228 usに対し乗算は約112→101 usに留まり，償却条件が厳しい                                                        |
| multiplication-matrix minpoly / modular reconstruction | 保留             | 現在のbounded degreeではincremental Krylovが小さな実装で十分な改善を出す                                                      |
| Gamma Bernoulli `B256` runtime eager生成               | 棄却             | exact Rational生成だけで約560 msのfirst-use tax。`B_2...B_128`固定表へ移行したが，より高次をruntime大量生成する案は採用しない |
| Gamma Stirling `maximumK>64`                           | 棄却             | 640 bit級でshiftは減るが，高次Bernoulli/Rationalと長いStirling和が勝ち約1.37 s→約2.7 sへ退行                                  |
| Gamma fixed-k **Rational-power**二分探索               | 棄却             | normalized Rational `x^(2k-1)` probeが高価。同じ判定をBigInt cross multiplicationへ再定式化した別方式を採用                   |

---

# 18. 次の候補

`Expr::Node` typed-node化はv1.5.3で採用済み。公開APIを維持した単独refactorとしてinternal regressionとrandom fuzzerを通し，同一環境の1024 MatrixでRSS約56.8%削減を確認した。

次の優先度は次のように考える。

1. paged packed storage導入後のcost balanceでblocked LU / QR等を再測定
2. pure numeric working kernelだけを対象にthreading 閾値を検討
3. parser AST側の巨大brace temporaryを必要に応じて追加監査
4. BigUInt / BigInt SBOは保守性を損なわないstorage abstractionとして独立benchmarkし，複雑さに見合う場合だけ採用
5. Toom-4 / higher Toom crossover，さらに巨大な整数ではFFT/NTT multiplication
6. Lehmer GCD
7. `log`のbit-burst / AGM 計算基盤
8. exact Cyclotomic FFTのmixed-radix / prime-length高速化と，128点超の2冪専用backend拡張。既存quotient計算基盤のcrossoverを実測してから判断する

Arrayについては，単一flat packed vectorを棄却し，immutable paged backing + stride viewを採用した。approximate Matrix algorithmは既存の専用連続working bufferを維持し，persistent Array storageと無理に統合しない。BigUInt SBOは今回明示的に見送る。

採用時にはこの文書へ「なぜ採用したか」「なぜ前案を棄却したか」「どのbenchmarkで判断したか」を追記する。

# 19. 特殊函数 performance-cliff 横断監査

個別benchmarkだけでは「既知の重い入力」を測る方向へ偏るため，`mmCal.Benchmarks --performance-cliffs [iterations]`を追加した。公開`KernelSession`から特殊函数のargument / precision sweepを評価し，wall timeと`EvaluationUsage`を同時に記録する。

重要なのは，**計算量制限付き境界とperformance cliffを分離した**ことである。numericからheld/errorへ切り替わった点は`BOUNDARY`として表示し，成功点同士の時間またはcertified workが急増した場合だけ`CLIFF`候補とする。これにより，`1F1`の161や`polylog`の99/100のように意図的に即時停止する点を「遅い」と誤判定しない。

2026-08-25，GCC Release / LTO off，2回測定のwarm側では以下が目立った。

- `1F1[1/2,5/4,z]`: `z=1`約0.8 ms，64約80 ms，128約619 ms，160約1.08 s。当時の`|z|<=160`境界は停止性を守っていたが，境界内部に既に大きなcliffがある。
- `polylog[2,z]`: `z=1/2`約0.8 msから`49/50`約139 msまで増加。境界直前ほど級数収束が支配的になる。
- 楕円積分: `m=9/10`で`F`約162 ms，`E`約206 ms，`Pi[n=9/10]`約307 ms。`Pi`はcharacteristic側の収束悪化が特に大きい。
- complex Fresnel: `1+I`約38 ms，`4+I`約166 ms，`6+I`約522 ms，`7+I`約1.97 s。`|z|<8`境界直前の級数計算基盤は明確な最適化候補。
- real zeta一般点: `s=3`約19 ms・約1.0万 refinementに対し，非整数側では約4.2万 refinementへ上がり，`s=3/2`約155 ms，`s=12/5`約321 ms級。極接近そのものよりEuler–Maclaurin側の一般経路への移行が大きい。
- complex zeta: 20桁の代表点で約0.45–0.55 s。`CertifiedRefinement`は約0.8–1.7万であり，単純なrefinement countより複素interval arithmetic / Euler–Maclaurin kernel自体のcostが大きい。
- complex `2F1` 解析接続代表は約0.98 s，complex `li[-2+I]`は約1.23 sで，既存fast path後も重い。

逆にGamma，ibeta，通常の`2F1` precision sweepは10→80桁で急激な4倍cliffを示さず，今回の監査では優先度を下げられる。

この結果から，特殊函数の次の最適化候補は概ね次の順である。

1. complex Fresnelの`|z|`大側にasymptotic / 解析接続 計算基盤を導入できるか検討
2. `1F1`大正実数でplain series以外の計算基盤を検討
3. elliptic `F/E/Pi`のnear-boundary seriesをAGM/Carlson symmetric forms等へ置換できるか検討
4. real/complex zetaのEuler–Maclaurin workとinterval演算を分解profiling
5. complex `2F1` 解析接続とcomplex `li`の内部合成costを再分解

ただし新計算基盤採用はexact-first/certified contractを崩さず，既存boundary/metamorphic fuzzerを通したうえでrunnerの前後比較で判断する。

# 20. complex Fresnel / large positive 1F1 の exact-majorant cliff 除去

performance-cliff runnerで最優先になったcomplex Fresnelと大きな正実数`1F1`を分解すると，どちらも主要因は単純な級数項数ではなく，**保証を付けるためのexact `Rational` bookkeepingが毎項肥大化していたこと**だった。

complex FresnelではMaclaurin項そのものは固定precisionの`ComplexInterval`だった一方，Taylor tailのmajorantだけをexact `Rational`で反復乗算していた。`|z|`が8へ近づくと，`pi`のdyadic上界を含むRationalの分子・分母が各項で成長し，20桁`7+I`で約1.97 sまでwall timeが増えていた。tail majorantに必要なのは厳密な値ではなく「真のmajorantを必ず含む上界」なので，これを外向き丸め`RealInterval`へ変更した。固定係数`z^4*pi^2/4`もcomplex intervalとして一度だけmaterializeし，term更新で再利用する。tail判定は収束域へ入った後も毎項行う必要がないため8項ごとに限定した。剰余を採用する際はinterval上端をRationalへ変換して従来と同じ明示的tail inflationを行うため，certified contractは変わらない。

大きな正実数`1F1`では従来，termとpartial sumを最後までexact `Rational`で保持していた。`z=160`では値が`e^z`級に成長するため，項数約1000よりもRational normalization / gcd / allocationが支配し，約1.08 sを要していた。`a>0`,`b>0`のmoderate パラメータかつ`16<=z<=160`では全項が正で相殺しないため，値の指数成長分を`2 ceil(z)` bitのguardとして先に確保し，term/sumを外向き丸め`RealInterval`へ移した。

また従来は`n>=6|z|`まで待ってから「以後の比<=1/2」としてtailを押さえていた。正パラメータでは，`N=n+1`以降について

```text
r_j = z (a+j) / ((b+j)(j+1))
    <= z (1+a/N) / (N+1) = Q
```

が成り立つ。`Q<1`なら，まだpartial sumへ加えていない`next`を含む残差全体を`|next|/(1-Q)`で上から押さえられる。この判定を16項ごとに行うことで，`z=160`のCertifiedRefinementは約964から約513へ減り，tail proof自身のRational costも小さくした。

2026-08-25 GCC Release / LTO offでの再測定は次の通り。絶対時間はCPU/compiler依存なので，採用判断は比率とwork量を主に見る。

| workload | 旧warm代表値 | 最適化後代表値 | 備考 |
| --- | ---: | ---: | --- |
| `N[fresnelc[7+I],20]` | 約1.97 s | 約14 ms | この時点では`\|z\|<8`境界を維持 |
| `N[1F1[1/2,5/4,64],20]` | 約80 ms | 約2.5 ms | 正実数interval fast path |
| `N[1F1[1/2,5/4,128],20]` | 約619 ms | 約4.5 ms | 同上 |
| `N[1F1[1/2,5/4,160],20]` | 約1.08 s | 約6 ms | この時点では`\|z\|<=160`境界を維持 |

precision sweepでも`1F1[...,160]`は10/20/40/80桁で概ね5/6/7/10 ms級，complex Fresnel `7+I`は概ね20/14/25/36 ms級で，旧秒級のalgorithmic cliffは消えた。complex Fresnelでは`6+I -> 7+I`の項数増加自体は残るが，wall timeは数十ms以下であり，exact-majorant representationが作っていた病的な崖とは性質が異なる。

なおpositive `1F1` fast pathは実験上`z=512`でも数十msでcertifyできたが，exact pointだけ境界を拡大するとfinite-precision/complex 計算基盤との意味論上の非対称性が生じる。今回は性能修正に限定し，公開計算量制限付き policyは`|z|<=160`のままとした。境界拡大は別途InformationEnclosureを含む横断監査後に判断する。

# 21. elliptic F/E/Pi の near-boundary series cliff 除去

performance-cliff監査では，20桁・`phi=1/2`で`m=9/10`の`ellipticF`が約162 ms，`ellipticE`が約206 ms，`n=9/10`の`ellipticPi[...,m=1/3]`が約307 msまで増えていた。当初はAGM / Carlson symmetric formsへの置換を候補としたが，profile上の第一原因はLegendre級数そのものではなく，**保証用のexact `Rational`状態を数百項にわたり反復更新していたこと**だった。

旧実装では，

```text
c_k
m^k
Pi用q_k
tail用r^k
```

をexact `Rational`として保持していた。`r=max(|m|,|n|)`が9/10へ近づくとtail証明まで数百項を要し，分子・分母bit長，GCD正規化，allocationが項数とともに増える。一方，級数和と三角函数値は既に`RealInterval`で保証されているため，これら中間量だけexact分数で保持する数学的必要はない。

そこで係数，パラメータ冪，Piのcombined coefficient，tail ratio powerを作業precisionの外向き`RealInterval`へ移した。tail判定で必要なのは真の剰余を上から押さえることであり，採用時にはinterval上端をRationalへ戻して従来どおり対称誤差を最終値へ加える。certified contractは変更しない。tail certificateは8項ごとに評価し，毎項の除算・上界materializationも削減した。

さらに振幅が`|phi|<=Pi/2`内にあることを保証できる場合，

```text
I_k(phi) = Integral[sin(t)^(2k), {t,0,phi}]
|I_k(phi)| <= |phi| sin(|phi|)^(2k)
```

を利用する。従来は一律`|I_k|<=|phi|`としていたため，`phi=1/2`でも収束率を`r=0.9`として扱っていた。新しいtail proofでは実効比を

```text
rho = r sin(phi)^2
```

まで下げられる。`phi=1/2`では`sin(phi)^2`が約0.23なので，`r=0.9`でも`rho`は約0.21となり，必要項数が大幅に減る。`|phi|<=Pi/2`を証明できない場合は従来の`rho=r`へ自動代替経路するため，広いamplitudeの意味論は変えない。

2026-08-26，GCC Release / LTO off，20桁での代表値は次の通り。

| workload                   |   監査時 | 最適化後 | CertifiedRefinement |
| -------------------------- | -------: | -------: | ------------------: |
| `ellipticF[1/2,9/10]`      | 約162 ms | 約2.3 ms |           799 → 102 |
| `ellipticE[1/2,9/10]`      | 約206 ms | 約2.6 ms |           799 → 102 |
| `ellipticPi[9/10,1/2,1/3]` | 約307 ms | 約2.8 ms |           862 → 110 |

`phi=3/2`では`sin(phi)^2`が1に近づくため，20桁でF/Eは約26 ms，Piは約28 msまで増える。これはcomplete caseへ近づくことで級数の実効収束率自体が悪化する残存costであり，旧Rational肥大化cliffとは別物である。performance-cliff runnerには`m/n=9/10`のamplitude sweepと10/20/40/80桁precision sweepを追加し，この残存領域を継続監視する。

Carlson symmetric formsはLegendre楕円積分の有力な将来計算基盤だが，今回の対象である`phi=1/2, m/n→9/10`の秒未満cliffは既存級数の証明表現を直すだけで数msまで落ちた。新計算基盤は実装・証明・branch/amplitude reductionの複雑さを伴うため，**今回のperformance修正には採用しない**。complete case近傍や将来`|m|>9/10`へ対応領域を広げる段階で，改めてcrossoverを測る。

# 22. zeta Euler–Maclaurin の certification / representation cliff 除去

performance-cliff監査では，real zetaの一般非整数点が20桁で約4.2万`CertifiedRefinement`，complex zetaが0.5秒前後を要していた。さらにcomplex側を80桁まで上げると約2.5秒へ跳ねる別のcliffも確認した。分解すると，Euler–Maclaurin公式自体より，**plannerと有限Dirichlet和で同じ数学量を繰り返し超越評価・exact Rational再構成していたこと**が支配的だった。

real plannerでは，各`k`についてrising factorial，`(2k)!`，`N^(-s-2k+1)`を先頭から再計算していた。これらを隣接`k`間の漸化更新へ変更し，`N^-2`を一度作ってpowerを更新する。有限Dirichlet和`sum n^-s`では，正整数底に対する完全乗法性

```text
(ab)^(-s) = a^(-s) b^(-s)    (a,b>0)
```

を利用し，素数底だけ`Log -> Exp`でcertifyし，合成数は既計算項のinterval積から構成する。整数指数は従来どおりexact power，半整数指数はexact整数冪とcertified `sqrt`へ分解し，例えば`n^(-3/2)=1/(n sqrt(n))`を毎回`Log/Exp`へ送らない。

complex plannerの高precision cliffはさらに明確で，Bernoulli remainder bound中の`(2Pi)^(-2k)`を高精度dyadic Piからexact `Rational`へ変換し，`rationalPower`で毎`k`巨大化させていた。80桁代表点ではwall timeのほぼ全てがplannerだけで消費されていた。必要なのは剰余の厳密値ではなく保証上界なので，

```text
rising bound
N^(1-sigma-2k)
(2Pi)^(-2k)
```

を固定working precisionの外向き`RealInterval`漸化へ変更した。採用判定時だけ上端をRationalへ戻し，従来と同じ明示的remainder inflationへ渡す。complex Euler–Maclaurin correction本体もrising factorial，factorial，`N^-2`を隣接項で更新する。

旧correction上限`K<=40`では80桁代表点が`N=64`であと数項足りず，`N=128,K=32`へ飛んで有限和の超越評価数を倍増させていた。interval planner化後は高次項の追加costが小さいため，上限を**40から48へ小幅に広げ**，同じ代表点を`N=64,K=44`で閉じる。これは対応領域の拡張ではなく，同じEuler–Maclaurin 計算基盤内のwork balance調整である。

2026-08-26 GCC Release / LTO offでの代表値は次の通り。

| workload            |              監査時 |   最適化後 |
| ------------------- | ------------------: | ---------: |
| `N[zeta[3/2],20]`   |            約155 ms | 約10–12 ms |
| `N[zeta[12/5],20]`  |        約275–321 ms | 約40–46 ms |
| `N[zeta[3/2+I],20]` |         約0.7–1.2 s |    約14 ms |
| `N[zeta[3/2+I],80]` | 第一段後でも約2.5 s |    約73 ms |

complex `zeta[3/2+I]`の10/20/40/80桁sweepは概ね11/14/29/73 msとなり，performance-cliff runnerの4倍判定を外れた。real一般点でも，従来約4.2万に張り付いていたrefinement量とwall timeの双方が大きく減った。

この最適化を導入した時点ではcertified 計算基盤領域を`Re(s)>1`のまま変更せず，主値分岐，`s=1` 極，`PrecisionInsufficient` / `CertifiedBackendUnsupported`の分類も維持した。その後のcapability監査でcritical stripと左半平面へcertified 解析接続を拡張したが，ここで導入した完全乗法性の再利用は正整数底に限定されたままであり，complex logarithmの分岐を跨ぐ変形は導入していない。

今回のFresnel，1F1，elliptic，zetaを通して共通していた教訓は，**certified算法では数学的なexactnessが必要な量と，保証上界だけ必要な量を分離すること**である。後者を無条件にnormalized `Rational`で保持すると，算法本来の収束より先にrepresentation costがperformance cliffを作る。今後の特殊函数監査でも，まずplanner / majorant / tail certificateの表現をprofileしてから新計算基盤導入を判断する。

# 23. complex 2F1 解析接続 / li の Arg・Gamma representation cliff 除去

performance-cliff監査で残っていた代表hotspotは，20桁のcomplex `2F1` 主値 `1/z` 解析接続が約0.96 s，`li[-2+I]`が約1.03 sであった。両者を分解すると，`2F1`のGamma個数や`li=Ei(Log(z))`という合成そのものより，**主値 complex Logが内部で使うcertified `atan`のexact Rational級数**が共通の支配要因だった。

`enclosePrincipalArgument`は象限を確定した後，`atan(x/y)`等へ帰着する。`atan`は`|x|<=1/2`へ縮約後，交代級数

```text
atan(x) = x - x^3/3 + x^5/5 - ...
```

で隣接部分和を使って厳密に囲っていた。この証明自体は簡潔だが，complex intervalの除算で得た端点は高precisionのdyadic Rationalである。これをexact Rationalのまま`x^(2k+1)`へ反復すると，分母bit長が項ごとに増え，非dyadic入力由来のArgだけが数百ms〜秒級へ膨らんでいた。`Log(-4.6-2I)`や`Ei(0.8+2.7I)`が典型であり，`Log(1+I)`のような単純比では問題が見えにくかった。

新実装では，入力端点`x`自体はexact Rationalとして受けるが，級数のpower / partial sum / next termは`precision+32` bitの外向き`RealInterval`で更新する。交代級数の真値はexactな隣接部分和の間にあるため，それぞれを包含する二つのintervalのhullも厳密なenclosureである。停止判定は次項intervalの上端が従来閾値以下であることを要求する。このためbranch/Argの証明規約を弱めず，normalized Rationalの分母成長だけを除去できる。

complex `Ei`ではtail majorantも同じ方針へ揃えた。項絶対値のexact Rational値は不要なので，majorantを固定working precisionのoutward `RealInterval`で伝播し，tail certificateは8項ごとに評価する。今回の`li`高速化の主因はatan側だが，Ei側にも同じrepresentation cliffを残さないための整理である。

`2F1` 解析接続ではさらに二つの無駄を除去した。

1. connection coefficientに現れる7個のGammaのうち，`a`，`b`，`b-a`，`a-b`等がprovably realならcomplex Gammaへ送らずreal Gamma 計算基盤を使い，結果だけComplexIntervalへliftする。
2. complex GammaのStirling remainder証明では，shift後の`Re(z)`下端が非dyadic入力由来の高precision dyadicになることがある。従来はそのRationalを`lower^(2n-1)`へexact冪乗し，80桁付近で数百msのplanner costを作っていた。remainderに必要なのは正の下界だけなので，`floor(lower)`を保守的な整数下界として用いる。これは真の`Re(z)`以下であるため剰余上界は安全側にしか動かず，巨大dyadic分母を完全に避けられる。

2026-08-26，GCC Release / LTO off，3回測定warm側では以下となった。

| workload                        |   監査時 |  最適化後 |    改善 |
| ------------------------------- | -------: | --------: | ------: |
| `N[2F1[3.4,5.6,4+I,4.6+2I],20]` | 約0.96 s |  約100 ms | 約9.6倍 |
| `N[li[-2+I],20]`                | 約1.03 s | 約13.5 ms |  約76倍 |
| `N[li[-2],20]`                  | 約575 ms |    約9 ms |  約64倍 |

precision sweepは次の通り。

```text
2F1 continuation : 10/20/40/80 digits ~= 77 / 100 / 177 / 485 ms
complex li        : 10/20/40/80 digits ~= 10 / 13.5 / 19 / 36.5 ms
```

`2F1` 解析接続は通常のunit-disk Gauss seriesより依然重い。これはGamma係数7個，内部`2F1` 2本，主値 power 2本を持つconnection formula自体の固定costであり，旧Rational/Arg cliffとは別物である。80桁でも40桁比4倍未満へ収まり，現時点ではより複雑なGamma-ratio専用計算基盤を導入する根拠は弱い。

公開`2F1` 計算量制限付き領域，主値 `1/z` connection条件，`li`の主値 `Log -> Ei`定義，分岐切断 / 極 / `PrecisionInsufficient`分類は変更していない。performance-cliff runnerには両者の10/20/40/80桁sweepを追加し，Arg/Gamma経路の再発を継続監視する。

### 2F1の旧9/10 work boundary撤去（2026-08-26）

初回performance-cliff監査では，`2F1`のGauss級数を`|z|<=9/10`へ制限していた。これは数学的な収束境界ではなく，当時のexact `Rational` certification stateで`z->1`へ近づくと分子・分母が肥大化することへのwork policyだった。その後，Gauss級数の項・tail証明をfixed-working-precisionのoutward intervalへ移したため，この閾値の根拠を再測定した。

GCC Release / LTO-offの20桁CLI監査では，`2F1[1/2,1/3,5/4,z]`が`z=19/20`で約0.07 s，`49/50`で約0.18 s，`99/100`で約0.41 sまでcertifyできた。`999/1000`では約5.3 sまで増えるため，unit circleへの接近には依然として本質的な収束costがある。しかしこれは連続的な性能悪化であり，`9/10`に算法上の不連続点はない。

したがって固定`9/10`閾値は撤去し，Gauss級数は真の収束域`|z|<1`をそのまま試行する。計算量制限付きはseries term capと共通`EvaluationBudget`へ委ねる。有限precisionのmagnitude intervalが`|z|=1`を跨ぐ場合は`PrecisionInsufficient`，`z=1`を証明し`Re(c-a-b)>0`ならGauss summationで閉じる。それ以外の`|z|=1`点はパラメータ依存のboundary formulaが必要なため現時点では`CertifiedBackendUnsupported`とする。`|z|>1`の主値 `1/z` 解析接続は従来どおりである。

これは，**過去のperformance workaroundとして導入したhard 閾値は，representation cliffを除去した後に必ず再測定し，数学的境界まで戻せるなら撤去する**という方針の例である。

## Lambert W 複素分岐 の高precision縮小写像

任意整数branch対応後，`N[lambertw[2,1],100]`で`Log_k(z)-Log(w)`のrectangle包含反復が100桁まで逐次`certified Log`を再評価し，約2.7秒まで伸びるcliffが見つかった。候補点の探索を線形収束のLog固定点反復から`w exp(w)-z=0`に対するNewton反復へ変更し，最終証明は従来どおり分岐付きLog写像のBanach disk包含で行う。候補探索値そのものを証明には使わず，disk半径がrequested precisionを満たす場合だけ早期returnする。

GCC Release / LTO-offでは100桁`lambertw[2,1]`が約2.7秒から約0.07秒へ低下した。低精度の近似値だけを早期返却しないよう，certified radiusにrequested-width gateを置くことが重要である。

# 24. Complex Root isolation / refinement の対称seed・再isolation cliff除去

`N[root[...,k,Complex],p]`を再監査すると，real Rootよりcomplex Rootで大きなperformance cliffが残っていた。特に疎で対称な多項式では，次数そのものより入力shapeに依存して病的な停滞が起きていた。

旧実装のDurand–Kerner候補生成は，全rootを同一半径かつ完全な等角度配置へ置いていた。このseedは一般多項式では十分だが，`x^n-a`のような回転対称性を持つ多項式ではiterationを対称な軌道へ閉じ込めることがある。例えば`root[{-2,0,0,0,1},1,Complex]`は10秒超，`x^8-2` / `x^10-2`も数秒級またはtimeoutとなった。

候補seedに決定論的な小さい角度・半径摂動を入れ，完全対称性だけを壊す。seedはroot candidateの探索にしか使わず，採用するroot diskは従来どおりexact Rational中心・半径に対するRouche判定で一意根性を証明するため，数学的意味論には影響しない。

さらに，exact `Root` Callは既に内部`algebraicValue`としてcanonical `AlgebraicNumber`とisolating diskを保持していたが，CertifiedEvaluatorはこれを無視して係数からComplexAlgebraicNumberを再構築していた。`N[Out[-1],p]`等ではexact評価済みdiskを直接再利用するよう変更した。

ComplexAlgebraicNumber::refinedも，旧実装では対象rootを1個だけ細分化したい場合でも`isolateComplexDisks`を再実行し，全rootのDurand–Kerner候補生成・Rouche証明・orderingをやり直していた。新実装は既存isolating diskの中心からNewtonで候補を磨き，

```text
new disk is certified unique by Rouche
AND
new disk is contained in old certified disk
```

を満たす場合だけ単根local refinementとして採用する。旧diskが一意根を含み，新diskがその内部に完全包含されるため，同じroot identityを保持できる。local proofが閉じない場合だけ従来の全根re-isolationへ切り替える。

local Rouché証明でも一つrepresentation cliffが見つかった。Newtonで得た高precision BigFloat中心をそのままexact Rationalへ変換すると，真値では0である成分が極小dyadicとして残る場合がある。`x^8-2`の純虚根ではこの「ほぼ0」の巨大分母がTaylor展開へ入り，単根refinementだけで秒級になっていた。証明中心はrootそのものである必要はなく，十分近い有理点でよいので，working precisionの半分へdyadic量子化し，その精度以下の微小成分は0へ正規化する。最終diskは改めてRoucheで証明するため，この正規化は候補中心の選択にしか影響しない。

2026-08-27，GCC Release / LTO off，`--algebraic-root-cliffs 2`の平均値は次の通り。

| case                   | initial isolation | refine 80 bit | refine 320 bit |
| ---------------------- | ----------------: | ------------: | -------------: |
| `x^2-2`                |            4.7 ms |        0.8 ms |         1.0 ms |
| `x^4-2`                |           13.9 ms |        3.2 ms |         9.7 ms |
| `x^5-2`                |           46.8 ms |        5.9 ms |        27.1 ms |
| `x^8-2`                |           68.7 ms |       10.1 ms |        47.3 ms |
| `x^10-2`               |          226.8 ms |       52.7 ms |       169.4 ms |
| `x^12-2`               |          332.3 ms |       95.2 ms |       286.8 ms |
| `x^16-2`               |          696.9 ms |      163.0 ms |       771.9 ms |
| generic dense degree 5 |           18.8 ms |        7.2 ms |        26.0 ms |
| close-root degree 4    |           28.4 ms |        0.3 ms |         2.1 ms |

旧`x^4-2`の10秒超pathologyは消え，少なくとも16次まではdegree増加に応じた連続的なcostへ戻った。一方，この時点の追加監査では`x^20-2`初期分離が約2.2秒，24次約4.4秒，32次約15.7秒まで増えたため，`maximumSupportedAlgebraicDegree=64`を一旦計算基盤 capability / resource safety capとして維持した。この判断は後続の高次数Rouché／candidate-generation監査で再評価され，現在は96へ更新されている。

この監査用に`mmCal.Benchmarks --algebraic-root-cliffs [iterations]`を追加した。初期全根分離と選択rootの80/320-bit local refinementを別列で測り，対称疎多項式，generic dense多項式，近接根caseを常設する。

# 25. `1F1` の固定 `|z|<=160` 境界撤去

Complex Root監査後にcertified特殊函数の残存hard limitを再点検したところ，`hypergeometric1F1`の`|z|<=160`は数学的境界ではなく，旧exact-`Rational`級数のrepresentation cliffを避けるために残っていたpolicy値であった。前段でlarge positive real `z`のterm / partial sumをoutward `RealInterval`へ移した後も，realの負側とcomplex側には固定閾値を外すための横断確認が残っていた。

今回，realのlarge-`|z|`経路を符号非依存のmajorantへ一般化した。`N>=2|a|,2|b|`なら，将来の級数項比は

```text
|t_(j+1)/t_j| <= 3 |z| / (N+1)
```

で一様に上から押さえられる。したがって右辺`Q`が1未満へ入った時点で，まだ加えていない次項を含むtailを`|next|/(1-Q)`で保証できる。term / sum自体はguard付きdyadic intervalで保持するため，負実数側の大きな相殺でもexact Rationalの分母肥大化へ戻らない。

complex 計算基盤も同じmajorantを使う。旧実装に残っていた暫定`|z|<=8000` cancellation-guard 閾値は撤去し，real/complexとも固定magnitude 閾値を持たない。代わりに，`2|a|`，`2|b|`，`3|z|`に基づく保守的なtail-certification開始点が最大250000 termsのalgorithm budget内へ入るかを事前判定し，series loopは共通`EvaluationBudget::CertifiedRefinement`でも停止する。これは巨大入力を無制限に計算することとは異なり，**入力値そのものではなく必要work量で計算量制限付きを決める**設計である。

finite-precision入力も同じ経路へ通すが，`InformationEnclosure`をpointへ縮退させない。例えば`N[hypergeometric1F1[1/2,5/4,N[161,5]],20]`は旧160境界を越えて数値を返せる一方，入力5桁から20桁を捏造せず，結果precisionは約2桁に留まる。

2026-08-27，GCC Release / LTO off，`--performance-cliffs 3`のwarm値は次の通りである。

| workload                    |      warm |
| --------------------------- | --------: |
| `N[1F1[1/2,5/4,64],20]`     |  約2.1 ms |
| `N[1F1[1/2,5/4,160],20]`    |  約6.0 ms |
| `N[1F1[1/2,5/4,161],20]`    |  約5.0 ms |
| `N[1F1[1/2,5/4,256],20]`    |  約9.7 ms |
| `N[1F1[1/2,5/4,512],20]`    | 約26.6 ms |
| `N[1F1[1/2,5/4,1000],20]`   |   約98 ms |
| `N[1F1[1/2,5/4,-512],20]`   |   約27 ms |
| `N[1F1[1/2,5/4,256+I],20]`  |   約53 ms |
| `N[1F1[1/2,5/4,-256+I],20]` |   約44 ms |

`z=512`のprecision sweepは10/20/40/80桁で約25/30/28/32 ms，complex `256+I`は約41/43/65/117 msであり，旧160位置にも高precision側にも4倍cliffは現れなかった。

この結果から，`|z|<=160`を別の固定値へ拡張する案は採らず，閾値自体を削除した。今後1F1で大引数の真のcliffが現れた場合はKummer変換やasymptotic 計算基盤を検討するが，現時点で固定magnitude limitを復活させる理由はない。

# 27. 実数 Ei / Si / Ci の固定96境界撤去と保証付き漸近計算基盤

`1F1`の固定境界撤去後，real `Ei/Si/Ci`に残っていた96も再監査した。interval Taylor化後の`97`は既に数ms～十数msで評価でき，96は旧exact-`Rational`級数の相殺・表現costを避けるための値になっていた。ただし単純にTaylor範囲だけを延長すると，`Ci[1000]`では`gamma+log(x)`と巨大Taylor和の相殺を支えるためEulerGammaへ約2000bit級を要求し，別計算基盤の上限へ到達する。したがって固定閾値を別の値へ移すのではなく，大引数用の解析的に異なる計算基盤を追加した。

positive real `Si/Ci`はDLMF 6.12の補助函数

```text
f(x) ~ 1/x (1 - 2!/x^2 + 4!/x^4 - ...)
g(x) ~ 1/x^2 (1 - 3!/x^2 + 5!/x^4 - ...)

Si(x) = Pi/2 - f(x) cos(x) - g(x) sin(x)
Ci(x) = f(x) sin(x) - g(x) cos(x)
```

を使う。正実軸では`f,g`の剰余は最初の未使用項以下で，しかも同符号であるため，最適打切り前のtermをoutward `RealInterval`で反復し，最初の未使用項を片側区間として足すだけでrigorous enclosureを作れる。項が再増大した時点で要求幅へ届いていなければTaylorへ切り替える。これによりEulerGammaを大引数Ciの値構成から完全に外せる。

positive real `Ei`は

```text
Ei(x) ~ exp(x)/x (1 + 1!/x + 2!/x^2 + ...)
```

を使う。DLMF 6.12.2の剰余上界`(1+chi(n+1)) * nextTerm`を直接実装した。整数argumentの`chi`は`chi(2)=2`, `chi(3)=3 Pi/4`, `chi(t+2)=chi(t)(t+2)/(t+1)`で保証区間を漸化更新でき，Gamma ratioを別途評価しない。`Ei(x)/(exp(x)/x)`は1程度なので，bracket側を相対精度相当までcertifyしてからprefactorを掛けることで巨大値でもsignificant-digits契約を保つ。負側は既存の`Ei(-x)=-E1(x)`保証付き漸近展開を維持する。

GCC Release / LTO-off，`--performance-cliffs 2`の20桁warm代表値は次の通り。

| workload         |     warm | CertifiedRefinement |
| ---------------- | -------: | ------------------: |
| `N[Ei[96],20]`   | 約6.0 ms |                 233 |
| `N[Ei[256],20]`  | 約1.8 ms |                 175 |
| `N[Ei[512],20]`  | 約1.5 ms |                 167 |
| `N[Si[1000],20]` | 約1.4 ms |                  73 |
| `N[Ci[1000],20]` | 約2.0 ms |                  73 |

さらに`N[Si[10000],20]`, `N[Ci[10000],20]`, `N[Ei[10000]/exp[10000],20]`も保証付きで短時間に完走する。旧96境界を別のmagnitude 閾値へ置換せず，**漸近剰余が要求精度へ届くか／Taylor term cap内か／共通EvaluationBudget内か**で計算基盤を選ぶ。

この監査では既存の`N[Ei[-64],20] -> 0.0`も発見した。負大引数Taylorの絶対tail基準が微小な真値に対して粗すぎたためで，`E1(x)>=exp(-x)/(x+1)`から相対精度相当のtail targetを作ることで`-2.46796855945...e-30`を保持するよう修正した。

# 28. complex Fresnel の固定 `|z|<8` 境界撤去

complex Fresnelのexact-majorant representation cliffを除去した後，旧`|z|<8`をcapability limitとして再監査した。境界だけ外してMaclaurinを試すと，`8+I`約0.02秒，`10+I`約0.08秒，`20+I`約1.2秒まではcertifyできるが，`32+I`では8秒を越える。したがって`8`自体は古いperformance workaroundだが，大きな軸近傍complex argumentには依然として振動級数の相殺cliffが存在する。

固定境界を別の値へ移す代わりに，DLMF 7.12の複素Fresnel補助函数`f/g`漸近展開を追加した。

```text
C(z) = 1/2 + f(z) sin(Pi z^2/2) - g(z) cos(Pi z^2/2)
S(z) = 1/2 - f(z) cos(Pi z^2/2) - g(z) sin(Pi z^2/2)
```

`f/g`は`1/(Pi z)`を核にした逆冪漸近級数で，`|arg z|<Pi/8`では剰余の絶対値を最初の未使用項で上から押さえられる。実装では角度を近似してwedge所属を推測せず，区間端点から

```text
Re(z) > 0
4 |Im(z)| <= Re(z)
```

を証明できる場合だけ漸近計算基盤を使う。`atan(1/4)<Pi/8`なのでDLMFの強い剰余条件の内部に確実に入る。さらにFresnel級数の90度回転対称性

```text
C(i z) = i C(z)
S(i z) = -i S(z)
```

を使い，4通りのquarter turnのどれかで同じwedgeへ入る軸近傍入力を共通計算基盤へ送る。これにより正負実軸・正負虚軸を別実装にしない。対角方向等でwedgeを証明できない場合は，Fresnel C/Sがentireであることを利用してguard付きMaclaurinへ切り替える。

`f/g`の各termは`ComplexInterval`で反復し，次の未使用項のabsolute upper boundをtail radiusとする。最終C/Sへはcertified complex `sin/cos`のmagnitudeを掛けた誤差上界まで含めてrequested targetへ届くことを確認してから採用する。漸近項が再増大し始めても要求幅へ届かない場合はseriesへ戻すため，発散漸近級数を無理に延長しない。

2026-08-27，GCC Release / LTO off，`--performance-cliffs 3`のwarm代表値は次の通り。

| workload                |       warm |
| ----------------------- | ---------: |
| `N[fresnelc[7+I],20]`   | 約15–24 ms |
| `N[fresnelc[8+I],20]`   |    約25 ms |
| `N[fresnelc[20+I],20]`  |    約11 ms |
| `N[fresnelc[32+I],20]`  |    約10 ms |
| `N[fresnelc[1+32I],20]` |    約12 ms |
| `N[fresnelc[8+8I],20]`  |    約13 ms |

`32+I`のprecision sweepは10/20/40/80桁で約8/10/16/35 msで，旧8境界の外側にも高precision側にも固定閾値を必要とするcliffはない。`8+8I`のような対角入力は漸近wedgeへ入らないが，series自体が速いためhybridとして相補的に機能する。

finite-precision入力でもInformationEnclosureは維持される。`N[fresnelc[N[8,5]+I],20]`は約2桁，`N[fresnelc[32+I*N[1,5]],20]`は約1桁に留まり，asymptotic dispatchがhidden exact centerを再取得しないことを確認した。

この結果，complex Fresnelの固定`|z|<8`は削除し，現在の計算量制限付き policyは**漸近剰余証明，Maclaurin term cap，共通EvaluationBudget**へ一本化した。

# 29. 楕円積分の9/10 capability boundary撤去

`ellipticF/E/Pi`には長く`|m|<=9/10`（`Pi`はさらに`|n|<=9/10`）という固定境界が残っていた。series内部のexact Rational肥大化をinterval化した後では0.9直前の計算は既に数ms級であり，この境界はstale performance workaroundになっていた。

今回，DLMF 19.25のLegendre形からCarlson symmetric formsへの変換と，DLMF 19.26のduplication formulaを用いるcertified `RF/RD/RJ` 計算基盤を追加した。Carlson残差は正実引数における単調性から上下を包含し，`RJ`の補正項で必要な`RC`は近接引数域で正則級数へ切り替えてinterval cancellationを避ける。exact `q Pi`振幅は`Pi`の近似値同士を除算せず，有理係数のままperiod還元する。

series/Carlson dispatchも固定0.9ではない。保証付き実効tail ratioが`<=9/10`ならseries fast pathを使う。したがって小振幅で`m,n`が1へ近づいても旧境界位置に人工的な性能段差を作らない。

2026-08-28，GCC Release / LTO offの20桁warm代表値：

| workload                    |      warm |
| --------------------------- | --------: |
| `ellipticF[1/2,0.90]`       |  約2.2 ms |
| `ellipticF[1/2,0.95]`       |  約2.6 ms |
| `ellipticF[1/2,0.99]`       |  約2.4 ms |
| `ellipticF[1/2,0.999]`      |  約2.5 ms |
| `ellipticPi[0.99,1/2,0.99]` |  約3.0 ms |
| `ellipticF[3/2,0.99]`       | 約26.7 ms |
| `ellipticE[3/2,0.99]`       | 約53.4 ms |
| `ellipticPi[0.99,3/2,0.99]` | 約75.2 ms |

Carlson経路の10/20/40/80桁precision sweepは，Fが約16/27/54/143 ms，Eが約32/53/115/317 ms，Piが約46/75/155/414 msで，precision増加に対して連続的に伸びる。

能力境界も固定パラメータ値から数学的条件へ移した。`m>1`/`n>1`でも，還元後の実積分路が分岐点/極より手前にあることを区間で証明できればcertifyする。特異点を跨ぐ場合や一般complex 解析接続は現計算基盤では`N::unsupported`を保持する。`ellipticE`のexact `m=1`は有限退化を専用処理する。

# 30. positive-order polylog の固定49/50境界撤去とDLMF 解析接続 / near-one 計算基盤

positive-order `polylog`には，旧exact-Rational級数のrepresentation cliffを避けるため`|z|<=49/50`という固定計算量制限付き境界が残っていた。まずこの境界を外して再測定したところ，数学的収束域`|z|<1`の内部でも，termを`z^k`の巨大exact Rationalとして保持することがunit point近傍の主なcostになっていた。

そこでdirect seriesは

```text
t_(k+1) = t_k z (k/(k+1))^s
```

のoutward `RealInterval` / `ComplexInterval` recurrenceへ変更した。tailは従来どおり将来項比の一様上界`|z|`から`|next|/(1-|z|)`で包含する。これによりexact Rationalの分子・分母成長を除去したが，`Li_2(0.999)`等では級数そのものの自然な収束遅延が残った。

`Li_2`にはDLMF 25.12.3，25.12.4，25.12.6の主値 connection formulaを追加した。実装は公式を無条件に適用せず，入力区間全体が対応する分岐切断を避けることを証明でき，かつ変換先のargument magnitudeが十分縮小する場合だけ採用する。

```text
Li_2(z) + Li_2(z/(z-1)) = -1/2 Log(1-z)^2
Li_2(z) + Li_2(1/z) = -Pi^2/6 - 1/2 Log(-z)^2
Li_2(x) + Li_2(1-x) = Pi^2/6 - log(x)log(1-x), 0<x<1
```

これにより正実unit point近傍，負実軸，`z=I`，分岐切断を避けた`|z|>1`等を小さいargumentへ移せる。exact `z>1`の正実軸は主値 分岐切断そのものであり，上側／下側の境界値をmmCalが勝手に選ぶべきではないため，従来どおり未評価に留める。

さらに`Li_3`以上では`z=1`自体が有限であるにもかかわらずdirect seriesだけでは`z=0.999`で秒級へ入った。DLMF 25.12.12の非整数order公式を`s -> n`の正整数極限へ取り，正実`0<z<1`，`mu=log(z)<0`に対して

```text
Li_n(exp(mu))
 = mu^(n-1)/(n-1)! (H_(n-1) - log(-mu))
   + sum_{k>=0, k!=n-1} zeta(n-k) mu^k/k!
```

をcertified near-one 計算基盤として追加した。負整数zeta値は

```text
zeta(1-2r) = -B_(2r)/(2r)
```

へexact化する。残差は

```text
|B_(2r)| = 2 (2r)! zeta(2r)/(2Pi)^(2r)
zeta(2r) < 2
```

から作るmajorantを用い，後続majorantの比を`(|mu|/(2Pi))^2`以下で押さえて幾何tailとして包含する。したがってDLMFの漸近式を「項が小さそうだから」打ち切るのではなく，明示的なremainder certificateを持つ。

このnear-one経路は能力境界ではない。GCC Release / LTO-offでのcrossover監査から，現在はorder 3～12かつ`|mu|<=1/20`の場合だけfast pathとして使い，それ以外は既存のunit-disk seriesへ切り替える。高orderでは`1/k^s`自体が強く効いてdirect seriesが既に速いためである。finite-precision実入力も同じ式を`RealInterval`のまま通し，midpointやhidden exact valueへ戻さない。

2026-08-28，GCC Release / LTO offの20桁代表値は次の通り。

| workload      |       旧監査 | 新計算基盤 |
| ------------- | -----------: | --------: |
| `Li_2(0.999)` |      約1.2 s |  約2.9 ms |
| `Li_3(0.999)` |     約0.96 s |   約11 ms |
| `Li_4(0.999)` |     約0.80 s |   約13 ms |
| `Li_8(0.999)` |     約0.24 s |   約24 ms |
| `Li_2(I)`     | series境界外 |   約35 ms |
| `Li_2(-2)`    | series境界外 |    約5 ms |
| `Li_2(2+I)`   | series境界外 |   約36 ms |

`Li_2(0.999)`の10/20/40/80桁は約2.3/2.9/4.2/7.3 msで連続的に伸びる。`Li_3(0.999)`も80桁でCLI込み約0.09 s程度である。`N[polylog[3,N[999/1000,8]],20]`は旧interval seriesで約2.9 sだったが，near-one interval経路では約0.01 sとなり，表示は`1.2004154`のままで入力8桁以上のprecisionを発明しない。

以上から固定`|z|<=49/50`は撤去し，現在の計算量制限付き contractは**unit-disk series term cap，共通`EvaluationBudget`，DLMF connection formulaの分岐証明，near-one Bernoulli remainder proof**で構成する。

# 31. complex Ei/Ci の固定512/128境界撤去

複素`Ei`には`|z|<=512`，複素`Ci`には`|z|<=128`という旧計算量制限付き境界が残っていた。単純に上限判定だけを削除すると，GCC Release / LTO offの20桁で`Ei[513+I]`は約0.7 s，`Ci[129+I]`は約1.5 sへ入り，`1000+I`ではcertified refinement budgetを使い切った。したがって旧境界はstaleな数字ではあったが，背後のTaylor cancellation cliffは実在していた。

大引数用にDLMF 6.12.1の

```text
E1(z) ~ exp(-z)/z (1 - 1!/z + 2!/z^2 - ...)
```

を`ComplexInterval`で反復するcertified kernelを追加した。`|arg z|<=Pi/2`では最初の未使用項，左半平面ではDLMFのsector boundに従い`csc(|arg z|)`を掛けたmajorantで剰余を包含する。漸近項が再増大する前にrequested whole-complex widthを証明できない場合は採用せず，収束seriesへ切り替える。

複素`Ei`は主値 identity

```text
Ei(z) = -E1(-z) + Log(z) - Log(-z)
```

でbranch correctionを主値 `Log`へ集約する。ただし正実軸はEiのcutではない一方，`E1(-z)`側では負実軸cutへ写る。有限precision虚部が0を跨ぐ正実軸近傍についてはreal `Ei(x)`をanchorにし，`Ei'(z)=exp(z)/z`を縦方向に積分したrigorous radius boundで入力幅をそのまま伝播する。これにより`N[Ei[1000+I*N[0,5]],20]`をhidden exact zeroへ決めず，resource limitにも落とさない。

複素`Ci`は右半平面で

```text
Ci(z) = -1/2 (E1(i z) + E1(-i z))
```

を使用し，左半平面は`Ci(z)=Ci(-z)+Log(z)-Log(-z)`で主値分岐を保つ。純虚軸は`Chi`へ退化させる。高precisionで漸近最適打切りが要求幅へ届かない`Ci[140+I]`等はseriesへ戻るため，series側のtail majorantも巨大exact Rationalからoutward `RealInterval` recurrenceへ変更した。

2026-08-28，GCC Release / LTO off，`--performance-cliffs 1`の代表値：

| workload           |  20-digit |
| ------------------ | --------: |
| `N[Ei[128+I],20]`  | 約13.2 ms |
| `N[Ei[512+I],20]`  |  約8.2 ms |
| `N[Ei[513+I],20]`  |  約9.4 ms |
| `N[Ei[1000+I],20]` |  約9.2 ms |
| `N[Ei[513I],20]`   | 約10.8 ms |
| `N[Ci[120+I],20]`  | 約49.7 ms |
| `N[Ci[128+I],20]`  | 約44.2 ms |
| `N[Ci[129+I],20]`  | 約45.7 ms |
| `N[Ci[140+I],20]`  | 約41.9 ms |
| `N[Ci[1000+I],20]` | 約18.4 ms |

`Ei[1000+I]`の10/20/40/80桁は約7.3/8.5/12.3/28.6 ms，`Ci[1000+I]`は約11.7/19.6/40.9/114.1 msで連続的に伸びる。series 代替経路が必要な`Ci[140+I]`も20/50/100桁で概ね0.04/0.28/0.36 sに収まり，高precision側へ新しい秒級cliffを作らない。

以上から固定`512/128`は撤去し，現在の計算量制限付き contractを**分岐切断 proof，E1漸近remainder certificate，guard付きseries term cap，共通`EvaluationBudget`**へ置換した。

## Lambert W `-1/e` 分岐点 の局所縮小写像

complex `LambertW`の一般branch 計算基盤は`w exp(w)-z=0`のNewton候補と`Log_k(z)-Log(w)`のBanach証明を使うが，`w=-1`では導函数が退化するため`-1/e`近傍だけ保証boxを閉じにくかった。DLMF 4.13.9_1の平方根局所構造を利用し，`u=W+1`, `q=e z+1`として `q=u^2 A(u)/2`，`A(0)=1` を用いる。反復は `p=sqrt(2q)` を一度固定して `u=±p/sqrt(A(u))` とし，`A`と`A'`をoutward complex interval級数で評価して写像包含と一様縮小率を証明する。

`sqrt(2q/A(u))`を反復ごとに直接評価する案は，区間sqrtのcut sideとrectangle dependencyで幅が約`1e-8`から縮まらないため不採用とした。`p`を固定し，`1/sqrt(A)`を`A-1`の保証付きbinomial級数で評価するとこのdependencyが消える。`W_0`は主値 `p`，`W_-1`は上側で`-p`，`W_1`は下側で`-p`へ接続する。候補生成にはDLMFのPuiseux先頭3項だけを使い，最終採用は局所写像の包含証明で行う。

極端な近傍では`z`全体を先に有限精度化してから`e z+1`を作ると分岐点との差がcancellationで失われるため，exact式の段階で`z+1/E`を簡約し，offsetを直接局所kernelへ渡す。停止判定も`W≈-1`全体の絶対幅ではなく`u=W+1`の各非零成分に対するrelative widthで行う。候補点はPuiseux seedから`q=u^2A(u)/2`のNewtonで精製し，候補精度に応じた小さい初期boxからBanach包含を開始する。Newton候補は証明値として採用せず，最終採用条件は従来どおり`T(B)⊂B`と一様縮小率である。

局所`A(u)`, `A'(u)`, `A(u)^(-1/2)`のtail majorantも巨大exact Rationalの累積を避け，現在項のoutward interval normと一様項比から幾何tailを直接証明する。2026-08-28のGCC Release / LTO-off runnerでは，100桁`-1/e+i*10^-2`約1.51 s，`-1/e+i*10^-12`約0.31 s，`-1/e-i*10^-160`約0.031 sである。実軸右側`-1/e+10^-12`もCLI起動込み約0.25 s（`W_-1`は約0.29 s）で，浅い近傍から極端な近傍まで秒級timeoutを作らず連続的に処理できる。

監査中にはLambert W外の共通性能バグも見つかった。`algebraicBinary`が左辺を非代数式と判定できる場合でも右辺の巨大complex Rationalを先に`AlgebraicNumber`へ変換していたため，`x-I/10^160`や`1/E-I/10^160`だけで秒級になることがあった。左側の変換失敗時に直ちにreturnする短絡評価へ変更し，不要な代数体構築を避けた。

# 32. Complex Root全根分離のNewton-polygon multi-radius候補

Complex Rootは既に，近似候補生成とexact証明を分離している。Durand-Kernerで中心候補を作り，各中心をNewtonで局所精製し，最後にexact RationalのRouché判定でdisk内の根がちょうど1個であることを証明する。残っていた弱点は初期配置であり，旧実装は全候補を1個のCauchy半径付近へ置いていた。根の絶対値が同程度なら問題になりにくいが，1本の多項式に桁違いのroot-radius群が共存すると候補移動が非常に大きくなる。

現在は非零係数について

```text
(k, log2(|a_k|))
```

の上側Newton polygonを作る。連続する頂点`k0 < k1`に対して`m=k1-k0`個の候補を

```text
r = (|a_k0 / a_k1|)^(1/m)
```

の近似半径を持つ円へ配置する。logarithmとradiusは候補生成専用の近似値であり，数学的証明には使わない。有限精度の凸包が不完全でも最終採用条件は従来のexact Rouché証明だけなので，影響は性能に限定され，誤ったRootをcertifyしない。定数項が0の場合や近似geometryを安全に作れない場合は旧Cauchy半径seedへ切り替える。

この初期値方針はFLINT/Arbの`acb_poly_find_roots`が文書化しているNewton-polygon circleによるdefault initial valuesと同系統である。mmCalは後段のexact Rational disk証明と`Re(z)+Pi Im(z)`による決定的orderingを独自に維持する。

常設benchmarkには

```text
(x^8 - 2^80)(x^8 - 2^-80)
= x^16 - (2^80 + 2^-80)x^8 + 1
```

を追加した。8根ずつが半径`2^-10`と`2^10`へ分かれるため，単一Cauchy半径seedの弱点を直接検出できる。2026-08-28 GCC Release / LTO offでは，旧単一Cauchy半径による初期全根分離が約10.79 sだったのに対し，Newton-polygon初期値では約0.36 sまで低下した。既存`x^16-2`は約0.91 sで従来水準を維持し，単一半径の対称familyを悪化させずscale-separated polynomialの病的挙動を除去している。

Newton-polygon seed導入時点ではComplex Rootのdegree 64 計算基盤上限を撤去しなかった。候補生成以外にもRouché certification，disk同士の分離確認，決定的orderingのcostが次数とともに増加していたためである。後続監査でRouché証明とDurand-Kerner停止条件を改善した結果，この上限は96へ再設定した。`--algebraic-root-cliffs`では従来の`x^n-2`系列，two-radius stress case，高次Solveを併記して継続監視する。

# 33. Advanced Integration の有限Fourier直接積分

`runAdvancedIntegrationTests`のwall timeをケース単位で監査すると，小次数の
`sin[u]^m cos[u]^n` gridが連続して0.5～1.1 s級となり，group全体の大半を占めていた。
重いのはderivative-back検証ではなく`integrate[...]`本体であり，従来経路は

```text
整数sin/cos冪を有限Fourier式へ展開
  -> 汎用simplify
  -> integrateCoreへ再投入
```

としていた。さらに有限Fourier familyへ到達する前にreverse-chain，Weierstrass，parts等の汎用候補探索も毎回通っていた。

現在は`trigonometric_polynomial`がExprだけでなく，

```text
共通argument
exact Rational coefficient
frequency
sin/cos種別
```

を持つ`TrigFourierExpansion`を返せるようにした。同一argumentの整数sin/cos冪で総次数3以上なら，汎用候補探索より前にこの構造を生成する。argumentの導函数がvariable非依存かつ非零と証明できる場合は，共通`du/dx`を1回だけ求め，

```text
c cos(k u) -> c sin(k u)/(k u')
c sin(k u) -> -c cos(k u)/(k u')
constant   -> constant*x
```

を直接構成する。session angle semanticsのscaleは従来の`TrigArgument`情報をそのまま使うため，Rad/Deg/Gradの意味論は変えない。

非線形argumentではこのfast pathを無理に使わない。例えば`sin[2x^2]^4`は共通`du/dx`が定数でないため直接積分を拒否し，従来どおりFourier Exprを汎用integratorへ戻す。これによりquadratic phaseがFresnelへ落ちる既存能力を維持する。平方や`sin[x]cos[x]`等の既存compact解も総次数3未満として従来経路を優先する。

2026-08-28，GCC Release / LTO offの代表値では，`runAdvancedIntegrationTests`は約14.9 sから約2.2 sへ低下した。小次数gridの個別`integrate`は概ね0.05～0.09 s以下となり，derivative-back側は引き続き数ms級である。256乗の有限Fourier testや非線形Fresnel compositionを含む全Advanced Integration回帰は維持している。

# 34. CalculusKnowledge のalgebraic field保持とLambert W disk certificate

`runCalculusKnowledgeTests`の監査では，主要hotspotが二つあった。

1. `(sqrt(2)+sqrt(3))^3`相当のalgebraic chained powerが約1.9 s。
2. `N[lambertw[-1,-1/E-I/10^8],20]`が約2.0～2.3 s。

前者では`sqrt(2)+sqrt(3)`のminimal polynomial

```text
x^4 - 10x^2 + 1
```

をQ上既約と現proverが証明できず，persistent `NumberFieldContext`を付与できないことが原因だった。small-prime reductionで「あるmod-p像が既約」ならQ上既約という十分条件は使っていたが，V4型の偶四次ではこの条件だけでは証明できない場合がある。

monic偶四次`x^4+b x^2+d`については，Q上のmonic quadratic factorizationを

```text
(x^2+p x+q)(x^2-p x+s)
```

と置くと`p(s-q)=0`である。したがって可約なのは，

```text
p=0  : b^2-4d がQの平方
q=s : d がQの平方で，いずれかの q=±sqrt(d) に対して 2q-b がQの平方
```

のいずれかの場合に限る。この完全なexact判定を追加し，上記四次を既約と証明できるようにした。

さらに既存field座標を持つreal `AlgebraicNumber`の整数冪は，binary exponentiation中に毎回visible `Root[minpoly,k]`へmaterializeせず，`AlgebraicElement`のQ(theta)座標のまま乗算・二乗・逆元を計算し，最後に一度だけminimal polynomialとroot isolationへ戻す。これにより`(sqrt(2)+sqrt(3))^2/^3`相当は約0.6/1.9 sから概ね0.14～0.17 sへ低下した。

後者のLambert Wでは，一般分岐の`Log_k(z)-Log(w)` 計算基盤がNewton候補の後にまずdisk Banach certificateを作っていた。存在・一意性そのものは既に証明できていたが，そのcertified diskの幅が内部guard targetへ数桁届かないだけで破棄し，高価なrectangle contractionを最大256回行っていた。

現在はdiskで`T(B) subset B`と縮小率`L<1`を証明できた時点でそのenclosureを返す。要求precisionへ十分狭いかは外側の`CertifiedEvaluator`が既に判定し，不足時だけ高precisionで再評価するため，計算基盤内部で同じtargetを二重に強制しない。証明条件自体は緩めていない。

この変更で`N[lambertw[-1,-1/E-I/10^8],20]`は約2.3 sから約0.05 sへ低下した。最終ケース間隔監査ではCalculusKnowledge内の最大値はcomplex Root isolation等の約0.31 s級となり，秒級の突出点は消えた。group全体は約7.8 sから約3.2 sへ低下した。

両最適化後も`runAdvancedIntegrationTests`と`runCalculusKnowledgeTests`は削減せず全件維持する。2026-08-28の最終GCC Release / LTO-offではinternal C++ regression 2820/2820 PASS，関連black-box 710/710 PASS，Certification Boundary Fuzzer 1000/1000 PASS（59 classification + 36 metamorphic）を確認した。

# 35. 一般高次 `solve` のmodular既約性証明と全根batch 正準化

`solve[x^16+x+1==0,x]`を監査すると，Complex Rootの候補生成やRouché証明そのものより，その前後に二つの構造的なperformance cliffがあった。

第一に，従来の`provenIrreducibleOverQ`は複数の小素数について「ある1個のmod-p像が元次数のまま既約ならQ上既約」という十分条件だけを試していた。`x^16+x+1`はQ上既約だが，旧prime集合では各mod-p像が分解するため証明できず，degree 16のKronecker完全探索へ代替経路していた。

現在はgood prime上で多項式がsquare-freeであることを確認したうえで，

```text
gcd(f, x^(p^k)-x)
```

のdegreeから各mod-p既約因子degreeの個数を復元する。Q上のproper factorがdegree `d`なら，leading degreeを失わないgood primeではそのmod-p像もdegree `d`を保ち，mod-p既約因子の一部を選んだdegree和として現れなければならない。そこで各primeについて`1..floor(n/2)`の可能factor-degree集合をsubset-sumで作り，複数prime間で共通部分を取る。共通候補degreeが空になった場合だけQ上既約と証明する。

`x^16+x+1`では代表的に，

```text
mod 2 : factor degrees 8 + 8   -> proper factor候補 {8}
mod 3 : factor degrees 1 + 15  -> proper factor候補 {1}
```

となり，共通候補がないためKronecker探索へ入る前にQ上既約と証明できる。この判定はprobable-prime型の推測ではなく，Gauss lemmaと有限体上のexact gcd/Frobenius計算に基づく一方向の証明である。square-free性やdegree保存を証明できないprimeは単に無視する。

第二に，Complex `solve`は最初に`ComplexAlgebraicNumber::isolateAll`で全根をcertifyした後，各rootをcanonical Rootへ変換するため`ComplexAlgebraicNumber::create`をroot数だけ呼んでいた。`create`は同じdefining polynomialの全根分離をもう一度行うため，一般degree `n`では実質的に全根分離を`1+n`回繰り返していた。

`ComplexAlgebraicNumber::canonicalizeAll`を追加し，既に得たall-root disk群を一括処理するようにした。元多項式が既約ならdiskとglobal root orderingをそのまま再利用する。degree 16以下の可約多項式ではfactorizationを一度だけ行い，各既約factorの根を一度ずつ分離して，元のordered isolating diskとのintersectionからminimal factorとfactor内root indexを対応付ける。対応を一意に証明できない場合だけ従来のper-root `create`へ切り替える。minimal-polynomial reductionの対象外であるより高次数では，従来`create`も元多項式を保持するため，証明済みall-root diskをそのまま再利用して無意味な再isolationを避ける。

2026-08-28，GCC Release / LTO offの代表値：

| workload | wall time |
| --- | ---: |
| `solve[x^6-3x^5-x^4+2x^3+2x^2-2x-1==0,x]` | 約0.046 s |
| `solve[(x^4-2)(x^4+1)==0,x]` | 約0.17--0.20 s |
| `solve[x^16+x+1==0,x]` | 約1.38--1.45 s |
| `solve[x^18+x+1==0,x]` | 約2.37 s |
| `solve[x^20+x+1==0,x]` | 約3.70 s |

監査時には`x^16+x+1`で候補生成が約70 ms，全16根のRouché certificationを含む初期全根分離が約1.5 sであり，根分離kernel自体は支配的な異常ではなかった。旧経路は既約性proof failure後のKronecker探索と，その後のper-root all-root再isolationが主因であった。`--algebraic-root-cliffs`には一般16次Solve，複数primeのfactor-degree証明を使う6次，batch 正準化を使う可約8次を常設し，Root isolationだけでなくSolve materializationまで継続監視する。


# 36. 高次Complex RootのRouché証明・候補反復・degree budget再監査

一般高次`solve`のbatch 正準化後，`solve[x^32-x+1==0,x]`はなお約25.5 sを要した。内部phaseを分けると，32次ではDurand-Kerner候補生成が約0.40 s，disk separation / orderingが約0.05 sである一方，全32根のRouché certificationが約11.05 sを占めていた。旧証明は各disk中心でexact Rational Taylor係数を構成し，巨大Rationalの正規化を全根で繰り返すため，高次数でrepresentation cliffを作っていた。

現在は各候補diskについて，directed BigFloat intervalで中心`c`における`p(c)`と`p'(c)`を外向き評価し，disk半径`r`と`M >= |c|+r`から

```text
|p(c)| + (r^2/2) max_{|z-c|<=r}|p''(z)| < |p'(c)| r
```

を証明するfast Rouché pathを先に試す。`max|p''|`は係数絶対値多項式を`M`でHorner評価して上から包含する。左辺上界と右辺下界はすべてdirected roundingで構成するため，BigFloat候補を真値として信用しない。この十分条件が閉じないclustered caseだけ，O(n^2) Horner translationで得たexact Rational Taylor係数による従来証明へ切り替える。

この変更で内部Rouché時間は32次で約11.05 sから約0.14 s，64次で約20.40 sから約0.47 sまで低下した。Rouchéが軽くなると次の支配項はDurand-Kerner候補生成になったため，候補補正が十分小さくなった後も最低degree回まで反復する旧条件を廃止し，4回以降にtiny correctionを検出した時点で候補生成を終了する。候補は証明値ではないため，早期終了がcorrectnessを弱めることはなく，disk certificationが閉じなければ従来どおりguard precisionを上げて候補生成から再試行する。

またdegree 16を超えるRootでは，algebraic-field arithmetic自体がcandidate degree 16で制限されているにもかかわらず，materialization時に各Rootへgenerator fieldを構築しようとして同じ高次既約性判定を繰り返していた。この無効なfield構築を高次数では省略した。

2026-08-28 GCC Release / LTO-off，最終`--algebraic-root-cliffs 1`の代表値：

| workload | wall time |
| --- | ---: |
| `x^16-2` all-root create | 約0.158 s |
| two-radius degree 16 create | 約0.223 s |
| `solve[x^16+x+1==0,x]` | 約0.087 s |
| `solve[x^32-x+1==0,x]` | 約0.361 s |
| `solve[x^64+x+1==0,x]` | 約1.562 s |
| `solve[x^65+x+1==0,x]` | 約1.521 s |

追加監査では65次以降にも64境界由来の算法的な不連続はなく，80次約2.9 s，96次約4.1 sまで連続的に評価できた。このため`maximumSupportedAlgebraicDegree`と既定`EvaluationLimits::maxAlgebraicDegree`を64から96へ再設定した。96は数学的定義域ではなく，現時点で実測済みのresource-safety boundaryである。

次の候補生成改善としてAberth–Ehrlich法を残す。比較する場合はNewton-polygon multi-radius seedを共通条件とし，通常根・clustered roots・scale-separated rootsで反復回数とwall timeをDurand-Kernerと比較する。採用してもcandidate generationだけを置換し，最終Rouché certification，disk separation，deterministic orderingの契約は維持する。


# 35. 高階symbolic derivativeの直接構成

`D[expr,{x,n}]`を一階微分の反復だけで処理すると，`LambertW`，`polylog`，`exp[q(x)]`で`cases`・積商則・共通指数因子が毎回再展開され，数学的には単純な高階導函数でも式木が急増する。2026-08-29から，direct variableかつ`n<=64`では函数族ごとのexact recurrenceを先に試す。Lambert WはDLMF 4.13.4_1--4.13.4_2の多項式`p_n(W)`，polylogは`theta=xD`とsigned Stirling number，二次`q`の`exp[q]`は`P_(n+1)=P'_n+q'P_n`を使う。

この経路は近似fast pathではなくexact Expr構成である。64を超えた場合や形が一致しない場合は従来の一般`D`へ代替経路し，公開order上限4096と微分意味論は変更しない。今回の作業では長時間benchmarkを行わず，代表式のcompact outputとtargeted compile/smoke testだけを確認した。

# 36. 有理函数積分のHermite reductionとalgebraic-log 代替経路

従来のexact rational integratorは，一次因子と既約二次因子へ分解できる場合には高速でcompactであったが，`x^3+x+1`のようなQ上既約な高次因子，およびその重複冪で停止していた。単純にComplex Rootへ全分解する経路を早期に置くと，`x/(1+x^4)`のようにreverse-chainで`atan[x^2]/2`へ閉じる式まで大きなRoot/Log和へ退行する。そのため高次algebraic pathはelementary探索後の代替経路とする。

重複度の抽出にはfactor engineを正しさの前提とせず，Q[x]上のYun square-free decompositionを使う。`gcd(Q,Q')`から`Q=product f_i^i`をexactに構成し，各square-free `f_i`の高い冪はHermite step

```text
A/f^k = D[B/f^(k-1)] + C/f^(k-1)
B = -A (f')^(-1)/(k-1) mod f
```

で1段ずつ下げる。`(f')^(-1) mod f`は拡張EuclidでQ[x]上exactに求めるため，数値rootは使用しない。最終的なsquare-free部分`P/Q`は，全complex root `r`を既存のcertified `ComplexAlgebraicNumber`で識別し，

```text
P(x)/Q(x) = sum_r P(r)/(Q'(r)(x-r))
```

から`sum_r P(r)/Q'(r) Log[x-r]`を構成する。residue自体もpersistent algebraic-number arithmeticでexactにmaterializeする。近似rootは候補や恒等式判定には使わない。

この代替経路は既存のdegree 12 specialized rational work budget内でのみ実行する。これは数学的定義域境界ではなく，全複素根isolation，algebraic residue materialization，出力項数を無制限化しないための計算量制限付き policyである。代表的なGCC Release確認では`integrate[1/(x^3+x+1),x]`約0.25 s，`1/(x^5+x+1)`約0.56 s，`1/(x^8+x+1)`約4.0 sであった。長時間benchmarkは実施していない。

将来の改善候補はRothstein–Trager / Lazard–Rioboo–Tragerである。現kernelは各rootを明示するためcorrectnessは閉じているが，等しいresidueや共役rootをgroupingして実`Log/atan`へまとめれば，式サイズとalgebraic materialization量を減らせる可能性がある。これは表現・性能改善であり，今回追加したHermite / square-free capabilityを置換する必要はない。
