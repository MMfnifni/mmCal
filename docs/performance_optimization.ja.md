# mmCal 高速化・算法選定記録

この文書は、v1.5.1で行った高速化について、**何を採用したか、何を比較したが棄却したか、なぜその判断をしたか**を残すための記録である。

速度だけでなく、mmCalの基本契約である

- exact arithmetic
- principal branch / definedness / domain
- directed rounding
- certified enclosure
- formatter round-trip safety

を壊さないことを採用条件とする。

benchmark値はCPU/compiler/allocator/cacheに依存する。以下は開発時のGCC系環境で得た代表値であり、普遍的な性能保証ではない。Visual Studio / MSVCを含む別環境では`mmCal.Benchmarks`で再測定する。

---

## 1. 採用判断の原則

高速化は次の順で判断する。

1. 旧実装と数学的意味論が一致すること
2. fixed-seed random / invariant testで差分がないこと
3. 境界値・extreme valueで正しさを維持すること
4. 実測で十分な利益があること
5. 小サイズを犠牲にする場合はthreshold dispatchで隔離できること
6. コード複雑性に見合わない微小改善は採用しないこと

したがって「理論上漸近的に速い」「有名libraryでも使われる」だけでは採用理由にしない。

---

# 2. BigInt multiplication

## 2.1 Schoolbook — 維持

小さいoperandでは二重loopのschoolbook multiplicationが最も低overheadだったためbase caseとして維持した。

## 2.2 Karatsuba — 採用

初期threshold sweepでは8～16 limbsからKaratsubaへ入れると明確に退行した。一方、32～48 limbs以降で利益が安定した。

v1.5.1既定値:

```text
Karatsuba crossover ≈ 48 limbs
1 limb = 32 bit
```

極端にunbalancedなoperandはschoolbookへ戻す。

代表値:

| operand | schoolbook | adaptive | speedup |
|---:|---:|---:|---:|
| 128 limbs | 約9.5 µs | 約6.9 µs | 約1.4x |
| 256 limbs | 約38 µs | 約21 µs | 約1.8x |
| 512 limbs | 約151 µs | 約63 µs | 約2.4x |

factorialにもproduct tree経由で波及し、代表測定では

```text
10000!  約4.4 ms → 約2.0 ms
20000!  約19.5 ms → 約6.8 ms
40000!  約88 ms → 約23 ms
```

程度まで改善した。

## 2.3 Toom-3 — 採用

Karatsubaよりさらに巨大なbalanced multiplicationではToom-3が有利になった。

v1.5.1の代表threshold:

```text
top-level Toom-3   ≈ 1280 limbs
recursive Toom-3   ≈ 448 limbs
```

512～1024 limbs付近ではToom-3のevaluation/interpolation overheadが勝つcaseがあるため、早過ぎるdispatchは避けた。

代表値ではKaratsuba-only比で

```text
4096 limbs  約1.17x
6144 limbs  約1.3x
```

程度の改善が得られた。

## 2.4 Toom-4 / FFT / NTT — 保留

v1.5.1では未導入。

Toom-3よりさらに巨大な領域では候補になるが、interactive用途でcrossoverが十分に現れるまで複雑化しない。将来は`mmCal.Benchmarks`でToom-4とFFT/NTTのthresholdを測って判断する。

---

# 3. Dedicated square

`x*x`は一般乗算と異なりcross termが対称なので、専用squareを採用した。

- small: symmetric schoolbook square
- large: Karatsuba square

代表値:

```text
512 limbs   一般乗算 約63 µs → square 約38 µs
1024 limbs  一般乗算 約189 µs → square 約111 µs
```

`pow`のrepeated squaringへ直接効く。

## Toom-3 square — 棄却

専用Toom-3 squareも実装して比較したが、現在のthreshold域ではKaratsuba squareより遅かった。

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

MSVC allocatorでは結果が変わる可能性があるため、将来再測定は可能。

---

# 5. Factorial

## Balanced product tree — 維持

既存の`productRange()`はbalanced treeで同程度の大きさのBigInt同士を掛けるため、adaptive multiplicationとの相性が良かった。

追加した軽量最適化:

- machine integer leafをdecimal string経由でparseせず`BigInt::fromUnsigned()`へ直接構築
- 1-limb multiplicationを`multiplySmall()`へ落とす

## Prime-Swing — 棄却

Prime-Swingも実装し、初版で見つかったprime exponent重複計算等も修正した上で比較した。

それでも代表値は

```text
320000!
product tree  約0.59 s
Prime-Swing   約1.5 s
```

で、現BigInt backendでは既存product treeが速かった。

したがって「一般に高級な算法だから」という理由では採用しなかった。

将来、prime処理や乗算backendが変われば再評価対象。

---

# 6. Division

## 6.1 Knuth normalized long division — 維持

小～中サイズでは低overheadで、Burnikel–Zieglerのbase caseとして優秀なため削除しない。

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

`x / 2^k`を一般divisionへ流さず、

```text
quotient  = x >> k
remainder = low k bits
```

とする。

4096 limbs級でms級から数µs級まで短縮した。

---

# 7. GCD

## Euclidean GCD — 維持

Burnikel–Zieglerの高速化をそのまま利用できるため、既存Euclidean `%`を維持した。

## Binary GCD — 棄却

Stein binary GCDを実装・比較したが、現BigIntでは約4～30倍遅いcaseがあった。

shift/subtraction回数と大きなtemporary処理が多く、現在の高速divisionを使うEuclidに勝てなかった。

次候補はLehmer GCD。

---

# 8. Decimal conversion

巨大factorialの計算本体を高速化すると、次にbinary-limb→decimal conversionが支配的になった。

## 8.1 `10^9` chunk — 採用

旧実装の`/10` 1桁ずつを、`/10^9` 9桁ずつへ変更。

`40000!`で約6秒から約0.5秒級まで改善。

parse側も9桁chunk化し、巨大decimal parseを大幅に短縮した。

## 8.2 Divide-and-conquer decimal conversion — 採用

さらに巨大な`10^(9*2^k)`でほぼ半分へ分割して再帰変換する。

`40000!`（約166,714桁）の代表値:

```text
旧 digit-wise        約6.0 s
10^9 chunk           約0.53 s
D&C                  約0.19 s
```

元実装比で約30倍級。

formatter固有overheadはこの規模でも小さく、主costはBigInt→decimal自体だった。

---

# 9. `tryToUint64`

巨大BigIntを一度10進文字列へ変換してから`from_chars`で失敗する経路を修正した。

64bit超がbit lengthで明白なら即時`nullopt`。

100000-bit級で十数msから測定限界近くまで短縮した。

---

# 10. BigFloat extreme exponent gap

旧加算は常にcommon exponentへexact alignmentし、`1 + 2^-5000000`でも巨大left shiftを作った。

v1.5.1では、要求precisionと符号から結果を一意に証明できるcaseだけfast pathへ入れる。

全rounding mode:

```text
NearestEven
TowardPositive
TowardNegative
TowardZero
```

を扱い、曖昧な場合は旧exact alignmentへfallbackする。

代表値:

```text
53-bit, 1 + 2^-5,000,000
約1.4 ms → sub-µs級
```

速度のためにsmall operandを単純破棄しているわけではない。

---

# 11. Pi

## Machin formula — 旧referenceへ降格

v1.5.0のMachin公式は保証構造が単純だったが、10000桁級で約12秒まで伸びた。

## Binary-splitting Chudnovsky — 採用

代表値:

```text
N[Pi,10000]
約12.3 s → 約0.4 s
```

高桁既知値と照合し、certified enclosureの契約を維持した。

---

# 12. `exp` / `E`

## 逐次RealInterval Taylor — 旧referenceへ降格

高桁でinterval object生成と巨大Rational中間値が支配的になった。

## Binary splitting + certified range reduction — 採用

小さいRationalはexact binary splitting、大きい分子・分母はfixed-precision interval binary splittingへ切り替える。

代表値:

```text
N[E,5000]  数秒級 → 約0.1 s級
```

一般Rationalでも深めのrange reductionを行い、中間項の成長を抑える。

---

# 13. `log`

## 逐次atanh級数 — 旧referenceへ降格

`N[log[2],3000]`が約6～7秒まで伸びた。

## Binary splitting + sqrt range reduction — 採用

- `log[2]`等の小係数: exact binary splitting
- 高bit mantissa: certified sqrtを複数回行って1近傍へ縮約し、interval binary splitting

代表値:

```text
N[log[2],3000]                 約6.7 s → 約0.1 s級
N[log[123456789/987654321],5000]  約22 s → 約3 s級
```

さらに高桁ではbit-burst / AGM系が将来候補。

---

# 14. 巨大Radianの三角函数

旧経路では巨大な生Radian Rationalを小区間へ落とさずTaylor評価し、

```text
N[sin[10^6],20]
```

が5秒timeout級になった。

## Certified argument reduction — 採用

`double fmod`は使わず、Piの保証区間から`x/(Pi/2)`の象限integerを一意に証明し、小区間へ縮約する。

point Rationalだけでなく、`10^6 sqrt[2]`のような保証区間入力にも適用する。

変更後は`10^6`～`10^12`級の代表入力が対話的時間へ戻った。

---

# 15. FFT plan cache

exact radix-2 FFTでは同じtransform sizeでbit-reversalとtwiddleを何度も構築するため、Evaluator/session内でplanをcacheする。

cache対象は入力結果ではなく、size依存のplan情報。

process-globalにはせずSymbol/session lifetimeを安全に保つ。

---

# 15.5. precision-aware `N` と certified FFT

## 旧経路 — exact FFT完成後に`N`を適用

従来の`N`は引数を通常評価してから呼ばれていたため、

```text
N[fft[data],16]
```

でもまず巨大なexact Fourier式を構築し、その後に各成分を近似していた。FFT本体が`Expr`のexact multiply/add/Simplifierをbutterflyごとに通るため、近似値しか要らない場合にもsymbolic costを全額支払っていた。

## precision-aware evaluation — 採用

`N`の第1引数を保持し、precisionを先に確定する。FFT dispatch時にprecision contextが存在すれば、`double`ではなく`ComplexInterval`/BigFloat端点で直接transformする。exactな`fft[...]`の経路は変更しない。

代表benchmark（同一GCC Release環境。当時の旧`N`仕様で16 fractional digits）:

```text
32 points   exact ~11.6 ms   certified ~2.7 ms
64 points   exact ~63.5 ms   certified ~6.2 ms
128 points  exact ~327.6 ms  certified ~13.1 ms
```

非2冪のcertified FFTではdirect DFTとBluesteinを比較し、65点ではdirectが約60 ms、127点ではBluesteinが約199 msでdirect約217 msを上回った。現在は96点未満をdirect、それ以上をBluesteinへ送る。これは数学定数ではなく現benchmark環境のpolicy値なので、MSVCでは再測定する。

採用理由:

- exact-first APIを維持したまま近似要求だけを高速化できる
- machine `double`を導入せず任意精度・外向き丸めを維持できる
- `N`のprecision伝播は将来ほかの高cost builtinにも再利用できる
- 非2冪のO(N^2)崖をapproximate pathではBluesteinで回避できる

exact FFT自体のsymbolic expression explosionは別問題であり、この変更では意図的に残している。

---

# 15.6. Array / Matrix Stage 1–2

## flat Array + exact Number backend — 採用

v1.5.2では行列専用のnested containerを増やさず，既存Arrayの`shape + row-major flat storage`を基盤とした。`MatrixView`はArrayをzero-copy参照し，Gaussian / Gauss-Jordan等で書換えが必要な場合だけflat `MatrixBuffer`へ複製する。これはv1.5.2 release時点の設計であり，Unreleasedではpersistent Arrayのphysical storageをimmutable paged backing + stride viewへ置換している。

exact Number行列ではpivot loopからExpr生成とSimplifier呼出しを外し，`Number`を直接累積・消去する。`dot`も全要素がNumberならcellごとの積和を`Number`だけで処理し，symbolicの場合のみExprを構築する。

一般symbolic `det` / `inverse`はLaplace/adjugate展開の仕事量を共有budgetで制限する。三角行列は次数に依存せず対角積へ落とし，疎行列はbudget内なら処理するが，dense高次行列は階乗級の式を作る前に未評価へ戻す。

## precision-aware certified Matrix — 採用

`N[det[A],p]`等はexact結果を完成してから近似せず，FFTと同じ`ApproximationContext`を受けて`ComplexInterval` backendへ直接dispatchする。expression→interval変換，decimalization，guard-digit増加はFFTと共通化した。

`matrixRank`はexact入力ではexact eliminationを優先する。近似backendではmachine epsilonを用いず、intervalが0を含むがexact zeroでもないpivotは`PrecisionInsufficient`としてguard precisionを増やす。full rank等を非零pivotから証明できる場合は確定するが、近似値だけからrank deficiencyを推測しない。

2026-08-13のRelease / LTO off計測例:

| size | exact `dot` | exact `det` | exact `rref` | `N[det,16]` |
|---:|---:|---:|---:|---:|
| 8 | — | 0.283 ms | 0.434 ms | 1.410 ms |
| 12 | — | 1.290 ms | 2.233 ms | 10.144 ms |
| 16 | 0.328 ms | 3.735 ms | 5.799 ms | 22.799 ms |
| 32 | 1.702 ms | — | — | — |
| 64 | 16.553 ms | — | — | — |

この段階のexact eliminationは通常Gaussian/Gauss-Jordanであり，整数/Rationalの中間分数膨張を抑えるBareiss/fraction-free eliminationはStage 3へ分離した。

# 15.7. Bareiss / fraction-free exact Matrix — 採用

Stage 3ではexact実数行列を行ごとの分母LCMで整数行列へliftし，`IntegerMatrixBuffer`上のBareiss eliminationを共通kernelとして追加した。整数行列はそのまま，Rational行列は各行を非零整数倍してから処理する。

- `det`: Bareissのfraction-free forward eliminationで計算し，Rational入力では行scale積を最後に一度だけ戻す。
- `rref`: Bareissでinteger echelon formまで進め，backward phaseだけRational正規化する。
- `matrixRank`: echelonのpivot数だけで決定し，RREF全体を構築しない。
- `inverse`: `B=D A` として `[B|D]` をfraction-free eliminationし，左側をidentityへ戻した右側を `A^-1` とする。
- exact complex: `Q(i)`等へ整数liftする専用環をまだ持たないため，従来`Number` Gaussian/Gauss-Jordanをfallbackとして保持する。

pivotは数値安定性のためではなく中間BigInt growthを抑えるため，候補中でbit lengthが小さい非零値を優先する。Bareissの各除算は`BigInt::divmod`で余り0を検証し，fraction-free invariantが壊れた場合は黙ってtruncationしない。

2026-08-13 Release / LTO off、同一benchmark入力でStage 2 Gaussianと比較:

| size | `det` Gaussian | `det` Bareiss | speedup | `rref` Gauss-Jordan | `rref` Bareiss | speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 0.256 ms | 0.047 ms | 5.4x | 0.447 ms | 0.058 ms | 7.7x |
| 12 | 1.241 ms | 0.109 ms | 11.4x | 2.058 ms | 0.146 ms | 14.1x |
| 16 | 3.570 ms | 0.385 ms | 9.3x | 5.798 ms | 0.375 ms | 15.5x |

`N[det[...],p]` / `N[inverse[...],p]`等はこのexact Bareiss結果を先に作らず，Stage 2で導入したFFT共通のprecision-aware certified Matrix backendへ直接dispatchする。したがってBareiss採用はexact pathの改善であり，`N`の近似経路を後退させない。

# 15.8. LU / Householder QR — 採用

分解処理は`linear_algebra/decomposition.*`へまとめ，`luDecomposition[A]`はrow-pivoted `P A = L U`，`qrDecomposition[A]`はHouseholder reflectorによる`A = Q R`を実装した。`N[...]`ではexact factorを先に構築せず，FFT/Matrixと共通の`ApproximationContext`からcertified `ComplexInterval` backendへ直接dispatchする。

Householderのapproximate kernelでは複数列を一度のrow-major走査で処理するcolumn-block版も試作した。block=1/8/16/32を複数回Release計測したが，8～24次で差は概ね数%以内かつ最速blockが安定せず，BigFloat/interval演算costが支配的だった。このため既定はblock=1相当とし，block kernelとbenchmarkのみ残した。

2026-08-13 Release / LTO offの代表値:

| size | exact `LU` | `N[LU,16]` | `N[QR,16]` |
|---:|---:|---:|---:|
| 8 | 0.290 ms | 1.811 ms | 7.062 ms |
| 12 | 1.315 ms | 5.355 ms | 21.218 ms |
| 16 | 3.548 ms | 10.562 ms | 46.718 ms |

一般exact Householder QRはradical式の膨張が速く，同benchmark系統で2×2が約1.2 ms，3×3が約59 ms，4×4では約18秒かつformatted outputが約677 KBまで増えた。したがって一般exact QRは3×3以下へpolicy制限し，上三角行列の`{I,A}` fast pathだけ任意次数を許す。4次以上の一般用途は`N[qrDecomposition[A],p]`を推奨する。

# 15.9. reduced SVD — 採用

数値SVDは条件数を二乗する`A^H A`を形成せず，Householder bidiagonalizationの後にone-sided Jacobiで列を直交化する。実数・複素数で同じprecision-aware policyを使い，候補factorはreconstructionとU/V orthogonalityを区間監査してから返す。exact SVDは自然に閉じるcaseへ限定する。

2026-08-13 Release / LTO offで`N[svd[A],16]`を複数回計測した代表値:

| size | time |
|---:|---:|
| 4×4 | 約3.8 ms |
| 8×8 | 約17.1 ms |
| 12×12 | 約46.6 ms |
| 16×16 | 約82.3 ms |

この範囲では急激な悪化はなく，概ね三次成長に沿う。支配costはJacobi反復とBigFloat演算であり，bidiagonalization側だけをcache block化しても寄与率が小さいため，SVD専用block policyは現時点で追加しない。QRのblock kernelは独立benchmarkとして残し，将来backendが変わった時に再測定する。

# 15.10. Eigen / complex Schur — 採用

一般固有値問題は`A^H A`等へ変形せず，Complex BigFloat上でHessenberg reduction → implicit shifted QR → complex Schur形へ進む。固有vectorが必要な場合はSchur三角行列からback substitutionし，Schur vectorを掛け戻す。exact pathは上三角/対角とdistinct-root exact Number 2×2を優先する。

反復停止精度と最終certificateを同じ桁へ置くと行列積で誤差余裕を使い切ることがfixed-seed 3×3で判明したため，内部QR停止精度は表示要求より10桁相当厳しくする。元入力の`ComplexInterval`に対する`A Q-Q T`および`A v-λv`を区間監査し，Schur vectorのunitarityも同時に監視する。非正規行列では固有値・固有vectorのcomponentwise enclosureを無条件には主張せず，Schur/eigenpair relationをcertification境界とする。

2026-08-13 Release / LTO off，random decimal Matrix（[-1,1]，小数10桁相当）:

| size | `N[eigenvalues,16]` | `N[eigensystem,16]` |
|---:|---:|---:|
| 4 | 13.4 ms | 14.0 ms |
| 8 | 68.1 ms | 77.0 ms |
| 16 | 364.8 ms | 466.3 ms |
| 32 | 2826 ms | 3658 ms |
| 64 | 19436 ms | >35 s（計測上限） |

# 15.11. large dense Matrix監査

添付のrandom matrix generatorと同じ「[-1,1]，小数10桁」という入力特性を再現する`--matrix-large <op> <size> [digits]`を`mmCal.Benchmarks`へ追加した。算法benchmarkではparser costを除くため，同じ範囲・10桁量子化をexact RationalとしてC++から直接構築する。PythonとRNG列そのものは一致させず，数値分布と桁幅を合わせる。一方CLI負荷はPython generatorと同じ出力形式を使って別測定する。

32/64次のrepresentative timing:

| op | 32×32 | 64×64 |
|---|---:|---:|
| `N[dot,16]` | 167 ms | 1.21 s |
| `N[det,16]` | 294 ms | 2.96 s |
| `N[inverse,16]` | 1.22 s | 8.68 s |
| `N[matrixRank,16]` | 400 ms | 4.16 s |
| `N[solveLinear,16]` | 549 ms | 5.42 s |
| `N[nullSpace,16]` | 401 ms | 4.01 s |
| `N[LU,16]` | 147 ms | 1.41 s |
| `N[QR,16]` | 1.47 s | 12.88 s |
| `N[SVD,16]` | 1.40 s | 11.60 s |
| `N[eigenvalues,16]` | 2.83 s | 19.44 s |

1024×1024では算法より先にrepresentation costが目立つ。C++から直接1,048,576個の10桁Rational Exprを構築したbenchmark processは入力だけで最大RSS約0.69 GB。`transpose`本体は約107 ms，`trace`本体は約60 msだった。Python形式の約14.16 MBテキストをCLIへ渡し`dimensions[...]`だけを評価した測定ではwall約10.9 s，最大RSS約1.99 GBだった。さらに1024次`N[dot,16]`は10秒上限で未完了（最大RSS約0.96 GB），`N[LU,16]`も10秒上限で未完了（最大RSS約1.59 GB）だったため，QR/SVD/Eigenの1024実走はメモリ圧迫を避けて中止した。

Unreleasedの`Expr::Node` typed-node refactor後，同一x86-64 GCC / Release / LTO offで旧variant sourceと新sourceを同じ`--matrix-large transpose 1024 16`へ掛けて再比較した。旧variant版は最大RSS `693312 KiB`（約677.1 MiB），typed-node版は`299668 KiB`（約292.6 MiB）で，約384.4 MiB / **56.8%削減**。単発`transpose` timingは181.7 ms→157.4 msだったが，timingはnoiseを含むため採用根拠はRSS削減と全regression維持を主とする。

第二段階ではpersistent `ArrayExpr`をfixed-size immutable pageへpackedし，shape / offset / stridesをbackingから分離した。最初に試した単一`vector<Rational>`方式は保存時のRSSは下がるものの，transposeでRational/BigIntを100万要素deep copyし，1024×1024 direct-packed transposeが約650～675 msへ退行したため棄却した。採用版は1024要素pageを`shared_ptr<const page>`で共有し，transposeをstride交換だけのview生成にした。

`ArrayBuilder`はpromotionを現在page内だけへ限定する。完成済みpageは不変なので，大規模Arrayの末尾にsymbolic値が出ても全量Generic化しない。1,048,576要素の最後だけ`x`にした専用測定はall-integer版とほぼ同じ約0.30 s / 約69.8 MiBで，最後のpageだけGenericだった。また矩形brace Lowererはleafを単一builderへ直接流し，numeric literalを一度`Expr` nodeへ包んでからpackする二重表現を避けた。

現`--matrix-large transpose 1024 16`はbenchmark fixtureも`ArrayBuilder`からexact Rationalを直接構築し，最大RSS `136576 KiB`（約133.4 MiB），transpose本体約0.059 ms。`trace 1024`は最大RSS約133.5 MiB，trace本体約10.95 msだった。P1 typed-nodeの292.6 MiBからさらに減っているが，fixture construction pathも同時に現実装へ合わせているため，この差はpersistent storage + builderの総合改善であって単一変更A/Bではない。

CLI側では13.63 MBの1024×1024・10桁decimal literalを`dimensions[...]`へ入力した測定がwall約4.55 s，最大RSS `441564 KiB`（約431 MiB）だった。v1.5.2監査の約14.16 MB / 10.9 s / 1.99 GBとは入力textが完全同一ではないため厳密比較ではないが，parse後段のExpr allocation削減が大きく効いていることを示す。

64次値から純粋なO(n^3)を仮定した1024次の粗い外挿でも，`N[LU]`約1.6時間，`N[dot]`約1.4時間，`N[det]`約3.4時間，`N[solveLinear]`約6.2時間，`N[inverse]`約9.9時間，`N[SVD]`約13時間，`N[QR]`約15時間，`N[eigenvalues]`約22時間となる。32→64の実測指数をそのまま延長すると約1～21時間程度へ揺れるため，これらは予測値であって1024実測ではない。cache・allocator・guard precision・反復回数により悪化し得る。

結論として，1024 dense自体はmachine double + BLASの世界では特別巨大な次数ではないが，mmCalのcertified arbitrary-precision dense算法にとっては依然stress領域である。一方，persistent exact Arrayのrepresentation固定費はtyped-node + paged packed backing + direct builderで大きく下がった。approximate SVD/Eigen等には既に連続working bufferがあるため，次はstorage改善後のcost balanceでblock化・threadingを再評価する。

# 16. `mmCal.Benchmarks`

v1.5.1でVisual Studio solutionへ独立Console projectとして追加した。

```text
mmCal
mmCal.Core
mmCal.Tests
mmCal.Benchmarks
```

通常testへbenchmark時間を混ぜず、CoreへProjectReferenceして次を行う。

- fixed-seed random BigInt division invariant
- decimal round-trip
- certified `exp(x)exp(-x)` invariant
- certified `log(x)+log(1/x)` invariant
- multiply / square / divide threshold benchmark
- factorial benchmark
- decimal parse/toString benchmark
- high-precision `Pi/exp/log` benchmark
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
```

threshold変更時は速度だけでなくrandom invariantを先に通す。

---

# 17. v1.5.1で意図的に採用しなかった一覧

| 候補 | 判断 | 理由 |
|---|---|---|
| Prime-Swing factorial | 棄却 | 現product treeより巨大factorialで遅い |
| binary GCD | 棄却 | Euclidean+B/Zより4～30倍遅いcase |
| Karatsuba vector pool | 棄却 | 5～10%程度退行 |
| Karatsuba depth scratch | 棄却 | 同様に管理costが勝る |
| Toom-3 square | 棄却 | Karatsuba squareより遅い |
| 低threshold Toom-3 | 棄却 | 512～1024 limbsでoverheadが勝つ |
| machine `fmod`による巨大trig縮約 | 棄却 | certified semanticsを失う |
| 全体をMachine/double化 | 方針として不採用 | exact-firstの意味論を変える |

---

# 18. 次の候補

`Expr::Node` typed-node化はUnreleasedで採用済み。公開APIを維持した単独refactorとしてinternal regressionとrandom fuzzerを通し，同一環境の1024 MatrixでRSS約56.8%削減を確認した。

次の優先度は次のように考える。

1. paged packed storage導入後のcost balanceでblocked LU / QR等を再測定
2. pure numeric working kernelだけを対象にthreading thresholdを検討
3. parser AST側の巨大brace temporaryを必要に応じて追加監査
4. BigUInt / BigInt SBOは保守性を損なわないstorage abstractionとして独立benchmarkし，複雑さに見合う場合だけ採用
5. Toom-4 / higher Toom crossover，さらに巨大な整数ではFFT/NTT multiplication
6. Lehmer GCD
7. `log`のbit-burst / AGM backend
8. exact FFTのCyclotomic backend

Arrayについては，単一flat packed vectorを棄却し，immutable paged backing + stride viewを採用した。approximate Matrix algorithmは既存の専用連続working bufferを維持し，persistent Array storageと無理に統合しない。BigUInt SBOは今回明示的に見送る。

採用時にはこの文書へ「なぜ採用したか」「なぜ前案を棄却したか」「どのbenchmarkで判断したか」を追記する。
