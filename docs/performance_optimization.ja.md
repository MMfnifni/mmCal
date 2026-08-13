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

代表benchmark（同一GCC Release環境、16 fractional digits）:

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
- fixed-seed certified FFT round-trip invariant

実行例:

```text
mmCal.Benchmarks
mmCal.Benchmarks --full
mmCal.Benchmarks --random-only
mmCal.Benchmarks --benchmark-only
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

優先度は次のように考える。

1. Toom-4 / higher Toom crossover測定
2. さらに巨大ならFFT/NTT integer multiplication
3. Lehmer GCD
4. `log`のbit-burst / AGM backend
5. Gamma / erf等の5000～10000桁横断benchmark
6. allocator/SBOは実測利益が出る場合のみ

採用時にはこの文書へ「なぜ採用したか」「なぜ前案を棄却したか」「どのbenchmarkで判断したか」を追記する。
