# mmCal v1.5.2 自作多倍長数値基盤 — 実装解剖

対象は主として次の層である。

```text
BigUInt
  ↓
BigInt
  ↓
Rational
  ↓
RealNumber / Number

BigInt + Rational
  ↓
BigFloat
  ↓
RealInterval / ComplexInterval
  ↓
CertifiedEvaluator
  ↓
DecimalApproximation
```

記号積分，Solver，Simplifier，特殊函数そのものの規則は扱わない。
ただし，それらが多倍長基盤をどう利用しているかを理解するために必要な範囲で，`N[...]`，保証区間，平方根などとの接続は説明する。

> 対象：mmCal v1.5.2
>
> この文書はAPIリファレンスではなく，mmCalが外部多倍長ライブラリを使わず，整数・有理数・任意精度2進浮動小数・保証区間をどのように積み上げているかを，実装を読みたい人向けに解剖する文書である。
> 「多倍長整数とは何を保存しているのか」「丸めはどこで発生するのか」「なぜ`BigFloat`だけでは保証にならないのか」「巨大行列でメモリを食うのはlimbなのかExprなのか」まで扱う。

## この文書の読み方

最初から通読してもよいが，興味別には次の順が読みやすい。

- **多倍長整数を自作したい**：3～15章。
- **任意精度浮動小数を自作したい**：18～26章。
- **誤差保証まで理解したい**：27～37章。
- **速さの理由を知りたい**：6，8，9，23，37，38章。
- **メモリがどこへ消えるか知りたい**：39章と46章。
- **実装で踏みやすい罠を見たい**：47～49章。

数学記号としての「任意精度」と，計算機上の「無限」を混同しないことが大前提である。mmCalの型は固定桁数を意味論上要求しないが，メモリ，`size_t`，`int64_t`指数，演算時間には当然有限の上限がある。

---

# 1. 設計の要点

mmCal の数値基盤は，大きく **exact 系**と**近似・保証付き評価系**に分かれる。

## 1.1 exact 系

```text
BigUInt      符号なし任意長整数の内部実装
BigInt       符号付き任意長整数
Rational     BigInt / BigInt による既約有理数
RealNumber   BigInt または Rational
Number       exact RealNumber または exact ComplexNumber
```

ここでは通常の整数・有理数演算に浮動小数点を使わない。

例えばソース入力

```text
0.1 + 0.2
```

は，概念的には

```text
1/10 + 1/5
```

として処理され，最終的に

```text
3/10
```

になる。

## 1.2 任意精度近似・保証付き評価系

```text
BigFloat              任意精度の有限2進浮動小数
RealInterval           真値を必ず含む実区間
ComplexInterval        実部・虚部をRealIntervalで囲う複素区間
CertifiedEvaluator     Expr全体を保証付き区間として評価
DecimalApproximation   確定した10進表示値と精度metadata
```

重要なのは，`BigFloat` 自体を「真値」とみなしていないことである。

`BigFloat` は指定方向へ丸められる**作業値**であり，mmCal は下端を負方向，上端を正方向へ丸めた `RealInterval` を構築することで，真値の包含を保証する。

したがって mmCal の高精度数値計算は，

```text
高精度の中心値を1個求める
```

という構造ではなく，基本的に

```text
lower <= true value <= upper
```

を維持しながら進む。

---

# 2. 関連ソース

主要ファイルは次の通り。

| 層 | ファイル | 主な責務 |
|---|---|---|
| unsigned integer | `src/numeric/detail/big_uint.hpp/.cpp` | 32bit limb，適応乗算・専用square，Knuth/Burnikel–Ziegler除算，基数変換 |
| signed integer | `src/numeric/big_int.hpp/.cpp` | 符号付きBigInt，符号付き除算 |
| integer algorithms | `src/numeric/integer_algorithms.hpp/.cpp` | gcd/lcm/pow/factorial/integer sqrt/cuberoot |
| rational | `src/numeric/rational.hpp/.cpp` | 既約有理数，演算前約分 |
| exact real tower | `src/numeric/real_number.hpp/.cpp` | BigInt/Rationalの統合 |
| exact number tower | `src/numeric/number.hpp/.cpp` | exact real/complexの統合 |
| binary scale | `src/numeric/detail/binary_scale.hpp/.cpp` | `floor(log2(p/q))` の整数演算による厳密決定 |
| arbitrary precision float | `src/numeric/big_float.hpp/.cpp` | dyadic BigFloat，4種丸め，四則 |
| rounding | `src/numeric/rounding_mode.hpp` | 最近接偶数・方向丸め |
| real interval | `src/approximation/real_interval.hpp/.cpp` | 外向き丸め区間演算 |
| complex interval | `src/approximation/complex_interval.hpp/.cpp` | 複素保証区間 |
| precision | `src/approximation/precision.hpp/.cpp` | 10進桁→安全な2進bit数 |
| approximation context | `src/approximation/approximation_context.hpp/.cpp` | 要求桁・guard桁・作業精度 |
| decimal result | `src/numeric/decimal_approximation.hpp/.cpp` | 10進丸めと保証区間metadata |
| certified evaluation | `src/approximation/certified_evaluator.hpp/.cpp` | Expr全体の区間評価 |

この周辺はv1.5.1～v1.5.2で高速算法・benchmark基盤が増えているため，行数そのものは仕様値として固定しない。

## 2.1 まず押さえる語彙

### limb

多倍長整数を一定幅に分割した1要素をlimbと呼ぶ。mmCalでは1 limb = 32 bitである。

```text
10進の「桁」ではない
2進32bitの「節」
```

したがって1000 decimal digitsの整数が1000 limbを持つわけではない。必要limb数は概ね

```text
ceil(bitLength / 32)
```

で決まる。

### exact

そのオブジェクトが表す値そのものに丸め誤差がないこと。`BigInt`，`Rational`，そして`BigFloat`の**保存されたdyadic値そのもの**はexactである。

### approximate

求めたい数学的真値と保存値が一致するとは限らないこと。例えば`BigFloat`で丸めた`1/3`はapproximateである。

### certified

真値を1個の近似値として信じるのではなく，

```text
lower <= true value <= upper
```

を演算全体で保証すること。mmCalでは`RealInterval` / `ComplexInterval`がこの責務を持つ。

この区別は重要である。`BigFloat`は「任意精度だから自動的に厳密」ではない。`BigFloat`が厳密なのは**保存している2進有理数そのもの**であって，本来求めたい`Pi`や`1/3`の真値ではない。

---

# 3. `BigUInt`: 最下層の符号なし多倍長整数

## 3.1 limb構造

`BigUInt` の内部基数は

```text
B = 2^32
```

である。

型は

```cpp
using limb_type = std::uint32_t;
using double_limb_type = std::uint64_t;
```

で，内部は

```cpp
std::vector<limb_type> limbs_;
```

を持つ。

limb は**下位から並ぶlittle-endian**で，

```text
limbs_[0] = 最下位32bit
limbs_[1] = 次の32bit
...
```

となる。

整数値 `x` は

```text
x = limbs[0]
  + limbs[1] * B
  + limbs[2] * B^2
  + ...
```

として表される。

例えば

```text
0x00000002_12345678
```

であれば概念的に

```text
limbs_[0] = 0x12345678
limbs_[1] = 0x00000002
```

となる。

## 3.2 zeroのcanonical representation

0 は

```text
limbs_.empty()
```

で表す。

```cpp
bool BigUInt::isZero() const noexcept {
    return limbs_.empty();
}
```

上位に0 limbを残さないことが重要なinvariantであり，演算後には `normalize()` が

```text
末尾の0 limbをpop_back
```

してcanonical化する。

したがって同じ整数に複数の内部表現を作らない。

## 3.3 bit length

最上位limbだけ `std::bit_width` を使う。

```text
bitLength
= (limbCount - 1) * 32
  + bit_width(highest_limb)
```

0のbit lengthは0。

この値は整数平方根，BigFloat変換，binary scale判定など多くの上位算法で使われる。

## 3.5 limbを手で追う

32bit limbを説明のため8bit limbへ縮めて考える。基数を`B=256`とすると，

```text
0x02_34
```

はlittle-endian limb列では

```text
{0x34, 0x02}
```

であり，値は

```text
0x34 + 0x02 * 256 = 564
```

となる。

mmCalの実装ではこれを8bitではなく32bitで行うだけである。内部基数が2の冪なので，bit shift，bit length，trailing zeroの計算が自然にlimb演算へ落ちる。

### なぜdecimal chunkで保存しないのか

10進入出力だけなら`10^9`を基数にした配列も考えられる。しかし内部基数`2^32`には，

- 32bit×32bitが64bitへちょうど収まる。
- shiftが自然。
- `std::bit_width` / `std::countr_zero`が直接使える。
- BigFloatの2進scaleと親和性が高い。

という利点がある。その代わり，人間向け10進表示にはbase conversionが必要になる。v1.5.1で10進divide-and-conquerを入れた理由はここにある。

## 3.4 trailing zero bits

下位limbから0を飛ばし，最初の非0 limbに `std::countr_zero` を使う。

これは後述する `BigFloat` の

```text
significand * 2^exponent
```

をcanonical化するために重要である。

---

# 4. `BigUInt` の加算

加算は通常のlimb単位carry propagation。

各位置について

```text
sum = lhs[i] + rhs[i] + carry
result[i] = sum mod 2^32
carry = floor(sum / 2^32)
```

を行う。

中間値は `uint64_t` なので，

```text
(2^32-1) + (2^32-1) + 1
= 2^33 - 1
```

を安全に保持できる。

計算量はlimb数を `n,m` とすると

```text
O(max(n,m))
```

である。

### self-addition対策

```cpp
x += x
```

では `resize()` によって右辺参照が無効化され得るため，実装は自己参照を検出して一度copyする。

---

# 5. `BigUInt` の減算

`BigUInt` は符号なしなので

```text
lhs < rhs
```

の場合は `underflow_error`。

各limbでborrowを伝播する。

```text
if lhs_i >= rhs_i + borrow:
    out_i = lhs_i - rhs_i - borrow
    borrow = 0
else:
    out_i = B + lhs_i - rhs_i - borrow
    borrow = 1
```

演算後に `normalize()` を行う。

```text
x - x
```

は空vector，すなわちcanonical zeroへ直接落とす。

計算量は

```text
O(n)
```

である。

---

# 6. `BigUInt` の乗算とsquare

v1.5.1で導入した適応dispatchをv1.5.2でも維持し，operand sizeと形状に応じてbackendを切り替える。

```text
small / unbalanced
    ↓
schoolbook
    ↓ 48 limbs前後
Karatsuba
    ↓ 1280 limbs前後（top-level）
Toom-3
```

1 limbは32bit。thresholdは数学的定数ではなく，GCC環境でのmicrobenchmarkから選んだ既定値であり，CPU/compiler/allocatorが変われば`mmCal.Benchmarks`で再測定する。

## 6.1 schoolbook

小さいoperandでは従来の二重loopが最速である。

```text
for i in lhs limbs:
    carry = 0
    for j in rhs limbs:
        k = i + j
        t = lhs[i] * rhs[j] + result[k] + carry
        result[k] = low32(t)
        carry = high32(t)
```

32bit limb × 32bit limbを64bit accumulatorで受けるため，

```text
(B-1)^2 + (B-1) + (B-1) = 2^64-1
```

まで`uint64_t`に収まる。

## 6.2 Karatsuba

十分大きく，かつ左右のサイズ差が極端でない場合に3回の再帰乗算へ分解する。

```text
x = x1 B^m + x0
y = y1 B^m + y0

z0 = x0 y0
z2 = x1 y1
z1 = (x0+x1)(y0+y1)-z0-z2

xy = z2 B^(2m) + z1 B^m + z0
```

実測では8～16 limbsから早期にKaratsubaへ入れるとoverheadで遅く，利益が安定するのはおおむね32～48 limbs以降だった。そのためv1.5.1の既定crossoverは48 limbs付近とした。

極端にunbalancedな積は分割効率が悪いためschoolbookへ戻す。

## 6.3 Toom-3

さらに巨大なbalanced operandは3分割し，`0, 1, -1, 2, infinity`の5点評価・補間で5回の再帰乗算へ落とす。

Top-levelでは小サイズでevaluation/interpolation overheadが勝つため，v1.5.1の既定thresholdは約1280 limbs。Toom再帰内部では既に分割overheadを払っているため，より低い約448 limbsを再帰thresholdとして使う。

代表測定ではKaratsuba-only比で4096 limbs級が約1.17倍，6144 limbs級が約1.3倍高速だった。

## 6.4 専用square

`x*x`は一般乗算へ流さない。

小サイズでは対称性を使い，対角項と上三角cross termだけを計算する。大サイズでは専用Karatsuba squareへ進む。

一般乗算と比較して512 limbsで約1.7倍，1024 limbsで約1.7倍程度の改善が得られ，`pow`のrepeated squaringにもそのまま波及する。

Toom-3専用squareも実装・比較したが，現在のthreshold域では通常のKaratsuba squareより遅かったため既定経路には採用していない。

## 6.5 採用しなかったworkspace化

Karatsuba再帰のtemporary `vector` allocation削減を狙い，

- pool型workspace
- 再帰depthごとのscratch型workspace

を比較した。しかしGCC環境では512～1024 limbs付近で最大約5～10%退行した。管理・resize・cache localityのcostがallocator削減を上回ったため，v1.5.1では採用しない。

これは将来MSVCやallocator特性が変わった場合の再測定候補である。

# 7. `BigUInt` のbit shift

`<<` / `>>` は

```text
32bit単位のlimb shift
+
limb内部のbit shift
```

に分解する。

例えば `bits` に対し

```text
limbShift = bits / 32
bitShift  = bits % 32
```

とする。

左shiftでは64bit中間値で上位carryを次のlimbへ流す。

右shiftでは

```text
current >> bitShift
```

と，一つ上のlimbから流れ込む

```text
high << (32-bitShift)
```

をORする。

右shiftで捨てられた下位bitは復元しない。

---

# 8. `BigUInt` の除算

v1.5.1で導入した**特殊case → Knuth base case → Burnikel–Ziegler**の段階dispatchをv1.5.2でも使う。

## 8.1 fast path

先に次を処理する。

```text
divisor == 0        → domain_error
dividend < divisor  → quotient=0, remainder=dividend
dividend == divisor → quotient=1, remainder=0
```

除数が1 limbなら`divideSmall()`。

さらに除数がexactな`2^k`なら一般長除算へ入れず，

```text
quotient  = dividend >> k
remainder = dividend の下位 k bit
```

で処理する。4096 limbs級では旧一般除算のms級から数µs級まで短縮された。

## 8.2 Knuth normalized long division

小～中サイズ，商が小さい場合，Burnikel–Ziegler再帰のbase caseには従来のKnuth型normalized long divisionを残す。

内部基数は`B=2^32`。

1. 除数最上位limbをleading-zero shiftで正規化
2. 上位2 limbから商digit `qhat` を推定
3. 次limbで過大推定を補正
4. `qhat * divisor` を減算し，必要なら1回add-back
5. remainderをde-normalize

この旧実装を捨てなかった理由は，小さいoperandでは再帰分割より定数costが小さく，Burnikel–Zieglerの良いbase caseになるためである。

## 8.3 Burnikel–Ziegler

巨大でbalancedなdivisionでは被除数・除数をblockへ分割し，2n/1n・3n/2n型の再帰divisionへ落とす。

v1.5.1のGCC benchmarkでは32 limbs前後から利益が安定したため，巨大balanced divisionへ適用する。商が小さいcaseはKnuthへ残す。

代表測定:

```text
1024-limb divisor/quotient: 約1.18 ms → 約0.12 ms
2048-limb divisor/quotient: 約5.0  ms → 約0.35 ms
```

これにより`%`, GCD, Rational正規化, `integerSqrt`, `integerCubeRoot`, divide-and-conquer decimal conversionにも波及する。

## 8.4 不変条件

どのbackendでも返す結果は必ず

```text
q*d + r == dividend
0 <= r < divisor
```

を満たす。

random testでは境界・不均衡サイズを含む巨大caseを固定seedで生成し，再構築不変条件を確認する。

# 9. `BigUInt` の文字列変換

## 9.1 parse

基数2..36に対応する。通常基数では

```text
value = value * radix + digit
```

をBigUInt上で進める。

10進についてはv1.5.1で`10^9` chunkを使い，9桁をまとめて取り込む。巨大10進文字列を1桁ずつ処理する旧経路より走査回数を大きく減らす。

`tryToUint64`も，明らかに64bitを超えるBigIntを一度10進文字列化して`from_chars`へ渡す旧経路を廃止し，bit lengthで即時棄却するfast pathを持つ。

## 9.2 toString(10)

v1.5.0初期の1桁ずつ`/10`する方式から，まず`10^9` chunkへ変更した。さらにv1.5.1では巨大値に**divide-and-conquer base conversion**を使う。

概念的には大きな

```text
10^(9 * 2^k)
```

をsquareで構築し，

```text
value = high * P + low
```

となるよう`divmod(value,P)`でほぼ半分へ分けて再帰変換する。小さい葉では`10^9` chunk変換へ戻る。

`40000!`（約166,714 decimal digits）の代表測定では，

```text
旧 /10 digit-wise          約6.0 s
10^9 chunk                 約0.53 s
divide-and-conquer         約0.19 s
```

まで短縮された。

formatterへBigIntを包む追加costはこの規模でもごく小さく，巨大整数表示の主costはほぼbinary-limb→decimal変換そのものだった。

## 9.3 他基数

2進・16進等はradixの性質に応じた既存経路を維持する。10進D&Cは表示上最も頻繁で，かつ旧実装のbottleneckが顕著だったため専用最適化としている。

# 10. `BigInt`: 符号付き任意精度整数

`BigInt` は2の補数多倍長ではなく，**sign-magnitude**である。

```cpp
bool negative_ = false;
detail::BigUInt magnitude_;
```

したがって値は

```text
negative ? -magnitude : magnitude
```

として保持する。

## 10.1 negative zeroを持たない

`normalizeSign()` により

```text
magnitude == 0
```

なら必ず

```text
negative = false
```

にする。

```text
-0
```

という別表現を内部に残さない。

## 10.2 `INT64_MIN`

`-INT64_MIN` は符号付き64bitでoverflowするため，constructorでは

```text
-(value + 1) + 1
```

の形でabsolute magnitudeを作る。

## 10.3 加算

符号が同じならmagnitude加算。

符号が異なるならmagnitudeを比較し，大きい側から小さい側を引き，結果符号を大きい側へ合わせる。

```text
(+a)+(+b) → +(a+b)
(-a)+(-b) → -(a+b)
(+a)+(-b) → magnitude比較後に差
```

## 10.4 乗算

magnitudeは `BigUInt` 乗算。

符号はXOR。

```text
negative = lhs.negative != rhs.negative
```

0になれば `normalizeSign()` でpositive zeroへ戻す。

## 10.5 除算規則

まずabsolute magnitude同士を `BigUInt::divmod` し，後で符号を付ける。

商:

```text
sign(dividend) XOR sign(divisor)
```

余り:

```text
sign(dividend)
```

したがってC++整数除算と同じく**0方向truncate**である。

例:

```text
 7 /  3 =  2, remainder  1
-7 /  3 = -2, remainder -1
 7 / -3 = -2, remainder  1
-7 / -3 =  2, remainder -1
```

常に

```text
q*d + r == dividend
|r| < |d|
```

を保つ。

## 10.6 負数の右shift

sign-magnitudeなので，

```text
-7 >> 1
```

は2の補数の算術shiftではない。

実装はabsolute magnitudeを右shiftして符号を戻すので，意味は**0方向切り捨て**である。

この仕様はbitwise integer型として使う場合に重要である。

## 10.7 符号のためにlimbを1個余計に使っているか

使っていない。`BigInt`は2の補数limb列ではなく，

```text
negative_ + magnitude_
```

というsign-magnitudeである。負数だから最上位へ`0xFFFFFFFF`を延々と並べるような符号拡張はない。

例えば`-1`もmagnitude側は正の`1`と同じ1 limbで，符号は`bool negative_`へ分離される。

ただしC++ objectとしてはalignment/paddingがあるため，「符号は1bitだからコストも1bit」という意味ではない。x86-64 GCCでのv1.5.2参考測定では，

```text
sizeof(BigUInt) = 24 bytes
sizeof(BigInt)  = 32 bytes
```

だった。これはABI依存の参考値であり仕様ではない。重要なのは，余分な**数値limb**を符号用に保持してはいないという点である。

---

# 11. 整数算法

`integer_algorithms.cpp` にはBigInt上の共通算法がまとまっている。

## 11.1 GCD

標準Euclid算法。

```text
a = abs(a)
b = abs(b)
while b != 0:
    a %= b
    swap(a,b)
return a
```

現在はbinary GCDではなく，`BigInt % BigInt`，つまり多倍長除算を使う。

## 11.2 LCM

```text
lcm(a,b) = abs((a/gcd(a,b))*b)
```

先にGCDで割ってから掛けるため，不要に巨大な中間積を作りにくい。

## 11.3 整数累乗

exponentiation by squaring。

```text
while exponent != 0:
    if exponent odd:
        result *= base
    exponent >>= 1
    if exponent != 0:
        base *= base
```

乗算回数は指数値そのものではなく概ね `O(log exponent)` 回。

## 11.4 factorial

既定実装は**balanced product tree**。

単純な

```text
1*2*3*...*n
```

という左結合にはせず，`productRange(first,last)`を中央で二分して同程度の大きさ同士を掛ける。

小区間はstraight loopとし，葉ではmachine整数を文字列化して再parseする旧経路を避け，`BigInt::fromUnsigned()`で直接構築する。また1-limb operandは一般multi-limb multiplicationより`multiplySmall()`を優先する。

Karatsuba/Toom-3と専用squareの導入後，product treeは巨大factorialでも大きく改善した。

### Prime-Swingを採用しなかった理由

Prime-Swing factorialも実装し，sieve，odd factorial，2の冪分離，balanced productまで比較した。しかし現在のBigInt backendでは既存product treeの方が速かった。

代表例:

```text
320000!:
product tree  約0.59 s
Prime-Swing   約1.5 s
```

初版Prime-Swingにはprime exponentの重複計算等の無駄があったためそこも最適化したが，最終的にも逆転しなかった。したがってv1.5.1では「高度そうだから」という理由だけで置換せず，実測で勝つproduct treeを既定とした。

Prime-Swing自体を一般に否定するものではなく，乗算backendやprime処理が変われば再評価可能である。

## 11.5 integer square root

非負BigIntに対して

```text
root = floor(sqrt(n))
remainder = n - root^2
```

を返す。

### 初期値

`bitLength(n)` から

```text
initialBit = (bitLength + 1) / 2
guess = 2^initialBit
```

とし，真のsqrtより上側から開始する。

### Newton反復

```text
next = (guess + n/guess) / 2
```

を繰り返し，

```text
next >= guess
```

になったところで停止する。

上側から単調に収束させるため，最終 `guess` が `floor(sqrt(n))` になる。

最後に

```text
remainder = n - guess*guess
```

を返す。

## 11.6 integer cube root

同様に

```text
root = floor(cuberoot(n))
```

をNewton法で求める。

初期値:

```text
2^ceil(bitLength/3)
```

実装上は

```text
initialBit = (bitLength + 2) / 3
```

Newton反復:

```text
next = (2*guess + n/(guess^2)) / 3
```

停止後，整数Newtonの停止点を厳密なfloorへ合わせるため

```text
while guess^3 > n:
    guess--
while (guess+1)^3 <= n:
    guess++
```

で補正する。

---

# 12. `Rational`: 任意精度有理数

内部は

```cpp
BigInt numerator_;
BigInt denominator_{1};
```

である。

常に次のcanonical invariantを保つ。

```text
denominator > 0
gcd(abs(numerator), denominator) = 1
0 は 0/1
```

## 12.1 decimal literalのexact parse

例えば

```text
12.500
```

は文字列から

```text
digits = 12500
fraction length = 3
```

として

```text
12500 / 10^3
```

を構築し，その後約分して

```text
25/2
```

にする。

基数指定も同じで，例えばbase 2なら

```text
101.01₂ = 10101₂ / 2^2 = 21/4
```

となる。

## 12.2 constructor normalize

一般の `Rational(a,b)` は

1. denominator=0を拒否
2. denominatorが負なら分子・分母双方を反転
3. zeroなら `0/1`
4. `gcd(abs(a),b)` で約分

を行う。

---

# 13. Rational加算の最適化

単純式

```text
a/b + c/d = (ad+bc)/(bd)
```

をそのまま使うと，巨大な `bd` を作った後で巨大GCDを計算することになる。

現行実装は，入力が既に既約であることを利用する。

```text
g = gcd(b,d)
lhsScale = d/g
rhsScale = b/g
n = a*lhsScale + c*rhsScale
r = gcd(abs(n), g)
n /= r
den = (b/r) * lhsScale
```

つまり最終分子と共通因子になり得る範囲を `g=gcd(b,d)` に限定し，巨大な最終分母全体とのGCD再計算を避ける。

これは多倍長演算では重要な最適化である。

---

# 14. Rational乗算の交差約分

```text
a/b * c/d
```

に対して積を作る前に

```text
leftCancel  = gcd(abs(a), d)
rightCancel = gcd(abs(c), b)
```

を求め，

```text
numerator
= (a/leftCancel) * (c/rightCancel)

denominator
= (b/rightCancel) * (d/leftCancel)
```

とする。

入力Rationalが既約なので，この2方向のcross-cancel後は結果も既約となり，巨大積を作った後の追加GCDを不要にしている。

---

# 15. Rational除算の事前約分

```text
a/b ÷ c/d = ad/bc
```

だが，先に

```text
numeratorCancel   = gcd(abs(a), abs(c))
denominatorCancel = gcd(b,d)
```

を取る。

その後

```text
numerator
= (a/numeratorCancel) * (d/denominatorCancel)

denominator
= (b/denominatorCancel) * (c/numeratorCancel)
```

とする。

これも中間BigIntサイズを抑えるための実装である。

---

# 16. `RealNumber` と exact numeric tower

`RealNumber` は

```cpp
std::variant<BigInt, Rational>
```

である。

Rationalの分母が1になった場合は `normalize()` により自動的にBigIntへ縮約する。

```text
6/3
```

をRational `2/1` として保持し続けず，

```text
BigInt(2)
```

へ戻す。

整数同士の加減乗はBigInt fast pathを使い，Rationalが必要になった時だけRationalへpromoteする。

例えば整数同士の除算は一旦

```text
Rational(lhs,rhs)
```

を作り，整数に閉じれば再度BigIntへnormalizeされる。

---

# 17. `Number`: exact real / complex

`Number` は

```text
RealNumber
or
ComplexNumber{RealNumber real, RealNumber imaginary}
```

を保持する。

つまり複素数の両成分も任意精度整数・有理数である。

```text
(3+4I)/5
```

なら内部的にもexactな

```text
3/5 + 4/5 I
```

を維持できる。

虚部がexact zeroなら `normalize()` によりRealNumberへ縮約する。

この層には `double` は入らない。

---

# 18. `BigFloat`: 任意精度2進浮動小数

`BigFloat` はIEEE floating pointをそのまま多倍長化した型ではない。

内部値は常に

```text
significand * 2^exponent
```

である。

```cpp
BigInt significand_;
std::int64_t exponent_;
std::size_t precisionBits_;
```

## 18.1 exponentは多倍長ではない

仮数はBigIntなので任意長だが，指数は

```text
int64_t
```

である。

したがって「任意精度」は仮数精度についてであり，指数範囲まで無限ではない。

指数overflow/underflowは明示的に検査し，例外にする。

## 18.2 finite only

現在のBigFloatは

- NaN
- +Infinity
- -Infinity

を内部値として持たない。

数学的domain errorやInfinityという数学記号は上位層で扱う。

---

# 19. BigFloatのcanonical dyadic representation

0以外ではsignificandから2の因子を全て除く。

例:

```text
12 * 2^-3
```

は

```text
12 = 3 * 2^2
```

なので

```text
3 * 2^-1
```

へnormalizeする。

実装は

```text
zeros = significand.trailingZeroBits()
significand >>= zeros
exponent += zeros
```

である。

これにより同一dyadic rationalに複数のrepresentationができにくい。

注意点として，`precisionBits_` は常にsignificandの実bit長そのものではない。

これは**その値を生成したときの目標有効bit数**として保持されるmetadataであり，normalizeでsignificandの末尾0が消えれば実bit長は小さくなり得る。

## 19.1 `precisionBits`は「そのbit数だけメモリを予約する」という意味ではない

ここは誤解しやすい。

```text
N[x,10000]
```

相当の作業で`precisionBits_`が大きくなっても，`BigFloat`が常にそのbit数ぶんの0-filled limb領域を抱えるわけではない。実際のpayloadは`BigInt significand_`であり，canonicalize後に必要なlimbだけを持つ。

例えばexactな`1`は高いprecision metadataを持っていても，significandそのものは`1`で済む。

逆に，非dyadicなRationalをp bitへ丸めれば，通常はp bit前後のsignificandが必要になる。

したがってBigFloatのメモリ量は，概ね

```text
object fixed cost + significand limb capacity
```

で決まり，`precisionBits_`の数値だけから一意には決まらない。

---

# 20. BigFloatの丸めmode

4種類ある。

```text
NearestEven
TowardZero
TowardPositive
TowardNegative
```

特に `RealInterval` では

```text
lower → TowardNegative
upper → TowardPositive
```

を使う。

## 20.1 exact remainderによる丸め判定

BigFloat変換では，machine floating pointのrounding flag等に依存しない。

整数除算から

```text
quotient
remainder
divisor
```

をexact BigIntとして持ち，そこから丸め方向を決める。

### TowardZero

余りがあってもincrementしない。

### TowardPositive

正数で余りがあればincrement。

### TowardNegative

負数で余りがあればmagnitudeをincrementし，より負側へ送る。

### NearestEven

```text
2*remainder ? divisor
```

をexact比較する。

```text
2r > d → increment
2r < d → keep
2r = d → tie
```

完全なtieなら保持側quotientの最下位bitを見て，偶数側を選ぶ。

したがって丸め判定そのものに `double` の誤差は入らない。

---

# 21. Rational → BigFloat変換

`BigFloat::fromRational()` は内部的に `fromPositiveRatio()` を使う。

正の比

```text
n/d
```

に対してまず

```text
e = floor(log2(n/d))
```

を厳密に求める。

ここで `log2()` の浮動小数函数は使わない。

## 21.1 `floorLog2PositiveRatio`

まず

```text
bitLength(n)
bitLength(d)
```

の差から候補指数を作る。

例えば `n` のbit長が大きい場合，

```text
k = bitLength(n)-bitLength(d)
```

として

```text
n ? d*2^k
```

をBigInt比較し，必要なら `k-1` に補正する。

これで

```text
2^e <= n/d < 2^(e+1)
```

を整数演算だけで決定する。

## 21.2 p-bit仮数の生成

要求precisionを `p` とする。

```text
m = (n/d) * 2^(p-1-e)
```

になるよう，nまたはdをbit shiftしてscaleする。

その後exact `divmod` を行い，

```text
quotient = floor(m)
remainder = exact remainder
```

を得る。

このremainderで前節の丸めを行う。

最終指数は

```text
outputExponent = exponentOffset + e - (p-1)
```

となる。

最後に `fromDyadic()` へ通してcanonicalizeする。

---

# 22. exact dyadic → BigFloat

入力が

```text
S * 2^e
```

というexact dyadicなら `fromDyadic()`。

`abs(S)` のbit長が要求precision以下なら値を一切丸めない。

超えていれば

```text
discardedBits = bitLength(S) - precisionBits
q = abs(S) >> discardedBits
remainder = abs(S) - (q << discardedBits)
divisor = 2^discardedBits
```

として丸めを判定する。

丸め後は捨てたbit数だけexponentへ移す。

```text
newExponent = e + discardedBits
```

---

# 23. BigFloat加算

通常caseでは指数を揃え，exact dyadic和を作ってから指定precisionへ丸める。

```text
commonExponent = min(lhs.exponent, rhs.exponent)
L' = L << (e1-common)
R' = R << (e2-common)
resultExactSignificand = L' + R'
```

この旧経路は近接値・cancellation・丸め境界を確実に扱えるため現在もbase pathとして残す。

## 23.1 extreme exponent-gap fast path

v1.5.1では，一方が要求precisionに比べ極端に小さく，かつ符号・距離からcancellationや丸め境界への影響を安全に判定できる場合だけ，巨大left shiftを作らず結果を直接決める。

重要なのは「小さいから捨てる」のではなく，丸めmodeごとに

```text
NearestEven
TowardPositive
TowardNegative
TowardZero
```

でdominant値そのもの，または直上/直下のrepresentable valueのどれになるかを証明して返すことである。

曖昧なcaseは必ず従来exact alignmentへfallbackする。

代表測定では53bit precisionの`1 + 2^-5,000,000`級が約1.4 msからsub-µs級へ短縮され，全4 rounding modeを旧実装と差分照合している。

# 24. BigFloat乗算

```text
(L*2^a) * (R*2^b)
= (L*R) * 2^(a+b)
```

なので，

```text
significand = lhs.significand * rhs.significand
exponent = lhs.exponent + rhs.exponent
```

をexactに作り，最後に `fromDyadic()` で指定precisionへ丸める。

仮数乗算の速度は基礎`BigInt`の適応schoolbook/Karatsuba/Toom-3 backendの影響を直接受ける。

---

# 25. BigFloat除算

```text
(L*2^a)/(R*2^b)
= (L/R) * 2^(a-b)
```

として，absolute significand ratioを `fromPositiveRatio()` へ渡す。

そこでexact BigInt `divmod` とremainderを使って指定方向へ丸める。

したがってBigFloat除算はmachine floating divisionに依存しない。

---

# 26. BigFloat比較

比較でも必要ならsignificandをscaleしてexact比較する。

ただし無条件に巨大shiftを作るわけではない。

指数gapが大きく，bitLengthだけで大小が確定する場合はその時点で返す。

例えば `lhs.exponent > rhs.exponent` の場合，gapがrhs magnitudeのbit長以上ならlhsの方が大きいことをshiftなしで確定できる。

このearly-outは，加算より比較の方が極端なscale差に強い理由である。

---

# 27. `RealInterval`: BigFloatを保証付き計算へ使う

`RealInterval` は

```text
[lower, upper]
```

という閉区間。

端点は `BigFloat` だが，各BigFloat自体はexactなdyadic rationalなので

```text
lower.toRational()
upper.toRational()
```

でexactな有理数へ戻せる。

invariantは

```text
lower <= true value <= upper
```

である。

## 27.1 Rationalからの区間化

exact Rational `x` をp bitへ変換するとき

```text
lower = BigFloat::fromRational(x,p,TowardNegative)
upper = BigFloat::fromRational(x,p,TowardPositive)
```

とする。

`x` がdyadic rationalなら双方が一致してpoint intervalになる。

例:

```text
1/2 → exact dyadic → [1/2,1/2]
1/3 → non-dyadic  → [lower,upper]
```

## 27.2 point intervalですら実装を雑にすると壊れる

`RealInterval::point(x)`は数学的には単に

```text
[x,x]
```

である。しかしC++実装ではv1.5.2開発中に，

```cpp
return RealInterval{value, std::move(value)};
```

という形がMSVCで問題を露出させた。同一objectを同一初期化式の中でcopyとmoveの双方へ使うと，引数評価順やmove後状態に依存し得る。負値で一方がmove後の0相当になれば，概念的に

```text
[0, negative]
```

となり，`lower > upper` invariantを破壊する。

現在は明示的に2端点へcopyしてからmoveする。

```cpp
BigFloat lower = value;
BigFloat upper = value;
return RealInterval{std::move(lower), std::move(upper)};
```

この一件は，多倍長や区間演算で怖いのは数学式だけではなく，**C++ object lifetimeと評価順も証明の一部**であることを示す好例である。

---

# 28. RealInterval四則

## 28.1 加算

```text
[a,b] + [c,d]
= [a+c, b+d]
```

ただし

```text
lower = add(a,c,TowardNegative)
upper = add(b,d,TowardPositive)
```

## 28.2 減算

```text
[a,b] - [c,d]
= [a-d, b-c]
```

同様に外向き丸め。

## 28.3 乗算

4隅

```text
ac, ad, bc, bd
```

をすべて計算する。

下端候補は全て負方向丸め，上端候補は全て正方向丸めし，

```text
lower = min(lower products)
upper = max(upper products)
```

とする。

## 28.4 除算

分母区間が0を含む場合は拒否。

0を跨がない場合は4隅の商を計算し，同様にmin/maxを取る。

---

# 29. ComplexInterval

複素近似も1個のcomplex BigFloatではなく，

```text
real:      RealInterval
imaginary: RealInterval
```

の直積として持つ。

虚部が「小さい」だけではrealへ縮約しない。

```text
imaginary == exact [0,0]
```

の場合だけ `isProvablyReal()` がtrueになる。

これはmmCalの「近いから同じとみなさない」というexact-first方針に一致する。

---

# 30. 10進要求桁から2進作業精度への変換

`N[...,n]` の `n` は10進小数部桁数。

内部 `BigFloat` は2進precisionなので変換が必要。

現在は

```text
log2(10) ≈ 3.321928...
```

に対して安全側の

```text
3.322 = 3322/1000
```

を使う。

つまり

```text
binaryBits = ceil(decimalDigits * 3322 / 1000)
```

相当をinteger arithmeticで計算する。

実際の作業桁は

```text
requested decimal digits + guard digits
```

である。

既定guardは8桁。

---

# 31. guard precisionの増加

`ApproximationContext` は

```text
decimalDigits
guardDigits
roundingMode
```

を持つ。

初期guardは8。

保証区間の幅がまだ広く，要求10進表示が一意に決まらなければguardを増やして再評価する。

Evaluatorの一般 `N` では

```text
growth = max(8, guard/2)
guard += growth
```

なので概ね

```text
8 → 16 → 24 → 36 → 54 → 81 → ...
```

と増える。

重要なのは，guard桁数そのものを「正しさの証拠」にしていないこと。

最終条件は**保証区間の両端が要求桁へ同じ10進丸め結果を持つこと**である。

---

# 32. `DecimalApproximation`

これは計算途中のfloating valueではなく，ユーザーへ返す**確定済み10進結果**である。

保持内容:

```text
text                      表示文字列
fractionalDigits          実際の小数部桁数
requestedFractionalDigits 要求桁数
rounded                   丸めが行われたか
origin                    ExactValue / CertifiedInterval
displayedValue            表示10進値そのもののexact Rational
certifiedLower            真値を含む下界 Rational
certifiedUpper            真値を含む上界 Rational
```

つまり

```text
"3.141592..."
```

だけを保存しているのではない。

表示値そのものもexact Rationalへ戻せ，さらにその値がどの保証区間から確定したかを保持する。

`accuracy`，`precision`，`rationalize` が文字列を再parseして精度を推測しなくてよいのはこのためである。

---

# 33. exact Rationalの10進化

`DecimalApproximation::fromReal()` では，まず分母から2と5を全て除いて

```text
残り == 1
```

なら有限10進と判定する。

有限なら割り切れるまでlong divisionする。

循環小数なら要求桁+guard digitを生成し，最近接偶数丸めする。

10進digit生成もBigIntだけで

```text
remainder *= 10
digit = remainder / denominator
remainder %= denominator
```

と進む。

---

# 34. 保証区間から10進結果を確定する方法

`fromCertifiedInterval(lower,upper,n)` では，lowerとupperをそれぞれ**同じn桁へnearest-evenで丸める**。

```text
round_n(lower)
round_n(upper)
```

が同じ文字列なら，区間内部の真値も同じ結果へ丸まるため，そのdecimalを確定できる。

異なるなら

```text
std::nullopt
```

を返し，上位Evaluatorが作業precisionを増やしてやり直す。

これがmmCalの `N[...]` の最終certification条件である。

## 34.1 v1.5.2のcompact decimal表示

内部のcertification metadataと，CLIで人間へ見せる末尾0は同じ情報ではない。v1.5.2ではcertified fixed-digit結果について，冗長な末尾0列を圧縮する。

```text
1.000000000000 → 1.0
1.500000000000 → 1.50
1.230000000000 → 1.230
```

ただし要求precisionや`certifiedLower/certifiedUpper`は内部に全て残る。したがって表示が`1.0`になっても，「1桁しか計算していない」わけではない。

一方，exact Rationalの有限小数はもともと必要以上に0埋めしない。

```text
N[1/2,10] → 0.5
```

これは表示policyであり，exact/approximate semanticsをdecimal literalへ逆輸入しない。mmCalでは`0.5`を再入力すればexact `1/2`であり，表示文字列だけを完全なapproximate round-trip syntaxとはみなしていない。

---

# 35. `N[...]` までの流れ

symbolic expressionに対する概略は次の通り。

```text
N[expr, n]
   ↓
ApproximationContext(decimalDigits=n, guard=8)
   ↓
workingBinaryBits = ceil((n+guard)*3.322)
   ↓
CertifiedEvaluator.enclose(expr, workingBinaryBits)
   ↓
RealInterval または ComplexInterval
   ↓
区間端点をexact Rationalへ戻す
   ↓
DecimalApproximation::fromCertifiedInterval(..., n)
   ↓
両端のn桁丸めが一致？
   ├─ yes → DecimalApproximationを返す
   └─ no  → guard増加 → 式全体を再評価
```

一方，入力が最初からexact `Number`，つまり整数・有理数・exact複素数なら，一般CertifiedEvaluatorを通さず直接 `DecimalApproximation::fromReal()` で10進化するfast pathがある。

## 35.1 v1.5.2のprecision-aware `N`

一般のscalar式では上記の`CertifiedEvaluator`が中心になる。ただしv1.5.2では，FFTや線形代数のように「exactな巨大中間式を作ってから近似すると本質的に遅い」builtinについて，`N`が要求precisionを**先に**確定してから専用backendへdispatchできる。

概念的には，

```text
旧: N[fft[data], p]
      ↓
    exact Fourier expressionを構築
      ↓
    最後に近似

現: N[fft[data], p]
      ↓
    pを先に確定
      ↓
    ComplexInterval / BigFloat FFTへ直接dispatch
```

Matrixも同様で，

```text
N[inverse[A], p]
N[qrDecomposition[A], p]
N[svd[A], p]
N[eigenvalues[A], p]
```

は，巨大なexact inverseやradical式を一度完成してから数値化する設計ではない。

これはexact-firstと矛盾しない。

```text
fft[exactData]      → exact path
N[fft[exactData],p] → certified approximate path
```

と，呼出し側が`N`で近似を明示しているからである。

### 不連続量は例外

`matrixRank`や`nullSpace`は微小摂動で結果の次元自体が変わる。このためexact入力ならexact eliminationで構造を確定できる場合を優先し，「近いから0」というepsilon判定はしない。precision-awareとは，何でも近似へ落とすことではない。

---

# 36. 多倍長整数を保証付き平方根へ接続する例

`certified_sqrt.cpp` は，BigInt基盤が上位数値評価へどう使われるか分かりやすい例である。

exact Rational

```text
x = p/q >= 0
```

のsqrtを求める。

## 36.1 2進scale正規化

まずexactに

```text
x = 4^k * z
1 <= z < 4
```

へ持っていく。

`floor(log2(x))` はmachine logを使わず，前述の `floorLog2PositiveRatio()` で決める。

## 36.2 fixed-point整数平方根へ変換

`sqrt(z)` をS bitの固定小数点整数 `m` として

```text
m = floor(sqrt(z) * 2^S)
```

としたい。

`z=p/q` なので

```text
Q = floor(p * 2^(2S) / q)
m = floor(sqrt(Q))
```

と変換できる。

ここで

- `p * 2^(2S)` → BigInt shift
- `/ q` → BigInt divmod
- `sqrt(Q)` → integerSqrt

で，主要部分を全てexact integer arithmeticとして処理する。

下界integerを `m` とし，真値がちょうどdyadic rootかを

```text
m^2 * q == p * 2^(2S)
```

でexact比較する。

exactでなければ上界を `m+1` とする。

最後に両者を指定precisionのBigFloatへ

```text
lower: TowardNegative
upper: TowardPositive
```

で丸め，RealIntervalを作る。

この構造は「多倍長整数を計算の土台にし，BigFloatは保証区間端点へ使う」というmmCalの設計をよく表している。

---

# 37. 高精度超越計算への接続

## 37.1 Pi: binary-splitting Chudnovsky

v1.5.0のMachin公式

```text
Pi = 16 atan(1/5) - 4 atan(1/239)
```

は保証構造が明快だったが，高桁で逐次Rational級数のcostが急増し，`N[Pi,10000]`が約12秒級になった。

v1.5.1では旧Machin実装をreferenceとして残し，既定を**binary-splitting Chudnovsky**へ変更した。巨大整数演算はKaratsuba/Toom/Burnikel–Ziegler等の下位backendを再利用し，最終的な区間化・丸め保証はRealInterval側の契約を維持する。

代表値:

```text
N[Pi,10000]: 約12.3 s → 約0.4 s
```

Chudnovsky内部の係数もterm indexが大きい場合にmachine 64bit overflowへ依存しないようBigIntで構築する。

## 37.2 exp / E

逐次RealInterval Taylorは高桁でinterval objectとRational中間値のcostが大きかったため，v1.5.1では**binary splitting + certified range reduction**へ移行した。

小さいexact Rationalではexact binary splitting，大きい分子・分母では固定precision interval binary splittingを使い，中間exact integerの異常膨張を避ける。

代表測定では`N[E,5000]`が数秒級から約0.1秒級まで改善した。

## 37.3 log

`log`も逐次atanh型級数からbinary splittingへ変更した。高bit mantissaではそのままexact splitすると中間整数が膨張するため，certified sqrtを複数回使って1近傍へ縮約し，interval binary splitting後にscaleを復元する。

`log[2]`のような小係数caseはexact binary splittingを維持する。

代表測定では`N[log[2],3000]`が約6～7秒から約0.1秒級へ改善した。

## 37.4 巨大Radianの三角函数

巨大なexact Rational radianをTaylorへ直接投入するとargument magnitudeに比例して実用不能になる。v1.5.1ではPiの保証区間を使い，`x/(Pi/2)`の象限integerが一意に確定するまでprecisionを確保してから，小区間へargument reductionする。

machine `fmod`や`double` Piへ落とさないため，certified semanticsを保つ。

point Rationalだけでなく，`sqrt[2]`を含むような保証区間全体にも縮約を拡張している。

# 38. 現行実装の計算量上の特徴

大まかには次の通り。threshold以下ではより単純な算法へ戻るため，表は巨大operand側の性格を示す。

| 演算 | v1.5.2現在の主要算法 | 備考 |
|---|---|---|
| BigUInt add/sub | linear carry/borrow | `O(n)` |
| BigUInt compare | 上位から比較 | `O(n)` worst |
| BigUInt shift | limb移動 + bit shift | `O(n)` |
| BigUInt multiply | schoolbook → Karatsuba → Toom-3 | size/shapeでdispatch |
| BigUInt square | symmetric schoolbook → Karatsuba square | `x*x`専用 |
| BigUInt divide | special path → Knuth → Burnikel–Ziegler | balanced huge divisionでBZ |
| decimal conversion | `10^9` chunk + divide-and-conquer | 巨大10進表示向け |
| gcd | Euclidean `%` | BZ除算の改善を間接利用 |
| pow | exponentiation by squaring | 専用squareを利用 |
| factorial | balanced product tree | adaptive multiplicationを利用 |
| integer sqrt/cbrt | Newton | BZ division / square改善が波及 |
| BigFloat add | exact alignment + safe exponent-gap fast path | directed rounding保持 |
| BigFloat mul | BigInt adaptive multiplication + rounding | exact significand積 |
| BigFloat div | BigInt divmod + rounding | BZが波及 |
| Pi | binary-splitting Chudnovsky | certified enclosureへ接続 |
| exp/log | binary splitting + range reduction | certified interval |
| trig huge radian | certified argument reduction | Pi enclosureを利用 |

v1.5.1で導入した乗算・除算・10進I/O・主要超越函数の高速化をv1.5.2でも維持しており，大きなquadratic/逐次bottleneckはかなり緩和されている。一方，さらに巨大な整数ではhigher Toom / FFT系，高桁`log`ではbit-burst/AGM系などが次候補になる。

# 39. メモリとサイズの実際の上限

「任意精度」は数学的に固定桁数を設けていないという意味であり，物理的に無限ではない。さらに，mmCalでは**limbそのものよりC++ objectの固定費が支配する領域**がある。

## 39.1 `BigUInt`のpayloadと固定費

`BigUInt`は

```cpp
std::vector<std::uint32_t> limbs_;
```

を持つ。したがって実際の最大長は，

- `vector::max_size()`，
- address space，
- available memory，
- 各算法の一時buffer，

で制約される。

0はempty vectorなのでheap payloadを必要としない。一方，1 limbの小整数でも通常の`std::vector`である以上，非zero payload用heap allocationが発生し得る。

また`normalize()`は上位0 limbを`pop_back()`するが，`vector::capacity()`を毎回縮めない。これは高速化には妥当であるが，一度巨大化した長寿命objectが小さくなっても高水位capacityを保持することはある。無条件`shrink_to_fit()`はallocator churnを増やすため，現在は採用していない。

## 39.2 BigIntの符号は数値limbを増やさない

前述の通りsign-magnitudeなので，負数の符号拡張limbはない。符号のコストはC++ object内の`bool`とpaddingであり，数値桁数に比例して増えるものではない。

## 39.3 BigFloatの指数範囲

仮数はBigIntだが，指数は`std::int64_t`である。したがって指数範囲は有限であり，overflow/underflowは明示的に検査する。

`precisionBits`は`std::size_t`だが，巨大shift，precision conversion，allocationには別の現実的上限がある。

## 39.4 x86-64 GCCでの参考`sizeof`

以下はv1.5.2 sourceをx86-64 GCCで測った**参考値**であり，MSVC ABIや将来実装の仕様ではない。

| 型 | `sizeof`参考値 | 主な理由 |
|---|---:|---|
| `BigUInt` | 24 B | `std::vector<uint32_t>` object |
| `BigInt` | 32 B | sign + padding + BigUInt |
| `Rational` | 64 B | BigInt × 2 |
| `BigFloat` | 48 B | BigInt + exponent + precision |
| `RealNumber` | 72 B | `variant<BigInt,Rational>` |
| `Number` | 152 B | real/complex variantの最大payload |
| `DecimalApproximation` | 248 B | text + metadata + Rational bounds |
| `ComplexDecimalApproximation` | 536 B | real/imag approximate metadata |
| `Expr` handle | 16 B | `shared_ptr` |
| `Expr::Node::Value`相当variant | 544 B | 最大alternativeをinline保持 |

ここから重要なことが分かる。

**「BigIntが符号1bitのために1 limb損している」ことは問題ではない。小さい数でも汎用C++ objectを何層も通る固定費の方が桁違いに大きい。**

## 39.5 `Expr::Node`が巨大dense Arrayで本丸になる理由

現在の`Expr::Node`は，

```cpp
std::variant<
    Number,
    DecimalApproximation,
    ComplexDecimalApproximation,
    bool,
    std::string,
    Symbol,
    ArrayExpr,
    ListExpr,
    CallExpr,
    std::shared_ptr<const SolutionSet>>
```

をinline保持する。`std::variant`は最大alternativeを入れられるだけの領域を全Nodeへ確保するため，小整数`1`のNodeでも`ComplexDecimalApproximation`級の箱を払う。

v1.5.2の1024×1024 dense Matrix監査では，添付形式と同等の10桁Rational Exprを1,048,576個C++から直接構築しただけで最大RSS約0.69 GB，約14.16 MBのテキストをCLIでparseし`dimensions[...]`を求めるだけでは最大RSS約1.99 GBだった。

このため次版ToDoでは，算法block化より前に，

1. `Expr::Node`をkind別typed nodeへ分離し巨大variant固定費を除去する。
2. BigUInt / BigIntへsmall-object optimizationを検討する。
3. numeric Array / approximate Matrixへpacked storageを検討する。
4. 巨大brace literalのparser/lowering allocationを削減する。

という順を候補にしている。

これは多倍長算術の意味論を変える最適化ではなく，**同じ値をもっと薄い器へ入れる**ためのrepresentation refactorである。

## 39.6 なぜBigUInt SBOを先にやらないのか

small-object optimizationで1～2 limb整数のheap allocationを消す価値は高い。しかし1024² Matrixの現状では，1要素あたり数byte～数十byteを節約する前に，数百byte級の`Expr::Node`固定費がある。

したがって費用対効果としては，

```text
Expr Node fixed cost
    ↓
small BigInt allocation
    ↓
packed numeric Array
```

の順に単独benchmarkする方が原因を切り分けやすい。

## 39.7 CertifiedEvaluatorの深さ制限

数値桁数とは別に，病的な深さのASTでC++ call stackを破壊しないため，certified expression depthには安全上限がある。これは多倍長値の桁数制限ではなく，式木評価の安全制限である。

# 40. 採用していない・まだ存在しないもの

v1.5.2現在でも「未実装」と「実装して比較したが棄却」を分ける。

## 40.1 比較したが既定採用しなかったもの

- Prime-Swing factorial — 現product treeより巨大factorialで遅かった
- binary GCD全面置換 — 現Euclidean `%`より約4～30倍遅いcaseがあった
- Karatsuba workspace pool — 管理costで約5～10%退行
- Karatsuba depth scratch — 同様に退行
- Toom-3専用square — Karatsuba squareより遅かった

これらはコードベースやCPU特性が変われば再評価可能であり，算法一般を否定しているわけではない。

## 40.2 まだ導入していないもの

- Toom-4 / higher Toom
- FFT/NTT integer multiplication
- Lehmer GCD
- bit-burst / AGM系の超高精度`log`
- small-buffer optimization for limbs
- limb-level custom allocator
- arbitrary-size BigFloat exponent
- MPFR/GMP/Boost.Multiprecision backend
- Exact/Certifiedと暗黙に混在するMachine/double backend

最後の項目は単なる未実装ではなく，通常意味論を速度のためにmachine精度へ落とさないという設計判断でもある。

# 41. 現行テストで確認している主な性質

## BigUInt

- 64bit境界超えのparse
- 基数2..36 round trip
- limb carry / borrow
- schoolbook/Karatsuba/Toom-3 multiplicationの境界・random照合
- 専用squareと一般乗算の一致
- shift境界 31/32/33/63/64bit等
- Knuth/Burnikel–Ziegler divisionのreconstruction invariant
- divisor top-limb全bit位置のnormalization
- 固定seedの巨大random division / decimal round-tripを`mmCal.Benchmarks`でも継続検証

## BigInt

- `INT64_MIN` constructor
- negative zero normalize
- signed arithmetic
- division truncate-toward-zero
- remainder sign = dividend sign
- `q*d+r == dividend`
- `int64_t`範囲のgrid cross-check

## Integer algorithms

- gcd/lcm
- `2^256`
- `100!`
- uint64 conversion境界
- 巨大整数sqrt
- sqrt invariantを0..5000で確認
- cbrt invariantを0..5000で確認

## Rational

- canonical reduction
- denominator sign normalization
- exact decimal parse
- base2/base16 fractional literal
- cross-cancellation
- 最適化演算と「巨大積を作ってconstructor normalizeするreference」の比較を160ケース

## BigFloat

- dyadic normalization
- 1/3の上下方向丸め
- negative directed rounding
- ties-to-even midpoint
- exact四則
- directed conversion property:
  - numerator -31..31
  - denominator 1..23
  - precision 1..18 bit
- directed arithmetic property:
  - dyadic grid -9..9
  - precision 1..12 bit

## RealInterval

- non-dyadic Rationalがpointにならない
- dyadic Rationalがpointになる
- add/sub/mul/divがexact resultを包含
- zero-containing denominator intervalを拒否

## DecimalApproximation

- terminating decimal
- repeating decimal
- nearest-even rounding
- integer carry
- fixed digits
- certified interval両端が同じ丸めならaccept
- rounding boundaryを跨ぐintervalはreject

---

# 42. v1.5.2以降の性能候補

v1.5.1でKaratsuba/Toom-3，専用square，Burnikel–Ziegler，10進D&C，BigFloat exponent-gap fast pathまで導入し，v1.5.2でもその基盤を維持しているため，次の候補は一段上になる。

## 42.1 higher multiplication

- Toom-4 / higher Toomのcrossover測定
- さらに巨大な整数向けFFT/NTT multiplication

ただし現在のinteractive用途ではToom-3までで十分な領域も広く，実測でcrossoverが現れるまで複雑化しない。

## 42.2 GCD

RationalはGCDを頻繁に使う。binary GCDは現backendでは退行したため，次に試すならLehmer GCD等が候補。

## 42.3 超高精度log/exp

binary splittingで数千桁は大幅改善したが，さらに高桁では`log`のbit-burst / AGM系，`exp`のrange reduction調整を比較する価値がある。

## 42.4 特殊函数横断benchmark

Gamma / erf等を5000～10000桁まで振り，逐次級数やinterval object生成が新たな崖にならないか`mmCal.Benchmarks`へ追加する。

性能候補は`performance_optimization.ja.md`に採用・棄却理由を残し，threshold変更時は固定seed正当性試験を先に通す。

## 42.5 representation最適化

v1.5.2の大行列監査で，算術algorithmよりrepresentation固定費が先に壁になる領域が確認された。次の候補は，

- `Expr::Node`巨大variantのkind別typed node化，
- BigUInt / BigInt small-object optimization，
- numeric Array / approximate Matrix packed storage，
- 巨大brace parser/loweringのallocation削減，

である。

これらは一括で変更しない。まずNodeだけ，次にBigInt SBOだけ，というように全regressionとRSS benchmarkを固定して採否を測る。

# 43. 実装上の強み

現行多倍長基盤の強みは，単に「大きな数字を扱える」ことではない。

## 43.1 exact integerとapproximationの境界が明確

`BigInt/Rational` と `BigFloat` を混同していない。

## 43.2 BigFloatを真値扱いしない

方向丸めをAPI引数で明示し，RealIntervalの包含保証へ使う。

## 43.3 BigFloat自身もexact dyadicとして検査できる

`toRational()` によりBigFloat端点をexact Rationalへ戻せる。

そのためテストでも

```text
lower <= exact <= upper
```

をexact arithmeticで直接検証できる。

## 43.4 Rationalの中間値膨張を意識している

加算のGCD縮小，乗除算の交差約分が既に入っている。

## 43.5 integer rootを上位certificationに再利用している

`sqrt` の保証付き評価をmachine sqrtへ丸投げせず，BigInt固定小数点問題へ落としている。

---

# 44. 実装上の弱点・今後の注意

## 44.1 thresholdは環境依存

Karatsuba / Toom-3 / Burnikel–ZieglerのcrossoverはCPU，compiler，allocator，cacheに依存する。v1.5.1値を普遍的な定数として扱わず，MSVC等では`mmCal.Benchmarks`で再測定する。

## 44.2 さらに巨大な整数

Toom-3より上のhigher Toom / FFT/NTTは未実装。数十万～数百万bit級を頻繁に扱う用途では次の構造的bottleneckになり得る。

## 44.3 GCD

Euclidean GCDはBZ除算の恩恵を受けるが，巨大Rationalの反復約分ではLehmer系の余地がある。binary GCDは実測退行したため採用していない。

## 44.4 高精度超越函数

`Pi`, `exp`, `log`, 巨大trigはv1.5.1で大幅改善したが，さらに高桁では`log`のbit-burst/AGM系や特殊函数固有のasymptotic/binary-splitting backendが候補になる。

## 44.5 directed roundingを壊さないこと

BigFloat backend最適化で最重要。NearestEvenだけ正しくても不十分で，

```text
TowardNegative
TowardPositive
```

が1 ulpでも内側へ入るとRealInterval全体の保証が壊れる。そのためexponent-gap fast pathも証明可能caseだけに限定し，曖昧なら旧exact alignmentへfallbackする。

## 44.6 performance testを正当性testと混同しない

速い結果が正しい証拠にはならない。`mmCal.Benchmarks`はthreshold sweepと固定seedrandom invariantを同じprojectへ置くが，通常のUnit/black-box regressionとは役割を分ける。

# 45. 既存基盤の要約

v1.5.2の多倍長・保証付き数値基盤は次のように整理できる。

```text
[ exact integer core ]
std::vector<uint32_t> limbs
        ↓
BigUInt
  multiply: schoolbook → Karatsuba → Toom-3
  square  : dedicated symmetric/Karatsuba path
  divide  : special → Knuth → Burnikel–Ziegler
  decimal : 10^9 chunk + divide-and-conquer
        ↓ sign-magnitude
BigInt
        ↓ normalized numerator/denominator
Rational
        ↓
RealNumber / exact Complex Number

[ arbitrary-precision working arithmetic ]
BigInt significand × 2^(int64 exponent)
        ↓ directed rounding
BigFloat
        ↓ outward rounding
RealInterval / ComplexInterval
        ↓ whole-expression enclosure
CertifiedEvaluator
        ↓ interval endpoints round identically
DecimalApproximation
```

上位の高精度函数も，

```text
Pi      → binary-splitting Chudnovsky
exp/log → binary splitting + certified range reduction
trig    → certified argument reduction
sqrt    → exact fixed-point integer root
```

のように，下位BigIntの高速化を再利用しつつ保証区間へ接続する。

v1.5.1～v1.5.2を通して重要なのは「高速算法を入れたこと」そのものではなく，**採用をbenchmarkで決め，速くならなかった算法は棄却し，Exact/Certifiedの意味論を変えないこと**である。

Prime-Swing，binary GCD，workspace化，Toom-3 squareは実際に試したが現環境では採用しなかった。逆にKaratsuba/Toom-3，専用square，Burnikel–Ziegler，decimal D&C，Chudnovsky，binary-splitting exp/logは実測利益と正当性試験の双方を確認して採用した。

今後backendを更新する場合も，最優先で守る境界は

```text
BigUInt/BigIntのexact arithmetic
BigFloatのdirected rounding
RealIntervalのoutward containment
```

の3点である。詳細な採用・棄却履歴と代表benchmarkは`performance_optimization.ja.md`を参照する。


# 46. 「自作多倍長」を層ごとに見る

多倍長を1個の巨大classとして実装すると，exact整数，符号，有理約分，丸め，区間保証が混ざる。mmCalは責務を分けている。

```text
BigUInt
  「桁列をどう足し，引き，掛け，割るか」

BigInt
  「符号をどう付けるか」

Rational
  「分子分母をどうcanonicalに保つか」

BigFloat
  「有限precisionへどちら向きに丸めるか」

RealInterval
  「真値を外へ逃がさないか」

DecimalApproximation
  「保証された内部結果をどう10進表示へ確定するか」
```

この分離により，例えばKaratsubaをToomへ差し替えても`integrate`や`solve`はlimbを知る必要がない。またBigFloatの高速化をしても，`TowardNegative/TowardPositive`契約を守る限りRealIntervalの証明構造は維持できる。

# 47. 実装者が踏みやすい罠

## 47.1 「任意精度floatだからexact」と思う

誤りである。`BigFloat::fromRational(1/3,p)`は有限dyadicへ丸める。保証が必要なら上下方向へ別々に丸め，intervalにする。

## 47.2 乗算だけ高速化する

巨大整数では，乗算を速くすると次にdivision，GCD，decimal conversionがbottleneckとして現れる。mmCalでもfactorial本体の高速化後，10進表示が支配的になった。

## 47.3 Rationalを毎回`ad+bc / bd`で作る

数学的には正しいが，中間値を巨大化させる。多倍長では「最終値が小さい」ことと「途中も小さい」ことは別問題である。

## 47.4 directed roundingをnearestで代用する

1 ulp内側へ入っただけでinterval certificateが壊れる。高速化のfast pathは，全rounding modeで正しいことを別々に証明する必要がある。

## 47.5 moveを数学的に無害だと思う

`[x,x]`を作るだけでも，同一C++ objectをcopy/move混在させれば評価順差で壊れ得る。数値証明の最下層ではobject lifetimeも正当性の一部である。

## 47.6 benchmarkで速いから採用する

random invariantや境界試験を先に通す。mmCalではPrime-Swing，binary GCD，workspace Karatsuba，Toom-3 squareなどを実装した上で，現backendでは遅かったため採用しなかった。

# 48. コードを読むならこの順

初見で`big_uint.cpp`の5万行級実装へ飛び込むより，次の順が理解しやすい。

1. `numeric/detail/big_uint.hpp` — limb APIとinvariantを見る。
2. `numeric/big_int.hpp` — sign-magnitudeの薄いwrapperであることを見る。
3. `numeric/rational.cpp` — cross-cancelとcanonicalizationを見る。
4. `numeric/big_float.hpp/.cpp` — dyadic表現とrounding modeを見る。
5. `approximation/real_interval.cpp` — outward roundingがどこで入るかを見る。
6. `numeric/decimal_approximation.cpp` — 最終10進表示の確定条件を見る。
7. `approximation/certified_evaluator.cpp` — 式全体をどうinterval化するかを見る。
8. `builtins/signal_processing.cpp`，`linear_algebra/*` — precision-aware `N`がexact式構築を避ける実例を見る。
9. `benchmarks/benchmark_main.cpp` — thresholdが理論ではなく実測で決められていることを見る。

# 49. 自作多倍長の設計チェックリスト

mmCal以外で同種の型を書く場合にも使える最低限の確認項目である。

- zero representationは1種類か。
- negative zeroは残らないか。
- limbの最大積＋carryが中間型へ収まるか。
- self-assignment / self-additionで参照無効化しないか。
- shift量0，limb幅境界，巨大shiftを扱えるか。
- divisionは`q*d+r==n`とremainder範囲を常に満たすか。
- signed divisionの丸め方向を仕様化しているか。
- Rationalは分母正・既約・zero canonicalを維持するか。
- BigFloatのrounding modeがAPI上明示されているか。
- nearest-evenのtie判定をexact remainderで行えるか。
- interval演算は必ず外向き丸めか。
- 0を含む区間によるdivisionを拒否するか。
- decimal表示と内部precision metadataを混同していないか。
- compiler差のある評価順やmove状態へ依存していないか。
- algorithm thresholdをCPU非依存の数学定数だと思っていないか。
- 大値だけでなく，小値のobject fixed costも測ったか。

# 50. v1.5.2時点での結論

mmCalの自作多倍長基盤は，

```text
32bit limb exact integer
    ↓
sign-magnitude BigInt
    ↓
canonical Rational
    ↓
exact dyadic BigFloat + directed rounding
    ↓
outward Real/Complex Interval
    ↓
certified decimal result
```

という一貫した層構造になっている。

性能面では，schoolbook/Karatsuba/Toom-3，専用square，Knuth/Burnikel–Ziegler，10進divide-and-conquer，Chudnovsky，binary-splitting exp/logまで入り，単純な「自作BigInt」の域はかなり越えている。

一方，v1.5.2の1024 dense Matrix監査で，次のbottleneckは算術algorithmだけではなく**representation**であることも明確になった。特に`Expr::Node`の巨大variant固定費は，small BigIntの符号やlimbより影響が大きい。

したがって次段の性能改善では，

```text
意味論を変えない
    ↓
器を薄くする
    ↓
その後に算法をさらに高度化する
```

という順が合理的である。

mmCalが守るべき核心は，最後まで次の3点である。

```text
BigUInt / BigInt のexact arithmetic
BigFloat の明示的directed rounding
RealInterval / ComplexInterval のoutward containment
```

ここを守る限り，内部表現や乗算algorithmは将来いくらでも交換できる。逆に，ここを曖昧にして得た高速化はmmCalの「exact-first，近似は明示的，解らないものは解らないと言う」という設計思想そのものを壊す。
