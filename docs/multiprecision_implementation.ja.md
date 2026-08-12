# mmCal 多倍長数値基盤 実装詳細

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
| unsigned integer | `src/numeric/detail/big_uint.hpp/.cpp` | 32bit limb，四則，shift，Knuth型除算，基数変換 |
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

この周辺だけで約3,400物理行ある。

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

# 6. `BigUInt` の乗算

現在は**schoolbook multiplication**である。

Karatsuba，Toom-Cook，FFT multiplication等は実装していない。

概念的には

```text
for i in lhs limbs:
    carry = 0
    for j in rhs limbs:
        k = i + j
        t = lhs[i] * rhs[j] + result[k] + carry
        result[k] = low32(t)
        carry = high32(t)
```

である。

32bit limb × 32bit limbを64bit accumulatorで受ける。

最悪値でも

```text
(B-1)^2 + (B-1) + (B-1)
= B^2 - 1
= 2^64 - 1
```

なので `uint64_t` にちょうど収まる。

limb数を `n,m` とすると計算量は

```text
O(nm)
```

同程度の長さなら

```text
O(n^2)
```

である。

これは現在のBigIntで巨大数乗算が重くなる主要因の一つである。

---

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

多倍長基盤の中でも重要な部分である。

## 8.1 fast path

次のケースを先に処理する。

```text
divisor == 0        → domain_error
dividend < divisor  → quotient=0, remainder=dividend
dividend == divisor → quotient=1, remainder=0
```

除数が1 limbなら `divideSmall()` を使う。

`divideSmall()` は上位limbから

```text
current = remainder * B + limb[i]
quotient[i] = current / divisor
remainder = current % divisor
```

と進む通常のlong divisionである。

## 8.2 multi-limb division

2 limb以上では**Knuth型のnormalized long division**を使う。

内部基数を

```text
B = 2^32
```

とする。

### Step 1: divisorの正規化

除数最上位limbのleading zero数

```cpp
normalizationShift = std::countl_zero(divisor.highestLimb)
```

だけ，被除数・除数の双方を左shiftする。

これにより除数最上位limbの最上位bitが1になり，商digit推定の条件を良くする。

### Step 2: 商1 limbの推定

除数長を `m`，対象位置を `j` とすると，被除数の上位2 limbから

```text
numerator = u[j+m] * B + u[j+m-1]
```

を作り，

```text
qhat = numerator / v[m-1]
rhat = numerator % v[m-1]
```

と推定する。

### Step 3: 次limbを使った補正

```text
qhat >= B
```

または

```text
qhat * v[m-2]
> B*rhat + u[j+m-2]
```

なら `qhat` を1減らして補正する。

### Step 4: `qhat * divisor` を減算

`subtractProduct()` で対象limb区間から

```text
qhat * divisor
```

を引く。

推定が1大きすぎてborrowが最上位まで抜けた場合は，

```text
qhat--
addBack(divisor)
```

で1回戻す。

### Step 5: remainderのde-normalize

最後にremainderを `normalizationShift` だけ右shiftして元のscaleへ戻す。

## 8.3 計算量

被除数 `n` limbs，除数 `m` limbsなら，現在のlong divisionは概ね

```text
O((n-m+1)m)
```

で，同程度の桁数なら `O(n^2)`。

Burnikel-Ziegler等の高速除算は現在ない。

## 8.4 テスト

現行テストでは，

- limb境界を跨ぐ除算
- 257bit / 129bitの2冪除算
- `q*d+r == dividend`
- `r < divisor`
- divisor最上位bit位置を0..31まで変えた正規化ケース
- 128件のrandom multi-limb input

を検証している。

random testでは別実装のbinary long divisionをreferenceとして比較している。

---

# 9. `BigUInt` の文字列変換

## 9.1 parse

基数2..36に対応。

各digitについて

```text
value = value * radix + digit
```

を繰り返す。

つまり decimal parse もmachine integerへ一旦収めず，最初から多倍長として構築する。

## 9.2 toString

逆に

```text
while value != 0:
    remainder = value % radix
    value /= radix
    output remainder
reverse(output)
```

とする。

実装は簡潔で正しいが，10進変換を `10^9` 等のchunkで処理する方式ではないため，非常に巨大な整数の文字列化については高速化余地がある。

---

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

単純に

```text
1*2*3*...*n
```

と左から掛け続けない。

`productRange(first,last)` でbalanced product treeを作る。

小区間，現在は

```text
last - first <= 15
```

だけstraight loopとし，それ以上は中央で二分する。

目的は，極端にサイズの違う巨大BigIntを順次掛け続ける形を避けること。

ただし基礎乗算自体はschoolbookなので，factorial全体が高度なprime-swing算法等になっているわけではない。

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

加算では指数を揃える。

```text
commonExponent = min(lhs.exponent, rhs.exponent)
```

上位指数側significandを左shiftして，双方をcommon exponentのexact整数へ揃える。

```text
lhs = L * 2^e1
rhs = R * 2^e2

common = min(e1,e2)
L' = L << (e1-common)
R' = R << (e2-common)
resultExactSignificand = L' + R'
```

そのexact和を `fromDyadic()` で指定precisionへ丸める。

### 現行実装上の注意

指数差が非常に大きい場合，加算は大きなleft shiftを作る。

つまり「小さい項はprecision上影響しない」と先に判定してsticky bitだけ処理する形式ではない。

これは正確で単純だが，極端なexponent gapではmemory/performance上の改善余地がある。

---

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

仮数乗算の速度は基礎 `BigInt`，すなわち現在はschoolbook `O(n^2)` の影響を直接受ける。

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

# 37. 高精度定数への接続例: Pi

Piの現行certified実装ではMachin公式

```text
Pi = 16 atan(1/5) - 4 atan(1/239)
```

を使う。

`atan(1/q)` の交代級数各項はexact `Rational` として構築され，それを `RealInterval::fromRational()` でBigFloat上下界へ変換する。

つまり，級数項自体をmachine doubleで生成しない。

交代級数の隣接部分和による数学的な包含と，BigFloat演算の外向き丸めの双方をRealIntervalへ吸収する。

Pi計算法として最速を目指す実装ではないが，BigInt/Rational/BigFloat/RealIntervalが一体として動くreference例になっている。

---

# 38. 現行実装の計算量上の特徴

大まかには次の通り。

| 演算 | 現行算法 | limb計算量の目安 |
|---|---|---:|
| BigUInt add/sub | linear carry/borrow | `O(n)` |
| BigUInt compare | 上位から比較 | `O(n)` worst |
| BigUInt shift | limb移動 + bit shift | `O(n)` |
| BigUInt multiply | schoolbook | `O(nm)` |
| BigUInt divide | normalized long division | `O((n-m+1)m)` |
| BigInt | BigUInt + sign処理 | 基礎演算に準ずる |
| gcd | Euclidean `%` | 除算コスト依存 |
| pow | exponentiation by squaring | `O(log exponent)`回の乗算 |
| factorial | balanced product tree | 乗算backend依存 |
| integer sqrt | Newton | 反復ごとに巨大除算・乗算 |
| integer cbrt | Newton + correction | 同上 |
| Rational add | GCD縮小付き | GCD/乗算依存 |
| Rational mul/div | 交差約分 | GCD/乗算依存 |
| BigFloat mul | BigInt積 + rounding | BigInt乗算依存 |
| BigFloat div | BigInt divmod + rounding | BigInt除算依存 |

現在，BigInt乗算にKaratsuba等がないため，数千～数万bitを大規模に扱う処理ではここが将来的なボトルネック候補になる。

---

# 39. メモリとサイズの実際の上限

「任意精度」は数学的に固定桁数を設けていないという意味であり，物理的に無限ではない。

## BigUInt

```text
std::vector<uint32_t>
```

なので，実際の最大長は

- `vector::max_size()`
- address space
- available memory

で制約される。

乗算・左shiftではoverflow前にサイズ検査がある。

## BigFloat

仮数はBigIntだが，指数は `int64_t`。

したがって指数範囲は有限。

要求precisionは `size_t` だが，巨大shiftやdecimal→binary precision変換でoverflowを明示的に検出する。

## CertifiedEvaluator

数値桁数とは別に，病的な深さのASTでC++ call stackを破壊しないため，certified expression depthに96段の安全上限がある。

これは多倍長値の桁数制限ではなく，式木の深さ制限である。

---

# 40. 現行実装で意図的に存在しないもの

この文書の基準ソースでは，次のような高度な多倍長最適化はまだない。

- Karatsuba multiplication
- Toom-Cook multiplication
- FFT/NTT integer multiplication
- Burnikel-Ziegler division
- Lehmer GCD
- binary GCDへの全面置換
- chunked decimal conversion (`10^9`単位等)
- small-buffer optimization for limbs
- limb-level custom allocator
- IEEE互換NaN/Infinityを持つBigFloat
- arbitrary-size BigFloat exponent
- sticky-bitベースの巨大exponent-gap加算fast path
- MPFR/GMP/Boost.Multiprecision backend

これは「未完成」というより，現行1.5系がまず意味論・exactness・certificationを優先し，基礎算法を自前で明快に保っている結果である。

将来高速化する場合も，上位のRational・BigFloat・RealIntervalの意味論を変えず，最下層backendだけを段階的に差し替えられる構造になっている。

---

# 41. 現行テストで確認している主な性質

## BigUInt

- 64bit境界超えのparse
- 基数2..36 round trip
- limb carry / borrow
- multi-limb multiplication
- shift境界 31/32/33/63/64bit等
- Knuth型divisionのreconstruction invariant
- divisor top-limb全bit位置のnormalization
- random 128ケースを別binary division referenceと照合

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

# 42. 多倍長基盤から見た性能上の優先候補

現在の意味論を変えずに性能を改善するなら，候補は概ね次の順になる。

## 42.1 multiplication thresholdの導入

小さい値:

```text
schoolbook
```

大きい値:

```text
Karatsuba
```

さらに巨大ならToom/FFT系，というdispatchが考えられる。

BigInt multiplicationが改善されれば，

- Rational
- factorial
- integer root
- BigFloat multiply
- certified transcendental

へ広く波及する。

## 42.2 decimal conversionのchunk化

現在の1 digitずつの

```text
*10
/10
```

を，例えば内部的に `10^9` chunkへすれば，巨大整数のI/O負荷を減らせる。

これは演算意味論へ影響しにくい。

## 42.3 GCD

RationalはGCDを頻繁に使うため，非常に巨大な値ではEuclidean `%` のコストが効く。

Lehmer GCD等は候補になる。

## 42.4 BigFloat addのexponent-gap fast path

現在はcommon exponentへexact shiftしてから丸める。

要求precisionに比べて一方が極端に小さい場合，guard/sticky情報だけを残して巨大shiftを避ける最適化余地がある。

ただしdirected roundingの正しさを崩しやすいため，これは整数backend高速化より慎重に扱うべきである。

---

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

## 44.1 大整数乗算はquadratic

現状最大の構造的性能限界。

## 44.2 除算もquadratic系

BigFloat，Rational GCD，Newton rootに広く効く。

## 44.3 入出力変換がdigit-wise

極端に巨大なdecimal I/Oには非効率。

## 44.4 BigFloat additionがexact alignment型

exponent gapが巨大な場合に大きなtemporary BigIntを作り得る。

## 44.5 `tryToUint64` がdecimal string経由

現在は

```text
BigInt -> toString() -> from_chars()
```

で変換する箇所がある。

頻繁に呼ぶhot pathになればlimbから直接判定・変換するAPIを追加する価値がある。

## 44.6 高速化でdirected roundingを壊さないこと

BigFloat backendを最適化するとき最も重要。

NearestEvenだけ正しくても不十分で，

```text
TowardNegative
TowardPositive
```

が1 ulpでも内側へ入るとRealIntervalの保証全体が壊れる。

---

# 45. まとめ

現行mmCalの多倍長基盤は，次のように整理できる。

```text
[ exact integer core ]
std::vector<uint32_t> limbs
        ↓
BigUInt
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

最下層のBigUIntは，

- base `2^32`
- 32bit limb / 64bit accumulator
- schoolbook multiplication
- Knuth型normalized long division

という比較的保守的で検証しやすい構成。

BigIntはsign-magnitude，Rationalは常時既約，BigFloatはcanonical dyadic，RealIntervalはdirected outward roundingという形で，それぞれの層に明確なinvariantがある。

mmCalにおける「多倍長」の本質は，BigIntだけではない。

```text
任意長整数
→ exact Rational
→ exact dyadic arbitrary-precision work value
→ certified interval
→ 一意に確定した10進表示
```

までが一続きの数値設計になっている。

現状の主要な改善余地は高速算法であり，意味論の基盤そのものは既に分離されている。したがって今後Karatsuba等を導入する場合も，最優先で守るべき境界は

```text
BigUInt/BigIntのexact arithmetic
BigFloatのdirected rounding
RealIntervalのoutward containment
```

の3点である。
