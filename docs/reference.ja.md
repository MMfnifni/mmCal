# mmCal 仕様・函数リファレンス

この文書は **mmCalの実装そのもの** を基準にした詳細仕様書である。
ユーザー向けの導入はルートの`README.ja.md`を参照する。

> 対象: **v1.5.4 開発版**

この文書はバージョンごとに常に変動するため，過去バージョンはgitより引っ張り出してください。

## 0. 不変の理念

### 設計理念

厳密に。全て自前で。近似は明示的に。解らないものは解らないと言う。

### 配布理念

唯一本の実行ファイルに。
Open source under the BSD 3-Clause License.

## 1. 現在の設計思想

mmCalは，入力を最初から`double`へ落とす電卓ではなく **exact-firstの小型CAS / 数値計算kernel** とする。
また整数・有理数・有限小数を可能な限り厳密値として扱い，代数数，記号式，微積分，方程式，複素数，行列，FFTまで同じexact-firstの体系上で処理する。

優先順位は次の通り。

1. 数学的厳密性
2. 主値分岐 / definedness / 定義域を失わないこと
3. exactとapproximateを混同しないこと
4. 共通知識をSolver・Simplifier・CertifiedEvaluatorで再利用すること
5. 上記を壊さない範囲で高速化すること

代表例:

```text
0.1 + 0.2
-> 3/10

sin[Pi/6]
-> 1/2

sqrt[-8]
-> 2I sqrt[2]

log[-1]
-> I Pi

N[Pi,30]
-> 3.14159265358979323846264338328
```

`Pi`や`sqrt[2]`は「内部に保存した小数」ではない。exactな数式として保持し，`N[...]`が指定されたときだけcertified numerical evaluationへ進む。

---

# 2. 数値・式の内部モデル

## 2.1 Integer

任意長整数`BigInt`。

```text
123456789012345678901234567890
-> 123456789012345678901234567890

100!
-> exact BigInt
```

固定64bit整数へ丸めない。

## 2.2 Rational

有限小数も最初からexact Rationalとしてparseする。

```text
0.1
-> 1/10

1.25
-> 5/4

0.1 + 0.2 == 0.3
-> True
```

したがってIEEE 754由来の`0.30000000000000004`問題は通常評価には存在しない。

## 2.3 exact complex

実部・虚部ともexactな`Number`として保持する。

```text
I^2
-> -1

(3 + 4I) / 5
-> 3/5 + 4/5I
```

## 2.4 Symbolic expression

exactに閉じない式はASTのまま保持する。

```text
sqrt[2]
-> sqrt[2]

sin[1]
-> sin[1]

gamma[1/3]
-> gamma[1/3]
```

ここで`sin[1]`の1は **1 radian**。既定角度はRadianである。

## 2.5 certified approximation

任意精度の作業値は`BigFloat`，証明付き区間は`RealInterval` / `ComplexInterval`。

`N[expr,p]`では，真値を含む区間の両端が同じ`p`有効桁の10進丸めへ入ることを確認してから`DecimalApproximation`を返す。0近傍で相対Precisionを定義できない場合でも，InformationEnclosureから絶対Accuracyを保証できるならzero-centered approximationを返せる。

現在の`DecimalApproximation`は表示文字列だけではなく，要求桁数，由来（exact入力 / certified interval），表示10進値そのもののexact Rationalに加えて，**CertifiedEnclosure**と**InformationEnclosure**の2種類のexact Rational区間を保持する。CertifiedEnclosureは真値包含を証明する区間，InformationEnclosureはその近似値から後続計算で利用してよい情報量を表す区間であり，常に`CertifiedEnclosure ⊆ InformationEnclosure`を満たす。`ComplexDecimalApproximation`も実部・虚部ごとに同じmetadataを保持する。`precision/accuracy/rationalize`はInformationEnclosureを直接使い，表示文字列を再parseして精度を推測しない。

```text
N[sqrt[2],30]
-> 1.41421356237309504880168872421
```

「差が小さくなったので終了」という経験的停止条件だけには依存しない。

---

# 3. predefined symbols

現在の保護されたpredefined symbol:

| 名前            | 意味                                                                                    |
| --------------- | --------------------------------------------------------------------------------------- |
| `Pi`            | 円周率。exact transcendental constant                                                   |
| `E`             | 自然対数の底。exact transcendental constant                                             |
| `Phi`           | 黄金比。exact algebraic constant                                                        |
| `I`             | 虚数単位                                                                                |
| `True`, `False` | Boolean                                                                                 |
| `Integer`       | 整数定義域                                                                              |
| `Rational`      | 有理数定義域                                                                            |
| `Real`          | 実数定義域                                                                              |
| `Complex`       | 複素数定義域                                                                            |
| `Infinity`      | 正の拡張実無限大。exact値の`precision/accuracy`もこれを返す                              |
| `ComplexInfinity` | 方向未定の無限大。非zeroと証明済みの値をexact zeroで割ったときに返す                  |
| `Indeterminate` | 数値を一意に定義できない式を表す，保護された非数値結果                                  |

旧版の`Tau`, `NA`, `ESP`は現在predefined constantではない。

---

# 4. 入力構文

## 4.1 function call

函数呼出には**角括弧 `[]` だけ**を使用する。丸括弧 `()` はgrouping専用であり，函数呼出delimiterにはしない。

```text
sin[Pi/6]
sqrt[2]
f[x]
```

既知の函数名に `sin(x)` のような旧丸括弧構文を使うとSyntaxError。通常identifierの `x(x+1)` は暗黙乗算として受理するが，Formatterは `x*(x+1)` と明示する。グルーピングは丸括弧を使い，単独の`[x+1]`はgroupではない。ユーザー函数の定義も `f[x] := ...` の形式に限定する。

## 4.2 配列

```text
{1,2,3}
{{1,2},{3,4}}
```

内部ではdense `ArrayExpr`として扱う。numeric値をimmutableなpacked pageへ保持し，shape / offset / stridesを別に持つため，transposeや一部reshape/sliceはbackingを共有できる。これは内部最適化であり，ユーザーからは通常のArrayとして見える。

## 4.3 変数・ユーザー函数

```text
x := 3
-> 3

f[t] := t^2 + 1
f[4]
-> 17
```

変数への`:=`は右辺を評価して保存する`Set`相当。
函数定義は本体を保持する`SetDelayed`相当。

既存定義を別内容で上書きすると，評価結果とは別にInfo diagnosticを返す。

```text
x := 2
-> 2

x := 4
INFO: x redefined (was 2)
-> 4
```

函数の同一arity再定義でも同様に以前の定義を表示する。同じ値を再代入しただけならInfoは出さない。

組込み名・predefined symbolは再定義できない。

現在の定義確認・削除は次を使う。

```text
Defs[]
UnDef[x]
UnDef[x,y,f]
```

`Defs[]`は現在のglobal user variableとuser function definitionを式の配列として返す。`UnDef[...]`は指定した名前について変数定義と函数定義を削除し，実際に変更した名前数を返す。

## 4.4 履歴

短縮記法は次のとおり。

```text
@       // In [-1]
@@      // In [-2]
@@@     // In [-3]
%       // Out[-1]
%%      // Out[-2]
%%%     // Out[-3]
```

正式な参照は`In [n]` / `Out[n]`を使う。`n > 0`は画面上の絶対入力番号，`n < 0`は相対参照，`n = 0`はTypeErrorとする。

```text
In [1]
Out[1]
In [-1]
Out[-1]
```

`In [n]`は対象入力のlowered Exprを取得した後，**現在のsession環境で通常評価する**。正の`In [n]`は絶対入力番号，負の`In [-n]`は現在評価中の入力slotを除いて過去の入力slotを数える。したがって`@` / `@@` / `@@@` / ... はそれぞれ`In [-1]` / `In [-2]` / `In [-3]` / ... を意味する。連続個数に固定上限はない。評価エラーになった入力でもparse/lowerまで成功していれば`In [-1]`で再評価できる。一方，Lexer / Parser / Lowererで式として成立しなかった入力は履歴へcommitせず，`In[n]`番号も消費しない。エラー表示にはその時点のpending入力番号を使うため，修正後の次入力は同じ`In[n]`から再開する。

`Out[n]`は保存済み結果snapshotを返し，再評価しない。正の`Out[n]`は絶対入力番号に対応するsnapshot，負の`Out[-n]`は**成功した出力だけ**を直前から数える。このため評価失敗を挟んでも常に`% == Out[-1]`，`%% == Out[-2]`，`%%% == Out[-3]`，... となる。`%`も連続個数に固定上限はない。
`D` / `integrate` / `limit`等のheld symbolic operator内でも，明示的な履歴参照はhistory参照として先に解決する。`Out[n]` / `%`は保存済みsnapshotを式へ貼り戻し，`In[n]` / `@`は保存入力を現在のsession環境で再評価する。したがって`integrate[Out[1],x]`や`integrate[2*In[1],x]`を未知函数`Out` / `In`として扱わない。

```text
In [1]> fft[{1,2,3}]
Out[1]> {6, ...}
In [2]> N[@,30]
Out[2]> {6, -1.50+0.866025403784...I, ...}
```

上の`N[@,30]`では，`Out[1]`を後から30桁化するのではなく，`In [1]`の`fft[...]`を30桁のapproximation contextで再評価できる。

絶対参照の例：

```text
In [1]> 1+1
Out[1]> 2
In [2]> 2+2
Out[2]> 4
In [3]> In [1]+In [2]
Out[3]> 6
In [4]> In [3]
Out[4]> 6
```

したがって`In [n]`は「過去入力を現在環境へ貼り戻して再実行する」意味である。過去入力が変数参照・代入・乱数等を含めば現在の定義やRNG stateを使う。生の入力ASTを表示する用途とは分離する。正の添字で現在評価中の入力自身を参照することは禁止し，自己再帰によるstack overflowを防ぐ。

評価エラーになった入力でもparse/lowerまで成功していれば絶対入力slot自体は残るが，対応する正の`Out[n]`は存在しない。

## 4.5 比較

```text
<  <=  >  >=  ==  !=
```

確定できれば`True/False`，symbolicに未確定ならPredicate式を保持する。

## 4.6 演算子優先順位

実装上の重要点:

1. postfix `!`
2. power `^` — 右結合
3. unary `+ -`
4. `* /` と暗黙乗算 — 同じtermレベルで左から処理
5. `+ -`
6. comparison
7. assignment `:=`

```text
2^3^2
-> 512

-2^2
-> -4
```

## 4.7 暗黙乗算

```text
2Pi
2(x+1)
(x+1)(x-1)
2 x
2exp[x]
2E^x
```

函数名と数値の直接結合は函数呼出とは解釈せず，暗黙乗算として扱う。数値直後の`e/E`は，その後ろに指数の数字が実際に続く場合だけ科学表記へ取り込む。

## 4.8 基数付き数値literal

現在のLowererは次を扱う。

```text
0b1010
0o17
0xFF
2#1010
16#FF
```

`base#digits`形式ではRational literalもparse可能な範囲で扱う。
`&` / `|` / `<<` / `>>`等の字句上のinfix bit演算子は未実装である。bit操作の函数APIは後述の`bitAnd` / `bitOr` / `bitXor` / `bitNot` / `bitShiftLeft` / `bitShiftRight`等を使う。

---

# 5. 条件分岐・数学的場合分け

```text
if[condition,trueExpr,falseExpr]
cases[value1 if condition1; value2 if condition2; ...]
```

`if[...]`は評価制御であり，conditionを先に評価して選択されたbranchだけを評価する。したがって非選択branch内のDomainErrorや乱数消費は発生しない。

`cases[...]`はscalar数学式としてのpiecewise表現である。各conditionがFalseと証明されたbranchは値を評価せず除去し，Trueと証明された分岐へ到達すればその値へ畳む。未確定conditionはBooleanへ捏造せず式中に保持する。`simplify`はassumptionで分岐を選択でき，`N`はpredicateを数値化せずbranch値だけを近似する。`integrate` / `limit`はconditionが微積分変数へ依存しない場合だけbranchごとに分配する。`D`はそれに加え，variable-dependentな不等式分岐のopen interiorを安全に微分できる場合は内部だけを微分し，閉境界点は微分可能性を別途証明できない限り未評価`D[...]`として残す。

```text
cases[1/x if x!=0; 0 if x==0]
-> cases[1/x if x != 0; 0 if x == 0]

simplify[cases[1/x if x!=0; 0 if x==0],x!=0]
-> 1/x
```

---

# 6. 角度仕様

**既定はRadian。**

```text
sin[Pi/6]
-> 1/2

asin[1/2]
-> Pi/6
```

明示単位:

```text
sin[30 Deg]
-> 1/2

sin[Pi/6 Rad]
-> 1/2

sin[100 Grad]
-> 1
```

単位接尾辞は数値literalだけでなく一般式へ適用できる。

```text
x Deg
Pi/6 Rad
(2x + 1) Grad
```

session既定角度は`angleMode`で確認・変更する。旧`:angle` commandは復活させない。

```text
angleMode[]
-> Rad

angleMode[Deg]
-> Deg

angleMode[Grad]
-> Grad
```

`angleMode`は表示設定ではなく数学的評価状態である。明示単位`Deg/Rad/Grad`は常にsession既定より優先する。Kernel APIの`setDefaultAngleUnit()`も同じ状態を操作する。

角度変換:

```text
DtoR[180] -> Pi
DtoG[90]  -> 100
RtoD[Pi]  -> 180
RtoG[Pi]  -> 200
GtoD[200] -> 180
GtoR[200] -> Pi
```

---

# 7. `N` — numerical approximation

```text
N[expr]
N[expr,p]
```

`p`は**有効10進桁数(significant decimal digits)**であり，既定は16桁。小数点以下の表示桁数ではない。固定小数表示は`:fix` / `--fix`が担当する。
Arrayへ再帰的に適用できるほか，`arg`などが返す明示角度単位では値の部分だけを近似し，単位は保持する。

`N`はprecision-aware evaluationの入口でもある。第2引数の要求精度を先に確定し，第1引数の評価中はそのprecision contextを保持する。通常builtinは従来どおりexact評価され，FFTなど明示的に対応したbuiltinだけが要求精度を受け取って直接certified 計算基盤へ降りる。したがってexact-firstの意味論を全体へ暗黙に変更しない。

whole-expressionをcertified数値として閉じられない場合でも，通常評価されるCall / Array / Listでは**数値閉包な部分だけ**を再帰的に近似する。自由symbolや未評価symbolic函数はexactのまま残す。`HoldAll` / `HoldFirst`等の評価属性を持つCallを勝手に再構築して保持規則を破らない。

```text
N[x+Pi,20]
-> 3.1415926535897932385+x

N[sin[x]+Pi,20]
-> 3.1415926535897932385+sin[x]

N[True,20]          -> True
N[Infinity,20]      -> Infinity
N[Indeterminate,20] -> Indeterminate
```

自由symbol，Boolean，`Infinity`，`ComplexInfinity`，`Indeterminate`等が意図的にexactのまま残ること自体はWarningではない。また内側の`D` / `limit` / `solve` / `rref`等が既に具体的なWarningを出した場合，外側`N`は重複した一般的な`N::unevaluated`を追加しない。数学的な値は存在するが保証付き計算が未実装な場合は`N::unsupported`としてWarningを出して式を保持し，真の定義域外と区別する。保護桁を増やしても有限精度入力の`InformationEnclosure`が狭まらず，分岐切断の側・極の除外・要求桁を確定できない場合は`N::precision`として有界回数で未評価へ戻す。したがってexactな`gamma[0]`や`2F1`の非正整数`c`はDomainErrorだが，`gamma[N[0,5]]`や`hypergeometric2F1[1,2,N[0,5],2]`は0から外れている可能性を入力情報が残すため`N::precision`である。同様に`log[-1+I*N[0,5]]`のように分岐切断の上下を決められない入力を一方の主値へ潰さない。complex Lambert Wで現縮小写像が保証boxを構成できない残存領域は`N::unsupported`として保持する。

```text
N[Pi,20]
-> 3.1415926535897932385
N[Phi,20]
-> 1.6180339887498948482
N[fft[{1,2,3,4}],20]
N[arg[-1],20]
-> 3.1415926535897932385 Rad

N[Pi*10^20,20]
-> 314159265358979323850
precision[N[Pi*10^20,20]]
-> 19
accuracy[N[Pi/10^20,20]]
-> 39
```

exact Rationalが有限10進になる場合，表示は必要以上に0埋めしない。例えば `N[1/2,10] -> 0.5` である。certified interval由来の有効桁結果では，要求桁に対応する末尾0の連続だけを圧縮し，最後に1個の0を残す。したがって内部の12桁保証が `1.000000000000` を確定していても表示は `1.0`，`1.500000000000` なら `1.50` とする。要求桁数，CertifiedEnclosure，InformationEnclosureはmetadataに全て保持し，表示上の0の個数を精度保証そのものとして扱わない。

## 7.1 CertifiedEnclosure / InformationEnclosure

`DecimalApproximation` / `ComplexDecimalApproximation`は，近似値ごとに2種類の区間を保持する。

- **CertifiedEnclosure** — 真値が必ず含まれることを計算基盤が証明した区間。内部guard桁により，ユーザーへ宣言した桁数より大幅に狭い場合がある。真値包含の検証と要求桁への一意丸め判定にはこちらを使う。
- **InformationEnclosure** — その近似値から後続計算で利用してよい情報量を表す区間。非zeroの`N[x,p]`で表示値`d`の10進指数を`e=floor(log10(|d|))`とすると，少なくとも`d ± 0.5*10^(e-p+1)`とCertifiedEnclosureの双方を包含する。したがって情報量は値のscaleに追従し，内部guard桁をユーザー可視のAccuracyとして後から回収しない。zero-centered approximationでは相対Precisionではなく，InformationEnclosureが直接absolute Accuracyを表す。

常に次を不変条件とする。

```text
CertifiedEnclosure ⊆ InformationEnclosure
```

`Infinity`は正の拡張実無限大，`ComplexInfinity`は実・複素方向を決められない無限大，`Indeterminate`は数値を一意に定義できないことを表す。いずれも有限代数symbolではなく保護atomである。現在は次の証明安全な不定形をcanonical化する。

| 形 | 結果 |
| --- | --- |
| `0/0`，`Infinity/Infinity`，`Infinity-Infinity`，`0*Infinity` | `Indeterminate` |
| `0^0`，`1^Infinity`，`Infinity^0`，`(-1)^Infinity`，`0^I` | `Indeterminate` |
| `abs[z]==1`を証明できるときの`z^Infinity` | `Indeterminate` |
| `a`の非zero性を証明できるときの`a/0` | `ComplexInfinity` |
| 任意の`x`に対する`x^(1/0)` | `Indeterminate` |

exact numeric baseの単位絶対値は実部・虚部のexact Rational演算で判定する。symbolic baseでは`simplify[z^Infinity,abs[z]==1]`のように明示的に証明できる場合だけ使い，単位円近傍のapproximationからexactな`abs[z]==1`を推測しない。`Indeterminate`は算術と登録済みscalar数学函数を伝播し，`N[Indeterminate,p]`でも保持され，`Indeterminate==Indeterminate`は`False`である。これは方向付き無限大を含む完全な代数ではなく，境界を明示したexceptional-value契約である。

InformationEnclosureは確率分布や統計的confidence intervalではない。また「真値がこの広い区間のどこにでもあり得る」と計算基盤が主張するものでもない。真値保証そのものはCertifiedEnclosureが担当し，InformationEnclosureは**現在の値から利用してよい情報量の契約**を表す。したがって計算基盤がより狭いCertifiedEnclosureを内部に持っていても，それだけを理由に既存近似値の情報量は増えない。

後続演算が近似値をexact zero / nonzeroと断定する場合，分岐切断の側を選ぶ場合，極を除外する場合，行列pivotや階数を確定する場合など，**意味論を離散的に決める判定にはInformationEnclosureだけを使う**。CertifiedEnclosureが内部で点としての零や一方のbranch側へ縮退していても，InformationEnclosureがその判定境界を跨ぐなら隠れたguard桁から結果を確定しない。したがって`N[0,5]^0`や`1/N[0,5]`はexactな`0^0` / `1/0`へ読み替えず，有限精度の曖昧性を保持する。

通常の`+ - * /`と単項`-`では2区間を独立に伝播する。exact `Number`は両方について同じpoint intervalとして混在できる。

複素近似値では実部・虚部を単純に別々のrelative Precisionへ変換して最小値を取らない。whole-complexのInformationEnclosureと値の大きさからPrecision/Accuracyを評価し，exact-zeroと証明済みのcomponentはeffective InformationEnclosure上でpoint `{0,0}`として扱う。したがって純虚数や非常に小さい純虚数でも，無関係な実部0が全体のprecisionを0へ落とさない。

```text
precision[N[I,20]]
-> 19

precision[N[I/10^100,20]]
-> 19
accuracy[N[I/10^100,20]]
-> 119
```

`+0` / `-0` / `*1` / `/1`および単項`-`は情報を消費する演算ではないため，近似値を新しい10進表示へ再量子化せず，既存のCertified/Information metadataをそのまま保持または符号反転する。exactな10の冪によるscale変更も相対Precisionを不要に1桁失わない。zero-centered結果はrelative PrecisionではなくInformationEnclosureのabsolute Accuracyから表示量子を選ぶ。

```text
precision[N[Pi,20]*10]
-> 19
accuracy[N[Pi,20]*10]
-> 18

precision[sin[N[Pi,20]]]
-> 0
accuracy[sin[N[Pi,20]]]
-> 19
```

演算後に保証precisionを`q`桁と判定した場合，丸め量子は最大`q+1` significant digitsの表示まで許す。これは「表示桁を増やして情報を発明する」という意味ではなく，InformationEnclosureが既に保証する半量子境界を余計に1 decade粗くしないための表示規則である。

`DecimalApproximation` / `ComplexDecimalApproximation`はcertified numerical evaluatorの第一級leafとして扱う。`sin` / `exp` / `log` / `sqrt` / 双曲線・逆函数・`gamma` / `erf` / `Ei` / `Si` / `Ci`に加え，対応領域の`FresnelC/S`，`1F1`，`2F1`，`zeta`，`polylog`等でも`ComplexInterval`上に両enclosureを独立伝播する。`log2` / `log10` / `fract`のようにprimitiveへrewriteされる函数もrewrite後に同じ経路へ入る。ordered comparison，`min` / `max`等の離散的判定は**InformationEnclosureだけで結論を証明できる場合**に限って確定し，内部guard桁をBoolean結果から漏らさない。exact Rational パラメータだけを受ける現行`1F1` / `2F1` / elliptic / `polylog`の一部計算基盤等は，approximate パラメータへ無理に拡張せず未評価に留める。

```text
N[Pi,20] + 1/3

Certified:    C(Pi) + {1/3}
Information:  I(Pi) + {1/3}
```

出力10進値の正当性はCertifiedEnclosureから決定し，出力として宣言できる桁数はInformationEnclosureを越えない範囲へ制限する。scale拡大や近接減算ではInformationEnclosureも演算されるため，`accuracy` / `precision`は自然に低下し得る。演算結果自身にも伝播後のInformationEnclosureを保存するので，複数回の演算を跨いでも単なる「要求桁数」へ情報を圧縮し直さない。

外側の`N`はInformationEnclosureを狭めて情報を発明しない。したがって

```text
N[N[Pi,20],100]
-> 3.1415926535897932385
```

は元の20桁保証を保持する。一方，より低い桁を要求した場合は表示丸めに対応するInformationEnclosureを追加して安全に情報を捨てられる。これらは`double`等のmachine arithmeticへ変換せず，両enclosureを`RealInterval` / `ComplexInterval`へ持ち上げて外向き丸めで計算する。

---

# 8. precision / accuracy / rationalize

## 8.1 `accuracy[x]`

`DecimalApproximation`について，真値に対する**保証可能な絶対10進桁数の整数下限**を返す。

```text
accuracy[N[1/3,20]]
-> 20
```

表示値`d`とInformationEnclosure `[iL,iU]`から

```text
max(|d-iL|, |d-iU|)
```

をabsolute error boundとして使う。`N[...,p]`生成時のInformationEnclosureには有効桁丸めに対応するscale依存の半量子が既に含まれるため，有限小数がCertifiedEnclosure上で真値と偶然完全一致していても，要求桁を越えた隠れた保護桁の情報をAccuracyとして回収しない。

exactな数・exact symbolic expressionは`Infinity`を返す。

```text
accuracy[1/3] -> Infinity
accuracy[Pi]  -> Infinity
```

## 8.2 `precision[x]`

同じabsolute error boundを，InformationEnclosureから得られる値絶対値の正の下限で割り，**保証可能な相対10進桁数の整数下限**を返す。

```text
precision[N[1/3,20]]
-> 19
```

これは要求した20有効桁を機械的に返す函数ではない。`1/3`近傍では有効20桁の丸め量子が`10^-20`なので，InformationEnclosureの相対的不確かさから保証できる整数桁数は19になる。InformationEnclosureが0を含む場合は値絶対値の正の下限を得られないため0を返す。exact expressionは`Infinity`。近接減算ではabsolute Accuracyを多く残したまま結果scaleだけが小さくなるため，Precisionだけが大きく落ちることがある。

## 8.3 `rationalize[x]`

近似値が持つInformationEnclosure内から，**分母が最小になるexact Rational**を求める。探索はexact Rational上のcontinued-fraction型interval recursionで行い，doubleへ変換しない。これにより，CertifiedEnclosureだけが保持している隠れた保護桁やhidden exact pointから，ユーザーへ宣言していない情報を`rationalize`で掘り返さない。

```text
rationalize[N[1/3,20]]
-> 1/3
```

`rationalize[x,tol]`は表示値を中心とする`[x-tol,x+tol]`内から最小分母Rationalを選ぶ。`tol`は非負exact real。

```text
rationalize[N[Pi,20],1/1000]
-> 201/64
```

`201/64`は`Pi`から0.001以内で，`355/113`より小さい分母を持つため，この仕様ではこちらが正しい。

`tol=0`は表示された有限10進値そのものをexact Rationalへ戻す。

```text
rationalize[N[1/3,20],0]
-> 33333333333333333333/100000000000000000000
```

mmCalではソースの`0.1`自体が最初から`1/10`なので，`rationalize[0.1]`は単に`1/10`のままである。Arrayや式内部のDecimalApproximationも再帰的にRational化する。

## 8.4 `explain[value]`

評価済みの値が**既に保持している情報だけ**を構造化して返す軽量introspection函数。対象を改めて`det` / `matrixRank` / LU / Eigen等へ掛けたり，Generic Arrayを全走査して性質を推論したりはしない。したがって`explain`自体のために高価な数学計算は開始しない。

```text
explain[{{1,2},{3,4}}]
-> {{"Kind","Array"},
    {"Domain","Integer"},
    {"Exactness","Exact"},
    {"ArrayRank",2},
    {"Dimensions",{2,2}},
    {"ElementCount",4},
    {"Rectangular",True},
    {"Empty",False},
    {"Matrix",True},
    {"Square",True},
    {"Order",2}}
```

通常の引数評価を先に行うため，`explain[1+2]`は`3`を説明する。履歴結果もそのまま対象にできる。

```text
explain[Out[25]]
explain[%]
```

返値はbrace構文で表示されるproperty/value pair列。`Dimensions`や各種Enclosureの値自体がArray/braceになり得るため，内部表現はdense `ArrayExpr`ではなく一般`ListExpr`を用いる。表示上は`{{"Kind",...},{...}}`であり，情報を文字列へ潰さない。

`Exactness`はbooleanではなく分類値を返す。現在の主な値は`"Exact"` / `"CertifiedApproximation"` / `"Unknown"`。近似値を単に`Exact=False`とは表現しない。

組込み数学定数・予約symbolは一般の未知symbolへ落とさない。`Pi/E/Phi`はMathRegistryに既に登録された数学metadataをO(1)で参照し，`Infinity`等はSymbolRegistryの予約意味を用いる。

```text
explain[Pi]
-> {{"Kind","Constant"},
    {"Domain","Real"},
    {"Exactness","Exact"},
    {"Name","Pi"},
    {"Real",True},
    {"Positive",True},
    {"Irrational",True},
    {"ArithmeticClass","Transcendental"}}

explain[Infinity]
-> {{"Kind","Constant"},
    {"Domain","ExtendedReal"},
    {"Exactness","Exact"},
    {"Name","Infinity"},
    {"Infinite",True},
    {"Finite",False},
    {"Sign","Positive"}}
```

`I`は通常評価でexact complex `Number`へloweringされるため，`explain[I]`は入力tokenではなく評価後の複素数値を説明する。

組込み函数symbolもBuiltinRegistry / MathRegistryの既存metadataだけをO(1)で参照して説明する。arity，held-argument規則，数学函数の定義域 / parity / period / 主値 inverse / real 値域等を確認できる。

```text
explain[sin]
-> {{"Kind","BuiltinFunction"},
    {"Domain","Function"},
    {"Exactness","Exact"},
    {"Name","sin"},
    {"Arity",1},
    {"ArgumentEvaluation","All"},
    {"FunctionDomain","ComplexToComplexRealPreserving"},
    {"Parity","Odd"},
    {"PeriodTurns",1},
    {"PrincipalInverse","asin"}, ...}

explain[table]
-> ... {"ArgumentEvaluation","HoldFirstAndTableIteratorSpec"} ...
```

これは函数を実行して性質を調べる機構ではなく，registryに既に登録済みのmetadataの照会である。

certified decimal approximationでは，要求有効桁数に加えてCertifiedEnclosureとInformationEnclosureを別々に確認できる。前者は真値保証，後者は後続計算で利用してよい情報量である。

```text
explain[N[Pi,20]]
-> {{"Kind","DecimalApproximation"},
    {"Domain","Real"},
    {"Exactness","CertifiedApproximation"},
    {"RequestedPrecisionDigits",20},
    ...
    {"CertifiedEnclosure",{certifiedLower,certifiedUpper}},
    {"InformationEnclosure",{informationLower,informationUpper}}}
```

Arrayではstorage metadataからO(1)で分かる`Domain` / `Exactness`と，shapeからほぼ無料で分かる`ArrayRank` / `Dimensions` / `ElementCount` / `Vector` / `Matrix` / `Square` / `Order` / `Empty`を返す。`det`，数学的`matrixRank`，invertibility，eigenvalue等は返さない。必要なら既存函数を明示的に呼ぶ。

整数では符号・zero・`BitLength`を返す。10進桁数は巨大整数で10進変換を必要とするため，軽量性を優先して自動計算しない。Rationalでは分子・分母のbit長を返す。

```text
explain[value,"internal"]
```

は開発・性能調査用で，`Representation`，Arrayの`Storage` / `Contiguous` / `StoredExpressions`，近似値の`ApproximationOrigin`等を追加する。**`"internal"`のproperty名・値は互換性保証対象ではない。** 未知modeはerrorとし，将来高コストな`"full"`相当を暗黙に実行しない。

---

# 9. 基本演算・代数

通常の演算子:

```text
+  -  *  /  ^  !
```

`pow[x,y]`は`Power[x,y]`の入力用alias，`fact[x]`はfactorialのaliasである。

主なexact simplification:

```text
sqrt[8]
-> 2 sqrt[2]

sqrt[2/3]
-> sqrt[6] / 3

sqrt[-8]
-> 2I sqrt[2]

sqrt[z]^2
-> z
```

ただし主値分岐を壊す変形はしない。

```text
sqrt[x^2]
-> sqrt[x ^ 2]       // xが何者か不明

simplify[sqrt[x^2], element[x,Real]]
-> abs[x]

simplify[sqrt[x^2], x >= 0]
-> x
```

---

# 10. 基本数学・複素函数

| 函数          | 概要                    | 例                               |
| ------------- | ----------------------- | -------------------------------- |
| `sqrt[x]`     | 主値平方根              | `sqrt[-4] -> 2I`                 |
| `cbrt[x]`     | 実立方根。実数定義域    | `cbrt[-8] -> -2`                 |
| `abs[z]`      | 絶対値・複素絶対値      | `abs[3+4I] -> 5`                 |
| `sign[z]`     | 実符号 / 複素`z/abs[z]` | `sign[3+4I] -> 3/5+4/5I`         |
| `re[z]`       | 実部                    | `re[3+4I] -> 3`                  |
| `im[z]`       | 虚部                    | `im[3+4I] -> 4`                  |
| `conj[z]`     | 複素共役                | `conj[3+4I] -> 3-4I`             |
| `arg[z]`      | 主値偏角                | `arg[-1] -> Pi Rad`              |
| `hypot[x,y]`  | exact `sqrt[x^2+y^2]`   | `hypot[3,4] -> 5`                |
| `cis[x]`      | `cos[x]+I sin[x]`       | `cis[Pi/3] -> 1/2 + I sqrt[3]/2` |
| `polar[r,t]`  | `r cis[t]`              | `polar[2,Pi/3]`                  |
| `nextpow2[x]` | 最小nで`2^n >= x`       | `nextpow2[9] -> 4`               |

互換alias:

```text
real -> re
imag -> im
mag  -> abs
unit,csgn -> sign
rect -> polar
```

---

# 11. 指数・対数

| 函数       | 仕様                        |
| ---------- | --------------------------- |
| `exp[x]`   | 複素平面全域で正則な指数函数 |
| `log[x]`   | 主値 natural logarithm |
| `log[b,x]` | 主値 `Log[x]/Log[b]`   |
| `log2[x]`  | `log[2,x]`の短縮形            |
| `log10[x]` | `log[10,x]`の短縮形           |
| `expm1[x]` | `exp[x]-1`を0近傍で安定評価 |
| `log1p[x]` | `log[1+x]`を0近傍で安定評価 |

`expm1` / `log1p`は，極小入力で`exp[x]-1`や`1+x`を要求精度のまま作って桁を失わない。入力の2進桁位置から必要な作業精度を増やし，最後に要求精度へ丸める。

例:

```text
exp[1]
-> E

log[E]
-> 1

log[-1]
-> I Pi

log[10,1000]
-> 3

log[1,10]
-> DomainError       // 底1
```

`log[0]`はInfinityへ置換せずDomainError。

---

# 12. 三角函数

実装済み:

```text
sin cos tan cot sec csc
asin acos atan atan2
```

既定Radなので:

```text
sin[Pi/6] -> 1/2
cos[Pi/3] -> 1/2
tan[Pi/4] -> 1
asin[1/2] -> Pi/6
atan2[1,-1] -> 3Pi/4
```

`tan`, `sec`, `cot`, `csc`の極はdefinednessとして扱い，有限値を捏造しない。

---

# 13. 双曲線函数

実装済み:

```text
sinh cosh tanh
asinh acosh atanh
csch sech coth
```

逆双曲線函数は主値の複素分岐を持ち，`MathRegistry`に分岐情報を保持する。

---

# 14. Cardinal・安定初等函数

```text
sinc[x]
cosc[x]
tanc[x]
sinhc[x]
tanhc[x]
expc[x]
```

可除特異点はexactに連続延長する。

```text
sinc[0]  -> 1
cosc[0]  -> 0
tanc[0]  -> 1
sinhc[0] -> 1
tanhc[0] -> 1
expc[0]  -> 1
```

三角cardinal函数は角度表現をRadian量へ正規化してから比を取る。有限precision入力のInformationEnclosureが0を跨ぐ場合でも，`sinc/cosc/sinhc/expc`は0近傍のTaylor多項式と剰余上界，`tanc/tanhc`は安全な連続延長形を用い，可除特異点を理由に未評価へ落とさない。

```text
sinc[Pi/2]
-> 2 / Pi

sinc[90 Deg]
-> 2 / Pi
```

---

# 15. 丸め・整数補助函数

```text
floor ceil trunc round frac
gcd lcm mod rem quotient
bitand bitor bitxor bitnot
bitshiftl bitshiftr bitlength bitcount bitget
fma clamp proj
```

例:

```text
floor[-3/2] -> -2
ceil[-3/2]  -> -1
trunc[-3/2] -> -1
round[5/2]  -> 2       // nearest-even
round[125,-1] -> 120
round[135,-1] -> 140
frac[-3/2]  -> 1/2

gcd[84,126,210] -> 42
lcm[6,8,9] -> 72

quotient[-5,3] -> -1
rem[-5,3]      -> -2
mod[-5,3]      -> 1

bitand[-1,5] -> 5
bitor[-8,3]  -> -5
bitxor[-1,5] -> -6
bitnot[5]    -> -6
bitshiftr[-3,1] -> -2
bitget[-2,100]  -> 1
```

`round[x,n]`は`10^-n`単位のnearest-even丸めで，負の`n`も許す。approximation入力ではInformationEnclosure全体が同じ丸め値へ入る場合だけ確定する。

bitwise函数は任意長BigIntに対して**無限長2の補数**として定義する。したがって負数は上位bitを1でsign-extensionする。`bitshiftr`はuser-facingには算術右shift，`bitcount`は無限個の1を持つ負数には定義せず非負整数だけを受ける。現段階ではlexerへ`& | << >>`等のinfix構文を追加せず，函数APIへ意味論を一意に集約している。

`fma[a,b,c]`はexact入力ではexactに`a*b+c`を返し，certified approximationを含む場合は中間DecimalApproximationへ一度丸めず一つのCertified/Information評価へ入れる。`clamp[x,lo,hi]`は実数用で，approximationではInformationEnclosureから境界順序を保証できる範囲だけ確定する。`proj[z]`は現在の有限exact/certified複素値に対して恒等写像であり，完全な複素Infinity/Riemann球面意味論はextended-real体系と同時に扱う。

`mod`はfloor quotient，`rem`はtruncate-toward-zero quotientに対応する。

---

# 16. 組合せ・軽量数論

```text
perm[n,r]
comb[n,r]
fib[n]
```

例:

```text
perm[10,3] -> 720
comb[10,3] -> 120
fib[100]   -> 354224848179261915075
```

`fib`はfast doubling O(log n)。
`comb`は対称性`r=min[r,n-r]`を利用。

軽量数論として次も実装している。

```text
isprime[n]
nextprime[n]
prevprime[n]
factorint[n]
totient[n]
```

```text
isprime[97]    -> True
nextprime[14]  -> 17
prevprime[14]  -> 13
factorint[360] -> {2, 2, 2, 3, 3, 5}
totient[9]     -> 6
```

`isprime`は`0 <= n <= 2^64-1`でdeterministicなstrong Miller-Rabin判定を行う。`factorint` / `totient`は同じ`uint64`範囲でdeterministic Pollard-Rho + exact prime verificationを使う。現在の証明計算基盤を超えるBigIntではprobable-primeを`True`として返さず，未評価を保持する。`factorint[-n]`は先頭に`-1`を付けたflat prime factor listを返し，`factorint[0]`はDomainErrorである。

---

# 17. 特殊函数

## 17.1 Gamma / LogGamma

```text
gamma[5]
-> 24

gamma[1/2]
-> sqrt[Pi]

gamma[-1/2]
-> -2 sqrt[Pi]

N[gamma[1/3],20]
-> 2.6789385347077476337
```

一般実数はStirling–Bernoulli + rigorous remainder，負実数はreflectionを使用。
非正整数極はDomainError。exact complex入力の`N`は右半平面へrecurrence shiftした後，`ComplexInterval`上のStirling展開と明示的剰余上界でcertifyする。

```text
N[gamma[1+I],20]
-> 0.49801566811835604271-0.15494982830181068512I
```

`lgamma[x]`は現在 **実軸上の`log[abs[gamma[x]]]`**。複素`LogGamma`とは分離している。

## 17.2 Erf

```text
erf[x]
erfc[x]
```

```text
erf[0]  -> 0
erfc[0] -> 1
N[erf[1],20] -> 0.84270079294971486934
N[erf[1+I],20]
-> 1.3161512816979476449+0.19045346923783468628I
```

複素入力は整級数を`ComplexInterval`上で評価し，明示的な剰余上界で要求桁を保証する。

## 17.3 Beta

現在は正実数域に限定。

```text
beta[2,3] -> 1/12
beta[1/2,1/2] -> Pi
betaln[1/2,1/2] -> log[Pi]
```

一般Gamma比へ無条件展開して極の相殺を壊さない。

## 17.4 Zeta / Digamma / Trigamma / regularized incomplete Beta

```text
zeta[s]
digamma[x]
trigamma[x]
ibeta[a,b,x]
```

`zeta`はRiemann zeta函数である。exactには代表値と自明零点を扱う。certified `N` 計算基盤は`s=1`を除く有限複素平面を対象とし，`Re[s]>=0`では補正次数ごとに剰余条件を検査するEuler-Maclaurin展開，左半平面ではfunctional equationを用いる。`s=1`だけが極であり，finite-precision入力区間が極を跨ぐ場合は`N::precision`とする。

```text
zeta[0]  -> -1/2
zeta[-2] -> 0
zeta[2]  -> Pi^2/6
N[zeta[3],20] -> 1.2020569031595942854
N[zeta[2+I],20]
-> 1.1503557032549026717-0.43753086591960788112I
N[zeta[1/2],20] -> -1.4603545088095868129
N[zeta[1/2+I],20]
-> 0.14393642707718906032-0.72209974353167308913I
```

`digamma[x]`は`d/dx lgamma[x]`，`trigamma[x]`は`d/dx digamma[x]`である。非正整数極はDomainError。現certified 計算基盤は正実数に加えてcomplex入力をrecurrenceで右半平面へ移し，Bernoulli漸近展開と明示的剰余上界を持つComplexIntervalで評価する。`trigamma[n]`の正整数値は`Pi^2/6`と有限二次調和和へexact還元する。

```text
N[digamma[1],20]  -> -0.57721566490153286061
N[trigamma[1],20] -> 1.6449340668482264365
N[digamma[1+I],30]  -> 0.0946503206224769772718784827219+1.07667404746858117413405079475I
N[trigamma[1+I],30] -> 0.463000096622763786298326518184-0.794233542759318865583013617157I
trigamma[2]        -> Pi^2/6-1
D[gamma[x],x]      -> digamma[x]gamma[x]
D[lgamma[x],x]     -> digamma[x]
D[digamma[x],x]    -> trigamma[x]
```

`ibeta[a,b,x]`は**正則化不完全Beta函数**`I_x(a,b)`である。real contractは`a>0`, `b>0`, `0<=x<=1`。正整数`a,b`とexact Rational `x`は有限binomial sumへexact還元する。certified `N`ではfinite-precisionを含む正の実数パラメータ `a,b`とreal `x`を区間のまま伝播し，`I_x(a,b)`が`a`に関して減少，`b`と`x`に関して増加することを利用してパラメータ enclosure全体を端点評価で包含する。

```text
ibeta[1,1,1/4] -> 1/4
ibeta[2,3,1/2] -> 11/16
N[ibeta[1/3,2/3,1/4],20] -> 0.53302858123542523627
N[ibeta[N[1/3,8],N[2/3,8],1/4],8] -> 0.53302858
```

polygamma高階と複素パラメータを持つ`ibeta`は引き続き未実装である。

## 17.5 generalized factorial系

```text
binom[x,n]
fallingfact[x,n]
risingfact[x,n]
```

現段階では非負整数次数へexact finite productを構成できる範囲が中心。

```text
binom[1/2,2] -> -1/8
fallingfact[5,3] -> 60
risingfact[5,3] -> 210
```

## 17.6 Fresnel C / S

mmCalでは標準Fresnel積分を

```text
fresnelc[x] = integral_0^x cos[Pi t^2/2] dt
fresnels[x] = integral_0^x sin[Pi t^2/2] dt
```

に対応する整な奇函数として扱う。exactに閉じない引数は記号式を保持し，`N`では実数に加えて複素入力も`ComplexInterval`上でcertifyする。複素経路に固定の`|z|`境界は設けず，中程度の引数や対角方向は整函数のMaclaurin級数，軸近傍の大引数はDLMF 7.12の`f/g`漸近展開と最初の未使用項による剰余保証を使う。90度回転対称性`C[i z]=i C[z]`，`S[i z]=-i S[z]`で任意の軸近傍を同じ保証wedgeへ写し，漸近証明が閉じない領域はMaclaurin級数へ切り替える。停止はterm capと共通`EvaluationBudget`で制御する。実数経路も大引数用の保証付き漸近算法を持つ。

```text
fresnelc[0] -> 0
fresnels[0] -> 0
N[fresnelc[1],20] -> 0.77989340037682282947
N[fresnels[1],20] -> 0.43825914739035476608
N[fresnelc[1+I],20] -> 2.5557937781024390246+2.5557937781024390246I
N[fresnels[1+I],20] -> -2.0618882191948404681+2.0618882191948404681I

D[fresnelc[x],x] -> cos[Pi x^2/2 Rad]
D[fresnels[x],x] -> sin[Pi x^2/2 Rad]
```

微分の位相には`Rad`を明示する。Fresnel函数の定義自体はsessionの既定角度単位に依存しないためである。

## 17.7 合流型超幾何函数 1F1

Kummerの合流型超幾何函数を

```text
hypergeometric1F1[a,b,z]
```

で表す。`z`について整函数であり，`b = 0,-1,-2,...`には一般にパラメータの極があるため，その場合を無条件に有限値へ簡約しない。現段階のexact評価は停止する級数，`z=0`，`a=b`等の安全に閉じる場合を扱う。`N`はexact Rational パラメータ `a,b`と実・複素`z`を，整級数＋明示的な剰余上界で`ComplexInterval`上にcertifyする。パラメータ自体をapproximateへ拡張することはまだしない。

```text
hypergeometric1F1[0,3,2] -> 1
hypergeometric1F1[-2,3,2] -> 0
hypergeometric1F1[2,2,1] -> E
N[hypergeometric1F1[1/6,7/6,1],20]
-> 1.1920688079818883008
N[hypergeometric1F1[1/2,5/4,1+I],20]
-> 1.2988692086674201067+0.72862412303674683434I
```

パラメータが微分変数に依存しないとき，

```text
D[hypergeometric1F1[a,b,z],z]
= a hypergeometric1F1[a+1,b+1,z]/b
```

を使う。積分器では，上側不完全Gammaによる局所式が主値分岐や原点のremovable holeを持つ場合に，原点を含めて正則な1F1表現を優先する。例えば，

```text
integrate[exp[-x^2],x]
-> erf[x]sqrt[Pi]/2

integrate[exp[-x^2],{x,0,Infinity}]
-> sqrt[Pi]/2

integrate[exp[-x^2],{x,-Infinity,Infinity}]
-> sqrt[Pi]

integrate[exp[x^6],x]
-> x hypergeometric1F1[1/6, 7/6, x^6]
```

より一般に正整数`n`について`exp[c x^n]`を同じ系列へ還元できる。

## 17.8 Gauss超幾何函数 2F1

Gaussの超幾何函数を

```text
hypergeometric2F1[a,b,c,z]
```

で表す。`c = 0,-1,-2,...`には一般にパラメータの極があり，`z`については主値分岐を採用する。現段階のexact評価は，上側パラメータが非正整数で停止する有限級数や，`a=0` / `b=0`など安全に閉じる場合を扱う。誤差保証付き`N` 計算基盤はexact 数値パラメータを`ComplexInterval`へ持ち上げ，`|z|<1`ではGauss級数とtail boundを使う。`z=1`では`Re(c-a-b)>0`を証明できる場合にGauss summation `Gamma[c] Gamma[c-a-b]/(Gamma[c-a] Gamma[c-b])`を用いる。`|z|>1`でも，パラメータ差の退化や分岐切断危険をexactに除外でき，`|1/z|<1`を証明できる場合は主値 `1/z` connection formulaへ送る。exactな実数`z>1`は主値の分岐切断上で定めた規約値として評価する。一方，`2+I*N[0,p]`のように有限精度情報が分岐切断の上下を同時に許す場合は，一方の値へ潰さず`N::precision`とする。証明できない境界・退化caseは未評価を保持する。

```text
hypergeometric2F1[-2,1,3,1/2] -> 17/24
hypergeometric2F1[0,2,3,x] -> 1
N[hypergeometric2F1[1/2,1/2,3/2,1/4],20]
-> 1.0471975511965977462
N[hypergeometric2F1[1/2,1/3,5/4,1/2+I/4],20]
-> 1.0768624682230010329+0.057258816434281232216I
N[hypergeometric2F1[3.4,5.6,4+I,4.6+2I],20]
-> 0.0046136876612922014955+0.0019659119401108294965I
```

パラメータが微分変数に依存しない場合，

```text
D[hypergeometric2F1[a,b,c,z],z]
= a b hypergeometric2F1[a+1,b+1,c+1,z]/c
```

を使う。積分器では，例えば

```text
integrate[sqrt[1+2x^3],x]
-> x hypergeometric2F1[-1/2, 1/3, 4/3, -2x^3]

integrate[1/(1+x^5),x]
-> x hypergeometric2F1[1, 1/5, 6/5, -x^5]
```

のようなbinomial-power familyへ利用する。一般の2F1を`Solve`で逆函数化する規則は持たない。大域単射性を証明できないためであり，停止級数やexact退化で既存代数式へ落ちた場合だけ通常のSolverへ渡す。

## 17.9 不完全楕円積分 F / E / Pi

mmCalではLegendre形の不完全楕円積分を

```text
ellipticF[phi,m]
ellipticE[phi,m]
ellipticPi[n,phi,m]
```

で表す。第2引数`m`はパラメータであり，振幅`phi`は**常にRadian**として解釈する。sessionの`Deg/Rad/Grad`設定には依存しない。主値分岐を採用し，一般複素パラメータの分岐切断や極を単純な「everywhere defined」には扱わない。

```text
ellipticF[x,0] -> x
ellipticE[x,0] -> x
ellipticPi[0,x,0] -> x

N[ellipticF[1/2,1/3],20]
-> 0.5068477562654311092
N[ellipticE[1/2,1/3],20]
-> 0.49331536201475850521
N[ellipticPi[1/5,1/2,1/3],20]
-> 0.51520338216141386085
N[ellipticF[1/2,99/100],20]
-> 0.52198775871658283077
N[ellipticE[Pi/2,1],20]
-> 1.0
N[ellipticF[1/2,2],20]
-> 0.55135887907967981413
```

保証付きtailの実効収束率が十分小さい場合は従来のguard付きLegendre級数をfast pathとして維持する。そこから外れる実数領域では振幅を`Pi`周期で還元し，Legendre形をcertified Carlson symmetric integrals `RF`，`RD`，`RJ`へ変換して評価する。固定`|m|<=9/10` / `|n|<=9/10` capability boundaryは撤去済みである。`m>1`または`n>1`でも，還元後の積分路全体が分岐点 / 極より手前にあることを証明できる局所実数値は扱う。周期を跨ぐ場合はその周期内に特異点がないことを証明できる必要がある。`ellipticE`のexact `m=1`は有限な実数退化を専用処理する。一般complex elliptic 解析接続は未対応である。振幅微分は

```text
D[ellipticF[phi,m],phi]
= 1/sqrt[1-m sin[phi Rad]^2]

D[ellipticE[phi,m],phi]
= sqrt[1-m sin[phi Rad]^2]

D[ellipticPi[n,phi,m],phi]
= 1/((1-n sin[phi Rad]^2)sqrt[1-m sin[phi Rad]^2])
```

である。したがって標準kernelは直接積分できる。

```text
integrate[1/sqrt[1-(1/3)sin[x]^2],x]
-> ellipticF[x, 1/3]

integrate[sqrt[1-(1/3)sin[x]^2],x]
-> ellipticE[x, 1/3]

integrate[1/((1-(1/5)sin[x]^2)sqrt[1-(1/3)sin[x]^2]),x]
-> ellipticPi[1/5, x, 1/3]

integrate[1/sqrt[1-x^4],x]
-> ellipticF[asin[x], -1]
```

最後のquartic reductionは正しい局所primitiveだが，現在の`fullSimplify`は`sin[asin[x]]`と主値 square rootの積を一般に安全な恒等式へ潰し切れない。そのためderivative-back harnessではResolutionOnlyとして監視し，証明器不足を理由に積分能力を削らない。一般の楕円函数方程式も，逆楕円函数族をまだ持たないため`Solve`は未解決を保持する。`m=0`等でexactに通常式へ退化した場合だけ既存Solverが解く。

## 17.10 Ei / Si / Ci / li / Polylogarithm

積分で頻出する主値 special functionsを次の名前で表す。

```text
Ei[x]
Si[x]
Ci[x]
li[x]
polylog[s,z]
```

`Ei`, `Ci`, `li`, `polylog`は一般に分岐を持つため，MathRegistryでは主値分岐として扱う。`Si`は整な奇函数である。`Ei/Si/Ci`は安全な複素級数領域まで`ComplexInterval`でcertifyする。`polylog[n,z]`は正整数orderの`|z|<1`を級数＋tail boundで扱うほか，`n=2`ではDLMF 25.12.3/25.12.4/25.12.6の主値 connection formulaを分岐切断から離れた領域で用いる。正実軸の`z≈1`ではorder 3～12に対する`mu=log(z)`正整数極限展開を保証付きfast pathとして使い，finite-precision実区間もpointへ縮退させず伝播する。分岐切断や収束境界を推測値で埋めない。

```text
Si[0] -> 0
Si[-1] -> -Si[1]
li[0] -> 0
polylog[0,z] -> z/(1-z)
polylog[1,z] -> -log[1-z]
polylog[2,1] -> Pi^2/6
polylog[2,-1] -> -Pi^2/12
polylog[3,1] -> zeta[3]
polylog[3,-1] -> -3zeta[3]/4

N[Ei[1],20] -> 1.8951178163559367555
N[Si[1],20] -> 0.94608307036718301494
N[Ci[1],20] -> 0.33740392290096813466
N[li[2],20] -> 1.0451637801174927848
N[polylog[2,1/2],20] -> 0.58224052646501250590
N[polylog[2,999/1000],20] -> 1.6370226052761177427
N[polylog[3,999/1000],20] -> 1.2004153539954643452
N[polylog[2,-2],20] -> -1.4367463668836809464
N[Ei[1+I],20] -> 1.7646259855638540684+2.3877698515105224193I
N[Si[1+I],20] -> 1.1042226582355817396+0.88245380500791774338I
N[Ci[1+I],20] -> 0.88217218055593632505+0.28724913351995593953I
N[polylog[2,1/2+I/4],20]
-> 0.54586750496407962676+0.33913769923976904082I
```

保証付き特殊函数の計算には，数学的な定義域とは別に計算量上限がある。現在の`1F1`は固定の`|z|`境界を持たず，将来項比の保証上界が最大250000項の級数budget内で収束域へ入る限り，実数・複素数とも評価を試みる。`2F1`のGauss級数は数学的な収束域そのものの`|z|<1`で使い，固定の内部閾値は設けず項数上限と共通`EvaluationBudget`で停止する。`z=1, Re(c-a-b)>0`はGauss summationで評価する。実数楕円積分`F/E/Pi`も旧`9/10`パラメータ境界を持たず，保証付き剰余項比`<=9/10`の級数高速経路とCarlson `RF/RD/RJ`計算基盤を，分岐・極の証明，duplication refinement上限，共通`EvaluationBudget`に従って使い分ける。複素`Ei/Ci`にも固定絶対値境界はない。中程度の引数はguard桁付き区間級数，大引数はDLMF 6.12の`E1`漸近展開と主値接続公式を使い，負実軸の分岐切断はInformationEnclosureから判定する。漸近剰余の証明が閉じない場合は級数へ切り替え，共通`EvaluationBudget`で計算量を制限する。正整数位数`polylog`の旧`|z|<=49/50`境界も撤去済みで，`|z|<1`の級数は最大1000000項と共通`EvaluationBudget`で制限する。`Li_2`は主値接続公式により，負実軸・単位円近傍・分岐切断を避けた一部`|z|>1`へ拡張するが，exactな正実分岐切断上では上下の境界値を勝手に選ばず未評価に留める。正実`z≈1`の高位数では`mu=log(z)`整数極限展開を高速経路として使う。実数`Ei/Si/Ci`も固定絶対値境界を設けず，中程度の引数のguard桁付きTaylor展開と大引数の保証付き漸近展開を使い分け，剰余証明・項数上限・共通`EvaluationBudget`で停止する。複素`2F1`は`|z|>1`を証明でき，パラメータ退化と分岐切断を安全に除外できる場合だけ主値の`1/z`接続公式を使う。

これらの閾値外で値が数学的に存在することと，現計算基盤がcertifyできることは別である。固定series領域・固定term上限・固定planner budgetを超えた場合はprecisionを無制限に増やさず，`N::unsupported`と未評価式へ戻す。`zeta`については例外的にcritical stripと左半平面までcertified 解析接続を持ち，有限複素平面では`s=1`だけを極として扱う。区間幅だけが計算基盤境界や分岐切断を跨いでいる場合はguard precisionを増やして再判定するが，top-level `N`は局所16回で必ず停止する。既存の有限precision入力区間そのものが境界を跨ぎ続ける等，要求桁を証明できないまま局所上限へ達した場合は`N::precision`を出して未評価式を保持し，global resource budget exhaustionへ読み替えない。

現在の微分Knowledgeは，引数・order パラメータが微分変数に依存しない範囲で

```text
D[Ei[x],x] -> exp[x]/x
D[Si[x],x] -> sinc[x Rad]
D[Ci[x],x] -> cos[x]/x
D[li[x],x] -> 1/log[x]
D[polylog[s,x],x] -> cases[polylog[s-1, x]/x if x != 0; 1 if x == 0]
```

を使う。`polylog[2,x]`は`polylog[1,x]`を`-log[1-x]`へexact退化させるため，

```text
D[polylog[2,x],x] -> cases[-log[1-x]/x if x != 0; 1 if x == 0]
```

まで閉じる。direct variableの高階微分`D[polylog[s,x],{x,n}]`は`n<=64`で，Euler作用素`theta=x D`とsigned Stirling numberを使ってgenericなnested `D[cases[...]]`を作らず直接構成する。`x=0`では級数係数から`n!/n^s`をexactに保持する。

積分器ではこの共有Knowledgeにより，

```text
integrate[exp[x]/x,x] -> Ei[x]
integrate[sin[x]/x,x] -> Si[x]
integrate[cos[x]/x,x] -> Ci[x]
integrate[1/log[x],x] -> li[x]
integrate[li[x],x] -> x li[x]-Ei[2log[x]]
integrate[log[1-x]/x,x] -> -polylog[2, x]
```

を返す。一般の`Ei/Si/Ci/li/polylog`方程式に主値 inverseを一個だけ返す`Solve`規則は持たない。大域単射性・分岐を証明できないためであり，`polylog[0,z]`や`polylog[1,z]`のように既存の代数函数・`log`へexact退化した場合だけ既存Solverへ渡す。

## 17.11 Lambert W

```text
lambertw[z]
lambertw[k,z]
```

`lambertw`は`w exp[w] == z`を満たすLambert W函数である。1引数形は主値分岐`k=0`，2引数形は整数分岐`k`を明示する。exactな記号函数として代表exact値と微分，Real `solve`の指数方程式分類に利用するとともに，`N`では実分岐と複素数の任意整数分岐を保証付きで評価する。

```text
lambertw[0] -> 0
lambertw[E] -> 1
lambertw[-1/E] -> -1
lambertw[-1,-1/E] -> -1
D[lambertw[x],x] -> exp[-lambertw[x]]/(1+lambertw[x])
D[lambertw[x],{x,2}] -> (-2-lambertw[x])exp[-2lambertw[x]]/(1+lambertw[x])^3
```

実軸では`k=0`と`k=-1`の実分岐をSolverが必要な範囲で区別し，実値領域では単調な実逆函数計算基盤を優先する。複素数の保証付き計算基盤は主値分岐のMaclaurin展開／縮小写像と，明示分岐`k`に対応する`Log[z]+2 Pi I k-Log[w]`の包含写像を用いる。任意の整数分岐番号を受け，近傍の分岐を推測せず指定された`k`を保証する。負実軸上の複素値も分岐切断の側をexactに決められる場合はこの経路で扱い，現行の縮小写像で保証区間を構成できない残存領域だけを`N::unsupported`として保持する。

```text
N[lambertw[1],20] -> 0.5671432904097838730
N[lambertw[-1,-1/10],20] -> -3.5771520639572972184
N[lambertw[1+I],20] -> 0.65696606923043640587+0.32545033941341502999I
N[lambertw[2,1],20] -> -2.4015851048680028842+10.776299516115070898I
N[lambertw[-1/E+I/10^8],20] -> -0.99983512787429915685+0.00016485400656056139308I
N[lambertw[-1,-1/E+I/10^8],20] -> -1.0001648721257003724-0.00016489025031827418034I
```

`-1/E`近傍の複素評価では，`u=W+1`と`q=E z+1`を用いる平方根局所座標へ切り替える。`W_0`は主値平方根側，`W_-1`は負実軸の上側から接続する局所枝，`W_1`は下側の対称枝として保証する。有限precision入力が非主値分岐の分岐切断側を決められない場合は，一方の枝を推測せず`N::precision`とする。

---

# 18. Array

## 18.1 Array基盤

Arrayはrankごとに別Value型を増やさず，共通のdense `ArrayExpr`を使う。physical storageをimmutable packed page，logical layoutをshape / offset / stridesへ分離している。
公開リテラルは従来どおり `{...}` / `{{...},...}` とし，Matrix演算へ渡せるArrayは常にdense rectangularである。

```text
dimensions[A]
arrayRank[A]
length[A]
at[A,i,...]
reshape[A,{d1,d2,...}]
identity[n]
zeros[rows,cols]
rows[A]
cols[A]
diag[A]
trace[A]
```

indexは0始まり。`at`はrank未満のprefix indexも受け取り，残り次元を保持したsubarrayを返す。全rank分を指定した場合だけscalarになる。さらにfinite `SolutionSet`にも同じ0始まりindexを使える。`at[solutions,i]`は第`i` branchだけを含む`SolutionSet`を返し，branchの条件，自由変数，multiplicity，solver変数domainを保持する。`at[solutions,i,x]`は第`i` branchにおける`x`のbinding右辺だけを返す。後者は条件metadataを返さないため，条件付き解を後段へ渡す場合は2引数版を使う。`Conditional` / `Universal` / `Unresolved`集合はexplicit branch indexを一意に定義できないため対象外である。

```text
dimensions[{{1,2,3},{4,5,6}}] -> {2, 3}
arrayRank[{{1,2},{3,4}}] -> 2
at[{{1,2},{3,4}},1] -> {3, 4}
at[{{1,2},{3,4}},1,0] -> 3
at[solve[x^2==1,x],0] -> {x == 1}
at[solve[x^2==1,x],1,x] -> -1
at[solve[{x+y==3,x*y==2},{x,y}],1,y] -> 1
reshape[{1,2,3,4},{2,2}] -> {{1, 2}, {3, 4}}
```

`{...}`はユーザー意味論として一般の有限brace containerである。child shapeが全て一致する矩形値は内部でdense `ArrayExpr`へ自動昇格し，`{{1,2},{3}}`や分解結果の`{Q,R}`のようにshapeが揃わない値は一般braceのまま保持する。numeric dense Arrayは内部でInteger / Rational / Number等のpacked pageを共有する場合があるが，storage種別はユーザー意味論へ露出しない。一般brace自体は正常な値であり，Matrix函数へ渡した時点で矩形性監査が入り，非矩形ならWarningを出して未評価保持する。`dimensions` / `arrayRank`は非矩形値では全childに共通するrectangular prefixだけを返す。

```text
dimensions[{{1,2},{3}}] -> {2}
arrayRank[{{1,2},{3}}] -> 1
length[{{1,2},{3}}] -> 2
at[{{1,2},{3}},0] -> {1, 2}
transpose[{{1,2},{3}}] -> Warning + unevaluated
```

`mget[A,row,col]` は互換aliasとして `at` と同じ0始まりindexを使う。

先頭側に0長次元を持つArrayはbrace literalだけではshapeを復元できないため，Formatterは必要な場合だけ `reshape` を使う。

```text
zeros[0,3]
-> reshape[{}, {0, 3}]

dimensions[zeros[0,3]]
-> {0,3}
```

評価後にArray要素がArrayへ変わる場合も，同一shapeなら自動的に一段flattenして共通Arrayへ正規化する。scalar/Array混在またはchild shape不一致はTypeError。

# 19. 集約函数

```text
sum
prod
min
max
mean
```

scalar variadicまたは1個のArrayを受ける。

```text
sum[1,2,3] -> 6
sum[{1,2,3}] -> 6
prod[] -> 1
sum[] -> 0
mean[1,2,4] -> 7/3
```

`min/max`は順序を証明できないsymbolic値を勝手に並べない。

```text
min[x,3]
-> min[x, 3]
```

`sum[f,{k,a,b}]`型の記号有限和はまだ未実装である。有限列を明示的に生成する場合は`table`と`sum`を組み合わせられる。

## 19.1 `range` / `table` / `map`

exactな有限列の生成と明示的なelement-wise適用には次を使う。

```text
range[n]
range[a,b]
range[a,b,increment]

table[expr,{i,n}]
table[expr,{i,a,b}]
table[expr,{i,a,b,increment}]

map[f,arrayOrBrace]
```

`range`の境界とincrementはexact real Integer / Rationalに限定する。終端はincrement方向で到達範囲に含まれる場合に含める。浮動小数incrementを暗黙に丸めて列を作らない。

```text
range[5] -> {1, 2, 3, 4, 5}
range[0,1,1/3] -> {0, 1/3, 2/3, 1}
range[5,1,-2] -> {5, 3, 1}
```

`table`はbodyを保持し，iterator変数だけを各反復のlocal scopeへ束縛して通常評価する。外側に同名の定義があっても反復終了後に復元され，nested `table`も独立scopeを持つ。

```text
table[i^2,{i,5}] -> {1, 4, 9, 16, 25}
table[i/2,{i,0,2,1/2}] -> {0, 1/4, 1/2, 3/4, 1}
```

`map[f,value]`はArrayまたは一般braceの**scalar leaf**へ`f[...]`を明示的に適用し，元のshape / ragged brace構造を保つ。`exp[A]`等を自動的にelement-wise化しないため，将来のmatrix function意味論と衝突しない。

```text
map[sin,{0,Pi/2,Pi}] -> {0, 1, 0}
```

---

# 20. 記述統計

統計函数は原則 **exact real data** を受け，Rationalで閉じる量はRationalのまま返す。
1個のrank-1 Arrayまたはscalar列を受けるものが多い。

## 20.1 順序統計

```text
median
mode
quantile
percentile
iqr
percentrank
```

`quantile`はHyndman–Fan Type 7。

```text
median[1,2,3,4] -> 5/2
quantile[1/4,1,2,3,4,5,6,7] -> 5/2
iqr[1,2,3,4] -> 3/2
```

`mode`が複数ならArrayを返し，全値が1回ずつなら空Array。

## 20.2 分散・標準偏差

```text
var      // population variance
vars     // sample unbiased variance
stddev   // sqrt[var]
stddevs  // sqrt[vars]
```

```text
var[1,2,3] -> 2/3
vars[1,2,3] -> 1
stddev[1,2,3] -> sqrt[6]/3
stddevs[1,2,3] -> 1
```

## 20.3 その他

```text
geomean harmmean rms
mad madR
skew kurtp kurts
cv stderr zscore
trimmean winsor winsorR
cov corr corrspearman
```

`cov/corr/corrspearman`の推奨形:

```text
cov[{1,2,3},{2,4,6}] -> 4/3
corr[{1,2,3},{2,4,6}] -> 1
```

互換用に偶数個scalarを前半/後半へ分割する形式も受ける。

---

# 21. Assumption / 定義域

```text
element[x,Real]
element[x,Integer]
```

`simplify/fullSimplify`の第2引数にはPredicate，Array，`And`相当の条件を渡せる。

```text
simplify[sqrt[x^2], element[x,Real]]
-> abs[x]

simplify[abs[x], x >= 0]
-> x
```

矛盾したassumptionはDomainError。

`element`は所属を証明できる場合の`True`だけでなく，排他的な数学知識から`False`も返す。非整数exact Rationalは`Integer`でないこと，既知irrational/transcendental定数は`Rational` / `Integer`でないこと，minimal polynomial次数>1を証明済みのalgebraic Rootは`Rational` / `Integer`でないことを利用する。

```text
element[1/2,Integer] -> False
element[Pi,Rational] -> False
element[Phi,Rational] -> False
element[root[{-2,0,1},2],Rational] -> False
```

証明不能は`False`へ落とさず未確定のまま保持する。

---

# 22. 式変形

```text
simplify[expr]
simplify[expr,assumptions]
fullSimplify[expr]
fullSimplify[expr,assumptions]
expand[expr]
factor[expr]
collect[expr,x]
```

例:

```text
simplify[sin[x]^2 + cos[x]^2]
-> 1

fullSimplify[x^2 + 2x + 1]
-> (1 + x)^2

expand[(x+1)^3]
-> x^3 + 3 x^2 + 3 x + 1

factor[x^2-1]
-> (x-1)(x+1)
```

`fullSimplify`はbounded candidate search。短い式を選ぶために定義域を変えてよいわけではない。

```text
fullSimplify[(x^2-1)/(x-1)]
-> 元のholeを保持

fullSimplify[(x^2-1)/(x-1), x != 1]
-> 1 + x

simplify[1/x-1/x]
-> 1/x-1/x

simplify[1/x-1/x, x != 0]
-> 0

simplify[x^0]
-> x^0

simplify[x^0, x != 0]
-> 1
```

`F-F -> 0`，`0*F -> 0`，`F/F -> 1`，`F^0 -> 1`等，結果が定数となって元の式のholeを消し得る規則は，`F`が現在のassumption下で定義済みと証明できる場合だけ適用する。mmCalでは`0^0 -> Indeterminate`であるため，symbolic `x^0`も`x != 0`を証明できない限り1へ潰さない。

同じ原則は「パラメータが値から消える」特殊函数・組合せ函数にも適用する。例えば`hypergeometric2F1[0,2,3,z]`や`polylog[s,0]`，`binom[x,0]`等は通常1/0へ退化するが，消える引数自体が未定義ならそのholeを捨てない。

```text
simplify[hypergeometric2F1[0,2,3,1/x]]
-> hypergeometric2F1[0, 2, 3, 1/x]

simplify[hypergeometric2F1[0,2,3,1/x], x != 0]
-> 1
```

`Power`のdefinednessはexact Rational指数まで区別する。正の非整数Rational指数ではbase=0を許す一方，負のRational指数はbase非零を要求する。`zeta[s]`は一般に「definedness不明」とせず，唯一の極 `s=1`を` s != 1 `として扱える。

## 22.1 `series` / `normal` / `toNormal`（v1.5.5 WIP）

```text
series[expr,{x,a,n}]
series[expr,{x,a,n},assumptions]
normal[seriesExpr]
toNormal[expr]
```

`series`は点`x=a`まわりの局所展開を，内部`seriesData[...]`として保持する。`normal`はトップレベルが`SeriesData`の場合だけ剰余次数を捨て，現在保持している打切り式へ戻す。`toNormal`は式木を再帰走査し，式中に入れ子になった対応済み構造を通常式へ戻す。現在はlist / array / call内の`SeriesData`に加え，finite/conditional `SolutionSet`のbinding右辺も再帰変換し，集合構造，branch条件，自由変数，multiplicity，domain metadataは保持する。通常の式・未対応構造はそのまま保持するため，既存`normal`のトップレベル限定挙動とは分離されている。TPSA基盤で定数，展開変数，和，差，積，除算，整数冪，`exp` / `log` / `sin` / `cos` / `sinh` / `cosh`，`tan/cot/sec/csc`，`tanh/coth/sech/csch`，`expm1/log1p`，`sinc/cosc/tanc`，`sinhc/tanhc/expc`，`log2/log10`，principal `sqrt` / exact有理冪を合成し，Taylor，有限principal partを持つLaurent，およびexact有理指数格子を持つPuiseux展開へ対応する。`log2/log10`は一般の`log[base,x] = log[x]/log[base]`として有限点および`+Infinity`のlog層へ接続する。analytic函数は高階`D`の反復ではなく係数漸化式で処理する。異なるPuiseux分母はLCM格子へexactに再配置する。記号的な先頭係数を逆数化する場合は非零性を証明できるときだけ展開し，principal `log`は展開中心が正の実数または非実数であることを証明できる場合に限る。direct trigは現在の角度modeと明示`Rad` / `Deg` / `Grad`を尊重する。分岐点上の非整数有理冪は，現在正の先頭係数を持つsimple zero/pole，または既にbranchを明示したPuiseux式に限定する。`sqrt[x^2]`のような高重複零点や負向きのprincipal branchを一意に証明できない形は推測せず未評価にする。

展開中心には`Infinity`も指定できる。ここで`Infinity`は複素無限遠一般ではなく実軸の`+Infinity`を意味し，内部では`t=1/x`として`t->0+`の局所展開へ写す。したがって`seriesData[x,Infinity,...]`の指数`r`は`(1/x)^r`，log係数層は`log[1/x]^k`を表す。例えば`series[1/(x+1),{x,Infinity,4}]`は`x^-1-x^-2+x^-3-x^-4+O[x^-5]`に対応し，`normal` / `toNormal`は`(1/x)^(-m)`のような中間形を残さず通常の`x`冪へ戻す。`D`では`dt/dx=-t^2`，`integrate`では`dx=-t^-2 dt`を係数演算へ反映し，`1/x`項の積分は`-log[1/x]`としてlog層へ閉じる。ただし打切り剰余がちょうど`O(1/x)`の場合，積分後の未知剰余がlog型になり得て現在の`O(t^r)`metadataだけでは表せないため未評価へ戻す。初期対応は有理函数，多項式成長，`exp[1/x]`，Puiseux冪，`log[1/x]`等である。`sin[x]`の無限振動，`exp[x]`のessential growth，直接の`log[x]`等は別のasymptotic providerを必要とするため未評価に保つ。

実軸`+Infinity`専用のlogarithmic asymptotic providerも備える。`A(x)~c(1/x)^r`の先頭係数`c`が正と証明できる場合に限り，`log[A(x)]`を`log[c]+r log[1/x]+log[1+h]`へ分解して既存TPSA/log層へ合成する。したがって`series[log[x],{x,Infinity,n}]`，`log[2x]`，`log[x+1]`，`log[x^2+1]`，`log[sqrt[x]+1]`，さらに`log[x]^k`や`log[x]/x^m`を同じ`SeriesData`で扱える。内部log基底は引き続き`log[1/x]`だが，`normal` / `toNormal`は実軸`+Infinity`の方向情報を使って通常の`log[x]`へ戻す。`1/log[x]`のようにlogの負冪を必要とするtransseries，先頭係数が負または複素でprincipal branchを追加判断する必要がある形は未評価に保つ。`sin[x]`の無限振動や`exp[x]`のessential growthも対象外である。

```text
series[log[x+1],{x,Infinity,4}]
-> seriesData[x, Infinity, {0, 1, -1/2, 1/3, -1/4}, 0, 5, 1, {{-1, 0, 0, 0, 0}}]

normal[series[log[x+1],{x,Infinity,3}]]
-> x^(-1)-x^(-2)/2+x^(-3)/3+log[x]

series[1/log[x],{x,Infinity,3}]
-> series[1/log[x], {x, Infinity, 3}]
```

特殊函数の局所Seriesはprimitive-composition方式で接続する。函数ごとの高階導函数表を持たず，既知の一階導函数を既存TPSAで展開し，係数積分して定数項へ元函数の中心値を戻す。`erf` / `Si` / `Ei` / `Ci`，`erfc` / `fresnelc` / `fresnels`，`li`，principal `asin` / `acos` / `atan`，`lambertw`の正則中心展開を接続している。`li`は`li'(z)=1/log[z]`を用い，逆三角函数は`asin'(z)=(1-z^2)^(-1/2)`，`acos'(z)=-(1-z^2)^(-1/2)`，`atan'(z)=1/(1+z^2)`をTPSAへ通す。`erfc`は`erfc[z]=1-erf[z]`と同じGaussian kernelを符号反転して共有し，Fresnel C/SはDLMF 7.2.7–7.2.8の`cos[Pi z^2/2]` / `sin[Pi z^2/2]`をTPSAへ通す。0まわりの係数はDLMF 7.6.1，7.6.4，7.6.6および6.6.5と一致する。`SeriesData`は`log(x-center)^k`の係数層も持ち，`Ei` / `Ci`の0における対数特異点もDLMF 6.6.1 / 6.6.6の局所級数としてexactに保持する。`log[x]`自身，`x log[x]`，`log[x]^k`，Puiseux因子との積，および`x^-1 log[x]^k`の積分も同じ表現で閉じる。さらに原点で引数が`A(t)=c t^r(1+h)`，`0<r<=1`，`c>0`と証明できる場合は`log A=log c+r log t+log(1+h)`へ分解し，`log` / `Ei` / `Ci`をTaylor/Puiseux零点へ合成する。principal branchの巻き数を変え得る`r>1`，負または複素の先頭係数は推測せず未評価にする。DLMF 6.2のprincipal cut上を正則中心とする要求も従来どおり未評価にする。`li`はprincipal `li(z)=Ei(Log(z))`の局所枝を推測しないため，現段階ではDLMF 6.2.8に沿う`x>1`を証明できる実中心，または非実中心だけを正則中心として扱う。`z=0`，`z=1`，`0<z<1`の実中心，負実軸中心は未評価に保つ。`asin` / `acos`はDLMF 4.23のprincipal cutを避け，実中心では`-1<a<1`を証明できる場合，または非実中心と証明できる場合だけTaylor/Puiseux合成を行う。`atan`は実中心，実部が非零と証明できる複素中心，または虚軸上で`-1<Im(a)<1`を証明できる中心を許可し，`±I`とそこから外側のprincipal cutは未評価にする。逆三角函数の戻り値はsessionのangle modeに従うため，Series係数にもRadian/Degree/Gradianの出力scaleを反映する。一方，`Si` / `Ci` / Fresnel C/Sの定義核は常にRadianであり，sessionの角度modeには依存しない。Lambert WはDLMF 4.13.4の導函数多項式`p_n(W)`を整数係数漸化式で生成し，正則中心の局所Taylor係数を既存TPSA/Puiseuxへ合成する。principal `W_0`は正則中心では`(-Infinity,-1/E]`をcutとして避け，明示整数branch `k!=0`は`(-Infinity,0]`を避ける。分岐点`z=-1/E`ではDLMF 4.13.9_1–4.13.9_2の`s=sqrt[E z+1]`展開を用い，exact中心`z=-1/E`における`W_0`と`W_-1`をsquare-root Puiseux級数として扱う。`d_n`は偶数次の有理係数と奇数次の`有理数*sqrt[2]`へ分離してexact漸化式で生成し，一般algebraic simplifierを係数生成ループへ持ち込まない。principal側は正の先頭係数からprincipal `sqrt`を証明できる方向だけを採用し，`W_-1`は同じ局所変数の符号を反転して反対sheetを選ぶ。他のbranch，cutへ向かう負の先頭方向，`sqrt[x^2]`を生む高重複接触は推測せず未評価に保つ。`gamma` / `lgamma`も同じ局所基盤へ接続している。DLMF 5.7の`log Gamma(1+z)`および半整数基底の係数を既存TPSAへ合成し，`Gamma(z+1)=z Gamma(z)`の対数微分を使ってexact整数・半整数中心へ移送する。係数生成中に`gamma`そのものを高階微分せず，`lgamma`も一度Gamma級数を作ってから`log`するのではなく局所対数増分を直接構成するため，高次での巨大な相殺式を避ける。`gamma`は正則な対応中心から複素方向やPuiseux引数へ合成できる。一方，現行`lgamma[x]`は実軸上の`log[abs[gamma[x]]]`であり複素`LogGamma`ではないため，Seriesでも実係数の局所方向だけを扱う。非正整数のGamma極，現基底でexact係数を構成できない`1/3`等の一般有理中心，複素方向の`lgamma`は未評価に保つ。`digamma` / `trigamma`も同じ局所基盤を使う。DLMF 5.7.4に一致する`digamma`係数は同じlog-Gamma係数を1階微分して生成し，`trigamma`はさらに1階微分するため，係数表を重複して持たない。基底中心1・1/2からexact整数・半整数中心への移送は`psi(z+1)=psi(z)+1/z`と`psi1(z+1)=psi1(z)-1/z^2`をSeries演算で適用する。正則な対応中心では複素方向とPuiseux引数へ合成できるが，非正整数極と現基底でexact係数を構成できない一般有理中心は未評価に保つ。`polylog[s,z]`の原点局所展開はDLMF 25.12.10の定義級数から直接構成する。order `s`が展開変数に依存しなければsymbolicでもよく，係数`n^(-s)`をexact expressionとして保持したままTPSA/Puiseux引数へ合成する。原点では`Li_s(z)=z+O[z^2]`の単純零点をvaluationへ公開するため，積・商やLaurent inversionにも接続できる。正整数orderでは非零正則中心にも対応する。`D^n Li_s(z)=z^(-n) sum_k s(n,k) Li_(s-k)(z)`をsigned Stirling数で直接Taylor係数化し，非正整数orderへ到達した項は`Li_0(z)=z/(1-z)`およびEulerian多項式`Li_{-m}(z)=z A_m(z)/(1-z)^(m+1)`へexact還元する。そのため高階導函数表を持たず，負order `polylog` headも係数へ残さない。principal cut `[1,Infinity)`を避けるため，実中心は`a<1`を証明できる場合，複素中心は非実と証明できる場合だけ展開する。`0<a<1`かつ正整数orderなら定義級数から`Li_s(a)>0`であることもSeries内部の非零証明として利用し，`1/polylog[s,a+x]`や負整数冪を安全にinverseへ接続する。非整数orderの非零中心，cut上の中心，正則性を証明できないsymbolic中心は未評価に保つ。特殊函数・逆三角函数・Lambert W・Gamma/psi系・polylogの既知valuationも積・商やLaurent inversionへ接続する。

```text
series[log[x],{x,0,4}]
-> seriesData[x, 0, {0, 0, 0, 0, 0}, 0, 5, 1, {{1, 0, 0, 0, 0}}]

series[Ei[x],{x,0,4}]
-> seriesData[x, 0, {-digamma[1], 1, 1/4, 1/18, 1/96}, 0, 5, 1, {{1, 0, 0, 0, 0}}]

series[Ci[x],{x,0,6}]
-> seriesData[x, 0, {-digamma[1], 0, -1/4, 0, 1/96, 0, -1/4320}, 0, 7, 1, {{1, 0, 0, 0, 0, 0, 0}}]

series[log[2*x],{x,0,4}]
-> seriesData[x, 0, {log[2], 0, 0, 0, 0}, 0, 5, 1, {{1, 0, 0, 0, 0}}]

series[log[sqrt[x]],{x,0,4}]
-> seriesData[x, 0, {0, 0, 0, 0, 0, 0, 0, 0, 0}, 0, 9, 2, {{1/2, 0, 0, 0, 0, 0, 0, 0, 0}}]

series[log[x^2],{x,0,4}]
-> series[log[x^2], {x, 0, 4}]
```

```text
series[lgamma[1+x],{x,0,5}]
-> seriesData[x, 0, {digamma[1], Pi^2/12, -zeta[3]/3, Pi^4/360, -zeta[5]/5}, 1, 6, 1]

series[gamma[1/2+x],{x,0,2}]
-> seriesData[x, 0, {sqrt[Pi], (digamma[1]-2log[2])sqrt[Pi], (Pi^2/2+(digamma[1]-2log[2])^2)sqrt[Pi]/2}, 0, 3, 1]

series[gamma[x],{x,0,3}]
-> series[gamma[x], {x, 0, 3}]

series[lgamma[1+I*x],{x,0,3}]
-> series[lgamma[I x+1], {x, 0, 3}]

series[digamma[1+x],{x,0,5}]
-> seriesData[x, 0, {digamma[1], Pi^2/6, -zeta[3], Pi^4/90, -zeta[5], zeta[6]}, 0, 6, 1]

series[trigamma[1/2+x],{x,0,3}]
-> seriesData[x, 0, {Pi^2/2, -14zeta[3], Pi^4/2, -124zeta[5]}, 0, 4, 1]

series[polylog[2,x],{x,0,6}]
-> seriesData[x, 0, {1, 1/4, 1/9, 1/16, 1/25, 1/36}, 1, 7, 1]

series[polylog[a,sqrt[x]],{x,0,2}]
-> seriesData[x, 0, {1, 2^(-a), 3^(-a), 4^(-a)}, 1, 5, 2]

series[polylog[2,1/2+x],{x,0,3}]
-> seriesData[x, 0, {polylog[2, 1/2], -2log[1/2], 2(1+log[1/2]), 4(-1-2log[1/2])/3}, 0, 4, 1]

series[1/polylog[2,1/2+x],{x,0,2}]
-> seriesData[x, 0, {1/polylog[2, 1/2], 2log[1/2]/polylog[2, 1/2]^2, -(2(1+log[1/2])/polylog[2, 1/2]-4log[1/2]^2/polylog[2, 1/2]^2)/polylog[2, 1/2]}, 0, 3, 1]
```

```text
series[(1+x)^3,{x,0,5}]
-> seriesData[x, 0, {1, 3, 3, 1, 0, 0}, 0, 6, 1]

normal[%]
-> x^3+3x^2+3x+1
```

```text
series[1/(1-x),{x,0,4}]
-> seriesData[x, 0, {1, 1, 1, 1, 1}, 0, 5, 1]

series[1/x,{x,0,3}]
-> seriesData[x, 0, {1, 0, 0, 0, 0}, -1, 4, 1]

normal[%]
-> x^(-1)
```

```text
toNormal[{series[(1+x)^2,{x,0,3}],series[log[x],{x,0,2}]}]
-> {x^2+2x+1, log[x]}

normal[{series[(1+x)^2,{x,0,3}]}]
-> {seriesData[x, 0, {1, 2, 1, 0}, 0, 4, 1]}
```

`SolutionSet`内部でもbindingだけを再帰変換し，解集合metadataは保持する。低コスト初等函数は既存TPSAへ正規化して展開する。

```text
series[tan[x],{x,0,5}]
-> seriesData[x, 0, {1, 0, 1/3, 0, 2/15}, 1, 6, 1]

series[log2[x],{x,Infinity,3}]
-> seriesData[x, Infinity, {0, 0, 0, 0}, 0, 4, 1, {{-1/log[2], 0, 0, 0}}]
```

`toNormal`は再帰変換であり，変換後の式へ再度適用しても結果は変わらない。`SolutionSet`ではbinding右辺だけを変換し，条件・自由変数・multiplicity・domainは変更しない。将来，通常形へ落とす特殊構造が増えた場合も同じfrontendへ追加できる。


```text
series[exp[x],{x,0,5}]
-> seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120}, 0, 6, 1]

series[1/sin[x],{x,0,5}]
-> seriesData[x, 0, {1, 0, 1/6, 0, 7/360, 0, 31/15120}, -1, 6, 1]

series[log[x],{x,I,3}]
-> seriesData[x, I, {I Pi/2, -I, 1/2, I/3}, 0, 4, 1]

series[sqrt[1+x],{x,0,6}]
-> seriesData[x, 0, {1, 1/2, -1/8, 1/16, -5/128, 7/256, -21/1024}, 0, 7, 1]

series[(1+x)^(3/2),{x,0,6}]
-> seriesData[x, 0, {1, 3/2, 3/8, -1/16, 3/128, -3/256, 7/1024}, 0, 7, 1]
```

```text
series[sqrt[x],{x,0,5}]
-> seriesData[x, 0, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 1, 11, 2]

series[sqrt[x]*(1+x),{x,0,4}]
-> seriesData[x, 0, {1, 0, 1, 0, 0, 0, 0, 0}, 1, 9, 2]

series[exp[sqrt[x]],{x,0,3}]
-> seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120, 1/720}, 0, 7, 2]
```

```text
series[erf[x],{x,0,7}]
-> seriesData[x, 0, {2/sqrt[Pi], 0, -2/sqrt[Pi]/3, 0, 1/(5sqrt[Pi]), 0, -1/(3sqrt[Pi])/7}, 1, 8, 1]

series[Si[x],{x,0,7}]
-> seriesData[x, 0, {1, 0, -1/18, 0, 1/600, 0, -1/35280}, 1, 8, 1]

series[Ei[1+x],{x,0,4}]
-> seriesData[x, 0, {Ei[1], E, 0, E/6, -E/12}, 0, 5, 1]

series[Ci[1+x],{x,0,4}]
-> seriesData[x, 0, {Ci[1], cos[1 Rad], (-cos[1 Rad]-sin[1 Rad])/2, (cos[1 Rad]/2+sin[1 Rad])/3, (-cos[1 Rad]/2-5sin[1 Rad]/6)/4}, 0, 5, 1]

series[erfc[x],{x,0,5}]
-> seriesData[x, 0, {1, -2/sqrt[Pi], 0, 2/(3sqrt[Pi]), 0, -1/sqrt[Pi]/5}, 0, 6, 1]

series[fresnelc[x],{x,0,5}]
-> seriesData[x, 0, {1, 0, 0, 0, -Pi^2/40}, 1, 6, 1]

series[fresnels[x],{x,0,7}]
-> seriesData[x, 0, {Pi/6, 0, 0, 0, -Pi Pi^2/336}, 3, 8, 1]

series[li[2+x],{x,0,2}]
-> seriesData[x, 0, {li[2], 1/log[2], -1/(2log[2]^2)/2}, 0, 3, 1]

series[li[a+x],{x,0,2},a>1]
-> seriesData[x, 0, {li[a], 1/log[a], -1/(a log[a]^2)/2}, 0, 3, 1]

series[asin[x],{x,0,7}]
-> seriesData[x, 0, {1, 0, 1/6, 0, 3/40, 0, 5/112}, 1, 8, 1]

series[atan[1+I+x],{x,0,3}]
-> seriesData[x, 0, {atan[1+I], 1/5-2I/5, -1/25+7I/25, -1/375-68I/375}, 0, 4, 1]

series[asin[a+x],{x,0,2},-1<a<1]
-> seriesData[x, 0, {asin[a], (1-a^2)^(-1/2), a*(1-a^2)^(-1/2)/(2(1-a^2))}, 0, 3, 1]

series[lambertw[x],{x,0,7}]
-> seriesData[x, 0, {1, -1, 3/2, -8/3, 125/24, -54/5, 16807/720}, 1, 8, 1]

series[lambertw[E+x],{x,0,4}]
-> seriesData[x, 0, {1, exp[-1]/2, -3*exp[-2]/16, 19*exp[-3]/192, -185*exp[-4]/3072}, 0, 5, 1]

series[1/lambertw[x],{x,0,5}]
-> seriesData[x, 0, {1, 1, -1/2, 2/3, -9/8, 32/15, -625/144}, -1, 6, 1]

series[lambertw[-1/E+x],{x,0,3}]
-> seriesData[x, 0, {-1, sqrt[2]sqrt[E], -2*E/3, 11*E sqrt[2]sqrt[E]/36, -43*exp[2]/135, 769*exp[2]sqrt[2]sqrt[E]/4320, -1768*E exp[2]/8505}, 0, 7, 2]

series[lambertw[-1,-1/E+x],{x,0,3}]
-> seriesData[x, 0, {-1, -sqrt[2]sqrt[E], -2*E/3, -11*E sqrt[2]sqrt[E]/36, -43*exp[2]/135, -769*exp[2]sqrt[2]sqrt[E]/4320, -1768*E exp[2]/8505}, 0, 7, 2]
```

通常REPLの複数行／自動layoutでは，例えば`series[sqrt[x]*(1+x),{x,0,4}]`を`x^(1/2) + x^(3/2) + O[x^(9/2)]`と表示する。`single`およびmachine-facing出力では`seriesData[...]`を保持する。

`series[log[x],{x,-1,n}]`のようなprincipal branchの分岐切断上を中心とする要求や，`series[sqrt[x^2],{x,0,n}]`のように単一のprincipal局所branchへ安全に落とせない要求は未評価のまま保持する。

`SeriesData`は展開変数自身について`D` / `integrate`で直接係数演算できる。通常係数だけでなく各`log(x-center)^k`層も積の微分則・部分積分漸化式で処理する。したがってTaylor/Laurent/Puiseux指数格子を保ったまま，`D[log[x]^k]`や`integrate[x^-1 log[x]^k,x]`もSeriesData内で閉じる。

```text
D[series[exp[x],{x,0,5}],x]
-> seriesData[x, 0, {1, 1, 1/2, 1/6, 1/24}, 0, 5, 1]

integrate[series[sqrt[x],{x,0,4}],x]
-> seriesData[x, 0, {2/3, 0, 0, 0, 0, 0, 0, 0}, 3, 11, 2]
```

`x^(-1)`項の積分はlog係数層へ移り，例えば`integrate[series[1/x,{x,0,4}],x]`は`log[x]`を表すSeriesDataを返す。`x^-1 log[x]^k`も`log[x]^(k+1)/(k+1)`へexactに移る。

`seriesData[variable,center,coefficients,minExponent,orderNumerator,exponentDenominator]`が従来の内部表現である。対数項を持つ場合だけ7引数目`logarithmicCoefficientLayers`を追加し，layer `k`は`log(x-center)^(k+1)`へ掛かる同一指数格子の係数列として解釈する。対数層が空なら6引数canonicalを維持する。通常利用では`series` / `normal`を使う。各係数の指数は`(minExponent+i)/exponentDenominator`，剰余次数は`orderNumerator/exponentDenominator`である。

---

# 23. 記号微分と数値微分

## 23.1 `D`

```text
D[expr,x]
D[expr,{x,n}]
D[expr,x,y,...]
```

`{x,n}`は非負整数`n`階微分，複数specは左から順に適用する。

```text
D[sin[x],{x,4}]
-> sin[x]

D[x^2 y^3,x,y]
-> 6 x y^2
```

HoldAllなので既存の変数値に置換せず式を微分する。

```text
D[x^3 + 2x,x]
-> 2 + 3 x^2

D[exp[x^2],x]
-> 2 x exp[x^2]
```

既定Radでは:

```text
D[sin[x],x]
-> cos[x]
```

`abs/sign/re/im/conj/arg`など一般複素変数で通常のholomorphic derivativeを持たないものは，偽の微分を返さず未評価`D[...]`を保持する。

0で可除特異点を連続延長して定義している安定化函数は，微分後もその点を失わない。例えば

```text
D[sinc[x],x]
-> cases[(x cos[x]-sin[x])/x^2 if x != 0; 0 if x == 0]

D[cosc[x],x]
-> cases[(-1+x sin[x]+cos[x])/x^2 if x != 0; 1/2 if x == 0]
```

これは第一級scalar値である`cases[...]`であり，評価制御`if[...]`でも`SolutionSet`でもない。`if[...]`は評価時に1 分岐を選択し，`cases[...]`は数学的分岐条件をsymbolicに保持する。`SolutionSet`は方程式の解集合専用である。

積分でも同じ微分知識を再利用する。

```text
D[integrate[f[x],x],x]
-> f[x]

D[integrate[t^2,{t,0,x}],x]
-> x^2
```

可変上下端や積分内部に微分変数が現れる一般形ではLeibniz ruleを形式的に構築する。分母が微分変数に依存しない商は，一般quotient ruleへ膨張させず`f'/c`を直接使う。

## 23.2 `diff`

```text
diff[expr,x,at]
diff[expr,x,at,digits]
```

有限差分の別公式ではなく，まず`D`でexactな導函数Exprを作り，それを指定点で`CertifiedEvaluator`へ渡す。有限精度入力では`CertifiedEnclosure`と`InformationEnclosure`を別々に評価し，結果の精度は入力が持つ情報量を超えない。被微分式の中にholdされた`N[...]`があっても同じ規則で扱う。入力の情報幅だけでは特異点や分岐切断の側を確定できない場合は，値を推測せずEvaluationErrorとする。
第4引数`digits`は小数部桁数であり，`N[...,p]`の有効桁数とは異なる。入力自身の有限精度がそれより低い場合は，`InformationEnclosure`が保証できる桁数を上限とする。

```text
diff[x^2,x,3]
-> 6.0
```

---

# 24. 記号積分とcertified数値積分

## 24.1 `integrate` — exact/symbolic積分

```text
integrate[expr,x]
integrate[expr,{x,a,b}]
integrate[expr,x,assumptions]
integrate[expr,{x,a,b},assumptions]
```

第1形式は不定積分，第2形式はexact/symbolic定積分。積分変数はbinderとしてholdされ，同名のglobal定義に置換されない。

不定積分は「加法定数を法として選んだ原始函数の代表元」を返す。したがって`+ C`は表示しない。

```text
integrate[x^2,x]
-> x^3 / 3

integrate[(2x+3)^5,x]
-> (3 + 2 x)^6 / 12

integrate[1/(2x+3),x]
-> log[3 + 2 x] / 2
```

現在の主なexact規則:

- 定数，`x`，任意の有限多項式
- affine baseの有理冪。指数`-1`はLogへ送る
- Rational係数の有理函数。1次/2次因子はexact partial fractionへ分解し，重複既約2次因子`(a x^2+b x+c)^k`も平方完成に基づくexact recurrenceで処理する。次数3以上を含む分母はQ[x]上のsquare-free decompositionとHermite reductionで重複度を落とし，残るsquare-free因子をcertified Complex `root[...,k,Complex]`とexact residue `P(r)/Q'(r)`によるalgebraic-log和へ分解する。elementary reverse-chainを先に試すため，`x/(1+x^4)`等の簡潔な`atan/asin`表現を優先する。specialized rational pathは現在degree 12まで。`x^m/(1+x^n)`型は必要に応じ`hypergeometric2F1` primitiveへも接続する
- 正のRational scaleを証明できる二次逆平方根型の`asin/asinh` primitive，およびexact Rational係数二次式`q(x)`の`sqrt[q(x)]` primitive
- `sin^m/cos^n`の有限Fourier reduction。積分器は正整数総次数256までを明示的に展開可能
- `sin[u]^(-n)` / `cos[u]^(-n)` (`1<=n<=256`) を `csc/sec` の標準漸化式で積分
- `tan/cot/sec/csc`の正整数冪 (`2<=n<=256`) を標準reduction formulaで積分
- 和・差・符号反転，積分変数に依存しない係数の線形性
- `exp/sin/cos/tan/cot/sec/csc`の安全な標準原始函数
- `sinh/cosh/tanh/coth/sech/csch`の安全な標準原始函数
- `log/log1p/expm1/sqrt/cbrt`。`log[x]/x`や`1/(x log[x])`は対数微分Knowledgeから認識
- `asin/acos/atan/asinh/acosh/atanh`
- `erf/erfc`
- `fresnelc/fresnels`。exact Rational係数の`sin/cos[a x^2+b x+c]`を平方完成して標準Fresnel積分へ還元。`Pi*x^2/2`の定義核も直接認識
- Gaussian family。`a`がexact positive Rationalの`exp[-a x^2]`は一般1F1より`erf`をpreferred canonical primitiveとして選び，improper endpointでも`erf`の±Infinity endpoint極限からexact Gaussian積分へ閉じる
- `hypergeometric1F1`。上記Gaussian以外の正整数`n>=2`の`exp[c x^n]`を原点で正則な1F1原始函数へ還元
- exactな逆chain rule。`f'(x) f(x)^p`はD後の偶然の式形に依存せず構造的にも認識
- 多項式×`exp/sin/cos/sinh/cosh`に対する有限回のintegration by parts
- `exp[a x+b] sin/cos[c x+d]`型を連立一次式としてexact積分
- 主値 `sqrt[x]`を含む有理的な形への`t=sqrt[x]`局所置換，および`sqrt[q(sqrt[x])]`の二次根号class
- 共通引数を持つ`R(sin(theta),cos(theta))`に対するbounded Weierstrass置換`t=tan(theta/2)`。変換後は既存exact有理積分器へ渡す
- `log[1+beta*x^n]/x`を`polylog[2,-beta*x^n]`へ還元するdilogarithm Knowledge
- boundedな積×短い和の分配。ただし全体がexact chain-ruleで一発に閉じる場合はchain-ruleを優先し，能力退行を防ぐ
- `x^n log[x]` (`n`が非負整数)

例:

```text
integrate[exp[2x+1],x]
-> exp[1 + 2 x] / 2

integrate[sin[2x],x]
-> -cos[2 x] / 2

integrate[tan[x],x]
-> -log[cos[x]]

integrate[sech[x],x]
-> atan[sinh[x]]

integrate[1/(1-x^2),x]
-> atanh[x]

integrate[1/(x^2-1),x]
-> -atanh[x]

integrate[(x+1)/(x+2),x]
-> x - log[2 + x]

integrate[(x+1)/(x^2+4),x]
-> log[4 + x^2] / 2 + atan[x/2] / 2

integrate[1/sqrt[4-x^2],x]
-> asin[x/2]

integrate[1/sqrt[x^2+4],x]
-> asinh[x/2]

integrate[sin[x]^2,x]
-> (x - sin[2 x] / 2) / 2

integrate[log[x]/x,x]
-> log[x]^2/2

integrate[sec[x]^3,x]
-> sec[x]tan[x]/2+log[sec[x]+tan[x]]/2

integrate[exp[x^6],x]
-> x hypergeometric1F1[1/6, 7/6, x^6]

integrate[sin[2x]^(-2),x]
-> -cot[2x]/2

integrate[cos[4x^2],x]
-> fresnelc[x sqrt[8/Pi]]/sqrt[8/Pi]

integrate[sin[2x^2]^4,x]
-> 3x/8+fresnelc[x sqrt[16/Pi]]/(8sqrt[16/Pi])-fresnelc[x sqrt[8/Pi]]/sqrt[8/Pi]/2

integrate[asin[x],x]
-> x asin[x] + sqrt[1 - x^2]

integrate[erf[x],x]
-> x erf[x] + exp[-x^2] / sqrt[Pi]

integrate[x exp[x],x]
-> x exp[x] - exp[x]

integrate[x^2 log[x],x]
-> x^3 log[x] / 3 - x^3 / 9

integrate[x^2+sin[x],x]
-> x^3 / 3 - cos[x]

integrate[E^x cos[x],x]
-> (cos[x] + sin[x]) / 2 * exp[x]

integrate[1/(x^3+1),{x,0,1}]
-> log[2] / 3 + Pi sqrt[3] / 9
```

今回追加した置換・接続の代表例:

```text
integrate[sqrt[4-x^2],x]
-> x sqrt[4-x^2]/2+2asin[x/2]

integrate[2*x*(1+x^2)^5,x]
-> (1+x^2)^6/6

integrate[x/(1+x^4),x]
-> atan[x^2]/2

integrate[cos[Pi*x^2/2],x]
-> fresnelc[x]

integrate[log[1+x^2]/x,x]
-> -polylog[2,-x^2]/2

integrate[1/(1+sin[x]),x]
-> -2/(1+tan[x/2])
```

角度単位も既存のAngleSemantics/Dと共有する。

```text
integrate[sin[x Deg],x]
-> -(180 / Pi cos[x Deg])
```

逆chain rule等で構造から候補原始函数を発見した場合は，既存の`D`をproof engineとして使い，候補の微分と元 integrand のexactな比例関係を証明してから採用する。単に数値点で一致した候補は採用しない。

branch/definednessを壊す**global simplification**は行わない。一方，原始函数は大域恒等式と同じ基準である必要はない。共通の解析領域上で正しい局所原始函数は，integrate専用の規則として採用できる。

例えば主値 square rootについて

```text
integrate[1/sqrt[x^2-1],x]
-> log[sqrt[x^2-1] + x]
```

を返すが，Simplifierへ

```text
sqrt[x^2-1] == sqrt[x-1] sqrt[x+1]
```

という危険な大域規則は追加しない。積分公式と代数的恒等式の正当性を分離している。

また和の線形性は部分評価する。

```text
integrate[x^2 + gamma[x],x]
-> x^3 / 3 + integrate[gamma[x],x]
WARN: integrate partially evaluated the expression; remaining subintegral(s) are outside the current symbolic rule set
```

解けない部分だけを保持し，既に求まった項まで巻き戻さない。

### 未評価理由の診断

未評価を単一のWARNへ潰さず，現在は次のdiagnostic codeを区別する。

- `integrate::unsupported` — 現在のsymbolic rule setに解法がない。**閉形式が存在しないことを意味しない**。
- `integrate::partial` — 一部は積分済みだが，残るsubintegralが現在のrule外。
- `integrate::conditionsRequired` — 定義域 / branch仮定不足で安全なprimitiveを選択できない。
- `integrate::noKnownClosedForm` — mmCalの現在の標準函数語彙で有限閉形式がない代表familyとして明示的に認識したもの。

例えば，

```text
integrate[gamma[x],x]
-> integrate[gamma[x],x]
WARN: mmCal has no implemented symbolic integration rule for this expression; this does not imply that no closed form exists

integrate[abs[x],x]
-> integrate[abs[x],x]
WARN: integrate needs additional domain or branch assumptions before it can choose a safe symbolic antiderivative

integrate[x^x,x]
-> integrate[x^x,x]
WARN: integrate recognized a family with no known finite closed form in mmCal's supported standard-function vocabulary; the integral remains unevaluated
```

`noKnownClosedForm`も「あらゆる数学的表現で不可能」という判定ではない。級数，新しい特殊函数，より広い函数classを許せば表現できる場合がある。mmCalが主張する範囲を現在サポートする有限標準函数語彙に限定する。

### nested square-root substitution

現在は次のclassも扱う。

```text
integrate[sqrt[x + sqrt[x]],x]
-> 2 (x + sqrt[x]) sqrt[x + sqrt[x]] / 3
   - ((1 + 2 sqrt[x]) sqrt[x + sqrt[x]] / 4
      - log[1 + 2 sqrt[x + sqrt[x]] + 2 sqrt[x]] / 8)
```

これは単発公式ではなく，`t=sqrt[x]`により`2 t sqrt[q(t)]`へ落ちる，`q`がexact Rational係数2次式で正leadingのclassを処理する。主値分岐上の局所置換則であり，一般のradical substitution探索ではない。

### exact/symbolic定積分

原始函数を安全に得られた場合，上下端へexact substitutionして差を取る。定義域条件がある函数は，可能ならCertifiedEvaluatorで**区間全体**をpreflightし，途中の極や分岐上の不成立を見逃さない。

```text
integrate[x^2,{x,0,1}]
-> 1/3

integrate[sin[x],{x,0,Pi}]
-> 2

integrate[1/x,{x,1,2}]
-> log[2]

integrate[log[x],{x,1,Pi}]
-> 1 + Pi log[Pi] - Pi
```

特異点を跨ぐ場合はendpoint代入だけで値を作らない。

```text
integrate[1/x,{x,-1,1}]
-> WARN + unevaluated integrate[...]

integrate[tan[x],{x,0,2}]
-> WARN + unevaluated integrate[...]
```

Cauchy 主値を自動的に意味することもない。

### assumptions と improper integral

第3引数にassumptionを渡せる。既存の`KnowledgeContext`へ統合され，`abs`や`sqrt[x^2]`のbranch-sensitive簡約に利用する。

```text
integrate[abs[x],x,x>=0]
-> x ^ 2 / 2

integrate[abs[x],x,x<=0]
-> -x ^ 2 / 2

integrate[sqrt[x^2],x,x>=0]
-> x ^ 2 / 2
```

endpointがInfinityまたは通常代入で定義されない場合は，原始函数の対応する片側/無限遠極限を使ってimproper integralを評価する。内部特異点がないことを証明できるclassだけを受理し，Cauchy 主値は推測しない。

```text
integrate[exp[-x],{x,0,Infinity}]
-> 1

integrate[1/x^2,{x,1,Infinity}]
-> 1

integrate[1/(1+x^2),{x,-Infinity,Infinity}]
-> Pi

integrate[1/sqrt[x],{x,0,1}]
-> 2

integrate[log[x],{x,0,1}]
-> -1

integrate[exp[-a*x],{x,0,Infinity},a>0]
-> 1/a

integrate[x^(s-1)*exp[-x],{x,0,Infinity},s>0]
-> gamma[s]

integrate[x^(a-1)*(1-x)^(b-1),{x,0,1},{a>0,b>0}]
-> beta[a,b]

integrate[1/(1+x^4),{x,0,Infinity}]
-> Pi/(2sqrt[2])

integrate[log[x]^2,{x,0,1}]
-> 2

integrate[sin[x]/x,{x,0,Infinity}]
-> Pi/2

integrate[cos[Pi*x^2/2],{x,0,Infinity}]
-> 1/2
```

これらのパラメータ付き familyは仮定から収束条件・実分岐を証明できた場合だけGamma/Beta/Mellin等へ還元する。finite intervalのsymbolic 極も区間外と証明できればexact endpoint ratioを使うが，区間内部の可能性が残れば評価しない。

途中に極を含む可能性を排除できない場合は未評価保持する。

```text
integrate[1/(x-2),{x,1,Infinity}]
-> WARN + unevaluated integrate[...]
```

## 24.2 `limit` — exact/symbolic limit

```text
limit[expr,x,a]
limit[expr,x,a,-1]   // 左極限
limit[expr,x,a,1]    // 右極限
limit[expr,{x,a,direction}]
```

第4引数またはbrace形式の第3要素は方向を表し，`-1`が左，`1`が右。省略時は二側極限。`limit[expr,{x,a,direction}]`は4引数形と等価である。方向は単なる表示指定ではなく一時的なassumption `x<a` / `x>a` としてKnowledgeContextへ渡される。

```text
limit[sin[x]/x,x,0]
-> 1

limit[(1-cos[x])/x^2,x,0]
-> 1/2

limit[1/x,x,0,1]
-> Infinity

limit[1/x,x,0,-1]
-> -Infinity

limit[abs[x]/x,x,0,1]
-> 1

limit[abs[x]/x,x,0,-1]
-> -1

limit[atan[x],x,Infinity]
-> Pi / 2

limit[exp[-x],x,Infinity]
-> 0

limit[sin[1/x],x,0]
-> Indeterminate

limit[sin[1/x],x,0,1]
-> Indeterminate

limit[x sin[1/x],x,0]
-> 0

limit[x*Ei[x],x,0,1]
-> 0

limit[sqrt[x]*log[x],x,0,1]
-> 0

limit[Ei[x]-log[x],x,0,1]
-> -digamma[1]

limit[Ci[x]-log[x],x,0,1]
-> -digamma[1]

limit[log[2*x]-log[x],x,0,1]
-> log[2]

limit[Ei[x],x,0]
-> -Infinity

limit[Ei[x],x,Infinity]
-> Infinity

limit[Ei[x],x,-Infinity]
-> 0

limit[Ci[x],x,0]
-> -Infinity

limit[Ci[x],x,Infinity]
-> 0

limit[Ci[x],x,-Infinity]
-> I Pi

limit[li[x],x,0]
-> 0

limit[li[x],x,1]
-> -Infinity

limit[li[x],{x,1,1}]
-> -Infinity

limit[li[x],x,Infinity]
-> Infinity

limit[li[x],x,-Infinity]
-> ComplexInfinity
```

`Ei` / `Ci` / `li`は主値分岐として扱う。実軸上では`Ei[x]`は`x -> 0`で`-Infinity`，`x -> -Infinity`で0へ収束する。`Ci[x]`は`x -> 0`で`-Infinity`，正の無限遠で0へ収束する一方，負の実軸は分岐切断上にあるため`x -> -Infinity`では分岐オフセットを保持して`I Pi`へ収束する。`li[x]`は原点で0へ収束する。負の実軸方向で`x -> -Infinity`とした場合は現行の方向表現では実数の`-Infinity`へ潰さず，方向未定の複素無限大`ComplexInfinity`を返す。

0/0型では既存`D`を使った反復l'Hopitalを安全弁付きで利用する。Rational functionの局所zero/極次数や無限遠次数比較はexactに処理する。これらの個別規則で決まらない有限点では，対応可能な場合に限って局所`SeriesData`を補助backendとして使い，先頭の非零項から有限定数または0を証明する。`Ei[x]-log[x]`等の`Infinity-Infinity`相殺もこの経路で扱う。負冪や定数次数の`log^k`が残る場合は発散方向を推測せず未評価へ戻すため，Series backendは既存limit kernelの置換ではない。`sin` / `cos` / `tan`について，実引数が一側または無限遠で`+/-Infinity`へ走ることをexactに証明できた場合は，周期振動により単一の極限値が存在しないことも証明済みとして`Indeterminate`を返す。有限点の二側極限では片側の無限振動だけでも不存在を確定できる。さらに，実Rational functionを引数に取る`sin` / `cos`は実軸上で絶対値1以下であることを使い，`x sin[1/x] -> 0`等の2因子積をsqueeze theoremでexactに閉じる。未解決二側極限を主値等へ潰すことはしない。

```text
limit[1/x,x,0]
-> WARN + limit[1 / x, x, 0]
```

## 24.3 `nintegrate` — certified数値積分

```text
nintegrate[expr,{x,a,b}]
nintegrate[expr,{x,a,b},digits]
```

```text
nintegrate[x^2,{x,0,1},12]
-> 0.333333333333
```

内部では区間を`x=a+(b-a)t`へ正規化し，exact Rational格子上のNewton–Cotesと高階導函数のcertified boundから求積誤差を包含する。

単なる`double` Simpsonの「近そうな値」ではなく，最終丸めが一意になった場合だけ返す。有限精度の積分区間・被積分函数については`CertifiedEnclosure`と`InformationEnclosure`を並行して積分し，入力以上の精度を作らない。holdされた`N[...]`も同じ情報量制約を保つ。
第3引数`digits`は小数部桁数を表す。したがって値の整数部が大きくても，入力情報が十分なら指定した小数部桁数まで確定する。有限精度入力が先に精度を制限する場合は，その情報量を越えて表示しない。

高階導函数を構築する前に元の被積分函数を区間全体でpreflightするため，明白な特異点を早期に拒否する。特異点を跨いで偶然相殺する処理は行わない。

---

# 25. Solver

```text
solve[equation,x]
solve[equation,x,domainOrConstraint]
solve[equation,domain]
solve[{equations...},{variables...}]
```

既定ambient 定義域は等式系でComplex。
ordered inequalityはRealまたはそのsubdomainで扱う。

`solve[equation,domain]`は`domain`が`Integer` / `Rational` / `Real` / `Complex`で，方程式中の未知user symbolを**ちょうど1個**に確定できる場合だけ変数を推定する短縮形である。未知数が0個または複数なら推測せずTypeErrorとする。明示変数形でも`Pi` / `Real`等の予約・builtin symbolをsolve変数として受理しない。

```text
solve[x^2 == 1,x]
-> {x == 1, x == -1}

solve[x^2 + 1 == 0,x,Real]
-> {}

solve[x^2 + 1 == 0,x,Complex]
-> {x == I, x == -I}

solve[x^2 < 4,x]
-> {x in Real if x > -2 && x < 2}

solve[{2x+3y==5,x-2y==9},{x,y}]
-> {{x==37/7, y==-13/7}}
```

`SolutionSet`はEmpty / Finite / Universal / Conditional / Unresolvedを区別する。
対応外の式を「解なし」と誤認しない。

### 多変数多項式・Gröbner basis

exact Rational係数の多変数多項式は`Q[x1,...,xn]`の専用ring representationへ変換する。公開接口は次である。

```text
groebnerBasis[polys,{x,y,...}]
groebnerBasis[polys,{x,y,...},Lex]
groebnerBasis[polys,{x,y,...},GrLex]
groebnerBasis[polys,{x,y,...},GrevLex]
polynomialReduce[f,G,{x,y,...}]
polynomialReduce[f,G,{x,y,...},order]
```

既定term orderは`GrevLex`。実装はmultivariate division，normal form，S-polynomial，Buchberger product/chain criteria，sugar pair selection，interreduction，reduced Gröbner basisをexact Rational arithmeticだけで行う。変数数や次数を固定しないが，critical pair，reduction回数，basis size，term数はboundedであり，資源爆発時は数学的結果を捏造しない。

```text
groebnerBasis[{x y-1,y^2-x},{x,y},Lex]
-> {x-y^2, y^3-1}

polynomialReduce[x^2+y^2,{x-y,y^2-1},{x,y},Lex]
-> {{x+y, 2}, 2}

polynomialReduce[x y-1,groebnerBasis[{x y-1,y^2-x},{x,y},Lex],{x,y},Lex]
-> {{y, 1}, 0}
```

`polynomialReduce[...]`内，または別の`groebnerBasis[...]`内へdirect `groebnerBasis[...]`結果を合成できる。この解決はheld polynomial変数を現在のsession bindingで一般評価せず，polynomial-ideal境界内だけで行う。

`solve[{...},{...}]`は非線形polynomial equalitiesを検出するとLex Gröbner eliminationを試す。矛盾idealは`{}`，zero-dimensional shape-position basisはunivariate exact Solve + back substitutionで列挙し，元方程式をexact polynomial remainderで検証する。shape-positionでない場合でも，univariate eliminantをexactに解け，その各rootで元systemをspecializeした低次元問題をすべて完全に解けるときはbounded recursionで全枝を列挙する。

positive-dimensional systemは一般多様体を推測してパラメータ化しない。ただし`f*g==0`をexact factor分解できる場合は`f==0` / `g==0`へ完全分岐し，`solve[{x*y==0},{x,y}] -> {x == 0 where y in Complex, y == 0 where x in Complex}`のように既存SolutionBranchの自由変数表現を使う。また複数方程式系でも，ある変数が**定数係数の一次式**としてexactに消去できる場合はそのbindingを残りのsystemへ代入し，低次元systemを再帰的に完全solveする。これにより`solve[{z-x-y==0,x^2+y^2==1},{x,y,z},Real]`は`x in [-1,1]`を保った2枝へ閉じる。残りが1方程式になり，あるsolver変数について次数1または2のsymbolic-coefficient polynomialとして完全に解ける場合は，他変数を自由パラメータとしてliftする。Complexでは`solve[{x*y==1},{x,y}] -> {y == 1/x where x in Complex if x != 0}`となる。Realでは係数をReal 自由パラメータとして扱い，一次式の分母非零条件と，二次式の主値 `sqrt`に必要な判別式／radicandの非負条件を分岐へ保持する。1自由パラメータの二次条件は一変数polynomial inequality solverでexact intervalへ正規化するため，`solve[{x^2+y^2==1},{x,y},Real] -> {y == sqrt[1-x^2] where x in Real if x >= -1 && x <= 1, y == -sqrt[1-x^2] where x in Real if x >= -1 && x <= 1}`となる。複数自由パラメータでは球面等の半代数Predicateをそのまま保持する。非定数係数の一次消去，高次function-field algebraic equation，一般多様体のrational parameterizationは引き続き`UnresolvedSolutionSet[...]`を保つ。

Solverは`HoldAll`の入力を一般Evaluatorへ流さず，専用の**solve-safe normalization**を通してから分類する。
この層はbuiltin aliasをcanonical headへ揃え，証明付きSimplifier rewriteだけを適用する。
そのため評価副作用を起こさずに`E^x`と`exp[x]`，`ln`と`log`，`log2` / `log10`等の表現差を解法能力差へ漏らさない。

分母zero，Logのdefinedness，rational-functionのhole/極等を可能な範囲でglobal conditionとして保持する。

Gammaの極集合のように現Predicateで完全表現できない条件は，不完全な条件を捏造せずunresolvedのまま扱う。

### 実軸global inverse solve

MathRegistryは主値 inverseだけでなく，実軸上のglobal injectivity・単調性・実値域をmetadataとして持つ。`solve`は**Realまたはそのsubdomain**で，globalに一対一であることを証明できる函数だけを安全に反転する。

```text
solve[exp[x]==2,x,Real]
-> {x == log[2]}

solve[log[x]==2,x,Real]
-> {x == exp[2]}

solve[sinh[3x]==2,x,Real]
-> {x == asinh[2] / 3}

solve[tanh[x]==1/2,x,Real]
-> {x == atanh[1/2]}

solve[tanh[x]==2,x,Real]
-> {}

solve[exp[x]==a,x,Real]
-> {x == log[a] if a in Real && a > 0}

solve[E^x==8,x,Real]
-> {x == log[8]}

solve[ln[x]==2,x,Real]
-> {x == exp[2]}

solve[log2[x]==3,x,Real]
-> {x == 8}
```

### Real-定義域 nonexistence / uniqueness proof

Real equalityには，具体的な逆函数解法の後段にexact proof layerを持つ。この層は「根を数値探索して見つからなかった」ことを非存在証明には使わない。現在は次を低コスト順に試す。

- `ValueFacts`がresidual `f(x)=lhs-rhs`を全実軸でstrict positive / strict negative / nonzeroと証明すれば`{}`。
- MathRegistryのreal 値域からtargetが像の外にあることを証明できれば`{}`。global injectivityと既知exact anchor `f(c)`が一致すれば，inverse builtinがなくてもargumentを`c`へ帰着する。
- `f`と`f'`がRealかつdefinedであるconnected 定義域 piece上で，`f'`がstrict signを持ち，exact anchorで`f(c)=0`を証明できれば，そのpieceのrootを唯一解として返す。また`f'>=0` / `f'<=0`しか証明できなくても，`f'=0`の完全解集合が有限またはIntegerパラメータ族として高々可算であり，実intervalを含まないことまで証明できればstrict 単調性へ昇格する。
- `f''`が全実軸でstrict signを持ち，`f'=0`の唯一のexact 臨界点を構成できれば，そこでのglobal minimum / maximumの符号から`{}`または接する唯一解を証明する。

一変数Real解析は，algebraicに証明できる定義域条件をconnected intervalsへ分解し，各intervalをexact 臨界点で再分割する。各pieceではstrict 単調性，左右endpoint limit，およびそれらから従う値域をexactに構成する。例えば`log[x]`は`(0,Infinity)`上で増加し値域は`(-Infinity,Infinity)`，`sqrt[x]`は`[0,Infinity)`上で増加し値域も`[0,Infinity)`，`atanh[x]`は`(-1,1)`上で増加し値域は`(-Infinity,Infinity)`と証明する。`1/(x^2-1)`のようなdisconnected 定義域は極を跨いで結合しない。

この定義域分解は現在，有限個のalgebraic境界とexactに解ける臨界点へ計算量制限付きで限定する。periodic 極の無限集合や，複数のcomplex-valued subexpressionが相殺して実数へ戻り得る式などで完全なreal-valued 定義域を証明できない場合は解析をUnknownに戻す。また，一意性を証明できてもexact rootを既存表現で構成できなければ`UnresolvedSolutionSet`を維持する。例えば`erf[x]==1/2`は実軸上で一意だが，現行mmCalにはinverse-erfまたは一般transcendental Root表現がないため未解決である。

`erf`のstrict 単調性はDLMF 7.10.1の`erf'(x)=2 exp(-x^2)/sqrt(Pi)>0`をMathRegistryのreal behavior knowledgeへ反映している。proof layerの構成は，Wolfram Languageの`Reduce` / `FunctionRange` / exact global optimizationが「解集合・値域・極値を数学的条件として扱う」設計を参考にするが，mmCalでは現在実装済みのexact certificateだけを採用する。

```text
solve[exp[x]==x,x,Real] -> {}
solve[exp[x]+x^2+1==0,x,Real] -> {}
solve[x+exp[x]-1==0,x,Real] -> {x == 0}
solve[exp[x]==x+1,x,Real] -> {x == 0}
solve[exp[x]-x+5==0,x,Real] -> {}
solve[log[x]-x-1==0,x,Real] -> {}
solve[log[x]-x+1==0,x,Real] -> {x == 1}
solve[1/x+x==0,x,Real] -> {}
solve[1/x-x==0,x,Real] -> {x == -1, x == 1}
solve[sin[x]==x,x,Real] -> {x == 0}
solve[erf[x]==0,x,Real] -> {x == 0}
solve[erf[x^2-1]==0,x,Real] -> {x == 1, x == -1}
solve[erf[x]==1,x,Real] -> {}
solve[erfc[x]==1,x,Real] -> {x == 0}
solve[erfc[x]==0,x,Real] -> {}
solve[erfc[x]==2,x,Real] -> {}
solve[erf[x]==1/2,x,Real] -> UnresolvedSolutionSet[x]
solve[exp[x]==x+2,x,Real] -> {x == -lambertw[-exp[-2]]-2, x == -lambertw[-1, -exp[-2]]-2}
solve[2^x==x,x,Real] -> {}
solve[(4/3)^x==x,x,Real] -> {x == -lambertw[-log[4/3]]/log[4/3], x == -lambertw[-1, -log[4/3]]/log[4/3]}
solve[x^x==1,x,Real] -> {x == 1}
solve[x^x==2,x,Real] -> {x == exp[lambertw[log[2]]]}
```

`exp[p x+q]==c x+d`および，positive constant baseを持つ`a^(m x+n)==c x+d`は，係数が実数で非零性をexactに証明できる場合にLambert W normal formへ変換する。実分岐は変換後argument `z` と`-1/E`の大小をcertifiedに比較して0/1/2本を完全分類し，分岐点では重複した`W_0/W_-1`を1根へ併合する。

主値 `x^x==r`は負実軸で一般に複素値となり，`0<r<1`では負の偶整数解が混在し得る。このため現段階では，負実根を完全に排除できる`r>1`，およびexactな`r=1,0,-1`，`r<-1`だけを分類する。例えば`solve[x^x==1/4,x,Real]`は不完全な正実分岐だけを返さず`UnresolvedSolutionSet[x]`を維持する。

Complex 定義域へこの証明を流用しない。たとえば`solve[exp[x]-x+5==0,x,Complex]`はReal非存在証明だけを理由に`{}`へしない。

### 絶対値方程式・不等式

Real relationでは`abs[u]`の値域`[0,Infinity)`を保ったままexactに変形する。ordered inequalityはReal semanticsなので，右辺の符号を証明した後に`u^2`のpolynomial inequalityへ帰着できる場合だけ既存solverへ渡す。等式`abs[u]==a`は**明示的なReal 定義域**でのみ`u==a`または`u==-a`へ分岐する。既定Complex 定義域の`abs[z]==a`は一般に円周などのlocusであり，有限2点へ潰さず`UnresolvedSolutionSet`を維持する。

```text
solve[abs[x]<2,x]
-> {x in Real if x > -2 && x < 2}

solve[abs[x-1]>=3,x]
-> {x in Real if x <= -2, x in Real if x >= 4}

solve[abs[2x-1]<=3,x]
-> {x in Real if x >= -1 && x <= 2}

solve[abs[x]==2,x,Real]
-> {x == 2, x == -2}

solve[abs[x]==-2,x,Real]
-> {}

solve[abs[x]!=2,x,Real]
-> {x in Real if x != 2 && x != -2}

solve[abs[x]==2,x]
-> UnresolvedSolutionSet[x]
```

### 実指数函数とLambert W

Real 定義域で`a>0`と指数が実数であることを証明できれば，主値 `a^u = exp[u log[a]]`は常に正である。
このため零方程式は数値探索なしに空集合へ確定する。

```text
solve[1.1^x == 0,x,Real] -> {}
solve[1.1^x == 0,Real]   -> {}
solve[2^x == 8,x,Real]   -> {x == 3}
solve[2^(2x+1) == 8,x,Real] -> {x == 1}
solve[2^x == -1,x,Real]  -> {}
```

右辺がsolve変数を含まないconstantで，`a>0`，`a!=1`，右辺`r>0`を証明できる場合，`a^u==r`を`u==log[a,r]`へ安全に反転して既存のpolynomial solverへ渡す。branch/定義域を証明できない場合はこの変形を行わない。

Lambert W正規化層は`u exp[u]==a`を直接認識する。Real 定義域では`a>=-1/E`の`W_0(a)`，さらに`-1/E<a<0`の`W_-1(a)`を条件付きで列挙し，分岐点で重複させない。`exp[-u]==u`のような式も同じnormal formへ送る。主値 `lambertw[u]==r`は`r>=-1`を証明できる実値域で`u==r exp[r]`へ反転する。一般複素分岐 familyを推測する規則ではない。

```text
solve[x*exp[x]==1,x,Real] -> {x == lambertw[1]}
solve[exp[-x]==x,x,Real] -> {x == lambertw[1]}
solve[x+log[x]==0,x,Real] -> {x == lambertw[1]}
solve[exp[x]+x==0,x,Real] -> {x == -lambertw[1]}
solve[lambertw[x]==1,x] -> {x == E}
solve[cosh[x]==2,x,Real] -> {x == acosh[2], x == -acosh[2]}
```

従来の`a^x==x^2`分類もこのLambert W knowledgeと併存する。`L=log[a]`として，常に存在する1根を主値分岐から構成し，`|L|<=2/E`がcertifiedに成立する場合だけ負引数側の`W_0` / `W_-1` 分岐を追加する。分岐点では両分岐が一致するため重複を返さない。

```text
solve[1.1^x == x^2,x,Real]
-> {x == -2lambertw[log[11/10]/2]/log[11/10],
    x == -2lambertw[-log[11/10]/2]/log[11/10],
    x == -2lambertw[-1,-log[11/10]/2]/log[11/10]}

N[solve[1.1^x == x^2,x,Real],20]
-> {x == -0.95548727594562198165,
    x == 1.0513800237472769374,
    x == 95.71683016840522274}
```

一般の`a^(b x+c)==P(x)`，Complex全branch，Lambert Wを含む不等式はまだ一般化しない。
branch条件を証明できない場合は条件付きbranchまたは`UnresolvedSolutionSet`を保持する。

### 主値 radical equation solve

`sqrt` / `cbrt`等式は，単純に両辺を冪乗して終わらせず，**変換後多項式の候補生成**と**元radicalの値域検証**を分離する。

主値 square rootでは

```text
sqrt[A] == B
    <=> A == B^2 かつ B が principal sqrt の像に属する
```

を使う。したがって平方で導入されたextraneous rootはexact 値域 proofで除外する。主値 sqrtの像は`Re(B)>0`または`Re(B)==0 && Im(B)>=0`であり，exact numberでは直接判定する。`B`がRealと証明できる場合は`B>=0`へ縮約する。現Predicateで複素symbolic パラメータの半平面条件を完全に表せない場合は，不完全な条件へ弱めず`UnresolvedSolutionSet`を維持する。

```text
solve[sqrt[x]==2,x]       -> {x == 4}
solve[sqrt[x]==-2,x]      -> {}
solve[sqrt[x]==I,x]       -> {x == -1}
solve[sqrt[x+1]==x-1,x]   -> {x == 3}
solve[sqrt[x^2]==2,x]     -> {x == 2, x == -2}
```

`cbrt`はreal cube rootであり，`cbrt[A]==B`は`A==B^3`かつ`B in Real`と同値である。変換後の候補が高次Rational polynomialになり，右辺のreal性からsolve変数のreal性まで従うaffine caseではReal algebraic isolationを利用する。

```text
solve[cbrt[x]==2,x]         -> {x == 8}
solve[cbrt[x]==-2,x]        -> {x == -8}
solve[cbrt[x+1]==x-1,x]     -> {x == root[{-2, 2, -3, 1}, 1]}
solve[cbrt[x]==a,x]         -> {x == a^3 if a in Real}
```

また，`cbrt[z]^3 -> z`は`z in Real`が証明できる場合だけSimplifierが適用する。これはreal `cbrt`の定義域 conditionを消さないための制約である。

### 実軸周期函数のパラメータ付き solution family

`sin/cos/tan`はglobal injectiveではないため主値 inverse一個へ潰さず，実軸で周期を保った整数パラメータ familyを返す。
formal パラメータは`where k in Integer`で局所的に束縛され，relation内に`k`が既に現れる場合は`k1`等のfresh nameを選ぶ。

```text
solve[sin[x]==0,x,Real]
-> {x==Pi k where k in Integer}

solve[cos[x]==0,x,Real]
-> {x==Pi/2+Pi k where k in Integer}

solve[tan[x]==1,x,Real]
-> {x==Pi/4+Pi k where k in Integer}
```

初版は`sin/cos/tan`の1引数函数で，引数がsolve変数に対する**exact非零一次係数を持つaffine式**の場合に限定する。
targetの実値域は既存Knowledgeで検証し，`sin[x]==2`等は空集合へ落とす。
非線形argumentやComplex全解を主値 inverseだけから捏造しない。

periodはsessionの角度modeに従う。初版のperiodic Solverはargument全体への明示`Rad` / `Deg` / `Grad` suffixをまだaffine polynomialとして正規化しないため，その形は未解決のまま保持する。

```text
angleMode[Deg]
solve[sin[x]==0,x,Real]
-> {x==180k where k in Integer}
```

Complex領域でも主値 inverseだけから全解を捏造しない。

### exact代数根 `root` / `AlgebraicNumber`

一般高次多項式の根はradicalへ無理に展開せず，実数・複素数ともexactな`root`表現で保持できる。

```text
root[{a0,a1,...,an},k]
root[{a0,a1,...,an},k,Complex]
```

係数列はexact Rationalを昇冪順に並べ，

```text
a0 + a1 x + ... + an x^n
```

を表す。2引数形は**異なる実根を昇順に並べた1-based第`k`根**である。
3引数`Complex`形は全複素根をexactに分離し，`Re(z)+Pi Im(z)`の昇順による決定的1-based indexを使う。
各rootの実部・虚部は代数数でPiは超越数なので，異なる代数根がこのordering keyで一致することはない。
root isolationそのものはこのorderingの数値近似へ依存せず，Rational center/radiusを持つ一意root diskをexactなRouché判定で証明した後に順序を確定する。

定義多項式はmonicかつsquare-freeへexactに正規化する。

```text
root[{4,0,-4,0,1},2]
-> root[{-2,0,1},2]

root[{4,0,4,0,1},1,Complex]
-> root[{2,0,1},1,Complex]
```

実根では`RealAlgebraicNumber`がRational Sturm列とRational isolating intervalを保持する。
複素根では`ComplexAlgebraicNumber`がexact Rational中心・半径のisolating diskを保持し，approximate root candidateは証明のための候補生成にだけ用いる。候補初期値は係数のNewton polygonから複数のroot-radius群を推定して配置し，最終的な根の存在・一意性はexact Rouché判定だけで決定する。
`N[root[...,k],p]`は既存のcertified interval/diskを再利用し，選択rootだけを局所精製してcertified approximationへ変換する。

```text
N[root[{-2,0,1},2],30]
-> 1.41421356237309504880168872421

N[root[{1,0,1},2,Complex],30]
-> I

solve[x^5-x+1==0,x,Real]
-> {x==root[{1,-1,0,0,0,1},1]}

solve[x^5-x+1==0,x]
-> {x==root[{1,-1,0,0,0,1},1,Complex], ...}
```

既存の線形・二次・binomial・Rational-root deflation等で自然なexact式に閉じる場合は従来Solverを優先する。
それでも閉じないRational係数多項式では，Real領域はSturm real Root，Complex/default領域はcertified complex Root isolationを代替経路として使う。
現在の定義多項式次数budgetは96である。

`AlgebraicNumber`はReal/Complex Rootを共通に扱い，bounded resultant arithmeticでRootとexact Rational/complex Rationalの`+ - * /`，および小さい整数冪をexact Rootへ閉じる。
演算後はresultantの候補多項式を作るだけでなく，operandのisolating interval/diskを演算して得た保証領域と照合し，正しいresult rootが一つに証明できた場合だけ簡約する。

```text
root[{-2,0,1},2]*root[{-2,0,1},2]
-> 2

root[{1,0,1},2,Complex]+I
-> 2I
```

Rootを個別に生成する際は，square-free monic多項式のまま止めず，**選択されたrootを含む有理既約因子をexactに証明できる場合だけminimal polynomialへ縮約**する。
現在のbounded factorizationは次数16以下を対象とし，小素数体上の既約性証明とexact Kronecker factor探索を組み合わせる。
実RootではSturm根数で対象因子への所属をexactに判定し，Complex Rootではcertified isolating diskと因子側root diskの一意対応を証明する。
証明できなければ元のsquare-free定義多項式を保持し，minimal polynomialを推測しない。

```text
root[{6,0,-5,0,1},1]
-> root[{-3,0,1},1]

root[{2,0,3,0,1},1,Complex]
-> root[{2,0,1},1,Complex]
```

`isolateAll`は同一多項式に対するroot enumeration契約を守るため元の多項式とglobal root indexを維持し，個々の`root[...]`生成またはSolve出力へ変換する境界でminimal-polynomial 正準化を行う。

Root同士のfield arithmeticでは，operandのminimal polynomialがQ上既約と証明でき，候補`theta=alpha+c beta`について最初のexact線形従属から得た多項式が次数`deg(alpha)deg(beta)`を持ちQ上既約と証明できた場合，`theta`を**primitive element**として採用する。
`Q(alpha,beta)=Q(theta)`を証明できた場合だけtensor-product basisからtheta power basisへexactに変換し，和・差・積・商のminimal polynomialをsimple extension内の線形従属から直接求める。
証明できない重なり拡大やbudget超過では従来のresultant＋isolating-region再同定へ切り替える。

```text
root[{-2,0,1},2]+root[{-3,0,0,1},1]
-> root[{1,-36,12,-6,-6,0,1},2]
```

primitive-element reductionでsimple extensionを証明できた場合は，その結果を一度`root[minpoly,k]`へ表示して終わらせず，内部では`NumberFieldContext`と`AlgebraicElement`を保持する。Contextはgeneratorのminimal polynomial，選択されたReal/Complex embedding，`theta^d mod m(theta)`のreductionをimmutable sharedで保持し，Elementはpower basis上のexact Rational座標を持つ。同一Context内の後続`+ - * /`はresultantやprimitive-element探索へ戻らず，係数演算と`mod m(theta)`だけで処理する。Q上既約性を証明できる個別Rootにもgenerator fieldを付与するため，同一Rootの連続演算も元のsimple extensionを再利用できる。

```text
root[{1,-1,0,0,0,1},1,Complex]^2
-> root[{-1,1,0,-2,0,1},3,Complex]
```

Rootのuser-visible canonical formは従来どおり`root[minpoly,k]`であり，field座標は表示・structural equalityへ露出しない。`Expr::rebuildCall`はCallの引数が構造的に不変な場合だけ内部Algebraic cacheを継承し，引数が変わった場合はcacheを破棄する。したがってSimplifier・置換・制約処理等を跨いでも，同じRoot値のfield lineageを安全に維持できる。

別々のexpression lineageで構築された`NumberFieldContext`についても，**同じembedded generator identity**（同じminimal polynomial・Root 定義域・root index）ならbounded weak internerで同一immutable Contextを共有する。cacheはContextを所有せず`weak_ptr`だけを保持し，expired entryを随時除去する。最大256 entryのLRUとし，cache miss/evictionは性能にだけ影響し数学的結果には影響しない。minimal polynomialが同じでも選択embeddingが異なるRootは共有せず，異なるprimitive generatorで表された同型体・subfield関係を推測して統合することもしない。

これにより，独立に構築された同一simple extensionの要素同士もpointer-level same-field fast pathへ入れる。例えば左右が別々にprimitive-element reductionされた次の式は，再度高次数のfield constructionへ戻らず同一体内の係数乗算でexact Rootへ閉じる。

```text
(root[{-2,0,1},2]+root[{-3,0,0,1},1])
*(root[{-2,0,1},2]-root[{-3,0,0,1},1])
-> root[{1,12,-6,1},1]
```

common-field construction自体も再利用する。primitive-element reductionが成功したRoot pairについて，compositumとなる`NumberFieldContext`と両operandのpower-basis embeddingを最大64 entryのbounded cacheへ保持する。cacheはRoot identityの順序反転も認識するため，`alpha+beta`の直後の`alpha-beta`のような演算で同じtensor-product / primitive-element探索を再実行しない。output fieldは`weak_ptr`で参照し，fieldが寿命を終えたentryは随時破棄する。

Real fieldでは`AlgebraicElement`の係数多項式をchosen generatorのcertified isolating interval上でexact Rational interval評価できる。その区間が結果minimal polynomialの根をただ1つ含むことをSturm列で証明し，区間より下の根数からcanonical root indexを直接決定する。従来のように結果多項式の全実根をisolateして候補を走査し，その後`root[minpoly,k]`生成時にもう一度全根isolationを行う必要はない。証明に失敗した場合だけ従来のisolating-region再同定へ切り替える。

さらにpersistent field representationが付いた値ではminimal polynomialの既約性が既に証明済みであるため，Real次数>1がRationalへ，Complex次数>2が`Q+iQ`へ退化しないことは次数だけで分かる。その場合`exactRationalParts`は192-bit refinementを行わず即座に非退化と判定する。これらはcanonical表現やexact semanticsを変更せず，既に得た証明を再利用する性能改善である。

同一`NumberFieldContext`内の除算で使う逆元も再利用する。power-basis座標`u`の逆元`u^-1 mod m(theta)`は最初のmiss時だけextended Euclidでexactに求め，fieldごと最大16 entryのthread-safe LRUへ座標対として保持する。`inverse(inverse(u))=u`なので逆方向も同時登録し，cache hitではRational coefficient vectorのコピーだけで済む。cache evictionは再計算を増やすだけで数学的結果には影響しない。定数座標`{q,0,...}`は`Q`からの埋め込みなので，多項式Euclidを回さず`{1/q,0,...}`を直接返す。multiplication matrixの常設cacheも検討したが，12次体の実測で行列構築約228 usに対し通常乗算112 usから行列-vector 101 us程度の短縮に留まり，十分な反復回数がないと償却できないため現段階では導入しない。

開発用には次を常設している。

```text
mmCal.Benchmarks --algebraic-field [iterations]
```

これは`(sqrt[2]+cuberoot[3])*(sqrt[2]-cuberoot[3])`相当のcompositum再利用について初回とwarm平均を同一session内で測定し，併せて12次simple extension上でreciprocalのfirst/warm，warm division，minimal polynomialのfirst/warmをmicrobenchmarkする。benchmark値はcompiler / build configuration / CPUに依存するため絶対性能保証ではなく，同一環境でのregression監視に用いる。

このfield表現はexact比較にも接続する。`==` / `!=`は同一`NumberFieldContext`ならpower-basis座標の完全一致で判定し，同じcanonical Root identityは即等値，同一多項式の別root indexまたは異なるQ上既約minimal polynomialはexactに不等と証明する。それだけで決まらないbounded caseでは差を既存primitive-element/resultant経路でexactに構成し，0かをfield座標またはcertified root isolationから判定する。証明不能やbudget超過を`False`へ落とすことはない。

実代数数の`< <= > >=`は数学的orderとして実装する。同一fieldでは`a-b`のpower-basis座標をchosen real embeddingのisolating interval上でexact Rational interval評価し，0から分離した符号で判定する。異なるfieldでは各Real Rootのcertified isolating intervalを細分化し，区間が分離すればその時点でexactに順序を確定する。必要なら差のexact algebraic constructionへ切り替える。Complex Rootの`Re(z)+Pi Im(z)`による決定的orderingはroot enumeration専用であり，ユーザー数学としての`<`ではないため，Complex algebraic valueの`< <= > >=`は未評価のまま保持する。

```text
root[{-2,0,1},2] > 1
-> True

root[{-2,0,1},2] != root[{-3,0,1},2]
-> True

root[{1,0,1},1,Complex] < root[{1,0,1},2,Complex]
-> root[{1,0,1},1,Complex] < root[{1,0,1},2,Complex]
```

表面構文がRootでなくてもexact algebraic valueと証明できる式は同じ計算基盤へ接続する。現在のbridge対象はcanonical `root[...]`，exact Rational / exact complex Rational，`sqrt[q]` / `cbrt[q]`（安全に実代数数として扱えるexact Rational引数），`Phi`，およびそれらのbounded `+ - * /`・小整数冪である。これにより表示形を強制的にRootへ書き換えず，内部比較・定義域証明・Solveだけが共通`AlgebraicNumber` viewを利用する。

```text
root[{-2,0,1},2] == sqrt[2]
-> True

root[{-2,0,0,1},1] == cbrt[2]
-> True

Phi == root[{-1,-1,1},2]
-> True

element[sqrt[2]+sqrt[3],Rational]
-> False

solve[x == sqrt[2],x,Rational]
-> {}

solve[x == sqrt[2],x,Real]
-> {x == sqrt[2]}
```

bridgeにはnode数と小整数冪のbudgetを設け，変換できない式は従来経路へ戻す。近似値からminimal polynomialを推測することはない。`0^0`等の未定義形はalgebraic bridgeへ入る前に`Indeterminate`となり，`a^0=1`はbaseのexact非zero性を証明できる場合だけ使う。

resultant・primitive-element双方の次数爆発を避けるため，現在の**algebraic-field候補次数budgetは16**。完全な任意次数Q因子分解，異なるprimitive generator / subfield関係を含む一般number-field merge・正準化，証明計算基盤がsimple extensionを構成できない重なり拡大の一般reduction，未証明・非minimal表現まで含む完全なcross-context equality，現在のrefinement / algebraic construction budgetを超えるordering，`rootApproximant`はまだ未実装である。したがって現在の`AlgebraicNumber`はpersistent field lineageとbounded exact comparisonを持つ代数体計算基盤であり，完全な代数体canonicalizerではない。

---

# 26. 線形代数

## 26.1 exact-first線形代数

canonical API:

```text
transpose[A]
conjugateTranspose[A]
dot[A,B]
det[A]
inverse[A]
rref[A]
matrixRank[A]
nullSpace[A]
solveLinear[A,b]
luDecomposition[A]
qrDecomposition[A]
svd[A]
conditionNumber[A]
pseudoInverse[A]
leastSquares[A,b]
eigenvalues[A]
eigenvectors[A]
eigensystem[A]
norm[v]
normalize[v]
trace[A]
```

`dot` はrank-1/rank-2を扱う。

```text
dot[{1,2,3},{4,5,6}] -> 32
dot[{{1,2},{3,4}},{5,6}] -> {17, 39}
dot[{5,6},{{1,2},{3,4}}] -> {23, 34}
dot[{{1,2},{3,4}},{{5,6},{7,8}}] -> {{19, 22}, {43, 50}}
```

Array同士の `*` は行列積にしない。`*` はscalar×Arrayだけを許し，行列積・vector contractionは明示的に `dot` を使う。同shapeの `+/-` はelement-wise。

exact実数/Rational行列は各行の分母を払って整数行列へliftし，Bareiss fraction-free eliminationを使う。これによりpivotごとのRational生成を避ける。exact complexはflat `Number` Gaussian 計算基盤へ切り替える。symbolic行列は非零性を証明できないpivotを勝手に選ばない。

```text
det[{{1,2},{3,4}}] -> -2
inverse[{{1,2},{3,4}}] -> {{-2, 1}, {3/2, -1/2}}
rref[{{1,2},{3,4}}] -> {{1, 0}, {0, 1}}
matrixRank[{{1,2},{2,4}}] -> 1
nullSpace[{{1,2},{2,4}}] -> {{-2, 1}}
solveLinear[{{2,1},{1,-1}},{5,1}] -> {2, 1}
```

`nullSpace[A]`はRREFのfree columnを昇順に取り，各自由変数を1としたcanonical basisを返す。返り値shapeは `{nullity, columns}` であり，full column rankでは `reshape[{}, {0,n}]` として空basisのvector次元を保持する。exact整数/RationalではBareiss forward eliminationを共有し，exact complexはGaussian 代替経路，symbolicではpivotの非零性を証明できる場合だけbasisを構成する。
既に有限precisionの要素を含む行列では，pivotだけでなく「その列にpivotが存在しないこと」もInformationEnclosureから証明できる場合だけnullityを確定する。したがって`nullSpace[N[{{Pi,1}},12]]`のようにpivot構造を入力情報から証明できる場合は評価するが，`N[0,p]`を含む零列からhidden CertifiedEnclosureのexact zeroを掘り返してnullityを決めることはしない。

`solveLinear[A,b]` は `A` を m×n 行列，`b` を長さmのvectorとして扱う。一意解が存在すれば長さnのvectorを返す。正方行列に限定せず，整合した過剰決定系もfull column rankなら解ける。不整合系，または自由変数が残る系は定義域エラー。一般parametric solutionはこの函数では捏造しない。

exact整数/Rational行列は行ごとに分母を除去して整数workspaceへliftし，規模・係数height・densityに応じてBareissまたは31-bit modular 計算基盤を自動選択する。modular `det`はHadamard上界までCRTして整数を一意復元し，modular `solveLinear`はCRT後のrational reconstruction候補を元の整数系でexact verificationした場合だけ返す。bad primeや復元不能時はBareissへ切り替える。modular inverse 計算基盤も存在するが，現GCC benchmarkで測定済み範囲ではBareissが優勢なためautomatic `inverse`はBareissを維持する。

`luDecomposition[A]` は現在正方行列を対象とし，shape `{3,n,n}` の `{P,L,U}` を返す。規約は `P A = L U`。certified approximate LUでは，非零を証明できた候補のうち`|pivot|^2`の区間下限が最大の行を選ぶpartial pivotingを使い，epsilon判定は行わない。row pivotingを行い，exact NumberではRational/complexをexactに保持する。三角symbolic行列は不要な除算を行わずそのまま分解でき，非零性を証明できないpivotが必要な一般symbolic行列は未評価に留める。factorはprefix indexingで取り出せる。

```text
lu = luDecomposition[A]
at[lu,0] -> P
at[lu,1] -> L
at[lu,2] -> U
```

`qrDecomposition[A]` は矩形m×nに対応するreduced QRで，`k=min(m,n)`として `Q:m×k`，`R:k×n` を一般brace `{Q,R}` で返し，規約は `A = Q R`。factor shapeが同じ正方caseでは内部的にdense Arrayへ自動最適化されるが，ユーザー構文は同じ`{Q,R}`である。exact実数行列は，各列をprimitive整数方向へ落としてからdivision-freeな射影除去を行うfraction-free直交化を使い，平方根を直交化途中には生成せずQ/Rを最終Expr化するときだけ導入する。full-rankの先頭列群ではGram行列の対称Bareiss分解（fraction-free LDLᵀ相当）から直交整数基底を復元し，rank落ち等では直接fraction-free直交化へ切り替える。従来の3×3 hard capは廃止した。上三角/上台形caseは引き続きfast pathを使う。`N[qrDecomposition[A],p]`はexact QRを先に展開せず，実/複素ともcertified interval Householder 計算基盤へ直接入る。 QR factorの列符号は数学的に一意ではないため，exact fraction-free 計算基盤とcertified Householder 計算基盤がcomponentwiseに同じ符号を返すことは契約に含めない。保証対象は`A=Q R`とQ列の直交・正規化であり，exact 計算基盤はprimitive整数方向から決まるdeterministicな符号を使う。

```text
qr = qrDecomposition[A]
at[qr,0] -> Q
at[qr,1] -> R
```

Householder適用には複数列を一度のrow-major走査で処理できるcolumn-block kernelも実装している。ただし外部BLASを使わない現計算基盤では8/16/24次の実測で一貫した高速化が得られなかったため，自動block化は採用せずunblocked相当を既定とする。block kernelとbenchmarkは今後のBigFloat/Matrix 計算基盤最適化用に残す。

`svd[A]`はreduced SVDを `{U,S,V}` で返す。m×n入力に対して`k=min(m,n)`，`U:m×k`，`S:k×k`，`V:n×k`。実数なら `A = U S Transpose[V]`，複素数なら `A = U S conjugateTranspose[V]`。一般数値計算基盤は条件数を二乗する`A^H A`を形成せず，Householder bidiagonalizationの後にone-sided Jacobiで列を直交化する。候補factorはreconstruction residualと`U^H U` / `V^H V`の直交性を区間演算で要求表示桁より厳しく監査し，証明できなければguard digitsを増やして再試行する。exact SVDは自然に閉じる実対角等へ限定する。重複特異値の部分空間ではsingular vector basisは一意ではないため，componentごとの「唯一の真値」を主張せず，再構成・直交性を保証する。

`conditionNumber[A]`は2-ノルム条件数 `σmax/σmin` を返す。exact行列で階数落ちを証明できれば`Infinity`，特異値が1個だけの非零矩形行列は`1`，実対角行列はexactな比を返す。一般のexact非対角行列では，特異値を閉形式へ無理に展開せず未評価に留め，`N[...]`で保証付きSVDへ送る。空行列の条件数はDomainErrorである。

`pseudoInverse[A]`はMoore–Penrose擬似逆行列を返す。exactな数値行列では階数分解 `A=FG` を作り，`A^+=G^H(GG^H)^-1(F^H F)^-1F^H` をexact算術で評価するため，階数落ちしたRational・複素行列も近似へ落とさない。0×n / n×0行列では形状を転置した空行列を返す。`N[pseudoInverse[A],p]`のようにexact行列から要求精度付きで入る場合は保証付きSVDを使う。一方，既に有限precisionの要素を含む行列は，現SVD 計算基盤が入力摂動全体を特異部分空間まで保証する段階ではないため，隠れたCertifiedEnclosureから階数や特異値を復元せず保守的に未評価へ戻す。

`leastSquares[A,b]`は `A^+ b` として最小ノルムの最小二乗解を返す。`b`の長さは`A`の行数と一致しなければならない。exactな数値入力では階数落ちを含めてexactに処理する。exact行列を外側`N`から数値化する経路は保証付きSVDを使うが，行列自体が既に有限precisionなら`pseudoInverse`と同じ理由で保守的に未評価へ戻す。

```text
conditionNumber[{{3,0},{0,4}}] -> 4/3
conditionNumber[{{1,2},{2,4}}] -> Infinity
pseudoInverse[{{1,2},{2,4}}] -> {{1/25, 2/25}, {2/25, 4/25}}
pseudoInverse[{{I,0},{0,2I}}] -> {{-I, 0}, {0, -I/2}}
leastSquares[{{1,0},{0,1},{1,1}},{1,2,4}] -> {4/3, 7/3}
dimensions[pseudoInverse[zeros[0,3]]] -> {3, 0}
```

`eigenvalues[A]` / `eigenvectors[A]` / `eigensystem[A]` は正方行列の固有値・固有vector・組を扱う。`eigenvectors`の各**列**が対応する固有vectorであり，`eigensystem[A]`は `{values,vectors}` を返す。exact pathは上三角行列の対角固有値，対角行列の標準基底，およびdistinct eigenvalueを持つexact Number 2×2を明示処理する。重根を持つ非対角2×2では不足する固有vectorを複製せず未評価に留める。一般行列の `N[...]` はComplex BigFloat Hessenberg reduction + implicit shifted QRからSchur形 `A Q ≈ Q T` を求め，Schur三角行列からback substitutionで固有vectorを構成する。元入力のcertified intervalに対するSchur relationと `A v ≈ λ v` residual，およびSchur vectorのunitarityを要求表示桁より厳しく区間監査し，証明できなければguard digitsを増やして再試行する。一般非正規行列では固有値・固有vectorは摂動に敏感であり，返した各componentが唯一の真値を個別区間包含するとは主張しない。保証対象は計算されたSchur/eigenpair relationである。近接重根・defective caseで独立固有vectorを安定に構成できない場合，`eigenvectors` / `eigensystem`は推測せず未評価に留める。

`conjugateTranspose[A]`はHermitian transposeであり，complex SVDの`V^H`や複素直交性の検証に使う。rank-1では成分の共役だけを行い，rank-2では転置と共役を同時に行う。

`norm` は複素vectorに対してHermitian normを使う。

```text
norm[{3,4}] -> 5
norm[{3+4I}] -> 5
normalize[{3,4}] -> {3/5, 4/5}
```

### precision-aware `N`

FFTと同じ `ApproximationContext` / certified interval変換を共有する。したがって例えば

```text
N[dot[A,B],100]
N[det[A],100]
N[inverse[A],100]
N[rref[A],100]
N[solveLinear[A,b],100]
N[luDecomposition[A],100]
N[qrDecomposition[A],100]
N[svd[A],100]
N[conditionNumber[A],100]
N[pseudoInverse[A],100]
N[leastSquares[A,b],100]
N[eigenvalues[A],100]
N[eigensystem[A],100]
N[norm[v],100]
```

は，巨大なexact中間式を完成させてから近似するのではなく，対応するBigFloat/interval 計算基盤へ要求精度を渡して直接評価できる。exact入力を外側`N`が直接数値計算基盤へ送る場合，入力区間は要求精度に応じて再精密化できる。これに対し，行列要素が既に`DecimalApproximation` / `ComplexDecimalApproximation`なら，値計算にはCertifiedEnclosure，pivot・零／非零・階数等の判定にはInformationEnclosureを使い，入力が宣言していない情報を復活させない。`solveLinear` / `inverse` / `rref` / `matrixRank` / `nullSpace`はInformationEnclosureからpivot構造を証明できる場合だけ結果を返し，証明不能なcaseをepsilonや内部の点値で補わない。`det` / `dot` / `norm`等の連続量はCertified/Information両区間を並行伝播し，相殺時には出力Precisionを自然に下げる。現在の`luDecomposition` / `qrDecomposition` / `svd` / `conditionNumber` / `pseudoInverse` / `leastSquares` / `eigen*`は，**既に有限precisionの行列**に対する摂動保証が未完成なため，その場合は保守的に未評価へ戻す。一方，exact行列に対する`N[...,p]`では従来どおり対応する保証付き数値計算基盤を利用する。

## 26.2 Vector補助函数と互換alias

Vector向けの補助函数には，独立したcanonical函数と互換aliasの両方がある。

canonical函数:

```text
madd
vadd vsub vscalar
vcross
inner outer
vproject vangle
vmanhattan veuclidean
vreflect vreflect_axis
vsum
grad divergence curl laplacian jacobian hessian
```

互換alias:

```text
matmul mmul vdot -> dot
rank mrank       -> matrixRank
mget             -> at
vnorm vlength    -> norm
vnormalize vunit -> normalize
vdistance distance -> veuclidean
cross              -> vcross
projection         -> vproject
gradient           -> grad
singularValueDecomposition -> svd
```

aliasは別算法を持たず，同じ`BuiltinId`へ束ねる。`vadd`等のcanonical vector helperは互換名ではなく，それ自体が公開APIである。`dot`は従来どおりbilinear contraction，`inner[a,b]`は第1引数を共役するHermitian内積であり，`norm[v]`は後者と整合する。`projection[a,b]` / `vproject[a,b]`は`b inner[b,a]/inner[b,b]`を用いる。

ベクトル解析はCartesian座標を明示して使う。

```text
grad[f,{x,y,z}]
divergence[{P,Q,R},{x,y,z}]
curl[{P,Q,R},{x,y,z}]
laplacian[f,{x,y,z}]
jacobian[{f1,f2,...},{x1,x2,...}]
hessian[f,{x1,x2,...}]
```

座標指定は重複のないsymbolからなるrank-1 Arrayでなければならない。`curl`は3次元Cartesian field専用である。曲線座標系のscale factorやmetricは暗黙に仮定しない。

---

# 27. Signal processing

```text
dft[v]
fft[v]
ifft[v]
convolve[a,b]
```

Fourier位相はsessionの角度既定に依存せず，内部で明示Radian。

```text
dft[{1,2,3,4}]
-> {10,-2+2I,-2,-2-2I}

ifft[fft[{1+I,2-I,3+2I,4-3I}]]
-> {1+I,2-I,3+2I,4-3I}

convolve[{1,2},{3,4}]
-> {3,10,8}
```

exact入力では2冪長FFTはradix-2 Cooley–Tukeyを使い，forwardの公開表現も従来どおり保持する。現在のdegree budget内（16～128点）の2冪長`ifft`では，入力にFFT由来のroot-of-unity式が含まれる場合だけ，その式を`Q[t]/Phi_n(t)`のRational power-basis座標へ再埋込みしてradix-2 inverseを行う。2冪円分体では`Phi_(2^m)(t)=t^(2^(m-1))+1`を利用し，twiddle乗算を係数shiftと符号反転で処理するため，`ifft[fft[v]]`の巨大なsymbolic展開を避けられる。純粋な数値spectrumや再埋込みを証明できない式は従来経路へ戻る。5点以上の非2冪長では，入力をexact Rational/Gaussian Rationalまたは同一cyclotomic quotientの式として証明付きで写せる場合に同じ`Q[t]/Phi_n(t)`座標上でexact変換する。Gaussian Rational入力では必要に応じconductorを`lcm(n,4)`へ拡張して`I`を同じcyclotomic fieldへ埋め込む。cyclotomic degreeが現在のbudget 64を超える場合，またはsymbolic入力をfield座標へ証明できない場合は従来のgeneric exact DFTへ切り替える。通常の`fft[...]`は引き続きexact-firstであり，machine `double`へ暗黙変換しない。

`N[fft[v],p]`では`N`が第1引数を先にexact展開せず，要求精度`p`をFFTへ伝播する。FFT側は`ComplexInterval`/BigFloat端点で直接butterflyを行い，各出力成分が要求桁へ一意に丸められることを証明してから`DecimalApproximation`を返す。近似入力を含む`fft[v]`ではCertifiedEnclosureとInformationEnclosureを別々に同じ変換へ通し，相殺で隠れたguard桁を復活させない。例えば5桁入力同士の差が入力情報より小さい成分は，高精度な微小値として露出せず0中心の低Precision値として保持する。

近似FFTでは2冪長をradix-2，十分大きい非2冪長をBluestein convolutionへ還元する。小さい非2冪はdirect DFTの定数項が小さいため，現在は384点未満をdirectとしている。これはGCCの強制算法比較とMSVCの`--full`結果を合わせて選んだ保守的なpolicy値であり，数学的境界ではない。`mmCal.Benchmarks --fft-threshold 1`で環境ごとに再測定できる。

---

# 28. 乱数

乱数だけはstateful builtin。
RNG stateはKernelSessionごとに独立する。

## 28.1 seed

```text
randSeed[42]
-> 42
```

同じseedなら同じ列へ戻る。

```text
randSeed[42]
a := rand[]
randSeed[42]
rand[] == a
-> True
```

`randSeed[]`はentropyから再seedし，その再現用整数seedを返す。

## 28.2 uniform real

```text
rand[]
rand[hi]
rand[lo,hi]
```

`rand[]`は53bit dyadic lattice上のexact Rational。
例えばseed 42の先頭sampleは:

```text
randSeed[42]
rand[]
-> 227930101193189/1125899906842624
```

範囲:

```text
rand[]       : [0,1)
rand[hi]     : [0,hi), hi >= 0
rand[lo,hi]  : [lo,hi), lo <= hi
```

## 28.3 integer

```text
randint[]
randint[a]
randint[a,b]
```

```text
randint[]    : {0,1}
randint[5]   : [0,5] inclusive
randint[-5]  : [-5,0] inclusive
randint[1,6] : [1,6] inclusive
```

BigInt範囲に対応。modulo biasを避けるrejection sampling。

## 28.4 choice

```text
choice[2,3,5,7]
choice[{2,3,5,7}]
```

rank-1 arrayまたはvariadic listから1要素。

## 28.5 normal

```text
randn[]
randn[mu]
randn[mu,sigma]
```

Box-Mullerをexact dyadic uniform sampleへ適用したsymbolic expressionを返す。

```text
randn[5,0]
-> 5

N[randn[],8]
-> 例: -0.80379286
```

`randn`の角度は明示Radで，sessionの角度既定に依存しない。

**暗号用途ではない。**

---

# 29. 主要alias

| alias                    | canonical    |
| ------------------------ | ------------ |
| `pow`                    | `Power`      |
| `fact`                   | `Factorial`  |
| `fract`                  | `frac`       |
| `ln`                     | `log`        |
| `real`                   | `re`         |
| `imag`                   | `im`         |
| `mag`                    | `abs`        |
| `unit`, `csgn`           | `sign`       |
| `rect`                   | `polar`      |
| `ave`                    | `mean`       |
| `matmul`, `mmul`, `vdot` | `dot`        |
| `mtranspose`             | `transpose`  |
| `mget`                   | `at`         |
| `singularValueDecomposition` | `svd` |
| `mdet`                   | `det`        |
| `minverse`               | `inverse`    |
| `rank`, `mrank`          | `matrixRank` |
| `mtrace`                 | `trace`      |
| `mrows`                  | `rows`       |
| `mcols`                  | `cols`       |
| `mdiag`                  | `diag`       |
| `vnorm`, `vlength`       | `norm`       |
| `vdistance`              | `veuclidean` |
| `vnormalize`, `vunit`    | `normalize`  |

aliasは別実装ではなく同一`BuiltinId`へ束ねる。数学metadataやSolver規則を二重管理しない。

mmCal 1.5.0では，Mathematica互換だけを目的とした大文字始まりalias（`Sin`, `ArcTan`, `Integrate`, `Solve`等）を削除した。数学函数はlowercase canonicalを原則とする。`D`, `N`, `In`, `Out`, `Exit`, `Clear`, `Defs`, `UnDef`は記号演算・Kernel操作の固有名として例外的に維持する。互換構文が必要になった場合は，default namespaceへaliasを増やすのではなく独立したimport/compatibility層として検討する。

---

# 30. 現在のsource-callable函数一覧

現在の開発treeでは **builtin/alias登録名274個 / sourceから呼出可能な名前254個**。内部headはsource-callable数に含めない。

```text
Clear, D, Defs, DtoG, DtoR, Exit, GtoD, GtoR, In, N,
Out, RtoD, RtoG, UnDef, abs, accuracy, acos, acosh, angleMode, arg,
arrayRank, asin, asinh, at, atan, atan2, atanh, ave, beta, betaln, binom, cbrt, cases,
ceil, choice, cis, collect, cols, comb, conditionNumber, conj, conjugateTranspose, convolve, corr, corrspearman,
cos, cosc, cosh, cot, coth, cov, cross, csc, csch, csgn, curl, cv,
det, dft, diag, digamma, diff, dimensions, distance, divergence, dot, eigenvalues, eigenvectors, eigensystem, element, erf, erfc, exp, explain, expand, expc,
Ei, Si, Ci, li, polylog, fresnelc, fresnels, hypergeometric1F1, hypergeometric2F1, ellipticF, ellipticE, ellipticPi,
expm1, fact, factor, factorint, fallingfact, fft, fib, floor, frac, fract, fullSimplify,
gamma, gcd, geomean, grad, gradient, groebnerBasis, harmmean, hessian, hypot, ibeta, identity, if, ifft, im, imag, inner,
integrate, inverse, iqr, isprime, jacobian, kurtp, kurts, laplacian, lcm, leastSquares, length, lgamma, lambertw, limit, ln, log,
log10, log1p, log2, mad, madR, madd, mag, map, matmul, max, mcols,
mdet, mdiag, matrixRank, mean, median, mget, min, minverse, mmul, mod, mode,
luDecomposition, mrank, mrows, mtrace, mtranspose, nextpow2, nextprime, nintegrate, norm, normal, toNormal, normalize, nullSpace, percentile, percentrank, perm, polar, prevprime,
outer, polynomialReduce, pow, precision, prod, projection, pseudoInverse, quantile, quotient, rand, randSeed, randint, randn, range, rank,
qrDecomposition, rationalize, re, real, rect, rem, reshape, risingfact, rms, root, round, rows, rref,
sec, sech, series, sign, simplify, sin, sinc, sinh, sinhc, skew, solve, solveLinear,
singularValueDecomposition, sqrt, stddev, stddevs, stderr, sum, svd, table, tan, tanc, tanh, tanhc, trace,
totient, transpose, trigamma, trimmean, trunc, unit, vadd, vangle, var, vars, vcross, vdistance,
vdot, veuclidean, vlength, vmanhattan, vnorm, vnormalize, vproject, vreflect, vreflect_axis, vscalar,
vsub, vsum, vunit, winsor, winsorR, zeros, zscore, zeta,
bitand, bitor, bitxor, bitnot, bitshiftl, bitshiftr, bitlength, bitcount, bitget, fma, clamp, proj
```

---

# 31. エラー / Warning方針

主なCalcError種別:

- Syntax
- 定義域
- Type
- Overflow
- Name
- Evaluation
- Internal
- ResourceLimit

評価自体は成功したがalgorithmic builtinが処理を完了できない場合は，結果Exprと別にWarningを返す。`D`, `solve`, `solveLinear`, `N`, `rref`, `matrixRank`, `nullSpace`, `luDecomposition`, `qrDecomposition`, `svd`, `conditionNumber`, `pseudoInverse`, `leastSquares`, `eigenvalues`, `eigenvectors`, `eigensystem`（互換alias `rank`を含む）, `integrate`に加え，`precision/accuracy/rationalize`が対象外入力を未評価保持する場合もWarningになる。

正常な状態変更の補足にはInfo diagnosticを使う。現在は変数・函数の再定義通知が対象。

```text
D[abs[x],x]
WARN: D could not fully evaluate the derivative; unevaluated D[...] remains
-> D[abs[x], x]
```

`sin[x]`のように記号函数として保持すること自体が正しい場合はWarningにしない。

数学的な特異値は，意味が定義できる場合は`ComplexInfinity` / `Indeterminate`等の明示的exceptional valueとして保持する。定義域違反や型違反など，値として表現できない失敗はErrorにする。

例:

```text
1/0
-> ComplexInfinity

0/0
-> Indeterminate

log[0]
-> DomainError

gamma[-2]
-> DomainError

randint[5,1]
-> DomainError
```

Parser/Evaluatorはsource spanとdocumentを保持し，函数定義経由のErrorにはcall traceを付与できる。

---

# 32. 性能方針

exact/certifiedはCPUのnative doubleより大幅に重い。
過去のmicrobenchmarkでは，対象によりdouble比で約100倍〜10万倍超の差がある。

それでも対話型CLIで数十µs〜数msの処理は実用上問題になりにくいため，通常意味論をdoubleへ落として速度を稼がない。

実施済み高速化例:

- Number real-real fast path
- Rational乗除算の重複GCD除去
- Rational加算の縮約範囲最小化
- Simplifier structural key再計算削減
- Add同類項索引
- BigInt cube root Newton法
- factorial balanced product tree
- certified Log 値域 reduction / log(2) enclosure共有
- FFT radix-2

将来`for/Plot`のように数千〜数百万回の評価を行う処理では，Exact/Certifiedとは別に明示的Machine evaluatorを追加する予定。

---

# 33. 現在未実装・保留の主な項目

代表:

- 一般parametric linear system
- `hilbert`（旧仕様の名称再確認）
- 工学函数，財務函数，単位変換
- 旧colon command `:defs`, `:unset`, `:undef` 等（函数版`Defs[]/UnDef[]`は実装済み）。`:angle`は`angleMode[]`へ置換し，frontend commandとして`:help` / `:fix` / `:layout` / `:status`を実装済み
- `for`, `plot`
- general Machine/double evaluation mode

今後追加検討:

今後は任意次数の完全Q因子分解・異なるprimitive generator間の一般number-field merge / 正準化・未証明表現まで含む完全なcross-context algebraic comparison・`rootApproximant`，より一般のパラメータ付き solution family，Machine evaluator / `for` / `plot`等を候補とする。

---

# 34. 現在のCLI

CLIはKernelの数学状態とfrontendの表示状態を分離する。起動時には次のオプションを指定できる。

```text
mmCal --fix 16 --angle deg
mmCal --angle rad
mmCal --angle grad --fix 8
mmCal --layout multi
mmCal --eval "expand[(x+1)^3]"
mmCal --batch < expressions.txt
```

- `--fix n`: 起動時の小数表示桁数上限。内部値は変更せず，末尾の不要な0は省略する
- `--angle deg|rad|grad`: 起動時の既定角度
- `--layout auto|single|multi`: 通常REPLだけの出力組版。非対話modeとは併用できない
- `--eval expr`: 1式だけを非対話評価する
- `--batch`: 標準入力を1行1式として同じsessionで非対話評価する。
- `--help`, `-h`: 使用法を表示

`--layout`は対話表示専用であり，`--eval` / `--batch`との同時指定は引数errorになる。既定`auto`はTTYでは端末幅と式構造からArray/List/`cases`/解集合の改行を選び，pipe/redirect時はcanonical 1行表示へ退避する。`single`は常に1行，`multi`は対象構造を複数行へ展開する。KernelのExpr，履歴，`formatExpr()`のcanonical表現は変更しない。

非対話modeはbanner，prompt，`Out[...]` label，終了挨拶を出さない。成功値だけをstdout，Warning / Errorをstderrへ出す。終了codeは成功`0`，引数error`2`，`SyntaxError` / `ResourceLimitError`は`3`，評価errorは`4`，`InternalError`は`5`である。batchは行単位error後も続行し，最大codeをprocessの終了codeとする。`--eval`と`--batch`は同時指定できない。

CLIは実行中のtop-level評価に限って`EvaluationCancellationToken`をconsole interruptへ接続する。WindowsではCtrl-C / Ctrl-Break，POSIXではSIGINTを受けるとcancelを要求し，長時間kernelが次のcancellation pollへ到達した時点で`ResourceLimitError: Evaluation cancelled by frontend`として停止する。非対話`--eval`ではこの分類に従い終了code `3`を返す。取消しはcooperativeであり，全subsystemを非同期強制停止するものではない。

Lexer / Parserはtoken，AST node，入れ子，演算子鎖，数値literal桁，函数引数，Array要素に独立budgetを持つ。上限超過は位置情報付き`ResourceLimitError`としてAST lowering前に停止する。

```text
In [1]> 1/3
Out[1]> 1/3
```

- 1行ごとにparse/evaluate
- promptは`In [n]>` / `Out[n]>`で固定し，余分な空白を入れない
- 数学函数としての終了は`Exit[]`。CLI互換commandとして`:quit` / `:exit`も受け付ける。裸の`exit` / `quit`は特別扱いしない
- `Clear[]`: user definitionsと全履歴を消し，次の入力番号を1へ戻す
- `Defs[]`, `UnDef[...]`: user definitionsの確認・削除
- 計算履歴 `@`, `%`, `%%`, ... および正負添字を持つ再評価型`In [n]`, snapshot型`Out[n]`

## 34.1 `:help`

```text
:help
:help sin
:help functions
:help constants
:help Pi
```

引数なしではREPL commandの短い一覧を表示する。`:help 函数名`は`BuiltinRegistry`をcanonical名，callable alias，arityの正本として使う。別の利用者向けcatalogが，登録済みの**全source-callable builtin**へ個別説明，明示的な入力規則，必要なnote，一つ以上の例を与える。`D` / `integrate` / `root` / `qrDecomposition` / `svd` / `solve`等の複数形式を持つ操作は，arityだけのplaceholderではなく複数usage・複数例を表示する。内部testは`BuiltinRegistry::sourceFunctionNames()`を全列挙し，詳細項目の欠落を検出する。

`:help Pi`と`:help constants` indexは保護された数学定数，Boolean値，数値定義域，`Infinity`，`Rad` / `Deg` / `Grad`単位symbolを扱う。`:help functions`はcallable canonical名とaliasをsortして一覧する。未知topicは成功したfrontend照会のままで，明確な近傍topicがあれば一件だけ決定的に提示する。例えば隣接転置`sdv`は`svd`を，探索用短縮`qr`はcallableな`qrDecomposition`を提案するが，`qr`を新規函数aliasとして登録はしない。

`:help`は式をparse/evaluateせず，成功・定数・未知名のいずれでも`In[n]`を進めず履歴へ入らない。未知名には函数・定数indexへの案内も返す。

## 34.2 `:fix` — presentation-only小数表示

```text
:fix 16
Display: Fixed(16)

In [1]> 1/3
Out[1]> 0.3333333333333333

:fix off
Display: Exact

In [2]> Out[1]
Out[2]> 1/3
```

`:fix n`は小数点以下最大`n`桁へ丸める**表示だけ**を変更する。末尾の不要な0は省略するため，例えば`:fix 5`で`31/10`は`3.1`と表示する。保存されるExpr，`Out[n]`，`precision/accuracy`の意味論は変更しない。exact値をMachine/doubleへ変換する機能ではない。`n`は現在0..1000。`:fix`のみなら現在の表示modeを表示する。

数値としてcertifyできる式全体は表示時だけ近似する。自由変数を含むsymbolic expressionはexact表記を維持する。

## 34.3 `:layout` — 対話REPLの組版

```text
:layout
Layout: Auto

:layout single
Layout: Single

:layout multi
Layout: Multi
```

`:layout`は通常REPLの**表示上の組版だけ**を切り替える。`single`は既存のcanonical 1行表現をそのまま使う。`multi`はArray/List/`cases`/有限・条件付き解集合を構造的に改行する。`auto`はTTY上で端末幅と式構造から選択し，stdoutがpipe/redirectなら1行表示へ退避する。式の内部表現，`Out[n]`，再parse可能なcanonical formatter，自動処理出力には影響しない。

例えば`multi`では，

```text
Out[1]> {
          {1, 2},
          {3, 4}
        }

Out[2]> cases[
          x^2 if x >= 0;
          -x if x < 0;
          0
        ]
```

`:layout`だけなら現在modeを表示する。

## 34.4 `:status`

```text
:status
Angle: Rad
Display: Exact
Layout: Auto
Evaluation: Exact-first
Definitions: 0
History: 0
```

`:status`も履歴へ入らないCLI commandである。数学状態の変更は`angleMode[...]`等のKernel函数，frontendの照会・presentation状態変更は`:help` / `:fix` / `:layout`等のCLI command，という境界を維持する。

## 34.5 `:quit` / `:exit`

```text
:quit
:exit
```

どちらも現在のCLI sessionを正常終了する互換commandであり，式としての`Exit[]`と同じ目的を持つ。parse/evaluateや履歴追加は行わない。裸の`quit` / `exit`は通常の入力として扱い，特別な終了commandにはしない。`--batch`でも認識し，その行で正常終了する。

## 34.6 console title

タイトルは補助情報として，例えば次の形式に更新する。

```text
mmCal <version> - Rad - Exact - Layout(Auto)
mmCal <version> - Deg - Fixed(16) - Layout(Multi)
```

- Windows: `SetConsoleTitleA`
- Linux/macOS: TTY時のみANSI OSC title sequence
- その他: no-op

タイトル変更失敗は計算Errorにしない。状態確認の正本は`:status`であり，terminalがtitleを上書きしても意味論には影響しない。

## 34.7 canonical formatter

通常`Out[n]`はAST dumpではなく，再parse可能なcompact数学表記とする。

```text
x^2+sin[x]
A-B+C
2(x+sqrt[x])sqrt[x+sqrt[x]]/3
```

- `+ -` / `+-`は出さず，負項を`-`として表示
- `A-(B-C)`のような加減算は表示時だけ`A-B+C`へflattenできる
- `+`, `-`, `*`, `/`, `^`の前後には不要な空白を置かない
- 比較演算子`==`, `!=`, `<`, `<=`, `>`, `>=`はrelationを読みやすくするため前後に1空白を置く
- implicit multiplicationは字句上安全な場合だけ連結する（`2x`, `2sqrt[x]`）。`2exp[x]`や`2E`のように指数表記と衝突する連結は`2*exp[x]`, `2*E`と明示する
- identifier同士など連結で別tokenになる場合は必要な空白を残す（`I Pi`, `x y`）
- 数字同士など曖昧になる場合は空白ではなく明示`*`を使う
- precedence/associativityを守り，format → parse → formatで意味が変わらないことを回帰テストする

canonical formatterは1行serialization契約の正本として維持する。対話REPLの`:layout`組版はこれとは別のpresentation layerであり，canonical表現を変更しない。内部構造を見せるdebug/full-form表示も通常formatterとは将来別機能に分離する。

---

# 35. 実装上の主要層

```text
Lexer / Parser
    ↓
Lowerer
    ↓
Expr AST
    ↓
Evaluator explicit task stack
    ↓
Builtin / user function / symbolic operation
    ↓
Simplifier + MathRegistry + KnowledgeContext
    ↓
Exact result
      or
CertifiedEvaluator -> interval -> DecimalApproximation
```

主要責務:

- `SymbolTable`: intern / identity
- `SymbolRegistry`: protected symbol / constant / 定義域 name
- `BuiltinRegistry`: name / alias / arity / Hold属性
- `MathRegistry`: 定義域 / parity / branch / definedness
- `ValueFacts`: conservative numeric-定義域/sign inference
- `KnowledgeContext`: permanent facts + assumptions
- `Simplifier`: safe local rewrite
- `FullSimplifier`: bounded candidate search
- `CertifiedEvaluator`: expression全体のinterval evaluation
- `SolutionSet`: Solverの解集合表現
- `RandomEngine`: session-local stateful PRNG

この分離を維持し，函数追加ごとにSolver・Simplifier・数値計算基盤へ同じ知識を重複記述しないことを基本方針とする。

## CertifiedEvaluator の安全限界

CertifiedEvaluatorは式を再帰的に区間評価するため，病的に深いASTについてはOSのstack overflowへ到達する前に評価対象外として扱う。現在の深さ上限は96段。通常のn-ary `Add` / `Multiply` の項数ではなく，ASTの入れ子深さに対する安全弁である。

`nintegrate` は高階導函数を作る前に被積分函数を区間全体でpreflightし，明白な特異点を先にDomainErrorへ落とす。
