# Mathematica 1.2.2f33 Enhanced 静的構造解析報告

**対象アーカイブ:** `Mathematica_1.22.zip`  
**実体バージョン:** `1.2.2f33 Enhanced`  
**解析方式:** Classic Mac OSファイル構造、MacBinary/resource fork、68k `CODE` resource、可読ASCIIシンボル、同梱Mathematica言語ソース（`.m`）の静的解析  
**SHA-256:** `b0cd32505d41b32f9b36715b7e5615b8abfd86724ef9fa81ff49742803dc5afd`

---

## 0. この報告書の範囲と確度

本報告書は、対象アーカイブを**実行せず**に解析した結果をまとめたものである。主な証拠は次の4種類である。

1. **MacBinaryヘッダとresource fork**
   - ファイルタイプ、creator、data/resource fork長、`vers`、`SIZE`、`CODE`等を直接解析した。
2. **68kコードresource中の可読シンボル**
   - `oIntegrate`, `oSolve`, `BuildDispatch`, `UseDispatch`, `NewBignum`等を抽出した。
   - これは内部機能の存在を示す強い証拠だが、可読名だけから低レベルアルゴリズムまで断定はしない。
3. **同梱 `.m` ソース**
   - パターン規則、パッケージ、起動処理、積分表、Series、統計、数値近似などはソースそのものを読める。
4. **Front End resourceと文字列**
   - Notebook/Cell、local/remote kernel、通信、印刷、PostScript等の構造を確認した。

本文では、必要に応じて次の意味で表現を使い分ける。

- **確認:** ファイルまたはソースに直接存在する。
- **強い推定:** 複数のバイナリ/resource/source証拠からほぼ確実だが、命令単位の逆アセンブルは未実施。
- **未確定:** 68k逆アセンブルまたは実行トレースが必要。

したがって本報告書は「完全な逆コンパイル」ではない。しかし、**システムの責務分割、評価言語の構造、主要CAS機能の実装層、起動構造、rule layerの設計**についてはかなり深いところまで確認できる。

---

# 1. アーカイブの実体

## 1.1 バージョン

`Mathematica.bin`, `Math A.bin`, `Math B.bin`, `Math C.bin`, `Mathematica Prefs.bin` の `vers` resourceはいずれも、

```text
1.2.2f33 Enhanced
Copyright Wolfram Research, Inc. 1988-90
```

を保持している。

したがってファイル名の `1.22` は、実体としては **Mathematica 1.2.2 build f33 Enhanced** を指す。

MacBinaryの作成・更新時刻も1988～1990年に集中しており、Classic Mac版の実体を保存したものと判断できる。

## 1.2 Classic Macのfork構造

トップレベルには次のような組が存在する。

```text
Mathematica            0 bytes
Mathematica.bin        678,656 bytes

Math A                 0 bytes
Math A.bin             634,112 bytes

Math B                 0 bytes
Math B.bin             768,256 bytes

Math C                 0 bytes
Math C.bin              60,416 bytes
```

0 byte側は「壊れたファイル」ではない。`Mathematica`, `Math A/B/C` は実質的にresource fork主体のClassic Macファイルであり、`.bin`側がMacBinary IIIとしてresource forkとFinder metadataを保存している。

確認したタイプコードは次の通り。

| ファイル | Macintosh type | creator | data fork | resource fork |
|---|---|---|---:|---:|
| `Mathematica.bin` | `APPL` | `OMEG` | 0 | 678,522 |
| `Math A.bin` | `OCDA` | `OMEG` | 0 | 633,884 |
| `Math B.bin` | `OCDB` | `OMEG` | 0 | 768,057 |
| `Math C.bin` | `OCDC` | `OMEG` | 0 | 60,275 |
| `Mathematica Help.bin` | `TEXT` | `OMEG` | 297,713 | 281,709 |
| `Mathematica Prefs.bin` | `OMPF` | `OMEG` | 0 | 11,528 |

`OCDA/OCDB/OCDC` はこのアプリケーション固有の外部コード容器と考えるのが自然である。

## 1.3 同梱ソースの規模

`Packages` 以下で可読な `.m` は **68ファイル、13,690行、531,029 bytes**。行数はClassic MacのCR改行をLFへ正規化して数えた。

ディレクトリ別では次の通り。

| 領域 | ファイル数 | 行数 | bytes |
|---|---:|---:|---:|
| `Algebra` | 4 | 484 | 14,380 |
| `Calculus` | 5 | 1,210 | 38,353 |
| `DataAnalysis` | 6 | 2,548 | 97,408 |
| `DiscreteMath` | 5 | 410 | 9,750 |
| `Examples` | 10 | 511 | 18,576 |
| `Geometry` | 2 | 242 | 5,663 |
| `Graphics` | 7 | 1,577 | 48,664 |
| `LinearAlgebra` | 2 | 72 | 1,405 |
| `Miscellaneous` | 2 | 569 | 13,023 |
| `NumberTheory` | 3 | 103 | 2,197 |
| `NumericalMath` | 4 | 1,604 | 56,748 |
| `StartUp` | 14 | 4,008 | 211,580 |
| `Utilities` | 2 | 56 | 1,277 |
| `init.m` | 1 | 86 | 4,648 |
| `sysinit.m` | 1 | 210 | 7,357 |

特に大きいものは、

- `StartUp/Series.m`: 1,130行
- `NumericalMath/Approximations.m`: 987行
- `Calculus/DefiniteIntegrate.m`: 951行
- `DataAnalysis/ContinuousDistributions.m`: 798行
- `StartUp/IntegralTables.m`: 798行
- `StartUp/info.m`: 654行
- `DataAnalysis/DiscreteDistributions.m`: 597行
- `DataAnalysis/ConfidenceIntervals.m`: 582行

である。

この分布だけでも、1.2.2の機能の多くが単なる68k builtinではなく、**Mathematica言語自身で書かれた上位レイヤ**だったことが分かる。

---

# 2. 全体アーキテクチャ

静的解析から見える全体像は、おおむね次のようになる。

```text
┌────────────────────────────────────────────┐
│ Macintosh Front End : Mathematica APPL     │
│                                            │
│ Notebook / Cells / Selection / Files       │
│ Help / Printing / PostScript / Image       │
│ Eval / MathTalk / Telecom / tcp            │
└──────────────────────┬─────────────────────┘
                       │
          local kernel │ remote kernel
                       │
        ┌──────────────▼──────────────┐
        │  Kernel code segments       │
        │  Math A / Math B / Math C   │
        │  111 CODE resources         │
        └──────────────┬──────────────┘
                       │
        ┌──────────────▼──────────────┐
        │ symbolic evaluation runtime │
        │ definitions / patterns      │
        │ attributes / rules          │
        │ numeric tower / CAS         │
        └──────────────┬──────────────┘
                       │
        ┌──────────────▼──────────────┐
        │ Mathematica-language layer  │
        │ StartUp/*.m                 │
        │ Algebra / Calculus / ...    │
        │ packages                    │
        └─────────────────────────────┘
```

重要なのは、**Front Endとkernelの分離がすでに明確**なことである。

---

# 3. Front End — `Mathematica.bin`

## 3.1 resource構成

`Mathematica.bin` のresource forkには **471 resources / 26 resource types** が存在する。

主なものは、

| type | count |
|---|---:|
| `CODE` | 59 |
| `STR ` | 159 |
| `STR#` | 11 |
| `DITL` | 48 |
| `DLOG` | 36 |
| `MENU` | 25 |
| `PICT` | 12 |
| `WIND` | 6 |
| `ICON` | 6 |
| `ICN#` | 9 |

である。

`CODE` resourceの総量は **512,142 bytes**。

## 3.2 CODE resourceの名前から見える責務

Front End側には次の名前付きCODE resourceがある。

```text
Main
About
Anim
Post2 / Post3 / Post
Files / Files2 / GetFile / PutFile
Cells
Cline
Complete
MathTalk
Eval
ViewPt
Find / Find2
Help
Image
Validate
Sel / Sel2
Notebook
Printing / PrintOps
Telecom
tcp
```

この構成は非常に明瞭である。

- `Notebook`, `Cells`, `Sel`: 文書・セル・選択範囲
- `Eval`, `MathTalk`: kernel評価との橋渡し
- `Telecom`, `tcp`: remote kernel通信
- `Post`, `Post2`, `Post3`, `ConvPS`: PostScript生成・変換
- `Printing`: 印刷
- `ViewPt`: 3D viewpoint UI
- `Complete`: 補完
- `Help`: help UI
- `Image`: 画像
- `Files`: ファイル管理

すなわちFront Endは単なる端末ではなく、**Notebook editor + graphics/printing + kernel session manager**である。

## 3.3 local / remote kernel

Front End内の文字列には明示的に、

```text
No Kernel
Local Kernel
Remote Kernel
```

があり、さらに、

```text
Do you really want to quit the local kernel?
Do you want to quit the remote kernel?
Waiting for connection to remote host.
Automatically start local kernel
```

等が存在する。

またCODE resourceに`MathTalk`, `Telecom`, `tcp`がある。

これは、1.2.2時点ですでに、

```text
Notebook Front End
   ├─ Local Kernel
   └─ Remote Kernel over network
```

という構造を持っていたことを直接裏付ける。

## 3.4 外部kernelファイルのロード

Front End内には、

```text
Please locate the file "Math A".
If you can't find this file, you won't be able to evaluate any expressions ...
```

という文字列が `Math A/B/C` それぞれについて存在する。

したがって `Math A/B/C` は付属データではなく、**ローカル評価器を成立させる必須コードセグメント**である。

## 3.5 メモリ要求

`SIZE` resourceは、

- preferred: **3,145,728 bytes = 3 MB**
- minimum: **524,288 bytes = 512 KB**

を指定している。

当時のMacintoshアプリとして、Front Endとローカルkernelをひとつのセッションで扱うための比較的大きなメモリモデルだったことが分かる。

---

# 4. Kernel segmentation — `Math A / B / C`

## 4.1 68k CODE resource

| ファイル | CODE数 | CODE総bytes |
|---|---:|---:|
| `Math A.bin` | 51 | 632,230 |
| `Math B.bin` | 56 | 766,308 |
| `Math C.bin` | 4 | 59,722 |
| **計** | **111** | **1,458,260** |

Front Endの512 KBとは別に、約1.46 MBの68kコードが外部kernelセグメントとして存在する。

## 4.2 `Math A`

可読シンボルから、主に次の領域が確認できる。

- polynomial / algebra:
  - `oFactor`, `oApart`, `oTogether`, `oResultant`
  - `oPolynomialGCD`, `oPolynomialDivision`, `oCollect`
- equation solving:
  - `oSolve`, `oMainSolve`, `oFullSolve`, `oSolveAlways`
  - `oEliminate`, `oAlgebraicRules`
- series / calculus:
  - `oSeries`, `oSeriesData`, `oSeriesCoefficient`
  - `oLimit`, `oResidue`, `oD`
- arbitrary precision:
  - `oN`, `oPrecision`, `oAccuracy`
  - `oSetPrecision`, `oSetAccuracy`
- elementary/special functions:
  - trig/hyperbolic/inverse trig
  - `Gamma`, `Beta`, `Zeta`, `PolyGamma`, `PolyLog`
- numerical/symbolic utility:
  - `NRoots`, `NSum`, `NProduct`
- plotting:
  - `Plot`, `Plot3`, `ParametricPlot`, `ContourPlot`, `DensityPlot`

大雑把には**symbolic algebra / solver / series / precision / plotting**の比重が高い。

## 4.3 `Math B`

`Math B` はさらに言語runtime色が強い。

### 評価・定義

```text
oSet
oSetDelayed
oUpSet
oTagSet
oUnset
oClear
oClearAll
oAttributes
oCondition
oHold
oRelease
```

### rule / traversal

```text
oReplaceAll
oMap
oMapAll
oMapAt
oApply
oScan
oCases
oSelect
oPosition
oMatchQ
oFixedPoint
oNest
```

### control flow

```text
oBlock
oFor
oWhile
oWhich
oSwitch
oCatch
oThrow
oCheck
```

### calculus

```text
oIntegrate
oNIntegrate
oDSolve
```

### linear algebra

```text
oDet
oInverse
oLinearSolve
oEigensystem
oSingularValues
oPseudoInverse
oRowReduce
oNullSpace
oMatrixPower
oMatrixExp
oLatticeReduce
```

### discrete / number theory

```text
FactorInteger
ExtendedGcd
JacobiSymbol
MoebiusMu
DivisorSigma
Divisors
PowerMod
ProbablePrimeQ
```

### special functions

Hypergeometric, Legendre, Jacobi, Hermite, Laguerre, elliptic functions等。

このため `Math B` は、**一般評価器・rule runtimeと多数の数学builtinが同居する主要kernelセグメント**とみられる。

## 4.4 `Math C`

`Math C` は約60 KBしかなく、可読operator名も少ない。

```text
oAiryAi
oBesselI
oBesselJ
oBesselK
oBesselY
oNumberQ
oPlus
oMinus
oSubtract
```

したがって、**Airy/Bessel系の専門数値コードを分離したセグメント**である可能性が高い。

`Plus`等が同居する理由は逆アセンブルなしには確定できないため、単純に「算術kernel」とは断定しない。

---

# 5. 起動系 — `sysinit.m` と `init.m`

## 5.1 `sysinit.m`

`Packages/sysinit.m` 自身が、kernel起動時に必ず評価され、これがないとkernelは起動しないと記している。

実際のboot sequenceは概ね次の通り。

```text
Begin["System`"]

General::noopen
General::writewarn

Needs["StartUp`ValueQ`"]
Needs["StartUp`Formats`"]
Needs["StartUp`Attributes`"]
Needs["StartUp`Digits`"]
Needs["StartUp`GroebnerBasis`"]
Needs["StartUp`InverseFunctions`"]
Needs["StartUp`LinearProgramming`"]

End[]
```

その後Front End専用の `FE`` contextを作り、

- `$Display = "stdout"`
- completion用 `FE`FC`
- template作成用 `FE`FT`
- list対応 `SetOptions`

等を定義する。

さらにMac版では`Quit`/`Exit`を直接kernel終了に使わせず、Front EndのFile menuを使うよう差し替えている。

最後に、

```text
$Path = Join[$Path, {StringJoin[First[$Path], ":StartUp"]}]
<<init.m
```

としてユーザー初期化へ移る。

## 5.2 `init.m`

`init.m`はユーザー向け初期化ファイルで、

- `$Path`追加
- default ViewPoint
- 任意のユーザー定義

を保存する場所になっている。

## 5.3 `.m`でありNotebookでもある

`sysinit.m`や`init.m`には、

```text
fontset = ...
:[font = ...]
```

といったFront Endのcell/style metadataがコメント領域に大量に含まれている。

つまり同じファイルが、

- kernelからは**Mathematica source**
- Front Endからは**整形されたNotebook-like文書**

として機能する。

さらに14個の`.m`には`.m.bin`も同梱され、data forkのsourceに加えてresource fork側のMac固有情報まで保存されている。

これは当時の「コードと文書の融合」が単なる理念ではなく、ファイル形式にも表れている。

---

# 6. 評価言語の中核構造

この版を理解する上で最も重要なのは、数学builtin一覧ではなく**評価runtime**である。

同梱source全体では、おおよそ次の構文が使われている。

| 構文 | 出現数 |
|---|---:|
| `:=` | 1,172 |
| `/;` | 454 |
| `:>` | 338 |
| `/:` | 461 |
| `Block[` | 380 |
| `Function[` | 63 |
| `While[` | 42 |
| `Dispatch[` | 10 |

数値は単純文字列countであり構文木解析ではないが、設計傾向を見るには十分である。

## 6.1 Set / SetDelayed / UpSet / TagSet

kernelには、

```text
oSet
oSetDelayed
oUpSet
oTagSet
oTagUnset
```

が存在する。

つまり定義は単なる「symbol tableへの値代入」ではなく、**式patternに対するruleをsymbolへ関連付ける仕組み**としてkernelの中心にある。

後世の`OwnValues`, `DownValues`等の公開API名は、このbuildの`info.m`には見当たらない。一方バイナリには、

```text
ValueCell
CreateValue
FreeValue
CopyValueList
```

があり、内部には明確なrule/value storageが存在する。

この版では、ユーザー側の introspection は主に `Definition[symbol]` として提供されている。

## 6.2 pattern engine

`info.m`には次のpattern primitivesが明記されている。

- `Blank`
- `BlankSequence`
- `BlankNullSequence`
- `Pattern`
- `PatternTest`
- `Optional`
- `Repeated`
- `RepeatedNull`
- `Condition`
- `Rule`
- `RuleDelayed`

バイナリ側にも、

```text
scanPattern
scanBlank
scanBlankSequence
CheckBlank
gmatch_pattern
pPattern
pBlank
CheckCondition
```

などがある。

したがってpattern matchingはlibrary emulationではなく**kernel nativeの評価機構**である。

## 6.3 Attributesがpattern matchingに介入する

`info.m`は次のattributeを明示する。

- `HoldFirst`
- `HoldRest`
- `HoldAll`
- `Flat`
- `OneIdentity`
- `Orderless`
- `Listable`
- `Protected`
- `Locked`
- `ReadProtected`
- `Constant`

特に`Flat`, `OneIdentity`, `Orderless`について、usage text自身が「pattern matchingで考慮される」と記している。

これは重要で、pattern matcherが単なるtree wildcard matchingではなく、

```text
head attributes
   ↓
associativity / commutativity / identity-like behavior
   ↓
pattern match
```

を扱うことを意味する。

CASとしての強さはこの層に大きく依存している。

## 6.4 Holdと評価制御

kernelには`oHold`, `oRelease`, `Condition`, `SetDelayed`等があり、source側でも`HoldFirst/HoldAll`が使用される。

例として、

```mathematica
Attributes[Movie] = {HoldFirst};
Attributes[EditDef] = HoldAll;
Attributes[ValueQ] = {HoldAll};
```

等が確認できる。

したがって「函数を呼ぶ前に全引数を評価する」という単純なcall semanticsではなく、**headごとにargument evaluation policyが変わる**。

## 6.5 Dispatch — rule tableのコンパイル

公開関数として、

```text
Dispatch[{lhs1->rhs1, lhs2->rhs2, ...}]
```

があり、usageには「optimized dispatch table representation」と明記されている。

さらにバイナリには、

```text
CreateDispatch
BuildDispatch
BuildDispatch1
UseDispatch
UpdateDispatch
```

が存在する。

`Algebra/Trigonometry.m`では実際に、

```mathematica
TrigCanonicalRel = Dispatch[TrigCanonicalRel]
TrigCanonical[e_] := e //. TrigCanonicalRel
```

という形で使われている。

これは非常に重要である。

**数学知識は可読rule listとして記述し、実行時には専用dispatch表へ変換する**という二段構成が、すでに1.2.2で成立している。

## 6.6 Block中心のスコープ

同梱sourceでは`Block[`が約380回使われる。

一方、少なくとも本buildのusage table、kernel operator文字列、同梱sourceには`Module`が確認できない。

したがってこの時代のパッケージ実装は、

- `BeginPackage`
- ``Private` context`
- `Block`

を中心に構成されており、後年一般化するlexical local symbol生成よりも、**context分離 + dynamic localization**が主要技法だったと考えられる。

---

# 7. Package / Contextシステム

同梱sourceには約40個の`BeginPackage[...]`が存在する。

典型構造は、

```text
BeginPackage["DataAnalysis`ContinuousDistributions`",
             "DataAnalysis`DescriptiveFunctions`",
             "NumericalMath`InverseStatisticalFunctions`"]

publicSymbol::usage = ...

Begin["`private`"]

implementation...

End[]
EndPackage[]
```

である。

この構造には、

1. public context
2. dependency context
3. private implementation context
4. usage metadata

が一体化している。

つまりpackage systemは単なるファイル読み込み規約ではなく、**symbol namespaceそのものがmodule system**である。

---

# 8. documentation / message metadataもsource

## 8.1 `info.m`

`StartUp/info.m`には、

- 648個の`::usage` assignment
- 647 unique symbols

を確認した。

graphicsからpattern、solver、numerics、I/O、control flowまで、大半のSystem symbolのusageがここへ集約されている。

## 8.2 `msg.m`

`StartUp/msg.m`には、

- 441 message assignments
- 439 unique `symbol::tag`

を確認した。

たとえば、

```text
Solve::ifun
Integrate::...
Set::...
Graphics::...
```

のような診断が数学コードと分離されている。

つまりこの版ですでに、

```text
machine code / algorithm
symbolic definitions
usage documentation
messages
```

を別レイヤへ分ける思想がある。

---

# 9. 数値システム

## 9.1 arbitrary precisionはkernel native

バイナリには次の内部名が確認できる。

```text
NewBignum
DoubleBignum
TimesBignumBignum
HashBignum

MInteger
MReal
MIntegerQ
MRealQ

bignumbits
bignumaccubits
bitstodigits
digitstobits

RaisePrecision
LowerPrecision
RaiseAccuracy
LowerAccuracy
```

公開側には、

```text
N
Precision
Accuracy
SetPrecision
SetAccuracy
WorkingPrecision
AccuracyGoal
Rationalize
Chop
```

が存在する。

従って任意精度計算は、上位packageで多桁decimalを模倣しているのではなく、**kernel内部のnumber representationと評価dispatchに組み込まれている**。

## 9.2 型組合せ別演算

バイナリには、

```text
PlusII
PlusIF
PlusIR

TimesII
TimesIF
TimesIR
TimesRR
TimesCC
```

のような名前が残る。

I/F/R/Cの正確な内部type対応は逆アセンブルしない限り断定しないが、少なくとも**operand型の組合せごとに算術経路を分けるdispatch**が存在した可能性が極めて高い。

## 9.3 ここからは分からないこと

静的文字列だけでは、

- bignum limb幅
- radix
- multiplication algorithm
- division algorithm
- GCD algorithm
- arbitrary precision transcendental functionの具体的手法
- guard digitsの厳密な管理規則

までは確定できない。

これらは68k CODE resourceの逆アセンブル対象である。

---

# 10. Algebra / Polynomial subsystem

kernel内にはかなり充実した多項式処理が存在する。

```text
Factor
FactorTerms
Apart
Together
Cancel
Coefficient / Coefficients
Exponent
PolynomialQ
PolynomialDivision
PolynomialGCD
PolynomialQuotient
PolynomialRemainder
Resultant
Variables
Decompose
```

さらにmodular / number theoretic pathも確認できる。

## 10.1 `GroebnerBasis` は薄いwrapper

興味深いことに、`StartUp/GroebnerBasis.m` は20行程度しかない。

実装は概念的に、

```text
GroebnerBasis
   ↓
AlgebraicRules[..., vars]
   ↓
AlgebraicRulesData
   ↓
内部fieldからbasisを取り出す
```

という構成である。

つまりGroebner basisを独立した巨大algorithmとして公開APIへ直結させるのではなく、**一般的なalgebraic rule生成器の結果を再利用している**。

`AlgebraicRulesData`はopaque capsuleに近い内部表現であり、`GroebnerBasis.m`ではその第6要素へアクセスしている。

内部データ形式の厳密な意味までは未確定だが、solverとGroebner処理が共通代数基盤を共有していることは強く示唆される。

---

# 11. Solve subsystem

## 11.1 kernel側の入口

`Math A`には、

```text
oMainSolve
oFullSolve
oSolve
oSolveAlways
SolveInit
TakeSolveOptions
TryLinearSolve
```

等がある。

`info.m`では`MainSolve`を、SolveとEliminateが呼び出すunderlying functionと明記している。

したがって構造は概ね、

```text
Solve / Eliminate
      ↓
   MainSolve
      ↓
 algebraic transformation
 elimination / directives
      ↓
 rules / logical result
```

と読める。

## 11.2 transcendental solutionは意図的に不完全性を通知

`msg.m`には、

```text
Solve::ifun =
"Warning: inverse functions are being used by Solve,
so some solutions may not be found."
```

がある。

また、

```text
Solve::tdep
Solve::dinv
```

等があり、

- transcendental dependence
- inverse functionsを適用できない複数argument依存

を区別している。

これはsolverが「何でも代数的に解く」のではなく、

```text
algebraic core
   ↓
必要なら inverse-function heuristic
   ↓
incompleteness warning
```

という層を持つことを示す。

## 11.3 inverse function table

`StartUp/InverseFunctions.m`には、

```text
Sin ↔ ArcSin
Cos ↔ ArcCos
Tan ↔ ArcTan
Log ↔ Exp
```

だけでなく、Jacobi elliptic functionsのinverse mappingまで定義されている。

つまりinverse knowledgeの一部はkernel hard-codeではなく**symbolic tableとして外出し**されている。

---

# 12. Integrate subsystem

同梱ソースで最も興味深い部分の一つである。

`StartUp/IntegralTables.m`は798行あり、

- `:>` 約253
- `/;` 約151

を含む。

単なる積分公式の羅列ではない。

## 12.1 rule groupを段階化

主なrule tableは、

```text
IntBase
ExpQuad
ExpLinear
ExpOther

AlgBase
AlgCase2
AlgCase3

PureTrigInt
TrigInt

IntOut
IntIn
IntToTan
IntFromTan
TrigExp

IntMatchIn
LogSimplify
ExpSimplify
IntMatchInt
IntMatchOut

IntAlg

TrigToComplexRel
ComplexToTrigRel
TrigCanonicalRel
```

である。

これは積分器が、

```text
normalization
→ family recognition
→ specialized integration rules
→ transformations
→ fallback pattern integrator
→ output normalization
```

という複数段のpipelineを持つことを示す。

## 12.2 rule tableの保護

冒頭に、

```mathematica
Attributes[Lock] = HoldAll
Lock[x_] := Attributes[x] = {Protected, Locked}
```

があり、各table構築後に`Lock[...]`している。

つまり数学規則は可変なMathematica expressionとして構築される一方、完成後は**Protected + Lockedで固定**する。

可読性と実行時安全性を両立する発想である。

## 12.3 三角積分

`PureTrigInt`には、

- 基本三角函数
- sec/csc
- tan/cot
- sin/cosの正負整数冪
- product patterns

が明示的に存在する。

たとえば正の偶数冪`Sin[X]^n`, `Cos[X]^n`にはdouble-factorialを使った一般ruleがあり、個々の指数を列挙していない。

## 12.4 affine argument lifting

特に興味深いruleが `TrigInt` にある。

概念的には、

```text
Int[f[a X + b]^n]
   ↓
Int[f[X]^n]
   ↓
内部trig rulesで処理
   ↓
X → a X+b
   ↓
1/a
```

とする。

ソースにはその直前に、歴史的に実に味わい深い、

```text
(* tags don't work arrgh! *)
```

というコメントまで残っている。

これは、個別patternを増殖させる代わりに、**argumentのaffine structureを一度標準形へ持ち上げて再利用する**設計である。

## 12.5 half-angle rationalization

`IntToTan` / `IntFromTan`では、

```text
T = Tan[X/2]
```

として、

```text
Sin[X] → 2T/(1+T^2)
Cos[X] → (1-T^2)/(1+T^2)
Tan[X] → ...
```

へ変換する。

つまり三角積分の一部は、**Weierstrass substitutionでrational problemへ落とす**経路を明示的に持つ。

## 12.6 特殊函数への着地

`IntMatchInt`には、

- `LogIntegral`
- `PolyLog`
- `ExpIntegralEi`
- `Erf`
- `SinIntegral`
- `CosIntegral`

等へ帰着する非初等積分ruleがある。

従ってintegratorは「elementary antiderivativeだけ」を対象にしていない。

## 12.7 recursive integration

積分規則の多くは右辺で再び`Int[...]`を呼ぶ。

例:

```text
x^n e^x
log[x]^n
x^m log[x]^n
exp(ax) sin(bx)^n
```

など。

つまりrule tableは単なるlookupではなく、**漸化式として実行されるterm-rewriting program**である。

## 12.8 performance hazardを作者自身が管理

ソースには、

```text
CAUTION: Patterns with a head of Plus ... can greatly slow the integrator
```

とあり、実際に一部ruleをコメントアウトしている。

さらに`TrigCanonicalRel`には、

```text
CAUTION: This rule can lead to an infinite loop
```

として無効化されたruleが複数ある。

rule systemの問題はすでに明確に認識されており、

- match cost
- rewrite loop
- canonical ordering

を人手で管理している。

## 12.9 branch-sensitiveなheuristic

`IntMatchIn`には、

```text
Log[a b] → Log[a] + Log[b]
Log[a^r] → r Log[a]
```

型の入力簡約が存在する。

複素数領域では一般にbranch-sensitiveであるため、この時代のintegratorが**coverageを広げるため比較的大胆なnormalizationを局所的に利用していた**ことが分かる。

---

# 13. Definite integration

`Calculus/DefiniteIntegrate.m`は951行で、冒頭から、

> still under development and has not been tested as fully as it might

と明記されている。

このsourceは初期Mathematicaの定積分器の構造をかなり直接見せる。

## 13.1 protocol

`PiecewiseIntegrate[f,{x,xmin,xmax}]` は、

1. `{ok, notok}`
2. `Fail`
3. それ以外（internal codeへ任せる）

という複数のreturn protocolを持つ。

つまりsymbolic layerとinternal integratorが**協調して仕事を分担**する。

## 13.2 基本pipeline

```text
Integrate[f,x] で不定積分
       ↓
不定積分に Integrate が残ったか
       ↓
symbolic endpointかnumeric endpointか
       ↓
pole / singularity探索
       ↓
区間を分割
       ↓
各区間でone-sided Limit
       ↓
合計
```

という構成。

## 13.3 singularity detection

`DefinitePoles`はexpression structureを調べ、

- polynomial
- rational
- Log
- Power
- その他のsubexpressions

を分岐してsingularity候補を探す。

多項式方程式を解く補助関数も独自に持つ。

## 13.4 negative cache

`DefiniteFailures`は空listから始まり、一度失敗した `{f,x,xmin,xmax}` を保存する。

次回同じ入力が来ると早期に`DefiniteFail`へ落とす。

これは**failure memoization / negative cache**である。

## 13.5 numeric comparisonへの依存

point sortingや区間判定の一部では、

```text
N[x] == N[y]
N[x] < N[y]
N[x] > N[y]
```

を使う。

またソース自身が、

```text
This is an O(n^2) sort.
```

と明記するsort routineを持つ。

この部分は非常に実務的で、理想化されたsymbolic orderingではなく、**singularity候補が少数であることを前提にした小規模heuristic code**である。

## 13.6 endpoint fallbackの危険を自覚

symbolic endpointsでは、

```text
Limit[int, x->xmax] - Limit[int, x->xmin]
```

へ戻る場合があるが、そのとき

> singularities may be missed

というwarningを出す。

つまりシステム自身がfallbackの数学的限界を認識し、message systemを通じて露出させる。

---

# 14. D / Limit / Series / Residue

## 14.1 Dはkernel builtin

`oD`がMath A/B双方に見え、Series sourceでも普通に`D`を利用する。

## 14.2 `SeriesData` はfirst-class symbolic object

`info.m`の定義は、

```text
SeriesData[x, x0, {a0,a1,...}, nmin, nmax, den]
```

である。

ここでpowerは、

```text
nmin/den, (nmin+1)/den, ..., nmax/den
```

となる。

つまり単なるTaylor polynomialではなく、**fractional powersを表現できるseries container**になっている。構造上はPuiseux型展開も扱える。

## 14.3 `Series.m` の規模

`StartUp/Series.m`は1,130行あり、

- `/:` 95回
- `:=` 138回
- `/;` 99回

を含む。

33種類以上のheadに対してUpValue形式でseries knowledgeを追加している。

対象例:

```text
AiryAi
BesselI/J/K/Y
Beta / BetaRegularized
ChebyshevT/U
ExpIntegralE/Ei
GammaRegularized
GegenbauerC
HermiteH
Hypergeometric...
LaguerreL
LegendreP/Q
LerchPhi
PolyGamma
PolyLog
SphericalHarmonicY
Zeta
SinIntegral
CosIntegral
```

## 14.4 「函数がSeriesを知る」設計

典型的には、

```text
SpecialFunction /:
    Series[SpecialFunction[g],{...}] := ...

SpecialFunction /:
    SpecialFunction[s_SeriesData] := ...
```

というUpValue形式。

これは中央の`Series`に全特殊函数のcaseを書き込むのではなく、**各symbolが自分のseries protocolを持つ**設計である。

## 14.5 微分方程式による級数生成

Airy/Besselなどでは`diffeq...`系内部helperが使われる。

つまり特殊函数の展開係数をformula tableだけで列挙するのではなく、**満たす微分方程式からseries recurrenceを構成する**経路がある。

## 14.6 InverseSeries

`InverseSeries`はNewton iterationを利用し、

```text
i = 2
while i < n:
    i = 2*i
    ...
```

と精度次数を倍増させる。

形式冪級数の逆函数計算に、**Newton doubling**をすでに利用している。

---

# 15. Special Functions

特殊函数はkernelとsymbolic layerへ分散している。

## Math Aで確認できるもの

```text
Gamma
Beta
Zeta
PolyGamma
PolyLog
LerchPhi
Erf
```

等。

## Math B

```text
Hypergeometric0F1
Hypergeometric1F1
Hypergeometric2F1
HypergeometricU
LegendreP/Q
Gegenbauer
Hermite
Laguerre
EllipticE/K
SphericalHarmonicY
ArithmeticGeometricMean
```

等。

## Math C

```text
AiryAi
BesselI/J/K/Y
```

。

さらに`StartUp/Series.m`, `StartUp/Elliptic.m`, `InverseFunctions.m`, `IntegralTables.m`が、

- symbolic identities
- series
- inverse relations
- integral relations
- branch-specific expressions

を付加する。

したがって特殊函数systemは、

```text
numeric/kernel implementation
        +
symbolic knowledge layer
```

という二層構成と見るのが妥当である。

---

# 16. Linear Algebra

kernel側で確認できる主なものは、

```text
Det
Inverse
LinearSolve
Eigensystem
SingularValues
PseudoInverse
RowReduce
NullSpace
MatrixPower
MatrixExp
LatticeReduce
Dot
Inner
Outer
Transpose
Diagonal
Minors
IdentityMatrix
```

。

`LinearAlgebra` packageはむしろ小さく、

- `Cross.m`
- `Vectors.m`

のみ。

これは基本的な行列演算、高度な固有値/SVD等を**kernelへかなり深く入れている**ことを示す。

`ZeroTest` optionのusageも存在し、symbolic elementを含むlinear algebraで「何を0とみなすか」を函数として差し替えられる設計になっている。

---

# 17. Optimization / Constrained computation

## 17.1 `ConstrainedMin / Max`

Math Aには、

```text
oConstrainedMin
oConstrainedMax
```

がある。

## 17.2 `LinearProgramming` はwrapper

`StartUp/LinearProgramming.m`は47行程度で、

```text
LinearProgramming
   ↓
ConstrainedMin
```

へ変換している。

これは特殊問題のpublic APIを、より一般的なoptimization primitiveへreductionする設計。

なお`AreWeValid`内に、

```text
length1 = Length[vector1]
...
If[Length1 != length4, ...]
```

という **`length1` / `Length1` の不一致**があり、静的には明らかなtypoに見える。

実際のdump済み定義で修正されている可能性までは否定できないが、同梱sourceそのものにはこの粗さが残る。

---

# 18. Numerical algorithms package

`NumericalMath/Approximations.m`は987行で、次を実装する。

- `Pade`
- `EconomizedRationalApproximation`
- `RationalInterpolation`
- `MiniMaxApproximation`
- `GeneralRationalInterpolation`
- `GeneralMiniMaxApproximation`

内部では、

- Chebyshev polynomial
- LinearSolve
- extrema location
- iterative refinement
- WorkingPrecision
- convergence failure

を組み合わせる。

これは「数値算法は全部native codeで高速化」という構成ではなく、**高級数値算法そのものをMathematica languageで組み立てる**思想を示す。

---

# 19. Runge–Kutta packageに見る性能意識

`NumericalMath/RungeKutta.m`はDormand–Prince 4(5)系のadaptive RKを実装する。

特にソースコメントが面白い。

ODE式を毎step、

```text
f /. Thread[vars -> vals]
```

で置換する代わりに、一度、

```text
Function[Release[vars], Release[f]]
```

へ変換し、`Apply`で評価する。

コメントでは明示的に、

> This saves time

としている。

つまり作者は、

- symbolic replacementは便利だが高コスト
- 評価回数が多い数値loopではfunction化する

という**評価器コストモデル**を意識してコードを書いている。

---

# 20. DataAnalysis — distributionを「式object」として実装

この領域は初期Mathematicaのexpression modelを理解する非常に良い例である。

`ContinuousDistributions.m`は、

```text
BetaDistribution[alpha,beta]
NormalDistribution[mu,sigma]
StudentTDistribution[n]
...
```

を単なるtagではなく**first-class symbolic object**として扱う。

そしてUpValueで、

```text
BetaDistribution /: Density[BetaDistribution[a,b], x] := ...
BetaDistribution /: Mean[BetaDistribution[a,b]] := ...
BetaDistribution /: Variance[...] := ...
BetaDistribution /: Quantile[...] := ...
BetaDistribution /: Random[...] := ...
```

というprotocolを持たせる。

これは現代的な用語なら、

```text
object = expression head + parameters
methods = tagged rewrite rules / upvalues
```

に近い。

継承class hierarchyを導入せず、**symbolic dispatchだけでobject protocolを作っている**点が重要。

## 20.1 Continuous distributions

17種:

- Beta
- Cauchy
- Chi
- ChiSquare
- Exponential
- ExtremeValue
- F
- Gamma
- Normal
- HalfNormal
- Laplace
- LogNormal
- Logistic
- Rayleigh
- StudentT
- Uniform
- Weibull

## 20.2 Discrete distributions

10種:

- Bernoulli
- BetaBinomial
- BetaPascal
- Binomial
- DiscreteUniform
- DiscreteWeibull
- Geometric
- Hypergeometric
- NegativeBinomial
- Poisson

## 20.3 protocol

`DescriptiveFunctions.m`側が、

```text
Density
CumulativeDensity
Mean
Variance
StandardDeviation
Skewness
Kurtosis
CharacteristicFunction
Quantile
```

という共通message/protocolを定義する。

これを各distribution headがUpValueで実装する。

これは非常に再利用性の高い設計である。

---

# 21. Discrete mathematics / number theory

package layerには、

- Gosper summation
- combinatorial functions
- permutation/cycle変換
- Clebsch–Gordan / Wigner / Racah
- tree expression tools
- continued fractions
- integer root processing
- number recognition

がある。

一方kernelにも、

```text
Prime
PrimeQ
FactorInteger
PowerMod
ExtendedGcd
JacobiSymbol
MoebiusMu
DivisorSigma
BernoulliB
StirlingS1/S2
PartitionsP/Q
Binomial
Multinomial
```

等がある。

ここでも**低レベル/頻用primitiveをkernel、特殊algorithmをpackage**という境界が見える。

---

# 22. Graphics architecture

## 22.1 kernel側

`Math A/B`には、

```text
Plot
Plot3
ParametricPlot
ContourPlot
DensityPlot
ListPlot
Graphics3D
SurfaceGraphics
Show
```

等がある。

`info.m`冒頭には`Graphics`, `Graphics3D`, primitives, optionsのusageが大量に並ぶ。

## 22.2 package側

`Graphics` directoryだけで1,577行あり、

- Animation
- color spaces
- log/polar plots
- bar/pie chart
- multiple plot
- ParametricPlot3D
- Polyhedra
- Shapes
- ThreeScript

を拡張する。

## 22.3 Front End側

一方Front EndのCODEには、

```text
Post / Post2 / Post3
ConvPS
Image
Printing
ViewPt
```

がある。

従って強い推定として、

```text
Kernel:
    symbolic Graphics expression生成
        ↓
Front End:
    display / PostScript / printing / viewpoint UI
```

という責務分割である。

---

# 23. Notebook / Cell model

Front End string/resourceには、

- Cell
- Notebook
- evaluation group
- cell group
- init/noninit cell
- cell form/type/height
- page break
- merge/ungroup
- style change
- undo/redo

の操作が大量に存在する。

Undo/Redo文字列だけでも、

```text
Undo Cut Cells
Undo Evaluation Group
Undo Group Open/Close
Undo Init/Noninit
Undo Cell Form Change
...
```

がある。

したがってNotebookは単なる「テキスト欄 + output欄」ではなく、

```text
document
  └─ cell groups
      ├─ input
      ├─ output
      ├─ message
      ├─ text
      └─ initialization metadata
```

を持つ構造化document editorとして実装されている。

---

# 24. Front End と kernel のプロトコル

`sysinit.m`はFront End contextに、

```text
FE`FC
FE`FlC
FE`FT
```

を定義する。

`FC`はcompletion、`FT`はfunction template生成用。

これらは`Names`, `ToExpression`, `::usage`, `Out`等の通常のkernel機能を使ってFront Endへ情報を返す。

つまりFront End専用機能の一部ですら、

```text
GUI request
   ↓
special FE` command
   ↓
ordinary symbolic kernel evaluation
   ↓
stdout/protocol
```

として構築されている。

これは「kernel APIを別物として大量に作る」より、**symbolic language自身をUI protocolにも再利用する**設計である。

---

# 25. Error / failure model

この版では失敗がひとつのmechanismに統一されているわけではない。

少なくとも、

- unevaluated expression
- `Fail`
- `Indeterminate`
- `Message`
- warning message
- `Check`
- `Throw/Catch`
- package固有sentinel (`DefiniteFail`等)

を使い分ける。

これは型付きexception systemではなく、**symbolic evaluatorに馴染む複数のfailure channel**である。

良い面は、数学的に未評価の式をそのまま値として保持できること。

悪い面は、packageごとにfailure protocolが異なりやすいことである。

---

# 26. 実装者が認識していた危険箇所

ソースコメントには、驚くほど率直に実装上の危険が残されている。

### rewrite performance

`IntegralTables.m`:

```text
Patterns with a head of Plus ... can greatly slow the integrator
```

### rewrite loop

```text
This rule can lead to an infinite loop
```

### definite integration

```text
still under development
singularities may be missed
```

### algorithmic complexity

```text
This is an O(n^2) sort.
```

### Series

```text
this is slow, but I don't know how else to do it
```

### minimax approximation

convergence failureを`Fail`で扱う多数のbranchを持つ。

これらは単なる「未完成さ」の証拠ではない。

むしろ、**汎用symbolic evaluator上でアルゴリズムを書くと、termination / match complexity / canonicalization / numerical fallbackが主要問題になる**ことを、開発者が非常に早い段階から経験していた証拠である。

---

# 27. 静的解析から推定できるEvaluatorの概念構造

バイナリには、

```text
eval
eval0
eval1
topeval
symeval
funceval
substeval
normeval
argeval
```

等がある。

pattern側には、

```text
gmatch_pattern
scanPattern
scanBlank
CheckCondition
```

。

dispatch/value側には、

```text
ValueCell
BuildDispatch
UseDispatch
UpdateDispatch
```

。

これらから、厳密なcall graphではないが概念的には、

```text
Input expression
      ↓
top-level evaluation
      ↓
head / attributes確認
      ↓
Hold policyに従ってargument評価
      ↓
symbol/value/rule lookup
      ↓
pattern matching
      ↓
condition test
      ↓
builtin handler または rule rhs
      ↓
必要なら再評価
      ↓
normal result / unevaluated result / message
```

という構造を想定するのが最も自然である。

**注意:** この順序の細部を確定するには68k逆アセンブルが必要。

---

# 28. この版で特に重要な「kernel vs language」境界

静的解析から見える最も本質的な設計は、機能ごとにnative codeかMathematicaかを固定するのではなく、次のように分けていること。

## kernelへ置くもの

- expression evaluation
- pattern matching
- attributes
- assignments/rules
- arbitrary precision number machinery
- polynomial primitives
- core Solve
- core Integrate
- core Series representation/operations
- linear algebra primitives
- plotting primitives
- many special functions
- I/O / control flow primitives

## Mathematica languageへ置くもの

- integral rule tables
- definite integration heuristics
- special-function series knowledge
- inverse function table
- GroebnerBasis public wrapper
- LinearProgramming wrapper
- statistics distributions
- confidence intervals
- numerical approximations
- Runge–Kutta
- graphics extensions
- units
- domain-specific discrete math

この境界の要点は、

> **「高速なprimitiveはkernel、数学知識とcompositionはlanguage」**

である。

ただし完全分離ではない。`Integrate`, `Series`, `Solve`のような大機能は明確に**hybrid**である。

---

# 29. `StartUp`を数学知識層として見る

`StartUp`は単なるstartup helperではない。

4,008行の中に、

- attributes helper
- formatting
- arbitrary digits helper
- elliptic identities
- Groebner wrapper
- **integral tables**
- inverse functions
- linear programming
- **series / special-function expansion**
- help usage
- messages

が入っている。

したがって`StartUp`は「起動時設定」というより、

**System contextへ注入される標準symbolic knowledge base**

に近い。

`sysinit.m`から直接loadされない`IntegralTables.m`や`Series.m`が存在する点については、build/dump時に取り込まれるか、visible `sysinit`外のbootstrapでloadされる可能性が高い。ここは実行traceまたはkernel dump解析なしには断定しない。

---

# 30. 現時点で未確定の内部

このアーカイブだけでも構造はかなり分かるが、次は逆アセンブル領域。

1. expression objectのメモリlayout
2. symbol table / hash table構造
3. `ValueCell`の具体構造
4. pattern matcherのbacktracking strategy
5. `Flat/Orderless/OneIdentity` matching algorithm
6. Dispatch tableの内部表現
7. bignum limb representation
8. arbitrary precision transcendental algorithms
9. polynomial factorizationの低レベルalgorithm
10. Solveの内部elimination strategy
11. Integrate builtin側と`IntegralTables`の正確なcall order
12. garbage collector / reference counting
13. Math A/B/C間のsegment call ABI
14. local/remote kernel protocol packet format
15. graphics expressionからFront End rendererへのwire format

これらは`CODE` resourceの68k disassemblyと、可能なら実行時traceで追うべき項目である。

---

# 31. 構造解析上の主要結論

Mathematica 1.2.2f33 Enhancedは、静的に見ても「巨大な数学函数集」ではない。

中心は次の五層である。

```text
1. Expression evaluator
2. Pattern / Rule / Attribute runtime
3. Native numerical + algebraic primitives
4. Mathematica-written symbolic knowledge
5. Notebook Front End / kernel communication
```

そして数学機能はこの上へ、

```text
hard-coded algorithm
+ rewrite rules
+ symbolic object protocol
+ package composition
```

として積み上げられる。

特に重要なのは、**Mathematica言語自体が製品の拡張言語ではなく、製品本体を実装する言語でもある**こと。

- trig simplifierがrule list
- integratorの大部分がrule list
- definite integralがMathematica source
- special-function seriesがUpValue
- distribution systemがUpValue
- numerical algorithmsがMathematica source
- Front End補完さえsymbolic kernel commandを使う

という事実が、そのことを直接示す。

---

# 32. 現代のCAS/数式処理系へ持ち帰る価値が高いもの

最後に、この構造から一般化して持ち帰る価値が高い設計を挙げる。

## 32.1 「数学函数」より先に、式を扱う一般runtimeを強くする

個々の`Sin`, `Integrate`, `Solve`を増やす前に、

- expression
- head
- attribute
- pattern
- rule
- hold
- condition
- traversal

を強くする。

これらが十分一般的なら、数学機能のかなりの部分を**データと規則として記述できる**。

## 32.2 可読ruleと実行用dispatchを分離する

`Dispatch`は非常に重要。

```text
human-readable rules
        ↓ compile
optimized matcher structure
        ↓
evaluation
```

とすれば、数学知識の保守性と速度を同時に取りやすい。

rule engineを導入するなら、最初からこの二段構造を考える価値がある。

## 32.3 rule tableはfamilyごとに分割する

`IntegralTables.m`の、

```text
IntBase
ExpQuad
ExpLinear
PureTrigInt
TrigInt
IntMatchInt
...
```

という構造は優秀。

巨大な一枚rule listより、

- preprocessing
- family-specific rules
- fallback
- postprocessing

へ分離した方が、性能・termination・debuggabilityを管理しやすい。

## 32.4 generic transformationでcoverageを稼ぐ

affine argument liftingのように、

```text
f[a x+b]
    ↓
f[x]用の既存知識を再利用
```

する。

個別caseを100個増やすより、**問題を既知のcanonical familyへ写像する変換**を1個作る方が強い。

## 32.5 数学objectは専用classでなく「式 + protocol」でも作れる

distribution implementationは好例。

```text
NormalDistribution[mu,sigma]
```

を単なるASTとして保持し、

```text
Mean[...]
Density[...]
Quantile[...]
```

へのruleをheadへ関連付ける。

この方式は、

- series
- distributions
- units
- transforms
- domains
- solution objects

などへ広く応用できる。

## 32.6 `SeriesData`のような中間表現をfirst-classにする

Taylor係数を毎函数が独自vectorで持つより、

```text
SeriesData[var, center, coefficients, min, max, denominator]
```

のような共通IRを持つ。

その上で各函数がprotocolを追加する。

これは、

- Series
- Limit
- Residue
- asymptotic expansion
- special functions

を同じ基盤に載せやすい。

## 32.7 数学知識と高速primitiveを分離する

全てをnative codeへ入れる必要はない。

kernel側は、

- exact primitives
- robust evaluator
- efficient matcher
- core data structures

を強くする。

高水準の数学知識はsource ruleとして置く。

この分離により、数学機能を**再コンパイルなしで読め、直せ、追加できる**。

## 32.8 failureをcacheする

`DefiniteFailures`は単純だが実用的。

「高価なsymbolic attemptが失敗した」という情報も計算結果である。

同じ式に同じ重い探索を繰り返させないnegative cacheは、CASでは非常に効く。

## 32.9 rewrite systemにはterminationとcostの設計を最初から入れる

1.2.2のsourceがすでに、

- infinite loop
- slow Plus pattern
- fixed point
- canonical ordering

に苦労している。

rule engineを作るなら、

- rule priority
- cost
- direction
- iteration cap
- visited-state
- canonical measure
- tracing

を後付けにしない方がよい。

## 32.10 繰り返し数値評価ではsymbolic substitutionを避ける

RungeKuttaの、

```text
式を毎回replace
```

ではなく、

```text
一度callable functionへ変換
```

する最適化は現在でも本質的。

symbolic representationとhot numerical loopの間に、**lowering / compilation / callableization**の境界を持つべきである。

## 32.11 documentationとmessageをcodeから分離する

`info.m` / `msg.m`方式は古いが発想は良い。

- usage
- diagnostics
- implementation

を分ければ、数学コード本体を読みやすく保てる。

## 32.12 Front Endをkernelの特権APIへ依存させすぎない

補完やtemplate取得を、

```text
Names
ToExpression
::usage
```

等の通常のsymbolic commandで実現している。

可能な限り、UIも**公開言語semanticsを再利用**すると、kernel/frontendの結合が弱くなる。

## 32.13 「全部を一つの層で解く」を避ける

この版の最も大きな教訓はこれである。

SolveもIntegrateもSeriesも、

```text
native coreだけ
```

でも、

```text
rulesだけ
```

でもない。

**native primitive + symbolic IR + rule knowledge + heuristic orchestration**

の組合せで作られている。

汎用CASでは、このhybrid構造が極めて合理的である。

---

# 付録A. package/source inventory

行数はCR改行を正規化後に計数。Public symbols欄は各ファイル内の`::usage` assignmentから先頭のみ抽出。

| file | lines | bytes | public symbols（抜粋） |
|---|---:|---:|---|
| `Algebra/CountRoots.m` | 54 | 1,546 | CountRoots |
| `Algebra/GosperSum.m` | 118 | 3,484 | GosperSum |
| `Algebra/ReIm.m` | 116 | 2,583 | — |
| `Algebra/Trigonometry.m` | 196 | 6,767 | TrigCanonical, TrigFactor, TrigReduce, TrigToComplex, ComplexToTrig |
| `Calculus/DefiniteIntegrate.m` | 951 | 31,764 | $PiecewiseIntegrate |
| `Calculus/InverseLaplace.m` | 25 | 631 | InverseLaplace |
| `Calculus/Laplace.m` | 46 | 1,049 | Laplace |
| `Calculus/ODE.m` | 5 | 104 | — |
| `Calculus/VectorAnalysis.m` | 183 | 4,805 | Coordinates, ScaleFactors, Cartesian, Cylindrical, Spherical, Parabolic ほか10 |
| `DataAnalysis/ConfidenceIntervals.m` | 582 | 23,701 | PopulationMean, PopulationMeanInterval, PopulationVariance, PopulationVarianceInterval, PopulationMeanDifference, PopulationMeanDifferenceInterval ほか20 |
| `DataAnalysis/ContinuousDistributions.m` | 798 | 31,739 | BetaDistribution, CauchyDistribution, ChiDistribution, ChiSquareDistribution, ExponentialDistribution, ExtremeValueDistribution ほか11 |
| `DataAnalysis/DataManipulation.m` | 270 | 9,256 | Column, ColumnTake, ColumnDrop, ColumnJoin, RowJoin, DropNonNumeric ほか14 |
| `DataAnalysis/DescriptiveFunctions.m` | 79 | 3,229 | Density, CumulativeDensity, Mean, Variance, StandardDeviation, Skewness ほか4 |
| `DataAnalysis/DescriptiveStatistics.m` | 222 | 5,865 | LocationReport, GeometricMean, HarmonicMean, RootMeanSquare, TrimmedMean, InterpolatedQuantile ほか14 |
| `DataAnalysis/DiscreteDistributions.m` | 597 | 23,618 | BernoulliDistribution, BetaBinomialDistribution, BetaPascalDistribution, BinomialDistribution, DiscreteUniformDistribution, DiscreteWeibullDistribution ほか4 |
| `DiscreteMath/ClebschGordan.m` | 113 | 2,726 | Clebsch, Wigner, Racah |
| `DiscreteMath/CombinatorialFunctions.m` | 59 | 1,284 | Subfactorial, CatalanNumber, Fibonacci, Hofstadter |
| `DiscreteMath/CombinatorialSimplification.m` | 35 | 827 | — |
| `DiscreteMath/Permutations.m` | 69 | 1,351 | PermutationQ, ToCycles, FromCycles, RandomPermutation |
| `DiscreteMath/Tree.m` | 134 | 3,562 | MakeTree, TreeFind, TreePlot, ExprPlot |
| `Examples/CellularAutomata.m` | 65 | 1,576 | UpdateCA, EvolveCA, ShowCA, NumberedRule, CenterSpot |
| `Examples/CollatzProblem.m` | 18 | 439 | Collatz |
| `Examples/CrystalStructure.m` | 76 | 3,661 | — |
| `Examples/EllipticCurves.m` | 112 | 6,924 | — |
| `Examples/Factor.m` | 75 | 1,516 | — |
| `Examples/FunctionalProgramming.m` | 19 | 355 | — |
| `Examples/ModularArithmetic.m` | 63 | 1,817 | — |
| `Examples/Mortgages.m` | 25 | 758 | — |
| `Examples/RingTheory.m` | 30 | 825 | — |
| `Examples/RungeKutta.m` | 28 | 705 | RungeKutta |
| `Geometry/Polytopes.m` | 182 | 4,090 | Vertices, Edges, Faces, Coordinates, Area, Inscribed ほか4 |
| `Geometry/Rotations.m` | 60 | 1,573 | RotationMatrix2D, Rotate2D, RotationMatrix3D, Rotate3D |
| `Graphics/Animation.m` | 201 | 7,906 | Animation, Animate, Movie, MoviePlot, MoviePlot3D, MovieDensityPlot ほか7 |
| `Graphics/Colors.m` | 137 | 3,887 | CMYColor, YIQColor, HSBColor, HLSColor |
| `Graphics/Graphics.m` | 501 | 14,297 | LinearScale, LogScale, UnitScale, PiScale, MultipleListPlot, TextListPlot ほか14 |
| `Graphics/ParametricPlot3D.m` | 164 | 4,758 | ParametricPlot3D, PointParametricPlot3D, SpaceCurve, PointSpaceCurve, SphericalPlot3D |
| `Graphics/Polyhedra.m` | 174 | 6,481 | Polyhedron, Vertices, Faces, Polyhedra, Icosahedron, Dodecahedron ほか4 |
| `Graphics/Shapes.m` | 204 | 6,736 | Shapes, Cylinder, Cone, Torus, Sphere, MoebiusStrip ほか7 |
| `Graphics/ThreeScript.m` | 196 | 4,599 | ThreeScript |
| `LinearAlgebra/Cross.m` | 18 | 351 | Cross |
| `LinearAlgebra/Vectors.m` | 54 | 1,054 | — |
| `Miscellaneous/PhysicalConstants.m` | 84 | 1,884 | — |
| `Miscellaneous/Units.m` | 485 | 11,139 | Convert, ConvertTemperature, SI, MKS, Meter, Kilogram ほか11 |
| `NumberTheory/ContinuedFractions.m` | 44 | 872 | ContinuedFraction, ContinuedFractionForm |
| `NumberTheory/IntegerRoots.m` | 28 | 759 | BreakRoots |
| `NumberTheory/Recognize.m` | 31 | 566 | Recognize |
| `NumericalMath/Approximations.m` | 987 | 36,711 | Pade, EconomizedRationalApproximation, RationalInterpolation, MiniMaxApproximation, GeneralRationalInterpolation, GeneralMiniMaxApproximation |
| `NumericalMath/InverseStatisticalFunctions.m` | 319 | 8,497 | Erfc, InverseErf, InverseErfc, InverseGammaRegularized, InverseBetaRegularized |
| `NumericalMath/ListIntegrate.m` | 78 | 2,443 | ListIntegrate |
| `NumericalMath/RungeKutta.m` | 220 | 9,097 | RungeKutta, PlotODESolution |
| `StartUp/Attributes.m` | 38 | 1,055 | ClearAttributes, SetAttributes |
| `StartUp/Digits.m` | 45 | 1,111 | Digits |
| `StartUp/Edit.m` | 389 | 7,352 | Edit, EditIn, EditDef, Recall |
| `StartUp/Elliptic.m` | 324 | 11,764 | EllipticF, JacobiSN, JacobiSD, JacobiSC, JacobiCS, JacobiCN ほか22 |
| `StartUp/Formats.m` | 25 | 637 | ComplexInfinity |
| `StartUp/GroebnerBasis.m` | 20 | 447 | GroebnerBasis |
| `StartUp/IntegralTables.m` | 798 | 27,713 | SinIntegral, CosIntegral |
| `StartUp/InverseFunctions.m` | 37 | 1,029 | — |
| `StartUp/LinearProgramming.m` | 47 | 1,308 | LinearProgramming |
| `StartUp/RunThrough.m` | 39 | 746 | RunThrough |
| `StartUp/Series.m` | 1,130 | 38,959 | InverseSeries |
| `StartUp/ValueQ.m` | 17 | 330 | ValueQ |
| `StartUp/info.m` | 654 | 88,361 | Graphics, Graphics3D, ContourGraphics, DensityGraphics, SurfaceGraphics, Show ほか642 |
| `StartUp/msg.m` | 445 | 30,768 | — |
| `Utilities/Record.m` | 9 | 102 | — |
| `Utilities/ShowTime.m` | 47 | 1,175 | ShowTime |
| `init.m` | 86 | 4,648 | — |
| `sysinit.m` | 210 | 7,357 | — |

---

# 付録B. resource inventory

## `Mathematica.bin`

- 471 resources
- 26 types
- 59 `CODE`
- `CODE`合計 512,142 bytes
- 代表CODE:
  - `Notebook`
  - `Cells`
  - `Eval`
  - `MathTalk`
  - `Telecom`
  - `tcp`
  - `Printing`
  - `Image`
  - `Help`
  - `ViewPt`
  - `Validate`

## `Math A.bin`

- 53 resources
- 51 `CODE`
- `CODE`合計 632,230 bytes

## `Math B.bin`

- 58 resources
- 56 `CODE`
- `CODE`合計 766,308 bytes

## `Math C.bin`

- 6 resources
- 4 `CODE`
- `CODE`合計 59,722 bytes

---

# 付録C. `o...`形式で抽出できたkernel operator識別子

`Math A/B/C`からASCII文字列として抽出し、`o[A-Z]...` に一致したものは**393 unique identifiers**。

これは「393 builtins」と同義ではない。debug/internal/operator wrapperを含む可能性がある。しかしkernel機能の分布を見る資料として有用である。


### Math A.bin — 160 identifiers

```text
oAbs, oAccuracy, oAlgebraicRules, oAlias, oApart, oAppend, oArcCos, oArcCosh, oArcCot, oArcCoth
oArcCsc, oArcCsch, oArcSec, oArcSech, oArcSin, oArcSinh, oArcTan, oArcTanh, oB, oBegin, oBeginPackage
oBeta, oBetaRegularized, oCancel, oCharacters, oChop, oCoefficient, oCoefficients, oCollect
oComposeSeries, oConstrainedMax, oConstrainedMin, oContext, oContourPlot, oCos, oCosh, oCot, oCoth
oCsc, oCsch, oCycloPoly, oD, oDebug, oDecompose, oDensityPlot, oDistribute, oDivide, oEliminate, oEnd
oEndAdd, oEndPackage, oErf, oEvenQ, oExp, oExpInt, oExponent, oFactor, oFactorTerms, oFirst, oFromASCII
oFullSolve, oGamma, oGammaRegularized, oGet, oHead, oHypInt, oImplies, oInput, oInsert, oIntegerQ
oInverseFunction, oJoin, oLast, oLeafCount, oLength, oLerchPhi, oLimit, oListPlot, oListPlot3, oLogInt
oLogicalExpand, oMainSolve, oMaxMemused, oMeminuse, oMemoryConstrained, oN, oNProduct, oNRoots, oNSum
oNeeds, oNormal, oNumDen, oNumerical, oOddQ, oOrderedQ, oP, oParametricPlot, oPermutations, oPlot
oPlot3, oPolyGamma, oPolyLog, oPolynomialDivision, oPolynomialGCD, oPolynomialQ, oPolynomialQuotient
oPolynomialRemainder, oPower, oPrecision, oPrepend, oPrime, oPrimeQ, oRandom, oRational, oRationalize
oRead, oRemove, oResidue, oRest, oResultant, oRoots, oRule, oSec, oSech, oSeedRandom, oSequenceLimit
oSeries, oSeriesCoefficient, oSeriesData, oSetAccuracy, oSetPrecision, oShare, oSign, oSignQ, oSimplify
oSin, oSinh, oSolve, oSolveAlways, oSort, oSplice, oSqrt, oStringJoin, oStringLength, oStringMatchQ, oT
oTake, oTan, oTanh, oTimes, oTiming, oToASCII, oToExpression, oToRules, oTogether, oUnAlias, oUnique
oVariables, oX, oZeta
```

### Math B.bin — 226 identifiers

```text
oAccumulate, oAddTo, oAnd, oApply, oArg, oArithmeticGeometricMean, oArray, oAtomQ, oAttributes, oB
oBernoulliB, oBinomial, oBlock, oCases, oCatch, oCeiling, oChebyshevT, oChebyshevU, oCheck, oClear
oClearAll, oClose, oComplement, oCompose, oCompoundExpression, oCondition, oConjugate, oConstruct
oContextToFilename, oCount, oD, oDFactorial, oDSolve, oDecrement, oDepth, oDet, oDiagonal, oDimensions
oDims, oDispatchTable, oDisplay, oDivideBy, oDivisorSigma, oDivisors, oDot, oEigensystem, oElapsedTime
oEllipticE, oEllipticExp, oEllipticK, oEllipticLog, oEncode, oEndPush, oEnvironment, oEqualQ, oEulerE
oEulerPhi, oExpand, oExpandAll, oExpandDenominator, oExpandNumerator, oExtendedGcd, oF, oFactorInteger
oFactorial, oFindMinimum, oFindRoot, oFit, oFixedPoint, oFlatten, oFloor, oFor, oFourier, oFreeofQ
oFullDepth, oGCD, oGegenbauer, oGoto, oGraphics3D, oGreater, oGreaterEqual, oHash, oHeadCompose
oHermite, oHold, oHypergeometric0F1, oHypergeometric0F1Regularized, oHypergeometric1F1
oHypergeometric1F1Regularized, oHypergeometric2F1, oHypergeometric2F1Regularized, oHypergeometricU
oIdentity, oIdentityMatrix, oInOut, oIncrement, oInequality, oInfo, oInner, oIntegrate
oInterpolatingPolynomial, oInterrupt, oIntersection, oInverse, oJacobiP, oJacobiSymbol, oL, oLaguerre
oLatticeReduce, oLcm, oLegendreP, oLegendreQ, oLess, oLessEqual, oLevel, oLinearSolve, oLog, oMap
oMapAll, oMapAt, oMatchQ, oMatrixExp, oMatrixPower, oMatrixQ, oMax, oMemberQ, oMessage, oMessageName
oMin, oMinors, oMod, oModDet, oMoebiusMu, oMultinomial, oNBernoulliB, oNIntegrate, oNameQ, oNames
oNest, oNestList, oNot, oNullSpace, oOff, oOpenTemporary, oOperate, oOptions, oOrder, oOuter, oOverflow
oPartition, oPartitionsP, oPartitionsQ, oPlotRange, oPochhammer, oPop, oPosition, oPowerMod
oPreDecrement, oPreIncrement, oPrint, oPrintForm, oProbablePrimeQ, oProduct, oPseudoInverse, oPush
oPut, oQuit, oQuotient, oRange, oReIm, oRelease, oReplaceAll, oResetMedium, oResultant2, oReverse
oRotate, oRound, oRowReduce, oRun, oSameQ, oSave, oScan, oSelect, oSet, oSetDelayed, oSetOptions, oShow
oSignature, oSingularValues, oSize, oSphericalHarmonicY, oStirlingS1, oStirlingS2, oSubtractFrom, oSum
oSurfaceGraphics, oSwitch, oTHn, oTable, oTagSet, oTagUnset, oThread, oThrough, oThrow, oTicks
oTimesBy, oToString, oTranspose, oTrigExpand, oTrueQ, oUnderflow, oUnequalQ, oUnion, oUnsameQ, oUnset
oUpSet, oUpdate, oValue, oVecTester, oVectorQ, oWhich, oWhile, oWrite, oWriteString, oXor, oZ
```

### Math C.bin — 9 identifiers

```text
oAiryAi, oBesselI, oBesselJ, oBesselK, oBesselY, oMinus, oNumberQ, oPlus, oSubtract
```


---

# 付録D. 解析に用いた特に重要なsource

今後さらに深掘りする場合の優先順。

1. `Packages/StartUp/IntegralTables.m`
   - rule-driven integrator本体の知識層。
2. `Packages/Calculus/DefiniteIntegrate.m`
   - symbolic/internal integrator境界、pole handling、failure cache。
3. `Packages/StartUp/Series.m`
   - `SeriesData`, special-function protocol, Newton series reversion。
4. `Packages/Algebra/Trigonometry.m`
   - `Dispatch`を使ったrule compilationの最も読みやすい実例。
5. `Packages/sysinit.m`
   - boot sequenceとFront End/kernel境界。
6. `Packages/StartUp/info.m`
   - public language semanticsの一覧。
7. `Packages/StartUp/msg.m`
   - failure semanticsと診断設計。
8. `Packages/DataAnalysis/ContinuousDistributions.m`
   - UpValueを使ったsymbolic object protocol。
9. `Packages/NumericalMath/Approximations.m`
   - 高水準数値algorithmを言語自身で実装する例。
10. `Packages/NumericalMath/RungeKutta.m`
   - symbolic evaluation costを意識した数値loop設計。

---

## 解析メモ

この報告書は、対象アーカイブの**static evidenceのみ**から作成した。  
特に68k machine code内部については、symbol stringsから責務はかなり推定できる一方、具体的algorithmを断定していない。

次段階で最も情報量が大きいのは、`Math A/B/C` の`CODE` resourceを個別抽出し、68k disassemblerへ掛けて、

```text
ValueCell
BuildDispatch / UseDispatch
gmatch_pattern
eval / topeval
NewBignum
oIntegrate
oMainSolve
oSeries
```

周辺のcall graphを復元することである。
