# mmCal v1.5.1 と取り巻く環境を鑑みたリファクタリング指針
## Refactoring mmCal v1.5.1 in the Context of Its Surrounding Ecosystem

**対象:** mmCalculator / mmCal v1.5.1  
**基準日:** 2026-08-12  
**文書種別:** アーキテクチャ評価・リファクタリング方針  
**Status:** Design proposal / not an implementation specification

---

# 日本語版(English version at the bottom)

## 0. 文書の目的

本書は、mmCal v1.5.1 時点の設計を、単独のコード品質だけではなく、以下の周辺事例を踏まえて評価し、今後の拡張に耐えるためのリファクタリング方針を定めるものである。

主に参照した系統は次の通り。

- **Eigenmath** — 小型・自己完結型 CAS を長期間維持した例
- **Qalculate! / libqalculate** — 電卓から汎用数学エンジンへ二十年以上成長した例
- **Mathematica 1.x** — Expression / Rule / Pattern を中核に据えた記号処理系の歴史的完成例
- **SymPy / Maxima / Giac** — 巨大化した CAS が specialized domain、assumption、simplifier、外部数値基盤をどう分離するかの参考例

本書の目的は「他製品へ似せる」ことではない。

目的は、mmCal が現在持つ

- **厳密計算優先（exact-first）**
- **定義域・主値・未定義点を壊さない**
- **不明な場合に嘘を返さない**
- **BigInt から誤差保証付き数値評価まで自己完結する**
- **CLI を中心に小さく理解可能な系として維持する**

という性格を守りつつ、今後の機能追加で設計が崩れないようにすることである。

---

# 1. 結論

mmCal v1.5.1 は、現時点で再度の全面再構築を必要としない。

v1.5.0 で行った Lexer / Parser / Lowerer / AST / Evaluator / Simplifier / Solver / CertifiedEvaluator / 数値モデルの再整理は有効であり、これを捨てる理由はない。

一方、v1.5.1 までに CAS としての能力が急速に増えたことで、いくつかの「小さな不具合」に見えるものが、実際には今後の拡張性を左右する**設計上の警告灯**になっている。

最優先で対処すべきものは次の五つである。

1. **Canonical Algebra Layer**
   - 数学的に同値な和・積・商・冪を共通の代数表現へ落とす。
2. **EvaluationContext + EvaluationBudget**
   - 深さだけではなく、rewrite 数、生成 node 数、整数 bit 数、候補数、working precision 等を統一管理する。
3. **AssumptionContext**
   - `x>0`, `element[x,Real]`, `x!=0` 等の知識を全アルゴリズムで共有する。
4. **Specialized Mathematical Views / IR**
   - generic AST だけで polynomial、rational function、linear system、series 等を処理しない。
5. **Declarative Rule Layer**
   - Simplify / Integrate 等の数学知識を C++ の `if` の集合だけにせず、条件・cost・verification policy 付きの規則として整理する。

これらは「新しい CAS を作り直す」変更ではない。

既存 AST の**横に追加する補助層**として導入するのが適切である。

---

# 2. v1.5.1 時点の位置付け

mmCal は現在、単純な関数電卓ではない。

最も近い分類は、

> **from-scratch exact-first symbolic calculator / compact CAS**

である。

主要な性格は以下。

### 2.1 自己完結した数値塔

概念的には次の層を持つ。

```text
BigUInt
  ↓
BigInt
  ↓
Rational
  ↓
exact Complex
  ↓
BigFloat
  ↓
RealInterval / ComplexInterval
  ↓
DecimalApproximation / ComplexDecimalApproximation
```

有限十進入力も原則 exact Rational として扱う。

```text
0.1
→ 1/10

0.1+0.2
→ 3/10

0.1+0.2==0.3
→ True
```

これは host language の `double` や binary floating point を入力意味論に持ち込まないという、mmCal の非常に強い特徴である。

### 2.2 誤差保証付き数値評価

`N[expr,n]` は単に working precision を n 桁へ設定して文字列化する機構ではない。

理想化すると、

```text
exact expression
  ↓
working precision を増加
  ↓
interval enclosure を計算
  ↓
要求した decimal rounding が一意に決まるか判定
  ↓
DecimalApproximation
```

という構造を持つ。

このため `accuracy`, `precision`, `rationalize` を、単なる表示 metadata ではなく真値との関係を持つ機能として構築できる。

### 2.3 記号計算

v1.5 系では概ね、

```text
Lexer
  ↓
Parser
  ↓
Lowerer
  ↓
AST
  ↓
Evaluator
  ├─ Simplifier
  ├─ Solver
  ├─ D
  ├─ integrate
  └─ CertifiedEvaluator
```

という分離が成立している。

### 2.4 CAS 機能

v1.5.1 時点では、

- exact arithmetic
- complex arithmetic
- arbitrary precision
- trigonometric / hyperbolic functions
- logarithm / exponential
- Gamma 系・Zeta 系等の特殊函数
- `D`
- `integrate`
- `solve`
- symbolic simplify / expand / factor / collect
- matrix 基礎機能
- statistics
- DFT / FFT / convolution
- special-function based integration の拡張

まで進み、「小型 CAS」の問題領域へ完全に入っている。

積分側でも、初等関数だけでなく、

- `Ei`
- `Si`, `Ci`
- `Shi`, `Chi`
- `fresnelc`, `fresnels`
- `polylog`
- incomplete gamma
- hypergeometric 系
- elliptic 系

へ接続する設計が進んでいる。

この段階では、単純に built-in を増やすより、**数学知識をどう整理して共有するか**が重要になる。

---

# 3. 周辺環境から得られる主要な教訓

## 3.1 Eigenmath — 「複雑にならない」ことを設計する

Eigenmath は 2026-08-12 の確認時点で GitHub 上に **3,420 commits** を持ちながら、現在も非常に小さな C ベース CAS として全体を追跡可能な構造を保っている。

中心は `struct atom` である。

現在の `defs.h` では、

- CONS
- kernel symbol
- user symbol
- Rational
- `double`
- string
- tensor

を一つの atom union で表現し、算術式を binary tree として保持することが明示されている。

また、

```text
STACKSIZE = 100000
BLOCKSIZE = 10000
MAXBLOCKS = 2000
MAXDIM = 24
```

のように、runtime の制約が非常に直接的である。

Evaluator も単純で、interrupt と evaluation depth を明示的に監視し、depth が 1000 を超えた場合は停止する。

Eigenmath の長所は、抽象化の巧妙さではない。

> **システム全体を一人の頭に収まる規模へ保つこと自体を設計原則にしている**

点にある。

### mmCal が学ぶべき点

- 中核の概念数を増やしすぎない。
- 「便利だから」という理由だけで一級 object を追加しない。
- runtime 全体を追えることを価値として扱う。
- 機能追加の前に「この機能は既存抽象化へ自然に載るか」を問う。

### 真似るべきでない点

mmCal は既に、

- arbitrary precision float
- certified interval
- exact/approximate separation
- richer solver semantics

を持つ。

したがって Eigenmath の `Rational + double + generic atom` 程度まで単純化することはできない。

Eigenmath は「天井を低く設定することで単純さを維持した」成功例であり、mmCal の最終形ではない。

---

## 3.2 Qalculate! — 「複雑になっても生き残る」ことを設計する

libqalculate は確認時点で GitHub 上に **2,028 commits** を持つ。

現在は、

- arbitrary precision rational / floating point
- complex
- infinity
- interval arithmetic
- uncertainty propagation
- symbolic simplification
- differentiation / integration
- equations / inequalities
- assumptions
- units
- physical constants
- matrices / vectors
- statistics
- CLI / GUI

まで持つ巨大な数学エンジンである。

現在の必須数値基盤は **GMP + MPFR**。

ソースは `Number`, `MathStructure`, `Calculator` を中心にしつつ、

```text
MathStructure-calculate.cc
MathStructure-decompose.cc
MathStructure-differentiate.cc
MathStructure-factor.cc
MathStructure-gcd.cc
MathStructure-integrate.cc
MathStructure-isolatex.cc
MathStructure-limit.cc
MathStructure-matrixvector.cc
MathStructure-polynomial.cc
...
```

のように機能別へ分割されている。

### 歴史上の重要点

#### 2004: core math code の書き直し

初期の段階で expression representation、simplification、calculation を含む core math code の rewrite が行われた。

これは、機能を増やし続ける前に、数学 kernel の意味論を整理した例である。

#### 2004–2006: library 化

CLI / GUI 内部の計算コードから、`libqalculate` という独立 math engine へ分離した。

この判断により、frontend の寿命と数学 kernel の寿命を分離できた。

#### 2017: CLN → GMP/MPFR

長年育った上位 CAS を保ったまま低層の数値 backend を交換した。

これは数値 backend と上位 symbolic 層の境界が重要であることを示す。

#### 2017: branch semantics の整理

`cbrt[-8]` と `(-8)^(1/3)` を分離し、

- real root function
- principal complex power

を別の意味論として扱うようになった。

#### 2017–2019: interval の目的分離

interval arithmetic は一度導入して終わりではなく、

- user-facing interval
- uncertainty propagation
- precision tracking

の目的差によって何度も改良された。

#### 近年: failure / termination が主要課題

近年の release history では、

- infinite loop
- segfault
- pathological equation expansion
- pole handling
- interval solution accuracy
- extremely large calculations
- assumption warnings

等の修正が継続している。

成熟 CAS では「函数を増やすこと」より、

> **止まるか、壊れないか、意味を誤らないか**

が支配的課題になる。

### mmCal が学ぶべき点

- backend boundary を作る。
- evaluation resource を統一管理する。
- branch semantics を函数ごとではなく全体原則として持つ。
- assumption を一級の共有機構にする。
- property / fuzz testing を早期に入れる。
- frontend と kernel を論理的に分離する。
- 「大きな中央 object」を作りすぎない。

### 真似るべきでない点

Qalculate! は二十年以上の実用機能が `MathStructure` 等へ集積した結果、中央 object が非常に大きくなっている。

現在の TODO にも "Sane and stable API" が残っている。

mmCal はこの規模になる前に、

- generic expression
- mathematical domain representation
- evaluation context
- numerical backend

の境界を明確にしておくべきである。

---

## 3.3 Mathematica 1.x — rule system の力と危険性

Mathematica 1.x の静的解析から得られる最大の教訓は、

> **数学知識を evaluator の C kernel だけへ埋め込まず、pattern / rule layer へ移せる**

ことである。

積分表、Series、inverse function table 等が、kernel の expression / rule runtime 上へ実装されている。

一方で古いソース自身に、

- pattern が integrator を大きく遅くする
- rule が infinite loop を起こし得る

という警告が残る。

したがって mmCal が rule layer を導入する場合、

**Mathematica型の unrestricted rewrite engine をそのまま導入すべきではない。**

必要なのは、

- domain predicate
- applicability
- cost
- termination guard
- verification policy

を持つ制御された rule table である。

---

# 4. mmCal v1.5.1 の設計上の負債

以下は「バグ一覧」ではない。

今後の規模拡大で問題化する可能性が高い、構造上の弱点である。

---

## 4.1 D1 — Canonical algebraic representation が不足している

既に症状が出ている。

数学的に同一の式が、

```text
(cos[x]+sin[x])/2*exp[x]
```

と

```text
(cos[x]+sin[x])exp[x]/2
```

のように異なる AST 構造として存在し得る。

同様に、

```text
a/b*c
a*c/b
a*(1/b)*c
(a*c)/b
c*a*b^(-1)
```

は数学的には同じ有理積を表す場合があるが、structural equality では別物になる。

### 影響

- `==` の証明能力低下
- simplify の規則増大
- integrate の候補検証失敗
- factor / collect の不安定化
- common subexpression の認識低下
- hash / memoization の効率低下
- pattern rule の重複
- formatter と internal form の癒着

### 評価

**最優先の設計課題。**

これは表示だけの問題ではない。

---

## 4.2 D2 — generic AST 一種類で数学 domain を処理し続ける危険

AST は構文・記号式の保存には適切である。

しかし次の算法では、generic tree は必ずしも最適な表現ではない。

- polynomial GCD
- polynomial factorization
- rational function
- partial fractions
- linear system
- matrix decomposition
- series
- algebraic number
- root isolation

例えば、

```text
x^100 + 2x + 1
```

を毎回 `Add(Power(Symbol,...),...)` として探索する必要はない。

一度

```text
Polynomial<Rational>
```

へ落とした方が、算法も計算量も明確になる。

### リスク

generic AST だけで進むと、

```cpp
if (isAdd(...))
if (isMul(...))
if (isPower(...))
if (looksLikePolynomial(...))
```

が各モジュールへ拡散する。

---

## 4.3 D3 — evaluation resource model が局所 guard に寄っている

CertifiedEvaluator には病的に深い AST を防ぐための depth guard がある。

これは必要だが、depth だけでは十分ではない。

危険な式は例えば、

```text
depth               = 20
rewrite steps       = 100000
generated terms     = 500000
integer bits        = 10000000
solver candidates   = 100000
working precision   = 1000000
```

のように浅くても爆発する。

### 必要なもの

共通の `EvaluationBudget`。

監視対象例:

- recursion depth
- AST node visits
- rewrite steps
- generated nodes
- generated terms
- BigInt bit length
- polynomial degree
- solver candidate count
- integral candidate count
- working precision
- interval refinement count

### 結果型

少なくとも、

```text
Completed
NotApplicable
Unresolved
BudgetExceeded
Cancelled
```

は区別すべきである。

---

## 4.4 D4 — exact-first semantics と execution strategy を分ける必要がある

mmCal の exact-first は守るべきである。

しかし、

> exact に意味を持つこと  
> と  
> 常に exact algorithm で実行すること

は同義ではない。

例:

- FFT
- 大規模 matrix
- plotting 用 sampling
- exploratory numeric solve
- massive statistics

では machine floating point が合理的な場合がある。

### 推奨分離

```text
Semantic layer
  Exact value / symbolic meaning

Execution policy
  ├─ Exact symbolic
  ├─ Certified arbitrary precision
  └─ Fast machine approximate
```

Machine evaluator を BigFloat の特殊 case として継ぎ足すべきではない。

独立した execution engine とする。

---

## 4.5 D5 — 自前数値 backend が上位層へ漏れる危険

BigInt / BigFloat / interval を自前で持つこと自体は欠点ではない。

むしろ mmCal の重要な価値である。

問題は、上位層が低層実装の詳細を直接知りすぎる場合である。

Qalculate! が長期運用中に CLN から GMP/MPFR へ移行できたことは、backend boundary の価値を示す。

mmCal でも外部ライブラリへ移行する必要はない。

むしろ、

> **自前 backend を長く維持するために boundary を作る**

べきである。

---

## 4.6 D6 — 数学知識が C++ recognizer へ散らばる危険

積分器が典型である。

機能が増えるほど、

```cpp
if (isSin(...))
if (isAffine(...))
if (isPower(...))
if (matchPolynomialTimesExp(...))
...
```

という recognizer が増える。

初期は読みやすいが、数百規則になると、

- 重複
- precedence conflict
- rule ordering
- infinite transformation
- coverage hole
- domain condition の漏れ

が起きる。

### 解決方向

完全な user-visible pattern language は不要。

まず、

```text
Rule
  pattern
  conditions
  transform
  cost
  verification policy
```

という C++ data structure へ整理するだけでよい。

---

## 4.7 D7 — assumption の知識共有を一級化する必要

現在既に、

```text
element[x,Real]
x>=0
x!=0
```

等の条件が simplify / solve 等へ影響する。

今後、

- simplify
- solve
- integrate
- limit
- series
- power
- log
- sqrt
- abs
- sign
- numerical evaluation

の全てが assumption を見る。

各モジュールが独自に、

```cpp
isKnownPositive(...)
isKnownReal(...)
isNonZero(...)
```

を実装すると矛盾する。

### 必要なモデル

三値論理を基本にする。

```text
True
False
Unknown
```

例:

```text
AssumptionContext

x:
  Real       = True
  Integer    = Unknown
  Positive   = True
  Zero       = False
  Finite     = True
  LowerBound = 0
  UpperBound = +Infinity
```

---

## 4.8 D8 — failure semantics を Solver 以外にも広げる必要

Solver は、

- Empty
- Finite
- Universal
- Conditional
- Unresolved

のように結果の性質を区別する方向へ進んでいる。

この考え方は他の subsystem へも広げるべきである。

例えば Integrate で、

```text
NoAntiderivativeKnown
UnsupportedFunction
DomainAmbiguous
BudgetExceeded
CandidateRejected
```

は異なる。

現状「未評価式を返す」だけでは、ユーザーにも開発者にも理由が分かりにくい。

内部 result object に reason を保持し、CLI では必要に応じて WARN として出すのがよい。

---

## 4.9 D9 — cache / memoization の前提となる stable canonical identity が弱い

CAS が大きくなると同じ subexpression を繰り返し処理する。

特に、

- `D`
- integrate candidate verification
- assumption queries
- polynomial conversion
- certified constants
- simplification

では cache が効く。

しかし数学的同値な式が別 AST なら cache hit 率が低い。

したがって cache 最適化より先に canonical algebra が必要である。

---

## 4.10 D10 — Kernel と CLI の論理分離を今のうちに維持する

UI を追加する必要はない。

しかし、

```text
Math Kernel
  ↑
CLI
```

という依存方向は守るべきである。

Qalculate! が早期に `libqalculate` を独立させたことは参考になる。

mmCal を共有ライブラリ化する必要はまだないが、

- CLI options
- history
- prompt
- formatting preferences

が mathematical evaluator へ侵入しないようにする。

`:fix` が presentation-only で exact Out 値を変更しない現在方針は、この意味で正しい。

---

## 4.11 D11 — property / fuzz testing が今後必須

機能数が増えると hand-written regression test だけでは不足する。

特に CAS では正しい関係式そのものを oracle として利用できる。

例:

```text
parse(format(parse(x))) ≡ parse(x)
simplify(x) ≡ x
expand(factor(p)) ≡ p
D(integrate(f,x),x) ≡ f
solve(f==0) の返した解を代入すると 0
certified interval contains reference value
rationalize(exact-derived approximation) recovers source when uniquely possible
```

random AST を大量生成して検査できる。

Qalculate! が大量の random expression testing を導入した歴史は、この段階の CAS に非常に参考になる。

---

# 5. リファクタリングの基本原則

## 5.1 やること

- 既存 AST は残す。
- 数値塔は残す。
- exact-first semantics は残す。
- CLI syntax は原則変更しない。
- specialized view を AST の横へ追加する。
- context / budget を shared infrastructure とする。
- rule は declarative に整理する。
- failure reason を内部で保持する。
- testability を architecture requirement とする。

## 5.2 やらないこと

- Mathematica 互換 evaluator への全面変更
- unrestricted pattern language
- UI / notebook 化
- GMP / MPFR への置換
- v1.5 の再度の全面 rewrite
- 全 built-in の class hierarchy 化
- 一つの巨大 `Expression` object へ全機能を押し込む
- 「速いから」という理由で exact semantics を machine float へ変える

---

# 6. 推奨ターゲットアーキテクチャ

```text
┌─────────────────────────────────────────────┐
│ CLI / Presentation                          │
│ prompt, history, :fix, formatter options   │
└──────────────────────┬──────────────────────┘
                       │
┌──────────────────────▼──────────────────────┐
│ Frontend                                    │
│ Lexer → Parser → Lowerer                    │
└──────────────────────┬──────────────────────┘
                       │
┌──────────────────────▼──────────────────────┐
│ Core Expression Layer                      │
│ AST / Symbol / Value                       │
└──────────┬───────────────────┬──────────────┘
           │                   │
           │                   │
┌──────────▼───────────┐ ┌────▼───────────────────┐
│ Canonical Algebra    │ │ Specialized Views / IR │
│ Add / Mul / Power    │ │ Polynomial             │
│ Rational monomials   │ │ RationalFunction       │
│ stable ordering      │ │ LinearSystem           │
└──────────┬───────────┘ │ Series                  │
           │             │ AlgebraicNumber        │
           │             └────┬───────────────────┘
           └──────────────┬────┘
                          │
┌─────────────────────────▼───────────────────┐
│ Evaluation Infrastructure                  │
│ EvaluationContext                          │
│ EvaluationBudget                           │
│ AssumptionContext                          │
│ Cancellation                               │
│ Diagnostics                                │
└─────────────────────────┬───────────────────┘
                          │
┌─────────────────────────▼───────────────────┐
│ Mathematical Algorithms                   │
│ Simplify / D / Integrate / Solve / Limit   │
│ Matrix / Polynomial / Special Functions    │
│ Declarative Rule Tables                    │
└─────────────────────────┬───────────────────┘
                          │
┌─────────────────────────▼───────────────────┐
│ Numerical Semantics                        │
│ Exact engine                               │
│ Certified arbitrary-precision engine       │
│ Machine approximate engine (future)        │
└─────────────────────────┬───────────────────┘
                          │
┌─────────────────────────▼───────────────────┐
│ Numeric Backend                            │
│ BigUInt / BigInt / Rational / BigFloat     │
│ RealInterval / ComplexInterval             │
└─────────────────────────────────────────────┘
```

重要なのは、上下関係を厳密に守ることである。

Specialized IR は AST を置き換えない。

**algorithm が必要なときだけ AST から view を構築する。**

---

# 7. Refactor R1 — Canonical Algebra Layer

最優先。

## 7.1 目的

以下を同一の数学表現として扱いやすくする。

```text
a*b
b*a

a/b*c
a*c/b

2*x/4
x/2

x*x*x
x^3
```

ただし branch-sensitive な Power まで無条件に代数化してはならない。

### 原則

- `Add` と `Mul` の associative / commutative な範囲を canonicalize。
- exact numeric coefficient を集約。
- factor exponent map を利用。
- denominator は負 exponent として扱えるが、definedness を失わない metadata を保持。
- noncommutative object を将来導入する余地を残す。
- branch-sensitive Power は conservative に扱う。

## 7.2 概念例

```cpp
struct CanonicalProduct {
    Rational coefficient;
    std::vector<FactorPower> factors;
};

struct FactorPower {
    Expr base;
    Rational exponent;
};
```

ただし実装では、

- integer exponent
- rational exponent
- generic exponent

を同列に扱うべきではない。

branch safety のため分類が必要。

## 7.3 `CanonicalAlgebraView`

既存 AST を破壊的に書き換えるより、

```cpp
auto view = CanonicalAlgebraView::tryCreate(expr, ctx);
```

のような非破壊 view が安全。

成功した範囲だけ canonical comparison に利用する。

## 7.4 stable ordering

和や積の順序は deterministic にする。

推奨 key:

1. exact numeric
2. symbol
3. power
4. function
5. compound expression

ただし表示順と内部順を完全に一致させる必要はない。

内部 canonical order と pretty printer は分ける。

---

# 8. Refactor R2 — EvaluationContext / EvaluationBudget

## 8.1 EvaluationContext

全 algorithm へばらばらの option を渡さない。

```cpp
struct EvaluationContext {
    AssumptionContext assumptions;
    EvaluationBudget* budget = nullptr;
    CancellationToken* cancellation = nullptr;
    EvaluationPolicy policy;
    Diagnostics* diagnostics = nullptr;
};
```

## 8.2 EvaluationPolicy

例:

```cpp
enum class ApproximationPolicy {
    ExactOnly,
    ExactPreferred,
    CertifiedApproximate,
    MachineApproximate
};
```

これは「値の意味」ではなく実行方針。

## 8.3 EvaluationBudget

```cpp
struct EvaluationBudget {
    std::size_t maxDepth;
    std::size_t maxNodeVisits;
    std::size_t maxRewriteSteps;
    std::size_t maxGeneratedNodes;
    std::size_t maxTerms;
    std::size_t maxCandidates;
    std::size_t maxIntegerBits;
    std::size_t maxWorkingPrecision;
};
```

全てを最初から hard limit にする必要はない。

まず counter を導入し、benchmark / fuzz で実態を観測してから default を決める。

## 8.4 budget はエラーではない

`BudgetExceeded` は数学的失敗ではない。

```text
能力的には解ける可能性がある
ただし現在の計算資源制限で停止した
```

という結果である。

これは「数学的にできない」と「実装能力がない」の区別にも使える。

---

# 9. Refactor R3 — AssumptionContext

## 9.1 共通 query API

```cpp
TruthValue isReal(const Expr&, const AssumptionContext&);
TruthValue isInteger(const Expr&, const AssumptionContext&);
TruthValue isPositive(const Expr&, const AssumptionContext&);
TruthValue isNonZero(const Expr&, const AssumptionContext&);
TruthValue isFinite(const Expr&, const AssumptionContext&);
```

```cpp
enum class TruthValue {
    False,
    Unknown,
    True
};
```

## 9.2 inference

例:

```text
x > 0
  ⇒ Real[x]
  ⇒ x != 0

x ∈ Integer
  ⇒ x ∈ Real

x > 3
  ⇒ x > 0
```

初期実装では完全な theorem prover は不要。

単純 closure だけでも大きい。

## 9.3 subsystem 間共有

同じ `AssumptionContext` を、

- Simplifier
- Solver
- Integrator
- D
- Limit
- CertifiedEvaluator

へ渡す。

`log[x^2] → 2log[x]` の可否と、`sqrt[x^2] → abs[x]` の可否を別々の独自判定で処理しない。

---

# 10. Refactor R4 — Specialized Mathematical Views

## 10.1 PolynomialView

```cpp
template<class Coeff>
struct Polynomial {
    Symbol variable;
    std::vector<Coeff> coefficients;
};
```

または sparse representation。

用途:

- degree
- GCD
- factor
- rational root
- derivative
- resultant
- solver

## 10.2 RationalFunctionView

```text
P(x) / Q(x)
```

を numerator / denominator polynomial として持つ。

用途:

- cancellation
- holes / poles
- partial fractions
- rational integration
- inequality solving

**重要:** cancellation 前の excluded point を metadata として保持する。

## 10.3 LinearSystemView

多変数一次式を、

```text
A x = b
```

へ落とす。

Solver が generic AST から毎回係数を拾わない。

## 10.4 SeriesData

将来 `Series`, `Limit`, asymptotic, special functions を実装するなら、Mathematica 1.x の `SeriesData` に相当する dedicated representation が有力。

例:

```cpp
struct SeriesData {
    Symbol variable;
    Expr center;
    Rational minExponent;
    Rational step;
    std::vector<Expr> coefficients;
    Rational order;
};
```

初期は integer exponent の Taylor series だけでもよい。

---

# 11. Refactor R5 — Declarative Rule Layer

## 11.1 目的

C++ の recognizer を消すことではない。

**数学規則の inventory と適用条件をデータ化すること**が目的。

## 11.2 Rule model

```cpp
struct RewriteRule {
    RuleId id;
    RuleDomain domain;
    Pattern pattern;
    ConditionFn condition;
    TransformFn transform;
    Cost cost;
    VerificationPolicy verification;
};
```

## 11.3 VerificationPolicy

例:

```text
None
Structural
SimplifyDifference
DifferentiateAndCompare
SubstituteAndCheck
IntervalCheck
```

積分規則なら、

```text
DifferentiateAndCompare
```

を既定にできる。

## 11.4 cost

rewrite には「正しいか」だけでなく、

- expression size
- tree depth
- expected simplification
- algorithm class

を cost として持たせる。

同じ式を行ったり来たりする rule を防ぐ。

## 11.5 rule phase

一つの巨大 rule pool にしない。

例:

```text
Normalize
Canonicalize
AlgebraicSimplify
TrigSimplify
IntegrationRecognition
IntegrationReduction
PostSimplify
```

phase 間で方向を決める。

---

# 12. Refactor R6 — Result / Failure Model

共通の algorithm result を導入する。

```cpp
enum class AlgorithmStatus {
    Success,
    NotApplicable,
    Unresolved,
    BudgetExceeded,
    Cancelled
};

template<class T>
struct AlgorithmResult {
    AlgorithmStatus status;
    std::optional<T> value;
    ReasonCode reason;
};
```

### Integrate の例

```text
Success
  primitive certified by D

NotApplicable
  this rule does not match

Unresolved
  no supported closed form found

BudgetExceeded
  candidate explosion

Cancelled
  user interrupt
```

この区別があれば、

> 「能力的にできなかったのか、数学的に閉じないのか分からない」

問題を大幅に改善できる。

数学的に「標準函数で表現不能」を一般に証明することは難しいため、

`NoKnownClosedForm` と `ProvenImpossible` を混同しないことも重要。

---

# 13. Refactor R7 — Numerical Backend Boundary

外部 dependency を導入する提案ではない。

## 13.1 目的

BigInt 等の representation detail を上位 CAS から隠す。

例えば、

```cpp
class Integer {
public:
    ...
};
```

の内部が limb vector であることを Solver が知る必要はない。

## 13.2 algorithm backend

BigUInt 内でも、

```text
multiply
  ├─ schoolbook
  └─ Karatsuba

divide
  ├─ classic
  └─ Burnikel-Ziegler
```

のように strategy が増えている。

この選択を arithmetic API の外へ漏らさない。

## 13.3 将来の利点

- Toom-Cook
- FFT/NTT multiplication
- Newton division
- faster sqrt
- alternate BigFloat kernel

を入れても symbolic 層を変更しなくてよい。

---

# 14. Refactor R8 — Machine Approximate Engine

これは後段でよい。

## 14.1 exact semantics を変更しない

```text
0.1
```

は今後も exact `1/10`。

machine engine は明示的な execution mode とする。

## 14.2 用途

- plot sampling
- large FFT
- large matrix
- statistical batch
- approximate initial guess
- numerical solver seed

## 14.3 CertifiedEvaluator との関係

```text
CertifiedEvaluator
  correctness / guaranteed rounding

MachineEvaluator
  throughput / exploratory computation
```

役割を混ぜない。

MachineEvaluator の結果を `DecimalApproximation` の certified result と同じ型で偽装しない。

---

# 15. テスト戦略

## 15.1 Regression

既存 700+ 系統のテストを維持。

## 15.2 Property test

### Parser / formatter

```text
format(parse(format(parse(x))))
==
format(parse(x))
```

### Simplifier

```text
simplify(x) ≡ x
```

ただし structural equality ではなく semantic check を使用。

### Algebra

```text
factor(expand(p)) ≡ p
```

### Calculus

```text
D(integrate(f,x),x) ≡ f
```

積分定数差を考慮。

### Solver

返された各有限解 `r` について、

```text
f(r) == 0
```

を確認。

### Certified numerics

interval が independent high-precision reference を包含することを検査。

## 15.3 Differential testing

開発用テストでは、

- Python `decimal`
- mpmath
- Boost.Multiprecision
- MPFR-based utility

等を**test oracle にのみ利用する**ことは検討価値がある。

本体 dependency に入れない。

## 15.4 fuzz generator

grammar-aware random AST generator を作る。

生成時に、

- depth
- node count
- function categories
- singularities
- exact / approximate
- complex branches

を制御可能にする。

---

# 16. 性能計測

リファクタリングは性能劣化を伴い得る。

したがって benchmark corpus を固定する。

カテゴリ:

```text
BigInt
Rational
BigFloat
Pi / exp / log
simplify
factor
D
integrate
solve
matrix
FFT
special functions
formatter
```

測るもの:

- wall time
- peak generated nodes
- rewrite count
- max BigInt bits
- max working precision
- allocation count

単に ms だけ測るより、**なぜ遅くなったか説明できる counter** を持つ。

---

# 17. リファクタリング順序

## Phase 0 — Baseline Freeze

- v1.5.1 の regression を固定。
- benchmark corpus を固定。
- current known failures を catalog 化。
- formatter round-trip corpus を固定。

**新機能追加を止める必要はないが、大規模機能追加より先に行う。**

## Phase 1 — Evaluation Infrastructure

導入:

- `EvaluationContext`
- `EvaluationBudget`
- `CancellationToken`
- `AlgorithmResult`

既存 algorithm を一つずつ対応。

この段階では動作を変えない。

## Phase 2 — Canonical Algebra

- Add normalization
- Mul normalization
- exact coefficient collection
- stable term ordering
- safe factor/exponent representation

最初の成功条件:

```text
(cos[x]+sin[x])/2*exp[x]
```

と数学的同値な標準形が consistent に比較できること。

## Phase 3 — AssumptionContext

- three-valued predicates
- basic implication
- bounds
- Simplifier / Solver から移行
- Integrate / CertifiedEvaluator へ展開

## Phase 4 — Mathematical Views

順序:

1. PolynomialView
2. RationalFunctionView
3. LinearSystemView
4. SeriesData

## Phase 5 — Declarative Rules

まず Integrate の一部だけを移す。

特に、

- affine lifting
- direct special-function primitives
- trig power reductions

等の明確な規則から始める。

全 Simplifier を一度に移行しない。

## Phase 6 — Property / Fuzz Infrastructure

random AST と property tests を CI へ。

BudgetExceeded / Cancelled の corner case も含める。

## Phase 7 — Machine Engine

必要性が十分高まってから追加。

---

# 18. 優先度表

| 項目 | 優先度 | 理由 |
|---|---:|---|
| Canonical Algebra | A+ | 既に equality / integrate 検証へ症状が出ている |
| EvaluationContext / Budget | A | CAS の停止性を共通管理する基盤 |
| AssumptionContext | A | domain / branch soundness の共有基盤 |
| PolynomialView | A- | solver / factor / integrate 全てへ効く |
| RationalFunctionView | A- | hole / pole を保った rational algorithms に必須 |
| AlgorithmResult / failure reason | A- | 「解けない理由」の区別に必要 |
| Declarative rule table | B+ | Integrate の成長前に入れたい |
| Fuzz / property tests | B+ | 今の規模から費用対効果が急増 |
| Numeric backend boundary | B | 長期保守性 |
| SeriesData | B | Series / Limit / special function 拡張時 |
| MachineEvaluator | B- | 性能用途。exact-first semantics より後 |
| Public stable library API | C | CLI 一本なら急ぐ必要なし |
| General pattern language | D | 複雑性が高く、現時点では不要 |
| GUI / Notebook | 非目標 | mmCal の価値と直接関係しない |

---

# 19. リスク

## 19.1 過剰抽象化

最も大きなリスク。

Eigenmath から学ぶべきなのは、抽象化を増やすこと自体に価値はないという点。

新 layer は、

> **現在複数箇所で同じ問題が実際に発生している場合だけ追加する**

こと。

## 19.2 canonicalization が意味を壊す

特に、

- branch cuts
- zero denominator
- `0^0`
- complex powers
- logarithm identities

で危険。

canonical algebra は「数学的に安全な範囲」を明示する。

## 19.3 rule engine が第二の evaluator になる

rule layer が独自 control flow を持ち始めると設計が二重化する。

Rule は、

- match
- condition
- transform
- verify

までに限定する。

一般 programming language へ育てない。

## 19.4 context object の肥大化

Qalculate!型の巨大 option object を避ける。

Context は小さな subsystem に分ける。

```text
EvaluationContext
  references
    AssumptionContext
    EvaluationBudget
    CancellationToken
    Diagnostics
```

とし、全設定値を一 struct に詰め込まない。

---

# 20. 成功条件

このリファクタリングは、新 function 数では評価しない。

成功条件は以下。

### Correctness

- exact semantics を維持。
- branch / domain regression なし。
- pole / hole regression なし。
- certified rounding regression なし。

### Architecture

- polynomial algorithm が generic AST traversal をほぼ行わない。
- assumption query が一箇所へ集約される。
- all major algorithms が budget を監視可能。
- rule applicability と verification が分離される。

### Reliability

- pathological input が crash / stack overflow せず、理由付きで停止。
- random expression test で infinite rewrite を検出可能。
- failure が `Unresolved` と `BudgetExceeded` で区別される。

### Maintainability

- 新しい積分規則の追加で既存 recognizer の if-chain を触る量が減る。
- 新しい polynomial algorithm が Parser / Formatter を知らない。
- 数値算法変更が Solver / Integrate へ伝播しない。

---

# 21. 最終方針

mmCal は Eigenmath のように「小さいまま」を最終目標にはできない。

既に数値塔・solver・integrator・certified numerics がその天井を越えている。

一方、Qalculate! のように中央の `MathStructure` 相当へすべてを集積すると、将来 API と意味論の整理が非常に重くなる。

したがって mmCal が取るべき道は中間である。

> **AST は小さく保つ。  
> 数学 domain は specialized view に分離する。  
> 数学知識は rule と algorithm に分離する。  
> exact semantics と execution policy を分離する。  
> 全ての計算へ共有 context / budget / assumptions を通す。**

mmCal の現在の最大の価値は、function catalog の大きさではない。

```text
exact-first
soundness-first
self-contained
```

という三点である。

リファクタリングの目的はこの三点を強化することであり、より巨大な CAS の外観を真似ることではない。

---

# 22. 持ち帰る設計原則

最後に、周辺環境から持ち帰るべきものを短く固定する。

### Eigenmath から

- 全体を理解可能な規模に保つ。
- runtime の制約を隠さない。
- 単純な design は大きな資産。
- 不要な一級概念を増やさない。

### Qalculate! から

- backend は交換可能な境界を持つ。
- assumption は後付け feature ではなく kernel infrastructure。
- interval の「数学 object」と「誤差追跡」は区別する。
- slow / infinite / crash の制御は機能追加と同等に重要。
- random expression testing は成熟 CAS で非常に有効。
- frontend と math kernel の寿命を分ける。

### Mathematica 1.x から

- 数学知識は hard-coded algorithm だけでなく rule layer に分離できる。
- generic affine lifting のような共通変換は個別公式より価値が高い。
- SeriesData のような専用 symbolic representation は強力。
- unrestricted rewrite は停止性と性能問題を必ず生む。

### mmCal が守るもの

- finite decimal の exact semantics
- principal branch と real-root function の区別
- definedness / hole / pole の保存
- candidate を証明してから採用する姿勢
- exact / certified / machine の意味を混ぜない
- 「不明」を正直に返す failure semantics

---

# 23. 参考情報

本書作成時に参照した主要一次情報:

- mmCal v1.5 系開発引き継ぎ資料および積分 catalog
- Eigenmath GitHub repository (`georgeweigt/eigenmath`)
  - repository history
  - `src/defs.h`
  - `src/eval.c`
- Qalculate! / libqalculate GitHub repository (`Qalculate/libqalculate`)
  - source tree
  - current TODO
  - project news / release history
- Mathematica 1.2.2f33 Enhanced 静的構造解析結果

Git commit 数や外部プロジェクトの現況は **2026-08-12 に確認した時点**の値であり、将来変化する。

---

---

# English Version

## 0. Purpose

This document evaluates the architecture of mmCal v1.5.1 not only in isolation, but in the context of several mature symbolic and calculator systems:

- **Eigenmath** — a long-lived example of a small, self-contained CAS;
- **Qalculate! / libqalculate** — a calculator that evolved into a large mathematical engine over more than two decades;
- **Mathematica 1.x** — an early, highly developed example of expression- and rule-oriented symbolic computation;
- **SymPy / Maxima / Giac** — useful reference points for specialized domains, assumptions, simplification, and numerical backends in larger CAS implementations.

The goal is not compatibility with any of them.

The goal is to preserve the core identity of mmCal:

- **exact-first semantics**;
- preservation of domains, principal branches, holes, and undefined points;
- refusal to return a false “complete” answer when the system cannot prove one;
- a self-contained stack from BigInt to certified numerical evaluation;
- a compact, CLI-oriented system that remains understandable;

while making the architecture robust enough for continued growth.

---

# 1. Executive conclusion

mmCal v1.5.1 does **not** require another full rewrite.

The v1.5 reorganization of the lexer, parser, lowerer, AST, evaluator, simplifier, solver, certified evaluator, and numeric model was directionally correct and should be preserved.

However, the rapid growth of CAS functionality has exposed several architectural warning signs. Some currently look like small simplification or formatting problems, but they are symptoms of issues that will become expensive as the system grows.

The five highest-priority refactorings are:

1. **Canonical Algebra Layer**
   - normalize mathematically equivalent additive and multiplicative forms into a stable algebraic representation;
2. **EvaluationContext + EvaluationBudget**
   - manage recursion, rewrite steps, generated nodes, integer growth, candidate growth, and working precision consistently;
3. **AssumptionContext**
   - share facts such as `x>0`, `element[x,Real]`, and `x!=0` across all algorithms;
4. **Specialized Mathematical Views / IR**
   - avoid using only the generic AST for polynomials, rational functions, linear systems, series, and related domains;
5. **Declarative Rule Layer**
   - represent mathematical rules with conditions, costs, and verification policies instead of accumulating only C++ recognizer chains.

These should be introduced **beside the existing AST**, not by replacing the v1.5 architecture.

---

# 2. Position of mmCal v1.5.1

The best short classification is:

> **from-scratch exact-first symbolic calculator / compact CAS**

## 2.1 Self-contained numeric tower

Conceptually:

```text
BigUInt
  ↓
BigInt
  ↓
Rational
  ↓
exact Complex
  ↓
BigFloat
  ↓
RealInterval / ComplexInterval
  ↓
DecimalApproximation / ComplexDecimalApproximation
```

Finite decimal literals are exact by default.

```text
0.1
→ 1/10

0.1+0.2
→ 3/10

0.1+0.2==0.3
→ True
```

This avoids importing the host language's binary floating-point semantics into the language itself.

## 2.2 Certified numerical evaluation

`N[expr,n]` is not merely “compute at n digits and print n digits.”

The intended architecture is closer to:

```text
exact expression
  ↓
increase working precision
  ↓
compute an interval enclosure
  ↓
prove that requested decimal rounding is unique
  ↓
DecimalApproximation
```

This gives `accuracy`, `precision`, and `rationalize` a meaningful connection to the true mathematical value.

## 2.3 Symbolic pipeline

The v1.5 design is approximately:

```text
Lexer
  ↓
Parser
  ↓
Lowerer
  ↓
AST
  ↓
Evaluator
  ├─ Simplifier
  ├─ Solver
  ├─ D
  ├─ integrate
  └─ CertifiedEvaluator
```

## 2.4 CAS breadth

By v1.5.1 the project includes substantial functionality in:

- exact arithmetic;
- complex arithmetic;
- arbitrary precision;
- elementary and transcendental functions;
- special functions;
- symbolic differentiation;
- symbolic integration;
- equation solving;
- simplification / expansion / factorization / collection;
- basic matrices;
- statistics;
- DFT / FFT / convolution;
- special-function-based antiderivatives.

Integration work is already expanding beyond elementary functions into:

- `Ei`;
- `Si`, `Ci`;
- `Shi`, `Chi`;
- `fresnelc`, `fresnels`;
- `polylog`;
- incomplete gamma;
- hypergeometric forms;
- elliptic forms.

At this point, organizing and sharing mathematical knowledge matters more than simply adding built-ins.

---

# 3. Lessons from the surrounding ecosystem

## 3.1 Eigenmath — design for remaining small

At the time of review on 2026-08-12, Eigenmath had **3,420 commits** on GitHub while retaining a remarkably small C-based architecture.

Its central representation is `struct atom`.

The current `defs.h` explicitly represents expressions as binary trees and provides atom variants for:

- cons cells;
- kernel symbols;
- user symbols;
- rational numbers;
- `double`;
- strings;
- tensors.

Runtime limits are direct and visible, including values such as:

```text
STACKSIZE = 100000
BLOCKSIZE = 10000
MAXBLOCKS = 2000
MAXDIM = 24
```

The evaluator directly checks for interrupts and excessive evaluation depth.

Eigenmath's primary architectural strength is therefore not elaborate abstraction.

It is:

> **keeping the whole system small enough to remain understandable.**

### Lessons for mmCal

- Limit the number of core concepts.
- Do not introduce a first-class abstraction merely because it is convenient.
- Treat whole-system comprehensibility as an architectural asset.
- Ask whether a new feature naturally fits existing abstractions before extending the kernel.

### What not to copy

mmCal already has requirements beyond Eigenmath's numerical model:

- arbitrary-precision floating point;
- certified intervals;
- richer exact/approximate semantics;
- richer solver result semantics.

The Eigenmath ceiling is therefore too low for mmCal.

---

## 3.2 Qalculate! — design for surviving complexity

At the time of review, libqalculate had **2,028 commits** on GitHub.

It now supports:

- arbitrary-precision rational and floating-point numbers;
- complex and infinite values;
- intervals;
- uncertainty propagation;
- symbolic simplification;
- differentiation and integration;
- equations and inequalities;
- assumptions;
- units and physical constants;
- matrices and vectors;
- statistics;
- multiple frontends.

Its required low-level numeric dependencies are currently **GMP and MPFR**.

The source tree reflects decades of growth, with files such as:

```text
MathStructure-calculate.cc
MathStructure-decompose.cc
MathStructure-differentiate.cc
MathStructure-factor.cc
MathStructure-gcd.cc
MathStructure-integrate.cc
MathStructure-isolatex.cc
MathStructure-limit.cc
MathStructure-matrixvector.cc
MathStructure-polynomial.cc
...
```

### Historically important events

#### 2004: core mathematics rewrite

The project rewrote core mathematical code, including expression representation, simplification, and calculation, relatively early in its lifetime.

This is a useful example of stabilizing semantics before uncontrolled feature growth.

#### 2004–2006: library separation

The mathematical engine became `libqalculate`, separated from CLI and GUI frontends.

This separated the lifetime of the mathematical kernel from the lifetime of any particular UI.

#### 2017: CLN → GMP/MPFR

The low-level numerical backend was replaced while preserving the upper symbolic system.

This demonstrates the value of a real backend boundary.

#### 2017: root and power semantics

The project explicitly distinguished:

- real root functions such as `cbrt(-8) = -2`;
- principal complex powers such as `(-8)^(1/3)`.

This is precisely the class of semantic distinction a mature CAS must make.

#### 2017–2019: interval semantics evolved

Intervals had to be refined over several releases because different roles emerged:

- user-visible intervals;
- uncertainty propagation;
- precision tracking.

An interval “type” alone was not enough.

#### Recent releases: failure and termination dominate

Recent release history repeatedly includes fixes for:

- infinite or nearly infinite loops;
- crashes;
- pathological equation transformations;
- poles in inequalities;
- interval solution accuracy;
- huge calculations;
- assumption warnings.

The mature-CAS lesson is clear:

> **termination, failure semantics, and semantic safety eventually become as important as mathematical breadth.**

### Lessons for mmCal

- Define backend boundaries.
- Make evaluation resources a shared concern.
- Treat branch semantics as a global design problem.
- Make assumptions kernel infrastructure.
- Introduce property and fuzz testing before the system becomes enormous.
- Keep frontend lifetime separate from math-kernel lifetime.
- Avoid concentrating everything into one giant central object.

---

## 3.3 Mathematica 1.x — the power and danger of rules

Static analysis of Mathematica 1.x demonstrates the architectural power of placing mathematical knowledge above a generic expression/rule runtime.

Integration tables, series knowledge, inverse-function tables, and related systems can live outside the lowest-level kernel.

The same historical sources also contain warnings about:

- expensive patterns;
- rules that can create infinite loops.

Therefore mmCal should not copy unrestricted Mathematica-style rewriting.

A mmCal rule layer should be controlled and include:

- domain predicates;
- applicability conditions;
- cost;
- termination controls;
- verification policies.

---

# 4. Architectural debt in mmCal v1.5.1

These items are not merely bugs. They are structural weaknesses likely to become more expensive with growth.

---

## 4.1 D1 — insufficient canonical algebraic representation

A visible symptom already exists.

Mathematically equivalent expressions such as:

```text
(cos[x]+sin[x])/2*exp[x]
```

and:

```text
(cos[x]+sin[x])exp[x]/2
```

can have different AST shapes.

Similarly:

```text
a/b*c
a*c/b
a*(1/b)*c
(a*c)/b
c*a*b^(-1)
```

may represent the same rational product while remaining structurally distinct.

### Consequences

- weaker semantic equality;
- more simplifier rules;
- failed antiderivative verification;
- unstable factoring and collection;
- weaker common-subexpression recognition;
- lower cache hit rates;
- duplicated pattern rules;
- coupling between formatter and internal form.

### Assessment

**Highest priority architectural issue.**

This is not merely a formatting problem.

---

## 4.2 D2 — risk of using one generic AST for every mathematical domain

The AST is appropriate for syntax and symbolic preservation.

It is not the best representation for every algorithm.

Domains that benefit from dedicated representations include:

- polynomials;
- rational functions;
- linear systems;
- matrices;
- series;
- algebraic numbers;
- root isolation.

For example:

```text
x^100 + 2x + 1
```

should not need to be repeatedly rediscovered as a polynomial by traversing generic `Add`, `Power`, and `Symbol` nodes.

Without specialized views, code tends to devolve into repeated checks such as:

```cpp
if (isAdd(...))
if (isMul(...))
if (isPower(...))
if (looksLikePolynomial(...))
```

throughout the system.

---

## 4.3 D3 — evaluation resources are still guarded locally

A recursion-depth guard is necessary, but insufficient.

A dangerous expression can be shallow and still create:

```text
rewrite steps       = 100000
generated terms     = 500000
integer bits        = 10000000
solver candidates   = 100000
working precision   = 1000000
```

A shared `EvaluationBudget` is needed.

Potential counters:

- recursion depth;
- AST node visits;
- rewrite steps;
- generated nodes;
- generated terms;
- integer bit length;
- polynomial degree;
- solver candidates;
- integration candidates;
- working precision;
- interval refinements.

Statuses should distinguish:

```text
Completed
NotApplicable
Unresolved
BudgetExceeded
Cancelled
```

---

## 4.4 D4 — exact-first semantics must be separated from execution strategy

Exact-first semantics should remain.

However:

> being exactly representable  
> is not the same as  
> always requiring an exact execution algorithm.

Large FFTs, large matrices, plotting samples, exploratory numerical solves, and bulk statistics may reasonably use machine floating point.

Recommended separation:

```text
Semantic layer
  Exact mathematical meaning

Execution policy
  ├─ Exact symbolic
  ├─ Certified arbitrary precision
  └─ Fast machine approximate
```

A future machine evaluator should not be implemented as a special case hidden inside BigFloat.

---

## 4.5 D5 — risk of leaking self-written numeric backends upward

Using self-written BigInt, BigFloat, and interval arithmetic is not itself a flaw.

It is an important feature of mmCal.

The risk appears if upper symbolic code depends on backend representation details.

The correct goal is:

> **create backend boundaries in order to preserve the self-written backend for the long term.**

Upper algorithms should not know limb layout, multiplication thresholds, or division implementation details.

---

## 4.6 D6 — mathematical knowledge may scatter into C++ recognizer chains

Integration is the clearest example.

As coverage grows, code naturally accumulates:

```cpp
if (isSin(...))
if (isAffine(...))
if (isPower(...))
if (matchPolynomialTimesExp(...))
...
```

At hundreds of rules this causes:

- duplication;
- rule ordering conflicts;
- transformation loops;
- coverage holes;
- forgotten domain conditions.

The answer is not necessarily a general pattern language.

A controlled C++ rule data model is enough:

```text
Rule
  pattern
  conditions
  transform
  cost
  verification policy
```

---

## 4.7 D7 — assumptions need to become first-class shared infrastructure

The system already uses facts such as:

```text
element[x,Real]
x>=0
x!=0
```

As functionality grows, assumptions affect:

- simplify;
- solve;
- integrate;
- limit;
- series;
- power;
- logarithms;
- roots;
- absolute values;
- signs;
- numerical evaluation.

Each subsystem must not invent its own definition of `isKnownPositive()`.

Use a shared three-valued logic:

```text
True
False
Unknown
```

Example conceptual state:

```text
x:
  Real       = True
  Integer    = Unknown
  Positive   = True
  Zero       = False
  Finite     = True
  LowerBound = 0
  UpperBound = +Infinity
```

---

## 4.8 D8 — solver-style failure semantics should spread to other subsystems

The solver already moves toward rich result categories.

The same idea should be generalized.

Integration failures such as:

```text
NoKnownAntiderivative
UnsupportedFunction
DomainAmbiguous
BudgetExceeded
CandidateRejected
```

are not equivalent.

Internally preserving a reason code lets the CLI explain failures when useful.

The system must also distinguish:

- “no known closed form in the implemented function set”;
- “not implemented”;
- “resource limit reached”;
- “proven impossible” — which is much stronger and usually unavailable.

---

## 4.9 D9 — caching depends on stable canonical identity

As the system grows, repeated work becomes common in:

- differentiation;
- antiderivative verification;
- assumption queries;
- polynomial conversion;
- certified constants;
- simplification.

Caching before canonical identity is stable produces poor hit rates and potentially inconsistent keys.

Canonical algebra should therefore precede aggressive memoization.

---

## 4.10 D10 — preserve the logical kernel/frontend boundary

There is no need to add a GUI.

However, the dependency direction should remain:

```text
Math Kernel
  ↑
CLI
```

Presentation settings such as `:fix`, prompt state, history, and display preferences must not alter mathematical values.

The current rule that `:fix` is presentation-only is architecturally sound.

---

## 4.11 D11 — property and fuzz testing are becoming mandatory

Handwritten regression tests alone will not scale.

CAS implementations have unusually strong mathematical properties that can serve as oracles.

Examples:

```text
parse(format(parse(x))) ≡ parse(x)
simplify(x) ≡ x
expand(factor(p)) ≡ p
D(integrate(f,x),x) ≡ f
returned finite solve roots satisfy the original equation
certified intervals contain independent reference values
```

Grammar-aware random AST generation should be introduced before the codebase becomes substantially larger.

---

# 5. Refactoring principles

## 5.1 Do

- Preserve the existing AST.
- Preserve the numeric tower.
- Preserve exact-first semantics.
- Preserve CLI syntax unless a change is mathematically necessary.
- Add specialized views beside the AST.
- Make context and budget shared infrastructure.
- Organize rules declaratively.
- Preserve internal failure reasons.
- Treat testability as an architectural requirement.

## 5.2 Do not

- Rewrite the evaluator as a Mathematica clone.
- Add an unrestricted pattern language.
- Turn the project into a notebook/GUI system.
- Replace the numeric stack with GMP/MPFR.
- repeat the v1.5 full rewrite.
- create a giant class hierarchy for every built-in.
- move all behavior into one enormous `Expression` object.
- weaken exact semantics merely for speed.

---

# 6. Recommended target architecture

```text
┌─────────────────────────────────────────────┐
│ CLI / Presentation                          │
│ prompt, history, :fix, formatter options   │
└──────────────────────┬──────────────────────┘
                       │
┌──────────────────────▼──────────────────────┐
│ Frontend                                    │
│ Lexer → Parser → Lowerer                    │
└──────────────────────┬──────────────────────┘
                       │
┌──────────────────────▼──────────────────────┐
│ Core Expression Layer                      │
│ AST / Symbol / Value                       │
└──────────┬───────────────────┬──────────────┘
           │                   │
┌──────────▼───────────┐ ┌────▼───────────────────┐
│ Canonical Algebra    │ │ Specialized Views / IR │
│ Add / Mul / Power    │ │ Polynomial             │
│ Rational monomials   │ │ RationalFunction       │
│ stable ordering      │ │ LinearSystem           │
└──────────┬───────────┘ │ Series                  │
           │             │ AlgebraicNumber        │
           │             └────┬───────────────────┘
           └──────────────┬────┘
                          │
┌─────────────────────────▼───────────────────┐
│ Evaluation Infrastructure                  │
│ EvaluationContext                          │
│ EvaluationBudget                           │
│ AssumptionContext                          │
│ Cancellation                               │
│ Diagnostics                                │
└─────────────────────────┬───────────────────┘
                          │
┌─────────────────────────▼───────────────────┐
│ Mathematical Algorithms                   │
│ Simplify / D / Integrate / Solve / Limit   │
│ Matrix / Polynomial / Special Functions    │
│ Declarative Rule Tables                    │
└─────────────────────────┬───────────────────┘
                          │
┌─────────────────────────▼───────────────────┐
│ Numerical Semantics                        │
│ Exact engine                               │
│ Certified arbitrary-precision engine       │
│ Machine approximate engine (future)        │
└─────────────────────────┬───────────────────┘
                          │
┌─────────────────────────▼───────────────────┐
│ Numeric Backend                            │
│ BigUInt / BigInt / Rational / BigFloat     │
│ RealInterval / ComplexInterval             │
└─────────────────────────────────────────────┘
```

The essential rule is dependency direction.

Specialized IR does not replace the AST.

An algorithm constructs a specialized view only when useful.

---

# 7. R1 — Canonical Algebra Layer

Highest priority.

## 7.1 Goal

Make equivalent forms easier to treat uniformly:

```text
a*b
b*a

a/b*c
a*c/b

2*x/4
x/2

x*x*x
x^3
```

Branch-sensitive powers must remain conservative.

## 7.2 Principles

- Canonicalize associative/commutative portions of `Add` and `Mul`.
- Collect exact numerical coefficients.
- Use a stable factor/exponent representation.
- Preserve excluded-point information when denominators are transformed.
- Leave room for future noncommutative objects.
- Treat rational and generic exponents differently when branch semantics require it.

## 7.3 Non-destructive view

Prefer:

```cpp
auto view = CanonicalAlgebraView::tryCreate(expr, ctx);
```

to immediately replacing every AST node.

This makes migration incremental.

---

# 8. R2 — EvaluationContext / EvaluationBudget

## 8.1 Context

```cpp
struct EvaluationContext {
    AssumptionContext assumptions;
    EvaluationBudget* budget = nullptr;
    CancellationToken* cancellation = nullptr;
    EvaluationPolicy policy;
    Diagnostics* diagnostics = nullptr;
};
```

## 8.2 Policy

```cpp
enum class ApproximationPolicy {
    ExactOnly,
    ExactPreferred,
    CertifiedApproximate,
    MachineApproximate
};
```

This describes execution policy, not the mathematical identity of the value.

## 8.3 Budget

```cpp
struct EvaluationBudget {
    std::size_t maxDepth;
    std::size_t maxNodeVisits;
    std::size_t maxRewriteSteps;
    std::size_t maxGeneratedNodes;
    std::size_t maxTerms;
    std::size_t maxCandidates;
    std::size_t maxIntegerBits;
    std::size_t maxWorkingPrecision;
};
```

Initially, counters can be observed before strict defaults are imposed.

`BudgetExceeded` is not a domain error. It is a resource outcome.

---

# 9. R3 — AssumptionContext

Shared API:

```cpp
TruthValue isReal(const Expr&, const AssumptionContext&);
TruthValue isInteger(const Expr&, const AssumptionContext&);
TruthValue isPositive(const Expr&, const AssumptionContext&);
TruthValue isNonZero(const Expr&, const AssumptionContext&);
TruthValue isFinite(const Expr&, const AssumptionContext&);
```

```cpp
enum class TruthValue {
    False,
    Unknown,
    True
};
```

A small implication closure is sufficient initially.

Examples:

```text
x > 0       ⇒ Real[x], x != 0
x ∈ Integer ⇒ x ∈ Real
x > 3       ⇒ x > 0
```

The same context should be shared by simplification, solving, integration, limits, and certified evaluation.

---

# 10. R4 — Specialized Mathematical Views

## 10.1 PolynomialView

Use dedicated polynomial representations for:

- degree;
- GCD;
- factoring;
- rational-root algorithms;
- derivatives;
- resultants;
- polynomial solving.

## 10.2 RationalFunctionView

Represent `P(x)/Q(x)` directly.

Use it for:

- cancellation;
- holes and poles;
- partial fractions;
- rational integration;
- rational inequalities.

Excluded points must survive cancellation.

## 10.3 LinearSystemView

Convert multivariable linear equations into:

```text
A x = b
```

rather than repeatedly extracting coefficients from generic ASTs.

## 10.4 SeriesData

For future `Series`, `Limit`, asymptotics, and special-function expansions, use a first-class truncated-series representation.

Start with integer-exponent Taylor series if necessary.

---

# 11. R5 — Declarative Rule Layer

The goal is not to remove all C++ algorithms.

The goal is to make the rule inventory explicit.

```cpp
struct RewriteRule {
    RuleId id;
    RuleDomain domain;
    Pattern pattern;
    ConditionFn condition;
    TransformFn transform;
    Cost cost;
    VerificationPolicy verification;
};
```

Possible verification policies:

```text
None
Structural
SimplifyDifference
DifferentiateAndCompare
SubstituteAndCheck
IntervalCheck
```

For integration, `DifferentiateAndCompare` is a natural default for many rule classes.

Rules should also be divided into phases:

```text
Normalize
Canonicalize
AlgebraicSimplify
TrigSimplify
IntegrationRecognition
IntegrationReduction
PostSimplify
```

This reduces rewrite loops and makes cost direction explicit.

---

# 12. R6 — Result and failure model

Recommended common status:

```cpp
enum class AlgorithmStatus {
    Success,
    NotApplicable,
    Unresolved,
    BudgetExceeded,
    Cancelled
};
```

Use reason codes internally.

For integration:

```text
Success
  certified primitive

NotApplicable
  rule does not match

Unresolved
  no supported result found

BudgetExceeded
  candidate or rewrite explosion

Cancelled
  user interruption
```

This directly helps distinguish:

- mathematical unknown;
- unsupported implementation;
- resource exhaustion;
- user cancellation.

---

# 13. R7 — Numerical backend boundary

This proposal does **not** advocate external numerical dependencies.

The goal is to stop upper layers from depending on representation details.

Within the backend, algorithm selection can continue to evolve:

```text
multiply
  ├─ schoolbook
  └─ Karatsuba

divide
  ├─ classic
  └─ Burnikel-Ziegler
```

Future algorithms such as Toom-Cook, NTT/FFT multiplication, or alternative BigFloat kernels should not require changes in Solve or Integrate.

---

# 14. R8 — Future machine evaluator

This can be delayed.

It must not change literal semantics:

```text
0.1
```

remains exact `1/10`.

Use a separate execution engine for:

- plotting;
- large FFTs;
- large matrices;
- batch statistics;
- numerical initial guesses.

Keep the roles distinct:

```text
CertifiedEvaluator
  correctness and guaranteed rounding

MachineEvaluator
  throughput and exploratory computation
```

Do not label machine results as certified `DecimalApproximation` values.

---

# 15. Test strategy

## 15.1 Regression

Preserve and grow the existing regression corpus.

## 15.2 Property testing

Parser / formatter:

```text
format(parse(format(parse(x))))
==
format(parse(x))
```

Simplifier:

```text
simplify(x) ≡ x
```

Algebra:

```text
factor(expand(p)) ≡ p
```

Calculus:

```text
D(integrate(f,x),x) ≡ f
```

Solver:

Every returned finite root should satisfy the original equation.

Certified numerics:

Every certified interval should contain an independent high-precision reference.

## 15.3 Differential testing

Development-only test oracles may use external libraries or tools without turning them into runtime dependencies.

## 15.4 Grammar-aware fuzzing

Generate random ASTs with controlled:

- depth;
- node counts;
- function categories;
- singularities;
- exact/approximate values;
- complex branches.

---

# 16. Performance measurement

Freeze a benchmark corpus covering:

```text
BigInt
Rational
BigFloat
Pi / exp / log
simplify
factor
D
integrate
solve
matrix
FFT
special functions
formatter
```

Track not only wall time, but also:

- generated nodes;
- rewrite count;
- maximum BigInt bits;
- maximum working precision;
- allocation count.

Performance counters should explain *why* a regression occurred.

---

# 17. Recommended migration order

## Phase 0 — Baseline freeze

- Freeze the v1.5.1 regression corpus.
- Freeze a benchmark corpus.
- Catalog known failures.
- Freeze formatter round-trip cases.

## Phase 1 — Evaluation infrastructure

Add:

- `EvaluationContext`;
- `EvaluationBudget`;
- `CancellationToken`;
- `AlgorithmResult`.

Do not change mathematical behavior yet.

## Phase 2 — Canonical algebra

Implement:

- additive normalization;
- multiplicative normalization;
- exact coefficient collection;
- stable ordering;
- safe factor/exponent handling.

## Phase 3 — AssumptionContext

Add:

- three-valued predicates;
- simple implication closure;
- bounds;
- migration of Simplifier and Solver;
- later migration of Integrate and CertifiedEvaluator.

## Phase 4 — Mathematical views

In order:

1. PolynomialView
2. RationalFunctionView
3. LinearSystemView
4. SeriesData

## Phase 5 — Declarative rules

Migrate a limited, well-understood subset of integration rules first.

Do not convert the entire simplifier at once.

## Phase 6 — Property and fuzz infrastructure

Add random-expression property testing to CI.

## Phase 7 — Machine engine

Add only when concrete performance use cases justify it.

---

# 18. Priority table

| Item | Priority | Reason |
|---|---:|---|
| Canonical Algebra | A+ | Existing equality/integration symptoms |
| EvaluationContext / Budget | A | Shared termination and resource model |
| AssumptionContext | A | Shared foundation for domain and branch safety |
| PolynomialView | A- | Benefits solve/factor/integrate |
| RationalFunctionView | A- | Required for safe hole/pole-aware rational algorithms |
| AlgorithmResult / failure reasons | A- | Distinguishes inability, unknown, and resource exhaustion |
| Declarative rule tables | B+ | Needed before integration rules proliferate further |
| Property / fuzz testing | B+ | High leverage at current project size |
| Numeric backend boundary | B | Long-term maintainability |
| SeriesData | B | Needed for Series/Limit/special functions |
| MachineEvaluator | B- | Performance layer, not semantic foundation |
| Stable public library API | C | Not urgent for a CLI-only project |
| General pattern language | D | High complexity, currently unnecessary |
| GUI / Notebook | Non-goal | Not central to mmCal's value |

---

# 19. Risks of the refactoring itself

## 19.1 Over-abstraction

The primary risk.

Only introduce a new layer when the same real problem already exists in multiple places.

## 19.2 Canonicalization may break semantics

Particularly dangerous around:

- branch cuts;
- zero denominators;
- `0^0`;
- complex powers;
- logarithmic identities.

Canonical algebra must define a deliberately conservative safe subset.

## 19.3 A rule system can become a second evaluator

Keep rules limited to:

- match;
- condition;
- transform;
- verify.

Do not evolve the rule system into a general programming language.

## 19.4 Context-object bloat

Avoid a giant option structure.

Prefer a small `EvaluationContext` that references separate services:

```text
AssumptionContext
EvaluationBudget
CancellationToken
Diagnostics
```

---

# 20. Success criteria

The refactoring should not be evaluated by the number of new functions.

## Correctness

- Exact semantics preserved.
- No branch/domain regressions.
- No hole/pole regressions.
- No certified-rounding regressions.

## Architecture

- Polynomial algorithms rarely traverse generic AST directly.
- Assumption queries are centralized.
- Major algorithms can observe a shared budget.
- Rule applicability is separate from rule verification.

## Reliability

- Pathological inputs terminate with reasons rather than crashing or overflowing the stack.
- Random-expression tests can detect rewrite loops.
- `Unresolved` and `BudgetExceeded` are distinct.

## Maintainability

- Adding an integration rule requires less modification of existing recognizer chains.
- New polynomial algorithms do not know Parser or Formatter details.
- Numeric backend changes do not propagate into Solver or Integrate.

---

# 21. Final architectural direction

mmCal cannot simply remain as small as Eigenmath; its numerical tower, solver, integrator, and certified evaluation have already exceeded that ceiling.

It should also avoid growing into a single giant `MathStructure`-style central object.

The recommended middle path is:

> **Keep the AST small.  
> Separate mathematical domains into specialized views.  
> Separate mathematical knowledge into algorithms and controlled rules.  
> Separate exact semantics from execution policy.  
> Pass shared contexts, budgets, and assumptions through all serious computations.**

The core value of mmCal is not the size of its function catalog.

It is:

```text
exact-first
soundness-first
self-contained
```

The refactoring should make those qualities easier to preserve as the project grows.

---

# 22. Design principles to carry forward

### From Eigenmath

- Keep the whole system understandable.
- Make runtime limits explicit.
- Simplicity is an architectural asset.
- Avoid unnecessary first-class concepts.

### From Qalculate!

- Maintain a real backend boundary.
- Assumptions are kernel infrastructure, not an optional feature.
- Distinguish intervals as mathematical objects from intervals used for error tracking.
- Termination, cancellation, and crash resistance are first-class features.
- Random-expression testing is highly valuable in a mature calculator/CAS.
- Keep frontend lifetime separate from mathematical-kernel lifetime.

### From Mathematica 1.x

- Mathematical knowledge can live in a rule layer rather than only in hard-coded algorithms.
- Generic transforms such as affine lifting can be more valuable than many individual formulas.
- Dedicated symbolic representations such as series objects are powerful.
- Unrestricted rewriting inevitably creates termination and performance problems.

### What mmCal should preserve

- exact semantics for finite decimal literals;
- distinction between real-root functions and principal complex powers;
- preservation of definedness, holes, and poles;
- verification of candidates before accepting them;
- strict separation of exact, certified, and machine results;
- honest `Unknown` / unresolved failure semantics.

---

# 23. Reference basis

Primary materials consulted for this document:

- mmCal v1.5 development handoff material and integration catalog;
- Eigenmath GitHub repository (`georgeweigt/eigenmath`)
  - repository history;
  - `src/defs.h`;
  - `src/eval.c`;
- Qalculate! / libqalculate GitHub repository (`Qalculate/libqalculate`)
  - source tree;
  - current TODO;
  - project news and release history;
- static structural analysis of Mathematica 1.2.2f33 Enhanced.

External repository commit counts and project status are observations made on **2026-08-12** and will naturally change over time.
