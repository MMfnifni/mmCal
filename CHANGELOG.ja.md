# Changelog

## v1.5.5 — Unreleased

### Series・漸近展開

- `SeriesData` / TPSAを核とする`series[...]`を追加した。exact Taylor / Laurent / Puiseux展開とlogarithmic係数層を同一表現で扱い，`D` / `integrate` / `normal`へ直接接続する。初等函数に加え，`erf` / `erfc` / `Si` / `Ei` / `Ci` / `li` / Fresnel函数，principal逆三角函数，Lambert W，`gamma` / `lgamma` / `digamma` / `trigamma`，`polylog`の対応可能な局所展開を係数漸化式とDLMFの公式からexactに構成する。`tan/cot/sec/csc`，双曲線系，`expm1/log1p`，cardinal函数，`log2/log10`も既存TPSAへの正規化で接続した。
- `series[expr,{x,Infinity,n}]`を追加した。`Infinity`は実軸の`+Infinity`として`t=1/x`の局所展開へ写し，有理函数，多項式成長，reciprocal composition，Puiseux，logarithmic asymptoticsを既存Series演算で扱う。正の先頭方向を証明できる`log[A(x)]`は`log[c]+r log[1/x]+log[1+h]`へ分解する。振動型，essential growth，logの負冪を要するtransseries等は推測せず未評価に保つ。
- 既存`limit`で決まらない有限点に限り，局所Seriesの先頭項と特異項相殺を利用する限定fallbackを追加した。`toNormal[expr]`はlist / array / call内の`SeriesData`を再帰的に通常形へ戻し，finite / conditional `SolutionSet`ではbinding右辺だけを変換して条件，自由変数，multiplicity，domainを保持する。
- `N[SeriesData]`が指数格子metadataまで有限precision化して後続の`normal` / `D`を壊す問題を修正した。係数とcenterだけを近似し，`minimumExponent` / `orderNumerator` / `exponentDenominator`はexact integerとして保持する。
- `+Infinity`展開で`1/(1/t)`等の未整理reciprocalがvaluationへ入り，`tan[1/x]`，`sec[1/x]`，`cot[1/x]`，`csc[1/x]`等が誤って`Division by zero`になる問題を修正した。局所変数への置換直後だけ安全なreciprocal正規化を行い，global Simplifierの意味論は変更しない。
- `gamma`の正則Seriesで，`lgamma`から指数Seriesへ戻す係数漸化式が各中間項でSimplifierを反復し，次数とともに式が急膨張する経路を軽量化した。係数ごとに生の和を構築してから一度だけ簡約し，低次数のcanonical出力を維持したまま秒級の性能段差を緩和した。
- affineなexact Rational引数の`Ei` / `li` Seriesへ係数漸化式を追加し，汎用Series積・逆数と反復Simplifierによる秒級の膨張を避けた。`li[2+sqrt[x]]`のようなPuiseux compositionも局所変数上のaffine形を利用して直接構成する。

### Array・Vector・ベクトル解析

- `Array / scalar`をelementwise scalar scalingとして追加した。exact zero除算でもscalarのexceptional-value意味論を成分ごとに保持し，`{0,1}/0 -> {Indeterminate, ComplexInfinity}`となる。`scalar / Array`と`Array / Array`は線形代数上の意味を暗黙に仮定せず，引き続き拒否する。
- Vector APIへ`inner[a,b]`，`outer[a,b]`，`distance[a,b]`，`manhattanDistance[a,b]`，`projection[a,b]`，`rejection[a,b]`を追加した。`dot`はbilinear contractionを維持し，`inner`は第1引数を共役するHermitian内積とする。`norm`，距離，射影，棄却成分，反射は同じHermitian意味論へ統合した。
- rank-2 Arrayの各行をvector集合として扱う`orthogonalQ`，`orthonormalQ`，`linearIndependentQ`，`gramSchmidt`を追加した。Hermitian内積に基づくexact-firstの直交化を行い，直交化中は非正規化basisを保持して根号を含む途中式の膨張を抑える。provably dependentな行は落とす一方，symbolicなzero/nonzero判定を証明できない場合はbasisを推測しない。
- 旧Vector convenience名を整理した。`vadd` / `vsub` / `vscalar` / `vsum`は公開登録から外し，Arrayの`+` / `-` / scalar multiplicationと`sum`へ一本化した。`vcross` / `vmanhattan` / `veuclidean` / `vproject` / `vangle` / `vreflect` / `vreflect_axis`も廃止し，`cross` / `manhattanDistance` / `distance` / `projection` / `vectorAngle` / `reflectNormal` / `reflectAxis`を正式名とする。
- `vectorAngle[a,b,assumptions]`を追加し，symbolic vectorの実数性を明示できるようにした。complex vectorへ幾何学的angleを勝手に一般化せず，実vectorだけを扱う。
- Cartesian座標の微分演算子として`grad`，`divergence`，`curl`，`laplacian`，`jacobian`，`hessian`，`directionalDerivative`を追加した。`curl`は2次元ではscalar curl，3次元ではvector curlを返す。`laplacian`はvector fieldへcomponentwiseに拡張した。円筒・球座標のscale factor等は暗黙に推測しない。
- `at`をfinite `SolutionSet`へ拡張した。`at[solutions,i]`は0始まりで1 branchだけを含む`SolutionSet`を返し，条件，自由変数，multiplicity，solver変数domainを保持する。`at[solutions,i,x]`は指定symbolのbinding右辺を取得する。非finite集合，範囲外index，未binding変数は明示的に拒否する。

### 記号計算・CLI・FFT

- 主値`sqrt`方程式で右辺がgeneric complex parameterの場合も完全な解集合を返せるようにした。`sqrt[x]==y`は`x==y^2`に加え，principal square rootの像条件`re[y] > 0`または`re[y] == 0 && im[y] >= 0`をSolutionBranchへ保持する。exact complex RHSは値域判定前にboundedなexact arithmeticで正規化し，`sqrt[x]==-2+3I`のような像外入力を偽の条件付きbranchへ残さない。
- Real `solve`でprincipal `asin` / `acos` / `atan` / `acosh`の逆向き方程式を値域条件付きで反転するようにした。三角逆函数の値域は現在のangle modeに従い，endpointはclosed/openを区別する。`log2[x]==y` / `log10[x]==y`も`y in Real`をSolutionBranchへ保持してsymbolic parameterを完全に扱う。
- principal逆函数のexact endpoint判定を強化し，`solve[asin[x]==-Pi/2,x,Real]`等で自明なrange条件を解に残さず，open endpointは逆函数本体を評価する前に空集合として確定する。
- complex大引数`Ci`の保証付き評価で，到達不能な漸近精度を事前判定して無駄な`E1`漸近試行とguard再試行を避けた。`N[Ci[140+I],100]`のようなbackend境界でfrontend timeoutへ落ちる性能崖を解消した。
- exact入力の`N[...]`では有限precision由来の情報幅が存在しないことを利用し，`CertifiedEnclosure`と`InformationEnclosure`のために同じ特殊函数backendを二重評価する固定費を除いた。finite-precision値，nested `N`，外部bindingを含む式は従来どおり別々のenclosureを伝播する。exact real `2F1`が`z=1`へ十分近い場合はDLMF 15.8.4の`1-z`接続をguard付きで利用し，exact Rationalな接続係数はRational専用Gamma / LogGamma backendへ直接送り，near-unit Gauss級数とgeneric interval Gammaの固定費を避ける。
- 保証付き特殊函数の実軸dispatchをさらに整理した。exact Rationalの`zeta` / `digamma` / `trigamma`は元の有理引数を保持してdyadic端点の二重評価を避け，負のfinite-precision `digamma`は実recurrenceで正実軸backendへ接続する。exactな半整数`Pi`倍の`ellipticF` / `ellipticE` / `ellipticPi`はperiod reductionで同じcomplete Carlson積分を重複評価せず，complete値の整数倍として構成する。
- exact有理函数積分のalgebraic-log経路を軽量化した。留数`P(r)/Q'(r)`を同じ`Root`上のexact式として保持し，中間`AlgebraicNumber`の反復canonicalizeを避ける。あわせて`x^(4m)+x^(2m)+1`型の低コスト因数分解を追加し，高次有理函数積分の秒級停滞を大幅に縮めた。
- 高階`D[exp[q(x)],{x,n}]`で二次以下の`q`を係数vector漸化式として構成し，nested `D`も内側をmaterializeしてから外側へ接続する。多項式×線形`exp/sin/cos/sinh/cosh`積分も有限係数漸化式へ移し，部分積分の深さ上限と高次数の性能崖を避けた。単純な異周波数`sin/cos`積は汎用探索前に積和公式へ送り，数学的に同等で式木が浅くなる場合は`sin[x]^2/2`より`-cos[2x]/4`のような平坦なprimitiveをcanonicalに採用する。高次数の純binomial radical / reciprocalは既存`2F1` primitiveへ早期dispatchし，`c/log[a*x+b]`はaffine chainを証明できる場合に`li`へ直接接続して，汎用候補探索の固定費を避ける。
- Real多項式`solve`ではreducibleな場合に因子ごとの実根分離を共有し，`(x^2-a)(x^2-b)`型の固定費を削減した。`x^(2m)+b x^m+c`でRational上の二次因数分解が可能な場合も低コストに分解し，既存fast pathへ接続する。有理根だけで実解集合が閉じる場合はRationalを直接返し，`x^4-1 -> {-1,1}`等の既存canonical出力を維持する。direct bindingは一般algebraic proofより先にdomain filterへ送り，`solve[x==Phi,x,Real|Rational]`等の不要なAlgebraicNumber構成を避ける。
- 異なる純二次`Root`同士の加減乗除は，既知のdegree-2/4 annihilating polynomialとexactな共役順序から結果根を直接構成し，一般resultant / primitive-element構成を回避する。`sqrt[2]±sqrt[3]`型のcanonical `root[...]`表現は維持する。
- `D[cases[...]]`が明示default枝の境界へ偽の導函数値を与える場合を修正し，高階`D`の加減算をcompactに正規化した。`simplify` / `fullSimplify`のexact有理係数線形結合と有理affine仮定推論も拡張した。
- 横断的なfrontend合成を修正した。`limit`は`D`，vector calculus，`solve`，`collect` / Gröbner操作，`Series` / `Normal`等の内側binder・制御変数を外側の点代入で破壊せず，安全にmaterializeできるものは先に評価する。変数依存`cases`の二側極限は左右を別々に評価し，一致時だけ値を返し，不一致は`Indeterminate`とする。brace/List値も成分ごとの極限へ接続した。
- `solve`のdefinednessとdispatchの非対称性を修正した。同一式方程式では有限Predicateで表せる定義条件を保持して`All`へ落とし，relation配列内でもscalar relationと同じtranscendental / radical dispatchを利用する。Real ambient domainを正規化へ渡し，`sqrt[x^2]==±x`，`abs[x]==±x`，`exp[log[x]]==x`，`log[exp[x]]==x`等を既存branch/domain知識から解けるようにした。
- principal逆函数の安全な向きの合成を整理し，有限式の`sin[asin[x]]` / `cos[acos[x]]` / `sinh[asinh[x]]` / `cosh[acosh[x]]`を縮約する一方，`tan[atan[x]]` / `tanh[atanh[x]]`は例外点を消さないdefinedness条件付きで扱う。逆向きの`asin[sin[x]]`等は一般には縮約しない。あわせて`x>1`等の仮定から`log/log2/log10`の符号・非零性を伝播し，Vectorのzero-directionエラーが委譲先`projection`ではなく公開函数自身の名前を報告するよう修正した。
- held symbolic frontend間の接続を補強し，`integrate[D[...]]`，`Series[D[...]]`，`solve[D[...]==...]`，`integrate[Normal[Series[...]]]`，nested `Normal[Series]`等を安全な限定materializationで評価する。この接続は通常評価される式木の内部，Vector Calculus，`simplify` / `fullSimplify` / `expand` / `factor` / `collect`にも適用する。式変形はheld引数を既存kernelへ直接渡し，通常のheld callを無差別に先行評価しないため，session definitionやbinder意味論を変えない。
- nested binder frontendの合成を補強した。外側`D` / `integrate` / `series` / `solve` / Vector Calculusから内側`limit` / `integrate`を先に閉じる際，外側binder変数を現在のsession definitionから一時的に保護する。`D[integrate[f,x],x]`は既存の微積分基本定理ruleを優先し，原始函数へ不必要に展開しない。`limit`の点代入後に残る`integrate[0,x]`や別変数のinner `limit`も安全な場合だけ再materializeする。
- REPLへ`:quit` / `:exit`と対話表示用`:layout`を追加した。
- 現在のdegree budget内の2冪exact `ifft[fft[v]]`で，FFT由来のroot-of-unity式を円分体座標へ再埋込みし，forward表現を維持したままinverseの式爆発を抑えた。
- 非2冪exact FFTの円分体再埋込みを補強し，`Q(ζ_12)`で現れる`sqrt[3] = ζ_12 + ζ_12^(-1)`を認識するようにした。6点Gaussian-integerの`ifft[fft[v]]`がgeneric Expr fallbackへ落ちて約1秒を要する性能崖を解消し，exact round-tripを数msで閉じる。
- Random Expression Fuzzerの単一case再現時は，PASSでもdepth，不変条件，生成式を表示するようにし，timeoutや性能崖の反例をfailure化せず直接診断できるようにした。

### 内部リファクタリング

- exact scalarのAST生成・判定，Rationalの整数冪，関係演算子の`BuiltinId`対応，`Cases/CaseBranch`構築を共通helperへ整理した。Series，積分，極限，Solver，Assumption等が同じexact値・relation意味論を共有し，追加時の同期漏れを減らす。
- exact / approximate線形代数の内部表現を整理した。exact `NumberMatrix`と行列四則wrapper，SVD / eigenで使うNearestEven BigFloat点演算を共通化し，shape処理・丸めmode・基本演算の重複実装を削除した。
- Vector / Matrix / statistics等のBuiltin family判定を共通predicateへ集約し，評価器ごとの巨大なID列挙の重複を減らした。一方，算法ごとに再試行単位が異なるcertified retry loopや局所的な一行aliasは，制御境界と可読性を保つため無理に抽象化していない。
- `PrecisionInsufficient`ではguard桁を増やして再試行し，`CertifiedBackendUnsupported`では再試行せずfallbackする等，保証付き線形代数の非自明な制御理由を日本語コメントで明示した。参照されなくなった積分helper等の死コードも削除した。公開API・数式意味論・Formatter契約は変更しない。

## v1.5.4 -2026-08-31

### Solver・記号計算

- 多変数polynomial `solve`を強化した。Gröbner消去と再帰的specializationを拡張し，対応可能なzero-dimensional系をexactに完全列挙するほか，`x y==0`や`x^2+y^2==1`等のpositive-dimensional系も，証明可能な範囲で自由パラメータとexact条件を持つ解集合として返す。一般多様体の不完全なparameterizationは引き続き推測しない。
- Real `solve`を拡張し，`u exp[u]=a`，`exp[p x+q]==c x+d`，正定数底の指数方程式，`x^x==r`の一部，`lambertw[u]==r`，`cosh[u]==r`等をLambert Wや既存inverseへ接続した。値域・単調性・凸性等による不存在／一意性証明も追加し，証明できた場合だけ空集合またはexact解を返す。
- `abs`および主値`radical`方程式をReal solverへ追加した。平方等で導入される偽根をexactに除外し，Complex上の円周等を有限個の解へ誤って潰さない。
- 数学的場合分け`cases[value if condition; ...]`を追加し，`simplify`，`N`，`D`，`integrate`，`limit`へ接続した。未知条件は保持し，false branchは評価しない。
- `groebnerBasis[...]` / `polynomialReduce[...]`を追加した。Lex / GrLex / GrevLexのexact Rational Gröbner基底を扱い，非線形多項式`solve`にも利用する。
- `simplify` / `fullSimplify`を定義性に配慮するよう強化し，`F/F -> 1`，`F^0 -> 1`，特殊函数の退化等は必要な定義条件を証明できる場合だけ適用する。exact複素数は演算後に虚部が0なら実数へ正規化する。

### 微分・積分・極限

- exact有理函数積分を一般化し，Hermite reductionとalgebraic `Root`を用いて高次・重複分母を含む有理函数まで扱える範囲を拡張した。既存の簡潔な`Log/atan/atanh`表現や初等的置換を優先する。
- 仮定付き定積分・広義積分を拡張し，Gamma/Beta/Mellin，Frullani型，対数moment，`1+x^q`型，Dirichlet/Fresnel等のexact条件付き評価を追加した。収束条件や内部特異点を証明できない場合は未評価を保持する。
- `D`の`cases[...]`処理と高階微分を強化した。変数依存の境界では非微分可能点へ偽の値を割り当てず，Lambert W，`polylog`，二次式指数函数等の高階微分をよりcompactなexact式で構成する。
- `limit[expr,{x,a,direction}]`を追加し，`Ei` / `Ci` / `li`や実函数の既知極限，周期振動，squeeze可能な極限を補強した。積分のderivative-back検証も共通定義域上の恒等性を確認するよう修正した。

### Algebraic Root・高次多項式

- 高次Complex `Root`と一般高次polynomial `solve`の性能段差を大幅に改善した。root isolation / refinementと既約性判定の重複計算を減らし，対称多項式やroot半径が大きく異なる場合の停滞も解消した。96次までの監査範囲を拡張している。
- AlgebraicNumber / number field内の演算を高速化し，中間`Root`生成を抑えてexact性を保ったまま高次冪等の負荷を低減した。

### 保証付き `N`・特殊函数

- 複素数を含む保証付き`N`を大幅に拡張した。`erf/erfc`，`Ei/Si/Ci`，Fresnel函数，`1F1`，`2F1`，`polylog`，`zeta`，`gamma`，`digamma/trigamma`，Lambert W，楕円積分等で主値枝と誤差保証を維持した評価範囲を広げた。
- 従来の固定的な大きさ境界を多数撤去し，実数`Ei/Si/Ci`，複素`Ei/Ci`，Fresnel，`1F1`，`2F1`，`polylog`，`ellipticF/E/Pi`等を収束性と資源上限に基づいて選択的に評価するようにした。`zeta`は`s=1`を除く有限複素平面へ，`ibeta`はfinite-precisionを含む正実パラメータへ対応範囲を拡張した。
- Lambert Wは`-1/e`近傍と任意整数branchの保証付き複素評価を追加した。finite-precision入力で分岐側を証明できない場合は推測せず`N::precision`を返す。
- `CertifiedEnclosure`と`InformationEnclosure`の役割を整理し，有限precision入力の隠れた保護桁から過剰な精度を捏造しないよう統一した。特異点・分岐切断・算法境界の曖昧性は，数学的なDomainErrorと区別して`N::precision`または`N::unsupported`として扱う。
- `N`の表示を簡潔化し，近零値や末尾0を整理しつつapproximation provenanceを保持する。nested `N`，`diff[...,digits]`，`nintegrate[...,digits]`も共通のprecision制約と情報量制約へ統一した。
- `acosh`の分岐切断近傍，`expm1/log1p`の0近傍，`2F1`の接続公式，複素`li`等の相殺・過剰refinementによる性能問題を修正した。

### Exact線形代数

- `conditionNumber`，`pseudoInverse`，`leastSquares`を追加した。exact行列ではrankをexactに判定し，外側`N[...]`では保証付きSVDを利用する。
- exact Integer/Rational行列の`det` / `solveLinear`へmodular + CRT経路を追加し，検証付きでBareissと自動選択する。`leastSquares`の中間丸めによる精度損失も修正した。

### 評価資源・CLI

- 評価ごとの共通`EvaluationBudget` / `EvaluationLimits`を追加した。巨大整数，Array/Matrix，Solver，積分，AlgebraicNumber，要求precision等を計算前または計算中に制限し，超過時は`ResourceLimitError`として報告する。
- Ctrl-C / Ctrl-Break / SIGINTによる協調的な評価取消しを追加した。
- 行単位自動化向け`--batch`を追加した。REPLの`:help`も全callable函数，定数，入力規則，例，近似候補を含む一覧へ拡張した。

### 検証・互換性

- public CLIだけを使うblack-box監査とRandom Expression Fuzzerを拡張し，`D` / `integrate` / `limit` / `N` / `cases` / Gröbner / Solve / exact FFT / 線形代数 / precision provenance等の性質検証を追加した。
- Formatter，履歴参照，`cases[...]`を跨ぐ微積分，主値分岐，可除特異点，DomainError分類等の回帰を修正した。公開構文とexact-first方針は維持する。

## v1.5.3 — 2026-08-16

v1.5.3は，v1.5.2の記号微積分，保証付き`N`，Array・線形代数を基礎に，**exact代数数体，Solverの意味論統一，高精度特殊函数，exact Cyclotomic FFT，内部表現の最適化**を接続した版である。exact-first・proof-onlyの方針を維持し，証明不能や未対応を推測値として確定しない。

### 代数数・数体

- 実数／複素数`root[...]`を共通`AlgebraicNumber`へ統合し，最小多項式簡約，原始元，終結式，保証付きの根再同定を有界計算で実装した。
- 不変な`NumberFieldContext` / `AlgebraicElement`を導入し，選択した埋め込みとRational冪基底座標を保持する。同一数体の四則演算，小整数冪，exact逆元を数体内で処理する。
- 数体・逆元・最小多項式等の再利用キャッシュを追加し，最小多項式はexact Krylov消去で導出する。exact代数数の等値判定と実数上の大小比較も追加した。
- `root[...]`，exact Rational/複素Rational，`sqrt` / `cbrt`，`Phi`等を共通の代数数評価へ接続し，比較・定義域判定・Solveで同一のexact値として扱う。表示は従来どおり`root[minpoly,k]`を基本とする。

### Solver・意味論

- `solve[equation,Integer|Rational|Real|Complex]`を追加し，未知user symbolが一意な場合だけ変数を推定する。protected symbolはsolve変数にしない。
- `E^x` / `exp[x]`，`ln` / `log`等の表現差を正規化し，証明可能な正の底の指数方程式を対数で反転する。exact記号`lambertw`と実数枝`k=0/-1`の保証付き`N`も追加した。
- 非整数・非Rationalであることを`ValueFacts`へ伝播し，`N[SolutionSet]`は解集合構造を保ったまま数値化できる右辺だけを近似する。比較演算子の表示も読みやすく整えた。

### 保証付き数値評価・特殊函数

- `DecimalApproximation` / `ComplexDecimalApproximation`を第一級の保証付き数値として扱い，`CertifiedEnclosure`と`InformationEnclosure`を後続演算へ伝播する。
- `N[expr,p]`を有効10進桁数として統一し，全体を保証できない場合もHold属性を壊さず数値部分だけを近似する。重複診断も抑制した。
- `zeta` / `digamma` / `trigamma` / regularized `ibeta`を追加し，Gamma/Beta系の高精度計算を高速化した。値は存在するが保証付き計算未対応の場合は`CertifiedBackendUnsupported`としてDomainErrorから分離する。

### Exact FFT・Array表現

- exact Cyclotomic FFTを正式化し，対応する非2冪exact入力を`Q[t]/Phi_n(t)`上で計算する。5/7/10/12点等の往復を巨大な`cis[...]`式なしでexactに閉じ，適用不能時は従来のexact DFTへ戻す。
- `Expr::Node`をkind別typed nodeへ変更し，公開`Expr` APIを保ったまま固定payloadを削減した。
- dense `ArrayExpr`をimmutableなpaged packed storage + shape/offset/strides viewへ変更した。transposeはzero-copy viewとなり，矩形数値braceは`ArrayBuilder`へ直接構築する。

### 函数・補助機能・診断

- `isprime` / `nextprime` / `prevprime` / `factorint` / `totient`を`uint64`全域でdeterministic exactに追加した。証明できないBigIntをprobable-primeだけで`True`にはしない。
- bit演算群，`round[x,n]`，`fma`，`clamp`，`proj`，`range` / `table` / `map`，`explain`，周期的な実数解族を追加した。
- 構文解析に失敗した入力は`In[n]`を消費せず，`Infinity-Infinity`や`0*Infinity`等も誤簡約しない。未使用のMachine/double shimも削除した。

### 性能・検証

- `mmCal.Benchmarks --random-expressions`へseed+caseで再現可能な並列の意味論fuzzerを追加し，長時間実行時の履歴蓄積も防止した。
- 代数数体，特殊函数，exact Cyclotomic FFT等のベンチマークを追加し，採用・棄却した最適化の実測根拠を`docs/performance_optimization.*`へ記録した。利益のないキャッシュや無条件な高精度展開は採用しない。

### 文書・互換性

- README / Reference / Architecture / Roadmap / 性能文書をv1.5.3へ同期し，意図的な未実装・探索上限はroadmapとintentional-limits memorandumへ分離した。
- 内部原始元やcache表現は利用者向けFormatterへ露出せず，代数数の標準表示は`root[minpoly,k]`を維持する。

## v1.5.2 — 2026-08-13

### 構文・REPL（破壊的変更を含む）

- 函数呼び出しを`name[...]`へ一本化し，`name(...)`を廃止した。`()`はgrouping専用で，既知函数へ旧形式を使うとSyntaxErrorとなる。
- `x(x+1)`は暗黙乗算として受理し，Formatterは`x*(x+1)`へ正規化する。
- 履歴参照を`In[n]` / `Out[n]`へ統一した。正添字は絶対番号，負添字は相対参照，`0`は無効で，`@`群は`In[-n]`，`%`群は`Out[-n]`の短縮形である。
- `D` / `integrate` / `limit`等のheld symbolic operator内でも`In[n]` / `Out[n]` / `%`を履歴参照として再帰的に解決し，`integrate[Out[1],x]`等が未知函数扱いになる不整合を修正した。
- 多変数polynomialの正次元代替経路に定数係数の一次変数消去を追加し，複数方程式を低次元systemへexactに落として既存のRealパラメータ-定義域 射影へ接続した。

### 記号微積分・積分規則

- `sin[u]^m cos[u]^n`の非負整数冪を有限Fourier多項式へ還元し，負整数冪は`csc/sec`漸化式，`tan/cot/sec/csc`整数冪は各還元公式で処理する。異周波数積のproduct-to-sumも追加した。
- 共通引数を持つ`R(sin[x],cos[x])`へ有界なWeierstrass置換`t=tan[x/2]`を追加し，exact有理函数積分へ接続した。
- inverse-chain / substitution認識と二次根号族を拡張し，`2x(1+x^2)^5`，`x/(1+x^4)`等を構造的に処理する。
- derivative-back監査を維持しつつ，証明器不足だけで正しい原始函数を捨てない`ResolutionOnly`を明文化した。積分失敗も`unsupported` / `partial` / `conditionsRequired` / `noKnownClosedForm`へ分類した。
- `docs/memorandum/integralCatalog.md`を用い，三角・有理・特殊函数・分岐に敏感な積分族を横断監査した。

### 特殊函数

- `fresnelc` / `fresnels`，`hypergeometric1F1` / `hypergeometric2F1`，`ellipticF` / `ellipticE` / `ellipticPi`を追加し，exact退化，微分，保証付き実数`N`，対応する積分規則へ接続した。
- `Ei` / `Si` / `Ci` / `li` / `polylog`を追加し，代表的なexact値・微分・保証付き実数`N`・積分を実装した。
- `integrate[exp[x^n],x]`等では主値枝を保ちやすい1F1表現を優先する。
- 存在しない一般的な特殊函数の逆函数をSolverで捏造せず，安全なexact退化で既存Solverへ落ちる場合だけ解く。

### 精度対応 `N`・FFT

- exact Cyclotomic FFTを追加し，5点以上の対応する非2冪exact Rational入力を`Q[t]/Phi_n(t)`上で処理する。`ifft[fft[v]]`を巨大なroot-of-unity式なしでexactに閉じ，次数上限や数体所属を証明できない場合は従来のexact DFTへ戻す。
- `N[expr,p]`を精度対応の評価入口へ拡張し，対応函数へ`ApproximationContext`を伝播する。`N[fft[data],p]`は巨大exact式を作らずBigFloat/`ComplexInterval`の保証付きFFTへ直接進む。
- 保証付き非2冪FFTは小サイズで直接DFT，大サイズでBluesteinを使う。FFT/Matrix間で区間変換，10進化，保護桁の精密化等を共有した。
- 保証付き近似値は保証情報を保ったまま末尾0列を短く表示し，`N`由来の有限小数・整数にはexact値との区別用に最低1個の0を残す。

### Array・線形代数

- `{...}`を一般有限brace containerへ拡張し，矩形要素はrow-majorのdense `ArrayExpr`へ自動最適化した。0長次元と空Arrayのshapeも保持する。
- `dimensions` / `arrayRank` / `length` / `at` / `reshape`を整備し，`MatrixView` / `MatrixBuffer`で不要な複製を減らした。主要APIは`dot` / `matrixRank` / `norm` / `normalize`へ統一し，旧名は互換aliasとして残す。
- integer/Rational行列へBareiss fraction-free eliminationを導入し，`det` / `rref` / `matrixRank` / `inverse` / `solveLinear` / `nullSpace`で共有した。`solveLinear`は一意解だけを返し，`nullSpace`は決定的な基底と空shapeを保持する。
- `luDecomposition[A]`を追加し，exactでは行ピボット付き`P A = L U`，`N[...]`では保証付きinterval partial pivotingを使う。
- `qrDecomposition[A]`を矩形reduced Householder QRとして追加した。一般exact QRは式爆発のため3×3以下に制限し，上三角／上台形には高速経路を使う。自動block化は実測で優位性が安定せず不採用とした。
- exact `qrDecomposition[A]`の一般経路をExpr上のHouseholder展開からfraction-free直交化へ置換した。直交化中は平方根・除算を生成せず，primitive整数vectorとGCD content reductionで進め，full-rank caseはGram行列の対称Bareiss（fraction-free LDLᵀ相当）をfast pathに使う。従来の3×3 hard capを撤去し，rank落ちは直接fraction-free経路へ切り替える。
- exact Matrixのperformance-cliff再監査で，`inverse` / unique `solveLinear`のBareiss後段をgeneric Rational RREFから共通分母BigInt back-substitutionへ置換した。full-column-rank `rref`もRational backward phaseを省略する。modular inverseは32×32・256-bitまで再計測して引き続きBareissより遅く，自動dispatchは有効化しない。
- 実・複素reduced `svd`を追加し，`A^H A`を形成せずHouseholder bidiagonalization + one-sided Jacobiを使う。`conjugateTranspose`も追加した。
- `eigenvalues` / `eigenvectors` / `eigensystem`を追加した。簡単なexact行列はexactのまま，一般`N[...]`はComplex BigFloat Hessenberg + shifted QR + Schur経路で処理する。
- 記号`det` / `inverse`には三角行列の高速経路と展開上限を設け，rank/null-space等の不連続量は任意epsilonではなくexact pivotまたは区間で証明できる場合だけ確定する。

### 性能・安定性・開発基盤

- internal test runnerは`--timings`指定時だけsuite別実時間を表示し，通常実行の出力契約は維持する。
- Bareissにより8～16次で`det`約5.4～11.4倍，`rref`約7.7～15.5倍の改善を実測した。`--matrix-large`も追加し，32/64次行列と1024次の保存・parse負荷を測定した。
- 1024×1024のRational行列ではExpr/Rational表現が先にメモリ・構文解析のボトルネックになることを確認し，型別`Expr::Node`，BigUInt/BigInt SBO，packed Array等を次版候補として記録した。
- MSVCで負の`RealInterval::point`が反転し得る評価順依存と`<stdexcept>`の推移的include依存を修正し，Matrix / FFTのfixed-seed invariantを拡充した。

### 文書・ライセンス

- README / Reference / Architecture / Roadmap / 性能文書をv1.5.2へ同期した。
- BSD 3-Clauseの著作権表示をプロジェクトmetadataと揃え，商標・ブランド利用は`TRADEMARKS.md` / `TRADEMARKS.ja.md`の別方針として明記した。

## v1.5.1 — 2026-08-12

v1.5.0のexact-first CAS基盤を維持しつつ，正規化，検証，多倍長整数，高精度数値評価，ベンチマーク基盤を重点的に改善した。

### 数式・CAS

- `Add`へAST全域の決定的な全順序を導入し，積・除算も定義性を保つ正規形へ整理した。`MathKnowledge`の非零知識と関係式の左右反転推論も強化した。
- Formatter/Parser 往復，積分derivative-back，Referenceとbuiltin registryの照合，極端値近似のテストを追加・拡充した。
- radix prefix衝突，Array隣接，負Rational表示を修正し，FFT plan/twiddleのtransform間cacheを追加した。

### 多倍長・高精度

- BigUInt乗算をschoolbook / Karatsuba / Toom-3の適応選択へ変更し，専用square，Burnikel–Ziegler除算，`2^k`除算，decimal divide-and-conquer変換を追加した。
- factorialのbalanced product treeを高速化し，`tryToUint64`の巨大値早期棄却，BigFloatの巨大指数差加減算の高速経路を追加した。
- `Pi`をbinary-splitting Chudnovskyへ，`exp/E`・`log`をbinary splitting + 保証付き範囲縮約へ変更し，巨大Radianの`sin/cos/tan`にも保証付き引数縮約を追加した。

### ベンチマーク・テスト

- Visual Studio solutionへ`mmCal.Benchmarks`を追加し，固定seedのrandom invariant，算法閾値，factorial，decimal I/O，高精度`Pi/exp/log`のベンチマークを常設化した。
- v1.5.1確定時点でinternal 1691 / 1691，ブラックボックス 1337 / 1337を確認した。

### 比較したが採用しなかったもの

- Prime-Swing factorial，binary GCD，Karatsuba workspace各種，Toom-3専用square／低閾値化，machine `fmod`による巨大trig縮約は実測上の利点が不足したため採用しなかった。

各判断の実測根拠は`docs/performance_optimization.ja.md`を参照。

---

## v1.5.0

旧版から数値モデル，Lexer/Parser/AST/Evaluator，Simplifier，Solver，CertifiedEvaluator，CLI，Formatter，テスト，文書をほぼ全面再構築し，mmCalをexact-firstなCLI電卓 / 小規模CASとして再定義した版である。
