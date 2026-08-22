# Unified EvaluationBudget

## 1. 目的

`EvaluationBudget`は，1回のtop-level評価が使用できる決定論的資源量を一つの文脈で管理する。数学的な定義域，解の有無，backendの対応範囲とは別概念である。

上限超過はすべて`ResourceLimitError`である。`DomainError`，空の`SolutionSet`，`Unresolved`，未評価式へ読み替えてはならない。

## 2. 所有期間

- `KernelSession::evaluate`は入力ごとに新しいbudgetを作る。
- 成功・失敗にかかわらず，次の入力へcounterを持ち越さない。
- 複数の`KernelSession`やfuzzer worker間でbudgetを共有しない。
- `EvaluationContext`がbudgetへの非所有pointerを持つ。深いsubsystemではscope付きthread-local参照を使うが，scope終了時に以前の参照を必ず復元する。
- Coreはwall clockを参照しない。最大入力byte数はcopy・tokenize前に決定論的に検査し，deadline / UI cancelはfrontendが`EvaluationCancellationToken`を発火して伝える。CLIはactiveなtop-level評価中だけCtrl-C / Ctrl-Break / SIGINTをtokenへ接続する。

## 3. 既定上限

| `EvaluationLimits` field | 既定値 | 意味 |
|---|---:|---|
| `maxInputBytes` | 16 MiB | source textをcopy・tokenizeする前の最大byte数 |
| `maxEvaluationSteps` | 10,000,000 | explicit evaluator machineが処理するtask数 |
| `maxDepth` | 1,024 | Expr・symbol解決・user函数の論理評価深度 |
| `maxGeneratedNodes` | 5,000,000 | 評価・探索が生成するExpr nodeの累積量 |
| `maxSimplificationCandidates` | 2,000,000 | Simplifier passとFullSimplify候補 |
| `maxSolverBranches` | 50,000 | Solver strategy・生成branch・constraint適用 |
| `maxIntegrationCandidates` | 100,000 | 積分再帰・置換・部分積分候補 |
| `maxCertifiedRefinements` | 1,000,000 | precision retry，保証付き級数項，argument reduction |
| `maxDenseArrayElements` | 10,000,000 | 評価対象となるdense Array要素の累積量 |
| `maxTemporaryMatrixElements` | 25,000,000 | Matrix作業buffer・symbolic minorの累積量 |
| `maxBigIntegerBits` | 8,000,000 | exact整数・Rational成分の最大bit長 |
| `maxRequestedPrecisionDigits` | 100,000 | `N[expr,p]`の要求10進桁数 |
| `maxAlgebraicDegree` | 64 | root / AlgebraicNumber構成次数 |
| `maxAlgebraicRefinements` | 100,000 | isolating interval / disk refinement回数 |

値は安全性の絶対保証だけでなく，既存の実用入力を不必要に拒否しない初期policyである。公開後は負荷測定とfailure corpusを根拠に調整し，変更時はchangelogへ記録する。

## 4. 局所上限との関係

共通budgetは既存の局所algorithm上限を撤去しない。例えば次を維持する。

- Parserのtoken / AST node / nesting / operator chain上限
- Simplifierの最大pass数とFullSimplifyの1回あたり候補数
- integrateの再帰深度と置換候補数
- symbolic determinantの展開数
- Matrix / SVD / Eigenのprecision retry数
- AlgebraicNumber内部の次数・分離上限

局所上限が先に尽き，正しい未評価fallbackを返す既存契約は維持する。一方，共通budgetが先に尽きた場合はfallbackへ偽装せず`ResourceLimitError`を送出する。

## 5. BigInt・精度・allocationの事前検査

整数冪とfactorialは，結果を作ってからbit長を測るだけではmemory防御にならない。この二経路では入力から結果bit長の下界を求め，上限超過が確実な場合は演算前に停止する。

通常のexact四則演算，packed Array，SolutionSet binding等は，生成結果の整数・Rational numerator / denominator・complex成分を走査し，実bit長を記録する。要求精度は`ApproximationContext`を構築する前に検査する。

`range` / `table`，`zeros` / `identity`，`dot`出力およびMatrix作業bufferは，vectorの`reserve` / `assign`より前に要素数を算出してbudgetへ請求する。結果を確保してから上限超過を報告する経路にはしない。入力のpacked / mixed Arrayについてもexact成分を走査し，数値literalのbit長を評価開始時に検査する。

## 6. APIとtelemetry

```cpp
mmcal::kernel::KernelSession session;
auto limits = session.evaluationLimits();
limits.maxEvaluationSteps = 100'000;
limits.maxBigIntegerBits = 1'000'000;
session.setEvaluationLimits(limits);

auto result = session.evaluate("factor[expand[(x+1)^20]]");
const auto& usage = session.lastEvaluationUsage();
```

deadlineを持つfrontendは時計をCoreへ渡さず，期限到達時にtokenを発火する。

```cpp
mmcal::evaluation::EvaluationCancellationToken cancellation;
// 別threadまたはfrontend event handler: cancellation.requestCancellation();
auto result = session.evaluate("longRunningExpression", cancellation);
```

`lastEvaluationUsage()`は直前の成功評価だけでなく，`ResourceLimitError`を含む失敗評価についても取得できる。fuzzerはこれを反例出力へ含める。`EvaluationUsage::modularPrimes`はexact modular Matrix backendが試行したprime image数を記録するtelemetry専用counterであり，独立したlimitではない。

代表負荷の計測には次を使う。

```text
mmCal.Benchmarks --budget-telemetry
```

代数式，modular determinant，modular solve，certified評価，積分について主要counterと`modular-primes`を表示する。CLI cancellationはcooperativeであり，長時間kernelが`checkEvaluationCancellation()`等のpollへ到達した時点で停止する。

`setEvaluationDepthLimit` / `evaluationDepthLimit`は互換APIとして残り，`EvaluationLimits::maxDepth`だけを変更・取得する。
