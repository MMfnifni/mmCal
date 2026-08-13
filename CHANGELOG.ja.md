# Changelog

## v1.5.2 — 開発中

### 構文（breaking change）

- 函数呼び出しを `name[...]` に一本化し、`name(...)` 呼び出しを廃止
- 丸括弧 `()` はgrouping専用とし、通常identifierの `x(x+1)` は暗黙乗算として扱う
- 既知の函数名に旧 `sin(x)` 構文を使った場合は、誤って乗算へ解釈せずSyntaxErrorを返す
- Formatterは函数を常に `name[...]`、identifierとgroupの積を `x*(...)` と明示してround-tripを一意化
- Parser/AST/Lowererから函数呼出delimiterの分岐を削除
- 履歴参照を`In[n]` / `Out[n]`へ整理し，負添字による相対参照を追加（`0`は無効）
- `@` / `@@` / ... を`In[-1]` / `In[-2]` / ... の短縮記法として追加し，`%` / `%%` / ... と`Out[-n]`の対応も明文化
- 負の`Out[-n]`は成功出力，負の`In[-n]`は入力slotを基準に数える


## v1.5.1 — 2026-08-12

v1.5.0のexact-first CAS基盤を維持しつつ、canonicalization、検証、巨大整数、高精度数値評価、benchmark基盤を重点的に改善した。

### 数式・CAS

- `Add`へAST全域の決定的strict total orderingを導入
- 積・除算をdefinednessを保つcanonical normal formへ整理
- `MathKnowledge`の非零知識と関係式左右反転推論を強化
- formatter/parser生成型round-trip testを追加し、radix prefix衝突、Array隣接、負Rationalの表記不安定を修正
- 積分rule/familyを横断するderivative-back harnessを追加
- Referenceとbuiltin registryの自動照合を追加
- approximation extreme-value testを拡充
- FFT plan/twiddleのtransform間cacheを追加

### 多倍長・高精度

- BigUInt乗算をschoolbook / Karatsuba / Toom-3の適応dispatchへ変更
- 専用squareを追加
- factorialはbalanced product treeを維持し、leaf構築・1-limb経路を高速化
- Burnikel–Ziegler divisionと`2^k`除算fast pathを追加
- decimal parse/toStringを`10^9` chunk + divide-and-conquer化
- `tryToUint64`の巨大値早期棄却を追加
- BigFloat extreme exponent-gap加減算をdirected roundingを保ったままfast-path化
- `Pi`をbinary-splitting Chudnovskyへ変更
- `exp/E`と`log`をbinary splitting + certified range reductionへ変更
- 巨大Radianの`sin/cos/tan`へcertified argument reductionを追加

### Benchmark / test

- Visual Studio solutionへ`mmCal.Benchmarks`を追加
- fixed-seed random invariant、threshold sweep、factorial、decimal I/O、高精度`Pi/exp/log`のbenchmarkを常設化
- v1.5.1確定時点: internal 1691 / 1691、black-box 1337 / 1337

### 比較したが採用しなかったもの

- Prime-Swing factorial
- binary GCD
- Karatsuba vector-pool workspace
- Karatsuba recursion-depth scratch workspace
- Toom-3専用square
- 低threshold Toom-3
- machine `fmod`による巨大trig縮約

各判断の実測根拠は`docs/performance_optimization.ja.md`を参照。

---

## v1.5.0

旧版から数値モデル、Lexer/Parser/AST/Evaluator、Simplifier、Solver、CertifiedEvaluator、CLI、formatter、tests、documentationをほぼ全面再構築し、mmCalをexact-first CLI calculator / compact CASとして再定義したrelease。
