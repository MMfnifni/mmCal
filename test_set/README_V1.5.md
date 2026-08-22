# mmCal V1.5 ブラックボックステスト

このディレクトリは，mmCal の public CLI を介して exact arithmetic，symbolic evaluation，certified numerical evaluation，calculus，special functions，matrix/statistics，definedness/error semantics を横断確認するブラックボックステスト集合である。

テストファイルは従来の34本から10本へ整理した。元の2012ケースは削除せず保持し，2026-08のcertified numerical auditで重要だった17ケースを追加して，現在は合計2029ケースである。

## 実行方法

`test_set` ディレクトリで実行する。

```powershell
py tester.py --exe "..\..\build\x64\Release\mmCal.exe"
```

実行ファイルを省略した場合は，`MMCAL_EXE` 環境変数と一般的なbuild出力先を順に探索する。

CLI自動処理contractは，通常の対話protocolとは別に検証する。

```powershell
py cli_automation_tests.py --exe "..\..\build\x64\Release\mmCal.exe"
```

このtestは `--eval`，`--batch`，`--bach` 互換alias，stdout/stderr分離，終了code，batchのsession保持とerror後継続，巨大Parser入力が `ResourceLimitError` で安全に停止することを確認する。

## ファイル構成

| ファイル | 主な対象 | ケース数 |
|---|---|---:|
| `test00_legacy_scalar.txt` | 旧fixed表示，syntax，実数，複素数 | 425 |
| `test01_legacy_numeric.txt` | statistics，matrix/linear algebra | 353 |
| `test02_legacy_errors.txt` | 旧error/exception回帰 | 165 |
| `test03_language_runtime.txt` | core syntax，angle mode，numeric lexer，history，iteration，variables | 161 |
| `test04_exact_arithmetic_algebra.txt` | BigInt/Rational，exact complex/algebraic，number theory，number fields，cyclotomic FFT | 224 |
| `test05_approximation_certification.txt` | `N`，Precision/Accuracy，explain，InformationEnclosure，branch/backend audit | 163 |
| `test06_calculus_solver.txt` | symbolic/numerical calculus，Solve，Lambert W | 175 |
| `test07_special_functions.txt` | gamma系，hypergeometric，elliptic，polylog等 | 108 |
| `test08_boundary_errors.txt` | 詳細なdomain/boundary/error message contract | 82 |
| `test09_semantics_audit.txt` | semantic coherence，definedness-preserving simplification，横断audit | 173 |

合計は2029ケースである。

## `# @session` と隔離性

物理ファイルを統合しても，従来の「テストファイルごとに新しいmmCal processを起動する」隔離性は失わない。`tester.py` は次のdirectiveを認識する。

```text
# @session
```

これ以降を新しいmmCal processで実行する。起動引数も同時に指定できる。

```text
# @session --fix 15 --angle deg
```

したがって，変数定義，`angleMode`，history，diagnostics等は隣のsessionへ漏れない。

従来形式の

```text
# @args --fix 15 --angle deg
```

も互換性のため引き続き利用可能である。ただし，そのsession内でテストケースが始まる前に置く必要がある。ケース開始後に別の起動条件へ切り替える場合は `# @session` を使う。

## テスト記法

```text
expression ==> expected
expression =>  expected
expression =>> expected
```

`==>` は正規化後のexact比較，`=>` はstrict numeric比較，`=>>` はloose numeric比較である。`Error` または具体的な `DomainError` 等も既存どおり利用できる。

## 追加したcertification回帰

今回の整理では，単なるファイル結合だけでなく，直近のcertified numerical auditで実際に発見された欠陥をblack-box側にも固定した。

- `N[Ci[120+I],20]` のwhole-complex precision certification
- finite-precision negative-real `digamma` / `trigamma` のcomplex recurrence
- finite-precision negative-real `Ci` のprincipal complex value
- negative/complex `li` のprincipal backend
- closed `N[...]` を含む `diff` / `nintegrate` のprecision provenance
- inverse-function branchから安全に離れた領域でのmean-value enclosure精度
- complex `z` に対するterminating `1F1` / `2F1` のexact polynomial化
- exact `2F1` principal-cut pointのcomplex continuation
- approximate terminating valueの `0.0` / `2.0` 表示contract

branch ambiguityそのものの `N::precision` / `N::unsupported` 分類は，public CLIの文字列表現だけでなく内部diagnostic codeを検査する必要があるため，C++ unit testと `mmCal.Benchmarks --certification-boundaries` が主たるoracleである。black-box testでは値・保持式・precision metadataとして安定して観測できるものを固定する。

## 性能監査

black-boxのPASS/FAILとは独立に，遅い物理test fileを確認する場合は次を使用する。

```powershell
py tester.py --exe "..\..\build\x64\Release\mmCal.exe" --timings 10
```

`--timings`だけなら上位10ファイル，`--timings N`なら上位Nファイルについて，同一物理ファイル内の全isolated sessionのmmCal process wall timeを合算して表示する。これはoptimization対象を探すためのprofiling signalであり，test結果やtimeout判定には影響しない。

個別ファイルだけ実行する場合は通常どおり指定する。

```powershell
py tester.py --exe "..\..\build\x64\Release\mmCal.exe" test05_approximation_certification.txt
```

## 旧ファイルからの対応

| 新ファイル | 統合した旧ファイル |
|---|---|
| `test00_legacy_scalar.txt` | `test0_display_fix`，`test1_Syntax`，`test2_real`，`test3_complex` |
| `test01_legacy_numeric.txt` | `test4_statistics`，`test5_matrix` |
| `test02_legacy_errors.txt` | `test6_error` |
| `test03_language_runtime.txt` | `test0_v1_5_core`，`test10_angleMode`，`test14_numeric_lexer_boundaries`，`test17_history_relative`，`test19_iteration`，`test7_vars` |
| `test04_exact_arithmetic_algebra.txt` | `test11_exact_integer_rational`，`test12_exact_symbolic_complex`，`test20_number_theory`，`test21_algebraic_comparison`，`test22_number_field_interning`，`test25_algebraic_expression_bridge`，`test26_cyclotomic_fft` |
| `test05_approximation_certification.txt` | `test13_exact_approximation`，`test18_explain`，`test30_n_precision_audit`，`test31_n_refinement_audit`，`test32_real_complex_backend_audit` |
| `test06_calculus_solver.txt` | `test8_calculus`，`test16_exact_calculus_solver`，`test23_exponential_lambert_solve` |
| `test07_special_functions.txt` | `test9_special_func` |
| `test08_boundary_errors.txt` | `test15_boundary_errors_detailed` |
| `test09_semantics_audit.txt` | `test24_semantic_coherence`，`test27_audit`，`test28_semantic_hardening`，`test29_simplifier_domain_audit` |
