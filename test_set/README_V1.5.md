# mmCal V1.5 ブラックボックステスト

## 実行方法

`test_set` ディレクトリで実行します。

```powershell
py tester.py --exe "..\..\build\x64\Release\mmCal.exe"
```

実行ファイルを省略した場合は、`MMCAL_EXE` 環境変数と一般的な build 出力先を順に探索します。

各テストファイル先頭の `# @args ...` で、そのファイルだけに適用する起動引数を指定できます。旧テストの多くは当時の挙動を再現するため `--fix 15 --angle deg` を指定しています。


## 性能監査

black-boxのPASS/FAILとは独立に，遅いtest fileを確認する場合は次を使用する。

```powershell
py tester.py --exe "..\..\build\x64\Release\mmCal.exe" --timings 10
```

`--timings`だけなら上位10 file，`--timings N`なら上位N fileについて，mmCal processのwall timeと1 testあたり平均を表示する。これはoptimization対象を探すためのprofiling signalであり，test結果やtimeout判定には影響しない。

## 追加Solver回帰

`test23_exponential_lambert_solve.txt`は，`solve[equation,domain]`の一意変数推定，positive-real exponentialの非零証明，Lambert Wによる`a^x==x^2`のReal branch，protected solve変数の拒否に加え，Real `lambertw`のcertified `N`と`N[solve[...]]`によるsolution binding右辺の数値化を固定する。

`test24_semantic_coherence.txt`は，Stage 7の意味論接続を固定する。Solve-safe normalizationによる`E^x` / `exp[x]` / constant-base Power / `ln` / `log2` / `log10`の整合，Integer/Rational membershipの否定知識，free symbolを保持するstructural `N`，および値は存在するがcertified Complex backend未実装なLambert WをDomainErrorへ誤分類しないことを横断的に確認する。

`test25_algebraic_expression_bridge.txt`は，Stage 7-6のexact algebraic expression bridgeを固定する。canonical `root[...]`，`sqrt` / `cbrt`由来radical，`Phi`，およびそれらのbounded四則演算を同じ`AlgebraicNumber` backendへ接続し，表現差を跨ぐ`== != < <= > >=`，Integer/Rational membership，およびdirect Solve bindingのdomain filteringを確認する。

`test26_cyclotomic_fft.txt`は，Stage 7-7のexact Cyclotomic FFT backendを固定する。5/7/10/12点の非2冪exact round-tripを`Q[t]/Phi_n(t)`上で閉じ，Gaussian Rational入力では必要に応じ`Q(zeta_lcm(n,4))`へ埋め込み，既存2冪radix-2経路を変えないことを確認する。
