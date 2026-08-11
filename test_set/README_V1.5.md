# mmCal V1.5 ブラックボックステスト

## 実行方法

`test_set` ディレクトリで実行します。

```powershell
py testerREPL.py --exe "..\..\build\x64\Release\mmCal.exe"
```

実行ファイルを省略した場合は、`MMCAL_EXE` 環境変数と一般的な build 出力先を順に探索します。

各テストファイル先頭の `# @args ...` で、そのファイルだけに適用する起動引数を指定できます。旧テストの多くは当時の挙動を再現するため `--fix 15 --angle deg` を指定しています。

