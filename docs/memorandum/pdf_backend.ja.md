# PDF backend memorandum

mmCal の初版 PDF backend は **PDF 1.4** を出力する。Plot 側は PDF を意識せず，`GraphicsScene` から backend 共通 dispatch を経由して PDF へ lower する。

## 初版の構造

```text
%PDF-1.4
1 0 obj  Catalog
2 0 obj  Pages
3 0 obj  Page
4 0 obj  Contents stream
5 0 obj  Helvetica
6 0 obj  Times-Roman
7 0 obj  Courier
8 0 obj  XMP Metadata stream
9 0 obj  Info Dictionary
xref
trailer
startxref
%%EOF
```

初版はデバッグ性を優先し，object stream，xref stream，content/XMP stream の圧縮を使わない。従来型の classic xref table を出力するため，PDF をテキストエディタで開いて content stream と metadata を直接追える。

## 座標

`GraphicsScene` は mm・y 上向き。PDF も y 上向きなので，content stream 冒頭で

```text
2.83464566929134 0 0 2.83464566929134 0 0 cm
```

として `1 mm = 72 / 25.4 pt` の CTM を設定する。以降の path 座標・線幅・文字寸法は mm の値をそのまま記述する。

canvas 外へ伸びる漸近線等は

```text
0 0 width height re W n
```

で明示的に clip する。

## path

- `GraphicsMoveTo` -> `m`
- `GraphicsLineTo` -> `l`
- `GraphicsCubicTo` -> `c`
- `GraphicsClosePath` -> `h`
- stroke -> `S`
- fill -> `f`
- fill + stroke -> `B`

PDF には quadratic Bézier operator がないため，`GraphicsQuadraticTo` は cubic へ exact degree elevation して `c` を出す。数学曲線の形状は変えない。

円 marker は PDF に circle primitive がないため，4 本の cubic Bézier で描画する。これは marker 描画用であり Plot の数学曲線特殊化とは独立である。

## metadata

version はルート `version.h` の `MMCAL_VERSION_STRING` から取得する。

Info Dictionary:

- Title: `mmCal [version] Plot`
- Creator: `mmCal [version]`
- Producer: `mmCal [version] PDF Plotter`
- CreationDate: auto
- ModDate: auto

XMP:

- `dc:title`: `mmCal [version] Plot`
- `dc:language`: `ja-JP`
- `pdf:Producer`: `mmCal [version] PDF Plotter`
- `xmp:CreatorTool`: `mmCal [version]`
- `xmp:CreateDate`, `xmp:ModifyDate`, `xmp:MetadataDate`: auto
- `xmpMM:DocumentID`, `xmpMM:InstanceID`: auto UUID
- `mmcal:Version`: `[version]`
- `mmcal:UserName`: OS user name

`Author`, `Subject`, `Keywords`, copyright/rights は初版では出力しない。PC 名，絶対ファイルパス，git commit/build number も記録しない。

`PdfRenderOptions::deterministic=true` では日時と ID を固定して byte-for-byte 回帰を可能にする。username はテスト用に override 可能。

## Export filename

明示 format を省略した場合は拡張子で推定する。

```text
Export[plot[...], "a.svg"] -> SVG
Export[plot[...], "a.eps"] -> EPS
Export[plot[...], "a.pdf"] -> PDF
```

format を明示し，渡された suffix が異なる場合は置換せず canonical suffix を追加する。

```text
Export[plot[...], "a.svg", "PDF"] -> "a.svg.pdf"
Export[plot[...], "a.dat", "EPS"] -> "a.dat.eps"
Export[plot[...], "a.pdf", "PDF"] -> "a.pdf"
```

Export の返値は実際に書き込んだ path とする。

## 初版の意図的制限

- Base 14 font による ASCII text のみ。日本語/Unicode text は font embedding と ToUnicode 実装まで明示的 unsupported。
- alpha transparency は未実装。PDF 1.4 自体は transparency を持つが，silent drop はせず unsupported とする。
- encryption, linearization, object stream, xref stream, stream compression は未使用。
