# mmCal EPS backendメモ

EPS（Encapsulated PostScript）は「1枚の図を他文書へ埋め込むためのPostScript」であり，SVGのようなXML木ではない。基本的には，ヘッダに図の外接矩形を宣言し，その後ろへPostScriptの描画命令を順番に並べる。

mmCalではPlot固有情報をEPSへ直接渡さず，次の共通経路を使う。

```text
Plot -> PlotScene -> GraphicsScene(mm, y-up)
                         |-> SVG backend
                         `-> EPS backend
```

`GraphicsScene`はbackend非依存で，座標，線幅，文字寸法をmmで保持する。EPS/PostScriptの座標はy上向きなので，SVGのようなy反転は不要である。

## 最小構造

mmCalが生成するEPSは概ね次の形になる。

```postscript
%!PS-Adobe-3.0 EPSF-3.0
%%BoundingBox: 0 0 426 284
%%HiResBoundingBox: 0 0 425.196850... 283.464566...
%%Creator: mmCal
%%LanguageLevel: 2
%%Pages: 1
%%EndComments

gsave
2.834645669... 2.834645669... scale

% canvas clip
newpath
0 0 moveto
150 0 lineto
150 100 lineto
0 100 lineto
closepath
clip
newpath

% drawing commands ...

grestore
%%EOF
```

`BoundingBox`はPostScript point単位で整数，`HiResBoundingBox`は実数で記録する。150 x 100 mmなら，

```text
150 mm = 425.196850... pt
100 mm = 283.464566... pt
```

となる。`BoundingBox`は図全体を確実に含むよう上端をceilして426 x 284 ptとする。

## mm座標

PostScriptの既定単位は1/72 inchなので，冒頭で

```postscript
2.834645669... 2.834645669... scale
```

を実行する。以後は`1 user unit = 1 mm`として扱える。したがってGraphicsSceneの座標や0.8 mmのcurve stroke等をそのまま出力できる。

## path命令

主な対応は次の通り。

```text
GraphicsMoveTo       -> moveto
GraphicsLineTo       -> lineto
GraphicsCubicTo      -> curveto
GraphicsClosePath    -> closepath
```

PostScriptにはquadratic Bezier命令が無い。`GraphicsQuadraticTo`は近似せず，厳密なdegree elevationでcubicへ変換する。

始点をP0，quadratic controlをQ，終点をP2とすると，cubic controlは

```text
C1 = P0 + 2/3 (Q - P0)
C2 = P2 + 2/3 (Q - P2)
```

であり，`P0, C1, C2, P2`のcubic Bezierは元のquadraticと完全に同じ曲線である。

## clip

SVGはroot viewportが暗黙のclipとして働くが，EPSにはそのようなcanvas viewportが自動では存在しない。そこで0..width x 0..height mmの矩形を明示的なclipping pathにする。

これにより`tan[x]`等でviewport外へ伸びる幾何がEPSを埋め込んだ先へ漏れない。

## stroke / fill

線幅，line cap，line join，RGB colorをPostScript graphics stateへ設定してから`stroke`する。circleは`arc`で構成する。

現EPS backendは透明度を黙って捨てない。alpha != 255のGraphicsSceneは`UnsupportedTransparency`として失敗させる。EPS/PostScript Level 2にはSVG/PDFのような通常のalpha transparencyが無いためである。

## text

初版はPlotの数値tick labelを対象とし，ASCII textをHelvetica/Times-Roman/Courierへ対応させる。文字列中の`(`, `)`, `\\`はPostScript stringとしてescapeする。

日本語・Unicode・font embeddingは初版EPSの対象外である。非ASCII文字を文字化けさせて出力せず，unsupportedとして扱う。

## semantic情報

SVGの`data-mmcal-kind`に相当する標準的なEPS属性は無いため，現在は

```postscript
% mmCal-semantic: curve 1
```

のようなコメントを描画命令の前へ残す。これは描画には影響しない。

## Export

```text
Export[plot[sin[x],{x,-Pi,Pi}], "plot.svg"]
Export[plot[sin[x],{x,-Pi,Pi}], "plot.eps"]
Export[plot[sin[x],{x,-Pi,Pi}], "out.dat", "EPS"]
```

拡張子推定は`.svg`と`.eps`を扱い，明示formatがある場合はそちらを優先する。
