## 目安

- 初等〜高校: I1, I6, D1, S1, S3の一部, S5の基本
- 大学基礎: I2–I4, I7, D2–D3, S2–S5, S7
- 大学上級: I5, I8, D4, S6, S8–S9
- 研究・CAS境界: R1–R2, X1 と prospective special-function 節

## I1 初等・基礎不定積分

- 段階: 初等〜高校・大学初年級
- 主分類: `E`
- 狙い: べき・指数・対数・三角・逆三角・双曲線の基本原始函数。

```text
  integrate[1,x]
  integrate[x,x]
  integrate[x^2,x]
  integrate[x^5,x]
  integrate[x^20,x]
  integrate[1/x,x]
  integrate[x^(-2),x]
  integrate[x^(1/2),x]
  integrate[x^(1/3),x]
  integrate[1/sqrt[x],x]
  integrate[1/cbrt[x],x]
  integrate[2x+3,x]
  integrate[x^3-4x+7,x]
  integrate[(2x+3)^5,x]
  integrate[(3x-1)^(-2),x]
  integrate[1/(2x+3),x]
  integrate[(1+x)^(-1/2),x]
  integrate[sqrt[1+x],x]
  integrate[cbrt[1+2x],x]
  integrate[exp[x],x]
  integrate[exp[2x+1],x]
  integrate[expm1[x],x]
  integrate[sin[x],x]
  integrate[cos[x],x]
  integrate[sin[2x],x]
  integrate[cos[3x+1],x]
  integrate[tan[x],x]
  integrate[cot[x],x]
  integrate[sec[x],x]
  integrate[csc[x],x]
  integrate[sinh[x],x]
  integrate[cosh[x],x]
  integrate[tanh[x],x]
  integrate[coth[x],x]
  integrate[sech[x],x]
  integrate[csch[x],x]
  integrate[log[x],x]
  integrate[log1p[x],x]
  integrate[asin[x],x]
  integrate[acos[x],x]
  integrate[atan[x],x]
  integrate[asinh[x],x]
  integrate[acosh[x],x]
  integrate[atanh[x],x]
  integrate[erf[x],x]
  integrate[erfc[x],x]
```

## I2 有理函数・部分分数・逆三角

- 段階: 大学初年級
- 主分類: `E/C`
- 狙い: 部分分数，既約二次因子，平方根型。hole と branch の監査にも使う。

```text
  integrate[1/(1+x),x]
  integrate[1/(1-x),x]
  integrate[1/(1+x^2),x]
  integrate[1/(1-x^2),x]
  integrate[1/(x^2-1),x]
  integrate[x/(1+x^2),x]
  integrate[x/(1+x^4),x]
  integrate[(x+1)/(x+2),x]
  integrate[(x+1)/(x^2+4),x]
  integrate[(2x+1)/(x^2+x+1),x]
  integrate[1/(x^2+2x+5),x]
  integrate[1/(x^2-2x+5),x]
  integrate[1/(x*(x+1)),x]
  integrate[1/(x*(x-1)),x]
  integrate[1/(x^2*(x+1)),x]
  integrate[1/((x-1)^2*(x+2)),x]
  integrate[(x^2+1)/(x^3-x),x]
  integrate[(2x^3+3x^2+1)/(x^2+1),x]
  integrate[1/(x^3+1),x]
  integrate[1/(x^4-1),x]
  integrate[1/(x^4+1),x]
  integrate[x^2/(x^4+1),x]
  integrate[(x^3+1)/(x^4+x^2+1),x]
  integrate[1/sqrt[4-x^2],x]
  integrate[1/sqrt[x^2+4],x]
  integrate[sqrt[4-x^2],x]
  integrate[sqrt[x^2+4],x]
  integrate[1/sqrt[x^2-1],x]
  integrate[x/sqrt[1-x^2],x]
  integrate[x/sqrt[1+x^2],x]
  integrate[1/(x*sqrt[x^2-1]),x]
```

## I3 三角・双曲線・reduction

- 段階: 大学初年級〜中級
- 主分類: `E/C`
- 狙い: 冪のreduction，Weierstrass置換，極を含む函数。

```text
  integrate[sin[x]^2,x]
  integrate[cos[x]^2,x]
  integrate[sin[x]^3,x]
  integrate[cos[x]^3,x]
  integrate[sin[x]^4,x]
  integrate[cos[x]^4,x]
  integrate[sin[x]^2*cos[x]^2,x]
  integrate[sin[2x]^(-2),x]
  integrate[cos[3x]^(-2),x]
  integrate[sec[x]^2,x]
  integrate[csc[x]^2,x]
  integrate[sec[x]^3,x]
  integrate[csc[x]^3,x]
  integrate[tan[x]^2,x]
  integrate[cot[x]^2,x]
  integrate[tan[x]^5,x]
  integrate[cot[x]^5,x]
  integrate[sec[x]^4,x]
  integrate[csc[x]^4,x]
  integrate[sin[x]*cos[x],x]
  integrate[sin[x]*cos[x]^3,x]
  integrate[sin[x]^3*cos[x]^2,x]
  integrate[1/(1+sin[x]),x]
  integrate[1/(1+cos[x]),x]
  integrate[1/(2+sin[x]),x]
  integrate[1/(2+cos[x]),x]
  integrate[1/(3+2sin[x]),x]
  integrate[1/(1+sin[x]+cos[x]),x]
  integrate[sinh[x]^2,x]
  integrate[cosh[x]^2,x]
  integrate[tanh[x]^2,x]
  integrate[sech[x]^2,x]
  integrate[csch[x]^2,x]
  integrate[sinh[x]*cosh[x],x]
  integrate[1/(1+cosh[x]),x]
```

## I4 置換・chain rule・部分積分

- 段階: 大学初年級〜中級
- 主分類: `E/S`
- 狙い: 逆chain rule，部分積分，多項式×超越函数。

```text
  integrate[2x*(1+x^2)^5,x]
  integrate[3x^2*exp[x^3],x]
  integrate[cos[x]*exp[sin[x]],x]
  integrate[sin[x]/(1+cos[x]),x]
  integrate[2x/(1+x^2),x]
  integrate[log[x]/x,x]
  integrate[1/(x*log[x]),x]
  integrate[log[1+x^2]/x,x]
  integrate[x*log[1+x^2],x]
  integrate[x*exp[x],x]
  integrate[x^2*exp[x],x]
  integrate[x^5*exp[x],x]
  integrate[x*sin[x],x]
  integrate[x*cos[x],x]
  integrate[x^2*sin[x],x]
  integrate[x^2*cos[x],x]
  integrate[x*sinh[x],x]
  integrate[x*cosh[x],x]
  integrate[x^2*log[x],x]
  integrate[x^7*log[x],x]
  integrate[exp[x]*sin[x],x]
  integrate[exp[x]*cos[x],x]
  integrate[exp[2x]*sin[3x],x]
  integrate[exp[2x]*cos[3x],x]
  integrate[E^x*cos[x],x]
  integrate[sqrt[x+sqrt[x]],x]
  integrate[1/(sqrt[x]*(1+sqrt[x])),x]
  integrate[sqrt[x]/(1+sqrt[x]),x]
```

## I5 Fresnel・Ei・Si・Ci・polylog・hypergeometric・elliptic

- 段階: 大学上級〜特殊函数
- 主分類: `S/U`
- 狙い: 特殊函数への還元と，特殊函数をさらに積分する能力境界。

```text
  integrate[exp[x]/x,x]
  integrate[sin[x]/x,x]
  integrate[cos[x]/x,x]
  integrate[1/log[x],x]
  integrate[li[x],x]
  integrate[log[1-x]/x,x]
  integrate[log[1+x]/x,x]
  integrate[log[1+x^2]/x,x]
  integrate[exp[-x^2],x]
  integrate[exp[x^2],x]
  integrate[exp[x^3],x]
  integrate[exp[x^4],x]
  integrate[exp[x^6],x]
  integrate[cos[Pi*x^2/2],x]
  integrate[sin[Pi*x^2/2],x]
  integrate[cos[4x^2],x]
  integrate[sin[2x^2],x]
  integrate[sin[2x^2]^4,x]
  integrate[sqrt[1+2x^3],x]
  integrate[1/(1+x^5),x]
  integrate[1/sqrt[1-x^4],x]
  integrate[1/sqrt[1-(1/3)*sin[x]^2],x]
  integrate[sqrt[1-(1/3)*sin[x]^2],x]
  integrate[1/((1-(1/5)*sin[x]^2)*sqrt[1-(1/3)*sin[x]^2]),x]
  integrate[hypergeometric1F1[a,b,x],x]
  integrate[hypergeometric2F1[a,b,c,x],x]
  integrate[polylog[2,x],x]
  integrate[polylog[3,x]/x,x]
  integrate[zeta[x],x]
  integrate[gamma[x],x]
  integrate[x^x,x]
```

## I6 定積分・基本

- 段階: 高校〜大学初年級
- 主分類: `E/S`
- 狙い: 有限区間，対称性，endpoint substitution，特殊値。

```text
  integrate[x^2,{x,0,1}]
  integrate[x^4+x^2+1,{x,1,3}]
  integrate[sin[x],{x,0,Pi}]
  integrate[cos[x],{x,0,Pi/2}]
  integrate[1/x,{x,1,2}]
  integrate[log[x],{x,1,Pi}]
  integrate[1/(1+x^2),{x,0,1}]
  integrate[1/(1+x^2),{x,-1,1}]
  integrate[1/(x^3+1),{x,0,1}]
  integrate[sqrt[1-x^2],{x,-1,1}]
  integrate[sqrt[1-x^2],{x,0,1}]
  integrate[1/sqrt[1-x^2],{x,-1,1}]
  integrate[sin[x]^2,{x,0,Pi}]
  integrate[cos[x]^4,{x,0,Pi/2}]
  integrate[sin[x]*cos[x],{x,0,Pi/2}]
  integrate[exp[x],{x,0,1}]
  integrate[exp[-x],{x,0,1}]
  integrate[erf[x],{x,0,1}]
  integrate[fresnelc[x],{x,0,1}]
  integrate[log[1+x]/x,{x,0,1}]
  integrate[log[1-x]/x,{x,0,1}]
  integrate[log[x]^2,{x,0,1}]
  integrate[x^2*log[x],{x,0,1}]
```

## I7 improper・無限区間・古典積分

- 段階: 大学中級〜解析
- 主分類: `S/C`
- 狙い: improper integral，条件収束，内部特異点。主値を暗黙にしない。

```text
  integrate[exp[-x],{x,0,Infinity}]
  integrate[1/x^2,{x,1,Infinity}]
  integrate[1/(1+x^2),{x,-Infinity,Infinity}]
  integrate[1/sqrt[x],{x,0,1}]
  integrate[log[x],{x,0,1}]
  integrate[exp[-x^2],{x,0,Infinity}]
  integrate[exp[-x^2],{x,-Infinity,Infinity}]
  integrate[exp[-x^4],{x,0,Infinity}]
  integrate[exp[-x^4],{x,-Infinity,Infinity}]
  integrate[sin[x]/x,{x,0,Infinity}]
  integrate[cos[x]/x,{x,1,Infinity}]
  integrate[sin[x]^2/x^2,{x,0,Infinity}]
  integrate[cos[Pi*x^2/2],{x,0,Infinity}]
  integrate[sin[Pi*x^2/2],{x,0,Infinity}]
  integrate[1/(1+x^4),{x,0,Infinity}]
  integrate[x/(1+x^4),{x,0,Infinity}]
  integrate[1/sqrt[1-x^4],{x,0,1}]
  integrate[log[sin[x]],{x,0,Pi/2}]
  integrate[log[cos[x]],{x,0,Pi/2}]
  integrate[atan[x]/x,{x,0,1}]
  integrate[log[x]/(1+x^2),{x,0,1}]
  integrate[1/(x-2),{x,1,Infinity}]
  integrate[1/x,{x,-1,1}]
  integrate[tan[x],{x,0,2}]
```

## I8 パラメータ・assumptions

- 段階: 大学中級〜解析
- 主分類: `C/S`
- 狙い: パラメータ依存収束条件，Gamma/Beta/Mellin型，assumption伝播。

```text
  integrate[x^n,x,n!=-1]
  integrate[x^n,{x,0,1},n>0]
  integrate[x^n,{x,0,a},{a>0,n>-1}]
  integrate[exp[-a*x],{x,0,Infinity},a>0]
  integrate[exp[-a*x^2],{x,-Infinity,Infinity},a>0]
  integrate[exp[-a*x^2],{x,0,Infinity},a>0]
  integrate[x*exp[-a*x^2],{x,0,Infinity},a>0]
  integrate[x^(s-1)*exp[-x],{x,0,Infinity},s>0]
  integrate[x^(a-1)*(1-x)^(b-1),{x,0,1},{a>0,b>0}]
  integrate[x^(a-1)/(1+x),{x,0,Infinity},{a>0,a<1}]
  integrate[x^(s-1)/(exp[x]-1),{x,0,Infinity},s>1]
  integrate[x^(s-1)/(exp[x]+1),{x,0,Infinity},s>0]
  integrate[(exp[-a*x]-exp[-b*x])/x,{x,0,Infinity},{a>0,b>0}]
  integrate[exp[-a*x]*cos[b*x],{x,0,Infinity},a>0]
  integrate[exp[-a*x]*sin[b*x],{x,0,Infinity},a>0]
  integrate[exp[-a*x^2]*cos[b*x],{x,-Infinity,Infinity},a>0]
  integrate[1/(x^2+a^2),{x,-Infinity,Infinity},a>0]
  integrate[1/sqrt[a^2-x^2],{x,-a,a},a>0]
  integrate[sqrt[a^2-x^2],{x,-a,a},a>0]
  integrate[abs[x],x,x>=0]
  integrate[abs[x],x,x<=0]
  integrate[sqrt[x^2],x,x>=0]
  integrate[sqrt[x^2],x,x<=0]
  integrate[1/(x-a),{x,0,1},{a<0}]
  integrate[1/(x-a),{x,0,1},{a>1}]
  integrate[1/(x-a),{x,0,1},{a>0,a<1}]
```

## D1 基礎微分

- 段階: 高校〜大学初年級
- 主分類: `E`
- 狙い: 基本則，合成函数，商，冪。

```text
  D[x,x]
  D[x^2,x]
  D[x^10,x]
  D[1/x,x]
  D[sqrt[x],x]
  D[cbrt[x],x]
  D[exp[x],x]
  D[log[x],x]
  D[sin[x],x]
  D[cos[x],x]
  D[tan[x],x]
  D[sinh[x],x]
  D[cosh[x],x]
  D[tanh[x],x]
  D[asin[x],x]
  D[acos[x],x]
  D[atan[x],x]
  D[asinh[x],x]
  D[acosh[x],x]
  D[atanh[x],x]
  D[x^x,x]
  D[exp[x^2],x]
  D[log[1+x^2],x]
  D[sin[x^2],x]
  D[sqrt[1+x^4],x]
  D[(x^2+1)/(x^2-1),x]
  D[(sin[x]+cos[x])^5,x]
```

## D2 高階微分

- 段階: 大学初年級〜上級
- 主分類: `E/S/U`
- 狙い: 高階微分，特殊函数，高次数での式膨張。

```text
  D[x^20,{x,5}]
  D[exp[x],{x,20}]
  D[sin[x],{x,4}]
  D[cos[x],{x,17}]
  D[exp[2x],{x,10}]
  D[1/x,{x,5}]
  D[log[x],{x,6}]
  D[x^x,{x,2}]
  D[exp[-x^2],{x,4}]
  D[sin[x^2],{x,3}]
  D[1/(1+x^2),{x,6}]
  D[lambertw[x],{x,2}]
  D[gamma[x],{x,2}]
  D[polylog[3,x],{x,3}]
```

## D3 多変数・偏微分

- 段階: 大学初年級〜中級
- 主分類: `E/S`
- 狙い: 偏微分，混合微分，パラメータ扱い。

```text
  D[x^2*y^3,x]
  D[x^2*y^3,y]
  D[x^2*y^3,x,y]
  D[x^2*y^3,x,y,y]
  D[exp[x*y],x]
  D[exp[x*y],x,y]
  D[sin[x*y+z],x,y,z]
  D[log[x^2+y^2],x]
  D[log[x^2+y^2],x,y]
  D[(x^2+y^2+z^2)^3,x,y,z]
  D[exp[x*y]*sin[z*x],x,y]
  D[hypergeometric1F1[a,b,x*y],x,y]
  D[hypergeometric2F1[a,b,c,x*y],x,y]
  D[x^y,x]
  D[x^y,y]
  D[x^y,x,y]
```

## D4 特殊函数・branch・未解決候補

- 段階: 大学上級〜研究境界
- 主分類: `S/C/U`
- 狙い: 特殊函数の微分知識，non-holomorphic函数，積分とのFTC接続。

```text
  D[gamma[x],x]
  D[lgamma[x],x]
  D[digamma[x],x]
  D[trigamma[x],x]
  D[erf[x],x]
  D[erfc[x],x]
  D[fresnelc[x],x]
  D[fresnels[x],x]
  D[Ei[x],x]
  D[Si[x],x]
  D[Ci[x],x]
  D[li[x],x]
  D[polylog[2,x],x]
  D[polylog[s,x],x]
  D[hypergeometric1F1[a,b,x],x]
  D[hypergeometric2F1[a,b,c,x],x]
  D[ellipticF[x,m],x]
  D[ellipticE[x,m],x]
  D[ellipticPi[n,x,m],x]
  D[lambertw[x],x]
  D[lambertw[-1,x],x]
  D[zeta[x],x]
  D[abs[x],x]
  D[sign[x],x]
  D[re[x],x]
  D[im[x],x]
  D[conj[x],x]
  D[arg[x],x]
  D[sinc[x],x]
  D[cosc[x],x]
  D[cases[x^2 if x>=0; -x if x<0],x]
  D[integrate[t^2,{t,0,x}],x]
  D[integrate[exp[-t^2],{t,0,x}],x]
  D[integrate[f[t],{t,0,x}],x]
```

## S1 線形・二次・有理方程式

- 段階: 中学〜大学初年級
- 主分類: `E/C`
- 狙い: 一次・二次・有理・根号。分母holeやextraneous root監査を含む。

```text
  solve[x==0,x]
  solve[x+1==0,x]
  solve[2x+3==0,x]
  solve[a*x+b==0,x]
  solve[x^2==1,x]
  solve[x^2==2,x]
  solve[x^2+1==0,x]
  solve[x^2-2x+1==0,x]
  solve[x^2+3x+2==0,x]
  solve[a*x^2+b*x+c==0,x]
  solve[(x-1)*(x+2)==0,x]
  solve[(x^2-1)/(x-1)==0,x]
  solve[1/(x-1)==0,x]
  solve[(x+1)/(x-2)==0,x]
  solve[x+1/x==2,x]
  solve[x+1/x==3,x]
  solve[x^2+1/x^2==2,x]
  solve[sqrt[x]==2,x]
  solve[sqrt[x+1]==x-1,x]
  solve[cbrt[x+1]==x-1,x]
```

## S2 高次多項式・代数数

- 段階: 高校〜代数学
- 主分類: `E/S`
- 狙い: cubic/quartic，高次Root，次数budget境界。

```text
  solve[x^3-1==0,x]
  solve[x^3-x==0,x]
  solve[x^3-2==0,x]
  solve[x^4-1==0,x]
  solve[x^4+1==0,x]
  solve[x^4-10x^2+9==0,x]
  solve[x^5-x+1==0,x]
  solve[x^5-x-1==0,x]
  solve[x^7+x+1==0,x]
  solve[x^10-2==0,x]
  solve[x^16+x+1==0,x]
  solve[x^32-x+1==0,x]
  solve[x^64+x+1==0,x]
  solve[x^65+x+1==0,x]
  solve[(x^2-2)*(x^3-3)==0,x]
  solve[x^6-5x^4+6x^2-1==0,x]
  solve[x^8+4x^6+6x^4+4x^2+2==0,x]
```

## S3 定義域・constraint・不等式

- 段階: 高校〜実解析
- 主分類: `C`
- 狙い: Real/Complex/Rational/Integer，ordered inequality，追加constraint。

```text
  solve[x^2+1==0,x,Real]
  solve[x^2+1==0,x,Complex]
  solve[x^2==2,x,Rational]
  solve[x^2==4,x,Integer]
  solve[x^2<4,x]
  solve[x^2<=4,x]
  solve[x^2>4,x]
  solve[x^2-5x+6>=0,x]
  solve[(x-1)/(x+2)>0,x]
  solve[1/x>0,x]
  solve[abs[x]<2,x]
  solve[abs[x-1]>=3,x]
  solve[sqrt[x]==2,x,Real]
  solve[log[x]==0,x,Real]
  solve[x^2==2,x,x>0]
  solve[x^2==2,x,x<0]
  solve[x^3==x,x,x>=0]
  solve[x^2==2,Real]
  solve[x^2+1==0,Real]
  solve[x^2+1==0,Complex]
  solve[x^2==4,Integer]
```

## S4 指数・対数・Lambert W

- 段階: 高校〜特殊函数
- 主分類: `S/U`
- 狙い: global inverse，Lambert W，閉形式を持たない超越方程式。

```text
  solve[exp[x]==1,x,Real]
  solve[exp[x]==2,x,Real]
  solve[exp[x]==0,x,Real]
  solve[log[x]==2,x,Real]
  solve[log[x]==0,x,Real]
  solve[log[x]==-1,x,Real]
  solve[log[x]==I*Pi,x]
  solve[E^x==8,x,Real]
  solve[2^x==8,x,Real]
  solve[2^(2x+1)==8,x,Real]
  solve[2^x==-1,x,Real]
  solve[1.1^x==0,x,Real]
  solve[1.1^x==x^2,x,Real]
  solve[2^x==x,x,Real]
  solve[2^x==x^2,x,Real]
  solve[x*exp[x]==1,x,Real]
  solve[x*exp[x]==a,x,Real]
  solve[x+log[x]==0,x,Real]
  solve[exp[x]+x==0,x,Real]
  solve[x^x==2,x,Real]
  solve[x^x==1,x,Real]
  solve[log[x]==x,x,Real]
  solve[exp[x]==x,x,Real]
  solve[exp[-x]==x,x,Real]
```

## S5 三角・双曲線

- 段階: 高校〜解析
- 主分類: `C/U`
- 狙い: 周期解族，値域，非線形超越方程式。

```text
  solve[sin[x]==0,x,Real]
  solve[cos[x]==0,x,Real]
  solve[tan[x]==1,x,Real]
  solve[sin[2x+1]==0,x,Real]
  solve[cos[3x-2]==1/2,x,Real]
  solve[tan[5x]==-1,x,Real]
  solve[sin[x]==1,x,Real]
  solve[sin[x]==2,x,Real]
  solve[cos[x]==-1,x,Real]
  solve[tan[x]==0,x,Real]
  solve[sin[x]==x,x,Real]
  solve[cos[x]==x,x,Real]
  solve[tan[x]==x,x,Real]
  solve[sinh[x]==2,x,Real]
  solve[cosh[x]==2,x,Real]
  solve[tanh[x]==1/2,x,Real]
  solve[tanh[x]==2,x,Real]
  solve[asinh[x]==2,x,Real]
  solve[atanh[x]==1/2,x,Real]
  solve[sinh[3x]==2,x,Real]
```

## S6 特殊函数方程式・高難度

- 段階: 大学上級〜研究境界
- 主分類: `U/R`
- 狙い: 一般に逆函数を一個返せない特殊函数方程式。

```text
  solve[lambertw[x]==1,x]
  solve[erf[x]==0,x,Real]
  solve[erf[x]==1/2,x,Real]
  solve[gamma[x]==1,x,Real]
  solve[gamma[x]==2,x,Real]
  solve[digamma[x]==0,x,Real]
  solve[zeta[x]==0,x,Real]
  solve[zeta[x]==0,x,Complex]
  solve[zeta[1/2+I*x]==0,x,Real]
  solve[zeta[x]==1,x,Real]
  solve[polylog[2,x]==0,x,Real]
  solve[polylog[2,x]==1,x,Real]
  solve[Ei[x]==0,x,Real]
  solve[Si[x]==1,x,Real]
  solve[Ci[x]==0,x,Real]
  solve[li[x]==0,x,Real]
  solve[fresnelc[x]==1/2,x,Real]
  solve[fresnels[x]==1/2,x,Real]
  solve[hypergeometric1F1[a,b,x]==0,x]
  solve[hypergeometric2F1[a,b,c,x]==0,x]
  solve[ellipticF[x,m]==1,x,Real]
```

## S7 多変数線形・多項式系

- 段階: 高校〜代数幾何
- 主分類: `E/R`
- 狙い: 線形系，zero-dimensional polynomial system，Groebner/elimination対象。

```text
  solve[{x+y==3,x-y==1},{x,y}]
  solve[{2x+3y==5,x-2y==9},{x,y}]
  solve[{x+y+z==6,x-y==0,y-z==0},{x,y,z}]
  solve[{x+y==1,2x+2y==2},{x,y}]
  solve[{x+y==1,2x+2y==3},{x,y}]
  solve[{x^2+y^2==1,y==x},{x,y}]
  solve[{x^2+y^2==1,y==0},{x,y}]
  solve[{x*y==1,x+y==3},{x,y}]
  solve[{x^2+y^2==5,x*y==2},{x,y}]
  solve[{x^2-y==0,y^2-x==0},{x,y}]
  solve[{x*y-1==0,y^2-x==0},{x,y}]
  solve[{x^2+y^2-1==0,x^2-y==0},{x,y}]
  solve[{x^3-3x^2-y+1==0,-x^2+y^2-1==0},{x,y}]
  solve[{x^3+y^3==2,x+y==2},{x,y}]
  solve[{x^4+y^4==1,x^2+y^2==1},{x,y}]
  solve[{x*y*z==1,x+y+z==3,x*y+y*z+z*x==3},{x,y,z}]
  solve[{x^2+y^2+z^2==1,x+y+z==0,x-y==0},{x,y,z}]
  solve[{x^2-y==0,x*y-z==0,z^2-x==0},{x,y,z}]
```

## S8 positive-dimensional・パラメータ・条件分岐

- 段階: 大学〜代数幾何
- 主分類: `C/U`
- 狙い: 無限解集合，パラメータ例外条件，conditional solution。

```text
  solve[{x+y==1,2x+2y==2},{x,y}]
  solve[{x*y==0},{x,y}]
  solve[{x^2+y^2==1},{x,y}]
  solve[{x*y==1},{x,y}]
  solve[{x^2-y^2==0},{x,y}]
  solve[a*x==1,x]
  solve[a*x==0,x]
  solve[x^2==a,x]
  solve[x^2+a*x+1==0,x]
  solve[exp[x]==a,x,Real]
  solve[tanh[x]==a,x,Real]
  solve[sin[x]==a,x,Real]
  solve[log[x]==a,x,Real]
```

## S9 整数・有理数定義域（単変数）

- 段階: 整数論
- 主分類: `C/U`
- 狙い: 定義域制約，離散解，現Solverの能力境界。

```text
  solve[x^2==2,x,Integer]
  solve[x^2==4,x,Integer]
  solve[x^3==8,x,Integer]
  solve[x^2-5x+6==0,x,Integer]
  solve[x^4-5x^2+4==0,x,Integer]
  solve[x^2==2,x,Rational]
  solve[3x==1,x,Rational]
  solve[x^3==2,x,Rational]
  solve[x^5-x==0,x,Integer]
  solve[2^x==8,x,Integer]
  solve[2^x==3,x,Integer]
  solve[fib[x]==55,x,Integer]
  solve[fact[x]==120,x,Integer]
```

## R1 特異点・branch・不連続・能力境界

- 段階: 解析・CAS意味論
- 主分類: `C/U`
- 狙い: 数学的に値がない／branch条件が不足／大域簡約が危険な例。

```text
  integrate[1/x,{x,-1,1}]
  integrate[1/(x-1)^2,{x,0,2}]
  integrate[log[x],{x,-1,1}]
  integrate[sqrt[x],{x,-1,1}]
  integrate[1/sqrt[x],{x,-1,1}]
  integrate[1/(x*log[x]),{x,1/2,2}]
  integrate[1/(1-x^2),{x,-2,2}]
  integrate[atanh[x],{x,-2,2}]
  integrate[abs[x],x]
  integrate[sign[x],x]
  integrate[arg[x],x]
  D[abs[x],x]
  D[conj[x],x]
  D[arg[x],x]
  solve[log[x]==log[-x],x]
  solve[sqrt[x^2]==x,x]
  solve[1/x==0,x]
  solve[(x^2-1)/(x-1)==x+1,x]
```

## R2 研究・特殊函数・解析数論方向

- 段階: 研究・計算数学
- 主分類: `R/U`
- 狙い: 解析数論，Abelian/hyperelliptic積分，特殊函数零点，高難度多項式系。

```text
  integrate[x^(s-1)/(exp[x]-1),{x,0,Infinity},s>1]
  integrate[log[x]^2/(1-x),{x,0,1}]
  integrate[log[x]^3/(1-x),{x,0,1}]
  integrate[sin[x]/(x*(exp[x]-1)),{x,0,Infinity}]
  integrate[cos[x^3],{x,0,Infinity}]
  integrate[sin[x^3],{x,0,Infinity}]
  integrate[cos[x^4],{x,0,Infinity}]
  integrate[1/sqrt[1-x^5],x]
  integrate[1/sqrt[1-x^7],x]
  integrate[1/sqrt[x^5-x+1],x]
  integrate[exp[-x^2]/(1+x^2),{x,0,Infinity}]
  integrate[log[1+x^2]/(1+x^2),{x,0,Infinity}]
  integrate[polylog[2,x]/x,{x,0,1}]
  integrate[polylog[3,x]/x,{x,0,1}]
  solve[zeta[1/2+I*x]==0,x,Real]
  solve[zeta[x]==0,x,Complex]
  solve[gamma[x]==x,x,Real]
  solve[lambertw[x]==x,x,Complex]
  solve[polylog[2,x]==x,x,Complex]
  solve[hypergeometric2F1[1/2,1/2,1,x]==0,x,Complex]
  solve[{x^5+y^5==1,x^2+y^2==1},{x,y}]
  solve[{x^3+y^3+z^3==33},{x,y,z}]
  solve[{x^4+y^4==z^4},{x,y,z}]
```

## X1 恒等・空集合・定義域hole・interface境界

- 段階: CAS意味論・negative test
- 主分類: `C/U`
- 狙い: Universal/Empty/Conditional，hole保存，symbolic derivative order，型・変数推定・予約symbol境界。

```text
  integrate[0,x]
  integrate[a,x]
  integrate[y,x]
  integrate[x*y,x]
  integrate[x^2+gamma[x],x]
  integrate[abs[x]+x^2,x]
  D[sin[x],{x,0}]
  D[x^3,{x,0}]
  D[x^3,{x,4}]
  D[sin[x],{x,n}]
  D[x^a,x]
  D[a^x,x]
  D[log[a,x],x]
  D[atan2[y,x],x]
  D[integrate[gamma[x],x],x]
  solve[0==0,x]
  solve[0==1,x]
  solve[x==x,x]
  solve[x!=x,x]
  solve[x^2>=0,x,Real]
  solve[x^2<0,x,Real]
  solve[1/x==1/x,x]
  solve[(x^2-1)/(x-1)==x+1,x]
  solve[x^2==1,x,x!=1]
  solve[x^2==a,x,a>0]
  solve[sin[x]>0,x,Real]
  solve[cos[x]>=0,x,Real]
  solve[tan[x]>0,x,Real]
  solve[x+y==1,Real]
  solve[2+2==4,Real]
  solve[x==1,Pi]
```

# 将来の特殊函数語彙を仮定した追加ターゲット

以下の head は現行 `reference.ja.md` の source-callable 一覧にない。よって **構文スタイル案／将来ターゲット** であり，現行登録名とはみなさない。Airy/Bessel/一般化超幾何函数を追加する段階で名前を確定してから正式 corpus へ昇格する。

## F1 Airy/Bessel/orthogonal-polynomial prospective heads

```text
integrate[cos[t^3/3+x*t],{t,0,Infinity}]
D[airyai[x],x]
D[airybi[x],x]
D[besselj[n,x],x]
D[bessely[n,x],x]
D[besseli[n,x],x]
D[besselk[n,x],x]
integrate[cos[z*cos[t]]*sin[t]^(2n+1),{t,0,Pi}]
solve[airyai[x]==0,x,Real]
solve[besselj[0,x]==0,x,Real]
solve[besselj[n,x]==0,x,Real]
```

## F2 generalized special-function / research prospective heads

```text
integrate[meijerg[a,b,c,d,x],x]
D[meijerg[a,b,c,d,x],x]
solve[meijerg[a,b,c,d,x]==0,x]
integrate[appellF1[a,b,bp,c,x,y],x]
D[appellF1[a,b,bp,c,x,y],x]
solve[appellF1[a,b,bp,c,x,y]==0,x]
```

# 推奨 expected-class

同じ入力群を将来 test case 化する際は，少なくとも次を別分類にする。

1. `ExactClosedForm` — exactな有限式へ閉じる。
2. `ConditionalClosedForm` — assumption/定義域 条件付きで閉じる。
3. `SpecialFunctionClosedForm` — 既知特殊函数へ閉じる。
4. `ParameterizedSolutionFamily` — 周期函数など整数パラメータ族。
5. `FiniteRootObject` — radicalでなく `root[...]` 等のexact object。
6. `UniversalSolutionSet` — 恒等的に成立。ただしdefinedness holeを落とさない。
7. `EmptySolutionSet` — 定義域上の真の空集合。
8. `UnresolvedSolutionSet` — 解の不存在を主張できず未解決。
9. `ConditionsRequired` — branch/定義域不足。
10. `NoKnownClosedForm` — 現標準函数語彙で有限閉形式を期待しない。
11. `DomainError` / `TypeError` — 数学的定義域外またはinterface contract違反。
