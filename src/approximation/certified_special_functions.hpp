#pragma once

#include "complex_interval.hpp"
#include "real_interval.hpp"
#include "numeric/rational.hpp"

#include <cstddef>

namespace mmcal::approximation {

// 実Gamma函数の包含区間。非正整数poleを含む入力区間は評価しない。
[[nodiscard]] RealInterval encloseGammaReal(
    const RealInterval& input,
    std::size_t precisionBits);

// exact Rational入力では元の有理値を失わず，argument shiftのrising factorialを
// exact productとして処理する。
[[nodiscard]] RealInterval encloseGammaRational(
    const numeric::Rational& input,
    std::size_t precisionBits);

// C/POSIXのlgammaと同じく、実軸上の log(|Gamma(x)|)。
[[nodiscard]] RealInterval encloseLogGammaReal(
    const RealInterval& input,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseLogGammaRational(
    const numeric::Rational& input,
    std::size_t precisionBits);

// Lambert W の実数branch。branch=0 は [-1/e,inf) 上のprincipal branch，
// branch=-1 は [-1/e,0) 上のlower real branchを単調性付きで直接囲う。
[[nodiscard]] RealInterval encloseLambertWReal(
    const RealInterval& input,
    int branch,
    std::size_t precisionBits);

// 複素branch。principal branchの小円板ではMaclaurin級数，その他は
// Log_k(z)-Log(w)の縮小写像を証明できる領域で任意整数branchを囲う。
[[nodiscard]] ComplexInterval encloseLambertWComplex(
    const ComplexInterval& input,
    const numeric::BigInt& branch,
    std::size_t precisionBits);

// -1/eからのoffsetをexact構造のまま受け取る局所backend。
// z=-1/e+offset として平方根branchをcertifyし，極端に小さいoffsetで
// e*z+1 の数値的cancellationにより情報を失うことを避ける。
[[nodiscard]] ComplexInterval encloseLambertWBranchPointOffset(
    const ComplexInterval& offset,
    const numeric::BigInt& branch,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseErfReal(
    const RealInterval& input,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseErfcReal(
    const RealInterval& input,
    std::size_t precisionBits);

// 複素特殊函数のcertified backend。principal branchを持つ函数は
// ComplexInterval上でbranch cutを跨がないことを証明できる場合だけ評価する。
[[nodiscard]] ComplexInterval encloseErfComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseErfcComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseGammaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);

// Fresnel C(x)=integral_0^x cos(Pi t^2/2)dt, S(x)=integral_0^x sin(Pi t^2/2)dt。
// 実軸ではTaylor/部分積分漸近展開を切り替え、剰余を明示的に包含する。
[[nodiscard]] RealInterval encloseFresnelCReal(
    const RealInterval& input,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseFresnelSReal(
    const RealInterval& input,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval encloseFresnelCComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval encloseFresnelSComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);

// Kummerの合流型超幾何函数 M(a,b,z)=1F1(a;b;z)。
// 現backendはexact Rationalのpoint引数を級数＋厳密tail boundで保証評価する。
[[nodiscard]] RealInterval encloseHypergeometric1F1Real(
    const numeric::Rational& a,
    const numeric::Rational& b,
    const numeric::Rational& z,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval encloseHypergeometric1F1Complex(
    const ComplexInterval& a,
    const ComplexInterval& b,
    const ComplexInterval& z,
    std::size_t precisionBits);

// Gauss 2F1。現backendはprincipal branchのうち |z|<1 のexact Rational pointを
// Gauss級数＋厳密な幾何tail boundで保証評価する。
[[nodiscard]] RealInterval encloseHypergeometric2F1Real(
    const numeric::Rational& a,
    const numeric::Rational& b,
    const numeric::Rational& c,
    const numeric::Rational& z,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval encloseHypergeometric2F1Complex(
    const ComplexInterval& a,
    const ComplexInterval& b,
    const ComplexInterval& c,
    const ComplexInterval& z,
    std::size_t precisionBits);

// Legendre不完全楕円積分。amplitudeは常にRadian。
// 小振幅ではm/n級数をfast pathとし，実効tail収束率が遅い領域は
// Carlson symmetric forms RF/RD/RJの保証付きduplicationへ送る。
[[nodiscard]] RealInterval encloseEllipticFReal(
    const numeric::Rational& phi,
    const numeric::Rational& m,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseEllipticEReal(
    const numeric::Rational& phi,
    const numeric::Rational& m,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseEllipticPiReal(
    const numeric::Rational& n,
    const numeric::Rational& phi,
    const numeric::Rational& m,
    std::size_t precisionBits);

// 有限precisionの実数入力用。級数fast pathとCarlson backendのどちらでも
// InformationEnclosureを保持したままbranch/pole条件を証明する。
[[nodiscard]] RealInterval encloseEllipticFReal(
    const RealInterval& phi,
    const RealInterval& m,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseEllipticEReal(
    const RealInterval& phi,
    const RealInterval& m,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseEllipticPiReal(
    const RealInterval& n,
    const RealInterval& phi,
    const RealInterval& m,
    std::size_t precisionBits);

// amplitudeがexactなq*Piである場合の専用経路。period番号をRationalのまま決定し，
// Piの独立な区間近似同士の除算による境界曖昧性を避ける。
[[nodiscard]] RealInterval encloseEllipticFRealPiMultiple(
    const numeric::Rational& piCoefficient,
    const RealInterval& m,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseEllipticERealPiMultiple(
    const numeric::Rational& piCoefficient,
    const RealInterval& m,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseEllipticPiRealPiMultiple(
    const RealInterval& n,
    const numeric::Rational& piCoefficient,
    const RealInterval& m,
    std::size_t precisionBits);

// 古典的積分函数。実数backendではprincipal real valueを保証区間で返す。
[[nodiscard]] RealInterval encloseExponentialIntegralEiReal(
    const RealInterval& input,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseSineIntegralSiReal(
    const RealInterval& input,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseCosineIntegralCiPositive(
    const RealInterval& input,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseLogarithmicIntegralLiPositive(
    const RealInterval& input,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval encloseExponentialIntegralEiComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseSineIntegralSiComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseCosineIntegralCiComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseLogarithmicIntegralLiComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);

// 正整数sを保証評価する。|z|<1の級数に加え，Li_2のprincipal connection formulaと，
// z≈1の正実数ではmu=log(z)整数極限展開を使う。有限precision実区間も幅を保って伝播する。
[[nodiscard]] RealInterval enclosePolylogReal(
    std::uint64_t order,
    const numeric::Rational& z,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval enclosePolylogComplex(
    std::uint64_t order,
    const ComplexInterval& z,
    std::size_t precisionBits);

// a,b>0 に対するBetaとlog Beta。Gammaの比ではなくlog-domainで評価する。
[[nodiscard]] RealInterval encloseBetaPositive(
    const RealInterval& a,
    const RealInterval& b,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseBetaRational(
    const numeric::Rational& a,
    const numeric::Rational& b,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseBetaLogPositive(
    const RealInterval& a,
    const RealInterval& b,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseBetaLogRational(
    const numeric::Rational& a,
    const numeric::Rational& b,
    std::size_t precisionBits);

// Riemann zeta。exact Rational s>1は元の有理値を保ったままEuler-Maclaurinへ送る。
[[nodiscard]] RealInterval encloseZetaRational(
    const numeric::Rational& input,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseZetaReal(
    const RealInterval& input,
    std::size_t precisionBits);

[[nodiscard]] ComplexInterval encloseZetaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);

// digamma/trigammaのexact Rational fast path。負の非整数はexact recurrenceで正側へ移す。
[[nodiscard]] RealInterval encloseDigammaRational(
    const numeric::Rational& input,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseTrigammaRational(
    const numeric::Rational& input,
    std::size_t precisionBits);

// digammaの実数backend。poleを跨がない負実数区間はrecurrenceで正実軸へ移す。
[[nodiscard]] RealInterval encloseDigammaReal(
    const RealInterval& input,
    std::size_t precisionBits);

// 正実数用の直接backend。
[[nodiscard]] RealInterval encloseDigammaPositive(
    const RealInterval& input,
    std::size_t precisionBits);
[[nodiscard]] RealInterval encloseTrigammaPositive(
    const RealInterval& input,
    std::size_t precisionBits);

// principal branchの複素digamma/trigamma。poleを含まず，右半平面への
// recurrence shiftと漸近剰余を保証できる区間だけ評価する。
[[nodiscard]] ComplexInterval encloseDigammaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);
[[nodiscard]] ComplexInterval encloseTrigammaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits);

// 正則化不完全Beta I_x(a,b)。a,bはexact positive Rational、xは[0,1]のcertified interval。
[[nodiscard]] RealInterval encloseIncompleteBetaRegularized(
    const RealInterval& a,
    const RealInterval& b,
    const RealInterval& x,
    std::size_t precisionBits);

[[nodiscard]] RealInterval encloseIncompleteBetaRegularized(
    const numeric::Rational& a,
    const numeric::Rational& b,
    const RealInterval& x,
    std::size_t precisionBits);

} // namespace mmcal::approximation
