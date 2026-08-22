// 三角函数の保証付き評価
#include "certified_trigonometry.hpp"
#include "certified_precision.hpp"

#include "approximation_context.hpp"
#include "certification_error.hpp"
#include "certified_constants.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "mathematics/math_ids.hpp"
#include "mathematics/trigonometric_reduction.hpp"
#include "numeric/big_int.hpp"
#include "numeric/detail/binary_scale.hpp"
#include "numeric/rational_rounding.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace mmcal::approximation {
namespace {

using mathematics::FunctionId;
using mathematics::ReducedTrigAngle;
using numeric::BigFloat;
using numeric::BigInt;
using numeric::DecimalApproximation;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}



[[nodiscard]] std::size_t nextGuardDigits(std::size_t guardDigits) {
    // guard桁数そのものを正しさの根拠にはしない。
    // 最終の10進丸めが区間両端で一致しなければ、作業precisionを増やして再試行する。
    const std::size_t growth = std::max<std::size_t>(8, guardDigits / 2);
    return checkedPrecisionAdd(guardDigits, growth, "Certified trigonometric precision is too large");
}

[[nodiscard]] BigFloat absoluteValue(const BigFloat& value) {
    return value.isNegative() ? -value : value;
}

[[nodiscard]] BigFloat absoluteUpperBound(const RealInterval& value) {
    const BigFloat lowerMagnitude = absoluteValue(value.lower());
    const BigFloat upperMagnitude = absoluteValue(value.upper());
    return lowerMagnitude > upperMagnitude ? lowerMagnitude : upperMagnitude;
}

[[nodiscard]] RealInterval symmetricInterval(const BigFloat& radius) {
    if (radius.isNegative())
        throw std::invalid_argument("Symmetric interval radius cannot be negative");
    return RealInterval{-radius, radius};
}

[[nodiscard]] CertifiedTrigEnclosure encloseTaylorPoint(
    FunctionId function,
    const Rational& x,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Trigonometric precision must be at least one bit");
    if (function != FunctionId::Sin && function != FunctionId::Cos)
        throw std::invalid_argument("Taylor trigonometric backend supports sin/cos only");

    // -------------------------------------------------------------------------
    // なぜBigFloatでTaylorを回してもcertifiedなのか
    // -------------------------------------------------------------------------
    // 以前はTaylorの各項と部分和を巨大なRationalとして完全exactに保持していた。
    // 証明は非常に分かりやすい一方、1000桁級では分母・分子が急激に巨大化する。
    // ここでは同じ数学的証明を RealInterval + directed rounding へ移す。
    //
    // 各 `term` と `sum` は単一BigFloatではなく区間である。乗算・除算・加算のたびに
    // lowerを -infinity、upperを +infinity へ丸めるので、真のTaylor項と真の部分和は常にその区間内に残る。
    // したがって途中丸めが何回起きても「真値を落とさない」。
    //
    // sin:
    //   t_0 = x
    //   t_(k+1) = -t_k*x^2 / ((2k+2)(2k+3))
    //
    // cos:
    //   t_0 = 1
    //   t_(k+1) = -t_k*x^2 / ((2k+1)(2k+2))
    //
    // TaylorのLagrange剰余について、sin/cosの任意階導関数の絶対値は実数全体で<=1。
    // よって現在の部分和の次のexact項を T とすれば
    //
    //   |remainder| <= |T|
    //
    // が成立する。`nextTerm` interval はそのexact Tを必ず含むので、区間両端の絶対値の大きい方 B は必ず |T| 以上である。
    // 従って
    //
    //   sum + [-B, +B]
    //
    // は、BigFloat丸め誤差とTaylor打切り誤差の両方を含んだcertified enclosureになる。
    //
    // 「nextTermが表示上0になった」「前回値から変わらない」といったmachine-float的な停止条件は一切使わない。
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Trigonometric working precision is too large");
    const std::size_t thresholdBits = checkedPrecisionAdd(
        precisionBits, 16, "Trigonometric precision is too large");

    // xはRationalだが、RealInterval::fromRationalが外向きに囲うので、dyadicでなくても入力の真値は失われない。
    // 現在の呼出経路ではPi区間の中点などdyadicが主である。
    const RealInterval xInterval = RealInterval::fromRational(x, workBits);
    const RealInterval xSquared = multiply(xInterval, xInterval, workBits);
    const BigFloat threshold = RealInterval::fromRational(
        binaryPrecisionThreshold(thresholdBits), workBits).upper();

    RealInterval term = function == FunctionId::Sin
        ? xInterval
        : RealInterval::fromRational(rational(1), workBits);
    RealInterval sum = term;

    BigInt firstFactor = function == FunctionId::Sin ? BigInt{2} : BigInt{1};
    BigInt secondFactor = function == FunctionId::Sin ? BigInt{3} : BigInt{2};
    std::size_t termsUsed = 1;

    for (;;) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        const Rational denominatorValue{firstFactor * secondFactor};
        const RealInterval denominator = RealInterval::fromRational(
            denominatorValue, workBits);

        RealInterval nextTerm = divide(
            multiply(term, xSquared, workBits), denominator, workBits);
        nextTerm = negate(nextTerm);

        // exactな次項はnextTerm区間内にある。従ってこのBはTaylor剰余の安全な上界。
        const BigFloat remainderBound = absoluteUpperBound(nextTerm);
        if (remainderBound <= threshold) {
            RealInterval enclosure = add(
                sum, symmetricInterval(remainderBound), workBits);
            return CertifiedTrigEnclosure{
                enclosure.roundedOutward(precisionBits),
                termsUsed,
                precisionBits
            };
        }

        sum = add(sum, nextTerm, workBits);
        term = std::move(nextTerm);
        if (termsUsed != std::numeric_limits<std::size_t>::max())
            ++termsUsed;

        firstFactor += BigInt{2};
        secondFactor += BigInt{2};
    }
}

[[nodiscard]] CertifiedTrigEnclosure encloseFunctionAtRadianInterval(
    FunctionId function,
    const RealInterval& radians,
    std::size_t precisionBits) {
    // Degree/Grad -> radian変換ではPiがintervalなので、入力角度も [a,b] になる。
    // 同じ不確定xをinterval Taylorへ直接何度も掛けるとdependency problemで区間が膨らむ。
    // そこで中点cだけをTaylorで評価し、入力区間半径rの影響は解析的に別途足す。
    //
    //   |sin'(x)| = |cos(x)| <= 1
    //   |cos'(x)| = |sin(x)| <= 1
    //
    // Mean Value Theoremより sin/cosはいずれも実数全体で1-Lipschitz:
    //
    //   |f(x)-f(c)| <= |x-c|
    //
    // x in [a,b], c=(a+b)/2, r=(b-a)/2 なら、f(c)のcertified intervalを
    // [-r,+r]だけ拡張すれば、単調性や象限境界を仮定せずf([a,b])を必ず包含できる。
    const Rational lower = radians.lower().toRational();
    const Rational upper = radians.upper().toRational();
    const Rational center = (lower + upper) / rational(2);
    const Rational radius = (upper - lower) / rational(2);

    const std::size_t innerBits = checkedPrecisionAdd(
        precisionBits, 16, "Trigonometric interval precision is too large");
    CertifiedTrigEnclosure point = encloseTaylorPoint(function, center, innerBits);
    const RealInterval radiusInterval = RealInterval::fromRational(radius, innerBits);
    const BigFloat radiusUpper = radiusInterval.upper();
    RealInterval expanded = add(
        point.interval,
        symmetricInterval(radiusUpper),
        innerBits);

    return CertifiedTrigEnclosure{
        expanded.roundedOutward(precisionBits),
        point.termsUsed,
        precisionBits
    };
}

[[nodiscard]] CertifiedTrigEnclosure tangentQuotient(
    const RealInterval& radians,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Tangent working precision is too large");
    const CertifiedTrigEnclosure sine = encloseFunctionAtRadianInterval(
        FunctionId::Sin, radians, workBits);
    const CertifiedTrigEnclosure cosine = encloseFunctionAtRadianInterval(
        FunctionId::Cos, radians, workBits);

    if (cosine.interval.containsZero())
        throw PrecisionInsufficient{
            "Tangent denominator cannot yet be proven nonzero"};

    return CertifiedTrigEnclosure{
        approximation::divide(sine.interval, cosine.interval, workBits)
            .roundedOutward(precisionBits),
        checkedPrecisionAdd(
            sine.termsUsed, cosine.termsUsed, "Tangent Taylor term count overflow"),
        precisionBits
    };
}

[[nodiscard]] unsigned quadrantModuloFour(const BigInt& value) {
    BigInt remainder = value % BigInt{4};
    if (remainder.isNegative())
        remainder += BigInt{4};
    if (remainder == BigInt{0}) return 0;
    if (remainder == BigInt{1}) return 1;
    if (remainder == BigInt{2}) return 2;
    return 3;
}

struct ReducedRadianInterval final {
    RealInterval remainder;
    unsigned quadrant = 0;
    std::size_t piTermsUsed = 0;
};

[[nodiscard]] ReducedRadianInterval reduceRadianArgument(
    const Rational& argument,
    std::size_t precisionBits) {
    if (argument.isZero()) {
        return ReducedRadianInterval{
            RealInterval::fromRational(argument, precisionBits), 0, 0};
    }

    // Payne-Hanek型に、x/(Pi/2) の最近整数kを保証区間から一意に確定する。
    // Rational x自体はexactなので、必要なPi精度は概ね|x|のbinary exponent + 出力精度。
    // 先にその分だけguardを積むことで、10^6等を低精度Piから何度も試す無駄を避ける。
    const Rational magnitude = argument.numerator().isNegative() ? -argument : argument;
    const auto binaryExponent = numeric::detail::floorLog2PositiveRatio(
        magnitude.numerator(), magnitude.denominator());
    const std::size_t magnitudeBits = binaryExponent > 0
        ? static_cast<std::size_t>(binaryExponent)
        : 0;
    std::size_t reductionBits = checkedPrecisionAdd(
        precisionBits, 40, "Trigonometric reduction precision is too large");
    reductionBits = checkedPrecisionAdd(
        reductionBits, magnitudeBits, "Trigonometric reduction precision is too large");

    for (;;) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        const CertifiedConstantResult pi = enclosePi(reductionBits);
        const Rational piLower = pi.interval.lower().toRational();
        const Rational piUpper = pi.interval.upper().toRational();
        const Rational twiceArgument = argument * rational(2);

        Rational quotientLower;
        Rational quotientUpper;
        if (argument.numerator().isNegative()) {
            quotientLower = twiceArgument / piLower;
            quotientUpper = twiceArgument / piUpper;
        }
        else {
            quotientLower = twiceArgument / piUpper;
            quotientUpper = twiceArgument / piLower;
        }

        const Rational half = rational(1, 2);
        const BigInt lowerNearest = numeric::floorToInteger(quotientLower + half);
        const BigInt upperNearest = numeric::floorToInteger(quotientUpper + half);
        if (lowerNearest != upperNearest) {
            const std::size_t growth = std::max<std::size_t>(32, reductionBits / 2);
            reductionBits = checkedPrecisionAdd(
                reductionBits, growth, "Trigonometric reduction precision is too large");
            continue;
        }

        const BigInt k = lowerNearest;
        const RealInterval x = RealInterval::fromRational(argument, reductionBits);
        const RealInterval halfPi = multiply(
            pi.interval,
            RealInterval::fromRational(half, reductionBits),
            reductionBits);
        const RealInterval kInterval = RealInterval::fromRational(
            Rational{k}, reductionBits);
        const RealInterval multiple = multiply(halfPi, kInterval, reductionBits);
        const RealInterval remainder = subtract(x, multiple, reductionBits);

        return ReducedRadianInterval{
            remainder, quadrantModuloFour(k), pi.termsUsed};
    }
}

[[nodiscard]] std::size_t rationalMagnitudeBits(const Rational& value) {
    if (value.isZero())
        return 0;
    const Rational magnitude = value.numerator().isNegative() ? -value : value;
    const auto exponent = numeric::detail::floorLog2PositiveRatio(
        magnitude.numerator(), magnitude.denominator());
    return exponent > 0 ? static_cast<std::size_t>(exponent) : 0;
}

[[nodiscard]] ReducedRadianInterval reduceRadianInterval(
    const RealInterval& argument,
    std::size_t precisionBits) {
    const Rational lower = argument.lower().toRational();
    const Rational upper = argument.upper().toRational();
    const std::size_t magnitudeBits = std::max(
        rationalMagnitudeBits(lower), rationalMagnitudeBits(upper));

    std::size_t reductionBits = checkedPrecisionAdd(
        precisionBits, 40, "Trigonometric reduction precision is too large");
    reductionBits = checkedPrecisionAdd(
        reductionBits, magnitudeBits, "Trigonometric reduction precision is too large");

    const CertifiedConstantResult pi = enclosePi(reductionBits);
    const Rational piLower = pi.interval.lower().toRational();
    const Rational piUpper = pi.interval.upper().toRational();
    const Rational twiceLower = lower * rational(2);
    const Rational twiceUpper = upper * rational(2);

    // Pi>0なので、x<0側だけ分母の大小による向きが反転する。
    const Rational quotientLower = lower.numerator().isNegative()
        ? twiceLower / piLower
        : twiceLower / piUpper;
    const Rational quotientUpper = upper.numerator().isNegative()
        ? twiceUpper / piUpper
        : twiceUpper / piLower;

    const Rational half = rational(1, 2);
    const BigInt lowerNearest = numeric::floorToInteger(quotientLower + half);
    const BigInt upperNearest = numeric::floorToInteger(quotientUpper + half);
    if (lowerNearest != upperNearest) {
        // 入力区間自体がPi/4境界を跨ぐ可能性がある。ここで勝手に象限を選ばず、
        // callerにworking precisionを上げて元式の区間を狭めてもらう。
        throw PrecisionInsufficient{
            "Trigonometric argument interval does not yet determine one quadrant"};
    }

    const BigInt k = lowerNearest;
    const RealInterval x = argument.roundedOutward(reductionBits);
    const RealInterval halfPi = multiply(
        pi.interval,
        RealInterval::fromRational(half, reductionBits),
        reductionBits);
    const RealInterval kInterval = RealInterval::fromRational(
        Rational{k}, reductionBits);
    const RealInterval remainder = subtract(
        x, multiply(halfPi, kInterval, reductionBits), reductionBits);

    return ReducedRadianInterval{
        remainder, quadrantModuloFour(k), pi.termsUsed};
}

[[nodiscard]] CertifiedTrigEnclosure mapReducedSinCos(
    FunctionId function,
    const ReducedRadianInterval& reduced,
    std::size_t precisionBits) {
    FunctionId baseFunction = FunctionId::Sin;
    bool negative = false;

    if (function == FunctionId::Sin) {
        switch (reduced.quadrant) {
        case 0: baseFunction = FunctionId::Sin; negative = false; break;
        case 1: baseFunction = FunctionId::Cos; negative = false; break;
        case 2: baseFunction = FunctionId::Sin; negative = true; break;
        default: baseFunction = FunctionId::Cos; negative = true; break;
        }
    }
    else {
        switch (reduced.quadrant) {
        case 0: baseFunction = FunctionId::Cos; negative = false; break;
        case 1: baseFunction = FunctionId::Sin; negative = true; break;
        case 2: baseFunction = FunctionId::Cos; negative = true; break;
        default: baseFunction = FunctionId::Sin; negative = false; break;
        }
    }

    CertifiedTrigEnclosure result = encloseFunctionAtRadianInterval(
        baseFunction, reduced.remainder, precisionBits);
    result.termsUsed = checkedPrecisionAdd(
        result.termsUsed,
        reduced.piTermsUsed,
        "Trigonometric term count overflow");
    if (negative)
        result.interval = negate(result.interval);
    return result;
}

[[nodiscard]] CertifiedTrigEnclosure encloseReducedRadian(
    FunctionId function,
    const Rational& argument,
    std::size_t precisionBits) {
    if (argument.isZero()) {
        if (function == FunctionId::Sin || function == FunctionId::Tan)
            return encloseTaylorPoint(FunctionId::Sin, argument, precisionBits);
        return encloseTaylorPoint(FunctionId::Cos, argument, precisionBits);
    }

    const ReducedRadianInterval reduced = reduceRadianArgument(argument, precisionBits);
    if (function == FunctionId::Sin || function == FunctionId::Cos)
        return mapReducedSinCos(function, reduced, precisionBits);

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Tangent working precision is too large");
    const CertifiedTrigEnclosure sine = mapReducedSinCos(
        FunctionId::Sin, reduced, workBits);
    const CertifiedTrigEnclosure cosine = mapReducedSinCos(
        FunctionId::Cos, reduced, workBits);
    if (cosine.interval.containsZero())
        throw PrecisionInsufficient{"Tangent denominator cannot yet be proven nonzero"};

    return CertifiedTrigEnclosure{
        approximation::divide(sine.interval, cosine.interval, workBits)
            .roundedOutward(precisionBits),
        checkedPrecisionAdd(
            sine.termsUsed, cosine.termsUsed,
            "Tangent Taylor term count overflow"),
        precisionBits};
}

[[nodiscard]] CertifiedTrigEnclosure encloseTurns(
    FunctionId function,
    const Rational& turns,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Trigonometric precision must be at least one bit");

    // 周期と象限はPiを近似する前にexactなRational turn上で処理する。
    // 例えば 1000000000000000001 Deg のような巨大Degreeでも、まず整数剰余だけで周期を除去できる。
    // これは巨大radian値から近似Piを何回引くか推測するより強い。
    const ReducedTrigAngle reduced = mathematics::reduceTrigTurns(function, turns);

    // tanのreference=1/4 turnはexact pole。区間精度不足とは区別して確定できる。
    if (function == FunctionId::Tan && reduced.referenceTurns == rational(1, 4))
        throw std::domain_error("tan is undefined where cos is zero");

    // referenceTurnsは[0,1/4]なので、対応radianは 2*t*Pi。Piもscaleもintervalとして外向き演算し、角度変換そのものの誤差も包含する。
    const CertifiedConstantResult pi = enclosePi(precisionBits);
    const Rational scale = reduced.referenceTurns * rational(2);
    const RealInterval scaleInterval = RealInterval::fromRational(scale, precisionBits);
    const RealInterval radians = multiply(pi.interval, scaleInterval, precisionBits);

    CertifiedTrigEnclosure result = function == FunctionId::Tan
        ? tangentQuotient(radians, precisionBits)
        : encloseFunctionAtRadianInterval(function, radians, precisionBits);
    result.termsUsed += pi.termsUsed;
    if (reduced.negative)
        result.interval = negate(result.interval);
    return result;
}

[[nodiscard]] std::optional<DecimalApproximation> tryCertifiedDecimal(
    const RealInterval& interval,
    std::size_t fractionalDigits) {
    return DecimalApproximation::fromCertifiedInterval(
        interval.lower().toRational(),
        interval.upper().toRational(),
        fractionalDigits);
}

template <class Encloser>
[[nodiscard]] CertifiedTrigResult approximateWithRetry(
    Encloser&& encloser,
    std::size_t fractionalDigits) {
    if (fractionalDigits == 0)
        throw std::invalid_argument("Trigonometric precision must be greater than zero");

    ApproximationContext context{fractionalDigits};
    for (;;) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        try {
            const CertifiedTrigEnclosure enclosure = encloser(context.workingBinaryBits());
            if (const auto decimal = tryCertifiedDecimal(enclosure.interval, fractionalDigits))
                return CertifiedTrigResult{*decimal, enclosure.termsUsed};
        }
        catch (const PrecisionInsufficient&) {
            // 現precisionではcosの非零性をまだ証明できないだけなので再試行する。
        }

        // 最終的な収束条件は「区間両端を要求桁へ丸めた結果が同一」であること。
        // 一致しなければ、丸め境界を跨いでいる可能性がまだ残るためprecisionを増やす。
        context.setGuardDigits(nextGuardDigits(context.guardDigits()));
    }
}

} // namespace

CertifiedTrigEnclosure encloseSinRadian(
    const Rational& argument,
    std::size_t precisionBits) {
    /*
    旧実装:
        return encloseTaylorPoint(FunctionId::Sin, argument, precisionBits);

    巨大な生radian値をそのままTaylorへ渡すと必要項数が|x|に比例して爆発する。
    Pi/2の最近整数倍をcertifiedに除去し、[-Pi/4,Pi/4]近傍へ落としてから評価する。
    */
    return encloseReducedRadian(FunctionId::Sin, argument, precisionBits);
}

CertifiedTrigEnclosure encloseCosRadian(
    const Rational& argument,
    std::size_t precisionBits) {
    /* 旧実装: return encloseTaylorPoint(FunctionId::Cos, argument, precisionBits); */
    return encloseReducedRadian(FunctionId::Cos, argument, precisionBits);
}

CertifiedTrigEnclosure encloseTanRadian(
    const Rational& argument,
    std::size_t precisionBits) {
    /*
    旧実装では巨大radian intervalを直接sin/cos Taylorへ渡していた。
    reduction後の同じ小区間からsin/cosを作り、pole判定を保ったまま商を取る。
    */
    return encloseReducedRadian(FunctionId::Tan, argument, precisionBits);
}

CertifiedTrigEnclosure encloseSinRadianInterval(
    const RealInterval& argument,
    std::size_t precisionBits) {
    /*
    旧実装ではpoint interval以外を常に巨大radianのままTaylor評価していた。
    小さい区間は従来Taylorを維持し、|x|が大きい場合だけPi/2の同一象限を
    区間全体で証明してから小区間へ縮約する。
    */
    if (argument.isPoint())
        return encloseReducedRadian(
            FunctionId::Sin, argument.lower().toRational(), precisionBits);

    const Rational lower = argument.lower().toRational();
    const Rational upper = argument.upper().toRational();
    if (lower >= rational(-4) && upper <= rational(4))
        return encloseFunctionAtRadianInterval(FunctionId::Sin, argument, precisionBits);

    return mapReducedSinCos(
        FunctionId::Sin, reduceRadianInterval(argument, precisionBits), precisionBits);
}

CertifiedTrigEnclosure encloseCosRadianInterval(
    const RealInterval& argument,
    std::size_t precisionBits) {
    if (argument.isPoint())
        return encloseReducedRadian(
            FunctionId::Cos, argument.lower().toRational(), precisionBits);

    const Rational lower = argument.lower().toRational();
    const Rational upper = argument.upper().toRational();
    if (lower >= rational(-4) && upper <= rational(4))
        return encloseFunctionAtRadianInterval(FunctionId::Cos, argument, precisionBits);

    return mapReducedSinCos(
        FunctionId::Cos, reduceRadianInterval(argument, precisionBits), precisionBits);
}

CertifiedTrigEnclosure encloseTanRadianInterval(
    const RealInterval& argument,
    std::size_t precisionBits) {
    if (argument.isPoint())
        return encloseReducedRadian(
            FunctionId::Tan, argument.lower().toRational(), precisionBits);

    const Rational lower = argument.lower().toRational();
    const Rational upper = argument.upper().toRational();
    if (lower >= rational(-4) && upper <= rational(4))
        return tangentQuotient(argument, precisionBits);

    const ReducedRadianInterval reduced = reduceRadianInterval(argument, precisionBits);
    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 24, "Tangent working precision is too large");
    const CertifiedTrigEnclosure sine = mapReducedSinCos(
        FunctionId::Sin, reduced, workBits);
    const CertifiedTrigEnclosure cosine = mapReducedSinCos(
        FunctionId::Cos, reduced, workBits);
    if (cosine.interval.containsZero())
        throw PrecisionInsufficient{"Tangent denominator cannot yet be proven nonzero"};
    return CertifiedTrigEnclosure{
        approximation::divide(sine.interval, cosine.interval, workBits)
            .roundedOutward(precisionBits),
        checkedPrecisionAdd(sine.termsUsed, cosine.termsUsed, "Tangent Taylor term count overflow"),
        precisionBits};
}

CertifiedTrigEnclosure encloseSinTurns(
    const Rational& turns,
    std::size_t precisionBits) {
    return encloseTurns(FunctionId::Sin, turns, precisionBits);
}

CertifiedTrigEnclosure encloseCosTurns(
    const Rational& turns,
    std::size_t precisionBits) {
    return encloseTurns(FunctionId::Cos, turns, precisionBits);
}

CertifiedTrigEnclosure encloseTanTurns(
    const Rational& turns,
    std::size_t precisionBits) {
    return encloseTurns(FunctionId::Tan, turns, precisionBits);
}

CertifiedTrigResult approximateSin(
    const RealNumber& argument,
    std::size_t fractionalDigits) {
    const Rational radians = argument.toRational();
    return approximateWithRetry(
        [&](std::size_t bits) { return encloseSinRadian(radians, bits); },
        fractionalDigits);
}

CertifiedTrigResult approximateCos(
    const RealNumber& argument,
    std::size_t fractionalDigits) {
    const Rational radians = argument.toRational();
    return approximateWithRetry(
        [&](std::size_t bits) { return encloseCosRadian(radians, bits); },
        fractionalDigits);
}

CertifiedTrigResult approximateTan(
    const RealNumber& argument,
    std::size_t fractionalDigits) {
    const Rational radians = argument.toRational();
    return approximateWithRetry(
        [&](std::size_t bits) { return encloseTanRadian(radians, bits); },
        fractionalDigits);
}

CertifiedTrigResult approximateSinTurns(
    const Rational& turns,
    std::size_t fractionalDigits) {
    return approximateWithRetry(
        [&](std::size_t bits) { return encloseSinTurns(turns, bits); },
        fractionalDigits);
}

CertifiedTrigResult approximateCosTurns(
    const Rational& turns,
    std::size_t fractionalDigits) {
    return approximateWithRetry(
        [&](std::size_t bits) { return encloseCosTurns(turns, bits); },
        fractionalDigits);
}

CertifiedTrigResult approximateTanTurns(
    const Rational& turns,
    std::size_t fractionalDigits) {
    return approximateWithRetry(
        [&](std::size_t bits) { return encloseTanTurns(turns, bits); },
        fractionalDigits);
}

} // namespace mmcal::approximation
